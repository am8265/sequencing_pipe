/// crude hack of sam.c - please see samtools-1.5/htslib-1.5/sam.c for who/why/where/when and the proper way to do things...
#include <strings.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <zlib.h>
#include <assert.h>
#include "htslib/sam.h"
#include "htslib/bgzf.h"
// #include "cram/cram.h"
#include "hts_internal.h"
#include "htslib/hfile.h"
#include "htslib/hts_endian.h"
#include "htslib/khash.h"
KHASH_DECLARE(s2i, kh_cstr_t, int64_t)
typedef khash_t(s2i) sdict_t;
#ifndef EOVERFLOW
#define EOVERFLOW ERANGE
#endif
#define KS_SEP_LINE  2 // line separator: "\n" (Unix) or "\r\n" (Windows)
#include <stdint.h>
#include <set>
#include <string>

///https://stackoverflow.com/questions/1471353/whats-the-c-equivalent-of-uint32-max
#define INT32_MAX ((uint32_t)-1)

static bool lazy_prob = false;
static int nasty1=0,nasty2=0,nasty3=0;

bam_hdr_t *bam_hdr_read_hack(BGZF *fp) {
    bam_hdr_t *h;
    char buf[4];
    int magic_len, has_EOF;
    int32_t i, name_len, num_names = 0;
    size_t bufsize;
    ssize_t bytes;
    // check EOF
    has_EOF = bgzf_check_EOF(fp);
    if (has_EOF < 0) { 
        perror("[W::bam_hdr_read] bgzf_check_EOF");
        lazy_prob = true;
    } else if (has_EOF == 0) { 
        hts_log_error("EOF marker is absent. The input is probably truncated");
        lazy_prob = true;
    }

    // read "BAM1"
    magic_len = bgzf_read(fp, buf, 4);
    if (magic_len != 4 || strncmp(buf, "BAM\1", 4)) { //////// 4241 4d01
        hts_log_error("Invalid BAM binary header");
        lazy_prob = true;
        // exit(1);
    }
    h = bam_hdr_init();
    if (!h) goto nomem;

    // read plain text and the number of reference sequences
    bytes = bgzf_read(fp, &h->l_text, 4); /////// l_text      length of the plain text in the header
// printf("just read %d bytes (length of header?!?)\n",bytes);
    if (bytes != 4) goto read_err;
    if (fp->is_be) ed_swap_4p(&h->l_text);

    bufsize = ((size_t) h->l_text) + 1;
    if (bufsize < h->l_text) goto nomem; // so large that adding 1 overflowed

    h->text = (char*)malloc(bufsize);
    if (!h->text) goto nomem;

    h->text[h->l_text] = 0; // make sure it is NULL terminated
    bytes = bgzf_read(fp, h->text, h->l_text);
// printf("just read %d bytes (the header itself)\n",bytes);
    if (bytes != h->l_text) goto read_err;

////// here we have the 'header' as appears in terms of parsing '^@[A-Z]{2}\t' lines but the 'actual' targets header is pulled below?!?
{
int i;
for (i =0; i<h->l_text; ++i) if(memcmp(h->text+i,"\n@RG",4)==0||memcmp(h->text+i,"\n@PG",4)==0) break;
char * x = h->text+i;
// printf("fucker=%s\n------\n",x);
char ent[16*1024], * entp=ent; // char ent[1024], * entp=ent; // silly long lines?!?

//// cannot have space between define name and param?!?
// std::string who, how;

//// ternary?!?
#define RAR(X,Y) if((X).empty()) (X)=(Y); \
else (X)+=std::string(";")+(Y); 

    // (who).insert((ent)); \
   
#define ARGH(x,entp,ent,who,how) if(*++(x)=='\t') *(entp)++=','; \
else if(*(x)==' ') *(entp)++='_'; \
else if(*(x)!='\n') *(entp)++=*(x); \
else { \
    *(entp)='\0'; \
    (entp)=(ent); \
    if((ent)[2]=='G') { \
        if((ent)[1]=='R') { RAR((who),(ent)); ++nasty1; \
        } else if((ent)[1]=='P') { \
            ++nasty2; \
            int j=0, z=0; \
            for (;j<strlen(ent);++j) { if(ent[j]=='C' && ent[j+1]=='L') break;} \
            if(ent[j-1]==',') ent[j-1]='\0'; \
            ent[j]='\0'; \
            RAR((how),(ent)); \
        } \
    }  \
} 

std::string who,how;
// std::set<std::string> who,how;
for (; i <h->l_text; ++i) ARGH(x,entp,ent,who,how);
ARGH(x,entp,ent,who,how);
puts(who.data());
puts(how.data());
}
    /*std::string rg;
    for(std::set<std::string>::iterator it = who.begin(); it != who.end(); it++ ) {
        if(rg.empty()) rg+=*it; // it->data();
        else rg+=","+*it;
    } */

/*
char * x = h->text+i;
// printf("<-\n->");
char rg[1024], * rgp = rg;
x+=4;
for(; *++x!='\n'; ++i) {
    if(*x!='\t') *rgp++=*x;
    else *rgp++=',';
    // printf("%c",*x);
}
*rgp='\0';
//// should check the reads all have the correct rg tags?!?
printf("ReadGroup=%s\n",rg); // printf("ReadGroup='%s'\n",rg);
// while(++i,*++x!='\n') printf("%c",*x);
// printf("%s",h->text+i);
/// can't be bothered with flag?!?
i+=6;
// printf("<-\n->");
int c=0;
for (; i<h->l_text; ++i) {
    if(c<3) printf("%c",*(h->text+i));
    if(h->text[i]=='\t') ++c;// ,puts(""); 
    else if(h->text[i]=='\n') c=0, puts("");
}
*/

    bytes = bgzf_read(fp, &h->n_targets, 4);
// printf("just read %d bytes (the number of target seq?!? - that is, now we're ready to read the 'real' header?!?!?)",bytes);
    if (bytes != 4) goto read_err; if (fp->is_be) ed_swap_4p(&h->n_targets); if (h->n_targets < 0) goto invalid;

    // read reference sequence names and lengths
    if(h->n_targets!=86) {
        fprintf(stderr,"what was this mapped to?!?\n"); // should check scf lengths etc.?!?
        lazy_prob = true;
    }
    // printf(" - i.e. refs=%d\n",h->n_targets);

    if (h->n_targets > 0) {
        h->target_name = (char**)calloc(h->n_targets, sizeof(char*));
        if (!h->target_name) goto nomem;
        h->target_len = (uint32_t*)calloc(h->n_targets, sizeof(uint32_t));
        if (!h->target_len) goto nomem;
    } else {
        h->target_name = NULL;
        h->target_len = NULL;
    }

    for (i = 0; i != h->n_targets; ++i) {

        bytes = bgzf_read(fp, &name_len, 4); if (bytes != 4) goto read_err; if (fp->is_be) ed_swap_4p(&name_len); if (name_len <= 0) goto invalid;
// printf("\tjust read %d bytes (length of this target name?!?)\n",bytes);

        h->target_name[i] = (char*)malloc(name_len); if (!h->target_name[i]) goto nomem;

        num_names++;

// printf("\tjust read %d bytes (name itself?!?)\n",bytes);
        bytes = bgzf_read(fp, h->target_name[i], name_len); if (bytes != name_len) goto read_err;

        if (h->target_name[i][name_len - 1] != '\0') { /* Fix missing NUL-termination.  Is this being too nice? We could alternatively bail out with an error. */
            char *new_name;
            if (name_len == INT32_MAX) goto invalid;
            new_name = (char*) realloc(h->target_name[i], name_len + 1);
            if (new_name == NULL) goto nomem;
            h->target_name[i] = new_name;
            h->target_name[i][name_len] = '\0';
        }

// printf("\t\tname=%s\n",h->target_name[i]);

        bytes = bgzf_read(fp, &h->target_len[i], 4); if (bytes != 4) goto read_err; if (fp->is_be) ed_swap_4p(&h->target_len[i]);
    }
    return h;
 nomem:
    hts_log_error("Out of memory");
    goto clean;
 read_err:
    if (bytes < 0) {
        hts_log_error("Error reading BGZF stream");
    } else {
        hts_log_error("Truncated BAM header");
    }
    goto clean;

 invalid:
    hts_log_error("Invalid BAM binary header");

 clean:
    if (h != NULL) {
        h->n_targets = num_names; // ensure we free only allocated target_names
        bam_hdr_destroy(h);
    }
    return NULL;
}

#define BAM_LIDX_SHIFT    14

int main(int argc, char** argv) {

    // const char *fnidx, int min_shift, int nthreads) {
    if(argc!=2) puts("blah"),exit(1);
    const char *fn = *(argv+1);

    htsFile *fpx;
    if ((fpx = hts_open(fn, "r")) == 0);

    // hts_idx_t *idx = bam_index(fp->fp.bgzf,BAM_LIDX_SHIFT);
    // static hts_idx_t *bam_index(BGZF *fp, 
    int min_shift = BAM_LIDX_SHIFT;

    BGZF *fp = fpx->fp.bgzf;

    int n_lvls, i, fmt, ret;
    bam1_t *b;
    hts_idx_t *idx;
    bam_hdr_t *h;

    ///// go through BAM1, header length, 'header' an then 'target' length/names...?!?
    h = bam_hdr_read_hack(fp);

    // if (h==1) return 1;

    if (min_shift > 0) {
        int64_t max_len = 0, s;
        for (i = 0; i < h->n_targets; ++i)
            if (max_len < h->target_len[i]) max_len = h->target_len[i];
        max_len += 256;
        for (n_lvls = 0, s = 1<<min_shift; max_len > s; ++n_lvls, s <<= 3);
        fmt = HTS_FMT_CSI;
    } else min_shift = 14, n_lvls = 5, fmt = HTS_FMT_BAI;

    // idx = hts_idx_init(h->n_targets, fmt, bgzf_tell(fp), min_shift, n_lvls);
    bam_hdr_destroy(h);
    b = bam_init1();
    int count=0;
///// can make this a dody md5sum?!?
    std::set<std::string> who;
    while ((ret = bam_read1(fp, b)) >= 0 && ++count<=100000) {
        char const * qn = bam_get_qname(b);
        // assert(memcmp(qn,"HWI-",4)==0);
        int y = 0, z =0;
        if((memcmp(qn,"HWI-",4)==0)) y = z = 4;
        for (int g = 0; y < strlen(qn); ++y ) {
            if(*(qn+y)==':') ++g;
            else if(g>=4) break;
        }
        char mc_rn_fc_ln[1024];
        memcpy(mc_rn_fc_ln,qn+z,y-z);
        mc_rn_fc_ln[y-z-1]='\0';
        who.insert(mc_rn_fc_ln);
        // puts(mc_rn_fc_ln);
        // printf("have %s : %d\n",mc_rn_fc_ln,y);
        // ret = hts_idx_push(idx, b->core.tid, b->core.pos, bam_endpos(b), bgzf_tell(fp), !(b->core.flag&BAM_FUNMAP));
        if (ret < 0) lazy_prob = true; // goto err; // unsorted
    }
    std::string rg;
    for(std::set<std::string>::iterator it = who.begin(); it != who.end(); it++ ) {
        if(rg.empty()) rg+=*it; // it->data();
        else rg+=","+*it;
        ++nasty3;
    }
    rg="SampleRg="+rg;
    if(rg.length()>250) printf("%s...\n",rg.substr(0,247).data());
    else puts(rg.data());
    if (ret < -1) goto err; // corrupted BAM file
    printf("%d:%d:%d\n",nasty1,nasty2,nasty3);
// if(nasty1!=nasty3) lazy_prob = true;
/////// remove this but make it check seqdb for FCID & LANE, make sure read names in header and body match too!?!?
if(nasty3>200) {
    puts("DODGY_OLD_BAM");
    // lazy_prob = true;
}
    // hts_idx_finish(idx, bgzf_tell(fp));
    bam_destroy1(b);
    if(lazy_prob) puts("ERROR");
    else puts("OKAY");
    return lazy_prob;
err:
    bam_destroy1(b);
    // hts_idx_destroy(idx);
    return 1;
}


