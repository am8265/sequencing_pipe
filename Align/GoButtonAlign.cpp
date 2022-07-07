#include <iostream>
#include <utility>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <stdlib.h>
#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <openssl/md5.h>
#include <stdint.h>
#include "rapidjson/document.h"     
#include "rapidjson/prettywriter.h" 
#include "rapidjson/error/en.h" // for stringify JSON
#include <curl/curl.h>
#include <string>
#include <EXTERN.h>
#include <perl.h>
#include <my_global.h>
#include <mysql.h>

#include <signal.h>
#include <execinfo.h>
// #ifdef _OPENMP // when was this used?
#include <omp.h>
// #endif
#undef min
#undef max

// #ifndef RELEASE_H
// #define RELEASE_H

// use vector of pairs of macro set funcname/func-pointer and shuffle
// put in those absurdly hacking merge recover things?!? 
// just add ethnicity here!
// add merge rescue object too?!?
// add reset to bcl BUT with lock issue telling them to run it ONLY if they're sure!?!? : update Flowcell set fc_status = 'registered' where FCillumID = 'HMJCJDSXX' and fc_status = 'sequenced' and fail = 0 ; select row_count() locked

// #define BASE_DIR "/nfs/seqscratch_ssd/informatics/"
// hc: never reference BASE_DIR, update anyway. 
#define BASE_DIR "/nfs/seqscratch_ssd/informatics/"
#define SCRIPT_DIR BASE_DIR "logs/merge/"
#define LOG_DIR "/nfs/central/home/dh2880/.logs/"
#define ALIGNSTATS "/nfs/goldstein/software/alignstats/alignstats"
// #define POSTMERGE "/nfs/seqscratch_ssd/informatics/logs/postmerge/"
#define POSTMERGE "/nfs/seqscratch_ssd/informatics/logs/postmerge/"

// int merge_and_release(int, char **);
// inline bool isregfile(const char* fn) { struct stat test; if (stat(fn, &test) != 0) { return false; }  return S_ISREG(test.st_mode); }                                   

using namespace std;


bool myfunction (float i,float j) { return (i<j); }

void finish_with_error(MYSQL *con) { // http://zetcode.com/db/mysqlc/
    fprintf(stderr, "%s\n", mysql_error(con));
    mysql_close(con);
    exit(1);        
}

struct FUNKY { std::string experiment_id, prepid, sample_name, sample_type, capture_kit, priority, is_external, end_point; };

#define FDP_OFFSET 11
// #define FDP_OFFSET 15



char const * REPLACE_INTO_DPS = 
  "replace into dragen_pipeline_step "
  "(pseudo_prepid,  pipeline_step_id,   version,    submit_time,            finish_time,            times_ran,  step_status) values "
  "(%s,             1,                  '0.0.1',    CURRENT_TIMESTAMP(),    CURRENT_TIMESTAMP(),    1,          '%s');";

inline bool isregfile(const char* fn) { struct stat test; if (stat(fn, &test) != 0) { return false; }  return S_ISREG(test.st_mode); }
inline bool isdir(const char* fn) { struct stat test; if (stat(fn, &test) != 0) { return false; }  return S_ISDIR(test.st_mode); }
inline void touchfile(const std::string& file) { FILE* fp = fopen(file.data(), "ab+"); fclose(fp); }       



namespace options { bool debug = false; }

// MUST STOP HAVING DPS/DQM ERRORS IF IT'S IN WALDB ALREADY. WIPE IT ALL AND SEND EMAIL - SHOULD DEPRECATE!?!?!?
enum status { 
    WGS = -5, INCONSISTENT = -6, MERGE_ERROR = -7, PROB = -8, DQM = -9, JUNK = -10, 
    // must split COMPONENT_ERROR into COMPONENT_RG_ERROR, COMPONENT_EOF_ERROR...
    MISSING_FCINFO = -11, RG_COUNT_ERROR = -12, DSP = -13, COMPONENT_ERROR = -14,
    MOVE = -15, NO_DIR = -16, NULL_COMPONENTS = -17,
    MERGE_RESCUE = -18, MERGE_RESCUE_ERROR = -19,
    SCRIPT_EXISTS = -20, CHECKPOINT_EXISTS = -21, RG_PU_VERSUS_READNAME_ERROR = -22, RG_LIMS_READNAME_MISMATCH_ERROR = -23,
    DT_TAGS = -24,
    MERGE_METRICS_RESCUE = -28, MERGE_METRICS_RESCUE_ERROR = -29,
    COMPONENT_RG_ERROR = -30, COMPONENT_EOF_ERROR = -31, FEW_MT = -32,
    GAFE_NO_DATA = -52, GAFE_ALREADY_HAS_DATA = -51, 
    FINAL_BAM_MISSING = -60, BAM_CHECKING_PROB = -61, BAM_IS_A_MESS = -62, BAM_LIMS_SAMPLE_NAME_MISMATCH = -63, 
    BAM_COUNT_MISMATCH = -64, NOT_A_BAM = -65, EXTERNAL_STATUS_SCREWED = -66, DB_SCREWED = -67, MERGE_METRICS_ERROR = -70
}; 

void tokenise(std::vector<std::string>& t, std::string const l, char s) {
// void tokenise(std::vector<std::string>& t, std::string const& l, char s) {

    using namespace std;
    size_t pos = 0;
    size_t lst = 0;

    while (pos != std::string::npos) {
        pos = l.find(s, pos + 1);
        size_t g = lst == 0 ? 0 : 1;
        string a = l.substr(lst + g, pos - lst - g);
        // cout << a << "\n";
        for (unsigned c=0 ;c<a.length(); ++c) {
            // cout << "["<<c<<"]" <<a[c] << "\n";
        }
        t.push_back(a);
        lst = pos;
    }
}

class Popen {
    public: 
    Popen(char const * cmd, int const z,char const * m) : _m(m), _ml(z), _bf(new char [_ml]) { 
        if(!(_in = popen(cmd,_m))) std::cout << "problem with command\n", exit(1); 
    }
    Popen();
    ~Popen() { fflush(_in); pclose(_in); delete [] _bf; }

    void write(std::string const & a) { 
    // void write(char const * a) {

        assert(_m[0]=='w');
        // bored of warning
        // fwrite(a.data(),sizeof(char),a.length(),_in);
        assert(fwrite(a.data(),sizeof(char),a.length(),_in));
        fflush(_in);

    }

    char* getline() const {
        assert(_m[0]=='r');
        char *hmm = _bf; int x = 0;
        for( x=0, hmm=_bf; (*hmm=fgetc(_in))!=EOF && *hmm++!='\n'; x++){ assert(x<_ml); }
        assert(!(_bf[0]=='\t'&&_bf[1]=='\t'));
        _bf[x]='\0';
        return _bf; 
    }

  private:
    char const  *_m;
    int const   _ml;
    char        *_bf;
    FILE        *_in;
};

int filesize(char const * p) {
    struct stat a;
    stat(p,&a);
    return a.st_size;
}

// implicit ctor in string or better yet use operator()?!?
template<typename A> struct Yum { // struct Yum {
    Yum(char const * cmd) : output(cmd) {} // Yum(char const * cmd, A a) {
    std::string operator()(A a){ // std::string operator()(char const * a){
        char tmp[1024]; sprintf(tmp,output.data(),a); 
        Popen ox(tmp,16*1024,"r"); 
        output=ox.getline(); 
        return output; 
    }
    void operator()(std::vector<std::string> & X,A a){ 
        char tmp[1024]; sprintf(tmp,output.data(),a); 
        Popen ox(tmp,16*1024,"r"); 
        char * z3;
        while( *( z3=ox.getline() ) != '\0') X.push_back(z3);
    }
    std::string output;
};

inline void Lazy(char const * a, char const * b) {
    char tmp[1024]; sprintf(tmp,a,b);
// std::cout << "LAZY : using\n'"<<tmp<<"'\n";
    if(system(tmp)) std::cout << "what '" << tmp << "'",exit(1);
}

// https://stackoverflow.com/questions/9317305/sending-an-email-from-a-c-c-program-in-linux
// should just open a socket to local port and send direct...?!?
int sendmail(const char *to, const char *from, const char *subject, const char *message, bool html = false) {

    int retval = -1;
    FILE *mailpipe = popen("/usr/sbin/sendmail -t", "w");
    if (mailpipe != NULL) {
        fprintf(mailpipe, "To: %s\n", to);
        fprintf(mailpipe, "From: %s\n", from);
        time_t t = time(0);
        struct tm * tm_s = localtime(&t);
        char bits[1024], n[256], blah[1024];
        strftime(blah,1024,"%c",tm_s);
        if(html) fprintf(mailpipe, "Content-Type: text/html\n");
        fprintf(mailpipe, "Subject: %s (%s)\n\n",subject,blah); // sick of merging subjects...?!?
        gethostname(n,256);
        sprintf(bits,"%s:%d : %s",n,getpid(),subject);
        // bored of warning
        // fwrite(message, 1, strlen(message), mailpipe);
        // fwrite(".\n", 1, 2, mailpipe);
        assert(fwrite(message, 1, strlen(message), mailpipe));
        assert(fwrite(".\n", 1, 2, mailpipe));
        pclose(mailpipe);
        retval = 0;
     }else{
         perror("Failed to invoke sendmail");
     }

     return retval;
}

namespace lazy {

    using namespace std;

    inline std::string GetFirstLinePopen(char const * const a) {

        Popen ns1(a,16*1024,"r"); 
        // cout << "GetFirstLinePopen=\""<<a<<"\"\n\n";

        // this was an error - duh! returning address 
        // return ns1.getline();
        return std::string(ns1.getline());
    }

    template<typename T1> inline std::string GetFirstLinePopen(char const * const a, T1 t1) {
        char tmp[1024]; sprintf(tmp,a,t1); return GetFirstLinePopen(tmp);
    }

    // this 'could' use next one down - i.e. just recurse instead of collapsing
    template<typename T1, typename T2> inline std::string GetFirstLinePopen(char const * const a, T1 t1, T2 t2) {
        char tmp[1024]; sprintf(tmp,a,t1,t2); return GetFirstLinePopen(tmp);
    }

    template<typename T1, typename T2, typename T3> inline std::string GetFirstLinePopen(char const * const a, T1 t1, T2 t2, T3 t3) {
        char tmp[1024]; sprintf(tmp,a,t1,t2,t3); return GetFirstLinePopen(tmp);
    }

    template<typename T1, typename T2, typename T3, typename T4> inline std::string GetFirstLinePopen(char const * a, T1 t1, T2 t2, T3 t3, T4 t4) {
        char tmp[1024]; sprintf(tmp,a,t1,t2,t3,t4); return GetFirstLinePopen(tmp);
    }

}

inline std::string get_single_line_output_as_string(char const * const a, char const * const b) {
    char tmp[1024];
    sprintf(tmp,a,b);
    Popen ns1(tmp,16*1024,"r"); 
    return ns1.getline();
}

inline bool run_line(char const * const a, char const * const b) {
    char tmp[1024];
    sprintf(tmp,a,b);
    if(system(tmp)!=0) {             std::cout << "what : " << tmp << "\n\n\n";               exit(1);                }
    return true;
}



namespace fastq { 

    static std::map<std::string,std::string> FCT;

    /* inline */ std::string fc(char const * const fn) {

        using namespace std;
        string ret;
        if(strlen(fn)<9) return ret;

        {
        for(unsigned int o=7;o<strlen(fn);++o){

            if(fn[o]=='X'
            // && fn[o+1]=='X'
            && ( ( (int)fn[o+1]>=(int)'A' && (int)fn[o+1]<=(int)'Z') || ( (int)fn[o+1]>=(int)'0' && (int)fn[o+1]<=(int)'9') )
            ) {

                bool fc=true;
                string lazy;

                for (unsigned int p=o-7;p<o+2;++p) {
                    // if( ( (int)fn[p]>=65&&(int)fn[p]<=90) || ( (int)fn[p]>=48&&(int)fn[p]<=57) ) {
                    if( ( (int)fn[p]>=(int)'A' && (int)fn[p]<=(int)'Z') || ( (int)fn[p]>=(int)'0' && (int)fn[p]<=(int)'9') ) {
                        lazy+=fn[p];
                    }else fc=false;
                }

                if((int)lazy[0]>=(int)'A' && (int)lazy[0]<=(int)'I') {
                }else fc=false;

                if(lazy.substr(0,2)=="EU")fc=false;
                if(lazy.substr(0,6)=="ALSNEU")fc=false;
                // assert((int)lazy[0]>=(int)'A' && (int)lazy[0]<=(int)'I');

                if(fc) {
                    return lazy;
                }

                cout << "\n";
            }
        }
        return ret;
        }
    }

    bool rncheck(char const * const rn,std::string &fcs, std::string &bc, std::string &mtype_from_fc, std::string &mtype_from_pfx) {

        using namespace std;

        if(FCT.size()==0){

            FCT.insert(make_pair("ACXX","HiSeq High-Output (8-lane) v3 flowcell (HiSeq 1000-2500)"));
            FCT.insert(make_pair("ANXX","HiSeq High-Output (8-lane) v4 flowcell (HiSeq 1500-2500)"));
            FCT.insert(make_pair("BGXX","High-Output NextSeq"));
            FCT.insert(make_pair("BGXY","High-Output NextSeq"));
            FCT.insert(make_pair("BGX2","High-Output NextSeq"));
            FCT.insert(make_pair("BGX3","SUDC Sample : Guessing it's High-Output NextSeq"));
            FCT.insert(make_pair("AFXX","Mid-Output NextSeq"));
            FCT.insert(make_pair("ADXY","Old_Perhaps"));
            FCT.insert(make_pair("BCX2","Old_Perhaps"));
            FCT.insert(make_pair("BBXX","HiSeq 4000 (8-lane) v1 flowcell"));
            FCT.insert(make_pair("BBXY","HiSeq 4000 (8-lane) v1 flowcell"));
            FCT.insert(make_pair("ABXX","HiSeq 2000"));
            FCT.insert(make_pair("DRXX","NovaSeq S1 flowcell"));
            FCT.insert(make_pair("DMXX","NovaSeq S2 flowcell"));
            FCT.insert(make_pair("BCX3","NovaSeq S2 flowcell (hack)"));
            FCT.insert(make_pair("DSXX","NovaSeq S4 flowcell"));
            FCT.insert(make_pair("DSXY","NovaSeq S4 flowcell (hack)"));
            FCT.insert(make_pair("BGX5","NovaSeq S4 flowcell (hack EGI)")); 
            FCT.insert(make_pair("BGX7","NovaSeq S4 flowcell (hack EGI)"));
	    FCT.insert(make_pair("DSX2","NovaSeq S4 flowcell (hack EGI)"));
            FCT.insert(make_pair("ALXX","HiSeqX (8-lane) flowcell"));
            FCT.insert(make_pair("CCXX","HiSeqX (8-lane) flowcell"));
            FCT.insert(make_pair("CCXY","HiSeqX (8-lane) flowcell"));
            FCT.insert(make_pair("AAXX","Genome Analyzer"));
            FCT.insert(make_pair("ADXX","HiSeq Rapid Run (2-lane) v1 (HiSeq 1500/2500)"));
            FCT.insert(make_pair("AGXX","High-Output NextSeq"));
            FCT.insert(make_pair("AMXX","HiSeq RR v2"));
            FCT.insert(make_pair("BCXX","HiSeq Rapid Run (2-lane) v1.5/v2 (HiSeq 1500/2500)"));
            FCT.insert(make_pair("BCXY","HiSeq Rapid Run (2-lane) v2 (HiSeq 1500/2500)"));
            FCT.insert(make_pair("DRXY","NovaSeq S4 flowcell"));
            FCT.insert(make_pair("DSX3","NovaSeq S4 flowcell"));
        }

        std::vector<std::vector<std::string> > info;
        std::vector<std::string> y;
        tokenise(y,rn,' ');
        for (unsigned int i=0;i<y.size();++i) {
            std::vector<std::string> x;
            tokenise(x,y[i],':');
            info.push_back(x);
        } 

        assert(info[0][0][0]=='@');
        info[0][0]=info[0][0].substr(1);//,info[0][0].length()-1);
        bool bcl2=false;
        if(info.size()==2) {
            assert(info[0].size()==7);
            assert(info[1].size()==4);
            ++bcl2;
            // ss << "THIS APPEARS TO BE ORIGINAL bcl2fastq...\n";
            bc=info[1][3];
        }

        int fcbcc=0;
        string bored;
        bool one8plus=false;
        for (unsigned int i=0;i<info.size();++i) { 
            for (unsigned int j=0;j<info[i].size();++j) { 
                string fcst;
                if(i==0 && !(fcst=fc(info[i][j].data())).empty() ) {
                    ++fcbcc;
                    fcs=fcst;
                    if(j==2) one8plus=true;
                    break;
                }
            }
        }
        assert(fcbcc<=1);

        std::stringstream ss;
        if(info[0][0].substr(0,5)=="HWUSI") { ss << "GA IIx\n";
        }else if(info[0][0].substr(0,5)=="HWI-M") { ss << "MiSeq\n"; 
        }else if(info[0][0].substr(0,5)=="HWI-C") { ss << "HiSeq (1500)\n";
        }else if(info[0][0].substr(0,5)=="HWI-S") { ss << "HiSeq (2000)\n";
        }else if(info[0][0].substr(0,5)=="HWI-D") { ss << "HiSeq (2500)\n";
        //////
        }else if(info[0][0].substr(0,5)=="HWI-E") { ss << "Old GAII/HiSeq?!?\n";
        //// should check all following are digits!?!
        }else if( ( info[0][0].substr(0,2)=="NB" || info[0][0].substr(0,2)=="NS" ) && ::isdigit(info[0][0][2]) ) { ss << "NextSeq\n";
        }else if( info[0][0].substr(0,2)=="MN" && ::isdigit(info[0][0][2]) ) { ss << "MiniSeq\n";
        }else if( info[0][0][0]=='D' && ::isdigit(info[0][0][1]) ) { ss << "HiSeq 2500\n";
        }else if( info[0][0][0]=='E' && ::isdigit(info[0][0][1]) ) { ss << "HiSeqX\n";
        }else if( info[0][0][0]=='J' && ::isdigit(info[0][0][1]) ) { ss << "HiSeq 3000\n";
        }else if( info[0][0][0]=='K' && ::isdigit(info[0][0][1]) ) { ss << "HiSeq 3000/4000\n";
        }else if( info[0][0][0]=='C' && ::isdigit(info[0][0][1]) ) { ss << "HiSeq 1500\n";
        }// else{ ss << "no_idea\n"; }

        mtype_from_pfx=ss.str();

        if(!fcs.empty()) {
            if(FCT.count(fcs.substr(5,4))==0) cout << "unknown chemistry : " << fcs << "\n", exit(1);
            else mtype_from_fc=FCT[fcs.substr(5,4)];
        }

        bool one4=true;
        if(fcs.empty()) {
            if(info[0].size()!=5) one4=false;
            if(info[0].size()==3) return false;
            if(one4) for(int o=1;o<4;++o){
                string &c=info[0][o];
                for(unsigned int g=1;g<c.length();++g){
                    if(!::isdigit(c[g])) {
                        one4=false;
                        break;
                    }
                }
            }
            if(info[0].size()<3) one4 =false;
            else if(strchr(info[0][4].data(),'#')==0) one4=false;
        }else one4=false;

        if(!one8plus){ 
            if(one4) cout << "classic 1.4\n";
            else if(info[0].size()==7 && info[0][2]=="FC") {
                cout << "WHAT:1: some form of nasty 1.4->1.8 hack?!? - eeeeew\n";
            } else if(info[0].size()==7 && ( info[0][2][0]=='h' || info[0][2][0]=='c' ) &&  info[0][2][7]=='x' ) {
                for(unsigned int y=0;y<info[0][2].length();++y) info[0][2][y]=::toupper(info[0][2][y]);
                cout <<"WHAT: looks like some whatwit has been at this?!? ="<<info[0][2] << "\n";
                fcs=info[0][2];
            } else if(info[0].size()==7) {
                cout << "WHAT:2: some other form of 1.4->1.8 hack?!? - yuck\n";
                // exit(1);
            }else if(!fcs.empty() && !one8plus){
                cout << "WHAT1:3: this is a strange format in which someone has either merged FC/Machine or happens to have a legit' FCID in machine name?!?\n";
                // exit(1);
            } else {
                cout << "WHAT:4: what in this?!?\n";
                // assert(0);
            }
        }

        return (!fcs.empty() && one8plus);
    }

}


#define SINLGE_CMD(N,X,Y,Z) (N).fill_in_name((X),(Y),sizeof((Y))); (Z).push_back((Y)); run_cmd((Z)); (Z).clear();
#define SINLGE_CMD2(N,X,Y,Z) { char b1[1024]; strcpy(b1,(Z)); SINLGE_CMD(N,b1,(X),(Y)) }

// #endif

#define HOW_OFTEN 10

#define LOCK_T "LOCKED_SAMPLES"
#define LOCK_Q "update dragen_sample_metadata set is_merged = %d where is_merged = %d and experiment_id = %s ; select row_count() " LOCK_T
#define LOCK_IT(Y,Z,A,X) \
    { NLIST updated = db::get_named_row("seqdb", (X), (Z), (Y), (A) ) ; \
    assert(updated[LOCK_T]=="1"); \
    cout << "got " << updated[LOCK_T] << "\n"; }

#define ALIGNSTATS "/nfs/goldstein/software/alignstats/alignstats"
#define VERIFYBAMID "/nfs/goldstein/software/bin/verifyBamID.20120620"

// #include "pull_info.h" // this is just wrong!?!

namespace seq {
    static char const * NETAPP_OUT_DIR = "/nfs/seqscratch_ssd/",
      * SEQ_RUN_DIR = "/nfs/hts/novaseq/sequence/Runs/";
}

// this is horrible!?!

/// namespace { // unnamed!?!

#define MIN_RGS 5

#define GET_FLOWCELL "select * from Flowcell where FCillumID = '%s' and fail = 0"

struct ARCHIVE { std::string data, sample_internal_name, lanenum, archive_dir, html, read1, read2; };

static std::map<std::string, std::map<std::string, int> > insane_frac;

inline void mkdirs(char const *bdir) {
    if(!isdir(bdir)) {
        // std::cout << "need to create dir " << bdir << "\n";
        char cmd[2048];
        sprintf(cmd,"mkdir -p %s",bdir);
        // assert(mkdir(bdir,0x664)==0);
        assert(system(cmd)==0);
    }
}

inline void link_file(char const *adir, char const *bdir, char const *file) {
    using namespace std;
    char fa[2048], fb[2048];
    mkdirs(bdir);
    sprintf(fa,"%s/%s",adir,file);
    sprintf(fb,"%s/%s",bdir,file);
    struct stat s1;
    int ret = stat( fa, &s1);
    assert(ret==0);
    assert(S_ISREG(s1.st_mode));
    if(s1.st_nlink==1) {
        cout << "we link " << fa << " to " << fb << "\n";
        assert(link(fa,fb)==0);
    }else if(s1.st_nlink==2) {
        cout << "no need to link " << fa << " to " << fb << "\n";
        assert(isregfile(fb));
    }else cout << "this is wrong\n";
}

std::string get_json(const char *path, char const *hmm1, char const *hmm2) {
// std::string get_json(char const *hmm1, char const *hmm2) {

    using namespace rapidjson;
    using namespace std;

    /* char *path1 = strdup(hmm1), *path2 = strdup(hmm2);
    dirname(path1),dirname(path2);
    assert(strcmp(path1,path2)==0);
    free(path2); */

    struct stat s1;
    stat( (string(path)+"/"+hmm1).data() ,&s1);

    struct stat s2;
    stat( (string(path)+"/"+hmm2).data() ,&s2);

    /* cout << "\t" << __FILE__ << ":" << __LINE__ << " (" << __func__ << ") : " << "file: " << hmm1 << ", mode1: "<<s1.st_mode<< "\n";
    cout << "\t" << __FILE__ << ":" << __LINE__ << " (" << __func__ << ") : " << "file: " << hmm1 << ", nlink1: "<<s1.st_nlink<< "\n";
    cout << "\t" << __FILE__ << ":" << __LINE__ << " (" << __func__ << ") : " << "file: " << hmm2 << ", mode2: "<<s2.st_mode<< "\n";
    cout << "\t" << __FILE__ << ":" << __LINE__ << " (" << __func__ << ") : " << "file: " << hmm2 << ", nlink2: "<<s2.st_nlink<< "\n"; */

    assert(S_ISREG(s1.st_mode));
    assert(S_ISREG(s2.st_mode));

    /* assert(s1.st_nlink==2);
    assert(s2.st_nlink==2); */

    // cout << "using " << hmm1 << " and " << hmm2 << "\n\n";
    
    char j[16*1024];
    sprintf(
      j," { \"fastq\" : { \"path\" : { \"scratch\" : \"%s\" }, "
      " \"type\" : \"pe\", "
      " \"r1\" : { \"basename\" : \"%s\", \"size\" : %llu, \"modification\" : %llu}, "
      " \"r2\" : { \"basename\" : \"%s\", \"size\" : %llu, \"modification\" : %llu} "
      " } }", 
      // this is wrong but getting things running with -pedantic
      path, hmm1, (long long)s1.st_size, (long long)s1.st_mtime, hmm2, (long long)s2.st_size, (long long)s2.st_mtime
    );

    Document doc; 
    char buffer[sizeof(j)];
    memcpy(buffer, j, sizeof(j));
    if (doc.ParseInsitu(buffer).HasParseError()) cout << "there's an issue with the json\n", exit(1);

    // static const char* kTypeNames[] = { "Null", "False", "True", "Object", "Array", "String", "Number" };

    // printf("Original JSON:\n %s\n", j);
    StringBuffer sb2;
    PrettyWriter<StringBuffer> writer2(sb2);
    doc.Accept(writer2); // Accept() traverses the DOM and generates Handler events.
    // puts(sb2.GetString());
    // return 
    string bored = sb2.GetString(), bored2;
    for (unsigned int i=0;i<bored.length();++i){
        if(bored[i]=='\n') {
        }else if(bored[i]=='"') {
            bored2+='\\';
            bored2+='"';
        }else bored2+=bored[i];

    }
    return bored2;
}


// json.data(), SE["sample_internal_name"].data(), SE["lanenum"].data(), fcid.data() ); se_archive_dir

static unsigned int restart_pos = 0; // static so value initialised anyway?!?
static int restart_int[] = { 10, 60, 300, 900, 1800 };

// not really required anymore - but also wasn't actually re-starting anyway
void sig_action_function(int, siginfo_t *info, void*) { 
    std::cout << "\n---\nMessage received from child : '" << (char*)info->si_value.sival_ptr<< "'\n---\n" << std::endl; 
    if (strcmp((char*)info->si_value.sival_ptr,"All done")==0) 
    std::cout << "\n---\nParent fork shutting down cos we're done\n",exit(0);
    else std::cout << "Parent fork shutting down\n",exit(1);
} 

typedef struct _sig_ucontext { unsigned long     uc_flags; struct ucontext   *uc_link; stack_t           uc_stack; struct sigcontext uc_mcontext; sigset_t          uc_sigmask; } sig_ucontext_t;

void crit_err_hdlr(int sig_num, siginfo_t * info, void * ucontext) {
    void *             array[50]; void *             caller_address; char **            messages; int                size, i; sig_ucontext_t *   uc;
    uc = (sig_ucontext_t *)ucontext;
    caller_address = (void *) uc->uc_mcontext.rip; // RIP: x86_64 specific
    fprintf(stderr, "signal %d (%s), address is %p from %p\n", 
    sig_num, strsignal(sig_num), info->si_addr, 
    (void *)caller_address);
    size = backtrace(array, 50);
    array[1] = caller_address;
    messages = backtrace_symbols(array, size);
    for (i = 1; i < size && messages != NULL; ++i) { fprintf(stderr, "[bt]: (%d) %s\n", i, messages[i]); }
    free(messages);
    exit(1);
}


namespace config {

char const * conf = "#================================================================================\n\
# Dragen 2.5 Configuration File\n\
#================================================================================\n\
# SAMPLE SETUP\n\
intermediate-results-dir = /staging/tmp\n\
ref-dir = /staging/REF/b37_decoy/\n\
fastq-file1 = {{FQ1}}\n\
fastq-file2 = {{FQ2}}\n\
fastq-offset = 33 		# For CASAVA1.8 samples\n\
enable-auto-multifile = true\n\
#================================================================================\n\
# OUTPUT\n\
output-file-prefix = {{RGSM}}.{{EXPT_ID}}.{{RGPU}}\n\
output-format = BAM\n\
output-directory = {{SCRATCH_DIR}}/ # ALIGNMENT/BUILD37/DRAGEN/{{ST}}/{{RGSM}}/\n\
#================================================================================\n\
# READ GROUP INFO\n\
RGID = {{RGID}} # Read group ID\n\
RGLB = {{RGLB}} # Library\n\
RGPL = {{RGPL}} # Platform/technology\n\
RGSM = {{RGSM}} # Sample Name\n\
RGPU = {{RGPU}} # Platform unit\n\
RGCN = {{RGCN}} # Sequencing Center\n\
RGDT = {{RGDT}} # Date the run was produced.  Format:ISO_8601\n\
#================================================================================\n\
# ALIGNMENT\n\
enable-map-align-output = true\n\
enable-bam-indexing = true\n\
enable-sort = true\n\
enable-duplicate-marking = true\n\
remove-duplicates = false\n\
enable-sampling = true 			# automatically detect paired-end parameters with aligner test\n\
enable-deterministic-sort = true 	# ensure sort order is completely repeatable at cost of a small decrease in speed\n\
#================================================================================\n\
[Aligner]\n\
match-score = 1 	# Score increment for matching reference nucleotide\n\
mismatch-pen = 4 	# Score penalty for a mismatch\n\
gap-open-pen = 6 	# Score penalty for opening a gap (insertion or deletion)\n\
gap-ext-pen = 1 	# Score penalty for gap extension\n\
unclip-score = 5 	# Score bonus for reaching each edge of the read\n\
global = 0 		# If alignment is global (N-W) rather than local (S-W)\n\
pe-orientation = 0 	# Expected paired-end orientation: 0=FR, 1=RF, 2=FF\n\
pe-max-penalty = 60 	# Maximum pairing score penalty, for unpaired or distant ends\n\
mapq-max = 60 		# Ceiling on reported MAPQ\n\
supp-aligns = 3 	# Maximum supplimentary (chimeric) alignments to report per read\n\
sec-aligns = 0 		# Maximum secondary (suboptimal) alignments to report per read\n\
supp-as-sec = 0 	# If supplementary alignments should be reported with secondary flag\n\
hard-clips = 6 		# Flags for hard clipping: 0 primary, 1 supplementary, 2 secondary\n\
unpaired-pen = 80 	# Penalty for unpaired alignments in Phred scale\n\
dedup-min-qual = 15 	# Minimum base quality for calculating read quality metric for deduplication\n\
no-unpaired = 0 	# If only properly paired alignments should be reported for paired reads\n\
#================================================================================\n\
[Mapper]\n\
seed-density = 0.5 	# Requested density of seeds from reads queried in the hash table\n\
edit-mode = 0 		# 0 = No edits, 1 = Chain len test, 2 = Paired chain len test, 3 = Edit all std seeds\n\
edit-seed-num = 6 	# For edit-mode 1 or 2: Requested number of seeds per read to allow editing on\n\
edit-read-len = 100 	# For edit-mode 1 or 2: Read length in which to try edit-seed-num edited seeds\n\
edit-chain-limit = 29 	# For edit-mode 1 or 2: Maximum seed chain length in a read to qualify for seed editing\n\
map-orientations = 0 	# 0=Normal, 1=No Rev Comp, 2=No Forward  (paired end requires Normal)\n\
#================================================================================\n\
# FIN\n\
#================================================================================";

}



namespace rarp {
    typedef std::vector<std::string> LIST;
    typedef std::map<std::string,std::string> NLIST;
    typedef std::vector<NLIST> NLISTS;
    typedef std::vector<LIST> LISTS;
}

namespace opts {
    float wgs_min = 29.3, // 30.0, /// this is put into capture but whatever?!?
      wes_min = 60.0;
    bool commit = false, force = false;
    char const * serv = "seqprod";

    static bool email = true;

    using std::cout;
    using std::cerr;

    /// this MUST be removed but need to change users ASAP
    class MysqlUser {
      public:

        char const * user() const { return _user; }
        char const * pass() const { return _pass; }
        char const * host() const { return _host; }
        char const * connstr() const { return _connstr; }
        char const * connstr_quick_hack() const { return _connstr_quick_hack; }

        MysqlUser() {
            char const * tmp = getenv("LIMS_USER");
            // char const * tmp = getenv("USER");
            if(!tmp) cerr << "Need LIMS_USER env\n",exit(1);
            strcpy(_user,tmp);
            tmp = getenv("LIMS_PASS");
            if(!tmp) cerr << "Need LIMS_PASS env\n",exit(1);
            strcpy(_pass,tmp);
            tmp = getenv("LIMS_HOST");
            if(!tmp) cerr << "Need LIMS_HOST env\n",exit(1);
            strcpy(_host,tmp);

            sprintf(_connstr,"mysql -u%s -p%s -h%s sequenceDB",_user,_pass,_host); 
            sprintf(_connstr_quick_hack,"%s -BN -e ",_connstr);

            // if(strcmp(_user,"dh2880")==0 || strcmp(_user,"dsth")==0 )
              // cout << "we have " << _user << "\n" << /* "we have " << _pass << "\n" << */ "we have " << _host << "\n";

            //// should do a test connection here to check it all works??!!?!?

            // return a string for value semanatics or just have a data member?!?
            //char const * LAZY_CONN_STR_SEQDB = "mysql -udh2880 -p -hseqprod.igm.cumc.columbia.edu sequenceDB",
            // whatever, return initialised member?!? clearly, dbname should be treated consistently... 

        }

      private:
        char _user[32], _pass[32], _host[32], _connstr[1024], _connstr_quick_hack[1024];
    };

    static MysqlUser myuser;

    char const * usage = 
" run       This is the main pipeline manager. Run a single instance of this.\n"
" bcl       Run BCL conversion of a flowcell.\n"
" pipe      Perform a number of auxiliary functions associated with FASTQ and Pipeline integrity. Run several instances of this.\n"
" load_atav Self-explanatory. Only ever run a single instance of this.\n";
}

// was done in a major hurry. clearly change innerds to use mysql C API?!?
namespace db {

    char const // * LAZY_CONN_STR_SEQDB = "mysql -udh2880 -p -hseqprod.igm.cumc.columbia.edu sequenceDB",
      * LAZY_CONN_STR_DRGDB = "", // mysql -udh2880 -p -hannodb06 WalDB",
      * LAZY_CONN_STR_PGM = "";// mysql -upipeline -p -h10.73.50.31 annodb_pgm" ;

    rarp::NLIST get_named_row(char const * w, char const * const q) {

        char const * d = strcmp(w,"seqdb")==0 ? opts::myuser.connstr() : strcmp(w,"drgdb")==0 ? LAZY_CONN_STR_DRGDB : strcmp(w,"pgmdb")==0 ? LAZY_CONN_STR_PGM : 0 ; assert(d);

        char tmp[2048]; sprintf(tmp,"%s -B -e \"%s\"",d,q); 

        Popen ns1(tmp,16*1024,"r"); 
        rarp::LIST hrow,row;
        char * z3=ns1.getline(); tokenise(hrow,z3,'\t');
        assert(hrow.size()>=1);
        z3=ns1.getline(); tokenise(row,z3,'\t');
        rarp::NLIST nrow;
        for(unsigned i=0; i<hrow.size(); ++i) {
            if(i==0 && hrow[i]=="" /* silly : */ && row[i]=="") continue;
            nrow[hrow[i]]=row[i];
        }
        return nrow;
    }

    template<typename A> inline rarp::NLIST get_named_row(char const * w, char const * const q, A a) {
        char preptq[16*1024]; 
        sprintf(preptq,q,a);
        return get_named_row(w,preptq); 
    }

    template<typename A, typename B> inline rarp::NLIST get_named_row(char const * w, char const * const q, A a, B b) {
        char preptq[16*1024]; 
        sprintf(preptq,q,a,b);
        return get_named_row(w,preptq); 
    }

    template<typename A, typename B, typename C> inline rarp::NLIST get_named_row(char const * w, char const * const q, A a, B b, C c) {
        char preptq[16*1024]; 
        sprintf(preptq,q,a,b,c);
        return get_named_row(w,preptq); 
    }

    template<typename A, typename B, typename C, typename D> inline rarp::NLIST get_named_row(char const * w, char const * const q, A a, B b, C c, D d) {
        char preptq[16*1024]; 
        sprintf(preptq,q,a,b,c,d); 
        return get_named_row(w,preptq); 
    }

    template<typename A, typename B, typename C, typename D, typename E> inline rarp::NLIST get_named_row(char const * w, char const * const q, A a, B b, C c, D d, E e) {
        char preptq[16*1024]; std::cout << "using5 " << w << ", "<< q <<", "<<a<<", "<<b <<", "<<c<< ", "<<d<<", "<<e<<"\n";
        sprintf(preptq,q,a,b,c,d,e); 
        return get_named_row(w,preptq); 
    }

    void get_named_rows(char const * w, rarp::NLISTS & nrows, char const * const q) {

        char const * d = strcmp(w,"seqdb")==0 ? opts::myuser.connstr() : strcmp(w,"drgdb")==0 ? LAZY_CONN_STR_DRGDB : strcmp(w,"pgmdb")==0 ? LAZY_CONN_STR_PGM : 0 ; assert(d);

        /// why the heck are we overflowing all of a sudden?!?
        char tmp[8*1024]; sprintf(tmp,"%s -B -e \"%s\"",d,q); 

        rarp::LIST hrow;
        Popen ns1(tmp,16*1024,"r"); char * z3=ns1.getline(); tokenise(hrow,z3,'\t');
        assert(hrow.size()>=1);
        while( *( z3=ns1.getline() ) != '\0') { 
            rarp::NLIST nrow; rarp::LIST row; tokenise(row,z3,'\t');
            for(unsigned i=0; i<hrow.size(); ++i) nrow[hrow[i]]=row[i];
            nrows.push_back(nrow);
        }
    }

    // use a variadic template?!?
    template<typename A> inline void get_named_rows(char const * w, rarp::NLISTS & nrows, char const * const q, A a) {
        char preptq[16*1024]; 
        sprintf(preptq,q,a);
        get_named_rows(w,nrows,preptq); 
    }

    template<typename A, typename B> inline void get_named_rows(char const * w, rarp::NLISTS & nrows, char const * const q, A a, B b) {
        char preptq[16*1024]; 
        sprintf(preptq,q,a,b);
        get_named_rows(w,nrows,preptq); 
    }

    template<typename A, typename B, typename C> inline void get_named_rows(char const * w, rarp::NLISTS & nrows, char const * const q, A a, B b, C c) {
        char preptq[16*1024]; 
        sprintf(preptq,q,a,b,c);
        get_named_rows(w,nrows,preptq); 
    }

    std::string get_core_query(std::string whater) { 
        return db::get_named_row("seqdb","select replace( concat( 'gaf_bin=', e.gafbin, '&irb_protocol__irb_protocol=', e.protocol, '&fund_code__fund_code=', e.fundcode, '&sub_project_name=', e.subproject ), ' ', '+' ) CORE_QUERY from SampleT s join Experiment e on s.sample_id=e.sample_id where id = %s",whater.data())["CORE_QUERY"]; 
    }

    std::string get_core_query_NEW(std::string whater) { 
                return db::get_named_row("seqdb","select concat( 'id=', e.subproject_id ) CORE_QUERY_NEW from SampleT s join Experiment e on s.sample_id=e.sample_id where id = %s",whater.data())["CORE_QUERY_NEW"]; 
    }


}

namespace checks {

    inline void check_rsync_checksum(std::string const & f, int m, long long s) {
        using namespace std;

        if(!isregfile(f.data())) cout << "file " << f << " is missing\n",exit(1);
        struct stat st;
        stat(f.data(),&st);

        if(st.st_mtime==m) { // cout << f << " modification matches ("<<m<<")\n";
        } else cout << f << " modification doesn't match ("<<m<<")\n",exit(1);

        if(st.st_size==s) { // cout << f << " size matches ("<<s<<")\n";
        } else cout << f << " size doesn't match ("<<s<<")\n",exit(1);
    }

}



namespace lists {

/*
#define SAMPLE_SCRATCH_DIR "{{SAMPLE_SCRATCH_DIR}}"
#define SAMPLE_ARCHIVE_DIR "{{SAMPLE_ARCHIVE_DIR}}"
#define SAMPLE_UNIQUE_NAME "{{SAMPLE_UNIQUE_NAME}}"
#define SAMPLE_NAME "{{SAMPLE_NAME}}"
*/

    #define FILL_IN_DIR(X,Y,Z,A) char X[(A)]; (Y).fill_in_name((Z),X,sizeof(X)); \
    { int j=0; \
    for (unsigned int i = 0; i < strlen((X))-1; ++i) { \
        while((X)[j]=='/'&&(X)[j+1]=='/') ++j; \
        (X)[i]=(X)[j++]; \
    } \
    if((X)[j-1]=='/') (X)[j-1]='\0'; \
    (X)[j]='\0'; }

    #define FILL_IN_DIR_2(R,Y,Z,A) { char X[(A)]; (Y).fill_in_name((Z),X,sizeof(X)); \
    int j=0; \
    for (unsigned int i = 0; i < strlen((X))-1; ++i) { \
        while((X)[j]=='/'&&(X)[j+1]=='/') ++j; \
        (X)[i]=(X)[j++]; \
    } \
    if((X)[j-1]=='/') (X)[j-1]='\0'; \
    (X)[j]='\0'; \
    (Y).add_name((R),X); }

    struct NAMES {

        // NAMES(char * a, char * b, char * c, char * d) : ssd(a), sad(

        void add_name(char const * a, char const * b) {
            names.insert(std::make_pair(a,b));
        }
        void show_names() {
            for(std::map<std::string,std::string>::iterator i = names.begin(); i!=names.end(); i++)
            std::cout << i->first << " = " << i->second << "\n";
        }

        void fill_in_name(char const * in, char * out, size_t l) {
            memset(out,0,l);
            int len = strlen(in);
            bool sub = false;
            char name[1024], * p = name, * outp=out;
            for (int i = 0; i < len; ++i) {
                if(in[i]=='{' && in[i+1]=='{') {
                    ++in;
                    sub = true;
                    memset(name,0,sizeof name);
                    p = name;
                } else if(in[i]=='}' && in[i+1]=='}') {
                    ++in;
                    if(names.count(name)==0) std::cout << "what the heck is " << name << "\n", exit(1);
                    for (unsigned i =0; i<names.find(name)->second.size(); ++i) {
                        *outp++=names.find(name)->second[i];
                    }
                    sub = false;
                } else if(sub) {
                    *p++=in[i];
                } else {// if( (in[i]!='{' && in[i]!='}') {
                    *outp++=in[i];
                }
            }
        }

        std::map<std::string,std::string> names;

    };

}

bool bam_check(char const * bamf, rarp::NLIST & x) {

    using std::cout;

    cout << "checking bam " << bamf << "\n";
    if(!isregfile(bamf)) return false;

    char cmd[2048];
    int blarp=250;
    sprintf(cmd,"samtools view -h %s hs37d5 2>&1 | head -n%d",bamf,blarp);

    // cout << "ARGH\n"<<cmd<<"\n";

    // Popen ps(cmd,16*1024,"r");
    Popen ps(cmd,32*1024,"r");
    int lc=0, hc=0;
    char * z; // , blah[2048];
    bool fr = false;
    std::vector<std::string> ars2;
    while( *(z=ps.getline())!='\0') {
        ++lc;
        if(strstr(z,"EOF marker is absent")!=0) {
            cout << "eof issue.\n";
            return false;
        }
        if(*z=='@'){
            ++hc;
            if(memcmp(z,"@RG",3)==0) { 
                std::cout << "got '" << z+7 << "'\n";
                std::cout << "wrt " << x["PU"] << "\n";
                if(memcmp(x["PU"].data(),z+7,x["PU"].length())!=0) {
                    cout << "rg tag issue.\n";
                    return false;
                }
            }
        }else if(!fr++) tokenise(ars2,z,'\t');
            // assert(memcmp(I
        // else strcpy(blah,z);
        // if(*z!='@') std::cout << "got '" << z << "'\n";
    }
    std::cout << "overall="<<lc<<"\n";
    if(ars2.size()==0){
        cout << "there were NO decoy reads\n";
        return false;
    }

    // cout << "using " << ars2[0] << "\n";
    std::string tmp = ars2[0];
    ars2.clear();
    tokenise(ars2,tmp,':');
    // std::cout << "should use function but for now " << ars2[2]  << "\n";
    if(x["PU"]!=ars2[2]+"."+ars2[3]) {
        cout << "readname issue\n";
        return false;
    }
    if(lc!=blarp) {
        cout << "warning this seems pretty poor coverage - decoy reads = " << (lc-hc) << "\n";
        if( (lc-hc) < 1 ) {
            cout << "why aren't there reads for decoy.\n";
            return false;
        }
    }

    return true;

}

void align(int /* argc */, char ** /* argv */) {

    using namespace std;
    using namespace rarp;

    char hostname[1024];
    gethostname(hostname,sizeof(hostname));
    cout << "running on " << hostname << "\n";
    char const * mping, * ering;
    if(strstr(hostname,"dragen1")!=0) mping="fastq_mapping_d1", ering="mapping_error_d1";
    else if(strstr(hostname,"dragen2")!=0) mping="fastq_mapping_d2", ering="mapping_error_d2";
    else if(strstr(hostname,"dragen4")!=0) mping="fastq_mapping_d4", ering="mapping_error_d4";
    else cout << "run me from a dragen machine\n", exit(1);

    { char arsv2[8*1024];
    // awful, needs sooo mnuch cleaning and factoring!?!?
    sprintf(arsv2,"update Lane set rg_status = '%s' where rg_status = '%s'; select row_count() updated ", ering,mping);
    NLIST u = db::get_named_row("seqdb",arsv2); 
    if(u["updated"]=="1") cout << "marked previous error as " << mping << "\n";
    else if(u["updated"]=="0") cout << "nothing to mark\n";//  as " << mping << "\n";
    else assert(0); }

    NLISTS fc;
    db::get_named_rows("seqdb",fc,"select concat(f.fcillumid,'.',l.lanenum) PU, sample_internal_name, exomekit capture_kit, sample_type, is_released, "
    "f.fcillumid,l.fcid, "
    " p.adapter_id, a.sequence adapterlet, " // " p.adapterlet, "
    " l.lanenum,f.machine,f.machinetype,f.chemver,rg_status, experiment_id, p_prepid pseudo_prepid, l.prepid, "
    "data, upper(sample_type) ST, libkit, DATE_FORMAT( dateread1, '%Y-%m-%dT%TZ' ) as dateread1 "
    "from Lane l join Flowcell f on l.fcid=f.fcid join prepT p on l.prepid=p.prepid " 
    " join adapter a on a.id=p.adapter_id "
    " where "
    " rg_status = 'fastq_ready' "
    " order by prepid "
    "limit 1");
        
        if(fc.size()==0) {
            exit(0);
        }

        assert(fc.size()==1);
        NLIST & x = fc[0];
        cout << "\n\n----------------\nwe will process " << fc[0]["data"] << "\n";
        for(NLIST::iterator it = x.begin(); it != x.end(); it++ ) cout << "["<<it->first<<"] "<<it->second<<"\n";

        rapidjson::Document doc; 

        char buffer[16*1024]; 
        memset(buffer,0,sizeof(buffer));
        memcpy(buffer, x["data"].data(), x["data"].length());
        assert(x["data"]!="NULL");
        cout << "parsing " << x["data"] << "\n";
        if (doc.Parse(buffer).HasParseError()) {
            cout << "there's an issue with the json\n";
            exit(1);
        }

        assert(doc.FindMember("fastq")->value.FindMember("path")->value.HasMember("scratch"));
        assert(doc.FindMember("fastq")->value.HasMember("type"));
        assert(doc.FindMember("fastq")->value.FindMember("r1")->value.HasMember("size"));
        assert(doc.FindMember("fastq")->value.FindMember("r2")->value.HasMember("size"));

        // assert(doc.FindMember("fastq")->value.FindMember("files")->value.HasMember("type"));
        cout << "checking " << doc.FindMember("fastq")->value.FindMember("type")->value.GetString() << "\n";
        assert(strcmp(doc.FindMember("fastq")->value.FindMember("type")->value.GetString(),"pe")==0);

        string fp_s = doc.FindMember("fastq")->value.FindMember("path")->value.FindMember("scratch")->value.GetString(),
          fq1_s = fp_s + "/" + doc.FindMember("fastq")->value.FindMember("r1")->value.FindMember("basename")->value.GetString(),
          fq2_s = fp_s + "/" + doc.FindMember("fastq")->value.FindMember("r2")->value.FindMember("basename")->value.GetString();

        int fq1m = doc.FindMember("fastq")->value.FindMember("r1")->value.FindMember("modification")->value.GetInt(),
          fq2m =   doc.FindMember("fastq")->value.FindMember("r2")->value.FindMember("modification")->value.GetInt();
        long long fq1s =   doc.FindMember("fastq")->value.FindMember("r1")->value.FindMember("size")->value.GetInt64(),
          fq2s =   doc.FindMember("fastq")->value.FindMember("r2")->value.FindMember("size")->value.GetInt64();

        // cout << "we run\n\nfq1= " << fq1_s << " ("<<fq1s <<"/"<<fq1s<<")\nfq2= " << fq2_s << " ("<<fq2m<<"/"<<fq2s<<")\n\n";
        checks::check_rsync_checksum(fq1_s,fq1m,fq1s);
        checks::check_rsync_checksum(fq2_s,fq2m,fq2s);

        string ckpnt = fp_s + "/" + x["sample_internal_name"]+"."+x["experiment_id"] + "." + x["PU"] + ".done", 
          bam = x["sample_internal_name"]+"."+x["experiment_id"]+"."+x["PU"]+".bam",
          bam_f = fp_s + "/" + bam;

        // cout << "running with " << ckpnt << "\n";

        {

            cout << "\n> running dragen\n\n";

            { char arsv2[8*1024];
            sprintf(arsv2,"update Lane set rg_status = '%s' " 
              " where rg_status = 'fastq_ready' "
              " and prepid = %s and lanenum = %s and fcid = '%s' ; select row_count() updated",
              mping,fc[0]["prepid"].data(),fc[0]["lanenum"].data(),fc[0]["fcid"].data()
            );
            NLIST u = db::get_named_row("seqdb",arsv2); assert(u["updated"]=="1"); }

            string cmd, ext_conf = fp_s+"/"+x["PU"]+".dragen.conf";// , bam;
            if(fc[0]["sample_type"]=="Exome"||fc[0]["sample_type"]=="Genome") {

                lists::NAMES n;
                n.add_name("RGSM",x["sample_internal_name"].data());
                n.add_name("RGID",/* whatever!?! */ x["PU"].data());
                n.add_name("FQ1",fq1_s.data());
                n.add_name("FQ2",fq2_s.data());
                n.add_name("EXPT_ID",x["experiment_id"].data());
                n.add_name("RGLB",( /* add in libkit without spaces?!? */ "EXPT_ID:"+x["experiment_id"]+";PREP_ID:"+x["prepid"]).data());
                n.add_name("RGPU",x["PU"].data());
                n.add_name("ST",x["ST"].data());
                n.add_name("SCRATCH_DIR",fp_s.data());
                // n.add_name("SCRATCH_DIR","/nfs/seqscratch_ssd/");
                n.add_name("RGPL","ILLUMINA");
                n.add_name("RGDT",x["dateread1"].data());
                n.add_name("RGCN","IGM");
                n.add_name("OUT_PFX",(x["sample_internal_name"]+"."+x["experiment_id"]+"."+x["PU"]).data());
                // n.add_name("OUT_DIR",(x["sample_name"]+"."+x["experiment_id"]+"."+x["PU"]).data());
                char cf[2*strlen(config::conf)];

                // cout << "WT?!? " << x["pseudo_prepid"] << " vs " << x["experiment_id"] << "\n";

                assert(x["pseudo_prepid"]==x["experiment_id"]);
                n.fill_in_name(config::conf,cf,sizeof(cf));
                ofstream fh(ext_conf.data());

                if(!fh) cout << "problem writting conf file " << ext_conf << "\n",exit(1);
                fh << cf;
                fh.close();
                cmd=" -v -c " + ext_conf + " ";

            }else if(fc[0]["sample_type"]=="RNAseq"){

                cmd=" --enable-rna true --enable-rna-gene-fusion true -a /nfs/goldsteindata/refDB/mouse/GRCm38.87/Mus_musculus.GRCm38.87.gtf --enable-map-align true --enable-sort true "
                "--intermediate-results-dir /staging/tmp "
                "-c /nfs/seqscratch_ssd/dsth/dragen-user-defaults_2c1de82e74231199af88dbda224d3e93.cfg "
                "-r /nfs/seqscratch_ssd/dsth/RNA_Enabled_RefHash/GRCm38 "
                "--output-dir " + fp_s + 
                " --output-file-prefix " + x["sample_internal_name"]+"."+x["experiment_id"]+"."+x["PU"]+ // "_TMP " +
                " -1 " + fq1_s + " -2 " + fq2_s;

            }else{ assert(0); }

            cmd = "time /opt/edico/bin/dragen -f " + cmd + " --watchdog-active-timeout 600 2>&1 | tee " + fp_s + "/" + x["PU"] + ".dragen.combined.log";
            // cout << "has dragen cmd:\n" << cmd << "\n";

            if(system(cmd.data())!=0 || !bam_check(bam_f.data(),x)) {

                char arsv2[8*1024];
                // cout << "we have " << fc[0]["data"] << "\n";
                if(fq1s<100000||fq2s<100000) { // if(fq1s<75000||fq2s<75000) { 

                    sprintf(arsv2,"update Lane set rg_status = 'rg_error', rg_metrics_status = 'low_yield', failr1 = 1, failr2 = 1 where prepid = %s and lanenum = %s and fcid = '%s' ; select row_count() updated",
                    fc[0]["prepid"].data(),fc[0]["lanenum"].data(),fc[0]["fcid"].data());
                    cout << "error rg out - will get blocked and it's rather suspect at best?!? : " << arsv2 << "\n";
                    NLIST u = db::get_named_row("seqdb",arsv2); assert(u["updated"]=="1"); 

                }else{

                    sprintf(arsv2,"update Lane set rg_status = '%s' where prepid = %s and lanenum = %s and fcid = '%s' ; select row_count() updated",
                    ering,fc[0]["prepid"].data(),fc[0]["lanenum"].data(),fc[0]["fcid"].data());
                    cout << "error this out : " << arsv2 << "\n";
                    NLIST u = db::get_named_row("seqdb",arsv2); assert(u["updated"]=="1"); 
                }

                // don't let it cascade down?!?
                exit(1);

            }

            struct stat s;
            stat(bam_f.data(),&s);
            cout << "bam mtime="<<s.st_mtime << "\n";
            cout << "bam size="<<s.st_size << "\n";
            if(isregfile(bam_f.data())) touchfile(ckpnt.data());

            { char j[16*1024];
            sprintf(j,"update Lane set data = json_remove(data, '$.bam') where fcid = '%s' and lanenum = %s and prepid = %s ; select row_count() updated ",
              x["fcid"].data(), x["lanenum"].data(), x["prepid"].data() );
            cout << "have:\n" << j << "\n";
            NLIST u = db::get_named_row("seqdb",j); }
            // for(NLIST::iterator it = u.begin(); it != u.end(); it++ ) cout << "\t["<<it->first<<"] "<<it->second<<"\n\n";
            
            { char j[16*1024];
            sprintf(j,"update Lane set data = json_merge(data, '{ \\\"bam\\\" : { \\\"path\\\" : { \\\"scratch\\\" : \\\"%s\\\" }, \\\"basename\\\" : \\\"%s\\\", \\\"size\\\" : %llu, \\\"modification\\\" : %llu} }') "
              "where fcid = '%s' and lanenum = %s and prepid = %s ; select row_count() updated", 
              fp_s.data(), bam.data(), (long long)s.st_size, (long long)s.st_mtime, x["fcid"].data(), x["lanenum"].data(), x["prepid"].data() );
            // cout << "have:\n" << j << "\n";
            NLIST u = db::get_named_row("seqdb",j); 
            // cout << "WE GOT\n";
            // for(NLIST::iterator it = u.begin(); it != u.end(); it++ ) cout << "\t["<<it->first<<"] "<<it->second<<"\n\n";
            assert(u["updated"]=="1"); }

        }

        { char buffer[16*1024]; 
        memset(buffer,0,sizeof(buffer));
        assert(isregfile(bam_f.data()));
        sprintf(buffer,"select data->'$.bam.basename' bam, data->'$.fastq.path.scratch' scratch, data->'$.fastq.path' fastq, data->'$.bam.size' bs, data->'$.bam.modification' bm from Lane "
          "where fcid = '%s' and lanenum = %s and prepid = %s ", x["fcid"].data(), x["lanenum"].data(), x["prepid"].data());
        NLIST u = db::get_named_row("seqdb",buffer); 
        cout << "have:\n" << buffer << "\n\n";
        // for(NLIST::iterator it = u.begin(); it != u.end(); it++ ) cout << "\t["<<it->first<<"] "<<it->second<<"\n";
        string bc;
        { string bct = u["scratch"]+"/"+u["bam"]; for (unsigned int y=0;y<bct.length();++y) if(bct[y]!='"') bc+=bct[y]; cout << "checking " << bct << "\n"; }
        // cout << "checking " << bc << "\n";
        assert(isregfile(bc.data()));

        { char whater[1024];
        sprintf(whater,"update Lane set rg_status = 'fastq_mapped' where rg_status = '%s' and fcid = '%s' and lanenum = %s and prepid = %s ; "
          "select row_count() updated ",
          mping,x["fcid"].data(), x["lanenum"].data(), x["prepid"].data());
        for(NLIST::iterator it = x.begin(); it != x.end(); it++ ) cout << "["<<it->first<<"] "<<it->second<<"\n";
        // cout << whater << "\n";
        u = db::get_named_row("seqdb",whater); assert(u["updated"]=="1"); }

        }

    return;
}

rarp::NLIST get_an_se_by_status(char const * y) {
    //// why?!?
    rarp::NLISTS fc;
    db::get_named_rows("seqdb",fc,"select concat(f.fcillumid,'.',l.lanenum) PU, sample_internal_name, exomekit capture_kit, sample_type, is_released, "
        "f.fcillumid,l.fcid, "
          " p.adapter_id, a.sequence adapterlet, " // " p.adapterlet, "
        " l.lanenum,f.machine,f.machinetype,f.chemver,rg_status, experiment_id, p_prepid pseudo_prepid, l.prepid, "
        "data, upper(sample_type) seq_type_upper, libkit, dateread1 "
        "from Lane l join Flowcell f on l.fcid=f.fcid "
        " join prepT p on l.prepid=p.prepid "
          " join adapter a on a.id=p.adapter_id "
        " where "
        " rg_status = '%s' order by prepid limit 1",y);
    return fc[0];
}

void /* rarp::NLISTS */ get_se_group_by_uniform_status(rarp::NLISTS & fc, char const * y) {

    // rarp::NLISTS fc;
    db::get_named_rows("seqdb",fc,"select l.id lane_id, concat(f.fcillumid,'.',l.lanenum) PU, sample_internal_name, exomekit capture_kit, sample_type, is_released, "
      "f.fcillumid,l.fcid, " 

          " p.adapter_id, a.sequence adapterlet, " // " p.adapterlet, "

        " data->'$.fastq' wtf_f , data->'$.bam' wtf_b, "
      " l.lanenum,f.machine,f.machinetype,f.chemver,rg_status, experiment_id, p_prepid pseudo_prepid, "
      "l.prepid, data, upper(sample_type) seq_type_upper, libkit, dateread1 from Lane l join Flowcell f on l.fcid=f.fcid "
      "join prepT p on l.prepid=p.prepid "
        " join adapter a on a.id=p.adapter_id "
      " where "
      " (l.fcid,l.prepid) = (select fcid,prepid from Lane group by prepid,fcid "
      "having count(rg_status not in ('%s') or null ) = 0 limit 1) ", y);

    // paranoid!?!
    for (unsigned int h=0;h<fc.size();++h) assert(fc[h]["rg_status"]==y);
    // return fc;
}


///// mkfifo - pipe markdups straight to flagstat - sambamba for both?!?

template<typename A> inline void do_void_thing(const char * const q, A a) {
    char preptq[16*1024], preptq2[16*1024]; 
    strcpy(preptq,opts::myuser.connstr_quick_hack());
    sprintf(preptq2,q,a);
    sprintf(preptq+strlen(opts::myuser.connstr_quick_hack()),"\"%s\"",preptq2);
    if(system(preptq)!=0) { std::cout << "what : " << preptq << "\n\n\n"; exit(1); } 
}

template<typename A, typename B> inline void do_void_thing(const char * const q, A a, B b) {
    char preptq[16*1024], preptq2[16*1024]; 
    strcpy(preptq,opts::myuser.connstr_quick_hack());
    sprintf(preptq2,q,a,b);
    sprintf(preptq+strlen(opts::myuser.connstr_quick_hack()),"\"%s\"",preptq2);
    if(system(preptq)!=0) { std::cout << "what : " << preptq << "\n\n\n"; exit(1); } 
}

template<typename A, typename B, typename C> inline void do_void_thing(const char * const q, A a, B b, C c) {
    char preptq[16*1024], preptq2[16*1024]; 
    strcpy(preptq,opts::myuser.connstr_quick_hack());
    sprintf(preptq2,q,a,b, c);
    sprintf(preptq+strlen(opts::myuser.connstr_quick_hack()),"\"%s\"",preptq2);
    if(system(preptq)!=0) { std::cout << "what : " << preptq << "\n\n\n"; exit(1); } 
}

/* template<typename A> */ inline void update_status(rarp::NLIST /* using operator[] : const */ & dsm, int b,std::string & e, std::string const & s) {
    char preptq[16*1024]; 
    strcpy(preptq,opts::myuser.connstr_quick_hack());
    sprintf(preptq+strlen(opts::myuser.connstr_quick_hack()),"\"update dragen_sample_metadata set is_merged = %d where pseudo_prepid = %s\"",b,dsm["pseudo_prepid"].data());
    std::stringstream whater;
    whater << "BLOCKING : " << dsm["sample_type"] << " : " << dsm["sample_name"] << " : " << dsm["pseudo_prepid"] << " : " << dsm["capture_kit"] << " : ";
    whater << "\nupdate dragen_sample_metadata set is_merged = " << b << " where pseudo_prepid = " << dsm["pseudo_prepid"] << "\n" << s << "\n";
    std::cout << whater.str();//<< b << "\n" <<s<< "\n\n";A
    e+=whater.str();
    fflush(stdout);
    std::cout << "RUNNING UPDATE\n";
    if(system(preptq)!=0) { std::cout << "what : " << preptq << "\n\n\n"; exit(1); } 
}

enum BORED { sample_name=0, seqtype, capture_kit, pseudo_prepid, seqscratch_drive, BAMS, merged};

// dynamic modules e.g. Data::Dumper...
EXTERN_C void xs_init (pTHX);
EXTERN_C void boot_DynaLoader (pTHX_ CV* cv);
EXTERN_C void xs_init(pTHX) {
	const char *file = __FILE__;
	dXSUB_SYS;
	newXS("DynaLoader::boot_DynaLoader", boot_DynaLoader, file);
}


int main(int argc, char **argv){

    using std::cout;
    using std::cerr;
    using std::string;

    // if(strcmp(getenv("USER"),"pipe")!=0) cout << "run me as pipe\n",exit(1);

    { char hostname[1024];
    gethostname(hostname,1024);
    cout << "running on " << hostname << " @ " << time(0) << "\n";
    if(strstr(hostname,"atav")==0 && strstr(hostname,"dragen")==0) {
    // if(memcmp(hostname,"atav",4)!=0 && memcmp(hostname,"dragen",6)!=0) {
        cout << "run me from atav/dragen\n";
        exit(1);
    } }

    if(system("which samtools >/dev/null 2>&1")!=0) cout << "samtools must be in PATH\n",exit(1);

    if(argc<2) cerr << "usage: " << *argv << " <mode>\n\n" << opts::usage, exit(1);

    //// argh : export PATH=$PATH:/nfs/goldstein/software/samtools/samtools-1.5
    
    // {char argh[1024]; gethostname(argh,sizeof(argh)); cout << ">>>>>>>>>> " << argh << " : " << *argv << ":" << *(argv+1) << " : " << time(0) << "\n"; }

    /*
    if( ( argc==2 && strcmp(1[argv],"s3")==0 / * strstr(1[argv],"s3")==0 * / ) ) {
        // || ( argc==3  && strcmp(1[argv],"retry" ) ) ) {
        // bcl_backup();
        return 0;
    }else if ( argc==2 && strcmp(1[argv],"clear")==0 ) {
        // bcl_wipe();
        return 0;
    }else if ( argc==2 && strcmp(1[argv],"weekly")==0 ) {
        lims_driven_summary();
        return 0;
    }else if ( argc>=2 && strcmp(1[argv],"diag")==0 ) {
        // diagnostic_projects(--argc,++argv);
        return 0;
    }
    */

    if(argc==2) {

        if(strcmp(*(argv+1),"align")==0){

            cout << "Align wrapper.\n";
            while(1) align(argc,argv);

        }
    }
   

    // }else cerr << "unknown mode '" << *(argv+1) << "'\n",exit(1);

    cout << "bye\n";

    return 0;


}

