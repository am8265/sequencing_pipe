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
#include <tinyxml2.h>
#include <myhtml/tree.h>
#include <signal.h>
#include <execinfo.h>
// #ifdef _OPENMP // when was this used?
#include <omp.h>
// #endif
#undef min
#undef max
#include "interop/model/metric_base/metric_set.h"
#include "interop/model/metrics/corrected_intensity_metric.h"
#include "interop/io/metric_file_stream.h"
#include "interop/util/statistics.h"
#include "interop/logic/summary/run_summary.h"
#include "interop/model/metric_base/metric_set.h"
#include "interop/model/metrics/tile_metric.h"
#include "interop/model/run_metrics.h"
// #ifndef RELEASE_H
// #define RELEASE_H

// use vector of pairs of macro set funcname/func-pointer and shuffle
// put in those absurdly hacking merge recover things?!? 
// just add ethnicity here!
// add merge rescue object too?!?
// add reset to bcl BUT with lock issue telling them to run it ONLY if they're sure!?!? : update Flowcell set fc_status = 'registered' where FCillumID = 'HMJCJDSXX' and fc_status = 'sequenced' and fail = 0 ; select row_count() locked

#define BASE_DIR "/nfs/seqscratch09/informatics/"
#define SCRIPT_DIR BASE_DIR "logs/merge/"
#define LOG_DIR "/nfs/central/home/dh2880/.logs/"
#define ALIGNSTATS "/nfs/goldstein/software/alignstats/alignstats"
#define POSTMERGE "/nfs/seqscratch09/informatics/logs/postmerge/"

// int merge_and_release(int, char **);
// inline bool isregfile(const char* fn) { struct stat test; if (stat(fn, &test) != 0) { return false; }  return S_ISREG(test.st_mode); }                                   

namespace opts {
    float wgs_min = 29.3, // 30.0, /// this is put into capture but whatever?!?
      wes_min = 65.0;
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

bool myfunction (float i,float j) { return (i<j); }

void finish_with_error(MYSQL *con) { // http://zetcode.com/db/mysqlc/
    fprintf(stderr, "%s\n", mysql_error(con));
    mysql_close(con);
    exit(1);        
}

struct FUNKY { std::string experiment_id, prepid, sample_name, sample_type, capture_kit, priority, is_external, end_point; };

#define FDP_OFFSET 11
// #define FDP_OFFSET 15

char const * Q1 = "select "
  "SUM(LNYIELD) l_lane_sum, "
  "sum(rg_metrics_capturemean) l_capmean_sum, "
  "avg(rg_metrics_contamination) l_contamination_mean, "
  "avg(rg_metrics_dups) l_dups_mean, "
"replace(concat('gaf_bin=',e.gafbin,'&irb_protocol__irb_protocol=',e.protocol,'&fund_code__fund_code=',e.fundcode,'&sub_project_name=',e.subproject),' ','+') CORE_QUERY, "
"e.sample_type e_sample_type, "
" e.id e_experiment_id, s.sample_internal_name, e.sample_id, p.sample_type, p.exomekit, "
  "w.prep_id w_experiment_id, w.sample_finished w_sample_finished, w.sample_failure w_sample_failure, "
"merge_metrics_capturemean,merge_metrics_capturemedian,merge_metrics_contamination,merge_metrics_dups,eligible_readgroups, "
" pool_count, name pool_name, "
" (select group_concat(concat(step_name,':',pipeline_step_id,':',step_status)) from dragen_pipeline_step dsp join dragen_pipeline_step_desc junk on junk.id=dsp.pipeline_step_id where dsp.pseudo_prepid=e.id) dsp_step_ids, "
" (select count(pipeline_step_id) from dragen_pipeline_step where pseudo_prepid=e.id) dsp_step_id_count, "
" (select count(distinct(pp.prepid)) from prepT pp where pp.experiment_id=e.id) sum_prepids, "
" (select count(distinct(pp.status)) from prepT pp where pp.experiment_id=e.id) sum_statuses_count, "
" (select group_concat(distinct(pp.status)) from prepT pp where pp.experiment_id=e.id) sum_statuses, "
" (select count(distinct(pp.poolid)) from prepT pp where pp.experiment_id=e.id) sum_pools, "
" (select count(distinct(ll.fcid)) from Lane ll where ll.prepid in (select ppp.prepid from prepT ppp where ppp.experiment_id=e.id)) sum_flowcells, "
" (select count(distinct(concat(rg_status,':',rg_metrics_status))) from Lane ll where ll.prepid in (select ppp.prepid from prepT ppp where ppp.experiment_id=e.id)) sum_rg_statuses_count, "
" ( select group_concat(distinct(concat(rg_status,':',rg_metrics_status))) "
" from Lane ll join Flowcell ff on ll.fcid=ff.fcid "
" where ff.fail=0 and ff.complete=1 and ff.pipelinecomplete=1 and ll.failr1 is null and ll.failr2 is null "
" and ll.prepid in (select ppp.prepid from prepT ppp where ppp.experiment_id=e.id) "
" ) sum_rg_statuses, "
"count(l.id) sum_eligible_readgroups, "
  "e.is_released  e_is_released, "
  " group_concat(concat(p.prepid,':',name,':',fcillumid,'.',lanenum),'---') summary, "
  " group_concat(rg_status) rg_statuses, " 
  " group_concat(rg_metrics_status) rg_metrics_statuses, group_concat(l.id order by l.id) NEW_RGS_STAMP, e.rgs CURRENT_RGS_STAMP, " 
  " dqm.experiment_id dqm_experiment_id, "
  " d.experiment_id dsm_experiment_id, d.is_merged d_ism "
  " from Lane l                      "
  " join    Flowcell f                  on l.fcid               =f.fcid "
  " join         prepT p                     on l.prepid                 =p.prepid "
  " join         Experiment e                   on p.experiment_id          =e.id"
  " join         SampleT s                   on s.sample_id          =e.sample_id "
  " left join    dragen_sample_metadata d    on d.pseudo_prepid      =p.experiment_id "
  " left join    dsth.sample w               on w.prep_id            =p.experiment_id "
  " left join    dragen_qc_metrics dqm       on dqm.pseudo_prepid    =p.experiment_id "
  " left join    pool on p.poolid=pool.id "
  // evil
  " where (failedprep = 0 or failedprep >= 100) "
  // don't care about 'complete' button : " and f.complete = 1 "
  " and pipelinecomplete = 1 "
  " and failr1 is null and failr2 is null "
  " and externaldata is null "
  " and e.sample_type in ('Genome', 'Exome') "
  " and e.is_released in ('not_released', 'release_rejected') "
  // 'could' impose flowcell and pool grouping here?!?
  " group by e.id order by e.sample_type, p.poolid; ";

char const * REPLACE_INTO_DPS = 
  "replace into dragen_pipeline_step "
  "(pseudo_prepid,  pipeline_step_id,   version,    submit_time,            finish_time,            times_ran,  step_status) values "
  "(%s,             1,                  '0.0.1',    CURRENT_TIMESTAMP(),    CURRENT_TIMESTAMP(),    1,          '%s');";

inline bool isregfile(const char* fn) { struct stat test; if (stat(fn, &test) != 0) { return false; }  return S_ISREG(test.st_mode); }
inline bool isdir(const char* fn) { struct stat test; if (stat(fn, &test) != 0) { return false; }  return S_ISDIR(test.st_mode); }
inline void touchfile(const std::string& file) { FILE* fp = fopen(file.data(), "ab+"); fclose(fp); }       

namespace rarp {
    typedef std::vector<std::string> LIST;
    typedef std::map<std::string,std::string> NLIST;
    typedef std::vector<NLIST> NLISTS;
    typedef std::vector<LIST> LISTS;
}

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

void get_named_table(rarp::NLISTS & nrows, char const * const q) {
    rarp::LIST hrow;
    Popen ns1(q,16*1024,"r"); char * z3=ns1.getline(); tokenise(hrow,z3,'\t');
    assert(hrow.size()>=1);
    while( *( z3=ns1.getline() ) != '\0') { 
        rarp::NLIST nrow; rarp::LIST row; tokenise(row,z3,'\t');
        for(unsigned i=0; i<hrow.size(); ++i) nrow[hrow[i]]=row[i];
        nrows.push_back(nrow);
    }
}

void get_table(rarp::LISTS & rows, char const * const q) {
    Popen ns1(q,16*1024,"r"); 
    char * z3;
    while( *( z3=ns1.getline() ) != '\0') { 
        rarp::LIST row; 
        tokenise(row,z3,'\t');
        rows.push_back(row);
    }
}

template<typename A> inline void get_table(rarp::LISTS & rows, char const * const q, A a) {
    char preptq[16*1024]; 
    sprintf(preptq,q,a);
    // std::cout << "get_table<1>= " << preptq << "\n";
    get_table(rows,preptq); 
}

template<typename A, typename B> inline void get_table(rarp::LISTS & rows, char const * const q, A a, B b) {
    char preptq[16*1024]; 
    sprintf(preptq,q,a,b);
    // std::cout << "get_table<2>= " << preptq << "\n";
    get_table(rows,preptq); 
}

template<typename A> inline void get_named_table(rarp::NLISTS & nrows, char const * const q, A a) {
    char preptq[16*1024]; 
    sprintf(preptq,q,a);
    get_named_table(nrows,preptq); 
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
            FCT.insert(make_pair("ALXX","HiSeqX (8-lane) flowcell"));
            FCT.insert(make_pair("CCXX","HiSeqX (8-lane) flowcell"));
            FCT.insert(make_pair("CCXY","HiSeqX (8-lane) flowcell"));
            FCT.insert(make_pair("AAXX","Genome Analyzer"));
            FCT.insert(make_pair("ADXX","HiSeq Rapid Run (2-lane) v1 (HiSeq 1500/2500)"));
            FCT.insert(make_pair("AGXX","High-Output NextSeq"));
            FCT.insert(make_pair("AMXX","HiSeq RR v2"));
            FCT.insert(make_pair("BCXX","HiSeq Rapid Run (2-lane) v1.5/v2 (HiSeq 1500/2500)"));
            FCT.insert(make_pair("BCXY","HiSeq Rapid Run (2-lane) v2 (HiSeq 1500/2500)"));
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

namespace lims {
    
    inline std::string get_uname(rarp::NLIST & entry) { return entry["dsm_sample_name"]+"."+entry["dsm_pseudo_prepid"]; }
    inline std::string get_archive(rarp::NLIST & entry) { return entry["qc_AlignSeqFileLoc"]+"/"+get_uname(entry)+"/"; }
    inline std::string get_scratch(rarp::NLIST & entry) { 
        std::string type = entry["dsm_sample_type"];
        for(unsigned i=0; i<type.length(); ++i) type[i]=::toupper(type[i]);                                                                      
        return "/nfs/"+entry["dsm_seqscratch_drive"]+"/ALIGNMENT/BUILD37/DRAGEN/"+type+"/"+get_uname(entry);
    }

}

namespace query {

    using namespace std;

    template <typename A> inline bool silly_update(char const * b, A a) {
        char q[2048];
        strcpy(q,b);
        strcat(q,"; select row_count() as affected");
        // std::cout << "using " << q << " and " << a << "\n";
        rarp::NLIST arsv2 = db::get_named_row("seqdb",q,a);
        // std::cout << "we modified " << arsv2["affected"] << " affected rows\n"; 
        return arsv2["affected"]=="1"; 
    }

}

struct Timey {

    Timey(struct tm * tm) : _tm(*tm) {
        _tm.tm_isdst=0; 
        // Timey(struct tm * tm) : _tm(memset(&_tm,tm,sizeof(struct tm))) {
        // printf("using struct tm\n");
        blah();
    }

    Timey(time_t t) : _tm(*localtime(&t)) /* use localtime or gmtime for this?!? */ {
        _tm.tm_isdst=0; 
        blah();
    }

    Timey(char const * a, char const * b,bool clean=false) {

        _tm.tm_isdst=0; 
        if(clean) {

            char silly[1024];
            strcpy(silly,a);

            for (unsigned p = 0, dop=0; p<strlen(silly); ++p ){           
                if(silly[p]=='.') ++dop;

                assert(dop<2);
                if(dop && silly[p]==' ') --dop;
                if(dop) silly[p]=' ';

            }
            a=silly;

        }

        strptime(a, b,&_tm);
        blah();
    }

    void blah() {
        _tm.tm_isdst=0; // not explicitly setting this, or zero-initialising gives non-determinstics behavirou w/-or-w/o 1h offset?!?
        assert(mktime(&_tm)>=1449263413); // ~end of 2015?!?
        assert(mktime(&_tm)<=1607116213); // ~end of 2020?!?
        std::cout.flush();
    }

    void check() { 
        // this seesm to force initialisation?!?
        assert(mktime(&_tm)==1512508213); 

    }

    // why?!?
    mutable char _buf[1024];
    struct tm _tm;

    // get re-write of value when invoking multiple times in cout... - i.e. undefined order?!?
    char const * epoch_time_as_string(char * const ext) const {
        memset(_buf,sizeof(_buf),0);
        strftime(_buf, sizeof(_buf), "%s", &_tm);
        strcpy(ext,_buf);
        return ext; // return _buf;
    }

    time_t epoch_time_as_time_t() { return mktime(&_tm); }

    char const * iso_time(char * const ext) const {
        memset(_buf,sizeof(_buf),0);
        strftime(_buf, sizeof(_buf), "%Y-%m-%dT%H:%M:%S", &_tm);
        strcpy(ext,_buf);
        return ext; // return _buf;
    }
};

namespace md5 {

    uint32_t grab_something(char const * const x) { // uint64_t grab_something(char const * const x) {
        unsigned char digest[MD5_DIGEST_LENGTH];
        MD5( (unsigned char*) x, strlen(x), (unsigned char*) &digest );  
        uint32_t tmp; // uint64_t tmp;
        memcpy(&tmp,&digest,sizeof(tmp));
        return tmp;
    }

    //https://stackoverflow.com/questions/7627723/how-to-create-a-md5-hash-of-a-string-in-c
    char *str2md5(const char *str, int length) {
        int n;
        MD5_CTX c;
        unsigned char digest[16];
        char *out = (char*)malloc(33);

        MD5_Init(&c);

        while (length > 0) {
            if (length > 512) {
                MD5_Update(&c, str, 512);
            } else {
                MD5_Update(&c, str, length);
            }
            length -= 512;
            str += 512;
        }

        MD5_Final(digest, &c);

        for (n = 0; n < 16; ++n) {
            snprintf(&(out[n*2]), 16*2, "%02x", (unsigned int)digest[n]);
        }

        return out;
    }

    /// clearly this should be called by str2md5
    /* error to return pointer to locally allocated mem!?! unsigned char * */ void str2md5bin(const char *str, int length,unsigned char * digest) {// , unsigned char &* digest) {
        // int n;
        MD5_CTX c;
        // unsigned char digest[16];
        // char *out = (char*)malloc(33);

        MD5_Init(&c);

        while (length > 0) {
            if (length > 512) {
                MD5_Update(&c, str, 512);
            } else {
                MD5_Update(&c, str, length);
            }
            length -= 512;
            str += 512;
        }

        MD5_Final(digest, &c);

        // return digest;
    }

}

namespace seq {

    static char const * ARCHIVE_DIR = "/nfs/archive/p2018/FASTQ/"; /* ,
      * NETAPP_OUT_DIR = "/nfs/seqscratch_ssd/",
      * SEQ_RUN_DIR = "/nfs/hts/novaseq/sequence/Runs/"; */

    #define CLEAN_NAME(X) \
    { int j=0; \
    for (unsigned int i = 0; i < strlen((X))-1; ++i) { \
        while((X)[j]=='/'&&(X)[j+1]=='/') ++j; \
        (X)[i]=(X)[j++]; \
    } \
    if((X)[j-1]=='/') (X)[j-1]='\0'; \
    (X)[j]='\0'; }

    std::string se_archive_dir(rarp::NLIST & SE) { 
        assert(SE.count("seq_type_upper"));
        assert(SE.count("sample_internal_name"));
        assert(SE.count("fcillumid"));
        std::string tmp = std::string(seq::ARCHIVE_DIR) + "/" + SE["seq_type_upper"] + "/" + SE["sample_internal_name"] + "/" + SE["fcillumid"] + "/";
        char Y[2048];
        strcpy(Y,tmp.data());
        CLEAN_NAME(Y);
        tmp=Y;
        return tmp;
    }

    #undef CLEAN_NAME
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

namespace pipey {

    void run_cmd(
      std::vector<std::string> const & cmds,
      std::string const & dir,
      std::string const & prefix
    ) {
        using namespace std;

        vector<string> run;
        for (vector<string>::const_iterator it = cmds.begin(); it != cmds.end(); it++) {
            string cf = dir + "/" + prefix  + "_" + md5::str2md5(it->data(),it->length());   
            cf+=".cp.txt";
            cout << "CKPT= " << cf << "\n";
            if(isregfile(cf.data())) {
                cout << "SKIPPING_COMPLETED= " << *it << "\n";
            }else{
                cout << "CMD= " << *it << "\n\n";
                // run.push_back(it->data());
                if(system(it->data())) cout << "there was a problem with " << *it << "\n\n",exit(1);
                touchfile(cf.data());
            }

            // how does this work without VLA?!?
            // std::thread t[run.size()];
            // for (int i = 0; i<run.size(); ++i) t[i]=std::thread(dummy, run[i], i);
            // for (int i = 0; i<run.size(); ++i) t[i].join();
        }

        // if(system(it->data())) cout << "there was a problem with " << *it << "\n\n",exit(1);
        // touchfile(cf.data());

        cout << "wait for series\n\n";
    }
}

#define SINLGE_CMD(N,X,Y,Z) (N).fill_in_name((X),(Y),sizeof((Y))); (Z).push_back((Y)); run_cmd((Z)); (Z).clear();
#define SINLGE_CMD2(N,X,Y,Z) { char b1[1024]; strcpy(b1,(Z)); SINLGE_CMD(N,b1,(X),(Y)) }

namespace core { 

    // static float min_wes, min_wgs;
    static float min_wes = opts::wes_min, min_wgs = opts::wgs_min;
    // static float min_wes = 75.0, min_wgs = 29.0;

    float get_min_wes() { return min_wes; }
    float get_min_wgs() { return min_wgs; }

    size_t process_response(char *ptr, size_t, size_t nmemb, void *userdata) {
    // size_t process_response(char *ptr, size_t size, size_t nmemb, void *userdata) {

        using std::cout;
        using std::string;
        // cout << "will copy " << size << " bytes of " << ptr << " in to " << userdata <<"\n";
        memcpy(userdata,ptr,nmemb);
        ///curl_easy_perform() failed: Failed writing received data to disk/application
        // cout << "WE GOT " << size << ", " << nmemb << ", "<<strlen((char*)userdata)<<"\n";
        return nmemb;
        /// return strlen((char*)userdata)+1;
        // cout << "will copy " << size << " bytes of " << ptr << " in to " << userdata <<"\n";
        // cout << "x="<<ptr <<"\n";
        // strcpy((char*)userdata,"3232");
        //
    }

class Core {

    bool _is_approved, _is_releasable;
    char * _useful_name;
    Core();
    float _min_wes, _min_wgs;

  public:

    // static float get_min_wes() { return min_wes; }
    // static float get_min_wgs() { return min_wgs; }

    Core(char const * query) : _useful_name(new char[2048]) { 

        using std::cout;
        using std::string;

        CURL *curl;
        CURLcode res;

        curl = curl_easy_init();
        if(!curl) cout << "had an issue\n", exit(1);

        char search[2048];
        memset(search,0,sizeof(search));
        strcpy(search,"https://core.igm.cumc.columbia.edu/core/api/subproject/?");
        strcat(search,query);

        //
        //Accept: application/json; indent=4
        //   curl_easy_setopt(curl, CURLOPT_POSTFIELDS, 
        #ifdef DEBUG 
        curl_easy_setopt(curl, CURLOPT_URL, search);
        struct curl_slist *list = NULL;
        list = curl_slist_append(list, "Accept: application/json; indent=4");
        curl_easy_setopt(curl, CURLOPT_HTTPHEADER, list);
        #else
        curl_easy_setopt(curl, CURLOPT_URL, search),
        #endif
        curl_easy_setopt(curl,  CURLOPT_SSL_VERIFYPEER, 0L), curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L), curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, process_response);
        // ///* example.com is redirected, so we tell libcurl to follow redirection */

        char rarp[4*2048];
        memset(rarp,0,sizeof(rarp));

        curl_easy_setopt(curl, CURLOPT_WRITEDATA, &rarp);

        res = curl_easy_perform(curl);
        if(res != CURLE_OK) fprintf(stderr, "curl_easy_perform() failed: %s\n",curl_easy_strerror(res));
        curl_easy_cleanup(curl);

        // cout << "rarp='"<<rarp<<"'\n";

        rapidjson::Document doc; 

        // char buffer[16*1024]; 
        // memset(buffer,0,sizeof(buffer));
        // memcpy(buffer, rgs[x]["data"].data(), rgs[x]["data"].length());
        
        #ifdef DEBUG
        // cout << getenv("USER") << "\n";
        cout << rarp << "\n";
        #endif

        if (doc.Parse(rarp).HasParseError()) cout << "there's an issue with the json\n", exit(1);
        // cout << "what1\n";
        // assert(doc.IsObject());
        assert(doc.IsArray());
        // assert(doc2.IsArray() && doc2.size()==1);
        assert(doc[0].IsObject());
        assert(doc[0].HasMember("id"));
        assert(doc[0].HasMember("approved"));

        int subproj=doc[0].FindMember("id")->value.GetInt();
        bool approved=doc[0].FindMember("approved")->value.GetBool();

        string project =doc[0].FindMember("project")->value.FindMember("project_name")->value.GetString(),
          sub_project_name =doc[0].FindMember("sub_project_name")->value.GetString();

        string sow="";
        if(!doc[0].FindMember("contract_uri")->value.IsNull()) sow=doc[0].FindMember("contract_uri")->value.GetString();

        /*
        memset(search,0,sizeof(search));
        strcpy(search,"https://core.igm.cumc.columbia.edu/core/api/project/?project_name=");
        strcat(search,project.data());
        // strcat(search,"Lymphatic+Anomaly");

        for (int y = 0 ; y < strlen(search); ++y ) if(search[y]==' ') search[y]='+';

        cout << "we use " << search << "\n";

        memset(rarp,0,sizeof(rarp));
        curl = curl_easy_init();

        curl_easy_setopt(curl, CURLOPT_URL, search),
        // curl_easy_setopt(curl, CURLOPT_URL, "https://core.igm.cumc.columbia.edu/core/api/project/?project_name=Lymphatic+Anomaly"),
        curl_easy_setopt(curl,  CURLOPT_SSL_VERIFYPEER, 0L), curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L); // /// * example.com is redirected, so we tell libcurl to follow redirection * /
        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, process_response);
        curl_easy_setopt(curl, CURLOPT_WRITEDATA, &rarp);
        res = curl_easy_perform(curl);
        if(res != CURLE_OK) fprintf(stderr, "curl_easy_perform() failed: %s\n",curl_easy_strerror(res));
        curl_easy_cleanup(curl);
        */

        assert(doc[0].FindMember("project")->value.IsObject()); // value.GetString();


        assert(doc[0].FindMember("project")->value.HasMember("id"));

        assert(doc[0].FindMember("project")->value.HasMember("project_name"));
        assert(doc[0].FindMember("project")->value.HasMember("principal_investigator"));
        /// only relevant to run_flowcell/metrics?!?
        assert(doc[0].FindMember("project")->value.HasMember("genome_max_contamination"));
        assert(doc[0].FindMember("project")->value.HasMember("exome_max_contamination"));

        //// just take arg of sample_type e.g. get_min("exome");
        assert(doc[0].FindMember("project")->value.HasMember("genome_target_min_coverage"));
        assert(doc[0].FindMember("project")->value.HasMember("exome_target_min_coverage"));
        assert(doc[0].FindMember("project")->value.HasMember("rnaseq_min_unique_fragments"));
        assert(doc[0].FindMember("project")->value.HasMember("project_type"));

        // cout << "sow="<< sow<< "\n";

        string who=doc[0].FindMember("project")->value.FindMember("principal_investigator")->value.GetString(),
        proj_type=doc[0].FindMember("project")->value.FindMember("project_type")->value.GetString();
        if(sow.empty()) sow="None";

        //cout << "\nsub_project_name="<< sub_project_name<< "("<<subproj<< ") of " << "project="<< project<< "(" 
        //<< doc[0].FindMember("project")->value.FindMember("id")->value.GetInt() << ") [type="<< proj_type << "/"<<sow <<"] (" << who << ")\n";

        //cout << "approved="<< approved<< "\n";

        _is_releasable=doc[0].FindMember("project")->value.FindMember("release_approved")->value.GetBool();
        _is_approved=approved;
        //cout << "_is_releasable= " << _is_releasable << "\n"
         // << "_is_approved= " << _is_approved << "\n";

        _min_wgs=
          ( doc[0].FindMember("project")->value.FindMember("genome_target_min_coverage")->value.IsNull() 
          ) ?
            // || doc[0].FindMember("project")->value.FindMember("genome_target_min_coverage")->value.GetFloat()==0.0 ) ?
          get_min_wgs():doc[0].FindMember("project")->value.FindMember("genome_target_min_coverage")->value.GetFloat();
        _min_wes=doc[0].FindMember("project")->value.FindMember("exome_target_min_coverage")->value.IsNull()?
          get_min_wes():doc[0].FindMember("project")->value.FindMember("exome_target_min_coverage")->value.GetFloat();

        //cout << "_min_wes= " << _min_wes << "\n"
        //  << "_min_wgs= " << _min_wgs << "\n";

        assert(doc[0].FindMember("project")->value.FindMember("rnaseq_min_unique_fragments")->value.IsNull());

        assert(doc[0].FindMember("project")->value.HasMember("genome_target_min_coverage"));
        assert(doc[0].FindMember("project")->value.HasMember("exome_target_min_coverage"));

        sprintf(_useful_name,"sub_project_name=%s(%d),project=%s(%d)",sub_project_name.data(),subproj,project.data(),doc[0].FindMember("project")->value.FindMember("id")->value.GetInt());

        // cout << "what2\n";
        // cout << "rarp='"<<rarp<<"'\n";

/*
        rapidjson::Document doc2; 
        if (doc2.Parse(rarp).HasParseError()) cout << "there's an issue with the json\n", exit(1);
        // cout << "what1\n";
        // assert(doc2.IsObject());
        assert(doc2.IsArray());// && doc2.size()==1);
        assert(doc2[0].IsObject());
        // assert(doc2[0].HasMember("genome_target_coverage"));
*/

        // assert(doc2[0].HasMember("exome_target_coverage"));
        // assert(doc2[0].HasMember("rnaseq_target_unique_fragments"));


    }

    ~Core() { delete[] _useful_name; }

    // float min_releasable(char const * why, int current) { /// tkae type and value or just type
    bool hits_coverage(std::string const & st,float cov) { 
        assert(st=="Genome"||st=="Exome");
        //std::cout << "we have " << st << " with cov= " << cov << " and wes_min= "<<_min_wes << " and wgs_min= " << _min_wgs << "\n";

        // if(cov>120) assert(0);

        if(st=="Genome") return cov>=_min_wgs;
        if(st=="Exome") return cov>=_min_wes;
        return false;
    }

    bool hits_legacy_coverage(std::string const & st,int /* yield */ ) { 
        assert(0);
        assert(st=="Genome"||st=="Exome");
        // if(st=="Genome") return (yield/4500) >=_min_wgs;
        // if(st=="Exome") return yield>cov>=_min_wes;
    }
    bool is_approved() const { return _is_approved; }
    bool is_releasable() const { return _is_releasable; }
    char const * useful_name() const { return _useful_name; }
};
     
    // static float min_wes = 75.0, min_wgs = 29.0;

}

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
namespace html {

    struct res_html {
        char  *html;
        size_t size;
    };

    struct res_html load_html_file(const char* filename) {
        FILE *fh = fopen(filename, "rb");
        if(fh == NULL) {
            fprintf(stderr, "Can't open html file: %s\n", filename);
            exit(EXIT_FAILURE);
        }
        if(fseek(fh, 0L, SEEK_END) != 0) {
            fprintf(stderr, "Can't set position (fseek) in file: %s\n", filename);
            exit(EXIT_FAILURE);
        }
        long size = ftell(fh);
        if(fseek(fh, 0L, SEEK_SET) != 0) {
            fprintf(stderr, "Can't set position (fseek) in file: %s\n", filename);
            exit(EXIT_FAILURE);
        }
        if(size <= 0) {
            fprintf(stderr, "Can't get file size or file is empty: %s\n", filename);
            exit(EXIT_FAILURE);
        }
        char *html = (char*)malloc(size + 1);
        if(html == NULL) {
            fprintf(stderr, "Can't allocate mem for html file: %s\n", filename);
            exit(EXIT_FAILURE);
        }
        size_t nread = fread(html, 1, size, fh);
        if ( (long)nread != size) {
            // fprintf(stderr, "could not read %ld bytes (" MyCORE_FMT_Z " bytes done)\n", size, nread);
            exit(EXIT_FAILURE);
        }
        fclose(fh);
        struct res_html res = {html, (size_t)size};
        return res;
    }

    void print_node_attr(myhtml_tree_node_t *node) {
        myhtml_tree_attr_t *attr = myhtml_node_attribute_first(node);
        printf("\n\t> atrribs :");
        while (attr) {
            const char *name = myhtml_attribute_key(attr, NULL);
            if(name) {
                printf("[a] %s", name);
                const char *value = myhtml_attribute_value(attr, NULL);
                if(value)
                    printf("[b]=\"%s\"", value);
            }
            attr = myhtml_attribute_next(attr);
        }
        puts("");
    }

    void print_tree2(myhtml_tree_t* tree, myhtml_tree_node_t *node, size_t inc,rarp::NLISTS & tbl) {
        using std::cout;
        while (node) {
            if(myhtml_node_tag_id(node)==MyHTML_TAG_H2){
                rarp::NLISTS out;
                char const * const ct = myhtml_node_text(node->child, NULL);
                // printf("i'm an h2 and my next's next should be a table = %s, ", myhtml_tag_name_by_id(tree, myhtml_node_tag_id(node->next->next), NULL) );
                // printf("and my text (via my child) is : '%s'\n", myhtml_node_text(node->child, NULL));
                rarp::LIST hdr;
                myhtml_tree_node_t *tr = node->next->next->child->next->child, *th=tr->child;
                assert(myhtml_node_tag_id(th->next)==MyHTML_TAG_TH); 
                while(th){
                    if(myhtml_node_tag_id(th)==MyHTML_TAG_TH) 
                    hdr.push_back(myhtml_node_text(th->child, NULL));
                    th = myhtml_node_next(th);
                }
                tr=tr->next->next;
                while(tr){
                    th=tr->child;
                    rarp::NLIST t;
                    int i=0;
                    while(th){

/*
cout << "this is a " << myhtml_node_tag_id(th) << " : th=" << MyHTML_TAG_TH << ", td="<< MyHTML_TAG_TD<<"\n";
cout << " we check the node has a child " << (th->child==0?"NULL":"it has a child") << "\n";
if(th->child!=0) {
cout << " child have text " << myhtml_node_text(th->child, NULL) << "\n";
}
fflush(stdout);
*/

                        // if(myhtml_node_tag_id(th)==MyHTML_TAG_TH || myhtml_node_tag_id(th)==MyHTML_TAG_TD) 
                        //   t[hdr[i++]]=myhtml_node_text(th->child, NULL);
                        if(myhtml_node_tag_id(th)==MyHTML_TAG_TH || myhtml_node_tag_id(th)==MyHTML_TAG_TD) {
                            if(th->child!=0) t[hdr[i++]]=myhtml_node_text(th->child, NULL);
                            else t[hdr[i++]]="NULL";
                        }
                        th = myhtml_node_next(th);
                    }
                    assert(t.size()==hdr.size());
                    out.push_back(t);
                    tr=tr->next->next; //tr = myhtml_node_next(tr);
                }
                cout << __FILE__ << ":" << __LINE__ << " ( " << __func__ << " ) TABLE " << ct << " is " << hdr.size() << " and " << out.size() << "\n";
                ///// seems we don't actually store 'Top Unknown Barcodes' even in terms of %?!?
                ///// thus 'should've' just used the comment remover and parsed out <table> </table> and then lines?!?
                if(strcmp(ct,"Lane Summary")==0) {
                    tbl=out;
                    return;
                /// already have this from interop?!?
                }else if(strcmp(ct,"Flowcell Summary")==0) {
                  for(rarp::NLIST::iterator it = out[0].begin(); it!=out[0].end(); it++){
                    cout << __FILE__ << ":" << __LINE__ << " ( " << __func__ << " ) ["<<it->first<<"] "<<it->second<<"\n";
                  }
                }
            } 
            print_tree2(tree, myhtml_node_child(node), (inc + 1),tbl); // print children
            node = myhtml_node_next(node);                         // puts("---- done with children ----");
        }
    }

    void /* int */ get_named_rows_from_fc_summary(rarp::NLISTS & tbl, char const * file) {
        struct res_html res = load_html_file(file);
        // basic init
        myhtml_t* myhtml = myhtml_create();
        myhtml_init(myhtml, MyHTML_OPTIONS_DEFAULT, 1, 0);
        // init tree
        myhtml_tree_t* tree = myhtml_tree_create();
        myhtml_tree_init(tree, myhtml);
        // pFLOWCELL_REPORT
        myhtml_parse(tree, MyENCODING_UTF_8, res.html, res.size);
        // print tree
        myhtml_tree_node_t *node = myhtml_tree_get_document(tree);
        print_tree2(tree, myhtml_node_child(node), 1,tbl);
        // release resources
        myhtml_tree_destroy(tree);
        myhtml_destroy(myhtml);
        free(res.html);
    }

}

/// namespace { // unnamed!?!

#define MIN_RGS 5

#define GET_FLOWCELL "select * from Flowcell where FCillumID = '%s' and fail = 0"

struct ARCHIVE { std::string data, sample_internal_name, lanenum, archive_dir, html, read1, read2; };

#define BCL_CMD "export LD_LIBRARY_PATH=/nfs/goldstein/software/bcl2fastq2-v2.20-x86_64/lib:/nfs/goldstein/software/libxml2-2.9.4-x86_64/lib:/nfs/goldstein/software/binutils-2.26.1/lib:/nfs/goldstein/software/boost_1_54_0/lib:/nfs/goldstein/software/gcc-4.9.3/lib64:$LD_LIBRARY_PATH ; \
  /nfs/goldstein/software/bcl2fastq2-v2.20-x86_64/bin/bcl2fastq \
  --runfolder-dir %s \
  --output-dir %s/%s \
  --sample-sheet %s \
  --barcode-mismatches %d"

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
        // cout << "we link " << fa << " to " << fb << "\n";
        assert(link(fa,fb)==0);
    }else if(s1.st_nlink==2) {
        // cout << "no need to link " << fa << " to " << fb << "\n";
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

int read_interop_file(
  const char* filename, 
  illumina::interop::model::metric_base::metric_set<illumina::interop::model::metrics::tile_metric>& tile_metric_set) {

    using namespace illumina::interop::model::metric_base;
    using namespace illumina::interop::model::metrics;
    using namespace illumina::interop::io;
    using namespace illumina::interop::util;
    using namespace illumina::interop;

    try { read_interop(filename, tile_metric_set); 
    }catch(const incomplete_file_exception&){ // Ignore incomplete files
    }catch(const bad_format_exception& ex) {
        std::cerr << "InterOp did not have the expected format: " << ex.what() << std::endl;
        return 1;
    }catch(const file_not_found_exception& ex){
        std::cerr << "Count not find InterOp file: " << ex.what() << std::endl;
        return 1;
    }
    return 0;
}

void set_status(rarp::NLISTS & rgs,char const * status,int userid) {

    for(unsigned int p=0; p<rgs.size(); ++p){
        char update[1024];
        sprintf(update,"update prepT set status = '%s', status_time=unix_timestamp() where prepid = %s ; select row_count() updated",status,rgs[p]["prepid"].data());
        // std::cout << "using " << update << "\n";
        rarp::NLIST x = db::get_named_row("seqdb",update); 
        
        // assert(x["updated"]=="1");
        if(x["updated"]!="1") {
            std::cout << "unable to update status\n"; // never seen this before?!?
            fflush(stdout);
            sleep(1);
        }

        // sprintf(update,"insert into statusT (status,status_time,prepid,userid,poolid,seqid,platename) values ('%s',unix_timestamp(),%s,%d,0,0,'') ; select row_count() inserted",
        sprintf(update,"replace into statusT (status,status_time,prepid,userid,poolid,seqid) values ('%s',unix_timestamp(),%s,%d,0,0) ; select row_count() inserted",
          status,rgs[p]["prepid"].data(),userid);

        x.clear(); 

        x = db::get_named_row("seqdb",update); 
        assert(x["inserted"]=="1");

    }
}

void update_fc_table(std::string const & run_dir, std::string const & fcid, 
  std::string const & archive_dir, tinyxml2::XMLElement /* so lazy */ *& rp) {

    using namespace std;

    // why on earch set LOCATION (seqsataloc) now?!?
    string ts_rta = get_single_line_output_as_string( "ls -l --full-time %s | awk '{print $6\" \"$7}' | cut -f1 -d. ",(run_dir+"/RTAComplete.txt").data()), 
      ts_r1 = get_single_line_output_as_string( "ls -l --full-time %s | awk '{print $6\" \"$7}' | cut -f1 -d. ",(run_dir+"/RunInfo.xml").data()); 

    rarp::NLISTS hmm;
    db::get_named_rows("seqdb",hmm,"select fc_insert,fc_update,dateread1,daterta,rtaver,hcsver from Flowcell where fcillumid = '%s'",fcid.data());
    assert(hmm.size()==1);
    rarp::NLIST & hff = hmm[0];

    if(hff["daterta"]!="0000-00-00 00:00:00") {
            cout << "Flow cell info already set\n\n", assert(hff["rtaver"]==rp->FirstChildElement( "RtaVersion" )->GetText()), 
              assert(hff["rtaver"]==rp->FirstChildElement( "RtaVersion" )->GetText()), assert(hff["daterta"]==ts_rta), assert(hff["dateread1"]==ts_r1);
    } else{
        char update[8*1024];
        sprintf(update,"update Flowcell set " 
          "rtaver='%s', hcsver='%s', "
          "dateread1='%s', daterta='%s', "
          "datebcl=unix_timestamp(), seqsataloc ='%s' "
          "where fcillumid = '%s'; select row_count() as affected",
          rp->FirstChildElement( "RtaVersion" )->GetText(), rp->FirstChildElement( "ApplicationVersion" )->GetText(),
          ts_r1.data(), ts_rta.data(), 
          archive_dir.data(), 
          fcid.data() 
        );
        // cout << "Update flowcell info " << update << "\n";
        rarp::NLIST x = db::get_named_row("seqdb",update);
        assert(x["affected"]=="1");
    }

    // for(rarp::NLIST::iterator it=hff.begin(); it!=hff.end(); it++) { cout << "\t[" << it->first << "] "<< it->second << "\n"; }
    // cout << "\n";

}

void interop_bits(std::string const & run_dir, std::string const & fcid, rarp::NLISTS & ns) {

    using namespace std;
    // cout << "pull metrics from " << run_dir << "\n\n";

    using namespace illumina::interop::model::metric_base;
    using namespace illumina::interop::model::metrics;
    using namespace illumina::interop::io;
    using namespace illumina::interop::util;
    using namespace illumina::interop;
    using namespace illumina::interop::logic::summary;
    using namespace illumina::interop::model::summary;

    illumina::interop::model::metrics::run_metrics rm;
    std::vector<unsigned char> valid_to_load;
    logic::utils::list_summary_metrics_to_load(valid_to_load); // only the ??!?
    rm.read(run_dir,valid_to_load);
    run_summary summary;
    // summarize_run_metrics(rm, summary);//, skip_median_calculation);

    // cout << "summary lane count (pre-init) = " << summary.lane_count() << "\n";
    summarize_run_metrics(rm,summary,true,true);

    assert(atoi(ns[0]["lane_num"].data())==(int)summary.lane_count()); // no risk here

    // int lane=0;

    // get rid of this and just use below...

    for (unsigned int r=0; r<summary.size();++r){

        for(unsigned int l=0;l<summary.lane_count();++l) {

            // cout << "read-summary (l,index,r)=" << r << ", lane="<< l << " tile-count="<<summary[r][l].tile_count() 
            if(r==0) {

                // cout << "THIS SHOULD BE MERGED\n\n";
                rarp::NLISTS hmm;
                db::get_named_rows("seqdb",hmm,"select distinct(clustden),clustdenstdev,clusterpf,clusterpfstdev from Lane l join Flowcell f on l.fcid=f.fcid where fcillumid = '%s' and lanenum = '%d'",fcid.data(),l+1);
                assert(hmm.size()==1);
                rarp::NLIST & hff = hmm[0];

                // this will prolly error out with unset data?!?
                float clustden =    (float) summary.at(r).at(l).density().mean()/1000,
                  clustdenstdev =   (float) summary.at(r).at(l).density().stddev()/1000,
                  clusterpf =       (float) summary.at(r).at(l).density_pf().mean()/1000,
                  clusterpfstdev =  (float) summary.at(r).at(l).density_pf().stddev()/1000;

                if(hff["clustden"]!="NULL") {
                    assert(int(clustden)==atoi(hff["clustden"].data()));
                    assert(int(clustdenstdev)==atoi(hff["clustdenstdev"].data()));
                    // cout << "getting mismatch for " << clusterpf << " and " << hff["clusterpf"].data() << "\n";
                    assert(int(clusterpf)==atoi(hff["clusterpf"].data()) || int(clusterpf+0.1)==atoi(hff["clusterpf"].data()));
                    // assert(int(clusterpf)==atoi(hff["clusterpf"].data()));
                    assert(int(clusterpfstdev)==atoi(hff["clusterpfstdev"].data()));
                } else {
                    char update[1024];
                    sprintf(
                      update,"update Lane l join Flowcell f on l.fcid=f.fcid set rg_status='sequenced',clustden=%0.2f,clustdenstdev=%0.2f,clusterpf=%0.2f,clusterpfstdev=%0.2f "
                      "where lanenum=%d and fcillumid='%s' ; select row_count() updated", clustden, clustdenstdev, clusterpf, clusterpfstdev, l+1, fcid.data() 
                    );
                    rarp::NLIST x = db::get_named_row("seqdb",update);
                    assert(x["updated"]!="0");
                }

                // for(rarp::NLIST::iterator it=hff.begin(); it!=hff.end(); it++) { cout << "\t[" << it->first << "] "<< it->second << "\n"; }
            }
        }
    }

    assert( (summary.size()==3 && ns[0]["leni2"]=="0") || summary.size()==4);

    tinyxml2::XMLDocument doc2;
    assert(isregfile((run_dir+"/RunInfo.xml").data()));
    doc2.LoadFile((run_dir+"/RunInfo.xml").data());

    int r1of=0, i1of=1, r2of= summary.size()==3?2:3;
    // assert(summary.size()==3||summary.size()==4);
    
    for(unsigned int l=0;l<summary.lane_count();++l) {

        ///////// THESE METRICS BELONG IN FLOWCELL NOT LANE!?!?!

        {
        rarp::NLISTS hmm;
        db::get_named_rows("seqdb",hmm,"select distinct(clustden),clustdenstdev,clusterpf,clusterpfstdev,perq30r1,perq30r2,perq30i1,errorr1,errorr2,percentalignr1,percentalignr2 from Lane l join Flowcell f on l.fcid=f.fcid where fcillumid = '%s' and lanenum = '%d'",fcid.data(),l+1);
        assert(hmm.size()==1);
        rarp::NLIST & hff = hmm[0];

        float perq30r1=(float) summary.at(r1of).at(l).percent_gt_q30(), 
          perq30r2=(float) summary.at(r2of).at(l).percent_gt_q30(), 
          perq30i1=(float) summary.at(i1of).at(l).percent_gt_q30(),
          errorr1          =(isnan(summary.at(r1of).at(l).error_rate().mean())?0.0:summary.at(r1of).at(l).error_rate().mean()),
          errorr2          =(isnan(summary.at(r2of).at(l).error_rate().mean())?0.0:summary.at(r2of).at(l).error_rate().mean()),
          percentalignr1   =summary.at(r1of).at(l).percent_aligned().mean(), 
          percentalignr2   =summary.at(r2of).at(l).percent_aligned().mean();

        if(hff["perq30r1"]!="NULL") { 
            assert(int(perq30r1)==atoi(hff["perq30r1"].data()));
            assert(int(perq30r2)==atoi(hff["perq30r2"].data()));
            assert(int(perq30i1)==atoi(hff["perq30i1"].data()));
            assert(int(errorr1)==atoi(hff["errorr1"].data()));
            assert(int(errorr2)==atoi(hff["errorr2"].data()));
            assert(int(percentalignr1)==atoi(hff["percentalignr1"].data()));
            assert(int(percentalignr2)==atoi(hff["percentalignr2"].data()));
        }else{

            char update[8*1024];
            sprintf(update,"update Lane l join Flowcell f on l.fcid=f.fcid set "
              "clustden=%0.2f,clustdenstdev=%0.2f," "clusterpf=%0.2f,clusterpfstdev=%0.2f, " ///// why the f' does the original set these separately from the rest?!?
              "perq30r1=%0.3f,perq30r2=%0.3f,perq30i1=%0.3f,"
              "errorr1=%0.3f,errorr2=%0.3f,"
              "percentalignr1=%0.3f,percentalignr2=%0.3f "
              "where lanenum=%d and fcillumid='%s' ; select row_count() updated ",
              (float) summary.at(r1of).at(l).density().mean()/1000, (float) summary.at(r1of).at(l).density().stddev()/1000,
              (float) summary.at(r1of).at(l).density_pf().mean()/1000, (float) summary.at(r1of).at(l).density_pf().stddev()/1000,
              (float) summary.at(r1of).at(l).percent_gt_q30(), (float) summary.at(r2of).at(l).percent_gt_q30(), (float) summary.at(i1of).at(l).percent_gt_q30(),
              (isnan(summary.at(r1of).at(l).error_rate().mean())?0.0:summary.at(r1of).at(l).error_rate().mean()),
              (isnan(summary.at(r2of).at(l).error_rate().mean())?0.0:summary.at(r2of).at(l).error_rate().mean()),
              summary.at(r1of).at(l).percent_aligned().mean(), summary.at(r2of).at(l).percent_aligned().mean(),
              l+1, fcid.data()
            );

            rarp::NLIST x = db::get_named_row("seqdb",update);
            /////// this should be lane not rg!?!?
            assert(x["updated"]!="0");
        }
        }
    }
}

void write_sample_sheet(
    std::string const & fcid, 
    std::string const & ssfn, 
    rarp::NLISTS & ns, 
    tinyxml2::XMLElement /* so lazy */ *& rp,
    rarp::NLISTS & rgs,
    bool legacy ) {

    using namespace std;

    rarp::NLIST & info = ns[0];
    stringstream ss;

    if(legacy) {
        ss << "[Header]\nIEMFileVersion,4\nInvestigator Name,"<<info["fullname"]<<"\nExperiment Name,"<<fcid<<"\nDate,"
          <<rp->FirstChildElement( "RunStartDate" )->GetText()<<"\nWorkflow,GenerateFASTQ\nApplication,PlatformType FASTQ Only\nAssay,LibType\n"
            "Description,"<<info["projects"]<<"\nChemistry,Default\n\n[Reads]\n"<<info["lenr1"]<<"\n"<<info["lenr2"]<<"\n\n"
            "[Data]\nLane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description\n";
    }else{
        ss << "[Header],,,,,\nWorkflow,GenerateFASTQ,,,,\n,,,,,\n[Settings],,,,,\n,,,,,\n[Data],,,,,\nLane,Sample_ID,Sample_Name,Sample_Project,index,index2\n"; 
    }

    assert(rgs.size()>MIN_RGS);

    for(unsigned int p=0; p<rgs.size(); ++p){

        assert(insane_frac.count(ns[0]["chemver"]));
        assert(insane_frac[ns[0]["chemver"]].count(rgs[p]["sample_type"]));
        // cout << "temp check " << rgs[p]["adapter"] << "\n"; // ", " << rgs[p]["old_adapter"] << "\n";

        if(legacy) {
            assert(strchr(rgs[p]["adapter"].data(),',')==0);
            char silly_l_frac_s[256];
            sprintf(silly_l_frac_s,"%0.4f",(float) 1.0/insane_frac[ns[0]["chemver"]][rgs[p]["sample_type"]]);
            ss << rgs[p]["lanenum"] << ","<<rgs[p]["sample_internal_name"] << ","<<rgs[p]["sample_internal_name"] <<",,,,"<< rgs[p]["adapter"] 
            <<","<< rgs[p]["project"] 
            <<","<<silly_l_frac_s<<"_"
            << rgs[p]["silly_desc"] << "pM\n";
        } else { 
            assert(strchr(rgs[p]["adapter"].data(),',')!=0);
            ss << rgs[p]["lanenum"] << ","<<rgs[p]["sample_internal_name"] << ","<<rgs[p]["sample_internal_name"] <<"," << rgs[p]["project"] << ","<< rgs[p]["adapter"] << "\n";
        }
            
    }

    ofstream ssf(ssfn.data());
    if(!ssf) cout << "couldn't open sample sheet\n",exit(1);
    else cout << "writting sample sheet\n";
    ssf << ss.str();
    ssf.close();

}

void print_info(tinyxml2::XMLElement /* so lazy */ *& rp) {

    printf( "\t%-30s\t%s\n","WorkflowType",rp->FirstChildElement( "WorkflowType" )->GetText());
    printf( "\t%-30s\t%s\n","RtaVersion",rp->FirstChildElement( "RtaVersion" )->GetText());
    printf( "\t%-30s\t%s\n","Side",rp->FirstChildElement( "Side" )->GetText());
    printf( "\t%-30s\t%s\n","ExperimentName",rp->FirstChildElement( "ExperimentName" )->GetText());
    printf( "\t%-30s\t%s\n","Read1NumberOfCycles",rp->FirstChildElement( "Read1NumberOfCycles" )->GetText());
    printf( "\t%-30s\t%s\n","Read2NumberOfCycles",rp->FirstChildElement( "Read2NumberOfCycles" )->GetText());
    printf( "\t%-30s\t%s\n","IndexRead1NumberOfCycles",rp->FirstChildElement( "IndexRead1NumberOfCycles" )->GetText());
    printf( "\t%-30s\t%s\n","IndexRead2NumberOfCycles",rp->FirstChildElement( "IndexRead2NumberOfCycles" )->GetText());
    printf( "\t%-30s\t%s\n","RunNumber",rp->FirstChildElement( "RunNumber" )->GetText());
    printf( "\t%-30s\t%s\n","InstrumentName/MachineName",rp->FirstChildElement( "InstrumentName" )->GetText());
    printf( "\t%-30s\t%s\n","Application",rp->FirstChildElement( "Application" )->GetText());
    printf( "\t%-30s\t%s\n","ApplicationVersion",rp->FirstChildElement( "ApplicationVersion" )->GetText());
    printf( "\t%-30s\t%s\n","RunStartDate/RunDate",rp->FirstChildElement( "RunStartDate" )->GetText());
    printf( "\t%-30s\t%s\n","RunId/RunFolder",rp->FirstChildElement( "RunId" )->GetText());
    printf( "\t%-30s\t%s\n","FlowCellSerialBarcode/FCID",rp->FirstChildElement("RfidsInfo")->FirstChildElement( "FlowCellSerialBarcode" )->GetText());

}

void update_actual_lane_fractions_from_html(rarp::NLISTS & hs,rarp::NLISTS & rgs,std::string const & fcid, int & ly, char const * fcyield) {

    using namespace std;

    stringstream tmp;
    tmp << "<style><th, td {padding: 5px;}\nth, td {text-align: center;}\nth ding: 5px;}\n</style>\n"
      << "<table border=\"1\" style=\"border:1px solid black;border-collapse:collapse;width:95%\">\n"
      << "<tr bgcolor=\"#d7d7d7\"><th><font color=\"#2996cc\">SampleID</font></th>"
      << "<th><font color=\"#2996cc\">Lane</font></th>"
      << "<th><font color=\"#2996cc\">PoolName</font></th>"
      << "<th><font color=\"#2996cc\">Sampletype</font></th>"
      << "<th><font color=\"#2996cc\">CaptureKit</font></th>"
      << "<th><font color=\"#2996cc\">Yield (Mb)</font></th></tr>\n";

    // for(int i=0; i<rgs.size(); ++i) { for(rarp::NLIST::iterator it=rgs[i].begin(); it!=rgs[i].end(); it++) cout << "["<<i<<"]["<<it->first<<"] "<<it->second<<"\n"; }
    map<string, map<string, int> > ordering;
    for(unsigned int i=0; i<rgs.size(); ++i) ordering[rgs[i]["lanenum"]][rgs[i]["sample_internal_name"]]=i;

    string cur_lane; // rgs[0]["lanenum"];
    for(unsigned int i=0; i<hs.size(); ++i) { // for(int i=0, i2=0; i<hs.size(); ++i,++i2) {

        int i2 = ordering[hs[i]["Lane"]][hs[i]["Sample"]];

        char bored[1024], *bp=bored;
        memset(bored,0,1024);
        for(unsigned int y=0;y<hs[i]["Yield (Mbases)"].length();++y) if(hs[i]["Yield (Mbases)"][y]!=',') *bp++=hs[i]["Yield (Mbases)"][y];
        // WHY ISN'T YIELD/SAMPLE DONE HERE?!? 
        ly+=atoi(bored);

        if(hs[i]["Sample"]=="Undetermined") {
            tmp << "<tr bgcolor=\"#d7d7d7\"><td>Undetermined</td><td>" << cur_lane << "</td><td>NA</td><td>NA</td><td>NA</td><td>"<<hs[i]["Yield (Mbases)"]<<"</td></tr>\n";
            continue;
        }
        cur_lane=rgs[i2]["lanenum"];

        assert(hs[i]["Sample"]==rgs[i2]["sample_internal_name"]);

        tmp << "<tr><td>"<<rgs[i2]["sample_internal_name"]<<"</td><td>"<<rgs[i2]["lanenum"]<<"</td><td>"<<rgs[i2]["pool_name"]<<"</td><td>"<<rgs[i2]["sample_type"]
          <<"</td><td>"<<rgs[i2]["capture_kit"]<<"</td><td>"<<hs[i]["Yield (Mbases)"]<<"</td></tr>\n";

        if(rgs[i2]["lnfractionact"]=="NULL"){
            cout << "-> Lane(3) info not been set\n";
            char update[8*1024];
            sprintf(update,"update Lane l join prepT p on l.prepid=p.prepid join Flowcell f on l.fcid=f.fcid set lnfractionact = %s, lnyield = %s "
              "where sample_internal_name = '%s' and lanenum=%s and fcillumid='%s'; select row_count() affected ", hs[i]["% of the"].data(), bored, 
              rgs[i2]["sample_internal_name"].data(), rgs[i2]["lanenum"].data(), fcid.data()
            ); 
            // cout << "\nusing " << update << "\n\n";
            rarp::NLIST x = db::get_named_row("seqdb",update);
            assert(x["affected"]=="1");
        }else{
            int rarp = 100.0*atof(rgs[i2]["lnfractionact"].data());
            assert(rarp==int(100.0*atof(hs[i]["% of the"].data())));
            assert(rgs[i2]["lnyield"]==bored);
        }
    }

    char subject[1024];
    sprintf(subject,"Flowcell %s lane/sample yield summaries (%d/%s)",fcid.data(),ly,fcyield);
    // { char what[2048]; getcwd(what,sizeof what); cout << "using " << what << "\n"; }
    // sprintf(fer2, "/nfs/goldstein/software/mutt-1.5.23/bin/mutt -e \"set content_type=text/html\" -s \"Flowcell %s lane/sample yield summaries (%d/%s)\" dh2880@cumc.columbia.edu < " FLOWCELL_REPORT,fcid.data(),ly,fcyield);
    // sendmail("nb2975@cumc.columbia.edu","nb2975@cumc.columbia.edu",subject,tmp.str().data(),true);

}
    
void move_fastq_to_sample_specific_dirs(
  rarp::NLISTS & /* hs */, 
  rarp::NLISTS & rgs, 
  std::string const & bcl_dir, 
  std::string const & bcl, 
  char const * netapp_out_dir_tmp, 
  std::string const & run_dir_base, 
  std::string const & fcid, 
  std::string const & archive_dir, 
  bool archive = 0
) {

    using namespace std;

    string blarp;
    int S=0;
    char s_[256];
    blarp=rgs[0]["lanenum"];
    std::set<string> proj_dirs_to_check;
    vector<struct ARCHIVE> ARSV2;
    
    map<string,int> S_NUM;
    int S_N=0, C_S_N;

    for(unsigned int p=0; p<rgs.size(); ++p){

        rarp::NLIST & SE = rgs[p];

        if(blarp!=SE["lanenum"]) {
            blarp=SE["lanenum"];
            S=0;
        }

        if(S_NUM.count(rgs[p]["sample_internal_name"])) C_S_N=S_NUM[rgs[p]["sample_internal_name"]];
        else C_S_N=S_NUM[rgs[p]["sample_internal_name"]]=++S_N;

        ++S;
        sprintf(s_,"%d",C_S_N);
        
        // cout << "\n[se]: sample_internal_name=" << rgs[p]["sample_internal_name"]<<", lanenum="<<rgs[p]["lanenum"] << ", p=" << p << ", S=" << S << ", C_S_N=" << C_S_N << ", S_N=" << S_N << "\n";

        string ssb= string(netapp_out_dir_tmp) + "/FASTQ/" + fcid +"/",
          SAMPLE_SPECIFIC_TMP_DIR = ssb + rgs[p]["seq_type_upper"] + "/" + rgs[p]["sample_internal_name"] + "/" + fcid,
          read1 = rgs[p]["sample_internal_name"] + "_S" + s_ + "_L00" + rgs[p]["lanenum"] + "_R1_001.fastq.gz",
          read2 = rgs[p]["sample_internal_name"] + "_S" + s_ + "_L00" + rgs[p]["lanenum"] + "_R2_001.fastq.gz",
          se_archive_dir = seq::se_archive_dir(SE);

        mkdirs(se_archive_dir.data());

        /*
        //// not necessary...
        if( isregfile( (se_archive_dir + "/" + read1 ).data() )
            && !isregfile( (SAMPLE_SPECIFIC_TMP_DIR + "/" + read1).data() )
            && !isregfile( (bcl_dir + rgs[p]["project"] + "/" + read1 ).data() ) ) {
            cout << "This is already archived AND cleaned up!?!?!\n";
            ///// check storage complete!?!
            exit(1);
            continue;
        }
        if(!isregfile( (bcl_dir + rgs[p]["project"]+ "/" + read1 ).data()) && !isregfile( (SAMPLE_SPECIFIC_TMP_DIR + "/" + read1 ).data() ) ) {
            { 
            rarp::LISTS check;
            ///// should use man 3 glob?!?
            get_table(check,("find "+bcl_dir+rgs[p]["project"]+"/ "+SAMPLE_SPECIFIC_TMP_DIR+"/ -maxdepth 1 -name '"+rgs[p]["sample_internal_name"] + "_S*_L00" + rgs[p]["lanenum"] + "_R1_001.fastq.gz'").data());
            cout << check.size() << "\n";
            assert(check.size()==1);
            char *b1 = strdup(check[0][0].data()), *b2=basename(b1), *b3=b2+rgs[p]["sample_internal_name"].length(), *b4=strchr(b3+2,'_');
            free(b1);
            assert(memcmp(b3,"_S",2)==0), assert(b4);
            char what[1024];
            memset(wt,0,1024);
            memcpy(wt,b3+2,b4-b3-2);
            cout << "we are swapping sample number " << wt << " for " << s_ << "\n";
            strcpy(s_,wt);
            read1 = rgs[p]["sample_internal_name"] + "_S" + s_ + "_L00" + rgs[p]["lanenum"] + "_R1_001.fastq.gz",
            read2 = rgs[p]["sample_internal_name"] + "_S" + s_ + "_L00" + rgs[p]["lanenum"] + "_R2_001.fastq.gz";
            }
        }
        */

        if(!archive) {

            string hdir = bcl + run_dir_base + "/Reports/html/" + fcid + "/" + rgs[p]["project"] +"/" + rgs[p]["sample_internal_name"] + "/all/",
            pdir = bcl_dir + rgs[p]["project"];

            proj_dirs_to_check.insert(pdir);

            // cout << "\nfor transfer(1)\t" << bcl_dir << rgs[p]["project"] << "\t" << SAMPLE_SPECIFIC_TMP_DIR << "\n";

            bool some_how_done = false;

            if(isregfile( (pdir + "/" + read1 ).data()) ) {

                // cout << "'may' not yet be moved\n";
                assert( isregfile( (pdir + "/" + read2 ).data() ) );

                // assert(SE["data"]=="NULL");

                link_file(hdir.data(),SAMPLE_SPECIFIC_TMP_DIR.data(),"laneBarcode.html");
                link_file(pdir.data(),SAMPLE_SPECIFIC_TMP_DIR.data(),read1.data());
                link_file(pdir.data(),SAMPLE_SPECIFIC_TMP_DIR.data(),read2.data());

                { string json = get_json( SAMPLE_SPECIFIC_TMP_DIR.data(), read1.data(), read2.data() );

                // should be in conditional but for now
                char update[8*1024]; 
                sprintf(update,"update Lane l join prepT p on l.prepid=p.prepid join Flowcell f on l.fcid=f.fcid set "
                  "status = 'Mapping', status_time = unix_timestamp(), rg_status = 'fastq_ready', step_status = 'in storage', "
                  "datestor=current_timestamp(), seqsataloc='%s', pipelinecomplete = 1, "
                  "data = '%s' "
                  "where sample_internal_name = '%s' and lanenum=%s and fcillumid='%s' ; select row_count() as affected ;", 
                  archive_dir.data(),json.data(), SE["sample_internal_name"].data(), SE["lanenum"].data(), fcid.data() 
                ); 

                // clearly this shouldn't be here!?!?
                // if(SE["data"]=="NULL") {
                // cout << "\nusing:\n" << update << "\n\n";
                rarp::NLISTS urows; // rarp::NLIST urow;
                db::get_named_rows("seqdb",urows,update); // get_named_row(update);
                // cout << "we got " << urows[0]["affected"] << "\n";
                // setting prepT.status as well as Lane...
                assert(urows[0]["affected"]!="0"); } 

                // assert(urows[0]["affected"]=="2");
                // }

                { struct stat s1;
                stat( (SAMPLE_SPECIFIC_TMP_DIR+"/"+read1).data() ,&s1);
                struct stat s2;
                stat( (SAMPLE_SPECIFIC_TMP_DIR+"/"+read2).data() ,&s2);
                assert(s1.st_nlink==2);
                assert(s2.st_nlink==2);
                char update[8*1024]; 
                // for(rarp::NLIST::iterator it=SE.begin(); it!=SE.end(); it++) cout << "["<<it->first<<"] "<<it->second<<"\n";
                sprintf(update,"select data, data->'$.fastq.r1.size' rs1, data->'$.fastq.r2.size' rs2, data->'$.fastq.r1.modification' rm1, data->'$.fastq.r2.modification' rm2 from "
                  "Lane where prepid = '%s' and lanenum=%s and fcid='%s'",SE["prepid"].data(), SE["lanenum"].data(), SE["fcid"].data() );
                // cout << "using " << update << "\n";
                rarp::NLIST x = db::get_named_row("seqdb",update); 
                // cout << "checking json - paranoia\n";
                // for(rarp::NLIST::iterator it=x.begin(); it!=x.end(); it++) cout << "["<<it->first<<"] "<<it->second<<"\n";
                // cout << s1.st_size << ", " << atoi(x["rs1"].data()) << "\n";

                assert(s1.st_size==atol(x["rs1"].data())), assert(s1.st_mtime==atol(x["rm1"].data())), // no need with mtime!?!
                  assert(s2.st_size==atol(x["rs2"].data())), assert(s2.st_mtime==atol(x["rm2"].data()));

                // cout << "checked json\n";
                // cout << "we remove " << pdir << "/" << read1 << "\n";
                assert(unlink( (pdir + "/" + read1).data() )==0);
                // cout << "we remove " << pdir << "/" << read2 << "\n";
                assert(unlink( (pdir + "/" + read2).data() )==0); }
                
                /* for (int i=0;i<urows.size();++i){ rarp::NLIST & urow = urows[i];
                    for(rarp::NLIST::iterator it=urow.begin(); it!=urow.end(); it++) {
                        cout << "[" << i << "] : ["<<it->first<<"] '"<<it->second<<"'\n"; } } */

                /* char cmd[2048]; ////// or use rename(full1,full2)
                sprintf(cmd,"mv %s/%s %s",pdir.data(),read2.data(),SAMPLE_SPECIFIC_TMP_DIR.data()); */

                // cout << "now we check files and update db\n";

            }else if(isregfile( (SAMPLE_SPECIFIC_TMP_DIR + "/" + read1 ).data() ) ) {
                assert(isregfile( (SAMPLE_SPECIFIC_TMP_DIR + "/" + read2 ).data() ) );
                assert(!isregfile( (pdir+ "/" + read1 ).data() ) );
                assert(!isregfile( (pdir+ "/" + read2 ).data() ) );
            }else {
                // cout << "INITIAL: " << pdir << "\n";
                // cout << "SSL: " << SAMPLE_SPECIFIC_TMP_DIR << "\n";
                // cout << "ARCHIVE: " << se_archive_dir + "/" + read1 << "\n";
               
                /* if( isregfile( (bcl_dir + "/StorageComplete").data() )) {
                    cout << "this was processed via legacy procedure adn already cleaned up?!?\n";
                    char arsv2[1024];
                    sprintf(arsv2,"update Flowcell set fc_status = 'legacy' where FCillumID = '%s' ; select row_count() locked",fcid.data());
                    rarp::NLIST x = get_named_row(arsv2); assert(x["locked"]=="1"); cout << "we have a lock on " << fcid << "\n";
                }else */

                // cout << "this doesn't exist " <<  pdir << "/" << read1 <<"\n";
                if( 
                    isregfile( (se_archive_dir + "/" + read1 ).data() )
                    && isregfile( ( se_archive_dir + "/" + read1 + ".md5sum").data() )
                ) {
                    // cout << "already archived and cleaned up!?!\n";
                    some_how_done=true;
                    // continue;
                }else{
                    cout << "This is plain wrong!?!\n\n";
                    cout << "check yield " << rgs[p]["lnyield"] <<"\n";
                    if(rgs[p]["lnyield"]=="0") {
                        // perhaps leave till next round?!?
                        // this should trigger an email - i.e. the pool/sample was left out?!?
                        // cout << "we should really make this use int and fail out the lane entry if <500?!?\n";
                        continue;
                    }else{
                        assert(0);
                        exit(1);
                    }
                }
            }
            
            // rsync uses timestamp so will completely re-run if bcl2fastq re-ran!?!?
            if(SE["data"]=="NULL") { 
                //// don't over-write as mapping data alread gone in!?!
                // cout << "MOFO '" << SE["data"] << "'\n";
                // assert(0);
                string json = get_json( SAMPLE_SPECIFIC_TMP_DIR.data(), read1.data(), read2.data() );
                char update[8*1024]; 
                sprintf(update,"update Lane l join prepT p on l.prepid=p.prepid join Flowcell f on l.fcid=f.fcid set "
                    "status = 'Mapping', status_time = unix_timestamp(), rg_status = 'fastq_ready', step_status = 'in storage', "
                    "datestor=current_timestamp(), seqsataloc='%s', pipelinecomplete = 1, "
                    "data = '%s' "
                    "where sample_internal_name = '%s' and lanenum=%s and fcillumid='%s' ; select row_count() as affected ;", 
                    archive_dir.data(),json.data(), SE["sample_internal_name"].data(), SE["lanenum"].data(), fcid.data() 
                ); 
                // cout << "\nusing:\n" << update << "\n\n";
                rarp::NLISTS urows; // rarp::NLIST urow;
                db::get_named_rows("seqdb",urows,update); // get_named_row(update);

            }

            if(!some_how_done) {

                struct stat s1; stat( (SAMPLE_SPECIFIC_TMP_DIR+"/"+read1).data() ,&s1); assert(s1.st_nlink==1);
                struct stat s2; stat( (SAMPLE_SPECIFIC_TMP_DIR+"/"+read2).data() ,&s2); assert(s2.st_nlink==1);
                char update[8*1024]; 
                // for(rarp::NLIST::iterator it=SE.begin(); it!=SE.end(); it++) cout << "["<<it->first<<"] "<<it->second<<"\n";
                sprintf(update,"select data, data->'$.fastq.r1.size' rs1, data->'$.fastq.r2.size' rs2, data->'$.fastq.r1.modification' rm1, data->'$.fastq.r2.modification' rm2 from "
                    "Lane where prepid = '%s' and lanenum=%s and fcid='%s'",SE["prepid"].data(), SE["lanenum"].data(), SE["fcid"].data() );
                // cout << "using " << update << "\n";
                rarp::NLIST x = db::get_named_row("seqdb",update); 
                // cout << "checking json - paranoia\n";
                // for(rarp::NLIST::iterator it=x.begin(); it!=x.end(); it++) cout << "["<<it->first<<"] "<<it->second<<"\n";
               
                assert(s1.st_size==atol(x["rs1"].data())), assert(s1.st_mtime==atol(x["rm1"].data())), // no need with mtime!?!
                    assert(s2.st_size==atol(x["rs2"].data())), assert(s2.st_mtime==atol(x["rm2"].data()));

                // cout << "checked json\n";
            } 

        }else{

            string r1path = SAMPLE_SPECIFIC_TMP_DIR + "/" + read1, r2path = SAMPLE_SPECIFIC_TMP_DIR + "/" + read2;
            mkdirs(se_archive_dir.data());

            // cout << "checking " << r1path<< ", " << r2path<< "\n";
            // bool some_how_done = false;

            if(!isregfile( r1path.data())||!isregfile(r2path.data())){

                ////// need to start using l.id?!?
                if(rgs[p]["lnyield"]=="0") {
                    // this should trigger an email - i.e. the pool/sample was left out?!?
                    continue;
                }

                assert(isregfile( (se_archive_dir + "/" + read1 ).data() ));
                assert(isregfile( (se_archive_dir + "/" + read1 +".md5sum").data() ));
                assert(isregfile( (se_archive_dir + "/" + read2 ).data() ));
                assert(isregfile( (se_archive_dir + "/" + read2 +".md5sum").data() ));

            }else{

                // cout << r1path << ", " << r2path << "\n";
                // cout << se_archive_dir << "\n\n";
                assert(isregfile( r1path.data() ));
                assert(isregfile( r2path.data() ));

                stringstream ss;

                // char update[8*1024]; sprintf(update,"update Lane l join prepT p on l.prepid=p.prepid join Flowcell f on l.fcid=f.fcid set "
                // json.data(), SE["sample_internal_name"].data(), SE["lanenum"].data(), fcid.data() ); se_archive_dir
                // json.data(), SE["sample_internal_name"].data(), SE["lanenum"].data(), fcid.data() ); se_archive_dir
        
                struct ARCHIVE a;
                a.data=SE["data"], a.sample_internal_name=SE["sample_internal_name"], a.lanenum=SE["lanenum"], a.archive_dir=se_archive_dir;

                // do this in popen?!?
                if(!isregfile( (se_archive_dir+"/laneBarcode.html").data() )) {
                    ss << "rsync -azpvv -L --no-owner " << SAMPLE_SPECIFIC_TMP_DIR << "/laneBarcode.html " << se_archive_dir << " > /dev/null";
                    a.html=ss.str(); ss.str(""), ss.clear();    
                }

                ss << "rsync -azpvv -L --no-owner " << r1path << " " << se_archive_dir << " > /dev/null";
                a.read1=ss.str(); ss.str(""), ss.clear();   

                ss << "rsync -azpvv -L --no-owner " << r2path << " " << se_archive_dir << " > /dev/null";
                a.read2=ss.str(); ss.str(""), ss.clear();   

                ARSV2.push_back(a);

            }
            // for(rarp::NLIST::iterator it=SE.begin(); it!=SE.end(); it++) { cout << "\t[" << it->first << "] "<< it->second << "\t"; } cout << "\n";
        }

    }

    if(!archive) for(std::set<string>::iterator it = proj_dirs_to_check.begin(); it != proj_dirs_to_check.end(); it++) {

        string gc = get_single_line_output_as_string( "ls %s/*.gz 2>/dev/null | wc -l",it->data());
        // cout << "we have archive_dir " << archive_dir << "\n";
        // cout <<"MUST CHECK PROJ DIR " << *it << " HAS COUNT OF " << gc << "\n";
        assert(gc=="0");
        // rarp::NLIST x = get_named_row("INSERT INTO statusT (sample_internal_name,STATUS_TIME,STATUS,sample_id,PREPID,USERID,POOLID,SEQID,PLATENAME) SELECT DISTINCT(pt.sample_internal_name),UNIX_TIMESTAMP(),'Archiving',pt.sample_id,pt.PREPID,'%d',0,0,' ' "

    }else {

        #ifdef _OPENMP
        omp_set_num_threads(8);
        #endif

        // #pragma omp parallel 

        for(unsigned int i = 0; i < ARSV2.size(); ++i ) { // for(std::set<string>::iterator it = proj_dirs_to_check.begin(); it != proj_dirs_to_check.end(); it++) {

            // cout << "run this " << *it << "\n";

            if(!ARSV2[i].html.empty()) {
                #pragma omp critical
                { cout << "we have " << ARSV2[i].data<<"\n\n";
                cout << "["<<i<<"] archive html " << ARSV2[i].html << "\n"; }
                if(system(ARSV2[i].html.data())) cout << "this is not good\n",exit(1);
            }

            // if(system(it->data())) cout << "this is not good\n",exit(1);
            
            #pragma omp critical
            { cout << "["<<i<<"] archive read1 " << ARSV2[i].read1 << "\n"; }
            if(system(ARSV2[i].read1.data())) cout << "this is not good\n",exit(1);

            #pragma omp critical
            { cout << "["<<i<<"] archive read2 " << ARSV2[i].read2 << "\n"; } 
            if(system(ARSV2[i].read2.data())) cout << "this is not good\n",exit(1);

            // move into main loop above since this is really just a warm up copy?!?
            char fer2[16*1024];
            sprintf(fer2,"update Lane l join prepT p on l.prepid=p.prepid join Flowcell f on l.fcid=f.fcid set data = json_set(data,'$.fastq.path.archive','%s') "
              "where sample_internal_name = '%s' and lanenum = %s and fcillumid='%s' ; select row_count() as affected ",ARSV2[i].archive_dir.data(),ARSV2[i].sample_internal_name.data(),ARSV2[i].lanenum.data(),fcid.data());

            rarp::NLIST x = db::get_named_row("seqdb",fer2);
            // cout << "We got " << x["affected"] << "\n";
            sprintf(fer2,"select Lane.data->'$.fastq.path.archive' archive from Lane join prepT on Lane.prepid=prepT.prepid join Flowcell on Flowcell.fcid=Lane.fcid where sample_internal_name = '%s' and lanenum = %s and fcillumid='%s'",
              ARSV2[i].sample_internal_name.data(),ARSV2[i].lanenum.data(),fcid.data()
            );
            x.clear(), x = db::get_named_row("seqdb",fer2);

            // cout << "We got " << x["archive"] << " and " << ARSV2[i].archive_dir << "\n";
            
            assert(x["archive"]=='"'+ARSV2[i].archive_dir+'"');

        }

        // do with bother with runparms.xml etc....

    }
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

int do_it(
  lists::NAMES & names,
  rarp::NLISTS & fc,
tinyxml2::XMLElement *& rp,
  std::string const & fcid,
  std::string const & archive_dir,
  std::string const & ssfn,
  std::string const & bcl_chkpnt,
  std::string const & bcl_dir,
  std::string const & bcl,
  std::string const & run_dir_base,
  std::string const & run_dir,
  int argc,
  int userid
);

int bcl(int argc, char **argv) {

    // FDP_CHECK;
    // bored::quickie();
    // integrate check_runs.cpp procedure but have it db driven - i.e. check if flowecell exists and automatically update the flowcell status!?!
    // use the run single command thing - i.e. with md5 for not rerunning of steps like bcl conversion!?!?!?!?

    using namespace std;
    using namespace rarp;

    int userid;
    { NLIST x = db::get_named_row("seqdb","select userid from users where netid='%s'",opts::myuser.user()); assert(x["userid"]!="");  userid=atoi(x["userid"].data()); }
    // { NLIST x = db::get_named_row("seqdb","select userid from users where netid='%s'",getenv("USER")); assert(x["userid"]!="");  userid=atoi(x["userid"].data()); }

    lists::NAMES names;

    // do we use this nonsense anymore?!?
    insane_frac["NSP"]["Exome"]=24,  insane_frac["NSP"]["Genome"]=3,  insane_frac["NSP"]["Custom_Capture"]=240,  insane_frac["NSP"]["RNAseq"]=24;
    insane_frac["NS1"]["Exome"]=48,  insane_frac["NS1"]["Genome"]=6,  insane_frac["NS1"]["Custom_Capture"]=480,  insane_frac["NS1"]["RNAseq"]=48;
    insane_frac["NS2"]["Exome"]=96,  insane_frac["NS2"]["Genome"]=12, insane_frac["NS2"]["Custom_Capture"]=960,  insane_frac["NS2"]["RNAseq"]=96;
    insane_frac["NS4"]["Exome"]=144, insane_frac["NS4"]["Genome"]=16, insane_frac["NS4"]["Custom_Capture"]=1440, insane_frac["NS4"]["RNAseq"]=144;
    
    string run_dir, run_parms, run_dir_base, fcid;

    NLISTS fc;

    #define MSG(X) cout << "USING " << #X << " = '" << X << "'\n"
    MSG(seq::NETAPP_OUT_DIR);
    MSG(seq::SEQ_RUN_DIR);
    MSG(seq::ARCHIVE_DIR);
    #undef MSG

    if(argc==1 || (argc==2 && strcmp(argv[1],"run")==0) ) {

        rarp::LISTS check;
        
        get_table(check,"find %s -maxdepth 2 -name CopyComplete.txt",seq::SEQ_RUN_DIR);

        assert(check.size()>0);
        cout << "there are " << check.size() << " flowcell dirs available\n";

        for (unsigned int i = 0; i < check.size() ; ++i) {

            check[i][0] = check[i][0].substr(0,check[i][0].length()-16)+"RunInfo.xml";
            string fc1 = get_single_line_output_as_string("perl -ne 'm{<Flowcell>(\\S+?)</Flowcell>} && print qq{$1\\n}' %s",check[i][0].data());

            fflush(stdout);

            rarp::NLIST x = db::get_named_row("seqdb",GET_FLOWCELL,fc1.data());
            // for(rarp::NLIST::iterator it=x.begin(); it!=x.end(); it++) cout << "["<<it->first<<"] "<<it->second<<"\n";
            
            if(x["fc_status"]=="registered" && x["fail"]=="0") {
            // if(x["complete"]=="0" && x["fail"]=="0") {
            
                cout << "flowcell '" << fc1 << "' can run\n";
                fcid=fc1;

                if(argc==2 && strcmp(argv[1],"run")==0) {
                    db::get_named_rows("seqdb",fc,GET_FLOWCELL,fcid.data());
                    for(rarp::NLIST::iterator it=x.begin(); it!=x.end(); it++) cout << "["<<it->first<<"] "<<it->second<<"\n";
                    break;
                }

            }
        } 

        if(fcid.empty()) cout << "there's nothing to do\n",exit(0); 
        
        //string msg = "FLOWCELL READY FOR BCL CONVERSION " + fcid;
        //sendmail("nb2975@cumc.columbia.edu","nb2975@cumc.columbia.edu",msg.data(),"As above.");

        if(argc==1 || strcmp(argv[1],"run")!=0) cout << "won't run\n", exit(0);
        // cout << "using " << fcid << "\n";

    }else if(argc==2) {

        fcid=*(argv+1); // bcl = string(NETAPP_OUT_DIR_TMP) + "/BCL/";
        db::get_named_rows("seqdb",fc,GET_FLOWCELL,fcid.data());

    }else cout << "usage: " << *argv << " <fcid>\n",exit(1);

    NLIST * flowcell = 0; // NLIST & flowcell = fc[0];

    { assert(fc.size()==1); for(unsigned int i=0; i<fc.size(); ++i) {
        cout << "-> " << __FILE__ << ":" << __LINE__ << " (" << __func__ << ") : seq.flowcell :\n";
        for(rarp::NLIST::iterator it=fc[i].begin(); it!=fc[i].end(); it++) cout << "["<<it->first<<"] "<<it->second<<"\n";
        cout << "\n";
    } flowcell = &fc[0]; }

    string archive_dir;

    { 
    long long space = atol(get_single_line_output_as_string("df %s | grep -P '\\d+\\s+\\d+\\s+\\d+' | awk '{print $3}'",seq::NETAPP_OUT_DIR).data());
    cout << "we have " << space << " (" << (float)space/(1024*1024*1024) << ") available\n";
    if( (float)space/(1024*1024*1024) < 10.0 ) cout << "insufficient space to run!\n",exit(1);
    else cout << "will run\n";
    }

    cout << " we have fcillumid=" << fcid << " and fcid=" << flowcell->operator[]("FCID") << "\n";

    { if(fc[0]["fc_status"]=="fastq_archived") cout << "this has been processed to completion\n",exit(0); // assert(fc[0]["fc_status"]=="registered");
    char arsv2[1024];
    sprintf(arsv2,"update Flowcell set fc_status = 'sequenced' where FCillumID = '%s' and fc_status = '%s' and fail = 0 ; select row_count() locked",fcid.data(),fc[0]["fc_status"].data());
    NLIST x = db::get_named_row("seqdb",arsv2);
    assert(x["locked"]=="1"); }
    // cout << "we have a lock on " << fcid << "\n"; } 
    
    if(fc[0]["fc_archive"]=="NULL") {

        cout << "updating archive dir\n";
        char arsv2[1024];
        // #define ARCHIVE_DIR "archive/p2018/FASTQ/"
        sprintf(arsv2,"update Flowcell set "
          "complete = 1, seqend = UNIX_TIMESTAMP(NOW()), "
          "fc_archive = '%s' where fcid = '%s' "
          "; select row_count() updated", seq::ARCHIVE_DIR, flowcell->operator[]("FCID").data() 
        ); 

        rarp::NLIST x = db::get_named_row("seqdb",arsv2);
        assert(x["updated"]=="1");

        rarp::NLISTS bored;
        db::get_named_rows("seqdb",bored,"select distinct fcm.seqID, fcm.poolID, fcm.prepID, pt.sample_internal_name, pt.pool_status "
          "from FCmembers fcm join prepT pt on pt.prepID=fcm.prepID where fcm.FCID = %s",(*flowcell)["FCID"].data());
        // cout << "we have " << bored.size() << "\n";
       
        // abolishing LIMS functionality
        for (unsigned int ii = 0 ; ii != bored.size(); ++ii){
            // cout << "["<<ii<<"]\n";
            // for(rarp::NLIST::iterator it=bored[ii].begin(); it!=bored[ii].end(); it++){ cout << "\t"<<it->first <<"\t"<<it->second <<"\n"; }
            // cout << "we modify " << ii << "\n";
            rarp::NLIST g = db::get_named_row("seqdb","update samplesTOrun set complete = 1, s2r_status = 'completed'  where seqID = %s ; select row_count() updated ",bored[ii]["seqID"].data());
            // cout << "\twe updated " << g["updated"] << "\n\n"; 
            assert(g["updated"]=="1");

            g = db::get_named_row("seqdb","insert into statusT set "
              "prepID = %s, status = 'Sequencing Complete', userID = %d, poolID = %s, status_time = UNIX_TIMESTAMP(NOW()), seqID = %s; select row_count() updated ",
              bored[ii]["prepID"].data(), userid, bored[ii]["poolID"].data(), bored[ii]["seqID"].data()
            );
            assert(g["updated"]=="1");
            if(bored[ii]["pool_status"]=="pool") {

                ////get pool members
                rarp::NLISTS bored2;
                db::get_named_rows("seqdb",bored2,"select sample_internal_name, prepID from prepT where pool_status = 'pool_members' and poolID = %s",bored[ii]["poolID"].data());

                for (unsigned int iii = 0 ; iii != bored2.size(); ++iii){

                    // cout << "["<<iii<<"]\n";
                    // for(rarp::NLIST::iterator it=bored2[iii].begin(); it!=bored2[iii].end(); it++){ cout << "\t"<<it->first <<"\t"<<it->second <<"\n"; }
                    g = db::get_named_row("seqdb","insert into statusT set "
                      "prepID = %s, status = 'Sequencing Complete', userID = %d, poolID = %s, status_time = UNIX_TIMESTAMP(NOW()), seqID = %s; select row_count() updated ",
                      bored2[iii]["prepID"].data(), userid, bored[ii]["poolID"].data(), bored[ii]["seqID"].data()
                    );
                    // cout << "\twe inserted " << g["updated"] << "\n\n"; 
                    assert(g["updated"]=="1");
                }

                rarp::NLIST g = db::get_named_row("seqdb","update prepT set status = 'Sequencing Complete', status_time = UNIX_TIMESTAMP(NOW()), "
                  "sequencing_complete = 1, sequencing_complete_time = UNIX_TIMESTAMP(NOW()) "
                  "where prepID = %s; select row_count() updated",bored[ii]["prepID"].data());
                // cout << "\twe updated prept " << g["updated"] << "\n\n"; 
                assert(g["updated"]=="1");
                
            }else{
                rarp::NLIST g = db::get_named_row("seqdb","update prepT set status = 'Sequencing Complete', status_time = UNIX_TIMESTAMP(NOW()), "
                "sequencing_complete = 1, sequencing_complete_time = UNIX_TIMESTAMP(NOW()) "
                "where prepID = %s; select row_count() updated",bored[ii]["prepID"].data());
                assert(g["updated"]=="1");
            }
        }

        archive_dir=seq::ARCHIVE_DIR;
        // #undef ARCHIVE_DIR
    }else{
        archive_dir=fc[0]["fc_archive"];
    }

    // why?!?
    assert( flowcell->find("complete")->second==flowcell->operator[]("complete") 
      && flowcell->operator[]("complete")==(*flowcell)["complete"] );

    if( (*flowcell)["complete"]=="1" ) cout << "-> this flowcell (" << fcid << ") is already completed\n";
    else cout << "-> this flowcell (" << fcid << ") has not been run\n";

    ////// hmm, didn't realise, should use this to clean stuff up in general?!?
    names.add_name("NETAPP_OUT_DIR",seq::NETAPP_OUT_DIR);
    FILL_IN_DIR(bcl,names,"{{NETAPP_OUT_DIR}}////BCL///",1024);
    FILL_IN_DIR_2("BCL_BASE_DIR",names,"{{NETAPP_OUT_DIR}}////BCL///",1024);
    // names.show_names();

    { LISTS check; char argh[1024];
    sprintf(argh,"find %s -type d -maxdepth 1 -name '*%s' 2>/dev/null", seq::SEQ_RUN_DIR, fcid.data() );

    get_table(check,argh);
    if(check.size()==0) cout << "cannot find run dir for " << fcid << " in " << seq::SEQ_RUN_DIR << "\n",exit(1);
    else if(check.size()>1) {
        cout << "run dir for " << fcid << " in " << seq::SEQ_RUN_DIR << " is ambiguous - please clean up partial runs\n";
        for (unsigned int i=0;i<check.size();++i) {
            cout <<"\t["<<i<<"] "<<check[i][0] << "\n";
        }
        exit(1);
    }

    run_dir=check[0][0]; 
    names.add_name("RUN_DIR", check[0][0].data());

    sprintf(argh,"find %s -type f -maxdepth 1 -iname 'runparameters.xml' 2>/dev/null",run_dir.data());

    check.clear();
    get_table(check,argh);
    run_parms=check[0][0]; 
    assert(check.size()==1); }
    // cout << "run params = " << run_parms << "\n\n";
    
    tinyxml2::XMLDocument doc;
    doc.LoadFile(run_parms.data());

    assert(run_parms.length()>=17);
    string machtype = run_parms.substr(run_parms.length()-17,17) == "RunParameters.xml" ?   "NovaSeq"
                    : run_parms.substr(run_parms.length()-17,17) == "runParameters.xml" ?   "HighSeq"
                    :                                                                       "";

    assert(run_parms!="");
    // cout << "-> using " << machtype << "\n\n";
    
    tinyxml2::XMLElement *rp = doc.FirstChildElement( "RunParameters" );
    // const char* title = doc.FirstChildElement( "RunParameters" )->FirstChildElement( "RtaVersion" )->GetText();
    // cout << "NAME="<<rp->FirstChildElement("RunStartDate")->Name()<<", VALUE="<<rp->FirstChildElement("RunStartDate")->GetText()<<"\n";

    //// check it's the same as end of run_dir?!?
    run_dir_base = string(rp->FirstChildElement("RunStartDate")->GetText())
      +"_" + rp->FirstChildElement("InstrumentName")->GetText()
      + rp->FirstChildElement( "Side" )->GetText()
      + "_" + fcid + "_Unaligned/";

    string ssfn = run_dir+"/Sample_Sheet_"+fcid+"_v1.0.csv";
    names.add_name("SS",ssfn.data());
    names.add_name("MM","1");
    names.add_name("BCL_OUT_DIR",(string(bcl) + "/" + run_dir_base).data());
    // +rp->FirstChildElement("RfidsInfo")->FirstChildElement( "FlowCellSerialBarcode" )->GetText();

    string bcl_dir = bcl+run_dir_base, bcl_chkpnt = bcl_dir+"bcl_complete";

    print_info(rp);

    // const char* title = doc.FirstChildElement( "RunParameters" )->FirstChildElement( "RtaVersion" )->GetText();
    if(isregfile((run_dir+"/RTAComplete.txt").data())) cout << "run finished\n";
    else cout << "run doesn't seem to have finished\n",exit(1);
    if(isregfile((run_dir+"/CopyComplete.txt").data())) cout << "copy finished\n";
    else cout << "copy doesn't seem to have finished : " << run_dir << "/CopyComplete.txt" ,exit(1);

    //sendmail("nb2975@cumc.columbia.edu","nb2975@cumc.columbia.edu",string("Running FC "+fcid).data(),"That's all.");

    // int run_it_in_fork(int argc, char ** argv){
    // cout << "must relocate lock to here!?!\n";

    struct sigaction act;
    memset (&act, '\0', sizeof(act));  
    act.sa_sigaction = sig_action_function;
    act.sa_flags = SA_SIGINFO;    
    sigaction(SIGUSR1, &act, 0); // sigaction(SIGALRM, &act, 0);  ///r could use pipe for IPC?!? - not really much point?!?
    sigaction(SIGUSR2, &act, 0); // sigaction(SIGALRM, &act, 0);
    struct sigaction act2;
    memset (&act2, '\0', sizeof(act2)); // wasn't initialising?!?
    act2.sa_sigaction = crit_err_hdlr;
    act2.sa_flags = SA_RESTART | SA_SIGINFO;
    if (sigaction(SIGSEGV, &act2, (struct sigaction *)NULL) != 0) {
        fprintf(stderr, "error setting signal handler for %d (%s)\n", SIGSEGV, strsignal(SIGSEGV));
        exit(1);
    }

    pid_t child_pid; // , wpid;
    int status = 0;
    while(1) { // for(;;) 
        child_pid = fork();
        switch (child_pid) { // if (child_pid==-1) {
            case -1:    std::cerr << "unable to fork!?! Exiting."<<endl;    exit(1);        break;
            case 0: { 
                cerr << "running within child process " << getpid() << "\n";

                do_it(names, fc, rp, fcid, archive_dir, ssfn, bcl_chkpnt, bcl_dir, bcl, run_dir_base, run_dir, argc, userid);

                // why?!? - why not just check status below and exit?!?
                // return 0;

                static char messageText[] = "all done";
                union sigval signal_value;
                signal_value.sival_ptr = messageText;
                cerr << "sending shut-down message to parent\n";
                sigqueue(getppid(), SIGUSR2, signal_value); // sigqueue(getppid(), SIGALRM, signal_value);
                sleep(1);
                return 0;
            } 
            default:        std::cerr << "generated child fork pid= " << child_pid << "\n";  break;
        }

        while ((wait(&status)) > 0) { // while ((wpid = wait(&status)) > 0) {

            if (status>0) cerr << "SEEMS WE EXCITED IN A MOST UNDIGNIFIED MANNER - SIGNAL : " << status << "\n";
            else {
                std::cerr << "child seems to have finished - should exit here?!?\n";
                // break;
            }

            exit(1);

            if(1) {
                cout << "not re-runnign right now\n";
                break;
            }

            if( restart_pos >= sizeof(restart_int)/sizeof(int) ) {
                cerr << "will not attempt again\n";
                break;
            }

            cerr << "[" << restart_pos << "] will restart in " << restart_int[restart_pos] << "\n";
            sleep(restart_int[restart_pos]);
            ++restart_pos;

        }
    }
    return 0;
// }
}

int do_it(
  lists::NAMES & names,
  rarp::NLISTS & fc,
tinyxml2::XMLElement *& rp,
  std::string const & fcid,
  std::string const & archive_dir,
  std::string const & ssfn,
  std::string const & bcl_chkpnt,
  std::string const & bcl_dir,
  std::string const & bcl,
  std::string const & run_dir_base,
  std::string const & run_dir,
  int /* argc */,
  int userid
) {

    using namespace std;
    using namespace rarp;

    NLISTS ns;

    db::get_named_rows("seqdb",ns,"select count(distinct lanenum) lane_num, fc_insert, fc_update, fc_status, group_concat(distinct lanenum separator '_') lanes, "
      "group_concat(distinct p.sample_internal_name separator ':') samples, group_concat(distinct gafbin separator ',') projects, "
      "machine instrument,complete,fail error_state, lenr1,lenr2,leni1,leni2,chemver, replace(name,' ','') fullname "
      "from SampleT s " "join Experiment e on s.sample_id=e.sample_id " "join prepT p on e.id=p.experiment_id "
      "join Lane l on p.prepid=l.prepid " "join Flowcell f on l.fcid=f.fcid "
      "join users u on f.userid=u.userid where fcillumid = '%s'",fcid.data()
    );

    // for(rarp::NLIST::iterator it=ns[0].begin(); it!=ns[0].end(); it++) { cout << "\t[" << it->first << "] "<< it->second << "\n"; }

    assert(ns.size()==1);
    cout << "\nlims machine = " << ns[0]["instrument"] << ", " << ns[0]["complete"] << ", " << ns[0]["error_state"]<<"\n";
    assert(ns[0]["fc_status"]=="sequenced"); assert(ns[0]["error_state"]=="0");
    cout << "\n";

    update_fc_table(run_dir, fcid, archive_dir, rp);

    interop_bits(run_dir, fcid, ns);

    { rarp::NLISTS rgs;
    db::get_named_rows("seqdb",rgs,"select lanenum, l.prepid, p.sample_type, upper(p.sample_type) seq_type_upper, p.sample_internal_name, lnfraction, lnfractionact, lnyield, data, rg_status, "
      "sequence adapter, pool.name pool_name, p.exomekit capture_kit, fcillumid, "
      "round(sr.picomoles,1) silly_desc, l.fcid, "
      "replace(e.gafbin,' ','') project "
      "from prepT p join Lane l on p.prepid=l.prepid join Flowcell f on l.fcid=f.fcid "
      "join adapter on adapter.id=p.adapter_id "
      "join pool on pool.id=p.poolid "
      "join samplesTOrun sr on sr.seqid=l.seqid " "join Experiment e on p.experiment_id=e.id " "join SampleT s on s.sample_id=e.sample_id "
      "where fcillumid = '%s' order by lanenum, prepid",
      fcid.data()
    ); // ,l+1);
    
    // cout << "there are " << rgs.size() << " read groups\n";
    assert(rgs.size()>MIN_RGS);

    // cout << "SHOULDN'T RE-WRITE SAMPLE SHEET!?!\n";
    // assert(ns[0]["leni2"]=="0");
    write_sample_sheet(fcid, ssfn, ns, rp, rgs, (ns[0]["leni2"]=="0"));

    { char blarp[16*1024];
    ///// submit as script?!?
    names.fill_in_name("#!/bin/bash\nexport LD_LIBRARY_PATH=/nfs/goldstein/software/bcl2fastq2-v2.20-x86_64/lib:"
      "/nfs/goldstein/software/libxml2-2.9.4-x86_64/lib:"
      "/nfs/goldstein/software/binutils-2.26.1/lib:"
      "/nfs/goldstein/software/boost_1_54_0/lib:"
      "/nfs/goldstein/software/gcc-4.9.3/lib64:"
      "$LD_LIBRARY_PATH\n"
      // "if [ -e \"{{BCL_OUT_DIR}}/bcl_complete\" ]; then\necho '
      "nohup " // wt?!?
      "/nfs/goldstein/software/bcl2fastq2-v2.20-x86_64/bin/bcl2fastq "
      "--runfolder-dir {{RUN_DIR}} "
      "--output-dir {{BCL_OUT_DIR}} "
      "--sample-sheet {{SS}} "
      "--barcode-mismatches {{MM}} 2>&1 | tee {{RUN_DIR}}/bcl_conversion_log.txt "
      " && touch {{BCL_OUT_DIR}}/bcl_complete\n",
      // "if [ $? -eq 0 ]; then\n touch {{BCL_OUT_DIR}}/bcl_complete\nelse\n touch {{BCL_OUT_DIR}}/bcl_error\nfi\n"
      blarp, sizeof(blarp));

    // cout << "BCL2FASTQ CMD:\n" << blarp << "\n";
    { ofstream arsv2((run_dir+"/bcl2fastq.sh").data());
    if(!arsv2) cout << "unable to write script\n",exit(1);
    arsv2 << blarp;
    arsv2.close();
    // cout << "using " << (run_dir+"/bcl2fastq.sh").data() << "\n";
    chmod((run_dir+"/bcl2fastq.sh").data(),00744);
    }

    if(isregfile(bcl_chkpnt.data())) {
        cout << "bcl conversion already run - skipping\n";
    }else{
        set_status(rgs,"BCL_Started",userid);
        // bus error happens no matter what when returning control from system call 
        // errrr, that sounds rather like an inode mem alignment error!?!
        // thread out the archiving to 8-12 threads?!?
        if(system((run_dir+"/bcl2fastq.sh").data())!=0) cout << "there was a problem\n",exit(1);
        set_status(rgs,"BCL_Completed",userid);
    } }

    { char arsv2[1024];
    sprintf(arsv2,"update Flowcell set fc_status = 'fastq_converted' where FCillumID = '%s' and fc_status = '%s' and fail = 0 ; select row_count() locked",fcid.data(),"sequenced");
    NLIST x = db::get_named_row("seqdb",arsv2); assert(x["locked"]=="1"); cout << "we have a lock on " << fcid << "\n"; }

    string html_report = bcl + run_dir_base + "/Reports/html/" + fcid + "/all/all/all/laneBarcode.html";
    rarp::NLISTS hs;
    html::get_named_rows_from_fc_summary(hs,html_report.data());
    int ln_yield=0;
    // float ln_yield=0.0;

    update_actual_lane_fractions_from_html(hs,rgs,fcid,ln_yield,fc[0]["fcYield"].data());

    // cout << "FINAL YIELD = " << ln_yield << " (" << fc[0]["fcYield"] << ")\n";
    if(ln_yield!=atoi(fc[0]["fcYield"].data())) {
        char argh[1024];
        sprintf(argh,"update Flowcell set fcyield = %d, casavaver = '%s' where FCillumid = '%s'; select row_count() updated",ln_yield,"v2.20.0.422",fcid.data());
        rarp::NLIST x = db::get_named_row("seqdb",argh);
        // cout << "using " << argh << "\n";
        assert(x["updated"]=="1");
    }else{
        if(atoi(fc[0]["fcYield"].data())!=ln_yield) cout << "wt?!?\n";
        cout << "total yield alaredy set\n";
        // assert(atoi(fc[0]["fcYield"].data())==ln_yield);
    }

    move_fastq_to_sample_specific_dirs(hs,rgs,bcl_dir,bcl,seq::NETAPP_OUT_DIR,run_dir_base,fcid,archive_dir);

    { string bored = archive_dir + "/summary/SAV/" + fcid + "_" + ns[0]["instrument"] + "_SAV.tar.gz";
    char argh[8*1024];
    sprintf(argh,"tar czvf %s %s/RunInfo.xml %s/RunParameters.xml %s/InterOp ",bored.data(),run_dir.data(),run_dir.data(),run_dir.data());
    // cout << "using\n"<<argh<<"\n";

    if(!isregfile(bored.data())) {
        if(system(argh)!=0) cout << "oh, what on earth...\n", exit(1);
    }else { 
        cout << "interop already archived... " << bored << "\n"; 
    } }

    move_fastq_to_sample_specific_dirs(hs,rgs,bcl_dir,bcl,seq::NETAPP_OUT_DIR,run_dir_base,fcid,archive_dir,true);

    { char arsv2[1024];
    sprintf(arsv2,"update Flowcell set fc_status = 'fastq_archived' where FCillumID = '%s' and fc_status = '%s' and fail = 0; select row_count() locked",fcid.data(),"fastq_converted");
    NLIST x = db::get_named_row("seqdb",arsv2); assert(x["locked"]=="1"); cout << "we have a lock on " << fcid << "\n"; } 

    cout << "bye\n";

    return 0; }

}

void submit_and_post_checks(bool);

namespace post_pipe {

static char const * list[] = { "bam", "realn.recal.bam", "realn.recal.bai", "analysisReady.annotated.vcf.gz", "analysisReady.annotated.vcf.gz.tbi", "g.vcf.gz", "g.vcf.gz.tbi", "coverage_bins", "pipeline_data.tar.gz" };

using namespace std;

void /* int */ check_pipeline_output(int argc , char ** argv) {
// void /* int */ check_pipeline_output(int /* argc */, char ** /* argv */) {
    // cout << "check_pipeline_output\n";
    --argc;
    ++argv;
    // this was extreme paranoia - had to implement checks in major hurry so double implemented and never removed...
    // why aren't we invoking directly as with pipeline?!?
    // cout << "run md5sum\n";
    submit_and_post_checks(true);

    /* 
    cleanup_pipeline_scratch(argc,argv); 
    char const * script = "PIPELINE_PREWIPE_10-20.pl";
    assert(isregfile(script));
    char cmd[1024];
    sprintf(cmd,"/usr/bin/perl %s",script);
    if(system(cmd)) std::cerr << "there's a problem\n", exit(1);
    */
    return;

}

// Timey StartTime(time(0));
// cout << "starting at " << StartTime.epoch_time_as_string(tt1) << ", " << StartTime.iso_time(tt2) << "\n";
// memset(timeyt,sizeof timeyt, 0); memcpy(timeyt,z+1,i-1);
// Timey O2(timeyt,"%a %b %d %H:%M:%S %Y",1);   
// int min = (StartTime.epoch_time_as_time_t()-O2.epoch_time_as_time_t())/60;

// this was separate as everything was done in a major hurry so was double implemented for caution and it requires sudo rules that were only avaialable on a single host
int cleanup_pipeline_scratch(int argc, char **argv) {
    
    using namespace std;

    // {char argh[1024]; gethostname(argh,sizeof(argh)); cout << ">>>>>>>>>> " << argh << " : " << *argv << " : " << time(0) << " : "; }

    rarp::NLISTS stuff;
    if(argc==1) {
        db::get_named_rows("seqdb",stuff,"select m.*,q.AlignSeqFileLoc FROM dragen_sample_metadata m, dragen_qc_metrics q "
          " where m.pseudo_prepid=q.pseudo_prepid and m.is_merged = 20 order by seqscratch_drive desc  limit 1");
    }else{
        db::get_named_rows("seqdb",stuff,"select m.*,q.AlignSeqFileLoc FROM dragen_sample_metadata m, dragen_qc_metrics q "
          " where m.pseudo_prepid=q.pseudo_prepid and m.pseudo_prepid = %d and is_merged >= 20 and is_merged <= 40 ",atoi(argv[1]));
    }

    if(stuff.size()==0) return 1; 

    rarp::NLIST & entry = stuff[0];
    if(argc==1){
    bool locked = query::silly_update("update dragen_sample_metadata set is_merged = 32 where is_merged = 20 and pseudo_prepid = %s",entry["pseudo_prepid"].data());
    assert(locked);
    }else query::silly_update("update dragen_sample_metadata set is_merged = 32 where pseudo_prepid = %s",entry["pseudo_prepid"].data());

    // for(rarp::NLIST::iterator i=stuff[0].begin(); i!=stuff[0].end(); ++i) { cout << " ["<<i->first<<"]= "<<i->second<<"\n";  }

    string type = entry["sample_type"];
    for(unsigned i=0; i<entry["sample_type"].length(); ++i) type[i]=::toupper(type[i]);
    string u_name = entry["sample_name"]+"."+entry["pseudo_prepid"],
      sdir = "/nfs/"+entry["seqscratch_drive"]+"/ALIGNMENT/BUILD37/DRAGEN/"+type+"/"+u_name+"/",
      adir = entry["AlignSeqFileLoc"]+"/"+u_name+"/",
      md5s = adir+"/md5sums.txt";

    // do we care?!?
    string who = sdir+"who.txt", sge = sdir+"sge_wrapper.log";

    cout << "\nwe have scratch dir = " << sdir << "\nwe have archive dir = " << adir << "\nmd5s = " << md5s << "\n\n";

    assert((isdir(sdir.data())));
    assert((isdir(adir.data())));
    assert((isregfile(md5s.data())));

    FILE * mf = fopen(md5s.data(),"r");
    assert(mf);
    char line[2048];
    // int lc=0;

    std::vector<string> files_to_protect;

    char const ** thing = entry["sample_type"]=="Genome_As_Fake_Exome" ? list : list+1;
    // char * lp;
    while(fgets(line,2048,mf)){
        char a_md5[1024],file[1024];
        memset(a_md5,0,sizeof(a_md5)); // can't be bothered to null-term later?!?
        memset(file,0,sizeof(file));
        // cout << "["<<lc++<<"] "<<line << " vs " << *thing<< "\n";
        int l = strlen(line);
        char * fn = 0;
        for(int i=0;i<l;++i){
            if(line[i++]==' '&&line[i]==' '){
                // cout << "have it at " << i << "\n";
                memcpy(a_md5,line,i-1);
                memcpy(file,line+i+1,l-i-2);
                fn = strrchr(file,'/');
                assert(fn);
                ++fn;
            }
        }

        files_to_protect.push_back(file);

        string tmp = u_name+"."+*thing++,s_file = sdir+fn;

        // cout << "a_md5= '" << a_md5 << "'\nfile       = '" << file << "'\nfn         = '" << fn 
            // << "'\n= '"<< tmp  << "'\ncheck2=    = '" << s_file << "'\n";

        assert(tmp==fn);

        char const * x = "md5sum %s | awk '{print $1}'";
        // Yum v(x);
        string s_md5 = Yum<char const *>(x)(s_file.data()); 
        
        // cout << "s_md5 = '"<< s_md5 << "'\n\n";
        assert(s_md5==a_md5);
        int a_size=filesize(file), s_size=filesize(s_file.data());
        // cout << "a_file " << file << " size = " << a_size << "\n";
        // cout << "s_file " << s_file << " size = " << s_size << "\n";
        assert(a_size==s_size);

    }

    query::silly_update("update dragen_sample_metadata set is_merged = 31 where is_merged = 32 and pseudo_prepid = %s",entry["pseudo_prepid"].data());

    { files_to_protect.push_back(adir+"md5sums.txt");
    files_to_protect.push_back(adir+"coverage.tar.gz");
    files_to_protect.push_back(adir+"gq.tar.gz");

    for(unsigned int o=0;o<files_to_protect.size();o++){
    // for(int o=0;o<(int)files_to_protect.size();o++){
        cout << "["<<o<<"] "<<files_to_protect[o]<<"\n"; 
        Lazy("sudo /bin/chown dragen.analysts %s",files_to_protect[o].data());
        Lazy("sudo /bin/chmod 440 %s",files_to_protect[o].data());
    }
    Lazy("sudo /bin/chown dragen.analysts %s",adir.data());
    // this was me being too nervous to change anything run on 80k samples
    // not sure if doing this...?!?
    Lazy("sudo /bin/chown dragen.analysts %s",adir.data());
    // Lazy("sudo /bin/chgrp -R analysts %s",adir.data());
    Lazy("sudo /bin/chmod 550 %s",adir.data()); }

    { std::vector<string> REMOVE_LIST;
    // cout << "mofo : (a)" << sdir << "\n";
    Yum<char const *>("find %s -type f")(REMOVE_LIST,sdir.data());
    // cout << "mofo : (b)" << REMOVE_LIST.size() << "\n";
    // for(int o=REMOVE_LIST.size()-1;o>0;--o){
    for(int o=REMOVE_LIST.size()-1;o>=0;--o){
        // cout << "we wipe ["<<o<<"] " << REMOVE_LIST[o] << "\n";
        assert(strchr(REMOVE_LIST[o].data(),' ')==NULL);
        assert(strchr(REMOVE_LIST[o].data(),'*')==NULL);
        int ret = unlink(REMOVE_LIST[o].data());
        if(ret!=0) cout << "we got = " << ret << "\n";
        // if(ret!=0) Lazy("echo '%s' >> /nfs/seqscratch09/dsth/play/could_not_wipe.txt",REMOVE_LIST[o].data());
    } }

    { std::vector<string> REMOVE_LIST;
    // cout << "mofo : (c)" << sdir << "\n";
    Yum<char const *>("find %s -type d")(REMOVE_LIST,sdir.data());
    // cout << "mofo : (d)" << REMOVE_LIST.size() << "\n";
    for(int o=REMOVE_LIST.size()-1;o>=0;--o){
        // cout << "we wipe ["<<o<<"] ";
        cout.flush();
        cout << REMOVE_LIST[o] << "\n";
        assert(strchr(REMOVE_LIST[o].data(),' ')==NULL);
        assert(strchr(REMOVE_LIST[o].data(),'*')==NULL);
        int ret = rmdir(REMOVE_LIST[o].data());
        if(ret!=0) cout << "we got = " << ret << "\n";
        // if(ret!=0) Lazy("echo '%s' >> /nfs/seqscratch09/dsth/play/could_not_wipe.txt",REMOVE_LIST[o].data());
    } }

    query::silly_update("update dragen_sample_metadata set is_merged = 30 where is_merged = 31 and pseudo_prepid = %s",entry["pseudo_prepid"].data());

    return 0;
}

} // post_pipe

#define RSYNC "rsync -azpvv -L --no-owner --no-perms "
#define RSYNC_CLEAN "rsync -azpvv -L --no-owner --no-perms --remove-source-files "

bool bam_check(char const *, rarp::NLIST &);

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

void metrics(int /* argc */, char ** /* argv */) {

    using namespace std;
    using namespace rarp;

    NLISTS fc;
    db::get_named_rows("seqdb",fc,"select id, rg_status,data->'$.bam' data ,concat(f.fcillumid,'.',l.lanenum) PU, lanenum, l.fcid, l.prepid, sample_internal_name, sample_type, exomekit "
        " from Lane l "
        " join Flowcell f on l.fcid=f.fcid "
        " join prepT p on l.prepid=p.prepid "
        " where data->'$.bam.path.scratch' is not null "
        " and rg_status in ('fastq_mapped','fastq_archiving','fastq_archived') "
        "  and rg_metrics_status = 'pending' "
        " order by rg_insert asc limit 1 ");

    if(fc.size()==0) {
        // cout << "There's nothing to do metrics for. Bye.\n";
        return;
        // exit(0);
    }

    assert(fc.size()==1);
    NLIST & x = fc[0];
    for(NLIST::iterator it = x.begin(); it != x.end(); it++ ) cout << "["<<it->first<<"] "<<it->second<<"\n";

    rapidjson::Document doc; 
    // this really shouldn't be legal?!? : char buffer[x["data"].length()];
    char buffer[16*1024]; 
    memset(buffer,0,sizeof(buffer));
    memcpy(buffer, x["data"].data(), x["data"].length());

    assert(x["data"]!="NULL");
    cout << "parsing " << x["data"] << "\n";
    if (doc.Parse(buffer).HasParseError()) {
    // if (doc.ParseInsitu(buffer).HasParseError()) {
        cout << "there's an issue with the json\n";
        exit(1);
    }

    { char arsv2[2048]; sprintf(arsv2,"update Lane set rg_metrics_status = 'running' where rg_metrics_status = 'pending' and prepid = %s and lanenum = %s and fcid = '%s' ; select row_count() updated",
x["prepid"].data(),x["lanenum"].data(),x["fcid"].data());
NLIST u = db::get_named_row("seqdb",arsv2); assert(u["updated"]=="1"); }

    assert(doc.FindMember("path")->value.HasMember("scratch"));
    assert(doc.HasMember("basename"));
    string base_metrics = doc.FindMember("path")->value.FindMember("scratch")->value.GetString();
    base_metrics += "/";
    base_metrics += doc.FindMember("basename")->value.GetString();
    assert(base_metrics.substr(base_metrics.length()-4,4)==".bam");
    string bam = base_metrics;
    base_metrics=base_metrics.substr(0,base_metrics.length()-4);
    base_metrics += ".mapping_metrics.csv";
    assert(isregfile(base_metrics.data()));

    float dups=-0.0, map=-0.0;
    long long useable=0, tot=0;

    { char const * TOT = "MAPPING/ALIGNING SUMMARY,,Total input reads,", * DUP ="MAPPING/ALIGNING SUMMARY,,Number of duplicate reads (marked),",
        * MAP = "MAPPING/ALIGNING SUMMARY,,Mapped reads,", * USE = "MAPPING/ALIGNING SUMMARY,,Number of unique & mapped reads (excl. dups),";
    char info[2048];
    FILE * file_to_check = fopen(base_metrics.data(),"r");
    assert(file_to_check);
    while(fgets(info,2048,file_to_check)) {
        if(memcmp(DUP,info,strlen(DUP))==0) {
            dups = atof(strrchr(info,',')+1);
        }else if(memcmp(MAP,info,strlen(MAP))==0) {
            map = atof(strrchr(info,',')+1);
        }else if(memcmp(TOT,info,strlen(TOT))==0) {
            tot = atol(strrchr(info,',')+1);
        }else if(memcmp(USE,info,strlen(USE))==0) {
            vector<string> x; tokenise(x,info,',');
            useable = atol(x[3].data());
        }
    }
    fclose(file_to_check); }

    string bored1 = base_metrics + ".target", bored2 = base_metrics + ".contam", bed;
    if(x["sample_type"]=="Genome") bed="/nfs/seqscratch_ssd/PIPELINE_DATA/ccds_regions.bed";
    else{
        cout << "need to grab the bed file!?!?\n";
        NLIST quero_chorar = db::get_named_row("seqdb","select count(1) count,region_file_lsrc bed from captureKit where prepT_name = '%s' and chr = 'all'",x["exomekit"].data());
        for(NLIST::iterator it = quero_chorar.begin(); it != quero_chorar.end(); it++ ) { cout << "["<<it->first<<"] "<<it->second<<"\n"; }
        assert(quero_chorar["count"]=="1");
        bed=quero_chorar["bed"];
    }

    cout << "USING " << bed << "\n";
    char bored[16*1024];
    //// if genome use ccds bed else use actual bed...?!?
    snprintf(bored,16*1024,ALIGNSTATS " -q 10 -i %s -t %s "
        " -o %s >/dev/null 2>&1 & "
        VERIFYBAMID " --vcf /nfs/seqscratch_ssd/PIPELINE_DATA/All.1kwgs.genomes.seqcap.exome.vcf.gz "
" --bam %s --ignoreRG --maxDepth 250 --precise --out %s >/dev/null 2>&1 & wait",
        bam.data(),bed.data(),
        bored1.data(), // base_metrics.data(),
        bam.data(), bored2.data()); // base_metrics.data());

    bored2+=".selfSM";

cout << "USING " << bored1 << "\n";
cout << "USING " << bored2 << "\n";
// cout << "TMP:\n"<<bored<<"\n";

/////// need to make sure the file isn't empty too!?!?
    if(isregfile(bored1.data())&&isregfile(bored2.data())) cout << "re-use files\n"; 
    else{
        // cout << "need to run\n";
        // cout << "using " << bored << "\n";
        if(system(bored)==0&&isregfile(bored1.data())&&isregfile(bored2.data())) { 
        // if(system(bored)==0) {
            cout << "YAY\n";
        }else{
            // separate from stalled!?!?
            char arsv2[2048]; 
            sprintf(arsv2,"update Lane set rg_metrics_status = 'error' where rg_metrics_status = 'none' and prepid = %s and lanenum = %s and fcid = '%s' ; select row_count() updated",
                x["prepid"].data(),x["lanenum"].data(),x["fcid"].data());
            exit(1);
        }
    }

    // cout << "dups="<<dups<<" map="<<map<<", useable="<<useable<<", tot="<<tot<<"\n";

    float ontarg=atof(get_single_line_output_as_string("grep CapReadsOnTargetOrBufferPct %s | awk '{print $2}'",bored1.data()).data()),
        soft=atof(get_single_line_output_as_string("grep \"SoftClippedReadsPct\" %s | awk '{print $2}'",bored1.data()).data()),
        capmean,capmedian;

    if(x["sample_type"]=="Genome") {
        capmean=atof(get_single_line_output_as_string("grep \"WgsCoverageMean\" %s | awk '{print $2}'",bored1.data()).data());
        capmedian=atof(get_single_line_output_as_string("grep \"WgsCoverageMedian\" %s | awk '{print $2}'",bored1.data()).data());
    }else{
        capmean=atof(get_single_line_output_as_string("grep \"CapCoverageMean\" %s | awk '{print $2}'",bored1.data()).data());
        capmedian=atof(get_single_line_output_as_string("grep \"CapCoverageMedian\" %s | awk '{print $2}'",bored1.data()).data());
    }

    string cs=get_single_line_output_as_string("cut -f7 %s | perl -pe 's{\n}{\t}'",bored2.data());
    cout << "soft="<<soft << ", ontarg="<< ontarg<<"\n";
    assert(memcmp(cs.data(),"FREEMIX",7)==0);
    float contam=atof(cs.substr(8).data());
    cout << "contam="<<contam<<"\n";

    if(contam<0.09){ 
        char arsv2[2048]; 
        sprintf(arsv2,"update Lane set rg_metrics_status = 'okay', "
            " rg_metrics_dups = %f, rg_metrics_contamination = %f, rg_metrics_ontarget = %f, rg_metrics_softclip = %f, "
            " rg_metrics_totalreads = %lld, rg_metrics_uniquemappedreads = %lld, "  
            " rg_metrics_capturemean = %f, rg_metrics_capturemedian = %f "
            " where prepid = %s and lanenum = %s and fcid = '%s' ; select row_count() updated",
            dups, contam, ontarg, soft, tot, useable,
            capmean, capmedian, 
            x["prepid"].data(),x["lanenum"].data(),x["fcid"].data()
        );
        NLIST u = db::get_named_row("seqdb",arsv2); assert(u["updated"]=="1"); 
    }else{
        char arsv2[2048]; 
        for(NLIST::iterator it = x.begin(); it != x.end(); it++ ) cout << "\t["<<it->first<<"] "<<it->second<<"\n\n";
        sprintf(arsv2,"update Lane set rg_metrics_status = 'contaminated', "
            " rg_metrics_dups = %f, rg_metrics_contamination = %f, rg_metrics_ontarget = %f, rg_metrics_softclip = %f, "
            " rg_metrics_totalreads = %lld, rg_metrics_uniquemappedreads = %lld, "
            " rg_metrics_capturemean = %f, rg_metrics_capturemedian = %f "
            " where prepid = %s and lanenum = %s and fcid = '%s' ; select row_count() updated",
            dups, contam, ontarg, soft, tot, useable,
            capmean, capmedian, 
            x["prepid"].data(),x["lanenum"].data(),x["fcid"].data()
        );
        cout << "using " << arsv2 << "\n";

        { NLIST u = db::get_named_row("seqdb",arsv2); 
        for(NLIST::iterator it = u.begin(); it != u.end(); it++ ) cout << "\t["<<it->first<<"] "<<it->second<<"\n";
        assert(u["updated"]=="1"); }

        sprintf(arsv2,"update prepT p set failedprep = 3, status = 'Contaminated' where prepid = %s ; select row_count() updated",x["prepid"].data());  
        cout << "using " << arsv2 << "\n";
        // u.clear();

        { NLIST u = db::get_named_row("seqdb",arsv2); 
        for(NLIST::iterator it = u.begin(); it != u.end(); it++ ) cout << "\t["<<it->first<<"] "<<it->second<<"\n";
        assert(u["updated"]=="1"); }
    }

    cout << "\nbye\n\n"; 

    return;
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

inline bool rsync_checksum_return(std::string const & f, int m, long long s) {
    using namespace std;
    if(!isregfile(f.data())) cout << "file " << f << " is missing\n";
    struct stat st;
    stat(f.data(),&st);
    if(st.st_mtime==m) cout << f << " modification matches ("<<m<<")\n";
    else {
        cout << f << " modification doesn't match ("<<m<<")\n";
        return false;
    }
    if(st.st_size==s) cout << f << " size matches ("<<s<<")\n";
    else {
        cout << f << " size doesn't match ("<<s<<")\n";
        return false;
    }
    cout << "seems okay\n";
    return true;
}

void protect(int /* argc */, char ** /* argv */) {

    using namespace std;
    using namespace rarp;

    // this is a quick hack to get around a sql issue that haven't had time or need to get around
    // trivial way to bypass issue would be to make the queries fc specific
    { NLIST check = db::get_named_row("seqdb","select count(1) fastq_ready_count from Lane where rg_status in ('fastq_ready','fastq_mapped','fastq_mapping_d1','fastq_mapping_d2','mapping_error_d1','mapping_error_d2')");
    assert(check.count("fastq_ready_count"));
    cout << "we have " << check["fastq_ready_count"] << "\n";
    if(check["fastq_ready_count"]!="0") {
        cout << "won't run yet\n";
        return;
    } else cout << "will run\n"; }

    /* 
    char hostname[1024];
    gethostname(hostname,1024);
    if(memcmp(hostname,"igm-atav",8)!=0 && memcmp(hostname,"atav",4)!=0) {
        cout << "run me from atav machine\n";
        return;
        //exit(1);
    } */

    NLISTS fc;
    get_se_group_by_uniform_status(fc,"fastq_archiving");
    for (unsigned int h=0;h<fc.size();++h) {
        cout << "["<<h<<"]\n";
        assert(fc[h]["rg_status"]=="fastq_archiving");
        for(NLIST::iterator it = fc[h].begin(); it != fc[h].end(); it++ ) {
            cout << "\t["<<it->first<<"] "<<it->second<<"\n";
        }
    }

    // = get_an_se_by_status("fastq_archiving");
    // NLIST x = get_an_se_by_status("fastq_archiving");
    if(fc.size()==0) { // if(x.size()==0) {
        return;
        // exit(0);
    }

    string adir;

for (unsigned int rgn=0;rgn<fc.size();++rgn) { // for (int rgn=0;rgn<(int)fc.size();++rgn) {

    /// so tired just gonna alias?!?
    rarp::NLIST & x=fc[rgn];

    cout << "\n\n----------------\nwe will process " << x["data"] << "\n";
    for(NLIST::iterator it = x.begin(); it != x.end(); it++ ) cout << "["<<it->first<<"] "<<it->second<<"\n";

    assert(x["data"]!="NULL");
    rapidjson::Document doc; 
    char buffer[16*1024]; 

    //////// NO WONDER IT'S NON-DETERMINISTIC!?!?!?!?!
    memset(buffer,0,sizeof(buffer));
    memcpy(buffer, x["data"].data(), x["data"].length());

    // this is clearly not a good idea - 'could' allow pulling of fastq_archiving/fastq_archived and ignore the latter if it's an issue but really need to investigate!?!
    rapidjson::ParseResult arsv3;
    if ( ! ( arsv3 = doc.Parse(buffer) ) ) {

        fprintf(stderr, "JSON parse error: %s (%zu)", GetParseError_En(arsv3.Code()), arsv3.Offset());

        cout << "there's an issue with the json\n";
        // cout << "INCREDIBLY NASTY\n";
        if (doc.Parse(buffer).HasParseError()) {
            cout << "still f'd\n";
            exit(1);
        }
    }

    { cout << "attempting lock\n";
    char arsv2[8*1024]; sprintf(arsv2,"update Lane set rg_status = 'protect_lock' where rg_status = 'fastq_archiving' and prepid = %s and lanenum = %s and fcid = '%s' ; select row_count() updated",x["prepid"].data(),x["lanenum"].data(),x["fcid"].data());
    NLIST u = db::get_named_row("seqdb",arsv2); assert(u["updated"]=="1"); }

    if(doc.FindMember("fastq")->value.FindMember("path")->value.HasMember("archive")) {
        cout << "right and proper\n";
    }else{
        char arsv2[8*1024]; sprintf(arsv2,"update Lane set rg_status = 'protect_error' where rg_status = 'protect_lock' and prepid = %s and fcid = '%s' ; select row_count() updated",x["prepid"].data(),x["fcid"].data());
        NLIST u = db::get_named_row("seqdb",arsv2); 
        // cout << "we got " << u["updated"] << " from '" << arsv2 << "'\n";
        assert(u["updated"]=="1");
        exit(1);
    }

    if(adir.empty()) adir = doc.FindMember("fastq")->value.FindMember("path")->value.FindMember("archive")->value.GetString();
    else assert(adir==doc.FindMember("fastq")->value.FindMember("path")->value.FindMember("archive")->value.GetString());

    assert(doc.FindMember("fastq")->value.FindMember("path")->value.HasMember("scratch"));
    assert(doc.FindMember("fastq")->value.HasMember("type"));
    assert(doc.FindMember("fastq")->value.FindMember("r1")->value.HasMember("size"));
    assert(doc.FindMember("fastq")->value.FindMember("r2")->value.HasMember("size"));

    // cout << "checking " << doc.FindMember("fastq")->value.FindMember("type")->value.GetString() << "\n";

    assert(strcmp(doc.FindMember("fastq")->value.FindMember("type")->value.GetString(),"pe")==0);

    string fp_a = doc.FindMember("fastq")->value.FindMember("path")->value.FindMember("archive")->value.GetString(),
      fq1_a     = fp_a + "/" + doc.FindMember("fastq")->value.FindMember("r1")->value.FindMember("basename")->value.GetString(),
      fq2_a     = fp_a + "/" + doc.FindMember("fastq")->value.FindMember("r2")->value.FindMember("basename")->value.GetString(),
      fp_s      = doc.FindMember("fastq")->value.FindMember("path")->value.FindMember("scratch")->value.GetString(),
      fq1_s     = fp_s + "/" + doc.FindMember("fastq")->value.FindMember("r1")->value.FindMember("basename")->value.GetString(),
      fq2_s     = fp_s + "/" + doc.FindMember("fastq")->value.FindMember("r2")->value.FindMember("basename")->value.GetString();

    int fq1m = doc.FindMember("fastq")->value.FindMember("r1")->value.FindMember("modification")->value.GetInt(),
      fq2m   =   doc.FindMember("fastq")->value.FindMember("r2")->value.FindMember("modification")->value.GetInt();

    long long fq1s =   doc.FindMember("fastq")->value.FindMember("r1")->value.FindMember("size")->value.GetInt64(),
      fq2s         =   doc.FindMember("fastq")->value.FindMember("r2")->value.FindMember("size")->value.GetInt64();

    // cout << "we run\n\nfq1= " << fq1_a << " ("<<fq1m <<"/"<<fq1s<<")\nfq2= " << fq2_a << " ("<<fq2m<<"/"<<fq2s<<")\n\n";

    checks::check_rsync_checksum(fq1_a,fq1m,fq1s);
    checks::check_rsync_checksum(fq2_a,fq2m,fq2s);

    string bam = x["sample_internal_name"]+"."+x["experiment_id"]+"."+x["PU"]+".bam", bam_f = fp_s + "/" + bam;

    int bm = doc.FindMember("bam")->value.FindMember("modification")->value.GetInt();
    long long bs =   doc.FindMember("bam")->value.FindMember("size")->value.GetInt64();
    checks::check_rsync_checksum(bam_f,bm,bs);

    assert(bam_check(bam_f.data(),x));

    string md5f1 = adir + "/" + doc.FindMember("fastq")->value.FindMember("r1")->value.FindMember("basename")->value.GetString() + ".md5sum",
      md5f2 = adir + "/" + doc.FindMember("fastq")->value.FindMember("r2")->value.FindMember("basename")->value.GetString() + ".md5sum";
    assert(isregfile(md5f1.data()));
    assert(isregfile(md5f2.data()));

    // assert(isregfile( (fp_s+"/laneBarcode.html").data() ));

    string cmds;

    if(rgn==0 && isregfile( (fp_s+"/laneBarcode.html").data() )) {
        cmds=RSYNC_CLEAN+fp_s+"/laneBarcode.html "+adir;
        cout << "using\n"<< cmds<<"\n";
        assert(system(cmds.data())==0);
        Lazy("sudo /bin/chown dragen.analysts %s", (fp_a+"/laneBarcode.html").data() );
        Lazy("sudo /bin/chmod 440 %s", (fp_a+"/laneBarcode.html").data() );
    }

    cmds=RSYNC_CLEAN+fq1_s+" "+adir;
    cout << "using\n"<< cmds<<"\n";

    // assert(system(cmds.data())==0);
    if(system(cmds.data())==0) {
        cout << "YAY\n";
    }else{
        // cout << "check for " << fq1_a << ", " << fq1s << ", " << fq1m << "\n";
        if(rsync_checksum_return(fq1_a,fq1m,fq1s)) {
            cout << "we will let this slide for now\n";
        }else{
            cout << "we need to make this a specific error state!?!\n";
            exit(1);
        }
    }

    Lazy("sudo /bin/chown dragen.analysts %s",fq1_a.data());
    Lazy("sudo /bin/chmod 440 %s",fq1_a.data());
    Lazy("sudo /bin/chown dragen.analysts %s",(fq1_a+".md5sum").data());
    Lazy("sudo /bin/chmod 440 %s",(fq1_a+".md5sum").data());

    cmds=RSYNC_CLEAN+fq2_s+" "+adir;
    // cout << "using\n"<< cmds<<"\n";

    // assert(system(cmds.data())==0);
    if(system(cmds.data())==0) {
        cout << "YAY\n";
    }else{
        // cout << "check for " << fq2_a << ", " << fq2s << ", " << fq2m << "\n";
        if(rsync_checksum_return(fq2_a,fq2m,fq2s)) {
            cout << "we will let this slide for now\n";
        }else{
            cout << "we need to make this a specific error state!?!\n";
            exit(1);
        }
    }

    /////// did we forget chgrp?!?
    Lazy("sudo /bin/chown dragen.analysts %s",fq2_a.data());
    Lazy("sudo /bin/chmod 440 %s",fq2_a.data());
    Lazy("sudo /bin/chown dragen.analysts %s",(fq2_a+".md5sum").data());
    Lazy("sudo /bin/chmod 440 %s",(fq2_a+".md5sum").data());

}

    Lazy("sudo /bin/chown dragen.analysts %s",adir.data());
    Lazy("sudo /bin/chmod 550 %s",adir.data());

    for (unsigned int rgn=0;rgn<fc.size();++rgn) {

        rarp::NLIST & x=fc[rgn];

        char arsv2[8*1024]; sprintf(arsv2,"update Lane set rg_status = 'fastq_archived', data = json_remove(data,'$.fastq.path.scratch') where rg_status = 'protect_lock' "
        "and prepid = %s and lanenum = %s and fcid = '%s' ; select row_count() updated",
        x["prepid"].data(),x["lanenum"].data(),x["fcid"].data()
        );

        NLIST u = db::get_named_row("seqdb",arsv2);  assert(u["updated"]=="1"); 

    }

    cout << "done with " << fc[0]["prepid"]<<"\n";

    return;
}

void archive(int /* argc */, char ** /* argv */) {

    using namespace std;
    using namespace rarp;

    NLISTS fc;
    db::get_named_rows("seqdb",fc,"select concat(f.fcillumid,'.',l.lanenum) PU, sample_internal_name, exomekit capture_kit, sample_type, is_released, "
        "f.fcillumid,l.fcid, "
        " p.adapter_id, a.sequence adapterlet, " // " p.adapterlet, "
        " l.lanenum,f.machine,f.machinetype,f.chemver,rg_status, experiment_id, p_prepid pseudo_prepid, l.prepid, "
        "data, upper(sample_type) seq_type_upper, libkit, dateread1 "
        "from Lane l join Flowcell f on l.fcid=f.fcid "
        " join prepT p on l.prepid=p.prepid "
        " join adapter a on a.id=p.adapter_id "
        " where "
        " rg_status = 'fastq_mapped' "
        " order by prepid "
        "limit 1");

    if(fc.size()==0) {
        return;
        // exit(0);
    }

    assert(fc.size()==1);
    NLIST & x = fc[0];

    cout << "\n\n----------------\nwe will process " << fc[0]["data"] << "\n";
    // for(NLIST::iterator it = x.begin(); it != x.end(); it++ ) cout << "["<<it->first<<"] "<<it->second<<"\n";

    rapidjson::Document doc; 
    char buffer[16*1024]; 
    memset(buffer,0,sizeof(buffer));
    memcpy(buffer, x["data"].data(), x["data"].length());

    assert(x["data"]!="NULL");
    cout << "parsing " << x["data"] << "\n";
    if (doc.ParseInsitu(buffer).HasParseError()) {
        cout << "there's an issue with the json\n";
        exit(1);
    }

    // 'should' we recheck md5sum?!?
    
    { char arsv2[8*1024];
    sprintf(arsv2,"update Lane set rg_status = 'archive_lock' where rg_status = 'fastq_mapped' and prepid = %s and lanenum = %s and fcid = '%s' ; select row_count() updated",fc[0]["prepid"].data(),fc[0]["lanenum"].data(),fc[0]["fcid"].data());
    NLIST u = db::get_named_row("seqdb",arsv2); 
    assert(u["updated"]=="1"); }

    string adir = doc.FindMember("fastq")->value.FindMember("path")->value.HasMember("archive") ? 
      doc.FindMember("fastq")->value.FindMember("path")->value.FindMember("archive")->value.GetString() :
      seq::se_archive_dir(x);

    // cout << "arcvhive dir should be " << adir << "\n";

    // this is silly and way too repetitive
    assert(doc.FindMember("fastq")->value.FindMember("path")->value.HasMember("scratch"));
    assert(doc.FindMember("fastq")->value.HasMember("type"));
    assert(doc.FindMember("fastq")->value.FindMember("r1")->value.HasMember("size"));
    assert(doc.FindMember("fastq")->value.FindMember("r2")->value.HasMember("size"));

    cout << "checking " << doc.FindMember("fastq")->value.FindMember("type")->value.GetString() << "\n";

    assert(strcmp(doc.FindMember("fastq")->value.FindMember("type")->value.GetString(),"pe")==0);
    string fp_s = doc.FindMember("fastq")->value.FindMember("path")->value.FindMember("scratch")->value.GetString(),
      fq1_s     = fp_s + "/" + doc.FindMember("fastq")->value.FindMember("r1")->value.FindMember("basename")->value.GetString(),
      fq2_s     = fp_s + "/" + doc.FindMember("fastq")->value.FindMember("r2")->value.FindMember("basename")->value.GetString();

    int fq1m = doc.FindMember("fastq")->value.FindMember("r1")->value.FindMember("modification")->value.GetInt(),
      fq2m   =   doc.FindMember("fastq")->value.FindMember("r2")->value.FindMember("modification")->value.GetInt();

    long long fq1s =   doc.FindMember("fastq")->value.FindMember("r1")->value.FindMember("size")->value.GetInt64(),
      fq2s         =   doc.FindMember("fastq")->value.FindMember("r2")->value.FindMember("size")->value.GetInt64();

    // cout << "we run\n\nfq1= " << fq1_s << " ("<<fq1s <<"/"<<fq1s<<")\nfq2= " << fq2_s << " ("<<fq2m<<"/"<<fq2s<<")\n\n";

    checks::check_rsync_checksum(fq1_s,fq1m,fq1s);
    checks::check_rsync_checksum(fq2_s,fq2m,fq2s);

    string ckpnt = fp_s + "/" + x["sample_internal_name"]+"."+x["experiment_id"] + "." + x["PU"] + ".done", 
        bam = x["sample_internal_name"]+"."+x["experiment_id"]+"."+x["PU"]+".bam",
        bam_f = fp_s + "/" + bam;

    // cout << "we should have bam file " << bam_f << " and checkpoint " << ckpnt << "\n";

    /// feeling so lazy so update to error briefly and then update to checked after
    assert(isregfile(bam_f.data()));
    assert(isregfile(ckpnt.data()));
    int bm = doc.FindMember("bam")->value.FindMember("modification")->value.GetInt();
    long long bs =   doc.FindMember("bam")->value.FindMember("size")->value.GetInt64();
    checks::check_rsync_checksum(bam_f,bm,bs);

    assert(bam_check(bam_f.data(),x));

    string cmds;
    if(isregfile( (fp_s+"/laneBarcode.html").data() )) {
        cout << "we have a barcode file\n";
        cmds=RSYNC;
        cmds+=fp_s+"/laneBarcode.html "+adir;
        cout << "using\n"<< cmds<<"\n";
        assert(system(cmds.data())==0);
    }else{
        cout << "there's no barcode file\n";
    }
    cmds=RSYNC+fq1_s+" "+adir;
    cout << "using\n"<< cmds<<"\n";
    assert(system(cmds.data())==0);
    cmds=RSYNC+fq2_s+" "+adir;
    cout << "using\n"<< cmds<<"\n";
    assert(system(cmds.data())==0);
    // cout << "do md5sum and lock it down\n";

    string fq_a = adir + "/" + doc.FindMember("fastq")->value.FindMember("r1")->value.FindMember("basename")->value.GetString();
    assert(chmod(fq_a.data(),S_IRUSR|S_IRGRP)==0);
    cmds="md5sum "+fq_a+" | tee "+fq_a+".md5sum";
    // cout << "using\n"<< cmds<<"\n";
    assert(strchr(fq_a.data(),' ')==0);
    if(!isregfile((fq_a+".md5sum").data())) assert(system(cmds.data())==0);

    fq_a = adir + "/" + doc.FindMember("fastq")->value.FindMember("r2")->value.FindMember("basename")->value.GetString();
    assert(chmod(fq_a.data(),S_IRUSR|S_IRGRP)==0);
    cmds="md5sum "+fq_a+" | tee "+fq_a+".md5sum";
    // cout << "using\n"<< cmds<<"\n";
    assert(strchr(fq_a.data(),' ')==0);
    if(!isregfile((fq_a+".md5sum").data())) assert(system(cmds.data())==0);
        
    { char arsv2[8*1024];
    sprintf(arsv2,"update Lane set rg_status = 'fastq_archiving' where rg_status = 'archive_lock' and prepid = %s and lanenum = %s and fcid = '%s' ; select row_count() updated",fc[0]["prepid"].data(),fc[0]["lanenum"].data(),fc[0]["fcid"].data());
    // cout << arsv2 << "\n";
    NLIST u = db::get_named_row("seqdb",arsv2); 
    assert(u["updated"]=="1"); }
    ///// we do any heavier stuff apart from mapping here to keep dragen free?!?

    return;

}

bool bam_check(char const * bamf, rarp::NLIST & x) {

    using std::cout;

    cout << "checking bam " << bamf << "\n";
    if(!isregfile(bamf)) return false;

    char cmd[2048];
    int blarp=250;
    sprintf(cmd,"samtools view -h %s hs37d5 2>&1 | head -n%d",bamf,blarp);

    // cout << "ARGH\n"<<cmd<<"\n";

    Popen ps(cmd,16*1024,"r");
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

bool check_final_bam( 
    rarp::NLIST /* allow operator[] : const */ & dsm,
    std::string const & rg_from_lims,
    std::string const & final_bam,
    bool igm,
    std::string & email_body,
    bool merge
) {

    using namespace std;

    stringstream issues;

    if(!isregfile(final_bam.data())) { 
        issues << "ERROR: THERE IS NO FINAL BAM " << final_bam
          << " (" << dsm["sample_type"] << ")'\n";
        update_status(dsm,(merge?MERGE_ERROR:FINAL_BAM_MISSING),email_body,issues.str());
        return false;
    }

    // if(atoi((get_single_line_output_as_string("samtools view -c %s MT",final_bam.data()).data()))<10){
    if( atoi(dsm["min_prepid"].data())>50000 && dsm["sample_type"]!="Custom_Capture" 
      && atoi((get_single_line_output_as_string("samtools view -c %s MT",final_bam.data()).data()))<10
      && atoi((get_single_line_output_as_string("samtools view -c %s 21",final_bam.data()).data()))<10000
    ){
        issues << "ERROR: (" << dsm["pseudo_prepid"] << ") cannot find many MT reads for final bam " << final_bam << "\n";

        // update_status(dsm,FEW_MT,email_body,issues.str());
        // return false;

        update_status(dsm,FEW_MT,email_body,issues.str());
        return false;
    }

    string h_rg, b_rg, pg, RG_from_bam;
    bool okay = false, error = false, warn = false;

    char rg_counts[1024];

{
    char tmp[10240]; sprintf(tmp,"/nfs/seqscratch09/dsth/vcf_thing/bam_thing %s",final_bam.data()); Popen ns10(tmp,16*1024,"r");
    char * z20;
    while( *( z20=ns10.getline() ) != '\0') { 
        //// don't use size of with string literals!?!?
        // cout << "HAVE: "<<z20<<"\n";
        if(memcmp(z20,"ERROR",strlen( "ERROR") )==0) error=true;
        if(memcmp(z20,"OKAY",strlen( "OKAY") )==0) okay=true;
        if(memcmp(z20,"DODGY_OLD_BAM",strlen( "DODGY_OLD_BAM" ))==0) warn=true;
        if(::isdigit(z20[0])) strcpy(rg_counts,z20);
        if(memcmp(z20,"SampleRg=",strlen( "SampleRg=" ))==0) b_rg=z20+9; 
        else if(memcmp(z20,"@RG,",strlen( "@RG," ) )==0) h_rg=z20+4;
        else if(memcmp(z20,"@PG,",strlen( "@PG," ) )==0) pg=z20+4;
        else std::cout << "bam: '" << z20 << "'\n";
    }
}

    if(warn) cout << "WARN: THE BAM LOOKS LIKE IT'S PRETTY OLD?!?\n";

    if(!okay&&!error) {
        issues << "ERROR: (" << final_bam << ") THERE'S A PROBLEM CHECKING THE BAM\n";
        // update_status(dsm["pseudo_prepid"],PROB,email_body,issues.str());
        update_status(dsm,(merge?MERGE_ERROR:BAM_CHECKING_PROB),email_body,issues.str());
        return false;
    } else if(error) {
        issues << "ERROR: (" << final_bam << ") THE BAM HAS ISSUES\n";
        // update_status(dsm["pseudo_prepid"],PROB,email_body,issues.str());
        update_status(dsm,(merge?MERGE_ERROR:BAM_IS_A_MESS),email_body,issues.str());
        return false;
    }else{
        std::cout << "HRG="<<h_rg<<"\nBRG="<<b_rg<<"\nPG="<<pg<<"\n";

        std::set<string> a_header_rg_set, b_body_rg_set;
        for (unsigned r =0 ; r<h_rg.length() ; ++r ) {
            if(h_rg.substr(r,3)=="PU:" ) {
                int rr=r;
                while(h_rg[++rr]!=',') {}
                std::cout << "have one @ " << r << " - " << rr << "\n";
                a_header_rg_set.insert(h_rg.substr(r+3,rr-r-3));
            }
            if(h_rg.substr(r,3)=="SM:" ) {
                int rr=r;
                // while(h_rg[++rr]!=',') {}
                while(h_rg[++rr]!=',' && h_rg[rr]!=';') {}
                std::cout << "bam_sample_name = " << h_rg.substr(r+3,rr-r-3) << "\n";

                if(dsm["sample_name"]!=h_rg.substr(r+3,rr-r-3)) {

                    issues << "ERROR: (" << dsm["pseudo_prepid"] << ") BAM SAMPLE NAME MISMATCH " << dsm["sample_name"] << " VS., " << h_rg.substr(r+3,rr-r-3) << "\n";

                    // update_status(dsm["pseudo_prepid"],PROB,email_body,issues.str());
                    // char lazy[1024];
                    rarp::NLISTS tmp;

                    db::get_named_rows("seqdb",tmp,"select s.sample_internal_name,s.sample_external_name from SampleT s join Experiment e on s.sample_id=e.sample_id "
                      "where e.id = %s",dsm["experiment_id"].data());
                      // get_named_rows(tmp,"select sample_internal_name,sample_external_name from SampleT where sample_id = %s",prept["sample_id"].data());
                    assert(tmp.size()==1);

                    cout << "this has is_external " << dsm["is_external"] << " sample_internal_name= " << tmp[0]["sample_internal_name"] << " and sample_external_name = " << tmp[0]["sample_external_name"] << "\n";

                    string clean_h_rg=h_rg.substr(r+3,rr-r-3);
                    { char seriously[1024], *s_p=seriously;
                    memset(seriously,0,sizeof seriously);
                    for(unsigned int i=0;i<clean_h_rg.length(); ++i) if(clean_h_rg[i]!='-' && clean_h_rg[i]!='_' && clean_h_rg[i]!='.') *s_p++=clean_h_rg[i];
                    // for(int i=0;i<dsm["sample_name"].length(); ++i) if(dsm["sample_name"][i]!='-') *s_p++=dsm["sample_name"][i];
                    cout <<"CHECKING DETOXED RG NAME '" << seriously << "'\n";
                    clean_h_rg=seriously;
                    }

                    /////// this has a ticket registered (newer external) and sample_external_name matches SM:(samplename)
                    if(atoi(dsm["is_external"].data())>1 && dsm["sample_name"]==clean_h_rg) {
                        cout << "let's this slide... (detox)?!?\n";
                        // assert(0);
                    }else if(atoi(dsm["is_external"].data())>1 && tmp[0]["sample_external_name"]==h_rg.substr(r+3,rr-r-3)) {
                        cout << "let's this slide... (sample_external_name)?!?\n";
                        // assert(0);
                    }else{

                        issues << "ERROR: (" << dsm["pseudo_prepid"] << ") BAM SAMPLE VS sample_external_name MISMATCH " << tmp[0]["sample_external_name"] << " VS., " << h_rg.substr(r+3,rr-r-3) << "\n";
                        rarp::NLIST t = db::get_named_row("seqdb","update prepT set status = 'Sample_Name_Mismatch \"%s\" vs. \"%s\"' where experiment_id = %s ; select row_count() updated",
                          tmp[0]["sample_external_name"].data(),h_rg.substr(r+3,rr-r-3).data(),dsm["experiment_id"].data()
                        );
                        cout << "we modified " << t["updated"]<< "\n";
                        // sprintf(silly,"update prepT set status = 'Sample_Name_Mismatch \"%s\" vs. \"%s\"' where experiment_id = %s",
                        update_status(dsm,(merge?MERGE_ERROR:BAM_LIMS_SAMPLE_NAME_MISMATCH),email_body,issues.str());
                        return false;
                    }
                }

            }
        }

        if(!warn) { // cout << "not old bam\n";
            for (unsigned r =0, rr =0 ; r<b_rg.length() ; ++r ) {
                if(b_rg[r]==',') rr=0;
                else if(b_rg[r]==':') {
                    ++rr;
                    if(rr==2) { 
                        unsigned rrr = r;
                        while(rrr<b_rg.length() && b_rg[++rrr]!=',') {}
                        b_body_rg_set.insert(b_rg.substr(r+1,rrr-r-1));
                    }
                }
            }
            string sa_header_rg_set_as_sorted_string, sb_body_rg_set_as_sorted_string;
            for (std::set<string>::iterator it = a_header_rg_set.begin(); it!=a_header_rg_set.end(); it++)
            if(sa_header_rg_set_as_sorted_string.empty()) sa_header_rg_set_as_sorted_string=*it;
            else sa_header_rg_set_as_sorted_string += ',' + *it;
            for (std::set<string>::iterator it = b_body_rg_set.begin(); it!=b_body_rg_set.end(); it++)
            if(sb_body_rg_set_as_sorted_string.empty()) sb_body_rg_set_as_sorted_string=*it;
            else sb_body_rg_set_as_sorted_string += ',' + *it;
            for (unsigned g =0 ; g < sb_body_rg_set_as_sorted_string.length(); ++g ) if(sb_body_rg_set_as_sorted_string[g]==':') sb_body_rg_set_as_sorted_string[g]='.';
            RG_from_bam = sb_body_rg_set_as_sorted_string;

            if(sa_header_rg_set_as_sorted_string==sb_body_rg_set_as_sorted_string) {
                cout << "HEADER/BODY RG MATCH= " << RG_from_bam << "\n";
            } else {

                // cout << "line="<<__LINE__<<"\n";
                std::string rn=get_single_line_output_as_string("samtools view %s | head -n1 | cut -f1",final_bam.data()),
                    fc, bc, mtype_from_fc, mtype_from_pfx;
                rn='@'+rn;
                cout << final_bam << "\n";
                cout << "igm="<<igm<<" " << rn <<"\n";
                bool is18=fastq::rncheck(rn.data(),fc, bc, mtype_from_fc, mtype_from_pfx);
                cout << "1.8="<<is18<<"\n";

                cout << "FCID= '"<<fc<<"'\nbc= "<<bc<<"\nmtype_from_fc= "<<mtype_from_fc<<"\nmtype_from_pfx= "<<mtype_from_pfx<<"\n";

                // if(igm||is18) {
                if(igm) {

                    rarp::LIST rgs; tokenise(rgs,rg_counts,':');
                    assert(rgs.size()==3);
                    int header_rg_count = atoi(rgs[0].data()), body_rg_sampling_count = atoi(rgs[2].data());

                    if(header_rg_count>=10 && header_rg_count==body_rg_sampling_count) {
                        cout << "\n\t> This is internal and PU is likley simply truncated there are so many RGs - need to print whole thing and/or just use FC?!?\n\n";
                        return true;
                    }else if(atoi(dsm["max_prepid"].data())<50000) {
                        cout << "\n\t> This is a mess...\n";
                        return true;
                    }else{

                        issues << "ERROR: (" << dsm["pseudo_prepid"] << ") BAM HEADER/BODY RG PU/READNAME MISMATCH : header= " 
                        << sa_header_rg_set_as_sorted_string << ", body= " << sb_body_rg_set_as_sorted_string << "\n"; 
                        // update_status(dsm["pseudo_prepid"],PROB,email_body,issues.str());
                        update_status(dsm,RG_PU_VERSUS_READNAME_ERROR,email_body,issues.str());
                        return false;

                    }

                }else if(is18) {

                    cout << "have rg_counts " << rg_counts << "\n";
                    rarp::LIST rgs; tokenise(rgs,rg_counts,':');
                    assert(rgs.size()==3);
                    int header_rg_count = atoi(rgs[0].data()), body_rg_sampling_count = atoi(rgs[2].data());

                    if(header_rg_count==body_rg_sampling_count) {
                        cout << "\n\t> This is external and while PU are misnamed distinct RG tags and readname/lane counts match\n\n";
                        return true;
                    }else if(header_rg_count>=5 && (double)header_rg_count/body_rg_sampling_count > 0.8) {
                        cout << "\n\t> This is external and while PU are misnamed distinct RG tags and readname/lane counts mismatch it's likely this is due to insufficient sampling with many RGs?!? - WE SHOULD HOWEVER ATTEMPT TO NAME MATCH WITH SET - I.E. INTERSECT VS. UNION COUNT SHOULD BE APPROPRIATE!?!\n\n";
                        return true;
                    }else if(!igm) { cout << "i'm sick of dealing with external so whatever...\n";
                        return true;
                    }else{
                        issues << "ERROR: (" << dsm["pseudo_prepid"] << ") BAM HEADER/BODY RG COUNT MISMATCH : " << rg_from_lims << " : " << RG_from_bam << "\n"; 
                        update_status(dsm,RG_COUNT_ERROR,email_body,issues.str());
                        return false;
                    }

                }else{
                    cout << "\n\t> This is external and old or strange so whatever!?!\n\n";
                    return true;
                }

            }
        }

        if( ( igm && rg_from_lims==RG_from_bam ) 
              || !igm 
              || atoi(dsm["max_prepid"].data())<50000 // have no FCBC for old samples - could/should pull from old bam?!?
              // || atoi(dsm["pseudo_prepid"].data())<30000 // have no FCBC for old samples - could/should pull from old bam?!?
        ) {

            cout << "have rg_counts " << rg_counts << "\n";
            rarp::LIST rgs; tokenise(rgs,rg_counts,':');
            assert(rgs.size()==3);

            if(atoi(rgs[0].data())==atoi(rgs[2].data())) {
                cout <<"WILL LET THIS " << (igm?"IGM":"NON-IGM") << " MERGED SAMPLE THROUGH\n";
                cout << "rg2= " << RG_from_bam << "\n";
                return true;
            }else{
                cout << "line="<<__LINE__<<"\n";

                std::string rn=get_single_line_output_as_string("samtools view %s | head -n1 | cut -f1",final_bam.data()),
                    fc, bc, mtype_from_fc, mtype_from_pfx;
                rn='@'+rn;
                cout << final_bam << "\n";
                cout << "igm="<<igm<<" " << rn <<"\n";
                bool is18=fastq::rncheck(rn.data(),fc, bc, mtype_from_fc, mtype_from_pfx);
                cout << "1.8="<<is18<<"\n";

                cout << "FCID= '"<<fc<<"'\nbc= "<<bc<<"\nmtype_from_fc= "<<mtype_from_fc<<"\nmtype_from_pfx= "<<mtype_from_pfx<<"\n";

                if(!igm) { cout << "i'm sick of dealing with external so whatever...\n";
                    return true;
                }else if(igm&&is18) {

                    issues << "ERROR: (" << dsm["pseudo_prepid"] << ") BAM HEADER/BODY RG COUNT MISMATCH : " << rg_from_lims << " : " << RG_from_bam << "\n"; 
                    // update_status(dsm["pseudo_prepid"],PROB,email_body,issues.str());
                    update_status(dsm,RG_COUNT_ERROR,email_body,issues.str());
                    // update_status(dsm["pseudo_prepid"],(merge?MERGE_ERROR:JUNK),email_body,issues.str());
                    return false;

                }else{
                    cout << "\n\t> This is external and old or strange so whatever!?!\n\n";
                    return true;
                }
            }

        } else {
            issues << "ERROR: (" << dsm["pseudo_prepid"] << ") SEQDB/BAM RG MISMATCH : " << rg_from_lims << " : " << RG_from_bam << "\n"; 
            update_status(dsm,RG_LIMS_READNAME_MISMATCH_ERROR,email_body,issues.str());
            // update_status(dsm["pseudo_prepid"],PROB,email_body,issues.str());
            // update_status(dsm["pseudo_prepid"],(merge?MERGE_ERROR:JUNK),email_body,issues.str());
            return false;
        }
    } 

    return false;
}

void release_merged_rgs(
    rarp::NLIST /* allow operator[] : const */ & dsm,
    std::string const & rg_from_lims,
    FILE * const log,
    std::string const & final_bam,
    std::string const & final_checkpoint,
    bool igm,
    std::string & email_body,
    std::string const & /* merge_intermediate */,
    std::string const & sh,
    std::stringstream & issues,
    std::string const & dir
) {

    using namespace std;
        
    // some checks got turned off cos of the way genomes_as_fake_exomes were done 
    // when out of space it seems that return may have been okay and things bulldozed through?!?

    // nasty!?!?
    if(get_single_line_output_as_string("grep -c -P 'Run| picard' %s.err",sh.substr(0,sh.length()-3).data())!="6") {
        issues << "ERROR: (" << dsm["pseudo_prepid"] << ") there appears to be an issue with the error file\n";
        update_status(dsm,MERGE_ERROR,email_body,issues.str());
        return;
    }

    if(get_single_line_output_as_string("grep -c 'Runtime' %s.err",sh.substr(0,sh.length()-3).data())!="2") {
        issues << "ERROR: (" << dsm["pseudo_prepid"] << ") there appears to be an issue with the error file\n";
        update_status(dsm,MERGE_ERROR,email_body,issues.str());
        return;
    }

    if(get_single_line_output_as_string("grep -c 'picard' %s.err",sh.substr(0,sh.length()-3).data())!="4") {
        issues << "ERROR: (" << dsm["pseudo_prepid"] << ") there appears to be an issue with the error file\n";
        update_status(dsm,MERGE_ERROR,email_body,issues.str());
        return;
    }

    if(get_single_line_output_as_string("grep -c 'done. Elapsed' %s.err",sh.substr(0,sh.length()-3).data())!="2") {
        issues << "ERROR: (" << dsm["pseudo_prepid"] << ") there appears to be an issue with the error file\n";
        update_status(dsm,MERGE_ERROR,email_body,issues.str());
        return;
    }

    if(dsm["sample_type"]!="Genome_As_Fake_Exome") {

        cout << "we check DT string as this is " << dsm["sample_type"] << "\n";

        if(get_single_line_output_as_string("grep -c 'DT tag value' %s.err",sh.substr(0,sh.length()-3).data())!="0") {
            issues << "ERROR: (" << dsm["pseudo_prepid"] << ") there appears to be an issue with the DT tags\n";
            update_status(dsm,DT_TAGS,email_body,issues.str());
            return;
        }else if(get_single_line_output_as_string("grep -c 'Error parsing SAM header' %s.err",sh.substr(0,sh.length()-3).data())!="0") {
        // }else if(get_single_line_output_as_string("grep -c 'Error parsing SAM header' %s.err",sh.substr(0,sh.length()-3).data())!="0") {
            issues << "ERROR: (" << dsm["pseudo_prepid"] << ") there appears to be an issue with the DT tags\n";
            update_status(dsm,DT_TAGS,email_body,issues.str());
            return;
        }

    }else{
        cout << "we don't check DT string as this is " << dsm["sample_type"] << "\n";
    }

    // really should check metrics file for library complexity 
    // we get lots of alignstats jobs dying. almost certainly mem - too lazy to check. means re-release gets old value
    // wipe origianl to avoid last result comign back through and force
    char bored[1024];
    // cout << "using " << sh << "\n";
    { sprintf(bored,POSTMERGE "/%s.%s.txt",dsm["experiment_id"].data(),dsm["sample_name"].data());
    // cout << "checking " << bored << "\n";
    if(!isregfile(bored)) {
        issues << "ERROR: (" << dsm["pseudo_prepid"] << ") there appears to be an issue with the metrics file\n";
        update_status(dsm,MERGE_METRICS_ERROR,email_body,issues.str());
        return;
    } }
    // cout << "we need to pull in the alignment stats here!?!?!!?!?\n";
    // cout << "if insufficient we store metrics BUT unrelease sample!?!?!?!\n";

    if( check_final_bam( dsm, rg_from_lims, final_bam, igm, email_body, true /*, prepidi, prept */) ) {

        // cout << "this seems fine\n";
        // this probably needs to be picu aware!?!??!
        // this is truelly awful!?!?
        string picard_dups=dir+dsm["sample_name"]+"."+dsm["experiment_id"]+".metrics_duplication.txt",
          merge_intermediate = dir+dsm["sample_name"]+"."+dsm["experiment_id"]+".merge.bam";

        if(isregfile(merge_intermediate.data())) {
            // we must apply a samtools quickcheck here!?!
            // cout << "we remove merge intermediate " << merge_intermediate << "\n";
            unlink(merge_intermediate.data());
        } // else cout << "there is no merge intermediate " << merge_intermediate << "\n";

        string core_q = db::get_core_query(dsm["experiment_id"]);

        assert(isregfile(bored));

        // cout << "core_q= " << core_q << "\n";
        core::Core core_info(core_q.data());

        float capmean=-0.0,capmedian=-0.0;

        // if touching core tables MUST use rowcount!!?!?!?!
        // don't undo prepT.is_released for billing reasons?!?
        if(dsm["sample_type"]=="Genome_As_Fake_Exome") {

            capmean=atof(get_single_line_output_as_string("grep \"WgsCoverageMean\" %s | awk '{print $2}'",bored).data()),
              capmedian=atof(get_single_line_output_as_string("grep \"WgsCoverageMedian\" %s | awk '{print $2}'",bored).data());
            // cout << "has coverage = " << capmean << "\n";

            if( !core_info.hits_coverage("Genome",capmean) ) { // capmean < ( core::get_min_wgs() 

                // if(capmean < ( 0.97 * opts::wgs_min ) ) {
                // do_void_thing("update dragen_sample_metadata set is_merged = 200001 where pseudo_prepid = %s",dsm["experiment_id"].data());
                do_void_thing("delete from dragen_sample_metadata where pseudo_prepid = %s",dsm["experiment_id"].data());
                do_void_thing("delete from dragen_pipeline_step where pseudo_prepid = %s",dsm["experiment_id"].data());
                do_void_thing("update Experiment set is_released = 'not_released', merge_metrics_capturemean = %0.2f where id = %s",capmean,dsm["experiment_id"].data());
                // do_void_thing("update Experiment set is_released = 'release_rejected' where id = %s",dsm["experiment_id"].data());
                do_void_thing("update prepT set status = 'Release_Rejected_UnderCov(%0.2f)', status_time = unix_timestamp() where experiment_id = %s",capmean,dsm["experiment_id"].data());
                return;
            }else{
                do_void_thing("update Experiment set merge_metrics_capturemean = %0.2f where id = %s",capmean,dsm["experiment_id"].data());
            }
        }else{
            capmean=atof(get_single_line_output_as_string("grep \"CapCoverageMean\" %s | awk '{print $2}'",bored).data()),
              capmedian=atof(get_single_line_output_as_string("grep \"CapCoverageMedian\" %s | awk '{print $2}'",bored).data());
            // cout << "has coverage = " << capmean << "\n";

            if( !core_info.hits_coverage("Exome",capmean) ) { // capmean < ( core::get_min_wes() ) 

                // if(capmean < ( 0.97 * opts::wes_min ) ) {
                do_void_thing("delete from dragen_sample_metadata where pseudo_prepid = %s",dsm["experiment_id"].data());
                do_void_thing("delete from dragen_pipeline_step where pseudo_prepid = %s",dsm["experiment_id"].data());
                // do_void_thing("update Experiment set is_released = 'release_rejected' where id = %s",dsm["experiment_id"].data());
                do_void_thing("update Experiment set is_released = 'not_released', merge_metrics_capturemean = %0.2f where id = %s",capmean,dsm["experiment_id"].data());
                // do_void_thing("update Experiment set is_released = 'not_released' where id = %s",dsm["experiment_id"].data());
                do_void_thing("update prepT set status = 'Release_Rejected_UnderCov(%0.2f)', status_time = unix_timestamp() where experiment_id = %s",capmean,dsm["experiment_id"].data());
                return;
            }else{
                do_void_thing("update Experiment set merge_metrics_capturemean = %0.2f where id = %s",capmean,dsm["experiment_id"].data());
            }
        }

        // simple way to avoid this is to replace into 
        rarp::NLIST ars2 = db::get_named_row("seqdb",
            "update dragen_pipeline_step set finish_time = CURRENT_TIMESTAMP(), step_status = 'completed' where pseudo_prepid = %s and pipeline_step_id = 1; select row_count() as affected",
            dsm["experiment_id"].data()
        );

        // for (rarp::NLIST::iterator it = ars2.begin(); it!=ars2.end(); it++ ) cout << "ars2["<< it->first <<"] "<< it->second << "\n"; 

        if(ars2["affected"]=="0") {
            /////// MUST WIPE DPS OR WILL CAUSE MORE HASSLES LATER!?!
            do_void_thing( REPLACE_INTO_DPS, dsm["experiment_id"].data(),"completed" );
        }else {
            cout << "running with " << dsm["experiment_id"].data() << "\n";
            assert(ars2["affected"]=="1");
        }
    
        do_void_thing("update dragen_sample_metadata set is_merged = 1 where pseudo_prepid = %s",dsm["experiment_id"].data());
        do_void_thing("update prepT set status = 'Released_to_Pipeline', status_time = unix_timestamp() where experiment_id = %s",dsm["experiment_id"].data());

        time_t t = time(0); struct tm * tm_s = localtime(&t); char blah[1024]; strftime(blah,1024,"%c",tm_s); fputs(blah,log); fputc('\n',log);

        // fputs(bored2,log); fputc('\n',log); fputs(bored1,log); fputc('\n',log);
        fputs(final_checkpoint.data(),log);
        fputc('\n',log);

    }

}

void merge_multiple(
    // rarp::LIST const & pieces,
    // LIST const & SEs,
    std::string const & final_bam,
    std::string const & final_checkpoint,
    FILE * const log,
    std::string & merge_cmd,
    std::string const & dir,
    std::string const & gg,
    std::string const & merge_intermediate,
    rarp::NLIST & dsm
) {

    using namespace std;

    merge_cmd = "/nfs/goldstein/software/jdk1.8.0_05/bin/java -Xms12g -Xmx16g -jar /nfs/goldstein/software/picard-tools-2.2.1/picard.jar MergeSamFiles " 
      // + merge_cmd + " O="+ merge_intermediate +" TMP_DIR=/nfs/seqscratch_ssd/MERGE_TEMP MAX_RECORDS_IN_RAM=10000000 VALIDATION_STRINGENCY=STRICT";
      + merge_cmd + " O="+ merge_intermediate +" TMP_DIR=/nfs/seqscratch_ssd/MERGE_TEMP MAX_RECORDS_IN_RAM=10000000 VALIDATION_STRINGENCY=LENIENT";

    // cout << "change this to sambamba and samtools and simply parse out flagstat in entryppoint checks\n";

    string markdups_cmd = "/nfs/goldstein/software/jdk1.8.0_05/bin/java -Xms12g -Xmx16g -jar /nfs/goldstein/software/picard-tools-2.2.1/picard.jar MarkDuplicates "
      " I="+ dir 
      +dsm["sample_name"] +"."
      +dsm["experiment_id"]
      +".merge.bam O="
      // have the DT issue - consider making stringency high - 'should' get proper marking via name sort prior to marking?!?
      // consider using VALIDATION_STRINGENCY=STRICT
      // + final_bam + " TMP_DIR=/nfs/seqscratch_ssd/MERGE_TEMP MAX_RECORDS_IN_RAM=10000000 VALIDATION_STRINGENCY=STRICT METRICS_FILE="+dir
      + final_bam + " TMP_DIR=/nfs/seqscratch_ssd/MERGE_TEMP MAX_RECORDS_IN_RAM=10000000 VALIDATION_STRINGENCY=LENIENT METRICS_FILE="+dir
      +dsm["sample_name"]+"."
      +dsm["experiment_id"]
      +".metrics_duplication.txt CREATE_INDEX=true REMOVE_DUPLICATES=false REMOVE_SEQUENCING_DUPLICATES=false";

    string go= (string(SCRIPT_DIR) +dsm["sample_name"]+"." +dsm["experiment_id"]+".out");
    string ge= (string(SCRIPT_DIR) +dsm["sample_name"]+"." +dsm["experiment_id"]+".err");

    // started so don't stomp...?!?
    if(isregfile(gg.data())) {
        // cout << "There's something really nasty going on " << gg << " but we aren't being so uptight anymore but will clean up\n";
        unlink(gg.data());
        // exit(1);
    }

    if(isregfile(go.data())) unlink(go.data());
    if(isregfile(ge.data())) unlink(ge.data());
    if(isregfile(final_checkpoint.data())) unlink(final_checkpoint.data());

    // splitting_release!?!
    char bored_alignstats[6*2048];
    { string bed; //horrible!?!
    if(dsm["sample_type"]=="Genome_As_Fake_Exome") bed="/nfs/seqscratch_ssd/PIPELINE_DATA/ccds_regions.bed"; // cout << "This is genomic so just use CCDS\n";
    else{
        rarp::NLIST quero_chorar = db::get_named_row("seqdb","select count(1) count,region_file_lsrc bed from captureKit where prepT_name = '%s' and chr = 'all'",dsm["capture_kit"].data());
        assert(quero_chorar["count"]=="1");
        bed=quero_chorar["bed"];
    }
    // nasty, nasty, nasty. need to make sure repeat invocations were newer one dies doesnt pick up last result. all needs factoring...!?!
    { char bored[1024];
    sprintf(bored,POSTMERGE"/%s.%s.txt",dsm["experiment_id"].data(),dsm["sample_name"].data());
    if(isregfile(bored)) unlink(bored); }
    //////////// add in CCDS check for exomes too!?! iff there is CCDS min value too have fudge factor to release when a sample is within ~5x of min requirement
    snprintf(
      bored_alignstats,sizeof(bored_alignstats),ALIGNSTATS " -q 10 -i %s -t %s "
      // bored_alignstats,sizeof(bored_alignstats),"/home/dh2880/tools/alignstats/alignstats -i %s -t %s "
      " -o " POSTMERGE "/%s.%s.txt >/dev/null 2>&1 ",
      final_bam.data(),bed.data(),dsm["experiment_id"].data(),dsm["sample_name"].data()
    ); }

    char FER2[16222];
    sprintf(FER2,"#!/bin/bash\n#$ -S /bin/bash\n#$ -cwd\n" 
        "#$ -V\n#$ -l mem_free=15G\n#$ -pe threaded 8\n" // 3\n"
        // "#$ -V\n#$ -l mem_free=15G\n#$ -p 99\n#$ -pe threaded 8\n" // 3\n"
        "#$ -M dh2880@cumc.columbia.edu\n#$ -o %s\n#$ -e %s\n"
        "#$ -N merge_mark_%s.%s\n"
        /* "set -e\n"
        "# this really isn't how we should be doing this!?!\n" */
        "%s\n%s\n%s\n/nfs/goldstein/software/samtools/samtools-1.5/samtools index %s\ntouch %s\n",
        go.data(), ge.data(),
        dsm["sample_name"].data(),dsm["experiment_id"].data(),
        merge_cmd.data(),markdups_cmd.data(),bored_alignstats,
        final_bam.data(),final_checkpoint.data()
    );


    FILE * qsub = fopen(gg.data(),"w");
    if(!qsub) {
        cout << "can't write to script file\n";
        exit(1);
    }
    // printf("will put='%s'\ninto='%s'\n",FER2,gg.data());
    int x=fputs(FER2,qsub);
    if(x==EOF) {
        printf("unable to write to script file (%d/%d)\n",x,(int)strlen(FER2));
        fclose(qsub);
        exit(1);
    }
    fclose(qsub);
    chmod(gg.data(),S_IRUSR|S_IWUSR|S_IXUSR);


    char FER3[1024], /* bored1[1024], */ bored2[1024];
    sprintf(FER3,"qsub %s",gg.data()); // { Popen ns1(FER3,16*1024,"r"); char * z2 = ns1.getline(); 

    do_void_thing ( "update dragen_sample_metadata set is_merged = -1 where pseudo_prepid = %s", dsm["experiment_id"].data() );
    do_void_thing( "update prepT set status = 'Bam_Merging', status_time = unix_timestamp() where p_prepid = %s", dsm["experiment_id"].data() );

    if(system(FER3)!=0) { 
        cout << "ouch: " << FER3 << "\n\n\n";
        exit(1);                
    }

    cout << "Merging Expt=" << dsm["experiment_id"] << " (" << dsm["sample_name"] << ")\n";

    do_void_thing( REPLACE_INTO_DPS, dsm["experiment_id"].data(),"started" );

    if(system(bored2)!=0) {             
        cout << "ouch: " << bored2 << "\n\n\n";               
        exit(1);                
    }

    fflush(stdout);

    time_t t = time(0); struct tm * tm_s = localtime(&t); char blah[1024]; strftime(blah,1024,"%c",tm_s); fputs(blah,log); fputc('\n',log);
    fputs(bored2,log);
    fputc('\n',log);

    // fputs(bored1,log);
    fputc('\n',log);
    fputs(FER3,log);
    fputc('\n',log);
                
}

void release_single_rg(
    FILE * const log,
    // rarp::LIST const & pieces,
    rarp::LIST const & SEs,
    std::string const & bits,
    std::string const & dir,
    std::string const & final_bam,
    std::string const & final_bai,
    std::string const & final_dups,
    std::string const & final_checkpoint,
    rarp::NLIST & dsm
) {

    using namespace std;

    // this is very bad!?! same as in merge other than for how we get values - don't go this!?!
    string core_q = db::get_core_query(dsm["experiment_id"]);
    cout << "core_q= " << core_q << "\n";
    core::Core core_info(core_q.data());

    { float capmean=-0.0; // ,capmedian=-0.0;
    if(dsm["sample_type"]=="Genome_As_Fake_Exome") {

        // cout << "this should NEVER happen\n";
        sendmail("nb2975@cumc.columbia.edu","nb2975@cumc.columbia.edu","Single RG WGS sample","Just grab from dragen metrics");
        // cout << "we use " << dir << "/" << SEs[0] << "\n";
        string mapping_metrics = dir + "/" + SEs[0];
        assert(mapping_metrics.substr(mapping_metrics.length()-4,4)==".bam");
        mapping_metrics=mapping_metrics.substr(0,mapping_metrics.length()-4);
        mapping_metrics+=".mapping_metrics.csv";
        assert(isregfile(mapping_metrics.data()));
        // cout << "we have " << mapping_metrics << "\n";
        capmean=atof(get_single_line_output_as_string("grep 'MAPPING/ALIGNING SUMMARY,,Average sequenced coverage over genome' %s "
          " | cut -f4 -d, ",mapping_metrics.data()).data());
        // cout << "has coverage = " << capmean << "\n";

        if( !core_info.hits_coverage("Genome",capmean) ) {

            do_void_thing("delete from dragen_sample_metadata where pseudo_prepid = %s",dsm["experiment_id"].data());
            do_void_thing("delete from dragen_pipeline_step where pseudo_prepid = %s",dsm["experiment_id"].data());
            do_void_thing("update Experiment set is_released = 'not_released', merge_metrics_capturemean = %0.2f where id = %s",capmean,dsm["experiment_id"].data());
            // do_void_thing("update dragen_sample_metadata set is_merged = 200001 where pseudo_prepid = %s",dsm["experiment_id"].data());
            // do_void_thing("update Experiment set is_released = 'release_rejected' where id = %s",dsm["experiment_id"].data());
            do_void_thing("update prepT set status = 'Release_Rejected_Single_RG_UnderCov(%0.2f)', status_time = unix_timestamp() where experiment_id = %s",capmean,dsm["experiment_id"].data());

            return;

        }else do_void_thing("update Experiment set merge_metrics_capturemean = %0.2f where id = %s",capmean,dsm["experiment_id"].data());

    }else {} 
   //sendmail("nb2975@cumc.columbia.edu","nb2975@cumc.columbia.edu","Single RG sample","Must implement capture metrics for this - though RG metrics should avoid premature release"); 
  }

    // "internal genomes should NEVER get here!?! 
    // if we start doing single RG exomes then need to implement coverage checks
    // "we seem to have some strange WGS so need to implement coverage checks here - for now ignore exome capture?!? "

    char bored3[1024],bored4[1024],bored5[1024];

    sprintf(bored3,"mv %s %s",(dir+SEs[0]).data(),final_bam.data());
    sprintf(bored4,"mv %s %s",(dir+SEs[0]+".bai").data(),final_bai.data());
    sprintf(bored5,"mv %s %s",(dir+"logs/"+SEs[0].substr(0,SEs[0].length()-4)+".dragen.out").data(),final_dups.data());

    time_t t = time(0); struct tm * tm_s = localtime(&t); char blah[1024]; strftime(blah,1024,"%c",tm_s); fputs(blah,log); fputc('\n',log);

    // clearly make inline/macro
    if(system(bored3)!=0) { cout << "issue : " << bored3 << "\n\n\n"; exit(1); } fputs(bored3,log); fputc('\n',log);
    if(system(bored4)!=0) { cout << "iseue : " << bored4 << "\n\n\n"; exit(1); } fputs(bored4,log); fputc('\n',log);

    if(!isregfile(final_dups.data())) {
        if(system(bored5)!=0) { cout << "issue (dup mv) : " << bored5 << "\n\n\n"; exit(1); } fputs(bored5,log); fputc('\n',log);
    }else{
        cout << "final dup file already exists?!?\n";
    }

    do_void_thing( REPLACE_INTO_DPS, dsm["experiment_id"].data(), "completed" );
    do_void_thing( "update dragen_sample_metadata set is_merged = 1 where pseudo_prepid = %s", dsm["experiment_id"].data() );
    do_void_thing( "update prepT set status = 'Released_to_Pipeline', status_time = unix_timestamp() where experiment_id = %s", dsm["experiment_id"].data() );

    fputs(bits.data(),log);

    touchfile(final_checkpoint.data()); 

}

int auto_merge();
void auto_release();

void auto_release() { // int argc, char **argv) {

    using namespace std;
    string email_body;

    /* {
        rarp::NLIST x = db::get_named_row("drgdb","desc sample");
        for (rarp::NLIST::iterator i =x.begin(); i!=x.end(); i++) {
            cout << "["<<i->first<<"]\t"<<i->second<<"\n";
        }
        x.clear();
        x= db::get_named_row("seqdb","desc dragen_sample_metadata");
        for (rarp::NLIST::iterator i =x.begin(); i!=x.end(); i++) {
            cout << "["<<i->first<<"]\t"<<i->second<<"\n";
        }
        rarp::NLISTS xs;
        db::get_named_rows("drgdb",xs,"desc sample");
        for (int y=0;y<xs.size();++y) {
            for (rarp::NLIST::iterator i =xs[y].begin(); i!=xs[y].end(); i++) {
                cout << "["<<i->first<<"]\t"<<i->second<<"\n";
            }
        }
        xs.clear();
        db::get_named_rows("seqdb",xs,"desc dragen_sample_metadata");
        for (int y=0;y<xs.size();++y) {
            for (rarp::NLIST::iterator i =xs[y].begin(); i!=xs[y].end(); i++) {
                cout << "["<<i->first<<"]\t"<<i->second<<"\n";
            }
        }
    } */

    cout << "\nCheck merging\n";

    char const * what = " in (0,-1) ";
    if(opts::force) { // argc==2 && strcmp(*(argv+1),"force")==0) {
        // cout << "will re-run all blocked experiment_ids\n";
        do_void_thing("update dragen_sample_metadata set is_merged = 0 where is_merged %s"," <= -3 ");
        // what = " <= -2 ";
        // ALL=true;
    }

    // time_t timey_start = time(0);

    // get rid of this junk
    FILE * log = fopen(SCRIPT_DIR "/SEs.txt","a");
    if(!log) cout << "problem with " << SCRIPT_DIR << "\n", exit(1);

    rarp::NLIST count = db::get_named_row("seqdb","select count(1) num from dragen_sample_metadata where is_merged %s",what);
// cout << "wt?!? = " << what << "\n";
    int county = atoi(count["num"].data());

    if(county==0) {
        cout << "there's no samples to merge\n";
        return;
    }

    cout << "Checking merge status of " << county << " samples\n";

    { rarp::NLISTS experiments;
    db::get_named_rows("seqdb",experiments,"select d.*,expt.sample_id expt_sample_id, expt.sample_type expt_sample_type,expt.exomeKit expt_exomeKit, "
      "min(p.prepid) min_prepid, max(p.prepid) max_prepid "
      "from dragen_sample_metadata d join Experiment expt on d.experiment_id=expt.id "
      "join prepT p on p.experiment_id=d.experiment_id "
      "where is_merged %s group by d.experiment_id",what);
    fflush(stdout);

    // Popen ns("",16*1024,"r"); 
cout << "will process " << experiments.size() << " experiments\n";
    
    for (unsigned int expt_i = 0; expt_i < experiments.size(); ++expt_i ) { 
    // while( *( z=ns.getline() ) != '\0') {  char * z = 0;

        bool run = true;

        stringstream issues;
        string d_sample_name, d_type, d_flavour;

        // rarp::LIST pieces; tokenise(pieces,z,'\t');
        rarp::NLIST & dsm = experiments[expt_i];

cout << "\n\n>>>>> " << expt_i << " : " << dsm["sample_name"] << "\n\n";
for (rarp::NLIST::iterator it = dsm.begin(); it!=dsm.end(); it++ ) cout << " dsm["<< it->first <<"] "<< it->second << "\n"; // for (unsigned i =0; i<dsm.size(); i++) 
        // fflush(stdout);

        if(dsm["component_bams"]=="NULL") {
            // issues << "there's no component bam list'\n";
            update_status(dsm,NULL_COMPONENTS,email_body,issues.str());
            continue; // run=false;
        }
        
        bool igm = false;

        /* if(prept["externalData"]=="NULL") igm=true;
        else if(prept["externalData"]!="1") {
            update_status(dsm,EXTERNAL_STATUS_SCREWED,email_body,issues.str());
            continue; 
        } */

        string dir = "/nfs/" + dsm["seqscratch_drive"] + "/ALIGNMENT/BUILD37/DRAGEN/";

        if(dsm["sample_type"]=="Genome") {

            // cout << "DSM IS Genome still!?!\n";
            // fflush(stdout);

            // for (rarp::NLIST::iterator it = dsm.begin(); it!=dsm.end(); it++ ) cout << "dsm["<< it->first <<"] "<< it->second << "\n"; 
            // cout << "need to update Genome\n";

            string olddir = dir + "GENOME/" + dsm["sample_name"] + "." + dsm["pseudo_prepid"];
            string newdir = dir + "GENOME_AS_FAKE_EXOME/" + dsm["sample_name"] + "." + dsm["pseudo_prepid"];

            if(!isdir(olddir.data())) {
                // cout << "CURRENT DIR IS MISSING!?! " << olddir << "\n";
                update_status(dsm,GAFE_NO_DATA,email_body,issues.str());
                continue; 
                // run=false;
            }

            if(isdir(newdir.data())) {
                // cout << "NEW DIR IS ALREADY PRESENT\n";
                assert(strchr(newdir.data(),' ')==0);
                assert(strchr(newdir.data(),'*')==0);
                string wipe_it = "rm -rf " + newdir + "/";
                // cout << "will wipe with " << wipe_it << "\n";
                if(system(wipe_it.data())!=0) {
                    cout << "unable to clear old dir - need to do that as these may be topups?!?\n";
                    exit(1);
                }
                // continue;

                /* update_status(dsm,GAFE_ALREADY_HAS_DATA,email_body,issues.str());
                continue; // run=false; */
            }

            if(rename(olddir.data(),newdir.data())) {
                // issues << "unable to move directory " << olddir << " to " << newdir << "\n";
                update_status(dsm,MOVE,email_body,issues.str());
                continue; // run=false;
            }

            do_void_thing("update dragen_sample_metadata set sample_type = 'Genome_As_Fake_Exome', capture_kit = 'Roche' where pseudo_prepid = %s",dsm["pseudo_prepid"].data());
            dsm["sample_type"]="Genome_As_Fake_Exome";
            dsm["capture_kit"]="Roche";

            if(!isdir(newdir.data())) cout << "NEW DIR IS NOT PRESENT\n",exit(1);

            // cout << "check " << olddir << "\nand " << newdir << "\n";

        }

        if(dsm["sample_type"]=="Genome_As_Fake_Exome") {
            // this was never the place to do this 
            // cout << "DSM IS Genome_As_Fake_Exome - not wiping dps&dqm!?!\n";
            // fflush(stdout);
            /*  NO, this needs to be applied before else you permanently block!?!
            // do_void_thing("delete from dragen_pipeline_step where pseudo_prepid = %s",dsm["pseudo_prepid"].data()); */
        }

        if(dsm["sample_type"]=="Exome") dir+="EXOME/";
        // else if(dsm["sample_type"]=="Genome") dir+="GENOME/";
        else if(dsm["sample_type"]=="Genome") {
            // dir+="GENOME/";
            // issues << "NOT RUNNING GENOME DIRECTLY ATM!?! : " << dir << "\n";
            update_status(dsm,WGS,email_body,issues.str());
            continue; // run=false;
        }else if(dsm["sample_type"]=="Genome_As_Fake_Exome") {
            dir+="GENOME_AS_FAKE_EXOME/";
            if(dsm["expt_sample_type"]!="Genome") {
                // issues << "NOT RUNNING INCONSISTENT SAMPLE!?! : " << dir << "\n";
                update_status(dsm,INCONSISTENT,email_body,issues.str());
                continue; // run=false;
            }

        }else if(dsm["sample_type"]=="Custom_Capture") dir+="CUSTOM_CAPTURE/";
        else {
            cout << "what is " << dsm["sample_type"] << "\n";
            exit(1);
        }

        dir+=dsm["sample_name"]+"."+dsm["pseudo_prepid"]+"/";

        if(!isdir(dir.data())) {
            // issues << "scratch dir " << dir << " doesn't exists\n";
            update_status(dsm,NO_DIR,email_body,issues.str());
            continue; // run=false;
        }

        string rg_from_lims;

        if(igm) {

            rarp::NLISTS nrows;
            // cout << "TYPE= "<< dsm["sample_type"] << "\n\n";

            /* static */ char const * GET_RGS_BY_EXPT = "select l.prepID,l.poolID,seqID,f.FCillumID,l.FCID, "                                                                               
                ///// this is about to go!?!
                "p.adapterLet, "
                "l.LaneNum,f.Machine,f.MachineType,f.ChemVer "
                "from Lane l join Flowcell f on f.fcid=l.fcid join prepT p on p.prepid=l.prepID where experiment_id = %s "
                " and f.complete = 1 and f.pipelineComplete = 1 and f.fail = 0 and l.failr1 is null and l.failr2 is null "
                " order by prepID asc ";

            db::get_named_rows("seqdb", nrows, GET_RGS_BY_EXPT, dsm["experiment_id"].data() );

            /* this is merged and adapterlet will be moved
            if(dsm["sample_type"]=="Exome" || dsm["sample_type"]=="Custom_Capture") {
                get_named_rows(     nrows,       query::pooled_rgs,      prept["prepID"].data()      );
            }else if(dsm["sample_type"]=="Genome_As_Fake_Exome") {
              get_named_rows(   nrows,     query::nonpooled_rgs,   prept["prepID"].data()  );
            } else cout << "what off\n", exit(1); */

            // cout << "FCBC: got " << nrows.size() << " registered lanes for internal bam\n";
            std::set<string> c;
            for(unsigned r=0; r<nrows.size(); ++r){
                // for (rarp::NLIST::iterator it = nrows[r].begin(); it!=nrows[r].end(); it++ ) { cout << "nrows["<< it->first <<"] "<< it->second << ", "; } cout << "\n";
                c.insert(nrows[r]["FCillumID"]+"."+nrows[r]["LaneNum"]);
            }

            if(nrows.size()==0) {

                if(atoi(dsm["min_prepid"].data())>30000) {

                    issues << "WON'T RUN NEWISH SAMPLE WITHOUT FCBC INFO\n";
                    // cout << "WON'T RUN NEWISH SAMPLE WITHOUT FCBC INFO\n";
                    update_status(dsm,MISSING_FCINFO,email_body,issues.str());

                    if(dsm["sample_type"]!="Exome") {
                        // cout << "THIS DOESN'T SEEM RIGHT!?! : MUST RE-CHECK JOIN FOR NON-HYB-POOLED FCBC\n";
                    }

                    continue; // run=false;

                }else{
                    // cout << "RUNNING OLD SAMPLE WITHOUT FCBC INFO\n";
                }

            }else{

                for (std::set<string>::iterator it = c.begin(); it!=c.end(); it++)
                if(rg_from_lims.empty()) rg_from_lims=*it;
                else rg_from_lims += ',' + *it;

            }

            // cout << "have rgs = " << rg_from_lims << "\n";
            // issues << "fcid= " << pieces2[4] << "\nlane="<< pieces2[7] << "\nbc= " << pieces2[6] << "\nmachine= " << pieces2[8] << "\nmtype= " << pieces2[9] << "\nversion= " << pieces2[10] << "\n";

        }

        if(dsm["component_bams"].length()>=5000) cout << "WARN: list length may overflow\n\n";

        if(dsm["component_bams"].substr(dsm["component_bams"].length()-4)!=".bam"){ 

            issues << "ERROR: this is f'd\n\n";
            update_status(dsm,NOT_A_BAM,email_body,issues.str());
            continue; // run=false;
        }

        // cout << "checking '" << dsm["sample_name"] << "'\n";

        string final_bam = dir +dsm["sample_name"]+"."+dsm["pseudo_prepid"]+".bam", 
          final_bai = dir +dsm["sample_name"]+"."+dsm["pseudo_prepid"]+".bai",
          final_dups = dir +"logs/"+dsm["sample_name"]+"."+dsm["pseudo_prepid"]+".dragen.out",
          final_checkpoint = dir +dsm["sample_name"]+"."+dsm["pseudo_prepid"]+".merged.marked.done",
          // vial re-declaration - about to get even more vial due to issue with flash space
          merge_intermediate = dir+dsm["sample_name"]+"."+dsm["pseudo_prepid"]+".merge.bam";

        string gg= (string(SCRIPT_DIR)+dsm["sample_name"]+"."+dsm["experiment_id"]+".sh");

        if(dsm["is_merged"]=="0"){// || ALL) { 

            // cout << "check " << dsm["is_merged"] << "\n";
            if(isregfile(gg.data())) unlink(gg.data());
            if(isregfile(final_checkpoint.data())) unlink(final_checkpoint.data());

            if(isregfile(final_checkpoint.data())) {
                // update_status(dsm,CHECKPOINT_EXISTS,email_body,issues.str());
                // continue; // run=false;
            }else if(isregfile(gg.data())) {
                issues << "There's something nasty going on " << gg << "\n";
                // update_status(dsm,SCRIPT_EXISTS,email_body,issues.str());
                // continue; // run=false;
            }

            assert(dsm["component_bams"].substr(dsm["component_bams"].length()-4)==".bam"); 

            rarp::LIST SEs; tokenise(SEs,dsm["component_bams"],',');

            string sn, fc;
            // cout << "there are " << SEs.size() << " SEs to process\n";

            char size_s[1024];
            sprintf(size_s,"%lu",SEs.size()); 
            string bits="", merge_cmd="";

            for (unsigned i2 = 0; i2<SEs.size() ; ++i2) {

                rarp::LIST lazy; tokenise(lazy,SEs[i2],'.');
                if(lazy.size()==6) {
                    issues << "ERROR: SE name doesn't conform '" << SEs[i2] << "'\n\n";
                    update_status(dsm,COMPONENT_ERROR,email_body,issues.str());
                    run=false; 
                    continue;
                }

                string se = dir+SEs[i2];

                if(!isregfile(se.data())) {
                    issues << "ERROR: SE bam has issue : " << se << "\n";
                    update_status(dsm,COMPONENT_ERROR,email_body,issues.str());
                    run=false; 
                    continue;
                }

                merge_cmd+=" I="+se;

                if(!isregfile((se+".bai").data())) {
                    issues << "ERROR: SE bai has issue : " << se << "\n";
                    update_status(dsm,COMPONENT_ERROR,email_body,issues.str());
                    run=false; 
                    continue;
                }

                // string index = se.substr(0,se.length()-3)+"bai";
                bits+=dsm["min_prepid"]+"\t"+dsm["experiment_id"]+"\t"+size_s+"\t" 
                  + dsm["sample_name"]+"\t"+dsm["sample_type"]+"\t"+dsm["capture_kit"]+"\t"+dsm["seqscratch_drive"]+"\t";

                { char arsv2[1024];
                FILE * f = fopen(se.data(),"r");
                if(!f) { cout << "wt?!?!\n"; 
                    run=false; 
                    continue;
                }else{ fseek(f,0L,SEEK_END); sprintf(arsv2,"%ld",ftell(f)); fclose(f); }
                bits+=SEs[i2]+"\t"+arsv2+"\t"; }

                bool okay = false, error = false, warn = false;
                char RGC[1024];

                { char tmp[10240]; sprintf(tmp,"/nfs/seqscratch09/dsth/vcf_thing/bam_thing %s",se.data()); Popen ns1(tmp,16*1024,"r");
                    char * z2;
                    while( *( z2=ns1.getline() ) != '\0') { 
                        if(memcmp(z2,"ERROR",strlen( "ERROR") )==0) error=true;
                        if(memcmp(z2,"OKAY",strlen( "OKAY") )==0) okay=true;
                        if(memcmp(z2,"DODGY_OLD_BAM",strlen( "DODGY_OLD_BAM" ))==0) warn=true;
                        if(::isdigit(z2[0])) strcpy(RGC,z2);
                        bits+=z2;
                        bits+="\t";
                    } }
                bits+="\n";

                if(!okay&&!error) {
                    issues << "ERROR: (" << se << ") THERE'S A PROBLEM CHECKING THE BAM\n"; run=false; 
                    update_status(dsm,COMPONENT_ERROR,email_body,issues.str());
                    run=false; 
                    continue;
                } else if(error) {
                    issues << "ERROR: (" << se << ") THE BAM HAS ISSUES\n"; run=false; 
                    update_status(dsm,COMPONENT_ERROR,email_body,issues.str());
                    run=false; 
                    continue;
                } if(warn) cout << "WARN: THE BAM LOOKS LIKE IT'S PRETTY OLD?!?\n";


                // cout << "have rg_counts " << RGC << "\n";
                rarp::LIST rgs; tokenise(rgs,RGC,':');
                assert(rgs.size()==3);
                if(rgs[0]!="1"||rgs[2]!="1") {

                    std::string rn=get_single_line_output_as_string("samtools view %s | head -n1 | cut -f1",se.data()),
                      fc, bc, mtype_from_fc, mtype_from_pfx;
                    rn='@'+rn;
                    bool is18=fastq::rncheck(rn.data(),fc, bc, mtype_from_fc, mtype_from_pfx);
                    // cout << "1.8="<<is18<<"\n";
                    int pid = atoi(dsm["max_prepid"].data());
                    // cout << "we have prepid " << pid << " and " << igm << "\n";
                    if(pid<50000) cout << "it's really old so whatever...\n";
                    else if(!igm) cout << "i'm sick of dealing with external so whatever...\n";
                    else if(igm&&is18) { // }else if(igm||is18) {
                        issues << "ERROR: (" << se << ") COMPONENT HAS RG PROBLEM " << RGC << "\n";
                        update_status(dsm,COMPONENT_RG_ERROR,email_body,issues.str());
                        run=false; 
                        continue;
                    }else cout << "it's external and dodgy format so whatever...\n"; 

                }

            } 

            // RUNIT:
            if(run) { 

                // cout << "WILL_RUN " << dir << "\n";
                // bool done = false;
                rarp::NLIST nc = db::get_named_row("seqdb", "select count(1) as step1count from dragen_pipeline_step where pipeline_step_id = 1 and pseudo_prepid = %s",dsm["experiment_id"].data() ),
                  nc2 = db::get_named_row("seqdb", "select count(1) as qccount from dragen_qc_metrics where pseudo_prepid = %s",dsm["experiment_id"].data() ); 

                if(nc.count("step1count")==0) cout << "issue : " << __LINE__ << "\n",exit(1);
                if(nc2.count("qccount")==0) cout << "issue : " << __LINE__ << "\n",exit(1);

                // pure evil 
                { string merge_intermediate = dir+dsm["sample_name"]+"."+dsm["pseudo_prepid"]+".merge.bam",
                  final_bam = dir +dsm["sample_name"]+"."+dsm["pseudo_prepid"]+".bam";
                if(isregfile(merge_intermediate.data())) unlink(merge_intermediate.data());
                if(isregfile(final_bam.data())) unlink(final_bam.data()); }

                // this is the one line thing?!?
                char tmp[1024]; sprintf(tmp,"ls %s/*.bam | grep -v cram2bam | wc -l",dir.data()); Popen ns1(tmp,16*1024,"r"); char * z2=ns1.getline();
                // cout <<"SE count = " << size_s << " there are " << z2 << " bams in dir\n";

                if(nc2["qccount"]!="0") {
                    issues << "ERROR: This sample seems to have finished pipeline before " << dsm["experiment_id"] << " got = '" << nc2["qccount"]
                      << "' (" << dsm["sample_type"] << ")\n";
                    update_status(dsm,DQM,email_body,issues.str());
                }else if(nc["step1count"]!="0") {
                    issues << "ERROR: This sample seems to have been mapped/merged before - has step 1 in dps!?! " << dsm["experiment_id"] << " got = '" << nc["step1count"]
                      << "' (" << dsm["sample_type"] << ")\n";
                    update_status(dsm,DSP,email_body,issues.str());
                    // continue;
                }else if(strcmp(size_s,z2)!=0) {
                    issues << "bam count doesn't match. won't run - if this appears to be okay reset to -1\n";
                    update_status(dsm,BAM_COUNT_MISMATCH,email_body,issues.str());
                    run=false; 
                    continue;
                }else if(SEs.size()==1) {
                    // cout << "there's nothing to do\n";
                    if( check_final_bam( dsm,rg_from_lims,dir+SEs[0],igm,email_body,false /*, prepidi, prept */ ) )
                       release_single_rg( log, /* pieces, */ SEs,bits,dir,final_bam,final_bai,final_dups,final_checkpoint , dsm);

                }else merge_multiple(final_bam,final_checkpoint,log,merge_cmd,dir,gg,merge_intermediate,dsm);
                // fflush(stdout);

            } else {
                issues << "ERROR: WILL_NOT_RUN_PROBLEM : " <<  dsm["experiment_id"] << " : " <<  dir << "\n";
                update_status(dsm,JUNK,email_body,issues.str());
            }

        }else{

            if(isregfile(final_checkpoint.data())) 
              release_merged_rgs( dsm,rg_from_lims, /* pieces, */ log,final_bam,final_checkpoint,igm,email_body,merge_intermediate,gg,issues,dir /* ,prepidi,prept */ );
            // else {
                // cout << "HMM: This seems to be running - must re-check qstat/when was started... - also make script exit with any issues at all and check qacct return?!?\n";
            // }
        }
        // fflush(stdout);
    }

    fclose(log);
    // cout << "stf-c?!?\n"; cout.flush();
    // cout << "stf-d?!?\n"; cout.flush();

    }

    if(!email_body.empty()) {
        // sendmail("dh2880@columbia.edu","Daniel Hughes <dh2880@columbia.edu>","Blocking problem bams",email_body.data());
        cout << "Blocking samples:\n"<<email_body<<"\n\n";
    }

    // exit(1);

    return;
}

// #endif

int auto_merge() { // int argc, char **argv) {
    
    using namespace std;

    string userids; 
    { rarp::NLISTS lims_userid;
    db::get_named_rows("seqdb",lims_userid,"SELECT userid FROM users WHERE NETID= '%s'",opts::myuser.user());
    // db::get_named_rows("seqdb",lims_userid,"SELECT userid FROM users WHERE NETID= '%s'",getenv("USER"));
    if(lims_userid.size()==1) userids=lims_userid[0]["userid"];
    else if(lims_userid.size()>1) {
        // cout << "you're userid is ambiguous. will use the first one\n";
        userids=lims_userid[0]["userid"];
    }else assert(lims_userid.size()==1); }

    // char const * serv = opts::serv;

    string file, st, capture_kit;
    int release_time = 0, rel = 0, will = 0, wont = 0;

    if(opts::commit) cout << "WARNING: Will release samples\n";
    else cout << "WARNING: This is a dry run\n";

    rarp::NLISTS cun2;
    db::get_named_rows("seqdb",cun2,Q1); // get_named_rows(cun2,Q1,argv[1]); 
    string message, cur_pool="";
    int pool_count=1;

    // cout << "there are " << cun2.size() << " samples to check for release eligibility.\n";

    std::map<int,FUNKY> list;
    std::vector<float> pool;
    for (unsigned int i = 0 ; i < cun2.size(); ++i ) {

        // cout << "\n-------------------\n\n";
        if(cun2[i]["l_capmean_sum"]=="0") {
            // cout << "I don't do legacy samples - use 'release' procedure or 'rg metrics' page\n";
            continue;
        }

        // cout << "STATUS_OF_RGS: previous= " << cun2[i]["CURRENT_RGS_STAMP"] << " new= " << cun2[i]["NEW_RGS_STAMP"] <<"\n";
        // cout << "> sample ["<<i<<"]\n";
        // this is basically pointless now that they want brute force release!?!
        if(cur_pool!=cun2[i]["pool_name"]) {

            cur_pool=cun2[i]["pool_name"];
            // cout << "previous pool count " << pool_count << "\n";
            // assert(pool_count==1);
            pool_count=atoi(cun2[i]["pool_count"].data());
            // cout << "\nNEW_POOL " << cur_pool << " of " << pool_count << "\n";

            if(pool.size()) {

                std::sort(pool.begin(),pool.end(),myfunction);
                // float a = pool[int(0.25*pool.size())], b = pool[int(0.75*pool.size())], c = b - a;

                // add in mean and median - perhaps ratio?!? but certainly anything that is outside of inner fences?!?
                // cout << "we have " << a << " - " << b << " : " << c << "\n"; // cout << " we have " << pool.size() << " : ";
                // cout << "lower_inner_fence=" << (a-(1.5*c)) << " - upper_inner_fence=" << (b+(1.5*c)) << "\n";
                // cout << "lower_inner_fence=" << (a-(1*c)) << " - upper_inner_fence=" << (b+(1*c)) << "\n";
                // for (unsigned int j = 0 ; j < pool.size() ; ++j) cout << pool[j] << ", "; cout << "\n";
                pool.clear();

            }

        }else{
            //cout << "same pool " << --pool_count << "\n";
        }

        pool.push_back(atof(cun2[i]["l_capmean_sum"].data()));

        // cout << "we check project info here\n";
        // cout << "running " << cun2[i]["CORE_QUERY"] << "\n";
        core::Core core_info(cun2[i]["CORE_QUERY"].data());
        string stamp_current = cun2[i]["CURRENT_RGS_STAMP"],
          stamp_new = string(core_info.is_releasable()?"+:":"-:")+cun2[i]["NEW_RGS_STAMP"];
        // cout << "stamp_current="<<stamp_current<<"\nstamp_new="<<stamp_new<<"\n\n";
if(strcmp(getenv("USER"),"dh2880")==0) {
          cout << "Sample info:\n";
          for (rarp::NLIST::iterator it=cun2[i].begin(); it!=cun2[i].end(); it++) { cout << " ["<<it->first << "]\t"<<it->second <<"\n"; }
}
        // for (;;) {
        // fflush(stdout);

        if(  cun2[i]["dsm_experiment_id"]!="NULL"  || cun2[i]["dqm_experiment_id"]!="NULL"  || cun2[i]["e_is_released"]!="not_released"  || cun2[i]["w_experiment_id"]!="NULL" ) {

            string borederer;
            { if(cun2[i]["dsm_experiment_id"]!="NULL") borederer+="sm";
            if(cun2[i]["dqm_experiment_id"]!="NULL") borederer+=":qm"; 
            if(cun2[i]["e_is_released"]!="not_released") borederer+=":eir";
            if(cun2[i]["w_experiment_id"]!="NULL") borederer+=":we";
            if(borederer[0]==':') borederer=borederer.substr(1,borederer.length()-1); }

            rel++;
            // cout << "dsm_experiment_id= " << cun2[i]["dsm_experiment_id"] << "\n";
            // cout << "dqm_experiment_id= "<< cun2[i]["dqm_experiment_id"]<< "\n";
            // cout << "e_is_released " << cun2[i]["dsm_experiment_id"] << "\n";
            // cout << "w_experiment_id= " << cun2[i]["w_experiment_id"] << "\n";
            // cout << "\n\n\tthis has been released already!?!\n\n";
            char msg[2048];
            sprintf(msg,"update prepT set status = 'This requires full deprecation and re-release (synch samples:%s)' where ( failedprep=0 or failedprep>=100) and experiment_id = %s ; select row_count() updated ",borederer.data(),cun2[i]["e_experiment_id"].data());
            // cout << "using " << msg << "\n";
            rarp::NLIST x = db::get_named_row("seqdb",msg);
            // cout << "got " << x["updated"]<<"\n";
            continue;
            assert(0);

        }else{

            // cout << "\n\n\tchecking yields!?!\n\n";
            float terrible_mb_measure=atof(cun2[i]["l_lane_sum"].data()), rg_mean_sums=atof(cun2[i]["l_capmean_sum"].data());

            // cout << "terrible_mb_measure= " << terrible_mb_measure << ", rg_mean_sums= " << rg_mean_sums << "\n\n";
            // << " : thresholds are " << opts::wes_min << "," << opts::wgs_min << "\n\n";
            // cout << "rg_mean_sums " << terrible_mb_measure << " : ";
            // assert(strstr(cu[i]["rg_metrics_statuses"].data(),"contaminated")==0);
            // assert(strstr(cu[i]["rg_metrics_statuses"].data(),"stalled")==0);
            // assert(strstr(cu[i]["rg_metrics_statuses"].data(),"pending")==0);
            // cout << "REL_SUM= ";

            // horrible, brute force release!?!
            // cout << "This is from " << core_info.useful_name() << "\n";
            if(
                cun2[i]["sum_rg_statuses"]!="fastq_archived:okay"
                && cun2[i]["sum_rg_statuses"]!="fastq_archived:contaminated,fastq_archived:okay"
            ) {
                // cout << "this isn't ready for release - i.e. all useable RGs should be fastq_archived & okay!?!\n";
            }else if(!core_info.is_releasable()) { // if(!core_info.is_approved()) {
                // cout << "this is not an approved project!?! - blocking " << cun2[i]["sample_internal_name"] << " - " << cun2[i]["e_experiment_id"] << "\n";
                char msg[2048];
                sprintf(msg,"update prepT set status = 'Error - %s is not an approved subproject' where ( failedprep=0 or failedprep>=100) and experiment_id = %s ; select row_count() updated ",core_info.useful_name(),cun2[i]["e_experiment_id"].data());
                // cout << "using " << msg << "\n";
                rarp::NLIST x = db::get_named_row("seqdb",msg);
                // cout << "updated="<<x["updated"]<<"\n";
                // exit(1);
            }else if( stamp_current==stamp_new ) {
                // cout << "SAME_RGS : previous= " << stamp_current << " new= " << stamp_new <<"\n";
                // cout << "we must check that approval status hasn't changed\n";
            }else if( core_info.hits_coverage(cun2[i]["sample_type"],rg_mean_sums)

              /* (cun2[i]["sample_type"]=="Exome" && (rg_mean_sums>opts::wes_min)) // (cun2[i]["sample_type"]=="Exome" && (rg_mean_sums>(opts::wes_min*1.05)))
              (cun2[i]["sample_type"]=="Exome" && (terrible_mb_measure>8000||rg_mean_sums>90.0)) ||  (cun2[i]["sample_type"]=="Genome" && (rg_mean_sums>opts::wgs_min))  */
          
              /* ||  ( cun2[i]["sample_type"]=="Genome" && cun2[i]["sample_internal_name"].substr(0,7)=="chapscz" && terrible_mb_measure>81000 ) */

              ) { 

                // *1.05))) ) { // terrible_mb_measure>115000 ) ) { // for this 30x is fine!?!||rg_mean_sums>82.5))
                // (cun2[i]["sample_type"]=="Genome" && (rg_mean_sums>(opts::wgs_min*1.05))) ) { // terrible_mb_measure>115000 ) ) { // for this 30x is fine!?!||rg_mean_sums>82.5))
                // (cun2[i]["sample_type"]=="Genome" && terrible_mb_measure>115000 ) ) { // for this 30x is fine!?!||rg_mean_sums>82.5))
                // cout << "\n\t\tWILL RELEASE " << cun2[i]["sample_type"] << " : " << rg_mean_sums << "\n\n";
                
                will++;
                FUNKY f;
                f.experiment_id=cun2[i]["e_experiment_id"], f.sample_name=cun2[i]["sample_internal_name"], f.sample_type=cun2[i]["sample_type"],
                  f.capture_kit=cun2[i]["exomekit"], f.priority="1", f.is_external="0",
                f.prepid=cun2[i]["prepid"], // f.sample_id=cun2[i]["sample_id"],
                f.end_point=cun2[i]["end_point"];
                list.insert(make_pair(atol(cun2[i]["e_experiment_id"].data()),f)); 

                { char arsv2[1024];
                sprintf(arsv2,"update Experiment set rgs = '%s' where id = %s ; select row_count() updated ",    stamp_new.data(), cun2[i]["e_experiment_id"].data());
                // cout << arsv2 << "\n";
                rarp::NLIST x = db::get_named_row("seqdb",arsv2);
                // cout << "updated="<<x["updated"]<<"\n";
                assert(x["updated"]!="0"); }

            }else{
                assert(
                    strstr(cun2[i]["sample_type"].data(),"Exome")==0
                    || strstr(cun2[i]["sample_type"].data(),"Genome")==0
                );

                wont++;
                // cout << "\n\t\tWON'T RELEASE '" << cun2[i]["sum_statuses"] << "'\n\n";
                char arsv2 [2048];
                if(terrible_mb_measure<500) sprintf(arsv2,"update prepT set status = 'Error - Requires HTS investigation' where ( failedprep=0 or failedprep>=100) and experiment_id = %s ; select row_count() updated ",cun2[i]["e_experiment_id"].data());
                else sprintf(arsv2,"update prepT set status = 'Not eligible for autorelease - Requires manual HTS release or topup (yield=%.02f/%.02f)' where ( failedprep=0 or failedprep>=100) and experiment_id = %s ; select row_count() updated ",terrible_mb_measure,rg_mean_sums,cun2[i]["e_experiment_id"].data());
                // cout << "using " << arsv2 << "\n";
                rarp::NLIST x = db::get_named_row("seqdb",arsv2);
                // cout << "updated="<<x["updated"]<<"\n";
                // assert(x["updated"]=="1" || x["updated"]=="0" || x["updated"]=="2");
                // for (rarp::NLIST::iterator it = cun[i].begin(); it!=cun[i].end(); it++) { cout << "["<<it->first << "] "<< it->second << "\n"; }

            }

            /* int yield=row[1]?atoi(row[1]):0;
            // cout << "\tyield= " << yield << "\n";
            if(st=="Genome" && yield<99000) cout << " > WARNING: WGS Sample has a low yield (" << yield << ")\n";
            else if(st=="Exome" && yield<10000) cout << " > WARNING: WES Sample has a low yield (" << yield << ")\n";
            else if(st=="Custom_Capture" && yield<1000) cout << " > WARNING: Custom Sample has a low yield (" << yield << ")\n";
            if(strcmp(row[9],"Storage")!=0) {
                cout << " > WARNING: Sample doesn't have status 'Storage' - current value is '" << row[9] << "'\n";
            }
            if(row[11]) {
                cout << " > WARNING: Sample has entries in DragenDB (" << row[11] << ")\n";
            } */

        }
    
    }

    /* auto-close tickets * presumably reason new RGs are greyed out is that they haven't go flowcell.complete - i.e. the silly button hasn't been pressed yet!?!? */

    //// wt - this is part of the legacy procedure. needs to be factored?!?
    MYSQL *con = mysql_init(NULL);
    if (con == NULL) fprintf(stderr, "mysql_init() failed\n"), exit(1);
    
    if (mysql_real_connect(con, opts::myuser.host(), opts::myuser.user(), opts::myuser.pass(), "sequenceDB", 0, NULL, 0) == NULL)
      finish_with_error(con);
    
    { stringstream x; x << "Will release " << list.size() << " samples.\n"; message+=x.str(); }
    
    // whatever, can't be bothered to clean up...?!?
    if(list.size()==0) return 0;

    int relc=0;
    for (std::map<int,FUNKY>::iterator it=list.begin(); it!=list.end(); it++) {

        assert(it->first==atoi(it->second.experiment_id.data()));

        { stringstream x;
        x << "Releasing expt_id=" << it->first;
        x << ", sample_name="<< it->second.sample_name;
        x << ", sample_type="<< it->second.sample_type;
        x << ", capture_kit="<< it->second.capture_kit;
        x << ", priority="<< it->second.priority;
        x << ", is_external= "<< it->second.is_external << "\n";
        // cout << x.str();
        message+=x.str(); }

        bool warnings=false;

        if (mysql_query(con, "START TRANSACTION")) finish_with_error(con);    

        char fer3[8*1024];
        sprintf(fer3,"INSERT INTO statusT (STATUS_TIME,    STATUS,     PREPID,     USERID,     POOLID,     SEQID)"
          " SELECT UNIX_TIMESTAMP(),   'Released to Bioinformatics Team', PREPID, %s, 0, 0 FROM prepT WHERE experiment_id= %s",
          userids.data(), it->second.experiment_id.data());
        if(mysql_query(con, fer3)) finish_with_error(con);    
        if(mysql_insert_id(con)==0) {
            // cout << "problem with statust insertion\n";
            warnings=true;
        }
        if(mysql_affected_rows(con)==0) { // if(mysql_affected_rows(con)!=1) {
            // cout << "problem with statust insertion(s)\n";
            warnings=true;
        }

        sprintf(
          fer3,
          "UPDATE prepT SET STATUS='Released to Bioinformatics Team', is_released = 1, released_time = ( UNIX_TIMESTAMP() + ( 24*3600*%d ) ) WHERE experiment_id= %s",
          release_time,it->second.experiment_id.data()
        );
        if (mysql_query(con, fer3)) finish_with_error(con);    

        // id = mysql_insert_id(con); // cout << "got id = " << id << "\n"; // cout << "there were " <<  mysql_affected_rows(con) << " affected rows\n";
        if(mysql_affected_rows(con)==0) cout << " > WARNING: Sample apparently already had status 'Released to Bioinformatics Team'\n";

        if(st!="RNAseq") {

            sprintf(fer3,"insert into dragen_sample_metadata "
              " (     experiment_id,    pseudo_prepid,  sample_name,    sample_type,    capture_kit,    priority,   is_external,    end_point, is_merged) "
              "values(%s,               %s,             '%s',           '%s',           '%s',           %s,         %s,             '%s', 80000)",
              it->second.experiment_id.data(), it->second.experiment_id.data(), it->second.sample_name.data(),
              it->second.sample_type.data(), it->second.capture_kit.data(), it->second.priority.data(),
              it->second.is_external.data(), it->second.end_point.data() );

            if(mysql_query(con, fer3)) finish_with_error(con);    
            // not autoincrementing pk...?!?
            // if(mysql_insert_id(con)==0) warnings=true;
            if(mysql_affected_rows(con)!=1) {
                cout << "unable to insert into dsm\n";
                warnings=true;
            }

            // id = mysql_insert_id(con);
            // assert(mysql_affected_rows(con)==1);
            // cout << "got id = " << id << "\n";
            // cout << "there were " <<  mysql_affected_rows(con) << " affected rows\n";
                  
        }

        sprintf(fer3,"update Experiment set is_released = 'released' where id=%s", it->second.experiment_id.data());
        if(mysql_query(con, fer3)) finish_with_error(con);    
        // duh : if(mysql_insert_id(con)==0) warnings=true;
        if(mysql_affected_rows(con)!=1) {
            cout << "unamble to update experiment\n";
            warnings=true;
        }

        if(opts::commit) {
            if(warnings) {
                cout << "ERROR: There was problem with the db - not commiting changes.\n";
                if (mysql_query(con, "ROLLBACK")) finish_with_error(con); 
            }else {
                if (mysql_query(con, "COMMIT")) finish_with_error(con);    
                ++relc;
            }
        }else{
            if (mysql_query(con, "ROLLBACK")) finish_with_error(con); 
        }

        // if (mysql_query(con, "COMMIT")) finish_with_error(con);    

    }

    { stringstream x; x << "Released " << relc << " samples.\n";
    cout << "already_released="<<rel<<", will_release="<<will<<", won't release="<<wont<<"\n";
    cout << x.str();
    message+=x.str(); }

    sendmail("nb2975@cumc.columbia.edu","nb2975@cumc.columbia.edu","Release summary",message.data());
    cout << "SUMMARY\n" << message << "\n\n";

    mysql_close(con);

    return 0;

}

// dynamic modules e.g. Data::Dumper...
EXTERN_C void xs_init (pTHX);
EXTERN_C void boot_DynaLoader (pTHX_ CV* cv);
EXTERN_C void xs_init(pTHX) {
	const char *file = __FILE__;
	dXSUB_SYS;
	newXS("DynaLoader::boot_DynaLoader", boot_DynaLoader, file);
}

static PerlInterpreter *my_perl;

void submit_and_post_checks(bool post_checks = false) {

    using std::cout;
    ///// should probalby pass in host,db,user,pword and list of sub?!?
    ///// no system call with pword and make changes to pword in single location - still horrible though!?!
    ///// need to build this path!?!
    // char const * mod[] = { BASE_DIR "/bin/GoButton.pm", "/home/dsth/Trabalho/DNA_PIPE/Refactor/GoButton.pm" };
    char const * mod[] = { "/nfs/goldstein/software/sequencing_pipe/master/sequencing_pipe/GoButton.pm", "/home/dsth/Trabalho/DNA_PIPE/Refactor/GoButton.pm" };
    // cout << "we have " << BASE_DIR << " and " << SCRIPT_DIR << "\n";
    char * my_argv[] = { (char*) "", (char*) (( strcmp(getenv("HOSTNAME"),"centos6")==0 )?*(mod+1):*mod) }; 
    my_perl = perl_alloc();
    perl_construct(my_perl);
    perl_parse(my_perl, xs_init, 2, my_argv, NULL);
    // char *noargs[] = { NULL };
    dSP;                            
    ENTER;                          
    SAVETMPS;                       
    PUSHMARK(SP);                   
    XPUSHs(sv_2mortal(newSVpv(opts::myuser.host(),strlen(opts::myuser.host()))));
    XPUSHs(sv_2mortal(newSVpv("sequenceDB",strlen("sequenceDB")))); 
    XPUSHs(sv_2mortal(newSVpv(opts::myuser.user(),strlen(opts::myuser.user()))));
    XPUSHs(sv_2mortal(newSVpv(opts::myuser.pass(),strlen(opts::myuser.pass()))));

    if(!post_checks) {
        // cout << "submit\n";
        sleep(3);

        XPUSHs(sv_2mortal(newSVpv("push",strlen("push")))); 

        time_t t = time(0);
        srand(t); int r = rand(), h = r%HOW_OFTEN==0; printf("> %d, %d -> %d\n",r,HOW_OFTEN,h);
        if(1) { // if(h) { 
            // silly but wanna make sure qstat is all caught up?!?...?!?
            // sleep(30);
            XPUSHs(sv_2mortal(newSVpv("push_back",strlen("push_back")))); 
        } 

    }else{
        // cout << "post_checks\n";
        sleep(3);
        XPUSHs(sv_2mortal(newSVpv("archive_md5",strlen("archive_md5")))); 
    }

    // XPUSHs(sv_2mortal(newSViv(p))); 
    PUTBACK;                        
    call_pv("GoButton::entry", G_SCALAR); // no return?!?
    SPAGAIN;                        
    PUTBACK;
    FREETMPS;                       
    LEAVE;
    perl_destruct(my_perl);
    perl_free(my_perl);
    cout << "done\n";
}

int release_manual(int, char **);// , opts::MysqlUser &);
void load(int, char **);

namespace email_bits {

    void lims_driven_summary(int, char **) {

        // bool email = true; // argc==2 && strcmp(argv[1],"email")==0;

        using namespace std;
        time_t since=1527843600, now = time(0);
        stringstream strtmp;
        for (time_t day = since, count = 0, capacity = 4; day < now - (2*24*3600); day += 24 * 3600,++count) {
            if(day>=1546336800)capacity=6;
            struct tm * tmp = localtime(&day);
            char silly[1024];
            strftime(silly,sizeof silly, "%FT%TZ", tmp);
            rarp::NLIST c = db::get_named_row("seqdb","select count(1) count,group_concat(concat(machine,':',chemver,':',fcillumid,':',date_format(dateread1,'%%d-%%H-%%i'),':',date_format(daterta,'%%d-%%H-%%i'),':',round((unix_timestamp(daterta)-unix_timestamp(dateread1))/3600,2) )) list from Flowcell where ( dateread1 < '%s' or seqtime < %llu ) and daterta > '%s'",silly,day,silly);
            // rarp::NLIST c = db::get_named_row("seqdb","select count(1) count,group_concat(concat(machine,':',chemver,':',fcillumid,':',date_format(dateread1,'%%d-%%H-%%i'),':',date_format(daterta,'%%d-%%H-%%i'),':',timediff(daterta,dateread1) )) list from Flowcell where ( dateread1 < '%s' or seqtime < %llu ) and daterta > '%s'",silly,day,silly);
            // int blarp = atoi(c["count"].data());
            // float util = 100 * blarp / capacity;
            strtmp
            << string(silly).substr(0,10) << " " 
            << ( tmp->tm_wday==0?"Sun": tmp->tm_wday==1?"Mon": tmp->tm_wday==2?"Tue": tmp->tm_wday==3?"Wed":
                    tmp->tm_wday==4?"Thu": tmp->tm_wday==5?"Fri":"Sat" ) << " " 
            << "(" << tmp->tm_yday << ")" 
            << " : " << c["count"] ;
            if(c["list"]!="NULL") strtmp << " ("<<c["list"] << ")";
            strtmp << "\n";
        }
        cout << strtmp.str();
        // sendmail("dh2880@columbia.edu,dsth@cantab.net,jb4393@cumc.columbia.edu","dh2880@columbia.edu","Flowcell Utilisation Summary",strtmp.str().data());
        if(opts::email) sendmail("igm-bioinfo@columbia.edu,jb4393@cumc.columbia.edu","igm-bioinfo@columbia.edu","Flowcell Utilisation Summary",strtmp.str().data());
        else sendmail("igm-bioinfo@columbia.edu","igm-bioinfo@columbia.edu","Flowcell Utilisation Summary",strtmp.str().data());
    }

    // #define DAYS "90"
    #define DAYS "14"

    int const update_int = 1000, MIN_DONE = 2;
    long long const AGE_LIM = 30;

    struct RUNINFO {

        std::string lims_status, machine, start_date, run_dir, fctype, cycle_rec, lims_machine, lims_machine_type, 
        lims_chemver, lims_readlens, lims_total_cycles, lims_date_bcl, lims_rundir, lims_rta_complete;
        int run_num, total_cycles, last_cycle, lims_fcid, lims_made;
        bool rta_complete, seq_complete, copy_complete, /* this is same as has RTA3.cfg * - i.e. pretty much same as has run dir...!?! : */ run_started; 
        long long lims_delete_bcl_time_set, lims_delete_bcl_time_complete, lims_s3_transfer_time_complete, lims_fc_made_seqtime,
        lims_rta_complete_epoch;
        char slot;

        RUNINFO() : run_num(0), total_cycles(0), last_cycle(0), lims_fcid(0), lims_made(0), 
        rta_complete(false), seq_complete(false), copy_complete(false), run_started(false), 
        lims_delete_bcl_time_set(0), lims_delete_bcl_time_complete(0), lims_s3_transfer_time_complete(0), lims_fc_made_seqtime(0),
        lims_rta_complete_epoch(0),
        slot('\0') {} // no need to initialise string explicity as class type will invoke default c-tor!?!?
    };

    namespace todo_move_to_shared_location {

    // #define AWS "LD_LIBRARY_PATH=/nfs/goldstein/software/python3.6.1-x86_64_shared/lib /home/dh2880/.local/bin/aws "
    #define RUNDIR "/nfs/hts/novaseq/sequence/Runs"

        // inline void check_aws_version() { 
        // inline std::string get_log_file(rarp::NLIST const & x, char const * which) {
        
        inline void check_base_rundir() { assert(isdir(RUNDIR)); } /// really just to make sure we're somehwere it's mounted??!

        inline std::string get_run_dir(rarp::NLIST const & x) {
            using namespace std;
            char tmp[1024];
            assert(x.count("rundir"));
            sprintf(tmp,"%s/%s",RUNDIR,x.find("rundir")->second.data());
            return tmp;
        }

        inline std::string get_run_dir(std::map<std::string,RUNINFO>::const_iterator x) {
            using namespace std;
            char tmp[1024];
            assert(!x->second.run_dir.empty());
            sprintf(tmp,"%s/%s",RUNDIR,x->second.run_dir.data());
            return tmp;
        }

        // inline bool in_a_hurry_to_delete_so_factor_later(rarp::NLIST const & fc_run,int value,char const * bit = "") {

    }

    static std::vector<std::string> nasty_lazy_message;
    // static std::string nasty_lazy_message;

    void get_list(std::vector<std::string> & dirs, char const * const q) {
        Popen ns1(q,16*1024,"r"); 
        char * z3;
        while( *( z3=ns1.getline() ) != '\0') { 
        dirs.push_back(z3);
        }
    }

    inline bool consumed(std::stringstream & silly,std::string const & line,std::string & file) {
        using namespace std;

        if(silly.eof()) {
            // cout << "line consumed!\n";
            return true;
        }else{
            assert( file.substr(file.length()-5,5)=="_Side" || file.substr(file.length()-4,4)=="/CIB" || file.substr(file.length()-7,7)=="/Adjust" );
            stringstream really(line);
            string argh;
            file="";
            for (int g=0; really>>argh ; ++g) { // really) { // while(g++,really>>argh) {
                if(g>=3) {
                    ///// clearly make this a ternary!?!
                    if(file.empty()) file+=argh; 
                    else {
                        file+=' ' + argh;
                    }
                }
            }
            return false;
        }
    }

    // inline void check_space(bool exitwithspace=true) {

    typedef std::pair<long long, std::string> WHATWHEN;


    void daily(int, char **) {

        // if(strcmp(getenv("USER"),"dh2880")!=0) exit(0);

        using namespace std;

        std::map<string,WHATWHEN> info;

        { char host[256]; gethostname(host,sizeof host); assert(strcmp(host,"igm-atav01.igm.cumc.columbia.edu")==0); }

        // bool email = true; // argc==2 && strcmp(argv[1],"email")==0;

        if(opts::email) {
            cout << "will email\n";
            time_t const now_t = time(0);
            char * now = ctime(&now_t);
            if(!now) cout << "unable to get time\n",exit(1);
            cout << "time is " << now;
            cout << "the hour is " << localtime(&now_t)->tm_hour << "\n";
            // assert(localtime(&now_t)->tm_hour==5);
        }

        std::map<std::string,RUNINFO> runs;

        vector<string> dirs;
        get_list(dirs,"find " RUNDIR " -ctime -" DAYS " -name RTA3.cfg -maxdepth 2 | sort | tac ");
        for (int i = 0 ; i<(int)dirs.size(); ++i ) {

            cout << "we have " << dirs[i] << "\n";
            char const * startish, *endish;

            if((startish=strstr(dirs[i].data(),"Runs/"))==0) {
                cout << "this is whated - we email!?!?\n", exit(1);
            }else if((endish=strchr(startish+5,'/'))==0) {
                cout << "this is whated - we email!?!?\n", exit(1);
            }

            char rundir[1024];
            memcpy(rundir,startish+5,endish-(startish+5)+1); // add in '/'?!?
            assert(rundir[endish-(startish+5)]=='/');
            rundir[endish-(startish+5)]='\0';
            cout << "rundir = " << rundir << "\n";

            vector<string> run_info;
            tokenise(run_info,rundir,'_');
            for(int j=0;j<4;++j) {
                cout << "\t"<<j<<"\t"<<run_info[j]<<"\n";
            }
            assert(run_info.size()==4);
            assert(run_info[3].length()==10);
            assert(run_info[3][0]=='A'||run_info[3][0]=='B');

            RUNINFO ri;
            ri.run_started=true; // RTA3.cfg is present...?!?
            ri.machine=run_info[1],ri.run_dir=rundir;
            ri.slot=run_info[3][0];
            int rn = atoi(run_info[2].data());
            assert(rn>0 && rn<100000);
            ri.run_num=rn;

            if(runs.count(run_info[3].substr(1,9))) {
                nasty_lazy_message.push_back(run_info[3].substr(1,9) + " has been seen before - " + dirs[i].substr(0,dirs[i].length()-8) + " is likely to be an aborted run and should be removed.");
                for(int y=0; y<29;++y) cout << "we have seen this flowcell before - ignore likely aborited run!?! - this should send an email asking them to clean-up their mess!?!\n";
                sendmail("nb2975@cumc.columbia.edu","nb2975@cumc.columbia.edu",(dirs[i].substr(0,dirs[i].length()-8) + " appears to be an aborted run.").data(),"Please check and if so remove it.\n");
                continue;
            }

            string cycle_rec = lazy::GetFirstLinePopen("grep NumCycles %s | cut -f2 -d= | perl -pe 's{ }{:}g' ",dirs[i].data());
            cout << "have cycle_rec= '"<<cycle_rec <<"'\n";

            int total_cycles=0;
            {
            vector<string> tmp;
            tokenise(tmp,cycle_rec,':');
            assert(tmp.size()==4); /// would be strange if they go back to single barcoding so block it for now?!?
            for (int t=0;t<4;++t) total_cycles+=atoi(tmp[t].data());
            }

            cout << "let's check a bit more?!? " << dirs[i] << "\n";
            assert(dirs[i].substr(dirs[i].length()-8,8)=="RTA3.cfg");

            ri.rta_complete = isregfile( (dirs[i].substr(0,dirs[i].length()-8)+"RTAComplete.txt").data() ), 
            ri.seq_complete = isregfile( (dirs[i].substr(0,dirs[i].length()-8)+"SequenceComplete.txt").data() ),
            ri.copy_complete = isregfile( (dirs[i].substr(0,dirs[i].length()-8)+"CopyComplete.txt").data() );

            string bored = string("ls -tr ")+dirs[i].substr(0,dirs[i].length()-8)+"Logs/*_Cycle* | tail -n1 | perl -ne 'm{Cycle(\\d+)} && print qq{$1\\n}'";
            // cout << "we get " << bored << "\n";
            bored = lazy::GetFirstLinePopen( bored.data() );
            cout << "we get " << bored << "\n";
            int cycle_cur = atoi(bored.data());
            cout << "cycle_cur = " << cycle_cur << "\n";
            cout << "total_cycles= " << total_cycles <<"\n";

            dirs[i]=dirs[i].substr(0,dirs[i].length()-8)+"RunParameters.xml";

            cout << "let's check a bit more?!? " << dirs[i] << "\n";
            assert(isregfile(dirs[i].data()));
            string fc_type = lazy::GetFirstLinePopen("grep FlowCellMode %s | perl -pe 'm{<FlowCellMode>(.*?)</FlowCellMode>} && print qq{$1\\n}'",dirs[i].data());

            ri.fctype=fc_type;
            ri.total_cycles=total_cycles,ri.last_cycle=cycle_cur;

            cout << "we have fc_type= '" << fc_type << "'\n";
            cout << "we add flowcell\n";

            ri.cycle_rec=cycle_rec;

            runs.insert(make_pair(run_info[3].substr(1,9),ri));

        }

        rarp::NLISTS lims_runs;
        db::get_named_rows("seqdb",lims_runs,"select FCID,FCillumID,rundir,made,Machine,MachineType,ChemVer,concat(LenR1,':',LenI1,':',LenI2,':',LenR2) readlens,fc_status, (LenR1+LenI1+LenI1+LenR2) total_cycles, "
        "ifnull(delete_bcl_time_set,0) delete_bcl_time_set, "
        "ifnull(delete_bcl_time_complete,0) delete_bcl_time_complete, "
        "ifnull(DateRTA,0) DateRTA, "
        "unix_timestamp(ifnull(DateRTA,0)) DateRTA_Epoch, "
        "ifnull(s3_transfer_time_complete,0) s3_transfer_time_complete, "
        "seqtime seqtime_aka_fc_made_time "
        "from Flowcell f where fc_insert>=now()-interval " DAYS " day ");

        assert(lims_runs.size()>0);

        for (int i = 0; i<(int)lims_runs.size(); ++i) {
            for (rarp::NLIST::iterator it=lims_runs[i].begin(); it!=lims_runs[i].end(); it++) {
            }

            runs[lims_runs[i]["FCillumID"]] .lims_status        =lims_runs[i]["fc_status"];
            runs[lims_runs[i]["FCillumID"]] .lims_machine       =lims_runs[i]["Machine"];
            runs[lims_runs[i]["FCillumID"]] .lims_machine_type  =lims_runs[i]["MachineType"];
            runs[lims_runs[i]["FCillumID"]] .lims_chemver       =lims_runs[i]["ChemVer"];
            runs[lims_runs[i]["FCillumID"]] .lims_readlens      =lims_runs[i]["readlens"];
            runs[lims_runs[i]["FCillumID"]] .lims_rundir        =lims_runs[i]["rundir"];
            runs[lims_runs[i]["FCillumID"]] .lims_rta_complete  =lims_runs[i]["DateRTA"];

            runs[lims_runs[i]["FCillumID"]].lims_fcid=atoi(lims_runs[i]["FCID"].data());
            runs[lims_runs[i]["FCillumID"]].lims_made=atoi(lims_runs[i]["made"].data());
            runs[lims_runs[i]["FCillumID"]].lims_total_cycles=atoi(lims_runs[i]["total_cycles"].data());

            runs[lims_runs[i]["FCillumID"]].lims_delete_bcl_time_set=atol(lims_runs[i]["delete_bcl_time_set"].data());
            runs[lims_runs[i]["FCillumID"]].lims_delete_bcl_time_complete=atol(lims_runs[i]["delete_bcl_time_complete"].data());
            runs[lims_runs[i]["FCillumID"]].lims_s3_transfer_time_complete=atol(lims_runs[i]["s3_transfer_time_complete"].data());
            runs[lims_runs[i]["FCillumID"]].lims_fc_made_seqtime=atol(lims_runs[i]["seqtime_aka_fc_made_time"].data());
            runs[lims_runs[i]["FCillumID"]].lims_rta_complete_epoch=atol(lims_runs[i]["DateRTA_Epoch"].data());

        }

        stringstream lazy_bored;

        for (std::map<std::string,RUNINFO>::iterator it=runs.begin(); it!=runs.end(); it++) {

            cout << "CHECKING FLOWCELL " << it->first << "\n";
            if(it->second.run_started) {

                if(it->second.lims_fcid<=0) {
                    cout << "SUMMARY : this flowcell has started but NOT been registered - we send an obnoxious email reminding them not to do this\n";
                    nasty_lazy_message.push_back(it->first + " has started but NOT been registered - we send an obnoxious email reminding them not to do this\n");
                    // sleep(3);
                    //sendmail("nb2975@cumc.columbia.edu","nb2975@cumc.columbia.edu",(it->first + " has not been registered.").data(),"\n");
                    continue;
                }else if(!it->second.lims_made){

                    cout << "SUMMARY : this flowcell wasn't made (need to set seqtime/made for search and machine,machinetype)\n";
                    nasty_lazy_message.push_back(it->first + " wasn't made - we'll make it ourselves - check sequence for exactly what it puts in (made,machine,machinetype,cbot=0,seqtype!?!?!?.");
                    cout << "we do everything 'make flowcell' used to do in one go and get rid of the damned thing - send a message stating we just made it\n";

                    string rundir = todo_move_to_shared_location::get_run_dir(it);
                    // assert(isdir(
                    struct stat get_ctime;
                    stat(rundir.data(),&get_ctime);
                    /////// 'should' we do chemver too?!?
                    string machine =        it->second.machine=="A00116" ? "N1" 
                                :        it->second.machine=="A00123" ? "N2"
                                :        it->second.machine=="A00736" ? "N3"
                                :                                       "";
                    assert(!machine.empty());
                    machine+=it->second.slot;

                    cout << "we check " << rundir << " : " << get_ctime.st_ctime << "\n";
                    cout << "machine " << machine << "\n";
                    char update[2048];
                    sprintf(update,"update Flowcell set seqtime = %lld, cbot = 0, made = 1, machinetype = 'NovaSeq6000', machine = '%s' where fcid = %d and FCillumID = '%s' ; select row_count() updated",
                    (long long)get_ctime.st_ctime,machine.data(),it->second.lims_fcid,it->first.data()); // need to cast to long long
                    cout << "USING " << update << "\n";
                    string up = db::get_named_row("seqdb",update)["updated"];
                    up = up=="0" ? up="Make flowcell error for " + it->first : "Make flowcell " + it->first + " ("+up+")";
                    sendmail("nb2975@cumc.columbia.edu","nb2975@cumc.columbia.edu",up.data(),update);
                    cout << "message " << up << "\n";

                }

                if(it->second.lims_status=="fastq_archived") {

                    cout << "this has finished bcl conversion so could consider wiping!?!?\n";

                    if(it->second.lims_delete_bcl_time_set && it->second.lims_delete_bcl_time_complete==0) {
                        nasty_lazy_message.push_back(it->first + " has been marked by hts for deletion - we check all relevant ism values - should require s3 time set too.");
                    }else if(it->second.lims_delete_bcl_time_set && it->second.lims_delete_bcl_time_complete<=0) {
                        cout << "A SAMPLE IS IN THEORY BEING DELETED\n";
                        nasty_lazy_message.push_back(it->first + " is in theory being deleted.");
                    }else if(!it->second.lims_delete_bcl_time_set && !it->second.lims_delete_bcl_time_complete) {
                        cout << "not even marked for deletion!?!\n";
                    }else{
                        cout << "in theory this was already wiped - clearly this CANNOT actually happen since it's run_started (aka has rta3.cfg) and yet is wiped in limsk\n";
                        assert(0);
                    }
                }else{

                    char done[256];
                    sprintf(done,"%s is in theory running/converting (%.1f%%) on %s",it->first.data(),((float) 100.0 * it->second.last_cycle / it->second.total_cycles),it->second.lims_machine.data());
                    nasty_lazy_message.push_back(done);
                    cout << "this hasn't 't yet completed bcl conversion so don't even consider cleaning up!?!?\n";

                }

                if(it->second.lims_rundir=="NULL") {
                    cout << "--- WE MUST UPDATE THE LIMS RUNDIR!?!? : " << it->second.run_dir.data() << " : " << it->second.lims_fcid << "\n";
                    string count = db::get_named_row("seqdb","update Flowcell set rundir = '%s' where fcid = '%d'; select row_count() updated ",it->second.run_dir.data(),it->second.lims_fcid)["updated"];
                    cout << "we updatd " << count << " entries\n";
                    assert(count=="1");
                }else{

                    cout << "we MUST check that the rundir hasn't changed!?!?!?!!?!? - that should also drive bcl conversion!?!?\n";
                    cout << "\tlims_status= "<<it->second.lims_status<<"\n";
                    cout << "\tlims_fcid= "<<it->second.lims_fcid<<"\n";
                    cout << "\tlims_machine= "<<it->second.lims_machine<<"\n";
                    cout << "\tlims_machine_type= "<<it->second.lims_machine_type<<"\n";
                    cout << "\tlims_chemver= "<<it->second.lims_chemver<<"\n";
                    cout << "\tlims_readlens= "<<it->second.lims_readlens<<"\n";
                    cout << "\tlims_total_cycles= "<<it->second.lims_total_cycles<<"\n";
                    cout << "\tlims_made= "<<it->second.lims_made<<"\n";
                    cout << "\tlims_rundir= "<<it->second.lims_rundir<<"\n";
                    cout << "\tlims_fc_made_seqtime= "<<it->second.lims_fc_made_seqtime<<"\n";
                    cout << "\tlims_s3_transfer_time_complete= "<<it->second.lims_s3_transfer_time_complete<<"\n";
                    cout << "\tlims_delete_bcl_time_set= "<<it->second.lims_delete_bcl_time_set<<"\n";
                    cout << "\tlims_delete_bcl_time_complete= "<<it->second.lims_delete_bcl_time_complete<<"\n";
                    cout << "\trun_dir= "<<it->second.run_dir<<"\n";
                    if(it->second.run_dir!=it->second.lims_rundir) {
                        sendmail("nb2975@cumc.columbia.edu","nb2975@cumc.columbia.edu",(it->first + " has a rundir mismatch!").data(),"Please investigate.\n");
                        nasty_lazy_message.push_back(string("ERROR: ") + it->first + " rundir mismatch!?! previous="+it->second.lims_rundir+", new="+it->second.run_dir);
                    }
                }

                if(it->second.copy_complete) {
                    cout << ">>> this copy complete - i.e. it's completely done\n";
                }else if(it->second.seq_complete) {
                    cout << ">>> this has seq complete\n";
                }else if(it->second.rta_complete) {
                    cout << ">>> this has rta complete - we really need a reminder to be sent if this was more than a few hours agao!?!?\n";
                }else{
                    cout << ">>> this in theory is running " << it->second.last_cycle << "\n";
                }

            }else{
                cout << "SUMMARY : this flowcell is registerd in lims BUT hasn't actually been started - we should probably send a reminder if it was more than a couple days ago\n";
                cout << "for now send me an email for testing purposes!?!\n";
                nasty_lazy_message.push_back(string("WARNING: ") + it->first + " is registerd in lims BUT hasn't actually been started - we should probably send a reminder if it was more than a couple days ago.");
            }

            cout << "we have " << it->first << "\n";
            cout << "\trun_started= "<<it->second.run_started<<"\n";
            cout << "\tmachine= "<<it->second.machine<<"\n";
            cout << "\tstart_date= "<<it->second.start_date<<"\n";
            cout << "\trun_dir= "<<it->second.run_dir<<"\n";
            cout << "\tfctype= "<<it->second.fctype<<"\n";
            cout << "\tcycle_rec= "<<it->second.cycle_rec<<"\n";
            cout << "\trun_num= "<<it->second.run_num<<"\n";
            cout << "\tslot= "<<it->second.slot<<"\n";
            cout << "\tlast_cycle= "<<it->second.last_cycle<<"\n";
            cout << "\ttotal_cycles= "<<it->second.total_cycles<<"\n";
            cout << "\trta_complete= "<<it->second.rta_complete<<"\n";
            cout << "\tseq_complete= "<<it->second.seq_complete<<"\n";
            cout << "\tcopy_complete= "<<it->second.copy_complete<<"\n";
            cout << "\tlims_status= "<<it->second.lims_status<<"\n";
            cout << "\tlims_fcid= "<<it->second.lims_fcid<<"\n";
            cout << "\tlims_machine= "<<it->second.lims_machine<<"\n";
            cout << "\tlims_machine_type= "<<it->second.lims_machine_type<<"\n";
            cout << "\tlims_chemver= "<<it->second.lims_chemver<<"\n";
            cout << "\tlims_readlens= "<<it->second.lims_readlens<<"\n";
            cout << "\tlims_total_cycles= "<<it->second.lims_total_cycles<<"\n";
            cout << "\tlims_made= "<<it->second.lims_made<<"\n";
            cout << "\tlims_rundir= "<<it->second.lims_rundir<<"\n";
            cout << "\tlims_rta_complete= "<<it->second.lims_rta_complete<<"\n";
            cout << "\tlims_rta_complete_epoch= "<<it->second.lims_rta_complete_epoch<<"\n";
            cout << "\tlims_s3_transfer_time_complete= "<<it->second.lims_s3_transfer_time_complete<<"\n";
            cout << "\tlims_delete_bcl_time_set= "<<it->second.lims_delete_bcl_time_set<<"\n";
            cout << "\tlims_delete_bcl_time_complete= "<<it->second.lims_delete_bcl_time_complete<<"\n";
            cout << "\tlims_fc_made_seqtime= "<<it->second.lims_fc_made_seqtime<<"\n";

            float done = ((float) 100.0 * it->second.last_cycle / it->second.total_cycles);
            cout << "\tdone= "<< done << "\n";

            string whatit;
            { rarp::NLISTS hmm;
            /// lexical ordering (i.e. case too) so RG, Sp,
            db::get_named_rows("seqdb",hmm,"select count(1) RGs,replace(e.subproject,' ','_') Sp,e.sample_type type,CurrProjLeader pm from SampleT s join Experiment e on s.sample_id=e.sample_id join prepT p on e.id=p.experiment_id join Lane l on p.prepid=l.prepid where fcid = %d group by e.subproject,e.sample_type order by RGs desc",it->second.lims_fcid);
            for (int x=0;x<(int)hmm.size();++x) {
                for(rarp::NLIST::iterator xit=hmm[x].begin(); xit!=hmm[x].end(); xit++) {
                    cout << "\t"<<xit->first<<"\t:\t"<<xit->second<<"\n";
                    whatit+=xit->second+':';
                }
                // most inefficient thing ever!?!?! - copy not avoid of terminate?!?
                if(whatit[whatit.length()-1]==':') whatit=whatit.substr(0,whatit.length()-1);
                whatit+=';';
            } }
            if(whatit[whatit.length()-1]==';') whatit=whatit.substr(0,whatit.length()-1);
            
            lazy_bored << "[" << it->first << "] " << it->second.lims_chemver << " : " << whatit << "\n";
            char msg[16*1024];
            if(!it->second.rta_complete) {
                sprintf(msg,"Running (%.1f%%) - %s : %s",done,it->second.lims_chemver.data(),it->first.data());
            }else{
                if(it->second.lims_rta_complete.substr(0,4)=="0000") { // this is terrible
                    //// NEED TO PUT IN COPY/RTA TIMESTAMP
                    sprintf(msg,"Run complete - Data converting - %s : %s",it->second.lims_chemver.data(),it->first.data());
                }else{
                    sprintf(msg,"Idle since %s - %s : %s",it->second.lims_rta_complete.substr(0,16).data(),it->second.lims_chemver.data(),it->first.data());
                }
            }

            strcat(msg,( string(" (") + whatit + ")").data());

            // clearly this is a terrible way to do this - should just order the qeury and emit the first of each machine!?!
            std::map<string,WHATWHEN>::iterator tit = info.find(it->second.lims_machine);
            if(tit!=info.end()) {
                if(tit->second.first<it->second.lims_fc_made_seqtime) {
                    tit->second.first=it->second.lims_fc_made_seqtime;
                    tit->second.second=msg;
                }
            }else{
                WHATWHEN x;
                x.first=it->second.lims_fc_made_seqtime; // x.first=pointless_now;
                x.second=msg;
                info.insert(make_pair(it->second.lims_machine,x));
            }

        }

        cout << "\n" << lazy_bored.str() << "\n";
            
        { stringstream body;// ,tmp;

        body << "<style><th, td {padding: 5px;}\nth, td {text-align: left;}\nth ding: 5px;}\n</style>\n"
            << "<table align=\"left\" cellpadding=\"5\" border=\"1\" style=\"border:1px solid black;border-collapse:collapse;width:95%\">\n"
            << "<tr bgcolor=\"#d7d7d7\">"
            << "<th><font color=\"#2996cc\">Slot</font></th>"
            << "<th><font color=\"#2996cc\">State</font></th></tr>\n";

        char const *ml[] = { "N1A", "N1B", "N2A", "N2B", "N3A", "N3B" };
        for (unsigned int i = 0; i<sizeof ml/sizeof ml[0]; ++i) {
            std::map<string,WHATWHEN>::const_iterator it = info.find(ml[i]);

            /* if(it!=info.end()) body << it->first << " : " << it->second.second << "\n";
            else body << ml[i] << " : Idle for longer than check interval (" << DAYS << " days)\n"; */

            body << "<tr bgcolor=\"#d7d7d7\"></td><td>";

            if(it!=info.end()) body << it->first << "</td><td>" << it->second.second << "\n";
            else body << ml[i] << "</td><td> Idle for longer than check interval (" << DAYS << " days)\n";

            body << "</td></tr>";

        }

        body << "</table>";

        cout << body.str();

        if(opts::email) sendmail("jb4393@cumc.columbia.edu,dg2875@cumc.columbia.edu,vsa2105@cumc.columbia.edu,igm-bioinfo@columbia.edu","igm-bioinfo@columbia.edu","NovaSeq occupancy report",body.str().data(),true);
        else sendmail("igm-bioinfo@columbia.edu","igm-bioinfo@columbia.edu","NovaSeq occupancy report",body.str().data(),true); 

        }

        cout << "\nfinal messages " << nasty_lazy_message.size() << "\n";
        for (int i=0; i<(int)nasty_lazy_message.size(); ++i) {
            cout << "["<<i+1<<"] "<<nasty_lazy_message[i] << "\n";
        }

        return;

    }
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

    if(strcmp(1[argv],"load_atav")==0) {
        load(--argc, ++argv);// , myuser); // duh, been a while, ptr decay
        return 0;
        exit(0);
    }else if(strcmp(1[argv],"legacy_release")==0) {
        release_manual(--argc, ++argv);
        return 0;
    }else if(strcmp(1[argv],"bcl")==0) {
        bcl(--argc, ++argv);
        return 0;
    }else if(argc==2) {

        if(strcmp(*(argv+1),"pipe")==0){

            --argc, ++argv;
            metrics(argc,argv);
            archive(argc,argv);
            metrics(argc,argv);
            protect(argc,argv);
            metrics(argc,argv);
            post_pipe::check_pipeline_output(argc,argv);
            metrics(argc,argv);
            post_pipe::cleanup_pipeline_scratch(argc,argv);
            cout << ">we're done\n";

        }else if(strcmp(*(argv+1),"align")==0){

            cout << "Align wrapper.\n";
            while(1) align(argc,argv);

        }else if ( strcmp(1[argv],"daily_summary")==0 ) {

            email_bits::daily(argc,argv);
            return 0;

        }else if ( strcmp(1[argv],"weekly_summary")==0 ) {

            email_bits::lims_driven_summary(argc,argv);
            return 0;

        }else if(strcmp(*(argv+1),"submit")==0) {

            if(!getenv("ATAV_USER")) cerr << "Need ATAV_USER env\n",exit(1);
            if(!getenv("ATAV_PASS")) cerr << "Need ATAV_PASS env\n",exit(1);
            if(!getenv("ATAV_HOST")) cerr << "Need ATAV_HOST env\n",exit(1);

            submit_and_post_checks();

        }else if(strcmp(*(argv+1),"run")==0) {

            opts::commit=true;
            auto_merge();
            auto_release();
            submit_and_post_checks();

        // }else if(strcmp(*(argv+1),"legacy_reset")==0) {
        // exit(0);

        }else cout << "what's '" << *(argv+1) << "\n";

    }else if(argc==3 && strcmp(*(argv+1),"deprecate")==0){

            // deprecate(*(argv+2));
           
    }

    // }else cerr << "unknown mode '" << *(argv+1) << "'\n",exit(1);
    
    cout << "bye\n";

    return 0;


}

#define LOAD_IT "export LUIGI_CONFIG_PATH=/nfs/goldstein/software/dragen_pipe/master/dragen/python/luigi.cfg ; \
  export PATH=/nfs/goldstein/software/git-2.5.0/bin:/nfs/goldstein/software/python2.7.7/bin/:$PATH ; \
  export PYTHONPATH=/nfs/goldstein/software/dragen_pipe/master/dragen/python:/nfs/goldstein/software/dragen_pipe/master/import/src:/nfs/central/home/dh2880/pipeline/gatk_wrappers_utils ; \
  /nfs/goldstein/software/python2.7.7/bin/luigi --module data_import_pipeline ImportSample \
  --sample-id  %s \
  --workers 16 \
  --database waldb_master \
  --dont-remove-tmp-dir-if-failure \
  --local-scheduler --run-locally "
  // --dont-load-data "

/* template <typename A> inline bool query::silly_update(char const * b, A a) {
    char q[2048];
    strcpy(q,b);
    strcat(q,"; select row_count() as affected");
    std::cout << "using " << q << "\n";
    NLIST ars = get_named_row(q,a);
    std::cout << "we modified " << ars["affected"] << " affected rows\n"; 
    return ars["affected"]=="1"; 
    // return ars["affected"]!="0"; 
}*/

// Timey StartTime(time(0));
// cout << "starting at " << StartTime.epoch_time_as_string(tt1) << ", " << StartTime.iso_time(tt2) << "\n";
// memset(timeyt,sizeof timeyt, 0); memcpy(timeyt,z+1,i-1);
// Timey O2(timeyt,"%a %b %d %H:%M:%S %Y",1);   
// int min = (StartTime.epoch_time_as_time_t()-O2.epoch_time_as_time_t())/60;

#define BORED_SQL "select \
  p.sample_type p_sample_type, p.exomeKit p_exomeKit, p.status p_status, FROM_UNIXTIME(p.status_time) p_status_time, p.externalData p_externalData, s.BroadPhenotype s_broadphenotype, \
  dsm.sample_name dsm_sample_name, dsm.sample_type dsm_sample_type, dsm.capture_kit dsm_capture_kit, dsm.pseudo_prepid dsm_pseudo_prepid, \
  dsm.priority dsm_priority, dsm.seqscratch_drive dsm_seqscratch_drive, dsm.is_merged dsm_is_merged, dsm.component_bams dsm_component_bams, \
  w.sample_name w_sample_name, w.sample_finished w_sample_finished, w.sample_failure w_sample_failure, w.prep_id w_pseudo_prepid, w.sample_id w_sampleid, \
  qc.AlignSeqFileLoc qc_archive_dir, qc.pseudo_prepid qc_pseudo_prepid \
  from sequenceDB.dragen_sample_metadata dsm \
  join prepT p on dsm.pseudo_prepid=p.p_prepid \
  left join dsth.sample w on dsm.pseudo_prepid=w.prep_id \
  left join sequenceDB.dragen_qc_metrics qc on dsm.pseudo_prepid=qc.pseudo_prepid \
  join Experiment ex on p.experiment_id=ex.id \
  join SampleT s on ex.sample_id=s.sample_id "


char const * Q1a = "select " // #define Q1 "select "
    " d.pseudo_prepid d_experiment_id , " 
    " SUM(LNYIELD) l_lane_sum , "         
    " p.experiment_id p_experiment_id, "  
    " s.sample_internal_name s_sample_name, "           
    " e.sample_type e_sample_type, "      
    " e.exomekit e_exomekit, "            
    " p.priority p_priority, "            
    " p.externalData p_externalData, "    
    " p.prepid p_prepid, "                
    " p.status p_status, "                
    " p.end_point, "                      
    " d.is_merged d_is_merged, "          
    " w.prep_id w_experiment_id, "                  
    " w.sample_finished w_sample_finished, "        
    " w.sample_failure w_sample_failure, "          
    " e.is_released e_is_released, "                
    " group_concat(concat(p.prepid,':',fcillumid,'.',lanenum)) summary, " // 16
    " dqm.experiment_id dqm_experiment_id, "    // 17
    " dps.experiment_id dps_experiment_id, "    // 18
    " group_concat(pipeline_step_id) steps, "   // 19
    " alignseqfileloc "                         // 20
 " from         Experiment e "
 " join         prepT p                     on e.id                 =p.experiment_id "
 " join         SampleT s                   on s.sample_id          =e.sample_id "
 " left join    Lane l                      on p.prepid             =l.prepid "
 // left join for external samples
 " left join    Flowcell f                  on l.fcid               =f.fcid "
 " left join    dragen_sample_metadata d    on d.pseudo_prepid      =p.experiment_id "
 " left join    dsth.sample w               on w.prep_id            =p.experiment_id "
 " left join    dragen_qc_metrics dqm       on dqm.pseudo_prepid    =p.experiment_id "
 " left join    dragen_pipeline_step dps    on dps.pseudo_prepid    =p.experiment_id"
 " where upper(e.sample_type) = '%s' and (failedprep=0 or failedprep>=100) and p.experiment_id=p.p_prepid ", 
* Q2 = " and e.exomekit = '%s' ";
// #define Q2 " and e.exomekit = '%s' "

int release_manual(int argc, char **argv) { // , opts::MysqlUser & myuser) {
    
    // 'could' have a release-by-flowcell procedure 

    using namespace std;

    { stringstream rarp, parp;
    char where[1024];
    gethostname(where,1024);
    if(
        memcmp(where,"igm-atav",8)!=0
        && memcmp(where,"atav",4)!=0
    ) cout << "Please run me from an atav\n",exit(1);

    rarp<<"Sample release : "<<getenv("USER")<<" @ " << where;
    assert(getcwd(where,1024)); // bored of ignoring return value...
    parp << "From : " << where << "\n";
    for (int o = 0; o < argc; ++o ) {
        parp << "[" << o << "] "<< argv[o] << "\n";
    }
    //sendmail("nb2975@cumc.columbia.edu","nb2975@cumc.columbia.edu",rarp.str().data(),parp.str().data()); 
   }

    string userids;
    // cout << "This wrong " << getenv("USER") << "\n";
    { rarp::NLISTS bored_beyond_words;
    db::get_named_rows("seqdb",bored_beyond_words,"SELECT userid FROM users WHERE NETID= '%s'",getenv("USER"));
    if(bored_beyond_words.size()!=1) cout << "unable to find userid for " << getenv("USER") << "\n",exit(1);
    userids=bored_beyond_words[0]["userid"]; }
    // cout << "using " << userids << "\n";
    // char const * serv = "seqprod";

    static bool commit = false;

    if(argc<4||argc>6) cerr << "usage: " << *argv << " mode=<release/check> file=<list_of_chgvid> sample_type=<sample_type> capture_kit=<capture_type> release_time=<1-14>\n", exit(1);

    // bool commit = false;
    string file, st, capture_kit;
    int release_time = 0;

    char const * of;
    for (int i = 1; i<argc; ++i ){

        if( (of = strchr(*(argv+i),'=')) ==0) cout << "args must be of form what=value\n",exit(1);
        ++of;
        if(memcmp(argv[i],"mode",4)==0) {
            if(strcmp(of,"check")==0) { 
            } else if(strcmp(of,"release")==0) commit=true;
        }
        if(memcmp(argv[i],"file",4)==0) file=of;
        if(memcmp(argv[i],"sample_type",11)==0) st=of;
        if(memcmp(argv[i],"capture_kit",11)==0) capture_kit=of;
        if(memcmp(argv[i],"release_time",12)==0) {
            // assert(*of=='-');
            release_time=atoi(of);
            if(release_time<-14||release_time>0) cout << "release_time offset must be within -14 and 0\n",exit(1);
        }
    }
    if(file.empty()||st.empty()) cout << "You must supply a file and sample_type\n", exit(1);

    if(st!="Genome"&&st!="Exome"&&st!="Custom_Capture"&&st!="RNAseq") 
      cout << "Unrecognised sample_type '"<<st<<"' - valids options are : Genome, Exome, Custom_Capture, RNAseq\n",exit(1);
    cout << "using sample_type '" << st << "'\n";

    if(capture_kit.empty()&&(st!="Genome"&&st!="RNAseq")) cout << "Sample_type " << st << " requires a capture_type argument\n", exit(1); 
    else if(!capture_kit.empty()&&(st=="Genome"||st=="RNAseq")) cout << "Please do not supply a capture_type for " << st << "\n",exit(1);

    cout << "using file '" << file << "'\n"; // cout << "using file '" << *(argv+2) << "'\n";
    if(!isregfile(file.data())) cerr << "There's a problem with your sample list file\n", exit(1);

    if(commit) cout << "WARNING: Will release samples\n";
    else cout << "WARNING: This is a dry run\n";

    std::map<int,FUNKY> list;
    // std::set<int> expt_ids;

    char line[1024];
    FILE *ifs = fopen(file.data(),"r");
    if(!ifs) cerr << "unable to read sample list file\n", exit(1);

    MYSQL *con = mysql_init(NULL);
    if (con == NULL) fprintf(stderr, "mysql_init() failed\n"), exit(1);
    
    if (mysql_real_connect(con, opts::myuser.host(), opts::myuser.user(), opts::myuser.pass(), "sequenceDB", 0, NULL, 0) == NULL)
      finish_with_error(con);
    
    string message;

    while(fgets(line,1024,ifs)){
        int l = strlen(line);
        --l;
        if(line[l]=='\n') line[l]='\0';
        for (int i = 0; i<l; ++i) {
            // if(!::isalpha(line[i])) {
            if( ! ( (line[i]>='a'&&line[i]<='z') || (line[i]>='A'&&line[i]<='Z') || (line[i]>='0'&&line[i]<='9') ) && line[i]!='-' ) 
              cout << "there are stray non-permitted characters in the file '" << line[i] << "'\n", exit(1);
        }

        char ars[2*2048], ars2[2*2048];
        if(st=="Genome") sprintf(ars,Q1a,"Genome");
        else if(st=="RNAseq") sprintf(ars,Q1a,"RNAseq");
        else {
            char tmp[2*2048];
            strcpy(tmp,Q1a);
            strcat(tmp,Q2);
            sprintf(ars,tmp,st.data(),capture_kit.data());
        }

        strcat(ars,"and s.sample_internal_name= '%s'");
        sprintf(ars2,ars,line);
            
        if(getenv("USER") && strcmp(getenv("USER"),"dh2880")==0) {
// cout << "\n"<<ars2<<"\n\n";
        }

        if (mysql_query(con,ars2)) finish_with_error(con);
        
        MYSQL_RES *result = mysql_store_result(con);
        if (result == NULL) finish_with_error(con);

        // int num_fields = mysql_num_fields(result);

        MYSQL_ROW row = mysql_fetch_row(result);
        MYSQL_FIELD *field = mysql_fetch_field(result); 

        fflush(stdout);

        assert(strcmp(field->name,"d_experiment_id")==0);

        cout << "Checking '"<< line <<"' : ";

        ////// null is null ptr!?!?
        if(!row[2]) { // if(result->row_count==0){
            char arsy[2048];
            sprintf(arsy,"expt. %s:%s:%s does not has any eligible (passed/deprecated) preps - skipping\n",line,
              st.data(),(!capture_kit.empty()?capture_kit.data():"NA")
            );
            cout << arsy; message += arsy;
            mysql_free_result(result);
            continue;
        }else if( row[0] || ( row[15] && strcmp(row[15],"released")==0 )  || row[12] || row[17] || row[18] ) {
            string blocks;

            if(row[15] && strcmp(row[15],"released")==0) blocks="expt";
            if(row[0]) { // stop : warning: suggest explicit braces to avoid ambiguous ‘else’
                if(blocks.empty()) blocks="dsm";
                else blocks+=";dsm";
            }
            if(row[18]) {
                if(blocks.empty()) blocks="dps";
                else blocks+=";dps";
            }
            if(row[17]) {
                if(blocks.empty()) blocks="dqm";
                else blocks+=";dqm";
            }
            if(row[12]) {
                if(blocks.empty()) blocks="dragendb";
                else blocks+=";dragendb";
            }

            char arsy[2048];
            sprintf(arsy,"expt. %s:%s:%s (expt_id:%s) has already been released to dragen pipeline - skipping (%s).\n",line,
              st.data(),(!capture_kit.empty()?capture_kit.data():"NA"),row[2], 
                blocks.data()
              );
            cout << arsy; message += arsy;

            mysql_free_result(result);
            continue;
        }else{
            char arsy[2048];
            sprintf(arsy,"expt. %s:%s:%s (expt_id:%s;release_state:%s) can be released (%s).\n",line,
              st.data(),(!capture_kit.empty()?capture_kit.data():"NA"),
              row[2],row[15],row[16]
            );
            cout << arsy; message += arsy;
            int yield=row[1]?atoi(row[1]):0;

            if(st=="Genome" && yield<99000) cout << " > WARNING: WGS Sample has a low yield (" << yield << ")\n";
            else if(st=="Exome" && yield<10000) cout << " > WARNING: WES Sample has a low yield (" << yield << ")\n";
            else if(st=="Custom_Capture" && yield<1000) cout << " > WARNING: Custom Sample has a low yield (" << yield << ")\n";

            if(strcmp(row[9],"Storage")!=0) {
                cout << " > WARNING: Sample doesn't have status 'Storage' - current value is '" << row[9] << "'\n";
            }

            if(row[11]) {
                cout << " > WARNING: Sample has entries in DragenDB (" << row[11] << ")\n";
            }

            field = mysql_fetch_field(result); 

            FUNKY f;
            f.experiment_id=row[2],     f.sample_name=row[3],   f.sample_type=row[4],   f.capture_kit=row[5], 
            f.priority=row[6], f.is_external=(row[7]?"1":"0"), f.prepid=row[8],
            // f.sample_id=row[10], 
            f.end_point=row[10];
              
            list.insert(make_pair(atoi(row[2]),f));

            mysql_free_result(result);

        }
    
    }

    { stringstream x;
    x << "Will release " << list.size() << " samples.\n";
    cout << x.str();
    message+=x.str(); }
    
    if(list.size()==0) {
        /// whatever, can't be bothered to clean up...?!?
        sendmail("nb2975@cumc.columbia.edu","nb2975@cumc.columbia.edu","Release summary",message.data());
        return 0;
    }

    int relc=0;
    for (std::map<int,FUNKY>::iterator it=list.begin(); it!=list.end(); it++) {

        assert(it->first==atoi(it->second.experiment_id.data()));

        { stringstream x;
        x << "Releasing expt_id=" << it->first;
        x << ", sample_name="<< it->second.sample_name;
        x << ", sample_type="<< it->second.sample_type;
        x << ", capture_kit="<< it->second.capture_kit;
        // x << ", sample_id="<< it->second.sample_id;
        x << ", prepid="<< it->second.prepid;
        x << ", priority="<< it->second.priority;
        x << ", is_external= "<< it->second.is_external << "\n";
        cout << x.str();
        message+=x.str(); }

        bool warnings=false;

        if (mysql_query(con, "START TRANSACTION")) finish_with_error(con);    

        char fuk[8*1024];
        sprintf(fuk,"INSERT INTO statusT (STATUS_TIME,    STATUS,     PREPID,     USERID,     POOLID,     SEQID)"
          " SELECT UNIX_TIMESTAMP(),   'Released to Bioinformatics Team', PREPID, %s, 0, 0 FROM prepT WHERE PREPID= %s",
          userids.data(), it->second.prepid.data());

        if(mysql_query(con, fuk)) finish_with_error(con);    
        if(mysql_insert_id(con)==0) warnings=true;
        if(mysql_affected_rows(con)!=1) warnings=true;

        sprintf(
          fuk,
          "UPDATE prepT SET STATUS='Released to Bioinformatics Team', is_released = 1, released_time = ( UNIX_TIMESTAMP() + ( 24*3600*%d ) ) WHERE experiment_id= %s",
          release_time,it->second.experiment_id.data()
        );
        if (mysql_query(con, fuk)) finish_with_error(con);    

        // id = mysql_insert_id(con);
        if(mysql_affected_rows(con)==0) cout << " > WARNING: Sample apparently already had status 'Released to Bioinformatics Team'\n";

        if(st!="RNAseq") {

            sprintf(fuk,"insert into dragen_sample_metadata "
              " (     experiment_id,    pseudo_prepid,  sample_name,    sample_type,    capture_kit,    priority,   is_external,    end_point, is_merged) "
              "values(%s,               %s,             '%s',           '%s',           '%s',           %s,         %s,             '%s', 80000)",
                    it->second.experiment_id.data(), it->second.experiment_id.data(), it->second.sample_name.data(),
                    it->second.sample_type.data(),   it->second.capture_kit.data(),   it->second.priority.data(),
                    it->second.is_external.data(),   it->second.end_point.data() );

            if(mysql_query(con, fuk)) finish_with_error(con);    
            if(mysql_affected_rows(con)!=1) warnings=true;
                  
        }

        sprintf(fuk,"update Experiment set is_released = 'released' where id=%s", it->second.experiment_id.data());
        if(mysql_query(con, fuk)) finish_with_error(con);    
        // duh : if(mysql_insert_id(con)==0) warnings=true;
        if(mysql_affected_rows(con)!=1) warnings=true;

        if(commit) {
            if(warnings) {
                cout << "ERROR: There was problem with the db - not commiting changes.\n";
                if (mysql_query(con, "ROLLBACK")) finish_with_error(con); 
            }else {
                if (mysql_query(con, "COMMIT")) finish_with_error(con);    
                ++relc;
            }
        }else{
            if (mysql_query(con, "ROLLBACK")) finish_with_error(con); 
        }

    }

    { stringstream x;
    x << "Released " << relc << " samples.\n";
    cout << x.str();
    message+=x.str(); }

    sendmail("nb2975@cumc.columbia.edu","nb2975@cumc.columbia.edu","Release summary",message.data());

    mysql_close(con);

    return 0;

}

void load(int argc, char **argv) {

    using namespace std;

    if(1) { // while(1) {

        // new list each time...
        rarp::NLISTS stuff;
        if(argc==1){
            db::get_named_rows("seqdb",stuff,
            BORED_SQL " where (dsm.pseudo_prepid <= 110000 || dsm.pseudo_prepid >= 120000) "
            " and dsm.is_merged = 30 "
            " order by dsm.end_point desc, p.p_prepid desc limit 1");
        }else{
            cout << "grabbing single sample " << argv[1] << "\n";
            db::get_named_rows("seqdb",stuff,BORED_SQL "where dsm.pseudo_prepid = %s group by p.experiment_id",*(argv+1));
        }

        // cout << "WILL PROCESS " << stuff.size() << " ENTRIES\n\n";
        cout.flush();

        if(stuff.size()!=1) return; // break;
        rarp::NLIST & entry = stuff[0];

        if(entry["w_sample_failure"]=="0" && entry["w_sample_finished"]=="1" && entry["p_status"]=="In DragenDB") {
            /* bool modified1 = */ query::silly_update("update dragen_sample_metadata set is_merged = 40 where pseudo_prepid = %s",entry["dsm_pseudo_prepid"].data());
            // cout << "we have 1 modified values = " << modified1 << "\n";
            return;
            // continue;
        }

        cout << "Sample\n";
        for(rarp::NLIST::iterator i=entry.begin(); i!=entry.end(); ++i) { 
            cout << " ["<<i->first<<"]= "<<i->second<<"\n"; 
        }

        assert(isdir(entry["qc_archive_dir"].data()));

        string type = entry["dsm_sample_type"];

        // string const & status = entry["p_status"], & is_merged = entry["dsm_is_merged"];

        for(unsigned i=0; i<entry["dsm_sample_type"].length(); ++i) type[i]=::toupper(type[i]);

        string u_name = entry["dsm_sample_name"]+"."+entry["dsm_pseudo_prepid"],
          dir = "/nfs/"+entry["dsm_seqscratch_drive"]+"/ALIGNMENT/BUILD37/DRAGEN/"+type+"/"+u_name+"/",
          who = dir+"who.txt", 
          sge = dir+"sge_wrapper.log";
          //the param is in waldb.cfg?!? : scratch_area=/nfs/seqscratch_ssd/dh2880/ANNOTATION/DB_LOADING/{sequencing_type}/{sample_name}.{prep_id}/ 

        // does it check for truncated files - i.e. if a chr job exits?!?
        // assert(!isdir(d_dir.data()));

        if( entry["w_sampleid"]=="NULL" && entry["w_sample_name"]=="NULL" && entry["w_sample_failure"]=="NULL"
              && entry["w_sample_finished"]=="NULL" && entry["w_pseudo_prepid"]=="NULL"
        ) {
            // cout << "entirely new...\n";
        }else if( 
              entry["w_sample_finished"]=="0"  && entry["w_sample_name"]==entry["dsm_sample_name"] // && entry["w_sample_failure"]=="0" 
        ) { // entry["w_sampleid"]=="NULL" &&  //&& entry["w_pseudo_prepid"]=="NULL"
            cout << "initialised... - what do we do?!?\n";
        }else if( entry["w_sample_failure"]=="0" && entry["w_sample_finished"]=="1" && entry["w_sample_name"]==entry["dsm_sample_name"]
        ) {
            bool modified1 = query::silly_update("update prepT set status = 'In DragenDB' where p_prepid = %s",entry["dsm_pseudo_prepid"].data());
            // cout << "we have 1 modified values = " << modified1 << "\n";
            assert(modified1==true);
            modified1 = query::silly_update("update dragen_sample_metadata set is_merged = 40 where pseudo_prepid = %s",entry["dsm_pseudo_prepid"].data());
            // cout << "we have 1 modified values = " << modified1 << "\n";
            assert(modified1==true);
            exit(1);
        }else {
            exit(1);
        }

        bool modified1;
        bool patching=false;
        if(entry["dsm_is_merged"]!="41") {
            patching=true;
            modified1 = query::silly_update("update dragen_sample_metadata set is_merged = 41 where pseudo_prepid = %s",entry["dsm_pseudo_prepid"].data());
            // cout << "we have 1 modified values = " << modified1 << "\n";
            assert(modified1==true);
        }

        modified1 = query::silly_update("replace dragen_pipeline_step (pseudo_prepid,pipeline_step_id, version, finish_time, times_ran, step_status) "
          " value ( %s, 32, '0.0.1', now(), 1, 'completed')",entry["dsm_pseudo_prepid"].data());
        // cout << "we have 1 modified values = " << modified1 << "\n";

        string w_sample_id;
        if(entry["w_sampleid"]=="NULL") {

            // cout << "NEW SAMPLE SO POPULATING WALDB\n";
            char YUCK[2048];
            sprintf(YUCK,"mysql -upipeline -pPipeLine001 -h annodb06 WalDB -e \" "
            " insert into sample (sample_name,sample_type,capture_kit,prep_id,priority) "
            " value ('%s','%s','%s',%s,%s ); "
            " select LAST_INSERT_ID()\" ",entry["dsm_sample_name"].data(),entry["dsm_sample_type"].data(),
            entry["dsm_capture_kit"].data(),entry["dsm_pseudo_prepid"].data(),entry["dsm_priority"].data());

            Popen p(YUCK,16*1024,"r");
            p.getline();
            w_sample_id = p.getline();
            assert(!w_sample_id.empty());
            assert(w_sample_id!="");

        }else {
            cout << "PRE-INITIALISED SAMPLE\n";
            cout.flush();
            w_sample_id=entry["w_sampleid"];
        }

        { char cmd[8*1024];
        sprintf(cmd,LOAD_IT,w_sample_id.data());
        // cout << "we will run " << cmd << "\n";
        if(system(cmd)) {

            // cout << "used " << cmd << "\n";
            bool modified1 = query::silly_update("update dragen_sample_metadata set is_merged = 44 where is_merged = 43 and pseudo_prepid = %s",entry["dsm_pseudo_prepid"].data());
            // cout << "we have 1 modified values = " << modified1 << "\n";
            if(!modified1) {
                modified1 = query::silly_update("update dragen_sample_metadata set is_merged = 45 where pseudo_prepid = %s",entry["dsm_pseudo_prepid"].data());
                assert(modified1=true);
            }

        }else{
            // cout<<"yay\n";

            { rarp::NLIST arsv2 = db::get_named_row("seqdb","select status from prepT where p_prepid = %s",entry["dsm_pseudo_prepid"].data()); 
            assert(arsv2["status"]=="In DragenDB"); }

            bool modified1 = query::silly_update("update dragen_sample_metadata set is_merged = 40 where is_merged = 43 and pseudo_prepid = %s",entry["dsm_pseudo_prepid"].data());
            // cout << "we have 1 modified values = " << modified1 << "\n";
            // cout << "used " << cmd << "\n";
            assert(modified1==true);

        }

}

        exit(1);

    }

    cout << "done\n";
}

