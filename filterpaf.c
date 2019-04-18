#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <execinfo.h>
#include <getopt.h>
#include <unistd.h>
#include <pthread.h>
#include <signal.h>
#include <errno.h> 
#include <sys/resource.h>
#include <sys/time.h>

#include "biorand.h"
#include "misc.h"
#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#include <string>
#include <iostream>
using namespace std;


/* ---------------paf file related parameters------------------------*/
//#define ORDERED 1 //if the mappings of a particular read are in sorted order
#define NUM_MAPPINGS 100 //max number of mappings per read
/* ------------------------------------------------------------------*/


/* ---------------split mapping related parameters--------------------*/
/*
                target gap
                   <--> 
ref  : ----------------------------
           |      |    |      |
read :     --------------------
                   <-->
                query gap 
*/


// profile for generic
//  LOW_THRESH_BASES_QUERY < query gaps < UPPER_THRESH_BASES_QUERY
#define LOW_THRESH_BASES_QUERY (opt->query_gap_min)    //lower threshold for a query gap 
#define UPPER_THRESH_BASES_QUERY (opt->query_gap_max)  //upper threshold for a query gap 
//  LOW_THRESH_BASES_TARGET < target gaps < UPPER_THRESH_BASES_TARGET
#define LOW_THRESH_BASES_TARGET (opt->target_gap_min)     //lower threshold for a target gap 
#define UPPER_THRESH_BASES_TARGET (opt->target_gap_max)  //upper threshold for a target gap 
// |query gap - target gap| < |target gap * GAP_DIFF_RATIO|
#define GAP_DIFF_RATIO (opt->gap_diff_ratio) 

// // profile for martin filters
// //  LOW_THRESH_BASES_QUERY < query gaps < UPPER_THRESH_BASES_QUERY
// #define LOW_THRESH_BASES_QUERY 200    //lower threshold for a query gap 
// #define UPPER_THRESH_BASES_QUERY 5000 //upper threshold for a query gap 
// //  LOW_THRESH_BASES_TARGET < target gaps < UPPER_THRESH_BASES_TARGET
// #define LOW_THRESH_BASES_TARGET 200    //lower threshold for a target gap 
// #define UPPER_THRESH_BASES_TARGET 5000 //upper threshold for a target gap 
// // |query gap - target gap| < |target gap * GAP_DIFF_RATIO|
// #define GAP_DIFF_RATIO 0.1

// // profile for insersions
// //  LOW_THRESH_BASES_QUERY < query gaps < UPPER_THRESH_BASES_QUERY
// #define LOW_THRESH_BASES_QUERY 200    //lower threshold for a query gap 
// #define UPPER_THRESH_BASES_QUERY 5000 //upper threshold for a query gap 
// //  LOW_THRESH_BASES_TARGET < target gaps < UPPER_THRESH_BASES_TARGET
// #define LOW_THRESH_BASES_TARGET -200    //lower threshold for a target gap 
// #define UPPER_THRESH_BASES_TARGET 200 //upper threshold for a target gap 
// // |query gap - target gap| < |target gap * GAP_DIFF_RATIO|
// #define GAP_DIFF_RATIO -1.0

/* ------------------------------------------------------------------*/


/* ------------Fastq phead quality related parameters----------------*/
//#define ABSOLUTE_THRESH 1 // for fastq quality drop, absolute drop or relative drop
#ifdef ABSOLUTE_THRESH
    #define QUAL_THRESH_UPPER 7 //upper threshold for the average phred qual
    #define QUAL_THRESH_LOWER 5
#else
    #define QUAL_THRESH (opt->qual_thresh)//2 //relative average phred score drop
    #define CHECK_QUAL_DROP_DEBUG 1
#endif
#define AVG_WINDOW_SIZE (opt->window_size) //500 //the window size for the average of quality scores outside the gap
/* ------------------------------------------------------------------*/


#define QUAL_SCORE(qual, i) ((qual.c_str()[(i)])-33)


typedef struct{
    string rid;
    int32_t query_start;
    int32_t query_end;
    int8_t strand;
    string tid;
    int32_t target_start;
    int32_t target_end;

    string qual;

}alignment_t;

//stats and global vars
typedef struct{
	
    //stats
    int32_t mapped_reads;
	int32_t mappings_total;
	int32_t split_mappings;
	int32_t split_mappings_same_gap;
	int32_t martian_mappings_with_qual_drop;

    //files 
    FILE* bedfile;

}stat_t;


typedef struct{

    /* ---------------split mapping related parameters--------------------*/
    int32_t query_gap_min;  //lower threshold for a query gap 
    int32_t query_gap_max;  //upper threshold for a query gap 
    int32_t target_gap_min; //lower threshold for a target gap 
    int32_t target_gap_max; //upper threshold for a target gap 
    float gap_diff_ratio;

    /* ------------Fastq phead quality related parameters----------------*/
    int32_t window_size;        //the window size for the average of quality scores outside the gap
    float qual_thresh;       //relative average phred score drop

}filterpaf_opt_t;







static struct option long_options[] = {
    {"profile", required_argument, 0, 'x'},         //0
    {"qmin",required_argument, 0, 0},                //1   
    {"qmax",required_argument, 0, 0},                //2  
    {"tmin",required_argument, 0, 0},                //3
    {"tmax",required_argument, 0, 0},               //4
    {"gap-diff",required_argument, 0, 0},            //5
    {"qual-thresh",required_argument, 0, 0},         //6
    {"w-size",required_argument, 0, 0},              //7
    // {"bam", required_argument, 0, 'b'},            //1
    // {"genome", required_argument, 0, 'g'},         //2
    //{"threads", required_argument, 0, 't'},        //3
    //{"batchsize", required_argument, 0, 'K'},      //4
    // {"print", no_argument, 0, 'p'},                //5
    // {"aaaaaa", no_argument, 0, 0},                 //6
    // {"help", no_argument, 0, 'h'},                 //7
    {0, 0, 0, 0}};


 static inline void print_usage(char **argv, filterpaf_opt_t* opt){
    fprintf(stderr,
    "Usage: %s %s [OPTIONS] reads.paf reads.fastq\n\n"
    "Options :\n"
    "   -x STR              profile (martian or insert) [martian]\n"
    "   --qmin INT          minimum query gap\n"
    "   --qmax INT          maximum query gap\n"
    "   --tmin INT          minimum target gap\n"
    "   --tmax INT          maximum target gap\n"
    "   --gap-diff FLOAT    gap difference ratio (-1.0 to disable)\n"
    "   --qual-thresh FLOAT relative phed quality drop in query\n"
    "   --w-size INT        window size outside the gap for relative phred score\n\n",
    argv[0],argv[1]);

    fprintf(stderr,"Definitions :\n");
    fprintf(stderr,  
                   "                           %d<target_gap<%d         \n"
                   "                                <-->                \n"
                   "      target(ref)  : ----------------------------   \n"
                   "                        |      |    |      |        \n"
                   "      query(read)  :    --------------------        \n"
                   "                                <-->                \n"
                   "                           %d<query_gap<%d          \n\n",
                   LOW_THRESH_BASES_TARGET,UPPER_THRESH_BASES_TARGET,
                   LOW_THRESH_BASES_QUERY,UPPER_THRESH_BASES_QUERY);  

    fprintf(stderr,"gap difference ratio : output only if |query_gap-target_gap|<|target_gap*gap-diff|, ignore if gap-diff<0.0 [current gap-diff %.1f]\n",GAP_DIFF_RATIO);
    fprintf(stderr,"with qual drop (relative drop of  %.1f, window outside gap %d)\n",QUAL_THRESH,AVG_WINDOW_SIZE);
 #ifdef ABSOLUTE_THRESH
    ERROR("%s",Absolute quality drop in effect, qual-thresh is ineffective);
 #endif
 
 }   



static inline void print_qual_score(alignment_t a){

    for (int k=0;k<(int)a.qual.length();k++){
        printf("%d ",QUAL_SCORE(a.qual,k));
    }
    cout << endl;
}

/* check if a split mappings */
static inline int check_if_split_mapping(alignment_t a, alignment_t b, filterpaf_opt_t *opt){
    //positive strand
    if(a.strand==0 && b.strand==0){
        if(a.query_start-b.query_end>LOW_THRESH_BASES_QUERY && a.query_start-b.query_end<UPPER_THRESH_BASES_QUERY  
            && a.target_start-b.target_end>LOW_THRESH_BASES_TARGET && a.target_start-b.target_end<UPPER_THRESH_BASES_TARGET){
            return 1;
        }
        else{
            return 0;
        }
    }
    //negative strand
    else if(a.strand==1 && b.strand==1){
        if(a.query_start-b.query_end>LOW_THRESH_BASES_QUERY && a.query_start-b.query_end<UPPER_THRESH_BASES_QUERY  
            && b.target_start-a.target_end>LOW_THRESH_BASES_TARGET && b.target_start-a.target_end<UPPER_THRESH_BASES_TARGET){
            return 1;    
        }
        else{
            return 0;
        }
    }
    else{
        assert(0);
    }
}

static inline int check_if_similar_gap(alignment_t a, alignment_t b, filterpaf_opt_t* opt){

    if(GAP_DIFF_RATIO<0){ //filter is disabled 
        return 1;
    }
    else{
        //positive stand
        if(a.strand==0 && b.strand==0){
            if( abs((a.query_start-b.query_end) - (a.target_start-b.target_end)) < abs(a.target_start-b.target_end)*GAP_DIFF_RATIO ){
                return 1;
            }
            else{
                return 0;
            }    

        }
        //negative strand
        else if(a.strand==1 && b.strand==1){
            if( abs((a.query_start-b.query_end) - (b.target_start-a.target_end)) < abs(b.target_start-a.target_end)*GAP_DIFF_RATIO ){
                return 1;
            }
            else{
                return 0;
            } 
        }
        else{
            assert(0);
        }
    }
}

/* check if there is a drop in the fastq phred quality scores */
//#define CHECK_QUAL_DROP_DEBUG 1
int32_t check_qual_drop(alignment_t a, alignment_t b,filterpaf_opt_t *opt){


    float sum,avg_prev,avg_split,avg_post;
    
    //check the average qual score for the area before the split
    sum=0;
    int32_t query_end_left = b.query_end-AVG_WINDOW_SIZE+1;

    for(int k=(query_end_left < 0 ? 0 : query_end_left);k<=b.query_end;k++){
        int32_t qual = QUAL_SCORE(a.qual,k);

        if(qual<0){
            cerr << a.qual << endl;
            cerr << k << "\t|" << a.qual.c_str()[(k)] << "|\t" << endl;;
        }
        assert(qual>=0);
        sum+=qual;
    }

    avg_prev = sum/(float)AVG_WINDOW_SIZE;
    if(avg_prev<0){
        cerr << b.rid  << "\t" << b.query_start << "\t" << b.query_end << "\t" << b.tid << "\t" << b.target_start << "\t" << b.target_end << endl;
        cerr << a.rid  << "\t" << a.query_start << "\t" << a.query_end << "\t" << a.tid << "\t" << a.target_start << "\t" << a.target_end << endl;
        cerr << "sum qual " << sum  << "\tavg qual " << avg_prev   << endl;
                                         
    }
    assert(avg_prev>=0);

#ifndef CHECK_QUAL_DROP_DEBUG 
    if(avg_prev<=QUAL_THRESH_LOWER){
        return 0;
    }
#endif

    //avg qual score for the split
    assert(b.query_end<a.query_start); 
    sum=0;    
    for(int k=b.query_end;k<a.query_start;k++){
        sum+=QUAL_SCORE(a.qual,k);
    }
    avg_split = sum/(float)(a.query_start-b.query_end);
    assert(avg_split>=0);
#ifndef CHECK_QUAL_DROP_DEBUG     
    if(avg_split>=QUAL_THRESH_UPPER){
        return 0;
    }
#endif

    //avg qual score after the split    
    sum=0;
    for(int k=a.query_start;k<a.query_start+AVG_WINDOW_SIZE && k<(int)a.qual.length() ;k++){
        sum+=QUAL_SCORE(a.qual,k);
    }              
    avg_post = sum/(float)AVG_WINDOW_SIZE;    
    assert(avg_post>=0); 
#ifndef CHECK_QUAL_DROP_DEBUG   
    if(avg_post<=QUAL_THRESH_LOWER){
        return 0;
    }
#endif

#ifdef CHECK_QUAL_DROP_DEBUG 
    printf("avg_prev %f\tavg_split %f\tavg_post %f\t",avg_prev,avg_split,avg_post);
#endif

#ifdef ABSOLUTE_THRESH
    if(avg_prev>QUAL_THRESH_UPPER && avg_split<QUAL_THRESH_LOWER && avg_post>QUAL_THRESH_UPPER){
#else
    if(avg_prev-avg_split>QUAL_THRESH && avg_post-avg_split>QUAL_THRESH){
#endif
        printf("qualdrop yes\n");
        return 1;    
    }
    else{
        printf("qualdrop no\n");
        return 0;
    }
}

void print_bed_entry(FILE *bedfile, alignment_t a, alignment_t b, filterpaf_opt_t *opt){
    //positive stand
    if(a.strand==0 && b.strand==0){
        fprintf(bedfile,"%s\t%d\t%d\t%s\t%d\t%c\n",a.tid.c_str(),b.target_end-AVG_WINDOW_SIZE,a.target_start+1+AVG_WINDOW_SIZE,
                            a.rid.c_str(),a.target_start-b.target_end, '+');  
    }
    //negative strand
    else if(a.strand==1 && b.strand==1){
         fprintf(bedfile,"%s\t%d\t%d\t%s\t%d\t%c\n",a.tid.c_str(),a.target_end-AVG_WINDOW_SIZE,
                            b.target_start+1+AVG_WINDOW_SIZE,a.rid.c_str(),b.target_start-a.target_end,'-');
                           
    }
    else{
        assert(0);
    }    
}

void print_custom(alignment_t a, alignment_t b, filterpaf_opt_t *opt){
        
    //positive stand
    if(a.strand==0 && b.strand==0){
        cout << b.rid  << "\t" << b.query_start << "\t" << b.query_end << "\t+\t" << b.tid << "\t" << b.target_start << "\t" << b.target_end << endl;
        cout << a.rid  << "\t" << a.query_start << "\t" << a.query_end << "\t+\t" << a.tid << "\t" << a.target_start << "\t" << a.target_end << endl;
        printf("readgap %d\tchrgap %d\n",a.query_start-b.query_end,a.target_start-b.target_end);
        print_qual_score(a);
        cout <<endl;
    }

    //negative strand
    else if(a.strand==1 && b.strand==1){
        cout << b.rid  << "\t" << b.query_start << "\t" << b.query_end << "\t-\t" << b.tid << "\t" << b.target_start << "\t" << b.target_end << endl;
        cout << a.rid  << "\t" << a.query_start << "\t" << a.query_end << "\t-\t" << a.tid << "\t" << a.target_start << "\t" << a.target_end << endl;
        printf("readgap %d\tchrgap %d\n",a.query_start-b.query_end,b.target_start-a.target_end);     
        print_qual_score(a);
        cout <<endl;     
    }
    else{
        assert(0);
    } 
}

void evaluate_mapping_pair(alignment_t a, alignment_t b, stat_t *stats, filterpaf_opt_t *opt){

    FILE *bedfile = stats->bedfile;

    //target ID should be the same
    if(a.tid==b.tid){

        //positive stand
        if(a.strand==0 && b.strand==0){
            
            if(check_if_split_mapping(a, b,opt)){
            //if(a.query_start-b.query_end>LOW_THRESH_BASES && a.query_start-b.query_end<UPPER_THRESH_BASES  && a.target_start-b.target_end>LOW_THRESH_BASES && a.target_start-b.target_end<UPPER_THRESH_BASES){
                    //cout << b.rid  << "\t" << b.query_start << "\t" << b.query_end << "\t+\t" << b.tid << "\t" << b.target_start << "\t" << b.target_end << endl;
                    //cout << a.rid  << "\t" << a.query_start << "\t" << a.query_end << "\t+\t" << a.tid << "\t" << a.target_start << "\t" << a.target_end << endl;
                    assert(a.qual==b.qual);
                    //printf("readgap %d\tchrgap %d\n",a.query_start-b.query_end,a.target_start-b.target_end);
                    print_custom(a,b,opt);

                    //print_qual_score(a);
                    if(LOW_THRESH_BASES_QUERY>=0) assert(a.query_start-b.query_end>=0 && a.target_start-b.target_end>=0);
                    //if( abs(abs(a.query_start-b.query_end) - abs(a.target_start-b.target_end)) < abs(a.target_start-b.target_end)*GAP_DIFF_RATIO ){
                    //if( abs((a.query_start-b.query_end) - (a.target_start-b.target_end)) < abs(a.target_start-b.target_end)*GAP_DIFF_RATIO ){
                    if(check_if_similar_gap(a,b,opt)){    
                        stats->split_mappings_same_gap++;
                        if(check_qual_drop(a,b,opt)){
                            //fprintf(stderr,"%s:%d-%d\t%s:%d-%d\n",a.rid.c_str(),b.query_end,a.query_start,a.tid.c_str(),b.target_end,a.target_start);
                            //fprintf(bedfile,"%s\t%d\t%d\t%s\t%d\t%c\n",a.tid.c_str(),b.target_end-AVG_WINDOW_SIZE,a.target_start+1+AVG_WINDOW_SIZE,
                            //a.rid.c_str(),a.target_start-b.target_end, '+');
                            print_bed_entry(bedfile,a,b,opt);
                            stats->martian_mappings_with_qual_drop++;
                        }                
                    }
                    //cout <<endl;
                    stats->split_mappings++;
            }    

        }

        //negative strand
        if(a.strand==1 && b.strand==1){
            if(check_if_split_mapping(a, b,opt)){
            //if(a.query_start-b.query_end>LOW_THRESH_BASES && a.query_start-b.query_end<UPPER_THRESH_BASES  && b.target_start-a.target_end>LOW_THRESH_BASES && b.target_start-a.target_end<UPPER_THRESH_BASES){
                    //cout << b.rid  << "\t" << b.query_start << "\t" << b.query_end << "\t-\t" << b.tid << "\t" << b.target_start << "\t" << b.target_end << endl;
                    //cout << a.rid  << "\t" << a.query_start << "\t" << a.query_end << "\t-\t" << a.tid << "\t" << a.target_start << "\t" << a.target_end << endl;
                    assert(a.qual==b.qual);
                    //printf("readgap %d\tchrgap %d\n",a.query_start-b.query_end,b.target_start-a.target_end);
                    print_custom(a,b,opt);
                    if(LOW_THRESH_BASES_QUERY>=0) assert(a.query_start-b.query_end>=0 && b.target_start-a.target_end>=0);
                    //print_qual_score(a);
                    //if( abs(abs(a.query_start-b.query_end) - abs(b.target_start-a.target_end)) < abs(b.target_start-a.target_end)*GAP_DIFF_RATIO ){
                    //if( abs((a.query_start-b.query_end) - (b.target_start-a.target_end)) < abs(b.target_start-a.target_end)*GAP_DIFF_RATIO ){
                    if(check_if_similar_gap(a,b,opt)){     
                        stats->split_mappings_same_gap++;
                        if(check_qual_drop(a,b,opt)){
                            //fprintf(bedfile,"%s\t%d\t%d\t%s\t%d\t%c\n",a.tid.c_str(),a.target_end-AVG_WINDOW_SIZE,b.target_start+1+AVG_WINDOW_SIZE,a.rid.c_str(),b.target_start-a.target_end,'-');
                            print_bed_entry(bedfile,a,b,opt);
                            stats->martian_mappings_with_qual_drop++;
                        }  
                    }
                    
                    //cout <<endl;    
                    stats->split_mappings++;
            } 
        }
    }
}


void set_profile(filterpaf_opt_t* opt,const char* profile){

    if(strcmp(profile,"martian")==0){
        opt->query_gap_min=200;
        opt->query_gap_max=5000;
        opt->target_gap_min=200;
        opt->target_gap_max=5000;
        opt->gap_diff_ratio=0.1;
    }
    else if(strcmp(profile,"insert")==0){
        opt->query_gap_min=200;
        opt->query_gap_max=5000;
        opt->target_gap_min=-200;
        opt->target_gap_max=200;
        opt->gap_diff_ratio=-1.0;
    }
    else{
        ERROR("Unkown profile %s",profile);
        exit(EXIT_FAILURE);
    }
}

void init_opt(filterpaf_opt_t* opt) {
    opt->window_size = 500;
    opt->qual_thresh = 2;
    set_profile(opt,"martian");
}

void print_final_stats(filterpaf_opt_t* opt, stat_t stats){

    fprintf(stderr,"[%s] definitions\n",__func__);
    fprintf(stderr,  
                   "                             target_gap             \n"
                   "                                <-->                \n"
                   "      target(ref)  : ----------------------------   \n"
                   "                        |      |    |      |        \n"
                   "      query(read)  :    --------------------        \n"
                   "                                <-->                \n"
                   "                              query_gap             \n");  

    fprintf(stderr," ");

    fprintf(stderr,"[%s] Mapped reads : %d, Total number of mappings : %d\n",__func__, stats.mapped_reads, stats.mappings_total);
    fprintf(stderr,"        Split mappings (query_gap [%d,%d], target_gap [%d-%d]) : %d\n",LOW_THRESH_BASES_QUERY,UPPER_THRESH_BASES_QUERY,LOW_THRESH_BASES_TARGET,UPPER_THRESH_BASES_TARGET,stats.split_mappings);
if(GAP_DIFF_RATIO>=0){    
    fprintf(stderr,"             - similar gap in both query and target ( |query_gap-target_gap| < |target_gap*%.1f| ) : %d\n",GAP_DIFF_RATIO,stats.split_mappings_same_gap);
}
#ifdef ABSOLUTE_THRESH
    fprintf(stderr,"                 - with qual drop (min qual %d, max qual %d, window outside gap %d) : %d\n",
            QUAL_THRESH_LOWER, QUAL_THRESH_UPPER,AVG_WINDOW_SIZE, martian_mappings_with_qual_drop);
#else
    fprintf(stderr,"                 - with qual drop (relative drop of %.1f, window outside gap %d) : %d\n",
                    QUAL_THRESH,AVG_WINDOW_SIZE, stats.martian_mappings_with_qual_drop);
#endif
}

void filterpaf(int argc, char* argv[]){
	

    const char* optstring = "x:";
    int longindex = 0;
    int32_t c = -1;

    optind=2;

    filterpaf_opt_t opt_s;
    filterpaf_opt_t* opt = &opt_s;;
    init_opt(opt);    
	
    while ((c = getopt_long(argc, argv, optstring, long_options, &longindex)) >=0) {
        if (c == 'x') {
            set_profile(opt,optarg);
        }
        else if (c == 0 && longindex == 1){
            opt->query_gap_min = atoi(optarg);
        }
        else if (c == 0 && longindex == 2){
            opt->query_gap_max = atoi(optarg);
        }
        else if (c == 0 && longindex == 3){
            opt->target_gap_min = atoi(optarg);
        }            
        else if (c == 0 && longindex == 4){
            opt->target_gap_min = atoi(optarg);
        }
        else if (c == 0 && longindex == 5){
            opt->gap_diff_ratio = atof(optarg);
        }
        else if (c == 0 && longindex == 6){
            opt->qual_thresh = atof(optarg);
        }
        else if (c == 0 && longindex == 7){
            opt->window_size = atoi(optarg);
        }                        
        // } else if (c == 'b') {
        //     bamfilename = optarg;
        // } else if (c == 'g') {
        //     fastafile = optarg;
        // // } else if (c == 'p') {
        //     // opt.flag |= F5C_PRINT_RAW;
        // } else if (c == 'K') {
        //     opt.batch_size = atoi(optarg);
        //     if (opt.batch_size < 1) {
        //         ERROR("Batch size should larger than 0. You entered %d",
        //               opt.batch_size);
        //         exit(EXIT_FAILURE);
        //     }

    }

    if (argc-optind < 2) {
        print_usage(argv,opt);
        exit(EXIT_FAILURE);
    }

    FILE* paffile= fopen(argv[optind],"r");
    F_CHK(paffile,argv[optind]);

    gzFile fp;
    kseq_t *seq;
    int l;
 
    fp = gzopen(argv[optind+1], "r");
    F_CHK(fp,argv[optind+1]);

    seq = kseq_init(fp);
    MALLOC_CHK(seq);

    //buffers for getline
    size_t bufferSize = 4096;
    char *buffer = (char *)malloc(sizeof(char)*(bufferSize)); 
    MALLOC_CHK(buffer);
    int readlinebytes=1;
    char *pch=NULL;


    //stats
    stat_t stats;
    memset(&stats,0, sizeof(stat_t));


    stats.bedfile = fopen("candidates.bed","w");
    F_CHK(stats.bedfile,"candidates.bed");

    alignment_t prev;
    prev.rid = "";

    alignment_t curr;
    curr.rid = "";

#ifndef ORDERED
    alignment_t mapping_of_reads[NUM_MAPPINGS];
    int32_t mapping_of_reads_n=0;
#endif

	while(1){
        readlinebytes=getline(&buffer, &bufferSize, paffile); 
        if(readlinebytes == -1){
            break;
        } 

        
        //read name
        pch = strtok (buffer,"\t\r\n"); assert(pch!=NULL);
        curr.rid = pch;

        //readlen
        pch = strtok (NULL,"\t\r\n"); assert(pch!=NULL);

        //query start
        pch = strtok (NULL,"\t\r\n"); assert(pch!=NULL);
        curr.query_start = atoi(pch);

        //query end
        pch = strtok (NULL,"\t\r\n"); assert(pch!=NULL);
        curr.query_end= atoi(pch);

        //relative strand
        pch = strtok (NULL,"\t\r\n"); assert(pch!=NULL);
        if(strcmp(pch,"+")==0){
            curr.strand=0;
        }
        else if(strcmp(pch,"-")==0){
            curr.strand=1;
        }
        else{
            assert(0);
        }


        //targetname
        pch = strtok (NULL,"\t\r\n"); assert(pch!=NULL);
        curr.tid = pch;

        //target len
        pch = strtok (NULL,"\t\r\n"); assert(pch!=NULL);

        //target start
        pch = strtok (NULL,"\t\r\n"); assert(pch!=NULL);
        curr.target_start = atoi(pch);

        //target end
        pch = strtok (NULL,"\t\r\n"); assert(pch!=NULL);
        curr.target_end= atoi(pch);

        stats.mappings_total++;

        if(curr.rid==prev.rid ){

            #ifndef ORDERED    
                mapping_of_reads[mapping_of_reads_n]=curr;
                mapping_of_reads_n++;
                if(mapping_of_reads_n>=NUM_MAPPINGS){
                    fprintf(stderr,"Too many mapping for read %s\n",curr.rid.c_str());
                    mapping_of_reads_n--;
                }

            #else
                if(curr.tid==prev.tid){
                    //only for the positive strand for the moment
                    if(curr.strand==0 && prev.strand==0){
                        if(curr.query_start-prev.query_end>LOW_THRESH_BASES && curr.query_start-prev.query_end<UPPER_THRESH_BASES  && curr.target_start-prev.target_end>LOW_THRESH_BASES && curr.target_start-prev.target_end<UPPER_THRESH_BASES){
                                cout << prev.rid  << "\t" << prev.query_start << "\t" << prev.query_end << "\t" << prev.tid << "\t" << prev.target_start << "\t" << prev.target_end << endl;
                                cout << curr.rid  << "\t" << curr.query_start << "\t" << curr.query_end << "\t" << curr.tid << "\t" << curr.target_start << "\t" << curr.target_end << endl << endl;
                                assert(curr.qual==prev.qual);
                                //cout << curr.qual << endl;;
                        }    

                    }
                }
            #endif

        }

        else{
    
            stats.mapped_reads++;

            #ifndef ORDERED    
                for(int32_t i=0;i<mapping_of_reads_n;i++){
                    for(int32_t j=0;j<mapping_of_reads_n;j++){
                        if(i!=j){
                            alignment_t a = mapping_of_reads[i];
                            alignment_t b = mapping_of_reads[j];
                            evaluate_mapping_pair(a,b,&stats,opt);
                        }
                    }
                }
            #endif

            while (1) {
                l = kseq_read(seq);
                assert(l>=0);
                //cout << "ffff\t" << seq->name.s << "\t" << curr.rid.c_str() << endl;
                if(strcmp(seq->name.s,curr.rid.c_str())==0){
                    //cout << seq->name.s << "\t" ;
                    //cout << "hhhhh\n";
                    curr.qual = seq->qual.s;
                    break;
                }
            }

        #ifndef ORDERED
            mapping_of_reads[0]=curr;
            mapping_of_reads_n=1;
        #endif

        }
        
        prev = curr;

    }

    print_final_stats(opt,stats);

    kseq_destroy(seq);
    gzclose(fp);
    fclose(paffile);
    fclose(stats.bedfile);
    free(buffer);
}
