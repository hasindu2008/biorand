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

#include <zlib.h>

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#include "biorand.h"
#include "misc.h"

//#define DEBUG_PRINT
#define STEP_BASES 200
#define QUAL_THRESH_UPPER 7
#define QUAL_THRESH_LOWER 5


#define QUAL_SCORE(seq, i) (((seq)->qual.s[(i)])-33)

static struct option long_options[] = {
    // {"reads", required_argument, 0, 'r'},          //0
    // {"bam", required_argument, 0, 'b'},            //1
    // {"genome", required_argument, 0, 'g'},         //2
    {"threads", required_argument, 0, 't'},        //3
    {"batchsize", required_argument, 0, 'K'},      //4
    // {"print", no_argument, 0, 'p'},                //5
    // {"aaaaaa", no_argument, 0, 0},                 //6
    // {"help", no_argument, 0, 'h'},                 //7
    // {"version", no_argument, 0, 'V'},              //8
    // {"min-mapq", required_argument, 0, 0},         //9
    // {"secondary", required_argument, 0, 0},        //10
    // {"kmer-model", required_argument, 0, 0},       //11
    // {"skip-unreadable", required_argument, 0, 0},  //12
    // {"print-events", required_argument, 0, 0},     //13
    // {"print-banded-aln", required_argument, 0, 0}, //14
    // {"print-scaling", required_argument, 0, 0},    //15
    // {"print-raw", required_argument, 0, 0},        //16
    // {"disable-cuda", required_argument, 0, 0},     //17
    // {"cuda-block-size",required_argument, 0, 0},   //18
    // {"debug-break",required_argument, 0, 0},   //19
    {0, 0, 0, 0}};

void sig_handler(int sig) {
    void* array[100];
    size_t size = backtrace(array, 100);
    ERROR("I regret to inform that a segmentation fault occurred. But at least "
          "it is better than a wrong answer%s",
          ".");
    fprintf(stderr,
            "[%s::DEBUG]\033[1;35m Here is the backtrace in case it is of any "
            "use:\n",
            __func__);
    backtrace_symbols_fd(&array[2], size - 1, STDERR_FILENO);
    fprintf(stderr, "\033[0m\n");
    exit(EXIT_FAILURE);
}

static inline int64_t mm_parse_num(const char* str) //taken from minimap2
{
    double x;
    char* p;
    x = strtod(str, &p);
    if (*p == 'G' || *p == 'g')
        x *= 1e9;
    else if (*p == 'M' || *p == 'm')
        x *= 1e6;
    else if (*p == 'K' || *p == 'k')
        x *= 1e3;
    return (int64_t)(x + .499);
}

static inline void yes_or_no(opt_t* opt, uint64_t flag, int long_idx,
                             const char* arg,
                             int yes_to_set) //taken from minimap2
{
    if (yes_to_set) {
        if (strcmp(arg, "yes") == 0 || strcmp(arg, "y") == 0) {
            opt->flag |= flag;
        } else if (strcmp(arg, "no") == 0 || strcmp(arg, "n") == 0) {
            opt->flag &= ~flag;
        } else {
            WARNING("option '--%s' only accepts 'yes' or 'no'.",
                    long_options[long_idx].name);
        }
    } else {
        if (strcmp(arg, "yes") == 0 || strcmp(arg, "y") == 0) {
            opt->flag &= ~flag;
        } else if (strcmp(arg, "no") == 0 || strcmp(arg, "n") == 0) {
            opt->flag |= flag;
        } else {
            WARNING("option '--%s' only accepts 'yes' or 'no'.",
                    long_options[long_idx].name);
        }
    }
}

void process_single(core_t* core, db_t* db, int32_t optind) {


    gzFile fp;
    kseq_t *seq;
    int l;
    char **argv=core->argv;
    double realtime0=core->realtime0;

    fprintf(stderr,"[%s::%.3f*%.2f] Processing %s (%d of %d)\n",__func__, realtime() - realtime0, cputime() / (realtime() - realtime0),core->argv[optind],optind,core->argc);
 
    fp = gzopen(argv[optind], "r");
    F_CHK(fp,argv[optind]);

    seq = kseq_init(fp);
    MALLOC_CHK(seq);

    int64_t mov_avg_c = 10000;
    float *mov_avg = (float *)malloc(sizeof(float)*mov_avg_c);
    MALLOC_CHK(mov_avg);


    char tmpfilename[1024];
    FILE *file;

    if(core->opt.num_thread>1){

        sprintf(tmpfilename,"%s.filt.fastq",core->argv[optind]);
        file=fopen(tmpfilename,"w");
        F_CHK(file,tmpfilename);
    }
    else{
        file=stdout;
    }


    while ((l = kseq_read(seq)) >= 0) {

        int8_t flag=0;
        assert(l==(int)strlen(seq->seq.s));
        assert(l==(int)strlen(seq->qual.s));
        int32_t i,j;
        double sum=0;

        if(l-STEP_BASES > mov_avg_c){
            mov_avg_c=l;
            mov_avg=(float *)realloc(mov_avg,sizeof(float)*mov_avg_c);
            MALLOC_CHK(mov_avg);
        }    

        for(i=0;i<l-STEP_BASES;i++){
            sum=0;
            for(j=i;j<i+STEP_BASES;j++){
                char quality = QUAL_SCORE(seq,j);
                if(quality<0 || quality>60){
                    //ERROR("Invalid quality score '%c'\n",seq->qual.s[i]);
                    assert(quality>=0 && quality<=60);
                }
                sum += quality;
            }
            mov_avg[i] = sum/STEP_BASES;
        }
        int32_t pos[1000];
        int32_t pos_index=0;

        for(i=0;i<l-STEP_BASES-STEP_BASES-STEP_BASES;i++){
            assert(i+STEP_BASES+STEP_BASES < l-STEP_BASES);
            if(mov_avg[i]>QUAL_THRESH_UPPER && mov_avg[i+STEP_BASES]<QUAL_THRESH_LOWER && mov_avg[i+STEP_BASES+STEP_BASES]>QUAL_THRESH_UPPER){
                flag=1;
                pos[pos_index]=i+STEP_BASES+STEP_BASES/2;
                pos_index++;
                assert(pos_index<1000);            
                i=i+STEP_BASES;
            }
        }

        // assert(sum>=0);
        // if(sum/(float)l>20){
        //     flag=1;
        //     fprintf(stderr,"avg qual %f\n",sum/(float)l);
        // }


        if(flag>0){
            fprintf(file,"@%s", seq->name.s);
            // if (seq->comment.l) {
            //     fprintf(file," %s\n", seq->comment.s);
            // }
            // else{
            //     fprintf(file,"\n");
            // }
            fprintf(file," qual_drop=");
            for(int j=0;j<pos_index;j++){
                fprintf(file,"%d,",pos[j]);
            }
            fprintf(file,"\n");
            fprintf(file,"%s\n", seq->seq.s);
            if (seq->qual.l){
                fprintf(file,"+\n%s\n", seq->qual.s);
            }     
        #ifdef DEBUG_PRINT    
            for(i=0;i<l;i++){
                fprintf(file,"%d ",QUAL_SCORE(seq,i));   
            }
            fprintf(file,"\n");
            fprintf(file,"pos : ");
            for(int j=0;j<pos_index;j++){
                fprintf(file,"%d ",pos[j]);
            }
            fprintf(file,"\n");
            for(i=0;i<l-STEP_BASES;i++){
                fprintf(file,"%f ",mov_avg[i]);
            }
            fprintf(file,"\n");
        #endif
        }

    

    }

    if(core->opt.num_thread>1){
        fclose(file);
    }


    free(mov_avg);

    if(l==-2){
        WARNING("File %s had a truncated quality string after read %s\n",argv[optind],seq->name.s);
    }
    //fprintf("return value: %d\n", l);
    kseq_destroy(seq);
    gzclose(fp);


}


void  filter(int argc, char* argv[]) {

    double realtime0 = realtime();

    const char* optstring = "r:b:g:t:K:hvp";
    int longindex = 0;
    int32_t c = -1;

    optind=2;

    opt_t opt;
    init_opt(&opt);
	
    while ((c = getopt_long(argc, argv, optstring, long_options, &longindex)) >=
           0) {
        // if (c == 'r') {
        //     fastqfile = optarg;
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
        // } 
        if (c == 't') {
            opt.num_thread = atoi(optarg);
            if (opt.num_thread < 1) {
                ERROR("Number of threads should larger than 0. You entered %d",
                      opt.num_thread);
                exit(EXIT_FAILURE);
            }
        }
        // } else if (c == 0 && longindex == 9) {
        //     opt.min_mapq =
        //         atoi(optarg); //check whether this is between 0 and 60
		// }
        // } else if (c == 0 && longindex == 10) { //consider secondary
            // yes_or_no(&opt, F5C_SECONDARY_YES, longindex, optarg, 1);
        // } else if (c == 0 && longindex == 11) {
            // opt.model_file = optarg;
        // } else if (c == 0 && longindex == 12) {
            // yes_or_no(&opt, F5C_SKIP_UNREADABLE, longindex, optarg, 1);
        // } else if (c == 0 && longindex == 13) {
            // yes_or_no(&opt, F5C_PRINT_EVENTS, longindex, optarg, 1);
        // } else if (c == 0 && longindex == 14) {
            // yes_or_no(&opt, F5C_PRINT_BANDED_ALN, longindex, optarg, 1);
        // } else if (c == 0 && longindex == 15) {
            // yes_or_no(&opt, F5C_PRINT_SCALING, longindex, optarg, 1);
        // } else if (c == 0 && longindex == 16) {
            // yes_or_no(&opt, F5C_PRINT_RAW, longindex, optarg, 1);
        // } else if (c == 0 && longindex == 17) {
// #ifdef HAVE_CUDA
            // yes_or_no(&opt, F5C_DISABLE_CUDA, longindex, optarg, 1);
// #else
            // WARNING("%s",
                    // "disable-cuda has no effect when compiled for the CPU");
// #endif
        // } else if(c == 0 && longindex == 18){
            // opt.cuda_block_size = atoi(optarg); //todo : warnining for cpu only mode, check limits
        // }else if(c == 0 && longindex == 19){
            // yes_or_no(&opt, F5C_DEBUG_BRK, longindex, optarg, 1);
        // }

    }

    if (argc-optind < 1) {
        fprintf(
            stderr,
            "Usage: %s %s [OPTIONS] reads.fastq\n",
            argv[0],argv[1]);
        exit(EXIT_FAILURE);
    }
 
    core_t* core = init_core(opt,realtime0);

    core->argc=argc-optind;
    core->argv=&argv[optind];

    if (core->opt.num_thread == 1) {
        int32_t i=0;
        while(argc-optind>0){
            process_single(core,0,i);
            i++;
            optind++;
        }
        

    } 
    else {
        pthread_db(core,0,process_single);    
    }
	

    free_core(core);

    return ;
}


int main(int argc, char* argv[]) {
    
    double realtime0 = realtime();

    signal(SIGSEGV, sig_handler);

    if(argc<2){
        fprintf(stderr,"Usage: %s [PROGRAM] [OPTIONS]\n",argv[0]);
        exit(EXIT_FAILURE);
    }
    else if(strcmp(argv[1],"filter")==0){
        filter(argc, argv);
    }
    else if(strcmp(argv[1],"filterpaf")==0){
        filterpaf(argc, argv);
    }
    else{
        fprintf(stderr,"Unknown program %s\nUsage: %s [PROGRAM] [OPTIONS]\n",argv[1],argv[0]);
        exit(EXIT_FAILURE);
    }




    fprintf(stderr, "[%s] CMD:", __func__);
    for (int i = 0; i < argc; ++i) {
        fprintf(stderr, " %s", argv[i]);
    }


    fprintf(stderr, "\n[%s] Real time: %.3f sec; CPU time: %.3f sec\n",
            __func__, realtime() - realtime0, cputime());
    // }

    // else {
    // fprintf(stderr,"Usage: %s [OPTIONS] -r reads.fa -b alignments.bam -g
    // genome.fa\n",argv[0]);
    // exit(EXIT_FAILURE);
    // }


    return 0;
}


