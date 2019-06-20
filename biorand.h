#ifndef F5C_H
#define F5C_H

// #include <htslib/faidx.h>
// #include <htslib/hts.h>
// #include <htslib/sam.h>

//option flags
// #define F5C_PRINT_RAW 0x001     
// #define F5C_SECONDARY_YES 0x002 
// #define F5C_SKIP_UNREADABLE   0x004 
// #define F5C_PRINT_EVENTS 0x008
// #define F5C_PRINT_BANDED_ALN 0x010
// #define F5C_PRINT_SCALING 0x020
// #define F5C_DISABLE_CUDA 0x040
// #define F5C_DEBUG_BRK 0x080

//flags for a read
// #define FAILED_CALIBRATION 0x001 //if the calibration failed
// #define FAILED_ALIGNMENT 0x002
// #define FAILED_QUALITY_CHK  0x004


#define WORK_STEAL 1
#define STEAL_THRESH 5

//#define IO_PROC_INTERLEAVE 1
//#define SECTIONAL_BENCHMARK 1   

typedef struct {
    uint32_t flag;
    int32_t batch_size;
    int32_t num_thread;
} opt_t;

typedef struct {
    // region string
    //char* region;

    // bam records
    // bam1_t** bam_rec;
    int32_t capacity_bam_rec; // will these overflow?
    int32_t n_bam_rec;

    // fasta cache //can optimise later by caching a common string for all
    // records in the batch
    char** fasta_cache;

    //read sequence //todo : optimise by grabbing it from bam seq. is it possible due to clipping?
    char** read;
    int32_t* read_len;

} db_t;

typedef struct {
    // options
    opt_t opt;

    //realtime0
    double realtime0;

    int argc;
    char **argv;

} core_t;

typedef struct {
    core_t* core;
    db_t* db;
    int32_t starti;
    int32_t endi;
    void (*func)(core_t*,db_t*,int);
#ifdef WORK_STEAL
    void *all_pthread_args;
#endif
} pthread_arg_t;

typedef struct {
    core_t* core;
    db_t* db;
    //conditional variable for notifying the processing to the output threads
    pthread_cond_t cond;
    pthread_mutex_t mutex;
} pthread_arg2_t;

core_t* init_core(opt_t opt,double realtime0);
void pthread_db(core_t* core, db_t* db, void (*func)(core_t*,db_t*,int));
void free_core(core_t* core);
void init_opt(opt_t* opt);
void filterpaf(int argc, char* argv[]);
void  filterfq(int argc, char* argv[]);
void  comparesam(int argc, char* argv[]);
void olp(int argc, char **argv);

#endif
