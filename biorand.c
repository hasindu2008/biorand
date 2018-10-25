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

core_t* init_core(opt_t opt,double realtime0) {
    core_t* core = (core_t*)malloc(sizeof(core_t));
    MALLOC_CHK(core);
    core->opt = opt;
    //realtime0
    core->realtime0=realtime0;
    return core;
}

void free_core(core_t* core) {
    free(core);
}

#ifdef WORK_STEAL
static inline int32_t steal_work(pthread_arg_t* all_args, int32_t n_threads)
{
	int32_t i, c_i = -1;
	int32_t k;
	for (i = 0; i < n_threads; ++i){
        pthread_arg_t args = all_args[i];
        //fprintf(stderr,"tid : %d, endi : %d, starti : %d\n",i,args.endi,args.starti);
		if (args.endi-args.starti > STEAL_THRESH) {
            //fprintf(stderr,"gap : %d\n",args.endi-args.starti);
            c_i = i;
            break;
        }
    }
    if(c_i<0){
        return -1;
    }
	k = __sync_fetch_and_add(&(all_args[c_i].starti), 1);
    //fprintf(stderr,"k : %d, end %d, start %d\n",k,all_args[c_i].endi,all_args[c_i].starti);
	return k >= all_args[c_i].endi ? -1 : k;
}
#endif

void* pthread_single(void* voidargs) {
    int32_t i;
    pthread_arg_t* args = (pthread_arg_t*)voidargs;
    db_t* db = args->db;
    core_t* core = args->core;

#ifndef WORK_STEAL
    for (i = args->starti; i < args->endi; i++) {
        args->func(core,db,i);
    }
#else
    pthread_arg_t* all_args = (pthread_arg_t*)(args->all_pthread_args);
    //adapted from kthread
	for (;;) {
		i = __sync_fetch_and_add(&args->starti, 1);
		if (i >= args->endi) {
            break;
        }
		args->func(core,db,i);
	}
	while ((i = steal_work(all_args,core->opt.num_thread)) >= 0){
		args->func(core,db,i);  
    }  
#endif

    //fprintf(stderr,"Thread %d done\n",(myargs->position)/THREADS);
    pthread_exit(0);
}



void pthread_db(core_t* core, db_t* db, void (*func)(core_t*,db_t*,int)){
    //create threads
    pthread_t tids[core->opt.num_thread];
    pthread_arg_t pt_args[core->opt.num_thread];
    int32_t t, ret;
    int32_t i = 0;
    int32_t num_thread = core->opt.num_thread;
    int32_t n = core->argc; 
    int32_t step = (n + num_thread - 1) / num_thread;
    
    //todo : check for higher num of threads than the data for efficiency
    //currently handled in a hacky way
    for (t = 0; t < num_thread; t++) {
        pt_args[t].core = core;
        pt_args[t].db = db;
        pt_args[t].starti = i;
        i += step;
        if (i > n) {
            pt_args[t].endi = n;
        } else {
            pt_args[t].endi = i;
        }
        pt_args[t].func=func;
    #ifdef WORK_STEAL    
        pt_args[t].all_pthread_args =  (void *)pt_args;
    #endif
        //fprintf(stderr,"t%d : %d-%d\n",t,pt_args[t].starti,pt_args[t].endi);
    }

    //create threads
    for (t = 0; t < num_thread; t++){
        ret = pthread_create(&tids[t], NULL, pthread_single,
                                (void*)(&pt_args[t]));
        NEG_CHK(ret);
    } 

    //pthread joining
    for (t = 0; t < num_thread ; t++) {
        int ret = pthread_join(tids[t], NULL);
        NEG_CHK(ret);
    }
}


void init_opt(opt_t* opt) {
    memset(opt, 0, sizeof(opt_t));
    opt->num_thread = 12;
}
