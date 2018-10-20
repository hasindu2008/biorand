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


#include "temp.h"
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

// db_t* init_db(core_t* core) {
    // db_t* db = (db_t*)(malloc(sizeof(db_t)));
    // MALLOC_CHK(db);

    // db->capacity_bam_rec = core->opt.batch_size;
    // db->n_bam_rec = 0;

    // db->bam_rec = (bam1_t**)(malloc(sizeof(bam1_t*) * db->capacity_bam_rec));
    // MALLOC_CHK(db->bam_rec);

    // int32_t i = 0;
    // for (i = 0; i < db->capacity_bam_rec; ++i) {
        // db->bam_rec[i] = bam_init1();
        // NULL_CHK(db->bam_rec[i]);
    // }

    // db->fasta_cache = (char**)(malloc(sizeof(char*) * db->capacity_bam_rec));
    // MALLOC_CHK(db->fasta_cache);
    // db->read = (char**)(malloc(sizeof(char*) * db->capacity_bam_rec));
    // MALLOC_CHK(db->read);
    // db->read_len = (int32_t*)(malloc(sizeof(int32_t) * db->capacity_bam_rec));
    // MALLOC_CHK(db->read_len);

    // db->f5 = (fast5_t**)malloc(sizeof(fast5_t*) * db->capacity_bam_rec);
    // MALLOC_CHK(db->f5);

    // db->et = (event_table*)malloc(sizeof(event_table) * db->capacity_bam_rec);
    // MALLOC_CHK(db->et);

    // db->scalings =
        // (scalings_t*)malloc(sizeof(scalings_t) * db->capacity_bam_rec);
    // MALLOC_CHK(db->scalings);

    // db->event_align_pairs =
        // (AlignedPair**)malloc(sizeof(AlignedPair*) * db->capacity_bam_rec);
    // MALLOC_CHK(db->event_align_pairs);
    // db->n_event_align_pairs =
        // (int32_t*)malloc(sizeof(int32_t) * db->capacity_bam_rec);
    // MALLOC_CHK(db->n_event_align_pairs);

    // db->event_alignment = (event_alignment_t**)malloc(
        // sizeof(event_alignment_t*) * db->capacity_bam_rec);
    // MALLOC_CHK(db->event_alignment);
    // db->n_event_alignment =
        // (int32_t*)malloc(sizeof(int32_t*) * db->capacity_bam_rec);
    // MALLOC_CHK(db->n_event_alignment);

    // db->events_per_base =
        // (double*)malloc(sizeof(double*) * db->capacity_bam_rec);
    // MALLOC_CHK(db->events_per_base);

    // db->base_to_event_map =
        // (index_pair_t**)malloc(sizeof(index_pair_t*) * db->capacity_bam_rec);
    // MALLOC_CHK(db->base_to_event_map);

    // db->read_stat_flag = (int32_t *)malloc(sizeof(int32_t) * db->capacity_bam_rec);
    // MALLOC_CHK(db->read_stat_flag);

    // return db;
// }

// int32_t load_db(core_t* core, db_t* db) {
    // // get bams
    // bam1_t* record;
    // int32_t result = 0;
    // db->n_bam_rec = 0;
    // int32_t i = 0;
    // while (db->n_bam_rec < db->capacity_bam_rec) {
        // record = db->bam_rec[db->n_bam_rec];
        // result = sam_itr_next(core->m_bam_fh, core->itr, record);

        // if (result < 0) {
            // break;
        // } else {
            // if ((record->core.flag & BAM_FUNMAP) == 0 &&
                // record->core.qual >=
                    // core->opt
                        // .min_mapq) { // remove secondraies? //need to use the user parameter
                // // printf("%s\t%d\n",bam_get_qname(db->bam_rec[db->n_bam_rec]),result);
                // db->n_bam_rec++;
            // }
        // }
    // }
    // // fprintf(stderr,"%s:: %d queries read\n",__func__,db->n_bam_rec);

    // // get ref sequences (can make efficient by taking the the start and end of
    // // the sorted bam)
    // for (i = 0; i < db->n_bam_rec; i++) {
        // bam1_t* record = db->bam_rec[i];
        // char* ref_name = core->m_hdr->target_name[record->core.tid];
        // // printf("refname : %s\n",ref_name);
        // int32_t ref_start_pos = record->core.pos;
        // int32_t ref_end_pos = bam_endpos(record);
        // assert(ref_end_pos >= ref_start_pos);

        // // Extract the reference sequence for this region
        // int32_t fetched_len = 0;
        // char* refseq =
            // faidx_fetch_seq(core->fai, ref_name, ref_start_pos, ref_end_pos,
                            // &fetched_len); // error handle?
        // db->fasta_cache[i] = refseq;
        // // printf("seq : %s\n",db->fasta_cache[i]);

        // // get the fast5

        // // Get the read type from the fast5 file
        // std::string qname = bam_get_qname(db->bam_rec[i]);
        // char* fast5_path =
            // (char*)malloc(core->readbb->get_signal_path(qname).size() +
                          // 10); // is +10 needed? do errorcheck
        // strcpy(fast5_path, core->readbb->get_signal_path(qname).c_str());

        // //fprintf(stderr,"readname : %s\n",qname.c_str());

        // hid_t hdf5_file = fast5_open(fast5_path);
        // if (hdf5_file >= 0) {
            // db->f5[i] = (fast5_t*)calloc(1, sizeof(fast5_t));
            // MALLOC_CHK(db->f5[i]);
            // fast5_read(hdf5_file, db->f5[i]); // todo : errorhandle
            // fast5_close(hdf5_file);
        // } else {
            // if (core->opt.flag & F5C_SKIP_UNREADABLE) {
                // WARNING("Fast5 file is unreadable and will be skipped: %s",
                        // fast5_path);
            // } else {
                // ERROR("Fast5 file could not be opened: %s", fast5_path);
                // exit(EXIT_FAILURE);
            // }
        // }

        // if (core->opt.flag & F5C_PRINT_RAW) {
            // printf(">%s\tPATH:%s\tLN:%llu\n", qname.c_str(), fast5_path,
                   // db->f5[i]->nsample);
            // uint32_t j = 0;
            // for (j = 0; j < db->f5[i]->nsample; j++) {
                // printf("%d\t", (int)db->f5[i]->rawptr[j]);
            // }
            // printf("\n");
        // }

        // free(fast5_path);

        // //get the read in ascci
        // db->read[i] =
            // (char*)malloc(core->readbb->get_read_sequence(qname).size() +
                          // 1); // is +1 needed? do errorcheck
        // strcpy(db->read[i], core->readbb->get_read_sequence(qname).c_str());
        // db->read_len[i] = strlen(db->read[i]);

        // db->read_stat_flag[i] = 0; //reset the flag
    // }
    // // fprintf(stderr,"%s:: %d fast5 read\n",__func__,db->n_bam_rec);

    // return db->n_bam_rec;
// }

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
    //adapted from ktherad
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







// void output_db(core_t* core, db_t* db) {
    // if (core->opt.flag & F5C_PRINT_EVENTS) {
        // int32_t i = 0;
        // for (i = 0; i < db->n_bam_rec; i++) {
            // printf(">%s\tLN:%d\tEVENTSTART:%d\tEVENTEND:%d\n",
                   // bam_get_qname(db->bam_rec[i]), (int)db->et[i].n,
                   // (int)db->et[i].start, (int)db->et[i].end);
            // uint32_t j = 0;
            // for (j = 0; j < db->et[i].n; j++) {
                // printf("{%d,%f,%f,%f,%d,%d}\t", (int)db->et[i].event[j].start,
                       // db->et[i].event[j].length, db->et[i].event[j].mean,
                       // db->et[i].event[j].stdv, (int)db->et[i].event[j].pos,
                       // (int)db->et[i].event[j].state);
            // }
            // printf("\n");
        // }
    // }
    // if (core->opt.flag & F5C_PRINT_BANDED_ALN) {
        // int32_t i = 0;
        // for (i = 0; i < db->n_bam_rec; i++) {
            // if((db->read_stat_flag[i]) & FAILED_ALIGNMENT){
                // continue;
            // }
            // printf(">%s\tN_ALGN_PAIR:%d\t{ref_os,read_pos}\n",
                   // bam_get_qname(db->bam_rec[i]),
                   // (int)db->n_event_align_pairs[i]);
            // AlignedPair* event_align_pairs = db->event_align_pairs[i];
            // int32_t j = 0;
            // for (j = 0; j < db->n_event_align_pairs[i]; j++) {
                // printf("{%d,%d}\t", event_align_pairs[j].ref_pos,
                       // event_align_pairs[j].read_pos);
            // }
            // printf("\n");
        // }
    // }

    // if (core->opt.flag & F5C_PRINT_SCALING) {
        // int32_t i = 0;
        // printf("read\tshift\tscale\tvar\n");

        // for (i = 0; i < db->n_bam_rec; i++) {
            // if((db->read_stat_flag[i])&(FAILED_ALIGNMENT|FAILED_CALIBRATION)){
                // continue;
            // }
            // printf("%s\t%.2lf\t%.2lf\t%.2lf\n", bam_get_qname(db->bam_rec[i]),
                   // db->scalings[i].shift, db->scalings[i].scale,
                   // db->scalings[i].var);
        // }
    // }
// }

// void free_db_tmp(db_t* db) {
    // int32_t i = 0;
    // for (i = 0; i < db->n_bam_rec; ++i) {
        // bam_destroy1(db->bam_rec[i]);
        // db->bam_rec[i] = bam_init1();
        // free(db->fasta_cache[i]);
        // free(db->read[i]);
        // free(db->f5[i]->rawptr);
        // free(db->f5[i]);
        // free(db->et[i].event);
        // free(db->event_align_pairs[i]);
        // free(db->base_to_event_map[i]);
    // }
// }

// void free_db(db_t* db) {
    // int32_t i = 0;
    // for (i = 0; i < db->capacity_bam_rec; ++i) {
        // bam_destroy1(db->bam_rec[i]);
    // }
    // free(db->bam_rec);
    // free(db->fasta_cache);
    // free(db->read);
    // free(db->read_len);
    // free(db->et);
    // free(db->f5);
    // free(db->scalings);
    // free(db->event_align_pairs);
    // free(db->n_event_align_pairs);
    // free(db->event_alignment);
    // free(db->n_event_alignment);
    // free(db->events_per_base);
    // free(db->base_to_event_map);
    // free(db->read_stat_flag);

    // free(db);
// }

void init_opt(opt_t* opt) {
    memset(opt, 0, sizeof(opt_t));
    //opt->min_mapq = 30;
    //opt->batch_size = 512;
    opt->num_thread = 12;
#ifndef HAVE_CUDA
    opt->flag |= F5C_DISABLE_CUDA;
#endif
    //opt->cuda_block_size=64;    
}
