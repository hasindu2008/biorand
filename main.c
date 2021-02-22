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

int main(int argc, char* argv[]) {

    double realtime0 = realtime();

    signal(SIGSEGV, sig_handler);

    if(argc<2){
        fprintf(stderr,"Usage: %s [PROGRAM] [OPTIONS]\n",argv[0]);
        fprintf(stderr,"PROGRAM : \n");
        fprintf(stderr,"\tfilterfq - apply martian filter for a fastq file\n");
        fprintf(stderr,"\tfilterpaf - apply martian filter for a paf file\n");
        fprintf(stderr,"\tcomparesam, - compare two sam files\n");
        fprintf(stderr,"\tolp, - find exact overlaps in reads in a fastq\n");
        fprintf(stderr,"\tidat, - read illumina idat methylation files and prints to stdout in tsv\n");        exit(EXIT_FAILURE);
    }
    else if(strcmp(argv[1],"filterfq")==0){
        filterfq(argc, argv);
    }
    else if(strcmp(argv[1],"filterpaf")==0){
        filterpaf(argc, argv);
    }
    else if(strcmp(argv[1],"comparesam")==0){
        comparesam(argc, argv);
    }
    else if(strcmp(argv[1],"olp")==0){
        olp(argc-1, &argv[1]);
    }
	else if(strcmp(argv[1],"idat")==0){
		idat(argc-1, &argv[1]);
	}
    else{
        fprintf(stderr,"Unknown program %s\nUsage: %s [PROGRAM] [OPTIONS]\n",argv[1],argv[0]);
        fprintf(stderr,"PROGRAM : \n");
        fprintf(stderr,"\tfilterfq - apply martian filter for a fastq file\n");
        fprintf(stderr,"\tfilterpaf - apply martian filter for a paf file\n");
        fprintf(stderr,"\tcomparesam, - compare two sam files\n");
        fprintf(stderr,"\tolp, - find exact overlaps in reads in a fastq\n");
        fprintf(stderr,"\tidat, - read illumina idat methylation files and prints to stdout in tsv\n");
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
