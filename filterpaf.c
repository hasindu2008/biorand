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
#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#include <string>
#include <iostream>
using namespace std;

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

//#define ORDERED 1
#define NUM_MAPPINGS 100

int filterpaf(int argc, char* argv[]){
	
    optind=2;

    if (argc-optind < 2) {
        fprintf(
            stderr,
            "Usage: %s %s [OPTIONS] reads.paf reads.fastq\n",
            argv[0],argv[1]);
        exit(EXIT_FAILURE);
    }

    FILE* paffile= fopen(argv[optind],"r");
    F_CHK(paffile,argv[optind]);

    gzFile fp;
    kseq_t *seq;
    int l;
 
    fp = gzopen(argv[optind+1], "r");
    F_CHK(fp,argv[optind]);

    seq = kseq_init(fp);
    MALLOC_CHK(seq);

    //buffers for getline
    size_t bufferSize = 4096;
    char *buffer = (char *)malloc(sizeof(char)*(bufferSize)); 
    MALLOC_CHK(buffer);
    int readlinebytes=1;
    char *pch=NULL;

	int32_t mapped_reads=0;
	int32_t mappings=0;
	int32_t martian_mappings=0;


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
        if(strcmp(pch,"+")){
            curr.strand=0;
        }
        else if(strcmp(pch,"-")){
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

        mappings++;

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
                        if(curr.query_start-prev.query_end>0 && curr.query_start-prev.query_end<5000  && curr.target_start-prev.target_end>0 && curr.target_start-prev.target_end<5000){
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
    
            mapped_reads++;

            #ifndef ORDERED    
                for(int32_t i=0;i<mapping_of_reads_n;i++){
                    for(int32_t j=0;j<mapping_of_reads_n;j++){
                        if(i!=j){
                            alignment_t a = mapping_of_reads[i];
                            alignment_t b = mapping_of_reads[j];
                            if(a.tid==b.tid){
                                //only for the positive strand for the moment
                                if(a.strand==0 && b.strand==0){
                                    if(a.query_start-b.query_end>0 && a.query_start-b.query_end<5000  && a.target_start-b.target_end>0 && a.target_start-b.target_end<5000){
                                            cout << b.rid  << "\t" << b.query_start << "\t" << b.query_end << "\t" << b.tid << "\t" << b.target_start << "\t" << b.target_end << endl;
                                            cout << a.rid  << "\t" << a.query_start << "\t" << a.query_end << "\t" << a.tid << "\t" << a.target_start << "\t" << a.target_end << endl << endl;
                                            assert(a.qual==b.qual);
                                            //cout << curr.qual << endl;;
                                            martian_mappings++;
                                    }    

                                }
                            }
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

    fprintf(stderr,"[%s] Mapped reads %d, Mappings %d, Candidate mappings %d\n",
            __func__, mapped_reads, mappings,martian_mappings);
    kseq_destroy(seq);
    gzclose(fp);
    fclose(paffile);
    free(buffer);
}
