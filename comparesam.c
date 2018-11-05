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

#include <string>
#include <iostream>
using namespace std;

#define NUM_MAPPINGS 1
#define F_UNMAPPED 0x04
#define F_SECONDARY 0x100
#define F_SUPPLEMENTARY 0x800

typedef struct{
    string rid;
    //int32_t query_start;
    //int32_t query_end;
    //int8_t strand;
    string tid;
    int32_t target_start;
    //int32_t target_end;
    int32_t mapq;
    int32_t score;
    int32_t flag;

    //string qual;

}alignment_t;


void comparesam(int argc, char* argv[]){
	
    optind=2;

    if (argc-optind < 2) {
        fprintf(
            stderr,
            "Usage: %s %s [OPTIONS] a.sam b.sam\n",
            argv[0],argv[1]);
        exit(EXIT_FAILURE);
    }

    FILE* samfile_a= fopen(argv[optind],"r");
    F_CHK(samfile_a,argv[optind]);
    FILE* samfile_b= fopen(argv[optind+1],"r");
    F_CHK(samfile_b,argv[optind+1]);


    /****************************sam A**********************************/
    //buffers for getline
    size_t bufferSize_a = 4096;
    char *buffer_a = (char *)malloc(sizeof(char)*(bufferSize_a)); 
    MALLOC_CHK(buffer_a);
    int readlinebytes_a=1;
    char *pch_a=NULL;

    int32_t num_entries_a=0;
	int32_t mapped_reads_a=0;
	int32_t mappings_a=0;

    alignment_t prev_a;
    prev_a.rid = "";

    alignment_t curr_a;
    curr_a.rid = "";

    alignment_t mapping_of_reads_a[NUM_MAPPINGS];
    int32_t mapping_of_reads_n_a=0;

    /****************************sam b********************************/
    //buffers for getline
    size_t bufferSize_b = 4096;
    char *buffer_b = (char *)malloc(sizeof(char)*(bufferSize_b)); 
    MALLOC_CHK(buffer_b);
    int readlinebytes_b=1;
    char *pch_b=NULL;

    int32_t num_entries_b=0;
    int32_t mapped_reads_b=0;
    int32_t mappings_b=0;

    alignment_t prev_b;
    prev_b.rid = "";

    alignment_t curr_b;
    curr_b.rid = "";

    alignment_t mapping_of_reads_b[NUM_MAPPINGS];
    int32_t mapping_of_reads_n_b=0;

	while(1){
        readlinebytes_a=getline(&buffer_a, &bufferSize_a, samfile_a); 
        if(readlinebytes_a == -1){
            break;
        } 

        //ignore header lines
        if(buffer_a[0]=='@' || buffer_a[0]=='\n' || buffer_a[0]=='\r'){
            continue;
        }
        
        //read name
        pch_a = strtok (buffer_a,"\t\r\n"); assert(pch_a!=NULL);
        curr_a.rid = pch_a;

        //flag
        pch_a = strtok (NULL,"\t\r\n"); assert(pch_a!=NULL);
        curr_a.flag = atoi(pch_a);

      
        //RNAME
        pch_a = strtok (NULL,"\t\r\n"); assert(pch_a!=NULL);
        curr_a.tid = pch_a;
        
        //POS
        pch_a = strtok (NULL,"\t\r\n"); assert(pch_a!=NULL);
        curr_a.target_start = atoi(pch_a);
        
        //MAPQ
        pch_a = strtok (NULL,"\t\r\n");  assert(pch_a!=NULL);
        
        //CIGAR
        pch_a = strtok (NULL,"\t\r\n");  assert(pch_a!=NULL);
        
        //RNEXT
        pch_a = strtok (NULL,"\t\r\n");  assert(pch_a!=NULL);
        
        //PNEXT
        pch_a = strtok (NULL,"\t\r\n");  assert(pch_a!=NULL);
        
        //TLEN
        pch_a = strtok (NULL,"\t\r\n");  assert(pch_a!=NULL);
        
        //SEQ
        pch_a = strtok (NULL,"\t\r\n");  assert(pch_a!=NULL);
        
        //QUAL
        pch_a = strtok (NULL,"\t\r\n");  assert(pch_a!=NULL);
        
        if(!(curr_a.flag & 0x04)){
            //NM
            pch_a = strtok (NULL,"\t\r\n");  assert(pch_a!=NULL);
            
            //ms
            pch_a = strtok (NULL,"\t\r\n");  assert(pch_a!=NULL);
            
            //alignment score
            pch_a = strtok (NULL,"\t\r\n");  assert(pch_a!=NULL);
            int32_t as;
            int32_t ret=sscanf(pch_a,"AS:i:%d",&as);
            assert(ret>0);
            mappings_a++;
        }

        num_entries_a++;


        if(curr_a.rid==prev_a.rid ){

            if((curr_a.flag&(F_UNMAPPED|F_SECONDARY|F_SUPPLEMENTARY))==0){           
                mapping_of_reads_a[mapping_of_reads_n_a]=curr_a;
                mapping_of_reads_n_a++;
                if(mapping_of_reads_n_a>=NUM_MAPPINGS){
                    fprintf(stderr,"Too many mapping for read %s\n",curr_a.rid.c_str());
                    mapping_of_reads_n_a--;
                }
            }   


        }

        else{
    
            if(!(curr_a.flag & 0x04)){
                mapped_reads_a++;
            }



            int8_t break_flag = 0;

            while(1){
                    readlinebytes_b=getline(&buffer_b, &bufferSize_b, samfile_b); 
                    if(readlinebytes_b == -1){
                        break_flag=-1;
                        break;
                    } 

                    //ignore header lines
                    if(buffer_b[0]=='@' || buffer_b[0]=='\n' || buffer_b[0]=='\r'){
                        continue;
                    }
                    
                    //read name
                    pch_b = strtok (buffer_b,"\t\r\n"); assert(pch_b!=NULL);
                    curr_b.rid = pch_b;

                    //flag
                    pch_b = strtok (NULL,"\t\r\n"); assert(pch_b!=NULL);
                    curr_b.flag = atoi(pch_b);

                
                    //RNAME
                    pch_b = strtok (NULL,"\t\r\n"); assert(pch_b!=NULL);
                    curr_b.tid = pch_b;
                    
                    //POS
                    pch_b = strtok (NULL,"\t\r\n"); assert(pch_b!=NULL);
                    curr_b.target_start = atoi(pch_b);
                    
                    //MAPQ
                    pch_b = strtok (NULL,"\t\r\n");  assert(pch_b!=NULL);
                    
                    //CIGAR
                    pch_b = strtok (NULL,"\t\r\n");  assert(pch_b!=NULL);
                    
                    //RNEXT
                    pch_b = strtok (NULL,"\t\r\n");  assert(pch_b!=NULL);
                    
                    //PNEXT
                    pch_b = strtok (NULL,"\t\r\n");  assert(pch_b!=NULL);
                    
                    //TLEN
                    pch_b = strtok (NULL,"\t\r\n");  assert(pch_b!=NULL);
                    
                    //SEQ
                    pch_b = strtok (NULL,"\t\r\n");  assert(pch_b!=NULL);
                    
                    //QUAL
                    pch_b = strtok (NULL,"\t\r\n");  assert(pch_b!=NULL);
                    
                    if(!(curr_b.flag & 0x04)){
                        //NM
                        pch_b = strtok (NULL,"\t\r\n");  assert(pch_b!=NULL);
                        
                        //MD
                        pch_b = strtok (NULL,"\t\r\n");  assert(pch_b!=NULL);
                        
                        //alignment score
                        pch_b = strtok (NULL,"\t\r\n");  assert(pch_b!=NULL);
                        int32_t as;
                        int32_t ret=sscanf(pch_b,"AS:i:%d",&as);
                        assert(ret>0);
                        mappings_b++;
                    }

                    num_entries_b++;


                    if(curr_b.rid==prev_b.rid ){

                        if((curr_b.flag&(F_UNMAPPED|F_SECONDARY|F_SUPPLEMENTARY))==0){  

                            mapping_of_reads_b[mapping_of_reads_n_b]=curr_b;
                            mapping_of_reads_n_b++;
                            if(mapping_of_reads_n_b>=NUM_MAPPINGS){
                                fprintf(stderr,"Too many mapping for read %s\n",curr_b.rid.c_str());
                                mapping_of_reads_n_b--;
                            }
                        }


                    }

                    else{
                        prev_b = curr_b;
                        break_flag=1;
                        break;

                    }
                    
                    prev_b = curr_b;

            }

            //compare
            for(int32_t i=0;i<mapping_of_reads_n_a;i++){
                for(int32_t j=0;j<mapping_of_reads_n_b;j++){
                        alignment_t a = mapping_of_reads_a[i];
                        alignment_t b = mapping_of_reads_b[j];
                        if(a.rid!=b.rid){
                            cout << a.rid << '\t' << b.rid <<endl;
                        }
                        assert(a.rid==b.rid);
                        if(a.tid==b.tid){
                           
                        }
                    
                }
            }           


            mapping_of_reads_a[0]=curr_a;
            mapping_of_reads_n_a=1;

            if(break_flag==1){
                 if(!(curr_b.flag & 0x04)){
                     mapped_reads_b++;
                 }
                
                mapping_of_reads_b[0]=curr_b;
                mapping_of_reads_n_b=1;
            }

        }
        
        prev_a = curr_a;
    // fprintf(stderr,"[%s] Num entries in a %d, mappings_a %d, Num_entries in b %d, mappings_b %d\n",
    //         __func__, num_entries_a, mappings_a,num_entries_b, mappings_b);
    }

    fprintf(stderr,"[%s] Num entries in a %d, mappings_a %d, mapped reads in a %d, Num_entries in b %d, mappings_b %d, mapped reads in b %d\n",
            __func__, num_entries_a, mappings_a,mapped_reads_a,num_entries_b, mappings_b,mapped_reads_b);

    fclose(samfile_a);
    fclose(samfile_b);
    free(buffer_a);
    free(buffer_b);
}
