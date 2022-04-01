#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <pthread.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include <time.h>

//cleans the disk cache by writing and reading back a lot of data
// gcc -Wall clean_fscache2.c -O2 -o clean_fscache -lpthread
 
#define FILE_TMP_PREFIX "/scratch/wv19/hg1112/tmp_file"
#define NUM_THERADS 64
#define GB_PER_THREAD 4

#define ERROR_PREFIX "[%s::ERROR]\033[1;31m "
#define DEBUG2_PREFIX "[%s::DEBUG]\033[1;35m " 

#define NO_COLOUR "\033[0m\n"
#define NEG_CHK(ret) neg_chk(ret, __func__, __FILE__, __LINE__ - 1)

static inline void neg_chk(int ret, const char* func, const char* file,
                           int line) {
    if (ret < 0) {
        fprintf(stderr, ERROR_PREFIX "Unexpected negative value : %s." NO_COLOUR,
                func, strerror(errno));
        fprintf(stderr, DEBUG2_PREFIX "Error occured at %s:%d." NO_COLOUR,
                func, file, line);

        exit(EXIT_FAILURE);
    }
}

#define F_CHK(ret, filename) f_chk((void*) ret, __func__, __FILE__, __LINE__ - 1, filename)

static inline void f_chk(void* ret, const char* func, const char* file,
                         int line, const char* fopen_f) {
    if (ret == NULL) {
        fprintf(stderr, ERROR_PREFIX "Failed to open %s : %s." NO_COLOUR,
                func, fopen_f, strerror(errno));
        fprintf(stderr, DEBUG2_PREFIX "Error occured at %s:%d." NO_COLOUR,
                func, file, line);

        exit(EXIT_FAILURE);
    }
}

/* argument wrapper for the multithreaded framework used for data processing */
typedef struct {
    int64_t write_gbs;
    int32_t thread_index;
	pid_t pid;
	time_t timestamp;
} pthread_arg_t;

void* pthread_single_read(void* voidargs) {
	pthread_arg_t* args = (pthread_arg_t*)voidargs;

	pid_t pid = args->pid;
	time_t timestamp = args->timestamp;
	
	char fname_r[1024];
	sprintf(fname_r, FILE_TMP_PREFIX "_%ld_%d_%d.tmp",timestamp,pid,args->thread_index);
	FILE *file_r = fopen(fname_r,"rb");
	F_CHK(file_r,fname_r);
	
	char fname_w[1024];
	sprintf(fname_w, FILE_TMP_PREFIX "2_%ld_%d_%d.tmp",timestamp,pid,args->thread_index);
	FILE *file_w = fopen(fname_w,"wb");
	F_CHK(file_w,fname_w);

	int64_t bytes =  (args->write_gbs)*1024*1024*1024/4;
	int64_t buff_size = 1024*1024;
	int32_t *buffer  = (int32_t *)malloc(buff_size*4);
	
	int64_t counts = bytes/buff_size;
 
	for(int64_t i=0; i<counts; i++){
		fread(buffer, 4, buff_size,file_r);
		fwrite(buffer, 4, buff_size,file_w);
	}
	
	free(buffer);
	
	fclose(file_r);
	fclose(file_w);
	
	remove(fname_r);
	remove(fname_w);
	
	fprintf(stderr,"%ld GBs done by thread %d\n",(bytes*4)/(1024*1024*1024),args->thread_index);
	
    //fprintf(stderr,"Thread %d done\n",(myargs->position)/THREADS);
    pthread_exit(0);
}

void* pthread_single_write(void* voidargs) {
	pthread_arg_t* args = (pthread_arg_t*)voidargs;
	pid_t pid = args->pid;
	time_t timestamp = args->timestamp;
	
	char fname[1024];
	sprintf(fname,FILE_TMP_PREFIX "_%ld_%d_%d.tmp",timestamp,pid,args->thread_index);
	FILE *file = fopen(fname,"wb");
	F_CHK(file,fname);

	int64_t bytes =  (args->write_gbs)*1024*1024*1024/4;
	int64_t buff_size = 1024*1024;
	int32_t *buffer  = (int32_t *)malloc(buff_size*4);
	
	int64_t counts = bytes/buff_size;
 
	for(int64_t i=0; i<counts; i++){
		for(int j=0;j<buff_size;j++){
				buffer[j]=j;
		}
		fwrite(buffer, 4, buff_size,file);
	}
	
	free(buffer);
	fclose(file);
	fprintf(stderr,"%ld GBs done by thread %d\n",(bytes*4)/(1024*1024*1024),args->thread_index);
	
    //fprintf(stderr,"Thread %d done\n",(myargs->position)/THREADS);
    pthread_exit(0);
}

void pthread_db(int num_thread, int64_t gb_per_thread, int write, pid_t pid, time_t timestamp){
    //create threads
    pthread_t tids[num_thread];
    pthread_arg_t pt_args[num_thread];
    int32_t t, ret;

    //set the data structures
    for (t = 0; t < num_thread; t++) {
        pt_args[t].thread_index = t;
		pt_args[t].write_gbs = gb_per_thread;
		pt_args[t].pid = pid;  
		pt_args[t].timestamp = timestamp;  
    }
	
    //create threads
    for(t = 0; t < num_thread; t++){
		if(write){
			ret = pthread_create(&tids[t], NULL, pthread_single_write,
                                (void*)(&pt_args[t]));
		}
		else{
			ret = pthread_create(&tids[t], NULL, pthread_single_read,
                                (void*)(&pt_args[t]));			
		}
        NEG_CHK(ret);
    }

    //pthread joining
    for (t = 0; t < num_thread; t++) {
        int ret = pthread_join(tids[t], NULL);
        NEG_CHK(ret);
    }
	
	
}

int main(){
	sync();
	
	int num_threads = NUM_THERADS;
	int64_t gb_per_thread = GB_PER_THREAD;
	
	pid_t pid = getpid();
	time_t timestamp = time(NULL);
	
	fprintf(stderr,"Writing tmp files: %s*\n",FILE_TMP_PREFIX);
	pthread_db(num_threads, gb_per_thread,1,pid,timestamp);
	sync();
	
	fprintf(stderr,"reading tmp files: %s*\n",FILE_TMP_PREFIX);
	pthread_db(num_threads, gb_per_thread,0,pid,timestamp);
	sync();	
	
	
	return 0;
}
