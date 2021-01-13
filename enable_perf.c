#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

//enables using perf and intel vtune 
//the resultant the binary can be added +s permission to be run from a non-root user; eg, sudo chown root:root enable_perf && sudo chmod +s enable_perf 

int main(){
	sync();
	int fd = open("/proc/sys/kernel/yama/ptrace_scope", O_WRONLY);
	if(fd<0){
		fprintf(stderr,"Opening ptrace scope failed\n");
		exit(1);
	}
	int ret=write(fd, "0", 1);
	if(ret!=1){
		fprintf(stderr,"Setting ptrace scope failed\n");
		exit(1);
	}
	close(fd);


	fd = open("/proc/sys/kernel/perf_event_paranoid", O_WRONLY);
	if(fd<0){
			fprintf(stderr,"Opening perf scope failed\n");
			exit(1);
	}
	ret=write(fd, "-1", 2);
	if(ret!=2){
			fprintf(stderr,"Setting perf scope failed\n");
			exit(1);
	}
    close(fd);

	return 0;
}
