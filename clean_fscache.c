#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h> 
#include <sys/types.h>
#include <sys/stat.h>


//cleans the disk cache of Linux O/S (see https://www.tecmint.com/clear-ram-memory-cache-buffer-and-swap-space-on-linux/)
//the resultant the binary can be added +s permission to be run from a non-root user; eg, sudo chown root:root clean_fscache && sudo chmod +s clean_fscache 

#define DROP_CACHE "/proc/sys/vm/drop_caches"

int main(){
	sync();
	int fd = open(DROP_CACHE, O_WRONLY);
	if(fd<0){
		fprintf(stderr,"Opening " DROP_CACHE " failed : %s\n",strerror(errno));
		exit(EXIT_FAILURE);
	}
	int ret=write(fd, "3", 1);
	if(ret!=1){
		fprintf(stderr,"Writing to " DROP_CACHE " failed : %s\n",strerror(errno));
		exit(EXIT_FAILURE);
	}
	close(fd);
	return 0;
}
