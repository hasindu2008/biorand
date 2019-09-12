#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <sys/time.h>
#include <assert.h>

#define MAX_READ_SIZE 125       //maximum possible read size in number of bases
#define MAX_READ_ID_SIZE 50     //maximum length of read ID
#define BIT_MWORD 64            //number of bits in a machine word
                
//#define DEBUG 1
#define LINEAR_SORT 1

typedef int seqid_t ;           //change this for large number of reads that will overflow int
typedef unsigned long long int MWORD;

//set based on commandline args
int READ_SIZE= 15;
int NUM_READS= 10;
int OVERLAP_MAX= 10;
int OVERLAP_MIN= 5;

//computed based on runtime info
int SR_WORD; //number of words that contain the converted bitvector  ((READ_SIZE+BIT_MWORD/2-1)/(BIT_MWORD/2)) i.e. ceil(READ_SIZE/(BIT_MWORD/2))
int SR_CONV_WORD; //number of words that contain the input ascii sequence ((READ_SIZE+7)/8) i.e. ceil(READ_SIZE/8)

//global counters
int overlap = 5; //the considered overlap at the meoment
long pairs = 0;  //counter to keep track of number of pairs

/*Die on error. Print the error and exit if the return value of the previous function NULL*/
#define errorCheckNULL(ret) ({\
    if (ret==NULL){ \
        fprintf(stderr,"Error at File %s line number %d : %s\n",__FILE__, __LINE__,strerror(errno));\
        exit(EXIT_FAILURE);\
    }\
    })

/*Die on error. Print the error and exit if the return value of the previous function is -1*/
#define errorCheck(ret) ({\
    if (ret<0){ \
        fprintf(stderr,"Error at File %s line number %d : %s\n",__FILE__, __LINE__,strerror(errno));\
        exit(EXIT_FAILURE);\
    }\
    })

//easy timers
#define INITIALISE_TIMERS   struct timeval start_time, stop_time; double computed_time; 
#define START_TIMER         gettimeofday(&start_time, NULL);
#define STOP_TIMER   	    gettimeofday(&stop_time, NULL); tv_sub(&stop_time, &start_time); \
                            computed_time=((stop_time.tv_sec)*1000.0 + (stop_time.tv_usec)/1000.0)/1000.0;
#define PRINT_TIME(...)     fprintf(stderr, __VA_ARGS__); fprintf(stderr, "(Time : %lf s)\n", computed_time);   

/**Function to substract two times taken using gettimeofday()*/
void tv_sub(struct  timeval *out, struct timeval *in){
	if ((out->tv_usec -= in->tv_usec) <0)
	{
		--out ->tv_sec;
		out ->tv_usec += 1000000;
	}
	out->tv_sec -= in->tv_sec;
}    
    
//structure for storing a read    
struct read{
    char read_id[MAX_READ_ID_SIZE+1] ; // this should be optimised
    MWORD sequence[(MAX_READ_SIZE+BIT_MWORD/2-1)/(BIT_MWORD/2)]; //READ in bitvector
    seqid_t seqid; //to do : concat this with reverse complement
    char reversecomplement; //0 for foward 1 for reverse
};    
    
    
void print_to_paf(FILE *output, struct read *targetReads, int location_reads, struct read *queryReads, int location_reads_left, int theoverlap){
 
// 1	string	Query sequence name                         : targetReads[location_reads].read_id
// 2	int	Query sequence length                           : READ_SIZE   
// 3	int	Query start (0-based)                           : query_start
// 4	int	Query end (0-based)                             : query_end
// 5	char	Relative strand: "+" or "-"                 : sign
// 6	string	Target sequence name                        : queryReads[location_reads_left].read_id
// 7	int	Target sequence length                          : READ_SIZE
// 8	int	Target start on original strand (0-based)       : target_start
// 9	int	Target end on original strand (0-based)         : target_end
// 10	int	Number of residue matches                       : overlap
// 11	int	Alignment block length                          : overlap
// 12	int	Mapping quality (0-255; 255 for missing)        : 255

    if(queryReads[location_reads_left].seqid < targetReads[location_reads].seqid){
    
        char sign;
        int query_start;
        int query_end;
        int target_start;
        int target_end;
        
        //is query or target reversed
        char target_reversed = targetReads[location_reads].reversecomplement;
        char query_reversed = queryReads[location_reads_left].reversecomplement;
        
        //both are forward
        if(query_reversed==0 && target_reversed==0){
            sign ='+';
            query_start = READ_SIZE-theoverlap;
            query_end = READ_SIZE-1;
            target_start = 0;
            target_end = overlap-1;
        }
        //both are reversed
        else if(query_reversed==1 && target_reversed==1){
            sign ='+';
            query_start = 0;
            query_end = theoverlap-1 ;       
            target_start = READ_SIZE-theoverlap;
            target_end = READ_SIZE-1;
        }
        //query forward other reverse
        else if(query_reversed==0 && target_reversed==1){
            sign ='-';
            query_start = READ_SIZE-theoverlap;
            query_end = READ_SIZE-1;
            target_start = READ_SIZE-theoverlap;
            target_end = READ_SIZE-1;      
            
        }
        //query reverse other forward
        else if(query_reversed==1 && target_reversed==0){
            sign ='-';
            query_start = 0;
            query_end = theoverlap-1 ;       
            target_start = 0;
            target_end = theoverlap-1 ;      
        }
        else{
            fprintf(stderr,"Thi sis impossible!\n");
            exit(EXIT_FAILURE);
        }
        
        
        fprintf(output,"%s\t%d\t%d\t%d\t%c\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n",queryReads[location_reads_left].read_id, READ_SIZE, query_start, query_end, sign, targetReads[location_reads].read_id, READ_SIZE, target_start, target_end,overlap,overlap,255 );
        pairs++;
        
    }
}    
  
//get the complement of a base  
char complement(char in){
    char out;
    switch(in){
        case 'A':
            out='T';
            break;
        case 'C':
            out='G';
            break;
        case 'G':
            out='C';
            break;
        case 'T':
            out='A';
            break;
        case 'N':
            fprintf(stderr,"A base of value N found\n");
            out='N';
            break;
        default:
            fprintf(stderr,"Some invalid base value %c\n",in);  
            exit(EXIT_FAILURE);    
    }
    return out;
        
}
 
//compute reverse complement of a string 
void reverse_complement(char *destination,char *source){
    int i;
    for(i=0;i<READ_SIZE;i++){
        destination[i] = complement(source[READ_SIZE-1-i]);
    }
}    
    
    
/********************************************** Bit vector stuff ************************************/

void Base2Bit(MWORD *bit_vector, MWORD *ascii_vector) // Bit level parallelized convert method
{
    int i;
    // convert 8 base in each word first
    for (i = 0; i < SR_CONV_WORD; i++){
        ascii_vector[i] &= 0x0606060606060606llu;
        ascii_vector[i] >>= 1;
        ascii_vector[i] |= (ascii_vector[i] << 10) | (ascii_vector[i] << 20) | (ascii_vector[i] << 30);
        ascii_vector[i] &= 0xFF000000FF000000;
        ascii_vector[i] = (ascii_vector[i] >> 8) | (ascii_vector[i] << 32);
        ascii_vector[i] &= 0xFFFF000000000000;
    }

    //merge every 4 words together
    
    for (i = 0; i < SR_WORD; i++){
        int idx = (i * 4);
        bit_vector[i] = ascii_vector[idx] | ascii_vector[idx + 1] >> 16 | ascii_vector[idx + 2] >> 32 | ascii_vector[idx + 3] >> 48;
    }
}

//print bitvectors which are left aligned
void Print_Bitvector(char *seq_ascii, MWORD *bit_vector){
    
    int printed = 0;
    int i,j;
    
    for (i = 0; i < SR_WORD; i++){
        MWORD temp = bit_vector[i];
        for (j = 0; j < BIT_MWORD/2; j++){
            if (printed == READ_SIZE){
                break;
            }
            switch (temp >> (BIT_MWORD - 2)){
                case 0: seq_ascii[printed]='A'; break;
                case 1: seq_ascii[printed]='C'; break;
                case 2: seq_ascii[printed]='T'; break;
                case 3: seq_ascii[printed]='G'; break;
            }
        temp <<= 2;
        printed++;
        }
    }
    seq_ascii[printed] = 0;

}

//print bitvectors which are right aligned
void Print_Bitvector_right_shifted(char *seq_ascii, MWORD *bit_vector){
    
    int printed = 0;
    int i,j;
    
    for (i = 0; i < SR_WORD; i++){
        MWORD temp = bit_vector[i];
        for (j = 0; j < BIT_MWORD/2; j++){
            if (printed == READ_SIZE){ 
                break;
            }
            switch (temp >> (BIT_MWORD - 2)){
                case 0: seq_ascii[printed]='A'; break;
                case 1: seq_ascii[printed]='C'; break;
                case 2: seq_ascii[printed]='T'; break;
                case 3: seq_ascii[printed]='G'; break;
            }
            temp <<= 2;
            printed++;
            if(i==0 && j < SR_WORD*BIT_MWORD/2-READ_SIZE){
                printed--;
            }
        }
    }
    seq_ascii[printed] = 0;
}


void Shift_Right(int shift, MWORD *bit_vector) // less than the size of a word
{
    int i;
    int shift_inverse = BIT_MWORD - shift;
    for (i = SR_WORD - 1; i > 0; i--)
    {
        bit_vector[i] >>= shift;
        bit_vector[i] |= bit_vector[i - 1] << shift_inverse;
        }
        bit_vector[0] >>= shift;
}
    
    
/***********************************************Sorters***********************************************/    
    
//compare function for sorting by right corner   
/*int sortbyLastBases(const void *a,const void *b){
    struct read *readA = (struct read *)a;
    struct read *readB = (struct read *)b;
    return strcmp(&(readA->sequence[READ_SIZE-overlap]),&(readB->sequence[READ_SIZE-overlap]));
    
}   */   
    
//compare function for sorting by right corner   
int sortbyLastBases(const void *a,const void *b){
    struct read *readA = (struct read *)a;
    struct read *readB = (struct read *)b;
    MWORD seqA;
    MWORD seqB;
    
    // this commented code is supports for a maximum of 32 bases overlap
    /*if(overlap<=BIT_MWORD/2){
        MWORD seqA = readA->sequence[SR_WORD-1];
        seqA = seqA << (BIT_MWORD-overlap*2);
        MWORD seqB = readB->sequence[SR_WORD-1];
        seqB = seqB << (BIT_MWORD-overlap*2);
        return ((seqA>seqB) - (seqA<seqB));
        
    }
    else{
        fprintf(stderr,"Please implement multiple machine words\n");
        exit(EXIT_FAILURE);
    }*/

    
    int N = ((overlap+BIT_MWORD/2-1)/ (BIT_MWORD/2)); //the number of machine words needed to be considered to contain the given overlap
    assert(N);
    
    int n = SR_WORD-N;  //index of the machine word that the overlap starts
    assert(n>=0 && n<SR_WORD);
    
    //the bit vectors are right aligned. So the first word needed to be shifted to get the only required part
    seqA = readA->sequence[n];
    seqA = seqA << ((N*BIT_MWORD/2-overlap)*2);
    seqB = readB->sequence[n];
    seqB = seqB << ((N*BIT_MWORD/2-overlap)*2);
    //if the first words is the lastword or if the first words itself are different
    if(seqA!=seqB || N==1){
        return ((seqA>seqB) - (seqA<seqB));    
    }
        
    
    //go through all remianing words
    int i=n+1; 
    assert(i<SR_WORD);
    do{
        seqA = readA->sequence[i];
        seqB = readB->sequence[i];
        if(i==SR_WORD-1){
            return ((seqA>seqB) - (seqA<seqB));
        }
        i++;    
    } while(seqA==seqB);
        
    return ((seqA>seqB) - (seqA<seqB));    
    
}    

//compare function for sorting by left corner
/*int sortbyFirstBases(const void *a,const void *b){
    struct read *readA = (struct read *)a;
    struct read *readB = (struct read *)b;
    return strncmp(readA->sequence,readB->sequence,OVERLAP_MAX);
    
}*/

//compare function for sorting by left corner
int sortbyFirstBases(const void *a,const void *b){
    struct read *readA = (struct read *)a;
    struct read *readB = (struct read *)b;
    
    MWORD seqA;
    MWORD seqB;
    
    int n = (OVERLAP_MAX+BIT_MWORD/2-1)/ (BIT_MWORD/2); //the number of machine words needed to be considered to contain the given overlap

    int i=0;
    do{
        seqA = readA->sequence[i];
        seqB = readB->sequence[i];
        if(i==n-1){
            return ((seqA>seqB) - (seqA<seqB));
        }
        i++;    
    } while(seqA==seqB);
        
    return ((seqA>seqB) - (seqA<seqB));
    
    
    //commented code is for a maximum of 32 bases of overlap
    /*if(OVERLAP_MAX<=BIT_MWORD/2){
        seqA = readA->sequence[0];
        //seqA = seqA >> (BIT_MWORD-OVERLAP_MAX*2);
        seqB = readB->sequence[0];
        //seqB = seqB >> (BIT_MWORD-OVERLAP_MAX*2);
        return ((seqA>seqB) - (seqA<seqB));
        
        
    }
    else{
        fprintf(stderr,"Please implement multiple machine words\n");
        exit(EXIT_FAILURE);
    }*/
    
    //return strncmp(readA->sequence,readB->sequence,OVERLAP_MAX);
    
} 



    
//compare two overlaps (one from left sligbed and the other from right aligned)    
int mycmp(MWORD *target, MWORD *query,int overlap){
    
    //targets left aligned
    //qeries are aligned by last bases
    
    MWORD seq_target;
    MWORD seq_query;
    
    //commented code is valid only for a maximum overlpa size of 32 bases
    /*if(overlap<=BIT_MWORD/2){
        seq_target = target[0];
        
        seq_query = query[SR_WORD-1];
        
        seq_target = seq_target >> (BIT_MWORD-overlap*2);

        seq_query = ((seq_query << (BIT_MWORD-overlap*2))>>(BIT_MWORD-overlap*2));   
    }
    else{
        fprintf(stderr,"Please implement multiple machine words\n");
        exit(EXIT_FAILURE);
    }    */
        
    int n = (overlap+BIT_MWORD/2-1)/ (BIT_MWORD/2); //the number of machine words needed to be considered to contain the given overlap

    int i=0;
    do{
        assert(SR_WORD-n+i < SR_WORD);
        seq_target = target[i];
        //since right aligned, need to take end part of the current words and starting part of the next word
        seq_query = (query[SR_WORD-n+i]<< ((n*BIT_MWORD/2-overlap)*2)) | (query[SR_WORD-n+i+1]>> ((overlap-(n-1)*BIT_MWORD/2)*2)); 
        if(i==n-1){ //if this is the last byte  things are a bit different
            seq_target =  seq_target >> ((n*BIT_MWORD/2-overlap)*2); //only the starting part of the left aligned read is neeeded
            seq_query = (query[SR_WORD-n+i] << ((n*BIT_MWORD/2-overlap)*2)) >> ((n*BIT_MWORD/2-overlap)*2); //only th eend part of the right aligned read is needed
            return ((seq_target>seq_query) - (seq_target<seq_query));
        }
        i++;    
    } while(seq_target==seq_query);
        
    return ((seq_target>seq_query) - (seq_target<seq_query));

}    


void linearsort(struct read **queryReads_addr,struct read **tempReads_addr){
    
    struct read *queryReads = *queryReads_addr;
    struct read *tempReads = *tempReads_addr;
    
    seqid_t count00=0;
    seqid_t count01=0;
    seqid_t count10=0;
    seqid_t count11=0;
    
    MWORD seqA;
    
    int i;
    for(i=0;i<NUM_READS*2;i++){
        
        int N = ((overlap+BIT_MWORD/2-1)/ (BIT_MWORD/2)); //the number of machine words needed to be considered to contain the given overlap
        assert(N);
        int n = SR_WORD-N;  //index of the machine word that the overlap starts
        assert(n>=0 && n<SR_WORD);
        
        //the bit vectors are right aligned. So the first word needed to be shifted to get the only required part
        seqA = queryReads[i].sequence[n];
        //seqA = seqA >> (BIT_MWORD-(N*BIT_MWORD-overlap*2)-2);
        seqA= (seqA << (N*BIT_MWORD-overlap*2))>>(BIT_MWORD-2);  //optimise
        
        switch(seqA){
            case 0:
                count00++;
                break;            
            case 1: 
                count01++;
                break;
            case 2: 
                count10++;
                break;               
            case 3: 
                count11++;
                break;
            default:
                fprintf(stderr,"Some improproper encoded base %llx",seqA);
                exit(EXIT_FAILURE);          
        }
        
    }

    //fprintf(stderr,"%d %d %d %d\n",count00,count01,count10,count11);
    assert(count00+count01+count10+count11==NUM_READS*2);    
    
    int i00=0;
    int i01=count00;
    int i10=count00+count01;
    int i11=count00+count01+count10;
    
    for(i=0;i<NUM_READS*2;i++){
        
        int N = ((overlap+BIT_MWORD/2-1)/ (BIT_MWORD/2)); //the number of machine words needed to be considered to contain the given overlap
        assert(N);
        int n = SR_WORD-N;  //index of the machine word that the overlap starts
        assert(n>=0 && n<SR_WORD);
        
        //the bit vectors are right aligned. So the first word needed to be shifted to get the only required part
        seqA = queryReads[i].sequence[n];
        //seqA = seqA >> (BIT_MWORD-(N*BIT_MWORD-overlap*2)-2);
        seqA= (seqA << (N*BIT_MWORD-overlap*2))>>(BIT_MWORD-2);  //optimise        
 
        switch(seqA){
            case 0:
                tempReads[i00]=queryReads[i];
                i00++;
                break;            
            case 1: 
                tempReads[i01]=queryReads[i];
                i01++;
                break;
            case 2: 
                tempReads[i10]=queryReads[i];
                i10++;
                break;               
            case 3: 
                tempReads[i11]=queryReads[i];
                i11++;
                break;
            default:
                fprintf(stderr,"Some improproper encoded base %llx",seqA);
                exit(EXIT_FAILURE);          
        }
 
    }
    //fprintf(stderr,"%d %d %d %d\n",i00,i01,i10,i11);
    assert(i00==count00);
    assert(i01==count00+count01);
    assert(i10==count00+count01+count10);
    assert(i11==NUM_READS*2);
    
    
    //change pointers
    *queryReads_addr = tempReads;
    *tempReads_addr = queryReads;
}
    
void olp(int argc, char **argv){
	 
    if(argc!=7){
       fprintf(stderr,"Please enter command as : ./biorand %s <in.fq> <out.paf> <read size> <num reads> <min overlap> <max overlap>\n",argv[0]);
       exit(EXIT_FAILURE);
    }    
      
    INITIALISE_TIMERS 
     
    //setting variables
    READ_SIZE= atoi(argv[3]);
    if(READ_SIZE>MAX_READ_SIZE){
        fprintf(stderr,"Max read size allowed is %d\n",MAX_READ_SIZE);
        exit(EXIT_FAILURE);
    }
    NUM_READS= atoi(argv[4]);
    OVERLAP_MIN= atoi(argv[5]);        
    OVERLAP_MAX= atoi(argv[6]);
    overlap = OVERLAP_MIN;
    fprintf(stderr,"Running for %d reads of read size %d bases\nOverlap min is %d bases and overlap max is %d bases\n",NUM_READS,READ_SIZE,OVERLAP_MIN,OVERLAP_MAX);
    

    if (BIT_MWORD != sizeof(MWORD)*8){
        fprintf(stderr,"size of a MWORD is stated as %u bits, but actually it is %u bits on this arhcitecture\n",(unsigned int)BIT_MWORD,(unsigned int)(sizeof(MWORD)*8 ));
        exit(EXIT_FAILURE);
    }    
    
    //calculate parameters
    SR_WORD = (READ_SIZE+BIT_MWORD/2-1)/(BIT_MWORD/2);
    SR_CONV_WORD = (READ_SIZE+7)/8;
    fprintf(stderr,"BIT_MWORD %u\nSR_WORD %u\nSR_CONV_WORD %u\n",BIT_MWORD,SR_WORD,SR_CONV_WORD);   
    
    //opening files 
    FILE *input = fopen(argv[1],"r");
    errorCheckNULL(input);
    
    FILE *output = fopen(argv[2],"w");
    errorCheckNULL(output);   
    
#ifdef DEBUG
    FILE *readsfile = fopen("readsfile.txt","w");
    errorCheckNULL(readsfile);
    FILE *sortedreadsfile = fopen("sortedreads.txt","w");
    errorCheckNULL(sortedreadsfile);
    FILE *sortedreadsleftfile = fopen("sortedreadsleft.txt","w");
    errorCheckNULL(sortedreadsleftfile);    
    FILE *mappings = fopen("mappings.txt","w");
    errorCheckNULL(mappings);        
    char sequence_ascii[MAX_READ_SIZE];
#endif
    
    
    //buffers for getline
    char *buffer = (char *)malloc(sizeof(char)*(READ_SIZE+2));  //READ+newline+nullcharacter
    errorCheckNULL(buffer);
    size_t bufferSize = READ_SIZE+2;
    //buffer for bit conversion (can be optimised)
    MWORD *seq_ascii = (MWORD *)malloc(sizeof(MWORD)*SR_CONV_WORD); errorCheckNULL(seq_ascii);
    
    //space for targetReads (sorted by left corner) (*2 for reverse complement)
    struct read *targetReads = (struct read *)malloc(sizeof(struct read)*NUM_READS*2);
    errorCheckNULL(targetReads);
    //space for queryReads (sorted by right corner) (*2 for reverse complement)
    struct read *queryReads = (struct read *)malloc(sizeof(struct read)*NUM_READS*2);
    errorCheckNULL(queryReads);

#ifdef LINEAR_SORT
     //temporary space
    struct read *tempReads = (struct read *)malloc(sizeof(struct read)*NUM_READS*2);
    errorCheckNULL(tempReads);   
#endif    
    
    
    
    seqid_t i;
    int readlinebytes;
    
    START_TIMER
    
    //Read fastq file and load the arrays
    for(i=0;i<NUM_READS;i++){
        
        /**********read Name*******************/
        readlinebytes=getline(&buffer, &bufferSize, input); 
        //file has ended
        if(readlinebytes == -1){
            fprintf(stderr,"file has ended prematurely\n");
            exit(EXIT_FAILURE);
        }  
        if(readlinebytes>MAX_READ_ID_SIZE){
            fprintf(stderr,"We expected a maximum readID of %d, but got %d\n",MAX_READ_ID_SIZE,readlinebytes);
            exit(EXIT_FAILURE);
        }
        if(readlinebytes<3){
            fprintf(stderr,"Malformed read ID. Read ID should be at least 1 letter\n");
            exit(EXIT_FAILURE);
        }
        if(buffer[0]!='@'){
            fprintf(stderr,"Malformed read ID. Read ID should start with @\n");
            exit(EXIT_FAILURE);
        }
        if(buffer[readlinebytes-1]=='\n'){  //check for \r as well and other funny characters
            buffer[readlinebytes-1]=0;
        }        
        strcpy(targetReads[i*2].read_id,&buffer[1]); //forward
        strcpy(targetReads[i*2+1].read_id,&buffer[1]); //reverse complement
        
        /**************sequence************************/
        readlinebytes=getline(&buffer, &bufferSize, input); 
        //file has ended
        if(readlinebytes == -1){
            fprintf(stderr,"file has ended prematurely\n");
            exit(EXIT_FAILURE);
        }           
        if(readlinebytes!=READ_SIZE+1){
            fprintf(stderr,"We expected a read size of %d, but got %d\n",READ_SIZE,readlinebytes);
            exit(EXIT_FAILURE);
        }
        if(buffer[readlinebytes-1]=='\n'){  //check for \r as well and other funny characters
            buffer[readlinebytes-1]=0;
        }
        //strcpy(targetReads[i*2].sequence,buffer);
        strncpy((char *)seq_ascii,buffer,READ_SIZE);  //copy
        Base2Bit(targetReads[i*2].sequence, seq_ascii);
        
        //reverse_complement(targetReads[i*2+1].sequence,buffer);
        reverse_complement((char *)seq_ascii,buffer); //optimise
        Base2Bit(targetReads[i*2+1].sequence, seq_ascii);
        
        //setting if forward or reversed complement
        targetReads[i*2].reversecomplement=0; //forward
        targetReads[i*2+1].reversecomplement=1; //reverse complement
        //setting the seq ID
        targetReads[i*2].seqid=i; 
        targetReads[i*2+1].seqid=i; 
        
        //+ sign : ignore
        readlinebytes=getline(&buffer, &bufferSize, input); 
        //file has ended
        if(readlinebytes == -1){
            fprintf(stderr,"file has ended prematurely\n");
            exit(EXIT_FAILURE);
        }  
        
  
        //Quality values : ignore
        readlinebytes=getline(&buffer, &bufferSize, input); 
        //file has ended
        if(readlinebytes == -1){
            fprintf(stderr,"file has ended prematurely\n");
            exit(EXIT_FAILURE);
        }          
    
    }
    
    STOP_TIMER
    PRINT_TIME("Reading files completed ... ")
	
#ifdef DEBUG
    for(i=0;i<NUM_READS*2;i++){
        Print_Bitvector(sequence_ascii,targetReads[i].sequence);
        fprintf(readsfile,"%s %s %d\n",targetReads[i].read_id, sequence_ascii, targetReads[i].reversecomplement);    
    } 
#endif   

    
    START_TIMER
    //copying the targetReads as query reads
    memcpy(queryReads,targetReads,sizeof(struct read)*NUM_READS*2);  
    //shift all bit vectors to the right
    for(i=0;i<NUM_READS*2;i++){
        Shift_Right( SR_WORD*BIT_MWORD-READ_SIZE*2,   queryReads[i].sequence);
    }
    STOP_TIMER
    PRINT_TIME("Mem copy finished ... ")
    

    
    //do the sorting  of target reads by left corner
    START_TIMER
    qsort(targetReads,NUM_READS*2,sizeof(struct read),sortbyFirstBases);
    STOP_TIMER
    PRINT_TIME("Sorted reads by left corner ... ")    


 
#ifdef DEBUG
    for(i=0;i<NUM_READS*2;i++){
        Print_Bitvector(sequence_ascii,targetReads[i].sequence);
        fprintf(sortedreadsfile,"%s %s %d\n",targetReads[i].read_id, sequence_ascii, targetReads[i].reversecomplement);    
    }
#endif
  

    int k;
    //perform the search for mappings
    for(k=0;k<OVERLAP_MAX-OVERLAP_MIN+1;k++){
        
        fprintf(stderr,"\nStep %d out of %d ... \n",k+1,OVERLAP_MAX-OVERLAP_MIN+1);
        pairs=0; //set counter to 0
    
    #ifdef DEBUG
        fprintf(mappings,"Overlap by %d bases\n",overlap);
    #endif    
    
        //sort query reads by right corner (in each loop extend sorting by one base)
        START_TIMER
#ifdef LINEAR_SORT
        if(k==0){
            qsort(queryReads,NUM_READS*2,sizeof(struct read),sortbyLastBases);
        }
        else{
            linearsort(&queryReads,&tempReads);
        }
#else          
        qsort(queryReads,NUM_READS*2,sizeof(struct read),sortbyLastBases);
#endif
        STOP_TIMER
        PRINT_TIME("Sorted reads by right corner by %d bases ... ",overlap)         
        
        START_TIMER
        
        seqid_t marker=0,j;  //marker stores the index for query reads
        //iterate through the target reads while finding mappins in the query reads
        for(i=0;i<NUM_READS*2;i++){
            
            int comp=1;
            //skip untill we come to a place where we have an equality or passed a probable equality
            while(comp>0 && marker<NUM_READS*2){
                comp=mycmp(targetReads[i].sequence,queryReads[marker].sequence,overlap);
                marker++;
            }
            j=marker; //copy the current marker
            marker--; //cancel the post increment in the loop 
            
            //if the last loop broke due to an equality need to scan further
            if(comp==0){
            #ifdef DEBUG
                char targetRead_ascii[MAX_READ_SIZE];
                char queryRead_ascii[MAX_READ_SIZE];
                Print_Bitvector(targetRead_ascii,targetReads[i].sequence);
                Print_Bitvector_right_shifted(queryRead_ascii,queryReads[marker].sequence);
                fprintf(mappings,"%c%s %c%s | %s(%d) %s(%d)\n",(queryReads[marker].reversecomplement==1 ? '-' : '+'), queryRead_ascii, \
                    (targetReads[i].reversecomplement==1 ? '-' : '+'),targetRead_ascii, queryReads[marker].read_id,marker, targetReads[i].read_id ,i);
            #endif                 
                print_to_paf(output, targetReads, i, queryReads, marker, overlap);
                while(j<NUM_READS*2){
                    comp=mycmp(targetReads[i].sequence,queryReads[j].sequence,overlap);
                    if(comp==0){
                    #ifdef DEBUG
                    Print_Bitvector(targetRead_ascii,targetReads[i].sequence);
                    Print_Bitvector_right_shifted(queryRead_ascii,queryReads[j].sequence);                    
                    fprintf(mappings,"%c%s %c%s | %s(%d) %s(%d)\n",(queryReads[j].reversecomplement==1 ? '-' : '+'), queryRead_ascii, \
                        (targetReads[i].reversecomplement==1 ? '-' : '+'),targetRead_ascii, queryReads[j].read_id,j, targetReads[i].read_id ,i);
                    #endif 
                        print_to_paf(output, targetReads, i, queryReads, j, overlap);
                    }
                    else if(comp<0){
                        break;
                    }
                    else{
                        fprintf(stderr,"cant happen. Probably sorting issue\n");
                        exit(EXIT_FAILURE);
                    }
                    j++;
                }            
            }
            
            
 
        }
        
        STOP_TIMER
        PRINT_TIME("%ld pairs found ... ",pairs)           

        
        
    #ifdef DEBUG
        fprintf(mappings,"\n");
    #endif         
        overlap++;

        
    }
    
    
#ifdef DEBUG
    for(i=0;i<NUM_READS*2;i++){
        Print_Bitvector_right_shifted(sequence_ascii,queryReads[i].sequence);
        fprintf(sortedreadsleftfile,"%s %s %d\n",queryReads[i].read_id,sequence_ascii,targetReads[i].reversecomplement);    
    }    
    
    fclose(readsfile);
    fclose(sortedreadsfile);
    fclose(sortedreadsleftfile);
    fclose(mappings);
#endif    

#ifdef LINEAR_SORT
    free(tempReads);   
#endif    
    
    free(buffer);
    free(targetReads);
    free(queryReads);
    fclose(input);
    fclose(output);
    
	//return 0;
}
