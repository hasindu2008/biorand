#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <getopt.h>
#include "misc.h"

#include <vector>
#include <string>

//https://github.com/HenrikBengtsson/illuminaio/blob/master/R/readIDAT_nonenc.R

//RunInfo and midblock not read

//#define READ_UNOWN_FILEDS 0

#define IDAT_N_SNPS_READ 1000
#define IDAT_ILLUMINA_ID 102
#define IDAT_SD 103
#define IDAT_MEAN 104
#define IDAT_N_BEADS 107
#define IDAT_RED_GREEN 400
#define IDAT_BARCODE 402
#define IDAT_CHIP_TYPE 403

typedef struct{
    int64_t nSNPsRead_offset;
    int64_t IlluminaID_offset;
    int64_t SD_offset;
    int64_t Mean_offset;
    int64_t NBeads_offset;
    int64_t RedGreen_offset;
    int64_t Barcode_offset;
    int64_t ChipType_offset;
} idat_fields_t;


int64_t read_idat_string(char **dest,FILE *file){
    // From [1] https://code.google.com/p/glu-genetics/source/browse/glu/lib/illumina.py#86:
    // String data are encoded as a sequence of one or more length
    // bytes followed by the specified number of data bytes.
    //
    // The lower 7 bits of each length byte encodes the bits that
    // comprise the length of the following byte string.  When the
    // most significant bit it set, then an additional length byte
    // follows with 7 additional high bits to be added to the current
    // length.  The following string lengths are accommodated by
    // increasing sequences of length bytes:
    //
    // length  maximum
    // bytes   length
    // ------  --------
    //   1       127 B
    //   2        16 KB
    //   3         2 MB
    //   4       256 MB
    //   5        32 GB
    //
    // While this seems like a sensible progression, there is some
    // uncertainty about this interpretation, since the longest of
    // string observed in the wild has been of length 6,264 with
    // two length bytes.
    //
    // EXAMPLES: by HB (2015-09-11)
    // Length   Len bytes  Iterations (n, m, k, shift)
    // 1        (1)
    // 127      (127)       -> n=128
    // 128      (128,1)     -> n=0  ,m=1  ,shift=7, k=128 -> n=128
    // 255      (255,1)     -> n=128,m=1  ,shift=7, k=128 -> n=255
    // 256      (128,2)     -> n=0  ,m=2  ,shift=7, k=256 -> n=256
    // 257      (129,2)     -> n=1  ,m=2  ,shift=7, k=256 -> n=257
    // 512      (128,4)     -> n=0  ,m=4  ,shift=7, k=512 -> n=512
    // 81921    (129,128,5) -> n=1  ,m=128,shift=7, k=0 -> n=1+0=1
    //                      -> n=0  ,m=5  ,shift=14,k=81920
    //                                              -> n=1+81920=81921

    // Parse the number of characters to read
    uint8_t m,n;

    size_t ret = fread(&m,1,1,file); assert(ret == 1);
    n = m % 128;

    int64_t shift =0;
    int64_t k =0;

    while (m / 128 == 1) {
        // Read next length byte ...
        ret = fread(&m,1,1,file); assert(ret == 1);

        // ... which represents the next 7 hi-bits
        shift = shift + 7;
        k = (m % 128) * 2^shift;

        // Total number of bytes to read
        n = n + k;
    }

    //fprintf(stderr,"String size %d\n",n);
    // Now read all bytes/characters
    char *str = *dest;
    str = (char *)malloc(sizeof(char)*(n+1));

    ret = fread(str,1,n,file); assert(ret == n);
    str[n] = 0;

    *dest = str;

    return n;
}


int idat(int argc, char **argv){

    if (argc < 2) {
        fprintf(
            stderr,
            "Usage: biorand %s [OPTIONS] a.idat\n",
            argv[0]);
        exit(EXIT_FAILURE);
    }

    size_t ret;
    int reti;

    FILE *file = fopen(argv[1],"rb");
    F_CHK(file,argv[1]);

    char magic[5];
    ret = fread(magic,1,4,file);
    assert(ret == 4);
    magic[4]=0;

    if(strcmp(magic,"IDAT")!=0){
        ERROR("Not an IDAT file. Unknown magic: %s",magic);
        exit(EXIT_FAILURE);
    }

    int32_t version, unknown0, nFields;
    ret = fread(&version,4,1,file); assert(ret == 1);


    if(version != 3){
        ERROR("Only IDAT files of version 3 are supported. Your file had version: %d",version);
        exit(EXIT_FAILURE);
    }

    ret = fread(&unknown0,4,1,file); assert(ret==1);
    ret = fread(&nFields,4,1,file); assert(ret==1);

    //fprintf(stderr,"%d fields\n",nFields);

    uint16_t field_code;
    int64_t field_offset;
    idat_fields_t idat_fields = {0};

    int64_t all_others_offset[nFields];
    int32_t all_others_offset_index=0;

    for(int i=0;i<nFields;i++){
        ret = fread(&field_code,2,1,file); assert(ret==1);
        ret = fread(&field_offset,8,1,file);assert(ret==1);
        //fprintf(stderr,"%d\t%ld\n",field_code, field_offset);

        switch (field_code) {
            case IDAT_N_SNPS_READ:
                idat_fields.nSNPsRead_offset = field_offset;
                break;
            case IDAT_ILLUMINA_ID:
                idat_fields.IlluminaID_offset = field_offset;
            case IDAT_SD:
                idat_fields.SD_offset = field_offset;
            case IDAT_MEAN:
                idat_fields.Mean_offset = field_offset;
            case IDAT_N_BEADS:
                idat_fields.NBeads_offset = field_offset;
            case IDAT_RED_GREEN:
                idat_fields.RedGreen_offset = field_offset;
            case IDAT_BARCODE:
                idat_fields.Barcode_offset = field_offset;
            case IDAT_CHIP_TYPE:
                idat_fields.ChipType_offset = field_offset;
            default:
                all_others_offset[all_others_offset_index] = field_offset;
                all_others_offset_index++;
                break;

        }

    }

    int32_t nSNPsRead;

    assert(idat_fields.nSNPsRead_offset > 0);
    assert(idat_fields.IlluminaID_offset > 0);
    assert(idat_fields.SD_offset);
    assert(idat_fields.Mean_offset);
    assert(idat_fields.NBeads_offset);
    assert(idat_fields.RedGreen_offset);
    assert(idat_fields.Barcode_offset);
    assert(idat_fields.ChipType_offset);

    reti=fseek(file,idat_fields.nSNPsRead_offset,SEEK_SET); assert(reti==0);
    ret = fread(&nSNPsRead, 4,1, file); assert(ret==1);

    int32_t *IlluminaID = (int32_t *)malloc(sizeof(int32_t)*nSNPsRead);
    reti=fseek(file,idat_fields.IlluminaID_offset,SEEK_SET); assert(reti==0);
    ret = fread(IlluminaID, 4,nSNPsRead, file); assert(ret==(size_t)nSNPsRead);

    uint16_t *SD = (uint16_t *)malloc(sizeof(int32_t)*nSNPsRead);
    reti=fseek(file,idat_fields.SD_offset,SEEK_SET); assert(reti==0);
    ret = fread(SD, 2,nSNPsRead, file); assert(ret==(size_t)nSNPsRead);

    uint16_t *Mean = (uint16_t *)malloc(sizeof(int32_t)*nSNPsRead);
    reti=fseek(file,idat_fields.Mean_offset,SEEK_SET); assert(reti==0);
    ret = fread(Mean, 2,nSNPsRead, file); assert(ret==(size_t)nSNPsRead);

    uint8_t *NBeads = (uint8_t *)malloc(sizeof(int32_t)*nSNPsRead);
    reti=fseek(file,idat_fields.NBeads_offset,SEEK_SET); assert(reti==0);
    ret = fread(NBeads, 1,nSNPsRead, file); assert(ret==(size_t)nSNPsRead);

    int32_t RedGreen;
    reti=fseek(file,idat_fields.RedGreen_offset,SEEK_SET); assert(reti==0);
    ret = fread(&RedGreen, 4,1, file); assert(ret==1);


    char *chip_type=NULL;
    reti=fseek(file,idat_fields.ChipType_offset,SEEK_SET); assert(reti==0);
    read_idat_string(&chip_type,file);
    assert(chip_type!=NULL);

    char *barcode=NULL;
    reti=fseek(file,idat_fields.Barcode_offset,SEEK_SET); assert(reti==0);
    read_idat_string(&barcode,file);
    assert(chip_type!=NULL);

    fprintf(stdout,"#nSNPsRead\t%d\n",nSNPsRead);
    fprintf(stdout,"#ChipType\t%s\n",chip_type);
    fprintf(stdout,"#Barcode\t%s\n",barcode);
    fprintf(stdout,"#RedGreen\t%d\n",RedGreen);

    #ifdef READ_UNOWN_FILEDS
        char *all_others[all_others_offset_index] = {NULL};
        for(int i=0;i<all_others_offset_index;i++){
            reti=fseek(file,all_others_offset[i],SEEK_SET); assert(reti==0);
            read_idat_string(&all_others[i],file);
            assert(chip_type!=NULL);
            fprintf(stdout,"#Unknown\t%s\n",all_others[i]);
        }
    #endif


    fprintf(stdout,"#IlluminaID\tSD\tMean\tNBeads\n");
    for(int i=0;i<nSNPsRead;i++){
        fprintf(stdout,"%d\t%d\t%d\t%d\n",IlluminaID[i],SD[i], Mean[i], NBeads[i]);
    }
    fprintf(stderr,"\n");

    free(IlluminaID);
    free(SD);
    free(Mean);
    free(NBeads);
    free(chip_type);
    free(barcode);
    #ifdef READ_UNOWN_FILEDS
        for(int i=0;i<all_others_offset_index;i++){
            if(all_others[i]){
                free(all_others[i]);
            }
        }
    #endif

    fclose(file);

    return 0;

}
