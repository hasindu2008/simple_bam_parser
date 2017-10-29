#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <inttypes.h>
#include <stdbool.h>
#include <assert.h>
#include <errno.h>
#include "htslib/sam.h"
#include "common.h"

    
int main(int argc,char** argv){
    
    //check args
    if(argc!=3){
        fprintf(stderr, "Usage %s in.bam out.bam\n",argv[0]);
        exit(EXIT_FAILURE);
    }
    
    //these come from htslib/sam.h
	samFile *in = NULL;
	samFile *out = NULL;
	bam1_t *b= NULL;
    bam_hdr_t *header = NULL;

    //open the BAM file (though called sam_open is opens bam files too :P)
    in = sam_open(argv[1], "r");       //for reading
    errorCheckNULL(in);
    out = sam_open(argv[2], "w");      //for writing 
    errorCheckNULL(out);
    
    //get the sam header.
    if ((header = sam_hdr_read(in)) == 0){
        fprintf(stderr,"No sam header?\n");
        exit(EXIT_FAILURE);
    }
    //write the SAM header
    int ret=sam_hdr_write(out,header);
    errorCheck(ret);

    
    //this must be the initialisation for the structure that stores a read (need to verify)
	b = bam_init1();
    
    //my structure for a read (see common.h)
    struct alignedRead* myread = (struct alignedRead*)malloc(sizeof(struct alignedRead));
    
    //repeat until all reads in the file are retreived
	while ( sam_read1(in, header, b) >= 0){ 
        //write the read
        ret=sam_write1(out,header,b);  
        errorCheck(ret);
	}
    
    //wrap up
    free(myread);    
	bam_destroy1(b);
	bam_hdr_destroy(header);
	sam_close(in);
	sam_close(out);
    
    return 0;
}