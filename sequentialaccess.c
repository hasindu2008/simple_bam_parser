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
    if(argc!=2){
        fprintf(stderr, "Usage %s file.bam \n",argv[0]);
        exit(EXIT_FAILURE);
    }
    
    //these come from htslib/sam.h
	samFile *in = NULL;
	bam1_t *b= NULL;
    bam_hdr_t *header = NULL;

    //open the BAM file (though called sam_open is opens bam files too :P)
    in = sam_open(argv[1], "r");
    errorCheckNULL(in);
    
    //get the sam header. Need to check how to parse this
    if ((header = sam_hdr_read(in)) == 0){
        fprintf(stderr,"No sam header?\n");
        exit(EXIT_FAILURE);
    }
    
    
    //should check what this initialisation is
	b = bam_init1();
    
    struct alignedRead* myread = (struct alignedRead*)malloc(sizeof(struct alignedRead));
    
    //repeat until all reads in the file
	while ( sam_read1(in, header, b) >= 0){
        getRead(myread, b);
        printRead(myread);  
	}
     free(myread);
    
    //wrap up
	bam_destroy1(b);
	bam_hdr_destroy(header);
	sam_close(in);
    
    return 0;
}