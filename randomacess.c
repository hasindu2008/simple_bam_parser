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
        fprintf(stderr, "Usage %s file.bam chr:start-stop\n",argv[0]);
        exit(EXIT_FAILURE);
    }
    
    //these come from htslib/sam.h
	hts_itr_t *iter=NULL;
	hts_idx_t *idx=NULL;
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
    
    //load the index file for BAM
	idx = sam_index_load(in, argv[1]);
	errorCheckNULL(idx);
    
    //the iterator for the BAM random access is probably initialised here. Note that we pass the region string to this function 
	iter  = sam_itr_querys(idx, header, argv[2]); 
	errorCheckNULL(iter);
    
    //should check what this initialisation is
	b = bam_init1();
    
    struct alignedRead* myread = (struct alignedRead*)malloc(sizeof(struct alignedRead));
    
    //repeat until all reads in the regions are retrieved
	while ( sam_itr_next(in, iter, b) >= 0){
        getRead(myread, b);
        printRead(myread);  
	}
     free(myread);
    
    //wrap up
	hts_idx_destroy(idx);
	hts_itr_destroy(iter);
	bam_destroy1(b);
	bam_hdr_destroy(header);
	sam_close(in);
    
    return 0;
}