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


#define MAX_READNAME_LEN 100  //The maximun size of a read name (qname size in bytes +1 )
#define MAX_READ_LEN 150      //maximum size of a read (number of bases +1)   
#define MAX_N_CIGAR 16        //no idea what this number of CIGAR ops mean at the moment  


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
    

struct alignedRead {
    
  char qname[MAX_READNAME_LEN];         //Query template NAME
  uint32_t flag;                        //bitwise FLAG
  int32_t chromID;                      //References sequence ID (chr1 is 0, chr2 is 1 and so on)
  uint32_t pos;                         //0-based leftmost mapping POSition
  uint8_t mapq;                         //MAPping Quality 
  uint32_t cigarOps[2*MAX_N_CIGAR];     //CIGAR ops            
  uint32_t mateChromID;                 //Ref. ID of the mate/next read
  uint32_t matePos;                     //Position of the mate/next read
  uint32_t tlen;                        //observed Template LENgth (not used)
  char seq[MAX_READ_LEN];               //segment SEQuence
  uint8_t qual[MAX_READ_LEN];           //quality string   

  uint32_t cigarLen;
  uint32_t rlen;                        //Length of SEQuence (l seq)
  uint32_t end;
  uint32_t insertSize;
 
  //need to get some other fields in the BAM

};    
   
//Is this efficient? Should try optimising   
char _getBase(uint8_t *s, int i){
    char* baseLookup="=ACMGRSVTWYHKDBN";
    return baseLookup[bam_seqi(s, i)];
}        
   
   
struct alignedRead* getRead(struct alignedRead*  theRead, bam1_t *b/*, int storeRgID, char** rgID*/){
    
        bam1_core_t *c = &(b->core);
        uint8_t *s = bam_get_seq(b);
        uint8_t *q = bam_get_qual(b);
        int32_t lenSeq = c->l_qseq;
        char *readname = bam_get_qname(b);
        
        int read_name_length=strlen(readname);
        if(MAX_READNAME_LEN< read_name_length+1){
            fprintf(stderr,"The maximum read name length is set to %d, but the actual read length is %d\n",MAX_READNAME_LEN,read_name_length + 1); 
            exit(EXIT_FAILURE);            
        }
        
        
        if (lenSeq == 0){
            fprintf(stderr,"The sequence length is 0. How come?\n"); 
            exit(EXIT_FAILURE);
        }
        
        if (q[0] == 0xff){
            fprintf(stderr,"The quality score is 255 for the first base. How come?\n"); 
            exit(EXIT_FAILURE);
        }
        
        //struct alignedRead* theRead = (struct alignedRead*)malloc(sizeof(struct alignedRead));

        
        if(MAX_READ_LEN<lenSeq+1){
            fprintf(stderr,"The maximum read length is set to %d, but the actual read length is %d\n",MAX_READ_LEN,lenSeq + 1); 
            exit(EXIT_FAILURE);
        }
        
        //char* seq             = (char*)malloc((lenSeq + 1) * sizeof(char));
        //uint8_t* qual            = (uint8_t*)malloc((lenSeq + 1) * sizeof(uint8_t));
        
        //assert (theRead != NULL);
        //assert (seq     != NULL);
        //assert (qual    != NULL);
        
        // Try to grab the read-group tag value
        /*uint8_t* v     = NULL;
        char* tempRgID = NULL;
        int lenRgID    = 0;
        if (storeRgID){
            v = bam_aux_get(b, "RG");
            tempRgID = bam_aux2Z(v);
            lenRgID = strlen(tempRgID);
            rgID[0] = (char*)(calloc(lenRgID + 1, sizeof(char)));
            strcpy(rgID[0], tempRgID);
        }*/
        
        int i = 0;
        for (i=0; i < lenSeq; i++){
            theRead->seq[i] = _getBase(s, i);
            theRead->qual[i] = q[i];
            assert (theRead->qual[i] <= 93);
            assert (theRead->qual[i] >= 0);
        }
        theRead->seq[lenSeq]  = '\0';
        theRead->qual[lenSeq] = '\0';
        
        int32_t readStart = c->pos; 
        
        if(MAX_N_CIGAR< (c->n_cigar)){
            fprintf(stderr,"The maximum number of cigar is set to %d, but the actual number of cigar is %d\n",MAX_N_CIGAR,c->n_cigar); 
            exit(EXIT_FAILURE);           
        }
        
        //uint32_t* cigarOps = (uint32_t*)malloc(2 * c->n_cigar * sizeof(uint32_t));
        //assert (cigarOps != NULL);
        uint32_t *cigar = bam_get_cigar(b);
        for (i=0 ;i < c->n_cigar; i++){
            uint32_t cigarFlag     = bam_cigar_op(cigar[i]);
            uint32_t cigarFlagLen  = bam_cigar_oplen(cigar[i]);
            theRead->cigarOps[2 * i]       = cigarFlag;
            theRead->cigarOps[(2 * i) + 1] = cigarFlagLen;
            
            // Soft-clipping of sequence at start of read changes the mapping
            // position. Recorded mapping pos is that of the first aligned (not soft-clipped)
            // base. I want to adjust this so that the read start refers to the first base in
            // the read.
            //if (i == 0 && cigarFlag == 4){
                //readStart -= cigarFlagLen;
            //}
        }
        
        strcpy(theRead->qname, readname);
        theRead->flag        = c->flag;
        theRead->chromID     = c->tid;
        theRead->pos         = readStart;        
        theRead->mapq        = c->qual;
        //theRead->cigarOps    = cigarOps;
        theRead->mateChromID = c->mtid; 
        theRead->matePos     = c->mpos;
        theRead->tlen        = 0;      
        //theRead->seq         = seq;
        //theRead->qual        = qual;


        theRead->cigarLen    = c->n_cigar;
        theRead->rlen        = lenSeq;
        theRead->end         = bam_endpos(b);
        theRead->insertSize  = c->isize;

        
        return theRead;
}     


void printRead(struct alignedRead* theRead){
    
    //print in sam like format
    //need to convert chromID to chromosome name
    //need to convert cigarOps to CIGAR string
    //need to convert qual to ASCII
    //need to get optional tags
    printf("%s\t%u\t%u\t%u\t%u\t%s\t%u\t%u\t%u\t%s\t%s\t%s\n",
            theRead->qname, theRead->flag, theRead->chromID, (theRead->pos+1), //+1 to convert to 1 based index
            theRead->mapq, "CIGAR", theRead->mateChromID, (theRead->matePos+1), 
            0, theRead->seq, "QUAL", "OPTIONAL_TAGS");    
    
}  


    
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