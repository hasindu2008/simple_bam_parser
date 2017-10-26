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
    
    
  char *qname;             //Query template NAME
  uint32_t flag;           //bitwise FLAG
  int32_t chromID;         //References sequence ID (chr1 is 0, chr2 is 1 and so on)
  uint32_t pos;            //0-based leftmost mapping POSition
  uint8_t mapq;            //MAPping Quality 
  uint32_t *cigarOps;      //CIGAR ops            
  uint32_t mateChromID;    //Ref. ID of the mate/next read
  uint32_t matePos;        //Position of the mate/next read
  uint32_t tlen;            //observed Template LENgth (not used)
  char *seq;               //segment SEQuence
  uint8_t *qual;           //quality string   

  uint32_t cigarLen;
  uint32_t rlen;            //Length of SEQuence (l seq)
  uint32_t end;
  uint32_t insertSize;
 
  //uint32_t *hash;

};    
   
char _getBase(uint8_t *s, int i){
    char* baseLookup="=ACMGRSVTWYHKDBN";
    return baseLookup[bam_seqi(s, i)];
}        
   
struct alignedRead* get(bam1_t *b /*, int storeRgID, char** rgID*/){
    
        bam1_core_t *c = &(b->core);
        uint8_t *s = bam_get_seq(b);
        uint8_t *q = bam_get_qual(b);
        int32_t lenSeq = c->l_qseq;
        
        if (lenSeq == 0){
            return NULL;
        }
        
        if (q[0] == 0xff){
            return NULL;
        }
        
        struct alignedRead* theRead = (struct alignedRead*)malloc(sizeof(struct alignedRead));
        char* seq             = (char*)malloc((lenSeq + 1) * sizeof(char));
        uint8_t* qual            = (uint8_t*)malloc((lenSeq + 1) * sizeof(uint8_t));
        
        assert (theRead != NULL);
        assert (seq     != NULL);
        assert (qual    != NULL);
        
        // Try to grab the read-group tag value
        uint8_t* v     = NULL;
        char* tempRgID = NULL;
        int lenRgID    = 0;
        /*if (storeRgID){
            v = bam_aux_get(b, "RG");
            tempRgID = bam_aux2Z(v);
            lenRgID = strlen(tempRgID);
            rgID[0] = (char*)(calloc(lenRgID + 1, sizeof(char)));
            strcpy(rgID[0], tempRgID);
        }*/
        
        int i = 0;
        for (i=0; i < lenSeq; i++){
            seq[i] = _getBase(s, i);
            qual[i] = q[i];
            assert (qual[i] <= 93);
            assert (qual[i] >= 0);
        }
        seq[lenSeq]  = '\0';
        qual[lenSeq] = '\0';
        
        int32_t readStart = c->pos;
        uint32_t* cigarOps = (uint32_t*)malloc(2 * c->n_cigar * sizeof(uint32_t));
        assert (cigarOps != NULL);
        uint32_t *cigar = bam_get_cigar(b);
        for (i=0 ;i < c->n_cigar; i++){
            uint32_t cigarFlag     = bam_cigar_op(cigar[i]);
            uint32_t cigarFlagLen  = bam_cigar_oplen(cigar[i]);
            cigarOps[2 * i]       = cigarFlag;
            cigarOps[(2 * i) + 1] = cigarFlagLen;
            
            // Soft-clipping of sequence at start of read changes the mapping
            // position. Recorded mapping pos is that of the first aligned (not soft-clipped)
            // base. I want to adjust this so that the read start refers to the first base in
            // the read.
            if (i == 0 && cigarFlag == 4){
                readStart -= cigarFlagLen;
            }
        }
        
        theRead->qname       = bam_get_qname(b);
        theRead->seq         = seq;
        theRead->qual        = qual;
        theRead->cigarOps    = cigarOps;
        //theRead->hash        = NULL;
        theRead->mateChromID = c->mtid;
        theRead->cigarLen    = c->n_cigar;
        theRead->chromID     = c->tid;
        theRead->rlen        = lenSeq;
        theRead->pos         = readStart;
        theRead->end         = bam_endpos(b);
        theRead->insertSize  = c->isize;
        theRead->matePos     = c->mpos;
        theRead->flag     = c->flag;
        theRead->mapq        = c->qual;
        
        //Read_SetUnCompressed(theRead);
        printf("%s\t%u\t%u\t%u\t%u\t%s\t%u\t%u\t%u\t%s\t%s\t%s\n",
                theRead->qname, theRead->flag, theRead->chromID, theRead->pos,
                theRead->mapq, "CIGAR", theRead->mateChromID, theRead->matePos, 
                0, theRead->seq, "QUAL", "OPTIONAL_TAGS");
        
        return theRead;
}       
    
int main(int argc,char** argv){
    
    if(argc!=3){
        fprintf(stderr, "Usage %s file.bam chr:start-stop\n",argv[0]);
        exit(EXIT_FAILURE);
    }
    
    
	hts_itr_t *iter=NULL;
	hts_idx_t *idx=NULL;
	samFile *in = NULL;
	bam1_t *b= NULL;
    bam_hdr_t *header = NULL;

    in = sam_open(argv[1], "r");
    errorCheckNULL(in);
    
    if ((header = sam_hdr_read(in)) == 0){
        fprintf(stderr,"No sam header?\n");
        exit(EXIT_FAILURE);
    }
    
	idx = sam_index_load(in,  argv[1]);
	errorCheckNULL(idx);
    
	iter  = sam_itr_querys(idx, header, argv[2]); 
	errorCheckNULL(iter);
    
    
    
	b = bam_init1();
    
	while ( sam_itr_next(in, iter, b) >= 0){
        struct alignedRead* myread= get(b);
        free(myread);
        
	}
    
	hts_itr_destroy(iter);
	bam_destroy1(b);
	bam_hdr_destroy(header);
	sam_close(in);
    
    return 0;
}