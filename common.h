
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

char _getBase(uint8_t *s, int i);
struct alignedRead* getRead(struct alignedRead*  theRead, bam1_t *b/*, int storeRgID, char** rgID*/);
void printRead(struct alignedRead* theRead);
   