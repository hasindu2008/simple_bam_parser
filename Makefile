CC     = gcc
LIBS = -Lhtslib -lz -lm -lbz2 -llzma -lpthread 
CFLAGS   = -g -Wall -O2
INC = -Ihtslib

all: randomacess.c htslib/libhts.a  sequentialaccess.c common.c common.h
	$(CC) randomacess.c common.c htslib/libhts.a $(CLFLAGS) $(INC) $(LIBS) -o randomacess
	$(CC) sequentialaccess.c common.c htslib/libhts.a $(CLFLAGS) $(INC) $(LIBS) -o sequentialaccess

clean: 
	rm randomacess sequentialaccess
