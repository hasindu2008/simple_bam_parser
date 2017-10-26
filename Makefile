CC     = gcc
LIBS = -Lhtslib -lz -lm -lbz2 -llzma -lpthread 
CFLAGS   = -g -Wall -O2
INC = -Ihtslib

all:
	$(CC) main.c htslib/libhts.a $(CLFLAGS) $(INC) $(LIBS) -o main

clean: 
	rm main
