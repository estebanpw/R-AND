CC=gcc
CXX=g++
CFLAGS=-O3 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -Wall #-DVERBOSE
BIN=../bin

all: DNArand



DNArand: DNArand.c
	$(CC) $(CFLAGS) alignmentFunctions.c -lm commonFunctions.c -lm DNArand.c -lm -o $(BIN)/DNArand

clean:
	rm -rf $(BIN)/DNArand
