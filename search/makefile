CC = gcc
CFLAGS = -Wall -m64 -D_LARGEFILE64_SOURCE=1 -D_FILE_OFFSET_BITS=64 -O9 -std=c99 -fPIC -fopenmp -DFM_COMP_64 #-DSA_64 -DCSALIB_SEARCH -I $(CUDA_LOCATION)/include/ -L /usr/local/cuda/lib64/ -lcudart -msse4.2 #-DVERBOSE_DBG

DEPEN = io.o csafm.o search.o preprocess.o results.o

.PHONY: clean

search.a: $(DEPEN)
	ar rc search.a $(DEPEN)

clean:
	rm -f *~ \#*\# .*.swp *.o *.a
