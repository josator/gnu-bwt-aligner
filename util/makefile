COMMONS_DIR = ../commons/

#Debian distro location
#CUDA_LOCATION=/usr/lib/nvidia-cuda-toolkit

CC = gcc
NVCC = nvcc
CFLAGS = -Wall -m64 -D_LARGEFILE64_SOURCE=1 -D_FILE_OFFSET_BITS=64 -O9 -fPIC -std=c99 -DCSALIB_SEARCH #-I $(CUDA_LOCATION)/include/ -L /usr/local/cuda/lib64/ -lcudart -msse4.2 #-DVERBOSE_DBG
#NVCCFLAGS = --compiler-options -Wall,-fopenmp,-m64 -O3 -Xptxas -v -arch=sm_13 -DFM_COMP_32 -D_LARGEFILE64_SOURCE=1 -D_FILE_OFFSET_BITS=64 #-DVERBOSE_DBG

.PHONY: all clean

all: ../bin/preprocess_dbwt_new ../bin/exact_search

../bin/preprocess_dbwt_new: preprocess_dbwt_new.c dbwt csalib search commons
	$(CC) $(CFLAGS) -o $@ preprocess_dbwt_new.c ../dbwt/dbwt.a ../csalib/csalib.a ../search/search.a $(COMMONS_DIR)/commons.a

../bin/exact_search: exact_search.c csalib search commons
	$(CC) $(CFLAGS) -o $@ exact_search.c ../csalib/csalib.a ../search/search.a $(COMMONS_DIR)/commons.a

dbwt:
	make -C ../dbwt/ dbwt.a

csalib:
	make -C ../csalib/ csalib.a

search:
	make -C ../search/ search.a

commons:
	make -C $(COMMONS_DIR) commons.a

################################UTIL###############################

clean:
	rm -f *~ \#*\# .*.swp *.o *.a
