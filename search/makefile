COMMONS_DIR=../commons/

#Debian distro location
#CUDA_LOCATION=/usr/lib/nvidia-cuda-toolkit

CC = gcc
#NVCC = nvcc
CFLAGS = -Wall -O3 -std=c99 -fopenmp -m64 -DCSALIB_SEARCH -DCSALIB_CSA -DFM_COMP_32 -D_LARGEFILE64_SOURCE=1 -D_FILE_OFFSET_BITS=64 #-I $(CUDA_LOCATION)/include/ -L /usr/local/cuda/lib64/ -lcudart -msse4.2 #-DVERBOSE_DBG
#NVCCFLAGS = --compiler-options -Wall,-fopenmp,-m64 -O3 -Xptxas -v -arch=sm_13 -DFM_COMP_32 -D_LARGEFILE64_SOURCE=1 -D_FILE_OFFSET_BITS=64 #-DVERBOSE_DBG

DEPEN = io.o csafm.o search.o preprocess.o results.o

.PHONY: all clean depen

depen: $(DEPEN)

io.o: io.c io.h runtime.h $(COMMONS_DIR)/commons.h $(COMMONS_DIR)/string_utils.h
	$(CC) $(CFLAGS) -c io.c

csafm.o: csafm.c csafm.h types.h $(COMMONS_DIR)/commons.h $(COMMONS_DIR)/string_utils.h
	$(CC) $(CFLAGS) -c csafm.c

search.o: search.c search.h runtime.h results.h $(COMMONS_DIR)/commons.h $(COMMONS_DIR)/string_utils.h
	$(CC) $(CFLAGS) -c search.c

preprocess.o: preprocess.c preprocess.h types.h $(COMMONS_DIR)/commons.h $(COMMONS_DIR)/string_utils.h
	$(CC) $(CFLAGS) -c preprocess.c

results.o: results.c results.h runtime.h $(COMMONS_DIR)/commons.h $(COMMONS_DIR)/string_utils.h
	$(CC) $(CFLAGS) -c results.c

clean:
	rm -f *~ \#*\# .*.swp *.o *.a
