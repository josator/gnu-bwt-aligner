COMMONS_DIR = ../commons/

#CUDA location
CUDA_LOCATION=/usr/local/cuda

CC = gcc
CFLAGS = -Wall -m64 -D_LARGEFILE64_SOURCE=1 -D_FILE_OFFSET_BITS=64 -O9 -fPIC -fopenmp -std=c99 -DFM_COMP_64 #-DCSALIB_SEARCH -DSA_64 -msse4.2 -DVERBOSE_DBG
CFLAGSCUDA = -I $(CUDA_LOCATION)/include/ -L $(CUDA_LOCATION)/lib64/ -lcudart

.PHONY: all gpu-tools dbwt csalib search gpu commons clean

all: ../bin/preprocess_dbwt_old ../bin/preprocess_dbwt_new ../bin/exact_search ../bin/inexact_search

gpu-tools: ../bin/search_gpu ../bin/optimize_speedup_vector

../bin/preprocess_dbwt_old: preprocess_dbwt_old.c dbwt csalib search commons
	$(CC) $(CFLAGS) -o $@ preprocess_dbwt_old.c ../search/search.a ../dbwt/dbwt.a ../csalib/csalib.a $(COMMONS_DIR)/commons.a

../bin/preprocess_dbwt_new: preprocess_dbwt_new.c dbwt csalib search commons
	$(CC) $(CFLAGS) -o $@ preprocess_dbwt_new.c ../search/search.a ../dbwt/dbwt.a ../csalib/csalib.a $(COMMONS_DIR)/commons.a

../bin/exact_search: exact_search.c csalib search commons
	$(CC) $(CFLAGS) -o $@ exact_search.c ../csalib/csalib.a ../search/search.a $(COMMONS_DIR)/commons.a

../bin/inexact_search: inexact_search.c csalib search commons
	$(CC) $(CFLAGS) -o $@ inexact_search.c ../csalib/csalib.a ../search/search.a $(COMMONS_DIR)/commons.a

../bin/search_gpu: search_gpu.c csalib search gpu commons
	$(CC) $(CFLAGS) -o $@ search_gpu.c ../csalib/csalib.a ../search/search.a ../gpu/gpu.a $(COMMONS_DIR)/commons.a $(CFLAGSCUDA)

../bin/optimize_speedup_vector: optimize_speedup_vector.c csalib search gpu commons
	$(CC) $(CFLAGS) -o $@ optimize_speedup_vector.c ../csalib/csalib.a ../search/search.a ../gpu/gpu.a $(COMMONS_DIR)/commons.a $(CFLAGSCUDA)

dbwt:
	make -C ../dbwt/ dbwt.a

csalib:
	make -C ../csalib/ csalib.a

search:
	make -C ../search/ search.a

gpu:
	make -C ../gpu/ gpu.a

commons:
	make -C $(COMMONS_DIR) commons.a

################################UTIL###############################

clean:
	rm -f *~ \#*\# .*.swp *.o *.a
