TARGETS = bin/preprocess_dbwt_new bin/preprocess_dbwt_old #bin/test_sadakane bin/inexact_search #bin/generateSRcomp bin/optimize_speedup_errors

COMMONS_DIR = commons/

#Debian distro location
#CUDA_LOCATION=/usr/lib/nvidia-cuda-toolkit

CC = gcc
NVCC = nvcc
CFLAGS = -Wall -O3 -std=c99 -fopenmp -m64 -DFM_COMP_32 -D_LARGEFILE64_SOURCE=1 -D_FILE_OFFSET_BITS=64 #-I $(CUDA_LOCATION)/include/ -L /usr/local/cuda/lib64/ -lcudart -msse4.2 #-DVERBOSE_DBG
#NVCCFLAGS = --compiler-options -Wall,-fopenmp,-m64 -O3 -Xptxas -v -arch=sm_13 -DFM_COMP_32 -D_LARGEFILE64_SOURCE=1 -D_FILE_OFFSET_BITS=64 #-DVERBOSE_DBG

COMMONS_OBJECTS = $(COMMONS_DIR)/commons.o $(COMMONS_DIR)/string_utils.o
DBWT_OBJECTS = dbwt/dbwt.o dbwt/sais.o dbwt/queue.o dbwt/utils.o
CSALIB_OBJECTS = csalib/csa.o csalib/mmap.o csalib/psi1.o csalib/diskbuf.o csalib/lf_dna.o csalib/lf_wt.o csalib/densearray.o csalib/huffman.o csalib/comparray.o csalib/psi2.o csalib/densearray2.o csalib/sparsearray2.o csalib/sparsearray.o csalib/cst.o csalib/lf_bit.o

.PHONY: all clean commons dbwt csalib

all: dbwt csalib $(TARGETS)

bin/preprocess_dbwt_old: $(COMMONS_OBJECTS) $(DBWT_OBJECTS) preprocess_dbwt_old.o preprocess.o csafm.o io.o search.o
	$(CC) $(CFLAGS) $(COMMONS_OBJECTS) $(DBWT_OBJECTS) preprocess.o csafm.o io.o search.o preprocess_dbwt_old.o -o bin/preprocess_dbwt_old

preprocess_dbwt_old.o: util/preprocess_dbwt_old.c types.h search.h preprocess.h csafm.h $(COMMONS_DIR)/commons.h $(COMMONS_DIR)/string_utils.h dbwt/dbwt.h
	$(CC) $(CFLAGS) -c util/preprocess_dbwt_old.c

bin/preprocess_dbwt_new: $(COMMONS_OBJECTS) $(DBWT_OBJECTS) $(CSALIB_OBJECTS) preprocess_dbwt_new.o preprocess.o csafm.o io.o search.o
	$(CC) $(CFLAGS) $(COMMONS_OBJECTS) $(DBWT_OBJECTS) $(CSALIB_OBJECTS) preprocess.o csafm.o io.o search.o preprocess_dbwt_new.o -o bin/preprocess_dbwt_new

preprocess_dbwt_new.o: util/preprocess_dbwt_new.c types.h search.h preprocess.h csafm.h $(COMMONS_DIR)/commons.h $(COMMONS_DIR)/string_utils.h dbwt/dbwt.h csalib/csa.h
	$(CC) $(CFLAGS) -c util/preprocess_dbwt_new.c

bin/test_sadakane: $(COMMONS_OBJECTS) $(DBWT_OBJECTS) test_sadakane.o preprocess.o csafm.o
	$(CC) $(CFLAGS) $(COMMONS_OBJECTS) $(DBWT_OBJECTS) preprocess.o csafm.o test_sadakane.o -o bin/test_sadakane

test_sadakane.o: util/test_sadakane.c types.h preprocess.h csafm.h $(COMMONS_DIR)/commons.h $(COMMONS_DIR)/string_utils.h
	$(CC) $(CFLAGS) -c util/test_sadakane.c

bin/search: search.o results.o csafm.o search.o io.o $(COMMONS_DIR)/string_utils.o
	$(CC) $(CFLAGS) results.o csafm.o search.o io.o search.o $(COMMONS_DIR)/string_utils.o -o bin/search

search.o: util/search.c gpu.cuh csafm.h types.h $(COMMONS_DIR)/commons.h $(COMMONS_DIR)/string_utils.h
	$(CC) $(CFLAGS) -c util/search.c

bin/inexact_search: inexact_search.o results.o csafm.o io.o search.o $(COMMONS_DIR)/string_utils.o
	$(CC) $(CFLAGS) results.o csafm.o io.o search.o inexact_search.o $(COMMONS_DIR)/string_utils.o -o bin/inexact_search

inexact_search.o: util/inexact_search.c gpu.cuh csafm.h types.h $(COMMONS_DIR)/commons.h $(COMMONS_DIR)/string_utils.h
	$(CC) $(CFLAGS) -c util/inexact_search.c

bin/generateSRcomp: generateSRcomp.o preprocess.o csafm.o results.o $(COMMONS_DIR)/string_utils.o
	$(CC) $(CFLAGS) preprocess.o csafm.o results.o generateSRcomp.o $(COMMONS_DIR)/string_utils.o -o bin/generateSRcomp

generateSRcomp.o: util/generateSRcomp.c preprocess.h csafm.h results.h types.h $(COMMONS_DIR)/commons.h $(COMMONS_DIR)/string_utils.h
	$(CC) $(CFLAGS) -c util/generateSRcomp.c


##########################OBJECTS##################################

io.o: io.c io.h types.h $(COMMONS_DIR)/commons.h $(COMMONS_DIR)/string_utils.h
	$(CC) $(CFLAGS) -c io.c

csafm.o: csafm.c csafm.h types.h $(COMMONS_DIR)/commons.h $(COMMONS_DIR)/string_utils.h
	$(CC) $(CFLAGS) -c csafm.c

search.o: search.c search.h types.h results.h $(COMMONS_DIR)/commons.h $(COMMONS_DIR)/string_utils.h
	$(CC) $(CFLAGS) -c search.c

preprocess.o: preprocess.c preprocess.h types.h $(COMMONS_DIR)/commons.h $(COMMONS_DIR)/string_utils.h
	$(CC) $(CFLAGS) -c preprocess.c

results.o: results.c results.h types.h $(COMMONS_DIR)/commons.h $(COMMONS_DIR)/string_utils.h
	$(CC) $(CFLAGS) -c results.c

string_utils.o: $(COMMONS_DIR)/types.h $(COMMONS_DIR)/string_utils.c $(COMMONS_DIR)/string_utils.h
	$(CC) $(CFLAGS) -c $(COMMONS_DIR)/string_utils.c

commons:
	make -C $(COMMONS_DIR) depen

dbwt:
	make -C dbwt/ depen

csalib:
	make -C csalib/ depen

################################UTIL###############################

clean:
	rm -f *~ \#*\# .*.swp *.o
	make -C $(COMMONS_DIR) clean
	make -C dbwt/ clean
	make -C csalib/ clean
