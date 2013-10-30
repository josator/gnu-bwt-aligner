TARGETS = bin/preprocess_dbwt_new bin/preprocess_dbwt_old #bin/test_sadakane bin/inexact_search bin/search_gpu #bin/generateSRcomp bin/optimize_speedup_errors

COMMONS_DIR = commons/

#Debian distro location
CUDA_LOCATION=/usr/lib/nvidia-cuda-toolkit

CC = g++
NVCC = nvcc
CFLAGS = -Wall -g -fopenmp -m64 -DFM_COMP_32 -D_LARGEFILE64_SOURCE=1 -D_FILE_OFFSET_BITS=64 -I $(CUDA_LOCATION)/include/ -L /usr/local/cuda/lib64/ -lcudart #-msse4.2 #-DVERBOSE_DBG
NVCCFLAGS = --compiler-options -Wall,-fopenmp,-m64 -O3 -Xptxas -v -arch=sm_13 -DFM_COMP_32 -D_LARGEFILE64_SOURCE=1 -D_FILE_OFFSET_BITS=64 #-DVERBOSE_DBG
#THRUST_FLAGS = -Xcompiler -fopenmp -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_OMP -lgomp

COMMONS_OBJECTS = $(COMMONS_DIR)/commons.o $(COMMONS_DIR)/string_utils.o
DBWT_OBJECTS = dbwt/dbwt.o dbwt/sais.o dbwt/queue.o dbwt/utils.o
CSALIB_OBJECTS = csalib/csa.o csalib/mmap.o csalib/psi1.o csalib/diskbuf.o csalib/lf_dna.o csalib/lf_wt.o csalib/densearray.o csalib/huffman.o csalib/comparray.o csalib/psi2.o csalib/densearray2.o csalib/sparsearray2.o csalib/sparsearray.o csalib/cst.o csalib/lf_bit.o

.PHONY: all clean commons dbwt csalib

all: dbwt csalib $(TARGETS)

bin/preprocess_dbwt_old: $(COMMONS_OBJECTS) $(DBWT_OBJECTS) preprocess_dbwt_old.o BW_preprocess.o BW_csafm.o BW_io.o BW_search.o
	$(CC) $(CFLAGS) $(COMMONS_OBJECTS) $(DBWT_OBJECTS) BW_preprocess.o BW_csafm.o BW_io.o BW_search.o preprocess_dbwt_old.o -o bin/preprocess_dbwt_old

preprocess_dbwt_old.o: util/preprocess_dbwt_old.c $(COMMONS_DIR)/BW_types.h BW_search.h BW_preprocess.h BW_csafm.h $(COMMONS_DIR)/commons.h $(COMMONS_DIR)/string_utils.h dbwt/dbwt.h
	$(CC) $(CFLAGS) -c util/preprocess_dbwt_old.c

bin/preprocess_dbwt_new: $(COMMONS_OBJECTS) $(DBWT_OBJECTS) $(CSALIB_OBJECTS) preprocess_dbwt_new.o BW_preprocess.o BW_csafm.o BW_io.o BW_search.o
	$(CC) $(CFLAGS) $(COMMONS_OBJECTS) $(DBWT_OBJECTS) $(CSALIB_OBJECTS) BW_preprocess.o BW_csafm.o BW_io.o BW_search.o preprocess_dbwt_new.o -o bin/preprocess_dbwt_new

preprocess_dbwt_new.o: util/preprocess_dbwt_new.c $(COMMONS_DIR)/BW_types.h BW_search.h BW_preprocess.h BW_csafm.h $(COMMONS_DIR)/commons.h $(COMMONS_DIR)/string_utils.h dbwt/dbwt.h csalib/csa.h
	$(CC) $(CFLAGS) -c util/preprocess_dbwt_new.c

bin/test_sadakane: $(COMMONS_OBJECTS) $(DBWT_OBJECTS) test_sadakane.o BW_preprocess.o BW_csafm.o
	$(CC) $(CFLAGS) $(COMMONS_OBJECTS) $(DBWT_OBJECTS) BW_preprocess.o BW_csafm.o test_sadakane.o -o bin/test_sadakane

test_sadakane.o: util/test_sadakane.c $(COMMONS_DIR)/BW_types.h BW_preprocess.h BW_csafm.h $(COMMONS_DIR)/commons.h $(COMMONS_DIR)/string_utils.h
	$(CC) $(CFLAGS) -c util/test_sadakane.c

bin/search: search.o BW_results.o BW_csafm.o BW_search.o BW_io.o $(COMMONS_DIR)/string_utils.o
	$(CC) $(CFLAGS) BW_results.o BW_csafm.o BW_search.o BW_io.o search.o $(COMMONS_DIR)/string_utils.o -o bin/search

search.o: util/search.c BW_gpu.cuh BW_csafm.h $(COMMONS_DIR)/BW_types.h $(COMMONS_DIR)/commons.h $(COMMONS_DIR)/string_utils.h
	$(CC) $(CFLAGS) -c util/search.c

bin/inexact_search: inexact_search.o BW_results.o BW_csafm.o BW_io.o BW_search.o $(COMMONS_DIR)/string_utils.o
	$(CC) $(CFLAGS) BW_results.o BW_csafm.o BW_io.o BW_search.o inexact_search.o $(COMMONS_DIR)/string_utils.o -o bin/inexact_search

inexact_search.o: util/inexact_search.c BW_gpu.cuh BW_csafm.h $(COMMONS_DIR)/BW_types.h $(COMMONS_DIR)/commons.h $(COMMONS_DIR)/string_utils.h
	$(CC) $(CFLAGS) -c util/inexact_search.c

bin/search_gpu: search_gpu.o BW_csafm.o BW_results.o BW_io.o BW_search.o BW_gpu.o $(COMMONS_DIR)/string_utils.o
	$(CC) $(CFLAGS) BW_csafm.o BW_results.o BW_io.o BW_search.o BW_gpu.o $(COMMONS_DIR)/string_utils.o search_gpu.o -o bin/search_gpu

search_gpu.o: util/search_gpu.c BW_gpu.cuh BW_csafm.h BW_results.h $(COMMONS_DIR)/BW_types.h $(COMMONS_DIR)/commons.h $(COMMONS_DIR)/string_utils.h
	$(CC) $(CFLAGS) -c util/search_gpu.c

bin/generateSRcomp: generateSRcomp.o BW_preprocess.o BW_csafm.o BW_results.o $(COMMONS_DIR)/string_utils.o
	$(CC) $(CFLAGS) BW_preprocess.o BW_csafm.o BW_results.o generateSRcomp.o $(COMMONS_DIR)/string_utils.o -o bin/generateSRcomp

generateSRcomp.o: util/generateSRcomp.c BW_preprocess.h BW_csafm.h BW_results.h $(COMMONS_DIR)/BW_types.h $(COMMONS_DIR)/commons.h $(COMMONS_DIR)/string_utils.h
	$(CC) $(CFLAGS) -c util/generateSRcomp.c

bin/optimize_speedup_errors: optimize_speedup_errors.o BW_search.o BW_csafm.o BW_results.o BW_gpu.o $(COMMONS_DIR)/string_utils.o
	$(CC) $(CFLAGS) BW_gpu.o BW_search.o BW_csafm.o BW_results.o $(COMMONS_DIR)/string_utils.o optimize_speedup_errors.o -o bin/optimize_speedup_errors

optimize_speedup_errors.o: util/optimize_speedup_errors.c $(COMMONS_DIR)/BW_types.h BW_csafm.h BW_results.h BW_search.h $(COMMONS_DIR)/commons.h $(COMMONS_DIR)/string_utils.h
	$(CC) $(CFLAGS) -c util/optimize_speedup_errors.c

##########################OBJECTS##################################

BW_io.o: BW_io.c BW_io.h $(COMMONS_DIR)/BW_types.h $(COMMONS_DIR)/commons.h $(COMMONS_DIR)/string_utils.h
	$(CC) $(CFLAGS) -c BW_io.c

BW_csafm.o: BW_csafm.c BW_csafm.h $(COMMONS_DIR)/BW_types.h $(COMMONS_DIR)/commons.h $(COMMONS_DIR)/string_utils.h
	$(CC) $(CFLAGS) -c BW_csafm.c

BW_search.o: BW_search.c BW_search.h $(COMMONS_DIR)/BW_types.h BW_results.h $(COMMONS_DIR)/commons.h $(COMMONS_DIR)/string_utils.h
	$(CC) $(CFLAGS) -c BW_search.c

BW_preprocess.o: BW_preprocess.c BW_preprocess.h $(COMMONS_DIR)/BW_types.h $(COMMONS_DIR)/commons.h $(COMMONS_DIR)/string_utils.h
	$(CC) $(CFLAGS) -c BW_preprocess.c

BW_gpu.o: BW_gpu.cu BW_gpu.cuh BW_results.h $(COMMONS_DIR)/BW_types.h $(COMMONS_DIR)/commons.h $(COMMONS_DIR)/string_utils.h
	$(NVCC) $(NVCCFLAGS) -c BW_gpu.cu

BW_results.o: BW_results.c BW_results.h $(COMMONS_DIR)/BW_types.h $(COMMONS_DIR)/commons.h $(COMMONS_DIR)/string_utils.h
	$(CC) $(CFLAGS) -c BW_results.c

string_utils.o: $(COMMONS_DIR)/BW_types.h $(COMMONS_DIR)/string_utils.c $(COMMONS_DIR)/string_utils.h
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
