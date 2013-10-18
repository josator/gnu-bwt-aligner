UTIL_DIR = util
BIN_DIR = bin

ifdef LIBS_HOME
	LIBS_ROOT = $(LIBS_HOME)
else
	LIBS_ROOT = .
endif

COMMONS_DIR = $(LIBS_ROOT)/commons

#Debian distro location
CUDA_LOCATION=/usr/lib/nvidia-cuda-toolkit

CC = g++
NVCC = nvcc
CFLAGS = -Wall -g -fopenmp -DFM_COMP_32 -D_LARGEFILE64_SOURCE=1 -D_FILE_OFFSET_BITS=64 -I . -I $(LIBS_ROOT)/common-libs/ -I $(CUDA_LOCATION)/include/ -L /usr/local/cuda/lib64/ -lcudart #-DVERBOSE_DBG
NVCCFLAGS = --compiler-options -Wall,-fopenmp -O3 -Xptxas -v -arch=sm_13 -DFM_COMP_32 -D_LARGEFILE64_SOURCE=1 -D_FILE_OFFSET_BITS=64 -I . -I $(LIBS_ROOT)/common-libs/ #-DVERBOSE_DBG
#THRUST_FLAGS = -Xcompiler -fopenmp -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_OMP -lgomp

TARGETS = $(BIN_DIR)/preprocess_sadakane #$(BIN_DIR)/preprocess $(BIN_DIR)/inexact_search #$(BIN_DIR)/search_gpu #$(BIN_DIR)/generateSRcomp $(BIN_DIR)/optimize_speedup_errors

DBWT_OBJECTS = dbwt.o sais.o queue.o utils.o

.PHONY: all clean dbwt

all: dbwt $(TARGETS)

$(BIN_DIR)/preprocess: $(DBWT_OBJECTS) preprocess.o BW_preprocess.o BW_csafm.o $(COMMONS_DIR)/string_utils.o
	$(CC) $(CFLAGS) $(DBWT_OBJECTS) BW_preprocess.o BW_csafm.o preprocess.o $(COMMONS_DIR)/string_utils.o -o $(BIN_DIR)/preprocess

preprocess.o: $(UTIL_DIR)/preprocess.c $(COMMONS_DIR)/BW_types.h BW_preprocess.h BW_csafm.h $(COMMONS_DIR)/commons.h $(COMMONS_DIR)/string_utils.h
	$(CC) $(CFLAGS) -c $(UTIL_DIR)/preprocess.c

$(BIN_DIR)/preprocess_sadakane: $(DBWT_OBJECTS) preprocess_sadakane.o BW_preprocess.o BW_csafm.o $(COMMONS_DIR)/string_utils.o
	$(CC) $(CFLAGS) $(DBWT_OBJECTS) BW_preprocess.o BW_csafm.o preprocess_sadakane.o $(COMMONS_DIR)/string_utils.o -o $(BIN_DIR)/preprocess_sadakane

preprocess_sadakane.o: $(UTIL_DIR)/preprocess_sadakane.c $(COMMONS_DIR)/BW_types.h BW_preprocess.h BW_csafm.h $(COMMONS_DIR)/commons.h $(COMMONS_DIR)/string_utils.h
	$(CC) $(CFLAGS) -c $(UTIL_DIR)/preprocess_sadakane.c

$(BIN_DIR)/search: search.o BW_results.o BW_csafm.o BW_search.o $(COMMONS_DIR)/string_utils.o
	$(CC) $(CFLAGS) BW_results.o BW_csafm.o BW_search.o search.o $(COMMONS_DIR)/string_utils.o -o $(BIN_DIR)/search

search.o: $(UTIL_DIR)/search.c BW_gpu.cuh BW_csafm.h $(COMMONS_DIR)/BW_types.h $(COMMONS_DIR)/commons.h $(COMMONS_DIR)/string_utils.h
	$(CC) $(CFLAGS) -c $(UTIL_DIR)/search.c

$(BIN_DIR)/inexact_search: inexact_search.o BW_results.o BW_csafm.o BW_search.o $(COMMONS_DIR)/string_utils.o
	$(CC) $(CFLAGS) BW_results.o BW_csafm.o BW_search.o inexact_search.o $(COMMONS_DIR)/string_utils.o -o $(BIN_DIR)/inexact_search

inexact_search.o: $(UTIL_DIR)/inexact_search.c BW_gpu.cuh BW_csafm.h $(COMMONS_DIR)/BW_types.h $(COMMONS_DIR)/commons.h $(COMMONS_DIR)/string_utils.h
	$(CC) $(CFLAGS) -c $(UTIL_DIR)/inexact_search.c

$(BIN_DIR)/search_gpu: search_gpu.o BW_csafm.o BW_results.o BW_search.o BW_gpu.o $(COMMONS_DIR)/string_utils.o
	$(CC) $(CFLAGS) BW_csafm.o BW_results.o BW_search.o BW_gpu.o $(COMMONS_DIR)/string_utils.o search_gpu.o -o $(BIN_DIR)/search_gpu

search_gpu.o: $(UTIL_DIR)/search_gpu.c BW_gpu.cuh BW_csafm.h BW_results.h $(COMMONS_DIR)/BW_types.h $(COMMONS_DIR)/commons.h $(COMMONS_DIR)/string_utils.h
	$(CC) $(CFLAGS) -c $(UTIL_DIR)/search_gpu.c

$(BIN_DIR)/generateSRcomp: generateSRcomp.o BW_preprocess.o BW_csafm.o BW_results.o $(COMMONS_DIR)/string_utils.o
	$(CC) $(CFLAGS) BW_preprocess.o BW_csafm.o BW_results.o generateSRcomp.o $(COMMONS_DIR)/string_utils.o -o $(BIN_DIR)/generateSRcomp

generateSRcomp.o: $(UTIL_DIR)/generateSRcomp.c BW_preprocess.h BW_csafm.h BW_results.h $(COMMONS_DIR)/BW_types.h $(COMMONS_DIR)/commons.h $(COMMONS_DIR)/string_utils.h
	$(CC) $(CFLAGS) -c $(UTIL_DIR)/generateSRcomp.c

$(BIN_DIR)/optimize_speedup_errors: optimize_speedup_errors.o BW_search.o BW_csafm.o BW_results.o BW_gpu.o $(COMMONS_DIR)/string_utils.o
	$(CC) $(CFLAGS) BW_gpu.o BW_search.o BW_csafm.o BW_results.o $(COMMONS_DIR)/string_utils.o optimize_speedup_errors.o -o $(BIN_DIR)/optimize_speedup_errors

optimize_speedup_errors.o: $(UTIL_DIR)/optimize_speedup_errors.c $(COMMONS_DIR)/BW_types.h BW_csafm.h BW_results.h BW_search.h $(COMMONS_DIR)/commons.h $(COMMONS_DIR)/string_utils.h
	$(CC) $(CFLAGS) -c $(UTIL_DIR)/optimize_speedup_errors.c

##########################OBJECTS##################################

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

dbwt:
		make -C dbwt/ depen

################################UTIL###############################

clean:
	rm -f *~ \#*\# .*.swp *.o
	make -C dbwt/ clean
