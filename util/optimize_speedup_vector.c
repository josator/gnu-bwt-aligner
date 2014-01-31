#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <cuda_runtime_api.h>

#include "../gpu/gpu.cuh"
#include "../search/search.h"
#include "../search/io.h"

#define TAM_BLOQUE_GPU 32

int MAXLINE;

float t_cpu=0, t_gpu=0;
int NUM_REP;
FILE *output_file;

int MAX_BUS_GPU;
int READS_PER_THREAD=0;

int main(int argc, char **argv) {
	char *h_Worig;
  uint8_t  *h_We, *d_We;
  uint64_t *h_nWe, *d_nWe;

  uint32_t *h_k, *h_l, *d_k, *d_l;
	intmax_t *h_k2, *h_l2;

	bwt_index backward;

  comp_matrix h_O, d_O;
  vector h_C, d_C, h_C1, d_C1;
  comp_vector h_S;

  cudaSetDevice(0);

  FILE *queries_file;

  cudaError_t error;

  if (argc!=8) {
    printf("Sintaxis:\n\t%s fichero_bus dir_entrada fichero_sal max_bus_gpu repeticiones max_length nucleotides\n", argv[0]);
    return 1;
  }

  timevars();
  init_replace_table(argv[7]);

  queries_file = fopen(argv[1], "r");
  if (!queries_file) {
    fprintf(stderr, "Error opening file: %s\n", argv[1]);
      return -1;
  }

  MAX_BUS_GPU = atoi(argv[4]);
	MAXLINE = atoi(argv[6]);

	tic("Cargando FM-Index");

  read_vector(&h_C, argv[2], "C");
  read_vector(&h_C1, argv[2], "C1");

  copy_vector_gpu(&d_C,   &h_C);
  copy_vector_gpu(&d_C1,  &h_C1);

  read_comp_matrix_gpu(&h_O, argv[2], "O");

  copy_comp_matrix_gpu(&d_O, &h_O);

  read_comp_vector(&h_S, argv[2], "S");

	backward.C  = h_C;
	backward.C1 = h_C1;
	backward.O  = h_O;
	backward.S  = h_S;

  toc();

  h_Worig = (char*)malloc(MAX_BUS_GPU * MAXLINE     * sizeof(char));
 
  cudaMallocHost((void**) &h_We, MAX_BUS_GPU * MAXLINE * sizeof(uint8_t));
  cudaMallocHost((void**) &h_nWe, MAX_BUS_GPU * sizeof(uint64_t));

  cudaMalloc((void**) &d_We,  MAX_BUS_GPU * MAXLINE * sizeof(uint8_t));
  manageCudaError();
  cudaMalloc((void**) &d_nWe, MAX_BUS_GPU * sizeof(uint64_t));
  manageCudaError();

  cudaMallocHost((void**) &h_k, MAX_BUS_GPU * MAXLINE * sizeof(uint32_t));
  cudaMallocHost((void**) &h_l, MAX_BUS_GPU * MAXLINE * sizeof(uint32_t));

  cudaMallocHost((void**) &h_k2, MAX_BUS_GPU * MAXLINE * sizeof(intmax_t));
  cudaMallocHost((void**) &h_l2, MAX_BUS_GPU * MAXLINE * sizeof(intmax_t));

  cudaMalloc((void**) &d_k, MAX_BUS_GPU * MAXLINE * sizeof(uint32_t));
  manageCudaError();
  cudaMalloc((void**) &d_l, MAX_BUS_GPU * MAXLINE * sizeof(uint32_t));
  manageCudaError();

  int TAM_BUS_GPU=0, NUM_BLOQUES_GPU=0;

  NUM_REP          = atoi(argv[5]);

  tic("Leer de disco");

  while(nextFASTAToken(queries_file, h_Worig + TAM_BUS_GPU * MAXLINE, h_We + TAM_BUS_GPU * MAXLINE, h_nWe + TAM_BUS_GPU)) {

    TAM_BUS_GPU++;

    if (TAM_BUS_GPU == MAX_BUS_GPU) break;

  }

  toc();

  NUM_BLOQUES_GPU = (TAM_BUS_GPU / TAM_BLOQUE_GPU);

  cudaThreadSynchronize();
  tic("CPU -> GPU");
  cudaMemcpy(d_We, h_We, TAM_BUS_GPU * MAXLINE * sizeof(uint8_t), cudaMemcpyHostToDevice);
  manageCudaError();
  cudaMemcpy(d_nWe,  h_nWe,  TAM_BUS_GPU * sizeof(uint64_t), cudaMemcpyHostToDevice);
  manageCudaError();
  cudaThreadSynchronize();
  toc();

	cudaThreadSynchronize();
  tic("GPU Kernel");
  BWExactSearchBackwardVectorGPUWrapper(NUM_BLOQUES_GPU, TAM_BLOQUE_GPU, d_We, d_nWe, MAXLINE, d_k, d_l, 0, d_O.siz-2, &d_C, &d_C1, &d_O);
  cudaThreadSynchronize();
  toc();

  cudaThreadSynchronize();
  tic("GPU -> CPU");
  cudaMemcpy(h_k, d_k, sizeof(uint32_t) * TAM_BUS_GPU * MAXLINE, cudaMemcpyDeviceToHost);
  manageCudaError();
  cudaMemcpy(h_l, d_l, sizeof(uint32_t) * TAM_BUS_GPU * MAXLINE, cudaMemcpyDeviceToHost);
  manageCudaError();  
  cudaThreadSynchronize();
  toc();

  tic("CPU kernel");
  for (int i=0; i<TAM_BUS_GPU; i++) {
    BWExactSearchVectorBackward(h_We + MAXLINE*i, 0, h_nWe[i]-1, 0, d_O.siz-2, h_k2 + MAXLINE*i, h_l2 + MAXLINE*i, &backward);
  }
  toc();

  /*
  for (int i=0; i<TAM_BUS_GPU; i++) {

    for (int j=0; j<h_nWe[i]; j++) {

      if (h_k[i*MAXLINE + j] != h_k2[i*MAXLINE + j]) {
	printf("Diferente %d %d\n", i, j);
	goto salir;
      }

    }

  }

  salir:
  */
  /*
  for (int i=0; i<h_nWe[0]; i++) {
    printf("%u ", h_k[i]);
  }
  printf("\n");

  printf("\n");

  for (int i=0; i<h_nWe[0]; i++) {
    printf("%u ", h_k2[i]);
  }
  printf("\n");
  */

  cudaFreeHost(h_k);
  cudaFree(d_k);
  cudaFreeHost(h_l);
  cudaFree(d_l);

  cudaFreeHost(h_We);
  cudaFreeHost(h_nWe);
  cudaFree(d_We);
  cudaFree(d_nWe);

  free(h_C.vector);
  cudaFree(d_C.vector);
  free(h_C1.vector);
  cudaFree(d_C1.vector);

  free_comp_matrix_gpu_host(NULL, &h_O);
  free_comp_matrix_gpu_device(NULL, &d_O);

  fclose(queries_file);

  return 0;

}
