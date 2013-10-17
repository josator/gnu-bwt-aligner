#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#include "commons/commons.h"
#include "BW_search.h"
#include "BW_csafm.h"
#include "BW_results.h"

int main(int argc, char **argv) {

	char *Worig;
  byte_vector W;
  vector C, C1;
  vector rC, rC1;
  comp_matrix O, Oi;
  comp_matrix rO, rOi;
  comp_vector S, Si;
  comp_vector R, Ri;

  results_list rl_prev, rl_next, rl_prev_i, rl_next_i, rl_final;
  uintmax_t read_index=0;

  VECTOR_TYPE RESULTS, FRAGSIZE;
  exome ex;

  FILE *queries_file, *output_file;

  if (argc!=7) {
    printf("Syntax:\n\t%s search_file input_dir output_file results_buffer frag_size nucleotides\n", argv[0]);
    return 1;
  }

  timevars();
  init_replace_table(argv[6]);

  queries_file = fopen(argv[1], "r");
  check_file_open(queries_file, argv[1]);
  output_file = fopen(argv[3], "w");
  check_file_open(output_file, argv[3]);

  RESULTS  = atoi(argv[4]);
  FRAGSIZE = atoi(argv[5]);

  if (FRAGSIZE <= 0) {
    fprintf(stderr, "Fragsize must be greater than 0\n");
    exit(1);
  }

  tic("Loading FM-Index");

  read_vector(&C, argv[2], "C");
  read_vector(&C1, argv[2], "C1");
  read_comp_matrix(&O, argv[2], "O");
  read_comp_matrix(&Oi, argv[2], "Oi");

  read_comp_vector(&S, argv[2], "Scomp");
  read_comp_vector(&Si, argv[2], "Scompi");
  read_comp_vector(&R, argv[2], "Rcomp");
  read_comp_vector(&Ri, argv[2], "Rcompi");

  reverse_strand_C(&rC, &C, &rC1, &C1);
  reverse_strand_O(&rO, &O);
  reverse_strand_O(&rOi, &Oi);

  toc();

  tic("Preparing search space");

  load_exome_file(&ex, argv[2]);
  new_results_list(&rl_prev, RESULTS); new_results_list(&rl_prev_i, RESULTS);
  new_results_list(&rl_next, RESULTS); new_results_list(&rl_next_i, RESULTS);
  new_results_list(&rl_final, RESULTS);

  Worig = (char *) malloc( MAXLINE * sizeof(char) );
  check_malloc(Worig, "main");
  W.vector = (uint8_t *) malloc( MAXLINE * sizeof(uint8_t) );
  check_malloc(W.vector, "main");

  toc();

  VECTOR_TYPE *k = (VECTOR_TYPE*)malloc(RESULTS * sizeof(VECTOR_TYPE));
  VECTOR_TYPE *l = (VECTOR_TYPE*)malloc(RESULTS * sizeof(VECTOR_TYPE));

  VECTOR_TYPE nW_aux;

  tic("Sequence mapping");

  while(nextFASTAToken(queries_file, Worig, W.vector, &nW_aux, NULL, NULL)) {

    if (!W.vector) {
      printf ("Error de lectura de la cadena a buscar\n");
      return -1;
    }

    W.n = nW_aux;

    rl_prev.read_index = read_index; rl_prev_i.read_index = read_index;
    rl_next.read_index = read_index; rl_next_i.read_index = read_index;
    rl_final.read_index = read_index;

    BWSearchCPU(W.vector, W.n, &C, &C1, &O, &Oi, &S, &R, &Si, &Ri, &rl_prev, &rl_next, &rl_prev_i, &rl_next_i, &rl_final, FRAGSIZE, 1);
    write_results(&rl_final, k, l, &ex, &S, &Si, &C, &O, &Oi, Worig, nW_aux, 1, output_file);

    read_index++;

  }

  toc();
  
  printf("Memory frees\n");

  free(W.vector);
  free(Worig);

  free(C.vector);
  free(C1.vector);
  free(rC.vector);
  free(rC1.vector);

  free_comp_matrix(&rO,&O);
  free_comp_matrix(&rOi,&Oi);

  free(S.vector);
  free(Si.vector);
  free(R.vector);
  free(Ri.vector);

	free(k);
  free(l);

  free(rl_prev.list);
  free(rl_next.list);
  free(rl_prev_i.list);
  free(rl_next_i.list);
  free(rl_final.list);

  fclose(queries_file);
  fclose(output_file);

  return 0;

}
