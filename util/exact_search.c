#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#include "../commons/commons.h"
#include "../search/search.h"
#include "../search/io.h"
#include "../csalib/csa.h"

int main(int argc, char **argv) {

	CSA SA;
	bwt_index backward;
	backward.csa = &SA;
	
	char *Worig;
  ref_vector W;

	results_list rl_prev, rl_next, rl_prev_i, rl_next_i, rl_final;
  uintmax_t read_index=0;

	uintmax_t RESULTS, FRAGSIZE;
  exome ex;

	FILE *queries_file, *output_file;
	char *fname[2];

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

	fname[0] = (char *) malloc(500 * sizeof(char));
	fname[0][0]='\0';
  strcat(fname[0], argv[2]);
  strcat(fname[0], "/backward.bwd");

	fname[1] = (char *) malloc(500 * sizeof(char));
	fname[1][0]='\0';
  strcat(fname[1], argv[2]);
  strcat(fname[1], "/backward.idx");

	csa_read(backward.csa, 2, fname);

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

  intmax_t *k = (intmax_t*)malloc(RESULTS * sizeof(intmax_t));
  intmax_t *l = (intmax_t*)malloc(RESULTS * sizeof(intmax_t));

  uintmax_t nW_aux;

  tic("Sequence mapping");

  while(nextFASTAToken(queries_file, Worig, W.vector, &nW_aux)) {

    if (!W.vector) {
      printf ("Error de lectura de la cadena a buscar\n");
      return -1;
    }

    W.n = nW_aux;

		rl_prev.read_index = read_index; rl_prev_i.read_index = read_index;
    rl_next.read_index = read_index; rl_next_i.read_index = read_index;
    rl_final.read_index = read_index;

		result r;
		init_result(&r, 0);
		change_result(&r, 0, size_SA(&backward)-1, W.n-1);
		bound_result(&r, 0, W.n-1);
		BWExactSearchBackward(W.vector, &backward, &r);
		//backward.csa->search(W.vector, W.n, backward.csa, &r_k, &r_l);
		printf("%ld - %ld : %ld\n", r.k, r.l);

		//BWSearchCPU(W.vector, W.n, &C, &C1, &O, &Oi, &S, &R, &Si, &Ri, &rl_prev, &rl_next, &rl_prev_i, &rl_next_i, &rl_final, FRAGSIZE, 1);
    //write_results(&rl_final, k, l, &ex, &S, &Si, &C, &O, &Oi, Worig, nW_aux, 1, output_file);

		read_index++;

	}

	toc();

	printf("Memory frees\n");

	free(W.vector);
  free(Worig);

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
