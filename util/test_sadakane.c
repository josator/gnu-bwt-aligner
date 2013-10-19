#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sys/time.h>

#include "BW_csafm.h"

int found;
unsigned int enW;
char search[MAXLINE];

unsigned int read_index=0;

int main(int argc, char **argv) {

  ref_vector B1, B2;

  check_syntax(argc, 2, "test_sadakane input_dir");

	read_ref_vector(&B1, argv[1], "B");
	read_ref_vector(&B2, argv[1], "B2");

  for (unsigned int i=0; i<B1.n; i++) {

    if (B1.vector[i] != B2.vector[i]) {
      printf("%u -> %u =! %u\n", i, B1.vector[i], B2.vector[i]);
    }

	}

	printf("Memory frees\n");

	free(B1.vector);
  free(B2.vector);

	return 0;

}
