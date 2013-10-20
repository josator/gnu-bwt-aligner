#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sys/time.h>

#include "commons/commons.h"
#include "commons/string_utils.h"

#include "BW_csafm.h"

int main(int argc, char **argv) {

  ref_vector B1, B2;

  check_syntax(argc, 3, "test_sadakane input_dir nucleotides");

  init_replace_table(argv[3]);

	read_ref_vector(&B1, argv[1], "B");
	read_ref_vector(&B2, argv[1], "B2");

  for (unsigned int i=0; i<B1.n; i++) {

//    if (B1.vector[i] != B2.vector[i]) {
      printf("%u\t-> %c - %c\n", i, rev_table[B1.vector[i]],rev_table[ B2.vector[i]]);
//    }

	}

	printf("Memory frees\n");

	free(B1.vector);
  free(B2.vector);

	return 0;

}
