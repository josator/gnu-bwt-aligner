#include "commons/commons.h"
#include "commons/string_utils.h"

#include "BW_preprocess.h"

int main(int argc, char **argv)
{

  ref_vector X, B, Bi;
  vector C, C1;
  comp_vector S, Si, Scomp, Scompi;
  comp_vector R, Ri, Rcomp, Rcompi;
  comp_matrix O, Oi;

  int s_ratio;

  exome ex;

	check_syntax(argc, 5, "preprocess ref_file output_dir s_ratio nucleotides");

  timevars();
  init_replace_table(argv[4]);

  s_ratio = atoi(argv[3]);

  encode_reference(&X, &ex, argv[1]);
  save_exome_file(&ex, argv[2]);
  save_ref_vector(&X, argv[2], "X");
  revstring((char *)X.vector, X.n-1);
  save_ref_vector(&X, argv[2], "Xi");
  revstring((char *)X.vector, X.n-1);

	tic("Calc. Suffix Array -> Sadakane direct SAIS");
  calculate_B_sadakane_SAIS(&B, &X);
	toc();

}
