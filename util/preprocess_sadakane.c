#include "commons/commons.h"
#include "commons/string_utils.h"

#include "BW_preprocess.h"

#include "dbwt/dbwt.h"

int main(int argc, char **argv)
{

  ref_vector X;

	check_syntax(argc, 5, "preprocess ref_file output_dir s_ratio nucleotides");

  timevars();
  init_replace_table(argv[4]);

  encode_reference(&X, NULL, false, argv[1]);

	tic("Calc. Suffix Array -> Sadakane direct SAIS");
  bwt(X.vector, X.n-1);
	toc();

}
