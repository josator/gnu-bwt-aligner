#include "../BW_io.h"
#include "../commons/commons.h"
#include "../dbwt/dbwt.h"
#include "../csalib/csa.h"

int main(int argc, char **argv)
{

	ref_vector X;

	exome ex;

	check_syntax(argc, 4, "preprocess ref_file output_dir nucleotides");

	timevars();
  init_replace_table(argv[3]);

	//encode_reference(&X, &ex, argv[1]);
	//save_exome_file(&ex, argv[2]);

	tic("Calc. Burrows-Wheeler Transform -> Sadakane direct SAIS");
	//direct_bwt(X.vector, X.n, argv[2], NULL);
	toc();

	tic("Calc. Suffix Array -> CSALIB DNA");
	csa_new_from_bwt_gnu_bwt_wrapper(argv[2], NULL);
	toc();

	return 0;

}
