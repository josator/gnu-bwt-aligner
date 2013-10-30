#include "../BW_preprocess.h"
#include "../BW_csafm.h"
#include "../BW_io.h"

int main(int argc, char **argv)
{

	ref_vector X, Xi, B, Bi;
	vector C, C1;
	comp_matrix O, Oi;
	comp_vector S, R, Si, Ri;

	exome ex;

	check_syntax(argc, 5, "preprocess ref_file output_dir s_ratio nucleotides");

	timevars();
  init_replace_table(argv[4]);

	SA_TYPE ratio = atoi(argv[3]);

	encode_reference(&X, &ex, false, argv[1]);
	save_ref_vector(&X, argv[2], "X");
	save_exome_file(&ex, false, argv[2]);
  print_vector(X.vector, X.n);

	tic("Calc. Suffix Array -> Sadakane direct SAIS");
	calculate_and_save_B(&X, argv[2], "B");
	toc();

	read_ref_vector(&B, argv[2], "B");
  print_vector(B.vector, B.n);

	calculate_C(&C, &C1, &B);
  print_vector(C.vector, C.n);
  print_vector(C1.vector, C1.n);
	save_vector(&C, argv[2], "C");
  save_vector(&C1,argv[2], "C1");

	calculate_O(&O, &B);
 	print_comp_matrix(O);
	save_comp_matrix(&O, argv[2], "O");

	calculate_S_and_R(&S, &R, &B, &C, &O, ratio);
	print_vector(S.vector, S.n);
	print_vector(R.vector, S.n);
	save_comp_vector(&S, argv[2], "S");
	save_comp_vector(&R, argv[2], "R");

	free(B.vector);
  free_comp_matrix(NULL, &O);
	free(S.vector);
	free(R.vector);

	read_ref_vector(&Xi, argv[2], "X");
	revstring(Xi.vector, Xi.n);
	save_ref_vector(&Xi, argv[2], "Xi");
  print_vector(Xi.vector, Xi.n);

	tic("Calc. Suffix Array -> Sadakane direct SAIS");
	calculate_and_save_B(&Xi, argv[2], "Bi");
	toc();

	read_ref_vector(&Bi, argv[2], "Bi");
  print_vector(Bi.vector, Bi.n);

	calculate_O(&Oi, &Bi);
 	print_comp_matrix(Oi);
	save_comp_matrix(&Oi, argv[2], "Oi");

	calculate_S_and_R(&Si, &Ri, &Bi, &C, &Oi, ratio);
	print_vector(Si.vector, Si.n);
	print_vector(Ri.vector, Si.n);
	save_comp_vector(&Si, argv[2], "Si");
	save_comp_vector(&Ri, argv[2], "Ri");

	free(Bi.vector);
  free_comp_matrix(NULL, &Oi);
	free(Si.vector);
	free(Ri.vector);

	free(C.vector);
  free(C1.vector);

	return 0;

}
