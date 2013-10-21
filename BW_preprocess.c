#include "BW_preprocess.h"

void encode_reference(ref_vector *X, exome *ex, bool duplicate, const char *ref_path) {

	FILE *ref_file;
	ref_file = fopen(ref_path, "r");
	check_file_open(ref_file, ref_path);

	size_t size, read;

	fseek(ref_file, 0, SEEK_END);
	read = ftell(ref_file);
	fseek(ref_file, 0, SEEK_SET);

  if (duplicate) size = 2*read+1;
	else size = read;

	X->vector = (REF_TYPE *) malloc( size * sizeof(REF_TYPE) );
	check_malloc(X->vector, ref_path);

	char *reference = (char *) X->vector;

	if (ex !=NULL) ex->size=0;

	SA_TYPE partial_length=0, total_length=0;

	while ( fgets(reference + total_length, MAXLINE, ref_file) ) {

		if ( (reference + total_length)[0] == '>') {

			if (ex!=NULL) {

				if (total_length == 0) {

					sscanf(reference + total_length, ">%s ", ex->chromosome + ex->size * IDMAX);
					ex->start[ex->size] = 0;

				} else {

					ex->end[ex->size] = partial_length - 1;
					partial_length=0;

					if (ex->size==0)
						ex->offset[0] = 0;
					else
						ex->offset[ex->size] = ex->offset[ex->size-1] + (ex->end[ex->size-1] - ex->start[ex->size-1] + 1);
					ex->size++;

					sscanf(reference + total_length, ">%s ", ex->chromosome + ex->size * IDMAX);
					ex->start[ex->size] = 0;

				}

			}

			continue;

		}

		size_t length = strlen(reference + total_length);
		if ((reference + total_length)[length-1]=='\n')
			length--;

		partial_length += length;
		total_length += length;

	}

	if (ex != NULL) {
		ex->end[ex->size] = partial_length - 1;
		partial_length=0;

		if (ex->size==0)
			ex->offset[0] = 0;
		else
			ex->offset[ex->size] = ex->offset[ex->size-1] + (ex->end[ex->size-1] - ex->start[ex->size-1] + 1);
		ex->size++;
	}

	encode_bases(X->vector, reference, total_length);
	X->n = total_length;
	X->dollar = 0;

	fclose(ref_file);

}

void calculate_R(comp_vector *R, comp_vector *S) {

	if (S->ratio != 1) {
		fprintf(stderr, "calculateR: The S vector must be uncompressed\n");
		exit(1);
	}

	R->siz   = S->siz;
	R->n     = S->n;
	R->ratio = S->ratio;

	R->vector = (SA_TYPE *) malloc(R->n * sizeof(SA_TYPE));

#pragma omp parallel for
	for(SA_TYPE i=0; i<S->n; i++) {
		R->vector[S->vector[i]] = i;
	}

}

void calculate_C(vector *C, vector *C1, ref_vector *B) {

	C->n  = nA;
	C1->n = nA;

	C->vector  = (SA_TYPE *) malloc(C->n  * sizeof(SA_TYPE));
	C1->vector = (SA_TYPE *) malloc(C1->n * sizeof(SA_TYPE));

	check_malloc(C->vector,  "calculateC");
	check_malloc(C1->vector, "calculateC");

	for (SA_TYPE i = 0; i < C->n; i++)
		C->vector[i] = 0;

	SA_TYPE dollar = 0;

	for (SA_TYPE i = 0; i<=B->n; i++) {
		if (i == B->dollar) {dollar = 1; continue;}
		if (B->vector[i - dollar] + 1 == nA) continue;
		C->vector[B->vector[i - dollar] + 1]++;
	}

	for (SA_TYPE i = 1; i < C->n; i++)
		C->vector[i] += C->vector[i-1];

	for (SA_TYPE i = 0; i < C->n; i++)
		C1->vector[i] = C->vector[i] + 1;

}

void calculate_O(comp_matrix *O, ref_vector *B) {

#if defined FM_COMP_32 || FM_COMP_64

	O->siz = B->n+2;          // Position 0 is -1 and the dollar, so I add two elements

	O->n_desp = nA;
	O->m_desp = (O->siz / FM_COMP_VALUE);
	if (O->siz % FM_COMP_VALUE) O->m_desp++;

	O->n_count = nA;
	O->m_count = O->m_desp;

	O->desp  = (SA_TYPE **) malloc(O->n_desp * sizeof(SA_TYPE *));
	check_malloc(O->desp, "calculateO");
	O->count = (FM_COMP_TYPE **) malloc(O->n_count * sizeof(FM_COMP_TYPE *));
	check_malloc(O->count, "calculateO");

#pragma omp parallel for
	for (SA_TYPE i=0; i<O->n_count; i++) {

		O->desp[i]  = (SA_TYPE *) malloc(O->m_desp * sizeof(SA_TYPE));
		check_malloc(O->desp[i], "calculateO");
		O->count[i] = (FM_COMP_TYPE *) malloc(O->m_count * sizeof(FM_COMP_TYPE));
		check_malloc(O->count[i], "calculateO");

		SA_TYPE k=0;
		SA_TYPE pos=0;
		SA_TYPE bit;

		O->desp[i][pos]  = 0;
		O->count[i][pos] = 0;

		/* printf("%d", 0); */

		SA_TYPE dollar = 0;

		for (SA_TYPE j=0; j<=B->n; j++) {

			bit = (j + 1) % FM_COMP_VALUE; //First column is -1 index

			if (!bit) {
				pos++;
				O->desp[i][pos]  = k;
				O->count[i][pos] = 0;          //Initialize to 0 bit-vector
			}

			if (j == B->dollar) { dollar = 1; continue;	}

			if (B->vector[j - dollar] == i) {
				k++;
				O->count[i][pos] |= ((FM_COMP_TYPE) 1) << bit;
			}

		}

	}

#else

	O->siz    = B->n+1;
	O->n_desp = nA;
	O->m_desp = O->siz;       // Position 0 is -1, so I add one element

	O->desp  = (SA_TYPE **) malloc(O->n_desp * sizeof(SA_TYPE *));
	check_malloc(O->desp, "calculateO");

#pragma omp parallel for
	for (SA_TYPE i=0; i<O->n_desp; i++) {

		O->desp[i]  = (SA_TYPE *) malloc(O->m_desp * sizeof(SA_TYPE));
		check_malloc(O->desp[i], "calculateO");

		SA_TYPE k=0;

		O->desp[i][0] = 0;     // First column is -1 index

		SA_TYPE dollar = 0;

		for (SA_TYPE j=0; j<=B->n; j++) {

			if (j == B->dollar) {
				dollar = 1;
			} else if (((SA_TYPE)B->vector[j - dollar]) == i) {
				k++;
			}

			O->desp[i][j + 1] = k; // First column is -1 index

		}

	}

#endif

}

void compress_SR(comp_vector *SR, comp_vector *SRcomp, SA_TYPE ratio) {

	SRcomp->n = (SR->siz / ratio);
	if (SR->siz % ratio) SRcomp->n++;

	SRcomp->siz = SR->siz;
	SRcomp->ratio = ratio;

	SRcomp->vector = (SA_TYPE *) malloc(SRcomp->n * sizeof(SA_TYPE));
	check_malloc(SRcomp->vector, "calculateO");

#pragma omp parallel for
	for (SA_TYPE i=0; i < SRcomp->siz; i+=ratio)
		SRcomp->vector[i/ratio] = SR->vector[i];

}
