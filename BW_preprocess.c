#include "BW_preprocess.h"

void encode_reference(ref_vector *X, exome *ex, const char *ref_path) {

	FILE *ref_file;
	ref_file = fopen(ref_path, "r");
	check_file_open(ref_file, ref_path);

	size_t size;

	fseek(ref_file, 0, SEEK_END);
	size = ftell(ref_file) + 2;
	fseek(ref_file, 0, SEEK_SET);

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
	X->vector[total_length] = (REF_TYPE) -1; // $ Symbol
	X->n = total_length + 1; // $ Symbol

	fclose(ref_file);
}

inline SA_TYPE ternary_quicksort_start(SA_TYPE *S, ref_vector *X, SA_TYPE *ranges, SA_TYPE start, SA_TYPE end, uintmax_t h) {

	if (h>1) {
		fprintf(stderr, "Ternary quicksort start only works in the two first iterations\n");
		exit(1);
	}

	if (end<=start) {
		fprintf(stderr, "Start and end range is wrong\n");
		exit(1);
	}

	SA_TYPE start_pivot_pos, start_pivot, end_pivot_pos, end_pivot, l, p, r, value;

	// Put $ suffix at the beginning and ?$ in the end and keep it until sorting ends, just to make second phase faster
	if (h==0) {
		end_pivot_pos = S[end];
		S[end] = S[end - 1];
		S[end - 1] = S[start];
		S[start] = end_pivot_pos;
		start++;
		end--;
	}

	//Select an efficient pivot and put it at the end
	for(value=start; value < end; value++) {
		if (X->vector[(S[value] + h) % X->n] == 1) {end_pivot_pos = S[value]; break;}
	}

	if (value < end) {
		S[value] = S[end];
		S[end] = end_pivot_pos;
	} else {
		printf("ERROR: Pivot not found\n");
		exit(1);
	}

	//Init pivots and perform first pass
	start_pivot_pos = S[start];
	start_pivot = X->vector[(start_pivot_pos + h) % X->n];

	end_pivot_pos = S[end];
	end_pivot = X->vector[(end_pivot_pos + h) % X->n];

	l = start;
	p = start+1;
	r = end;

	while (p != r) {

		value = X->vector[(S[p] + h) % X->n];

		if (value == end_pivot) { //Ternary quicksort (if = does nothing)
			p++;
		} else if (value > end_pivot) {
			S[r] = S[p];
			r--;
			S[p] = S[r];
		} else {
			S[l] = S[p];
			l++;
			S[p] = S[l];
			p++;
		}

	}

	if (start_pivot < end_pivot) {
		S[l] = start_pivot_pos;
		S[r] = end_pivot_pos;
		if (r < end  ) r++;
	} else if (start_pivot > end_pivot) {
		S[r] = start_pivot_pos;
		S[l] = end_pivot_pos;
		if (l > start) l--;
	} else {
		S[r] = start_pivot_pos;
		S[l] = end_pivot_pos;
		if (l > start) l--;
		if (r < end  ) r++;
	}

	ranges[0] = l - start + 1;
	ranges[1] = r - l - 1;

	if (nA==4) {
		//Init pivots for second pass
		start = r;

		start_pivot_pos = S[start];
		start_pivot = X->vector[(start_pivot_pos + h) % X->n];

		end_pivot_pos = S[end];
		end_pivot = X->vector[(end_pivot_pos + h) % X->n];

		l = start;
		p = start+1;
		r = end;

		while (p != r) {

			value = X->vector[(S[p] + h) % X->n];

			if (value == end_pivot) { //Ternary quicksort (if = does nothing)
				p++;
			} else if (value > end_pivot)  {
				S[r] = S[p];
				r--;
				S[p] = S[r];
			} else {
				S[l] = S[p];
				l++;
				S[p] = S[l];
				p++;
			}

		}

		if (start_pivot < end_pivot) {
			S[l] = start_pivot_pos;
			S[r] = end_pivot_pos;
			if (r < end  ) r++;
		} else if (start_pivot > end_pivot) {
			S[r] = start_pivot_pos;
			S[l] = end_pivot_pos;
			if (l > start) l--;
		} else {
			S[r] = start_pivot_pos;
			S[l] = end_pivot_pos;
			if (l > start) l--;
			if (r < end  ) r++;
		}

		if (end_pivot==2) {
			ranges[2] = r - l;
			ranges[3] = end - r + 1;
		} else {
			ranges[2] = l - start + 1;
			ranges[3] = r - l;
		}

	} else if (nA==3) {
		ranges[2] = end - r + 1;
	} else if (nA==2) {
	}

	//Put ?$ on place
	if (h==0) {
		end++;
		value=1;
		int i;
		for (i=0; i < X->vector[(S[end] + h) % X->n]; i++) {
			value += ranges[i];
		}

		end_pivot_pos=S[value];
		end_pivot=X->vector[(S[value] + h) % X->n];
		S[value]=S[end];

		for (; i <nA; i++) {
			value += ranges[i];
			start_pivot_pos = S[value];
			S[value]=end_pivot_pos;
			end_pivot_pos=start_pivot_pos;
		}

		ranges[end_pivot]++;

	} else {

		end_pivot=-1;

	}

	return end_pivot;

}

inline void ternary_quicksort(SA_TYPE *S, SA_TYPE *V, S_SA_TYPE *L, unsigned int start, unsigned int end, size_t h) {

	SA_TYPE stack[512*512];
	intmax_t i=0;

	SA_TYPE start_pivot_pos, start_pivot, end_pivot_pos, end_pivot, left, right, l, p, r, value;

	stack[i++] = start;
	stack[i++] = end;

	if (end<=start) {
		fprintf(stderr, "Start and end range is wrong\n");
		exit(1);
	}

	while (i > 0) {

		right = stack[--i];
		left  = stack[--i];

		if (h>256) {
			value = left + (random() % (right-left));
			end_pivot_pos = S[value];
			S[value] = S[right];
			S[right] = end_pivot_pos;
		}

		//Init pivots
		start_pivot_pos = S[left];
		start_pivot = V[start_pivot_pos];

		end_pivot_pos = S[right];
		end_pivot = V[end_pivot_pos];

		l = left;
		p = left+1;
		r = right;

		while (p != r) {

			value = V[S[p]];

			if (value == end_pivot) { //Ternary quicksort (if = does nothing)
				p++;
			} else if (value > end_pivot) {
				S[r] = S[p];
				r--;
				S[p] = S[r];
			} else {
				S[l] = S[p];
				l++;
				S[p] = S[l];
				p++;
			}

		}

		if (start_pivot < end_pivot) {
			S[l] = start_pivot_pos;
			S[r] = end_pivot_pos;
			if (r < right ) r++;
		} else if (start_pivot > end_pivot) {
			S[l] = end_pivot_pos;
			S[r] = start_pivot_pos;
			if (l > left) l--;
		} else {
			S[r] = start_pivot_pos;
			S[l] = end_pivot_pos;
			if (l > left ) l--;
			if (r < right) r++;
		}

		if (l > left) {
			stack[i++] = left;
			stack[i++] = l;
		}

		if (r < right) {
			stack[i++] = r;
			stack[i++] = right;
		}

	}

}

void calculate_S(comp_vector *S, ref_vector *X) {

	vector V;
	SA_TYPE ranges[nA];
	SA_TYPE offsets[nA];
	SA_TYPE ranges2[nA*nA];
	SA_TYPE offsets2[nA*nA];

	S->siz    = X->n;
	S->n      = S->siz;
	S->ratio  = 1;
	S->vector = (SA_TYPE *) malloc(S->n * sizeof(SA_TYPE));

	for (SA_TYPE i=0; i < X->n; i++) { //TODO: Paralelizar bucle con OpenMP
		S->vector[i] = i;
	}

	print_vector(S->vector, S->n);

	uintmax_t h = 0;
	SA_TYPE dollar;

	printf("**** Suborder: %lu ****\n", h);
	fflush(stdout);

	dollar = ternary_quicksort_start(S->vector, X, ranges, 0, S->siz-1, h);

	for (SA_TYPE i=0, group=1; i < nA; i++) {
		offsets[i] = group;
		group += ranges[i];
	}

	print_vector(S->vector, S->n);
	print_vector(ranges, nA);
	print_vector(offsets, nA);

	h = 1;

	printf("**** Suborder: %lu ****\n", h);
	fflush(stdout);

	S_SA_TYPE **L = (S_SA_TYPE **) malloc(nA * nA * sizeof(S_SA_TYPE *));

#pragma omp parallel for
	for (SA_TYPE r=0; r < nA; r++) {

		if (r==dollar) // Manage ?$
			ternary_quicksort_start(S->vector + offsets[r] + 1, X, ranges2 + nA*r, 0, ranges[r]-2, h);
		else
			ternary_quicksort_start(S->vector + offsets[r]    , X, ranges2 + nA*r, 0, ranges[r]-1, h);

		//Initialize L vectors 
		SA_TYPE aux_desp;
		for (SA_TYPE j=0; j < nA; j++) {

			aux_desp=0;
			if(r==0 && j==0)      aux_desp++;
			if(r==dollar && j==0) aux_desp++;

			L[r*nA + j] = (S_SA_TYPE *) malloc((ranges2[r*nA + j] + aux_desp) * sizeof(S_SA_TYPE));

			if (aux_desp) L[r*nA + j][0] = -aux_desp;
			L[r*nA + j][aux_desp] = (ranges2[r*nA + j]==1)?-1:ranges2[r*nA + j];

			ranges2[r*nA + j] += aux_desp;

			//printIntVector(L[i*nA + j], ranges2[i*nA + j]);

		}

	}

	for (SA_TYPE i=0, group=0; i < nA * nA; i++) {
		offsets2[i] = group;
		group += ranges2[i];
	}

	free(X->vector);

	V.n = X->n;
	V.vector  = (SA_TYPE *) malloc(V.n * sizeof(SA_TYPE));

	print_vector(S->vector, S->n);

	print_vector(ranges2, nA*nA);
	print_vector(offsets2, nA*nA);

	bool end;

	do {

		end = false;
		h = h * 2;

		printf("**** Suborder: %lu ****\n", h);
		fflush(stdout);

#pragma omp parallel for
		for(SA_TYPE r = 0; r < nA*nA; r++) {

			SA_TYPE offset=offsets2[r], offset_last;
			for(SA_TYPE i = 0; i<ranges2[r]; i += imaxabs(L[r][i])) {

				offset_last=offset;
				offset += imaxabs(L[r][i]);

				if (L[r][i]<0) {
					for(SA_TYPE j = offset_last; j < offset; j++)
						V.vector[(S->vector[j] + X->n - h) % X->n] = j;
				} else {
					for(SA_TYPE j = offset_last; j < offset; j++)
						V.vector[(S->vector[j] + X->n - h) % X->n] = offset - 1;
				}

			}

		}

		print_vector(V.vector, V.n);

#pragma omp parallel for
		for(SA_TYPE r = 0; r < nA*nA; r++) {

			intmax_t pos_neg=-1, pos_pos=-1;
			SA_TYPE pos_size=0;
			SA_TYPE offset=offsets2[r], offset_last, increment;
			for(SA_TYPE i = 0; i<ranges2[r]; i += increment) {

				increment = imaxabs(L[r][i]);
				offset_last=offset;
				offset += increment;

				if (L[r][i]<0) {

					if (pos_neg==-1)
						pos_neg=i;
					else
						L[r][pos_neg] += L[r][i];

					continue;

				}

				//printf("++ Sorting segment %ju - %ju\n", offset_last, offset - 1);
				ternary_quicksort(S->vector + offset_last, V.vector, L[r] + i, 0, increment - 1, h);

				//printIntVector(L[r], ranges2[r]);

				pos_pos=-1, pos_size=0;

				SA_TYPE j;
				for(j = i; j < i + increment - 1; j++) {

					if (V.vector[S->vector[offsets2[r] + j]]==V.vector[S->vector[offsets2[r] + j+1]]) {

						pos_neg = -1;

						if (pos_pos==-1) {
							pos_pos = j;
							pos_size = 2;
						} else {
							pos_size++;
						}

					} else {

						if (pos_pos != -1) {
							L[r][pos_pos] = pos_size;
							pos_pos = -1;
							continue;
						}

						if(pos_neg == -1) {
							L[r][j] = -1;
							pos_neg = j;
						} else {
							L[r][pos_neg]--;
						}

						//printf("Condicion restar %lu,%lu %d\n", j, i, pos_neg != -1);

					}

				}

				//printf(":%u %u\n", V.vector[S->vector[j-1]], V.vector[S->vector[j]]);

				if (V.vector[S->vector[offsets2[r] + j-1]]!=V.vector[S->vector[offsets2[r] + j]]) {

					if(pos_neg == -1) {
						L[r][j] = -1;
						pos_neg = j;
					} else {
						L[r][pos_neg]--;
					}

				} else {

					if (pos_pos != -1) {
						L[r][pos_pos] = pos_size;
						pos_pos = -1;
					}

				}

			}

			//printIntVector(L[r], ranges2[r]);
			if(L[r][0] != (S_SA_TYPE) -ranges2[r]) {
				//printf("Ranges => %d != %d\n", L[r][0], (S_SA_TYPE) -ranges2[r]);
				end=true;
			}

		}

		print_vector(S->vector, S->n);

	} while(end);

	free(V.vector);
	for(SA_TYPE i=0; i<nA*nA; i++)
		free(L[i]);

	printf("**** Finished ****\n");
	fflush(stdout);
}

void calculate_B_sadakane_SAIS(ref_vector *B, ref_vector *X) {
		bwt(X->vector, X->n);
}

void calculate_B(ref_vector *B, ref_vector *X, comp_vector *S) {

	B->n = X->n;
	B->vector = (REF_TYPE *) malloc(B->n * sizeof(REF_TYPE));

	SA_TYPE nBlast = B->n - 1;

#pragma omp parallel for
	for(SA_TYPE i=0; i<B->n; i++) {
		B->vector[i] = X->vector[(S->vector[i] + nBlast) % S->n];
	}

#ifdef VERBOSE_DBG
	for(SA_TYPE i=0; i<B->n; i++) {
		for(SA_TYPE j=0; j<B->n; j++) {

			if(X->vector[(S->vector[i]+j) % S->n] == (REF_TYPE) -1)
				printf("-");
			else
				printf("%d", X->vector[(S->vector[i]+j) % S->n]);

		}
		printf("\n");
	}
#endif

}

void calculate_R(comp_vector *R, comp_vector *S) {

	if (S->ratio != 1) {
		fprintf(stderr, "calculateR: The S vector must be uncompressed\n");
		exit(1);
	}

	R->siz   = S->siz;
	R->n     = S->n;
	R->ratio = S->ratio;

	R->vector = (SA_TYPE *) malloc(S->n * sizeof(SA_TYPE));

#pragma omp parallel for
	for(SA_TYPE i=0; i<S->n; i++) {
		R->vector[S->vector[i]] = i;
	}

}

void calculate_C(vector *C, vector *C1, ref_vector *B) {

	SA_TYPE i;

	C->n  = nA;
	C1->n = nA;

	C->vector  = (SA_TYPE *) malloc(C->n  * sizeof(SA_TYPE));
	C1->vector = (SA_TYPE *) malloc(C1->n * sizeof(SA_TYPE));

	check_malloc(C->vector,  "calculateC");
	check_malloc(C1->vector, "calculateC");

	for (i = 0; i < C->n; i++)
		C->vector[i] = 0;

	for (i = 0; i<B->n; i++) {
		if (B->vector[i] == (REF_TYPE) -1) continue;
		C->vector[B->vector[i] + 1]++;
	}

	for (i = 1; i < C->n; i++)
		C->vector[i] += C->vector[i-1];

	for (i = 0; i < C->n; i++)
		C1->vector[i] = C->vector[i] + 1;

}

void calculate_O(comp_matrix *O, ref_vector *B) {

#if defined FM_COMP_32 || FM_COMP_64

	O->siz = B->n+1;          // Position 0 is -1, so I add one element

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
	for (REF_TYPE i=0; i<O->n_count; i++) {

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

		for (SA_TYPE j=0; j<B->n; j++) {

			bit = (j+1) % FM_COMP_VALUE; //First column is -1 index

			if (!bit) {
				pos++;
				O->desp[i][pos]  = k;
				O->count[i][pos] = 0;          //Initialize to 0 bit-vector
			}

			if (B->vector[j] == i) {
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
		for (SA_TYPE j=0; j<B->n; j++) {
			if (((SA_TYPE)B->vector[j]) == i)
				k++;
			O->desp[i][j+1] = k; // First column is -1 index
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
