#ifndef _SEARCH_RUNTIME_
#define _SEARCH_RUNTIME_

#if defined CSALIB_SEARCH
#include "../csalib/csa.h"
#else
#include "types.h"
#include "csafm.h"
#endif

typedef struct {

#if defined CSALIB_SEARCH
	CSA *csa;
#else
	vector *C, *C1;
	comp_matrix *O;
	comp_vector *S, *R;
#endif

} bwt_index;

#if defined CSALIB_SEARCH

#define BWiterationVariables() uintmax_t k0, l0;

#define BWiteration(k_in,l_in, k_out, l_out, b, index)\
	do {\
		k0 = (index)->csa->K[(b)+1];\
		l0 = (index)->csa->K[(b)+2]-1;\
		(k_out) = (index)->csa->psi_succ((index)->csa,(k_in),k0,l0);\
		(l_out) = (index)->csa->psi_pred((index)->csa,(l_in),k0,l0);\
	} while (0);
//printf("k-> %lu, l-> %lu\n", (k_out), (l_out));

#else

#define BWiterationVariables();

#define BWiteration(k_in,l_in, k_out, l_out, b, index)\
	do {\
		(k_out) = (index)->C1->vector[(b)] + get_O((b), (k_in)  , (index)->O);\
		(l_out) = (index)->C->vector[(b)]  + get_O((b), (l_in)+1, (index)->O);\
	} while (0);
//printf("k-> %lu, l-> %lu, O(k) -> %u, O(l) -> %u, C -> %u, C1 -> %u\n", (k_out), (l_out), getOcompValue((b), (k_in), (O)), getOcompValue((b), (l_in)+1, (O)), (C)->vector[(b)], (C1)->vector[(b)]);

#endif

#if defined CSALIB_SEARCH
#define size_SA(index) ((index)->csa->n+1)
#else
#define size_SA(index) ((index)->S->siz)
#endif

inline uintmax_t get_SA(uintmax_t m, bwt_index *index) {

#if defined CSALIB_SEARCH

	return csa_lookup(index->csa, m);

#else

	if (index->S->ratio==1) {

		return index->S->vector[m];

	} else {

		SA_TYPE i,j;
		uint8_t b_aux;

		i=m; j=0;

		while (i % index->S->ratio) {

			b_aux = get_B_from_O(i, index->O);

			if (b_aux == (uint8_t) -1) {

				i=0;

			} else {

				i = index->C->vector[b_aux] + get_O(b_aux, i+1/*0 is -1*/, index->O);

			}

			j++;

		}

		return (index->S->vector[i / index->S->ratio] + j) % (index->O->siz-1);

	}

#endif

}

inline uintmax_t get_ISA(uintmax_t m, bwt_index *index) {

#if defined CSALIB_SEARCH

		return csa_inverse(index->csa, m);

#else

	if (index->R->ratio==1) {

		return index->R->vector[m];

	} else {

		SA_TYPE i, j, k;
		uint8_t b_aux;

		i = (index->R->ratio - (m % index->R->ratio)) % index->R->ratio;
		k = m + i;

		if (k < index->R->siz) {
			j = index->R->vector[k / index->R->ratio];
		} else {
			j = index->R->vector[0];
			i = index->R->siz - m;
		}

		while (i) {

			b_aux = get_B_from_O(j, index->O);

			if (b_aux == (uint8_t) -1) {

				j=0;

			} else {

				j = index->C->vector[b_aux] + get_O(b_aux, j+1/*0 is -1*/, index->O);

			}

			i--;

		}

		return j;

	}

#endif

}

#endif
