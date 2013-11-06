#ifndef _SEARCH_RUNTIME_
#define _SEARCH_RUNTIME_

#if defined CSALIB_SEARCH
#include "../csalib/csa.h"
#else
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

#define BWiteration(k_in,l_in, k_out, l_out, b, index)\
	do {\
		(k_out) = k_in;\
		(l_out) = l_in;\
		(index)->csa->searchsub((b), (index)->csa, &(k_out), &(l_out));\
	} while (0);
//printf("k-> %lu, l-> %lu\n", (k_out), (l_out));

#else

#define BWiteration(k_in,l_in, k_out, l_out, b, index)\
	do {\
		(k_out) = (index)->C1->vector[(b)] + get_O((b), (k_in) , (index)->O);\
		(l_out) = (index)->C->vector[(b)] + get_O((b), (l_in)+1, (index)->O);\
	} while (0);
//printf("k-> %lu, l-> %lu, O(k) -> %u, O(l) -> %u, C -> %u, C1 -> %u\n", (k_out), (l_out), get_O((b), (k_in), (index)->O), get_O((b), (l_in)+1, (index)->O), (index)->C->vector[(b)], (index)->C1->vector[(b)]);
#endif

#if defined CSALIB_SEARCH
#define size_SA(index) ((index)->csa->n+1)
#define get_SA(m, index) (index)->csa->lookup((index)->csa, (m))
#define get_ISA(m, index) (index)->csa->inverse((index)->csa, (m))
#else
#define size_SA(index) ((index)->S->siz)
#define get_SA(m, index) getScompValue((m), (index)->S, (index)->C, (index)->O)
#define get_ISA(m, index) getRcompValue((m), (index)->R, (index)->C, (index)->O)
#endif

#endif
