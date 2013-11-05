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
#define get_SA(m, index) (index)->csa->lookup((index)->csa, (m))
#define get_ISA(m, index) (index)->csa->inverse((index)->csa, (m))
#else
#define size_SA(index) ((index)->S->siz)
#define get_SA(m, index) getScompValue((m), (index)->S, (index)->C, (index)->O)
#define get_ISA(m, index) getRcompValue((m), (index)->R, (index)->C, (index)->O)
#endif

#endif
