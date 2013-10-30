#if defined CSALIB_SEARCH
#include "csalib/csa.h"
#else
#include "commons/BW_types.h"
#endif

typedef struct {

#if defined CSALIB_SEARCH
	CSA csa;
#else
	vector C, C1;
  comp_matrix O;
#endif

} fm_index;

typedef struct {

#if defined CSALIB_SEARCH
	CSA csa;
#else
  comp_vector S, R;
#endif

} suffix_array;

#if defined FM_COMP_32 || FM_COMP_64

#define BWiteration(k_in,l_in, k_out, l_out, b, C, C1, O)\
	do {\
		(k_out) = (C1)->vector[(b)] + getOcompValue((b), (k_in)  , (O));\
		(l_out) = (C)->vector[(b)]  + getOcompValue((b), (l_in)+1, (O));\
	} while (0);
//printf("k-> %lu, l-> %lu, O(k) -> %u, O(l) -> %u, C -> %u, C1 -> %u\n", (k_out), (l_out), getOcompValue((b), (k_in), (O)), getOcompValue((b), (l_in)+1, (O)), (C)->vector[(b)], (C1)->vector[(b)]);
#else

#define BWiteration(k_in, l_in, k_out, l_out, b, C, C1, O)\
	do {\
		(k_out) = (C1)->vector[(b)] + (O)->desp[(b)][(k_in)];\
		(l_out) = (C)->vector[(b)]  + (O)->desp[(b)][(l_in)+1];\
	} while(0);
//printf("k-> %lu, l-> %lu, O(k) -> %u, O(l) -> %u, C -> %u, C1 -> %u\n", k_out, l_out, (O)->desp[(b)][(k_in)], (O)->desp[(b)][(l_in)+1], (C)->vector[(b)], (C1)->vector[(b)]);
#endif
