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
