#ifndef _SEARCH_TYPES_
#define _SEARCH_TYPES_

#include <inttypes.h>
#include <stdbool.h>

#if   defined SA_64

typedef uint64_t SA_TYPE;
typedef int64_t S_SA_TYPE;

#elif defined SA_32

typedef uint32_t SA_TYPE;
typedef int32_t S_SA_TYPE;

#elif defined SA_16

typedef uint16_t SA_TYPE;
typedef int16_t S_SA_TYPE;

#elif defined SA_8

typedef uint8_t SA_TYPE;
typedef int8_t S_SA_TYPE;

#else

typedef uint32_t SA_TYPE;
typedef int32_t S_SA_TYPE;

#endif

#if   defined FM_COMP_64

typedef uint64_t FM_COMP_TYPE;
#define FM_COMP_VALUE 64

#elif defined FM_COMP_32

typedef uint32_t FM_COMP_TYPE;
#define FM_COMP_VALUE 32

#endif

typedef struct {

	SA_TYPE **desp; // nA
	SA_TYPE siz; //Real number of columns of the uncompressed matrix
  SA_TYPE n_desp;
  SA_TYPE m_desp;

#if defined FM_COMP_32 || FM_COMP_64
  FM_COMP_TYPE **count; // nA
  SA_TYPE n_count;
  SA_TYPE m_count;
#endif

} comp_matrix;

typedef struct {

	SA_TYPE *vector;
  SA_TYPE n;

} vector;

typedef struct {

	SA_TYPE *vector;

  SA_TYPE siz; //Real size of the uncompressed vector
	SA_TYPE n;
  SA_TYPE ratio;

} comp_vector;

typedef struct {

	uint8_t *vector;
  SA_TYPE n;
	SA_TYPE dollar; //Position ending with the $ symbol (the first in the reference)

} ref_vector;

#endif