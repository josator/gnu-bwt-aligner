#ifndef _BW_TYPES_
#define _BW_TYPES_

#include <inttypes.h>

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

typedef uint8_t REF_TYPE;

typedef struct {

	SA_TYPE siz;

	SA_TYPE **desp; // nA
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

  SA_TYPE siz;

	SA_TYPE *vector;
	SA_TYPE n;
  SA_TYPE ratio;

} comp_vector;

typedef struct {

	REF_TYPE *vector;
  SA_TYPE n;

} ref_vector;

#endif
