#ifndef _SEARCH_CSAFM_
#define _SEARCH_CSAFM_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>

#include "commons/commons.h"
#include "commons/string_utils.h"

#include "commons/BW_types.h"

void reverse_strand_C(vector *r_C, vector *s_C, vector *r_C1, vector *s_C1);
void reverse_strand_O(comp_matrix *r_O, comp_matrix *s_O);

void read_vector(vector *vector, const char *directory, const char *name);
void read_comp_vector(comp_vector *vector, const char *directory, const char *name);
void read_ref_vector(ref_vector *vector, const char *directory, const char *name);
void read_comp_matrix(comp_matrix *matrix, const char *directory, const char *name);

void save_vector(vector *vector, const char *directory, const char *name);
void save_comp_vector(comp_vector *vector, const char *directory, const char *name);
void save_ref_vector(ref_vector *vector, const char *directory, const char *name);
void save_comp_matrix(comp_matrix *matrix, const char *directory, const char *name);

void free_comp_matrix(comp_matrix *reverse, comp_matrix *strand);

#ifdef VERBOSE_DBG

#define print_matrix(M,n,m)\
{\
  printf("Matrix " #M ":\n");\
  for (SA_TYPE i_=0; i_<((SA_TYPE) (n)); i_++) {\
    printf("%ju: ", (uintmax_t) i_);\
    for (SA_TYPE j_=0; j_<((SA_TYPE) (m)); j_++) {\
      printf("%ju ", (uintmax_t) (M)[i_][j_]);\
    }\
    printf("\n");\
  }\
}

#define print_comp_matrix_count(M,n,m,siz)\
{\
	SA_TYPE bit, bbb;\
	printf("Matrix " #M ":\n");\
	for (SA_TYPE i_=0; i_<((SA_TYPE) (n)); i_++) {\
		printf("%ju: ", (uintmax_t) i_);\
		for (SA_TYPE j_=0; j_<((SA_TYPE) (siz)); j_++) {\
			bbb = j_ / FM_COMP_VALUE;\
			bit  = j_ % FM_COMP_VALUE;\
			printf("%ju ", (uintmax_t) (((M)[i_][bbb] >> bit) & 1));\
			if (bit+1 == FM_COMP_VALUE) printf("\t");\
		}\
		printf("\n");\
	}\
}

#if defined FM_COMP_32 || FM_COMP_64

#define print_comp_matrix(O)\
{\
  print_comp_matrix_count((O).count, (O).n_count, (O).m_count, (O).siz);\
  print_matrix((O).desp, (O).n_desp, (O).m_desp);\
}

#else

#define print_comp_matrix(O)\
  print_matrix((O).desp, (O).n_desp, (O).m_desp);
 
#endif

#define print_vector(V,n)\
{\
  printf("Vector " #V ":\n");\
  for (SA_TYPE i_=0; i_<((SA_TYPE) (n)); i_++) {\
		printf("%ju ", (uintmax_t) (V)[i_]);\
  }\
  printf("\n");\
}

#define print_string(S)\
	printf("String " #S ":\n%s\n", S);

#else

#define print_matrix(M,n,m);
#define print_comp_matrix_count_32(M,n,m,siz);
#define print_comp_matrix_count_64(M,n,m,siz);
#define print_comp_matrix(O);
#define print_vector(V,n);
#define print_string(S);

#endif

#ifdef __SSE4_2__
#include <smmintrin.h>
#endif

#if   defined FM_COMP_32

#ifdef __SSE4_2__
#define popcount(x) _mm_popcnt_u32(x)
#else
inline SA_TYPE popcount(uint32_t i) {
  i = i - ((i >> 1) & 0x55555555);
  i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
  return (((i + (i >> 4)) & 0xF0F0F0F) * 0x1010101) >> 24;
}
#endif

#elif defined FM_COMP_64

#ifdef __SSE4_2__
#define popcount(x) _mm_popcnt_u64(x)
#else
inline SA_TYPE popcount(uint64_t i) {
  i = i - ((i >> 1) & 0x5555555555555555);
  i = (i & 0x3333333333333333) + ((i >> 2) & 0x3333333333333333);
  return (((i + (i >> 4)) & 0xF0F0F0F0F0F0F0F) * 0x101010101010101) >> 56;
}
#endif

#endif

#if defined FM_COMP_32 || FM_COMP_64

inline SA_TYPE getOcompValue(SA_TYPE n, SA_TYPE m, comp_matrix *O) {
  SA_TYPE pos, desp;
  pos  = m / FM_COMP_VALUE;
  desp = m % FM_COMP_VALUE;
  return O->desp[n][pos] + popcount( O->count[n][pos] << (FM_COMP_VALUE - (desp + 1)) );
}

#endif

inline REF_TYPE getBfromO(SA_TYPE m, comp_matrix *O) {

#if defined FM_COMP_32 || FM_COMP_64

  m++;

  SA_TYPE pos, desp;
  pos  = m / FM_COMP_VALUE;
  desp = m % FM_COMP_VALUE;

  for(REF_TYPE i=0; i<O->n_count; i++) {
    if ( ( (O->count[i][pos] >> desp) & ((FM_COMP_TYPE) 1)) != 0) return i;
  }

#else

  for(REF_TYPE i=0; i<O->n_desp; i++)
    if ( (O->desp[i][m] < O->desp[i][m+1]) ) return i;

#endif

	return (REF_TYPE) -1;

}

inline SA_TYPE getScompValue(SA_TYPE m, comp_vector *Scomp, vector *C, comp_matrix *O) {

  SA_TYPE i,j;
  REF_TYPE b_aux;
  
  i=m; j=0;
  
  while (i % Scomp->ratio) {
    
    b_aux = getBfromO(i, O);
    
    if (b_aux == (REF_TYPE) -1) {
      
      i=0;

    } else {

#if defined FM_COMP_32 || FM_COMP_64
      i = C->vector[b_aux] + getOcompValue(b_aux, i+1, O);
#else
      i = C->vector[b_aux] + O->desp[b_aux][i+1/*0 is -1*/];
#endif

    }
    
    j++;

  }

  return (Scomp->vector[i / Scomp->ratio] + j) % (O->siz-1);

}

inline SA_TYPE getScompValueB(SA_TYPE m, comp_vector *Scomp, vector *C, comp_matrix *O, ref_vector *B) {
  
  SA_TYPE i, j;
  REF_TYPE b_aux;
  
  i=m; j=0;

  while (i % Scomp->ratio) {
 
    if (i == B->dollar) {

      i=0;

		} else {

			if (i > B->dollar)
    		b_aux = B->vector[i-1];
			else
				b_aux = B->vector[i];

#if defined FM_COMP_32 || FM_COMP_64
      i = C->vector[b_aux] + getOcompValue(b_aux, i+1, O);
#else
      i = C->vector[b_aux] + O->desp[b_aux][i+1/*0 is -1*/];
#endif

    }

    j++;

  }

  return (Scomp->vector[i / Scomp->ratio] + j) % (O->siz-1);

}

inline SA_TYPE getRcompValue(SA_TYPE m, comp_vector *Rcomp, vector *C, comp_matrix *O) {

  SA_TYPE i, j, k;
  REF_TYPE b_aux;

  i = (Rcomp->ratio - (m % Rcomp->ratio)) % Rcomp->ratio;
  k = m + i;

  if (k < Rcomp->siz) {
    j = Rcomp->vector[k / Rcomp->ratio];
  } else {
    j = Rcomp->vector[0];
    i = Rcomp->siz - m;
  }

  while (i) {

    b_aux = getBfromO(j, O);

    if (b_aux == (REF_TYPE) -1) {

      j=0;

    } else {

#if defined FM_COMP_32 || FM_COMP_64
			j = C->vector[b_aux] + getOcompValue(b_aux, j+1, O);
#else
			j = C->vector[b_aux] + O->desp[b_aux][j+1/*0 is -1*/];
#endif

		}

		i--;

	}

	return j;

}

inline SA_TYPE getRcompValueB(SA_TYPE m, comp_vector *Rcomp, vector *C, comp_matrix *O, ref_vector *B) {
  SA_TYPE i, j, k;
  REF_TYPE b_aux;

  i = (Rcomp->ratio - (m % Rcomp->ratio)) % Rcomp->ratio;
  k = m + i;

  if(k < Rcomp->siz) {
    j = Rcomp->vector[k / Rcomp->ratio];
  } else {
    j = Rcomp->vector[0];
    i = Rcomp->siz - m;
  }

  while (i) {


    if (j == B->dollar) {

      j=0;

    } else {

			if (j > B->dollar)
				b_aux = B->vector[j-1];
			else
				b_aux = B->vector[j];

#if defined FM_COMP_32 || FM_COMP_64
      j = C->vector[b_aux] + getOcompValue(b_aux, j+1, O);
#else
      j = C->vector[b_aux] + O->desp[b_aux][j+1/*0 is -1*/];
#endif

    }

    i--;

  }

  return j;

}

//Data structure for chromosome or exome separation positions in the reference

#define INDEX_EXOME 24000
#define IDMAX 100

typedef struct {
  char chromosome[INDEX_EXOME*IDMAX];
  SA_TYPE start[INDEX_EXOME];
  SA_TYPE end[INDEX_EXOME];
  SA_TYPE offset[INDEX_EXOME];
  SA_TYPE size;
} exome;

void load_exome_file(exome *ex, const char *directory);
void save_exome_file(exome *ex, const char *directory);

#endif
