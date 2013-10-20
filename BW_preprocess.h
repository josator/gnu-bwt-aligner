#ifndef _BW_PREPROCESS_
#define _BW_PREPROCESS_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "commons/commons.h"
#include "commons/string_utils.h"

#include "BW_csafm.h"

void encode_reference(ref_vector *X, exome *ex, bool duplicate, const char *ref_path);

void calculate_S(comp_vector *S, ref_vector *X);
void calculate_B(ref_vector *B, ref_vector *X, comp_vector *S);
void calculate_R(comp_vector *R, comp_vector *S);
void calculate_C(vector *C, vector *C1, ref_vector *B);
void calculate_O(comp_matrix *O, ref_vector *B);

void compress_SR(comp_vector *SRcomp, comp_vector *SR, SA_TYPE ratio);

void calculateBWTdebug(ref_vector *B, comp_vector *S, ref_vector *X, bool reverse);

#endif
