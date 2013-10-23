#ifndef STRING_UTILS_H
#define STRING_UTILS_H

#include <ctype.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <inttypes.h>
#include <wchar.h>

#include "BW_types.h"

extern REF_TYPE nA;
extern REF_TYPE AA, CC, GG, TT;

extern REF_TYPE table[128];
extern char *rev_table;

/**
 *  @brief Inits table for nucleotide coding/decoding 
 *  @return void
 * 
 *  Inits table[128] for nucleotide coding/decoding 
 */
void init_table();
void init_replace_table(const char *str);

/**
 *  @brief Encodes a sequence of plain nucleotides
 *  @param dest pointer to destination with encoded nucleotides
 *  @param src pointer to char with plain nucleotides
 *  @param length length of the nucleotide sequence
 * 
 *  Encodes a sequence of plain nucleotides
 */
void encode_bases(REF_TYPE* dest, char* src, uintmax_t length);

/**
 *  @brief Decodes a sequence of plain nucleotides
 *  @param dest pointer to destination with plain nucleotides
 *  @param src pointer to char with encoded nucleotides
 *  @param length length of the nucleotide sequence
 * 
 *  Encodes a sequence of plain nucleotides
 */
void decode_bases(char* dest, REF_TYPE* src, uintmax_t length);

void revstring(REF_TYPE *X, uintmax_t nX);

#endif
