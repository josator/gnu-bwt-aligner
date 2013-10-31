#ifndef STRING_UTILS_H
#define STRING_UTILS_H

#include <ctype.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <inttypes.h>
#include <wchar.h>

extern uint8_t nA;
extern uint8_t AA, CC, GG, TT;

extern uint8_t table[128];
extern char rev_table[4];
extern char reserve[4];

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
void encode_bases(uint8_t* dest, char* src, uintmax_t length);

/**
 *  @brief Decodes a sequence of plain nucleotides
 *  @param dest pointer to destination with plain nucleotides
 *  @param src pointer to char with encoded nucleotides
 *  @param length length of the nucleotide sequence
 * 
 *  Encodes a sequence of plain nucleotides
 */
void decode_bases(char* dest, uint8_t* src, uintmax_t length);

void revstring(uint8_t *X, uintmax_t nX);

void duplicate_reverse(uint8_t *X, uintmax_t nX);

#endif
