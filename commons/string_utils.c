#include "string_utils.h"

REF_TYPE nA;
REF_TYPE AA = -1, CC = -1, GG = -1, TT = -1;

REF_TYPE table[128];
char *rev_table;

void init_replace_table(const char *str) {

  if (str == NULL) {

    nA = 4;
    AA = 0; CC = 1; GG = 2; TT = 3;

  	rev_table = (char *) malloc(nA * sizeof(char));

    table['a'] = AA; table['A'] = AA;
    table['c'] = CC; table['C'] = CC;
    table['t'] = TT; table['T'] = TT;
    table['g'] = GG; table['G'] = GG;
    table['n'] = AA; table['N'] = AA;

    rev_table[AA] = 'A';
    rev_table[CC] = 'C';
    rev_table[GG] = 'G';
    rev_table[TT] = 'T';

  } else {

    nA = strlen(str);

    rev_table = (char *) malloc(nA * sizeof(char));

    for (REF_TYPE i = 0; i < nA; i++) {
      rev_table[i] = toupper(str[i]);

      table[toupper(str[i])] = i;
      table[tolower(str[i])] = i;

      if      (toupper(str[i]) == 'A') AA = i;
      else if (toupper(str[i]) == 'C') CC = i;
      else if (toupper(str[i]) == 'G') GG = i;
      else if (toupper(str[i]) == 'T') TT = i;
    }

  }

}

void init_table() {
  init_replace_table(NULL);
}

void encode_bases(REF_TYPE* dest, char* src, uintmax_t length) {
  for (uintmax_t i=0; i<length; i++)
    dest[i] = table[(uintmax_t)src[i]];
}

void decode_bases(char* dest, REF_TYPE* src, uintmax_t length) {
  for (uintmax_t i=0; i<length; i++)
    dest[i] = rev_table[(uintmax_t)src[i]];
  dest[length] = '\0';
}

uintmax_t comp4basesInByte(REF_TYPE *X, uintmax_t nX, uint8_t *Y) {

  uintmax_t desp=2, iaux=0;
  uint8_t aux=0;

  if (nX>0) aux = X[0];

  for (uintmax_t i=1; i<nX; i++) {

    if (desp==8) {
      desp=0;
      Y[iaux] = aux;
      iaux++ ;
      aux = aux & 0;
    }

    aux = aux | (X[i] << desp);
    desp = desp +2;

  }

  Y[iaux] = aux;

	uintmax_t nY = nX / 4;
  if (nX % 4) nY++;

  return nY;
}

void revstring(char *X, uintmax_t nX) {

  char tmp;
  uintmax_t i, j;

  if (nX <= 1) return;

  for (i=0, j=nX-1; i<=j; i++, j--) {
    tmp = X[i];
    X[i] = X[j];
    X[j] = tmp;
  }

}
