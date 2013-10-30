#ifndef _SEARCH_IO_
#define _SEARCH_IO_

#include <string.h>
#include <stdbool.h>

#include "commons/commons.h"
#include "commons/string_utils.h"
#include "commons/BW_types.h"

void load_exome_file(exome *ex, const char *directory);
void save_exome_file(exome *ex, bool reverse, const char *directory);

void encode_reference(ref_vector *X, exome *ex, bool reverse, const char *ref_path);
bool nextFASTAToken(FILE *queries_file, char *uncoded, REF_TYPE *coded, SA_TYPE *nquery);

inline SA_TYPE binsearch(SA_TYPE *array, SA_TYPE size, SA_TYPE key) {

  if( !array ) return 0;

  SA_TYPE *p = array;
  SA_TYPE w;

  while( size > 0 ) {

    w=size/2;
    
    if ( p[w] <= key ) {
      p+=w+1;
      size-=w+1;
    } else {
      size=w;
    }

  }

  return p - array;
}

#endif

/*
This file is part of GNU BWT Aligner.

		Copyright José Salavert Torres, Ignacio Blanquer Espert
							Universitat Politècnica of València
							 - Institut de Instrumentació per a la Imatge Molecular (GRyCAP)
							with partial support of Bull (S.A) Informatique
							2010 - 2013

    GNU BWT Aligner is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    GNU BWT Aligner is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with GNU BWT Aligner.  If not, see <http://www.gnu.org/licenses/>.
*/
