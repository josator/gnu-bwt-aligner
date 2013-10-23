#include "BW_io.h"

void load_exome_file(exome *ex, const char *directory) {

  FILE *fp;

  char path[500];
  path[0]='\0';
  strcat(path, directory);
  strcat(path, "/index");

  fp  = fopen(path,  "r");
  check_file_open(fp, path);

  char line[MAXLINE];

  ex->offset[0]=0;
  ex->size = 0;

  while (fgets(line, MAXLINE, fp) ) {

    if (line[0]=='>') {

			int j;
			char c = 0;

			for(j=0; j<IDMAX-1; j++) {
				c = line[j+1];
				if (c==' ') break;
				ex->chromosome[ex->size*IDMAX+j] = c;
			}

			ex->chromosome[ex->size*IDMAX+j] = '\0';

			sscanf(line + j + 2, "%ju %ju %*s", &ex->start[ex->size], &ex->end[ex->size]);
			ex->size++;
			ex->offset[ex->size] = ex->offset[ex->size-1] + (ex->end[ex->size-1] - ex->start[ex->size-1]+1);

		}

	}

	fclose(fp);

}

void save_exome_file(exome *ex, const char *directory) {

  FILE *fp;

  char path[500];
  path[0]='\0';
  strcat(path, directory);
  strcat(path, "/index");

  fp  = fopen(path, "w");
  check_file_open(fp, path);

  for(SA_TYPE i=0; i<ex->size; i++) {
    fprintf(fp, ">%s %ju %ju\n", ex->chromosome + i*IDMAX, (uintmax_t) ex->start[i], (uintmax_t) ex->end[i]);
  }

}

void encode_reference(ref_vector *X, exome *ex, const char *ref_path) {

	FILE *ref_file;
	ref_file = fopen(ref_path, "r");
	check_file_open(ref_file, ref_path);

	size_t size;

	fseek(ref_file, 0, SEEK_END);
	size = ftell(ref_file) + 1; //Valgrind errors on dbwt
	fseek(ref_file, 0, SEEK_SET);

	X->vector = (REF_TYPE *) malloc( size * sizeof(REF_TYPE) );
	check_malloc(X->vector, ref_path);

	char *reference = (char *) X->vector;

	if (ex !=NULL) ex->size=0;

	SA_TYPE partial_length=0, total_length=0;

	while ( fgets(reference + total_length, MAXLINE, ref_file) ) {

		if ( (reference + total_length)[0] == '>') {

			if (ex!=NULL) {

				if (total_length == 0) {

					sscanf(reference + total_length, ">%s ", ex->chromosome + ex->size * IDMAX);
					ex->start[ex->size] = 0;

				} else {

					ex->end[ex->size] = partial_length - 1;
					partial_length=0;

					if (ex->size==0)
						ex->offset[0] = 0;
					else
						ex->offset[ex->size] = ex->offset[ex->size-1] + (ex->end[ex->size-1] - ex->start[ex->size-1] + 1);
					ex->size++;

					sscanf(reference + total_length, ">%s ", ex->chromosome + ex->size * IDMAX);
					ex->start[ex->size] = 0;

				}

			}

			continue;

		}

		size_t length = strlen(reference + total_length);
		if ((reference + total_length)[length-1]=='\n')
			length--;

		partial_length += length;
		total_length += length;

	}

	if (ex != NULL) {
		ex->end[ex->size] = partial_length - 1;
		partial_length=0;

		if (ex->size==0)
			ex->offset[0] = 0;
		else
			ex->offset[ex->size] = ex->offset[ex->size-1] + (ex->end[ex->size-1] - ex->start[ex->size-1] + 1);
		ex->size++;
	}

	encode_bases(X->vector, reference, total_length);
	X->n = total_length;
	X->dollar = 0;
	X->vector[X->n] = 0; //Valgrind errors on dbwt

	fclose(ref_file);

}

bool nextFASTAToken(FILE *queries_file, char *uncoded, REF_TYPE *coded, SA_TYPE *nquery) {

	char line[MAXLINE];
	size_t length=0;

	*nquery=0;
	uncoded[0]='\0';

	while ( fgets(line, MAXLINE, queries_file) ) {

		if (line[0] == '>') {
			if (*nquery) break;
			else continue;
		}

		length=strlen(line);
		if (line[length-1]=='\n')
			length--;

		uncoded[*nquery] = '\0';
		strncpy(uncoded + *nquery, line , length);

		*nquery += length;

	}

	if (*nquery) {

		encode_bases(coded, uncoded, *nquery);

		return true;

	} else {

		return false;

	}

}
