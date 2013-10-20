#include "BW_csafm.h"

void free_comp_matrix(comp_matrix *reverse, comp_matrix *strand) {

	for (SA_TYPE i=0; i<strand->n_desp; i++) {
		free(strand->desp[i]);
#if defined FM_COMP_32 || FM_COMP_64
		free(strand->count[i]);
#endif
	}

	free(strand->desp);
	if (reverse != NULL) free(reverse->desp);
#if defined FM_COMP_32 || FM_COMP_64
	free(strand->count);
	if (reverse != NULL) free(reverse->count);
#endif

}

//TODO: Arreglar para mapear el complementario reverso cuando trabajo con
// el análisis con menos de 4 bases. ¿¿¿Esto está hecho???
void reverse_strand_C(vector *r_C, vector *s_C, vector *r_C1, vector *s_C1) {

  r_C->n  = s_C->n; r_C1->n = s_C1->n;

	r_C->vector  = (SA_TYPE *)malloc(r_C->n  * sizeof(SA_TYPE));
	check_malloc(r_C->vector,  "reverseStrandC r_C");
	r_C1->vector = (SA_TYPE *)malloc(r_C1->n * sizeof(SA_TYPE));
	check_malloc(r_C1->vector, "reverseStrandC r_C1");

  if (AA != -1 && TT !=-1) {
    r_C->vector[AA] = s_C->vector[TT]; r_C1->vector[AA] = s_C1->vector[TT];
    r_C->vector[TT] = s_C->vector[AA]; r_C1->vector[TT] = s_C1->vector[AA];
  } else if (AA != -1) {
    r_C->vector[AA] = s_C->vector[AA]; r_C1->vector[AA] = s_C1->vector[AA];
  } else if (TT != -1) {
    r_C->vector[TT] = s_C->vector[TT]; r_C1->vector[TT] = s_C1->vector[TT];
  }

  if (CC != -1 && GG !=-1) {
    r_C->vector[CC] = s_C->vector[GG]; r_C1->vector[CC] = s_C1->vector[GG];
    r_C->vector[GG] = s_C->vector[CC]; r_C1->vector[GG] = s_C1->vector[CC];
  } else if (CC != -1) {
    r_C->vector[CC] = s_C->vector[CC]; r_C1->vector[CC] = s_C1->vector[CC];
  } else if (GG != -1) {
    r_C->vector[GG] = s_C->vector[GG]; r_C1->vector[GG] = s_C1->vector[GG];
  }

}

void reverse_strand_O(comp_matrix *r_O, comp_matrix *s_O) {

  r_O->siz = s_O->siz;

  r_O->n_desp = s_O->n_desp;
  r_O->m_desp = s_O->m_desp;

	r_O->desp = (SA_TYPE **) malloc(r_O->n_desp * sizeof(SA_TYPE *));
	check_malloc(r_O->desp, "reverse_strand_O");

  if (AA != -1 && TT !=-1) {
    r_O->desp[AA] = s_O->desp[TT];
    r_O->desp[TT] = s_O->desp[AA];
  } else if (AA != -1) {
    r_O->desp[AA] = s_O->desp[AA];
  } else if (TT != -1) {
    r_O->desp[TT] = s_O->desp[TT];
  }

  if (CC != -1 && GG !=-1) {
    r_O->desp[CC] = s_O->desp[GG];
    r_O->desp[GG] = s_O->desp[CC];
  } else if (CC != -1) {
    r_O->desp[CC] = s_O->desp[CC];
  } else if (GG != -1) {
    r_O->desp[GG] = s_O->desp[GG];
  }

#if defined FM_COMP_32 || FM_COMP_64

  r_O->n_count = s_O->n_count;
  r_O->m_count = s_O->m_count;

	r_O->count = (FM_COMP_TYPE **) malloc(r_O->n_count * sizeof(FM_COMP_TYPE *));
	check_malloc(r_O->count, "reverse_strand_O");

  if (AA != -1 && TT !=-1) {
    r_O->count[AA] = s_O->count[TT];
    r_O->count[TT] = s_O->count[AA];
  } else if (AA != -1) {
    r_O->count[AA] = s_O->count[AA];
  } else if (TT != -1) {
    r_O->count[TT] = s_O->count[TT];
  }

if (CC != -1 && GG !=-1) {
    r_O->count[CC] = s_O->count[GG];
    r_O->count[GG] = s_O->count[CC];
  } else if (CC != -1) {
    r_O->count[CC] = s_O->count[CC];
  } else if (GG != -1) {
    r_O->count[GG] = s_O->count[GG];
  }

#endif

}

void read_vector(vector *vector, const char *directory, const char *name) {

  size_t err=0;
  FILE *fp;

  char path[500];
  path[0]='\0';
  strcat(path, directory);
  strcat(path, "/");
  strcat(path, name);
  strcat(path, ".vec");

  fp  = fopen(path,  "rb+");
  check_file_open(fp, path);

  err = fread(&vector->n, sizeof(SA_TYPE),  1, fp);
  check_file_read(err, 1, path);

  vector->vector = (SA_TYPE *) malloc(vector->n * sizeof(SA_TYPE));
  check_malloc(vector->vector, path);

  err = fread(vector->vector, sizeof(SA_TYPE), vector->n, fp);
  check_file_read(err, vector->n, path);

	fclose(fp);
}

void read_comp_vector(comp_vector *vector, const char *directory, const char *name) {

  size_t err=0;
  FILE *fp;

  char path[500];
  path[0]='\0';
  strcat(path, directory);
  strcat(path, "/");
  strcat(path, name);
  strcat(path, ".vec");

  fp  = fopen(path,  "rb+");
  check_file_open(fp, path);

  err = fread(&vector->n, sizeof(SA_TYPE),  1, fp);
  check_file_read(err, 1, path);

  err = fread(&vector->siz, sizeof(SA_TYPE),  1, fp);
  check_file_read(err, 1, path);

  err = fread(&vector->ratio, sizeof(SA_TYPE),  1, fp);
  check_file_read(err, 1, path);

  vector->vector = (SA_TYPE *) malloc(vector->n * sizeof(SA_TYPE));
  check_malloc(vector->vector, path);

  err = fread(vector->vector, sizeof(SA_TYPE), vector->n, fp);
  check_file_read(err, vector->n, path);

  fclose(fp);
}

void read_ref_vector(ref_vector *vector, const char *directory, const char *name) {

  size_t err=0;
  FILE *fp;

  char path[500];

  path[0]='\0';
  strcat(path, directory);
  strcat(path, "/");
  strcat(path, name);
  strcat(path, ".vec");

  fp  = fopen(path,  "rb+");
  check_file_open(fp, path);

  err = fread(&vector->n, sizeof(SA_TYPE),  1, fp);
  check_file_read(err, 1, path);

  vector->vector = (REF_TYPE *) malloc(vector->n * sizeof(REF_TYPE));
  check_malloc(vector->vector, path);

	err = fread(vector->vector, sizeof(REF_TYPE), vector->n, fp);
  check_file_read(err, vector->n, path);

  fclose(fp);
}

void read_comp_matrix(comp_matrix *matrix, const char *directory, const char *name) {

  size_t err=0;
  FILE *fp;

  char path[500];

  path[0]='\0';
  strcat(path, directory);
  strcat(path, "/");
  strcat(path, name);
  strcat(path, ".desp");

  fp  = fopen(path,  "rb+");
  check_file_open(fp, path);

  err = fread(&matrix->siz,      sizeof(SA_TYPE),  1, fp);
  check_file_read(err, 1, path);

  err = fread(&matrix->n_desp,   sizeof(SA_TYPE),  1, fp);
  check_file_read(err, 1, path);

  err = fread(&matrix->m_desp,   sizeof(SA_TYPE),  1, fp);
  check_file_read(err, 1, path);

	matrix->desp = (SA_TYPE **) malloc(matrix->n_desp * sizeof(SA_TYPE *));
	check_malloc(matrix->desp, path);

  for (SA_TYPE i=0; i<matrix->n_desp; i++) {
    matrix->desp[i] = (SA_TYPE *) malloc(matrix->m_desp * sizeof(SA_TYPE));
    check_malloc(matrix->desp[i], path);

    err = fread(matrix->desp[i], sizeof(SA_TYPE), matrix->m_desp, fp);
    check_file_read(err, matrix->m_desp, path);
  }

	fclose(fp);

#if defined FM_COMP_32 || FM_COMP_64

  path[0]='\0';
  strcat(path, directory);
  strcat(path, "/");
  strcat(path, name);
  strcat(path, ".count");

  fp  = fopen(path,  "rb+");
  check_file_open(fp, path);

  err = fread(&matrix->n_count,   sizeof(SA_TYPE),  1, fp);
  check_file_read(err, 1, path);

  err = fread(&matrix->m_count,   sizeof(SA_TYPE),  1, fp);
  check_file_read(err, 1, path);

	matrix->count = (FM_COMP_TYPE **) malloc(matrix->n_count * sizeof(FM_COMP_TYPE *));
	check_malloc(matrix->count, path);

  for (SA_TYPE i=0; i<matrix->n_count; i++){
    matrix->count[i] = (FM_COMP_TYPE *) malloc(matrix->m_count * sizeof(FM_COMP_TYPE));
    check_malloc(matrix->count[i], path);
    err = fread(matrix->count[i], sizeof(FM_COMP_TYPE), matrix->m_count, fp);
    check_file_read(err, matrix->m_count, path);
  }

  fclose(fp);

#endif

}

void save_vector(vector *vector, const char *directory, const char *name) {

  size_t err=0;
  FILE *fp;

  char path[500];

  path[0]='\0';
  strcat(path, directory);
  strcat(path, "/");
  strcat(path, name);
  strcat(path, ".vec");

  fp  = fopen(path,  "wb+");
  check_file_open(fp, path);

  err = fwrite(&vector->n,     sizeof(SA_TYPE), 1, fp);
  check_file_write(err, 1, path);

  err = fwrite(vector->vector, sizeof(SA_TYPE), vector->n, fp);
  check_file_write(err, vector->n, path);

  fclose(fp);

}

void save_comp_vector(comp_vector *vector, const char *directory, const char *name) {

  size_t err=0;
  FILE *fp;

  char path[500];

  path[0]='\0';
  strcat(path, directory);
  strcat(path, "/");
  strcat(path, name);
  strcat(path, ".vec");

  fp  = fopen(path,  "wb+");
  check_file_open(fp, path);

  err = fwrite(&vector->n,     sizeof(SA_TYPE), 1, fp);
  check_file_write(err, 1, path);

  err = fwrite(&vector->siz,   sizeof(SA_TYPE), 1, fp);
  check_file_write(err, 1, path);

  err = fwrite(&vector->ratio, sizeof(SA_TYPE), 1, fp);
  check_file_write(err, 1, path);

  err = fwrite(vector->vector, sizeof(SA_TYPE), vector->n, fp);
  check_file_write(err, vector->n, path);

  fclose(fp);

}

void save_ref_vector(ref_vector *vector, const char *directory, const char *name) {

  size_t err=0;
  FILE *fp;

  char path[500];

  path[0]='\0';
  strcat(path, directory);
  strcat(path, "/");
  strcat(path, name);
  strcat(path, ".vec");

  fp  = fopen(path,  "wb+");
  check_file_open(fp, path);

  err = fwrite(&vector->n,     sizeof(SA_TYPE), 1, fp);
  check_file_write(err, 1, path);

	err = fwrite(vector->vector, sizeof(REF_TYPE), vector->n, fp);
  check_file_write(err, vector->n, path);

	fclose(fp);

}

void save_comp_matrix(comp_matrix *matrix, const char *directory, const char *name) {

  size_t err=0;
  FILE *fp;

  char path[500];

  path[0]='\0';
  strcat(path, directory);
  strcat(path, "/");
  strcat(path, name);
  strcat(path, ".desp");

  fp  = fopen(path,  "wb+");
  check_file_open(fp, path);

  err = fwrite(&matrix->siz,    sizeof(SA_TYPE), 1, fp);
  check_file_write(err, 1, path);

  err = fwrite(&matrix->n_desp, sizeof(SA_TYPE), 1, fp);
  check_file_write(err, 1, path);

  err = fwrite(&matrix->m_desp, sizeof(SA_TYPE), 1, fp);
  check_file_write(err, 1, path);

  for (SA_TYPE i=0; i<matrix->n_desp; i++) {
    err = fwrite(matrix->desp[i], sizeof(SA_TYPE), matrix->m_desp, fp);
    check_file_write(err, matrix->m_desp, path);
  }

  fclose(fp);

#if defined FM_COMP_32 || FM_COMP_64

  path[0]='\0';
  strcat(path, directory);
  strcat(path, "/");
  strcat(path, name);
  strcat(path, ".count");

  fp  = fopen(path,  "wb+");
  check_file_open(fp, path);

  err = fwrite(&matrix->n_count, sizeof(SA_TYPE), 1, fp);
  check_file_write(err, 1, path);

  err = fwrite(&matrix->m_count, sizeof(SA_TYPE), 1, fp);
  check_file_write(err, 1, path);

  for (SA_TYPE i=0; i<matrix->n_count; i++) {
    err = fwrite(matrix->count[i], sizeof(FM_COMP_TYPE), matrix->m_count, fp);
    check_file_write(err, matrix->m_count, path);
  }

  fclose(fp);

#endif

}

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
