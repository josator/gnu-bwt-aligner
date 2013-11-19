#ifndef COMMONS_H
#define COMMONS_H

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <stdbool.h>

#define MAXLINE      2000
#define MAXLINECOMP  500

//The meaning of Insertion and Deletion may be swapped to the meaning of the SAM file format
#define MATCH	  0
#define DELETION  1
#define MISMATCH  2
#define INSERTION 3

#define timevars() struct timeval t1, t2;

#ifndef min
#define min(x,y) (((x)<(y))?(x):(y))
#endif

#ifndef max
#define max(x,y) (((x)>(y))?(x):(y))
#endif

#define tic(msg)\
  printf("--->> " msg " \n");\
  fflush(stdout);\
  gettimeofday(&t1, NULL);

#define toc()\
  gettimeofday(&t2, NULL);\
  printf("<<--- complete in %.0f usecs\n", (t2.tv_sec-t1.tv_sec)*1e6+(t2.tv_usec-t1.tv_usec));\
  fflush(stdout);

#define start_timer(start_time)\
  gettimeofday(&start_time, NULL);

#define stop_timer(start_time, end_time, timer)\
  gettimeofday(&end_time, NULL);\
  timer = timer + (end_time.tv_sec-start_time.tv_sec)*1e6 + (end_time.tv_usec-start_time.tv_usec);

#define check_syntax(argc, n, params)\
  if ((argc)!=(n)) {\
    fprintf(stderr, "Syntax:\n\t%s\n", (params));\
    exit(1);\
  }

#define check_malloc(D, path)\
  if ((D)==NULL) {\
    fprintf(stderr, "Data structure " #D " in %s is too large\n", (path));\
    exit(1);\
  }

#define check_file_open(fp, path)\
  if (!(fp)) {\
    fprintf(stderr, "Error opening file: %s\n", (path));\
    exit(1);\
  }

#define check_file_read(err, nmemb, path)\
  if ((err) != (nmemb)) {\
    fprintf(stderr, "Error reading file '%s'\n", (path));\
    exit(1);\
  }

#define check_file_write(err, nmemb, path)\
  if ((err) != (nmemb)) {\
    fprintf(stderr, "Error writing file '%s'\n", (path));\
    exit(1);\
  }

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
#define print_vector(V,n);
#define print_string(S);

#endif

void *mymalloc(size_t n);
void *myrealloc(void *ptr, size_t next, size_t last);
void myfree(void *p, size_t s);
void report_mem(const char *s);

extern size_t cur_alloc, max_alloc;

#endif
