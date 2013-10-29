#include "commons.h"

size_t cur_alloc=0, max_alloc=0;

void *mymalloc(size_t n)
{
  void *p;

  p = malloc(n);
  if (p == NULL) {
    printf("malloc failed.\n");
    exit(1);
  }
  cur_alloc += n;
  if (cur_alloc > max_alloc) {
    max_alloc = cur_alloc;
    //printf("allocated %ld\n",max_alloc);
  }

  if (n == 0) {
    printf("warning: 0 bytes allocated p=%p\n",p);
  }

  return p;
}

void *myrealloc(void *ptr, size_t next, size_t last)
{
  void *p;

  p = realloc(ptr, next);
  if (next > 0 && p == NULL) {
    printf("realloc failed. ptr=%p new=%zd old=%zd\n",ptr,next,last);
    exit(1);
  }
  cur_alloc += next - last;
  if (cur_alloc > max_alloc) {
    max_alloc = cur_alloc;
    //printf("allocated %ld\n",max_alloc);
  }

  return p;
}

void report_mem(const char *s)
{
  puts(s);
  printf("allocated total %zu   max %zu\n",cur_alloc,max_alloc);
}


void myfree(void *p, size_t s)
{
  free(p);
  cur_alloc -= s;
}
