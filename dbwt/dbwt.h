#ifndef _DBWT_DBWT_
#define _DBWT_DBWT_

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>

//#include "mman.h"
#include "utils.h"
#include "queue.h"
#include "sais.h"

#include "../commons/BW_types.h"

#define HSIZ 67777
//#define HSIZ 375559
#define HBSIZ (256*2)

typedef struct {
  int rest[HSIZ];
  int head[HSIZ];
  uchar *buf;
  int bufsiz;
} htbl;

void bwt(uchar *T, long n);

#endif

/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/