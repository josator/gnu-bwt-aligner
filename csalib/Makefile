all: csa.a mkcsa csa unbwt mkcst approx

#CC = gcc43
CC = gcc

#CFLAGS = -m64 -D_LARGEFILE64_SOURCE=1 -D_FILE_OFFSET_BITS=64 -msse4.2 -O9 -fPIC
CFLAGS = -m64 -D_LARGEFILE64_SOURCE=1 -D_FILE_OFFSET_BITS=64 -O9 -fPIC

CFILES = csa.c mmap.c psi1.c diskbuf.c lf_dna.c lf_wt.c densearray.c huffman.c comparray.c psi2.c densearray2.c sparsearray2.c sparsearray.c cst.c lf_bit.c
OFILES = csa.o mmap.o psi1.o diskbuf.o lf_dna.o lf_wt.o densearray.o huffman.o comparray.o psi2.o densearray2.o sparsearray2.o sparsearray.o cst.o lf_bit.o
HFILES = typedef.h csa6.h mman.h psi1.h diskbuf.h lf_dna.h lf_wt.h densearray.h huffman.h comparray.h psi2.h densearray2.h sparsearray2.h sparsearray.h lf_bit.h cst.h

csa.a: $(OFILES)
	ar rc csa.a $(OFILES)

mkcsa: mkarray2.c $(OFILES)
	$(CC)  $(CFLAGS) -o mkcsa mkarray2.c csa.a
csa: test_csa.c $(OFILES)
	$(CC)  $(CFLAGS) -o csa test_csa.c csa.a
unbwt: unbwt.c $(OFILES)
	$(CC)  $(CFLAGS) -o unbwt unbwt.c csa.a
mkcst: mkcst.c
	$(CC)  $(CFLAGS) -o mkcst mkcst.c
approx: test_approx.c approx.c csa.a
	$(CC)  $(CFLAGS) -o approx test_approx.c approx.c csa.a

clean:
	rm -f $(OFILES) csa.a mkcsa csa unbwt mkcst approx
