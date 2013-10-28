all: bin/csa.a bin/mkcsa bin/csa bin/unbwt bin/mkcst bin/approx

#CC = gcc43
CC = gcc

#CFLAGS = -m64 -D_LARGEFILE64_SOURCE=1 -D_FILE_OFFSET_BITS=64 -msse4.2 -O9 -fPIC
CFLAGS = -m64 -D_LARGEFILE64_SOURCE=1 -D_FILE_OFFSET_BITS=64 -O9 -fPIC

CFILES = csa.c mmap.c psi1.c diskbuf.c lf_dna.c lf_wt.c densearray.c huffman.c comparray.c psi2.c densearray2.c sparsearray2.c sparsearray.c cst.c lf_bit.c
OFILES = csa.o mmap.o psi1.o diskbuf.o lf_dna.o lf_wt.o densearray.o huffman.o comparray.o psi2.o densearray2.o sparsearray2.o sparsearray.o cst.o lf_bit.o
HFILES = typedef.h csa6.h mman.h psi1.h diskbuf.h lf_dna.h lf_wt.h densearray.h huffman.h comparray.h psi2.h densearray2.h sparsearray2.h sparsearray.h lf_bit.h cst.h

bin/csa.a: $(OFILES)
	ar rc bin/csa.a $(OFILES)

bin/mkcsa: mkarray2.c $(OFILES)
	$(CC)  $(CFLAGS) -o bin/mkcsa mkarray2.c bin/csa.a
bin/csa: test_csa.c $(OFILES)
	$(CC)  $(CFLAGS) -o bin/csa test_csa.c bin/csa.a
bin/unbwt: unbwt.c $(OFILES)
	$(CC)  $(CFLAGS) -o bin/unbwt unbwt.c bin/csa.a
bin/mkcst: mkcst.c
	$(CC)  $(CFLAGS) -o bin/mkcst mkcst.c
bin/approx: test_approx.c approx.c bin/csa.a
	$(CC)  $(CFLAGS) -o bin/approx test_approx.c approx.c bin/csa.a

clean:
	rm -f $(OFILES) bin/csa.a bin/mkcsa bin/csa bin/unbwt bin/mkcst bin/approx
