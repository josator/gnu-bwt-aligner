CC = gcc
CFLAGS = -m64 -D_LARGEFILE64_SOURCE=1 -D_FILE_OFFSET_BITS=64 -O3 -std=c99 -fPIC
DEPEN = commons.o string_utils.o

.PHONY: clean

commons.a: $(DEPEN)
	ar rc commons.a $(DEPEN)

clean:
	rm -f *~ \#*\# .*.swp *.o *.a
