CC = gcc
CFLAGS = -m64 -D_LARGEFILE64_SOURCE=1 -D_FILE_OFFSET_BITS=64 -O3 -std=c99
DEPEN = commons.o string_utils.o

.PHONY: clean

depen: $(DEPEN)

commons.o: commons.c commons.h
	$(CC) $(CFLAGS) -c commons.c

string_utils.o: string_utils.c string_utils.h
	$(CC) $(CFLAGS) -c string_utils.c

clean:
	rm -f *~ \#*\# .*.swp *.o *.a
