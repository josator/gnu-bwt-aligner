CC = clang

.PHONY: clean

bin/dbwt: dbwt.c queue.c queue.h sais.c utils.c utils.h
	$(CC) -o bin/dbwt dbwt.c queue.c sais.c utils.c -D_LARGEFILE64_SOURCE=1 -D_FILE_OFFSET_BITS=64 -O3

clean:
	rm -f *~ \#*\# bwt *.o dbwt 
