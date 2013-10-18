CC = g++
CFLAGS = -D_LARGEFILE64_SOURCE=1 -D_FILE_OFFSET_BITS=64 -O3
DEPEN = queue.o utils.o sais.o dbwt.o
OBJECTS = $(DEPEN) main.o

.PHONY: clean all depen

all: bin/dbwt

depen: $(DEPEN)
	cp $(DEPEN) ../

bin/dbwt: $(OBJECTS) queue.h utils.h sais.h dbwt.h
	$(CC) $(CFLAGS) $(OBJECTS) -o bin/dbwt

queue.o: queue.c queue.h utils.h
	$(CC) $(CFLAGS) -c queue.c

utils.o: utils.c utils.h
	$(CC) $(CFLAGS) -c utils.c

sais.o: sais.c sais.h
	$(CC) $(CFLAGS) -c sais.c

dbwt.o: dbwt.c dbwt.h queue.h utils.h sais.h
	$(CC) $(CFLAGS) -c dbwt.c

main.o: main.c dbwt.h
	$(CC) $(CFLAGS) -c main.c

clean:
	rm -f *~ \#*\# bwt *.o bin/dbwt
