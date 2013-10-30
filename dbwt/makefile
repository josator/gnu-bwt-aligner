COMMONS_DIR = ../commons/

CC = g++
CFLAGS = -m64 -D_LARGEFILE64_SOURCE=1 -D_FILE_OFFSET_BITS=64 -g
DEPEN = queue.o utils.o sais.o dbwt.o
OBJECTS = $(DEPEN) main.o

.PHONY: clean all depen

all: commons bin/dbwt

depen: $(DEPEN)

bin/dbwt: $(OBJECTS) queue.h utils.h sais.h dbwt.h
	$(CC) $(CFLAGS) $(COMMONS_DIR)/commons.o $(OBJECTS) -o bin/dbwt

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

commons:
	make -C $(COMMONS_DIR) depen

clean:
	rm -f *~ \#*\# *.o
