CC = g++
CFLAGS = -m64 -g
DEPEN = commons.o string_utils.o

.PHONY: clean

depen: $(DEPEN)

commons.o: commons.c commons.h
	$(CC) $(CFLAGS) -c commons.c

string_utils.o: string_utils.c string_utils.h
	$(CC) $(CFLAGS) -c string_utils.c

clean:
	rm -f *~ \#*\# bwt *.o
