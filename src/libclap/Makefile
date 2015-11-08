CC=g++
CFLAGS=-Wall -I./ -Wfatal-errors -O2
CFLAGS+= -std=c++11

all: clap.o

clap.o: clap.hpp clap.cpp
	$(CC) $(CFLAGS) -c clap.cpp

clean:
	rm -f clap.o
