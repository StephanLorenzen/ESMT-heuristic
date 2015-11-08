CC=g++
# Add -DNDEBUG and -O2 in CFLAGS before testing
CFLAGS=-Wall -I./ -Wfatal-errors # -O2 -DNDEBUG

MAIN=esmt-heuristic

OBJECTS = $(subst cpp,o,$(wildcard *.cpp))
OBJECTS+= $(subst cpp,o,$(wildcard steiner/*.cpp))
OBJECTS+= $(subst cpp,o,$(wildcard steiner/*/*.cpp))
OBJECTS+= $(subst cpp,o,$(wildcard steiner/*/*/*.cpp))
OBJECTS+= $(subst cpp,o,$(wildcard test/*.cpp))

CFLAGS+= -std=c++11


all: $(MAIN) $(OBJECTS)

%.o: %.cpp
	g++ $(CFLAGS) -c -o $@ $^

libclap/clap.o: libclap/clap.cpp
	$(MAKE) -C libclap

$(MAIN): $(OBJECTS) libclap/clap.o
	$(CC) -o $(MAIN) $(OBJECTS) libclap/clap.o

clean:
	rm -f $(MAIN) $(OBJECTS) libclap/clap.o

run:
	./$(MAIN)
