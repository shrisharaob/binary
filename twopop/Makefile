CC=g++
CFLAGS=-Ofast -O4 -Wall
DEPS=globals.h

default: a.out

a.out: main.cpp $(DEPS) 
	$(CC) -std=c++0x $< $(CFLAGS)
debug: main.cpp globals.h
	g++ -g -std=c++0x $< -O0 -Wall
clean:
	rm *.o *.out
.PHONY: clean
