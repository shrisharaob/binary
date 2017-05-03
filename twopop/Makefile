
default: a.out

a.out: main.c globals.h
	g++ -std=c++0x $< -Ofast -O4 -Wall
debug: main.c globals.h
	g++ -g -std=c++0x $< -O0 -Wall
clean:
	rm *.o *.out
.PHONY: clean
