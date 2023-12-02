CC=g++
CFLAG= -Wall -I. -O0 -pthread -std=c++17

TARGETS=samplesort samplesortCPPThreads samplesortPThreads

all: $(TARGETS)

samplesort: samplesort.o microtime.o
	$(CC) -g3 -o $@ $^

samplesort.o: samplesort.cpp microtime.h samplesort.h
	$(CC) $(CFLAG) -g3 -c $<

samplesortCPPThreads: samplesortCPPThreads.o microtime.o
	$(CC) -g3 -o $@ $^

samplesortCPPThreads.o: samplesortCPPThreads.cpp microtime.h samplesortCPPThreads.hpp
	$(CC) $(CFLAG) -g3 -c $<

samplesortPThreads: samplesortPThreads.o microtime.o
	$(CC) -g3 -o $@ $^

samplesortPThreads.o: samplesortPThreads.cpp microtime.h samplesortPThreads.hpp
	$(CC) $(CFLAG) -g3 -c $<

microtime.o: microtime.cpp microtime.h
	$(CC) $(CFLAG) -g3 -c $<

clean:
	rm -f *.o *~ core $(TARGETS)
