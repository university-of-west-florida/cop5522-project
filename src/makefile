CC=g++
CFLAG= -Wall -I. -O0

TARGETS=samplesort open_mp_samplesort

all: $(TARGETS)

samplesort: samplesort.o microtime.o
	$(CC) -g3 -o $@ $^

open_mp_samplesort: open_mp.o microtime.o
	$(CC) $(CFLAG) -fopenmp -g3 -o $@ $^ -lm  # Add -fopenmp and -lm for math library

samplesort.o: samplesort.cpp microtime.h
	$(CC) $(CFLAG) -g3 -c $<

microtime.o: microtime.cpp microtime.h
	$(CC) $(CFLAG) -g3 -c $<

open_mp.o: open_mp.cpp
	$(CC) $(CFLAG) -fopenmp -g3 -c $<  # Add -fopenmp flag

clean:
	rm -f *.o *~ core $(TARGETS)
