CC=g++
CFLAG= -Wall -I. -O0 -fopenmp -lm

TARGETS=samplesort open_mp_samplesort dw_open_mp_samplesort

all: $(TARGETS)

samplesort: samplesort.o microtime.o
	$(CC) $(CFLAG) -g3 -o $@ $^

samplesort.o: samplesort.cpp microtime.h
	$(CC) $(CFLAG) -g3 -c $<

dw_open_mp_samplesort: dw_open_mp_samplesort.o microtime.o
	$(CC) $(CFLAG) -g3 -o $@ $^

dw_open_mp_samplesort.o: dw_open_mp_samplesort.cpp microtime.h
	$(CC) $(CFLAG) -g3 -c $<

open_mp_samplesort: open_mp.o microtime.o
	$(CC) $(CFLAG) -fopenmp -g3 -o $@ $^ -lm  # Add -fopenmp and -lm for math library

open_mp.o: open_mp.cpp
	$(CC) $(CFLAG) -fopenmp -g3 -c $<  # Add -fopenmp flag

microtime.o: microtime.cpp microtime.h
	$(CC) $(CFLAG) -g3 -c $<

clean:
	rm -f *.o *~ core $(TARGETS)
