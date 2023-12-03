CC=g++
CFLAG= -Wall -I. -O0
MPICC=mpic++
OMPCFLAG= -Wall -I. -O0 -fopenmp -lm
PTHREADCFLAG= -Wall -I. -O0 -pthread -std=c++17

TARGETS=samplesort open_mp_samplesort dw_open_mp_samplesort mpi_samplesort samplesortCPPThreads samplesortPThreads

all: $(TARGETS)

samplesort: samplesort.o microtime.o
	$(CC) $(CFLAG) -g3 -o $@ $^

samplesort.o: samplesort.cpp microtime.h samplesort.h
	$(CC) $(CFLAG) -g3 -c $<

dw_open_mp_samplesort: dw_open_mp_samplesort.o microtime.o
	$(CC) $(OMPCFLAG) -g3 -o $@ $^

dw_open_mp_samplesort.o: dw_open_mp_samplesort.cpp microtime.h
	$(CC) $(OMPCFLAG) -g3 -c $<

open_mp_samplesort: open_mp.o microtime.o
	$(CC) $(CFLAG) -fopenmp -g3 -o $@ $^ -lm  # Add -fopenmp and -lm for math library

open_mp.o: open_mp.cpp
	$(CC) $(CFLAG) -fopenmp -g3 -c $<  # Add -fopenmp flag

samplesortCPPThreads: samplesortCPPThreads.o microtime.o
	$(CC) $(PTHREADCFLAG) -g3 -o $@ $^

samplesortCPPThreads.o: samplesortCPPThreads.cpp microtime.h samplesortCPPThreads.hpp
	$(CC) $(PTHREADCFLAG) -g3 -c $<

samplesortPThreads: samplesortPThreads.o microtime.o
	$(CC) $(PTHREADCFLAG) -g3 -o $@ $^

samplesortPThreads.o: samplesortPThreads.cpp microtime.h samplesortPThreads.hpp
	$(CC) $(PTHREADCFLAG) -g3 -c $<

microtime.o: microtime.cpp microtime.h
	$(CC) $(CFLAG) -g3 -c $<

mpi_samplesort: mpi_samplesort.o mpi_microtime.o
	$(MPICC) -g3 -o $@ $^

mpi_samplesort.o: mpi_samplesort.cpp microtime.h
	$(MPICC) $(CFLAG) -g3 -c $<

mpi_microtime.o: microtime.cpp microtime.h
	$(MPICC) $(CFLAG) -g3 -c -o $@ $<

clean:
	rm -f *.o *~ core $(TARGETS)
