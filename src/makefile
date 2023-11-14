CC=g++
CFLAG= -Wall -I. -O0

TARGETS=samplesort # add your target here

all: $(TARGETS)

samplesort: samplesort.o microtime.o
	$(CC) -g3 -o $@ $^

samplesort.o: samplesort.cpp microtime.h
	$(CC) $(CFLAG) -g3 -c $<

microtime.o: microtime.cpp microtime.h
	$(CC) $(CFLAG) -g3 -c $<

clean:
	rm -f *.o *~ core $(TARGETS)
