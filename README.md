# Prerequisites
GCC

# Build
Run the "Make" command in cop5522-project/src

# SampleSort Serial and MPI Command
./samplesort numNums randomNumBound samplesPerBucket

# Parallel Command
./nameOfProgram numOfProcessors numNums randomNumBound samplesPerBucket

# PThreads Command
./nameOfProgram numNums randomNumBound samplesPerBucket [seed numOfProcessors]

numOfProcessors - The number of processors/threads allocated to use in the parallel framework. Optional in Pthreads and defaults to number of processors in the system
numNums - The number of numbers to sort (They will be automatically generated)  
randomNumBound - The bound of the generated integers (All numbers will fall within [0, randomNumBound])  
samplesPerBucket - When initially pulling samples to approximate evenly-spaced splitters, 
                   this is how many samples are pulled per bucket (See references). Can just use value of 1
seed (optional) - Used to seed the random number generator

# References
http://blog.s-schoener.com/2021-05-24-parallel-sorting/  
https://en.wikipedia.org/wiki/Samplesort  

![Alt text](image.png)

# Timing data automation scripts

`src/script` or `src/script_v2`
Runs `make clean`, `make`, and each program with hardcoded parameters to measure parallelization efficiency.

# dw_openMPVersions and MPIVersions folders
Files from previous attempts used for reference.
