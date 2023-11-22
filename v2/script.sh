make clean
make
mpiexec -n 8 ./samplesort 409 800 5
#mpiexec -n 8 ./samplesort 80000 70000
