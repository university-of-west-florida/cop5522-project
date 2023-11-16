make clean
make
rm -f output1.txt output2.txt

n_array=(1000 5000 10000 20000 30000 40000)
Tcount_array=(1 2 4 8 16 20)

#: '
for n in ${n_array[@]}
do
	echo "unoptimized, n=$n" >> output1.txt
	./unoptimized $n $n >> output1.txt
done
#'

for n in ${n_array[@]}
do
	for Tcount in ${Tcount_array[@]}
	do
		export OMP_NUM_THREADS=$Tcount

		for k in {1..20}
		do
			echo "opt1, n=$n, thread_count=$Tcount, run_num=$k" >> output2.txt
			./opt1 $n $n >> output2.txt

			echo "opt2, n=$n, thread_count=$Tcount, run_num=$k" >> output2.txt
			./opt2 $n $n >> output2.txt

			echo "optimized, n=$n, thread_count=$Tcount, run_num=$k" >> output2.txt
			./optimized $n $n >> output2.txt
		done
		echo "matrix size = $n , thread count = $Tcount is complete"
	done
done
