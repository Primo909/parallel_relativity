#!/bin/bash
#file="./bin/parallel_main.exe"

file="bin/parallel_main.exe"
exe="parallel_main.exe"

echo "Making executable"
make $exe > /dev/null
run_simulation () {
echo "control,N,size,time"
for cores in 1 2 3 4
do
	for N in 1000 2000 3000 4000
	do
		for i in {1..5}; do mpiexec -np $cores $file $N 1 0 | grep control; done
	done
done
}

run_simulation | tee Data/time.csv
python3 pyt/time_evaluation.py
