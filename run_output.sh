#!/bin/bash
file="bin/parallel_main.exe"
exe="parallel_main.exe"

echo "Making executable"
make $exe > /dev/null

run_simulation () {
echo "control,N,size,time"

for cores in 1 2 3 4
do
	for N in 3000
	do
		for i in {1..5}; do mpiexec -np $cores $file $N 0 0 | grep control; echo  "no output"; 
        mpiexec -np $cores $file $N 1 0 | grep control; echo  "output"; done
	done
done
}

echo "Running the simulation on 4 cores"
run_simulation | tee Data/time.csv
