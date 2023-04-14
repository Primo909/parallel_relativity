#!/bin/bash
#file="./bin/parallel_main.exe"

file="bin/parallel_main.exe"
exe="parallel_main.exe"
N="1000"
echo "Making executable"
make $exe > /dev/null

run_simulation () {
echo "control,N,size,time"
for cores in 1 2 3 4
do
	for i in {1..5} 
	do 
		make cleanFIG > /dev/null
		mpiexec -np $cores $file $N $1 0 | grep control
	done
done
}

echo "Running the simulstion on 4 cores"
run_simulation 1 | tee Data/time_yes_save.csv
run_simulation 0 | tee Data/time_no_save.csv
python3 pyt/compare.py
