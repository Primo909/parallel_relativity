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
		for i in {1..2}; do mpiexec -np $cores $file $N 0 0 | grep control; done
	done
done
}

echo "Running the simulstion on 4 cores"
run_simulation | tee Data/time.csv
python3 pyt/time_evaluation.py