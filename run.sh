#file="./bin/parallel_main.exe"
file="bin/parallel_main.exe"
echo "control,N,size,time"
for cores in 1 2 3 4 5 6 7 8
do
	for N in 500 1000 2000 3000 4000 5000
	do
		for i in {1..8}; do mpiexec -np $cores $file $N | grep control; done
	done
done
