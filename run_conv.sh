file="bin/parallel_main.exe"

echo "control,N,time"

for N in 1000 2000 3000 4000 5000
do
	for i in {1..5}; do mpiexec -np 4 $file $N 0 1 ""Data/pointConvTest"$N".dat"" | grep control; done
done

