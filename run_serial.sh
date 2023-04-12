file="bin/main.exe"

echo "control,N,time"

for N in 1000 2000 3000 4000 5000
do
	for i in {1..5}; do $file $N | grep control; done
done

