file="parallel_main.exe"
echo "Cleaning Data folder"
make cleanFIG > /dev/null
echo "Making executable"
make $file > /dev/null
echo "Running simulation"
mpiexec -np 4 bin/$file 3000 1 0 > /dev/null
python3 pyt/plots.py
