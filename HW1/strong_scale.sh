N=65536
iterMax=100
for i in 1 2 4 8 16
do
     mpirun-openmpi-mp -np $i ./jacobi-mpi $N $iterMax
done
