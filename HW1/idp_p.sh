N=8
iterMax=5
for i in 1 2 4 8
do
  echo "======Test for $i Processor(s)===="
  mpirun-openmpi-mp -np $i ./jacobi-mpi $N $iterMax
#  echo "=================================="
done

