# N=1024
iterMax=10
for i in 1 2 4 8
do
  N=$((65536*$i))
#  echo "$N"
  mpirun-openmpi-mp -np $i ./jacobi-mpi $N $iterMax
done
