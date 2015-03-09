EXECUTABLES = int_ring int_ring_vector jacobi-mpi
COMPILER = mpicc-openmpi-mp
FLAGS = -O3 -Wall

all: $(EXECUTABLES)

int_ring: int_ring.c
	$(COMPILER) $(FLAGS) int_ring.c -o int_ring

int_ring_vector: int_ring_vector.c
	$(COMPILER) $(FLAGS) int_ring_vector.c -o int_ring_vector

jacobi-mpi: jacobi-mpi.c
	$(COMPILER) $(FLAGS) jacobi-mpi.c -o jacobi-mpi

clean: $(EXECUTABLES)
	rm $(EXECUTABLES)
