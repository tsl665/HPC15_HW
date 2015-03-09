EXECUTABLES = int_ring int_ring_vector
COMPILER = mpicc-openmpi-mp
FLAGS = -O3 -Wall

all: $(EXECUTABLES)

int_ring: int_ring.c
	$(COMPILER) $(FLAGS) int_ring.c -o int_ring

int_ring_vector: int_ring_vector.c
	$(COMPILER) $(FLAGS) int_ring_vector.c -o int_ring_vector

clean: $(EXECUTABLES)
	rm $(EXECUTABLES)
