EXECUTABLES = int_ring
COMPILER = mpicc-openmpi-mp
FLAGS = -O3 -Wall

all: $(EXECUTABLES)

int_ring: int_ring.c
	$(COMPILER) $(FLAGS) int_ring.c -o int_ring

clean: $(EXECUTABLES)
	rm $(EXECUTABLES)
