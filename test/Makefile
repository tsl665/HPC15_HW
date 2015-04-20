EXECUTABLES = try_mpi_gather
COMPILER = mpicc-openmpi-mp
FLAGS = -O3 -Wall

all: $(EXECUTABLES)

try_mpi_gather: try_mpi_gather.c
	$(COMPILER) $(FLAGS) try_mpi_gather.c -o try_mpi_gather

clean: $(EXECUTABLES)
	rm $(EXECUTABLES)
