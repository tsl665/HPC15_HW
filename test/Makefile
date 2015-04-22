EXECUTABLES = try_mpi_gather try_mpi_gatherv
COMPILER = mpicc-openmpi-mp
FLAGS = -O3 -Wall

all: $(EXECUTABLES)

try_mpi_gather: try_mpi_gather.c
	$(COMPILER) $(FLAGS) try_mpi_gather.c -o try_mpi_gather

try_mpi_gatherv: try_mpi_gatherv.c
	$(COMPILER) $(FLAGS) try_mpi_gatherv.c -o try_mpi_gatherv

clean: $(EXECUTABLES)
	rm $(EXECUTABLES)
