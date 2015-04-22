EXECUTABLES = ssort_gatherv
COMPILER = mpicc-openmpi-mp
FLAGS = -O3 -Wall

all: $(EXECUTABLES)



#ssort: ssort.c
#	$(COMPILER) $(FLAGS) ssort.c -o ssort

ssort_gatherv: ssort_gatherv.c
	$(COMPILER) $(FLAGS) ssort_gatherv.c -o ssort_gatherv

clean: $(EXECUTABLES)
	rm $(EXECUTABLES)
