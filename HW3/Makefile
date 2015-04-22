EXECUTABLES = jacobi RBGS
COMPILER = gcc-mp-4.9
FLAGS = -O3 -Wall
OPTION = -fopenmp -o

all: $(EXECUTABLES)

jacobi: Jacobi.c
	$(COMPILER) Jacobi.c $(OPTION) jacobi

RBGS: RBGS.c
	$(COMPILER) RBGS.c $(OPTION) RBGS

clean: $(EXECUTABLES)
	rm $(EXECUTABLES)
