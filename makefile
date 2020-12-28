CC=mpic++
CR=mpirun

CPU_COUNT=5

all: main.cpp
	${CC} -std=c++17 -fopenmp -g main.cpp -o output.o
	mpirun -n ${CPU_COUNT} ./output.o

release: main.cpp
	${CC} -O3 -std=c++17 -fopenmp main.cpp -o output.o
	mpirun -n ${CPU_COUNT} ./output.o
