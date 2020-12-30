CC=mpic++
CR=mpirun

CPU_COUNT=4

all: main.cpp
	${CC} -std=c++17 -g main.cpp -o output.o
	mpirun -n ${CPU_COUNT} ./output.o

release: main.cpp
	${CC} -O3 -std=c++17 main.cpp -o output.o
	mpirun -n ${CPU_COUNT} ./output.o
