#!/bin/sh
#SBATCH -p short
#SBATCH -n4
#SBATCH -C beta

start_execution_loop () {
  for i in 1 2 3 4
  do
    echo "---loop $i start---"
    mpirun "$1"
    echo "---loop $i end---"
  done
}

echo "Script is starting..."

echo "-----------------Running scenarios with -n 2-----------------"
echo "Initializing scenario with 3 objects"
mpic++ -std=c++17 -g main_3.cpp -o output_3.o
start_execution_loop output_3.o
printf "\n\n"

echo "Initializing scenario with 4 objects"
mpic++ -std=c++17 -g main_4.cpp -o output_4.o
start_execution_loop output_4.o
printf "\n\n"

echo "Initializing scenario with 5 objects"
mpic++ -std=c++17 -g main_5.cpp -o output_5.o
start_execution_loop output_5.o
printf "\n\n"
echo "-----------------Stopping-----------------"

echo "Script is stopping..."