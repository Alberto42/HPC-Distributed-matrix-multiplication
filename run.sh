#!/usr/bin/env bash

#mpiexec -n 2 ../cmake-build-debug/matrixmul -f /home/albert/HPC/MPI/ac370756/exported_tests/sparse05_00010_001 -s 00085 -c 1 -e 1
#mpiexec -n 2 ../cmake-build-debug/matrixmul -f /home/albert/HPC/MPI/ac370756/exported_tests/sparse05_00010_000 -s 00042 -c 2 -e 1 -g -1
mpiexec -n 1 ../cmake-build-debug/matrixmul -f tests/1 -s 42 -c 1 -e 1 -v -i