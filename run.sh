#!/usr/bin/env bash

mpiexec -n 1 matrixmul -f tests/1 -s 42 -c 1 -e 42
#./matrixmul -f tests/1 -s 42 -c 1 -e 42