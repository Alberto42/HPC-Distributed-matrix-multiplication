//
// Created by albert on 22.05.19.
//

#include <mpi.h>
#include <src/utils.h>
#include <src/parseInput.h>
#include <assert.h>
#include "MatmulAlgorithm.h"


void MatmulAlgorithm::init(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myProcessRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
    initLogger(myProcessRank);
}

void MatmulAlgorithm::calcGroups() {
    assert(numProcesses % spec.c == 0);
    numberOfGroups = numProcesses / spec.c;
}
