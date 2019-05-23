//
// Created by albert on 22.05.19.
//

#include <mpi.h>
#include <src/utils.h>
#include <src/parseInput.h>
#include <assert.h>
#include <src/matrices/CSCMatrix.h>
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

void MatmulAlgorithm::extendA(CSCMatrix *fullMatrixA, int numProcesses) {
    assert(fullMatrixA->n == fullMatrixA->m);
    assert(numProcesses <= fullMatrixA->n);
    int n = fullMatrixA->n;
    int tmp = n / numProcesses;
    if (tmp * numProcesses < n) {
        int targetSize = (tmp + 1) * numProcesses;
        fullMatrixA->n = targetSize;
        fullMatrixA->m = targetSize;
        int *newExtents = new int[targetSize + 1];
        for (int i = 0; i < targetSize + 1; i++) {
            if (i < n + 1)
                newExtents[i] = fullMatrixA->extents[i];
            else
                newExtents[i] = fullMatrixA->extents[n];
        }
        delete[] fullMatrixA->extents;
        fullMatrixA->extents = newExtents;
    }
}