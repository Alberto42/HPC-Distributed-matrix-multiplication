//
// Created by albert on 22.05.19.
//

#include <mpi.h>
#include <src/utils.h>
#include <src/parseInput.h>
#include <assert.h>
#include <src/matrices/CSCMatrix.h>
#include <src/const.h>
#include "MatmulAlgorithm.h"


void MatmulAlgorithm::init(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myProcessRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
    initLogger(myProcessRank);
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

void MatmulAlgorithm::scatterAAmongGroups(CSCMatrix &fullMatrixA, CSCMatrix &localA) {
    vector<CSCMatrix> blocks;
    if (myProcessRank == 0) {
        blocks = fullMatrixA.split(numberOfBlocks);
        localA = blocks[0];
        maxANonzeros = fullMatrixA.count;
        for (int i = 1; i < numberOfBlocks; i++) {
            blocks[i].sendSync(i, INITIAL_SCATTER_TAG);
            MPI_Send(&fullMatrixA.count, 1, MPI_INT, i, SEND_A_NONZEROS_TAG, MPI_COMM_WORLD);
            MPI_Send(&n, 1, MPI_INT, i, SEND_N_TAG, MPI_COMM_WORLD);
            MPI_Send(&nBeforeExtending, 1, MPI_INT, i, SEND_N_BEFORE_EXTENDING_TAG, MPI_COMM_WORLD);
        }
    } else if (myProcessRank < numberOfBlocks) {
        localA.receiveSync(0, INITIAL_SCATTER_TAG);
        MPI_Status status;
        MPI_Recv(&maxANonzeros, 1, MPI_INT, 0, SEND_A_NONZEROS_TAG, MPI_COMM_WORLD, &status);
        MPI_Recv(&n, 1, MPI_INT, 0, SEND_N_TAG, MPI_COMM_WORLD, &status);
        MPI_Recv(&nBeforeExtending, 1, MPI_INT, 0, SEND_N_BEFORE_EXTENDING_TAG, MPI_COMM_WORLD, &status);
    }
}

void MatmulAlgorithm::replicateA(CSCMatrix &localA) {
    MPI_Bcast((void *) &localA, sizeof(CSCMatrix), MPI_BYTE, 0, myGroup);
    if (myProcessRank >= numberOfBlocks) {
        localA.nonzeros = new double[localA.count];
        localA.extents = new int[localA.m + 1];
        localA.indices = new int[localA.count];
    }
    MPI_Bcast((void *) localA.nonzeros, localA.count, MPI_DOUBLE, 0, myGroup);
    MPI_Bcast((void *) localA.extents, localA.m + 1, MPI_INT, 0, myGroup);
    MPI_Bcast((void *) localA.indices, localA.count, MPI_INT, 0, myGroup);

    MPI_Bcast(&maxANonzeros, 1, MPI_INT, 0, myGroup);
    MPI_Bcast(&n, 1, MPI_INT, 0, myGroup);
    MPI_Bcast(&nBeforeExtending, 1, MPI_INT, 0, myGroup);
}