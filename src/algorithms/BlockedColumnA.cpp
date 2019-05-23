//
// Created by albert on 22.05.19.
//

#include "BlockedColumnA.h"

#include <iostream>
#include <string>
#include <mpi.h>
#include <fstream>
#include <sys/stat.h>
#include <algorithm>
#include <src/matrices/CSCMatrix.h>
#include <src/utils.h>
#include <src/matrices/DenseMatrix.h>
#include "src/parseInput.h"
#include "src/algorithms/BlockedColumnA.h"
#include "src/algorithms/BlockedInnerABC.h"
#include "src/const.h"

using namespace std;

void BlockedColumnA::scatterAAmongGroups(CSCMatrix &fullMatrixA, CSCMatrix &localAPencil) {
    vector<CSCMatrix> pencils;
    if (myProcessRank == 0) {
        pencils = fullMatrixA.split(numberOfGroups);
        localAPencil = pencils[0];
        maxANonzeros = fullMatrixA.count;
        for (int i = 1; i < numberOfGroups; i++) {
            pencils[i].sendSync(i, INITIAL_SCATTER_TAG);
            MPI_Send(&fullMatrixA.count, 1, MPI_INT, i, SEND_A_NONZEROS_TAG, MPI_COMM_WORLD);
            MPI_Send(&n, 1, MPI_INT, i, SEND_N_TAG, MPI_COMM_WORLD);
            MPI_Send(&nBeforeExtending, 1, MPI_INT, i, SEND_N_BEFORE_EXTENDING_TAG, MPI_COMM_WORLD);
        }
    } else if (myProcessRank < numberOfGroups) {
        localAPencil.receiveSync(0, INITIAL_SCATTER_TAG);
        MPI_Status status;
        MPI_Recv(&maxANonzeros, 1, MPI_INT, 0, SEND_A_NONZEROS_TAG, MPI_COMM_WORLD, &status);
        MPI_Recv(&n, 1, MPI_INT, 0, SEND_N_TAG, MPI_COMM_WORLD, &status);
        MPI_Recv(&nBeforeExtending, 1, MPI_INT, 0, SEND_N_BEFORE_EXTENDING_TAG, MPI_COMM_WORLD, &status);
    }
}

void BlockedColumnA::replicateAPencils(CSCMatrix &localAPencil) {
    MPI_Bcast((void *) &localAPencil, sizeof(CSCMatrix), MPI_BYTE, 0, myGroup);
    if (myProcessRank >= numberOfGroups) {
        localAPencil.nonzeros = new double[localAPencil.count];
        localAPencil.extents = new int[localAPencil.m + 1];
        localAPencil.indices = new int[localAPencil.count];
    }
    MPI_Bcast((void *) localAPencil.nonzeros, localAPencil.count, MPI_DOUBLE, 0, myGroup);
    MPI_Bcast((void *) localAPencil.extents, localAPencil.m + 1, MPI_INT, 0, myGroup);
    MPI_Bcast((void *) localAPencil.indices, localAPencil.count, MPI_INT, 0, myGroup);

    MPI_Bcast(&maxANonzeros, 1, MPI_INT, 0, myGroup);
    MPI_Bcast(&n, 1, MPI_INT, 0, myGroup);
    MPI_Bcast(&nBeforeExtending, 1, MPI_INT, 0, myGroup);
}

void BlockedColumnA::createMPICommunicators() {
    MPI_Comm_split(MPI_COMM_WORLD, myProcessRank % numberOfGroups, myProcessRank, &myGroup);
}

void BlockedColumnA::sparseTimesDense(const CSCMatrix &A, DenseMatrix &B, DenseMatrix &result) {
    for (int i = 1; i < A.m + 1; i++) {
        int extentBegin = A.extents[i - 1] - A.offset;
        int extentEnd = A.extents[i] - A.offset;
        int colA = i - 1 + A.shiftHorizontal;

        for (int j = extentBegin; j < extentEnd; j++) {
            int rowA = A.indices[j];
            double valueA = A.nonzeros[j];
            int rowB = colA;
            int colBBegin = B.shift;
            int colBEnd = B.shift + B.m;

            for (int colB = colBBegin; colB < colBEnd; colB++) {
                double valueB = B.get(rowB, colB - colBBegin);
                result.add(rowA, colB - colBBegin, valueA * valueB);
            }
        }
    }
}

void BlockedColumnA::shift(CSCMatrix *&localAPencil, CSCMatrix *&localAPencilTmp) {
    MPI_Request requests[8];
    MPI_Status statuses[8];

    localAPencil->sendAsync((myProcessRank + 1) % numProcesses, SHIFT_TAGS, requests);

    int src = (myProcessRank - 1 + numProcesses) % numProcesses;

    MPI_Irecv((void *) localAPencilTmp, sizeof(CSCMatrix), MPI_BYTE, src, SHIFT_TAGS[0], MPI_COMM_WORLD, requests + 4);

    double *nonzeros = new double[maxANonzeros];
    int *extents = new int[localAPencil->m + 1];
    int *indices = new int[maxANonzeros];

    MPI_Irecv(nonzeros, maxANonzeros, MPI_DOUBLE, src, SHIFT_TAGS[1],
              MPI_COMM_WORLD, requests + 5);
    MPI_Irecv(extents, localAPencil->m + 1, MPI_INT, src, SHIFT_TAGS[2], MPI_COMM_WORLD,
              requests + 6);
    MPI_Irecv(indices, maxANonzeros, MPI_INT, src, SHIFT_TAGS[3], MPI_COMM_WORLD,
              requests + 7);

    MPI_Waitall(8, requests, statuses);

    localAPencilTmp->nonzeros = nonzeros;
    localAPencilTmp->extents = extents;
    localAPencilTmp->indices = indices;

    swap(localAPencil, localAPencilTmp);

    delete[] localAPencilTmp->nonzeros;
    delete[] localAPencilTmp->extents;
    delete[] localAPencilTmp->indices;
}

DenseMatrix *BlockedColumnA::gatherResultVerbose(DenseMatrix *localCPencil) {
    MPI_Datatype dtDenseMatrix;
    const size_t localCPencilSize = sizeof(DenseMatrix) + n * localCPencil->m * sizeof(double);
    MPI_Type_contiguous(localCPencilSize, MPI_BYTE, &dtDenseMatrix);
    MPI_Type_commit(&dtDenseMatrix);

    DenseMatrix *receiverCMatrices = nullptr;
    if (myProcessRank == 0) {
        receiverCMatrices = (DenseMatrix *) malloc(numProcesses * localCPencilSize);
    }
    MPI_Gather((void *) localCPencil, 1, dtDenseMatrix, receiverCMatrices, 1, dtDenseMatrix, 0, MPI_COMM_WORLD);
    return receiverCMatrices;
}

int BlockedColumnA::gatherResultGreater(DenseMatrix *localCPencil) {
    int localGreaterCount = 0;
    for (int row = 0; row < nBeforeExtending; row++) {
        for (int col = 0; localCPencil->shift + col < nBeforeExtending && col < localCPencil->m; col++) {
            localGreaterCount += (localCPencil->get(row, col) >= spec.g);
        }
    }
    int *buff = nullptr;
    if (myProcessRank == 0) {
        buff = new int[numProcesses];
    }
    MPI_Gather(&localGreaterCount, 1, MPI_INT, buff, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (myProcessRank == 0) {
        int totalCount = 0;
        for (int i = 0; i < numProcesses; i++)
            totalCount += buff[i];
        return totalCount;
    } else
        return -1;

}

void BlockedColumnA::printResult(DenseMatrix *receiverCMatrices) {
    if (myProcessRank == 0) {
        cout << nBeforeExtending << " " << nBeforeExtending << endl;

        for (int row = 0; row < nBeforeExtending; row++) {
            for (int i = 0; i < numProcesses; i++) {
                DenseMatrix &m = *getIthMatrix(receiverCMatrices, i);
                for (int colInM = 0; colInM < m.m && m.m * i + colInM < nBeforeExtending; colInM++) {
                    cout.precision(5);
                    cout << "   " << fixed << m.get(row, colInM) << " ";
                }
            }
            cout << endl;
        }
    }
}

void BlockedColumnA::assignCMatrixToBMatrix(DenseMatrix *localBPencil, DenseMatrix *localCPencil) {
    assert(localBPencil->size() == localCPencil->size());
    memcpy(localBPencil, localCPencil, localBPencil->size());
    for (int i = 0; i < localCPencil->n * localCPencil->m; i++)
        localCPencil->matrix[i] = 0;

}



void BlockedColumnA::columnAAlgorithm(int argc, char **argv) {

    CSCMatrix fullMatrixA, *localAPencil, *localAPencilTmp;

    localAPencil = new CSCMatrix();
    localAPencilTmp = new CSCMatrix();
    DenseMatrix *localBPencil, *localCPencil;
    int greaterCount;
    DenseMatrix *receiverCMatrices;
    double startTime, endTime;

    init(argc, argv);
    calcGroups();
    createMPICommunicators();

    log("start algorithm ------------------------------------------------------------");
    if (myProcessRank == 0) {
        ifstream ifs = ifstream(spec.file);

        ifs >> fullMatrixA;

        nBeforeExtending = fullMatrixA.n;
        extendA(&fullMatrixA, numProcesses);
        n = fullMatrixA.n;

    }

    scatterAAmongGroups(fullMatrixA, *localAPencil);


    MPI_Barrier(MPI_COMM_WORLD);
    if (myProcessRank == 0) {
        startTime = MPI_Wtime();
    }

    log("replicateAPencils");
    replicateAPencils(*localAPencil);

    const int pencilBCWidth = n / numProcesses;
    const int BCOffset = myProcessRank * pencilBCWidth;
    localBPencil = makeDenseMatrix(myProcessRank, numProcesses, n, spec.seed, nBeforeExtending);
    localCPencil = makeDenseMatrix(n, pencilBCWidth, BCOffset);

    log("main loop");
    for (int j = 0; j < spec.exponent; j++) {
        for (int i = 0; i < numberOfGroups; i++) {
            log("sparseTimesDense");
            sparseTimesDense(*localAPencil, *localBPencil, *localCPencil);
            if (i == numberOfGroups - 1 && j == spec.exponent - 1)
                break;
            log("shift");
            shift(localAPencil, localAPencilTmp);
        }
        if (j != spec.exponent - 1)
            assignCMatrixToBMatrix(localBPencil, localCPencil);
    }

    log("gatherResult");

    if (spec.verbose) {
        receiverCMatrices = gatherResultVerbose(localCPencil);
    } else {
        greaterCount = gatherResultGreater(localCPencil);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (myProcessRank == 0) {
        endTime = MPI_Wtime();
    }
    log("printResult");
    if (spec.verbose) {
        printResult(receiverCMatrices);
    } else {
        if (myProcessRank == 0)
            cout << greaterCount << endl;
    }

    MPI_Finalize();
    log("after finalize");
}