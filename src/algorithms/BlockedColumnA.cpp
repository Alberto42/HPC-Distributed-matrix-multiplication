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

void BlockedColumnA::createMPICommunicators() {
    MPI_Comm_split(MPI_COMM_WORLD, myProcessRank % numberOfBlocks, myProcessRank, &myGroup);
}

void BlockedColumnA::sparseTimesDense(const CSCMatrix &A, DenseMatrix &B, DenseMatrix &result) {
    for (int i = 1; i < A.m + 1; i++) {
        int extentBegin = A.extents[i - 1] - A.offset;
        int extentEnd = A.extents[i] - A.offset;
        int colA = i - 1 + A.shift;

        for (int j = extentBegin; j < extentEnd; j++) {
            int rowA = A.indices[j];
            double valueA = A.nonzeros[j];
            int rowB = colA;
            int colBBegin = B.shiftHorizontal;
            int colBEnd = B.shiftHorizontal + B.m;

            for (int colB = colBBegin; colB < colBEnd; colB++) {
                double valueB = B.get(rowB, colB - colBBegin);
                result.add(rowA, colB - colBBegin, valueA * valueB);
            }
        }
    }
}

void BlockedColumnA::calcGroups() {
    assert(numProcesses % spec.c == 0);
    numberOfBlocks = numProcesses / spec.c;
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

    log("replicateA");
    replicateA(*localAPencil);

    const int pencilBCWidth = n / numProcesses;
    const int BCShift = myProcessRank * pencilBCWidth;
    localBPencil = makeDenseMatrix(myProcessRank, numProcesses, n, spec.seed, nBeforeExtending);
    localCPencil = makeDenseMatrix(n, pencilBCWidth, BCShift, 0);

    log("main loop");
    for (int j = 0; j < spec.exponent; j++) {
        for (int i = 0; i < numberOfBlocks; i++) {
            log("sparseTimesDense");
            sparseTimesDense(*localAPencil, *localBPencil, *localCPencil);
            if (i == numberOfBlocks - 1 && j == spec.exponent - 1)
                break;
            log("shift");
            shift(localAPencil, localAPencilTmp, MPI_COMM_WORLD);
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