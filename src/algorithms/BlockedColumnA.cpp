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
#include <cstring>
#include <chrono>

using namespace std;
using namespace std::chrono;

void BlockedColumnA::createMPICommunicators() {
    MPI_Comm_split(MPI_COMM_WORLD, myProcessRank % numberOfBlocks, myProcessRank, &myGroup);
}

void BlockedColumnA::sparseTimesDense(const CSCMatrix &A, DenseMatrix &B, DenseMatrix &result) {

    steady_clock::time_point begin = steady_clock::now();
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
    steady_clock::time_point end = steady_clock::now();
    sparseTimeDenseTotalTime += timeDiffInMs(begin,end);
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
    steady_clock::time_point assignStartTime = steady_clock::now();
    assert(localBPencil->size() == localCPencil->size());
    memcpy(localBPencil, localCPencil, localBPencil->size());
    for (int i = 0; i < localCPencil->n * localCPencil->m; i++) {
        localCPencil->matrix[i] = 0;
    }
    steady_clock::time_point assignEndTime = steady_clock::now();
    assignTotalTime += timeDiffInMs(assignStartTime, assignEndTime);

}


void BlockedColumnA::columnAAlgorithm(int argc, char **argv) {

    CSCMatrix fullMatrixA, *localAPencil, *localAPencilTmp;

    localAPencil = new CSCMatrix();
    localAPencilTmp = new CSCMatrix();
    DenseMatrix *localBPencil, *localCPencil;
    long long greaterCount;
    DenseMatrix *receiverCMatrices;

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

    log("scatterAAmongGroups");
    scatterAAmongGroups(fullMatrixA, *localAPencil);


    MPI_Barrier(MPI_COMM_WORLD);
    startTime = steady_clock::now();

    log("replicateA");
    replicateA(*localAPencil);

    const int pencilBCWidth = n / numProcesses;
    const int BCShift = myProcessRank * pencilBCWidth;
    localBPencil = makeDenseMatrix(myProcessRank, numProcesses, n, spec.seed, nBeforeExtending);
    localCPencil = makeDenseMatrix(n, pencilBCWidth, BCShift, 0);

    log("main loop");
    for (int j = 0; j < spec.exponent; j++) {
        for (int i = 0; i < numberOfBlocks; i++) {
            sparseTimesDense(*localAPencil, *localBPencil, *localCPencil);
            if (i == numberOfBlocks - 1 && j == spec.exponent - 1)
                break;
            shift(localAPencil, localAPencilTmp, MPI_COMM_WORLD);
        }
        if (j != spec.exponent - 1) {
            log("assignCMatrixToBMatrix");
            assignCMatrixToBMatrix(localBPencil, localCPencil);
        }
    }

    log("gatherResult");
    gatherStartTime = steady_clock::now();
    if (spec.verbose) {
        receiverCMatrices = gatherResultVerbose(localCPencil);
    } else {
        greaterCount = gatherResultGreater(localCPencil);
    }
    gatherEndTime = steady_clock::now();

    MPI_Barrier(MPI_COMM_WORLD);
    endTime = steady_clock::now();

    log("printResult");
    if (spec.verbose) {
        printResult(receiverCMatrices);
    } else {
        stream() << "resultGather: " << greaterCount << endl;
        if (myProcessRank == 0)
            cout << greaterCount << endl;
    }

    MPI_Finalize();
    log("after finalize");
    log("");
    stream() << "sparseTimesDenseTotal: " << sparseTimeDenseTotalTime << endl;
    stream() << "shiftTotal: " << shiftTotalTime << endl;
    stream() << "assignTotal: " << assignTotalTime<< endl << endl;
    stream() << "replicateATotal: " << timeDiffInMs(replicateStartTime, replicateEndTime) << endl;
    stream() << "gatherTotal: " << timeDiffInMs(gatherStartTime, gatherEndTime) << endl;
    stream() << "totalAlgorithmTime: " << timeDiffInMs(startTime, endTime) << endl;
    delete[] localAPencil->nonzeros;
    delete[] localAPencil->extents;
    delete[] localAPencil->indices;
//    delete[] localAPencil;
//    delete[] localAPencilTmp;
}