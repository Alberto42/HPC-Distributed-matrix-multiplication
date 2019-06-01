//
// Created by albert on 22.05.19.
//

#include <src/matrices/CSRMatrix.h>
#include <src/matrices/DenseMatrix.h>
#include <src/utils.h>
#include <src/parseInput.h>
#include <src/const.h>
#include "BlockedInnerABC.h"
#include "MatmulAlgorithm.h"
#include <cstring>

using namespace std::chrono;

void BlockedInnerABC::innerABCAlgorithm(int argc, char **argv) {
    CSRMatrix fullMatrixA, *localA, *localATmp;
    DenseMatrix *localBPencil, *localCPencil, *fullMatrixC;
    long long greaterCount;

    localA = new CSRMatrix();
    localATmp = new CSRMatrix();

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
    scatterAAmongGroups(fullMatrixA, *localA);

    if (myProcessRank < numberOfBlocks) {
        localBPencil = makeDenseMatrix(myProcessRank, numberOfBlocks, n, spec.seed, nBeforeExtending);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    startTime = steady_clock::now();

    log("replicateA");
    replicateStartTime = steady_clock::now();
    replicateA(*localA);

    const int blockCWidth = n / numberOfBlocks;
    const int CShiftHorizontal = (myProcessRank % numberOfBlocks) * blockCWidth;

    log("replicateB");
    if (myProcessRank >= numberOfBlocks) {
        localBPencil = makeDenseMatrix(n, blockCWidth, CShiftHorizontal, 0);
    }
    MPI_Bcast(localBPencil, localBPencil->size(), MPI_BYTE, 0, groupDenseReplicate);
    replicateEndTime = steady_clock::now();

    localCPencil = makeDenseMatrix(n, blockCWidth, CShiftHorizontal, 0);

    log("main loop");
    for (int j = 0; j < spec.exponent; j++) {
        for (int i = 0; i < numberOfBlocks / spec.c; i++) {
            sparseTimesDense(*localA, *localBPencil, *localCPencil);
            if (i == numberOfBlocks / spec.c - 1 && j == spec.exponent - 1)
                break;
            shift(localA, localATmp, groupShift);
        }
        if (j != spec.exponent - 1)
            gatherResultAfterMultiplicationAndAssign(localBPencil, localCPencil, groupDenseReplicate);
    }

    gatherStartTime = steady_clock::now();
    if (spec.verbose) {
        fullMatrixC = gatherResultVerbose(localCPencil);
    } else {
        gatherResultAfterMultiplicationAndAssign(localBPencil, localCPencil, groupDenseReplicate);
        greaterCount = gatherResultGreater(localBPencil) / spec.c;
    }
    gatherEndTime = steady_clock::now();

    MPI_Barrier(MPI_COMM_WORLD);
    endTime = steady_clock::now();

    log("printResult");
    if (myProcessRank == 0) {
        if (spec.verbose) {
            printResult(fullMatrixC);
        } else {
            cout << greaterCount << endl;
        }
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

}

void BlockedInnerABC::calcGroups() {
    assert(numProcesses % (spec.c * spec.c) == 0);
    numberOfBlocks = numProcesses / spec.c;
    myRowBlock = ((myProcessRank / numberOfBlocks) * (numberOfBlocks / spec.c) + (myProcessRank % numberOfBlocks)) %
                 numberOfBlocks;
}

void BlockedInnerABC::createMPICommunicators() {
    int color = myRowBlock;
    MPI_Comm_split(MPI_COMM_WORLD, color, myProcessRank, &myGroup);
    MPI_Comm_split(MPI_COMM_WORLD, myProcessRank % numberOfBlocks, myProcessRank, &groupDenseReplicate);
    MPI_Comm_split(MPI_COMM_WORLD, myProcessRank / numberOfBlocks, myProcessRank, &groupShift);
    MPI_Comm_rank(groupShift, &groupShiftRank);

}

void BlockedInnerABC::sparseTimesDense(const CSRMatrix &A, DenseMatrix &B, DenseMatrix &result) {
    steady_clock::time_point begin = steady_clock::now();
    for (int i = 1; i < A.m + 1; i++) {
        int extentBegin = A.extents[i - 1] - A.offset;
        int extentEnd = A.extents[i] - A.offset;
        int rowA = i - 1 + A.shift;

        for (int j = extentBegin; j < extentEnd; j++) {
            int colA = A.indices[j];
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

void BlockedInnerABC::shift(CSRMatrix *&localA, CSRMatrix *&localATmp, MPI_Comm comm) {
    steady_clock::time_point begin = steady_clock::now();

    MPI_Request requests[8];
    MPI_Status statuses[8];

    int dest = (groupShiftRank + 1) % numberOfBlocks;
    int src = (groupShiftRank - 1 + numberOfBlocks) % numberOfBlocks;

    localA->CSCMatrix::sendAsync(dest, SHIFT_TAGS, requests, comm);


    MPI_Irecv((void *) localATmp, sizeof(CSCMatrix), MPI_BYTE, src, SHIFT_TAGS[0], comm, requests + 4);

    double *nonzeros = new double[maxANonzeros];
    int *extents = new int[localA->m + 1];
    int *indices = new int[maxANonzeros];

    MPI_Irecv(nonzeros, maxANonzeros, MPI_DOUBLE, src, SHIFT_TAGS[1],
              comm, requests + 5);
    MPI_Irecv(extents, localA->m + 1, MPI_INT, src, SHIFT_TAGS[2], comm,
              requests + 6);
    MPI_Irecv(indices, maxANonzeros, MPI_INT, src, SHIFT_TAGS[3], comm,
              requests + 7);

    MPI_Waitall(8, requests, statuses);

    localATmp->nonzeros = nonzeros;
    localATmp->extents = extents;
    localATmp->indices = indices;

    swap(localA, localATmp);

    delete[] localATmp->nonzeros;
    delete[] localATmp->extents;
    delete[] localATmp->indices;

    steady_clock::time_point end = steady_clock::now();
    shiftTotalTime += timeDiffInMs(begin,end);
}

void BlockedInnerABC::gatherResultAfterMultiplicationAndAssign(DenseMatrix *localBPencil, DenseMatrix *localCPencil,
                                                               MPI_Comm comm) {
    steady_clock::time_point assignStartTime = steady_clock::now();

    assert(localBPencil->n == localCPencil->n);
    assert(localBPencil->m == localCPencil->m);
    MPI_Datatype dtLocalCBlock;
    const size_t localCSize = sizeof(DenseMatrix) + localCPencil->n * localCPencil->m * sizeof(double);
    MPI_Type_contiguous(localCSize, MPI_BYTE, &dtLocalCBlock);
    MPI_Type_commit(&dtLocalCBlock);
    DenseMatrix *receiverCMatrices = nullptr;

    int numProcessesInGroup;
    MPI_Comm_size(comm, &numProcessesInGroup);

    receiverCMatrices = (DenseMatrix *) malloc(numProcessesInGroup * localCSize);

    MPI_Allgather((void *) localCPencil, 1, dtLocalCBlock, receiverCMatrices, 1, dtLocalCBlock, comm);

    for (int i = 0; i < numProcessesInGroup; i++) {
        DenseMatrix *m = getIthMatrix(receiverCMatrices, i);
        for (int localRow = 0; localRow < m->n; localRow++) {
            for (int localCol = 0; localCol < m->m; localCol++) {
                if (i == 0)
                    localBPencil->set(localRow, localCol, 0);
                localBPencil->add(localRow, localCol, m->get(localRow, localCol));
            }
        }
    }

    for (int row = 0; row < localCPencil->n; row++) {
        for (int col = 0; col < localCPencil->m; col++) {
            localCPencil->set(row, col, 0);
        }
    }

    delete[] receiverCMatrices;

    steady_clock::time_point assignEndTime = steady_clock::now();
    assignTotalTime += timeDiffInMs(assignStartTime, assignEndTime);
}

DenseMatrix *BlockedInnerABC::gatherResultVerbose(DenseMatrix *localCPencil) {
    MPI_Datatype dtLocalCBlock;
    const size_t localCSize = sizeof(DenseMatrix) + localCPencil->n * localCPencil->m * sizeof(double);
    MPI_Type_contiguous(localCSize, MPI_BYTE, &dtLocalCBlock);
    MPI_Type_commit(&dtLocalCBlock);
    DenseMatrix *receiverCMatrices = nullptr;
    if (myProcessRank == 0) {
        receiverCMatrices = (DenseMatrix *) malloc(numProcesses * localCSize);
    }
    MPI_Gather((void *) localCPencil, 1, dtLocalCBlock, receiverCMatrices, 1, dtLocalCBlock, 0, MPI_COMM_WORLD);

    DenseMatrix *fullC = nullptr;
    if (myProcessRank == 0) {
        fullC = makeDenseMatrix(n, n, 0, 0);
        for (int i = 0; i < numProcesses; i++) {
            DenseMatrix *m = getIthMatrix(receiverCMatrices, i);
            for (int localRow = 0; localRow < m->n; localRow++) {
                for (int localCol = 0; localCol < m->m; localCol++) {
                    int row = localRow + m->shiftVertical;
                    int col = localCol + m->shiftHorizontal;
                    fullC->add(row, col, m->get(localRow, localCol));
                }
            }
        }
    }
    return fullC;
}

void BlockedInnerABC::printResult(DenseMatrix *fullC) {
    if (myProcessRank == 0) {
        cout << nBeforeExtending << " " << nBeforeExtending << endl;

        for (int row = 0; row < nBeforeExtending; row++) {
            for (int col = 0; col < nBeforeExtending; col++) {
                cout.precision(5);
                cout << "   " << fixed << fullC->get(row, col) << " ";
            }
            cout << endl;
        }
    }
}


