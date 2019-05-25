//
// Created by albert on 22.05.19.
//

#include <src/matrices/CSRMatrix.h>
#include <src/matrices/DenseMatrix.h>
#include <src/utils.h>
#include <src/parseInput.h>
#include "BlockedInnerABC.h"
#include "MatmulAlgorithm.h"

void BlockedInnerABC::innerABCAlgorithm(int argc,char **argv){
    CSRMatrix fullMatrixA, *localA;
    DenseMatrix *localBPencil, *localCBlock;
    int greaterCount;

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
    scatterAAmongGroups(fullMatrixA, *localA);
    if (myProcessRank < numberOfBlocks) {
        localBPencil = makeDenseMatrix(myProcessRank, numberOfBlocks, n, spec.seed, nBeforeExtending);
    }

    log("replicateB");
    MPI_Bcast(localBPencil, localBPencil->size(), MPI_BYTE, 0, groupDenseReplicate);

    log("replicateA");
    replicateA(*localA);

    const int blockCWidth = n / numberOfBlocks;
    const int blockCLength = blockCWidth * (numberOfBlocks/spec.c);
    const int CShiftHorizontal = (myProcessRank % numberOfBlocks) * blockCWidth;

    localCBlock = makeDenseMatrix(blockCWidth, blockCLength,CShiftHorizontal, myRowBlock*blockCWidth);

    log("main loop");
    for (int j = 0; j < spec.exponent; j++) {
        for (int i = 0; i < numberOfBlocks/spec.c; i++) {
            log("sparseTimesDense");
            sparseTimesDense(*localA, *localBPencil, *localCBlock);
            if (i == numberOfBlocks - 1 && j == spec.exponent - 1)
                break;
            log("shift");
//            shift(localAPencil, localAPencilTmp);
        }
    }


}

void BlockedInnerABC::calcGroups() {
    assert(numProcesses % (spec.c * spec.c) == 0);
    numberOfBlocks = numProcesses / spec.c;
    myRowBlock = (numProcesses % numberOfBlocks)+(numProcesses / numberOfBlocks) * (numberOfBlocks/spec.c);
}

void BlockedInnerABC::createMPICommunicators() {
    int color = myRowBlock;
    MPI_Comm_split(MPI_COMM_WORLD, color, myProcessRank, &myGroup);
    MPI_Comm_split(MPI_COMM_WORLD, myProcessRank % numberOfBlocks, myProcessRank, &groupDenseReplicate);

}

void BlockedInnerABC::sparseTimesDense(const CSRMatrix&A, DenseMatrix &B, DenseMatrix &result) {
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
            int rowABegin = B.shiftVertical;

            for (int colB = colBBegin; colB < colBEnd; colB++) {
                double valueB = B.get(rowB, colB - colBBegin);
                result.add(rowA - rowABegin, colB - colBBegin, valueA * valueB);
            }
        }
    }
}


