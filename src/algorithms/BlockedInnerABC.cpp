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
    CSRMatrix fullMatrixA, localA;
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
    scatterAAmongGroups(fullMatrixA, localA);
    if (myProcessRank < numberOfBlocks) {
        localBPencil = makeDenseMatrix(myProcessRank, numberOfBlocks, n, spec.seed, nBeforeExtending);
    }

    log("replicateB");
    MPI_Bcast(localBPencil, localBPencil->size(), MPI_BYTE, 0, groupDenseReplicate);

    log("replicateA");
    replicateA(localA);

    const int blockCWidth = n / numberOfBlocks;
    const int CShiftHorizontal = (myProcessRank % numberOfBlocks) * blockCWidth;

    localCBlock = makeDenseMatrix(blockCWidth, blockCWidth,n,CShiftHorizontal, myRowBlock*blockCWidth);


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


