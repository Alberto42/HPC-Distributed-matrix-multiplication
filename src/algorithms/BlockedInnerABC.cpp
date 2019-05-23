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

    log("replicateAPencils");



}

void BlockedInnerABC::createMPICommunicators() {

}

void BlockedInnerABC::replicateAPencils(CSCMatrix &localAPencil) {
    MPI_Bcast((void *) &localAPencil, sizeof(CSCMatrix), MPI_BYTE, 0, myGroup);
    if (myProcessRank >= numberOfBlocks) {
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
