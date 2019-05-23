//
// Created by albert on 22.05.19.
//

#ifndef AC370756_COLUMNA_H
#define AC370756_COLUMNA_H

#include <mpi.h>
#include <src/matrices/CSCMatrix.h>
#include <src/matrices/DenseMatrix.h>
#include "MatmulAlgorithm.h"

using namespace std;


class BlockedColumnA : public MatmulAlgorithm {
    int n, nBeforeExtending;
    int maxANonzeros;
    MPI_Comm myGroup;
public:
    void columnAAlgorithm(int argc, char **argv);

    void scatterAAmongGroups(CSCMatrix &fullMatrixA, CSCMatrix &localAPencil);

    void replicateAPencils(CSCMatrix &localAPencil);

    void createMPICommunicators();

    void sparseTimesDense(const CSCMatrix &A, DenseMatrix &B, DenseMatrix &result);

    void shift(CSCMatrix *&localAPencil, CSCMatrix *&localAPencilTmp);

    DenseMatrix *gatherResultVerbose(DenseMatrix *localCPencil);

    int gatherResultGreater(DenseMatrix *localCPencil);

    void printResult(DenseMatrix *receiverCMatrices);

    void assignCMatrixToBMatrix(DenseMatrix *localBPencil, DenseMatrix *localCPencil);

    void extendA(CSCMatrix *fullMatrixA, int numProcesses);
};


#endif //AC370756_COLUMNA_H
