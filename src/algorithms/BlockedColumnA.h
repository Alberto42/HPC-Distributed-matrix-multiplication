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
public:
    void columnAAlgorithm(int argc, char **argv);

    void createMPICommunicators();

    void calcGroups();

    void sparseTimesDense(const CSCMatrix &A, DenseMatrix &B, DenseMatrix &result);

    DenseMatrix *gatherResultVerbose(DenseMatrix *localCPencil);

    void printResult(DenseMatrix *receiverCMatrices);

    void assignCMatrixToBMatrix(DenseMatrix *localBPencil, DenseMatrix *localCPencil);
};


#endif //AC370756_COLUMNA_H
