//
// Created by albert on 22.05.19.
//

#ifndef AC370756_BLOCKEDINNERABC_H
#define AC370756_BLOCKEDINNERABC_H


#include "MatmulAlgorithm.h"
#include "../matrices/CSRMatrix.h"

class BlockedInnerABC : public MatmulAlgorithm{
public:
    MPI_Comm groupDenseReplicate;
    MPI_Comm groupShift;
    int myRowBlock;
    int groupShiftRank;
    void innerABCAlgorithm(int argc,char **argv);

    void createMPICommunicators();

    void calcGroups();

    void sparseTimesDense(const CSRMatrix &A, DenseMatrix &B, DenseMatrix &result);

    void shift(CSRMatrix *&localA, CSRMatrix *&localATmp, MPI_Comm comm);

    DenseMatrix * gatherResultVerbose(DenseMatrix *localCPencil);

    void printResult(DenseMatrix *fullC);

    void assignCMatrixToBMatrix(DenseMatrix *localBPencil, DenseMatrix *localCPencil);

    void gatherResultAfterMultiplicationAndAssign(DenseMatrix *localBPencil, DenseMatrix *localCPencil, MPI_Comm comm);
};


#endif //AC370756_BLOCKEDINNERABC_H
