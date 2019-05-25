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

    DenseMatrix * gatherResultVerbose(DenseMatrix *localCBlock);

    void printResult(DenseMatrix *fullC);

    void assignCMatrixToBMatrix(DenseMatrix *localBPencil, DenseMatrix *localCBlock);
};


#endif //AC370756_BLOCKEDINNERABC_H
