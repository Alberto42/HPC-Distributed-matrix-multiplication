//
// Created by albert on 22.05.19.
//

#ifndef AC370756_MATMULALGORITHM_H
#define AC370756_MATMULALGORITHM_H


class MatmulAlgorithm {
public:

    int myProcessRank;
    int numProcesses;
    int numberOfBlocks;
    int n, nBeforeExtending;
    int maxANonzeros;
    MPI_Comm myGroup;

    void init(int argc, char **argv);

    void extendA(CSCMatrix *fullMatrixA, int numProcesses);

    void scatterAAmongGroups(CSCMatrix &fullMatrixA, CSCMatrix &localA);

    void replicateA(CSCMatrix &localA);

    void shift(CSCMatrix *&localAPencil, CSCMatrix *&localAPencilTmp, MPI_Comm comm);

    long long gatherResultGreater(DenseMatrix *localCPencil);
};


#endif //AC370756_MATMULALGORITHM_H
