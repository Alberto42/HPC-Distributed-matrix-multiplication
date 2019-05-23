//
// Created by albert on 22.05.19.
//

#ifndef AC370756_MATMULALGORITHM_H
#define AC370756_MATMULALGORITHM_H


class MatmulAlgorithm {
public:

    int myProcessRank;
    int numProcesses;
    int numberOfGroups;
    int n, nBeforeExtending;
    void init(int argc, char **argv);

    void calcGroups();

    void extendA(CSCMatrix *fullMatrixA, int numProcesses);
};


#endif //AC370756_MATMULALGORITHM_H
