//
// Created by albert on 22.05.19.
//

#ifndef AC370756_MATMULALGORITHM_H
#define AC370756_MATMULALGORITHM_H

#include <chrono>

using namespace std::chrono;

class MatmulAlgorithm {
public:

    int myProcessRank;
    int numProcesses;
    int numberOfBlocks;
    int n, nBeforeExtending;
    int maxANonzeros;
    MPI_Comm myGroup;
    double sparseTimeDenseTotalTime, shiftTotalTime;
    steady_clock::time_point startTime, endTime, replicateStartTime, replicateEndTime, gatherStartTime, gatherEndTime;

    void init(int argc, char **argv);

    void extendA(CSCMatrix *fullMatrixA, int numProcesses);

    void scatterAAmongGroups(CSCMatrix &fullMatrixA, CSCMatrix &localA);

    void replicateA(CSCMatrix &localA);

    void shift(CSCMatrix *&localAPencil, CSCMatrix *&localAPencilTmp, MPI_Comm comm);

    long long gatherResultGreater(DenseMatrix *localCPencil);
};


#endif //AC370756_MATMULALGORITHM_H
