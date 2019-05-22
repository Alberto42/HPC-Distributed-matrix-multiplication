//
// Created by albert on 23.05.19.
//

#ifndef AC370756_SPARSEMATRIX_H
#define AC370756_SPARSEMATRIX_H


#include <vector>
#include <fstream>
#include <assert.h>
#include <mpi.h>

class SparseMatrix {
public:
    SparseMatrix();

    SparseMatrix(double *nonzeros, int *extents, int *indices, int n, int m, int count, int maxNonzeroInRow, int offset,
                 int shiftHorizontal);

    double *nonzeros;
    int *extents, *indices;
    int n, m, count, maxNonzeroInRow, offset, shiftHorizontal;

    size_t getSize();

    void sendSync(int dest, const int *tags);

    void sendAsync(int dest, const int *tags, MPI_Request *req);

    void receiveSync(int src, const int *tags);
};


#endif //AC370756_SPARSEMATRIX_H
