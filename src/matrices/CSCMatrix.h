//
// Created by albert on 17.05.19.
//

#ifndef AC370756_CSCMATRIX_H
#define AC370756_CSCMATRIX_H

#include <vector>
#include <fstream>
#include <assert.h>
#include <mpi.h>

using namespace std;

class CSCMatrix {
public:
    double *nonzeros;
    int *extents, *indices;
    int n, m, count, maxNonzeroInRow, offset, shift;

    CSCMatrix();

    CSCMatrix(double *nonzeros, int *extents, int *indices, int n, int m, int count, int maxNonzeroInRow, int offset,
              int shift);

    void sendAsync(int dest, const int *tags, MPI_Request *req);

    void sendSync(int dest, const int *tags);

    void receiveSync(int src, const int *tags);

    void receiveAsync(int src, const int *tags, MPI_Request *req, int m, int maxNonzeros);

    vector<CSCMatrix> split(int pencilsCount);

    friend ostream &operator<<(ostream &os, const CSCMatrix &matrix);
};

CSCMatrix operator>>(istream &stream, CSCMatrix &matrix);


#endif //AC370756_CSCMATRIX_H
