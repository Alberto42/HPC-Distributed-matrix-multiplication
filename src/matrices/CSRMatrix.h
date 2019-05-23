//
// Created by albert on 23.05.19.
//

#ifndef AC370756_CSRMATRIX_H
#define AC370756_CSRMATRIX_H


#include <cstddef>
#include <vector>
#include <fstream>
#include <assert.h>
#include <mpi.h>
#include <ostream>
#include "CSCMatrix.h"

using namespace std;
class CSRMatrix : public CSCMatrix {
public:

    vector<CSRMatrix> split(int peacesCount);

    friend ostream &operator<<(ostream &os, const CSRMatrix &matrix);

    CSRMatrix();

    CSRMatrix(double *nonzeros, int *extents, int *indices, int n, int m, int count, int maxNonzeroInRow, int offset,
            int shiftVertical);


    double *nonzeros;
    int *extents, *indices;
    int n, m, count, maxNonzeroInRow, offset, shiftVertical;

    void sendSync(int dest, const int *tags);

    void sendAsync(int dest, const int *tags, MPI_Request *req);

    void receiveSync(int src, const int *tags);
};

CSRMatrix operator>>(istream &stream, CSRMatrix &matrix);


#endif //AC370756_CSRMATRIX_H
