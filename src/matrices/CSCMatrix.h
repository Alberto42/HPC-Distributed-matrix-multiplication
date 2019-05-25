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

class CSCMatrix{
public:

    vector<CSCMatrix> split(int peacesCount);

    friend ostream &operator<<(ostream &os, const CSCMatrix &matrix);

    CSCMatrix();

    CSCMatrix(double *nonzeros, int *extents, int *indices, int n, int m, int count, int maxNonzeroInRow, int offset,
                 int shiftHorizontal);

    double *nonzeros;
    int *extents, *indices;
    int n, m, count, maxNonzeroInRow, offset, shift;

    void sendSync(int dest, const int *tags);

    void sendAsync(int dest, const int *tags, MPI_Request *req, MPI_Comm comm);

    void receiveSync(int src, const int *tags);
};

CSCMatrix operator>>(istream &stream, CSCMatrix &matrix);


#endif //AC370756_CSCMATRIX_H
