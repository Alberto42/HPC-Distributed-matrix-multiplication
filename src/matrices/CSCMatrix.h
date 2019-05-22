//
// Created by albert on 17.05.19.
//

#ifndef AC370756_CSCMATRIX_H
#define AC370756_CSCMATRIX_H

#include <vector>
#include <fstream>
#include <assert.h>
#include <mpi.h>
#include "SparseMatrix.h"

using namespace std;

class CSCMatrix : public SparseMatrix{
public:

    using SparseMatrix::SparseMatrix;

    vector<CSCMatrix> split(int pencilsCount);
    size_t getSize();

    friend ostream &operator<<(ostream &os, const CSCMatrix &matrix);
};

CSCMatrix operator>>(istream &stream, CSCMatrix &matrix);


#endif //AC370756_CSCMATRIX_H
