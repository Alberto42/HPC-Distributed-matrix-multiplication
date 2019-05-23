//
// Created by albert on 19.05.19.
//

#include "DenseMatrix.h"
#include "src/densematgen.h"
#include <cstdlib>
#include "src/utils.h"


DenseMatrix::DenseMatrix() {}

void DenseMatrix::set(int row, int col, double value) {
    matrix[row * m + col] = value;
}

void DenseMatrix::add(int row, int col, double value) {
    set(row, col, get(row, col) + value);
}

double DenseMatrix::get(const int row, const int col) {
    return matrix[row * m + col];
}

size_t DenseMatrix::size() {
    return sizeof(DenseMatrix) + n * m * sizeof(double);
}

ostream &operator<<(ostream &os, DenseMatrix &matrix) {
    os << "n: " << matrix.n << " m: " << matrix.m << " shift: " << matrix.shiftHorizontal << endl;

    for(int row=0;row<matrix.n;row++) {
        for(int col=0;col<matrix.m;col++) {
            os << matrix.get(row,col) << " ";
        }
        os << endl;
    }
    return os;
}

DenseMatrix *makeDenseMatrix(int pencilNumber, int numProcesses, int n, int seed, int nBeforeExtending) {
    int columnsInPeace = n / numProcesses;
    DenseMatrix *d = (DenseMatrix *) malloc(sizeof(DenseMatrix) + n * columnsInPeace * sizeof(double));

    int colRangeBegin = pencilNumber * columnsInPeace;
    int colRangeEnd = colRangeBegin + columnsInPeace;

    d->n = n;
    d->m = columnsInPeace;
    d->shiftHorizontal = colRangeBegin;

    for (int row = 0; row < n; row++) {
        for (int col = colRangeBegin; col < colRangeEnd; col++) {
            int colIdxLocal = col - colRangeBegin;
            if (row < nBeforeExtending && colIdxLocal < nBeforeExtending)
                d->set(row, colIdxLocal, generate_double(seed, row, col));
            else d->set(row, colIdxLocal, 0);
        }
    }
    return d;
}

DenseMatrix *makeDenseMatrix(int n, int m, int shiftHorizontal, int shiftVertical) {
    DenseMatrix *d = (DenseMatrix *) malloc(sizeof(DenseMatrix) + n * m * sizeof(double));
    d->n = n;
    d->m = m;
    d->shiftHorizontal = shiftHorizontal;
    d->shiftVertical = shiftVertical;
    for (int i = 0; i < n * m; i++)
        d->matrix[i] = 0;
    return d;
}

DenseMatrix *getIthMatrix(DenseMatrix* first,int i) {
    size_t len = first->size();
    return (DenseMatrix*)( (char*)first + len*i );
}