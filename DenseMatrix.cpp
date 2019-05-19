//
// Created by albert on 19.05.19.
//

#include "DenseMatrix.h"
#include "densematgen.h"
#include <cstdlib>
#include "utils.h"


DenseMatrix::DenseMatrix() {}

void DenseMatrix::set(int row, int col, double value) {
    matrix[row * n + col] = value;
}

void DenseMatrix::add(int row, int col, double value) {
    set(row, col, get(row, col) + value);
}

double DenseMatrix::get(const int row, const int col) {
    return matrix[row * n + col];
}

DenseMatrix *makeDenseMatrix(int pencilNumber, int numProcesses, int n, int seed) {
    int columnsInPeace = n / numProcesses;
    DenseMatrix *d = (DenseMatrix *) malloc(sizeof(DenseMatrix) + n * columnsInPeace * sizeof(double));

    int colRangeBegin = pencilNumber * columnsInPeace;
    int colRangeEnd = colRangeBegin + columnsInPeace;

    d->n = n;
    d->m = columnsInPeace;
    d->shift = colRangeBegin;

    for (int row = 0; row < n; row++) {
        for (int col = colRangeBegin; col < colRangeEnd; col++) {
            d->set(row, col, generate_double(seed, row, col));
        }
    }
    return d;
}

DenseMatrix *makeDenseMatrix(int n, int m, int shift) {
    DenseMatrix *d = (DenseMatrix *) malloc(sizeof(DenseMatrix) + n * m * sizeof(double));
    d->n = n;
    d->m = m;
    d->shift = shift;
    for (int i = 0; i < n * m; i++)
        d->matrix[i] = 0;
    return d;
}