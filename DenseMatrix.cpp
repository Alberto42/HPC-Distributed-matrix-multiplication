//
// Created by albert on 19.05.19.
//

#include "DenseMatrix.h"
#include "densematgen.h"

DenseMatrix::DenseMatrix(int pencilNumber, int numProcesses, int n, int seed) {
    int columnsInPeace = m / numProcesses;
    int colRangeBegin = pencilNumber * columnsInPeace;
    int colRangeEnd = colRangeBegin + columnsInPeace;
    matrix = new double *[n];

    for (int row = 0; row < n; row++) {
        matrix[row] = new double[columnsInPeace];
        for (int col = colRangeBegin; col < colRangeEnd; col++) {
            matrix[row][col] = generate_double(seed, row, col);
        }
    }
    this->n = n;
    this->m = columnsInPeace;
    shift = colRangeBegin;

}

DenseMatrix::DenseMatrix() {}

DenseMatrix::DenseMatrix(int n, int m) {
    matrix = new double *[n];

    for (int row = 0; row < n; row++) {
        matrix[row] = new double[m];
        for (int col = 0; col < n; col++) {
            matrix[row][col] = 0;
        }
    }
}

void DenseMatrix::add(int row, int col, double value) {
    matrix[row][col] += value;
}