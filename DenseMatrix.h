//
// Created by albert on 19.05.19.
//

#ifndef AC370756_DENSEMATRIX_H
#define AC370756_DENSEMATRIX_H


class DenseMatrix {
public:
    int n, m, shift;
    double matrix[0];

    DenseMatrix();

    void set(int row, int col, double value);

    double get(const int row, const int col);

    void add(int row, int col, double value);
};

DenseMatrix *makeDenseMatrix(int pencilNumber, int numProcesses, int n, int seed);

DenseMatrix *makeDenseMatrix(int n, int m, int shift);


#endif //AC370756_DENSEMATRIX_H
