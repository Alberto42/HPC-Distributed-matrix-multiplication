//
// Created by albert on 19.05.19.
//

#ifndef AC370756_DENSEMATRIX_H
#define AC370756_DENSEMATRIX_H


class DenseMatrix {
public:
    double **matrix;
    int n, m, shift;

    DenseMatrix(int pencilNumber, int numProcesses, int n, int seed);

    DenseMatrix();

    DenseMatrix(int n, int m);

    void add(int row, int col, double value);
};


#endif //AC370756_DENSEMATRIX_H
