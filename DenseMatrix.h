//
// Created by albert on 19.05.19.
//

#ifndef AC370756_DENSEMATRIX_H
#define AC370756_DENSEMATRIX_H


class DenseMatrix {
    double **matrix;
    int n,m;
public:
    DenseMatrix(int pencilNumber, int numProcesses, int n, int seed);
    DenseMatrix();
};


#endif //AC370756_DENSEMATRIX_H
