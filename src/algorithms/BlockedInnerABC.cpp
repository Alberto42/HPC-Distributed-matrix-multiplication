//
// Created by albert on 22.05.19.
//

#include <src/matrices/CSRMatrix.h>
#include <src/matrices/DenseMatrix.h>
#include "BlockedInnerABC.h"
#include "MatmulAlgorithm.h"

void BlockedInnerABC::innerABCAlgorithm(int argc,char **argv){
    CSRMatrix fullMatrixA;
    DenseMatrix *localBPencil, *localCBlock;
    int greaterCount;

    init(argc, argv);

}