//
// Created by albert on 23.05.19.
//

#include <algorithm>
#include <sys/stat.h>
#include <fstream>
#include <mpi.h>
#include <string>
#include <boost/algorithm/string/replace.hpp>
#include <boost/program_options.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/options_description.hpp>
#include <iostream>
#include "src/utils.h"
#include "CSCMatrix.h"
#include "SparseMatrix.h"


SparseMatrix::SparseMatrix() { offset = 0; }

SparseMatrix::SparseMatrix(double *nonzeros, int *extents, int *indices, int n, int m, int count, int maxNonzeroInRow,
                     int offset, int shiftHorizontal)
        : nonzeros(nonzeros), extents(extents), indices(indices), n(n), m(m), count(count),
          maxNonzeroInRow(maxNonzeroInRow), offset(offset), shiftHorizontal(shiftHorizontal) {}

void SparseMatrix::sendSync(int dest, const int *tags) {
    MPI_Send((const void *) this, getSize(), MPI_BYTE, dest, tags[0], MPI_COMM_WORLD);
    MPI_Send((const void *) nonzeros, count, MPI_DOUBLE, dest, tags[1], MPI_COMM_WORLD);
    MPI_Send((const void *) extents, m + 1, MPI_INT, dest, tags[2], MPI_COMM_WORLD);
    MPI_Send((const void *) indices, count, MPI_INT, dest, tags[3], MPI_COMM_WORLD);
}

void SparseMatrix::sendAsync(int dest, const int *tags, MPI_Request *req) {
    MPI_Isend((const void *) this, getSize(), MPI_BYTE, dest, tags[0], MPI_COMM_WORLD, req);
    MPI_Isend((const void *) nonzeros, count, MPI_DOUBLE, dest, tags[1], MPI_COMM_WORLD, req + 1);
    MPI_Isend((const void *) extents, m + 1, MPI_INT, dest, tags[2], MPI_COMM_WORLD, req + 2);
    MPI_Isend((const void *) indices, count, MPI_INT, dest, tags[3], MPI_COMM_WORLD, req + 3);
}

void SparseMatrix::receiveSync(int src, const int *tags) {
    MPI_Status status;
    MPI_Recv((void *) this, getSize(), MPI_BYTE, src, tags[0], MPI_COMM_WORLD, &status);
    nonzeros = new double[count];
    extents = new int[m + 1];
    indices = new int[count];
    MPI_Recv((void *) nonzeros, count, MPI_DOUBLE, src, tags[1],
             MPI_COMM_WORLD, &status);
    MPI_Recv((void *) extents, m + 1, MPI_INT, src, tags[2], MPI_COMM_WORLD,
             &status);
    MPI_Recv((void *) indices, count, MPI_INT, src, tags[3], MPI_COMM_WORLD,
             &status);
}

size_t SparseMatrix::getSize() {
    return sizeof(CSCMatrix);
//    return 0;
}
