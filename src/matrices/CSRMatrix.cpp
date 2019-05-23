//
// Created by albert on 23.05.19.
//

#include "CSRMatrix.h"
#include "src/utils.h"
#include <iostream>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <string>
#include <mpi.h>
#include <fstream>
#include <sys/stat.h>
#include <algorithm>

using namespace std;


vector<CSRMatrix> CSRMatrix::split(int peacesCount) {
    assert(m % peacesCount == 0);
    int rowsInPeace = m / peacesCount;
    vector<CSRMatrix> result;
    for (int rowRangeBegin = 0; rowRangeBegin < m; rowRangeBegin += rowsInPeace) {
        int rowRangeEnd = rowRangeBegin + rowsInPeace;
        int first = extents[rowRangeBegin], second = extents[rowRangeEnd];
        CSRMatrix nextMatrix(
                nonzeros + first,
                extents + rowRangeBegin,
                indices + first,
                n,
                rowsInPeace,
                second - first,
                -1,
                first,
                rowRangeBegin
                );
        result.push_back(nextMatrix);
    }
    return result;
}

CSRMatrix operator>>(istream &stream, CSRMatrix &matrix) {
    stream >> matrix.n >> matrix.m >> matrix.count >> matrix.maxNonzeroInRow;
    matrix.nonzeros = new double[matrix.count];
    matrix.extents = new int[matrix.m + 1];
    matrix.indices = new int[matrix.count];

    for (int i = 0; i < matrix.count; i++) {
        stream >> matrix.nonzeros[i];
    }
    for (int i = 0; i < matrix.n + 1; i++) {
        stream >> matrix.extents[i];
    }
    for (int i = 0; i < matrix.count; i++) {
        stream >> matrix.indices[i];
    }

    return matrix;
}

CSRMatrix::CSRMatrix() { offset = 0; }

CSRMatrix::CSRMatrix(double *nonzeros, int *extents, int *indices, int n, int m, int count, int maxNonzeroInRow,
                     int offset, int shiftVertical)
        : nonzeros(nonzeros), extents(extents), indices(indices), n(n), m(m), count(count),
          maxNonzeroInRow(maxNonzeroInRow), offset(offset), shiftVertical(shiftVertical) {}

void CSRMatrix::sendSync(int dest, const int *tags) {
    size_t size = sizeof(CSRMatrix);
    MPI_Send((const void *) this, size, MPI_BYTE, dest, tags[0], MPI_COMM_WORLD);
    MPI_Send((const void *) nonzeros, count, MPI_DOUBLE, dest, tags[1], MPI_COMM_WORLD);
    MPI_Send((const void *) extents, m + 1, MPI_INT, dest, tags[2], MPI_COMM_WORLD);
    MPI_Send((const void *) indices, count, MPI_INT, dest, tags[3], MPI_COMM_WORLD);
}

void CSRMatrix::sendAsync(int dest, const int *tags, MPI_Request *req) {
    size_t size = sizeof(CSRMatrix);
    MPI_Isend((const void *) this, size, MPI_BYTE, dest, tags[0], MPI_COMM_WORLD, req);
    MPI_Isend((const void *) nonzeros, count, MPI_DOUBLE, dest, tags[1], MPI_COMM_WORLD, req + 1);
    MPI_Isend((const void *) extents, m + 1, MPI_INT, dest, tags[2], MPI_COMM_WORLD, req + 2);
    MPI_Isend((const void *) indices, count, MPI_INT, dest, tags[3], MPI_COMM_WORLD, req + 3);
}

void CSRMatrix::receiveSync(int src, const int *tags) {
    size_t size = sizeof(CSRMatrix);
    MPI_Status status;
    MPI_Recv((void *) this, size, MPI_BYTE, src, tags[0], MPI_COMM_WORLD, &status);
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

ostream &operator<<(ostream &os, const CSRMatrix &matrix) {
    os << "nonzeros: " << matrix.nonzeros << " extents: " << matrix.extents << " indices: " << matrix.indices << " n: "
       << matrix.n << " m: " << matrix.m << " count: " << matrix.count << " maxNonzeroInRow: " << matrix.maxNonzeroInRow
       << " offset: " << matrix.offset <<  " shiftVertical: "
       << matrix.shiftVertical;
    printArray<double>(matrix.nonzeros, matrix.count, os);
    printArray<int>(matrix.extents, matrix.m + 1, os);
    printArray<int>(matrix.indices, matrix.count, os);
    return os;
}