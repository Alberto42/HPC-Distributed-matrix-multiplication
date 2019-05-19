//
// Created by albert on 17.05.19.
//

#include "CSCMatrix.h"
#include "utils.h"
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

CSCMatrix::CSCMatrix() { offset = 0; }

CSCMatrix::CSCMatrix(double *nonzeros, int *extents, int *indices, int n, int m, int count, int maxNonzeroInRow,
                     int offset, int shift)
        : nonzeros(nonzeros), extents(extents), indices(indices), n(n), m(m), count(count),
          maxNonzeroInRow(maxNonzeroInRow), offset(offset), shift(shift) {}

vector<CSCMatrix> CSCMatrix::split(int pencilsCount) {
    assert(m % pencilsCount == 0);
    int columnsInPeace = m / pencilsCount;
    vector<CSCMatrix> result;
    for (int colRangeBegin = 0; colRangeBegin < m; colRangeBegin += columnsInPeace) {
        int colRangeEnd = colRangeBegin + columnsInPeace;
        int first = extents[colRangeBegin], second = extents[colRangeEnd];
        CSCMatrix nextMatrix(
                nonzeros + first,
                extents + colRangeBegin,
                indices + first,
                n,
                columnsInPeace,
                second - first,
                -1,
                first,
                colRangeBegin);
        result.push_back(nextMatrix);
    }
    return result;
}

ostream &operator<<(ostream &os, const CSCMatrix &matrix) {
    os << "nonzeros: " << matrix.nonzeros << " extents: " << matrix.extents << " indices: " << matrix.indices
       << " n: " << matrix.n << " m: " << matrix.m << " count: " << matrix.count << " maxNonzeroInRow: "
       << matrix.maxNonzeroInRow << " offset: " << matrix.offset << endl;
    printArray<double>(matrix.nonzeros, matrix.count, os);
    printArray<int>(matrix.extents, matrix.m + 1, os);
    printArray<int>(matrix.indices, matrix.count, os);
    return os;
}

CSCMatrix operator>>(istream &stream, CSCMatrix &matrix) {
    stream >> matrix.n >> matrix.m >> matrix.count >> matrix.maxNonzeroInRow;
    double *csrNonzeros;
    int *csrExtends, *csrIndices;
    csrNonzeros = new double[matrix.count];
    csrExtends = new int[matrix.n + 1];
    csrIndices = new int[matrix.count];
    matrix.nonzeros = new double[matrix.count];
    matrix.extents = new int[matrix.m + 1];
    matrix.indices = new int[matrix.count];

    for (int i = 0; i < matrix.count; i++) {
        stream >> csrNonzeros[i];
    }
    for (int i = 0; i < matrix.n + 1; i++) {
        stream >> csrExtends[i];
    }
    for (int i = 0; i < matrix.count; i++) {
        stream >> csrIndices[i];
    }
    vector<tuple<int, int, double> > tmp;
    for (int i = 1; i < matrix.n + 1; i++) {
        int l = csrExtends[i - 1], r = csrExtends[i];
        int rowIdx = i - 1;
        for (int j = l; j < r; j++) {
            int colIdx = csrIndices[j];
            double nonzero = csrNonzeros[j];
            tmp.emplace_back(make_tuple(colIdx, rowIdx, nonzero));
        }
    }
    sort(tmp.begin(), tmp.end());
    assert((int) tmp.size() == matrix.count);
    for (unsigned i = 0; i < tmp.size(); i++) {
        matrix.nonzeros[i] = get<2>(tmp[i]);
        matrix.indices[i] = get<1>(tmp[i]);
    }
    int a = 0;
    matrix.extents[0] = 0;
    for (int colIdx = 0; colIdx < matrix.m; colIdx++) {
        while (get<0>(tmp[a]) == colIdx)
            a++;
        matrix.extents[colIdx + 1] = a;
    }

    delete[] csrNonzeros;
    delete[] csrExtends;
    delete[] csrIndices;

    return matrix;
}

void CSCMatrix::sendSync(int dest, const int* tags) {
    MPI_Send( (const void *) this, sizeof(CSCMatrix), MPI_BYTE, dest, tags[0], MPI_COMM_WORLD );
    MPI_Send( (const void *) nonzeros, count, MPI_DOUBLE, dest, tags[1], MPI_COMM_WORLD );
    MPI_Send( (const void *) extents, m + 1, MPI_INT, dest, tags[2], MPI_COMM_WORLD );
    MPI_Send( (const void *) indices, count, MPI_INT, dest, tags[3], MPI_COMM_WORLD );
}

void CSCMatrix::sendAsync(int dest, const int* tags, MPI_Request *req) {
    MPI_Isend( (const void *) this, sizeof(CSCMatrix), MPI_BYTE, dest, tags[0], MPI_COMM_WORLD, req);
    MPI_Isend( (const void *) nonzeros, count, MPI_DOUBLE, dest, tags[1], MPI_COMM_WORLD, req+1);
    MPI_Isend( (const void *) extents, m + 1, MPI_INT, dest, tags[2], MPI_COMM_WORLD, req+2);
    MPI_Isend( (const void *) indices, count, MPI_INT, dest, tags[3], MPI_COMM_WORLD, req+3);
}
void CSCMatrix::receiveSync(int src, const int* tags) {
    MPI_Status status;
    MPI_Recv((void *) this, sizeof(CSCMatrix), MPI_BYTE, src, tags[0], MPI_COMM_WORLD, &status);
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

void CSCMatrix::receiveAsync(int src, const int* tags, MPI_Request *req) {
    MPI_Irecv((void *) this, sizeof(CSCMatrix), MPI_BYTE, src, tags[0], MPI_COMM_WORLD, req);
    nonzeros = new double[count];
    extents = new int[m + 1];
    indices = new int[count];
    MPI_Irecv((void *) nonzeros, count, MPI_DOUBLE, src, tags[1],
             MPI_COMM_WORLD, req+1);
    MPI_Irecv((void *) extents, m + 1, MPI_INT, src, tags[2], MPI_COMM_WORLD,
             req+2);
    MPI_Irecv((void *) indices, count, MPI_INT, src, tags[3], MPI_COMM_WORLD,
             req+3);
}