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
#include "utils.h"
#include "parseInput.h"
#include "CSCMatrix.h"
#include "const.h"
#include "DenseMatrix.h"

namespace po = boost::program_options;
using namespace std;

int myProcessRank;
int numProcesses;
int groupId, numberOfGroups, processesPerGroup, sizeOfGroup;
int n;
int maxANonzeros;
MPI_Comm myGroup;

void scatterAAmongGroups(CSCMatrix &fullMatrixA, CSCMatrix &localAPencil) {
    vector<CSCMatrix> pencils;
    if (myProcessRank == 0) {
        pencils = fullMatrixA.split(numberOfGroups);
        localAPencil = pencils[0];
        maxANonzeros = fullMatrixA.count;
        for (int i = 1; i < numberOfGroups; i++) {
            pencils[i].sendSync(i, INITIAL_SCATTER_TAG);
            MPI_Send(&fullMatrixA.count,1,MPI_INT,i,SEND_A_NONZEROS_TAG,MPI_COMM_WORLD);
        }
    } else if (myProcessRank < numberOfGroups) {
        localAPencil.receiveSync(0, INITIAL_SCATTER_TAG);
        MPI_Status status;
        MPI_Recv(&maxANonzeros,1,MPI_INT,0,SEND_A_NONZEROS_TAG,MPI_COMM_WORLD, &status);
    }
    n = localAPencil.n;

}

void init(int argc, char **argv) {
    initSpec(argc, argv);
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myProcessRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
    initLogger(myProcessRank);
}

void calcGroups() {
    assert(numProcesses % spec.c == 0);
    numberOfGroups = numProcesses / spec.c;

    groupId = myProcessRank % numberOfGroups;
    processesPerGroup = spec.c;
}

void replicateAPencils(CSCMatrix &localAPencil) {
    MPI_Bcast((void *) &localAPencil, sizeof(CSCMatrix), MPI_BYTE, 0, myGroup);
    MPI_Bcast((void *) localAPencil.nonzeros, localAPencil.count, MPI_DOUBLE, 0, myGroup);
    MPI_Bcast((void *) localAPencil.extents, localAPencil.m + 1, MPI_INT, 0, myGroup);
    MPI_Bcast((void *) localAPencil.indices, localAPencil.count, MPI_INT, 0, myGroup);
    MPI_Bcast(&maxANonzeros,1,MPI_INT, 0, myGroup);
}

void createMPICommunicators() {
    MPI_Comm_split(MPI_COMM_WORLD, myProcessRank % numberOfGroups, myProcessRank, &myGroup);
}

void sparseTimesDense(const CSCMatrix &A, DenseMatrix &B, DenseMatrix &result) {
    for (int i = 1; i < A.m + 1; i++) {
        int extentBegin = A.extents[i - 1] - A.offset;
        int extentEnd = A.extents[i] - A.offset;
        int colA = i - 1 + A.shift;

        for (int j = extentBegin; j < extentEnd; j++) {
            int rowA = A.indices[j];
            double valueA = A.nonzeros[j];
            int rowB = colA;
            int colBBegin = B.shift;
            int colBEnd = B.shift + B.m;

            for (int colB = colBBegin; colB < colBEnd; colB++) {
                double valueB = B.get(rowB, colB);
                result.add(rowA, colB-colBBegin, valueA * valueB);
            }

        }
    }
}

void shift(CSCMatrix *&localAPencil, CSCMatrix *&localAPencilTmp) {
    MPI_Request requests[8];
    MPI_Status statuses[8];

    localAPencil->sendAsync((myProcessRank + 1) % numProcesses, SHIFT_TAGS, requests);

    int src = (myProcessRank - 1 + numProcesses) % numProcesses;

    MPI_Irecv((void *) localAPencilTmp, sizeof(CSCMatrix), MPI_BYTE, src, SHIFT_TAGS[0], MPI_COMM_WORLD, requests + 4);

    double *nonzeros = new double[maxANonzeros];
    int *extents = new int[localAPencil->m + 1];
    int *indices = new int[maxANonzeros];

    MPI_Irecv(nonzeros, maxANonzeros, MPI_DOUBLE, src, SHIFT_TAGS[1],
              MPI_COMM_WORLD, requests + 5);
    MPI_Irecv(extents, localAPencil->m + 1, MPI_INT, src, SHIFT_TAGS[2], MPI_COMM_WORLD,
              requests + 6);
    MPI_Irecv(indices, maxANonzeros, MPI_INT, src, SHIFT_TAGS[3], MPI_COMM_WORLD,
              requests + 7);

    MPI_Waitall(8, requests, statuses);

    localAPencilTmp->nonzeros = nonzeros;
    localAPencilTmp->extents = extents;
    localAPencilTmp->indices = indices;

    delete[] localAPencil->nonzeros;
    delete[] localAPencil->extents;
    delete[] localAPencil->indices;

    swap(localAPencil, localAPencilTmp);

    delete[] localAPencilTmp->nonzeros;
    delete[] localAPencilTmp->extents;
    delete[] localAPencilTmp->indices;
}

DenseMatrix *gatherResult(DenseMatrix *localCPencil) {
    MPI_Datatype dtDenseMatrix;
    const size_t localCPencilSize = sizeof(DenseMatrix) + n * localCPencil->m * sizeof(double);
    MPI_Type_contiguous(localCPencilSize, MPI_BYTE, &dtDenseMatrix);
    MPI_Type_commit(&dtDenseMatrix);

    DenseMatrix *receiverCMatrices = nullptr;
    if (myProcessRank == 0) {
        receiverCMatrices = (DenseMatrix *) malloc(numProcesses * localCPencilSize);
    }
    MPI_Gather((void *) localCPencil, 1, dtDenseMatrix, receiverCMatrices, 1, dtDenseMatrix, 0, MPI_COMM_WORLD);
    return receiverCMatrices;
}

void printResult(DenseMatrix *receiverCMatrices) {
    if (myProcessRank == 0) {
        cout << n << " " << n << endl;

        for (int row = 0; row < n; row++) {
            for (int i = 0; i < numProcesses; i++) {
                DenseMatrix &m = *getIthMatrix(receiverCMatrices, i);
                for (int colInM = 0; colInM < m.m; colInM++) {
                    cout.precision(5);
                    cout << "   " << fixed << m.get(row, colInM) << " ";
                }
            }
            cout << endl;
        }
    }
}

void columnAAlgorithm(int argc, char **argv) {

    CSCMatrix fullMatrixA, *localAPencil, *localAPencilTmp;

    localAPencil = new CSCMatrix();
    localAPencilTmp = new CSCMatrix();
    DenseMatrix *localBPencil, *localCPencil;
    double startTime, endTime;

    init(argc, argv);
    calcGroups();
    createMPICommunicators();

    log("start algorithm ------------------------------------------------------------");
    if (myProcessRank == 0) {
        ifstream ifs = ifstream(spec.file);

        ifs >> fullMatrixA;
    }

    log("Prepare matrices");
    scatterAAmongGroups(fullMatrixA, *localAPencil);

    n = localAPencil->n;
    const int pencilBCWidth = n / numProcesses;
    const int BCShift = myProcessRank * pencilBCWidth;
    localBPencil = makeDenseMatrix(myProcessRank, numProcesses, n, spec.seed);
    localCPencil = makeDenseMatrix(n, pencilBCWidth, BCShift);
    MPI_Barrier(MPI_COMM_WORLD);

    if (myProcessRank == 0) {
        startTime = MPI_Wtime();
    }

    log("replicateAPencils");
    replicateAPencils(*localAPencil);

    log("main loop");
    for (int i = 0; i < numberOfGroups; i++) {
        log("sparseTimesDense");
        sparseTimesDense(*localAPencil, *localBPencil, *localCPencil);
        if (i == numberOfGroups - 1)
            break;
        log("shift");
        shift(localAPencil, localAPencilTmp);
    }

    log("gatherResult");
    DenseMatrix *receiverCMatrices = gatherResult(localCPencil);

    MPI_Barrier(MPI_COMM_WORLD);
    if (myProcessRank == 0) {
        endTime = MPI_Wtime();
    }
    log("printResult");
    printResult(receiverCMatrices);

    MPI_Finalize();
    log("after finalize");
}

int main(int argc, char **argv) {
    columnAAlgorithm(argc, argv);
}