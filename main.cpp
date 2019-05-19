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
int groupId, groupsCount, processesPerGroup;
int n;
MPI_Comm myGroup;

void scatterAAmongGroups(CSCMatrix &fullMatrixA, CSCMatrix &localAPencil) {
    vector<CSCMatrix> pencils;
    if (myProcessRank == 0) {
        pencils = fullMatrixA.split(processesPerGroup);
        localAPencil = pencils[0];
        for (int i = 1; i < processesPerGroup; i++) {
            pencils[i].send(i,INITIAL_SCATTER_TAG);
        }
    } else if (myProcessRank < processesPerGroup) {
        localAPencil.receive(0,INITIAL_SCATTER_TAG);
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
    groupsCount = spec.c;
    groupId = myProcessRank % numProcesses;
    processesPerGroup = numProcesses / groupsCount;
}

void replicateAPencils(CSCMatrix &localAPencil) {
    MPI_Bcast((void *) &localAPencil, sizeof(CSCMatrix), MPI_BYTE, 0, myGroup);
    MPI_Bcast((void *) localAPencil.nonzeros, localAPencil.count, MPI_DOUBLE, 0, myGroup);
    MPI_Bcast((void *) localAPencil.extents, localAPencil.m + 1, MPI_INT, 0, myGroup);
    MPI_Bcast((void *) localAPencil.indices, localAPencil.count, MPI_INT, 0, myGroup);
}

void createMPICommunicators() {
    MPI_Comm_split(MPI_COMM_WORLD, myProcessRank % processesPerGroup, myProcessRank, &myGroup);
}

void sparseTimesDense(const CSCMatrix &A, const DenseMatrix &B, DenseMatrix &result) {
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
                double valueB = B.matrix[rowB][colB];
                result.add(rowA, colB, valueA * valueB);
            }

        }
    }
}
void columnAAlgorithm(int argc, char **argv) {

    CSCMatrix fullMatrixA, *localAPencil, *localAPencilTmp;
    localAPencil = new CSCMatrix();
    localAPencilTmp = new CSCMatrix();
    DenseMatrix localBPencil, localCPencil;
    double startTime, endTime;

    init(argc, argv);
    calcGroups();

    if (myProcessRank == 0) {
        ifstream ifs = ifstream(spec.file);

        ifs >> fullMatrixA;
    }
    scatterAAmongGroups(fullMatrixA, *localAPencil);
    localBPencil = DenseMatrix(myProcessRank, numProcesses, n, spec.seed);
    MPI_Barrier(MPI_COMM_WORLD);

    if (myProcessRank == 0) {
        startTime = MPI_Wtime();
    }

    replicateAPencils(*localAPencil);

    for(int i=0;i<groupsCount;i++) {
        sparseTimesDense(*localAPencil, localBPencil, localCPencil);
    }

    // gather results
    // print results

    MPI_Barrier(MPI_COMM_WORLD);
    if (myProcessRank == 0) {
        endTime = MPI_Wtime();
    }

    MPI_Finalize();
}

int main(int argc, char **argv) {
    columnAAlgorithm(argc, argv);
}