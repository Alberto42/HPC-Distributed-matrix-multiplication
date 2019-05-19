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
            MPI_Send(
                    (const void *) &pencils[i],
                    sizeof(CSCMatrix),
                    MPI_BYTE,
                    i,
                    INITIAL_SCATTER_TAG1,
                    MPI_COMM_WORLD
            );
            MPI_Send(
                    (const void *) pencils[i].nonzeros,
                    pencils[i].count,
                    MPI_DOUBLE,
                    i,
                    INITIAL_SCATTER_TAG2,
                    MPI_COMM_WORLD
            );
            MPI_Send(
                    (const void *) pencils[i].extents,
                    pencils[i].m + 1,
                    MPI_INT,
                    i,
                    INITIAL_SCATTER_TAG3,
                    MPI_COMM_WORLD
            );
            MPI_Send(
                    (const void *) pencils[i].indices,
                    pencils[i].count,
                    MPI_INT,
                    i,
                    INITIAL_SCATTER_TAG4,
                    MPI_COMM_WORLD
            );
        }
    } else if (myProcessRank < processesPerGroup){
        MPI_Status status;
        MPI_Recv((void *) &localAPencil, sizeof(CSCMatrix), MPI_BYTE, 0, INITIAL_SCATTER_TAG1, MPI_COMM_WORLD, &status);
        localAPencil.nonzeros = new double[localAPencil.count];
        localAPencil.extents = new int[localAPencil.m + 1];
        localAPencil.indices = new int[localAPencil.count];
        MPI_Recv((void *) localAPencil.nonzeros, localAPencil.count, MPI_DOUBLE, 0, INITIAL_SCATTER_TAG2,
                 MPI_COMM_WORLD, &status);
        MPI_Recv((void *) localAPencil.extents, localAPencil.m + 1, MPI_INT, 0, INITIAL_SCATTER_TAG3, MPI_COMM_WORLD,
                 &status);
        MPI_Recv((void *) localAPencil.indices, localAPencil.count, MPI_INT, 0, INITIAL_SCATTER_TAG4, MPI_COMM_WORLD,
                 &status);
    }
    n = localAPencil.n;

}
void init(int argc,char **argv) {
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
    MPI_Bcast((void *)&localAPencil, sizeof(CSCMatrix),MPI_BYTE,0,myGroup);
    MPI_Bcast((void *)localAPencil.nonzeros, localAPencil.count,MPI_DOUBLE,0,myGroup);
    MPI_Bcast((void *)localAPencil.extents, localAPencil.m+1,MPI_INT,0,myGroup);
    MPI_Bcast((void *)localAPencil.indices, localAPencil.count, MPI_INT,0,myGroup);
}
void createMPICommunicators() {
    MPI_Comm_split(MPI_COMM_WORLD,myProcessRank % processesPerGroup,myProcessRank,&myGroup);
}
void sparseTimesDense(const CSCMatrix &A, const DenseMatrix &B, DenseMatrix &result) {
    for(int i=1;i < A.m+1;i++) {
        int extentBegin = A.extents[i-1]-A.offset;
        int extentEnd = A.extents[i]-A.offset;
        int colA = i-1 + A.shift;
        for(int j=extentBegin;j < extentEnd; j++) {
            int rowA = A.indices[j];
            double valueA = A.nonzeros[j];
            int rowB = colA;
            int colBBegin = B.shift;
            int colBEnd = B.shift + B.m;
            for(int colB = colBBegin;colB < colBEnd;colB++) {
                double valueB = B.matrix[rowB][colB];
                result.add(rowA, colB, valueA * valueB);
            }

        }
    }
}
void columnAAlgorithm(int argc, char **argv) {

    CSCMatrix fullMatrixA, localAPencil;
    DenseMatrix localBPencil, localCPencil;
    double startTime, endTime;

    init(argc, argv);
    calcGroups();

    if (myProcessRank == 0) {
        ifstream ifs = ifstream(spec.file);

        ifs >> fullMatrixA;
    }
    scatterAAmongGroups(fullMatrixA, localAPencil);
    localBPencil = DenseMatrix(myProcessRank, numProcesses, n, spec.seed);
    MPI_Barrier(MPI_COMM_WORLD);

    if (myProcessRank == 0) {
        startTime = MPI_Wtime();
    }

    replicateAPencils(localAPencil);

    // multiply, shift
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