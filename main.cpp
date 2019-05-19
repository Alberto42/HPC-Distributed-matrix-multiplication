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

namespace po = boost::program_options;
using namespace std;

int myProcessNo;
int numProcesses;
int groupId, groupsCount, processesPerGroup;

void scatterAAmongGroups(CSCMatrix &fullMatrixA, CSCMatrix &localAColumn) {
    vector<CSCMatrix> pencils;
    if (myProcessNo == 0) {
        pencils = fullMatrixA.split(processesPerGroup);
        localAColumn = pencils[0];
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
    } else {
        MPI_Status status;
        MPI_Recv((void *) &localAColumn, sizeof(CSCMatrix), MPI_BYTE, 0, INITIAL_SCATTER_TAG1, MPI_COMM_WORLD, &status);
        localAColumn.nonzeros = new double[localAColumn.count];
        localAColumn.extents = new int[localAColumn.m + 1];
        localAColumn.indices = new int[localAColumn.count];
        MPI_Recv((void *) localAColumn.nonzeros, localAColumn.count, MPI_DOUBLE, 0, INITIAL_SCATTER_TAG2,
                 MPI_COMM_WORLD, &status);
        MPI_Recv((void *) localAColumn.extents, localAColumn.m + 1, MPI_INT, 0, INITIAL_SCATTER_TAG3, MPI_COMM_WORLD,
                 &status);
        MPI_Recv((void *) localAColumn.indices, localAColumn.count, MPI_INT, 0, INITIAL_SCATTER_TAG4, MPI_COMM_WORLD,
                 &status);
    }

}
void init(int argc,char **argv) {
    initSpec(argc, argv);
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myProcessNo);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
    initLogger(myProcessNo);
}
void calcGroups() {
    assert(numProcesses % spec.c == 0);
    groupsCount = spec.c;
    groupId = myProcessNo % numProcesses;
    processesPerGroup = numProcesses / groupsCount;
}

void sparseTimesDense(int argc, char *argv[]) {

    CSCMatrix fullMatrixA, localAColumn;

    init(argc, argv);
    calcGroups();

    if (myProcessNo == 0) {
        ifstream ifs = ifstream(spec.file);

        ifs >> fullMatrixA;
    }
    scatterAAmongGroups(fullMatrixA, localAColumn);
    log(localAColumn);
    // generate appropriate B submatrices
    // synchronize
    // start timer
    // replicate matrices inside groups
    // multiply, shift
    // gather results
    // print results

    MPI_Finalize();
}

int main(int argc, char **argv) {
    sparseTimesDense(argc, argv);
}