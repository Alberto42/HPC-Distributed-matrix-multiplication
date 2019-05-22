#include <iostream>
#include <string>
#include <mpi.h>
#include <fstream>
#include <sys/stat.h>
#include <algorithm>
#include "src/parseInput.h"
#include "src/algorithms/BlockedColumnA.h"
#include "src/algorithms/BlockedInnerABC.h"

using namespace std;

int main(int argc, char **argv) {
    initSpec(argc, argv);
    if (spec.i) {
        BlockedInnerABC blockedInnerAbc;
        blockedInnerAbc.innerABCAlgorithm(argc, argv);
    } else {
        BlockedColumnA blockedColumnA;
        blockedColumnA.columnAAlgorithm(argc, argv);
    }
}