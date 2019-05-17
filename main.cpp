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

namespace po = boost::program_options;
using namespace std;

struct ProgramSpec {
    string file;
    int seed, c, exponent, g;
    bool verbose, i, m;

    ProgramSpec() : g(-1), verbose(false), i(false), m(false) {}

};

template<class T>
void printArray(T* array,int size, ostream& stream = cout) {
    for(int i=0;i<size;i++) {
        stream<< array[i] << " ";
    }
    stream << endl;
}

class Logger {
public:
    ofstream outfile;

    Logger(int myProcessNo) {
        mkdir("logger", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        string name = "logger/" + to_string(myProcessNo);
        outfile = ofstream(name);
    }

    Logger() {}

    virtual ~Logger() {
        outfile.close();

    }

    void log(const string &s) {
        outfile << s << endl;
    }
};

Logger *logger;

class CSCMatrix {
public:
    double *nonzeros;
    int *extents, *indices;
    int n, m, count, maxNonzeroInRow, offset;

    CSCMatrix() { offset = 0; }

    CSCMatrix(double *nonzeros, int *extents, int *indices, int n, int m, int count, int maxNonzeroInRow, int offset)
            : nonzeros(nonzeros), extents(extents), indices(indices), n(n), m(m), count(count),
              maxNonzeroInRow(maxNonzeroInRow), offset(offset) {}

    vector<CSCMatrix> split(int pencilsCount) {
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
                    first);
            result.push_back(nextMatrix);
        }
        return result;
    }
};

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
void parseArgs(int argc, char **argv, ProgramSpec &s) {
    po::options_description desc{"Options"};
    try {
        desc.add_options()
                ("f", po::value<std::string>(), "Sparse matrix file")
                ("s", po::value<int>(), "Seed")
                ("c", po::value<int>())
                ("e", po::value<int>())
                ("g", po::value<int>())
                ("v", "verbose")
                ("i", "i")
                ("m", "m")
                ("help", "help_message");

        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).options(desc)
                          .style(po::command_line_style::default_style |
                                 po::command_line_style::allow_long_disguise)
                          .run(), vm);
        po::notify(vm);
        if (vm.count("help")) {
            std::cout << desc << '\n';
            exit(0);
        }
        s.file = vm["f"].as<string>();
        s.seed = vm["s"].as<int>();
        s.c = vm["c"].as<int>();
        s.exponent = vm["e"].as<int>();
        if (vm.find("g") != vm.end()) {
            s.g = vm["g"].as<int>();
        }
        if (vm.find("v") != vm.end()) {
            s.verbose = true;
        }
        if (vm.find("i") != vm.end()) {
            s.i = true;
        }
        if (vm.find("m") != vm.end()) {
            s.m = true;
        }
    }
    catch (const po::error &ex) {
        std::cerr << ex.what() << '\n';
    }
    catch (...) {
        std::stringstream stream;
        stream << desc;
        string helpMsg = stream.str();
        boost::algorithm::replace_all(helpMsg, "--", "-");
        cout << helpMsg << endl;
    }
}
int myProcessNo;
int numProcesses;
int groupId, groupsCount, processesPerGroup;
const int INITIAL_SCATTER_TAG1 = 1;
const int INITIAL_SCATTER_TAG2 = 2;
const int INITIAL_SCATTER_TAG3 = 3;

const int INITIAL_SCATTER_TAG4 = 4;

void scatterAAmongGroups(CSCMatrix &fullMatrixA, CSCMatrix &localAColumn) {
    vector<CSCMatrix> pencils;
    if (myProcessNo == 0) {
        pencils = fullMatrixA.split(groupsCount);
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
                    (const void *) &pencils[i].nonzeros,
                    pencils[i].count,
                    MPI_DOUBLE,
                    i,
                    INITIAL_SCATTER_TAG2,
                    MPI_COMM_WORLD
            );
            MPI_Send(
                    (const void *) &pencils[i].extents,
                    pencils[i].m + 1,
                    MPI_INT,
                    i,
                    INITIAL_SCATTER_TAG3,
                    MPI_COMM_WORLD
            );
            MPI_Send(
                    (const void *) &pencils[i].indices,
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
        MPI_Recv((void *) localAColumn.nonzeros, localAColumn.count, MPI_DOUBLE, 0, INITIAL_SCATTER_TAG2,
                 MPI_COMM_WORLD, &status);
        MPI_Recv((void *) localAColumn.extents, localAColumn.m + 1, MPI_INT, 0, INITIAL_SCATTER_TAG3, MPI_COMM_WORLD,
                 &status);
        MPI_Recv((void *) localAColumn.indices, localAColumn.count, MPI_INT, 0, INITIAL_SCATTER_TAG4, MPI_COMM_WORLD,
                 &status);
    }

}
void sparseTimesDense(int argc, char *argv[]) {

    CSCMatrix fullMatrixA, localAColumn;

    ProgramSpec s;
    parseArgs(argc, argv, s);

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myProcessNo);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);

    assert(numProcesses % s.c == 0);
    groupsCount = s.c;
    groupId = myProcessNo % numProcesses;
    processesPerGroup = numProcesses / groupsCount;

    logger = new Logger(myProcessNo);

    if (myProcessNo == 0) {
        ifstream ifs = ifstream(s.file);

        ifs >> fullMatrixA;
    }
    scatterAAmongGroups(fullMatrixA, localAColumn);
    // generate appropriate B submatrices
    // synchronize
    // start timer
    // replicate matrices inside groups
    // multiply, shift
    // gather results
    // print results
}

int main(int argc, char **argv) {
    sparseTimesDense(argc, argv);
}