#include <iostream>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <string>

namespace po = boost::program_options;
using namespace std;

struct ProgramSpec {
    string file;
    int seed,c,exponent,g;
    bool verbose,i,m;

    ProgramSpec() : g(-1), verbose(false), i(false), m(false) {}

};
class CSCMatrix {
public:
    double *nonzeros;
    int *extents, *indices;
    int n,m, count,maxNonzeroInRow;
};
CSCMatrix operator>>(istream& stream, CSCMatrix& matrix) {
    stream >> matrix.n >> matrix.m >> matrix.count >> matrix.maxNonzeroInRow;
    double *csrNonzeros;
    int *csrExtends, *csrIndices;
    csrNonzeros = new double[matrix.count];
    csrExtends = new int[matrix.n+1];
    csrIndices = new int[matrix.count];
    matrix.nonzeros = new double[matrix.count];
    matrix.extents = new int[matrix.m+1];
    matrix.indices = new int[matrix.count];

    for(int i=0;i<matrix.count;i++) {
        cin >> csrNonzeros[i];
    }
    for(int i=0;i<matrix.n+1;i++) {
        cin >> csrExtends[i];
    }
    for(int i=0;i<matrix.count;i++) {
        cin >> csrIndices[i];
    }
    vector< tuple<int,int,double> > tmp;
    for(int i=1;i<matrix.n+1;i++) {
        int l=csrExtends[i-1],r=csrExtends[i];
        int rowIdx = i-1;
        for (int j = l;j < r;j++) {
            int colIdx = csrIndices[j];
            double nonzero = csrNonzeros[j];
            tmp.emplace_back(make_tuple(colIdx, rowIdx, nonzero));
        }
    }
    assert((int)tmp.size() == matrix.count);
    for(unsigned i=0;i<tmp.size();i++) {
        matrix.nonzeros[i] = get<2>(tmp[i]);
        matrix.indices[i] = get<0>(tmp[i]);
    }
    int a=0;
    matrix.extents[0] = 0;
    for(int colIdx=0;colIdx < matrix.m;colIdx++) {
        while(get<1>(tmp[a]) == colIdx)
            a++;
        matrix.extents[colIdx+1] = a;
    }

    delete [] csrNonzeros;
    delete [] csrExtends;
    delete [] csrIndices;

    return matrix;
}
void parseArgs(int argc, const char **argv, ProgramSpec &s) {
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
                ("help", "help_message")
                ;

        po::variables_map vm;
        po::store (po::command_line_parser (argc, argv).options (desc)
                                .style (po::command_line_style::default_style |
                                        po::command_line_style::allow_long_disguise)
                                .run (), vm);
        po::notify(vm);
        if (vm.count("help")) {
            std::cout << desc << '\n';
            exit(0);
        }
        s.file = vm["f"].as<string>();
        s.seed = vm["s"].as<int>();
        s.c = vm["c"].as<int>();
        s.exponent = vm["e"].as<int>();
        s.g = vm["g"].as<int>();
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
        string helpMsg = stream.str ();
        boost::algorithm::replace_all (helpMsg, "--", "-");
        cout << helpMsg << endl;
    }
}
int main(int argc,const char **argv) {
    ProgramSpec s;
    parseArgs(argc, argv, s);
    std::cout << "Hello, World!" << std::endl;
}