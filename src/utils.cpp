//
// Created by albert on 17.05.19.
//

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


Logger::Logger(int myProcessNo) {
    mkdir("logger", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    string name = "logger/" + to_string(myProcessNo);
    outfile = ofstream(name);
}

Logger::Logger() {}

Logger::~Logger() {
    outfile.close();

}


ofstream &Logger::stream() {
    return outfile;
}

namespace utils {
    Logger *logger;
}

void initLogger(int myProcessNo) {
    utils::logger = new Logger(myProcessNo);
    utils::logger->startTime = clock();
}

ofstream &stream() {
    return utils::logger->stream();
}
