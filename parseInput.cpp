//
// Created by albert on 17.05.19.
//

#include "utils.h"
#include <algorithm>
#include <sys/stat.h>
#include <fstream>
#include <mpi.h>
#include <string>
#include <boost/algorithm/string/replace.hpp>
#include <boost/program_options.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/options_description.hpp>
#include <iostream>
#include "parseInput.h"

ProgramSpec::ProgramSpec() : g(-1), verbose(false), i(false), m(false) {}

void parseArgs(int argc, char **argv, ProgramSpec &s) {
    boost::program_options::options_description desc{"Options"};
    try {
        desc.add_options()
                ("f", boost::program_options::value<std::string>(), "Sparse matrix file")
                ("s", boost::program_options::value<int>(), "Seed")
                ("c", boost::program_options::value<int>())
                ("e", boost::program_options::value<int>())
                ("g", boost::program_options::value<int>())
                ("v", "verbose")
                ("i", "i")
                ("m", "m")
                ("help", "help_message");

        boost::program_options::variables_map vm;
        boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(desc)
                          .style(boost::program_options::command_line_style::default_style |
                                 boost::program_options::command_line_style::allow_long_disguise)
                          .run(), vm);
        boost::program_options::notify(vm);
        if (vm.count("help")) {
            std::cout << desc << '\n';
            exit(0);
        }
        s.file = vm["f"].as<std::__cxx11::string>();
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
    catch (const boost::program_options::error &ex) {
        std::cerr << ex.what() << '\n';
    }
    catch (...) {
        std::stringstream stream;
        stream << desc;
        std::__cxx11::string helpMsg = stream.str();
        boost::algorithm::replace_all(helpMsg, "--", "-");
        std::cout << helpMsg << std::endl;
    }
}