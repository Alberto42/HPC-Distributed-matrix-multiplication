//
// Created by albert on 17.05.19.
//

#ifndef AC370756_PARSEINPUT_H
#define AC370756_PARSEINPUT_H


struct ProgramSpec {
    std::__cxx11::string file;
    int seed, c, exponent, g;
    bool verbose, i, m;

    ProgramSpec();

};

void parseArgs(int argc, char **argv, ProgramSpec &s);


#endif //AC370756_PARSEINPUT_H
