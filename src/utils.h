//
// Created by albert on 17.05.19.
//

#ifndef AC370756_UTILS_H
#define AC370756_UTILS_H

#include <iostream>
#include <fstream>

using namespace std;

template<class T>
void printArray(T *array, int size, std::ostream &stream) {
    for (int i = 0; i < size; i++) {
        stream << array[i] << " ";
    }
    stream << std::endl;
}

class Logger {
public:
    ofstream outfile;

    Logger(int myProcessNo);

    Logger();

    virtual ~Logger();

    template<class T>
    void log(T &s);

    ofstream &stream();
};

void initLogger(int myProcessNo);

template<class T>
void Logger::log(T &s) {
    outfile << s << endl;
}

namespace utils {
    extern Logger *logger;
}

template<class T>
void log(T &s) {
    utils::logger->log(s);
}

ofstream &stream();


#endif //AC370756_UTILS_H
