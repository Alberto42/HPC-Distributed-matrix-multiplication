//
// Created by albert on 17.05.19.
//

#ifndef AC370756_UTILS_H
#define AC370756_UTILS_H

#include <iostream>
#include <fstream>
#include <chrono>

using namespace std;
using namespace std::chrono;

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
    std::chrono::steady_clock::time_point startTime;

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
    std::chrono::steady_clock::time_point current = std::chrono::steady_clock::now();
    outfile << s << " (Time since start: " << std::chrono::duration_cast<std::chrono::microseconds>(current - startTime).count() / 1000.0 << " ms)" << endl;
}

namespace utils {
    extern Logger *logger;
}

template<class T>
void log(T &s) {
    utils::logger->log(s);
}

ofstream &stream();

double timeDiffInMs(steady_clock::time_point a, steady_clock::time_point b);


#endif //AC370756_UTILS_H
