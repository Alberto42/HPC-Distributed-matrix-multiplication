//
// Created by albert on 17.05.19.
//

#ifndef AC370756_UTILS_H
#define AC370756_UTILS_H

#include <fstream>
#include <iostream>

using namespace std;
template<class T>
void printArray(T* array,int size, std::ostream& stream) {
    for(int i=0;i<size;i++) {
        stream<< array[i] << " ";
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
    ofstream& stream();
};

template<class T>
void Logger::log(T &s) {
    outfile << s << endl;
}

#endif //AC370756_UTILS_H
