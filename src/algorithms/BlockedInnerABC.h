//
// Created by albert on 22.05.19.
//

#ifndef AC370756_BLOCKEDINNERABC_H
#define AC370756_BLOCKEDINNERABC_H


#include "MatmulAlgorithm.h"

class BlockedInnerABC : public MatmulAlgorithm{
public:
    MPI_Comm groupDenseReplicate;
    int myRowBlock;
    void innerABCAlgorithm(int argc,char **argv);

    void createMPICommunicators();

    void calcGroups();
};


#endif //AC370756_BLOCKEDINNERABC_H
