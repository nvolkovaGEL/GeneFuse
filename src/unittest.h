#ifndef UNITDISTANCE_H
#define UNITDISTANCE_H

#include <stdio.h>
#include <stdlib.h>
#include <string>

using namespace std;

class UnitTest{
public:
    UnitTest();
    void run();
    bool report(bool result, string message);
};

#endif