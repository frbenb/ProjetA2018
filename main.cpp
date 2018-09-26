#include <string>
#include <iostream>
//#include <dirent.h>

// REMOVE
#include <fstream>
#include <string>

#include "Mesh.h"

#define NDIMS 2 // CHECK move

using namespace std;

int main(){

    double crap1, crap2, crap3, crap7, crap8;
    double crap4, crap5, crap6, crap9, crap10;


    ifstream fpin;
    fpin.open("input");

    fpin >> crap4;
    fpin >> crap4 >> crap1 >> crap5;
    fpin >> crap1 >> crap4 >> crap2 >> crap5 >> crap3 >> crap6 >> crap7 >> crap9;
    fpin >> crap4;
    fpin >> crap4 >> crap5 >> crap6 >> crap9 >> crap10;

    fpin >> crap1 >> crap2 >> crap3 >> crap7 >> crap8;

    cout << crap3 << endl;

    fpin.close();

    return 0;
}