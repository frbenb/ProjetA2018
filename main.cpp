

#include <string>
#include <iostream>
//#include <dirent.h>

#include "Shape.h"
#include "Triangle.h"

#define NDIMS 2 // CHECK move

using namespace std;

int main(){
    cout << "Hello world!" << endl;

    float* p0 = new float[NDIMS];
    float* p1 = new float[NDIMS];
    float* p2 = new float[NDIMS];

    cout << "new points" << endl;

    p0[0] = 1;
    p0[1] = 2;

    p1[0] = 3;
    p1[1] = 4;

    p2[0] = 5;
    p2[1] = 6;

    cout << "points filled" << endl;

    Triangle* tri = new Triangle(p0, p1, p2);

    cout << "triangle created" << endl;

    tri->print();

    cout << "printed" << endl;

    return 0;
}