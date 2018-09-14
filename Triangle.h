#ifndef TRIANGLE_H_
#define TRIANGLE_H_

#include "Shape.h"
using namespace std;

class Triangle : public Shape {
    private:

    public:
        Triangle(float* p0, float* p1, float* p2);
        float** get();
        void print();
};

#endif