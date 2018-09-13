#ifndef TRIANGLE
#define TRIANGLE

#include "Shape.h"
using namespace std;

class Triangle : public Shape {
    private:

    public:
        Triangle(float* p1, float* p2, float* p3);
        float** get();
        void print();
};

#endif