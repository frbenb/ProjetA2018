#ifndef QUAD_H_
#define QUAD_H_

#include "Shape.h"
using namespace std;

class Quad : public Shape {
    private:

    public:
        Quad(float* p0, float* p1, float* p2, float* p3);
        float** get();
        void print();
};

#endif