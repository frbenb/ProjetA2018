#ifndef LINE_H_
#define LINE_H_

#include "Shape.h"
using namespace std;

class Line : public Shape {
    private:

    public:
        Line(float* p1, float* p2);
        float** get();
        void print();
};

#endif