#ifndef SHAPE_H_
#define SHAPE_H_

using namespace std;

class Shape {
    protected: // CHECK change to protected
        float ** points_;
        unsigned int npoints_;

    public:
        Shape();
        ~Shape(); // CHECK make virtual
        virtual void print() = 0;
        //virtual float** get() const = 0;
};

#endif