#ifndef SHAPE
#define SHAPE

using namespace std;

class Shape {
    public: // CHECK change to protected
        float ** points_;
        unsigned int npoints_;

    //public:
        Shape();
        ~Shape(); // CHECK make virtual
        //virtual float** get() const = 0;
};

#endif