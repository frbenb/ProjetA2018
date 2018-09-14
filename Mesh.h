#ifndef MESH_H_
#define MESH_H_

#include "Shape.h"
#include <string>

using namespace std;

class Mesh {
    private:
        Shape** shapes_; // CHECK if should be single ptr
        unsigned int nshapes_;

    public:
        Mesh();
        ~Mesh(); // CHECK make virtual
        //virtual float** get() const = 0;
        void print();
        void read_su2(string filename);
        void save_to_tecplot(string filename) const;
};

#endif