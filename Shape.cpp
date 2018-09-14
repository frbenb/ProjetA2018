#include "Shape.h"

#define NDIMS 2 // CHECK move

Shape::Shape() : points_(nullptr), npoints_(0){

}

Shape::~Shape(){
    for (unsigned int i = 0; i < npoints_; i++)
        delete [] points_[i];

    delete [] points_;
    points_ = nullptr;
}

void Shape::print(){}