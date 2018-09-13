#include "Triangle.h"

Triangle::Triangle(float* p0, float* p1, float* p2){
    npoints_ = 3;
    float** points_ = new float*[npoints_];

    for (unsigned int i = 0; i < npoints_; i++)
        points_[i] = new float[3];

    for (unsigned int i = 0; i < 3; i++)
        points_[0][i] = p0[i];
    for (unsigned int i = 0; i < 3; i++)
        points_[1][i] = p1[i];
    for (unsigned int i = 0; i < 3; i++)
        points_[2][i] = p2[i];
}

float** Triangle::get(){

}

void Triangle::print(){
    
}