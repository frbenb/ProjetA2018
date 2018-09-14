#include "Quad.h"

// REMOVE
#include <iostream>

using namespace std;

#define NDIMS 2 // CHECK move

Quad::Quad(float* p0, float* p1, float* p2, float* p3){
    this->npoints_ = 4;
    this->points_  = new float*[this->npoints_];

    for (unsigned int i = 0; i < this->npoints_; i++)
        this->points_[i] = new float[NDIMS];

    for (unsigned int i = 0; i < NDIMS; i++)
        this->points_[0][i] = p0[i];
    for (unsigned int i = 0; i < NDIMS; i++)
        this->points_[1][i] = p1[i];
    for (unsigned int i = 0; i < NDIMS; i++)
        this->points_[2][i] = p2[i];
    for (unsigned int i = 0; i < NDIMS; i++)
        this->points_[3][i] = p3[i];
}

float** Quad::get(){
    cout << "gotten" << endl;
    return this->points_;
}

void Quad::print(){
    cout << "Quad with points:" << endl;
    for (unsigned int i = 0; i < this->npoints_; i++){
        for (unsigned int j = 0; j < NDIMS; j++){
            cout << " " << this->points_[i][j];
        }
        cout << endl;
    }
}