#include "Line.h"

// REMOVE
#include <iostream>

using namespace std;

#define NDIMS 2 // CHECK move

Line::Line(float* p0, float* p1){
    this->npoints_ = 2;
    this->points_  = new float*[this->npoints_];

    for (unsigned int i = 0; i < this->npoints_; i++)
        this->points_[i] = new float[NDIMS];

    for (unsigned int i = 0; i < NDIMS; i++)
        this->points_[0][i] = p0[i];
    for (unsigned int i = 0; i < NDIMS; i++)
        this->points_[1][i] = p1[i];
}

float** Line::get(){
    cout << "gotten" << endl;
    return this->points_;
}

void Line::print(){
    cout << "Line with points: " << endl;
    for (unsigned int i = 0; i < this->npoints_; i++){
        for (unsigned int j = 0; j < NDIMS; j++){
            cout << " " << this->points_[i][j];
        }
        cout << endl;
    }
}