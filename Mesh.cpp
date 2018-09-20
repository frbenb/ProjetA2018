#include "Mesh.h"
#include "Shape.h"
#include <iostream>
#include <fstream>
#include <string>
#include "Triangle.h"
#include "Quad.h"
#include "Line.h"

using namespace std;

#define NDIMS 2 // CHECK move

Mesh::Mesh() : shapes_(nullptr), nshapes_(0){}

void Mesh::read_su2(string filename){
    if (shapes_ != nullptr){
        for (unsigned int i = 0; i < nshapes_; i++){
            delete [] shapes_[i];
        }
        shapes_ = nullptr;
        nshapes_ = 0;
    }

    cout << "Filename: " << filename << endl;

    ifstream meshfile;
    meshfile.open(filename);

    if (!meshfile.is_open()){
        cout << "File \"" << filename << "\" could not be opened." << endl;
        return;
    }

    string token;
    unsigned int Ndime, npoints, shape_id;

    meshfile >> token >> Ndime;
    if ((token.compare("NDIME=") == 0) && (Ndime != NDIMS)){
        cout << "Wrong number of dimentions. Is " << Ndime << " in file, but " << NDIMS << " in code." << endl;
        return;
    }

    meshfile >> token >> npoints;
    if (token.compare("NPOIN=") != 0){
        cout << "Can't fine NPOIN= token.";
        return;
    }

    float** points = new float*[npoints];

    float point_coord;
    unsigned int point_index[8]; // Max number of points

    for (unsigned int points_counter = 0; points_counter < npoints; points_counter++){
        points[points_counter] = new float[NDIMS];
        for (unsigned int dim = 0; dim < NDIMS; dim++){
            meshfile >> point_coord;
            points[points_counter][dim] = point_coord;
        }
    }

    meshfile >> token >> nshapes_;
    if (token.compare("NELEM=") != 0){
        cout << "Can't fine NELEM= token.";
        return;
    }

    shapes_  = new Shape*[nshapes_];
    Shape* shape_toadd;

    for (unsigned int shape_counter = 0; shape_counter < nshapes_; shape_counter++){
        meshfile >> shape_id;
        
        switch (shape_id){
            case 3:
                meshfile >> point_index[0] >> point_index[1];
                shape_toadd = new Line(points[point_index[0]], points[point_index[1]]);
                shapes_[shape_counter] = shape_toadd;
                break;

            case 5:
                meshfile >> point_index[0] >> point_index[1] >> point_index[2];
                shape_toadd = new Triangle(points[point_index[0]], points[point_index[1]], points[point_index[2]]);
                shapes_[shape_counter] = shape_toadd;
                break;

            case 9:
                meshfile >> point_index[0] >> point_index[1] >> point_index[2] >> point_index[3];
                shape_toadd = new Quad(points[point_index[0]], points[point_index[1]], points[point_index[2]], points[point_index[3]]);
                shapes_[shape_counter] = shape_toadd;
                break;

            case 10:
                cout << "10: Tetrahedral not implemented!" << endl;
                break;

            case 12:
                cout << "12: Hexahedral not implemented!" << endl;
                break;

            case 13:
                cout << "13: Prism not implemented!" << endl;
                break;

            case 14:
                cout << "14: Pyramid not implemented!" << endl;
                break;
        }
    }

    meshfile.close();

    for (unsigned int i = 0; i < npoints; i++){
        delete [] points[i];
    }
}

Mesh::~Mesh(){
    for (unsigned int i = 0; i < nshapes_; i++){
        delete [] shapes_[i];
    }

    shapes_ = nullptr;
    nshapes_ = 0;
}

void Mesh::print(){

    for (unsigned int i = 0; i < nshapes_; i++){
        shapes_[i]->print();
    }
}

void write_tecplot(string filename){

    cout << "Saving to : " << filename << endl;

    ofstream tecplot_file;
    tecplot_file.open(filename);

    tecplot_file << "VARIABLES=\"X\",\"Y\",\"RO\",\"U\",\"V\",\"P\"\n";
    tecplot_file << "ZONE T=\"FLOW_FIELD\" i=" << _imax << " j=" << _jmax << "\n"; // CHECK change to imax and jmax

    for (int j = 2; j <= _jmax+1; j++){
        for (int i = 2; i <= _imax+1; i++){
            tecplot_file << _x[i][j] << " " << _y[i][j] << " " << _rocv[i][j] << " " << _uucv[i][j] << " " << _vvcv[i][j] << " " << _ppcv[i][j] << "\n";
        }
    }

    tecplot_file.close();
    cout << "Saved to : " << filename << endl;
}

void read_tecplot(string filename){
    unsigned int nblocks, imax, jmax, i, j;
    double point_coord;

    cout << "Filename: " << filename << endl;

    ifstream meshfile;
    meshfile.open(filename);

    meshfile >> nblocks;
    if (nblocks > 1){
        cout << "Input meshes should only have 1 block" << endl;
    }

    meshfile >> imax >> jmax;

    for (j=2;j<=_rjmax+1;j++){
        for (i=2;i<=_rimax+1;i++){
            meshfile >> point_coord;
            _x[i][j]=point_coord/_cmac; // cmac should be defined
        }
    }
    for (j=2;j<=_rjmax+1;j++){
        for (i=2;i<=_rimax+1;i++){
            meshfile >> point_coord;
            _y[i][j]=point_coord/_cmac;
        }
    }

    meshfile.close();
}