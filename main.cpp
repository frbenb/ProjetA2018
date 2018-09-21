

#include <string>
#include <iostream>
//#include <dirent.h>

#include "Mesh.h"

#define NDIMS 2 // CHECK move

using namespace std;

int main(){

    Mesh* mesh = new Mesh();
    mesh->read_tecplot("./naca0012_129x129_1B_JAMESON.x");
    //mesh->print();

    return 0;
}