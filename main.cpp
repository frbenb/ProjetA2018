#include <string>
#include <iostream>
//#include <dirent.h>

// REMOVE
#include <fstream>
#include <string>

#include "Mesh.h"

#define NDIMS 2 // CHECK move

using namespace std;

int main(){
    InitialSystem* init = new InitialSystem();

    init->readctrl("input");

    Mesh* mesh = new Mesh(init);
    
    mesh->read_tecplot();

    mesh->write_tecplot("test.x");

    return 0;
}