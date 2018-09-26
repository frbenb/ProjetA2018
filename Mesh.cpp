#include "Mesh.h"
#include "InitialSystem.h"
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

#define NDIMS 2 // CHECK move

Mesh::Mesh(unsigned int imax, unsigned int jmax, unsigned int itl, unsigned int itu, InitialSystem* NSC) : 
                imax_(imax), jmax_(jmax), itl_(itl), itu_(itu), imaxGhost_(imax+2), jmaxGhost_(jmax+2), 
                nbKnots_((imax+3)*(jmax+3)), numberOfCells_((imax+2)*(jmax+2)), 
                rimax_(imax), rjmax_(jmax), inci_(jmax + 3), incj_(1), initSyst_(NSC){

    // himax is imaxGhost_
    // hjmax is jmaxGhost_

    x_ = new double*[imaxGhost_+1];
    y_ = new double*[imaxGhost_+1];
    cellArea_ = new double*[imaxGhost_+1];
    normal_i_x_ = new double*[imaxGhost_+1];
    normal_i_y_ = new double*[imaxGhost_+1];
    normal_j_x_ = new double*[imaxGhost_+1];
    normal_j_y_ = new double*[imaxGhost_+1];
    rho_ = new double*[imaxGhost_+1];
    u_ = new double*[imaxGhost_+1];
    v_ = new double*[imaxGhost_+1];
    p_ = new double*[imaxGhost_+1];
    rho_nodes_ = new double*[imaxGhost_+1];
    u_nodes_ = new double*[imaxGhost_+1];
    v_nodes_ = new double*[imaxGhost_+1];
    p_nodes_ = new double*[imaxGhost_+1];
    rho_0_ = new double*[imaxGhost_+1];
    u_0_ = new double*[imaxGhost_+1];
    v_0_ = new double*[imaxGhost_+1];
    p_0_ = new double*[imaxGhost_+1];
    speci_ = new double*[imaxGhost_+1];
    specj_ = new double*[imaxGhost_+1];
    deltaT_ = new double*[imaxGhost_+1];
    residualInviscid_rho_ = new double*[imaxGhost_+1];
    residualInviscid_u_ = new double*[imaxGhost_+1];
    residualInviscid_v_ = new double*[imaxGhost_+1];
    residualInviscid_p_ = new double*[imaxGhost_+1];
    residualDissip_rho_ = new double*[imaxGhost_+1];
    residualDissip_u_ = new double*[imaxGhost_+1];
    residualDissip_v_ = new double*[imaxGhost_+1];
    residualDissip_p_ = new double*[imaxGhost_+1];
    tmp_rho_ = new double*[imaxGhost_+1];
    tmp_u_ = new double*[imaxGhost_+1];
    tmp_v_ = new double*[imaxGhost_+1];
    tmp_p_ = new double*[imaxGhost_+1];

    for (unsigned int i = 0; i < imaxGhost_+1; i++){
        x_[i] = new double[jmaxGhost_+1];
        y_[i] = new double[jmaxGhost_+1];
        cellArea_[i] = new double[jmaxGhost_+1];
        normal_i_x_[i] = new double[jmaxGhost_+1];
        normal_i_y_[i] = new double[jmaxGhost_+1];
        normal_j_x_[i] = new double[jmaxGhost_+1];
        normal_j_y_[i] = new double[jmaxGhost_+1];
        rho_[i] = new double[jmaxGhost_+1];
        u_[i] = new double[jmaxGhost_+1];
        v_[i] = new double[jmaxGhost_+1];
        p_[i] = new double[jmaxGhost_+1];
        rho_nodes_[i] = new double[jmaxGhost_+1];
        u_nodes_[i] = new double[jmaxGhost_+1];
        v_nodes_[i] = new double[jmaxGhost_+1];
        p_nodes_[i] = new double[jmaxGhost_+1];
        rho_0_[i] = new double[jmaxGhost_+1];
        u_0_[i] = new double[jmaxGhost_+1];
        v_0_[i] = new double[jmaxGhost_+1];
        p_0_[i] = new double[jmaxGhost_+1];
        speci_[i] = new double[jmaxGhost_+1];
        specj_[i] = new double[jmaxGhost_+1];
        deltaT_[i] = new double[jmaxGhost_+1];
        residualInviscid_rho_[i] = new double[jmaxGhost_+1];
        residualInviscid_u_[i] = new double[jmaxGhost_+1];
        residualInviscid_v_[i] = new double[jmaxGhost_+1];
        residualInviscid_p_[i] = new double[jmaxGhost_+1];
        residualDissip_rho_[i] = new double[jmaxGhost_+1];
        residualDissip_u_[i] = new double[jmaxGhost_+1];
        residualDissip_v_[i] = new double[jmaxGhost_+1];
        residualDissip_p_[i] = new double[jmaxGhost_+1];
        tmp_rho_[i] = new double[jmaxGhost_+1];
        tmp_u_[i] = new double[jmaxGhost_+1];
        tmp_v_[i] = new double[jmaxGhost_+1];
        tmp_p_[i] = new double[jmaxGhost_+1];
    }
    // fill other in loop, check that all in big comment above are dealt with, 
    // fill destructor, fix read functions, finish readCtrl, don't forget for read stuff to delete ry.
}

void Mesh::read_su2(string filename){
    /*if (shapes_ != nullptr){
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
    }*/
}

Mesh::~Mesh(){
    unsigned int i;
    if (x_ != nullptr){
        for (i = 0; i < imaxGhost_+1; i++){
            delete [] x_[i];
        }
        delete [] x_;
        x_ = nullptr;
    }
    if (y_ != nullptr){
        for (i = 0; i < imaxGhost_+1; i++){
            delete [] y_[i];
        }
        delete [] y_;
        y_ = nullptr;
    }
    if (cellArea_ != nullptr){
        for (i = 0; i < imaxGhost_+1; i++){
            delete [] cellArea_[i];
        }
        delete [] cellArea_;
        cellArea_ = nullptr;
    }
    if (normal_i_x_ != nullptr){
        for (i = 0; i < imaxGhost_+1; i++){
            delete [] normal_i_x_[i];
        }
        delete [] normal_i_x_;
        normal_i_x_ = nullptr;
    }
    if (normal_i_y_ != nullptr){
        for (i = 0; i < imaxGhost_+1; i++){
            delete [] normal_i_y_[i];
        }
        delete [] normal_i_y_;
        normal_i_y_ = nullptr;
    }
    if (normal_j_x_ != nullptr){
        for (i = 0; i < imaxGhost_+1; i++){
            delete [] normal_j_x_[i];
        }
        delete [] normal_j_x_;
        normal_j_x_ = nullptr;
    }
    if (normal_j_y_ != nullptr){
        for (i = 0; i < imaxGhost_+1; i++){
            delete [] normal_j_y_[i];
        }
        delete [] normal_j_y_;
        normal_j_y_ = nullptr;
    }
    if (rho_ != nullptr){
        for (i = 0; i < imaxGhost_+1; i++){
            delete [] rho_[i];
        }
        delete [] rho_;
        rho_ = nullptr;
    }
    if (u_ != nullptr){
        for (i = 0; i < imaxGhost_+1; i++){
            delete [] u_[i];
        }
        delete [] u_;
        u_ = nullptr;
    }
    if (v_ != nullptr){
        for (i = 0; i < imaxGhost_+1; i++){
            delete [] v_[i];
        }
        delete [] v_;
        v_ = nullptr;
    }
    if (p_ != nullptr){
        for (i = 0; i < imaxGhost_+1; i++){
            delete [] p_[i];
        }
        delete [] p_;
        p_ = nullptr;
    }
    if (rho_nodes_ != nullptr){
        for (i = 0; i < imaxGhost_+1; i++){
            delete [] rho_nodes_[i];
        }
        delete [] rho_nodes_;
        rho_nodes_ = nullptr;
    }
    if (u_nodes_ != nullptr){
        for (i = 0; i < imaxGhost_+1; i++){
            delete [] u_nodes_[i];
        }
        delete [] u_nodes_;
        u_nodes_ = nullptr;
    }
    if (v_nodes_ != nullptr){
        for (i = 0; i < imaxGhost_+1; i++){
            delete [] v_nodes_[i];
        }
        delete [] v_nodes_;
        v_nodes_ = nullptr;
    }
    if (p_nodes_ != nullptr){
        for (i = 0; i < imaxGhost_+1; i++){
            delete [] p_nodes_[i];
        }
        delete [] p_nodes_;
        p_nodes_ = nullptr;
    }
    if (rho_0_ != nullptr){
        for (i = 0; i < imaxGhost_+1; i++){
            delete [] rho_0_[i];
        }
        delete [] rho_0_;
        rho_0_ = nullptr;
    }
    if (u_0_ != nullptr){
        for (i = 0; i < imaxGhost_+1; i++){
            delete [] u_0_[i];
        }
        delete [] u_0_;
        u_0_ = nullptr;
    }
    if (v_0_ != nullptr){
        for (i = 0; i < imaxGhost_+1; i++){
            delete [] v_0_[i];
        }
        delete [] v_0_;
        v_0_ = nullptr;
    }
    if (p_0_ != nullptr){
        for (i = 0; i < imaxGhost_+1; i++){
            delete [] p_0_[i];
        }
        delete [] p_0_;
        p_0_ = nullptr;
    }
    if (speci_ != nullptr){
        for (i = 0; i < imaxGhost_+1; i++){
            delete [] speci_[i];
        }
        delete [] speci_;
        speci_ = nullptr;
    }
    if (specj_ != nullptr){
        for (i = 0; i < imaxGhost_+1; i++){
            delete [] specj_[i];
        }
        delete [] specj_;
        specj_ = nullptr;
    }
    if (deltaT_ != nullptr){
        for (i = 0; i < imaxGhost_+1; i++){
            delete [] deltaT_[i];
        }
        delete [] deltaT_;
        deltaT_ = nullptr;
    }
    if (residualInviscid_rho_ != nullptr){
        for (i = 0; i < imaxGhost_+1; i++){
            delete [] residualInviscid_rho_[i];
        }
        delete [] residualInviscid_rho_;
        residualInviscid_rho_ = nullptr;
    }
    if (residualInviscid_u_ != nullptr){
        for (i = 0; i < imaxGhost_+1; i++){
            delete [] residualInviscid_u_[i];
        }
        delete [] residualInviscid_u_;
        residualInviscid_u_ = nullptr;
    }
    if (residualInviscid_v_ != nullptr){
        for (i = 0; i < imaxGhost_+1; i++){
            delete [] residualInviscid_v_[i];
        }
        delete [] residualInviscid_v_;
        residualInviscid_v_ = nullptr;
    }
    if (residualInviscid_p_ != nullptr){
        for (i = 0; i < imaxGhost_+1; i++){
            delete [] residualInviscid_p_[i];
        }
        delete [] residualInviscid_p_;
        residualInviscid_p_ = nullptr;
    }
    if (residualDissip_rho_ != nullptr){
        for (i = 0; i < imaxGhost_+1; i++){
            delete [] residualDissip_rho_[i];
        }
        delete [] residualDissip_rho_;
        residualDissip_rho_ = nullptr;
    }
    if (residualDissip_u_ != nullptr){
        for (i = 0; i < imaxGhost_+1; i++){
            delete [] residualDissip_u_[i];
        }
        delete [] residualDissip_u_;
        residualDissip_u_ = nullptr;
    }
    if (residualDissip_v_ != nullptr){
        for (i = 0; i < imaxGhost_+1; i++){
            delete [] residualDissip_v_[i];
        }
        delete [] residualDissip_v_;
        residualDissip_v_ = nullptr;
    }
    if (residualDissip_p_ != nullptr){
        for (i = 0; i < imaxGhost_+1; i++){
            delete [] residualDissip_p_[i];
        }
        delete [] residualDissip_p_;
        residualDissip_p_ = nullptr;
    }
    if (tmp_rho_ != nullptr){
        for (i = 0; i < imaxGhost_+1; i++){
            delete [] tmp_rho_[i];
        }
        delete [] tmp_rho_;
        tmp_rho_ = nullptr;
    }
    if (tmp_u_ != nullptr){
        for (i = 0; i < imaxGhost_+1; i++){
            delete [] tmp_u_[i];
        }
        delete [] tmp_u_;
        tmp_u_ = nullptr;
    }
    if (tmp_v_ != nullptr){
        for (i = 0; i < imaxGhost_+1; i++){
            delete [] tmp_v_[i];
        }
        delete [] tmp_v_;
        tmp_v_ = nullptr;
    }
    if (tmp_p_ != nullptr){
        for (i = 0; i < imaxGhost_+1; i++){
            delete [] tmp_p_[i];
        }
        delete [] tmp_p_;
        tmp_p_ = nullptr;
    }
}

void Mesh::write_tecplot(string filename){

    cout << "Saving to : " << filename << endl;

    ofstream tecplot_file;
    tecplot_file.open(filename);

    tecplot_file << "VARIABLES=\"X\",\"Y\",\"RO\",\"U\",\"V\",\"P\"\n";
    tecplot_file << "ZONE T=\"FLOW_FIELD\" i=" << imax_ << " j=" << jmax_ << "\n"; // CHECK change to imax and jmax

    for (unsigned int j = 2; j <= jmax_+1; j++){
        for (unsigned int i = 2; i <= imax_+1; i++){
            tecplot_file << x_[i][j] << " " << y_[i][j] << " " << rho_nodes_[i][j] << " " << u_nodes_[i][j] << " " << v_nodes_[i][j] << " " << p_nodes_[i][j] << "\n";
        }
    }

    tecplot_file.close();
    cout << "Saved to : " << filename << endl;
}

void Mesh::read_tecplot(string filename){
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

    for (j = 2; j <= rjmax_+1; j++){
        for (i = 2; i <= rimax_+1; i++){
            meshfile >> point_coord;
            x_[i][j]=point_coord/initSyst_->getCmac(); // cmac should be defined in NSC
        }
    }
    for (j = 2; j <= rjmax_+1; j++){
        for (i = 2; i <= rimax_+1; i++){
            meshfile >> point_coord;
            y_[i][j]=point_coord/initSyst_->getCmac(); // cmac should be defined in NSC
        }
    }

    meshfile.close();
}

const void Mesh::print(){
    cout << endl;
}