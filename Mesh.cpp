#include "Mesh.h"
#include "InitialSystem.h"
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

#define NDIMS 2 // CHECK move

Mesh::Mesh() : imax_(0), jmax_(0), imaxGhost_(0), jmaxGhost_(0), nbKnots_(0), numberOfCells_(0), 
                rimax_(0), rjmax_(0), inci_(0), incj_(0), x_(nullptr), y_(nullptr), cellArea_(nullptr),
                normal_i_x_(nullptr), normal_i_y_(nullptr), normal_j_x_(nullptr), normal_j_y_(nullptr), 
                rho_(nullptr), u_(nullptr), v_(nullptr), p_(nullptr), rho_nodes_(nullptr), 
                u_nodes_(nullptr), v_nodes_(nullptr), p_nodes_(nullptr), rho_0_(nullptr), u_0_(nullptr), 
                v_0_(nullptr), p_0_(nullptr), residualInviscid_rho_(nullptr), residualInviscid_u_(nullptr),
                residualInviscid_v_(nullptr), residualInviscid_p_(nullptr), residualDissip_rho_(nullptr),
                residualDissip_u_(nullptr), residualDissip_v_(nullptr), residualDissip_p_(nullptr),
                tmp_rho_(nullptr), tmp_u_(nullptr), tmp_v_(nullptr), tmp_p_(nullptr), deltaT_(nullptr),
                speci_(nullptr), specj_(nullptr){}

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
    /*for (unsigned int i = 0; i < nshapes_; i++){
        delete [] shapes_[i];
    }

    shapes_ = nullptr;
    nshapes_ = 0;*/
}

/*void Mesh::print(){

    for (unsigned int i = 0; i < nshapes_; i++){
        shapes_[i]->print();
    }
}*/

void Mesh::write_tecplot(string filename){

    cout << "Saving to : " << filename << endl;

    ofstream tecplot_file;
    tecplot_file.open(filename);

    tecplot_file << "VARIABLES=\"X\",\"Y\",\"RO\",\"U\",\"V\",\"P\"\n";
    tecplot_file << "ZONE T=\"FLOW_FIELD\" i=" << imax_ << " j=" << jmax_ << "\n"; // CHECK change to imax and jmax

    for (int j = 2; j <= jmax_+1; j++){
        for (int i = 2; i <= imax_+1; i++){
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
            x_[i][j]=point_coord/NSC_->getCmac(); // cmac should be defined in NSC
        }
    }
    for (j = 2; j <= rjmax_+1; j++){
        for (i=2; i <= rimax_+1; i++){
            meshfile >> point_coord;
            y_[i][j]=point_coord/NSC_->getCmac(); // cmac should be defined in NSC
        }
    }

    meshfile.close();
}
void Mesh::iterate_pseudo_timestep(int level, int nstage)
{
    int istage;

    //Computation of the timestep.
    timestep();
    save_w0();

    for (istage=0;istage<nstage;istage++)
    {
        // TO CODE Here
        spectral_radius(level);

        residual(level, NSC_->rk_beta[istage], istage, NSC_->dissip_);
        
        update_solution(); // Implementation needed.
        update_boundary(); // Implementation needed.

    }


}

void Mesh::update_solution()
{
    //TBD.
}

void Mesh::update_boundary()
{
    //TBD
}

void Mesh::timestep()
{
    unsigned int i, j;
    double **dt, **area, **speci, **specj;

    dt = deltaT_;
    area = cellArea_;
    speci = speci_;
    specj = specj_;

    for(j=2;j<= rjmax_;j++)
    {
        for(i=2; i<=rimax_;i++)
        {
            
           dt[i][j]= NSC_->cfl_ * area[i][j]/(speci[i][j]+specj[i][j]);
        }

    }
    
}

void Mesh::save_w0()
{
    unsigned int i,j;
    double g, **ro, **uu, **vv, **pp;
    double **ro0, **ru0, **rv0, **re0; // Conservative variables.

    g = NSC_->gamma_;

    ro = rho_;
    uu = u_;
    vv = v_;
    pp = p_;

    ro0 = rho_0_;
    ru0 = u_0_;
    rv0 = v_0_;
    re0 = p_0_;

    for (j=0; j<=jmaxGhost_;j++)
    {
        for(i=0; i<=imaxGhost_;i++)
        {
            ro0[i][j] = ro[i][j];
            ru0[i][j] = ro[i][j]*uu[i][j];
            rv0[i][j] = ro[i][j]*vv[i][j];
            re0[i][j] = 0.5*ro[i][j] * (uu[i][j] * uu[i][j] + vv[i][j] * vv[i][j]) + 1./(g-1.)*pp[i][j];
        }

    }

    return;

}


void Mesh::spectral_radius(int level)
{
    unsigned int i,j;
    double **ro,**uu,**vv,**pp,g,**six,**siy,sx,sy,u_dot_n,cc;

    g=NSC_->gamma_;
  
    ro=rho_;
    uu=u_;
    vv=v_;
    pp=p_;

    /* in i */
    six=normal_i_x_;
    siy=normal_i_y_;

    for (j=1;j<=rjmax_+1;j++)
    {
        for (i=1;i<=rimax_+1;i++)
        {
            sx=0.5*(six[i][j]+six[i+1][j]);
            sy=0.5*(siy[i][j]+siy[i+1][j]);
            u_dot_n=uu[i][j]*sx+vv[i][j]*sy;
            cc=g*pp[i][j]/ro[i][j];

            speci_[i][j]= abs(u_dot_n)+sqrt(cc*(sx*sx+sy*sy));
        }
    }
  
    /* in j */
    six=normal_j_x_;
    siy=normal_j_y_;

    for (j=1;j<=rjmax_+1;j++)
    {
        for (i=1;i<=rimax_+1;i++)
        {
            sx=0.5*(six[i][j]+six[i][j+1]);
            sy=0.5*(siy[i][j]+siy[i][j+1]);
            u_dot_n=uu[i][j]*sx+vv[i][j]*sy;
            cc=g*pp[i][j]/ro[i][j];

            specj_[i][j]=abs(u_dot_n)+sqrt(cc*(sx*sx+sy*sy));
        }
    }

    return;


}

void Mesh::residual(int level, double beta, int istage,int dissip)
{
    unsigned int i,j;

    for (j=0;j<=jmaxGhost_;j++)
    {
        for (i=0;i<=imaxGhost_;i++)
        {
            residualInviscid_rho_[i][j]=0.;
            residualInviscid_u_[i][j]=0.;
            residualInviscid_v_[i][j]=0.;
            residualInviscid_p_[i][j]=0.;
        }
    }
  if (istage==0)
  {
    for (j=0;j<=jmaxGhost_;j++)
    {
      for (i=0;i<=imaxGhost_;i++)
      {
        residualDissip_rho_[i][j]=0.;
        residualDissip_u_[i][j]=0.;
        residualDissip_v_[i][j]=0.;
        residualDissip_p_[i][j]=0.;
      }
    }
  }

  eflux(level);
  
  if(beta>NSC_->epsilon_)
  {
    if (dissip==1)dflux(level,beta);
    if (dissip==2)dflux2(level,beta);
  }
  
  for (j=0;j<=jmaxGhost_;j++)
  {
    for (i=0;i<=imaxGhost_;i++)
    {
        residualInviscid_rho_[i][j]+=residualDissip_rho_[i][j];
        residualInviscid_u_[i][j]+=residualDissip_u_[i][j];
        residualInviscid_v_[i][j]+=residualDissip_v_[i][j];
        residualInviscid_p_[i][j]+=residualDissip_p_[i][j];
    }
  }

  return;

}

void Mesh::eflux(int level)
{
    //TBD.
}

void Mesh::dflux(int level, int beta)
{
    //TBD.

}

void Mesh::dflux2(int level, int beta)
{
    //TBD.

}

