#include "Mesh.h"
#include "InitialSystem.h"
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

#define NDIMS 2 // CHECK move

Mesh::Mesh(InitialSystem* NSC) : 
                imax_(NSC->imax_), jmax_(NSC->jmax_), itl_(NSC->itl_), itu_(NSC->itu_), imaxGhost_(NSC->imax_+2), jmaxGhost_(NSC->jmax_+2), 
                nbKnots_((NSC->imax_+3)*(NSC->jmax_+3)), numberOfCells_((NSC->imax_+2)*(NSC->jmax_+2)), 
                rimax_(NSC->imax_), rjmax_(NSC->jmax_), inci_(NSC->jmax_ + 3), incj_(1), NSC_(NSC){

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

void Mesh::read_tecplot(){
    unsigned int nblocks, imax, jmax, i, j;
    double point_coord;
    string filename = NSC_->meshfilename_;

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
        for (i = 2; i <= rimax_+1; i++){
            meshfile >> point_coord;
            y_[i][j]=point_coord/NSC_->getCmac(); // cmac should be defined in NSC
        }
    }

    meshfile.close();
}

void Mesh::iterate_pseudo_timestep(int nstage)
{
    int istage;

    //Computation of the timestep.
    timestep();
    save_w0();

    for (istage=0;istage<nstage;istage++)
    {
        // TO CODE Here
        spectral_radius();

        residual(NSC_->rk_beta[istage], istage, NSC_->dissip_);
        
        update_solution(NSC_->rk_alpha[istage]); // Implementation needed.
        update_boundary(); // Implementation needed.

    }


}

void Mesh::update_solution(float alfa)
{
  int i,j;
  double g,ronew,runew,rvnew,renew,**ro,**uu,**vv,**pp;
  double **ro0,**ru0,**rv0,**re0,**dt;
  double **Ri_ro,**Ri_uu,**Ri_vv,**Ri_pp;
  
  g=NSC_->gamma_;

  ro = rho_; ro0=rho_0_; Ri_ro=residualInviscid_rho_;
  uu= u_; ru0=u_0_; Ri_uu=residualInviscid_u_;
  vv= v_; rv0=v_0_; Ri_vv=residualInviscid_v_;
  pp= p_; re0=p_0_; Ri_pp=residualInviscid_p_; 
  dt=deltaT_;

  for (j=2;j<=rjmax_;j++)
  {
    for (i=2;i<=rimax_;i++)
    {
       ronew=ro0[i][j]-alfa*dt[i][j]*Ri_ro[i][j];
       runew=ru0[i][j]-alfa*dt[i][j]*Ri_uu[i][j];
       rvnew=rv0[i][j]-alfa*dt[i][j]*Ri_vv[i][j];
       renew=re0[i][j]-alfa*dt[i][j]*Ri_pp[i][j];
       
       ro[i][j]=ronew;
       uu[i][j]=runew/ronew;
       vv[i][j]=rvnew/ronew;
       pp[i][j]=(g-1.)*(renew-0.5*(runew*runew+rvnew*rvnew)/ronew);
    }
  }
}

void Mesh::update_boundary()
{
    int i,j,himax,hjmax,rimax;
    double **ro,**uu,**vv,**pp,**sx,**sy,robc,uubc,vvbc,ppbc;
    double ro1,uu1,vv1,pp1,ssx,ssy,ss,un1;
    double g,gm1,cfree,chav_in,el,R4e,R4f,R4,chav_out,R5e,R5f,R5,
            unbc,ccbc,dun,uubc_inlet,vvbc_inlet,ssbc_inlet,
	       uubc_outlet,vvbc_outlet,ssbc_outlet,ela,elb,cc1,unf,ssbc,cc2;
 

  himax=imaxGhost_;
  hjmax=jmaxGhost_;
  rimax=rimax_;

  ro=rho_;
  uu=u_;
  vv=v_;
  pp=p_;
  
  g=NSC_->gamma_;
  gm1=g-1.;
  cfree=sqrt(g*NSC_->pInfini_/NSC_->rhoInfini_);

  /* jmin halos: wall*/
  sx=normal_j_x_;
  sy=normal_j_y_;
  for (i=2;i<=rimax;i++)
  {
    ro1=ro[i][2];
    uu1=uu[i][2];
    vv1=vv[i][2];
    pp1=pp[i][2];
    ssx=sx[i][2];
    ssy=sy[i][2];
    ss=sqrt(ssx*ssx+ssy*ssy);
    ssx/=ss; ssy/=ss;
    ssx*=-1.; ssy*=-1.;
    un1=uu1*ssx+vv1*ssy;

    robc=ro1;
    uubc=uu1-un1*ssx;
    vvbc=vv1-un1*ssy;
    ppbc=pp1;


    ro[i][1]=robc;
    uu[i][1]=2.*uubc-uu1;
    vv[i][1]=2.*vvbc-vv1;
    pp[i][1]=ppbc;

    ro[i][0]=robc;
    uu[i][0]=uu[i][1];
    vv[i][0]=vv[i][1];
    pp[i][0]=pp[i][1];

  }
  /*jmax halos: far-field */
  
  for (i=2;i<=rimax;i++)
  {

    ro1=ro[i][hjmax-2];
    uu1=uu[i][hjmax-2];
    vv1=vv[i][hjmax-2];
    pp1=pp[i][hjmax-2];
    cc1=sqrt(g*pp1/ro1);
    
    ssx=sx[i][hjmax-1];
    ssy=sy[i][hjmax-1];
    ss=sqrt(ssx*ssx+ssy*ssy);
    ssx/=ss; ssy/=ss;
    ssx*=-1.; ssy*=-1.;
    un1=uu1*ssx+vv1*ssy;
    unf=NSC_->uInfini_*ssx+NSC_->vInfini_*ssy;
    chav_in=unf+cfree;
    el=sign(chav_in);
    R4e=un1+2.*cc1/gm1;
    R4f=unf+2.*cfree/gm1;
    R4=0.5*((1+el)*R4f+(1.-el)*R4e);
    chav_out=un1-cc1;
    el=sign(chav_out);
    R5e=un1-2.*cc1/gm1;
    R5f=unf-2.*cfree/gm1;
    R5=0.5*((1+el)*R5f+(1.-el)*R5e);
    unbc=0.5*(R4+R5);
    ccbc=0.25*(R4-R5)*gm1;
    el=sign(unbc);
    dun=unbc-unf;
    uubc_inlet=NSC_->uInfini_+dun*ssx;
    vvbc_inlet=NSC_->vInfini_+dun*ssy;
    ssbc_inlet=NSC_->pInfini_/pow(NSC_->rhoInfini_,g);
    dun=unbc-un1;
    uubc_outlet=uu1+dun*ssx;
    vvbc_outlet=vv1+dun*ssy;
    ssbc_outlet=pp1/pow(ro1,g);
    
    ela=0.5*(1.+el);
    elb=0.5*(1.-el);
    uubc=ela*uubc_inlet+elb*uubc_outlet;
    vvbc=ela*vvbc_inlet+elb*vvbc_outlet;
    ssbc=ela*ssbc_inlet+elb*ssbc_outlet;
    cc2=ccbc*ccbc;
    robc=cc2/g/ssbc;
    robc=pow(robc,1./gm1);
    ppbc=robc*cc2/g;

    ro[i][hjmax-1]=2.*robc- ro1;
    uu[i][hjmax-1]=2.*uubc- uu1;
    vv[i][hjmax-1]=2.*vvbc- vv1;
    pp[i][hjmax-1]=2.*ppbc- pp1;

    ro[i][hjmax]=2.*ro[i][hjmax-1]- ro1;
    uu[i][hjmax]=2.*uu[i][hjmax-1]- uu1;
    vv[i][hjmax]=2.*vv[i][hjmax-1]- vv1;
    pp[i][hjmax]=2.*pp[i][hjmax-1]- pp1;

  }

  /* imin and imax halos: wall+connecting BC*/
  for (j=0;j<=hjmax;j++)
  {
    ro[      0][j]=ro[rimax-1][j]; pp[      0][j]=pp[rimax-1][j];
    ro[      1][j]=ro[rimax  ][j]; pp[      1][j]=pp[rimax  ][j];
    ro[himax-1][j]=ro[      2][j]; pp[himax-1][j]=pp[	   2][j];
    ro[himax  ][j]=ro[      3][j]; pp[himax  ][j]=pp[	   3][j];
    
    uu[      0][j]=uu[rimax-1][j]; vv[      0][j]=vv[rimax-1][j];
    uu[      1][j]=uu[rimax  ][j]; vv[      1][j]=vv[rimax  ][j];
    uu[himax-1][j]=uu[      2][j]; vv[himax-1][j]=vv[	   2][j];
    uu[himax  ][j]=uu[      3][j]; vv[himax  ][j]=vv[	   3][j];
  }
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


void Mesh::spectral_radius()
{
    unsigned int i,j;
    double **ro,**uu,**vv,**pp,g,**six,**siy,sx,sy,u_dot_n,cc;

    g=NSC_->getGamma();
  
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

void Mesh::residual(double beta, int istage,int dissip)
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

  eflux();
  
  if(beta>NSC_->epsilon_)
  {
    if (dissip==1)dflux(beta);
    if (dissip==2)dflux2(beta);
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


void Mesh::monitor_convergence()
{
    int i,j,rimax,rjmax;
    double **pp,**sx,**sy,rms,ppbc,cpbc,cl,cd,dynhead,cmac,alpha,clwind,cdwind;

    rimax=rimax_;
    rjmax=rjmax_;

    pp=p_;

    /* jmin halos: wall*/
    sx=normal_j_x_;
    sy=normal_j_y_;
    cmac=NSC_->cmac_/NSC_->cmac_;
    dynhead=0.5*NSC_->gamma_*NSC_->mach_*NSC_->mach_;
    alpha=NSC_->alpha_*NSC_->pi_/180.;
    cl=0; cd=0.; 
    j=2;

    for (i=2;i<=rimax;i++)
    {
        ppbc=0.5*(pp[i][j]+pp[i][j-1]);
        cpbc=(ppbc-1.)/dynhead;

        if (NSC_->itertot_ == NSC_->nbiter_ - 1) fprintf(NSC_->file_cp_,
            "%f %f %f\n",x_[i][j]*NSC_->cmac_,y_[i][j]*NSC_->cmac_,cpbc);
    cl += -ppbc*sy[i][j];
    cd += -ppbc*sx[i][j];
    
  }
  cl=cl/(dynhead*cmac);
  cd=cd/(dynhead*cmac);
  
  clwind=cl*cos(alpha) - cd*sin(alpha);
  cdwind=cl*sin(alpha) + cd*cos(alpha);

  rms=0.;
  for (j=2;j<=rjmax;j++){
    for (i=2;i<=rimax;i++){
       rms+= residualInviscid_rho_[i][j]*residualInviscid_rho_[i][j];}}

  rms=1./((rjmax-1)*(rimax-1))*sqrt(rms);
  if (NSC_->itertot_==0) NSC_->rms0_=rms;

  printf("%d   %f   %f  %f\n",NSC_->itertot_,log10(rms)-log10(NSC_->rms0_),clwind,cdwind);
  fprintf(NSC_->file_conv_,"%d %f %f %f\n",NSC_->itertot_,log10(rms)-log10(NSC_->rms0_),clwind,cdwind);

}

void Mesh::initial_flow_parameters()
{
  double alpha,mach,g;
    
  printf("in initial_flow_parameters..........................................\n");

  g = NSC_->gamma_;

  alpha=NSC_->alpha_ * NSC_->pi_/180;
  mach=NSC_->mach_;
  
  NSC_->rhoInfini_=1.0;
  NSC_->uInfini_=mach*sqrt(g)*cos(alpha);
  NSC_->vInfini_=mach*sqrt(g)*sin(alpha);
  NSC_->pInfini_=1.0;

  printf("in initial_flow_parameters..........................................DONE\n");
}

void Mesh::initial_system()
{
  printf("in initial_system..........................................\n");
  NSC_->pi_=4.*atan(1.);
  NSC_->gamma_=1.4;
  NSC_->epsilon_=1.0e-28;
  NSC_->file_conv_=fopen("conv","w");
  NSC_->file_cp_=fopen("cp","w");
  printf("in initial_system..........................................DONE\n");
}

void Mesh::mesh4halos()
{
  int i,j,himax,hjmax;
  double **x,**y;
  

  printf("in mesh4halos..........................................\n");
  x=x;
  y=y;
  himax=imaxGhost_;
  hjmax=jmaxGhost_;

  /* imin and imax halos */
  
  for (j=2;j<hjmax;j++)
  {
    x[    0][j]=x[himax-3][j]; y[    0][j]=y[himax-3][j];
    x[    1][j]=x[himax-2][j]; y[    1][j]=y[himax-2][j];
    x[himax][j]=x[      3][j]; y[himax][j]=y[	   3][j];
  }
  /* jmin & jmax halos */
  
  for (i=0;i<=himax;i++)
  {
    x[i][    0]=x[i][      2]; y[i][	0]=y[i][      2];
    x[i][    1]=x[i][      2]; y[i][	1]=y[i][      2];
    x[i][hjmax]=x[i][hjmax-1]; y[i][hjmax]=y[i][hjmax-1];
  }
  
  printf("in mesh4halos..........................................DONE\n");
  return;
}


void Mesh::initial_field()
{
  int i,j;
  double **ro,**uu,**vv,**pp;
  

  printf("in initial_field..........................................\n");

  ro=rho_;
  uu=u_;
  vv=v_;
  pp=p_;

  for (j=0;j<=jmaxGhost_;j++)
  {
    for (i=0;i<=imaxGhost_;i++)
    {
      ro[i][j]=NSC_->rhoInfini_;
      uu[i][j]=NSC_->uInfini_;
      vv[i][j]=NSC_->vInfini_;
      pp[i][j]=NSC_->pInfini_;
    }
  }
  printf("in initial_field..........................................DONE\n");
}

void Mesh::eflux()
{
    unsigned int i, j;
    double gamma_; //Taken from NSC object
    double rhoTemp_, uTemp_, vTemp_, pTemp_; //Temporary variables
    double leftFlux_rho_, leftFlux_u_, leftFlux_v_, leftFlux_p_; //Flux on cell's left
    double rightFlux_rho, rightFlux_u_, rightFlux_v_, rightFlux_p_; //Flux on cell's right
    double sx_, sy_, un_, qq_; //???

    gamma_ = NSC_->getGamma();

    for (j=0; j<=jmax_; j++)
    {
        for (i=0; i<=imax_; i++)
        {
            tmp_rho_[i][j] = 0;
            tmp_u_[i][j] = 0;
            tmp_v_[i][j] = 0;
            tmp_p_[i][j] = 0;
        }
    }

    /* i direction */
    for(j = 2; j<=rjmax_; j++)
    {
        for(i=2; i<=rimax_+1; i++)
        {
            rhoTemp_ = rho_[i-1][j];
            uTemp_ = u_[i-1][j];
            uTemp_ = v_[i-1][j];
            pTemp_ = p_[i-1][j];
            sx_ = normal_i_x_[i][j];
            sy_ = normal_i_y_[i][j];
            un_ = uTemp_ * sx_ + uTemp_ * sy_;
            qq_ = uTemp_ * uTemp_ + uTemp_ * uTemp_;
        
            leftFlux_rho_ = rhoTemp_ * un_;
            leftFlux_u_ = rhoTemp_ * un_ * uTemp_ + pTemp_ * sx_;
            leftFlux_v_ = rhoTemp_ * un_ * vTemp_ + pTemp_ * sy_;
            leftFlux_p_ = rhoTemp_ * un_ * (0.5 * qq_ + gamma_/(gamma_ - 1.) * pTemp_/rhoTemp_);

            rhoTemp_ = rho_[i][j];
            uTemp_ = u_[i][j];
            vTemp_ = v_[i][j];
            pTemp_ = p_[i][j];
            sx_ = normal_i_x_[i][j];
            sy_ = normal_i_y_[i][j];
            un_ = uTemp_ * sx_ + vTemp_ * sy_;
            qq_ = uTemp_ * uTemp_ + vTemp_ * vTemp_;
            
            rightFlux_rho = rhoTemp_ * un_;
            rightFlux_u_ = rhoTemp_ * un_ * uTemp_ + pTemp_ * sx_;
            rightFlux_v_ = rhoTemp_ * un_ * vTemp_ + pTemp_ * sy_;
            rightFlux_p_ = rhoTemp_ * un_ * (0.5 * qq_ + gamma_/(gamma_ - 1.) * pTemp_/rhoTemp_);
            
            tmp_rho_[i][j] = 0.5 * (leftFlux_rho_ + rightFlux_rho);
            tmp_u_[i][j] = 0.5 * (leftFlux_u_ + rightFlux_u_);
            tmp_v_[i][j] = 0.5 * (leftFlux_v_ + rightFlux_v_);
            tmp_p_[i][j] = 0.5 * (leftFlux_p_ + rightFlux_p_);
        }
    }

    for(j=2; j<=rjmax_; j++)
    {
        for(i=2; i<=rimax_; i++)
        {
            residualInviscid_rho_[i][j] -= tmp_rho_[i][j]; 
            residualInviscid_rho_[i-1][j] += tmp_rho_[i][j];

            residualInviscid_u_[i][j] -= tmp_u_[i][j]; 
            residualInviscid_u_[i-1][j] += tmp_u_[i][j];

            residualInviscid_v_[i][j] -= tmp_v_[i][j]; 
            residualInviscid_v_[i-1][j] += tmp_v_[i][j];

            residualInviscid_p_[i][j] -= tmp_p_[i][j]; 
            residualInviscid_p_[i-1][j] += tmp_p_[i][j];
        }
    }

    /* j direction */
    for(i=2; i<=rimax_; i++)
    {
        for(j=2; j<=rjmax_+1; j++)
        {
            rhoTemp_ = rho_[i][j-1];
            uTemp_ = u_[i][j-1];
            vTemp_ = v_[i][j-1];
            pTemp_ = p_[i][j-1];
            sx_ = normal_j_x_[i][j];
            sy_ = normal_j_y_[i][j];
            un_ = uTemp_ * sx_ + vTemp_ * sy_;
            qq_ = uTemp_ * uTemp_ + vTemp_ * vTemp_;
            
            leftFlux_rho_ = rhoTemp_ * un_;
            leftFlux_u_ = rhoTemp_ * un_ * uTemp_ + pTemp_ * sx_;
            leftFlux_v_ = rhoTemp_ * un_ * vTemp_ + pTemp_ * sy_;
            leftFlux_p_ = rhoTemp_ * un_ * (0.5 * qq_ + gamma_/(gamma_ - 1.) * pTemp_/rhoTemp_);

            rhoTemp_ = rho_[i][j];
            uTemp_ = u_[i][j];
            vTemp_ = v_[i][j];
            pTemp_ = p_[i][j];
            sx_ = normal_j_x_[i][j];
            sy_ = normal_j_y_[i][j];
            un_ = uTemp_ * sx_ + vTemp_ * sy_;
            qq_ = uTemp_ * uTemp_ + vTemp_ * vTemp_;
            
            rightFlux_rho = rhoTemp_ * un_;
            rightFlux_u_ = rhoTemp_ * un_ * uTemp_ + pTemp_ * sx_;
            rightFlux_v_ = rhoTemp_ * un_ * vTemp_ + pTemp_ * sy_;
            rightFlux_p_ = rhoTemp_ * un_ * (0.5 * qq_ + gamma_/(gamma_ - 1.) * pTemp_/rhoTemp_);
            
            tmp_rho_[i][j] = 0.5 * (leftFlux_rho_ + rightFlux_rho);
            tmp_u_[i][j] = 0.5 * (leftFlux_u_ + rightFlux_u_);
            tmp_v_[i][j] = 0.5 * (leftFlux_v_ + rightFlux_v_);
            tmp_p_[i][j]= 0.5 * (leftFlux_p_ + rightFlux_p_);           
        }
    }

    /* j direction BC */
    j = 2;
    for(i=2; i<=rimax_; i++)
    {
        rhoTemp_ = 0.5 * (rho_[i][j] + rho_[i][j-1]);
        uTemp_ = 0.5 * (u_[i][j] + u_[i][j-1]);
        vTemp_ = 0.5 * (v_[i][j] + v_[i][j-1]);
        pTemp_ = 0.5 * (p_[i][j] + p_[i][j-1]);
        sx_ = normal_j_x_[i][j];
        sy_ = normal_j_y_[i][j];
        un_ = uTemp_ * sx_ + vTemp_ * sy_;
        qq_ = uTemp_ * uTemp_ + vTemp_ * vTemp_;

        rightFlux_rho = rhoTemp_ * un_;
        rightFlux_u_ = rhoTemp_ * un_ * uTemp_ + pTemp_ * sx_;
        rightFlux_v_ = rhoTemp_ * un_ * vTemp_ + pTemp_ * sy_;
        rightFlux_p_ = rhoTemp_ * un_ * (0.5 * qq_ + gamma_/(gamma_ - 1.) * pTemp_/rhoTemp_);

        tmp_rho_[i][j] = rightFlux_rho;
        tmp_u_[i][j] = rightFlux_u_;
        tmp_v_[i][j] = rightFlux_v_;
        tmp_p_[i][j] = rightFlux_p_;
    }

  j = rjmax_ + 1;
  for(i=2; i<=rimax_; i++)
  {
    rhoTemp_ = 0.5 * (rho_[i][j] + rho_[i][j-1]);
    uTemp_ = 0.5 * (u_[i][j] + u_[i][j-1]);
    vTemp_ = 0.5 * (v_[i][j] + v_[i][j-1]);
    pTemp_ = 0.5 * (p_[i][j] + p_[i][j-1]);
    sx_ = normal_j_x_[i][j];
    sy_ = normal_j_y_[i][j];
    un_ = uTemp_ * sx_ + vTemp_ * sy_;
    qq_ = uTemp_ * uTemp_ + vTemp_ * vTemp_;

    rightFlux_rho = rhoTemp_ * un_;
    rightFlux_u_ = rhoTemp_ * un_ * uTemp_ + pTemp_ * sx_;
    rightFlux_v_ = rhoTemp_ * un_ * vTemp_ + pTemp_ * sy_;
    rightFlux_p_ = rhoTemp_ * un_ * (0.5 * qq_ + gamma_/(gamma_ - 1.) * pTemp_/rhoTemp_);

    tmp_rho_[i][j] = rightFlux_rho;
    tmp_u_[i][j] = rightFlux_u_;
    tmp_v_[i][j] = rightFlux_v_;
    tmp_p_[i][j] = rightFlux_p_;
  }

  for(i=2; i<=rimax_; i++)
  {
    for (j=2; j<=rjmax_ + 1; j++)
    {
      residualInviscid_rho_[i][j] -= tmp_rho_[i][j]; 
      residualInviscid_rho_[i][j-1] += tmp_rho_[i][j];

      residualInviscid_u_[i][j] -= tmp_u_[i][j]; 
      residualInviscid_u_[i][j-1] += tmp_u_[i][j];
      
      residualInviscid_v_[i][j] -= tmp_v_[i][j]; 
      residualInviscid_v_[i][j-1] += tmp_v_[i][j];

      residualInviscid_p_[i][j] -= tmp_p_[i][j]; 
      residualInviscid_p_[i][j-1] += tmp_p_[i][j];
    }
  }
  
  for(i=2; i<=rimax_; i++)
  {
    for(j=2; j<=rjmax_; j++)
    {
      residualInviscid_rho_[i][j] = residualInviscid_rho_[i][j] / cellArea_[i][j];
      residualInviscid_u_[i][j] = residualInviscid_u_[i][j] / cellArea_[i][j];
      residualInviscid_v_[i][j] = residualInviscid_v_[i][j] / cellArea_[i][j];
      residualInviscid_p_[i][j] = residualInviscid_p_[i][j] / cellArea_[i][j];
    }
  }
}

void Mesh::metric()
{
  int i,j;
  double **x,**y,x1,y1,x2,y2;
  

  printf("in metric..........................................\n");
  
  /* area */
  for (j=1;j<=rjmax_+1;j++)
  {
    for (i=1;i<=rimax_+1;i++)
    {
      x1=x[i+1][j+1]-x[i][j];
      y1=y[i+1][j+1]-y[i][j];
      x2=x[i][j+1]-x[i+1][j];
      y2=y[i][j+1]-y[i+1][j];
      cellArea_[i][j]=abs(0.5*(x1*y2-x2*y1));
    }
  }
  
  /* i direction */
  for (j=1;j<=rjmax_+1;j++)
  {
    for (i=1;i<=rimax_+2;i++)
    {
      normal_i_x_[i][j]=   y[  i][j+1]-y[i][j];
      normal_i_y_[i][j]= -(x[  i][j+1]-x[i][j]);
    }
  }
  
  /* j direction */
  for (j=1;j<=rjmax_+2;j++)
  {
    for (i=1;i<=rimax_+1;i++)
    {
      normal_j_x_[i][j]=  (y[i][j  ]-y[i+1][j]);
      normal_j_y_[i][j]= -(x[i][j  ]-x[i+1][j]);
    }
  }

  printf("in metric..........................................DONE\n");
}


void Mesh::transCC_to_CV()
{

  int i,j,inci,incj,ia0,ia1,ia2,ia3;
  double *ro,*uu,*vv,*pp,*rocv,*uucv,*vvcv,*ppcv;

  printf("in transCC_to_CV..........................................\n");
  
  inci=inci_;
  incj=incj_;

  ro=&(rho_[0][0]); rocv=&(rho_nodes_[0][0]);
  uu=&(u_[0][0]); uucv=&(u_nodes_[0][0]);
  vv=&(v_[0][0]); vvcv=&(v_nodes_[0][0]);
  pp=&(p_[0][0]); ppcv=&(p_nodes_[0][0]);

  for (j=2;j<=rjmax_+1;j++)
  {
    for (i=2;i<=rimax_+1;i++)
    {
       ia0= i   *inci+ j   *incj;
       ia1=(i-1)*inci+ j   *incj;
       ia2=(i-1)*inci+(j-1)*incj;
       ia3= i   *inci+(j-1)*incj;
       
       rocv[ia0]=0.25*(ro[ia0]+ro[ia1]+ro[ia2]+ro[ia3]);
       uucv[ia0]=0.25*(uu[ia0]+uu[ia1]+uu[ia2]+uu[ia3]);
       vvcv[ia0]=0.25*(vv[ia0]+vv[ia1]+vv[ia2]+vv[ia3]);
       ppcv[ia0]=0.25*(pp[ia0]+pp[ia1]+pp[ia2]+pp[ia3]);
    
    }
  printf("in transCC_to_CV..........................................DONE\n");
  return;
}

}

void Mesh::tridiagonal(int il,int iu,double *b,double *d,double *a, double *c)
{
  int lp,i,j;
  double r;
  
  /* taken from Anderson, Tannehill & Pletcher */

  lp=il+1;  
  for (i=lp;i<=iu;i++)
  {
    r=b[i]/d[i-1];
    d[i]=d[i]-r*a[i-1];
    c[i]=c[i]-r*c[i-1];
  }
  
  c[iu]=c[iu]/d[iu];
  for (i=lp;i<=iu;i++)
  {
    j=iu-i+il;
    c[j]=(c[j]-a[j]*c[j+1])/d[j];
  }
  
  return;
  
}

void Mesh::dflux( int beta)
{
    //TBD.

}

void Mesh::dflux2(int beta)
{
    //TBD.

}
