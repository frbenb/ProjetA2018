#ifndef EFLUX_H
#define EFLUX_H

class Eflux
{
public:
    void eflux();
    void ~eflux();

private:
    int i, j;
    float gamma_;
    float rho_, u_, v_, p_; 

    // float sx,sy,un,qq; 
    // float fl0,fl1,fl2,fl3,fr0,fr1,fr2,fr3;
    
    double **rho_, **u_ ,**v_, **p_; //primitives variables cell centered
    double **residualInviscid_rho_, **residualInviscid_u_;
    double **residualInviscid_v_, **residualInviscid_p_;


    float **flux0,**flux1,**flux2,**flux3;
    float **snx,**sny;

    Mesh *mesh;
    InitialSystem * initialSystem_;
}
#endif