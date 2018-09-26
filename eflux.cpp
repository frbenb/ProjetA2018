#include "eflux.h"

void Eflux::eflux()
        : {}











g=nsc.gamma;
  
  ro=mesh->ro; Ri_ro=mesh->Ri_ro;
  uu=mesh->uu; Ri_uu=mesh->Ri_uu;
  vv=mesh->vv; Ri_vv=mesh->Ri_vv;
  pp=mesh->pp; Ri_pp=mesh->Ri_pp;

  // Coordonnees du flux dans le mesh
  flux0   = mesh->tmp_ro; // rho
  flux1   = mesh->tmp_uu; // vitesse u
  flux2   = mesh->tmp_vv; // vitesse v
  flux3   = mesh->tmp_pp; // Pression

  // Initialisation des coordonnees du flux???
  for (j=0;j<=mesh->hjmax;j++)
    for (i=0;i<=mesh->himax;i++)
    {
     flux0[i][j]=0.;
     flux1[i][j]=0.;
     flux2[i][j]=0.;
     flux3[i][j]=0.;
    }

  /* i direction */
  
  // C'est quoi snx et sny???
  snx=mesh->six;
  sny=mesh->siy;
  for (j=2;j<=mesh->rjmax;j++)
  {
    for (i=2;i<=mesh->rimax+1;i++)
    {
      r=ro[i-1][j];
      u=uu[i-1][j];
      v=vv[i-1][j];
      p=pp[i-1][j];
      sx=snx[i][j];
      sy=sny[i][j];
      un=u*sx+v*sy;
      qq=u*u+v*v;
      
      fl0=r*un;
      fl1=r*un*u+p*sx;
      fl2=r*un*v+p*sy;
      fl3=r*un*(0.5*qq+g/(g-1.)*p/r);

      r=ro[i][j];
      u=uu[i][j];
      v=vv[i][j];
      p=pp[i][j];
      sx=snx[i][j];
      sy=sny[i][j];
      un=u*sx+v*sy;
      qq=u*u+v*v;
      
      fr0=r*un;
      fr1=r*un*u+p*sx;
      fr2=r*un*v+p*sy;
      fr3=r*un*(0.5*qq+g/(g-1.)*p/r);
      
      flux0[i][j]=0.5*(fl0+fr0);
      flux1[i][j]=0.5*(fl1+fr1);
      flux2[i][j]=0.5*(fl2+fr2);
      flux3[i][j]=0.5*(fl3+fr3);
    }
  }
  for (j=2;j<=mesh->rjmax;j++)
  {
    for (i=2;i<=mesh->rimax+1;i++)
    {      
      Ri_ro[i][j] -= flux0[i][j]; Ri_ro[i-1][j] += flux0[i][j];
      Ri_uu[i][j] -= flux1[i][j]; Ri_uu[i-1][j] += flux1[i][j];
      Ri_vv[i][j] -= flux2[i][j]; Ri_vv[i-1][j] += flux2[i][j];
      Ri_pp[i][j] -= flux3[i][j]; Ri_pp[i-1][j] += flux3[i][j];
    }
  }
  
  /* j direction */
  
  snx=mesh->sjx;
  sny=mesh->sjy;
  for (i=2;i<=mesh->rimax;i++)
  {
    for (j=2;j<=mesh->rjmax+1;j++)
    {
      r=ro[i][j-1];
      u=uu[i][j-1];
      v=vv[i][j-1];
      p=pp[i][j-1];
      sx=snx[i][j];
      sy=sny[i][j];
      un=u*sx+v*sy;
      qq=u*u+v*v;
      
      fl0=r*un;
      fl1=r*un*u+p*sx;
      fl2=r*un*v+p*sy;
      fl3=r*un*(0.5*qq+g/(g-1.)*p/r);

      r=ro[i][j];
      u=uu[i][j];
      v=vv[i][j];
      p=pp[i][j];
      sx=snx[i][j];
      sy=sny[i][j];
      un=u*sx+v*sy;
      qq=u*u+v*v;
      
      fr0=r*un;
      fr1=r*un*u+p*sx;
      fr2=r*un*v+p*sy;
      fr3=r*un*(0.5*qq+g/(g-1.)*p/r);
      
      flux0[i][j]=0.5*(fl0+fr0);
      flux1[i][j]=0.5*(fl1+fr1);
      flux2[i][j]=0.5*(fl2+fr2);
      flux3[i][j]=0.5*(fl3+fr3);
    }
  }
  
  /* j direction BC */
  j=2;
  for (i=2;i<=mesh->rimax;i++)
  {
    r=0.5*(ro[i][j]+ro[i][j-1]);
    u=0.5*(uu[i][j]+uu[i][j-1]);
    v=0.5*(vv[i][j]+vv[i][j-1]);
    p=0.5*(pp[i][j]+pp[i][j-1]);
    sx=snx[i][j];
    sy=sny[i][j];
    un=u*sx+v*sy;
    qq=u*u+v*v;

    fr0=r*un;
    fr1=r*un*u+p*sx;
    fr2=r*un*v+p*sy;
    fr3=r*un*(0.5*qq+g/(g-1.)*p/r);

    flux0[i][j]=fr0;
    flux1[i][j]=fr1;
    flux2[i][j]=fr2;
    flux3[i][j]=fr3;
  }
  j=mesh->rjmax+1;
  for (i=2;i<=mesh->rimax;i++)
  {
    r=0.5*(ro[i][j]+ro[i][j-1]);
    u=0.5*(uu[i][j]+uu[i][j-1]);
    v=0.5*(vv[i][j]+vv[i][j-1]);
    p=0.5*(pp[i][j]+pp[i][j-1]);
    sx=snx[i][j];
    sy=sny[i][j];
    un=u*sx+v*sy;
    qq=u*u+v*v;

    fr0=r*un;
    fr1=r*un*u+p*sx;
    fr2=r*un*v+p*sy;
    fr3=r*un*(0.5*qq+g/(g-1.)*p/r);

    flux0[i][j]=fr0;
    flux1[i][j]=fr1;
    flux2[i][j]=fr2;
    flux3[i][j]=fr3;
  }
  for (i=2;i<=mesh->rimax;i++)
  {
    for (j=2;j<=mesh->rjmax+1;j++)
    {
      Ri_ro[i][j] -= flux0[i][j]; Ri_ro[i][j-1] += flux0[i][j];
      Ri_uu[i][j] -= flux1[i][j]; Ri_uu[i][j-1] += flux1[i][j];
      Ri_vv[i][j] -= flux2[i][j]; Ri_vv[i][j-1] += flux2[i][j];
      Ri_pp[i][j] -= flux3[i][j]; Ri_pp[i][j-1] += flux3[i][j];
    }
  }
  
  for (i=2;i<=mesh->rimax;i++)
  {
    for (j=2;j<=mesh->rjmax;j++)
    {
      Ri_ro[i][j] = Ri_ro[i][j]/mesh->area[i][j];
      Ri_uu[i][j] = Ri_uu[i][j]/mesh->area[i][j];
      Ri_vv[i][j] = Ri_vv[i][j]/mesh->area[i][j];
      Ri_pp[i][j] = Ri_pp[i][j]/mesh->area[i][j];
    }
  }
}
