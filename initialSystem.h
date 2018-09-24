#ifndef INITIALSYSTEM_H
#define INITIALSYSTEM_H


class InitialSystem
{
public:
	InitialSystem(); //Default Constructor.

  //Get methods:
  float getPi() const;
  float getGamma() const;
  float getEpsilon() const;
  unsigned int getDissip() const;
  unsigned int getNbiter() const;
  unsigned int getRungekutta() const;

  float getMach() const;
  float getAlpha() const;
  float getReynolds() const;
  float getTinf() const;
  float getXref() const;
  float getYref() const;
  float getCmac() const; //cmac?
  float getRhoInfini() const;
  float getUInfini() const;
  float getVInfini() const;
  float getPInfini() const;

  float getClTargert() const; //cltarget?
  float getDcl() const; //dcl?

  float getrms0() const; //rms0?


  //Attributes.
public: 

  float pi_;
  float gamma_;
  float epsilon_; 
  float cmac_;
	/* constants from "input file" */
  int dissip_;
  int nbiter_;
  int rungekutta_;

  /* flow & geometry properties */
	float mach_;
  float alpha_;
  float reynolds_; //from "input file"
  float tinf_;
	float xref_;
  float yref_; //from "input file"
  float rhoInfini_;
  float uInfini_;
  float vInfini_;
  float pInfini_; //rho_free, u_free, v_free, p_free
  
  /*constant cl run */
  float cltarget_;
  float dcl_; //from "input file"

  /*convergence*/
  float rms0_;

  //Sover
  float cfl_;

  //Variables of Runge Kutta
  double rk_alfa[5];
  double rk_beta[5];


  private:

};

#endif