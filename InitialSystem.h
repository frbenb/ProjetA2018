#ifndef INITIALSYSTEM_H_
#define INITIALSYSTEM_H_

#define MAX_MGLEVEL 5

#include <string>

using namespace std;

class InitialSystem
{
public:
	InitialSystem(); //Constructor
  ~InitialSystem();

  void add_NSC_reference();

  void readctrl(string controlFileName);

  //Get methods:
  float getPi() const;
  float getGamma() const;
  float getEpsilon() const;
  unsigned int getDissip() const;
  unsigned int getNbiter() const;

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

  unsigned int imax_, jmax_, itl_, itu_;
  string meshfilename_;

private:
  string ctrlfilename_, title_;
  float pi_, gamma_, epsilon_; 

	/* constants from "input file" */
  unsigned int dissip_, nbiter_;
  unsigned int niter_[MAX_MGLEVEL]; // Number of iterations per run
  unsigned int rungekutta_[MAX_MGLEVEL];
  unsigned int itccfl_[MAX_MGLEVEL]; //iterate timestep

  /* flow & geometry properties */
	float mach_, alpha_, reynolds_; //from "input file"
  float tinf_;
	float xref_, yref_, cmac_; //from "input file"
  float rhoInfini_, uInfini_, vInfini_, pInfini_; //rho_free, u_free, v_free, p_free
  
  /*constant cl run */
  float cltarget_;
  float dcl_; //from "input file"

  /*convergence*/
  float rms0_;
};

#endif