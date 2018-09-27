/*This class includes the following functions of NSCODE:
  - initial_system
  - initial_flow_parameter
and the input file.*/

#ifndef INITIALSYSTEM_H_
#define INITIALSYSTEM_H_

#define MAX_MGLEVEL 5

#include <string>
#include <math.h>
#include <fstream>

using namespace std;

class InitialSystem
{
public:
	InitialSystem(); //Constructor
  ~InitialSystem();

  void readctrl(string controlFileName);
  void rungeKuttaInit();


  //Get methods:
  double getPi() const;
  double getGamma() const;
  double getEpsilon() const;
  unsigned int getDissip() const;
  unsigned int getNbiter() const;

  double getMach() const;
  double getAlpha() const;
  double getReynolds() const;
  double getTinf() const;
  double getXref() const;
  double getYref() const;
  double getCmac() const; //cmac?
  double getRhoInfini() const;
  double getUInfini() const;
  double getVInfini() const;
  double getPInfini() const;

  double getClTargert() const; //cltarget?
  double getDcl() const; //dcl?

  double getrms0() const; //rms0?

  double getCfl() const;

public:
  unsigned int imax_, jmax_, itl_, itu_;
  string meshfilename_;
  double rk_beta_[5];
  double rk_alpha_[5];
  unsigned int dissip_;

private:
  string ctrlfilename_, title_;
  double pi_, gamma_, epsilon_; 

	/* constants from "input file" */
  unsigned int nbiter_;
  unsigned int niter_[MAX_MGLEVEL]; // Number of iterations per run
  unsigned int rungekutta_[MAX_MGLEVEL];
  unsigned int itccfl_[MAX_MGLEVEL]; //iterate timestep

  /* flow & geometry properties */
	double mach_, alpha_, reynolds_; //from "input file"
  double tinf_;
	double xref_, yref_, cmac_; //from "input file"
  double rhoInfini_, uInfini_, vInfini_, pInfini_; //rho_free, u_free, v_free, p_free
  
  /*constant cl run */
  double cltarget_;
  double dcl_; //from "input file"

  /*convergence*/
  double rms0_;

  double cfl_;
  
  int itertot_;

  FILE* file_cp_;
  FILE* file_conv_;
};

#endif