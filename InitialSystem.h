/*This class includes the following functions of NSCODE:
  - initial_system
  - initial_flow_parameter
and the input file.*/

#ifndef INITIALSYSTEM_H_
#define INITIALSYSTEM_H_

#include <math.h>

class InitialSystem
{
public:
	InitialSystem(); //Constructor
  ~InitialSystem();

  //Get methods:
  double getPi() const;
  double getGamma() const;
  double getEpsilon() const;
  unsigned int getDissip() const;
  unsigned int getNbiter() const;
  unsigned int getRungekutta() const;

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

  double pi_, gamma_, epsilon_; 

	/* constants from "input file" */
  unsigned int dissip_, nbiter_, rungekutta_;

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

  double rk_beta[5];
  double rk_alpha[5];


};

#endif