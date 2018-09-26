/*This class includes the following functions of NSCODE:
  - initial_system
  - initial_flow_parameter
and the input file.*/

#ifndef INITIALSYSTEM_H_
#define INITIALSYSTEM_H_

class InitialSystem
{
public:
	InitialSystem(); //Constructor
  ~InitialSystem();

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

private:
  float pi_, gamma_, epsilon_; 

	/* constants from "input file" */
  unsigned int dissip_, nbiter_, rungekutta_;

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