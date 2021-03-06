#include "InitialSystem.h"
#include <string>
#include <fstream>

#define MAX_MGLEVEL 5

using namespace std;

InitialSystem::InitialSystem() : pi_(3.1416), gamma_(1.4), epsilon_(1.0e-28), dissip_(1), nbiter_(300),
		mach_(0.80), alpha_(1.25), reynolds_(11e6), tinf_(0), xref_(0.25), yref_(0), cmac_(1.0),
		rhoInfini_(1.0), uInfini_(0), vInfini_(0), pInfini_(1.0), 
		cltarget_(0), dcl_(0.001), rms0_(0), imax_(0), jmax_(0), itl_(0), itu_(0){
    for (unsigned int i = 0; i < MAX_MGLEVEL; i++){
		rungekutta_[i] = 5;
	}
	uInfini_ = mach_*sqrt(gamma_)*cos(alpha_);
	vInfini_ = mach_*sqrt(gamma_)*sin(alpha_);
}

InitialSystem::~InitialSystem(){

	if(file_cp_ != nullptr)
	{
    	delete file_cp_;
		file_cp_ = nullptr;
	}

	if(file_conv_ != nullptr)
	{
    	delete file_conv_;
		file_conv_ = nullptr;
	}
}

void InitialSystem::readctrl(string controlFileName){
	unsigned int i;
	double fdum, fdum1;
	string line, str, str1, str2, unusedString;
	double unusedInput;

	//printf("reading %s\n",ctrlfilename_);

	ifstream fpin;
    fpin.open(controlFileName);

	fpin >> title_;
	fpin >> str >> unusedInput >> str1;

	fpin >> imax_ >> str >> jmax_ >> str1 >> itl_ >> str2 >> itu_ >> unusedString; // last one is not in NSCODE, maybe remove

	fpin >> meshfilename_;

	// Or, you know, skip the line itself.
	fpin >> unusedString;
	fpin >> unusedString;
	fpin >> unusedString;
	fpin >> unusedString;
	fpin >> unusedString;


	fpin >> mach_ >> alpha_ >> cltarget_ >> dcl_ >> reynolds_;
	// if cltarget is "no" and not a number, set to "NO_CLRUN"
	// ..whatever that is.

	tinf_ = 300.;

	// Or, you know, skip the line itself.
	fpin >> unusedString;
	fpin >> unusedString;
	fpin >> unusedString;

	fpin >> xref_ >> yref_ >> cmac_;
	
	// Or, you know, skip the line itself.
	fpin >> unusedString;
	fpin >> unusedString;
	fpin >> unusedString;

	fpin >> dissip_ >> fdum >> fdum1; // vis2 vis4

	// Or, you kno, skip the line itself.
	fpin >> unusedString;
	fpin >> unusedString;
	fpin >> unusedString;
	fpin >> unusedString;

	fpin >> fdum; // residual smoothing

	// Or, you know, skip the line itself.
	fpin >> unusedString;

	fpin >> nbiter_;

	// Or, you kno, skip the line itself.
	fpin >> unusedString;
	fpin >> unusedString;
	fpin >> unusedString;
	fpin >> unusedString;
	fpin >> unusedString;
 
	// coarsest_level_=0; // useless

	for (i=0;i<nbiter_;i++){
		fpin >> fdum >> niter_[i] >> fdum1 >> rungekutta_[i] >> itccfl_[i];
	}

	fpin.close();
}

double InitialSystem::getPi() const
{
	return pi_;
}

double InitialSystem::getGamma() const
{
	return gamma_;
}

double InitialSystem::getEpsilon() const
{
	return epsilon_;
}

unsigned int InitialSystem::getDissip() const
{
	return dissip_;
}

unsigned int InitialSystem::getNbiter() const
{
	return nbiter_;
}

double InitialSystem::getMach() const
{
	return mach_;
}

double InitialSystem::getAlpha() const
{
	return alpha_;
}

double InitialSystem::getReynolds() const
{
	return reynolds_;
}

double InitialSystem::getTinf() const
{
	return tinf_;
}

double InitialSystem::getXref() const
{
	return xref_;
}

double InitialSystem::getYref() const
{
	return yref_;
}

double InitialSystem::getCmac() const
{
	return cmac_;
}

double InitialSystem::getRhoInfini() const
{
	return rhoInfini_;
}

double InitialSystem::getUInfini() const
{
	return uInfini_;
}

double InitialSystem::getVInfini() const
{
	return vInfini_;
}

double InitialSystem::getPInfini() const
{
	return pInfini_;
}

double InitialSystem::getClTargert() const
{
	return cltarget_;
}

double InitialSystem::getDcl() const
{
	return dcl_;
}

double InitialSystem::getrms0() const
{
	return rms0_;
}

void InitialSystem::rungeKuttaInit() 
{
	rk_alpha_[0]=0.25; 
    rk_alpha_[1]=0.1666667; 
    rk_alpha_[2]=0.375; 
    rk_alpha_[3]=0.5; 
    rk_alpha_[4]=1.0; 
    rk_beta_[0]=1.0; 
    rk_beta_[1]=0.0; 
    rk_beta_[2]=0.56; 
    rk_beta_[3]=0.0; 
    rk_beta_[4]=0.44;
}