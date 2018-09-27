#include "InitialSystem.h"
#include <string>
#include <fstream>

using namespace std;

InitialSystem::InitialSystem() : pi_(3.1416), gamma_(1.4), epsilon_(1.0e-28), dissip_(1), nbiter_(300), rungekutta_(5),
		mach_(0.80), alpha_(1.25), reynolds_(11e6), tinf_(0), xref_(0.25), yref_(0), cmac_(1.0),
		rhoInfini_(0), uInfini_(0), vInfini_(0), pInfini_(0), 
		cltarget_(0), dcl_(0.001), rms0_(0), imax_(0), jmax_(0), itl_(0), itu_(0){}



InitialSystem::~InitialSystem(){

}

void InitialSystem::readctrl(string controlFileName){
	unsigned int nread, imax, jmax, itl, itu, i;
	double fdum, fdum1, fdum2, fdum3, fdum4;
	string line, str, str1, str2, unusedString;
	double unusedInput;

	printf("reading %s\n",ctrlfilename_);

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

	fpin >> niter_;

	// Or, you kno, skip the line itself.
	fpin >> unusedString;
	fpin >> unusedString;
	fpin >> unusedString;
	fpin >> unusedString;
	fpin >> unusedString;
 
	// coarsest_level_=0; // useless

	for (i=0;i<niter_;i++){
		fpin >> fdum >> nbiter_[i] >> fdum1 >> rungekutta_[i] >> itccfl_[i];
	}

	fpin.close();
}

float InitialSystem::getPi() const
{
	return pi_;
}

float InitialSystem::getGamma() const
{
	return gamma_;
}

float InitialSystem::getEpsilon() const
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

unsigned int InitialSystem::getRungekutta() const
{
	return rungekutta_;
}

float InitialSystem::getMach() const
{
	return mach_;
}

float InitialSystem::getAlpha() const
{
	return alpha_;
}

float InitialSystem::getReynolds() const
{
	return reynolds_;
}

float InitialSystem::getTinf() const
{
	return tinf_;
}

float InitialSystem::getXref() const
{
	return xref_;
}

float InitialSystem::getYref() const
{
	return yref_;
}

float InitialSystem::getCmac() const
{
	return cmac_;
}

float InitialSystem::getRhoInfini() const
{
	return rhoInfini_;
}

float InitialSystem::getUInfini() const
{
	return uInfini_;
}

float InitialSystem::getVInfini() const
{
	return vInfini_;
}

float InitialSystem::getPInfini() const
{
	return pInfini_;
}

float InitialSystem::getClTargert() const
{
	return cltarget_;
}

float InitialSystem::getDcl() const
{
	return dcl_;
}

float InitialSystem::getrms0() const
{
	return rms0_;
}

