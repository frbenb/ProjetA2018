#include "InitialSystem.h"

InitialSystem::InitialSystem() : pi_(3.1416), gamma_(1.4), epsilon_(1.0e-28), dissip_(1), nbiter_(300), rungekutta_(5),
		mach_(0.80), alpha_(1.25), reynolds_(11e6), tinf_(0), xref_(0.25), yref_(0), cmac_(1.0),
		rhoInfini_(1.0), uInfini_(0), vInfini_(0), pInfini_(1.0), 
		cltarget_(0), dcl_(0.001), rms0_(0)
		{
			uInfini_ = mach_*sqrt(gamma_)*cos(alpha_);
			vInfini_ = mach_*sqrt(gamma_)*sin(alpha_);
		}

InitialSystem::~InitialSystem(){

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

