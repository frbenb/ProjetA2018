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

	if(file_cp_ != nullptr)
	{
		file_cp_ = nullptr;
	}
	delete file_cp_;

	if(file_conv_ != nullptr)
	{
		file_conv_ = nullptr;
	}
	delete file_conv_;

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

unsigned int InitialSystem::getRungekutta() const
{
	return rungekutta_;
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

