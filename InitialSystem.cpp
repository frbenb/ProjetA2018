#include "InitialSystem.h"
#include "Mesh.h"
#include <string>
#include <fstream>

using namespace std;

InitialSystem::InitialSystem() : pi_(3.1416), gamma_(1.4), epsilon_(1.0e-28), dissip_(1), nbiter_(300), rungekutta_(5),
		mach_(0.80), alpha_(1.25), reynolds_(11e6), tinf_(0), xref_(0.25), yref_(0), cmac_(1.0),
		rhoInfini_(0), uInfini_(0), vInfini_(0), pInfini_(0), 
		cltarget_(0), dcl_(0.001), rms0_(0), mesh_(nullptr){}

InitialSystem::~InitialSystem(){

}

void InitialSystem::readctrl(string controlFileName){
	unsigned int nread, imax, jmax, itl, itu, level, i;
	float fdum, fdum1, fdum2, fdum3, fdum4;
	string line, str, str1, str2, unusedString;
	double unusedInput;

	printf("reading %s\n",ctrlfilename_);

	ifstream fpin;
    fpin.open(controlFileName);

	fpin >> title_;
	fpin >> str >> unusedInput >> str1;

	level=0; // ideally remove

	fpin >> imax >> str >> jmax >> str1 >> itl >> str2 >> itu >> unusedString; // last one is not in NSCODE, maybe remove
	
	mesh_ = new Mesh(imax,jmax,itl,itu, this);

	fpin >> meshfilename_;

	// Or, you know, skip the line itself.
	fpin >> unusedString;
	fpin >> unusedString;
	fpin >> unusedString;
	fpin >> unusedString;
	fpin >> unusedString;

	fpin >>



	fgets(line,STLEN,fpin);
	nread=sscanf(line,"%f %f %f %f %f",&fdum,&fdum1,&fdum2,&fdum3,&fdum4);
	printf("nread is %d\n",nread);
	if (nread==5)
	{
		nsc.mach=fdum;
		nsc.alpha=fdum1;
		nsc.cltarget=fdum2;
		nsc.dcl=fdum3;
		nsc.reynolds=fdum4;
	}
	else
	{
		nread=sscanf(line,"%f %f %s %f %f",&fdum,&fdum1,str,&fdum3,&fdum4);
		printf("nread is %d\n",nread);
		if (nread!=5) printerror("mach alpha cl dcl reynolds not well specified");
		nsc.mach=fdum;
		nsc.alpha=fdum1;
		nsc.cltarget=NO_CLRUN;
		nsc.dcl=fdum3;
		nsc.reynolds=fdum4;
	}
	nsc.tinf=300.;
	printf("%f %f %f %f %f %f\n",nsc.mach,nsc.alpha,nsc.cltarget,nsc.dcl,nsc.reynolds,nsc.tinf);  
	fgets(line,STLEN,fpin);
	fgets(line,STLEN,fpin);
	sscanf(line,"%f %f %f",&fdum,&fdum1,&fdum2);
	nsc.xref=fdum;	  
	nsc.yref=fdum1;	  
	nsc.cmac=fdum2;  
	printf("%f %f %f\n",nsc.xref,nsc.yref,nsc.cmac);
	
	
	/* solver */
	fgets(line,STLEN,fpin);
	fgets(line,STLEN,fpin);
	sscanf(line,"%d %f %f",&nsc.dissip,&fdum,&fdum1);
	nsc.vis2=fdum;
	nsc.vis4=fdum1;
	printf("dissip %d vis2 %f vis4 %f\n",nsc.dissip,nsc.vis2,nsc.vis4);
	nsc.vis2=nsc.vis2/2.;  
	nsc.vis4=nsc.vis4/32.;  
	fgets(line,STLEN,fpin);
	fgets(line,STLEN,fpin);
	sscanf(line,"%f",&fdum);
	nsc.ressmoo=fdum;
	printf("ressmoo %f\n",nsc.ressmoo);  
	fgets(line,STLEN,fpin);
	fgets(line,STLEN,fpin);
	sscanf(line,"%d",&nsc.nitc);
	fgets(line,STLEN,fpin);
	printf("   itc level  iter  mglevel  rk    cfl\n");  
	nsc.coarsest_level=0;
	for (i=0;i<nsc.nitc;i++)
	{
		fgets(line,STLEN,fpin);
		sscanf(line,"%d %d %d  %d %f",&nsc.itclevel[i],&nsc.niter[i],&nsc.mglevel[i],&nsc.rklevel[i],&fdum);
		nsc.itccfl[i]=fdum;
		if (nsc.itclevel[i]> nsc.coarsest_level) nsc.coarsest_level=nsc.itclevel[i];
		printf("%5d %5d %5d %5d  %5d %10.5f\n",i,nsc.itclevel[i],nsc.niter[i],nsc.mglevel[i],nsc.rklevel[i],nsc.itccfl[i]);
		printf("%5d\n",nsc.coarsest_level);
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

