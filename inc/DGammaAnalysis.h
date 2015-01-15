#ifndef __DGammaAnalysis_h__
#define __DGammaAnalysis_h__

#include <iostream>
#include <fstream>
using namespace std;
#include <cstdio>
#include <string> 

#include "PPhysics.h"

//This may not actually have all of the required things and may need to change defintion of variables
//E.g. change names to things we need

class	DGammaAnalysis : public PPhysics
{
private:

	Double_t time;
	TH1*    PTheta;
	TH1*    PEk;
	TH1*    PEg;
	TH2* 	EpEg;	//2D histogram of Proton energy against photon energy, could do this as 3D? Will this even work?
	//Need to define something that will fill this histogram?
	TH2*    EpTp; //2D Histogram of Proton energy vs angle theta
	Int_t 	N_P;
		
protected:

public:
    DGammaAnalysis();
    virtual ~DGammaAnalysis();

    virtual Bool_t	Init(const char* configfile);	
    virtual Bool_t	File(const char* file_in, const char* file_out);    
    virtual void 	Analyse();
    virtual void	Reconstruct();
    void	PostReconstruction();

    void	DefineHistograms();
    Bool_t	WriteHistograms(TFile* pfile);
    Bool_t	WriteHistograms() {return WriteHistograms(HistFile);}
};
#endif
