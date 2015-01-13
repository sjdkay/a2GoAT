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
	TH1* 	time_pi0;
	TH1* 	time_pi0_cuts;	

	TH1* 	MM_prompt_pi0;
	TH1* 	MM_random_pi0;
	TH1* 	MM_pi0;
	
	TH1* 	MM_prompt_pi0_n_2g;
	TH1* 	MM_random_pi0_n_2g;
	TH1* 	MM_pi0_n_2g;		

	Int_t 	N_pi0;
		
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
