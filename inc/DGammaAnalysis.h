#ifndef __DGammaAnalysis_h__
#define __DGammaAnalysis_h__

#include <iostream>
#include <fstream>
using namespace std;
#include <cstdio>
#include <string> 

#include "PPhysics.h"

class	DGammaAnalysis : public PPhysics
{
private:

	Double_t time;
	Double_t N_PD;
	Double_t DiffD;
	Double_t B;

	TH1* 	time_proton;
	TH1* 	time_proton_cuts;	
	TH1* 	prompt_proton;
	TH1* 	random_proton;
	TH1* 	proton;
	TH1*    Eg_Epsum;
	TH1*    Eg_EpsumPrompt;
	TH1*    Eg_EpsumRandom;
	TH1*	PEpPrompt;
	TH1*    PEpRandom;
	TH1*    PEgPrompt;
	TH1*    PEgRandom;
	TH1*    PEpTotPrompt;
	TH1*    PEpTotRandom;	

	TH1*    PTheta;
	TH1*    PPhi;
	TH1*    PEp;
	TH1*    PEg;
	TH1*    PVX;
	TH1*    PVY;
	TH1*    PVZ;
	TH1*    nTAPS;
	TH1*    PEpTot;
	TH1*    PPhiDiff;
	
	TH2* 	EpEg;	
	TH2*    EpTp;
	TH2*    EpdE;	
	TH2*    TpPp;
	TH2*    TpdE;
	TH2*    EpPVZ;
	TH2*    dE1_dE2;
	TH2*    Ep1_Ep2;
	TH2*    PVZ1_PVZ2;
	TH2*    PTheta1_PTheta2;
	TH2*    Ep1dE1;
	TH2*    Ep2dE2;
	
	
	Int_t 	N_P;
	Int_t   N_P2;
	Int_t   Diff;
	Int_t   k;

	TLorentzVector GV1;
	TLorentzVector GV2;
	TLorentzVector GV1B;
	TLorentzVector GV2B;
	TLorentzVector sum;
	TLorentzVector sumB;
	TLorentzVector boostvector;
	TVector3 b;
	
		
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
