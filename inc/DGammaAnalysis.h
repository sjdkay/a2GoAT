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
	TH1*    PTheta;
	TH1*    PPhi;
	//TH1*    PThetaLow;
	//TH1*    PPhiLow;
        //TH1*    PThetaHigh;
	//TH1*    PPhiHigh;
	TH1*    PEp;
	TH1*    PEg;
	TH1*    PVX;
	TH1*    PVY;
	TH1*    PVZ;
	//TH1*    PVXLow;
	//TH1*    PVYLow;
	//TH1*    PVZLow;
	//TH1*    PVXHigh;
	//TH1*    PVYHigh;
	//TH1*    PVZHigh;
	TH1*    nTAPS;
	TH1*    PEpTot;
	TH1*    PEpTot_2;
	TH2* 	EpEg;	
	TH2*    EpTp;
	TH2*    EpdE;
	//TH2*    EpdELow;
	//TH2*    EpdEHigh;	
	TH2*    TpPp;
	TH2*    TpdE;
	TH2*    EpPVZ;
	TH2*    dE1_dE2;
	TH2*    Ep1_Ep2;
	TH2*    PVZ1_PVZ2;
	TH2*    PTheta1_PTheta2;
	TH2*    dE1_dE2_2;
	TH2*    Ep1_Ep2_2;
	TH2*    PVZ1_PVZ2_2;
	TH2*    PTheta1_PTheta2_2;	
	//TH2*    TpdELow;
        //TH2*    TpdEHigh;
	Int_t 	N_P;
	Int_t   j;
		
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
