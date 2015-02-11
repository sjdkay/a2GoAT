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
	Double_t Theta1;
	Double_t Theta2;
	Double_t Theta1B;
	Double_t Theta2B;
	Double_t P1CalcTheta;
	Double_t P1ThetaDiff;
	Double_t P2CalcTheta;
	Double_t P2ThetaDiff;
	Double_t ThetaNReal1;
	Double_t ThetaNReal2;
	Double_t l1;
	Double_t l2;

	TH1* 	time_proton;
	TH1* 	time_proton_cuts;	
	TH1* 	prompt_proton;
	TH1* 	random_proton;
	TH1*    Eg_EpsumPrompt;
	TH1*    Eg_EpsumRandom;
	TH1*    PEgPrompt;
	TH1*    PEgRandom;
	TH1*    PThetaCMPrompt;
	TH1*    PThetaCMRandom;

	TH1* 	proton;
	TH1*    Eg_Epsum;
	TH1*    PTheta;
	TH1*    PThetaCM;
	TH1*    PEp;
	TH1*    PEg;
	TH1*    PEpTot;
	TH1*    P1CDiff;
	TH1*    P2CDiff;
	TH1*    MM;
	TH1*    MM2;
	
	TH2*    EpdE;

	Int_t 	N_P;
	Int_t   N_P2;
	Int_t   N_Part;
	Int_t   Diff;
	Int_t   k;
	Int_t   a;

	TLorentzVector GV1;
	TLorentzVector GV2;
	TLorentzVector Gamma;
	TLorentzVector Deut;
	TLorentzVector boostvector;
	TVector3 b;

	char 	        cutfilename[256];
	char		cutname[256];	
	TFile* 		CutFile;
	TCutG* 		Cut;
	TCutG* 		Cut_proton;
	TCutG* 		Cut_proton_4Sig;
	TCutG* 		Cut_proton_Full;	
	TCutG*		Cut_pion;
	TCutG*		Cut_neutron;
	TCutG* 		Cut_CB_proton;
	TCutG* 		Cut_CB_proton_4Sig;
	TCutG* 		Cut_CB_proton_Full;
	TCutG* 		Cut_CB_pion;
	TCutG* 		Cut_CB_neutron;
		
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
    TCutG*	OpenCutFile(Char_t* filename, Char_t* cutname);
};
#endif
