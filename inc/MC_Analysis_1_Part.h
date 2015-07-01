#ifndef __MC_Analysis_1_Part_h__
#define __MC_Analysis_1_Part_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
using namespace std;
#include "GTreeManager.h"
#include "PPhysics.h"
#include "TCutG.h"
#include "TObject.h"
#include "TGraph.h"
#include "TRandom2.h"

class	MC_Analysis_1_Part : public PPhysics
{
private:

  TH1*	TaggerAccScal;

  Int_t NP;
  Int_t NPi;
  Int_t NRoo;
  Int_t NTag;
  Int_t NTrack;
  Int_t PIDElement;

  double_t Time;
  double_t TaggerTime;
  double_t EGamma;
  double_t B;
  double_t Theta1;
  double_t Phi1;
  double_t ThetaB;
  double_t d;
  double_t z1;
  double_t E1;
  double_t dE1;
  double_t dE1Corr;

  TLorentzVector GV1;
  TLorentzVector GV1B;
  TLorentzVector Gamma;
  TLorentzVector Deut;
  TLorentzVector boostvector;
  TVector3 b;
  TVector3 Gamma3;

  TRandom2 rGen;

  GH1* time;
  GH1* time_cut;
  GH1* CM150;
  GH1* CM250;
  GH1* CM350;
  GH1* CM450;
  GH1* CM550;
  GH1* Ek;
  GH1* Eg;
  GH1* Theta;
  GH1* ThetaCM;

  GH1* Z_Vert;

  GH2* E_dE;
  GH2* ThetaPhi;
  GH2* ThetaPhiIn;
  GH2* ThetaPhiOut;
  GH2* PhidEFix;

  char cutfilename[256];
  char cutname[256];
  TFile* CutFile;
  TCutG* Cut;
  TCutG* Cut_proton;
  TCutG* Cut_proton_4Sig;
  TCutG* Cut_proton_Full;
  TCutG* Cut_pion;
  TCutG* Cut_neutron;
  TCutG* Cut_CB_proton;
  TCutG* Cut_CB_proton_4Sig;
  TCutG* Cut_CB_proton_Full;
  TCutG* Cut_CB_pion;
  TCutG* Cut_CB_neutron;

protected:
    virtual Bool_t  Start();
    virtual void    ProcessEvent();
    virtual void    ProcessScalerRead();
    virtual Bool_t  Write();

public:

    MC_Analysis_1_Part();
    virtual ~MC_Analysis_1_Part();
    virtual Bool_t  Init();
    TCutG* OpenCutFile(Char_t* filename, Char_t* cutname);
    Int_t GetEvent();
    TLorentzVector InitialVect();
    Double_t InitialProp();
    Double_t MCSmearing();
    Double_t PNProp(Int_t ProtonParticleNumber);
    TLorentzVector PNVect(Int_t ProtonParticleNumber);
    TVector3 DefineAxes(TVector3 ProtonVector, TVector3 ReconstuctedNeutronVector);
    Double_t LabBoost();
    Double_t LabScatter();
    Double_t nFrameScatter();
    void FillHists();

};
#endif
