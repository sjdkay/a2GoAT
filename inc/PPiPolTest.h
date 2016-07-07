#ifndef __PPiPolTest_h__
#define __PPiPolTest_h__

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
#include "TMath.h"

class	PPiPolTest : public PPhysics
{
private:

  TH1*	TaggerAccScal;

  Int_t NP;
  Int_t NPi;
  Int_t NPi0;
  Int_t NRoo;
  Int_t NPhotons;
  Int_t NTag;
  Int_t NTrack;
  Int_t MCTrueID1;
  Int_t MCTrueID2;
  Int_t PIDHits1;
  Int_t MWPCHits1;
  Int_t CBHits1;
  Int_t PIDHits2;
  Int_t MWPCHits2;
  Int_t CBHits2;

  double_t i;
  double_t k;
  double_t Time;
  double_t TaggerTime;
  double_t EGamma;
  double_t Mn;
  double_t Mp;
  double_t MpiC;
  double_t Mpi0;
  double_t MN;
  double_t MPi;
  double_t d;
  double_t ThetaNTrue;
  double_t ThetaPiTrue;
  double_t Theta1;
  double_t Phi1;
  double_t z1;
  double_t E1;
  double_t EPiTrue;
  double_t ENTrue;
  double_t EgTrue;
  double_t dE1;

  TLorentzVector Gamma;
  TLorentzVector TrueGamma;
  TLorentzVector Deut;
  TLorentzVector NucleonTrueVect;
  TLorentzVector PionTrueVect;
  TLorentzVector GV1;

  TVector3 PionTrue3Vect;
  TVector3 NucleonTrue3Vect;

  TRandom2 rGen;

  Bool_t MCData;
  Bool_t Nucleon1;
  Bool_t Nucleon2;
  Bool_t Pion1;
  Bool_t Pion2;
  Bool_t NeutralPion;

  GH1* time;
  GH1* time_cut;
  GH1* TrueThetaPi;
  GH1* TrueThetaN;

  GH2* TrueThetaPiThetaN;
  GH2* ThetaPiEPiTrue;
  GH2* ThetaPiEPi;
  GH2* EPiEgTrue;
  GH2* EPiEg;

  GH3* TrueEgThetaPiEPi;
  GH3* TrueEgThetaNEN;
  GH3* EgThetaEPi;
  GH3* EgThetaEPiTrue;

  char cutfilename[256];
  char cutname[256];
  TFile* CutFile;
  TCutG* Cut;
  TCutG* Cut_proton;
  TCutG* Cut_pion;
  TCutG* Cut_ROI;
  TCutG* Cut_neutron;
  TCutG* Cut_CB_proton;
  TCutG* Cut_CB_pion;
  TCutG* Cut_CB_ROI;

protected:
    virtual Bool_t  Start();
    virtual void    ProcessEvent();
    virtual void    ProcessScalerRead();
    virtual Bool_t  Write();

public:

    PPiPolTest();
    virtual ~PPiPolTest();
    virtual Bool_t  Init();
    TCutG* OpenCutFile(Char_t* filename, Char_t* cutname);
    Bool_t MCDataCheck();
    Int_t GetEvent();
    Int_t MCTrueID();
    Bool_t ParticleAssignment();
    Double_t GetTrueTheta();
    Int_t GetDetectorHitInfoChargedPi();
    Int_t GetDetectorHitInfoNeutralPi();
    TLorentzVector GetTrackVector();
    Double_t GetTrackInfo();
    Double_t MCSmearing();
    TLorentzVector GetTrackVectorPi0();
    Double_t GetTrackInfoPi0();
    Double_t MCSmearingPi0();
    void FillHists();

};
#endif
