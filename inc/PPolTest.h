#ifndef __PPolTest_h__
#define __PPolTest_h__

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

class	PPolTest : public PPhysics
{
private:

  TH1*	TaggerAccScal;

  Int_t NP;
  Int_t NPi;
  Int_t NRoo;
  Int_t NTag;
  Int_t NTrack;
  Int_t i;
  Int_t k;
  Int_t MCTrueID1;
  Int_t MCTrueID2;

  double_t Time;
  double_t TaggerTime;
  double_t EGamma;
  double_t Mn;
  double_t Mp;
  double_t d;
  double_t ThetanTrue;
  double_t ThetapTrue;

  TLorentzVector Gamma;
  TLorentzVector Deut;
  TLorentzVector NeutronTrueVect;
  TLorentzVector ProtonTrueVect;

  TVector3 ProtonTrue3Vect;
  TVector3 NeutronTrue3Vect;

  TRandom2 rGen;

  Bool_t MCData;
  Bool_t Proton1;
  Bool_t Proton2;
  Bool_t Neutron1;
  Bool_t Neutron2;

  GH1*	time;
  GH1*	time_cut;

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

    PPolTest();
    virtual ~PPolTest();
    virtual Bool_t  Init();
    TCutG* OpenCutFile(Char_t* filename, Char_t* cutname);
    Bool_t MCDataCheck();
    Int_t GetEvent();
    Int_t MCTrueID();
    Bool_t ParticleAssignment();
    Double_t GetTrueTheta();
    void FillHists();

};
#endif
