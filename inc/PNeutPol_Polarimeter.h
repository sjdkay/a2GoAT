#ifndef __PNeutPol_Polarimeter_h__
#define __PNeutPol_Polarimeter_h__

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

class	PNeutPol_Polarimeter : public PPhysics
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

  Int_t Detectors1;
  Int_t Detectors2;
  Int_t DetectorsSum;

  double_t Time;
  double_t TaggerTime;
  double_t EGamma;
  double_t Mn;
  double_t Mp;
  double_t mm1;
  double_t mm2;
  double_t mmp;
  double_t mmn;
  double_t mmadj;
  double_t GVpUnadjE;
  double_t B;
  double_t Theta1;
  double_t Theta2;
  double_t Phi1;
  double_t Phi2;
  double_t mm1Diff;
  double_t mm2Diff;
  double_t d;
  double_t lrec;
  double_t l;
  double_t z1;
  double_t z2;
  double_t zdiff;
  double_t ln;
  double_t Zp;
  double_t Zn;
  double_t zWCRec;
  double_t zWC;
  double_t zWCRec1;
  double_t zWC1;
  double_t zWCRec2;
  double_t zWC2;
  double_t Thetap;
  double_t Thetan;
  double_t ThetanRec;
  double_t Phip;
  double_t Phin;
  double_t PhinRec;
  double_t ThetaWCn;
  double_t PhiDiff;
  double_t E1;
  double_t E2;
  double_t Ep;
  double_t En;
  double_t EnVectCalc;
  double_t EnKinCalc;
  double_t dE1;
  double_t dE2;
  double_t dEp;
  double_t dEn;
  double_t ScattX;
  double_t ScattY;
  double_t ScattZ;
  double_t ScattThetaLab;
  double_t ScattTheta;
  double_t ScattPhi;

  Bool_t nBanana;
  Bool_t MCData;
  Bool_t Proton1;
  Bool_t Proton2;

  TLorentzVector GV1;
  TLorentzVector GV2;
  TLorentzVector GVp;
  TLorentzVector GVn;
  TLorentzVector Gamma;
  TLorentzVector Deut;
  TVector3 b;
  TVector3 Gamma3;
  TVector3 GVp3;
  TVector3 GVn3;
  TVector3 GVn3Rec;

  TRandom2 rGen;

  GH1*	time;
  GH1*	time_cut;
  GH1*  Zp_Vert;
  GH1*  Zn_Vert;
  GH1*  Ekp;
  GH1*  Ekn;
  GH1*  EkSum;
  GH1*  Eg;
  GH1*  ThetaProt;
  GH1*  ThetaNeut;
  GH1*  PhiProt;
  GH1*  PhiNeut;

  GH2* E_dE;
  GH2* E_dE_p;
  GH2* E_dE_n;
  GH2* EkEg;
  GH2* EkEg_p;
  GH2* EkEg_n;

  TLorentzVector MCTrueVect1;
  TLorentzVector MCTrueVect2;

  Int_t MCTrueID1;
  Int_t MCTrueID2;

  double_t MCTheta1;
  double_t MCTheta2;
  double_t MCE1;
  double_t MCE2;
  double_t MCTheta1True;
  double_t MCTheta2True;
  double_t MCE1True;
  double_t MCE2True;

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

    PNeutPol_Polarimeter();
    virtual ~PNeutPol_Polarimeter();
    virtual Bool_t  Init();
    TCutG* OpenCutFile(Char_t* filename, Char_t* cutname);
    void MCHists();
    Bool_t MCDataCheck();
    Int_t GetEvent();
    TLorentzVector InitialVect();
    Double_t InitialProp();
    Int_t DetectorCheck();
    Int_t MCTrueID();
    TLorentzVector MCTrueVectors();
    Double_t MCSmearing();
    Double_t PNProp(Int_t ProtonParticleNumber);
    TLorentzVector PNVect(Int_t ProtonParticleNumber);
    Double_t LabAngles();
    void FillHists();

};
#endif
