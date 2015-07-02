#ifndef __PNeutPol_h__
#define __PNeutPol_h__

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

class	PNeutPol : public PPhysics
{
private:

  TH1*	TaggerAccScal;

  Int_t NP;
  Int_t NPi;
  Int_t NRoo;
  Int_t NTag;
  Int_t NTrack;

  double_t Time;
  double_t TaggerTime;
  double_t EGamma;
  double_t mm1;
  double_t mm2;
  double_t mmp;
  double_t mmn;
  double_t mmadj;
  double_t GVpUnadjE;
  double_t B;
  double_t Theta1;
  double_t Theta2;
  double_t ThetapB;
  double_t ThetanB;
  double_t PhipB;
  double_t PhinB;
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
  double_t Thetap;
  double_t Thetan;
  double_t ThetanCalc;
  double_t Phip;
  double_t Phin;
  double_t PhinCalc;
  double_t ThetaWCn;
  double_t E1;
  double_t E2;
  double_t Ep;
  double_t En;
  double_t dE1;
  double_t dE2;
  double_t dEp;
  double_t dEn;
  double_t dE1Corr;
  double_t dE2Corr;
  double_t ScattX;
  double_t ScattY;
  double_t ScattZ;
  double_t ScattThetaLab;
  double_t ScattTheta;
  double_t ScattPhi;

  TLorentzVector GV1;
  TLorentzVector GV2;
  TLorentzVector GVpB;
  TLorentzVector GVnB;
  TLorentzVector GVp;
  TLorentzVector GVn;
  TLorentzVector GVnCalc;
  TLorentzVector Gamma;
  TLorentzVector Deut;
  TLorentzVector boostvector;
  TVector3 b;
  TVector3 Gamma3;
  TVector3 GVp3;
  TVector3 GVn3;
  TVector3 GVnCalc3;
  TVector3 fX;
  TVector3 fY;
  TVector3 fZ;
  TVector3 nFrame;

  TRandom2 rGen;

  GH1*	time;
  GH1*	time_cut;

  GH1* MM_Proton;
  GH1* MM_Neutron;
  GH1* Ek;
  GH1* Eg;
  GH1* EgESum;
  GH1* EkSum;
  GH1* Theta;
  GH1* ThetaCM;
  GH1* Z_Vert;
  GH1* Zp_Vert;
  GH1* Zn_Vert;
  GH1* Z_WireChamber;
  GH1* Z_WireChamberRec;
  GH1* Z_WireChamberDifference;
  GH1* ThetaCMProton;
  GH1* ThetaCMNeutron;
  GH1* PhiCMProton;
  GH1* PhiCMNeutron;
  GH1* ThetaWCProton;
  GH1* ThetaWCNeutron;
  GH1* CM150;
  GH1* CM250;
  GH1* CM350;
  GH1* CM450;
  GH1* CM550;
  GH1* ThetaScLab;
  GH1* ThetaSc;
  GH1* PhiSc;
  GH1* PhiScCut;
  GH1* PhiSc125;
  GH1* PhiSc175;
  GH1* PhiSc225;
  GH1* PhiSc275;
  GH1* PhiSc325;
  GH1* PhiSc375;
  GH1* PhiSc425;
  GH1* PhiSc475;
  GH1* PhiSc525;
  GH1* PhiSc575;
  GH1* PhiSc125Cut;
  GH1* PhiSc175Cut;
  GH1* PhiSc225Cut;
  GH1* PhiSc275Cut;
  GH1* PhiSc325Cut;
  GH1* PhiSc375Cut;
  GH1* PhiSc425Cut;
  GH1* PhiSc475Cut;
  GH1* PhiSc525Cut;
  GH1* PhiSc575Cut;

  GH1* PhiSc_In;
  GH1* PhiSc125_In;
  GH1* PhiSc175_In;
  GH1* PhiSc225_In;
  GH1* PhiSc275_In;
  GH1* PhiSc325_In;
  GH1* PhiSc375_In;
  GH1* PhiSc425_In;
  GH1* PhiSc475_In;
  GH1* PhiSc525_In;
  GH1* PhiSc575_In;

  GH1* PhiSc_Out;
  GH1* PhiSc125_Out;
  GH1* PhiSc175_Out;
  GH1* PhiSc225_Out;
  GH1* PhiSc275_Out;
  GH1* PhiSc325_Out;
  GH1* PhiSc375_Out;
  GH1* PhiSc425_Out;
  GH1* PhiSc475_Out;
  GH1* PhiSc525_Out;
  GH1* PhiSc575_Out;

  GH1* PhiScCut_In;
  GH1* PhiSc125Cut_In;
  GH1* PhiSc175Cut_In;
  GH1* PhiSc225Cut_In;
  GH1* PhiSc275Cut_In;
  GH1* PhiSc325Cut_In;
  GH1* PhiSc375Cut_In;
  GH1* PhiSc425Cut_In;
  GH1* PhiSc475Cut_In;
  GH1* PhiSc525Cut_In;
  GH1* PhiSc575Cut_In;

  GH1* PhiScCut_Out;
  GH1* PhiSc125Cut_Out;
  GH1* PhiSc175Cut_Out;
  GH1* PhiSc225Cut_Out;
  GH1* PhiSc275Cut_Out;
  GH1* PhiSc325Cut_Out;
  GH1* PhiSc375Cut_Out;
  GH1* PhiSc425Cut_Out;
  GH1* PhiSc475Cut_Out;
  GH1* PhiSc525Cut_Out;
  GH1* PhiSc575Cut_Out;


  GH2* E_dE;
  GH2* E_dE_Proton;
  GH2* E_dE_Neutron;
  GH2* ThetanCM_PhinCM;
  GH2* Thetan_Vs_Phin_Lab;
  GH2* Theta_Vs_Phi;
  GH2* PhiSc_dEn;
  GH2* mmE;
  GH2* PhidEFixp;
  GH2* PhidECorrFixp;

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

    PNeutPol();
    virtual ~PNeutPol();
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
    Double_t WCVertex();
    Double_t nFrameScatter();
    void FillHists();

};
#endif
