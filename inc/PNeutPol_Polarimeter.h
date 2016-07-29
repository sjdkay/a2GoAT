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

  Int_t PIDEle1;
  Int_t PIDEle2;
  double_t PIDElePhi1;
  double_t PIDElePhi2;
  double_t PIDEleRec1;
  double_t PIDEleRec2;
  double_t PIDElement;
  double_t PIDDiffMeas1;
  double_t PIDDiffMeas2;
  double_t PIDDiffRecMeas1;
  double_t PIDDiffRecMeas2;

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
  double_t Phi1Rec;
  double_t Phi2Rec;
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
  double_t zWCRec1;
  double_t zWC1;
  double_t zWCRec2;
  double_t zWC2;
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
  TLorentzVector GV1Rec;
  TLorentzVector GV2Rec;
  TLorentzVector GVpB;
  TLorentzVector GVnB;
  TLorentzVector GVp;
  TLorentzVector GVn;
  TLorentzVector GVpCalc;
  TLorentzVector GVnCalc;
  TLorentzVector Gamma;
  TLorentzVector Deut;
  TLorentzVector boostvector;
  TVector3 GV1_3;
  TVector3 GV2_3;
  TVector3 GV1Rec_3;
  TVector3 GV2Rec_3;
  TVector3 b;
  TVector3 Gamma3;
  TVector3 GVp3;
  TVector3 GVn3;
  TVector3 GVpCalc3;
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
  GH1* ENeutronVectCalc;
  GH1* ENeutronKinCalc;
  GH1* ENeutronDiff;
  GH1* Theta;
  GH1* ThetaCM;
  GH1* Z_Vert;
  GH1* Zp_Vert;
  GH1* Zn_Vert;
  GH1* Z_WireChamber;
  GH1* Z_WireChamberRec;
  GH1* Z_WireChamberDifference;
  GH1* Z1_WireChamber;
  GH1* Z1_WireChamberRec;
  GH1* Z1_WireChamberDifference;
  GH1* Z2_WireChamber;
  GH1* Z2_WireChamberRec;
  GH1* Z2_WireChamberDifference;
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

  GH2* E_dE;
  GH2* E_dE_Proton;
  GH2* E_dE_Neutron;
  GH2* Thetap_Ep;
  GH2* Thetan_En;

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

  GH1* EgMC_In;
  GH1* ThetapMC_In;
  GH1* ThetanMC_In;
  GH1* ThetanMC_Rec_In;
  GH1* EpMC_In;
  GH1* EnMC_In;
  GH1* EgMC_Out;
  GH1* ThetapMC_Out;
  GH1* ThetanMC_Out;
  GH1* ThetanMC_Rec_Out;
  GH1* EpMC_Out;
  GH1* EnMC_Out;
  GH2* Thetan_dE_MC_In;
  GH2* ThetanRec_dE_MC_In;
  GH2* Phin_dE_MC_In;
  GH2* Thetan_dE_MC_Out;
  GH2* ThetanRec_dE_MC_Out;
  GH2* Phin_dE_MC_Out;
  GH2* MCThetap_Ep;
  GH2* MCThetan_En;
  GH2* MCThetap_Ep_True;
  GH2* MCThetan_En_True;

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
    void ParticleSelection();
    void AlternativeParticleSelection();
    Bool_t MCDataCheck();
    Int_t GetEvent();
    TLorentzVector InitialVect();
    Double_t InitialProp();
    Int_t DetectorElements();
    TLorentzVector ReconstructVectors();
    Int_t ReconstructDetElements();
    Double_t PIDElementsFromPhi(Double_t PhiVal);
    Int_t MCTrueID();
    TLorentzVector MCTrueVectors();
    Double_t MCSmearing();
    Double_t PNProp(Int_t ProtonParticleNumber);
    TLorentzVector PNVect(Int_t ProtonParticleNumber);
    TVector3 DefineAxes(TVector3 ProtonVector, TVector3 ReconstuctedNeutronVector);
    Double_t NeutronEnergy();
    Double_t LabBoost();
    Double_t LabScatter();
    Double_t WCVertex(TVector3 MeasuredVector, TVector3 ReconstructedVector, double_t ReconstructorZ, double_t MeasuredZ);
    Double_t nFrameScatter();
    void FillHists();

};
#endif
