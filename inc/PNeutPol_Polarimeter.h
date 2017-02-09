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
  Int_t EventCounter;
  Int_t EventCounterTrackCut;
  Int_t EventCounterZCut;
  Int_t EventCounterCoplanarCut;

  Int_t Detectors1;
  Int_t Detectors2;
  Int_t DetectorsSum;

  double_t Time;
  double_t TaggerTime;
  double_t EGamma;
  double_t Mn;
  double_t Mp;
  double_t Md;
  double_t z1;
  double_t z2;
  double_t ln;
  double_t Zp;
  double_t Zn;
  double_t WC1pX;
  double_t WC1pY;
  double_t WC1pZ;
  double_t WC1nX;
  double_t WC1nY;
  double_t WC1nZ;
  double_t WCThetap;
  double_t WCThetapRad;
  double_t WCThetan;
  double_t WCPhip;
  double_t WCPhipRad;
  double_t WCPhin;
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
  double_t PhiWCDiff;
  double_t Ep;
  double_t En;
  double_t EnVectCalc;
  double_t EnKinCalc;
  double_t dEp;
  double_t dEn;
  double_t ScattX;
  double_t ScattY;
  double_t ScattZ;
  double_t ScattThetaLab;
  double_t ScattTheta;
  double_t ScattPhi;
  double_t KinEp;
  double_t KinEpWC;
  double_t KinEpDiff;
  double_t KinEDiff;
  double_t EpDiff;
  double_t EpCorr;
  double_t KinEpMB;
  double_t EpTot;
  double_t Pp;
  double_t Ppx;
  double_t Ppy;
  double_t Ppz;
  double_t MMpKin;
  double_t MMpEpCorr;
  double_t MMpKinMB;
  double_t OpeningAngle;
  double_t ThetanDiff;
  double_t PhinDiff;

  Bool_t nBanana;
  Bool_t Proton1;
  Bool_t Proton2;
  Bool_t BeamHelicity;

  TLorentzVector GVp;
  TLorentzVector GVn;
  TLorentzVector Gamma;
  TLorentzVector Deut;
  TLorentzVector P4Vect;
  TLorentzVector N4Vect;
  TLorentzVector RecKinProton;
  TLorentzVector RecKinNeutron;
  TLorentzVector RecProtonEpCorr;
  TLorentzVector RecNeutronEpCorr;
  TLorentzVector RecKinMBProton;
  TLorentzVector RecKinMBNeutron;
  TVector3 GVp3;
  TVector3 GVn3;
  TVector3 GVn3Rec;
  TVector3 WC3Vectp;
  TVector3 WC3Vectn;
  TVector3 P3Vect;
  TVector3 N3Vect;
  TVector3 RecProtonEpCorr3;
  TVector3 RecNeutronEpCorr3;

  TH1D*	time;
  TH1D*	time_cut;

  GH1* Zp_Vert;
  GH1* Zn_Vert;
  GH1* Ekp;
  GH1* Ekn;
  GH1* EkSum;
  GH1* Eg;
  GH1* ThetaProt;
  GH1* ThetaNeut;
  GH1* PhiProt;
  GH1* PhiNeut;
  GH1* WCPhiDifference;
  GH1* WCThetaProt;
  GH1* WCThetaNeut;
  GH1* WCPhiProt;
  GH1* WCPhiNeut;
  GH1* EpKin;
  GH1* EpCorrected;
  GH1* EpKinEpCorrDiff;
  GH1* EpEpCorrDiff;

  GH1* WCXp;
  GH1* WCYp;
  GH1* WCZp;
  GH1* WCXn;
  GH1* WCYn;
  GH1* WCZn;
  GH1* MMp;
  GH1* MMpEpCorrected;
  GH1* OAngle;

  GH1* MMp200300;
  GH1* MMp300400;
  GH1* MMp400500;
  GH1* MMp500600;
  GH1* MMp600700;
  GH1* MMp700800;
  GH1* MMp800900;
  GH1* EgCut;
  GH1* MMpEpCorrectedCut;
  GH1* OAngleCut;
  GH1* OAngleCut200400;
  GH1* ScattFrameTheta;
  GH1* ScattFramePhi;

  GH1* ThetaSc;
  GH1* PhiSc;
  GH1* PhiSc275;
  GH1* PhiSc325;
  GH1* PhiSc375;
  GH1* PhiSc425;
  GH1* PhiSc475;
  GH1* PhiSc525;
  GH1* PhiSc575;
  GH1* PhiSc625;
  GH1* PhiSc675;
  GH1* PhiSc725;
  GH1* PhiSc775;
  GH1* PhiSc825;
  GH1* PhiSc875;

  GH1* PhiScNegHel;
  GH1* PhiScPosHel;

  GH1* PhiSc275NegHel;
  GH1* PhiSc325NegHel;
  GH1* PhiSc375NegHel;
  GH1* PhiSc425NegHel;
  GH1* PhiSc475NegHel;
  GH1* PhiSc525NegHel;
  GH1* PhiSc575NegHel;
  GH1* PhiSc625NegHel;
  GH1* PhiSc675NegHel;
  GH1* PhiSc725NegHel;
  GH1* PhiSc775NegHel;
  GH1* PhiSc825NegHel;
  GH1* PhiSc875NegHel;

  GH1* PhiSc275PosHel;
  GH1* PhiSc325PosHel;
  GH1* PhiSc375PosHel;
  GH1* PhiSc425PosHel;
  GH1* PhiSc475PosHel;
  GH1* PhiSc525PosHel;
  GH1* PhiSc575PosHel;
  GH1* PhiSc625PosHel;
  GH1* PhiSc675PosHel;
  GH1* PhiSc725PosHel;
  GH1* PhiSc775PosHel;
  GH1* PhiSc825PosHel;
  GH1* PhiSc875PosHel;

  GH1* PhiSc275NegHelTheta010;
  GH1* PhiSc325NegHelTheta010;
  GH1* PhiSc375NegHelTheta010;
  GH1* PhiSc425NegHelTheta010;
  GH1* PhiSc475NegHelTheta010;
  GH1* PhiSc525NegHelTheta010;
  GH1* PhiSc575NegHelTheta010;
  GH1* PhiSc625NegHelTheta010;
  GH1* PhiSc675NegHelTheta010;
  GH1* PhiSc725NegHelTheta010;
  GH1* PhiSc775NegHelTheta010;
  GH1* PhiSc825NegHelTheta010;
  GH1* PhiSc875NegHelTheta010;

  GH1* PhiSc275PosHelTheta010;
  GH1* PhiSc325PosHelTheta010;
  GH1* PhiSc375PosHelTheta010;
  GH1* PhiSc425PosHelTheta010;
  GH1* PhiSc475PosHelTheta010;
  GH1* PhiSc525PosHelTheta010;
  GH1* PhiSc575PosHelTheta010;
  GH1* PhiSc625PosHelTheta010;
  GH1* PhiSc675PosHelTheta010;
  GH1* PhiSc725PosHelTheta010;
  GH1* PhiSc775PosHelTheta010;
  GH1* PhiSc825PosHelTheta010;
  GH1* PhiSc875PosHelTheta010;

  GH1* PhiSc275NegHelTheta1020;
  GH1* PhiSc325NegHelTheta1020;
  GH1* PhiSc375NegHelTheta1020;
  GH1* PhiSc425NegHelTheta1020;
  GH1* PhiSc475NegHelTheta1020;
  GH1* PhiSc525NegHelTheta1020;
  GH1* PhiSc575NegHelTheta1020;
  GH1* PhiSc625NegHelTheta1020;
  GH1* PhiSc675NegHelTheta1020;
  GH1* PhiSc725NegHelTheta1020;
  GH1* PhiSc775NegHelTheta1020;
  GH1* PhiSc825NegHelTheta1020;
  GH1* PhiSc875NegHelTheta1020;

  GH1* PhiSc275PosHelTheta1020;
  GH1* PhiSc325PosHelTheta1020;
  GH1* PhiSc375PosHelTheta1020;
  GH1* PhiSc425PosHelTheta1020;
  GH1* PhiSc475PosHelTheta1020;
  GH1* PhiSc525PosHelTheta1020;
  GH1* PhiSc575PosHelTheta1020;
  GH1* PhiSc625PosHelTheta1020;
  GH1* PhiSc675PosHelTheta1020;
  GH1* PhiSc725PosHelTheta1020;
  GH1* PhiSc775PosHelTheta1020;
  GH1* PhiSc825PosHelTheta1020;
  GH1* PhiSc875PosHelTheta1020;

  GH1* PhiSc275NegHelTheta2030;
  GH1* PhiSc325NegHelTheta2030;
  GH1* PhiSc375NegHelTheta2030;
  GH1* PhiSc425NegHelTheta2030;
  GH1* PhiSc475NegHelTheta2030;
  GH1* PhiSc525NegHelTheta2030;
  GH1* PhiSc575NegHelTheta2030;
  GH1* PhiSc625NegHelTheta2030;
  GH1* PhiSc675NegHelTheta2030;
  GH1* PhiSc725NegHelTheta2030;
  GH1* PhiSc775NegHelTheta2030;
  GH1* PhiSc825NegHelTheta2030;
  GH1* PhiSc875NegHelTheta2030;

  GH1* PhiSc275PosHelTheta2030;
  GH1* PhiSc325PosHelTheta2030;
  GH1* PhiSc375PosHelTheta2030;
  GH1* PhiSc425PosHelTheta2030;
  GH1* PhiSc475PosHelTheta2030;
  GH1* PhiSc525PosHelTheta2030;
  GH1* PhiSc575PosHelTheta2030;
  GH1* PhiSc625PosHelTheta2030;
  GH1* PhiSc675PosHelTheta2030;
  GH1* PhiSc725PosHelTheta2030;
  GH1* PhiSc775PosHelTheta2030;
  GH1* PhiSc825PosHelTheta2030;
  GH1* PhiSc875PosHelTheta2030;

  GH2* E_dE;
  GH2* E_dE_Cut;
  GH2* KinEp_dE;
  GH2* KinEp_dE_GoodCut;
  GH2* ThetaScPhiSc;
  GH2* OAnglePhiDiffLab;
  GH2* ThetaDiffPhiDiffLab;

  char cutfilename[256];
  char cutname[256];
  TFile* CutFile;
  TCutG* Cut;
  TCutG* Cut_proton;
  TCutG* Cut_pion;
  TCutG* Cut_protonKinGood;
  TCutG* Cut_protonKinBad;
  TCutG* Cut_CB_proton;
  TCutG* Cut_CB_pion;
  TCutG* Cut_CB_protonKinGood;
  TCutG* Cut_CB_protonKinBad;

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
    TLorentzVector Proton4VectorKin(Double_t KinE, Double_t Theta, Double_t Phi);
    TLorentzVector Neutron4VectorKin(TLorentzVector ProtonKinVector);
    Double_t LabAngles();
    void FillHists();

};
#endif
