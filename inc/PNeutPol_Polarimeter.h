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
  Int_t PIDEle;
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
  double_t GVpUnadjE;
  double_t Theta1;
  double_t Theta2;
  double_t Phi1;
  double_t Phi2;
  double_t WC1X1;
  double_t WC1Y1;
  double_t WC1Z1;
  double_t WC1X2;
  double_t WC1Y2;
  double_t WC1Z2;
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
  double_t WC1pX;
  double_t WC1pY;
  double_t WC1pZ;
  double_t WC1nX;
  double_t WC1nY;
  double_t WC1nZ;
  double_t WCThetap;
  double_t WCThetan;
  double_t WCPhip;
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
  double_t KinEp;
  double_t KinEpWC;
  double_t KinEpDiff;
  double_t KinEDiff;
  double_t KinEpMB;
  double_t KinEpMB2;
  double_t EpTot;
  double_t Pp;
  double_t Ppx;
  double_t Ppy;
  double_t Ppz;
  double_t MMpKin;
  double_t MMpKinMB;
  double_t MMpKinMB2;

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
  TLorentzVector P4Vect;
  TLorentzVector N4Vect;
  TLorentzVector RecKinProton;
  TLorentzVector RecKinNeutron;
  TLorentzVector RecKinMBProton;
  TLorentzVector RecKinMBNeutron;
  TLorentzVector RecKinMBProton2;
  TLorentzVector RecKinMBNeutron2;
  TVector3 GVp3;
  TVector3 GVn3;
  TVector3 GVn3Rec;
  TVector3 WC3Vectp;
  TVector3 WC3Vectn;

  TRandom2 rGen;

  TH1D*	time;
  TH1D*	time_cut;

  GH1*  Zp_Vert;
  GH1*  Zn_Vert;
  GH1*  Ekp;
  GH1*  Ekn;
  GH1*  EkSum;
  GH1*  Eg;
  GH1*  EgCut;
  GH1*  ThetaProt;
  GH1*  ThetaNeut;
  GH1*  PhiProt;
  GH1*  PhiNeut;
  GH1*  WCPhiDifference;
  GH1*  WCThetaProt;
  GH1*  WCThetaNeut;
  GH1*  WCPhiProt;
  GH1*  WCPhiNeut;
  GH1*  EpKin;
  GH1*  EpKinMB;

  GH1*  WCXp;
  GH1*  WCYp;
  GH1*  WCZp;
  GH1*  WCXn;
  GH1*  WCYn;
  GH1*  WCZn;
  GH1*  MMp;
  GH1*  MMpMB;
  GH1*  MMpMB2;
  GH1*  MMpCut;
  GH1*  MMpMBCut;
  GH1*  MMpMB2Cut;

  GH2* E_dE;
  GH2* E_dE_Cut;
  GH2* EpKinEpKinMBDiffPTheta;
  GH2* MMpKinEKin;
  GH2* MMpKinEKinMB;
  GH2* MMpKinEKinCut;
  GH2* MMpKinEKinMBCut;
  GH2* MMpKinTheta;
  GH2* MMpKinThetaMB;

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
    Int_t GetEvent();
    TLorentzVector InitialVect();
    Double_t InitialProp();
    Int_t DetectorCheck();
    Double_t PNProp(Int_t ProtonParticleNumber);
    TLorentzVector PNVect(Int_t ProtonParticleNumber);
    TVector3 WC3Vectors(Double_t WCpX, Double_t WCpY, Double_t WCpZ, Double_t WCnX, Double_t WCnY, Double_t WCnZ);
    Double_t WCAngles(TVector3 MWPCp3Vector, TVector3 MWPCn3Vector);
    TLorentzVector Proton4VectorKin(Double_t KinE, Double_t Theta, Double_t Phi);
    TLorentzVector Neutron4VectorKin(TLorentzVector ProtonKinVector);
    Double_t LabAngles();
    void FillHists();

};
#endif
