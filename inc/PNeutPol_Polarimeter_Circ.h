#ifndef __PNeutPol_Polarimeter_Circ_h__
#define __PNeutPol_Polarimeter_Circ_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <APLCON.hpp>
using namespace std;
#include "GTreeManager.h"
#include "PPhysics.h"
#include "TCutG.h"
#include "TObject.h"
#include "TGraph.h"
#include "TRandom2.h"
#include "TMath.h"

class	PNeutPol_Polarimeter_Circ : public PPhysics
{
private:

  TH1*	TaggerAccScal;

  Int_t NP;
  Int_t NPi;
  Int_t NRoo;
  Int_t NTag;
  Int_t NTrack;
  Int_t EventNum;

  Int_t Detectors1;
  Int_t Detectors2;
  Int_t DetectorsSum;

  double_t Time;
  double_t TaggerTime;
  double_t Timep;
  double_t Timen;
  double_t EGamma;
  double_t Mn;
  double_t Mp;
  double_t Md;
  double_t Mpi;
  double_t z1;
  double_t z2;
  double_t ln;
  double_t Xp;
  double_t Yp;
  double_t Zp;
  double_t Xn;
  double_t Yn;
  double_t Zn;
  double_t Thp;
  double_t ThpRad;
  double_t Thn;
  double_t Php;
  double_t PhpRad;
  double_t Phn;
  double_t WC1pX;
  double_t WC1pY;
  double_t WC1pZ;
  double_t WC1nX;
  double_t WC1nY;
  double_t WC1nZ;
  double_t WC2nX;
  double_t WC2nY;
  double_t WC2nZ;
  double_t WCZnRec;
  double_t ThetapCM;
  double_t CosThetapCM;
  double_t Thetan;
  double_t ThetapRec;
  double_t ThetanRec;
  double_t ThetaPiRec;
  double_t Pp;
  double_t Pn;
  double_t PpKin;
  double_t PhipRec;
  double_t PhinRec;
  double_t PhiPiRec;
  double_t ThetaWCn;
  double_t PhiDiff;
  double_t ThetaPiRecDiff;
  double_t ThetapRecDiff;
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
  double_t KinEpPi;
  double_t KinPi;
  double_t KinEDiff;
  double_t EpDiff;
  double_t EpCorr;
  double_t EpTot;
  double_t Ppx;
  double_t Ppy;
  double_t Ppz;
  double_t MMpEpCorr;
  double_t OpeningAngle;
  double_t ThetanDiff;
  double_t PhinDiff;
  double_t Ncor1;
  double_t Ncor2;
  double_t Ncor3;
  double_t NcorR;
  double_t NcorRR;
  double_t ThetanCorr;
  double_t DOCA;
  double_t POCAx;
  double_t POCAy;
  double_t POCAz;
  double_t r;

  Bool_t nBanana;
  Bool_t Proton1;
  Bool_t Proton2;
  Bool_t BeamHelicity;

  TLorentzVector B;
  TLorentzVector GVp;
  TLorentzVector GVn;
  TLorentzVector GVpCorr;
  TLorentzVector GVpCorrB;
  TLorentzVector GVnCorr;
  TLorentzVector Gamma;
  TLorentzVector Deut;
  TLorentzVector Neut;
  TLorentzVector P4Vect;
  TLorentzVector N4Vect;
  TLorentzVector N4VectCorr;
  TLorentzVector Pi4Vect;
  TLorentzVector RecKinProton;
  TLorentzVector RecKinNeutron;
  TLorentzVector RecKinProtonPi;
  TLorentzVector RecKinPion;
  TLorentzVector RecKinPionP;
  TLorentzVector RecKinPPi;
  TLorentzVector RecProtonEpCorr;
  TLorentzVector RecNeutronEpCorr;
  TVector3 b;
  TVector3 pVertex;
  TVector3 nVertex;
  TVector3 DOCAVertex1;
  TVector3 DOCAVertex2;
  TVector3 POCA;
  TVector3 GVpCorr3;
  TVector3 GVnCorr3;
  TVector3 GVn3Unit;
  TVector3 GVnCorr3Unit;
  TVector3 GVn3Rec;
  TVector3 WC3Vectp;
  TVector3 WC13Vectn;
  TVector3 WC23Vectn;
  TVector3 P3Vect;
  TVector3 N3Vect;
  TVector3 N3VectUnit;
  TVector3 RecProtonEpCorr3;
  TVector3 RecNeutronEpCorr3;

  TH1D*	time;
  TH1D*	time_cut;

  GH1* EkSum;
  GH1* Eg;
  GH1* ThetaProt655705;
  GH1* PhiDifference;
  GH1* EpKin;
  GH1* EpCorrected;
  GH1* EpKinEpCorrDiff;
  GH1* EpEpCorrDiff;

  GH1* MMpEpCorrected;
  GH1* OAngle;
  GH1* WCZnRecon;

  GH1* MMp200300;
  GH1* MMp300400;
  GH1* MMp400500;
  GH1* MMp500600;
  GH1* MMp600700;
  GH1* MMp700800;
  GH1* MMp800900;

  GH1* ZpDist;
  GH1* ZpPhiScatNeg180;
  GH1* ZpPhiScat0;
  GH1* ZpPhiScatPos180;

  GH1* ThetaSc;
  GH1* PhiSc;

  GH1* PhiSc315NegHelCM1;
  GH1* PhiSc345NegHelCM1;
  GH1* PhiSc375NegHelCM1;
  GH1* PhiSc405NegHelCM1;
  GH1* PhiSc435NegHelCM1;
  GH1* PhiSc465NegHelCM1;
  GH1* PhiSc495NegHelCM1;
  GH1* PhiSc525NegHelCM1;
  GH1* PhiSc555NegHelCM1;
  GH1* PhiSc585NegHelCM1;
  GH1* PhiSc615NegHelCM1;
  GH1* PhiSc645NegHelCM1;
  GH1* PhiSc675NegHelCM1;
  GH1* PhiSc315NegHelCM2;
  GH1* PhiSc345NegHelCM2;
  GH1* PhiSc375NegHelCM2;
  GH1* PhiSc405NegHelCM2;
  GH1* PhiSc435NegHelCM2;
  GH1* PhiSc465NegHelCM2;
  GH1* PhiSc495NegHelCM2;
  GH1* PhiSc525NegHelCM2;
  GH1* PhiSc555NegHelCM2;
  GH1* PhiSc585NegHelCM2;
  GH1* PhiSc615NegHelCM2;
  GH1* PhiSc645NegHelCM2;
  GH1* PhiSc675NegHelCM2;
  GH1* PhiSc315NegHelCM3;
  GH1* PhiSc345NegHelCM3;
  GH1* PhiSc375NegHelCM3;
  GH1* PhiSc405NegHelCM3;
  GH1* PhiSc435NegHelCM3;
  GH1* PhiSc465NegHelCM3;
  GH1* PhiSc495NegHelCM3;
  GH1* PhiSc525NegHelCM3;
  GH1* PhiSc555NegHelCM3;
  GH1* PhiSc585NegHelCM3;
  GH1* PhiSc615NegHelCM3;
  GH1* PhiSc645NegHelCM3;
  GH1* PhiSc675NegHelCM3;
  GH1* PhiSc315NegHelCM4;
  GH1* PhiSc345NegHelCM4;
  GH1* PhiSc375NegHelCM4;
  GH1* PhiSc405NegHelCM4;
  GH1* PhiSc435NegHelCM4;
  GH1* PhiSc465NegHelCM4;
  GH1* PhiSc495NegHelCM4;
  GH1* PhiSc525NegHelCM4;
  GH1* PhiSc555NegHelCM4;
  GH1* PhiSc585NegHelCM4;
  GH1* PhiSc615NegHelCM4;
  GH1* PhiSc645NegHelCM4;
  GH1* PhiSc675NegHelCM4;
  GH1* PhiSc315NegHelCM5;
  GH1* PhiSc345NegHelCM5;
  GH1* PhiSc375NegHelCM5;
  GH1* PhiSc405NegHelCM5;
  GH1* PhiSc435NegHelCM5;
  GH1* PhiSc465NegHelCM5;
  GH1* PhiSc495NegHelCM5;
  GH1* PhiSc525NegHelCM5;
  GH1* PhiSc555NegHelCM5;
  GH1* PhiSc585NegHelCM5;
  GH1* PhiSc615NegHelCM5;
  GH1* PhiSc645NegHelCM5;
  GH1* PhiSc675NegHelCM5;
  GH1* PhiSc315NegHelCM6;
  GH1* PhiSc345NegHelCM6;
  GH1* PhiSc375NegHelCM6;
  GH1* PhiSc405NegHelCM6;
  GH1* PhiSc435NegHelCM6;
  GH1* PhiSc465NegHelCM6;
  GH1* PhiSc495NegHelCM6;
  GH1* PhiSc525NegHelCM6;
  GH1* PhiSc555NegHelCM6;
  GH1* PhiSc585NegHelCM6;
  GH1* PhiSc615NegHelCM6;
  GH1* PhiSc645NegHelCM6;
  GH1* PhiSc675NegHelCM6;
  GH1* PhiSc315NegHelCM7;
  GH1* PhiSc345NegHelCM7;
  GH1* PhiSc375NegHelCM7;
  GH1* PhiSc405NegHelCM7;
  GH1* PhiSc435NegHelCM7;
  GH1* PhiSc465NegHelCM7;
  GH1* PhiSc495NegHelCM7;
  GH1* PhiSc525NegHelCM7;
  GH1* PhiSc555NegHelCM7;
  GH1* PhiSc585NegHelCM7;
  GH1* PhiSc615NegHelCM7;
  GH1* PhiSc645NegHelCM7;
  GH1* PhiSc675NegHelCM7;
  GH1* PhiSc315NegHelCM8;
  GH1* PhiSc345NegHelCM8;
  GH1* PhiSc375NegHelCM8;
  GH1* PhiSc405NegHelCM8;
  GH1* PhiSc435NegHelCM8;
  GH1* PhiSc465NegHelCM8;
  GH1* PhiSc495NegHelCM8;
  GH1* PhiSc525NegHelCM8;
  GH1* PhiSc555NegHelCM8;
  GH1* PhiSc585NegHelCM8;
  GH1* PhiSc615NegHelCM8;
  GH1* PhiSc645NegHelCM8;
  GH1* PhiSc675NegHelCM8;

  GH1* PhiSc315PosHelCM1;
  GH1* PhiSc345PosHelCM1;
  GH1* PhiSc375PosHelCM1;
  GH1* PhiSc405PosHelCM1;
  GH1* PhiSc435PosHelCM1;
  GH1* PhiSc465PosHelCM1;
  GH1* PhiSc495PosHelCM1;
  GH1* PhiSc525PosHelCM1;
  GH1* PhiSc555PosHelCM1;
  GH1* PhiSc585PosHelCM1;
  GH1* PhiSc615PosHelCM1;
  GH1* PhiSc645PosHelCM1;
  GH1* PhiSc675PosHelCM1;
  GH1* PhiSc315PosHelCM2;
  GH1* PhiSc345PosHelCM2;
  GH1* PhiSc375PosHelCM2;
  GH1* PhiSc405PosHelCM2;
  GH1* PhiSc435PosHelCM2;
  GH1* PhiSc465PosHelCM2;
  GH1* PhiSc495PosHelCM2;
  GH1* PhiSc525PosHelCM2;
  GH1* PhiSc555PosHelCM2;
  GH1* PhiSc585PosHelCM2;
  GH1* PhiSc615PosHelCM2;
  GH1* PhiSc645PosHelCM2;
  GH1* PhiSc675PosHelCM2;
  GH1* PhiSc315PosHelCM3;
  GH1* PhiSc345PosHelCM3;
  GH1* PhiSc375PosHelCM3;
  GH1* PhiSc405PosHelCM3;
  GH1* PhiSc435PosHelCM3;
  GH1* PhiSc465PosHelCM3;
  GH1* PhiSc495PosHelCM3;
  GH1* PhiSc525PosHelCM3;
  GH1* PhiSc555PosHelCM3;
  GH1* PhiSc585PosHelCM3;
  GH1* PhiSc615PosHelCM3;
  GH1* PhiSc645PosHelCM3;
  GH1* PhiSc675PosHelCM3;
  GH1* PhiSc315PosHelCM4;
  GH1* PhiSc345PosHelCM4;
  GH1* PhiSc375PosHelCM4;
  GH1* PhiSc405PosHelCM4;
  GH1* PhiSc435PosHelCM4;
  GH1* PhiSc465PosHelCM4;
  GH1* PhiSc495PosHelCM4;
  GH1* PhiSc525PosHelCM4;
  GH1* PhiSc555PosHelCM4;
  GH1* PhiSc585PosHelCM4;
  GH1* PhiSc615PosHelCM4;
  GH1* PhiSc645PosHelCM4;
  GH1* PhiSc675PosHelCM4;
  GH1* PhiSc315PosHelCM5;
  GH1* PhiSc345PosHelCM5;
  GH1* PhiSc375PosHelCM5;
  GH1* PhiSc405PosHelCM5;
  GH1* PhiSc435PosHelCM5;
  GH1* PhiSc465PosHelCM5;
  GH1* PhiSc495PosHelCM5;
  GH1* PhiSc525PosHelCM5;
  GH1* PhiSc555PosHelCM5;
  GH1* PhiSc585PosHelCM5;
  GH1* PhiSc615PosHelCM5;
  GH1* PhiSc645PosHelCM5;
  GH1* PhiSc675PosHelCM5;
  GH1* PhiSc315PosHelCM6;
  GH1* PhiSc345PosHelCM6;
  GH1* PhiSc375PosHelCM6;
  GH1* PhiSc405PosHelCM6;
  GH1* PhiSc435PosHelCM6;
  GH1* PhiSc465PosHelCM6;
  GH1* PhiSc495PosHelCM6;
  GH1* PhiSc525PosHelCM6;
  GH1* PhiSc555PosHelCM6;
  GH1* PhiSc585PosHelCM6;
  GH1* PhiSc615PosHelCM6;
  GH1* PhiSc645PosHelCM6;
  GH1* PhiSc675PosHelCM6;
  GH1* PhiSc315PosHelCM7;
  GH1* PhiSc345PosHelCM7;
  GH1* PhiSc375PosHelCM7;
  GH1* PhiSc405PosHelCM7;
  GH1* PhiSc435PosHelCM7;
  GH1* PhiSc465PosHelCM7;
  GH1* PhiSc495PosHelCM7;
  GH1* PhiSc525PosHelCM7;
  GH1* PhiSc555PosHelCM7;
  GH1* PhiSc585PosHelCM7;
  GH1* PhiSc615PosHelCM7;
  GH1* PhiSc645PosHelCM7;
  GH1* PhiSc675PosHelCM7;
  GH1* PhiSc315PosHelCM8;
  GH1* PhiSc345PosHelCM8;
  GH1* PhiSc375PosHelCM8;
  GH1* PhiSc405PosHelCM8;
  GH1* PhiSc435PosHelCM8;
  GH1* PhiSc465PosHelCM8;
  GH1* PhiSc495PosHelCM8;
  GH1* PhiSc525PosHelCM8;
  GH1* PhiSc555PosHelCM8;
  GH1* PhiSc585PosHelCM8;
  GH1* PhiSc615PosHelCM8;
  GH1* PhiSc645PosHelCM8;
  GH1* PhiSc675PosHelCM8;

  GH2* E_dE;
  GH2* KinEp_dE;
  GH2* ThetaScPhiSc;
  GH2* EpCorr_KinEp;
  GH2* PhiDiffThetaSc;
  GH2* PhinDiffWCZRec;
  GH2* ThetaDiffPhiDiff;

  GH2* DeutKinPiKin;

  GH1* ThetanDist;
  GH1* ThetanRecDist;
  GH1* ThetanDiffDist;
  GH2* ThetanDiffZp;

  GH1* ThetanCorrDist;
  GH1* ThetanCorrDiffDist;
  GH1* ThetanCorrRecDiffDist;
  GH2* ThetanCorrDiffZp;

  GH1* ThetaRecPiDiff;
  GH2* ThetanThetaRecPi;
  GH2* ThetanThetaRecPiDiff;

  GH1* ThetaRecPDiff;
  GH2* ThetanThetaRecP;
  GH2* ThetanThetaRecPDiff;

  GH1* ClosestApproach;
  GH1* POCAr;
  GH1* ScatterVertexZ;
  GH2* ScatterVertexZr;
  GH2* ScatterVertexXY;
  GH3* ScatterVertex;

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

      // lightweight structure for linking to fitter
    struct FitParticle{
        void SetFromVector(const TLorentzVector& p_) {
            Ek = p_.E()-p_.M();
            Theta = p_.Theta();
            Phi = p_.Phi();
        }

        static TLorentzVector Make(const std::vector<double>& EkThetaPhi,
        const Double_t m);
        static TLorentzVector Make(const FitParticle& p, const Double_t m) {
            return Make(std::vector<double>{p.Ek, p.Theta, p.Phi}, m);
        }

    std::vector<double*> Link() {
        return {std::addressof(Ek),
        std::addressof(Theta),
        std::addressof(Phi)};
    }
        std::vector<double*> LinkSigma() {
        return {std::addressof(Ek_Sigma),
        std::addressof(Theta_Sigma),
        std::addressof(Phi_Sigma)};
    }

    std::vector<APLCON::Variable_Settings_t> LinkSettings()
    {
        return{Ek_Setting, Theta_Setting, Phi_Setting};
    }

    void Smear(std::vector<double> unc , int particle);

    void APLCONSettings();

    double Ek;
    double Ek_Sigma;
    APLCON::Variable_Settings_t Ek_Setting;
    double Theta;
    double Theta_Sigma;
    APLCON::Variable_Settings_t Theta_Setting;
    double Phi;
    double Phi_Sigma;
    APLCON::Variable_Settings_t Phi_Setting;

    bool isCB;

    private:
        static std::default_random_engine generator;

    };

    FitParticle beamF;
    FitParticle protonF;
    FitParticle neutronF;
    //APLCON kinfit("EMcons", settings);

public:

    PNeutPol_Polarimeter_Circ();
    virtual ~PNeutPol_Polarimeter_Circ();
    virtual Bool_t  Init();
    TCutG* OpenCutFile(Char_t* filename, Char_t* cutname);
    TLorentzVector CNeutron4VectorCorr(Double_t ZVert, TLorentzVector n4Vector, Double_t nE, Double_t MagP, Double_t nMass, Double_t nPhi);
    TLorentzVector CProton4VectorKin(Double_t KinE, Double_t Theta, Double_t Phi);
    TLorentzVector CNeutron4VectorKin(TLorentzVector ProtonKinVector);
    TLorentzVector CPion4VectorKin(TLorentzVector ProtonKinVector);
    Double_t LabAngles();
    void FillHists();

};
#endif
