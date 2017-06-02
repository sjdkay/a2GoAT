#ifndef __PNeutPol_Polarimeter_Lin_h__
#define __PNeutPol_Polarimeter_Lin_h__

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

class	PNeutPol_Polarimeter_Lin : public PPhysics
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

  Bool_t nBanana;
  Bool_t Proton1;
  Bool_t Proton2;

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
  TVector3 GVpCorr3;
  TVector3 GVnCorr3;
  TVector3 GVn3Rec;
  TVector3 WC13Vectn;
  TVector3 WC23Vectn;
  TVector3 P3Vect;
  TVector3 N3Vect;
  TVector3 RecProtonEpCorr3;
  TVector3 RecNeutronEpCorr3;

  TH1D*	time;
  TH1D*	time_cut;

  GH1* EkSum;
  GH1* Eg;
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
  GH1* OAngle200400;

  GH1* ZpDist;
  GH1* ZpPhiScatNeg180;
  GH1* ZpPhiScat0;
  GH1* ZpPhiScatPos180;

  GH1* ThetaSc;
  GH1* PhiSc;
  GH1* PhiSc410;
  GH1* PhiSc430;
  GH1* PhiSc450;
  GH1* PhiSc470;
  GH1* PhiSc490;
  GH1* PhiSc510;
  GH1* PhiSc530;
  GH1* PhiSc550;
  GH1* PhiSc570;
  GH1* PhiSc590;
  GH1* PhiSc610;
  GH1* PhiSc630;

  GH1* Phip410CM1;
  GH1* Phip430CM1;
  GH1* Phip450CM1;
  GH1* Phip470CM1;
  GH1* Phip490CM1;
  GH1* Phip510CM1;
  GH1* Phip530CM1;
  GH1* Phip550CM1;
  GH1* Phip570CM1;
  GH1* Phip590CM1;
  GH1* Phip610CM1;
  GH1* Phip630CM1;

  GH1* Phip410CM2;
  GH1* Phip430CM2;
  GH1* Phip450CM2;
  GH1* Phip470CM2;
  GH1* Phip490CM2;
  GH1* Phip510CM2;
  GH1* Phip530CM2;
  GH1* Phip550CM2;
  GH1* Phip570CM2;
  GH1* Phip590CM2;
  GH1* Phip610CM2;
  GH1* Phip630CM2;

  GH1* Phip410CM3;
  GH1* Phip430CM3;
  GH1* Phip450CM3;
  GH1* Phip470CM3;
  GH1* Phip490CM3;
  GH1* Phip510CM3;
  GH1* Phip530CM3;
  GH1* Phip550CM3;
  GH1* Phip570CM3;
  GH1* Phip590CM3;
  GH1* Phip610CM3;
  GH1* Phip630CM3;

  GH1* Phip410CM4;
  GH1* Phip430CM4;
  GH1* Phip450CM4;
  GH1* Phip470CM4;
  GH1* Phip490CM4;
  GH1* Phip510CM4;
  GH1* Phip530CM4;
  GH1* Phip550CM4;
  GH1* Phip570CM4;
  GH1* Phip590CM4;
  GH1* Phip610CM4;
  GH1* Phip630CM4;

  GH1* Phip410CM5;
  GH1* Phip430CM5;
  GH1* Phip450CM5;
  GH1* Phip470CM5;
  GH1* Phip490CM5;
  GH1* Phip510CM5;
  GH1* Phip530CM5;
  GH1* Phip550CM5;
  GH1* Phip570CM5;
  GH1* Phip590CM5;
  GH1* Phip610CM5;
  GH1* Phip630CM5;

  GH1* Phip410CM6;
  GH1* Phip430CM6;
  GH1* Phip450CM6;
  GH1* Phip470CM6;
  GH1* Phip490CM6;
  GH1* Phip510CM6;
  GH1* Phip530CM6;
  GH1* Phip550CM6;
  GH1* Phip570CM6;
  GH1* Phip590CM6;
  GH1* Phip610CM6;
  GH1* Phip630CM6;

  GH1* Phip410CM7;
  GH1* Phip430CM7;
  GH1* Phip450CM7;
  GH1* Phip470CM7;
  GH1* Phip490CM7;
  GH1* Phip510CM7;
  GH1* Phip530CM7;
  GH1* Phip550CM7;
  GH1* Phip570CM7;
  GH1* Phip590CM7;
  GH1* Phip610CM7;
  GH1* Phip630CM7;

  GH1* Phip410CM8;
  GH1* Phip430CM8;
  GH1* Phip450CM8;
  GH1* Phip470CM8;
  GH1* Phip490CM8;
  GH1* Phip510CM8;
  GH1* Phip530CM8;
  GH1* Phip550CM8;
  GH1* Phip570CM8;
  GH1* Phip590CM8;
  GH1* Phip610CM8;
  GH1* Phip630CM8;

  GH1* Phip410CM9;
  GH1* Phip430CM9;
  GH1* Phip450CM9;
  GH1* Phip470CM9;
  GH1* Phip490CM9;
  GH1* Phip510CM9;
  GH1* Phip530CM9;
  GH1* Phip550CM9;
  GH1* Phip570CM9;
  GH1* Phip590CM9;
  GH1* Phip610CM9;
  GH1* Phip630CM9;

  GH1* Phip410CM10;
  GH1* Phip430CM10;
  GH1* Phip450CM10;
  GH1* Phip470CM10;
  GH1* Phip490CM10;
  GH1* Phip510CM10;
  GH1* Phip530CM10;
  GH1* Phip550CM10;
  GH1* Phip570CM10;
  GH1* Phip590CM10;
  GH1* Phip610CM10;
  GH1* Phip630CM10;

  GH2* E_dE;
  GH2* E_dE_Cut;
  GH2* KinEp_dE;
  GH2* KinEp_dE_GoodCut;
  GH2* ThetaScPhiSc;
  GH2* E_KinEp;
  GH2* E_KinEpCut;
  GH2* PhinDiffWCZRec;
  GH2* PhinDiffWCZRec_KinCut;

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

  GH2* DeutKinPiKin;

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

public:

    PNeutPol_Polarimeter_Lin();
    virtual ~PNeutPol_Polarimeter_Lin();
    virtual Bool_t  Init();
    TCutG* OpenCutFile(Char_t* filename, Char_t* cutname);
    TLorentzVector LNeutron4VectorCorr(Double_t ZVert, TLorentzVector n4Vector, Double_t nE, Double_t MagP, Double_t nMass, Double_t nPhi);
    TLorentzVector LProton4VectorKin(Double_t KinE, Double_t Theta, Double_t Phi);
    TLorentzVector LNeutron4VectorKin(TLorentzVector ProtonKinVector);
    TLorentzVector LPion4VectorKin(TLorentzVector ProtonKinVector);
    Double_t LabAngles();
    void FillHists();

};
#endif
