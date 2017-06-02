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
  TVector3 GVpCorr3;
  TVector3 GVnCorr3;
  TVector3 GVn3Rec;
  TVector3 WC3Vectp;
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

  GH2* E_dE;
  GH2* KinEp_dE;
  GH2* ThetaScPhiSc;
  GH2* E_KinEp;
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
