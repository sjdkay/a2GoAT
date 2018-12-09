#ifndef __PNeutPol_Polarimeter_Lin_NoScatt_h__
#define __PNeutPol_Polarimeter_Lin_NoScatt_h__

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

class	PNeutPol_Polarimeter_Lin_NoScatt : public PPhysics
{
private:

    TH1*	TaggerAccScal;

    Int_t k;
    Int_t NP;
    Int_t NPi;
    Int_t NRoo;
    Int_t NTag;
    Int_t NTrack;
    Int_t EventNum;

    Int_t Detectors1;
    Int_t Detectors2;
    Int_t DetectorsSum;
    Int_t EventNumber;

    double_t PromptLow;
    double_t PromptHigh;
    double_t RandomLow1;
    double_t RandomHigh1;
    double_t RandomLow2;
    double_t RandomHigh2;
    double_t PvRratio;

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
    double_t WCZnRec;
    double_t ThetapCM;
    double_t ThetanCM;
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
    double_t PhinCorr;

    Bool_t MCData;
    Bool_t Proton1;
    Bool_t Proton2;
    Bool_t BeamHelicity;
    Bool_t ExecuteCut;

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
    TVector3 P3Vect;
    TVector3 N3Vect;
    TVector3 RecProtonEpCorr3;
    TVector3 RecNeutronEpCorr3;

    TLorentzVector pKin;
    TLorentzVector nKin;
    TLorentzVector pKinB;
    TLorentzVector nKinB;
    TVector3 pKin3;
    TVector3 nKin3;

    TH1D*	time;
    TH1D*	time_cut;

    double_t ELow;
    double_t EHigh;
    double_t EWidth;
    double_t ELow2;
    double_t EHigh2;
    double_t EWidth2;
    double_t CosThetaLow;
    double_t CosThetaHigh;
    double_t CosThetaWidth;
    double_t CosThetaLow2;
    double_t CosThetaHigh2;
    double_t CosThetaWidth2;

    GH1* Eg;
    GH1* PhiDiffDist;
    GH1* OAngle;
    GH1* MMpEpCorrected;
    GH1* ZpDist;
    GH1* ThetanDist;
    GH1* ThetanDiffDist;

    GH2* E_dEKin;
    GH2* DeutKinPiKin;

    TH1D* Eg2;
    TH1D* EgPrompt;
    TH1D* EgRandom;
    TH2D* EdE;
    TH2D* EdEPrompt;
    TH2D* EdERandom;
    TH1D* Thetap;
    TH1D* ThetapPrompt;
    TH1D* ThetapRandom;

    TH1D* PhipPosHel[7][5];
    TH1D* PhipPosHelPrompt[7][5];
    TH1D* PhipPosHelRandom[7][5];
    TH1D* PhipNegHel[7][5];
    TH1D* PhipNegHelPrompt[7][5];
    TH1D* PhipNegHelRandom[7][5];

    TH1D* PhipSet[21][20];
    TH1D* PhipSetPrompt[21][20];
    TH1D* PhipSetRandom[21][20];

    TH1D* MMpEgamma[12];
    TH1D* MMpEgammaPrompt[12];
    TH1D* MMpEgammaRandom[12];

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

    PNeutPol_Polarimeter_Lin_NoScatt();
    virtual ~PNeutPol_Polarimeter_Lin_NoScatt();
    virtual Bool_t  Init();
    TCutG* OpenCutFile(Char_t* filename, Char_t* cutname);
    TLorentzVector LNeutron4VectorCorr(Double_t ZVert, TLorentzVector n4Vector, Double_t nE, Double_t MagP, Double_t nMass, Double_t nPhi);
    TLorentzVector LProton4VectorKin(Double_t KinE, Double_t Theta, Double_t Phi);
    TLorentzVector LNeutron4VectorKin(TLorentzVector ProtonKinVector);
    TLorentzVector LPion4VectorKin(TLorentzVector ProtonKinVector);
    Double_t LabAngles();
    void FillHists();
    void BGSub();

};
#endif
