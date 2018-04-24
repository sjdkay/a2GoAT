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
    double_t ln;
    double_t Thp;
    double_t ThpRad;
    double_t Thn;
    double_t Php;
    double_t PhpRad;
    double_t Phn;
    double_t PhnRad;
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
    double_t ScattThetaRad;
    double_t ScattPhiRad;
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
    double_t DOCA;
    double_t POCAx;
    double_t POCAy;
    double_t POCAz;
    double_t r;
    double_t tdif;
    double_t MWPC0pE;
    double_t MWPC1pE;
    double_t MWPC0nE;
    double_t MWPC1nE;
    double_t AngularDiffr;
    double_t MWPCPLp;
    double_t Wgt;

    double_t ELow;
    double_t EHigh;
    double_t EWidth;
    double_t CosThetaLow;
    double_t CosThetaHigh;
    double_t CosThetaWidth;
    double_t PhiScLow;
    double_t PhiScHigh;
    double_t PhiScWidth;

    double RtmpMass,RtmpMom;
    double chi2;
    TLorentzVector* Pp1;
    TLorentzVector* Pp2;
    TLorentzVector* Pbeam;
    TLorentzVector* Pp1C;
    TLorentzVector* Pp2C;
    TLorentzVector* PbeamC;
    TLorentzVector* PtargetC;

    double_t WC1Phip;
    double_t WC2Phip;
    double_t WC1Phin;
    double_t WC2Phin;
    double_t WCPhiDiffp;
    double_t WCPhiDiffn;

    Bool_t MCData;
    Bool_t Proton1;
    Bool_t Proton2;
    Bool_t BeamHelicity;

    TLorentzVector B;
    TLorentzVector GVp;
    TLorentzVector GVn;
    TLorentzVector GVpCorr;
    TLorentzVector GVnCorr;
    TLorentzVector Gamma;
    TLorentzVector Deut;
    TLorentzVector P4Vect;
    TLorentzVector N4Vect;
    TLorentzVector N4VectCorr;
    TLorentzVector RecNeutronEpCorr;

    TLorentzVector pKin;
    TLorentzVector nKin;
    TLorentzVector pKinB;
    TLorentzVector nKinB;
    TVector3 pKin3;
    TVector3 nKin3;

    TVector3 b;
    TVector3 pVertex;
    TVector3 nVertex;
    TVector3 DOCAVertex1;
    TVector3 DOCAVertex2;
    TVector3 POCA;
    TVector3 GVpCorr3;
    TVector3 GVnCorr3;
    TVector3 GVn3;
    TVector3 GVn3Rec;
    TVector3 WC13Vectp;
    TVector3 WC23Vectp;
    TVector3 WC13Vectn;
    TVector3 WC23Vectn;
    TVector3 WCDirp;
    TVector3 RecoilVector;

    // Variables for calculating P vertex from MWPC hit positions
    double_t num0;
    double_t denum0;
    double_t tsk;
    TVector3 MWPCpDir0;
    TVector3 MWPCpDir1;
    TVector3 MWPCpPerp;
    TVector3 MWPCpCalcVertex;
    double_t ZpMWPC;

    TH1D*	time;
    TH1D*	time_cut;

    TH1D* Eg;
    TH1D* ThetaSc;
    TH1D* RPoca;
    TH1D* Thetap;
    TH1D* Phip;
    TH1D* EProton;
    TH1D* MMProton;
    TH1D* PhiSc;
    TH1D* PhiSc320;
    TH1D* PhiSc360;
    TH1D* PhiSc400;
    TH1D* PhiSc440;
    TH1D* PhiSc480;
    TH1D* PhiSc520;
    TH1D* PhiSc560;
    TH1D* PhiSc600;
    TH1D* PhiSc640;
    TH1D* PhiSc680;

    TH2D* EdE;
    TH2D* EdEPrompt;
    TH2D* EdERandom;

    TH2D* ThetaScThetap;
    TH2D* ThetaScThetapPrompt;
    TH2D* ThetaScThetapRandom;

    TH2D* PhiScThetap;
    TH2D* PhiScThetapPrompt;
    TH2D* PhiScThetapRandom;

    TH2D* rPocaThetap;
    TH2D* rPocaThetapPrompt;
    TH2D* rPocaThetapRandom;

    TH2D* ThetaDiffPhiDiff;
    TH2D* ThetaDiffPhiDiffPrompt;
    TH2D* ThetaDiffPhiDiffRandom;

    TH1D* EgPrompt;
    TH1D* ThetaScPrompt;
    TH1D* RPocaPrompt;
    TH1D* ThetapPrompt;
    TH1D* PhipPrompt;
    TH1D* EProtonPrompt;
    TH1D* MMProtonPrompt;
    TH1D* EgRandom;
    TH1D* ThetaScRandom;
    TH1D* RPocaRandom;
    TH1D* ThetapRandom;
    TH1D* PhipRandom;
    TH1D* EProtonRandom;
    TH1D* MMProtonRandom;
    TH1D* PhiScPrompt;
    TH1D* PhiScRandom;
    TH1D* PhiSc320Prompt;
    TH1D* PhiSc360Prompt;
    TH1D* PhiSc400Prompt;
    TH1D* PhiSc440Prompt;
    TH1D* PhiSc480Prompt;
    TH1D* PhiSc520Prompt;
    TH1D* PhiSc560Prompt;
    TH1D* PhiSc600Prompt;
    TH1D* PhiSc640Prompt;
    TH1D* PhiSc680Prompt;
    TH1D* PhiSc320Random;
    TH1D* PhiSc360Random;
    TH1D* PhiSc400Random;
    TH1D* PhiSc440Random;
    TH1D* PhiSc480Random;
    TH1D* PhiSc520Random;
    TH1D* PhiSc560Random;
    TH1D* PhiSc600Random;
    TH1D* PhiSc640Random;
    TH1D* PhiSc680Random;

    TH1D* MM_Eg_CM_Tot[7][5];
    TH1D* MM_Eg_CM_Tot_Prompt[7][5];
    TH1D* MM_Eg_CM_Tot_Random[7][5];
    TH1D* MM_Eg_CM_Binned[7][5][10];
    TH1D* MM_Eg_CM_Binned_Prompt[7][5][10];
    TH1D* MM_Eg_CM_Binned_Random[7][5][10];

    TH1D* PhiScPosHel[7][5];
    TH1D* PhiScPosHelPrompt[7][5];
    TH1D* PhiScPosHelRandom[7][5];
    TH1D* PhiScNegHel[7][5];
    TH1D* PhiScNegHelPrompt[7][5];
    TH1D* PhiScNegHelRandom[7][5];

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
    Double_t LabAngles();
    void FillHists();
    void FillHistsPreMMCut();
    void BGSub();

};
#endif
