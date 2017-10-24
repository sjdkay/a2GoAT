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

    //APLCON DOCA

    double_t DOCAApl;
    double_t POCAxApl;
    double_t POCAyApl;
    double_t POCAzApl;
    TVector3 DOCAVertex1Apl;
    TVector3 DOCAVertex2Apl;
    TVector3 POCAApl;
    double_t rApl;
    GH1* ClosestApproachApl;
    GH1* POCArApl;
    GH2* POCArPOCArAPL;

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

       GH1* Eg;
    GH1* PhiDet;
    GH1* PhiRec;
    GH1* ThetaSc;
    GH1* PhiSc;
    GH1* MMpEpCorrected;
    GH1* ZpDist;
    GH1* ThetanDist;
    GH1* ThetanCMDist;

    GH2* DeutKinPiKin;
    GH2* E_dE;
    GH2* ThetaScPhiSc;
    GH2* PhiPp1Phip;
    GH2* EdEMWPCp;
    GH2* EdEMWPCn;

    GH1* ClosestApproach;
    GH1* POCAr;
    GH1* ScatterVertexZ;
    GH2* ScatterVertexZr;
    GH2* ScatterVertexXY;
    GH3* ScatterVertex;
    GH2* POCArPhiSc;

    GH1* ThetapCorrDiff;
    GH1* PhipCorrDiff;
    GH1* ThetaDiff;
    GH1* PhiDiff;
    GH2* PhiDiffThetaDiff;

    GH2* PhiScEg;
    GH2* PhiScEp;
    GH2* PhiScThetan;

    GH2* EMWPCnPhiSc;

    GH1* ZpVertexDiff;

    GH1* MMp200300;
    GH1* MMp300400;
    GH1* MMp400500;
    GH1* MMp500600;
    GH1* MMp600700;
    GH1* MMp700800;
    GH1* MMp800900;

    GH1* PhiSc320;
    GH1* PhiSc360;
    GH1* PhiSc400;
    GH1* PhiSc440;
    GH1* PhiSc480;
    GH1* PhiSc520;
    GH1* PhiSc560;
    GH1* PhiSc600;
    GH1* PhiSc640;
    GH1* PhiSc680;

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

        GH1* PhiSc265NegHelCM1;
    GH1* PhiSc335NegHelCM1;
    GH1* PhiSc405NegHelCM1;
    GH1* PhiSc475NegHelCM1;
    GH1* PhiSc545NegHelCM1;
    GH1* PhiSc615NegHelCM1;
    GH1* PhiSc685NegHelCM1;
    GH1* PhiSc265NegHelCM2;
    GH1* PhiSc335NegHelCM2;
    GH1* PhiSc405NegHelCM2;
    GH1* PhiSc475NegHelCM2;
    GH1* PhiSc545NegHelCM2;
    GH1* PhiSc615NegHelCM2;
    GH1* PhiSc685NegHelCM2;
    GH1* PhiSc265NegHelCM3;
    GH1* PhiSc335NegHelCM3;
    GH1* PhiSc405NegHelCM3;
    GH1* PhiSc475NegHelCM3;
    GH1* PhiSc545NegHelCM3;
    GH1* PhiSc615NegHelCM3;
    GH1* PhiSc685NegHelCM3;
    GH1* PhiSc265NegHelCM4;
    GH1* PhiSc335NegHelCM4;
    GH1* PhiSc405NegHelCM4;
    GH1* PhiSc475NegHelCM4;
    GH1* PhiSc545NegHelCM4;
    GH1* PhiSc615NegHelCM4;
    GH1* PhiSc685NegHelCM4;
    GH1* PhiSc265NegHelCM5;
    GH1* PhiSc335NegHelCM5;
    GH1* PhiSc405NegHelCM5;
    GH1* PhiSc475NegHelCM5;
    GH1* PhiSc545NegHelCM5;
    GH1* PhiSc615NegHelCM5;
    GH1* PhiSc685NegHelCM5;
    GH1* PhiSc265NegHelCM6;
    GH1* PhiSc335NegHelCM6;
    GH1* PhiSc405NegHelCM6;
    GH1* PhiSc475NegHelCM6;
    GH1* PhiSc545NegHelCM6;
    GH1* PhiSc615NegHelCM6;
    GH1* PhiSc685NegHelCM6;
    GH1* PhiSc265NegHelCM7;
    GH1* PhiSc335NegHelCM7;
    GH1* PhiSc405NegHelCM7;
    GH1* PhiSc475NegHelCM7;
    GH1* PhiSc545NegHelCM7;
    GH1* PhiSc615NegHelCM7;
    GH1* PhiSc685NegHelCM7;
    GH1* PhiSc265NegHelCM8;
    GH1* PhiSc335NegHelCM8;
    GH1* PhiSc405NegHelCM8;
    GH1* PhiSc475NegHelCM8;
    GH1* PhiSc545NegHelCM8;
    GH1* PhiSc615NegHelCM8;
    GH1* PhiSc685NegHelCM8;

    GH1* PhiSc265PosHelCM1;
    GH1* PhiSc335PosHelCM1;
    GH1* PhiSc405PosHelCM1;
    GH1* PhiSc475PosHelCM1;
    GH1* PhiSc545PosHelCM1;
    GH1* PhiSc615PosHelCM1;
    GH1* PhiSc685PosHelCM1;
    GH1* PhiSc265PosHelCM2;
    GH1* PhiSc335PosHelCM2;
    GH1* PhiSc405PosHelCM2;
    GH1* PhiSc475PosHelCM2;
    GH1* PhiSc545PosHelCM2;
    GH1* PhiSc615PosHelCM2;
    GH1* PhiSc685PosHelCM2;
    GH1* PhiSc265PosHelCM3;
    GH1* PhiSc335PosHelCM3;
    GH1* PhiSc405PosHelCM3;
    GH1* PhiSc475PosHelCM3;
    GH1* PhiSc545PosHelCM3;
    GH1* PhiSc615PosHelCM3;
    GH1* PhiSc685PosHelCM3;
    GH1* PhiSc265PosHelCM4;
    GH1* PhiSc335PosHelCM4;
    GH1* PhiSc405PosHelCM4;
    GH1* PhiSc475PosHelCM4;
    GH1* PhiSc545PosHelCM4;
    GH1* PhiSc615PosHelCM4;
    GH1* PhiSc685PosHelCM4;
    GH1* PhiSc265PosHelCM5;
    GH1* PhiSc335PosHelCM5;
    GH1* PhiSc405PosHelCM5;
    GH1* PhiSc475PosHelCM5;
    GH1* PhiSc545PosHelCM5;
    GH1* PhiSc615PosHelCM5;
    GH1* PhiSc685PosHelCM5;
    GH1* PhiSc265PosHelCM6;
    GH1* PhiSc335PosHelCM6;
    GH1* PhiSc405PosHelCM6;
    GH1* PhiSc475PosHelCM6;
    GH1* PhiSc545PosHelCM6;
    GH1* PhiSc615PosHelCM6;
    GH1* PhiSc685PosHelCM6;
    GH1* PhiSc265PosHelCM7;
    GH1* PhiSc335PosHelCM7;
    GH1* PhiSc405PosHelCM7;
    GH1* PhiSc475PosHelCM7;
    GH1* PhiSc545PosHelCM7;
    GH1* PhiSc615PosHelCM7;
    GH1* PhiSc685PosHelCM7;
    GH1* PhiSc265PosHelCM8;
    GH1* PhiSc335PosHelCM8;
    GH1* PhiSc405PosHelCM8;
    GH1* PhiSc475PosHelCM8;
    GH1* PhiSc545PosHelCM8;
    GH1* PhiSc615PosHelCM8;
    GH1* PhiSc685PosHelCM8;

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
