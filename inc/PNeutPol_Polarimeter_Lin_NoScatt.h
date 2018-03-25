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

    Bool_t MCData;
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

    GH1* Eg;
    GH1* OAngle;
    GH1* MMpEpCorrected;
    GH1* ZpDist;
    GH1* ThetanDist;

    GH2* E_dE;
    GH2* DeutKinPiKin;

    GH1* MMp200300;
    GH1* MMp300400;
    GH1* MMp400500;
    GH1* MMp500600;
    GH1* MMp600700;
    GH1* MMp700800;
    GH1* MMp800900;

    GH1* Phip415CM1;
    GH1* Phip425CM1;
    GH1* Phip435CM1;
    GH1* Phip445CM1;
    GH1* Phip455CM1;
    GH1* Phip465CM1;
    GH1* Phip475CM1;
    GH1* Phip485CM1;
    GH1* Phip495CM1;
    GH1* Phip505CM1;
    GH1* Phip515CM1;
    GH1* Phip525CM1;
    GH1* Phip535CM1;
    GH1* Phip545CM1;
    GH1* Phip555CM1;
    GH1* Phip565CM1;
    GH1* Phip575CM1;
    GH1* Phip585CM1;
    GH1* Phip595CM1;
    GH1* Phip605CM1;
    GH1* Phip615CM1;

    GH1* Phip415CM2;
    GH1* Phip425CM2;
    GH1* Phip435CM2;
    GH1* Phip445CM2;
    GH1* Phip455CM2;
    GH1* Phip465CM2;
    GH1* Phip475CM2;
    GH1* Phip485CM2;
    GH1* Phip495CM2;
    GH1* Phip505CM2;
    GH1* Phip515CM2;
    GH1* Phip525CM2;
    GH1* Phip535CM2;
    GH1* Phip545CM2;
    GH1* Phip555CM2;
    GH1* Phip565CM2;
    GH1* Phip575CM2;
    GH1* Phip585CM2;
    GH1* Phip595CM2;
    GH1* Phip605CM2;
    GH1* Phip615CM2;

    GH1* Phip415CM3;
    GH1* Phip425CM3;
    GH1* Phip435CM3;
    GH1* Phip445CM3;
    GH1* Phip455CM3;
    GH1* Phip465CM3;
    GH1* Phip475CM3;
    GH1* Phip485CM3;
    GH1* Phip495CM3;
    GH1* Phip505CM3;
    GH1* Phip515CM3;
    GH1* Phip525CM3;
    GH1* Phip535CM3;
    GH1* Phip545CM3;
    GH1* Phip555CM3;
    GH1* Phip565CM3;
    GH1* Phip575CM3;
    GH1* Phip585CM3;
    GH1* Phip595CM3;
    GH1* Phip605CM3;
    GH1* Phip615CM3;

    GH1* Phip415CM4;
    GH1* Phip425CM4;
    GH1* Phip435CM4;
    GH1* Phip445CM4;
    GH1* Phip455CM4;
    GH1* Phip465CM4;
    GH1* Phip475CM4;
    GH1* Phip485CM4;
    GH1* Phip495CM4;
    GH1* Phip505CM4;
    GH1* Phip515CM4;
    GH1* Phip525CM4;
    GH1* Phip535CM4;
    GH1* Phip545CM4;
    GH1* Phip555CM4;
    GH1* Phip565CM4;
    GH1* Phip575CM4;
    GH1* Phip585CM4;
    GH1* Phip595CM4;
    GH1* Phip605CM4;
    GH1* Phip615CM4;

    GH1* Phip415CM5;
    GH1* Phip425CM5;
    GH1* Phip435CM5;
    GH1* Phip445CM5;
    GH1* Phip455CM5;
    GH1* Phip465CM5;
    GH1* Phip475CM5;
    GH1* Phip485CM5;
    GH1* Phip495CM5;
    GH1* Phip505CM5;
    GH1* Phip515CM5;
    GH1* Phip525CM5;
    GH1* Phip535CM5;
    GH1* Phip545CM5;
    GH1* Phip555CM5;
    GH1* Phip565CM5;
    GH1* Phip575CM5;
    GH1* Phip585CM5;
    GH1* Phip595CM5;
    GH1* Phip605CM5;
    GH1* Phip615CM5;

    GH1* Phip415CM6;
    GH1* Phip425CM6;
    GH1* Phip435CM6;
    GH1* Phip445CM6;
    GH1* Phip455CM6;
    GH1* Phip465CM6;
    GH1* Phip475CM6;
    GH1* Phip485CM6;
    GH1* Phip495CM6;
    GH1* Phip505CM6;
    GH1* Phip515CM6;
    GH1* Phip525CM6;
    GH1* Phip535CM6;
    GH1* Phip545CM6;
    GH1* Phip555CM6;
    GH1* Phip565CM6;
    GH1* Phip575CM6;
    GH1* Phip585CM6;
    GH1* Phip595CM6;
    GH1* Phip605CM6;
    GH1* Phip615CM6;

    GH1* Phip415CM7;
    GH1* Phip425CM7;
    GH1* Phip435CM7;
    GH1* Phip445CM7;
    GH1* Phip455CM7;
    GH1* Phip465CM7;
    GH1* Phip475CM7;
    GH1* Phip485CM7;
    GH1* Phip495CM7;
    GH1* Phip505CM7;
    GH1* Phip515CM7;
    GH1* Phip525CM7;
    GH1* Phip535CM7;
    GH1* Phip545CM7;
    GH1* Phip555CM7;
    GH1* Phip565CM7;
    GH1* Phip575CM7;
    GH1* Phip585CM7;
    GH1* Phip595CM7;
    GH1* Phip605CM7;
    GH1* Phip615CM7;

    GH1* Phip415CM8;
    GH1* Phip425CM8;
    GH1* Phip435CM8;
    GH1* Phip445CM8;
    GH1* Phip455CM8;
    GH1* Phip465CM8;
    GH1* Phip475CM8;
    GH1* Phip485CM8;
    GH1* Phip495CM8;
    GH1* Phip505CM8;
    GH1* Phip515CM8;
    GH1* Phip525CM8;
    GH1* Phip535CM8;
    GH1* Phip545CM8;
    GH1* Phip555CM8;
    GH1* Phip565CM8;
    GH1* Phip575CM8;
    GH1* Phip585CM8;
    GH1* Phip595CM8;
    GH1* Phip605CM8;
    GH1* Phip615CM8;

    GH1* Phip415CM9;
    GH1* Phip425CM9;
    GH1* Phip435CM9;
    GH1* Phip445CM9;
    GH1* Phip455CM9;
    GH1* Phip465CM9;
    GH1* Phip475CM9;
    GH1* Phip485CM9;
    GH1* Phip495CM9;
    GH1* Phip505CM9;
    GH1* Phip515CM9;
    GH1* Phip525CM9;
    GH1* Phip535CM9;
    GH1* Phip545CM9;
    GH1* Phip555CM9;
    GH1* Phip565CM9;
    GH1* Phip575CM9;
    GH1* Phip585CM9;
    GH1* Phip595CM9;
    GH1* Phip605CM9;
    GH1* Phip615CM9;

    GH1* Phip415CM10;
    GH1* Phip425CM10;
    GH1* Phip435CM10;
    GH1* Phip445CM10;
    GH1* Phip455CM10;
    GH1* Phip465CM10;
    GH1* Phip475CM10;
    GH1* Phip485CM10;
    GH1* Phip495CM10;
    GH1* Phip505CM10;
    GH1* Phip515CM10;
    GH1* Phip525CM10;
    GH1* Phip535CM10;
    GH1* Phip545CM10;
    GH1* Phip555CM10;
    GH1* Phip565CM10;
    GH1* Phip575CM10;
    GH1* Phip585CM10;
    GH1* Phip595CM10;
    GH1* Phip605CM10;
    GH1* Phip615CM10;

    TH1D* Eg2;
    TH1D* EgPrompt;
    TH1D* EgRandom;

    TH1D* Phip;
    TH1D* Phip320;
    TH1D* Phip360;
    TH1D* Phip400;
    TH1D* Phip440;
    TH1D* Phip480;
    TH1D* Phip520;
    TH1D* Phip560;
    TH1D* Phip600;
    TH1D* Phip640;
    TH1D* Phip680;

    TH1D* PhipPrompt;
    TH1D* PhipRandom;
    TH1D* Phip320Prompt;
    TH1D* Phip360Prompt;
    TH1D* Phip400Prompt;
    TH1D* Phip440Prompt;
    TH1D* Phip480Prompt;
    TH1D* Phip520Prompt;
    TH1D* Phip560Prompt;
    TH1D* Phip600Prompt;
    TH1D* Phip640Prompt;
    TH1D* Phip680Prompt;
    TH1D* Phip320Random;
    TH1D* Phip360Random;
    TH1D* Phip400Random;
    TH1D* Phip440Random;
    TH1D* Phip480Random;
    TH1D* Phip520Random;
    TH1D* Phip560Random;
    TH1D* Phip600Random;
    TH1D* Phip640Random;
    TH1D* Phip680Random;

    TH1D* Phip265NegHelCM1;
    TH1D* Phip335NegHelCM1;
    TH1D* Phip405NegHelCM1;
    TH1D* Phip475NegHelCM1;
    TH1D* Phip545NegHelCM1;
    TH1D* Phip615NegHelCM1;
    TH1D* Phip685NegHelCM1;
    TH1D* Phip265NegHelCM2;
    TH1D* Phip335NegHelCM2;
    TH1D* Phip405NegHelCM2;
    TH1D* Phip475NegHelCM2;
    TH1D* Phip545NegHelCM2;
    TH1D* Phip615NegHelCM2;
    TH1D* Phip685NegHelCM2;
    TH1D* Phip265NegHelCM3;
    TH1D* Phip335NegHelCM3;
    TH1D* Phip405NegHelCM3;
    TH1D* Phip475NegHelCM3;
    TH1D* Phip545NegHelCM3;
    TH1D* Phip615NegHelCM3;
    TH1D* Phip685NegHelCM3;
    TH1D* Phip265NegHelCM4;
    TH1D* Phip335NegHelCM4;
    TH1D* Phip405NegHelCM4;
    TH1D* Phip475NegHelCM4;
    TH1D* Phip545NegHelCM4;
    TH1D* Phip615NegHelCM4;
    TH1D* Phip685NegHelCM4;
    TH1D* Phip265NegHelCM5;
    TH1D* Phip335NegHelCM5;
    TH1D* Phip405NegHelCM5;
    TH1D* Phip475NegHelCM5;
    TH1D* Phip545NegHelCM5;
    TH1D* Phip615NegHelCM5;
    TH1D* Phip685NegHelCM5;
    TH1D* Phip265NegHelCM6;
    TH1D* Phip335NegHelCM6;
    TH1D* Phip405NegHelCM6;
    TH1D* Phip475NegHelCM6;
    TH1D* Phip545NegHelCM6;
    TH1D* Phip615NegHelCM6;
    TH1D* Phip685NegHelCM6;
    TH1D* Phip265NegHelCM7;
    TH1D* Phip335NegHelCM7;
    TH1D* Phip405NegHelCM7;
    TH1D* Phip475NegHelCM7;
    TH1D* Phip545NegHelCM7;
    TH1D* Phip615NegHelCM7;
    TH1D* Phip685NegHelCM7;
    TH1D* Phip265NegHelCM8;
    TH1D* Phip335NegHelCM8;
    TH1D* Phip405NegHelCM8;
    TH1D* Phip475NegHelCM8;
    TH1D* Phip545NegHelCM8;
    TH1D* Phip615NegHelCM8;
    TH1D* Phip685NegHelCM8;

    TH1D* Phip265PosHelCM1;
    TH1D* Phip335PosHelCM1;
    TH1D* Phip405PosHelCM1;
    TH1D* Phip475PosHelCM1;
    TH1D* Phip545PosHelCM1;
    TH1D* Phip615PosHelCM1;
    TH1D* Phip685PosHelCM1;
    TH1D* Phip265PosHelCM2;
    TH1D* Phip335PosHelCM2;
    TH1D* Phip405PosHelCM2;
    TH1D* Phip475PosHelCM2;
    TH1D* Phip545PosHelCM2;
    TH1D* Phip615PosHelCM2;
    TH1D* Phip685PosHelCM2;
    TH1D* Phip265PosHelCM3;
    TH1D* Phip335PosHelCM3;
    TH1D* Phip405PosHelCM3;
    TH1D* Phip475PosHelCM3;
    TH1D* Phip545PosHelCM3;
    TH1D* Phip615PosHelCM3;
    TH1D* Phip685PosHelCM3;
    TH1D* Phip265PosHelCM4;
    TH1D* Phip335PosHelCM4;
    TH1D* Phip405PosHelCM4;
    TH1D* Phip475PosHelCM4;
    TH1D* Phip545PosHelCM4;
    TH1D* Phip615PosHelCM4;
    TH1D* Phip685PosHelCM4;
    TH1D* Phip265PosHelCM5;
    TH1D* Phip335PosHelCM5;
    TH1D* Phip405PosHelCM5;
    TH1D* Phip475PosHelCM5;
    TH1D* Phip545PosHelCM5;
    TH1D* Phip615PosHelCM5;
    TH1D* Phip685PosHelCM5;
    TH1D* Phip265PosHelCM6;
    TH1D* Phip335PosHelCM6;
    TH1D* Phip405PosHelCM6;
    TH1D* Phip475PosHelCM6;
    TH1D* Phip545PosHelCM6;
    TH1D* Phip615PosHelCM6;
    TH1D* Phip685PosHelCM6;
    TH1D* Phip265PosHelCM7;
    TH1D* Phip335PosHelCM7;
    TH1D* Phip405PosHelCM7;
    TH1D* Phip475PosHelCM7;
    TH1D* Phip545PosHelCM7;
    TH1D* Phip615PosHelCM7;
    TH1D* Phip685PosHelCM7;
    TH1D* Phip265PosHelCM8;
    TH1D* Phip335PosHelCM8;
    TH1D* Phip405PosHelCM8;
    TH1D* Phip475PosHelCM8;
    TH1D* Phip545PosHelCM8;
    TH1D* Phip615PosHelCM8;
    TH1D* Phip685PosHelCM8;

    TH1D* Phip265NegHelCM1Prompt;
    TH1D* Phip335NegHelCM1Prompt;
    TH1D* Phip405NegHelCM1Prompt;
    TH1D* Phip475NegHelCM1Prompt;
    TH1D* Phip545NegHelCM1Prompt;
    TH1D* Phip615NegHelCM1Prompt;
    TH1D* Phip685NegHelCM1Prompt;
    TH1D* Phip265NegHelCM2Prompt;
    TH1D* Phip335NegHelCM2Prompt;
    TH1D* Phip405NegHelCM2Prompt;
    TH1D* Phip475NegHelCM2Prompt;
    TH1D* Phip545NegHelCM2Prompt;
    TH1D* Phip615NegHelCM2Prompt;
    TH1D* Phip685NegHelCM2Prompt;
    TH1D* Phip265NegHelCM3Prompt;
    TH1D* Phip335NegHelCM3Prompt;
    TH1D* Phip405NegHelCM3Prompt;
    TH1D* Phip475NegHelCM3Prompt;
    TH1D* Phip545NegHelCM3Prompt;
    TH1D* Phip615NegHelCM3Prompt;
    TH1D* Phip685NegHelCM3Prompt;
    TH1D* Phip265NegHelCM4Prompt;
    TH1D* Phip335NegHelCM4Prompt;
    TH1D* Phip405NegHelCM4Prompt;
    TH1D* Phip475NegHelCM4Prompt;
    TH1D* Phip545NegHelCM4Prompt;
    TH1D* Phip615NegHelCM4Prompt;
    TH1D* Phip685NegHelCM4Prompt;
    TH1D* Phip265NegHelCM5Prompt;
    TH1D* Phip335NegHelCM5Prompt;
    TH1D* Phip405NegHelCM5Prompt;
    TH1D* Phip475NegHelCM5Prompt;
    TH1D* Phip545NegHelCM5Prompt;
    TH1D* Phip615NegHelCM5Prompt;
    TH1D* Phip685NegHelCM5Prompt;
    TH1D* Phip265NegHelCM6Prompt;
    TH1D* Phip335NegHelCM6Prompt;
    TH1D* Phip405NegHelCM6Prompt;
    TH1D* Phip475NegHelCM6Prompt;
    TH1D* Phip545NegHelCM6Prompt;
    TH1D* Phip615NegHelCM6Prompt;
    TH1D* Phip685NegHelCM6Prompt;
    TH1D* Phip265NegHelCM7Prompt;
    TH1D* Phip335NegHelCM7Prompt;
    TH1D* Phip405NegHelCM7Prompt;
    TH1D* Phip475NegHelCM7Prompt;
    TH1D* Phip545NegHelCM7Prompt;
    TH1D* Phip615NegHelCM7Prompt;
    TH1D* Phip685NegHelCM7Prompt;
    TH1D* Phip265NegHelCM8Prompt;
    TH1D* Phip335NegHelCM8Prompt;
    TH1D* Phip405NegHelCM8Prompt;
    TH1D* Phip475NegHelCM8Prompt;
    TH1D* Phip545NegHelCM8Prompt;
    TH1D* Phip615NegHelCM8Prompt;
    TH1D* Phip685NegHelCM8Prompt;

    TH1D* Phip265PosHelCM1Prompt;
    TH1D* Phip335PosHelCM1Prompt;
    TH1D* Phip405PosHelCM1Prompt;
    TH1D* Phip475PosHelCM1Prompt;
    TH1D* Phip545PosHelCM1Prompt;
    TH1D* Phip615PosHelCM1Prompt;
    TH1D* Phip685PosHelCM1Prompt;
    TH1D* Phip265PosHelCM2Prompt;
    TH1D* Phip335PosHelCM2Prompt;
    TH1D* Phip405PosHelCM2Prompt;
    TH1D* Phip475PosHelCM2Prompt;
    TH1D* Phip545PosHelCM2Prompt;
    TH1D* Phip615PosHelCM2Prompt;
    TH1D* Phip685PosHelCM2Prompt;
    TH1D* Phip265PosHelCM3Prompt;
    TH1D* Phip335PosHelCM3Prompt;
    TH1D* Phip405PosHelCM3Prompt;
    TH1D* Phip475PosHelCM3Prompt;
    TH1D* Phip545PosHelCM3Prompt;
    TH1D* Phip615PosHelCM3Prompt;
    TH1D* Phip685PosHelCM3Prompt;
    TH1D* Phip265PosHelCM4Prompt;
    TH1D* Phip335PosHelCM4Prompt;
    TH1D* Phip405PosHelCM4Prompt;
    TH1D* Phip475PosHelCM4Prompt;
    TH1D* Phip545PosHelCM4Prompt;
    TH1D* Phip615PosHelCM4Prompt;
    TH1D* Phip685PosHelCM4Prompt;
    TH1D* Phip265PosHelCM5Prompt;
    TH1D* Phip335PosHelCM5Prompt;
    TH1D* Phip405PosHelCM5Prompt;
    TH1D* Phip475PosHelCM5Prompt;
    TH1D* Phip545PosHelCM5Prompt;
    TH1D* Phip615PosHelCM5Prompt;
    TH1D* Phip685PosHelCM5Prompt;
    TH1D* Phip265PosHelCM6Prompt;
    TH1D* Phip335PosHelCM6Prompt;
    TH1D* Phip405PosHelCM6Prompt;
    TH1D* Phip475PosHelCM6Prompt;
    TH1D* Phip545PosHelCM6Prompt;
    TH1D* Phip615PosHelCM6Prompt;
    TH1D* Phip685PosHelCM6Prompt;
    TH1D* Phip265PosHelCM7Prompt;
    TH1D* Phip335PosHelCM7Prompt;
    TH1D* Phip405PosHelCM7Prompt;
    TH1D* Phip475PosHelCM7Prompt;
    TH1D* Phip545PosHelCM7Prompt;
    TH1D* Phip615PosHelCM7Prompt;
    TH1D* Phip685PosHelCM7Prompt;
    TH1D* Phip265PosHelCM8Prompt;
    TH1D* Phip335PosHelCM8Prompt;
    TH1D* Phip405PosHelCM8Prompt;
    TH1D* Phip475PosHelCM8Prompt;
    TH1D* Phip545PosHelCM8Prompt;
    TH1D* Phip615PosHelCM8Prompt;
    TH1D* Phip685PosHelCM8Prompt;

    TH1D* Phip265NegHelCM1Random;
    TH1D* Phip335NegHelCM1Random;
    TH1D* Phip405NegHelCM1Random;
    TH1D* Phip475NegHelCM1Random;
    TH1D* Phip545NegHelCM1Random;
    TH1D* Phip615NegHelCM1Random;
    TH1D* Phip685NegHelCM1Random;
    TH1D* Phip265NegHelCM2Random;
    TH1D* Phip335NegHelCM2Random;
    TH1D* Phip405NegHelCM2Random;
    TH1D* Phip475NegHelCM2Random;
    TH1D* Phip545NegHelCM2Random;
    TH1D* Phip615NegHelCM2Random;
    TH1D* Phip685NegHelCM2Random;
    TH1D* Phip265NegHelCM3Random;
    TH1D* Phip335NegHelCM3Random;
    TH1D* Phip405NegHelCM3Random;
    TH1D* Phip475NegHelCM3Random;
    TH1D* Phip545NegHelCM3Random;
    TH1D* Phip615NegHelCM3Random;
    TH1D* Phip685NegHelCM3Random;
    TH1D* Phip265NegHelCM4Random;
    TH1D* Phip335NegHelCM4Random;
    TH1D* Phip405NegHelCM4Random;
    TH1D* Phip475NegHelCM4Random;
    TH1D* Phip545NegHelCM4Random;
    TH1D* Phip615NegHelCM4Random;
    TH1D* Phip685NegHelCM4Random;
    TH1D* Phip265NegHelCM5Random;
    TH1D* Phip335NegHelCM5Random;
    TH1D* Phip405NegHelCM5Random;
    TH1D* Phip475NegHelCM5Random;
    TH1D* Phip545NegHelCM5Random;
    TH1D* Phip615NegHelCM5Random;
    TH1D* Phip685NegHelCM5Random;
    TH1D* Phip265NegHelCM6Random;
    TH1D* Phip335NegHelCM6Random;
    TH1D* Phip405NegHelCM6Random;
    TH1D* Phip475NegHelCM6Random;
    TH1D* Phip545NegHelCM6Random;
    TH1D* Phip615NegHelCM6Random;
    TH1D* Phip685NegHelCM6Random;
    TH1D* Phip265NegHelCM7Random;
    TH1D* Phip335NegHelCM7Random;
    TH1D* Phip405NegHelCM7Random;
    TH1D* Phip475NegHelCM7Random;
    TH1D* Phip545NegHelCM7Random;
    TH1D* Phip615NegHelCM7Random;
    TH1D* Phip685NegHelCM7Random;
    TH1D* Phip265NegHelCM8Random;
    TH1D* Phip335NegHelCM8Random;
    TH1D* Phip405NegHelCM8Random;
    TH1D* Phip475NegHelCM8Random;
    TH1D* Phip545NegHelCM8Random;
    TH1D* Phip615NegHelCM8Random;
    TH1D* Phip685NegHelCM8Random;

    TH1D* Phip265PosHelCM1Random;
    TH1D* Phip335PosHelCM1Random;
    TH1D* Phip405PosHelCM1Random;
    TH1D* Phip475PosHelCM1Random;
    TH1D* Phip545PosHelCM1Random;
    TH1D* Phip615PosHelCM1Random;
    TH1D* Phip685PosHelCM1Random;
    TH1D* Phip265PosHelCM2Random;
    TH1D* Phip335PosHelCM2Random;
    TH1D* Phip405PosHelCM2Random;
    TH1D* Phip475PosHelCM2Random;
    TH1D* Phip545PosHelCM2Random;
    TH1D* Phip615PosHelCM2Random;
    TH1D* Phip685PosHelCM2Random;
    TH1D* Phip265PosHelCM3Random;
    TH1D* Phip335PosHelCM3Random;
    TH1D* Phip405PosHelCM3Random;
    TH1D* Phip475PosHelCM3Random;
    TH1D* Phip545PosHelCM3Random;
    TH1D* Phip615PosHelCM3Random;
    TH1D* Phip685PosHelCM3Random;
    TH1D* Phip265PosHelCM4Random;
    TH1D* Phip335PosHelCM4Random;
    TH1D* Phip405PosHelCM4Random;
    TH1D* Phip475PosHelCM4Random;
    TH1D* Phip545PosHelCM4Random;
    TH1D* Phip615PosHelCM4Random;
    TH1D* Phip685PosHelCM4Random;
    TH1D* Phip265PosHelCM5Random;
    TH1D* Phip335PosHelCM5Random;
    TH1D* Phip405PosHelCM5Random;
    TH1D* Phip475PosHelCM5Random;
    TH1D* Phip545PosHelCM5Random;
    TH1D* Phip615PosHelCM5Random;
    TH1D* Phip685PosHelCM5Random;
    TH1D* Phip265PosHelCM6Random;
    TH1D* Phip335PosHelCM6Random;
    TH1D* Phip405PosHelCM6Random;
    TH1D* Phip475PosHelCM6Random;
    TH1D* Phip545PosHelCM6Random;
    TH1D* Phip615PosHelCM6Random;
    TH1D* Phip685PosHelCM6Random;
    TH1D* Phip265PosHelCM7Random;
    TH1D* Phip335PosHelCM7Random;
    TH1D* Phip405PosHelCM7Random;
    TH1D* Phip475PosHelCM7Random;
    TH1D* Phip545PosHelCM7Random;
    TH1D* Phip615PosHelCM7Random;
    TH1D* Phip685PosHelCM7Random;
    TH1D* Phip265PosHelCM8Random;
    TH1D* Phip335PosHelCM8Random;
    TH1D* Phip405PosHelCM8Random;
    TH1D* Phip475PosHelCM8Random;
    TH1D* Phip545PosHelCM8Random;
    TH1D* Phip615PosHelCM8Random;
    TH1D* Phip685PosHelCM8Random;

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
