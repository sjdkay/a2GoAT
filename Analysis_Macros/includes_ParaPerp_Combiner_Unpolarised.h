#include "TROOT.h"
#include <TChain.h>
#include "TFile.h"

#include "TObject.h"
#include "TApplication.h"
#include "TVector3.h"
#include "TVector2.h"
#include "TString.h"
#include "TChain.h"
#include "TGraphErrors.h"
#include "TLorentzVector.h"

#include "TStyle.h"
#include "TColor.h"
#include "TLegend.h"

#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TMath.h"

#include "TFile.h"
#include "TTree.h"

#include <stdio.h>
#include <iostream>
#include <unistd.h>

#include <fstream>

#include <string.h>

using namespace std;
using namespace TMath;

TH1D* time_Para;
TH1D* time_cut_Para;
TH1D* Eg_Para;
TH1D* WCPhiDifference_Para;
TH1D* EpKin_Para;
TH1D* EpCorrected_Para;
TH1D* OAngle_Para;
TH1D* WCZnRecon_Para;
TH1D* Theta_Scattered_Para;
TH1D* Phi_Scattered_Para;
TH1D* EpKinEpCorrDiff_Para;
TH1D* EpEpCorrDiff_Para;
TH1D* MMpEpCorrected_Para;
TH1D* ZpDist_Para;
TH1D* ZpPhiScatNeg180_Para;
TH1D* ZpPhiScat0_Para;
TH1D* ZpPhiScatPos180_Para;
TH1D* MMp200300_Para;
TH1D* MMp300400_Para;
TH1D* MMp400500_Para;
TH1D* MMp500600_Para;
TH1D* MMp600700_Para;
TH1D* MMp700800_Para;
TH1D* MMp800900_Para;
TH1D* Phi_Scattered_410MeV_Para;
TH1D* Phi_Scattered_430MeV_Para;
TH1D* Phi_Scattered_450MeV_Para;
TH1D* Phi_Scattered_470MeV_Para;
TH1D* Phi_Scattered_490MeV_Para;
TH1D* Phi_Scattered_510MeV_Para;
TH1D* Phi_Scattered_530MeV_Para;
TH1D* Phi_Scattered_550MeV_Para;
TH1D* Phi_Scattered_570MeV_Para;
TH1D* Phi_Scattered_590MeV_Para;
TH1D* Phi_Scattered_610MeV_Para;
TH1D* Phi_Scattered_630MeV_Para;

TH1D* ThetanDist_Para;
TH1D* ThetanRecDist_Para;
TH1D* ThetanDiffDist_Para;
TH2D* ThetanDiffZp_Para;
TH1D* ThetanCorrDist_Para;
TH1D* ThetanCorrDiffDist_Para;
TH1D* ThetanCorrRecDiffDist_Para;
TH2D* ThetanCorrDiffZp_Para;

TH1D* ThetaRecPiDiff_Para;
TH2D* ThetanThetaRecPi_Para;
TH2D* ThetanThetaRecPiDiff_Para;
TH1D* ThetaRecPDiff_Para;
TH2D* ThetanThetaRecP_Para;
TH2D* ThetanThetaRecPDiff_Para;
TH2D* E_dE_Para;
TH2D* KinEp_dE_Para;
TH2D* ThetaScPhiSc_Para;
TH2D* E_KinEp_Para;
TH2D* PhinDiffWCZRec_Para;

TH1D* time_Perp;
TH1D* time_cut_Perp;
TH1D* Eg_Perp;
TH1D* WCPhiDifference_Perp;
TH1D* EpKin_Perp;
TH1D* EpCorrected_Perp;
TH1D* OAngle_Perp;
TH1D* WCZnRecon_Perp;
TH1D* Theta_Scattered_Perp;
TH1D* Phi_Scattered_Perp;
TH1D* EpKinEpCorrDiff_Perp;
TH1D* EpEpCorrDiff_Perp;
TH1D* MMpEpCorrected_Perp;
TH1D* ZpDist_Perp;
TH1D* ZpPhiScatNeg180_Perp;
TH1D* ZpPhiScat0_Perp;
TH1D* ZpPhiScatPos180_Perp;
TH1D* MMp200300_Perp;
TH1D* MMp300400_Perp;
TH1D* MMp400500_Perp;
TH1D* MMp500600_Perp;
TH1D* MMp600700_Perp;
TH1D* MMp700800_Perp;
TH1D* MMp800900_Perp;
TH1D* Phi_Scattered_410MeV_Perp;
TH1D* Phi_Scattered_430MeV_Perp;
TH1D* Phi_Scattered_450MeV_Perp;
TH1D* Phi_Scattered_470MeV_Perp;
TH1D* Phi_Scattered_490MeV_Perp;
TH1D* Phi_Scattered_510MeV_Perp;
TH1D* Phi_Scattered_530MeV_Perp;
TH1D* Phi_Scattered_550MeV_Perp;
TH1D* Phi_Scattered_570MeV_Perp;
TH1D* Phi_Scattered_590MeV_Perp;
TH1D* Phi_Scattered_610MeV_Perp;
TH1D* Phi_Scattered_630MeV_Perp;

TH1D* ThetanDist_Perp;
TH1D* ThetanRecDist_Perp;
TH1D* ThetanDiffDist_Perp;
TH2D* ThetanDiffZp_Perp;
TH1D* ThetanCorrDist_Perp;
TH1D* ThetanCorrDiffDist_Perp;
TH1D* ThetanCorrRecDiffDist_Perp;
TH2D* ThetanCorrDiffZp_Perp;

TH1D* ThetaRecPiDiff_Perp;
TH2D* ThetanThetaRecPi_Perp;
TH2D* ThetanThetaRecPiDiff_Perp;
TH1D* ThetaRecPDiff_Perp;
TH2D* ThetanThetaRecP_Perp;
TH2D* ThetanThetaRecPDiff_Perp;
TH2D* E_dE_Perp;
TH2D* KinEp_dE_Perp;
TH2D* ThetaScPhiSc_Perp;
TH2D* E_KinEp_Perp;
TH2D* PhinDiffWCZRec_Perp;

TH1D* time;
TH1D* time_cut;
TH1D* Eg;
TH1D* WCPhiDifference;
TH1D* EpKin;
TH1D* EpCorrected;
TH1D* OAngle;
TH1D* WCZnRecon;
TH1D* ThetaSc;
TH1D* PhiSc;
TH1D* EpKinEpCorrDiff;
TH1D* EpEpCorrDiff;
TH1D* MMpEpCorrected;
TH1D* ZpDist;
TH1D* ZpPhiScatNeg180;
TH1D* ZpPhiScat0;
TH1D* ZpPhiScatPos180;
TH1D* MMp200300;
TH1D* MMp300400;
TH1D* MMp400500;
TH1D* MMp500600;
TH1D* MMp600700;
TH1D* MMp700800;
TH1D* MMp800900;
TH1D* PhiSc410;
TH1D* PhiSc430;
TH1D* PhiSc450;
TH1D* PhiSc470;
TH1D* PhiSc490;
TH1D* PhiSc510;
TH1D* PhiSc530;
TH1D* PhiSc550;
TH1D* PhiSc570;
TH1D* PhiSc590;
TH1D* PhiSc610;
TH1D* PhiSc630;

TH1D* ThetanDist;
TH1D* ThetanRecDist;
TH1D* ThetanDiffDist;
TH2D* ThetanDiffZp;
TH1D* ThetanCorrDist;
TH1D* ThetanCorrDiffDist;
TH1D* ThetanCorrRecDiffDist;
TH2D* ThetanCorrDiffZp;

TH1D* ThetaRecPiDiff;
TH2D* ThetanThetaRecPi;
TH2D* ThetanThetaRecPiDiff;
TH1D* ThetaRecPDiff;
TH2D* ThetanThetaRecP;
TH2D* ThetanThetaRecPDiff;
TH2D* E_dE;
TH2D* KinEp_dE;
TH2D* ThetaScPhiSc;
TH2D* E_KinEp;
TH2D* PhinDiffWCZRec;

Double_t NPara;
Double_t NPerp;
Double_t ScaleFactor;
Double_t ScaleFactorErr;

