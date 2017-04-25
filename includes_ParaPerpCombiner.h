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
TH1D* MMpEpCorrectedCut_Para;
TH1D* OAngleCut_Para;
TH1D* EgCut_Para;
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
TH1D* Phip_410MeVCM1_Para;
TH1D* Phip_430MeVCM1_Para;
TH1D* Phip_450MeVCM1_Para;
TH1D* Phip_470MeVCM1_Para;  
TH1D* Phip_490MeVCM1_Para;  
TH1D* Phip_510MeVCM1_Para;  
TH1D* Phip_530MeVCM1_Para;  
TH1D* Phip_550MeVCM1_Para;  
TH1D* Phip_570MeVCM1_Para;  
TH1D* Phip_590MeVCM1_Para;  
TH1D* Phip_610MeVCM1_Para;  
TH1D* Phip_630MeVCM1_Para;  
TH1D* Phin_410MeVCM1_Para;
TH1D* Phin_430MeVCM1_Para;
TH1D* Phin_450MeVCM1_Para;
TH1D* Phin_470MeVCM1_Para;  
TH1D* Phin_490MeVCM1_Para;  
TH1D* Phin_510MeVCM1_Para;  
TH1D* Phin_530MeVCM1_Para;  
TH1D* Phin_550MeVCM1_Para;  
TH1D* Phin_570MeVCM1_Para;  
TH1D* Phin_590MeVCM1_Para;  
TH1D* Phin_610MeVCM1_Para;  
TH1D* Phin_630MeVCM1_Para;
TH1D* Phip_410MeVCM2_Para;
TH1D* Phip_430MeVCM2_Para;
TH1D* Phip_450MeVCM2_Para;
TH1D* Phip_470MeVCM2_Para;  
TH1D* Phip_490MeVCM2_Para;  
TH1D* Phip_510MeVCM2_Para;  
TH1D* Phip_530MeVCM2_Para;  
TH1D* Phip_550MeVCM2_Para;  
TH1D* Phip_570MeVCM2_Para;  
TH1D* Phip_590MeVCM2_Para;  
TH1D* Phip_610MeVCM2_Para;  
TH1D* Phip_630MeVCM2_Para;  
TH1D* Phin_410MeVCM2_Para;
TH1D* Phin_430MeVCM2_Para;
TH1D* Phin_450MeVCM2_Para;
TH1D* Phin_470MeVCM2_Para;  
TH1D* Phin_490MeVCM2_Para;  
TH1D* Phin_510MeVCM2_Para;  
TH1D* Phin_530MeVCM2_Para;  
TH1D* Phin_550MeVCM2_Para;  
TH1D* Phin_570MeVCM2_Para;  
TH1D* Phin_590MeVCM2_Para;  
TH1D* Phin_610MeVCM2_Para;  
TH1D* Phin_630MeVCM2_Para;
TH1D* Phip_410MeVCM3_Para;
TH1D* Phip_430MeVCM3_Para;
TH1D* Phip_450MeVCM3_Para;
TH1D* Phip_470MeVCM3_Para;  
TH1D* Phip_490MeVCM3_Para;  
TH1D* Phip_510MeVCM3_Para;  
TH1D* Phip_530MeVCM3_Para;  
TH1D* Phip_550MeVCM3_Para;  
TH1D* Phip_570MeVCM3_Para;  
TH1D* Phip_590MeVCM3_Para;  
TH1D* Phip_610MeVCM3_Para;  
TH1D* Phip_630MeVCM3_Para;  
TH1D* Phin_410MeVCM3_Para;
TH1D* Phin_430MeVCM3_Para;
TH1D* Phin_450MeVCM3_Para;
TH1D* Phin_470MeVCM3_Para;  
TH1D* Phin_490MeVCM3_Para;  
TH1D* Phin_510MeVCM3_Para;  
TH1D* Phin_530MeVCM3_Para;  
TH1D* Phin_550MeVCM3_Para;  
TH1D* Phin_570MeVCM3_Para;  
TH1D* Phin_590MeVCM3_Para;  
TH1D* Phin_610MeVCM3_Para;  
TH1D* Phin_630MeVCM3_Para;
TH1D* Phip_410MeVCM4_Para;
TH1D* Phip_430MeVCM4_Para;
TH1D* Phip_450MeVCM4_Para;
TH1D* Phip_470MeVCM4_Para;  
TH1D* Phip_490MeVCM4_Para;  
TH1D* Phip_510MeVCM4_Para;  
TH1D* Phip_530MeVCM4_Para;  
TH1D* Phip_550MeVCM4_Para;  
TH1D* Phip_570MeVCM4_Para;  
TH1D* Phip_590MeVCM4_Para;  
TH1D* Phip_610MeVCM4_Para;  
TH1D* Phip_630MeVCM4_Para;  
TH1D* Phin_410MeVCM4_Para;
TH1D* Phin_430MeVCM4_Para;
TH1D* Phin_450MeVCM4_Para;
TH1D* Phin_470MeVCM4_Para;  
TH1D* Phin_490MeVCM4_Para;  
TH1D* Phin_510MeVCM4_Para;  
TH1D* Phin_530MeVCM4_Para;  
TH1D* Phin_550MeVCM4_Para;  
TH1D* Phin_570MeVCM4_Para;  
TH1D* Phin_590MeVCM4_Para;  
TH1D* Phin_610MeVCM4_Para;  
TH1D* Phin_630MeVCM4_Para;
TH1D* Phip_410MeVCM5_Para;
TH1D* Phip_430MeVCM5_Para;
TH1D* Phip_450MeVCM5_Para;
TH1D* Phip_470MeVCM5_Para;  
TH1D* Phip_490MeVCM5_Para;  
TH1D* Phip_510MeVCM5_Para;  
TH1D* Phip_530MeVCM5_Para;  
TH1D* Phip_550MeVCM5_Para;  
TH1D* Phip_570MeVCM5_Para;  
TH1D* Phip_590MeVCM5_Para;  
TH1D* Phip_610MeVCM5_Para;  
TH1D* Phip_630MeVCM5_Para;  
TH1D* Phin_410MeVCM5_Para;
TH1D* Phin_430MeVCM5_Para;
TH1D* Phin_450MeVCM5_Para;
TH1D* Phin_470MeVCM5_Para;  
TH1D* Phin_490MeVCM5_Para;  
TH1D* Phin_510MeVCM5_Para;  
TH1D* Phin_530MeVCM5_Para;  
TH1D* Phin_550MeVCM5_Para;  
TH1D* Phin_570MeVCM5_Para;  
TH1D* Phin_590MeVCM5_Para;  
TH1D* Phin_610MeVCM5_Para;  
TH1D* Phin_630MeVCM5_Para;
TH1D* Phip_410MeVCM6_Para;
TH1D* Phip_430MeVCM6_Para;
TH1D* Phip_450MeVCM6_Para;
TH1D* Phip_470MeVCM6_Para;  
TH1D* Phip_490MeVCM6_Para;  
TH1D* Phip_510MeVCM6_Para;  
TH1D* Phip_530MeVCM6_Para;  
TH1D* Phip_550MeVCM6_Para;  
TH1D* Phip_570MeVCM6_Para;  
TH1D* Phip_590MeVCM6_Para;  
TH1D* Phip_610MeVCM6_Para;  
TH1D* Phip_630MeVCM6_Para;  
TH1D* Phin_410MeVCM6_Para;
TH1D* Phin_430MeVCM6_Para;
TH1D* Phin_450MeVCM6_Para;
TH1D* Phin_470MeVCM6_Para;  
TH1D* Phin_490MeVCM6_Para;  
TH1D* Phin_510MeVCM6_Para;  
TH1D* Phin_530MeVCM6_Para;  
TH1D* Phin_550MeVCM6_Para;  
TH1D* Phin_570MeVCM6_Para;  
TH1D* Phin_590MeVCM6_Para;  
TH1D* Phin_610MeVCM6_Para;  
TH1D* Phin_630MeVCM6_Para;
TH1D* ThetaRecPiDiff_Para;
TH2D* ThetanThetaRecPi_Para;
TH2D* ThetanRetaRecPiDiff_Para;
TH1D* ThetaRecPDiff_Para;
TH2D* ThetanThetaRecP_Para;
TH2D* ThetanRetaRecPDiff_Para;
TH2D* E_dE_Para;
TH2D* E_dE_KinCut_Para;
TH2D* KinEp_dE_Para;
TH2D* ThetaScPhiSc_Para;
TH2D* E_KinEP_Para;
TH2D* PhinDiffWCZRec_Para;
TH2D PhinDiffWCZRec_KinCut_Para;

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
TH1D* MMpEpCorrectedCut_Perp;
TH1D* OAngleCut_Perp;
TH1D* EgCut_Perp;
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
TH1D* Phip_410MeVCM1_Perp;
TH1D* Phip_430MeVCM1_Perp;
TH1D* Phip_450MeVCM1_Perp;
TH1D* Phip_470MeVCM1_Perp;  
TH1D* Phip_490MeVCM1_Perp;  
TH1D* Phip_510MeVCM1_Perp;  
TH1D* Phip_530MeVCM1_Perp;  
TH1D* Phip_550MeVCM1_Perp;  
TH1D* Phip_570MeVCM1_Perp;  
TH1D* Phip_590MeVCM1_Perp;  
TH1D* Phip_610MeVCM1_Perp;  
TH1D* Phip_630MeVCM1_Perp;  
TH1D* Phin_410MeVCM1_Perp;
TH1D* Phin_430MeVCM1_Perp;
TH1D* Phin_450MeVCM1_Perp;
TH1D* Phin_470MeVCM1_Perp;  
TH1D* Phin_490MeVCM1_Perp;  
TH1D* Phin_510MeVCM1_Perp;  
TH1D* Phin_530MeVCM1_Perp;  
TH1D* Phin_550MeVCM1_Perp;  
TH1D* Phin_570MeVCM1_Perp;  
TH1D* Phin_590MeVCM1_Perp;  
TH1D* Phin_610MeVCM1_Perp;  
TH1D* Phin_630MeVCM1_Perp;
TH1D* Phip_410MeVCM2_Perp;
TH1D* Phip_430MeVCM2_Perp;
TH1D* Phip_450MeVCM2_Perp;
TH1D* Phip_470MeVCM2_Perp;  
TH1D* Phip_490MeVCM2_Perp;  
TH1D* Phip_510MeVCM2_Perp;  
TH1D* Phip_530MeVCM2_Perp;  
TH1D* Phip_550MeVCM2_Perp;  
TH1D* Phip_570MeVCM2_Perp;  
TH1D* Phip_590MeVCM2_Perp;  
TH1D* Phip_610MeVCM2_Perp;  
TH1D* Phip_630MeVCM2_Perp;  
TH1D* Phin_410MeVCM2_Perp;
TH1D* Phin_430MeVCM2_Perp;
TH1D* Phin_450MeVCM2_Perp;
TH1D* Phin_470MeVCM2_Perp;  
TH1D* Phin_490MeVCM2_Perp;  
TH1D* Phin_510MeVCM2_Perp;  
TH1D* Phin_530MeVCM2_Perp;  
TH1D* Phin_550MeVCM2_Perp;  
TH1D* Phin_570MeVCM2_Perp;  
TH1D* Phin_590MeVCM2_Perp;  
TH1D* Phin_610MeVCM2_Perp;  
TH1D* Phin_630MeVCM2_Perp;
TH1D* Phip_410MeVCM3_Perp;
TH1D* Phip_430MeVCM3_Perp;
TH1D* Phip_450MeVCM3_Perp;
TH1D* Phip_470MeVCM3_Perp;  
TH1D* Phip_490MeVCM3_Perp;  
TH1D* Phip_510MeVCM3_Perp;  
TH1D* Phip_530MeVCM3_Perp;  
TH1D* Phip_550MeVCM3_Perp;  
TH1D* Phip_570MeVCM3_Perp;  
TH1D* Phip_590MeVCM3_Perp;  
TH1D* Phip_610MeVCM3_Perp;  
TH1D* Phip_630MeVCM3_Perp;  
TH1D* Phin_410MeVCM3_Perp;
TH1D* Phin_430MeVCM3_Perp;
TH1D* Phin_450MeVCM3_Perp;
TH1D* Phin_470MeVCM3_Perp;  
TH1D* Phin_490MeVCM3_Perp;  
TH1D* Phin_510MeVCM3_Perp;  
TH1D* Phin_530MeVCM3_Perp;  
TH1D* Phin_550MeVCM3_Perp;  
TH1D* Phin_570MeVCM3_Perp;  
TH1D* Phin_590MeVCM3_Perp;  
TH1D* Phin_610MeVCM3_Perp;  
TH1D* Phin_630MeVCM3_Perp;
TH1D* Phip_410MeVCM4_Perp;
TH1D* Phip_430MeVCM4_Perp;
TH1D* Phip_450MeVCM4_Perp;
TH1D* Phip_470MeVCM4_Perp;  
TH1D* Phip_490MeVCM4_Perp;  
TH1D* Phip_510MeVCM4_Perp;  
TH1D* Phip_530MeVCM4_Perp;  
TH1D* Phip_550MeVCM4_Perp;  
TH1D* Phip_570MeVCM4_Perp;  
TH1D* Phip_590MeVCM4_Perp;  
TH1D* Phip_610MeVCM4_Perp;  
TH1D* Phip_630MeVCM4_Perp;  
TH1D* Phin_410MeVCM4_Perp;
TH1D* Phin_430MeVCM4_Perp;
TH1D* Phin_450MeVCM4_Perp;
TH1D* Phin_470MeVCM4_Perp;  
TH1D* Phin_490MeVCM4_Perp;  
TH1D* Phin_510MeVCM4_Perp;  
TH1D* Phin_530MeVCM4_Perp;  
TH1D* Phin_550MeVCM4_Perp;  
TH1D* Phin_570MeVCM4_Perp;  
TH1D* Phin_590MeVCM4_Perp;  
TH1D* Phin_610MeVCM4_Perp;  
TH1D* Phin_630MeVCM4_Perp;
TH1D* Phip_410MeVCM5_Perp;
TH1D* Phip_430MeVCM5_Perp;
TH1D* Phip_450MeVCM5_Perp;
TH1D* Phip_470MeVCM5_Perp;  
TH1D* Phip_490MeVCM5_Perp;  
TH1D* Phip_510MeVCM5_Perp;  
TH1D* Phip_530MeVCM5_Perp;  
TH1D* Phip_550MeVCM5_Perp;  
TH1D* Phip_570MeVCM5_Perp;  
TH1D* Phip_590MeVCM5_Perp;  
TH1D* Phip_610MeVCM5_Perp;  
TH1D* Phip_630MeVCM5_Perp;  
TH1D* Phin_410MeVCM5_Perp;
TH1D* Phin_430MeVCM5_Perp;
TH1D* Phin_450MeVCM5_Perp;
TH1D* Phin_470MeVCM5_Perp;  
TH1D* Phin_490MeVCM5_Perp;  
TH1D* Phin_510MeVCM5_Perp;  
TH1D* Phin_530MeVCM5_Perp;  
TH1D* Phin_550MeVCM5_Perp;  
TH1D* Phin_570MeVCM5_Perp;  
TH1D* Phin_590MeVCM5_Perp;  
TH1D* Phin_610MeVCM5_Perp;  
TH1D* Phin_630MeVCM5_Perp;
TH1D* Phip_410MeVCM6_Perp;
TH1D* Phip_430MeVCM6_Perp;
TH1D* Phip_450MeVCM6_Perp;
TH1D* Phip_470MeVCM6_Perp;  
TH1D* Phip_490MeVCM6_Perp;  
TH1D* Phip_510MeVCM6_Perp;  
TH1D* Phip_530MeVCM6_Perp;  
TH1D* Phip_550MeVCM6_Perp;  
TH1D* Phip_570MeVCM6_Perp;  
TH1D* Phip_590MeVCM6_Perp;  
TH1D* Phip_610MeVCM6_Perp;  
TH1D* Phip_630MeVCM6_Perp;  
TH1D* Phin_410MeVCM6_Perp;
TH1D* Phin_430MeVCM6_Perp;
TH1D* Phin_450MeVCM6_Perp;
TH1D* Phin_470MeVCM6_Perp;  
TH1D* Phin_490MeVCM6_Perp;  
TH1D* Phin_510MeVCM6_Perp;  
TH1D* Phin_530MeVCM6_Perp;  
TH1D* Phin_550MeVCM6_Perp;  
TH1D* Phin_570MeVCM6_Perp;  
TH1D* Phin_590MeVCM6_Perp;  
TH1D* Phin_610MeVCM6_Perp;  
TH1D* Phin_630MeVCM6_Perp;
TH1D* ThetaRecPiDiff_Perp;
TH2D* ThetanThetaRecPi_Perp;
TH2D* ThetanRetaRecPiDiff_Perp;
TH1D* ThetaRecPDiff_Perp;
TH2D* ThetanThetaRecP_Perp;
TH2D* ThetanRetaRecPDiff_Perp;
TH2D* E_dE_Perp;
TH2D* E_dE_KinCut_Perp;
TH2D* KinEp_dE_Perp;
TH2D* ThetaScPhiSc_Perp;
TH2D* E_KinEP_Perp;
TH2D* PhinDiffWCZRec_Perp;
TH2D PhinDiffWCZRec_KinCut_Perp;
