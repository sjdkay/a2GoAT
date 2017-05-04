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

TH1D* time_Unpol;
TH1D* time_cut_Unpol;
TH1D* Eg_Unpol;
TH1D* WCPhiDifference_Unpol;
TH1D* EpKin_Unpol;
TH1D* EpCorrected_Unpol;
TH1D* OAngle_Unpol;
TH1D* WCZnRecon_Unpol;
TH1D* Theta_Scattered_Unpol;
TH1D* Phi_Scattered_Unpol;
TH1D* EpKinEpCorrDiff_Unpol;
TH1D* EpEpCorrDiff_Unpol;
TH1D* MMpEpCorrected_Unpol;
TH1D* ZpDist_Unpol;
TH1D* ZpPhiScatNeg180_Unpol;
TH1D* ZpPhiScat0_Unpol;
TH1D* ZpPhiScatPos180_Unpol;
TH1D* MMp200300_Unpol;
TH1D* MMp300400_Unpol;
TH1D* MMp400500_Unpol;
TH1D* MMp500600_Unpol;
TH1D* MMp600700_Unpol;
TH1D* MMp700800_Unpol;
TH1D* MMp800900_Unpol;
TH1D* Phi_Scattered_410MeV_Unpol;
TH1D* Phi_Scattered_430MeV_Unpol;
TH1D* Phi_Scattered_450MeV_Unpol;
TH1D* Phi_Scattered_470MeV_Unpol;
TH1D* Phi_Scattered_490MeV_Unpol;
TH1D* Phi_Scattered_510MeV_Unpol;
TH1D* Phi_Scattered_530MeV_Unpol;
TH1D* Phi_Scattered_550MeV_Unpol;
TH1D* Phi_Scattered_570MeV_Unpol;
TH1D* Phi_Scattered_590MeV_Unpol;
TH1D* Phi_Scattered_610MeV_Unpol;
TH1D* Phi_Scattered_630MeV_Unpol;
TH1D* ThetaRecPiDiff_Unpol;
TH2D* ThetanThetaRecPi_Unpol;
TH2D* ThetanThetaRecPiDiff_Unpol;
TH1D* ThetaRecPDiff_Unpol;
TH2D* ThetanThetaRecP_Unpol;
TH2D* ThetanThetaRecPDiff_Unpol;
TH2D* E_dE_Unpol;
TH2D* KinEp_dE_Unpol;
TH2D* ThetaScPhiSc_Unpol;
TH2D* E_KinEp_Unpol;
TH2D* PhinDiffWCZRec_Unpol;

Double_t NPara;
Double_t NPerp;
Double_t ScaleFactor;
Double_t ScaleFactorErr;

