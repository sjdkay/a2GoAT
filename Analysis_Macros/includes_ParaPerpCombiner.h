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
TH1D* PhiDifference_Para;
TH1D* EpKin_Para;
TH1D* EpCorrected_Para;
TH1D* EpKinEpCorrDiff_Para;
TH1D* EpEpCorrDiff_Para;

TH1D* MMpEpCorrected_Para;
TH1D* OAngle_Para;
TH1D* WCZnRecon_Para;

TH1D* MMp200300_Para;
TH1D* MMp300400_Para;
TH1D* MMp400500_Para;
TH1D* MMp500600_Para;
TH1D* MMp600700_Para;
TH1D* MMp700800_Para;
TH1D* MMp800900_Para;

TH1D* ZpDist_Para;

TH1D* Phip_435MeVCM1_Para;
TH1D* Phip_455MeVCM1_Para;
TH1D* Phip_475MeVCM1_Para;
TH1D* Phip_495MeVCM1_Para;
TH1D* Phip_515MeVCM1_Para;
TH1D* Phip_535MeVCM1_Para;
TH1D* Phip_555MeVCM1_Para;
TH1D* Phip_575MeVCM1_Para;
TH1D* Phip_595MeVCM1_Para;
TH1D* Phip_615MeVCM1_Para;

TH1D* Phip_435MeVCM2_Para;
TH1D* Phip_455MeVCM2_Para;
TH1D* Phip_475MeVCM2_Para;
TH1D* Phip_495MeVCM2_Para;
TH1D* Phip_515MeVCM2_Para;
TH1D* Phip_535MeVCM2_Para;
TH1D* Phip_555MeVCM2_Para;
TH1D* Phip_575MeVCM2_Para;
TH1D* Phip_595MeVCM2_Para;
TH1D* Phip_615MeVCM2_Para;

TH1D* Phip_435MeVCM3_Para;
TH1D* Phip_455MeVCM3_Para;
TH1D* Phip_475MeVCM3_Para;
TH1D* Phip_495MeVCM3_Para;
TH1D* Phip_515MeVCM3_Para;
TH1D* Phip_535MeVCM3_Para;
TH1D* Phip_555MeVCM3_Para;
TH1D* Phip_575MeVCM3_Para;
TH1D* Phip_595MeVCM3_Para;
TH1D* Phip_615MeVCM3_Para;

TH1D* Phip_435MeVCM4_Para;
TH1D* Phip_455MeVCM4_Para;
TH1D* Phip_475MeVCM4_Para;
TH1D* Phip_495MeVCM4_Para;
TH1D* Phip_515MeVCM4_Para;
TH1D* Phip_535MeVCM4_Para;
TH1D* Phip_555MeVCM4_Para;
TH1D* Phip_575MeVCM4_Para;
TH1D* Phip_595MeVCM4_Para;
TH1D* Phip_615MeVCM4_Para;

TH1D* Phip_435MeVCM5_Para;
TH1D* Phip_455MeVCM5_Para;
TH1D* Phip_475MeVCM5_Para;
TH1D* Phip_495MeVCM5_Para;
TH1D* Phip_515MeVCM5_Para;
TH1D* Phip_535MeVCM5_Para;
TH1D* Phip_555MeVCM5_Para;
TH1D* Phip_575MeVCM5_Para;
TH1D* Phip_595MeVCM5_Para;
TH1D* Phip_615MeVCM5_Para;

TH1D* Phip_435MeVCM6_Para;
TH1D* Phip_455MeVCM6_Para;
TH1D* Phip_475MeVCM6_Para;
TH1D* Phip_495MeVCM6_Para;
TH1D* Phip_515MeVCM6_Para;
TH1D* Phip_535MeVCM6_Para;
TH1D* Phip_555MeVCM6_Para;
TH1D* Phip_575MeVCM6_Para;
TH1D* Phip_595MeVCM6_Para;
TH1D* Phip_615MeVCM6_Para;

TH1D* Phip_435MeVCM7_Para;
TH1D* Phip_455MeVCM7_Para;
TH1D* Phip_475MeVCM7_Para;
TH1D* Phip_495MeVCM7_Para;
TH1D* Phip_515MeVCM7_Para;
TH1D* Phip_535MeVCM7_Para;
TH1D* Phip_555MeVCM7_Para;
TH1D* Phip_575MeVCM7_Para;
TH1D* Phip_595MeVCM7_Para;
TH1D* Phip_615MeVCM7_Para;

TH1D* Phip_435MeVCM8_Para;
TH1D* Phip_455MeVCM8_Para;
TH1D* Phip_475MeVCM8_Para;
TH1D* Phip_495MeVCM8_Para;
TH1D* Phip_515MeVCM8_Para;
TH1D* Phip_535MeVCM8_Para;
TH1D* Phip_555MeVCM8_Para;
TH1D* Phip_575MeVCM8_Para;
TH1D* Phip_595MeVCM8_Para;
TH1D* Phip_615MeVCM8_Para;

TH1D* Phip_435MeVCM9_Para;
TH1D* Phip_455MeVCM9_Para;
TH1D* Phip_475MeVCM9_Para;
TH1D* Phip_495MeVCM9_Para;
TH1D* Phip_515MeVCM9_Para;
TH1D* Phip_535MeVCM9_Para;
TH1D* Phip_555MeVCM9_Para;
TH1D* Phip_575MeVCM9_Para;
TH1D* Phip_595MeVCM9_Para;
TH1D* Phip_615MeVCM9_Para;

TH1D* Phip_435MeVCM10_Para;
TH1D* Phip_455MeVCM10_Para;
TH1D* Phip_475MeVCM10_Para;
TH1D* Phip_495MeVCM10_Para;
TH1D* Phip_515MeVCM10_Para;
TH1D* Phip_535MeVCM10_Para;
TH1D* Phip_555MeVCM10_Para;
TH1D* Phip_575MeVCM10_Para;
TH1D* Phip_595MeVCM10_Para;
TH1D* Phip_615MeVCM10_Para;

TH1D* ThetanDist_Para;
TH1D* ThetanRecDist_Para;
TH1D* ThetanDiffDist_Para;
TH2D* ThetanDiffZp_Para;
TH1D* ThetanCorrDist_Para;
TH1D* ThetanCorrDiffDist_Para;
TH1D* ThetanCorrRecDiffDist_Para;
TH2D* ThetanCorrDiffZp_Para;

TH2D* E_dE_Para;
TH2D* KinEp_dE_Para;
TH2D* E_KinEp_Para;

TH1D* ThetaRecPiDiff_Para;
TH2D* ThetanThetaRecPi_Para;
TH2D* ThetanThetaRecPiDiff_Para;
TH1D* ThetaRecPDiff_Para;
TH2D* ThetanThetaRecP_Para;
TH2D* ThetanThetaRecPDiff_Para;

TH2D* DeutKinPiKin_Para;

TH1D* time_Perp;
TH1D* time_cut_Perp;

TH1D* Eg_Perp;
TH1D* PhiDifference_Perp;
TH1D* EpKin_Perp;
TH1D* EpCorrected_Perp;
TH1D* EpKinEpCorrDiff_Perp;
TH1D* EpEpCorrDiff_Perp;

TH1D* MMpEpCorrected_Perp;
TH1D* OAngle_Perp;
TH1D* WCZnRecon_Perp;

TH1D* MMp200300_Perp;
TH1D* MMp300400_Perp;
TH1D* MMp400500_Perp;
TH1D* MMp500600_Perp;
TH1D* MMp600700_Perp;
TH1D* MMp700800_Perp;
TH1D* MMp800900_Perp;

TH1D* ZpDist_Perp;

TH1D* Phip_435MeVCM1_Perp;
TH1D* Phip_455MeVCM1_Perp;
TH1D* Phip_475MeVCM1_Perp;
TH1D* Phip_495MeVCM1_Perp;
TH1D* Phip_515MeVCM1_Perp;
TH1D* Phip_535MeVCM1_Perp;
TH1D* Phip_555MeVCM1_Perp;
TH1D* Phip_575MeVCM1_Perp;
TH1D* Phip_595MeVCM1_Perp;
TH1D* Phip_615MeVCM1_Perp;

TH1D* Phip_435MeVCM2_Perp;
TH1D* Phip_455MeVCM2_Perp;
TH1D* Phip_475MeVCM2_Perp;
TH1D* Phip_495MeVCM2_Perp;
TH1D* Phip_515MeVCM2_Perp;
TH1D* Phip_535MeVCM2_Perp;
TH1D* Phip_555MeVCM2_Perp;
TH1D* Phip_575MeVCM2_Perp;
TH1D* Phip_595MeVCM2_Perp;
TH1D* Phip_615MeVCM2_Perp;

TH1D* Phip_435MeVCM3_Perp;
TH1D* Phip_455MeVCM3_Perp;
TH1D* Phip_475MeVCM3_Perp;
TH1D* Phip_495MeVCM3_Perp;
TH1D* Phip_515MeVCM3_Perp;
TH1D* Phip_535MeVCM3_Perp;
TH1D* Phip_555MeVCM3_Perp;
TH1D* Phip_575MeVCM3_Perp;
TH1D* Phip_595MeVCM3_Perp;
TH1D* Phip_615MeVCM3_Perp;

TH1D* Phip_435MeVCM4_Perp;
TH1D* Phip_455MeVCM4_Perp;
TH1D* Phip_475MeVCM4_Perp;
TH1D* Phip_495MeVCM4_Perp;
TH1D* Phip_515MeVCM4_Perp;
TH1D* Phip_535MeVCM4_Perp;
TH1D* Phip_555MeVCM4_Perp;
TH1D* Phip_575MeVCM4_Perp;
TH1D* Phip_595MeVCM4_Perp;
TH1D* Phip_615MeVCM4_Perp;

TH1D* Phip_435MeVCM5_Perp;
TH1D* Phip_455MeVCM5_Perp;
TH1D* Phip_475MeVCM5_Perp;
TH1D* Phip_495MeVCM5_Perp;
TH1D* Phip_515MeVCM5_Perp;
TH1D* Phip_535MeVCM5_Perp;
TH1D* Phip_555MeVCM5_Perp;
TH1D* Phip_575MeVCM5_Perp;
TH1D* Phip_595MeVCM5_Perp;
TH1D* Phip_615MeVCM5_Perp;

TH1D* Phip_435MeVCM6_Perp;
TH1D* Phip_455MeVCM6_Perp;
TH1D* Phip_475MeVCM6_Perp;
TH1D* Phip_495MeVCM6_Perp;
TH1D* Phip_515MeVCM6_Perp;
TH1D* Phip_535MeVCM6_Perp;
TH1D* Phip_555MeVCM6_Perp;
TH1D* Phip_575MeVCM6_Perp;
TH1D* Phip_595MeVCM6_Perp;
TH1D* Phip_615MeVCM6_Perp;

TH1D* Phip_435MeVCM7_Perp;
TH1D* Phip_455MeVCM7_Perp;
TH1D* Phip_475MeVCM7_Perp;
TH1D* Phip_495MeVCM7_Perp;
TH1D* Phip_515MeVCM7_Perp;
TH1D* Phip_535MeVCM7_Perp;
TH1D* Phip_555MeVCM7_Perp;
TH1D* Phip_575MeVCM7_Perp;
TH1D* Phip_595MeVCM7_Perp;
TH1D* Phip_615MeVCM7_Perp;

TH1D* Phip_435MeVCM8_Perp;
TH1D* Phip_455MeVCM8_Perp;
TH1D* Phip_475MeVCM8_Perp;
TH1D* Phip_495MeVCM8_Perp;
TH1D* Phip_515MeVCM8_Perp;
TH1D* Phip_535MeVCM8_Perp;
TH1D* Phip_555MeVCM8_Perp;
TH1D* Phip_575MeVCM8_Perp;
TH1D* Phip_595MeVCM8_Perp;
TH1D* Phip_615MeVCM8_Perp;

TH1D* Phip_435MeVCM9_Perp;
TH1D* Phip_455MeVCM9_Perp;
TH1D* Phip_475MeVCM9_Perp;
TH1D* Phip_495MeVCM9_Perp;
TH1D* Phip_515MeVCM9_Perp;
TH1D* Phip_535MeVCM9_Perp;
TH1D* Phip_555MeVCM9_Perp;
TH1D* Phip_575MeVCM9_Perp;
TH1D* Phip_595MeVCM9_Perp;
TH1D* Phip_615MeVCM9_Perp;

TH1D* Phip_435MeVCM10_Perp;
TH1D* Phip_455MeVCM10_Perp;
TH1D* Phip_475MeVCM10_Perp;
TH1D* Phip_495MeVCM10_Perp;
TH1D* Phip_515MeVCM10_Perp;
TH1D* Phip_535MeVCM10_Perp;
TH1D* Phip_555MeVCM10_Perp;
TH1D* Phip_575MeVCM10_Perp;
TH1D* Phip_595MeVCM10_Perp;
TH1D* Phip_615MeVCM10_Perp;

TH1D* ThetanDist_Perp;
TH1D* ThetanRecDist_Perp;
TH1D* ThetanDiffDist_Perp;
TH2D* ThetanDiffZp_Perp;
TH1D* ThetanCorrDist_Perp;
TH1D* ThetanCorrDiffDist_Perp;
TH1D* ThetanCorrRecDiffDist_Perp;
TH2D* ThetanCorrDiffZp_Perp;

TH2D* DeutKinPiKin_Perp;

TH2D* E_dE_Perp;
TH2D* KinEp_dE_Perp;
TH2D* E_KinEp_Perp;

TH1D* ThetaRecPiDiff_Perp;
TH2D* ThetanThetaRecPi_Perp;
TH2D* ThetanThetaRecPiDiff_Perp;
TH1D* ThetaRecPDiff_Perp;
TH2D* ThetanThetaRecP_Perp;
TH2D* ThetanThetaRecPDiff_Perp;
