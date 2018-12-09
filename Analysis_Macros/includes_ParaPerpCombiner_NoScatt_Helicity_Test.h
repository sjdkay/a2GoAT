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

TH1D* Phip_425MeVCM1_Para_NegHel;
TH1D* Phip_435MeVCM1_Para_NegHel;
TH1D* Phip_445MeVCM1_Para_NegHel;
TH1D* Phip_455MeVCM1_Para_NegHel;
TH1D* Phip_465MeVCM1_Para_NegHel;
TH1D* Phip_475MeVCM1_Para_NegHel;
TH1D* Phip_485MeVCM1_Para_NegHel;
TH1D* Phip_495MeVCM1_Para_NegHel;
TH1D* Phip_505MeVCM1_Para_NegHel;
TH1D* Phip_515MeVCM1_Para_NegHel;
TH1D* Phip_525MeVCM1_Para_NegHel;
TH1D* Phip_535MeVCM1_Para_NegHel;
TH1D* Phip_545MeVCM1_Para_NegHel;
TH1D* Phip_555MeVCM1_Para_NegHel;
TH1D* Phip_565MeVCM1_Para_NegHel;
TH1D* Phip_575MeVCM1_Para_NegHel;
TH1D* Phip_585MeVCM1_Para_NegHel;
TH1D* Phip_595MeVCM1_Para_NegHel;
TH1D* Phip_605MeVCM1_Para_NegHel;
TH1D* Phip_615MeVCM1_Para_NegHel;

TH1D* Phip_425MeVCM2_Para_NegHel;
TH1D* Phip_435MeVCM2_Para_NegHel;
TH1D* Phip_445MeVCM2_Para_NegHel;
TH1D* Phip_455MeVCM2_Para_NegHel;
TH1D* Phip_465MeVCM2_Para_NegHel;
TH1D* Phip_475MeVCM2_Para_NegHel;
TH1D* Phip_485MeVCM2_Para_NegHel;
TH1D* Phip_495MeVCM2_Para_NegHel;
TH1D* Phip_505MeVCM2_Para_NegHel;
TH1D* Phip_515MeVCM2_Para_NegHel;
TH1D* Phip_525MeVCM2_Para_NegHel;
TH1D* Phip_535MeVCM2_Para_NegHel;
TH1D* Phip_545MeVCM2_Para_NegHel;
TH1D* Phip_555MeVCM2_Para_NegHel;
TH1D* Phip_565MeVCM2_Para_NegHel;
TH1D* Phip_575MeVCM2_Para_NegHel;
TH1D* Phip_585MeVCM2_Para_NegHel;
TH1D* Phip_595MeVCM2_Para_NegHel;
TH1D* Phip_605MeVCM2_Para_NegHel;
TH1D* Phip_615MeVCM2_Para_NegHel;

TH1D* Phip_425MeVCM3_Para_NegHel;
TH1D* Phip_435MeVCM3_Para_NegHel;
TH1D* Phip_445MeVCM3_Para_NegHel;
TH1D* Phip_455MeVCM3_Para_NegHel;
TH1D* Phip_465MeVCM3_Para_NegHel;
TH1D* Phip_475MeVCM3_Para_NegHel;
TH1D* Phip_485MeVCM3_Para_NegHel;
TH1D* Phip_495MeVCM3_Para_NegHel;
TH1D* Phip_505MeVCM3_Para_NegHel;
TH1D* Phip_515MeVCM3_Para_NegHel;
TH1D* Phip_525MeVCM3_Para_NegHel;
TH1D* Phip_535MeVCM3_Para_NegHel;
TH1D* Phip_545MeVCM3_Para_NegHel;
TH1D* Phip_555MeVCM3_Para_NegHel;
TH1D* Phip_565MeVCM3_Para_NegHel;
TH1D* Phip_575MeVCM3_Para_NegHel;
TH1D* Phip_585MeVCM3_Para_NegHel;
TH1D* Phip_595MeVCM3_Para_NegHel;
TH1D* Phip_605MeVCM3_Para_NegHel;
TH1D* Phip_615MeVCM3_Para_NegHel;

TH1D* Phip_425MeVCM4_Para_NegHel;
TH1D* Phip_435MeVCM4_Para_NegHel;
TH1D* Phip_445MeVCM4_Para_NegHel;
TH1D* Phip_455MeVCM4_Para_NegHel;
TH1D* Phip_465MeVCM4_Para_NegHel;
TH1D* Phip_475MeVCM4_Para_NegHel;
TH1D* Phip_485MeVCM4_Para_NegHel;
TH1D* Phip_495MeVCM4_Para_NegHel;
TH1D* Phip_505MeVCM4_Para_NegHel;
TH1D* Phip_515MeVCM4_Para_NegHel;
TH1D* Phip_525MeVCM4_Para_NegHel;
TH1D* Phip_535MeVCM4_Para_NegHel;
TH1D* Phip_545MeVCM4_Para_NegHel;
TH1D* Phip_555MeVCM4_Para_NegHel;
TH1D* Phip_565MeVCM4_Para_NegHel;
TH1D* Phip_575MeVCM4_Para_NegHel;
TH1D* Phip_585MeVCM4_Para_NegHel;
TH1D* Phip_595MeVCM4_Para_NegHel;
TH1D* Phip_605MeVCM4_Para_NegHel;
TH1D* Phip_615MeVCM4_Para_NegHel;

TH1D* Phip_425MeVCM5_Para_NegHel;
TH1D* Phip_435MeVCM5_Para_NegHel;
TH1D* Phip_445MeVCM5_Para_NegHel;
TH1D* Phip_455MeVCM5_Para_NegHel;
TH1D* Phip_465MeVCM5_Para_NegHel;
TH1D* Phip_475MeVCM5_Para_NegHel;
TH1D* Phip_485MeVCM5_Para_NegHel;
TH1D* Phip_495MeVCM5_Para_NegHel;
TH1D* Phip_505MeVCM5_Para_NegHel;
TH1D* Phip_515MeVCM5_Para_NegHel;
TH1D* Phip_525MeVCM5_Para_NegHel;
TH1D* Phip_535MeVCM5_Para_NegHel;
TH1D* Phip_545MeVCM5_Para_NegHel;
TH1D* Phip_555MeVCM5_Para_NegHel;
TH1D* Phip_565MeVCM5_Para_NegHel;
TH1D* Phip_575MeVCM5_Para_NegHel;
TH1D* Phip_585MeVCM5_Para_NegHel;
TH1D* Phip_595MeVCM5_Para_NegHel;
TH1D* Phip_605MeVCM5_Para_NegHel;
TH1D* Phip_615MeVCM5_Para_NegHel;

TH1D* Phip_425MeVCM6_Para_NegHel;
TH1D* Phip_435MeVCM6_Para_NegHel;
TH1D* Phip_445MeVCM6_Para_NegHel;
TH1D* Phip_455MeVCM6_Para_NegHel;
TH1D* Phip_465MeVCM6_Para_NegHel;
TH1D* Phip_475MeVCM6_Para_NegHel;
TH1D* Phip_485MeVCM6_Para_NegHel;
TH1D* Phip_495MeVCM6_Para_NegHel;
TH1D* Phip_505MeVCM6_Para_NegHel;
TH1D* Phip_515MeVCM6_Para_NegHel;
TH1D* Phip_525MeVCM6_Para_NegHel;
TH1D* Phip_535MeVCM6_Para_NegHel;
TH1D* Phip_545MeVCM6_Para_NegHel;
TH1D* Phip_555MeVCM6_Para_NegHel;
TH1D* Phip_565MeVCM6_Para_NegHel;
TH1D* Phip_575MeVCM6_Para_NegHel;
TH1D* Phip_585MeVCM6_Para_NegHel;
TH1D* Phip_595MeVCM6_Para_NegHel;
TH1D* Phip_605MeVCM6_Para_NegHel;
TH1D* Phip_615MeVCM6_Para_NegHel;

TH1D* Phip_425MeVCM7_Para_NegHel;
TH1D* Phip_435MeVCM7_Para_NegHel;
TH1D* Phip_445MeVCM7_Para_NegHel;
TH1D* Phip_455MeVCM7_Para_NegHel;
TH1D* Phip_465MeVCM7_Para_NegHel;
TH1D* Phip_475MeVCM7_Para_NegHel;
TH1D* Phip_485MeVCM7_Para_NegHel;
TH1D* Phip_495MeVCM7_Para_NegHel;
TH1D* Phip_505MeVCM7_Para_NegHel;
TH1D* Phip_515MeVCM7_Para_NegHel;
TH1D* Phip_525MeVCM7_Para_NegHel;
TH1D* Phip_535MeVCM7_Para_NegHel;
TH1D* Phip_545MeVCM7_Para_NegHel;
TH1D* Phip_555MeVCM7_Para_NegHel;
TH1D* Phip_565MeVCM7_Para_NegHel;
TH1D* Phip_575MeVCM7_Para_NegHel;
TH1D* Phip_585MeVCM7_Para_NegHel;
TH1D* Phip_595MeVCM7_Para_NegHel;
TH1D* Phip_605MeVCM7_Para_NegHel;
TH1D* Phip_615MeVCM7_Para_NegHel;

TH1D* Phip_425MeVCM8_Para_NegHel;
TH1D* Phip_435MeVCM8_Para_NegHel;
TH1D* Phip_445MeVCM8_Para_NegHel;
TH1D* Phip_455MeVCM8_Para_NegHel;
TH1D* Phip_465MeVCM8_Para_NegHel;
TH1D* Phip_475MeVCM8_Para_NegHel;
TH1D* Phip_485MeVCM8_Para_NegHel;
TH1D* Phip_495MeVCM8_Para_NegHel;
TH1D* Phip_505MeVCM8_Para_NegHel;
TH1D* Phip_515MeVCM8_Para_NegHel;
TH1D* Phip_525MeVCM8_Para_NegHel;
TH1D* Phip_535MeVCM8_Para_NegHel;
TH1D* Phip_545MeVCM8_Para_NegHel;
TH1D* Phip_555MeVCM8_Para_NegHel;
TH1D* Phip_565MeVCM8_Para_NegHel;
TH1D* Phip_575MeVCM8_Para_NegHel;
TH1D* Phip_585MeVCM8_Para_NegHel;
TH1D* Phip_595MeVCM8_Para_NegHel;
TH1D* Phip_605MeVCM8_Para_NegHel;
TH1D* Phip_615MeVCM8_Para_NegHel;

TH1D* Phip_425MeVCM9_Para_NegHel;
TH1D* Phip_435MeVCM9_Para_NegHel;
TH1D* Phip_445MeVCM9_Para_NegHel;
TH1D* Phip_455MeVCM9_Para_NegHel;
TH1D* Phip_465MeVCM9_Para_NegHel;
TH1D* Phip_475MeVCM9_Para_NegHel;
TH1D* Phip_485MeVCM9_Para_NegHel;
TH1D* Phip_495MeVCM9_Para_NegHel;
TH1D* Phip_505MeVCM9_Para_NegHel;
TH1D* Phip_515MeVCM9_Para_NegHel;
TH1D* Phip_525MeVCM9_Para_NegHel;
TH1D* Phip_535MeVCM9_Para_NegHel;
TH1D* Phip_545MeVCM9_Para_NegHel;
TH1D* Phip_555MeVCM9_Para_NegHel;
TH1D* Phip_565MeVCM9_Para_NegHel;
TH1D* Phip_575MeVCM9_Para_NegHel;
TH1D* Phip_585MeVCM9_Para_NegHel;
TH1D* Phip_595MeVCM9_Para_NegHel;
TH1D* Phip_605MeVCM9_Para_NegHel;
TH1D* Phip_615MeVCM9_Para_NegHel;

TH1D* Phip_425MeVCM10_Para_NegHel;
TH1D* Phip_435MeVCM10_Para_NegHel;
TH1D* Phip_445MeVCM10_Para_NegHel;
TH1D* Phip_455MeVCM10_Para_NegHel;
TH1D* Phip_465MeVCM10_Para_NegHel;
TH1D* Phip_475MeVCM10_Para_NegHel;
TH1D* Phip_485MeVCM10_Para_NegHel;
TH1D* Phip_495MeVCM10_Para_NegHel;
TH1D* Phip_505MeVCM10_Para_NegHel;
TH1D* Phip_515MeVCM10_Para_NegHel;
TH1D* Phip_525MeVCM10_Para_NegHel;
TH1D* Phip_535MeVCM10_Para_NegHel;
TH1D* Phip_545MeVCM10_Para_NegHel;
TH1D* Phip_555MeVCM10_Para_NegHel;
TH1D* Phip_565MeVCM10_Para_NegHel;
TH1D* Phip_575MeVCM10_Para_NegHel;
TH1D* Phip_585MeVCM10_Para_NegHel;
TH1D* Phip_595MeVCM10_Para_NegHel;
TH1D* Phip_605MeVCM10_Para_NegHel;
TH1D* Phip_615MeVCM10_Para_NegHel;

TH1D* Phip_425MeVCM1_Para_PosHel;
TH1D* Phip_435MeVCM1_Para_PosHel;
TH1D* Phip_445MeVCM1_Para_PosHel;
TH1D* Phip_455MeVCM1_Para_PosHel;
TH1D* Phip_465MeVCM1_Para_PosHel;
TH1D* Phip_475MeVCM1_Para_PosHel;
TH1D* Phip_485MeVCM1_Para_PosHel;
TH1D* Phip_495MeVCM1_Para_PosHel;
TH1D* Phip_505MeVCM1_Para_PosHel;
TH1D* Phip_515MeVCM1_Para_PosHel;
TH1D* Phip_525MeVCM1_Para_PosHel;
TH1D* Phip_535MeVCM1_Para_PosHel;
TH1D* Phip_545MeVCM1_Para_PosHel;
TH1D* Phip_555MeVCM1_Para_PosHel;
TH1D* Phip_565MeVCM1_Para_PosHel;
TH1D* Phip_575MeVCM1_Para_PosHel;
TH1D* Phip_585MeVCM1_Para_PosHel;
TH1D* Phip_595MeVCM1_Para_PosHel;
TH1D* Phip_605MeVCM1_Para_PosHel;
TH1D* Phip_615MeVCM1_Para_PosHel;

TH1D* Phip_425MeVCM2_Para_PosHel;
TH1D* Phip_435MeVCM2_Para_PosHel;
TH1D* Phip_445MeVCM2_Para_PosHel;
TH1D* Phip_455MeVCM2_Para_PosHel;
TH1D* Phip_465MeVCM2_Para_PosHel;
TH1D* Phip_475MeVCM2_Para_PosHel;
TH1D* Phip_485MeVCM2_Para_PosHel;
TH1D* Phip_495MeVCM2_Para_PosHel;
TH1D* Phip_505MeVCM2_Para_PosHel;
TH1D* Phip_515MeVCM2_Para_PosHel;
TH1D* Phip_525MeVCM2_Para_PosHel;
TH1D* Phip_535MeVCM2_Para_PosHel;
TH1D* Phip_545MeVCM2_Para_PosHel;
TH1D* Phip_555MeVCM2_Para_PosHel;
TH1D* Phip_565MeVCM2_Para_PosHel;
TH1D* Phip_575MeVCM2_Para_PosHel;
TH1D* Phip_585MeVCM2_Para_PosHel;
TH1D* Phip_595MeVCM2_Para_PosHel;
TH1D* Phip_605MeVCM2_Para_PosHel;
TH1D* Phip_615MeVCM2_Para_PosHel;

TH1D* Phip_425MeVCM3_Para_PosHel;
TH1D* Phip_435MeVCM3_Para_PosHel;
TH1D* Phip_445MeVCM3_Para_PosHel;
TH1D* Phip_455MeVCM3_Para_PosHel;
TH1D* Phip_465MeVCM3_Para_PosHel;
TH1D* Phip_475MeVCM3_Para_PosHel;
TH1D* Phip_485MeVCM3_Para_PosHel;
TH1D* Phip_495MeVCM3_Para_PosHel;
TH1D* Phip_505MeVCM3_Para_PosHel;
TH1D* Phip_515MeVCM3_Para_PosHel;
TH1D* Phip_525MeVCM3_Para_PosHel;
TH1D* Phip_535MeVCM3_Para_PosHel;
TH1D* Phip_545MeVCM3_Para_PosHel;
TH1D* Phip_555MeVCM3_Para_PosHel;
TH1D* Phip_565MeVCM3_Para_PosHel;
TH1D* Phip_575MeVCM3_Para_PosHel;
TH1D* Phip_585MeVCM3_Para_PosHel;
TH1D* Phip_595MeVCM3_Para_PosHel;
TH1D* Phip_605MeVCM3_Para_PosHel;
TH1D* Phip_615MeVCM3_Para_PosHel;

TH1D* Phip_425MeVCM4_Para_PosHel;
TH1D* Phip_435MeVCM4_Para_PosHel;
TH1D* Phip_445MeVCM4_Para_PosHel;
TH1D* Phip_455MeVCM4_Para_PosHel;
TH1D* Phip_465MeVCM4_Para_PosHel;
TH1D* Phip_475MeVCM4_Para_PosHel;
TH1D* Phip_485MeVCM4_Para_PosHel;
TH1D* Phip_495MeVCM4_Para_PosHel;
TH1D* Phip_505MeVCM4_Para_PosHel;
TH1D* Phip_515MeVCM4_Para_PosHel;
TH1D* Phip_525MeVCM4_Para_PosHel;
TH1D* Phip_535MeVCM4_Para_PosHel;
TH1D* Phip_545MeVCM4_Para_PosHel;
TH1D* Phip_555MeVCM4_Para_PosHel;
TH1D* Phip_565MeVCM4_Para_PosHel;
TH1D* Phip_575MeVCM4_Para_PosHel;
TH1D* Phip_585MeVCM4_Para_PosHel;
TH1D* Phip_595MeVCM4_Para_PosHel;
TH1D* Phip_605MeVCM4_Para_PosHel;
TH1D* Phip_615MeVCM4_Para_PosHel;

TH1D* Phip_425MeVCM5_Para_PosHel;
TH1D* Phip_435MeVCM5_Para_PosHel;
TH1D* Phip_445MeVCM5_Para_PosHel;
TH1D* Phip_455MeVCM5_Para_PosHel;
TH1D* Phip_465MeVCM5_Para_PosHel;
TH1D* Phip_475MeVCM5_Para_PosHel;
TH1D* Phip_485MeVCM5_Para_PosHel;
TH1D* Phip_495MeVCM5_Para_PosHel;
TH1D* Phip_505MeVCM5_Para_PosHel;
TH1D* Phip_515MeVCM5_Para_PosHel;
TH1D* Phip_525MeVCM5_Para_PosHel;
TH1D* Phip_535MeVCM5_Para_PosHel;
TH1D* Phip_545MeVCM5_Para_PosHel;
TH1D* Phip_555MeVCM5_Para_PosHel;
TH1D* Phip_565MeVCM5_Para_PosHel;
TH1D* Phip_575MeVCM5_Para_PosHel;
TH1D* Phip_585MeVCM5_Para_PosHel;
TH1D* Phip_595MeVCM5_Para_PosHel;
TH1D* Phip_605MeVCM5_Para_PosHel;
TH1D* Phip_615MeVCM5_Para_PosHel;

TH1D* Phip_425MeVCM6_Para_PosHel;
TH1D* Phip_435MeVCM6_Para_PosHel;
TH1D* Phip_445MeVCM6_Para_PosHel;
TH1D* Phip_455MeVCM6_Para_PosHel;
TH1D* Phip_465MeVCM6_Para_PosHel;
TH1D* Phip_475MeVCM6_Para_PosHel;
TH1D* Phip_485MeVCM6_Para_PosHel;
TH1D* Phip_495MeVCM6_Para_PosHel;
TH1D* Phip_505MeVCM6_Para_PosHel;
TH1D* Phip_515MeVCM6_Para_PosHel;
TH1D* Phip_525MeVCM6_Para_PosHel;
TH1D* Phip_535MeVCM6_Para_PosHel;
TH1D* Phip_545MeVCM6_Para_PosHel;
TH1D* Phip_555MeVCM6_Para_PosHel;
TH1D* Phip_565MeVCM6_Para_PosHel;
TH1D* Phip_575MeVCM6_Para_PosHel;
TH1D* Phip_585MeVCM6_Para_PosHel;
TH1D* Phip_595MeVCM6_Para_PosHel;
TH1D* Phip_605MeVCM6_Para_PosHel;
TH1D* Phip_615MeVCM6_Para_PosHel;

TH1D* Phip_425MeVCM7_Para_PosHel;
TH1D* Phip_435MeVCM7_Para_PosHel;
TH1D* Phip_445MeVCM7_Para_PosHel;
TH1D* Phip_455MeVCM7_Para_PosHel;
TH1D* Phip_465MeVCM7_Para_PosHel;
TH1D* Phip_475MeVCM7_Para_PosHel;
TH1D* Phip_485MeVCM7_Para_PosHel;
TH1D* Phip_495MeVCM7_Para_PosHel;
TH1D* Phip_505MeVCM7_Para_PosHel;
TH1D* Phip_515MeVCM7_Para_PosHel;
TH1D* Phip_525MeVCM7_Para_PosHel;
TH1D* Phip_535MeVCM7_Para_PosHel;
TH1D* Phip_545MeVCM7_Para_PosHel;
TH1D* Phip_555MeVCM7_Para_PosHel;
TH1D* Phip_565MeVCM7_Para_PosHel;
TH1D* Phip_575MeVCM7_Para_PosHel;
TH1D* Phip_585MeVCM7_Para_PosHel;
TH1D* Phip_595MeVCM7_Para_PosHel;
TH1D* Phip_605MeVCM7_Para_PosHel;
TH1D* Phip_615MeVCM7_Para_PosHel;

TH1D* Phip_425MeVCM8_Para_PosHel;
TH1D* Phip_435MeVCM8_Para_PosHel;
TH1D* Phip_445MeVCM8_Para_PosHel;
TH1D* Phip_455MeVCM8_Para_PosHel;
TH1D* Phip_465MeVCM8_Para_PosHel;
TH1D* Phip_475MeVCM8_Para_PosHel;
TH1D* Phip_485MeVCM8_Para_PosHel;
TH1D* Phip_495MeVCM8_Para_PosHel;
TH1D* Phip_505MeVCM8_Para_PosHel;
TH1D* Phip_515MeVCM8_Para_PosHel;
TH1D* Phip_525MeVCM8_Para_PosHel;
TH1D* Phip_535MeVCM8_Para_PosHel;
TH1D* Phip_545MeVCM8_Para_PosHel;
TH1D* Phip_555MeVCM8_Para_PosHel;
TH1D* Phip_565MeVCM8_Para_PosHel;
TH1D* Phip_575MeVCM8_Para_PosHel;
TH1D* Phip_585MeVCM8_Para_PosHel;
TH1D* Phip_595MeVCM8_Para_PosHel;
TH1D* Phip_605MeVCM8_Para_PosHel;
TH1D* Phip_615MeVCM8_Para_PosHel;

TH1D* Phip_425MeVCM9_Para_PosHel;
TH1D* Phip_435MeVCM9_Para_PosHel;
TH1D* Phip_445MeVCM9_Para_PosHel;
TH1D* Phip_455MeVCM9_Para_PosHel;
TH1D* Phip_465MeVCM9_Para_PosHel;
TH1D* Phip_475MeVCM9_Para_PosHel;
TH1D* Phip_485MeVCM9_Para_PosHel;
TH1D* Phip_495MeVCM9_Para_PosHel;
TH1D* Phip_505MeVCM9_Para_PosHel;
TH1D* Phip_515MeVCM9_Para_PosHel;
TH1D* Phip_525MeVCM9_Para_PosHel;
TH1D* Phip_535MeVCM9_Para_PosHel;
TH1D* Phip_545MeVCM9_Para_PosHel;
TH1D* Phip_555MeVCM9_Para_PosHel;
TH1D* Phip_565MeVCM9_Para_PosHel;
TH1D* Phip_575MeVCM9_Para_PosHel;
TH1D* Phip_585MeVCM9_Para_PosHel;
TH1D* Phip_595MeVCM9_Para_PosHel;
TH1D* Phip_605MeVCM9_Para_PosHel;
TH1D* Phip_615MeVCM9_Para_PosHel;

TH1D* Phip_425MeVCM10_Para_PosHel;
TH1D* Phip_435MeVCM10_Para_PosHel;
TH1D* Phip_445MeVCM10_Para_PosHel;
TH1D* Phip_455MeVCM10_Para_PosHel;
TH1D* Phip_465MeVCM10_Para_PosHel;
TH1D* Phip_475MeVCM10_Para_PosHel;
TH1D* Phip_485MeVCM10_Para_PosHel;
TH1D* Phip_495MeVCM10_Para_PosHel;
TH1D* Phip_505MeVCM10_Para_PosHel;
TH1D* Phip_515MeVCM10_Para_PosHel;
TH1D* Phip_525MeVCM10_Para_PosHel;
TH1D* Phip_535MeVCM10_Para_PosHel;
TH1D* Phip_545MeVCM10_Para_PosHel;
TH1D* Phip_555MeVCM10_Para_PosHel;
TH1D* Phip_565MeVCM10_Para_PosHel;
TH1D* Phip_575MeVCM10_Para_PosHel;
TH1D* Phip_585MeVCM10_Para_PosHel;
TH1D* Phip_595MeVCM10_Para_PosHel;
TH1D* Phip_605MeVCM10_Para_PosHel;
TH1D* Phip_615MeVCM10_Para_PosHel;

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

TH1D* Phip_425MeVCM1_Perp_NegHel;
TH1D* Phip_435MeVCM1_Perp_NegHel;
TH1D* Phip_445MeVCM1_Perp_NegHel;
TH1D* Phip_455MeVCM1_Perp_NegHel;
TH1D* Phip_465MeVCM1_Perp_NegHel;
TH1D* Phip_475MeVCM1_Perp_NegHel;
TH1D* Phip_485MeVCM1_Perp_NegHel;
TH1D* Phip_495MeVCM1_Perp_NegHel;
TH1D* Phip_505MeVCM1_Perp_NegHel;
TH1D* Phip_515MeVCM1_Perp_NegHel;
TH1D* Phip_525MeVCM1_Perp_NegHel;
TH1D* Phip_535MeVCM1_Perp_NegHel;
TH1D* Phip_545MeVCM1_Perp_NegHel;
TH1D* Phip_555MeVCM1_Perp_NegHel;
TH1D* Phip_565MeVCM1_Perp_NegHel;
TH1D* Phip_575MeVCM1_Perp_NegHel;
TH1D* Phip_585MeVCM1_Perp_NegHel;
TH1D* Phip_595MeVCM1_Perp_NegHel;
TH1D* Phip_605MeVCM1_Perp_NegHel;
TH1D* Phip_615MeVCM1_Perp_NegHel;

TH1D* Phip_425MeVCM2_Perp_NegHel;
TH1D* Phip_435MeVCM2_Perp_NegHel;
TH1D* Phip_445MeVCM2_Perp_NegHel;
TH1D* Phip_455MeVCM2_Perp_NegHel;
TH1D* Phip_465MeVCM2_Perp_NegHel;
TH1D* Phip_475MeVCM2_Perp_NegHel;
TH1D* Phip_485MeVCM2_Perp_NegHel;
TH1D* Phip_495MeVCM2_Perp_NegHel;
TH1D* Phip_505MeVCM2_Perp_NegHel;
TH1D* Phip_515MeVCM2_Perp_NegHel;
TH1D* Phip_525MeVCM2_Perp_NegHel;
TH1D* Phip_535MeVCM2_Perp_NegHel;
TH1D* Phip_545MeVCM2_Perp_NegHel;
TH1D* Phip_555MeVCM2_Perp_NegHel;
TH1D* Phip_565MeVCM2_Perp_NegHel;
TH1D* Phip_575MeVCM2_Perp_NegHel;
TH1D* Phip_585MeVCM2_Perp_NegHel;
TH1D* Phip_595MeVCM2_Perp_NegHel;
TH1D* Phip_605MeVCM2_Perp_NegHel;
TH1D* Phip_615MeVCM2_Perp_NegHel;

TH1D* Phip_425MeVCM3_Perp_NegHel;
TH1D* Phip_435MeVCM3_Perp_NegHel;
TH1D* Phip_445MeVCM3_Perp_NegHel;
TH1D* Phip_455MeVCM3_Perp_NegHel;
TH1D* Phip_465MeVCM3_Perp_NegHel;
TH1D* Phip_475MeVCM3_Perp_NegHel;
TH1D* Phip_485MeVCM3_Perp_NegHel;
TH1D* Phip_495MeVCM3_Perp_NegHel;
TH1D* Phip_505MeVCM3_Perp_NegHel;
TH1D* Phip_515MeVCM3_Perp_NegHel;
TH1D* Phip_525MeVCM3_Perp_NegHel;
TH1D* Phip_535MeVCM3_Perp_NegHel;
TH1D* Phip_545MeVCM3_Perp_NegHel;
TH1D* Phip_555MeVCM3_Perp_NegHel;
TH1D* Phip_565MeVCM3_Perp_NegHel;
TH1D* Phip_575MeVCM3_Perp_NegHel;
TH1D* Phip_585MeVCM3_Perp_NegHel;
TH1D* Phip_595MeVCM3_Perp_NegHel;
TH1D* Phip_605MeVCM3_Perp_NegHel;
TH1D* Phip_615MeVCM3_Perp_NegHel;

TH1D* Phip_425MeVCM4_Perp_NegHel;
TH1D* Phip_435MeVCM4_Perp_NegHel;
TH1D* Phip_445MeVCM4_Perp_NegHel;
TH1D* Phip_455MeVCM4_Perp_NegHel;
TH1D* Phip_465MeVCM4_Perp_NegHel;
TH1D* Phip_475MeVCM4_Perp_NegHel;
TH1D* Phip_485MeVCM4_Perp_NegHel;
TH1D* Phip_495MeVCM4_Perp_NegHel;
TH1D* Phip_505MeVCM4_Perp_NegHel;
TH1D* Phip_515MeVCM4_Perp_NegHel;
TH1D* Phip_525MeVCM4_Perp_NegHel;
TH1D* Phip_535MeVCM4_Perp_NegHel;
TH1D* Phip_545MeVCM4_Perp_NegHel;
TH1D* Phip_555MeVCM4_Perp_NegHel;
TH1D* Phip_565MeVCM4_Perp_NegHel;
TH1D* Phip_575MeVCM4_Perp_NegHel;
TH1D* Phip_585MeVCM4_Perp_NegHel;
TH1D* Phip_595MeVCM4_Perp_NegHel;
TH1D* Phip_605MeVCM4_Perp_NegHel;
TH1D* Phip_615MeVCM4_Perp_NegHel;

TH1D* Phip_425MeVCM5_Perp_NegHel;
TH1D* Phip_435MeVCM5_Perp_NegHel;
TH1D* Phip_445MeVCM5_Perp_NegHel;
TH1D* Phip_455MeVCM5_Perp_NegHel;
TH1D* Phip_465MeVCM5_Perp_NegHel;
TH1D* Phip_475MeVCM5_Perp_NegHel;
TH1D* Phip_485MeVCM5_Perp_NegHel;
TH1D* Phip_495MeVCM5_Perp_NegHel;
TH1D* Phip_505MeVCM5_Perp_NegHel;
TH1D* Phip_515MeVCM5_Perp_NegHel;
TH1D* Phip_525MeVCM5_Perp_NegHel;
TH1D* Phip_535MeVCM5_Perp_NegHel;
TH1D* Phip_545MeVCM5_Perp_NegHel;
TH1D* Phip_555MeVCM5_Perp_NegHel;
TH1D* Phip_565MeVCM5_Perp_NegHel;
TH1D* Phip_575MeVCM5_Perp_NegHel;
TH1D* Phip_585MeVCM5_Perp_NegHel;
TH1D* Phip_595MeVCM5_Perp_NegHel;
TH1D* Phip_605MeVCM5_Perp_NegHel;
TH1D* Phip_615MeVCM5_Perp_NegHel;

TH1D* Phip_425MeVCM6_Perp_NegHel;
TH1D* Phip_435MeVCM6_Perp_NegHel;
TH1D* Phip_445MeVCM6_Perp_NegHel;
TH1D* Phip_455MeVCM6_Perp_NegHel;
TH1D* Phip_465MeVCM6_Perp_NegHel;
TH1D* Phip_475MeVCM6_Perp_NegHel;
TH1D* Phip_485MeVCM6_Perp_NegHel;
TH1D* Phip_495MeVCM6_Perp_NegHel;
TH1D* Phip_505MeVCM6_Perp_NegHel;
TH1D* Phip_515MeVCM6_Perp_NegHel;
TH1D* Phip_525MeVCM6_Perp_NegHel;
TH1D* Phip_535MeVCM6_Perp_NegHel;
TH1D* Phip_545MeVCM6_Perp_NegHel;
TH1D* Phip_555MeVCM6_Perp_NegHel;
TH1D* Phip_565MeVCM6_Perp_NegHel;
TH1D* Phip_575MeVCM6_Perp_NegHel;
TH1D* Phip_585MeVCM6_Perp_NegHel;
TH1D* Phip_595MeVCM6_Perp_NegHel;
TH1D* Phip_605MeVCM6_Perp_NegHel;
TH1D* Phip_615MeVCM6_Perp_NegHel;

TH1D* Phip_425MeVCM7_Perp_NegHel;
TH1D* Phip_435MeVCM7_Perp_NegHel;
TH1D* Phip_445MeVCM7_Perp_NegHel;
TH1D* Phip_455MeVCM7_Perp_NegHel;
TH1D* Phip_465MeVCM7_Perp_NegHel;
TH1D* Phip_475MeVCM7_Perp_NegHel;
TH1D* Phip_485MeVCM7_Perp_NegHel;
TH1D* Phip_495MeVCM7_Perp_NegHel;
TH1D* Phip_505MeVCM7_Perp_NegHel;
TH1D* Phip_515MeVCM7_Perp_NegHel;
TH1D* Phip_525MeVCM7_Perp_NegHel;
TH1D* Phip_535MeVCM7_Perp_NegHel;
TH1D* Phip_545MeVCM7_Perp_NegHel;
TH1D* Phip_555MeVCM7_Perp_NegHel;
TH1D* Phip_565MeVCM7_Perp_NegHel;
TH1D* Phip_575MeVCM7_Perp_NegHel;
TH1D* Phip_585MeVCM7_Perp_NegHel;
TH1D* Phip_595MeVCM7_Perp_NegHel;
TH1D* Phip_605MeVCM7_Perp_NegHel;
TH1D* Phip_615MeVCM7_Perp_NegHel;

TH1D* Phip_425MeVCM8_Perp_NegHel;
TH1D* Phip_435MeVCM8_Perp_NegHel;
TH1D* Phip_445MeVCM8_Perp_NegHel;
TH1D* Phip_455MeVCM8_Perp_NegHel;
TH1D* Phip_465MeVCM8_Perp_NegHel;
TH1D* Phip_475MeVCM8_Perp_NegHel;
TH1D* Phip_485MeVCM8_Perp_NegHel;
TH1D* Phip_495MeVCM8_Perp_NegHel;
TH1D* Phip_505MeVCM8_Perp_NegHel;
TH1D* Phip_515MeVCM8_Perp_NegHel;
TH1D* Phip_525MeVCM8_Perp_NegHel;
TH1D* Phip_535MeVCM8_Perp_NegHel;
TH1D* Phip_545MeVCM8_Perp_NegHel;
TH1D* Phip_555MeVCM8_Perp_NegHel;
TH1D* Phip_565MeVCM8_Perp_NegHel;
TH1D* Phip_575MeVCM8_Perp_NegHel;
TH1D* Phip_585MeVCM8_Perp_NegHel;
TH1D* Phip_595MeVCM8_Perp_NegHel;
TH1D* Phip_605MeVCM8_Perp_NegHel;
TH1D* Phip_615MeVCM8_Perp_NegHel;

TH1D* Phip_425MeVCM9_Perp_NegHel;
TH1D* Phip_435MeVCM9_Perp_NegHel;
TH1D* Phip_445MeVCM9_Perp_NegHel;
TH1D* Phip_455MeVCM9_Perp_NegHel;
TH1D* Phip_465MeVCM9_Perp_NegHel;
TH1D* Phip_475MeVCM9_Perp_NegHel;
TH1D* Phip_485MeVCM9_Perp_NegHel;
TH1D* Phip_495MeVCM9_Perp_NegHel;
TH1D* Phip_505MeVCM9_Perp_NegHel;
TH1D* Phip_515MeVCM9_Perp_NegHel;
TH1D* Phip_525MeVCM9_Perp_NegHel;
TH1D* Phip_535MeVCM9_Perp_NegHel;
TH1D* Phip_545MeVCM9_Perp_NegHel;
TH1D* Phip_555MeVCM9_Perp_NegHel;
TH1D* Phip_565MeVCM9_Perp_NegHel;
TH1D* Phip_575MeVCM9_Perp_NegHel;
TH1D* Phip_585MeVCM9_Perp_NegHel;
TH1D* Phip_595MeVCM9_Perp_NegHel;
TH1D* Phip_605MeVCM9_Perp_NegHel;
TH1D* Phip_615MeVCM9_Perp_NegHel;

TH1D* Phip_425MeVCM10_Perp_NegHel;
TH1D* Phip_435MeVCM10_Perp_NegHel;
TH1D* Phip_445MeVCM10_Perp_NegHel;
TH1D* Phip_455MeVCM10_Perp_NegHel;
TH1D* Phip_465MeVCM10_Perp_NegHel;
TH1D* Phip_475MeVCM10_Perp_NegHel;
TH1D* Phip_485MeVCM10_Perp_NegHel;
TH1D* Phip_495MeVCM10_Perp_NegHel;
TH1D* Phip_505MeVCM10_Perp_NegHel;
TH1D* Phip_515MeVCM10_Perp_NegHel;
TH1D* Phip_525MeVCM10_Perp_NegHel;
TH1D* Phip_535MeVCM10_Perp_NegHel;
TH1D* Phip_545MeVCM10_Perp_NegHel;
TH1D* Phip_555MeVCM10_Perp_NegHel;
TH1D* Phip_565MeVCM10_Perp_NegHel;
TH1D* Phip_575MeVCM10_Perp_NegHel;
TH1D* Phip_585MeVCM10_Perp_NegHel;
TH1D* Phip_595MeVCM10_Perp_NegHel;
TH1D* Phip_605MeVCM10_Perp_NegHel;
TH1D* Phip_615MeVCM10_Perp_NegHel;

TH1D* Phip_425MeVCM1_Perp_PosHel;
TH1D* Phip_435MeVCM1_Perp_PosHel;
TH1D* Phip_445MeVCM1_Perp_PosHel;
TH1D* Phip_455MeVCM1_Perp_PosHel;
TH1D* Phip_465MeVCM1_Perp_PosHel;
TH1D* Phip_475MeVCM1_Perp_PosHel;
TH1D* Phip_485MeVCM1_Perp_PosHel;
TH1D* Phip_495MeVCM1_Perp_PosHel;
TH1D* Phip_505MeVCM1_Perp_PosHel;
TH1D* Phip_515MeVCM1_Perp_PosHel;
TH1D* Phip_525MeVCM1_Perp_PosHel;
TH1D* Phip_535MeVCM1_Perp_PosHel;
TH1D* Phip_545MeVCM1_Perp_PosHel;
TH1D* Phip_555MeVCM1_Perp_PosHel;
TH1D* Phip_565MeVCM1_Perp_PosHel;
TH1D* Phip_575MeVCM1_Perp_PosHel;
TH1D* Phip_585MeVCM1_Perp_PosHel;
TH1D* Phip_595MeVCM1_Perp_PosHel;
TH1D* Phip_605MeVCM1_Perp_PosHel;
TH1D* Phip_615MeVCM1_Perp_PosHel;

TH1D* Phip_425MeVCM2_Perp_PosHel;
TH1D* Phip_435MeVCM2_Perp_PosHel;
TH1D* Phip_445MeVCM2_Perp_PosHel;
TH1D* Phip_455MeVCM2_Perp_PosHel;
TH1D* Phip_465MeVCM2_Perp_PosHel;
TH1D* Phip_475MeVCM2_Perp_PosHel;
TH1D* Phip_485MeVCM2_Perp_PosHel;
TH1D* Phip_495MeVCM2_Perp_PosHel;
TH1D* Phip_505MeVCM2_Perp_PosHel;
TH1D* Phip_515MeVCM2_Perp_PosHel;
TH1D* Phip_525MeVCM2_Perp_PosHel;
TH1D* Phip_535MeVCM2_Perp_PosHel;
TH1D* Phip_545MeVCM2_Perp_PosHel;
TH1D* Phip_555MeVCM2_Perp_PosHel;
TH1D* Phip_565MeVCM2_Perp_PosHel;
TH1D* Phip_575MeVCM2_Perp_PosHel;
TH1D* Phip_585MeVCM2_Perp_PosHel;
TH1D* Phip_595MeVCM2_Perp_PosHel;
TH1D* Phip_605MeVCM2_Perp_PosHel;
TH1D* Phip_615MeVCM2_Perp_PosHel;

TH1D* Phip_425MeVCM3_Perp_PosHel;
TH1D* Phip_435MeVCM3_Perp_PosHel;
TH1D* Phip_445MeVCM3_Perp_PosHel;
TH1D* Phip_455MeVCM3_Perp_PosHel;
TH1D* Phip_465MeVCM3_Perp_PosHel;
TH1D* Phip_475MeVCM3_Perp_PosHel;
TH1D* Phip_485MeVCM3_Perp_PosHel;
TH1D* Phip_495MeVCM3_Perp_PosHel;
TH1D* Phip_505MeVCM3_Perp_PosHel;
TH1D* Phip_515MeVCM3_Perp_PosHel;
TH1D* Phip_525MeVCM3_Perp_PosHel;
TH1D* Phip_535MeVCM3_Perp_PosHel;
TH1D* Phip_545MeVCM3_Perp_PosHel;
TH1D* Phip_555MeVCM3_Perp_PosHel;
TH1D* Phip_565MeVCM3_Perp_PosHel;
TH1D* Phip_575MeVCM3_Perp_PosHel;
TH1D* Phip_585MeVCM3_Perp_PosHel;
TH1D* Phip_595MeVCM3_Perp_PosHel;
TH1D* Phip_605MeVCM3_Perp_PosHel;
TH1D* Phip_615MeVCM3_Perp_PosHel;

TH1D* Phip_425MeVCM4_Perp_PosHel;
TH1D* Phip_435MeVCM4_Perp_PosHel;
TH1D* Phip_445MeVCM4_Perp_PosHel;
TH1D* Phip_455MeVCM4_Perp_PosHel;
TH1D* Phip_465MeVCM4_Perp_PosHel;
TH1D* Phip_475MeVCM4_Perp_PosHel;
TH1D* Phip_485MeVCM4_Perp_PosHel;
TH1D* Phip_495MeVCM4_Perp_PosHel;
TH1D* Phip_505MeVCM4_Perp_PosHel;
TH1D* Phip_515MeVCM4_Perp_PosHel;
TH1D* Phip_525MeVCM4_Perp_PosHel;
TH1D* Phip_535MeVCM4_Perp_PosHel;
TH1D* Phip_545MeVCM4_Perp_PosHel;
TH1D* Phip_555MeVCM4_Perp_PosHel;
TH1D* Phip_565MeVCM4_Perp_PosHel;
TH1D* Phip_575MeVCM4_Perp_PosHel;
TH1D* Phip_585MeVCM4_Perp_PosHel;
TH1D* Phip_595MeVCM4_Perp_PosHel;
TH1D* Phip_605MeVCM4_Perp_PosHel;
TH1D* Phip_615MeVCM4_Perp_PosHel;

TH1D* Phip_425MeVCM5_Perp_PosHel;
TH1D* Phip_435MeVCM5_Perp_PosHel;
TH1D* Phip_445MeVCM5_Perp_PosHel;
TH1D* Phip_455MeVCM5_Perp_PosHel;
TH1D* Phip_465MeVCM5_Perp_PosHel;
TH1D* Phip_475MeVCM5_Perp_PosHel;
TH1D* Phip_485MeVCM5_Perp_PosHel;
TH1D* Phip_495MeVCM5_Perp_PosHel;
TH1D* Phip_505MeVCM5_Perp_PosHel;
TH1D* Phip_515MeVCM5_Perp_PosHel;
TH1D* Phip_525MeVCM5_Perp_PosHel;
TH1D* Phip_535MeVCM5_Perp_PosHel;
TH1D* Phip_545MeVCM5_Perp_PosHel;
TH1D* Phip_555MeVCM5_Perp_PosHel;
TH1D* Phip_565MeVCM5_Perp_PosHel;
TH1D* Phip_575MeVCM5_Perp_PosHel;
TH1D* Phip_585MeVCM5_Perp_PosHel;
TH1D* Phip_595MeVCM5_Perp_PosHel;
TH1D* Phip_605MeVCM5_Perp_PosHel;
TH1D* Phip_615MeVCM5_Perp_PosHel;

TH1D* Phip_425MeVCM6_Perp_PosHel;
TH1D* Phip_435MeVCM6_Perp_PosHel;
TH1D* Phip_445MeVCM6_Perp_PosHel;
TH1D* Phip_455MeVCM6_Perp_PosHel;
TH1D* Phip_465MeVCM6_Perp_PosHel;
TH1D* Phip_475MeVCM6_Perp_PosHel;
TH1D* Phip_485MeVCM6_Perp_PosHel;
TH1D* Phip_495MeVCM6_Perp_PosHel;
TH1D* Phip_505MeVCM6_Perp_PosHel;
TH1D* Phip_515MeVCM6_Perp_PosHel;
TH1D* Phip_525MeVCM6_Perp_PosHel;
TH1D* Phip_535MeVCM6_Perp_PosHel;
TH1D* Phip_545MeVCM6_Perp_PosHel;
TH1D* Phip_555MeVCM6_Perp_PosHel;
TH1D* Phip_565MeVCM6_Perp_PosHel;
TH1D* Phip_575MeVCM6_Perp_PosHel;
TH1D* Phip_585MeVCM6_Perp_PosHel;
TH1D* Phip_595MeVCM6_Perp_PosHel;
TH1D* Phip_605MeVCM6_Perp_PosHel;
TH1D* Phip_615MeVCM6_Perp_PosHel;

TH1D* Phip_425MeVCM7_Perp_PosHel;
TH1D* Phip_435MeVCM7_Perp_PosHel;
TH1D* Phip_445MeVCM7_Perp_PosHel;
TH1D* Phip_455MeVCM7_Perp_PosHel;
TH1D* Phip_465MeVCM7_Perp_PosHel;
TH1D* Phip_475MeVCM7_Perp_PosHel;
TH1D* Phip_485MeVCM7_Perp_PosHel;
TH1D* Phip_495MeVCM7_Perp_PosHel;
TH1D* Phip_505MeVCM7_Perp_PosHel;
TH1D* Phip_515MeVCM7_Perp_PosHel;
TH1D* Phip_525MeVCM7_Perp_PosHel;
TH1D* Phip_535MeVCM7_Perp_PosHel;
TH1D* Phip_545MeVCM7_Perp_PosHel;
TH1D* Phip_555MeVCM7_Perp_PosHel;
TH1D* Phip_565MeVCM7_Perp_PosHel;
TH1D* Phip_575MeVCM7_Perp_PosHel;
TH1D* Phip_585MeVCM7_Perp_PosHel;
TH1D* Phip_595MeVCM7_Perp_PosHel;
TH1D* Phip_605MeVCM7_Perp_PosHel;
TH1D* Phip_615MeVCM7_Perp_PosHel;

TH1D* Phip_425MeVCM8_Perp_PosHel;
TH1D* Phip_435MeVCM8_Perp_PosHel;
TH1D* Phip_445MeVCM8_Perp_PosHel;
TH1D* Phip_455MeVCM8_Perp_PosHel;
TH1D* Phip_465MeVCM8_Perp_PosHel;
TH1D* Phip_475MeVCM8_Perp_PosHel;
TH1D* Phip_485MeVCM8_Perp_PosHel;
TH1D* Phip_495MeVCM8_Perp_PosHel;
TH1D* Phip_505MeVCM8_Perp_PosHel;
TH1D* Phip_515MeVCM8_Perp_PosHel;
TH1D* Phip_525MeVCM8_Perp_PosHel;
TH1D* Phip_535MeVCM8_Perp_PosHel;
TH1D* Phip_545MeVCM8_Perp_PosHel;
TH1D* Phip_555MeVCM8_Perp_PosHel;
TH1D* Phip_565MeVCM8_Perp_PosHel;
TH1D* Phip_575MeVCM8_Perp_PosHel;
TH1D* Phip_585MeVCM8_Perp_PosHel;
TH1D* Phip_595MeVCM8_Perp_PosHel;
TH1D* Phip_605MeVCM8_Perp_PosHel;
TH1D* Phip_615MeVCM8_Perp_PosHel;

TH1D* Phip_425MeVCM9_Perp_PosHel;
TH1D* Phip_435MeVCM9_Perp_PosHel;
TH1D* Phip_445MeVCM9_Perp_PosHel;
TH1D* Phip_455MeVCM9_Perp_PosHel;
TH1D* Phip_465MeVCM9_Perp_PosHel;
TH1D* Phip_475MeVCM9_Perp_PosHel;
TH1D* Phip_485MeVCM9_Perp_PosHel;
TH1D* Phip_495MeVCM9_Perp_PosHel;
TH1D* Phip_505MeVCM9_Perp_PosHel;
TH1D* Phip_515MeVCM9_Perp_PosHel;
TH1D* Phip_525MeVCM9_Perp_PosHel;
TH1D* Phip_535MeVCM9_Perp_PosHel;
TH1D* Phip_545MeVCM9_Perp_PosHel;
TH1D* Phip_555MeVCM9_Perp_PosHel;
TH1D* Phip_565MeVCM9_Perp_PosHel;
TH1D* Phip_575MeVCM9_Perp_PosHel;
TH1D* Phip_585MeVCM9_Perp_PosHel;
TH1D* Phip_595MeVCM9_Perp_PosHel;
TH1D* Phip_605MeVCM9_Perp_PosHel;
TH1D* Phip_615MeVCM9_Perp_PosHel;

TH1D* Phip_425MeVCM10_Perp_PosHel;
TH1D* Phip_435MeVCM10_Perp_PosHel;
TH1D* Phip_445MeVCM10_Perp_PosHel;
TH1D* Phip_455MeVCM10_Perp_PosHel;
TH1D* Phip_465MeVCM10_Perp_PosHel;
TH1D* Phip_475MeVCM10_Perp_PosHel;
TH1D* Phip_485MeVCM10_Perp_PosHel;
TH1D* Phip_495MeVCM10_Perp_PosHel;
TH1D* Phip_505MeVCM10_Perp_PosHel;
TH1D* Phip_515MeVCM10_Perp_PosHel;
TH1D* Phip_525MeVCM10_Perp_PosHel;
TH1D* Phip_535MeVCM10_Perp_PosHel;
TH1D* Phip_545MeVCM10_Perp_PosHel;
TH1D* Phip_555MeVCM10_Perp_PosHel;
TH1D* Phip_565MeVCM10_Perp_PosHel;
TH1D* Phip_575MeVCM10_Perp_PosHel;
TH1D* Phip_585MeVCM10_Perp_PosHel;
TH1D* Phip_595MeVCM10_Perp_PosHel;
TH1D* Phip_605MeVCM10_Perp_PosHel;
TH1D* Phip_615MeVCM10_Perp_PosHel;

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
