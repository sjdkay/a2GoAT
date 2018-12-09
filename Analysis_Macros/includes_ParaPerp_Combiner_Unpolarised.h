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

TH1D* PhiDet_Para;
TH1D* PhiRec_Para;
TH1D* Theta_Scattered_Para;
TH1D* Phi_Scattered_Para;
TH1D* MMpEpCorrected_Para;
TH1D* ZpDist_Para;
TH1D* ThetanDist_Para;
TH1D* ThetanCMDist_Para;
TH2D* E_dE_Para;
TH2D* ThetaScPhiSc_Para;
TH2D* EdEMWPCp_Para;
TH2D* EdEMWPCn_Para;
TH1D* ClosestApproach_Para;
TH1D* POCAr_Para;
TH1D* ScatterVertexZ_Para;
TH2D* ScatterVertexZr_Para;
TH2D* ScatterVertexXY_Para;
TH3D* ScatterVertex_Para;
TH2D* POCArPhiSc_Para;
TH1D* ThetapCorrDiff_Para;
TH1D* PhipCorrDiff_Para;
TH1D* ThetaDiff_Para;
TH1D* PhiDiff_Para;
TH2D* PhiDiffThetaDiff_Para;
TH2D* PhiScEg_Para;
TH2D* PhiScEp_Para;
TH2D* PhiScThetan_Para;
TH2D* EMWPCnPhiSc_Para;

TH1D* PhiSc320_Para;
TH1D* PhiSc360_Para;
TH1D* PhiSc400_Para;
TH1D* PhiSc440_Para;
TH1D* PhiSc480_Para;
TH1D* PhiSc520_Para;
TH1D* PhiSc560_Para;
TH1D* PhiSc600_Para;
TH1D* PhiSc640_Para;
TH1D* PhiSc680_Para;

TH1D* MMp200300_Para;
TH1D* MMp300400_Para;
TH1D* MMp400500_Para;
TH1D* MMp500600_Para;
TH1D* MMp600700_Para;
TH1D* MMp700800_Para;
TH1D* MMp800900_Para;

TH1D* Phi_Scattered_265MeV_NegHelCM1_Para;
TH1D* Phi_Scattered_335MeV_NegHelCM1_Para;
TH1D* Phi_Scattered_405MeV_NegHelCM1_Para;
TH1D* Phi_Scattered_475MeV_NegHelCM1_Para;
TH1D* Phi_Scattered_545MeV_NegHelCM1_Para;
TH1D* Phi_Scattered_615MeV_NegHelCM1_Para;
TH1D* Phi_Scattered_685MeV_NegHelCM1_Para;

TH1D* Phi_Scattered_265MeV_NegHelCM2_Para;
TH1D* Phi_Scattered_335MeV_NegHelCM2_Para;
TH1D* Phi_Scattered_405MeV_NegHelCM2_Para;
TH1D* Phi_Scattered_475MeV_NegHelCM2_Para;
TH1D* Phi_Scattered_545MeV_NegHelCM2_Para;
TH1D* Phi_Scattered_615MeV_NegHelCM2_Para;
TH1D* Phi_Scattered_685MeV_NegHelCM2_Para;

TH1D* Phi_Scattered_265MeV_NegHelCM3_Para;
TH1D* Phi_Scattered_335MeV_NegHelCM3_Para;
TH1D* Phi_Scattered_405MeV_NegHelCM3_Para;
TH1D* Phi_Scattered_475MeV_NegHelCM3_Para;
TH1D* Phi_Scattered_545MeV_NegHelCM3_Para;
TH1D* Phi_Scattered_615MeV_NegHelCM3_Para;
TH1D* Phi_Scattered_685MeV_NegHelCM3_Para;

TH1D* Phi_Scattered_265MeV_NegHelCM4_Para;
TH1D* Phi_Scattered_335MeV_NegHelCM4_Para;
TH1D* Phi_Scattered_405MeV_NegHelCM4_Para;
TH1D* Phi_Scattered_475MeV_NegHelCM4_Para;
TH1D* Phi_Scattered_545MeV_NegHelCM4_Para;
TH1D* Phi_Scattered_615MeV_NegHelCM4_Para;
TH1D* Phi_Scattered_685MeV_NegHelCM4_Para;

TH1D* Phi_Scattered_265MeV_NegHelCM5_Para;
TH1D* Phi_Scattered_335MeV_NegHelCM5_Para;
TH1D* Phi_Scattered_405MeV_NegHelCM5_Para;
TH1D* Phi_Scattered_475MeV_NegHelCM5_Para;
TH1D* Phi_Scattered_545MeV_NegHelCM5_Para;
TH1D* Phi_Scattered_615MeV_NegHelCM5_Para;
TH1D* Phi_Scattered_685MeV_NegHelCM5_Para;

TH1D* Phi_Scattered_265MeV_NegHelCM6_Para;
TH1D* Phi_Scattered_335MeV_NegHelCM6_Para;
TH1D* Phi_Scattered_405MeV_NegHelCM6_Para;
TH1D* Phi_Scattered_475MeV_NegHelCM6_Para;
TH1D* Phi_Scattered_545MeV_NegHelCM6_Para;
TH1D* Phi_Scattered_615MeV_NegHelCM6_Para;
TH1D* Phi_Scattered_685MeV_NegHelCM6_Para;

TH1D* Phi_Scattered_265MeV_NegHelCM7_Para;
TH1D* Phi_Scattered_335MeV_NegHelCM7_Para;
TH1D* Phi_Scattered_405MeV_NegHelCM7_Para;
TH1D* Phi_Scattered_475MeV_NegHelCM7_Para;
TH1D* Phi_Scattered_545MeV_NegHelCM7_Para;
TH1D* Phi_Scattered_615MeV_NegHelCM7_Para;
TH1D* Phi_Scattered_685MeV_NegHelCM7_Para;

TH1D* Phi_Scattered_265MeV_NegHelCM8_Para;
TH1D* Phi_Scattered_335MeV_NegHelCM8_Para;
TH1D* Phi_Scattered_405MeV_NegHelCM8_Para;
TH1D* Phi_Scattered_475MeV_NegHelCM8_Para;
TH1D* Phi_Scattered_545MeV_NegHelCM8_Para;
TH1D* Phi_Scattered_615MeV_NegHelCM8_Para;
TH1D* Phi_Scattered_685MeV_NegHelCM8_Para;

TH1D* Phi_Scattered_265MeV_TH1D* PosHelCM1_Para;
TH1D* Phi_Scattered_335MeV_TH1D* PosHelCM1_Para;
TH1D* Phi_Scattered_405MeV_TH1D* PosHelCM1_Para;
TH1D* Phi_Scattered_475MeV_TH1D* PosHelCM1_Para;
TH1D* Phi_Scattered_545MeV_TH1D* PosHelCM1_Para;
TH1D* Phi_Scattered_615MeV_TH1D* PosHelCM1_Para;
TH1D* Phi_Scattered_685MeV_TH1D* PosHelCM1_Para;

TH1D* Phi_Scattered_265MeV_TH1D* PosHelCM2_Para;
TH1D* Phi_Scattered_335MeV_TH1D* PosHelCM2_Para;
TH1D* Phi_Scattered_405MeV_TH1D* PosHelCM2_Para;
TH1D* Phi_Scattered_475MeV_TH1D* PosHelCM2_Para;
TH1D* Phi_Scattered_545MeV_TH1D* PosHelCM2_Para;
TH1D* Phi_Scattered_615MeV_TH1D* PosHelCM2_Para;
TH1D* Phi_Scattered_685MeV_TH1D* PosHelCM2_Para;

TH1D* Phi_Scattered_265MeV_TH1D* PosHelCM3_Para;
TH1D* Phi_Scattered_335MeV_TH1D* PosHelCM3_Para;
TH1D* Phi_Scattered_405MeV_TH1D* PosHelCM3_Para;
TH1D* Phi_Scattered_475MeV_TH1D* PosHelCM3_Para;
TH1D* Phi_Scattered_545MeV_TH1D* PosHelCM3_Para;
TH1D* Phi_Scattered_615MeV_TH1D* PosHelCM3_Para;
TH1D* Phi_Scattered_685MeV_TH1D* PosHelCM3_Para;

TH1D* Phi_Scattered_265MeV_TH1D* PosHelCM4_Para;
TH1D* Phi_Scattered_335MeV_TH1D* PosHelCM4_Para;
TH1D* Phi_Scattered_405MeV_TH1D* PosHelCM4_Para;
TH1D* Phi_Scattered_475MeV_TH1D* PosHelCM4_Para;
TH1D* Phi_Scattered_545MeV_TH1D* PosHelCM4_Para;
TH1D* Phi_Scattered_615MeV_TH1D* PosHelCM4_Para;
TH1D* Phi_Scattered_685MeV_TH1D* PosHelCM4_Para;

TH1D* Phi_Scattered_265MeV_TH1D* PosHelCM5_Para;
TH1D* Phi_Scattered_335MeV_TH1D* PosHelCM5_Para;
TH1D* Phi_Scattered_405MeV_TH1D* PosHelCM5_Para;
TH1D* Phi_Scattered_475MeV_TH1D* PosHelCM5_Para;
TH1D* Phi_Scattered_545MeV_TH1D* PosHelCM5_Para;
TH1D* Phi_Scattered_615MeV_TH1D* PosHelCM5_Para;
TH1D* Phi_Scattered_685MeV_TH1D* PosHelCM5_Para;

TH1D* Phi_Scattered_265MeV_TH1D* PosHelCM6_Para;
TH1D* Phi_Scattered_335MeV_TH1D* PosHelCM6_Para;
TH1D* Phi_Scattered_405MeV_TH1D* PosHelCM6_Para;
TH1D* Phi_Scattered_475MeV_TH1D* PosHelCM6_Para;
TH1D* Phi_Scattered_545MeV_TH1D* PosHelCM6_Para;
TH1D* Phi_Scattered_615MeV_TH1D* PosHelCM6_Para;
TH1D* Phi_Scattered_685MeV_TH1D* PosHelCM6_Para;

TH1D* Phi_Scattered_265MeV_TH1D* PosHelCM7_Para;
TH1D* Phi_Scattered_335MeV_TH1D* PosHelCM7_Para;
TH1D* Phi_Scattered_405MeV_TH1D* PosHelCM7_Para;
TH1D* Phi_Scattered_475MeV_TH1D* PosHelCM7_Para;
TH1D* Phi_Scattered_545MeV_TH1D* PosHelCM7_Para;
TH1D* Phi_Scattered_615MeV_TH1D* PosHelCM7_Para;
TH1D* Phi_Scattered_685MeV_TH1D* PosHelCM7_Para;

TH1D* Phi_Scattered_265MeV_TH1D* PosHelCM8_Para;
TH1D* Phi_Scattered_335MeV_TH1D* PosHelCM8_Para;
TH1D* Phi_Scattered_405MeV_TH1D* PosHelCM8_Para;
TH1D* Phi_Scattered_475MeV_TH1D* PosHelCM8_Para;
TH1D* Phi_Scattered_545MeV_TH1D* PosHelCM8_Para;
TH1D* Phi_Scattered_615MeV_TH1D* PosHelCM8_Para;
TH1D* Phi_Scattered_685MeV_TH1D* PosHelCM8_Para;

TH1D* time_Perp;
TH1D* time_cut_Perp;
TH1D* Eg_Perp;

TH1D* PhiDet_Perp;
TH1D* PhiRec_Perp;
TH1D* Theta_Scattered_Perp;
TH1D* Phi_Scattered_Perp;
TH1D* MMpEpCorrected_Perp;
TH1D* ZpDist_Perp;
TH1D* ThetanDist_Perp;
TH1D* ThetanCMDist_Perp;
TH2D* E_dE_Perp;
TH2D* ThetaScPhiSc_Perp;
TH2D* EdEMWPCp_Perp;
TH2D* EdEMWPCn_Perp;
TH1D* ClosestApproach_Perp;
TH1D* POCAr_Perp;
TH1D* ScatterVertexZ_Perp;
TH2D* ScatterVertexZr_Perp;
TH2D* ScatterVertexXY_Perp;
TH3D* ScatterVertex_Perp;
TH2D* POCArPhiSc_Perp;
TH1D* ThetapCorrDiff_Perp;
TH1D* PhipCorrDiff_Perp;
TH1D* ThetaDiff_Perp;
TH1D* PhiDiff_Perp;
TH2D* PhiDiffThetaDiff_Perp;
TH2D* PhiScEg_Perp;
TH2D* PhiScEp_Perp;
TH2D* PhiScThetan_Perp;
TH2D* EMWPCnPhiSc_Perp;

TH1D* PhiSc320_Perp;
TH1D* PhiSc360_Perp;
TH1D* PhiSc400_Perp;
TH1D* PhiSc440_Perp;
TH1D* PhiSc480_Perp;
TH1D* PhiSc520_Perp;
TH1D* PhiSc560_Perp;
TH1D* PhiSc600_Perp;
TH1D* PhiSc640_Perp;
TH1D* PhiSc680_Perp;

TH1D* MMp200300_Perp;
TH1D* MMp300400_Perp;
TH1D* MMp400500_Perp;
TH1D* MMp500600_Perp;
TH1D* MMp600700_Perp;
TH1D* MMp700800_Perp;
TH1D* MMp800900_Perp;

TH1D* Phi_Scattered_265MeV_NegHelCM1_Perp;
TH1D* Phi_Scattered_335MeV_NegHelCM1_Perp;
TH1D* Phi_Scattered_405MeV_NegHelCM1_Perp;
TH1D* Phi_Scattered_475MeV_NegHelCM1_Perp;
TH1D* Phi_Scattered_545MeV_NegHelCM1_Perp;
TH1D* Phi_Scattered_615MeV_NegHelCM1_Perp;
TH1D* Phi_Scattered_685MeV_NegHelCM1_Perp;

TH1D* Phi_Scattered_265MeV_NegHelCM2_Perp;
TH1D* Phi_Scattered_335MeV_NegHelCM2_Perp;
TH1D* Phi_Scattered_405MeV_NegHelCM2_Perp;
TH1D* Phi_Scattered_475MeV_NegHelCM2_Perp;
TH1D* Phi_Scattered_545MeV_NegHelCM2_Perp;
TH1D* Phi_Scattered_615MeV_NegHelCM2_Perp;
TH1D* Phi_Scattered_685MeV_NegHelCM2_Perp;

TH1D* Phi_Scattered_265MeV_NegHelCM3_Perp;
TH1D* Phi_Scattered_335MeV_NegHelCM3_Perp;
TH1D* Phi_Scattered_405MeV_NegHelCM3_Perp;
TH1D* Phi_Scattered_475MeV_NegHelCM3_Perp;
TH1D* Phi_Scattered_545MeV_NegHelCM3_Perp;
TH1D* Phi_Scattered_615MeV_NegHelCM3_Perp;
TH1D* Phi_Scattered_685MeV_NegHelCM3_Perp;

TH1D* Phi_Scattered_265MeV_NegHelCM4_Perp;
TH1D* Phi_Scattered_335MeV_NegHelCM4_Perp;
TH1D* Phi_Scattered_405MeV_NegHelCM4_Perp;
TH1D* Phi_Scattered_475MeV_NegHelCM4_Perp;
TH1D* Phi_Scattered_545MeV_NegHelCM4_Perp;
TH1D* Phi_Scattered_615MeV_NegHelCM4_Perp;
TH1D* Phi_Scattered_685MeV_NegHelCM4_Perp;

TH1D* Phi_Scattered_265MeV_NegHelCM5_Perp;
TH1D* Phi_Scattered_335MeV_NegHelCM5_Perp;
TH1D* Phi_Scattered_405MeV_NegHelCM5_Perp;
TH1D* Phi_Scattered_475MeV_NegHelCM5_Perp;
TH1D* Phi_Scattered_545MeV_NegHelCM5_Perp;
TH1D* Phi_Scattered_615MeV_NegHelCM5_Perp;
TH1D* Phi_Scattered_685MeV_NegHelCM5_Perp;

TH1D* Phi_Scattered_265MeV_NegHelCM6_Perp;
TH1D* Phi_Scattered_335MeV_NegHelCM6_Perp;
TH1D* Phi_Scattered_405MeV_NegHelCM6_Perp;
TH1D* Phi_Scattered_475MeV_NegHelCM6_Perp;
TH1D* Phi_Scattered_545MeV_NegHelCM6_Perp;
TH1D* Phi_Scattered_615MeV_NegHelCM6_Perp;
TH1D* Phi_Scattered_685MeV_NegHelCM6_Perp;

TH1D* Phi_Scattered_265MeV_NegHelCM7_Perp;
TH1D* Phi_Scattered_335MeV_NegHelCM7_Perp;
TH1D* Phi_Scattered_405MeV_NegHelCM7_Perp;
TH1D* Phi_Scattered_475MeV_NegHelCM7_Perp;
TH1D* Phi_Scattered_545MeV_NegHelCM7_Perp;
TH1D* Phi_Scattered_615MeV_NegHelCM7_Perp;
TH1D* Phi_Scattered_685MeV_NegHelCM7_Perp;

TH1D* Phi_Scattered_265MeV_NegHelCM8_Perp;
TH1D* Phi_Scattered_335MeV_NegHelCM8_Perp;
TH1D* Phi_Scattered_405MeV_NegHelCM8_Perp;
TH1D* Phi_Scattered_475MeV_NegHelCM8_Perp;
TH1D* Phi_Scattered_545MeV_NegHelCM8_Perp;
TH1D* Phi_Scattered_615MeV_NegHelCM8_Perp;
TH1D* Phi_Scattered_685MeV_NegHelCM8_Perp;

TH1D* Phi_Scattered_265MeV_TH1D* PosHelCM1_Perp;
TH1D* Phi_Scattered_335MeV_TH1D* PosHelCM1_Perp;
TH1D* Phi_Scattered_405MeV_TH1D* PosHelCM1_Perp;
TH1D* Phi_Scattered_475MeV_TH1D* PosHelCM1_Perp;
TH1D* Phi_Scattered_545MeV_TH1D* PosHelCM1_Perp;
TH1D* Phi_Scattered_615MeV_TH1D* PosHelCM1_Perp;
TH1D* Phi_Scattered_685MeV_TH1D* PosHelCM1_Perp;

TH1D* Phi_Scattered_265MeV_TH1D* PosHelCM2_Perp;
TH1D* Phi_Scattered_335MeV_TH1D* PosHelCM2_Perp;
TH1D* Phi_Scattered_405MeV_TH1D* PosHelCM2_Perp;
TH1D* Phi_Scattered_475MeV_TH1D* PosHelCM2_Perp;
TH1D* Phi_Scattered_545MeV_TH1D* PosHelCM2_Perp;
TH1D* Phi_Scattered_615MeV_TH1D* PosHelCM2_Perp;
TH1D* Phi_Scattered_685MeV_TH1D* PosHelCM2_Perp;

TH1D* Phi_Scattered_265MeV_TH1D* PosHelCM3_Perp;
TH1D* Phi_Scattered_335MeV_TH1D* PosHelCM3_Perp;
TH1D* Phi_Scattered_405MeV_TH1D* PosHelCM3_Perp;
TH1D* Phi_Scattered_475MeV_TH1D* PosHelCM3_Perp;
TH1D* Phi_Scattered_545MeV_TH1D* PosHelCM3_Perp;
TH1D* Phi_Scattered_615MeV_TH1D* PosHelCM3_Perp;
TH1D* Phi_Scattered_685MeV_TH1D* PosHelCM3_Perp;

TH1D* Phi_Scattered_265MeV_TH1D* PosHelCM4_Perp;
TH1D* Phi_Scattered_335MeV_TH1D* PosHelCM4_Perp;
TH1D* Phi_Scattered_405MeV_TH1D* PosHelCM4_Perp;
TH1D* Phi_Scattered_475MeV_TH1D* PosHelCM4_Perp;
TH1D* Phi_Scattered_545MeV_TH1D* PosHelCM4_Perp;
TH1D* Phi_Scattered_615MeV_TH1D* PosHelCM4_Perp;
TH1D* Phi_Scattered_685MeV_TH1D* PosHelCM4_Perp;

TH1D* Phi_Scattered_265MeV_TH1D* PosHelCM5_Perp;
TH1D* Phi_Scattered_335MeV_TH1D* PosHelCM5_Perp;
TH1D* Phi_Scattered_405MeV_TH1D* PosHelCM5_Perp;
TH1D* Phi_Scattered_475MeV_TH1D* PosHelCM5_Perp;
TH1D* Phi_Scattered_545MeV_TH1D* PosHelCM5_Perp;
TH1D* Phi_Scattered_615MeV_TH1D* PosHelCM5_Perp;
TH1D* Phi_Scattered_685MeV_TH1D* PosHelCM5_Perp;

TH1D* Phi_Scattered_265MeV_TH1D* PosHelCM6_Perp;
TH1D* Phi_Scattered_335MeV_TH1D* PosHelCM6_Perp;
TH1D* Phi_Scattered_405MeV_TH1D* PosHelCM6_Perp;
TH1D* Phi_Scattered_475MeV_TH1D* PosHelCM6_Perp;
TH1D* Phi_Scattered_545MeV_TH1D* PosHelCM6_Perp;
TH1D* Phi_Scattered_615MeV_TH1D* PosHelCM6_Perp;
TH1D* Phi_Scattered_685MeV_TH1D* PosHelCM6_Perp;

TH1D* Phi_Scattered_265MeV_TH1D* PosHelCM7_Perp;
TH1D* Phi_Scattered_335MeV_TH1D* PosHelCM7_Perp;
TH1D* Phi_Scattered_405MeV_TH1D* PosHelCM7_Perp;
TH1D* Phi_Scattered_475MeV_TH1D* PosHelCM7_Perp;
TH1D* Phi_Scattered_545MeV_TH1D* PosHelCM7_Perp;
TH1D* Phi_Scattered_615MeV_TH1D* PosHelCM7_Perp;
TH1D* Phi_Scattered_685MeV_TH1D* PosHelCM7_Perp;

TH1D* Phi_Scattered_265MeV_TH1D* PosHelCM8_Perp;
TH1D* Phi_Scattered_335MeV_TH1D* PosHelCM8_Perp;
TH1D* Phi_Scattered_405MeV_TH1D* PosHelCM8_Perp;
TH1D* Phi_Scattered_475MeV_TH1D* PosHelCM8_Perp;
TH1D* Phi_Scattered_545MeV_TH1D* PosHelCM8_Perp;
TH1D* Phi_Scattered_615MeV_TH1D* PosHelCM8_Perp;
TH1D* Phi_Scattered_685MeV_TH1D* PosHelCM8_Perp;

TH1D* time;
TH1D* time_cut;
TH1D* Eg;

TH1D* PhiDet;
TH1D* PhiRec;
TH1D* Theta_Scattered;
TH1D* Phi_Scattered;
TH1D* MMpEpCorrected;
TH1D* ZpDist;
TH1D* ThetanDist;
TH1D* ThetanCMDist;
TH2D* E_dE;
TH2D* ThetaScPhiSc;
TH2D* EdEMWPCp;
TH2D* EdEMWPCn;
TH1D* ClosestApproach;
TH1D* POCAr;
TH1D* ScatterVertexZ;
TH2D* ScatterVertexZr;
TH2D* ScatterVertexXY;
TH3D* ScatterVertex;
TH2D* POCArPhiSc;
TH1D* ThetapCorrDiff;
TH1D* PhipCorrDiff;
TH1D* ThetaDiff;
TH1D* PhiDiff;
TH2D* PhiDiffThetaDiff;
TH2D* PhiScEg;
TH2D* PhiScEp;
TH2D* PhiScThetan;
TH2D* EMWPCnPhiSc;

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

TH1D* MMp200300;
TH1D* MMp300400;
TH1D* MMp400500;
TH1D* MMp500600;
TH1D* MMp600700;
TH1D* MMp700800;
TH1D* MMp800900;

TH1D* Phi_Scattered_265MeV_NegHelCM1;
TH1D* Phi_Scattered_335MeV_NegHelCM1;
TH1D* Phi_Scattered_405MeV_NegHelCM1;
TH1D* Phi_Scattered_475MeV_NegHelCM1;
TH1D* Phi_Scattered_545MeV_NegHelCM1;
TH1D* Phi_Scattered_615MeV_NegHelCM1;
TH1D* Phi_Scattered_685MeV_NegHelCM1;

TH1D* Phi_Scattered_265MeV_NegHelCM2;
TH1D* Phi_Scattered_335MeV_NegHelCM2;
TH1D* Phi_Scattered_405MeV_NegHelCM2;
TH1D* Phi_Scattered_475MeV_NegHelCM2;
TH1D* Phi_Scattered_545MeV_NegHelCM2;
TH1D* Phi_Scattered_615MeV_NegHelCM2;
TH1D* Phi_Scattered_685MeV_NegHelCM2;

TH1D* Phi_Scattered_265MeV_NegHelCM3;
TH1D* Phi_Scattered_335MeV_NegHelCM3;
TH1D* Phi_Scattered_405MeV_NegHelCM3;
TH1D* Phi_Scattered_475MeV_NegHelCM3;
TH1D* Phi_Scattered_545MeV_NegHelCM3;
TH1D* Phi_Scattered_615MeV_NegHelCM3;
TH1D* Phi_Scattered_685MeV_NegHelCM3;

TH1D* Phi_Scattered_265MeV_NegHelCM4;
TH1D* Phi_Scattered_335MeV_NegHelCM4;
TH1D* Phi_Scattered_405MeV_NegHelCM4;
TH1D* Phi_Scattered_475MeV_NegHelCM4;
TH1D* Phi_Scattered_545MeV_NegHelCM4;
TH1D* Phi_Scattered_615MeV_NegHelCM4;
TH1D* Phi_Scattered_685MeV_NegHelCM4;

TH1D* Phi_Scattered_265MeV_NegHelCM5;
TH1D* Phi_Scattered_335MeV_NegHelCM5;
TH1D* Phi_Scattered_405MeV_NegHelCM5;
TH1D* Phi_Scattered_475MeV_NegHelCM5;
TH1D* Phi_Scattered_545MeV_NegHelCM5;
TH1D* Phi_Scattered_615MeV_NegHelCM5;
TH1D* Phi_Scattered_685MeV_NegHelCM5;

TH1D* Phi_Scattered_265MeV_NegHelCM6;
TH1D* Phi_Scattered_335MeV_NegHelCM6;
TH1D* Phi_Scattered_405MeV_NegHelCM6;
TH1D* Phi_Scattered_475MeV_NegHelCM6;
TH1D* Phi_Scattered_545MeV_NegHelCM6;
TH1D* Phi_Scattered_615MeV_NegHelCM6;
TH1D* Phi_Scattered_685MeV_NegHelCM6;

TH1D* Phi_Scattered_265MeV_NegHelCM7;
TH1D* Phi_Scattered_335MeV_NegHelCM7;
TH1D* Phi_Scattered_405MeV_NegHelCM7;
TH1D* Phi_Scattered_475MeV_NegHelCM7;
TH1D* Phi_Scattered_545MeV_NegHelCM7;
TH1D* Phi_Scattered_615MeV_NegHelCM7;
TH1D* Phi_Scattered_685MeV_NegHelCM7;

TH1D* Phi_Scattered_265MeV_NegHelCM8;
TH1D* Phi_Scattered_335MeV_NegHelCM8;
TH1D* Phi_Scattered_405MeV_NegHelCM8;
TH1D* Phi_Scattered_475MeV_NegHelCM8;
TH1D* Phi_Scattered_545MeV_NegHelCM8;
TH1D* Phi_Scattered_615MeV_NegHelCM8;
TH1D* Phi_Scattered_685MeV_NegHelCM8;

TH1D* Phi_Scattered_265MeV_TH1D* PosHelCM1;
TH1D* Phi_Scattered_335MeV_TH1D* PosHelCM1;
TH1D* Phi_Scattered_405MeV_TH1D* PosHelCM1;
TH1D* Phi_Scattered_475MeV_TH1D* PosHelCM1;
TH1D* Phi_Scattered_545MeV_TH1D* PosHelCM1;
TH1D* Phi_Scattered_615MeV_TH1D* PosHelCM1;
TH1D* Phi_Scattered_685MeV_TH1D* PosHelCM1;

TH1D* Phi_Scattered_265MeV_TH1D* PosHelCM2;
TH1D* Phi_Scattered_335MeV_TH1D* PosHelCM2;
TH1D* Phi_Scattered_405MeV_TH1D* PosHelCM2;
TH1D* Phi_Scattered_475MeV_TH1D* PosHelCM2;
TH1D* Phi_Scattered_545MeV_TH1D* PosHelCM2;
TH1D* Phi_Scattered_615MeV_TH1D* PosHelCM2;
TH1D* Phi_Scattered_685MeV_TH1D* PosHelCM2;

TH1D* Phi_Scattered_265MeV_TH1D* PosHelCM3;
TH1D* Phi_Scattered_335MeV_TH1D* PosHelCM3;
TH1D* Phi_Scattered_405MeV_TH1D* PosHelCM3;
TH1D* Phi_Scattered_475MeV_TH1D* PosHelCM3;
TH1D* Phi_Scattered_545MeV_TH1D* PosHelCM3;
TH1D* Phi_Scattered_615MeV_TH1D* PosHelCM3;
TH1D* Phi_Scattered_685MeV_TH1D* PosHelCM3;

TH1D* Phi_Scattered_265MeV_TH1D* PosHelCM4;
TH1D* Phi_Scattered_335MeV_TH1D* PosHelCM4;
TH1D* Phi_Scattered_405MeV_TH1D* PosHelCM4;
TH1D* Phi_Scattered_475MeV_TH1D* PosHelCM4;
TH1D* Phi_Scattered_545MeV_TH1D* PosHelCM4;
TH1D* Phi_Scattered_615MeV_TH1D* PosHelCM4;
TH1D* Phi_Scattered_685MeV_TH1D* PosHelCM4;

TH1D* Phi_Scattered_265MeV_TH1D* PosHelCM5;
TH1D* Phi_Scattered_335MeV_TH1D* PosHelCM5;
TH1D* Phi_Scattered_405MeV_TH1D* PosHelCM5;
TH1D* Phi_Scattered_475MeV_TH1D* PosHelCM5;
TH1D* Phi_Scattered_545MeV_TH1D* PosHelCM5;
TH1D* Phi_Scattered_615MeV_TH1D* PosHelCM5;
TH1D* Phi_Scattered_685MeV_TH1D* PosHelCM5;

TH1D* Phi_Scattered_265MeV_TH1D* PosHelCM6;
TH1D* Phi_Scattered_335MeV_TH1D* PosHelCM6;
TH1D* Phi_Scattered_405MeV_TH1D* PosHelCM6;
TH1D* Phi_Scattered_475MeV_TH1D* PosHelCM6;
TH1D* Phi_Scattered_545MeV_TH1D* PosHelCM6;
TH1D* Phi_Scattered_615MeV_TH1D* PosHelCM6;
TH1D* Phi_Scattered_685MeV_TH1D* PosHelCM6;

TH1D* Phi_Scattered_265MeV_TH1D* PosHelCM7;
TH1D* Phi_Scattered_335MeV_TH1D* PosHelCM7;
TH1D* Phi_Scattered_405MeV_TH1D* PosHelCM7;
TH1D* Phi_Scattered_475MeV_TH1D* PosHelCM7;
TH1D* Phi_Scattered_545MeV_TH1D* PosHelCM7;
TH1D* Phi_Scattered_615MeV_TH1D* PosHelCM7;
TH1D* Phi_Scattered_685MeV_TH1D* PosHelCM7;

TH1D* Phi_Scattered_265MeV_TH1D* PosHelCM8;
TH1D* Phi_Scattered_335MeV_TH1D* PosHelCM8;
TH1D* Phi_Scattered_405MeV_TH1D* PosHelCM8;
TH1D* Phi_Scattered_475MeV_TH1D* PosHelCM8;
TH1D* Phi_Scattered_545MeV_TH1D* PosHelCM8;
TH1D* Phi_Scattered_615MeV_TH1D* PosHelCM8;
TH1D* Phi_Scattered_685MeV_TH1D* PosHelCM8;

Double_t NPara;
Double_t NPerp;
Double_t ScaleFactor;
Double_t ScaleFactorErr;

