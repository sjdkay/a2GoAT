#include "./includes_ParaPerp_Combiner_Unpolarised.h"

void ParaPerp_Combiner_Unpolarised() {

    TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/Para/Physics_Total_Para_28_30_11_17.root"); // Open latest Para file

//    TH1D* time_Para = (TH1D*)f->Get("time")->Clone();
//    time_Para->SetName("time_Para");
//    TH1D* time_cut_Para = (TH1D*)f->Get("time_cut")->Clone();
//    time_cut_Para->SetName("time_cut_Para");
//    TH1D* Eg_Para = (TH1D*)f->Get("Eg")->Clone();
//    Eg_Para->SetName("Eg_Para");
//    TH1D* PhiDet_Para = (TH1D*)f->Get("PhiDet")->Clone();
//    PhiDet_Para->SetName("PhiDet_Para");
//    TH1D* PhiRec_Para = (TH1D*)f->Get("PhiRec")->Clone();
//    PhiRec_Para->SetName("PhiRec_Para");
//    TH1D* Theta_Scattered_Para = (TH1D*)f->Get("Theta_Scattered")->Clone();
//    Theta_Scattered_Para->SetName("Theta_Scattered_Para");
//    TH1D* Phi_Scattered_Para = (TH1D*)f->Get("Phi_Scattered")->Clone();
//    Phi_Scattered_Para->SetName("Phi_Scattered_Para");
//    TH1D* MMpEpCorrected_Para = (TH1D*)f->Get("MMpEpCorrected")->Clone();
//    MMpEpCorrected_Para->SetName("MMpEpCorrected_Para");
//    TH1D* ZpDist_Para = (TH1D*)f->Get("ZpDist")->Clone();
//    ZpDist_Para->SetName("ZpDist_Para");
//    TH1D* ThetanDist_Para = (TH1D*)f->Get("ThetanDist")->Clone();
//    ThetanDist_Para->SetName("ThetanDist_Para");
//    TH1D* ThetanCMDist_Para = (TH1D*)f->Get("ThetanCMDist")->Clone();
//    ThetanCMDist_Para->SetName("ThetanCMDist_Para");
//    TH2D* E_dE_Para = (TH2D*)f->Get("E_dE")->Clone();
//    E_dE_Para->SetName("E_dE_Para");
//    TH2D* ThetaScPhiSc_Para = (TH2D*)f->Get("ThetaScPhiSc")->Clone();
//    ThetaScPhiSc_Para->SetName("ThetaScPhiSc_Para");
//    TH2D* EdEMWPCp_Para = (TH2D*)f->Get("EdEMWPCp")->Clone();
//    EdEMWPCp_Para->SetName("EdEMWPCp_Para");
//    TH2D* EdEMWPCn_Para = (TH2D*)f->Get("EdEMWPCn")->Clone();
//    EdEMWPCn_Para->SetName("EdEMWPCn_Para");
//    TH1D* ClosestApproach_Para = (TH1D*)f->Get("ClosestApproach")->Clone();
//    ClosestApproach_Para->SetName("ClosestApproach_Para");
//    TH1D* POCAr_Para = (TH1D*)f->Get("POCAr")->Clone();
//    POCAr_Para->SetName("POCAr_Para");
//    TH1D* ScatterVertexZ_Para = (TH1D*)f->Get("ScatterVertexZ")->Clone();
//    ScatterVertexZ_Para->SetName("ScatterVertexZ_Para");
//    TH2D* ScatterVertexZr_Para = (TH2D*)f->Get("ScatterVertexZr")->Clone();
//    ScatterVertexZr_Para->SetName("ScatterVertexZr_Para");
//    TH2D* ScatterVertexXY_Para = (TH2D*)f->Get("ScatterVertexXY")->Clone();
//    ScatterVertexXY_Para->SetName("ScatterVertexXY_Para");
//    TH3D* ScatterVertex_Para = (TH3D*)f->Get("ScatterVertex")->Clone();
//    ScatterVertex_Para->SetName("ScatterVertex_Para");
//    TH2D* POCArPhiSc_Para = (TH2D*)f->Get("POCArPhiSc")->Clone();
//    POCArPhiSc_Para->SetName("POCArPhiSc_Para");
//    TH1D* ThetapCorrDiff_Para = (TH1D*)f->Get("ThetapCorrDiff")->Clone();
//    ThetapCorrDiff_Para->SetName("ThetapCorrDiff_Para");
//    TH1D* PhipCorrDiff_Para = (TH1D*)f->Get("PhipCorrDiff")->Clone();
//    PhipCorrDiff_Para->SetName("PhipCorrDiff_Para");
//    TH1D* ThetaDiff_Para = (TH1D*)f->Get("ThetaDiff")->Clone();
//    ThetaDiff_Para->SetName("ThetaDiff_Para");
//    TH1D* PhiDiff_Para = (TH1D*)f->Get("PhiDiff")->Clone();
//    PhiDiff_Para->SetName("PhiDiff_Para");
//    TH2D* PhiDiffThetaDiff_Para = (TH2D*)f->Get("PhiDiffThetaDiff")->Clone();
//    PhiDiffThetaDiff_Para->SetName("PhiDiffThetaDiff_Para");
//    TH2D* PhiScEg_Para = (TH2D*)f->Get("PhiScEg")->Clone();
//    PhiScEg_Para->SetName("PhiScEg_Para");
//    TH2D* PhiScEp_Para = (TH2D*)f->Get("PhiScEp")->Clone();
//    PhiScEp_Para->SetName("PhiScEp_Para");
//    TH2D* PhiScThetan_Para = (TH2D*)f->Get("PhiScThetan")->Clone();
//    PhiScThetan_Para->SetName("PhiScThetan_Para");
//    TH2D* EMWPCnPhiSc_Para = (TH2D*)f->Get("EMWPCnPhiSc")->Clone();
//    EMWPCnPhiSc_Para->SetName("EMWPCnPhiSc_Para");
//
//    TH1D* MMp200300_Para = (TH1D*)f->Get("MMp200300")->Clone();
//    MMp200300_Para->SetName("MMp200300_Para");
//    TH1D* MMp300400_Para = (TH1D*)f->Get("MMp300400")->Clone();
//    MMp300400_Para->SetName("MMp300400_Para");
//    TH1D* MMp400500_Para = (TH1D*)f->Get("MMp400500")->Clone();
//    MMp400500_Para->SetName("MMp400500_Para");
//    TH1D* MMp500600_Para = (TH1D*)f->Get("MMp500600")->Clone();
//    MMp500600_Para->SetName("MMp500600_Para");
//    TH1D* MMp600700_Para = (TH1D*)f->Get("MMp600700")->Clone();
//    MMp600700_Para->SetName("MMp600700_Para");
//    TH1D* MMp700800_Para = (TH1D*)f->Get("MMp700800")->Clone();
//    MMp700800_Para->SetName("MMp700800_Para");
//    TH1D* MMp800900_Para = (TH1D*)f->Get("MMp800900")->Clone();
//    MMp800900_Para->SetName("MMp800900_Para");

    TH1D* PhiSc320_Para = (TH1D*)f->Get("PhiSc320")->Clone();
    PhiSc320_Para->SetName("PhiSc320_Para");
    TH1D* PhiSc360_Para = (TH1D*)f->Get("PhiSc360")->Clone();
    PhiSc360_Para->SetName("PhiSc360_Para");
    TH1D* PhiSc400_Para = (TH1D*)f->Get("PhiSc400")->Clone();
    PhiSc400_Para->SetName("PhiSc400_Para");
    TH1D* PhiSc440_Para = (TH1D*)f->Get("PhiSc440")->Clone();
    PhiSc440_Para->SetName("PhiSc440_Para");
    TH1D* PhiSc480_Para = (TH1D*)f->Get("PhiSc480")->Clone();
    PhiSc480_Para->SetName("PhiSc480_Para");
    TH1D* PhiSc520_Para = (TH1D*)f->Get("PhiSc520")->Clone();
    PhiSc520_Para->SetName("PhiSc520_Para");
    TH1D* PhiSc560_Para = (TH1D*)f->Get("PhiSc560")->Clone();
    PhiSc560_Para->SetName("PhiSc560_Para");
    TH1D* PhiSc600_Para = (TH1D*)f->Get("PhiSc600")->Clone();
    PhiSc600_Para->SetName("PhiSc600_Para");
    TH1D* PhiSc640_Para = (TH1D*)f->Get("PhiSc640")->Clone();
    PhiSc640_Para->SetName("PhiSc640_Para");
    TH1D* PhiSc680_Para = (TH1D*)f->Get("PhiSc680")->Clone();
    PhiSc680_Para->SetName("PhiSc680_Para");
//
//    Phi_Scattered_265MeV_NegHelCM1_Para = (TH1D*)f->Get("Phi_Scattered_265MeV_NegHelCM1")->Clone();
//    Phi_Scattered_265MeV_NegHelCM1_Para->SetName("Phi_Scattered_265MeV_NegHelCM1_Para");
//    Phi_Scattered_335MeV_NegHelCM1_Para = (TH1D*)f->Get("Phi_Scattered_335MeV_NegHelCM1")->Clone();
//    Phi_Scattered_335MeV_NegHelCM1_Para->SetName("Phi_Scattered_335MeV_NegHelCM1_Para");
//    Phi_Scattered_405MeV_NegHelCM1_Para = (TH1D*)f->Get("Phi_Scattered_405MeV_NegHelCM1")->Clone();
//    Phi_Scattered_405MeV_NegHelCM1_Para->SetName("Phi_Scattered_405MeV_NegHelCM1_Para");
//    Phi_Scattered_475MeV_NegHelCM1_Para = (TH1D*)f->Get("Phi_Scattered_475MeV_NegHelCM1")->Clone();
//    Phi_Scattered_475MeV_NegHelCM1_Para->SetName("Phi_Scattered_475MeV_NegHelCM1_Para");
//    Phi_Scattered_545MeV_NegHelCM1_Para = (TH1D*)f->Get("Phi_Scattered_545MeV_NegHelCM1")->Clone();
//    Phi_Scattered_545MeV_NegHelCM1_Para->SetName("Phi_Scattered_545MeV_NegHelCM1_Para");
//    Phi_Scattered_615MeV_NegHelCM1_Para = (TH1D*)f->Get("Phi_Scattered_615MeV_NegHelCM1")->Clone();
//    Phi_Scattered_615MeV_NegHelCM1_Para->SetName("Phi_Scattered_615MeV_NegHelCM1_Para");
//    Phi_Scattered_685MeV_NegHelCM1_Para = (TH1D*)f->Get("Phi_Scattered_685MeV_NegHelCM1")->Clone();
//    Phi_Scattered_685MeV_NegHelCM1_Para->SetName("Phi_Scattered_685MeV_NegHelCM1_Para");
//
//    Phi_Scattered_265MeV_NegHelCM2_Para = (TH1D*)f->Get("Phi_Scattered_265MeV_NegHelCM2")->Clone();
//    Phi_Scattered_265MeV_NegHelCM2_Para->SetName("Phi_Scattered_265MeV_NegHelCM2_Para");
//    Phi_Scattered_335MeV_NegHelCM2_Para = (TH1D*)f->Get("Phi_Scattered_335MeV_NegHelCM2")->Clone();
//    Phi_Scattered_335MeV_NegHelCM2_Para->SetName("Phi_Scattered_335MeV_NegHelCM2_Para");
//    Phi_Scattered_405MeV_NegHelCM2_Para = (TH1D*)f->Get("Phi_Scattered_405MeV_NegHelCM2")->Clone();
//    Phi_Scattered_405MeV_NegHelCM2_Para->SetName("Phi_Scattered_405MeV_NegHelCM2_Para");
//    Phi_Scattered_475MeV_NegHelCM2_Para = (TH1D*)f->Get("Phi_Scattered_475MeV_NegHelCM2")->Clone();
//    Phi_Scattered_475MeV_NegHelCM2_Para->SetName("Phi_Scattered_475MeV_NegHelCM2_Para");
//    Phi_Scattered_545MeV_NegHelCM2_Para = (TH1D*)f->Get("Phi_Scattered_545MeV_NegHelCM2")->Clone();
//    Phi_Scattered_545MeV_NegHelCM2_Para->SetName("Phi_Scattered_545MeV_NegHelCM2_Para");
//    Phi_Scattered_615MeV_NegHelCM2_Para = (TH1D*)f->Get("Phi_Scattered_615MeV_NegHelCM2")->Clone();
//    Phi_Scattered_615MeV_NegHelCM2_Para->SetName("Phi_Scattered_615MeV_NegHelCM2_Para");
//    Phi_Scattered_685MeV_NegHelCM2_Para = (TH1D*)f->Get("Phi_Scattered_685MeV_NegHelCM2")->Clone();
//    Phi_Scattered_685MeV_NegHelCM2_Para->SetName("Phi_Scattered_685MeV_NegHelCM2_Para");
//
//    Phi_Scattered_265MeV_NegHelCM3_Para = (TH1D*)f->Get("Phi_Scattered_265MeV_NegHelCM3")->Clone();
//    Phi_Scattered_265MeV_NegHelCM3_Para->SetName("Phi_Scattered_265MeV_NegHelCM3_Para");
//    Phi_Scattered_335MeV_NegHelCM3_Para = (TH1D*)f->Get("Phi_Scattered_335MeV_NegHelCM3")->Clone();
//    Phi_Scattered_335MeV_NegHelCM3_Para->SetName("Phi_Scattered_335MeV_NegHelCM3_Para");
//    Phi_Scattered_405MeV_NegHelCM3_Para = (TH1D*)f->Get("Phi_Scattered_405MeV_NegHelCM3")->Clone();
//    Phi_Scattered_405MeV_NegHelCM3_Para->SetName("Phi_Scattered_405MeV_NegHelCM3_Para");
//    Phi_Scattered_475MeV_NegHelCM3_Para = (TH1D*)f->Get("Phi_Scattered_475MeV_NegHelCM3")->Clone();
//    Phi_Scattered_475MeV_NegHelCM3_Para->SetName("Phi_Scattered_475MeV_NegHelCM3_Para");
//    Phi_Scattered_545MeV_NegHelCM3_Para = (TH1D*)f->Get("Phi_Scattered_545MeV_NegHelCM3")->Clone();
//    Phi_Scattered_545MeV_NegHelCM3_Para->SetName("Phi_Scattered_545MeV_NegHelCM3_Para");
//    Phi_Scattered_615MeV_NegHelCM3_Para = (TH1D*)f->Get("Phi_Scattered_615MeV_NegHelCM3")->Clone();
//    Phi_Scattered_615MeV_NegHelCM3_Para->SetName("Phi_Scattered_615MeV_NegHelCM3_Para");
//    Phi_Scattered_685MeV_NegHelCM3_Para = (TH1D*)f->Get("Phi_Scattered_685MeV_NegHelCM3")->Clone();
//    Phi_Scattered_685MeV_NegHelCM3_Para->SetName("Phi_Scattered_685MeV_NegHelCM3_Para");
//
//    Phi_Scattered_265MeV_NegHelCM4_Para = (TH1D*)f->Get("Phi_Scattered_265MeV_NegHelCM4")->Clone();
//    Phi_Scattered_265MeV_NegHelCM4_Para->SetName("Phi_Scattered_265MeV_NegHelCM4_Para");
//    Phi_Scattered_335MeV_NegHelCM4_Para = (TH1D*)f->Get("Phi_Scattered_335MeV_NegHelCM4")->Clone();
//    Phi_Scattered_335MeV_NegHelCM4_Para->SetName("Phi_Scattered_335MeV_NegHelCM4_Para");
//    Phi_Scattered_405MeV_NegHelCM4_Para = (TH1D*)f->Get("Phi_Scattered_405MeV_NegHelCM4")->Clone();
//    Phi_Scattered_405MeV_NegHelCM4_Para->SetName("Phi_Scattered_405MeV_NegHelCM4_Para");
//    Phi_Scattered_475MeV_NegHelCM4_Para = (TH1D*)f->Get("Phi_Scattered_475MeV_NegHelCM4")->Clone();
//    Phi_Scattered_475MeV_NegHelCM4_Para->SetName("Phi_Scattered_475MeV_NegHelCM4_Para");
//    Phi_Scattered_545MeV_NegHelCM4_Para = (TH1D*)f->Get("Phi_Scattered_545MeV_NegHelCM4")->Clone();
//    Phi_Scattered_545MeV_NegHelCM4_Para->SetName("Phi_Scattered_545MeV_NegHelCM4_Para");
//    Phi_Scattered_615MeV_NegHelCM4_Para = (TH1D*)f->Get("Phi_Scattered_615MeV_NegHelCM4")->Clone();
//    Phi_Scattered_615MeV_NegHelCM4_Para->SetName("Phi_Scattered_615MeV_NegHelCM4_Para");
//    Phi_Scattered_685MeV_NegHelCM4_Para = (TH1D*)f->Get("Phi_Scattered_685MeV_NegHelCM4")->Clone();
//    Phi_Scattered_685MeV_NegHelCM4_Para->SetName("Phi_Scattered_685MeV_NegHelCM4_Para");
//
//    Phi_Scattered_265MeV_NegHelCM5_Para = (TH1D*)f->Get("Phi_Scattered_265MeV_NegHelCM5")->Clone();
//    Phi_Scattered_265MeV_NegHelCM5_Para->SetName("Phi_Scattered_265MeV_NegHelCM5_Para");
//    Phi_Scattered_335MeV_NegHelCM5_Para = (TH1D*)f->Get("Phi_Scattered_335MeV_NegHelCM5")->Clone();
//    Phi_Scattered_335MeV_NegHelCM5_Para->SetName("Phi_Scattered_335MeV_NegHelCM5_Para");
//    Phi_Scattered_405MeV_NegHelCM5_Para = (TH1D*)f->Get("Phi_Scattered_405MeV_NegHelCM5")->Clone();
//    Phi_Scattered_405MeV_NegHelCM5_Para->SetName("Phi_Scattered_405MeV_NegHelCM5_Para");
//    Phi_Scattered_475MeV_NegHelCM5_Para = (TH1D*)f->Get("Phi_Scattered_475MeV_NegHelCM5")->Clone();
//    Phi_Scattered_475MeV_NegHelCM5_Para->SetName("Phi_Scattered_475MeV_NegHelCM5_Para");
//    Phi_Scattered_545MeV_NegHelCM5_Para = (TH1D*)f->Get("Phi_Scattered_545MeV_NegHelCM5")->Clone();
//    Phi_Scattered_545MeV_NegHelCM5_Para->SetName("Phi_Scattered_545MeV_NegHelCM5_Para");
//    Phi_Scattered_615MeV_NegHelCM5_Para = (TH1D*)f->Get("Phi_Scattered_615MeV_NegHelCM5")->Clone();
//    Phi_Scattered_615MeV_NegHelCM5_Para->SetName("Phi_Scattered_615MeV_NegHelCM5_Para");
//    Phi_Scattered_685MeV_NegHelCM5_Para = (TH1D*)f->Get("Phi_Scattered_685MeV_NegHelCM5")->Clone();
//    Phi_Scattered_685MeV_NegHelCM5_Para->SetName("Phi_Scattered_685MeV_NegHelCM5_Para");
//
//    Phi_Scattered_265MeV_NegHelCM6_Para = (TH1D*)f->Get("Phi_Scattered_265MeV_NegHelCM6")->Clone();
//    Phi_Scattered_265MeV_NegHelCM6_Para->SetName("Phi_Scattered_265MeV_NegHelCM6_Para");
//    Phi_Scattered_335MeV_NegHelCM6_Para = (TH1D*)f->Get("Phi_Scattered_335MeV_NegHelCM6")->Clone();
//    Phi_Scattered_335MeV_NegHelCM6_Para->SetName("Phi_Scattered_335MeV_NegHelCM6_Para");
//    Phi_Scattered_405MeV_NegHelCM6_Para = (TH1D*)f->Get("Phi_Scattered_405MeV_NegHelCM6")->Clone();
//    Phi_Scattered_405MeV_NegHelCM6_Para->SetName("Phi_Scattered_405MeV_NegHelCM6_Para");
//    Phi_Scattered_475MeV_NegHelCM6_Para = (TH1D*)f->Get("Phi_Scattered_475MeV_NegHelCM6")->Clone();
//    Phi_Scattered_475MeV_NegHelCM6_Para->SetName("Phi_Scattered_475MeV_NegHelCM6_Para");
//    Phi_Scattered_545MeV_NegHelCM6_Para = (TH1D*)f->Get("Phi_Scattered_545MeV_NegHelCM6")->Clone();
//    Phi_Scattered_545MeV_NegHelCM6_Para->SetName("Phi_Scattered_545MeV_NegHelCM6_Para");
//    Phi_Scattered_615MeV_NegHelCM6_Para = (TH1D*)f->Get("Phi_Scattered_615MeV_NegHelCM6")->Clone();
//    Phi_Scattered_615MeV_NegHelCM6_Para->SetName("Phi_Scattered_615MeV_NegHelCM6_Para");
//    Phi_Scattered_685MeV_NegHelCM6_Para = (TH1D*)f->Get("Phi_Scattered_685MeV_NegHelCM6")->Clone();
//    Phi_Scattered_685MeV_NegHelCM6_Para->SetName("Phi_Scattered_685MeV_NegHelCM6_Para");
//
//    Phi_Scattered_265MeV_NegHelCM7_Para = (TH1D*)f->Get("Phi_Scattered_265MeV_NegHelCM7")->Clone();
//    Phi_Scattered_265MeV_NegHelCM7_Para->SetName("Phi_Scattered_265MeV_NegHelCM7_Para");
//    Phi_Scattered_335MeV_NegHelCM7_Para = (TH1D*)f->Get("Phi_Scattered_335MeV_NegHelCM7")->Clone();
//    Phi_Scattered_335MeV_NegHelCM7_Para->SetName("Phi_Scattered_335MeV_NegHelCM7_Para");
//    Phi_Scattered_405MeV_NegHelCM7_Para = (TH1D*)f->Get("Phi_Scattered_405MeV_NegHelCM7")->Clone();
//    Phi_Scattered_405MeV_NegHelCM7_Para->SetName("Phi_Scattered_405MeV_NegHelCM7_Para");
//    Phi_Scattered_475MeV_NegHelCM7_Para = (TH1D*)f->Get("Phi_Scattered_475MeV_NegHelCM7")->Clone();
//    Phi_Scattered_475MeV_NegHelCM7_Para->SetName("Phi_Scattered_475MeV_NegHelCM7_Para");
//    Phi_Scattered_545MeV_NegHelCM7_Para = (TH1D*)f->Get("Phi_Scattered_545MeV_NegHelCM7")->Clone();
//    Phi_Scattered_545MeV_NegHelCM7_Para->SetName("Phi_Scattered_545MeV_NegHelCM7_Para");
//    Phi_Scattered_615MeV_NegHelCM7_Para = (TH1D*)f->Get("Phi_Scattered_615MeV_NegHelCM7")->Clone();
//    Phi_Scattered_615MeV_NegHelCM7_Para->SetName("Phi_Scattered_615MeV_NegHelCM7_Para");
//    Phi_Scattered_685MeV_NegHelCM7_Para = (TH1D*)f->Get("Phi_Scattered_685MeV_NegHelCM7")->Clone();
//    Phi_Scattered_685MeV_NegHelCM7_Para->SetName("Phi_Scattered_685MeV_NegHelCM7_Para");
//
//    Phi_Scattered_265MeV_NegHelCM8_Para = (TH1D*)f->Get("Phi_Scattered_265MeV_NegHelCM8")->Clone();
//    Phi_Scattered_265MeV_NegHelCM8_Para->SetName("Phi_Scattered_265MeV_NegHelCM8_Para");
//    Phi_Scattered_335MeV_NegHelCM8_Para = (TH1D*)f->Get("Phi_Scattered_335MeV_NegHelCM8")->Clone();
//    Phi_Scattered_335MeV_NegHelCM8_Para->SetName("Phi_Scattered_335MeV_NegHelCM8_Para");
//    Phi_Scattered_405MeV_NegHelCM8_Para = (TH1D*)f->Get("Phi_Scattered_405MeV_NegHelCM8")->Clone();
//    Phi_Scattered_405MeV_NegHelCM8_Para->SetName("Phi_Scattered_405MeV_NegHelCM8_Para");
//    Phi_Scattered_475MeV_NegHelCM8_Para = (TH1D*)f->Get("Phi_Scattered_475MeV_NegHelCM8")->Clone();
//    Phi_Scattered_475MeV_NegHelCM8_Para->SetName("Phi_Scattered_475MeV_NegHelCM8_Para");
//    Phi_Scattered_545MeV_NegHelCM8_Para = (TH1D*)f->Get("Phi_Scattered_545MeV_NegHelCM8")->Clone();
//    Phi_Scattered_545MeV_NegHelCM8_Para->SetName("Phi_Scattered_545MeV_NegHelCM8_Para");
//    Phi_Scattered_615MeV_NegHelCM8_Para = (TH1D*)f->Get("Phi_Scattered_615MeV_NegHelCM8")->Clone();
//    Phi_Scattered_615MeV_NegHelCM8_Para->SetName("Phi_Scattered_615MeV_NegHelCM8_Para");
//    Phi_Scattered_685MeV_NegHelCM8_Para = (TH1D*)f->Get("Phi_Scattered_685MeV_NegHelCM8")->Clone();
//    Phi_Scattered_685MeV_NegHelCM8_Para->SetName("Phi_Scattered_685MeV_NegHelCM8_Para");
//
//    Phi_Scattered_265MeV_PosHelCM1_Para = (TH1D*)f->Get("Phi_Scattered_265MeV_PosHelCM1")->Clone();
//    Phi_Scattered_265MeV_PosHelCM1_Para->SetName("Phi_Scattered_265MeV_PosHelCM1_Para");
//    Phi_Scattered_335MeV_PosHelCM1_Para = (TH1D*)f->Get("Phi_Scattered_335MeV_PosHelCM1")->Clone();
//    Phi_Scattered_335MeV_PosHelCM1_Para->SetName("Phi_Scattered_335MeV_PosHelCM1_Para");
//    Phi_Scattered_405MeV_PosHelCM1_Para = (TH1D*)f->Get("Phi_Scattered_405MeV_PosHelCM1")->Clone();
//    Phi_Scattered_405MeV_PosHelCM1_Para->SetName("Phi_Scattered_405MeV_PosHelCM1_Para");
//    Phi_Scattered_475MeV_PosHelCM1_Para = (TH1D*)f->Get("Phi_Scattered_475MeV_PosHelCM1")->Clone();
//    Phi_Scattered_475MeV_PosHelCM1_Para->SetName("Phi_Scattered_475MeV_PosHelCM1_Para");
//    Phi_Scattered_545MeV_PosHelCM1_Para = (TH1D*)f->Get("Phi_Scattered_545MeV_PosHelCM1")->Clone();
//    Phi_Scattered_545MeV_PosHelCM1_Para->SetName("Phi_Scattered_545MeV_PosHelCM1_Para");
//    Phi_Scattered_615MeV_PosHelCM1_Para = (TH1D*)f->Get("Phi_Scattered_615MeV_PosHelCM1")->Clone();
//    Phi_Scattered_615MeV_PosHelCM1_Para->SetName("Phi_Scattered_615MeV_PosHelCM1_Para");
//    Phi_Scattered_685MeV_PosHelCM1_Para = (TH1D*)f->Get("Phi_Scattered_685MeV_PosHelCM1")->Clone();
//    Phi_Scattered_685MeV_PosHelCM1_Para->SetName("Phi_Scattered_685MeV_PosHelCM1_Para");
//
//    Phi_Scattered_265MeV_PosHelCM2_Para = (TH1D*)f->Get("Phi_Scattered_265MeV_PosHelCM2")->Clone();
//    Phi_Scattered_265MeV_PosHelCM2_Para->SetName("Phi_Scattered_265MeV_PosHelCM2_Para");
//    Phi_Scattered_335MeV_PosHelCM2_Para = (TH1D*)f->Get("Phi_Scattered_335MeV_PosHelCM2")->Clone();
//    Phi_Scattered_335MeV_PosHelCM2_Para->SetName("Phi_Scattered_335MeV_PosHelCM2_Para");
//    Phi_Scattered_405MeV_PosHelCM2_Para = (TH1D*)f->Get("Phi_Scattered_405MeV_PosHelCM2")->Clone();
//    Phi_Scattered_405MeV_PosHelCM2_Para->SetName("Phi_Scattered_405MeV_PosHelCM2_Para");
//    Phi_Scattered_475MeV_PosHelCM2_Para = (TH1D*)f->Get("Phi_Scattered_475MeV_PosHelCM2")->Clone();
//    Phi_Scattered_475MeV_PosHelCM2_Para->SetName("Phi_Scattered_475MeV_PosHelCM2_Para");
//    Phi_Scattered_545MeV_PosHelCM2_Para = (TH1D*)f->Get("Phi_Scattered_545MeV_PosHelCM2")->Clone();
//    Phi_Scattered_545MeV_PosHelCM2_Para->SetName("Phi_Scattered_545MeV_PosHelCM2_Para");
//    Phi_Scattered_615MeV_PosHelCM2_Para = (TH1D*)f->Get("Phi_Scattered_615MeV_PosHelCM2")->Clone();
//    Phi_Scattered_615MeV_PosHelCM2_Para->SetName("Phi_Scattered_615MeV_PosHelCM2_Para");
//    Phi_Scattered_685MeV_PosHelCM2_Para = (TH1D*)f->Get("Phi_Scattered_685MeV_PosHelCM2")->Clone();
//    Phi_Scattered_685MeV_PosHelCM2_Para->SetName("Phi_Scattered_685MeV_PosHelCM2_Para");
//
//    Phi_Scattered_265MeV_PosHelCM3_Para = (TH1D*)f->Get("Phi_Scattered_265MeV_PosHelCM3")->Clone();
//    Phi_Scattered_265MeV_PosHelCM3_Para->SetName("Phi_Scattered_265MeV_PosHelCM3_Para");
//    Phi_Scattered_335MeV_PosHelCM3_Para = (TH1D*)f->Get("Phi_Scattered_335MeV_PosHelCM3")->Clone();
//    Phi_Scattered_335MeV_PosHelCM3_Para->SetName("Phi_Scattered_335MeV_PosHelCM3_Para");
//    Phi_Scattered_405MeV_PosHelCM3_Para = (TH1D*)f->Get("Phi_Scattered_405MeV_PosHelCM3")->Clone();
//    Phi_Scattered_405MeV_PosHelCM3_Para->SetName("Phi_Scattered_405MeV_PosHelCM3_Para");
//    Phi_Scattered_475MeV_PosHelCM3_Para = (TH1D*)f->Get("Phi_Scattered_475MeV_PosHelCM3")->Clone();
//    Phi_Scattered_475MeV_PosHelCM3_Para->SetName("Phi_Scattered_475MeV_PosHelCM3_Para");
//    Phi_Scattered_545MeV_PosHelCM3_Para = (TH1D*)f->Get("Phi_Scattered_545MeV_PosHelCM3")->Clone();
//    Phi_Scattered_545MeV_PosHelCM3_Para->SetName("Phi_Scattered_545MeV_PosHelCM3_Para");
//    Phi_Scattered_615MeV_PosHelCM3_Para = (TH1D*)f->Get("Phi_Scattered_615MeV_PosHelCM3")->Clone();
//    Phi_Scattered_615MeV_PosHelCM3_Para->SetName("Phi_Scattered_615MeV_PosHelCM3_Para");
//    Phi_Scattered_685MeV_PosHelCM3_Para = (TH1D*)f->Get("Phi_Scattered_685MeV_PosHelCM3")->Clone();
//    Phi_Scattered_685MeV_PosHelCM3_Para->SetName("Phi_Scattered_685MeV_PosHelCM3_Para");
//
//    Phi_Scattered_265MeV_PosHelCM4_Para = (TH1D*)f->Get("Phi_Scattered_265MeV_PosHelCM4")->Clone();
//    Phi_Scattered_265MeV_PosHelCM4_Para->SetName("Phi_Scattered_265MeV_PosHelCM4_Para");
//    Phi_Scattered_335MeV_PosHelCM4_Para = (TH1D*)f->Get("Phi_Scattered_335MeV_PosHelCM4")->Clone();
//    Phi_Scattered_335MeV_PosHelCM4_Para->SetName("Phi_Scattered_335MeV_PosHelCM4_Para");
//    Phi_Scattered_405MeV_PosHelCM4_Para = (TH1D*)f->Get("Phi_Scattered_405MeV_PosHelCM4")->Clone();
//    Phi_Scattered_405MeV_PosHelCM4_Para->SetName("Phi_Scattered_405MeV_PosHelCM4_Para");
//    Phi_Scattered_475MeV_PosHelCM4_Para = (TH1D*)f->Get("Phi_Scattered_475MeV_PosHelCM4")->Clone();
//    Phi_Scattered_475MeV_PosHelCM4_Para->SetName("Phi_Scattered_475MeV_PosHelCM4_Para");
//    Phi_Scattered_545MeV_PosHelCM4_Para = (TH1D*)f->Get("Phi_Scattered_545MeV_PosHelCM4")->Clone();
//    Phi_Scattered_545MeV_PosHelCM4_Para->SetName("Phi_Scattered_545MeV_PosHelCM4_Para");
//    Phi_Scattered_615MeV_PosHelCM4_Para = (TH1D*)f->Get("Phi_Scattered_615MeV_PosHelCM4")->Clone();
//    Phi_Scattered_615MeV_PosHelCM4_Para->SetName("Phi_Scattered_615MeV_PosHelCM4_Para");
//    Phi_Scattered_685MeV_PosHelCM4_Para = (TH1D*)f->Get("Phi_Scattered_685MeV_PosHelCM4")->Clone();
//    Phi_Scattered_685MeV_PosHelCM4_Para->SetName("Phi_Scattered_685MeV_PosHelCM4_Para");
//
//    Phi_Scattered_265MeV_PosHelCM5_Para = (TH1D*)f->Get("Phi_Scattered_265MeV_PosHelCM5")->Clone();
//    Phi_Scattered_265MeV_PosHelCM5_Para->SetName("Phi_Scattered_265MeV_PosHelCM5_Para");
//    Phi_Scattered_335MeV_PosHelCM5_Para = (TH1D*)f->Get("Phi_Scattered_335MeV_PosHelCM5")->Clone();
//    Phi_Scattered_335MeV_PosHelCM5_Para->SetName("Phi_Scattered_335MeV_PosHelCM5_Para");
//    Phi_Scattered_405MeV_PosHelCM5_Para = (TH1D*)f->Get("Phi_Scattered_405MeV_PosHelCM5")->Clone();
//    Phi_Scattered_405MeV_PosHelCM5_Para->SetName("Phi_Scattered_405MeV_PosHelCM5_Para");
//    Phi_Scattered_475MeV_PosHelCM5_Para = (TH1D*)f->Get("Phi_Scattered_475MeV_PosHelCM5")->Clone();
//    Phi_Scattered_475MeV_PosHelCM5_Para->SetName("Phi_Scattered_475MeV_PosHelCM5_Para");
//    Phi_Scattered_545MeV_PosHelCM5_Para = (TH1D*)f->Get("Phi_Scattered_545MeV_PosHelCM5")->Clone();
//    Phi_Scattered_545MeV_PosHelCM5_Para->SetName("Phi_Scattered_545MeV_PosHelCM5_Para");
//    Phi_Scattered_615MeV_PosHelCM5_Para = (TH1D*)f->Get("Phi_Scattered_615MeV_PosHelCM5")->Clone();
//    Phi_Scattered_615MeV_PosHelCM5_Para->SetName("Phi_Scattered_615MeV_PosHelCM5_Para");
//    Phi_Scattered_685MeV_PosHelCM5_Para = (TH1D*)f->Get("Phi_Scattered_685MeV_PosHelCM5")->Clone();
//    Phi_Scattered_685MeV_PosHelCM5_Para->SetName("Phi_Scattered_685MeV_PosHelCM5_Para");
//
//    Phi_Scattered_265MeV_PosHelCM6_Para = (TH1D*)f->Get("Phi_Scattered_265MeV_PosHelCM6")->Clone();
//    Phi_Scattered_265MeV_PosHelCM6_Para->SetName("Phi_Scattered_265MeV_PosHelCM6_Para");
//    Phi_Scattered_335MeV_PosHelCM6_Para = (TH1D*)f->Get("Phi_Scattered_335MeV_PosHelCM6")->Clone();
//    Phi_Scattered_335MeV_PosHelCM6_Para->SetName("Phi_Scattered_335MeV_PosHelCM6_Para");
//    Phi_Scattered_405MeV_PosHelCM6_Para = (TH1D*)f->Get("Phi_Scattered_405MeV_PosHelCM6")->Clone();
//    Phi_Scattered_405MeV_PosHelCM6_Para->SetName("Phi_Scattered_405MeV_PosHelCM6_Para");
//    Phi_Scattered_475MeV_PosHelCM6_Para = (TH1D*)f->Get("Phi_Scattered_475MeV_PosHelCM6")->Clone();
//    Phi_Scattered_475MeV_PosHelCM6_Para->SetName("Phi_Scattered_475MeV_PosHelCM6_Para");
//    Phi_Scattered_545MeV_PosHelCM6_Para = (TH1D*)f->Get("Phi_Scattered_545MeV_PosHelCM6")->Clone();
//    Phi_Scattered_545MeV_PosHelCM6_Para->SetName("Phi_Scattered_545MeV_PosHelCM6_Para");
//    Phi_Scattered_615MeV_PosHelCM6_Para = (TH1D*)f->Get("Phi_Scattered_615MeV_PosHelCM6")->Clone();
//    Phi_Scattered_615MeV_PosHelCM6_Para->SetName("Phi_Scattered_615MeV_PosHelCM6_Para");
//    Phi_Scattered_685MeV_PosHelCM6_Para = (TH1D*)f->Get("Phi_Scattered_685MeV_PosHelCM6")->Clone();
//    Phi_Scattered_685MeV_PosHelCM6_Para->SetName("Phi_Scattered_685MeV_PosHelCM6_Para");
//
//    Phi_Scattered_265MeV_PosHelCM7_Para = (TH1D*)f->Get("Phi_Scattered_265MeV_PosHelCM7")->Clone();
//    Phi_Scattered_265MeV_PosHelCM7_Para->SetName("Phi_Scattered_265MeV_PosHelCM7_Para");
//    Phi_Scattered_335MeV_PosHelCM7_Para = (TH1D*)f->Get("Phi_Scattered_335MeV_PosHelCM7")->Clone();
//    Phi_Scattered_335MeV_PosHelCM7_Para->SetName("Phi_Scattered_335MeV_PosHelCM7_Para");
//    Phi_Scattered_405MeV_PosHelCM7_Para = (TH1D*)f->Get("Phi_Scattered_405MeV_PosHelCM7")->Clone();
//    Phi_Scattered_405MeV_PosHelCM7_Para->SetName("Phi_Scattered_405MeV_PosHelCM7_Para");
//    Phi_Scattered_475MeV_PosHelCM7_Para = (TH1D*)f->Get("Phi_Scattered_475MeV_PosHelCM7")->Clone();
//    Phi_Scattered_475MeV_PosHelCM7_Para->SetName("Phi_Scattered_475MeV_PosHelCM7_Para");
//    Phi_Scattered_545MeV_PosHelCM7_Para = (TH1D*)f->Get("Phi_Scattered_545MeV_PosHelCM7")->Clone();
//    Phi_Scattered_545MeV_PosHelCM7_Para->SetName("Phi_Scattered_545MeV_PosHelCM7_Para");
//    Phi_Scattered_615MeV_PosHelCM7_Para = (TH1D*)f->Get("Phi_Scattered_615MeV_PosHelCM7")->Clone();
//    Phi_Scattered_615MeV_PosHelCM7_Para->SetName("Phi_Scattered_615MeV_PosHelCM7_Para");
//    Phi_Scattered_685MeV_PosHelCM7_Para = (TH1D*)f->Get("Phi_Scattered_685MeV_PosHelCM7")->Clone();
//    Phi_Scattered_685MeV_PosHelCM7_Para->SetName("Phi_Scattered_685MeV_PosHelCM7_Para");
//
//    Phi_Scattered_265MeV_PosHelCM8_Para = (TH1D*)f->Get("Phi_Scattered_265MeV_PosHelCM8")->Clone();
//    Phi_Scattered_265MeV_PosHelCM8_Para->SetName("Phi_Scattered_265MeV_PosHelCM8_Para");
//    Phi_Scattered_335MeV_PosHelCM8_Para = (TH1D*)f->Get("Phi_Scattered_335MeV_PosHelCM8")->Clone();
//    Phi_Scattered_335MeV_PosHelCM8_Para->SetName("Phi_Scattered_335MeV_PosHelCM8_Para");
//    Phi_Scattered_405MeV_PosHelCM8_Para = (TH1D*)f->Get("Phi_Scattered_405MeV_PosHelCM8")->Clone();
//    Phi_Scattered_405MeV_PosHelCM8_Para->SetName("Phi_Scattered_405MeV_PosHelCM8_Para");
//    Phi_Scattered_475MeV_PosHelCM8_Para = (TH1D*)f->Get("Phi_Scattered_475MeV_PosHelCM8")->Clone();
//    Phi_Scattered_475MeV_PosHelCM8_Para->SetName("Phi_Scattered_475MeV_PosHelCM8_Para");
//    Phi_Scattered_545MeV_PosHelCM8_Para = (TH1D*)f->Get("Phi_Scattered_545MeV_PosHelCM8")->Clone();
//    Phi_Scattered_545MeV_PosHelCM8_Para->SetName("Phi_Scattered_545MeV_PosHelCM8_Para");
//    Phi_Scattered_615MeV_PosHelCM8_Para = (TH1D*)f->Get("Phi_Scattered_615MeV_PosHelCM8")->Clone();
//    Phi_Scattered_615MeV_PosHelCM8_Para->SetName("Phi_Scattered_615MeV_PosHelCM8_Para");
//    Phi_Scattered_685MeV_PosHelCM8_Para = (TH1D*)f->Get("Phi_Scattered_685MeV_PosHelCM8")->Clone();
//    Phi_Scattered_685MeV_PosHelCM8_Para->SetName("Phi_Scattered_685MeV_PosHelCM8_Para");

    ////////////////////////////////////////////////////
    ////////////////////////////////////////////////////
    ///////////////// PARA DONE ////////////////////////
    ////////////////////////////////////////////////////
    ////////////////////////////////////////////////////

    TFile *f1 = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/Perp/Physics_Total_Perp_28_30_11_17.root"); // Open latest Para file

//    TH1D* time_Perp = (TH1D*)f->Get("time")->Clone();
//    time_Perp->SetName("time_Perp");
//    TH1D* time_cut_Perp = (TH1D*)f->Get("time_cut")->Clone();
//    time_cut_Perp->SetName("time_cut_Perp");
//    TH1D* Eg_Perp = (TH1D*)f->Get("Eg")->Clone();
//    Eg_Perp->SetName("Eg_Perp");
//    TH1D* PhiDet_Perp = (TH1D*)f->Get("PhiDet")->Clone();
//    PhiDet_Perp->SetName("PhiDet_Perp");
//    TH1D* PhiRec_Perp = (TH1D*)f->Get("PhiRec")->Clone();
//    PhiRec_Perp->SetName("PhiRec_Perp");
//    TH1D* Theta_Scattered_Perp = (TH1D*)f->Get("Theta_Scattered")->Clone();
//    Theta_Scattered_Perp->SetName("Theta_Scattered_Perp");
//    TH1D* Phi_Scattered_Perp = (TH1D*)f->Get("Phi_Scattered")->Clone();
//    Phi_Scattered_Perp->SetName("Phi_Scattered_Perp");
//    TH1D* MMpEpCorrected_Perp = (TH1D*)f->Get("MMpEpCorrected")->Clone();
//    MMpEpCorrected_Perp->SetName("MMpEpCorrected_Perp");
//    TH1D* ZpDist_Perp = (TH1D*)f->Get("ZpDist")->Clone();
//    ZpDist_Perp->SetName("ZpDist_Perp");
//    TH1D* ThetanDist_Perp = (TH1D*)f->Get("ThetanDist")->Clone();
//    ThetanDist_Perp->SetName("ThetanDist_Perp");
//    TH1D* ThetanCMDist_Perp = (TH1D*)f->Get("ThetanCMDist")->Clone();
//    ThetanCMDist_Perp->SetName("ThetanCMDist_Perp");
//    TH2D* E_dE_Perp = (TH2D*)f->Get("E_dE")->Clone();
//    E_dE_Perp->SetName("E_dE_Perp");
//    TH2D* ThetaScPhiSc_Perp = (TH2D*)f->Get("ThetaScPhiSc")->Clone();
//    ThetaScPhiSc_Perp->SetName("ThetaScPhiSc_Perp");
//    TH2D* EdEMWPCp_Perp = (TH2D*)f->Get("EdEMWPCp")->Clone();
//    EdEMWPCp_Perp->SetName("EdEMWPCp_Perp");
//    TH2D* EdEMWPCn_Perp = (TH2D*)f->Get("EdEMWPCn")->Clone();
//    EdEMWPCn_Perp->SetName("EdEMWPCn_Perp");
//    TH1D* ClosestApproach_Perp = (TH1D*)f->Get("ClosestApproach")->Clone();
//    ClosestApproach_Perp->SetName("ClosestApproach_Perp");
//    TH1D* POCAr_Perp = (TH1D*)f->Get("POCAr")->Clone();
//    POCAr_Perp->SetName("POCAr_Perp");
//    TH1D* ScatterVertexZ_Perp = (TH1D*)f->Get("ScatterVertexZ")->Clone();
//    ScatterVertexZ_Perp->SetName("ScatterVertexZ_Perp");
//    TH2D* ScatterVertexZr_Perp = (TH2D*)f->Get("ScatterVertexZr")->Clone();
//    ScatterVertexZr_Perp->SetName("ScatterVertexZr_Perp");
//    TH2D* ScatterVertexXY_Perp = (TH2D*)f->Get("ScatterVertexXY")->Clone();
//    ScatterVertexXY_Perp->SetName("ScatterVertexXY_Perp");
//    TH3D* ScatterVertex_Perp = (TH3D*)f->Get("ScatterVertex")->Clone();
//    ScatterVertex_Perp->SetName("ScatterVertex_Perp");
//    TH2D* POCArPhiSc_Perp = (TH2D*)f->Get("POCArPhiSc")->Clone();
//    POCArPhiSc_Perp->SetName("POCArPhiSc_Perp");
//    TH1D* ThetapCorrDiff_Perp = (TH1D*)f->Get("ThetapCorrDiff")->Clone();
//    ThetapCorrDiff_Perp->SetName("ThetapCorrDiff_Perp");
//    TH1D* PhipCorrDiff_Perp = (TH1D*)f->Get("PhipCorrDiff")->Clone();
//    PhipCorrDiff_Perp->SetName("PhipCorrDiff_Perp");
//    TH1D* ThetaDiff_Perp = (TH1D*)f->Get("ThetaDiff")->Clone();
//    ThetaDiff_Perp->SetName("ThetaDiff_Perp");
//    TH1D* PhiDiff_Perp = (TH1D*)f->Get("PhiDiff")->Clone();
//    PhiDiff_Perp->SetName("PhiDiff_Perp");
//    TH2D* PhiDiffThetaDiff_Perp = (TH2D*)f->Get("PhiDiffThetaDiff")->Clone();
//    PhiDiffThetaDiff_Perp->SetName("PhiDiffThetaDiff_Perp");
//    TH2D* PhiScEg_Perp = (TH2D*)f->Get("PhiScEg")->Clone();
//    PhiScEg_Perp->SetName("PhiScEg_Perp");
//    TH2D* PhiScEp_Perp = (TH2D*)f->Get("PhiScEp")->Clone();
//    PhiScEp_Perp->SetName("PhiScEp_Perp");
//    TH2D* PhiScThetan_Perp = (TH2D*)f->Get("PhiScThetan")->Clone();
//    PhiScThetan_Perp->SetName("PhiScThetan_Perp");
//    TH2D* EMWPCnPhiSc_Perp = (TH2D*)f->Get("EMWPCnPhiSc")->Clone();
//    EMWPCnPhiSc_Perp->SetName("EMWPCnPhiSc_Perp");
//
//    TH1D* MMp200300_Perp = (TH1D*)f->Get("MMp200300")->Clone();
//    MMp200300_Perp->SetName("MMp200300_Perp");
//    TH1D* MMp300400_Perp = (TH1D*)f->Get("MMp300400")->Clone();
//    MMp300400_Perp->SetName("MMp300400_Perp");
//    TH1D* MMp400500_Perp = (TH1D*)f->Get("MMp400500")->Clone();
//    MMp400500_Perp->SetName("MMp400500_Perp");
//    TH1D* MMp500600_Perp = (TH1D*)f->Get("MMp500600")->Clone();
//    MMp500600_Perp->SetName("MMp500600_Perp");
//    TH1D* MMp600700_Perp = (TH1D*)f->Get("MMp600700")->Clone();
//    MMp600700_Perp->SetName("MMp600700_Perp");
//    TH1D* MMp700800_Perp = (TH1D*)f->Get("MMp700800")->Clone();
//    MMp700800_Perp->SetName("MMp700800_Perp");
//    TH1D* MMp800900_Perp = (TH1D*)f->Get("MMp800900")->Clone();
//    MMp800900_Perp->SetName("MMp800900_Perp");

    TH1D* PhiSc320_Perp = (TH1D*)f->Get("PhiSc320")->Clone();
    PhiSc320_Perp->SetName("PhiSc320_Perp");
    TH1D* PhiSc360_Perp = (TH1D*)f->Get("PhiSc360")->Clone();
    PhiSc360_Perp->SetName("PhiSc360_Perp");
    TH1D* PhiSc400_Perp = (TH1D*)f->Get("PhiSc400")->Clone();
    PhiSc400_Perp->SetName("PhiSc400_Perp");
    TH1D* PhiSc440_Perp = (TH1D*)f->Get("PhiSc440")->Clone();
    PhiSc440_Perp->SetName("PhiSc440_Perp");
    TH1D* PhiSc480_Perp = (TH1D*)f->Get("PhiSc480")->Clone();
    PhiSc480_Perp->SetName("PhiSc480_Perp");
    TH1D* PhiSc520_Perp = (TH1D*)f->Get("PhiSc520")->Clone();
    PhiSc520_Perp->SetName("PhiSc520_Perp");
    TH1D* PhiSc560_Perp = (TH1D*)f->Get("PhiSc560")->Clone();
    PhiSc560_Perp->SetName("PhiSc560_Perp");
    TH1D* PhiSc600_Perp = (TH1D*)f->Get("PhiSc600")->Clone();
    PhiSc600_Perp->SetName("PhiSc600_Perp");
    TH1D* PhiSc640_Perp = (TH1D*)f->Get("PhiSc640")->Clone();
    PhiSc640_Perp->SetName("PhiSc640_Perp");
    TH1D* PhiSc680_Perp = (TH1D*)f->Get("PhiSc680")->Clone();
    PhiSc680_Perp->SetName("PhiSc680_Perp");

//    Phi_Scattered_265MeV_NegHelCM1_Perp = (TH1D*)f->Get("Phi_Scattered_265MeV_NegHelCM1")->Clone();
//    Phi_Scattered_265MeV_NegHelCM1_Perp->SetName("Phi_Scattered_265MeV_NegHelCM1_Perp");
//    Phi_Scattered_335MeV_NegHelCM1_Perp = (TH1D*)f->Get("Phi_Scattered_335MeV_NegHelCM1")->Clone();
//    Phi_Scattered_335MeV_NegHelCM1_Perp->SetName("Phi_Scattered_335MeV_NegHelCM1_Perp");
//    Phi_Scattered_405MeV_NegHelCM1_Perp = (TH1D*)f->Get("Phi_Scattered_405MeV_NegHelCM1")->Clone();
//    Phi_Scattered_405MeV_NegHelCM1_Perp->SetName("Phi_Scattered_405MeV_NegHelCM1_Perp");
//    Phi_Scattered_475MeV_NegHelCM1_Perp = (TH1D*)f->Get("Phi_Scattered_475MeV_NegHelCM1")->Clone();
//    Phi_Scattered_475MeV_NegHelCM1_Perp->SetName("Phi_Scattered_475MeV_NegHelCM1_Perp");
//    Phi_Scattered_545MeV_NegHelCM1_Perp = (TH1D*)f->Get("Phi_Scattered_545MeV_NegHelCM1")->Clone();
//    Phi_Scattered_545MeV_NegHelCM1_Perp->SetName("Phi_Scattered_545MeV_NegHelCM1_Perp");
//    Phi_Scattered_615MeV_NegHelCM1_Perp = (TH1D*)f->Get("Phi_Scattered_615MeV_NegHelCM1")->Clone();
//    Phi_Scattered_615MeV_NegHelCM1_Perp->SetName("Phi_Scattered_615MeV_NegHelCM1_Perp");
//    Phi_Scattered_685MeV_NegHelCM1_Perp = (TH1D*)f->Get("Phi_Scattered_685MeV_NegHelCM1")->Clone();
//    Phi_Scattered_685MeV_NegHelCM1_Perp->SetName("Phi_Scattered_685MeV_NegHelCM1_Perp");
//
//    Phi_Scattered_265MeV_NegHelCM2_Perp = (TH1D*)f->Get("Phi_Scattered_265MeV_NegHelCM2")->Clone();
//    Phi_Scattered_265MeV_NegHelCM2_Perp->SetName("Phi_Scattered_265MeV_NegHelCM2_Perp");
//    Phi_Scattered_335MeV_NegHelCM2_Perp = (TH1D*)f->Get("Phi_Scattered_335MeV_NegHelCM2")->Clone();
//    Phi_Scattered_335MeV_NegHelCM2_Perp->SetName("Phi_Scattered_335MeV_NegHelCM2_Perp");
//    Phi_Scattered_405MeV_NegHelCM2_Perp = (TH1D*)f->Get("Phi_Scattered_405MeV_NegHelCM2")->Clone();
//    Phi_Scattered_405MeV_NegHelCM2_Perp->SetName("Phi_Scattered_405MeV_NegHelCM2_Perp");
//    Phi_Scattered_475MeV_NegHelCM2_Perp = (TH1D*)f->Get("Phi_Scattered_475MeV_NegHelCM2")->Clone();
//    Phi_Scattered_475MeV_NegHelCM2_Perp->SetName("Phi_Scattered_475MeV_NegHelCM2_Perp");
//    Phi_Scattered_545MeV_NegHelCM2_Perp = (TH1D*)f->Get("Phi_Scattered_545MeV_NegHelCM2")->Clone();
//    Phi_Scattered_545MeV_NegHelCM2_Perp->SetName("Phi_Scattered_545MeV_NegHelCM2_Perp");
//    Phi_Scattered_615MeV_NegHelCM2_Perp = (TH1D*)f->Get("Phi_Scattered_615MeV_NegHelCM2")->Clone();
//    Phi_Scattered_615MeV_NegHelCM2_Perp->SetName("Phi_Scattered_615MeV_NegHelCM2_Perp");
//    Phi_Scattered_685MeV_NegHelCM2_Perp = (TH1D*)f->Get("Phi_Scattered_685MeV_NegHelCM2")->Clone();
//    Phi_Scattered_685MeV_NegHelCM2_Perp->SetName("Phi_Scattered_685MeV_NegHelCM2_Perp");
//
//    Phi_Scattered_265MeV_NegHelCM3_Perp = (TH1D*)f->Get("Phi_Scattered_265MeV_NegHelCM3")->Clone();
//    Phi_Scattered_265MeV_NegHelCM3_Perp->SetName("Phi_Scattered_265MeV_NegHelCM3_Perp");
//    Phi_Scattered_335MeV_NegHelCM3_Perp = (TH1D*)f->Get("Phi_Scattered_335MeV_NegHelCM3")->Clone();
//    Phi_Scattered_335MeV_NegHelCM3_Perp->SetName("Phi_Scattered_335MeV_NegHelCM3_Perp");
//    Phi_Scattered_405MeV_NegHelCM3_Perp = (TH1D*)f->Get("Phi_Scattered_405MeV_NegHelCM3")->Clone();
//    Phi_Scattered_405MeV_NegHelCM3_Perp->SetName("Phi_Scattered_405MeV_NegHelCM3_Perp");
//    Phi_Scattered_475MeV_NegHelCM3_Perp = (TH1D*)f->Get("Phi_Scattered_475MeV_NegHelCM3")->Clone();
//    Phi_Scattered_475MeV_NegHelCM3_Perp->SetName("Phi_Scattered_475MeV_NegHelCM3_Perp");
//    Phi_Scattered_545MeV_NegHelCM3_Perp = (TH1D*)f->Get("Phi_Scattered_545MeV_NegHelCM3")->Clone();
//    Phi_Scattered_545MeV_NegHelCM3_Perp->SetName("Phi_Scattered_545MeV_NegHelCM3_Perp");
//    Phi_Scattered_615MeV_NegHelCM3_Perp = (TH1D*)f->Get("Phi_Scattered_615MeV_NegHelCM3")->Clone();
//    Phi_Scattered_615MeV_NegHelCM3_Perp->SetName("Phi_Scattered_615MeV_NegHelCM3_Perp");
//    Phi_Scattered_685MeV_NegHelCM3_Perp = (TH1D*)f->Get("Phi_Scattered_685MeV_NegHelCM3")->Clone();
//    Phi_Scattered_685MeV_NegHelCM3_Perp->SetName("Phi_Scattered_685MeV_NegHelCM3_Perp");
//
//    Phi_Scattered_265MeV_NegHelCM4_Perp = (TH1D*)f->Get("Phi_Scattered_265MeV_NegHelCM4")->Clone();
//    Phi_Scattered_265MeV_NegHelCM4_Perp->SetName("Phi_Scattered_265MeV_NegHelCM4_Perp");
//    Phi_Scattered_335MeV_NegHelCM4_Perp = (TH1D*)f->Get("Phi_Scattered_335MeV_NegHelCM4")->Clone();
//    Phi_Scattered_335MeV_NegHelCM4_Perp->SetName("Phi_Scattered_335MeV_NegHelCM4_Perp");
//    Phi_Scattered_405MeV_NegHelCM4_Perp = (TH1D*)f->Get("Phi_Scattered_405MeV_NegHelCM4")->Clone();
//    Phi_Scattered_405MeV_NegHelCM4_Perp->SetName("Phi_Scattered_405MeV_NegHelCM4_Perp");
//    Phi_Scattered_475MeV_NegHelCM4_Perp = (TH1D*)f->Get("Phi_Scattered_475MeV_NegHelCM4")->Clone();
//    Phi_Scattered_475MeV_NegHelCM4_Perp->SetName("Phi_Scattered_475MeV_NegHelCM4_Perp");
//    Phi_Scattered_545MeV_NegHelCM4_Perp = (TH1D*)f->Get("Phi_Scattered_545MeV_NegHelCM4")->Clone();
//    Phi_Scattered_545MeV_NegHelCM4_Perp->SetName("Phi_Scattered_545MeV_NegHelCM4_Perp");
//    Phi_Scattered_615MeV_NegHelCM4_Perp = (TH1D*)f->Get("Phi_Scattered_615MeV_NegHelCM4")->Clone();
//    Phi_Scattered_615MeV_NegHelCM4_Perp->SetName("Phi_Scattered_615MeV_NegHelCM4_Perp");
//    Phi_Scattered_685MeV_NegHelCM4_Perp = (TH1D*)f->Get("Phi_Scattered_685MeV_NegHelCM4")->Clone();
//    Phi_Scattered_685MeV_NegHelCM4_Perp->SetName("Phi_Scattered_685MeV_NegHelCM4_Perp");
//
//    Phi_Scattered_265MeV_NegHelCM5_Perp = (TH1D*)f->Get("Phi_Scattered_265MeV_NegHelCM5")->Clone();
//    Phi_Scattered_265MeV_NegHelCM5_Perp->SetName("Phi_Scattered_265MeV_NegHelCM5_Perp");
//    Phi_Scattered_335MeV_NegHelCM5_Perp = (TH1D*)f->Get("Phi_Scattered_335MeV_NegHelCM5")->Clone();
//    Phi_Scattered_335MeV_NegHelCM5_Perp->SetName("Phi_Scattered_335MeV_NegHelCM5_Perp");
//    Phi_Scattered_405MeV_NegHelCM5_Perp = (TH1D*)f->Get("Phi_Scattered_405MeV_NegHelCM5")->Clone();
//    Phi_Scattered_405MeV_NegHelCM5_Perp->SetName("Phi_Scattered_405MeV_NegHelCM5_Perp");
//    Phi_Scattered_475MeV_NegHelCM5_Perp = (TH1D*)f->Get("Phi_Scattered_475MeV_NegHelCM5")->Clone();
//    Phi_Scattered_475MeV_NegHelCM5_Perp->SetName("Phi_Scattered_475MeV_NegHelCM5_Perp");
//    Phi_Scattered_545MeV_NegHelCM5_Perp = (TH1D*)f->Get("Phi_Scattered_545MeV_NegHelCM5")->Clone();
//    Phi_Scattered_545MeV_NegHelCM5_Perp->SetName("Phi_Scattered_545MeV_NegHelCM5_Perp");
//    Phi_Scattered_615MeV_NegHelCM5_Perp = (TH1D*)f->Get("Phi_Scattered_615MeV_NegHelCM5")->Clone();
//    Phi_Scattered_615MeV_NegHelCM5_Perp->SetName("Phi_Scattered_615MeV_NegHelCM5_Perp");
//    Phi_Scattered_685MeV_NegHelCM5_Perp = (TH1D*)f->Get("Phi_Scattered_685MeV_NegHelCM5")->Clone();
//    Phi_Scattered_685MeV_NegHelCM5_Perp->SetName("Phi_Scattered_685MeV_NegHelCM5_Perp");
//
//    Phi_Scattered_265MeV_NegHelCM6_Perp = (TH1D*)f->Get("Phi_Scattered_265MeV_NegHelCM6")->Clone();
//    Phi_Scattered_265MeV_NegHelCM6_Perp->SetName("Phi_Scattered_265MeV_NegHelCM6_Perp");
//    Phi_Scattered_335MeV_NegHelCM6_Perp = (TH1D*)f->Get("Phi_Scattered_335MeV_NegHelCM6")->Clone();
//    Phi_Scattered_335MeV_NegHelCM6_Perp->SetName("Phi_Scattered_335MeV_NegHelCM6_Perp");
//    Phi_Scattered_405MeV_NegHelCM6_Perp = (TH1D*)f->Get("Phi_Scattered_405MeV_NegHelCM6")->Clone();
//    Phi_Scattered_405MeV_NegHelCM6_Perp->SetName("Phi_Scattered_405MeV_NegHelCM6_Perp");
//    Phi_Scattered_475MeV_NegHelCM6_Perp = (TH1D*)f->Get("Phi_Scattered_475MeV_NegHelCM6")->Clone();
//    Phi_Scattered_475MeV_NegHelCM6_Perp->SetName("Phi_Scattered_475MeV_NegHelCM6_Perp");
//    Phi_Scattered_545MeV_NegHelCM6_Perp = (TH1D*)f->Get("Phi_Scattered_545MeV_NegHelCM6")->Clone();
//    Phi_Scattered_545MeV_NegHelCM6_Perp->SetName("Phi_Scattered_545MeV_NegHelCM6_Perp");
//    Phi_Scattered_615MeV_NegHelCM6_Perp = (TH1D*)f->Get("Phi_Scattered_615MeV_NegHelCM6")->Clone();
//    Phi_Scattered_615MeV_NegHelCM6_Perp->SetName("Phi_Scattered_615MeV_NegHelCM6_Perp");
//    Phi_Scattered_685MeV_NegHelCM6_Perp = (TH1D*)f->Get("Phi_Scattered_685MeV_NegHelCM6")->Clone();
//    Phi_Scattered_685MeV_NegHelCM6_Perp->SetName("Phi_Scattered_685MeV_NegHelCM6_Perp");
//
//    Phi_Scattered_265MeV_NegHelCM7_Perp = (TH1D*)f->Get("Phi_Scattered_265MeV_NegHelCM7")->Clone();
//    Phi_Scattered_265MeV_NegHelCM7_Perp->SetName("Phi_Scattered_265MeV_NegHelCM7_Perp");
//    Phi_Scattered_335MeV_NegHelCM7_Perp = (TH1D*)f->Get("Phi_Scattered_335MeV_NegHelCM7")->Clone();
//    Phi_Scattered_335MeV_NegHelCM7_Perp->SetName("Phi_Scattered_335MeV_NegHelCM7_Perp");
//    Phi_Scattered_405MeV_NegHelCM7_Perp = (TH1D*)f->Get("Phi_Scattered_405MeV_NegHelCM7")->Clone();
//    Phi_Scattered_405MeV_NegHelCM7_Perp->SetName("Phi_Scattered_405MeV_NegHelCM7_Perp");
//    Phi_Scattered_475MeV_NegHelCM7_Perp = (TH1D*)f->Get("Phi_Scattered_475MeV_NegHelCM7")->Clone();
//    Phi_Scattered_475MeV_NegHelCM7_Perp->SetName("Phi_Scattered_475MeV_NegHelCM7_Perp");
//    Phi_Scattered_545MeV_NegHelCM7_Perp = (TH1D*)f->Get("Phi_Scattered_545MeV_NegHelCM7")->Clone();
//    Phi_Scattered_545MeV_NegHelCM7_Perp->SetName("Phi_Scattered_545MeV_NegHelCM7_Perp");
//    Phi_Scattered_615MeV_NegHelCM7_Perp = (TH1D*)f->Get("Phi_Scattered_615MeV_NegHelCM7")->Clone();
//    Phi_Scattered_615MeV_NegHelCM7_Perp->SetName("Phi_Scattered_615MeV_NegHelCM7_Perp");
//    Phi_Scattered_685MeV_NegHelCM7_Perp = (TH1D*)f->Get("Phi_Scattered_685MeV_NegHelCM7")->Clone();
//    Phi_Scattered_685MeV_NegHelCM7_Perp->SetName("Phi_Scattered_685MeV_NegHelCM7_Perp");
//
//    Phi_Scattered_265MeV_NegHelCM8_Perp = (TH1D*)f->Get("Phi_Scattered_265MeV_NegHelCM8")->Clone();
//    Phi_Scattered_265MeV_NegHelCM8_Perp->SetName("Phi_Scattered_265MeV_NegHelCM8_Perp");
//    Phi_Scattered_335MeV_NegHelCM8_Perp = (TH1D*)f->Get("Phi_Scattered_335MeV_NegHelCM8")->Clone();
//    Phi_Scattered_335MeV_NegHelCM8_Perp->SetName("Phi_Scattered_335MeV_NegHelCM8_Perp");
//    Phi_Scattered_405MeV_NegHelCM8_Perp = (TH1D*)f->Get("Phi_Scattered_405MeV_NegHelCM8")->Clone();
//    Phi_Scattered_405MeV_NegHelCM8_Perp->SetName("Phi_Scattered_405MeV_NegHelCM8_Perp");
//    Phi_Scattered_475MeV_NegHelCM8_Perp = (TH1D*)f->Get("Phi_Scattered_475MeV_NegHelCM8")->Clone();
//    Phi_Scattered_475MeV_NegHelCM8_Perp->SetName("Phi_Scattered_475MeV_NegHelCM8_Perp");
//    Phi_Scattered_545MeV_NegHelCM8_Perp = (TH1D*)f->Get("Phi_Scattered_545MeV_NegHelCM8")->Clone();
//    Phi_Scattered_545MeV_NegHelCM8_Perp->SetName("Phi_Scattered_545MeV_NegHelCM8_Perp");
//    Phi_Scattered_615MeV_NegHelCM8_Perp = (TH1D*)f->Get("Phi_Scattered_615MeV_NegHelCM8")->Clone();
//    Phi_Scattered_615MeV_NegHelCM8_Perp->SetName("Phi_Scattered_615MeV_NegHelCM8_Perp");
//    Phi_Scattered_685MeV_NegHelCM8_Perp = (TH1D*)f->Get("Phi_Scattered_685MeV_NegHelCM8")->Clone();
//    Phi_Scattered_685MeV_NegHelCM8_Perp->SetName("Phi_Scattered_685MeV_NegHelCM8_Perp");
//
//    Phi_Scattered_265MeV_PosHelCM1_Perp = (TH1D*)f->Get("Phi_Scattered_265MeV_PosHelCM1")->Clone();
//    Phi_Scattered_265MeV_PosHelCM1_Perp->SetName("Phi_Scattered_265MeV_PosHelCM1_Perp");
//    Phi_Scattered_335MeV_PosHelCM1_Perp = (TH1D*)f->Get("Phi_Scattered_335MeV_PosHelCM1")->Clone();
//    Phi_Scattered_335MeV_PosHelCM1_Perp->SetName("Phi_Scattered_335MeV_PosHelCM1_Perp");
//    Phi_Scattered_405MeV_PosHelCM1_Perp = (TH1D*)f->Get("Phi_Scattered_405MeV_PosHelCM1")->Clone();
//    Phi_Scattered_405MeV_PosHelCM1_Perp->SetName("Phi_Scattered_405MeV_PosHelCM1_Perp");
//    Phi_Scattered_475MeV_PosHelCM1_Perp = (TH1D*)f->Get("Phi_Scattered_475MeV_PosHelCM1")->Clone();
//    Phi_Scattered_475MeV_PosHelCM1_Perp->SetName("Phi_Scattered_475MeV_PosHelCM1_Perp");
//    Phi_Scattered_545MeV_PosHelCM1_Perp = (TH1D*)f->Get("Phi_Scattered_545MeV_PosHelCM1")->Clone();
//    Phi_Scattered_545MeV_PosHelCM1_Perp->SetName("Phi_Scattered_545MeV_PosHelCM1_Perp");
//    Phi_Scattered_615MeV_PosHelCM1_Perp = (TH1D*)f->Get("Phi_Scattered_615MeV_PosHelCM1")->Clone();
//    Phi_Scattered_615MeV_PosHelCM1_Perp->SetName("Phi_Scattered_615MeV_PosHelCM1_Perp");
//    Phi_Scattered_685MeV_PosHelCM1_Perp = (TH1D*)f->Get("Phi_Scattered_685MeV_PosHelCM1")->Clone();
//    Phi_Scattered_685MeV_PosHelCM1_Perp->SetName("Phi_Scattered_685MeV_PosHelCM1_Perp");
//
//    Phi_Scattered_265MeV_PosHelCM2_Perp = (TH1D*)f->Get("Phi_Scattered_265MeV_PosHelCM2")->Clone();
//    Phi_Scattered_265MeV_PosHelCM2_Perp->SetName("Phi_Scattered_265MeV_PosHelCM2_Perp");
//    Phi_Scattered_335MeV_PosHelCM2_Perp = (TH1D*)f->Get("Phi_Scattered_335MeV_PosHelCM2")->Clone();
//    Phi_Scattered_335MeV_PosHelCM2_Perp->SetName("Phi_Scattered_335MeV_PosHelCM2_Perp");
//    Phi_Scattered_405MeV_PosHelCM2_Perp = (TH1D*)f->Get("Phi_Scattered_405MeV_PosHelCM2")->Clone();
//    Phi_Scattered_405MeV_PosHelCM2_Perp->SetName("Phi_Scattered_405MeV_PosHelCM2_Perp");
//    Phi_Scattered_475MeV_PosHelCM2_Perp = (TH1D*)f->Get("Phi_Scattered_475MeV_PosHelCM2")->Clone();
//    Phi_Scattered_475MeV_PosHelCM2_Perp->SetName("Phi_Scattered_475MeV_PosHelCM2_Perp");
//    Phi_Scattered_545MeV_PosHelCM2_Perp = (TH1D*)f->Get("Phi_Scattered_545MeV_PosHelCM2")->Clone();
//    Phi_Scattered_545MeV_PosHelCM2_Perp->SetName("Phi_Scattered_545MeV_PosHelCM2_Perp");
//    Phi_Scattered_615MeV_PosHelCM2_Perp = (TH1D*)f->Get("Phi_Scattered_615MeV_PosHelCM2")->Clone();
//    Phi_Scattered_615MeV_PosHelCM2_Perp->SetName("Phi_Scattered_615MeV_PosHelCM2_Perp");
//    Phi_Scattered_685MeV_PosHelCM2_Perp = (TH1D*)f->Get("Phi_Scattered_685MeV_PosHelCM2")->Clone();
//    Phi_Scattered_685MeV_PosHelCM2_Perp->SetName("Phi_Scattered_685MeV_PosHelCM2_Perp");
//
//    Phi_Scattered_265MeV_PosHelCM3_Perp = (TH1D*)f->Get("Phi_Scattered_265MeV_PosHelCM3")->Clone();
//    Phi_Scattered_265MeV_PosHelCM3_Perp->SetName("Phi_Scattered_265MeV_PosHelCM3_Perp");
//    Phi_Scattered_335MeV_PosHelCM3_Perp = (TH1D*)f->Get("Phi_Scattered_335MeV_PosHelCM3")->Clone();
//    Phi_Scattered_335MeV_PosHelCM3_Perp->SetName("Phi_Scattered_335MeV_PosHelCM3_Perp");
//    Phi_Scattered_405MeV_PosHelCM3_Perp = (TH1D*)f->Get("Phi_Scattered_405MeV_PosHelCM3")->Clone();
//    Phi_Scattered_405MeV_PosHelCM3_Perp->SetName("Phi_Scattered_405MeV_PosHelCM3_Perp");
//    Phi_Scattered_475MeV_PosHelCM3_Perp = (TH1D*)f->Get("Phi_Scattered_475MeV_PosHelCM3")->Clone();
//    Phi_Scattered_475MeV_PosHelCM3_Perp->SetName("Phi_Scattered_475MeV_PosHelCM3_Perp");
//    Phi_Scattered_545MeV_PosHelCM3_Perp = (TH1D*)f->Get("Phi_Scattered_545MeV_PosHelCM3")->Clone();
//    Phi_Scattered_545MeV_PosHelCM3_Perp->SetName("Phi_Scattered_545MeV_PosHelCM3_Perp");
//    Phi_Scattered_615MeV_PosHelCM3_Perp = (TH1D*)f->Get("Phi_Scattered_615MeV_PosHelCM3")->Clone();
//    Phi_Scattered_615MeV_PosHelCM3_Perp->SetName("Phi_Scattered_615MeV_PosHelCM3_Perp");
//    Phi_Scattered_685MeV_PosHelCM3_Perp = (TH1D*)f->Get("Phi_Scattered_685MeV_PosHelCM3")->Clone();
//    Phi_Scattered_685MeV_PosHelCM3_Perp->SetName("Phi_Scattered_685MeV_PosHelCM3_Perp");
//
//    Phi_Scattered_265MeV_PosHelCM4_Perp = (TH1D*)f->Get("Phi_Scattered_265MeV_PosHelCM4")->Clone();
//    Phi_Scattered_265MeV_PosHelCM4_Perp->SetName("Phi_Scattered_265MeV_PosHelCM4_Perp");
//    Phi_Scattered_335MeV_PosHelCM4_Perp = (TH1D*)f->Get("Phi_Scattered_335MeV_PosHelCM4")->Clone();
//    Phi_Scattered_335MeV_PosHelCM4_Perp->SetName("Phi_Scattered_335MeV_PosHelCM4_Perp");
//    Phi_Scattered_405MeV_PosHelCM4_Perp = (TH1D*)f->Get("Phi_Scattered_405MeV_PosHelCM4")->Clone();
//    Phi_Scattered_405MeV_PosHelCM4_Perp->SetName("Phi_Scattered_405MeV_PosHelCM4_Perp");
//    Phi_Scattered_475MeV_PosHelCM4_Perp = (TH1D*)f->Get("Phi_Scattered_475MeV_PosHelCM4")->Clone();
//    Phi_Scattered_475MeV_PosHelCM4_Perp->SetName("Phi_Scattered_475MeV_PosHelCM4_Perp");
//    Phi_Scattered_545MeV_PosHelCM4_Perp = (TH1D*)f->Get("Phi_Scattered_545MeV_PosHelCM4")->Clone();
//    Phi_Scattered_545MeV_PosHelCM4_Perp->SetName("Phi_Scattered_545MeV_PosHelCM4_Perp");
//    Phi_Scattered_615MeV_PosHelCM4_Perp = (TH1D*)f->Get("Phi_Scattered_615MeV_PosHelCM4")->Clone();
//    Phi_Scattered_615MeV_PosHelCM4_Perp->SetName("Phi_Scattered_615MeV_PosHelCM4_Perp");
//    Phi_Scattered_685MeV_PosHelCM4_Perp = (TH1D*)f->Get("Phi_Scattered_685MeV_PosHelCM4")->Clone();
//    Phi_Scattered_685MeV_PosHelCM4_Perp->SetName("Phi_Scattered_685MeV_PosHelCM4_Perp");
//
//    Phi_Scattered_265MeV_PosHelCM5_Perp = (TH1D*)f->Get("Phi_Scattered_265MeV_PosHelCM5")->Clone();
//    Phi_Scattered_265MeV_PosHelCM5_Perp->SetName("Phi_Scattered_265MeV_PosHelCM5_Perp");
//    Phi_Scattered_335MeV_PosHelCM5_Perp = (TH1D*)f->Get("Phi_Scattered_335MeV_PosHelCM5")->Clone();
//    Phi_Scattered_335MeV_PosHelCM5_Perp->SetName("Phi_Scattered_335MeV_PosHelCM5_Perp");
//    Phi_Scattered_405MeV_PosHelCM5_Perp = (TH1D*)f->Get("Phi_Scattered_405MeV_PosHelCM5")->Clone();
//    Phi_Scattered_405MeV_PosHelCM5_Perp->SetName("Phi_Scattered_405MeV_PosHelCM5_Perp");
//    Phi_Scattered_475MeV_PosHelCM5_Perp = (TH1D*)f->Get("Phi_Scattered_475MeV_PosHelCM5")->Clone();
//    Phi_Scattered_475MeV_PosHelCM5_Perp->SetName("Phi_Scattered_475MeV_PosHelCM5_Perp");
//    Phi_Scattered_545MeV_PosHelCM5_Perp = (TH1D*)f->Get("Phi_Scattered_545MeV_PosHelCM5")->Clone();
//    Phi_Scattered_545MeV_PosHelCM5_Perp->SetName("Phi_Scattered_545MeV_PosHelCM5_Perp");
//    Phi_Scattered_615MeV_PosHelCM5_Perp = (TH1D*)f->Get("Phi_Scattered_615MeV_PosHelCM5")->Clone();
//    Phi_Scattered_615MeV_PosHelCM5_Perp->SetName("Phi_Scattered_615MeV_PosHelCM5_Perp");
//    Phi_Scattered_685MeV_PosHelCM5_Perp = (TH1D*)f->Get("Phi_Scattered_685MeV_PosHelCM5")->Clone();
//    Phi_Scattered_685MeV_PosHelCM5_Perp->SetName("Phi_Scattered_685MeV_PosHelCM5_Perp");
//
//    Phi_Scattered_265MeV_PosHelCM6_Perp = (TH1D*)f->Get("Phi_Scattered_265MeV_PosHelCM6")->Clone();
//    Phi_Scattered_265MeV_PosHelCM6_Perp->SetName("Phi_Scattered_265MeV_PosHelCM6_Perp");
//    Phi_Scattered_335MeV_PosHelCM6_Perp = (TH1D*)f->Get("Phi_Scattered_335MeV_PosHelCM6")->Clone();
//    Phi_Scattered_335MeV_PosHelCM6_Perp->SetName("Phi_Scattered_335MeV_PosHelCM6_Perp");
//    Phi_Scattered_405MeV_PosHelCM6_Perp = (TH1D*)f->Get("Phi_Scattered_405MeV_PosHelCM6")->Clone();
//    Phi_Scattered_405MeV_PosHelCM6_Perp->SetName("Phi_Scattered_405MeV_PosHelCM6_Perp");
//    Phi_Scattered_475MeV_PosHelCM6_Perp = (TH1D*)f->Get("Phi_Scattered_475MeV_PosHelCM6")->Clone();
//    Phi_Scattered_475MeV_PosHelCM6_Perp->SetName("Phi_Scattered_475MeV_PosHelCM6_Perp");
//    Phi_Scattered_545MeV_PosHelCM6_Perp = (TH1D*)f->Get("Phi_Scattered_545MeV_PosHelCM6")->Clone();
//    Phi_Scattered_545MeV_PosHelCM6_Perp->SetName("Phi_Scattered_545MeV_PosHelCM6_Perp");
//    Phi_Scattered_615MeV_PosHelCM6_Perp = (TH1D*)f->Get("Phi_Scattered_615MeV_PosHelCM6")->Clone();
//    Phi_Scattered_615MeV_PosHelCM6_Perp->SetName("Phi_Scattered_615MeV_PosHelCM6_Perp");
//    Phi_Scattered_685MeV_PosHelCM6_Perp = (TH1D*)f->Get("Phi_Scattered_685MeV_PosHelCM6")->Clone();
//    Phi_Scattered_685MeV_PosHelCM6_Perp->SetName("Phi_Scattered_685MeV_PosHelCM6_Perp");
//
//    Phi_Scattered_265MeV_PosHelCM7_Perp = (TH1D*)f->Get("Phi_Scattered_265MeV_PosHelCM7")->Clone();
//    Phi_Scattered_265MeV_PosHelCM7_Perp->SetName("Phi_Scattered_265MeV_PosHelCM7_Perp");
//    Phi_Scattered_335MeV_PosHelCM7_Perp = (TH1D*)f->Get("Phi_Scattered_335MeV_PosHelCM7")->Clone();
//    Phi_Scattered_335MeV_PosHelCM7_Perp->SetName("Phi_Scattered_335MeV_PosHelCM7_Perp");
//    Phi_Scattered_405MeV_PosHelCM7_Perp = (TH1D*)f->Get("Phi_Scattered_405MeV_PosHelCM7")->Clone();
//    Phi_Scattered_405MeV_PosHelCM7_Perp->SetName("Phi_Scattered_405MeV_PosHelCM7_Perp");
//    Phi_Scattered_475MeV_PosHelCM7_Perp = (TH1D*)f->Get("Phi_Scattered_475MeV_PosHelCM7")->Clone();
//    Phi_Scattered_475MeV_PosHelCM7_Perp->SetName("Phi_Scattered_475MeV_PosHelCM7_Perp");
//    Phi_Scattered_545MeV_PosHelCM7_Perp = (TH1D*)f->Get("Phi_Scattered_545MeV_PosHelCM7")->Clone();
//    Phi_Scattered_545MeV_PosHelCM7_Perp->SetName("Phi_Scattered_545MeV_PosHelCM7_Perp");
//    Phi_Scattered_615MeV_PosHelCM7_Perp = (TH1D*)f->Get("Phi_Scattered_615MeV_PosHelCM7")->Clone();
//    Phi_Scattered_615MeV_PosHelCM7_Perp->SetName("Phi_Scattered_615MeV_PosHelCM7_Perp");
//    Phi_Scattered_685MeV_PosHelCM7_Perp = (TH1D*)f->Get("Phi_Scattered_685MeV_PosHelCM7")->Clone();
//    Phi_Scattered_685MeV_PosHelCM7_Perp->SetName("Phi_Scattered_685MeV_PosHelCM7_Perp");
//
//    Phi_Scattered_265MeV_PosHelCM8_Perp = (TH1D*)f->Get("Phi_Scattered_265MeV_PosHelCM8")->Clone();
//    Phi_Scattered_265MeV_PosHelCM8_Perp->SetName("Phi_Scattered_265MeV_PosHelCM8_Perp");
//    Phi_Scattered_335MeV_PosHelCM8_Perp = (TH1D*)f->Get("Phi_Scattered_335MeV_PosHelCM8")->Clone();
//    Phi_Scattered_335MeV_PosHelCM8_Perp->SetName("Phi_Scattered_335MeV_PosHelCM8_Perp");
//    Phi_Scattered_405MeV_PosHelCM8_Perp = (TH1D*)f->Get("Phi_Scattered_405MeV_PosHelCM8")->Clone();
//    Phi_Scattered_405MeV_PosHelCM8_Perp->SetName("Phi_Scattered_405MeV_PosHelCM8_Perp");
//    Phi_Scattered_475MeV_PosHelCM8_Perp = (TH1D*)f->Get("Phi_Scattered_475MeV_PosHelCM8")->Clone();
//    Phi_Scattered_475MeV_PosHelCM8_Perp->SetName("Phi_Scattered_475MeV_PosHelCM8_Perp");
//    Phi_Scattered_545MeV_PosHelCM8_Perp = (TH1D*)f->Get("Phi_Scattered_545MeV_PosHelCM8")->Clone();
//    Phi_Scattered_545MeV_PosHelCM8_Perp->SetName("Phi_Scattered_545MeV_PosHelCM8_Perp");
//    Phi_Scattered_615MeV_PosHelCM8_Perp = (TH1D*)f->Get("Phi_Scattered_615MeV_PosHelCM8")->Clone();
//    Phi_Scattered_615MeV_PosHelCM8_Perp->SetName("Phi_Scattered_615MeV_PosHelCM8_Perp");
//    Phi_Scattered_685MeV_PosHelCM8_Perp = (TH1D*)f->Get("Phi_Scattered_685MeV_PosHelCM8")->Clone();
//    Phi_Scattered_685MeV_PosHelCM8_Perp->SetName("Phi_Scattered_685MeV_PosHelCM8_Perp");

    ////////////////////////////////////////////////////
    ////////////////////////////////////////////////////
    //// ALL HISTS CLONED - NEED TO SCALE AND MERGE ////
    ////////////////////////////////////////////////////
    ////////////////////////////////////////////////////

    NPara = Eg_Para->GetEntries();
    NPerp = Eg_Perp->GetEntries();
    ScaleFactor = NPara/NPerp;
    ScaleFactorErr = sqrt( (NPara/((TMath::Power(NPerp,2)))) + (((TMath::Power(NPara,2)))/(TMath::Power(NPerp,3))) ); // Error Propagation of above formula, see notebook 5

    // Now scale all Perp histograms and merge to form new ones, suffix Unpol

//    TList *list1 = new TList;
//    list1->Add(time_Para);
//    time_Perp->Scale(ScaleFactor);
//    list1->Add(time_Perp);
//    time = new TH1D("time", "time", 1400, -700, 700);
//    time->Merge(list1);
//
//    TList *list2 = new TList;
//    list2->Add(time_cut_Para);
//    time_cut_Perp->Scale(ScaleFactor);
//    list2->Add(time_cut_Perp);
//    time_cut = new TH1D("time_cut", "time_cut", 1400, -700, 700);
//    time_cut->Merge(list2);
//
//    TList *list3 = new TList;
//    list3->Add(Eg_Para);
//    Eg_Perp->Scale(ScaleFactor);
//    list3->Add(Eg_Perp);
//    Eg = new TH1D("Eg", "E_{#gamma} Distribution", 200, 100, 1600);
//    Eg->Merge(list3);
//
//    TList *list4 = new TList;
//    list4->Add(PhiDet_Para);
//    PhiDet_Perp->Scale(ScaleFactor);
//    list4->Add(PhiDet_Perp);
//    PhiDet = new TH1D("PhiDet", "#phi_{Det} distribution for Neutrons", 360, -180, 180);
//    PhiDet->Merge(list4);
//
//    TList *list5 = new TList;
//    list5->Add(PhiRec_Para);
//    PhiRec_Perp->Scale(ScaleFactor);
//    list5->Add(PhiRec_Perp);
//    PhiRec = new TH1D("PhiRec", "#phi_{Rec} distribution for Neutrons", 360, -180, 180);
//    PhiRec->Merge(list5);
//
//    TList *list6 = new TList;
//    list6->Add(Theta_Scattered_Para);
//    Theta_Scattered_Perp->Scale(ScaleFactor);
//    list6->Add(Theta_Scattered_Perp);
//    Theta_Scattered = new TH1D("Theta_Scattered", "#theta_{sc} Proton Distribution", 200, 0, 4);
//    Theta_Scattered->Merge(list6);
//
//    TList *list7 = new TList;
//    list7->Add(Phi_Scattered_Para);
//    Phi_Scattered_Perp->Scale(ScaleFactor);
//    list7->Add(Phi_Scattered_Perp);
//    Phi_Scattered = new TH1D("Phi_Scattered", "#phi_{sc} Proton Distribution", 200, -4, 4);
//    Phi_Scattered->Merge(list7);
//
//    TList *list8 = new TList;
//    list8->Add(MMpEpCorrected_Para);
//    MMpEpCorrected_Perp->Scale(ScaleFactor);
//    list8->Add(MMpEpCorrected_Perp);
//    MMpEpCorrected = new TH1D("MMpEpCorrected", "Missing mass seen by Proton (E Loss Corrected)", 400, 0, 2000);
//    MMpEpCorrected->Merge(list8);
//
//    TList *list9 = new TList;
//    list9->Add(ZpDist_Para);
//    ZpDist_Perp->Scale(ScaleFactor);
//    list9->Add(ZpDist_Perp);
//    ZpDist = new TH1D("ZpDist", "Proton Pseudo Z Vertex Distribution", 200, -400, 400);
//    ZpDist->Merge(list9);
//
//    TList *list10 = new TList;
//    list10->Add(ThetanDist_Para);
//    ThetanDist_Perp->Scale(ScaleFactor);
//    list10->Add(ThetanDist_Perp);
//    ThetanDist = new TH1D("ThetanDist", "#theta_{n} Distribution", 200, 0, 180);
//    ThetanDist->Merge(list10);
//
//    TList *list11 = new TList;
//    list11->Add(ThetanCMDist_Para);
//    ThetanCMDist_Perp->Scale(ScaleFactor);
//    list11->Add(ThetanCMDist_Perp);
//    ThetanCMDist = new TH1D("ThetanCMDist", "#theta_{nCM} Distribution", 200, 0, 180);
//    ThetanCMDist->Merge(list11);
//
//    TList *list13 = new TList;
//    list13->Add(E_dE_Para);
//    E_dE_Perp->Scale(ScaleFactor);
//    list13->Add(E_dE_Perp);
//    E_dE = new TH2D("E_dE", "EdE Plot With E Loss Adjustment", 100, 0, 500, 100, 0, 5);
//    E_dE->Merge(list13);
//
//    TList *list14 = new TList;
//    list14->Add(ThetaScPhiSc_Para);
//    ThetaScPhiSc_Perp->Scale(ScaleFactor);
//    list14->Add(ThetaScPhiSc_Perp);
//    ThetaScPhiSc = new TH2D("ThetaScPhiSc", "#Phi_{Sc} as a function of #theta_{Sc}", 100, 0, 1, 100, -4, 4);
//    ThetaScPhiSc->Merge(list14);
//
//    TList *list15 = new TList;
//    list15->Add(EdEMWPCp_Para);
//    EdEMWPCp_Perp->Scale(ScaleFactor);
//    list15->Add(EdEMWPCp_Perp);
//    EdEMWPCp = new TH2D("EdEMWPCp", "EdEMWPC0 Plot for Proton Track", 200, 0, 500, 200, 0, 400);
//    EdEMWPCp->Merge(list15);
//
//    TList *list16 = new TList;
//    list16->Add(EdEMWPCn_Para);
//    EdEMWPCn_Perp->Scale(ScaleFactor);
//    list16->Add(EdEMWPCn_Perp);
//    EdEMWPCn = new TH2D("EdEMWPCn", "EdEMWPC0 Plot for Neutron Track", 200, 0, 500, 200, 0, 400);
//    EdEMWPCn->Merge(list16);
//
//    TList *list17 = new TList;
//    list17->Add(ClosestApproach_Para);
//    ClosestApproach_Perp->Scale(ScaleFactor);
//    list17->Add(ClosestApproach_Perp);
//    ClosestApproach = new TH1D("ClosestApproach", "DOCA of n and p' vectors", 200, -200, 200);
//    ClosestApproach->Merge(list17);
//
//    TList *list18 = new TList;
//    list18->Add(POCAr_Para);
//    POCAr_Perp->Scale(ScaleFactor);
//    list18->Add(POCAr_Perp);
//    POCAr = new TH1D("POCAr", "r_{POCA}", 200, 0, 150);
//    POCAr->Merge(list18);
//
//    TList *list19 = new TList;
//    list19->Add(ScatterVertexZ_Para);
//    ScatterVertexZ_Perp->Scale(ScaleFactor);
//    list19->Add(ScatterVertexZ_Perp);
//    ScatterVertexZ = new TH1D("ScatterVertexZ", "Z_{POCA}", 200, -200, 200);
//    ScatterVertexZ->Merge(list19);
//
//    TList *list20 = new TList;
//    list20->Add(ScatterVertexZr_Para);
//    ScatterVertexZr_Perp->Scale(ScaleFactor);
//    list20->Add(ScatterVertexZr_Perp);
//    ScatterVertexZr = new TH2D("ScatterVertexZr", "Z_{POCA} vs r_{POCA}", 200, -150, 150, 200, 0, 150);
//    ScatterVertexZr->Merge(list20);
//
//    TList *list21 = new TList;
//    list21->Add(ScatterVertexXY_Para);
//    ScatterVertexXY_Perp->Scale(ScaleFactor);
//    list21->Add(ScatterVertexXY_Perp);
//    ScatterVertexXY = new TH2D("ScatterVertexXY", "XY Vertex Point of Scatter from DOCA Method", 100, -100, 100, 100, -100, 100);
//    ScatterVertexXY->Merge(list21);
//
//    TList *list22 = new TList;
//    list22->Add(ScatterVertex_Para);
//    ScatterVertex_Perp->Scale(ScaleFactor);
//    list22->Add(ScatterVertex_Perp);
//    ScatterVertex = new TH3D("ScatterVertex", "Vertex Point of Scatter from DOCA Method", 100, -80, 80, 100, -80, 80, 100, -200, 200);
//    ScatterVertex->Merge(list22);
//
//    TList *list23 = new TList;
//    list23->Add(POCArPhiSc_Para);
//    POCArPhiSc_Perp->Scale(ScaleFactor);
//    list23->Add(POCArPhiSc_Perp);
//    POCArPhiSc = new TH2D("POCArPhiSc", "#phi_{Sc} as a Function of r_{POCA}", 100, 0, 200, 100, -4, 4);
//    POCArPhiSc->Merge(list23);
//
//    TList *list24 = new TList;
//    list24->Add(ThetapCorrDiff_Para);
//    ThetapCorrDiff_Perp->Scale(ScaleFactor);
//    list24->Add(ThetapCorrDiff_Perp);
//    ThetapCorrDiff = new TH1D("ThetapCorrDiff", "Difference Between #theta_{p} and #theta_{pCorr}", 200, -180, 180);
//    ThetapCorrDiff->Merge(list24);
//
//    TList *list25 = new TList;
//    list25->Add(PhipCorrDiff_Para);
//    PhipCorrDiff_Perp->Scale(ScaleFactor);
//    list25->Add(PhipCorrDiff_Perp);
//    PhipCorrDiff = new TH1D("PhipCorrDiff", "Difference Between #phi_{p} and #phi_{pCorr}", 200, -180, 180);
//    PhipCorrDiff->Merge(list25);
//
//    TList *list26 = new TList;
//    list26->Add(ThetaDiff_Para);
//    ThetaDiff_Perp->Scale(ScaleFactor);
//    list26->Add(ThetaDiff_Perp);
//    ThetaDiff = new TH1D("ThetaDiff", "Difference Between #theta_{Det} and #theta_{Rec}", 200, -90, 90);
//    ThetaDiff->Merge(list26);
//
//    TList *list27 = new TList;
//    list27->Add(PhiDiff_Para);
//    PhiDiff_Perp->Scale(ScaleFactor);
//    list27->Add(PhiDiff_Perp);
//    PhiDiff = new TH1D("PhiDiff", "Difference Between #phi_{Det} and #phi_{Rec}", 360, -180, 180);
//    PhiDiff->Merge(list27);
//
//    TList *list28 = new TList;
//    list28->Add(PhiDiffThetaDiff_Para);
//    PhiDiffThetaDiff_Perp->Scale(ScaleFactor);
//    list28->Add(PhiDiffThetaDiff_Perp);
//    PhiDiffThetaDiff = new TH2D("PhiDiffThetaDiff", "#phi_{Diff} as a Function of #theta_{Diff}", 100, -100, 100, 100, -100, 100);
//    PhiDiffThetaDiff->Merge(list28);
//
//    TList *list29 = new TList;
//    list29->Add(PhiScEg_Para);
//    PhiScEg_Perp->Scale(ScaleFactor);
//    list29->Add(PhiScEg_Perp);
//    PhiScEg = new TH2D("PhiScEg", "#phi_{Sc} as a Function of E_{#gamma}", 150, 200, 800, 150, -4, 4);
//    PhiScEg->Merge(list29);
//
//    TList *list30 = new TList;
//    list30->Add(PhiScEp_Para);
//    PhiScEp_Perp->Scale(ScaleFactor);
//    list30->Add(PhiScEp_Perp);
//    PhiScEp = new TH2D("PhiScEp", "#phi_{Sc} as a Function of E_{p}", 150, 100, 500, 150, -4, 4);
//    PhiScEp->Merge(list30);
//
//    TList *list31 = new TList;
//    list31->Add(PhiScThetan_Para);
//    PhiScThetan_Perp->Scale(ScaleFactor);
//    list31->Add(PhiScThetan_Perp);
//    PhiScThetan = new TH2D("PhiScThetan", "#phi_{Sc} as a Function of #theta_{n}", 150, 0, 4, 150, -4, 4);
//    PhiScThetan->Merge(list31);
//
//    TList *list32 = new TList;
//    list32->Add(EMWPCnPhiSc_Para);
//    EMWPCnPhiSc_Perp->Scale(ScaleFactor);
//    list32->Add(EMWPCnPhiSc_Perp);
//    EMWPCnPhiSc = new TH2D("EMWPCnPhiSc", "#phi_{Sc} as a Function of MWPC E_{Sum}", 200, 0, 750, 200, 0, 4);
//    EMWPCnPhiSc->Merge(list32);

    TList *list33 = new TList;
    list33->Add(PhiSc320_Para);
    PhiSc320_Perp->Scale(ScaleFactor);
    list33->Add(PhiSc320_Perp);
    PhiSc320 = new TH1D("PhiSc320", "#phi_{Sc} (300-340MeV)", 24, -4, 4);
    PhiSc320->Merge(list33);

    TList *list34 = new TList;
    list34->Add(PhiSc360_Para);
    PhiSc360_Perp->Scale(ScaleFactor);
    list34->Add(PhiSc360_Perp);
    PhiSc360 = new TH1D("PhiSc360", "#phi_{Sc} (340-380MeV)", 24, -4, 4);
    PhiSc360->Merge(list34);

    TList *list35 = new TList;
    list35->Add(PhiSc400_Para);
    PhiSc400_Perp->Scale(ScaleFactor);
    list35->Add(PhiSc400_Perp);
    PhiSc400 = new TH1D("PhiSc400", "#phi_{Sc} (380-420MeV)", 24, -4, 4);
    PhiSc400->Merge(list35);

    TList *list36 = new TList;
    list36->Add(PhiSc440_Para);
    PhiSc440_Perp->Scale(ScaleFactor);
    list36->Add(PhiSc440_Perp);
    PhiSc440 = new TH1D("PhiSc440", "#phi_{Sc} (420-460MeV)", 24, -4, 4);
    PhiSc440->Merge(list36);

    TList *list37 = new TList;
    list37->Add(PhiSc480_Para);
    PhiSc480_Perp->Scale(ScaleFactor);
    list37->Add(PhiSc480_Perp);
    PhiSc480 = new TH1D("PhiSc480", "#phi_{Sc} (460-500MeV)", 24, -4, 4);
    PhiSc480->Merge(list37);

    TList *list38 = new TList;
    list38->Add(PhiSc520_Para);
    PhiSc520_Perp->Scale(ScaleFactor);
    list38->Add(PhiSc520_Perp);
    PhiSc520 = new TH1D("PhiSc520", "#phi_{Sc} (500-540MeV)", 24, -4, 4);
    PhiSc520->Merge(list38);

    TList *list39 = new TList;
    list39->Add(PhiSc560_Para);
    PhiSc560_Perp->Scale(ScaleFactor);
    list39->Add(PhiSc560_Perp);
    PhiSc560 = new TH1D("PhiSc560", "#phi_{Sc} (540-580MeV)", 24, -4, 4);
    PhiSc560->Merge(list39);

    TList *list40 = new TList;
    list40->Add(PhiSc600_Para);
    PhiSc600_Perp->Scale(ScaleFactor);
    list40->Add(PhiSc600_Perp);
    PhiSc600 = new TH1D("PhiSc600", "#phi_{Sc} (580-620MeV)", 24, -4, 4);
    PhiSc600->Merge(list40);

    TList *list41 = new TList;
    list41->Add(PhiSc640_Para);
    PhiSc640_Perp->Scale(ScaleFactor);
    list41->Add(PhiSc640_Perp);
    PhiSc640 = new TH1D("PhiSc640", "#phi_{Sc} (620-660MeV)", 24, -4, 4);
    PhiSc640->Merge(list41);

    TList *list42 = new TList;
    list42->Add(PhiSc680_Para);
    PhiSc680_Perp->Scale(ScaleFactor);
    list42->Add(PhiSc680_Perp);
    PhiSc680 = new TH1D("PhiSc680", "#phi_{Sc} (660-700MeV)", 24, -4, 4);
    PhiSc680->Merge(list42);

//    TList *list43 = new TList;
//    list43->Add(MMp200300_Para);
//    MMp200300_Perp->Scale(ScaleFactor);
//    list43->Add(MMp200300_Perp);
//    MMp200300 = new TH1D("MMp200300", "Missing mass as seen by Proton (200-300MeV E_{#gamma})", 400, 0, 2000);
//    MMp200300->Merge(list43);
//
//    TList *list44 = new TList;
//    list44->Add(MMp300400_Para);
//    MMp300400_Perp->Scale(ScaleFactor);
//    list44->Add(MMp300400_Perp);
//    MMp300400 = new TH1D("MMp300400", "Missing mass as seen by Proton (300-400MeV E_{#gamma})", 400, 0, 2000);
//    MMp300400->Merge(list44);
//
//    TList *list45 = new TList;
//    list45->Add(MMp400500_Para);
//    MMp400500_Perp->Scale(ScaleFactor);
//    list45->Add(MMp400500_Perp);
//    MMp400500 = new TH1D("MMp400500", "Missing mass as seen by Proton (400-500MeV E_{#gamma})", 400, 0, 2000);
//    MMp400500->Merge(list45);
//
//    TList *list46 = new TList;
//    list46->Add(MMp500600_Para);
//    MMp500600_Perp->Scale(ScaleFactor);
//    list46->Add(MMp500600_Perp);
//    MMp500600 = new TH1D("MMp500600", "Missing mass as seen by Proton (500-600MeV E_{#gamma})", 400, 0, 2000);
//    MMp500600->Merge(list46);
//
//    TList *list47 = new TList;
//    list47->Add(MMp600700_Para);
//    MMp600700_Perp->Scale(ScaleFactor);
//    list47->Add(MMp600700_Perp);
//    MMp600700 = new TH1D("MMp600700", "Missing mass as seen by Proton (600-700MeV E_{#gamma})", 400, 0, 2000);
//    MMp600700->Merge(list47);
//
//    TList *list48 = new TList;
//    list48->Add(MMp700800_Para);
//    MMp700800_Perp->Scale(ScaleFactor);
//    list48->Add(MMp700800_Perp);
//    MMp700800 = new TH1D("MMp700800", "Missing mass as seen by Proton (700-800MeV E_{#gamma})", 400, 0, 2000);
//    MMp700800->Merge(list48);
//
//    TList *list49 = new TList;
//    list49->Add(MMp800900_Para);
//    MMp800900_Perp->Scale(ScaleFactor);
//    list49->Add(MMp800900_Perp);
//    MMp800900 = new TH1D("MMp800900", "Missing mass as seen by Proton (800-900MeV E_{#gamma})", 400, 0, 2000);
//    MMp800900->Merge(list49);
//
//    /////////////////////////////////////
    ////// All Neg Hel PhiSc Dists //////
    /////////////////////////////////////

//    TList *list50 = new TList;
//    list50->Add(Phi_Scattered_265MeV_NegHelCM1_Para);
//    Phi_Scattered_265MeV_NegHelCM1_Perp->Scale(ScaleFactor);
//    list50->Add(Phi_Scattered_265MeV_NegHelCM1_Perp);
//    Phi_Scattered_265MeV_NegHelCM1 = new TH1D( "Phi_Scattered_265MeV_NegHelCM1", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_265MeV_NegHelCM1->Merge(list50);
//
//    TList *list51 = new TList;
//    list51->Add(Phi_Scattered_335MeV_NegHelCM1_Para);
//    Phi_Scattered_335MeV_NegHelCM1_Perp->Scale(ScaleFactor);
//    list51->Add(Phi_Scattered_335MeV_NegHelCM1_Perp);
//    Phi_Scattered_335MeV_NegHelCM1 = new TH1D( "Phi_Scattered_335MeV_NegHelCM1", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_335MeV_NegHelCM1->Merge(list51);
//
//    TList *list52 = new TList;
//    list52->Add(Phi_Scattered_405MeV_NegHelCM1_Para);
//    Phi_Scattered_405MeV_NegHelCM1_Perp->Scale(ScaleFactor);
//    list52->Add(Phi_Scattered_405MeV_NegHelCM1_Perp);
//    Phi_Scattered_405MeV_NegHelCM1 = new TH1D( "Phi_Scattered_405MeV_NegHelCM1", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_405MeV_NegHelCM1->Merge(list52);
//
//    TList *list53 = new TList;
//    list53->Add(Phi_Scattered_475MeV_NegHelCM1_Para);
//    Phi_Scattered_475MeV_NegHelCM1_Perp->Scale(ScaleFactor);
//    list53->Add(Phi_Scattered_475MeV_NegHelCM1_Perp);
//    Phi_Scattered_475MeV_NegHelCM1 = new TH1D( "Phi_Scattered_475MeV_NegHelCM1", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_475MeV_NegHelCM1->Merge(list53);
//
//    TList *list54 = new TList;
//    list54->Add(Phi_Scattered_545MeV_NegHelCM1_Para);
//    Phi_Scattered_545MeV_NegHelCM1_Perp->Scale(ScaleFactor);
//    list54->Add(Phi_Scattered_545MeV_NegHelCM1_Perp);
//    Phi_Scattered_545MeV_NegHelCM1 = new TH1D( "Phi_Scattered_545MeV_NegHelCM1", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_545MeV_NegHelCM1->Merge(list54);
//
//    TList *list55 = new TList;
//    list55->Add(Phi_Scattered_615MeV_NegHelCM1_Para);
//    Phi_Scattered_615MeV_NegHelCM1_Perp->Scale(ScaleFactor);
//    list55->Add(Phi_Scattered_615MeV_NegHelCM1_Perp);
//    Phi_Scattered_615MeV_NegHelCM1 = new TH1D( "Phi_Scattered_615MeV_NegHelCM1", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_615MeV_NegHelCM1->Merge(list55);
//
//    TList *list56 = new TList;
//    list56->Add(Phi_Scattered_685MeV_NegHelCM1_Para);
//    Phi_Scattered_685MeV_NegHelCM1_Perp->Scale(ScaleFactor);
//    list56->Add(Phi_Scattered_685MeV_NegHelCM1_Perp);
//    Phi_Scattered_685MeV_NegHelCM1 = new TH1D( "Phi_Scattered_685MeV_NegHelCM1", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_685MeV_NegHelCM1->Merge(list56);
//
//    TList *list57 = new TList;
//    list57->Add(Phi_Scattered_265MeV_NegHelCM2_Para);
//    Phi_Scattered_265MeV_NegHelCM2_Perp->Scale(ScaleFactor);
//    list57->Add(Phi_Scattered_265MeV_NegHelCM2_Perp);
//    Phi_Scattered_265MeV_NegHelCM2 = new TH1D( "Phi_Scattered_265MeV_NegHelCM2", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.75-0.5) for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_265MeV_NegHelCM2->Merge(list57);
//
//    TList *list58 = new TList;
//    list58->Add(Phi_Scattered_335MeV_NegHelCM2_Para);
//    Phi_Scattered_335MeV_NegHelCM2_Perp->Scale(ScaleFactor);
//    list58->Add(Phi_Scattered_335MeV_NegHelCM2_Perp);
//    Phi_Scattered_335MeV_NegHelCM2 = new TH1D( "Phi_Scattered_335MeV_NegHelCM2", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.75-0.5) for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_335MeV_NegHelCM2->Merge(list58);
//
//    TList *list59 = new TList;
//    list59->Add(Phi_Scattered_405MeV_NegHelCM2_Para);
//    Phi_Scattered_405MeV_NegHelCM2_Perp->Scale(ScaleFactor);
//    list59->Add(Phi_Scattered_405MeV_NegHelCM2_Perp);
//    Phi_Scattered_405MeV_NegHelCM2 = new TH1D( "Phi_Scattered_405MeV_NegHelCM2", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.75-0.5) for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_405MeV_NegHelCM2->Merge(list59);
//
//    TList *list60 = new TList;
//    list60->Add(Phi_Scattered_475MeV_NegHelCM2_Para);
//    Phi_Scattered_475MeV_NegHelCM2_Perp->Scale(ScaleFactor);
//    list60->Add(Phi_Scattered_475MeV_NegHelCM2_Perp);
//    Phi_Scattered_475MeV_NegHelCM2 = new TH1D( "Phi_Scattered_475MeV_NegHelCM2", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.75-0.5) for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_475MeV_NegHelCM2->Merge(list60);
//
//    TList *list61 = new TList;
//    list61->Add(Phi_Scattered_545MeV_NegHelCM2_Para);
//    Phi_Scattered_545MeV_NegHelCM2_Perp->Scale(ScaleFactor);
//    list61->Add(Phi_Scattered_545MeV_NegHelCM2_Perp);
//    Phi_Scattered_545MeV_NegHelCM2 = new TH1D( "Phi_Scattered_545MeV_NegHelCM2", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.75-0.5) for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_545MeV_NegHelCM2->Merge(list61);
//
//    TList *list62 = new TList;
//    list62->Add(Phi_Scattered_615MeV_NegHelCM2_Para);
//    Phi_Scattered_615MeV_NegHelCM2_Perp->Scale(ScaleFactor);
//    list62->Add(Phi_Scattered_615MeV_NegHelCM2_Perp);
//    Phi_Scattered_615MeV_NegHelCM2 = new TH1D( "Phi_Scattered_615MeV_NegHelCM2", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.75-0.5) for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_615MeV_NegHelCM2->Merge(list62);
//
//    TList *list63 = new TList;
//    list63->Add(Phi_Scattered_685MeV_NegHelCM2_Para);
//    Phi_Scattered_685MeV_NegHelCM2_Perp->Scale(ScaleFactor);
//    list63->Add(Phi_Scattered_685MeV_NegHelCM2_Perp);
//    Phi_Scattered_685MeV_NegHelCM2 = new TH1D( "Phi_Scattered_685MeV_NegHelCM2", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.75-0.5) for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_685MeV_NegHelCM2->Merge(list63);
//
//    TList *list64 = new TList;
//    list64->Add(Phi_Scattered_265MeV_NegHelCM3_Para);
//    Phi_Scattered_265MeV_NegHelCM3_Perp->Scale(ScaleFactor);
//    list64->Add(Phi_Scattered_265MeV_NegHelCM3_Perp);
//    Phi_Scattered_265MeV_NegHelCM3 = new TH1D( "Phi_Scattered_265MeV_NegHelCM3", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.5-0.25 for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_265MeV_NegHelCM3->Merge(list64);
//
//    TList *list65 = new TList;
//    list65->Add(Phi_Scattered_335MeV_NegHelCM3_Para);
//    Phi_Scattered_335MeV_NegHelCM3_Perp->Scale(ScaleFactor);
//    list65->Add(Phi_Scattered_335MeV_NegHelCM3_Perp);
//    Phi_Scattered_335MeV_NegHelCM3 = new TH1D( "Phi_Scattered_335MeV_NegHelCM3", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.5-0.25 for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_335MeV_NegHelCM3->Merge(list65);
//
//    TList *list66 = new TList;
//    list66->Add(Phi_Scattered_405MeV_NegHelCM3_Para);
//    Phi_Scattered_405MeV_NegHelCM3_Perp->Scale(ScaleFactor);
//    list66->Add(Phi_Scattered_405MeV_NegHelCM3_Perp);
//    Phi_Scattered_405MeV_NegHelCM3 = new TH1D( "Phi_Scattered_405MeV_NegHelCM3", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.5-0.25 for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_405MeV_NegHelCM3->Merge(list66);
//
//    TList *list67 = new TList;
//    list67->Add(Phi_Scattered_475MeV_NegHelCM3_Para);
//    Phi_Scattered_475MeV_NegHelCM3_Perp->Scale(ScaleFactor);
//    list67->Add(Phi_Scattered_475MeV_NegHelCM3_Perp);
//    Phi_Scattered_475MeV_NegHelCM3 = new TH1D( "Phi_Scattered_475MeV_NegHelCM3", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.5-0.25 for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_475MeV_NegHelCM3->Merge(list67);
//
//    TList *list68 = new TList;
//    list68->Add(Phi_Scattered_545MeV_NegHelCM3_Para);
//    Phi_Scattered_545MeV_NegHelCM3_Perp->Scale(ScaleFactor);
//    list68->Add(Phi_Scattered_545MeV_NegHelCM3_Perp);
//    Phi_Scattered_545MeV_NegHelCM3 = new TH1D( "Phi_Scattered_545MeV_NegHelCM3", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.5-0.25 for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_545MeV_NegHelCM3->Merge(list68);
//
//    TList *list69 = new TList;
//    list69->Add(Phi_Scattered_615MeV_NegHelCM3_Para);
//    Phi_Scattered_615MeV_NegHelCM3_Perp->Scale(ScaleFactor);
//    list69->Add(Phi_Scattered_615MeV_NegHelCM3_Perp);
//    Phi_Scattered_615MeV_NegHelCM3 = new TH1D( "Phi_Scattered_615MeV_NegHelCM3", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.5-0.25 for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_615MeV_NegHelCM3->Merge(list69);
//
//    TList *list70 = new TList;
//    list70->Add(Phi_Scattered_685MeV_NegHelCM3_Para);
//    Phi_Scattered_685MeV_NegHelCM3_Perp->Scale(ScaleFactor);
//    list70->Add(Phi_Scattered_685MeV_NegHelCM3_Perp);
//    Phi_Scattered_685MeV_NegHelCM3 = new TH1D( "Phi_Scattered_685MeV_NegHelCM3", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.5-0.25 for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_685MeV_NegHelCM3->Merge(list70);
//
//    TList *list71 = new TList;
//    list71->Add(Phi_Scattered_265MeV_NegHelCM4_Para);
//    Phi_Scattered_265MeV_NegHelCM4_Perp->Scale(ScaleFactor);
//    list71->Add(Phi_Scattered_265MeV_NegHelCM4_Perp);
//    Phi_Scattered_265MeV_NegHelCM4 = new TH1D( "Phi_Scattered_265MeV_NegHelCM4", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.25-0.0 for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_265MeV_NegHelCM4->Merge(list71);
//
//    TList *list72 = new TList;
//    list72->Add(Phi_Scattered_335MeV_NegHelCM4_Para);
//    Phi_Scattered_335MeV_NegHelCM4_Perp->Scale(ScaleFactor);
//    list72->Add(Phi_Scattered_335MeV_NegHelCM4_Perp);
//    Phi_Scattered_335MeV_NegHelCM4 = new TH1D( "Phi_Scattered_335MeV_NegHelCM4", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.25-0.0 for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_335MeV_NegHelCM4->Merge(list72);
//
//    TList *list73 = new TList;
//    list73->Add(Phi_Scattered_405MeV_NegHelCM4_Para);
//    Phi_Scattered_405MeV_NegHelCM4_Perp->Scale(ScaleFactor);
//    list73->Add(Phi_Scattered_405MeV_NegHelCM4_Perp);
//    Phi_Scattered_405MeV_NegHelCM4 = new TH1D( "Phi_Scattered_405MeV_NegHelCM4", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.25-0.0 for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_405MeV_NegHelCM4->Merge(list73);
//
//    TList *list74 = new TList;
//    list74->Add(Phi_Scattered_475MeV_NegHelCM4_Para);
//    Phi_Scattered_475MeV_NegHelCM4_Perp->Scale(ScaleFactor);
//    list74->Add(Phi_Scattered_475MeV_NegHelCM4_Perp);
//    Phi_Scattered_475MeV_NegHelCM4 = new TH1D( "Phi_Scattered_475MeV_NegHelCM4", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.25-0.0 for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_475MeV_NegHelCM4->Merge(list74);
//
//    TList *list75 = new TList;
//    list75->Add(Phi_Scattered_545MeV_NegHelCM4_Para);
//    Phi_Scattered_545MeV_NegHelCM4_Perp->Scale(ScaleFactor);
//    list75->Add(Phi_Scattered_545MeV_NegHelCM4_Perp);
//    Phi_Scattered_545MeV_NegHelCM4 = new TH1D( "Phi_Scattered_545MeV_NegHelCM4", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.25-0.0 for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_545MeV_NegHelCM4->Merge(list75);
//
//    TList *list76 = new TList;
//    list76->Add(Phi_Scattered_615MeV_NegHelCM4_Para);
//    Phi_Scattered_615MeV_NegHelCM4_Perp->Scale(ScaleFactor);
//    list76->Add(Phi_Scattered_615MeV_NegHelCM4_Perp);
//    Phi_Scattered_615MeV_NegHelCM4 = new TH1D( "Phi_Scattered_615MeV_NegHelCM4", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.25-0.0 for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_615MeV_NegHelCM4->Merge(list76);
//
//    TList *list77 = new TList;
//    list77->Add(Phi_Scattered_685MeV_NegHelCM4_Para);
//    Phi_Scattered_685MeV_NegHelCM4_Perp->Scale(ScaleFactor);
//    list77->Add(Phi_Scattered_685MeV_NegHelCM4_Perp);
//    Phi_Scattered_685MeV_NegHelCM4 = new TH1D( "Phi_Scattered_685MeV_NegHelCM4", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.25-0.0 for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_685MeV_NegHelCM4->Merge(list77);
//
//    TList *list78 = new TList;
//    list78->Add(Phi_Scattered_265MeV_NegHelCM5_Para);
//    Phi_Scattered_265MeV_NegHelCM5_Perp->Scale(ScaleFactor);
//    list78->Add(Phi_Scattered_265MeV_NegHelCM5_Perp);
//    Phi_Scattered_265MeV_NegHelCM5 = new TH1D( "Phi_Scattered_265MeV_NegHelCM5", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25) for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_265MeV_NegHelCM5->Merge(list78);
//
//    TList *list79 = new TList;
//    list79->Add(Phi_Scattered_335MeV_NegHelCM5_Para);
//    Phi_Scattered_335MeV_NegHelCM5_Perp->Scale(ScaleFactor);
//    list79->Add(Phi_Scattered_335MeV_NegHelCM5_Perp);
//    Phi_Scattered_335MeV_NegHelCM5 = new TH1D( "Phi_Scattered_335MeV_NegHelCM5", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25) for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_335MeV_NegHelCM5->Merge(list79);
//
//    TList *list80 = new TList;
//    list80->Add(Phi_Scattered_405MeV_NegHelCM5_Para);
//    Phi_Scattered_405MeV_NegHelCM5_Perp->Scale(ScaleFactor);
//    list80->Add(Phi_Scattered_405MeV_NegHelCM5_Perp);
//    Phi_Scattered_405MeV_NegHelCM5 = new TH1D( "Phi_Scattered_405MeV_NegHelCM5", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25) for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_405MeV_NegHelCM5->Merge(list80);
//
//    TList *list81 = new TList;
//    list81->Add(Phi_Scattered_475MeV_NegHelCM5_Para);
//    Phi_Scattered_475MeV_NegHelCM5_Perp->Scale(ScaleFactor);
//    list81->Add(Phi_Scattered_475MeV_NegHelCM5_Perp);
//    Phi_Scattered_475MeV_NegHelCM5 = new TH1D( "Phi_Scattered_475MeV_NegHelCM5", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25) for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_475MeV_NegHelCM5->Merge(list81);
//
//    TList *list82 = new TList;
//    list82->Add(Phi_Scattered_545MeV_NegHelCM5_Para);
//    Phi_Scattered_545MeV_NegHelCM5_Perp->Scale(ScaleFactor);
//    list82->Add(Phi_Scattered_545MeV_NegHelCM5_Perp);
//    Phi_Scattered_545MeV_NegHelCM5 = new TH1D( "Phi_Scattered_545MeV_NegHelCM5", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25) for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_545MeV_NegHelCM5->Merge(list82);
//
//    TList *list83 = new TList;
//    list83->Add(Phi_Scattered_615MeV_NegHelCM5_Para);
//    Phi_Scattered_615MeV_NegHelCM5_Perp->Scale(ScaleFactor);
//    list83->Add(Phi_Scattered_615MeV_NegHelCM5_Perp);
//    Phi_Scattered_615MeV_NegHelCM5 = new TH1D( "Phi_Scattered_615MeV_NegHelCM5", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25) for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_615MeV_NegHelCM5->Merge(list83);
//
//    TList *list84 = new TList;
//    list84->Add(Phi_Scattered_685MeV_NegHelCM5_Para);
//    Phi_Scattered_685MeV_NegHelCM5_Perp->Scale(ScaleFactor);
//    list84->Add(Phi_Scattered_685MeV_NegHelCM5_Perp);
//    Phi_Scattered_685MeV_NegHelCM5 = new TH1D( "Phi_Scattered_685MeV_NegHelCM5", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25) for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_685MeV_NegHelCM5->Merge(list84);
//
//    TList *list85 = new TList;
//    list85->Add(Phi_Scattered_265MeV_NegHelCM6_Para);
//    Phi_Scattered_265MeV_NegHelCM6_Perp->Scale(ScaleFactor);
//    list85->Add(Phi_Scattered_265MeV_NegHelCM6_Perp);
//    Phi_Scattered_265MeV_NegHelCM6 = new TH1D( "Phi_Scattered_265MeV_NegHelCM6", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5) for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_265MeV_NegHelCM6->Merge(list85);
//
//    TList *list86 = new TList;
//    list86->Add(Phi_Scattered_335MeV_NegHelCM6_Para);
//    Phi_Scattered_335MeV_NegHelCM6_Perp->Scale(ScaleFactor);
//    list86->Add(Phi_Scattered_335MeV_NegHelCM6_Perp);
//    Phi_Scattered_335MeV_NegHelCM6 = new TH1D( "Phi_Scattered_335MeV_NegHelCM6", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5) for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_335MeV_NegHelCM6->Merge(list86);
//
//    TList *list87 = new TList;
//    list87->Add(Phi_Scattered_405MeV_NegHelCM6_Para);
//    Phi_Scattered_405MeV_NegHelCM6_Perp->Scale(ScaleFactor);
//    list87->Add(Phi_Scattered_405MeV_NegHelCM6_Perp);
//    Phi_Scattered_405MeV_NegHelCM6 = new TH1D( "Phi_Scattered_405MeV_NegHelCM6", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5) for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_405MeV_NegHelCM6->Merge(list87);
//
//    TList *list88 = new TList;
//    list88->Add(Phi_Scattered_475MeV_NegHelCM6_Para);
//    Phi_Scattered_475MeV_NegHelCM6_Perp->Scale(ScaleFactor);
//    list88->Add(Phi_Scattered_475MeV_NegHelCM6_Perp);
//    Phi_Scattered_475MeV_NegHelCM6 = new TH1D( "Phi_Scattered_475MeV_NegHelCM6", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5) for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_475MeV_NegHelCM6->Merge(list88);
//
//    TList *list89 = new TList;
//    list89->Add(Phi_Scattered_545MeV_NegHelCM6_Para);
//    Phi_Scattered_545MeV_NegHelCM6_Perp->Scale(ScaleFactor);
//    list89->Add(Phi_Scattered_545MeV_NegHelCM6_Perp);
//    Phi_Scattered_545MeV_NegHelCM6 = new TH1D( "Phi_Scattered_545MeV_NegHelCM6", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5) for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_545MeV_NegHelCM6->Merge(list89);
//
//    TList *list90 = new TList;
//    list90->Add(Phi_Scattered_615MeV_NegHelCM6_Para);
//    Phi_Scattered_615MeV_NegHelCM6_Perp->Scale(ScaleFactor);
//    list90->Add(Phi_Scattered_615MeV_NegHelCM6_Perp);
//    Phi_Scattered_615MeV_NegHelCM6 = new TH1D( "Phi_Scattered_615MeV_NegHelCM6", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5) for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_615MeV_NegHelCM6->Merge(list90);
//
//    TList *list91 = new TList;
//    list91->Add(Phi_Scattered_685MeV_NegHelCM6_Para);
//    Phi_Scattered_685MeV_NegHelCM6_Perp->Scale(ScaleFactor);
//    list91->Add(Phi_Scattered_685MeV_NegHelCM6_Perp);
//    Phi_Scattered_685MeV_NegHelCM6 = new TH1D( "Phi_Scattered_685MeV_NegHelCM6", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5) for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_685MeV_NegHelCM6->Merge(list91);
//
//    TList *list92 = new TList;
//    list92->Add(Phi_Scattered_265MeV_NegHelCM7_Para);
//    Phi_Scattered_265MeV_NegHelCM7_Perp->Scale(ScaleFactor);
//    list92->Add(Phi_Scattered_265MeV_NegHelCM7_Perp);
//    Phi_Scattered_265MeV_NegHelCM7 = new TH1D( "Phi_Scattered_265MeV_NegHelCM7", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75) for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_265MeV_NegHelCM7->Merge(list92);
//
//    TList *list93 = new TList;
//    list93->Add(Phi_Scattered_335MeV_NegHelCM7_Para);
//    Phi_Scattered_335MeV_NegHelCM7_Perp->Scale(ScaleFactor);
//    list93->Add(Phi_Scattered_335MeV_NegHelCM7_Perp);
//    Phi_Scattered_335MeV_NegHelCM7 = new TH1D( "Phi_Scattered_335MeV_NegHelCM7", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75) for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_335MeV_NegHelCM7->Merge(list93);
//
//    TList *list94 = new TList;
//    list94->Add(Phi_Scattered_405MeV_NegHelCM7_Para);
//    Phi_Scattered_405MeV_NegHelCM7_Perp->Scale(ScaleFactor);
//    list94->Add(Phi_Scattered_405MeV_NegHelCM7_Perp);
//    Phi_Scattered_405MeV_NegHelCM7 = new TH1D( "Phi_Scattered_405MeV_NegHelCM7", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75) for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_405MeV_NegHelCM7->Merge(list94);
//
//    TList *list95 = new TList;
//    list95->Add(Phi_Scattered_475MeV_NegHelCM7_Para);
//    Phi_Scattered_475MeV_NegHelCM7_Perp->Scale(ScaleFactor);
//    list95->Add(Phi_Scattered_475MeV_NegHelCM7_Perp);
//    Phi_Scattered_475MeV_NegHelCM7 = new TH1D( "Phi_Scattered_475MeV_NegHelCM7", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75) for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_475MeV_NegHelCM7->Merge(list95);
//
//    TList *list96 = new TList;
//    list96->Add(Phi_Scattered_545MeV_NegHelCM7_Para);
//    Phi_Scattered_545MeV_NegHelCM7_Perp->Scale(ScaleFactor);
//    list96->Add(Phi_Scattered_545MeV_NegHelCM7_Perp);
//    Phi_Scattered_545MeV_NegHelCM7 = new TH1D( "Phi_Scattered_545MeV_NegHelCM7", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75) for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_545MeV_NegHelCM7->Merge(list96);
//
//    TList *list97 = new TList;
//    list97->Add(Phi_Scattered_615MeV_NegHelCM7_Para);
//    Phi_Scattered_615MeV_NegHelCM7_Perp->Scale(ScaleFactor);
//    list97->Add(Phi_Scattered_615MeV_NegHelCM7_Perp);
//    Phi_Scattered_615MeV_NegHelCM7 = new TH1D( "Phi_Scattered_615MeV_NegHelCM7", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75) for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_615MeV_NegHelCM7->Merge(list97);
//
//    TList *list98 = new TList;
//    list98->Add(Phi_Scattered_685MeV_NegHelCM7_Para);
//    Phi_Scattered_685MeV_NegHelCM7_Perp->Scale(ScaleFactor);
//    list98->Add(Phi_Scattered_685MeV_NegHelCM7_Perp);
//    Phi_Scattered_685MeV_NegHelCM7 = new TH1D( "Phi_Scattered_685MeV_NegHelCM7", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75) for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_685MeV_NegHelCM7->Merge(list98);
//
//    TList *list99 = new TList;
//    list99->Add(Phi_Scattered_265MeV_NegHelCM8_Para);
//    Phi_Scattered_265MeV_NegHelCM8_Perp->Scale(ScaleFactor);
//    list99->Add(Phi_Scattered_265MeV_NegHelCM8_Perp);
//    Phi_Scattered_265MeV_NegHelCM8 = new TH1D( "Phi_Scattered_265MeV_NegHelCM8", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}-0.57-(-1.0) for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_265MeV_NegHelCM8->Merge(list99);
//
//    TList *list100 = new TList;
//    list100->Add(Phi_Scattered_335MeV_NegHelCM8_Para);
//    Phi_Scattered_335MeV_NegHelCM8_Perp->Scale(ScaleFactor);
//    list100->Add(Phi_Scattered_335MeV_NegHelCM8_Perp);
//    Phi_Scattered_335MeV_NegHelCM8 = new TH1D( "Phi_Scattered_335MeV_NegHelCM8", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}-0.57-(-1.0) for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_335MeV_NegHelCM8->Merge(list100);
//
//    TList *list101 = new TList;
//    list101->Add(Phi_Scattered_405MeV_NegHelCM8_Para);
//    Phi_Scattered_405MeV_NegHelCM8_Perp->Scale(ScaleFactor);
//    list101->Add(Phi_Scattered_405MeV_NegHelCM8_Perp);
//    Phi_Scattered_405MeV_NegHelCM8 = new TH1D( "Phi_Scattered_405MeV_NegHelCM8", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}-0.57-(-1.0) for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_405MeV_NegHelCM8->Merge(list101);
//
//    TList *list102 = new TList;
//    list102->Add(Phi_Scattered_475MeV_NegHelCM8_Para);
//    Phi_Scattered_475MeV_NegHelCM8_Perp->Scale(ScaleFactor);
//    list102->Add(Phi_Scattered_475MeV_NegHelCM8_Perp);
//    Phi_Scattered_475MeV_NegHelCM8 = new TH1D( "Phi_Scattered_475MeV_NegHelCM8", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}-0.57-(-1.0) for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_475MeV_NegHelCM8->Merge(list102);
//
//    TList *list103 = new TList;
//    list103->Add(Phi_Scattered_545MeV_NegHelCM8_Para);
//    Phi_Scattered_545MeV_NegHelCM8_Perp->Scale(ScaleFactor);
//    list103->Add(Phi_Scattered_545MeV_NegHelCM8_Perp);
//    Phi_Scattered_545MeV_NegHelCM8 = new TH1D( "Phi_Scattered_545MeV_NegHelCM8", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}-0.57-(-1.0) for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_545MeV_NegHelCM8->Merge(list103);
//
//    TList *list104 = new TList;
//    list104->Add(Phi_Scattered_615MeV_NegHelCM8_Para);
//    Phi_Scattered_615MeV_NegHelCM8_Perp->Scale(ScaleFactor);
//    list104->Add(Phi_Scattered_615MeV_NegHelCM8_Perp);
//    Phi_Scattered_615MeV_NegHelCM8 = new TH1D( "Phi_Scattered_615MeV_NegHelCM8", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}-0.57-(-1.0) for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_615MeV_NegHelCM8->Merge(list104);
//
//    TList *list105 = new TList;
//    list105->Add(Phi_Scattered_685MeV_NegHelCM8_Para);
//    Phi_Scattered_685MeV_NegHelCM8_Perp->Scale(ScaleFactor);
//    list105->Add(Phi_Scattered_685MeV_NegHelCM8_Perp);
//    Phi_Scattered_685MeV_NegHelCM8 = new TH1D( "Phi_Scattered_685MeV_NegHelCM8", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}-0.57-(-1.0) for -ve Helicity", 20, -4, 4);
//    Phi_Scattered_685MeV_NegHelCM8->Merge(list105);
//
//    /////////////////////////////////////
//    ////// All Pos Hel PhiSc Dists //////
//    /////////////////////////////////////
//
//    TList *list50a = new TList;
//    list50a->Add(Phi_Scattered_265MeV_PosHelCM1_Para);
//    Phi_Scattered_265MeV_PosHelCM1_Perp->Scale(ScaleFactor);
//    list50a->Add(Phi_Scattered_265MeV_PosHelCM1_Perp);
//    Phi_Scattered_265MeV_PosHelCM1 = new TH1D( "Phi_Scattered_265MeV_PosHelCM1", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_265MeV_PosHelCM1->Merge(list50a);
//
//    TList *list51a = new TList;
//    list51a->Add(Phi_Scattered_335MeV_PosHelCM1_Para);
//    Phi_Scattered_335MeV_PosHelCM1_Perp->Scale(ScaleFactor);
//    list51a->Add(Phi_Scattered_335MeV_PosHelCM1_Perp);
//    Phi_Scattered_335MeV_PosHelCM1 = new TH1D( "Phi_Scattered_335MeV_PosHelCM1", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_335MeV_PosHelCM1->Merge(list51a);
//
//    TList *list52a = new TList;
//    list52a->Add(Phi_Scattered_405MeV_PosHelCM1_Para);
//    Phi_Scattered_405MeV_PosHelCM1_Perp->Scale(ScaleFactor);
//    list52a->Add(Phi_Scattered_405MeV_PosHelCM1_Perp);
//    Phi_Scattered_405MeV_PosHelCM1 = new TH1D( "Phi_Scattered_405MeV_PosHelCM1", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_405MeV_PosHelCM1->Merge(list52a);
//
//    TList *list53a = new TList;
//    list53a->Add(Phi_Scattered_475MeV_PosHelCM1_Para);
//    Phi_Scattered_475MeV_PosHelCM1_Perp->Scale(ScaleFactor);
//    list53a->Add(Phi_Scattered_475MeV_PosHelCM1_Perp);
//    Phi_Scattered_475MeV_PosHelCM1 = new TH1D( "Phi_Scattered_475MeV_PosHelCM1", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_475MeV_PosHelCM1->Merge(list53a);
//
//    TList *list54a = new TList;
//    list54a->Add(Phi_Scattered_545MeV_PosHelCM1_Para);
//    Phi_Scattered_545MeV_PosHelCM1_Perp->Scale(ScaleFactor);
//    list54a->Add(Phi_Scattered_545MeV_PosHelCM1_Perp);
//    Phi_Scattered_545MeV_PosHelCM1 = new TH1D( "Phi_Scattered_545MeV_PosHelCM1", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_545MeV_PosHelCM1->Merge(list54a);
//
//    TList *list55a = new TList;
//    list55a->Add(Phi_Scattered_615MeV_PosHelCM1_Para);
//    Phi_Scattered_615MeV_PosHelCM1_Perp->Scale(ScaleFactor);
//    list55a->Add(Phi_Scattered_615MeV_PosHelCM1_Perp);
//    Phi_Scattered_615MeV_PosHelCM1 = new TH1D( "Phi_Scattered_615MeV_PosHelCM1", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_615MeV_PosHelCM1->Merge(list55a);
//
//    TList *list56a = new TList;
//    list56a->Add(Phi_Scattered_685MeV_PosHelCM1_Para);
//    Phi_Scattered_685MeV_PosHelCM1_Perp->Scale(ScaleFactor);
//    list56a->Add(Phi_Scattered_685MeV_PosHelCM1_Perp);
//    Phi_Scattered_685MeV_PosHelCM1 = new TH1D( "Phi_Scattered_685MeV_PosHelCM1", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_685MeV_PosHelCM1->Merge(list56a);
//
//    TList *list57a = new TList;
//    list57a->Add(Phi_Scattered_265MeV_PosHelCM2_Para);
//    Phi_Scattered_265MeV_PosHelCM2_Perp->Scale(ScaleFactor);
//    list57a->Add(Phi_Scattered_265MeV_PosHelCM2_Perp);
//    Phi_Scattered_265MeV_PosHelCM2 = new TH1D( "Phi_Scattered_265MeV_PosHelCM2", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.75-0.5) for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_265MeV_PosHelCM2->Merge(list57a);
//
//    TList *list58a = new TList;
//    list58a->Add(Phi_Scattered_335MeV_PosHelCM2_Para);
//    Phi_Scattered_335MeV_PosHelCM2_Perp->Scale(ScaleFactor);
//    list58a->Add(Phi_Scattered_335MeV_PosHelCM2_Perp);
//    Phi_Scattered_335MeV_PosHelCM2 = new TH1D( "Phi_Scattered_335MeV_PosHelCM2", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.75-0.5) for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_335MeV_PosHelCM2->Merge(list58a);
//
//    TList *list59a = new TList;
//    list59a->Add(Phi_Scattered_405MeV_PosHelCM2_Para);
//    Phi_Scattered_405MeV_PosHelCM2_Perp->Scale(ScaleFactor);
//    list59a->Add(Phi_Scattered_405MeV_PosHelCM2_Perp);
//    Phi_Scattered_405MeV_PosHelCM2 = new TH1D( "Phi_Scattered_405MeV_PosHelCM2", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.75-0.5) for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_405MeV_PosHelCM2->Merge(list59a);
//
//    TList *list60a = new TList;
//    list60a->Add(Phi_Scattered_475MeV_PosHelCM2_Para);
//    Phi_Scattered_475MeV_PosHelCM2_Perp->Scale(ScaleFactor);
//    list60a->Add(Phi_Scattered_475MeV_PosHelCM2_Perp);
//    Phi_Scattered_475MeV_PosHelCM2 = new TH1D( "Phi_Scattered_475MeV_PosHelCM2", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.75-0.5) for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_475MeV_PosHelCM2->Merge(list60a);
//
//    TList *list61a = new TList;
//    list61a->Add(Phi_Scattered_545MeV_PosHelCM2_Para);
//    Phi_Scattered_545MeV_PosHelCM2_Perp->Scale(ScaleFactor);
//    list61a->Add(Phi_Scattered_545MeV_PosHelCM2_Perp);
//    Phi_Scattered_545MeV_PosHelCM2 = new TH1D( "Phi_Scattered_545MeV_PosHelCM2", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.75-0.5) for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_545MeV_PosHelCM2->Merge(list61a);
//
//    TList *list62a = new TList;
//    list62a->Add(Phi_Scattered_615MeV_PosHelCM2_Para);
//    Phi_Scattered_615MeV_PosHelCM2_Perp->Scale(ScaleFactor);
//    list62a->Add(Phi_Scattered_615MeV_PosHelCM2_Perp);
//    Phi_Scattered_615MeV_PosHelCM2 = new TH1D( "Phi_Scattered_615MeV_PosHelCM2", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.75-0.5) for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_615MeV_PosHelCM2->Merge(list62a);
//
//    TList *list63a = new TList;
//    list63a->Add(Phi_Scattered_685MeV_PosHelCM2_Para);
//    Phi_Scattered_685MeV_PosHelCM2_Perp->Scale(ScaleFactor);
//    list63a->Add(Phi_Scattered_685MeV_PosHelCM2_Perp);
//    Phi_Scattered_685MeV_PosHelCM2 = new TH1D( "Phi_Scattered_685MeV_PosHelCM2", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.75-0.5) for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_685MeV_PosHelCM2->Merge(list63a);
//
//    TList *list64a = new TList;
//    list64a->Add(Phi_Scattered_265MeV_PosHelCM3_Para);
//    Phi_Scattered_265MeV_PosHelCM3_Perp->Scale(ScaleFactor);
//    list64a->Add(Phi_Scattered_265MeV_PosHelCM3_Perp);
//    Phi_Scattered_265MeV_PosHelCM3 = new TH1D( "Phi_Scattered_265MeV_PosHelCM3", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.5-0.25 for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_265MeV_PosHelCM3->Merge(list64a);
//
//    TList *list65a = new TList;
//    list65a->Add(Phi_Scattered_335MeV_PosHelCM3_Para);
//    Phi_Scattered_335MeV_PosHelCM3_Perp->Scale(ScaleFactor);
//    list65a->Add(Phi_Scattered_335MeV_PosHelCM3_Perp);
//    Phi_Scattered_335MeV_PosHelCM3 = new TH1D( "Phi_Scattered_335MeV_PosHelCM3", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.5-0.25 for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_335MeV_PosHelCM3->Merge(list65a);
//
//    TList *list66a = new TList;
//    list66a->Add(Phi_Scattered_405MeV_PosHelCM3_Para);
//    Phi_Scattered_405MeV_PosHelCM3_Perp->Scale(ScaleFactor);
//    list66a->Add(Phi_Scattered_405MeV_PosHelCM3_Perp);
//    Phi_Scattered_405MeV_PosHelCM3 = new TH1D( "Phi_Scattered_405MeV_PosHelCM3", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.5-0.25 for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_405MeV_PosHelCM3->Merge(list66a);
//
//    TList *list67a = new TList;
//    list67a->Add(Phi_Scattered_475MeV_PosHelCM3_Para);
//    Phi_Scattered_475MeV_PosHelCM3_Perp->Scale(ScaleFactor);
//    list67a->Add(Phi_Scattered_475MeV_PosHelCM3_Perp);
//    Phi_Scattered_475MeV_PosHelCM3 = new TH1D( "Phi_Scattered_475MeV_PosHelCM3", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.5-0.25 for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_475MeV_PosHelCM3->Merge(list67a);
//
//    TList *list68a = new TList;
//    list68a->Add(Phi_Scattered_545MeV_PosHelCM3_Para);
//    Phi_Scattered_545MeV_PosHelCM3_Perp->Scale(ScaleFactor);
//    list68a->Add(Phi_Scattered_545MeV_PosHelCM3_Perp);
//    Phi_Scattered_545MeV_PosHelCM3 = new TH1D( "Phi_Scattered_545MeV_PosHelCM3", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.5-0.25 for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_545MeV_PosHelCM3->Merge(list68a);
//
//    TList *list69a = new TList;
//    list69a->Add(Phi_Scattered_615MeV_PosHelCM3_Para);
//    Phi_Scattered_615MeV_PosHelCM3_Perp->Scale(ScaleFactor);
//    list69a->Add(Phi_Scattered_615MeV_PosHelCM3_Perp);
//    Phi_Scattered_615MeV_PosHelCM3 = new TH1D( "Phi_Scattered_615MeV_PosHelCM3", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.5-0.25 for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_615MeV_PosHelCM3->Merge(list69a);
//
//    TList *list70a = new TList;
//    list70a->Add(Phi_Scattered_685MeV_PosHelCM3_Para);
//    Phi_Scattered_685MeV_PosHelCM3_Perp->Scale(ScaleFactor);
//    list70a->Add(Phi_Scattered_685MeV_PosHelCM3_Perp);
//    Phi_Scattered_685MeV_PosHelCM3 = new TH1D( "Phi_Scattered_685MeV_PosHelCM3", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.5-0.25 for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_685MeV_PosHelCM3->Merge(list70a);
//
//    TList *list71a = new TList;
//    list71a->Add(Phi_Scattered_265MeV_PosHelCM4_Para);
//    Phi_Scattered_265MeV_PosHelCM4_Perp->Scale(ScaleFactor);
//    list71a->Add(Phi_Scattered_265MeV_PosHelCM4_Perp);
//    Phi_Scattered_265MeV_PosHelCM4 = new TH1D( "Phi_Scattered_265MeV_PosHelCM4", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.25-0.0 for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_265MeV_PosHelCM4->Merge(list71a);
//
//    TList *list72a = new TList;
//    list72a->Add(Phi_Scattered_335MeV_PosHelCM4_Para);
//    Phi_Scattered_335MeV_PosHelCM4_Perp->Scale(ScaleFactor);
//    list72a->Add(Phi_Scattered_335MeV_PosHelCM4_Perp);
//    Phi_Scattered_335MeV_PosHelCM4 = new TH1D( "Phi_Scattered_335MeV_PosHelCM4", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.25-0.0 for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_335MeV_PosHelCM4->Merge(list72a);
//
//    TList *list73a = new TList;
//    list73a->Add(Phi_Scattered_405MeV_PosHelCM4_Para);
//    Phi_Scattered_405MeV_PosHelCM4_Perp->Scale(ScaleFactor);
//    list73a->Add(Phi_Scattered_405MeV_PosHelCM4_Perp);
//    Phi_Scattered_405MeV_PosHelCM4 = new TH1D( "Phi_Scattered_405MeV_PosHelCM4", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.25-0.0 for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_405MeV_PosHelCM4->Merge(list73a);
//
//    TList *list74a = new TList;
//    list74a->Add(Phi_Scattered_475MeV_PosHelCM4_Para);
//    Phi_Scattered_475MeV_PosHelCM4_Perp->Scale(ScaleFactor);
//    list74a->Add(Phi_Scattered_475MeV_PosHelCM4_Perp);
//    Phi_Scattered_475MeV_PosHelCM4 = new TH1D( "Phi_Scattered_475MeV_PosHelCM4", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.25-0.0 for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_475MeV_PosHelCM4->Merge(list74a);
//
//    TList *list75a = new TList;
//    list75a->Add(Phi_Scattered_545MeV_PosHelCM4_Para);
//    Phi_Scattered_545MeV_PosHelCM4_Perp->Scale(ScaleFactor);
//    list75a->Add(Phi_Scattered_545MeV_PosHelCM4_Perp);
//    Phi_Scattered_545MeV_PosHelCM4 = new TH1D( "Phi_Scattered_545MeV_PosHelCM4", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.25-0.0 for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_545MeV_PosHelCM4->Merge(list75a);
//
//    TList *list76a = new TList;
//    list76a->Add(Phi_Scattered_615MeV_PosHelCM4_Para);
//    Phi_Scattered_615MeV_PosHelCM4_Perp->Scale(ScaleFactor);
//    list76a->Add(Phi_Scattered_615MeV_PosHelCM4_Perp);
//    Phi_Scattered_615MeV_PosHelCM4 = new TH1D( "Phi_Scattered_615MeV_PosHelCM4", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.25-0.0 for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_615MeV_PosHelCM4->Merge(list76a);
//
//    TList *list77a = new TList;
//    list77a->Add(Phi_Scattered_685MeV_PosHelCM4_Para);
//    Phi_Scattered_685MeV_PosHelCM4_Perp->Scale(ScaleFactor);
//    list77a->Add(Phi_Scattered_685MeV_PosHelCM4_Perp);
//    Phi_Scattered_685MeV_PosHelCM4 = new TH1D( "Phi_Scattered_685MeV_PosHelCM4", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.25-0.0 for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_685MeV_PosHelCM4->Merge(list77a);
//
//    TList *list78a = new TList;
//    list78a->Add(Phi_Scattered_265MeV_PosHelCM5_Para);
//    Phi_Scattered_265MeV_PosHelCM5_Perp->Scale(ScaleFactor);
//    list78a->Add(Phi_Scattered_265MeV_PosHelCM5_Perp);
//    Phi_Scattered_265MeV_PosHelCM5 = new TH1D( "Phi_Scattered_265MeV_PosHelCM5", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25) for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_265MeV_PosHelCM5->Merge(list78a);
//
//    TList *list79a = new TList;
//    list79a->Add(Phi_Scattered_335MeV_PosHelCM5_Para);
//    Phi_Scattered_335MeV_PosHelCM5_Perp->Scale(ScaleFactor);
//    list79a->Add(Phi_Scattered_335MeV_PosHelCM5_Perp);
//    Phi_Scattered_335MeV_PosHelCM5 = new TH1D( "Phi_Scattered_335MeV_PosHelCM5", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25) for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_335MeV_PosHelCM5->Merge(list79a);
//
//    TList *list80a = new TList;
//    list80a->Add(Phi_Scattered_405MeV_PosHelCM5_Para);
//    Phi_Scattered_405MeV_PosHelCM5_Perp->Scale(ScaleFactor);
//    list80a->Add(Phi_Scattered_405MeV_PosHelCM5_Perp);
//    Phi_Scattered_405MeV_PosHelCM5 = new TH1D( "Phi_Scattered_405MeV_PosHelCM5", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25) for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_405MeV_PosHelCM5->Merge(list80a);
//
//    TList *list81a = new TList;
//    list81a->Add(Phi_Scattered_475MeV_PosHelCM5_Para);
//    Phi_Scattered_475MeV_PosHelCM5_Perp->Scale(ScaleFactor);
//    list81a->Add(Phi_Scattered_475MeV_PosHelCM5_Perp);
//    Phi_Scattered_475MeV_PosHelCM5 = new TH1D( "Phi_Scattered_475MeV_PosHelCM5", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25) for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_475MeV_PosHelCM5->Merge(list81a);
//
//    TList *list82a = new TList;
//    list82a->Add(Phi_Scattered_545MeV_PosHelCM5_Para);
//    Phi_Scattered_545MeV_PosHelCM5_Perp->Scale(ScaleFactor);
//    list82a->Add(Phi_Scattered_545MeV_PosHelCM5_Perp);
//    Phi_Scattered_545MeV_PosHelCM5 = new TH1D( "Phi_Scattered_545MeV_PosHelCM5", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25) for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_545MeV_PosHelCM5->Merge(list82a);
//
//    TList *list83a = new TList;
//    list83a->Add(Phi_Scattered_615MeV_PosHelCM5_Para);
//    Phi_Scattered_615MeV_PosHelCM5_Perp->Scale(ScaleFactor);
//    list83a->Add(Phi_Scattered_615MeV_PosHelCM5_Perp);
//    Phi_Scattered_615MeV_PosHelCM5 = new TH1D( "Phi_Scattered_615MeV_PosHelCM5", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25) for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_615MeV_PosHelCM5->Merge(list83a);
//
//    TList *list84a = new TList;
//    list84a->Add(Phi_Scattered_685MeV_PosHelCM5_Para);
//    Phi_Scattered_685MeV_PosHelCM5_Perp->Scale(ScaleFactor);
//    list84a->Add(Phi_Scattered_685MeV_PosHelCM5_Perp);
//    Phi_Scattered_685MeV_PosHelCM5 = new TH1D( "Phi_Scattered_685MeV_PosHelCM5", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25) for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_685MeV_PosHelCM5->Merge(list84a);
//
//    TList *list85a = new TList;
//    list85a->Add(Phi_Scattered_265MeV_PosHelCM6_Para);
//    Phi_Scattered_265MeV_PosHelCM6_Perp->Scale(ScaleFactor);
//    list85a->Add(Phi_Scattered_265MeV_PosHelCM6_Perp);
//    Phi_Scattered_265MeV_PosHelCM6 = new TH1D( "Phi_Scattered_265MeV_PosHelCM6", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5) for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_265MeV_PosHelCM6->Merge(list85a);
//
//    TList *list86a = new TList;
//    list86a->Add(Phi_Scattered_335MeV_PosHelCM6_Para);
//    Phi_Scattered_335MeV_PosHelCM6_Perp->Scale(ScaleFactor);
//    list86a->Add(Phi_Scattered_335MeV_PosHelCM6_Perp);
//    Phi_Scattered_335MeV_PosHelCM6 = new TH1D( "Phi_Scattered_335MeV_PosHelCM6", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5) for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_335MeV_PosHelCM6->Merge(list86a);
//
//    TList *list87a = new TList;
//    list87a->Add(Phi_Scattered_405MeV_PosHelCM6_Para);
//    Phi_Scattered_405MeV_PosHelCM6_Perp->Scale(ScaleFactor);
//    list87a->Add(Phi_Scattered_405MeV_PosHelCM6_Perp);
//    Phi_Scattered_405MeV_PosHelCM6 = new TH1D( "Phi_Scattered_405MeV_PosHelCM6", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5) for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_405MeV_PosHelCM6->Merge(list87a);
//
//    TList *list88a = new TList;
//    list88a->Add(Phi_Scattered_475MeV_PosHelCM6_Para);
//    Phi_Scattered_475MeV_PosHelCM6_Perp->Scale(ScaleFactor);
//    list88a->Add(Phi_Scattered_475MeV_PosHelCM6_Perp);
//    Phi_Scattered_475MeV_PosHelCM6 = new TH1D( "Phi_Scattered_475MeV_PosHelCM6", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5) for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_475MeV_PosHelCM6->Merge(list88a);
//
//    TList *list89a = new TList;
//    list89a->Add(Phi_Scattered_545MeV_PosHelCM6_Para);
//    Phi_Scattered_545MeV_PosHelCM6_Perp->Scale(ScaleFactor);
//    list89a->Add(Phi_Scattered_545MeV_PosHelCM6_Perp);
//    Phi_Scattered_545MeV_PosHelCM6 = new TH1D( "Phi_Scattered_545MeV_PosHelCM6", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5) for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_545MeV_PosHelCM6->Merge(list89a);
//
//    TList *list90a = new TList;
//    list90a->Add(Phi_Scattered_615MeV_PosHelCM6_Para);
//    Phi_Scattered_615MeV_PosHelCM6_Perp->Scale(ScaleFactor);
//    list90a->Add(Phi_Scattered_615MeV_PosHelCM6_Perp);
//    Phi_Scattered_615MeV_PosHelCM6 = new TH1D( "Phi_Scattered_615MeV_PosHelCM6", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5) for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_615MeV_PosHelCM6->Merge(list90a);
//
//    TList *list91a = new TList;
//    list91a->Add(Phi_Scattered_685MeV_PosHelCM6_Para);
//    Phi_Scattered_685MeV_PosHelCM6_Perp->Scale(ScaleFactor);
//    list91a->Add(Phi_Scattered_685MeV_PosHelCM6_Perp);
//    Phi_Scattered_685MeV_PosHelCM6 = new TH1D( "Phi_Scattered_685MeV_PosHelCM6", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5) for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_685MeV_PosHelCM6->Merge(list91a);
//
//    TList *list92a = new TList;
//    list92a->Add(Phi_Scattered_265MeV_PosHelCM7_Para);
//    Phi_Scattered_265MeV_PosHelCM7_Perp->Scale(ScaleFactor);
//    list92a->Add(Phi_Scattered_265MeV_PosHelCM7_Perp);
//    Phi_Scattered_265MeV_PosHelCM7 = new TH1D( "Phi_Scattered_265MeV_PosHelCM7", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75) for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_265MeV_PosHelCM7->Merge(list92a);
//
//    TList *list93a = new TList;
//    list93a->Add(Phi_Scattered_335MeV_PosHelCM7_Para);
//    Phi_Scattered_335MeV_PosHelCM7_Perp->Scale(ScaleFactor);
//    list93a->Add(Phi_Scattered_335MeV_PosHelCM7_Perp);
//    Phi_Scattered_335MeV_PosHelCM7 = new TH1D( "Phi_Scattered_335MeV_PosHelCM7", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75) for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_335MeV_PosHelCM7->Merge(list93a);
//
//    TList *list94a = new TList;
//    list94a->Add(Phi_Scattered_405MeV_PosHelCM7_Para);
//    Phi_Scattered_405MeV_PosHelCM7_Perp->Scale(ScaleFactor);
//    list94a->Add(Phi_Scattered_405MeV_PosHelCM7_Perp);
//    Phi_Scattered_405MeV_PosHelCM7 = new TH1D( "Phi_Scattered_405MeV_PosHelCM7", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75) for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_405MeV_PosHelCM7->Merge(list94a);
//
//    TList *list95a = new TList;
//    list95a->Add(Phi_Scattered_475MeV_PosHelCM7_Para);
//    Phi_Scattered_475MeV_PosHelCM7_Perp->Scale(ScaleFactor);
//    list95a->Add(Phi_Scattered_475MeV_PosHelCM7_Perp);
//    Phi_Scattered_475MeV_PosHelCM7 = new TH1D( "Phi_Scattered_475MeV_PosHelCM7", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75) for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_475MeV_PosHelCM7->Merge(list95a);
//
//    TList *list96a = new TList;
//    list96a->Add(Phi_Scattered_545MeV_PosHelCM7_Para);
//    Phi_Scattered_545MeV_PosHelCM7_Perp->Scale(ScaleFactor);
//    list96a->Add(Phi_Scattered_545MeV_PosHelCM7_Perp);
//    Phi_Scattered_545MeV_PosHelCM7 = new TH1D( "Phi_Scattered_545MeV_PosHelCM7", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75) for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_545MeV_PosHelCM7->Merge(list96a);
//
//    TList *list97a = new TList;
//    list97a->Add(Phi_Scattered_615MeV_PosHelCM7_Para);
//    Phi_Scattered_615MeV_PosHelCM7_Perp->Scale(ScaleFactor);
//    list97a->Add(Phi_Scattered_615MeV_PosHelCM7_Perp);
//    Phi_Scattered_615MeV_PosHelCM7 = new TH1D( "Phi_Scattered_615MeV_PosHelCM7", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75) for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_615MeV_PosHelCM7->Merge(list97a);
//
//    TList *list98a = new TList;
//    list98a->Add(Phi_Scattered_685MeV_PosHelCM7_Para);
//    Phi_Scattered_685MeV_PosHelCM7_Perp->Scale(ScaleFactor);
//    list98a->Add(Phi_Scattered_685MeV_PosHelCM7_Perp);
//    Phi_Scattered_685MeV_PosHelCM7 = new TH1D( "Phi_Scattered_685MeV_PosHelCM7", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75) for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_685MeV_PosHelCM7->Merge(list98a);
//
//    TList *list99a = new TList;
//    list99a->Add(Phi_Scattered_265MeV_PosHelCM8_Para);
//    Phi_Scattered_265MeV_PosHelCM8_Perp->Scale(ScaleFactor);
//    list99a->Add(Phi_Scattered_265MeV_PosHelCM8_Perp);
//    Phi_Scattered_265MeV_PosHelCM8 = new TH1D( "Phi_Scattered_265MeV_PosHelCM8", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}-0.57-(-1.0) for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_265MeV_PosHelCM8->Merge(list99a);
//
//    TList *list100a = new TList;
//    list100a->Add(Phi_Scattered_335MeV_PosHelCM8_Para);
//    Phi_Scattered_335MeV_PosHelCM8_Perp->Scale(ScaleFactor);
//    list100a->Add(Phi_Scattered_335MeV_PosHelCM8_Perp);
//    Phi_Scattered_335MeV_PosHelCM8 = new TH1D( "Phi_Scattered_335MeV_PosHelCM8", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}-0.57-(-1.0) for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_335MeV_PosHelCM8->Merge(list100a);
//
//    TList *list101a = new TList;
//    list101a->Add(Phi_Scattered_405MeV_PosHelCM8_Para);
//    Phi_Scattered_405MeV_PosHelCM8_Perp->Scale(ScaleFactor);
//    list101a->Add(Phi_Scattered_405MeV_PosHelCM8_Perp);
//    Phi_Scattered_405MeV_PosHelCM8 = new TH1D( "Phi_Scattered_405MeV_PosHelCM8", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}-0.57-(-1.0) for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_405MeV_PosHelCM8->Merge(list101a);
//
//    TList *list102a = new TList;
//    list102a->Add(Phi_Scattered_475MeV_PosHelCM8_Para);
//    Phi_Scattered_475MeV_PosHelCM8_Perp->Scale(ScaleFactor);
//    list102a->Add(Phi_Scattered_475MeV_PosHelCM8_Perp);
//    Phi_Scattered_475MeV_PosHelCM8 = new TH1D( "Phi_Scattered_475MeV_PosHelCM8", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}-0.57-(-1.0) for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_475MeV_PosHelCM8->Merge(list102a);
//
//    TList *list103a = new TList;
//    list103a->Add(Phi_Scattered_545MeV_PosHelCM8_Para);
//    Phi_Scattered_545MeV_PosHelCM8_Perp->Scale(ScaleFactor);
//    list103a->Add(Phi_Scattered_545MeV_PosHelCM8_Perp);
//    Phi_Scattered_545MeV_PosHelCM8 = new TH1D( "Phi_Scattered_545MeV_PosHelCM8", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}-0.57-(-1.0) for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_545MeV_PosHelCM8->Merge(list103a);
//
//    TList *list104a = new TList;
//    list104a->Add(Phi_Scattered_615MeV_PosHelCM8_Para);
//    Phi_Scattered_615MeV_PosHelCM8_Perp->Scale(ScaleFactor);
//    list104a->Add(Phi_Scattered_615MeV_PosHelCM8_Perp);
//    Phi_Scattered_615MeV_PosHelCM8 = new TH1D( "Phi_Scattered_615MeV_PosHelCM8", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}-0.57-(-1.0) for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_615MeV_PosHelCM8->Merge(list104a);
//
//    TList *list105a = new TList;
//    list105a->Add(Phi_Scattered_685MeV_PosHelCM8_Para);
//    Phi_Scattered_685MeV_PosHelCM8_Perp->Scale(ScaleFactor);
//    list105a->Add(Phi_Scattered_685MeV_PosHelCM8_Perp);
//    Phi_Scattered_685MeV_PosHelCM8 = new TH1D( "Phi_Scattered_685MeV_PosHelCM8", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}-0.57-(-1.0) for +ve Helicity", 20, -4, 4);
//    Phi_Scattered_685MeV_PosHelCM8->Merge(list105a);

    TFile f2("ParaPerp_Total_28_Combined_Unpolarised.root", "RECREATE");

//    time->Write();
//    time_cut->Write();
//    Eg->Write();
//
//    PhiDet->Write();
//    PhiRec->Write();
//    Theta_Scattered->Write();
//    Phi_Scattered->Write();
//    MMpEpCorrected->Write();
//    ZpDist->Write();
//    ThetanDist->Write();
//    ThetanCMDist->Write();
//    E_dE->Write();
//    ThetaScPhiSc->Write();
//    EdEMWPCp->Write();
//    EdEMWPCn->Write();
//    ClosestApproach->Write();
//    POCAr->Write();
//    ScatterVertexZ->Write();
//    ScatterVertexZr->Write();
//    ScatterVertexXY->Write();
//    ScatterVertex->Write();
//    POCArPhiSc->Write();
//    ThetapCorrDiff->Write();
//    PhipCorrDiff->Write();
//    ThetaDiff->Write();
//    PhiDiff->Write();
//    PhiDiffThetaDiff->Write();
//    PhiScEg->Write();
//    PhiScEp->Write();
//    PhiScThetan->Write();
//    EMWPCnPhiSc->Write();

    PhiSc320->Write();
    PhiSc360->Write();
    PhiSc400->Write();
    PhiSc440->Write();
    PhiSc480->Write();
    PhiSc520->Write();
    PhiSc560->Write();
    PhiSc600->Write();
    PhiSc640->Write();
    PhiSc680->Write();
//
//    MMp200300->Write();
//    MMp300400->Write();
//    MMp400500->Write();
//    MMp500600->Write();
//    MMp600700->Write();
//    MMp700800->Write();
//    MMp800900->Write();
//
//    Phi_Scattered_265MeV_NegHelCM1->Write();
//    Phi_Scattered_335MeV_NegHelCM1->Write();
//    Phi_Scattered_405MeV_NegHelCM1->Write();
//    Phi_Scattered_475MeV_NegHelCM1->Write();
//    Phi_Scattered_545MeV_NegHelCM1->Write();
//    Phi_Scattered_615MeV_NegHelCM1->Write();
//    Phi_Scattered_685MeV_NegHelCM1->Write();
//
//    Phi_Scattered_265MeV_NegHelCM2->Write();
//    Phi_Scattered_335MeV_NegHelCM2->Write();
//    Phi_Scattered_405MeV_NegHelCM2->Write();
//    Phi_Scattered_475MeV_NegHelCM2->Write();
//    Phi_Scattered_545MeV_NegHelCM2->Write();
//    Phi_Scattered_615MeV_NegHelCM2->Write();
//    Phi_Scattered_685MeV_NegHelCM2->Write();
//
//    Phi_Scattered_265MeV_NegHelCM3->Write();
//    Phi_Scattered_335MeV_NegHelCM3->Write();
//    Phi_Scattered_405MeV_NegHelCM3->Write();
//    Phi_Scattered_475MeV_NegHelCM3->Write();
//    Phi_Scattered_545MeV_NegHelCM3->Write();
//    Phi_Scattered_615MeV_NegHelCM3->Write();
//    Phi_Scattered_685MeV_NegHelCM3->Write();
//
//    Phi_Scattered_265MeV_NegHelCM4->Write();
//    Phi_Scattered_335MeV_NegHelCM4->Write();
//    Phi_Scattered_405MeV_NegHelCM4->Write();
//    Phi_Scattered_475MeV_NegHelCM4->Write();
//    Phi_Scattered_545MeV_NegHelCM4->Write();
//    Phi_Scattered_615MeV_NegHelCM4->Write();
//    Phi_Scattered_685MeV_NegHelCM4->Write();
//
//    Phi_Scattered_265MeV_NegHelCM5->Write();
//    Phi_Scattered_335MeV_NegHelCM5->Write();
//    Phi_Scattered_405MeV_NegHelCM5->Write();
//    Phi_Scattered_475MeV_NegHelCM5->Write();
//    Phi_Scattered_545MeV_NegHelCM5->Write();
//    Phi_Scattered_615MeV_NegHelCM5->Write();
//    Phi_Scattered_685MeV_NegHelCM5->Write();
//
//    Phi_Scattered_265MeV_NegHelCM6->Write();
//    Phi_Scattered_335MeV_NegHelCM6->Write();
//    Phi_Scattered_405MeV_NegHelCM6->Write();
//    Phi_Scattered_475MeV_NegHelCM6->Write();
//    Phi_Scattered_545MeV_NegHelCM6->Write();
//    Phi_Scattered_615MeV_NegHelCM6->Write();
//    Phi_Scattered_685MeV_NegHelCM6->Write();
//
//    Phi_Scattered_265MeV_NegHelCM7->Write();
//    Phi_Scattered_335MeV_NegHelCM7->Write();
//    Phi_Scattered_405MeV_NegHelCM7->Write();
//    Phi_Scattered_475MeV_NegHelCM7->Write();
//    Phi_Scattered_545MeV_NegHelCM7->Write();
//    Phi_Scattered_615MeV_NegHelCM7->Write();
//    Phi_Scattered_685MeV_NegHelCM7->Write();
//
//    Phi_Scattered_265MeV_NegHelCM8->Write();
//    Phi_Scattered_335MeV_NegHelCM8->Write();
//    Phi_Scattered_405MeV_NegHelCM8->Write();
//    Phi_Scattered_475MeV_NegHelCM8->Write();
//    Phi_Scattered_545MeV_NegHelCM8->Write();
//    Phi_Scattered_615MeV_NegHelCM8->Write();
//    Phi_Scattered_685MeV_NegHelCM8->Write();
//
//    Phi_Scattered_265MeV_PosHelCM1->Write();
//    Phi_Scattered_335MeV_PosHelCM1->Write();
//    Phi_Scattered_405MeV_PosHelCM1->Write();
//    Phi_Scattered_475MeV_PosHelCM1->Write();
//    Phi_Scattered_545MeV_PosHelCM1->Write();
//    Phi_Scattered_615MeV_PosHelCM1->Write();
//    Phi_Scattered_685MeV_PosHelCM1->Write();
//
//    Phi_Scattered_265MeV_PosHelCM2->Write();
//    Phi_Scattered_335MeV_PosHelCM2->Write();
//    Phi_Scattered_405MeV_PosHelCM2->Write();
//    Phi_Scattered_475MeV_PosHelCM2->Write();
//    Phi_Scattered_545MeV_PosHelCM2->Write();
//    Phi_Scattered_615MeV_PosHelCM2->Write();
//    Phi_Scattered_685MeV_PosHelCM2->Write();
//
//    Phi_Scattered_265MeV_PosHelCM3->Write();
//    Phi_Scattered_335MeV_PosHelCM3->Write();
//    Phi_Scattered_405MeV_PosHelCM3->Write();
//    Phi_Scattered_475MeV_PosHelCM3->Write();
//    Phi_Scattered_545MeV_PosHelCM3->Write();
//    Phi_Scattered_615MeV_PosHelCM3->Write();
//    Phi_Scattered_685MeV_PosHelCM3->Write();
//
//    Phi_Scattered_265MeV_PosHelCM4->Write();
//    Phi_Scattered_335MeV_PosHelCM4->Write();
//    Phi_Scattered_405MeV_PosHelCM4->Write();
//    Phi_Scattered_475MeV_PosHelCM4->Write();
//    Phi_Scattered_545MeV_PosHelCM4->Write();
//    Phi_Scattered_615MeV_PosHelCM4->Write();
//    Phi_Scattered_685MeV_PosHelCM4->Write();
//
//    Phi_Scattered_265MeV_PosHelCM5->Write();
//    Phi_Scattered_335MeV_PosHelCM5->Write();
//    Phi_Scattered_405MeV_PosHelCM5->Write();
//    Phi_Scattered_475MeV_PosHelCM5->Write();
//    Phi_Scattered_545MeV_PosHelCM5->Write();
//    Phi_Scattered_615MeV_PosHelCM5->Write();
//    Phi_Scattered_685MeV_PosHelCM5->Write();
//
//    Phi_Scattered_265MeV_PosHelCM6->Write();
//    Phi_Scattered_335MeV_PosHelCM6->Write();
//    Phi_Scattered_405MeV_PosHelCM6->Write();
//    Phi_Scattered_475MeV_PosHelCM6->Write();
//    Phi_Scattered_545MeV_PosHelCM6->Write();
//    Phi_Scattered_615MeV_PosHelCM6->Write();
//    Phi_Scattered_685MeV_PosHelCM6->Write();
//
//    Phi_Scattered_265MeV_PosHelCM7->Write();
//    Phi_Scattered_335MeV_PosHelCM7->Write();
//    Phi_Scattered_405MeV_PosHelCM7->Write();
//    Phi_Scattered_475MeV_PosHelCM7->Write();
//    Phi_Scattered_545MeV_PosHelCM7->Write();
//    Phi_Scattered_615MeV_PosHelCM7->Write();
//    Phi_Scattered_685MeV_PosHelCM7->Write();
//
//    Phi_Scattered_265MeV_PosHelCM8->Write();
//    Phi_Scattered_335MeV_PosHelCM8->Write();
//    Phi_Scattered_405MeV_PosHelCM8->Write();
//    Phi_Scattered_475MeV_PosHelCM8->Write();
//    Phi_Scattered_545MeV_PosHelCM8->Write();
//    Phi_Scattered_615MeV_PosHelCM8->Write();
//    Phi_Scattered_685MeV_PosHelCM8->Write();

    f2.Write();

}
