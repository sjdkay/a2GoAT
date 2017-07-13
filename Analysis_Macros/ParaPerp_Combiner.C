#include "./includes_ParaPerpCombiner.h"

void ParaPerp_Combiner(){

  TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/Para/Physics_Total_Para_19_12_7_17.root"); // Open latest Para file

    TH1D* time_Para = (TH1D*)f->Get("time")->Clone();
    time_Para->SetName("time_Para");
    TH1D* time_cut_Para = (TH1D*)f->Get("time_cut")->Clone();
    time_cut_Para->SetName("time_cut_Para");

    TH1D* Eg_Para = (TH1D*)f->Get("Eg")->Clone();
    Eg_Para->SetName("Eg_Para");
    TH1D* PhiDifference_Para = (TH1D*)f->Get("PhiDifference")->Clone();
    PhiDifference_Para->SetName("PhiDifference_Para");
    TH1D* EpKin_Para = (TH1D*)f->Get("EpKin")->Clone();
    EpKin_Para->SetName("EpKin_Para");
    TH1D* EpCorrected_Para = (TH1D*)f->Get("EpCorrected")->Clone();
    EpCorrected_Para->SetName("EpCorrected_Para");
    TH1D* EpKinEpCorrDiff_Para = (TH1D*)f->Get("EpKinEpCorrDiff")->Clone();
    EpKinEpCorrDiff_Para->SetName("EpKinEpCorrDiff_Para");
    TH1D* EpEpCorrDiff_Para = (TH1D*)f->Get("EpEpCorrDiff")->Clone();
    EpEpCorrDiff_Para->SetName("EpEpCorrDiff_Para");

    TH1D* MMpEpCorrected_Para = (TH1D*)f->Get("MMpEpCorrected")->Clone();
    MMpEpCorrected_Para->SetName("MMpEpCorrected_Para");
    TH1D* OAngle_Para = (TH1D*)f->Get("OAngle")->Clone();
    OAngle_Para->SetName("OAngle_Para");
    TH1D* WCZnRecon_Para = (TH1D*)f->Get("WCZnRecon")->Clone();
    WCZnRecon_Para->SetName("WCZnRecon_Para");

    TH1D* MMp200300_Para = (TH1D*)f->Get("MMp200300")->Clone();
    MMp200300_Para->SetName("MMp200300_Para");
    TH1D* MMp300400_Para = (TH1D*)f->Get("MMp300400")->Clone();
    MMp300400_Para->SetName("MMp300400_Para");
    TH1D* MMp400500_Para = (TH1D*)f->Get("MMp400500")->Clone();
    MMp400500_Para->SetName("MMp400500_Para");
    TH1D* MMp500600_Para = (TH1D*)f->Get("MMp500600")->Clone();
    MMp500600_Para->SetName("MMp500600_Para");
    TH1D* MMp600700_Para = (TH1D*)f->Get("MMp600700")->Clone();
    MMp600700_Para->SetName("MMp600700_Para");
    TH1D* MMp700800_Para = (TH1D*)f->Get("MMp700800")->Clone();
    MMp700800_Para->SetName("MMp700800_Para");
    TH1D* MMp800900_Para = (TH1D*)f->Get("MMp800900")->Clone();
    MMp800900_Para->SetName("MMp800900_Para");

    TH1D* ZpDist_Para = (TH1D*)f->Get("ZpDist")->Clone();
    ZpDist_Para->SetName("ZpDist_Para");

    TH1D* Phip_430MeVCM1_Para = (TH1D*)f->Get("Phip_430MeVCM1")->Clone();
    Phip_430MeVCM1_Para->SetName("Phip_430MeVCM1_Para");
    TH1D* Phip_450MeVCM1_Para = (TH1D*)f->Get("Phip_450MeVCM1")->Clone();
    Phip_450MeVCM1_Para->SetName("Phip_450MeVCM1_Para");
    TH1D* Phip_470MeVCM1_Para = (TH1D*)f->Get("Phip_470MeVCM1")->Clone();
    Phip_470MeVCM1_Para->SetName("Phip_470MeVCM1_Para");
    TH1D* Phip_490MeVCM1_Para = (TH1D*)f->Get("Phip_490MeVCM1")->Clone();
    Phip_490MeVCM1_Para->SetName("Phip_490MeVCM1_Para");
    TH1D* Phip_510MeVCM1_Para = (TH1D*)f->Get("Phip_510MeVCM1")->Clone();
    Phip_510MeVCM1_Para->SetName("Phip_510MeVCM1_Para");
    TH1D* Phip_530MeVCM1_Para = (TH1D*)f->Get("Phip_530MeVCM1")->Clone();
    Phip_530MeVCM1_Para->SetName("Phip_530MeVCM1_Para");
    TH1D* Phip_550MeVCM1_Para = (TH1D*)f->Get("Phip_550MeVCM1")->Clone();
    Phip_550MeVCM1_Para->SetName("Phip_550MeVCM1_Para");
    TH1D* Phip_570MeVCM1_Para = (TH1D*)f->Get("Phip_570MeVCM1")->Clone();
    Phip_570MeVCM1_Para->SetName("Phip_570MeVCM1_Para");
    TH1D* Phip_590MeVCM1_Para = (TH1D*)f->Get("Phip_590MeVCM1")->Clone();
    Phip_590MeVCM1_Para->SetName("Phip_590MeVCM1_Para");
    TH1D* Phip_610MeVCM1_Para = (TH1D*)f->Get("Phip_610MeVCM1")->Clone();
    Phip_610MeVCM1_Para->SetName("Phip_610MeVCM1_Para");

    TH1D* Phip_430MeVCM2_Para = (TH1D*)f->Get("Phip_430MeVCM2")->Clone();
    Phip_430MeVCM2_Para->SetName("Phip_430MeVCM2_Para");
    TH1D* Phip_450MeVCM2_Para = (TH1D*)f->Get("Phip_450MeVCM2")->Clone();
    Phip_450MeVCM2_Para->SetName("Phip_450MeVCM2_Para");
    TH1D* Phip_470MeVCM2_Para = (TH1D*)f->Get("Phip_470MeVCM2")->Clone();
    Phip_470MeVCM2_Para->SetName("Phip_470MeVCM2_Para");
    TH1D* Phip_490MeVCM2_Para = (TH1D*)f->Get("Phip_490MeVCM2")->Clone();
    Phip_490MeVCM2_Para->SetName("Phip_490MeVCM2_Para");
    TH1D* Phip_510MeVCM2_Para = (TH1D*)f->Get("Phip_510MeVCM2")->Clone();
    Phip_510MeVCM2_Para->SetName("Phip_510MeVCM2_Para");
    TH1D* Phip_530MeVCM2_Para = (TH1D*)f->Get("Phip_530MeVCM2")->Clone();
    Phip_530MeVCM2_Para->SetName("Phip_530MeVCM2_Para");
    TH1D* Phip_550MeVCM2_Para = (TH1D*)f->Get("Phip_550MeVCM2")->Clone();
    Phip_550MeVCM2_Para->SetName("Phip_550MeVCM2_Para");
    TH1D* Phip_570MeVCM2_Para = (TH1D*)f->Get("Phip_570MeVCM2")->Clone();
    Phip_570MeVCM2_Para->SetName("Phip_570MeVCM2_Para");
    TH1D* Phip_590MeVCM2_Para = (TH1D*)f->Get("Phip_590MeVCM2")->Clone();
    Phip_590MeVCM2_Para->SetName("Phip_590MeVCM2_Para");
    TH1D* Phip_610MeVCM2_Para = (TH1D*)f->Get("Phip_610MeVCM2")->Clone();
    Phip_610MeVCM2_Para->SetName("Phip_610MeVCM2_Para");

    TH1D* Phip_430MeVCM3_Para = (TH1D*)f->Get("Phip_430MeVCM3")->Clone();
    Phip_430MeVCM3_Para->SetName("Phip_430MeVCM3_Para");
    TH1D* Phip_450MeVCM3_Para = (TH1D*)f->Get("Phip_450MeVCM3")->Clone();
    Phip_450MeVCM3_Para->SetName("Phip_450MeVCM3_Para");
    TH1D* Phip_470MeVCM3_Para = (TH1D*)f->Get("Phip_470MeVCM3")->Clone();
    Phip_470MeVCM3_Para->SetName("Phip_470MeVCM3_Para");
    TH1D* Phip_490MeVCM3_Para = (TH1D*)f->Get("Phip_490MeVCM3")->Clone();
    Phip_490MeVCM3_Para->SetName("Phip_490MeVCM3_Para");
    TH1D* Phip_510MeVCM3_Para = (TH1D*)f->Get("Phip_510MeVCM3")->Clone();
    Phip_510MeVCM3_Para->SetName("Phip_510MeVCM3_Para");
    TH1D* Phip_530MeVCM3_Para = (TH1D*)f->Get("Phip_530MeVCM3")->Clone();
    Phip_530MeVCM3_Para->SetName("Phip_530MeVCM3_Para");
    TH1D* Phip_550MeVCM3_Para = (TH1D*)f->Get("Phip_550MeVCM3")->Clone();
    Phip_550MeVCM3_Para->SetName("Phip_550MeVCM3_Para");
    TH1D* Phip_570MeVCM3_Para = (TH1D*)f->Get("Phip_570MeVCM3")->Clone();
    Phip_570MeVCM3_Para->SetName("Phip_570MeVCM3_Para");
    TH1D* Phip_590MeVCM3_Para = (TH1D*)f->Get("Phip_590MeVCM3")->Clone();
    Phip_590MeVCM3_Para->SetName("Phip_590MeVCM3_Para");
    TH1D* Phip_610MeVCM3_Para = (TH1D*)f->Get("Phip_610MeVCM3")->Clone();
    Phip_610MeVCM3_Para->SetName("Phip_610MeVCM3_Para");

    TH1D* Phip_430MeVCM4_Para = (TH1D*)f->Get("Phip_430MeVCM4")->Clone();
    Phip_430MeVCM4_Para->SetName("Phip_430MeVCM4_Para");
    TH1D* Phip_450MeVCM4_Para = (TH1D*)f->Get("Phip_450MeVCM4")->Clone();
    Phip_450MeVCM4_Para->SetName("Phip_450MeVCM4_Para");
    TH1D* Phip_470MeVCM4_Para = (TH1D*)f->Get("Phip_470MeVCM4")->Clone();
    Phip_470MeVCM4_Para->SetName("Phip_470MeVCM4_Para");
    TH1D* Phip_490MeVCM4_Para = (TH1D*)f->Get("Phip_490MeVCM4")->Clone();
    Phip_490MeVCM4_Para->SetName("Phip_490MeVCM4_Para");
    TH1D* Phip_510MeVCM4_Para = (TH1D*)f->Get("Phip_510MeVCM4")->Clone();
    Phip_510MeVCM4_Para->SetName("Phip_510MeVCM4_Para");
    TH1D* Phip_530MeVCM4_Para = (TH1D*)f->Get("Phip_530MeVCM4")->Clone();
    Phip_530MeVCM4_Para->SetName("Phip_530MeVCM4_Para");
    TH1D* Phip_550MeVCM4_Para = (TH1D*)f->Get("Phip_550MeVCM4")->Clone();
    Phip_550MeVCM4_Para->SetName("Phip_550MeVCM4_Para");
    TH1D* Phip_570MeVCM4_Para = (TH1D*)f->Get("Phip_570MeVCM4")->Clone();
    Phip_570MeVCM4_Para->SetName("Phip_570MeVCM4_Para");
    TH1D* Phip_590MeVCM4_Para = (TH1D*)f->Get("Phip_590MeVCM4")->Clone();
    Phip_590MeVCM4_Para->SetName("Phip_590MeVCM4_Para");
    TH1D* Phip_610MeVCM4_Para = (TH1D*)f->Get("Phip_610MeVCM4")->Clone();
    Phip_610MeVCM4_Para->SetName("Phip_610MeVCM4_Para");

    TH1D* Phip_430MeVCM5_Para = (TH1D*)f->Get("Phip_430MeVCM5")->Clone();
    Phip_430MeVCM5_Para->SetName("Phip_430MeVCM5_Para");
    TH1D* Phip_450MeVCM5_Para = (TH1D*)f->Get("Phip_450MeVCM5")->Clone();
    Phip_450MeVCM5_Para->SetName("Phip_450MeVCM5_Para");
    TH1D* Phip_470MeVCM5_Para = (TH1D*)f->Get("Phip_470MeVCM5")->Clone();
    Phip_470MeVCM5_Para->SetName("Phip_470MeVCM5_Para");
    TH1D* Phip_490MeVCM5_Para = (TH1D*)f->Get("Phip_490MeVCM5")->Clone();
    Phip_490MeVCM5_Para->SetName("Phip_490MeVCM5_Para");
    TH1D* Phip_510MeVCM5_Para = (TH1D*)f->Get("Phip_510MeVCM5")->Clone();
    Phip_510MeVCM5_Para->SetName("Phip_510MeVCM5_Para");
    TH1D* Phip_530MeVCM5_Para = (TH1D*)f->Get("Phip_530MeVCM5")->Clone();
    Phip_530MeVCM5_Para->SetName("Phip_530MeVCM5_Para");
    TH1D* Phip_550MeVCM5_Para = (TH1D*)f->Get("Phip_550MeVCM5")->Clone();
    Phip_550MeVCM5_Para->SetName("Phip_550MeVCM5_Para");
    TH1D* Phip_570MeVCM5_Para = (TH1D*)f->Get("Phip_570MeVCM5")->Clone();
    Phip_570MeVCM5_Para->SetName("Phip_570MeVCM5_Para");
    TH1D* Phip_590MeVCM5_Para = (TH1D*)f->Get("Phip_590MeVCM5")->Clone();
    Phip_590MeVCM5_Para->SetName("Phip_590MeVCM5_Para");
    TH1D* Phip_610MeVCM5_Para = (TH1D*)f->Get("Phip_610MeVCM5")->Clone();
    Phip_610MeVCM5_Para->SetName("Phip_610MeVCM5_Para");

    TH1D* Phip_430MeVCM6_Para = (TH1D*)f->Get("Phip_430MeVCM6")->Clone();
    Phip_430MeVCM6_Para->SetName("Phip_430MeVCM6_Para");
    TH1D* Phip_450MeVCM6_Para = (TH1D*)f->Get("Phip_450MeVCM6")->Clone();
    Phip_450MeVCM6_Para->SetName("Phip_450MeVCM6_Para");
    TH1D* Phip_470MeVCM6_Para = (TH1D*)f->Get("Phip_470MeVCM6")->Clone();
    Phip_470MeVCM6_Para->SetName("Phip_470MeVCM6_Para");
    TH1D* Phip_490MeVCM6_Para = (TH1D*)f->Get("Phip_490MeVCM6")->Clone();
    Phip_490MeVCM6_Para->SetName("Phip_490MeVCM6_Para");
    TH1D* Phip_510MeVCM6_Para = (TH1D*)f->Get("Phip_510MeVCM6")->Clone();
    Phip_510MeVCM6_Para->SetName("Phip_510MeVCM6_Para");
    TH1D* Phip_530MeVCM6_Para = (TH1D*)f->Get("Phip_530MeVCM6")->Clone();
    Phip_530MeVCM6_Para->SetName("Phip_530MeVCM6_Para");
    TH1D* Phip_550MeVCM6_Para = (TH1D*)f->Get("Phip_550MeVCM6")->Clone();
    Phip_550MeVCM6_Para->SetName("Phip_550MeVCM6_Para");
    TH1D* Phip_570MeVCM6_Para = (TH1D*)f->Get("Phip_570MeVCM6")->Clone();
    Phip_570MeVCM6_Para->SetName("Phip_570MeVCM6_Para");
    TH1D* Phip_590MeVCM6_Para = (TH1D*)f->Get("Phip_590MeVCM6")->Clone();
    Phip_590MeVCM6_Para->SetName("Phip_590MeVCM6_Para");
    TH1D* Phip_610MeVCM6_Para = (TH1D*)f->Get("Phip_610MeVCM6")->Clone();
    Phip_610MeVCM6_Para->SetName("Phip_610MeVCM6_Para");

    TH1D* Phip_430MeVCM7_Para = (TH1D*)f->Get("Phip_430MeVCM7")->Clone();
    Phip_430MeVCM7_Para->SetName("Phip_430MeVCM7_Para");
    TH1D* Phip_450MeVCM7_Para = (TH1D*)f->Get("Phip_450MeVCM7")->Clone();
    Phip_450MeVCM7_Para->SetName("Phip_450MeVCM7_Para");
    TH1D* Phip_470MeVCM7_Para = (TH1D*)f->Get("Phip_470MeVCM7")->Clone();
    Phip_470MeVCM7_Para->SetName("Phip_470MeVCM7_Para");
    TH1D* Phip_490MeVCM7_Para = (TH1D*)f->Get("Phip_490MeVCM7")->Clone();
    Phip_490MeVCM7_Para->SetName("Phip_490MeVCM7_Para");
    TH1D* Phip_510MeVCM7_Para = (TH1D*)f->Get("Phip_510MeVCM7")->Clone();
    Phip_510MeVCM7_Para->SetName("Phip_510MeVCM7_Para");
    TH1D* Phip_530MeVCM7_Para = (TH1D*)f->Get("Phip_530MeVCM7")->Clone();
    Phip_530MeVCM7_Para->SetName("Phip_530MeVCM7_Para");
    TH1D* Phip_550MeVCM7_Para = (TH1D*)f->Get("Phip_550MeVCM7")->Clone();
    Phip_550MeVCM7_Para->SetName("Phip_550MeVCM7_Para");
    TH1D* Phip_570MeVCM7_Para = (TH1D*)f->Get("Phip_570MeVCM7")->Clone();
    Phip_570MeVCM7_Para->SetName("Phip_570MeVCM7_Para");
    TH1D* Phip_590MeVCM7_Para = (TH1D*)f->Get("Phip_590MeVCM7")->Clone();
    Phip_590MeVCM7_Para->SetName("Phip_590MeVCM7_Para");
    TH1D* Phip_610MeVCM7_Para = (TH1D*)f->Get("Phip_610MeVCM7")->Clone();
    Phip_610MeVCM7_Para->SetName("Phip_610MeVCM7_Para");

    TH1D* Phip_430MeVCM8_Para = (TH1D*)f->Get("Phip_430MeVCM8")->Clone();
    Phip_430MeVCM8_Para->SetName("Phip_430MeVCM8_Para");
    TH1D* Phip_450MeVCM8_Para = (TH1D*)f->Get("Phip_450MeVCM8")->Clone();
    Phip_450MeVCM8_Para->SetName("Phip_450MeVCM8_Para");
    TH1D* Phip_470MeVCM8_Para = (TH1D*)f->Get("Phip_470MeVCM8")->Clone();
    Phip_470MeVCM8_Para->SetName("Phip_470MeVCM8_Para");
    TH1D* Phip_490MeVCM8_Para = (TH1D*)f->Get("Phip_490MeVCM8")->Clone();
    Phip_490MeVCM8_Para->SetName("Phip_490MeVCM8_Para");
    TH1D* Phip_510MeVCM8_Para = (TH1D*)f->Get("Phip_510MeVCM8")->Clone();
    Phip_510MeVCM8_Para->SetName("Phip_510MeVCM8_Para");
    TH1D* Phip_530MeVCM8_Para = (TH1D*)f->Get("Phip_530MeVCM8")->Clone();
    Phip_530MeVCM8_Para->SetName("Phip_530MeVCM8_Para");
    TH1D* Phip_550MeVCM8_Para = (TH1D*)f->Get("Phip_550MeVCM8")->Clone();
    Phip_550MeVCM8_Para->SetName("Phip_550MeVCM8_Para");
    TH1D* Phip_570MeVCM8_Para = (TH1D*)f->Get("Phip_570MeVCM8")->Clone();
    Phip_570MeVCM8_Para->SetName("Phip_570MeVCM8_Para");
    TH1D* Phip_590MeVCM8_Para = (TH1D*)f->Get("Phip_590MeVCM8")->Clone();
    Phip_590MeVCM8_Para->SetName("Phip_590MeVCM8_Para");
    TH1D* Phip_610MeVCM8_Para = (TH1D*)f->Get("Phip_610MeVCM8")->Clone();
    Phip_610MeVCM8_Para->SetName("Phip_610MeVCM8_Para");

    TH1D* Phip_430MeVCM9_Para = (TH1D*)f->Get("Phip_430MeVCM9")->Clone();
    Phip_430MeVCM9_Para->SetName("Phip_430MeVCM9_Para");
    TH1D* Phip_450MeVCM9_Para = (TH1D*)f->Get("Phip_450MeVCM9")->Clone();
    Phip_450MeVCM9_Para->SetName("Phip_450MeVCM9_Para");
    TH1D* Phip_470MeVCM9_Para = (TH1D*)f->Get("Phip_470MeVCM9")->Clone();
    Phip_470MeVCM9_Para->SetName("Phip_470MeVCM9_Para");
    TH1D* Phip_490MeVCM9_Para = (TH1D*)f->Get("Phip_490MeVCM9")->Clone();
    Phip_490MeVCM9_Para->SetName("Phip_490MeVCM9_Para");
    TH1D* Phip_510MeVCM9_Para = (TH1D*)f->Get("Phip_510MeVCM9")->Clone();
    Phip_510MeVCM9_Para->SetName("Phip_510MeVCM9_Para");
    TH1D* Phip_530MeVCM9_Para = (TH1D*)f->Get("Phip_530MeVCM9")->Clone();
    Phip_530MeVCM9_Para->SetName("Phip_530MeVCM9_Para");
    TH1D* Phip_550MeVCM9_Para = (TH1D*)f->Get("Phip_550MeVCM9")->Clone();
    Phip_550MeVCM9_Para->SetName("Phip_550MeVCM9_Para");
    TH1D* Phip_570MeVCM9_Para = (TH1D*)f->Get("Phip_570MeVCM9")->Clone();
    Phip_570MeVCM9_Para->SetName("Phip_570MeVCM9_Para");
    TH1D* Phip_590MeVCM9_Para = (TH1D*)f->Get("Phip_590MeVCM9")->Clone();
    Phip_590MeVCM9_Para->SetName("Phip_590MeVCM9_Para");
    TH1D* Phip_610MeVCM9_Para = (TH1D*)f->Get("Phip_610MeVCM9")->Clone();
    Phip_610MeVCM9_Para->SetName("Phip_610MeVCM9_Para");

    TH1D* Phip_430MeVCM10_Para = (TH1D*)f->Get("Phip_430MeVCM10")->Clone();
    Phip_430MeVCM10_Para->SetName("Phip_430MeVCM10_Para");
    TH1D* Phip_450MeVCM10_Para = (TH1D*)f->Get("Phip_450MeVCM10")->Clone();
    Phip_450MeVCM10_Para->SetName("Phip_450MeVCM10_Para");
    TH1D* Phip_470MeVCM10_Para = (TH1D*)f->Get("Phip_470MeVCM10")->Clone();
    Phip_470MeVCM10_Para->SetName("Phip_470MeVCM10_Para");
    TH1D* Phip_490MeVCM10_Para = (TH1D*)f->Get("Phip_490MeVCM10")->Clone();
    Phip_490MeVCM10_Para->SetName("Phip_490MeVCM10_Para");
    TH1D* Phip_510MeVCM10_Para = (TH1D*)f->Get("Phip_510MeVCM10")->Clone();
    Phip_510MeVCM10_Para->SetName("Phip_510MeVCM10_Para");
    TH1D* Phip_530MeVCM10_Para = (TH1D*)f->Get("Phip_530MeVCM10")->Clone();
    Phip_530MeVCM10_Para->SetName("Phip_530MeVCM10_Para");
    TH1D* Phip_550MeVCM10_Para = (TH1D*)f->Get("Phip_550MeVCM10")->Clone();
    Phip_550MeVCM10_Para->SetName("Phip_550MeVCM10_Para");
    TH1D* Phip_570MeVCM10_Para = (TH1D*)f->Get("Phip_570MeVCM10")->Clone();
    Phip_570MeVCM10_Para->SetName("Phip_570MeVCM10_Para");
    TH1D* Phip_590MeVCM10_Para = (TH1D*)f->Get("Phip_590MeVCM10")->Clone();
    Phip_590MeVCM10_Para->SetName("Phip_590MeVCM10_Para");
    TH1D* Phip_610MeVCM10_Para = (TH1D*)f->Get("Phip_610MeVCM10")->Clone();
    Phip_610MeVCM10_Para->SetName("Phip_610MeVCM10_Para");

    TH1D* ThetanDist_Para = (TH1D*)f->Get("ThetanDist")->Clone();
    ThetanDist_Para->SetName("ThetanDist_Para");
    TH1D* ThetanRecDist_Para = (TH1D*)f->Get("ThetanRecDist")->Clone();
    ThetanRecDist_Para->SetName("ThetanRecDist_Para");
    TH1D* ThetanDiffDist_Para = (TH1D*)f->Get("ThetanDiffDist")->Clone();
    ThetanDiffDist_Para->SetName("ThetanDiffDist_Para");
    TH2D* ThetanDiffZp_Para = (TH2D*)f->Get("ThetanDiffZp")->Clone();
    ThetanDiffZp_Para->SetName("ThetanDiffZp_Para");

    TH1D* ThetanCorrDist_Para = (TH1D*)f->Get("ThetanCorrDist")->Clone();
    ThetanCorrDist_Para->SetName("ThetanCorrDist_Para");
    TH1D* ThetanCorrDiffDist_Para = (TH1D*)f->Get("ThetanCorrDiffDist")->Clone();
    ThetanCorrDiffDist_Para->SetName("ThetanCorrDiffDist_Para");
    TH1D* ThetanCorrRecDiffDist_Para = (TH1D*)f->Get("ThetanCorrRecDiffDist")->Clone();
    ThetanCorrRecDiffDist_Para->SetName("ThetanCorrRecDiffDist_Para");
    TH2D* ThetanCorrDiffZp_Para = (TH2D*)f->Get("ThetanCorrDiffZp")->Clone();
    ThetanCorrDiffZp_Para->SetName("ThetanCorrDiffZp_Para");

    TH2D* E_dE_Para = (TH2D*)f->Get("E_dE")->Clone();
    E_dE_Para->SetName("E_dE_Para");
    TH2D* KinEp_dE_Para = (TH2D*)f->Get("KinEp_dE")->Clone();
    KinEp_dE_Para->SetName("KinEp_dE_Para");
    TH2D* E_KinEp_Para = (TH2D*)f->Get("E_KinEp")->Clone();
    E_KinEp_Para->SetName("E_KinEp_Para");

    TH1D* ThetaRecPiDiff_Para = (TH1D*)f->Get("ThetaRecPiDiff")->Clone();
    ThetaRecPiDiff_Para->SetName("ThetaRecPiDiff_Para");
    TH2D* ThetanThetaRecPi_Para = (TH2D*)f->Get("ThetanThetaRecPi")->Clone();
    ThetanThetaRecPi_Para->SetName("ThetanThetaRecPi_Para");
    TH2D* ThetanThetaRecPiDiff_Para = (TH2D*)f->Get("ThetanThetaRecPiDiff")->Clone();
    ThetanThetaRecPiDiff_Para->SetName("ThetanThetaRecPiDiff_Para");
    TH1D* ThetaRecPDiff_Para = (TH1D*)f->Get("ThetaRecPDiff")->Clone();
    ThetaRecPDiff_Para->SetName("ThetaRecPDiff_Para");
    TH2D* ThetanThetaRecP_Para = (TH2D*)f->Get("ThetanThetaRecP")->Clone();
    ThetanThetaRecP_Para->SetName("ThetanThetaRecP_Para");
    TH2D* ThetanThetaRecPDiff_Para = (TH2D*)f->Get("ThetanThetaRecPDiff")->Clone();
    ThetanThetaRecPDiff_Para->SetName("ThetanThetaRecPDiff_Para");

    TH2D* DeutKinPiKin_Para = (TH2D*)f->Get("DeutKinPiKin")->Clone();
    DeutKinPiKin_Para->SetName("DeutKinPiKin_Para");

    ////////////////////////////////////////////////////
    ////////////////////////////////////////////////////
    ///////////////// PARA DONE ////////////////////////
    ////////////////////////////////////////////////////
    ////////////////////////////////////////////////////

    TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/Perp/Physics_Total_Perp_19_12_7_17.root"); // Open latest Perp file

    TH1D* time_Perp = (TH1D*)f->Get("time")->Clone();
    time_Perp->SetName("time_Perp");
    TH1D* time_cut_Perp = (TH1D*)f->Get("time_cut")->Clone();
    time_cut_Perp->SetName("time_cut_Perp");

    TH1D* Eg_Perp = (TH1D*)f->Get("Eg")->Clone();
    Eg_Perp->SetName("Eg_Perp");
    TH1D* PhiDifference_Perp = (TH1D*)f->Get("PhiDifference")->Clone();
    PhiDifference_Perp->SetName("PhiDifference_Perp");
    TH1D* EpKin_Perp = (TH1D*)f->Get("EpKin")->Clone();
    EpKin_Perp->SetName("EpKin_Perp");
    TH1D* EpCorrected_Perp = (TH1D*)f->Get("EpCorrected")->Clone();
    EpCorrected_Perp->SetName("EpCorrected_Perp");
    TH1D* EpKinEpCorrDiff_Perp = (TH1D*)f->Get("EpKinEpCorrDiff")->Clone();
    EpKinEpCorrDiff_Perp->SetName("EpKinEpCorrDiff_Perp");
    TH1D* EpEpCorrDiff_Perp = (TH1D*)f->Get("EpEpCorrDiff")->Clone();
    EpEpCorrDiff_Perp->SetName("EpEpCorrDiff_Perp");

    TH1D* MMpEpCorrected_Perp = (TH1D*)f->Get("MMpEpCorrected")->Clone();
    MMpEpCorrected_Perp->SetName("MMpEpCorrected_Perp");
    TH1D* OAngle_Perp = (TH1D*)f->Get("OAngle")->Clone();
    OAngle_Perp->SetName("OAngle_Perp");
    TH1D* WCZnRecon_Perp = (TH1D*)f->Get("WCZnRecon")->Clone();
    WCZnRecon_Perp->SetName("WCZnRecon_Perp");

    TH1D* MMp200300_Perp = (TH1D*)f->Get("MMp200300")->Clone();
    MMp200300_Perp->SetName("MMp200300_Perp");
    TH1D* MMp300400_Perp = (TH1D*)f->Get("MMp300400")->Clone();
    MMp300400_Perp->SetName("MMp300400_Perp");
    TH1D* MMp400500_Perp = (TH1D*)f->Get("MMp400500")->Clone();
    MMp400500_Perp->SetName("MMp400500_Perp");
    TH1D* MMp500600_Perp = (TH1D*)f->Get("MMp500600")->Clone();
    MMp500600_Perp->SetName("MMp500600_Perp");
    TH1D* MMp600700_Perp = (TH1D*)f->Get("MMp600700")->Clone();
    MMp600700_Perp->SetName("MMp600700_Perp");
    TH1D* MMp700800_Perp = (TH1D*)f->Get("MMp700800")->Clone();
    MMp700800_Perp->SetName("MMp700800_Perp");
    TH1D* MMp800900_Perp = (TH1D*)f->Get("MMp800900")->Clone();
    MMp800900_Perp->SetName("MMp800900_Perp");

    TH1D* ZpDist_Perp = (TH1D*)f->Get("ZpDist")->Clone();
    ZpDist_Perp->SetName("ZpDist_Perp");

    TH1D* Phip_430MeVCM1_Perp = (TH1D*)f->Get("Phip_430MeVCM1")->Clone();
    Phip_430MeVCM1_Perp->SetName("Phip_430MeVCM1_Perp");
    TH1D* Phip_450MeVCM1_Perp = (TH1D*)f->Get("Phip_450MeVCM1")->Clone();
    Phip_450MeVCM1_Perp->SetName("Phip_450MeVCM1_Perp");
    TH1D* Phip_470MeVCM1_Perp = (TH1D*)f->Get("Phip_470MeVCM1")->Clone();
    Phip_470MeVCM1_Perp->SetName("Phip_470MeVCM1_Perp");
    TH1D* Phip_490MeVCM1_Perp = (TH1D*)f->Get("Phip_490MeVCM1")->Clone();
    Phip_490MeVCM1_Perp->SetName("Phip_490MeVCM1_Perp");
    TH1D* Phip_510MeVCM1_Perp = (TH1D*)f->Get("Phip_510MeVCM1")->Clone();
    Phip_510MeVCM1_Perp->SetName("Phip_510MeVCM1_Perp");
    TH1D* Phip_530MeVCM1_Perp = (TH1D*)f->Get("Phip_530MeVCM1")->Clone();
    Phip_530MeVCM1_Perp->SetName("Phip_530MeVCM1_Perp");
    TH1D* Phip_550MeVCM1_Perp = (TH1D*)f->Get("Phip_550MeVCM1")->Clone();
    Phip_550MeVCM1_Perp->SetName("Phip_550MeVCM1_Perp");
    TH1D* Phip_570MeVCM1_Perp = (TH1D*)f->Get("Phip_570MeVCM1")->Clone();
    Phip_570MeVCM1_Perp->SetName("Phip_570MeVCM1_Perp");
    TH1D* Phip_590MeVCM1_Perp = (TH1D*)f->Get("Phip_590MeVCM1")->Clone();
    Phip_590MeVCM1_Perp->SetName("Phip_590MeVCM1_Perp");
    TH1D* Phip_610MeVCM1_Perp = (TH1D*)f->Get("Phip_610MeVCM1")->Clone();
    Phip_610MeVCM1_Perp->SetName("Phip_610MeVCM1_Perp");

    TH1D* Phip_430MeVCM2_Perp = (TH1D*)f->Get("Phip_430MeVCM2")->Clone();
    Phip_430MeVCM2_Perp->SetName("Phip_430MeVCM2_Perp");
    TH1D* Phip_450MeVCM2_Perp = (TH1D*)f->Get("Phip_450MeVCM2")->Clone();
    Phip_450MeVCM2_Perp->SetName("Phip_450MeVCM2_Perp");
    TH1D* Phip_470MeVCM2_Perp = (TH1D*)f->Get("Phip_470MeVCM2")->Clone();
    Phip_470MeVCM2_Perp->SetName("Phip_470MeVCM2_Perp");
    TH1D* Phip_490MeVCM2_Perp = (TH1D*)f->Get("Phip_490MeVCM2")->Clone();
    Phip_490MeVCM2_Perp->SetName("Phip_490MeVCM2_Perp");
    TH1D* Phip_510MeVCM2_Perp = (TH1D*)f->Get("Phip_510MeVCM2")->Clone();
    Phip_510MeVCM2_Perp->SetName("Phip_510MeVCM2_Perp");
    TH1D* Phip_530MeVCM2_Perp = (TH1D*)f->Get("Phip_530MeVCM2")->Clone();
    Phip_530MeVCM2_Perp->SetName("Phip_530MeVCM2_Perp");
    TH1D* Phip_550MeVCM2_Perp = (TH1D*)f->Get("Phip_550MeVCM2")->Clone();
    Phip_550MeVCM2_Perp->SetName("Phip_550MeVCM2_Perp");
    TH1D* Phip_570MeVCM2_Perp = (TH1D*)f->Get("Phip_570MeVCM2")->Clone();
    Phip_570MeVCM2_Perp->SetName("Phip_570MeVCM2_Perp");
    TH1D* Phip_590MeVCM2_Perp = (TH1D*)f->Get("Phip_590MeVCM2")->Clone();
    Phip_590MeVCM2_Perp->SetName("Phip_590MeVCM2_Perp");
    TH1D* Phip_610MeVCM2_Perp = (TH1D*)f->Get("Phip_610MeVCM2")->Clone();
    Phip_610MeVCM2_Perp->SetName("Phip_610MeVCM2_Perp");

    TH1D* Phip_430MeVCM3_Perp = (TH1D*)f->Get("Phip_430MeVCM3")->Clone();
    Phip_430MeVCM3_Perp->SetName("Phip_430MeVCM3_Perp");
    TH1D* Phip_450MeVCM3_Perp = (TH1D*)f->Get("Phip_450MeVCM3")->Clone();
    Phip_450MeVCM3_Perp->SetName("Phip_450MeVCM3_Perp");
    TH1D* Phip_470MeVCM3_Perp = (TH1D*)f->Get("Phip_470MeVCM3")->Clone();
    Phip_470MeVCM3_Perp->SetName("Phip_470MeVCM3_Perp");
    TH1D* Phip_490MeVCM3_Perp = (TH1D*)f->Get("Phip_490MeVCM3")->Clone();
    Phip_490MeVCM3_Perp->SetName("Phip_490MeVCM3_Perp");
    TH1D* Phip_510MeVCM3_Perp = (TH1D*)f->Get("Phip_510MeVCM3")->Clone();
    Phip_510MeVCM3_Perp->SetName("Phip_510MeVCM3_Perp");
    TH1D* Phip_530MeVCM3_Perp = (TH1D*)f->Get("Phip_530MeVCM3")->Clone();
    Phip_530MeVCM3_Perp->SetName("Phip_530MeVCM3_Perp");
    TH1D* Phip_550MeVCM3_Perp = (TH1D*)f->Get("Phip_550MeVCM3")->Clone();
    Phip_550MeVCM3_Perp->SetName("Phip_550MeVCM3_Perp");
    TH1D* Phip_570MeVCM3_Perp = (TH1D*)f->Get("Phip_570MeVCM3")->Clone();
    Phip_570MeVCM3_Perp->SetName("Phip_570MeVCM3_Perp");
    TH1D* Phip_590MeVCM3_Perp = (TH1D*)f->Get("Phip_590MeVCM3")->Clone();
    Phip_590MeVCM3_Perp->SetName("Phip_590MeVCM3_Perp");
    TH1D* Phip_610MeVCM3_Perp = (TH1D*)f->Get("Phip_610MeVCM3")->Clone();
    Phip_610MeVCM3_Perp->SetName("Phip_610MeVCM3_Perp");

    TH1D* Phip_430MeVCM4_Perp = (TH1D*)f->Get("Phip_430MeVCM4")->Clone();
    Phip_430MeVCM4_Perp->SetName("Phip_430MeVCM4_Perp");
    TH1D* Phip_450MeVCM4_Perp = (TH1D*)f->Get("Phip_450MeVCM4")->Clone();
    Phip_450MeVCM4_Perp->SetName("Phip_450MeVCM4_Perp");
    TH1D* Phip_470MeVCM4_Perp = (TH1D*)f->Get("Phip_470MeVCM4")->Clone();
    Phip_470MeVCM4_Perp->SetName("Phip_470MeVCM4_Perp");
    TH1D* Phip_490MeVCM4_Perp = (TH1D*)f->Get("Phip_490MeVCM4")->Clone();
    Phip_490MeVCM4_Perp->SetName("Phip_490MeVCM4_Perp");
    TH1D* Phip_510MeVCM4_Perp = (TH1D*)f->Get("Phip_510MeVCM4")->Clone();
    Phip_510MeVCM4_Perp->SetName("Phip_510MeVCM4_Perp");
    TH1D* Phip_530MeVCM4_Perp = (TH1D*)f->Get("Phip_530MeVCM4")->Clone();
    Phip_530MeVCM4_Perp->SetName("Phip_530MeVCM4_Perp");
    TH1D* Phip_550MeVCM4_Perp = (TH1D*)f->Get("Phip_550MeVCM4")->Clone();
    Phip_550MeVCM4_Perp->SetName("Phip_550MeVCM4_Perp");
    TH1D* Phip_570MeVCM4_Perp = (TH1D*)f->Get("Phip_570MeVCM4")->Clone();
    Phip_570MeVCM4_Perp->SetName("Phip_570MeVCM4_Perp");
    TH1D* Phip_590MeVCM4_Perp = (TH1D*)f->Get("Phip_590MeVCM4")->Clone();
    Phip_590MeVCM4_Perp->SetName("Phip_590MeVCM4_Perp");
    TH1D* Phip_610MeVCM4_Perp = (TH1D*)f->Get("Phip_610MeVCM4")->Clone();
    Phip_610MeVCM4_Perp->SetName("Phip_610MeVCM4_Perp");

    TH1D* Phip_430MeVCM5_Perp = (TH1D*)f->Get("Phip_430MeVCM5")->Clone();
    Phip_430MeVCM5_Perp->SetName("Phip_430MeVCM5_Perp");
    TH1D* Phip_450MeVCM5_Perp = (TH1D*)f->Get("Phip_450MeVCM5")->Clone();
    Phip_450MeVCM5_Perp->SetName("Phip_450MeVCM5_Perp");
    TH1D* Phip_470MeVCM5_Perp = (TH1D*)f->Get("Phip_470MeVCM5")->Clone();
    Phip_470MeVCM5_Perp->SetName("Phip_470MeVCM5_Perp");
    TH1D* Phip_490MeVCM5_Perp = (TH1D*)f->Get("Phip_490MeVCM5")->Clone();
    Phip_490MeVCM5_Perp->SetName("Phip_490MeVCM5_Perp");
    TH1D* Phip_510MeVCM5_Perp = (TH1D*)f->Get("Phip_510MeVCM5")->Clone();
    Phip_510MeVCM5_Perp->SetName("Phip_510MeVCM5_Perp");
    TH1D* Phip_530MeVCM5_Perp = (TH1D*)f->Get("Phip_530MeVCM5")->Clone();
    Phip_530MeVCM5_Perp->SetName("Phip_530MeVCM5_Perp");
    TH1D* Phip_550MeVCM5_Perp = (TH1D*)f->Get("Phip_550MeVCM5")->Clone();
    Phip_550MeVCM5_Perp->SetName("Phip_550MeVCM5_Perp");
    TH1D* Phip_570MeVCM5_Perp = (TH1D*)f->Get("Phip_570MeVCM5")->Clone();
    Phip_570MeVCM5_Perp->SetName("Phip_570MeVCM5_Perp");
    TH1D* Phip_590MeVCM5_Perp = (TH1D*)f->Get("Phip_590MeVCM5")->Clone();
    Phip_590MeVCM5_Perp->SetName("Phip_590MeVCM5_Perp");
    TH1D* Phip_610MeVCM5_Perp = (TH1D*)f->Get("Phip_610MeVCM5")->Clone();
    Phip_610MeVCM5_Perp->SetName("Phip_610MeVCM5_Perp");

    TH1D* Phip_430MeVCM6_Perp = (TH1D*)f->Get("Phip_430MeVCM6")->Clone();
    Phip_430MeVCM6_Perp->SetName("Phip_430MeVCM6_Perp");
    TH1D* Phip_450MeVCM6_Perp = (TH1D*)f->Get("Phip_450MeVCM6")->Clone();
    Phip_450MeVCM6_Perp->SetName("Phip_450MeVCM6_Perp");
    TH1D* Phip_470MeVCM6_Perp = (TH1D*)f->Get("Phip_470MeVCM6")->Clone();
    Phip_470MeVCM6_Perp->SetName("Phip_470MeVCM6_Perp");
    TH1D* Phip_490MeVCM6_Perp = (TH1D*)f->Get("Phip_490MeVCM6")->Clone();
    Phip_490MeVCM6_Perp->SetName("Phip_490MeVCM6_Perp");
    TH1D* Phip_510MeVCM6_Perp = (TH1D*)f->Get("Phip_510MeVCM6")->Clone();
    Phip_510MeVCM6_Perp->SetName("Phip_510MeVCM6_Perp");
    TH1D* Phip_530MeVCM6_Perp = (TH1D*)f->Get("Phip_530MeVCM6")->Clone();
    Phip_530MeVCM6_Perp->SetName("Phip_530MeVCM6_Perp");
    TH1D* Phip_550MeVCM6_Perp = (TH1D*)f->Get("Phip_550MeVCM6")->Clone();
    Phip_550MeVCM6_Perp->SetName("Phip_550MeVCM6_Perp");
    TH1D* Phip_570MeVCM6_Perp = (TH1D*)f->Get("Phip_570MeVCM6")->Clone();
    Phip_570MeVCM6_Perp->SetName("Phip_570MeVCM6_Perp");
    TH1D* Phip_590MeVCM6_Perp = (TH1D*)f->Get("Phip_590MeVCM6")->Clone();
    Phip_590MeVCM6_Perp->SetName("Phip_590MeVCM6_Perp");
    TH1D* Phip_610MeVCM6_Perp = (TH1D*)f->Get("Phip_610MeVCM6")->Clone();
    Phip_610MeVCM6_Perp->SetName("Phip_610MeVCM6_Perp");

    TH1D* Phip_430MeVCM7_Perp = (TH1D*)f->Get("Phip_430MeVCM7")->Clone();
    Phip_430MeVCM7_Perp->SetName("Phip_430MeVCM7_Perp");
    TH1D* Phip_450MeVCM7_Perp = (TH1D*)f->Get("Phip_450MeVCM7")->Clone();
    Phip_450MeVCM7_Perp->SetName("Phip_450MeVCM7_Perp");
    TH1D* Phip_470MeVCM7_Perp = (TH1D*)f->Get("Phip_470MeVCM7")->Clone();
    Phip_470MeVCM7_Perp->SetName("Phip_470MeVCM7_Perp");
    TH1D* Phip_490MeVCM7_Perp = (TH1D*)f->Get("Phip_490MeVCM7")->Clone();
    Phip_490MeVCM7_Perp->SetName("Phip_490MeVCM7_Perp");
    TH1D* Phip_510MeVCM7_Perp = (TH1D*)f->Get("Phip_510MeVCM7")->Clone();
    Phip_510MeVCM7_Perp->SetName("Phip_510MeVCM7_Perp");
    TH1D* Phip_530MeVCM7_Perp = (TH1D*)f->Get("Phip_530MeVCM7")->Clone();
    Phip_530MeVCM7_Perp->SetName("Phip_530MeVCM7_Perp");
    TH1D* Phip_550MeVCM7_Perp = (TH1D*)f->Get("Phip_550MeVCM7")->Clone();
    Phip_550MeVCM7_Perp->SetName("Phip_550MeVCM7_Perp");
    TH1D* Phip_570MeVCM7_Perp = (TH1D*)f->Get("Phip_570MeVCM7")->Clone();
    Phip_570MeVCM7_Perp->SetName("Phip_570MeVCM7_Perp");
    TH1D* Phip_590MeVCM7_Perp = (TH1D*)f->Get("Phip_590MeVCM7")->Clone();
    Phip_590MeVCM7_Perp->SetName("Phip_590MeVCM7_Perp");
    TH1D* Phip_610MeVCM7_Perp = (TH1D*)f->Get("Phip_610MeVCM7")->Clone();
    Phip_610MeVCM7_Perp->SetName("Phip_610MeVCM7_Perp");

    TH1D* Phip_430MeVCM8_Perp = (TH1D*)f->Get("Phip_430MeVCM8")->Clone();
    Phip_430MeVCM8_Perp->SetName("Phip_430MeVCM8_Perp");
    TH1D* Phip_450MeVCM8_Perp = (TH1D*)f->Get("Phip_450MeVCM8")->Clone();
    Phip_450MeVCM8_Perp->SetName("Phip_450MeVCM8_Perp");
    TH1D* Phip_470MeVCM8_Perp = (TH1D*)f->Get("Phip_470MeVCM8")->Clone();
    Phip_470MeVCM8_Perp->SetName("Phip_470MeVCM8_Perp");
    TH1D* Phip_490MeVCM8_Perp = (TH1D*)f->Get("Phip_490MeVCM8")->Clone();
    Phip_490MeVCM8_Perp->SetName("Phip_490MeVCM8_Perp");
    TH1D* Phip_510MeVCM8_Perp = (TH1D*)f->Get("Phip_510MeVCM8")->Clone();
    Phip_510MeVCM8_Perp->SetName("Phip_510MeVCM8_Perp");
    TH1D* Phip_530MeVCM8_Perp = (TH1D*)f->Get("Phip_530MeVCM8")->Clone();
    Phip_530MeVCM8_Perp->SetName("Phip_530MeVCM8_Perp");
    TH1D* Phip_550MeVCM8_Perp = (TH1D*)f->Get("Phip_550MeVCM8")->Clone();
    Phip_550MeVCM8_Perp->SetName("Phip_550MeVCM8_Perp");
    TH1D* Phip_570MeVCM8_Perp = (TH1D*)f->Get("Phip_570MeVCM8")->Clone();
    Phip_570MeVCM8_Perp->SetName("Phip_570MeVCM8_Perp");
    TH1D* Phip_590MeVCM8_Perp = (TH1D*)f->Get("Phip_590MeVCM8")->Clone();
    Phip_590MeVCM8_Perp->SetName("Phip_590MeVCM8_Perp");
    TH1D* Phip_610MeVCM8_Perp = (TH1D*)f->Get("Phip_610MeVCM8")->Clone();
    Phip_610MeVCM8_Perp->SetName("Phip_610MeVCM8_Perp");

    TH1D* Phip_430MeVCM9_Perp = (TH1D*)f->Get("Phip_430MeVCM9")->Clone();
    Phip_430MeVCM9_Perp->SetName("Phip_430MeVCM9_Perp");
    TH1D* Phip_450MeVCM9_Perp = (TH1D*)f->Get("Phip_450MeVCM9")->Clone();
    Phip_450MeVCM9_Perp->SetName("Phip_450MeVCM9_Perp");
    TH1D* Phip_470MeVCM9_Perp = (TH1D*)f->Get("Phip_470MeVCM9")->Clone();
    Phip_470MeVCM9_Perp->SetName("Phip_470MeVCM9_Perp");
    TH1D* Phip_490MeVCM9_Perp = (TH1D*)f->Get("Phip_490MeVCM9")->Clone();
    Phip_490MeVCM9_Perp->SetName("Phip_490MeVCM9_Perp");
    TH1D* Phip_510MeVCM9_Perp = (TH1D*)f->Get("Phip_510MeVCM9")->Clone();
    Phip_510MeVCM9_Perp->SetName("Phip_510MeVCM9_Perp");
    TH1D* Phip_530MeVCM9_Perp = (TH1D*)f->Get("Phip_530MeVCM9")->Clone();
    Phip_530MeVCM9_Perp->SetName("Phip_530MeVCM9_Perp");
    TH1D* Phip_550MeVCM9_Perp = (TH1D*)f->Get("Phip_550MeVCM9")->Clone();
    Phip_550MeVCM9_Perp->SetName("Phip_550MeVCM9_Perp");
    TH1D* Phip_570MeVCM9_Perp = (TH1D*)f->Get("Phip_570MeVCM9")->Clone();
    Phip_570MeVCM9_Perp->SetName("Phip_570MeVCM9_Perp");
    TH1D* Phip_590MeVCM9_Perp = (TH1D*)f->Get("Phip_590MeVCM9")->Clone();
    Phip_590MeVCM9_Perp->SetName("Phip_590MeVCM9_Perp");
    TH1D* Phip_610MeVCM9_Perp = (TH1D*)f->Get("Phip_610MeVCM9")->Clone();
    Phip_610MeVCM9_Perp->SetName("Phip_610MeVCM9_Perp");

    TH1D* Phip_430MeVCM10_Perp = (TH1D*)f->Get("Phip_430MeVCM10")->Clone();
    Phip_430MeVCM10_Perp->SetName("Phip_430MeVCM10_Perp");
    TH1D* Phip_450MeVCM10_Perp = (TH1D*)f->Get("Phip_450MeVCM10")->Clone();
    Phip_450MeVCM10_Perp->SetName("Phip_450MeVCM10_Perp");
    TH1D* Phip_470MeVCM10_Perp = (TH1D*)f->Get("Phip_470MeVCM10")->Clone();
    Phip_470MeVCM10_Perp->SetName("Phip_470MeVCM10_Perp");
    TH1D* Phip_490MeVCM10_Perp = (TH1D*)f->Get("Phip_490MeVCM10")->Clone();
    Phip_490MeVCM10_Perp->SetName("Phip_490MeVCM10_Perp");
    TH1D* Phip_510MeVCM10_Perp = (TH1D*)f->Get("Phip_510MeVCM10")->Clone();
    Phip_510MeVCM10_Perp->SetName("Phip_510MeVCM10_Perp");
    TH1D* Phip_530MeVCM10_Perp = (TH1D*)f->Get("Phip_530MeVCM10")->Clone();
    Phip_530MeVCM10_Perp->SetName("Phip_530MeVCM10_Perp");
    TH1D* Phip_550MeVCM10_Perp = (TH1D*)f->Get("Phip_550MeVCM10")->Clone();
    Phip_550MeVCM10_Perp->SetName("Phip_550MeVCM10_Perp");
    TH1D* Phip_570MeVCM10_Perp = (TH1D*)f->Get("Phip_570MeVCM10")->Clone();
    Phip_570MeVCM10_Perp->SetName("Phip_570MeVCM10_Perp");
    TH1D* Phip_590MeVCM10_Perp = (TH1D*)f->Get("Phip_590MeVCM10")->Clone();
    Phip_590MeVCM10_Perp->SetName("Phip_590MeVCM10_Perp");
    TH1D* Phip_610MeVCM10_Perp = (TH1D*)f->Get("Phip_610MeVCM10")->Clone();
    Phip_610MeVCM10_Perp->SetName("Phip_610MeVCM10_Perp");

    TH1D* ThetanDist_Perp = (TH1D*)f->Get("ThetanDist")->Clone();
    ThetanDist_Perp->SetName("ThetanDist_Perp");
    TH1D* ThetanRecDist_Perp = (TH1D*)f->Get("ThetanRecDist")->Clone();
    ThetanRecDist_Perp->SetName("ThetanRecDist_Perp");
    TH1D* ThetanDiffDist_Perp = (TH1D*)f->Get("ThetanDiffDist")->Clone();
    ThetanDiffDist_Perp->SetName("ThetanDiffDist_Perp");
    TH2D* ThetanDiffZp_Perp = (TH2D*)f->Get("ThetanDiffZp")->Clone();
    ThetanDiffZp_Perp->SetName("ThetanDiffZp_Perp");

    TH1D* ThetanCorrDist_Perp = (TH1D*)f->Get("ThetanCorrDist")->Clone();
    ThetanCorrDist_Perp->SetName("ThetanCorrDist_Perp");
    TH1D* ThetanCorrDiffDist_Perp = (TH1D*)f->Get("ThetanCorrDiffDist")->Clone();
    ThetanCorrDiffDist_Perp->SetName("ThetanCorrDiffDist_Perp");
    TH1D* ThetanCorrRecDiffDist_Perp = (TH1D*)f->Get("ThetanCorrRecDiffDist")->Clone();
    ThetanCorrRecDiffDist_Perp->SetName("ThetanCorrRecDiffDist_Perp");
    TH2D* ThetanCorrDiffZp_Perp = (TH2D*)f->Get("ThetanCorrDiffZp")->Clone();
    ThetanCorrDiffZp_Perp->SetName("ThetanCorrDiffZp_Perp");

    TH2D* E_dE_Perp = (TH2D*)f->Get("E_dE")->Clone();
    E_dE_Perp->SetName("E_dE_Perp");
    TH2D* KinEp_dE_Perp = (TH2D*)f->Get("KinEp_dE")->Clone();
    KinEp_dE_Perp->SetName("KinEp_dE_Perp");
    TH2D* E_KinEp_Perp = (TH2D*)f->Get("E_KinEp")->Clone();
    E_KinEp_Perp->SetName("E_KinEp_Perp");

    TH1D* ThetaRecPiDiff_Perp = (TH1D*)f->Get("ThetaRecPiDiff")->Clone();
    ThetaRecPiDiff_Perp->SetName("ThetaRecPiDiff_Perp");
    TH2D* ThetanThetaRecPi_Perp = (TH2D*)f->Get("ThetanThetaRecPi")->Clone();
    ThetanThetaRecPi_Perp->SetName("ThetanThetaRecPi_Perp");
    TH2D* ThetanThetaRecPiDiff_Perp = (TH2D*)f->Get("ThetanThetaRecPiDiff")->Clone();
    ThetanThetaRecPiDiff_Perp->SetName("ThetanThetaRecPiDiff_Perp");
    TH1D* ThetaRecPDiff_Perp = (TH1D*)f->Get("ThetaRecPDiff")->Clone();
    ThetaRecPDiff_Perp->SetName("ThetaRecPDiff_Perp");
    TH2D* ThetanThetaRecP_Perp = (TH2D*)f->Get("ThetanThetaRecP")->Clone();
    ThetanThetaRecP_Perp->SetName("ThetanThetaRecP_Perp");
    TH2D* ThetanThetaRecPDiff_Perp = (TH2D*)f->Get("ThetanThetaRecPDiff")->Clone();
    ThetanThetaRecPDiff_Perp->SetName("ThetanThetaRecPDiff_Perp");

    TH2D* DeutKinPiKin_Perp = (TH2D*)f->Get("DeutKinPiKin")->Clone();
    DeutKinPiKin_Perp->SetName("DeutKinPiKin_Perp");

    ////////////////////////////////////////////////////
    ////////////////////////////////////////////////////
    ///////////////// PERP DONE ////////////////////////
    ////////////////////////////////////////////////////
    ////////////////////////////////////////////////////

    TFile f2("ParaPerp_Total_19_Combined.root", "RECREATE");

    time_Para->Write();
    time_cut_Para->Write();

    Eg_Para->Write();
    PhiDifference_Para->Write();
    EpKin_Para->Write();
    EpCorrected_Para->Write();
    EpKinEpCorrDiff_Para->Write();
    EpEpCorrDiff_Para->Write();

    MMpEpCorrected_Para->Write();
    OAngle_Para->Write();
    WCZnRecon_Para->Write();

    ZpDist_Para->Write();

    MMp200300_Para->Write();
    MMp300400_Para->Write();
    MMp400500_Para->Write();
    MMp500600_Para->Write();
    MMp600700_Para->Write();
    MMp700800_Para->Write();
    MMp800900_Para->Write();

    Phip_430MeVCM1_Para->Write();
    Phip_450MeVCM1_Para->Write();
    Phip_470MeVCM1_Para->Write();
    Phip_490MeVCM1_Para->Write();
    Phip_510MeVCM1_Para->Write();
    Phip_530MeVCM1_Para->Write();
    Phip_550MeVCM1_Para->Write();
    Phip_570MeVCM1_Para->Write();
    Phip_590MeVCM1_Para->Write();
    Phip_610MeVCM1_Para->Write();

    Phip_430MeVCM2_Para->Write();
    Phip_450MeVCM2_Para->Write();
    Phip_470MeVCM2_Para->Write();
    Phip_490MeVCM2_Para->Write();
    Phip_510MeVCM2_Para->Write();
    Phip_530MeVCM2_Para->Write();
    Phip_550MeVCM2_Para->Write();
    Phip_570MeVCM2_Para->Write();
    Phip_590MeVCM2_Para->Write();
    Phip_610MeVCM2_Para->Write();

    Phip_430MeVCM3_Para->Write();
    Phip_450MeVCM3_Para->Write();
    Phip_470MeVCM3_Para->Write();
    Phip_490MeVCM3_Para->Write();
    Phip_510MeVCM3_Para->Write();
    Phip_530MeVCM3_Para->Write();
    Phip_550MeVCM3_Para->Write();
    Phip_570MeVCM3_Para->Write();
    Phip_590MeVCM3_Para->Write();
    Phip_610MeVCM3_Para->Write();

    Phip_430MeVCM4_Para->Write();
    Phip_450MeVCM4_Para->Write();
    Phip_470MeVCM4_Para->Write();
    Phip_490MeVCM4_Para->Write();
    Phip_510MeVCM4_Para->Write();
    Phip_530MeVCM4_Para->Write();
    Phip_550MeVCM4_Para->Write();
    Phip_570MeVCM4_Para->Write();
    Phip_590MeVCM4_Para->Write();
    Phip_610MeVCM4_Para->Write();

    Phip_430MeVCM5_Para->Write();
    Phip_450MeVCM5_Para->Write();
    Phip_470MeVCM5_Para->Write();
    Phip_490MeVCM5_Para->Write();
    Phip_510MeVCM5_Para->Write();
    Phip_530MeVCM5_Para->Write();
    Phip_550MeVCM5_Para->Write();
    Phip_570MeVCM5_Para->Write();
    Phip_590MeVCM5_Para->Write();
    Phip_610MeVCM5_Para->Write();

    Phip_430MeVCM6_Para->Write();
    Phip_450MeVCM6_Para->Write();
    Phip_470MeVCM6_Para->Write();
    Phip_490MeVCM6_Para->Write();
    Phip_510MeVCM6_Para->Write();
    Phip_530MeVCM6_Para->Write();
    Phip_550MeVCM6_Para->Write();
    Phip_570MeVCM6_Para->Write();
    Phip_590MeVCM6_Para->Write();
    Phip_610MeVCM6_Para->Write();

    Phip_430MeVCM7_Para->Write();
    Phip_450MeVCM7_Para->Write();
    Phip_470MeVCM7_Para->Write();
    Phip_490MeVCM7_Para->Write();
    Phip_510MeVCM7_Para->Write();
    Phip_530MeVCM7_Para->Write();
    Phip_550MeVCM7_Para->Write();
    Phip_570MeVCM7_Para->Write();
    Phip_590MeVCM7_Para->Write();
    Phip_610MeVCM7_Para->Write();

    Phip_430MeVCM8_Para->Write();
    Phip_450MeVCM8_Para->Write();
    Phip_470MeVCM8_Para->Write();
    Phip_490MeVCM8_Para->Write();
    Phip_510MeVCM8_Para->Write();
    Phip_530MeVCM8_Para->Write();
    Phip_550MeVCM8_Para->Write();
    Phip_570MeVCM8_Para->Write();
    Phip_590MeVCM8_Para->Write();
    Phip_610MeVCM8_Para->Write();

    Phip_430MeVCM9_Para->Write();
    Phip_450MeVCM9_Para->Write();
    Phip_470MeVCM9_Para->Write();
    Phip_490MeVCM9_Para->Write();
    Phip_510MeVCM9_Para->Write();
    Phip_530MeVCM9_Para->Write();
    Phip_550MeVCM9_Para->Write();
    Phip_570MeVCM9_Para->Write();
    Phip_590MeVCM9_Para->Write();
    Phip_610MeVCM9_Para->Write();

    Phip_430MeVCM10_Para->Write();
    Phip_450MeVCM10_Para->Write();
    Phip_470MeVCM10_Para->Write();
    Phip_490MeVCM10_Para->Write();
    Phip_510MeVCM10_Para->Write();
    Phip_530MeVCM10_Para->Write();
    Phip_550MeVCM10_Para->Write();
    Phip_570MeVCM10_Para->Write();
    Phip_590MeVCM10_Para->Write();
    Phip_610MeVCM10_Para->Write();

    ThetanDist_Para->Write();
    ThetanRecDist_Para->Write();
    ThetanDiffDist_Para->Write();
    ThetanDiffZp_Para->Write();
    ThetanCorrDist_Para->Write();
    ThetanCorrDiffDist_Para->Write();
    ThetanCorrRecDiffDist_Para->Write();
    ThetanCorrDiffZp_Para->Write();

    E_dE_Para->Write();;
    KinEp_dE_Para->Write();
    E_KinEp_Para->Write();

    ThetaRecPiDiff_Para->Write();
    ThetanThetaRecPi_Para->Write();
    ThetanThetaRecPiDiff_Para->Write();
    ThetaRecPDiff_Para->Write();
    ThetanThetaRecP_Para->Write();
    ThetanThetaRecPDiff_Para->Write();

    DeutKinPiKin_Para->Write();

    ////////////////////////////////////////////////////
    ////////////////////////////////////////////////////
    ///////////////// PARA DONE ////////////////////////
    ////////////////////////////////////////////////////
    ////////////////////////////////////////////////////

    time_Perp->Write();
    time_cut_Perp->Write();

    Eg_Perp->Write();
    PhiDifference_Perp->Write();
    EpKin_Perp->Write();
    EpCorrected_Perp->Write();
    EpKinEpCorrDiff_Perp->Write();
    EpEpCorrDiff_Perp->Write();

    MMpEpCorrected_Perp->Write();
    OAngle_Perp->Write();
    WCZnRecon_Perp->Write();

    ZpDist_Perp->Write();

    MMp200300_Perp->Write();
    MMp300400_Perp->Write();
    MMp400500_Perp->Write();
    MMp500600_Perp->Write();
    MMp600700_Perp->Write();
    MMp700800_Perp->Write();
    MMp800900_Perp->Write();

    Phip_430MeVCM1_Perp->Write();
    Phip_450MeVCM1_Perp->Write();
    Phip_470MeVCM1_Perp->Write();
    Phip_490MeVCM1_Perp->Write();
    Phip_510MeVCM1_Perp->Write();
    Phip_530MeVCM1_Perp->Write();
    Phip_550MeVCM1_Perp->Write();
    Phip_570MeVCM1_Perp->Write();
    Phip_590MeVCM1_Perp->Write();
    Phip_610MeVCM1_Perp->Write();

    Phip_430MeVCM2_Perp->Write();
    Phip_450MeVCM2_Perp->Write();
    Phip_470MeVCM2_Perp->Write();
    Phip_490MeVCM2_Perp->Write();
    Phip_510MeVCM2_Perp->Write();
    Phip_530MeVCM2_Perp->Write();
    Phip_550MeVCM2_Perp->Write();
    Phip_570MeVCM2_Perp->Write();
    Phip_590MeVCM2_Perp->Write();
    Phip_610MeVCM2_Perp->Write();

    Phip_430MeVCM3_Perp->Write();
    Phip_450MeVCM3_Perp->Write();
    Phip_470MeVCM3_Perp->Write();
    Phip_490MeVCM3_Perp->Write();
    Phip_510MeVCM3_Perp->Write();
    Phip_530MeVCM3_Perp->Write();
    Phip_550MeVCM3_Perp->Write();
    Phip_570MeVCM3_Perp->Write();
    Phip_590MeVCM3_Perp->Write();
    Phip_610MeVCM3_Perp->Write();

    Phip_430MeVCM4_Perp->Write();
    Phip_450MeVCM4_Perp->Write();
    Phip_470MeVCM4_Perp->Write();
    Phip_490MeVCM4_Perp->Write();
    Phip_510MeVCM4_Perp->Write();
    Phip_530MeVCM4_Perp->Write();
    Phip_550MeVCM4_Perp->Write();
    Phip_570MeVCM4_Perp->Write();
    Phip_590MeVCM4_Perp->Write();
    Phip_610MeVCM4_Perp->Write();

    Phip_430MeVCM5_Perp->Write();
    Phip_450MeVCM5_Perp->Write();
    Phip_470MeVCM5_Perp->Write();
    Phip_490MeVCM5_Perp->Write();
    Phip_510MeVCM5_Perp->Write();
    Phip_530MeVCM5_Perp->Write();
    Phip_550MeVCM5_Perp->Write();
    Phip_570MeVCM5_Perp->Write();
    Phip_590MeVCM5_Perp->Write();
    Phip_610MeVCM5_Perp->Write();

    Phip_430MeVCM6_Perp->Write();
    Phip_450MeVCM6_Perp->Write();
    Phip_470MeVCM6_Perp->Write();
    Phip_490MeVCM6_Perp->Write();
    Phip_510MeVCM6_Perp->Write();
    Phip_530MeVCM6_Perp->Write();
    Phip_550MeVCM6_Perp->Write();
    Phip_570MeVCM6_Perp->Write();
    Phip_590MeVCM6_Perp->Write();
    Phip_610MeVCM6_Perp->Write();

    Phip_430MeVCM7_Perp->Write();
    Phip_450MeVCM7_Perp->Write();
    Phip_470MeVCM7_Perp->Write();
    Phip_490MeVCM7_Perp->Write();
    Phip_510MeVCM7_Perp->Write();
    Phip_530MeVCM7_Perp->Write();
    Phip_550MeVCM7_Perp->Write();
    Phip_570MeVCM7_Perp->Write();
    Phip_590MeVCM7_Perp->Write();
    Phip_610MeVCM7_Perp->Write();

    Phip_430MeVCM8_Perp->Write();
    Phip_450MeVCM8_Perp->Write();
    Phip_470MeVCM8_Perp->Write();
    Phip_490MeVCM8_Perp->Write();
    Phip_510MeVCM8_Perp->Write();
    Phip_530MeVCM8_Perp->Write();
    Phip_550MeVCM8_Perp->Write();
    Phip_570MeVCM8_Perp->Write();
    Phip_590MeVCM8_Perp->Write();
    Phip_610MeVCM8_Perp->Write();

    Phip_430MeVCM9_Perp->Write();
    Phip_450MeVCM9_Perp->Write();
    Phip_470MeVCM9_Perp->Write();
    Phip_490MeVCM9_Perp->Write();
    Phip_510MeVCM9_Perp->Write();
    Phip_530MeVCM9_Perp->Write();
    Phip_550MeVCM9_Perp->Write();
    Phip_570MeVCM9_Perp->Write();
    Phip_590MeVCM9_Perp->Write();
    Phip_610MeVCM9_Perp->Write();

    Phip_430MeVCM10_Perp->Write();
    Phip_450MeVCM10_Perp->Write();
    Phip_470MeVCM10_Perp->Write();
    Phip_490MeVCM10_Perp->Write();
    Phip_510MeVCM10_Perp->Write();
    Phip_530MeVCM10_Perp->Write();
    Phip_550MeVCM10_Perp->Write();
    Phip_570MeVCM10_Perp->Write();
    Phip_590MeVCM10_Perp->Write();
    Phip_610MeVCM10_Perp->Write();

    ThetanDist_Perp->Write();
    ThetanRecDist_Perp->Write();
    ThetanDiffDist_Perp->Write();
    ThetanDiffZp_Perp->Write();
    ThetanCorrDist_Perp->Write();
    ThetanCorrDiffDist_Perp->Write();
    ThetanCorrRecDiffDist_Perp->Write();
    ThetanCorrDiffZp_Perp->Write();

    E_dE_Perp->Write();;
    KinEp_dE_Perp->Write();
    E_KinEp_Perp->Write();

    ThetaRecPiDiff_Perp->Write();
    ThetanThetaRecPi_Perp->Write();
    ThetanThetaRecPiDiff_Perp->Write();
    ThetaRecPDiff_Perp->Write();
    ThetanThetaRecP_Perp->Write();
    ThetanThetaRecPDiff_Perp->Write();

    DeutKinPiKin_Perp->Write();

    f2.Write();
}
