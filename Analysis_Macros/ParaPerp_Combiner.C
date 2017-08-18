#include "./includes_ParaPerpCombiner.h"

void ParaPerp_Combiner(){

  TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/Para/Physics_Total_Para_20_16_08_17.root"); // Open latest Para file

    TH1D* time_Para = (TH1D*)f->Get("time")->Clone();
    time_Para->SetName("time_Para");
    TH1D* time_cut_Para = (TH1D*)f->Get("time_cut")->Clone();
    time_cut_Para->SetName("time_cut_Para");

    TH1D* Eg_Para = (TH1D*)f->Get("Eg")->Clone();
    Eg_Para->SetName("Eg_Para");

    TH1D* MMpEpCorrected_Para = (TH1D*)f->Get("MMpEpCorrected")->Clone();
    MMpEpCorrected_Para->SetName("MMpEpCorrected_Para");
    TH1D* OAngle_Para = (TH1D*)f->Get("OAngle")->Clone();
    OAngle_Para->SetName("OAngle_Para");

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

    TH1D* ThetanDist_Para = (TH1D*)f->Get("ThetanDist")->Clone();
    ThetanDist_Para->SetName("ThetanDist_Para");

    TH2D* E_dE_Para = (TH2D*)f->Get("E_dE")->Clone();
    E_dE_Para->SetName("E_dE_Para");

    TH2D* DeutKinPiKin_Para = (TH2D*)f->Get("DeutKinPiKin")->Clone();
    DeutKinPiKin_Para->SetName("DeutKinPiKin_Para");

    ////////////////////////////////////////////////////
    ////////////////////////////////////////////////////
    ///////////////// PARA DONE ////////////////////////
    ////////////////////////////////////////////////////
    ////////////////////////////////////////////////////

    TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/Perp/Physics_Total_Perp_20_16_08_17.root"); // Open latest Perp file

    TH1D* time_Perp = (TH1D*)f->Get("time")->Clone();
    time_Perp->SetName("time_Perp");
    TH1D* time_cut_Perp = (TH1D*)f->Get("time_cut")->Clone();
    time_cut_Perp->SetName("time_cut_Perp");

    TH1D* Eg_Perp = (TH1D*)f->Get("Eg")->Clone();
    Eg_Perp->SetName("Eg_Perp");

    TH1D* MMpEpCorrected_Perp = (TH1D*)f->Get("MMpEpCorrected")->Clone();
    MMpEpCorrected_Perp->SetName("MMpEpCorrected_Perp");
    TH1D* OAngle_Perp = (TH1D*)f->Get("OAngle")->Clone();
    OAngle_Perp->SetName("OAngle_Perp");

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

    TH1D* ThetanDist_Perp = (TH1D*)f->Get("ThetanDist")->Clone();
    ThetanDist_Perp->SetName("ThetanDist_Perp");

    TH2D* E_dE_Perp = (TH2D*)f->Get("E_dE")->Clone();
    E_dE_Perp->SetName("E_dE_Perp");

    TH2D* DeutKinPiKin_Perp = (TH2D*)f->Get("DeutKinPiKin")->Clone();
    DeutKinPiKin_Perp->SetName("DeutKinPiKin_Perp");

    ////////////////////////////////////////////////////
    ////////////////////////////////////////////////////
    ///////////////// PERP DONE ////////////////////////
    ////////////////////////////////////////////////////
    ////////////////////////////////////////////////////

    TFile f2("ParaPerp_Total_20_Combined.root", "RECREATE");

    time_Para->Write();
    time_cut_Para->Write();

    Eg_Para->Write();
    MMpEpCorrected_Para->Write();
    OAngle_Para->Write();

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

    ThetanDist_Para->Write();

    E_dE_Para->Write();;
    DeutKinPiKin_Para->Write();

    ////////////////////////////////////////////////////
    ////////////////////////////////////////////////////
    ///////////////// PARA DONE ////////////////////////
    ////////////////////////////////////////////////////
    ////////////////////////////////////////////////////

    time_Perp->Write();
    time_cut_Perp->Write();

    Eg_Perp->Write();
    MMpEpCorrected_Perp->Write();
    OAngle_Perp->Write();

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

    ThetanDist_Perp->Write();

    E_dE_Perp->Write();;
    DeutKinPiKin_Perp->Write();

    f2.Write();
}
