#include "./includes_ParaPerp_Combiner_Unpolarised.h"

void ParaPerp_Combiner_Unpolarised() {

    TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/Para/Physics_Total_Para_12_15_5_17.root"); // Open latest Para file

    TH1D* time_Para = (TH1D*)f->Get("time")->Clone();
    time_Para->SetName("time_Para");
    TH1D* time_cut_Para = (TH1D*)f->Get("time_cut")->Clone();
    time_cut_Para->SetName("time_cut_Para");
    TH1D* Eg_Para = (TH1D*)f->Get("Eg")->Clone();
    Eg_Para->SetName("Eg_Para");
    TH1D* WCPhiDifference_Para = (TH1D*)f->Get("WCPhiDifference")->Clone();
    WCPhiDifference_Para->SetName("WCPhiDifference_Para");
    TH1D* EpKin_Para = (TH1D*)f->Get("EpKin")->Clone();
    EpKin_Para->SetName("EpKin_Para");
    TH1D* EpCorrected_Para = (TH1D*)f->Get("EpCorrected")->Clone();
    EpCorrected_Para->SetName("EpCorrected_Para");
    TH1D* OAngle_Para = (TH1D*)f->Get("OAngle")->Clone();
    OAngle_Para->SetName("OAngle_Para");
    TH1D* WCZnRecon_Para = (TH1D*)f->Get("WCZnRecon")->Clone();
    WCZnRecon_Para->SetName("WCZnRecon_Para");
    TH1D* Theta_Scattered_Para = (TH1D*)f->Get("Theta_Scattered")->Clone();
    Theta_Scattered_Para->SetName("Theta_Scattered_Para");
    TH1D* Phi_Scattered_Para = (TH1D*)f->Get("Phi_Scattered")->Clone();
    Phi_Scattered_Para->SetName("Phi_Scattered_Para");
    TH1D* EpKinEpCorrDiff_Para = (TH1D*)f->Get("EpKinEpCorrDiff")->Clone();
    EpKinEpCorrDiff_Para->SetName("EpKinEpCorrDiff_Para");
    TH1D* EpEpCorrDiff_Para = (TH1D*)f->Get("EpEpCorrDiff")->Clone();
    EpEpCorrDiff_Para->SetName("EpEpCorrDiff_Para");
    TH1D* MMpEpCorrected_Para = (TH1D*)f->Get("MMpEpCorrected")->Clone();
    MMpEpCorrected_Para->SetName("MMpEpCorrected_Para");
    TH1D* ZpDist_Para = (TH1D*)f->Get("ZpDist")->Clone();
    ZpDist_Para->SetName("ZpDist_Para");
    TH1D* ZpPhiScatNeg180_Para = (TH1D*)f->Get("ZpPhiScatNeg180")->Clone();
    ZpPhiScatNeg180_Para->SetName("ZpPhiScatNeg180_Para");
    TH1D* ZpPhiScat0_Para = (TH1D*)f->Get("ZpPhiScat0")->Clone();
    ZpPhiScat0_Para->SetName("ZpPhiScat0_Para");
    TH1D* ZpPhiScatPos180_Para = (TH1D*)f->Get("ZpPhiScatPos180")->Clone();
    ZpPhiScatPos180_Para->SetName("ZpPhiScatPos180_Para");
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
    MMp800900_Para->SetName("MM8600900_Para");
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
    TH2D* E_dE_Para = (TH2D*)f->Get("E_dE")->Clone();
    E_dE_Para->SetName("E_dE_Para");
    TH2D* KinEp_dE_Para = (TH2D*)f->Get("KinEp_dE")->Clone();
    KinEp_dE_Para->SetName("KinEp_dE_Para");
    TH2D* ThetaScPhiSc_Para = (TH2D*)f->Get("ThetaScPhiSc")->Clone();
    ThetaScPhiSc_Para->SetName("ThetaScPhiSc_Para");
    TH2D* E_KinEp_Para = (TH2D*)f->Get("E_KinEp")->Clone();
    E_KinEp_Para->SetName("E_KinEp_Para");
    TH2D* PhinDiffWCZRec_Para = (TH2D*)f->Get("PhinDiffWCZRec")->Clone();
    PhinDiffWCZRec_Para->SetName(" PhinDiffWCZRec_Para");

    TH1D* Phi_Scattered_410MeV_Para = (TH1D*)f->Get("Phi_Scattered_410MeV")->Clone();
    Phi_Scattered_410MeV_Para->SetName("Phi_Scattered_410MeV_Para");
    TH1D* Phi_Scattered_430MeV_Para = (TH1D*)f->Get("Phi_Scattered_430MeV")->Clone();
    Phi_Scattered_430MeV_Para->SetName("Phi_Scattered_430MeV_Para");
    TH1D* Phi_Scattered_450MeV_Para = (TH1D*)f->Get("Phi_Scattered_450MeV")->Clone();
    Phi_Scattered_450MeV_Para->SetName("Phi_Scattered_450MeV_Para");
    TH1D* Phi_Scattered_470MeV_Para = (TH1D*)f->Get("Phi_Scattered_470MeV")->Clone();
    Phi_Scattered_470MeV_Para->SetName("Phi_Scattered_470MeV_Para");
    TH1D* Phi_Scattered_490MeV_Para = (TH1D*)f->Get("Phi_Scattered_490MeV")->Clone();
    Phi_Scattered_490MeV_Para->SetName("Phi_Scattered_490MeV_Para");
    TH1D* Phi_Scattered_510MeV_Para = (TH1D*)f->Get("Phi_Scattered_510MeV")->Clone();
    Phi_Scattered_510MeV_Para->SetName("Phi_Scattered_510MeV_Para");
    TH1D* Phi_Scattered_530MeV_Para = (TH1D*)f->Get("Phi_Scattered_530MeV")->Clone();
    Phi_Scattered_530MeV_Para->SetName("Phi_Scattered_530MeV_Para");
    TH1D* Phi_Scattered_550MeV_Para = (TH1D*)f->Get("Phi_Scattered_550MeV")->Clone();
    Phi_Scattered_550MeV_Para->SetName("Phi_Scattered_550MeV_Para");
    TH1D* Phi_Scattered_570MeV_Para = (TH1D*)f->Get("Phi_Scattered_570MeV")->Clone();
    Phi_Scattered_570MeV_Para->SetName("Phi_Scattered_570MeV_Para");
    TH1D* Phi_Scattered_590MeV_Para = (TH1D*)f->Get("Phi_Scattered_590MeV")->Clone();
    Phi_Scattered_590MeV_Para->SetName("Phi_Scattered_590MeV_Para");
    TH1D* Phi_Scattered_610MeV_Para = (TH1D*)f->Get("Phi_Scattered_610MeV")->Clone();
    Phi_Scattered_610MeV_Para->SetName("Phi_Scattered_610MeV_Para");
    TH1D* Phi_Scattered_630MeV_Para = (TH1D*)f->Get("Phi_Scattered_630MeV")->Clone();
    Phi_Scattered_630MeV_Para->SetName("Phi_Scattered_630MeV_Para");

    ////////////////////////////////////////////////////
    ////////////////////////////////////////////////////
    ///////////////// PARA DONE ////////////////////////
    ////////////////////////////////////////////////////
    ////////////////////////////////////////////////////

    TFile *f1 = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/Perp/Physics_Total_Perp_12_15_5_17.root"); // Open latest Para file

    TH1D* time_Perp = (TH1D*)f1->Get("time")->Clone();
    time_Perp->SetName("time_Perp");
    TH1D* time_cut_Perp = (TH1D*)f1->Get("time_cut")->Clone();
    time_cut_Perp->SetName("time_cut_Perp");
    TH1D* Eg_Perp = (TH1D*)f1->Get("Eg")->Clone();
    Eg_Perp->SetName("Eg_Perp");
    TH1D* WCPhiDifference_Perp = (TH1D*)f1->Get("WCPhiDifference")->Clone();
    WCPhiDifference_Perp->SetName("WCPhiDifference_Perp");
    TH1D* EpKin_Perp = (TH1D*)f1->Get("EpKin")->Clone();
    EpKin_Perp->SetName("EpKin_Perp");
    TH1D* EpCorrected_Perp = (TH1D*)f1->Get("EpCorrected")->Clone();
    EpCorrected_Perp->SetName("EpCorrected_Perp");
    TH1D* OAngle_Perp = (TH1D*)f1->Get("OAngle")->Clone();
    OAngle_Perp->SetName("OAngle_Perp");
    TH1D* WCZnRecon_Perp = (TH1D*)f1->Get("WCZnRecon")->Clone();
    WCZnRecon_Perp->SetName("WCZnRecon_Perp");
    TH1D* Theta_Scattered_Perp = (TH1D*)f1->Get("Theta_Scattered")->Clone();
    Theta_Scattered_Perp->SetName("Theta_Scattered_Perp");
    TH1D* Phi_Scattered_Perp = (TH1D*)f1->Get("Phi_Scattered")->Clone();
    Phi_Scattered_Perp->SetName("Phi_Scattered_Perp");
    TH1D* EpKinEpCorrDiff_Perp = (TH1D*)f1->Get("EpKinEpCorrDiff")->Clone();
    EpKinEpCorrDiff_Perp->SetName("EpKinEpCorrDiff_Perp");
    TH1D* EpEpCorrDiff_Perp = (TH1D*)f1->Get("EpEpCorrDiff")->Clone();
    EpEpCorrDiff_Perp->SetName("EpEpCorrDiff_Perp");
    TH1D* MMpEpCorrected_Perp = (TH1D*)f1->Get("MMpEpCorrected")->Clone();
    MMpEpCorrected_Perp->SetName("MMpEpCorrected_Perp");
    TH1D* ZpDist_Perp = (TH1D*)f1->Get("ZpDist")->Clone();
    ZpDist_Perp->SetName("ZpDist_Perp");
    TH1D* ZpPhiScatNeg180_Perp = (TH1D*)f1->Get("ZpPhiScatNeg180")->Clone();
    ZpPhiScatNeg180_Perp->SetName("ZpPhiScatNeg180_Perp");
    TH1D* ZpPhiScat0_Perp = (TH1D*)f1->Get("ZpPhiScat0")->Clone();
    ZpPhiScat0_Perp->SetName("ZpPhiScat0_Perp");
    TH1D* ZpPhiScatPos180_Perp = (TH1D*)f1->Get("ZpPhiScatPos180")->Clone();
    ZpPhiScatPos180_Perp->SetName("ZpPhiScatPos180_Perp");
    TH1D* MMp200300_Perp = (TH1D*)f1->Get("MMp200300")->Clone();
    MMp200300_Perp->SetName("MMp200300_Perp");
    TH1D* MMp300400_Perp = (TH1D*)f1->Get("MMp300400")->Clone();
    MMp300400_Perp->SetName("MMp300400_Perp");
    TH1D* MMp400500_Perp = (TH1D*)f1->Get("MMp400500")->Clone();
    MMp400500_Perp->SetName("MMp400500_Perp");
    TH1D* MMp500600_Perp = (TH1D*)f1->Get("MMp500600")->Clone();
    MMp500600_Perp->SetName("MMp500600_Perp");
    TH1D* MMp600700_Perp = (TH1D*)f1->Get("MMp600700")->Clone();
    MMp600700_Perp->SetName("MMp600700_Perp");
    TH1D* MMp700800_Perp = (TH1D*)f1->Get("MMp700800")->Clone();
    MMp700800_Perp->SetName("MMp700800_Perp");
    TH1D* MMp800900_Perp = (TH1D*)f1->Get("MMp800900")->Clone();
    MMp800900_Perp->SetName("MM8600900_Perp");
    TH1D* ThetaRecPiDiff_Perp = (TH1D*)f1->Get("ThetaRecPiDiff")->Clone();
    ThetaRecPiDiff_Perp->SetName("ThetaRecPiDiff_Perp");
    TH2D* ThetanThetaRecPi_Perp = (TH2D*)f1->Get("ThetanThetaRecPi")->Clone();
    ThetanThetaRecPi_Perp->SetName("ThetanThetaRecPi_Perp");
    TH2D* ThetanThetaRecPiDiff_Perp = (TH2D*)f1->Get("ThetanThetaRecPiDiff")->Clone();
    ThetanThetaRecPiDiff_Perp->SetName("ThetanThetaRecPiDiff_Perp");
    TH1D* ThetaRecPDiff_Perp = (TH1D*)f1->Get("ThetaRecPDiff")->Clone();
    ThetaRecPDiff_Perp->SetName("ThetaRecPDiff_Perp");
    TH2D* ThetanThetaRecP_Perp = (TH2D*)f1->Get("ThetanThetaRecP")->Clone();
    ThetanThetaRecP_Perp->SetName("ThetanThetaRecP_Perp");
    TH2D* ThetanThetaRecPDiff_Perp = (TH2D*)f1->Get("ThetanThetaRecPDiff")->Clone();
    ThetanThetaRecPDiff_Perp->SetName("ThetanThetaRecPDiff_Perp");
    TH2D* E_dE_Perp = (TH2D*)f1->Get("E_dE")->Clone();
    E_dE_Perp->SetName("E_dE_Perp");
    TH2D* KinEp_dE_Perp = (TH2D*)f1->Get("KinEp_dE")->Clone();
    KinEp_dE_Perp->SetName("KinEp_dE_Perp");
    TH2D* ThetaScPhiSc_Perp = (TH2D*)f1->Get("ThetaScPhiSc")->Clone();
    ThetaScPhiSc_Perp->SetName("ThetaScPhiSc_Perp");
    TH2D* E_KinEp_Perp = (TH2D*)f1->Get("E_KinEp")->Clone();
    E_KinEp_Perp->SetName("E_KinEp_Perp");
    TH2D* PhinDiffWCZRec_Perp = (TH2D*)f1->Get("PhinDiffWCZRec")->Clone();
    PhinDiffWCZRec_Perp->SetName(" PhinDiffWCZRec_Perp");

    TH1D* Phi_Scattered_410MeV_Perp = (TH1D*)f1->Get("Phi_Scattered_410MeV")->Clone();
    Phi_Scattered_410MeV_Perp->SetName("Phi_Scattered_410MeV_Perp");
    TH1D* Phi_Scattered_430MeV_Perp = (TH1D*)f1->Get("Phi_Scattered_430MeV")->Clone();
    Phi_Scattered_430MeV_Perp->SetName("Phi_Scattered_430MeV_Perp");
    TH1D* Phi_Scattered_450MeV_Perp = (TH1D*)f1->Get("Phi_Scattered_450MeV")->Clone();
    Phi_Scattered_450MeV_Perp->SetName("Phi_Scattered_450MeV_Perp");
    TH1D* Phi_Scattered_470MeV_Perp = (TH1D*)f1->Get("Phi_Scattered_470MeV")->Clone();
    Phi_Scattered_470MeV_Perp->SetName("Phi_Scattered_470MeV_Perp");
    TH1D* Phi_Scattered_490MeV_Perp = (TH1D*)f1->Get("Phi_Scattered_490MeV")->Clone();
    Phi_Scattered_490MeV_Perp->SetName("Phi_Scattered_490MeV_Perp");
    TH1D* Phi_Scattered_510MeV_Perp = (TH1D*)f1->Get("Phi_Scattered_510MeV")->Clone();
    Phi_Scattered_510MeV_Perp->SetName("Phi_Scattered_510MeV_Perp");
    TH1D* Phi_Scattered_530MeV_Perp = (TH1D*)f1->Get("Phi_Scattered_530MeV")->Clone();
    Phi_Scattered_530MeV_Perp->SetName("Phi_Scattered_530MeV_Perp");
    TH1D* Phi_Scattered_550MeV_Perp = (TH1D*)f1->Get("Phi_Scattered_550MeV")->Clone();
    Phi_Scattered_550MeV_Perp->SetName("Phi_Scattered_550MeV_Perp");
    TH1D* Phi_Scattered_570MeV_Perp = (TH1D*)f1->Get("Phi_Scattered_570MeV")->Clone();
    Phi_Scattered_570MeV_Perp->SetName("Phi_Scattered_570MeV_Perp");
    TH1D* Phi_Scattered_590MeV_Perp = (TH1D*)f1->Get("Phi_Scattered_590MeV")->Clone();
    Phi_Scattered_590MeV_Perp->SetName("Phi_Scattered_590MeV_Perp");
    TH1D* Phi_Scattered_610MeV_Perp = (TH1D*)f1->Get("Phi_Scattered_610MeV")->Clone();
    Phi_Scattered_610MeV_Perp->SetName("Phi_Scattered_610MeV_Perp");
    TH1D* Phi_Scattered_630MeV_Perp = (TH1D*)f1->Get("Phi_Scattered_630MeV")->Clone();
    Phi_Scattered_630MeV_Perp->SetName("Phi_Scattered_630MeV_Perp");


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

    TList *list1 = new TList;
    list1->Add(time_Para);
    time_Perp->Scale(ScaleFactor);
    list1->Add(time_Perp);
    time = new TH1D("time", "time", 1400, -700, 700);
    time->Merge(list1);

    TList *list2 = new TList;
    list2->Add(time_cut_Para);
    time_cut_Perp->Scale(ScaleFactor);
    list2->Add(time_cut_Perp);
    time_cut = new TH1D("time_cut", "time_cut", 1400, -700, 700);
    time_cut->Merge(list2);

    TList *list3 = new TList;
    list3->Add(Eg_Para);
    Eg_Perp->Scale(ScaleFactor);
    list3->Add(Eg_Perp);
    Eg = new TH1D("Eg", "Eg", 200, 100, 1600);
    Eg->Merge(list3);

    TList *list4 = new TList;
    list4->Add(WCPhiDifference_Para);
    WCPhiDifference_Perp->Scale(ScaleFactor);
    list4->Add(WCPhiDifference_Perp);
    WCPhiDifference = new TH1D("WCPhiDifference", "WCPhiDifference", 180, 0, 360);
    WCPhiDifference->Merge(list4);

    TList *list5 = new TList;
    list5->Add(EpKin_Para);
    EpKin_Perp->Scale(ScaleFactor);
    list5->Add(EpKin_Perp);
    EpKin = new TH1D("EpKin", "EpKin", 100, 0, 500);
    EpKin->Merge(list5);

    TList *list6 = new TList;
    list6->Add(EpCorrected_Para);
    EpCorrected_Perp->Scale(ScaleFactor);
    list6->Add(EpCorrected_Perp);
    EpCorrected = new TH1D("EpCorrected", "EpCorrected", 100, 0, 500);
    EpCorrected->Merge(list6);

    TList *list7 = new TList;
    list7->Add(OAngle_Para);
    OAngle_Perp->Scale(ScaleFactor);
    list7->Add(OAngle_Perp);
    OAngle = new TH1D("OAngle", "OAngle", 180, 0, 180);
    OAngle->Merge(list7);

    TList *list8 = new TList;
    list8->Add(WCZnRecon_Para);
    WCZnRecon_Perp->Scale(ScaleFactor);
    list8->Add(WCZnRecon_Perp);
    WCZnRecon = new TH1D("WCZnRecon", "WCZnRecon", 200, 0, 400);
    WCZnRecon->Merge(list8);

    TList *list9 = new TList;
    list9->Add(Theta_Scattered_Para);
    Theta_Scattered_Perp->Scale(ScaleFactor);
    list9->Add(Theta_Scattered_Perp);
    ThetaSc = new TH1D("Theta_Scattered", "Theta_Scattered", 180, 0, 180);
    ThetaSc->Merge(list9);

    TList *list10 = new TList;
    list10->Add(Phi_Scattered_Para);
    Phi_Scattered_Perp->Scale(ScaleFactor);
    list10->Add(Phi_Scattered_Perp);
    PhiSc = new TH1D("Phi_Scattered", "Phi_Scattered", 90, -180, 180);
    PhiSc->Merge(list10);

    TList *list11 = new TList;
    list11->Add(EpKinEpCorrDiff_Para);
    EpKinEpCorrDiff_Perp->Scale(ScaleFactor);
    list11->Add(EpKinEpCorrDiff_Perp);
    EpKinEpCorrDiff = new TH1D("EpKinEpCorrDiff", "EpKinEpCorrDiff", 300, -300, 300);
    EpKinEpCorrDiff->Merge(list11);

    TList *list12 = new TList;
    list12->Add(EpEpCorrDiff_Para);
    EpEpCorrDiff_Perp->Scale(ScaleFactor);
    list12->Add(EpEpCorrDiff_Perp);
    EpEpCorrDiff = new TH1D("EpEpCorrDiff", "EpEpCorrDiff", 200, 0, 200);
    EpEpCorrDiff->Merge(list12);

    TList *list13 = new TList;
    list13->Add(MMpEpCorrected_Para);
    MMpEpCorrected_Perp->Scale(ScaleFactor);
    list13->Add(MMpEpCorrected_Perp);
    MMpEpCorrected = new TH1D("MMpEpCorrected", "MMpEpCorrected", 400, 0, 2000);
    MMpEpCorrected->Merge(list13);

    TList *list14 = new TList;
    list14->Add(ZpDist_Para);
    ZpDist_Perp->Scale(ScaleFactor);
    list14->Add(ZpDist_Perp);
    ZpDist = new TH1D("ZpDist", "ZpDist", 200, -400, 400);
    ZpDist->Merge(list14);

    TList *list15 = new TList;
    list15->Add(ZpPhiScatNeg180_Para);
    ZpPhiScatNeg180_Perp->Scale(ScaleFactor);
    list15->Add(ZpPhiScatNeg180_Perp);
    ZpPhiScatNeg180 = new TH1D("ZpPhiScatNeg180", "ZpPhiScatNeg180", 200, -200, 200);
    ZpPhiScatNeg180->Merge(list15);

    TList *list16 = new TList;
    list16->Add(ZpPhiScat0_Para);
    ZpPhiScat0_Perp->Scale(ScaleFactor);
    list16->Add(ZpPhiScat0_Perp);
    ZpPhiScat0 = new TH1D("ZpPhiScat0", "ZpPhiScat0", 200, -200, 200);
    ZpPhiScat0->Merge(list16);

    TList *list17 = new TList;
    list17->Add(ZpPhiScatPos180_Para);
    ZpPhiScatPos180_Perp->Scale(ScaleFactor);
    list17->Add(ZpPhiScatPos180_Perp);
    ZpPhiScatPos180 = new TH1D("ZpPhiScatPos180", "ZpPhiScatPos180", 200, -200, 200);
    ZpPhiScatPos180->Merge(list17);

    TList *list18 = new TList;
    list18->Add(MMp200300_Para);
    MMp200300_Perp->Scale(ScaleFactor);
    list18->Add(MMp200300_Perp);
    MMp200300 = new TH1D("MMp200300", "MMp200300", 400, 0, 2000);
    MMp200300->Merge(list18);

    TList *list19 = new TList;
    list19->Add(MMp300400_Para);
    MMp300400_Perp->Scale(ScaleFactor);
    list19->Add(MMp300400_Perp);
    MMp300400 = new TH1D("MMp300400", "MMp300400", 400, 0, 2000);
    MMp300400->Merge(list19);

    TList *list20 = new TList;
    list20->Add(MMp400500_Para);
    MMp400500_Perp->Scale(ScaleFactor);
    list20->Add(MMp400500_Perp);
    MMp400500 = new TH1D("MMp400500", "MMp400500", 400, 0, 2000);
    MMp400500->Merge(list20);

    TList *list21 = new TList;
    list21->Add(MMp500600_Para);
    MMp500600_Perp->Scale(ScaleFactor);
    list21->Add(MMp500600_Perp);
    MMp500600 = new TH1D("MMp500600", "MMp500600", 400, 0, 2000);
    MMp500600->Merge(list21);

    TList *list22 = new TList;
    list22->Add(MMp600700_Para);
    MMp600700_Perp->Scale(ScaleFactor);
    list22->Add(MMp600700_Perp);
    MMp600700 = new TH1D("MMp600700", "MMp600700", 400, 0, 2000);
    MMp600700->Merge(list22);

    TList *list23 = new TList;
    list23->Add(MMp700800_Para);
    MMp700800_Perp->Scale(ScaleFactor);
    list23->Add(MMp700800_Perp);
    MMp700800 = new TH1D("MMp700800", "MMp700800", 400, 0, 2000);
    MMp700800->Merge(list23);

    TList *list24 = new TList;
    list24->Add(MMp800900_Para);
    MMp800900_Perp->Scale(ScaleFactor);
    list24->Add(MMp800900_Perp);
    MMp800900 = new TH1D("MMp800900", "MMp800900", 400, 0, 2000);
    MMp800900->Merge(list24);

    TList *list25 = new TList;
    list25->Add(Phi_Scattered_410MeV_Para);
    Phi_Scattered_410MeV_Perp->Scale(ScaleFactor);
    list25->Add(Phi_Scattered_410MeV_Perp);
    PhiSc410 = new TH1D("PhiSc410", "PhiSc410", 2, -180, 180);
    PhiSc410->Merge(list25);

    TList *list26 = new TList;
    list26->Add(Phi_Scattered_430MeV_Para);
    Phi_Scattered_430MeV_Perp->Scale(ScaleFactor);
    list26->Add(Phi_Scattered_430MeV_Perp);
    PhiSc430 = new TH1D("PhiSc430", "PhiSc430", 2, -180, 180);
    PhiSc430->Merge(list26);

    TList *list27 = new TList;
    list27->Add(Phi_Scattered_450MeV_Para);
    Phi_Scattered_450MeV_Perp->Scale(ScaleFactor);
    list27->Add(Phi_Scattered_450MeV_Perp);
    PhiSc450 = new TH1D("PhiSc450", "PhiSc450", 2, -180, 180);
    PhiSc450->Merge(list27);

    TList *list28 = new TList;
    list28->Add(Phi_Scattered_470MeV_Para);
    Phi_Scattered_470MeV_Perp->Scale(ScaleFactor);
    list28->Add(Phi_Scattered_470MeV_Perp);
    PhiSc470 = new TH1D("PhiSc470", "PhiSc470", 2, -180, 180);
    PhiSc470->Merge(list28);

    TList *list28 = new TList;
    list28->Add(Phi_Scattered_490MeV_Para);
    Phi_Scattered_490MeV_Perp->Scale(ScaleFactor);
    list28->Add(Phi_Scattered_490MeV_Perp);
    PhiSc490 = new TH1D("PhiSc490", "PhiSc490", 2, -180, 180);
    PhiSc490->Merge(list28);

    TList *list29 = new TList;
    list29->Add(Phi_Scattered_510MeV_Para);
    Phi_Scattered_510MeV_Perp->Scale(ScaleFactor);
    list29->Add(Phi_Scattered_510MeV_Perp);
    PhiSc510 = new TH1D("PhiSc510", "PhiSc510", 2, -180, 180);
    PhiSc510->Merge(list29);

    TList *list30 = new TList;
    list30->Add(Phi_Scattered_530MeV_Para);
    Phi_Scattered_530MeV_Perp->Scale(ScaleFactor);
    list30->Add(Phi_Scattered_530MeV_Perp);
    PhiSc530 = new TH1D("PhiSc530", "PhiSc530", 2, -180, 180);
    PhiSc530->Merge(list30);

    TList *list31 = new TList;
    list31->Add(Phi_Scattered_550MeV_Para);
    Phi_Scattered_550MeV_Perp->Scale(ScaleFactor);
    list31->Add(Phi_Scattered_550MeV_Perp);
    PhiSc550 = new TH1D("PhiSc550", "PhiSc550", 2, -180, 180);
    PhiSc550->Merge(list31);

    TList *list32 = new TList;
    list32->Add(Phi_Scattered_570MeV_Para);
    Phi_Scattered_570MeV_Perp->Scale(ScaleFactor);
    list32->Add(Phi_Scattered_570MeV_Perp);
    PhiSc570 = new TH1D("PhiSc570", "PhiSc570", 2, -180, 180);
    PhiSc570->Merge(list32);

    TList *list33 = new TList;
    list33->Add(Phi_Scattered_590MeV_Para);
    Phi_Scattered_590MeV_Perp->Scale(ScaleFactor);
    list33->Add(Phi_Scattered_590MeV_Perp);
    PhiSc590 = new TH1D("PhiSc590", "PhiSc590", 2, -180, 180);
    PhiSc590->Merge(list33);

    TList *list34 = new TList;
    list34->Add(Phi_Scattered_610MeV_Para);
    Phi_Scattered_610MeV_Perp->Scale(ScaleFactor);
    list34->Add(Phi_Scattered_610MeV_Perp);
    PhiSc610 = new TH1D("PhiSc610", "PhiSc610", 2, -180, 180);
    PhiSc610->Merge(list34);

    TList *list35 = new TList;
    list35->Add(Phi_Scattered_630MeV_Para);
    Phi_Scattered_630MeV_Perp->Scale(ScaleFactor);
    list35->Add(Phi_Scattered_630MeV_Perp);
    PhiSc630 = new TH1D("PhiSc630", "PhiSc630", 2, -180, 180);
    PhiSc630->Merge(list35);

    TList *list36 = new TList;
    list36->Add(ThetaRecPiDiff_Para);
    ThetaRecPiDiff_Perp->Scale(ScaleFactor);
    list36->Add(ThetaRecPiDiff_Perp);
    ThetaRecPiDiff = new TH1D("ThetaRecPiDiff", "ThetaRecPiDiff", 200, 0, 180);
    ThetaRecPiDiff->Merge(list36);

    TList *list37 = new TList;
    list37->Add(ThetanThetaRecPi_Para);
    ThetanThetaRecPi_Perp->Scale(ScaleFactor);
    list37->Add(ThetanThetaRecPi_Perp);
    ThetanThetaRecPi = new TH2D("ThetanThetaRecPi", "ThetanThetaRecPi", 100, 0, 180, 100, 0, 180);
    ThetanThetaRecPi->Merge(list37);

    TList *list38 = new TList;
    list38->Add(ThetanThetaRecPiDiff_Para);
    ThetanThetaRecPiDiff_Perp->Scale(ScaleFactor);
    list38->Add(ThetanThetaRecPiDiff_Perp);
    ThetanThetaRecPiDiff = new TH2D("ThetanThetaRecPiDiff", "ThetanThetaRecPiDiff", 100, 0, 180, 100, 0, 180);
    ThetanThetaRecPiDiff->Merge(list38);

    TList *list39 = new TList;
    list39->Add(ThetaRecPDiff_Para);
    ThetaRecPDiff_Perp->Scale(ScaleFactor);
    list39->Add(ThetaRecPDiff_Perp);
    ThetaRecPDiff = new TH1D("ThetaRecPDiff", "ThetaRecPDiff", 200, 0, 180);
    ThetaRecPDiff->Merge(list39);

    TList *list40 = new TList;
    list40->Add(ThetanThetaRecP_Para);
    ThetanThetaRecP_Perp->Scale(ScaleFactor);
    list40->Add(ThetanThetaRecP_Perp);
    ThetanThetaRecP = new TH2D("ThetanThetaRecP", "ThetanThetaRecP", 100, 0, 180, 100, 0, 180);
    ThetanThetaRecP->Merge(list40);

    TList *list41 = new TList;
    list41->Add(ThetanThetaRecPDiff_Para);
    ThetanThetaRecPDiff_Perp->Scale(ScaleFactor);
    list41->Add(ThetanThetaRecPDiff_Perp);
    ThetanThetaRecPDiff = new TH2D("ThetanThetaRecPDiff", "ThetanThetaRecPDiff", 100, 0, 180, 100, 0, 180);
    ThetanThetaRecPDiff->Merge(list41);

    TList *list42 = new TList;
    list42->Add(E_dE_Para);
    E_dE_Perp->Scale(ScaleFactor);
    list42->Add(E_dE_Perp);
    E_dE = new TH2D("E_dE", "E_dE", 100, 0, 500, 100, 0, 5);
    E_dE->Merge(list42);

    TList *list43 = new TList;
    list43->Add(KinEp_dE_Para);
    KinEp_dE_Perp->Scale(ScaleFactor);
    list43->Add(KinEp_dE_Perp);
    KinEp_dE = new TH2D("KinEp_dE", "KinEp_dE", 100, 0, 500, 100, 0, 5);
    KinEp_dE->Merge(list43);

    TList *list44 = new TList;
    list44->Add(ThetaScPhiSc_Para);
    ThetaScPhiSc_Perp->Scale(ScaleFactor);
    list44->Add(ThetaScPhiSc_Perp);
    ThetaScPhiSc = new TH2D("ThetaScPhiSc", "ThetaScPhiSc", 100, 0, 180, 100, -180, 180);
    ThetaScPhiSc->Merge(list44);

    TList *list45 = new TList;
    list45->Add(E_KinEp_Para);
    E_KinEp_Perp->Scale(ScaleFactor);
    list45->Add(E_KinEp_Perp);
    E_KinEp = new TH2D("E_KinEp", "E_KinEp", 100, 0, 500, 100, 0, 500);
    E_KinEp->Merge(list45);

    TList *list46 = new TList;
    list46->Add(PhinDiffWCZRec_Para);
    PhinDiffWCZRec_Perp->Scale(ScaleFactor);
    list46->Add(PhinDiffWCZRec_Perp);
    PhinDiffWCZRec = new TH2D("PhinDiffWCZRec", "PhinDiffWCZRec", 100, 0, 200, 100, 0, 180);
    PhinDiffWCZRec->Merge(list46);

    TFile f2("ParaPerp_Total_12_Combined_Unpolarised.root", "RECREATE");

    time->Write();
    time_cut->Write();
    Eg->Write();
    WCPhiDifference->Write();
    EpKin->Write();
    EpCorrected->Write();
    OAngle->Write();
    WCZnRecon->Write();
    ThetaSc->Write();
    PhiSc->Write();
    EpKinEpCorrDiff->Write();
    EpEpCorrDiff->Write();
    MMpEpCorrected->Write();
    ZpDist->Write();
    ZpPhiScatNeg180->Write();
    ZpPhiScat0->Write();
    ZpPhiScatPos180->Write();
    MMp200300->Write();
    MMp300400->Write();
    MMp400500->Write();
    MMp500600->Write();
    MMp600700->Write();
    MMp700800->Write();
    MMp800900->Write();
    PhiSc410->Write();
    PhiSc430->Write();
    PhiSc450->Write();
    PhiSc470->Write();
    PhiSc490->Write();
    PhiSc510->Write();
    PhiSc530->Write();
    PhiSc550->Write();
    PhiSc570->Write();
    PhiSc590->Write();
    PhiSc610->Write();
    PhiSc630->Write();
    ThetaRecPiDiff->Write();
    ThetanThetaRecPi->Write();
    ThetanThetaRecPiDiff->Write();
    ThetaRecPDiff->Write();
    ThetanThetaRecP->Write();
    ThetanThetaRecPDiff->Write();
    E_dE->Write();
    KinEp_dE->Write();
    ThetaScPhiSc->Write();
    E_KinEp->Write();
    PhinDiffWCZRec->Write();

    f2.Write();

}
