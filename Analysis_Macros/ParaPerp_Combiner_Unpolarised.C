#include "./includes_ParaPerp_Combiner_Unpolarised.h"

void ParaPerp_Combiner_Unpolarised() {

    TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/Para/Physics_Total_Para_10_28_4_17.root"); // Open latest Para file

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

    ////////////////////////////////////////////////////
    ////////////////////////////////////////////////////
    ///////////////// PARA DONE ////////////////////////
    ////////////////////////////////////////////////////
    ////////////////////////////////////////////////////

    TFile *f1 = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/Perp/Physics_Total_Perp_10_28_4_17.root"); // Open latest Para file

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
    time_Unpol = new TH1D("time_Unpol", "time_Unpol", 2, -1, 1);
    time_Unpol->Merge(list1);

    TFile f2("ParaPerp_Total_10_Combined_Unpolarised.root", "RECREATE");

    time_Unpol->Write();

    f2.Write();

}
