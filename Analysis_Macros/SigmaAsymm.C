#include "./includes_SigmaAsymm.h"

void SigmaAsymm(){

    TF1 *CosFunc = new TF1("CosFit", "[0]*cos((2*x*TMath::DegToRad())+acos(0))");
    CosFunc->SetParNames("Amplitude");

    double pCosAmp[10][10]; // Format of array is Theta bin (x) by Egamma bin (y), 10 CosTheta bins, 10 20MeV Egamma bins
    double pCosAmpErr[10][10];
    double pCosA435;
    double pCosAErr435;
    double pCosA455;
    double pCosAErr455;
    double pCosA475;
    double pCosAErr475;
    double pCosA495;
    double pCosAErr495;
    double pCosA515;
    double pCosAErr515;
    double pCosA535;
    double pCosAErr535;
    double pCosA555;
    double pCosAErr555;
    double pCosA575;
    double pCosAErr575;
    double pCosA595;
    double pCosAErr595;
    double pCosA615;
    double pCosAErr615;

    TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/ParaPerp_Total_17_Combined.root"); // Open the latest PTotal combined file to load histograms from
    NPara = Eg_Para->GetEntries();
    NPerp = Eg_Perp->GetEntries();
    ScaleFactor = NPara/NPerp;
    ScaleFactorErr = sqrt( (NPara/((TMath::Power(NPerp,2)))) + (((TMath::Power(NPara,2)))/(TMath::Power(NPerp,3))) ); // Error Propagation of above formula, see notebook 5

    ///////////////////////////////////////////
    ////////////////  CM1  ////////////////////
    ///////////////////////////////////////////

    ParaPerpAsymmPhip_435MeVCM1 = Phip_435MeVCM1_Para->GetAsymmetry(Phip_435MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_435MeVCM1->SetName("ParaPerpAsymmPhip435MeVCM1");
    ParaPerpAsymmPhip_435MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 435 #pm 10MeV (CosTheta1-0.8)");
    ParaPerpAsymmPhip_435MeVCM1->Fit("CosFit", "Q");
    pCosAmp[0][0] = CosFit->GetParameter(0);
    pCosAmpErr[0][0] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_455MeVCM1 = Phip_455MeVCM1_Para->GetAsymmetry(Phip_455MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_455MeVCM1->SetName("ParaPerpAsymmPhip455MeVCM1");
    ParaPerpAsymmPhip_455MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 455 #pm 10MeV (CosTheta1-0.8)");
    ParaPerpAsymmPhip_455MeVCM1->Fit("CosFit", "Q");
    pCosAmp[0][1] = CosFit->GetParameter(0);
    pCosAmpErr[0][1] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_475MeVCM1 = Phip_475MeVCM1_Para->GetAsymmetry(Phip_475MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_475MeVCM1->SetName("ParaPerpAsymmPhip475MeVCM1");
    ParaPerpAsymmPhip_475MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 475 #pm 10MeV (CosTheta1-0.8)");
    ParaPerpAsymmPhip_475MeVCM1->Fit("CosFit", "Q");
    pCosAmp[0][2] = CosFit->GetParameter(0);
    pCosAmpErr[0][2] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_495MeVCM1 = Phip_495MeVCM1_Para->GetAsymmetry(Phip_495MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_495MeVCM1->SetName("ParaPerpAsymmPhip495MeVCM1");
    ParaPerpAsymmPhip_495MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 495 #pm 10MeV (CosTheta1-0.8)");
    ParaPerpAsymmPhip_495MeVCM1->Fit("CosFit", "Q");
    pCosAmp[0][3] = CosFit->GetParameter(0);
    pCosAmpErr[0][3] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_515MeVCM1 = Phip_515MeVCM1_Para->GetAsymmetry(Phip_515MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_515MeVCM1->SetName("ParaPerpAsymmPhip515MeVCM1");
    ParaPerpAsymmPhip_515MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 515 #pm 10MeV (CosTheta1-0.8)");
    ParaPerpAsymmPhip_515MeVCM1->Fit("CosFit", "Q");
    pCosAmp[0][4] = CosFit->GetParameter(0);
    pCosAmpErr[0][4] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_535MeVCM1 = Phip_535MeVCM1_Para->GetAsymmetry(Phip_535MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_535MeVCM1->SetName("ParaPerpAsymmPhip535MeVCM1");
    ParaPerpAsymmPhip_535MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 535 #pm 10MeV (CosTheta1-0.8)");
    ParaPerpAsymmPhip_535MeVCM1->Fit("CosFit", "Q");
    pCosAmp[0][5] = CosFit->GetParameter(0);
    pCosAmpErr[0][5] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_555MeVCM1 = Phip_555MeVCM1_Para->GetAsymmetry(Phip_555MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_555MeVCM1->SetName("ParaPerpAsymmPhip555MeVCM1");
    ParaPerpAsymmPhip_555MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 555 #pm 10MeV (CosTheta1-0.8)");
    ParaPerpAsymmPhip_555MeVCM1->Fit("CosFit", "Q");
    pCosAmp[0][6] = CosFit->GetParameter(0);
    pCosAmpErr[0][6] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_575MeVCM1 = Phip_575MeVCM1_Para->GetAsymmetry(Phip_575MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_575MeVCM1->SetName("ParaPerpAsymmPhip575MeVCM1");
    ParaPerpAsymmPhip_575MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 575 #pm 10MeV (CosTheta1-0.8)");
    ParaPerpAsymmPhip_575MeVCM1->Fit("CosFit", "Q");
    pCosAmp[0][7] = CosFit->GetParameter(0);
    pCosAmpErr[0][7] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_595MeVCM1 = Phip_595MeVCM1_Para->GetAsymmetry(Phip_595MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_595MeVCM1->SetName("ParaPerpAsymmPhip595MeVCM1");
    ParaPerpAsymmPhip_595MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 595 #pm 10MeV (CosTheta1-0.8)");
    ParaPerpAsymmPhip_595MeVCM1->Fit("CosFit", "Q");
    pCosAmp[0][8] = CosFit->GetParameter(0);
    pCosAmpErr[0][8] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_615MeVCM1 = Phip_615MeVCM1_Para->GetAsymmetry(Phip_615MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_615MeVCM1->SetName("ParaPerpAsymmPhip615MeVCM1");
    ParaPerpAsymmPhip_615MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 615 #pm 10MeV (CosTheta1-0.8)");
    ParaPerpAsymmPhip_615MeVCM1->Fit("CosFit", "Q");
    pCosAmp[0][9] = CosFit->GetParameter(0);
    pCosAmpErr[0][9] = CosFit->GetParError(0);

    ///////////////////////////////////////////
    ////////////////  CM2  ////////////////////
    ///////////////////////////////////////////

    ParaPerpAsymmPhip_435MeVCM2 = Phip_435MeVCM2_Para->GetAsymmetry(Phip_435MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_435MeVCM2->SetName("ParaPerpAsymmPhip435MeVCM2");
    ParaPerpAsymmPhip_435MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 435 #pm 10MeV (CosTheta0.8-0.6)");
    ParaPerpAsymmPhip_435MeVCM2->Fit("CosFit", "Q");
    pCosAmp[1][0] = CosFit->GetParameter(0);
    pCosAmpErr[1][0] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_455MeVCM2 = Phip_455MeVCM2_Para->GetAsymmetry(Phip_455MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_455MeVCM2->SetName("ParaPerpAsymmPhip455MeVCM2");
    ParaPerpAsymmPhip_455MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 455 #pm 10MeV (CosTheta0.8-0.6)");
    ParaPerpAsymmPhip_455MeVCM2->Fit("CosFit", "Q");
    pCosAmp[1][1] = CosFit->GetParameter(0);
    pCosAmpErr[1][1] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_475MeVCM2 = Phip_475MeVCM2_Para->GetAsymmetry(Phip_475MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_475MeVCM2->SetName("ParaPerpAsymmPhip475MeVCM2");
    ParaPerpAsymmPhip_475MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 475 #pm 10MeV (CosTheta0.8-0.6)");
    ParaPerpAsymmPhip_475MeVCM2->Fit("CosFit", "Q");
    pCosAmp[1][2] = CosFit->GetParameter(0);
    pCosAmpErr[1][2] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_495MeVCM2 = Phip_495MeVCM2_Para->GetAsymmetry(Phip_495MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_495MeVCM2->SetName("ParaPerpAsymmPhip495MeVCM2");
    ParaPerpAsymmPhip_495MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 495 #pm 10MeV (CosTheta0.8-0.6)");
    ParaPerpAsymmPhip_495MeVCM2->Fit("CosFit", "Q");
    pCosAmp[1][3] = CosFit->GetParameter(0);
    pCosAmpErr[1][3] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_515MeVCM2 = Phip_515MeVCM2_Para->GetAsymmetry(Phip_515MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_515MeVCM2->SetName("ParaPerpAsymmPhip515MeVCM2");
    ParaPerpAsymmPhip_515MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 515 #pm 10MeV (CosTheta0.8-0.6)");
    ParaPerpAsymmPhip_515MeVCM2->Fit("CosFit", "Q");
    pCosAmp[1][4] = CosFit->GetParameter(0);
    pCosAmpErr[1][4] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_535MeVCM2 = Phip_535MeVCM2_Para->GetAsymmetry(Phip_535MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_535MeVCM2->SetName("ParaPerpAsymmPhip535MeVCM2");
    ParaPerpAsymmPhip_535MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 535 #pm 10MeV (CosTheta0.8-0.6)");
    ParaPerpAsymmPhip_535MeVCM2->Fit("CosFit", "Q");
    pCosAmp[1][5] = CosFit->GetParameter(0);
    pCosAmpErr[1][5] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_555MeVCM2 = Phip_555MeVCM2_Para->GetAsymmetry(Phip_555MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_555MeVCM2->SetName("ParaPerpAsymmPhip555MeVCM2");
    ParaPerpAsymmPhip_555MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 555 #pm 10MeV (CosTheta0.8-0.6)");
    ParaPerpAsymmPhip_555MeVCM2->Fit("CosFit", "Q");
    pCosAmp[1][6] = CosFit->GetParameter(0);
    pCosAmpErr[1][6] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_575MeVCM2 = Phip_575MeVCM2_Para->GetAsymmetry(Phip_575MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_575MeVCM2->SetName("ParaPerpAsymmPhip575MeVCM2");
    ParaPerpAsymmPhip_575MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 575 #pm 10MeV (CosTheta0.8-0.6)");
    ParaPerpAsymmPhip_575MeVCM2->Fit("CosFit", "Q");
    pCosAmp[1][7] = CosFit->GetParameter(0);
    pCosAmpErr[1][7] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_595MeVCM2 = Phip_595MeVCM2_Para->GetAsymmetry(Phip_595MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_595MeVCM2->SetName("ParaPerpAsymmPhip595MeVCM2");
    ParaPerpAsymmPhip_595MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 595 #pm 10MeV (CosTheta0.8-0.6)");
    ParaPerpAsymmPhip_595MeVCM2->Fit("CosFit", "Q");
    pCosAmp[1][8] = CosFit->GetParameter(0);
    pCosAmpErr[1][8] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_615MeVCM2 = Phip_615MeVCM2_Para->GetAsymmetry(Phip_615MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_615MeVCM2->SetName("ParaPerpAsymmPhip615MeVCM2");
    ParaPerpAsymmPhip_615MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 615 #pm 10MeV (CosTheta0.8-0.6)");
    ParaPerpAsymmPhip_615MeVCM2->Fit("CosFit", "Q");
    pCosAmp[1][9] = CosFit->GetParameter(0);
    pCosAmpErr[1][9] = CosFit->GetParError(0);

    ///////////////////////////////////////////
    ////////////////  CM3  ////////////////////
    ///////////////////////////////////////////

    ParaPerpAsymmPhip_435MeVCM3 = Phip_435MeVCM3_Para->GetAsymmetry(Phip_435MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_435MeVCM3->SetName("ParaPerpAsymmPhip435MeVCM3");
    ParaPerpAsymmPhip_435MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 435 #pm 10MeV (CosTheta0.6-0.4)");
    ParaPerpAsymmPhip_435MeVCM3->Fit("CosFit", "Q");
    pCosAmp[2][0] = CosFit->GetParameter(0);
    pCosAmpErr[2][0] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_455MeVCM3 = Phip_455MeVCM3_Para->GetAsymmetry(Phip_455MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_455MeVCM3->SetName("ParaPerpAsymmPhip455MeVCM3");
    ParaPerpAsymmPhip_455MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 455 #pm 10MeV (CosTheta0.6-0.4)");
    ParaPerpAsymmPhip_455MeVCM3->Fit("CosFit", "Q");
    pCosAmp[2][1] = CosFit->GetParameter(0);
    pCosAmpErr[2][1] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_475MeVCM3 = Phip_475MeVCM3_Para->GetAsymmetry(Phip_475MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_475MeVCM3->SetName("ParaPerpAsymmPhip475MeVCM3");
    ParaPerpAsymmPhip_475MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 475 #pm 10MeV (CosTheta0.6-0.4)");
    ParaPerpAsymmPhip_475MeVCM3->Fit("CosFit", "Q");
    pCosAmp[2][2] = CosFit->GetParameter(0);
    pCosAmpErr[2][2] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_495MeVCM3 = Phip_495MeVCM3_Para->GetAsymmetry(Phip_495MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_495MeVCM3->SetName("ParaPerpAsymmPhip495MeVCM3");
    ParaPerpAsymmPhip_495MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 495 #pm 10MeV (CosTheta0.6-0.4)");
    ParaPerpAsymmPhip_495MeVCM3->Fit("CosFit", "Q");
    pCosAmp[2][3] = CosFit->GetParameter(0);
    pCosAmpErr[2][3] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_515MeVCM3 = Phip_515MeVCM3_Para->GetAsymmetry(Phip_515MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_515MeVCM3->SetName("ParaPerpAsymmPhip515MeVCM3");
    ParaPerpAsymmPhip_515MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 515 #pm 10MeV (CosTheta0.6-0.4)");
    ParaPerpAsymmPhip_515MeVCM3->Fit("CosFit", "Q");
    pCosAmp[2][4] = CosFit->GetParameter(0);
    pCosAmpErr[2][4] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_535MeVCM3 = Phip_535MeVCM3_Para->GetAsymmetry(Phip_535MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_535MeVCM3->SetName("ParaPerpAsymmPhip535MeVCM3");
    ParaPerpAsymmPhip_535MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 535 #pm 10MeV (CosTheta0.6-0.4)");
    ParaPerpAsymmPhip_535MeVCM3->Fit("CosFit", "Q");
    pCosAmp[2][5] = CosFit->GetParameter(0);
    pCosAmpErr[2][5] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_555MeVCM3 = Phip_555MeVCM3_Para->GetAsymmetry(Phip_555MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_555MeVCM3->SetName("ParaPerpAsymmPhip555MeVCM3");
    ParaPerpAsymmPhip_555MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 555 #pm 10MeV (CosTheta0.6-0.4)");
    ParaPerpAsymmPhip_555MeVCM3->Fit("CosFit", "Q");
    pCosAmp[2][6] = CosFit->GetParameter(0);
    pCosAmpErr[2][6] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_575MeVCM3 = Phip_575MeVCM3_Para->GetAsymmetry(Phip_575MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_575MeVCM3->SetName("ParaPerpAsymmPhip575MeVCM3");
    ParaPerpAsymmPhip_575MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 575 #pm 10MeV (CosTheta0.6-0.4)");
    ParaPerpAsymmPhip_575MeVCM3->Fit("CosFit", "Q");
    pCosAmp[2][7] = CosFit->GetParameter(0);
    pCosAmpErr[2][7] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_595MeVCM3 = Phip_595MeVCM3_Para->GetAsymmetry(Phip_595MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_595MeVCM3->SetName("ParaPerpAsymmPhip595MeVCM3");
    ParaPerpAsymmPhip_595MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 595 #pm 10MeV (CosTheta0.6-0.4)");
    ParaPerpAsymmPhip_595MeVCM3->Fit("CosFit", "Q");
    pCosAmp[2][8] = CosFit->GetParameter(0);
    pCosAmpErr[2][8] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_615MeVCM3 = Phip_615MeVCM3_Para->GetAsymmetry(Phip_615MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_615MeVCM3->SetName("ParaPerpAsymmPhip615MeVCM3");
    ParaPerpAsymmPhip_615MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 615 #pm 10MeV (CosTheta0.6-0.4)");
    ParaPerpAsymmPhip_615MeVCM3->Fit("CosFit", "Q");
    pCosAmp[2][9] = CosFit->GetParameter(0);
    pCosAmpErr[2][9] = CosFit->GetParError(0);

    ///////////////////////////////////////////
    ////////////////  CM4  ////////////////////
    ///////////////////////////////////////////

    ParaPerpAsymmPhip_435MeVCM4 = Phip_435MeVCM4_Para->GetAsymmetry(Phip_435MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_435MeVCM4->SetName("ParaPerpAsymmPhip435MeVCM4");
    ParaPerpAsymmPhip_435MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 435 #pm 10MeV (CosTheta0.4-0.2)");
    ParaPerpAsymmPhip_435MeVCM4->Fit("CosFit", "Q");
    pCosAmp[3][0] = CosFit->GetParameter(0);
    pCosAmpErr[3][0] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_455MeVCM4 = Phip_455MeVCM4_Para->GetAsymmetry(Phip_455MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_455MeVCM4->SetName("ParaPerpAsymmPhip455MeVCM4");
    ParaPerpAsymmPhip_455MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 455 #pm 10MeV (CosTheta0.4-0.2)");
    ParaPerpAsymmPhip_455MeVCM4->Fit("CosFit", "Q");
    pCosAmp[3][1] = CosFit->GetParameter(0);
    pCosAmpErr[3][1] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_475MeVCM4 = Phip_475MeVCM4_Para->GetAsymmetry(Phip_475MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_475MeVCM4->SetName("ParaPerpAsymmPhip475MeVCM4");
    ParaPerpAsymmPhip_475MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 475 #pm 10MeV (CosTheta0.4-0.2)");
    ParaPerpAsymmPhip_475MeVCM4->Fit("CosFit", "Q");
    pCosAmp[3][2] = CosFit->GetParameter(0);
    pCosAmpErr[3][2] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_495MeVCM4 = Phip_495MeVCM4_Para->GetAsymmetry(Phip_495MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_495MeVCM4->SetName("ParaPerpAsymmPhip495MeVCM4");
    ParaPerpAsymmPhip_495MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 495 #pm 10MeV (CosTheta0.4-0.2)");
    ParaPerpAsymmPhip_495MeVCM4->Fit("CosFit", "Q");
    pCosAmp[3][3] = CosFit->GetParameter(0);
    pCosAmpErr[3][3] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_515MeVCM4 = Phip_515MeVCM4_Para->GetAsymmetry(Phip_515MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_515MeVCM4->SetName("ParaPerpAsymmPhip515MeVCM4");
    ParaPerpAsymmPhip_515MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 515 #pm 10MeV (CosTheta0.4-0.2)");
    ParaPerpAsymmPhip_515MeVCM4->Fit("CosFit", "Q");
    pCosAmp[3][4] = CosFit->GetParameter(0);
    pCosAmpErr[3][4] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_535MeVCM4 = Phip_535MeVCM4_Para->GetAsymmetry(Phip_535MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_535MeVCM4->SetName("ParaPerpAsymmPhip535MeVCM4");
    ParaPerpAsymmPhip_535MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 535 #pm 10MeV (CosTheta0.4-0.2)");
    ParaPerpAsymmPhip_535MeVCM4->Fit("CosFit", "Q");
    pCosAmp[3][5] = CosFit->GetParameter(0);
    pCosAmpErr[3][5] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_555MeVCM4 = Phip_555MeVCM4_Para->GetAsymmetry(Phip_555MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_555MeVCM4->SetName("ParaPerpAsymmPhip555MeVCM4");
    ParaPerpAsymmPhip_555MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 555 #pm 10MeV (CosTheta0.4-0.2)");
    ParaPerpAsymmPhip_555MeVCM4->Fit("CosFit", "Q");
    pCosAmp[3][6] = CosFit->GetParameter(0);
    pCosAmpErr[3][6] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_575MeVCM4 = Phip_575MeVCM4_Para->GetAsymmetry(Phip_575MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_575MeVCM4->SetName("ParaPerpAsymmPhip575MeVCM4");
    ParaPerpAsymmPhip_575MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 575 #pm 10MeV (CosTheta0.4-0.2)");
    ParaPerpAsymmPhip_575MeVCM4->Fit("CosFit", "Q");
    pCosAmp[3][7] = CosFit->GetParameter(0);
    pCosAmpErr[3][7] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_595MeVCM4 = Phip_595MeVCM4_Para->GetAsymmetry(Phip_595MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_595MeVCM4->SetName("ParaPerpAsymmPhip595MeVCM4");
    ParaPerpAsymmPhip_595MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 595 #pm 10MeV (CosTheta0.4-0.2)");
    ParaPerpAsymmPhip_595MeVCM4->Fit("CosFit", "Q");
    pCosAmp[3][8] = CosFit->GetParameter(0);
    pCosAmpErr[3][8] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_615MeVCM4 = Phip_615MeVCM4_Para->GetAsymmetry(Phip_615MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_615MeVCM4->SetName("ParaPerpAsymmPhip615MeVCM4");
    ParaPerpAsymmPhip_615MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 615 #pm 10MeV (CosTheta0.4-0.2)");
    ParaPerpAsymmPhip_615MeVCM4->Fit("CosFit", "Q");
    pCosAmp[3][9] = CosFit->GetParameter(0);
    pCosAmpErr[3][9] = CosFit->GetParError(0);

    ///////////////////////////////////////////
    ////////////////  CM5  ////////////////////
    ///////////////////////////////////////////

    ParaPerpAsymmPhip_435MeVCM5 = Phip_435MeVCM5_Para->GetAsymmetry(Phip_435MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_435MeVCM5->SetName("ParaPerpAsymmPhip435MeVCM5");
    ParaPerpAsymmPhip_435MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 435 #pm 10MeV (CosTheta0.2-0.0)");
    ParaPerpAsymmPhip_435MeVCM5->Fit("CosFit", "Q");
    pCosAmp[4][0] = CosFit->GetParameter(0);
    pCosAmpErr[4][0] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_455MeVCM5 = Phip_455MeVCM5_Para->GetAsymmetry(Phip_455MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_455MeVCM5->SetName("ParaPerpAsymmPhip455MeVCM5");
    ParaPerpAsymmPhip_455MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 455 #pm 10MeV (CosTheta0.2-0.0)");
    ParaPerpAsymmPhip_455MeVCM5->Fit("CosFit", "Q");
    pCosAmp[4][1] = CosFit->GetParameter(0);
    pCosAmpErr[4][1] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_475MeVCM5 = Phip_475MeVCM5_Para->GetAsymmetry(Phip_475MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_475MeVCM5->SetName("ParaPerpAsymmPhip475MeVCM5");
    ParaPerpAsymmPhip_475MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 475 #pm 10MeV (CosTheta0.2-0.0)");
    ParaPerpAsymmPhip_475MeVCM5->Fit("CosFit", "Q");
    pCosAmp[4][2] = CosFit->GetParameter(0);
    pCosAmpErr[4][2] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_495MeVCM5 = Phip_495MeVCM5_Para->GetAsymmetry(Phip_495MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_495MeVCM5->SetName("ParaPerpAsymmPhip495MeVCM5");
    ParaPerpAsymmPhip_495MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 495 #pm 10MeV (CosTheta0.2-0.0)");
    ParaPerpAsymmPhip_495MeVCM5->Fit("CosFit", "Q");
    pCosAmp[4][3] = CosFit->GetParameter(0);
    pCosAmpErr[4][3] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_515MeVCM5 = Phip_515MeVCM5_Para->GetAsymmetry(Phip_515MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_515MeVCM5->SetName("ParaPerpAsymmPhip515MeVCM5");
    ParaPerpAsymmPhip_515MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 515 #pm 10MeV (CosTheta0.2-0.0)");
    ParaPerpAsymmPhip_515MeVCM5->Fit("CosFit", "Q");
    pCosAmp[4][4] = CosFit->GetParameter(0);
    pCosAmpErr[4][4] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_535MeVCM5 = Phip_535MeVCM5_Para->GetAsymmetry(Phip_535MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_535MeVCM5->SetName("ParaPerpAsymmPhip535MeVCM5");
    ParaPerpAsymmPhip_535MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 535 #pm 10MeV (CosTheta0.2-0.0)");
    ParaPerpAsymmPhip_535MeVCM5->Fit("CosFit", "Q");
    pCosAmp[4][5] = CosFit->GetParameter(0);
    pCosAmpErr[4][5] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_555MeVCM5 = Phip_555MeVCM5_Para->GetAsymmetry(Phip_555MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_555MeVCM5->SetName("ParaPerpAsymmPhip555MeVCM5");
    ParaPerpAsymmPhip_555MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 555 #pm 10MeV (CosTheta0.2-0.0)");
    ParaPerpAsymmPhip_555MeVCM5->Fit("CosFit", "Q");
    pCosAmp[4][6] = CosFit->GetParameter(0);
    pCosAmpErr[4][6] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_575MeVCM5 = Phip_575MeVCM5_Para->GetAsymmetry(Phip_575MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_575MeVCM5->SetName("ParaPerpAsymmPhip575MeVCM5");
    ParaPerpAsymmPhip_575MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 575 #pm 10MeV (CosTheta0.2-0.0)");
    ParaPerpAsymmPhip_575MeVCM5->Fit("CosFit", "Q");
    pCosAmp[4][7] = CosFit->GetParameter(0);
    pCosAmpErr[4][7] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_595MeVCM5 = Phip_595MeVCM5_Para->GetAsymmetry(Phip_595MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_595MeVCM5->SetName("ParaPerpAsymmPhip595MeVCM5");
    ParaPerpAsymmPhip_595MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 595 #pm 10MeV (CosTheta0.2-0.0)");
    ParaPerpAsymmPhip_595MeVCM5->Fit("CosFit", "Q");
    pCosAmp[4][8] = CosFit->GetParameter(0);
    pCosAmpErr[4][8] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_615MeVCM5 = Phip_615MeVCM5_Para->GetAsymmetry(Phip_615MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_615MeVCM5->SetName("ParaPerpAsymmPhip615MeVCM5");
    ParaPerpAsymmPhip_615MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 615 #pm 10MeV (CosTheta0.2-0.0)");
    ParaPerpAsymmPhip_615MeVCM5->Fit("CosFit", "Q");
    pCosAmp[4][9] = CosFit->GetParameter(0);
    pCosAmpErr[4][9] = CosFit->GetParError(0);

    ///////////////////////////////////////////
    ////////////////  CM6  ////////////////////
    ///////////////////////////////////////////

    ParaPerpAsymmPhip_435MeVCM6 = Phip_435MeVCM6_Para->GetAsymmetry(Phip_435MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_435MeVCM6->SetName("ParaPerpAsymmPhip435MeVCM6");
    ParaPerpAsymmPhip_435MeVCM6->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 435 #pm 10MeV (CosTheta0.0-(-0.2))");
    ParaPerpAsymmPhip_435MeVCM6->Fit("CosFit", "Q");
    pCosAmp[5][0] = CosFit->GetParameter(0);
    pCosAmpErr[5][0] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_455MeVCM6 = Phip_455MeVCM6_Para->GetAsymmetry(Phip_455MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_455MeVCM6->SetName("ParaPerpAsymmPhip455MeVCM6");
    ParaPerpAsymmPhip_455MeVCM6->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 455 #pm 10MeV (CosTheta0.0-(-0.2))");
    ParaPerpAsymmPhip_455MeVCM6->Fit("CosFit", "Q");
    pCosAmp[5][1] = CosFit->GetParameter(0);
    pCosAmpErr[5][1] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_475MeVCM6 = Phip_475MeVCM6_Para->GetAsymmetry(Phip_475MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_475MeVCM6->SetName("ParaPerpAsymmPhip475MeVCM6");
    ParaPerpAsymmPhip_475MeVCM6->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 475 #pm 10MeV (CosTheta0.0-(-0.2))");
    ParaPerpAsymmPhip_475MeVCM6->Fit("CosFit", "Q");
    pCosAmp[5][2] = CosFit->GetParameter(0);
    pCosAmpErr[5][2] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_495MeVCM6 = Phip_495MeVCM6_Para->GetAsymmetry(Phip_495MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_495MeVCM6->SetName("ParaPerpAsymmPhip495MeVCM6");
    ParaPerpAsymmPhip_495MeVCM6->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 495 #pm 10MeV (CosTheta0.0-(-0.2))");
    ParaPerpAsymmPhip_495MeVCM6->Fit("CosFit", "Q");
    pCosAmp[5][3] = CosFit->GetParameter(0);
    pCosAmpErr[5][3] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_515MeVCM6 = Phip_515MeVCM6_Para->GetAsymmetry(Phip_515MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_515MeVCM6->SetName("ParaPerpAsymmPhip515MeVCM6");
    ParaPerpAsymmPhip_515MeVCM6->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 515 #pm 10MeV (CosTheta0.0-(-0.2))");
    ParaPerpAsymmPhip_515MeVCM6->Fit("CosFit", "Q");
    pCosAmp[5][4] = CosFit->GetParameter(0);
    pCosAmpErr[5][4] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_535MeVCM6 = Phip_535MeVCM6_Para->GetAsymmetry(Phip_535MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_535MeVCM6->SetName("ParaPerpAsymmPhip535MeVCM6");
    ParaPerpAsymmPhip_535MeVCM6->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 535 #pm 10MeV (CosTheta0.0-(-0.2))");
    ParaPerpAsymmPhip_535MeVCM6->Fit("CosFit", "Q");
    pCosAmp[5][5] = CosFit->GetParameter(0);
    pCosAmpErr[5][5] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_555MeVCM6 = Phip_555MeVCM6_Para->GetAsymmetry(Phip_555MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_555MeVCM6->SetName("ParaPerpAsymmPhip555MeVCM6");
    ParaPerpAsymmPhip_555MeVCM6->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 555 #pm 10MeV (CosTheta0.0-(-0.2))");
    ParaPerpAsymmPhip_555MeVCM6->Fit("CosFit", "Q");
    pCosAmp[5][6] = CosFit->GetParameter(0);
    pCosAmpErr[5][6] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_575MeVCM6 = Phip_575MeVCM6_Para->GetAsymmetry(Phip_575MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_575MeVCM6->SetName("ParaPerpAsymmPhip575MeVCM6");
    ParaPerpAsymmPhip_575MeVCM6->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 575 #pm 10MeV (CosTheta0.0-(-0.2))");
    ParaPerpAsymmPhip_575MeVCM6->Fit("CosFit", "Q");
    pCosAmp[5][7] = CosFit->GetParameter(0);
    pCosAmpErr[5][7] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_595MeVCM6 = Phip_595MeVCM6_Para->GetAsymmetry(Phip_595MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_595MeVCM6->SetName("ParaPerpAsymmPhip595MeVCM6");
    ParaPerpAsymmPhip_595MeVCM6->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 595 #pm 10MeV (CosTheta0.0-(-0.2))");
    ParaPerpAsymmPhip_595MeVCM6->Fit("CosFit", "Q");
    pCosAmp[5][8] = CosFit->GetParameter(0);
    pCosAmpErr[5][8] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_615MeVCM6 = Phip_615MeVCM6_Para->GetAsymmetry(Phip_615MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_615MeVCM6->SetName("ParaPerpAsymmPhip615MeVCM6");
    ParaPerpAsymmPhip_615MeVCM6->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 615 #pm 10MeV (CosTheta0.0-(-0.2))");
    ParaPerpAsymmPhip_615MeVCM6->Fit("CosFit", "Q");
    pCosAmp[5][9] = CosFit->GetParameter(0);
    pCosAmpErr[5][9] = CosFit->GetParError(0);

    ///////////////////////////////////////////
    ////////////////  CM7  ////////////////////
    ///////////////////////////////////////////

    ParaPerpAsymmPhip_435MeVCM7 = Phip_435MeVCM7_Para->GetAsymmetry(Phip_435MeVCM7_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_435MeVCM7->SetName("ParaPerpAsymmPhip435MeVCM7");
    ParaPerpAsymmPhip_435MeVCM7->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 435 #pm 10MeV (CosTheta-0.2-(-0.4))");
    ParaPerpAsymmPhip_435MeVCM7->Fit("CosFit", "Q");
    pCosAmp[6][0] = CosFit->GetParameter(0);
    pCosAmpErr[6][0] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_455MeVCM7 = Phip_455MeVCM7_Para->GetAsymmetry(Phip_455MeVCM7_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_455MeVCM7->SetName("ParaPerpAsymmPhip455MeVCM7");
    ParaPerpAsymmPhip_455MeVCM7->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 455 #pm 10MeV (CosTheta-0.2-(-0.4))");
    ParaPerpAsymmPhip_455MeVCM7->Fit("CosFit", "Q");
    pCosAmp[6][1] = CosFit->GetParameter(0);
    pCosAmpErr[6][1] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_475MeVCM7 = Phip_475MeVCM7_Para->GetAsymmetry(Phip_475MeVCM7_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_475MeVCM7->SetName("ParaPerpAsymmPhip475MeVCM7");
    ParaPerpAsymmPhip_475MeVCM7->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 475 #pm 10MeV (CosTheta-0.2-(-0.4))");
    ParaPerpAsymmPhip_475MeVCM7->Fit("CosFit", "Q");
    pCosAmp[6][2] = CosFit->GetParameter(0);
    pCosAmpErr[6][2] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_495MeVCM7 = Phip_495MeVCM7_Para->GetAsymmetry(Phip_495MeVCM7_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_495MeVCM7->SetName("ParaPerpAsymmPhip495MeVCM7");
    ParaPerpAsymmPhip_495MeVCM7->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 495 #pm 10MeV (CosTheta-0.2-(-0.4))");
    ParaPerpAsymmPhip_495MeVCM7->Fit("CosFit", "Q");
    pCosAmp[6][3] = CosFit->GetParameter(0);
    pCosAmpErr[6][3] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_515MeVCM7 = Phip_515MeVCM7_Para->GetAsymmetry(Phip_515MeVCM7_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_515MeVCM7->SetName("ParaPerpAsymmPhip515MeVCM7");
    ParaPerpAsymmPhip_515MeVCM7->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 515 #pm 10MeV (CosTheta-0.2-(-0.4))");
    ParaPerpAsymmPhip_515MeVCM7->Fit("CosFit", "Q");
    pCosAmp[6][4] = CosFit->GetParameter(0);
    pCosAmpErr[6][4] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_535MeVCM7 = Phip_535MeVCM7_Para->GetAsymmetry(Phip_535MeVCM7_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_535MeVCM7->SetName("ParaPerpAsymmPhip535MeVCM7");
    ParaPerpAsymmPhip_535MeVCM7->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 535 #pm 10MeV (CosTheta-0.2-(-0.4))");
    ParaPerpAsymmPhip_535MeVCM7->Fit("CosFit", "Q");
    pCosAmp[6][5] = CosFit->GetParameter(0);
    pCosAmpErr[6][5] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_555MeVCM7 = Phip_555MeVCM7_Para->GetAsymmetry(Phip_555MeVCM7_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_555MeVCM7->SetName("ParaPerpAsymmPhip555MeVCM7");
    ParaPerpAsymmPhip_555MeVCM7->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 555 #pm 10MeV (CosTheta-0.2-(-0.4))");
    ParaPerpAsymmPhip_555MeVCM7->Fit("CosFit", "Q");
    pCosAmp[6][6] = CosFit->GetParameter(0);
    pCosAmpErr[6][6] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_575MeVCM7 = Phip_575MeVCM7_Para->GetAsymmetry(Phip_575MeVCM7_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_575MeVCM7->SetName("ParaPerpAsymmPhip575MeVCM7");
    ParaPerpAsymmPhip_575MeVCM7->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 575 #pm 10MeV (CosTheta-0.2-(-0.4))");
    ParaPerpAsymmPhip_575MeVCM7->Fit("CosFit", "Q");
    pCosAmp[6][7] = CosFit->GetParameter(0);
    pCosAmpErr[6][7] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_595MeVCM7 = Phip_595MeVCM7_Para->GetAsymmetry(Phip_595MeVCM7_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_595MeVCM7->SetName("ParaPerpAsymmPhip595MeVCM7");
    ParaPerpAsymmPhip_595MeVCM7->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 595 #pm 10MeV (CosTheta-0.2-(-0.4))");
    ParaPerpAsymmPhip_595MeVCM7->Fit("CosFit", "Q");
    pCosAmp[6][8] = CosFit->GetParameter(0);
    pCosAmpErr[6][8] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_615MeVCM7 = Phip_615MeVCM7_Para->GetAsymmetry(Phip_615MeVCM7_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_615MeVCM7->SetName("ParaPerpAsymmPhip615MeVCM7");
    ParaPerpAsymmPhip_615MeVCM7->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 615 #pm 10MeV (CosTheta-0.2-(-0.4))");
    ParaPerpAsymmPhip_615MeVCM7->Fit("CosFit", "Q");
    pCosAmp[6][9] = CosFit->GetParameter(0);
    pCosAmpErr[6][9] = CosFit->GetParError(0);

    ///////////////////////////////////////////
    ////////////////  CM8  ////////////////////
    ///////////////////////////////////////////

    ParaPerpAsymmPhip_435MeVCM8 = Phip_435MeVCM8_Para->GetAsymmetry(Phip_435MeVCM8_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_435MeVCM8->SetName("ParaPerpAsymmPhip435MeVCM8");
    ParaPerpAsymmPhip_435MeVCM8->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 435 #pm 10MeV (CosTheta-0.4-(-0.6))");
    ParaPerpAsymmPhip_435MeVCM8->Fit("CosFit", "Q");
    pCosAmp[7][0] = CosFit->GetParameter(0);
    pCosAmpErr[7][0] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_455MeVCM8 = Phip_455MeVCM8_Para->GetAsymmetry(Phip_455MeVCM8_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_455MeVCM8->SetName("ParaPerpAsymmPhip455MeVCM8");
    ParaPerpAsymmPhip_455MeVCM8->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 455 #pm 10MeV (CosTheta-0.4-(-0.6))");
    ParaPerpAsymmPhip_455MeVCM8->Fit("CosFit", "Q");
    pCosAmp[7][1] = CosFit->GetParameter(0);
    pCosAmpErr[7][1] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_475MeVCM8 = Phip_475MeVCM8_Para->GetAsymmetry(Phip_475MeVCM8_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_475MeVCM8->SetName("ParaPerpAsymmPhip475MeVCM8");
    ParaPerpAsymmPhip_475MeVCM8->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 475 #pm 10MeV (CosTheta-0.4-(-0.6))");
    ParaPerpAsymmPhip_475MeVCM8->Fit("CosFit", "Q");
    pCosAmp[7][2] = CosFit->GetParameter(0);
    pCosAmpErr[7][2] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_495MeVCM8 = Phip_495MeVCM8_Para->GetAsymmetry(Phip_495MeVCM8_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_495MeVCM8->SetName("ParaPerpAsymmPhip495MeVCM8");
    ParaPerpAsymmPhip_495MeVCM8->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 495 #pm 10MeV (CosTheta-0.4-(-0.6))");
    ParaPerpAsymmPhip_495MeVCM8->Fit("CosFit", "Q");
    pCosAmp[7][3] = CosFit->GetParameter(0);
    pCosAmpErr[7][3] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_515MeVCM8 = Phip_515MeVCM8_Para->GetAsymmetry(Phip_515MeVCM8_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_515MeVCM8->SetName("ParaPerpAsymmPhip515MeVCM8");
    ParaPerpAsymmPhip_515MeVCM8->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 515 #pm 10MeV (CosTheta-0.4-(-0.6))");
    ParaPerpAsymmPhip_515MeVCM8->Fit("CosFit", "Q");
    pCosAmp[7][4] = CosFit->GetParameter(0);
    pCosAmpErr[7][4] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_535MeVCM8 = Phip_535MeVCM8_Para->GetAsymmetry(Phip_535MeVCM8_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_535MeVCM8->SetName("ParaPerpAsymmPhip535MeVCM8");
    ParaPerpAsymmPhip_535MeVCM8->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 535 #pm 10MeV (CosTheta-0.4-(-0.6))");
    ParaPerpAsymmPhip_535MeVCM8->Fit("CosFit", "Q");
    pCosAmp[7][5] = CosFit->GetParameter(0);
    pCosAmpErr[7][5] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_555MeVCM8 = Phip_555MeVCM8_Para->GetAsymmetry(Phip_555MeVCM8_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_555MeVCM8->SetName("ParaPerpAsymmPhip555MeVCM8");
    ParaPerpAsymmPhip_555MeVCM8->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 555 #pm 10MeV (CosTheta-0.4-(-0.6))");
    ParaPerpAsymmPhip_555MeVCM8->Fit("CosFit", "Q");
    pCosAmp[7][6] = CosFit->GetParameter(0);
    pCosAmpErr[7][6] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_575MeVCM8 = Phip_575MeVCM8_Para->GetAsymmetry(Phip_575MeVCM8_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_575MeVCM8->SetName("ParaPerpAsymmPhip575MeVCM8");
    ParaPerpAsymmPhip_575MeVCM8->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 575 #pm 10MeV (CosTheta-0.4-(-0.6))");
    ParaPerpAsymmPhip_575MeVCM8->Fit("CosFit", "Q");
    pCosAmp[7][7] = CosFit->GetParameter(0);
    pCosAmpErr[7][7] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_595MeVCM8 = Phip_595MeVCM8_Para->GetAsymmetry(Phip_595MeVCM8_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_595MeVCM8->SetName("ParaPerpAsymmPhip595MeVCM8");
    ParaPerpAsymmPhip_595MeVCM8->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 595 #pm 10MeV (CosTheta-0.4-(-0.6))");
    ParaPerpAsymmPhip_595MeVCM8->Fit("CosFit", "Q");
    pCosAmp[7][8] = CosFit->GetParameter(0);
    pCosAmpErr[7][8] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_615MeVCM8 = Phip_615MeVCM8_Para->GetAsymmetry(Phip_615MeVCM8_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_615MeVCM8->SetName("ParaPerpAsymmPhip615MeVCM8");
    ParaPerpAsymmPhip_615MeVCM8->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 615 #pm 10MeV (CosTheta-0.4-(-0.6))");
    ParaPerpAsymmPhip_615MeVCM8->Fit("CosFit", "Q");
    pCosAmp[7][9] = CosFit->GetParameter(0);
    pCosAmpErr[7][9] = CosFit->GetParError(0);

    ///////////////////////////////////////////
    ////////////////  CM9  ////////////////////
    ///////////////////////////////////////////

    ParaPerpAsymmPhip_435MeVCM9 = Phip_435MeVCM9_Para->GetAsymmetry(Phip_435MeVCM9_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_435MeVCM9->SetName("ParaPerpAsymmPhip435MeVCM9");
    ParaPerpAsymmPhip_435MeVCM9->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 435 #pm 10MeV (CosTheta-0.6-(-0.8))");
    ParaPerpAsymmPhip_435MeVCM9->Fit("CosFit", "Q");
    pCosAmp[8][0] = CosFit->GetParameter(0);
    pCosAmpErr[8][0] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_455MeVCM9 = Phip_455MeVCM9_Para->GetAsymmetry(Phip_455MeVCM9_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_455MeVCM9->SetName("ParaPerpAsymmPhip455MeVCM9");
    ParaPerpAsymmPhip_455MeVCM9->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 455 #pm 10MeV (CosTheta-0.6-(-0.8))");
    ParaPerpAsymmPhip_455MeVCM9->Fit("CosFit", "Q");
    pCosAmp[8][1] = CosFit->GetParameter(0);
    pCosAmpErr[8][1] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_475MeVCM9 = Phip_475MeVCM9_Para->GetAsymmetry(Phip_475MeVCM9_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_475MeVCM9->SetName("ParaPerpAsymmPhip475MeVCM9");
    ParaPerpAsymmPhip_475MeVCM9->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 475 #pm 10MeV (CosTheta-0.6-(-0.8))");
    ParaPerpAsymmPhip_475MeVCM9->Fit("CosFit", "Q");
    pCosAmp[8][2] = CosFit->GetParameter(0);
    pCosAmpErr[8][2] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_495MeVCM9 = Phip_495MeVCM9_Para->GetAsymmetry(Phip_495MeVCM9_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_495MeVCM9->SetName("ParaPerpAsymmPhip495MeVCM9");
    ParaPerpAsymmPhip_495MeVCM9->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 495 #pm 10MeV (CosTheta-0.6-(-0.8))");
    ParaPerpAsymmPhip_495MeVCM9->Fit("CosFit", "Q");
    pCosAmp[8][3] = CosFit->GetParameter(0);
    pCosAmpErr[8][3] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_515MeVCM9 = Phip_515MeVCM9_Para->GetAsymmetry(Phip_515MeVCM9_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_515MeVCM9->SetName("ParaPerpAsymmPhip515MeVCM9");
    ParaPerpAsymmPhip_515MeVCM9->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 515 #pm 10MeV (CosTheta-0.6-(-0.8))");
    ParaPerpAsymmPhip_515MeVCM9->Fit("CosFit", "Q");
    pCosAmp[8][4] = CosFit->GetParameter(0);
    pCosAmpErr[8][4] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_535MeVCM9 = Phip_535MeVCM9_Para->GetAsymmetry(Phip_535MeVCM9_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_535MeVCM9->SetName("ParaPerpAsymmPhip535MeVCM9");
    ParaPerpAsymmPhip_535MeVCM9->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 535 #pm 10MeV (CosTheta-0.6-(-0.8))");
    ParaPerpAsymmPhip_535MeVCM9->Fit("CosFit", "Q");
    pCosAmp[8][5] = CosFit->GetParameter(0);
    pCosAmpErr[8][5] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_555MeVCM9 = Phip_555MeVCM9_Para->GetAsymmetry(Phip_555MeVCM9_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_555MeVCM9->SetName("ParaPerpAsymmPhip555MeVCM9");
    ParaPerpAsymmPhip_555MeVCM9->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 555 #pm 10MeV (CosTheta-0.6-(-0.8))");
    ParaPerpAsymmPhip_555MeVCM9->Fit("CosFit", "Q");
    pCosAmp[8][6] = CosFit->GetParameter(0);
    pCosAmpErr[8][6] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_575MeVCM9 = Phip_575MeVCM9_Para->GetAsymmetry(Phip_575MeVCM9_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_575MeVCM9->SetName("ParaPerpAsymmPhip575MeVCM9");
    ParaPerpAsymmPhip_575MeVCM9->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 575 #pm 10MeV (CosTheta-0.6-(-0.8))");
    ParaPerpAsymmPhip_575MeVCM9->Fit("CosFit", "Q");
    pCosAmp[8][7] = CosFit->GetParameter(0);
    pCosAmpErr[8][7] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_595MeVCM9 = Phip_595MeVCM9_Para->GetAsymmetry(Phip_595MeVCM9_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_595MeVCM9->SetName("ParaPerpAsymmPhip595MeVCM9");
    ParaPerpAsymmPhip_595MeVCM9->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 595 #pm 10MeV (CosTheta-0.6-(-0.8))");
    ParaPerpAsymmPhip_595MeVCM9->Fit("CosFit", "Q");
    pCosAmp[8][8] = CosFit->GetParameter(0);
    pCosAmpErr[8][8] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_615MeVCM9 = Phip_615MeVCM9_Para->GetAsymmetry(Phip_615MeVCM9_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_615MeVCM9->SetName("ParaPerpAsymmPhip615MeVCM9");
    ParaPerpAsymmPhip_615MeVCM9->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 615 #pm 10MeV (CosTheta-0.6-(-0.8))");
    ParaPerpAsymmPhip_615MeVCM9->Fit("CosFit", "Q");
    pCosAmp[8][9] = CosFit->GetParameter(0);
    pCosAmpErr[8][9] = CosFit->GetParError(0);

    ///////////////////////////////////////////
    ////////////////  CM10  ///////////////////
    ///////////////////////////////////////////

    ParaPerpAsymmPhip_435MeVCM10 = Phip_435MeVCM10_Para->GetAsymmetry(Phip_435MeVCM10_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_435MeVCM10->SetName("ParaPerpAsymmPhip435MeVCM10");
    ParaPerpAsymmPhip_435MeVCM10->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 435 #pm 10MeV (CosTheta-0.8-(-1.0))");
    ParaPerpAsymmPhip_435MeVCM10->Fit("CosFit", "Q");
    pCosAmp[9][0] = CosFit->GetParameter(0);
    pCosAmpErr[9][0] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_455MeVCM10 = Phip_455MeVCM10_Para->GetAsymmetry(Phip_455MeVCM10_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_455MeVCM10->SetName("ParaPerpAsymmPhip455MeVCM10");
    ParaPerpAsymmPhip_455MeVCM10->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 455 #pm 10MeV (CosTheta-0.8-(-1.0))");
    ParaPerpAsymmPhip_455MeVCM10->Fit("CosFit", "Q");
    pCosAmp[9][1] = CosFit->GetParameter(0);
    pCosAmpErr[9][1] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_475MeVCM10 = Phip_475MeVCM10_Para->GetAsymmetry(Phip_475MeVCM10_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_475MeVCM10->SetName("ParaPerpAsymmPhip475MeVCM10");
    ParaPerpAsymmPhip_475MeVCM10->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 475 #pm 10MeV (CosTheta-0.8-(-1.0))");
    ParaPerpAsymmPhip_475MeVCM10->Fit("CosFit", "Q");
    pCosAmp[9][2] = CosFit->GetParameter(0);
    pCosAmpErr[9][2] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_495MeVCM10 = Phip_495MeVCM10_Para->GetAsymmetry(Phip_495MeVCM10_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_495MeVCM10->SetName("ParaPerpAsymmPhip495MeVCM10");
    ParaPerpAsymmPhip_495MeVCM10->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 495 #pm 10MeV (CosTheta-0.8-(-1.0))");
    ParaPerpAsymmPhip_495MeVCM10->Fit("CosFit", "Q");
    pCosAmp[9][3] = CosFit->GetParameter(0);
    pCosAmpErr[9][3] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_515MeVCM10 = Phip_515MeVCM10_Para->GetAsymmetry(Phip_515MeVCM10_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_515MeVCM10->SetName("ParaPerpAsymmPhip515MeVCM10");
    ParaPerpAsymmPhip_515MeVCM10->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 515 #pm 10MeV (CosTheta-0.8-(-1.0))");
    ParaPerpAsymmPhip_515MeVCM10->Fit("CosFit", "Q");
    pCosAmp[9][4] = CosFit->GetParameter(0);
    pCosAmpErr[9][4] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_535MeVCM10 = Phip_535MeVCM10_Para->GetAsymmetry(Phip_535MeVCM10_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_535MeVCM10->SetName("ParaPerpAsymmPhip535MeVCM10");
    ParaPerpAsymmPhip_535MeVCM10->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 535 #pm 10MeV (CosTheta-0.8-(-1.0))");
    ParaPerpAsymmPhip_535MeVCM10->Fit("CosFit", "Q");
    pCosAmp[9][5] = CosFit->GetParameter(0);
    pCosAmpErr[9][5] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_555MeVCM10 = Phip_555MeVCM10_Para->GetAsymmetry(Phip_555MeVCM10_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_555MeVCM10->SetName("ParaPerpAsymmPhip555MeVCM10");
    ParaPerpAsymmPhip_555MeVCM10->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 555 #pm 10MeV (CosTheta-0.8-(-1.0))");
    ParaPerpAsymmPhip_555MeVCM10->Fit("CosFit", "Q");
    pCosAmp[9][6] = CosFit->GetParameter(0);
    pCosAmpErr[9][6] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_575MeVCM10 = Phip_575MeVCM10_Para->GetAsymmetry(Phip_575MeVCM10_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_575MeVCM10->SetName("ParaPerpAsymmPhip575MeVCM10");
    ParaPerpAsymmPhip_575MeVCM10->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 575 #pm 10MeV (CosTheta-0.8-(-1.0))");
    ParaPerpAsymmPhip_575MeVCM10->Fit("CosFit", "Q");
    pCosAmp[9][7] = CosFit->GetParameter(0);
    pCosAmpErr[9][7] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_595MeVCM10 = Phip_595MeVCM10_Para->GetAsymmetry(Phip_595MeVCM10_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_595MeVCM10->SetName("ParaPerpAsymmPhip595MeVCM10");
    ParaPerpAsymmPhip_595MeVCM10->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 595 #pm 10MeV (CosTheta-0.8-(-1.0))");
    ParaPerpAsymmPhip_595MeVCM10->Fit("CosFit", "Q");
    pCosAmp[9][8] = CosFit->GetParameter(0);
    pCosAmpErr[9][8] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_615MeVCM10 = Phip_615MeVCM10_Para->GetAsymmetry(Phip_615MeVCM10_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_615MeVCM10->SetName("ParaPerpAsymmPhip615MeVCM10");
    ParaPerpAsymmPhip_615MeVCM10->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 615 #pm 10MeV (CosTheta-0.8-(-1.0))");
    ParaPerpAsymmPhip_615MeVCM10->Fit("CosFit", "Q");
    pCosAmp[9][9] = CosFit->GetParameter(0);
    pCosAmpErr[9][9] = CosFit->GetParError(0);

    TFile f1("ParaPerpAsymm_Total_17.root", "RECREATE");

    ParaPerpAsymmPhip_435MeVCM1->Write();
    ParaPerpAsymmPhip_455MeVCM1->Write();
    ParaPerpAsymmPhip_475MeVCM1->Write();
    ParaPerpAsymmPhip_495MeVCM1->Write();
    ParaPerpAsymmPhip_515MeVCM1->Write();
    ParaPerpAsymmPhip_535MeVCM1->Write();
    ParaPerpAsymmPhip_555MeVCM1->Write();
    ParaPerpAsymmPhip_575MeVCM1->Write();
    ParaPerpAsymmPhip_595MeVCM1->Write();
    ParaPerpAsymmPhip_615MeVCM1->Write();

    ParaPerpAsymmPhip_435MeVCM2->Write();
    ParaPerpAsymmPhip_455MeVCM2->Write();
    ParaPerpAsymmPhip_475MeVCM2->Write();
    ParaPerpAsymmPhip_495MeVCM2->Write();
    ParaPerpAsymmPhip_515MeVCM2->Write();
    ParaPerpAsymmPhip_535MeVCM2->Write();
    ParaPerpAsymmPhip_555MeVCM2->Write();
    ParaPerpAsymmPhip_575MeVCM2->Write();
    ParaPerpAsymmPhip_595MeVCM2->Write();
    ParaPerpAsymmPhip_615MeVCM2->Write();

    ParaPerpAsymmPhip_435MeVCM3->Write();
    ParaPerpAsymmPhip_455MeVCM3->Write();
    ParaPerpAsymmPhip_475MeVCM3->Write();
    ParaPerpAsymmPhip_495MeVCM3->Write();
    ParaPerpAsymmPhip_515MeVCM3->Write();
    ParaPerpAsymmPhip_535MeVCM3->Write();
    ParaPerpAsymmPhip_555MeVCM3->Write();
    ParaPerpAsymmPhip_575MeVCM3->Write();
    ParaPerpAsymmPhip_595MeVCM3->Write();
    ParaPerpAsymmPhip_615MeVCM3->Write();

    ParaPerpAsymmPhip_435MeVCM4->Write();
    ParaPerpAsymmPhip_455MeVCM4->Write();
    ParaPerpAsymmPhip_475MeVCM4->Write();
    ParaPerpAsymmPhip_495MeVCM4->Write();
    ParaPerpAsymmPhip_515MeVCM4->Write();
    ParaPerpAsymmPhip_535MeVCM4->Write();
    ParaPerpAsymmPhip_555MeVCM4->Write();
    ParaPerpAsymmPhip_575MeVCM4->Write();
    ParaPerpAsymmPhip_595MeVCM4->Write();
    ParaPerpAsymmPhip_615MeVCM4->Write();

    ParaPerpAsymmPhip_435MeVCM5->Write();
    ParaPerpAsymmPhip_455MeVCM5->Write();
    ParaPerpAsymmPhip_475MeVCM5->Write();
    ParaPerpAsymmPhip_495MeVCM5->Write();
    ParaPerpAsymmPhip_515MeVCM5->Write();
    ParaPerpAsymmPhip_535MeVCM5->Write();
    ParaPerpAsymmPhip_555MeVCM5->Write();
    ParaPerpAsymmPhip_575MeVCM5->Write();
    ParaPerpAsymmPhip_595MeVCM5->Write();
    ParaPerpAsymmPhip_615MeVCM5->Write();

    ParaPerpAsymmPhip_435MeVCM6->Write();
    ParaPerpAsymmPhip_455MeVCM6->Write();
    ParaPerpAsymmPhip_475MeVCM6->Write();
    ParaPerpAsymmPhip_495MeVCM6->Write();
    ParaPerpAsymmPhip_515MeVCM6->Write();
    ParaPerpAsymmPhip_535MeVCM6->Write();
    ParaPerpAsymmPhip_555MeVCM6->Write();
    ParaPerpAsymmPhip_575MeVCM6->Write();
    ParaPerpAsymmPhip_595MeVCM6->Write();
    ParaPerpAsymmPhip_615MeVCM6->Write();

    ParaPerpAsymmPhip_435MeVCM7->Write();
    ParaPerpAsymmPhip_455MeVCM7->Write();
    ParaPerpAsymmPhip_475MeVCM7->Write();
    ParaPerpAsymmPhip_495MeVCM7->Write();
    ParaPerpAsymmPhip_515MeVCM7->Write();
    ParaPerpAsymmPhip_535MeVCM7->Write();
    ParaPerpAsymmPhip_555MeVCM7->Write();
    ParaPerpAsymmPhip_575MeVCM7->Write();
    ParaPerpAsymmPhip_595MeVCM7->Write();
    ParaPerpAsymmPhip_615MeVCM7->Write();

    ParaPerpAsymmPhip_435MeVCM8->Write();
    ParaPerpAsymmPhip_455MeVCM8->Write();
    ParaPerpAsymmPhip_475MeVCM8->Write();
    ParaPerpAsymmPhip_495MeVCM8->Write();
    ParaPerpAsymmPhip_515MeVCM8->Write();
    ParaPerpAsymmPhip_535MeVCM8->Write();
    ParaPerpAsymmPhip_555MeVCM8->Write();
    ParaPerpAsymmPhip_575MeVCM8->Write();
    ParaPerpAsymmPhip_595MeVCM8->Write();
    ParaPerpAsymmPhip_615MeVCM8->Write();

    ParaPerpAsymmPhip_435MeVCM9->Write();
    ParaPerpAsymmPhip_455MeVCM9->Write();
    ParaPerpAsymmPhip_475MeVCM9->Write();
    ParaPerpAsymmPhip_495MeVCM9->Write();
    ParaPerpAsymmPhip_515MeVCM9->Write();
    ParaPerpAsymmPhip_535MeVCM9->Write();
    ParaPerpAsymmPhip_555MeVCM9->Write();
    ParaPerpAsymmPhip_575MeVCM9->Write();
    ParaPerpAsymmPhip_595MeVCM9->Write();
    ParaPerpAsymmPhip_615MeVCM9->Write();

    ParaPerpAsymmPhip_435MeVCM10->Write();
    ParaPerpAsymmPhip_455MeVCM10->Write();
    ParaPerpAsymmPhip_475MeVCM10->Write();
    ParaPerpAsymmPhip_495MeVCM10->Write();
    ParaPerpAsymmPhip_515MeVCM10->Write();
    ParaPerpAsymmPhip_535MeVCM10->Write();
    ParaPerpAsymmPhip_555MeVCM10->Write();
    ParaPerpAsymmPhip_575MeVCM10->Write();
    ParaPerpAsymmPhip_595MeVCM10->Write();
    ParaPerpAsymmPhip_615MeVCM10->Write();

    //Define new tree to store parameters in
    TTree* tree = new TTree("Parameter_Values", "Tree_of_Values");

    // Define branches to store parameters, (Branch Name, Variable, Type of Variable)
    tree->Branch("pCosAmp435", &pCosA435, "pCosA435/D");
    tree->Branch("pCosAmpErr435", &pCosAErr435, "pCosAErr435/D");
    tree->Branch("pCosAmp455", &pCosA455, "pCosA455/D");
    tree->Branch("pCosAmpErr455", &pCosAErr455, "pCosAErr455/D");
    tree->Branch("pCosAmp475", &pCosA475, "pCosA475/D");
    tree->Branch("pCosAmpErr475", &pCosAErr475, "pCosAErr475/D");
    tree->Branch("pCosAmp495", &pCosA495, "pCosA495/D");
    tree->Branch("pCosAmpErr495", &pCosAErr495, "pCosAErr495/D");
    tree->Branch("pCosAmp515", &pCosA515, "pCosA515/D");
    tree->Branch("pCosAmpErr515", &pCosAErr515, "pCosAErr515/D");
    tree->Branch("pCosAmp535", &pCosA535, "pCosA535/D");
    tree->Branch("pCosAmpErr535", &pCosAErr535, "pCosAErr535/D");
    tree->Branch("pCosAmp555", &pCosA555, "pCosA555/D");
    tree->Branch("pCosAmpErr555", &pCosAErr555, "pCosAErr555/D");
    tree->Branch("pCosAmp575", &pCosA575, "pCosA575/D");
    tree->Branch("pCosAmpErr575", &pCosAErr575, "pCosAErr575/D");
    tree->Branch("pCosAmp595", &pCosA595, "pCosA595/D");
    tree->Branch("pCosAmpErr595", &pCosAErr595, "pCosAErr595/D");
    tree->Branch("pCosAmp615", &pCosA615, "pCosA615/D");
    tree->Branch("pCosAmpErr615", &pCosAErr615, "pCosAErr615/D");

    // Fill branches (and hence tree) with corresponding parameters from above
    for (Int_t m = 0; m < 10; m++){
        pCosA435 = pCosAmp[m][0];
        pCosAErr435 = pCosAmpErr[m][0];
        pCosA455 = pCosAmp[m][1];
        pCosAErr455 = pCosAmpErr[m][1];
        pCosA475 = pCosAmp[m][2];
        pCosAErr475 = pCosAmpErr[m][2];
        pCosA495 = pCosAmp[m][3];
        pCosAErr495 = pCosAmpErr[m][3];
        pCosA515 = pCosAmp[m][4];
        pCosAErr515 = pCosAmpErr[m][4];
        pCosA535 = pCosAmp[m][5];
        pCosAErr535 = pCosAmpErr[m][5];
        pCosA555 = pCosAmp[m][6];
        pCosAErr555 = pCosAmpErr[m][6];
        pCosA575 = pCosAmp[m][7];
        pCosAErr575 = pCosAmpErr[m][7];
        pCosA595 = pCosAmp[m][8];
        pCosAErr595 = pCosAmpErr[m][8];
        pCosA615 = pCosAmp[m][9];
        pCosAErr615 = pCosAmpErr[m][9];
        tree->Fill();
    }

    f1.Write();

}
