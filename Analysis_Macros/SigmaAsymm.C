#include "./includes_SigmaAsymm.h"

void SigmaAsymm(){

    TF1 *CosFunc = new TF1("CosFit", "[0]*cos((2*x*TMath::DegToRad())+acos(0))");
    CosFunc->SetParNames("Amplitude");

    double pCosAmp[5][10]; // Format of array is Theta bin (x) by Egamma bin (y), 10 CosTheta bins, 10 20MeV Egamma bins
    double pCosAmpErr[5][10];
    double pCosA430;
    double pCosAErr430;
    double pCosA450;
    double pCosAErr450;
    double pCosA470;
    double pCosAErr470;
    double pCosA490;
    double pCosAErr490;
    double pCosA510;
    double pCosAErr510;
    double pCosA530;
    double pCosAErr530;
    double pCosA550;
    double pCosAErr550;
    double pCosA570;
    double pCosAErr570;
    double pCosA590;
    double pCosAErr590;
    double pCosA610;
    double pCosAErr610;

    TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/ParaPerp_Total_20_Combined.root"); // Open the latest PTotal combined file to load histograms from
    NPara = Eg_Para->GetEntries();
    NPerp = Eg_Perp->GetEntries();
    ScaleFactor = NPara/NPerp;
    ScaleFactorErr = sqrt( (NPara/((TMath::Power(NPerp,2)))) + (((TMath::Power(NPara,2)))/(TMath::Power(NPerp,3))) ); // Error Propagation of above formula, see notebook 5

    ///////////////////////////////////////////
    ////////////////  CM1  ////////////////////
    ///////////////////////////////////////////

    ParaPerpAsymmPhip_430MeVCM1 = Phip_430MeVCM1_Para->GetAsymmetry(Phip_430MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_430MeVCM1->SetName("ParaPerpAsymmPhip430MeVCM1");
    ParaPerpAsymmPhip_430MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 435 #pm 10MeV (CosTheta1-0.6)");
    ParaPerpAsymmPhip_430MeVCM1->Fit("CosFit", "Q");
    pCosAmp[0][0] = CosFit->GetParameter(0);
    pCosAmpErr[0][0] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_450MeVCM1 = Phip_450MeVCM1_Para->GetAsymmetry(Phip_450MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_450MeVCM1->SetName("ParaPerpAsymmPhip450MeVCM1");
    ParaPerpAsymmPhip_450MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 455 #pm 10MeV (CosTheta1-0.6)");
    ParaPerpAsymmPhip_450MeVCM1->Fit("CosFit", "Q");
    pCosAmp[0][1] = CosFit->GetParameter(0);
    pCosAmpErr[0][1] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_470MeVCM1 = Phip_470MeVCM1_Para->GetAsymmetry(Phip_470MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_470MeVCM1->SetName("ParaPerpAsymmPhip470MeVCM1");
    ParaPerpAsymmPhip_470MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 475 #pm 10MeV (CosTheta1-0.6)");
    ParaPerpAsymmPhip_470MeVCM1->Fit("CosFit", "Q");
    pCosAmp[0][2] = CosFit->GetParameter(0);
    pCosAmpErr[0][2] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_490MeVCM1 = Phip_490MeVCM1_Para->GetAsymmetry(Phip_490MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_490MeVCM1->SetName("ParaPerpAsymmPhip490MeVCM1");
    ParaPerpAsymmPhip_490MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 495 #pm 10MeV (CosTheta1-0.6)");
    ParaPerpAsymmPhip_490MeVCM1->Fit("CosFit", "Q");
    pCosAmp[0][3] = CosFit->GetParameter(0);
    pCosAmpErr[0][3] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_510MeVCM1 = Phip_510MeVCM1_Para->GetAsymmetry(Phip_510MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_510MeVCM1->SetName("ParaPerpAsymmPhip510MeVCM1");
    ParaPerpAsymmPhip_510MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 515 #pm 10MeV (CosTheta1-0.6)");
    ParaPerpAsymmPhip_510MeVCM1->Fit("CosFit", "Q");
    pCosAmp[0][4] = CosFit->GetParameter(0);
    pCosAmpErr[0][4] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_530MeVCM1 = Phip_530MeVCM1_Para->GetAsymmetry(Phip_530MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_530MeVCM1->SetName("ParaPerpAsymmPhip530MeVCM1");
    ParaPerpAsymmPhip_530MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 535 #pm 10MeV (CosTheta1-0.6)");
    ParaPerpAsymmPhip_530MeVCM1->Fit("CosFit", "Q");
    pCosAmp[0][5] = CosFit->GetParameter(0);
    pCosAmpErr[0][5] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_550MeVCM1 = Phip_550MeVCM1_Para->GetAsymmetry(Phip_550MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_550MeVCM1->SetName("ParaPerpAsymmPhip550MeVCM1");
    ParaPerpAsymmPhip_550MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 555 #pm 10MeV (CosTheta1-0.6)");
    ParaPerpAsymmPhip_550MeVCM1->Fit("CosFit", "Q");
    pCosAmp[0][6] = CosFit->GetParameter(0);
    pCosAmpErr[0][6] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_570MeVCM1 = Phip_570MeVCM1_Para->GetAsymmetry(Phip_570MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_570MeVCM1->SetName("ParaPerpAsymmPhip570MeVCM1");
    ParaPerpAsymmPhip_570MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 575 #pm 10MeV (CosTheta1-0.6)");
    ParaPerpAsymmPhip_570MeVCM1->Fit("CosFit", "Q");
    pCosAmp[0][7] = CosFit->GetParameter(0);
    pCosAmpErr[0][7] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_590MeVCM1 = Phip_590MeVCM1_Para->GetAsymmetry(Phip_590MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_590MeVCM1->SetName("ParaPerpAsymmPhip590MeVCM1");
    ParaPerpAsymmPhip_590MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 595 #pm 10MeV (CosTheta1-0.6)");
    ParaPerpAsymmPhip_590MeVCM1->Fit("CosFit", "Q");
    pCosAmp[0][8] = CosFit->GetParameter(0);
    pCosAmpErr[0][8] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_610MeVCM1 = Phip_610MeVCM1_Para->GetAsymmetry(Phip_610MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_610MeVCM1->SetName("ParaPerpAsymmPhip610MeVCM1");
    ParaPerpAsymmPhip_610MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 615 #pm 10MeV (CosTheta1-0.6)");
    ParaPerpAsymmPhip_610MeVCM1->Fit("CosFit", "Q");
    pCosAmp[0][9] = CosFit->GetParameter(0);
    pCosAmpErr[0][9] = CosFit->GetParError(0);

    ///////////////////////////////////////////
    ////////////////  CM2  ////////////////////
    ///////////////////////////////////////////

    ParaPerpAsymmPhip_430MeVCM2 = Phip_430MeVCM2_Para->GetAsymmetry(Phip_430MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_430MeVCM2->SetName("ParaPerpAsymmPhip430MeVCM2");
    ParaPerpAsymmPhip_430MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 435 #pm 10MeV (CosTheta0.6-0.2)");
    ParaPerpAsymmPhip_430MeVCM2->Fit("CosFit", "Q");
    pCosAmp[1][0] = CosFit->GetParameter(0);
    pCosAmpErr[1][0] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_450MeVCM2 = Phip_450MeVCM2_Para->GetAsymmetry(Phip_450MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_450MeVCM2->SetName("ParaPerpAsymmPhip450MeVCM2");
    ParaPerpAsymmPhip_450MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 455 #pm 10MeV (CosTheta0.6-0.2)");
    ParaPerpAsymmPhip_450MeVCM2->Fit("CosFit", "Q");
    pCosAmp[1][1] = CosFit->GetParameter(0);
    pCosAmpErr[1][1] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_470MeVCM2 = Phip_470MeVCM2_Para->GetAsymmetry(Phip_470MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_470MeVCM2->SetName("ParaPerpAsymmPhip470MeVCM2");
    ParaPerpAsymmPhip_470MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 475 #pm 10MeV (CosTheta0.6-0.2)");
    ParaPerpAsymmPhip_470MeVCM2->Fit("CosFit", "Q");
    pCosAmp[1][2] = CosFit->GetParameter(0);
    pCosAmpErr[1][2] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_490MeVCM2 = Phip_490MeVCM2_Para->GetAsymmetry(Phip_490MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_490MeVCM2->SetName("ParaPerpAsymmPhip490MeVCM2");
    ParaPerpAsymmPhip_490MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 495 #pm 10MeV (CosTheta0.6-0.2)");
    ParaPerpAsymmPhip_490MeVCM2->Fit("CosFit", "Q");
    pCosAmp[1][3] = CosFit->GetParameter(0);
    pCosAmpErr[1][3] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_510MeVCM2 = Phip_510MeVCM2_Para->GetAsymmetry(Phip_510MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_510MeVCM2->SetName("ParaPerpAsymmPhip510MeVCM2");
    ParaPerpAsymmPhip_510MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 515 #pm 10MeV (CosTheta0.6-0.2)");
    ParaPerpAsymmPhip_510MeVCM2->Fit("CosFit", "Q");
    pCosAmp[1][4] = CosFit->GetParameter(0);
    pCosAmpErr[1][4] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_530MeVCM2 = Phip_530MeVCM2_Para->GetAsymmetry(Phip_530MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_530MeVCM2->SetName("ParaPerpAsymmPhip530MeVCM2");
    ParaPerpAsymmPhip_530MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 535 #pm 10MeV (CosTheta0.6-0.2)");
    ParaPerpAsymmPhip_530MeVCM2->Fit("CosFit", "Q");
    pCosAmp[1][5] = CosFit->GetParameter(0);
    pCosAmpErr[1][5] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_550MeVCM2 = Phip_550MeVCM2_Para->GetAsymmetry(Phip_550MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_550MeVCM2->SetName("ParaPerpAsymmPhip550MeVCM2");
    ParaPerpAsymmPhip_550MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 555 #pm 10MeV (CosTheta0.6-0.2)");
    ParaPerpAsymmPhip_550MeVCM2->Fit("CosFit", "Q");
    pCosAmp[1][6] = CosFit->GetParameter(0);
    pCosAmpErr[1][6] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_570MeVCM2 = Phip_570MeVCM2_Para->GetAsymmetry(Phip_570MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_570MeVCM2->SetName("ParaPerpAsymmPhip570MeVCM2");
    ParaPerpAsymmPhip_570MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 575 #pm 10MeV (CosTheta0.6-0.2)");
    ParaPerpAsymmPhip_570MeVCM2->Fit("CosFit", "Q");
    pCosAmp[1][7] = CosFit->GetParameter(0);
    pCosAmpErr[1][7] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_590MeVCM2 = Phip_590MeVCM2_Para->GetAsymmetry(Phip_590MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_590MeVCM2->SetName("ParaPerpAsymmPhip590MeVCM2");
    ParaPerpAsymmPhip_590MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 595 #pm 10MeV (CosTheta0.6-0.2)");
    ParaPerpAsymmPhip_590MeVCM2->Fit("CosFit", "Q");
    pCosAmp[1][8] = CosFit->GetParameter(0);
    pCosAmpErr[1][8] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_610MeVCM2 = Phip_610MeVCM2_Para->GetAsymmetry(Phip_610MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_610MeVCM2->SetName("ParaPerpAsymmPhip610MeVCM2");
    ParaPerpAsymmPhip_610MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 615 #pm 10MeV (CosTheta0.6-0.2)");
    ParaPerpAsymmPhip_610MeVCM2->Fit("CosFit", "Q");
    pCosAmp[1][9] = CosFit->GetParameter(0);
    pCosAmpErr[1][9] = CosFit->GetParError(0);

    ///////////////////////////////////////////
    ////////////////  CM3  ////////////////////
    ///////////////////////////////////////////

    ParaPerpAsymmPhip_430MeVCM3 = Phip_430MeVCM3_Para->GetAsymmetry(Phip_430MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_430MeVCM3->SetName("ParaPerpAsymmPhip430MeVCM3");
    ParaPerpAsymmPhip_430MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 435 #pm 10MeV (CosTheta0.2-(-0.2))");
    ParaPerpAsymmPhip_430MeVCM3->Fit("CosFit", "Q");
    pCosAmp[2][0] = CosFit->GetParameter(0);
    pCosAmpErr[2][0] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_450MeVCM3 = Phip_450MeVCM3_Para->GetAsymmetry(Phip_450MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_450MeVCM3->SetName("ParaPerpAsymmPhip450MeVCM3");
    ParaPerpAsymmPhip_450MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 455 #pm 10MeV (CosTheta0.2-(-0.2))");
    ParaPerpAsymmPhip_450MeVCM3->Fit("CosFit", "Q");
    pCosAmp[2][1] = CosFit->GetParameter(0);
    pCosAmpErr[2][1] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_470MeVCM3 = Phip_470MeVCM3_Para->GetAsymmetry(Phip_470MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_470MeVCM3->SetName("ParaPerpAsymmPhip470MeVCM3");
    ParaPerpAsymmPhip_470MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 475 #pm 10MeV (CosTheta0.2-(-0.2))");
    ParaPerpAsymmPhip_470MeVCM3->Fit("CosFit", "Q");
    pCosAmp[2][2] = CosFit->GetParameter(0);
    pCosAmpErr[2][2] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_490MeVCM3 = Phip_490MeVCM3_Para->GetAsymmetry(Phip_490MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_490MeVCM3->SetName("ParaPerpAsymmPhip490MeVCM3");
    ParaPerpAsymmPhip_490MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 495 #pm 10MeV (CosTheta0.2-(-0.2))");
    ParaPerpAsymmPhip_490MeVCM3->Fit("CosFit", "Q");
    pCosAmp[2][3] = CosFit->GetParameter(0);
    pCosAmpErr[2][3] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_510MeVCM3 = Phip_510MeVCM3_Para->GetAsymmetry(Phip_510MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_510MeVCM3->SetName("ParaPerpAsymmPhip510MeVCM3");
    ParaPerpAsymmPhip_510MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 515 #pm 10MeV (CosTheta0.2-(-0.2))");
    ParaPerpAsymmPhip_510MeVCM3->Fit("CosFit", "Q");
    pCosAmp[2][4] = CosFit->GetParameter(0);
    pCosAmpErr[2][4] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_530MeVCM3 = Phip_530MeVCM3_Para->GetAsymmetry(Phip_530MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_530MeVCM3->SetName("ParaPerpAsymmPhip530MeVCM3");
    ParaPerpAsymmPhip_530MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 535 #pm 10MeV (CosTheta0.2-(-0.2))");
    ParaPerpAsymmPhip_530MeVCM3->Fit("CosFit", "Q");
    pCosAmp[2][5] = CosFit->GetParameter(0);
    pCosAmpErr[2][5] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_550MeVCM3 = Phip_550MeVCM3_Para->GetAsymmetry(Phip_550MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_550MeVCM3->SetName("ParaPerpAsymmPhip550MeVCM3");
    ParaPerpAsymmPhip_550MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 555 #pm 10MeV (CosTheta0.2-(-0.2))");
    ParaPerpAsymmPhip_550MeVCM3->Fit("CosFit", "Q");
    pCosAmp[2][6] = CosFit->GetParameter(0);
    pCosAmpErr[2][6] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_570MeVCM3 = Phip_570MeVCM3_Para->GetAsymmetry(Phip_570MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_570MeVCM3->SetName("ParaPerpAsymmPhip570MeVCM3");
    ParaPerpAsymmPhip_570MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 575 #pm 10MeV (CosTheta0.2-(-0.2))");
    ParaPerpAsymmPhip_570MeVCM3->Fit("CosFit", "Q");
    pCosAmp[2][7] = CosFit->GetParameter(0);
    pCosAmpErr[2][7] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_590MeVCM3 = Phip_590MeVCM3_Para->GetAsymmetry(Phip_590MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_590MeVCM3->SetName("ParaPerpAsymmPhip590MeVCM3");
    ParaPerpAsymmPhip_590MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 595 #pm 10MeV (CosTheta0.2-(-0.2))");
    ParaPerpAsymmPhip_590MeVCM3->Fit("CosFit", "Q");
    pCosAmp[2][8] = CosFit->GetParameter(0);
    pCosAmpErr[2][8] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_610MeVCM3 = Phip_610MeVCM3_Para->GetAsymmetry(Phip_610MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_610MeVCM3->SetName("ParaPerpAsymmPhip610MeVCM3");
    ParaPerpAsymmPhip_610MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 615 #pm 10MeV (CosTheta0.2-(-0.2))");
    ParaPerpAsymmPhip_610MeVCM3->Fit("CosFit", "Q");
    pCosAmp[2][9] = CosFit->GetParameter(0);
    pCosAmpErr[2][9] = CosFit->GetParError(0);

    ///////////////////////////////////////////
    ////////////////  CM4  ////////////////////
    ///////////////////////////////////////////

    ParaPerpAsymmPhip_430MeVCM4 = Phip_430MeVCM4_Para->GetAsymmetry(Phip_430MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_430MeVCM4->SetName("ParaPerpAsymmPhip430MeVCM4");
    ParaPerpAsymmPhip_430MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 435 #pm 10MeV (CosTheta-0.2-(-0.6))");
    ParaPerpAsymmPhip_430MeVCM4->Fit("CosFit", "Q");
    pCosAmp[3][0] = CosFit->GetParameter(0);
    pCosAmpErr[3][0] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_450MeVCM4 = Phip_450MeVCM4_Para->GetAsymmetry(Phip_450MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_450MeVCM4->SetName("ParaPerpAsymmPhip450MeVCM4");
    ParaPerpAsymmPhip_450MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 455 #pm 10MeV (CosTheta-0.2-(-0.6))");
    ParaPerpAsymmPhip_450MeVCM4->Fit("CosFit", "Q");
    pCosAmp[3][1] = CosFit->GetParameter(0);
    pCosAmpErr[3][1] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_470MeVCM4 = Phip_470MeVCM4_Para->GetAsymmetry(Phip_470MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_470MeVCM4->SetName("ParaPerpAsymmPhip470MeVCM4");
    ParaPerpAsymmPhip_470MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 475 #pm 10MeV (CosTheta-0.2-(-0.6))");
    ParaPerpAsymmPhip_470MeVCM4->Fit("CosFit", "Q");
    pCosAmp[3][2] = CosFit->GetParameter(0);
    pCosAmpErr[3][2] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_490MeVCM4 = Phip_490MeVCM4_Para->GetAsymmetry(Phip_490MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_490MeVCM4->SetName("ParaPerpAsymmPhip490MeVCM4");
    ParaPerpAsymmPhip_490MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 495 #pm 10MeV (CosTheta-0.2-(-0.6))");
    ParaPerpAsymmPhip_490MeVCM4->Fit("CosFit", "Q");
    pCosAmp[3][3] = CosFit->GetParameter(0);
    pCosAmpErr[3][3] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_510MeVCM4 = Phip_510MeVCM4_Para->GetAsymmetry(Phip_510MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_510MeVCM4->SetName("ParaPerpAsymmPhip510MeVCM4");
    ParaPerpAsymmPhip_510MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 515 #pm 10MeV (CosTheta-0.2-(-0.6))");
    ParaPerpAsymmPhip_510MeVCM4->Fit("CosFit", "Q");
    pCosAmp[3][4] = CosFit->GetParameter(0);
    pCosAmpErr[3][4] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_530MeVCM4 = Phip_530MeVCM4_Para->GetAsymmetry(Phip_530MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_530MeVCM4->SetName("ParaPerpAsymmPhip530MeVCM4");
    ParaPerpAsymmPhip_530MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 535 #pm 10MeV (CosTheta-0.2-(-0.6))");
    ParaPerpAsymmPhip_530MeVCM4->Fit("CosFit", "Q");
    pCosAmp[3][5] = CosFit->GetParameter(0);
    pCosAmpErr[3][5] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_550MeVCM4 = Phip_550MeVCM4_Para->GetAsymmetry(Phip_550MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_550MeVCM4->SetName("ParaPerpAsymmPhip550MeVCM4");
    ParaPerpAsymmPhip_550MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 555 #pm 10MeV (CosTheta-0.2-(-0.6))");
    ParaPerpAsymmPhip_550MeVCM4->Fit("CosFit", "Q");
    pCosAmp[3][6] = CosFit->GetParameter(0);
    pCosAmpErr[3][6] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_570MeVCM4 = Phip_570MeVCM4_Para->GetAsymmetry(Phip_570MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_570MeVCM4->SetName("ParaPerpAsymmPhip570MeVCM4");
    ParaPerpAsymmPhip_570MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 575 #pm 10MeV (CosTheta-0.2-(-0.6))");
    ParaPerpAsymmPhip_570MeVCM4->Fit("CosFit", "Q");
    pCosAmp[3][7] = CosFit->GetParameter(0);
    pCosAmpErr[3][7] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_590MeVCM4 = Phip_590MeVCM4_Para->GetAsymmetry(Phip_590MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_590MeVCM4->SetName("ParaPerpAsymmPhip590MeVCM4");
    ParaPerpAsymmPhip_590MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 595 #pm 10MeV (CosTheta-0.2-(-0.6))");
    ParaPerpAsymmPhip_590MeVCM4->Fit("CosFit", "Q");
    pCosAmp[3][8] = CosFit->GetParameter(0);
    pCosAmpErr[3][8] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_610MeVCM4 = Phip_610MeVCM4_Para->GetAsymmetry(Phip_610MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_610MeVCM4->SetName("ParaPerpAsymmPhip610MeVCM4");
    ParaPerpAsymmPhip_610MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 615 #pm 10MeV (CosTheta-0.2-(-0.6))");
    ParaPerpAsymmPhip_610MeVCM4->Fit("CosFit", "Q");
    pCosAmp[3][9] = CosFit->GetParameter(0);
    pCosAmpErr[3][9] = CosFit->GetParError(0);

    ///////////////////////////////////////////
    ////////////////  CM5  ////////////////////
    ///////////////////////////////////////////

    ParaPerpAsymmPhip_430MeVCM5 = Phip_430MeVCM5_Para->GetAsymmetry(Phip_430MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_430MeVCM5->SetName("ParaPerpAsymmPhip430MeVCM5");
    ParaPerpAsymmPhip_430MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 435 #pm 10MeV (CosTheta-0.6-(-1))");
    ParaPerpAsymmPhip_430MeVCM5->Fit("CosFit", "Q");
    pCosAmp[4][0] = CosFit->GetParameter(0);
    pCosAmpErr[4][0] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_450MeVCM5 = Phip_450MeVCM5_Para->GetAsymmetry(Phip_450MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_450MeVCM5->SetName("ParaPerpAsymmPhip450MeVCM5");
    ParaPerpAsymmPhip_450MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 455 #pm 10MeV (CosTheta-0.6-(-1))");
    ParaPerpAsymmPhip_450MeVCM5->Fit("CosFit", "Q");
    pCosAmp[4][1] = CosFit->GetParameter(0);
    pCosAmpErr[4][1] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_470MeVCM5 = Phip_470MeVCM5_Para->GetAsymmetry(Phip_470MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_470MeVCM5->SetName("ParaPerpAsymmPhip470MeVCM5");
    ParaPerpAsymmPhip_470MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 475 #pm 10MeV (CosTheta-0.6-(-1))");
    ParaPerpAsymmPhip_470MeVCM5->Fit("CosFit", "Q");
    pCosAmp[4][2] = CosFit->GetParameter(0);
    pCosAmpErr[4][2] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_490MeVCM5 = Phip_490MeVCM5_Para->GetAsymmetry(Phip_490MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_490MeVCM5->SetName("ParaPerpAsymmPhip490MeVCM5");
    ParaPerpAsymmPhip_490MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 495 #pm 10MeV (CosTheta-0.6-(-1))");
    ParaPerpAsymmPhip_490MeVCM5->Fit("CosFit", "Q");
    pCosAmp[4][3] = CosFit->GetParameter(0);
    pCosAmpErr[4][3] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_510MeVCM5 = Phip_510MeVCM5_Para->GetAsymmetry(Phip_510MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_510MeVCM5->SetName("ParaPerpAsymmPhip510MeVCM5");
    ParaPerpAsymmPhip_510MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 515 #pm 10MeV (CosTheta-0.6-(-1))");
    ParaPerpAsymmPhip_510MeVCM5->Fit("CosFit", "Q");
    pCosAmp[4][4] = CosFit->GetParameter(0);
    pCosAmpErr[4][4] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_530MeVCM5 = Phip_530MeVCM5_Para->GetAsymmetry(Phip_530MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_530MeVCM5->SetName("ParaPerpAsymmPhip530MeVCM5");
    ParaPerpAsymmPhip_530MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 535 #pm 10MeV (CosTheta-0.6-(-1))");
    ParaPerpAsymmPhip_530MeVCM5->Fit("CosFit", "Q");
    pCosAmp[4][5] = CosFit->GetParameter(0);
    pCosAmpErr[4][5] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_550MeVCM5 = Phip_550MeVCM5_Para->GetAsymmetry(Phip_550MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_550MeVCM5->SetName("ParaPerpAsymmPhip550MeVCM5");
    ParaPerpAsymmPhip_550MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 555 #pm 10MeV (CosTheta-0.6-(-1))");
    ParaPerpAsymmPhip_550MeVCM5->Fit("CosFit", "Q");
    pCosAmp[4][6] = CosFit->GetParameter(0);
    pCosAmpErr[4][6] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_570MeVCM5 = Phip_570MeVCM5_Para->GetAsymmetry(Phip_570MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_570MeVCM5->SetName("ParaPerpAsymmPhip570MeVCM5");
    ParaPerpAsymmPhip_570MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 575 #pm 10MeV (CosTheta-0.6-(-1))");
    ParaPerpAsymmPhip_570MeVCM5->Fit("CosFit", "Q");
    pCosAmp[4][7] = CosFit->GetParameter(0);
    pCosAmpErr[4][7] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_590MeVCM5 = Phip_590MeVCM5_Para->GetAsymmetry(Phip_590MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_590MeVCM5->SetName("ParaPerpAsymmPhip590MeVCM5");
    ParaPerpAsymmPhip_590MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 595 #pm 10MeV (CosTheta-0.6-(-1))");
    ParaPerpAsymmPhip_590MeVCM5->Fit("CosFit", "Q");
    pCosAmp[4][8] = CosFit->GetParameter(0);
    pCosAmpErr[4][8] = CosFit->GetParError(0);

    ParaPerpAsymmPhip_610MeVCM5 = Phip_610MeVCM5_Para->GetAsymmetry(Phip_610MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
    ParaPerpAsymmPhip_610MeVCM5->SetName("ParaPerpAsymmPhip610MeVCM5");
    ParaPerpAsymmPhip_610MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for E_{#gamma} 615 #pm 10MeV (CosTheta-0.6-(-1))");
    ParaPerpAsymmPhip_610MeVCM5->Fit("CosFit", "Q");
    pCosAmp[4][9] = CosFit->GetParameter(0);
    pCosAmpErr[4][9] = CosFit->GetParError(0);

    TFile f1("ParaPerpAsymm_Total_20.root", "RECREATE");

    ParaPerpAsymmPhip_430MeVCM1->Write();
    ParaPerpAsymmPhip_450MeVCM1->Write();
    ParaPerpAsymmPhip_470MeVCM1->Write();
    ParaPerpAsymmPhip_490MeVCM1->Write();
    ParaPerpAsymmPhip_510MeVCM1->Write();
    ParaPerpAsymmPhip_530MeVCM1->Write();
    ParaPerpAsymmPhip_550MeVCM1->Write();
    ParaPerpAsymmPhip_570MeVCM1->Write();
    ParaPerpAsymmPhip_590MeVCM1->Write();
    ParaPerpAsymmPhip_610MeVCM1->Write();

    ParaPerpAsymmPhip_430MeVCM2->Write();
    ParaPerpAsymmPhip_450MeVCM2->Write();
    ParaPerpAsymmPhip_470MeVCM2->Write();
    ParaPerpAsymmPhip_490MeVCM2->Write();
    ParaPerpAsymmPhip_510MeVCM2->Write();
    ParaPerpAsymmPhip_530MeVCM2->Write();
    ParaPerpAsymmPhip_550MeVCM2->Write();
    ParaPerpAsymmPhip_570MeVCM2->Write();
    ParaPerpAsymmPhip_590MeVCM2->Write();
    ParaPerpAsymmPhip_610MeVCM2->Write();

    ParaPerpAsymmPhip_430MeVCM3->Write();
    ParaPerpAsymmPhip_450MeVCM3->Write();
    ParaPerpAsymmPhip_470MeVCM3->Write();
    ParaPerpAsymmPhip_490MeVCM3->Write();
    ParaPerpAsymmPhip_510MeVCM3->Write();
    ParaPerpAsymmPhip_530MeVCM3->Write();
    ParaPerpAsymmPhip_550MeVCM3->Write();
    ParaPerpAsymmPhip_570MeVCM3->Write();
    ParaPerpAsymmPhip_590MeVCM3->Write();
    ParaPerpAsymmPhip_610MeVCM3->Write();

    ParaPerpAsymmPhip_430MeVCM4->Write();
    ParaPerpAsymmPhip_450MeVCM4->Write();
    ParaPerpAsymmPhip_470MeVCM4->Write();
    ParaPerpAsymmPhip_490MeVCM4->Write();
    ParaPerpAsymmPhip_510MeVCM4->Write();
    ParaPerpAsymmPhip_530MeVCM4->Write();
    ParaPerpAsymmPhip_550MeVCM4->Write();
    ParaPerpAsymmPhip_570MeVCM4->Write();
    ParaPerpAsymmPhip_590MeVCM4->Write();
    ParaPerpAsymmPhip_610MeVCM4->Write();

    ParaPerpAsymmPhip_430MeVCM5->Write();
    ParaPerpAsymmPhip_450MeVCM5->Write();
    ParaPerpAsymmPhip_470MeVCM5->Write();
    ParaPerpAsymmPhip_490MeVCM5->Write();
    ParaPerpAsymmPhip_510MeVCM5->Write();
    ParaPerpAsymmPhip_530MeVCM5->Write();
    ParaPerpAsymmPhip_550MeVCM5->Write();
    ParaPerpAsymmPhip_570MeVCM5->Write();
    ParaPerpAsymmPhip_590MeVCM5->Write();
    ParaPerpAsymmPhip_610MeVCM5->Write();

    //Define new tree to store parameters in
    TTree* tree = new TTree("Parameter_Values", "Tree_of_Values");

    // Define branches to store parameters, (Branch Name, Variable, Type of Variable)
    tree->Branch("pCosAmp430", &pCosA430, "pCosA430/D");
    tree->Branch("pCosAmpErr430", &pCosAErr430, "pCosAErr430/D");
    tree->Branch("pCosAmp450", &pCosA450, "pCosA450/D");
    tree->Branch("pCosAmpErr450", &pCosAErr450, "pCosAErr450/D");
    tree->Branch("pCosAmp470", &pCosA470, "pCosA470/D");
    tree->Branch("pCosAmpErr470", &pCosAErr470, "pCosAErr470/D");
    tree->Branch("pCosAmp490", &pCosA490, "pCosA490/D");
    tree->Branch("pCosAmpErr490", &pCosAErr490, "pCosAErr490/D");
    tree->Branch("pCosAmp510", &pCosA510, "pCosA510/D");
    tree->Branch("pCosAmpErr510", &pCosAErr510, "pCosAErr510/D");
    tree->Branch("pCosAmp530", &pCosA530, "pCosA530/D");
    tree->Branch("pCosAmpErr530", &pCosAErr530, "pCosAErr530/D");
    tree->Branch("pCosAmp550", &pCosA550, "pCosA550/D");
    tree->Branch("pCosAmpErr550", &pCosAErr550, "pCosAErr550/D");
    tree->Branch("pCosAmp570", &pCosA570, "pCosA570/D");
    tree->Branch("pCosAmpErr570", &pCosAErr570, "pCosAErr570/D");
    tree->Branch("pCosAmp590", &pCosA590, "pCosA590/D");
    tree->Branch("pCosAmpErr590", &pCosAErr590, "pCosAErr590/D");
    tree->Branch("pCosAmp610", &pCosA610, "pCosA610/D");
    tree->Branch("pCosAmpErr610", &pCosAErr610, "pCosAErr610/D");

    // Fill branches (and hence tree) with corresponding parameters from above
    for (Int_t m = 0; m < 5; m++){
        pCosA430 = pCosAmp[m][0];
        pCosAErr430 = pCosAmpErr[m][0];
        pCosA450 = pCosAmp[m][1];
        pCosAErr450 = pCosAmpErr[m][1];
        pCosA470 = pCosAmp[m][2];
        pCosAErr470 = pCosAmpErr[m][2];
        pCosA490 = pCosAmp[m][3];
        pCosAErr490 = pCosAmpErr[m][3];
        pCosA510 = pCosAmp[m][4];
        pCosAErr510 = pCosAmpErr[m][4];
        pCosA530 = pCosAmp[m][5];
        pCosAErr530 = pCosAmpErr[m][5];
        pCosA550 = pCosAmp[m][6];
        pCosAErr550 = pCosAmpErr[m][6];
        pCosA570 = pCosAmp[m][7];
        pCosAErr570 = pCosAmpErr[m][7];
        pCosA590 = pCosAmp[m][8];
        pCosAErr590 = pCosAmpErr[m][8];
        pCosA610 = pCosAmp[m][9];
        pCosAErr610 = pCosAmpErr[m][9];
        tree->Fill();
    }

    f1.Write();

}
