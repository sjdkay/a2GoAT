#include "./includes.h"

// define a function with 4 parameters
Double_t fitf(Double_t *x,Double_t *par)
{
    Double_t fitval = 0;
    fitval =  ((par[0]*sin(x[0]*TMath::DegToRad()))/(1 + (par[1]*cos(x[0]*TMath::DegToRad()))));
    return fitval;
}

void CxAsymm(){

    double InitialSinAmp[8][7];
    double InitialSinAmpErr[8][7];
    double SinAmp[8][7]; // Format of arrays is Theta Bin Defined by i, Energy bin defined by j
    double SinAmpErr[8][7];
    double CosAmp[8][7]; // Format of arrays is Theta Bin Defined by i, Energy bin defined by j
    double CosAmpErr[8][7];
    Int_t i;
    double ISinAm335;
    double ISinAmErr335;
    double SinAm335;
    double SinAmErr335;
    double CosAm335;
    double CosAmErr335;
    double ISinAm405;
    double ISinAmErr405;
    double SinAm405;
    double SinAmErr405;
    double CosAm405;
    double CosAmErr405;
    double ISinAm475;
    double ISinAmErr475;
    double SinAm475;
    double SinAmErr475;
    double CosAm475;
    double CosAmErr475;
    double ISinAm545;
    double ISinAmErr545;
    double SinAm545;
    double SinAmErr545;
    double CosAm545;
    double CosAmErr545;
    double ISinAm615;
    double ISinAmErr615;
    double SinAm615;
    double SinAmErr615;
    double CosAm615;
    double CosAmErr615;
    double ISinAm685;
    double ISinAmErr685;
    double SinAm685;
    double SinAmErr685;
    double CosAm685;
    double CosAmErr685;

    Double_t Cx335[8], Cx405[8], Cx475[8], Cx545[8], Cx615[8], Cx685[8];
    Double_t CxErr335[8], CxErr405[8], CxErr475[8], CxErr545[8], CxErr615[8], CxErr685[8];

    TF1 *AsymmFunc = new TF1("AsymmFit",  fitf, -130.0, 130.0, 2); //Give a name and range to the fitting funcion
    AsymmFunc->SetParNames("SinAmp", "CosAmp"); //Name the parameters
    AsymmFunc->SetParameter(0, 0);
    TF1 *SinFunc = new TF1("SinFit", "[0]*sin(x*TMath::DegToRad())", -130, 130);
    SinFunc->SetParNames("InitialSinAmp");
    TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/Amo/Physics_Total_59_20_7_17.root"); // Open the latest PTotal file to load histograms from

    ///////////////////////////////////////////
    //////////////////  CM1  //////////////////
    ///////////////////////////////////////////

    PhiSc335AsymmCM1 = Phi_Scattered_335MeV_NegHelCM1->GetAsymmetry(Phi_Scattered_335MeV_PosHelCM1);
    PhiSc335AsymmCM1->SetName("PhiSc335AsymmCM1");
    PhiSc335AsymmCM1->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 300-370MeV (Cos #theta_{CM} 1-0.75)");
    PhiSc335AsymmCM1->Fit("SinFit");
    InitialSinAmp[0][0] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[0][0] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc335AsymmCM1->Fit("AsymmFit");
    SinAmp[0][0] = AsymmFit->GetParameter(0);
    SinAmpErr[0][0] = AsymmFit->GetParError(0);
    CosAmp[0][0] = AsymmFit->GetParameter(1);
    CosAmpErr[0][0] = AsymmFit->GetParError(1);

    PhiSc405AsymmCM1 = Phi_Scattered_405MeV_NegHelCM1->GetAsymmetry(Phi_Scattered_405MeV_PosHelCM1);
    PhiSc405AsymmCM1->SetName("PhiSc405AsymmCM1");
    PhiSc405AsymmCM1->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 370-440MeV (Cos #theta_{CM} 1-0.75)");
    PhiSc405AsymmCM1->Fit("SinFit");
    InitialSinAmp[0][1] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[0][1] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc405AsymmCM1->Fit("AsymmFit");
    SinAmp[0][1] = AsymmFit->GetParameter(0);
    SinAmpErr[0][1] = AsymmFit->GetParError(0);
    CosAmp[0][1] = AsymmFit->GetParameter(1);
    CosAmpErr[0][1] = AsymmFit->GetParError(1);

    PhiSc475AsymmCM1 = Phi_Scattered_475MeV_NegHelCM1->GetAsymmetry(Phi_Scattered_475MeV_PosHelCM1);
    PhiSc475AsymmCM1->SetName("PhiSc475AsymmCM1");
    PhiSc475AsymmCM1->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 440-510MeV (Cos #theta_{CM} 1-0.75)");
    PhiSc475AsymmCM1->Fit("SinFit");
    InitialSinAmp[0][2] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[0][2] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc475AsymmCM1->Fit("AsymmFit");
    SinAmp[0][2] = AsymmFit->GetParameter(0);
    SinAmpErr[0][2] = AsymmFit->GetParError(0);
    CosAmp[0][2] = AsymmFit->GetParameter(1);
    CosAmpErr[0][2] = AsymmFit->GetParError(1);

    PhiSc545AsymmCM1 = Phi_Scattered_545MeV_NegHelCM1->GetAsymmetry(Phi_Scattered_545MeV_PosHelCM1);
    PhiSc545AsymmCM1->SetName("PhiSc545AsymmCM1");
    PhiSc545AsymmCM1->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 510-580MeV (Cos #theta_{CM} 1-0.75)");
    PhiSc545AsymmCM1->Fit("SinFit");
    InitialSinAmp[0][3] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[0][3] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc545AsymmCM1->Fit("AsymmFit");
    SinAmp[0][3] = AsymmFit->GetParameter(0);
    SinAmpErr[0][3] = AsymmFit->GetParError(0);
    CosAmp[0][3] = AsymmFit->GetParameter(1);
    CosAmpErr[0][3] = AsymmFit->GetParError(1);

    PhiSc615AsymmCM1 = Phi_Scattered_615MeV_NegHelCM1->GetAsymmetry(Phi_Scattered_615MeV_PosHelCM1);
    PhiSc615AsymmCM1->SetName("PhiSc615AsymmCM1");
    PhiSc615AsymmCM1->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 580-650MeV (Cos #theta_{CM} 1-0.75)");
    PhiSc615AsymmCM1->Fit("SinFit");
    InitialSinAmp[0][4] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[0][4] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc615AsymmCM1->Fit("AsymmFit");
    SinAmp[0][4] = AsymmFit->GetParameter(0);
    SinAmpErr[0][4] = AsymmFit->GetParError(0);
    CosAmp[0][4] = AsymmFit->GetParameter(1);
    CosAmpErr[0][4] = AsymmFit->GetParError(1);

    PhiSc685AsymmCM1 = Phi_Scattered_685MeV_NegHelCM1->GetAsymmetry(Phi_Scattered_685MeV_PosHelCM1);
    PhiSc685AsymmCM1->SetName("PhiSc685AsymmCM1");
    PhiSc685AsymmCM1->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 650-720MeV (Cos #theta_{CM} 1-0.75)");
    PhiSc685AsymmCM1->Fit("SinFit");
    InitialSinAmp[0][5] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[0][5] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc685AsymmCM1->Fit("AsymmFit");
    SinAmp[0][5] = AsymmFit->GetParameter(0);
    SinAmpErr[0][5] = AsymmFit->GetParError(0);
    CosAmp[0][5] = AsymmFit->GetParameter(1);
    CosAmpErr[0][5] = AsymmFit->GetParError(1);

    ///////////////////////////////////////////
    //////////////////  CM2  //////////////////
    ///////////////////////////////////////////

    PhiSc335AsymmCM2 = Phi_Scattered_335MeV_NegHelCM2->GetAsymmetry(Phi_Scattered_335MeV_PosHelCM2);
    PhiSc335AsymmCM2->SetName("PhiSc335AsymmCM2");
    PhiSc335AsymmCM2->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 300-370MeV (Cos #theta_{CM} 0.75-0.5)");
    PhiSc335AsymmCM2->Fit("SinFit");
    InitialSinAmp[1][0] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[1][0] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc335AsymmCM2->Fit("AsymmFit");
    SinAmp[1][0] = AsymmFit->GetParameter(0);
    SinAmpErr[1][0] = AsymmFit->GetParError(0);
    CosAmp[1][0] = AsymmFit->GetParameter(1);
    CosAmpErr[1][0] = AsymmFit->GetParError(1);

    PhiSc405AsymmCM2 = Phi_Scattered_405MeV_NegHelCM2->GetAsymmetry(Phi_Scattered_405MeV_PosHelCM2);
    PhiSc405AsymmCM2->SetName("PhiSc405AsymmCM2");
    PhiSc405AsymmCM2->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 370-440MeV (Cos #theta_{CM} 0.75-0.5)");
    PhiSc405AsymmCM2->Fit("SinFit");
    InitialSinAmp[1][1] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[1][1] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc405AsymmCM2->Fit("AsymmFit");
    SinAmp[1][1] = AsymmFit->GetParameter(0);
    SinAmpErr[1][1] = AsymmFit->GetParError(0);
    CosAmp[1][1] = AsymmFit->GetParameter(1);
    CosAmpErr[1][1] = AsymmFit->GetParError(1);

    PhiSc475AsymmCM2 = Phi_Scattered_475MeV_NegHelCM2->GetAsymmetry(Phi_Scattered_475MeV_PosHelCM2);
    PhiSc475AsymmCM2->SetName("PhiSc475AsymmCM2");
    PhiSc475AsymmCM2->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 440-510MeV (Cos #theta_{CM} 0.75-0.5)");
    PhiSc475AsymmCM2->Fit("SinFit");
    InitialSinAmp[1][2] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[1][2] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc475AsymmCM2->Fit("AsymmFit");
    SinAmp[1][2] = AsymmFit->GetParameter(0);
    SinAmpErr[1][2] = AsymmFit->GetParError(0);
    CosAmp[1][2] = AsymmFit->GetParameter(1);
    CosAmpErr[1][2] = AsymmFit->GetParError(1);

    PhiSc545AsymmCM2 = Phi_Scattered_545MeV_NegHelCM2->GetAsymmetry(Phi_Scattered_545MeV_PosHelCM2);
    PhiSc545AsymmCM2->SetName("PhiSc545AsymmCM2");
    PhiSc545AsymmCM2->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 510-580MeV (Cos #theta_{CM} 0.75-0.5)");
    PhiSc545AsymmCM2->Fit("SinFit");
    InitialSinAmp[1][3] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[1][3] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc545AsymmCM2->Fit("AsymmFit");
    SinAmp[1][3] = AsymmFit->GetParameter(0);
    SinAmpErr[1][3] = AsymmFit->GetParError(0);
    CosAmp[1][3] = AsymmFit->GetParameter(1);
    CosAmpErr[1][3] = AsymmFit->GetParError(1);

    PhiSc615AsymmCM2 = Phi_Scattered_615MeV_NegHelCM2->GetAsymmetry(Phi_Scattered_615MeV_PosHelCM2);
    PhiSc615AsymmCM2->SetName("PhiSc615AsymmCM2");
    PhiSc615AsymmCM2->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 580-650MeV (Cos #theta_{CM} 0.75-0.5)");
    PhiSc615AsymmCM2->Fit("SinFit");
    InitialSinAmp[1][4] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[1][4] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc615AsymmCM2->Fit("AsymmFit");
    SinAmp[1][4] = AsymmFit->GetParameter(0);
    SinAmpErr[1][4] = AsymmFit->GetParError(0);
    CosAmp[1][4] = AsymmFit->GetParameter(1);
    CosAmpErr[1][4] = AsymmFit->GetParError(1);

    PhiSc685AsymmCM2 = Phi_Scattered_685MeV_NegHelCM2->GetAsymmetry(Phi_Scattered_685MeV_PosHelCM2);
    PhiSc685AsymmCM2->SetName("PhiSc685AsymmCM2");
    PhiSc685AsymmCM2->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 650-720MeV (Cos #theta_{CM} 0.75-0.5)");
    PhiSc685AsymmCM2->Fit("SinFit");
    InitialSinAmp[1][5] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[1][5] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc685AsymmCM2->Fit("AsymmFit");
    SinAmp[1][5] = AsymmFit->GetParameter(0);
    SinAmpErr[1][5] = AsymmFit->GetParError(0);
    CosAmp[1][5] = AsymmFit->GetParameter(1);
    CosAmpErr[1][5] = AsymmFit->GetParError(1);

    ///////////////////////////////////////////
    //////////////////  CM3  //////////////////
    ///////////////////////////////////////////

    PhiSc335AsymmCM3 = Phi_Scattered_335MeV_NegHelCM3->GetAsymmetry(Phi_Scattered_335MeV_PosHelCM3);
    PhiSc335AsymmCM3->SetName("PhiSc335AsymmCM3");
    PhiSc335AsymmCM3->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 300-370MeV (Cos #theta_{CM} 0.5-0.25)");
    PhiSc335AsymmCM3->Fit("SinFit");
    InitialSinAmp[2][0] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[2][0] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc335AsymmCM3->Fit("AsymmFit");
    SinAmp[2][0] = AsymmFit->GetParameter(0);
    SinAmpErr[2][0] = AsymmFit->GetParError(0);
    CosAmp[2][0] = AsymmFit->GetParameter(1);
    CosAmpErr[2][0] = AsymmFit->GetParError(1);

    PhiSc405AsymmCM3 = Phi_Scattered_405MeV_NegHelCM3->GetAsymmetry(Phi_Scattered_405MeV_PosHelCM3);
    PhiSc405AsymmCM3->SetName("PhiSc405AsymmCM3");
    PhiSc405AsymmCM3->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 370-440MeV (Cos #theta_{CM} 0.5-0.25)");
    PhiSc405AsymmCM3->Fit("SinFit");
    InitialSinAmp[2][1] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[2][1] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc405AsymmCM3->Fit("AsymmFit");
    SinAmp[2][1] = AsymmFit->GetParameter(0);
    SinAmpErr[2][1] = AsymmFit->GetParError(0);
    CosAmp[2][1] = AsymmFit->GetParameter(1);
    CosAmpErr[2][1] = AsymmFit->GetParError(1);

    PhiSc475AsymmCM3 = Phi_Scattered_475MeV_NegHelCM3->GetAsymmetry(Phi_Scattered_475MeV_PosHelCM3);
    PhiSc475AsymmCM3->SetName("PhiSc475AsymmCM3");
    PhiSc475AsymmCM3->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 440-510MeV (Cos #theta_{CM} 0.5-0.25)");
    PhiSc475AsymmCM3->Fit("SinFit");
    InitialSinAmp[2][2] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[2][2] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc475AsymmCM3->Fit("AsymmFit");
    SinAmp[2][2] = AsymmFit->GetParameter(0);
    SinAmpErr[2][2] = AsymmFit->GetParError(0);
    CosAmp[2][2] = AsymmFit->GetParameter(1);
    CosAmpErr[2][2] = AsymmFit->GetParError(1);

    PhiSc545AsymmCM3 = Phi_Scattered_545MeV_NegHelCM3->GetAsymmetry(Phi_Scattered_545MeV_PosHelCM3);
    PhiSc545AsymmCM3->SetName("PhiSc545AsymmCM3");
    PhiSc545AsymmCM3->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 510-580MeV (Cos #theta_{CM} 0.5-0.25)");
    PhiSc545AsymmCM3->Fit("SinFit");
    InitialSinAmp[2][3] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[2][3] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc545AsymmCM3->Fit("AsymmFit");
    SinAmp[2][3] = AsymmFit->GetParameter(0);
    SinAmpErr[2][3] = AsymmFit->GetParError(0);
    CosAmp[2][3] = AsymmFit->GetParameter(1);
    CosAmpErr[2][3] = AsymmFit->GetParError(1);

    PhiSc615AsymmCM3 = Phi_Scattered_615MeV_NegHelCM3->GetAsymmetry(Phi_Scattered_615MeV_PosHelCM3);
    PhiSc615AsymmCM3->SetName("PhiSc615AsymmCM3");
    PhiSc615AsymmCM3->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 580-650MeV (Cos #theta_{CM} 0.5-0.25)");
    PhiSc615AsymmCM3->Fit("SinFit");
    InitialSinAmp[2][4] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[2][4] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc615AsymmCM3->Fit("AsymmFit");
    SinAmp[2][4] = AsymmFit->GetParameter(0);
    SinAmpErr[2][4] = AsymmFit->GetParError(0);
    CosAmp[2][4] = AsymmFit->GetParameter(1);
    CosAmpErr[2][4] = AsymmFit->GetParError(1);

    PhiSc685AsymmCM3 = Phi_Scattered_685MeV_NegHelCM3->GetAsymmetry(Phi_Scattered_685MeV_PosHelCM3);
    PhiSc685AsymmCM3->SetName("PhiSc685AsymmCM3");
    PhiSc685AsymmCM3->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 650-720MeV (Cos #theta_{CM} 0.5-0.25)");
    PhiSc685AsymmCM3->Fit("SinFit");
    InitialSinAmp[2][5] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[2][5] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc685AsymmCM3->Fit("AsymmFit");
    SinAmp[2][5] = AsymmFit->GetParameter(0);
    SinAmpErr[2][5] = AsymmFit->GetParError(0);
    CosAmp[2][5] = AsymmFit->GetParameter(1);
    CosAmpErr[2][5] = AsymmFit->GetParError(1);

    ///////////////////////////////////////////
    //////////////////  CM4  //////////////////
    ///////////////////////////////////////////

    PhiSc335AsymmCM4 = Phi_Scattered_335MeV_NegHelCM4->GetAsymmetry(Phi_Scattered_335MeV_PosHelCM4);
    PhiSc335AsymmCM4->SetName("PhiSc335AsymmCM4");
    PhiSc335AsymmCM4->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 300-370MeV (Cos #theta_{CM} 0.25-0.0)");
    PhiSc335AsymmCM4->Fit("SinFit");
    InitialSinAmp[3][0] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[3][0] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc335AsymmCM4->Fit("AsymmFit");
    SinAmp[3][0] = AsymmFit->GetParameter(0);
    SinAmpErr[3][0] = AsymmFit->GetParError(0);
    CosAmp[3][0] = AsymmFit->GetParameter(1);
    CosAmpErr[3][0] = AsymmFit->GetParError(1);

    PhiSc405AsymmCM4 = Phi_Scattered_405MeV_NegHelCM4->GetAsymmetry(Phi_Scattered_405MeV_PosHelCM4);
    PhiSc405AsymmCM4->SetName("PhiSc405AsymmCM4");
    PhiSc405AsymmCM4->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 370-440MeV (Cos #theta_{CM} 0.25-0.0)");
    PhiSc405AsymmCM4->Fit("SinFit");
    InitialSinAmp[3][1] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[3][1] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc405AsymmCM4->Fit("AsymmFit");
    SinAmp[3][1] = AsymmFit->GetParameter(0);
    SinAmpErr[3][1] = AsymmFit->GetParError(0);
    CosAmp[3][1] = AsymmFit->GetParameter(1);
    CosAmpErr[3][1] = AsymmFit->GetParError(1);

    PhiSc475AsymmCM4 = Phi_Scattered_475MeV_NegHelCM4->GetAsymmetry(Phi_Scattered_475MeV_PosHelCM4);
    PhiSc475AsymmCM4->SetName("PhiSc475AsymmCM4");
    PhiSc475AsymmCM4->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 440-510MeV (Cos #theta_{CM} 0.25-0.0)");
    PhiSc475AsymmCM4->Fit("SinFit");
    InitialSinAmp[3][2] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[3][2] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc475AsymmCM4->Fit("AsymmFit");
    SinAmp[3][2] = AsymmFit->GetParameter(0);
    SinAmpErr[3][2] = AsymmFit->GetParError(0);
    CosAmp[3][2] = AsymmFit->GetParameter(1);
    CosAmpErr[3][2] = AsymmFit->GetParError(1);

    PhiSc545AsymmCM4 = Phi_Scattered_545MeV_NegHelCM4->GetAsymmetry(Phi_Scattered_545MeV_PosHelCM4);
    PhiSc545AsymmCM4->SetName("PhiSc545AsymmCM4");
    PhiSc545AsymmCM4->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 510-580MeV (Cos #theta_{CM} 0.25-0.0)");
    PhiSc545AsymmCM4->Fit("SinFit");
    InitialSinAmp[3][3] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[3][3] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc545AsymmCM4->Fit("AsymmFit");
    SinAmp[3][3] = AsymmFit->GetParameter(0);
    SinAmpErr[3][3] = AsymmFit->GetParError(0);
    CosAmp[3][3] = AsymmFit->GetParameter(1);
    CosAmpErr[3][3] = AsymmFit->GetParError(1);

    PhiSc615AsymmCM4 = Phi_Scattered_615MeV_NegHelCM4->GetAsymmetry(Phi_Scattered_615MeV_PosHelCM4);
    PhiSc615AsymmCM4->SetName("PhiSc615AsymmCM4");
    PhiSc615AsymmCM4->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 580-650MeV (Cos #theta_{CM} 0.25-0.0)");
    PhiSc615AsymmCM4->Fit("SinFit");
    InitialSinAmp[3][4] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[3][4] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc615AsymmCM4->Fit("AsymmFit");
    SinAmp[3][4] = AsymmFit->GetParameter(0);
    SinAmpErr[3][4] = AsymmFit->GetParError(0);
    CosAmp[3][4] = AsymmFit->GetParameter(1);
    CosAmpErr[3][4] = AsymmFit->GetParError(1);

    PhiSc685AsymmCM4 = Phi_Scattered_685MeV_NegHelCM4->GetAsymmetry(Phi_Scattered_685MeV_PosHelCM4);
    PhiSc685AsymmCM4->SetName("PhiSc685AsymmCM4");
    PhiSc685AsymmCM4->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 650-720MeV (Cos #theta_{CM} 0.25-0.0)");
    PhiSc685AsymmCM4->Fit("SinFit");
    InitialSinAmp[3][5] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[3][5] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc685AsymmCM4->Fit("AsymmFit");
    SinAmp[3][5] = AsymmFit->GetParameter(0);
    SinAmpErr[3][5] = AsymmFit->GetParError(0);
    CosAmp[3][5] = AsymmFit->GetParameter(1);
    CosAmpErr[3][5] = AsymmFit->GetParError(1);

    ///////////////////////////////////////////
    //////////////////  CM5  //////////////////
    ///////////////////////////////////////////

    PhiSc335AsymmCM5 = Phi_Scattered_335MeV_NegHelCM5->GetAsymmetry(Phi_Scattered_335MeV_PosHelCM5);
    PhiSc335AsymmCM5->SetName("PhiSc335AsymmCM5");
    PhiSc335AsymmCM5->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 300-370MeV (Cos #theta_{CM} 0.0-(-0.25))");
    PhiSc335AsymmCM5->Fit("SinFit");
    InitialSinAmp[4][0] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[4][0] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc335AsymmCM5->Fit("AsymmFit");
    SinAmp[4][0] = AsymmFit->GetParameter(0);
    SinAmpErr[4][0] = AsymmFit->GetParError(0);
    CosAmp[4][0] = AsymmFit->GetParameter(1);
    CosAmpErr[4][0] = AsymmFit->GetParError(1);

    PhiSc405AsymmCM5 = Phi_Scattered_405MeV_NegHelCM5->GetAsymmetry(Phi_Scattered_405MeV_PosHelCM5);
    PhiSc405AsymmCM5->SetName("PhiSc405AsymmCM5");
    PhiSc405AsymmCM5->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 370-440MeV (Cos #theta_{CM} 0.0-(-0.25))");
    PhiSc405AsymmCM5->Fit("SinFit");
    InitialSinAmp[4][1] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[4][1] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc405AsymmCM5->Fit("AsymmFit");
    SinAmp[4][1] = AsymmFit->GetParameter(0);
    SinAmpErr[4][1] = AsymmFit->GetParError(0);
    CosAmp[4][1] = AsymmFit->GetParameter(1);
    CosAmpErr[4][1] = AsymmFit->GetParError(1);

    PhiSc475AsymmCM5 = Phi_Scattered_475MeV_NegHelCM5->GetAsymmetry(Phi_Scattered_475MeV_PosHelCM5);
    PhiSc475AsymmCM5->SetName("PhiSc475AsymmCM5");
    PhiSc475AsymmCM5->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 440-510MeV (Cos #theta_{CM} 0.0-(-0.25))");
    PhiSc475AsymmCM5->Fit("SinFit");
    InitialSinAmp[4][2] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[4][2] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc475AsymmCM5->Fit("AsymmFit");
    SinAmp[4][2] = AsymmFit->GetParameter(0);
    SinAmpErr[4][2] = AsymmFit->GetParError(0);
    CosAmp[4][2] = AsymmFit->GetParameter(1);
    CosAmpErr[4][2] = AsymmFit->GetParError(1);

    PhiSc545AsymmCM5 = Phi_Scattered_545MeV_NegHelCM5->GetAsymmetry(Phi_Scattered_545MeV_PosHelCM5);
    PhiSc545AsymmCM5->SetName("PhiSc545AsymmCM5");
    PhiSc545AsymmCM5->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 510-580MeV (Cos #theta_{CM} 0.0-(-0.25))");
    PhiSc545AsymmCM5->Fit("SinFit");
    InitialSinAmp[4][3] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[4][3] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc545AsymmCM5->Fit("AsymmFit");
    SinAmp[4][3] = AsymmFit->GetParameter(0);
    SinAmpErr[4][3] = AsymmFit->GetParError(0);
    CosAmp[4][3] = AsymmFit->GetParameter(1);
    CosAmpErr[4][3] = AsymmFit->GetParError(1);

    PhiSc615AsymmCM5 = Phi_Scattered_615MeV_NegHelCM5->GetAsymmetry(Phi_Scattered_615MeV_PosHelCM5);
    PhiSc615AsymmCM5->SetName("PhiSc615AsymmCM5");
    PhiSc615AsymmCM5->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 580-650MeV (Cos #theta_{CM} 0.0-(-0.25))");
    PhiSc615AsymmCM5->Fit("SinFit");
    InitialSinAmp[4][4] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[4][4] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc615AsymmCM5->Fit("AsymmFit");
    SinAmp[4][4] = AsymmFit->GetParameter(0);
    SinAmpErr[4][4] = AsymmFit->GetParError(0);
    CosAmp[4][4] = AsymmFit->GetParameter(1);
    CosAmpErr[4][4] = AsymmFit->GetParError(1);

    PhiSc685AsymmCM5 = Phi_Scattered_685MeV_NegHelCM5->GetAsymmetry(Phi_Scattered_685MeV_PosHelCM5);
    PhiSc685AsymmCM5->SetName("PhiSc685AsymmCM5");
    PhiSc685AsymmCM5->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 650-720MeV (Cos #theta_{CM} 0.0-(-0.25))");
    PhiSc685AsymmCM5->Fit("SinFit");
    InitialSinAmp[4][5] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[4][5] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc685AsymmCM5->Fit("AsymmFit");
    SinAmp[4][5] = AsymmFit->GetParameter(0);
    SinAmpErr[4][5] = AsymmFit->GetParError(0);
    CosAmp[4][5] = AsymmFit->GetParameter(1);
    CosAmpErr[4][5] = AsymmFit->GetParError(1);

    ///////////////////////////////////////////
    //////////////////  CM6  //////////////////
    ///////////////////////////////////////////

    PhiSc335AsymmCM6 = Phi_Scattered_335MeV_NegHelCM6->GetAsymmetry(Phi_Scattered_335MeV_PosHelCM6);
    PhiSc335AsymmCM6->SetName("PhiSc335AsymmCM6");
    PhiSc335AsymmCM6->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 300-370MeV (Cos #theta_{CM} -0.25-(-0.5))");
    PhiSc335AsymmCM6->Fit("SinFit");
    InitialSinAmp[5][0] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[5][0] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc335AsymmCM6->Fit("AsymmFit");
    SinAmp[5][0] = AsymmFit->GetParameter(0);
    SinAmpErr[5][0] = AsymmFit->GetParError(0);
    CosAmp[5][0] = AsymmFit->GetParameter(1);
    CosAmpErr[5][0] = AsymmFit->GetParError(1);

    PhiSc405AsymmCM6 = Phi_Scattered_405MeV_NegHelCM6->GetAsymmetry(Phi_Scattered_405MeV_PosHelCM6);
    PhiSc405AsymmCM6->SetName("PhiSc405AsymmCM6");
    PhiSc405AsymmCM6->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 370-440MeV (Cos #theta_{CM} -0.25-(-0.5))");
    PhiSc405AsymmCM6->Fit("SinFit");
    InitialSinAmp[5][1] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[5][1] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc405AsymmCM6->Fit("AsymmFit");
    SinAmp[5][1] = AsymmFit->GetParameter(0);
    SinAmpErr[5][1] = AsymmFit->GetParError(0);
    CosAmp[5][1] = AsymmFit->GetParameter(1);
    CosAmpErr[5][1] = AsymmFit->GetParError(1);

    PhiSc475AsymmCM6 = Phi_Scattered_475MeV_NegHelCM6->GetAsymmetry(Phi_Scattered_475MeV_PosHelCM6);
    PhiSc475AsymmCM6->SetName("PhiSc475AsymmCM6");
    PhiSc475AsymmCM6->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 440-510MeV (Cos #theta_{CM} -0.25-(-0.5))");
    PhiSc475AsymmCM6->Fit("SinFit");
    InitialSinAmp[5][2] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[5][2] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc475AsymmCM6->Fit("AsymmFit");
    SinAmp[5][2] = AsymmFit->GetParameter(0);
    SinAmpErr[5][2] = AsymmFit->GetParError(0);
    CosAmp[5][2] = AsymmFit->GetParameter(1);
    CosAmpErr[5][2] = AsymmFit->GetParError(1);

    PhiSc545AsymmCM6 = Phi_Scattered_545MeV_NegHelCM6->GetAsymmetry(Phi_Scattered_545MeV_PosHelCM6);
    PhiSc545AsymmCM6->SetName("PhiSc545AsymmCM6");
    PhiSc545AsymmCM6->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 510-580MeV (Cos #theta_{CM} -0.25-(-0.5))");
    PhiSc545AsymmCM6->Fit("SinFit");
    InitialSinAmp[5][3] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[5][3] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc545AsymmCM6->Fit("AsymmFit");
    SinAmp[5][3] = AsymmFit->GetParameter(0);
    SinAmpErr[5][3] = AsymmFit->GetParError(0);
    CosAmp[5][3] = AsymmFit->GetParameter(1);
    CosAmpErr[5][3] = AsymmFit->GetParError(1);

    PhiSc615AsymmCM6 = Phi_Scattered_615MeV_NegHelCM6->GetAsymmetry(Phi_Scattered_615MeV_PosHelCM6);
    PhiSc615AsymmCM6->SetName("PhiSc615AsymmCM6");
    PhiSc615AsymmCM6->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 580-650MeV (Cos #theta_{CM} -0.25-(-0.5))");
    PhiSc615AsymmCM6->Fit("SinFit");
    InitialSinAmp[5][4] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[5][4] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc615AsymmCM6->Fit("AsymmFit");
    SinAmp[5][4] = AsymmFit->GetParameter(0);
    SinAmpErr[5][4] = AsymmFit->GetParError(0);
    CosAmp[5][4] = AsymmFit->GetParameter(1);
    CosAmpErr[5][4] = AsymmFit->GetParError(1);

    PhiSc685AsymmCM6 = Phi_Scattered_685MeV_NegHelCM6->GetAsymmetry(Phi_Scattered_685MeV_PosHelCM6);
    PhiSc685AsymmCM6->SetName("PhiSc685AsymmCM6");
    PhiSc685AsymmCM6->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 650-720MeV (Cos #theta_{CM} -0.25-(-0.5))");
    PhiSc685AsymmCM6->Fit("SinFit");
    InitialSinAmp[5][5] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[5][5] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc685AsymmCM6->Fit("AsymmFit");
    SinAmp[5][5] = AsymmFit->GetParameter(0);
    SinAmpErr[5][5] = AsymmFit->GetParError(0);
    CosAmp[5][5] = AsymmFit->GetParameter(1);
    CosAmpErr[5][5] = AsymmFit->GetParError(1);

    ///////////////////////////////////////////
    //////////////////  CM7  //////////////////
    ///////////////////////////////////////////

    PhiSc335AsymmCM7 = Phi_Scattered_335MeV_NegHelCM7->GetAsymmetry(Phi_Scattered_335MeV_PosHelCM7);
    PhiSc335AsymmCM7->SetName("PhiSc335AsymmCM7");
    PhiSc335AsymmCM7->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 300-370MeV (Cos #theta_{CM} -0.5-(-0.75))");
    PhiSc335AsymmCM7->Fit("SinFit");
    InitialSinAmp[6][0] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[6][0] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc335AsymmCM7->Fit("AsymmFit");
    SinAmp[6][0] = AsymmFit->GetParameter(0);
    SinAmpErr[6][0] = AsymmFit->GetParError(0);
    CosAmp[6][0] = AsymmFit->GetParameter(1);
    CosAmpErr[6][0] = AsymmFit->GetParError(1);

    PhiSc405AsymmCM7 = Phi_Scattered_405MeV_NegHelCM7->GetAsymmetry(Phi_Scattered_405MeV_PosHelCM7);
    PhiSc405AsymmCM7->SetName("PhiSc405AsymmCM7");
    PhiSc405AsymmCM7->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 370-440MeV (Cos #theta_{CM} -0.5-(-0.75))");
    PhiSc405AsymmCM7->Fit("SinFit");
    InitialSinAmp[6][1] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[6][1] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc405AsymmCM7->Fit("AsymmFit");
    SinAmp[6][1] = AsymmFit->GetParameter(0);
    SinAmpErr[6][1] = AsymmFit->GetParError(0);
    CosAmp[6][1] = AsymmFit->GetParameter(1);
    CosAmpErr[6][1] = AsymmFit->GetParError(1);

    PhiSc475AsymmCM7 = Phi_Scattered_475MeV_NegHelCM7->GetAsymmetry(Phi_Scattered_475MeV_PosHelCM7);
    PhiSc475AsymmCM7->SetName("PhiSc475AsymmCM7");
    PhiSc475AsymmCM7->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 440-510MeV (Cos #theta_{CM} -0.5-(-0.75))");
    PhiSc475AsymmCM7->Fit("SinFit");
    InitialSinAmp[6][2] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[6][2] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc475AsymmCM7->Fit("AsymmFit");
    SinAmp[6][2] = AsymmFit->GetParameter(0);
    SinAmpErr[6][2] = AsymmFit->GetParError(0);
    CosAmp[6][2] = AsymmFit->GetParameter(1);
    CosAmpErr[6][2] = AsymmFit->GetParError(1);

    PhiSc545AsymmCM7 = Phi_Scattered_545MeV_NegHelCM7->GetAsymmetry(Phi_Scattered_545MeV_PosHelCM7);
    PhiSc545AsymmCM7->SetName("PhiSc545AsymmCM7");
    PhiSc545AsymmCM7->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 510-580MeV (Cos #theta_{CM} -0.5-(-0.75))");
    PhiSc545AsymmCM7->Fit("SinFit");
    InitialSinAmp[6][3] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[6][3] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc545AsymmCM7->Fit("AsymmFit");
    SinAmp[6][3] = AsymmFit->GetParameter(0);
    SinAmpErr[6][3] = AsymmFit->GetParError(0);
    CosAmp[6][3] = AsymmFit->GetParameter(1);
    CosAmpErr[6][3] = AsymmFit->GetParError(1);

    PhiSc615AsymmCM7 = Phi_Scattered_615MeV_NegHelCM7->GetAsymmetry(Phi_Scattered_615MeV_PosHelCM7);
    PhiSc615AsymmCM7->SetName("PhiSc615AsymmCM7");
    PhiSc615AsymmCM7->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 580-650MeV (Cos #theta_{CM} -0.5-(-0.75))");
    PhiSc615AsymmCM7->Fit("SinFit");
    InitialSinAmp[6][4] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[6][4] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc615AsymmCM7->Fit("AsymmFit");
    SinAmp[6][4] = AsymmFit->GetParameter(0);
    SinAmpErr[6][4] = AsymmFit->GetParError(0);
    CosAmp[6][4] = AsymmFit->GetParameter(1);
    CosAmpErr[6][4] = AsymmFit->GetParError(1);

    PhiSc685AsymmCM7 = Phi_Scattered_685MeV_NegHelCM7->GetAsymmetry(Phi_Scattered_685MeV_PosHelCM7);
    PhiSc685AsymmCM7->SetName("PhiSc685AsymmCM7");
    PhiSc685AsymmCM7->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 650-720MeV (Cos #theta_{CM} -0.5-(-0.75))");
    PhiSc685AsymmCM7->Fit("SinFit");
    InitialSinAmp[6][5] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[6][5] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc685AsymmCM7->Fit("AsymmFit");
    SinAmp[6][5] = AsymmFit->GetParameter(0);
    SinAmpErr[6][5] = AsymmFit->GetParError(0);
    CosAmp[6][5] = AsymmFit->GetParameter(1);
    CosAmpErr[6][5] = AsymmFit->GetParError(1);

    ///////////////////////////////////////////
    //////////////////  CM8  //////////////////
    ///////////////////////////////////////////

    PhiSc335AsymmCM8 = Phi_Scattered_335MeV_NegHelCM8->GetAsymmetry(Phi_Scattered_335MeV_PosHelCM8);
    PhiSc335AsymmCM8->SetName("PhiSc335AsymmCM8");
    PhiSc335AsymmCM8->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 300-370MeV (Cos #theta_{CM} -0.75-(-1.0))");
    PhiSc335AsymmCM8->Fit("SinFit");
    InitialSinAmp[7][0] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[7][0] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc335AsymmCM8->Fit("AsymmFit");
    SinAmp[7][0] = AsymmFit->GetParameter(0);
    SinAmpErr[7][0] = AsymmFit->GetParError(0);
    CosAmp[7][0] = AsymmFit->GetParameter(1);
    CosAmpErr[7][0] = AsymmFit->GetParError(1);

    PhiSc405AsymmCM8 = Phi_Scattered_405MeV_NegHelCM8->GetAsymmetry(Phi_Scattered_405MeV_PosHelCM8);
    PhiSc405AsymmCM8->SetName("PhiSc405AsymmCM8");
    PhiSc405AsymmCM8->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 370-440MeV (Cos #theta_{CM} -0.75-(-1.0))");
    PhiSc405AsymmCM8->Fit("SinFit");
    InitialSinAmp[7][1] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[7][1] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc405AsymmCM8->Fit("AsymmFit");
    SinAmp[7][1] = AsymmFit->GetParameter(0);
    SinAmpErr[7][1] = AsymmFit->GetParError(0);
    CosAmp[7][1] = AsymmFit->GetParameter(1);
    CosAmpErr[7][1] = AsymmFit->GetParError(1);

    PhiSc475AsymmCM8 = Phi_Scattered_475MeV_NegHelCM8->GetAsymmetry(Phi_Scattered_475MeV_PosHelCM8);
    PhiSc475AsymmCM8->SetName("PhiSc475AsymmCM8");
    PhiSc475AsymmCM8->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 440-510MeV (Cos #theta_{CM} -0.75-(-1.0))");
    PhiSc475AsymmCM8->Fit("SinFit");
    InitialSinAmp[7][2] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[7][2] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc475AsymmCM8->Fit("AsymmFit");
    SinAmp[7][2] = AsymmFit->GetParameter(0);
    SinAmpErr[7][2] = AsymmFit->GetParError(0);
    CosAmp[7][2] = AsymmFit->GetParameter(1);
    CosAmpErr[7][2] = AsymmFit->GetParError(1);

    PhiSc545AsymmCM8 = Phi_Scattered_545MeV_NegHelCM8->GetAsymmetry(Phi_Scattered_545MeV_PosHelCM8);
    PhiSc545AsymmCM8->SetName("PhiSc545AsymmCM8");
    PhiSc545AsymmCM8->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 510-580MeV (Cos #theta_{CM} -0.75-(-1.0))");
    PhiSc545AsymmCM8->Fit("SinFit");
    InitialSinAmp[7][3] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[7][3] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc545AsymmCM8->Fit("AsymmFit");
    SinAmp[7][3] = AsymmFit->GetParameter(0);
    SinAmpErr[7][3] = AsymmFit->GetParError(0);
    CosAmp[7][3] = AsymmFit->GetParameter(1);
    CosAmpErr[7][3] = AsymmFit->GetParError(1);

    PhiSc615AsymmCM8 = Phi_Scattered_615MeV_NegHelCM8->GetAsymmetry(Phi_Scattered_615MeV_PosHelCM8);
    PhiSc615AsymmCM8->SetName("PhiSc615AsymmCM8");
    PhiSc615AsymmCM8->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 580-650MeV (Cos #theta_{CM} -0.75-(-1.0))");
    PhiSc615AsymmCM8->Fit("SinFit");
    InitialSinAmp[7][4] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[7][4] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc615AsymmCM8->Fit("AsymmFit");
    SinAmp[7][4] = AsymmFit->GetParameter(0);
    SinAmpErr[7][4] = AsymmFit->GetParError(0);
    CosAmp[7][4] = AsymmFit->GetParameter(1);
    CosAmpErr[7][4] = AsymmFit->GetParError(1);

    PhiSc685AsymmCM8 = Phi_Scattered_685MeV_NegHelCM8->GetAsymmetry(Phi_Scattered_685MeV_PosHelCM8);
    PhiSc685AsymmCM8->SetName("PhiSc685AsymmCM8");
    PhiSc685AsymmCM8->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 650-720MeV (Cos #theta_{CM} -0.75-(-1.0))");
    PhiSc685AsymmCM8->Fit("SinFit");
    InitialSinAmp[7][5] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[7][5] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc685AsymmCM8->Fit("AsymmFit");
    SinAmp[7][5] = AsymmFit->GetParameter(0);
    SinAmpErr[7][5] = AsymmFit->GetParError(0);
    CosAmp[7][5] = AsymmFit->GetParameter(1);
    CosAmpErr[7][5] = AsymmFit->GetParError(1);

    // Define new file to store fit parameters
    TFile f1("AsymmFits_PTotal59.root", "RECREATE");

    PhiSc335AsymmCM1->Write();
    PhiSc405AsymmCM1->Write();
    PhiSc475AsymmCM1->Write();
    PhiSc545AsymmCM1->Write();
    PhiSc615AsymmCM1->Write();
    PhiSc685AsymmCM1->Write();

    PhiSc335AsymmCM2->Write();
    PhiSc405AsymmCM2->Write();
    PhiSc475AsymmCM2->Write();
    PhiSc545AsymmCM2->Write();
    PhiSc615AsymmCM2->Write();
    PhiSc685AsymmCM2->Write();

    PhiSc335AsymmCM3->Write();
    PhiSc405AsymmCM3->Write();
    PhiSc475AsymmCM3->Write();
    PhiSc545AsymmCM3->Write();
    PhiSc615AsymmCM3->Write();
    PhiSc685AsymmCM3->Write();

    PhiSc335AsymmCM4->Write();
    PhiSc405AsymmCM4->Write();
    PhiSc475AsymmCM4->Write();
    PhiSc545AsymmCM4->Write();
    PhiSc615AsymmCM4->Write();
    PhiSc685AsymmCM4->Write();

    PhiSc335AsymmCM5->Write();
    PhiSc405AsymmCM5->Write();
    PhiSc475AsymmCM5->Write();
    PhiSc545AsymmCM5->Write();
    PhiSc615AsymmCM5->Write();
    PhiSc685AsymmCM5->Write();

    PhiSc335AsymmCM6->Write();
    PhiSc405AsymmCM6->Write();
    PhiSc475AsymmCM6->Write();
    PhiSc545AsymmCM6->Write();
    PhiSc615AsymmCM6->Write();
    PhiSc685AsymmCM6->Write();

    PhiSc335AsymmCM7->Write();
    PhiSc405AsymmCM7->Write();
    PhiSc475AsymmCM7->Write();
    PhiSc545AsymmCM7->Write();
    PhiSc615AsymmCM7->Write();
    PhiSc685AsymmCM7->Write();

    PhiSc335AsymmCM8->Write();
    PhiSc405AsymmCM8->Write();
    PhiSc475AsymmCM8->Write();
    PhiSc545AsymmCM8->Write();
    PhiSc615AsymmCM8->Write();
    PhiSc685AsymmCM8->Write();

    //Define new tree to store parameters in
    TTree* tree = new TTree("Parameter_Values", "Tree_of_Values");

    // Define branches to store parameters, (Branch Name, Variable, Type of Variable)

    tree->Branch("InitialSinAmp335", &ISinAm335, "ISinAm335/D");
    tree->Branch("InitialSinAmpErr335", &ISinAmErr335, "ISinAmErr335/D");
    tree->Branch("CosAmp335", &CosAm335, "CosAm335/D");
    tree->Branch("CosAmpErr335", &CosAmErr335, "CosAmErr335/D");
    tree->Branch("SinAmp335", &SinAm335, "SinAm335/D");
    tree->Branch("SinAmpErr335", &SinAmErr335, "SinAmErr335/D");
    tree->Branch("InitialSinAmp405", &ISinAm405, "ISinAm405/D");
    tree->Branch("InitialSinAmpErr405", &ISinAmErr405, "ISinAmErr405/D");
    tree->Branch("CosAmp405", &CosAm405, "CosAm405/D");
    tree->Branch("CosAmpErr405", &CosAmErr405, "CosAmErr405/D");
    tree->Branch("SinAmp405", &SinAm405, "SinAm405/D");
    tree->Branch("SinAmpErr405", &SinAmErr405, "SinAmErr405/D");
    tree->Branch("InitialSinAmp475", &ISinAm475, "ISinAm475/D");
    tree->Branch("InitialSinAmpErr475", &ISinAmErr475, "ISinAmErr475/D");
    tree->Branch("CosAmp475", &CosAm475, "CosAm475/D");
    tree->Branch("CosAmpErr475", &CosAmErr475, "CosAmErr475/D");
    tree->Branch("SinAmp475", &SinAm475, "SinAm475/D");
    tree->Branch("SinAmpErr475", &SinAmErr475, "SinAmErr475/D");
    tree->Branch("InitialSinAmp545", &ISinAm545, "ISinAm545/D");
    tree->Branch("InitialSinAmpErr545", &ISinAmErr545, "ISinAmErr545/D");
    tree->Branch("CosAmp545", &CosAm545, "CosAm545/D");
    tree->Branch("CosAmpErr545", &CosAmErr545, "CosAmErr545/D");
    tree->Branch("SinAmp545", &SinAm545, "SinAm545/D");
    tree->Branch("SinAmpErr545", &SinAmErr545, "SinAmErr545/D");
    tree->Branch("InitialSinAmp615", &ISinAm615, "ISinAm615/D");
    tree->Branch("InitialSinAmpErr615", &ISinAmErr615, "ISinAmErr615/D");
    tree->Branch("CosAmp615", &CosAm615, "CosAm615/D");
    tree->Branch("CosAmpErr615", &CosAmErr615, "CosAmErr615/D");
    tree->Branch("SinAmp615", &SinAm615, "SinAm615/D");
    tree->Branch("SinAmpErr615", &SinAmErr615, "SinAmErr615/D");
    tree->Branch("InitialSinAmp685", &ISinAm685, "ISinAm685/D");
    tree->Branch("InitialSinAmpErr685", &ISinAmErr685, "ISinAmErr685/D");
    tree->Branch("CosAmp685", &CosAm685, "CosAm685/D");
    tree->Branch("CosAmpErr685", &CosAmErr685, "CosAmErr685/D");
    tree->Branch("SinAmp685", &SinAm685, "SinAm685/D");
    tree->Branch("SinAmpErr685", &SinAmErr685, "SinAmErr685/D");

    for(Int_t m = 0; m < 8; m++){

        ISinAm335 = InitialSinAmp[m][0];
        ISinAmErr335 = InitialSinAmpErr[m][0];
        SinAm335 = SinAmp[m][0];
        SinAmErr335 = SinAmpErr[m][0];
        CosAm335 = CosAmp[m][0];
        CosAmErr335 = CosAmpErr[m][0];
        ISinAm405 = InitialSinAmp[m][1];
        ISinAmErr405 = InitialSinAmpErr[m][1];
        SinAm405 = SinAmp[m][1];
        SinAmErr405 = SinAmpErr[m][1];
        CosAm405 = CosAmp[m][1];
        CosAmErr405 = CosAmpErr[m][1];
        ISinAm475 = InitialSinAmp[m][2];
        ISinAmErr475 = InitialSinAmpErr[m][2];
        SinAm475 = SinAmp[m][2];
        SinAmErr475 = SinAmpErr[m][2];
        CosAm475 = CosAmp[m][2];
        CosAmErr475 = CosAmpErr[m][2];
        ISinAm545 = InitialSinAmp[m][3];
        ISinAmErr545 = InitialSinAmpErr[m][3];
        SinAm545 = SinAmp[m][3];
        SinAmErr545 = SinAmpErr[m][3];
        CosAm545 = CosAmp[m][3];
        CosAmErr545 = CosAmpErr[m][3];
        ISinAm615 = InitialSinAmp[m][4];
        ISinAmErr615 = InitialSinAmpErr[m][4];
        SinAm615 = SinAmp[m][4];
        SinAmErr615 = SinAmpErr[m][4];
        CosAm615 = CosAmp[m][4];
        CosAmErr615 = CosAmpErr[m][4];
        ISinAm685 = InitialSinAmp[m][5];
        ISinAmErr685 = InitialSinAmpErr[m][5];
        SinAm685 = SinAmp[m][5];
        SinAmErr685 = SinAmpErr[m][5];
        CosAm685 = CosAmp[m][5];
        CosAmErr685 = CosAmpErr[m][5];

        tree->Fill();

    }

    f1.Write();
    f1.Close();

    TFile *f2= TFile::Open("/scratch/Mainz_Software/a2GoAT/CircPol_Aug16.root");

    // Calculate values of Cx for each angular and energy bin
    for (Int_t n = 0; n < 8; n++){

        Cx335[n] = SinAmp[n][0]/(0.1*(Graph->Eval(335,0)));
        CxErr335[n] = SinAmpErr[n][0]/((0.1)*(Graph->Eval(335,0)));
        Cx405[n] = SinAmp[n][1]/(0.1*(Graph->Eval(405,0)));
        CxErr405[n] = SinAmpErr[n][1]/((0.1)*(Graph->Eval(405,0)));
        Cx475[n] = SinAmp[n][2]/(0.1*(Graph->Eval(475,0)));
        CxErr475[n] = SinAmpErr[n][2]/((0.1)*(Graph->Eval(475,0)));
        Cx545[n] = SinAmp[n][3]/(0.1*(Graph->Eval(545,0)));
        CxErr545[n] = SinAmpErr[n][3]/((0.1)*(Graph->Eval(545,0)));
        Cx615[n] = SinAmp[n][4]/(0.1*(Graph->Eval(615,0)));
        CxErr615[n] = SinAmpErr[n][4]/((0.1)*(Graph->Eval(615,0)));
        Cx685[n] = SinAmp[n][5]/(0.1*(Graph->Eval(685,0)));
        CxErr685[n] = SinAmpErr[n][5]/((0.1)*(Graph->Eval(685,0)));

    }

    TFile f3("Cx_Plots_59.root", "RECREATE");

    Float_t xMin = -1;
    Float_t xMax = 1;
    Float_t yMin = -5;
    Float_t yMax = 5;
    Double_t x[8] = {0.875, 0.625, 0.375, 0.125, -0.125, -0.375, -0.625, -0.875}; // Need to adjust
    Double_t ex[8] = {0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125}; // Need to adjust

    TCanvas *canvas = new TCanvas("canvas1","canvas1", 1920, 1080);
    TPad *pad1 = new TPad("pad1","",0,0,1,1);
    pad1->Draw();
    pad1->cd();

    pad1->SetTickx(1);
    pad1->SetTicky(1);
    pad1->SetGridx(1);
    pad1->SetGridy(1);
    TH1F  *hr;
    hr1 = canvas1->DrawFrame(xMin, -1,xMax, 1);
    hr1->SetTitle("C_{x} as fn of Cos#theta_{CM} (E_{#gamma} 300-370MeV)");

    gr1 = new TGraphErrors(8, x, Cx335, ex, CxErr335);
    gr1->SetMarkerColor(2);
    gr1->SetMarkerStyle(5);
    gr1->SetMarkerSize(2);
    gr1->SetTitle("C_{x} as fn of Cos#theta_{CM} (E_{#gamma} 300-370MeV)");
    gr1->GetXaxis()->SetTitle("Cos#theta_{CM}");
    gr1->GetYaxis()->SetTitle("C_{x}");
    gr1->SetName("Cx335");
    gr1->Draw("ep");

    TCanvas *canvas2 = new TCanvas("canvas2","canvas2", 1920, 1080);
    TPad *pad2 = new TPad("pad2","",0,0,1,1);
    pad2->Draw();
    pad2->cd();

    pad2->SetTickx(1);
    pad2->SetTicky(1);
    pad2->SetGridx(1);
    pad2->SetGridy(1);
    TH1F  *hr2;
    hr2 = canvas2->DrawFrame(xMin, -1,xMax, 1);
    hr2->SetTitle("C_{x} as fn of Cos#theta_{CM} (E_{#gamma} 370-440MeV)");

    gr2 = new TGraphErrors(8, x, Cx405, ex, CxErr405);
    gr2->SetMarkerColor(2);
    gr2->SetMarkerStyle(5);
    gr2->SetMarkerSize(2);
    gr2->SetTitle("C_{x} as fn of Cos#theta_{CM} (E_{#gamma} 370-440MeV)");
    gr2->GetXaxis()->SetTitle("Cos#theta_{CM}");
    gr2->GetYaxis()->SetTitle("C_{x}");
    gr2->SetName("Cx405");
    gr2->Draw("ep");

    TCanvas *canvas3 = new TCanvas("canvas3","canvas3", 1920, 1080);
    TPad *pad3 = new TPad("pad3","",0,0,1,1);
    pad3->Draw();
    pad3->cd();

    pad3->SetTickx(1);
    pad3->SetTicky(1);
    pad3->SetGridx(1);
    pad3->SetGridy(1);
    TH1F  *hr3;
    hr3 = canvas3->DrawFrame(xMin, -1,xMax, 1);
    hr3->SetTitle("C_{x} as fn of Cos#theta_{CM} (E_{#gamma} 440-510MeV)");

    gr3 = new TGraphErrors(8, x, Cx475, ex, CxErr475);
    gr3->SetMarkerColor(2);
    gr3->SetMarkerStyle(5);
    gr3->SetMarkerSize(2);
    gr3->SetTitle("C_{x} as fn of Cos#theta_{CM} (E_{#gamma} 440-510MeV)");
    gr3->GetXaxis()->SetTitle("Cos#theta_{CM}");
    gr3->GetYaxis()->SetTitle("C_{x}");
    gr3->SetName("Cx475");
    gr3->Draw("ep");

    TCanvas *canvas4 = new TCanvas("canvas4","canvas4", 1920, 1080);
    TPad *pad4 = new TPad("pad4","",0,0,1,1);
    pad4->Draw();
    pad4->cd();

    pad4->SetTickx(1);
    pad4->SetTicky(1);
    pad4->SetGridx(1);
    pad4->SetGridy(1);
    TH1F  *hr4;
    hr4 = canvas4->DrawFrame(xMin, -1,xMax, 1);
    hr4->SetTitle("C_{x} as fn of Cos#theta_{CM} (E_{#gamma} 510-580MeV)");

    gr4 = new TGraphErrors(8, x, Cx545, ex, CxErr545);
    gr4->SetMarkerColor(2);
    gr4->SetMarkerStyle(5);
    gr4->SetMarkerSize(2);
    gr4->SetTitle("C_{x} as fn of Cos#theta_{CM} (E_{#gamma} 510-580MeV)");
    gr4->GetXaxis()->SetTitle("Cos#theta_{CM}");
    gr4->GetYaxis()->SetTitle("C_{x}");
    gr4->SetName("Cx545");
    gr4->Draw("ep");

    TCanvas *canvas5 = new TCanvas("canvas5","canvas5", 1920, 1080);
    TPad *pad5 = new TPad("pad5","",0,0,1,1);
    pad5->Draw();
    pad5->cd();

    pad5->SetTickx(1);
    pad5->SetTicky(1);
    pad5->SetGridx(1);
    pad5->SetGridy(1);
    TH1F  *hr5;
    hr5 = canvas5->DrawFrame(xMin, -1,xMax, 1);
    hr5->SetTitle("C_{x} as fn of Cos#theta_{CM} (E_{#gamma} 580-650MeV)");

    gr5 = new TGraphErrors(8, x, Cx615, ex, CxErr615);
    gr5->SetMarkerColor(2);
    gr5->SetMarkerStyle(5);
    gr5->SetMarkerSize(2);
    gr5->SetTitle("C_{x} as fn of Cos#theta_{CM} (E_{#gamma} 580-650MeV)");
    gr5->GetXaxis()->SetTitle("Cos#theta_{CM}");
    gr5->GetYaxis()->SetTitle("C_{x}");
    gr5->SetName("Cx615");
    gr5->Draw("ep");

    TCanvas *canvas6 = new TCanvas("canvas6","canvas6", 1920, 1080);
    TPad *pad6 = new TPad("pad6","",0,0,1,1);
    pad6->Draw();
    pad6->cd();

    pad6->SetTickx(1);
    pad6->SetTicky(1);
    pad6->SetGridx(1);
    pad6->SetGridy(1);
    TH1F  *hr6;
    hr6 = canvas6->DrawFrame(xMin, -1,xMax, 1);
    hr6->SetTitle("C_{x} as fn of Cos#theta_{CM} (E_{#gamma} 650-720MeV)");

    gr6 = new TGraphErrors(8, x, Cx685, ex, CxErr685);
    gr6->SetMarkerColor(2);
    gr6->SetMarkerStyle(5);
    gr6->SetMarkerSize(2);
    gr6->SetTitle("C_{x} as fn of Cos#theta_{CM} (E_{#gamma} 650-720MeV)");
    gr6->GetXaxis()->SetTitle("Cos#theta_{CM}");
    gr6->GetYaxis()->SetTitle("C_{x}");
    gr6->SetName("Cx685");
    gr6->Draw("ep");

    canvas1->Write();
    gr1->Write();
    canvas2->Write();
    gr2->Write();
    canvas3->Write();;
    gr3->Write();
    canvas4->Write();
    gr4->Write();
    canvas5->Write();
    gr5->Write();
    canvas6->Write();
    gr6->Write();

    TCanvas *canvas7 = new TCanvas("canvas7","canvas7", 1920, 1080);
    canvas7->Divide(3,2);
    canvas7->cd(1);
    pad1->Draw();
    canvas7->cd(2);
    pad2->Draw();
    canvas7->cd(3);
    pad3->Draw();
    canvas7->cd(4);
    pad4->Draw();
    canvas7->cd(5);
    pad5->Draw();
    canvas7->cd(6);
    pad6->Draw();
    canvas7->Write();

    f3.Write();

}
