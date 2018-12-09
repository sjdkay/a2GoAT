#include "./includes.h"

// define a function with 4 parameters
Double_t fitf(Double_t *x,Double_t *par)
{
    Double_t fitval = 0;
    fitval =  ((par[0]*sin(x[0]))/(1 + (par[1]*cos(x[0]))));
    return fitval;
}

void CxAsymm() {

    double InitialSinAmp[8][8];
    double InitialSinAmpErr[8][8];
    double SinAmp[8][8]; // Format of arrays is Theta Bin Defined by i, Energy bin defined by j
    double SinAmpErr[8][8];
    double CosAmp[8][8]; // Format of arrays is Theta Bin Defined by i, Energy bin defined by j
    double CosAmpErr[8][8];
    Int_t i;
    double ISinAm265;
    double ISinAmErr265;
    double SinAm265;
    double SinAmErr265;
    double CosAm265;
    double CosAmErr265;
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

    Double_t Cx265[8], Cx335[8], Cx405[8], Cx475[8], Cx545[8], Cx615[8], Cx685[8];
    Double_t CxErr265[8], CxErr335[8], CxErr405[8], CxErr475[8], CxErr545[8], CxErr615[8], CxErr685[8];
    Double_t Pn265[8], Pn335[8], Pn405[8], Pn475[8], Pn545[8], Pn615[8], Pn685[8];
    Double_t PnErr265[8], PnErr335[8], PnErr405[8], PnErr475[8], PnErr545[8], PnErr615[8], PnErr685[8];

    TF1 *AsymmFunc = new TF1("AsymmFit",  fitf, -3.0, 3.0, 2); //Give a name and range to the fitting funcion
    AsymmFunc->SetParNames("SinAmp", "CosAmp"); //Name the parameters
    AsymmFunc->SetParameter(0, 0);
    TF1 *SinFunc = new TF1("SinFit", "[0]*sin(x*TMath::DegToRad())", -3, 3);
    SinFunc->SetParNames("InitialSinAmp");
    TFile *f = new TFile("/scratch/Mainz_Software/a2GoAT/Physics_Total_105_5_4_18.root"); // Open the latest PTotal file to load histograms from
    TF1 *Pn90CM = new TF1("Pn90CM", "1.64576-2.95484*(x/1000)+0.684577*(x/1000)**2-0.65*90**2/4/((x-560)**2+90**2/4)+(5.32305-35.3819*(x/1000)+70.145*(x/1000)**2-44.2899*(x/1000)**3)",300,700);

    ///////////////////////////////////////////
    //////////////////  CM1  //////////////////
    ///////////////////////////////////////////

    PhiSc265AsymmCM1 = Phi_Scattered_265MeV_NegHelCM1->GetAsymmetry(Phi_Scattered_265MeV_PosHelCM1);
    PhiSc265AsymmCM1->SetName("PhiSc265AsymmCM1");
    PhiSc265AsymmCM1->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 230-300MeV (Cos #theta_{CM} 1-0.75)");
    PhiSc265AsymmCM1->Fit("SinFit", "Q");
    InitialSinAmp[0][0] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[0][0] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc265AsymmCM1->Fit("AsymmFit", "QLL");
    SinAmp[0][0] = AsymmFit->GetParameter(0);
    SinAmpErr[0][0] = AsymmFit->GetParError(0);
    CosAmp[0][0] = AsymmFit->GetParameter(1);
    CosAmpErr[0][0] = AsymmFit->GetParError(1);

    PhiSc335AsymmCM1 = Phi_Scattered_335MeV_NegHelCM1->GetAsymmetry(Phi_Scattered_335MeV_PosHelCM1);
    PhiSc335AsymmCM1->SetName("PhiSc335AsymmCM1");
    PhiSc335AsymmCM1->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 300-370MeV (Cos #theta_{CM} 1-0.75)");
    PhiSc335AsymmCM1->Fit("SinFit", "Q");
    InitialSinAmp[0][1] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[0][1] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc335AsymmCM1->Fit("AsymmFit", "QLL");
    SinAmp[0][1] = AsymmFit->GetParameter(0);
    SinAmpErr[0][1] = AsymmFit->GetParError(0);
    CosAmp[0][1] = AsymmFit->GetParameter(1);
    CosAmpErr[0][1] = AsymmFit->GetParError(1);

    PhiSc405AsymmCM1 = Phi_Scattered_405MeV_NegHelCM1->GetAsymmetry(Phi_Scattered_405MeV_PosHelCM1);
    PhiSc405AsymmCM1->SetName("PhiSc405AsymmCM1");
    PhiSc405AsymmCM1->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 370-440MeV (Cos #theta_{CM} 1-0.75)");
    PhiSc405AsymmCM1->Fit("SinFit", "Q");
    InitialSinAmp[0][2] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[0][2] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc405AsymmCM1->Fit("AsymmFit", "QLL");
    SinAmp[0][2] = AsymmFit->GetParameter(0);
    SinAmpErr[0][2] = AsymmFit->GetParError(0);
    CosAmp[0][2] = AsymmFit->GetParameter(1);
    CosAmpErr[0][2] = AsymmFit->GetParError(1);

    PhiSc475AsymmCM1 = Phi_Scattered_475MeV_NegHelCM1->GetAsymmetry(Phi_Scattered_475MeV_PosHelCM1);
    PhiSc475AsymmCM1->SetName("PhiSc475AsymmCM1");
    PhiSc475AsymmCM1->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 440-510MeV (Cos #theta_{CM} 1-0.75)");
    PhiSc475AsymmCM1->Fit("SinFit", "Q");
    InitialSinAmp[0][3] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[0][3] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc475AsymmCM1->Fit("AsymmFit", "QLL");
    SinAmp[0][3] = AsymmFit->GetParameter(0);
    SinAmpErr[0][3] = AsymmFit->GetParError(0);
    CosAmp[0][3] = AsymmFit->GetParameter(1);
    CosAmpErr[0][3] = AsymmFit->GetParError(1);

    PhiSc545AsymmCM1 = Phi_Scattered_545MeV_NegHelCM1->GetAsymmetry(Phi_Scattered_545MeV_PosHelCM1);
    PhiSc545AsymmCM1->SetName("PhiSc545AsymmCM1");
    PhiSc545AsymmCM1->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 510-580MeV (Cos #theta_{CM} 1-0.75)");
    PhiSc545AsymmCM1->Fit("SinFit", "Q");
    InitialSinAmp[0][4] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[0][4] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc545AsymmCM1->Fit("AsymmFit", "QLL");
    SinAmp[0][4] = AsymmFit->GetParameter(0);
    SinAmpErr[0][4] = AsymmFit->GetParError(0);
    CosAmp[0][4] = AsymmFit->GetParameter(1);
    CosAmpErr[0][4] = AsymmFit->GetParError(1);

    PhiSc615AsymmCM1 = Phi_Scattered_615MeV_NegHelCM1->GetAsymmetry(Phi_Scattered_615MeV_PosHelCM1);
    PhiSc615AsymmCM1->SetName("PhiSc615AsymmCM1");
    PhiSc615AsymmCM1->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 580-650MeV (Cos #theta_{CM} 1-0.75)");
    PhiSc615AsymmCM1->Fit("SinFit", "Q");
    InitialSinAmp[0][5] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[0][5] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc615AsymmCM1->Fit("AsymmFit", "QLL");
    SinAmp[0][5] = AsymmFit->GetParameter(0);
    SinAmpErr[0][5] = AsymmFit->GetParError(0);
    CosAmp[0][5] = AsymmFit->GetParameter(1);
    CosAmpErr[0][5] = AsymmFit->GetParError(1);

    PhiSc685AsymmCM1 = Phi_Scattered_685MeV_NegHelCM1->GetAsymmetry(Phi_Scattered_685MeV_PosHelCM1);
    PhiSc685AsymmCM1->SetName("PhiSc685AsymmCM1");
    PhiSc685AsymmCM1->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 650-720MeV (Cos #theta_{CM} 1-0.75)");
    PhiSc685AsymmCM1->Fit("SinFit", "Q");
    InitialSinAmp[0][6] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[0][6] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc685AsymmCM1->Fit("AsymmFit", "QLL");
    SinAmp[0][6] = AsymmFit->GetParameter(0);
    SinAmpErr[0][6] = AsymmFit->GetParError(0);
    CosAmp[0][6] = AsymmFit->GetParameter(1);
    CosAmpErr[0][6] = AsymmFit->GetParError(1);

    ///////////////////////////////////////////
    //////////////////  CM2  //////////////////
    ///////////////////////////////////////////

    PhiSc265AsymmCM2 = Phi_Scattered_265MeV_NegHelCM2->GetAsymmetry(Phi_Scattered_265MeV_PosHelCM2);
    PhiSc265AsymmCM2->SetName("PhiSc265AsymmCM2");
    PhiSc265AsymmCM2->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 230-300MeV (Cos #theta_{CM} 0.75-0.5)");
    PhiSc265AsymmCM2->Fit("SinFit", "Q");
    InitialSinAmp[1][0] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[1][0] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc265AsymmCM2->Fit("AsymmFit", "QLL");
    SinAmp[1][0] = AsymmFit->GetParameter(0);
    SinAmpErr[1][0] = AsymmFit->GetParError(0);
    CosAmp[1][0] = AsymmFit->GetParameter(1);
    CosAmpErr[1][0] = AsymmFit->GetParError(1);

    PhiSc335AsymmCM2 = Phi_Scattered_335MeV_NegHelCM2->GetAsymmetry(Phi_Scattered_335MeV_PosHelCM2);
    PhiSc335AsymmCM2->SetName("PhiSc335AsymmCM2");
    PhiSc335AsymmCM2->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 300-370MeV (Cos #theta_{CM} 0.75-0.5)");
    PhiSc335AsymmCM2->Fit("SinFit", "Q");
    InitialSinAmp[1][1] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[1][1] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc335AsymmCM2->Fit("AsymmFit", "QLL");
    SinAmp[1][1] = AsymmFit->GetParameter(0);
    SinAmpErr[1][1] = AsymmFit->GetParError(0);
    CosAmp[1][1] = AsymmFit->GetParameter(1);
    CosAmpErr[1][1] = AsymmFit->GetParError(1);

    PhiSc405AsymmCM2 = Phi_Scattered_405MeV_NegHelCM2->GetAsymmetry(Phi_Scattered_405MeV_PosHelCM2);
    PhiSc405AsymmCM2->SetName("PhiSc405AsymmCM2");
    PhiSc405AsymmCM2->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 370-440MeV (Cos #theta_{CM} 0.75-0.5)");
    PhiSc405AsymmCM2->Fit("SinFit", "Q");
    InitialSinAmp[1][2] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[1][2] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc405AsymmCM2->Fit("AsymmFit", "QLL");
    SinAmp[1][2] = AsymmFit->GetParameter(0);
    SinAmpErr[1][2] = AsymmFit->GetParError(0);
    CosAmp[1][2] = AsymmFit->GetParameter(1);
    CosAmpErr[1][2] = AsymmFit->GetParError(1);

    PhiSc475AsymmCM2 = Phi_Scattered_475MeV_NegHelCM2->GetAsymmetry(Phi_Scattered_475MeV_PosHelCM2);
    PhiSc475AsymmCM2->SetName("PhiSc475AsymmCM2");
    PhiSc475AsymmCM2->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 440-510MeV (Cos #theta_{CM} 0.75-0.5)");
    PhiSc475AsymmCM2->Fit("SinFit", "Q");
    InitialSinAmp[1][3] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[1][3] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc475AsymmCM2->Fit("AsymmFit", "QLL");
    SinAmp[1][3] = AsymmFit->GetParameter(0);
    SinAmpErr[1][3] = AsymmFit->GetParError(0);
    CosAmp[1][3] = AsymmFit->GetParameter(1);
    CosAmpErr[1][3] = AsymmFit->GetParError(1);

    PhiSc545AsymmCM2 = Phi_Scattered_545MeV_NegHelCM2->GetAsymmetry(Phi_Scattered_545MeV_PosHelCM2);
    PhiSc545AsymmCM2->SetName("PhiSc545AsymmCM2");
    PhiSc545AsymmCM2->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 510-580MeV (Cos #theta_{CM} 0.75-0.5)");
    PhiSc545AsymmCM2->Fit("SinFit", "Q");
    InitialSinAmp[1][4] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[1][4] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc545AsymmCM2->Fit("AsymmFit", "QLL");
    SinAmp[1][4] = AsymmFit->GetParameter(0);
    SinAmpErr[1][4] = AsymmFit->GetParError(0);
    CosAmp[1][4] = AsymmFit->GetParameter(1);
    CosAmpErr[1][4] = AsymmFit->GetParError(1);

    PhiSc615AsymmCM2 = Phi_Scattered_615MeV_NegHelCM2->GetAsymmetry(Phi_Scattered_615MeV_PosHelCM2);
    PhiSc615AsymmCM2->SetName("PhiSc615AsymmCM2");
    PhiSc615AsymmCM2->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 580-650MeV (Cos #theta_{CM} 0.75-0.5)");
    PhiSc615AsymmCM2->Fit("SinFit", "Q");
    InitialSinAmp[1][5] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[1][5] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc615AsymmCM2->Fit("AsymmFit", "QLL");
    SinAmp[1][5] = AsymmFit->GetParameter(0);
    SinAmpErr[1][5] = AsymmFit->GetParError(0);
    CosAmp[1][5] = AsymmFit->GetParameter(1);
    CosAmpErr[1][5] = AsymmFit->GetParError(1);

    PhiSc685AsymmCM2 = Phi_Scattered_685MeV_NegHelCM2->GetAsymmetry(Phi_Scattered_685MeV_PosHelCM2);
    PhiSc685AsymmCM2->SetName("PhiSc685AsymmCM2");
    PhiSc685AsymmCM2->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 650-720MeV (Cos #theta_{CM} 0.75-0.5)");
    PhiSc685AsymmCM2->Fit("SinFit", "Q");
    InitialSinAmp[1][6] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[1][6] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc685AsymmCM2->Fit("AsymmFit", "QLL");
    SinAmp[1][6] = AsymmFit->GetParameter(0);
    SinAmpErr[1][6] = AsymmFit->GetParError(0);
    CosAmp[1][6] = AsymmFit->GetParameter(1);
    CosAmpErr[1][6] = AsymmFit->GetParError(1);

    ///////////////////////////////////////////
    //////////////////  CM3  //////////////////
    ///////////////////////////////////////////

    PhiSc265AsymmCM3 = Phi_Scattered_265MeV_NegHelCM3->GetAsymmetry(Phi_Scattered_265MeV_PosHelCM3);
    PhiSc265AsymmCM3->SetName("PhiSc265AsymmCM3");
    PhiSc265AsymmCM3->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 230-300MeV (Cos #theta_{CM} 0.5-0.25)");
    PhiSc265AsymmCM3->Fit("SinFit", "Q");
    InitialSinAmp[2][0] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[2][0] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc265AsymmCM3->Fit("AsymmFit", "QLL");
    SinAmp[2][0] = AsymmFit->GetParameter(0);
    SinAmpErr[2][0] = AsymmFit->GetParError(0);
    CosAmp[2][0] = AsymmFit->GetParameter(1);
    CosAmpErr[2][0] = AsymmFit->GetParError(1);

    PhiSc335AsymmCM3 = Phi_Scattered_335MeV_NegHelCM3->GetAsymmetry(Phi_Scattered_335MeV_PosHelCM3);
    PhiSc335AsymmCM3->SetName("PhiSc335AsymmCM3");
    PhiSc335AsymmCM3->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 300-370MeV (Cos #theta_{CM} 0.5-0.25)");
    PhiSc335AsymmCM3->Fit("SinFit", "Q");
    InitialSinAmp[2][1] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[2][1] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc335AsymmCM3->Fit("AsymmFit", "QLL");
    SinAmp[2][1] = AsymmFit->GetParameter(0);
    SinAmpErr[2][1] = AsymmFit->GetParError(0);
    CosAmp[2][1] = AsymmFit->GetParameter(1);
    CosAmpErr[2][1] = AsymmFit->GetParError(1);

    PhiSc405AsymmCM3 = Phi_Scattered_405MeV_NegHelCM3->GetAsymmetry(Phi_Scattered_405MeV_PosHelCM3);
    PhiSc405AsymmCM3->SetName("PhiSc405AsymmCM3");
    PhiSc405AsymmCM3->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 370-440MeV (Cos #theta_{CM} 0.5-0.25)");
    PhiSc405AsymmCM3->Fit("SinFit", "Q");
    InitialSinAmp[2][2] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[2][2] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc405AsymmCM3->Fit("AsymmFit", "QLL");
    SinAmp[2][2] = AsymmFit->GetParameter(0);
    SinAmpErr[2][2] = AsymmFit->GetParError(0);
    CosAmp[2][2] = AsymmFit->GetParameter(1);
    CosAmpErr[2][2] = AsymmFit->GetParError(1);

    PhiSc475AsymmCM3 = Phi_Scattered_475MeV_NegHelCM3->GetAsymmetry(Phi_Scattered_475MeV_PosHelCM3);
    PhiSc475AsymmCM3->SetName("PhiSc475AsymmCM3");
    PhiSc475AsymmCM3->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 440-510MeV (Cos #theta_{CM} 0.5-0.25)");
    PhiSc475AsymmCM3->Fit("SinFit", "Q");
    InitialSinAmp[2][3] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[2][3] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc475AsymmCM3->Fit("AsymmFit", "QLL");
    SinAmp[2][3] = AsymmFit->GetParameter(0);
    SinAmpErr[2][3] = AsymmFit->GetParError(0);
    CosAmp[2][3] = AsymmFit->GetParameter(1);
    CosAmpErr[2][3] = AsymmFit->GetParError(1);

    PhiSc545AsymmCM3 = Phi_Scattered_545MeV_NegHelCM3->GetAsymmetry(Phi_Scattered_545MeV_PosHelCM3);
    PhiSc545AsymmCM3->SetName("PhiSc545AsymmCM3");
    PhiSc545AsymmCM3->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 510-580MeV (Cos #theta_{CM} 0.5-0.25)");
    PhiSc545AsymmCM3->Fit("SinFit", "Q");
    InitialSinAmp[2][4] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[2][4] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc545AsymmCM3->Fit("AsymmFit", "QLL");
    SinAmp[2][4] = AsymmFit->GetParameter(0);
    SinAmpErr[2][4] = AsymmFit->GetParError(0);
    CosAmp[2][4] = AsymmFit->GetParameter(1);
    CosAmpErr[2][4] = AsymmFit->GetParError(1);

    PhiSc615AsymmCM3 = Phi_Scattered_615MeV_NegHelCM3->GetAsymmetry(Phi_Scattered_615MeV_PosHelCM3);
    PhiSc615AsymmCM3->SetName("PhiSc615AsymmCM3");
    PhiSc615AsymmCM3->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 580-650MeV (Cos #theta_{CM} 0.5-0.25)");
    PhiSc615AsymmCM3->Fit("SinFit", "Q");
    InitialSinAmp[2][5] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[2][5] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc615AsymmCM3->Fit("AsymmFit", "QLL");
    SinAmp[2][5] = AsymmFit->GetParameter(0);
    SinAmpErr[2][5] = AsymmFit->GetParError(0);
    CosAmp[2][5] = AsymmFit->GetParameter(1);
    CosAmpErr[2][5] = AsymmFit->GetParError(1);

    PhiSc685AsymmCM3 = Phi_Scattered_685MeV_NegHelCM3->GetAsymmetry(Phi_Scattered_685MeV_PosHelCM3);
    PhiSc685AsymmCM3->SetName("PhiSc685AsymmCM3");
    PhiSc685AsymmCM3->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 650-720MeV (Cos #theta_{CM} 0.5-0.25)");
    PhiSc685AsymmCM3->Fit("SinFit", "Q");
    InitialSinAmp[2][6] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[2][6] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc685AsymmCM3->Fit("AsymmFit", "QLL");
    SinAmp[2][6] = AsymmFit->GetParameter(0);
    SinAmpErr[2][6] = AsymmFit->GetParError(0);
    CosAmp[2][6] = AsymmFit->GetParameter(1);
    CosAmpErr[2][6] = AsymmFit->GetParError(1);

    ///////////////////////////////////////////
    //////////////////  CM4  //////////////////
    ///////////////////////////////////////////

    PhiSc265AsymmCM4 = Phi_Scattered_265MeV_NegHelCM4->GetAsymmetry(Phi_Scattered_265MeV_PosHelCM4);
    PhiSc265AsymmCM4->SetName("PhiSc265AsymmCM4");
    PhiSc265AsymmCM4->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 230-300MeV (Cos #theta_{CM} 0.25-0.0)");
    PhiSc265AsymmCM4->Fit("SinFit", "Q");
    InitialSinAmp[3][0] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[3][0] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc265AsymmCM4->Fit("AsymmFit", "QLL");
    SinAmp[3][0] = AsymmFit->GetParameter(0);
    SinAmpErr[3][0] = AsymmFit->GetParError(0);
    CosAmp[3][0] = AsymmFit->GetParameter(1);
    CosAmpErr[3][0] = AsymmFit->GetParError(1);

    PhiSc335AsymmCM4 = Phi_Scattered_335MeV_NegHelCM4->GetAsymmetry(Phi_Scattered_335MeV_PosHelCM4);
    PhiSc335AsymmCM4->SetName("PhiSc335AsymmCM4");
    PhiSc335AsymmCM4->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 300-370MeV (Cos #theta_{CM} 0.25-0.0)");
    PhiSc335AsymmCM4->Fit("SinFit", "Q");
    InitialSinAmp[3][1] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[3][1] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc335AsymmCM4->Fit("AsymmFit", "QLL");
    SinAmp[3][1] = AsymmFit->GetParameter(0);
    SinAmpErr[3][1] = AsymmFit->GetParError(0);
    CosAmp[3][1] = AsymmFit->GetParameter(1);
    CosAmpErr[3][1] = AsymmFit->GetParError(1);

    PhiSc405AsymmCM4 = Phi_Scattered_405MeV_NegHelCM4->GetAsymmetry(Phi_Scattered_405MeV_PosHelCM4);
    PhiSc405AsymmCM4->SetName("PhiSc405AsymmCM4");
    PhiSc405AsymmCM4->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 370-440MeV (Cos #theta_{CM} 0.25-0.0)");
    PhiSc405AsymmCM4->Fit("SinFit", "Q");
    InitialSinAmp[3][2] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[3][2] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc405AsymmCM4->Fit("AsymmFit", "QLL");
    SinAmp[3][2] = AsymmFit->GetParameter(0);
    SinAmpErr[3][2] = AsymmFit->GetParError(0);
    CosAmp[3][2] = AsymmFit->GetParameter(1);
    CosAmpErr[3][2] = AsymmFit->GetParError(1);

    PhiSc475AsymmCM4 = Phi_Scattered_475MeV_NegHelCM4->GetAsymmetry(Phi_Scattered_475MeV_PosHelCM4);
    PhiSc475AsymmCM4->SetName("PhiSc475AsymmCM4");
    PhiSc475AsymmCM4->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 440-510MeV (Cos #theta_{CM} 0.25-0.0)");
    PhiSc475AsymmCM4->Fit("SinFit", "Q");
    InitialSinAmp[3][3] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[3][3] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc475AsymmCM4->Fit("AsymmFit", "QLL");
    SinAmp[3][3] = AsymmFit->GetParameter(0);
    SinAmpErr[3][3] = AsymmFit->GetParError(0);
    CosAmp[3][3] = AsymmFit->GetParameter(1);
    CosAmpErr[3][3] = AsymmFit->GetParError(1);

    PhiSc545AsymmCM4 = Phi_Scattered_545MeV_NegHelCM4->GetAsymmetry(Phi_Scattered_545MeV_PosHelCM4);
    PhiSc545AsymmCM4->SetName("PhiSc545AsymmCM4");
    PhiSc545AsymmCM4->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 510-580MeV (Cos #theta_{CM} 0.25-0.0)");
    PhiSc545AsymmCM4->Fit("SinFit", "Q");
    InitialSinAmp[3][4] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[3][4] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc545AsymmCM4->Fit("AsymmFit", "QLL");
    SinAmp[3][4] = AsymmFit->GetParameter(0);
    SinAmpErr[3][4] = AsymmFit->GetParError(0);
    CosAmp[3][4] = AsymmFit->GetParameter(1);
    CosAmpErr[3][4] = AsymmFit->GetParError(1);

    PhiSc615AsymmCM4 = Phi_Scattered_615MeV_NegHelCM4->GetAsymmetry(Phi_Scattered_615MeV_PosHelCM4);
    PhiSc615AsymmCM4->SetName("PhiSc615AsymmCM4");
    PhiSc615AsymmCM4->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 580-650MeV (Cos #theta_{CM} 0.25-0.0)");
    PhiSc615AsymmCM4->Fit("SinFit", "Q");
    InitialSinAmp[3][5] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[3][5] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc615AsymmCM4->Fit("AsymmFit", "QLL");
    SinAmp[3][5] = AsymmFit->GetParameter(0);
    SinAmpErr[3][5] = AsymmFit->GetParError(0);
    CosAmp[3][5] = AsymmFit->GetParameter(1);
    CosAmpErr[3][5] = AsymmFit->GetParError(1);

    PhiSc685AsymmCM4 = Phi_Scattered_685MeV_NegHelCM4->GetAsymmetry(Phi_Scattered_685MeV_PosHelCM4);
    PhiSc685AsymmCM4->SetName("PhiSc685AsymmCM4");
    PhiSc685AsymmCM4->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 650-720MeV (Cos #theta_{CM} 0.25-0.0)");
    PhiSc685AsymmCM4->Fit("SinFit", "Q");
    InitialSinAmp[3][6] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[3][6] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc685AsymmCM4->Fit("AsymmFit", "QLL");
    SinAmp[3][6] = AsymmFit->GetParameter(0);
    SinAmpErr[3][6] = AsymmFit->GetParError(0);
    CosAmp[3][6] = AsymmFit->GetParameter(1);
    CosAmpErr[3][6] = AsymmFit->GetParError(1);

    ///////////////////////////////////////////
    //////////////////  CM5  //////////////////
    ///////////////////////////////////////////

    PhiSc265AsymmCM5 = Phi_Scattered_265MeV_NegHelCM5->GetAsymmetry(Phi_Scattered_265MeV_PosHelCM5);
    PhiSc265AsymmCM5->SetName("PhiSc265AsymmCM5");
    PhiSc265AsymmCM5->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 230-300MeV (Cos #theta_{CM} 0.0-(-0.25))");
    PhiSc265AsymmCM5->Fit("SinFit", "Q");
    InitialSinAmp[4][0] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[4][0] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc265AsymmCM5->Fit("AsymmFit", "QLL");
    SinAmp[4][0] = AsymmFit->GetParameter(0);
    SinAmpErr[4][0] = AsymmFit->GetParError(0);
    CosAmp[4][0] = AsymmFit->GetParameter(1);
    CosAmpErr[4][0] = AsymmFit->GetParError(1);

    PhiSc335AsymmCM5 = Phi_Scattered_335MeV_NegHelCM5->GetAsymmetry(Phi_Scattered_335MeV_PosHelCM5);
    PhiSc335AsymmCM5->SetName("PhiSc335AsymmCM5");
    PhiSc335AsymmCM5->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 300-370MeV (Cos #theta_{CM} 0.0-(-0.25))");
    PhiSc335AsymmCM5->Fit("SinFit", "Q");
    InitialSinAmp[4][1] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[4][1] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc335AsymmCM5->Fit("AsymmFit", "QLL");
    SinAmp[4][1] = AsymmFit->GetParameter(0);
    SinAmpErr[4][1] = AsymmFit->GetParError(0);
    CosAmp[4][1] = AsymmFit->GetParameter(1);
    CosAmpErr[4][1] = AsymmFit->GetParError(1);

    PhiSc405AsymmCM5 = Phi_Scattered_405MeV_NegHelCM5->GetAsymmetry(Phi_Scattered_405MeV_PosHelCM5);
    PhiSc405AsymmCM5->SetName("PhiSc405AsymmCM5");
    PhiSc405AsymmCM5->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 370-440MeV (Cos #theta_{CM} 0.0-(-0.25))");
    PhiSc405AsymmCM5->Fit("SinFit", "Q");
    InitialSinAmp[4][2] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[4][2] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc405AsymmCM5->Fit("AsymmFit", "QLL");
    SinAmp[4][2] = AsymmFit->GetParameter(0);
    SinAmpErr[4][2] = AsymmFit->GetParError(0);
    CosAmp[4][2] = AsymmFit->GetParameter(1);
    CosAmpErr[4][2] = AsymmFit->GetParError(1);

    PhiSc475AsymmCM5 = Phi_Scattered_475MeV_NegHelCM5->GetAsymmetry(Phi_Scattered_475MeV_PosHelCM5);
    PhiSc475AsymmCM5->SetName("PhiSc475AsymmCM5");
    PhiSc475AsymmCM5->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 440-510MeV (Cos #theta_{CM} 0.0-(-0.25))");
    PhiSc475AsymmCM5->Fit("SinFit", "Q");
    InitialSinAmp[4][3] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[4][3] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc475AsymmCM5->Fit("AsymmFit", "QLL");
    SinAmp[4][3] = AsymmFit->GetParameter(0);
    SinAmpErr[4][3] = AsymmFit->GetParError(0);
    CosAmp[4][3] = AsymmFit->GetParameter(1);
    CosAmpErr[4][3] = AsymmFit->GetParError(1);

    PhiSc545AsymmCM5 = Phi_Scattered_545MeV_NegHelCM5->GetAsymmetry(Phi_Scattered_545MeV_PosHelCM5);
    PhiSc545AsymmCM5->SetName("PhiSc545AsymmCM5");
    PhiSc545AsymmCM5->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 510-580MeV (Cos #theta_{CM} 0.0-(-0.25))");
    PhiSc545AsymmCM5->Fit("SinFit", "Q");
    InitialSinAmp[4][4] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[4][4] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc545AsymmCM5->Fit("AsymmFit", "QLL");
    SinAmp[4][4] = AsymmFit->GetParameter(0);
    SinAmpErr[4][4] = AsymmFit->GetParError(0);
    CosAmp[4][4] = AsymmFit->GetParameter(1);
    CosAmpErr[4][4] = AsymmFit->GetParError(1);

    PhiSc615AsymmCM5 = Phi_Scattered_615MeV_NegHelCM5->GetAsymmetry(Phi_Scattered_615MeV_PosHelCM5);
    PhiSc615AsymmCM5->SetName("PhiSc615AsymmCM5");
    PhiSc615AsymmCM5->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 580-650MeV (Cos #theta_{CM} 0.0-(-0.25))");
    PhiSc615AsymmCM5->Fit("SinFit", "Q");
    InitialSinAmp[4][5] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[4][5] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc615AsymmCM5->Fit("AsymmFit", "QLL");
    SinAmp[4][5] = AsymmFit->GetParameter(0);
    SinAmpErr[4][5] = AsymmFit->GetParError(0);
    CosAmp[4][5] = AsymmFit->GetParameter(1);
    CosAmpErr[4][5] = AsymmFit->GetParError(1);

    PhiSc685AsymmCM5 = Phi_Scattered_685MeV_NegHelCM5->GetAsymmetry(Phi_Scattered_685MeV_PosHelCM5);
    PhiSc685AsymmCM5->SetName("PhiSc685AsymmCM5");
    PhiSc685AsymmCM5->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 650-720MeV (Cos #theta_{CM} 0.0-(-0.25))");
    PhiSc685AsymmCM5->Fit("SinFit", "Q");
    InitialSinAmp[4][6] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[4][6] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc685AsymmCM5->Fit("AsymmFit", "QLL");
    SinAmp[4][6] = AsymmFit->GetParameter(0);
    SinAmpErr[4][6] = AsymmFit->GetParError(0);
    CosAmp[4][6] = AsymmFit->GetParameter(1);
    CosAmpErr[4][6] = AsymmFit->GetParError(1);

    ///////////////////////////////////////////
    //////////////////  CM6  //////////////////
    ///////////////////////////////////////////

    PhiSc265AsymmCM6 = Phi_Scattered_265MeV_NegHelCM6->GetAsymmetry(Phi_Scattered_265MeV_PosHelCM6);
    PhiSc265AsymmCM6->SetName("PhiSc265AsymmCM6");
    PhiSc265AsymmCM6->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 230-300MeV (Cos #theta_{CM} -0.25-(-0.5))");
    PhiSc265AsymmCM6->Fit("SinFit", "Q");
    InitialSinAmp[5][0] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[5][0] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc265AsymmCM6->Fit("AsymmFit", "QLL");
    SinAmp[5][0] = AsymmFit->GetParameter(0);
    SinAmpErr[5][0] = AsymmFit->GetParError(0);
    CosAmp[5][0] = AsymmFit->GetParameter(1);
    CosAmpErr[5][0] = AsymmFit->GetParError(1);

    PhiSc335AsymmCM6 = Phi_Scattered_335MeV_NegHelCM6->GetAsymmetry(Phi_Scattered_335MeV_PosHelCM6);
    PhiSc335AsymmCM6->SetName("PhiSc335AsymmCM6");
    PhiSc335AsymmCM6->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 300-370MeV (Cos #theta_{CM} -0.25-(-0.5))");
    PhiSc335AsymmCM6->Fit("SinFit", "Q");
    InitialSinAmp[5][1] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[5][1] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc335AsymmCM6->Fit("AsymmFit", "QLL");
    SinAmp[5][1] = AsymmFit->GetParameter(0);
    SinAmpErr[5][1] = AsymmFit->GetParError(0);
    CosAmp[5][1] = AsymmFit->GetParameter(1);
    CosAmpErr[5][1] = AsymmFit->GetParError(1);

    PhiSc405AsymmCM6 = Phi_Scattered_405MeV_NegHelCM6->GetAsymmetry(Phi_Scattered_405MeV_PosHelCM6);
    PhiSc405AsymmCM6->SetName("PhiSc405AsymmCM6");
    PhiSc405AsymmCM6->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 370-440MeV (Cos #theta_{CM} -0.25-(-0.5))");
    PhiSc405AsymmCM6->Fit("SinFit", "Q");
    InitialSinAmp[5][2] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[5][2] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc405AsymmCM6->Fit("AsymmFit", "QLL");
    SinAmp[5][2] = AsymmFit->GetParameter(0);
    SinAmpErr[5][2] = AsymmFit->GetParError(0);
    CosAmp[5][2] = AsymmFit->GetParameter(1);
    CosAmpErr[5][2] = AsymmFit->GetParError(1);

    PhiSc475AsymmCM6 = Phi_Scattered_475MeV_NegHelCM6->GetAsymmetry(Phi_Scattered_475MeV_PosHelCM6);
    PhiSc475AsymmCM6->SetName("PhiSc475AsymmCM6");
    PhiSc475AsymmCM6->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 440-510MeV (Cos #theta_{CM} -0.25-(-0.5))");
    PhiSc475AsymmCM6->Fit("SinFit", "Q");
    InitialSinAmp[5][3] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[5][3] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc475AsymmCM6->Fit("AsymmFit", "QLL");
    SinAmp[5][3] = AsymmFit->GetParameter(0);
    SinAmpErr[5][3] = AsymmFit->GetParError(0);
    CosAmp[5][3] = AsymmFit->GetParameter(1);
    CosAmpErr[5][3] = AsymmFit->GetParError(1);

    PhiSc545AsymmCM6 = Phi_Scattered_545MeV_NegHelCM6->GetAsymmetry(Phi_Scattered_545MeV_PosHelCM6);
    PhiSc545AsymmCM6->SetName("PhiSc545AsymmCM6");
    PhiSc545AsymmCM6->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 510-580MeV (Cos #theta_{CM} -0.25-(-0.5))");
    PhiSc545AsymmCM6->Fit("SinFit", "Q");
    InitialSinAmp[5][4] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[5][4] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc545AsymmCM6->Fit("AsymmFit", "QLL");
    SinAmp[5][4] = AsymmFit->GetParameter(0);
    SinAmpErr[5][4] = AsymmFit->GetParError(0);
    CosAmp[5][4] = AsymmFit->GetParameter(1);
    CosAmpErr[5][4] = AsymmFit->GetParError(1);

    PhiSc615AsymmCM6 = Phi_Scattered_615MeV_NegHelCM6->GetAsymmetry(Phi_Scattered_615MeV_PosHelCM6);
    PhiSc615AsymmCM6->SetName("PhiSc615AsymmCM6");
    PhiSc615AsymmCM6->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 580-650MeV (Cos #theta_{CM} -0.25-(-0.5))");
    PhiSc615AsymmCM6->Fit("SinFit", "Q");
    InitialSinAmp[5][5] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[5][5] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc615AsymmCM6->Fit("AsymmFit", "QLL");
    SinAmp[5][5] = AsymmFit->GetParameter(0);
    SinAmpErr[5][5] = AsymmFit->GetParError(0);
    CosAmp[5][5] = AsymmFit->GetParameter(1);
    CosAmpErr[5][5] = AsymmFit->GetParError(1);

    PhiSc685AsymmCM6 = Phi_Scattered_685MeV_NegHelCM6->GetAsymmetry(Phi_Scattered_685MeV_PosHelCM6);
    PhiSc685AsymmCM6->SetName("PhiSc685AsymmCM6");
    PhiSc685AsymmCM6->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 650-720MeV (Cos #theta_{CM} -0.25-(-0.5))");
    PhiSc685AsymmCM6->Fit("SinFit", "Q");
    InitialSinAmp[5][6] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[5][6] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc685AsymmCM6->Fit("AsymmFit", "QLL");
    SinAmp[5][6] = AsymmFit->GetParameter(0);
    SinAmpErr[5][6] = AsymmFit->GetParError(0);
    CosAmp[5][6] = AsymmFit->GetParameter(1);
    CosAmpErr[5][6] = AsymmFit->GetParError(1);

    ///////////////////////////////////////////
    //////////////////  CM7  //////////////////
    ///////////////////////////////////////////

    PhiSc265AsymmCM7 = Phi_Scattered_265MeV_NegHelCM7->GetAsymmetry(Phi_Scattered_265MeV_PosHelCM7);
    PhiSc265AsymmCM7->SetName("PhiSc265AsymmCM7");
    PhiSc265AsymmCM7->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 230-300MeV (Cos #theta_{CM} -0.5-(-0.75))");
    PhiSc265AsymmCM7->Fit("SinFit", "Q");
    InitialSinAmp[6][0] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[6][0] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc265AsymmCM7->Fit("AsymmFit", "QLL");
    SinAmp[6][0] = AsymmFit->GetParameter(0);
    SinAmpErr[6][0] = AsymmFit->GetParError(0);
    CosAmp[6][0] = AsymmFit->GetParameter(1);
    CosAmpErr[6][0] = AsymmFit->GetParError(1);

    PhiSc335AsymmCM7 = Phi_Scattered_335MeV_NegHelCM7->GetAsymmetry(Phi_Scattered_335MeV_PosHelCM7);
    PhiSc335AsymmCM7->SetName("PhiSc335AsymmCM7");
    PhiSc335AsymmCM7->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 300-370MeV (Cos #theta_{CM} -0.5-(-0.75))");
    PhiSc335AsymmCM7->Fit("SinFit", "Q");
    InitialSinAmp[6][1] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[6][1] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc335AsymmCM7->Fit("AsymmFit", "QLL");
    SinAmp[6][1] = AsymmFit->GetParameter(0);
    SinAmpErr[6][1] = AsymmFit->GetParError(0);
    CosAmp[6][1] = AsymmFit->GetParameter(1);
    CosAmpErr[6][1] = AsymmFit->GetParError(1);

    PhiSc405AsymmCM7 = Phi_Scattered_405MeV_NegHelCM7->GetAsymmetry(Phi_Scattered_405MeV_PosHelCM7);
    PhiSc405AsymmCM7->SetName("PhiSc405AsymmCM7");
    PhiSc405AsymmCM7->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 370-440MeV (Cos #theta_{CM} -0.5-(-0.75))");
    PhiSc405AsymmCM7->Fit("SinFit", "Q");
    InitialSinAmp[6][2] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[6][2] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc405AsymmCM7->Fit("AsymmFit", "QLL");
    SinAmp[6][2] = AsymmFit->GetParameter(0);
    SinAmpErr[6][2] = AsymmFit->GetParError(0);
    CosAmp[6][2] = AsymmFit->GetParameter(1);
    CosAmpErr[6][2] = AsymmFit->GetParError(1);

    PhiSc475AsymmCM7 = Phi_Scattered_475MeV_NegHelCM7->GetAsymmetry(Phi_Scattered_475MeV_PosHelCM7);
    PhiSc475AsymmCM7->SetName("PhiSc475AsymmCM7");
    PhiSc475AsymmCM7->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 440-510MeV (Cos #theta_{CM} -0.5-(-0.75))");
    PhiSc475AsymmCM7->Fit("SinFit", "Q");
    InitialSinAmp[6][3] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[6][3] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc475AsymmCM7->Fit("AsymmFit", "QLL");
    SinAmp[6][3] = AsymmFit->GetParameter(0);
    SinAmpErr[6][3] = AsymmFit->GetParError(0);
    CosAmp[6][3] = AsymmFit->GetParameter(1);
    CosAmpErr[6][3] = AsymmFit->GetParError(1);

    PhiSc545AsymmCM7 = Phi_Scattered_545MeV_NegHelCM7->GetAsymmetry(Phi_Scattered_545MeV_PosHelCM7);
    PhiSc545AsymmCM7->SetName("PhiSc545AsymmCM7");
    PhiSc545AsymmCM7->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 510-580MeV (Cos #theta_{CM} -0.5-(-0.75))");
    PhiSc545AsymmCM7->Fit("SinFit", "Q");
    InitialSinAmp[6][4] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[6][4] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc545AsymmCM7->Fit("AsymmFit", "QLL");
    SinAmp[6][4] = AsymmFit->GetParameter(0);
    SinAmpErr[6][4] = AsymmFit->GetParError(0);
    CosAmp[6][4] = AsymmFit->GetParameter(1);
    CosAmpErr[6][4] = AsymmFit->GetParError(1);

    PhiSc615AsymmCM7 = Phi_Scattered_615MeV_NegHelCM7->GetAsymmetry(Phi_Scattered_615MeV_PosHelCM7);
    PhiSc615AsymmCM7->SetName("PhiSc615AsymmCM7");
    PhiSc615AsymmCM7->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 580-650MeV (Cos #theta_{CM} -0.5-(-0.75))");
    PhiSc615AsymmCM7->Fit("SinFit", "Q");
    InitialSinAmp[6][5] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[6][5] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc615AsymmCM7->Fit("AsymmFit", "QLL");
    SinAmp[6][5] = AsymmFit->GetParameter(0);
    SinAmpErr[6][5] = AsymmFit->GetParError(0);
    CosAmp[6][5] = AsymmFit->GetParameter(1);
    CosAmpErr[6][5] = AsymmFit->GetParError(1);

    PhiSc685AsymmCM7 = Phi_Scattered_685MeV_NegHelCM7->GetAsymmetry(Phi_Scattered_685MeV_PosHelCM7);
    PhiSc685AsymmCM7->SetName("PhiSc685AsymmCM7");
    PhiSc685AsymmCM7->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 650-720MeV (Cos #theta_{CM} -0.5-(-0.75))");
    PhiSc685AsymmCM7->Fit("SinFit", "Q");
    InitialSinAmp[6][6] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[6][6] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc685AsymmCM7->Fit("AsymmFit", "QLL");
    SinAmp[6][6] = AsymmFit->GetParameter(0);
    SinAmpErr[6][6] = AsymmFit->GetParError(0);
    CosAmp[6][6] = AsymmFit->GetParameter(1);
    CosAmpErr[6][6] = AsymmFit->GetParError(1);

     ///////////////////////////////////////////
    //////////////////  CM8  //////////////////
    ///////////////////////////////////////////

    PhiSc265AsymmCM8 = Phi_Scattered_265MeV_NegHelCM8->GetAsymmetry(Phi_Scattered_265MeV_PosHelCM8);
    PhiSc265AsymmCM8->SetName("PhiSc265AsymmCM8");
    PhiSc265AsymmCM8->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 230-300MeV (Cos #theta_{CM} -0.75-(-1.0))");
    PhiSc265AsymmCM8->Fit("SinFit", "Q");
    InitialSinAmp[7][0] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[7][0] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc265AsymmCM8->Fit("AsymmFit", "QLL");
    SinAmp[7][0] = AsymmFit->GetParameter(0);
    SinAmpErr[7][0] = AsymmFit->GetParError(0);
    CosAmp[7][0] = AsymmFit->GetParameter(1);
    CosAmpErr[7][0] = AsymmFit->GetParError(1);

    PhiSc335AsymmCM8 = Phi_Scattered_335MeV_NegHelCM8->GetAsymmetry(Phi_Scattered_335MeV_PosHelCM8);
    PhiSc335AsymmCM8->SetName("PhiSc335AsymmCM8");
    PhiSc335AsymmCM8->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 300-370MeV (Cos #theta_{CM} -0.75-(-1.0))");
    PhiSc335AsymmCM8->Fit("SinFit", "Q");
    InitialSinAmp[7][1] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[7][1] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc335AsymmCM8->Fit("AsymmFit", "QLL");
    SinAmp[7][1] = AsymmFit->GetParameter(0);
    SinAmpErr[7][1] = AsymmFit->GetParError(0);
    CosAmp[7][1] = AsymmFit->GetParameter(1);
    CosAmpErr[7][1] = AsymmFit->GetParError(1);

    PhiSc405AsymmCM8 = Phi_Scattered_405MeV_NegHelCM8->GetAsymmetry(Phi_Scattered_405MeV_PosHelCM8);
    PhiSc405AsymmCM8->SetName("PhiSc405AsymmCM8");
    PhiSc405AsymmCM8->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 370-440MeV (Cos #theta_{CM} -0.75-(-1.0))");
    PhiSc405AsymmCM8->Fit("SinFit", "Q");
    InitialSinAmp[7][2] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[7][2] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc405AsymmCM8->Fit("AsymmFit", "QLL");
    SinAmp[7][2] = AsymmFit->GetParameter(0);
    SinAmpErr[7][2] = AsymmFit->GetParError(0);
    CosAmp[7][2] = AsymmFit->GetParameter(1);
    CosAmpErr[7][2] = AsymmFit->GetParError(1);

    PhiSc475AsymmCM8 = Phi_Scattered_475MeV_NegHelCM8->GetAsymmetry(Phi_Scattered_475MeV_PosHelCM8);
    PhiSc475AsymmCM8->SetName("PhiSc475AsymmCM8");
    PhiSc475AsymmCM8->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 440-510MeV (Cos #theta_{CM} -0.75-(-1.0))");
    PhiSc475AsymmCM8->Fit("SinFit", "Q");
    InitialSinAmp[7][3] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[7][3] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc475AsymmCM8->Fit("AsymmFit", "QLL");
    SinAmp[7][3] = AsymmFit->GetParameter(0);
    SinAmpErr[7][3] = AsymmFit->GetParError(0);
    CosAmp[7][3] = AsymmFit->GetParameter(1);
    CosAmpErr[7][3] = AsymmFit->GetParError(1);

    PhiSc545AsymmCM8 = Phi_Scattered_545MeV_NegHelCM8->GetAsymmetry(Phi_Scattered_545MeV_PosHelCM8);
    PhiSc545AsymmCM8->SetName("PhiSc545AsymmCM8");
    PhiSc545AsymmCM8->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 510-580MeV (Cos #theta_{CM} -0.75-(-1.0))");
    PhiSc545AsymmCM8->Fit("SinFit", "Q");
    InitialSinAmp[7][4] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[7][4] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc545AsymmCM8->Fit("AsymmFit", "QLL");
    SinAmp[7][4] = AsymmFit->GetParameter(0);
    SinAmpErr[7][4] = AsymmFit->GetParError(0);
    CosAmp[7][4] = AsymmFit->GetParameter(1);
    CosAmpErr[7][4] = AsymmFit->GetParError(1);

    PhiSc615AsymmCM8 = Phi_Scattered_615MeV_NegHelCM8->GetAsymmetry(Phi_Scattered_615MeV_PosHelCM8);
    PhiSc615AsymmCM8->SetName("PhiSc615AsymmCM8");
    PhiSc615AsymmCM8->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 580-650MeV (Cos #theta_{CM} -0.75-(-1.0))");
    PhiSc615AsymmCM8->Fit("SinFit", "Q");
    InitialSinAmp[7][5] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[7][5] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc615AsymmCM8->Fit("AsymmFit", "QLL");
    SinAmp[7][5] = AsymmFit->GetParameter(0);
    SinAmpErr[7][5] = AsymmFit->GetParError(0);
    CosAmp[7][5] = AsymmFit->GetParameter(1);
    CosAmpErr[7][5] = AsymmFit->GetParError(1);

    PhiSc685AsymmCM8 = Phi_Scattered_685MeV_NegHelCM8->GetAsymmetry(Phi_Scattered_685MeV_PosHelCM8);
    PhiSc685AsymmCM8->SetName("PhiSc685AsymmCM8");
    PhiSc685AsymmCM8->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 650-720MeV (Cos #theta_{CM} -0.75-(-1.0))");
    PhiSc685AsymmCM8->Fit("SinFit", "Q");
    InitialSinAmp[7][6] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[7][6] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -1, 1);
    AsymmFunc->SetParLimits(1, -1, 1);
    PhiSc685AsymmCM8->Fit("AsymmFit", "QLL");
    SinAmp[7][6] = AsymmFit->GetParameter(0);
    SinAmpErr[7][6] = AsymmFit->GetParError(0);
    CosAmp[7][6] = AsymmFit->GetParameter(1);
    CosAmpErr[7][6] = AsymmFit->GetParError(1);

    // Define new file to store fit parameters
    TFile f1("AsymmFits_PTotal_AmoOnly_105.root", "RECREATE");

    PhiSc265AsymmCM1->Write();
    PhiSc335AsymmCM1->Write();
    PhiSc405AsymmCM1->Write();
    PhiSc475AsymmCM1->Write();
    PhiSc545AsymmCM1->Write();
    PhiSc615AsymmCM1->Write();
    PhiSc685AsymmCM1->Write();

    PhiSc265AsymmCM2->Write();
    PhiSc335AsymmCM2->Write();
    PhiSc405AsymmCM2->Write();
    PhiSc475AsymmCM2->Write();
    PhiSc545AsymmCM2->Write();
    PhiSc615AsymmCM2->Write();
    PhiSc685AsymmCM2->Write();

    PhiSc265AsymmCM3->Write();
    PhiSc335AsymmCM3->Write();
    PhiSc405AsymmCM3->Write();
    PhiSc475AsymmCM3->Write();
    PhiSc545AsymmCM3->Write();
    PhiSc615AsymmCM3->Write();
    PhiSc685AsymmCM3->Write();

    PhiSc265AsymmCM4->Write();
    PhiSc335AsymmCM4->Write();
    PhiSc405AsymmCM4->Write();
    PhiSc475AsymmCM4->Write();
    PhiSc545AsymmCM4->Write();
    PhiSc615AsymmCM4->Write();
    PhiSc685AsymmCM4->Write();

    PhiSc265AsymmCM5->Write();
    PhiSc335AsymmCM5->Write();
    PhiSc405AsymmCM5->Write();
    PhiSc475AsymmCM5->Write();
    PhiSc545AsymmCM5->Write();
    PhiSc615AsymmCM5->Write();
    PhiSc685AsymmCM5->Write();

    PhiSc265AsymmCM6->Write();
    PhiSc335AsymmCM6->Write();
    PhiSc405AsymmCM6->Write();
    PhiSc475AsymmCM6->Write();
    PhiSc545AsymmCM6->Write();
    PhiSc615AsymmCM6->Write();
    PhiSc685AsymmCM6->Write();

    PhiSc265AsymmCM7->Write();
    PhiSc335AsymmCM7->Write();
    PhiSc405AsymmCM7->Write();
    PhiSc475AsymmCM7->Write();
    PhiSc545AsymmCM7->Write();
    PhiSc615AsymmCM7->Write();
    PhiSc685AsymmCM7->Write();

    PhiSc265AsymmCM8->Write();
    PhiSc335AsymmCM8->Write();
    PhiSc405AsymmCM8->Write();
    PhiSc475AsymmCM8->Write();
    PhiSc545AsymmCM8->Write();
    PhiSc615AsymmCM8->Write();
    PhiSc685AsymmCM8->Write();

    //Define new tree to store parameters in
    TTree* tree = new TTree("Parameter_Values", "Tree_of_Values");

    // Define branches to store parameters, (Branch Name, Variable, Type of Variable)

    tree->Branch("InitialSinAmp265", &ISinAm265, "ISinAm265/D");
    tree->Branch("InitialSinAmpErr265", &ISinAmErr265, "ISinAmErr265/D");
    tree->Branch("CosAmp265", &CosAm265, "CosAm265/D");
    tree->Branch("CosAmpErr265", &CosAmErr265, "CosAmErr265/D");
    tree->Branch("SinAmp265", &SinAm265, "SinAm265/D");
    tree->Branch("SinAmpErr265", &SinAmErr265, "SinAmErr265/D");
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

        ISinAm265 = InitialSinAmp[m][0];
        ISinAmErr265 = InitialSinAmpErr[m][0];
        SinAm265 = SinAmp[m][0];
        SinAmErr265 = SinAmpErr[m][0];
        CosAm265 = CosAmp[m][0];
        CosAmErr265 = CosAmpErr[m][1];
        ISinAm335 = InitialSinAmp[m][1];
        ISinAmErr335 = InitialSinAmpErr[m][1];
        SinAm335 = SinAmp[m][1];
        SinAmErr335 = SinAmpErr[m][1];
        CosAm335 = CosAmp[m][1];
        CosAmErr335 = CosAmpErr[m][1];
        ISinAm405 = InitialSinAmp[m][2];
        ISinAmErr405 = InitialSinAmpErr[m][2];
        SinAm405 = SinAmp[m][2];
        SinAmErr405 = SinAmpErr[m][2];
        CosAm405 = CosAmp[m][2];
        CosAmErr405 = CosAmpErr[m][2];
        ISinAm475 = InitialSinAmp[m][3];
        ISinAmErr475 = InitialSinAmpErr[m][3];
        SinAm475 = SinAmp[m][3];
        SinAmErr475 = SinAmpErr[m][3];
        CosAm475 = CosAmp[m][3];
        CosAmErr475 = CosAmpErr[m][3];
        ISinAm545 = InitialSinAmp[m][4];
        ISinAmErr545 = InitialSinAmpErr[m][4];
        SinAm545 = SinAmp[m][4];
        SinAmErr545 = SinAmpErr[m][4];
        CosAm545 = CosAmp[m][4];
        CosAmErr545 = CosAmpErr[m][4];
        ISinAm615 = InitialSinAmp[m][5];
        ISinAmErr615 = InitialSinAmpErr[m][5];
        SinAm615 = SinAmp[m][5];
        SinAmErr615 = SinAmpErr[m][5];
        CosAm615 = CosAmp[m][5];
        CosAmErr615 = CosAmpErr[m][5];
        ISinAm685 = InitialSinAmp[m][6];
        ISinAmErr685 = InitialSinAmpErr[m][6];
        SinAm685 = SinAmp[m][6];
        SinAmErr685 = SinAmpErr[m][6];
        CosAm685 = CosAmp[m][6];
        CosAmErr685 = CosAmpErr[m][6];

        tree->Fill();

    }

    f1.Write();
    f1.Close();

    TFile *f2= TFile::Open("/scratch/Mainz_Software/a2GoAT/CircPol_Aug16.root");

    // Calculate values of Cx for each angular and energy bin
    for (Int_t n = 0; n < 8; n++){

        Cx265[n] = SinAmp[n][0]/(0.1*(Graph->Eval(265,0)));
        CxErr265[n] = SinAmpErr[n][0]/((0.1)*(Graph->Eval(265,0)));
        Cx335[n] = SinAmp[n][1]/(0.1*(Graph->Eval(335,0)));
        CxErr335[n] = SinAmpErr[n][1]/((0.1)*(Graph->Eval(335,0)));
        Cx405[n] = SinAmp[n][2]/(0.1*(Graph->Eval(405,0)));
        CxErr405[n] = SinAmpErr[n][2]/((0.1)*(Graph->Eval(405,0)));
        Cx475[n] = SinAmp[n][3]/(0.1*(Graph->Eval(475,0)));
        CxErr475[n] = SinAmpErr[n][3]/((0.1)*(Graph->Eval(475,0)));
        Cx545[n] = SinAmp[n][4]/(0.1*(Graph->Eval(545,0)));
        CxErr545[n] = SinAmpErr[n][4]/((0.1)*(Graph->Eval(545,0)));
        Cx615[n] = SinAmp[n][5]/(0.1*(Graph->Eval(615,0)));
        CxErr615[n] = SinAmpErr[n][5]/((0.1)*(Graph->Eval(615,0)));
        Cx685[n] = SinAmp[n][6]/(0.1*(Graph->Eval(685,0)));
        CxErr685[n] = SinAmpErr[n][6]/((0.1)*(Graph->Eval(685,0)));

        // NOTE - This is to show amplitude only! Need to divide by Ay too!
        Pn265[n] = CosAmp[n][0];
        PnErr265[n] = CosAmpErr[n][0];
        Pn335[n] = CosAmp[n][1];
        PnErr335[n] = CosAmpErr[n][1];
        Pn405[n] = CosAmp[n][2];
        PnErr405[n] = CosAmpErr[n][2];
        Pn475[n] = CosAmp[n][3];
        PnErr475[n] = CosAmpErr[n][3];
        Pn545[n] = CosAmp[n][4];
        PnErr545[n] = CosAmpErr[n][4];
        Pn615[n] = CosAmp[n][5];
        PnErr615[n] = CosAmpErr[n][5];
        Pn685[n] = CosAmp[n][6];
        PnErr685[n] = CosAmpErr[n][6];

    }

    TFile f3("Cx_Plots_AmoOnly_105.root", "RECREATE");

    Float_t xMin = -1;
    Float_t xMax = 1;
    Float_t yMin = -5;
    Float_t yMax = 5;
    Double_t x[8] = {0.875, 0.625, 0.375, 0.125, -0.125, -0.375, -0.625, -0.875}; // Need to adjust
    Double_t ex[8] = {0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125}; // Need to adjust

    TCanvas *canvas = new TCanvas("canvas","canvas", 1920, 1080);
    TPad *pad = new TPad("pad","",0,0,1,1);
    pad->Draw();
    pad->cd();

    pad->SetTickx(1);
    pad->SetTicky(1);
    pad->SetGridx(1);
    pad->SetGridy(1);
    TH1F  *hr;
    hr = canvas->DrawFrame(xMin, -1,xMax, 1);
    hr->SetTitle("C_{x} as fn of Cos#theta_{CM} (E_{#gamma} 230-300MeV)");

    gr = new TGraphErrors(8, x, Cx265, ex, CxErr265);
    gr->SetMarkerColor(2);
    gr->SetMarkerStyle(5);
    gr->SetMarkerSize(2);
    gr->SetTitle("C_{x} as fn of Cos#theta_{CM} (E_{#gamma} 230-300MeV)");
    gr->GetXaxis()->SetTitle("Cos#theta_{CM}");
    gr->GetYaxis()->SetTitle("C_{x}");
    gr->SetName("Cx265");
    gr->Draw("ep");

    TCanvas *canvas1 = new TCanvas("canvas1","canvas1", 1920, 1080);
    TPad *pad1 = new TPad("pad1","",0,0,1,1);
    pad1->Draw();
    pad1->cd();

    pad1->SetTickx(1);
    pad1->SetTicky(1);
    pad1->SetGridx(1);
    pad1->SetGridy(1);
    TH1F  *hr1;
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

    canvas->Write();
    gr->Write();
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
    canvas7->Divide(4,2);
    canvas7->cd(1);
    pad->Draw();
    canvas7->cd(2);
    pad1->Draw();
    canvas7->cd(3);
    pad2->Draw();
    canvas7->cd(4);
    pad3->Draw();
    canvas7->cd(5);
    pad4->Draw();
    canvas7->cd(6);
    pad5->Draw();
    canvas7->cd(7);
    pad6->Draw();
    canvas7->Write();

    TCanvas *canvas8 = new TCanvas("canvas8","canvas8", 1920, 1080);
    TPad *pad7 = new TPad("pad7","",0,0,1,1);
    pad7->Draw();
    pad7->cd();

    pad7->SetTickx(1);
    pad7->SetTicky(1);
    pad7->SetGridx(1);
    pad7->SetGridy(1);
    TH1F  *hr7;
    hr7 = canvas->DrawFrame(xMin, -10, xMax, 10);
    hr7->SetTitle("P_{n} as fn of Cos#theta_{CM} (E_{#gamma} 230-300MeV)");

    gr7 = new TGraphErrors(8, x, Pn265, ex, PnErr265);
    gr7->SetMarkerColor(2);
    gr7->SetMarkerStyle(5);
    gr7->SetMarkerSize(2);
    gr7->SetTitle("P_{n} as fn of Cos#theta_{CM} (E_{#gamma} 230-300MeV)");
    gr7->GetXaxis()->SetTitle("Cos#theta_{CM}");
    gr7->GetYaxis()->SetTitle("P_{n}");
    gr7->SetName("Pn265");
    gr7->Draw("ep");

    TCanvas *canvas9 = new TCanvas("canvas9","canvas9", 1920, 1080);
    TPad *pad8 = new TPad("pad8","",0,0,1,1);
    pad8->Draw();
    pad8->cd();

    pad8->SetTickx(1);
    pad8->SetTicky(1);
    pad8->SetGridx(1);
    pad8->SetGridy(1);
    TH1F  *hr8;
    hr8 = canvas9->DrawFrame(xMin, -10, xMax, 10);
    hr8->SetTitle("P_{n} as fn of Cos#theta_{CM} (E_{#gamma} 300-370MeV)");

    gr8 = new TGraphErrors(8, x, Pn335, ex, PnErr335);
    gr8->SetMarkerColor(2);
    gr8->SetMarkerStyle(5);
    gr8->SetMarkerSize(2);
    gr8->SetTitle("P_{n} as fn of Cos#theta_{CM} (E_{#gamma} 300-370MeV)");
    gr8->GetXaxis()->SetTitle("Cos#theta_{CM}");
    gr8->GetYaxis()->SetTitle("P_{n}");
    gr8->SetName("Pn335");
    gr8->Draw("ep");

    TCanvas *canvas10 = new TCanvas("canvas10","canvas10", 1920, 1080);
    TPad *pad9 = new TPad("pad9","",0,0,1,1);
    pad9->Draw();
    pad9->cd();

    pad9->SetTickx(1);
    pad9->SetTicky(1);
    pad9->SetGridx(1);
    pad9->SetGridy(1);
    TH1F  *hr9;
    hr9 = canvas10->DrawFrame(xMin, -10, xMax, 10);
    hr9->SetTitle("P_{n} as fn of Cos#theta_{CM} (E_{#gamma} 370-440MeV)");

    gr9 = new TGraphErrors(8, x, Pn405, ex, PnErr405);
    gr9->SetMarkerColor(2);
    gr9->SetMarkerStyle(5);
    gr9->SetMarkerSize(2);
    gr9->SetTitle("P_{n} as fn of Cos#theta_{CM} (E_{#gamma} 370-440MeV)");
    gr9->GetXaxis()->SetTitle("Cos#theta_{CM}");
    gr9->GetYaxis()->SetTitle("P_{n}");
    gr9->SetName("Pn405");
    gr9->Draw("ep");

    TCanvas *canvas11 = new TCanvas("canvas11","canvas11", 1920, 1080);
    TPad *pad10 = new TPad("pad10","",0,0,1,1);
    pad10->Draw();
    pad10->cd();

    pad10->SetTickx(1);
    pad10->SetTicky(1);
    pad10->SetGridx(1);
    pad10->SetGridy(1);
    TH1F  *hr10;
    hr10 = canvas11->DrawFrame(xMin, -10, xMax, 10);
    hr10->SetTitle("P_{n} as fn of Cos#theta_{CM} (E_{#gamma} 440-510MeV)");

    gr10 = new TGraphErrors(8, x, Pn475, ex, PnErr475);
    gr10->SetMarkerColor(2);
    gr10->SetMarkerStyle(5);
    gr10->SetMarkerSize(2);
    gr10->SetTitle("P_{n} as fn of Cos#theta_{CM} (E_{#gamma} 440-510MeV)");
    gr10->GetXaxis()->SetTitle("Cos#theta_{CM}");
    gr10->GetYaxis()->SetTitle("P_{n}");
    gr10->SetName("Pn475");
    gr10->Draw("ep");

    TCanvas *canvas12 = new TCanvas("canvas12","canvas12", 1920, 1080);
    TPad *pad11 = new TPad("pad11","",0,0,1,1);
    pad11->Draw();
    pad11->cd();

    pad11->SetTickx(1);
    pad11->SetTicky(1);
    pad11->SetGridx(1);
    pad11->SetGridy(1);
    TH1F  *hr11;
    hr11 = canvas12->DrawFrame(xMin, -10, xMax, 10);
    hr11->SetTitle("P_{n} as fn of Cos#theta_{CM} (E_{#gamma} 510-580MeV)");

    gr11 = new TGraphErrors(8, x, Pn545, ex, PnErr545);
    gr11->SetMarkerColor(2);
    gr11->SetMarkerStyle(5);
    gr11->SetMarkerSize(2);
    gr11->SetTitle("P_{n} as fn of Cos#theta_{CM} (E_{#gamma} 510-580MeV)");
    gr11->GetXaxis()->SetTitle("Cos#theta_{CM}");
    gr11->GetYaxis()->SetTitle("P_{n}");
    gr11->SetName("Pn545");
    gr11->Draw("ep");

    TCanvas *canvas13 = new TCanvas("canvas13","canvas13", 1920, 1080);
    TPad *pad12 = new TPad("pad12","",0,0,1,1);
    pad12->Draw();
    pad12->cd();

    pad12->SetTickx(1);
    pad12->SetTicky(1);
    pad12->SetGridx(1);
    pad12->SetGridy(1);
    TH1F  *hr12;
    hr12 = canvas13->DrawFrame(xMin, -10, xMax, 10);
    hr12->SetTitle("P_{n} as fn of Cos#theta_{CM} (E_{#gamma} 580-650MeV)");

    gr12 = new TGraphErrors(8, x, Pn615, ex, PnErr615);
    gr12->SetMarkerColor(2);
    gr12->SetMarkerStyle(5);
    gr12->SetMarkerSize(2);
    gr12->SetTitle("P_{n} as fn of Cos#theta_{CM} (E_{#gamma} 580-650MeV)");
    gr12->GetXaxis()->SetTitle("Cos#theta_{CM}");
    gr12->GetYaxis()->SetTitle("P_{n}");
    gr12->SetName("Pn615");
    gr12->Draw("ep");

    TCanvas *canvas14 = new TCanvas("canvas14","canvas14", 1920, 1080);
    TPad *pad13 = new TPad("pad13","",0,0,1,1);
    pad13->Draw();
    pad13->cd();

    pad13->SetTickx(1);
    pad13->SetTicky(1);
    pad13->SetGridx(1);
    pad13->SetGridy(1);
    TH1F  *hr13;
    hr13 = canvas14->DrawFrame(xMin, -10, xMax, 10);
    hr13->SetTitle("P_{n} as fn of Cos#theta_{CM} (E_{#gamma} 650-720MeV)");

    gr13 = new TGraphErrors(8, x, Pn685, ex, PnErr685);
    gr13->SetMarkerColor(2);
    gr13->SetMarkerStyle(5);
    gr13->SetMarkerSize(2);
    gr13->SetTitle("P_{n} as fn of Cos#theta_{CM} (E_{#gamma} 650-720MeV)");
    gr13->GetXaxis()->SetTitle("Cos#theta_{CM}");
    gr13->GetYaxis()->SetTitle("P_{n}");
    gr13->SetName("Pn685");
    gr13->Draw("ep");

    TCanvas *canvas15 = new TCanvas("canvas15","canvas15", 1920, 1080);
    canvas15->Divide(4,2);
    canvas15->cd(1);
    pad7->Draw();
    canvas15->cd(2);
    pad8->Draw();
    canvas15->cd(3);
    pad9->Draw();
    canvas15->cd(4);
    pad10->Draw();
    canvas15->cd(5);
    pad11->Draw();
    canvas15->cd(6);
    pad12->Draw();
    canvas15->cd(7);
    pad13->Draw();
    canvas15->Write();

    f3.Write();

}
