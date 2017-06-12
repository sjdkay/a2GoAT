#include "./includes.h"

// define a function with 4 parameters
Double_t fitf(Double_t *x,Double_t *par)
{
    Double_t fitval = 0;
    fitval =  par[0] + ((par[1]*sin(x[0]*TMath::DegToRad()))/(1 + (par[2]*cos(x[0]*TMath::DegToRad()))));
    return fitval;
}

void AsymmNew(){

    double InitialSinAmp[8][13];
    double InitialSinAmpErr[8][13];
    double Offset[8][13]; // Format of arrays is Theta Bin Defined by i, Energy bin defined by j
    double OffsetErr[8][13];
    double SinAmp[8][13]; // Format of arrays is Theta Bin Defined by i, Energy bin defined by j
    double SinAmpErr[8][13];
    double CosAmp[8][13]; // Format of arrays is Theta Bin Defined by i, Energy bin defined by j
    double CosAmpErr[8][13];
    Int_t i;
    double ISinAm315;
    double ISinAmErr315;
    double Offs315;
    double OffsErr315;
    double SinAm315;
    double SinAmErr315;
    double CosAm315;
    double CosAmErr315;
    double ISinAm345;
    double ISinAmErr345;
    double Offs345;
    double OffsErr345;
    double SinAm345;
    double SinAmErr345;
    double CosAm345;
    double CosAmErr345;
    double ISinAm375;
    double ISinAmErr375;
    double Offs375;
    double OffsErr375;
    double SinAm375;
    double SinAmErr375;
    double CosAm375;
    double CosAmErr375;
    double ISinAm405;
    double ISinAmErr405;
    double Offs405;
    double OffsErr405;
    double SinAm405;
    double SinAmErr405;
    double CosAm405;
    double CosAmErr405;
    double ISinAm435;
    double ISinAmErr435;
    double Offs435;
    double OffsErr435;
    double SinAm435;
    double SinAmErr435;
    double CosAm435;
    double CosAmErr435;
    double ISinAm465;
    double ISinAmErr465;
    double Offs465;
    double OffsErr465;
    double SinAm465;
    double SinAmErr465;
    double CosAm465;
    double CosAmErr465;
    double ISinAm495;
    double ISinAmErr495;
    double Offs495;
    double OffsErr495;
    double SinAm495;
    double SinAmErr495;
    double CosAm495;
    double CosAmErr495;
    double ISinAm525;
    double ISinAmErr525;
    double Offs525;
    double OffsErr525;
    double SinAm525;
    double SinAmErr525;
    double CosAm525;
    double CosAmErr525;
    double ISinAm555;
    double ISinAmErr555;
    double Offs555;
    double OffsErr555;
    double SinAm555;
    double SinAmErr555;
    double CosAm555;
    double CosAmErr555;
    double ISinAm585;
    double ISinAmErr585;
    double Offs585;
    double OffsErr585;
    double SinAm585;
    double SinAmErr585;
    double CosAm585;
    double CosAmErr585;
    double ISinAm615;
    double ISinAmErr615;
    double Offs615;
    double OffsErr615;
    double SinAm615;
    double SinAmErr615;
    double CosAm615;
    double CosAmErr615;
    double ISinAm645;
    double ISinAmErr645;
    double Offs645;
    double OffsErr645;
    double SinAm645;
    double SinAmErr645;
    double CosAm645;
    double CosAmErr645;
    double ISinAm675;
    double ISinAmErr675;
    double Offs675;
    double OffsErr675;
    double SinAm675;
    double SinAmErr675;
    double CosAm675;
    double CosAmErr675;
    Double_t Cx315[8], Cx345[8], Cx375[8], Cx405[8], Cx435[8], Cx465[8], Cx495[8], Cx525[8], Cx555[8], Cx585[8], Cx615[8], Cx645[8], Cx675[8];
    Double_t CxErr315[8], CxErr345[8], CxErr375[8], CxErr405[8], CxErr435[8], CxErr465[8], CxErr495[8], CxErr525[8], CxErr555[8], CxErr585[8], CxErr615[8], CxErr645[8], CxErr675[8];

    TF1 *AsymmFunc = new TF1("AsymmFit",  fitf, -130.0, 130.0, 3); //Give a name and range to the fitting funcion
    AsymmFunc->SetParNames("Offset", "SinAmp", "CosAmp"); //Name the parameters
    AsymmFunc->SetParameter(0, 0);
    TF1 *SinFunc = new TF1("SinFit", "[0]*sin(x*TMath::DegToRad())", -130, 130);
    SinFunc->SetParNames("InitialSinAmp");
    TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/Amo/Physics_Total_48_12_6_17.root"); // Open the latest PTotal file to load histograms from

    ///////////////////////////////////////////
    //////////////////  CM1  //////////////////
    ///////////////////////////////////////////

    PhiSc315AsymmCM1 = PhiSc315NegHelCM1->GetAsymmetry(PhiSc315PosHelCM1);
    PhiSc315AsymmCM1->SetName("PhiSc315AsymmCM1");
    PhiSc315AsymmCM1->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 300-330MeV (Cos #theta_{CM} 1-0.75)")
    PhiSc315AsymmCM1->Fit("SinFit", "Q");
    InitialSinAmp[0][0] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[0][0] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc315AsymmCM1->Fit("AsymmFit", "Q");
    Offset[0][0] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[0][0] = AsymmFit->GetParError(0);
    SinAmp[0][0] = AsymmFit->GetParameter(1);
    SinAmpErr[0][0] = AsymmFit->GetParError(1);
    CosAmp[0][0] = AsymmFit->GetParameter(2);
    CosAmpErr[0][0] = AsymmFit->GetParError(2);

    PhiSc345AsymmCM1 = PhiSc345NegHelCM1->GetAsymmetry(PhiSc345PosHelCM1);
    PhiSc345AsymmCM1->SetName("PhiSc345AsymmCM1");
    PhiSc345AsymmCM1->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 330-360MeV (Cos #theta_{CM} 1-0.75)")
    PhiSc345AsymmCM1->Fit("SinFit", "Q");
    InitialSinAmp[0][1] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[0][1] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc345AsymmCM1->Fit("AsymmFit", "Q");
    Offset[0][1] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[0][1] = AsymmFit->GetParError(0);
    SinAmp[0][1] = AsymmFit->GetParameter(1);
    SinAmpErr[0][1] = AsymmFit->GetParError(1);
    CosAmp[0][1] = AsymmFit->GetParameter(2);
    CosAmpErr[0][1] = AsymmFit->GetParError(2);

    PhiSc375AsymmCM1 = PhiSc375NegHelCM1->GetAsymmetry(PhiSc375PosHelCM1);
    PhiSc375AsymmCM1->SetName("PhiSc375AsymmCM1");
    PhiSc375AsymmCM1->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 360-390MeV (Cos #theta_{CM} 1-0.75)")
    PhiSc375AsymmCM1->Fit("SinFit", "Q");
    InitialSinAmp[0][2] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[0][2] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc375AsymmCM1->Fit("AsymmFit", "Q");
    Offset[0][2] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[0][2] = AsymmFit->GetParError(0);
    SinAmp[0][2] = AsymmFit->GetParameter(1);
    SinAmpErr[0][2] = AsymmFit->GetParError(1);
    CosAmp[0][2] = AsymmFit->GetParameter(2);
    CosAmpErr[0][2] = AsymmFit->GetParError(2);

    PhiSc405AsymmCM1 = PhiSc405NegHelCM1->GetAsymmetry(PhiSc405PosHelCM1);
    PhiSc405AsymmCM1->SetName("PhiSc405AsymmCM1");
    PhiSc405AsymmCM1->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 390-420MeV (Cos #theta_{CM} 1-0.75)")
    PhiSc405AsymmCM1->Fit("SinFit", "Q");
    InitialSinAmp[0][3] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[0][3] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc405AsymmCM1->Fit("AsymmFit", "Q");
    Offset[0][3] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[0][3] = AsymmFit->GetParError(0);
    SinAmp[0][3] = AsymmFit->GetParameter(1);
    SinAmpErr[0][3] = AsymmFit->GetParError(1);
    CosAmp[0][3] = AsymmFit->GetParameter(2);
    CosAmpErr[0][3] = AsymmFit->GetParError(2);

    PhiSc435AsymmCM1 = PhiSc435NegHelCM1->GetAsymmetry(PhiSc435PosHelCM1);
    PhiSc435AsymmCM1->SetName("PhiSc435AsymmCM1");
    PhiSc435AsymmCM1->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 420-450MeV (Cos #theta_{CM} 1-0.75)")
    PhiSc435AsymmCM1->Fit("SinFit", "Q");
    InitialSinAmp[0][4] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[0][4] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc435AsymmCM1->Fit("AsymmFit", "Q");
    Offset[0][4] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[0][4] = AsymmFit->GetParError(0);
    SinAmp[0][4] = AsymmFit->GetParameter(1);
    SinAmpErr[0][4] = AsymmFit->GetParError(1);
    CosAmp[0][4] = AsymmFit->GetParameter(2);
    CosAmpErr[0][4] = AsymmFit->GetParError(2);

    PhiSc465AsymmCM1 = PhiSc465NegHelCM1->GetAsymmetry(PhiSc465PosHelCM1);
    PhiSc465AsymmCM1->SetName("PhiSc465AsymmCM1");
    PhiSc465AsymmCM1->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 450-480MeV (Cos #theta_{CM} 1-0.75)")
    PhiSc465AsymmCM1->Fit("SinFit", "Q");
    InitialSinAmp[0][5] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[0][5] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc465AsymmCM1->Fit("AsymmFit", "Q");
    Offset[0][5] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[0][5] = AsymmFit->GetParError(0);
    SinAmp[0][5] = AsymmFit->GetParameter(1);
    SinAmpErr[0][5] = AsymmFit->GetParError(1);
    CosAmp[0][5] = AsymmFit->GetParameter(2);
    CosAmpErr[0][5] = AsymmFit->GetParError(2);

    PhiSc495AsymmCM1 = PhiSc495NegHelCM1->GetAsymmetry(PhiSc495PosHelCM1);
    PhiSc495AsymmCM1->SetName("PhiSc495AsymmCM1");
    PhiSc495AsymmCM1->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 480-510MeV (Cos #theta_{CM} 1-0.75)")
    PhiSc495AsymmCM1->Fit("SinFit", "Q");
    InitialSinAmp[0][6] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[0][6] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc495AsymmCM1->Fit("AsymmFit", "Q");
    Offset[0][6] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[0][6] = AsymmFit->GetParError(0);
    SinAmp[0][6] = AsymmFit->GetParameter(1);
    SinAmpErr[0][6] = AsymmFit->GetParError(1);
    CosAmp[0][6] = AsymmFit->GetParameter(2);
    CosAmpErr[0][6] = AsymmFit->GetParError(2);

    PhiSc525AsymmCM1 = PhiSc525NegHelCM1->GetAsymmetry(PhiSc525PosHelCM1);
    PhiSc525AsymmCM1->SetName("PhiSc525AsymmCM1");
    PhiSc525AsymmCM1->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 510-5400MeV (Cos #theta_{CM} 1-0.75)")
    PhiSc525AsymmCM1->Fit("SinFit", "Q");
    InitialSinAmp[0][7] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[0][7] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc525AsymmCM1->Fit("AsymmFit", "Q");
    Offset[0][7] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[0][7] = AsymmFit->GetParError(0);
    SinAmp[0][7] = AsymmFit->GetParameter(1);
    SinAmpErr[0][7] = AsymmFit->GetParError(1);
    CosAmp[0][7] = AsymmFit->GetParameter(2);
    CosAmpErr[0][7] = AsymmFit->GetParError(2);

    PhiSc555AsymmCM1 = PhiSc555NegHelCM1->GetAsymmetry(PhiSc555PosHelCM1);
    PhiSc555AsymmCM1->SetName("PhiSc555AsymmCM1");
    PhiSc555AsymmCM1->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 540-570MeV (Cos #theta_{CM} 1-0.75)")
    PhiSc555AsymmCM1->Fit("SinFit", "Q");
    InitialSinAmp[0][8] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[0][8] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc555AsymmCM1->Fit("AsymmFit", "Q");
    Offset[0][8] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[0][8] = AsymmFit->GetParError(0);
    SinAmp[0][8] = AsymmFit->GetParameter(1);
    SinAmpErr[0][8] = AsymmFit->GetParError(1);
    CosAmp[0][8] = AsymmFit->GetParameter(2);
    CosAmpErr[0][8] = AsymmFit->GetParError(2);

    PhiSc585AsymmCM1 = PhiSc585NegHelCM1->GetAsymmetry(PhiSc585PosHelCM1);
    PhiSc585AsymmCM1->SetName("PhiSc585AsymmCM1");
    PhiSc585AsymmCM1->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 570-600MeV (Cos #theta_{CM} 1-0.75)")
    PhiSc585AsymmCM1->Fit("SinFit", "Q");
    InitialSinAmp[0][9] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[0][9] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc585AsymmCM1->Fit("AsymmFit", "Q");
    Offset[0][9] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[0][9] = AsymmFit->GetParError(0);
    SinAmp[0][9] = AsymmFit->GetParameter(1);
    SinAmpErr[0][9] = AsymmFit->GetParError(1);
    CosAmp[0][9] = AsymmFit->GetParameter(2);
    CosAmpErr[0][9] = AsymmFit->GetParError(2);

    PhiSc615AsymmCM1 = PhiSc615NegHelCM1->GetAsymmetry(PhiSc615PosHelCM1);
    PhiSc615AsymmCM1->SetName("PhiSc615AsymmCM1");
    PhiSc615AsymmCM1->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 600-630MeV (Cos #theta_{CM} 1-0.75)")
    PhiSc615AsymmCM1->Fit("SinFit", "Q");
    InitialSinAmp[0][10] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[0][10] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc615AsymmCM1->Fit("AsymmFit", "Q");
    Offset[0][10] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[0][10] = AsymmFit->GetParError(0);
    SinAmp[0][10] = AsymmFit->GetParameter(1);
    SinAmpErr[0][10] = AsymmFit->GetParError(1);
    CosAmp[0][10] = AsymmFit->GetParameter(2);
    CosAmpErr[0][10] = AsymmFit->GetParError(2);

    PhiSc645AsymmCM1 = PhiSc645NegHelCM1->GetAsymmetry(PhiSc645PosHelCM1);
    PhiSc645AsymmCM1->SetName("PhiSc645AsymmCM1");
    PhiSc645AsymmCM1->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 630-660MeV (Cos #theta_{CM} 1-0.75)")
    PhiSc645AsymmCM1->Fit("SinFit", "Q");
    InitialSinAmp[0][11] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[0][11] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc645AsymmCM1->Fit("AsymmFit", "Q");
    Offset[0][11] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[0][11] = AsymmFit->GetParError(0);
    SinAmp[0][11] = AsymmFit->GetParameter(1);
    SinAmpErr[0][11] = AsymmFit->GetParError(1);
    CosAmp[0][11] = AsymmFit->GetParameter(2);
    CosAmpErr[0][11] = AsymmFit->GetParError(2);

    PhiSc675AsymmCM1 = PhiSc675NegHelCM1->GetAsymmetry(PhiSc675PosHelCM1);
    PhiSc675AsymmCM1->SetName("PhiSc675AsymmCM1");
    PhiSc675AsymmCM1->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 660-690MeV (Cos #theta_{CM} 1-0.75)")
    PhiSc675AsymmCM1->Fit("SinFit", "Q");
    InitialSinAmp[0][12] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[0][12] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc675AsymmCM1->Fit("AsymmFit", "Q");
    Offset[0][12] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[0][12] = AsymmFit->GetParError(0);
    SinAmp[0][12] = AsymmFit->GetParameter(1);
    SinAmpErr[0][12] = AsymmFit->GetParError(1);
    CosAmp[0][12] = AsymmFit->GetParameter(2);
    CosAmpErr[0][12] = AsymmFit->GetParError(2);

    ///////////////////////////////////////////
    //////////////////  CM2  //////////////////
    ///////////////////////////////////////////

    PhiSc315AsymmCM2 = PhiSc315NegHelCM2->GetAsymmetry(PhiSc315PosHelCM2);
    PhiSc315AsymmCM2->SetName("PhiSc315AsymmCM2");
    PhiSc315AsymmCM2->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 300-330MeV (Cos #theta_{CM} 0.75-0.5)")
    PhiSc315AsymmCM2->Fit("SinFit", "Q");
    InitialSinAmp[1][0] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[1][0] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc315AsymmCM2->Fit("AsymmFit", "Q");
    Offset[1][0] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[1][0] = AsymmFit->GetParError(0);
    SinAmp[1][0] = AsymmFit->GetParameter(1);
    SinAmpErr[1][0] = AsymmFit->GetParError(1);
    CosAmp[1][0] = AsymmFit->GetParameter(2);
    CosAmpErr[1][0] = AsymmFit->GetParError(2);

    PhiSc345AsymmCM2 = PhiSc345NegHelCM2->GetAsymmetry(PhiSc345PosHelCM2);
    PhiSc345AsymmCM2->SetName("PhiSc345AsymmCM2");
    PhiSc345AsymmCM2->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 330-360MeV (Cos #theta_{CM} 0.75-0.5)")
    PhiSc345AsymmCM2->Fit("SinFit", "Q");
    InitialSinAmp[1][1] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[1][1] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc345AsymmCM2->Fit("AsymmFit", "Q");
    Offset[1][1] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[1][1] = AsymmFit->GetParError(0);
    SinAmp[1][1] = AsymmFit->GetParameter(1);
    SinAmpErr[1][1] = AsymmFit->GetParError(1);
    CosAmp[1][1] = AsymmFit->GetParameter(2);
    CosAmpErr[1][1] = AsymmFit->GetParError(2);

    PhiSc375AsymmCM2 = PhiSc375NegHelCM2->GetAsymmetry(PhiSc375PosHelCM2);
    PhiSc375AsymmCM2->SetName("PhiSc375AsymmCM2");
    PhiSc375AsymmCM2->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 360-390MeV (Cos #theta_{CM} 0.75-0.5)")
    PhiSc375AsymmCM2->Fit("SinFit", "Q");
    InitialSinAmp[1][2] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[1][2] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc375AsymmCM2->Fit("AsymmFit", "Q");
    Offset[1][2] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[1][2] = AsymmFit->GetParError(0);
    SinAmp[1][2] = AsymmFit->GetParameter(1);
    SinAmpErr[1][2] = AsymmFit->GetParError(1);
    CosAmp[1][2] = AsymmFit->GetParameter(2);
    CosAmpErr[1][2] = AsymmFit->GetParError(2);

    PhiSc405AsymmCM2 = PhiSc405NegHelCM2->GetAsymmetry(PhiSc405PosHelCM2);
    PhiSc405AsymmCM2->SetName("PhiSc405AsymmCM2");
    PhiSc405AsymmCM2->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 390-420MeV (Cos #theta_{CM} 0.75-0.5)")
    PhiSc405AsymmCM2->Fit("SinFit", "Q");
    InitialSinAmp[1][3] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[1][3] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc405AsymmCM2->Fit("AsymmFit", "Q");
    Offset[1][3] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[1][3] = AsymmFit->GetParError(0);
    SinAmp[1][3] = AsymmFit->GetParameter(1);
    SinAmpErr[1][3] = AsymmFit->GetParError(1);
    CosAmp[1][3] = AsymmFit->GetParameter(2);
    CosAmpErr[1][3] = AsymmFit->GetParError(2);

    PhiSc435AsymmCM2 = PhiSc435NegHelCM2->GetAsymmetry(PhiSc435PosHelCM2);
    PhiSc435AsymmCM2->SetName("PhiSc435AsymmCM2");
    PhiSc435AsymmCM2->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 420-450MeV (Cos #theta_{CM} 0.75-0.5)")
    PhiSc435AsymmCM2->Fit("SinFit", "Q");
    InitialSinAmp[1][4] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[1][4] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc435AsymmCM2->Fit("AsymmFit", "Q");
    Offset[1][4] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[1][4] = AsymmFit->GetParError(0);
    SinAmp[1][4] = AsymmFit->GetParameter(1);
    SinAmpErr[1][4] = AsymmFit->GetParError(1);
    CosAmp[1][4] = AsymmFit->GetParameter(2);
    CosAmpErr[1][4] = AsymmFit->GetParError(2);

    PhiSc465AsymmCM2 = PhiSc465NegHelCM2->GetAsymmetry(PhiSc465PosHelCM2);
    PhiSc465AsymmCM2->SetName("PhiSc465AsymmCM2");
    PhiSc465AsymmCM2->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 450-480MeV (Cos #theta_{CM} 0.75-0.5)")
    PhiSc465AsymmCM2->Fit("SinFit", "Q");
    InitialSinAmp[1][5] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[1][5] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc465AsymmCM2->Fit("AsymmFit", "Q");
    Offset[1][5] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[1][5] = AsymmFit->GetParError(0);
    SinAmp[1][5] = AsymmFit->GetParameter(1);
    SinAmpErr[1][5] = AsymmFit->GetParError(1);
    CosAmp[1][5] = AsymmFit->GetParameter(2);
    CosAmpErr[1][5] = AsymmFit->GetParError(2);

    PhiSc495AsymmCM2 = PhiSc495NegHelCM2->GetAsymmetry(PhiSc495PosHelCM2);
    PhiSc495AsymmCM2->SetName("PhiSc495AsymmCM2");
    PhiSc495AsymmCM2->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 480-510MeV (Cos #theta_{CM} 0.75-0.5)")
    PhiSc495AsymmCM2->Fit("SinFit", "Q");
    InitialSinAmp[1][6] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[1][6] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc495AsymmCM2->Fit("AsymmFit", "Q");
    Offset[1][6] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[1][6] = AsymmFit->GetParError(0);
    SinAmp[1][6] = AsymmFit->GetParameter(1);
    SinAmpErr[1][6] = AsymmFit->GetParError(1);
    CosAmp[1][6] = AsymmFit->GetParameter(2);
    CosAmpErr[1][6] = AsymmFit->GetParError(2);

    PhiSc525AsymmCM2 = PhiSc525NegHelCM2->GetAsymmetry(PhiSc525PosHelCM2);
    PhiSc525AsymmCM2->SetName("PhiSc525AsymmCM2");
    PhiSc525AsymmCM2->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 510-5400MeV (Cos #theta_{CM} 0.75-0.5)")
    PhiSc525AsymmCM2->Fit("SinFit", "Q");
    InitialSinAmp[1][7] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[1][7] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc525AsymmCM2->Fit("AsymmFit", "Q");
    Offset[1][7] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[1][7] = AsymmFit->GetParError(0);
    SinAmp[1][7] = AsymmFit->GetParameter(1);
    SinAmpErr[1][7] = AsymmFit->GetParError(1);
    CosAmp[1][7] = AsymmFit->GetParameter(2);
    CosAmpErr[1][7] = AsymmFit->GetParError(2);

    PhiSc555AsymmCM2 = PhiSc555NegHelCM2->GetAsymmetry(PhiSc555PosHelCM2);
    PhiSc555AsymmCM2->SetName("PhiSc555AsymmCM2");
    PhiSc555AsymmCM2->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 540-570MeV (Cos #theta_{CM} 0.75-0.5)")
    PhiSc555AsymmCM2->Fit("SinFit", "Q");
    InitialSinAmp[1][8] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[1][8] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc555AsymmCM2->Fit("AsymmFit", "Q");
    Offset[1][8] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[1][8] = AsymmFit->GetParError(0);
    SinAmp[1][8] = AsymmFit->GetParameter(1);
    SinAmpErr[1][8] = AsymmFit->GetParError(1);
    CosAmp[1][8] = AsymmFit->GetParameter(2);
    CosAmpErr[1][8] = AsymmFit->GetParError(2);

    PhiSc585AsymmCM2 = PhiSc585NegHelCM2->GetAsymmetry(PhiSc585PosHelCM2);
    PhiSc585AsymmCM2->SetName("PhiSc585AsymmCM2");
    PhiSc585AsymmCM2->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 570-600MeV (Cos #theta_{CM} 0.75-0.5)")
    PhiSc585AsymmCM2->Fit("SinFit", "Q");
    InitialSinAmp[1][9] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[1][9] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc585AsymmCM2->Fit("AsymmFit", "Q");
    Offset[1][9] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[1][9] = AsymmFit->GetParError(0);
    SinAmp[1][9] = AsymmFit->GetParameter(1);
    SinAmpErr[1][9] = AsymmFit->GetParError(1);
    CosAmp[1][9] = AsymmFit->GetParameter(2);
    CosAmpErr[1][9] = AsymmFit->GetParError(2);

    PhiSc615AsymmCM2 = PhiSc615NegHelCM2->GetAsymmetry(PhiSc615PosHelCM2);
    PhiSc615AsymmCM2->SetName("PhiSc615AsymmCM2");
    PhiSc615AsymmCM2->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 600-630MeV (Cos #theta_{CM} 0.75-0.5)")
    PhiSc615AsymmCM2->Fit("SinFit", "Q");
    InitialSinAmp[1][10] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[1][10] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc615AsymmCM2->Fit("AsymmFit", "Q");
    Offset[1][10] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[1][10] = AsymmFit->GetParError(0);
    SinAmp[1][10] = AsymmFit->GetParameter(1);
    SinAmpErr[1][10] = AsymmFit->GetParError(1);
    CosAmp[1][10] = AsymmFit->GetParameter(2);
    CosAmpErr[1][10] = AsymmFit->GetParError(2);

    PhiSc645AsymmCM2 = PhiSc645NegHelCM2->GetAsymmetry(PhiSc645PosHelCM2);
    PhiSc645AsymmCM2->SetName("PhiSc645AsymmCM2");
    PhiSc645AsymmCM2->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 630-660MeV (Cos #theta_{CM} 0.75-0.5)")
    PhiSc645AsymmCM2->Fit("SinFit", "Q");
    InitialSinAmp[1][11] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[1][11] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc645AsymmCM2->Fit("AsymmFit", "Q");
    Offset[1][11] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[1][11] = AsymmFit->GetParError(0);
    SinAmp[1][11] = AsymmFit->GetParameter(1);
    SinAmpErr[1][11] = AsymmFit->GetParError(1);
    CosAmp[1][11] = AsymmFit->GetParameter(2);
    CosAmpErr[1][11] = AsymmFit->GetParError(2);

    PhiSc675AsymmCM2 = PhiSc675NegHelCM2->GetAsymmetry(PhiSc675PosHelCM2);
    PhiSc675AsymmCM2->SetName("PhiSc675AsymmCM2");
    PhiSc675AsymmCM2->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 660-690MeV (Cos #theta_{CM} 0.75-0.5)")
    PhiSc675AsymmCM2->Fit("SinFit", "Q");
    InitialSinAmp[1][12] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[1][12] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc675AsymmCM2->Fit("AsymmFit", "Q");
    Offset[1][12] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[1][12] = AsymmFit->GetParError(0);
    SinAmp[1][12] = AsymmFit->GetParameter(1);
    SinAmpErr[1][12] = AsymmFit->GetParError(1);
    CosAmp[1][12] = AsymmFit->GetParameter(2);
    CosAmpErr[1][12] = AsymmFit->GetParError(2);

    ///////////////////////////////////////////
    //////////////////  CM3  //////////////////
    ///////////////////////////////////////////

    PhiSc315AsymmCM3 = PhiSc315NegHelCM3->GetAsymmetry(PhiSc315PosHelCM3);
    PhiSc315AsymmCM3->SetName("PhiSc315AsymmCM3");
    PhiSc315AsymmCM3->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 300-330MeV (Cos #theta_{CM} 0.5-0.25)")
    PhiSc315AsymmCM3->Fit("SinFit", "Q");
    InitialSinAmp[2][0] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[2][0] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc315AsymmCM3->Fit("AsymmFit", "Q");
    Offset[2][0] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[2][0] = AsymmFit->GetParError(0);
    SinAmp[2][0] = AsymmFit->GetParameter(1);
    SinAmpErr[2][0] = AsymmFit->GetParError(1);
    CosAmp[2][0] = AsymmFit->GetParameter(2);
    CosAmpErr[2][0] = AsymmFit->GetParError(2);

    PhiSc345AsymmCM3 = PhiSc345NegHelCM3->GetAsymmetry(PhiSc345PosHelCM3);
    PhiSc345AsymmCM3->SetName("PhiSc345AsymmCM3");
    PhiSc345AsymmCM3->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 330-360MeV (Cos #theta_{CM} 0.5-0.25)")
    PhiSc345AsymmCM3->Fit("SinFit", "Q");
    InitialSinAmp[2][1] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[2][1] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc345AsymmCM3->Fit("AsymmFit", "Q");
    Offset[2][1] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[2][1] = AsymmFit->GetParError(0);
    SinAmp[2][1] = AsymmFit->GetParameter(1);
    SinAmpErr[2][1] = AsymmFit->GetParError(1);
    CosAmp[2][1] = AsymmFit->GetParameter(2);
    CosAmpErr[2][1] = AsymmFit->GetParError(2);

    PhiSc375AsymmCM3 = PhiSc375NegHelCM3->GetAsymmetry(PhiSc375PosHelCM3);
    PhiSc375AsymmCM3->SetName("PhiSc375AsymmCM3");
    PhiSc375AsymmCM3->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 360-390MeV (Cos #theta_{CM} 0.5-0.25)")
    PhiSc375AsymmCM3->Fit("SinFit", "Q");
    InitialSinAmp[2][2] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[2][2] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc375AsymmCM3->Fit("AsymmFit", "Q");
    Offset[2][2] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[2][2] = AsymmFit->GetParError(0);
    SinAmp[2][2] = AsymmFit->GetParameter(1);
    SinAmpErr[2][2] = AsymmFit->GetParError(1);
    CosAmp[2][2] = AsymmFit->GetParameter(2);
    CosAmpErr[2][2] = AsymmFit->GetParError(2);

    PhiSc405AsymmCM3 = PhiSc405NegHelCM3->GetAsymmetry(PhiSc405PosHelCM3);
    PhiSc405AsymmCM3->SetName("PhiSc405AsymmCM3");
    PhiSc405AsymmCM3->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 390-420MeV (Cos #theta_{CM} 0.5-0.25)")
    PhiSc405AsymmCM3->Fit("SinFit", "Q");
    InitialSinAmp[2][3] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[2][3] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc405AsymmCM3->Fit("AsymmFit", "Q");
    Offset[2][3] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[2][3] = AsymmFit->GetParError(0);
    SinAmp[2][3] = AsymmFit->GetParameter(1);
    SinAmpErr[2][3] = AsymmFit->GetParError(1);
    CosAmp[2][3] = AsymmFit->GetParameter(2);
    CosAmpErr[2][3] = AsymmFit->GetParError(2);

    PhiSc435AsymmCM3 = PhiSc435NegHelCM3->GetAsymmetry(PhiSc435PosHelCM3);
    PhiSc435AsymmCM3->SetName("PhiSc435AsymmCM3");
    PhiSc435AsymmCM3->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 420-450MeV (Cos #theta_{CM} 0.5-0.25)")
    PhiSc435AsymmCM3->Fit("SinFit", "Q");
    InitialSinAmp[2][4] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[2][4] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc435AsymmCM3->Fit("AsymmFit", "Q");
    Offset[2][4] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[2][4] = AsymmFit->GetParError(0);
    SinAmp[2][4] = AsymmFit->GetParameter(1);
    SinAmpErr[2][4] = AsymmFit->GetParError(1);
    CosAmp[2][4] = AsymmFit->GetParameter(2);
    CosAmpErr[2][4] = AsymmFit->GetParError(2);

    PhiSc465AsymmCM3 = PhiSc465NegHelCM3->GetAsymmetry(PhiSc465PosHelCM3);
    PhiSc465AsymmCM3->SetName("PhiSc465AsymmCM3");
    PhiSc465AsymmCM3->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 450-480MeV (Cos #theta_{CM} 0.5-0.25)")
    PhiSc465AsymmCM3->Fit("SinFit", "Q");
    InitialSinAmp[2][5] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[2][5] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc465AsymmCM3->Fit("AsymmFit", "Q");
    Offset[2][5] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[2][5] = AsymmFit->GetParError(0);
    SinAmp[2][5] = AsymmFit->GetParameter(1);
    SinAmpErr[2][5] = AsymmFit->GetParError(1);
    CosAmp[2][5] = AsymmFit->GetParameter(2);
    CosAmpErr[2][5] = AsymmFit->GetParError(2);

    PhiSc495AsymmCM3 = PhiSc495NegHelCM3->GetAsymmetry(PhiSc495PosHelCM3);
    PhiSc495AsymmCM3->SetName("PhiSc495AsymmCM3");
    PhiSc495AsymmCM3->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 480-510MeV (Cos #theta_{CM} 0.5-0.25)")
    PhiSc495AsymmCM3->Fit("SinFit", "Q");
    InitialSinAmp[2][6] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[2][6] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc495AsymmCM3->Fit("AsymmFit", "Q");
    Offset[2][6] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[2][6] = AsymmFit->GetParError(0);
    SinAmp[2][6] = AsymmFit->GetParameter(1);
    SinAmpErr[2][6] = AsymmFit->GetParError(1);
    CosAmp[2][6] = AsymmFit->GetParameter(2);
    CosAmpErr[2][6] = AsymmFit->GetParError(2);

    PhiSc525AsymmCM3 = PhiSc525NegHelCM3->GetAsymmetry(PhiSc525PosHelCM3);
    PhiSc525AsymmCM3->SetName("PhiSc525AsymmCM3");
    PhiSc525AsymmCM3->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 510-5400MeV (Cos #theta_{CM} 0.5-0.25)")
    PhiSc525AsymmCM3->Fit("SinFit", "Q");
    InitialSinAmp[2][7] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[2][7] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc525AsymmCM3->Fit("AsymmFit", "Q");
    Offset[2][7] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[2][7] = AsymmFit->GetParError(0);
    SinAmp[2][7] = AsymmFit->GetParameter(1);
    SinAmpErr[2][7] = AsymmFit->GetParError(1);
    CosAmp[2][7] = AsymmFit->GetParameter(2);
    CosAmpErr[2][7] = AsymmFit->GetParError(2);

    PhiSc555AsymmCM3 = PhiSc555NegHelCM3->GetAsymmetry(PhiSc555PosHelCM3);
    PhiSc555AsymmCM3->SetName("PhiSc555AsymmCM3");
    PhiSc555AsymmCM3->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 540-570MeV (Cos #theta_{CM} 0.5-0.25)")
    PhiSc555AsymmCM3->Fit("SinFit", "Q");
    InitialSinAmp[2][8] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[2][8] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc555AsymmCM3->Fit("AsymmFit", "Q");
    Offset[2][8] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[2][8] = AsymmFit->GetParError(0);
    SinAmp[2][8] = AsymmFit->GetParameter(1);
    SinAmpErr[2][8] = AsymmFit->GetParError(1);
    CosAmp[2][8] = AsymmFit->GetParameter(2);
    CosAmpErr[2][8] = AsymmFit->GetParError(2);

    PhiSc585AsymmCM3 = PhiSc585NegHelCM3->GetAsymmetry(PhiSc585PosHelCM3);
    PhiSc585AsymmCM3->SetName("PhiSc585AsymmCM3");
    PhiSc585AsymmCM3->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 570-600MeV (Cos #theta_{CM} 0.5-0.25)")
    PhiSc585AsymmCM3->Fit("SinFit", "Q");
    InitialSinAmp[2][9] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[2][9] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc585AsymmCM3->Fit("AsymmFit", "Q");
    Offset[2][9] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[2][9] = AsymmFit->GetParError(0);
    SinAmp[2][9] = AsymmFit->GetParameter(1);
    SinAmpErr[2][9] = AsymmFit->GetParError(1);
    CosAmp[2][9] = AsymmFit->GetParameter(2);
    CosAmpErr[2][9] = AsymmFit->GetParError(2);

    PhiSc615AsymmCM3 = PhiSc615NegHelCM3->GetAsymmetry(PhiSc615PosHelCM3);
    PhiSc615AsymmCM3->SetName("PhiSc615AsymmCM3");
    PhiSc615AsymmCM3->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 600-630MeV (Cos #theta_{CM} 0.5-0.25)")
    PhiSc615AsymmCM3->Fit("SinFit", "Q");
    InitialSinAmp[2][10] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[2][10] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc615AsymmCM3->Fit("AsymmFit", "Q");
    Offset[2][10] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[2][10] = AsymmFit->GetParError(0);
    SinAmp[2][10] = AsymmFit->GetParameter(1);
    SinAmpErr[2][10] = AsymmFit->GetParError(1);
    CosAmp[2][10] = AsymmFit->GetParameter(2);
    CosAmpErr[2][10] = AsymmFit->GetParError(2);

    PhiSc645AsymmCM3 = PhiSc645NegHelCM3->GetAsymmetry(PhiSc645PosHelCM3);
    PhiSc645AsymmCM3->SetName("PhiSc645AsymmCM3");
    PhiSc645AsymmCM3->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 630-660MeV (Cos #theta_{CM} 0.5-0.25)")
    PhiSc645AsymmCM3->Fit("SinFit", "Q");
    InitialSinAmp[2][11] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[2][11] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc645AsymmCM3->Fit("AsymmFit", "Q");
    Offset[2][11] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[2][11] = AsymmFit->GetParError(0);
    SinAmp[2][11] = AsymmFit->GetParameter(1);
    SinAmpErr[2][11] = AsymmFit->GetParError(1);
    CosAmp[2][11] = AsymmFit->GetParameter(2);
    CosAmpErr[2][11] = AsymmFit->GetParError(2);

    PhiSc675AsymmCM3 = PhiSc675NegHelCM3->GetAsymmetry(PhiSc675PosHelCM3);
    PhiSc675AsymmCM3->SetName("PhiSc675AsymmCM3");
    PhiSc675AsymmCM3->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 660-690MeV (Cos #theta_{CM} 0.5-0.25)")
    PhiSc675AsymmCM3->Fit("SinFit", "Q");
    InitialSinAmp[2][12] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[2][12] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc675AsymmCM3->Fit("AsymmFit", "Q");
    Offset[2][12] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[2][12] = AsymmFit->GetParError(0);
    SinAmp[2][12] = AsymmFit->GetParameter(1);
    SinAmpErr[2][12] = AsymmFit->GetParError(1);
    CosAmp[2][12] = AsymmFit->GetParameter(2);
    CosAmpErr[2][12] = AsymmFit->GetParError(2);

    ///////////////////////////////////////////
    //////////////////  CM4  //////////////////
    ///////////////////////////////////////////

    PhiSc315AsymmCM4 = PhiSc315NegHelCM4->GetAsymmetry(PhiSc315PosHelCM4);
    PhiSc315AsymmCM4->SetName("PhiSc315AsymmCM4");
    PhiSc315AsymmCM4->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 300-330MeV (Cos #theta_{CM} 0.25-0.0)")
    PhiSc315AsymmCM4->Fit("SinFit", "Q");
    InitialSinAmp[3][0] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[3][0] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc315AsymmCM4->Fit("AsymmFit", "Q");
    Offset[3][0] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[3][0] = AsymmFit->GetParError(0);
    SinAmp[3][0] = AsymmFit->GetParameter(1);
    SinAmpErr[3][0] = AsymmFit->GetParError(1);
    CosAmp[3][0] = AsymmFit->GetParameter(2);
    CosAmpErr[3][0] = AsymmFit->GetParError(2);

    PhiSc345AsymmCM4 = PhiSc345NegHelCM4->GetAsymmetry(PhiSc345PosHelCM4);
    PhiSc345AsymmCM4->SetName("PhiSc345AsymmCM4");
    PhiSc345AsymmCM4->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 330-360MeV (Cos #theta_{CM} 0.25-0.0)")
    PhiSc345AsymmCM4->Fit("SinFit", "Q");
    InitialSinAmp[3][1] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[3][1] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc345AsymmCM4->Fit("AsymmFit", "Q");
    Offset[3][1] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[3][1] = AsymmFit->GetParError(0);
    SinAmp[3][1] = AsymmFit->GetParameter(1);
    SinAmpErr[3][1] = AsymmFit->GetParError(1);
    CosAmp[3][1] = AsymmFit->GetParameter(2);
    CosAmpErr[3][1] = AsymmFit->GetParError(2);

    PhiSc375AsymmCM4 = PhiSc375NegHelCM4->GetAsymmetry(PhiSc375PosHelCM4);
    PhiSc375AsymmCM4->SetName("PhiSc375AsymmCM4");
    PhiSc375AsymmCM4->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 360-390MeV (Cos #theta_{CM} 0.25-0.0)")
    PhiSc375AsymmCM4->Fit("SinFit", "Q");
    InitialSinAmp[3][2] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[3][2] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc375AsymmCM4->Fit("AsymmFit", "Q");
    Offset[3][2] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[3][2] = AsymmFit->GetParError(0);
    SinAmp[3][2] = AsymmFit->GetParameter(1);
    SinAmpErr[3][2] = AsymmFit->GetParError(1);
    CosAmp[3][2] = AsymmFit->GetParameter(2);
    CosAmpErr[3][2] = AsymmFit->GetParError(2);

    PhiSc405AsymmCM4 = PhiSc405NegHelCM4->GetAsymmetry(PhiSc405PosHelCM4);
    PhiSc405AsymmCM4->SetName("PhiSc405AsymmCM4");
    PhiSc405AsymmCM4->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 390-420MeV (Cos #theta_{CM} 0.25-0.0)")
    PhiSc405AsymmCM4->Fit("SinFit", "Q");
    InitialSinAmp[3][3] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[3][3] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc405AsymmCM4->Fit("AsymmFit", "Q");
    Offset[3][3] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[3][3] = AsymmFit->GetParError(0);
    SinAmp[3][3] = AsymmFit->GetParameter(1);
    SinAmpErr[3][3] = AsymmFit->GetParError(1);
    CosAmp[3][3] = AsymmFit->GetParameter(2);
    CosAmpErr[3][3] = AsymmFit->GetParError(2);

    PhiSc435AsymmCM4 = PhiSc435NegHelCM4->GetAsymmetry(PhiSc435PosHelCM4);
    PhiSc435AsymmCM4->SetName("PhiSc435AsymmCM4");
    PhiSc435AsymmCM4->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 420-450MeV (Cos #theta_{CM} 0.25-0.0)")
    PhiSc435AsymmCM4->Fit("SinFit", "Q");
    InitialSinAmp[3][4] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[3][4] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc435AsymmCM4->Fit("AsymmFit", "Q");
    Offset[3][4] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[3][4] = AsymmFit->GetParError(0);
    SinAmp[3][4] = AsymmFit->GetParameter(1);
    SinAmpErr[3][4] = AsymmFit->GetParError(1);
    CosAmp[3][4] = AsymmFit->GetParameter(2);
    CosAmpErr[3][4] = AsymmFit->GetParError(2);

    PhiSc465AsymmCM4 = PhiSc465NegHelCM4->GetAsymmetry(PhiSc465PosHelCM4);
    PhiSc465AsymmCM4->SetName("PhiSc465AsymmCM4");
    PhiSc465AsymmCM4->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 450-480MeV (Cos #theta_{CM} 0.25-0.0)")
    PhiSc465AsymmCM4->Fit("SinFit", "Q");
    InitialSinAmp[3][5] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[3][5] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc465AsymmCM4->Fit("AsymmFit", "Q");
    Offset[3][5] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[3][5] = AsymmFit->GetParError(0);
    SinAmp[3][5] = AsymmFit->GetParameter(1);
    SinAmpErr[3][5] = AsymmFit->GetParError(1);
    CosAmp[3][5] = AsymmFit->GetParameter(2);
    CosAmpErr[3][5] = AsymmFit->GetParError(2);

    PhiSc495AsymmCM4 = PhiSc495NegHelCM4->GetAsymmetry(PhiSc495PosHelCM4);
    PhiSc495AsymmCM4->SetName("PhiSc495AsymmCM4");
    PhiSc495AsymmCM4->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 480-510MeV (Cos #theta_{CM} 0.25-0.0)")
    PhiSc495AsymmCM4->Fit("SinFit", "Q");
    InitialSinAmp[3][6] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[3][6] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc495AsymmCM4->Fit("AsymmFit", "Q");
    Offset[3][6] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[3][6] = AsymmFit->GetParError(0);
    SinAmp[3][6] = AsymmFit->GetParameter(1);
    SinAmpErr[3][6] = AsymmFit->GetParError(1);
    CosAmp[3][6] = AsymmFit->GetParameter(2);
    CosAmpErr[3][6] = AsymmFit->GetParError(2);

    PhiSc525AsymmCM4 = PhiSc525NegHelCM4->GetAsymmetry(PhiSc525PosHelCM4);
    PhiSc525AsymmCM4->SetName("PhiSc525AsymmCM4");
    PhiSc525AsymmCM4->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 510-5400MeV (Cos #theta_{CM} 0.25-0.0)")
    PhiSc525AsymmCM4->Fit("SinFit", "Q");
    InitialSinAmp[3][7] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[3][7] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc525AsymmCM4->Fit("AsymmFit", "Q");
    Offset[3][7] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[3][7] = AsymmFit->GetParError(0);
    SinAmp[3][7] = AsymmFit->GetParameter(1);
    SinAmpErr[3][7] = AsymmFit->GetParError(1);
    CosAmp[3][7] = AsymmFit->GetParameter(2);
    CosAmpErr[3][7] = AsymmFit->GetParError(2);

    PhiSc555AsymmCM4 = PhiSc555NegHelCM4->GetAsymmetry(PhiSc555PosHelCM4);
    PhiSc555AsymmCM4->SetName("PhiSc555AsymmCM4");
    PhiSc555AsymmCM4->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 540-570MeV (Cos #theta_{CM} 0.25-0.0)")
    PhiSc555AsymmCM4->Fit("SinFit", "Q");
    InitialSinAmp[3][8] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[3][8] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc555AsymmCM4->Fit("AsymmFit", "Q");
    Offset[3][8] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[3][8] = AsymmFit->GetParError(0);
    SinAmp[3][8] = AsymmFit->GetParameter(1);
    SinAmpErr[3][8] = AsymmFit->GetParError(1);
    CosAmp[3][8] = AsymmFit->GetParameter(2);
    CosAmpErr[3][8] = AsymmFit->GetParError(2);

    PhiSc585AsymmCM4 = PhiSc585NegHelCM4->GetAsymmetry(PhiSc585PosHelCM4);
    PhiSc585AsymmCM4->SetName("PhiSc585AsymmCM4");
    PhiSc585AsymmCM4->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 570-600MeV (Cos #theta_{CM} 0.25-0.0)")
    PhiSc585AsymmCM4->Fit("SinFit", "Q");
    InitialSinAmp[3][9] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[3][9] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc585AsymmCM4->Fit("AsymmFit", "Q");
    Offset[3][9] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[3][9] = AsymmFit->GetParError(0);
    SinAmp[3][9] = AsymmFit->GetParameter(1);
    SinAmpErr[3][9] = AsymmFit->GetParError(1);
    CosAmp[3][9] = AsymmFit->GetParameter(2);
    CosAmpErr[3][9] = AsymmFit->GetParError(2);

    PhiSc615AsymmCM4 = PhiSc615NegHelCM4->GetAsymmetry(PhiSc615PosHelCM4);
    PhiSc615AsymmCM4->SetName("PhiSc615AsymmCM4");
    PhiSc615AsymmCM4->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 600-630MeV (Cos #theta_{CM} 0.25-0.0)")
    PhiSc615AsymmCM4->Fit("SinFit", "Q");
    InitialSinAmp[3][10] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[3][10] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc615AsymmCM4->Fit("AsymmFit", "Q");
    Offset[3][10] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[3][10] = AsymmFit->GetParError(0);
    SinAmp[3][10] = AsymmFit->GetParameter(1);
    SinAmpErr[3][10] = AsymmFit->GetParError(1);
    CosAmp[3][10] = AsymmFit->GetParameter(2);
    CosAmpErr[3][10] = AsymmFit->GetParError(2);

    PhiSc645AsymmCM4 = PhiSc645NegHelCM4->GetAsymmetry(PhiSc645PosHelCM4);
    PhiSc645AsymmCM4->SetName("PhiSc645AsymmCM4");
    PhiSc645AsymmCM4->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 630-660MeV (Cos #theta_{CM} 0.25-0.0)")
    PhiSc645AsymmCM4->Fit("SinFit", "Q");
    InitialSinAmp[3][11] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[3][11] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc645AsymmCM4->Fit("AsymmFit", "Q");
    Offset[3][11] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[3][11] = AsymmFit->GetParError(0);
    SinAmp[3][11] = AsymmFit->GetParameter(1);
    SinAmpErr[3][11] = AsymmFit->GetParError(1);
    CosAmp[3][11] = AsymmFit->GetParameter(2);
    CosAmpErr[3][11] = AsymmFit->GetParError(2);

    PhiSc675AsymmCM4 = PhiSc675NegHelCM4->GetAsymmetry(PhiSc675PosHelCM4);
    PhiSc675AsymmCM4->SetName("PhiSc675AsymmCM4");
    PhiSc675AsymmCM4->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 660-690MeV (Cos #theta_{CM} 0.25-0.0)")
    PhiSc675AsymmCM4->Fit("SinFit", "Q");
    InitialSinAmp[3][12] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[3][12] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc675AsymmCM4->Fit("AsymmFit", "Q");
    Offset[3][12] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[3][12] = AsymmFit->GetParError(0);
    SinAmp[3][12] = AsymmFit->GetParameter(1);
    SinAmpErr[3][12] = AsymmFit->GetParError(1);
    CosAmp[3][12] = AsymmFit->GetParameter(2);
    CosAmpErr[3][12] = AsymmFit->GetParError(2);

    ///////////////////////////////////////////
    //////////////////  CM5  //////////////////
    ///////////////////////////////////////////

    PhiSc315AsymmCM5 = PhiSc315NegHelCM5->GetAsymmetry(PhiSc315PosHelCM5);
    PhiSc315AsymmCM5->SetName("PhiSc315AsymmCM5");
    PhiSc315AsymmCM5->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 300-330MeV (Cos #theta_{CM} 0.0-(-0.25))")
    PhiSc315AsymmCM5->Fit("SinFit", "Q");
    InitialSinAmp[4][0] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[4][0] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc315AsymmCM5->Fit("AsymmFit", "Q");
    Offset[4][0] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[4][0] = AsymmFit->GetParError(0);
    SinAmp[4][0] = AsymmFit->GetParameter(1);
    SinAmpErr[4][0] = AsymmFit->GetParError(1);
    CosAmp[4][0] = AsymmFit->GetParameter(2);
    CosAmpErr[4][0] = AsymmFit->GetParError(2);

    PhiSc345AsymmCM5 = PhiSc345NegHelCM5->GetAsymmetry(PhiSc345PosHelCM5);
    PhiSc345AsymmCM5->SetName("PhiSc345AsymmCM5");
    PhiSc345AsymmCM5->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 330-360MeV (Cos #theta_{CM} 0.0-(-0.25))")
    PhiSc345AsymmCM5->Fit("SinFit", "Q");
    InitialSinAmp[4][1] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[4][1] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc345AsymmCM5->Fit("AsymmFit", "Q");
    Offset[4][1] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[4][1] = AsymmFit->GetParError(0);
    SinAmp[4][1] = AsymmFit->GetParameter(1);
    SinAmpErr[4][1] = AsymmFit->GetParError(1);
    CosAmp[4][1] = AsymmFit->GetParameter(2);
    CosAmpErr[4][1] = AsymmFit->GetParError(2);

    PhiSc375AsymmCM5 = PhiSc375NegHelCM5->GetAsymmetry(PhiSc375PosHelCM5);
    PhiSc375AsymmCM5->SetName("PhiSc375AsymmCM5");
    PhiSc375AsymmCM5->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 360-390MeV (Cos #theta_{CM} 0.0-(-0.25))")
    PhiSc375AsymmCM5->Fit("SinFit", "Q");
    InitialSinAmp[4][2] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[4][2] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc375AsymmCM5->Fit("AsymmFit", "Q");
    Offset[4][2] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[4][2] = AsymmFit->GetParError(0);
    SinAmp[4][2] = AsymmFit->GetParameter(1);
    SinAmpErr[4][2] = AsymmFit->GetParError(1);
    CosAmp[4][2] = AsymmFit->GetParameter(2);
    CosAmpErr[4][2] = AsymmFit->GetParError(2);

    PhiSc405AsymmCM5 = PhiSc405NegHelCM5->GetAsymmetry(PhiSc405PosHelCM5);
    PhiSc405AsymmCM5->SetName("PhiSc405AsymmCM5");
    PhiSc405AsymmCM5->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 390-420MeV (Cos #theta_{CM} 0.0-(-0.25))")
    PhiSc405AsymmCM5->Fit("SinFit", "Q");
    InitialSinAmp[4][3] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[4][3] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc405AsymmCM5->Fit("AsymmFit", "Q");
    Offset[4][3] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[4][3] = AsymmFit->GetParError(0);
    SinAmp[4][3] = AsymmFit->GetParameter(1);
    SinAmpErr[4][3] = AsymmFit->GetParError(1);
    CosAmp[4][3] = AsymmFit->GetParameter(2);
    CosAmpErr[4][3] = AsymmFit->GetParError(2);

    PhiSc435AsymmCM5 = PhiSc435NegHelCM5->GetAsymmetry(PhiSc435PosHelCM5);
    PhiSc435AsymmCM5->SetName("PhiSc435AsymmCM5");
    PhiSc435AsymmCM5->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 420-450MeV (Cos #theta_{CM} 0.0-(-0.25))")
    PhiSc435AsymmCM5->Fit("SinFit", "Q");
    InitialSinAmp[4][4] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[4][4] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc435AsymmCM5->Fit("AsymmFit", "Q");
    Offset[4][4] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[4][4] = AsymmFit->GetParError(0);
    SinAmp[4][4] = AsymmFit->GetParameter(1);
    SinAmpErr[4][4] = AsymmFit->GetParError(1);
    CosAmp[4][4] = AsymmFit->GetParameter(2);
    CosAmpErr[4][4] = AsymmFit->GetParError(2);

    PhiSc465AsymmCM5 = PhiSc465NegHelCM5->GetAsymmetry(PhiSc465PosHelCM5);
    PhiSc465AsymmCM5->SetName("PhiSc465AsymmCM5");
    PhiSc465AsymmCM5->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 450-480MeV (Cos #theta_{CM} 0.0-(-0.25))")
    PhiSc465AsymmCM5->Fit("SinFit", "Q");
    InitialSinAmp[4][5] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[4][5] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc465AsymmCM5->Fit("AsymmFit", "Q");
    Offset[4][5] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[4][5] = AsymmFit->GetParError(0);
    SinAmp[4][5] = AsymmFit->GetParameter(1);
    SinAmpErr[4][5] = AsymmFit->GetParError(1);
    CosAmp[4][5] = AsymmFit->GetParameter(2);
    CosAmpErr[4][5] = AsymmFit->GetParError(2);

    PhiSc495AsymmCM5 = PhiSc495NegHelCM5->GetAsymmetry(PhiSc495PosHelCM5);
    PhiSc495AsymmCM5->SetName("PhiSc495AsymmCM5");
    PhiSc495AsymmCM5->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 480-510MeV (Cos #theta_{CM} 0.0-(-0.25))")
    PhiSc495AsymmCM5->Fit("SinFit", "Q");
    InitialSinAmp[4][6] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[4][6] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc495AsymmCM5->Fit("AsymmFit", "Q");
    Offset[4][6] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[4][6] = AsymmFit->GetParError(0);
    SinAmp[4][6] = AsymmFit->GetParameter(1);
    SinAmpErr[4][6] = AsymmFit->GetParError(1);
    CosAmp[4][6] = AsymmFit->GetParameter(2);
    CosAmpErr[4][6] = AsymmFit->GetParError(2);

    PhiSc525AsymmCM5 = PhiSc525NegHelCM5->GetAsymmetry(PhiSc525PosHelCM5);
    PhiSc525AsymmCM5->SetName("PhiSc525AsymmCM5");
    PhiSc525AsymmCM5->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 510-5400MeV (Cos #theta_{CM} 0.0-(-0.25))")
    PhiSc525AsymmCM5->Fit("SinFit", "Q");
    InitialSinAmp[4][7] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[4][7] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc525AsymmCM5->Fit("AsymmFit", "Q");
    Offset[4][7] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[4][7] = AsymmFit->GetParError(0);
    SinAmp[4][7] = AsymmFit->GetParameter(1);
    SinAmpErr[4][7] = AsymmFit->GetParError(1);
    CosAmp[4][7] = AsymmFit->GetParameter(2);
    CosAmpErr[4][7] = AsymmFit->GetParError(2);

    PhiSc555AsymmCM5 = PhiSc555NegHelCM5->GetAsymmetry(PhiSc555PosHelCM5);
    PhiSc555AsymmCM5->SetName("PhiSc555AsymmCM5");
    PhiSc555AsymmCM5->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 540-570MeV (Cos #theta_{CM} 0.0-(-0.25))")
    PhiSc555AsymmCM5->Fit("SinFit", "Q");
    InitialSinAmp[4][8] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[4][8] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc555AsymmCM5->Fit("AsymmFit", "Q");
    Offset[4][8] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[4][8] = AsymmFit->GetParError(0);
    SinAmp[4][8] = AsymmFit->GetParameter(1);
    SinAmpErr[4][8] = AsymmFit->GetParError(1);
    CosAmp[4][8] = AsymmFit->GetParameter(2);
    CosAmpErr[4][8] = AsymmFit->GetParError(2);

    PhiSc585AsymmCM5 = PhiSc585NegHelCM5->GetAsymmetry(PhiSc585PosHelCM5);
    PhiSc585AsymmCM5->SetName("PhiSc585AsymmCM5");
    PhiSc585AsymmCM5->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 570-600MeV (Cos #theta_{CM} 0.0-(-0.25))")
    PhiSc585AsymmCM5->Fit("SinFit", "Q");
    InitialSinAmp[4][9] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[4][9] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc585AsymmCM5->Fit("AsymmFit", "Q");
    Offset[4][9] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[4][9] = AsymmFit->GetParError(0);
    SinAmp[4][9] = AsymmFit->GetParameter(1);
    SinAmpErr[4][9] = AsymmFit->GetParError(1);
    CosAmp[4][9] = AsymmFit->GetParameter(2);
    CosAmpErr[4][9] = AsymmFit->GetParError(2);

    PhiSc615AsymmCM5 = PhiSc615NegHelCM5->GetAsymmetry(PhiSc615PosHelCM5);
    PhiSc615AsymmCM5->SetName("PhiSc615AsymmCM5");
    PhiSc615AsymmCM5->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 600-630MeV (Cos #theta_{CM} 0.0-(-0.25))")
    PhiSc615AsymmCM5->Fit("SinFit", "Q");
    InitialSinAmp[4][10] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[4][10] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc615AsymmCM5->Fit("AsymmFit", "Q");
    Offset[4][10] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[4][10] = AsymmFit->GetParError(0);
    SinAmp[4][10] = AsymmFit->GetParameter(1);
    SinAmpErr[4][10] = AsymmFit->GetParError(1);
    CosAmp[4][10] = AsymmFit->GetParameter(2);
    CosAmpErr[4][10] = AsymmFit->GetParError(2);

    PhiSc645AsymmCM5 = PhiSc645NegHelCM5->GetAsymmetry(PhiSc645PosHelCM5);
    PhiSc645AsymmCM5->SetName("PhiSc645AsymmCM5");
    PhiSc645AsymmCM5->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 630-660MeV (Cos #theta_{CM} 0.0-(-0.25))")
    PhiSc645AsymmCM5->Fit("SinFit", "Q");
    InitialSinAmp[4][11] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[4][11] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc645AsymmCM5->Fit("AsymmFit", "Q");
    Offset[4][11] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[4][11] = AsymmFit->GetParError(0);
    SinAmp[4][11] = AsymmFit->GetParameter(1);
    SinAmpErr[4][11] = AsymmFit->GetParError(1);
    CosAmp[4][11] = AsymmFit->GetParameter(2);
    CosAmpErr[4][11] = AsymmFit->GetParError(2);

    PhiSc675AsymmCM5 = PhiSc675NegHelCM5->GetAsymmetry(PhiSc675PosHelCM5);
    PhiSc675AsymmCM5->SetName("PhiSc675AsymmCM5");
    PhiSc675AsymmCM5->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 660-690MeV (Cos #theta_{CM} 0.0-(-0.25))")
    PhiSc675AsymmCM5->Fit("SinFit", "Q");
    InitialSinAmp[4][12] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[4][12] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc675AsymmCM5->Fit("AsymmFit", "Q");
    Offset[4][12] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[4][12] = AsymmFit->GetParError(0);
    SinAmp[4][12] = AsymmFit->GetParameter(1);
    SinAmpErr[4][12] = AsymmFit->GetParError(1);
    CosAmp[4][12] = AsymmFit->GetParameter(2);
    CosAmpErr[4][12] = AsymmFit->GetParError(2);

    ///////////////////////////////////////////
    //////////////////  CM6  //////////////////
    ///////////////////////////////////////////

    PhiSc315AsymmCM6 = PhiSc315NegHelCM6->GetAsymmetry(PhiSc315PosHelCM6);
    PhiSc315AsymmCM6->SetName("PhiSc315AsymmCM6");
    PhiSc315AsymmCM6->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 300-330MeV (Cos #theta_{CM} -0.25-(-0.5))")
    PhiSc315AsymmCM6->Fit("SinFit", "Q");
    InitialSinAmp[5][0] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[5][0] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc315AsymmCM6->Fit("AsymmFit", "Q");
    Offset[5][0] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[5][0] = AsymmFit->GetParError(0);
    SinAmp[5][0] = AsymmFit->GetParameter(1);
    SinAmpErr[5][0] = AsymmFit->GetParError(1);
    CosAmp[5][0] = AsymmFit->GetParameter(2);
    CosAmpErr[5][0] = AsymmFit->GetParError(2);

    PhiSc345AsymmCM6 = PhiSc345NegHelCM6->GetAsymmetry(PhiSc345PosHelCM6);
    PhiSc345AsymmCM6->SetName("PhiSc345AsymmCM6");
    PhiSc345AsymmCM6->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 330-360MeV (Cos #theta_{CM} -0.25-(-0.5))")
    PhiSc345AsymmCM6->Fit("SinFit", "Q");
    InitialSinAmp[5][1] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[5][1] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc345AsymmCM6->Fit("AsymmFit", "Q");
    Offset[5][1] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[5][1] = AsymmFit->GetParError(0);
    SinAmp[5][1] = AsymmFit->GetParameter(1);
    SinAmpErr[5][1] = AsymmFit->GetParError(1);
    CosAmp[5][1] = AsymmFit->GetParameter(2);
    CosAmpErr[5][1] = AsymmFit->GetParError(2);

    PhiSc375AsymmCM6 = PhiSc375NegHelCM6->GetAsymmetry(PhiSc375PosHelCM6);
    PhiSc375AsymmCM6->SetName("PhiSc375AsymmCM6");
    PhiSc375AsymmCM6->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 360-390MeV (Cos #theta_{CM} -0.25-(-0.5))")
    PhiSc375AsymmCM6->Fit("SinFit", "Q");
    InitialSinAmp[5][2] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[5][2] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc375AsymmCM6->Fit("AsymmFit", "Q");
    Offset[5][2] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[5][2] = AsymmFit->GetParError(0);
    SinAmp[5][2] = AsymmFit->GetParameter(1);
    SinAmpErr[5][2] = AsymmFit->GetParError(1);
    CosAmp[5][2] = AsymmFit->GetParameter(2);
    CosAmpErr[5][2] = AsymmFit->GetParError(2);

    PhiSc405AsymmCM6 = PhiSc405NegHelCM6->GetAsymmetry(PhiSc405PosHelCM6);
    PhiSc405AsymmCM6->SetName("PhiSc405AsymmCM6");
    PhiSc405AsymmCM6->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 390-420MeV (Cos #theta_{CM} -0.25-(-0.5))")
    PhiSc405AsymmCM6->Fit("SinFit", "Q");
    InitialSinAmp[5][3] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[5][3] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc405AsymmCM6->Fit("AsymmFit", "Q");
    Offset[5][3] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[5][3] = AsymmFit->GetParError(0);
    SinAmp[5][3] = AsymmFit->GetParameter(1);
    SinAmpErr[5][3] = AsymmFit->GetParError(1);
    CosAmp[5][3] = AsymmFit->GetParameter(2);
    CosAmpErr[5][3] = AsymmFit->GetParError(2);

    PhiSc435AsymmCM6 = PhiSc435NegHelCM6->GetAsymmetry(PhiSc435PosHelCM6);
    PhiSc435AsymmCM6->SetName("PhiSc435AsymmCM6");
    PhiSc435AsymmCM6->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 420-450MeV (Cos #theta_{CM} -0.25-(-0.5))")
    PhiSc435AsymmCM6->Fit("SinFit", "Q");
    InitialSinAmp[5][4] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[5][4] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc435AsymmCM6->Fit("AsymmFit", "Q");
    Offset[5][4] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[5][4] = AsymmFit->GetParError(0);
    SinAmp[5][4] = AsymmFit->GetParameter(1);
    SinAmpErr[5][4] = AsymmFit->GetParError(1);
    CosAmp[5][4] = AsymmFit->GetParameter(2);
    CosAmpErr[5][4] = AsymmFit->GetParError(2);

    PhiSc465AsymmCM6 = PhiSc465NegHelCM6->GetAsymmetry(PhiSc465PosHelCM6);
    PhiSc465AsymmCM6->SetName("PhiSc465AsymmCM6");
    PhiSc465AsymmCM6->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 450-480MeV (Cos #theta_{CM} -0.25-(-0.5))")
    PhiSc465AsymmCM6->Fit("SinFit", "Q");
    InitialSinAmp[5][5] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[5][5] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc465AsymmCM6->Fit("AsymmFit", "Q");
    Offset[5][5] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[5][5] = AsymmFit->GetParError(0);
    SinAmp[5][5] = AsymmFit->GetParameter(1);
    SinAmpErr[5][5] = AsymmFit->GetParError(1);
    CosAmp[5][5] = AsymmFit->GetParameter(2);
    CosAmpErr[5][5] = AsymmFit->GetParError(2);

    PhiSc495AsymmCM6 = PhiSc495NegHelCM6->GetAsymmetry(PhiSc495PosHelCM6);
    PhiSc495AsymmCM6->SetName("PhiSc495AsymmCM6");
    PhiSc495AsymmCM6->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 480-510MeV (Cos #theta_{CM} -0.25-(-0.5))")
    PhiSc495AsymmCM6->Fit("SinFit", "Q");
    InitialSinAmp[5][6] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[5][6] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc495AsymmCM6->Fit("AsymmFit", "Q");
    Offset[5][6] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[5][6] = AsymmFit->GetParError(0);
    SinAmp[5][6] = AsymmFit->GetParameter(1);
    SinAmpErr[5][6] = AsymmFit->GetParError(1);
    CosAmp[5][6] = AsymmFit->GetParameter(2);
    CosAmpErr[5][6] = AsymmFit->GetParError(2);

    PhiSc525AsymmCM6 = PhiSc525NegHelCM6->GetAsymmetry(PhiSc525PosHelCM6);
    PhiSc525AsymmCM6->SetName("PhiSc525AsymmCM6");
    PhiSc525AsymmCM6->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 510-5400MeV (Cos #theta_{CM} -0.25-(-0.5))")
    PhiSc525AsymmCM6->Fit("SinFit", "Q");
    InitialSinAmp[5][7] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[5][7] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc525AsymmCM6->Fit("AsymmFit", "Q");
    Offset[5][7] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[5][7] = AsymmFit->GetParError(0);
    SinAmp[5][7] = AsymmFit->GetParameter(1);
    SinAmpErr[5][7] = AsymmFit->GetParError(1);
    CosAmp[5][7] = AsymmFit->GetParameter(2);
    CosAmpErr[5][7] = AsymmFit->GetParError(2);

    PhiSc555AsymmCM6 = PhiSc555NegHelCM6->GetAsymmetry(PhiSc555PosHelCM6);
    PhiSc555AsymmCM6->SetName("PhiSc555AsymmCM6");
    PhiSc555AsymmCM6->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 540-570MeV (Cos #theta_{CM} -0.25-(-0.5))")
    PhiSc555AsymmCM6->Fit("SinFit", "Q");
    InitialSinAmp[5][8] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[5][8] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc555AsymmCM6->Fit("AsymmFit", "Q");
    Offset[5][8] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[5][8] = AsymmFit->GetParError(0);
    SinAmp[5][8] = AsymmFit->GetParameter(1);
    SinAmpErr[5][8] = AsymmFit->GetParError(1);
    CosAmp[5][8] = AsymmFit->GetParameter(2);
    CosAmpErr[5][8] = AsymmFit->GetParError(2);

    PhiSc585AsymmCM6 = PhiSc585NegHelCM6->GetAsymmetry(PhiSc585PosHelCM6);
    PhiSc585AsymmCM6->SetName("PhiSc585AsymmCM6");
    PhiSc585AsymmCM6->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 570-600MeV (Cos #theta_{CM} -0.25-(-0.5))")
    PhiSc585AsymmCM6->Fit("SinFit", "Q");
    InitialSinAmp[5][9] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[5][9] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc585AsymmCM6->Fit("AsymmFit", "Q");
    Offset[5][9] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[5][9] = AsymmFit->GetParError(0);
    SinAmp[5][9] = AsymmFit->GetParameter(1);
    SinAmpErr[5][9] = AsymmFit->GetParError(1);
    CosAmp[5][9] = AsymmFit->GetParameter(2);
    CosAmpErr[5][9] = AsymmFit->GetParError(2);

    PhiSc615AsymmCM6 = PhiSc615NegHelCM6->GetAsymmetry(PhiSc615PosHelCM6);
    PhiSc615AsymmCM6->SetName("PhiSc615AsymmCM6");
    PhiSc615AsymmCM6->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 600-630MeV (Cos #theta_{CM} -0.25-(-0.5))")
    PhiSc615AsymmCM6->Fit("SinFit", "Q");
    InitialSinAmp[5][10] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[5][10] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc615AsymmCM6->Fit("AsymmFit", "Q");
    Offset[5][10] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[5][10] = AsymmFit->GetParError(0);
    SinAmp[5][10] = AsymmFit->GetParameter(1);
    SinAmpErr[5][10] = AsymmFit->GetParError(1);
    CosAmp[5][10] = AsymmFit->GetParameter(2);
    CosAmpErr[5][10] = AsymmFit->GetParError(2);

    PhiSc645AsymmCM6 = PhiSc645NegHelCM6->GetAsymmetry(PhiSc645PosHelCM6);
    PhiSc645AsymmCM6->SetName("PhiSc645AsymmCM6");
    PhiSc645AsymmCM6->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 630-660MeV (Cos #theta_{CM} -0.25-(-0.5))")
    PhiSc645AsymmCM6->Fit("SinFit", "Q");
    InitialSinAmp[5][11] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[5][11] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc645AsymmCM6->Fit("AsymmFit", "Q");
    Offset[5][11] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[5][11] = AsymmFit->GetParError(0);
    SinAmp[5][11] = AsymmFit->GetParameter(1);
    SinAmpErr[5][11] = AsymmFit->GetParError(1);
    CosAmp[5][11] = AsymmFit->GetParameter(2);
    CosAmpErr[5][11] = AsymmFit->GetParError(2);

    PhiSc675AsymmCM6 = PhiSc675NegHelCM6->GetAsymmetry(PhiSc675PosHelCM6);
    PhiSc675AsymmCM6->SetName("PhiSc675AsymmCM6");
    PhiSc675AsymmCM6->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 660-690MeV (Cos #theta_{CM} -0.25-(-0.5))")
    PhiSc675AsymmCM6->Fit("SinFit", "Q");
    InitialSinAmp[5][12] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[5][12] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc675AsymmCM6->Fit("AsymmFit", "Q");
    Offset[5][12] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[5][12] = AsymmFit->GetParError(0);
    SinAmp[5][12] = AsymmFit->GetParameter(1);
    SinAmpErr[5][12] = AsymmFit->GetParError(1);
    CosAmp[5][12] = AsymmFit->GetParameter(2);
    CosAmpErr[5][12] = AsymmFit->GetParError(2);

    ///////////////////////////////////////////
    //////////////////  CM7 //////////////////
    ///////////////////////////////////////////

    PhiSc315AsymmCM7 = PhiSc315NegHelCM7->GetAsymmetry(PhiSc315PosHelCM7);
    PhiSc315AsymmCM7->SetName("PhiSc315AsymmCM7");
    PhiSc315AsymmCM7->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 300-330MeV (Cos #theta_{CM} -0.5-(-0.75))")
    PhiSc315AsymmCM7->Fit("SinFit", "Q");
    InitialSinAmp[6][0] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[6][0] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc315AsymmCM7->Fit("AsymmFit", "Q");
    Offset[6][0] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[6][0] = AsymmFit->GetParError(0);
    SinAmp[6][0] = AsymmFit->GetParameter(1);
    SinAmpErr[6][0] = AsymmFit->GetParError(1);
    CosAmp[6][0] = AsymmFit->GetParameter(2);
    CosAmpErr[6][0] = AsymmFit->GetParError(2);

    PhiSc345AsymmCM7 = PhiSc345NegHelCM7->GetAsymmetry(PhiSc345PosHelCM7);
    PhiSc345AsymmCM7->SetName("PhiSc345AsymmCM7");
    PhiSc345AsymmCM7->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 330-360MeV (Cos #theta_{CM} -0.5-(-0.75))")
    PhiSc345AsymmCM7->Fit("SinFit", "Q");
    InitialSinAmp[6][1] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[6][1] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc345AsymmCM7->Fit("AsymmFit", "Q");
    Offset[6][1] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[6][1] = AsymmFit->GetParError(0);
    SinAmp[6][1] = AsymmFit->GetParameter(1);
    SinAmpErr[6][1] = AsymmFit->GetParError(1);
    CosAmp[6][1] = AsymmFit->GetParameter(2);
    CosAmpErr[6][1] = AsymmFit->GetParError(2);

    PhiSc375AsymmCM7 = PhiSc375NegHelCM7->GetAsymmetry(PhiSc375PosHelCM7);
    PhiSc375AsymmCM7->SetName("PhiSc375AsymmCM7");
    PhiSc375AsymmCM7->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 360-390MeV (Cos #theta_{CM} -0.5-(-0.75))")
    PhiSc375AsymmCM7->Fit("SinFit", "Q");
    InitialSinAmp[6][2] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[6][2] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc375AsymmCM7->Fit("AsymmFit", "Q");
    Offset[6][2] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[6][2] = AsymmFit->GetParError(0);
    SinAmp[6][2] = AsymmFit->GetParameter(1);
    SinAmpErr[6][2] = AsymmFit->GetParError(1);
    CosAmp[6][2] = AsymmFit->GetParameter(2);
    CosAmpErr[6][2] = AsymmFit->GetParError(2);

    PhiSc405AsymmCM7 = PhiSc405NegHelCM7->GetAsymmetry(PhiSc405PosHelCM7);
    PhiSc405AsymmCM7->SetName("PhiSc405AsymmCM7");
    PhiSc405AsymmCM7->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 390-420MeV (Cos #theta_{CM} -0.5-(-0.75))")
    PhiSc405AsymmCM7->Fit("SinFit", "Q");
    InitialSinAmp[6][3] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[6][3] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc405AsymmCM7->Fit("AsymmFit", "Q");
    Offset[6][3] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[6][3] = AsymmFit->GetParError(0);
    SinAmp[6][3] = AsymmFit->GetParameter(1);
    SinAmpErr[6][3] = AsymmFit->GetParError(1);
    CosAmp[6][3] = AsymmFit->GetParameter(2);
    CosAmpErr[6][3] = AsymmFit->GetParError(2);

    PhiSc435AsymmCM7 = PhiSc435NegHelCM7->GetAsymmetry(PhiSc435PosHelCM7);
    PhiSc435AsymmCM7->SetName("PhiSc435AsymmCM7");
    PhiSc435AsymmCM7->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 420-450MeV (Cos #theta_{CM} -0.5-(-0.75))")
    PhiSc435AsymmCM7->Fit("SinFit", "Q");
    InitialSinAmp[6][4] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[6][4] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc435AsymmCM7->Fit("AsymmFit", "Q");
    Offset[6][4] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[6][4] = AsymmFit->GetParError(0);
    SinAmp[6][4] = AsymmFit->GetParameter(1);
    SinAmpErr[6][4] = AsymmFit->GetParError(1);
    CosAmp[6][4] = AsymmFit->GetParameter(2);
    CosAmpErr[6][4] = AsymmFit->GetParError(2);

    PhiSc465AsymmCM7 = PhiSc465NegHelCM7->GetAsymmetry(PhiSc465PosHelCM7);
    PhiSc465AsymmCM7->SetName("PhiSc465AsymmCM7");
    PhiSc465AsymmCM7->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 450-480MeV (Cos #theta_{CM} -0.5-(-0.75))")
    PhiSc465AsymmCM7->Fit("SinFit", "Q");
    InitialSinAmp[6][5] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[6][5] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc465AsymmCM7->Fit("AsymmFit", "Q");
    Offset[6][5] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[6][5] = AsymmFit->GetParError(0);
    SinAmp[6][5] = AsymmFit->GetParameter(1);
    SinAmpErr[6][5] = AsymmFit->GetParError(1);
    CosAmp[6][5] = AsymmFit->GetParameter(2);
    CosAmpErr[6][5] = AsymmFit->GetParError(2);

    PhiSc495AsymmCM7 = PhiSc495NegHelCM7->GetAsymmetry(PhiSc495PosHelCM7);
    PhiSc495AsymmCM7->SetName("PhiSc495AsymmCM7");
    PhiSc495AsymmCM7->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 480-510MeV (Cos #theta_{CM} -0.5-(-0.75))")
    PhiSc495AsymmCM7->Fit("SinFit", "Q");
    InitialSinAmp[6][6] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[6][6] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc495AsymmCM7->Fit("AsymmFit", "Q");
    Offset[6][6] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[6][6] = AsymmFit->GetParError(0);
    SinAmp[6][6] = AsymmFit->GetParameter(1);
    SinAmpErr[6][6] = AsymmFit->GetParError(1);
    CosAmp[6][6] = AsymmFit->GetParameter(2);
    CosAmpErr[6][6] = AsymmFit->GetParError(2);

    PhiSc525AsymmCM7 = PhiSc525NegHelCM7->GetAsymmetry(PhiSc525PosHelCM7);
    PhiSc525AsymmCM7->SetName("PhiSc525AsymmCM7");
    PhiSc525AsymmCM7->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 510-5400MeV (Cos #theta_{CM} -0.5-(-0.75))")
    PhiSc525AsymmCM7->Fit("SinFit", "Q");
    InitialSinAmp[6][7] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[6][7] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc525AsymmCM7->Fit("AsymmFit", "Q");
    Offset[6][7] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[6][7] = AsymmFit->GetParError(0);
    SinAmp[6][7] = AsymmFit->GetParameter(1);
    SinAmpErr[6][7] = AsymmFit->GetParError(1);
    CosAmp[6][7] = AsymmFit->GetParameter(2);
    CosAmpErr[6][7] = AsymmFit->GetParError(2);

    PhiSc555AsymmCM7 = PhiSc555NegHelCM7->GetAsymmetry(PhiSc555PosHelCM7);
    PhiSc555AsymmCM7->SetName("PhiSc555AsymmCM7");
    PhiSc555AsymmCM7->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 540-570MeV (Cos #theta_{CM} -0.5-(-0.75))")
    PhiSc555AsymmCM7->Fit("SinFit", "Q");
    InitialSinAmp[6][8] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[6][8] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc555AsymmCM7->Fit("AsymmFit", "Q");
    Offset[6][8] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[6][8] = AsymmFit->GetParError(0);
    SinAmp[6][8] = AsymmFit->GetParameter(1);
    SinAmpErr[6][8] = AsymmFit->GetParError(1);
    CosAmp[6][8] = AsymmFit->GetParameter(2);
    CosAmpErr[6][8] = AsymmFit->GetParError(2);

    PhiSc585AsymmCM7 = PhiSc585NegHelCM7->GetAsymmetry(PhiSc585PosHelCM7);
    PhiSc585AsymmCM7->SetName("PhiSc585AsymmCM7");
    PhiSc585AsymmCM7->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 570-600MeV (Cos #theta_{CM} -0.5-(-0.75))")
    PhiSc585AsymmCM7->Fit("SinFit", "Q");
    InitialSinAmp[6][9] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[6][9] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc585AsymmCM7->Fit("AsymmFit", "Q");
    Offset[6][9] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[6][9] = AsymmFit->GetParError(0);
    SinAmp[6][9] = AsymmFit->GetParameter(1);
    SinAmpErr[6][9] = AsymmFit->GetParError(1);
    CosAmp[6][9] = AsymmFit->GetParameter(2);
    CosAmpErr[6][9] = AsymmFit->GetParError(2);

    PhiSc615AsymmCM7 = PhiSc615NegHelCM7->GetAsymmetry(PhiSc615PosHelCM7);
    PhiSc615AsymmCM7->SetName("PhiSc615AsymmCM7");
    PhiSc615AsymmCM7->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 600-630MeV (Cos #theta_{CM} -0.5-(-0.75))")
    PhiSc615AsymmCM7->Fit("SinFit", "Q");
    InitialSinAmp[6][10] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[6][10] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc615AsymmCM7->Fit("AsymmFit", "Q");
    Offset[6][10] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[6][10] = AsymmFit->GetParError(0);
    SinAmp[6][10] = AsymmFit->GetParameter(1);
    SinAmpErr[6][10] = AsymmFit->GetParError(1);
    CosAmp[6][10] = AsymmFit->GetParameter(2);
    CosAmpErr[6][10] = AsymmFit->GetParError(2);

    PhiSc645AsymmCM7 = PhiSc645NegHelCM7->GetAsymmetry(PhiSc645PosHelCM7);
    PhiSc645AsymmCM7->SetName("PhiSc645AsymmCM7");
    PhiSc645AsymmCM7->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 630-660MeV (Cos #theta_{CM} -0.5-(-0.75))")
    PhiSc645AsymmCM7->Fit("SinFit", "Q");
    InitialSinAmp[6][11] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[6][11] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc645AsymmCM7->Fit("AsymmFit", "Q");
    Offset[6][11] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[6][11] = AsymmFit->GetParError(0);
    SinAmp[6][11] = AsymmFit->GetParameter(1);
    SinAmpErr[6][11] = AsymmFit->GetParError(1);
    CosAmp[6][11] = AsymmFit->GetParameter(2);
    CosAmpErr[6][11] = AsymmFit->GetParError(2);

    PhiSc675AsymmCM7 = PhiSc675NegHelCM7->GetAsymmetry(PhiSc675PosHelCM7);
    PhiSc675AsymmCM7->SetName("PhiSc675AsymmCM7");
    PhiSc675AsymmCM7->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 660-690MeV (Cos #theta_{CM} -0.5-(-0.75))")
    PhiSc675AsymmCM7->Fit("SinFit", "Q");
    InitialSinAmp[6][12] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[6][12] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc675AsymmCM7->Fit("AsymmFit", "Q");
    Offset[6][12] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[6][12] = AsymmFit->GetParError(0);
    SinAmp[6][12] = AsymmFit->GetParameter(1);
    SinAmpErr[6][12] = AsymmFit->GetParError(1);
    CosAmp[6][12] = AsymmFit->GetParameter(2);
    CosAmpErr[6][12] = AsymmFit->GetParError(2);

    ///////////////////////////////////////////
    //////////////////  CM8 //////////////////
    ///////////////////////////////////////////

    PhiSc315AsymmCM8 = PhiSc315NegHelCM8->GetAsymmetry(PhiSc315PosHelCM8);
    PhiSc315AsymmCM8->SetName("PhiSc315AsymmCM8");
    PhiSc315AsymmCM8->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 300-330MeV (Cos #theta_{CM} -0.75-(-1.0))")
    PhiSc315AsymmCM8->Fit("SinFit", "Q");
    InitialSinAmp[7][0] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[7][0] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc315AsymmCM8->Fit("AsymmFit", "Q");
    Offset[7][0] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[7][0] = AsymmFit->GetParError(0);
    SinAmp[7][0] = AsymmFit->GetParameter(1);
    SinAmpErr[7][0] = AsymmFit->GetParError(1);
    CosAmp[7][0] = AsymmFit->GetParameter(2);
    CosAmpErr[7][0] = AsymmFit->GetParError(2);

    PhiSc345AsymmCM8 = PhiSc345NegHelCM8->GetAsymmetry(PhiSc345PosHelCM8);
    PhiSc345AsymmCM8->SetName("PhiSc345AsymmCM8");
    PhiSc345AsymmCM8->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 330-360MeV (Cos #theta_{CM} -0.75-(-1.0))")
    PhiSc345AsymmCM8->Fit("SinFit", "Q");
    InitialSinAmp[7][1] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[7][1] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc345AsymmCM8->Fit("AsymmFit", "Q");
    Offset[7][1] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[7][1] = AsymmFit->GetParError(0);
    SinAmp[7][1] = AsymmFit->GetParameter(1);
    SinAmpErr[7][1] = AsymmFit->GetParError(1);
    CosAmp[7][1] = AsymmFit->GetParameter(2);
    CosAmpErr[7][1] = AsymmFit->GetParError(2);

    PhiSc375AsymmCM8 = PhiSc375NegHelCM8->GetAsymmetry(PhiSc375PosHelCM8);
    PhiSc375AsymmCM8->SetName("PhiSc375AsymmCM8");
    PhiSc375AsymmCM8->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 360-390MeV (Cos #theta_{CM} -0.75-(-1.0))")
    PhiSc375AsymmCM8->Fit("SinFit", "Q");
    InitialSinAmp[7][2] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[7][2] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc375AsymmCM8->Fit("AsymmFit", "Q");
    Offset[7][2] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[7][2] = AsymmFit->GetParError(0);
    SinAmp[7][2] = AsymmFit->GetParameter(1);
    SinAmpErr[7][2] = AsymmFit->GetParError(1);
    CosAmp[7][2] = AsymmFit->GetParameter(2);
    CosAmpErr[7][2] = AsymmFit->GetParError(2);

    PhiSc405AsymmCM8 = PhiSc405NegHelCM8->GetAsymmetry(PhiSc405PosHelCM8);
    PhiSc405AsymmCM8->SetName("PhiSc405AsymmCM8");
    PhiSc405AsymmCM8->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 390-420MeV (Cos #theta_{CM} -0.75-(-1.0))")
    PhiSc405AsymmCM8->Fit("SinFit", "Q");
    InitialSinAmp[7][3] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[7][3] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc405AsymmCM8->Fit("AsymmFit", "Q");
    Offset[7][3] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[7][3] = AsymmFit->GetParError(0);
    SinAmp[7][3] = AsymmFit->GetParameter(1);
    SinAmpErr[7][3] = AsymmFit->GetParError(1);
    CosAmp[7][3] = AsymmFit->GetParameter(2);
    CosAmpErr[7][3] = AsymmFit->GetParError(2);

    PhiSc435AsymmCM8 = PhiSc435NegHelCM8->GetAsymmetry(PhiSc435PosHelCM8);
    PhiSc435AsymmCM8->SetName("PhiSc435AsymmCM8");
    PhiSc435AsymmCM8->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 420-450MeV (Cos #theta_{CM} -0.75-(-1.0))")
    PhiSc435AsymmCM8->Fit("SinFit", "Q");
    InitialSinAmp[7][4] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[7][4] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc435AsymmCM8->Fit("AsymmFit", "Q");
    Offset[7][4] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[7][4] = AsymmFit->GetParError(0);
    SinAmp[7][4] = AsymmFit->GetParameter(1);
    SinAmpErr[7][4] = AsymmFit->GetParError(1);
    CosAmp[7][4] = AsymmFit->GetParameter(2);
    CosAmpErr[7][4] = AsymmFit->GetParError(2);

    PhiSc465AsymmCM8 = PhiSc465NegHelCM8->GetAsymmetry(PhiSc465PosHelCM8);
    PhiSc465AsymmCM8->SetName("PhiSc465AsymmCM8");
    PhiSc465AsymmCM8->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 450-480MeV (Cos #theta_{CM} -0.75-(-1.0))")
    PhiSc465AsymmCM8->Fit("SinFit", "Q");
    InitialSinAmp[7][5] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[7][5] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc465AsymmCM8->Fit("AsymmFit", "Q");
    Offset[7][5] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[7][5] = AsymmFit->GetParError(0);
    SinAmp[7][5] = AsymmFit->GetParameter(1);
    SinAmpErr[7][5] = AsymmFit->GetParError(1);
    CosAmp[7][5] = AsymmFit->GetParameter(2);
    CosAmpErr[7][5] = AsymmFit->GetParError(2);

    PhiSc495AsymmCM8 = PhiSc495NegHelCM8->GetAsymmetry(PhiSc495PosHelCM8);
    PhiSc495AsymmCM8->SetName("PhiSc495AsymmCM8");
    PhiSc495AsymmCM8->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 480-510MeV (Cos #theta_{CM} -0.75-(-1.0))")
    PhiSc495AsymmCM8->Fit("SinFit", "Q");
    InitialSinAmp[7][6] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[7][6] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc495AsymmCM8->Fit("AsymmFit", "Q");
    Offset[7][6] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[7][6] = AsymmFit->GetParError(0);
    SinAmp[7][6] = AsymmFit->GetParameter(1);
    SinAmpErr[7][6] = AsymmFit->GetParError(1);
    CosAmp[7][6] = AsymmFit->GetParameter(2);
    CosAmpErr[7][6] = AsymmFit->GetParError(2);

    PhiSc525AsymmCM8 = PhiSc525NegHelCM8->GetAsymmetry(PhiSc525PosHelCM8);
    PhiSc525AsymmCM8->SetName("PhiSc525AsymmCM8");
    PhiSc525AsymmCM8->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 510-5400MeV (Cos #theta_{CM} -0.75-(-1.0))")
    PhiSc525AsymmCM8->Fit("SinFit", "Q");
    InitialSinAmp[7][7] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[7][7] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc525AsymmCM8->Fit("AsymmFit", "Q");
    Offset[7][7] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[7][7] = AsymmFit->GetParError(0);
    SinAmp[7][7] = AsymmFit->GetParameter(1);
    SinAmpErr[7][7] = AsymmFit->GetParError(1);
    CosAmp[7][7] = AsymmFit->GetParameter(2);
    CosAmpErr[7][7] = AsymmFit->GetParError(2);

    PhiSc555AsymmCM8 = PhiSc555NegHelCM8->GetAsymmetry(PhiSc555PosHelCM8);
    PhiSc555AsymmCM8->SetName("PhiSc555AsymmCM8");
    PhiSc555AsymmCM8->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 540-570MeV (Cos #theta_{CM} -0.75-(-1.0))")
    PhiSc555AsymmCM8->Fit("SinFit", "Q");
    InitialSinAmp[7][8] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[7][8] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc555AsymmCM8->Fit("AsymmFit", "Q");
    Offset[7][8] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[7][8] = AsymmFit->GetParError(0);
    SinAmp[7][8] = AsymmFit->GetParameter(1);
    SinAmpErr[7][8] = AsymmFit->GetParError(1);
    CosAmp[7][8] = AsymmFit->GetParameter(2);
    CosAmpErr[7][8] = AsymmFit->GetParError(2);

    PhiSc585AsymmCM8 = PhiSc585NegHelCM8->GetAsymmetry(PhiSc585PosHelCM8);
    PhiSc585AsymmCM8->SetName("PhiSc585AsymmCM8");
    PhiSc585AsymmCM8->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 570-600MeV (Cos #theta_{CM} -0.75-(-1.0))")
    PhiSc585AsymmCM8->Fit("SinFit", "Q");
    InitialSinAmp[7][9] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[7][9] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc585AsymmCM8->Fit("AsymmFit", "Q");
    Offset[7][9] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[7][9] = AsymmFit->GetParError(0);
    SinAmp[7][9] = AsymmFit->GetParameter(1);
    SinAmpErr[7][9] = AsymmFit->GetParError(1);
    CosAmp[7][9] = AsymmFit->GetParameter(2);
    CosAmpErr[7][9] = AsymmFit->GetParError(2);

    PhiSc615AsymmCM8 = PhiSc615NegHelCM8->GetAsymmetry(PhiSc615PosHelCM8);
    PhiSc615AsymmCM8->SetName("PhiSc615AsymmCM8");
    PhiSc615AsymmCM8->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 600-630MeV (Cos #theta_{CM} -0.75-(-1.0))")
    PhiSc615AsymmCM8->Fit("SinFit", "Q");
    InitialSinAmp[7][10] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[7][10] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc615AsymmCM8->Fit("AsymmFit", "Q");
    Offset[7][10] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[7][10] = AsymmFit->GetParError(0);
    SinAmp[7][10] = AsymmFit->GetParameter(1);
    SinAmpErr[7][10] = AsymmFit->GetParError(1);
    CosAmp[7][10] = AsymmFit->GetParameter(2);
    CosAmpErr[7][10] = AsymmFit->GetParError(2);

    PhiSc645AsymmCM8 = PhiSc645NegHelCM8->GetAsymmetry(PhiSc645PosHelCM8);
    PhiSc645AsymmCM8->SetName("PhiSc645AsymmCM8");
    PhiSc645AsymmCM8->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 630-660MeV (Cos #theta_{CM} -0.75-(-1.0))")
    PhiSc645AsymmCM8->Fit("SinFit", "Q");
    InitialSinAmp[7][11] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[7][11] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc645AsymmCM8->Fit("AsymmFit", "Q");
    Offset[7][11] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[7][11] = AsymmFit->GetParError(0);
    SinAmp[7][11] = AsymmFit->GetParameter(1);
    SinAmpErr[7][11] = AsymmFit->GetParError(1);
    CosAmp[7][11] = AsymmFit->GetParameter(2);
    CosAmpErr[7][11] = AsymmFit->GetParError(2);

    PhiSc675AsymmCM8 = PhiSc675NegHelCM8->GetAsymmetry(PhiSc675PosHelCM8);
    PhiSc675AsymmCM8->SetName("PhiSc675AsymmCM8");
    PhiSc675AsymmCM8->SetTitle("-ve/+ve Helicity Asymmetry for #phi_{Sc} for E_{#gamma} 660-690MeV (Cos #theta_{CM} -0.75-(-1.0))")
    PhiSc675AsymmCM8->Fit("SinFit", "Q");
    InitialSinAmp[7][12] = SinFunc->GetParameter(0);// Fit sine function to histogram
    InitialSinAmpErr[7][12] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
    PhiSc675AsymmCM8->Fit("AsymmFit", "Q");
    Offset[7][12] = AsymmFit->GetParameter(0); // Add values of the fit to an array
    OffsetErr[7][12] = AsymmFit->GetParError(0);
    SinAmp[7][12] = AsymmFit->GetParameter(1);
    SinAmpErr[7][12] = AsymmFit->GetParError(1);
    CosAmp[7][12] = AsymmFit->GetParameter(2);
    CosAmpErr[7][12] = AsymmFit->GetParError(2);

    // Define new file to store fit parameters
    TFile f1("AsymmFits_PTotal48.root", "RECREATE");

    PhiSc315AsymmCM1->Write();
    PhiSc345AsymmCM1->Write();
    PhiSc375AsymmCM1->Write();
    PhiSc405AsymmCM1->Write();
    PhiSc435AsymmCM1->Write();
    PhiSc465AsymmCM1->Write();
    PhiSc495AsymmCM1->Write();
    PhiSc525AsymmCM1->Write();
    PhiSc555AsymmCM1->Write();
    PhiSc585AsymmCM1->Write();
    PhiSc615AsymmCM1->Write();
    PhiSc645AsymmCM1->Write();
    PhiSc675AsymmCM1->Write();

    PhiSc315AsymmCM2->Write();
    PhiSc345AsymmCM2->Write();
    PhiSc375AsymmCM2->Write();
    PhiSc405AsymmCM2->Write();
    PhiSc435AsymmCM2->Write();
    PhiSc465AsymmCM2->Write();
    PhiSc495AsymmCM2->Write();
    PhiSc525AsymmCM2->Write();
    PhiSc555AsymmCM2->Write();
    PhiSc585AsymmCM2->Write();
    PhiSc615AsymmCM2->Write();
    PhiSc645AsymmCM2->Write();
    PhiSc675AsymmCM2->Write();

    PhiSc315AsymmCM3->Write();
    PhiSc345AsymmCM3->Write();
    PhiSc375AsymmCM3->Write();
    PhiSc405AsymmCM3->Write();
    PhiSc435AsymmCM3->Write();
    PhiSc465AsymmCM3->Write();
    PhiSc495AsymmCM3->Write();
    PhiSc525AsymmCM3->Write();
    PhiSc555AsymmCM3->Write();
    PhiSc585AsymmCM3->Write();
    PhiSc615AsymmCM3->Write();
    PhiSc645AsymmCM3->Write();
    PhiSc675AsymmCM3->Write();

    PhiSc315AsymmCM4->Write();
    PhiSc345AsymmCM4->Write();
    PhiSc375AsymmCM4->Write();
    PhiSc405AsymmCM4->Write();
    PhiSc435AsymmCM4->Write();
    PhiSc465AsymmCM4->Write();
    PhiSc495AsymmCM4->Write();
    PhiSc525AsymmCM4->Write();
    PhiSc555AsymmCM4->Write();
    PhiSc585AsymmCM4->Write();
    PhiSc615AsymmCM4->Write();
    PhiSc645AsymmCM4->Write();
    PhiSc675AsymmCM4->Write();

    PhiSc315AsymmCM5->Write();
    PhiSc345AsymmCM5->Write();
    PhiSc375AsymmCM5->Write();
    PhiSc405AsymmCM5->Write();
    PhiSc435AsymmCM5->Write();
    PhiSc465AsymmCM5->Write();
    PhiSc495AsymmCM5->Write();
    PhiSc525AsymmCM5->Write();
    PhiSc555AsymmCM5->Write();
    PhiSc585AsymmCM5->Write();
    PhiSc615AsymmCM5->Write();
    PhiSc645AsymmCM5->Write();
    PhiSc675AsymmCM5->Write();

    PhiSc315AsymmCM6->Write();
    PhiSc345AsymmCM6->Write();
    PhiSc375AsymmCM6->Write();
    PhiSc405AsymmCM6->Write();
    PhiSc435AsymmCM6->Write();
    PhiSc465AsymmCM6->Write();
    PhiSc495AsymmCM6->Write();
    PhiSc525AsymmCM6->Write();
    PhiSc555AsymmCM6->Write();
    PhiSc585AsymmCM6->Write();
    PhiSc615AsymmCM6->Write();
    PhiSc645AsymmCM6->Write();
    PhiSc675AsymmCM6->Write();

    PhiSc315AsymmCM7->Write();
    PhiSc345AsymmCM7->Write();
    PhiSc375AsymmCM7->Write();
    PhiSc405AsymmCM7->Write();
    PhiSc435AsymmCM7->Write();
    PhiSc465AsymmCM7->Write();
    PhiSc495AsymmCM7->Write();
    PhiSc525AsymmCM7->Write();
    PhiSc555AsymmCM7->Write();
    PhiSc585AsymmCM7->Write();
    PhiSc615AsymmCM7->Write();
    PhiSc645AsymmCM7->Write();
    PhiSc675AsymmCM7->Write();

    PhiSc315AsymmCM8->Write();
    PhiSc345AsymmCM8->Write();
    PhiSc375AsymmCM8->Write();
    PhiSc405AsymmCM8->Write();
    PhiSc435AsymmCM8->Write();
    PhiSc465AsymmCM8->Write();
    PhiSc495AsymmCM8->Write();
    PhiSc525AsymmCM8->Write();
    PhiSc555AsymmCM8->Write();
    PhiSc585AsymmCM8->Write();
    PhiSc615AsymmCM8->Write();
    PhiSc645AsymmCM8->Write();
    PhiSc675AsymmCM8->Write();

    //Define new tree to store parameters in
    TTree* tree = new TTree("Parameter_Values", "Tree_of_Values");

    // Define branches to store parameters, (Branch Name, Variable, Type of Variable)

    tree->Branch("InitialSinAmp315", &ISinAm315, "ISinAm315/D");
    tree->Branch("InitialSinAmpErr315", &ISinAmErr315, "ISinAmErr315/D");
    tree->Branch("CosAmp315", &CosAm315, "CosAm315/D");
    tree->Branch("CosAmpErr315", &CosAmErr315, "CosAmErr315/D");
    tree->Branch("SinAmp315", &SinAm315, "SinAm315/D");
    tree->Branch("SinAmpErr315", &SinAmErr315, "SinAmErr315/D");
    tree->Branch("Offset315", &Offs315, "Offs315/D");
    tree->Branch("OffsetErr315", &OffsErr315, "OffsErr315/D");
    tree->Branch("InitialSinAmp345", &ISinAm345, "ISinAm345/D");
    tree->Branch("InitialSinAmpErr345", &ISinAmErr345, "ISinAmErr345/D");
    tree->Branch("CosAmp345", &CosAm345, "CosAm345/D");
    tree->Branch("CosAmpErr345", &CosAmErr345, "CosAmErr345/D");
    tree->Branch("SinAmp345", &SinAm345, "SinAm345/D");
    tree->Branch("SinAmpErr345", &SinAmErr345, "SinAmErr345/D");
    tree->Branch("Offset345", &Offs345, "Offs345/D");
    tree->Branch("OffsetErr345", &OffsErr345, "OffsErr345/D");
    tree->Branch("InitialSinAmp375", &ISinAm375, "ISinAm375/D");
    tree->Branch("InitialSinAmpErr375", &ISinAmErr375, "ISinAmErr375/D");
    tree->Branch("CosAmp375", &CosAm375, "CosAm375/D");
    tree->Branch("CosAmpErr375", &CosAmErr375, "CosAmErr375/D");
    tree->Branch("SinAmp375", &SinAm375, "SinAm375/D");
    tree->Branch("SinAmpErr375", &SinAmErr375, "SinAmErr375/D");
    tree->Branch("Offset375", &Offs375, "Offs375/D");
    tree->Branch("OffsetErr375", &OffsErr375, "OffsErr375/D");
    tree->Branch("InitialSinAmp405", &ISinAm405, "ISinAm405/D");
    tree->Branch("InitialSinAmpErr405", &ISinAmErr405, "ISinAmErr405/D");
    tree->Branch("CosAmp405", &CosAm405, "CosAm405/D");
    tree->Branch("CosAmpErr405", &CosAmErr405, "CosAmErr405/D");
    tree->Branch("SinAmp405", &SinAm405, "SinAm405/D");
    tree->Branch("SinAmpErr405", &SinAmErr405, "SinAmErr405/D");
    tree->Branch("Offset405", &Offs405, "Offs405/D");
    tree->Branch("OffsetErr405", &OffsErr405, "OffsErr405/D");
    tree->Branch("InitialSinAmp435", &ISinAm435, "ISinAm435/D");
    tree->Branch("InitialSinAmpErr435", &ISinAmErr435, "ISinAmErr435/D");
    tree->Branch("CosAmp435", &CosAm435, "CosAm435/D");
    tree->Branch("CosAmpErr435", &CosAmErr435, "CosAmErr435/D");
    tree->Branch("SinAmp435", &SinAm435, "SinAm435/D");
    tree->Branch("SinAmpErr435", &SinAmErr435, "SinAmErr435/D");
    tree->Branch("Offset435", &Offs435, "Offs435/D");
    tree->Branch("OffsetErr435", &OffsErr435, "OffsErr435/D");
    tree->Branch("InitialSinAmp465", &ISinAm465, "ISinAm465/D");
    tree->Branch("InitialSinAmpErr465", &ISinAmErr465, "ISinAmErr465/D");
    tree->Branch("CosAmp465", &CosAm465, "CosAm465/D");
    tree->Branch("CosAmpErr465", &CosAmErr465, "CosAmErr465/D");
    tree->Branch("SinAmp465", &SinAm465, "SinAm465/D");
    tree->Branch("SinAmpErr465", &SinAmErr465, "SinAmErr465/D");
    tree->Branch("Offset465", &Offs465, "Offs465/D");
    tree->Branch("OffsetErr465", &OffsErr465, "OffsErr465/D");
    tree->Branch("InitialSinAmp495", &ISinAm495, "ISinAm495/D");
    tree->Branch("InitialSinAmpErr495", &ISinAmErr495, "ISinAmErr495/D");
    tree->Branch("CosAmp495", &CosAm495, "CosAm495/D");
    tree->Branch("CosAmpErr495", &CosAmErr495, "CosAmErr495/D");
    tree->Branch("SinAmp495", &SinAm495, "SinAm495/D");
    tree->Branch("SinAmpErr495", &SinAmErr495, "SinAmErr495/D");
    tree->Branch("Offset495", &Offs495, "Offs495/D");
    tree->Branch("OffsetErr495", &OffsErr495, "OffsErr495/D");
    tree->Branch("InitialSinAmp525", &ISinAm525, "ISinAm525/D");
    tree->Branch("InitialSinAmpErr525", &ISinAmErr525, "ISinAmErr525/D");
    tree->Branch("CosAmp525", &CosAm525, "CosAm525/D");
    tree->Branch("CosAmpErr525", &CosAmErr525, "CosAmErr525/D");
    tree->Branch("SinAmp525", &SinAm525, "SinAm525/D");
    tree->Branch("SinAmpErr525", &SinAmErr525, "SinAmErr525/D");
    tree->Branch("Offset525", &Offs525, "Offs525/D");
    tree->Branch("OffsetErr525", &OffsErr525, "OffsErr525/D");
    tree->Branch("InitialSinAmp555", &ISinAm555, "ISinAm555/D");
    tree->Branch("InitialSinAmpErr555", &ISinAmErr555, "ISinAmErr555/D");
    tree->Branch("CosAmp555", &CosAm555, "CosAm555/D");
    tree->Branch("CosAmpErr555", &CosAmErr555, "CosAmErr555/D");
    tree->Branch("SinAmp555", &SinAm555, "SinAm555/D");
    tree->Branch("SinAmpErr555", &SinAmErr555, "SinAmErr555/D");
    tree->Branch("Offset555", &Offs555, "Offs555/D");
    tree->Branch("OffsetErr555", &OffsErr555, "OffsErr555/D");
    tree->Branch("InitialSinAmp585", &ISinAm585, "ISinAm585/D");
    tree->Branch("InitialSinAmpErr585", &ISinAmErr585, "ISinAmErr585/D");
    tree->Branch("CosAmp585", &CosAm585, "CosAm585/D");
    tree->Branch("CosAmpErr585", &CosAmErr585, "CosAmErr585/D");
    tree->Branch("SinAmp585", &SinAm585, "SinAm585/D");
    tree->Branch("SinAmpErr585", &SinAmErr585, "SinAmErr585/D");
    tree->Branch("Offset585", &Offs585, "Offs585/D");
    tree->Branch("OffsetErr585", &OffsErr585, "OffsErr585/D");
    tree->Branch("InitialSinAmp615", &ISinAm615, "ISinAm615/D");
    tree->Branch("InitialSinAmpErr615", &ISinAmErr615, "ISinAmErr615/D");
    tree->Branch("CosAmp615", &CosAm615, "CosAm615/D");
    tree->Branch("CosAmpErr615", &CosAmErr615, "CosAmErr615/D");
    tree->Branch("SinAmp615", &SinAm615, "SinAm615/D");
    tree->Branch("SinAmpErr615", &SinAmErr615, "SinAmErr615/D");
    tree->Branch("Offset615", &Offs615, "Offs615/D");
    tree->Branch("OffsetErr615", &OffsErr615, "OffsErr615/D");
    tree->Branch("InitialSinAmp645", &ISinAm645, "ISinAm645/D");
    tree->Branch("InitialSinAmpErr645", &ISinAmErr645, "ISinAmErr645/D");
    tree->Branch("CosAmp645", &CosAm645, "CosAm645/D");
    tree->Branch("CosAmpErr645", &CosAmErr645, "CosAmErr645/D");
    tree->Branch("SinAmp645", &SinAm645, "SinAm645/D");
    tree->Branch("SinAmpErr645", &SinAmErr645, "SinAmErr645/D");
    tree->Branch("Offset645", &Offs645, "Offs645/D");
    tree->Branch("OffsetErr645", &OffsErr645, "OffsErr645/D");
    tree->Branch("InitialSinAmp675", &ISinAm675, "ISinAm675/D");
    tree->Branch("InitialSinAmpErr675", &ISinAmErr675, "ISinAmErr675/D");
    tree->Branch("CosAmp675", &CosAm675, "CosAm675/D");
    tree->Branch("CosAmpErr675", &CosAmErr675, "CosAmErr675/D");
    tree->Branch("SinAmp675", &SinAm675, "SinAm675/D");
    tree->Branch("SinAmpErr675", &SinAmErr675, "SinAmErr675/D");
    tree->Branch("Offset675", &Offs675, "Offs675/D");
    tree->Branch("OffsetErr675", &OffsErr675, "OffsErr675/D");

    for(Int_t m = 0; m < 8; m++){

        ISinAm315 = InitialSinAmp[m][0];
        ISinAmErr315 = InitialSinAmpErr[m][0];
        Offs315 = Offset[m][0];
        OffsErr315 = OffsetErr[m][0];
        SinAm315 = SinAmp[m][0];
        SinAmErr315 = SinAmpErr[m][0];
        CosAm315 = CosAmp[m][0];
        CosAmErr315 = CosAmpErr[m][0];
        ISinAm345 = InitialSinAmp[m][1];
        ISinAmErr345 = InitialSinAmpErr[m][1];
        Offs345 = Offset[m][1];
        OffsErr345 = OffsetErr[m][1];
        SinAm345 = SinAmp[m][1];
        SinAmErr345 = SinAmpErr[m][1];
        CosAm345 = CosAmp[m][1];
        CosAmErr345 = CosAmpErr[m][1];
        ISinAm375 = InitialSinAmp[m][2];
        ISinAmErr375 = InitialSinAmpErr[m][2];
        Offs375 = Offset[m][2];
        OffsErr375 = OffsetErr[m][2];
        SinAm375 = SinAmp[m][2];
        SinAmErr375 = SinAmpErr[m][2];
        CosAm375 = CosAmp[m][2];
        CosAmErr375 = CosAmpErr[m][2];
        ISinAm405 = InitialSinAmp[m][3];
        ISinAmErr405 = InitialSinAmpErr[m][3];
        Offs405 = Offset[m][3];
        OffsErr405 = OffsetErr[m][3];
        SinAm405 = SinAmp[m][3];
        SinAmErr405 = SinAmpErr[m][3];
        CosAm405 = CosAmp[m][3];
        CosAmErr405 = CosAmpErr[m][3];
        ISinAm435 = InitialSinAmp[m][4];
        ISinAmErr435 = InitialSinAmpErr[m][4];
        Offs435 = Offset[m][4];
        OffsErr435 = OffsetErr[m][4];
        SinAm435 = SinAmp[m][4];
        SinAmErr435 = SinAmpErr[m][4];
        CosAm435 = CosAmp[m][4];
        CosAmErr435 = CosAmpErr[m][4];
        ISinAm465 = InitialSinAmp[m][5];
        ISinAmErr465 = InitialSinAmpErr[m][5];
        Offs465 = Offset[m][5];
        OffsErr465 = OffsetErr[m][5];
        SinAm465 = SinAmp[m][5];
        SinAmErr465 = SinAmpErr[m][5];
        CosAm465 = CosAmp[m][5];
        CosAmErr465 = CosAmpErr[m][5];
        ISinAm495 = InitialSinAmp[m][6];
        ISinAmErr495 = InitialSinAmpErr[m][6];
        Offs495 = Offset[m][6];
        OffsErr495 = OffsetErr[m][6];
        SinAm495 = SinAmp[m][6];
        SinAmErr495 = SinAmpErr[m][6];
        CosAm495 = CosAmp[m][6];
        CosAmErr495 = CosAmpErr[m][6];
        ISinAm525 = InitialSinAmp[m][7];
        ISinAmErr525 = InitialSinAmpErr[m][7];
        Offs525 = Offset[m][7];
        OffsErr525 = OffsetErr[m][7];
        SinAm525 = SinAmp[m][7];
        SinAmErr525 = SinAmpErr[m][7];
        CosAm525 = CosAmp[m][7];
        CosAmErr525 = CosAmpErr[m][7];
        ISinAm555 = InitialSinAmp[m][8];
        ISinAmErr555 = InitialSinAmpErr[m][8];
        Offs555 = Offset[m][8];
        OffsErr555 = OffsetErr[m][8];
        SinAm555 = SinAmp[m][8];
        SinAmErr555 = SinAmpErr[m][8];
        CosAm555 = CosAmp[m][8];
        CosAmErr555 = CosAmpErr[m][8];
        ISinAm585 = InitialSinAmp[m][9];
        ISinAmErr585 = InitialSinAmpErr[m][9];
        Offs585 = Offset[m][9];
        OffsErr585 = OffsetErr[m][9];
        SinAm585 = SinAmp[m][9];
        SinAmErr585 = SinAmpErr[m][9];
        CosAm585 = CosAmp[m][9];
        CosAmErr585 = CosAmpErr[m][9];
        ISinAm615 = InitialSinAmp[m][10];
        ISinAmErr615 = InitialSinAmpErr[m][10];
        Offs615 = Offset[m][10];
        OffsErr615 = OffsetErr[m][10];
        SinAm615 = SinAmp[m][10];
        SinAmErr615 = SinAmpErr[m][10];
        CosAm615 = CosAmp[m][10];
        CosAmErr615 = CosAmpErr[m][10];
        ISinAm645 = InitialSinAmp[m][11];
        ISinAmErr645 = InitialSinAmpErr[m][11];
        Offs645 = Offset[m][11];
        OffsErr645 = OffsetErr[m][11];
        SinAm645 = SinAmp[m][11];
        SinAmErr645 = SinAmpErr[m][11];
        CosAm645 = CosAmp[m][11];
        CosAmErr645 = CosAmpErr[m][11];
        ISinAm675 = InitialSinAmp[m][12];
        ISinAmErr675 = InitialSinAmpErr[m][12];
        Offs675 = Offset[m][12];
        OffsErr675 = OffsetErr[m][12];
        SinAm675 = SinAmp[m][12];
        SinAmErr675 = SinAmpErr[m][12];
        CosAm675 = CosAmp[m][12];
        CosAmErr675 = CosAmpErr[m][12];

        tree->Fill();

    }

    f1.Write();
    f1.Close();

    TFile *f2= TFile::Open("/scratch/Mainz_Software/a2GoAT/CircPol_Aug16.root");

    // Calculate values of Cx for each angular and energy bin
    for (Int_t n = 0; n < 8; n++){

        Cx315[n] = SinAmp[n][0]/(Graph->Eval(315,0));
        CxErr315[n] = SinAmpErr[n][0]/((0.1)*(Graph->Eval(315,0)));
        Cx345[n] = SinAmp[n][1]/(Graph->Eval(345,0));
        CxErr345[n] = SinAmpErr[n][1]/((0.1)*(Graph->Eval(345,0)));
        Cx375[n] = SinAmp[n][2]/(Graph->Eval(375,0));
        CxErr375[n] = SinAmpErr[n][2]/((0.1)*(Graph->Eval(375,0)));
        Cx405[n] = SinAmp[n][3]/(Graph->Eval(405,0));
        CxErr405[n] = SinAmpErr[n][3]/((0.1)*(Graph->Eval(405,0)));
        Cx435[n] = SinAmp[n][4]/(Graph->Eval(435,0));
        CxErr435[n] = SinAmpErr[n][4]/((0.1)*(Graph->Eval(435,0)));
        Cx465[n] = SinAmp[n][5]/(Graph->Eval(465,0));
        CxErr465[n] = SinAmpErr[n][5]/((0.1)*(Graph->Eval(465,0)));
        Cx495[n] = SinAmp[n][6]/(Graph->Eval(495,0));
        CxErr495[n] = SinAmpErr[n][6]/((0.1)*(Graph->Eval(495,0)));
        Cx525[n] = SinAmp[n][7]/(Graph->Eval(525,0));
        CxErr525[n] = SinAmpErr[n][7]/((0.1)*(Graph->Eval(525,0)));
        Cx555[n] = SinAmp[n][8]/(Graph->Eval(555,0));
        CxErr555[n] = SinAmpErr[n][8]/((0.1)*(Graph->Eval(555,0)));
        Cx585[n] = SinAmp[n][9]/(Graph->Eval(585,0));
        CxErr585[n] = SinAmpErr[n][9]/((0.1)*(Graph->Eval(585,0)));
        Cx615[n] = SinAmp[n][10]/(Graph->Eval(615,0));
        CxErr615[n] = SinAmpErr[n][10]/((0.1)*(Graph->Eval(615,0)));
        Cx645[n] = SinAmp[n][11]/(Graph->Eval(645,0));
        CxErr645[n] = SinAmpErr[n][11]/((0.1)*(Graph->Eval(645,0)));
        Cx675[n] = SinAmp[n][12]/(Graph->Eval(675,0));
        CxErr675[n] = SinAmpErr[n][12]/((0.1)*(Graph->Eval(675,0)));

    }

    TFile f3("Cx_Plots_Multi_48.root", "RECREATE");

    Float_t xMin = -1;
    Float_t xMax = 1;
    Float_t yMin = -5;
    Float_t yMax = 5;
    Double_t x[8] = {0.875, 0.625, 0.375, 0.125, -0.125, -0.375, -0.625, -0.875}; // Need to adjust
    Double_t ex[8] = {0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125}; // Need to adjust

    TCanvas *canvas = new TCanvas("canvas","canvas", 1920, 1080);
    TPad *pad1 = new TPad("pad1","",0,0,1,1);
    pad1->Draw();
    pad1->cd();

    pad1->SetTickx(1);
    pad1->SetTicky(1);
    pad1->SetGridx(1);
    pad1->SetGridy(1);
    TH1F  *hr;
    hr = canvas->DrawFrame(xMin, -1,xMax, 1);
    hr->SetTitle("#C_{x} as fn of Cos#theta_{CM} (E_{#gamma} 300-330MeV)");

    gr1 = new TGraphErrors(8, x, Cx315, ex, CxErr315);
    gr1->SetMarkerColor(2);
    gr1->SetMarkerStyle(5);
    gr1->SetMarkerSize(2);
    gr1->SetTitle("C_{x} as fn of Cos#theta_{CM} (E_{#gamma} 300-330MeV)");
    gr1->GetXaxis()->SetTitle("Cos#theta_{CM}");
    gr1->GetYaxis()->SetTitle("C_{x}");
    gr1->Draw("ep");

    TCanvas *canvas1 = new TCanvas("canvas1","canvas1", 1920, 1080);
    TPad *pad2 = new TPad("pad2","",0,0,1,1);
    pad2->Draw();
    pad2->cd();

    pad2->SetTickx(1);
    pad2->SetTicky(1);
    pad2->SetGridx(1);
    pad2->SetGridy(1);
    TH1F  *hr1;
    hr1 = canvas1->DrawFrame(xMin, -1,xMax, 1);
    hr1->SetTitle("#C_{x} as fn of Cos#theta_{CM} (E_{#gamma} 330-360MeV)");

    gr2 = new TGraphErrors(8, x, Cx345, ex, CxErr345);
    gr2->SetMarkerColor(2);
    gr2->SetMarkerStyle(5);
    gr2->SetMarkerSize(2);
    gr2->SetTitle("C_{x} as fn of Cos#theta_{CM} (E_{#gamma} 330-360MeV)");
    gr2->GetXaxis()->SetTitle("Cos#theta_{CM}");
    gr2->GetYaxis()->SetTitle("C_{x}");
    gr2->Draw("ep");

    TCanvas *canvas2 = new TCanvas("canvas2","canvas2", 1920, 1080);
    TPad *pad3 = new TPad("pad3","",0,0,1,1);
    pad3->Draw();
    pad3->cd();

    pad3->SetTickx(1);
    pad3->SetTicky(1);
    pad3->SetGridx(1);
    pad3->SetGridy(1);
    TH1F  *hr2;
    hr2 = canvas2->DrawFrame(xMin, -1,xMax, 1);
    hr2->SetTitle("#C_{x} as fn of Cos#theta_{CM} (E_{#gamma} 360-390MeV)");

    gr3 = new TGraphErrors(8, x, Cx345, ex, CxErr345);
    gr3->SetMarkerColor(2);
    gr3->SetMarkerStyle(5);
    gr3->SetMarkerSize(2);
    gr3->SetTitle("C_{x} as fn of Cos#theta_{CM} (E_{#gamma} 360-390MeV)");
    gr3->GetXaxis()->SetTitle("Cos#theta_{CM}");
    gr3->GetYaxis()->SetTitle("C_{x}");
    gr3->Draw("ep");

    canvas->Write();
    gr1->SetName("Cx_315");
    Cx_315->Write();
    canvas1->Write();
    gr2->SetName("Cx_345");
    Cx_345->Write();
    canvas2->Write();
    gr3->SetName("Cx_375");
    Cx_375->Write();

    f3.Write();

}
