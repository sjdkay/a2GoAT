#include "./includes.h"

// define a function with 4 parameters
Double_t fitf(Double_t *x,Double_t *par)
{
    Double_t fitval = 0;
    fitval = (1 + (par[0]*cos(x[0]))) + par[1];
    return fitval;
}

void CxAsymm_Denom() {

    double CosAmp[8][8]; // Format of arrays is Theta Bin Defined by i, Energy bin defined by j
    double CosAmpErr[8][8];
    Int_t i;
    double CosAm265;
    double CosAmErr265;
    double CosAm335;
    double CosAmErr335;
    double CosAm405;
    double CosAmErr405;
    double CosAm475;
    double CosAmErr475;
    double CosAm545;
    double CosAmErr545;
    double CosAm615;
    double CosAmErr615;
    double CosAm685;
    double CosAmErr685;
    Int_t nBins;

    Double_t SumAm265[8], SumAm335[8], SumAm405[8], SumAm475[8], SumAm545[8], SumAm615[8], SumAm685[8];
    Double_t SumAmErr265[8], SumAmErr335[8], SumAmErr405[8], SumAmErr475[8], SumAmErr545[8], SumAmErr615[8], SumAmErr685[8];


    TF1 *Func = new TF1("Fit",  fitf, -3.0, 3.0, 2); //Give a name and range to the fitting funcion
    Func->SetParNames("CosAmp"); //Name the parameters
    Func->SetParLimits(0, -1000, 1000);
    Func->SetParLimits(1, -1000, 1000);
    TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/Amo/Physics_Total_108_6_4_18.root"); // Open the latest PTotal file to load histograms from
    nBins = Phi_Scattered_265MeV_NegHelCM1->GetNbinsX();

    ///////////////////////////////////////////
    //////////////////  CM1  //////////////////
    ///////////////////////////////////////////

    PhiSc265SumCM1 = new TH1D("PhiSc265SumCM1", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 230-300MeV (Cos #theta_{CM} 1-0.75)", nBins, -4, 4);
    PhiSc265SumCM1->Add(Phi_Scattered_265MeV_NegHelCM1, Phi_Scattered_265MeV_PosHelCM1);
    PhiSc265SumCM1->Fit("Fit", "Q");
    CosAmp[0][0] = Fit->GetParameter(1);
    CosAmpErr[0][0] = Fit->GetParError(1);

    PhiSc335SumCM1 = new TH1D("PhiSc335SumCM1", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 300-370MeV (Cos #theta_{CM} 1-0.75)", nBins, -4, 4);
    PhiSc335SumCM1->Add(Phi_Scattered_335MeV_NegHelCM1, Phi_Scattered_335MeV_PosHelCM1);
    PhiSc335SumCM1->Fit("Fit", "Q");
    CosAmp[0][1] = Fit->GetParameter(1);
    CosAmpErr[0][1] = Fit->GetParError(1);

    PhiSc405SumCM1 = new TH1D("PhiSc405SumCM1", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 370-440MeV (Cos #theta_{CM} 1-0.75)", nBins, -4, 4);
    PhiSc405SumCM1->Add(Phi_Scattered_405MeV_NegHelCM1, Phi_Scattered_405MeV_PosHelCM1);
    PhiSc405SumCM1->Fit("Fit", "Q");
    CosAmp[0][2] = Fit->GetParameter(1);
    CosAmpErr[0][2] = Fit->GetParError(1);

    PhiSc475SumCM1 = new TH1D("PhiSc475SumCM1", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 440-510MeV (Cos #theta_{CM} 1-0.75)", nBins, -4, 4);
    PhiSc475SumCM1->Add(Phi_Scattered_475MeV_NegHelCM1, Phi_Scattered_475MeV_PosHelCM1);
    PhiSc475SumCM1->Fit("Fit", "Q");
    CosAmp[0][3] = Fit->GetParameter(1);
    CosAmpErr[0][3] = Fit->GetParError(1);

    PhiSc545SumCM1 = new TH1D("PhiSc545SumCM1", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 510-580MeV (Cos #theta_{CM} 1-0.75)", nBins, -4, 4);
    PhiSc545SumCM1->Add(Phi_Scattered_545MeV_NegHelCM1, Phi_Scattered_545MeV_PosHelCM1);
    PhiSc545SumCM1->Fit("Fit", "Q");
    CosAmp[0][4] = Fit->GetParameter(1);
    CosAmpErr[0][4] = Fit->GetParError(1);

    PhiSc615SumCM1 = new TH1D("PhiSc615SumCM1", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 580-650MeV (Cos #theta_{CM} 1-0.75)", nBins, -4, 4);
    PhiSc615SumCM1->Add(Phi_Scattered_615MeV_NegHelCM1, Phi_Scattered_615MeV_PosHelCM1);
    PhiSc615SumCM1->Fit("Fit", "Q");
    CosAmp[0][5] = Fit->GetParameter(1);
    CosAmpErr[0][5] = Fit->GetParError(1);

    PhiSc685SumCM1 = new TH1D("PhiSc685SumCM1", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 650-720MeV (Cos #theta_{CM} 1-0.75)", nBins, -4, 4);
    PhiSc685SumCM1->Add(Phi_Scattered_685MeV_NegHelCM1, Phi_Scattered_685MeV_PosHelCM1);
    PhiSc685SumCM1->Fit("Fit", "Q");
    CosAmp[0][6] = Fit->GetParameter(1);
    CosAmpErr[0][6] = Fit->GetParError(1);

    ///////////////////////////////////////////
    //////////////////  CM2  //////////////////
    ///////////////////////////////////////////

    PhiSc265SumCM2 = new TH1D("PhiSc265SumCM2", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 230-300MeV (Cos #theta_{CM} 0.75-0.5)", nBins, -4, 4);
    PhiSc265SumCM2->Add(Phi_Scattered_265MeV_NegHelCM2, Phi_Scattered_265MeV_PosHelCM2);
    PhiSc265SumCM2->Fit("Fit", "Q");
    CosAmp[1][0] = Fit->GetParameter(1);
    CosAmpErr[1][0] = Fit->GetParError(1);

    PhiSc335SumCM2 = new TH1D("PhiSc335SumCM2", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 300-370MeV (Cos #theta_{CM} 0.75-0.5)", nBins, -4, 4);
    PhiSc335SumCM2->Add(Phi_Scattered_335MeV_NegHelCM2, Phi_Scattered_335MeV_PosHelCM2);
    PhiSc335SumCM2->Fit("Fit", "Q");
    CosAmp[1][1] = Fit->GetParameter(1);
    CosAmpErr[1][1] = Fit->GetParError(1);

    PhiSc405SumCM2 = new TH1D("PhiSc405SumCM2", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 370-440MeV (Cos #theta_{CM} 0.75-0.5)", nBins, -4, 4);
    PhiSc405SumCM2->Add(Phi_Scattered_405MeV_NegHelCM2, Phi_Scattered_405MeV_PosHelCM2);
    PhiSc405SumCM2->Fit("Fit", "Q");
    CosAmp[1][2] = Fit->GetParameter(1);
    CosAmpErr[1][2] = Fit->GetParError(1);

    PhiSc475SumCM2 = new TH1D("PhiSc475SumCM2", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 440-510MeV (Cos #theta_{CM} 0.75-0.5)", nBins, -4, 4);
    PhiSc475SumCM2->Add(Phi_Scattered_475MeV_NegHelCM2, Phi_Scattered_475MeV_PosHelCM2);
    PhiSc475SumCM2->Fit("Fit", "Q");
    CosAmp[1][3] = Fit->GetParameter(1);
    CosAmpErr[1][3] = Fit->GetParError(1);

    PhiSc545SumCM2 = new TH1D("PhiSc545SumCM2", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 510-580MeV (Cos #theta_{CM} 0.75-0.5)", nBins, -4, 4);
    PhiSc545SumCM2->Add(Phi_Scattered_545MeV_NegHelCM2, Phi_Scattered_545MeV_PosHelCM2);
    PhiSc545SumCM2->Fit("Fit", "Q");
    CosAmp[1][4] = Fit->GetParameter(1);
    CosAmpErr[1][4] = Fit->GetParError(1);

    PhiSc615SumCM2 = new TH1D("PhiSc615SumCM2", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 580-650MeV (Cos #theta_{CM} 0.75-0.5)", nBins, -4, 4);
    PhiSc615SumCM2->Add(Phi_Scattered_615MeV_NegHelCM2, Phi_Scattered_615MeV_PosHelCM2);
    PhiSc615SumCM2->Fit("Fit", "Q");
    CosAmp[1][5] = Fit->GetParameter(1);
    CosAmpErr[1][5] = Fit->GetParError(1);

    PhiSc685SumCM2 = new TH1D("PhiSc685SumCM2", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 650-720MeV (Cos #theta_{CM} 0.75-0.5)", nBins, -4, 4);
    PhiSc685SumCM2->Add(Phi_Scattered_685MeV_NegHelCM2, Phi_Scattered_685MeV_PosHelCM2);
    PhiSc685SumCM2->Fit("Fit", "Q");
    CosAmp[1][6] = Fit->GetParameter(1);
    CosAmpErr[1][6] = Fit->GetParError(1);

    ///////////////////////////////////////////
    //////////////////  CM3  //////////////////
    ///////////////////////////////////////////

    PhiSc265SumCM3 = new TH1D("PhiSc265SumCM3", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 230-300MeV (Cos #theta_{CM} 0.5-0.25)", nBins, -4, 4);
    PhiSc265SumCM3->Add(Phi_Scattered_265MeV_NegHelCM3, Phi_Scattered_265MeV_PosHelCM3);
    PhiSc265SumCM3->Fit("Fit", "Q");
    CosAmp[2][0] = Fit->GetParameter(1);
    CosAmpErr[2][0] = Fit->GetParError(1);

    PhiSc335SumCM3 = new TH1D("PhiSc335SumCM3", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 300-370MeV (Cos #theta_{CM} 0.5-0.25)", nBins, -4, 4);
    PhiSc335SumCM3->Add(Phi_Scattered_335MeV_NegHelCM3, Phi_Scattered_335MeV_PosHelCM3);
    PhiSc335SumCM3->Fit("Fit", "Q");
    CosAmp[2][1] = Fit->GetParameter(1);
    CosAmpErr[2][1] = Fit->GetParError(1);

    PhiSc405SumCM3 = new TH1D("PhiSc405SumCM3", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 370-440MeV (Cos #theta_{CM} 0.5-0.25)", nBins, -4, 4);
    PhiSc405SumCM3->Add(Phi_Scattered_405MeV_NegHelCM3, Phi_Scattered_405MeV_PosHelCM3);
    PhiSc405SumCM3->Fit("Fit", "Q");
    CosAmp[2][2] = Fit->GetParameter(1);
    CosAmpErr[2][2] = Fit->GetParError(1);

    PhiSc475SumCM3 = new TH1D("PhiSc475SumCM3", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 440-510MeV (Cos #theta_{CM} 0.5-0.25)", nBins, -4, 4);
    PhiSc475SumCM3->Add(Phi_Scattered_475MeV_NegHelCM3, Phi_Scattered_475MeV_PosHelCM3);
    PhiSc475SumCM3->Fit("Fit", "Q");
    CosAmp[2][3] = Fit->GetParameter(1);
    CosAmpErr[2][3] = Fit->GetParError(1);

    PhiSc545SumCM3 = new TH1D("PhiSc545SumCM3", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 510-580MeV (Cos #theta_{CM} 0.5-0.25)", nBins, -4, 4);
    PhiSc545SumCM3->Add(Phi_Scattered_545MeV_NegHelCM3, Phi_Scattered_545MeV_PosHelCM3);
    PhiSc545SumCM3->Fit("Fit", "Q");
    CosAmp[2][4] = Fit->GetParameter(1);
    CosAmpErr[2][4] = Fit->GetParError(1);

    PhiSc615SumCM3 = new TH1D("PhiSc615SumCM3", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 580-650MeV (Cos #theta_{CM} 0.5-0.25)", nBins, -4, 4);
    PhiSc615SumCM3->Add(Phi_Scattered_615MeV_NegHelCM3, Phi_Scattered_615MeV_PosHelCM3);
    PhiSc615SumCM3->Fit("Fit", "Q");
    CosAmp[2][5] = Fit->GetParameter(1);
    CosAmpErr[2][5] = Fit->GetParError(1);

    PhiSc685SumCM3 = new TH1D("PhiSc685SumCM3", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 650-720MeV (Cos #theta_{CM} 0.5-0.25)", nBins, -4, 4);
    PhiSc685SumCM3->Add(Phi_Scattered_685MeV_NegHelCM3, Phi_Scattered_685MeV_PosHelCM3);
    PhiSc685SumCM3->Fit("Fit", "Q");
    CosAmp[2][6] = Fit->GetParameter(1);
    CosAmpErr[2][6] = Fit->GetParError(1);

    ///////////////////////////////////////////
    //////////////////  CM4  //////////////////
    ///////////////////////////////////////////

    PhiSc265SumCM4 = new TH1D("PhiSc265SumCM4", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 230-300MeV (Cos #theta_{CM} 0.25-0.0)", nBins, -4, 4);
    PhiSc265SumCM4->Add(Phi_Scattered_265MeV_NegHelCM4, Phi_Scattered_265MeV_PosHelCM4);
    PhiSc265SumCM4->Fit("Fit", "Q");
    CosAmp[3][0] = Fit->GetParameter(1);
    CosAmpErr[3][0] = Fit->GetParError(1);

    PhiSc335SumCM4 = new TH1D("PhiSc335SumCM4", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 300-370MeV (Cos #theta_{CM} 0.25-0.0)", nBins, -4, 4);
    PhiSc335SumCM4->Add(Phi_Scattered_335MeV_NegHelCM4, Phi_Scattered_335MeV_PosHelCM4);
    PhiSc335SumCM4->Fit("Fit", "Q");
    CosAmp[3][1] = Fit->GetParameter(1);
    CosAmpErr[3][1] = Fit->GetParError(1);

    PhiSc405SumCM4 = new TH1D("PhiSc405SumCM4", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 370-440MeV (Cos #theta_{CM} 0.25-0.0)", nBins, -4, 4);
    PhiSc405SumCM4->Add(Phi_Scattered_405MeV_NegHelCM4, Phi_Scattered_405MeV_PosHelCM4);
    PhiSc405SumCM4->Fit("Fit", "Q");
    CosAmp[3][2] = Fit->GetParameter(1);
    CosAmpErr[3][2] = Fit->GetParError(1);

    PhiSc475SumCM4 = new TH1D("PhiSc475SumCM4", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 440-510MeV (Cos #theta_{CM} 0.25-0.0)", nBins, -4, 4);
    PhiSc475SumCM4->Add(Phi_Scattered_475MeV_NegHelCM4, Phi_Scattered_475MeV_PosHelCM4);
    PhiSc475SumCM4->Fit("Fit", "Q");
    CosAmp[3][3] = Fit->GetParameter(1);
    CosAmpErr[3][3] = Fit->GetParError(1);

    PhiSc545SumCM4 = new TH1D("PhiSc545SumCM4", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 510-580MeV (Cos #theta_{CM} 0.25-0.0)", nBins, -4, 4);
    PhiSc545SumCM4->Add(Phi_Scattered_545MeV_NegHelCM4, Phi_Scattered_545MeV_PosHelCM4);
    PhiSc545SumCM4->Fit("Fit", "Q");
    CosAmp[3][4] = Fit->GetParameter(1);
    CosAmpErr[3][4] = Fit->GetParError(1);

    PhiSc615SumCM4 = new TH1D("PhiSc615SumCM4", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 580-650MeV (Cos #theta_{CM} 0.25-0.0)", nBins, -4, 4);
    PhiSc615SumCM4->Add(Phi_Scattered_615MeV_NegHelCM4, Phi_Scattered_615MeV_PosHelCM4);
    PhiSc615SumCM4->Fit("Fit", "Q");
    CosAmp[3][5] = Fit->GetParameter(1);
    CosAmpErr[3][5] = Fit->GetParError(1);

    PhiSc685SumCM4 = new TH1D("PhiSc685SumCM4", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 650-720MeV (Cos #theta_{CM} 0.25-0.0)", nBins, -4, 4);
    PhiSc685SumCM4->Add(Phi_Scattered_685MeV_NegHelCM4, Phi_Scattered_685MeV_PosHelCM4);
    PhiSc685SumCM4->Fit("Fit", "Q");
    CosAmp[3][6] = Fit->GetParameter(1);
    CosAmpErr[3][6] = Fit->GetParError(1);

    ///////////////////////////////////////////
    //////////////////  CM5  //////////////////
    ///////////////////////////////////////////

    PhiSc265SumCM5 = new TH1D("PhiSc265SumCM5", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 230-300MeV (Cos #theta_{CM} 0.0-(-0.25))", nBins, -4, 4);
    PhiSc265SumCM5->Add(Phi_Scattered_265MeV_NegHelCM5, Phi_Scattered_265MeV_PosHelCM5);
    PhiSc265SumCM5->Fit("Fit", "Q");
    CosAmp[4][0] = Fit->GetParameter(1);
    CosAmpErr[4][0] = Fit->GetParError(1);

    PhiSc335SumCM5 = new TH1D("PhiSc335SumCM5", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 300-370MeV (Cos #theta_{CM} 0.0-(-0.25))", nBins, -4, 4);
    PhiSc335SumCM5->Add(Phi_Scattered_335MeV_NegHelCM5, Phi_Scattered_335MeV_PosHelCM5);
    PhiSc335SumCM5->Fit("Fit", "Q");
    CosAmp[4][1] = Fit->GetParameter(1);
    CosAmpErr[4][1] = Fit->GetParError(1);

    PhiSc405SumCM5 = new TH1D("PhiSc405SumCM5", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 370-440MeV (Cos #theta_{CM} 0.0-(-0.25))", nBins, -4, 4);
    PhiSc405SumCM5->Add(Phi_Scattered_405MeV_NegHelCM5, Phi_Scattered_405MeV_PosHelCM5);
    PhiSc405SumCM5->Fit("Fit", "Q");
    CosAmp[4][2] = Fit->GetParameter(1);
    CosAmpErr[4][2] = Fit->GetParError(1);

    PhiSc475SumCM5 = new TH1D("PhiSc475SumCM5", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 440-510MeV (Cos #theta_{CM} 0.0-(-0.25))", nBins, -4, 4);
    PhiSc475SumCM5->Add(Phi_Scattered_475MeV_NegHelCM5, Phi_Scattered_475MeV_PosHelCM5);
    PhiSc475SumCM5->Fit("Fit", "Q");
    CosAmp[4][3] = Fit->GetParameter(1);
    CosAmpErr[4][3] = Fit->GetParError(1);

    PhiSc545SumCM5 = new TH1D("PhiSc545SumCM5", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 510-580MeV (Cos #theta_{CM} 0.0-(-0.25))", nBins, -4, 4);
    PhiSc545SumCM5->Add(Phi_Scattered_545MeV_NegHelCM5, Phi_Scattered_545MeV_PosHelCM5);
    PhiSc545SumCM5->Fit("Fit", "Q");
    CosAmp[4][4] = Fit->GetParameter(1);
    CosAmpErr[4][4] = Fit->GetParError(1);

    PhiSc615SumCM5 = new TH1D("PhiSc615SumCM5", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 580-650MeV (Cos #theta_{CM} 0.0-(-0.25))", nBins, -4, 4);
    PhiSc615SumCM5->Add(Phi_Scattered_615MeV_NegHelCM5, Phi_Scattered_615MeV_PosHelCM5);
    PhiSc615SumCM5->Fit("Fit", "Q");
    CosAmp[4][5] = Fit->GetParameter(1);
    CosAmpErr[4][5] = Fit->GetParError(1);

    PhiSc685SumCM5 = new TH1D("PhiSc685SumCM5", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 650-720MeV (Cos #theta_{CM} 0.0-(-0.25))", nBins, -4, 4);
    PhiSc685SumCM5->Add(Phi_Scattered_685MeV_NegHelCM5, Phi_Scattered_685MeV_PosHelCM5);
    PhiSc685SumCM5->Fit("Fit", "Q");
    CosAmp[4][6] = Fit->GetParameter(1);
    CosAmpErr[4][6] = Fit->GetParError(1);

    ///////////////////////////////////////////
    //////////////////  CM6  //////////////////
    ///////////////////////////////////////////

    PhiSc265SumCM6 = new TH1D("PhiSc265SumCM6", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 230-300MeV (Cos #theta_{CM} -0.25-(-0.5))", nBins, -4, 4);
    PhiSc265SumCM6->Add(Phi_Scattered_265MeV_NegHelCM6, Phi_Scattered_265MeV_PosHelCM6);
    PhiSc265SumCM6->Fit("Fit", "Q");
    CosAmp[5][0] = Fit->GetParameter(1);
    CosAmpErr[5][0] = Fit->GetParError(1);

    PhiSc335SumCM6 = new TH1D("PhiSc335SumCM6", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 300-370MeV (Cos #theta_{CM} -0.25-(-0.5))", nBins, -4, 4);
    PhiSc335SumCM6->Add(Phi_Scattered_335MeV_NegHelCM6, Phi_Scattered_335MeV_PosHelCM6);
    PhiSc335SumCM6->Fit("Fit", "Q");
    CosAmp[5][1] = Fit->GetParameter(1);
    CosAmpErr[5][1] = Fit->GetParError(1);

    PhiSc405SumCM6 = new TH1D("PhiSc405SumCM6", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 370-440MeV (Cos #theta_{CM} -0.25-(-0.5))", nBins, -4, 4);
    PhiSc405SumCM6->Add(Phi_Scattered_405MeV_NegHelCM6, Phi_Scattered_405MeV_PosHelCM6);
    PhiSc405SumCM6->Fit("Fit", "Q");
    CosAmp[5][2] = Fit->GetParameter(1);
    CosAmpErr[5][2] = Fit->GetParError(1);

    PhiSc475SumCM6 = new TH1D("PhiSc475SumCM6", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 440-510MeV (Cos #theta_{CM} -0.25-(-0.5))", nBins, -4, 4);
    PhiSc475SumCM6->Add(Phi_Scattered_475MeV_NegHelCM6, Phi_Scattered_475MeV_PosHelCM6);
    PhiSc475SumCM6->Fit("Fit", "Q");
    CosAmp[5][3] = Fit->GetParameter(1);
    CosAmpErr[5][3] = Fit->GetParError(1);

    PhiSc545SumCM6 = new TH1D("PhiSc545SumCM6", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 510-580MeV (Cos #theta_{CM} -0.25-(-0.5))", nBins, -4, 4);
    PhiSc545SumCM6->Add(Phi_Scattered_545MeV_NegHelCM6, Phi_Scattered_545MeV_PosHelCM6);
    PhiSc545SumCM6->Fit("Fit", "Q");
    CosAmp[5][4] = Fit->GetParameter(1);
    CosAmpErr[5][4] = Fit->GetParError(1);

    PhiSc615SumCM6 = new TH1D("PhiSc615SumCM6", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 580-650MeV (Cos #theta_{CM} -0.25-(-0.5))", nBins, -4, 4);
    PhiSc615SumCM6->Add(Phi_Scattered_615MeV_NegHelCM6, Phi_Scattered_615MeV_PosHelCM6);
    PhiSc615SumCM6->Fit("Fit", "Q");
    CosAmp[5][5] = Fit->GetParameter(1);
    CosAmpErr[5][5] = Fit->GetParError(1);

    PhiSc685SumCM6 = new TH1D("PhiSc685SumCM6", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 650-720MeV (Cos #theta_{CM} -0.25-(-0.5))", nBins, -4, 4);
    PhiSc685SumCM6->Add(Phi_Scattered_685MeV_NegHelCM6, Phi_Scattered_685MeV_PosHelCM6);
    PhiSc685SumCM6->Fit("Fit", "Q");
    CosAmp[5][6] = Fit->GetParameter(1);
    CosAmpErr[5][6] = Fit->GetParError(1);

    ///////////////////////////////////////////
    //////////////////  CM7  //////////////////
    ///////////////////////////////////////////

    PhiSc265SumCM7 = new TH1D("PhiSc265SumCM7", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 230-300MeV (Cos #theta_{CM} -0.5-(-0.75))", nBins, -4, 4);
    PhiSc265SumCM7->Add(Phi_Scattered_265MeV_NegHelCM7, Phi_Scattered_265MeV_PosHelCM7);
    PhiSc265SumCM7->Fit("Fit", "Q");
    CosAmp[6][0] = Fit->GetParameter(1);
    CosAmpErr[6][0] = Fit->GetParError(1);

    PhiSc335SumCM7 = new TH1D("PhiSc335SumCM7", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 300-370MeV (Cos #theta_{CM} -0.5-(-0.75))", nBins, -4, 4);
    PhiSc335SumCM7->Add(Phi_Scattered_335MeV_NegHelCM7, Phi_Scattered_335MeV_PosHelCM7);
    PhiSc335SumCM7->Fit("Fit", "Q");
    CosAmp[6][1] = Fit->GetParameter(1);
    CosAmpErr[6][1] = Fit->GetParError(1);

    PhiSc405SumCM7 = new TH1D("PhiSc405SumCM7", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 370-440MeV (Cos #theta_{CM} -0.5-(-0.75))", nBins, -4, 4);
    PhiSc405SumCM7->Add(Phi_Scattered_405MeV_NegHelCM7, Phi_Scattered_405MeV_PosHelCM7);
    PhiSc405SumCM7->Fit("Fit", "Q");
    CosAmp[6][2] = Fit->GetParameter(1);
    CosAmpErr[6][2] = Fit->GetParError(1);

    PhiSc475SumCM7 = new TH1D("PhiSc475SumCM7", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 440-510MeV (Cos #theta_{CM} -0.5-(-0.75))", nBins, -4, 4);
    PhiSc475SumCM7->Add(Phi_Scattered_475MeV_NegHelCM7, Phi_Scattered_475MeV_PosHelCM7);
    PhiSc475SumCM7->Fit("Fit", "Q");
    CosAmp[6][3] = Fit->GetParameter(1);
    CosAmpErr[6][3] = Fit->GetParError(1);

    PhiSc545SumCM7 = new TH1D("PhiSc545SumCM7", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 510-580MeV (Cos #theta_{CM} -0.5-(-0.75))", nBins, -4, 4);
    PhiSc545SumCM7->Add(Phi_Scattered_545MeV_NegHelCM7, Phi_Scattered_545MeV_PosHelCM7);
    PhiSc545SumCM7->Fit("Fit", "Q");
    CosAmp[6][4] = Fit->GetParameter(1);
    CosAmpErr[6][4] = Fit->GetParError(1);

    PhiSc615SumCM7 = new TH1D("PhiSc615SumCM7", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 580-650MeV (Cos #theta_{CM} -0.5-(-0.75))", nBins, -4, 4);
    PhiSc615SumCM7->Add(Phi_Scattered_615MeV_NegHelCM7, Phi_Scattered_615MeV_PosHelCM7);
    PhiSc615SumCM7->Fit("Fit", "Q");
    CosAmp[6][5] = Fit->GetParameter(1);
    CosAmpErr[6][5] = Fit->GetParError(1);

    PhiSc685SumCM7 = new TH1D("PhiSc685SumCM7", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 650-720MeV (Cos #theta_{CM} -0.5-(-0.75))", nBins, -4, 4);
    PhiSc685SumCM7->Add(Phi_Scattered_685MeV_NegHelCM7, Phi_Scattered_685MeV_PosHelCM7);
    PhiSc685SumCM7->Fit("Fit", "Q");
    CosAmp[6][6] = Fit->GetParameter(1);
    CosAmpErr[6][6] = Fit->GetParError(1);

    ///////////////////////////////////////////
    //////////////////  CM8  //////////////////
    ///////////////////////////////////////////

    PhiSc265SumCM8 = new TH1D("PhiSc265SumCM8", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 230-300MeV (Cos #theta_{CM} -0.75-(-1.0))", nBins, -4, 4);
    PhiSc265SumCM8->Add(Phi_Scattered_265MeV_NegHelCM8, Phi_Scattered_265MeV_PosHelCM8);
    PhiSc265SumCM8->Fit("Fit", "Q");
    CosAmp[7][0] = Fit->GetParameter(1);
    CosAmpErr[7][0] = Fit->GetParError(1);

    PhiSc335SumCM8 = new TH1D("PhiSc335SumCM8", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 300-370MeV (Cos #theta_{CM} -0.75-(-1.0))", nBins, -4, 4);
    PhiSc335SumCM8->Add(Phi_Scattered_335MeV_NegHelCM8, Phi_Scattered_335MeV_PosHelCM8);
    PhiSc335SumCM8->Fit("Fit", "Q");
    CosAmp[7][1] = Fit->GetParameter(1);
    CosAmpErr[7][1] = Fit->GetParError(1);

    PhiSc405SumCM8 = new TH1D("PhiSc405SumCM8", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 370-440MeV (Cos #theta_{CM} -0.75-(-1.0))", nBins, -4, 4);
    PhiSc405SumCM8->Add(Phi_Scattered_405MeV_NegHelCM8, Phi_Scattered_405MeV_PosHelCM8);
    PhiSc405SumCM8->Fit("Fit", "Q");
    CosAmp[7][2] = Fit->GetParameter(1);
    CosAmpErr[7][2] = Fit->GetParError(1);

    PhiSc475SumCM8 = new TH1D("PhiSc475SumCM8", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 440-510MeV (Cos #theta_{CM} -0.75-(-1.0))", nBins, -4, 4);
    PhiSc475SumCM8->Add(Phi_Scattered_475MeV_NegHelCM8, Phi_Scattered_475MeV_PosHelCM8);
    PhiSc475SumCM8->Fit("Fit", "Q");
    CosAmp[7][3] = Fit->GetParameter(1);
    CosAmpErr[7][3] = Fit->GetParError(1);

    PhiSc545SumCM8 = new TH1D("PhiSc545SumCM8", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 510-580MeV (Cos #theta_{CM} -0.75-(-1.0))", nBins, -4, 4);
    PhiSc545SumCM8->Add(Phi_Scattered_545MeV_NegHelCM8, Phi_Scattered_545MeV_PosHelCM8);
    PhiSc545SumCM8->Fit("Fit", "Q");
    CosAmp[7][4] = Fit->GetParameter(1);
    CosAmpErr[7][4] = Fit->GetParError(1);

    PhiSc615SumCM8 = new TH1D("PhiSc615SumCM8", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 580-650MeV (Cos #theta_{CM} -0.75-(-1.0))", nBins, -4, 4);
    PhiSc615SumCM8->Add(Phi_Scattered_615MeV_NegHelCM8, Phi_Scattered_615MeV_PosHelCM8);
    PhiSc615SumCM8->Fit("Fit", "Q");
    CosAmp[7][5] = Fit->GetParameter(1);
    CosAmpErr[7][5] = Fit->GetParError(1);

    PhiSc685SumCM8 = new TH1D("PhiSc685SumCM8", "-ve/+ve Helicity Sum for #phi_{Sc} for E_{#gamma} 650-720MeV (Cos #theta_{CM} -0.75-(-1.0))", nBins, -4, 4);
    PhiSc685SumCM8->Add(Phi_Scattered_685MeV_NegHelCM8, Phi_Scattered_685MeV_PosHelCM8);
    PhiSc685SumCM8->Fit("Fit", "Q");
    CosAmp[7][6] = Fit->GetParameter(1);
    CosAmpErr[7][6] = Fit->GetParError(1);


    // Define new file to store fit parameters
    TFile f1("SumFits_PTotal_108_v1.root", "RECREATE");

    PhiSc265SumCM1->Write();
    PhiSc335SumCM1->Write();
    PhiSc405SumCM1->Write();
    PhiSc475SumCM1->Write();
    PhiSc545SumCM1->Write();
    PhiSc615SumCM1->Write();
    PhiSc685SumCM1->Write();

    PhiSc265SumCM2->Write();
    PhiSc335SumCM2->Write();
    PhiSc405SumCM2->Write();
    PhiSc475SumCM2->Write();
    PhiSc545SumCM2->Write();
    PhiSc615SumCM2->Write();
    PhiSc685SumCM2->Write();

    PhiSc265SumCM3->Write();
    PhiSc335SumCM3->Write();
    PhiSc405SumCM3->Write();
    PhiSc475SumCM3->Write();
    PhiSc545SumCM3->Write();
    PhiSc615SumCM3->Write();
    PhiSc685SumCM3->Write();

    PhiSc265SumCM4->Write();
    PhiSc335SumCM4->Write();
    PhiSc405SumCM4->Write();
    PhiSc475SumCM4->Write();
    PhiSc545SumCM4->Write();
    PhiSc615SumCM4->Write();
    PhiSc685SumCM4->Write();

    PhiSc265SumCM5->Write();
    PhiSc335SumCM5->Write();
    PhiSc405SumCM5->Write();
    PhiSc475SumCM5->Write();
    PhiSc545SumCM5->Write();
    PhiSc615SumCM5->Write();
    PhiSc685SumCM5->Write();

    PhiSc265SumCM6->Write();
    PhiSc335SumCM6->Write();
    PhiSc405SumCM6->Write();
    PhiSc475SumCM6->Write();
    PhiSc545SumCM6->Write();
    PhiSc615SumCM6->Write();
    PhiSc685SumCM6->Write();

    PhiSc265SumCM7->Write();
    PhiSc335SumCM7->Write();
    PhiSc405SumCM7->Write();
    PhiSc475SumCM7->Write();
    PhiSc545SumCM7->Write();
    PhiSc615SumCM7->Write();
    PhiSc685SumCM7->Write();

    PhiSc265SumCM8->Write();
    PhiSc335SumCM8->Write();
    PhiSc405SumCM8->Write();
    PhiSc475SumCM8->Write();
    PhiSc545SumCM8->Write();
    PhiSc615SumCM8->Write();
    PhiSc685SumCM8->Write();

    //Define new tree to store parameters in
    TTree* tree = new TTree("Denominator_Values", "Tree_of_Values");

    // Define branches to store parameters, (Branch Name, Variable, Type of Variable)

    tree->Branch("CosAmp265", &CosAm265, "CosAm265/D");
    tree->Branch("CosAmpErr265", &CosAmErr265, "CosAmErr265/D");
    tree->Branch("CosAmp335", &CosAm335, "CosAm335/D");
    tree->Branch("CosAmpErr335", &CosAmErr335, "CosAmErr335/D");
    tree->Branch("CosAmp405", &CosAm405, "CosAm405/D");
    tree->Branch("CosAmpErr405", &CosAmErr405, "CosAmErr405/D");
    tree->Branch("CosAmp475", &CosAm475, "CosAm475/D");
    tree->Branch("CosAmpErr475", &CosAmErr475, "CosAmErr475/D");
    tree->Branch("CosAmp545", &CosAm545, "CosAm545/D");
    tree->Branch("CosAmpErr545", &CosAmErr545, "CosAmErr545/D");
    tree->Branch("CosAmp615", &CosAm615, "CosAm615/D");
    tree->Branch("CosAmpErr615", &CosAmErr615, "CosAmErr615/D");
    tree->Branch("CosAmp685", &CosAm685, "CosAm685/D");
    tree->Branch("CosAmpErr685", &CosAmErr685, "CosAmErr685/D");

    for(Int_t m = 0; m < 8; m++){

        CosAm265 = CosAmp[m][0];
        CosAmErr265 = CosAmpErr[m][1];
        CosAm335 = CosAmp[m][1];
        CosAmErr335 = CosAmpErr[m][1];
        CosAm405 = CosAmp[m][2];
        CosAmErr405 = CosAmpErr[m][2];
        CosAm475 = CosAmp[m][3];
        CosAmErr475 = CosAmpErr[m][3];
        CosAm545 = CosAmp[m][4];
        CosAmErr545 = CosAmpErr[m][4];
        CosAm615 = CosAmp[m][5];
        CosAmErr615 = CosAmpErr[m][5];
        CosAm685 = CosAmp[m][6];
        CosAmErr685 = CosAmpErr[m][6];

        SumAm265[m] = CosAmp[m][0];
        SumAmErr265[m] = CosAmpErr[m][0];
        SumAm335[m] = CosAmp[m][1];
        SumAmErr335[m] = CosAmpErr[m][1];
        SumAm405[m] = CosAmp[m][2];
        SumAmErr405[m] = CosAmpErr[m][2];
        SumAm475[m] = CosAmp[m][3];
        SumAmErr475[m] = CosAmpErr[m][3];
        SumAm545[m] = CosAmp[m][4];
        SumAmErr545[m] = CosAmpErr[m][4];
        SumAm615[m] = CosAmp[m][5];
        SumAmErr615[m] = CosAmpErr[m][5];
        SumAm685[m] = CosAmp[m][6];
        SumAmErr685[m] = CosAmpErr[m][6];

        tree->Fill();

    }

    f1.Write();
    f1.Close();

    TFile f3("Cx_Denom_Plots_PTotal_108_v1.root", "RECREATE");

    Float_t xMin = -1;
    Float_t xMax = 1;
    Float_t yMin = -1000;
    Float_t yMax = 1000;
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
    hr = canvas->DrawFrame(xMin, yMin, xMax, yMax);
    hr->SetTitle("#Sigma^{#pm}(Cos#theta_{CM}) E_{#gamma} 230-300 MeV");

    gr = new TGraphErrors(8, x, SumAm265, ex, SumAmErr265);
    gr->SetMarkerColor(2);
    gr->SetMarkerStyle(5);
    gr->SetMarkerSize(2);
    gr->SetTitle("#Sigma^{#pm}(Cos#theta_{CM}) E_{#gamma} 230-300 MeV");
    gr->GetXaxis()->SetTitle("Cos#theta_{CM}");
    gr->GetYaxis()->SetTitle("#Sigma^{#pm} Amp");
    gr->SetName("SumAm265");
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
    hr1 = canvas1->DrawFrame(xMin, yMin, xMax, yMax);
    hr1->SetTitle("#Sigma^{#pm}(Cos#theta_{CM}) E_{#gamma} 300-370 MeV");

    gr1 = new TGraphErrors(8, x, SumAm335, ex, SumAmErr335);
    gr1->SetMarkerColor(2);
    gr1->SetMarkerStyle(5);
    gr1->SetMarkerSize(2);
    gr1->SetTitle("#Sigma^{#pm}(Cos#theta_{CM}) E_{#gamma} 300-370 MeV");
    gr1->GetXaxis()->SetTitle("Cos#theta_{CM}");
    gr1->GetYaxis()->SetTitle("#Sigma^{#pm} Amp");
    gr1->SetName("SumAm335");
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
    hr2 = canvas2->DrawFrame(xMin, yMin, xMax, yMax);
    hr2->SetTitle("#Sigma^{#pm}(Cos#theta_{CM}) E_{#gamma} 370-440 MeV");

    gr2 = new TGraphErrors(8, x, SumAm405, ex, SumAmErr405);
    gr2->SetMarkerColor(2);
    gr2->SetMarkerStyle(5);
    gr2->SetMarkerSize(2);
    gr2->SetTitle("#Sigma^{#pm}(Cos#theta_{CM}) E_{#gamma} 370-440 MeV");
    gr2->GetXaxis()->SetTitle("Cos#theta_{CM}");
    gr2->GetYaxis()->SetTitle("#Sigma^{#pm} Amp");
    gr2->SetName("SumAm405");
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
    hr3 = canvas3->DrawFrame(xMin, yMin, xMax, yMax);
    hr3->SetTitle("#Sigma^{#pm}(Cos#theta_{CM}) E_{#gamma} 440-510 MeV");

    gr3 = new TGraphErrors(8, x, SumAm475, ex, SumAmErr475);
    gr3->SetMarkerColor(2);
    gr3->SetMarkerStyle(5);
    gr3->SetMarkerSize(2);
    gr3->SetTitle("#Sigma^{#pm}(Cos#theta_{CM}) E_{#gamma} 440-510 MeV");
    gr3->GetXaxis()->SetTitle("Cos#theta_{CM}");
    gr3->GetYaxis()->SetTitle("#Sigma^{#pm} Amp");
    gr3->SetName("SumAm475");
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
    hr4 = canvas4->DrawFrame(xMin, yMin, xMax, yMax);
    hr4->SetTitle("#Sigma^{#pm}(Cos#theta_{CM}) (E_{#gamma} 510-580 MeV");

    gr4 = new TGraphErrors(8, x, SumAm545, ex, SumAmErr545);
    gr4->SetMarkerColor(2);
    gr4->SetMarkerStyle(5);
    gr4->SetMarkerSize(2);
    gr4->SetTitle("#Sigma^{#pm}(Cos#theta_{CM}) E_{#gamma} 510-580 MeV");
    gr4->GetXaxis()->SetTitle("Cos#theta_{CM}");
    gr4->GetYaxis()->SetTitle("#Sigma^{#pm} Amp");
    gr4->SetName("SumAm545");
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
    hr5 = canvas5->DrawFrame(xMin, yMin, xMax, yMax);
    hr5->SetTitle("#Sigma^{#pm}(Cos#theta_{CM}) E_{#gamma} 580-650 MeV");

    gr5 = new TGraphErrors(8, x, SumAm615, ex, SumAmErr615);
    gr5->SetMarkerColor(2);
    gr5->SetMarkerStyle(5);
    gr5->SetMarkerSize(2);
    gr5->SetTitle("#Sigma^{#pm}(Cos#theta_{CM}) E_{#gamma} 580-650 MeV");
    gr5->GetXaxis()->SetTitle("Cos#theta_{CM}");
    gr5->GetYaxis()->SetTitle("#Sigma^{#pm} Amp");
    gr5->SetName("SumAm615");
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
    hr6 = canvas6->DrawFrame(xMin, yMin, xMax, yMax);
    hr6->SetTitle("#Sigma^{#pm}(Cos#theta_{CM}) E_{#gamma} 650-720 MeV");

    gr6 = new TGraphErrors(8, x, SumAm685, ex, SumAmErr685);
    gr6->SetMarkerColor(2);
    gr6->SetMarkerStyle(5);
    gr6->SetMarkerSize(2);
    gr6->SetTitle("#Sigma^{#pm}(Cos#theta_{CM}) E_{#gamma} 650-720 MeV");
    gr6->GetXaxis()->SetTitle("Cos#theta_{CM}");
    gr6->GetYaxis()->SetTitle("#Sigma^{#pm} Amp");
    gr6->SetName("SumAm685");
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

    f3.Write();

}
