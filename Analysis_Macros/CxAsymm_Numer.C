#include "./includes.h"

// define a function with 4 parameters
Double_t fitf(Double_t *x,Double_t *par)
{
    Double_t fitval = 0;
    fitval = (par[0]*sin(x[0])) + par[1];
    return fitval;
}

void CxAsymm_Numer() {

    double SinAmp[8][8]; // Format of arrays is Theta Bin Defined by i, Energy bin defined by j
    double SinAmpErr[8][8];
    Int_t i;
    double SinAm265;
    double SinAmErr265;
    double SinAm335;
    double SinAmErr335;
    double SinAm405;
    double SinAmErr405;
    double SinAm475;
    double SinAmErr475;
    double SinAm545;
    double SinAmErr545;
    double SinAm615;
    double SinAmErr615;
    double SinAm685;
    double SinAmErr685;
    Int_t nBins;

    Double_t DiffAm265[8], DiffAm335[8], DiffAm405[8], DiffAm475[8], DiffAm545[8], DiffAm615[8], DiffAm685[8];
    Double_t DiffAmErr265[8], DiffAmErr335[8], DiffAmErr405[8], DiffAmErr475[8], DiffAmErr545[8], DiffAmErr615[8], DiffAmErr685[8];

    TF1 *Func = new TF1("Fit",  fitf, -3.0, 3.0, 2); //Give a name and range to the fitting funcion
    Func->SetParNames("SinAmp"); //Name the parameters
    Func->SetParLimits(0, -1000, 1000);
    TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/Amo/Physics_Total_108_6_4_18.root"); // Open the latest PTotal file to load histograms from
    nBins = Phi_Scattered_265MeV_NegHelCM1->GetNbinsX();

    ///////////////////////////////////////////
    //////////////////  CM1  //////////////////
    ///////////////////////////////////////////

    PhiSc265DiffCM1 = new TH1D("PhiSc265DiffCM1", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 230-300MeV (Cos #theta_{CM} 1-0.75)", nBins, -4, 4);
    PhiSc265DiffCM1->Add(Phi_Scattered_265MeV_PosHelCM1);
    PhiSc265DiffCM1->Add(Phi_Scattered_265MeV_NegHelCM1, -1);
    PhiSc265DiffCM1->Fit("Fit", "Q");
    SinAmp[0][0] = Fit->GetParameter(0);
    SinAmpErr[0][0] = Fit->GetParError(0);

    PhiSc335DiffCM1 = new TH1D("PhiSc335DiffCM1", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 300-370MeV (Cos #theta_{CM} 1-0.75)", nBins, -4, 4);
    PhiSc335DiffCM1->Add(Phi_Scattered_335MeV_PosHelCM1);
    PhiSc335DiffCM1->Add(Phi_Scattered_335MeV_NegHelCM1, -1);
    PhiSc335DiffCM1->Fit("Fit", "Q");
    SinAmp[0][1] = Fit->GetParameter(0);
    SinAmpErr[0][1] = Fit->GetParError(0);

    PhiSc405DiffCM1 = new TH1D("PhiSc405DiffCM1", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 370-440MeV (Cos #theta_{CM} 1-0.75)", nBins, -4, 4);
    PhiSc405DiffCM1->Add(Phi_Scattered_405MeV_NegHelCM1, Phi_Scattered_405MeV_PosHelCM1);
    PhiSc405DiffCM1->Fit("Fit", "Q");
    SinAmp[0][2] = Fit->GetParameter(0);
    SinAmpErr[0][2] = Fit->GetParError(0);

    PhiSc475DiffCM1 = new TH1D("PhiSc475DiffCM1", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 440-510MeV (Cos #theta_{CM} 1-0.75)", nBins, -4, 4);
    PhiSc475DiffCM1->Add(Phi_Scattered_475MeV_PosHelCM1);
    PhiSc475DiffCM1->Add(Phi_Scattered_475MeV_NegHelCM1, -1);
    PhiSc475DiffCM1->Fit("Fit", "Q");
    SinAmp[0][3] = Fit->GetParameter(0);
    SinAmpErr[0][3] = Fit->GetParError(0);

    PhiSc545DiffCM1 = new TH1D("PhiSc545DiffCM1", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 510-580MeV (Cos #theta_{CM} 1-0.75)", nBins, -4, 4);
    PhiSc545DiffCM1->Add(Phi_Scattered_545MeV_PosHelCM1);
    PhiSc545DiffCM1->Add(Phi_Scattered_545MeV_NegHelCM1, -1);
    PhiSc545DiffCM1->Fit("Fit", "Q");
    SinAmp[0][4] = Fit->GetParameter(0);
    SinAmpErr[0][4] = Fit->GetParError(0);

    PhiSc615DiffCM1 = new TH1D("PhiSc615DiffCM1", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 580-650MeV (Cos #theta_{CM} 1-0.75)", nBins, -4, 4);
    PhiSc615DiffCM1->Add(Phi_Scattered_615MeV_PosHelCM1);
    PhiSc615DiffCM1->Add(Phi_Scattered_615MeV_NegHelCM1, -1);
    PhiSc615DiffCM1->Fit("Fit", "Q");
    SinAmp[0][5] = Fit->GetParameter(0);
    SinAmpErr[0][5] = Fit->GetParError(0);

    PhiSc685DiffCM1 = new TH1D("PhiSc685DiffCM1", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 650-720MeV (Cos #theta_{CM} 1-0.75)", nBins, -4, 4);
    PhiSc685DiffCM1->Add(Phi_Scattered_685MeV_PosHelCM1);
    PhiSc685DiffCM1->Add(Phi_Scattered_685MeV_NegHelCM1, -1);
    PhiSc685DiffCM1->Fit("Fit", "Q");
    SinAmp[0][6] = Fit->GetParameter(0);
    SinAmpErr[0][6] = Fit->GetParError(0);

    ///////////////////////////////////////////
    //////////////////  CM2  //////////////////
    ///////////////////////////////////////////

    PhiSc265DiffCM2 = new TH1D("PhiSc265DiffCM2", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 230-300MeV (Cos #theta_{CM} 0.75-0.5)", nBins, -4, 4);
    PhiSc265DiffCM2->Add(Phi_Scattered_265MeV_PosHelCM2);
    PhiSc265DiffCM2->Add(Phi_Scattered_265MeV_NegHelCM2, -1);
    PhiSc265DiffCM2->Fit("Fit", "Q");
    SinAmp[1][0] = Fit->GetParameter(0);
    SinAmpErr[1][0] = Fit->GetParError(0);

    PhiSc335DiffCM2 = new TH1D("PhiSc335DiffCM2", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 300-370MeV (Cos #theta_{CM} 0.75-0.5)", nBins, -4, 4);
    PhiSc335DiffCM2->Add(Phi_Scattered_335MeV_PosHelCM2);
    PhiSc335DiffCM2->Add(Phi_Scattered_335MeV_NegHelCM2, -1);
    PhiSc335DiffCM2->Fit("Fit", "Q");
    SinAmp[1][1] = Fit->GetParameter(0);
    SinAmpErr[1][1] = Fit->GetParError(0);

    PhiSc405DiffCM2 = new TH1D("PhiSc405DiffCM2", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 370-440MeV (Cos #theta_{CM} 0.75-0.5)", nBins, -4, 4);
    PhiSc405DiffCM2->Add(Phi_Scattered_405MeV_NegHelCM2, Phi_Scattered_405MeV_PosHelCM2);
    PhiSc405DiffCM2->Fit("Fit", "Q");
    SinAmp[1][2] = Fit->GetParameter(0);
    SinAmpErr[1][2] = Fit->GetParError(0);

    PhiSc475DiffCM2 = new TH1D("PhiSc475DiffCM2", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 440-510MeV (Cos #theta_{CM} 0.75-0.5)", nBins, -4, 4);
    PhiSc475DiffCM2->Add(Phi_Scattered_475MeV_PosHelCM2);
    PhiSc475DiffCM2->Add(Phi_Scattered_475MeV_NegHelCM2, -1);
    PhiSc475DiffCM2->Fit("Fit", "Q");
    SinAmp[1][3] = Fit->GetParameter(0);
    SinAmpErr[1][3] = Fit->GetParError(0);

    PhiSc545DiffCM2 = new TH1D("PhiSc545DiffCM2", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 510-580MeV (Cos #theta_{CM} 0.75-0.5)", nBins, -4, 4);
    PhiSc545DiffCM2->Add(Phi_Scattered_545MeV_PosHelCM2);
    PhiSc545DiffCM2->Add(Phi_Scattered_545MeV_NegHelCM2, -1);
    PhiSc545DiffCM2->Fit("Fit", "Q");
    SinAmp[1][4] = Fit->GetParameter(0);
    SinAmpErr[1][4] = Fit->GetParError(0);

    PhiSc615DiffCM2 = new TH1D("PhiSc615DiffCM2", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 580-650MeV (Cos #theta_{CM} 0.75-0.5)", nBins, -4, 4);
    PhiSc615DiffCM2->Add(Phi_Scattered_615MeV_PosHelCM2);
    PhiSc615DiffCM2->Add(Phi_Scattered_615MeV_NegHelCM2, -1);
    PhiSc615DiffCM2->Fit("Fit", "Q");
    SinAmp[1][5] = Fit->GetParameter(0);
    SinAmpErr[1][5] = Fit->GetParError(0);

    PhiSc685DiffCM2 = new TH1D("PhiSc685DiffCM2", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 650-720MeV (Cos #theta_{CM} 0.75-0.5)", nBins, -4, 4);
    PhiSc685DiffCM2->Add(Phi_Scattered_685MeV_PosHelCM2);
    PhiSc685DiffCM2->Add(Phi_Scattered_685MeV_NegHelCM2, -1);
    PhiSc685DiffCM2->Fit("Fit", "Q");
    SinAmp[1][6] = Fit->GetParameter(0);
    SinAmpErr[1][6] = Fit->GetParError(0);

    ///////////////////////////////////////////
    //////////////////  CM3  //////////////////
    ///////////////////////////////////////////

    PhiSc265DiffCM3 = new TH1D("PhiSc265DiffCM3", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 230-300MeV (Cos #theta_{CM} 0.5-0.25)", nBins, -4, 4);
    PhiSc265DiffCM3->Add(Phi_Scattered_265MeV_PosHelCM3);
    PhiSc265DiffCM3->Add(Phi_Scattered_265MeV_NegHelCM3, -1);
    PhiSc265DiffCM3->Fit("Fit", "Q");
    SinAmp[2][0] = Fit->GetParameter(0);
    SinAmpErr[2][0] = Fit->GetParError(0);

    PhiSc335DiffCM3 = new TH1D("PhiSc335DiffCM3", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 300-370MeV (Cos #theta_{CM} 0.5-0.25)", nBins, -4, 4);
    PhiSc335DiffCM3->Add(Phi_Scattered_335MeV_PosHelCM3);
    PhiSc335DiffCM3->Add(Phi_Scattered_335MeV_NegHelCM3, -1);
    PhiSc335DiffCM3->Fit("Fit", "Q");
    SinAmp[2][1] = Fit->GetParameter(0);
    SinAmpErr[2][1] = Fit->GetParError(0);

    PhiSc405DiffCM3 = new TH1D("PhiSc405DiffCM3", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 370-440MeV (Cos #theta_{CM} 0.5-0.25)", nBins, -4, 4);
    PhiSc405DiffCM3->Add(Phi_Scattered_405MeV_NegHelCM3, Phi_Scattered_405MeV_PosHelCM3);
    PhiSc405DiffCM3->Fit("Fit", "Q");
    SinAmp[2][2] = Fit->GetParameter(0);
    SinAmpErr[2][2] = Fit->GetParError(0);

    PhiSc475DiffCM3 = new TH1D("PhiSc475DiffCM3", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 440-510MeV (Cos #theta_{CM} 0.5-0.25)", nBins, -4, 4);
    PhiSc475DiffCM3->Add(Phi_Scattered_475MeV_PosHelCM3);
    PhiSc475DiffCM3->Add(Phi_Scattered_475MeV_NegHelCM3, -1);
    PhiSc475DiffCM3->Fit("Fit", "Q");
    SinAmp[2][3] = Fit->GetParameter(0);
    SinAmpErr[2][3] = Fit->GetParError(0);

    PhiSc545DiffCM3 = new TH1D("PhiSc545DiffCM3", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 510-580MeV (Cos #theta_{CM} 0.5-0.25)", nBins, -4, 4);
    PhiSc545DiffCM3->Add(Phi_Scattered_545MeV_PosHelCM3);
    PhiSc545DiffCM3->Add(Phi_Scattered_545MeV_NegHelCM3, -1);
    PhiSc545DiffCM3->Fit("Fit", "Q");
    SinAmp[2][4] = Fit->GetParameter(0);
    SinAmpErr[2][4] = Fit->GetParError(0);

    PhiSc615DiffCM3 = new TH1D("PhiSc615DiffCM3", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 580-650MeV (Cos #theta_{CM} 0.5-0.25)", nBins, -4, 4);
    PhiSc615DiffCM3->Add(Phi_Scattered_615MeV_PosHelCM3);
    PhiSc615DiffCM3->Add(Phi_Scattered_615MeV_NegHelCM3, -1);
    PhiSc615DiffCM3->Fit("Fit", "Q");
    SinAmp[2][5] = Fit->GetParameter(0);
    SinAmpErr[2][5] = Fit->GetParError(0);

    PhiSc685DiffCM3 = new TH1D("PhiSc685DiffCM3", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 650-720MeV (Cos #theta_{CM} 0.5-0.25)", nBins, -4, 4);
    PhiSc685DiffCM3->Add(Phi_Scattered_685MeV_PosHelCM3);
    PhiSc685DiffCM3->Add(Phi_Scattered_685MeV_NegHelCM3, -1);
    PhiSc685DiffCM3->Fit("Fit", "Q");
    SinAmp[2][6] = Fit->GetParameter(0);
    SinAmpErr[2][6] = Fit->GetParError(0);

    ///////////////////////////////////////////
    //////////////////  CM4  //////////////////
    ///////////////////////////////////////////

    PhiSc265DiffCM4 = new TH1D("PhiSc265DiffCM4", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 230-300MeV (Cos #theta_{CM} 0.25-0.0)", nBins, -4, 4);
    PhiSc265DiffCM4->Add(Phi_Scattered_265MeV_PosHelCM4);
    PhiSc265DiffCM4->Add(Phi_Scattered_265MeV_NegHelCM4, -1);
    PhiSc265DiffCM4->Fit("Fit", "Q");
    SinAmp[3][0] = Fit->GetParameter(0);
    SinAmpErr[3][0] = Fit->GetParError(0);

    PhiSc335DiffCM4 = new TH1D("PhiSc335DiffCM4", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 300-370MeV (Cos #theta_{CM} 0.25-0.0)", nBins, -4, 4);
    PhiSc335DiffCM4->Add(Phi_Scattered_335MeV_PosHelCM4);
    PhiSc335DiffCM4->Add(Phi_Scattered_335MeV_NegHelCM4, -1);
    PhiSc335DiffCM4->Fit("Fit", "Q");
    SinAmp[3][1] = Fit->GetParameter(0);
    SinAmpErr[3][1] = Fit->GetParError(0);

    PhiSc405DiffCM4 = new TH1D("PhiSc405DiffCM4", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 370-440MeV (Cos #theta_{CM} 0.25-0.0)", nBins, -4, 4);
    PhiSc405DiffCM4->Add(Phi_Scattered_405MeV_NegHelCM4, Phi_Scattered_405MeV_PosHelCM4);
    PhiSc405DiffCM4->Fit("Fit", "Q");
    SinAmp[3][2] = Fit->GetParameter(0);
    SinAmpErr[3][2] = Fit->GetParError(0);

    PhiSc475DiffCM4 = new TH1D("PhiSc475DiffCM4", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 440-510MeV (Cos #theta_{CM} 0.25-0.0)", nBins, -4, 4);
    PhiSc475DiffCM4->Add(Phi_Scattered_475MeV_PosHelCM4);
    PhiSc475DiffCM4->Add(Phi_Scattered_475MeV_NegHelCM4, -1);
    PhiSc475DiffCM4->Fit("Fit", "Q");
    SinAmp[3][3] = Fit->GetParameter(0);
    SinAmpErr[3][3] = Fit->GetParError(0);

    PhiSc545DiffCM4 = new TH1D("PhiSc545DiffCM4", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 510-580MeV (Cos #theta_{CM} 0.25-0.0)", nBins, -4, 4);
    PhiSc545DiffCM4->Add(Phi_Scattered_545MeV_PosHelCM4);
    PhiSc545DiffCM4->Add(Phi_Scattered_545MeV_NegHelCM4, -1);
    PhiSc545DiffCM4->Fit("Fit", "Q");
    SinAmp[3][4] = Fit->GetParameter(0);
    SinAmpErr[3][4] = Fit->GetParError(0);

    PhiSc615DiffCM4 = new TH1D("PhiSc615DiffCM4", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 580-650MeV (Cos #theta_{CM} 0.25-0.0)", nBins, -4, 4);
    PhiSc615DiffCM4->Add(Phi_Scattered_615MeV_PosHelCM4);
    PhiSc615DiffCM4->Add(Phi_Scattered_615MeV_NegHelCM4, -1);
    PhiSc615DiffCM4->Fit("Fit", "Q");
    SinAmp[3][5] = Fit->GetParameter(0);
    SinAmpErr[3][5] = Fit->GetParError(0);

    PhiSc685DiffCM4 = new TH1D("PhiSc685DiffCM4", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 650-720MeV (Cos #theta_{CM} 0.25-0.0)", nBins, -4, 4);
    PhiSc685DiffCM4->Add(Phi_Scattered_685MeV_PosHelCM4);
    PhiSc685DiffCM4->Add(Phi_Scattered_685MeV_NegHelCM4, -1);
    PhiSc685DiffCM4->Fit("Fit", "Q");
    SinAmp[3][6] = Fit->GetParameter(0);
    SinAmpErr[3][6] = Fit->GetParError(0);

    ///////////////////////////////////////////
    //////////////////  CM5  //////////////////
    ///////////////////////////////////////////

    PhiSc265DiffCM5 = new TH1D("PhiSc265DiffCM5", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 230-300MeV (Cos #theta_{CM} 0.0-(-0.25))", nBins, -4, 4);
    PhiSc265DiffCM5->Add(Phi_Scattered_265MeV_PosHelCM5);
    PhiSc265DiffCM5->Add(Phi_Scattered_265MeV_NegHelCM5, -1);
    PhiSc265DiffCM5->Fit("Fit", "Q");
    SinAmp[4][0] = Fit->GetParameter(0);
    SinAmpErr[4][0] = Fit->GetParError(0);

    PhiSc335DiffCM5 = new TH1D("PhiSc335DiffCM5", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 300-370MeV (Cos #theta_{CM} 0.0-(-0.25))", nBins, -4, 4);
    PhiSc335DiffCM5->Add(Phi_Scattered_335MeV_PosHelCM5);
    PhiSc335DiffCM5->Add(Phi_Scattered_335MeV_NegHelCM5, -1);
    PhiSc335DiffCM5->Fit("Fit", "Q");
    SinAmp[4][1] = Fit->GetParameter(0);
    SinAmpErr[4][1] = Fit->GetParError(0);

    PhiSc405DiffCM5 = new TH1D("PhiSc405DiffCM5", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 370-440MeV (Cos #theta_{CM} 0.0-(-0.25))", nBins, -4, 4);
    PhiSc405DiffCM5->Add(Phi_Scattered_405MeV_NegHelCM5, Phi_Scattered_405MeV_PosHelCM5);
    PhiSc405DiffCM5->Fit("Fit", "Q");
    SinAmp[4][2] = Fit->GetParameter(0);
    SinAmpErr[4][2] = Fit->GetParError(0);

    PhiSc475DiffCM5 = new TH1D("PhiSc475DiffCM5", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 440-510MeV (Cos #theta_{CM} 0.0-(-0.25))", nBins, -4, 4);
    PhiSc475DiffCM5->Add(Phi_Scattered_475MeV_PosHelCM5);
    PhiSc475DiffCM5->Add(Phi_Scattered_475MeV_NegHelCM5, -1);
    PhiSc475DiffCM5->Fit("Fit", "Q");
    SinAmp[4][3] = Fit->GetParameter(0);
    SinAmpErr[4][3] = Fit->GetParError(0);

    PhiSc545DiffCM5 = new TH1D("PhiSc545DiffCM5", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 510-580MeV (Cos #theta_{CM} 0.0-(-0.25))", nBins, -4, 4);
    PhiSc545DiffCM5->Add(Phi_Scattered_545MeV_PosHelCM5);
    PhiSc545DiffCM5->Add(Phi_Scattered_545MeV_NegHelCM5, -1);
    PhiSc545DiffCM5->Fit("Fit", "Q");
    SinAmp[4][4] = Fit->GetParameter(0);
    SinAmpErr[4][4] = Fit->GetParError(0);

    PhiSc615DiffCM5 = new TH1D("PhiSc615DiffCM5", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 580-650MeV (Cos #theta_{CM} 0.0-(-0.25))", nBins, -4, 4);
    PhiSc615DiffCM5->Add(Phi_Scattered_615MeV_PosHelCM5);
    PhiSc615DiffCM5->Add(Phi_Scattered_615MeV_NegHelCM5, -1);
    PhiSc615DiffCM5->Fit("Fit", "Q");
    SinAmp[4][5] = Fit->GetParameter(0);
    SinAmpErr[4][5] = Fit->GetParError(0);

    PhiSc685DiffCM5 = new TH1D("PhiSc685DiffCM5", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 650-720MeV (Cos #theta_{CM} 0.0-(-0.25))", nBins, -4, 4);
    PhiSc685DiffCM5->Add(Phi_Scattered_685MeV_PosHelCM5);
    PhiSc685DiffCM5->Add(Phi_Scattered_685MeV_NegHelCM5, -1);
    PhiSc685DiffCM5->Fit("Fit", "Q");
    SinAmp[4][6] = Fit->GetParameter(0);
    SinAmpErr[4][6] = Fit->GetParError(0);

    ///////////////////////////////////////////
    //////////////////  CM6  //////////////////
    ///////////////////////////////////////////

    PhiSc265DiffCM6 = new TH1D("PhiSc265DiffCM6", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 230-300MeV (Cos #theta_{CM} -0.25-(-0.5))", nBins, -4, 4);
    PhiSc265DiffCM6->Add(Phi_Scattered_265MeV_PosHelCM6);
    PhiSc265DiffCM6->Add(Phi_Scattered_265MeV_NegHelCM6, -1);
    PhiSc265DiffCM6->Fit("Fit", "Q");
    SinAmp[5][0] = Fit->GetParameter(0);
    SinAmpErr[5][0] = Fit->GetParError(0);

    PhiSc335DiffCM6 = new TH1D("PhiSc335DiffCM6", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 300-370MeV (Cos #theta_{CM} -0.25-(-0.5))", nBins, -4, 4);
    PhiSc335DiffCM6->Add(Phi_Scattered_335MeV_PosHelCM6);
    PhiSc335DiffCM6->Add(Phi_Scattered_335MeV_NegHelCM6, -1);
    PhiSc335DiffCM6->Fit("Fit", "Q");
    SinAmp[5][1] = Fit->GetParameter(0);
    SinAmpErr[5][1] = Fit->GetParError(0);

    PhiSc405DiffCM6 = new TH1D("PhiSc405DiffCM6", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 370-440MeV (Cos #theta_{CM} -0.25-(-0.5))", nBins, -4, 4);
    PhiSc405DiffCM6->Add(Phi_Scattered_405MeV_NegHelCM6, Phi_Scattered_405MeV_PosHelCM6);
    PhiSc405DiffCM6->Fit("Fit", "Q");
    SinAmp[5][2] = Fit->GetParameter(0);
    SinAmpErr[5][2] = Fit->GetParError(0);

    PhiSc475DiffCM6 = new TH1D("PhiSc475DiffCM6", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 440-510MeV (Cos #theta_{CM} -0.25-(-0.5))", nBins, -4, 4);
    PhiSc475DiffCM6->Add(Phi_Scattered_475MeV_PosHelCM6);
    PhiSc475DiffCM6->Add(Phi_Scattered_475MeV_NegHelCM6, -1);
    PhiSc475DiffCM6->Fit("Fit", "Q");
    SinAmp[5][3] = Fit->GetParameter(0);
    SinAmpErr[5][3] = Fit->GetParError(0);

    PhiSc545DiffCM6 = new TH1D("PhiSc545DiffCM6", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 510-580MeV (Cos #theta_{CM} -0.25-(-0.5))", nBins, -4, 4);
    PhiSc545DiffCM6->Add(Phi_Scattered_545MeV_PosHelCM6);
    PhiSc545DiffCM6->Add(Phi_Scattered_545MeV_NegHelCM6, -1);
    PhiSc545DiffCM6->Fit("Fit", "Q");
    SinAmp[5][4] = Fit->GetParameter(0);
    SinAmpErr[5][4] = Fit->GetParError(0);

    PhiSc615DiffCM6 = new TH1D("PhiSc615DiffCM6", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 580-650MeV (Cos #theta_{CM} -0.25-(-0.5))", nBins, -4, 4);
    PhiSc615DiffCM6->Add(Phi_Scattered_615MeV_PosHelCM6);
    PhiSc615DiffCM6->Add(Phi_Scattered_615MeV_NegHelCM6, -1);
    PhiSc615DiffCM6->Fit("Fit", "Q");
    SinAmp[5][5] = Fit->GetParameter(0);
    SinAmpErr[5][5] = Fit->GetParError(0);

    PhiSc685DiffCM6 = new TH1D("PhiSc685DiffCM6", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 650-720MeV (Cos #theta_{CM} -0.25-(-0.5))", nBins, -4, 4);
    PhiSc685DiffCM6->Add(Phi_Scattered_685MeV_PosHelCM6);
    PhiSc685DiffCM6->Add(Phi_Scattered_685MeV_NegHelCM6, -1);
    PhiSc685DiffCM6->Fit("Fit", "Q");
    SinAmp[5][6] = Fit->GetParameter(0);
    SinAmpErr[5][6] = Fit->GetParError(0);

    ///////////////////////////////////////////
    //////////////////  CM7  //////////////////
    ///////////////////////////////////////////

    PhiSc265DiffCM7 = new TH1D("PhiSc265DiffCM7", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 230-300MeV (Cos #theta_{CM} -0.5-(-0.75))", nBins, -4, 4);
    PhiSc265DiffCM7->Add(Phi_Scattered_265MeV_PosHelCM7);
    PhiSc265DiffCM7->Add(Phi_Scattered_265MeV_NegHelCM7, -1);
    PhiSc265DiffCM7->Fit("Fit", "Q");
    SinAmp[6][0] = Fit->GetParameter(0);
    SinAmpErr[6][0] = Fit->GetParError(0);

    PhiSc335DiffCM7 = new TH1D("PhiSc335DiffCM7", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 300-370MeV (Cos #theta_{CM} -0.5-(-0.75))", nBins, -4, 4);
    PhiSc335DiffCM7->Add(Phi_Scattered_335MeV_PosHelCM7);
    PhiSc335DiffCM7->Add(Phi_Scattered_335MeV_NegHelCM7, -1);
    PhiSc335DiffCM7->Fit("Fit", "Q");
    SinAmp[6][1] = Fit->GetParameter(0);
    SinAmpErr[6][1] = Fit->GetParError(0);

    PhiSc405DiffCM7 = new TH1D("PhiSc405DiffCM7", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 370-440MeV (Cos #theta_{CM} -0.5-(-0.75))", nBins, -4, 4);
    PhiSc405DiffCM7->Add(Phi_Scattered_405MeV_NegHelCM7, Phi_Scattered_405MeV_PosHelCM7);
    PhiSc405DiffCM7->Fit("Fit", "Q");
    SinAmp[6][2] = Fit->GetParameter(0);
    SinAmpErr[6][2] = Fit->GetParError(0);

    PhiSc475DiffCM7 = new TH1D("PhiSc475DiffCM7", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 440-510MeV (Cos #theta_{CM} -0.5-(-0.75))", nBins, -4, 4);
    PhiSc475DiffCM7->Add(Phi_Scattered_475MeV_PosHelCM7);
    PhiSc475DiffCM7->Add(Phi_Scattered_475MeV_NegHelCM7, -1);
    PhiSc475DiffCM7->Fit("Fit", "Q");
    SinAmp[6][3] = Fit->GetParameter(0);
    SinAmpErr[6][3] = Fit->GetParError(0);

    PhiSc545DiffCM7 = new TH1D("PhiSc545DiffCM7", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 510-580MeV (Cos #theta_{CM} -0.5-(-0.75))", nBins, -4, 4);
    PhiSc545DiffCM7->Add(Phi_Scattered_545MeV_PosHelCM7);
    PhiSc545DiffCM7->Add(Phi_Scattered_545MeV_NegHelCM7, -1);
    PhiSc545DiffCM7->Fit("Fit", "Q");
    SinAmp[6][4] = Fit->GetParameter(0);
    SinAmpErr[6][4] = Fit->GetParError(0);

    PhiSc615DiffCM7 = new TH1D("PhiSc615DiffCM7", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 580-650MeV (Cos #theta_{CM} -0.5-(-0.75))", nBins, -4, 4);
    PhiSc615DiffCM7->Add(Phi_Scattered_615MeV_PosHelCM7);
    PhiSc615DiffCM7->Add(Phi_Scattered_615MeV_NegHelCM7, -1);
    PhiSc615DiffCM7->Fit("Fit", "Q");
    SinAmp[6][5] = Fit->GetParameter(0);
    SinAmpErr[6][5] = Fit->GetParError(0);

    PhiSc685DiffCM7 = new TH1D("PhiSc685DiffCM7", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 650-720MeV (Cos #theta_{CM} -0.5-(-0.75))", nBins, -4, 4);
    PhiSc685DiffCM7->Add(Phi_Scattered_685MeV_PosHelCM7);
    PhiSc685DiffCM7->Add(Phi_Scattered_685MeV_NegHelCM7, -1);
    PhiSc685DiffCM7->Fit("Fit", "Q");
    SinAmp[6][6] = Fit->GetParameter(0);
    SinAmpErr[6][6] = Fit->GetParError(0);

    ///////////////////////////////////////////
    //////////////////  CM8  //////////////////
    ///////////////////////////////////////////

    PhiSc265DiffCM8 = new TH1D("PhiSc265DiffCM8", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 230-300MeV (Cos #theta_{CM} -0.75-(-1.0))", nBins, -4, 4);
    PhiSc265DiffCM8->Add(Phi_Scattered_265MeV_PosHelCM8);
    PhiSc265DiffCM8->Add(Phi_Scattered_265MeV_NegHelCM8, -1);
    PhiSc265DiffCM8->Fit("Fit", "Q");
    SinAmp[7][0] = Fit->GetParameter(0);
    SinAmpErr[7][0] = Fit->GetParError(0);

    PhiSc335DiffCM8 = new TH1D("PhiSc335DiffCM8", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 300-370MeV (Cos #theta_{CM} -0.75-(-1.0))", nBins, -4, 4);
    PhiSc335DiffCM8->Add(Phi_Scattered_335MeV_PosHelCM8);
    PhiSc335DiffCM8->Add(Phi_Scattered_335MeV_NegHelCM8, -1);
    PhiSc335DiffCM8->Fit("Fit", "Q");
    SinAmp[7][1] = Fit->GetParameter(0);
    SinAmpErr[7][1] = Fit->GetParError(0);

    PhiSc405DiffCM8 = new TH1D("PhiSc405DiffCM8", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 370-440MeV (Cos #theta_{CM} -0.75-(-1.0))", nBins, -4, 4);
    PhiSc405DiffCM8->Add(Phi_Scattered_405MeV_NegHelCM8, Phi_Scattered_405MeV_PosHelCM8);
    PhiSc405DiffCM8->Fit("Fit", "Q");
    SinAmp[7][2] = Fit->GetParameter(0);
    SinAmpErr[7][2] = Fit->GetParError(0);

    PhiSc475DiffCM8 = new TH1D("PhiSc475DiffCM8", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 440-510MeV (Cos #theta_{CM} -0.75-(-1.0))", nBins, -4, 4);
    PhiSc475DiffCM8->Add(Phi_Scattered_475MeV_PosHelCM8);
    PhiSc475DiffCM8->Add(Phi_Scattered_475MeV_NegHelCM8, -1);
    PhiSc475DiffCM8->Fit("Fit", "Q");
    SinAmp[7][3] = Fit->GetParameter(0);
    SinAmpErr[7][3] = Fit->GetParError(0);

    PhiSc545DiffCM8 = new TH1D("PhiSc545DiffCM8", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 510-580MeV (Cos #theta_{CM} -0.75-(-1.0))", nBins, -4, 4);
    PhiSc545DiffCM8->Add(Phi_Scattered_545MeV_PosHelCM8);
    PhiSc545DiffCM8->Add(Phi_Scattered_545MeV_NegHelCM8, -1);
    PhiSc545DiffCM8->Fit("Fit", "Q");
    SinAmp[7][4] = Fit->GetParameter(0);
    SinAmpErr[7][4] = Fit->GetParError(0);

    PhiSc615DiffCM8 = new TH1D("PhiSc615DiffCM8", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 580-650MeV (Cos #theta_{CM} -0.75-(-1.0))", nBins, -4, 4);
    PhiSc615DiffCM8->Add(Phi_Scattered_615MeV_PosHelCM8);
    PhiSc615DiffCM8->Add(Phi_Scattered_615MeV_NegHelCM8, -1);
    PhiSc615DiffCM8->Fit("Fit", "Q");
    SinAmp[7][5] = Fit->GetParameter(0);
    SinAmpErr[7][5] = Fit->GetParError(0);

    PhiSc685DiffCM8 = new TH1D("PhiSc685DiffCM8", "-ve/+ve Helicity Diff for #phi_{Sc} for E_{#gamma} 650-720MeV (Cos #theta_{CM} -0.75-(-1.0))", nBins, -4, 4);
    PhiSc685DiffCM8->Add(Phi_Scattered_685MeV_PosHelCM8);
    PhiSc685DiffCM8->Add(Phi_Scattered_685MeV_NegHelCM8, -1);
    PhiSc685DiffCM8->Fit("Fit", "Q");
    SinAmp[7][6] = Fit->GetParameter(0);
    SinAmpErr[7][6] = Fit->GetParError(0);

    // Define new file to store fit parameters
    TFile f1("DiffFits_PTotal_108_v1.root", "RECREATE");

    PhiSc265DiffCM1->Write();
    PhiSc335DiffCM1->Write();
    PhiSc405DiffCM1->Write();
    PhiSc475DiffCM1->Write();
    PhiSc545DiffCM1->Write();
    PhiSc615DiffCM1->Write();
    PhiSc685DiffCM1->Write();

    PhiSc265DiffCM2->Write();
    PhiSc335DiffCM2->Write();
    PhiSc405DiffCM2->Write();
    PhiSc475DiffCM2->Write();
    PhiSc545DiffCM2->Write();
    PhiSc615DiffCM2->Write();
    PhiSc685DiffCM2->Write();

    PhiSc265DiffCM3->Write();
    PhiSc335DiffCM3->Write();
    PhiSc405DiffCM3->Write();
    PhiSc475DiffCM3->Write();
    PhiSc545DiffCM3->Write();
    PhiSc615DiffCM3->Write();
    PhiSc685DiffCM3->Write();

    PhiSc265DiffCM4->Write();
    PhiSc335DiffCM4->Write();
    PhiSc405DiffCM4->Write();
    PhiSc475DiffCM4->Write();
    PhiSc545DiffCM4->Write();
    PhiSc615DiffCM4->Write();
    PhiSc685DiffCM4->Write();

    PhiSc265DiffCM5->Write();
    PhiSc335DiffCM5->Write();
    PhiSc405DiffCM5->Write();
    PhiSc475DiffCM5->Write();
    PhiSc545DiffCM5->Write();
    PhiSc615DiffCM5->Write();
    PhiSc685DiffCM5->Write();

    PhiSc265DiffCM6->Write();
    PhiSc335DiffCM6->Write();
    PhiSc405DiffCM6->Write();
    PhiSc475DiffCM6->Write();
    PhiSc545DiffCM6->Write();
    PhiSc615DiffCM6->Write();
    PhiSc685DiffCM6->Write();

    PhiSc265DiffCM7->Write();
    PhiSc335DiffCM7->Write();
    PhiSc405DiffCM7->Write();
    PhiSc475DiffCM7->Write();
    PhiSc545DiffCM7->Write();
    PhiSc615DiffCM7->Write();
    PhiSc685DiffCM7->Write();

    PhiSc265DiffCM8->Write();
    PhiSc335DiffCM8->Write();
    PhiSc405DiffCM8->Write();
    PhiSc475DiffCM8->Write();
    PhiSc545DiffCM8->Write();
    PhiSc615DiffCM8->Write();
    PhiSc685DiffCM8->Write();

    //Define new tree to store parameters in
    TTree* tree = new TTree("Numerator_Values", "Tree_of_Values");

    // Define branches to store parameters, (Branch Name, Variable, Type of Variable)

    tree->Branch("SinAmp265", &SinAm265, "SinAm265/D");
    tree->Branch("SinAmpErr265", &SinAmErr265, "SinAmErr265/D");
    tree->Branch("SinAmp335", &SinAm335, "SinAm335/D");
    tree->Branch("SinAmpErr335", &SinAmErr335, "SinAmErr335/D");
    tree->Branch("SinAmp405", &SinAm405, "SinAm405/D");
    tree->Branch("SinAmpErr405", &SinAmErr405, "SinAmErr405/D");
    tree->Branch("SinAmp475", &SinAm475, "SinAm475/D");
    tree->Branch("SinAmpErr475", &SinAmErr475, "SinAmErr475/D");
    tree->Branch("SinAmp545", &SinAm545, "SinAm545/D");
    tree->Branch("SinAmpErr545", &SinAmErr545, "SinAmErr545/D");
    tree->Branch("SinAmp615", &SinAm615, "SinAm615/D");
    tree->Branch("SinAmpErr615", &SinAmErr615, "SinAmErr615/D");
    tree->Branch("SinAmp685", &SinAm685, "SinAm685/D");
    tree->Branch("SinAmpErr685", &SinAmErr685, "SinAmErr685/D");

    for(Int_t m = 0; m < 8; m++){

        SinAm265 = SinAmp[m][0];
        SinAmErr265 = SinAmpErr[m][1];
        SinAm335 = SinAmp[m][1];
        SinAmErr335 = SinAmpErr[m][1];
        SinAm405 = SinAmp[m][2];
        SinAmErr405 = SinAmpErr[m][2];
        SinAm475 = SinAmp[m][3];
        SinAmErr475 = SinAmpErr[m][3];
        SinAm545 = SinAmp[m][4];
        SinAmErr545 = SinAmpErr[m][4];
        SinAm615 = SinAmp[m][5];
        SinAmErr615 = SinAmpErr[m][5];
        SinAm685 = SinAmp[m][6];
        SinAmErr685 = SinAmpErr[m][6];

        DiffAm265[m] = SinAmp[m][0];
        DiffAmErr265[m] = SinAmpErr[m][0];
        DiffAm335[m] = SinAmp[m][1];
        DiffAmErr335[m] = SinAmpErr[m][1];
        DiffAm405[m] = SinAmp[m][2];
        DiffAmErr405[m] = SinAmpErr[m][2];
        DiffAm475[m] = SinAmp[m][3];
        DiffAmErr475[m] = SinAmpErr[m][3];
        DiffAm545[m] = SinAmp[m][4];
        DiffAmErr545[m] = SinAmpErr[m][4];
        DiffAm615[m] = SinAmp[m][5];
        DiffAmErr615[m] = SinAmpErr[m][5];
        DiffAm685[m] = SinAmp[m][6];
        DiffAmErr685[m] = SinAmpErr[m][6];

        tree->Fill();

    }

    f1.Write();
    f1.Close();

    TFile f3("Cx_Numer_Plots_PTotal_108_v1.root", "RECREATE");

    Float_t xMin = -1;
    Float_t xMax = 1;
    Float_t yMin = -10;
    Float_t yMax = 10;
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
    hr->SetTitle("#Delta^{#pm}(Cos#theta_{CM}) E_{#gamma} 230-300 MeV");

    gr = new TGraphErrors(8, x, DiffAm265, ex, DiffAmErr265);
    gr->SetMarkerColor(2);
    gr->SetMarkerStyle(5);
    gr->SetMarkerSize(2);
    gr->SetTitle("#Delta^{#pm}(Cos#theta_{CM}) E_{#gamma} 230-300 MeV");
    gr->GetXaxis()->SetTitle("Cos#theta_{CM}");
    gr->GetYaxis()->SetTitle("#Delta^{#pm} Amp");
    gr->SetName("DiffAm265");
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
    hr1->SetTitle("#Delta^{#pm}(Cos#theta_{CM}) E_{#gamma} 300-370 MeV");

    gr1 = new TGraphErrors(8, x, DiffAm335, ex, DiffAmErr335);
    gr1->SetMarkerColor(2);
    gr1->SetMarkerStyle(5);
    gr1->SetMarkerSize(2);
    gr1->SetTitle("#Delta^{#pm}(Cos#theta_{CM}) E_{#gamma} 300-370 MeV");
    gr1->GetXaxis()->SetTitle("Cos#theta_{CM}");
    gr1->GetYaxis()->SetTitle("#Delta^{#pm} Amp");
    gr1->SetName("DiffAm335");
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
    hr2->SetTitle("#Delta^{#pm}(Cos#theta_{CM}) E_{#gamma} 370-440 MeV");

    gr2 = new TGraphErrors(8, x, DiffAm405, ex, DiffAmErr405);
    gr2->SetMarkerColor(2);
    gr2->SetMarkerStyle(5);
    gr2->SetMarkerSize(2);
    gr2->SetTitle("#Delta^{#pm}(Cos#theta_{CM}) E_{#gamma} 370-440 MeV");
    gr2->GetXaxis()->SetTitle("Cos#theta_{CM}");
    gr2->GetYaxis()->SetTitle("#Delta^{#pm} Amp");
    gr2->SetName("DiffAm405");
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
    hr3->SetTitle("#Delta^{#pm}(Cos#theta_{CM}) E_{#gamma} 440-510 MeV");

    gr3 = new TGraphErrors(8, x, DiffAm475, ex, DiffAmErr475);
    gr3->SetMarkerColor(2);
    gr3->SetMarkerStyle(5);
    gr3->SetMarkerSize(2);
    gr3->SetTitle("#Delta^{#pm}(Cos#theta_{CM}) E_{#gamma} 440-510 MeV");
    gr3->GetXaxis()->SetTitle("Cos#theta_{CM}");
    gr3->GetYaxis()->SetTitle("#Delta^{#pm} Amp");
    gr3->SetName("DiffAm475");
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
    hr4->SetTitle("#Delta^{#pm}(Cos#theta_{CM}) E_{#gamma} 510-580 MeV");

    gr4 = new TGraphErrors(8, x, DiffAm545, ex, DiffAmErr545);
    gr4->SetMarkerColor(2);
    gr4->SetMarkerStyle(5);
    gr4->SetMarkerSize(2);
    gr4->SetTitle("#Delta^{#pm}(Cos#theta_{CM}) E_{#gamma} 510-580 MeV");
    gr4->GetXaxis()->SetTitle("Cos#theta_{CM}");
    gr4->GetYaxis()->SetTitle("#Delta^{#pm} Amp");
    gr4->SetName("DiffAm545");
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
    hr5->SetTitle("#Delta^{#pm}(Cos#theta_{CM}) E_{#gamma} 580-650 MeV");

    gr5 = new TGraphErrors(8, x, DiffAm615, ex, DiffAmErr615);
    gr5->SetMarkerColor(2);
    gr5->SetMarkerStyle(5);
    gr5->SetMarkerSize(2);
    gr5->SetTitle("#Delta^{#pm}(Cos#theta_{CM}) E_{#gamma} 580-650 MeV");
    gr5->GetXaxis()->SetTitle("Cos#theta_{CM}");
    gr5->GetYaxis()->SetTitle("#Delta^{#pm} Amp");
    gr5->SetName("DiffAm615");
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
    hr6->SetTitle("#Delta^{#pm}(Cos#theta_{CM}) E_{#gamma} 650-720 MeV");

    gr6 = new TGraphErrors(8, x, DiffAm685, ex, DiffAmErr685);
    gr6->SetMarkerColor(2);
    gr6->SetMarkerStyle(5);
    gr6->SetMarkerSize(2);
    gr6->SetTitle("#Delta^{#pm}(Cos#theta_{CM}) E_{#gamma} 650-720 MeV");
    gr6->GetXaxis()->SetTitle("Cos#theta_{CM}");
    gr6->GetYaxis()->SetTitle("#Delta^{#pm} Amp");
    gr6->SetName("DiffAm685");
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
