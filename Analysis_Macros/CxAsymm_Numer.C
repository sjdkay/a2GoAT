#include "./includes.h"

// define a function with 4 parameters
Double_t fitf(Double_t *x,Double_t *par)
{
    Double_t fitval = 0;
    fitval = (par[0]*sin(x[0])) + par[1];
    return fitval;
}

void CxAsymm_Numer() {

    double SinAmp[7][8]; // Format of arrays is Energy by i, Theta by j
    double SinAmpErr[7][8];
    double offset[7][8];
    double offsetErr[7][8];
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
    double Offset265;
    double OffsetErr265;
    double Offset335;
    double OffsetErr335;
    double Offset405;
    double OffsetErr405;
    double Offset475;
    double OffsetErr475;
    double Offset545;
    double OffsetErr545;
    double Offset615;
    double OffsetErr615;
    double Offset685;
    double OffsetErr685;
    Int_t nBins;
    TH1F* PosHelHist[7][8];
    TH1F* NegHelHist[7][8];
    TH1F* FitHist[7][8];
    char PosHelName[40];
    char NegHelName[40];
    char NewName[40];
    char Title[60];
    TGraphErrors* P0Plots[7];
    char P0name[60];
    char P0title[60];
    TGraphErrors* P1Plots[7];
    char P1name[60];
    char P1title[60];

    TF1 *Func = new TF1("Fit",  fitf, -3.0, 3.0, 2); //Give a name and range to the fitting funcion
    Func->SetParNames("SinAmp"); //Name the parameters
    Func->SetParLimits(0, -1000, 1000);
    TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/Amo/Physics_Total_118_01_05_18.root"); // Open the latest PTotal file to load histograms from
    nBins = Phi_Scattered_265MeV_NegHelCM1->GetNbinsX();

    for(Int_t i = 0; i < 7; i++){ // Energy
        for(Int_t j = 0; j < 8; j++){ // Theta
            sprintf(PosHelName, "Phi_Scattered_%iMeV_PosHelCM%i", 265+(i*70) , j+1);
            sprintf(NegHelName, "Phi_Scattered_%iMeV_NegHelCM%i", 265+(i*70) , j+1);
            sprintf(NewName, "FitHist_%iMeV_CM%i", 265+(i*70) , j+1);
            sprintf(Title, "%iMeV_CM%i", 265+(i*70) , j+1);
            PosHelHist[i][j] = (TH1F*)f->Get(PosHelName);
            NegHelHist[i][j] = (TH1F*)f->Get(NegHelName);
            FitHist[i][j] = new TH1F(NewName, Title, nBins, -4, 4);
            FitHist[i][j]->Add(PosHelHist[i][j]);
            FitHist[i][j]->Add(NegHelHist[i][j], -1);
            FitHist[i][j]->Fit("Fit", "Q");
            SinAmp[i][j] = Fit->GetParameter(0);
            SinAmpErr[i][j] = Fit->GetParError(0);
            offset[i][j] = Fit->GetParameter(1);
            offsetErr[i][j] = Fit->GetParError(1);
        }
    }

    // Define new file to store fit parameters
    TFile f1("DiffFits_PTotal_118_v1.root", "RECREATE");

    for(Int_t i = 0; i < 7; i++){
        for(Int_t j = 0; j < 8; j++){
            FitHist[i][j]->Write();
        }
    }

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
    tree->Branch("Offset265", &Offset265, "Offset265/D");
    tree->Branch("OffsetErr265", &OffsetErr265, "OffsetErr265/D");
    tree->Branch("Offset335", &Offset335, "Offset335/D");
    tree->Branch("OffsetErr335", &OffsetErr335, "OffsetErr335/D");
    tree->Branch("Offset405", &Offset405, "Offset405/D");
    tree->Branch("OffsetErr405", &OffsetErr405, "OffsetErr405/D");
    tree->Branch("Offset475", &Offset475, "Offset475/D");
    tree->Branch("OffsetErr475", &OffsetErr475, "OffsetErr475/D");
    tree->Branch("Offset545", &Offset545, "Offset545/D");
    tree->Branch("OffsetErr545", &OffsetErr545, "OffsetErr545/D");
    tree->Branch("Offset615", &Offset615, "Offset615/D");
    tree->Branch("OffsetErr615", &OffsetErr615, "OffsetErr615/D");
    tree->Branch("Offset685", &Offset685, "Offset685/D");
    tree->Branch("OffsetErr685", &OffsetErr685, "OffsetErr685/D");

    for(Int_t m = 0; m < 8; m++){

        SinAm265 = SinAmp[0][m];
        SinAmErr265 = SinAmpErr[0][m];
        SinAm335 = SinAmp[1][m];
        SinAmErr335 = SinAmpErr[1][m];
        SinAm405 = SinAmp[2][m];
        SinAmErr405 = SinAmpErr[2][m];
        SinAm475 = SinAmp[3][m];
        SinAmErr475 = SinAmpErr[3][m];
        SinAm545 = SinAmp[4][m];
        SinAmErr545 = SinAmpErr[4][m];
        SinAm615 = SinAmp[5][m];
        SinAmErr615 = SinAmpErr[5][m];
        SinAm685 = SinAmp[6][m];
        SinAmErr685 = SinAmpErr[6][m];

        Offset265 = offset[0][m];
        OffsetErr265 = offsetErr[0][m];
        Offset335 = offset[1][m];
        OffsetErr335 = offsetErr[1][m];
        Offset405 = offset[2][m];
        OffsetErr405 = offsetErr[2][m];
        Offset475 = offset[3][m];
        OffsetErr475 = offsetErr[3][m];
        Offset545 = offset[4][m];
        OffsetErr545 = offsetErr[4][m];
        Offset615 = offset[5][m];
        OffsetErr615 = offsetErr[5][m];
        Offset685 = offset[6][m];
        OffsetErr685 = offsetErr[6][m];

        tree->Fill();

    }

    f1.Write();
    f1.Close();

    TFile f2("Cx_Numer_Plots_PTotal_108_v1.root", "RECREATE");

    Double_t x[8] = {0.875, 0.625, 0.375, 0.125, -0.125, -0.375, -0.625, -0.875}; // Need to adjust
    Double_t ex[8] = {0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125}; // Need to adjust

    for(Int_t k = 0; k < 7; k++){

        sprintf(P0name, "P0_%i", 265+(k*70));
        sprintf(P0title, "P0 E_{#gamma} %i #pm 35 MeV", 265+(k*70));
        P0Plots[k] = new TGraphErrors(8 , x, SinAmp[k], ex, SinAmpErr[k]);
        P0Plots[k]->SetMarkerColor(4);
        P0Plots[k]->SetLineColor(4);
        P0Plots[k]->SetMarkerStyle(8);
        P0Plots[k]->SetMarkerSize(1);
        P0Plots[k]->GetXaxis()->SetTitle("Cos#theta_{CM}");
        P0Plots[k]->GetXaxis()->SetRangeUser(-1, 1);
        P0Plots[k]->GetYaxis()->SetTitle("P_{0}");
        P0Plots[k]->SetName(P0name);
        P0Plots[k]->SetTitle(P0title);

        sprintf(P1name, "P1_%i", 265+(k*70));
        sprintf(P1title, "P1 E_{#gamma} %i #pm 35 MeV", 265+(k*70));
        P1Plots[k] = new TGraphErrors(8 , x, offset[k], ex, offsetErr[k]);
        P1Plots[k]->SetMarkerColor(4);
        P1Plots[k]->SetLineColor(4);
        P1Plots[k]->SetMarkerStyle(8);
        P1Plots[k]->SetMarkerSize(1);
        P1Plots[k]->GetXaxis()->SetTitle("Cos#theta_{CM}");
        P1Plots[k]->GetXaxis()->SetRangeUser(-1, 1);
        P1Plots[k]->GetYaxis()->SetTitle("P_{1}");
        P1Plots[k]->SetName(P1name);
        P1Plots[k]->SetTitle(P1title);

        P0Plots[k]->Write();
        P1Plots[k]->Write();

    }

    TCanvas *canvas = new TCanvas("canvas","canvas", 1920, 1080);
    canvas->Divide(4,2);
    for(int i = 1 ; i < 8 ; i++){
        canvas->cd(i);
        P0Plots[i-1]->Draw("AEP");
    }

    TCanvas *canvas1 = new TCanvas("canvas1","canvas1", 1920, 1080);
    canvas1->Divide(4,2);
    for(int i = 1 ; i < 8 ; i++){
        canvas1->cd(i);
        P1Plots[i-1]->Draw("AEP");
    }

    canvas->Write();
    canvas1->Write();

    f2.Write();

}
