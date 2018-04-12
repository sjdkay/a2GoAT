#include "./includes.h"

// define a function with 4 parameters
Double_t fitf(Double_t *x,Double_t *par)
{
    Double_t fitval = 0;
    fitval = (1 + (par[0]*cos(x[0]))) + par[1];
    return fitval;
}

void CxAsymm_Denom() {

    double CosAmp[7][8]; // Format of arrays is Energy by i, Theta by j
    double CosAmpErr[7][8];
    double offset[7][8];
    double offsetErr[7][8];
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
    Func->SetParNames("CosAmp"); //Name the parameters
    Func->SetParLimits(0, -1000, 1000);
    TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/Amo/Physics_Total_108_6_4_18.root"); // Open the latest PTotal file to load histograms from
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
            FitHist[i][j]->Add(NegHelHist[i][j]);
            FitHist[i][j]->Fit("Fit", "Q");
            CosAmp[i][j] = Fit->GetParameter(0);
            CosAmpErr[i][j] = Fit->GetParError(0);
            offset[i][j] = Fit->GetParameter(1);
            offsetErr[i][j] = Fit->GetParError(1);
        }
    }

    // Define new file to store fit parameters
    TFile f1("SumFits_PTotal_108_v1.root", "RECREATE");

    for(Int_t i = 0; i < 7; i++){
        for(Int_t j = 0; j < 8; j++){
            FitHist[i][j]->Write();
        }
    }

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

        CosAm265 = CosAmp[0][m];
        CosAmErr265 = CosAmpErr[0][m];
        CosAm335 = CosAmp[1][m];
        CosAmErr335 = CosAmpErr[1][m];
        CosAm405 = CosAmp[2][m];
        CosAmErr405 = CosAmpErr[2][m];
        CosAm475 = CosAmp[3][m];
        CosAmErr475 = CosAmpErr[3][m];
        CosAm545 = CosAmp[4][m];
        CosAmErr545 = CosAmpErr[4][m];
        CosAm615 = CosAmp[5][m];
        CosAmErr615 = CosAmpErr[5][m];
        CosAm685 = CosAmp[6][m];
        CosAmErr685 = CosAmpErr[6][m];

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

    TFile f2("Cx_Denom_Plots_PTotal_108_v1.root", "RECREATE");

    Double_t x[8] = {0.875, 0.625, 0.375, 0.125, -0.125, -0.375, -0.625, -0.875}; // Need to adjust
    Double_t ex[8] = {0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125}; // Need to adjust

    for(Int_t k = 0; k < 7; k++){

        sprintf(P0name, "P0_%i", 265+(k*70));
        sprintf(P0title, "P0 E_{#gamma} %i #pm 35 MeV", 265+(k*70));
        P0Plots[k] = new TGraphErrors(8 , x, CosAmp[k], ex, CosAmpErr[k]);
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
