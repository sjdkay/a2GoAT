#include "./includes.h"

// define a function with 4 parameters
Double_t fitf(Double_t *x,Double_t *par)
{
    Double_t fitval = 0;
//    fitval =  ((par[0]*sin(x[0]))/(1 + (par[1]*cos(x[0]))));
    fitval =  par[0]*sin(x[0]);

    return fitval;
}

void CxAsymm_V2() {

    char PosHelHistName[60];
    char NegHelHistName[60];
    char AsymmHistName[60];
    char ISinAmBranchName[60];
    char ISinAmErBranchName[60];
    char SinAmBranchName[60];
    char SinAmErBranchName[60];
    char CosAmBranchName[60];
    char CosAmErBranchName[60];
    char name[60];
    char title[60];

    double InitialSinAmp[7][5];
    double InitialSinAmpErr[7][5];
    double SinAmp[7][5]; // Format of arrays is Theta Bin Defined by i, Energy bin defined by j
    double SinAmpErr[7][5];
    double CosAmp[7][5]; // Format of arrays is Theta Bin Defined by i, Energy bin defined by j
    double CosAmpErr[7][5];
    TGraphErrors* CxPlots[7];
    TGraphErrors* PnPlots[7];

    TH1F* AsymmHists[7][5];

    Double_t Cx[7][5];
    Double_t CxErr[7][5];
    Double_t Pn[7][5];
    Double_t PnErr[7][5];

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

    TF1 *AsymmFunc = new TF1("AsymmFit",  fitf, -3.0, 3.0, 1); //Give a name and range to the fitting funcion
//    AsymmFunc->SetParNames("SinAmp", "CosAmp"); //Name the parameters
    AsymmFunc->SetParNames("SinAmp"); //Name the parameters
    AsymmFunc->SetParameter(0, 0);
    TF1 *SinFunc = new TF1("SinFit", "[0]*sin(x*TMath::DegToRad())", -3, 3);
    SinFunc->SetParNames("InitialSinAmp");
    TFile *f = new TFile("/scratch/Mainz_Software/a2GoAT/Physics_Total_Amo122_Lin32_Combined.root"); // Open the latest PTotal file to load histograms from
    TF1 *Pn90CM = new TF1("Pn90CM", "1.64576-2.95484*(x/1000)+0.684577*(x/1000)**2-0.65*90**2/4/((x-560)**2+90**2/4)+(5.32305-35.3819*(x/1000)+70.145*(x/1000)**2-44.2899*(x/1000)**3)",300,700);

    for(Int_t i = 0; i < 6; i++){ // Energy
        for(Int_t j = 0; j < 3; j++){ // Theta
            sprintf(PosHelHistName, "PhiSc%iPosHelCM%i", 200+(i*100), j+1);
            sprintf(NegHelHistName, "PhiSc%iNegHelCM%i", 200+(i*100), j+1);
            sprintf(AsymmHistName, "CxAsymm%iCM%i", 200+(i*100), j+1);
            AsymmHists[i][j] = (TH1F*) (((TH1F*)f->Get(PosHelHistName))->GetAsymmetry(((TH1F*)f->Get(NegHelHistName)))));
            AsymmHists[i][j]->SetName(AsymmHistName);
            AsymmHists[i][j]->Fit("SinFit", "Q");
//            AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
//            AsymmFunc->SetParError(1, SinFunc->GetParError(0));
//            AsymmFunc->SetParLimits(0, -1, 1);
//            AsymmFunc->SetParLimits(1, -1, 1);
            AsymmHists[i][j]->Fit("AsymmFit", "QM");
            SinAmp[i][j] = AsymmFit->GetParameter(0);
            SinAmpErr[i][j] = AsymmFit->GetParError(0);
//            CosAmp[i][j] = AsymmFit->GetParameter(1);
//            CosAmpErr[i][j] = AsymmFit->GetParError(1);
        }
    }

    TFile f1("AsymmFits_PTotal_122_32_V2.root", "RECREATE");

    for(Int_t i = 0; i < 6; i++){ // Energy
        for(Int_t j = 0; j < 3; j++){ // Theta
            AsymmHists[i][j]->Write();
        }
    }

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

    //Fill branches (and hence tree) with corresponding parameters from above
//    for(Int_t A = 0; A < 7; A++){
//        sprintf(ISinAmBranchName, "ISinAm%i", 265+(A*70));
//        sprintf(ISinAmErBranchName, "ISinAmErr%i", 265+(A*70));
//        sprintf(SinAmBranchName, "SinAm%i", 265+(A*70));
//        sprintf(SinAmErBranchName, "SinAmErr%i", 265+(A*70));
//        sprintf(CosAmBranchName , "CosAm%i", 265+(A*70));
//        sprintf(CosAmErBranchName, "CosAmErr%i", 265+(A*70));
//        for(Int_t m = 0; m < 5; m++){
//            ((TBranch*)f1->Get(ISinAmBranchName)) = InitialSinAmp[A][m];
//            ((TBranch*)f1->Get(ISinAmErBranchName)) = InitialSinAmpErr[A][m];
//            ((TBranch*)f1->Get(SinAmErBranchName)) = SinAmp[A][m];
//            ((TBranch*)f1->Get(SinAmErBranchName)) = SinAmpErr[A][m];
//            ((TBranch*)f1->Get(CosAmErBranchName)) = CosAmp[A][m];
//            ((TBranch*)f1->Get(CosAmErBranchName)) = CosAmpErr[A][m];
//
//            tree->Fill();
//        }
//    }

    f1.Write();
    f1.Close();

    TFile *f2= TFile::Open("/scratch/Mainz_Software/a2GoAT/CircPol_Aug16.root");

    for (Int_t n = 0; n < 6; n++){
        for(Int_t k = 0; k < 3; k++){
        Double_t EPoint = 250 + (n*100);
        Cx[n][k] = SinAmp[n][k]/(0.1*(Graph->Eval(EPoint ,0)));
        CxErr[n][k] = SinAmpErr[n][k]/(0.1*(Graph->Eval(EPoint ,0)));
        }
    }

    TFile f3("Cx_Plots_122_32_V2.root", "RECREATE");

    Float_t xMin = -1;
    Float_t xMax = 1;
    double x[3] = {0.66, 0, -0.66};
    double ex[3] = {0.33, 0.33, 0.33};

    for(Int_t i = 0 ; i < 6 ; i++)
    {
        sprintf(name, "Cx_%i", 200+(i*100));
        sprintf(title, "C_{x}(Cos#theta_{CM}) E_{#gamma} %i #pm 50 MeV", 200+(i*100));
        CxPlots[i] = new TGraphErrors(3 , x, Cx[i], ex, CxErr[i]);
        CxPlots[i]->SetMarkerColor(4);
        CxPlots[i]->SetLineColor(4);
        CxPlots[i]->SetMarkerStyle(8);
        CxPlots[i]->SetMarkerSize(1);
        CxPlots[i]->GetXaxis()->SetTitle("Cos#theta_{CM}");
        CxPlots[i]->GetXaxis()->SetRangeUser(-1, 1);
//        CxPlots[i]->GetYaxis()->SetRangeUser(-1, 1);
        CxPlots[i]->GetYaxis()->SetTitle("C_{x}");
        CxPlots[i]->SetName(name);
        CxPlots[i]->SetTitle(title);

//        sprintf(name, "Pn_%i", 200+(i*100));
//        sprintf(title, "P_{n}(Cos#theta_{CM}) E_{#gamma} %i #pm 50 MeV", 200+(i*100));
//        PnPlots[i] = new TGraphErrors(3 , x, CosAmp[i], ex, CosAmpErr[i]);
//        PnPlots[i]->SetMarkerColor(4);
//        PnPlots[i]->SetLineColor(4);
//        PnPlots[i]->SetMarkerStyle(8);
//        PnPlots[i]->SetMarkerSize(1);
//        PnPlots[i]->GetXaxis()->SetTitle("Cos#theta_{CM}");
//        PnPlots[i]->GetXaxis()->SetRangeUser(-1, 1);
//        PnPlots[i]->GetYaxis()->SetRangeUser(-5, 5);
//        PnPlots[i]->GetYaxis()->SetTitle("P_{n}");
//        PnPlots[i]->SetName(name);
//        PnPlots[i]->SetTitle(title);

        CxPlots[i]->Write();
//        PnPlots[i]->Write();
    }

    TCanvas *canvas20 = new TCanvas("canvas20","canvas20", 1920, 1080);
    canvas20->Divide(4, 2);
    for(int i = 1 ; i < 7 ; i++){
        canvas20->cd(i);
        CxPlots[i-1]->Draw("AEP");
    }
    canvas20->Write();
//
//    TCanvas *canvas21 = new TCanvas("canvas21","canvas21", 1920, 1080);
//    canvas21->Divide(4, 2);
//    for(int i = 1 ; i < 7 ; i++){
//        canvas21->cd(i);
//        PnPlots[i-1]->Draw("AEP");
//    }
//
//    canvas21->Write();

    f3.Write();

}
