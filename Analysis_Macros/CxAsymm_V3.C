#include "./includes.h"

// define a function with 4 parameters
Double_t fitf(Double_t *x,Double_t *par)
{
    Double_t fitval = 0;
    fitval =  ((par[0]*sin(x[0]))/(1 + (par[1]*cos(x[0]))));
    return fitval;
}

void CxAsymm_V3() {

    char PosHelHistName[60];
    char NegHelHistName[60];
    char AsymmHistName[60];
    char NeutronEThetaScName[60];
    char name[60];
    char title[60];

    TH1F* AsymmHists[4][4];
    TGraphErrors* CxPlots[4];

    double SinAmps[4][4]; // To store and get Cx amps
    double CosAmps[2][4]; // To store and get P amps
    double SinAmpErrs[4][4]; // To store and get Cx amps
    double CosAmpErrs[2][4]; // To store and get P amps
    double MeanX[4];
    double MeanY[4];

    Double_t Cx[4][4]; // Fixed P Chi2 fit, Fixed P LL fit, Variable p Chi2 fit, Variable P LL fit
    Double_t CxErr[4][4];

    Double_t PVal;
    Double_t ECentre;
    Double_t AmpVal;
    Double_t AEff[4];
    Double_t CorrFac[4];
    Double_t AEffCorr[4];

    TF1 *AsymmFunc = new TF1("AsymmFit",  fitf, -3.0, 3.0, 2); //Give a name and range to the fitting funcion
    AsymmFunc->SetParNames("SinAmp", "CosAmp"); //Name the parameters
    AsymmFunc->SetParameter(0, 0);
    TF1 *SinFunc = new TF1("SinFit", "[0]*sin(x)", -3, 3);
    SinFunc->SetParNames("InitialSinAmp");
    TFile *f = new TFile("/scratch/Mainz_Software/a2GoAT/Physics_Total_Amo137_Lin37_Combined.root"); // Open the latest PTotal file to load histograms from
    TFile *fAy = new TFile ("/scratch/Mainz_Software/a2GoAT/npAy.root");
    TF1 *Pn90CM = new TF1("Pn90CM", "1.64576-2.95484*(x/1000)+0.684577*(x/1000)**2-0.65*90**2/4/((x-560)**2+90**2/4)+(5.32305-35.3819*(x/1000)+70.145*(x/1000)**2-44.2899*(x/1000)**3)",300,700);

    for(Int_t i = 0; i < 2; i++){ // Fit version
        for(Int_t j = 0; j < 4; j++){ // Energy

            sprintf(PosHelHistName, "PhiSc%iPosHelCM2", 400+(j*100));
            sprintf(NegHelHistName, "PhiSc%iNegHelCM2", 400+(j*100));
            sprintf(AsymmHistName, "CxAsymm%iCM2Fit%i", 400+(j*100), i+1);
            sprintf(NeutronEThetaScName, "NeutronEThetaSc%iCM2", 400+(j*100));
            AsymmHists[i][j] = (TH1F*) (((TH1F*)f->Get(PosHelHistName))->GetAsymmetry(((TH1F*)f->Get(NegHelHistName)))));
            AsymmHists[i][j]->SetName(AsymmHistName);

            MeanX[j] = ((TH2D*)f->Get(NeutronEThetaScName))->GetMean(1);
            MeanY[j] = ((TH2D*)f->Get(NeutronEThetaScName))->GetMean(2);

            ECentre = 450+(j*100); // Get centre of energy bin
            PVal = (Pn90CM->Eval(ECentre, 0)); // Get PValue for energy bin to fix parameter with
            AEff[j] = ((TH2F*)fAy->Get("nppnAy"))->Interpolate(MeanX[j], MeanY[j]);
            CorrFac[j] = (1+exp(1.81572-(0.0139530*MeanX[j])));
            AEffCorr[j] = AEff[j] * CorrFac[j]; //Analysing power for 12C based on correction factor from Mikhail
            AmpVal = PVal*AEff[j];

            if (i == 0){
                AsymmFit->SetLineColor(4);
                AsymmFit->FixParameter(1, 0.);
                AsymmHists[i][j]->Fit("AsymmFit", "Q");
                SinAmps[i][j] = AsymmFit->GetParameter(0);
                SinAmpErrs[i][j] = AsymmFit->GetParError(0);
            }

            if(i == 1){
                AsymmFit->SetLineColor(4);
                AsymmFit->SetLineStyle(10);
                AsymmFit->FixParameter(1, AmpVal);
                AsymmHists[i][j]->Fit("AsymmFit", "Q");
                SinAmps[i][j] = AsymmFit->GetParameter(0);
                SinAmpErrs[i][j] = AsymmFit->GetParError(0);
            }

        }
    }

    TFile f1("CM2_AsymmFits_PTotal_137_37_V1.root", "RECREATE");

    for(Int_t i = 0; i < 2; i++){ // Fit version
        for(Int_t j = 0; j < 4; j++){ // Energy
            AsymmHists[i][j]->Write();
        }
    }

    TCanvas *canvas20 = new TCanvas("canvas20","canvas20", 1920, 1080);
    canvas20->Divide(2, 2);
    for(int i = 1 ; i < 5 ; i++){
        canvas20->cd(i);
        AsymmHists[0][i-1]->Draw("EP");
        for(int j = 1 ; j < 2 ; j++){
            AsymmHists[j][i-1]->Draw("SAMEEP");
        }
    }

    canvas20->Write();

    f1.Write();
    f1.Close();

    TFile *f2= TFile::Open("/scratch/Mainz_Software/a2GoAT/CircPol_Aug16.root");

    for (Int_t m = 0; m < 2; m++){ // Fit Version
        for (Int_t n = 0; n < 4; n++){ //E

            Double_t EPoint = 450 + (n*100);

            Cx[m][n] = SinAmps[m][n]/(AEff[n]*(Graph->Eval(EPoint ,0)));
            CxErr[m][n] = fabs(SinAmpErrs[m][n]/(AEff[n]*(Graph->Eval(EPoint ,0))));
        }
    }


    TFile f3("CM2_Cx_Plots_137_37_V1.root", "RECREATE");

    double x[4] = {400, 500, 600, 700};
    double ex[4] = {50, 50, 50, 50};

    for(Int_t i = 0 ; i < 2 ; i++)
    {
        sprintf(name, "Cx_V%i", (i+1));
        CxPlots[i] = new TGraphErrors(4 , x, Cx[i], ex, CxErr[i]);
        CxPlots[i]->SetName(name);
        CxPlots[i]->SetTitle("C_{x}(E_{#gamma})");
        CxPlots[i]->SetMarkerColor(4);
        CxPlots[i]->SetLineColor(4);
        CxPlots[i]->SetMarkerStyle(8);
        CxPlots[i]->SetMarkerSize(1);
        CxPlots[i]->GetXaxis()->SetTitle("E_{#gamma}");
        CxPlots[i]->GetXaxis()->SetRangeUser(300, 800);
        CxPlots[i]->GetYaxis()->SetRangeUser(-1, 1);
        CxPlots[i]->GetYaxis()->SetTitle("C_{x}");
        CxPlots[i]->Write();
    }

    TCanvas *canvas = new TCanvas("canvas","canvas", 1920, 1080);
    CxPlots[0]->Draw("AEP");
    CxPlots[1]->SetMarkerColor(2);
    CxPlots[1]->SetLineColor(2);
    CxPlots[1]->Draw("SAMEEP");

    leg = new TLegend(0.75, 0.75, 0.9, 0.9);
    leg->AddEntry(CxPlots[0], "Fixed P = 0, #chi^{2}", "lepz");
    leg->AddEntry(CxPlots[1], "Fixed P, #chi^{2}", "lepz");

    leg->Draw("SAME");

    canvas->Write();

    f3.Write();

}
