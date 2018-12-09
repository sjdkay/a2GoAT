#include "./includes.h"

// define a function with 3 parameters
Double_t fitf(Double_t *x,Double_t *par)
{
    Double_t fitval = 0;
    fitval =  par[0] * (1+(par[1]*cos(x[0])));
    return fitval;
}

void P_Fitting(){

    // Define a bunch of arrays to be used later
    double x[12] = {200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750};
    double ex[12] = {25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25};
    double xLiu[19]={172, 200, 248, 299, 228, 237, 249, 296, 346, 401, 446, 305, 355, 407, 460, 341, 352, 401, 436};
    double xLiuErr[19]={13, 7, 8, 13, 7, 7, 8 ,10, 12, 15, 19, 9, 11, 13, 16, 10, 10, 12, 14};
    double xIkeda[19]={449.2, 501.7, 551.1, 603.5, 403.3, 446.6, 503.3, 547.6, 599, 511.4, 549, 597.8, 452.3, 503.3, 549.2, 597.9, 646.2, 500.6, 548.6};
    double PrevDataLiu[19]={-0.41, -0.25, -0.21, -0.01, -0.01, -0.04, -0.05, -0.13, -0.22, -0.47, -0.45, -0.21, -0.35, -0.52, -0.59, -0.23, -0.23, -0.44, -0.58};
    double PrevDataIkeda[19]={-0.34, -0.24, -0.34, 0.1, -0.44, -0.54, -0.61, -0.62, -0.42, -0.77, -0.51, -0.58, -0.35, -0.46, -0.54, -0.59, -0.38, -0.35, -0.55};
    double PrevDataLiuErr[19]={0.14, 0.12, 0.13, 0.13, 0.07, 0.11, 0.08, 0.07, 0.08, 0.08, 0.19, 0.06, 0.09, 0.09, 0.13, 0.07, 0.1, 0.06, 0.08};
    double PrevDataIkedaErr[19]={0.08, 0.09, 0.1, 0.13, 0.09, 0.1, 0.11, 0.11, 0.11, 0.13, 0.08, 0.12, 0.1, 0.11, 0.15, 0.24, 0.2};
    double Amp[12];
    double AmpErr[12];
    double Y_Off[12];
    double Y_OffErr[12];
    double AEff[12];
    double Pol[12];
    double PolErr[12];
    double Y_OffCorr[12];
    double Y_OffCorrErr[12];
    double P1;
    double P2;
    double P3;
    double Phi;
    double PolyVal;
    double F;
    double BinValue;
    double AdjBinValue;
    char name[60];
    char NeutronEThetaScName[60];
    double MeanX[6][3];
    double MeanY[6][3];

    TH1F* FitHists[12];

    TF1 *CosFit = new TF1("CosFit",  fitf, -3.0, 3.0, 2); //Give a name and range to the fitting funcion

    CosFit->SetParNames("Y_Offset", "Amplitdue"); //Name the parameters
    TF1 *Pn90CM = new TF1("Pn90CM", "1.64576-2.95484*(x/1000)+0.684577*(x/1000)**2-0.65*90**2/4/((x-560)**2+90**2/4)+(5.32305-35.3819*(x/1000)+70.145*(x/1000)**2-44.2899*(x/1000)**3)",300,700);

    // Add all relevant histograms to a list that can then be looped over
    TFile *f = new TFile("Physics_Total_Amo140_Lin40_Combined.root"); // Open the latest PTotal file to load histograms from
    TFile *fAy = new TFile ("/scratch/Mainz_Software/a2GoAT/npAy.root");
    TFile f1("Py_140.root", "RECREATE");

    for (Int_t k = 0; k < 12; k++){
        double ECent = (200 + (k*50))/2; // Assume Neutron takes HALF energy as it is going at roughly 90 deg into Polarimeter
        sprintf(name, "PhiSc%i", 200+(k*50));
        FitHists[k] = (TH1F*) ((TH1F*)f->Get(name));
        FitHists[k]->Fit("CosFit", "Q");
        Amp[k] = CosFit->GetParameter(1); // Add values of the fit to an array
        AmpErr[k] = CosFit->GetParError(1);
        Y_Off[k]  = CosFit->GetParameter(0);
        Y_OffErr[k] = CosFit->GetParError(0);
        AEff[k] = ((TH2F*)fAy->Get("nppnAy"))->Interpolate(ECent, 22.5);
        Pol[k] = Amp[k]/AEff[k];
        PolErr[k] = AmpErr[k]/AEff[k];
        FitHists[k]->Write();
        cout << Amp[k] << endl;
    }

    gr = new TGraphErrors(12, x, Pol, ex, PolErr);
    gr->SetMarkerColor(2);
    gr->SetMarkerStyle(5);
    gr->SetMarkerSize(2);
    gr->GetXaxis()->SetRangeUser(150, 800);
    gr->GetYaxis()->SetRangeUser(-30, 10);
    gr->SetTitle("P_{y}(E_{#gamma}) (#theta_{nCM} = 90)");
    gr->GetXaxis()->SetTitle("E_{#gamma}");
    gr->GetYaxis()->SetTitle("p_{y}");
    gr->SetName("py200750");
    gr->GetXaxis()->SetLabelSize(0.06);
    gr->GetXaxis()->SetTitleSize(0.06);
    gr->GetXaxis()->SetTitleOffset(0.7);
    gr->GetXaxis()->CenterTitle();
    gr->GetYaxis()->SetLabelSize(0.06);
    gr->GetYaxis()->SetTitleSize(0.08);
    gr->GetYaxis()->SetTitleOffset(0.5);
    gr->GetYaxis()->CenterTitle();
    gr->Draw("AEP");
    Pn90CM->Draw("SAME");

    gr->Write();
    f1.Write();
}
