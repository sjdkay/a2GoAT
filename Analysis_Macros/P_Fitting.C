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
    double x[10] = {320, 360, 400, 440, 480, 520, 560, 600, 640, 680};
    double ex[10] = {20, 20, 20, 20, 20, 20, 20, 20, 20, 20};
    double xLiu[19]={172, 200, 248, 299, 228, 237, 249, 296, 346, 401, 446, 305, 355, 407, 460, 341, 352, 401, 436};
    double xLiuErr[19]={13, 7, 8, 13, 7, 7, 8 ,10, 12, 15, 19, 9, 11, 13, 16, 10, 10, 12, 14};
    double xIkeda[19]={449.2, 501.7, 551.1, 603.5, 403.3, 446.6, 503.3, 547.6, 599, 511.4, 549, 597.8, 452.3, 503.3, 549.2, 597.9, 646.2, 500.6, 548.6};
    double PrevDataLiu[19]={-0.41, -0.25, -0.21, -0.01, -0.01, -0.04, -0.05, -0.13, -0.22, -0.47, -0.45, -0.21, -0.35, -0.52, -0.59, -0.23, -0.23, -0.44, -0.58};
    double PrevDataIkeda[19]={-0.34, -0.24, -0.34, 0.1, -0.44, -0.54, -0.61, -0.62, -0.42, -0.77, -0.51, -0.58, -0.35, -0.46, -0.54, -0.59, -0.38, -0.35, -0.55};
    double PrevDataLiuErr[19]={0.14, 0.12, 0.13, 0.13, 0.07, 0.11, 0.08, 0.07, 0.08, 0.08, 0.19, 0.06, 0.09, 0.09, 0.13, 0.07, 0.1, 0.06, 0.08};
    double PrevDataIkedaErr[19]={0.08, 0.09, 0.1, 0.13, 0.09, 0.1, 0.11, 0.11, 0.11, 0.13, 0.08, 0.12, 0.1, 0.11, 0.15, 0.24, 0.2};
    double Amp[10];
    double AmpErr[10];
    double Y_Off[10];
    double Y_OffErr[10];
    double APow = 0.1;
    double Pol[10];
    double PolErr[10];
    double Y_OffCorr[10];
    double Y_OffCorrErr[10];
    double P1;
    double P2;
    double P3;
    double Phi;
    double PolyVal;
    double F;
    double BinValue;
    double AdjBinValue;
    char name[60];

    TF1 *CosFit = new TF1("CosFit",  fitf, -3.0, 3.0, 2); //Give a name and range to the fitting funcion

    CosFit->SetParNames("Y_Offset", "Amplitdue"); //Name the parameters
    TF1 *Pn90CM = new TF1("Pn90CM", "1.64576-2.95484*(x/1000)+0.684577*(x/1000)**2-0.65*90**2/4/((x-560)**2+90**2/4)+(5.32305-35.3819*(x/1000)+70.145*(x/1000)**2-44.2899*(x/1000)**3)",300,700);

    // Add all relevant histograms to a list that can then be looped over
    TFile *f = new TFile("Physics_Total_119_03_5_18.root"); // Open the latest PTotal file to load histograms from
    TText *warn = new TText(0, 0 ,"PRELIMINARY"); // Preliminary warning label text
    TList *PhiScList = new TList;

    for (Int_t k = 0; k < 10; k++){
        sprintf(name, "PhiSc_%iMeV", 320+(k*40));
        PhiScList->Add((TH1F*)f->Get(name));
    }

    TList *PhiScFitList = new TList;

    for (Int_t k = 0; k < 10; k++){

        TH1* hist = PhiScList->At(k);
        hist->Fit("CosFit", "Q");
        Amp[k] = CosFit->GetParameter(1); // Add values of the fit to an array
        AmpErr[k] = CosFit->GetParError(1);
        Y_Off[k]  = CosFit->GetParameter(0);
        Y_OffErr[k] = CosFit->GetParError(0);
        Pol[k] = Amp[k]/APow;
        PolErr[k] = AmpErr[k]/APow;
        PhiScFitList->Add(hist);
    }

    TFile f1("Py_119.root", "RECREATE");

    TCanvas *canvas = new TCanvas("canvas","canvas", 1920, 1080);
    TPad *pad = new TPad("pad","",0,0,1,1);
    pad->Draw();
    pad->cd();

    TH1F  *hr;
    hr = canvas->DrawFrame(300, -1, 700, 1);
    hr->SetTitle("P_{y} as fn of E_{#gamma} (#theta_{nCM} = 90)");

    gr = new TGraphErrors(10, x, Pol, ex, PolErr);
    gr->SetMarkerColor(2);
    gr->SetMarkerStyle(5);
    gr->SetMarkerSize(2);
    gr->SetTitle("P_{y} as fn of E_{#gamma} (#theta_{nCM} = 90)");
    gr->GetXaxis()->SetTitle("E_{#gamma}");
    gr->GetYaxis()->SetTitle("P_{y}");
    gr->SetName("P300700");
    gr->Draw("ep");
    Pn90CM->Draw("SAME");

    gr->Write();
    canvas->Write();
    PhiScFitList->Write();

    f1.Write();
}
