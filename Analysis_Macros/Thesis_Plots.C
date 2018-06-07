#include "./includes.h"

void Thesis_Plots() {

    gStyle->SetTitleFontSize(0.06);

    TFile *f = new TFile("/scratch/Mainz_Software/a2GoAT/Cx_143_148_Combined.root"); // Open the latest PTotal file to load histograms from
    TFile f1("Cx_143_148_Plots_Combo.root", "RECREATE");

    char name[60];
    TGraphErrors* CxPlots3Bins[3][8]; // As fn of Angular bin for fixed energies
    TGraphErrors* CxPlots5Bins[5][8]; // As fn of Angular bin for fixed energies

    for(Int_t i = 0; i < 3; i++){ // Fit version
        for(Int_t j = 0; j < 8; j++){ // Energy
                sprintf(name, "Cx_%i_V%i_3Theta", 250+(j*100), i+1);
                CxPlots3Bins[i][j] = (TGraphErrors*) ((TGraphErrors*)f->Get(name));
        }
    }

    for(Int_t i = 0; i < 3; i++){ // Fit version
        for(Int_t j = 0; j < 8; j++){ // Energy
                sprintf(name, "Cx_%i_V%i_5Theta", 250+(j*100), i+1);
                CxPlots5Bins[i][j] = (TGraphErrors*) ((TGraphErrors*)f->Get(name));
                if(i == 0) {CxPlots5Bins[i][j]->SetMarkerStyle(24);}
                else if (i == 1) {CxPlots5Bins[i][j]->SetMarkerStyle(25);}
                else if (i == 2) {CxPlots5Bins[i][j]->SetMarkerStyle(26);}
        }
    }

    TCanvas *canvas = new TCanvas("canvas","canvas", 2560, 1440);
    canvas->Divide(5, 2, 0.000000001, 0.01);

    for(Int_t i = 0; i < 8; i++){
        canvas->cd(i+1);
        CxPlots3Bins[0][i]->Draw("AEP");
        CxPlots5Bins[0][i]->Draw("SAMEEP");

    }

    canvas->cd(9);

    leg = new TLegend(0.3, 0.3, 0.6, 0.6);
    leg->AddEntry(CxPlots3Bins[0][0], "Fixed p_{y} = 0, 3#theta Bins", "lepz");
    leg->AddEntry(CxPlots5Bins[0][0], "Fixed p_{y} = 0, 5#theta Bins", "lepz");

    leg->Draw("SAME");

    canvas->Write();

    double GlisterCx250x[6] = {cos(20*TMath::DegToRad()), cos(30*TMath::DegToRad()), cos(50*TMath::DegToRad()), cos(70*TMath::DegToRad()), cos(90*TMath::DegToRad()), cos(100*TMath::DegToRad())};
    double GlisterCx250y[6] = {0., 0., -0.1, -0.2,  -0.2, -0.3};
    double GlisterCx250yErr[6] = {0.01, 0.01, 0.01, 0.01, 0.02, 0.15};
    double GlisterCx350x[8] = {cos(20*TMath::DegToRad()), cos(30*TMath::DegToRad()), cos(40*TMath::DegToRad()), cos(50*TMath::DegToRad()), cos(70*TMath::DegToRad()), cos(90*TMath::DegToRad()), cos(110*TMath::DegToRad()), cos(120*TMath::DegToRad())};
    double GlisterCx350y[8] = {0.05, 0., -0.02, -0.05, -0.2, -0.2, 0.15, -0.05};
    double GlisterCx350yErr[8] = {0.01, 0.01, 0.01, 0.01, 0.01, 0.02, 0.1, 0.1};

    GlisterCx250 = new TGraphErrors(6, GlisterCx250x, GlisterCx250y, 0, GlisterCx250yErr);
    GlisterCx250->SetMarkerColor(618);
    GlisterCx250->SetLineColor(617);
    GlisterCx250->SetMarkerStyle(33);
    GlisterCx250->SetMarkerSize(2);
    GlisterCx250->SetName("GlisterCx250");
    GlisterCx250->Write();

    GlisterCx350 = new TGraphErrors(8, GlisterCx350x, GlisterCx350y, 0, GlisterCx350yErr);
    GlisterCx350->SetMarkerColor(618);
    GlisterCx350->SetLineColor(618);
    GlisterCx350->SetMarkerStyle(33);
    GlisterCx350->SetMarkerSize(2);
    GlisterCx350->SetName("GlisterCx350");
    GlisterCx350->Write();

    TCanvas *canvas6 = new TCanvas("canvas6", "canvas6", 2560, 1440);

    CxPlots3Bins[0][1]->Draw("AEP");
    CxPlots5Bins[0][1]->Draw("SAMEEP");
    GlisterCx350->Draw("SAMEEP");

    leg2 = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg2->AddEntry(CxPlots3Bins[0][1], "Fixed p_{y} = 0, 3#theta Bins", "lepz");
    leg2->AddEntry(CxPlots5Bins[0][1], "Fixed p_{y} = 0, 5#theta Bins", "lepz");
    leg2->AddEntry(GlisterCx350, "J. Glister PLB11 Results", "lepz");
    leg2->Draw("SAME");

    canvas6->Write();

    f1.Write();

}
