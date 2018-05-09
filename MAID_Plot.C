#include "./includes.h"

void MAID_Plot(){

    // Values for MAID analysis of photon pi0 p/n at Q^2 = 0 (real photon), polarisation = 1
    // W = sqrt(s) = 1396 for 570 MeV photon in gamma N system

    double x[19] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180};
    double ex[19] = {5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5};
    double Pyp[19] = {0, 0.1374, 0.2727, 0.4017, 0.5166, 0.6083, 0.6698, 0.6991, 0.7006, 0.6819, 0.6498, 0.5251, 0.3611, 0.0642, -0.2138, -0.2599, -0.1509, 0};
    double Pyn[19] = {0.0, 0.0493, 0.0982, 0.1463, 0.1932, 0.2378, 0.2791, 0.3167, 0.3522, 0.3885, 0.4293, 0.4777, 0.5351, 0.5989, 0.6575, 0.6818, 0.6146, 0.3857, 0.0};

    PypPlot = new TGraphErrors(19 , x, Pyp, ex, 0);
    PynPlot = new TGraphErrors(19 , x, Pyn, ex, 0);
    PypPlot->SetLineColor(2);
    PypPlot->SetMarkerColor(2);
    PypPlot->SetMarkerStyle(8);
    PynPlot->SetLineColor(4);
    PynPlot->SetMarkerColor(4);
    PynPlot->SetMarkerStyle(8);
    PypPlot->GetXaxis()->SetRangeUser(0, 180);
    PypPlot->GetYaxis()->SetRangeUser(-1, 1);
    PypPlot->SetTitle("P_{y'}(#theta) for  W = 1396MeV #gammaN #rightarrow N#pi^{0} MAID");
    PypPlot->GetXaxis()->SetTitle("#theta");
    PypPlot->GetYaxis()->SetTitle("P_{y'}");

    leg = new TLegend(0.8, 0.8, 0.9, 0.9);
    leg->AddEntry(PypPlot, "#gamma p #rightarrow p #pi^{0}", "lepz");
    leg->AddEntry(PynPlot, "#gamma n #rightarrow n #pi^{0}", "lepz");

    TCanvas *canvas1 = new TCanvas("canvas1","canvas1", 1920, 1080);
    TPad *pad1 = new TPad("pad1","",0,0,1,1);
    pad1->Draw();
    pad1->cd();

    PypPlot->Draw("AEP");
    PynPlot->Draw("SAMEEP");
    leg->Draw("SAME");

}
