#include "./includes_Sigma.h"

void Sigma(){

  TFile *f1= TFile::Open("/scratch/Mainz_Software/a2GoAT/ParaPerpAsymm_Total_9.root");
  TTree *t1 = (TTree*)f1->Get("Parameter_Values");

  Double_t pValuesCM1[12], pValuesCM2[12], pValuesCM3[12], pValuesCM4[12], pValuesCM5[12], pValuesCM6[12];
  Double_t nValuesCM1[12], nValuesCM2[12], nValuesCM3[12], nValuesCM4[12], nValuesCM5[12], nValuesCM6[12];
  Double_t pErrValuesCM1[12], pErrValuesCM2[12], pErrValuesCM3[12], pErrValuesCM4[12], pErrValuesCM5[12], pErrValuesCM6[12];
  Double_t nErrValuesCM1[12], nErrValuesCM2[12], nErrValuesCM3[12], nErrValuesCM4[12], nErrValuesCM5[12], nErrValuesCM6[12];
  Double_t pSigmaValuesCM1[12], pSigmaValuesCM2[12], pSigmaValuesCM3[12], pSigmaValuesCM4[12], pSigmaValuesCM5[12], pSigmaValuesCM6[12];
  Double_t nSigmaValuesCM1[12], nSigmaValuesCM2[12], nSigmaValuesCM3[12], nSigmaValuesCM4[12], nSigmaValuesCM5[12], nSigmaValuesCM6[12];
  Double_t pSigmaErrValuesCM1[12], pSigmaErrValuesCM2[12], pSigmaErrValuesCM3[12], pSigmaErrValuesCM4[12], pSigmaErrValuesCM5[12], pSigmaErrValuesCM6[12];
  Double_t nSigmaErrValuesCM1[12], nSigmaErrValuesCM2[12], nSigmaErrValuesCM3[12], nSigmaErrValuesCM4[12], nSigmaErrValuesCM5[12], nSigmaErrValuesCM6[12];

  Double_t pCosAmpCM1, pCosAmpErrCM1, nCosAmpCM1, nCosAmpErrCM1;
  Double_t pCosAmpCM2, pCosAmpErrCM2, nCosAmpCM2, nCosAmpErrCM2;
  Double_t pCosAmpCM3, pCosAmpErrCM3, nCosAmpCM3, nCosAmpErrCM3;
  Double_t pCosAmpCM4, pCosAmpErrCM4, nCosAmpCM4, nCosAmpErrCM4;
  Double_t pCosAmpCM5, pCosAmpErrCM5, nCosAmpCM5, nCosAmpErrCM5;
  Double_t pCosAmpCM6, pCosAmpErrCM6, nCosAmpCM6, nCosAmpErrCM6;

  // Set branch addresses to get values from
  t1->SetBranchAddress("pCosAmpCM1", &pCosAmpCM1);
  t1->SetBranchAddress("pCosAmpErrCM1", &pCosAmpErrCM1);
  t1->SetBranchAddress("nCosAmpCM1", &nCosAmpCM1);
  t1->SetBranchAddress("nCosAmpErrCM1", &nCosAmpErrCM1);
  t1->SetBranchAddress("pCosAmpCM2", &pCosAmpCM2);
  t1->SetBranchAddress("pCosAmpErrCM2", &pCosAmpErrCM2);
  t1->SetBranchAddress("nCosAmpCM2", &nCosAmpCM2);
  t1->SetBranchAddress("nCosAmpErrCM2", &nCosAmpErrCM2);
  t1->SetBranchAddress("pCosAmpCM3", &pCosAmpCM3);
  t1->SetBranchAddress("pCosAmpErrCM3", &pCosAmpErrCM3);
  t1->SetBranchAddress("nCosAmpCM3", &nCosAmpCM3);
  t1->SetBranchAddress("nCosAmpErrCM3", &nCosAmpErrCM3);
  t1->SetBranchAddress("pCosAmpCM4", &pCosAmpCM4);
  t1->SetBranchAddress("pCosAmpErrCM4", &pCosAmpErrCM4);
  t1->SetBranchAddress("nCosAmpCM4", &nCosAmpCM4);
  t1->SetBranchAddress("nCosAmpErrCM4", &nCosAmpErrCM4);
  t1->SetBranchAddress("pCosAmpCM5", &pCosAmpCM5);
  t1->SetBranchAddress("pCosAmpErrCM5", &pCosAmpErrCM5);
  t1->SetBranchAddress("nCosAmpCM5", &nCosAmpCM5);
  t1->SetBranchAddress("nCosAmpErrCM5", &nCosAmpErrCM5);
  t1->SetBranchAddress("pCosAmpCM6", &pCosAmpCM6);
  t1->SetBranchAddress("pCosAmpErrCM6", &pCosAmpErrCM6);
  t1->SetBranchAddress("nCosAmpCM6", &nCosAmpCM6);
  t1->SetBranchAddress("nCosAmpErrCM6", &nCosAmpErrCM6);

  // Load values from tree and asign values back into an array
  for (Int_t k = 0; k < 12; k++){
    Parameter_Values->GetEntry(k);
    pValuesCM1[k] = pCosAmpCM1;
    pErrValuesCM1[k] = pCosAmpErrCM1;
    nValuesCM1[k] = nCosAmpCM1;
    nErrValuesCM1[k] = nCosAmpErrCM1;
    pValuesCM2[k] = pCosAmpCM2;
    pErrValuesCM2[k] = pCosAmpErrCM2;
    nValuesCM2[k] = nCosAmpCM2;
    nErrValuesCM2[k] = nCosAmpErrCM2;
    pValuesCM3[k] = pCosAmpCM3;
    pErrValuesCM3[k] = pCosAmpErrCM3;
    nValuesCM3[k] = nCosAmpCM3;
    nErrValuesCM3[k] = nCosAmpErrCM3;
    pValuesCM4[k] = pCosAmpCM4;
    pErrValuesCM4[k] = pCosAmpErrCM4;
    nValuesCM4[k] = nCosAmpCM4;
    nErrValuesCM4[k] = nCosAmpErrCM4;
    pValuesCM5[k] = pCosAmpCM5;
    pErrValuesCM5[k] = pCosAmpErrCM5;
    nValuesCM5[k] = nCosAmpCM5;
    nErrValuesCM5[k] = nCosAmpErrCM5;
    pValuesCM5[k] = pCosAmpCM6;
    pErrValuesCM5[k] = pCosAmpErrCM6;
    nValuesCM5[k] = nCosAmpCM6;
    nErrValuesCM5[k] = nCosAmpErrCM6;
  }

  TFile *f2= TFile::Open("/scratch/Mainz_Software/a2GoAT/LinPol_Aug16.root"); // Open linear polarisation plot

  // Calculate values of sigma for each angular and energy bin
  for (Int_t i = 0; i < 12; i++){
    
    pSigmaValuesCM1[i] = pValuesCM1[i]/(Graph->Eval((410+(i*20)),0));
    pSigmaErrValuesCM1[i] = pErrValuesCM1[i]/(Graph->Eval((410+(i*20)),0));
    pSigmaValuesCM2[i] = pValuesCM2[i]/(Graph->Eval((410+(i*20)),0));
    pSigmaErrValuesCM2[i] = pErrValuesCM2[i]/(Graph->Eval((410+(i*20)),0));
    pSigmaValuesCM3[i] = pValuesCM3[i]/(Graph->Eval((410+(i*20)),0));
    pSigmaErrValuesCM3[i] = pErrValuesCM3[i]/(Graph->Eval((410+(i*20)),0));
    pSigmaValuesCM4[i] = pValuesCM4[i]/(Graph->Eval((410+(i*20)),0));
    pSigmaErrValuesCM4[i] = pErrValuesCM4[i]/(Graph->Eval((410+(i*20)),0));
    pSigmaValuesCM5[i] = pValuesCM5[i]/(Graph->Eval((410+(i*20)),0));
    pSigmaErrValuesCM5[i] = pErrValuesCM5[i]/(Graph->Eval((410+(i*20)),0));
    pSigmaValuesCM6[i] = pValuesCM6[i]/(Graph->Eval((410+(i*20)),0));
    pSigmaErrValuesCM6[i] = pErrValuesCM6[i]/(Graph->Eval((410+(i*20)),0));

    nSigmaValuesCM1[i] = nValuesCM1[i]/(Graph->Eval((410+(i*20)),0));
    nSigmaErrValuesCM1[i] = nErrValuesCM1[i]/(Graph->Eval((410+(i*20)),0));
    nSigmaValuesCM2[i] = nValuesCM2[i]/(Graph->Eval((410+(i*20)),0));
    nSigmaErrValuesCM2[i] = nErrValuesCM2[i]/(Graph->Eval((410+(i*20)),0));
    nSigmaValuesCM3[i] = nValuesCM3[i]/(Graph->Eval((410+(i*20)),0));
    nSigmaErrValuesCM3[i] = nErrValuesCM3[i]/(Graph->Eval((410+(i*20)),0));
    nSigmaValuesCM4[i] = nValuesCM4[i]/(Graph->Eval((410+(i*20)),0));
    nSigmaErrValuesCM4[i] = nErrValuesCM4[i]/(Graph->Eval((410+(i*20)),0));
    nSigmaValuesCM5[i] = nValuesCM5[i]/(Graph->Eval((410+(i*20)),0));
    nSigmaErrValuesCM5[i] = nErrValuesCM5[i]/(Graph->Eval((410+(i*20)),0));
    nSigmaValuesCM6[i] = nValuesCM6[i]/(Graph->Eval((410+(i*20)),0));
    nSigmaErrValuesCM6[i] = nErrValuesCM6[i]/(Graph->Eval((410+(i*20)),0));

  }

  TFile f3("Sigma_Plots.root", "RECREATE");

  Float_t xMin = 400;
  Float_t xMax = 650;
  Float_t yMin = -5;
  Float_t yMax = 5;
  Double_t x[12] = {410, 430, 450, 470, 490, 510, 530, 550, 570, 590, 610, 630};
  Double_t ex[12] = {10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10};

  TCanvas *canvas = new TCanvas("canvas","canvas", 1920, 1080);
  TPad *pad1 = new TPad("pad1","",0,0,1,1);
  pad1->Draw();
  pad1->cd();
    
  pad1->SetTickx(1);
  pad1->SetTicky(1);
  pad1->SetGridx(1);
  pad1->SetGridy(1);
  TH1F  *hr;
  hr = canvas->DrawFrame(xMin,yMin,xMax,yMax);

  gr = new TGraphErrors(12, x, pSigmaValuesCM1 , ex, pSigmaErrValuesCM1);
  gr->SetMarkerColor(2);
  gr->SetMarkerStyle(21);
  gr->SetMarkerSize(2);
  gr2 = new TGraphErrors(12, x, nSigmaValuesCM1, ex, nSigmaErrValuesCM1);
  gr2->SetMarkerColor(4);
  gr2->SetMarkerStyle(22);
  gr2->SetMarkerSize(2);
  leg = new TLegend(0.75, 0.75, 0.95, 0.95); // Add legend to plot
  leg->AddEntry(gr, "Proton", "ep");
  leg->AddEntry(gr2, "Neutron", "ep");
  gr->Draw("E1P");
  gr2->Draw("ESAMEP");
  leg->Draw("Same");

  TCanvas *canvas1 = new TCanvas("canvas1","canvas1", 1920, 1080);
  TPad *pad2 = new TPad("pad2","",0,0,1,1);
  pad2->Draw();
  pad2->cd();

  pad2->SetTickx(1);
  pad2->SetTicky(1);
  pad2->SetGridx(1);
  pad2->SetGridy(1);
  TH1F  *hr1;
  hr1 = canvas1->DrawFrame(xMin,yMin,xMax,yMax);

  gr3 = new TGraphErrors(12, x, pSigmaValuesCM2 , ex, pSigmaErrValuesCM2);
  gr3->SetMarkerColor(2);
  gr3->SetMarkerStyle(21);
  gr3->SetMarkerSize(2);
  gr4 = new TGraphErrors(12, x, nSigmaValuesCM2, ex, nSigmaErrValuesCM2);
  gr4->SetMarkerColor(4);
  gr4->SetMarkerStyle(22);
  gr4->SetMarkerSize(2);
  leg2 = new TLegend(0.75, 0.75, 0.95, 0.95); // Add legend to plot
  leg2->AddEntry(gr3, "Proton", "ep");
  leg2->AddEntry(gr4, "Neutron", "ep");
  gr3->Draw("E1P");
  gr4->Draw("ESAMEP");
  leg2->Draw("Same");

  TCanvas *canvas2 = new TCanvas("canvas2","canvas2", 1920, 1080);
  TPad *pad3 = new TPad("pad3","",0,0,1,1);
  pad3->Draw();
  pad3->cd();
    
  pad3->SetTickx(1);
  pad3->SetTicky(1);
  pad3->SetGridx(1);
  pad3->SetGridy(1);
  TH1F  *hr2;
  hr2 = canvas2->DrawFrame(xMin,yMin,xMax,yMax);

  gr5 = new TGraphErrors(12, x, pSigmaValuesCM3 , ex, pSigmaErrValuesCM3);
  gr5->SetMarkerColor(2);
  gr5->SetMarkerStyle(21);
  gr5->SetMarkerSize(2);
  gr6 = new TGraphErrors(12, x, nSigmaValuesCM3, ex, nSigmaErrValuesCM3);
  gr6->SetMarkerColor(4);
  gr6->SetMarkerStyle(22);
  gr6->SetMarkerSize(2);
  leg3 = new TLegend(0.75, 0.75, 0.95, 0.95); // Add legend to plot
  leg3->AddEntry(gr5, "Proton", "ep");
  leg3->AddEntry(gr6, "Neutron", "ep");
  gr5->Draw("E1P");
  gr6->Draw("ESAMEP");
  leg3->Draw("Same");

  TCanvas *canvas3 = new TCanvas("canvas3","canvas3", 1920, 1080);
  TPad *pad4 = new TPad("pad4","",0,0,1,1);
  pad4->Draw();
  pad4->cd();
    
  pad4->SetTickx(1);
  pad4->SetTicky(1);
  pad4->SetGridx(1);
  pad4->SetGridy(1);
  TH1F  *hr3;
  hr3 = canvas3->DrawFrame(xMin,yMin,xMax,yMax);

  gr7 = new TGraphErrors(12, x, pSigmaValuesCM4 , ex, pSigmaErrValuesCM4);
  gr7->SetMarkerColor(2);
  gr7->SetMarkerStyle(21);
  gr7->SetMarkerSize(2);
  gr8 = new TGraphErrors(12, x, nSigmaValuesCM4, ex, nSigmaErrValuesCM4);
  gr8->SetMarkerColor(4);
  gr8->SetMarkerStyle(22);
  gr8->SetMarkerSize(2);
  leg4 = new TLegend(0.75, 0.75, 0.95, 0.95); // Add legend to plot
  leg4->AddEntry(gr7, "Proton", "ep");
  leg4->AddEntry(gr8, "Neutron", "ep");
  gr7->Draw("E1P");
  gr8->Draw("ESAMEP");
  leg4->Draw("Same");

  TCanvas *canvas4 = new TCanvas("canvas4","canvas4", 1920, 1080);
  TPad *pad5 = new TPad("pad5","",0,0,1,1);
  pad5->Draw();
  pad5->cd();

  pad5->SetTickx(1);
  pad5->SetTicky(1);
  pad5->SetGridx(1);
  pad5->SetGridy(1);
  TH1F  *hr4;
  hr4 = canvas4->DrawFrame(xMin,yMin,xMax,yMax);

  gr9 = new TGraphErrors(12, x, pSigmaValuesCM5 , ex, pSigmaErrValuesCM5);
  gr9->SetMarkerColor(2);
  gr9->SetMarkerStyle(21);
  gr9->SetMarkerSize(2);
  gr10 = new TGraphErrors(12, x, nSigmaValuesCM5, ex, nSigmaErrValuesCM5);
  gr10->SetMarkerColor(4);
  gr10->SetMarkerStyle(22);
  gr10->SetMarkerSize(2);
  leg5 = new TLegend(0.75, 0.75, 0.95, 0.95); // Add legend to plot
  leg5->AddEntry(gr9, "Proton", "ep");
  leg5->AddEntry(gr10, "Neutron", "ep");
  gr9->Draw("E1P");
  gr10->Draw("ESAMEP");
  leg5->Draw("Same");

  TCanvas *canvas5 = new TCanvas("canvas5","canvas5", 1920, 1080);
  TPad *pad6 = new TPad("pad6","",0,0,1,1);
  pad6->Draw();
  pad6->cd();
    
  pad6->SetTickx(1);
  pad6->SetTicky(1);
  pad6->SetGridx(1);
  pad6->SetGridy(1);
  TH1F  *hr5;
  hr5 = canvas5->DrawFrame(xMin,yMin,xMax,yMax);

  gr11 = new TGraphErrors(12, x, pSigmaValuesCM6 , ex, pSigmaErrValuesCM6);
  gr11->SetMarkerColor(2);
  gr11->SetMarkerStyle(21);
  gr11->SetMarkerSize(2);
  gr12 = new TGraphErrors(12, x, nSigmaValuesCM6, ex, nSigmaErrValuesCM6);
  gr12->SetMarkerColor(4);
  gr12->SetMarkerStyle(22);
  gr12->SetMarkerSize(2);
  leg6 = new TLegend(0.75, 0.75, 0.95, 0.95); // Add legend to plot
  leg6->AddEntry(gr11, "Proton", "ep");
  leg6->AddEntry(gr12, "Neutron", "ep");
  gr11->Draw("E1P");
  gr12->Draw("ESAMEP");
  leg6->Draw("Same");

  canvas->Write();
  canvas1->Write();
  canvas2->Write();
  canvas3->Write();
  canvas4->Write();
  canvas5->Write();
  f3.Write();

}
