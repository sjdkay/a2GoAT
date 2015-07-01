#include "./includes.h"

void macroTGraph(){

  double Theta[14] = {30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160};
  double CM150[14] = {4.4, 4.72, 4.94, 5.31, 5.02, 5.03, 4.44, 4.31, 4.14, 3.8, 3.49, 3.19, 2.96, 1.15};
  double CM250[14] = {5.15, 5.6, 5.88, 5.89, 6.07, 6.23, 5.93, 5.5, 5.37, 4.92, 4.62, 4.07, 4.15, 3.82};
  double CM350[14] = {3.08, 2.9, 3.1, 3.22, 3.12, 3.32, 3.11, 3.08, 2.83, 2.34, 2.34, 2.1, 2.1, 1.86};
  double CM450[14] = {1.01, 1.02, 1.04, 1.1, 1.22, 1.21, 1.14, 1.07, 0.99, 0.78, 0.8, 0.68, 0.59, 0.59};
  double CM550[13] = {0.62, 0.68, 0.66, 0.64, 0.62, 0.65, 0.67, 0.56, 0.43, 0.41, 0.33, 0.27, 0.2};
  double CM150Err[14] = {0.17, 0.15, 0.13, 0.13, 0.12, 0.12, 0.11, 0.11, 0.12, 0.12, 0.12, 0.13, 0.15, 0.14};
  double CM250Err[14] = {0.3, 0.24, 0.21, 0.19, 0.19, 0.18, 0.18, 0.17, 0.18, 0.18, 0.19, 0.2, 0.24, 0.33};
  double CM350Err[14] = {0.25, 0.16, 0.14, 0.13, 0.12, 0.12, 0.12, 0.12, 0.12, 0.11, 0.12, 0.13, 0.16, 0.2};
  double CM450Err[14] = {0.12, 0.07, 0.06, 0.06, 0.06, 0.06, 0.06, 0.06, 0.05, 0.05, 0.06, 0.06, 0.07, 0.09};
  double CM550Err[13] = {0.07, 0.06, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.04, 0.04, 0.04, 0.05, 0.05}; 
  TFile *f = new TFile("../PhysicsTotal11_24_03_15.root");
  TText *warn = new TText(0, 0 , "PRELIMINARY");

  TCanvas *canvas = new TCanvas("canvas","canvas",1000,10,550,400);
  
  f -> cd();

  TPad *pad1 = new TPad("pad1","",0,0,1,1);
  pad1->Draw();
  pad1->cd();
  
  pad1->SetTickx(1);
  pad1->SetTicky(1);
  pad1->SetGridx(1);
  pad1->SetGridy(1);
  
  TH1F  *hr;
  Char_t hrTitle[64];
  
  Float_t xMin = 0.0;
  Float_t xMax = 180;
  Float_t yMin = 0;
  Float_t yMax = 900;
  
  strcpy(hrTitle,"CM Angular Distribution for 150 +/- 50 MeV; CM Angle");
  hr = canvas->DrawFrame(xMin,yMin,xMax,yMax);
  hr->SetTitle(hrTitle);
  const Int_t nBins = 14;
  
  Float_t x[nBins];
  Float_t y[nBins];
  Float_t Ey[nBins];

 for( Int_t i = 0 ; i < nBins ; i++ ){
  
   x[i] = Theta[i];
   y[i] = CM150[i]*150;
   Ey[i] = CM150Err[i]*150;
  }

  TGraphErrors *gr;
  
  gr  = new TGraphErrors(nBins,x,y,0,Ey);    
  gr->SetMarkerColor(4);
  gr->SetMarkerSize(0.8);
  gr->Draw("P");
  gr->Draw("SAMEE1");
  CM_150MeV->SetMarkerStyle(1);
  CM_150MeV->SetLineColor(2);
  CM_150MeV->Rebin(2);
  CM_150MeV->Draw("EHISTSAMEP");

  warn->Draw();

  canvas->SaveAs("./CM150_Hist_and_Graph.pdf");

TCanvas *canvas = new TCanvas("canvas","canvas",1000,10,550,400);
  
  f -> cd();

  TPad *pad1 = new TPad("pad1","",0,0,1,1);
  pad1->Draw();
  pad1->cd();
  
  pad1->SetTickx(1);
  pad1->SetTicky(1);
  pad1->SetGridx(1);
  pad1->SetGridy(1);
  
  TH1F  *hr;
  Char_t hrTitle[64];
  
  Float_t xMin = 0.0;
  Float_t xMax = 180;
  Float_t yMin = 0;
  Float_t yMax = 1100;
  
  strcpy(hrTitle,"CM Angular Distribution for 250 +/- 50 MeV; CM Angle");
  hr = canvas->DrawFrame(xMin,yMin,xMax,yMax);
  hr->SetTitle(hrTitle);
  const Int_t nBins = 14;
  
  Float_t x[nBins];
  Float_t y[nBins];
  Float_t Ey[nBins];

 for( Int_t i = 0 ; i < nBins ; i++ ){
  
   x[i] = Theta[i];
   y[i] = CM250[i]*150;
   Ey[i] = CM250Err[i]*150;
  }

  TGraphErrors *gr;
  
  gr  = new TGraphErrors(nBins,x,y,0,Ey);    
  gr->SetMarkerColor(4);
  gr->SetMarkerSize(0.8);
  gr->Draw("P");
  gr->Draw("SAMEE1");
  CM_250MeV->SetMarkerStyle(1);
  CM_250MeV->SetLineColor(2);
  CM_250MeV->Rebin(2);
  CM_250MeV->Draw("EHISTSAMEP");

  warn->Draw();

  canvas->SaveAs("./CM250_Hist_and_Graph.pdf");

TCanvas *canvas = new TCanvas("canvas","canvas",1000,10,550,400);
  
  f -> cd();

  TPad *pad1 = new TPad("pad1","",0,0,1,1);
  pad1->Draw();
  pad1->cd();
  
  pad1->SetTickx(1);
  pad1->SetTicky(1);
  pad1->SetGridx(1);
  pad1->SetGridy(1);
  
  TH1F  *hr;
  Char_t hrTitle[64];
  
  Float_t xMin = 0.0;
  Float_t xMax = 180;
  Float_t yMin = 0;
  Float_t yMax = 500;
  
  strcpy(hrTitle,"CM Angular Distribution for 350 +/- 50 MeV; CM Angle");
  hr = canvas->DrawFrame(xMin,yMin,xMax,yMax);
  hr->SetTitle(hrTitle);
  const Int_t nBins = 14;
  
  Float_t x[nBins];
  Float_t y[nBins];
  Float_t Ey[nBins];

 for( Int_t i = 0 ; i < nBins ; i++ ){
  
   x[i] = Theta[i];
   y[i] = CM350[i]*125;
   Ey[i] = CM350Err[i]*125;
  }

  TGraphErrors *gr;
  
  gr  = new TGraphErrors(nBins,x,y,0,Ey);    
  gr->SetMarkerColor(4);
  gr->SetMarkerSize(0.8);
  gr->Draw("P");
  gr->Draw("SAMEE1");
  CM_350MeV->SetMarkerStyle(1);
  CM_350MeV->SetLineColor(2);
  CM_350MeV->Rebin(2);
  CM_350MeV->Draw("EHISTSAMEP");

  warn->Draw();

  canvas->SaveAs("./CM350_Hist_and_Graph.pdf");


TCanvas *canvas = new TCanvas("canvas","canvas",1000,10,550,400);
  
  f -> cd();

  TPad *pad1 = new TPad("pad1","",0,0,1,1);
  pad1->Draw();
  pad1->cd();
  
  pad1->SetTickx(1);
  pad1->SetTicky(1);
  pad1->SetGridx(1);
  pad1->SetGridy(1);
  
  TH1F  *hr;
  Char_t hrTitle[64];
  
  Float_t xMin = 0.0;
  Float_t xMax = 180;
  Float_t yMin = 0;
  Float_t yMax = 200;
  
  strcpy(hrTitle,"CM Angular Distribution for 450 +/- 50 MeV; CM Angle");
  hr = canvas->DrawFrame(xMin,yMin,xMax,yMax);
  hr->SetTitle(hrTitle);
  const Int_t nBins = 14;
  
  Float_t x[nBins];
  Float_t y[nBins];
  Float_t Ey[nBins];

 for( Int_t i = 0 ; i < nBins ; i++ ){
  
   x[i] = Theta[i];
   y[i] = CM450[i]*125;
   Ey[i] = CM450Err[i]*125;
  }

  TGraphErrors *gr;
  
  gr  = new TGraphErrors(nBins,x,y,0,Ey);    
  gr->SetMarkerColor(4);
  gr->SetMarkerSize(0.8);
  gr->Draw("P");
  gr->Draw("SAMEE1");
  CM_450MeV->SetMarkerStyle(1);
  CM_450MeV->SetLineColor(2);
  CM_450MeV->Rebin(2);
  CM_450MeV->Draw("EHISTSAMEP");

  warn->Draw();

  canvas->SaveAs("./CM450_Hist_and_Graph.pdf");

TCanvas *canvas = new TCanvas("canvas","canvas",1000,10,550,400);
  
  f -> cd();

  TPad *pad1 = new TPad("pad1","",0,0,1,1);
  pad1->Draw();
  pad1->cd();
  
  pad1->SetTickx(1);
  pad1->SetTicky(1);
  pad1->SetGridx(1);
  pad1->SetGridy(1);
  
  TH1F  *hr;
  Char_t hrTitle[64];
  
  Float_t xMin = 0.0;
  Float_t xMax = 180;
  Float_t yMin = 0;
  Float_t yMax = 60;
  
  strcpy(hrTitle,"CM Angular Distribution for 550 +/- 50 MeV; CM Angle");
  hr = canvas->DrawFrame(xMin,yMin,xMax,yMax);
  hr->SetTitle(hrTitle);
  const Int_t nBins = 13;
  
  Float_t x[nBins];
  Float_t y[nBins];
  Float_t Ey[nBins];

 for( Int_t i = 0 ; i < nBins ; i++ ){
  
   x[i] = Theta[i];
   y[i] = CM550[i]*75;
   Ey[i] = CM550Err[i]*75;
  }

  TGraphErrors *gr;
  
  gr  = new TGraphErrors(nBins,x,y,0,Ey);    
  gr->SetMarkerColor(4);
  gr->SetMarkerSize(0.8);
  gr->Draw("P");
  gr->Draw("SAMEE1");
  CM_550MeV->SetMarkerStyle(1);
  CM_550MeV->SetLineColor(2);
  CM_550MeV->Draw("EHISTSAMEP");

  warn->Draw(); 

  canvas->SaveAs("./CM550_Hist_and_Graph.pdf");

TCanvas *canvas = new TCanvas("canvas","canvas",1000,10,550,400);
  
  f -> cd();

  TPad *pad1 = new TPad("pad1","",0,0,1,1);
  pad1->Draw();
  pad1->cd();
  
  pad1->SetTickx(1);
  pad1->SetTicky(1);
  pad1->SetGridx(1);
  pad1->SetGridy(1);
  
  TH1F  *hr;
  Char_t hrTitle[64];
  
  Float_t xMin = 0.0;
  Float_t xMax = 180;
  Float_t yMin = 0;
  Float_t yMax = 60;
  
  strcpy(hrTitle,"CM Angular Distribution for 550 +/- 50 MeV; CM Angle");
  hr = canvas->DrawFrame(xMin,yMin,xMax,yMax);
  hr->SetTitle(hrTitle);
  const Int_t nBins = 13;
  
  Float_t x[nBins];
  Float_t y[nBins];
  Float_t Ey[nBins];

 for( Int_t i = 0 ; i < nBins ; i++ ){
  
   x[i] = Theta[i];
   y[i] = CM550[i]*75;
   Ey[i] = CM550Err[i]*75;
  }

  TGraphErrors *gr;
  
  gr  = new TGraphErrors(nBins,x,y,0,Ey);    
  gr->SetMarkerColor(4);
  gr->SetMarkerSize(0.8);
  CM_550MeV->SetMarkerStyle(1);
  CM_550MeV->SetLineColor(2);

  leg = new TLegend(0,0,1,1);
  leg->AddEntry(gr,"Previous Results","ep");
  leg->AddEntry(CM_550MeV,"Current Data","lep");
  leg->Draw();

  canvas->SaveAs("./legend.pdf");
}
