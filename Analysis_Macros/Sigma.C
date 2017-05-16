#include "./includes_Sigma.h"

void Sigma(){

  TFile *f1= TFile::Open("/scratch/Mainz_Software/a2GoAT/ParaPerpAsymm_Total_12.root");
  TTree *t1 = (TTree*)f1->Get("Parameter_Values");

  Double_t pValues410[10], pValues430[10], pValues450[10], pValues470[10], pValues490[10], pValues510[10], pValues530[10], pValues550[10], pValues570[10], pValues590[10], pValues610[10], pValues630[10];
  Double_t pErrValues410[10], pErrValues430[10], pErrValues450[10], pErrValues470[10], pErrValues490[10], pErrValues510[10], pErrValues530[10], pErrValues550[10], pErrValues570[10], pErrValues590[10], pErrValues610[10], pErrValues630[10];
  Double_t pSigmaValues410[10], pSigmaValues430[10], pSigmaValues450[10], pSigmaValues470[10], pSigmaValues490[10], pSigmaValues510[10], pSigmaValues530[10], pSigmaValues550[10], pSigmaValues570[10], pSigmaValues590[10], pSigmaValues610[10], pSigmaValues630[10];
  Double_t pSigmaErrValues410[10], pSigmaErrValues430[10], pSigmaErrValues450[10], pSigmaErrValues470[10], pSigmaErrValues490[10], pSigmaErrValues510[10], pSigmaErrValues530[10], pSigmaErrValues550[10], pSigmaErrValues570[10], pSigmaErrValues590[10], pSigmaErrValues610[10], pSigmaErrValues630[10];

  Double_t pCosAmp410, pCosAmpErr410;
  Double_t pCosAmp430, pCosAmpErr430;
  Double_t pCosAmp450, pCosAmpErr450;
  Double_t pCosAmp470, pCosAmpErr470;
  Double_t pCosAmp490, pCosAmpErr490;
  Double_t pCosAmp510, pCosAmpErr510;
  Double_t pCosAmp530, pCosAmpErr530;
  Double_t pCosAmp550, pCosAmpErr550;
  Double_t pCosAmp570, pCosAmpErr570;
  Double_t pCosAmp590, pCosAmpErr590;
  Double_t pCosAmp610, pCosAmpErr610;
  Double_t pCosAmp630, pCosAmpErr630;

  // Set branch addresses to get values from
  t1->SetBranchAddress("pCosAmp410", &pCosAmp410);
  t1->SetBranchAddress("pCosAmpErr410", &pCosAmpErr410);
  t1->SetBranchAddress("pCosAmp430", &pCosAmp430);
  t1->SetBranchAddress("pCosAmpErr430", &pCosAmpErr430);
  t1->SetBranchAddress("pCosAmp450", &pCosAmp450);
  t1->SetBranchAddress("pCosAmpErr450", &pCosAmpErr450);
  t1->SetBranchAddress("pCosAmp470", &pCosAmp470);
  t1->SetBranchAddress("pCosAmpErr470", &pCosAmpErr470);
  t1->SetBranchAddress("pCosAmp490", &pCosAmp490);
  t1->SetBranchAddress("pCosAmpErr490", &pCosAmpErr490);
  t1->SetBranchAddress("pCosAmp510", &pCosAmp510);
  t1->SetBranchAddress("pCosAmpErr510", &pCosAmpErr510);
  t1->SetBranchAddress("pCosAmp530", &pCosAmp530);
  t1->SetBranchAddress("pCosAmpErr530", &pCosAmpErr530);
  t1->SetBranchAddress("pCosAmp550", &pCosAmp550);
  t1->SetBranchAddress("pCosAmpErr550", &pCosAmpErr550);
  t1->SetBranchAddress("pCosAmp570", &pCosAmp570);
  t1->SetBranchAddress("pCosAmpErr570", &pCosAmpErr570);
  t1->SetBranchAddress("pCosAmp590", &pCosAmp590);
  t1->SetBranchAddress("pCosAmpErr590", &pCosAmpErr590);
  t1->SetBranchAddress("pCosAmp610", &pCosAmp610);
  t1->SetBranchAddress("pCosAmpErr610", &pCosAmpErr610);
  t1->SetBranchAddress("pCosAmp630", &pCosAmp630);
  t1->SetBranchAddress("pCosAmpErr630", &pCosAmpErr630);

  // Load values from tree and asign values back into an array
  for (Int_t k = 0; k < 10; k++){
    Parameter_Values->GetEntry(k);
    pValues410[k] = pCosAmp410;
    pErrValues410[k] = pCosAmpErr410;
    pValues430[k] = pCosAmp430;
    pErrValues430[k] = pCosAmpErr430;
    pValues450[k] = pCosAmp450;
    pErrValues450[k] = pCosAmpErr450;
    pValues470[k] = pCosAmp470;
    pErrValues470[k] = pCosAmpErr470;
    pValues490[k] = pCosAmp490;
    pErrValues490[k] = pCosAmpErr490;
    pValues510[k] = pCosAmp510;
    pErrValues510[k] = pCosAmpErr510;
    pValues530[k] = pCosAmp530;
    pErrValues530[k] = pCosAmpErr530;
    pValues550[k] = pCosAmp550;
    pErrValues550[k] = pCosAmpErr550;
    pValues570[k] = pCosAmp570;
    pErrValues570[k] = pCosAmpErr570;
    pValues590[k] = pCosAmp590;
    pErrValues590[k] = pCosAmpErr590;
    pValues610[k] = pCosAmp610;
    pErrValues610[k] = pCosAmpErr610;
    pValues630[k] = pCosAmp630;
    pErrValues630[k] = pCosAmpErr630;

  }

  TFile *f2= TFile::Open("/scratch/Mainz_Software/a2GoAT/LinPol_Aug16.root"); // Open linear polarisation plot

  // Calculate values of sigma for each angular and energy bin
  for (Int_t i = 0; i < 10; i++){

    pSigmaValues410[i] = pValues410[i]/(Graph->Eval(410,0));
    pSigmaErrValues410[i] = pErrValues410[i]/(Graph->Eval(410,0));
    pSigmaValues430[i] = pValues430[i]/(Graph->Eval(430,0));
    pSigmaErrValues430[i] = pErrValues430[i]/(Graph->Eval(430,0));
    pSigmaValues450[i] = pValues450[i]/(Graph->Eval(450,0));
    pSigmaErrValues450[i] = pErrValues450[i]/(Graph->Eval(450,0));
    pSigmaValues470[i] = pValues470[i]/(Graph->Eval(470,0));
    pSigmaErrValues470[i] = pErrValues470[i]/(Graph->Eval(470,0));
    pSigmaValues490[i] = pValues490[i]/(Graph->Eval(490,0));
    pSigmaErrValues490[i] = pErrValues490[i]/(Graph->Eval(490,0));
    pSigmaValues510[i] = pValues510[i]/(Graph->Eval(510,0));
    pSigmaErrValues510[i] = pErrValues510[i]/(Graph->Eval(510,0));
    pSigmaValues530[i] = pValues530[i]/(Graph->Eval(530,0));
    pSigmaErrValues530[i] = pErrValues530[i]/(Graph->Eval(530,0));
    pSigmaValues550[i] = pValues550[i]/(Graph->Eval(550,0));
    pSigmaErrValues550[i] = pErrValues550[i]/(Graph->Eval(550,0));
    pSigmaValues570[i] = pValues570[i]/(Graph->Eval(570,0));
    pSigmaErrValues570[i] = pErrValues570[i]/(Graph->Eval(570,0));
    pSigmaValues590[i] = pValues590[i]/(Graph->Eval(590,0));
    pSigmaErrValues590[i] = pErrValues590[i]/(Graph->Eval(590,0));
    pSigmaValues610[i] = pValues610[i]/(Graph->Eval(610,0));
    pSigmaErrValues610[i] = pErrValues610[i]/(Graph->Eval(610,0));
    pSigmaValues630[i] = pValues630[i]/(Graph->Eval(630,0));
    pSigmaErrValues630[i] = pErrValues630[i]/(Graph->Eval(630,0));

  }

  TFile f3("Sigma_Plots_12.root", "RECREATE");

  Float_t xMin = -1;
  Float_t xMax = 1;
  Float_t yMin = -5;
  Float_t yMax = 5;
  Double_t x[10] = {0.9, 0.7, 0.5, 0.3, 0.1, -0.1, -0.3, -0.5, -0.7, -0.9}; // Need to adjust
  Double_t ex[10] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1}; // Need to adjust

  TCanvas *canvas = new TCanvas("canvas","canvas", 1920, 1080);
  TPad *pad1 = new TPad("pad1","",0,0,1,1);
  pad1->Draw();
  pad1->cd();

  pad1->SetTickx(1);
  pad1->SetTicky(1);
  pad1->SetGridx(1);
  pad1->SetGridy(1);
  TH1F  *hr;
  hr = canvas->DrawFrame(xMin,-2.5,xMax,2.5);
  hr->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 400-420MeV)");

  gr1 = new TGraphErrors(10, x, pSigmaValues410, ex, pSigmaErrValues410);
  gr1->SetMarkerColor(2);
  gr1->SetMarkerStyle(5);
  gr1->SetMarkerSize(2);
  gr1->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 400-420MeV)");
  gr1->GetXaxis()->SetTitle("Cos#theta_{CM}");
  gr1->GetYaxis()->SetTitle("#Sigma");
  gr1->Draw("ep");

  TCanvas *canvas1 = new TCanvas("canvas1","canvas1", 1920, 1080);
  TPad *pad2 = new TPad("pad2","",0,0,1,1);
  pad2->Draw();
  pad2->cd();

  pad2->SetTickx(1);
  pad2->SetTicky(1);
  pad2->SetGridx(1);
  pad2->SetGridy(1);
  TH1F  *hr1;
  hr1 = canvas1->DrawFrame(xMin,-2.5,xMax,2.5);
  hr1->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 420-440MeV)");

  gr2 = new TGraphErrors(10, x, pSigmaValues430 , ex, pSigmaErrValues430);
  gr2->SetMarkerColor(2);
  gr2->SetMarkerStyle(5);
  gr2->SetMarkerSize(2);
  gr2->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 420-440MeV)");
  gr2->GetXaxis()->SetTitle("Cos#theta_{CM}");
  gr2->GetYaxis()->SetTitle("#Sigma");
  gr2->Draw("ep");

  TCanvas *canvas2 = new TCanvas("canvas2","canvas2", 1920, 1080);
  TPad *pad3 = new TPad("pad3","",0,0,1,1);
  pad3->Draw();
  pad3->cd();

  pad3->SetTickx(1);
  pad3->SetTicky(1);
  pad3->SetGridx(1);
  pad3->SetGridy(1);
  TH1F  *hr2;
  hr2 = canvas2->DrawFrame(xMin,-2.5,xMax,2.5);
  hr2->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 440-460MeV)");

  gr3 = new TGraphErrors(10, x, pSigmaValues450 , ex, pSigmaErrValues450);
  gr3->SetMarkerColor(2);
  gr3->SetMarkerStyle(5);
  gr3->SetMarkerSize(2);
  gr3->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 440-460MeV)");
  gr3->GetXaxis()->SetTitle("Cos#theta_{CM}");
  gr3->GetYaxis()->SetTitle("#Sigma");
  gr3->Draw("ep");

  TCanvas *canvas3 = new TCanvas("canvas3","canvas3", 1920, 1080);
  TPad *pad4 = new TPad("pad4","",0,0,1,1);
  pad4->Draw();
  pad4->cd();

  pad4->SetTickx(1);
  pad4->SetTicky(1);
  pad4->SetGridx(1);
  pad4->SetGridy(1);
  TH1F  *hr3;
  hr3 = canvas3->DrawFrame(xMin,-2.5,xMax,2.5);
  hr3->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 460-480MeV)");

  gr4 = new TGraphErrors(10, x, pSigmaValues470 , ex, pSigmaErrValues470);
  gr4->SetMarkerColor(2);
  gr4->SetMarkerStyle(5);
  gr4->SetMarkerSize(2);
  gr4->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 460-480MeV)");
  gr4->GetXaxis()->SetTitle("Cos#theta_{CM}");
  gr4->GetYaxis()->SetTitle("#Sigma");
  gr4->Draw("ep");

  TCanvas *canvas4 = new TCanvas("canvas4","canvas4", 1920, 1080);
  TPad *pad5 = new TPad("pad5","",0,0,1,1);
  pad5->Draw();
  pad5->cd();

  pad5->SetTickx(1);
  pad5->SetTicky(1);
  pad5->SetGridx(1);
  pad5->SetGridy(1);
  TH1F  *hr4;
  hr4 = canvas4->DrawFrame(xMin,-2.5,xMax,2.5);
  hr4->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 480-500MeV)");

  gr5 = new TGraphErrors(10, x, pSigmaValues490, ex, pSigmaErrValues490);
  gr5->SetMarkerColor(2);
  gr5->SetMarkerStyle(5);
  gr5->SetMarkerSize(2);
  gr5->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 480-500MeV)");
  gr5->GetXaxis()->SetTitle("Cos#theta_{CM}");
  gr5->GetYaxis()->SetTitle("#Sigma");
  gr5->Draw("ep");

  TCanvas *canvas5 = new TCanvas("canvas5","canvas5", 1920, 1080);
  TPad *pad6 = new TPad("pad6","",0,0,1,1);
  pad6->Draw();
  pad6->cd();

  pad6->SetTickx(1);
  pad6->SetTicky(1);
  pad6->SetGridx(1);
  pad6->SetGridy(1);
  TH1F  *hr5;
  hr5 = canvas5->DrawFrame(xMin,-2.5,xMax,2.5);
  hr5->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 500-520MeV)");

  gr6 = new TGraphErrors(10, x, pSigmaValues510, ex, pSigmaErrValues510);
  gr6->SetMarkerColor(2);
  gr6->SetMarkerStyle(5);
  gr6->SetMarkerSize(2);
  gr6->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 500-520MeV)");
  gr6->GetXaxis()->SetTitle("Cos#theta_{CM}");
  gr6->GetYaxis()->SetTitle("#Sigma");
  gr6->Draw("ep");

  TCanvas *canvas6 = new TCanvas("canvas6","canvas6", 1920, 1080);
  TPad *pad7 = new TPad("pad7","",0,0,1,1);
  pad7->Draw();
  pad7->cd();

  pad7->SetTickx(1);
  pad7->SetTicky(1);
  pad7->SetGridx(1);
  pad7->SetGridy(1);
  TH1F  *hr6;
  hr6 = canvas6->DrawFrame(xMin,-1,xMax,1);
  hr6->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 520-540MeV)");

  gr7 = new TGraphErrors(10, x, pSigmaValues530, ex, pSigmaErrValues530);
  gr7->SetMarkerColor(2);
  gr7->SetMarkerStyle(5);
  gr7->SetMarkerSize(2);
  gr7->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 520-540MeV)");
  gr7->GetXaxis()->SetTitle("Cos#theta_{CM}");
  gr7->GetYaxis()->SetTitle("#Sigma");
  gr7->Draw("ep");

  TCanvas *canvas7 = new TCanvas("canvas7","canvas7", 1920, 1080);
  TPad *pad8 = new TPad("pad8","",0,0,1,1);
  pad8->Draw();
  pad8->cd();

  pad8->SetTickx(1);
  pad8->SetTicky(1);
  pad8->SetGridx(1);
  pad8->SetGridy(1);
  TH1F  *hr7;
  hr7 = canvas7->DrawFrame(xMin,-2.5,xMax,2.5);
  hr7->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 540-560MeV)");

  gr8 = new TGraphErrors(10, x, pSigmaValues550, ex, pSigmaErrValues550);
  gr8->SetMarkerColor(2);
  gr8->SetMarkerStyle(5);
  gr8->SetMarkerSize(2);
  gr8->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 540-560MeV)");
  gr8->GetXaxis()->SetTitle("Cos#theta_{CM}");
  gr8->GetYaxis()->SetTitle("#Sigma");
  gr8->Draw("ep");

  TCanvas *canvas8 = new TCanvas("canvas8","canvas8", 1920, 1080);
  TPad *pad9 = new TPad("pad9","",0,0,1,1);
  pad9->Draw();
  pad9->cd();

  pad9->SetTickx(1);
  pad9->SetTicky(1);
  pad9->SetGridx(1);
  pad9->SetGridy(1);
  TH1F  *hr8;
  hr8 = canvas8->DrawFrame(xMin,-2.5,xMax,2.5);
  hr8->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 560-580MeV)");

  gr9 = new TGraphErrors(10, x, pSigmaValues570, ex, pSigmaErrValues570);
  gr9->SetMarkerColor(2);
  gr9->SetMarkerStyle(5);
  gr9->SetMarkerSize(2);
  gr9->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 560-580MeV)");
  gr9->GetXaxis()->SetTitle("Cos#theta_{CM}");
  gr9->GetYaxis()->SetTitle("#Sigma");
  gr9->Draw("ep");

  TCanvas *canvas9 = new TCanvas("canvas9","canvas9", 1920, 1080);
  TPad *pad10 = new TPad("pad10","",0,0,1,1);
  pad10->Draw();
  pad10->cd();

  pad10->SetTickx(1);
  pad10->SetTicky(1);
  pad10->SetGridx(1);
  pad10->SetGridy(1);
  TH1F  *hr9;
  hr9 = canvas9->DrawFrame(xMin,-1.5,xMax,1.5);
  hr9->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 580-600MeV)");

  gr10 = new TGraphErrors(10, x, pSigmaValues590, ex, pSigmaErrValues590);
  gr10->SetMarkerColor(2);
  gr10->SetMarkerStyle(5);
  gr10->SetMarkerSize(2);
  gr10->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 580-600MeV)");
  gr10->GetXaxis()->SetTitle("Cos#theta_{CM}");
  gr10->GetYaxis()->SetTitle("#Sigma");
  gr10->Draw("ep");

  TCanvas *canvas10 = new TCanvas("canvas10","canvas10", 1920, 1080);
  TPad *pad11 = new TPad("pad11","",0,0,1,1);
  pad11->Draw();
  pad11->cd();

  pad11->SetTickx(1);
  pad11->SetTicky(1);
  pad11->SetGridx(1);
  pad11->SetGridy(1);
  TH1F  *hr10;
  hr10 = canvas10->DrawFrame(xMin,-1.5,xMax,1.5);
  hr10->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 600-620MeV)");

  gr11 = new TGraphErrors(10, x, pSigmaValues610, ex, pSigmaErrValues610);
  gr11->SetMarkerColor(2);
  gr11->SetMarkerStyle(5);
  gr11->SetMarkerSize(2);
  gr11->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 600-620MeV)");
  gr11->GetXaxis()->SetTitle("Cos#theta_{CM}");
  gr11->GetYaxis()->SetTitle("#Sigma");
  gr11->Draw("ep");

  TCanvas *canvas11 = new TCanvas("canvas11","canvas11", 1920, 1080);
  TPad *pad12 = new TPad("pad12","",0,0,1,1);
  pad12->Draw();
  pad12->cd();

  pad12->SetTickx(1);
  pad12->SetTicky(1);
  pad12->SetGridx(1);
  pad12->SetGridy(1);
  TH1F  *hr11;
  hr11 = canvas11->DrawFrame(xMin,-1.5,xMax,1.5);
  hr11->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 620-640MeV)");

  gr12 = new TGraphErrors(10, x, pSigmaValues630, ex, pSigmaErrValues630);
  gr12->SetMarkerColor(2);
  gr12->SetMarkerStyle(5);
  gr12->SetMarkerSize(2);
  gr12->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 620-640MeV)");
  gr12->GetXaxis()->SetTitle("Cos#theta_{CM}");
  gr12->GetYaxis()->SetTitle("#Sigma");
  gr12->Draw("ep");

  //TCanvas *canvas12 = new TCanvas("canvas12","canvas12", 1920, 1080);
  //canvas12->Divide(4,3);
  //canvas12->cd(1);
  //pad1->Draw();
  //canvas12->cd(2);
  //pad2->Draw();
  //canvas12->cd(3);
  //pad3->Draw();
  //canvas12->cd(4);
  //pad4->Draw();
  //canvas12->cd(5);
  //pad5->Draw();
  //canvas12->cd(6);
  //pad6->Draw();
  //canvas12->cd(7);
  //pad7->Draw();
  //canvas12->cd(8);
  //pad8->Draw();
  //canvas12->cd(9);
  //pad9->Draw();
  //canvas12->cd(10);
  //pad10->Draw();
  //canvas12->cd(11);
  //pad11->Draw();
  //canvas12->cd(12);
  //pad12->Draw();

  canvas->Write();
  canvas1->Write();
  canvas2->Write();
  canvas3->Write();
  canvas4->Write();
  canvas5->Write();
  canvas6->Write();
  canvas7->Write();
  canvas8->Write();
  canvas9->Write();
  canvas10->Write();
  canvas11->Write();
  gr1->Write();
  gr2->Write();
  gr3->Write();
  gr4->Write();
  gr5->Write();
  gr6->Write();
  gr7->Write();
  gr8->Write();
  gr9->Write();
  gr10->Write();
  gr11->Write();
  gr12->Write();
  //canvas12->Write();
  f3.Write();

}
