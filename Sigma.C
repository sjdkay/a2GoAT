#include "./includes_Sigma.h"

void Sigma(){

  TFile *f1= TFile::Open("/scratch/Mainz_Software/a2GoAT/ParaPerpAsymm_Total_10.root");
  TTree *t1 = (TTree*)f1->Get("Parameter_Values");

  Double_t pValues410[10], pValues430[10], pValues450[10], pValues470[10], pValues490[10], pValues510[10], pValues530[10], pValues550[10], pValues570[10], pValues590[10], pValues610[10], pValues630[10];
  Double_t nValues410[10], nValues430[10], nValues450[10], nValues470[10], nValues490[10], nValues510[10], nValues530[10], nValues550[10], nValues570[10], nValues590[10], nValues610[10], nValues630[10];
  Double_t nErrValues410[10], nErrValues430[10], nErrValues450[10], nErrValues470[10], nErrValues490[10], nErrValues510[10], nErrValues530[10], nErrValues550[10], nErrValues570[10], nErrValues590[10], nErrValues610[10], nErrValues630[10];
  Double_t nErrValues410[10], nErrValues430[10], nErrValues450[10], nErrValues470[10], nErrValues490[10], nErrValues510[10], nErrValues530[10], nErrValues550[10], nErrValues570[10], nErrValues590[10], nErrValues610[10], nErrValues630[10];
  Double_t pSigmaValues410[10], pSigmaValues430[10], pSigmaValues450[10], pSigmaValues470[10], pSigmaValues490[10], pSigmaValues510[10], pSigmaValues530[10], pSigmaValues550[10], pSigmaValues570[10], pSigmaValues590[10], pSigmaValues610[10], pSigmaValues630[10];
  Double_t nSigmaValues410[10], nSigmaValues430[10], nSigmaValues450[10], nSigmaValues470[10], nSigmaValues490[10], nSigmaValues510[10], nSigmaValues530[10], nSigmaValues550[10], nSigmaValues570[10], nSigmaValues590[10], nSigmaValues610[10], nSigmaValues630[10];
  Double_t pSigmaErrValues410[10], pSigmaErrValues430[10], pSigmaErrValues450[10], pSigmaErrValues470[10], pSigmaErrValues490[10], pSigmaErrValues510[10], pSigmaErrValues530[10], pSigmaErrValues550[10], pSigmaErrValues570[10], pSigmaErrValues590[10], pSigmaErrValues610[10], pSigmaErrValues630[10];
  Double_t nSigmaErrValues410[10], nSigmaErrValues430[10], nSigmaErrValues450[10], nSigmaErrValues470[10], nSigmaErrValues490[10], nSigmaErrValues510[10], nSigmaErrValues530[10], nSigmaErrValues550[10], nSigmaErrValues570[10], nSigmaErrValues590[10], nSigmaErrValues610[10], nSigmaErrValues630[10];

  Double_t pCosAmp410, pCosAmpErr410, nCosAmp410, nCosAmpErr410;
  Double_t pCosAmp430, pCosAmpErr430, nCosAmp430, nCosAmpErr430;
  Double_t pCosAmp450, pCosAmpErr450, nCosAmp450, nCosAmpErr450;
  Double_t pCosAmp470, pCosAmpErr470, nCosAmp470, nCosAmpErr470;
  Double_t pCosAmp490, pCosAmpErr490, nCosAmp490, nCosAmpErr490;
  Double_t pCosAmp510, pCosAmpErr510, nCosAmp510, nCosAmpErr510;
  Double_t pCosAmp530, pCosAmpErr530, nCosAmp530, nCosAmpErr530;
  Double_t pCosAmp550, pCosAmpErr550, nCosAmp550, nCosAmpErr550;
  Double_t pCosAmp570, pCosAmpErr570, nCosAmp570, nCosAmpErr570;
  Double_t pCosAmp590, pCosAmpErr590, nCosAmp590, nCosAmpErr590;
  Double_t pCosAmp610, pCosAmpErr610, nCosAmp610, nCosAmpErr610;
  Double_t pCosAmp630, pCosAmpErr630, nCosAmp630, nCosAmpErr630;

  // Set branch addresses to get values from
  t1->SetBranchAddress("pCosAmp410", &pCosAmp410);
  t1->SetBranchAddress("pCosAmpErr410", &pCosAmpErr410);
  t1->SetBranchAddress("nCosAmp410", &nCosAmp410);
  t1->SetBranchAddress("nCosAmpErr410", &nCosAmpErr410);
  t1->SetBranchAddress("pCosAmp430", &pCosAmp430);
  t1->SetBranchAddress("pCosAmpErr430", &pCosAmpErr430);
  t1->SetBranchAddress("nCosAmp430", &nCosAmp430);
  t1->SetBranchAddress("nCosAmpErr430", &nCosAmpErr430);
  t1->SetBranchAddress("pCosAmp450", &pCosAmp450);
  t1->SetBranchAddress("pCosAmpErr450", &pCosAmpErr450);
  t1->SetBranchAddress("nCosAmp450", &nCosAmp450);
  t1->SetBranchAddress("nCosAmpErr450", &nCosAmpErr450);
  t1->SetBranchAddress("pCosAmp470", &pCosAmp470);
  t1->SetBranchAddress("pCosAmpErr470", &pCosAmpErr470);
  t1->SetBranchAddress("nCosAmp470", &nCosAmp470);
  t1->SetBranchAddress("nCosAmpErr470", &nCosAmpErr470);
  t1->SetBranchAddress("pCosAmp490", &pCosAmp490);
  t1->SetBranchAddress("pCosAmpErr490", &pCosAmpErr490);
  t1->SetBranchAddress("nCosAmp490", &nCosAmp490);
  t1->SetBranchAddress("nCosAmpErr490", &nCosAmpErr490);
  t1->SetBranchAddress("pCosAmp510", &pCosAmp510);
  t1->SetBranchAddress("pCosAmpErr510", &pCosAmpErr510);
  t1->SetBranchAddress("nCosAmp510", &nCosAmp510);
  t1->SetBranchAddress("nCosAmpErr510", &nCosAmpErr510);
  t1->SetBranchAddress("pCosAmp530", &pCosAmp530);
  t1->SetBranchAddress("pCosAmpErr530", &pCosAmpErr530);
  t1->SetBranchAddress("nCosAmp530", &nCosAmp530);
  t1->SetBranchAddress("nCosAmpErr530", &nCosAmpErr530);
  t1->SetBranchAddress("pCosAmp550", &pCosAmp550);
  t1->SetBranchAddress("pCosAmpErr550", &pCosAmpErr550);
  t1->SetBranchAddress("nCosAmp550", &nCosAmp550);
  t1->SetBranchAddress("nCosAmpErr550", &nCosAmpErr550);
  t1->SetBranchAddress("pCosAmp570", &pCosAmp570);
  t1->SetBranchAddress("pCosAmpErr570", &pCosAmpErr570);
  t1->SetBranchAddress("nCosAmp570", &nCosAmp570);
  t1->SetBranchAddress("nCosAmpErr570", &nCosAmpErr570);
  t1->SetBranchAddress("pCosAmp590", &pCosAmp590);
  t1->SetBranchAddress("pCosAmpErr590", &pCosAmpErr590);
  t1->SetBranchAddress("nCosAmp590", &nCosAmp590);
  t1->SetBranchAddress("nCosAmpErr590", &nCosAmpErr590);
  t1->SetBranchAddress("pCosAmp610", &pCosAmp610);
  t1->SetBranchAddress("pCosAmpErr610", &pCosAmpErr610);
  t1->SetBranchAddress("nCosAmp610", &nCosAmp610);
  t1->SetBranchAddress("nCosAmpErr610", &nCosAmpErr610);
  t1->SetBranchAddress("pCosAmp630", &pCosAmp630);
  t1->SetBranchAddress("pCosAmpErr630", &pCosAmpErr630);
  t1->SetBranchAddress("nCosAmp630", &nCosAmp630);
  t1->SetBranchAddress("nCosAmpErr630", &nCosAmpErr630);

  // Load values from tree and asign values back into an array
  for (Int_t k = 0; k < 10; k++){
    Parameter_Values->GetEntry(k);
    pValues410[k] = pCosAmp410;
    pErrValues410[k] = pCosAmpErr410;
    nValues410[k] = nCosAmp410;
    nErrValues410[k] = nCosAmpErr410;
    pValues430[k] = pCosAmp430;
    pErrValues430[k] = pCosAmpErr430;
    nValues430[k] = nCosAmp430;
    nErrValues430[k] = nCosAmpErr430;
    pValues450[k] = pCosAmp450;
    pErrValues450[k] = pCosAmpErr450;
    nValues450[k] = nCosAmp450;
    nErrValues450[k] = nCosAmpErr450;
    pValues470[k] = pCosAmp470;
    pErrValues470[k] = pCosAmpErr470;
    nValues470[k] = nCosAmp470;
    nErrValues470[k] = nCosAmpErr470;
    pValues490[k] = pCosAmp490;
    pErrValues490[k] = pCosAmpErr490;
    nValues490[k] = nCosAmp490;
    nErrValues490[k] = nCosAmpErr490;
    pValues510[k] = pCosAmp510;
    pErrValues510[k] = pCosAmpErr510;
    nValues510[k] = nCosAmp510;
    nErrValues510[k] = nCosAmpErr510;
    pValues530[k] = pCosAmp530;
    pErrValues530[k] = pCosAmpErr530;
    nValues530[k] = nCosAmp530;
    nErrValues530[k] = nCosAmpErr530;
    pValues550[k] = pCosAmp550;
    pErrValues550[k] = pCosAmpErr550;
    nValues550[k] = nCosAmp550;
    nErrValues550[k] = nCosAmpErr550;
    pValues570[k] = pCosAmp570;
    pErrValues570[k] = pCosAmpErr570;
    nValues570[k] = nCosAmp570;
    nErrValues570[k] = nCosAmpErr570;
    pValues590[k] = pCosAmp590;
    pErrValues590[k] = pCosAmpErr590;
    nValues590[k] = nCosAmp590;
    nErrValues590[k] = nCosAmpErr590;
    pValues610[k] = pCosAmp610;
    pErrValues610[k] = pCosAmpErr610;
    nValues610[k] = nCosAmp610;
    nErrValues610[k] = nCosAmpErr610;
    pValues630[k] = pCosAmp630;
    pErrValues630[k] = pCosAmpErr630;
    nValues630[k] = nCosAmp630;
    nErrValues630[k] = nCosAmpErr630;

  }

  TFile *f2= TFile::Open("/scratch/Mainz_Software/a2GoAT/LinPol_Aug16.root"); // Open linear polarisation plot

  // Calculate values of sigma for each angular and energy bin
  for (Int_t i = 0; i < 10; i++){

    pSigmaValues410[i] = pValues410[i]/(Graph->Eval(410,0));
    pSigmaErrValues410[i] = pErrValues410[i]/(Graph->Eval(410,0)); // Double check error formula for this as new way around now
    pSigmaValues430[i] = pValues430[i]/(Graph->Eval(430,0));
    pSigmaErrValues430[i] = pErrValues430[i]/(Graph->Eval(430,0)); // Double check error formula for this as new way around now
    pSigmaValues450[i] = pValues450[i]/(Graph->Eval(450,0));
    pSigmaErrValues450[i] = pErrValues450[i]/(Graph->Eval(450,0)); // Double check error formula for this as new way around now
    pSigmaValues470[i] = pValues470[i]/(Graph->Eval(470,0));
    pSigmaErrValues470[i] = pErrValues470[i]/(Graph->Eval(470,0)); // Double check error formula for this as new way around now
    pSigmaValues490[i] = pValues490[i]/(Graph->Eval(490,0));
    pSigmaErrValues490[i] = pErrValues490[i]/(Graph->Eval(490,0)); // Double check error formula for this as new way around now
    pSigmaValues510[i] = pValues510[i]/(Graph->Eval(510,0));
    pSigmaErrValues510[i] = pErrValues510[i]/(Graph->Eval(510,0)); // Double check error formula for this as new way around now
    pSigmaValues530[i] = pValues530[i]/(Graph->Eval(530,0));
    pSigmaErrValues530[i] = pErrValues530[i]/(Graph->Eval(530,0)); // Double check error formula for this as new way around now
    pSigmaValues550[i] = pValues550[i]/(Graph->Eval(550,0));
    pSigmaErrValues550[i] = pErrValues550[i]/(Graph->Eval(550,0)); // Double check error formula for this as new way around now
    pSigmaValues570[i] = pValues570[i]/(Graph->Eval(570,0));
    pSigmaErrValues570[i] = pErrValues570[i]/(Graph->Eval(570,0)); // Double check error formula for this as new way around now
    pSigmaValues590[i] = pValues590[i]/(Graph->Eval(590,0));
    pSigmaErrValues590[i] = pErrValues590[i]/(Graph->Eval(590,0)); // Double check error formula for this as new way around now
    pSigmaValues610[i] = pValues610[i]/(Graph->Eval(610,0));
    pSigmaErrValues610[i] = pErrValues610[i]/(Graph->Eval(610,0)); // Double check error formula for this as new way around now
    pSigmaValues630[i] = pValues630[i]/(Graph->Eval(630,0));
    pSigmaErrValues630[i] = pErrValues630[i]/(Graph->Eval(630,0)); // Double check error formula for this as new way around now

    nSigmaValues410[i] = nValues410[i]/(Graph->Eval(410,0));
    nSigmaErrValues410[i] = nErrValues410[i]/(Graph->Eval(410,0)); // Double check error formula for this as new way around now
    nSigmaValues430[i] = nValues430[i]/(Graph->Eval(430,0));
    nSigmaErrValues430[i] = nErrValues430[i]/(Graph->Eval(430,0)); // Double check error formula for this as new way around now
    nSigmaValues450[i] = nValues450[i]/(Graph->Eval(450,0));
    nSigmaErrValues450[i] = nErrValues450[i]/(Graph->Eval(450,0)); // Double check error formula for this as new way around now
    nSigmaValues470[i] = nValues470[i]/(Graph->Eval(470,0));
    nSigmaErrValues470[i] = nErrValues470[i]/(Graph->Eval(470,0)); // Double check error formula for this as new way around now
    nSigmaValues490[i] = nValues490[i]/(Graph->Eval(490,0));
    nSigmaErrValues490[i] = nErrValues490[i]/(Graph->Eval(490,0)); // Double check error formula for this as new way around now
    nSigmaValues510[i] = nValues510[i]/(Graph->Eval(510,0));
    nSigmaErrValues510[i] = nErrValues510[i]/(Graph->Eval(510,0)); // Double check error formula for this as new way around now
    nSigmaValues530[i] = nValues530[i]/(Graph->Eval(530,0));
    nSigmaErrValues530[i] = nErrValues530[i]/(Graph->Eval(530,0)); // Double check error formula for this as new way around now
    nSigmaValues550[i] = nValues550[i]/(Graph->Eval(550,0));
    nSigmaErrValues550[i] = nErrValues550[i]/(Graph->Eval(550,0)); // Double check error formula for this as new way around now
    nSigmaValues570[i] = nValues570[i]/(Graph->Eval(570,0));
    nSigmaErrValues570[i] = nErrValues570[i]/(Graph->Eval(570,0)); // Double check error formula for this as new way around now
    nSigmaValues590[i] = nValues590[i]/(Graph->Eval(590,0));
    nSigmaErrValues590[i] = nErrValues590[i]/(Graph->Eval(590,0)); // Double check error formula for this as new way around now
    nSigmaValues610[i] = nValues610[i]/(Graph->Eval(610,0));
    nSigmaErrValues610[i] = nErrValues610[i]/(Graph->Eval(610,0)); // Double check error formula for this as new way around now
    nSigmaValues630[i] = nValues630[i]/(Graph->Eval(630,0));
    nSigmaErrValues630[i] = nErrValues630[i]/(Graph->Eval(630,0)); // Double check error formula for this as new way around now

  }

  TFile f3("Sigma_Plots.root", "RECREATE");

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
  hr = canvas->DrawFrame(xMin,yMin,xMax,yMax);
  hr->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 400-420MeV)");

  gr1 = new TGraphErrors(10, x, pSigmaValues410, ex, pSigmaErrValues410);
  gr1->SetMarkerColor(2);
  gr1->SetMarkerStyle(21);
  gr1->SetMarkerSize(2);
  gr2 = new TGraphErrors(10, x, nSigmaValues410, ex, nSigmaErrValues410);
  gr2->SetMarkerColor(4);
  gr2->SetMarkerStyle(22);
  gr2->SetMarkerSize(2);
  TMultiGraph *mg1 = new TMultiGraph();
  mg1->Add(gr1, "ep");
  mg1->Add(gr2, "ep");
  mg1->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 400-420MeV)");
  mg1->GetXaxis()->SetTitle("Cos#theta_{CM}");
  mg1->GetYaxis()->SetTitle("#Sigma");
  leg = new TLegend(0.75, 0.75, 0.95, 0.95); // Add legend to plot
  leg->AddEntry(gr1, "Proton", "ep");
  leg->AddEntry(gr2, "Neutron", "ep");
  mg1->Draw();
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
  hr1 = canvas1->DrawFrame(xMin,yMin/2,xMax,yMax/2);
  hr1->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 420-440MeV)");

  gr3 = new TGraphErrors(10, x, pSigmaValues430 , ex, pSigmaErrValues430);
  gr3->SetMarkerColor(2);
  gr3->SetMarkerStyle(21);
  gr3->SetMarkerSize(2);
  gr4 = new TGraphErrors(10, x, nSigmaValues430, ex, nSigmaErrValues430);
  gr4->SetMarkerColor(4);
  gr4->SetMarkerStyle(22);
  gr4->SetMarkerSize(2);
  TMultiGraph *mg2 = new TMultiGraph();
  mg2->Add(gr3, "ep");
  mg2->Add(gr4, "ep");
  mg2->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 420-440MeV)");
  mg2->GetXaxis()->SetTitle("Cos#theta_{CM}");
  mg2->GetYaxis()->SetTitle("#Sigma");
  leg2 = new TLegend(0.75, 0.75, 0.95, 0.95); // Add legend to plot
  leg2->AddEntry(gr3, "Proton", "ep");
  leg2->AddEntry(gr4, "Neutron", "ep");
  mg2->Draw();
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
  hr2 = canvas2->DrawFrame(xMin,yMin/2,xMax,yMax/2);
  hr2->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 440-460MeV)");

  gr5 = new TGraphErrors(10, x, pSigmaValues450 , ex, pSigmaErrValues450);
  gr5->SetMarkerColor(2);
  gr5->SetMarkerStyle(21);
  gr5->SetMarkerSize(2);
  gr6 = new TGraphErrors(10, x, nSigmaValues450, ex, nSigmaErrValues450);
  gr6->SetMarkerColor(4);
  gr6->SetMarkerStyle(22);
  gr6->SetMarkerSize(2);
  TMultiGraph *mg3 = new TMultiGraph();
  mg3->Add(gr5, "ep");
  mg3->Add(gr6, "ep");
  mg3->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 440-460MeV)");
  mg3->GetXaxis()->SetTitle("Cos#theta_{CM}");
  mg3->GetYaxis()->SetTitle("#Sigma");
  leg3 = new TLegend(0.75, 0.75, 0.95, 0.95); // Add legend to plot
  leg3->AddEntry(gr5, "Proton", "ep");
  leg3->AddEntry(gr6, "Neutron", "ep");
  mg3->Draw();
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
  hr3 = canvas3->DrawFrame(xMin,yMin/2,xMax,yMax/2);
  hr3->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 460-480MeV)");

  gr7 = new TGraphErrors(10, x, pSigmaValues470 , ex, pSigmaErrValues470);
  gr7->SetMarkerColor(2);
  gr7->SetMarkerStyle(21);
  gr7->SetMarkerSize(2);
  gr8 = new TGraphErrors(10, x, nSigmaValues470, ex, nSigmaErrValues470);
  gr8->SetMarkerColor(4);
  gr8->SetMarkerStyle(22);
  gr8->SetMarkerSize(2);
  TMultiGraph *mg4 = new TMultiGraph();
  mg4->Add(gr7, "ep");
  mg4->Add(gr8, "ep");
  mg4->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 460-480MeV)");
  mg4->GetXaxis()->SetTitle("Cos#theta_{CM}");
  mg4->GetYaxis()->SetTitle("#Sigma");
  leg4 = new TLegend(0.75, 0.75, 0.95, 0.95); // Add legend to plot
  leg4->AddEntry(gr7, "Proton", "ep");
  leg4->AddEntry(gr8, "Neutron", "ep");
  mg4->Draw();
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
  hr4 = canvas4->DrawFrame(xMin,yMin/2,xMax,yMax/2);
  hr4->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 480-500MeV)");

  gr9 = new TGraphErrors(10, x, pSigmaValues490, ex, pSigmaErrValues490);
  gr9->SetMarkerColor(2);
  gr9->SetMarkerStyle(21);
  gr9->SetMarkerSize(2);
  gr10 = new TGraphErrors(10, x, nSigmaValues490, ex, nSigmaErrValues490);
  gr10->SetMarkerColor(4);
  gr10->SetMarkerStyle(22);
  gr10->SetMarkerSize(2);
  TMultiGraph *mg5 = new TMultiGraph();
  mg5->Add(gr9, "ep");
  mg5->Add(gr10, "ep");
  mg5->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 480-500MeV)");
  mg5->GetXaxis()->SetTitle("Cos#theta_{CM}");
  mg5->GetYaxis()->SetTitle("#Sigma");
  leg5 = new TLegend(0.75, 0.75, 0.95, 0.95); // Add legend to plot
  leg5->AddEntry(gr9, "Proton", "ep");
  leg5->AddEntry(gr10, "Neutron", "ep");
  mg5->Draw();
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
  hr5->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 500-520MeV)");

  gr11 = new TGraphErrors(10, x, pSigmaValues510, ex, pSigmaErrValues510);
  gr11->SetMarkerColor(2);
  gr11->SetMarkerStyle(21);
  gr11->SetMarkerSize(2);
  gr12 = new TGraphErrors(10, x, nSigmaValues510, ex, nSigmaErrValues510);
  gr12->SetMarkerColor(4);
  gr12->SetMarkerStyle(22);
  gr12->SetMarkerSize(2);
  TMultiGraph *mg6 = new TMultiGraph();
  mg6->Add(gr11, "ep");
  mg6->Add(gr12, "ep");
  mg6->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 500-520MeV)");
  mg6->GetXaxis()->SetTitle("Cos#theta_{CM}");
  mg6->GetYaxis()->SetTitle("#Sigma");
  leg6 = new TLegend(0.75, 0.75, 0.95, 0.95); // Add legend to plot
  leg6->AddEntry(gr11, "Proton", "ep");
  leg6->AddEntry(gr12, "Neutron", "ep");
  mg6->Draw();
  leg6->Draw("Same");

  TCanvas *canvas6 = new TCanvas("canvas6","canvas6", 1920, 1080);
  TPad *pad7 = new TPad("pad7","",0,0,1,1);
  pad7->Draw();
  pad7->cd();

  pad7->SetTickx(1);
  pad7->SetTicky(1);
  pad7->SetGridx(1);
  pad7->SetGridy(1);
  TH1F  *hr6;
  hr6 = canvas6->DrawFrame(xMin,yMin,xMax,yMax);
  hr6->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 520-540MeV)");

  gr13 = new TGraphErrors(10, x, pSigmaValues510, ex, pSigmaErrValues510);
  gr13->SetMarkerColor(2);
  gr13->SetMarkerStyle(21);
  gr13->SetMarkerSize(2);
  gr14 = new TGraphErrors(10, x, nSigmaValues510, ex, nSigmaErrValues510);
  gr14->SetMarkerColor(4);
  gr14->SetMarkerStyle(22);
  gr14->SetMarkerSize(2);
  TMultiGraph *mg7 = new TMultiGraph();
  mg7->Add(gr13, "ep");
  mg7->Add(gr14, "ep");
  mg7->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 520-540MeV)");
  mg7->GetXaxis()->SetTitle("Cos#theta_{CM}");
  mg7->GetYaxis()->SetTitle("#Sigma");
  leg7 = new TLegend(0.75, 0.75, 0.95, 0.95); // Add legend to plot
  leg7->AddEntry(gr13, "Proton", "ep");
  leg7->AddEntry(gr14, "Neutron", "ep");
  mg7->Draw();
  leg7->Draw("Same");

  TCanvas *canvas7 = new TCanvas("canvas7","canvas7", 1920, 1080);
  TPad *pad8 = new TPad("pad8","",0,0,1,1);
  pad8->Draw();
  pad8->cd();

  pad8->SetTickx(1);
  pad8->SetTicky(1);
  pad8->SetGridx(1);
  pad8->SetGridy(1);
  TH1F  *hr7;
  hr7 = canvas7->DrawFrame(xMin,yMin,xMax,yMax);
  hr7->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 540-560MeV)");

  gr15 = new TGraphErrors(10, x, pSigmaValues510, ex, pSigmaErrValues510);
  gr15->SetMarkerColor(2);
  gr15->SetMarkerStyle(21);
  gr15->SetMarkerSize(2);
  gr16 = new TGraphErrors(10, x, nSigmaValues510, ex, nSigmaErrValues510);
  gr16->SetMarkerColor(4);
  gr16->SetMarkerStyle(22);
  gr16->SetMarkerSize(2);
  TMultiGraph *mg8 = new TMultiGraph();
  mg8->Add(gr15, "ep");
  mg8->Add(gr16, "ep");
  mg8->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 540-560MeV)");
  mg8->GetXaxis()->SetTitle("Cos#theta_{CM}");
  mg8->GetYaxis()->SetTitle("#Sigma");
  leg8 = new TLegend(0.75, 0.75, 0.95, 0.95); // Add legend to plot
  leg8->AddEntry(gr15, "Proton", "ep");
  leg8->AddEntry(gr16, "Neutron", "ep");
  mg8->Draw();
  leg8->Draw("Same");

  TCanvas *canvas8 = new TCanvas("canvas8","canvas8", 1920, 1080);
  TPad *pad9 = new TPad("pad9","",0,0,1,1);
  pad9->Draw();
  pad9->cd();

  pad9->SetTickx(1);
  pad9->SetTicky(1);
  pad9->SetGridx(1);
  pad9->SetGridy(1);
  TH1F  *hr8;
  hr8 = canvas7->DrawFrame(xMin,yMin,xMax,yMax);
  hr8->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 560-580MeV)");

  gr17 = new TGraphErrors(10, x, pSigmaValues510, ex, pSigmaErrValues510);
  gr17->SetMarkerColor(2);
  gr17->SetMarkerStyle(21);
  gr17->SetMarkerSize(2);
  gr18 = new TGraphErrors(10, x, nSigmaValues510, ex, nSigmaErrValues510);
  gr18->SetMarkerColor(4);
  gr18->SetMarkerStyle(22);
  gr18->SetMarkerSize(2);
  TMultiGraph *mg9 = new TMultiGraph();
  mg9->Add(gr17, "ep");
  mg9->Add(gr18, "ep");
  mg9->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 560-580MeV)");
  mg9->GetXaxis()->SetTitle("Cos#theta_{CM}");
  mg9->GetYaxis()->SetTitle("#Sigma");
  leg9 = new TLegend(0.75, 0.75, 0.95, 0.95); // Add legend to plot
  leg9->AddEntry(gr17, "Proton", "ep");
  leg9->AddEntry(gr18, "Neutron", "ep");
  mg9->Draw();
  leg9->Draw("Same");

  TCanvas *canvas9 = new TCanvas("canvas9","canvas9", 1920, 1080);
  TPad *pad10 = new TPad("pad10","",0,0,1,1);
  pad10->Draw();
  pad10->cd();

  pad10->SetTickx(1);
  pad10->SetTicky(1);
  pad10->SetGridx(1);
  pad10->SetGridy(1);
  TH1F  *hr9;
  hr9 = canvas7->DrawFrame(xMin,yMin,xMax,yMax);
  hr9->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 580-600MeV)");

  gr19 = new TGraphErrors(10, x, pSigmaValues510, ex, pSigmaErrValues510);
  gr19->SetMarkerColor(2);
  gr19->SetMarkerStyle(21);
  gr19->SetMarkerSize(2);
  gr20 = new TGraphErrors(10, x, nSigmaValues510, ex, nSigmaErrValues510);
  gr20->SetMarkerColor(4);
  gr20->SetMarkerStyle(22);
  gr20->SetMarkerSize(2);
  TMultiGraph *mg10 = new TMultiGraph();
  mg10->Add(gr19, "ep");
  mg10->Add(gr20, "ep");
  mg10->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 580-600MeV)");
  mg10->GetXaxis()->SetTitle("Cos#theta_{CM}");
  mg10->GetYaxis()->SetTitle("#Sigma");
  leg10 = new TLegend(0.75, 0.75, 0.95, 0.95); // Add legend to plot
  leg10->AddEntry(gr19, "Proton", "ep");
  leg10->AddEntry(gr20, "Neutron", "ep");
  mg10->Draw();
  leg10->Draw("Same");

  TCanvas *canvas10 = new TCanvas("canvas10","canvas10", 1920, 1080);
  TPad *pad11 = new TPad("pad11","",0,0,1,1);
  pad11->Draw();
  pad11->cd();

  pad11->SetTickx(1);
  pad11->SetTicky(1);
  pad11->SetGridx(1);
  pad11->SetGridy(1);
  TH1F  *hr10;
  hr10 = canvas7->DrawFrame(xMin,yMin,xMax,yMax);
  hr10->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 600-620MeV)");

  gr21 = new TGraphErrors(10, x, pSigmaValues510, ex, pSigmaErrValues510);
  gr21->SetMarkerColor(2);
  gr21->SetMarkerStyle(21);
  gr21->SetMarkerSize(2);
  gr22 = new TGraphErrors(10, x, nSigmaValues510, ex, nSigmaErrValues510);
  gr22->SetMarkerColor(4);
  gr22->SetMarkerStyle(22);
  gr22->SetMarkerSize(2);
  TMultiGraph *mg11 = new TMultiGraph();
  mg11->Add(gr21, "ep");
  mg11->Add(gr22, "ep");
  mg11->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 600-620MeV)");
  mg11->GetXaxis()->SetTitle("Cos#theta_{CM}");
  mg11->GetYaxis()->SetTitle("#Sigma");
  leg11 = new TLegend(0.75, 0.75, 0.95, 0.95); // Add legend to plot
  leg11->AddEntry(gr21, "Proton", "ep");
  leg11->AddEntry(gr22, "Neutron", "ep");
  mg11->Draw();
  leg11->Draw("Same");

  TCanvas *canvas11 = new TCanvas("canvas11","canvas11", 1920, 1080);
  TPad *pad12 = new TPad("pad12","",0,0,1,1);
  pad12->Draw();
  pad12->cd();

  pad12->SetTickx(1);
  pad12->SetTicky(1);
  pad12->SetGridx(1);
  pad12->SetGridy(1);
  TH1F  *hr11;
  hr11 = canvas7->DrawFrame(xMin,yMin,xMax,yMax);
  hr11->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 620-640MeV)");

  gr23 = new TGraphErrors(10, x, pSigmaValues510, ex, pSigmaErrValues510);
  gr23->SetMarkerColor(2);
  gr23->SetMarkerStyle(21);
  gr23->SetMarkerSize(2);
  gr24 = new TGraphErrors(10, x, nSigmaValues510, ex, nSigmaErrValues510);
  gr24->SetMarkerColor(4);
  gr24->SetMarkerStyle(22);
  gr24->SetMarkerSize(2);
  TMultiGraph *mg12 = new TMultiGraph();
  mg12->Add(gr23, "ep");
  mg12->Add(gr24, "ep");
  mg12->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 620-640MeV)");
  mg12->GetXaxis()->SetTitle("Cos#theta_{CM}");
  mg12->GetYaxis()->SetTitle("#Sigma");
  leg12 = new TLegend(0.75, 0.75, 0.95, 0.95); // Add legend to plot
  leg12->AddEntry(gr21, "Proton", "ep");
  leg12->AddEntry(gr22, "Neutron", "ep");
  mg12->Draw();
  leg12->Draw("Same");

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
  //canvas12->Write();
  f3.Write();

}
