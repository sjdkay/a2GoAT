#include "./includes_Sigma_NoScatt.h"

void Sigma_NoScatt(){

  TFile *f1= TFile::Open("/scratch/Mainz_Software/a2GoAT/ParaPerpAsymm_NoScatt_Total_1.root");
  TTree *t1 = (TTree*)f1->Get("Parameter_Values");

  Double_t pValues425[10], pValues435[10], pValues445[10], pValues455[10], pValues465[10], pValues475[10], pValues485[10], pValues495[10], pValues505[10], pValues515[10], pValues525[10], pValues535[10], pValues545[10], pValues555[10], pValues565[10], pValues575[10], pValues585[10], pValues595[10], pValues605[10], pValues615[10];
  Double_t pErrValues425[10], pErrValues435[10], pErrValues445[10], pErrValues455[10], pErrValues465[10], pErrValues475[10], pErrValues485[10], pErrValues495[10], pErrValues505[10], pErrValues515[10], pErrValues525[10], pErrValues535[10], pErrValues545[10], pErrValues555[10], pErrValues565[10], pErrValues575[10], pErrValues585[10], pErrValues595[10], pErrValues605[10], pErrValues615[10];
  Double_t pSigmaValues425[10], pSigmaValues435[10], pSigmaValues445[10], pSigmaValues455[10], pSigmaValues465[10], pSigmaValues475[10], pSigmaValues485[10], pSigmaValues495[10], pSigmaValues505[10], pSigmaValues515[10], pSigmaValues525[10], pSigmaValues535[10], pSigmaValues545[10], pSigmaValues555[10], pSigmaValues565[10], pSigmaValues575[10], pSigmaValues585[10], pSigmaValues595[10], pSigmaValues605[10], pSigmaValues615[10];
  Double_t pSigmaErrValues425[10], pSigmaErrValues435[10], pSigmaErrValues445[10], pSigmaErrValues455[10], pSigmaErrValues465[10], pSigmaErrValues475[10], pSigmaErrValues485[10], pSigmaErrValues495[10], pSigmaErrValues505[10], pSigmaErrValues515[10], pSigmaErrValues525[10], pSigmaErrValues535[10], pSigmaErrValues545[10], pSigmaErrValues555[10], pSigmaErrValues565[10], pSigmaErrValues575[10], pSigmaErrValues585[10], pSigmaErrValues595[10], pSigmaErrValues605[10], pSigmaErrValues615[10];;

  Double_t pCosAmp425, pCosAmpErr425;
  Double_t pCosAmp435, pCosAmpErr435;
  Double_t pCosAmp445, pCosAmpErr445;
  Double_t pCosAmp455, pCosAmpErr455;
  Double_t pCosAmp465, pCosAmpErr465;
  Double_t pCosAmp475, pCosAmpErr475;
  Double_t pCosAmp485, pCosAmpErr485;
  Double_t pCosAmp495, pCosAmpErr495;
  Double_t pCosAmp505, pCosAmpErr505;
  Double_t pCosAmp515, pCosAmpErr515;
  Double_t pCosAmp525, pCosAmpErr525;
  Double_t pCosAmp535, pCosAmpErr535;
  Double_t pCosAmp545, pCosAmpErr545;
  Double_t pCosAmp555, pCosAmpErr555;
  Double_t pCosAmp565, pCosAmpErr565;
  Double_t pCosAmp575, pCosAmpErr575;
  Double_t pCosAmp585, pCosAmpErr585;
  Double_t pCosAmp595, pCosAmpErr595;
  Double_t pCosAmp605, pCosAmpErr605;
  Double_t pCosAmp615, pCosAmpErr615;

  // Set branch addresses to get values from
  t1->SetBranchAddress("pCosAmp425", &pCosAmp425);
  t1->SetBranchAddress("pCosAmpErr425", &pCosAmpErr425);
  t1->SetBranchAddress("pCosAmp435", &pCosAmp435);
  t1->SetBranchAddress("pCosAmpErr435", &pCosAmpErr435);
  t1->SetBranchAddress("pCosAmp445", &pCosAmp445);
  t1->SetBranchAddress("pCosAmpErr445", &pCosAmpErr445);
  t1->SetBranchAddress("pCosAmp455", &pCosAmp455);
  t1->SetBranchAddress("pCosAmpErr455", &pCosAmpErr455);
  t1->SetBranchAddress("pCosAmp465", &pCosAmp465);
  t1->SetBranchAddress("pCosAmpErr465", &pCosAmpErr465);
  t1->SetBranchAddress("pCosAmp475", &pCosAmp475);
  t1->SetBranchAddress("pCosAmpErr475", &pCosAmpErr475);
  t1->SetBranchAddress("pCosAmp485", &pCosAmp485);
  t1->SetBranchAddress("pCosAmpErr485", &pCosAmpErr485);
  t1->SetBranchAddress("pCosAmp495", &pCosAmp495);
  t1->SetBranchAddress("pCosAmpErr495", &pCosAmpErr495);
  t1->SetBranchAddress("pCosAmp505", &pCosAmp505);
  t1->SetBranchAddress("pCosAmpErr505", &pCosAmpErr505);
  t1->SetBranchAddress("pCosAmp515", &pCosAmp515);
  t1->SetBranchAddress("pCosAmpErr515", &pCosAmpErr515);
  t1->SetBranchAddress("pCosAmp525", &pCosAmp525);
  t1->SetBranchAddress("pCosAmpErr525", &pCosAmpErr525);
  t1->SetBranchAddress("pCosAmp535", &pCosAmp535);
  t1->SetBranchAddress("pCosAmpErr535", &pCosAmpErr535);
  t1->SetBranchAddress("pCosAmp545", &pCosAmp545);
  t1->SetBranchAddress("pCosAmpErr545", &pCosAmpErr545);
  t1->SetBranchAddress("pCosAmp555", &pCosAmp555);
  t1->SetBranchAddress("pCosAmpErr555", &pCosAmpErr555);
  t1->SetBranchAddress("pCosAmp565", &pCosAmp565);
  t1->SetBranchAddress("pCosAmpErr565", &pCosAmpErr565);
  t1->SetBranchAddress("pCosAmp575", &pCosAmp575);
  t1->SetBranchAddress("pCosAmpErr575", &pCosAmpErr575);
  t1->SetBranchAddress("pCosAmp585", &pCosAmp585);
  t1->SetBranchAddress("pCosAmpErr585", &pCosAmpErr585);
  t1->SetBranchAddress("pCosAmp595", &pCosAmp595);
  t1->SetBranchAddress("pCosAmpErr595", &pCosAmpErr595);
  t1->SetBranchAddress("pCosAmp605", &pCosAmp605);
  t1->SetBranchAddress("pCosAmpErr605", &pCosAmpErr605);
  t1->SetBranchAddress("pCosAmp615", &pCosAmp615);
  t1->SetBranchAddress("pCosAmpErr615", &pCosAmpErr615);

  // Load values from tree and asign values back into an array
  for (Int_t k = 0; k < 10; k++){
    Parameter_Values->GetEntry(k);
    pValues425[k] = pCosAmp425;
    pErrValues425[k] = pCosAmpErr425;
    pValues435[k] = pCosAmp435;
    pErrValues435[k] = pCosAmpErr435;
    pValues445[k] = pCosAmp445;
    pErrValues445[k] = pCosAmpErr445;
    pValues455[k] = pCosAmp455;
    pErrValues455[k] = pCosAmpErr455;
    pValues465[k] = pCosAmp465;
    pErrValues465[k] = pCosAmpErr465;
    pValues475[k] = pCosAmp475;
    pErrValues475[k] = pCosAmpErr475;
    pValues485[k] = pCosAmp485;
    pErrValues485[k] = pCosAmpErr485;
    pValues495[k] = pCosAmp495;
    pErrValues495[k] = pCosAmpErr495;
    pValues505[k] = pCosAmp505;
    pErrValues505[k] = pCosAmpErr505;
    pValues515[k] = pCosAmp515;
    pErrValues515[k] = pCosAmpErr515;
    pValues525[k] = pCosAmp525;
    pErrValues525[k] = pCosAmpErr525;
    pValues535[k] = pCosAmp535;
    pErrValues535[k] = pCosAmpErr535;
    pValues545[k] = pCosAmp545;
    pErrValues545[k] = pCosAmpErr545;
    pValues555[k] = pCosAmp555;
    pErrValues555[k] = pCosAmpErr555;
    pValues565[k] = pCosAmp565;
    pErrValues565[k] = pCosAmpErr565;
    pValues575[k] = pCosAmp575;
    pErrValues575[k] = pCosAmpErr575;
    pValues585[k] = pCosAmp585;
    pErrValues585[k] = pCosAmpErr585;
    pValues595[k] = pCosAmp595;
    pErrValues595[k] = pCosAmpErr595;
    pValues605[k] = pCosAmp605;
    pErrValues605[k] = pCosAmpErr605;
    pValues615[k] = pCosAmp615;
    pErrValues615[k] = pCosAmpErr615;

  }

  TFile *f2= TFile::Open("/scratch/Mainz_Software/a2GoAT/LinPol_Aug16.root"); // Open linear polarisation plot

  // Calculate values of sigma for each angular and energy bin
  for (Int_t i = 0; i < 10; i++){

    pSigmaValues425[i] = pValues425[i]/(Graph->Eval(425,0));
    pSigmaErrValues425[i] = pErrValues425[i]/(Graph->Eval(425,0));
    pSigmaValues435[i] = pValues435[i]/(Graph->Eval(435,0));
    pSigmaErrValues435[i] = pErrValues435[i]/(Graph->Eval(435,0));
    pSigmaValues445[i] = pValues445[i]/(Graph->Eval(445,0));
    pSigmaErrValues445[i] = pErrValues445[i]/(Graph->Eval(445,0));
    pSigmaValues455[i] = pValues455[i]/(Graph->Eval(455,0));
    pSigmaErrValues455[i] = pErrValues455[i]/(Graph->Eval(455,0));
    pSigmaValues465[i] = pValues465[i]/(Graph->Eval(465,0));
    pSigmaErrValues465[i] = pErrValues465[i]/(Graph->Eval(465,0));
    pSigmaValues475[i] = pValues475[i]/(Graph->Eval(475,0));
    pSigmaErrValues475[i] = pErrValues475[i]/(Graph->Eval(475,0));
    pSigmaValues485[i] = pValues485[i]/(Graph->Eval(485,0));
    pSigmaErrValues485[i] = pErrValues485[i]/(Graph->Eval(485,0));
    pSigmaValues495[i] = pValues495[i]/(Graph->Eval(495,0));
    pSigmaErrValues495[i] = pErrValues495[i]/(Graph->Eval(495,0));
    pSigmaValues505[i] = pValues505[i]/(Graph->Eval(505,0));
    pSigmaErrValues505[i] = pErrValues505[i]/(Graph->Eval(505,0));
    pSigmaValues515[i] = pValues515[i]/(Graph->Eval(515,0));
    pSigmaErrValues515[i] = pErrValues515[i]/(Graph->Eval(515,0));
    pSigmaValues525[i] = pValues525[i]/(Graph->Eval(525,0));
    pSigmaErrValues525[i] = pErrValues525[i]/(Graph->Eval(525,0));
    pSigmaValues545[i] = pValues545[i]/(Graph->Eval(545,0));
    pSigmaErrValues545[i] = pErrValues545[i]/(Graph->Eval(545,0));
    pSigmaValues555[i] = pValues555[i]/(Graph->Eval(555,0));
    pSigmaErrValues555[i] = pErrValues555[i]/(Graph->Eval(555,0));
    pSigmaValues565[i] = pValues565[i]/(Graph->Eval(565,0));
    pSigmaErrValues565[i] = pErrValues565[i]/(Graph->Eval(565,0));
    pSigmaValues575[i] = pValues575[i]/(Graph->Eval(575,0));
    pSigmaErrValues575[i] = pErrValues575[i]/(Graph->Eval(575,0));
    pSigmaValues585[i] = pValues585[i]/(Graph->Eval(585,0));
    pSigmaErrValues585[i] = pErrValues585[i]/(Graph->Eval(585,0));
    pSigmaValues595[i] = pValues595[i]/(Graph->Eval(595,0));
    pSigmaErrValues595[i] = pErrValues595[i]/(Graph->Eval(595,0));
    pSigmaValues605[i] = pValues605[i]/(Graph->Eval(605,0));
    pSigmaErrValues605[i] = pErrValues605[i]/(Graph->Eval(605,0));
    pSigmaValues615[i] = pValues615[i]/(Graph->Eval(615,0));
    pSigmaErrValues615[i] = pErrValues615[i]/(Graph->Eval(615,0));

  }

  TFile f3("Sigma_Plots_NoScatt_1.root", "RECREATE");

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
  hr->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 420-430MeV)");

  gr1 = new TGraphErrors(10, x, pSigmaValues425, ex, pSigmaErrValues425);
  gr1->SetMarkerColor(2);
  gr1->SetMarkerStyle(5);
  gr1->SetMarkerSize(2);
  gr1->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 420-430MeV)");
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
  hr1->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 430-440MeV)");

  gr2 = new TGraphErrors(10, x, pSigmaValues435 , ex, pSigmaErrValues435);
  gr2->SetMarkerColor(2);
  gr2->SetMarkerStyle(5);
  gr2->SetMarkerSize(2);
  gr2->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 430-440MeV)");
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
  hr2->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 440-450MeV)");

  gr3 = new TGraphErrors(10, x, pSigmaValues445 , ex, pSigmaErrValues445);
  gr3->SetMarkerColor(2);
  gr3->SetMarkerStyle(5);
  gr3->SetMarkerSize(2);
  gr3->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 440-450MeV)");
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
  hr3->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 450-460MeV)");

  gr4 = new TGraphErrors(10, x, pSigmaValues455 , ex, pSigmaErrValues455);
  gr4->SetMarkerColor(2);
  gr4->SetMarkerStyle(5);
  gr4->SetMarkerSize(2);
  gr4->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 450-460MeV)");
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
  hr4->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 460-470MeV)");

  gr5 = new TGraphErrors(10, x, pSigmaValues465, ex, pSigmaErrValues465);
  gr5->SetMarkerColor(2);
  gr5->SetMarkerStyle(5);
  gr5->SetMarkerSize(2);
  gr5->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 460-470MeV)");
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
  hr5->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 470-480MeV)");

  gr6 = new TGraphErrors(10, x, pSigmaValues475, ex, pSigmaErrValues475);
  gr6->SetMarkerColor(2);
  gr6->SetMarkerStyle(5);
  gr6->SetMarkerSize(2);
  gr6->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 470-480MeV)");
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
  hr6->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 480-490MeV)");

  gr7 = new TGraphErrors(10, x, pSigmaValues485, ex, pSigmaErrValues485);
  gr7->SetMarkerColor(2);
  gr7->SetMarkerStyle(5);
  gr7->SetMarkerSize(2);
  gr7->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 480-490MeV)");
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
  hr7->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 490-500MeV)");

  gr8 = new TGraphErrors(10, x, pSigmaValues495, ex, pSigmaErrValues495);
  gr8->SetMarkerColor(2);
  gr8->SetMarkerStyle(5);
  gr8->SetMarkerSize(2);
  gr8->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 490-500MeV)");
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
  hr8->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 500-510MeV)");

  gr9 = new TGraphErrors(10, x, pSigmaValues505, ex, pSigmaErrValues505);
  gr9->SetMarkerColor(2);
  gr9->SetMarkerStyle(5);
  gr9->SetMarkerSize(2);
  gr9->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 500-510MeV)");
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
  hr9->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 510-520MeV)");

  gr10 = new TGraphErrors(10, x, pSigmaValues515, ex, pSigmaErrValues515);
  gr10->SetMarkerColor(2);
  gr10->SetMarkerStyle(5);
  gr10->SetMarkerSize(2);
  gr10->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 510-520MeV)");
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
  hr10->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 520-530MeV)");

  gr11 = new TGraphErrors(10, x, pSigmaValues525, ex, pSigmaErrValues525);
  gr11->SetMarkerColor(2);
  gr11->SetMarkerStyle(5);
  gr11->SetMarkerSize(2);
  gr11->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 520-530MeV)");
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
  hr11->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 530-540MeV)");

  gr12 = new TGraphErrors(10, x, pSigmaValues535, ex, pSigmaErrValues535);
  gr12->SetMarkerColor(2);
  gr12->SetMarkerStyle(5);
  gr12->SetMarkerSize(2);
  gr12->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 530-540MeV)");
  gr12->GetXaxis()->SetTitle("Cos#theta_{CM}");
  gr12->GetYaxis()->SetTitle("#Sigma");
  gr12->Draw("ep");

  TCanvas *canvas12 = new TCanvas("canvas12","canvas12", 1920, 1080);
  TPad *pad13 = new TPad("pad13","",0,0,1,1);
  pad13->Draw();
  pad13->cd();

  pad13->SetTickx(1);
  pad13->SetTicky(1);
  pad13->SetGridx(1);
  pad13->SetGridy(1);
  TH1F  *hr12;
  hr12 = canvas12->DrawFrame(xMin,-1.5,xMax,1.5);
  hr12->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 540-550MeV)");

  gr13 = new TGraphErrors(10, x, pSigmaValues545, ex, pSigmaErrValues545);
  gr13->SetMarkerColor(2);
  gr13->SetMarkerStyle(5);
  gr13->SetMarkerSize(2);
  gr13->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 540-550MeV)");
  gr13->GetXaxis()->SetTitle("Cos#theta_{CM}");
  gr13->GetYaxis()->SetTitle("#Sigma");
  gr13->Draw("ep");

  TCanvas *canvas13 = new TCanvas("canvas13","canvas13", 1920, 1080);
  TPad *pad14 = new TPad("pad14","",0,0,1,1);
  pad14->Draw();
  pad14->cd();

  pad14->SetTickx(1);
  pad14->SetTicky(1);
  pad14->SetGridx(1);
  pad14->SetGridy(1);
  TH1F  *hr13;
  hr13 = canvas13->DrawFrame(xMin,-1.5,xMax,1.5);
  hr13->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 550-560MeV)");

  gr14 = new TGraphErrors(10, x, pSigmaValues555, ex, pSigmaErrValues555);
  gr14->SetMarkerColor(2);
  gr14->SetMarkerStyle(5);
  gr14->SetMarkerSize(2);
  gr14->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 550-560MeV)");
  gr14->GetXaxis()->SetTitle("Cos#theta_{CM}");
  gr14->GetYaxis()->SetTitle("#Sigma");
  gr14->Draw("ep");

  TCanvas *canvas14 = new TCanvas("canvas14","canvas14", 1920, 1080);
  TPad *pad15 = new TPad("pad15","",0,0,1,1);
  pad15->Draw();
  pad15->cd();

  pad15->SetTickx(1);
  pad15->SetTicky(1);
  pad15->SetGridx(1);
  pad15->SetGridy(1);
  TH1F  *hr14;
  hr14 = canvas14->DrawFrame(xMin,-1.5,xMax,1.5);
  hr14->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 560-570MeV)");

  gr15 = new TGraphErrors(10, x, pSigmaValues565, ex, pSigmaErrValues565);
  gr15->SetMarkerColor(2);
  gr15->SetMarkerStyle(5);
  gr15->SetMarkerSize(2);
  gr15->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 560-570MeV)");
  gr15->GetXaxis()->SetTitle("Cos#theta_{CM}");
  gr15->GetYaxis()->SetTitle("#Sigma");
  gr15->Draw("ep");

  TCanvas *canvas15 = new TCanvas("canvas15","canvas15", 1920, 1080);
  TPad *pad16 = new TPad("pad16","",0,0,1,1);
  pad16->Draw();
  pad16->cd();

  pad16->SetTickx(1);
  pad16->SetTicky(1);
  pad16->SetGridx(1);
  pad16->SetGridy(1);
  TH1F  *hr15;
  hr15 = canvas15->DrawFrame(xMin,-1.5,xMax,1.5);
  hr15->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 570-580MeV)");

  gr16 = new TGraphErrors(10, x, pSigmaValues575, ex, pSigmaErrValues575);
  gr16->SetMarkerColor(2);
  gr16->SetMarkerStyle(5);
  gr16->SetMarkerSize(2);
  gr16->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 570-580MeV)");
  gr16->GetXaxis()->SetTitle("Cos#theta_{CM}");
  gr16->GetYaxis()->SetTitle("#Sigma");
  gr16->Draw("ep");

  TCanvas *canvas16 = new TCanvas("canvas16","canvas16", 1920, 1080);
  TPad *pad17 = new TPad("pad17","",0,0,1,1);
  pad17->Draw();
  pad17->cd();

  pad17->SetTickx(1);
  pad17->SetTicky(1);
  pad17->SetGridx(1);
  pad17->SetGridy(1);
  TH1F  *hr16;
  hr16 = canvas16->DrawFrame(xMin,-1.5,xMax,1.5);
  hr16->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 580-590MeV)");

  gr17 = new TGraphErrors(10, x, pSigmaValues585, ex, pSigmaErrValues585);
  gr17->SetMarkerColor(2);
  gr17->SetMarkerStyle(5);
  gr17->SetMarkerSize(2);
  gr17->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 580-590MeV)");
  gr17->GetXaxis()->SetTitle("Cos#theta_{CM}");
  gr17->GetYaxis()->SetTitle("#Sigma");
  gr17->Draw("ep");

  TCanvas *canvas17 = new TCanvas("canvas17","canvas17", 1920, 1080);
  TPad *pad18 = new TPad("pad18","",0,0,1,1);
  pad18->Draw();
  pad18->cd();

  pad18->SetTickx(1);
  pad18->SetTicky(1);
  pad18->SetGridx(1);
  pad18->SetGridy(1);
  TH1F  *hr17;
  hr17 = canvas17->DrawFrame(xMin,-1.5,xMax,1.5);
  hr17->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 590-600MeV)");

  gr18 = new TGraphErrors(10, x, pSigmaValues595, ex, pSigmaErrValues595);
  gr18->SetMarkerColor(2);
  gr18->SetMarkerStyle(5);
  gr18->SetMarkerSize(2);
  gr18->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 590-600MeV)");
  gr18->GetXaxis()->SetTitle("Cos#theta_{CM}");
  gr18->GetYaxis()->SetTitle("#Sigma");
  gr18->Draw("ep");

  TCanvas *canvas18 = new TCanvas("canvas18","canvas18", 1920, 1080);
  TPad *pad19 = new TPad("pad19","",0,0,1,1);
  pad19->Draw();
  pad19->cd();

  pad19->SetTickx(1);
  pad19->SetTicky(1);
  pad19->SetGridx(1);
  pad19->SetGridy(1);
  TH1F  *hr18;
  hr18 = canvas18->DrawFrame(xMin,-1.5,xMax,1.5);
  hr18->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 600-610MeV)");

  gr19 = new TGraphErrors(10, x, pSigmaValues605, ex, pSigmaErrValues605);
  gr19->SetMarkerColor(2);
  gr19->SetMarkerStyle(5);
  gr19->SetMarkerSize(2);
  gr19->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 600-610MeV)");
  gr19->GetXaxis()->SetTitle("Cos#theta_{CM}");
  gr19->GetYaxis()->SetTitle("#Sigma");
  gr19->Draw("ep");

  TCanvas *canvas19 = new TCanvas("canvas19","canvas19", 1920, 1080);
  TPad *pad20 = new TPad("pad20","",0,0,1,1);
  pad20->Draw();
  pad20->cd();

  pad20->SetTickx(1);
  pad20->SetTicky(1);
  pad20->SetGridx(1);
  pad20->SetGridy(1);
  TH1F  *hr19;
  hr19 = canvas19->DrawFrame(xMin,-1.5,xMax,1.5);
  hr19->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 610-620MeV)");

  gr20 = new TGraphErrors(10, x, pSigmaValues615, ex, pSigmaErrValues615);
  gr20->SetMarkerColor(2);
  gr20->SetMarkerStyle(5);
  gr20->SetMarkerSize(2);
  gr20->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 610-620MeV)");
  gr20->GetXaxis()->SetTitle("Cos#theta_{CM}");
  gr20->GetYaxis()->SetTitle("#Sigma");
  gr20->Draw("ep");

  //TCanvas *canvas20 = new TCanvas("canvas20","canvas20", 1920, 1080);
  //canvas20->Divide(5,4);
  //canvas20->cd(1);
  //pad1->Draw();
  //canvas20->cd(2);
  //pad2->Draw();
  //canvas20->cd(3);
  //pad3->Draw();
  //canvas20->cd(4);
  //pad4->Draw();
  //canvas20->cd(5);
  //pad5->Draw();
  //canvas20->cd(6);
  //pad6->Draw();
  //canvas20->cd(7);
  //pad7->Draw();
  //canvas20->cd(8);
  //pad8->Draw();
  //canvas20->cd(9);
  //pad9->Draw();
  //canvas20->cd(10);
  //pad10->Draw();
  //canvas20->cd(11);
  //pad11->Draw();
  //canvas20->cd(12);
  //pad12->Draw();
  //canvas20->cd(13);
  //pad13->Draw();
  //canvas20->cd(14);
  //pad14->Draw();
  //canvas20->cd(15);
  //pad15->Draw();
  //canvas20->cd(16);
  //pad16->Draw();
  //canvas20->cd(17);
  //pad17->Draw();
  //canvas20->cd(18);
  //pad18->Draw();
  //canvas20->cd(19);
  //pad19->Draw();
  //canvas20->cd(20);
  //pad20->Draw();

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
  canvas12->Write();
  canvas13->Write();
  canvas14->Write();
  canvas15->Write();
  canvas16->Write();
  canvas17->Write();
  canvas18->Write();
  canvas19->Write();
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
  gr13->Write();
  gr14->Write();
  gr15->Write();
  gr16->Write();
  gr17->Write();
  gr18->Write();
  gr19->Write();
  gr20->Write();
  //canvas20->Write();
  f3.Write();

}
