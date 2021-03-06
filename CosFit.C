#include "./includes.h"

// define a function with 3 parameters
Double_t fitf(Double_t *x,Double_t *par)
{
    Double_t fitval = 0;
    fitval =  par[0] * (1+(par[1]*cos(x[0]*((TMath::Pi())/180))));
    return fitval;
}

void CosFit(){

  double x[10] = {125, 175, 225, 275, 325, 375, 425, 475, 525, 575};
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
  double APow = 0.5;
  double Pol[10];
  double PolErr[10];
  double Y_OffCorr[10];
  double Y_OffCorrErr[10];

  TF1 *CosFit = new TF1("CosFit",  fitf, -180.0, 180.0, 2);
  CosFit->SetParLimits(0, -1000, 1000);
  CosFit->SetParLimits(1, -1, 1);
  CosFit->SetParNames("Y_Offset", "Amplitdue");

  TFile *f = new TFile("PhysicsTotal32_03_06_15.root");
  TText *warn = new TText(0, 0 ,"PRELIMINARY");

  for(Int_t i = 0; i < 10; i++){

    if(i==0){
        Char_t* Title = "PhiScatt in Scattered Proton Frame at 125 +/- 25 MeV; PhiScatt125MeV";
        TH1D *hist = (TH1D*)f->Get("Phi_Scattered_125MeV");
        Char_t* GraphPDF = "./CosFit_125MeV.pdf";
        Char_t* GraphPNG = "./CosFit_125MeV.png";
        RebinVal = 2;
        Float_t yMax = 200;
      }

    if(i==1){
        Char_t* Title = "PhiScatt in Scattered Proton Frame at 175 +/- 25 MeV; PhiScatt175MeV";
        TH1D *hist = (TH1D*)f->Get("Phi_Scattered_175MeV");
        Char_t* GraphPDF = "./CosFit_175MeV.pdf";
        Char_t* GraphPNG = "./CosFit_175MeV.png";
        RebinVal = 2;
        Float_t yMax = 250;
      }

    if(i==2){
        Char_t* Title = "PhiScatt in Scattered Proton Frame at 225 +/- 25 MeV; PhiScatt225MeV";
        TH1D *hist = (TH1D*)f->Get("Phi_Scattered_225MeV");
        Char_t* GraphPDF = "./CosFit_225MeV.pdf";
        Char_t* GraphPNG = "./CosFit_225MeV.png";
        RebinVal = 5;
        Float_t yMax = 550;
      }

    if(i==3){
        Char_t* Title = "PhiScatt in Scattered Proton Frame at 275 +/- 25 MeV; PhiScatt275MeV";
        TH1D *hist = (TH1D*)f->Get("Phi_Scattered_275MeV");
        Char_t* GraphPDF = "./CosFit_275MeV.pdf";
        Char_t* GraphPNG = "./CosFit_275MeV.png";
        RebinVal = 2;
        Float_t yMax = 250;
      }

    if(i==4){
        Char_t* Title = "PhiScatt in Scattered Proton Frame at 325 +/- 25 MeV; PhiScatt325MeV";
        TH1D *hist = (TH1D*)f->Get("Phi_Scattered_325MeV");
        Char_t* GraphPDF = "./CosFit_325MeV.pdf";
        Char_t* GraphPNG = "./CosFit_325MeV.png";
        RebinVal = 5;
        Float_t yMax = 250;
      }

    if(i==5){
        Char_t* Title = "PhiScatt in Scattered Proton Frame at 375 +/- 25 MeV; PhiScatt375MeV";
        TH1D *hist = (TH1D*)f->Get("Phi_Scattered_375MeV");
        Char_t* GraphPDF = "./CosFit_375MeV.pdf";
        Char_t* GraphPNG = "./CosFit_375MeV.png";
        RebinVal = 10;
        Float_t yMax = 250;
      }

    if(i==6){
        Char_t* Title = "PhiScatt in Scattered Proton Frame at 425 +/- 25 MeV; PhiScatt425MeV";
        TH1D *hist = (TH1D*)f->Get("Phi_Scattered_425MeV");
        Char_t* GraphPDF = "./CosFit_425MeV.pdf";
        Char_t* GraphPNG = "./CosFit_425MeV.png";
        RebinVal = 20;
        Float_t yMax = 200;
      }

    if(i==7){
        Char_t* Title = "PhiScatt in Scattered Proton Frame at 475 +/- 25 MeV; PhiScatt475MeV";
        TH1D *hist = (TH1D*)f->Get("Phi_Scattered_475MeV");
        Char_t* GraphPDF = "./CosFit_475MeV.pdf";
        Char_t* GraphPNG = "./CosFit_475MeV.png";
        RebinVal = 10;
        Float_t yMax = 80;
      }

    if(i==8){
        Char_t* Title = "PhiScatt in Scattered Proton Frame at 525 +/- 25 MeV; PhiScatt525MeV";
        TH1D *hist = (TH1D*)f->Get("Phi_Scattered_525MeV");
        Char_t* GraphPDF = "./CosFit_525MeV.pdf";
        Char_t* GraphPNG = "./CosFit_525MeV.png";
        RebinVal = 20;
        Float_t yMax = 60;
      }

    if(i==9){
        Char_t* Title = "PhiScatt in Scattered Proton Frame at 575 +/- 25 MeV; PhiScatt575MeV";
        TH1D *hist = (TH1D*)f->Get("Phi_Scattered_575MeV");
        Char_t* GraphPDF = "./CosFit_575MeV.pdf";
        Char_t* GraphPNG = "./CosFit_575MeV.png";
        RebinVal = 10;
        Float_t yMax = 30;
      }

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

    Float_t xMin = -180;
    Float_t xMax = 180;
    Float_t yMin = 0;
    if (i == 9) Float_t yMin = -10;
    if (i == 0) Float_t yMin = -75;

    strcpy(hrTitle, Title);
    hr = canvas->DrawFrame(xMin,yMin,xMax,yMax);
    hr->SetTitle(hrTitle);

    hist->SetMarkerStyle(1);
    hist->SetLineColor(2);
    hist->Rebin(RebinVal);
    hist->Draw("EHISTSAMES");
    hist->Fit("CosFit");
    CosFit->SetLineColor(4);
    CosFit->Draw("SAMES");
    warn->Draw("SAMES");
    gStyle->SetOptFit(0111);
    gPad->Update();

    Amp[i] = CosFit->GetParameter(1);
    AmpErr[i] = CosFit->GetParError(1);
    Y_Off[i]  = CosFit->GetParameter(0);
    Y_OffErr[i] = CosFit->GetParError(0);

    BinWidth = RebinVal;

    Pol[i] = Amp[i]/APow;
    PolErr[i] = AmpErr[i]/APow;
    Y_OffCorr[i] = Y_Off[i]/BinWidth;
    Y_OffCorrErr[i] = Y_OffCorr[i]/BinWidth;

    canvas->SaveAs(filename = GraphPDF);
    canvas->SaveAs(filename = GraphPNG);

  }

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

  xMin = 0;
  xMax = 700;
  yMin = -2;
  yMax = 1;

  strcpy(hrTitle, "Polarisation as a Function of Photon Energy");
  hr = canvas->DrawFrame(xMin,yMin,xMax,yMax);
  hr->SetTitle(hrTitle);

  TGraphErrors *PolEGamma;
  TGraphErrors *LiuData;
  TGraphErrors *IkedaData;

  PolEGamma  = new TGraphErrors(10,x,Pol,0,PolErr);
  LiuData = new TGraphErrors(19, xLiu, PrevDataLiu, xLiuErr, PrevDataLiuErr);
  IkedaData = new TGraphErrors(19, xIkeda, PrevDataIkeda,0 , PrevDataIkedaErr);
  PolEGamma->SetMarkerColor(1);
  PolEGamma->SetMarkerStyle(5);
  PolEGamma->SetMarkerSize(0.5);
  PolEGamma->Draw("P");
  PolEGamma->Draw("E1");
  LiuData->SetMarkerColor(2);
  LiuData->SetMarkerStyle(21);
  LiuData->SetMarkerSize(0.5);
  LiuData->Draw("ESAMEP");
  IkedaData->SetMarkerColor(4);
  IkedaData->SetMarkerStyle(22);
  IkedaData->SetMarkerSize(0.5);
  IkedaData->Draw("ESAMEP");
  warn->Draw("SAME");
  leg = new TLegend(0.75, 0.75, 0.95, 0.95);
  leg->AddEntry(PolEGamma, "Current Data", "ep");
  leg->AddEntry(LiuData, "Liu 1968", "ep");
  leg->AddEntry(IkedaData, "Ikeda 1980", "ep");
  leg->Draw("Same");
  canvas->SaveAs("./PolEGamma.pdf");
  canvas->SaveAs("./PolEGamma.png");

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

  xMin = 0;
  xMax = 600;
  yMin = 0;
  yMax = 75;


  strcpy(hrTitle, "Y_Offset as a Function of Photon Energy");
  hr = canvas->DrawFrame(xMin,yMin,xMax,yMax);
  hr->SetTitle(hrTitle);

  TGraphErrors *Y_OffsetEGamma;

  Y_OffsetEGamma = new TGraphErrors(10,x,Y_OffCorr,0,Y_OffCorrErr);
  Y_OffsetEGamma->SetMarkerColor(4);
  Y_OffsetEGamma->SetMarkerStyle(5);
  Y_OffsetEGamma->SetMarkerSize(0.5);
  Y_OffsetEGamma->Draw("P");
  Y_OffsetEGamma->Draw("E1");
  warn->Draw("SAME");
  canvas->SaveAs("./Y_OffsetEGamma.pdf");
  canvas->SaveAs("./Y_OffsetEGamma.png");

}
