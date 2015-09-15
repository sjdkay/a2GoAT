#include "./includes.h"

// define a function with 3 parameters
Double_t fitf(Double_t *x,Double_t *par)
{
    Double_t fitval = 0;
    fitval =  par[0] * (1+(par[1]*cos(x[0]*((TMath::Pi())/180))));
    return fitval;
}

void CosFitMCOutCut(){

  double x[10] = {125, 175, 225, 275, 325, 375, 425, 475, 525, 575};
  double xLiu[10]={200, 248, 296, 346, 401, 355, 407, 460, 401, 436};
  double xLiuErr[10]={7, 8, 10, 12, 15, 11, 13, 16, 12, 14};
  double xIkeda[8]={403.3, 446.6, 503.3, 547.6, 599, 511.4, 549, 597.8};
  double PrevDataLiu[10]={-0.25, -0.21, -0.13, -0.22, -0.47, -0.35, -0.52, -0.59, -0.44, -0.58};
  double PrevDataIkeda[8]={-0.44, -0.54, -0.61, -0.62, -0.42, -0.77, -0.51, -0.58};
  double PrevDataLiuErr[10]={0.12, 0.21, 0.07, 0.08, 0.08, 0.09, 0.09, 0.13, 0.06, 0.08};
  double PrevDataIkedaErr[8]={0.09, 0.1, 0.11, 0.11, 0.11, 0.13, 0.08, 0.12};
  double Amp[10];
  double AmpErr[10];
  double Y_Off[10];
  double Y_OffErr[10];
  double APow = 0.5;
  double PolCut[10];
  double PolErrCut[10];
  double Y_OffCorr[10];
  double Y_OffCorrErr[10];

  TF1 *CosFit = new TF1("CosFit",  fitf, -180.0, 180.0, 2);
  CosFit->SetParLimits(0, 0, 1000);
  CosFit->SetParLimits(1, -1, 1);
  CosFit->SetParNames("Y_Offset", "Amplitdue");

  TFile *f = new TFile("Physics_10e7_6_12_08_15_2.root");
  TText *warn = new TText(0, 0 ,"PRELIMINARY");

  for(Int_t i = 0; i < 10; i++){

    if(i==0){
        Char_t* Title = "PhiScatt in Scattered Proton Frame at 125 +/- 25 MeV; PhiScattCut125MeV";
        TH1D *hist = (TH1D*)f->Get("Phi_Scattered_125MeV_Cut_Out");
        Char_t* GraphPDF = "./CosFit_Cut_125MeV_Out.pdf";
        Char_t* GraphPNG = "./CosFit_Cut_125MeV_Out.png";
        RebinVal = 1;
        Float_t yMax = 50;
      }

    if(i==1){
        Char_t* Title = "PhiScatt in Scattered Proton Frame at 175 +/- 25 MeV; PhiScattCut175MeV";
        TH1D *hist = (TH1D*)f->Get("Phi_Scattered_175MeV_Cut_Out");
        Char_t* GraphPDF = "./CosFit_Cut_175MeV_Out.pdf";
        Char_t* GraphPNG = "./CosFit_Cut_175MeV_Out.png";
        RebinVal = 1;
        Float_t yMax = 50;
      }

    if(i==2){
        Char_t* Title = "PhiScatt in Scattered Proton Frame at 225 +/- 25 MeV; PhiScattCut225MeV";
        TH1D *hist = (TH1D*)f->Get("Phi_Scattered_225MeV_Cut_Out");
        Char_t* GraphPDF = "./CosFit_Cut_225MeV_Out.pdf";
        Char_t* GraphPNG = "./CosFit_Cut_225MeV_Out.png";
        RebinVal = 1;
        Float_t yMax = 30;
      }

    if(i==3){
        Char_t* Title = "PhiScatt in Scattered Proton Frame at 275 +/- 25 MeV; PhiScattCut275MeV";
        TH1D *hist = (TH1D*)f->Get("Phi_Scattered_275MeV_Cut_Out");
        Char_t* GraphPDF = "./CosFit_Cut_275MeV_Out.pdf";
        Char_t* GraphPNG = "./CosFit_Cut_275MeV_Out.png";
        RebinVal = 4;
        Float_t yMax = 80;
      }

    if(i==4){
        Char_t* Title = "PhiScatt in Scattered Proton Frame at 325 +/- 25 MeV; PhiScattCut325MeV";
        TH1D *hist = (TH1D*)f->Get("Phi_Scattered_325MeV_Cut_Out");
        Char_t* GraphPDF = "./CosFit_Cut_325MeV_Out.pdf";
        Char_t* GraphPNG = "./CosFit_Cut_325MeV_Out.png";
        RebinVal = 4;
        Float_t yMax = 80;
      }

    if(i==5){
        Char_t* Title = "PhiScatt in Scattered Proton Frame at 375 +/- 25 MeV; PhiScattCut375MeV";
        TH1D *hist = (TH1D*)f->Get("Phi_Scattered_375MeV_Cut_Out");
        Char_t* GraphPDF = "./CosFit_Cut_375MeV_Out.pdf";
        Char_t* GraphPNG = "./CosFit_Cut_375MeV_Out.png";
        RebinVal = 4;
        Float_t yMax = 60;
      }

    if(i==6){
        Char_t* Title = "PhiScatt in Scattered Proton Frame at 425 +/- 25 MeV; PhiScattCut425MeV";
        TH1D *hist = (TH1D*)f->Get("Phi_Scattered_425MeV_Cut_Out");
        Char_t* GraphPDF = "./CosFit_Cut_425MeV_Out.pdf";
        Char_t* GraphPNG = "./CosFit_Cut_425MeV_Out.png";
        RebinVal = 4;
        Float_t yMax = 30;
      }

    if(i==7){
        Char_t* Title = "PhiScatt in Scattered Proton Frame at 475 +/- 25 MeV; PhiScattCut475MeV";
        TH1D *hist = (TH1D*)f->Get("Phi_Scattered_475MeV_Cut_Out");
        Char_t* GraphPDF = "./CosFit_Cut_475MeV_Out.pdf";
        Char_t* GraphPNG = "./CosFit_Cut_475MeV_Out.png";
        RebinVal = 4;
        Float_t yMax = 20;
      }

    if(i==8){
        Char_t* Title = "PhiScatt in Scattered Proton Frame at 525 +/- 25 MeV; PhiScattCut525MeV";
        TH1D *hist = (TH1D*)f->Get("Phi_Scattered_525MeV_Cut_Out");
        Char_t* GraphPDF = "./CosFit_Cut_525MeV_Out.pdf";
        Char_t* GraphPNG = "./CosFit_Cut_525MeV_Out.png";
        RebinVal = 4;
        Float_t yMax = 20;
      }

    if(i==9){
        Char_t* Title = "PhiScatt in Scattered Proton Frame at 575 +/- 25 MeV; PhiScattCut575MeV";
        TH1D *hist = (TH1D*)f->Get("Phi_Scattered_575MeV_Cut_Out");
        Char_t* GraphPDF = "./CosFit_Cut_575MeV_Out.pdf";
        Char_t* GraphPNG = "./CosFit_Cut_575MeV_Out.png";
        RebinVal = 4;
        Float_t yMax = 20;
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
    //if (i == 9) Float_t yMin = -10;
    //if (i == 0) Float_t yMin = -75;

    strcpy(hrTitle, Title);
    hr = canvas->DrawFrame(xMin,yMin,xMax,yMax);
    hr->SetTitle(hrTitle);

    hist->SetMarkerStyle(1);
    hist->SetLineColor(2);
    hist->Rebin(RebinVal);
    hist->Draw("EHISTSAMES");
    hist->Fit("CosFit", "LL");
    CosFit->SetLineColor(4);
    CosFit->Draw("SAMES");
    warn->Draw("SAMES");
    gStyle->SetOptFit(0111);
    gPad->Update();

    Amp[i] = CosFit->GetParameter(1);
    AmpErr[i] = CosFit->GetParError(1);
    Y_Off[i]  = CosFit->GetParameter(0);
    Y_OffErr[i] = CosFit->GetParError(0);

    BinWidth = RebinVal*10;

    PolCut[i] = Amp[i]/APow;
    PolErrCut[i] = AmpErr[i]/APow;
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

  strcpy(hrTitle, "Polarisation as a Function of Photon Energy with Theta Neutron Cut");
  hr = canvas->DrawFrame(xMin,yMin,xMax,yMax);
  hr->SetTitle(hrTitle);

  TGraphErrors *PolEGammaCut;
  TGraphErrors *LiuData;
  TGraphErrors *IkedaData;

  PolEGammaCut  = new TGraphErrors(10,x,PolCut,0,PolErrCut);
  LiuData = new TGraphErrors(8, xLiu, PrevDataLiu, xLiuErr, PrevDataLiuErr);
  IkedaData = new TGraphErrors(10, xIkeda, PrevDataIkeda,0 , PrevDataIkedaErr);
  PolEGammaCut->SetMarkerColor(1);
  PolEGammaCut->SetMarkerStyle(5);
  PolEGammaCut->SetMarkerSize(0.5);
  PolEGammaCut->Draw("P");
  PolEGammaCut->Draw("E1");
  LiuData->SetMarkerColor(2);
  LiuData->SetMarkerStyle(21);
  LiuData->SetMarkerSize(0.5);
  LiuData->Draw("ESAMEP");
  IkedaData->SetMarkerColor(4);
  IkedaData->SetMarkerStyle(22);
  IkedaData->SetMarkerSize(0.5);
  IkedaData->Draw("ESAMEP");
  leg = new TLegend(0.75, 0.75, 0.95, 0.95);
  leg->AddEntry(PolEGammaCut, "Current Data", "ep");
  leg->AddEntry(LiuData, "Liu 1968", "ep");
  leg->AddEntry(IkedaData, "Ikeda 1980", "ep");
  leg->Draw("Same");
  warn->Draw("SAME");
  canvas->SaveAs("./PolEGamma_MCOutCut.pdf");
  canvas->SaveAs("./PolEGamma_MCOutCut.png");

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
  yMax = 30;

  strcpy(hrTitle, "Y_Offset as a Function of Photon Energy with Theta Neutron Cut");
  hr = canvas->DrawFrame(xMin,yMin,xMax,yMax);
  hr->SetTitle(hrTitle);

  TGraphErrors *Y_OffsetEGammaCut;

  Y_OffsetEGammaCut  = new TGraphErrors(10,x,Y_OffCorr,0,Y_OffCorrErr);
  Y_OffsetEGammaCut->SetMarkerColor(4);
  Y_OffsetEGammaCut->SetMarkerStyle(5);
  Y_OffsetEGammaCut->SetMarkerSize(0.5);
  Y_OffsetEGammaCut->Draw("P");
  Y_OffsetEGammaCut->Draw("E1");
  warn->Draw("SAME");
  canvas->SaveAs("./Y_OffsetEGamma_MCOutCut.pdf");
  canvas->SaveAs("./Y_OffsetEGamma_MCOutCut.png");

}
