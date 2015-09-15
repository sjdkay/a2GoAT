#include "./includes.h"


void RatioMCIn(){

  Int_t nBins;
  double LowMax;
  double HighMax;
  double AvgMax;
  double ratio[10][36];
  double P1[10];
  double P1Err[10];
  double P2[10];
  double P2Err[10];
  double P3[10];
  double P3Err[10];
  double Par1;
  double Par1Err;
  double Par2;
  double Par2Err;
  double Par3;
  double Par3Err;

  TFile *f = new TFile("Physics_10e7_6_12_08_15_2.root");
  TText *warn = new TText(0, 0 ,"PRELIMINARY");

  for(Int_t i = 0; i < 10; i++){

    if(i==0){
        Char_t* Title = "Ratio of Angle bin to maximum for 125 +/- 25 MeV; Ratio125MeV_In";
        TH1D *hist = (TH1D*)f->Get("Phi_Scattered_125MeV_In");
        Char_t* RatioPDF = "./Ratio_125MeV_In.pdf";
        Char_t* RatioPNG = "./Ratio_125MeV_In.png";
        RebinVal = 1;
      }

    if(i==1){

        Char_t* Title = "Ratio of Angle bin to maximum for 175 +/- 25 MeV; Ratio175MeV_In";
        TH1D *hist = (TH1D*)f->Get("Phi_Scattered_175MeV_In");
        Char_t* RatioPDF = "./Ratio_175MeV_In.pdf";
        Char_t* RatioPNG = "./Ratio_175MeV_In.png";
        RebinVal = 1;
      }

    if(i==2){

        Char_t* Title = "Ratio of Angle bin to maximum for 225 +/- 25 MeV; Ratio225MeV_In";
        TH1D *hist = (TH1D*)f->Get("Phi_Scattered_225MeV_In");
        Char_t* RatioPDF = "./Ratio_225MeV_In.pdf";
        Char_t* RatioPNG = "./Ratio_225MeV_In.png";
        RebinVal = 1;
      }

    if(i==3){
        Char_t* Title = "Ratio of Angle bin to maximum for 275 +/- 25 MeV; Ratio275MeV_In";
        TH1D *hist = (TH1D*)f->Get("Phi_Scattered_275MeV_In");
        Char_t* RatioPDF = "./Ratio_275MeV_In.pdf";
        Char_t* RatioPNG = "./Ratio_275MeV_In.png";
        RebinVal = 1;
      }

    if(i==4){
        Char_t* Title = "Ratio of Angle bin to maximum for 325 +/- 25 MeV; Ratio325MeV_In";
        TH1D *hist = (TH1D*)f->Get("Phi_Scattered_325MeV_In");
        Char_t* RatioPDF = "./Ratio_325MeV_In.pdf";
        Char_t* RatioPNG = "./Ratio_325MeV_In.png";
        RebinVal = 1;
      }

    if(i==5){
        Char_t* Title = "Ratio of Angle bin to maximum for 375 +/- 25 MeV; Ratio375MeV_In";
        TH1D *hist = (TH1D*)f->Get("Phi_Scattered_375MeV_In");
        Char_t* RatioPDF = "./Ratio_375MeV_In.pdf";
        Char_t* RatioPNG = "./Ratio_375MeV_In.png";
        RebinVal = 2;
      }

    if(i==6){
        Char_t* Title = "Ratio of Angle bin to maximum for 425 +/- 25 MeV; Ratio425MeV_In";
        TH1D *hist = (TH1D*)f->Get("Phi_Scattered_425MeV_In");
        Char_t* RatioPDF = "./Ratio_425MeV_In.pdf";
        Char_t* RatioPNG = "./Ratio_425MeV_In.png";
        RebinVal = 2;
      }

    if(i==7){
        Char_t* Title = "Ratio of Angle bin to maximum for 475 +/- 25 MeV; Ratio475MeV_In";
        TH1D *hist = (TH1D*)f->Get("Phi_Scattered_475MeV_In");
        Char_t* RatioPDF = "./Ratio_475MeV_In.pdf";
        Char_t* RatioPNG = "./Ratio_475MeV_In.png";
        RebinVal = 2;
      }

    if(i==8){
        Char_t* Title = "Ratio of Angle bin to maximum for 525 +/- 25 MeV; Ratio525MeV_In";
        TH1D *hist = (TH1D*)f->Get("Phi_Scattered_525MeV_In");
        Char_t* RatioPDF = "./Ratio_525MeV_In.pdf";
        Char_t* RatioPNG = "./Ratio_525MeV_In.png";
        RebinVal = 2;
      }

    if(i==9){
        Char_t* Title = "Ratio of Angle bin to maximum for 575 +/- 25 MeV; Ratio575MeV_In";
        TH1D *hist = (TH1D*)f->Get("Phi_Scattered_575MeV_In");
        Char_t* RatioPDF = "./Ratio_575MeV_In.pdf";
        Char_t* RatioPNG = "./Ratio_575MeV_In.png";
        RebinVal = 2;
      }

    BinWidth = RebinVal*10; // Default bin size is 10 degrees so x by 10

    TCanvas *canvas = new TCanvas("canvas","canvas",1000,10,550,400); // Create new canvas for ratio plot

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
    Float_t yMax = 2;

    strcpy(hrTitle, Title);
    hr = canvas->DrawFrame(xMin,yMin,xMax,yMax);
    hr->SetTitle(hrTitle);

    hist->Rebin(RebinVal);
    nBins = hist->GetSize() - 2; // -2 as otherwise under/overflow included
    LowMax = hist->GetBinContent(1); // Get value of leftmost bin
    HighMax = hist->GetBinContent(nBins); // Get value of final bin
    AvgMax = (LowMax + HighMax)/2; // Avg out maxima

    for (Int_t j = 0; j < (nBins); j++){

    ratio[i][j] = hist->GetBinContent(j+1)/AvgMax;

    }

    TGraph *gr;

    if (RebinVal == 1) {

        double x1[36];
        double y1[36];

        for (Int_t k = 0; k < 36; k++){

            x1[k] = ((-180 + (BinWidth/2)) + (k*BinWidth));
            y1[k] = ratio[i][k];
            //cout << x1[k] << "    "  << y1[k] << endl;

        }

        gr  = new TGraph(nBins, x1, y1);

    }

    if (RebinVal == 2) {

        double x2[18];
        double y2[18];

        for (Int_t k = 0; k < 18; k++){

            x2[k] = ((-180 + (BinWidth/2)) + (k*BinWidth));
            y2[k] = ratio[i][k];
            //cout << x2[k] << "    "  << y2[k] << endl;

        }

        gr  = new TGraph(nBins, x2, y2);

    }

    gr->SetLineColor(4);
    gr->SetMarkerStyle(2);
    gr->SetMarkerColor(4);
    gr->SetMarkerSize(1.5);
    gr->Draw("SAMEP");
    gr->Fit("pol2");
    pol2->SetLineColor(2);
    pol2->Draw("SAMES");
    gPad->Update();

    double Param1 = pol2->GetParameter(0);
    double Param1Error = pol2->GetParError(0);
    double Param2 = pol2->GetParameter(1);
    double Param2Error = pol2->GetParError(1);
    double Param3 = pol2->GetParameter(2);
    double Param3Error = pol2->GetParError(2);

    P1[i] = Param1;
    P1Err[i] = Param1Error;
    P2[i] = Param2;
    P2Err[i] = Param2Error;
    P3[i] = Param3;
    P3Err[i]= Param3Error;

    canvas->SaveAs(filename = RatioPDF);
    canvas->SaveAs(filename = RatioPNG);

  }

    TFile f1("MCInParameters.root", "RECREATE");

    TTree* tree = new TTree("Parameter_Values", "Tree_of_Values");

    tree->Branch("Par1", &Par1, "Par1/D");
    tree->Branch("Par1Err", &Par1Err, "Par1Err/D");
    tree->Branch("Par2", &Par2, "Par2/D");
    tree->Branch("Par2Err", &Par2Err, "Par2Err/D");
    tree->Branch("Par3", &Par3, "Par3/D");
    tree->Branch("Par3Err", &Par3Err, "Par3Err/D");

    for (Int_t m = 0; m < 10; m++){

        Par1 = P1[m];
        Par1Err = P1Err[m];
        Par2 = P2[m];
        Par2Err = P2Err[m];
        Par3 = P3[m];
        Par3Err = P3Err[m];
        tree->Fill();

    }

    f1.Write();

}
