#include "./includes_Sigma.h"

// define a legendre polynomial, need to define some limits?
Double_t legendre(Double_t *x,Double_t *par)
{
    Double_t fitval = 0;
    fitval = (1-TMath::Power(x[0],2)*(par[0]*3+par[1]*15*x[0]+par[2]*15.0/2*(7*TMath::Power(x[0],2)-1)+par[3]*105.0/2*x[0]*(3*TMath::Power(x[0],2)-1)+par[4]*105.0/8*(33*TMath::Power(x[0],4)-18*TMath::Power(x[0],2)+1)+par[5]*63.0/8*x*(143*TMath::Power(x[0],4)-110*TMath::Power(x[0],2)+15)+par[6]*315.0/16*(143*TMath::Power(x[0],6)-143*TMath::Power(x[0],4)+33*TMath::Power(x[0],2)-1)));
    return fitval;
}

void Sigma(){
    TFile *MBData = TFile::Open("/scratch/Mainz_Software/a2GoAT/Sig_res_St.root");

    TFile *f1= TFile::Open("/scratch/Mainz_Software/a2GoAT/ParaPerpAsymm_Total_20.root");
    TTree *t1 = (TTree*)f1->Get("Parameter_Values");
    TF1 *LegendreFunc = new TF1("LegendreFit",  legendre, -1, 1, 7); //Give a name and range to the fitting funcion

    Double_t pValues430[5], pValues450[5], pValues470[5], pValues490[5], pValues510[5], pValues530[5], pValues550[5], pValues570[5],, pValues590[5], pValues610[5];
    Double_t pErrValues430[5], pErrValues450[5], pErrValues470[5],pErrValues490[5], pErrValues510[5], pErrValues530[5], pErrValues550[5], pErrValues570[5], pErrValues580[5], pErrValues590[5], pErrValues610[5];
    Double_t pSigmaValues430[5], pSigmaValues450[5], pSigmaValues470[5],, pSigmaValues490[5], pSigmaValues510[5], pSigmaValues530[5], pSigmaValues550[5], pSigmaValues570[5], pSigmaValues590[5], pSigmaValues610[5];
    Double_t pSigmaErrValues430[5], pSigmaErrValues450[5], pSigmaErrValues470[5], pSigmaErrValues490[5], pSigmaErrValues510[5], pSigmaErrValues530[5], pSigmaErrValues550[5], pSigmaErrValues570[5], pSigmaErrValues590[5], pSigmaErrValues610[5];;

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

    // Set branch addresses to get values from
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

    // Load values from tree and asign values back into an array
    for (Int_t k = 0; k < 5; k++){
        Parameter_Values->GetEntry(k);
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
    }

    TFile *f2= TFile::Open("/scratch/Mainz_Software/a2GoAT/LinPol_Aug16.root"); // Open linear polarisation plot

    // Calculate values of sigma for each angular and energy bin
    for (Int_t i = 0; i < 5; i++){

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

    }

    TFile f3("Sigma_Plots_20.root", "RECREATE");

    Float_t xMin = -1;
    Float_t xMax = 1;
    Float_t yMin = -5;
    Float_t yMax = 5;
    Double_t x[5] = {0.8, 0.4, 0.0, -0.4, -0.8};
    Double_t ex[5] = {0.2, 0.2, 0.2, 0.2, 0.2};

    TCanvas *canvas1 = new TCanvas("canvas1","canvas1", 1920, 1080);
    TPad *pad1 = new TPad("pad1","",0,0,1,1);
    pad1->Draw();
    pad1->cd();

    pad1->SetTickx(1);
    pad1->SetTicky(1);
    pad1->SetGridx(1);
    pad1->SetGridy(1);
    TH1F  *hr1;
    hr1 = canvas1->DrawFrame(xMin,-1 ,xMax, 1 );
    hr1->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 430 #pm 10 MeV)");

    gr1 = new TGraphErrors(5, x, pSigmaValues430 , ex, pSigmaErrValues430);
    gr1->Fit("LegendreFit");
    gr1->SetMarkerColor(2);
    gr1->SetLineColor(2);
    gr1->SetMarkerStyle(5);
    gr1->SetMarkerSize(2);
    gr1->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 430 #pm 10 MeV)");
    gr1->GetXaxis()->SetTitle("Cos#theta_{CM}");
    gr1->GetYaxis()->SetTitle("#Sigma");
    gr1->SetName("Sigma430_S");
    gr1->Draw("ep");

    TCanvas *canvas2 = new TCanvas("canvas2","canvas2", 1920, 1080);
    TPad *pad2 = new TPad("pad2","",0,0,1,1);
    pad2->Draw();
    pad2->cd();

    pad2->SetTickx(1);
    pad2->SetTicky(1);
    pad2->SetGridx(1);
    pad2->SetGridy(1);
    TH1F  *hr2;
    hr2 = canvas2->DrawFrame(xMin,-1,xMax,1);
    hr2->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 450 #pm 10 MeV)");

    gr2 = new TGraphErrors(5, x, pSigmaValues450 , ex, pSigmaErrValues450);
    gr2->Fit("LegendreFit");
    gr2->SetMarkerColor(2);
    gr2->SetLineColor(2);
    gr2->SetMarkerStyle(5);
    gr2->SetMarkerSize(2);
    gr2->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 450 #pm 10 MeV)");
    gr2->GetXaxis()->SetTitle("Cos#theta_{CM}");
    gr2->GetYaxis()->SetTitle("#Sigma");
    gr2->SetName("Sigma450_S");
    gr2->Draw("ep");

    TCanvas *canvas3 = new TCanvas("canvas3","canvas3", 1920, 1080);
    TPad *pad3 = new TPad("pad3","",0,0,1,1);
    pad3->Draw();
    pad3->cd();

    pad3->SetTickx(1);
    pad3->SetTicky(1);
    pad3->SetGridx(1);
    pad3->SetGridy(1);
    TH1F  *hr3;
    hr3 = canvas3->DrawFrame(xMin,-1,xMax,1);
    hr3->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 470 #pm 10 MeV)");

    gr3 = new TGraphErrors(5, x, pSigmaValues470, ex, pSigmaErrValues470);
    gr3->Fit("LegendreFit");
    gr3->SetMarkerColor(2);
    gr3->SetLineColor(2);
    gr3->SetMarkerStyle(5);
    gr3->SetMarkerSize(2);
    gr3->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 470 #pm 10 MeV)");
    gr3->GetXaxis()->SetTitle("Cos#theta_{CM}");
    gr3->GetYaxis()->SetTitle("#Sigma");
    gr3->SetName("Sigma470_S");
    gr3->Draw("ep");


    TCanvas *canvas4 = new TCanvas("canvas4","canvas4", 1920, 1080);
    TPad *pad4 = new TPad("pad4","",0,0,1,1);
    pad4->Draw();
    pad4->cd();

    pad4->SetTickx(1);
    pad4->SetTicky(1);
    pad4->SetGridx(1);
    pad4->SetGridy(1);
    TH1F  *hr4;
    hr4 = canvas4->DrawFrame(xMin,-1,xMax,1);
    hr4->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 490 #pm 10 MeV)");

    gr4 = new TGraphErrors(5, x, pSigmaValues490, ex, pSigmaErrValues490);
    gr4->Fit("LegendreFit");
    gr4->SetMarkerColor(2);
    gr4->SetLineColor(2);
    gr4->SetMarkerStyle(5);
    gr4->SetMarkerSize(2);
    gr4->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 490 #pm 10 MeV)");
    gr4->GetXaxis()->SetTitle("Cos#theta_{CM}");
    gr4->GetYaxis()->SetTitle("#Sigma");
    gr4->SetName("Sigma490_S");
    gr4->Draw("ep");

    TCanvas *canvas5 = new TCanvas("canvas5","canvas5", 1920, 1080);
    TPad *pad5 = new TPad("pad5","",0,0,1,1);
    pad5->Draw();
    pad5->cd();

    pad5->SetTickx(1);
    pad5->SetTicky(1);
    pad5->SetGridx(1);
    pad5->SetGridy(1);
    TH1F  *hr5;
    hr5 = canvas5->DrawFrame(xMin,-1,xMax,1);
    hr5->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 510 #pm 10 MeV)");

    gr5 = new TGraphErrors(5, x, pSigmaValues510, ex, pSigmaErrValues510);
    gr5->Fit("LegendreFit");
    gr5->SetMarkerColor(2);
    gr5->SetLineColor(2);
    gr5->SetMarkerStyle(5);
    gr5->SetMarkerSize(2);
    gr5->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 510 #pm 10 MeV)");
    gr5->GetXaxis()->SetTitle("Cos#theta_{CM}");
    gr5->GetYaxis()->SetTitle("#Sigma");
    gr5->SetName("Sigma510_S");
    gr5->Draw("ep");

    TCanvas *canvas6 = new TCanvas("canvas6","canvas6", 1920, 1080);
    TPad *pad6 = new TPad("pad6","",0,0,1,1);
    pad6->Draw();
    pad6->cd();

    pad6->SetTickx(1);
    pad6->SetTicky(1);
    pad6->SetGridx(1);
    pad6->SetGridy(1);
    TH1F  *hr6;
    hr6 = canvas6->DrawFrame(xMin,-1,xMax,1);
    hr6->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 530 #pm 10 MeV)");

    gr6 = new TGraphErrors(5, x, pSigmaValues530, ex, pSigmaErrValues530);
    gr6->Fit("LegendreFit");
    gr6->SetMarkerColor(2);
    gr6->SetLineColor(2);
    gr6->SetMarkerStyle(5);
    gr6->SetMarkerSize(2);
    gr6->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 530 #pm 10 MeV)");
    gr6->GetXaxis()->SetTitle("Cos#theta_{CM}");
    gr6->GetYaxis()->SetTitle("#Sigma");
    gr6->SetName("Sigma530_S");
    gr6->Draw("ep");

    TCanvas *canvas7 = new TCanvas("canvas7","canvas7", 1920, 1080);
    TPad *pad7 = new TPad("pad7","",0,0,1,1);
    pad7->Draw();
    pad7->cd();

    pad7->SetTickx(1);
    pad7->SetTicky(1);
    pad7->SetGridx(1);
    pad7->SetGridy(1);
    TH1F  *hr7;
    hr7 = canvas7->DrawFrame(xMin,-1,xMax,1);
    hr7->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 550 #pm 10 MeV)");

    gr7 = new TGraphErrors(5, x, pSigmaValues550, ex, pSigmaErrValues550);
    gr7->Fit("LegendreFit");
    gr7->SetMarkerColor(2);
    gr7->SetLineColor(2);
    gr7->SetMarkerStyle(5);
    gr7->SetMarkerSize(2);
    gr7->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 550 #pm 10 MeV)");
    gr7->GetXaxis()->SetTitle("Cos#theta_{CM}");
    gr7->GetYaxis()->SetTitle("#Sigma");
    gr7->SetName("Sigma550_S");
    gr7->Draw("ep");

    TCanvas *canvas8 = new TCanvas("canvas8","canvas8", 1920, 1080);
    TPad *pad8 = new TPad("pad8","",0,0,1,1);
    pad8->Draw();
    pad8->cd();

    pad8->SetTickx(1);
    pad8->SetTicky(1);
    pad8->SetGridx(1);
    pad8->SetGridy(1);
    TH1F  *hr8;
    hr8 = canvas8->DrawFrame(xMin,-1,xMax,1);
    hr8->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 570 #pm 10 MeV)");

    gr8 = new TGraphErrors(5, x, pSigmaValues570, ex, pSigmaErrValues570);
    gr8->Fit("LegendreFit");
    gr8->SetMarkerColor(2);
    gr8->SetLineColor(2);
    gr8->SetMarkerStyle(5);
    gr8->SetMarkerSize(2);
    gr8->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 570 #pm 10 MeV)");
    gr8->GetXaxis()->SetTitle("Cos#theta_{CM}");
    gr8->GetYaxis()->SetTitle("#Sigma");
    gr8->SetName("Sigma570_S");
    gr8->Draw("ep");

    TCanvas *canvas9 = new TCanvas("canvas9","canvas9", 1920, 1080);
    TPad *pad9 = new TPad("pad9","",0,0,1,1);
    pad9->Draw();
    pad9->cd();

    pad9->SetTickx(1);
    pad9->SetTicky(1);
    pad9->SetGridx(1);
    pad9->SetGridy(1);
    TH1F  *hr9;
    hr9 = canvas9->DrawFrame(xMin,-1,xMax,1);
    hr9->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 590 #pm 10 MeV)");

    gr9 = new TGraphErrors(5, x, pSigmaValues590, ex, pSigmaErrValues590);
    gr9->Fit("LegendreFit");
    gr9->SetMarkerColor(2);
    gr9->SetLineColor(2);
    gr9->SetMarkerStyle(5);
    gr9->SetMarkerSize(2);
    gr9->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 590 #pm 10 MeV)");
    gr9->GetXaxis()->SetTitle("Cos#theta_{CM}");
    gr9->GetYaxis()->SetTitle("#Sigma");
    gr9->SetName("Sigma590_S");
    gr9->Draw("ep");

    TCanvas *canvas10 = new TCanvas("canvas10","canvas10", 1920, 1080);
    TPad *pad10 = new TPad("pad10","",0,0,1,1);
    pad10->Draw();
    pad10->cd();

    pad10->SetTickx(1);
    pad10->SetTicky(1);
    pad10->SetGridx(1);
    pad10->SetGridy(1);
    TH1F  *hr10;
    hr10 = canvas10->DrawFrame(xMin,-1,xMax,1);
    hr10->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 610 #pm 10 MeV)");

    gr10 = new TGraphErrors(5, x, pSigmaValues610, ex, pSigmaErrValues610);
    gr10->Fit("LegendreFit");
    gr10->SetMarkerColor(2);
    gr10->SetLineColor(2);
    gr10->SetMarkerStyle(5);
    gr10->SetMarkerSize(2);
    gr10->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 610 #pm 10 MeV)");
    gr10->GetXaxis()->SetTitle("Cos#theta_{CM}");
    gr10->GetYaxis()->SetTitle("#Sigma");
    gr10->SetName("Sigma610_S");
    gr10->Draw("ep");

    canvas1->Write();
    gr1->Write();
    canvas2->Write();
    gr2->Write();
    canvas3->Write();
    gr3->Write();
    canvas4->Write();
    gr4->Write();
    canvas5->Write();
    gr5->Write();
    canvas6->Write();
    gr6->Write();
    canvas7->Write();
    gr7->Write();
    canvas8->Write();
    gr8->Write();
    canvas9->Write();
    gr9->Write();
    canvas10->Write();
    gr10->Write();

    TCanvas *canvas11 = new TCanvas("canvas11","canvas11", 1920, 1080);
    canvas11->Divide(4,3);
    canvas11->cd(1);
    pad1->Draw();
    canvas11->cd(2);
    pad2->Draw();
    canvas11->cd(3);
    pad3->Draw();
    canvas11->cd(4);
    pad4->Draw();
    canvas11->cd(5);
    pad5->Draw();
    canvas11->cd(6);
    pad6->Draw();
    canvas11->cd(7);
    pad7->Draw();
    canvas11->cd(8);
    pad8->Draw();
    canvas11->cd(9);
    pad9->Draw();
    canvas11->cd(10);
    pad10->Draw();

    canvas11->Write();
    f3.Write();
}
