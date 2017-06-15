#include "./includes_Sigma.h"

void Sigma(){
TFile *MBData = TFile::Open("/scratch/Mainz_Software/a2GoAT/Sig_res_St.root");

    TFile *f1= TFile::Open("/scratch/Mainz_Software/a2GoAT/ParaPerpAsymm_Total_16.root");
    TTree *t1 = (TTree*)f1->Get("Parameter_Values");

    Double_t pValues435[10], pValues455[10], pValues475[10], pValues495[10], pValues515[10], pValues535[10], pValues555[10], pValues575[10],, pValues595[10], pValues615[10];
    Double_t pErrValues435[10], pErrValues455[10], pErrValues475[10],pErrValues495[10], pErrValues515[10], pErrValues535[10], pErrValues555[10], pErrValues575[10], pErrValues585[10], pErrValues595[10], pErrValues615[10];
    Double_t pSigmaValues435[10], pSigmaValues455[10], pSigmaValues475[10],, pSigmaValues495[10], pSigmaValues515[10], pSigmaValues535[10], pSigmaValues555[10], pSigmaValues575[10], pSigmaValues595[10], pSigmaValues615[10];
    Double_t pSigmaErrValues435[10], pSigmaErrValues455[10], pSigmaErrValues475[10], pSigmaErrValues495[10], pSigmaErrValues515[10], pSigmaErrValues535[10], pSigmaErrValues555[10], pSigmaErrValues575[10], pSigmaErrValues595[10], pSigmaErrValues615[10];;

    Double_t pCosAmp435, pCosAmpErr435;
    Double_t pCosAmp455, pCosAmpErr455;
    Double_t pCosAmp475, pCosAmpErr475;
    Double_t pCosAmp495, pCosAmpErr495;
    Double_t pCosAmp515, pCosAmpErr515;
    Double_t pCosAmp535, pCosAmpErr535;
    Double_t pCosAmp555, pCosAmpErr555;
    Double_t pCosAmp575, pCosAmpErr575;
    Double_t pCosAmp595, pCosAmpErr595;
    Double_t pCosAmp615, pCosAmpErr615;

    // Set branch addresses to get values from
    t1->SetBranchAddress("pCosAmp435", &pCosAmp435);
    t1->SetBranchAddress("pCosAmpErr435", &pCosAmpErr435);
    t1->SetBranchAddress("pCosAmp455", &pCosAmp455);
    t1->SetBranchAddress("pCosAmpErr455", &pCosAmpErr455);
    t1->SetBranchAddress("pCosAmp475", &pCosAmp475);
    t1->SetBranchAddress("pCosAmpErr475", &pCosAmpErr475);
    t1->SetBranchAddress("pCosAmp495", &pCosAmp495);
    t1->SetBranchAddress("pCosAmpErr495", &pCosAmpErr495);
    t1->SetBranchAddress("pCosAmp515", &pCosAmp515);
    t1->SetBranchAddress("pCosAmpErr515", &pCosAmpErr515);
    t1->SetBranchAddress("pCosAmp535", &pCosAmp535);
    t1->SetBranchAddress("pCosAmpErr535", &pCosAmpErr535);
    t1->SetBranchAddress("pCosAmp555", &pCosAmp555);
    t1->SetBranchAddress("pCosAmpErr555", &pCosAmpErr555);
    t1->SetBranchAddress("pCosAmp575", &pCosAmp575);
    t1->SetBranchAddress("pCosAmpErr575", &pCosAmpErr575);
    t1->SetBranchAddress("pCosAmp595", &pCosAmp595);
    t1->SetBranchAddress("pCosAmpErr595", &pCosAmpErr595);
    t1->SetBranchAddress("pCosAmp615", &pCosAmp615);
    t1->SetBranchAddress("pCosAmpErr615", &pCosAmpErr615);

    // Load values from tree and asign values back into an array
    for (Int_t k = 0; k < 10; k++){
        Parameter_Values->GetEntry(k);
        pValues435[k] = pCosAmp435;
        pErrValues435[k] = pCosAmpErr435;
        pValues455[k] = pCosAmp455;
        pErrValues455[k] = pCosAmpErr455;
        pValues475[k] = pCosAmp475;
        pErrValues475[k] = pCosAmpErr475;
        pValues495[k] = pCosAmp495;
        pErrValues495[k] = pCosAmpErr495;
        pValues515[k] = pCosAmp515;
        pErrValues515[k] = pCosAmpErr515;
        pValues535[k] = pCosAmp535;
        pErrValues535[k] = pCosAmpErr535;
        pValues555[k] = pCosAmp555;
        pErrValues555[k] = pCosAmpErr555;
        pValues575[k] = pCosAmp575;
        pErrValues575[k] = pCosAmpErr575;
        pValues595[k] = pCosAmp595;
        pErrValues595[k] = pCosAmpErr595;
        pValues615[k] = pCosAmp615;
        pErrValues615[k] = pCosAmpErr615;
    }

    TFile *f2= TFile::Open("/scratch/Mainz_Software/a2GoAT/LinPol_Aug16.root"); // Open linear polarisation plot

    // Calculate values of sigma for each angular and energy bin
    for (Int_t i = 0; i < 10; i++){

        pSigmaValues435[i] = pValues435[i]/(Graph->Eval(435,0));
        pSigmaErrValues435[i] = pErrValues435[i]/(Graph->Eval(435,0));
        pSigmaValues455[i] = pValues455[i]/(Graph->Eval(455,0));
        pSigmaErrValues455[i] = pErrValues455[i]/(Graph->Eval(455,0));
        pSigmaValues475[i] = pValues475[i]/(Graph->Eval(475,0));
        pSigmaErrValues475[i] = pErrValues475[i]/(Graph->Eval(475,0));
        pSigmaValues495[i] = pValues495[i]/(Graph->Eval(495,0));
        pSigmaErrValues495[i] = pErrValues495[i]/(Graph->Eval(495,0));
        pSigmaValues515[i] = pValues515[i]/(Graph->Eval(515,0));
        pSigmaErrValues515[i] = pErrValues515[i]/(Graph->Eval(515,0));
        pSigmaValues535[i] = pValues535[i]/(Graph->Eval(535,0));
        pSigmaErrValues535[i] = pErrValues535[i]/(Graph->Eval(535,0));
        pSigmaValues555[i] = pValues555[i]/(Graph->Eval(555,0));
        pSigmaErrValues555[i] = pErrValues555[i]/(Graph->Eval(555,0));
        pSigmaValues575[i] = pValues575[i]/(Graph->Eval(575,0));
        pSigmaErrValues575[i] = pErrValues575[i]/(Graph->Eval(575,0));
        pSigmaValues595[i] = pValues595[i]/(Graph->Eval(595,0));
        pSigmaErrValues595[i] = pErrValues595[i]/(Graph->Eval(595,0));
        pSigmaValues615[i] = pValues615[i]/(Graph->Eval(615,0));
        pSigmaErrValues615[i] = pErrValues615[i]/(Graph->Eval(615,0));

    }

    TFile f3("Sigma_Plots_16.root", "RECREATE");

    Float_t xMin = -1;
    Float_t xMax = 1;
    Float_t yMin = -5;
    Float_t yMax = 5;
    Double_t x[10] = {0.9, 0.7, 0.5, 0.3, 0.1, -0.1, -0.3, -0.5, -0.7, -0.9}; // Need to adjust
    Double_t ex[10] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1}; // Need to adjust

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
    hr1->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 435 #pm 10 MeV)");

    gr1 = new TGraphErrors(10, x, pSigmaValues435 , ex, pSigmaErrValues435);
    gr1->SetMarkerColor(2);
    gr1->SetMarkerStyle(5);
    gr1->SetMarkerSize(2);
    gr1->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 435 #pm 10 MeV)");
    gr1->GetXaxis()->SetTitle("Cos#theta_{CM}");
    gr1->GetYaxis()->SetTitle("#Sigma");
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
    hr2->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 455 #pm 10 MeV)");

    gr2 = new TGraphErrors(10, x, pSigmaValues455 , ex, pSigmaErrValues455);
    gr2->SetMarkerColor(2);
    gr2->SetMarkerStyle(5);
    gr2->SetMarkerSize(2);
    gr2->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 455 #pm 10 MeV)");
    gr2->GetXaxis()->SetTitle("Cos#theta_{CM}");
    gr2->GetYaxis()->SetTitle("#Sigma");
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
    hr3->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 475 #pm 10 MeV)");

    gr3 = new TGraphErrors(10, x, pSigmaValues475, ex, pSigmaErrValues475);
    gr3->SetMarkerColor(2);
    gr3->SetMarkerStyle(5);
    gr3->SetMarkerSize(2);
    gr3->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 475 #pm 10 MeV)");
    gr3->GetXaxis()->SetTitle("Cos#theta_{CM}");
    gr3->GetYaxis()->SetTitle("#Sigma");
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
    hr4->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 495 #pm 10 MeV)");

    gr4 = new TGraphErrors(10, x, pSigmaValues495, ex, pSigmaErrValues495);
    gr4->SetMarkerColor(2);
    gr4->SetMarkerStyle(5);
    gr4->SetMarkerSize(2);
    gr4->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 495 #pm 10 MeV)");
    gr4->GetXaxis()->SetTitle("Cos#theta_{CM}");
    gr4->GetYaxis()->SetTitle("#Sigma");
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
    hr5->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 515 #pm 10 MeV)");

    gr5 = new TGraphErrors(10, x, pSigmaValues515, ex, pSigmaErrValues515);
    gr5->SetMarkerColor(2);
    gr5->SetMarkerStyle(5);
    gr5->SetMarkerSize(2);
    gr5->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 515 #pm 10 MeV)");
    gr5->GetXaxis()->SetTitle("Cos#theta_{CM}");
    gr5->GetYaxis()->SetTitle("#Sigma");
    gr5->Draw("ep");

    TCanvas *canvas6 = new TCanvas("canvas6","canvas6", 1920, 1080);
    TPad *pad6 = new TPad("pad6","",0,0,1,1);
    pad6->Draw();
    pad->cd();

    pad6->SetTickx(1);
    pad6->SetTicky(1);
    pad6->SetGridx(1);
    pad6->SetGridy(1);
    TH1F  *hr6;
    hr6 = canvas6->DrawFrame(xMin,-1,xMax,1);
    hr6->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 535 #pm 10 MeV)");

    gr6 = new TGraphErrors(10, x, pSigmaValues535, ex, pSigmaErrValues535);
    gr6->SetMarkerColor(2);
    gr6->SetMarkerStyle(5);
    gr6->SetMarkerSize(2);
    gr6->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 535 #pm 10 MeV)");
    gr6->GetXaxis()->SetTitle("Cos#theta_{CM}");
    gr6->GetYaxis()->SetTitle("#Sigma");
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
    hr7->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 555 #pm 10 MeV)");

    gr7 = new TGraphErrors(10, x, pSigmaValues555, ex, pSigmaErrValues555);
    gr7->SetMarkerColor(2);
    gr7->SetMarkerStyle(5);
    gr7->SetMarkerSize(2);
    gr7->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 555 #pm 10 MeV)");
    gr7->GetXaxis()->SetTitle("Cos#theta_{CM}");
    gr7->GetYaxis()->SetTitle("#Sigma");
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
    hr8->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 575 #pm 10 MeV)");

    gr8 = new TGraphErrors(10, x, pSigmaValues575, ex, pSigmaErrValues575);
    gr8->SetMarkerColor(2);
    gr8->SetMarkerStyle(5);
    gr8->SetMarkerSize(2);
    gr8->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 575 #pm 10 MeV)");
    gr8->GetXaxis()->SetTitle("Cos#theta_{CM}");
    gr8->GetYaxis()->SetTitle("#Sigma");
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
    hr9->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 595 #pm 10 MeV)");

    gr9 = new TGraphErrors(10, x, pSigmaValues595, ex, pSigmaErrValues595);
    gr9->SetMarkerColor(2);
    gr9->SetMarkerStyle(5);
    gr9->SetMarkerSize(2);
    gr9->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 595 #pm 10 MeV)");
    gr9->GetXaxis()->SetTitle("Cos#theta_{CM}");
    gr9->GetYaxis()->SetTitle("#Sigma");
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
    hr10->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 615 #pm 10 MeV)");

    gr10 = new TGraphErrors(10, x, pSigmaValues615, ex, pSigmaErrValues615);
    gr10->SetMarkerColor(2);
    gr10->SetMarkerStyle(5);
    gr10->SetMarkerSize(2);
    gr10->SetTitle("#Sigma as fn of Cos#theta_{CM} (E_{#gamma} 615 #pm 10 MeV)");
    gr10->GetXaxis()->SetTitle("Cos#theta_{CM}");
    gr10->GetYaxis()->SetTitle("#Sigma");
    gr10->Draw("ep");

    canvas1->Write();
    canvas2->Write();
    canvas3->Write();
    canvas4->Write();
    canvas5->Write();
    canvas6->Write();
    canvas7->Write();
    canvas8->Write();
    canvas9->Write();

    TCanvas *canvas11 = new TCanvas("canvas11","canvas11", 1920, 1080);
    canvas11->Divide(5,2);
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

    canvas11->Write();
    f3.Write();
}
