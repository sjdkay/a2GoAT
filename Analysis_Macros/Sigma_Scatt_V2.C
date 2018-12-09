#include "./includes_Sigma_Scatt_V2.h"

void Sigma_Scatt_V2(){

    TFile *MBData = TFile::Open("/scratch/Mainz_Software/a2GoAT/Sig_res_St.root");
    TGraphErrors* SigmaPlots[10];

    char name[21];
    char title[60];
    char MBname[20];
    char name2[60];
    char title2[60];

    TH1F* MBHist[20];

    for (Int_t i = 0; i < 20; i++){
        sprintf(MBname, "hslC%i", 43+i);
        MBHist[i] = (TH1F*)MBData->Get(MBname);
    }

    TFile *f1= TFile::Open("/scratch/Mainz_Software/a2GoAT/Results/ParaPerpAsymm_S35.root");
    TTree *t1 = (TTree*)f1->Get("Parameter_Values");

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


    // Load values from tree and asign values back into an array
    for (Int_t k = 0; k < 5; k++){
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

    }

    double x412[7] = { cos(35.0*TMath::DegToRad()), cos(55.0*TMath::DegToRad()), cos(75.0*TMath::DegToRad()), cos(95.0*TMath::DegToRad()), cos(115.0*TMath::DegToRad()), cos(135.0*TMath::DegToRad()), cos(155.0*TMath::DegToRad())};
    double ex412[7] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
    double y412[7] = {0.353, -0.124, -0.188, -0.474, -0.355, -0.301, 0.003};
    double ey412[7] = {sqrt(0.114**2+0.033**2), sqrt(0.072**2+0.014**2), sqrt(0.072**2+0.02**2), sqrt(0.061**2+0.049**2), sqrt(0.058**2+0.032**2), sqrt(0.073**2+0.027**2), sqrt(0.103**2+0.015**2)};
    e412 =  new TGraphErrors(7, x412, y412, ex412, ey412);
    e412->SetMarkerStyle(24);
    e412->SetMarkerColor(1);
    e412->SetMarkerSize(1.5);
    e412->SetLineColor(1);
    e412->SetLineWidth(3);

    double x435[2] = {cos(75.0*TMath::DegToRad()), cos(90.0*TMath::DegToRad())};
    double ex435[2] = {0.1, 0.1};
    double y435[2] = {-0.2, -0.23};
    double ey435[2] = {0.05, 0.04};
    e435 = new TGraphErrors(2, x435, y435, ex435, ey435);
    e435->SetMarkerStyle(24);
    e435->SetMarkerColor(1);
    e435->SetMarkerSize(1.5);
    e435->SetLineColor(1);
    e435->SetLineWidth(3);

    double x455[2] = {cos(75.0*TMath::DegToRad()), cos(90.0*TMath::DegToRad())};
    double ex455[2] = {0.1, 0.1};
    double y455[2] = {-0.25, -0.23};
    double ey455[2] = {0.04, 0.04};
    e455 = new TGraphErrors(2, x455, y455, ex455, ey455);
    e455->SetMarkerStyle(24);
    e455->SetMarkerColor(1);
    e455->SetMarkerSize(1.5);
    e455->SetLineColor(1);
    e455->SetLineWidth(3);

    double x465[2] = {cos(45.0*TMath::DegToRad()), cos(60.0*TMath::DegToRad())};
    double ex465[2] = {0.1, 0.1};
    double y465[2] = {0.04, -0.04};
    double ey465[2] = {0.04, 0.04};
    e465 = new TGraphErrors(2, x465, y465, ex465, ey465);
    e465->SetMarkerStyle(24);
    e465->SetMarkerColor(1);
    e465->SetMarkerSize(1.5);
    e465->SetLineColor(1);
    e465->SetLineWidth(3);

    double x475[3] = {cos(45.0*TMath::DegToRad()), cos(75.0*TMath::DegToRad()), cos(90.0*TMath::DegToRad())};
    double ex475[3] = {0.1, 0.1, 0.1};
    double y475[3] = {-0.01, -0.15, -0.09};
    double ey475[3] = {0.05, 0.04, 0.04};
    e475 = new TGraphErrors(3, x475, y475, ex475, ey475);
    e475->SetMarkerStyle(24);
    e475->SetMarkerColor(1);
    e475->SetMarkerSize(1.5);
    e475->SetLineColor(1);
    e475->SetLineWidth(3);

    double x485[2] = {cos(60.0*TMath::DegToRad()), cos(75.0*TMath::DegToRad())};
    double ex485[2] = {0.1, 0.1};
    double y485[2] = {0.02, -0.22};
    double ey485[2] = {0.05, 0.06};
    e485 = new TGraphErrors(2, x485, y485, ex485, ey485);
    e485->SetMarkerStyle(24);
    e485->SetMarkerColor(1);
    e485->SetMarkerSize(1.5);
    e485->SetLineColor(1);
    e485->SetLineWidth(3);

    double x505[3] = {cos(45.0*TMath::DegToRad()), cos(60.0*TMath::DegToRad()), cos(90.0*TMath::DegToRad())};
    double ex505[3] = {0.1, 0.1, 0.1};
    double y505[3] = {0.07, 0.02, -0.15};
    double ey505[3] = {0.05, 0.06, 0.05};
    e505 = new TGraphErrors(3, x505, y505, ex505, ey505);
    e505->SetMarkerStyle(24);
    e505->SetMarkerColor(1);
    e505->SetMarkerSize(1.5);
    e505->SetLineColor(1);
    e505->SetLineWidth(3);

    double x515[2] = {cos(45.0*TMath::DegToRad()), cos(75.0*TMath::DegToRad())};
    double ex515[2] = {0.1, 0.1};
    double y515[2] = {0.07, -0.15};
    double ey515[2] = {0.05, 0.05};
    e515 = new TGraphErrors(2, x515, y515, ex515, ey515);
    e515->SetMarkerStyle(24);
    e515->SetMarkerColor(1);
    e515->SetMarkerSize(1.5);
    e515->SetLineColor(1);
    e515->SetLineWidth(3);

    double x525[3] = {cos(45.0*TMath::DegToRad()), cos(60.0*TMath::DegToRad()), cos(90.0*TMath::DegToRad())};
    double ex525[3] = {0.1, 0.1, 0.1};
    double y525[3] = {0.06, 0.05, -0.08};
    double ey525[3] = {0.05, 0.06, 0.05};
    e525 = new TGraphErrors(3, x525, y525, ex525, ey525);
    e525->SetMarkerStyle(24);
    e525->SetMarkerColor(1);
    e525->SetMarkerSize(1.5);
    e525->SetLineColor(1);
    e525->SetLineWidth(3);

    double x545[2] = {cos(45.0*TMath::DegToRad()), cos(60.0*TMath::DegToRad())};
    double ex545[2] = {0.1, 0.1};
    double y545[2] = {0.12, 0.16};
    double ey545[2] = {0.05, 0.05};
    e545 = new TGraphErrors(2, x545, y545, ex545, ey545);
    e545->SetMarkerStyle(24);
    e545->SetMarkerColor(1);
    e545->SetMarkerSize(1.5);
    e545->SetLineColor(1);
    e545->SetLineWidth(3);

    double x575[3] = {cos(45.0*TMath::DegToRad()), cos(75.0*TMath::DegToRad()), cos(90.0*TMath::DegToRad())};
    double ex575[3] = {0.1, 0.1, 0.1};
    double y575[3] = {0.19, 0.12, 0.02};
    double ey575[3] = {0.04, 0.05, 0.06};
    e575 = new TGraphErrors(3, x575, y575, ex575, ey575);
    e575->SetMarkerStyle(24);
    e575->SetMarkerColor(1);
    e575->SetMarkerSize(1.5);
    e575->SetLineColor(1);
    e575->SetLineWidth(3);

    double x595[3] = {cos(60.0*TMath::DegToRad()), cos(75.0*TMath::DegToRad()), cos(90.0*TMath::DegToRad())};
    double ex595[3] = {0.1, 0.1, 0.1};
    double y595[3] = {0.23, 0.16, 0.18};
    double ey595[3] = {0.05, 0.06, 0.07};
    e595 = new TGraphErrors(3, x595, y595, ex595, ey595);
    e595->SetMarkerStyle(24);
    e595->SetMarkerColor(1);
    e595->SetMarkerSize(1.5);
    e595->SetLineColor(1);
    e595->SetLineWidth(3);

    double x615[4] = {cos(45.0*TMath::DegToRad()), cos(60.0*TMath::DegToRad()), cos(75.0*TMath::DegToRad()), cos(90.0*TMath::DegToRad())};
    double ex615[4] = {0.1, 0.1, 0.1, 0.1};
    double y615[4] = {0.23, 0.22, 0.2, 0.15};
    double ey615[4] = {0.05, 0.05, 0.06, 0.07};
    e615 = new TGraphErrors(4, x615, y615, ex615, ey615);
    e615->SetMarkerStyle(24);
    e615->SetMarkerColor(1);
    e615->SetMarkerSize(1.5);
    e615->SetLineColor(1);
    e615->SetLineWidth(3);

    TFile *f2= TFile::Open("/scratch/Mainz_Software/a2GoAT/LinPol_Aug16.root"); // Open linear polarisation plot

    // Calculate values of sigma for each angular and energy bin
    // Exclude bins at edges!

    for (Int_t i = 0; i < 5; i++){ // Need to re-write this in a less terrible way...

        SigmaValues[0][i] = pValues415[i]/(Graph->Eval(415,0));
        SigmaErrValues[0][i] = pErrValues415[i]/(Graph->Eval(415,0));
        SigmaValues[1][i] = pValues435[i]/(Graph->Eval(435,0));
        SigmaErrValues[1][i] = pErrValues435[i]/(Graph->Eval(435,0));
        SigmaValues[2][i] = pValues455[i]/(Graph->Eval(455,0));
        SigmaErrValues[2][i] = pErrValues455[i]/(Graph->Eval(455,0));
        SigmaValues[3][i] = pValues470[i]/(Graph->Eval(470,0));
        SigmaErrValues[3][i] = pErrValues470[i]/(Graph->Eval(470,0));
        SigmaValues[4][i] = pValues490[i]/(Graph->Eval(490,0));
        SigmaErrValues[4][i] = pErrValues490[i]/(Graph->Eval(490,0));
        SigmaValues[5][i] = pValues510[i]/(Graph->Eval(510,0));
        SigmaErrValues[5][i] = pErrValues510[i]/(Graph->Eval(510,0));
        SigmaValues[6][i] = pValues530[i]/(Graph->Eval(530,0));
        SigmaErrValues[6][i] = pErrValues530[i]/(Graph->Eval(530,0));
        SigmaValues[7][i] = pValues550[i]/(Graph->Eval(550,0));
        SigmaErrValues[7][i] = pErrValues550[i]/(Graph->Eval(550,0));
        SigmaValues[8][i] = pValues570[i]/(Graph->Eval(570,0));
        SigmaErrValues[8][i] = pErrValues570[i]/(Graph->Eval(570,0));
        SigmaValues[9][i] = pValues590[i]/(Graph->Eval(590,0));
        SigmaErrValues[9][i] = pErrValues590[i]/(Graph->Eval(590,0));

    }

    TFile f3("Sigma_Plots_S35.root", "RECREATE");

    Float_t xMin = -1;
    Float_t xMax = 1;
    double x[5] = {0.8, 0.4, 0, -0.4, -0.8};
    double ex[5] = {0.2, 0.2, 0.2, 0.2, 0.2};

    for(Int_t i = 0 ; i < 10 ; i++)
    {
        sprintf(name, "Sigma_%i", 415+(i*20));
        sprintf(title, "#Sigma(Cos#theta_{CM}) E_{#gamma} %i #pm 10 MeV", 415+(i*20));
        SigmaPlots[i] = new TGraphErrors(5 , x, SigmaValues[i], ex, SigmaErrValues[i]);
        SigmaPlots[i]->SetMarkerColor(4);
        SigmaPlots[i]->SetLineColor(4);
        SigmaPlots[i]->SetMarkerStyle(8);
        SigmaPlots[i]->SetMarkerSize(1);
        SigmaPlots[i]->GetXaxis()->SetTitle("Cos#theta_{CM}");
        SigmaPlots[i]->GetXaxis()->SetRangeUser(-1, 1);
        SigmaPlots[i]->GetYaxis()->SetRangeUser(-1, 1);
        SigmaPlots[i]->GetYaxis()->SetTitle("#Sigma");
        SigmaPlots[i]->SetName(name);
        SigmaPlots[i]->SetTitle(title);

        SigmaPlots[i]->Write();

    }

    TCanvas *canvas20 = new TCanvas("canvas20","canvas20", 1920, 1080);
    canvas20->Divide(5,2);
    for(int i = 1 ; i < 11 ; i++){
        canvas20->cd(i);
        SigmaPlots[i-1]->Draw("AEP");
    }

    canvas20->Write();

    f3.Write();
}
