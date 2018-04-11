#include "./includes_Sigma_NoScatt.h"

// define a legendre polynomial, 6 parameters. Need to set limits? This is different polynomial?
//Double_t legendre(Double_t *x,Double_t *par)
//{
//    Double_t fitval = 0;
//    fitval = (1-TMath::Power(x[0],2)*(par[0]*3+par[1]*15*x[0]+par[2]*15.0/2*(7*TMath::Power(x[0],2)-1)+par[3]*105.0/2*x[0]*(3*TMath::Power(x[0],2)-1)+par[4]*105.0/8*(33*TMath::Power(x[0],4)-18*TMath::Power(x[0],2)+1)+par[5]*63.0/8*x*(143*TMath::Power(x[0],4)-110*TMath::Power(x[0],2)+15)+par[6]*315.0/16*(143*TMath::Power(x[0],6)-143*TMath::Power(x[0],4)+33*TMath::Power(x[0],2)-1)));
//    return fitval;
//}

void Sigma_NoScatt_V2(){

    TFile *MBData = TFile::Open("/scratch/Mainz_Software/a2GoAT/Sig_res_St.root");
    TGraphErrors* SigmaPlots[21];
    char name[21];
    char title[60];

    TFile *f1= TFile::Open("/scratch/Mainz_Software/a2GoAT/Results/ParaPerpAsymm_NS14.root");
    TTree *t1 = (TTree*)f1->Get("Parameter_Values");
    //TF1 *LegPol = new TF1("LegPol","(1-x**2)*([0]*3+[1]*15*x+[2]*15.0/2*(7*x**2-1)+[3]*105.0/2*x*(3*x**2-1)+[4]*105.0/8*(33*x**4-18*x**2+1)+[5]*63.0/8*x*(143*x**4-110*x**2+15)+[6]*315.0/16*(143*x**6-143*x**4+33*x**2-1))",-1,1);
    //TF1* LegFunc = new TF1("LegPol", legendre, -1, 1, 7);
    TF1 *LegPol = new TF1("LegPol", "[0]+[1]*x+[2]*(0.5*(3*x**2-1))+[3]*(0.5*(5*x**3-3*x))+[4]*(0.125*(35*x**4-30*x**2+3))+[5]*(1.0/8.0*(63*x**5-70*x**3+15*x))+[6]*(1.0/16*(231*x**6-315*x**4+105*x**2-5))", -1, 1);
    LegPol->SetParLimits(0,-1,1);
    LegPol->SetParLimits(1,-1,1);
    LegPol->SetParLimits(2,-1,1);
    LegPol->SetParLimits(3,-1,1);
    LegPol->SetParLimits(4,-1,1);
    LegPol->SetParLimits(5,-1,1);
    LegPol->SetParLimits(6,-1,1);
    LegPol->SetLineColor(4);
    LegPol->SetLineWidth(3);

    Double_t pValues415[20], pValues425[20], pValues435[20], pValues445[20], pValues455[20], pValues465[20], pValues475[20], pValues485[20], pValues495[20], pValues505[20], pValues515[20], pValues525[20], pValues535[20], pValues545[20], pValues555[20], pValues565[20], pValues575[20], pValues585[20], pValues595[20], pValues605[20], pValues615[20];
    Double_t pErrValues415[20], pErrValues425[20], pErrValues435[20], pErrValues445[20], pErrValues455[20], pErrValues465[20], pErrValues475[20], pErrValues485[20], pErrValues495[20], pErrValues505[20], pErrValues515[20], pErrValues525[20], pErrValues535[20], pErrValues545[20], pErrValues555[20], pErrValues565[20], pErrValues575[20], pErrValues585[20], pErrValues595[20], pErrValues605[20], pErrValues615[20];
    Double_t LegPar[7][21];
    Double_t LegParErr[7][21];
    Double_t SigmaValues[21][18];
    Double_t SigmaErrValues[21][18];
    double p0, p0Err, p1, p1Err, p2, p2Err, p3, p3Err, p4, p4Err, p5, p5Err, p6, p6Err;

    Double_t pCosAmp415, pCosAmpErr415;
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
    t1->SetBranchAddress("pCosAmp415", &pCosAmp415);
    t1->SetBranchAddress("pCosAmpErr415", &pCosAmpErr415);
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
    for (Int_t k = 0; k < 20; k++){
        Parameter_Values->GetEntry(k);
        pValues415[k] = pCosAmp415;
        pErrValues415[k] = pCosAmpErr415;
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
    // Exclude bins at edges!

    for (Int_t i = 0; i < 18; i++){

        SigmaValues[0][i] = pValues415[i+1]/(Graph->Eval(425,0));
        SigmaErrValues[0][i] = pErrValues415[i+1]/(Graph->Eval(425,0));
        SigmaValues[1][i] = pValues425[i+1]/(Graph->Eval(425,0));
        SigmaErrValues[1][i] = pErrValues425[i+1]/(Graph->Eval(425,0));
        SigmaValues[2][i] = pValues435[i+1]/(Graph->Eval(435,0));
        SigmaErrValues[2][i] = pErrValues435[i+1]/(Graph->Eval(435,0));
        SigmaValues[3][i] = pValues445[i+1]/(Graph->Eval(445,0));
        SigmaErrValues[3][i] = pErrValues445[i+1]/(Graph->Eval(445,0));
        SigmaValues[4][i] = pValues455[i+1]/(Graph->Eval(455,0));
        SigmaErrValues[4][i] = pErrValues455[i+1]/(Graph->Eval(455,0));
        SigmaValues[5][i] = pValues465[i+1]/(Graph->Eval(465,0));
        SigmaErrValues[5][i] = pErrValues465[i+1]/(Graph->Eval(465,0));
        SigmaValues[6][i] = pValues475[i+1]/(Graph->Eval(475,0));
        SigmaErrValues[6][i] = pErrValues475[i+1]/(Graph->Eval(475,0));
        SigmaValues[7][i] = pValues485[i+1]/(Graph->Eval(485,0));
        SigmaErrValues[7][i] = pErrValues485[i+1]/(Graph->Eval(485,0));
        SigmaValues[8][i] = pValues495[i+1]/(Graph->Eval(495,0));
        SigmaErrValues[8][i] = pErrValues495[i+1]/(Graph->Eval(495,0));
        SigmaValues[9][i] = pValues505[i+1]/(Graph->Eval(505,0));
        SigmaErrValues[9][i] = pErrValues505[i+1]/(Graph->Eval(505,0));
        SigmaValues[10][i] = pValues515[i+1]/(Graph->Eval(515,0));
        SigmaErrValues[10][i] = pErrValues515[i+1]/(Graph->Eval(515,0));
        SigmaValues[11][i] = pValues525[i+1]/(Graph->Eval(525,0));
        SigmaErrValues[11][i] = pErrValues525[i+1]/(Graph->Eval(525,0));
        SigmaValues[12][i] = pValues535[i+1]/(Graph->Eval(535,0));
        SigmaErrValues[12][i] = pErrValues535[i+1]/(Graph->Eval(535,0));
        SigmaValues[13][i] = pValues545[i+1]/(Graph->Eval(545,0));
        SigmaErrValues[13][i] = pErrValues545[i+1]/(Graph->Eval(545,0));
        SigmaValues[14][i] = pValues555[i+1]/(Graph->Eval(555,0));
        SigmaErrValues[14][i] = pErrValues555[i+1]/(Graph->Eval(555,0));
        SigmaValues[15][i] = pValues565[i+1]/(Graph->Eval(565,0));
        SigmaErrValues[15][i] = pErrValues565[i+1]/(Graph->Eval(565,0));
        SigmaValues[16][i] = pValues575[i+1]/(Graph->Eval(575,0));
        SigmaErrValues[16][i] = pErrValues575[i+1]/(Graph->Eval(575,0));
        SigmaValues[17][i] = pValues585[i+1]/(Graph->Eval(585,0));
        SigmaErrValues[17][i] = pErrValues585[i+1]/(Graph->Eval(585,0));
        SigmaValues[18][i] = pValues595[i+1]/(Graph->Eval(595,0));
        SigmaErrValues[18][i] = pErrValues595[i+1]/(Graph->Eval(595,0));
        SigmaValues[19][i] = pValues605[i+1]/(Graph->Eval(605,0));
        SigmaErrValues[19][i] = pErrValues605[i+1]/(Graph->Eval(605,0));
        SigmaValues[20][i] = pValues615[i+1]/(Graph->Eval(615,0));
        SigmaErrValues[20][i] = pErrValues615[i+1]/(Graph->Eval(615,0));

    }

    TFile f3("Sigma_Plots_NS14_V2.root", "RECREATE");

    Float_t xMin = -1;
    Float_t xMax = 1;
    double x[18] = {0.85, 0.75, 0.65, 0.55, 0.45, 0.35, 0.25, 0.15, 0.05, -0.05, -0.15, -0.25, -0.35, -0.45, -0.55, -0.65, -0.75, -0.85}; // Need to adjust
    double ex[18] = {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05}; // Need to adjust

    for(int i = 0 ; i < 21 ; i++)
    {
        sprintf(name, "Sigma_%i", 415+(i*10));
        sprintf(title, "#Sigma(Cos#theta_{CM}) for E_{#gamma} %i #pm 10 MeV", 415+(i*10));
        SigmaPlots[i] = new TGraphErrors(18 , x, SigmaValues[i], ex, SigmaErrValues[i]);
        if (i != 14)SigmaPlots[i]->Fit("LegPol", "B");
        SigmaPlots[i]->SetMarkerColor(4);
        SigmaPlots[i]->SetLineColor(4);
        SigmaPlots[i]->SetMarkerStyle(8);
        SigmaPlots[i]->SetMarkerSize(1);
        SigmaPlots[i]->GetXaxis()->SetTitle("Cos#theta_{CM}");
        SigmaPlots[i]->GetYaxis()->SetTitle("#Sigma");
        SigmaPlots[i]->SetName(name);
        SigmaPlots[i]->SetTitle(title);

        for(Int_t j = 0; j < 7; j++){
            LegPar[j][i] = LegPol->GetParameter(j);
            LegParErr[j][i] = LegPol->GetParError(j);
        }

        SigmaPlots[i]->Write();

    }

    TCanvas *canvas20 = new TCanvas("canvas20","canvas20", 1920, 1080);
    canvas20->Divide(5,5);
    for(int i = 1 ; i < 22 ; i++){
        canvas20->cd(i);
        SigmaPlots[i-1]->Draw("AEP");
    }

    canvas20->Write();

    //Define new tree to store parameters in
    TTree* p0tree = new TTree("Parameter0_Values", "Tree_of_Values");
    p0tree->Branch("p0", &p0, "p0/D");
    p0tree->Branch("p0Err", &p0Err, "p0Err/D");
    for (Int_t m = 0; m < 21; m++){
        p0 = LegPar[0][m];
        p0Err = LegParErr[0][m];
        p0tree->Fill();
    }

    //Define new tree to store parameters in
    TTree* p1tree = new TTree("Parameter1_Values", "Tree_of_Values");
    p1tree->Branch("p1", &p1, "p1/D");
    p1tree->Branch("p1Err", &p1Err, "p1Err/D");
    for (Int_t m = 0; m < 21; m++){
        p1 = LegPar[1][m];
        p1Err = LegParErr[1][m];
        p1tree->Fill();
    }

    //Define new tree to store parameters in
    TTree* p2tree = new TTree("Parameter2_Values", "Tree_of_Values");
    p2tree->Branch("p2", &p2, "p2/D");
    p2tree->Branch("p2Err", &p2Err, "p2Err/D");
    for (Int_t m = 0; m < 21; m++){
        p2 = LegPar[2][m];
        p2Err = LegParErr[2][m];
        p2tree->Fill();
    }

    //Define new tree to store parameters in
    TTree* p3tree = new TTree("Parameter3_Values", "Tree_of_Values");
    p3tree->Branch("p3", &p3, "p3/D");
    p3tree->Branch("p3Err", &p3Err, "p3Err/D");
    for (Int_t m = 0; m < 21; m++){
        p3 = LegPar[3][m];
        p3Err = LegParErr[3][m];
        p3tree->Fill();
    }

    //Define new tree to store parameters in
    TTree* p4tree = new TTree("Parameter4_Values", "Tree_of_Values");
    p4tree->Branch("p4", &p4, "p4/D");
    p4tree->Branch("p4Err", &p4Err, "p4Err/D");
    for (Int_t m = 0; m < 21; m++){
        p4 = LegPar[4][m];
        p4Err = LegParErr[4][m];
        p4tree->Fill();
    }

    //Define new tree to store parameters in
    TTree* p5tree = new TTree("Parameter5_Values", "Tree_of_Values");
    p5tree->Branch("p5", &p5, "p5/D");
    p5tree->Branch("p5Err", &p5Err, "p5Err/D");
    for (Int_t m = 0; m < 21; m++){
        p5 = LegPar[5][m];
        p5Err = LegParErr[5][m];
        p5tree->Fill();
    }

    //Define new tree to store parameters in
    TTree* p6tree = new TTree("Parameter6_Values", "Tree_of_Values");
    p6tree->Branch("p6", &p6, "p6/D");
    p6tree->Branch("p6Err", &p6Err, "p6Err/D");
    for (Int_t m = 0; m < 21; m++){
        p6 = LegPar[6][m];
        p6Err = LegParErr[6][m];
        p6tree->Fill();
    }

    f3.Write();
}
