#include "./includes_Sigma_NoScatt_V2.h"

void Sigma_NoScatt_V2(){

    TFile *MBData = TFile::Open("/scratch/Mainz_Software/a2GoAT/Sig_res_St.root");
    TGraphErrors* SigmaPlots[21];
    TGraphErrors* ParameterPlots[6];
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

    TFile *f1= TFile::Open("/scratch/Mainz_Software/a2GoAT/Results/ParaPerpAsymm_NS17.root");
    TTree *t1 = (TTree*)f1->Get("Parameter_Values");
    // Old version
    //TF1 *LegPol = new TF1("LegPol", "[0]+[1]*x+[2]*(0.5*(3*x**2-1))+[3]*(0.5*(5*x**3-3*x))+[4]*(0.125*(35*x**4-30*x**2+3))+[5]*(1.0/8.0*(63*x**5-70*x**3+15*x))+[6]*(1.0/16*(231*x**6-315*x**4+105*x**2-5))", -1, 1);
    // New Version, 8 parameters
    TF1 *LegPol = new TF1("LegPol", "(1-x**2)*([0]*3+[1]*15*x+[2]*15.0/2*(7*x**2-1)+[3]*105.0/2*x*(3*x**2-1)+[4]*105.0/8*(33*x**4-18*x**2+1)+[5]*63.0/8*x*(143*x**4-110*x**2+15)+[6]*315.0/16*(143*x**6-143*x**4+33*x**2-1)+[7]*495.0/16*(221*x**7-273*x**5+91*x**3-7*x))", -1, 1);
    LegPol->SetLineColor(4);
    LegPol->SetLineWidth(2);
    LegPol->SetParLimits(6, 0.0, 0.0);
    LegPol->SetParLimits(7, 0.0, 0.0);
    LegPol->FixParameter(6, 0.0); // These seem to be ignored?
    LegPol->FixParameter(7, 0.0);

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

    TFile f3("Sigma_Plots_NS17.root", "RECREATE");

    Float_t xMin = -1;
    Float_t xMax = 1;
    double x[18] = {0.85, 0.75, 0.65, 0.55, 0.45, 0.35, 0.25, 0.15, 0.05, -0.05, -0.15, -0.25, -0.35, -0.45, -0.55, -0.65, -0.75, -0.85}; // Need to adjust
    double ex[18] = {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05}; // Need to adjust

    for(Int_t i = 0 ; i < 21 ; i++)
    {
        sprintf(name, "Sigma_%i", 415+(i*10));
        sprintf(title, "#Sigma(Cos#theta_{CM}) E_{#gamma} %i #pm 10 MeV", 415+(i*10));
        SigmaPlots[i] = new TGraphErrors(18 , x, SigmaValues[i], ex, SigmaErrValues[i]);

        if(i == 0){
            LegPol->SetParLimits(0, -1, 1);
            LegPol->SetParLimits(1, -1, 1);
            LegPol->SetParLimits(2, -1, 1);
            LegPol->SetParLimits(3, -1, 1);
            LegPol->SetParLimits(4, -1, 1);
            LegPol->SetParLimits(5, -1, 1);
        }

        if (i != 14) {
            SigmaPlots[i]->Fit("LegPol", "M");
            cout << "NDOF " << LegPol->GetNDF() << "   " << "Chi2 " << LegPol->GetChisquare() << "   " << "Chi2/DoF " << LegPol->GetChisquare()/LegPol->GetNDF() << endl;
        }

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

        for(Int_t j = 0; j < 8; j++){
            LegPar[j][i] = LegPol->GetParameter(j);
            LegParErr[j][i] = LegPol->GetParError(j);
        }

        SigmaPlots[i]->Write();

    }

    double x2[21] = {415, 425, 435, 445, 455, 465, 475, 485, 495, 505, 515, 525, 535, 545, 555, 565, 575, 585, 595, 605, 615}; // Need to adjust
    double ex2[21] = {5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5}; // Need to adjust

    for(Int_t k = 0; k < 6; k++){
        sprintf(name2, "P%i", k);
        sprintf(title2, "P^{2}_{%i}", k+2);
        ParameterPlots[k] = new TGraphErrors(21, x2, LegPar[k], ex2, LegParErr[k]);
        ParameterPlots[k]->SetMarkerColor(4);
        ParameterPlots[k]->SetLineColor(4);
        ParameterPlots[k]->SetMarkerStyle(8);
        ParameterPlots[k]->SetMarkerSize(1);
        ParameterPlots[k]->GetXaxis()->SetTitle("E_{#gamma}");
        ParameterPlots[k]->GetXaxis()->SetRangeUser(410, 620);
        ParameterPlots[k]->GetYaxis()->SetTitle(title2);
        ParameterPlots[k]->SetName(name2);
        ParameterPlots[k]->SetTitle(title2);
        ParameterPlots[k]->Write();
    }

    TCanvas *canvas20 = new TCanvas("canvas20","canvas20", 1920, 1080);
    canvas20->Divide(5,5);
    for(int i = 1 ; i < 22 ; i++){
        canvas20->cd(i);
        SigmaPlots[i-1]->Draw("AEP");
    }

    canvas20->Write();

    leg = new TLegend(0.1, 0.1, 0.9, 0.9);
    leg->AddEntry(MBHist[1], "APLCON Analysis", "lepz");
    leg->AddEntry(SigmaPlots[0], "Alternative Analysis", "lepz");

    TCanvas *canvas21 = new TCanvas("canvas21","canvas21", 1920, 1080);
    canvas21->Divide(5,5);
    for(int i = 1 ; i < 22 ; i++){
        canvas21->cd(i);
        if (i != 1){
            SigmaPlots[i-1]->Draw("AEP");
            MBHist[i-2]->Draw("SAMEEP");
        }
        elseif (i == 1){
            SigmaPlots[i-1]->Draw("AEP");
        }
    }
    canvas21->cd(22);
    leg->Draw();

    canvas21->Write();

    leg2 = new TLegend(0.1, 0.1, 0.9, 0.9);
    leg2->AddEntry(e435, "Previous Results", "lepz");
    leg2->AddEntry(SigmaPlots[0], "This Analysis", "lepz");

    TCanvas *canvas22 = new TCanvas("canvas22","canvas22", 1920, 1080);
    canvas22->Divide(5,3);
    canvas22->cd(1);
    SigmaPlots[0]->Draw("AEP");
    e412->Draw("sameepz");
    canvas22->cd(2);
    SigmaPlots[2]->Draw("AEP");
    e435->Draw("sameepz");
    canvas22->cd(3);
    SigmaPlots[4]->Draw("AEP");
    e455->Draw("sameepz");
    canvas22->cd(4);
    SigmaPlots[5]->Draw("AEP");
    e465->Draw("sameepz");
    canvas22->cd(5);
    SigmaPlots[6]->Draw("AEP");
    e475->Draw("sameepz");
    canvas22->cd(6);
    SigmaPlots[7]->Draw("AEP");
    e485->Draw("sameepz");
    canvas22->cd(7);
    SigmaPlots[9]->Draw("AEP");
    e505->Draw("sameepz");
    canvas22->cd(8);
    SigmaPlots[10]->Draw("AEP");
    e515->Draw("sameepz");
    canvas22->cd(9);
    SigmaPlots[11]->Draw("AEP");
    e525->Draw("sameepz");
    canvas22->cd(10);
    SigmaPlots[13]->Draw("AEP");
    e545->Draw("sameepz");
    canvas22->cd(11);
    SigmaPlots[16]->Draw("AEP");
    e575->Draw("sameepz");
    canvas22->cd(12);
    SigmaPlots[18]->Draw("AEP");
    e595->Draw("sameepz");
    canvas22->cd(13);
    SigmaPlots[20]->Draw("AEP");
    e615->Draw("sameepz");
    canvas22->cd(14);
    leg2->Draw();

    canvas22->Write();

    TCanvas *canvas23 = new TCanvas("canvas23","canvas23", 1920, 1080);
    canvas23->Divide(3,2);
    for(int i = 1 ; i < 7 ; i++){
        canvas23->cd(i);
        ParameterPlots[i-1]->Draw("AEP");
    }

    canvas23->Write();

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

    //Define new tree to store parameters in
    TTree* p7tree = new TTree("Parameter7_Values", "Tree_of_Values");
    p7tree->Branch("p7", &p7, "p7/D");
    p7tree->Branch("p7Err", &p7Err, "p7Err/D");
    for (Int_t m = 0; m < 21; m++){
        p7 = LegPar[7][m];
        p7Err = LegParErr[7][m];
        p7tree->Fill();
    }

    f3.Write();
}
