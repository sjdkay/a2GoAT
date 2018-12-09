#include "./includes_Sigma_NoScatt_V3.h"

void Sigma_NoScatt_V3(){

    gStyle->SetTitleFontSize(0.06);

    char HistName[60];
    char NewHistName[60];
    TH1F* PhipPara[21][20];
    TH1F* PhipPerp[21][20];
    double EStart = 415; // What is the initial EGamma bin? Change this depending on value
    Double_t* SystematicArray[21];
    double Systematic[21][18];
    double Sigma90[20];
    double SigmaErr90[20];
    TGraphErrors* Sigma90DegCM;
    TGraphErrors* Sigma90DegCMAdamian;
    TGraphErrors* Sigma90DegCMGorbenko;

    TF1 *CosFunc = new TF1("CosFit", "[0]*cos((2*x*TMath::DegToRad())+acos(0))");
    CosFunc->SetParNames("Amplitude");

    TH1F* AsymmHists[21][20];

    char ParaHistName[60];
    char PerpHistName[60];
    char AsymmHistName[60];
    char GraphName[60];
    char ScGraphName[60];
    char ScGraphNameAdj[60];
    char FitName[60];
    char FitErrName[60];
    char FitGraphName[60];

    TF1 *LegPol = new TF1("LegPol", "(1-x**2)*([0]*3+[1]*15*x+[2]*15.0/2*(7*x**2-1)+[3]*105.0/2*x*(3*x**2-1)+[4]*105.0/8*(33*x**4-18*x**2+1)+[5]*63.0/8*x*(143*x**4-110*x**2+15)+[6]*315.0/16*(143*x**6-143*x**4+33*x**2-1)+[7]*495.0/16*(221*x**7-273*x**5+91*x**3-7*x))", -1, 1);
    LegPol->SetLineColor(4);
    LegPol->SetLineWidth(2);
    LegPol->SetParLimits(6, 0.0, 0.0);
    LegPol->SetParLimits(7, 0.0, 0.0);
    LegPol->FixParameter(6, 0.0); // These seem to be ignored?
    LegPol->FixParameter(7, 0.0);

    TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/Para/NoScatt/Physics_Total_Para_NoScatt_18_26_4_18_V2.root"); // Open latest Para file

    TH1D* Eg_Para = (TH1D*)f->Get("Eg2")->Clone();
    Eg_Para->SetName("Eg_Para");

    for(Int_t i = 0; i < 21; i++){ // Energy
        for(Int_t j = 0; j < 20; j++){ // Theta
            sprintf(HistName, "Phip_%iMeVCM%i", EStart+(i*10) , j+1);
            sprintf(NewHistName, "Phip_%iMeVCM%i_Para", EStart+(i*10) , j+1);
            PhipPara[i][j] = ((TH1F*)f->Get(HistName));
            PhipPara[i][j]->SetName(NewHistName);
        }
    }


    /////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////
    ////////////////////////// Para Done ////////////////////////////
    /////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////

    TFile *f1 = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/Perp/NoScatt/Physics_Total_Perp_NoScatt_18_26_4_18_V2.root"); // Open latest Para file

    TH1D* Eg_Perp = (TH1D*)f1->Get("Eg2")->Clone();
    Eg_Perp->SetName("Eg_Perp");

    for(Int_t i = 0; i < 21; i++){ // Energy
        for(Int_t j = 0; j < 20; j++){ // Theta
            sprintf(HistName, "Phip_%iMeVCM%i", EStart+(i*10) , j+1);
            sprintf(NewHistName, "Phip_%iMeVCM%i_Perp", EStart+(i*10) , j+1);
            PhipPerp[i][j] = ((TH1F*)f1->Get(HistName));
            PhipPerp[i][j]->SetName(NewHistName);
        }
    }

    TFile f2("ParaPerp_NS18_Combined_V2.root", "RECREATE");

    Eg_Para->Write();
    Eg_Perp->Write();

    for(Int_t i = 0; i < 21; i++){ // Energy
        for(Int_t j = 0; j < 20; j++){ // Theta
            PhipPara[i][j]->Write();
            PhipPerp[i][j]->Write();
        }
    }

    f2.Write();

    NPara = Eg_Para->GetEntries();
    NPerp = Eg_Perp->GetEntries();
    ScaleFactor = NPara/NPerp;
    ScaleFactorErr = sqrt( (NPara/((TMath::Power(NPerp,2)))) + (((TMath::Power(NPara,2)))/(TMath::Power(NPerp,3))) ); // Error Propagation of above formula, see notebook 5

    for(Int_t i = 0; i < 21; i++){ // Energy
        for(Int_t j = 0; j < 20; j++){ // Theta
            sprintf(ParaHistName, "Phip_%iMeVCM%i_Para", EStart+(i*10) , j+1);
            sprintf(PerpHistName, "Phip_%iMeVCM%i_Perp", EStart+(i*10) , j+1);
            sprintf(AsymmHistName, "Sigma_%iMeV_CM%i_Hist", EStart+(i*10) , j+1);
            AsymmHists[i][j] = (TH1F*) (((TH1F*)f2.Get(ParaHistName))->GetAsymmetry(((TH1F*)f2.Get(PerpHistName)), ScaleFactor, ScaleFactorErr)));
            AsymmHists[i][j]->SetName(AsymmHistName);
            AsymmHists[i][j]-> Fit("CosFit", "Q");
            cout << "NDOF " << CosFit->GetNDF() << "   " << "Chi2 " << CosFit->GetChisquare() << "   " << "Chi2/DoF " << CosFit->GetChisquare()/CosFit->GetNDF() << endl;
            pCosAmp[i][j] = CosFit->GetParameter(0);
            pCosAmpErr[i][j] = CosFit->GetParError(0);
        }
    }

    TFile f3("ParaPerpAsymm_NS18_V2.root", "RECREATE");
    // In other version fill a tree here to read for next file - Don't bother now?
    // Add back in later maybe?

    for(Int_t i = 0; i < 21; i++){ // Energy
        for(Int_t j = 0; j < 20; j++){ // Theta
            AsymmHists[i][j]->Write();
        }
    }

    f3.Write();

    TFile *MBData = TFile::Open("/scratch/Mainz_Software/a2GoAT/Sig_res_St.root");
    TGraphErrors* SigmaPlots[21];
    TGraphErrors* SigmaSystPlots[21];
    TGraphErrors* SigmaSystDiffPlots[21];
    TGraphErrors* ParameterPlots[6];
    TGraphErrors* SigmaScPlots[10];
    char name[21];
    char title[60];
    char MBname[20];
    char name2[60];
    char title2[60];

    TH1F* MBHist[20];

    for (Int_t i = 0; i < 20; i++){ // Get Mikhail plots
        sprintf(MBname, "hslC%i", 43+i);
        MBHist[i] = (TH1F*)MBData->Get(MBname);
    }

    // Get previous data plots

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

    TFile *f4= TFile::Open("/scratch/Mainz_Software/a2GoAT/LinPol_Aug16.root"); // Open linear polarisation plot

    // Calculate values of sigma for each angular and energy bin
    // Exclude bins at edges!

    for(Int_t i = 0; i < 21; i++){ // Energy
        for(Int_t j = 0; j < 18; j++){ // Theta
            Double_t EValue = EStart + (i*10);
            Double_t PolVal = Graph->Eval(EValue, 0);
            Double_t PolErrVal = (3./100.)*PolVal; // 3% error on Polarisation assumed here
            SigmaValues[i][j] = pCosAmp[i][j+1]/(PolVal);
            SigmaErrValues[i][j] = sqrt(((pCosAmpErr[i][j+1]/(PolVal))**2) + (((PolErrVal*pCosAmp[i][j+1])/(PolVal**2))**2));
        }
    }

    // Open file to grap systematic errors for 1 sigma vs 2 sigma MM cut
    TFile *fSyst = TFile::Open("/scratch/Mainz_Software/a2GoAT/Sigma_Systematic_18_17_V5.root");

    double xSc[5] = {0.8, 0.4, 0, -0.4, -0.8};
    double exSc[5] = {0.2, 0.2, 0.2, 0.2, 0.2};

    for(Int_t i = 0; i < 21; i ++){ // Egamma value
        sprintf(GraphName, "SigmaSyst_%i", EStart+(i*10));
        SigmaSystDiffPlots[i] = ((TGraphErrors*)fSyst->Get(GraphName));
        Mean = SigmaSystDiffPlots[i]->GetMean(2);
        RMS = SigmaSystDiffPlots[i]->GetRMS(2);
        SystVal = Mean + RMS;
        for(Int_t j = 0; j < 18; j ++){
            Systematic[i][j] = -1+(fabs(SystVal));
            if(i == 14) Systematic[i][j] = -1;
        }
    }

    // Open file to get latest values for Sigma from scattered data
    TFile *fScatt = TFile::Open("/scratch/Mainz_Software/a2GoAT/Sigma_Plots_S48_V2.root");

    for(Int_t i = 0; i < 10; i ++){ // Egamma value
        sprintf(ScGraphName, "Sigma_%i", EStart+(i*20));
        sprintf(ScGraphNameAdj, "SigmaSc_%i", EStart+(i*20)-5);
        SigmaScPlots[i] = ((TGraphErrors*)fScatt->Get(ScGraphName));
        SigmaScPlots[i]->SetName(ScGraphNameAdj);
        SigmaScPlots[i]->SetLineColor(618);
        SigmaScPlots[i]->SetMarkerColor(618);
        SigmaScPlots[i]->SetMarkerStyle(33);
        SigmaScPlots[i]->SetMarkerSize(1.5);
        SigmaScPlots[i]->GetXaxis()->SetLabelSize(0.06);
        SigmaScPlots[i]->GetXaxis()->SetTitleSize(0.06);
        SigmaScPlots[i]->GetXaxis()->SetTitleOffset(0.7);
        SigmaScPlots[i]->GetXaxis()->CenterTitle();
        SigmaScPlots[i]->GetYaxis()->SetLabelSize(0.06);
        SigmaScPlots[i]->GetYaxis()->SetTitleSize(0.08);
        SigmaScPlots[i]->GetYaxis()->SetTitleOffset(0.5);
        SigmaScPlots[i]->GetYaxis()->CenterTitle();
    }

    TFile f5("Sigma_Plots_NS18_S48_V3_2.root", "RECREATE");

    Float_t xMin = -1;
    Float_t xMax = 1;
    double x[18] = {0.85, 0.75, 0.65, 0.55, 0.45, 0.35, 0.25, 0.15, 0.05, -0.05, -0.15, -0.25, -0.35, -0.45, -0.55, -0.65, -0.75, -0.85}; // Need to adjust
    double ex[18] = {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05}; // Need to adjust

    for(Int_t i = 0 ; i < 21 ; i++)
    {
        sprintf(name, "Sigma_%i", EStart+(i*10));
        sprintf(title, "#Sigma(Cos#theta_{CM}) E_{#gamma} %i #pm 5 MeV", EStart+(i*10));
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
        SigmaPlots[i]->GetXaxis()->SetLabelSize(0.05);
        SigmaPlots[i]->GetXaxis()->SetTitleSize(0.06);
        SigmaPlots[i]->GetXaxis()->SetTitleOffset(0.7);
        SigmaPlots[i]->GetXaxis()->CenterTitle();
        SigmaPlots[i]->GetYaxis()->SetLabelSize(0.06);
        SigmaPlots[i]->GetYaxis()->SetTitleSize(0.08);
        SigmaPlots[i]->GetYaxis()->SetTitleOffset(0.5);
        SigmaPlots[i]->GetYaxis()->CenterTitle();

        for(Int_t j = 0; j < 8; j++){
            LegPar[j][i] = LegPol->GetParameter(j);
            LegParErr[j][i] = LegPol->GetParError(j);
        }
        SigmaPlots[i]->Write();
    }

    for(Int_t i = 0; i < 10; i ++){ // Egamma value
        SigmaScPlots[i]->Write();
    }


    double x2[21];
    for(Int_t i = 0 ; i < 21 ; i++){
        x2[i] = EStart + (i*10);
    }
    double ex2[21] = {5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5};

    for(Int_t k = 0; k < 6; k++){
        sprintf(name2, "P%i", k);
        sprintf(title2, "P^{2}_{%i}", k+2);
        ParameterPlots[k] = new TGraphErrors(21, x2, LegPar[k], ex2, LegParErr[k]);
        ParameterPlots[k]->SetMarkerColor(4);
        ParameterPlots[k]->SetLineColor(4);
        ParameterPlots[k]->SetMarkerStyle(8);
        ParameterPlots[k]->SetMarkerSize(1);
        ParameterPlots[k]->GetXaxis()->SetTitle("E_{#gamma}");
        ParameterPlots[k]->GetXaxis()->SetRangeUser(400, 620);
        ParameterPlots[k]->GetYaxis()->SetTitle(title2);
        ParameterPlots[k]->GetXaxis()->SetLabelSize(0.04);
        ParameterPlots[k]->GetXaxis()->SetTitleSize(0.06);
        ParameterPlots[k]->GetXaxis()->SetTitleOffset(0.7);
        ParameterPlots[k]->GetXaxis()->CenterTitle();
        ParameterPlots[k]->GetYaxis()->SetLabelSize(0.02);
        ParameterPlots[k]->GetYaxis()->SetTitleSize(0.06);
        ParameterPlots[k]->GetYaxis()->SetTitleOffset(0.75);
        ParameterPlots[k]->GetYaxis()->CenterTitle();
        ParameterPlots[k]->SetName(name2);
        ParameterPlots[k]->RemovePoint(14);
        ParameterPlots[k]->SetTitle(title2);
        ParameterPlots[k]->Write();
    }

    // Draw all plots for Sigma without systematics
    TCanvas *canvas20 = new TCanvas("canvas20","canvas20", 2560, 1440);
    canvas20->Divide(5,4, 0.0000001, 0.001);
    for(int i = 0 ; i < 21 ; i++){
        if (i < 14){
            canvas20->cd(i+1);
            SigmaPlots[i]->Draw("AEP");
        }
        if (i > 14){
            canvas20->cd(i);
            SigmaPlots[i]->Draw("AEP");
        }
    }

    canvas20->Write();

    // Draw all plots for Sigma with systematics
    TCanvas *canvas20a = new TCanvas("canvas20a","canvas20a", 2560, 1440);
    canvas20a->Divide(2,2, 0.0000001, 0.001);
    for(int i = 0 ; i < 4 ; i++){
        canvas20a->cd(i+1);
        SigmaPlots[i]->Draw("AEP");
        SigmaSystPlots[i] = new TGraphErrors(18 , x, Systematic[i], ex, 0);
        SigmaSystPlots[i]->SetFillColor(25);
        SigmaSystPlots[i]->Draw("SAMEB1");
    }
    canvas20a->Write();

    TCanvas *canvas20b = new TCanvas("canvas20b","canvas20b", 2560, 1440);
    canvas20b->Divide(2,2, 0.0000001, 0.001);
    for(int i = 4 ; i < 8 ; i++){
        canvas20b->cd(i-3);
        SigmaPlots[i]->Draw("AEP");
        SigmaSystPlots[i] = new TGraphErrors(18 , x, Systematic[i], ex, 0);
        SigmaSystPlots[i]->SetFillColor(25);
        SigmaSystPlots[i]->Draw("SAMEB1");
    }
    canvas20b->Write();

    TCanvas *canvas20c = new TCanvas("canvas20c","canvas20c", 2560, 1440);
    canvas20c->Divide(2,2, 0.0000001, 0.001);
    for(int i = 8 ; i < 12 ; i++){
        canvas20c->cd(i-7);
        SigmaPlots[i]->Draw("AEP");
        SigmaSystPlots[i] = new TGraphErrors(18 , x, Systematic[i], ex, 0);
        SigmaSystPlots[i]->SetFillColor(25);
        SigmaSystPlots[i]->Draw("SAMEB1");
    }
    canvas20c->Write();

    TCanvas *canvas20d = new TCanvas("canvas20d","canvas20d", 2560, 1440);
    canvas20d->Divide(2,2, 0.0000001, 0.001);
    for(int i = 0 ; i < 4 ; i++){
        canvas20d->cd(i+1);
        if(i==0){
            SigmaPlots[12]->Draw("AEP");
            SigmaSystPlots[12] = new TGraphErrors(18 , x, Systematic[12], ex, 0);
            SigmaSystPlots[12]->SetFillColor(25);
            SigmaSystPlots[12]->Draw("SAMEB1");
        }
        if(i==1){
            SigmaPlots[13]->Draw("AEP");
            SigmaSystPlots[13] = new TGraphErrors(18 , x, Systematic[12], ex, 0);
            SigmaSystPlots[13]->SetFillColor(25);
            SigmaSystPlots[13]->Draw("SAMEB1");
        }
        if(i==2){
            SigmaPlots[15]->Draw("AEP");
            SigmaSystPlots[15] = new TGraphErrors(18 , x, Systematic[15], ex, 0);
            SigmaSystPlots[15]->SetFillColor(25);
            SigmaSystPlots[15]->Draw("SAMEB1");
        }
        if(i==3){
            SigmaPlots[16]->Draw("AEP");
            SigmaSystPlots[16] = new TGraphErrors(18 , x, Systematic[16], ex, 0);
            SigmaSystPlots[16]->SetFillColor(25);
            SigmaSystPlots[16]->Draw("SAMEB1");
        }
    }
    canvas20d->Write();

    TCanvas *canvas20e = new TCanvas("canvas20e","canvas20e", 2560, 1440);
    canvas20e->Divide(2,2, 0.0000001, 0.001);
    for(int i = 17 ; i < 21 ; i++){
        canvas20e->cd(i-16);
        SigmaPlots[i]->Draw("AEP");
        SigmaSystPlots[i] = new TGraphErrors(18 , x, Systematic[i], ex, 0);
        SigmaSystPlots[i]->SetFillColor(25);
        SigmaSystPlots[i]->Draw("SAMEB1");
    }
    canvas20e->Write();

    // Draw all plots for Sigma without systematics
    TCanvas *canvas20f = new TCanvas("canvas20f","canvas20f", 2560, 1440);
    canvas20f->Divide(5,4, 0.0000001, 0.001);
    for(int i = 0 ; i < 20 ; i++){
        if (i < 14){
            canvas20f->cd(i+1);
            SigmaPlots[i]->Draw("AEP");
            if(i % 2 == kFALSE){
                SigmaScPlots[i/2]->Draw("SAMEEP");
            }
            else if (i % 2){
                SigmaScPlots[(i-1)/2]->Draw("SAMEEP");
            }
        }

        if (i > 14){
            canvas20f->cd(i);
            SigmaPlots[i]->Draw("AEP");
            if(i % 2 == kFALSE){
                SigmaScPlots[i/2]->Draw("SAMEEP");
            }
            else if (i % 2){
                SigmaScPlots[(i-1)/2]->Draw("SAMEEP");
            }
        }
    }

    legend = new TLegend(0.1, 0.1, 0.9, 0.9);
    legend->AddEntry(SigmaPlots[1], "Non-Scattered Events", "lepz");
    legend->AddEntry(SigmaScPlots[1], "Scattered Events", "lepz");
    canvas20f->cd(20);
    legend->Draw();

    canvas20f->Write();

    leg = new TLegend(0.1, 0.1, 0.9, 0.9);
    leg->AddEntry(MBHist[1], "APLCON Analysis", "lepz");
    leg->AddEntry(SigmaPlots[0], "Alternative Analysis", "lepz");

    TCanvas *canvas21 = new TCanvas("canvas21","canvas21", 2560, 1440);
    canvas21->Divide(5,5, 0.0000001, 0.001);
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

    TCanvas *canvas22 = new TCanvas("canvas22","canvas22", 2560, 1440);
    canvas22->Divide(5,3, 0.0000001, 0.001);
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

    TCanvas *canvas22a = new TCanvas("canvas22a","canvas22a", 2560, 1440);
    canvas22a->Divide(5,3, 0.0000001, 0.001);
    canvas22a->cd(1);
    SigmaPlots[0]->Draw("AEP");
    e412->Draw("sameepz");
    SigmaSystPlots[0]->Draw("SAMEB1");
    canvas22a->cd(2);
    SigmaPlots[2]->Draw("AEP");
    e435->Draw("sameepz");
    SigmaSystPlots[2]->Draw("SAMEB1");
    canvas22a->cd(3);
    SigmaPlots[4]->Draw("AEP");
    e455->Draw("sameepz");
    SigmaSystPlots[4]->Draw("SAMEB1");
    canvas22a->cd(4);
    SigmaPlots[5]->Draw("AEP");
    e465->Draw("sameepz");
    SigmaSystPlots[5]->Draw("SAMEB1");
    canvas22a->cd(5);
    SigmaPlots[6]->Draw("AEP");
    e475->Draw("sameepz");
    SigmaSystPlots[6]->Draw("SAMEB1");
    canvas22a->cd(6);
    SigmaPlots[7]->Draw("AEP");
    e485->Draw("sameepz");
    SigmaSystPlots[7]->Draw("SAMEB1");
    canvas22a->cd(7);
    SigmaPlots[9]->Draw("AEP");
    e505->Draw("sameepz");
    SigmaSystPlots[9]->Draw("SAMEB1");
    canvas22a->cd(8);
    SigmaPlots[10]->Draw("AEP");
    e515->Draw("sameepz");
    SigmaSystPlots[10]->Draw("SAMEB1");
    canvas22a->cd(9);
    SigmaPlots[11]->Draw("AEP");
    e525->Draw("sameepz");
    SigmaSystPlots[11]->Draw("SAMEB1");
    canvas22a->cd(10);
    SigmaPlots[13]->Draw("AEP");
    e545->Draw("sameepz");
    SigmaSystPlots[13]->Draw("SAMEB1");
    canvas22a->cd(11);
    SigmaPlots[16]->Draw("AEP");
    e575->Draw("sameepz");
    SigmaSystPlots[16]->Draw("SAMEB1");
    canvas22a->cd(12);
    SigmaPlots[18]->Draw("AEP");
    e595->Draw("sameepz");
    SigmaSystPlots[18]->Draw("SAMEB1");
    canvas22a->cd(13);
    SigmaPlots[20]->Draw("AEP");
    e615->Draw("sameepz");
    SigmaSystPlots[20]->Draw("SAMEB1");
    canvas22a->cd(14);
    leg2->Draw();

    canvas22a->Write();

    // Average bins at pm 0.05 Cos Theta to get a bin at 90, plot these as funtion of EGamma

    double Adamian90x[12] = {414, 433, 451, 476, 504, 528, 552, 575, 598, 620, 644, 672};
    double Adamian90y[12] = {-0.27, -0.23, -0.23, -0.09, -0.15, -0.08, 0.0, 0.02, 0.18, 0.15, 0.27, 0.35};
    double Adamian90yErr[12] = {0.04, 0.04, 0.04, 0.04, 0.05, 0.05, 0.05, 0.06, 0.07, 0.07, 0.07, 0.06};
    double Gorbenko90x[6] = {400, 400, 450, 500, 550, 600};
    double Gorbenko90y[6] = {-0.267, -0.257, -0.164, -0.150, 0.0, 0.126};
    double Gorbenko90yErr[6] = {0.061, 0.033, 0.028, 0.057, 0.094, 0.142};

    for(Int_t i = 0; i < 21; i++){
        if(i < 14){
            Sigma90[i] = (SigmaValues[i][8] + SigmaValues[i][9])/2;
            SigmaErr90[i] = 0.5*(sqrt( (SigmaErrValues[i][8]**2) + (SigmaErrValues[i][9]**2)));
        }
        else if (i==14) continue;
        else if (i > 14){
            Sigma90[i-1] = (SigmaValues[i][8] + SigmaValues[i][9])/2;
            SigmaErr90[i-1] = 0.5*(sqrt( (SigmaErrValues[i][8]**2) + (SigmaErrValues[i][9]**2)));
        }
    }

    double x3[20] = {415, 425, 435, 445, 455, 465, 475, 485, 495, 505, 515, 525, 535, 545, 565, 575, 585, 595, 605, 615};
    double ex3[20] = {5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5};

    Sigma90DegCM = new TGraphErrors(20, x3, Sigma90, ex3, SigmaErr90);
    Sigma90DegCM->SetMarkerColor(4);
    Sigma90DegCM->SetLineColor(4);
    Sigma90DegCM->SetMarkerStyle(21);
    Sigma90DegCM->SetMarkerSize(2);
    Sigma90DegCM->GetXaxis()->SetTitle("E_{#gamma}");
    Sigma90DegCM->GetXaxis()->SetRangeUser(400, 700);
    Sigma90DegCM->GetYaxis()->SetTitle("#Sigma(90^{#circ})");
    Sigma90DegCM->GetXaxis()->SetLabelSize(0.04);
    Sigma90DegCM->GetXaxis()->SetTitleSize(0.06);
    Sigma90DegCM->GetXaxis()->SetTitleOffset(0.7);
    Sigma90DegCM->GetXaxis()->CenterTitle();
    Sigma90DegCM->GetYaxis()->SetLabelSize(0.04);
    Sigma90DegCM->GetYaxis()->SetTitleSize(0.06);
    Sigma90DegCM->GetYaxis()->SetTitleOffset(0.75);
    Sigma90DegCM->GetYaxis()->CenterTitle();
    Sigma90DegCM->SetName("Sigma90DegCM");
    Sigma90DegCM->SetTitle("#Sigma(E_{#gamma})");
    Sigma90DegCM->Write();

    Sigma90DegCMAdamian = new TGraphErrors(12, Adamian90x, Adamian90y, 0, Adamian90yErr);
    Sigma90DegCMAdamian->SetMarkerColor(1);
    Sigma90DegCMAdamian->SetLineColor(1);
    Sigma90DegCMAdamian->SetMarkerStyle(20);
    Sigma90DegCMAdamian->SetMarkerSize(2);
    Sigma90DegCMAdamian->SetName("Sigma90DegCMAdamian");
    Sigma90DegCMAdamian->Write();

    Sigma90DegCMGorbenko = new TGraphErrors(6, Gorbenko90x, Gorbenko90y, 0, Gorbenko90yErr);
    Sigma90DegCMGorbenko->SetMarkerColor(2);
    Sigma90DegCMGorbenko->SetLineColor(2);
    Sigma90DegCMGorbenko->SetMarkerStyle(22);
    Sigma90DegCMGorbenko->SetMarkerSize(2);
    Sigma90DegCMGorbenko->SetName("Sigma90DegCMGorbenko");
    Sigma90DegCMGorbenko->Write();

    leg3 = new TLegend(0.7, 0.2, 0.9, 0.4);
    leg3->AddEntry(Sigma90DegCM, "This Work", "lepz");
    leg3->AddEntry(Sigma90DegCMAdamian, "Adamian91", "lepz");
    leg3->AddEntry(Sigma90DegCMGorbenko, "Gorbenko82", "lepz");

    TCanvas *canvas22b = new TCanvas("canvas22b","canvas22b", 2560, 1440);
    Sigma90DegCM->Draw("AEP");
    Sigma90DegCMAdamian->Draw("SAMEEP");
    Sigma90DegCMGorbenko->Draw("SAMEEP");
    leg3->Draw("SAME");
    canvas22b->Write();

    TCanvas *canvas23 = new TCanvas("canvas23","canvas23", 2560, 1440);
    canvas23->Divide(3,2, 0.0000001, 0.001);
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

// Couting loop for results to be inputted into LaTeX if needed
//    for(Int_t i = 0; i < 21; i++){ // Energy
//        Double_t EValue = EStart + (i*10);
//        cout << std::setprecision(3) << EValue << " #pm " << 5 << "\t";
//        for (Int_t j = 0; j < 18; j++){
//                if (j!= 17) cout << std::setprecision(3) << std::fixed << SigmaValues[i][j] << " #pm " << SigmaErrValues[i][j] << "\t";
//                else if (j == 17) cout << SigmaValues[i][j] << " #pm " << SigmaErrValues[i][j] << endl;
//            }
//        }

    double SigmaEDepPar[6][7] = {{2.99203,-11.6753,22.4077,-27.858,23.5892,-12.9431,3.63412},
                                {-1.67285,5.72189,-9.9412,11.3262,-8.83277,4.43907,-1.10018},
                                {0.648495,-1.96551,3.04773,-3.14691,2.31897,-1.16499,0.309782},
                                {-0.884155,3.09574,-5.41028,6.06124,-4.58708,2.23079,-0.544889},
                                {-0.134468,0.560091,-1.0877,1.2847,-0.972672,0.445549,-0.0954476},
                                {-0.460013,1.41646,-2.17844,2.15514,-1.45138,0.633165,-0.139968}}; // Par 0 is P22 e.t.c.
    double SigmaEDepParErr[6][7] = {{0.0186005,0.0245129,0.0207498,0.0193895,0.0178788,0.0180634,0.0120939},
                                    {0.0102835,0.0134829,0.0113874,0.0105955,0.0097473,0.00981033,0.0065908},
                                    {0.00698575,0.00911475,0.00767362,0.00710237,0.0065157,0.00653866,0.00441064},
                                    {0.00542185,0.00703949,0.00593895,0.00550025,0.00505964,0.00508222,0.00344486},
                                    {0.00423913,0.00553335,0.0046663,0.00433196,0.00399039,0.0040243,0.00272008},
                                    {0.00353928,0.00460135,0.00388868,0.00360703,0.00332274,0.00333902,0.00226802}};

    // Define Energy dependence fits of 7 gaussians to the Sigma values
    TF1 **EDepP2 = new TF1*[6];
    double EDepP2FitValsX[6][21];
    double EDepP2FitValsY[6][21];
    // Define Energy dependence fits of 7 gaussians to the Sigma values for inclusion of errors
    TF1 **EDepErrP2 = new TF1*[6];
    double EDepErrP2FitValsY[6][21];
    TGraphErrors* EDepP2Graphs[6];

    for(Int_t i = 0; i < 6; i++){
        sprintf(FitName, "EDepP2%i", i+2);
        sprintf(FitErrName, "EDepP2%iErr", i+2);
        EDepP2[i] = new TF1 (FitName,"[0]*exp(-1*(x-420)**2/2/60**2)+[1]*exp(-1*(x-450)**2/2/60**2)+[2]*exp(-1*(x-480)**2/2/60**2)+[3]*exp(-1*(x-510)**2/2/60**2)+[4]*exp(-1*(x-540)**2/2/60**2)+[5]*exp(-1*(x-570)**2/2/60**2)+[6]*exp(-1*(x-600)**2/2/60**2)" , 410, 620);
        EDepErrP2[i] = new TF1 (FitErrName,"[0]*exp(-1*(x-420)**2/2/60**2)+[1]*exp(-1*(x-450)**2/2/60**2)+[2]*exp(-1*(x-480)**2/2/60**2)+[3]*exp(-1*(x-510)**2/2/60**2)+[4]*exp(-1*(x-540)**2/2/60**2)+[5]*exp(-1*(x-570)**2/2/60**2)+[6]*exp(-1*(x-600)**2/2/60**2)" , 410, 620);
        for(Int_t j =0 ; j<7; j++){
            EDepP2[i]->SetParameter(j, SigmaEDepPar[i][j]);
            EDepErrP2[i]->SetParameter(j, SigmaEDepPar[i][j] + 0.5*SigmaEDepParErr[i][j]);
        }
        for(Int_t k = 0; k < 21; k++){
            EDepP2FitValsX[i][k] = 415+(k*10);
            EDepP2FitValsY[i][k] = EDepP2[i]->Eval(415+(k*10));
            EDepErrP2FitValsY[i][k] = fabs(EDepP2[i]->Eval(415+(k*10)) - EDepErrP2[i]->Eval(415+(k*15)))/2;
        }
    }

    for(int i = 0; i < 6; i++){
        sprintf(FitGraphName, "EDepP2%i_Graph", i+2);
        EDepP2Graphs[i] = new TGraphErrors(21 , EDepP2FitValsX[i], EDepP2FitValsY[i], ex2, EDepErrP2FitValsY[i]);
        EDepP2Graphs[i]->SetFillColor(2);
        EDepP2Graphs[i]->SetFillStyle(3005);
        EDepP2Graphs[i]->SetName(FitGraphName);
        EDepP2Graphs[i]->Write();
    }

    TCanvas *canvas24 = new TCanvas("canvas24","canvas24", 2560, 1440);
    canvas24->Divide(3,2, 0.0000001, 0.001);
    for(int i = 1; i < 7; i++){
        canvas24->cd(i);
        ParameterPlots[i-1]->Draw("AEP");
        EDepP2[i-1]->Draw("SAME");
    }
    canvas24->Write();

    TCanvas *canvas25 = new TCanvas("canvas25","canvas25", 2560, 1440);
    canvas25->Divide(3,2, 0.0000001, 0.001);
    for(int i = 1; i < 7; i++){
        canvas25->cd(i);
        ParameterPlots[i-1]->Draw("AEP");
        EDepP2[i-1]->Draw("SAME");
        EDepP2Graphs[i-1]->Draw("SAMEE3");
    }
    canvas25->Write();

    double SigmaEDepPar2[6][4] = {{-0.108754,-0.00356053,0.0745446,-0.0212027},
                                   {0.0238134,-0.00919569,0.0305643,-0.0102451},
                                    {0.00121507,-0.00177259,0.0044306,-0.00475281},
                                    {0.00533593,-0.000667558,0.00358099,-0.00188308},
                                    {0.0029606,-0.00206478,0.00386589,-0.00364641},
                                    {0.00146578,-0.000521466,0.0013484,-0.00115159}}; // Par 0 is P22 e.t.c.
    double SigmaEDepParErr2[6][4] = {{0.0117154,0.0145275,0.00995927,0.0083138},
                                    {0.00659844,0.00815144,0.00555855,0.00462971},
                                    {0.0045234,0.00556333,0.00377952,0.00312741},
                                    {0.00359515,0.00441473,0.00302053,0.00248303},
                                    {0.00273471,0.00336765,0.00231762,0.00190911},
                                    {0.00235227,0.00289976,0.0020076,0.00164492}};


    // Define Energy dependence fits of 3 gaussians + BW to the Sigma values
    TF1 **EDepP2_2 = new TF1*[6];
    double EDepP2FitValsX2[6][21];
    double EDepP2FitValsY2[6][21];
    // Define Energy dependence fits of 3 gaussians + BW to the Sigma values for inclusion of errors
    TF1 **EDepErrP2_2 = new TF1*[6];
    double EDepErrP2FitValsY2[6][21];
    TGraphErrors* EDepP2Graphs2[6];

    for(Int_t i = 0; i < 6; i++){
        sprintf(FitName, "EDepP2%i_2", i+2);
        sprintf(FitErrName, "EDepP2%iErr_2", i+2);
        EDepP2_2[i] = new TF1 (FitName, "[0]*exp(-1*(x-420)*(x-420)/2/100/100)+[1]*exp(-1*(x-520)*(x-520)/2/100/100)+[2]*exp(-1*(x-620)*(x-620)/2/100/100)+[3]*70*70/4/((x-570)*(x-570)+70*70/4)", 410, 620);
        EDepErrP2_2[i] = new TF1 (FitErrName, "[0]*exp(-1*(x-420)*(x-420)/2/100/100)+[1]*exp(-1*(x-520)*(x-520)/2/100/100)+[2]*exp(-1*(x-620)*(x-620)/2/100/100)+[3]*70*70/4/((x-570)*(x-570)+70*70/4)", 410, 620);
        for(Int_t j =0 ; j < 4; j++){
            EDepP2_2[i]->SetParameter(j, SigmaEDepPar2[i][j]);
            EDepErrP2_2[i]->SetParameter(j, SigmaEDepPar2[i][j] + 0.5*SigmaEDepParErr2[i][j]);
        }
        for(Int_t k = 0; k < 21; k++){
            EDepP2FitValsX2[i][k] = 415+(k*10);
            EDepP2FitValsY2[i][k] = EDepP2_2[i]->Eval(415+(k*10));
            EDepErrP2FitValsY2[i][k] = fabs(EDepP2_2[i]->Eval(415+(k*10)) - EDepErrP2_2[i]->Eval(415+(k*15)))/2;
        }
    }

    for(int i = 0; i < 6; i++){
        sprintf(FitGraphName, "EDepP2%i_Graph2", i+2);
        EDepP2Graphs2[i] = new TGraphErrors(21 , EDepP2FitValsX2[i], EDepP2FitValsY2[i], ex2, EDepErrP2FitValsY2[i]);
        EDepP2Graphs2[i]->SetFillColor(2);
        EDepP2Graphs2[i]->SetFillStyle(3005);
        EDepP2Graphs2[i]->SetName(FitGraphName);
        EDepP2Graphs2[i]->Write();
    }

    double SigmaEDepPar3[6][3] = {{-0.163267,0.0763012,0.0745446},
                                    {0.00272942,0.02196,0.0305643},
                                    {0.00109954,-0.000781213,0.0044306},
                                    {0.00319168,0.00246247,0.00358099},
                                    {0.00232062,-0.000615963,0.00386589},
                                    {0.00109522,0.000139162,0.0013484}}; // Par 0 is P22 e.t.c.
    double SigmaEDepParErr3[6][3] = {{0.00249515,0.00171847,0.00721249},
                                    {0.00137007,0.000936932,1.41421},
                                    {0.000924276,0.000626941,1.41421},
                                    {0.00071391,0.00048596,1.41421},
                                    {0.000561895,0.000383696,1.41421},
                                    {0.000466879,0.000318795,1.41421}};

    // Define Energy dependence fits of 3 gaussiansto the Sigma values
    TF1 **EDepP2_3 = new TF1*[6];
    double EDepP2FitValsX3[6][21];
    double EDepP2FitValsY3[6][21];
    // Define Energy dependence fits of 3 gaussians to the Sigma values for inclusion of errors
    TF1 **EDepErrP2_3 = new TF1*[6];
    double EDepErrP2FitValsY3[6][21];
    TGraphErrors* EDepP2Graphs3[6];

    for(Int_t i = 0; i < 6; i++){
        sprintf(FitName, "EDepP2%i_3", i+2);
        sprintf(FitErrName, "EDepP2%iErr_3", i+2);
        EDepP2_3[i] = new TF1 (FitName, "[0]*exp(-1*(x-420)*(x-420)/2/100/100)+[1]*exp(-1*(x-520)*(x-520)/2/100/100)+[2]*exp(-1*(x-620)*(x-620)/2/100/100)", 410, 620);
        EDepErrP2_3[i] = new TF1 (FitErrName, "[0]*exp(-1*(x-420)*(x-420)/2/100/100)+[1]*exp(-1*(x-520)*(x-520)/2/100/100)+[2]*exp(-1*(x-620)*(x-620)/2/100/100)", 410, 620);
        EDepP2_3[i]->SetLineStyle(9);
        for(Int_t j =0 ; j < 3; j++){
            EDepP2_3[i]->SetParameter(j, SigmaEDepPar2[i][j]);
            EDepErrP2_3[i]->SetParameter(j, SigmaEDepPar2[i][j] + 0.5*SigmaEDepParErr2[i][j]);
        }
        for(Int_t k = 0; k < 21; k++){
            EDepP2FitValsX3[i][k] = 415+(k*10);
            EDepP2FitValsY3[i][k] = EDepP2_3[i]->Eval(415+(k*10));
            EDepErrP2FitValsY3[i][k] = fabs(EDepP2_3[i]->Eval(415+(k*10)) - EDepErrP2_3[i]->Eval(415+(k*15)))/2;
        }
    }

    for(int i = 0; i < 6; i++){
        sprintf(FitGraphName, "EDepP2%i_Graph3", i+2);
        EDepP2Graphs3[i] = new TGraphErrors(21 , EDepP2FitValsX3[i], EDepP2FitValsY3[i], ex2, 0);
        EDepP2Graphs3[i]->SetFillColor(2);
        EDepP2Graphs3[i]->SetLineStyle(9);
        EDepP2Graphs3[i]->SetName(FitGraphName);
        EDepP2Graphs3[i]->Write();
    }

    TCanvas *canvas26 = new TCanvas("canvas26","canvas26", 2560, 1440);
    canvas26->Divide(3,2, 0.0000001, 0.001);
    for(int i = 1; i < 7; i++){
        canvas26->cd(i);
        ParameterPlots[i-1]->Draw("AEP");
        EDepP2_2[i-1]->Draw("SAME");
    }
    canvas26->Write();

    TCanvas *canvas27 = new TCanvas("canvas27","canvas27", 2560, 1440);
    canvas27->Divide(3,2, 0.0000001, 0.001);
    for(int i = 1; i < 7; i++){
        canvas27->cd(i);
        ParameterPlots[i-1]->Draw("AEP");
        EDepP2_2[i-1]->Draw("SAME");
        EDepP2Graphs2[i-1]->Draw("SAMEE3");
    }
    canvas27->Write();

    TCanvas *canvas27a = new TCanvas("canvas27a","canvas27a", 2560, 1440);
    canvas27a->Divide(3,2, 0.0000001, 0.001);
    for(int i = 1; i < 7; i++){
        canvas27a->cd(i);
        ParameterPlots[i-1]->Draw("AEP");
        EDepP2_2[i-1]->Draw("SAME");
        EDepP2_3[i-1]->Draw("SAME");
        EDepP2Graphs2[i-1]->Draw("SAMEE3");
    }
    canvas27a->Write();

    TF1 **LegPolEDep1 = new TF1*[21];
    TF1 **LegPolEDep2 = new TF1*[21];
    for(Int_t i = 0; i < 21; i++){
        LegPolEDep1[i] = new TF1(Form("LegPolEDep1_%i", 415+(i*10)), "(1-x**2)*([0]*3+[1]*15*x+[2]*15.0/2*(7*x**2-1)+[3]*105.0/2*x*(3*x**2-1)+[4]*105.0/8*(33*x**4-18*x**2+1)+[5]*63.0/8*x*(143*x**4-110*x**2+15)+[6]*315.0/16*(143*x**6-143*x**4+33*x**2-1)+[7]*495.0/16*(221*x**7-273*x**5+91*x**3-7*x))", -1, 1);
        LegPolEDep1[i]->SetLineColor(807);
        LegPolEDep2[i] = new TF1(Form("LegPolEDep2_%i", 415+(i*10)), "(1-x**2)*([0]*3+[1]*15*x+[2]*15.0/2*(7*x**2-1)+[3]*105.0/2*x*(3*x**2-1)+[4]*105.0/8*(33*x**4-18*x**2+1)+[5]*63.0/8*x*(143*x**4-110*x**2+15)+[6]*315.0/16*(143*x**6-143*x**4+33*x**2-1)+[7]*495.0/16*(221*x**7-273*x**5+91*x**3-7*x))", -1, 1);
        LegPolEDep2[i]->SetLineColor(1);
        for(Int_t j = 0; j < 6; j++){
            LegPolEDep1[i]->SetParameter(j, EDepP2[j]->Eval(415+(i*10)));
            LegPolEDep2[i]->SetParameter(j, EDepP2_2[j]->Eval(415+(i*10)));
        }
        LegPolEDep1[i]->SetParameter(6, 0);
        LegPolEDep2[i]->SetParameter(6, 0);
        LegPolEDep1[i]->SetParameter(7, 0);
        LegPolEDep2[i]->SetParameter(7, 0);
    }

    TCanvas *canvas28 = new TCanvas("canvas28","canvas28", 2560, 1440);
    canvas28->Divide(5,4, 0.0000001, 0.001);
    for(int i = 0 ; i < 21 ; i++){
        if (i < 14){
            canvas28->cd(i+1);
            SigmaPlots[i]->Draw("AEP");
            LegPolEDep1[i]->Draw("SAME");
            LegPolEDep2[i]->Draw("SAME");
        }
        if (i > 14){
            canvas28->cd(i);
            SigmaPlots[i]->Draw("AEP");
            LegPolEDep1[i]->Draw("SAME");
            LegPolEDep2[i]->Draw("SAME");
        }
    }

    canvas28->Write();

    leg4 = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg4->AddEntry(SigmaPlots[7], "E Indep Fit", "l");
    leg4->AddEntry(LegPolEDep1[7], "E Dep Fit (7 Gaus)", "l");
    leg4->AddEntry(LegPolEDep2[7], "E Dep Fit (3 Gaus + BW)", "l");

    TCanvas *canvas29 = new TCanvas("canvas29","canvas29", 2560, 1440);
    SigmaPlots[7]->Draw("AEP");
    LegPolEDep1[7]->Draw("SAME");
    LegPolEDep2[7]->Draw("SAME");
    leg4->Draw("SAME");
    canvas29->Write();

    f5.Write();
     // Save Sigma Values to .dat file
    ofstream outfile1("Sigma_NS18.dat");

    Int_t Dummy1=0;
    for(Int_t i = 0;i < 21;i++){
        for(Int_t j =0 ; j<18; j++){
            if(i != 14){
            outfile1 << Dummy1 <<"   "<< 415+(i*10) <<"   "<< 0.85 - (j*0.1) <<"   "<< SigmaValues[i][j] <<"   "<< SigmaErrValues[i][j] <<endl;
            Dummy1+=1;
            }
        }
    }

}
