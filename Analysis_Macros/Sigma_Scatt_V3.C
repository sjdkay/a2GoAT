#include "./includes_Sigma_Scatt_V3.h"

void Sigma_Scatt_V3(){

    char HistName[60];
    char NewHistName[60];
    TH1F* PhipPara[10][5];
    TH1F* PhipPerp[10][5];
    double EStart = 415; // What is the initial EGamma bin? Change this depending on value

    TF1 *CosFunc = new TF1("CosFit", "[0]*cos((2*x*TMath::DegToRad())+acos(0))");
    CosFunc->SetParNames("Amplitude");

    TH1F* AsymmHists[10][5];

    char ParaHistName[60];
    char PerpHistName[60];
    char AsymmHistName[60];
    char GraphName[60];
    char ScGraphName[60];
    char ScGraphNameAdj[60];

    TF1 *LegPol = new TF1("LegPol", "(1-x**2)*([0]*3+[1]*15*x+[2]*15.0/2*(7*x**2-1)+[3]*105.0/2*x*(3*x**2-1)+[4]*105.0/8*(33*x**4-18*x**2+1)+[5]*63.0/8*x*(143*x**4-110*x**2+15)+[6]*315.0/16*(143*x**6-143*x**4+33*x**2-1)+[7]*495.0/16*(221*x**7-273*x**5+91*x**3-7*x))", -1, 1);
    LegPol->SetLineColor(4);
    LegPol->SetLineWidth(2);
    LegPol->SetParLimits(6, 0.0, 0.0);
    LegPol->SetParLimits(7, 0.0, 0.0);
    LegPol->FixParameter(6, 0.0); // These seem to be ignored?
    LegPol->FixParameter(7, 0.0);

    TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/Para/Physics_Total_Para_43_25_5_18.root"); // Open latest Para file

    TH1D* Eg_Para = (TH1D*)f->Get("Eg")->Clone();
    Eg_Para->SetName("Eg_Para");

    for(Int_t i = 0; i < 10; i++){ // Energy
        for(Int_t j = 0; j < 5; j++){ // Theta
            sprintf(HistName, "Phip_%iMeVCM%i", EStart+(i*20) , j+1);
            sprintf(NewHistName, "Phip_%iMeVCM%i_Para", EStart+(i*20) , j+1);
            PhipPara[i][j] = ((TH1F*)f->Get(HistName));
            PhipPara[i][j]->SetName(NewHistName);
        }
    }


    /////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////
    ////////////////////////// Para Done ////////////////////////////
    /////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////

    TFile *f1 = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/Perp/Physics_Total_Perp_43_25_5_18.root"); // Open latest Para file

    TH1D* Eg_Perp = (TH1D*)f1->Get("Eg")->Clone();
    Eg_Perp->SetName("Eg_Perp");

    for(Int_t i = 0; i < 10; i++){ // Energy
        for(Int_t j = 0; j < 5; j++){ // Theta
            sprintf(HistName, "Phip_%iMeVCM%i", EStart+(i*20) , j+1);
            sprintf(NewHistName, "Phip_%iMeVCM%i_Perp", EStart+(i*20) , j+1);
            PhipPerp[i][j] = ((TH1F*)f1->Get(HistName));
            PhipPerp[i][j]->SetName(NewHistName);
        }
    }

    TFile f2("ParaPerp_S43_Combined.root", "RECREATE");

    Eg_Para->Write();
    Eg_Perp->Write();

    for(Int_t i = 0; i < 10; i++){ // Energy
        for(Int_t j = 0; j < 5; j++){ // Theta
            PhipPara[i][j]->Write();
            PhipPerp[i][j]->Write();
        }
    }

    f2.Write();

    NPara = Eg_Para->GetEntries();
    NPerp = Eg_Perp->GetEntries();
    ScaleFactor = NPara/NPerp;
    ScaleFactorErr = sqrt( (NPara/((TMath::Power(NPerp,2)))) + (((TMath::Power(NPara,2)))/(TMath::Power(NPerp,3))) ); // Error Propagation of above formula, see notebook 5

    for(Int_t i = 0; i < 10; i++){ // Energy
        for(Int_t j = 0; j < 5; j++){ // Theta
            sprintf(ParaHistName, "Phip_%iMeVCM%i_Para", EStart+(i*20) , j+1);
            sprintf(PerpHistName, "Phip_%iMeVCM%i_Perp", EStart+(i*20) , j+1);
            sprintf(AsymmHistName, "Sigma_%iMeV_CM%i_Hist", EStart+(i*20) , j+1);
            AsymmHists[i][j] = (TH1F*) (((TH1F*)f2.Get(ParaHistName))->GetAsymmetry(((TH1F*)f2.Get(PerpHistName)), ScaleFactor, ScaleFactorErr)));
            AsymmHists[i][j]->SetName(AsymmHistName);
            AsymmHists[i][j]-> Fit("CosFit", "Q");
            cout << "NDOF " << CosFit->GetNDF() << "   " << "Chi2 " << CosFit->GetChisquare() << "   " << "Chi2/DoF " << CosFit->GetChisquare()/CosFit->GetNDF() << endl;
            pCosAmp[i][j] = CosFit->GetParameter(0);
            pCosAmpErr[i][j] = CosFit->GetParError(0);
        }
    }

    TFile f3("ParaPerpAsymm_S43.root", "RECREATE");

    for(Int_t i = 0; i < 10; i++){ // Energy
        for(Int_t j = 0; j < 5; j++){ // Theta
            AsymmHists[i][j]->Write();
        }
    }

    f3.Write();

    TGraphErrors* SigmaPlots[10];
    char name[21];
    char title[60];
    char name2[60];
    char title2[60];

    TFile *f4= TFile::Open("/scratch/Mainz_Software/a2GoAT/LinPol_Aug16.root"); // Open linear polarisation plot

    // Calculate values of sigma for each angular and energy bin

    for(Int_t i = 0; i < 10; i++){
        for(Int_t j = 0; j < 5; j++){
            Double_t EValue = EStart + (i*10);
            Double_t PolVal = Graph->Eval(EValue, 0);
            Double_t PolErrVal = (3./100.)*PolVal; // 3% error on Polarisation assumed here
            SigmaValues[i][j] = pCosAmp[i][j+1]/(PolVal);
            SigmaErrValues[i][j] = sqrt(((pCosAmpErr[i][j+1]/(PolVal))**2) + (((PolErrVal*pCosAmp[i][j+1])/(PolVal**2))**2));
        }
    }


    TFile f5("Sigma_Plots_S43.root", "RECREATE");

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
    canvas20->Divide(4,3);
    for(int i = 1 ; i < 11 ; i++){
        canvas20->cd(i);
        SigmaPlots[i-1]->Draw("AEP");
    }

    canvas20->Write();

    f5.Write();
}
