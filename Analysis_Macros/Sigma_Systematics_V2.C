#include "./includes_Sigma_Systematics_V2.h"

void Sigma_Systematics_V2(){

    char HistName[60];
    char NewHistName[60];
    TH1F* PhipPara[21][20];
    TH1F* PhipPerp[21][20];
    TH1F* PhipPara2[21][20];
    TH1F* PhipPerp2[21][20];
    double EStart = 415; // What is the initial EGamma bin? Change this depending on value

    TF1 *CosFunc = new TF1("CosFit", "[0]*cos((2*x*TMath::DegToRad())+acos(0))");
    CosFunc->SetParNames("Amplitude");

    TH1F* AsymmHists[21][20];
    TH1F* AsymmHists2[21][20];

    char ParaHistName[60];
    char PerpHistName[60];
    char AsymmHistName[60];
    char name[21];
    char title[60];

    TGraphErrors* SigmaSystPlots[21];

    TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/Para/NoScatt/Physics_Total_Para_NoScatt_18_26_4_18.root"); // Open latest Para file
    TFile *f1 = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/Perp/NoScatt/Physics_Total_Perp_NoScatt_18_26_4_18.root"); // Open latest Para file

    TH1D* Eg_Para = (TH1D*)f->Get("Eg2")->Clone();
    Eg_Para->SetName("Eg_Para");

    TH1D* Eg_Perp = (TH1D*)f1->Get("Eg2")->Clone();
    Eg_Perp->SetName("Eg_Perp");

    NPara = Eg_Para->GetEntries();
    NPerp = Eg_Perp->GetEntries();
    ScaleFactor = NPara/NPerp;
    ScaleFactorErr = sqrt( (NPara/((TMath::Power(NPerp,2)))) + (((TMath::Power(NPara,2)))/(TMath::Power(NPerp,3))) ); // Error Propagation of above formula, see notebook 5

    for(Int_t i = 0; i < 21; i++){ // Energy
        for(Int_t j = 0; j < 20; j++){ // Theta
            sprintf(ParaHistName, "Phip_%iMeVCM%i", EStart+(i*10) , j+1);
            sprintf(PerpHistName, "Phip_%iMeVCM%i", EStart+(i*10) , j+1);
            sprintf(AsymmHistName, "Sigma_%iMeV_CM%i_Hist", EStart+(i*10) , j+1);
            AsymmHists[i][j] = (TH1F*) (((TH1F*)f->Get(ParaHistName))->GetAsymmetry(((TH1F*)f1->Get(PerpHistName)), ScaleFactor, ScaleFactorErr)));
            AsymmHists[i][j]->SetName(AsymmHistName);
            AsymmHists[i][j]-> Fit("CosFit", "Q");
            cout << "NDOF " << CosFit->GetNDF() << "   " << "Chi2 " << CosFit->GetChisquare() << "   " << "Chi2/DoF " << CosFit->GetChisquare()/CosFit->GetNDF() << endl;
            pCosAmp[i][j] = CosFit->GetParameter(0);
            pCosAmpErr[i][j] = CosFit->GetParError(0);
        }
    }

    TFile *f2 = TFile::Open("/scratch/Mainz_Software/a2GoAT/LinPol_Aug16.root"); // Open linear polarisation plot

    for(Int_t i = 0; i < 21; i++){
        for(Int_t j = 0; j < 18; j++){
            Double_t EValue = EStart + (i*10);
            Double_t PolVal = Graph->Eval(EValue, 0);
            Double_t PolErrVal = (3./100.)*PolVal; // 3% error on Polarisation assumed here
            SigmaValues[0][i][j] = pCosAmp[i][j+1]/(PolVal);
            SigmaErrValues[0][i][j] = sqrt(((pCosAmpErr[i][j+1]/(PolVal))**2) + (((PolErrVal*pCosAmp[i][j+1])/(PolVal**2))**2));
        }
    }

    /////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////
    ////////////////////////// 1Sig File Done ///////////////////////
    /////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////


    TFile *f3 = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/Para/NoScatt/Physics_Total_Para_NoScatt_17_26_4_18.root"); // Open latest Para file
    TFile *f4 = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/Perp/NoScatt/Physics_Total_Perp_NoScatt_17_26_4_18.root"); // Open latest Para file

    TH1D* Eg_Para2 = (TH1D*)f3->Get("Eg2")->Clone();
    Eg_Para2->SetName("Eg_Para2");

    TH1D* Eg_Perp2 = (TH1D*)f4->Get("Eg2")->Clone();
    Eg_Perp2->SetName("Eg_Perp2");

    NPara2 = Eg_Para2->GetEntries();
    NPerp2 = Eg_Perp2->GetEntries();
    ScaleFactor2 = NPara2/NPerp2;
    ScaleFactorErr2 = sqrt( (NPara2/((TMath::Power(NPerp2,2)))) + (((TMath::Power(NPara2,2)))/(TMath::Power(NPerp2,3))) ); // Error Propagation of above formula, see notebook 5

    for(Int_t i = 0; i < 21; i++){ // Energy
        for(Int_t j = 0; j < 20; j++){ // Theta
            sprintf(ParaHistName, "Phip_%iMeVCM%i", EStart+(i*10) , j+1);
            sprintf(PerpHistName, "Phip_%iMeVCM%i", EStart+(i*10) , j+1);
            sprintf(AsymmHistName, "Sigma_%iMeV_CM%i_Hist", EStart+(i*10) , j+1);
            AsymmHists2[i][j] = (TH1F*) (((TH1F*)f3->Get(ParaHistName))->GetAsymmetry(((TH1F*)f4->Get(PerpHistName)), ScaleFactor2, ScaleFactorErr2)));
            AsymmHists2[i][j]->SetName(AsymmHistName);
            AsymmHists2[i][j]-> Fit("CosFit", "Q");
            cout << "NDOF " << CosFit->GetNDF() << "   " << "Chi2 " << CosFit->GetChisquare() << "   " << "Chi2/DoF " << CosFit->GetChisquare()/CosFit->GetNDF() << endl;
            pCosAmp2[i][j] = CosFit->GetParameter(0);
            pCosAmpErr2[i][j] = CosFit->GetParError(0);
        }
    }

    TFile *f5 = TFile::Open("/scratch/Mainz_Software/a2GoAT/LinPol_Aug16.root"); // Open linear polarisation plot

    for(Int_t i = 0; i < 21; i++){
        for(Int_t j = 0; j < 18; j++){
            Double_t EValue = EStart + (i*10);
            Double_t PolVal = Graph->Eval(EValue, 0);
            Double_t PolErrVal = (3./100.)*PolVal; // 3% error on Polarisation assumed here
            SigmaValues[1][i][j] = pCosAmp2[i][j+1]/(PolVal);
            SigmaErrValues[1][i][j] = sqrt(((pCosAmpErr2[i][j+1]/(PolVal))**2) + (((PolErrVal*pCosAmp2[i][j+1])/(PolVal**2))**2));
        }
    }

    TFile f6("Sigma_Systematic_18_17_V2.root", "RECREATE");

    TF1 *Line = new TF1("Line", "[0]", -1, 1);
    double x[18] = {0.85, 0.75, 0.65, 0.55, 0.45, 0.35, 0.25, 0.15, 0.05, -0.05, -0.15, -0.25, -0.35, -0.45, -0.55, -0.65, -0.75, -0.85}; // Need to adjust
    double ex[18] = {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05}; // Need to adjust

    for(Int_t i = 0; i < 21; i ++){ // Egamma value
        for(Int_t j = 0; j < 18; j ++){ // CM value
            SystValues[i][j] = SigmaValues[0][i][j] - SigmaValues[1][i][j];
            SystErrValues[i][j] = sqrt( ((SigmaErrValues[0][i][j])**2) + ((SigmaErrValues[1][i][j])**2) );
        }

        sprintf(title, "#Sigma_{2#sigma} - #Sigma_{1#sigma} E_{#gamma} %i #pm 10 MeV", EStart+(i*10));
        SigmaSystPlots[i] = new TGraphErrors(18 , x, SystValues[i] , ex, SystErrValues[i]) ;
        SigmaSystPlots[i]->Fit("Line", "M");
        SigmaSystPlots[i]->SetMarkerColor(4);
        SigmaSystPlots[i]->SetLineColor(4);
        SigmaSystPlots[i]->SetMarkerStyle(8);
        SigmaSystPlots[i]->SetMarkerSize(1);
        SigmaSystPlots[i]->GetXaxis()->SetTitle("Cos#theta_{CM}");
        SigmaSystPlots[i]->GetXaxis()->SetRangeUser(-1, 1);
        SigmaSystPlots[i]->GetYaxis()->SetTitle("#Sigma_{2#sigma} - #Sigma_{1#sigma}");
        SigmaSystPlots[i]->SetName(name);
        SigmaSystPlots[i]->SetTitle(title);
        SigmaSystPlots[i]->Write();

        LinePar[i] = Line->GetParameter(0);
        LineParErr[i] = Line->GetParError(0);

        cout << (LinePar[i])/(LineParErr[i]) << endl;
    }

    TCanvas *canvas20 = new TCanvas("canvas20","canvas20", 1920, 1080);
    canvas20->Divide(5,4);
    for(int i = 1 ; i < 22 ; i++){
        if(i == 15) continue;
        if (i < 15) canvas20->cd(i);
        else if (i > 15) canvas20->cd(i-1);
        SigmaSystPlots[i-1]->Draw("AEP");
    }

    canvas20->Write();

    f6.Write();

}

