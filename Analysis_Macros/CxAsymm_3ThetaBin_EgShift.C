#include "./includes.h"

// define a function with 4 parameters
Double_t fitf(Double_t *x,Double_t *par)
{
    Double_t fitval = 0;
    fitval =  ((par[0]*sin(x[0]))/(1 + (par[1]*cos(x[0]))));
    return fitval;
}

void CxAsymm_3ThetaBin_EgShift() {

    gStyle->SetTitleFontSize(0.06);
    gStyle->SetOptStat(0);
    char FileName[200];
    char PosHelHistName[60];
    char NegHelHistName[60];
    char AsymmHistName[60];
    char NeutronEThetaScName[60];
    char name[60];
    char title[60];
    char name2[60];
    char title2[60];
    char name3[60];
    char title3[60];
    char name4[60];
    char title4[60];
    char AsymmHistTitle[60];

    TH1F* AsymmHists[8][3][8][3];
    TGraphErrors* CxPlots[8][3][8]; // As fn of Angular bin for fixed energies
    TGraphErrors* Cx2Plots[8][3][3]; // As fn of energy bin for fixed angular bins
    TGraphErrors* pyPlots[8][8]; // As fn of Angular bin for fixed energies
    TGraphErrors* py2Plots[8][8]; // As fn of Angular bin for fixed energies

    double SinAmps[8][3][8][3]; // To store and get Cx amps
    double SinAmpErrs[8][3][8][3];// To store and get Cx amps
    double CosAmps[8][8][3];
    double CosAmpErrs[8][8][3];
    double MeanX[8][8][3];
    double MeanY[8][8][3];

    Double_t Cx[8][3][8][3]; // Fixed P Chi2 fit, Fixed P LL fit, Variable p Chi2 fit, Variable P LL fit
    Double_t CxErr[8][3][8][3];
    Double_t py[8][8][3];
    Double_t pyErr[8][8][3];

    Double_t Cx2[8][3][8]; // Fixed P Chi2 fit, Fixed P LL fit, Variable p Chi2 fit, Variable P LL fit
    Double_t Cx2Err[8][3][8];

    Double_t py2[8][3][8];
    Double_t py2Err[8][3][8];

    Double_t PVal;
    Double_t ECentre;
    Double_t EPoint;
    Double_t AmpVal;
    Double_t AEff[8][8][3];
    Double_t CorrFac[8][8][3];
    Double_t AEffCorr[8][8][3];

    TF1 *AsymmFunc = new TF1("AsymmFit",  fitf, -3.0, 3.0, 2); //Give a name and range to the fitting funcion
    AsymmFunc->SetParNames("SinAmp", "CosAmp"); //Name the parameters
    AsymmFunc->SetParameter(0, 0);
    TFile *f = new TFile("/scratch/Mainz_Software/a2GoAT/Physics_Total_Amo151_Lin51_Combined_EgShift.root"); // Open the latest PTotal file to load histograms from
    TFile *fAy = new TFile ("/scratch/Mainz_Software/a2GoAT/npAy.root");
    TF1 *Pn90CM = new TF1("Pn90CM", "1.64576-2.95484*(x/1000)+0.684577*(x/1000)**2-0.65*90**2/4/((x-560)**2+90**2/4)+(5.32305-35.3819*(x/1000)+70.145*(x/1000)**2-44.2899*(x/1000)**3)",200,1000);

    for(Int_t l = 0; l < 8; l++){ // EgShift
        for(Int_t i = 0; i < 3; i++){ // Fit version
            for(Int_t j = 0; j < 8; j++){ // Energy
                if (l < 4) {ECentre = 260 + (j*100)+(l*10);}
                else{ECentre = 270 + (j*100)+(l*10);} // Get centre of energy bin
                PVal = (Pn90CM->Eval(ECentre, 0)); // Get PValue for energy bin to fix parameter with
                for(Int_t k = 0; k < 3; k++)){ // CosThetaBin

                    if (l < 4){
                        sprintf(PosHelHistName, "PhiSc%iPosHelCM%i", 260 + (j*100)+(l*10), k+1);
                        sprintf(NegHelHistName, "PhiSc%iNegHelCM%i", 260 + (j*100)+(l*10), k+1);
                        sprintf(AsymmHistName, "CxAsymm%iCM%iFit%i", 260 + (j*100)+(l*10), k+1, i+1);
                        sprintf(AsymmHistTitle, "#phi_{Sc} ASymmetry (%i #pm 50) MeV CM%i", 260 + (j*100)+(l*10), k+1);
                        sprintf(NeutronEThetaScName, "NeutronEThetaSc%iCM%i", 260 + (j*100)+(l*10), k+1);
                    }

                    else{
                        sprintf(PosHelHistName, "PhiSc%iPosHelCM%i", 270 + (j*100)+(l*10), k+1);
                        sprintf(NegHelHistName, "PhiSc%iNegHelCM%i", 270 + (j*100)+(l*10), k+1);
                        sprintf(AsymmHistName, "CxAsymm%iCM%iFit%i", 270 + (j*100)+(l*10), k+1, i+1);
                        sprintf(AsymmHistTitle, "#phi_{Sc} ASymmetry (%i #pm 50) MeV CM%i", 270 + (j*100)+(l*10), k+1);
                        sprintf(NeutronEThetaScName, "NeutronEThetaSc%iCM%i", 270 + (j*100)+(l*10), k+1);
                    }

                    AsymmHists[l][i][j][k] = (TH1F*) (((TH1F*)f->Get(PosHelHistName))->GetAsymmetry(((TH1F*)f->Get(NegHelHistName)))));
                    AsymmHists[l][i][j][k]->SetName(AsymmHistName);
                    AsymmHists[l][i][j][k]->SetTitle(AsymmHistTitle);
                    AsymmHists[l][i][j][k]->GetXaxis()->SetTitle("#phi_{Sc}/Rad");
                    AsymmHists[l][i][j][k]->GetYaxis()->SetTitle("Asymmetry");
                    AsymmHists[l][i][j][k]->GetXaxis()->SetLabelSize(0.06);
                    AsymmHists[l][i][j][k]->GetXaxis()->SetTitleSize(0.06);
                    AsymmHists[l][i][j][k]->GetXaxis()->SetTitleOffset(0.7);
                    AsymmHists[l][i][j][k]->GetXaxis()->CenterTitle();
                    AsymmHists[l][i][j][k]->GetYaxis()->SetLabelSize(0.06);
                    AsymmHists[l][i][j][k]->GetYaxis()->SetTitleSize(0.06);
                    AsymmHists[l][i][j][k]->GetYaxis()->SetTitleOffset(0.7);
                    AsymmHists[l][i][j][k]->GetYaxis()->CenterTitle();

                    MeanX[l][j][k] = ((TH2D*)f->Get(NeutronEThetaScName))->GetMean(1);
                    MeanY[l][j][k] = ((TH2D*)f->Get(NeutronEThetaScName))->GetMean(2);

                    AEff[l][j][k] = ((TH2F*)fAy->Get("nppnAy"))->Interpolate(MeanX[l][j][k], MeanY[l][j][k]);
                    CorrFac[l][j][k] = (1+exp(1.81572-(0.0139530*MeanX[l][j][k])));
                    AEffCorr[l][j][k] = AEff[l][j][k] * CorrFac[l][j][k]; //Analysing power for 12C based on correction factor from Mikhail
                    AmpVal = PVal*AEff[l][j][k];

                    if (i == 0){
                        AsymmFit->SetLineColor(4);
                        AsymmFit->FixParameter(1, 0.);
                        AsymmHists[l][i][j][k]->Fit("AsymmFit", "Q");
                        AsymmHists[l][i][j][k]->SetLineColor(4);
                        SinAmps[l][i][j][k] = AsymmFit->GetParameter(0);
                        SinAmpErrs[l][i][j][k] = AsymmFit->GetParError(0);
                    }

                    if(i == 1){
                        AsymmFit->SetLineColor(2);
                        AsymmFit->FixParameter(1, AmpVal);
                        AsymmHists[l][i][j][k]->SetLineColor(2);
                        AsymmHists[l][i][j][k]->Fit("AsymmFit", "Q");
                        SinAmps[l][i][j][k] = AsymmFit->GetParameter(0);
                        SinAmpErrs[l][i][j][k] = AsymmFit->GetParError(0);
                    }

                    if(i == 2){
                        AsymmFit->SetLineColor(1);
                        AsymmFit->ReleaseParameter(1);
                        AsymmHists[l][i][j][k]->SetLineColor(1);
                        AsymmFit->SetParLimits(1, -1*fabs(AEff[l][j][k]), fabs(AEff[l][j][k]));
                        AsymmHists[l][i][j][k]->Fit("AsymmFit", "Q");
                        SinAmps[l][i][j][k] = AsymmFit->GetParameter(0);
                        SinAmpErrs[l][i][j][k] = AsymmFit->GetParError(0);
                        CosAmps[l][j][k] = AsymmFit->GetParameter(1);
                        CosAmpErrs[l][j][k] = AsymmFit->GetParError(1);
                    }
                }
            }
        }
    }

    TFile f1("AsymmFits_PTotal_151_51_EgShift_V1.root", "RECREATE");

    for(Int_t l = 0; l < 8; l++){
        for(Int_t i = 0; i < 3; i++){ // Fit version
            for(Int_t j = 0; j < 8; j++){ // Energy
                for(Int_t k = 0; k < 3; k++)){
                    AsymmHists[l][i][j][k]->Write();
                }
            }
        }
    }

    f1.Write();
    f1.Close();

    TFile *f2= TFile::Open("/scratch/Mainz_Software/a2GoAT/CircPol_Aug16.root");

    for(Int_t k = 0; k < 8; k++){
        for (Int_t m = 0; m < 3; m++){ // Fit Version
            for (Int_t n = 0; n < 8; n++){ //E
                if (k < 4) {EPoint = 260 + (n*100)+(k*10);}
                else{EPoint = 270 + (n*100)+(k*10);} // Get centre of energy bin
                for(Int_t l = 0; l < 3; l++)){ // Theta CM bin
                    Cx[k][m][n][l] = SinAmps[k][m][n][l]/(AEff[k][n][l]*(Graph->Eval(EPoint ,0)));
                    CxErr[k][m][n][l] = fabs(SinAmpErrs[k][m][n][l]/(AEff[k][n][l]*(Graph->Eval(EPoint ,0))));
                    if(m == 2){
                        py[k][n][l] = CosAmps[k][n][l]/AEff[k][n][l];
                        pyErr[k][n][l] = fabs(CosAmpErrs[k][n][l]/AEff[k][n][l]);
                    }
                }
            }
        }
    }

    TFile f3("Cx_Plots_151_51_V1.root", "RECREATE");

    double x[3] = {2./3., 0., -2./3.};
    double ex[3] = {1./3., 1./3., 1./3.};

    for(Int_t h = 0; h < 8; h++){
        for(Int_t i = 0; i < 3; i++){
            for(Int_t j = 0; j < 8; j++){

                if (h < 4){
                    sprintf(name, "Cx_%i_V%i_3Theta", 260 + (j*100)+(h*10), i+1);
                    sprintf(title, "C_{x'}(Cos#theta_{CM}) E_{#gamma} %i #pm 50 MeV", 260 + (j*100)+(h*10), i+1);
                }

                else{
                    sprintf(name, "Cx_%i_V%i_3Theta", 270 + (j*100)+(h*10), i+1);
                    sprintf(title, "C_{x'}(Cos#theta_{CM}) E_{#gamma} %i #pm 50 MeV", 270 + (j*100)+(h*10), i+1);
                }

                CxPlots[h][i][j] = new TGraphErrors(3 , x, Cx[h][i][j], ex, CxErr[h][i][j]);
                CxPlots[h][i][j]->SetName(name);
                CxPlots[h][i][j]->SetTitle(title);
                CxPlots[h][i][j]->SetMarkerSize(1.2);
                if(i == 0){
                    CxPlots[h][i][j]->SetMarkerColor(4);
                    CxPlots[h][i][j]->SetMarkerStyle(20);
                    CxPlots[h][i][j]->SetMarkerSize(1.5);
                    CxPlots[h][i][j]->SetLineColor(4);
                }

                else if (i == 1){
                    CxPlots[h][i][j]->SetMarkerColor(2);
                    CxPlots[h][i][j]->SetMarkerStyle(21);
                    CxPlots[h][i][j]->SetLineColor(2);
                }

                else if (i == 2){
                    CxPlots[h][i][j]->SetMarkerColor(1);
                    CxPlots[h][i][j]->SetMarkerStyle(22);
                    CxPlots[h][i][j]->SetLineColor(1);
                }
                CxPlots[h][i][j]->GetXaxis()->SetTitle("Cos#theta_{CM}");
                CxPlots[h][i][j]->GetXaxis()->SetRangeUser(-1, 1);
                CxPlots[h][i][j]->GetYaxis()->SetRangeUser(-2, 2);
                CxPlots[h][i][j]->GetYaxis()->SetTitle("C_{x'}");
                CxPlots[h][i][j]->GetXaxis()->SetLabelSize(0.06);
                CxPlots[h][i][j]->GetXaxis()->SetTitleSize(0.06);
                CxPlots[h][i][j]->GetXaxis()->SetTitleOffset(0.7);
                CxPlots[h][i][j]->GetXaxis()->CenterTitle();
                CxPlots[h][i][j]->GetYaxis()->SetLabelSize(0.06);
                CxPlots[h][i][j]->GetYaxis()->SetTitleSize(0.08);
                CxPlots[h][i][j]->GetYaxis()->SetTitleOffset(0.5);
                CxPlots[h][i][j]->GetYaxis()->CenterTitle();
                CxPlots[h][i][j]->Write();
            }
        }
    }

    double x2[8];
    double ex2[8] = {50, 50, 50, 50, 50, 50, 50, 50};

    for(Int_t h = 0; h < 8; h++){
        for(Int_t i = 0; i < 3; i++){
            for(Int_t k = 0; k < 3; k++){
                for(Int_t j = 0; j < 8; j++){
                    Cx2[h][i][j] = Cx[h][i][j][k];
                    Cx2Err[h][i][j] = CxErr[h][i][j][k];
                    if(h < 4){x2[j] = 260 + (j*100)+(h*10);}
                    else{x2[j] = 270 + (j*100)+(h*10);}
                }
                sprintf(title2, "C_{x'}(E_{#gamma}) #theta_{CM} %i - %i", k*60 ,60 + (k*60));
                if (h < 4){sprintf(name2, "Cx_%iMeVCent_CM%i_V%i_3Theta", 60+(h*10), k+1, i+1);}
                else{sprintf(name2, "Cx_%iMeVCent_CM%i_V%i_3Theta", -30+(h*10), k+1, i+1);}
                Cx2Plots[h][i][k] = new TGraphErrors(8 , x2, Cx2[h][i], ex2, Cx2Err[h][i]);
                Cx2Plots[h][i][k]->SetName(name2);
                Cx2Plots[h][i][k]->SetTitle(title2);
                Cx2Plots[h][i][k]->SetMarkerSize(1.2);
                if(i == 0){
                    Cx2Plots[h][i][k]->SetMarkerColor(4);
                    Cx2Plots[h][i][k]->SetMarkerStyle(20);
                    Cx2Plots[h][i][k]->SetMarkerSize(1.5);
                    Cx2Plots[h][i][k]->SetLineColor(4);
                }

                else if (i == 1){
                    Cx2Plots[h][i][k]->SetMarkerColor(2);
                    Cx2Plots[h][i][k]->SetMarkerStyle(21);
                    Cx2Plots[h][i][k]->SetLineColor(2);
                }
                else if (i == 2){
                    Cx2Plots[h][i][k]->SetMarkerColor(1);
                    Cx2Plots[h][i][k]->SetMarkerStyle(22);
                    Cx2Plots[h][i][k]->SetLineColor(1);
                }
                Cx2Plots[h][i][k]->GetXaxis()->SetTitle("E_{#gamma}");
                Cx2Plots[h][i][k]->GetXaxis()->SetRangeUser(100, 1200);
                Cx2Plots[h][i][k]->GetYaxis()->SetRangeUser(-2, 2);
                Cx2Plots[h][i][k]->GetYaxis()->SetTitle("C_{x'}");
                Cx2Plots[h][i][k]->GetXaxis()->SetLabelSize(0.06);
                Cx2Plots[h][i][k]->GetXaxis()->SetTitleSize(0.06);
                Cx2Plots[h][i][k]->GetXaxis()->SetTitleOffset(0.7);
                Cx2Plots[h][i][k]->GetXaxis()->CenterTitle();
                Cx2Plots[h][i][k]->GetYaxis()->SetLabelSize(0.06);
                Cx2Plots[h][i][k]->GetYaxis()->SetTitleSize(0.08);
                Cx2Plots[h][i][k]->GetYaxis()->SetTitleOffset(0.5);
                Cx2Plots[h][i][k]->GetYaxis()->CenterTitle();
                Cx2Plots[h][i][k]->Write();
            }
        }
    }

    leg2 = new TLegend(0.2, 0.2, 0.8, 0.8);
    leg2->AddEntry(Cx2Plots[0][0][0], "Fixed p_{y} = 0", "lepz");
    leg2->AddEntry(Cx2Plots[0][1][0], "Fixed p_{y}", "lepz");
    leg2->AddEntry(Cx2Plots[0][2][0], "Variable p_{y}", "lepz");

    TCanvas *canvas2 = new TCanvas("canvas2","canvas2", 2560, 1440);
    canvas2->Divide(2, 2, 0.000000001, 0.01);

    for(Int_t i = 0; i < 3; i++){
        canvas2->cd(i+1);
        for(Int_t j = 0; j < 8; j++){
            if(j == 0){Cx2Plots[j][0][i]->Draw("AEP");}
            else{Cx2Plots[j][0][i]->Draw("SAMEEP");}
            Cx2Plots[j][1][i]->Draw("SAMEEP");
            Cx2Plots[j][2][i]->Draw("SAMEEP");
        }
    }

    canvas2->cd(4);
    leg2->Draw();
    canvas2->Write();
//
//    for (Int_t n = 0; n < 8; n++){ //E
//        for(Int_t l = 0; l < 3; l++)){ // Theta CM bin
//            py2[l][n] = py[n][l];
//            py2Err[l][n] = pyErr[n][l];
//        }
//    }
//
//    for(Int_t j = 0; j < 8; j++){
//        sprintf(name3, "py_%i", 250+(j*100));
//        sprintf(title3, "p_{y}(Cos#theta_{CM}) E_{#gamma} %i #pm 50 MeV", 250+(j*100));
//        pyPlots[j] = new TGraphErrors(3 , x, py[j], ex, pyErr[j]);
//        pyPlots[j]->SetName(name3);
//        pyPlots[j]->SetTitle(title3);
//        pyPlots[j]->SetMarkerSize(1.2);
//        pyPlots[j]->SetMarkerColor(1);
//        pyPlots[j]->SetMarkerStyle(22);
//        pyPlots[j]->SetMarkerSize(1.2);
//        pyPlots[j]->SetLineColor(1);
//        pyPlots[j]->GetXaxis()->SetTitle("Cos#theta_{CM}");
//        pyPlots[j]->GetXaxis()->SetRangeUser(-1, 1);
//        pyPlots[j]->GetYaxis()->SetRangeUser(-2, 2);
//        pyPlots[j]->GetYaxis()->SetTitle("p_{y}");
//        pyPlots[j]->GetXaxis()->SetLabelSize(0.06);
//        pyPlots[j]->GetXaxis()->SetTitleSize(0.06);
//        pyPlots[j]->GetXaxis()->SetTitleOffset(0.7);
//        pyPlots[j]->GetXaxis()->CenterTitle();
//        pyPlots[j]->GetYaxis()->SetLabelSize(0.06);
//        pyPlots[j]->GetYaxis()->SetTitleSize(0.08);
//        pyPlots[j]->GetYaxis()->SetTitleOffset(0.5);
//        pyPlots[j]->GetYaxis()->CenterTitle();
//        pyPlots[j]->Write();
//    }
//
//    for(Int_t j = 0; j < 3; j++){
//        sprintf(name4, "py_CM%i", j+1);
//        sprintf(title4, "p_{y}(E_{#gamma}) #theta_{CM} %i - %i", j*60 ,60 + (j*60));
//        py2Plots[j] = new TGraphErrors(8 , x2, py2[j], ex2, py2Err[j]);
//        py2Plots[j]->SetName(name4);
//        py2Plots[j]->SetTitle(title4);
//        py2Plots[j]->SetMarkerSize(1.2);
//        py2Plots[j]->SetMarkerColor(1);
//        py2Plots[j]->SetMarkerStyle(22);
//        py2Plots[j]->SetMarkerSize(1.2);
//        py2Plots[j]->SetLineColor(1);
//        py2Plots[j]->GetXaxis()->SetTitle("E_{#gamma}");
//        py2Plots[j]->GetXaxis()->SetRangeUser(200, 1000);
//        py2Plots[j]->GetYaxis()->SetRangeUser(-2, 2);
//        py2Plots[j]->GetYaxis()->SetTitle("p_{y}");
//        py2Plots[j]->GetXaxis()->SetLabelSize(0.06);
//        py2Plots[j]->GetXaxis()->SetTitleSize(0.06);
//        py2Plots[j]->GetXaxis()->SetTitleOffset(0.7);
//        py2Plots[j]->GetXaxis()->CenterTitle();
//        py2Plots[j]->GetYaxis()->SetLabelSize(0.06);
//        py2Plots[j]->GetYaxis()->SetTitleSize(0.08);
//        py2Plots[j]->GetYaxis()->SetTitleOffset(0.5);
//        py2Plots[j]->GetYaxis()->CenterTitle();
//        py2Plots[j]->Write();
//    }
//
//    TCanvas *canvas4 = new TCanvas("canvas4", "canvas4", 2560, 1440);
//    canvas4->Divide(3, 3, 0.000000001, 0.01);
//
//    for(Int_t i = 0; i < 8; i++){
//        canvas4->cd(i+1);
//        pyPlots[i]->Draw("AEP");
//    }
//
//    canvas4->Write();
//
//    TCanvas *canvas5 = new TCanvas("canvas5", "canvas5", 2560, 1440);
//    canvas5->Divide(2, 2, 0.000000001, 0.01);
//
//    for(Int_t i = 0; i < 3; i++){
//        canvas5->cd(i+1);
//        py2Plots[i]->Draw("AEP");
//    }
//
//    canvas5->Write();
//
//    double GlisterCx250x[6] = {cos(20*TMath::DegToRad()), cos(30*TMath::DegToRad()), cos(50*TMath::DegToRad()), cos(70*TMath::DegToRad()), cos(90*TMath::DegToRad()), cos(100*TMath::DegToRad())};
//    double GlisterCx250y[6] = {0., 0., -0.1, -0.2,  -0.2, -0.3};
//    double GlisterCx250yErr[6] = {0.01, 0.01, 0.01, 0.01, 0.02, 0.15};
//    double GlisterCx350x[8] = {cos(20*TMath::DegToRad()), cos(30*TMath::DegToRad()), cos(40*TMath::DegToRad()), cos(50*TMath::DegToRad()), cos(70*TMath::DegToRad()), cos(90*TMath::DegToRad()), cos(110*TMath::DegToRad()), cos(120*TMath::DegToRad())};
//    double GlisterCx350y[8] = {0.05, 0., -0.02, -0.05, -0.2, -0.2, 0.15, -0.05};
//    double GlisterCx350yErr[8] = {0.01, 0.01, 0.01, 0.01, 0.01, 0.02, 0.1, 0.1};
//
//    GlisterCx250 = new TGraphErrors(6, GlisterCx250x, GlisterCx250y, 0, GlisterCx250yErr);
//    GlisterCx250->SetMarkerColor(618);
//    GlisterCx250->SetLineColor(617);
//    GlisterCx250->SetMarkerStyle(33);
//    GlisterCx250->SetMarkerSize(2);
//    GlisterCx250->SetName("GlisterCx250");
//    GlisterCx250->Write();
//
//    GlisterCx350 = new TGraphErrors(8, GlisterCx350x, GlisterCx350y, 0, GlisterCx350yErr);
//    GlisterCx350->SetMarkerColor(618);
//    GlisterCx350->SetLineColor(618);
//    GlisterCx350->SetMarkerStyle(33);
//    GlisterCx350->SetMarkerSize(2);
//    GlisterCx350->SetName("GlisterCx350");
//    GlisterCx350->Write();
//
//    TCanvas *canvas6 = new TCanvas("canvas6", "canvas6", 2560, 1440);
//    canvas6->Divide(2, 1, 0.000000001, 0.01);
//
//    for(Int_t i = 0; i < 2; i++){
//        canvas6->cd(i+1);
//        CxPlots[0][i]->Draw("AEP");
//        CxPlots[1][i]->Draw("SAMEEP");
//        CxPlots[2][i]->Draw("SAMEEP");
//        if(i == 0) GlisterCx250->Draw("SAMEEP");
//        else if (i == 1) GlisterCx350->Draw("SAMEEP");
//    }
//
//    canvas6->Write();
//
////    Couting loop for results to be inputted into LaTeX if needed
////    for(Int_t k = 0; k < 3; k++){
////        for(Int_t i = 0; i < 8; i++){ // Energy
////            Double_t EValue = 250 + (i*100);
////            cout << std::setprecision(3) << EValue << " #pm " << 50 << "\t";
////            PVal = (Pn90CM->Eval(EValue, 0));
////            for (Int_t j = 0; j < 3; j++){
////                    if(k == 0) {cout << 0 << "\t";}
////                    else if (k == 1){cout << PVal*AEff[i][j] << "\t";}
////                    else if (k == 2){cout << py[i][j] << "\t";}
////                    cout << std::setprecision(3) << std::fixed << Cx[k][i][j]<< " #pm " << CxErr[k][i][j] << "\t";
////            }
////            cout << endl;
////        }
////    }
//
    f3.Write();

}
