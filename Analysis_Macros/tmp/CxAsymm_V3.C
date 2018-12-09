#include "./includes.h"

// define a function with 4 parameters
Double_t fitf(Double_t *x,Double_t *par)
{
    Double_t fitval = 0;
    fitval =  ((par[0]*sin(x[0]))/(1 + (par[1]*cos(x[0]))));
    return fitval;
}

void CxAsymm_V3() {

    gStyle->SetTitleFontSize(0.06);

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

    TH1F* AsymmHists[3][6][3];
    TGraphErrors* CxPlots[3][6]; // As fn of Angular bin for fixed energies
    TGraphErrors* Cx2Plots[3][3]; // As fn of energy bin for fixed angular bins
    TGraphErrors* pyPlots[6]; // As fn of Angular bin for fixed energies
    TGraphErrors* py2Plots[3]; // As fn of Angular bin for fixed energies

    double SinAmps[3][6][3]; // To store and get Cx amps
    double SinAmpErrs[3][6][3];// To store and get Cx amps
    double CosAmps[6][3];
    double CosAmpErrs[6][3];
    double MeanX[6][3];
    double MeanY[6][3];

    Double_t Cx[3][6][3]; // Fixed P Chi2 fit, Fixed P LL fit, Variable p Chi2 fit, Variable P LL fit
    Double_t CxErr[3][6][3];
    Double_t py[6][3];
    Double_t pyErr[6][3];

    Double_t Cx2[3][6]; // Fixed P Chi2 fit, Fixed P LL fit, Variable p Chi2 fit, Variable P LL fit
    Double_t Cx2Err[3][6];

    Double_t py2[3][6];
    Double_t py2Err[3][6];

    Double_t PVal;
    Double_t ECentre;
    Double_t AmpVal;
    Double_t AEff[6][3];
    Double_t CorrFac[6][3];
    Double_t AEffCorr[6][3];

    TF1 *AsymmFunc = new TF1("AsymmFit",  fitf, -3.0, 3.0, 2); //Give a name and range to the fitting funcion
    AsymmFunc->SetParNames("SinAmp", "CosAmp"); //Name the parameters
    AsymmFunc->SetParameter(0, 0);
    TFile *f = new TFile("/scratch/Mainz_Software/a2GoAT/Physics_Total_Amo142_Lin42_Combined.root"); // Open the latest PTotal file to load histograms from
    TFile *fAy = new TFile ("/scratch/Mainz_Software/a2GoAT/npAy.root");
    TF1 *Pn90CM = new TF1("Pn90CM", "1.64576-2.95484*(x/1000)+0.684577*(x/1000)**2-0.65*90**2/4/((x-560)**2+90**2/4)+(5.32305-35.3819*(x/1000)+70.145*(x/1000)**2-44.2899*(x/1000)**3)",300,700);

    for(Int_t i = 0; i < 3; i++){ // Fit version
        for(Int_t j = 0; j < 6; j++){ // Energy
            ECentre = 250+(j*100); // Get centre of energy bin
            PVal = (Pn90CM->Eval(ECentre, 0)); // Get PValue for energy bin to fix parameter with
            for(Int_t k = 0; k < 3; k++)){ // CosThetaBin

                sprintf(PosHelHistName, "PhiSc%iPosHelCM%i", 250+(j*100), k+1);
                sprintf(NegHelHistName, "PhiSc%iNegHelCM%i", 250+(j*100), k+1);
                sprintf(AsymmHistName, "CxAsymm%iCM%iFit%i", 250+(j*100), k+1, i+1);
                sprintf(NeutronEThetaScName, "NeutronEThetaSc%iCM%i", 200+(j*100), k+1);
                AsymmHists[i][j][k] = (TH1F*) (((TH1F*)f->Get(PosHelHistName))->GetAsymmetry(((TH1F*)f->Get(NegHelHistName)))));
                AsymmHists[i][j][k]->SetName(AsymmHistName);

                MeanX[j][k] = ((TH2D*)f->Get(NeutronEThetaScName))->GetMean(1);
                MeanY[j][k] = ((TH2D*)f->Get(NeutronEThetaScName))->GetMean(2);

                AEff[j][k] = ((TH2F*)fAy->Get("nppnAy"))->Interpolate(MeanX[j][k], MeanY[j][k]);
                CorrFac[j][k] = (1+exp(1.81572-(0.0139530*MeanX[j][k])));
                AEffCorr[j][k] = AEff[j][k] * CorrFac[j][k]; //Analysing power for 12C based on correction factor from Mikhail
                AmpVal = PVal*AEff[j][k];

                if (i == 0){
                    AsymmFit->SetLineColor(4);
                    AsymmFit->FixParameter(1, 0.);
                    AsymmHists[i][j][k]->Fit("AsymmFit", "Q");
                    SinAmps[i][j][k] = AsymmFit->GetParameter(0);
                    SinAmpErrs[i][j][k] = AsymmFit->GetParError(0);
                }

                if(i == 1){
                    AsymmFit->SetLineColor(4);
                    AsymmFit->SetLineStyle(10);
                    AsymmFit->FixParameter(1, AmpVal);
                    AsymmHists[i][j][k]->Fit("AsymmFit", "Q");
                    SinAmps[i][j][k] = AsymmFit->GetParameter(0);
                    SinAmpErrs[i][j][k] = AsymmFit->GetParError(0);
                }

                if(i == 2){
                    AsymmFit->SetLineColor(4);
                    AsymmFit->SetLineStyle(10);
                    AsymmFit->ReleaseParameter(1);
                    AsymmFit->SetParLimits(1, -1*fabs(AEff[j][k]), fabs(AEff[j][k]));
                    AsymmHists[i][j][k]->Fit("AsymmFit", "Q");
                    SinAmps[i][j][k] = AsymmFit->GetParameter(0);
                    SinAmpErrs[i][j][k] = AsymmFit->GetParError(0);
                    CosAmps[j][k] = AsymmFit->GetParameter(1);
                    CosAmpErrs[j][k] = AsymmFit->GetParError(1);
                    cout << AEff[j][k] << "   " << AsymmFit->GetParameter(1)/AEff[j][k] << endl;
                }


            }
        }
    }

    TFile f1("AsymmFits_PTotal_142_42_V4.root", "RECREATE");

    for(Int_t i = 0; i < 3; i++){ // Fit version
        for(Int_t j = 0; j < 6; j++){ // Energy
            for(Int_t k = 0; k < 3; k++)){
                AsymmHists[i][j][k]->Write();
            }
        }
    }

//    TCanvas *canvas20 = new TCanvas("canvas20","canvas20", 1920, 1080);
//    canvas20->Divide(2, 2);
//    for(int i = 1 ; i < 5 ; i++){
//        canvas20->cd(i);
//        AsymmHists[0][i-1]->Draw("EP");
//        for(int j = 1 ; j < 2 ; j++){
//            AsymmHists[j][i-1]->Draw("SAMEEP");
//        }
//    }
//
//    canvas20->Write();

    f1.Write();
    f1.Close();

    TFile *f2= TFile::Open("/scratch/Mainz_Software/a2GoAT/CircPol_Aug16.root");

    for (Int_t m = 0; m < 3; m++){ // Fit Version
        for (Int_t n = 0; n < 6; n++){ //E
            Double_t EPoint = 250 + (n*100);
            for(Int_t l = 0; l < 3; l++)){ // Theta CM bin
                Cx[m][n][l] = SinAmps[m][n][l]/(AEff[n][l]*(Graph->Eval(EPoint ,0)));
                CxErr[m][n][l] = fabs(SinAmpErrs[m][n][l]/(AEff[n][l]*(Graph->Eval(EPoint ,0))));
                if(m == 2){
                    py[n][l] = CosAmps[n][l]/AEff[n][l];
                    pyErr[n][l] = fabs(CosAmpErrs[n][l]/AEff[n][l]);
                }
            }
        }
    }

    TFile f3("Cx_Plots_142_42_V4.root", "RECREATE");

    double x[3] = {2./3., 0., -2./3.};
    double ex[3] = {1./3., 1./3., 1./3. };

    for(Int_t i = 0; i < 3; i++){
        for(Int_t j = 0; j < 6; j++){
            sprintf(name, "Cx_%i_V%i", 250+(j*100), i+1);
            sprintf(title, "Cx(Cos#theta_{CM}) E_{#gamma} %i #pm 50 MeV", 250+(j*100), i+1);
            CxPlots[i][j] = new TGraphErrors(3 , x, Cx[i][j], ex, CxErr[i][j]);
            CxPlots[i][j]->SetName(name);
            CxPlots[i][j]->SetTitle(title);
            CxPlots[i][j]->SetMarkerSize(1.2);
            if(i == 0){
                CxPlots[i][j]->SetMarkerColor(4);
                CxPlots[i][j]->SetMarkerStyle(20);
                CxPlots[i][j]->SetMarkerSize(1.5);
                CxPlots[i][j]->SetLineColor(4);
            }

            else if (i == 1){
                CxPlots[i][j]->SetMarkerColor(2);
                CxPlots[i][j]->SetMarkerStyle(21);
                CxPlots[i][j]->SetLineColor(2);
            }

            else if (i == 2){
                CxPlots[i][j]->SetMarkerColor(1);
                CxPlots[i][j]->SetMarkerStyle(22);
                CxPlots[i][j]->SetLineColor(1);
            }
            CxPlots[i][j]->GetXaxis()->SetTitle("Cos#theta_{CM}");
            CxPlots[i][j]->GetXaxis()->SetRangeUser(-1, 1);
            CxPlots[i][j]->GetYaxis()->SetRangeUser(-2, 2);
            CxPlots[i][j]->GetYaxis()->SetTitle("C_{x}");
            CxPlots[i][j]->GetXaxis()->SetLabelSize(0.06);
            CxPlots[i][j]->GetXaxis()->SetTitleSize(0.06);
            CxPlots[i][j]->GetXaxis()->SetTitleOffset(0.7);
            CxPlots[i][j]->GetXaxis()->CenterTitle();
            CxPlots[i][j]->GetYaxis()->SetLabelSize(0.06);
            CxPlots[i][j]->GetYaxis()->SetTitleSize(0.08);
            CxPlots[i][j]->GetYaxis()->SetTitleOffset(0.5);
            CxPlots[i][j]->GetYaxis()->CenterTitle();
            CxPlots[i][j]->Write();
        }
    }

    TCanvas *canvas = new TCanvas("canvas","canvas", 2560, 1440);
    canvas->Divide(4, 2, 0.000000001, 0.01);

    for(Int_t i = 0; i < 6; i++){
        canvas->cd(i+1);
        CxPlots[0][i]->Draw("AEP");
        CxPlots[1][i]->Draw("SAMEEP");
        CxPlots[2][i]->Draw("SAMEEP");
    }

    canvas->cd(7);

    leg = new TLegend(0.1, 0.1, 0.9, 0.9);
    leg->AddEntry(CxPlots[0][0], "Fixed py = 0, #chi^{2}", "lepz");
    leg->AddEntry(CxPlots[1][0], "Fixed py, #chi^{2}", "lepz");
    leg->AddEntry(CxPlots[2][0], "Variable py, #chi^{2}", "lepz");

    leg->Draw("SAME");

    canvas->Write();

    double x2[6] = {250, 350, 450, 550, 650, 750};
    double ex2[6] = {50, 50, 50, 50, 50, 50};

    for(Int_t i = 0; i < 3; i++){
        for(Int_t k = 0; k < 3; k++){
            for(Int_t j = 0; j < 6; j++){
                Cx2[i][j] = Cx[i][j][k];
                Cx2Err[i][j] = CxErr[i][j][k];
            }
            sprintf(name2, "Cx_CM%i_V%i", k+1, i+1);
            sprintf(title2, "Cx(E_{#gamma}) #theta_{CM} %i - %i", k*60 ,60 + (k*60));
            Cx2Plots[i][k] = new TGraphErrors(6 , x2, Cx2[i], ex2, Cx2Err[i]);
            Cx2Plots[i][k]->SetName(name2);
            Cx2Plots[i][k]->SetTitle(title2);
            Cx2Plots[i][k]->SetMarkerSize(1.2);
            if(i == 0){
                Cx2Plots[i][k]->SetMarkerColor(4);
                Cx2Plots[i][k]->SetMarkerStyle(20);
                Cx2Plots[i][k]->SetMarkerSize(1.5);
                Cx2Plots[i][k]->SetLineColor(4);
            }

            else if (i == 1){
                Cx2Plots[i][k]->SetMarkerColor(2);
                Cx2Plots[i][k]->SetMarkerStyle(21);
                Cx2Plots[i][k]->SetLineColor(2);
            }
            else if (i == 2){
                Cx2Plots[i][k]->SetMarkerColor(1);
                Cx2Plots[i][k]->SetMarkerStyle(22);
                Cx2Plots[i][k]->SetLineColor(1);
            }
            Cx2Plots[i][k]->GetXaxis()->SetTitle("E_{#gamma}");
            Cx2Plots[i][k]->GetXaxis()->SetRangeUser(100, 800);
            Cx2Plots[i][k]->GetYaxis()->SetRangeUser(-2, 2);
            Cx2Plots[i][k]->GetYaxis()->SetTitle("C_{x}");
            Cx2Plots[i][k]->GetXaxis()->SetLabelSize(0.06);
            Cx2Plots[i][k]->GetXaxis()->SetTitleSize(0.06);
            Cx2Plots[i][k]->GetXaxis()->SetTitleOffset(0.7);
            Cx2Plots[i][k]->GetXaxis()->CenterTitle();
            Cx2Plots[i][k]->GetYaxis()->SetLabelSize(0.06);
            Cx2Plots[i][k]->GetYaxis()->SetTitleSize(0.08);
            Cx2Plots[i][k]->GetYaxis()->SetTitleOffset(0.5);
            Cx2Plots[i][k]->GetYaxis()->CenterTitle();
            Cx2Plots[i][k]->Write();
        }
    }

    leg2 = new TLegend(0.0, 0.06, 0.99, 1.0);
    leg2->AddEntry(Cx2Plots[0][0], "Fixed p_{y} = 0, #chi^{2}", "lepz");
    leg2->AddEntry(Cx2Plots[1][0], "Fixed p_{y}, #chi^{2}", "lepz");
    leg2->AddEntry(Cx2Plots[2][0], "Variable p_{y}, #chi^{2}", "lepz");

    TCanvas *canvas2 = new TCanvas("canvas2","canvas2", 2560, 1440);
    canvas2->Divide(2, 2, 0.000000001, 0.01);

    for(Int_t i = 0; i < 3; i++){
        canvas2->cd(i+1);
        Cx2Plots[0][i]->Draw("AEP");
        Cx2Plots[1][i]->Draw("SAMEEP");
        Cx2Plots[2][i]->Draw("SAMEEP");
    }

    canvas2->cd(4);
    leg2->Draw();
    canvas2->Write();

    TCanvas *canvas3 = new TCanvas("canvas3", "canvas3", 2560, 1440);
    canvas3->Divide(2, 2, 0.000000001, 0.01);

    for(Int_t i = 0; i < 3; i++){
        canvas3->cd(i+1);
        Cx2Plots[2][i]->Draw("AEP");
    }

    for(Int_t i = 0; i < 3; i++){
        for(Int_t k = 0; k < 3; k++){
            for(Int_t j = 0; j < 6; j++){
                Cx2[i][j] = Cx[i][j][k];
                Cx2Err[i][j] = CxErr[i][j][k];
            }
            sprintf(name2, "Cx_CM%i_V%i", k+1, i+1);
            sprintf(title2, "Cx(E_{#gamma}) #theta_{CM} %i - %i", k*60 ,60 + (k*60));
            Cx2Plots[i][k] = new TGraphErrors(6 , x2, Cx2[i], ex2, Cx2Err[i]);
            Cx2Plots[i][k]->SetName(name2);
            Cx2Plots[i][k]->SetTitle(title2);
            Cx2Plots[i][k]->SetMarkerSize(1.2);
            if(i == 0){
                Cx2Plots[i][k]->SetMarkerColor(4);
                Cx2Plots[i][k]->SetMarkerStyle(20);
                Cx2Plots[i][k]->SetMarkerSize(1.5);
                Cx2Plots[i][k]->SetLineColor(4);
            }

            else if (i == 1){
                Cx2Plots[i][k]->SetMarkerColor(2);
                Cx2Plots[i][k]->SetMarkerStyle(21);
                Cx2Plots[i][k]->SetLineColor(2);
            }
            else if (i == 2){
                Cx2Plots[i][k]->SetMarkerColor(1);
                Cx2Plots[i][k]->SetMarkerStyle(22);
                Cx2Plots[i][k]->SetLineColor(1);
            }
            Cx2Plots[i][k]->GetXaxis()->SetTitle("E_{#gamma}");
            Cx2Plots[i][k]->GetXaxis()->SetRangeUser(100, 800);
            Cx2Plots[i][k]->GetYaxis()->SetRangeUser(-2, 2);
            Cx2Plots[i][k]->GetYaxis()->SetTitle("C_{x}");
            Cx2Plots[i][k]->GetXaxis()->SetLabelSize(0.06);
            Cx2Plots[i][k]->GetXaxis()->SetTitleSize(0.06);
            Cx2Plots[i][k]->GetXaxis()->SetTitleOffset(0.7);
            Cx2Plots[i][k]->GetXaxis()->CenterTitle();
            Cx2Plots[i][k]->GetYaxis()->SetLabelSize(0.06);
            Cx2Plots[i][k]->GetYaxis()->SetTitleSize(0.08);
            Cx2Plots[i][k]->GetYaxis()->SetTitleOffset(0.5);
            Cx2Plots[i][k]->GetYaxis()->CenterTitle();
            Cx2Plots[i][k]->Write();
        }
    }

    canvas3->Write();

    for (Int_t n = 0; n < 6; n++){ //E
        for(Int_t l = 0; l < 3; l++)){ // Theta CM bin
            py2[l][n] = py[n][l];
            py2Err[l][n] = pyErr[n][l];
        }
    }

    for(Int_t j = 0; j < 6; j++){
        sprintf(name3, "py_%i", 250+(j*100));
        sprintf(title3, "p_{y}(Cos#theta_{CM}) E_{#gamma} %i #pm 50 MeV", 250+(j*100));
        pyPlots[j] = new TGraphErrors(3 , x, py[j], ex, pyErr[j]);
        pyPlots[j]->SetName(name3);
        pyPlots[j]->SetTitle(title3);
        pyPlots[j]->SetMarkerSize(1.2);
        pyPlots[j]->SetMarkerColor(1);
        pyPlots[j]->SetMarkerStyle(22);
        pyPlots[j]->SetMarkerSize(1.2);
        pyPlots[j]->SetLineColor(1);
        pyPlots[j]->GetXaxis()->SetTitle("Cos#theta_{CM}");
        pyPlots[j]->GetXaxis()->SetRangeUser(-1, 1);
        pyPlots[j]->GetYaxis()->SetRangeUser(-2, 2);
        pyPlots[j]->GetYaxis()->SetTitle("p_{y}");
        pyPlots[j]->GetXaxis()->SetLabelSize(0.06);
        pyPlots[j]->GetXaxis()->SetTitleSize(0.06);
        pyPlots[j]->GetXaxis()->SetTitleOffset(0.7);
        pyPlots[j]->GetXaxis()->CenterTitle();
        pyPlots[j]->GetYaxis()->SetLabelSize(0.06);
        pyPlots[j]->GetYaxis()->SetTitleSize(0.08);
        pyPlots[j]->GetYaxis()->SetTitleOffset(0.5);
        pyPlots[j]->GetYaxis()->CenterTitle();
        pyPlots[j]->Write();
    }

    for(Int_t j = 0; j < 3; j++){
        sprintf(name4, "py_CM%i", j+1);
        sprintf(title4, "p_{y}(E_{#gamma}) #theta_{CM} %i - %i", j*60 ,60 + (j*60));
        py2Plots[j] = new TGraphErrors(6 , x2, py2[j], ex2, py2Err[j]);
        py2Plots[j]->SetName(name4);
        py2Plots[j]->SetTitle(title4);
        py2Plots[j]->SetMarkerSize(1.2);
        py2Plots[j]->SetMarkerColor(1);
        py2Plots[j]->SetMarkerStyle(22);
        py2Plots[j]->SetMarkerSize(1.2);
        py2Plots[j]->SetLineColor(1);
        py2Plots[j]->GetXaxis()->SetTitle("E_{#gamma}");
        py2Plots[j]->GetXaxis()->SetRangeUser(200, 800);
        py2Plots[j]->GetYaxis()->SetRangeUser(-2, 2);
        py2Plots[j]->GetYaxis()->SetTitle("p_{y}");
        py2Plots[j]->GetXaxis()->SetLabelSize(0.06);
        py2Plots[j]->GetXaxis()->SetTitleSize(0.06);
        py2Plots[j]->GetXaxis()->SetTitleOffset(0.7);
        py2Plots[j]->GetXaxis()->CenterTitle();
        py2Plots[j]->GetYaxis()->SetLabelSize(0.06);
        py2Plots[j]->GetYaxis()->SetTitleSize(0.08);
        py2Plots[j]->GetYaxis()->SetTitleOffset(0.5);
        py2Plots[j]->GetYaxis()->CenterTitle();
        py2Plots[j]->Write();
    }

    TCanvas *canvas4 = new TCanvas("canvas4", "canvas4", 2560, 1440);
    canvas4->Divide(3, 2, 0.000000001, 0.01);

    for(Int_t i = 0; i < 6; i++){
        canvas4->cd(i+1);
        pyPlots[i]->Draw("AEP");
    }

    canvas4->Write();

    TCanvas *canvas5 = new TCanvas("canvas5", "canvas5", 2560, 1440);
    canvas5->Divide(2, 2, 0.000000001, 0.01);

    for(Int_t i = 0; i < 3; i++){
        canvas5->cd(i+1);
        py2Plots[i]->Draw("AEP");
    }

    canvas5->Write();

    f3.Write();

}
