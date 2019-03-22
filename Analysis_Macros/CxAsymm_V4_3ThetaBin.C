#include "./includes.h"

// Version 4 utilises the 2D histograms of AEff as a fn of PhiSc for each bin
// AEff is extracted from this and averaged to calculate the value of Cx'

// define a function with 4 parameters
Double_t fitf(Double_t *x,Double_t *par)
{
    Double_t fitval = 0;
    fitval =  ((par[0]*sin(x[0]))/(1 + (par[1]*cos(x[0]))));
    return fitval;
}

void CxAsymm_V4_3ThetaBin() {

    gStyle->SetTitleFontSize(0.06);
    gStyle->SetOptStat(0);

    char PosHelHistName[60];
    char NegHelHistName[60];
    char PosHelProfName[60];
    char NegHelProfName[60];
    char PosHelProjName[60];
    char NegHelProjName[60];
    char AsymmHistName[60];
    char NeutronEThetaScName[60];;
    char AsymmHistTitle[60];

    TProfile* PhiScAEffXProfile[9][2][8][3]; // Need each helicity so 0 is -ve, 1 is +ve for first index
    TH1D* PhiScAEffXProjection[9][2][8][3];
    TH1D* AsymmHists[9][8][3];
    TGraphErrors* CxEg[10][3]; // As fn of Angular bin for fixed energies
    double SinAmps[9][8][3]; // To store and get Cx amps
    double SinAmpErrs[9][8][3];// To store and get Cx amps
    Double_t Cx[9][8][3];
    Double_t CxErr[9][8][3];

    Double_t NumerSum1;
    Double_t NumerSum2;
    Double_t DenomSum1;
    Double_t DenomSum2;
    Double_t AEff[9][8][3];
    Double_t AEffErr[9][8][3];
    Double_t CorrFac[9][8][3];
    Double_t AEffCorr[9][8][3];

    Double_t Pm = 0.7666; // Max beam polarisation and energy
    Double_t Em = 1557;
    Double_t CircPol[9][8];
    Double_t CircPolErr [9][8];

    TF1 *AsymmFunc = new TF1("AsymmFit",  fitf, -3.0, 3.0, 2); //Give a name and range to the fitting funcion
    AsymmFunc->SetParNames("SinAmp", "CosAmp"); //Name the parameters
    AsymmFunc->SetParameter(0, 0);
    TF1 *Pn90CM = new TF1("Pn90CM", "1.64576-2.95484*(x/1000)+0.684577*(x/1000)**2-0.65*90**2/4/((x-560)**2+90**2/4)+(5.32305-35.3819*(x/1000)+70.145*(x/1000)**2-44.2899*(x/1000)**3)",200,1000);
    TFile *f = new TFile("/home/sjdkay/work/A2/a2GoAT/Physics_Total_153_53/Physics_Total_Amo153_Lin53_Combined_EgShift_V2.root"); // Open the latest PTotal file to load histograms from

    // PROFILE of PhiScAEff gives the averaged AEff values for each PhiSc bin, e.g. will give a 1D hist of this
    // PhiScAEff440PosHelCM2->ProfileX("Test2", 0, -1, "ode")
    // PROJECTION of PhiScAEff gives the normal 1D hist of PhiSc again, so this is what we need to take as the asymmetry
    // PhiScAEff440PosHelCM2->ProjectionX("Test", 0, -1, "ode")
    // Options o keep axes, d, draws (so don't have this) and e gets the error
    // See https://root.cern.ch/doc/master/classTProfile.html#a1ff9340284c73ce8762ab6e7dc0e6725 for error output options
    // Default version probably fine averages over all bins y
    // So for each Asymmetry should probably determine a WEIGHTED average value for AEff

    TFile Results_File("AsymmFits_PTotal_153_53_V4_TEST.root", "RECREATE");
    Double_t EVals[9] = {250, 260, 270, 280, 290, 310, 320, 330, 340};

    TDirectory *Asymm_Hists = Results_File->mkdir("Asymmetry Histograms"); Asymm_Hists->cd();

    for(Int_t i = 0; i < 9; i++){ // Centre of energy bin
        for(Int_t j = 0; j < 8; j++){ // Energy
            CircPol[i][j] = (Pm*(EVals[i]+(j*100))*(Em+((1/3)*(Em-(EVals[i]+(j*100))))))/((TMath::Power(Em,2))+(TMath::Power((Em-(EVals[i]+(j*100))),2))-((2/3)*Em*(Em-(EVals[i]+(j*100)))));
            CircPolErr[i][j] = 0.; // For now, average between edges or something?
            for(Int_t k = 0; k < 3; k++)){ // CosThetaBin
            sprintf(PosHelHistName, "PhiScAEff%iPosHelCM%i", EVals[i]+(j*100), k+1);
            sprintf(NegHelHistName, "PhiScAEff%iNegHelCM%i", EVals[i]+(j*100), k+1);
            sprintf(PosHelProfName, "AEff%iPosHelCM%i", EVals[i]+(j*100), k+1);
            sprintf(NegHelProfName, "AEff%iNegHelCM%i", EVals[i]+(j*100), k+1);
            sprintf(PosHelProjName, "PhiScProj%iPosHelCM%i", EVals[i]+(j*100), k+1);
            sprintf(NegHelProjName, "PhiScProj%iNegHelCM%i", EVals[i]+(j*100), k+1);
            PhiScAEffXProfile[i][0][j][k] = ((TH2F*)f->Get(NegHelHistName))->ProfileX(NegHelProfName, 0, -1, "oe");
            PhiScAEffXProfile[i][1][j][k] = ((TH2F*)f->Get(PosHelHistName))->ProfileX(PosHelProfName, 0, -1, "oe");
            PhiScAEffXProjection[i][0][j][k] = (((TH2F*)f->Get(NegHelHistName))->ProjectionX(NegHelProjName, 0, -1, "oe"));
            PhiScAEffXProjection[i][1][j][k] = (((TH2F*)f->Get(PosHelHistName))->ProjectionX(PosHelProjName, 0, -1, "oe"));

            Int_t NBinsProf = PhiScAEffXProfile[i][0][j][k]->GetNbinsX();
            DenomSum1 = 0;
            DenomSum2 = 0;
            NumerSum1 = 0;
            NumerSum2 = 0;
            for(Int_t l = 1; l < NBinsProf; l++){ // Calculate weighted average of AEff for each Energy/CosTheta bin for each helicity
                DenomSum1 =+ (1)/(Power(PhiScAEffXProfile[i][0][j][k]->GetBinError(l), 2));
                DenomSum2 =+ (1)/(Power(PhiScAEffXProfile[i][1][j][k]->GetBinError(l), 2));
                NumerSum1 =+ (PhiScAEffXProfile[i][0][j][k]->GetBinContent(l))/(Power(PhiScAEffXProfile[i][0][j][k]->GetBinError(l), 2));
                NumerSum2 =+ (PhiScAEffXProfile[i][1][j][k]->GetBinContent(l))/(Power(PhiScAEffXProfile[i][1][j][k]->GetBinError(l), 2));
            }

            AEff[i][j][k] = (NumerSum1+NumerSum2)/(DenomSum1+DenomSum2);
            AEffErr[i][j][k] = sqrt(1/(DenomSum1+DenomSum2));

            // Will need to modify this, need to plot PhiScAeffCorr in main analysis routine
            //CorrFac[i][j][k] = (1+exp(1.81572-(0.0139530*MeanX[j][k])));
            //AEffCorr[i][j][k] = AEff[i][j][k] * CorrFac[i][j][k]; //Analysing power for 12C based on correction factor from Mikhail

            sprintf(AsymmHistName, "CxAsymm%iCM%i", EVals[i]+(j*100), k+1);
            sprintf(AsymmHistTitle, "#phi_{Sc} Asymmetry (%i #pm 50) MeV CM%i",EVals[i]+(j*100), k+1);

            AsymmHists[i][j][k] = (TH1D*) ((PhiScAEffXProjection[i][1][j][k]->GetAsymmetry(PhiScAEffXProjection[i][0][j][k])));
            AsymmHists[i][j][k]->SetName(AsymmHistName);
            AsymmHists[i][j][k]->SetTitle(AsymmHistTitle);
            AsymmHists[i][j][k]->GetXaxis()->SetTitle("#phi_{Sc}/Rad");
            AsymmHists[i][j][k]->GetYaxis()->SetTitle("Asymmetry");
            AsymmHists[i][j][k]->GetXaxis()->SetLabelSize(0.06);
            AsymmHists[i][j][k]->GetXaxis()->SetTitleSize(0.06);
            AsymmHists[i][j][k]->GetXaxis()->SetTitleOffset(0.7);
            AsymmHists[i][j][k]->GetXaxis()->CenterTitle();
            AsymmHists[i][j][k]->GetYaxis()->SetLabelSize(0.06);
            AsymmHists[i][j][k]->GetYaxis()->SetTitleSize(0.06);
            AsymmHists[i][j][k]->GetYaxis()->SetTitleOffset(0.7);
            AsymmHists[i][j][k]->GetYaxis()->CenterTitle();

            AsymmFit->SetLineColor(4);
            AsymmFit->FixParameter(1, 0.);
            AsymmHists[i][j][k]->Fit("AsymmFit", "Q");
            AsymmHists[i][j][k]->SetLineColor(4);
            SinAmps[i][j][k] = AsymmFit->GetParameter(0);
            SinAmpErrs[i][j][k] = AsymmFit->GetParError(0);
            Cx[i][j][k] = SinAmps[i][j][k]/(AEff[i][j][k]*CircPol[i][j]);
            CxErr[i][j][k] = sqrt( ((SinAmpErrs[i][j][k]/(CircPol[i][j]*AEff[i][j][k]))**2) + (((AEffErr[i][j][k]*SinAmps[i][j][k])/(((AEff[i][j][k])**2)*CircPol[i][j]))**2) + (((CircPolErr[i][j]*SinAmps[i][j][k])/(((CircPol[i][j])**2)*AEff[i][j][k]))**2));
            AsymmHists[i][j][k]->Write();
            }
        }
    }

    TDirectory *Cx_Plots_Individual = Results_File.mkdir("Cx Plots Individual"); Cx_Plots_Individual->cd();

    for(Int_t i = 0; i < 9; i++){ // Centre of energy bin
        for(Int_t k = 0; k < 3; k++)){ // CosThetaBin
            CxEg[i][k] = new TGraphErrors(8);
            for(Int_t j = 0; j < 8; j++){
                CxEg[i][k]->SetPoint(j, EVals[i]+(j*100), Cx[i][j][k]);
                CxEg[i][k]->SetPointError(j, 50, CxErr[i][j][k]);
            }
            CxEg[i][k]->SetName(Form("Cx_CM%i_V%i", k+1, i+1));
            CxEg[i][k]->SetTitle(Form("C_{x'}(E_{#gamma}) #theta_{CM} %i - %i^{#circ}", k*60 ,60 + (k*60)));
            CxEg[i][k]->SetMarkerSize(1.2);
            CxEg[i][k]->SetMarkerColor(4);
            CxEg[i][k]->SetMarkerStyle(20);
            CxEg[i][k]->SetMarkerSize(1.5);
            CxEg[i][k]->SetLineColor(4);
            CxEg[i][k]->GetXaxis()->SetTitle("E_{#gamma}/MeV");
            CxEg[i][k]->GetXaxis()->SetRangeUser(100, 1200);
            CxEg[i][k]->GetYaxis()->SetRangeUser(-2, 2);
            CxEg[i][k]->GetYaxis()->SetTitle("C_{x'}");
            CxEg[i][k]->GetXaxis()->SetLabelSize(0.04);
            CxEg[i][k]->GetXaxis()->SetTitleSize(0.04);
            CxEg[i][k]->GetXaxis()->SetTitleOffset(0.9);
            CxEg[i][k]->GetXaxis()->CenterTitle();
            CxEg[i][k]->GetYaxis()->SetLabelSize(0.04);
            CxEg[i][k]->GetYaxis()->SetTitleSize(0.04);
            CxEg[i][k]->GetYaxis()->SetTitleOffset(0.5);
            CxEg[i][k]->GetYaxis()->CenterTitle();
            CxEg[i][k]->Write();
        }
    }

    TDirectory *Cx_Plots = Results_File.mkdir("Cx Plots"); Cx_Plots->cd();

    for(Int_t k = 0; k < 3; k++)){ // CosThetaBin
        CxEg[9][k] = new TGraphErrors(72);
            Int_t BinNo = 0;
            for(Int_t i = 0; i < 9; i++){
                for(Int_t j = 0; j < 8; j++){
                CxEg[9][k]->SetPoint(BinNo, EVals[i]+(j*100), Cx[i][j][k]);
                CxEg[9][k]->SetPointError(BinNo, 50, CxErr[i][j][k]);
                BinNo += 1;
            }
        }
        CxEg[9][k]->SetName(Form("Cx_CM%i", k+1));
        CxEg[9][k]->SetTitle(Form("C_{x'}(E_{#gamma}) #theta_{CM} %i - %i^{#circ}", k*60 ,60 + (k*60)));
        CxEg[9][k]->SetMarkerSize(1.2);
        CxEg[9][k]->SetMarkerColor(4);
        CxEg[9][k]->SetMarkerStyle(20);
        CxEg[9][k]->SetMarkerSize(1.5);
        CxEg[9][k]->SetLineColor(4);
        CxEg[9][k]->GetXaxis()->SetTitle("E_{#gamma}/MeV");
        CxEg[9][k]->GetXaxis()->SetRangeUser(100, 1200);
        CxEg[9][k]->GetYaxis()->SetRangeUser(-2, 2);
        CxEg[9][k]->GetYaxis()->SetTitle("C_{x'}");
        CxEg[9][k]->GetXaxis()->SetLabelSize(0.04);
        CxEg[9][k]->GetXaxis()->SetTitleSize(0.04);
        CxEg[9][k]->GetXaxis()->SetTitleOffset(0.9);
        CxEg[9][k]->GetXaxis()->CenterTitle();
        CxEg[9][k]->GetYaxis()->SetLabelSize(0.04);
        CxEg[9][k]->GetYaxis()->SetTitleSize(0.04);
        CxEg[9][k]->GetYaxis()->SetTitleOffset(0.5);
        CxEg[9][k]->GetYaxis()->CenterTitle();
        CxEg[9][k]->Write();
    }

    TDirectory *Results_250MeV = Results_File.mkdir("Results_250MeV"); Results_250MeV->cd(); // Results for normal binning
    // Plot of Cx(Eg) for each CM range with fixed energy binning and then canvas with all 5 on one
    for(Int_t k = 0; k < 3; k++)){ // CosThetaBin
        TCanvas *C1 = new TCanvas(Form("CxEg_50MeVCent_CM%i", k+1), Form("C_{x'}(E_{#gamma}) #theta_{CM} %i - %i^{#circ}", k*36 ,36 + (k*36)), 2560, 1440);
        CxEg[0][k]->Draw("AEP");
        C1->Write();
    }
    TCanvas *C2= new TCanvas("CxEg_50MeVCent", "CxEg_50MeVCent" , 2560, 1440);
    C2->Divide(3, 2, 0.000000001, 0.01);
    for(Int_t k = 0; k < 3; k++){
        C2->cd(k+1);
        CxEg[0][k]->Draw("AEP");
    }
    C2->Write();

    TDirectory *Results_MovingBin = Results_File.mkdir("Results_MovingBin"); Results_MovingBin->cd(); // Results for normal binning
    // Plot of Cx(Eg) for each CM range with moving energy binning and then canvas with all 5 on one
    for(Int_t k = 0; k < 3; k++)){ // CosThetaBin
        TCanvas *C3 = new TCanvas(Form("CxEg_CM%i", k+1), Form("CxEg_CM%i", k+1), 2560, 1440);
        CxEg[9][k]->Draw("AEP");
        C3->Write();
    }
    TCanvas *C4 = new TCanvas("CxEg", "CxEg" , 2560, 1440);
    C4->Divide(2, 2, 0.000000001, 0.01);
    for(Int_t k = 0; k < 3; k++){
        C4->cd(k+1);
        CxEg[9][k]->Draw("AEP");
    }
    C4->Write();

    Results_File.Write();
    Results_File.Close();

}
