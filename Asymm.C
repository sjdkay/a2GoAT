#include "./includes.h"

void Asymm(){

     double ParanBin;
     double PerpnBins;
     double ParaBinValues[90];
     double PerpBinValues[90];
     double DiffBinValues[90];
     double SumBinValues[90];
     double AsymmetryValues[90];

     TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_15_11_16/Para/Physics_Total_Para_1_16_1_17.root"); // Open the latest PTotal file to load histograms from
     ParanBins = Phi_Scattered->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < ParanBins; i++){
     	 ParaBinValues[i] = Phi_Scattered->GetBinContent(i+1);
     }

     TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_15_11_16/Perp/Physics_Total_Perp_1_16_1_17.root"); // Open the latest PTotal file to load histograms from
     PerpnBins = Phi_Scattered->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < PerpnBins; i++){
     	 PerpBinValues[i] = Phi_Scattered->GetBinContent(i+1);
     }

     // Define our new histogram to plot the asymmetry here
     TCanvas *canvas = new TCanvas("canvas","canvas", 1920, 1080);
     TH1D *AsymmHist = new TH1D("AsymmHist", "Asymmetry in PhiSc Between Para and Perp", 90, -180, 180);

     for (Int_t i = 0; i < PerpnBins; i++){
     	 DiffBinValues[i] = ParaBinValues[i] - PerpBinValues[i];
	 SumBinValues[i] = ParaBinValues[i] + PerpBinValues[i];
	 AsymmetryValues[i] = DiffBinValues[i]/SumBinValues[i];
	 AsymmHist->SetBinContent(i+1, AsymmetryValues[i]);
     }

     //Save our asymmetry histogram
     AsymmHist->Draw();
     canvas->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_15_11_16/Para_Perp_Asymm_1.png");
     canvas->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_15_11_16/Para_Perp_Asymm_1.pdf");
}
