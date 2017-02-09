#include "./includes.h"

void Asymm(){

     double NegHelnBin;
     double PosHelnBins;
     double NegHelBinValues[90];
     double PosHelBinValues[90];
     double DiffBinValues[90];
     double SumBinValues[90];
     double AsymmetryValues[90];

     TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/Amo/Physics_Total_7_3_2_17.root"); // Open the latest PTotal file to load histograms from
     NegHelnBins = PhiScNegHel->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < NegHelnBins; i++){
     	 NegHelBinValues[i] = PhiScNegHel->GetBinContent(i+1);
     }

     PosHelnBins = PhiScPosHel->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < PosHelnBins; i++){
     	 PosHelBinValues[i] = PhiScPosHel->GetBinContent(i+1);
     }

     // Define our new histogram to plot the asymmetry here
     TCanvas *canvas = new TCanvas("canvas","canvas", 1920, 1080);
     TH1D *AsymmHist = new TH1D("AsymmHist", "Asymmetry in PhiSc Between -ve and +ve Helicity", 90, -180, 180);

     for (Int_t i = 0; i < PosHelnBins; i++){
     	 DiffBinValues[i] = NegHelBinValues[i] - PosHelBinValues[i];
	 SumBinValues[i] = NegHelBinValues[i] + PosHelBinValues[i];
	 AsymmetryValues[i] = DiffBinValues[i]/SumBinValues[i];
	 AsymmHist->SetBinContent(i+1, AsymmetryValues[i]);
     }

     //Save our asymmetry histogram
     AsymmHist->Draw();
     canvas->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_Asymm_7.png");
     canvas->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_Asymm_7.pdf");

     // ###############################################################################################################

     NegHelnBins = Phi_Scattered_275MeV_NegHel->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < NegHelnBins; i++){
     	 NegHelBinValues[i] = Phi_Scattered_275MeV_NegHel->GetBinContent(i+1);
     }

     PosHelnBins = Phi_Scattered_275MeV_PosHel->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < PosHelnBins; i++){
     	 PosHelBinValues[i] = Phi_Scattered_275MeV_PosHel->GetBinContent(i+1);
     }

     // Define our new histogram to plot the asymmetry here
     TCanvas *canvas1 = new TCanvas("canvas1","canvas1", 1920, 1080);
     TH1D *AsymmHist275 = new TH1D("AsymmHist275", "Asymmetry in PhiSc Between -ve and +ve Helicity (275pm25MeV)", 90, -180, 180);

     for (Int_t i = 0; i < PosHelnBins; i++){
     	 DiffBinValues[i] = NegHelBinValues[i] - PosHelBinValues[i];
	 SumBinValues[i] = NegHelBinValues[i] + PosHelBinValues[i];
	 AsymmetryValues[i] = DiffBinValues[i]/SumBinValues[i];
	 AsymmHist275->SetBinContent(i+1, AsymmetryValues[i]);
     }

     //Save our asymmetry histogram
     AsymmHist275->Draw();
     canvas1->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_275MeV_Asymm_7.png");
     canvas1->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_275MeV_Asymm_7.pdf");

     // ###############################################################################################################

     NegHelnBins = Phi_Scattered_325MeV_NegHel->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < NegHelnBins; i++){
     	 NegHelBinValues[i] = Phi_Scattered_325MeV_NegHel->GetBinContent(i+1);
     }

     PosHelnBins = Phi_Scattered_325MeV_PosHel->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < PosHelnBins; i++){
     	 PosHelBinValues[i] = Phi_Scattered_325MeV_PosHel->GetBinContent(i+1);
     }

     // Define our new histogram to plot the asymmetry here
     TCanvas *canvas2 = new TCanvas("canvas2","canvas2", 1920, 1080);
     TH1D *AsymmHist325 = new TH1D("AsymmHist325", "Asymmetry in PhiSc Between -ve and +ve Helicity (325pm25MeV)", 90, -180, 180);

     for (Int_t i = 0; i < PosHelnBins; i++){
     	 DiffBinValues[i] = NegHelBinValues[i] - PosHelBinValues[i];
	 SumBinValues[i] = NegHelBinValues[i] + PosHelBinValues[i];
	 AsymmetryValues[i] = DiffBinValues[i]/SumBinValues[i];
	 AsymmHist325->SetBinContent(i+1, AsymmetryValues[i]);
     }

     //Save our asymmetry histogram
     AsymmHist325->Draw();
     canvas2->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_325MeV_Asymm_7.png");
     canvas2->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_325MeV_Asymm_7.pdf");

     // ###############################################################################################################

     NegHelnBins = Phi_Scattered_375MeV_NegHel->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < NegHelnBins; i++){
     	 NegHelBinValues[i] = Phi_Scattered_375MeV_NegHel->GetBinContent(i+1);
     }

     PosHelnBins = Phi_Scattered_375MeV_PosHel->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < PosHelnBins; i++){
     	 PosHelBinValues[i] = Phi_Scattered_375MeV_PosHel->GetBinContent(i+1);
     }

     // Define our new histogram to plot the asymmetry here
     TCanvas *canvas3 = new TCanvas("canvas3","canvas3", 1920, 1080);
     TH1D *AsymmHist375 = new TH1D("AsymmHist375", "Asymmetry in PhiSc Between -ve and +ve Helicity (375pm25MeV)", 90, -180, 180);

     for (Int_t i = 0; i < PosHelnBins; i++){
     	 DiffBinValues[i] = NegHelBinValues[i] - PosHelBinValues[i];
	 SumBinValues[i] = NegHelBinValues[i] + PosHelBinValues[i];
	 AsymmetryValues[i] = DiffBinValues[i]/SumBinValues[i];
	 AsymmHist375->SetBinContent(i+1, AsymmetryValues[i]);
     }

     //Save our asymmetry histogram
     AsymmHist375->Draw();
     canvas3->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_375MeV_Asymm_7.png");
     canvas3->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_375MeV_Asymm_7.pdf");

     // ###############################################################################################################

     NegHelnBins = Phi_Scattered_425MeV_NegHel->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < NegHelnBins; i++){
     	 NegHelBinValues[i] = Phi_Scattered_425MeV_NegHel->GetBinContent(i+1);
     }

     PosHelnBins = Phi_Scattered_425MeV_PosHel->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < PosHelnBins; i++){
     	 PosHelBinValues[i] = Phi_Scattered_425MeV_PosHel->GetBinContent(i+1);
     }

     // Define our new histogram to plot the asymmetry here
     TCanvas *canvas4 = new TCanvas("canvas4","canvas4", 1920, 1080);
     TH1D *AsymmHist425 = new TH1D("AsymmHist425", "Asymmetry in PhiSc Between -ve and +ve Helicity (425pm25MeV)", 90, -180, 180);

     for (Int_t i = 0; i < PosHelnBins; i++){
     	 DiffBinValues[i] = NegHelBinValues[i] - PosHelBinValues[i];
	 SumBinValues[i] = NegHelBinValues[i] + PosHelBinValues[i];
	 AsymmetryValues[i] = DiffBinValues[i]/SumBinValues[i];
	 AsymmHist425->SetBinContent(i+1, AsymmetryValues[i]);
     }

     //Save our asymmetry histogram
     AsymmHist425->Draw();
     canvas4->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_425MeV_Asymm_7.png");
     canvas4->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_425MeV_Asymm_7.pdf");
    
     // ###############################################################################################################

     NegHelnBins = Phi_Scattered_475MeV_NegHel->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < NegHelnBins; i++){
     	 NegHelBinValues[i] = Phi_Scattered_475MeV_NegHel->GetBinContent(i+1);
     }

     PosHelnBins = Phi_Scattered_475MeV_PosHel->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < PosHelnBins; i++){
     	 PosHelBinValues[i] = Phi_Scattered_475MeV_PosHel->GetBinContent(i+1);
     }

     // Define our new histogram to plot the asymmetry here
     TCanvas *canvas5 = new TCanvas("canvas5","canvas5", 1920, 1080);
     TH1D *AsymmHist475 = new TH1D("AsymmHist475", "Asymmetry in PhiSc Between -ve and +ve Helicity (475pm25MeV)", 90, -180, 180);

     for (Int_t i = 0; i < PosHelnBins; i++){
     	 DiffBinValues[i] = NegHelBinValues[i] - PosHelBinValues[i];
	 SumBinValues[i] = NegHelBinValues[i] + PosHelBinValues[i];
	 AsymmetryValues[i] = DiffBinValues[i]/SumBinValues[i];
	 AsymmHist475->SetBinContent(i+1, AsymmetryValues[i]);
     }

     //Save our asymmetry histogram
     AsymmHist475->Draw();
     canvas5->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_475MeV_Asymm_7.png");
     canvas5->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_475MeV_Asymm_7.pdf");

     // ###############################################################################################################

     NegHelnBins = Phi_Scattered_525MeV_NegHel->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < NegHelnBins; i++){
     	 NegHelBinValues[i] = Phi_Scattered_525MeV_NegHel->GetBinContent(i+1);
     }

     PosHelnBins = Phi_Scattered_525MeV_PosHel->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < PosHelnBins; i++){
     	 PosHelBinValues[i] = Phi_Scattered_525MeV_PosHel->GetBinContent(i+1);
     }

     // Define our new histogram to plot the asymmetry here
     TCanvas *canvas6 = new TCanvas("canvas6","canvas6", 1920, 1080);
     TH1D *AsymmHist525 = new TH1D("AsymmHist525", "Asymmetry in PhiSc Between -ve and +ve Helicity (525pm25MeV)", 90, -180, 180);

     for (Int_t i = 0; i < PosHelnBins; i++){
     	 DiffBinValues[i] = NegHelBinValues[i] - PosHelBinValues[i];
	 SumBinValues[i] = NegHelBinValues[i] + PosHelBinValues[i];
	 AsymmetryValues[i] = DiffBinValues[i]/SumBinValues[i];
	 AsymmHist525->SetBinContent(i+1, AsymmetryValues[i]);
     }

     //Save our asymmetry histogram
     AsymmHist525->Draw();
     canvas6->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_525MeV_Asymm_7.png");
     canvas6->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_525MeV_Asymm_7.pdf");

     // ###############################################################################################################

     NegHelnBins = Phi_Scattered_575MeV_NegHel->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < NegHelnBins; i++){
     	 NegHelBinValues[i] = Phi_Scattered_575MeV_NegHel->GetBinContent(i+1);
     }

     PosHelnBins = Phi_Scattered_575MeV_PosHel->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < PosHelnBins; i++){
     	 PosHelBinValues[i] = Phi_Scattered_575MeV_PosHel->GetBinContent(i+1);
     }

     // Define our new histogram to plot the asymmetry here
     TCanvas *canvas7 = new TCanvas("canvas7","canvas7", 1920, 1080);
     TH1D *AsymmHist575 = new TH1D("AsymmHist575", "Asymmetry in PhiSc Between -ve and +ve Helicity (575pm25MeV)", 90, -180, 180);

     for (Int_t i = 0; i < PosHelnBins; i++){
     	 DiffBinValues[i] = NegHelBinValues[i] - PosHelBinValues[i];
	 SumBinValues[i] = NegHelBinValues[i] + PosHelBinValues[i];
	 AsymmetryValues[i] = DiffBinValues[i]/SumBinValues[i];
	 AsymmHist575->SetBinContent(i+1, AsymmetryValues[i]);
     }

     //Save our asymmetry histogram
     AsymmHist575->Draw();
     canvas7->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_575MeV_Asymm_7.png");
     canvas7->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_575MeV_Asymm_7.pdf");

     // ###############################################################################################################

     NegHelnBins = Phi_Scattered_625MeV_NegHel->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < NegHelnBins; i++){
     	 NegHelBinValues[i] = Phi_Scattered_625MeV_NegHel->GetBinContent(i+1);
     }

     PosHelnBins = Phi_Scattered_625MeV_PosHel->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < PosHelnBins; i++){
     	 PosHelBinValues[i] = Phi_Scattered_625MeV_PosHel->GetBinContent(i+1);
     }

     // Define our new histogram to plot the asymmetry here
     TCanvas *canvas8 = new TCanvas("canvas8","canvas8", 1920, 1080);
     TH1D *AsymmHist625 = new TH1D("AsymmHist625", "Asymmetry in PhiSc Between -ve and +ve Helicity (625pm25MeV)", 90, -180, 180);

     for (Int_t i = 0; i < PosHelnBins; i++){
     	 DiffBinValues[i] = NegHelBinValues[i] - PosHelBinValues[i];
	 SumBinValues[i] = NegHelBinValues[i] + PosHelBinValues[i];
	 AsymmetryValues[i] = DiffBinValues[i]/SumBinValues[i];
	 AsymmHist625->SetBinContent(i+1, AsymmetryValues[i]);
     }

     //Save our asymmetry histogram
     AsymmHist625->Draw();
     canvas8->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_625MeV_Asymm_7.png");
     canvas8->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_625MeV_Asymm_7.pdf");

     // ###############################################################################################################

     NegHelnBins = Phi_Scattered_675MeV_NegHel->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < NegHelnBins; i++){
     	 NegHelBinValues[i] = Phi_Scattered_675MeV_NegHel->GetBinContent(i+1);
     }

     PosHelnBins = Phi_Scattered_675MeV_PosHel->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < PosHelnBins; i++){
     	 PosHelBinValues[i] = Phi_Scattered_675MeV_PosHel->GetBinContent(i+1);
     }

     // Define our new histogram to plot the asymmetry here
     TCanvas *canvas9 = new TCanvas("canvas9","canvas9", 1920, 1080);
     TH1D *AsymmHist675 = new TH1D("AsymmHist675", "Asymmetry in PhiSc Between -ve and +ve Helicity (675pm25MeV)", 90, -180, 180);

     for (Int_t i = 0; i < PosHelnBins; i++){
     	 DiffBinValues[i] = NegHelBinValues[i] - PosHelBinValues[i];
	 SumBinValues[i] = NegHelBinValues[i] + PosHelBinValues[i];
	 AsymmetryValues[i] = DiffBinValues[i]/SumBinValues[i];
	 AsymmHist675->SetBinContent(i+1, AsymmetryValues[i]);
     }

     //Save our asymmetry histogram
     AsymmHist675->Draw();
     canvas9->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_675MeV_Asymm_7.png");
     canvas9->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_675MeV_Asymm_7.pdf");

     // ###############################################################################################################

     NegHelnBins = Phi_Scattered_725MeV_NegHel->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < NegHelnBins; i++){
     	 NegHelBinValues[i] = Phi_Scattered_725MeV_NegHel->GetBinContent(i+1);
     }

     PosHelnBins = Phi_Scattered_725MeV_PosHel->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < PosHelnBins; i++){
     	 PosHelBinValues[i] = Phi_Scattered_725MeV_PosHel->GetBinContent(i+1);
     }

     // Define our new histogram to plot the asymmetry here
     TCanvas *canvas10 = new TCanvas("canvas10","canvas10", 1920, 1080);
     TH1D *AsymmHist725 = new TH1D("AsymmHist725", "Asymmetry in PhiSc Between -ve and +ve Helicity (725pm25MeV)", 90, -180, 180);

     for (Int_t i = 0; i < PosHelnBins; i++){
     	 DiffBinValues[i] = NegHelBinValues[i] - PosHelBinValues[i];
	 SumBinValues[i] = NegHelBinValues[i] + PosHelBinValues[i];
	 AsymmetryValues[i] = DiffBinValues[i]/SumBinValues[i];
	 AsymmHist725->SetBinContent(i+1, AsymmetryValues[i]);
     }

     //Save our asymmetry histogram
     AsymmHist725->Draw();
     canvas10->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_725MeV_Asymm_7.png");
     canvas10->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_725MeV_Asymm_7.pdf");

     // ###############################################################################################################

     NegHelnBins = Phi_Scattered_775MeV_NegHel->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < NegHelnBins; i++){
     	 NegHelBinValues[i] = Phi_Scattered_775MeV_NegHel->GetBinContent(i+1);
     }

     PosHelnBins = Phi_Scattered_775MeV_PosHel->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < PosHelnBins; i++){
     	 PosHelBinValues[i] = Phi_Scattered_775MeV_PosHel->GetBinContent(i+1);
     }

     // Define our new histogram to plot the asymmetry here
     TCanvas *canvas11 = new TCanvas("canvas11","canvas11", 1920, 1080);
     TH1D *AsymmHist775 = new TH1D("AsymmHist775", "Asymmetry in PhiSc Between -ve and +ve Helicity (775pm25MeV)", 90, -180, 180);

     for (Int_t i = 0; i < PosHelnBins; i++){
     	 DiffBinValues[i] = NegHelBinValues[i] - PosHelBinValues[i];
	 SumBinValues[i] = NegHelBinValues[i] + PosHelBinValues[i];
	 AsymmetryValues[i] = DiffBinValues[i]/SumBinValues[i];
	 AsymmHist775->SetBinContent(i+1, AsymmetryValues[i]);
     }

     //Save our asymmetry histogram
     AsymmHist775->Draw();
     canvas11->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_775MeV_Asymm_7.png");
     canvas11->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_775MeV_Asymm_7.pdf");
     
     // ###############################################################################################################

     NegHelnBins = Phi_Scattered_825MeV_NegHel->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < NegHelnBins; i++){
     	 NegHelBinValues[i] = Phi_Scattered_825MeV_NegHel->GetBinContent(i+1);
     }

     PosHelnBins = Phi_Scattered_825MeV_PosHel->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < PosHelnBins; i++){
     	 PosHelBinValues[i] = Phi_Scattered_825MeV_PosHel->GetBinContent(i+1);
     }

     // Define our new histogram to plot the asymmetry here
     TCanvas *canvas12 = new TCanvas("canvas12","canvas12", 1920, 1080);
     TH1D *AsymmHist825 = new TH1D("AsymmHist825", "Asymmetry in PhiSc Between -ve and +ve Helicity (825pm25MeV)", 90, -180, 180);

     for (Int_t i = 0; i < PosHelnBins; i++){
     	 DiffBinValues[i] = NegHelBinValues[i] - PosHelBinValues[i];
	 SumBinValues[i] = NegHelBinValues[i] + PosHelBinValues[i];
	 AsymmetryValues[i] = DiffBinValues[i]/SumBinValues[i];
	 AsymmHist825->SetBinContent(i+1, AsymmetryValues[i]);
     }

     //Save our asymmetry histogram
     AsymmHist825->Draw();
     canvas12->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_825MeV_Asymm_7.png");
     canvas12->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_825MeV_Asymm_7.pdf");

     // ###############################################################################################################

     NegHelnBins = Phi_Scattered_875MeV_NegHel->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < NegHelnBins; i++){
     	 NegHelBinValues[i] = Phi_Scattered_875MeV_NegHel->GetBinContent(i+1);
     }

     PosHelnBins = Phi_Scattered_875MeV_PosHel->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < PosHelnBins; i++){
     	 PosHelBinValues[i] = Phi_Scattered_875MeV_PosHel->GetBinContent(i+1);
     }

     // Define our new histogram to plot the asymmetry here
     TCanvas *canvas13 = new TCanvas("canvas13","canvas13", 1920, 1080);
     TH1D *AsymmHist875 = new TH1D("AsymmHist875", "Asymmetry in PhiSc Between -ve and +ve Helicity (875pm25MeV)", 90, -180, 180);

     for (Int_t i = 0; i < PosHelnBins; i++){
     	 DiffBinValues[i] = NegHelBinValues[i] - PosHelBinValues[i];
	 SumBinValues[i] = NegHelBinValues[i] + PosHelBinValues[i];
	 AsymmetryValues[i] = DiffBinValues[i]/SumBinValues[i];
	 AsymmHist875->SetBinContent(i+1, AsymmetryValues[i]);
     }

     //Save our asymmetry histogram
     AsymmHist875->Draw();
     canvas13->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_875MeV_Asymm_7.png");
     canvas13->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_875MeV_Asymm_7.pdf");

     // #######################################################################################################################
     // #######################################################################################################################
     // 0-10 Degree Plots
     // #######################################################################################################################
     // #######################################################################################################################

     NegHelnBins = Phi_Scattered_275MeV_NegHelTheta010->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < NegHelnBins; i++){
     	 NegHelBinValues[i] = Phi_Scattered_275MeV_NegHelTheta010->GetBinContent(i+1);
     }

     PosHelnBins = Phi_Scattered_275MeV_PosHelTheta010->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < PosHelnBins; i++){
     	 PosHelBinValues[i] = Phi_Scattered_275MeV_PosHelTheta010->GetBinContent(i+1);
     }

     // Define our new histogram to plot the asymmetry here
     TCanvas *canvas14 = new TCanvas("canvas14","canvas14", 1920, 1080);
     TH1D *AsymmHist275Theta010 = new TH1D("AsymmHist275Theta010", "Asymmetry in PhiSc Between -ve and +ve Helicity (275pm25MeV, ThetaSc 0-10)", 90, -180, 180);

     for (Int_t i = 0; i < PosHelnBins; i++){
     	 DiffBinValues[i] = NegHelBinValues[i] - PosHelBinValues[i];
	 SumBinValues[i] = NegHelBinValues[i] + PosHelBinValues[i];
	 AsymmetryValues[i] = DiffBinValues[i]/SumBinValues[i];
	 AsymmHist275Theta010->SetBinContent(i+1, AsymmetryValues[i]);
     }

     //Save our asymmetry histogram
     AsymmHist275Theta010->Draw();
     canvas14->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_275MeV_Theta010_Asymm_7.png");
     canvas14->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_275MeV_Theta010_Asymm_7.pdf");

     // ###############################################################################################################

     NegHelnBins = Phi_Scattered_325MeV_NegHelTheta010->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < NegHelnBins; i++){
     	 NegHelBinValues[i] = Phi_Scattered_325MeV_NegHelTheta010->GetBinContent(i+1);
     }

     PosHelnBins = Phi_Scattered_325MeV_PosHelTheta010->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < PosHelnBins; i++){
     	 PosHelBinValues[i] = Phi_Scattered_325MeV_PosHelTheta010->GetBinContent(i+1);
     }

     // Define our new histogram to plot the asymmetry here
     TCanvas *canvas15 = new TCanvas("canvas15","canvas15", 1920, 1080);
     TH1D *AsymmHist325Theta010 = new TH1D("AsymmHist325Theta010", "Asymmetry in PhiSc Between -ve and +ve Helicity (325pm25MeV, ThetaSc 0-10)", 90, -180, 180);

     for (Int_t i = 0; i < PosHelnBins; i++){
     	 DiffBinValues[i] = NegHelBinValues[i] - PosHelBinValues[i];
	 SumBinValues[i] = NegHelBinValues[i] + PosHelBinValues[i];
	 AsymmetryValues[i] = DiffBinValues[i]/SumBinValues[i];
	 AsymmHist325Theta010->SetBinContent(i+1, AsymmetryValues[i]);
     }

     //Save our asymmetry histogram
     AsymmHist325Theta010->Draw();
     canvas15->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_325MeV_Theta010_Asymm_7.png");
     canvas15->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_325MeV_Theta010_Asymm_7.pdf");

     // ###############################################################################################################

     NegHelnBins = Phi_Scattered_375MeV_NegHelTheta010->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < NegHelnBins; i++){
     	 NegHelBinValues[i] = Phi_Scattered_375MeV_NegHelTheta010->GetBinContent(i+1);
     }

     PosHelnBins = Phi_Scattered_375MeV_PosHelTheta010->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < PosHelnBins; i++){
     	 PosHelBinValues[i] = Phi_Scattered_375MeV_PosHelTheta010->GetBinContent(i+1);
     }

     // Define our new histogram to plot the asymmetry here
     TCanvas *canvas16 = new TCanvas("canvas16","canvas16", 1920, 1080);
     TH1D *AsymmHist375Theta010 = new TH1D("AsymmHist375Theta010", "Asymmetry in PhiSc Between -ve and +ve Helicity (375pm25MeV, ThetaSc 0-10)", 90, -180, 180);

     for (Int_t i = 0; i < PosHelnBins; i++){
     	 DiffBinValues[i] = NegHelBinValues[i] - PosHelBinValues[i];
	 SumBinValues[i] = NegHelBinValues[i] + PosHelBinValues[i];
	 AsymmetryValues[i] = DiffBinValues[i]/SumBinValues[i];
	 AsymmHist375Theta010->SetBinContent(i+1, AsymmetryValues[i]);
     }

     //Save our asymmetry histogram
     AsymmHist375Theta010->Draw();
     canvas16->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_375MeV_Theta010_Asymm_7.png");
     canvas16->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_375MeV_Theta010_Asymm_7.pdf");

     // ###############################################################################################################

     NegHelnBins = Phi_Scattered_425MeV_NegHelTheta010->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < NegHelnBins; i++){
     	 NegHelBinValues[i] = Phi_Scattered_425MeV_NegHelTheta010->GetBinContent(i+1);
     }

     PosHelnBins = Phi_Scattered_425MeV_PosHelTheta010->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < PosHelnBins; i++){
     	 PosHelBinValues[i] = Phi_Scattered_425MeV_PosHelTheta010->GetBinContent(i+1);
     }

     // Define our new histogram to plot the asymmetry here
     TCanvas *canvas17 = new TCanvas("canvas17","canvas17", 1920, 1080);
     TH1D *AsymmHist425Theta010 = new TH1D("AsymmHist425Theta010", "Asymmetry in PhiSc Between -ve and +ve Helicity (425pm25MeV, ThetaSc 0-10)", 90, -180, 180);

     for (Int_t i = 0; i < PosHelnBins; i++){
     	 DiffBinValues[i] = NegHelBinValues[i] - PosHelBinValues[i];
	 SumBinValues[i] = NegHelBinValues[i] + PosHelBinValues[i];
	 AsymmetryValues[i] = DiffBinValues[i]/SumBinValues[i];
	 AsymmHist425Theta010->SetBinContent(i+1, AsymmetryValues[i]);
     }

     //Save our asymmetry histogram
     AsymmHist425Theta010->Draw();
     canvas17->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_425MeV_Theta010_Asymm_7.png");
     canvas17->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_425MeV_Theta010_Asymm_7.pdf");

     // ###############################################################################################################

     NegHelnBins = Phi_Scattered_475MeV_NegHelTheta010->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < NegHelnBins; i++){
     	 NegHelBinValues[i] = Phi_Scattered_475MeV_NegHelTheta010->GetBinContent(i+1);
     }

     PosHelnBins = Phi_Scattered_475MeV_PosHelTheta010->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < PosHelnBins; i++){
     	 PosHelBinValues[i] = Phi_Scattered_475MeV_PosHelTheta010->GetBinContent(i+1);
     }

     // Define our new histogram to plot the asymmetry here
     TCanvas *canvas18 = new TCanvas("canvas18","canvas18", 1920, 1080);
     TH1D *AsymmHist475Theta010 = new TH1D("AsymmHist475Theta010", "Asymmetry in PhiSc Between -ve and +ve Helicity (475pm25MeV, ThetaSc 0-10)", 90, -180, 180);

     for (Int_t i = 0; i < PosHelnBins; i++){
     	 DiffBinValues[i] = NegHelBinValues[i] - PosHelBinValues[i];
	 SumBinValues[i] = NegHelBinValues[i] + PosHelBinValues[i];
	 AsymmetryValues[i] = DiffBinValues[i]/SumBinValues[i];
	 AsymmHist475Theta010->SetBinContent(i+1, AsymmetryValues[i]);
     }

     //Save our asymmetry histogram
     AsymmHist475Theta010->Draw();
     canvas18->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_475MeV_Theta010_Asymm_7.png");
     canvas18->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_475MeV_Theta010_Asymm_7.pdf");

     // ###############################################################################################################

     NegHelnBins = Phi_Scattered_525MeV_NegHelTheta010->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < NegHelnBins; i++){
     	 NegHelBinValues[i] = Phi_Scattered_525MeV_NegHelTheta010->GetBinContent(i+1);
     }

     PosHelnBins = Phi_Scattered_525MeV_PosHelTheta010->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < PosHelnBins; i++){
     	 PosHelBinValues[i] = Phi_Scattered_525MeV_PosHelTheta010->GetBinContent(i+1);
     }

     // Define our new histogram to plot the asymmetry here
     TCanvas *canvas19 = new TCanvas("canvas19","canvas19", 1920, 1080);
     TH1D *AsymmHist525Theta010 = new TH1D("AsymmHist525Theta010", "Asymmetry in PhiSc Between -ve and +ve Helicity (525pm25MeV, ThetaSc 0-10)", 90, -180, 180);

     for (Int_t i = 0; i < PosHelnBins; i++){
     	 DiffBinValues[i] = NegHelBinValues[i] - PosHelBinValues[i];
	 SumBinValues[i] = NegHelBinValues[i] + PosHelBinValues[i];
	 AsymmetryValues[i] = DiffBinValues[i]/SumBinValues[i];
	 AsymmHist525Theta010->SetBinContent(i+1, AsymmetryValues[i]);
     }

     //Save our asymmetry histogram
     AsymmHist525Theta010->Draw();
     canvas19->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_525MeV_Theta010_Asymm_7.png");
     canvas19->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_525MeV_Theta010_Asymm_7.pdf");

     // ###############################################################################################################

     NegHelnBins = Phi_Scattered_575MeV_NegHelTheta010->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < NegHelnBins; i++){
     	 NegHelBinValues[i] = Phi_Scattered_575MeV_NegHelTheta010->GetBinContent(i+1);
     }

     PosHelnBins = Phi_Scattered_575MeV_PosHelTheta010->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < PosHelnBins; i++){
     	 PosHelBinValues[i] = Phi_Scattered_575MeV_PosHelTheta010->GetBinContent(i+1);
     }

     // Define our new histogram to plot the asymmetry here
     TCanvas *canvas20 = new TCanvas("canvas20","canvas20", 1920, 1080);
     TH1D *AsymmHist575Theta010 = new TH1D("AsymmHist575Theta010", "Asymmetry in PhiSc Between -ve and +ve Helicity (575pm25MeV, ThetaSc 0-10)", 90, -180, 180);

     for (Int_t i = 0; i < PosHelnBins; i++){
     	 DiffBinValues[i] = NegHelBinValues[i] - PosHelBinValues[i];
	 SumBinValues[i] = NegHelBinValues[i] + PosHelBinValues[i];
	 AsymmetryValues[i] = DiffBinValues[i]/SumBinValues[i];
	 AsymmHist575Theta010->SetBinContent(i+1, AsymmetryValues[i]);
     }

     //Save our asymmetry histogram
     AsymmHist575Theta010->Draw();
     canvas20->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_575MeV_Theta010_Asymm_7.png");
     canvas20->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_575MeV_Theta010_Asymm_7.pdf");

     // ###############################################################################################################

     NegHelnBins = Phi_Scattered_625MeV_NegHelTheta010->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < NegHelnBins; i++){
     	 NegHelBinValues[i] = Phi_Scattered_625MeV_NegHelTheta010->GetBinContent(i+1);
     }

     PosHelnBins = Phi_Scattered_625MeV_PosHelTheta010->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < PosHelnBins; i++){
     	 PosHelBinValues[i] = Phi_Scattered_625MeV_PosHelTheta010->GetBinContent(i+1);
     }

     // Define our new histogram to plot the asymmetry here
     TCanvas *canvas21 = new TCanvas("canvas21","canvas21", 1920, 1080);
     TH1D *AsymmHist625Theta010 = new TH1D("AsymmHist625Theta010", "Asymmetry in PhiSc Between -ve and +ve Helicity (625pm25MeV, ThetaSc 0-10)", 90, -180, 180);

     for (Int_t i = 0; i < PosHelnBins; i++){
     	 DiffBinValues[i] = NegHelBinValues[i] - PosHelBinValues[i];
	 SumBinValues[i] = NegHelBinValues[i] + PosHelBinValues[i];
	 AsymmetryValues[i] = DiffBinValues[i]/SumBinValues[i];
	 AsymmHist625Theta010->SetBinContent(i+1, AsymmetryValues[i]);
     }

     //Save our asymmetry histogram
     AsymmHist625Theta010->Draw();
     canvas21->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_625MeV_Theta010_Asymm_7.png");
     canvas21->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_625MeV_Theta010_Asymm_7.pdf");

     // ###############################################################################################################

     NegHelnBins = Phi_Scattered_675MeV_NegHelTheta010->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < NegHelnBins; i++){
     	 NegHelBinValues[i] = Phi_Scattered_675MeV_NegHelTheta010->GetBinContent(i+1);
     }

     PosHelnBins = Phi_Scattered_675MeV_PosHelTheta010->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < PosHelnBins; i++){
     	 PosHelBinValues[i] = Phi_Scattered_675MeV_PosHelTheta010->GetBinContent(i+1);
     }

     // Define our new histogram to plot the asymmetry here
     TCanvas *canvas22 = new TCanvas("canvas22","canvas22", 1920, 1080);
     TH1D *AsymmHist675Theta010 = new TH1D("AsymmHist675Theta010", "Asymmetry in PhiSc Between -ve and +ve Helicity (675pm25MeV, ThetaSc 0-10)", 90, -180, 180);

     for (Int_t i = 0; i < PosHelnBins; i++){
     	 DiffBinValues[i] = NegHelBinValues[i] - PosHelBinValues[i];
	 SumBinValues[i] = NegHelBinValues[i] + PosHelBinValues[i];
	 AsymmetryValues[i] = DiffBinValues[i]/SumBinValues[i];
	 AsymmHist675Theta010->SetBinContent(i+1, AsymmetryValues[i]);
     }

     //Save our asymmetry histogram
     AsymmHist675Theta010->Draw();
     canvas22->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_675MeV_Theta010_Asymm_7.png");
     canvas22->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_675MeV_Theta010_Asymm_7.pdf");
    
     // ###############################################################################################################

     NegHelnBins = Phi_Scattered_725MeV_NegHelTheta010->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < NegHelnBins; i++){
     	 NegHelBinValues[i] = Phi_Scattered_725MeV_NegHelTheta010->GetBinContent(i+1);
     }

     PosHelnBins = Phi_Scattered_725MeV_PosHelTheta010->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < PosHelnBins; i++){
     	 PosHelBinValues[i] = Phi_Scattered_725MeV_PosHelTheta010->GetBinContent(i+1);
     }

     // Define our new histogram to plot the asymmetry here
     TCanvas *canvas23 = new TCanvas("canvas23","canvas23", 1920, 1080);
     TH1D *AsymmHist725Theta010 = new TH1D("AsymmHist725Theta010", "Asymmetry in PhiSc Between -ve and +ve Helicity (725pm25MeV, ThetaSc 0-10)", 90, -180, 180);

     for (Int_t i = 0; i < PosHelnBins; i++){
     	 DiffBinValues[i] = NegHelBinValues[i] - PosHelBinValues[i];
	 SumBinValues[i] = NegHelBinValues[i] + PosHelBinValues[i];
	 AsymmetryValues[i] = DiffBinValues[i]/SumBinValues[i];
	 AsymmHist725Theta010->SetBinContent(i+1, AsymmetryValues[i]);
     }

     //Save our asymmetry histogram
     AsymmHist725Theta010->Draw();
     canvas23->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_725MeV_Theta010_Asymm_7.png");
     canvas23->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_725MeV_Theta010_Asymm_7.pdf");

     // ###############################################################################################################

     NegHelnBins = Phi_Scattered_775MeV_NegHelTheta010->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < NegHelnBins; i++){
     	 NegHelBinValues[i] = Phi_Scattered_775MeV_NegHelTheta010->GetBinContent(i+1);
     }

     PosHelnBins = Phi_Scattered_775MeV_PosHelTheta010->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < PosHelnBins; i++){
     	 PosHelBinValues[i] = Phi_Scattered_775MeV_PosHelTheta010->GetBinContent(i+1);
     }

     // Define our new histogram to plot the asymmetry here
     TCanvas *canvas24 = new TCanvas("canvas24","canvas24", 1920, 1080);
     TH1D *AsymmHist775Theta010 = new TH1D("AsymmHist775Theta010", "Asymmetry in PhiSc Between -ve and +ve Helicity (775pm25MeV, ThetaSc 0-10)", 90, -180, 180);

     for (Int_t i = 0; i < PosHelnBins; i++){
     	 DiffBinValues[i] = NegHelBinValues[i] - PosHelBinValues[i];
	 SumBinValues[i] = NegHelBinValues[i] + PosHelBinValues[i];
	 AsymmetryValues[i] = DiffBinValues[i]/SumBinValues[i];
	 AsymmHist775Theta010->SetBinContent(i+1, AsymmetryValues[i]);
     }

     //Save our asymmetry histogram
     AsymmHist775Theta010->Draw();
     canvas24->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_775MeV_Theta010_Asymm_7.png");
     canvas24->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_775MeV_Theta010_Asymm_7.pdf");

     // ###############################################################################################################

     NegHelnBins = Phi_Scattered_825MeV_NegHelTheta010->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < NegHelnBins; i++){
     	 NegHelBinValues[i] = Phi_Scattered_825MeV_NegHelTheta010->GetBinContent(i+1);
     }

     PosHelnBins = Phi_Scattered_825MeV_PosHelTheta010->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < PosHelnBins; i++){
     	 PosHelBinValues[i] = Phi_Scattered_825MeV_PosHelTheta010->GetBinContent(i+1);
     }

     // Define our new histogram to plot the asymmetry here
     TCanvas *canvas25 = new TCanvas("canvas25","canvas25", 1920, 1080);
     TH1D *AsymmHist825Theta010 = new TH1D("AsymmHist825Theta010", "Asymmetry in PhiSc Between -ve and +ve Helicity (825pm25MeV, ThetaSc 0-10)", 90, -180, 180);

     for (Int_t i = 0; i < PosHelnBins; i++){
     	 DiffBinValues[i] = NegHelBinValues[i] - PosHelBinValues[i];
	 SumBinValues[i] = NegHelBinValues[i] + PosHelBinValues[i];
	 AsymmetryValues[i] = DiffBinValues[i]/SumBinValues[i];
	 AsymmHist825Theta010->SetBinContent(i+1, AsymmetryValues[i]);
     }

     //Save our asymmetry histogram
     AsymmHist825Theta010->Draw();
     canvas25->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_825MeV_Theta010_Asymm_7.png");
     canvas25->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_825MeV_Theta010_Asymm_7.pdf");

     // ###############################################################################################################

     NegHelnBins = Phi_Scattered_875MeV_NegHelTheta010->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < NegHelnBins; i++){
     	 NegHelBinValues[i] = Phi_Scattered_875MeV_NegHelTheta010->GetBinContent(i+1);
     }

     PosHelnBins = Phi_Scattered_875MeV_PosHelTheta010->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < PosHelnBins; i++){
     	 PosHelBinValues[i] = Phi_Scattered_875MeV_PosHelTheta010->GetBinContent(i+1);
     }

     // Define our new histogram to plot the asymmetry here
     TCanvas *canvas26 = new TCanvas("canvas26","canvas26", 1920, 1080);
     TH1D *AsymmHist875Theta010 = new TH1D("AsymmHist875Theta010", "Asymmetry in PhiSc Between -ve and +ve Helicity (875pm25MeV, ThetaSc 0-10)", 90, -180, 180);

     for (Int_t i = 0; i < PosHelnBins; i++){
     	 DiffBinValues[i] = NegHelBinValues[i] - PosHelBinValues[i];
	 SumBinValues[i] = NegHelBinValues[i] + PosHelBinValues[i];
	 AsymmetryValues[i] = DiffBinValues[i]/SumBinValues[i];
	 AsymmHist875Theta010->SetBinContent(i+1, AsymmetryValues[i]);
     }

     //Save our asymmetry histogram
     AsymmHist875Theta010->Draw();
     canvas26->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_875MeV_Theta010_Asymm_7.png");
     canvas26->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_875MeV_Theta010_Asymm_7.pdf");

     // #######################################################################################################################
     // #######################################################################################################################
     // 10-20 Degree Plots
     // #######################################################################################################################
     // #######################################################################################################################

     NegHelnBins = Phi_Scattered_275MeV_NegHelTheta1020->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < NegHelnBins; i++){
     	 NegHelBinValues[i] = Phi_Scattered_275MeV_NegHelTheta1020->GetBinContent(i+1);
     }

     PosHelnBins = Phi_Scattered_275MeV_PosHelTheta1020->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < PosHelnBins; i++){
     	 PosHelBinValues[i] = Phi_Scattered_275MeV_PosHelTheta1020->GetBinContent(i+1);
     }

     // Define our new histogram to plot the asymmetry here
     TCanvas *canvas27 = new TCanvas("canvas27","canvas27", 1920, 1080);
     TH1D *AsymmHist275Theta1020 = new TH1D("AsymmHist275Theta1020", "Asymmetry in PhiSc Between -ve and +ve Helicity (275pm25MeV, ThetaSc 10-20)", 90, -180, 180);

     for (Int_t i = 0; i < PosHelnBins; i++){
     	 DiffBinValues[i] = NegHelBinValues[i] - PosHelBinValues[i];
	 SumBinValues[i] = NegHelBinValues[i] + PosHelBinValues[i];
	 AsymmetryValues[i] = DiffBinValues[i]/SumBinValues[i];
	 AsymmHist275Theta1020->SetBinContent(i+1, AsymmetryValues[i]);
     }

     //Save our asymmetry histogram
     AsymmHist275Theta1020->Draw();
     canvas27->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_275MeV_Theta1020_Asymm_7.png");
     canvas27->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_275MeV_Theta1020_Asymm_7.pdf");

     // ###############################################################################################################

     NegHelnBins = Phi_Scattered_325MeV_NegHelTheta1020->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < NegHelnBins; i++){
     	 NegHelBinValues[i] = Phi_Scattered_325MeV_NegHelTheta1020->GetBinContent(i+1);
     }

     PosHelnBins = Phi_Scattered_325MeV_PosHelTheta1020->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < PosHelnBins; i++){
     	 PosHelBinValues[i] = Phi_Scattered_325MeV_PosHelTheta1020->GetBinContent(i+1);
     }

     // Define our new histogram to plot the asymmetry here
     TCanvas *canvas28 = new TCanvas("canvas28","canvas28", 1920, 1080);
     TH1D *AsymmHist325Theta1020 = new TH1D("AsymmHist325Theta1020", "Asymmetry in PhiSc Between -ve and +ve Helicity (325pm25MeV, ThetaSc 10-20)", 90, -180, 180);

     for (Int_t i = 0; i < PosHelnBins; i++){
     	 DiffBinValues[i] = NegHelBinValues[i] - PosHelBinValues[i];
	 SumBinValues[i] = NegHelBinValues[i] + PosHelBinValues[i];
	 AsymmetryValues[i] = DiffBinValues[i]/SumBinValues[i];
	 AsymmHist325Theta1020->SetBinContent(i+1, AsymmetryValues[i]);
     }

     //Save our asymmetry histogram
     AsymmHist325Theta1020->Draw();
     canvas28->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_325MeV_Theta1020_Asymm_7.png");
     canvas28->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_325MeV_Theta1020_Asymm_7.pdf");

     // ###############################################################################################################

     NegHelnBins = Phi_Scattered_375MeV_NegHelTheta1020->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < NegHelnBins; i++){
     	 NegHelBinValues[i] = Phi_Scattered_375MeV_NegHelTheta1020->GetBinContent(i+1);
     }

     PosHelnBins = Phi_Scattered_375MeV_PosHelTheta1020->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < PosHelnBins; i++){
     	 PosHelBinValues[i] = Phi_Scattered_375MeV_PosHelTheta1020->GetBinContent(i+1);
     }

     // Define our new histogram to plot the asymmetry here
     TCanvas *canvas29 = new TCanvas("canvas29","canvas29", 1920, 1080);
     TH1D *AsymmHist375Theta1020 = new TH1D("AsymmHist375Theta1020", "Asymmetry in PhiSc Between -ve and +ve Helicity (375pm25MeV, ThetaSc 10-20)", 90, -180, 180);

     for (Int_t i = 0; i < PosHelnBins; i++){
     	 DiffBinValues[i] = NegHelBinValues[i] - PosHelBinValues[i];
	 SumBinValues[i] = NegHelBinValues[i] + PosHelBinValues[i];
	 AsymmetryValues[i] = DiffBinValues[i]/SumBinValues[i];
	 AsymmHist375Theta1020->SetBinContent(i+1, AsymmetryValues[i]);
     }

     //Save our asymmetry histogram
     AsymmHist375Theta1020->Draw();
     canvas29->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_375MeV_Theta1020_Asymm_7.png");
     canvas29->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_375MeV_Theta1020_Asymm_7.pdf");

     // ###############################################################################################################

     NegHelnBins = Phi_Scattered_425MeV_NegHelTheta1020->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < NegHelnBins; i++){
     	 NegHelBinValues[i] = Phi_Scattered_425MeV_NegHelTheta1020->GetBinContent(i+1);
     }

     PosHelnBins = Phi_Scattered_425MeV_PosHelTheta1020->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < PosHelnBins; i++){
     	 PosHelBinValues[i] = Phi_Scattered_425MeV_PosHelTheta1020->GetBinContent(i+1);
     }

     // Define our new histogram to plot the asymmetry here
     TCanvas *canvas30 = new TCanvas("canvas30","canvas30", 1920, 1080);
     TH1D *AsymmHist425Theta1020 = new TH1D("AsymmHist425Theta1020", "Asymmetry in PhiSc Between -ve and +ve Helicity (425pm25MeV, ThetaSc 10-20)", 90, -180, 180);

     for (Int_t i = 0; i < PosHelnBins; i++){
     	 DiffBinValues[i] = NegHelBinValues[i] - PosHelBinValues[i];
	 SumBinValues[i] = NegHelBinValues[i] + PosHelBinValues[i];
	 AsymmetryValues[i] = DiffBinValues[i]/SumBinValues[i];
	 AsymmHist425Theta1020->SetBinContent(i+1, AsymmetryValues[i]);
     }

     //Save our asymmetry histogram
     AsymmHist425Theta1020->Draw();
     canvas30->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_425MeV_Theta1020_Asymm_7.png");
     canvas30->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_425MeV_Theta1020_Asymm_7.pdf");

     // ###############################################################################################################

     NegHelnBins = Phi_Scattered_475MeV_NegHelTheta1020->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < NegHelnBins; i++){
     	 NegHelBinValues[i] = Phi_Scattered_475MeV_NegHelTheta1020->GetBinContent(i+1);
     }

     PosHelnBins = Phi_Scattered_475MeV_PosHelTheta1020->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < PosHelnBins; i++){
     	 PosHelBinValues[i] = Phi_Scattered_475MeV_PosHelTheta1020->GetBinContent(i+1);
     }

     // Define our new histogram to plot the asymmetry here
     TCanvas *canvas31 = new TCanvas("canvas31","canvas31", 1920, 1080);
     TH1D *AsymmHist475Theta1020 = new TH1D("AsymmHist475Theta1020", "Asymmetry in PhiSc Between -ve and +ve Helicity (475pm25MeV, ThetaSc 10-20)", 90, -180, 180);

     for (Int_t i = 0; i < PosHelnBins; i++){
     	 DiffBinValues[i] = NegHelBinValues[i] - PosHelBinValues[i];
	 SumBinValues[i] = NegHelBinValues[i] + PosHelBinValues[i];
	 AsymmetryValues[i] = DiffBinValues[i]/SumBinValues[i];
	 AsymmHist475Theta1020->SetBinContent(i+1, AsymmetryValues[i]);
     }

     //Save our asymmetry histogram
     AsymmHist475Theta1020->Draw();
     canvas31->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_475MeV_Theta1020_Asymm_7.png");
     canvas31->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_475MeV_Theta1020_Asymm_7.pdf");

     // ###############################################################################################################

     NegHelnBins = Phi_Scattered_525MeV_NegHelTheta1020->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < NegHelnBins; i++){
     	 NegHelBinValues[i] = Phi_Scattered_525MeV_NegHelTheta1020->GetBinContent(i+1);
     }

     PosHelnBins = Phi_Scattered_525MeV_PosHelTheta1020->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < PosHelnBins; i++){
     	 PosHelBinValues[i] = Phi_Scattered_525MeV_PosHelTheta1020->GetBinContent(i+1);
     }

     // Define our new histogram to plot the asymmetry here
     TCanvas *canvas32 = new TCanvas("canvas32","canvas32", 1920, 1080);
     TH1D *AsymmHist525Theta1020 = new TH1D("AsymmHist525Theta1020", "Asymmetry in PhiSc Between -ve and +ve Helicity (525pm25MeV, ThetaSc 10-20)", 90, -180, 180);

     for (Int_t i = 0; i < PosHelnBins; i++){
     	 DiffBinValues[i] = NegHelBinValues[i] - PosHelBinValues[i];
	 SumBinValues[i] = NegHelBinValues[i] + PosHelBinValues[i];
	 AsymmetryValues[i] = DiffBinValues[i]/SumBinValues[i];
	 AsymmHist525Theta1020->SetBinContent(i+1, AsymmetryValues[i]);
     }

     //Save our asymmetry histogram
     AsymmHist525Theta1020->Draw();
     canvas32->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_525MeV_Theta1020_Asymm_7.png");
     canvas32->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_525MeV_Theta1020_Asymm_7.pdf");

     // ###############################################################################################################

     NegHelnBins = Phi_Scattered_575MeV_NegHelTheta1020->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < NegHelnBins; i++){
     	 NegHelBinValues[i] = Phi_Scattered_575MeV_NegHelTheta1020->GetBinContent(i+1);
     }

     PosHelnBins = Phi_Scattered_575MeV_PosHelTheta1020->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < PosHelnBins; i++){
     	 PosHelBinValues[i] = Phi_Scattered_575MeV_PosHelTheta1020->GetBinContent(i+1);
     }

     // Define our new histogram to plot the asymmetry here
     TCanvas *canvas33 = new TCanvas("canvas33","canvas33", 1920, 1080);
     TH1D *AsymmHist575Theta1020 = new TH1D("AsymmHist575Theta1020", "Asymmetry in PhiSc Between -ve and +ve Helicity (575pm25MeV, ThetaSc 10-20)", 90, -180, 180);

     for (Int_t i = 0; i < PosHelnBins; i++){
     	 DiffBinValues[i] = NegHelBinValues[i] - PosHelBinValues[i];
	 SumBinValues[i] = NegHelBinValues[i] + PosHelBinValues[i];
	 AsymmetryValues[i] = DiffBinValues[i]/SumBinValues[i];
	 AsymmHist575Theta1020->SetBinContent(i+1, AsymmetryValues[i]);
     }

     //Save our asymmetry histogram
     AsymmHist575Theta1020->Draw();
     canvas33->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_575MeV_Theta1020_Asymm_7.png");
     canvas33->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_575MeV_Theta1020_Asymm_7.pdf");

     // ###############################################################################################################

     NegHelnBins = Phi_Scattered_625MeV_NegHelTheta1020->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < NegHelnBins; i++){
     	 NegHelBinValues[i] = Phi_Scattered_625MeV_NegHelTheta1020->GetBinContent(i+1);
     }

     PosHelnBins = Phi_Scattered_625MeV_PosHelTheta1020->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < PosHelnBins; i++){
     	 PosHelBinValues[i] = Phi_Scattered_625MeV_PosHelTheta1020->GetBinContent(i+1);
     }

     // Define our new histogram to plot the asymmetry here
     TCanvas *canvas34 = new TCanvas("canvas34","canvas34", 1920, 1080);
     TH1D *AsymmHist625Theta1020 = new TH1D("AsymmHist625Theta1020", "Asymmetry in PhiSc Between -ve and +ve Helicity (625pm25MeV, ThetaSc 10-20)", 90, -180, 180);

     for (Int_t i = 0; i < PosHelnBins; i++){
     	 DiffBinValues[i] = NegHelBinValues[i] - PosHelBinValues[i];
	 SumBinValues[i] = NegHelBinValues[i] + PosHelBinValues[i];
	 AsymmetryValues[i] = DiffBinValues[i]/SumBinValues[i];
	 AsymmHist625Theta1020->SetBinContent(i+1, AsymmetryValues[i]);
     }

     //Save our asymmetry histogram
     AsymmHist625Theta1020->Draw();
     canvas34->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_625MeV_Theta1020_Asymm_7.png");
     canvas34->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_625MeV_Theta1020_Asymm_7.pdf");

     // ###############################################################################################################

     NegHelnBins = Phi_Scattered_675MeV_NegHelTheta1020->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < NegHelnBins; i++){
     	 NegHelBinValues[i] = Phi_Scattered_675MeV_NegHelTheta1020->GetBinContent(i+1);
     }

     PosHelnBins = Phi_Scattered_675MeV_PosHelTheta1020->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < PosHelnBins; i++){
     	 PosHelBinValues[i] = Phi_Scattered_675MeV_PosHelTheta1020->GetBinContent(i+1);
     }

     // Define our new histogram to plot the asymmetry here
     TCanvas *canvas35 = new TCanvas("canvas35","canvas35", 1920, 1080);
     TH1D *AsymmHist675Theta1020 = new TH1D("AsymmHist675Theta1020", "Asymmetry in PhiSc Between -ve and +ve Helicity (675pm25MeV, ThetaSc 10-20)", 90, -180, 180);

     for (Int_t i = 0; i < PosHelnBins; i++){
     	 DiffBinValues[i] = NegHelBinValues[i] - PosHelBinValues[i];
	 SumBinValues[i] = NegHelBinValues[i] + PosHelBinValues[i];
	 AsymmetryValues[i] = DiffBinValues[i]/SumBinValues[i];
	 AsymmHist675Theta1020->SetBinContent(i+1, AsymmetryValues[i]);
     }

     //Save our asymmetry histogram
     AsymmHist675Theta1020->Draw();
     canvas35->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_675MeV_Theta1020_Asymm_7.png");
     canvas35->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_675MeV_Theta1020_Asymm_7.pdf");
    
     // ###############################################################################################################

     NegHelnBins = Phi_Scattered_725MeV_NegHelTheta1020->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < NegHelnBins; i++){
     	 NegHelBinValues[i] = Phi_Scattered_725MeV_NegHelTheta1020->GetBinContent(i+1);
     }

     PosHelnBins = Phi_Scattered_725MeV_PosHelTheta1020->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < PosHelnBins; i++){
     	 PosHelBinValues[i] = Phi_Scattered_725MeV_PosHelTheta1020->GetBinContent(i+1);
     }

     // Define our new histogram to plot the asymmetry here
     TCanvas *canvas36 = new TCanvas("canvas36","canvas36", 1920, 1080);
     TH1D *AsymmHist725Theta1020 = new TH1D("AsymmHist725Theta1020", "Asymmetry in PhiSc Between -ve and +ve Helicity (725pm25MeV, ThetaSc 10-20)", 90, -180, 180);

     for (Int_t i = 0; i < PosHelnBins; i++){
     	 DiffBinValues[i] = NegHelBinValues[i] - PosHelBinValues[i];
	 SumBinValues[i] = NegHelBinValues[i] + PosHelBinValues[i];
	 AsymmetryValues[i] = DiffBinValues[i]/SumBinValues[i];
	 AsymmHist725Theta1020->SetBinContent(i+1, AsymmetryValues[i]);
     }

     //Save our asymmetry histogram
     AsymmHist725Theta1020->Draw();
     canvas36->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_725MeV_Theta1020_Asymm_7.png");
     canvas36->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_725MeV_Theta1020_Asymm_7.pdf");

     // ###############################################################################################################

     NegHelnBins = Phi_Scattered_775MeV_NegHelTheta1020->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < NegHelnBins; i++){
     	 NegHelBinValues[i] = Phi_Scattered_775MeV_NegHelTheta1020->GetBinContent(i+1);
     }

     PosHelnBins = Phi_Scattered_775MeV_PosHelTheta1020->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < PosHelnBins; i++){
     	 PosHelBinValues[i] = Phi_Scattered_775MeV_PosHelTheta1020->GetBinContent(i+1);
     }

     // Define our new histogram to plot the asymmetry here
     TCanvas *canvas37 = new TCanvas("canvas37","canvas37", 1920, 1080);
     TH1D *AsymmHist775Theta1020 = new TH1D("AsymmHist775Theta1020", "Asymmetry in PhiSc Between -ve and +ve Helicity (775pm25MeV, ThetaSc 10-20)", 90, -180, 180);

     for (Int_t i = 0; i < PosHelnBins; i++){
     	 DiffBinValues[i] = NegHelBinValues[i] - PosHelBinValues[i];
	 SumBinValues[i] = NegHelBinValues[i] + PosHelBinValues[i];
	 AsymmetryValues[i] = DiffBinValues[i]/SumBinValues[i];
	 AsymmHist775Theta1020->SetBinContent(i+1, AsymmetryValues[i]);
     }

     //Save our asymmetry histogram
     AsymmHist775Theta1020->Draw();
     canvas37->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_775MeV_Theta1020_Asymm_7.png");
     canvas37->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_775MeV_Theta1020_Asymm_7.pdf");

     // ###############################################################################################################

     NegHelnBins = Phi_Scattered_825MeV_NegHelTheta1020->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < NegHelnBins; i++){
     	 NegHelBinValues[i] = Phi_Scattered_825MeV_NegHelTheta1020->GetBinContent(i+1);
     }

     PosHelnBins = Phi_Scattered_825MeV_PosHelTheta1020->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < PosHelnBins; i++){
     	 PosHelBinValues[i] = Phi_Scattered_825MeV_PosHelTheta1020->GetBinContent(i+1);
     }

     // Define our new histogram to plot the asymmetry here
     TCanvas *canvas38 = new TCanvas("canvas38","canvas38", 1920, 1080);
     TH1D *AsymmHist825Theta1020 = new TH1D("AsymmHist825Theta1020", "Asymmetry in PhiSc Between -ve and +ve Helicity (825pm25MeV, ThetaSc 10-20)", 90, -180, 180);

     for (Int_t i = 0; i < PosHelnBins; i++){
     	 DiffBinValues[i] = NegHelBinValues[i] - PosHelBinValues[i];
	 SumBinValues[i] = NegHelBinValues[i] + PosHelBinValues[i];
	 AsymmetryValues[i] = DiffBinValues[i]/SumBinValues[i];
	 AsymmHist825Theta1020->SetBinContent(i+1, AsymmetryValues[i]);
     }

     //Save our asymmetry histogram
     AsymmHist825Theta1020->Draw();
     canvas38->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_825MeV_Theta1020_Asymm_7.png");
     canvas38->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_825MeV_Theta1020_Asymm_7.pdf");

     // ###############################################################################################################

     NegHelnBins = Phi_Scattered_875MeV_NegHelTheta1020->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < NegHelnBins; i++){
     	 NegHelBinValues[i] = Phi_Scattered_875MeV_NegHelTheta1020->GetBinContent(i+1);
     }

     PosHelnBins = Phi_Scattered_875MeV_PosHelTheta1020->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < PosHelnBins; i++){
     	 PosHelBinValues[i] = Phi_Scattered_875MeV_PosHelTheta1020->GetBinContent(i+1);
     }

     // Define our new histogram to plot the asymmetry here
     TCanvas *canvas39 = new TCanvas("canvas39","canvas39", 1920, 1080);
     TH1D *AsymmHist875Theta1020 = new TH1D("AsymmHist875Theta1020", "Asymmetry in PhiSc Between -ve and +ve Helicity (875pm25MeV, ThetaSc 10-20)", 90, -180, 180);

     for (Int_t i = 0; i < PosHelnBins; i++){
     	 DiffBinValues[i] = NegHelBinValues[i] - PosHelBinValues[i];
	 SumBinValues[i] = NegHelBinValues[i] + PosHelBinValues[i];
	 AsymmetryValues[i] = DiffBinValues[i]/SumBinValues[i];
	 AsymmHist875Theta1020->SetBinContent(i+1, AsymmetryValues[i]);
     }

     //Save our asymmetry histogram
     AsymmHist875Theta1020->Draw();
     canvas39->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_875MeV_Theta1020_Asymm_7.png");
     canvas39->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_875MeV_Theta1020_Asymm_7.pdf");

     // #######################################################################################################################
     // #######################################################################################################################
     // 20-30 Degree Plots
     // #######################################################################################################################
     // #######################################################################################################################

     NegHelnBins = Phi_Scattered_275MeV_NegHelTheta2030->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < NegHelnBins; i++){
     	 NegHelBinValues[i] = Phi_Scattered_275MeV_NegHelTheta2030->GetBinContent(i+1);
     }

     PosHelnBins = Phi_Scattered_275MeV_PosHelTheta2030->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < PosHelnBins; i++){
     	 PosHelBinValues[i] = Phi_Scattered_275MeV_PosHelTheta2030->GetBinContent(i+1);
     }

     // Define our new histogram to plot the asymmetry here
     TCanvas *canvas40 = new TCanvas("canvas40","canvas40", 1920, 1080);
     TH1D *AsymmHist275Theta2030 = new TH1D("AsymmHist275Theta2030", "Asymmetry in PhiSc Between -ve and +ve Helicity (275pm25MeV, ThetaSc 20-30)", 90, -180, 180);

     for (Int_t i = 0; i < PosHelnBins; i++){
     	 DiffBinValues[i] = NegHelBinValues[i] - PosHelBinValues[i];
	 SumBinValues[i] = NegHelBinValues[i] + PosHelBinValues[i];
	 AsymmetryValues[i] = DiffBinValues[i]/SumBinValues[i];
	 AsymmHist275Theta2030->SetBinContent(i+1, AsymmetryValues[i]);
     }

     //Save our asymmetry histogram
     AsymmHist275Theta2030->Draw();
     canvas40->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_275MeV_Theta2030_Asymm_7.png");
     canvas40->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_275MeV_Theta2030_Asymm_7.pdf");

     // ###############################################################################################################

     NegHelnBins = Phi_Scattered_325MeV_NegHelTheta2030->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < NegHelnBins; i++){
     	 NegHelBinValues[i] = Phi_Scattered_325MeV_NegHelTheta2030->GetBinContent(i+1);
     }

     PosHelnBins = Phi_Scattered_325MeV_PosHelTheta2030->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < PosHelnBins; i++){
     	 PosHelBinValues[i] = Phi_Scattered_325MeV_PosHelTheta2030->GetBinContent(i+1);
     }

     // Define our new histogram to plot the asymmetry here
     TCanvas *canvas41 = new TCanvas("canvas41","canvas41", 1920, 1080);
     TH1D *AsymmHist325Theta2030 = new TH1D("AsymmHist325Theta2030", "Asymmetry in PhiSc Between -ve and +ve Helicity (325pm25MeV, ThetaSc 20-30)", 90, -180, 180);

     for (Int_t i = 0; i < PosHelnBins; i++){
     	 DiffBinValues[i] = NegHelBinValues[i] - PosHelBinValues[i];
	 SumBinValues[i] = NegHelBinValues[i] + PosHelBinValues[i];
	 AsymmetryValues[i] = DiffBinValues[i]/SumBinValues[i];
	 AsymmHist325Theta2030->SetBinContent(i+1, AsymmetryValues[i]);
     }

     //Save our asymmetry histogram
     AsymmHist325Theta2030->Draw();
     canvas41->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_325MeV_Theta2030_Asymm_7.png");
     canvas41->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_325MeV_Theta2030_Asymm_7.pdf");

     // ###############################################################################################################

     NegHelnBins = Phi_Scattered_375MeV_NegHelTheta2030->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < NegHelnBins; i++){
     	 NegHelBinValues[i] = Phi_Scattered_375MeV_NegHelTheta2030->GetBinContent(i+1);
     }

     PosHelnBins = Phi_Scattered_375MeV_PosHelTheta2030->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < PosHelnBins; i++){
     	 PosHelBinValues[i] = Phi_Scattered_375MeV_PosHelTheta2030->GetBinContent(i+1);
     }

     // Define our new histogram to plot the asymmetry here
     TCanvas *canvas42 = new TCanvas("canvas42","canvas42", 1920, 1080);
     TH1D *AsymmHist375Theta2030 = new TH1D("AsymmHist375Theta2030", "Asymmetry in PhiSc Between -ve and +ve Helicity (375pm25MeV, ThetaSc 20-30)", 90, -180, 180);

     for (Int_t i = 0; i < PosHelnBins; i++){
     	 DiffBinValues[i] = NegHelBinValues[i] - PosHelBinValues[i];
	 SumBinValues[i] = NegHelBinValues[i] + PosHelBinValues[i];
	 AsymmetryValues[i] = DiffBinValues[i]/SumBinValues[i];
	 AsymmHist375Theta2030->SetBinContent(i+1, AsymmetryValues[i]);
     }

     //Save our asymmetry histogram
     AsymmHist375Theta2030->Draw();
     canvas42->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_375MeV_Theta2030_Asymm_7.png");
     canvas42->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_375MeV_Theta2030_Asymm_7.pdf");

     // ###############################################################################################################

     NegHelnBins = Phi_Scattered_425MeV_NegHelTheta2030->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < NegHelnBins; i++){
     	 NegHelBinValues[i] = Phi_Scattered_425MeV_NegHelTheta2030->GetBinContent(i+1);
     }

     PosHelnBins = Phi_Scattered_425MeV_PosHelTheta2030->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < PosHelnBins; i++){
     	 PosHelBinValues[i] = Phi_Scattered_425MeV_PosHelTheta2030->GetBinContent(i+1);
     }

     // Define our new histogram to plot the asymmetry here
     TCanvas *canvas43 = new TCanvas("canvas43","canvas43", 1920, 1080);
     TH1D *AsymmHist425Theta2030 = new TH1D("AsymmHist425Theta2030", "Asymmetry in PhiSc Between -ve and +ve Helicity (425pm25MeV, ThetaSc 20-30)", 90, -180, 180);

     for (Int_t i = 0; i < PosHelnBins; i++){
     	 DiffBinValues[i] = NegHelBinValues[i] - PosHelBinValues[i];
	 SumBinValues[i] = NegHelBinValues[i] + PosHelBinValues[i];
	 AsymmetryValues[i] = DiffBinValues[i]/SumBinValues[i];
	 AsymmHist425Theta2030->SetBinContent(i+1, AsymmetryValues[i]);
     }

     //Save our asymmetry histogram
     AsymmHist425Theta2030->Draw();
     canvas43->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_425MeV_Theta2030_Asymm_7.png");
     canvas43->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_425MeV_Theta2030_Asymm_7.pdf");

     // ###############################################################################################################

     NegHelnBins = Phi_Scattered_475MeV_NegHelTheta2030->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < NegHelnBins; i++){
     	 NegHelBinValues[i] = Phi_Scattered_475MeV_NegHelTheta2030->GetBinContent(i+1);
     }

     PosHelnBins = Phi_Scattered_475MeV_PosHelTheta2030->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < PosHelnBins; i++){
     	 PosHelBinValues[i] = Phi_Scattered_475MeV_PosHelTheta2030->GetBinContent(i+1);
     }

     // Define our new histogram to plot the asymmetry here
     TCanvas *canvas44 = new TCanvas("canvas44","canvas44", 1920, 1080);
     TH1D *AsymmHist475Theta2030 = new TH1D("AsymmHist475Theta2030", "Asymmetry in PhiSc Between -ve and +ve Helicity (475pm25MeV, ThetaSc 20-30)", 90, -180, 180);

     for (Int_t i = 0; i < PosHelnBins; i++){
     	 DiffBinValues[i] = NegHelBinValues[i] - PosHelBinValues[i];
	 SumBinValues[i] = NegHelBinValues[i] + PosHelBinValues[i];
	 AsymmetryValues[i] = DiffBinValues[i]/SumBinValues[i];
	 AsymmHist475Theta2030->SetBinContent(i+1, AsymmetryValues[i]);
     }

     //Save our asymmetry histogram
     AsymmHist475Theta2030->Draw();
     canvas44->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_475MeV_Theta2030_Asymm_7.png");
     canvas44->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_475MeV_Theta2030_Asymm_7.pdf");

     // ###############################################################################################################

     NegHelnBins = Phi_Scattered_525MeV_NegHelTheta2030->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < NegHelnBins; i++){
     	 NegHelBinValues[i] = Phi_Scattered_525MeV_NegHelTheta2030->GetBinContent(i+1);
     }

     PosHelnBins = Phi_Scattered_525MeV_PosHelTheta2030->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < PosHelnBins; i++){
     	 PosHelBinValues[i] = Phi_Scattered_525MeV_PosHelTheta2030->GetBinContent(i+1);
     }

     // Define our new histogram to plot the asymmetry here
     TCanvas *canvas45 = new TCanvas("canvas45","canvas45", 1920, 1080);
     TH1D *AsymmHist525Theta2030 = new TH1D("AsymmHist525Theta2030", "Asymmetry in PhiSc Between -ve and +ve Helicity (525pm25MeV, ThetaSc 20-30)", 90, -180, 180);

     for (Int_t i = 0; i < PosHelnBins; i++){
     	 DiffBinValues[i] = NegHelBinValues[i] - PosHelBinValues[i];
	 SumBinValues[i] = NegHelBinValues[i] + PosHelBinValues[i];
	 AsymmetryValues[i] = DiffBinValues[i]/SumBinValues[i];
	 AsymmHist525Theta2030->SetBinContent(i+1, AsymmetryValues[i]);
     }

     //Save our asymmetry histogram
     AsymmHist525Theta2030->Draw();
     canvas45->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_525MeV_Theta2030_Asymm_7.png");
     canvas45->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_525MeV_Theta2030_Asymm_7.pdf");

     // ###############################################################################################################

     NegHelnBins = Phi_Scattered_575MeV_NegHelTheta2030->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < NegHelnBins; i++){
     	 NegHelBinValues[i] = Phi_Scattered_575MeV_NegHelTheta2030->GetBinContent(i+1);
     }

     PosHelnBins = Phi_Scattered_575MeV_PosHelTheta2030->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < PosHelnBins; i++){
     	 PosHelBinValues[i] = Phi_Scattered_575MeV_PosHelTheta2030->GetBinContent(i+1);
     }

     // Define our new histogram to plot the asymmetry here
     TCanvas *canvas46 = new TCanvas("canvas46","canvas46", 1920, 1080);
     TH1D *AsymmHist575Theta2030 = new TH1D("AsymmHist575Theta2030", "Asymmetry in PhiSc Between -ve and +ve Helicity (575pm25MeV, ThetaSc 20-30)", 90, -180, 180);

     for (Int_t i = 0; i < PosHelnBins; i++){
     	 DiffBinValues[i] = NegHelBinValues[i] - PosHelBinValues[i];
	 SumBinValues[i] = NegHelBinValues[i] + PosHelBinValues[i];
	 AsymmetryValues[i] = DiffBinValues[i]/SumBinValues[i];
	 AsymmHist575Theta2030->SetBinContent(i+1, AsymmetryValues[i]);
     }

     //Save our asymmetry histogram
     AsymmHist575Theta2030->Draw();
     canvas46->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_575MeV_Theta2030_Asymm_7.png");
     canvas46->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_575MeV_Theta2030_Asymm_7.pdf");

     // ###############################################################################################################

     NegHelnBins = Phi_Scattered_625MeV_NegHelTheta2030->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < NegHelnBins; i++){
     	 NegHelBinValues[i] = Phi_Scattered_625MeV_NegHelTheta2030->GetBinContent(i+1);
     }

     PosHelnBins = Phi_Scattered_625MeV_PosHelTheta2030->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < PosHelnBins; i++){
     	 PosHelBinValues[i] = Phi_Scattered_625MeV_PosHelTheta2030->GetBinContent(i+1);
     }

     // Define our new histogram to plot the asymmetry here
     TCanvas *canvas47 = new TCanvas("canvas47","canvas47", 1920, 1080);
     TH1D *AsymmHist625Theta2030 = new TH1D("AsymmHist625Theta2030", "Asymmetry in PhiSc Between -ve and +ve Helicity (625pm25MeV, ThetaSc 20-30)", 90, -180, 180);

     for (Int_t i = 0; i < PosHelnBins; i++){
     	 DiffBinValues[i] = NegHelBinValues[i] - PosHelBinValues[i];
	 SumBinValues[i] = NegHelBinValues[i] + PosHelBinValues[i];
	 AsymmetryValues[i] = DiffBinValues[i]/SumBinValues[i];
	 AsymmHist625Theta2030->SetBinContent(i+1, AsymmetryValues[i]);
     }

     //Save our asymmetry histogram
     AsymmHist625Theta2030->Draw();
     canvas47->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_625MeV_Theta2030_Asymm_7.png");
     canvas47->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_625MeV_Theta2030_Asymm_7.pdf");

     // ###############################################################################################################

     NegHelnBins = Phi_Scattered_675MeV_NegHelTheta2030->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < NegHelnBins; i++){
     	 NegHelBinValues[i] = Phi_Scattered_675MeV_NegHelTheta2030->GetBinContent(i+1);
     }

     PosHelnBins = Phi_Scattered_675MeV_PosHelTheta2030->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < PosHelnBins; i++){
     	 PosHelBinValues[i] = Phi_Scattered_675MeV_PosHelTheta2030->GetBinContent(i+1);
     }

     // Define our new histogram to plot the asymmetry here
     TCanvas *canvas48 = new TCanvas("canvas48","canvas48", 1920, 1080);
     TH1D *AsymmHist675Theta2030 = new TH1D("AsymmHist675Theta2030", "Asymmetry in PhiSc Between -ve and +ve Helicity (675pm25MeV, ThetaSc 20-30)", 90, -180, 180);

     for (Int_t i = 0; i < PosHelnBins; i++){
     	 DiffBinValues[i] = NegHelBinValues[i] - PosHelBinValues[i];
	 SumBinValues[i] = NegHelBinValues[i] + PosHelBinValues[i];
	 AsymmetryValues[i] = DiffBinValues[i]/SumBinValues[i];
	 AsymmHist675Theta2030->SetBinContent(i+1, AsymmetryValues[i]);
     }

     //Save our asymmetry histogram
     AsymmHist675Theta2030->Draw();
     canvas48->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_675MeV_Theta2030_Asymm_7.png");
     canvas48->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_675MeV_Theta2030_Asymm_7.pdf");
    
     // ###############################################################################################################

     NegHelnBins = Phi_Scattered_725MeV_NegHelTheta2030->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < NegHelnBins; i++){
     	 NegHelBinValues[i] = Phi_Scattered_725MeV_NegHelTheta2030->GetBinContent(i+1);
     }

     PosHelnBins = Phi_Scattered_725MeV_PosHelTheta2030->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < PosHelnBins; i++){
     	 PosHelBinValues[i] = Phi_Scattered_725MeV_PosHelTheta2030->GetBinContent(i+1);
     }

     // Define our new histogram to plot the asymmetry here
     TCanvas *canvas49 = new TCanvas("canvas49","canvas49", 1920, 1080);
     TH1D *AsymmHist725Theta2030 = new TH1D("AsymmHist725Theta2030", "Asymmetry in PhiSc Between -ve and +ve Helicity (725pm25MeV, ThetaSc 20-30)", 90, -180, 180);

     for (Int_t i = 0; i < PosHelnBins; i++){
     	 DiffBinValues[i] = NegHelBinValues[i] - PosHelBinValues[i];
	 SumBinValues[i] = NegHelBinValues[i] + PosHelBinValues[i];
	 AsymmetryValues[i] = DiffBinValues[i]/SumBinValues[i];
	 AsymmHist725Theta2030->SetBinContent(i+1, AsymmetryValues[i]);
     }

     //Save our asymmetry histogram
     AsymmHist725Theta2030->Draw();
     canvas49->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_725MeV_Theta2030_Asymm_7.png");
     canvas49->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_725MeV_Theta2030_Asymm_7.pdf");

     // ###############################################################################################################

     NegHelnBins = Phi_Scattered_775MeV_NegHelTheta2030->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < NegHelnBins; i++){
     	 NegHelBinValues[i] = Phi_Scattered_775MeV_NegHelTheta2030->GetBinContent(i+1);
     }

     PosHelnBins = Phi_Scattered_775MeV_PosHelTheta2030->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < PosHelnBins; i++){
     	 PosHelBinValues[i] = Phi_Scattered_775MeV_PosHelTheta2030->GetBinContent(i+1);
     }

     // Define our new histogram to plot the asymmetry here
     TCanvas *canvas50 = new TCanvas("canvas50","canvas50", 1920, 1080);
     TH1D *AsymmHist775Theta2030 = new TH1D("AsymmHist775Theta2030", "Asymmetry in PhiSc Between -ve and +ve Helicity (775pm25MeV, ThetaSc 20-30)", 90, -180, 180);

     for (Int_t i = 0; i < PosHelnBins; i++){
     	 DiffBinValues[i] = NegHelBinValues[i] - PosHelBinValues[i];
	 SumBinValues[i] = NegHelBinValues[i] + PosHelBinValues[i];
	 AsymmetryValues[i] = DiffBinValues[i]/SumBinValues[i];
	 AsymmHist775Theta2030->SetBinContent(i+1, AsymmetryValues[i]);
     }

     //Save our asymmetry histogram
     AsymmHist775Theta2030->Draw();
     canvas50->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_775MeV_Theta2030_Asymm_7.png");
     canvas50->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_775MeV_Theta2030_Asymm_7.pdf");

     // ###############################################################################################################

     NegHelnBins = Phi_Scattered_825MeV_NegHelTheta2030->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < NegHelnBins; i++){
     	 NegHelBinValues[i] = Phi_Scattered_825MeV_NegHelTheta2030->GetBinContent(i+1);
     }

     PosHelnBins = Phi_Scattered_825MeV_PosHelTheta2030->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < PosHelnBins; i++){
     	 PosHelBinValues[i] = Phi_Scattered_825MeV_PosHelTheta2030->GetBinContent(i+1);
     }

     // Define our new histogram to plot the asymmetry here
     TCanvas *canvas51 = new TCanvas("canvas51","canvas51", 1920, 1080);
     TH1D *AsymmHist825Theta2030 = new TH1D("AsymmHist825Theta2030", "Asymmetry in PhiSc Between -ve and +ve Helicity (825pm25MeV, ThetaSc 20-30)", 90, -180, 180);

     for (Int_t i = 0; i < PosHelnBins; i++){
     	 DiffBinValues[i] = NegHelBinValues[i] - PosHelBinValues[i];
	 SumBinValues[i] = NegHelBinValues[i] + PosHelBinValues[i];
	 AsymmetryValues[i] = DiffBinValues[i]/SumBinValues[i];
	 AsymmHist825Theta2030->SetBinContent(i+1, AsymmetryValues[i]);
     }

     //Save our asymmetry histogram
     AsymmHist825Theta2030->Draw();
     canvas51->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_825MeV_Theta2030_Asymm_7.png");
     canvas51->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_825MeV_Theta2030_Asymm_7.pdf");

     // ###############################################################################################################

     NegHelnBins = Phi_Scattered_875MeV_NegHelTheta2030->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < NegHelnBins; i++){
     	 NegHelBinValues[i] = Phi_Scattered_875MeV_NegHelTheta2030->GetBinContent(i+1);
     }

     PosHelnBins = Phi_Scattered_875MeV_PosHelTheta2030->GetSize() - 2; // -2 as otherwise under/overflow included
     for (Int_t i = 0; i < PosHelnBins; i++){
     	 PosHelBinValues[i] = Phi_Scattered_875MeV_PosHelTheta2030->GetBinContent(i+1);
     }

     // Define our new histogram to plot the asymmetry here
     TCanvas *canvas52 = new TCanvas("canvas52","canvas52", 1920, 1080);
     TH1D *AsymmHist875Theta2030 = new TH1D("AsymmHist875Theta2030", "Asymmetry in PhiSc Between -ve and +ve Helicity (875pm25MeV, ThetaSc 20-30)", 90, -180, 180);

     for (Int_t i = 0; i < PosHelnBins; i++){
     	 DiffBinValues[i] = NegHelBinValues[i] - PosHelBinValues[i];
	 SumBinValues[i] = NegHelBinValues[i] + PosHelBinValues[i];
	 AsymmetryValues[i] = DiffBinValues[i]/SumBinValues[i];
	 AsymmHist875Theta2030->SetBinContent(i+1, AsymmetryValues[i]);
     }

     //Save our asymmetry histogram
     AsymmHist875Theta2030->Draw();
     canvas52->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_875MeV_Theta2030_Asymm_7.png");
     canvas52->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_875MeV_Theta2030_Asymm_7.pdf");

     // ###############################################################################################################

}
