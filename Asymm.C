#include "./includes.h"

void Asymm(){

     double NegHelnBin;
     double PosHelnBins;
     double NegHelBinValues[90];
     double PosHelBinValues[90];
     double DiffBinValues[90];
     double SumBinValues[90];
     double AsymmetryValues[90];

     TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/Amo/Physics_Total_3_27_1_17.root"); // Open the latest PTotal file to load histograms from
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
     canvas->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_Asymm_3.png");
     canvas->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_Asymm_3.pdf");

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
     canvas1->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_275MeV_Asymm_3.png");
     canvas1->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_275MeV_Asymm_3.pdf");

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
     canvas2->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_325MeV_Asymm_3.png");
     canvas2->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_325MeV_Asymm_3.pdf");

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
     canvas3->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_375MeV_Asymm_3.png");
     canvas3->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_375MeV_Asymm_3.pdf");

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
     canvas4->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_425MeV_Asymm_3.png");
     canvas4->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_425MeV_Asymm_3.pdf");
    
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
     canvas5->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_475MeV_Asymm_3.png");
     canvas5->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_475MeV_Asymm_3.pdf");

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
     canvas6->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_525MeV_Asymm_3.png");
     canvas6->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_525MeV_Asymm_3.pdf");

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
     canvas7->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_575MeV_Asymm_3.png");
     canvas7->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_575MeV_Asymm_3.pdf");

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
     canvas8->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_625MeV_Asymm_3.png");
     canvas8->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_625MeV_Asymm_3.pdf");

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
     canvas9->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_675MeV_Asymm_3.png");
     canvas9->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_675MeV_Asymm_3.pdf");

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
     canvas10->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_725MeV_Asymm_3.png");
     canvas10->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_725MeV_Asymm_3.pdf");

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
     canvas11->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_775MeV_Asymm_3.png");
     canvas11->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_775MeV_Asymm_3.pdf");
     
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
     canvas12->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_825MeV_Asymm_3.png");
     canvas12->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_825MeV_Asymm_3.pdf");

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
     canvas13->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_875MeV_Asymm_3.png");
     canvas13->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_875MeV_Asymm_3.pdf");

}
