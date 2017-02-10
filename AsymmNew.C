#include "./includes.h"

void AsymmNew(){

     double NegHelnBins;
     double PosHelnBins;
     double NegHelBinValues[90];
     double PosHelBinValues[90];
     double DiffBinValues[90];
     double SumBinValues[90];
     double AsymmetryValues[90];
     double Amp[4][13]; // Format of arrays is Theta Bin Defined by i, Energy bin defined by j
     double AmpErr[4][13];
     Int_t i;
     double AmpAllTheta;
     double AmpAllThetaErr;
     double AmpTheta010;
     double AmpTheta010Err;
     double AmpTheta1020;
     double AmpTheta1020Err;
     double AmpTheta2030;
     double AmpTheta2030Err;
     
     TF1 *SinFunc = new TF1("SinFit",  "[0]*sin(x*TMath::DegToRad())", -180.0, 180.0); //Give a name and range to the fitting funcion
     SinFunc->SetParLimits(0, -1, 1);
     SinFunc->SetParNames("Amplitdue"); //Name the parameters
     SinFunc->SetParameter(0, 0);
     TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/Amo/Physics_Total_8_10_2_17.root"); // Open the latest PTotal file to load histograms from

     for(Int_t j=0; j < 13; j++){
       
       i=0;
       
       if (j==0){
	         Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (275pm25MeV)"; // Set title of output graph
		 TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_275MeV_NegHel"); // Select the correct histogram
		 TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_275MeV_PosHel"); // Select the correct histogram
		 Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_275MeV_Asymm_8.png"; // Name the output images
		 Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_275MeV_Asymm_8.pdf";
       }

       if (j==1){
	         Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (325pm25MeV)"; // Set title of output graph
		 TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_325MeV_NegHel"); // Select the correct histogram
		 TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_325MeV_PosHel"); // Select the correct histogram
		 Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_325MeV_Asymm_8.png"; // Name the output images
		 Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_325MeV_Asymm_8.pdf";
       }
       
       if (j==2){
	         Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (375pm25MeV)"; // Set title of output graph
		 TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_375MeV_NegHel"); // Select the correct histogram
		 TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_375MeV_PosHel"); // Select the correct histogram
		 Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_375MeV_Asymm_8.png"; // Name the output images
		 Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_375MeV_Asymm_8.pdf";
       }

       if (j==3){
	         Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (425pm25MeV)"; // Set title of output graph
		 TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_425MeV_PosHel"); // Select the correct histogram
		 TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_425MeV_PosHel"); // Select the correct histogram
		 Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_425MeV_Asymm_8.png"; // Name the output images
		 Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_425MeV_Asymm_8.pdf";
       }

       if (j==4){
	         Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (475pm25MeV)"; // Set title of output graph
		 TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_475MeV_NegHel"); // Select the correct histogram
		 TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_475MeV_PosHel"); // Select the correct histogram
		 Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_475MeV_Asymm_8.png"; // Name the output images
		 Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_475MeV_Asymm_8.pdf";
       }

       if (j==5){
	         Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (525pm25MeV)"; // Set title of output graph
		 TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_525MeV_NegHel"); // Select the correct histogram
		 TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_525MeV_PosHel"); // Select the correct histogram
		 Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_525MeV_Asymm_8.png"; // Name the output images
		 Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_525MeV_Asymm_8.pdf";
       }

       if (j==6){
	         Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (575pm25MeV)"; // Set title of output graph
		 TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_575MeV_NegHel"); // Select the correct histogram
		 TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_575MeV_PosHel"); // Select the correct histogram
		 Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_575MeV_Asymm_8.png"; // Name the output images
		 Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_575MeV_Asymm_8.pdf";
       }

       if (j==7){
	         Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (625pm25MeV)"; // Set title of output graph
		 TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_625MeV_NegHel"); // Select the correct histogram
		 TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_625MeV_PosHel"); // Select the correct histogram
		 Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_625MeV_Asymm_8.png"; // Name the output images
		 Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_625MeV_Asymm_8.pdf";
       }

       if (j==8){
	         Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (675pm25MeV)"; // Set title of output graph
		 TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_675MeV_NegHel"); // Select the correct histogram
		 TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_675MeV_PosHel"); // Select the correct histogram
		 Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_675MeV_Asymm_8.png"; // Name the output images
		 Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_675MeV_Asymm_8.pdf";
       }

       if (j==9){
	         Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (725pm25MeV)"; // Set title of output graph
		 TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_725MeV_NegHel"); // Select the correct histogram
		 TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_725MeV_PosHel"); // Select the correct histogram
		 Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_725MeV_Asymm_8.png"; // Name the output images
		 Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_725MeV_Asymm_8.pdf";
       }

       if (j==10){
	         Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (775pm25MeV)"; // Set title of output graph
		 TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_775MeV_NegHel"); // Select the correct histogram
		 TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_775MeV_PosHel"); // Select the correct histogram
		 Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_775MeV_Asymm_8.png"; // Name the output images
		 Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_775MeV_Asymm_8.pdf";
       }

      if (j==11){
	         Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (825pm25MeV)"; // Set title of output graph
		 TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_825MeV_NegHel"); // Select the correct histogram
		 TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_825MeV_PosHel"); // Select the correct histogram
		 Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_825MeV_Asymm_8.png"; // Name the output images
		 Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_825MeV_Asymm_8.pdf";
       }

      if (j==12){
	         Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (875pm25MeV)"; // Set title of output graph
		 TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_875MeV_NegHel"); // Select the correct histogram
		 TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_875MeV_PosHel"); // Select the correct histogram
		 Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_875MeV_Asymm_8.png"; // Name the output images
		 Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_875MeV_Asymm_8.pdf";
       }

       TCanvas *canvas = new TCanvas("canvas","canvas", 1920, 1080);
       TH1D *AsymmHist = new TH1D("AsymmHist", "", 90, -180, 180);

       NegHelnBins = histNeg->GetSize() - 2; // -2 as otherwise under/overflow included
       for (Int_t k = 0; k < NegHelnBins; k++){
	 NegHelBinValues[k] = histNeg->GetBinContent(k+1);
       }

       PosHelnBins = histPos->GetSize() - 2; // -2 as otherwise under/overflow included
       for (Int_t k = 0; k < PosHelnBins; k++){
	 PosHelBinValues[k] = histPos->GetBinContent(k+1);
       }

       for (Int_t k = 0; k < PosHelnBins; k++){
     	 DiffBinValues[k] = NegHelBinValues[k] - PosHelBinValues[k];
	 SumBinValues[k] = NegHelBinValues[k] + PosHelBinValues[k];
	 AsymmetryValues[k] = DiffBinValues[k]/SumBinValues[k];
	 AsymmHist->SetBinContent(k+1, AsymmetryValues[k]);
       }

       TPad *pad1 = new TPad("pad1","",0,0,1,1);
       pad1->Draw();
       pad1->cd();

       pad1->SetTickx(1);
       pad1->SetTicky(1);
       pad1->SetGridx(1);
       pad1->SetGridy(1);

       TH1F  *hr;
       Char_t hrTitle[64];
       
       Float_t xMin = -180;
       Float_t xMax = 180;
       Float_t yMin = -2;
       Float_t yMax = 2;
       
       strcpy(hrTitle, Title);
       hr = canvas->DrawFrame(xMin,yMin,xMax,yMax);
       hr->SetTitle(hrTitle);
       
       AsymmHist->SetMarkerStyle(1); // Style options for graph
       AsymmHist->SetLineColor(2);
       AsymmHist->Rebin(2);

       AsymmHist->Draw("HISTSAMES"); // Draw the histogram
       AsymmHist->Fit("SinFit"); // Fit sine function to histogram
       SinFit->SetLineColor(4);
       SinFit->Draw("SAMES"); // Draw the resulting fit
       gStyle->SetOptFit(0111);
       gPad->Update(); // Refresh plot

       Amp[i][j] = SinFit->GetParameter(0); // Add values of the fit to an array
       AmpErr[i][j] = SinFit->GetParError(0);

       canvas->SaveAs(filename = GraphPDF);
       canvas->SaveAs(filename = GraphPNG);
       
     }

     // #######################################################################################################################
     // #######################################################################################################################
     // 0-10 Degree Plots
     // #######################################################################################################################
     // #######################################################################################################################

     for(Int_t j=0; j < 13; j++){
       
       i=1;
       
       if (j==0){
	         Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (275pm25MeV, ThetaSc 0-10)"; // Set title of output graph
		 TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_275MeV_NegHelTheta010"); // Select the correct histogram
		 TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_275MeV_PosHelTheta010"); // Select the correct histogram
		 Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_275MeV_Theta010_Asymm_8.png"; 
		 Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_275MeV_Theta010_Asymm_8.pdf";
       }

       if (j==1){
	         Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (325pm25MeV, ThetaSc 0-10)"; // Set title of output graph
		 TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_325MeV_NegHelTheta010"); // Select the correct histogram
		 TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_325MeV_PosHelTheta010"); // Select the correct histogram
		 Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_325MeV_Theta010_Asymm_8.png"; 
		 Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_325MeV_Theta010_Asymm_8.pdf";
       }
       
       if (j==2){
	         Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (375pm25MeV, ThetaSc 0-10)"; // Set title of output graph
		 TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_375MeV_NegHelTheta010"); // Select the correct histogram
		 TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_375MeV_PosHelTheta010"); // Select the correct histogram
		 Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_375MeV_Theta010_Asymm_8.png"; 
		 Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_375MeV_Theta010_Asymm_8.pdf";
       }

       if (j==3){
	         Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (425pm25MeV, ThetaSc 0-10)"; // Set title of output graph
		 TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_425MeV_PosHelTheta010"); // Select the correct histogram
		 TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_425MeV_PosHelTheta010"); // Select the correct histogram
		 Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_425MeV_Theta010_Asymm_8.png"; 
		 Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_425MeV_Theta010_Asymm_8.pdf";
       }

       if (j==4){
	         Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (475pm25MeV, ThetaSc 0-10)"; // Set title of output graph
		 TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_475MeV_NegHelTheta010"); // Select the correct histogram
		 TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_475MeV_PosHelTheta010"); // Select the correct histogram
		 Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_475MeV_Theta010_Asymm_8.png"; 
		 Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_475MeV_Theta010_Asymm_8.pdf";
       }

       if (j==5){
	         Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (525pm25MeV, ThetaSc 0-10)"; // Set title of output graph
		 TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_525MeV_NegHelTheta010"); // Select the correct histogram
		 TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_525MeV_PosHelTheta010"); // Select the correct histogram
		 Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_525MeV_Theta010_Asymm_8.png"; 
		 Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_525MeV_Theta010_Asymm_8.pdf";
       }

       if (j==6){
	         Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (575pm25MeV, ThetaSc 0-10)"; // Set title of output graph
		 TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_575MeV_NegHelTheta010"); // Select the correct histogram
		 TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_575MeV_PosHelTheta010"); // Select the correct histogram
		 Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_575MeV_Theta010_Asymm_8.png"; 
		 Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_575MeV_Theta010_Asymm_8.pdf";
       }

       if (j==7){
	         Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (625pm25MeV, ThetaSc 0-10)"; // Set title of output graph
		 TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_625MeV_NegHelTheta010"); // Select the correct histogram
		 TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_625MeV_PosHelTheta010"); // Select the correct histogram
		 Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_625MeV_Theta010_Asymm_8.png"; 
		 Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_625MeV_Theta010_Asymm_8.pdf";
       }

       if (j==8){
	         Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (675pm25MeV, ThetaSc 0-10)"; // Set title of output graph
		 TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_675MeV_NegHelTheta010"); // Select the correct histogram
		 TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_675MeV_PosHelTheta010"); // Select the correct histogram
		 Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_675MeV_Theta010_Asymm_8.png"; 
		 Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_675MeV_Theta010_Asymm_8.pdf";
       }

       if (j==9){
	         Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (725pm25MeV, ThetaSc 0-10)"; // Set title of output graph
		 TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_725MeV_NegHelTheta010"); // Select the correct histogram
		 TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_725MeV_PosHelTheta010"); // Select the correct histogram
		 Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_725MeV_Theta010_Asymm_8.png"; 
		 Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_725MeV_Theta010_Asymm_8.pdf";
       }

       if (j==10){
	         Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (775pm25MeV, ThetaSc 0-10)"; // Set title of output graph
		 TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_775MeV_NegHelTheta010"); // Select the correct histogram
		 TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_775MeV_PosHelTheta010"); // Select the correct histogram
		 Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_775MeV_Theta010_Asymm_8.png"; 
		 Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_775MeV_Theta010_Asymm_8.pdf";
       }

      if (j==11){
	         Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (825pm25MeV, ThetaSc 0-10)"; // Set title of output graph
		 TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_825MeV_NegHelTheta010"); // Select the correct histogram
		 TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_825MeV_PosHelTheta010"); // Select the correct histogram
		 Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_825MeV_Theta010_Asymm_8.png"; 
		 Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_825MeV_Theta010_Asymm_8.pdf";
       }

      if (j==12){
	         Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (875pm25MeV, ThetaSc 0-10)"; // Set title of output graph
		 TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_875MeV_NegHelTheta010"); // Select the correct histogram
		 TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_875MeV_PosHelTheta010"); // Select the correct histogram
		 Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_875MeV_Theta010_Asymm_8.png";
		 Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_875MeV_Theta010_Asymm_8.pdf";
       }

       TCanvas *canvas = new TCanvas("canvas","canvas", 1920, 1080);
       TH1D *AsymmHist = new TH1D("AsymmHist", "", 90, -180, 180);

       NegHelnBins = histNeg->GetSize() - 2; // -2 as otherwise under/overflow included
       for (Int_t k = 0; k < NegHelnBins; k++){
	 NegHelBinValues[k] = histNeg->GetBinContent(k+1);
       }

       PosHelnBins = histPos->GetSize() - 2; // -2 as otherwise under/overflow included
       for (Int_t k = 0; k < PosHelnBins; k++){
	 PosHelBinValues[k] = histPos->GetBinContent(k+1);
       }

       for (Int_t k = 0; k < PosHelnBins; k++){
     	 DiffBinValues[k] = NegHelBinValues[k] - PosHelBinValues[k];
	 SumBinValues[k] = NegHelBinValues[k] + PosHelBinValues[k];
	 AsymmetryValues[k] = DiffBinValues[k]/SumBinValues[k];
	 AsymmHist->SetBinContent(k+1, AsymmetryValues[k]);
       }

       TPad *pad1 = new TPad("pad1","",0,0,1,1);
       pad1->Draw();
       pad1->cd();

       pad1->SetTickx(1);
       pad1->SetTicky(1);
       pad1->SetGridx(1);
       pad1->SetGridy(1);

       TH1F  *hr;
       Char_t hrTitle[64];
       
       Float_t xMin = -180;
       Float_t xMax = 180;
       Float_t yMin = -2;
       Float_t yMax = 2;
       
       strcpy(hrTitle, Title);
       hr = canvas->DrawFrame(xMin,yMin,xMax,yMax);
       hr->SetTitle(hrTitle);
       
       AsymmHist->SetMarkerStyle(1); // Style options for graph
       AsymmHist->SetLineColor(2);

       AsymmHist->Draw("HISTSAMES"); // Draw the histogram
       AsymmHist->Fit("SinFit"); // Fit sine function to histogram
       SinFit->SetLineColor(4);
       SinFit->Draw("SAMES"); // Draw the resulting fit
       gStyle->SetOptFit(0111);
       gPad->Update(); // Refresh plot

       Amp[i][j] = SinFit->GetParameter(0); // Add values of the fit to an array
       AmpErr[i][j] = SinFit->GetParError(0);

       canvas->SaveAs(filename = GraphPDF);
       canvas->SaveAs(filename = GraphPNG);
       
     }

     // #######################################################################################################################
     // #######################################################################################################################
     // 10-20 Degree Plots
     // #######################################################################################################################
     // #######################################################################################################################

     for(Int_t j=0; j < 13; j++){
       
       i=2;
       
       if (j==0){
	         Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (275pm25MeV, ThetaSc 10-20)"; // Set title of output graph
		 TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_275MeV_NegHelTheta1020"); // Select the correct histogram
		 TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_275MeV_PosHelTheta1020"); // Select the correct histogram
		 Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_275MeV_Theta1020_Asymm_8.png"; 
		 Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_275MeV_Theta1020_Asymm_8.pdf";
       }

       if (j==1){
	         Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (325pm25MeV, ThetaSc 10-20)"; // Set title of output graph
		 TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_325MeV_NegHelTheta1020"); // Select the correct histogram
		 TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_325MeV_PosHelTheta1020"); // Select the correct histogram
		 Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_325MeV_Theta1020_Asymm_8.png"; 
		 Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_325MeV_Theta1020_Asymm_8.pdf";
       }
       
       if (j==2){
	         Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (375pm25MeV, ThetaSc 10-20)"; // Set title of output graph
		 TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_375MeV_NegHelTheta1020"); // Select the correct histogram
		 TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_375MeV_PosHelTheta1020"); // Select the correct histogram
		 Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_375MeV_Theta1020_Asymm_8.png"; 
		 Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_375MeV_Theta1020_Asymm_8.pdf";
       }

       if (j==3){
	         Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (425pm25MeV, ThetaSc 10-20)"; // Set title of output graph
		 TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_425MeV_PosHelTheta1020"); // Select the correct histogram
		 TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_425MeV_PosHelTheta1020"); // Select the correct histogram
		 Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_425MeV_Theta1020_Asymm_8.png"; 
		 Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_425MeV_Theta1020_Asymm_8.pdf";
       }

       if (j==4){
	         Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (475pm25MeV, ThetaSc 10-20)"; // Set title of output graph
		 TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_475MeV_NegHelTheta1020"); // Select the correct histogram
		 TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_475MeV_PosHelTheta1020"); // Select the correct histogram
		 Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_475MeV_Theta1020_Asymm_8.png"; 
		 Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_475MeV_Theta1020_Asymm_8.pdf";
       }

       if (j==5){
	         Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (525pm25MeV, ThetaSc 10-20)"; // Set title of output graph
		 TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_525MeV_NegHelTheta1020"); // Select the correct histogram
		 TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_525MeV_PosHelTheta1020"); // Select the correct histogram
		 Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_525MeV_Theta1020_Asymm_8.png"; 
		 Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_525MeV_Theta1020_Asymm_8.pdf";
       }

       if (j==6){
	         Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (575pm25MeV, ThetaSc 10-20)"; // Set title of output graph
		 TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_575MeV_NegHelTheta1020"); // Select the correct histogram
		 TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_575MeV_PosHelTheta1020"); // Select the correct histogram
		 Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_575MeV_Theta1020_Asymm_8.png"; 
		 Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_575MeV_Theta1020_Asymm_8.pdf";
       }

       if (j==7){
	         Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (625pm25MeV, ThetaSc 10-20)"; // Set title of output graph
		 TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_625MeV_NegHelTheta1020"); // Select the correct histogram
		 TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_625MeV_PosHelTheta1020"); // Select the correct histogram
		 Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_625MeV_Theta1020_Asymm_8.png"; 
		 Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_625MeV_Theta1020_Asymm_8.pdf";
       }

       if (j==8){
	         Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (675pm25MeV, ThetaSc 10-20)"; // Set title of output graph
		 TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_675MeV_NegHelTheta1020"); // Select the correct histogram
		 TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_675MeV_PosHelTheta1020"); // Select the correct histogram
		 Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_675MeV_Theta1020_Asymm_8.png"; 
		 Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_675MeV_Theta1020_Asymm_8.pdf";
       }

       if (j==9){
	         Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (725pm25MeV, ThetaSc 10-20)"; // Set title of output graph
		 TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_725MeV_NegHelTheta1020"); // Select the correct histogram
		 TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_725MeV_PosHelTheta1020"); // Select the correct histogram
		 Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_725MeV_Theta1020_Asymm_8.png"; 
		 Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_725MeV_Theta1020_Asymm_8.pdf";
       }

       if (j==10){
	         Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (775pm25MeV, ThetaSc 10-20)"; // Set title of output graph
		 TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_775MeV_NegHelTheta1020"); // Select the correct histogram
		 TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_775MeV_PosHelTheta1020"); // Select the correct histogram
		 Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_775MeV_Theta1020_Asymm_8.png"; 
		 Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_775MeV_Theta1020_Asymm_8.pdf";
       }

      if (j==11){
	         Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (825pm25MeV, ThetaSc 10-20)"; // Set title of output graph
		 TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_825MeV_NegHelTheta1020"); // Select the correct histogram
		 TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_825MeV_PosHelTheta1020"); // Select the correct histogram
		 Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_825MeV_Theta1020_Asymm_8.png"; 
		 Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_825MeV_Theta1020_Asymm_8.pdf";
       }

      if (j==12){
	         Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (875pm25MeV, ThetaSc 10-20)"; // Set title of output graph
		 TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_875MeV_NegHelTheta1020"); // Select the correct histogram
		 TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_875MeV_PosHelTheta1020"); // Select the correct histogram
		 Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_875MeV_Theta1020_Asymm_8.png";
		 Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_875MeV_Theta1020_Asymm_8.pdf";
       }

       TCanvas *canvas = new TCanvas("canvas","canvas", 1920, 1080);
       TH1D *AsymmHist = new TH1D("AsymmHist", "", 90, -180, 180);

       NegHelnBins = histNeg->GetSize() - 2; // -2 as otherwise under/overflow included
       for (Int_t k = 0; k < NegHelnBins; k++){
	 NegHelBinValues[k] = histNeg->GetBinContent(k+1);
       }

       PosHelnBins = histPos->GetSize() - 2; // -2 as otherwise under/overflow included
       for (Int_t k = 0; k < PosHelnBins; k++){
	 PosHelBinValues[k] = histPos->GetBinContent(k+1);
       }

       for (Int_t k = 0; k < PosHelnBins; k++){
     	 DiffBinValues[k] = NegHelBinValues[k] - PosHelBinValues[k];
	 SumBinValues[k] = NegHelBinValues[k] + PosHelBinValues[k];
	 AsymmetryValues[k] = DiffBinValues[k]/SumBinValues[k];
	 AsymmHist->SetBinContent(k+1, AsymmetryValues[k]);
       }

       TPad *pad1 = new TPad("pad1","",0,0,1,1);
       pad1->Draw();
       pad1->cd();

       pad1->SetTickx(1);
       pad1->SetTicky(1);
       pad1->SetGridx(1);
       pad1->SetGridy(1);

       TH1F  *hr;
       Char_t hrTitle[64];
       
       Float_t xMin = -180;
       Float_t xMax = 180;
       Float_t yMin = -2;
       Float_t yMax = 2;
       
       strcpy(hrTitle, Title);
       hr = canvas->DrawFrame(xMin,yMin,xMax,yMax);
       hr->SetTitle(hrTitle);
       
       AsymmHist->SetMarkerStyle(1); // Style options for graph
       AsymmHist->SetLineColor(2);

       AsymmHist->Draw("HISTSAMES"); // Draw the histogram
       AsymmHist->Fit("SinFit"); // Fit sine function to histogram
       SinFit->SetLineColor(4);
       SinFit->Draw("SAMES"); // Draw the resulting fit
       gStyle->SetOptFit(0111);
       gPad->Update(); // Refresh plot

       Amp[i][j] = SinFit->GetParameter(0); // Add values of the fit to an array
       AmpErr[i][j] = SinFit->GetParError(0);

       canvas->SaveAs(filename = GraphPDF);
       canvas->SaveAs(filename = GraphPNG);
       
     }

     // #######################################################################################################################
     // #######################################################################################################################
     // 20-30 Degree Plots
     // #######################################################################################################################
     // #######################################################################################################################

     for(Int_t j=0; j < 13; j++){
       
       i=3;
       
       if (j==0){
	         Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (275pm25MeV, ThetaSc 20-30)"; // Set title of output graph
		 TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_275MeV_NegHelTheta2030"); // Select the correct histogram
		 TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_275MeV_PosHelTheta2030"); // Select the correct histogram
		 Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_275MeV_Theta2030_Asymm_8.png"; 
		 Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_275MeV_Theta2030_Asymm_8.pdf";
       }

       if (j==1){
	         Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (325pm25MeV, ThetaSc 20-30)"; // Set title of output graph
		 TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_325MeV_NegHelTheta2030"); // Select the correct histogram
		 TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_325MeV_PosHelTheta2030"); // Select the correct histogram
		 Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_325MeV_Theta2030_Asymm_8.png"; 
		 Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_325MeV_Theta2030_Asymm_8.pdf";
       }
       
       if (j==2){
	         Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (375pm25MeV, ThetaSc 20-30)"; // Set title of output graph
		 TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_375MeV_NegHelTheta2030"); // Select the correct histogram
		 TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_375MeV_PosHelTheta2030"); // Select the correct histogram
		 Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_375MeV_Theta2030_Asymm_8.png"; 
		 Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_375MeV_Theta2030_Asymm_8.pdf";
       }

       if (j==3){
	         Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (425pm25MeV, ThetaSc 20-30)"; // Set title of output graph
		 TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_425MeV_PosHelTheta2030"); // Select the correct histogram
		 TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_425MeV_PosHelTheta2030"); // Select the correct histogram
		 Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_425MeV_Theta2030_Asymm_8.png"; 
		 Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_425MeV_Theta2030_Asymm_8.pdf";
       }

       if (j==4){
	         Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (475pm25MeV, ThetaSc 20-30)"; // Set title of output graph
		 TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_475MeV_NegHelTheta2030"); // Select the correct histogram
		 TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_475MeV_PosHelTheta2030"); // Select the correct histogram
		 Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_475MeV_Theta2030_Asymm_8.png"; 
		 Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_475MeV_Theta2030_Asymm_8.pdf";
       }

       if (j==5){
	         Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (525pm25MeV, ThetaSc 20-30)"; // Set title of output graph
		 TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_525MeV_NegHelTheta2030"); // Select the correct histogram
		 TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_525MeV_PosHelTheta2030"); // Select the correct histogram
		 Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_525MeV_Theta2030_Asymm_8.png"; 
		 Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_525MeV_Theta2030_Asymm_8.pdf";
       }

       if (j==6){
	         Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (575pm25MeV, ThetaSc 20-30)"; // Set title of output graph
		 TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_575MeV_NegHelTheta2030"); // Select the correct histogram
		 TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_575MeV_PosHelTheta2030"); // Select the correct histogram
		 Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_575MeV_Theta2030_Asymm_8.png"; 
		 Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_575MeV_Theta2030_Asymm_8.pdf";
       }

       if (j==7){
	         Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (625pm25MeV, ThetaSc 20-30)"; // Set title of output graph
		 TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_625MeV_NegHelTheta2030"); // Select the correct histogram
		 TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_625MeV_PosHelTheta2030"); // Select the correct histogram
		 Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_625MeV_Theta2030_Asymm_8.png"; 
		 Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_625MeV_Theta2030_Asymm_8.pdf";
       }

       if (j==8){
	         Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (675pm25MeV, ThetaSc 20-30)"; // Set title of output graph
		 TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_675MeV_NegHelTheta2030"); // Select the correct histogram
		 TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_675MeV_PosHelTheta2030"); // Select the correct histogram
		 Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_675MeV_Theta2030_Asymm_8.png"; 
		 Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_675MeV_Theta2030_Asymm_8.pdf";
       }

       if (j==9){
	         Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (725pm25MeV, ThetaSc 20-30)"; // Set title of output graph
		 TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_725MeV_NegHelTheta2030"); // Select the correct histogram
		 TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_725MeV_PosHelTheta2030"); // Select the correct histogram
		 Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_725MeV_Theta2030_Asymm_8.png"; 
		 Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_725MeV_Theta2030_Asymm_8.pdf";
       }

       if (j==10){
	         Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (775pm25MeV, ThetaSc 20-30)"; // Set title of output graph
		 TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_775MeV_NegHelTheta2030"); // Select the correct histogram
		 TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_775MeV_PosHelTheta2030"); // Select the correct histogram
		 Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_775MeV_Theta2030_Asymm_8.png"; 
		 Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_775MeV_Theta2030_Asymm_8.pdf";
       }

      if (j==11){
	         Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (825pm25MeV, ThetaSc 20-30)"; // Set title of output graph
		 TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_825MeV_NegHelTheta2030"); // Select the correct histogram
		 TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_825MeV_PosHelTheta2030"); // Select the correct histogram
		 Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_825MeV_Theta2030_Asymm_8.png"; 
		 Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_825MeV_Theta2030_Asymm_8.pdf";
       }

      if (j==12){
	         Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (875pm25MeV, ThetaSc 20-30)"; // Set title of output graph
		 TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_875MeV_NegHelTheta2030"); // Select the correct histogram
		 TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_875MeV_PosHelTheta2030"); // Select the correct histogram
		 Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_875MeV_Theta2030_Asymm_8.png";
		 Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_875MeV_Theta2030_Asymm_8.pdf";
       }

       TCanvas *canvas = new TCanvas("canvas","canvas", 1920, 1080);
       TH1D *AsymmHist = new TH1D("AsymmHist", "", 90, -180, 180);

       NegHelnBins = histNeg->GetSize() - 2; // -2 as otherwise under/overflow included
       for (Int_t k = 0; k < NegHelnBins; k++){
	 NegHelBinValues[k] = histNeg->GetBinContent(k+1);
       }

       PosHelnBins = histPos->GetSize() - 2; // -2 as otherwise under/overflow included
       for (Int_t k = 0; k < PosHelnBins; k++){
	 PosHelBinValues[k] = histPos->GetBinContent(k+1);
       }

       for (Int_t k = 0; k < PosHelnBins; k++){
     	 DiffBinValues[k] = NegHelBinValues[k] - PosHelBinValues[k];
	 SumBinValues[k] = NegHelBinValues[k] + PosHelBinValues[k];
	 AsymmetryValues[k] = DiffBinValues[k]/SumBinValues[k];
	 AsymmHist->SetBinContent(k+1, AsymmetryValues[k]);
       }

       TPad *pad1 = new TPad("pad1","",0,0,1,1);
       pad1->Draw();
       pad1->cd();

       pad1->SetTickx(1);
       pad1->SetTicky(1);
       pad1->SetGridx(1);
       pad1->SetGridy(1);

       TH1F  *hr;
       Char_t hrTitle[64];
       
       Float_t xMin = -180;
       Float_t xMax = 180;
       Float_t yMin = -2;
       Float_t yMax = 2;
       
       strcpy(hrTitle, Title);
       hr = canvas->DrawFrame(xMin,yMin,xMax,yMax);
       hr->SetTitle(hrTitle);
       
       AsymmHist->SetMarkerStyle(1); // Style options for graph
       AsymmHist->SetLineColor(2);

       AsymmHist->Draw("HISTSAMES"); // Draw the histogram
       AsymmHist->Fit("SinFit", "LL"); // Fit sine function to histogram
       SinFit->SetLineColor(4);
       SinFit->Draw("SAMES"); // Draw the resulting fit
       gStyle->SetOptFit(0111);
       gPad->Update(); // Refresh plot

       Amp[i][j] = SinFit->GetParameter(0); // Add values of the fit to an array
       AmpErr[i][j] = SinFit->GetParError(0);

       canvas->SaveAs(filename = GraphPDF);
       canvas->SaveAs(filename = GraphPNG);
       
     }

    // Define new file to store fit parameters
    TFile f1("SinPhiFitValues.root", "RECREATE");

    //Define new tree to store parameters in
    TTree* tree = new TTree("Parameter_Values", "Tree_of_Values");

    // Define branches to store parameters, (Branch Name, Variable, Type of Variable)
    tree->Branch("AmpAllTheta", &AmpAllTheta, "AmpAllTheta/D");
    tree->Branch("AmpAllThetaErr", &AmpAllThetaErr, "AmpAllThetaErr/D");
    tree->Branch("AmpTheta010", &AmpTheta010, "AmpTheta010/D");
    tree->Branch("AmpTheta010Err", &AmpTheta010Err, "AmpTheta010Err/D");
    tree->Branch("AmpTheta1020", &AmpTheta1020, "AmpTheta1020/D");
    tree->Branch("AmpTheta1020Err", &AmpTheta1020Err, "AmpTheta1020Err/D");
    tree->Branch("AmpTheta2030", &AmpTheta2030, "AmpTheta2030/D");
    tree->Branch("AmpTheta2030Err", &AmpTheta2030Err, "AmpTheta2030Err/D");

    // Fill branches (and hence tree) with corresponding parameters from above
    for (Int_t m = 0; m < 13; m++){
      
      AmpAllTheta = Amp[0][m];
      AmpAllThetaErr = AmpErr[0][m];
      AmpTheta010 = Amp[1][m];
      AmpTheta010Err = AmpErr[1][m];
      AmpTheta1020 = Amp[2][m];
      AmpTheta1020Err = AmpErr[2][m];
      AmpTheta2030 = Amp[3][m];
      AmpTheta2030Err = AmpErr[3][m];
      tree->Fill();

    }

    f1.Write();

}
