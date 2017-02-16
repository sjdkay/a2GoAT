#include "./includes.h"

// define a function with 4 parameters
Double_t fitf(Double_t *x,Double_t *par)
{
    Double_t fitval = 0;
    fitval =  par[0] + ((par[1]*sin(x[0]*TMath::DegToRad()))/(1 + (par[2]*cos(x[0]*TMath::DegToRad()))));
    return fitval;
}

void AsymmNew(){

     double Offset[4][13]; // Format of arrays is Theta Bin Defined by i, Energy bin defined by j
     double OffsetErr[4][13];
     double SinAmp[4][13]; // Format of arrays is Theta Bin Defined by i, Energy bin defined by j
     double SinAmpErr[4][13];
     double CosAmp[4][13]; // Format of arrays is Theta Bin Defined by i, Energy bin defined by j
     double CosAmpErr[4][13];
     Int_t i;
     double OffsetAllTheta;
     double OffsetAllThetaErr;
     double OffsetTheta010;
     double OffsetTheta010Err;
     double OffsetTheta1020;
     double OffsetTheta1020Err;
     double OffsetTheta2030;
     double OffsetTheta2030Err;
     double SinAmpAllTheta;
     double SinAmpAllThetaErr;
     double SinAmpTheta010;
     double SinAmpTheta010Err;
     double SinAmpTheta1020;
     double SinAmpTheta1020Err;
     double SinAmpTheta2030;
     double SinAmpTheta2030Err;
     double CosAmpAllTheta;
     double CosAmpAllThetaErr;
     double CosAmpTheta010;
     double CosAmpTheta010Err;
     double CosAmpTheta1020;
     double CosAmpTheta1020Err;
     double CosAmpTheta2030;
     double CosAmpTheta2030Err;
     
     TF1 *AsymmFunc = new TF1("AsymmFit",  fitf, -180.0, 180.0, 3); //Give a name and range to the fitting funcion
     AsymmFunc->SetParNames("Offset", "SinAmp", "CosAmp"); //Name the parameters
     AsymmFunc->SetParameter(0, 0);
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

       AsymmHist = histNeg->GetAsymmetry(histPos);

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
       AsymmHist->Fit("AsymmFit", "LL"); // Fit sine function to histogram
       AsymmFit->SetLineColor(4);
       AsymmFit->Draw("SAMES"); // Draw the resulting fit
       gStyle->SetOptFit(11);
       gStyle->SetStatW(0.1);
       gStyle->SetStatH(0.1);
       gPad->Update(); // Refresh plot

       Offset[i][j] = AsymmFit->GetParameter(0); // Add values of the fit to an array
       OffsetErr[i][j] = AsymmFit->GetParError(0);
       SinAmp[i][j] = AsymmFit->GetParameter(1);
       SinAmpErr[i][j] = AsymmFit->GetParError(1);
       CosAmp[i][j] = AsymmFit->GetParameter(2);
       CosAmpErr[i][j] = AsymmFit->GetParError(2);

       canvas->SaveAs(filename = GraphPDF);
       canvas->SaveAs(filename = GraphPNG);
       canvas->Close();
       
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

       AsymmHist = histNeg->GetAsymmetry(histPos);

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
       AsymmHist->Fit("AsymmFit", "LL"); // Fit sine function to histogram
       AsymmFit->SetLineColor(4);
       AsymmFit->Draw("SAMES"); // Draw the resulting fit
       gStyle->SetOptFit(11);
       gStyle->SetStatW(0.1);
       gStyle->SetStatH(0.1);
       gPad->Update(); // Refresh plot

       Offset[i][j] = AsymmFit->GetParameter(0); // Add values of the fit to an array
       OffsetErr[i][j] = AsymmFit->GetParError(0);
       SinAmp[i][j] = AsymmFit->GetParameter(1);
       SinAmpErr[i][j] = AsymmFit->GetParError(1);
       CosAmp[i][j] = AsymmFit->GetParameter(2);
       CosAmpErr[i][j] = AsymmFit->GetParError(2);

       canvas->SaveAs(filename = GraphPDF);
       canvas->SaveAs(filename = GraphPNG);
       canvas->Close();
       
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

       AsymmHist = histNeg->GetAsymmetry(histPos);

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
       AsymmHist->Fit("AsymmFit", "LL"); // Fit sine function to histogram
       AsymmFit->SetLineColor(4);
       AsymmFit->Draw("SAMES"); // Draw the resulting fit
       gStyle->SetOptFit(11);
       gStyle->SetStatW(0.1);
       gStyle->SetStatH(0.1);
       gPad->Update(); // Refresh plot

       Offset[i][j] = AsymmFit->GetParameter(0); // Add values of the fit to an array
       OffsetErr[i][j] = AsymmFit->GetParError(0);
       SinAmp[i][j] = AsymmFit->GetParameter(1);
       SinAmpErr[i][j] = AsymmFit->GetParError(1);
       CosAmp[i][j] = AsymmFit->GetParameter(2);
       CosAmpErr[i][j] = AsymmFit->GetParError(2);

       canvas->SaveAs(filename = GraphPDF);
       canvas->SaveAs(filename = GraphPNG);
       canvas->Close();
       
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

       AsymmHist = histNeg->GetAsymmetry(histPos);

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
       AsymmHist->Fit("AsymmFit", "LL"); // Fit sine function to histogram
       AsymmFit->SetLineColor(4);
       AsymmFit->Draw("SAMES"); // Draw the resulting fit
       gStyle->SetOptFit(11);
       gStyle->SetStatW(0.1);
       gStyle->SetStatH(0.1);
       gPad->Update(); // Refresh plot

       Offset[i][j] = AsymmFit->GetParameter(0); // Add values of the fit to an array
       OffsetErr[i][j] = AsymmFit->GetParError(0);
       SinAmp[i][j] = AsymmFit->GetParameter(1);
       SinAmpErr[i][j] = AsymmFit->GetParError(1);
       CosAmp[i][j] = AsymmFit->GetParameter(2);
       CosAmpErr[i][j] = AsymmFit->GetParError(2);

       canvas->SaveAs(filename = GraphPDF);
       canvas->SaveAs(filename = GraphPNG);
       canvas->Close();
       
     }

    // Define new file to store fit parameters
    TFile f1("SinPhiFitValues.root", "RECREATE");

    //Define new tree to store parameters in
    TTree* tree = new TTree("Parameter_Values", "Tree_of_Values");

    // Define branches to store parameters, (Branch Name, Variable, Type of Variable)
    tree->Branch("OffsetAllTheta", &OffsetAllTheta, "OffsetAllTheta/D");
    tree->Branch("OffsetAllThetaErr", &OffsetAllThetaErr, "OffsetAllThetaErr/D");
    tree->Branch("OffsetTheta010", &OffsetTheta010, "OffsetTheta010/D");
    tree->Branch("OffsetTheta010Err", &OffsetTheta010Err, "OffsetTheta010Err/D");
    tree->Branch("OffsetTheta1020", &OffsetTheta1020, "OffsetTheta1020/D");
    tree->Branch("OffsetTheta1020Err", &OffsetTheta1020Err, "OffsetTheta1020Err/D");
    tree->Branch("OffsetTheta2030", &OffsetTheta2030, "OffsetTheta2030/D");
    tree->Branch("OffsetTheta2030Err", &OffsetTheta2030Err, "OffsetTheta2030Err/D");
    tree->Branch("SinAmpAllTheta", &SinAmpAllTheta, "SinAmpAllTheta/D");
    tree->Branch("SinAmpAllThetaErr", &SinAmpAllThetaErr, "SinAmpAllThetaErr/D");
    tree->Branch("SinAmpTheta010", &SinAmpTheta010, "SinAmpTheta010/D");
    tree->Branch("SinAmpTheta010Err", &SinAmpTheta010Err, "SinAmpTheta010Err/D");
    tree->Branch("SinAmpTheta1020", &SinAmpTheta1020, "SinAmpTheta1020/D");
    tree->Branch("SinAmpTheta1020Err", &SinAmpTheta1020Err, "SinAmpTheta1020Err/D");
    tree->Branch("SinAmpTheta2030", &SinAmpTheta2030, "SinAmpTheta2030/D");
    tree->Branch("SinAmpTheta2030Err", &SinAmpTheta2030Err, "SinAmpTheta2030Err/D");
    tree->Branch("CosAmpAllTheta", &CosAmpAllTheta, "CosAmpAllTheta/D");
    tree->Branch("CosAmpAllThetaErr", &CosAmpAllThetaErr, "CosAmpAllThetaErr/D");
    tree->Branch("CosAmpTheta010", &CosAmpTheta010, "CosAmpTheta010/D");
    tree->Branch("CosAmpTheta010Err", &CosAmpTheta010Err, "CosAmpTheta010Err/D");
    tree->Branch("CosAmpTheta1020", &CosAmpTheta1020, "CosAmpTheta1020/D");
    tree->Branch("CosAmpTheta1020Err", &CosAmpTheta1020Err, "CosAmpTheta1020Err/D");
    tree->Branch("CosAmpTheta2030", &CosAmpTheta2030, "CosAmpTheta2030/D");
    tree->Branch("CosAmpTheta2030Err", &CosAmpTheta2030Err, "CosAmpTheta2030Err/D");

    // Fill branches (and hence tree) with corresponding parameters from above
    for (Int_t m = 0; m < 13; m++){
      
      OffsetAllTheta = Offset[0][m];
      OffsetAllThetaErr = OffsetErr[0][m];
      OffsetTheta010 = Offset[1][m];
      OffsetTheta010Err = OffsetErr[1][m];
      OffsetTheta1020 = Offset[2][m];
      OffsetTheta1020Err = OffsetErr[2][m];
      OffsetTheta2030 = Offset[3][m];
      OffsetTheta2030Err = OffsetErr[3][m];
      SinAmpAllTheta = SinAmp[0][m];
      SinAmpAllThetaErr = SinAmpErr[0][m];
      SinAmpTheta010 = SinAmp[1][m];
      SinAmpTheta010Err = SinAmpErr[1][m];
      SinAmpTheta1020 = SinAmp[2][m];
      SinAmpTheta1020Err = SinAmpErr[2][m];
      SinAmpTheta2030 = SinAmp[3][m];
      SinAmpTheta2030Err = SinAmpErr[3][m];
      CosAmpAllTheta = CosAmp[0][m];
      CosAmpAllThetaErr = CosAmpErr[0][m];
      CosAmpTheta010 = CosAmp[1][m];
      CosAmpTheta010Err = CosAmpErr[1][m];
      CosAmpTheta1020 = CosAmp[2][m];
      CosAmpTheta1020Err = CosAmpErr[2][m];
      CosAmpTheta2030 = CosAmp[3][m];
      CosAmpTheta2030Err = CosAmpErr[3][m];
      tree->Fill();

    }

    f1.Write();

}
