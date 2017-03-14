#include "./includes.h"

// define a function with 4 parameters
Double_t fitf(Double_t *x,Double_t *par)
{
    Double_t fitval = 0;
    fitval =  par[0] + ((par[1]*sin(x[0]*TMath::DegToRad()))/(1 + (par[2]*cos(x[0]*TMath::DegToRad()))));
    return fitval;
}

void AsymmNew(){

  double InitialSinAmp[4][13];
  double InitialSinAmpErr[4][13];
  double Offset[4][13]; // Format of arrays is Theta Bin Defined by i, Energy bin defined by j
  double OffsetErr[4][13];
  double SinAmp[4][13]; // Format of arrays is Theta Bin Defined by i, Energy bin defined by j
  double SinAmpErr[4][13];
  double CosAmp[4][13]; // Format of arrays is Theta Bin Defined by i, Energy bin defined by j
  double CosAmpErr[4][13];
  Int_t i;
  double ISinAm;
  double ISinAmErr;
  double Offs;
  double OffsErr;
  double SinAm;
  double SinAmErr;
  double CosAm;
  double CosAmErr;
     
  TF1 *AsymmFunc = new TF1("AsymmFit",  fitf, -130.0, 130.0, 3); //Give a name and range to the fitting funcion
  AsymmFunc->SetParNames("Offset", "SinAmp", "CosAmp"); //Name the parameters
  AsymmFunc->SetParameter(0, 0);
  TF1 *SinFunc = new TF1("SinFit", "[0]*sin(x*TMath::DegToRad())", -130, 130);
  SinFunc->SetParNames("InitialSinAmp");
  TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/Amo/Physics_Total_16_13_3_17.root"); // Open the latest PTotal file to load histograms from
  
  for(Int_t j=0; j < 13; j++){
    
    i=0;
    
    if (j==0){
      Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (275pm25MeV)"; // Set title of output graph
      TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_275MeV_NegHel"); // Select the correct histogram
      TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_275MeV_PosHel"); // Select the correct histogram
      Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_275MeV_Asymm_16.png"; // Name the output images
      Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_275MeV_Asymm_16.pdf";
      Char_t* GraphROOT = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_275MeV_Asymm_16.root";
    }
    
    if (j==1){
      Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (325pm25MeV)"; // Set title of output graph
      TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_325MeV_NegHel"); // Select the correct histogram
      TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_325MeV_PosHel"); // Select the correct histogram
      Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_325MeV_Asymm_16.png"; // Name the output images
      Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_325MeV_Asymm_16.pdf";
      Char_t* GraphROOT = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_325MeV_Asymm_16.root";
    }
    
    if (j==2){
      Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (375pm25MeV)"; // Set title of output graph
      TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_375MeV_NegHel"); // Select the correct histogram
      TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_375MeV_PosHel"); // Select the correct histogram
      Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_375MeV_Asymm_16.png"; // Name the output images
      Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_375MeV_Asymm_16.pdf";
      Char_t* GraphROOT = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_375MeV_Asymm_16.root";
    }
    
    if (j==3){
      Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (425pm25MeV)"; // Set title of output graph
      TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_425MeV_NegHel"); // Select the correct histogram
      TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_425MeV_PosHel"); // Select the correct histogram
      Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_425MeV_Asymm_16.png"; // Name the output images
      Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_425MeV_Asymm_16.pdf";
      Char_t* GraphROOT = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_425MeV_Asymm_16.root";

    }
    
    if (j==4){
      Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (475pm25MeV)"; // Set title of output graph
      TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_475MeV_NegHel"); // Select the correct histogram
      TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_475MeV_PosHel"); // Select the correct histogram
      Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_475MeV_Asymm_16.png"; // Name the output images
      Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_475MeV_Asymm_16.pdf";
      Char_t* GraphROOT = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_475MeV_Asymm_16.root";

    }
    
    if (j==5){
      Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (525pm25MeV)"; // Set title of output graph
      TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_525MeV_NegHel"); // Select the correct histogram
      TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_525MeV_PosHel"); // Select the correct histogram
      Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_525MeV_Asymm_16.png"; // Name the output images
      Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_525MeV_Asymm_16.pdf";
      Char_t* GraphROOT = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_525MeV_Asymm_16.root";
    }
    
    if (j==6){
      Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (575pm25MeV)"; // Set title of output graph
      TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_575MeV_NegHel"); // Select the correct histogram
      TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_575MeV_PosHel"); // Select the correct histogram
      Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_575MeV_Asymm_16.png"; // Name the output images
      Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_575MeV_Asymm_16.pdf";
      Char_t* GraphROOT = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_575MeV_Asymm_16.root";
    }
    
    if (j==7){
      Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (625pm25MeV)"; // Set title of output graph
      TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_625MeV_NegHel"); // Select the correct histogram
      TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_625MeV_PosHel"); // Select the correct histogram
      Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_625MeV_Asymm_16.png"; // Name the output images
      Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_625MeV_Asymm_16.pdf";
      Char_t* GraphROOT = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_625MeV_Asymm_16.root";
    }
    
    if (j==8){
      Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (675pm25MeV)"; // Set title of output graph
      TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_675MeV_NegHel"); // Select the correct histogram
      TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_675MeV_PosHel"); // Select the correct histogram
      Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_675MeV_Asymm_16.png"; // Name the output images
      Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_675MeV_Asymm_16.pdf";
      Char_t* GraphROOT = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_675MeV_Asymm_16.root";
    }
    
    if (j==9){
      Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (725pm25MeV)"; // Set title of output graph
      TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_725MeV_NegHel"); // Select the correct histogram
      TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_725MeV_PosHel"); // Select the correct histogram
      Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_725MeV_Asymm_16.png"; // Name the output images
      Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_725MeV_Asymm_16.pdf";
      Char_t* GraphROOT = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_725MeV_Asymm_16.root";

    }
    
    if (j==10){
      Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (775pm25MeV)"; // Set title of output graph
      TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_775MeV_NegHel"); // Select the correct histogram
      TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_775MeV_PosHel"); // Select the correct histogram
      Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_775MeV_Asymm_16.png"; // Name the output images
      Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_775MeV_Asymm_16.pdf";
      Char_t* GraphROOT = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_775MeV_Asymm_16.root";
    }
    
    if (j==11){
      Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (825pm25MeV)"; // Set title of output graph
      TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_825MeV_NegHel"); // Select the correct histogram
      TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_825MeV_PosHel"); // Select the correct histogram
      Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_825MeV_Asymm_16.png"; // Name the output images
      Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_825MeV_Asymm_16.pdf";
      Char_t* GraphROOT = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_825MeV_Asymm_16.root";
    }
    
    if (j==12){
      Char_t* Title = "Asymmetry in PhiSc Between -ve and +ve Helicity (875pm25MeV)"; // Set title of output graph
      TH1D *histNeg = (TH1D*)f->Get("Phi_Scattered_875MeV_NegHel"); // Select the correct histogram
      TH1D *histPos = (TH1D*)f->Get("Phi_Scattered_875MeV_PosHel"); // Select the correct histogram
      Char_t* GraphPDF = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_875MeV_Asymm_16.png"; // Name the output images
      Char_t* GraphPNG = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_875MeV_Asymm_16.pdf";
      Char_t* GraphROOT = "/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/NegPosHel_875MeV_Asymm_16.root";
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
    Float_t yMin = -1;
    Float_t yMax = 1;
    
    strcpy(hrTitle, Title);
    hr = canvas->DrawFrame(xMin,yMin,xMax,yMax);
    hr->SetTitle(hrTitle);
    
    AsymmHist->SetMarkerStyle(1); // Style options for graph
    AsymmHist->SetLineColor(2);
    AsymmHist->Draw("EHISTSAMES"); // Draw the histogram
    AsymmHist->Fit("SinFit", "LL");
    SinFit->SetLineColor(1);
    SinFit->Draw("SAMES");
    InitialSinAmp[i][j] = SinFunc->GetParameter(0);
    InitialSinAmpErr[i][j] = SinFunc->GetParError(0);
    AsymmFunc->SetParameter(1, SinFunc->GetParameter(0));
    AsymmFunc->SetParError(1, SinFunc->GetParError(0));
    AsymmFunc->SetParLimits(0, -2, 2);
    AsymmFunc->SetParLimits(1, -1, 1);
    AsymmFunc->SetParLimits(2, -1, 1);
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
    canvas->SaveAs(filename = GraphROOT);
    canvas->Close();
    
  }
  
  // Define new file to store fit parameters
  TFile f1("SinPhiFitValues.root", "RECREATE");
  
  //Define new tree to store parameters in
  TTree* tree = new TTree("Parameter_Values", "Tree_of_Values");
  
  // Define branches to store parameters, (Branch Name, Variable, Type of Variable)
  tree->Branch("InitialSinAmp", &ISinAm, "ISinAm/D");
  tree->Branch("InitialSinAmpErr", &ISinAmErr, "ISinAmErr/D");
  tree->Branch("Offset", &Offs, "Offs/D");
  tree->Branch("OffsetErr", &OffsErr, "OffsErr/D");
  tree->Branch("SinAmp", &SinAm, "SinAm/D");
  tree->Branch("SinAmpErr", &SinAmErr, "SinAmErr/D");
  tree->Branch("CosAmp", &CosAm, "CosAm/D");
  tree->Branch("CosAmpErr", &CosAmErr, "CosAmErr/D");
  
  // Fill branches (and hence tree) with corresponding parameters from above
  for (Int_t m = 0; m < 13; m++){
    
    ISinAm = InitialSinAmp[0][m];
    ISinAmErr = InitialSinAmpErr[0][m];
    Offs = Offset[0][m];
    OffsErr = OffsetErr[0][m];
    SinAm = SinAmp[0][m];
    SinAmErr = SinAmpErr[0][m];
    CosAm = CosAmp[0][m];
    CosAmErr = CosAmpErr[0][m];
    
    tree->Fill();
    
  }
  
  f1.Write();
  
}
