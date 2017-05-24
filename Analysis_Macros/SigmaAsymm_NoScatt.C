#include "./includes_SigmaAsymm_NoScatt.h"

void SigmaAsymm_NoScatt(){

  TF1 *CosFunc = new TF1("CosFit", "[0]*cos((2*x*TMath::DegToRad())+acos(0))");
  CosFunc->SetParNames("Amplitude");

  double pCosAmp[10][20]; // Format of array is Theta bin (x) by Egamma bin (y), 6 theta bins of 30, 12 20MeV Egamma bins
  double pCosAmpErr[10][20];
  double pCosA425;
  double pCosAErr425;
  double pCosA435;
  double pCosAErr435;
  double pCosA445;
  double pCosAErr445;
  double pCosA455;
  double pCosAErr455;
  double pCosA465;
  double pCosAErr465;
  double pCosA475;
  double pCosAErr475;
  double pCosA485;
  double pCosAErr485;
  double pCosA495;
  double pCosAErr495;
  double pCosA505;
  double pCosAErr505;
  double pCosA515;
  double pCosAErr515;
  double pCosA525;
  double pCosAErr525;
  double pCosA535;
  double pCosAErr535;
  double pCosA545;
  double pCosAErr545;
  double pCosA555;
  double pCosAErr555;
  double pCosA565;
  double pCosAErr565;
  double pCosA575;
  double pCosAErr575;
  double pCosA585;
  double pCosAErr585;
  double pCosA595;
  double pCosAErr595;
  double pCosA605;
  double pCosAErr605;
  double pCosA615;
  double pCosAErr615;

  TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/ParaPerp_NoScatt_Total_7_Combined.root"); // Open the latest PTotal combined file to load histograms from
  NPara = Eg_Para->GetEntries();
  NPerp = Eg_Perp->GetEntries();
  ScaleFactor = NPara/NPerp;
  ScaleFactorErr = sqrt( (NPara/((TMath::Power(NPerp,2)))) + (((TMath::Power(NPara,2)))/(TMath::Power(NPerp,3))) ); // Error Propagation of above formula, see notebook 5

  ParaPerpAsymmPhip_425MeVCM1 = Phip_425MeVCM1_Para->GetAsymmetry(Phip_425MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_425MeVCM1->SetName("ParaPerpAsymmPhip425MeVCM1");
  ParaPerpAsymmPhip_425MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for 420-430MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_425MeVCM1->Fit("CosFit", "Q");
  pCosAmp[0][0] = CosFit->GetParameter(0);
  pCosAmpErr[0][0] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_435MeVCM1 = Phip_435MeVCM1_Para->GetAsymmetry(Phip_435MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_435MeVCM1->SetName("ParaPerpAsymmPhip435MeVCM1");
  ParaPerpAsymmPhip_435MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for 430-440MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_435MeVCM1->Fit("CosFit", "Q");
  pCosAmp[0][1] = CosFit->GetParameter(0);
  pCosAmpErr[0][1] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_445MeVCM1 = Phip_445MeVCM1_Para->GetAsymmetry(Phip_445MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_445MeVCM1->SetName("ParaPerpAsymmPhip445MeVCM1");
  ParaPerpAsymmPhip_445MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for 440-450MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_445MeVCM1->Fit("CosFit", "Q");
  pCosAmp[0][2] = CosFit->GetParameter(0);
  pCosAmpErr[0][2] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_455MeVCM1 = Phip_455MeVCM1_Para->GetAsymmetry(Phip_455MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_455MeVCM1->SetName("ParaPerpAsymmPhip455MeVCM1");
  ParaPerpAsymmPhip_455MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for 450-460MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_455MeVCM1->Fit("CosFit", "Q");
  pCosAmp[0][3] = CosFit->GetParameter(0);
  pCosAmpErr[0][3] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_465MeVCM1 = Phip_465MeVCM1_Para->GetAsymmetry(Phip_465MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_465MeVCM1->SetName("ParaPerpAsymmPhip465MeVCM1");
  ParaPerpAsymmPhip_465MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for 460-470MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_465MeVCM1->Fit("CosFit", "Q");
  pCosAmp[0][4] = CosFit->GetParameter(0);
  pCosAmpErr[0][4] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_475MeVCM1 = Phip_475MeVCM1_Para->GetAsymmetry(Phip_475MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_475MeVCM1->SetName("ParaPerpAsymmPhip475MeVCM1");
  ParaPerpAsymmPhip_475MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for 470-480MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_475MeVCM1->Fit("CosFit", "Q");
  pCosAmp[0][5] = CosFit->GetParameter(0);
  pCosAmpErr[0][5] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_485MeVCM1 = Phip_485MeVCM1_Para->GetAsymmetry(Phip_485MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_485MeVCM1->SetName("ParaPerpAsymmPhip485MeVCM1");
  ParaPerpAsymmPhip_485MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for 480-490MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_485MeVCM1->Fit("CosFit", "Q");
  pCosAmp[0][6] = CosFit->GetParameter(0);
  pCosAmpErr[0][6] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_495MeVCM1 = Phip_495MeVCM1_Para->GetAsymmetry(Phip_495MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_495MeVCM1->SetName("ParaPerpAsymmPhip495MeVCM1");
  ParaPerpAsymmPhip_495MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for 490-500MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_495MeVCM1->Fit("CosFit", "Q");
  pCosAmp[0][7] = CosFit->GetParameter(0);
  pCosAmpErr[0][7] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_505MeVCM1 = Phip_505MeVCM1_Para->GetAsymmetry(Phip_505MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_505MeVCM1->SetName("ParaPerpAsymmPhip505MeVCM1");
  ParaPerpAsymmPhip_505MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for 500-510MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_505MeVCM1->Fit("CosFit", "Q");
  pCosAmp[0][8] = CosFit->GetParameter(0);
  pCosAmpErr[0][8] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_515MeVCM1 = Phip_515MeVCM1_Para->GetAsymmetry(Phip_515MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_515MeVCM1->SetName("ParaPerpAsymmPhip515MeVCM1");
  ParaPerpAsymmPhip_515MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for 510-520MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_515MeVCM1->Fit("CosFit", "Q");
  pCosAmp[0][9] = CosFit->GetParameter(0);
  pCosAmpErr[0][9] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_525MeVCM1 = Phip_525MeVCM1_Para->GetAsymmetry(Phip_525MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_525MeVCM1->SetName("ParaPerpAsymmPhip525MeVCM1");
  ParaPerpAsymmPhip_525MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for 520-530MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_525MeVCM1->Fit("CosFit", "Q");
  pCosAmp[0][10] = CosFit->GetParameter(0);
  pCosAmpErr[0][10] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_535MeVCM1 = Phip_535MeVCM1_Para->GetAsymmetry(Phip_535MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_535MeVCM1->SetName("ParaPerpAsymmPhip535MeVCM1");
  ParaPerpAsymmPhip_535MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for 530-540MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_535MeVCM1->Fit("CosFit", "Q");
  pCosAmp[0][11] = CosFit->GetParameter(0);
  pCosAmpErr[0][11] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_545MeVCM1 = Phip_545MeVCM1_Para->GetAsymmetry(Phip_545MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_545MeVCM1->SetName("ParaPerpAsymmPhip545MeVCM1");
  ParaPerpAsymmPhip_545MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for 540-550MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_545MeVCM1->Fit("CosFit", "Q");
  pCosAmp[0][12] = CosFit->GetParameter(0);
  pCosAmpErr[0][12] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_555MeVCM1 = Phip_555MeVCM1_Para->GetAsymmetry(Phip_555MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_555MeVCM1->SetName("ParaPerpAsymmPhip555MeVCM1");
  ParaPerpAsymmPhip_555MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for 550-560MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_555MeVCM1->Fit("CosFit", "Q");
  pCosAmp[0][13] = CosFit->GetParameter(0);
  pCosAmpErr[0][13] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_565MeVCM1 = Phip_565MeVCM1_Para->GetAsymmetry(Phip_565MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_565MeVCM1->SetName("ParaPerpAsymmPhip565MeVCM1");
  ParaPerpAsymmPhip_565MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for 560-570MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_565MeVCM1->Fit("CosFit", "Q");
  pCosAmp[0][14] = CosFit->GetParameter(0);
  pCosAmpErr[0][14] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_575MeVCM1 = Phip_575MeVCM1_Para->GetAsymmetry(Phip_575MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_575MeVCM1->SetName("ParaPerpAsymmPhip575MeVCM1");
  ParaPerpAsymmPhip_575MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for 570-580MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_575MeVCM1->Fit("CosFit", "Q");
  pCosAmp[0][15] = CosFit->GetParameter(0);
  pCosAmpErr[0][15] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_585MeVCM1 = Phip_585MeVCM1_Para->GetAsymmetry(Phip_585MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_585MeVCM1->SetName("ParaPerpAsymmPhip585MeVCM1");
  ParaPerpAsymmPhip_585MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for 580-590MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_585MeVCM1->Fit("CosFit", "Q");
  pCosAmp[0][16] = CosFit->GetParameter(0);
  pCosAmpErr[0][16] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_595MeVCM1 = Phip_595MeVCM1_Para->GetAsymmetry(Phip_595MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_595MeVCM1->SetName("ParaPerpAsymmPhip595MeVCM1");
  ParaPerpAsymmPhip_595MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for 590-600MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_595MeVCM1->Fit("CosFit", "Q");
  pCosAmp[0][17] = CosFit->GetParameter(0);
  pCosAmpErr[0][17] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_605MeVCM1 = Phip_605MeVCM1_Para->GetAsymmetry(Phip_605MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_605MeVCM1->SetName("ParaPerpAsymmPhip605MeVCM1");
  ParaPerpAsymmPhip_605MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for 600-610MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_605MeVCM1->Fit("CosFit", "Q");
  pCosAmp[0][18] = CosFit->GetParameter(0);
  pCosAmpErr[0][18] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_615MeVCM1 = Phip_615MeVCM1_Para->GetAsymmetry(Phip_615MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_615MeVCM1->SetName("ParaPerpAsymmPhip615MeVCM1");
  ParaPerpAsymmPhip_615MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for 610-620MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_615MeVCM1->Fit("CosFit", "Q");
  pCosAmp[0][19] = CosFit->GetParameter(0);
  pCosAmpErr[0][19] = CosFit->GetParError(0);

  ///////////////////////////////////////////
  //////////////////CM2 p ///////////////////
  ///////////////////////////////////////////

  ParaPerpAsymmPhip_425MeVCM2 = Phip_425MeVCM2_Para->GetAsymmetry(Phip_425MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_425MeVCM2->SetName("ParaPerpAsymmPhip425MeVCM2");
  ParaPerpAsymmPhip_425MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for 420-430MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_425MeVCM2->Fit("CosFit", "Q");
  pCosAmp[1][0] = CosFit->GetParameter(0);
  pCosAmpErr[1][0] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_435MeVCM2 = Phip_435MeVCM2_Para->GetAsymmetry(Phip_435MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_435MeVCM2->SetName("ParaPerpAsymmPhip435MeVCM2");
  ParaPerpAsymmPhip_435MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for 430-440MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_435MeVCM2->Fit("CosFit", "Q");
  pCosAmp[1][1] = CosFit->GetParameter(0);
  pCosAmpErr[1][1] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_445MeVCM2 = Phip_445MeVCM2_Para->GetAsymmetry(Phip_445MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_445MeVCM2->SetName("ParaPerpAsymmPhip445MeVCM2");
  ParaPerpAsymmPhip_445MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for 440-450MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_445MeVCM2->Fit("CosFit", "Q");
  pCosAmp[1][2] = CosFit->GetParameter(0);
  pCosAmpErr[1][2] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_455MeVCM2 = Phip_455MeVCM2_Para->GetAsymmetry(Phip_455MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_455MeVCM2->SetName("ParaPerpAsymmPhip455MeVCM2");
  ParaPerpAsymmPhip_455MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for 450-460MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_455MeVCM2->Fit("CosFit", "Q");
  pCosAmp[1][3] = CosFit->GetParameter(0);
  pCosAmpErr[1][3] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_465MeVCM2 = Phip_465MeVCM2_Para->GetAsymmetry(Phip_465MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_465MeVCM2->SetName("ParaPerpAsymmPhip465MeVCM2");
  ParaPerpAsymmPhip_465MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for 460-470MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_465MeVCM2->Fit("CosFit", "Q");
  pCosAmp[1][4] = CosFit->GetParameter(0);
  pCosAmpErr[1][4] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_475MeVCM2 = Phip_475MeVCM2_Para->GetAsymmetry(Phip_475MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_475MeVCM2->SetName("ParaPerpAsymmPhip475MeVCM2");
  ParaPerpAsymmPhip_475MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for 470-480MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_475MeVCM2->Fit("CosFit", "Q");
  pCosAmp[1][5] = CosFit->GetParameter(0);
  pCosAmpErr[1][5] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_485MeVCM2 = Phip_485MeVCM2_Para->GetAsymmetry(Phip_485MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_485MeVCM2->SetName("ParaPerpAsymmPhip485MeVCM2");
  ParaPerpAsymmPhip_485MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for 480-490MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_485MeVCM2->Fit("CosFit", "Q");
  pCosAmp[1][6] = CosFit->GetParameter(0);
  pCosAmpErr[1][6] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_495MeVCM2 = Phip_495MeVCM2_Para->GetAsymmetry(Phip_495MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_495MeVCM2->SetName("ParaPerpAsymmPhip495MeVCM2");
  ParaPerpAsymmPhip_495MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for 490-500MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_495MeVCM2->Fit("CosFit", "Q");
  pCosAmp[1][7] = CosFit->GetParameter(0);
  pCosAmpErr[1][7] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_505MeVCM2 = Phip_505MeVCM2_Para->GetAsymmetry(Phip_505MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_505MeVCM2->SetName("ParaPerpAsymmPhip505MeVCM2");
  ParaPerpAsymmPhip_505MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for 500-510MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_505MeVCM2->Fit("CosFit", "Q");
  pCosAmp[1][8] = CosFit->GetParameter(0);
  pCosAmpErr[1][8] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_515MeVCM2 = Phip_515MeVCM2_Para->GetAsymmetry(Phip_515MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_515MeVCM2->SetName("ParaPerpAsymmPhip515MeVCM2");
  ParaPerpAsymmPhip_515MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for 510-520MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_515MeVCM2->Fit("CosFit", "Q");
  pCosAmp[1][9] = CosFit->GetParameter(0);
  pCosAmpErr[1][9] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_525MeVCM2 = Phip_525MeVCM2_Para->GetAsymmetry(Phip_525MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_525MeVCM2->SetName("ParaPerpAsymmPhip525MeVCM2");
  ParaPerpAsymmPhip_525MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for 520-530MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_525MeVCM2->Fit("CosFit", "Q");
  pCosAmp[1][10] = CosFit->GetParameter(0);
  pCosAmpErr[1][10] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_535MeVCM2 = Phip_535MeVCM2_Para->GetAsymmetry(Phip_535MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_535MeVCM2->SetName("ParaPerpAsymmPhip535MeVCM2");
  ParaPerpAsymmPhip_535MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for 530-540MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_535MeVCM2->Fit("CosFit", "Q");
  pCosAmp[1][11] = CosFit->GetParameter(0);
  pCosAmpErr[1][11] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_545MeVCM2 = Phip_545MeVCM2_Para->GetAsymmetry(Phip_545MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_545MeVCM2->SetName("ParaPerpAsymmPhip545MeVCM2");
  ParaPerpAsymmPhip_545MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for 540-550MeV Photon Energy (CosTheta0.8-0.6)");
  ParaPerpAsymmPhip_545MeVCM2->Fit("CosFit", "Q");
  pCosAmp[1][12] = CosFit->GetParameter(0);
  pCosAmpErr[1][12] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_555MeVCM2 = Phip_555MeVCM2_Para->GetAsymmetry(Phip_555MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_555MeVCM2->SetName("ParaPerpAsymmPhip555MeVCM2");
  ParaPerpAsymmPhip_555MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for 550-560MeV Photon Energy (CosTheta0.8-0.6)");
  ParaPerpAsymmPhip_555MeVCM2->Fit("CosFit", "Q");
  pCosAmp[1][13] = CosFit->GetParameter(0);
  pCosAmpErr[1][13] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_565MeVCM2 = Phip_565MeVCM2_Para->GetAsymmetry(Phip_565MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_565MeVCM2->SetName("ParaPerpAsymmPhip565MeVCM2");
  ParaPerpAsymmPhip_565MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for 560-570MeV Photon Energy (CosTheta0.8-0.6)");
  ParaPerpAsymmPhip_565MeVCM2->Fit("CosFit", "Q");
  pCosAmp[1][14] = CosFit->GetParameter(0);
  pCosAmpErr[1][14] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_575MeVCM2 = Phip_575MeVCM2_Para->GetAsymmetry(Phip_575MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_575MeVCM2->SetName("ParaPerpAsymmPhip575MeVCM2");
  ParaPerpAsymmPhip_575MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for 570-580MeV Photon Energy (CosTheta0.8-0.6)");
  ParaPerpAsymmPhip_575MeVCM2->Fit("CosFit", "Q");
  pCosAmp[1][15] = CosFit->GetParameter(0);
  pCosAmpErr[1][15] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_585MeVCM2 = Phip_585MeVCM2_Para->GetAsymmetry(Phip_585MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_585MeVCM2->SetName("ParaPerpAsymmPhip585MeVCM2");
  ParaPerpAsymmPhip_585MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for 580-590MeV Photon Energy (CosTheta0.8-0.6)");
  ParaPerpAsymmPhip_585MeVCM2->Fit("CosFit", "Q");
  pCosAmp[1][16] = CosFit->GetParameter(0);
  pCosAmpErr[1][16] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_595MeVCM2 = Phip_595MeVCM2_Para->GetAsymmetry(Phip_595MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_595MeVCM2->SetName("ParaPerpAsymmPhip595MeVCM2");
  ParaPerpAsymmPhip_595MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for 590-600MeV Photon Energy (CosTheta0.8-0.6)");
  ParaPerpAsymmPhip_595MeVCM2->Fit("CosFit", "Q");
  pCosAmp[1][17] = CosFit->GetParameter(0);
  pCosAmpErr[1][17] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_605MeVCM2 = Phip_605MeVCM2_Para->GetAsymmetry(Phip_605MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_605MeVCM2->SetName("ParaPerpAsymmPhip605MeVCM2");
  ParaPerpAsymmPhip_605MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for 600-610MeV Photon Energy (CosTheta0.8-0.6)");
  ParaPerpAsymmPhip_605MeVCM2->Fit("CosFit", "Q");
  pCosAmp[1][18] = CosFit->GetParameter(0);
  pCosAmpErr[1][18] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_615MeVCM2 = Phip_615MeVCM2_Para->GetAsymmetry(Phip_615MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_615MeVCM2->SetName("ParaPerpAsymmPhip615MeVCM2");
  ParaPerpAsymmPhip_615MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for 610-620MeV Photon Energy (CosTheta0.8-0.6)");
  ParaPerpAsymmPhip_615MeVCM2->Fit("CosFit", "Q");
  pCosAmp[1][19] = CosFit->GetParameter(0);
  pCosAmpErr[1][19] = CosFit->GetParError(0);

  ///////////////////////////////////////////
  //////////////////CM3 p ///////////////////
  ///////////////////////////////////////////

  ParaPerpAsymmPhip_425MeVCM3 = Phip_425MeVCM3_Para->GetAsymmetry(Phip_425MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_425MeVCM3->SetName("ParaPerpAsymmPhip425MeVCM3");
  ParaPerpAsymmPhip_425MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for 420-430MeV Photon Energy (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_425MeVCM3->Fit("CosFit", "Q");
  pCosAmp[2][0] = CosFit->GetParameter(0);
  pCosAmpErr[2][0] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_435MeVCM3 = Phip_435MeVCM3_Para->GetAsymmetry(Phip_435MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_435MeVCM3->SetName("ParaPerpAsymmPhip435MeVCM3");
  ParaPerpAsymmPhip_435MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for 430-440MeV Photon Energy (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_435MeVCM3->Fit("CosFit", "Q");
  pCosAmp[2][1] = CosFit->GetParameter(0);
  pCosAmpErr[2][1] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_445MeVCM3 = Phip_445MeVCM3_Para->GetAsymmetry(Phip_445MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_445MeVCM3->SetName("ParaPerpAsymmPhip445MeVCM3");
  ParaPerpAsymmPhip_445MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for 440-450MeV Photon Energy (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_445MeVCM3->Fit("CosFit", "Q");
  pCosAmp[2][2] = CosFit->GetParameter(0);
  pCosAmpErr[2][2] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_455MeVCM3 = Phip_455MeVCM3_Para->GetAsymmetry(Phip_455MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_455MeVCM3->SetName("ParaPerpAsymmPhip455MeVCM3");
  ParaPerpAsymmPhip_455MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for 450-460MeV Photon Energy (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_455MeVCM3->Fit("CosFit", "Q");
  pCosAmp[2][3] = CosFit->GetParameter(0);
  pCosAmpErr[2][3] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_465MeVCM3 = Phip_465MeVCM3_Para->GetAsymmetry(Phip_465MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_465MeVCM3->SetName("ParaPerpAsymmPhip465MeVCM3");
  ParaPerpAsymmPhip_465MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for 460-470MeV Photon Energy (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_465MeVCM3->Fit("CosFit", "Q");
  pCosAmp[2][4] = CosFit->GetParameter(0);
  pCosAmpErr[2][4] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_475MeVCM3 = Phip_475MeVCM3_Para->GetAsymmetry(Phip_475MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_475MeVCM3->SetName("ParaPerpAsymmPhip475MeVCM3");
  ParaPerpAsymmPhip_475MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for 470-480MeV Photon Energy (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_475MeVCM3->Fit("CosFit", "Q");
  pCosAmp[2][5] = CosFit->GetParameter(0);
  pCosAmpErr[2][5] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_485MeVCM3 = Phip_485MeVCM3_Para->GetAsymmetry(Phip_485MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_485MeVCM3->SetName("ParaPerpAsymmPhip485MeVCM3");
  ParaPerpAsymmPhip_485MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for 480-490MeV Photon Energy (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_485MeVCM3->Fit("CosFit", "Q");
  pCosAmp[2][6] = CosFit->GetParameter(0);
  pCosAmpErr[2][6] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_495MeVCM3 = Phip_495MeVCM3_Para->GetAsymmetry(Phip_495MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_495MeVCM3->SetName("ParaPerpAsymmPhip495MeVCM3");
  ParaPerpAsymmPhip_495MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for 490-500MeV Photon Energy (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_495MeVCM3->Fit("CosFit", "Q");
  pCosAmp[2][7] = CosFit->GetParameter(0);
  pCosAmpErr[2][7] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_505MeVCM3 = Phip_505MeVCM3_Para->GetAsymmetry(Phip_505MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_505MeVCM3->SetName("ParaPerpAsymmPhip505MeVCM3");
  ParaPerpAsymmPhip_505MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for 500-510MeV Photon Energy (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_505MeVCM3->Fit("CosFit", "Q");
  pCosAmp[2][8] = CosFit->GetParameter(0);
  pCosAmpErr[2][8] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_515MeVCM3 = Phip_515MeVCM3_Para->GetAsymmetry(Phip_515MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_515MeVCM3->SetName("ParaPerpAsymmPhip515MeVCM3");
  ParaPerpAsymmPhip_515MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for 510-520MeV Photon Energy (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_515MeVCM3->Fit("CosFit", "Q");
  pCosAmp[2][9] = CosFit->GetParameter(0);
  pCosAmpErr[2][9] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_525MeVCM3 = Phip_525MeVCM3_Para->GetAsymmetry(Phip_525MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_525MeVCM3->SetName("ParaPerpAsymmPhip525MeVCM3");
  ParaPerpAsymmPhip_525MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for 520-530MeV Photon Energy (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_525MeVCM3->Fit("CosFit", "Q");
  pCosAmp[2][10] = CosFit->GetParameter(0);
  pCosAmpErr[2][10] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_535MeVCM3 = Phip_535MeVCM3_Para->GetAsymmetry(Phip_535MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_535MeVCM3->SetName("ParaPerpAsymmPhip535MeVCM3");
  ParaPerpAsymmPhip_535MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for 530-540MeV Photon Energy (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_535MeVCM3->Fit("CosFit", "Q");
  pCosAmp[2][11] = CosFit->GetParameter(0);
  pCosAmpErr[2][11] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_545MeVCM3 = Phip_545MeVCM3_Para->GetAsymmetry(Phip_545MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_545MeVCM3->SetName("ParaPerpAsymmPhip545MeVCM3");
  ParaPerpAsymmPhip_545MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for 540-550MeV Photon Energy (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_545MeVCM3->Fit("CosFit", "Q");
  pCosAmp[2][12] = CosFit->GetParameter(0);
  pCosAmpErr[2][12] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_555MeVCM3 = Phip_555MeVCM3_Para->GetAsymmetry(Phip_555MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_555MeVCM3->SetName("ParaPerpAsymmPhip555MeVCM3");
  ParaPerpAsymmPhip_555MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for 550-560MeV Photon Energy (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_555MeVCM3->Fit("CosFit", "Q");
  pCosAmp[2][13] = CosFit->GetParameter(0);
  pCosAmpErr[2][13] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_565MeVCM3 = Phip_565MeVCM3_Para->GetAsymmetry(Phip_565MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_565MeVCM3->SetName("ParaPerpAsymmPhip565MeVCM3");
  ParaPerpAsymmPhip_565MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for 560-570MeV Photon Energy (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_565MeVCM3->Fit("CosFit", "Q");
  pCosAmp[2][14] = CosFit->GetParameter(0);
  pCosAmpErr[2][14] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_575MeVCM3 = Phip_575MeVCM3_Para->GetAsymmetry(Phip_575MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_575MeVCM3->SetName("ParaPerpAsymmPhip575MeVCM3");
  ParaPerpAsymmPhip_575MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for 570-580MeV Photon Energy (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_575MeVCM3->Fit("CosFit", "Q");
  pCosAmp[2][15] = CosFit->GetParameter(0);
  pCosAmpErr[2][15] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_585MeVCM3 = Phip_585MeVCM3_Para->GetAsymmetry(Phip_585MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_585MeVCM3->SetName("ParaPerpAsymmPhip585MeVCM3");
  ParaPerpAsymmPhip_585MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for 580-590MeV Photon Energy (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_585MeVCM3->Fit("CosFit", "Q");
  pCosAmp[2][16] = CosFit->GetParameter(0);
  pCosAmpErr[2][16] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_595MeVCM3 = Phip_595MeVCM3_Para->GetAsymmetry(Phip_595MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_595MeVCM3->SetName("ParaPerpAsymmPhip595MeVCM3");
  ParaPerpAsymmPhip_595MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for 590-600MeV Photon Energy (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_595MeVCM3->Fit("CosFit", "Q");
  pCosAmp[2][17] = CosFit->GetParameter(0);
  pCosAmpErr[2][17] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_605MeVCM3 = Phip_605MeVCM3_Para->GetAsymmetry(Phip_605MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_605MeVCM3->SetName("ParaPerpAsymmPhip605MeVCM3");
  ParaPerpAsymmPhip_605MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for 600-610MeV Photon Energy (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_605MeVCM3->Fit("CosFit", "Q");
  pCosAmp[2][18] = CosFit->GetParameter(0);
  pCosAmpErr[2][18] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_615MeVCM3 = Phip_615MeVCM3_Para->GetAsymmetry(Phip_615MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_615MeVCM3->SetName("ParaPerpAsymmPhip615MeVCM3");
  ParaPerpAsymmPhip_615MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for 610-620MeV Photon Energy (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_615MeVCM3->Fit("CosFit", "Q");
  pCosAmp[2][19] = CosFit->GetParameter(0);
  pCosAmpErr[2][19] = CosFit->GetParError(0);

  ///////////////////////////////////////////
  //////////////////CM4 p ///////////////////
  ///////////////////////////////////////////

  ParaPerpAsymmPhip_425MeVCM4 = Phip_425MeVCM4_Para->GetAsymmetry(Phip_425MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_425MeVCM4->SetName("ParaPerpAsymmPhip425MeVCM4");
  ParaPerpAsymmPhip_425MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for 420-430MeV Photon Energy (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_425MeVCM4->Fit("CosFit", "Q");
  pCosAmp[3][0] = CosFit->GetParameter(0);
  pCosAmpErr[3][0] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_435MeVCM4 = Phip_435MeVCM4_Para->GetAsymmetry(Phip_435MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_435MeVCM4->SetName("ParaPerpAsymmPhip435MeVCM4");
  ParaPerpAsymmPhip_435MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for 430-440MeV Photon Energy (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_435MeVCM4->Fit("CosFit", "Q");
  pCosAmp[3][1] = CosFit->GetParameter(0);
  pCosAmpErr[3][1] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_445MeVCM4 = Phip_445MeVCM4_Para->GetAsymmetry(Phip_445MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_445MeVCM4->SetName("ParaPerpAsymmPhip445MeVCM4");
  ParaPerpAsymmPhip_445MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for 440-450MeV Photon Energy (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_445MeVCM4->Fit("CosFit", "Q");
  pCosAmp[3][2] = CosFit->GetParameter(0);
  pCosAmpErr[3][2] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_455MeVCM4 = Phip_455MeVCM4_Para->GetAsymmetry(Phip_455MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_455MeVCM4->SetName("ParaPerpAsymmPhip455MeVCM4");
  ParaPerpAsymmPhip_455MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for 450-460MeV Photon Energy (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_455MeVCM4->Fit("CosFit", "Q");
  pCosAmp[3][3] = CosFit->GetParameter(0);
  pCosAmpErr[3][3] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_465MeVCM4 = Phip_465MeVCM4_Para->GetAsymmetry(Phip_465MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_465MeVCM4->SetName("ParaPerpAsymmPhip465MeVCM4");
  ParaPerpAsymmPhip_465MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for 460-470MeV Photon Energy (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_465MeVCM4->Fit("CosFit", "Q");
  pCosAmp[3][4] = CosFit->GetParameter(0);
  pCosAmpErr[3][4] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_475MeVCM4 = Phip_475MeVCM4_Para->GetAsymmetry(Phip_475MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_475MeVCM4->SetName("ParaPerpAsymmPhip475MeVCM4");
  ParaPerpAsymmPhip_475MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for 470-480MeV Photon Energy (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_475MeVCM4->Fit("CosFit", "Q");
  pCosAmp[3][5] = CosFit->GetParameter(0);
  pCosAmpErr[3][5] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_485MeVCM4 = Phip_485MeVCM4_Para->GetAsymmetry(Phip_485MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_485MeVCM4->SetName("ParaPerpAsymmPhip485MeVCM4");
  ParaPerpAsymmPhip_485MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for 480-490MeV Photon Energy (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_485MeVCM4->Fit("CosFit", "Q");
  pCosAmp[3][6] = CosFit->GetParameter(0);
  pCosAmpErr[3][6] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_495MeVCM4 = Phip_495MeVCM4_Para->GetAsymmetry(Phip_495MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_495MeVCM4->SetName("ParaPerpAsymmPhip495MeVCM4");
  ParaPerpAsymmPhip_495MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for 490-500MeV Photon Energy (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_495MeVCM4->Fit("CosFit", "Q");
  pCosAmp[3][7] = CosFit->GetParameter(0);
  pCosAmpErr[3][7] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_505MeVCM4 = Phip_505MeVCM4_Para->GetAsymmetry(Phip_505MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_505MeVCM4->SetName("ParaPerpAsymmPhip505MeVCM4");
  ParaPerpAsymmPhip_505MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for 500-510MeV Photon Energy (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_505MeVCM4->Fit("CosFit", "Q");
  pCosAmp[3][8] = CosFit->GetParameter(0);
  pCosAmpErr[3][8] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_515MeVCM4 = Phip_515MeVCM4_Para->GetAsymmetry(Phip_515MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_515MeVCM4->SetName("ParaPerpAsymmPhip515MeVCM4");
  ParaPerpAsymmPhip_515MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for 510-520MeV Photon Energy (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_515MeVCM4->Fit("CosFit", "Q");
  pCosAmp[3][9] = CosFit->GetParameter(0);
  pCosAmpErr[3][9] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_525MeVCM4 = Phip_525MeVCM4_Para->GetAsymmetry(Phip_525MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_525MeVCM4->SetName("ParaPerpAsymmPhip525MeVCM4");
  ParaPerpAsymmPhip_525MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for 520-530MeV Photon Energy (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_525MeVCM4->Fit("CosFit", "Q");
  pCosAmp[3][10] = CosFit->GetParameter(0);
  pCosAmpErr[3][10] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_535MeVCM4 = Phip_535MeVCM4_Para->GetAsymmetry(Phip_535MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_535MeVCM4->SetName("ParaPerpAsymmPhip535MeVCM4");
  ParaPerpAsymmPhip_535MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for 530-540MeV Photon Energy (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_535MeVCM4->Fit("CosFit", "Q");
  pCosAmp[3][11] = CosFit->GetParameter(0);
  pCosAmpErr[3][11] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_545MeVCM4 = Phip_545MeVCM4_Para->GetAsymmetry(Phip_545MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_545MeVCM4->SetName("ParaPerpAsymmPhip545MeVCM4");
  ParaPerpAsymmPhip_545MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for 540-550MeV Photon Energy (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_545MeVCM4->Fit("CosFit", "Q");
  pCosAmp[3][12] = CosFit->GetParameter(0);
  pCosAmpErr[3][12] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_555MeVCM4 = Phip_555MeVCM4_Para->GetAsymmetry(Phip_555MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_555MeVCM4->SetName("ParaPerpAsymmPhip555MeVCM4");
  ParaPerpAsymmPhip_555MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for 550-560MeV Photon Energy (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_555MeVCM4->Fit("CosFit", "Q");
  pCosAmp[3][13] = CosFit->GetParameter(0);
  pCosAmpErr[3][13] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_565MeVCM4 = Phip_565MeVCM4_Para->GetAsymmetry(Phip_565MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_565MeVCM4->SetName("ParaPerpAsymmPhip565MeVCM4");
  ParaPerpAsymmPhip_565MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for 560-570MeV Photon Energy (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_565MeVCM4->Fit("CosFit", "Q");
  pCosAmp[3][14] = CosFit->GetParameter(0);
  pCosAmpErr[3][14] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_575MeVCM4 = Phip_575MeVCM4_Para->GetAsymmetry(Phip_575MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_575MeVCM4->SetName("ParaPerpAsymmPhip575MeVCM4");
  ParaPerpAsymmPhip_575MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for 570-580MeV Photon Energy (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_575MeVCM4->Fit("CosFit", "Q");
  pCosAmp[3][15] = CosFit->GetParameter(0);
  pCosAmpErr[3][15] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_585MeVCM4 = Phip_585MeVCM4_Para->GetAsymmetry(Phip_585MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_585MeVCM4->SetName("ParaPerpAsymmPhip585MeVCM4");
  ParaPerpAsymmPhip_585MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for 580-590MeV Photon Energy (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_585MeVCM4->Fit("CosFit", "Q");
  pCosAmp[3][16] = CosFit->GetParameter(0);
  pCosAmpErr[3][16] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_595MeVCM4 = Phip_595MeVCM4_Para->GetAsymmetry(Phip_595MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_595MeVCM4->SetName("ParaPerpAsymmPhip595MeVCM4");
  ParaPerpAsymmPhip_595MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for 590-600MeV Photon Energy (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_595MeVCM4->Fit("CosFit", "Q");
  pCosAmp[3][17] = CosFit->GetParameter(0);
  pCosAmpErr[3][17] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_605MeVCM4 = Phip_605MeVCM4_Para->GetAsymmetry(Phip_605MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_605MeVCM4->SetName("ParaPerpAsymmPhip605MeVCM4");
  ParaPerpAsymmPhip_605MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for 600-610MeV Photon Energy (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_605MeVCM4->Fit("CosFit", "Q");
  pCosAmp[3][18] = CosFit->GetParameter(0);
  pCosAmpErr[3][18] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_615MeVCM4 = Phip_615MeVCM4_Para->GetAsymmetry(Phip_615MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_615MeVCM4->SetName("ParaPerpAsymmPhip615MeVCM4");
  ParaPerpAsymmPhip_615MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for 610-620MeV Photon Energy (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_615MeVCM4->Fit("CosFit", "Q");
  pCosAmp[3][19] = CosFit->GetParameter(0);
  pCosAmpErr[3][19] = CosFit->GetParError(0);

  ///////////////////////////////////////////
  //////////////////CM5 p ///////////////////
  ///////////////////////////////////////////

  ParaPerpAsymmPhip_425MeVCM5 = Phip_425MeVCM5_Para->GetAsymmetry(Phip_425MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_425MeVCM5->SetName("ParaPerpAsymmPhip425MeVCM5");
  ParaPerpAsymmPhip_425MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for 420-430MeV Photon Energy (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_425MeVCM5->Fit("CosFit", "Q");
  pCosAmp[4][0] = CosFit->GetParameter(0);
  pCosAmpErr[4][0] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_435MeVCM5 = Phip_435MeVCM5_Para->GetAsymmetry(Phip_435MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_435MeVCM5->SetName("ParaPerpAsymmPhip435MeVCM5");
  ParaPerpAsymmPhip_435MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for 430-440MeV Photon Energy (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_435MeVCM5->Fit("CosFit", "Q");
  pCosAmp[4][1] = CosFit->GetParameter(0);
  pCosAmpErr[4][1] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_445MeVCM5 = Phip_445MeVCM5_Para->GetAsymmetry(Phip_445MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_445MeVCM5->SetName("ParaPerpAsymmPhip445MeVCM5");
  ParaPerpAsymmPhip_445MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for 440-450MeV Photon Energy (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_445MeVCM5->Fit("CosFit", "Q");
  pCosAmp[4][2] = CosFit->GetParameter(0);
  pCosAmpErr[4][2] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_455MeVCM5 = Phip_455MeVCM5_Para->GetAsymmetry(Phip_455MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_455MeVCM5->SetName("ParaPerpAsymmPhip455MeVCM5");
  ParaPerpAsymmPhip_455MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for 450-460MeV Photon Energy (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_455MeVCM5->Fit("CosFit", "Q");
  pCosAmp[4][3] = CosFit->GetParameter(0);
  pCosAmpErr[4][3] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_465MeVCM5 = Phip_465MeVCM5_Para->GetAsymmetry(Phip_465MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_465MeVCM5->SetName("ParaPerpAsymmPhip465MeVCM5");
  ParaPerpAsymmPhip_465MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for 460-470MeV Photon Energy (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_465MeVCM5->Fit("CosFit", "Q");
  pCosAmp[4][4] = CosFit->GetParameter(0);
  pCosAmpErr[4][4] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_475MeVCM5 = Phip_475MeVCM5_Para->GetAsymmetry(Phip_475MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_475MeVCM5->SetName("ParaPerpAsymmPhip475MeVCM5");
  ParaPerpAsymmPhip_475MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for 470-480MeV Photon Energy (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_475MeVCM5->Fit("CosFit", "Q");
  pCosAmp[4][5] = CosFit->GetParameter(0);
  pCosAmpErr[4][5] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_485MeVCM5 = Phip_485MeVCM5_Para->GetAsymmetry(Phip_485MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_485MeVCM5->SetName("ParaPerpAsymmPhip485MeVCM5");
  ParaPerpAsymmPhip_485MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for 480-490MeV Photon Energy (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_485MeVCM5->Fit("CosFit", "Q");
  pCosAmp[4][6] = CosFit->GetParameter(0);
  pCosAmpErr[4][6] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_495MeVCM5 = Phip_495MeVCM5_Para->GetAsymmetry(Phip_495MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_495MeVCM5->SetName("ParaPerpAsymmPhip495MeVCM5");
  ParaPerpAsymmPhip_495MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for 490-500MeV Photon Energy (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_495MeVCM5->Fit("CosFit", "Q");
  pCosAmp[4][7] = CosFit->GetParameter(0);
  pCosAmpErr[4][7] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_505MeVCM5 = Phip_505MeVCM5_Para->GetAsymmetry(Phip_505MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_505MeVCM5->SetName("ParaPerpAsymmPhip505MeVCM5");
  ParaPerpAsymmPhip_505MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for 500-510MeV Photon Energy (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_505MeVCM5->Fit("CosFit", "Q");
  pCosAmp[4][8] = CosFit->GetParameter(0);
  pCosAmpErr[4][8] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_515MeVCM5 = Phip_515MeVCM5_Para->GetAsymmetry(Phip_515MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_515MeVCM5->SetName("ParaPerpAsymmPhip515MeVCM5");
  ParaPerpAsymmPhip_515MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for 510-520MeV Photon Energy (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_515MeVCM5->Fit("CosFit", "Q");
  pCosAmp[4][9] = CosFit->GetParameter(0);
  pCosAmpErr[4][9] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_525MeVCM5 = Phip_525MeVCM5_Para->GetAsymmetry(Phip_525MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_525MeVCM5->SetName("ParaPerpAsymmPhip525MeVCM5");
  ParaPerpAsymmPhip_525MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for 520-530MeV Photon Energy (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_525MeVCM5->Fit("CosFit", "Q");
  pCosAmp[4][10] = CosFit->GetParameter(0);
  pCosAmpErr[4][10] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_535MeVCM5 = Phip_535MeVCM5_Para->GetAsymmetry(Phip_535MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_535MeVCM5->SetName("ParaPerpAsymmPhip535MeVCM5");
  ParaPerpAsymmPhip_535MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for 530-540MeV Photon Energy (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_535MeVCM5->Fit("CosFit", "Q");
  pCosAmp[4][11] = CosFit->GetParameter(0);
  pCosAmpErr[4][11] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_545MeVCM5 = Phip_545MeVCM5_Para->GetAsymmetry(Phip_545MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_545MeVCM5->SetName("ParaPerpAsymmPhip545MeVCM5");
  ParaPerpAsymmPhip_545MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for 540-550MeV Photon Energy (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_545MeVCM5->Fit("CosFit", "Q");
  pCosAmp[4][12] = CosFit->GetParameter(0);
  pCosAmpErr[4][12] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_555MeVCM5 = Phip_555MeVCM5_Para->GetAsymmetry(Phip_555MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_555MeVCM5->SetName("ParaPerpAsymmPhip555MeVCM5");
  ParaPerpAsymmPhip_555MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for 550-560MeV Photon Energy (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_555MeVCM5->Fit("CosFit", "Q");
  pCosAmp[4][13] = CosFit->GetParameter(0);
  pCosAmpErr[4][13] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_565MeVCM5 = Phip_565MeVCM5_Para->GetAsymmetry(Phip_565MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_565MeVCM5->SetName("ParaPerpAsymmPhip565MeVCM5");
  ParaPerpAsymmPhip_565MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for 560-570MeV Photon Energy (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_565MeVCM5->Fit("CosFit", "Q");
  pCosAmp[4][14] = CosFit->GetParameter(0);
  pCosAmpErr[4][14] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_575MeVCM5 = Phip_575MeVCM5_Para->GetAsymmetry(Phip_575MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_575MeVCM5->SetName("ParaPerpAsymmPhip575MeVCM5");
  ParaPerpAsymmPhip_575MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for 570-580MeV Photon Energy (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_575MeVCM5->Fit("CosFit", "Q");
  pCosAmp[4][15] = CosFit->GetParameter(0);
  pCosAmpErr[4][15] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_585MeVCM5 = Phip_585MeVCM5_Para->GetAsymmetry(Phip_585MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_585MeVCM5->SetName("ParaPerpAsymmPhip585MeVCM5");
  ParaPerpAsymmPhip_585MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for 580-590MeV Photon Energy (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_585MeVCM5->Fit("CosFit", "Q");
  pCosAmp[4][16] = CosFit->GetParameter(0);
  pCosAmpErr[4][16] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_595MeVCM5 = Phip_595MeVCM5_Para->GetAsymmetry(Phip_595MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_595MeVCM5->SetName("ParaPerpAsymmPhip595MeVCM5");
  ParaPerpAsymmPhip_595MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for 590-600MeV Photon Energy (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_595MeVCM5->Fit("CosFit", "Q");
  pCosAmp[4][17] = CosFit->GetParameter(0);
  pCosAmpErr[4][17] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_605MeVCM5 = Phip_605MeVCM5_Para->GetAsymmetry(Phip_605MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_605MeVCM5->SetName("ParaPerpAsymmPhip605MeVCM5");
  ParaPerpAsymmPhip_605MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for 600-610MeV Photon Energy (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_605MeVCM5->Fit("CosFit", "Q");
  pCosAmp[4][18] = CosFit->GetParameter(0);
  pCosAmpErr[4][18] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_615MeVCM5 = Phip_615MeVCM5_Para->GetAsymmetry(Phip_615MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_615MeVCM5->SetName("ParaPerpAsymmPhip615MeVCM5");
  ParaPerpAsymmPhip_615MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for 610-620MeV Photon Energy (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_615MeVCM5->Fit("CosFit", "Q");
  pCosAmp[4][19] = CosFit->GetParameter(0);
  pCosAmpErr[4][19] = CosFit->GetParError(0);

  ///////////////////////////////////////////
  //////////////////CM6 p ///////////////////
  ///////////////////////////////////////////

  ParaPerpAsymmPhip_425MeVCM6 = Phip_425MeVCM6_Para->GetAsymmetry(Phip_425MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_425MeVCM6->SetName("ParaPerpAsymmPhip425MeVCM6");
  ParaPerpAsymmPhip_425MeVCM6->SetTitle("Proton Para/Perp Phi Asymmetry for 420-430MeV Photon Energy (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_425MeVCM6->Fit("CosFit", "Q");
  pCosAmp[5][0] = CosFit->GetParameter(0);
  pCosAmpErr[5][0] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_435MeVCM6 = Phip_435MeVCM6_Para->GetAsymmetry(Phip_435MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_435MeVCM6->SetName("ParaPerpAsymmPhip435MeVCM6");
  ParaPerpAsymmPhip_435MeVCM6->SetTitle("Proton Para/Perp Phi Asymmetry for 430-440MeV Photon Energy (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_435MeVCM6->Fit("CosFit", "Q");
  pCosAmp[5][1] = CosFit->GetParameter(0);
  pCosAmpErr[5][1] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_445MeVCM6 = Phip_445MeVCM6_Para->GetAsymmetry(Phip_445MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_445MeVCM6->SetName("ParaPerpAsymmPhip445MeVCM6");
  ParaPerpAsymmPhip_445MeVCM6->SetTitle("Proton Para/Perp Phi Asymmetry for 440-450MeV Photon Energy (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_445MeVCM6->Fit("CosFit", "Q");
  pCosAmp[5][2] = CosFit->GetParameter(0);
  pCosAmpErr[5][2] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_455MeVCM6 = Phip_455MeVCM6_Para->GetAsymmetry(Phip_455MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_455MeVCM6->SetName("ParaPerpAsymmPhip455MeVCM6");
  ParaPerpAsymmPhip_455MeVCM6->SetTitle("Proton Para/Perp Phi Asymmetry for 450-460MeV Photon Energy (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_455MeVCM6->Fit("CosFit", "Q");
  pCosAmp[5][3] = CosFit->GetParameter(0);
  pCosAmpErr[5][3] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_465MeVCM6 = Phip_465MeVCM6_Para->GetAsymmetry(Phip_465MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_465MeVCM6->SetName("ParaPerpAsymmPhip465MeVCM6");
  ParaPerpAsymmPhip_465MeVCM6->SetTitle("Proton Para/Perp Phi Asymmetry for 460-470MeV Photon Energy (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_465MeVCM6->Fit("CosFit", "Q");
  pCosAmp[5][4] = CosFit->GetParameter(0);
  pCosAmpErr[5][4] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_475MeVCM6 = Phip_475MeVCM6_Para->GetAsymmetry(Phip_475MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_475MeVCM6->SetName("ParaPerpAsymmPhip475MeVCM6");
  ParaPerpAsymmPhip_475MeVCM6->SetTitle("Proton Para/Perp Phi Asymmetry for 470-480MeV Photon Energy (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_475MeVCM6->Fit("CosFit", "Q");
  pCosAmp[5][5] = CosFit->GetParameter(0);
  pCosAmpErr[5][5] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_485MeVCM6 = Phip_485MeVCM6_Para->GetAsymmetry(Phip_485MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_485MeVCM6->SetName("ParaPerpAsymmPhip485MeVCM6");
  ParaPerpAsymmPhip_485MeVCM6->SetTitle("Proton Para/Perp Phi Asymmetry for 480-490MeV Photon Energy (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_485MeVCM6->Fit("CosFit", "Q");
  pCosAmp[5][6] = CosFit->GetParameter(0);
  pCosAmpErr[5][6] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_495MeVCM6 = Phip_495MeVCM6_Para->GetAsymmetry(Phip_495MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_495MeVCM6->SetName("ParaPerpAsymmPhip495MeVCM6");
  ParaPerpAsymmPhip_495MeVCM6->SetTitle("Proton Para/Perp Phi Asymmetry for 490-500MeV Photon Energy (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_495MeVCM6->Fit("CosFit", "Q");
  pCosAmp[5][7] = CosFit->GetParameter(0);
  pCosAmpErr[5][7] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_505MeVCM6 = Phip_505MeVCM6_Para->GetAsymmetry(Phip_505MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_505MeVCM6->SetName("ParaPerpAsymmPhip505MeVCM6");
  ParaPerpAsymmPhip_505MeVCM6->SetTitle("Proton Para/Perp Phi Asymmetry for 500-510MeV Photon Energy (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_505MeVCM6->Fit("CosFit", "Q");
  pCosAmp[5][8] = CosFit->GetParameter(0);
  pCosAmpErr[5][8] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_515MeVCM6 = Phip_515MeVCM6_Para->GetAsymmetry(Phip_515MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_515MeVCM6->SetName("ParaPerpAsymmPhip515MeVCM6");
  ParaPerpAsymmPhip_515MeVCM6->SetTitle("Proton Para/Perp Phi Asymmetry for 510-520MeV Photon Energy (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_515MeVCM6->Fit("CosFit", "Q");
  pCosAmp[5][9] = CosFit->GetParameter(0);
  pCosAmpErr[5][9] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_525MeVCM6 = Phip_525MeVCM6_Para->GetAsymmetry(Phip_525MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_525MeVCM6->SetName("ParaPerpAsymmPhip525MeVCM6");
  ParaPerpAsymmPhip_525MeVCM6->SetTitle("Proton Para/Perp Phi Asymmetry for 520-530MeV Photon Energy (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_525MeVCM6->Fit("CosFit", "Q");
  pCosAmp[5][10] = CosFit->GetParameter(0);
  pCosAmpErr[5][10] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_535MeVCM6 = Phip_535MeVCM6_Para->GetAsymmetry(Phip_535MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_535MeVCM6->SetName("ParaPerpAsymmPhip535MeVCM6");
  ParaPerpAsymmPhip_535MeVCM6->SetTitle("Proton Para/Perp Phi Asymmetry for 530-540MeV Photon Energy (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_535MeVCM6->Fit("CosFit", "Q");
  pCosAmp[5][11] = CosFit->GetParameter(0);
  pCosAmpErr[5][11] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_545MeVCM6 = Phip_545MeVCM6_Para->GetAsymmetry(Phip_545MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_545MeVCM6->SetName("ParaPerpAsymmPhip545MeVCM6");
  ParaPerpAsymmPhip_545MeVCM6->SetTitle("Proton Para/Perp Phi Asymmetry for 540-550MeV Photon Energy (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_545MeVCM6->Fit("CosFit", "Q");
  pCosAmp[5][12] = CosFit->GetParameter(0);
  pCosAmpErr[5][12] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_555MeVCM6 = Phip_555MeVCM6_Para->GetAsymmetry(Phip_555MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_555MeVCM6->SetName("ParaPerpAsymmPhip555MeVCM6");
  ParaPerpAsymmPhip_555MeVCM6->SetTitle("Proton Para/Perp Phi Asymmetry for 550-560MeV Photon Energy (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_555MeVCM6->Fit("CosFit", "Q");
  pCosAmp[5][13] = CosFit->GetParameter(0);
  pCosAmpErr[5][13] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_565MeVCM6 = Phip_565MeVCM6_Para->GetAsymmetry(Phip_565MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_565MeVCM6->SetName("ParaPerpAsymmPhip565MeVCM6");
  ParaPerpAsymmPhip_565MeVCM6->SetTitle("Proton Para/Perp Phi Asymmetry for 560-570MeV Photon Energy (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_565MeVCM6->Fit("CosFit", "Q");
  pCosAmp[5][14] = CosFit->GetParameter(0);
  pCosAmpErr[5][14] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_575MeVCM6 = Phip_575MeVCM6_Para->GetAsymmetry(Phip_575MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_575MeVCM6->SetName("ParaPerpAsymmPhip575MeVCM6");
  ParaPerpAsymmPhip_575MeVCM6->SetTitle("Proton Para/Perp Phi Asymmetry for 570-580MeV Photon Energy (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_575MeVCM6->Fit("CosFit", "Q");
  pCosAmp[5][15] = CosFit->GetParameter(0);
  pCosAmpErr[5][15] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_585MeVCM6 = Phip_585MeVCM6_Para->GetAsymmetry(Phip_585MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_585MeVCM6->SetName("ParaPerpAsymmPhip585MeVCM6");
  ParaPerpAsymmPhip_585MeVCM6->SetTitle("Proton Para/Perp Phi Asymmetry for 580-590MeV Photon Energy (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_585MeVCM6->Fit("CosFit", "Q");
  pCosAmp[5][16] = CosFit->GetParameter(0);
  pCosAmpErr[5][16] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_595MeVCM6 = Phip_595MeVCM6_Para->GetAsymmetry(Phip_595MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_595MeVCM6->SetName("ParaPerpAsymmPhip595MeVCM6");
  ParaPerpAsymmPhip_595MeVCM6->SetTitle("Proton Para/Perp Phi Asymmetry for 590-600MeV Photon Energy (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_595MeVCM6->Fit("CosFit", "Q");
  pCosAmp[5][17] = CosFit->GetParameter(0);
  pCosAmpErr[5][17] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_605MeVCM6 = Phip_605MeVCM6_Para->GetAsymmetry(Phip_605MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_605MeVCM6->SetName("ParaPerpAsymmPhip605MeVCM6");
  ParaPerpAsymmPhip_605MeVCM6->SetTitle("Proton Para/Perp Phi Asymmetry for 600-610MeV Photon Energy (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_605MeVCM6->Fit("CosFit", "Q");
  pCosAmp[5][18] = CosFit->GetParameter(0);
  pCosAmpErr[5][18] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_615MeVCM6 = Phip_615MeVCM6_Para->GetAsymmetry(Phip_615MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_615MeVCM6->SetName("ParaPerpAsymmPhip615MeVCM6");
  ParaPerpAsymmPhip_615MeVCM6->SetTitle("Proton Para/Perp Phi Asymmetry for 610-620MeV Photon Energy (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_615MeVCM6->Fit("CosFit", "Q");
  pCosAmp[5][19] = CosFit->GetParameter(0);
  pCosAmpErr[5][19] = CosFit->GetParError(0);

  ///////////////////////////////////////////
  //////////////////CM7 p ///////////////////
  ///////////////////////////////////////////

  ParaPerpAsymmPhip_425MeVCM7 = Phip_425MeVCM7_Para->GetAsymmetry(Phip_425MeVCM7_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_425MeVCM7->SetName("ParaPerpAsymmPhip425MeVCM7");
  ParaPerpAsymmPhip_425MeVCM7->SetTitle("Proton Para/Perp Phi Asymmetry for 420-430MeV Photon Energy (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_425MeVCM7->Fit("CosFit", "Q");
  pCosAmp[6][0] = CosFit->GetParameter(0);
  pCosAmpErr[6][0] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_435MeVCM7 = Phip_435MeVCM7_Para->GetAsymmetry(Phip_435MeVCM7_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_435MeVCM7->SetName("ParaPerpAsymmPhip435MeVCM7");
  ParaPerpAsymmPhip_435MeVCM7->SetTitle("Proton Para/Perp Phi Asymmetry for 430-440MeV Photon Energy (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_435MeVCM7->Fit("CosFit", "Q");
  pCosAmp[6][1] = CosFit->GetParameter(0);
  pCosAmpErr[6][1] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_445MeVCM7 = Phip_445MeVCM7_Para->GetAsymmetry(Phip_445MeVCM7_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_445MeVCM7->SetName("ParaPerpAsymmPhip445MeVCM7");
  ParaPerpAsymmPhip_445MeVCM7->SetTitle("Proton Para/Perp Phi Asymmetry for 440-450MeV Photon Energy (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_445MeVCM7->Fit("CosFit", "Q");
  pCosAmp[6][2] = CosFit->GetParameter(0);
  pCosAmpErr[6][2] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_455MeVCM7 = Phip_455MeVCM7_Para->GetAsymmetry(Phip_455MeVCM7_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_455MeVCM7->SetName("ParaPerpAsymmPhip455MeVCM7");
  ParaPerpAsymmPhip_455MeVCM7->SetTitle("Proton Para/Perp Phi Asymmetry for 450-460MeV Photon Energy (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_455MeVCM7->Fit("CosFit", "Q");
  pCosAmp[6][3] = CosFit->GetParameter(0);
  pCosAmpErr[6][3] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_465MeVCM7 = Phip_465MeVCM7_Para->GetAsymmetry(Phip_465MeVCM7_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_465MeVCM7->SetName("ParaPerpAsymmPhip465MeVCM7");
  ParaPerpAsymmPhip_465MeVCM7->SetTitle("Proton Para/Perp Phi Asymmetry for 460-470MeV Photon Energy (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_465MeVCM7->Fit("CosFit", "Q");
  pCosAmp[6][4] = CosFit->GetParameter(0);
  pCosAmpErr[6][4] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_475MeVCM7 = Phip_475MeVCM7_Para->GetAsymmetry(Phip_475MeVCM7_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_475MeVCM7->SetName("ParaPerpAsymmPhip475MeVCM7");
  ParaPerpAsymmPhip_475MeVCM7->SetTitle("Proton Para/Perp Phi Asymmetry for 470-480MeV Photon Energy (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_475MeVCM7->Fit("CosFit", "Q");
  pCosAmp[6][5] = CosFit->GetParameter(0);
  pCosAmpErr[6][5] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_485MeVCM7 = Phip_485MeVCM7_Para->GetAsymmetry(Phip_485MeVCM7_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_485MeVCM7->SetName("ParaPerpAsymmPhip485MeVCM7");
  ParaPerpAsymmPhip_485MeVCM7->SetTitle("Proton Para/Perp Phi Asymmetry for 480-490MeV Photon Energy (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_485MeVCM7->Fit("CosFit", "Q");
  pCosAmp[6][6] = CosFit->GetParameter(0);
  pCosAmpErr[6][6] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_495MeVCM7 = Phip_495MeVCM7_Para->GetAsymmetry(Phip_495MeVCM7_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_495MeVCM7->SetName("ParaPerpAsymmPhip495MeVCM7");
  ParaPerpAsymmPhip_495MeVCM7->SetTitle("Proton Para/Perp Phi Asymmetry for 490-500MeV Photon Energy (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_495MeVCM7->Fit("CosFit", "Q");
  pCosAmp[6][7] = CosFit->GetParameter(0);
  pCosAmpErr[6][7] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_505MeVCM7 = Phip_505MeVCM7_Para->GetAsymmetry(Phip_505MeVCM7_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_505MeVCM7->SetName("ParaPerpAsymmPhip505MeVCM7");
  ParaPerpAsymmPhip_505MeVCM7->SetTitle("Proton Para/Perp Phi Asymmetry for 500-510MeV Photon Energy (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_505MeVCM7->Fit("CosFit", "Q");
  pCosAmp[6][8] = CosFit->GetParameter(0);
  pCosAmpErr[6][8] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_515MeVCM7 = Phip_515MeVCM7_Para->GetAsymmetry(Phip_515MeVCM7_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_515MeVCM7->SetName("ParaPerpAsymmPhip515MeVCM7");
  ParaPerpAsymmPhip_515MeVCM7->SetTitle("Proton Para/Perp Phi Asymmetry for 510-520MeV Photon Energy (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_515MeVCM7->Fit("CosFit", "Q");
  pCosAmp[6][9] = CosFit->GetParameter(0);
  pCosAmpErr[6][9] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_525MeVCM7 = Phip_525MeVCM7_Para->GetAsymmetry(Phip_525MeVCM7_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_525MeVCM7->SetName("ParaPerpAsymmPhip525MeVCM7");
  ParaPerpAsymmPhip_525MeVCM7->SetTitle("Proton Para/Perp Phi Asymmetry for 520-530MeV Photon Energy (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_525MeVCM7->Fit("CosFit", "Q");
  pCosAmp[6][10] = CosFit->GetParameter(0);
  pCosAmpErr[6][10] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_535MeVCM7 = Phip_535MeVCM7_Para->GetAsymmetry(Phip_535MeVCM7_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_535MeVCM7->SetName("ParaPerpAsymmPhip535MeVCM7");
  ParaPerpAsymmPhip_535MeVCM7->SetTitle("Proton Para/Perp Phi Asymmetry for 530-540MeV Photon Energy (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_535MeVCM7->Fit("CosFit", "Q");
  pCosAmp[6][11] = CosFit->GetParameter(0);
  pCosAmpErr[6][11] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_545MeVCM7 = Phip_545MeVCM7_Para->GetAsymmetry(Phip_545MeVCM7_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_545MeVCM7->SetName("ParaPerpAsymmPhip545MeVCM7");
  ParaPerpAsymmPhip_545MeVCM7->SetTitle("Proton Para/Perp Phi Asymmetry for 540-550MeV Photon Energy (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_545MeVCM7->Fit("CosFit", "Q");
  pCosAmp[6][12] = CosFit->GetParameter(0);
  pCosAmpErr[6][12] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_555MeVCM7 = Phip_555MeVCM7_Para->GetAsymmetry(Phip_555MeVCM7_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_555MeVCM7->SetName("ParaPerpAsymmPhip555MeVCM7");
  ParaPerpAsymmPhip_555MeVCM7->SetTitle("Proton Para/Perp Phi Asymmetry for 550-560MeV Photon Energy (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_555MeVCM7->Fit("CosFit", "Q");
  pCosAmp[6][13] = CosFit->GetParameter(0);
  pCosAmpErr[6][13] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_565MeVCM7 = Phip_565MeVCM7_Para->GetAsymmetry(Phip_565MeVCM7_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_565MeVCM7->SetName("ParaPerpAsymmPhip565MeVCM7");
  ParaPerpAsymmPhip_565MeVCM7->SetTitle("Proton Para/Perp Phi Asymmetry for 560-570MeV Photon Energy (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_565MeVCM7->Fit("CosFit", "Q");
  pCosAmp[6][14] = CosFit->GetParameter(0);
  pCosAmpErr[6][14] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_575MeVCM7 = Phip_575MeVCM7_Para->GetAsymmetry(Phip_575MeVCM7_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_575MeVCM7->SetName("ParaPerpAsymmPhip575MeVCM7");
  ParaPerpAsymmPhip_575MeVCM7->SetTitle("Proton Para/Perp Phi Asymmetry for 570-580MeV Photon Energy (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_575MeVCM7->Fit("CosFit", "Q");
  pCosAmp[6][15] = CosFit->GetParameter(0);
  pCosAmpErr[6][15] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_585MeVCM7 = Phip_585MeVCM7_Para->GetAsymmetry(Phip_585MeVCM7_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_585MeVCM7->SetName("ParaPerpAsymmPhip585MeVCM7");
  ParaPerpAsymmPhip_585MeVCM7->SetTitle("Proton Para/Perp Phi Asymmetry for 580-590MeV Photon Energy (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_585MeVCM7->Fit("CosFit", "Q");
  pCosAmp[6][16] = CosFit->GetParameter(0);
  pCosAmpErr[6][16] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_595MeVCM7 = Phip_595MeVCM7_Para->GetAsymmetry(Phip_595MeVCM7_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_595MeVCM7->SetName("ParaPerpAsymmPhip595MeVCM7");
  ParaPerpAsymmPhip_595MeVCM7->SetTitle("Proton Para/Perp Phi Asymmetry for 590-600MeV Photon Energy (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_595MeVCM7->Fit("CosFit", "Q");
  pCosAmp[6][17] = CosFit->GetParameter(0);
  pCosAmpErr[6][17] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_605MeVCM7 = Phip_605MeVCM7_Para->GetAsymmetry(Phip_605MeVCM7_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_605MeVCM7->SetName("ParaPerpAsymmPhip605MeVCM7");
  ParaPerpAsymmPhip_605MeVCM7->SetTitle("Proton Para/Perp Phi Asymmetry for 600-610MeV Photon Energy (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_605MeVCM7->Fit("CosFit", "Q");
  pCosAmp[6][18] = CosFit->GetParameter(0);
  pCosAmpErr[6][18] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_615MeVCM7 = Phip_615MeVCM7_Para->GetAsymmetry(Phip_615MeVCM7_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_615MeVCM7->SetName("ParaPerpAsymmPhip615MeVCM7");
  ParaPerpAsymmPhip_615MeVCM7->SetTitle("Proton Para/Perp Phi Asymmetry for 610-620MeV Photon Energy (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_615MeVCM7->Fit("CosFit", "Q");
  pCosAmp[6][19] = CosFit->GetParameter(0);
  pCosAmpErr[6][19] = CosFit->GetParError(0);

  ///////////////////////////////////////////
  //////////////////CM8 p ///////////////////
  ///////////////////////////////////////////

  ParaPerpAsymmPhip_425MeVCM8 = Phip_425MeVCM8_Para->GetAsymmetry(Phip_425MeVCM8_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_425MeVCM8->SetName("ParaPerpAsymmPhip425MeVCM8");
  ParaPerpAsymmPhip_425MeVCM8->SetTitle("Proton Para/Perp Phi Asymmetry for 420-430MeV Photon Energy (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_425MeVCM8->Fit("CosFit", "Q");
  pCosAmp[7][0] = CosFit->GetParameter(0);
  pCosAmpErr[7][0] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_435MeVCM8 = Phip_435MeVCM8_Para->GetAsymmetry(Phip_435MeVCM8_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_435MeVCM8->SetName("ParaPerpAsymmPhip435MeVCM8");
  ParaPerpAsymmPhip_435MeVCM8->SetTitle("Proton Para/Perp Phi Asymmetry for 430-440MeV Photon Energy (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_435MeVCM8->Fit("CosFit", "Q");
  pCosAmp[7][1] = CosFit->GetParameter(0);
  pCosAmpErr[7][1] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_445MeVCM8 = Phip_445MeVCM8_Para->GetAsymmetry(Phip_445MeVCM8_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_445MeVCM8->SetName("ParaPerpAsymmPhip445MeVCM8");
  ParaPerpAsymmPhip_445MeVCM8->SetTitle("Proton Para/Perp Phi Asymmetry for 440-450MeV Photon Energy (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_445MeVCM8->Fit("CosFit", "Q");
  pCosAmp[7][2] = CosFit->GetParameter(0);
  pCosAmpErr[7][2] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_455MeVCM8 = Phip_455MeVCM8_Para->GetAsymmetry(Phip_455MeVCM8_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_455MeVCM8->SetName("ParaPerpAsymmPhip455MeVCM8");
  ParaPerpAsymmPhip_455MeVCM8->SetTitle("Proton Para/Perp Phi Asymmetry for 450-460MeV Photon Energy (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_455MeVCM8->Fit("CosFit", "Q");
  pCosAmp[7][3] = CosFit->GetParameter(0);
  pCosAmpErr[7][3] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_465MeVCM8 = Phip_465MeVCM8_Para->GetAsymmetry(Phip_465MeVCM8_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_465MeVCM8->SetName("ParaPerpAsymmPhip465MeVCM8");
  ParaPerpAsymmPhip_465MeVCM8->SetTitle("Proton Para/Perp Phi Asymmetry for 460-470MeV Photon Energy (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_465MeVCM8->Fit("CosFit", "Q");
  pCosAmp[7][4] = CosFit->GetParameter(0);
  pCosAmpErr[7][4] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_475MeVCM8 = Phip_475MeVCM8_Para->GetAsymmetry(Phip_475MeVCM8_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_475MeVCM8->SetName("ParaPerpAsymmPhip475MeVCM8");
  ParaPerpAsymmPhip_475MeVCM8->SetTitle("Proton Para/Perp Phi Asymmetry for 470-480MeV Photon Energy (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_475MeVCM8->Fit("CosFit", "Q");
  pCosAmp[7][5] = CosFit->GetParameter(0);
  pCosAmpErr[7][5] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_485MeVCM8 = Phip_485MeVCM8_Para->GetAsymmetry(Phip_485MeVCM8_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_485MeVCM8->SetName("ParaPerpAsymmPhip485MeVCM8");
  ParaPerpAsymmPhip_485MeVCM8->SetTitle("Proton Para/Perp Phi Asymmetry for 480-490MeV Photon Energy (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_485MeVCM8->Fit("CosFit", "Q");
  pCosAmp[7][6] = CosFit->GetParameter(0);
  pCosAmpErr[7][6] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_495MeVCM8 = Phip_495MeVCM8_Para->GetAsymmetry(Phip_495MeVCM8_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_495MeVCM8->SetName("ParaPerpAsymmPhip495MeVCM8");
  ParaPerpAsymmPhip_495MeVCM8->SetTitle("Proton Para/Perp Phi Asymmetry for 490-500MeV Photon Energy (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_495MeVCM8->Fit("CosFit", "Q");
  pCosAmp[7][7] = CosFit->GetParameter(0);
  pCosAmpErr[7][7] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_505MeVCM8 = Phip_505MeVCM8_Para->GetAsymmetry(Phip_505MeVCM8_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_505MeVCM8->SetName("ParaPerpAsymmPhip505MeVCM8");
  ParaPerpAsymmPhip_505MeVCM8->SetTitle("Proton Para/Perp Phi Asymmetry for 500-510MeV Photon Energy (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_505MeVCM8->Fit("CosFit", "Q");
  pCosAmp[7][8] = CosFit->GetParameter(0);
  pCosAmpErr[7][8] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_515MeVCM8 = Phip_515MeVCM8_Para->GetAsymmetry(Phip_515MeVCM8_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_515MeVCM8->SetName("ParaPerpAsymmPhip515MeVCM8");
  ParaPerpAsymmPhip_515MeVCM8->SetTitle("Proton Para/Perp Phi Asymmetry for 510-520MeV Photon Energy (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_515MeVCM8->Fit("CosFit", "Q");
  pCosAmp[7][9] = CosFit->GetParameter(0);
  pCosAmpErr[7][9] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_525MeVCM8 = Phip_525MeVCM8_Para->GetAsymmetry(Phip_525MeVCM8_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_525MeVCM8->SetName("ParaPerpAsymmPhip525MeVCM8");
  ParaPerpAsymmPhip_525MeVCM8->SetTitle("Proton Para/Perp Phi Asymmetry for 520-530MeV Photon Energy (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_525MeVCM8->Fit("CosFit", "Q");
  pCosAmp[7][10] = CosFit->GetParameter(0);
  pCosAmpErr[7][10] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_535MeVCM8 = Phip_535MeVCM8_Para->GetAsymmetry(Phip_535MeVCM8_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_535MeVCM8->SetName("ParaPerpAsymmPhip535MeVCM8");
  ParaPerpAsymmPhip_535MeVCM8->SetTitle("Proton Para/Perp Phi Asymmetry for 530-540MeV Photon Energy (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_535MeVCM8->Fit("CosFit", "Q");
  pCosAmp[7][11] = CosFit->GetParameter(0);
  pCosAmpErr[7][11] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_545MeVCM8 = Phip_545MeVCM8_Para->GetAsymmetry(Phip_545MeVCM8_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_545MeVCM8->SetName("ParaPerpAsymmPhip545MeVCM8");
  ParaPerpAsymmPhip_545MeVCM8->SetTitle("Proton Para/Perp Phi Asymmetry for 540-550MeV Photon Energy (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_545MeVCM8->Fit("CosFit", "Q");
  pCosAmp[7][12] = CosFit->GetParameter(0);
  pCosAmpErr[7][12] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_555MeVCM8 = Phip_555MeVCM8_Para->GetAsymmetry(Phip_555MeVCM8_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_555MeVCM8->SetName("ParaPerpAsymmPhip555MeVCM8");
  ParaPerpAsymmPhip_555MeVCM8->SetTitle("Proton Para/Perp Phi Asymmetry for 550-560MeV Photon Energy (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_555MeVCM8->Fit("CosFit", "Q");
  pCosAmp[7][13] = CosFit->GetParameter(0);
  pCosAmpErr[7][13] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_565MeVCM8 = Phip_565MeVCM8_Para->GetAsymmetry(Phip_565MeVCM8_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_565MeVCM8->SetName("ParaPerpAsymmPhip565MeVCM8");
  ParaPerpAsymmPhip_565MeVCM8->SetTitle("Proton Para/Perp Phi Asymmetry for 560-570MeV Photon Energy (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_565MeVCM8->Fit("CosFit", "Q");
  pCosAmp[7][14] = CosFit->GetParameter(0);
  pCosAmpErr[7][14] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_575MeVCM8 = Phip_575MeVCM8_Para->GetAsymmetry(Phip_575MeVCM8_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_575MeVCM8->SetName("ParaPerpAsymmPhip575MeVCM8");
  ParaPerpAsymmPhip_575MeVCM8->SetTitle("Proton Para/Perp Phi Asymmetry for 570-580MeV Photon Energy (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_575MeVCM8->Fit("CosFit", "Q");
  pCosAmp[7][15] = CosFit->GetParameter(0);
  pCosAmpErr[7][15] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_585MeVCM8 = Phip_585MeVCM8_Para->GetAsymmetry(Phip_585MeVCM8_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_585MeVCM8->SetName("ParaPerpAsymmPhip585MeVCM8");
  ParaPerpAsymmPhip_585MeVCM8->SetTitle("Proton Para/Perp Phi Asymmetry for 580-590MeV Photon Energy (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_585MeVCM8->Fit("CosFit", "Q");
  pCosAmp[7][16] = CosFit->GetParameter(0);
  pCosAmpErr[7][16] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_595MeVCM8 = Phip_595MeVCM8_Para->GetAsymmetry(Phip_595MeVCM8_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_595MeVCM8->SetName("ParaPerpAsymmPhip595MeVCM8");
  ParaPerpAsymmPhip_595MeVCM8->SetTitle("Proton Para/Perp Phi Asymmetry for 590-600MeV Photon Energy (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_595MeVCM8->Fit("CosFit", "Q");
  pCosAmp[7][17] = CosFit->GetParameter(0);
  pCosAmpErr[7][17] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_605MeVCM8 = Phip_605MeVCM8_Para->GetAsymmetry(Phip_605MeVCM8_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_605MeVCM8->SetName("ParaPerpAsymmPhip605MeVCM8");
  ParaPerpAsymmPhip_605MeVCM8->SetTitle("Proton Para/Perp Phi Asymmetry for 600-610MeV Photon Energy (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_605MeVCM8->Fit("CosFit", "Q");
  pCosAmp[7][18] = CosFit->GetParameter(0);
  pCosAmpErr[7][18] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_615MeVCM8 = Phip_615MeVCM8_Para->GetAsymmetry(Phip_615MeVCM8_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_615MeVCM8->SetName("ParaPerpAsymmPhip615MeVCM8");
  ParaPerpAsymmPhip_615MeVCM8->SetTitle("Proton Para/Perp Phi Asymmetry for 610-620MeV Photon Energy (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_615MeVCM8->Fit("CosFit", "Q");
  pCosAmp[7][19] = CosFit->GetParameter(0);
  pCosAmpErr[7][19] = CosFit->GetParError(0);

  ///////////////////////////////////////////
  //////////////////CM9 p ///////////////////
  ///////////////////////////////////////////

  ParaPerpAsymmPhip_425MeVCM9 = Phip_425MeVCM9_Para->GetAsymmetry(Phip_425MeVCM9_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_425MeVCM9->SetName("ParaPerpAsymmPhip425MeVCM9");
  ParaPerpAsymmPhip_425MeVCM9->SetTitle("Proton Para/Perp Phi Asymmetry for 420-430MeV Photon Energy (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_425MeVCM9->Fit("CosFit", "Q");
  pCosAmp[8][0] = CosFit->GetParameter(0);
  pCosAmpErr[8][0] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_435MeVCM9 = Phip_435MeVCM9_Para->GetAsymmetry(Phip_435MeVCM9_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_435MeVCM9->SetName("ParaPerpAsymmPhip435MeVCM9");
  ParaPerpAsymmPhip_435MeVCM9->SetTitle("Proton Para/Perp Phi Asymmetry for 430-440MeV Photon Energy (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_435MeVCM9->Fit("CosFit", "Q");
  pCosAmp[8][1] = CosFit->GetParameter(0);
  pCosAmpErr[8][1] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_445MeVCM9 = Phip_445MeVCM9_Para->GetAsymmetry(Phip_445MeVCM9_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_445MeVCM9->SetName("ParaPerpAsymmPhip445MeVCM9");
  ParaPerpAsymmPhip_445MeVCM9->SetTitle("Proton Para/Perp Phi Asymmetry for 440-450MeV Photon Energy (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_445MeVCM9->Fit("CosFit", "Q");
  pCosAmp[8][2] = CosFit->GetParameter(0);
  pCosAmpErr[8][2] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_455MeVCM9 = Phip_455MeVCM9_Para->GetAsymmetry(Phip_455MeVCM9_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_455MeVCM9->SetName("ParaPerpAsymmPhip455MeVCM9");
  ParaPerpAsymmPhip_455MeVCM9->SetTitle("Proton Para/Perp Phi Asymmetry for 450-460MeV Photon Energy (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_455MeVCM9->Fit("CosFit", "Q");
  pCosAmp[8][3] = CosFit->GetParameter(0);
  pCosAmpErr[8][3] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_465MeVCM9 = Phip_465MeVCM9_Para->GetAsymmetry(Phip_465MeVCM9_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_465MeVCM9->SetName("ParaPerpAsymmPhip465MeVCM9");
  ParaPerpAsymmPhip_465MeVCM9->SetTitle("Proton Para/Perp Phi Asymmetry for 460-470MeV Photon Energy (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_465MeVCM9->Fit("CosFit", "Q");
  pCosAmp[8][4] = CosFit->GetParameter(0);
  pCosAmpErr[8][4] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_475MeVCM9 = Phip_475MeVCM9_Para->GetAsymmetry(Phip_475MeVCM9_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_475MeVCM9->SetName("ParaPerpAsymmPhip475MeVCM9");
  ParaPerpAsymmPhip_475MeVCM9->SetTitle("Proton Para/Perp Phi Asymmetry for 470-480MeV Photon Energy (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_475MeVCM9->Fit("CosFit", "Q");
  pCosAmp[8][5] = CosFit->GetParameter(0);
  pCosAmpErr[8][5] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_485MeVCM9 = Phip_485MeVCM9_Para->GetAsymmetry(Phip_485MeVCM9_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_485MeVCM9->SetName("ParaPerpAsymmPhip485MeVCM9");
  ParaPerpAsymmPhip_485MeVCM9->SetTitle("Proton Para/Perp Phi Asymmetry for 480-490MeV Photon Energy (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_485MeVCM9->Fit("CosFit", "Q");
  pCosAmp[8][6] = CosFit->GetParameter(0);
  pCosAmpErr[8][6] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_495MeVCM9 = Phip_495MeVCM9_Para->GetAsymmetry(Phip_495MeVCM9_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_495MeVCM9->SetName("ParaPerpAsymmPhip495MeVCM9");
  ParaPerpAsymmPhip_495MeVCM9->SetTitle("Proton Para/Perp Phi Asymmetry for 490-500MeV Photon Energy (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_495MeVCM9->Fit("CosFit", "Q");
  pCosAmp[8][7] = CosFit->GetParameter(0);
  pCosAmpErr[8][7] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_505MeVCM9 = Phip_505MeVCM9_Para->GetAsymmetry(Phip_505MeVCM9_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_505MeVCM9->SetName("ParaPerpAsymmPhip505MeVCM9");
  ParaPerpAsymmPhip_505MeVCM9->SetTitle("Proton Para/Perp Phi Asymmetry for 500-510MeV Photon Energy (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_505MeVCM9->Fit("CosFit", "Q");
  pCosAmp[8][8] = CosFit->GetParameter(0);
  pCosAmpErr[8][8] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_515MeVCM9 = Phip_515MeVCM9_Para->GetAsymmetry(Phip_515MeVCM9_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_515MeVCM9->SetName("ParaPerpAsymmPhip515MeVCM9");
  ParaPerpAsymmPhip_515MeVCM9->SetTitle("Proton Para/Perp Phi Asymmetry for 510-520MeV Photon Energy (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_515MeVCM9->Fit("CosFit", "Q");
  pCosAmp[8][9] = CosFit->GetParameter(0);
  pCosAmpErr[8][9] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_525MeVCM9 = Phip_525MeVCM9_Para->GetAsymmetry(Phip_525MeVCM9_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_525MeVCM9->SetName("ParaPerpAsymmPhip525MeVCM9");
  ParaPerpAsymmPhip_525MeVCM9->SetTitle("Proton Para/Perp Phi Asymmetry for 520-530MeV Photon Energy (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_525MeVCM9->Fit("CosFit", "Q");
  pCosAmp[8][10] = CosFit->GetParameter(0);
  pCosAmpErr[8][10] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_535MeVCM9 = Phip_535MeVCM9_Para->GetAsymmetry(Phip_535MeVCM9_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_535MeVCM9->SetName("ParaPerpAsymmPhip535MeVCM9");
  ParaPerpAsymmPhip_535MeVCM9->SetTitle("Proton Para/Perp Phi Asymmetry for 530-540MeV Photon Energy (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_535MeVCM9->Fit("CosFit", "Q");
  pCosAmp[8][11] = CosFit->GetParameter(0);
  pCosAmpErr[8][11] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_545MeVCM9 = Phip_545MeVCM9_Para->GetAsymmetry(Phip_545MeVCM9_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_545MeVCM9->SetName("ParaPerpAsymmPhip545MeVCM9");
  ParaPerpAsymmPhip_545MeVCM9->SetTitle("Proton Para/Perp Phi Asymmetry for 540-550MeV Photon Energy (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_545MeVCM9->Fit("CosFit", "Q");
  pCosAmp[8][12] = CosFit->GetParameter(0);
  pCosAmpErr[8][12] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_555MeVCM9 = Phip_555MeVCM9_Para->GetAsymmetry(Phip_555MeVCM9_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_555MeVCM9->SetName("ParaPerpAsymmPhip555MeVCM9");
  ParaPerpAsymmPhip_555MeVCM9->SetTitle("Proton Para/Perp Phi Asymmetry for 550-560MeV Photon Energy (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_555MeVCM9->Fit("CosFit", "Q");
  pCosAmp[8][13] = CosFit->GetParameter(0);
  pCosAmpErr[8][13] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_565MeVCM9 = Phip_565MeVCM9_Para->GetAsymmetry(Phip_565MeVCM9_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_565MeVCM9->SetName("ParaPerpAsymmPhip565MeVCM9");
  ParaPerpAsymmPhip_565MeVCM9->SetTitle("Proton Para/Perp Phi Asymmetry for 560-570MeV Photon Energy (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_565MeVCM9->Fit("CosFit", "Q");
  pCosAmp[8][14] = CosFit->GetParameter(0);
  pCosAmpErr[8][14] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_575MeVCM9 = Phip_575MeVCM9_Para->GetAsymmetry(Phip_575MeVCM9_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_575MeVCM9->SetName("ParaPerpAsymmPhip575MeVCM9");
  ParaPerpAsymmPhip_575MeVCM9->SetTitle("Proton Para/Perp Phi Asymmetry for 570-580MeV Photon Energy (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_575MeVCM9->Fit("CosFit", "Q");
  pCosAmp[8][15] = CosFit->GetParameter(0);
  pCosAmpErr[8][15] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_585MeVCM9 = Phip_585MeVCM9_Para->GetAsymmetry(Phip_585MeVCM9_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_585MeVCM9->SetName("ParaPerpAsymmPhip585MeVCM9");
  ParaPerpAsymmPhip_585MeVCM9->SetTitle("Proton Para/Perp Phi Asymmetry for 580-590MeV Photon Energy (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_585MeVCM9->Fit("CosFit", "Q");
  pCosAmp[8][16] = CosFit->GetParameter(0);
  pCosAmpErr[8][16] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_595MeVCM9 = Phip_595MeVCM9_Para->GetAsymmetry(Phip_595MeVCM9_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_595MeVCM9->SetName("ParaPerpAsymmPhip595MeVCM9");
  ParaPerpAsymmPhip_595MeVCM9->SetTitle("Proton Para/Perp Phi Asymmetry for 590-600MeV Photon Energy (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_595MeVCM9->Fit("CosFit", "Q");
  pCosAmp[8][17] = CosFit->GetParameter(0);
  pCosAmpErr[8][17] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_605MeVCM9 = Phip_605MeVCM9_Para->GetAsymmetry(Phip_605MeVCM9_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_605MeVCM9->SetName("ParaPerpAsymmPhip605MeVCM9");
  ParaPerpAsymmPhip_605MeVCM9->SetTitle("Proton Para/Perp Phi Asymmetry for 600-610MeV Photon Energy (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_605MeVCM9->Fit("CosFit", "Q");
  pCosAmp[8][18] = CosFit->GetParameter(0);
  pCosAmpErr[8][18] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_615MeVCM9 = Phip_615MeVCM9_Para->GetAsymmetry(Phip_615MeVCM9_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_615MeVCM9->SetName("ParaPerpAsymmPhip615MeVCM9");
  ParaPerpAsymmPhip_615MeVCM9->SetTitle("Proton Para/Perp Phi Asymmetry for 610-620MeV Photon Energy (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_615MeVCM9->Fit("CosFit", "Q");
  pCosAmp[8][19] = CosFit->GetParameter(0);
  pCosAmpErr[8][19] = CosFit->GetParError(0);

  ///////////////////////////////////////////
  //////////////////CM10 p ///////////////////
  ///////////////////////////////////////////

  ParaPerpAsymmPhip_425MeVCM10 = Phip_425MeVCM10_Para->GetAsymmetry(Phip_425MeVCM10_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_425MeVCM10->SetName("ParaPerpAsymmPhip425MeVCM10");
  ParaPerpAsymmPhip_425MeVCM10->SetTitle("Proton Para/Perp Phi Asymmetry for 420-430MeV Photon Energy (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_425MeVCM10->Fit("CosFit", "Q");
  pCosAmp[9][0] = CosFit->GetParameter(0);
  pCosAmpErr[9][0] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_435MeVCM10 = Phip_435MeVCM10_Para->GetAsymmetry(Phip_435MeVCM10_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_435MeVCM10->SetName("ParaPerpAsymmPhip435MeVCM10");
  ParaPerpAsymmPhip_435MeVCM10->SetTitle("Proton Para/Perp Phi Asymmetry for 430-440MeV Photon Energy (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_435MeVCM10->Fit("CosFit", "Q");
  pCosAmp[9][1] = CosFit->GetParameter(0);
  pCosAmpErr[9][1] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_445MeVCM10 = Phip_445MeVCM10_Para->GetAsymmetry(Phip_445MeVCM10_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_445MeVCM10->SetName("ParaPerpAsymmPhip445MeVCM10");
  ParaPerpAsymmPhip_445MeVCM10->SetTitle("Proton Para/Perp Phi Asymmetry for 440-450MeV Photon Energy (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_445MeVCM10->Fit("CosFit", "Q");
  pCosAmp[9][2] = CosFit->GetParameter(0);
  pCosAmpErr[9][2] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_455MeVCM10 = Phip_455MeVCM10_Para->GetAsymmetry(Phip_455MeVCM10_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_455MeVCM10->SetName("ParaPerpAsymmPhip455MeVCM10");
  ParaPerpAsymmPhip_455MeVCM10->SetTitle("Proton Para/Perp Phi Asymmetry for 450-460MeV Photon Energy (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_455MeVCM10->Fit("CosFit", "Q");
  pCosAmp[9][3] = CosFit->GetParameter(0);
  pCosAmpErr[9][3] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_465MeVCM10 = Phip_465MeVCM10_Para->GetAsymmetry(Phip_465MeVCM10_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_465MeVCM10->SetName("ParaPerpAsymmPhip465MeVCM10");
  ParaPerpAsymmPhip_465MeVCM10->SetTitle("Proton Para/Perp Phi Asymmetry for 460-470MeV Photon Energy (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_465MeVCM10->Fit("CosFit", "Q");
  pCosAmp[9][4] = CosFit->GetParameter(0);
  pCosAmpErr[9][4] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_475MeVCM10 = Phip_475MeVCM10_Para->GetAsymmetry(Phip_475MeVCM10_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_475MeVCM10->SetName("ParaPerpAsymmPhip475MeVCM10");
  ParaPerpAsymmPhip_475MeVCM10->SetTitle("Proton Para/Perp Phi Asymmetry for 470-480MeV Photon Energy (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_475MeVCM10->Fit("CosFit", "Q");
  pCosAmp[9][5] = CosFit->GetParameter(0);
  pCosAmpErr[9][5] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_485MeVCM10 = Phip_485MeVCM10_Para->GetAsymmetry(Phip_485MeVCM10_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_485MeVCM10->SetName("ParaPerpAsymmPhip485MeVCM10");
  ParaPerpAsymmPhip_485MeVCM10->SetTitle("Proton Para/Perp Phi Asymmetry for 480-490MeV Photon Energy (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_485MeVCM10->Fit("CosFit", "Q");
  pCosAmp[9][6] = CosFit->GetParameter(0);
  pCosAmpErr[9][6] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_495MeVCM10 = Phip_495MeVCM10_Para->GetAsymmetry(Phip_495MeVCM10_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_495MeVCM10->SetName("ParaPerpAsymmPhip495MeVCM10");
  ParaPerpAsymmPhip_495MeVCM10->SetTitle("Proton Para/Perp Phi Asymmetry for 490-500MeV Photon Energy (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_495MeVCM10->Fit("CosFit", "Q");
  pCosAmp[9][7] = CosFit->GetParameter(0);
  pCosAmpErr[9][7] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_505MeVCM10 = Phip_505MeVCM10_Para->GetAsymmetry(Phip_505MeVCM10_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_505MeVCM10->SetName("ParaPerpAsymmPhip505MeVCM10");
  ParaPerpAsymmPhip_505MeVCM10->SetTitle("Proton Para/Perp Phi Asymmetry for 500-510MeV Photon Energy (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_505MeVCM10->Fit("CosFit", "Q");
  pCosAmp[9][8] = CosFit->GetParameter(0);
  pCosAmpErr[9][8] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_515MeVCM10 = Phip_515MeVCM10_Para->GetAsymmetry(Phip_515MeVCM10_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_515MeVCM10->SetName("ParaPerpAsymmPhip515MeVCM10");
  ParaPerpAsymmPhip_515MeVCM10->SetTitle("Proton Para/Perp Phi Asymmetry for 510-520MeV Photon Energy (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_515MeVCM10->Fit("CosFit", "Q");
  pCosAmp[9][9] = CosFit->GetParameter(0);
  pCosAmpErr[9][9] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_525MeVCM10 = Phip_525MeVCM10_Para->GetAsymmetry(Phip_525MeVCM10_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_525MeVCM10->SetName("ParaPerpAsymmPhip525MeVCM10");
  ParaPerpAsymmPhip_525MeVCM10->SetTitle("Proton Para/Perp Phi Asymmetry for 520-530MeV Photon Energy (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_525MeVCM10->Fit("CosFit", "Q");
  pCosAmp[9][10] = CosFit->GetParameter(0);
  pCosAmpErr[9][10] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_535MeVCM10 = Phip_535MeVCM10_Para->GetAsymmetry(Phip_535MeVCM10_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_535MeVCM10->SetName("ParaPerpAsymmPhip535MeVCM10");
  ParaPerpAsymmPhip_535MeVCM10->SetTitle("Proton Para/Perp Phi Asymmetry for 530-540MeV Photon Energy (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_535MeVCM10->Fit("CosFit", "Q");
  pCosAmp[9][11] = CosFit->GetParameter(0);
  pCosAmpErr[9][11] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_545MeVCM10 = Phip_545MeVCM10_Para->GetAsymmetry(Phip_545MeVCM10_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_545MeVCM10->SetName("ParaPerpAsymmPhip545MeVCM10");
  ParaPerpAsymmPhip_545MeVCM10->SetTitle("Proton Para/Perp Phi Asymmetry for 540-550MeV Photon Energy (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_545MeVCM10->Fit("CosFit", "Q");
  pCosAmp[9][12] = CosFit->GetParameter(0);
  pCosAmpErr[9][12] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_555MeVCM10 = Phip_555MeVCM10_Para->GetAsymmetry(Phip_555MeVCM10_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_555MeVCM10->SetName("ParaPerpAsymmPhip555MeVCM10");
  ParaPerpAsymmPhip_555MeVCM10->SetTitle("Proton Para/Perp Phi Asymmetry for 550-560MeV Photon Energy (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_555MeVCM10->Fit("CosFit", "Q");
  pCosAmp[9][13] = CosFit->GetParameter(0);
  pCosAmpErr[9][13] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_565MeVCM10 = Phip_565MeVCM10_Para->GetAsymmetry(Phip_565MeVCM10_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_565MeVCM10->SetName("ParaPerpAsymmPhip565MeVCM10");
  ParaPerpAsymmPhip_565MeVCM10->SetTitle("Proton Para/Perp Phi Asymmetry for 560-570MeV Photon Energy (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_565MeVCM10->Fit("CosFit", "Q");
  pCosAmp[9][14] = CosFit->GetParameter(0);
  pCosAmpErr[9][14] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_575MeVCM10 = Phip_575MeVCM10_Para->GetAsymmetry(Phip_575MeVCM10_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_575MeVCM10->SetName("ParaPerpAsymmPhip575MeVCM10");
  ParaPerpAsymmPhip_575MeVCM10->SetTitle("Proton Para/Perp Phi Asymmetry for 570-580MeV Photon Energy (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_575MeVCM10->Fit("CosFit", "Q");
  pCosAmp[9][15] = CosFit->GetParameter(0);
  pCosAmpErr[9][15] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_585MeVCM10 = Phip_585MeVCM10_Para->GetAsymmetry(Phip_585MeVCM10_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_585MeVCM10->SetName("ParaPerpAsymmPhip585MeVCM10");
  ParaPerpAsymmPhip_585MeVCM10->SetTitle("Proton Para/Perp Phi Asymmetry for 580-590MeV Photon Energy (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_585MeVCM10->Fit("CosFit", "Q");
  pCosAmp[9][16] = CosFit->GetParameter(0);
  pCosAmpErr[9][16] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_595MeVCM10 = Phip_595MeVCM10_Para->GetAsymmetry(Phip_595MeVCM10_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_595MeVCM10->SetName("ParaPerpAsymmPhip595MeVCM10");
  ParaPerpAsymmPhip_595MeVCM10->SetTitle("Proton Para/Perp Phi Asymmetry for 590-600MeV Photon Energy (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_595MeVCM10->Fit("CosFit", "Q");
  pCosAmp[9][17] = CosFit->GetParameter(0);
  pCosAmpErr[9][17] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_605MeVCM10 = Phip_605MeVCM10_Para->GetAsymmetry(Phip_605MeVCM10_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_605MeVCM10->SetName("ParaPerpAsymmPhip605MeVCM10");
  ParaPerpAsymmPhip_605MeVCM10->SetTitle("Proton Para/Perp Phi Asymmetry for 600-610MeV Photon Energy (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_605MeVCM10->Fit("CosFit", "Q");
  pCosAmp[9][18] = CosFit->GetParameter(0);
  pCosAmpErr[9][18] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_615MeVCM10 = Phip_615MeVCM10_Para->GetAsymmetry(Phip_615MeVCM10_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_615MeVCM10->SetName("ParaPerpAsymmPhip615MeVCM10");
  ParaPerpAsymmPhip_615MeVCM10->SetTitle("Proton Para/Perp Phi Asymmetry for 610-620MeV Photon Energy (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_615MeVCM10->Fit("CosFit", "Q");
  pCosAmp[9][19] = CosFit->GetParameter(0);
  pCosAmpErr[9][19] = CosFit->GetParError(0);

  TFile f1("ParaPerpAsymm_NoScatt_Total_7.root", "RECREATE");

  ParaPerpAsymmPhip_425MeVCM1->Write();
  ParaPerpAsymmPhip_435MeVCM1->Write();
  ParaPerpAsymmPhip_445MeVCM1->Write();
  ParaPerpAsymmPhip_455MeVCM1->Write();
  ParaPerpAsymmPhip_465MeVCM1->Write();
  ParaPerpAsymmPhip_475MeVCM1->Write();
  ParaPerpAsymmPhip_485MeVCM1->Write();
  ParaPerpAsymmPhip_495MeVCM1->Write();
  ParaPerpAsymmPhip_505MeVCM1->Write();
  ParaPerpAsymmPhip_515MeVCM1->Write();
  ParaPerpAsymmPhip_525MeVCM1->Write();
  ParaPerpAsymmPhip_535MeVCM1->Write();
  ParaPerpAsymmPhip_545MeVCM1->Write();
  ParaPerpAsymmPhip_555MeVCM1->Write();
  ParaPerpAsymmPhip_565MeVCM1->Write();
  ParaPerpAsymmPhip_575MeVCM1->Write();
  ParaPerpAsymmPhip_585MeVCM1->Write();
  ParaPerpAsymmPhip_595MeVCM1->Write();
  ParaPerpAsymmPhip_605MeVCM1->Write();
  ParaPerpAsymmPhip_615MeVCM1->Write();

  ParaPerpAsymmPhip_425MeVCM2->Write();
  ParaPerpAsymmPhip_435MeVCM2->Write();
  ParaPerpAsymmPhip_445MeVCM2->Write();
  ParaPerpAsymmPhip_455MeVCM2->Write();
  ParaPerpAsymmPhip_465MeVCM2->Write();
  ParaPerpAsymmPhip_475MeVCM2->Write();
  ParaPerpAsymmPhip_485MeVCM2->Write();
  ParaPerpAsymmPhip_495MeVCM2->Write();
  ParaPerpAsymmPhip_505MeVCM2->Write();
  ParaPerpAsymmPhip_515MeVCM2->Write();
  ParaPerpAsymmPhip_525MeVCM2->Write();
  ParaPerpAsymmPhip_535MeVCM2->Write();
  ParaPerpAsymmPhip_545MeVCM2->Write();
  ParaPerpAsymmPhip_555MeVCM2->Write();
  ParaPerpAsymmPhip_565MeVCM2->Write();
  ParaPerpAsymmPhip_575MeVCM2->Write();
  ParaPerpAsymmPhip_585MeVCM2->Write();
  ParaPerpAsymmPhip_595MeVCM2->Write();
  ParaPerpAsymmPhip_605MeVCM2->Write();
  ParaPerpAsymmPhip_615MeVCM2->Write();

  ParaPerpAsymmPhip_425MeVCM3->Write();
  ParaPerpAsymmPhip_435MeVCM3->Write();
  ParaPerpAsymmPhip_445MeVCM3->Write();
  ParaPerpAsymmPhip_455MeVCM3->Write();
  ParaPerpAsymmPhip_465MeVCM3->Write();
  ParaPerpAsymmPhip_475MeVCM3->Write();
  ParaPerpAsymmPhip_485MeVCM3->Write();
  ParaPerpAsymmPhip_495MeVCM3->Write();
  ParaPerpAsymmPhip_505MeVCM3->Write();
  ParaPerpAsymmPhip_515MeVCM3->Write();
  ParaPerpAsymmPhip_525MeVCM3->Write();
  ParaPerpAsymmPhip_535MeVCM3->Write();
  ParaPerpAsymmPhip_545MeVCM3->Write();
  ParaPerpAsymmPhip_555MeVCM3->Write();
  ParaPerpAsymmPhip_565MeVCM3->Write();
  ParaPerpAsymmPhip_575MeVCM3->Write();
  ParaPerpAsymmPhip_585MeVCM3->Write();
  ParaPerpAsymmPhip_595MeVCM3->Write();
  ParaPerpAsymmPhip_605MeVCM3->Write();
  ParaPerpAsymmPhip_615MeVCM3->Write();

  ParaPerpAsymmPhip_425MeVCM4->Write();
  ParaPerpAsymmPhip_435MeVCM4->Write();
  ParaPerpAsymmPhip_445MeVCM4->Write();
  ParaPerpAsymmPhip_455MeVCM4->Write();
  ParaPerpAsymmPhip_465MeVCM4->Write();
  ParaPerpAsymmPhip_475MeVCM4->Write();
  ParaPerpAsymmPhip_485MeVCM4->Write();
  ParaPerpAsymmPhip_495MeVCM4->Write();
  ParaPerpAsymmPhip_505MeVCM4->Write();
  ParaPerpAsymmPhip_515MeVCM4->Write();
  ParaPerpAsymmPhip_525MeVCM4->Write();
  ParaPerpAsymmPhip_535MeVCM4->Write();
  ParaPerpAsymmPhip_545MeVCM4->Write();
  ParaPerpAsymmPhip_555MeVCM4->Write();
  ParaPerpAsymmPhip_565MeVCM4->Write();
  ParaPerpAsymmPhip_575MeVCM4->Write();
  ParaPerpAsymmPhip_585MeVCM4->Write();
  ParaPerpAsymmPhip_595MeVCM4->Write();
  ParaPerpAsymmPhip_605MeVCM4->Write();
  ParaPerpAsymmPhip_615MeVCM4->Write();

  ParaPerpAsymmPhip_425MeVCM5->Write();
  ParaPerpAsymmPhip_435MeVCM5->Write();
  ParaPerpAsymmPhip_445MeVCM5->Write();
  ParaPerpAsymmPhip_455MeVCM5->Write();
  ParaPerpAsymmPhip_465MeVCM5->Write();
  ParaPerpAsymmPhip_475MeVCM5->Write();
  ParaPerpAsymmPhip_485MeVCM5->Write();
  ParaPerpAsymmPhip_495MeVCM5->Write();
  ParaPerpAsymmPhip_505MeVCM5->Write();
  ParaPerpAsymmPhip_515MeVCM5->Write();
  ParaPerpAsymmPhip_525MeVCM5->Write();
  ParaPerpAsymmPhip_535MeVCM5->Write();
  ParaPerpAsymmPhip_545MeVCM5->Write();
  ParaPerpAsymmPhip_555MeVCM5->Write();
  ParaPerpAsymmPhip_565MeVCM5->Write();
  ParaPerpAsymmPhip_575MeVCM5->Write();
  ParaPerpAsymmPhip_585MeVCM5->Write();
  ParaPerpAsymmPhip_595MeVCM5->Write();
  ParaPerpAsymmPhip_605MeVCM5->Write();
  ParaPerpAsymmPhip_615MeVCM5->Write();

  ParaPerpAsymmPhip_425MeVCM6->Write();
  ParaPerpAsymmPhip_435MeVCM6->Write();
  ParaPerpAsymmPhip_445MeVCM6->Write();
  ParaPerpAsymmPhip_455MeVCM6->Write();
  ParaPerpAsymmPhip_465MeVCM6->Write();
  ParaPerpAsymmPhip_475MeVCM6->Write();
  ParaPerpAsymmPhip_485MeVCM6->Write();
  ParaPerpAsymmPhip_495MeVCM6->Write();
  ParaPerpAsymmPhip_505MeVCM6->Write();
  ParaPerpAsymmPhip_515MeVCM6->Write();
  ParaPerpAsymmPhip_525MeVCM6->Write();
  ParaPerpAsymmPhip_535MeVCM6->Write();
  ParaPerpAsymmPhip_545MeVCM6->Write();
  ParaPerpAsymmPhip_555MeVCM6->Write();
  ParaPerpAsymmPhip_565MeVCM6->Write();
  ParaPerpAsymmPhip_575MeVCM6->Write();
  ParaPerpAsymmPhip_585MeVCM6->Write();
  ParaPerpAsymmPhip_595MeVCM6->Write();
  ParaPerpAsymmPhip_605MeVCM6->Write();
  ParaPerpAsymmPhip_615MeVCM6->Write();

  ParaPerpAsymmPhip_425MeVCM7->Write();
  ParaPerpAsymmPhip_435MeVCM7->Write();
  ParaPerpAsymmPhip_445MeVCM7->Write();
  ParaPerpAsymmPhip_455MeVCM7->Write();
  ParaPerpAsymmPhip_465MeVCM7->Write();
  ParaPerpAsymmPhip_475MeVCM7->Write();
  ParaPerpAsymmPhip_485MeVCM7->Write();
  ParaPerpAsymmPhip_495MeVCM7->Write();
  ParaPerpAsymmPhip_505MeVCM7->Write();
  ParaPerpAsymmPhip_515MeVCM7->Write();
  ParaPerpAsymmPhip_525MeVCM7->Write();
  ParaPerpAsymmPhip_535MeVCM7->Write();
  ParaPerpAsymmPhip_545MeVCM7->Write();
  ParaPerpAsymmPhip_555MeVCM7->Write();
  ParaPerpAsymmPhip_565MeVCM7->Write();
  ParaPerpAsymmPhip_575MeVCM7->Write();
  ParaPerpAsymmPhip_585MeVCM7->Write();
  ParaPerpAsymmPhip_595MeVCM7->Write();
  ParaPerpAsymmPhip_605MeVCM7->Write();
  ParaPerpAsymmPhip_615MeVCM7->Write();

  ParaPerpAsymmPhip_425MeVCM8->Write();
  ParaPerpAsymmPhip_435MeVCM8->Write();
  ParaPerpAsymmPhip_445MeVCM8->Write();
  ParaPerpAsymmPhip_455MeVCM8->Write();
  ParaPerpAsymmPhip_465MeVCM8->Write();
  ParaPerpAsymmPhip_475MeVCM8->Write();
  ParaPerpAsymmPhip_485MeVCM8->Write();
  ParaPerpAsymmPhip_495MeVCM8->Write();
  ParaPerpAsymmPhip_505MeVCM8->Write();
  ParaPerpAsymmPhip_515MeVCM8->Write();
  ParaPerpAsymmPhip_525MeVCM8->Write();
  ParaPerpAsymmPhip_535MeVCM8->Write();
  ParaPerpAsymmPhip_545MeVCM8->Write();
  ParaPerpAsymmPhip_555MeVCM8->Write();
  ParaPerpAsymmPhip_565MeVCM8->Write();
  ParaPerpAsymmPhip_575MeVCM8->Write();
  ParaPerpAsymmPhip_585MeVCM8->Write();
  ParaPerpAsymmPhip_595MeVCM8->Write();
  ParaPerpAsymmPhip_605MeVCM8->Write();
  ParaPerpAsymmPhip_615MeVCM8->Write();

  ParaPerpAsymmPhip_425MeVCM9->Write();
  ParaPerpAsymmPhip_435MeVCM9->Write();
  ParaPerpAsymmPhip_445MeVCM9->Write();
  ParaPerpAsymmPhip_455MeVCM9->Write();
  ParaPerpAsymmPhip_465MeVCM9->Write();
  ParaPerpAsymmPhip_475MeVCM9->Write();
  ParaPerpAsymmPhip_485MeVCM9->Write();
  ParaPerpAsymmPhip_495MeVCM9->Write();
  ParaPerpAsymmPhip_505MeVCM9->Write();
  ParaPerpAsymmPhip_515MeVCM9->Write();
  ParaPerpAsymmPhip_525MeVCM9->Write();
  ParaPerpAsymmPhip_535MeVCM9->Write();
  ParaPerpAsymmPhip_545MeVCM9->Write();
  ParaPerpAsymmPhip_555MeVCM9->Write();
  ParaPerpAsymmPhip_565MeVCM9->Write();
  ParaPerpAsymmPhip_575MeVCM9->Write();
  ParaPerpAsymmPhip_585MeVCM9->Write();
  ParaPerpAsymmPhip_595MeVCM9->Write();
  ParaPerpAsymmPhip_605MeVCM9->Write();
  ParaPerpAsymmPhip_615MeVCM9->Write();

  ParaPerpAsymmPhip_425MeVCM10->Write();
  ParaPerpAsymmPhip_435MeVCM10->Write();
  ParaPerpAsymmPhip_445MeVCM10->Write();
  ParaPerpAsymmPhip_455MeVCM10->Write();
  ParaPerpAsymmPhip_465MeVCM10->Write();
  ParaPerpAsymmPhip_475MeVCM10->Write();
  ParaPerpAsymmPhip_485MeVCM10->Write();
  ParaPerpAsymmPhip_495MeVCM10->Write();
  ParaPerpAsymmPhip_505MeVCM10->Write();
  ParaPerpAsymmPhip_515MeVCM10->Write();
  ParaPerpAsymmPhip_525MeVCM10->Write();
  ParaPerpAsymmPhip_535MeVCM10->Write();
  ParaPerpAsymmPhip_545MeVCM10->Write();
  ParaPerpAsymmPhip_555MeVCM10->Write();
  ParaPerpAsymmPhip_565MeVCM10->Write();
  ParaPerpAsymmPhip_575MeVCM10->Write();
  ParaPerpAsymmPhip_585MeVCM10->Write();
  ParaPerpAsymmPhip_595MeVCM10->Write();
  ParaPerpAsymmPhip_605MeVCM10->Write();
  ParaPerpAsymmPhip_615MeVCM10->Write();

  //Define new tree to store parameters in
  TTree* tree = new TTree("Parameter_Values", "Tree_of_Values");

  // Define branches to store parameters, (Branch Name, Variable, Type of Variable)
  tree->Branch("pCosAmp425", &pCosA425, "pCosA425/D");
  tree->Branch("pCosAmpErr425", &pCosAErr425, "pCosAErr425/D");
  tree->Branch("pCosAmp435", &pCosA435, "pCosA435/D");
  tree->Branch("pCosAmpErr435", &pCosAErr435, "pCosAErr435/D");
  tree->Branch("pCosAmp445", &pCosA445, "pCosA445/D");
  tree->Branch("pCosAmpErr445", &pCosAErr445, "pCosAErr445/D");
  tree->Branch("pCosAmp455", &pCosA455, "pCosA455/D");
  tree->Branch("pCosAmpErr455", &pCosAErr455, "pCosAErr455/D");
  tree->Branch("pCosAmp465", &pCosA465, "pCosA465/D");
  tree->Branch("pCosAmpErr465", &pCosAErr465, "pCosAErr465/D");
  tree->Branch("pCosAmp475", &pCosA475, "pCosA475/D");
  tree->Branch("pCosAmpErr475", &pCosAErr475, "pCosAErr475/D");
  tree->Branch("pCosAmp485", &pCosA485, "pCosA485/D");
  tree->Branch("pCosAmpErr485", &pCosAErr485, "pCosAErr485/D");
  tree->Branch("pCosAmp495", &pCosA495, "pCosA495/D");
  tree->Branch("pCosAmpErr495", &pCosAErr495, "pCosAErr495/D");
  tree->Branch("pCosAmp505", &pCosA505, "pCosA505/D");
  tree->Branch("pCosAmpErr505", &pCosAErr505, "pCosAErr505/D");
  tree->Branch("pCosAmp515", &pCosA515, "pCosA515/D");
  tree->Branch("pCosAmpErr515", &pCosAErr515, "pCosAErr515/D");
  tree->Branch("pCosAmp525", &pCosA525, "pCosA525/D");
  tree->Branch("pCosAmpErr525", &pCosAErr525, "pCosAErr525/D");
  tree->Branch("pCosAmp535", &pCosA535, "pCosA535/D");
  tree->Branch("pCosAmpErr535", &pCosAErr535, "pCosAErr535/D");
  tree->Branch("pCosAmp545", &pCosA545, "pCosA545/D");
  tree->Branch("pCosAmpErr545", &pCosAErr545, "pCosAErr545/D");
  tree->Branch("pCosAmp555", &pCosA555, "pCosA555/D");
  tree->Branch("pCosAmpErr555", &pCosAErr555, "pCosAErr555/D");
  tree->Branch("pCosAmp565", &pCosA565, "pCosA565/D");
  tree->Branch("pCosAmpErr565", &pCosAErr565, "pCosAErr565/D");
  tree->Branch("pCosAmp575", &pCosA575, "pCosA575/D");
  tree->Branch("pCosAmpErr575", &pCosAErr575, "pCosAErr575/D");
  tree->Branch("pCosAmp585", &pCosA585, "pCosA585/D");
  tree->Branch("pCosAmpErr585", &pCosAErr585, "pCosAErr585/D");
  tree->Branch("pCosAmp595", &pCosA595, "pCosA595/D");
  tree->Branch("pCosAmpErr595", &pCosAErr595, "pCosAErr595/D");
  tree->Branch("pCosAmp605", &pCosA605, "pCosA605/D");
  tree->Branch("pCosAmpErr605", &pCosAErr605, "pCosAErr605/D");
  tree->Branch("pCosAmp615", &pCosA615, "pCosA615/D");
  tree->Branch("pCosAmpErr615", &pCosAErr615, "pCosAErr615/D");

  // Fill branches (and hence tree) with corresponding parameters from above
  for (Int_t m = 0; m < 10; m++){
    pCosA425 = pCosAmp[m][0];
    pCosAErr425 = pCosAmpErr[m][0];
    pCosA435 = pCosAmp[m][1];
    pCosAErr435 = pCosAmpErr[m][1];
    pCosA445 = pCosAmp[m][2];
    pCosAErr445 = pCosAmpErr[m][2];
    pCosA455 = pCosAmp[m][3];
    pCosAErr455 = pCosAmpErr[m][3];
    pCosA465 = pCosAmp[m][4];
    pCosAErr465 = pCosAmpErr[m][4];
    pCosA475 = pCosAmp[m][5];
    pCosAErr475 = pCosAmpErr[m][5];
    pCosA485 = pCosAmp[m][6];
    pCosAErr485 = pCosAmpErr[m][6];
    pCosA495 = pCosAmp[m][7];
    pCosAErr495 = pCosAmpErr[m][7];
    pCosA505 = pCosAmp[m][8];
    pCosAErr505 = pCosAmpErr[m][8];
    pCosA515 = pCosAmp[m][9];
    pCosAErr515 = pCosAmpErr[m][9];
    pCosA525 = pCosAmp[m][10];
    pCosAErr525 = pCosAmpErr[m][10];
    pCosA535 = pCosAmp[m][11];
    pCosAErr535 = pCosAmpErr[m][11];
    pCosA545 = pCosAmp[m][12];
    pCosAErr545= pCosAmpErr[m][12];
    pCosA555 = pCosAmp[m][13];
    pCosAErr555 = pCosAmpErr[m][13];
    pCosA565 = pCosAmp[m][14];
    pCosAErr565 = pCosAmpErr[m][14];
    pCosA575 = pCosAmp[m][15];
    pCosAErr575 = pCosAmpErr[m][15];
    pCosA585 = pCosAmp[m][16];
    pCosAErr585 = pCosAmpErr[m][16];
    pCosA595 = pCosAmp[m][17];
    pCosAErr595 = pCosAmpErr[m][17];
    pCosA605 = pCosAmp[m][18];
    pCosAErr605 = pCosAmpErr[m][18];
    pCosA615 = pCosAmp[m][19];
    pCosAErr615 = pCosAmpErr[m][19];

    tree->Fill();
  }

  f1.Write();

}
