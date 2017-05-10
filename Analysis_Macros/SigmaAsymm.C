#include "./includes_SigmaAsymm.h"

void SigmaAsymm(){

  TF1 *CosFunc = new TF1("CosFit", "[0]*cos((2*x*TMath::DegToRad())+acos(0))");
  CosFunc->SetParNames("Amplitude");

  double pCosAmp[10][12]; // Format of array is Theta bin (x) by Egamma bin (y), 6 theta bins of 30, 12 20MeV Egamma bins
  double pCosAmpErr[10][12];
  double pCosA410;
  double pCosAErr410;
  double pCosA430;
  double pCosAErr430;
  double pCosA450;
  double pCosAErr450;
  double pCosA470;
  double pCosAErr470;
  double pCosA490;
  double pCosAErr490;
  double pCosA510;
  double pCosAErr510;
  double pCosA530;
  double pCosAErr530;
  double pCosA550;
  double pCosAErr550;
  double pCosA570;
  double pCosAErr570;
  double pCosA590;
  double pCosAErr590;
  double pCosA610;
  double pCosAErr610;
  double pCosA630;
  double pCosAErr630;

  TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/ParaPerp_Total_10_Combined.root"); // Open the latest PTotal combined file to load histograms from
  NPara = Eg_Para->GetEntries();
  NPerp = Eg_Perp->GetEntries();
  ScaleFactor = NPara/NPerp;
  ScaleFactorErr = sqrt( (NPara/((TMath::Power(NPerp,2)))) + (((TMath::Power(NPara,2)))/(TMath::Power(NPerp,3))) ); // Error Propagation of above formula, see notebook 5

  ParaPerpAsymmPhip_410MeVCM1 = Phip_410MeVCM1_Para->GetAsymmetry(Phip_410MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_410MeVCM1->SetName("ParaPerpAsymmPhip410MeVCM1");
  ParaPerpAsymmPhip_410MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for 400-420MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_410MeVCM1->Fit("CosFit", "LL");
  pCosAmp[0][0] = CosFit->GetParameter(0);
  pCosAmpErr[0][0] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_430MeVCM1 = Phip_430MeVCM1_Para->GetAsymmetry(Phip_430MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_430MeVCM1->SetName("ParaPerpAsymmPhip430MeVCM1");
  ParaPerpAsymmPhip_430MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for 420-440MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_430MeVCM1->Fit("CosFit", "LL");
  pCosAmp[0][1] = CosFit->GetParameter(0);
  pCosAmpErr[0][1] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_450MeVCM1 = Phip_450MeVCM1_Para->GetAsymmetry(Phip_450MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_450MeVCM1->SetName("ParaPerpAsymmPhip450MeVCM1");
  ParaPerpAsymmPhip_450MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for 440-460MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_450MeVCM1->Fit("CosFit", "LL");
  pCosAmp[0][2] = CosFit->GetParameter(0);
  pCosAmpErr[0][2] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_470MeVCM1 = Phip_470MeVCM1_Para->GetAsymmetry(Phip_470MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_470MeVCM1->SetName("ParaPerpAsymmPhip470MeVCM1");
  ParaPerpAsymmPhip_470MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for 460-480MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_470MeVCM1->Fit("CosFit", "LL");
  pCosAmp[0][3] = CosFit->GetParameter(0);
  pCosAmpErr[0][3] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_490MeVCM1 = Phip_490MeVCM1_Para->GetAsymmetry(Phip_490MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_490MeVCM1->SetName("ParaPerpAsymmPhip490MeVCM1");
  ParaPerpAsymmPhip_490MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for 480-500MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_490MeVCM1->Fit("CosFit", "LL");
  pCosAmp[0][4] = CosFit->GetParameter(0);
  pCosAmpErr[0][4] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_510MeVCM1 = Phip_510MeVCM1_Para->GetAsymmetry(Phip_510MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_510MeVCM1->SetName("ParaPerpAsymmPhip510MeVCM1");
  ParaPerpAsymmPhip_510MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for 500-520MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_510MeVCM1->Fit("CosFit", "LL");
  pCosAmp[0][5] = CosFit->GetParameter(0);
  pCosAmpErr[0][5] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_530MeVCM1 = Phip_530MeVCM1_Para->GetAsymmetry(Phip_530MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_530MeVCM1->SetName("ParaPerpAsymmPhip530MeVCM1");
  ParaPerpAsymmPhip_530MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for 520-540MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_530MeVCM1->Fit("CosFit", "LL");
  pCosAmp[0][6] = CosFit->GetParameter(0);
  pCosAmpErr[0][6] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_550MeVCM1 = Phip_550MeVCM1_Para->GetAsymmetry(Phip_550MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_550MeVCM1->SetName("ParaPerpAsymmPhip550MeVCM1");
  ParaPerpAsymmPhip_550MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for 540-560MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_550MeVCM1->Fit("CosFit", "LL");
  pCosAmp[0][7] = CosFit->GetParameter(0);
  pCosAmpErr[0][7] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_570MeVCM1 = Phip_570MeVCM1_Para->GetAsymmetry(Phip_570MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_570MeVCM1->SetName("ParaPerpAsymmPhip570MeVCM1");
  ParaPerpAsymmPhip_570MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for 560-580MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_570MeVCM1->Fit("CosFit", "LL");
  pCosAmp[0][8] = CosFit->GetParameter(0);
  pCosAmpErr[0][8] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_590MeVCM1 = Phip_590MeVCM1_Para->GetAsymmetry(Phip_590MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_590MeVCM1->SetName("ParaPerpAsymmPhip590MeVCM1");
  ParaPerpAsymmPhip_590MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for 580-600MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_590MeVCM1->Fit("CosFit", "LL");
  pCosAmp[0][9] = CosFit->GetParameter(0);
  pCosAmpErr[0][9] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_610MeVCM1 = Phip_610MeVCM1_Para->GetAsymmetry(Phip_610MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_610MeVCM1->SetName("ParaPerpAsymmPhip610MeVCM1");
  ParaPerpAsymmPhip_610MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for 600-620MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_610MeVCM1->Fit("CosFit", "LL");
  pCosAmp[0][10] = CosFit->GetParameter(0);
  pCosAmpErr[0][10] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_630MeVCM1 = Phip_630MeVCM1_Para->GetAsymmetry(Phip_630MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_630MeVCM1->SetName("ParaPerpAsymmPhip630MeVCM1");
  ParaPerpAsymmPhip_630MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for 620-640MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_630MeVCM1->Fit("CosFit", "LL");
  pCosAmp[0][11] = CosFit->GetParameter(0);
  pCosAmpErr[0][11] = CosFit->GetParError(0);

  ///////////////////////////////////////////
  //////////////////CM2 p ///////////////////
  ///////////////////////////////////////////

  ParaPerpAsymmPhip_410MeVCM2 = Phip_410MeVCM2_Para->GetAsymmetry(Phip_410MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_410MeVCM2->SetName("ParaPerpAsymmPhip410MeVCM2");
  ParaPerpAsymmPhip_410MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for 400-420MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_410MeVCM2->Fit("CosFit", "LL");
  pCosAmp[1][0] = CosFit->GetParameter(0);
  pCosAmpErr[1][0] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_430MeVCM2 = Phip_430MeVCM2_Para->GetAsymmetry(Phip_430MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_430MeVCM2->SetName("ParaPerpAsymmPhip430MeVCM2");
  ParaPerpAsymmPhip_430MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for 420-440MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_430MeVCM2->Fit("CosFit", "LL");
  pCosAmp[1][1] = CosFit->GetParameter(0);
  pCosAmpErr[1][1] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_450MeVCM2 = Phip_450MeVCM2_Para->GetAsymmetry(Phip_450MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_450MeVCM2->SetName("ParaPerpAsymmPhip450MeVCM2");
  ParaPerpAsymmPhip_450MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for 440-460MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_450MeVCM2->Fit("CosFit", "LL");
  pCosAmp[1][2] = CosFit->GetParameter(0);
  pCosAmpErr[1][2] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_470MeVCM2 = Phip_470MeVCM2_Para->GetAsymmetry(Phip_470MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_470MeVCM2->SetName("ParaPerpAsymmPhip470MeVCM2");
  ParaPerpAsymmPhip_470MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for 460-480MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_470MeVCM2->Fit("CosFit", "LL");
  pCosAmp[1][3] = CosFit->GetParameter(0);
  pCosAmpErr[1][3] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_490MeVCM2 = Phip_490MeVCM2_Para->GetAsymmetry(Phip_490MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_490MeVCM2->SetName("ParaPerpAsymmPhip490MeVCM2");
  ParaPerpAsymmPhip_490MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for 480-500MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_490MeVCM2->Fit("CosFit", "LL");
  pCosAmp[1][4] = CosFit->GetParameter(0);
  pCosAmpErr[1][4] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_510MeVCM2 = Phip_510MeVCM2_Para->GetAsymmetry(Phip_510MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_510MeVCM2->SetName("ParaPerpAsymmPhip510MeVCM2");
  ParaPerpAsymmPhip_510MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for 500-520MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_510MeVCM2->Fit("CosFit", "LL");
  pCosAmp[1][5] = CosFit->GetParameter(0);
  pCosAmpErr[1][5] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_530MeVCM2 = Phip_530MeVCM2_Para->GetAsymmetry(Phip_530MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_530MeVCM2->SetName("ParaPerpAsymmPhip530MeVCM2");
  ParaPerpAsymmPhip_530MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for 520-540MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_530MeVCM2->Fit("CosFit", "LL");
  pCosAmp[1][6] = CosFit->GetParameter(0);
  pCosAmpErr[1][6] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_550MeVCM2 = Phip_550MeVCM2_Para->GetAsymmetry(Phip_550MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_550MeVCM2->SetName("ParaPerpAsymmPhip550MeVCM2");
  ParaPerpAsymmPhip_550MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for 540-560MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_550MeVCM2->Fit("CosFit", "LL");
  pCosAmp[1][7] = CosFit->GetParameter(0);
  pCosAmpErr[1][7] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_570MeVCM2 = Phip_570MeVCM2_Para->GetAsymmetry(Phip_570MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_570MeVCM2->SetName("ParaPerpAsymmPhip570MeVCM2");
  ParaPerpAsymmPhip_570MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for 560-580MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_570MeVCM2->Fit("CosFit", "LL");
  pCosAmp[1][8] = CosFit->GetParameter(0);
  pCosAmpErr[1][8] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_590MeVCM2 = Phip_590MeVCM2_Para->GetAsymmetry(Phip_590MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_590MeVCM2->SetName("ParaPerpAsymmPhip590MeVCM2");
  ParaPerpAsymmPhip_590MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for 580-600MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_590MeVCM2->Fit("CosFit", "LL");
  pCosAmp[1][9] = CosFit->GetParameter(0);
  pCosAmpErr[1][9] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_610MeVCM2 = Phip_610MeVCM2_Para->GetAsymmetry(Phip_610MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_610MeVCM2->SetName("ParaPerpAsymmPhip610MeVCM2");
  ParaPerpAsymmPhip_610MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for 600-620MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_610MeVCM2->Fit("CosFit", "LL");
  pCosAmp[1][10] = CosFit->GetParameter(0);
  pCosAmpErr[1][10] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_630MeVCM2 = Phip_630MeVCM2_Para->GetAsymmetry(Phip_630MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_630MeVCM2->SetName("ParaPerpAsymmPhip630MeVCM2");
  ParaPerpAsymmPhip_630MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for 620-640MeV Photon Energy (CosTheta1-0.8)");
  ParaPerpAsymmPhip_630MeVCM2->Fit("CosFit", "LL");
  pCosAmp[1][11] = CosFit->GetParameter(0);
  pCosAmpErr[1][11] = CosFit->GetParError(0);

  ///////////////////////////////////////////
  //////////////////CM3 p ///////////////////
  ///////////////////////////////////////////

  ParaPerpAsymmPhip_410MeVCM3 = Phip_410MeVCM3_Para->GetAsymmetry(Phip_410MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_410MeVCM3->SetName("ParaPerpAsymmPhip410MeVCM3");
  ParaPerpAsymmPhip_410MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for 400-420MeV Photon Energy (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_410MeVCM3->Fit("CosFit", "LL");
  pCosAmp[2][0] = CosFit->GetParameter(0);
  pCosAmpErr[2][0] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_430MeVCM3 = Phip_430MeVCM3_Para->GetAsymmetry(Phip_430MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_430MeVCM3->SetName("ParaPerpAsymmPhip430MeVCM3");
  ParaPerpAsymmPhip_430MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for 420-440MeV Photon Energy (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_430MeVCM3->Fit("CosFit", "LL");
  pCosAmp[2][1] = CosFit->GetParameter(0);
  pCosAmpErr[2][1] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_450MeVCM3 = Phip_450MeVCM3_Para->GetAsymmetry(Phip_450MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_450MeVCM3->SetName("ParaPerpAsymmPhip450MeVCM3");
  ParaPerpAsymmPhip_450MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for 440-460MeV Photon Energy (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_450MeVCM3->Fit("CosFit", "LL");
  pCosAmp[2][2] = CosFit->GetParameter(0);
  pCosAmpErr[2][2] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_470MeVCM3 = Phip_470MeVCM3_Para->GetAsymmetry(Phip_470MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_470MeVCM3->SetName("ParaPerpAsymmPhip470MeVCM3");
  ParaPerpAsymmPhip_470MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for 460-480MeV Photon Energy (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_470MeVCM3->Fit("CosFit", "LL");
  pCosAmp[2][3] = CosFit->GetParameter(0);
  pCosAmpErr[2][3] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_490MeVCM3 = Phip_490MeVCM3_Para->GetAsymmetry(Phip_490MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_490MeVCM3->SetName("ParaPerpAsymmPhip490MeVCM3");
  ParaPerpAsymmPhip_490MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for 480-500MeV Photon Energy (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_490MeVCM3->Fit("CosFit", "LL");
  pCosAmp[2][4] = CosFit->GetParameter(0);
  pCosAmpErr[2][4] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_510MeVCM3 = Phip_510MeVCM3_Para->GetAsymmetry(Phip_510MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_510MeVCM3->SetName("ParaPerpAsymmPhip510MeVCM3");
  ParaPerpAsymmPhip_510MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for 500-520MeV Photon Energy (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_510MeVCM3->Fit("CosFit", "LL");
  pCosAmp[2][5] = CosFit->GetParameter(0);
  pCosAmpErr[2][5] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_530MeVCM3 = Phip_530MeVCM3_Para->GetAsymmetry(Phip_530MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_530MeVCM3->SetName("ParaPerpAsymmPhip530MeVCM3");
  ParaPerpAsymmPhip_530MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for 520-540MeV Photon Energy (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_530MeVCM3->Fit("CosFit", "LL");
  pCosAmp[2][6] = CosFit->GetParameter(0);
  pCosAmpErr[2][6] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_550MeVCM3 = Phip_550MeVCM3_Para->GetAsymmetry(Phip_550MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_550MeVCM3->SetName("ParaPerpAsymmPhip550MeVCM3");
  ParaPerpAsymmPhip_550MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for 540-560MeV Photon Energy (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_550MeVCM3->Fit("CosFit", "LL");
  pCosAmp[2][7] = CosFit->GetParameter(0);
  pCosAmpErr[2][7] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_570MeVCM3 = Phip_570MeVCM3_Para->GetAsymmetry(Phip_570MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_570MeVCM3->SetName("ParaPerpAsymmPhip570MeVCM3");
  ParaPerpAsymmPhip_570MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for 560-580MeV Photon Energy (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_570MeVCM3->Fit("CosFit", "LL");
  pCosAmp[2][8] = CosFit->GetParameter(0);
  pCosAmpErr[2][8] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_590MeVCM3 = Phip_590MeVCM3_Para->GetAsymmetry(Phip_590MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_590MeVCM3->SetName("ParaPerpAsymmPhip590MeVCM3");
  ParaPerpAsymmPhip_590MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for 580-600MeV Photon Energy (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_590MeVCM3->Fit("CosFit", "LL");
  pCosAmp[2][9] = CosFit->GetParameter(0);
  pCosAmpErr[2][9] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_610MeVCM3 = Phip_610MeVCM3_Para->GetAsymmetry(Phip_610MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_610MeVCM3->SetName("ParaPerpAsymmPhip610MeVCM3");
  ParaPerpAsymmPhip_610MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for 600-620MeV Photon Energy (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_610MeVCM3->Fit("CosFit", "LL");
  pCosAmp[2][10] = CosFit->GetParameter(0);
  pCosAmpErr[2][10] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_630MeVCM3 = Phip_630MeVCM3_Para->GetAsymmetry(Phip_630MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_630MeVCM3->SetName("ParaPerpAsymmPhip630MeVCM3");
  ParaPerpAsymmPhip_630MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for 620-640MeV Photon Energy (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_630MeVCM3->Fit("CosFit", "LL");
  pCosAmp[2][11] = CosFit->GetParameter(0);
  pCosAmpErr[2][11] = CosFit->GetParError(0);

  ///////////////////////////////////////////
  //////////////////CM4 p ///////////////////
  ///////////////////////////////////////////

  ParaPerpAsymmPhip_410MeVCM4 = Phip_410MeVCM4_Para->GetAsymmetry(Phip_410MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_410MeVCM4->SetName("ParaPerpAsymmPhip410MeVCM4");
  ParaPerpAsymmPhip_410MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for 400-420MeV Photon Energy (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_410MeVCM4->Fit("CosFit", "LL");
  pCosAmp[3][0] = CosFit->GetParameter(0);
  pCosAmpErr[3][0] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_430MeVCM4 = Phip_430MeVCM4_Para->GetAsymmetry(Phip_430MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_430MeVCM4->SetName("ParaPerpAsymmPhip430MeVCM4");
  ParaPerpAsymmPhip_430MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for 420-440MeV Photon Energy (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_430MeVCM4->Fit("CosFit", "LL");
  pCosAmp[3][1] = CosFit->GetParameter(0);
  pCosAmpErr[3][1] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_450MeVCM4 = Phip_450MeVCM4_Para->GetAsymmetry(Phip_450MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_450MeVCM4->SetName("ParaPerpAsymmPhip450MeVCM4");
  ParaPerpAsymmPhip_450MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for 440-460MeV Photon Energy (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_450MeVCM4->Fit("CosFit", "LL");
  pCosAmp[3][2] = CosFit->GetParameter(0);
  pCosAmpErr[3][2] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_470MeVCM4 = Phip_470MeVCM4_Para->GetAsymmetry(Phip_470MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_470MeVCM4->SetName("ParaPerpAsymmPhip470MeVCM4");
  ParaPerpAsymmPhip_470MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for 460-480MeV Photon Energy (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_470MeVCM4->Fit("CosFit", "LL");
  pCosAmp[3][3] = CosFit->GetParameter(0);
  pCosAmpErr[3][3] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_490MeVCM4 = Phip_490MeVCM4_Para->GetAsymmetry(Phip_490MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_490MeVCM4->SetName("ParaPerpAsymmPhip490MeVCM4");
  ParaPerpAsymmPhip_490MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for 480-500MeV Photon Energy (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_490MeVCM4->Fit("CosFit", "LL");
  pCosAmp[3][4] = CosFit->GetParameter(0);
  pCosAmpErr[3][4] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_510MeVCM4 = Phip_510MeVCM4_Para->GetAsymmetry(Phip_510MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_510MeVCM4->SetName("ParaPerpAsymmPhip510MeVCM4");
  ParaPerpAsymmPhip_510MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for 500-520MeV Photon Energy (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_510MeVCM4->Fit("CosFit", "LL");
  pCosAmp[3][5] = CosFit->GetParameter(0);
  pCosAmpErr[3][5] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_530MeVCM4 = Phip_530MeVCM4_Para->GetAsymmetry(Phip_530MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_530MeVCM4->SetName("ParaPerpAsymmPhip530MeVCM4");
  ParaPerpAsymmPhip_530MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for 520-540MeV Photon Energy (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_530MeVCM4->Fit("CosFit", "LL");
  pCosAmp[3][6] = CosFit->GetParameter(0);
  pCosAmpErr[3][6] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_550MeVCM4 = Phip_550MeVCM4_Para->GetAsymmetry(Phip_550MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_550MeVCM4->SetName("ParaPerpAsymmPhip550MeVCM4");
  ParaPerpAsymmPhip_550MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for 540-560MeV Photon Energy (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_550MeVCM4->Fit("CosFit", "LL");
  pCosAmp[3][7] = CosFit->GetParameter(0);
  pCosAmpErr[3][7] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_570MeVCM4 = Phip_570MeVCM4_Para->GetAsymmetry(Phip_570MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_570MeVCM4->SetName("ParaPerpAsymmPhip570MeVCM4");
  ParaPerpAsymmPhip_570MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for 560-580MeV Photon Energy (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_570MeVCM4->Fit("CosFit", "LL");
  pCosAmp[3][8] = CosFit->GetParameter(0);
  pCosAmpErr[3][8] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_590MeVCM4 = Phip_590MeVCM4_Para->GetAsymmetry(Phip_590MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_590MeVCM4->SetName("ParaPerpAsymmPhip590MeVCM4");
  ParaPerpAsymmPhip_590MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for 580-600MeV Photon Energy (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_590MeVCM4->Fit("CosFit", "LL");
  pCosAmp[3][9] = CosFit->GetParameter(0);
  pCosAmpErr[3][9] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_610MeVCM4 = Phip_610MeVCM4_Para->GetAsymmetry(Phip_610MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_610MeVCM4->SetName("ParaPerpAsymmPhip610MeVCM4");
  ParaPerpAsymmPhip_610MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for 600-620MeV Photon Energy (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_610MeVCM4->Fit("CosFit", "LL");
  pCosAmp[3][10] = CosFit->GetParameter(0);
  pCosAmpErr[3][10] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_630MeVCM4 = Phip_630MeVCM4_Para->GetAsymmetry(Phip_630MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_630MeVCM4->SetName("ParaPerpAsymmPhip630MeVCM4");
  ParaPerpAsymmPhip_630MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for 620-640MeV Photon Energy (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_630MeVCM4->Fit("CosFit", "LL");
  pCosAmp[3][11] = CosFit->GetParameter(0);
  pCosAmpErr[3][11] = CosFit->GetParError(0);

  ///////////////////////////////////////////
  //////////////////CM5 p ///////////////////
  ///////////////////////////////////////////

  ParaPerpAsymmPhip_410MeVCM5 = Phip_410MeVCM5_Para->GetAsymmetry(Phip_410MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_410MeVCM5->SetName("ParaPerpAsymmPhip410MeVCM5");
  ParaPerpAsymmPhip_410MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for 400-420MeV Photon Energy (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_410MeVCM5->Fit("CosFit", "LL");
  pCosAmp[4][0] = CosFit->GetParameter(0);
  pCosAmpErr[4][0] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_430MeVCM5 = Phip_430MeVCM5_Para->GetAsymmetry(Phip_430MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_430MeVCM5->SetName("ParaPerpAsymmPhip430MeVCM5");
  ParaPerpAsymmPhip_430MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for 420-440MeV Photon Energy (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_430MeVCM5->Fit("CosFit", "LL");
  pCosAmp[4][1] = CosFit->GetParameter(0);
  pCosAmpErr[4][1] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_450MeVCM5 = Phip_450MeVCM5_Para->GetAsymmetry(Phip_450MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_450MeVCM5->SetName("ParaPerpAsymmPhip450MeVCM5");
  ParaPerpAsymmPhip_450MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for 440-460MeV Photon Energy (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_450MeVCM5->Fit("CosFit", "LL");
  pCosAmp[4][2] = CosFit->GetParameter(0);
  pCosAmpErr[4][2] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_470MeVCM5 = Phip_470MeVCM5_Para->GetAsymmetry(Phip_470MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_470MeVCM5->SetName("ParaPerpAsymmPhip470MeVCM5");
  ParaPerpAsymmPhip_470MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for 460-480MeV Photon Energy (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_470MeVCM5->Fit("CosFit", "LL");
  pCosAmp[4][3] = CosFit->GetParameter(0);
  pCosAmpErr[4][3] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_490MeVCM5 = Phip_490MeVCM5_Para->GetAsymmetry(Phip_490MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_490MeVCM5->SetName("ParaPerpAsymmPhip490MeVCM5");
  ParaPerpAsymmPhip_490MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for 480-500MeV Photon Energy (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_490MeVCM5->Fit("CosFit", "LL");
  pCosAmp[4][4] = CosFit->GetParameter(0);
  pCosAmpErr[4][4] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_510MeVCM5 = Phip_510MeVCM5_Para->GetAsymmetry(Phip_510MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_510MeVCM5->SetName("ParaPerpAsymmPhip510MeVCM5");
  ParaPerpAsymmPhip_510MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for 500-520MeV Photon Energy (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_510MeVCM5->Fit("CosFit", "LL");
  pCosAmp[4][5] = CosFit->GetParameter(0);
  pCosAmpErr[4][5] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_530MeVCM5 = Phip_530MeVCM5_Para->GetAsymmetry(Phip_530MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_530MeVCM5->SetName("ParaPerpAsymmPhip530MeVCM5");
  ParaPerpAsymmPhip_530MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for 520-540MeV Photon Energy (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_530MeVCM5->Fit("CosFit", "LL");
  pCosAmp[4][6] = CosFit->GetParameter(0);
  pCosAmpErr[4][6] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_550MeVCM5 = Phip_550MeVCM5_Para->GetAsymmetry(Phip_550MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_550MeVCM5->SetName("ParaPerpAsymmPhip550MeVCM5");
  ParaPerpAsymmPhip_550MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for 540-560MeV Photon Energy (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_550MeVCM5->Fit("CosFit", "LL");
  pCosAmp[4][7] = CosFit->GetParameter(0);
  pCosAmpErr[4][7] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_570MeVCM5 = Phip_570MeVCM5_Para->GetAsymmetry(Phip_570MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_570MeVCM5->SetName("ParaPerpAsymmPhip570MeVCM5");
  ParaPerpAsymmPhip_570MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for 560-580MeV Photon Energy (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_570MeVCM5->Fit("CosFit", "LL");
  pCosAmp[4][8] = CosFit->GetParameter(0);
  pCosAmpErr[4][8] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_590MeVCM5 = Phip_590MeVCM5_Para->GetAsymmetry(Phip_590MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_590MeVCM5->SetName("ParaPerpAsymmPhip590MeVCM5");
  ParaPerpAsymmPhip_590MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for 580-600MeV Photon Energy (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_590MeVCM5->Fit("CosFit", "LL");
  pCosAmp[4][9] = CosFit->GetParameter(0);
  pCosAmpErr[4][9] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_610MeVCM5 = Phip_610MeVCM5_Para->GetAsymmetry(Phip_610MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_610MeVCM5->SetName("ParaPerpAsymmPhip610MeVCM5");
  ParaPerpAsymmPhip_610MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for 600-620MeV Photon Energy (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_610MeVCM5->Fit("CosFit", "LL");
  pCosAmp[4][10] = CosFit->GetParameter(0);
  pCosAmpErr[4][10] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_630MeVCM5 = Phip_630MeVCM5_Para->GetAsymmetry(Phip_630MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_630MeVCM5->SetName("ParaPerpAsymmPhip630MeVCM5");
  ParaPerpAsymmPhip_630MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for 620-640MeV Photon Energy (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_630MeVCM5->Fit("CosFit", "LL");
  pCosAmp[4][11] = CosFit->GetParameter(0);
  pCosAmpErr[4][11] = CosFit->GetParError(0);

  ///////////////////////////////////////////
  //////////////////CM6 p ///////////////////
  ///////////////////////////////////////////

  ParaPerpAsymmPhip_410MeVCM6 = Phip_410MeVCM6_Para->GetAsymmetry(Phip_410MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_410MeVCM6->SetName("ParaPerpAsymmPhip410MeVCM6");
  ParaPerpAsymmPhip_410MeVCM6->SetTitle("Proton Para/Perp Phi Asymmetry for 400-420MeV Photon Energy (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_410MeVCM6->Fit("CosFit", "LL");
  pCosAmp[5][0] = CosFit->GetParameter(0);
  pCosAmpErr[5][0] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_430MeVCM6 = Phip_430MeVCM6_Para->GetAsymmetry(Phip_430MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_430MeVCM6->SetName("ParaPerpAsymmPhip430MeVCM6");
  ParaPerpAsymmPhip_430MeVCM6->SetTitle("Proton Para/Perp Phi Asymmetry for 420-440MeV Photon Energy (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_430MeVCM6->Fit("CosFit", "LL");
  pCosAmp[5][1] = CosFit->GetParameter(0);
  pCosAmpErr[5][1] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_450MeVCM6 = Phip_450MeVCM6_Para->GetAsymmetry(Phip_450MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_450MeVCM6->SetName("ParaPerpAsymmPhip450MeVCM6");
  ParaPerpAsymmPhip_450MeVCM6->SetTitle("Proton Para/Perp Phi Asymmetry for 440-460MeV Photon Energy (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_450MeVCM6->Fit("CosFit", "LL");
  pCosAmp[5][2] = CosFit->GetParameter(0);
  pCosAmpErr[5][2] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_470MeVCM6 = Phip_470MeVCM6_Para->GetAsymmetry(Phip_470MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_470MeVCM6->SetName("ParaPerpAsymmPhip470MeVCM6");
  ParaPerpAsymmPhip_470MeVCM6->SetTitle("Proton Para/Perp Phi Asymmetry for 460-480MeV Photon Energy (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_470MeVCM6->Fit("CosFit", "LL");
  pCosAmp[5][3] = CosFit->GetParameter(0);
  pCosAmpErr[5][3] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_490MeVCM6 = Phip_490MeVCM6_Para->GetAsymmetry(Phip_490MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_490MeVCM6->SetName("ParaPerpAsymmPhip490MeVCM6");
  ParaPerpAsymmPhip_490MeVCM6->SetTitle("Proton Para/Perp Phi Asymmetry for 480-500MeV Photon Energy (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_490MeVCM6->Fit("CosFit", "LL");
  pCosAmp[5][4] = CosFit->GetParameter(0);
  pCosAmpErr[5][4] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_510MeVCM6 = Phip_510MeVCM6_Para->GetAsymmetry(Phip_510MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_510MeVCM6->SetName("ParaPerpAsymmPhip510MeVCM6");
  ParaPerpAsymmPhip_510MeVCM6->SetTitle("Proton Para/Perp Phi Asymmetry for 500-520MeV Photon Energy (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_510MeVCM6->Fit("CosFit", "LL");
  pCosAmp[5][5] = CosFit->GetParameter(0);
  pCosAmpErr[5][5] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_530MeVCM6 = Phip_530MeVCM6_Para->GetAsymmetry(Phip_530MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_530MeVCM6->SetName("ParaPerpAsymmPhip530MeVCM6");
  ParaPerpAsymmPhip_530MeVCM6->SetTitle("Proton Para/Perp Phi Asymmetry for 520-540MeV Photon Energy (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_530MeVCM6->Fit("CosFit", "LL");
  pCosAmp[5][6] = CosFit->GetParameter(0);
  pCosAmpErr[5][6] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_550MeVCM6 = Phip_550MeVCM6_Para->GetAsymmetry(Phip_550MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_550MeVCM6->SetName("ParaPerpAsymmPhip550MeVCM6");
  ParaPerpAsymmPhip_550MeVCM6->SetTitle("Proton Para/Perp Phi Asymmetry for 540-560MeV Photon Energy (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_550MeVCM6->Fit("CosFit", "LL");
  pCosAmp[5][7] = CosFit->GetParameter(0);
  pCosAmpErr[5][7] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_570MeVCM6 = Phip_570MeVCM6_Para->GetAsymmetry(Phip_570MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_570MeVCM6->SetName("ParaPerpAsymmPhip570MeVCM6");
  ParaPerpAsymmPhip_570MeVCM6->SetTitle("Proton Para/Perp Phi Asymmetry for 560-580MeV Photon Energy (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_570MeVCM6->Fit("CosFit", "LL");
  pCosAmp[5][8] = CosFit->GetParameter(0);
  pCosAmpErr[5][8] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_590MeVCM6 = Phip_590MeVCM6_Para->GetAsymmetry(Phip_590MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_590MeVCM6->SetName("ParaPerpAsymmPhip590MeVCM6");
  ParaPerpAsymmPhip_590MeVCM6->SetTitle("Proton Para/Perp Phi Asymmetry for 580-600MeV Photon Energy (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_590MeVCM6->Fit("CosFit", "LL");
  pCosAmp[5][9] = CosFit->GetParameter(0);
  pCosAmpErr[5][9] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_610MeVCM6 = Phip_610MeVCM6_Para->GetAsymmetry(Phip_610MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_610MeVCM6->SetName("ParaPerpAsymmPhip610MeVCM6");
  ParaPerpAsymmPhip_610MeVCM6->SetTitle("Proton Para/Perp Phi Asymmetry for 600-620MeV Photon Energy (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_610MeVCM6->Fit("CosFit", "LL");
  pCosAmp[5][10] = CosFit->GetParameter(0);
  pCosAmpErr[5][10] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_630MeVCM6 = Phip_630MeVCM6_Para->GetAsymmetry(Phip_630MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_630MeVCM6->SetName("ParaPerpAsymmPhip630MeVCM6");
  ParaPerpAsymmPhip_630MeVCM6->SetTitle("Proton Para/Perp Phi Asymmetry for 620-640MeV Photon Energy (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_630MeVCM6->Fit("CosFit", "LL");
  pCosAmp[5][11] = CosFit->GetParameter(0);
  pCosAmpErr[5][11] = CosFit->GetParError(0);

  ///////////////////////////////////////////
  //////////////////CM7 p ///////////////////
  ///////////////////////////////////////////

  ParaPerpAsymmPhip_410MeVCM7 = Phip_410MeVCM7_Para->GetAsymmetry(Phip_410MeVCM7_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_410MeVCM7->SetName("ParaPerpAsymmPhip410MeVCM7");
  ParaPerpAsymmPhip_410MeVCM7->SetTitle("Proton Para/Perp Phi Asymmetry for 400-420MeV Photon Energy (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_410MeVCM7->Fit("CosFit", "LL");
  pCosAmp[6][0] = CosFit->GetParameter(0);
  pCosAmpErr[6][0] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_430MeVCM7 = Phip_430MeVCM7_Para->GetAsymmetry(Phip_430MeVCM7_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_430MeVCM7->SetName("ParaPerpAsymmPhip430MeVCM7");
  ParaPerpAsymmPhip_430MeVCM7->SetTitle("Proton Para/Perp Phi Asymmetry for 420-440MeV Photon Energy (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_430MeVCM7->Fit("CosFit", "LL");
  pCosAmp[6][1] = CosFit->GetParameter(0);
  pCosAmpErr[6][1] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_450MeVCM7 = Phip_450MeVCM7_Para->GetAsymmetry(Phip_450MeVCM7_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_450MeVCM7->SetName("ParaPerpAsymmPhip450MeVCM7");
  ParaPerpAsymmPhip_450MeVCM7->SetTitle("Proton Para/Perp Phi Asymmetry for 440-460MeV Photon Energy (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_450MeVCM7->Fit("CosFit", "LL");
  pCosAmp[6][2] = CosFit->GetParameter(0);
  pCosAmpErr[6][2] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_470MeVCM7 = Phip_470MeVCM7_Para->GetAsymmetry(Phip_470MeVCM7_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_470MeVCM7->SetName("ParaPerpAsymmPhip470MeVCM7");
  ParaPerpAsymmPhip_470MeVCM7->SetTitle("Proton Para/Perp Phi Asymmetry for 460-480MeV Photon Energy (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_470MeVCM7->Fit("CosFit", "LL");
  pCosAmp[6][3] = CosFit->GetParameter(0);
  pCosAmpErr[6][3] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_490MeVCM7 = Phip_490MeVCM7_Para->GetAsymmetry(Phip_490MeVCM7_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_490MeVCM7->SetName("ParaPerpAsymmPhip490MeVCM7");
  ParaPerpAsymmPhip_490MeVCM7->SetTitle("Proton Para/Perp Phi Asymmetry for 480-500MeV Photon Energy (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_490MeVCM7->Fit("CosFit", "LL");
  pCosAmp[6][4] = CosFit->GetParameter(0);
  pCosAmpErr[6][4] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_510MeVCM7 = Phip_510MeVCM7_Para->GetAsymmetry(Phip_510MeVCM7_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_510MeVCM7->SetName("ParaPerpAsymmPhip510MeVCM7");
  ParaPerpAsymmPhip_510MeVCM7->SetTitle("Proton Para/Perp Phi Asymmetry for 500-520MeV Photon Energy (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_510MeVCM7->Fit("CosFit", "LL");
  pCosAmp[6][5] = CosFit->GetParameter(0);
  pCosAmpErr[6][5] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_530MeVCM7 = Phip_530MeVCM7_Para->GetAsymmetry(Phip_530MeVCM7_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_530MeVCM7->SetName("ParaPerpAsymmPhip530MeVCM7");
  ParaPerpAsymmPhip_530MeVCM7->SetTitle("Proton Para/Perp Phi Asymmetry for 520-540MeV Photon Energy (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_530MeVCM7->Fit("CosFit", "LL");
  pCosAmp[6][6] = CosFit->GetParameter(0);
  pCosAmpErr[6][6] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_550MeVCM7 = Phip_550MeVCM7_Para->GetAsymmetry(Phip_550MeVCM7_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_550MeVCM7->SetName("ParaPerpAsymmPhip550MeVCM7");
  ParaPerpAsymmPhip_550MeVCM7->SetTitle("Proton Para/Perp Phi Asymmetry for 540-560MeV Photon Energy (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_550MeVCM7->Fit("CosFit", "LL");
  pCosAmp[6][7] = CosFit->GetParameter(0);
  pCosAmpErr[6][7] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_570MeVCM7 = Phip_570MeVCM7_Para->GetAsymmetry(Phip_570MeVCM7_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_570MeVCM7->SetName("ParaPerpAsymmPhip570MeVCM7");
  ParaPerpAsymmPhip_570MeVCM7->SetTitle("Proton Para/Perp Phi Asymmetry for 560-580MeV Photon Energy (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_570MeVCM7->Fit("CosFit", "LL");
  pCosAmp[6][8] = CosFit->GetParameter(0);
  pCosAmpErr[6][8] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_590MeVCM7 = Phip_590MeVCM7_Para->GetAsymmetry(Phip_590MeVCM7_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_590MeVCM7->SetName("ParaPerpAsymmPhip590MeVCM7");
  ParaPerpAsymmPhip_590MeVCM7->SetTitle("Proton Para/Perp Phi Asymmetry for 580-600MeV Photon Energy (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_590MeVCM7->Fit("CosFit", "LL");
  pCosAmp[6][9] = CosFit->GetParameter(0);
  pCosAmpErr[6][9] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_610MeVCM7 = Phip_610MeVCM7_Para->GetAsymmetry(Phip_610MeVCM7_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_610MeVCM7->SetName("ParaPerpAsymmPhip610MeVCM7");
  ParaPerpAsymmPhip_610MeVCM7->SetTitle("Proton Para/Perp Phi Asymmetry for 600-620MeV Photon Energy (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_610MeVCM7->Fit("CosFit", "LL");
  pCosAmp[6][10] = CosFit->GetParameter(0);
  pCosAmpErr[6][10] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_630MeVCM7 = Phip_630MeVCM7_Para->GetAsymmetry(Phip_630MeVCM7_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_630MeVCM7->SetName("ParaPerpAsymmPhip630MeVCM7");
  ParaPerpAsymmPhip_630MeVCM7->SetTitle("Proton Para/Perp Phi Asymmetry for 620-640MeV Photon Energy (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_630MeVCM7->Fit("CosFit", "LL");
  pCosAmp[6][11] = CosFit->GetParameter(0);
  pCosAmpErr[6][11] = CosFit->GetParError(0);

  ///////////////////////////////////////////
  //////////////////CM8 p ///////////////////
  ///////////////////////////////////////////

  ParaPerpAsymmPhip_410MeVCM8 = Phip_410MeVCM8_Para->GetAsymmetry(Phip_410MeVCM8_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_410MeVCM8->SetName("ParaPerpAsymmPhip410MeVCM8");
  ParaPerpAsymmPhip_410MeVCM8->SetTitle("Proton Para/Perp Phi Asymmetry for 400-420MeV Photon Energy (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_410MeVCM8->Fit("CosFit", "LL");
  pCosAmp[7][0] = CosFit->GetParameter(0);
  pCosAmpErr[7][0] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_430MeVCM8 = Phip_430MeVCM8_Para->GetAsymmetry(Phip_430MeVCM8_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_430MeVCM8->SetName("ParaPerpAsymmPhip430MeVCM8");
  ParaPerpAsymmPhip_430MeVCM8->SetTitle("Proton Para/Perp Phi Asymmetry for 420-440MeV Photon Energy (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_430MeVCM8->Fit("CosFit", "LL");
  pCosAmp[7][1] = CosFit->GetParameter(0);
  pCosAmpErr[7][1] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_450MeVCM8 = Phip_450MeVCM8_Para->GetAsymmetry(Phip_450MeVCM8_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_450MeVCM8->SetName("ParaPerpAsymmPhip450MeVCM8");
  ParaPerpAsymmPhip_450MeVCM8->SetTitle("Proton Para/Perp Phi Asymmetry for 440-460MeV Photon Energy (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_450MeVCM8->Fit("CosFit", "LL");
  pCosAmp[7][2] = CosFit->GetParameter(0);
  pCosAmpErr[7][2] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_470MeVCM8 = Phip_470MeVCM8_Para->GetAsymmetry(Phip_470MeVCM8_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_470MeVCM8->SetName("ParaPerpAsymmPhip470MeVCM8");
  ParaPerpAsymmPhip_470MeVCM8->SetTitle("Proton Para/Perp Phi Asymmetry for 460-480MeV Photon Energy (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_470MeVCM8->Fit("CosFit", "LL");
  pCosAmp[7][3] = CosFit->GetParameter(0);
  pCosAmpErr[7][3] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_490MeVCM8 = Phip_490MeVCM8_Para->GetAsymmetry(Phip_490MeVCM8_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_490MeVCM8->SetName("ParaPerpAsymmPhip490MeVCM8");
  ParaPerpAsymmPhip_490MeVCM8->SetTitle("Proton Para/Perp Phi Asymmetry for 480-500MeV Photon Energy (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_490MeVCM8->Fit("CosFit", "LL");
  pCosAmp[7][4] = CosFit->GetParameter(0);
  pCosAmpErr[7][4] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_510MeVCM8 = Phip_510MeVCM8_Para->GetAsymmetry(Phip_510MeVCM8_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_510MeVCM8->SetName("ParaPerpAsymmPhip510MeVCM8");
  ParaPerpAsymmPhip_510MeVCM8->SetTitle("Proton Para/Perp Phi Asymmetry for 500-520MeV Photon Energy (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_510MeVCM8->Fit("CosFit", "LL");
  pCosAmp[7][5] = CosFit->GetParameter(0);
  pCosAmpErr[7][5] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_530MeVCM8 = Phip_530MeVCM8_Para->GetAsymmetry(Phip_530MeVCM8_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_530MeVCM8->SetName("ParaPerpAsymmPhip530MeVCM8");
  ParaPerpAsymmPhip_530MeVCM8->SetTitle("Proton Para/Perp Phi Asymmetry for 520-540MeV Photon Energy (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_530MeVCM8->Fit("CosFit", "LL");
  pCosAmp[7][6] = CosFit->GetParameter(0);
  pCosAmpErr[7][6] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_550MeVCM8 = Phip_550MeVCM8_Para->GetAsymmetry(Phip_550MeVCM8_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_550MeVCM8->SetName("ParaPerpAsymmPhip550MeVCM8");
  ParaPerpAsymmPhip_550MeVCM8->SetTitle("Proton Para/Perp Phi Asymmetry for 540-560MeV Photon Energy (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_550MeVCM8->Fit("CosFit", "LL");
  pCosAmp[7][7] = CosFit->GetParameter(0);
  pCosAmpErr[7][7] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_570MeVCM8 = Phip_570MeVCM8_Para->GetAsymmetry(Phip_570MeVCM8_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_570MeVCM8->SetName("ParaPerpAsymmPhip570MeVCM8");
  ParaPerpAsymmPhip_570MeVCM8->SetTitle("Proton Para/Perp Phi Asymmetry for 560-580MeV Photon Energy (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_570MeVCM8->Fit("CosFit", "LL");
  pCosAmp[7][8] = CosFit->GetParameter(0);
  pCosAmpErr[7][8] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_590MeVCM8 = Phip_590MeVCM8_Para->GetAsymmetry(Phip_590MeVCM8_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_590MeVCM8->SetName("ParaPerpAsymmPhip590MeVCM8");
  ParaPerpAsymmPhip_590MeVCM8->SetTitle("Proton Para/Perp Phi Asymmetry for 580-600MeV Photon Energy (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_590MeVCM8->Fit("CosFit", "LL");
  pCosAmp[7][9] = CosFit->GetParameter(0);
  pCosAmpErr[7][9] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_610MeVCM8 = Phip_610MeVCM8_Para->GetAsymmetry(Phip_610MeVCM8_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_610MeVCM8->SetName("ParaPerpAsymmPhip610MeVCM8");
  ParaPerpAsymmPhip_610MeVCM8->SetTitle("Proton Para/Perp Phi Asymmetry for 600-620MeV Photon Energy (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_610MeVCM8->Fit("CosFit", "LL");
  pCosAmp[7][10] = CosFit->GetParameter(0);
  pCosAmpErr[7][10] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_630MeVCM8 = Phip_630MeVCM8_Para->GetAsymmetry(Phip_630MeVCM8_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_630MeVCM8->SetName("ParaPerpAsymmPhip630MeVCM8");
  ParaPerpAsymmPhip_630MeVCM8->SetTitle("Proton Para/Perp Phi Asymmetry for 620-640MeV Photon Energy (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_630MeVCM8->Fit("CosFit", "LL");
  pCosAmp[7][11] = CosFit->GetParameter(0);
  pCosAmpErr[7][11] = CosFit->GetParError(0);

  ///////////////////////////////////////////
  //////////////////CM9 p ///////////////////
  ///////////////////////////////////////////

  ParaPerpAsymmPhip_410MeVCM9 = Phip_410MeVCM9_Para->GetAsymmetry(Phip_410MeVCM9_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_410MeVCM9->SetName("ParaPerpAsymmPhip410MeVCM9");
  ParaPerpAsymmPhip_410MeVCM9->SetTitle("Proton Para/Perp Phi Asymmetry for 400-420MeV Photon Energy (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_410MeVCM9->Fit("CosFit", "LL");
  pCosAmp[8][0] = CosFit->GetParameter(0);
  pCosAmpErr[8][0] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_430MeVCM9 = Phip_430MeVCM9_Para->GetAsymmetry(Phip_430MeVCM9_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_430MeVCM9->SetName("ParaPerpAsymmPhip430MeVCM9");
  ParaPerpAsymmPhip_430MeVCM9->SetTitle("Proton Para/Perp Phi Asymmetry for 420-440MeV Photon Energy (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_430MeVCM9->Fit("CosFit", "LL");
  pCosAmp[8][1] = CosFit->GetParameter(0);
  pCosAmpErr[8][1] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_450MeVCM9 = Phip_450MeVCM9_Para->GetAsymmetry(Phip_450MeVCM9_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_450MeVCM9->SetName("ParaPerpAsymmPhip450MeVCM9");
  ParaPerpAsymmPhip_450MeVCM9->SetTitle("Proton Para/Perp Phi Asymmetry for 440-460MeV Photon Energy (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_450MeVCM9->Fit("CosFit", "LL");
  pCosAmp[8][2] = CosFit->GetParameter(0);
  pCosAmpErr[8][2] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_470MeVCM9 = Phip_470MeVCM9_Para->GetAsymmetry(Phip_470MeVCM9_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_470MeVCM9->SetName("ParaPerpAsymmPhip470MeVCM9");
  ParaPerpAsymmPhip_470MeVCM9->SetTitle("Proton Para/Perp Phi Asymmetry for 460-480MeV Photon Energy (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_470MeVCM9->Fit("CosFit", "LL");
  pCosAmp[8][3] = CosFit->GetParameter(0);
  pCosAmpErr[8][3] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_490MeVCM9 = Phip_490MeVCM9_Para->GetAsymmetry(Phip_490MeVCM9_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_490MeVCM9->SetName("ParaPerpAsymmPhip490MeVCM9");
  ParaPerpAsymmPhip_490MeVCM9->SetTitle("Proton Para/Perp Phi Asymmetry for 480-500MeV Photon Energy (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_490MeVCM9->Fit("CosFit", "LL");
  pCosAmp[8][4] = CosFit->GetParameter(0);
  pCosAmpErr[8][4] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_510MeVCM9 = Phip_510MeVCM9_Para->GetAsymmetry(Phip_510MeVCM9_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_510MeVCM9->SetName("ParaPerpAsymmPhip510MeVCM9");
  ParaPerpAsymmPhip_510MeVCM9->SetTitle("Proton Para/Perp Phi Asymmetry for 500-520MeV Photon Energy (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_510MeVCM9->Fit("CosFit", "LL");
  pCosAmp[8][5] = CosFit->GetParameter(0);
  pCosAmpErr[8][5] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_530MeVCM9 = Phip_530MeVCM9_Para->GetAsymmetry(Phip_530MeVCM9_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_530MeVCM9->SetName("ParaPerpAsymmPhip530MeVCM9");
  ParaPerpAsymmPhip_530MeVCM9->SetTitle("Proton Para/Perp Phi Asymmetry for 520-540MeV Photon Energy (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_530MeVCM9->Fit("CosFit", "LL");
  pCosAmp[8][6] = CosFit->GetParameter(0);
  pCosAmpErr[8][6] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_550MeVCM9 = Phip_550MeVCM9_Para->GetAsymmetry(Phip_550MeVCM9_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_550MeVCM9->SetName("ParaPerpAsymmPhip550MeVCM9");
  ParaPerpAsymmPhip_550MeVCM9->SetTitle("Proton Para/Perp Phi Asymmetry for 540-560MeV Photon Energy (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_550MeVCM9->Fit("CosFit", "LL");
  pCosAmp[8][7] = CosFit->GetParameter(0);
  pCosAmpErr[8][7] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_570MeVCM9 = Phip_570MeVCM9_Para->GetAsymmetry(Phip_570MeVCM9_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_570MeVCM9->SetName("ParaPerpAsymmPhip570MeVCM9");
  ParaPerpAsymmPhip_570MeVCM9->SetTitle("Proton Para/Perp Phi Asymmetry for 560-580MeV Photon Energy (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_570MeVCM9->Fit("CosFit", "LL");
  pCosAmp[8][8] = CosFit->GetParameter(0);
  pCosAmpErr[8][8] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_590MeVCM9 = Phip_590MeVCM9_Para->GetAsymmetry(Phip_590MeVCM9_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_590MeVCM9->SetName("ParaPerpAsymmPhip590MeVCM9");
  ParaPerpAsymmPhip_590MeVCM9->SetTitle("Proton Para/Perp Phi Asymmetry for 580-600MeV Photon Energy (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_590MeVCM9->Fit("CosFit", "LL");
  pCosAmp[8][9] = CosFit->GetParameter(0);
  pCosAmpErr[8][9] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_610MeVCM9 = Phip_610MeVCM9_Para->GetAsymmetry(Phip_610MeVCM9_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_610MeVCM9->SetName("ParaPerpAsymmPhip610MeVCM9");
  ParaPerpAsymmPhip_610MeVCM9->SetTitle("Proton Para/Perp Phi Asymmetry for 600-620MeV Photon Energy (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_610MeVCM9->Fit("CosFit", "LL");
  pCosAmp[8][10] = CosFit->GetParameter(0);
  pCosAmpErr[8][10] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_630MeVCM9 = Phip_630MeVCM9_Para->GetAsymmetry(Phip_630MeVCM9_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_630MeVCM9->SetName("ParaPerpAsymmPhip630MeVCM9");
  ParaPerpAsymmPhip_630MeVCM9->SetTitle("Proton Para/Perp Phi Asymmetry for 620-640MeV Photon Energy (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_630MeVCM9->Fit("CosFit", "LL");
  pCosAmp[8][11] = CosFit->GetParameter(0);
  pCosAmpErr[8][11] = CosFit->GetParError(0);

  ///////////////////////////////////////////
  //////////////////CM10 p ///////////////////
  ///////////////////////////////////////////

  ParaPerpAsymmPhip_410MeVCM10 = Phip_410MeVCM10_Para->GetAsymmetry(Phip_410MeVCM10_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_410MeVCM10->SetName("ParaPerpAsymmPhip410MeVCM10");
  ParaPerpAsymmPhip_410MeVCM10->SetTitle("Proton Para/Perp Phi Asymmetry for 400-420MeV Photon Energy (CosTheta-0.8-(-1.))");
  ParaPerpAsymmPhip_410MeVCM10->Fit("CosFit", "LL");
  pCosAmp[9][0] = CosFit->GetParameter(0);
  pCosAmpErr[9][0] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_430MeVCM10 = Phip_430MeVCM10_Para->GetAsymmetry(Phip_430MeVCM10_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_430MeVCM10->SetName("ParaPerpAsymmPhip430MeVCM10");
  ParaPerpAsymmPhip_430MeVCM10->SetTitle("Proton Para/Perp Phi Asymmetry for 420-440MeV Photon Energy (CosTheta-0.8-(-1.))");
  ParaPerpAsymmPhip_430MeVCM10->Fit("CosFit", "LL");
  pCosAmp[9][1] = CosFit->GetParameter(0);
  pCosAmpErr[9][1] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_450MeVCM10 = Phip_450MeVCM10_Para->GetAsymmetry(Phip_450MeVCM10_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_450MeVCM10->SetName("ParaPerpAsymmPhip450MeVCM10");
  ParaPerpAsymmPhip_450MeVCM10->SetTitle("Proton Para/Perp Phi Asymmetry for 440-460MeV Photon Energy (CosTheta-0.8-(-1.))");
  ParaPerpAsymmPhip_450MeVCM10->Fit("CosFit", "LL");
  pCosAmp[9][2] = CosFit->GetParameter(0);
  pCosAmpErr[9][2] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_470MeVCM10 = Phip_470MeVCM10_Para->GetAsymmetry(Phip_470MeVCM10_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_470MeVCM10->SetName("ParaPerpAsymmPhip470MeVCM10");
  ParaPerpAsymmPhip_470MeVCM10->SetTitle("Proton Para/Perp Phi Asymmetry for 460-480MeV Photon Energy (CosTheta-0.8-(-1.))");
  ParaPerpAsymmPhip_470MeVCM10->Fit("CosFit", "LL");
  pCosAmp[9][3] = CosFit->GetParameter(0);
  pCosAmpErr[9][3] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_490MeVCM10 = Phip_490MeVCM10_Para->GetAsymmetry(Phip_490MeVCM10_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_490MeVCM10->SetName("ParaPerpAsymmPhip490MeVCM10");
  ParaPerpAsymmPhip_490MeVCM10->SetTitle("Proton Para/Perp Phi Asymmetry for 480-500MeV Photon Energy (CosTheta-0.8-(-1.))");
  ParaPerpAsymmPhip_490MeVCM10->Fit("CosFit", "LL");
  pCosAmp[9][4] = CosFit->GetParameter(0);
  pCosAmpErr[9][4] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_510MeVCM10 = Phip_510MeVCM10_Para->GetAsymmetry(Phip_510MeVCM10_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_510MeVCM10->SetName("ParaPerpAsymmPhip510MeVCM10");
  ParaPerpAsymmPhip_510MeVCM10->SetTitle("Proton Para/Perp Phi Asymmetry for 500-520MeV Photon Energy (CosTheta-0.8-(-1.))");
  ParaPerpAsymmPhip_510MeVCM10->Fit("CosFit", "LL");
  pCosAmp[9][5] = CosFit->GetParameter(0);
  pCosAmpErr[9][5] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_530MeVCM10 = Phip_530MeVCM10_Para->GetAsymmetry(Phip_530MeVCM10_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_530MeVCM10->SetName("ParaPerpAsymmPhip530MeVCM10");
  ParaPerpAsymmPhip_530MeVCM10->SetTitle("Proton Para/Perp Phi Asymmetry for 520-540MeV Photon Energy (CosTheta-0.8-(-1.))");
  ParaPerpAsymmPhip_530MeVCM10->Fit("CosFit", "LL");
  pCosAmp[9][6] = CosFit->GetParameter(0);
  pCosAmpErr[9][6] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_550MeVCM10 = Phip_550MeVCM10_Para->GetAsymmetry(Phip_550MeVCM10_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_550MeVCM10->SetName("ParaPerpAsymmPhip550MeVCM10");
  ParaPerpAsymmPhip_550MeVCM10->SetTitle("Proton Para/Perp Phi Asymmetry for 540-560MeV Photon Energy (CosTheta-0.8-(-1.))");
  ParaPerpAsymmPhip_550MeVCM10->Fit("CosFit", "LL");
  pCosAmp[9][7] = CosFit->GetParameter(0);
  pCosAmpErr[9][7] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_570MeVCM10 = Phip_570MeVCM10_Para->GetAsymmetry(Phip_570MeVCM10_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_570MeVCM10->SetName("ParaPerpAsymmPhip570MeVCM10");
  ParaPerpAsymmPhip_570MeVCM10->SetTitle("Proton Para/Perp Phi Asymmetry for 560-580MeV Photon Energy (CosTheta-0.8-(-1.))");
  ParaPerpAsymmPhip_570MeVCM10->Fit("CosFit", "LL");
  pCosAmp[9][8] = CosFit->GetParameter(0);
  pCosAmpErr[9][8] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_590MeVCM10 = Phip_590MeVCM10_Para->GetAsymmetry(Phip_590MeVCM10_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_590MeVCM10->SetName("ParaPerpAsymmPhip590MeVCM10");
  ParaPerpAsymmPhip_590MeVCM10->SetTitle("Proton Para/Perp Phi Asymmetry for 580-600MeV Photon Energy (CosTheta-0.8-(-1.))");
  ParaPerpAsymmPhip_590MeVCM10->Fit("CosFit", "LL");
  pCosAmp[9][9] = CosFit->GetParameter(0);
  pCosAmpErr[9][9] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_610MeVCM10 = Phip_610MeVCM10_Para->GetAsymmetry(Phip_610MeVCM10_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_610MeVCM10->SetName("ParaPerpAsymmPhip610MeVCM10");
  ParaPerpAsymmPhip_610MeVCM10->SetTitle("Proton Para/Perp Phi Asymmetry for 600-620MeV Photon Energy (CosTheta-0.8-(-1.))");
  ParaPerpAsymmPhip_610MeVCM10->Fit("CosFit", "LL");
  pCosAmp[9][10] = CosFit->GetParameter(0);
  pCosAmpErr[9][10] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_630MeVCM10 = Phip_630MeVCM10_Para->GetAsymmetry(Phip_630MeVCM10_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_630MeVCM10->SetName("ParaPerpAsymmPhip630MeVCM10");
  ParaPerpAsymmPhip_630MeVCM10->SetTitle("Proton Para/Perp Phi Asymmetry for 620-640MeV Photon Energy (CosTheta-0.8-(-1.))");
  ParaPerpAsymmPhip_630MeVCM10->Fit("CosFit", "LL");
  pCosAmp[9][11] = CosFit->GetParameter(0);
  pCosAmpErr[9][11] = CosFit->GetParError(0);

  TFile f1("ParaPerpAsymm_Total_10_v2.root", "RECREATE");

  ParaPerpAsymmPhip_410MeVCM1->Write();
  ParaPerpAsymmPhip_430MeVCM1->Write();
  ParaPerpAsymmPhip_450MeVCM1->Write();
  ParaPerpAsymmPhip_470MeVCM1->Write();
  ParaPerpAsymmPhip_490MeVCM1->Write();
  ParaPerpAsymmPhip_510MeVCM1->Write();
  ParaPerpAsymmPhip_530MeVCM1->Write();
  ParaPerpAsymmPhip_550MeVCM1->Write();
  ParaPerpAsymmPhip_570MeVCM1->Write();
  ParaPerpAsymmPhip_590MeVCM1->Write();
  ParaPerpAsymmPhip_610MeVCM1->Write();
  ParaPerpAsymmPhip_630MeVCM1->Write();

  ParaPerpAsymmPhip_410MeVCM2->Write();
  ParaPerpAsymmPhip_430MeVCM2->Write();
  ParaPerpAsymmPhip_450MeVCM2->Write();
  ParaPerpAsymmPhip_470MeVCM2->Write();
  ParaPerpAsymmPhip_490MeVCM2->Write();
  ParaPerpAsymmPhip_510MeVCM2->Write();
  ParaPerpAsymmPhip_530MeVCM2->Write();
  ParaPerpAsymmPhip_550MeVCM2->Write();
  ParaPerpAsymmPhip_570MeVCM2->Write();
  ParaPerpAsymmPhip_590MeVCM2->Write();
  ParaPerpAsymmPhip_610MeVCM2->Write();
  ParaPerpAsymmPhip_630MeVCM2->Write();

  ParaPerpAsymmPhip_410MeVCM3->Write();
  ParaPerpAsymmPhip_430MeVCM3->Write();
  ParaPerpAsymmPhip_450MeVCM3->Write();
  ParaPerpAsymmPhip_470MeVCM3->Write();
  ParaPerpAsymmPhip_490MeVCM3->Write();
  ParaPerpAsymmPhip_510MeVCM3->Write();
  ParaPerpAsymmPhip_530MeVCM3->Write();
  ParaPerpAsymmPhip_550MeVCM3->Write();
  ParaPerpAsymmPhip_570MeVCM3->Write();
  ParaPerpAsymmPhip_590MeVCM3->Write();
  ParaPerpAsymmPhip_610MeVCM3->Write();
  ParaPerpAsymmPhip_630MeVCM3->Write();

  ParaPerpAsymmPhip_410MeVCM4->Write();
  ParaPerpAsymmPhip_430MeVCM4->Write();
  ParaPerpAsymmPhip_450MeVCM4->Write();
  ParaPerpAsymmPhip_470MeVCM4->Write();
  ParaPerpAsymmPhip_490MeVCM4->Write();
  ParaPerpAsymmPhip_510MeVCM4->Write();
  ParaPerpAsymmPhip_530MeVCM4->Write();
  ParaPerpAsymmPhip_550MeVCM4->Write();
  ParaPerpAsymmPhip_570MeVCM4->Write();
  ParaPerpAsymmPhip_590MeVCM4->Write();
  ParaPerpAsymmPhip_610MeVCM4->Write();
  ParaPerpAsymmPhip_630MeVCM4->Write();

  ParaPerpAsymmPhip_410MeVCM5->Write();
  ParaPerpAsymmPhip_430MeVCM5->Write();
  ParaPerpAsymmPhip_450MeVCM5->Write();
  ParaPerpAsymmPhip_470MeVCM5->Write();
  ParaPerpAsymmPhip_490MeVCM5->Write();
  ParaPerpAsymmPhip_510MeVCM5->Write();
  ParaPerpAsymmPhip_530MeVCM5->Write();
  ParaPerpAsymmPhip_550MeVCM5->Write();
  ParaPerpAsymmPhip_570MeVCM5->Write();
  ParaPerpAsymmPhip_590MeVCM5->Write();
  ParaPerpAsymmPhip_610MeVCM5->Write();
  ParaPerpAsymmPhip_630MeVCM5->Write();

  ParaPerpAsymmPhip_410MeVCM6->Write();
  ParaPerpAsymmPhip_430MeVCM6->Write();
  ParaPerpAsymmPhip_450MeVCM6->Write();
  ParaPerpAsymmPhip_470MeVCM6->Write();
  ParaPerpAsymmPhip_490MeVCM6->Write();
  ParaPerpAsymmPhip_510MeVCM6->Write();
  ParaPerpAsymmPhip_530MeVCM6->Write();
  ParaPerpAsymmPhip_550MeVCM6->Write();
  ParaPerpAsymmPhip_570MeVCM6->Write();
  ParaPerpAsymmPhip_590MeVCM6->Write();
  ParaPerpAsymmPhip_610MeVCM6->Write();
  ParaPerpAsymmPhip_630MeVCM6->Write();

  ParaPerpAsymmPhip_410MeVCM7->Write();
  ParaPerpAsymmPhip_430MeVCM7->Write();
  ParaPerpAsymmPhip_450MeVCM7->Write();
  ParaPerpAsymmPhip_470MeVCM7->Write();
  ParaPerpAsymmPhip_490MeVCM7->Write();
  ParaPerpAsymmPhip_510MeVCM7->Write();
  ParaPerpAsymmPhip_530MeVCM7->Write();
  ParaPerpAsymmPhip_550MeVCM7->Write();
  ParaPerpAsymmPhip_570MeVCM7->Write();
  ParaPerpAsymmPhip_590MeVCM7->Write();
  ParaPerpAsymmPhip_610MeVCM7->Write();
  ParaPerpAsymmPhip_630MeVCM7->Write();

  ParaPerpAsymmPhip_410MeVCM8->Write();
  ParaPerpAsymmPhip_430MeVCM8->Write();
  ParaPerpAsymmPhip_450MeVCM8->Write();
  ParaPerpAsymmPhip_470MeVCM8->Write();
  ParaPerpAsymmPhip_490MeVCM8->Write();
  ParaPerpAsymmPhip_510MeVCM8->Write();
  ParaPerpAsymmPhip_530MeVCM8->Write();
  ParaPerpAsymmPhip_550MeVCM8->Write();
  ParaPerpAsymmPhip_570MeVCM8->Write();
  ParaPerpAsymmPhip_590MeVCM8->Write();
  ParaPerpAsymmPhip_610MeVCM8->Write();
  ParaPerpAsymmPhip_630MeVCM8->Write();

  ParaPerpAsymmPhip_410MeVCM9->Write();
  ParaPerpAsymmPhip_430MeVCM9->Write();
  ParaPerpAsymmPhip_450MeVCM9->Write();
  ParaPerpAsymmPhip_470MeVCM9->Write();
  ParaPerpAsymmPhip_490MeVCM9->Write();
  ParaPerpAsymmPhip_510MeVCM9->Write();
  ParaPerpAsymmPhip_530MeVCM9->Write();
  ParaPerpAsymmPhip_550MeVCM9->Write();
  ParaPerpAsymmPhip_570MeVCM9->Write();
  ParaPerpAsymmPhip_590MeVCM9->Write();
  ParaPerpAsymmPhip_610MeVCM9->Write();
  ParaPerpAsymmPhip_630MeVCM9->Write();

  ParaPerpAsymmPhip_410MeVCM10->Write();
  ParaPerpAsymmPhip_430MeVCM10->Write();
  ParaPerpAsymmPhip_450MeVCM10->Write();
  ParaPerpAsymmPhip_470MeVCM10->Write();
  ParaPerpAsymmPhip_490MeVCM10->Write();
  ParaPerpAsymmPhip_510MeVCM10->Write();
  ParaPerpAsymmPhip_530MeVCM10->Write();
  ParaPerpAsymmPhip_550MeVCM10->Write();
  ParaPerpAsymmPhip_570MeVCM10->Write();
  ParaPerpAsymmPhip_590MeVCM10->Write();
  ParaPerpAsymmPhip_610MeVCM10->Write();
  ParaPerpAsymmPhip_630MeVCM10->Write();


  //Define new tree to store parameters in
  TTree* tree = new TTree("Parameter_Values", "Tree_of_Values");

  // Define branches to store parameters, (Branch Name, Variable, Type of Variable)
  tree->Branch("pCosAmp410", &pCosA410, "pCosA410/D");
  tree->Branch("pCosAmpErr410", &pCosAErr410, "pCosAErr410/D");
  tree->Branch("pCosAmp430", &pCosA430, "pCosA430/D");
  tree->Branch("pCosAmpErr430", &pCosAErr430, "pCosAErr430/D");
  tree->Branch("pCosAmp450", &pCosA450, "pCosA450/D");
  tree->Branch("pCosAmpErr450", &pCosAErr450, "pCosAErr450/D");
  tree->Branch("pCosAmp470", &pCosA470, "pCosA470/D");
  tree->Branch("pCosAmpErr470", &pCosAErr470, "pCosAErr470/D");
  tree->Branch("pCosAmp490", &pCosA490, "pCosA490/D");
  tree->Branch("pCosAmpErr490", &pCosAErr490, "pCosAErr490/D");
  tree->Branch("pCosAmp510", &pCosA510, "pCosA510/D");
  tree->Branch("pCosAmpErr510", &pCosAErr510, "pCosAErr510/D");
  tree->Branch("pCosAmp530", &pCosA530, "pCosA530/D");
  tree->Branch("pCosAmpErr530", &pCosAErr530, "pCosAErr530/D");
  tree->Branch("pCosAmp550", &pCosA550, "pCosA550/D");
  tree->Branch("pCosAmpErr550", &pCosAErr550, "pCosAErr550/D");
  tree->Branch("pCosAmp570", &pCosA570, "pCosA570/D");
  tree->Branch("pCosAmpErr570", &pCosAErr570, "pCosAErr570/D");
  tree->Branch("pCosAmp590", &pCosA590, "pCosA590/D");
  tree->Branch("pCosAmpErr590", &pCosAErr590, "pCosAErr590/D");
  tree->Branch("pCosAmp610", &pCosA610, "pCosA610/D");
  tree->Branch("pCosAmpErr610", &pCosAErr610, "pCosAErr610/D");
  tree->Branch("pCosAmp630", &pCosA630, "pCosA630/D");
  tree->Branch("pCosAmpErr630", &pCosAErr630, "pCosAErr630/D");

  // Fill branches (and hence tree) with corresponding parameters from above
  for (Int_t m = 0; m < 10; m++){
    pCosA410 = pCosAmp[m][0];
    pCosAErr410 = pCosAmpErr[m][0];
    pCosA430 = pCosAmp[m][1];
    pCosAErr430 = pCosAmpErr[m][1];
    pCosA450 = pCosAmp[m][2];
    pCosAErr450 = pCosAmpErr[m][2];
    pCosA470 = pCosAmp[m][3];
    pCosAErr470 = pCosAmpErr[m][3];
    pCosA490 = pCosAmp[m][4];
    pCosAErr490 = pCosAmpErr[m][4];
    pCosA510 = pCosAmp[m][5];
    pCosAErr510 = pCosAmpErr[m][5];
    pCosA530 = pCosAmp[m][6];
    pCosAErr530 = pCosAmpErr[m][6];
    pCosA550 = pCosAmp[m][7];
    pCosAErr550 = pCosAmpErr[m][7];
    pCosA570 = pCosAmp[m][8];
    pCosAErr570= pCosAmpErr[m][8];
    pCosA590 = pCosAmp[m][9];
    pCosAErr590= pCosAmpErr[m][9];
    pCosA610 = pCosAmp[m][10];
    pCosAErr610= pCosAmpErr[m][10];
    pCosA630 = pCosAmp[m][11];
    pCosAErr630= pCosAmpErr[m][11];
    tree->Fill();
  }

  f1.Write();

}
