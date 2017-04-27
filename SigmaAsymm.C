#include "./includes_SigmaAsymm.h"

void SigmaAsymm(){

  TF1 *CosFunc = new TF1("CosFit", "[0]*cos(2*x*TMath::DegToRad())", -130, 130);
  CosFunc->SetParNames("Amplitude");

  double pCosAmp[6][12]; // Format of array is Theta bin (x) by Egamma bin (y), 6 theta bins of 30, 12 20MeV Egamma bins
  double pCosAmpErr[6][12];
  double nCosAmp[6][12]; // Format of array is Theta bin (x) by Egamma bin (y), 6 theta bins of 30, 12 20MeV Egamma bins
  double nCosAmpErr[6][12];
  double pCosA;
  double pCosAErr;
  double nCosA;
  double nCosAErr;

  TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/ParaPerp_Total_9_Combined.root"); // Open the latest PTotal combined file to load histograms from
  NPara = Eg_Para->GetEntries();
  NPerp = Eg_Perp->GetEntries();
  ScaleFactor = NPara/NPerp;
  ScaleFactorErr = sqrt( (NPara/((TMath::Power(NPerp,2)))) + (((TMath::Power(NPara,2)))/(TMath::Power(NPerp,3))) ); // Error Propagation of above formula, see notebook 5
  
  ParaPerpAsymmPhip_410MeVCM1 = Phip_410MeVCM1_Para->GetAsymmetry(Phip_410MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_410MeVCM1->SetName("ParaPerpAsymmPhip410MeVCM1");
  ParaPerpAsymmPhip_410MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for 400-420MeV Photon Energy (ThetaCM 0-30)");
  ParaPerpAsymmPhip_410MeVCM1->Fit("CosFit", "LL");
  pCosAmp[0][0] = CosFit->GetParameter(0);
  pCosAmpErr[0][0] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_430MeVCM1 = Phip_430MeVCM1_Para->GetAsymmetry(Phip_430MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_430MeVCM1->SetName("ParaPerpAsymmPhip430MeVCM1");
  ParaPerpAsymmPhip_430MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for 420-440MeV Photon Energy (ThetaCM 0-30)");
  ParaPerpAsymmPhip_430MeVCM1->Fit("CosFit", "LL");
  pCosAmp[0][1] = CosFit->GetParameter(0);
  pCosAmpErr[0][1] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_450MeVCM1 = Phip_450MeVCM1_Para->GetAsymmetry(Phip_450MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_450MeVCM1->SetName("ParaPerpAsymmPhip450MeVCM1");
  ParaPerpAsymmPhip_450MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for 440-460MeV Photon Energy (ThetaCM 0-30)");
  ParaPerpAsymmPhip_450MeVCM1->Fit("CosFit", "LL");
  pCosAmp[0][2] = CosFit->GetParameter(0);
  pCosAmpErr[0][2] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_470MeVCM1 = Phip_470MeVCM1_Para->GetAsymmetry(Phip_470MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_470MeVCM1->SetName("ParaPerpAsymmPhip470MeVCM1");
  ParaPerpAsymmPhip_470MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for 460-480MeV Photon Energy (ThetaCM 0-30)");
  ParaPerpAsymmPhip_470MeVCM1->Fit("CosFit", "LL");
  pCosAmp[0][3] = CosFit->GetParameter(0);
  pCosAmpErr[0][3] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_490MeVCM1 = Phip_490MeVCM1_Para->GetAsymmetry(Phip_490MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_490MeVCM1->SetName("ParaPerpAsymmPhip490MeVCM1");
  ParaPerpAsymmPhip_490MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for 480-500MeV Photon Energy (ThetaCM 0-30)");
  ParaPerpAsymmPhip_490MeVCM1->Fit("CosFit", "LL");
  pCosAmp[0][4] = CosFit->GetParameter(0);
  pCosAmpErr[0][4] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_510MeVCM1 = Phip_510MeVCM1_Para->GetAsymmetry(Phip_510MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_510MeVCM1->SetName("ParaPerpAsymmPhip510MeVCM1");
  ParaPerpAsymmPhip_510MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for 500-520MeV Photon Energy (ThetaCM 0-30)");
  ParaPerpAsymmPhip_510MeVCM1->Fit("CosFit", "LL");
  pCosAmp[0][5] = CosFit->GetParameter(0);
  pCosAmpErr[0][5] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_530MeVCM1 = Phip_530MeVCM1_Para->GetAsymmetry(Phip_530MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_530MeVCM1->SetName("ParaPerpAsymmPhip530MeVCM1");
  ParaPerpAsymmPhip_530MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for 520-540MeV Photon Energy (ThetaCM 0-30)");
  ParaPerpAsymmPhip_530MeVCM1->Fit("CosFit", "LL");
  pCosAmp[0][6] = CosFit->GetParameter(0);
  pCosAmpErr[0][6] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_550MeVCM1 = Phip_550MeVCM1_Para->GetAsymmetry(Phip_550MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_550MeVCM1->SetName("ParaPerpAsymmPhip550MeVCM1");
  ParaPerpAsymmPhip_550MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for 540-560MeV Photon Energy (ThetaCM 0-30)");
  ParaPerpAsymmPhip_550MeVCM1->Fit("CosFit", "LL");
  pCosAmp[0][7] = CosFit->GetParameter(0);
  pCosAmpErr[0][7] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_570MeVCM1 = Phip_570MeVCM1_Para->GetAsymmetry(Phip_570MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_570MeVCM1->SetName("ParaPerpAsymmPhip570MeVCM1");
  ParaPerpAsymmPhip_570MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for 560-580MeV Photon Energy (ThetaCM 0-30)");
  ParaPerpAsymmPhip_570MeVCM1->Fit("CosFit", "LL");
  pCosAmp[0][8] = CosFit->GetParameter(0);
  pCosAmpErr[0][8] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_590MeVCM1 = Phip_590MeVCM1_Para->GetAsymmetry(Phip_590MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_590MeVCM1->SetName("ParaPerpAsymmPhip590MeVCM1");
  ParaPerpAsymmPhip_590MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for 580-600MeV Photon Energy (ThetaCM 0-30)");
  ParaPerpAsymmPhip_590MeVCM1->Fit("CosFit", "LL");
  pCosAmp[0][9] = CosFit->GetParameter(0);
  pCosAmpErr[0][9] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_610MeVCM1 = Phip_610MeVCM1_Para->GetAsymmetry(Phip_610MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_610MeVCM1->SetName("ParaPerpAsymmPhip610MeVCM1");
  ParaPerpAsymmPhip_610MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for 600-620MeV Photon Energy (ThetaCM 0-30)");
  ParaPerpAsymmPhip_610MeVCM1->Fit("CosFit", "LL");
  pCosAmp[0][10] = CosFit->GetParameter(0);
  pCosAmpErr[0][10] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_630MeVCM1 = Phip_630MeVCM1_Para->GetAsymmetry(Phip_630MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_630MeVCM1->SetName("ParaPerpAsymmPhip630MeVCM1");
  ParaPerpAsymmPhip_630MeVCM1->SetTitle("Proton Para/Perp Phi Asymmetry for 620-640MeV Photon Energy (ThetaCM 0-30)");
  ParaPerpAsymmPhip_630MeVCM1->Fit("CosFit", "LL");
  pCosAmp[0][11] = CosFit->GetParameter(0);
  pCosAmpErr[0][11] = CosFit->GetParError(0);

  ///////////////////////////////////////////
  //////////////////CM2 p ///////////////////
  ///////////////////////////////////////////

  ParaPerpAsymmPhip_410MeVCM2 = Phip_410MeVCM2_Para->GetAsymmetry(Phip_410MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_410MeVCM2->SetName("ParaPerpAsymmPhip410MeVCM2");
  ParaPerpAsymmPhip_410MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for 400-420MeV Photon Energy (ThetaCM 30-60)");
  ParaPerpAsymmPhip_410MeVCM2->Fit("CosFit", "LL");
  pCosAmp[1][0] = CosFit->GetParameter(0);
  pCosAmpErr[1][0] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_430MeVCM2 = Phip_430MeVCM2_Para->GetAsymmetry(Phip_430MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_430MeVCM2->SetName("ParaPerpAsymmPhip430MeVCM2");
  ParaPerpAsymmPhip_430MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for 420-440MeV Photon Energy (ThetaCM 30-60)");
  ParaPerpAsymmPhip_430MeVCM2->Fit("CosFit", "LL");
  pCosAmp[1][1] = CosFit->GetParameter(0);
  pCosAmpErr[1][1] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_450MeVCM2 = Phip_450MeVCM2_Para->GetAsymmetry(Phip_450MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_450MeVCM2->SetName("ParaPerpAsymmPhip450MeVCM2");
  ParaPerpAsymmPhip_450MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for 440-460MeV Photon Energy (ThetaCM 30-60)");
  ParaPerpAsymmPhip_450MeVCM2->Fit("CosFit", "LL");
  pCosAmp[1][2] = CosFit->GetParameter(0);
  pCosAmpErr[1][2] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_470MeVCM2 = Phip_470MeVCM2_Para->GetAsymmetry(Phip_470MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_470MeVCM2->SetName("ParaPerpAsymmPhip470MeVCM2");
  ParaPerpAsymmPhip_470MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for 460-480MeV Photon Energy (ThetaCM 30-60)");
  ParaPerpAsymmPhip_470MeVCM2->Fit("CosFit", "LL");
  pCosAmp[1][3] = CosFit->GetParameter(0);
  pCosAmpErr[1][3] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_490MeVCM2 = Phip_490MeVCM2_Para->GetAsymmetry(Phip_490MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_490MeVCM2->SetName("ParaPerpAsymmPhip490MeVCM2");
  ParaPerpAsymmPhip_490MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for 480-500MeV Photon Energy (ThetaCM 30-60)");
  ParaPerpAsymmPhip_490MeVCM2->Fit("CosFit", "LL");
  pCosAmp[1][4] = CosFit->GetParameter(0);
  pCosAmpErr[1][4] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_510MeVCM2 = Phip_510MeVCM2_Para->GetAsymmetry(Phip_510MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_510MeVCM2->SetName("ParaPerpAsymmPhip510MeVCM2");
  ParaPerpAsymmPhip_510MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for 500-520MeV Photon Energy (ThetaCM 30-60)");
  ParaPerpAsymmPhip_510MeVCM2->Fit("CosFit", "LL");
  pCosAmp[1][5] = CosFit->GetParameter(0);
  pCosAmpErr[1][5] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_530MeVCM2 = Phip_530MeVCM2_Para->GetAsymmetry(Phip_530MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_530MeVCM2->SetName("ParaPerpAsymmPhip530MeVCM2");
  ParaPerpAsymmPhip_530MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for 520-540MeV Photon Energy (ThetaCM 30-60)");
  ParaPerpAsymmPhip_530MeVCM2->Fit("CosFit", "LL");
  pCosAmp[1][6] = CosFit->GetParameter(0);
  pCosAmpErr[1][6] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_550MeVCM2 = Phip_550MeVCM2_Para->GetAsymmetry(Phip_550MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_550MeVCM2->SetName("ParaPerpAsymmPhip550MeVCM2");
  ParaPerpAsymmPhip_550MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for 540-560MeV Photon Energy (ThetaCM 30-60)");
  ParaPerpAsymmPhip_550MeVCM2->Fit("CosFit", "LL");
  pCosAmp[1][7] = CosFit->GetParameter(0);
  pCosAmpErr[1][7] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_570MeVCM2 = Phip_570MeVCM2_Para->GetAsymmetry(Phip_570MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_570MeVCM2->SetName("ParaPerpAsymmPhip570MeVCM2");
  ParaPerpAsymmPhip_570MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for 560-580MeV Photon Energy (ThetaCM 30-60)");
  ParaPerpAsymmPhip_570MeVCM2->Fit("CosFit", "LL");
  pCosAmp[1][8] = CosFit->GetParameter(0);
  pCosAmpErr[1][8] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_590MeVCM2 = Phip_590MeVCM2_Para->GetAsymmetry(Phip_590MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_590MeVCM2->SetName("ParaPerpAsymmPhip590MeVCM2");
  ParaPerpAsymmPhip_590MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for 580-600MeV Photon Energy (ThetaCM 30-60)");
  ParaPerpAsymmPhip_590MeVCM2->Fit("CosFit", "LL");
  pCosAmp[1][9] = CosFit->GetParameter(0);
  pCosAmpErr[1][9] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_610MeVCM2 = Phip_610MeVCM2_Para->GetAsymmetry(Phip_610MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_610MeVCM2->SetName("ParaPerpAsymmPhip610MeVCM2");
  ParaPerpAsymmPhip_610MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for 600-620MeV Photon Energy (ThetaCM 30-60)");
  ParaPerpAsymmPhip_610MeVCM2->Fit("CosFit", "LL");
  pCosAmp[1][10] = CosFit->GetParameter(0);
  pCosAmpErr[1][10] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_630MeVCM2 = Phip_630MeVCM2_Para->GetAsymmetry(Phip_630MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_630MeVCM2->SetName("ParaPerpAsymmPhip630MeVCM2");
  ParaPerpAsymmPhip_630MeVCM2->SetTitle("Proton Para/Perp Phi Asymmetry for 620-640MeV Photon Energy (ThetaCM 30-60)");
  ParaPerpAsymmPhip_630MeVCM2->Fit("CosFit", "LL");
  pCosAmp[1][11] = CosFit->GetParameter(0);
  pCosAmpErr[1][11] = CosFit->GetParError(0);

  ///////////////////////////////////////////
  //////////////////CM3 p ///////////////////
  ///////////////////////////////////////////

  ParaPerpAsymmPhip_410MeVCM3 = Phip_410MeVCM3_Para->GetAsymmetry(Phip_410MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_410MeVCM3->SetName("ParaPerpAsymmPhip410MeVCM3");
  ParaPerpAsymmPhip_410MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for 400-420MeV Photon Energy (ThetaCM 60-90)");
  ParaPerpAsymmPhip_410MeVCM3->Fit("CosFit", "LL");
  pCosAmp[2][0] = CosFit->GetParameter(0);
  pCosAmpErr[2][0] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_430MeVCM3 = Phip_430MeVCM3_Para->GetAsymmetry(Phip_430MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_430MeVCM3->SetName("ParaPerpAsymmPhip430MeVCM3");
  ParaPerpAsymmPhip_430MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for 420-440MeV Photon Energy (ThetaCM 60-90)");
  ParaPerpAsymmPhip_430MeVCM3->Fit("CosFit", "LL");
  pCosAmp[2][1] = CosFit->GetParameter(0);
  pCosAmpErr[2][1] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_450MeVCM3 = Phip_450MeVCM3_Para->GetAsymmetry(Phip_450MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_450MeVCM3->SetName("ParaPerpAsymmPhip450MeVCM3");
  ParaPerpAsymmPhip_450MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for 440-460MeV Photon Energy (ThetaCM 60-90)");
  ParaPerpAsymmPhip_450MeVCM3->Fit("CosFit", "LL");
  pCosAmp[2][2] = CosFit->GetParameter(0);
  pCosAmpErr[2][2] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_470MeVCM3 = Phip_470MeVCM3_Para->GetAsymmetry(Phip_470MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_470MeVCM3->SetName("ParaPerpAsymmPhip470MeVCM3");
  ParaPerpAsymmPhip_470MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for 460-480MeV Photon Energy (ThetaCM 60-90)");
  ParaPerpAsymmPhip_470MeVCM3->Fit("CosFit", "LL");
  pCosAmp[2][3] = CosFit->GetParameter(0);
  pCosAmpErr[2][3] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_490MeVCM3 = Phip_490MeVCM3_Para->GetAsymmetry(Phip_490MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_490MeVCM3->SetName("ParaPerpAsymmPhip490MeVCM3");
  ParaPerpAsymmPhip_490MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for 480-500MeV Photon Energy (ThetaCM 60-90)");
  ParaPerpAsymmPhip_490MeVCM3->Fit("CosFit", "LL");
  pCosAmp[2][4] = CosFit->GetParameter(0);
  pCosAmpErr[2][4] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_510MeVCM3 = Phip_510MeVCM3_Para->GetAsymmetry(Phip_510MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_510MeVCM3->SetName("ParaPerpAsymmPhip510MeVCM3");
  ParaPerpAsymmPhip_510MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for 500-520MeV Photon Energy (ThetaCM 60-90)");
  ParaPerpAsymmPhip_510MeVCM3->Fit("CosFit", "LL");
  pCosAmp[2][5] = CosFit->GetParameter(0);
  pCosAmpErr[2][5] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_530MeVCM3 = Phip_530MeVCM3_Para->GetAsymmetry(Phip_530MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_530MeVCM3->SetName("ParaPerpAsymmPhip530MeVCM3");
  ParaPerpAsymmPhip_530MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for 520-540MeV Photon Energy (ThetaCM 60-90)");
  ParaPerpAsymmPhip_530MeVCM3->Fit("CosFit", "LL");
  pCosAmp[2][6] = CosFit->GetParameter(0);
  pCosAmpErr[2][6] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_550MeVCM3 = Phip_550MeVCM3_Para->GetAsymmetry(Phip_550MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_550MeVCM3->SetName("ParaPerpAsymmPhip550MeVCM3");
  ParaPerpAsymmPhip_550MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for 540-560MeV Photon Energy (ThetaCM 60-90)");
  ParaPerpAsymmPhip_550MeVCM3->Fit("CosFit", "LL");
  pCosAmp[2][7] = CosFit->GetParameter(0);
  pCosAmpErr[2][7] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_570MeVCM3 = Phip_570MeVCM3_Para->GetAsymmetry(Phip_570MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_570MeVCM3->SetName("ParaPerpAsymmPhip570MeVCM3");
  ParaPerpAsymmPhip_570MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for 560-580MeV Photon Energy (ThetaCM 60-90)");
  ParaPerpAsymmPhip_570MeVCM3->Fit("CosFit", "LL");
  pCosAmp[2][8] = CosFit->GetParameter(0);
  pCosAmpErr[2][8] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_590MeVCM3 = Phip_590MeVCM3_Para->GetAsymmetry(Phip_590MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_590MeVCM3->SetName("ParaPerpAsymmPhip590MeVCM3");
  ParaPerpAsymmPhip_590MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for 580-600MeV Photon Energy (ThetaCM 60-90)");
  ParaPerpAsymmPhip_590MeVCM3->Fit("CosFit", "LL");
  pCosAmp[2][9] = CosFit->GetParameter(0);
  pCosAmpErr[2][9] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_610MeVCM3 = Phip_610MeVCM3_Para->GetAsymmetry(Phip_610MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_610MeVCM3->SetName("ParaPerpAsymmPhip610MeVCM3");
  ParaPerpAsymmPhip_610MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for 600-620MeV Photon Energy (ThetaCM 60-90)");
  ParaPerpAsymmPhip_610MeVCM3->Fit("CosFit", "LL");
  pCosAmp[2][10] = CosFit->GetParameter(0);
  pCosAmpErr[2][10] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_630MeVCM3 = Phip_630MeVCM3_Para->GetAsymmetry(Phip_630MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_630MeVCM3->SetName("ParaPerpAsymmPhip630MeVCM3");
  ParaPerpAsymmPhip_630MeVCM3->SetTitle("Proton Para/Perp Phi Asymmetry for 620-640MeV Photon Energy (ThetaCM 60-90)");
  ParaPerpAsymmPhip_630MeVCM3->Fit("CosFit", "LL");
  pCosAmp[2][11] = CosFit->GetParameter(0);
  pCosAmpErr[2][11] = CosFit->GetParError(0);

  ///////////////////////////////////////////
  //////////////////CM4 p ///////////////////
  ///////////////////////////////////////////

  ParaPerpAsymmPhip_410MeVCM4 = Phip_410MeVCM4_Para->GetAsymmetry(Phip_410MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_410MeVCM4->SetName("ParaPerpAsymmPhip410MeVCM4");
  ParaPerpAsymmPhip_410MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for 400-420MeV Photon Energy (ThetaCM 90-120)");
  ParaPerpAsymmPhip_410MeVCM4->Fit("CosFit", "LL");
  pCosAmp[3][0] = CosFit->GetParameter(0);
  pCosAmpErr[3][0] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_430MeVCM4 = Phip_430MeVCM4_Para->GetAsymmetry(Phip_430MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_430MeVCM4->SetName("ParaPerpAsymmPhip430MeVCM4");
  ParaPerpAsymmPhip_430MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for 420-440MeV Photon Energy (ThetaCM 90-120)");
  ParaPerpAsymmPhip_430MeVCM4->Fit("CosFit", "LL");
  pCosAmp[3][1] = CosFit->GetParameter(0);
  pCosAmpErr[3][1] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_450MeVCM4 = Phip_450MeVCM4_Para->GetAsymmetry(Phip_450MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_450MeVCM4->SetName("ParaPerpAsymmPhip450MeVCM4");
  ParaPerpAsymmPhip_450MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for 440-460MeV Photon Energy (ThetaCM 90-120)");
  ParaPerpAsymmPhip_450MeVCM4->Fit("CosFit", "LL");
  pCosAmp[3][2] = CosFit->GetParameter(0);
  pCosAmpErr[3][2] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_470MeVCM4 = Phip_470MeVCM4_Para->GetAsymmetry(Phip_470MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_470MeVCM4->SetName("ParaPerpAsymmPhip470MeVCM4");
  ParaPerpAsymmPhip_470MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for 460-480MeV Photon Energy (ThetaCM 90-120)");
  ParaPerpAsymmPhip_470MeVCM4->Fit("CosFit", "LL");
  pCosAmp[3][3] = CosFit->GetParameter(0);
  pCosAmpErr[3][3] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_490MeVCM4 = Phip_490MeVCM4_Para->GetAsymmetry(Phip_490MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_490MeVCM4->SetName("ParaPerpAsymmPhip490MeVCM4");
  ParaPerpAsymmPhip_490MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for 480-500MeV Photon Energy (ThetaCM 90-120)");
  ParaPerpAsymmPhip_490MeVCM4->Fit("CosFit", "LL");
  pCosAmp[3][4] = CosFit->GetParameter(0);
  pCosAmpErr[3][4] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_510MeVCM4 = Phip_510MeVCM4_Para->GetAsymmetry(Phip_510MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_510MeVCM4->SetName("ParaPerpAsymmPhip510MeVCM4");
  ParaPerpAsymmPhip_510MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for 500-520MeV Photon Energy (ThetaCM 90-120)");
  ParaPerpAsymmPhip_510MeVCM4->Fit("CosFit", "LL");
  pCosAmp[3][5] = CosFit->GetParameter(0);
  pCosAmpErr[3][5] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_530MeVCM4 = Phip_530MeVCM4_Para->GetAsymmetry(Phip_530MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_530MeVCM4->SetName("ParaPerpAsymmPhip530MeVCM4");
  ParaPerpAsymmPhip_530MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for 520-540MeV Photon Energy (ThetaCM 90-120)");
  ParaPerpAsymmPhip_530MeVCM4->Fit("CosFit", "LL");
  pCosAmp[3][6] = CosFit->GetParameter(0);
  pCosAmpErr[3][6] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_550MeVCM4 = Phip_550MeVCM4_Para->GetAsymmetry(Phip_550MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_550MeVCM4->SetName("ParaPerpAsymmPhip550MeVCM4");
  ParaPerpAsymmPhip_550MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for 540-560MeV Photon Energy (ThetaCM 90-120)");
  ParaPerpAsymmPhip_550MeVCM4->Fit("CosFit", "LL");
  pCosAmp[3][7] = CosFit->GetParameter(0);
  pCosAmpErr[3][7] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_570MeVCM4 = Phip_570MeVCM4_Para->GetAsymmetry(Phip_570MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_570MeVCM4->SetName("ParaPerpAsymmPhip570MeVCM4");
  ParaPerpAsymmPhip_570MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for 560-580MeV Photon Energy (ThetaCM 90-120)");
  ParaPerpAsymmPhip_570MeVCM4->Fit("CosFit", "LL");
  pCosAmp[3][8] = CosFit->GetParameter(0);
  pCosAmpErr[3][8] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_590MeVCM4 = Phip_590MeVCM4_Para->GetAsymmetry(Phip_590MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_590MeVCM4->SetName("ParaPerpAsymmPhip590MeVCM4");
  ParaPerpAsymmPhip_590MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for 580-600MeV Photon Energy (ThetaCM 90-120)");
  ParaPerpAsymmPhip_590MeVCM4->Fit("CosFit", "LL");
  pCosAmp[3][9] = CosFit->GetParameter(0);
  pCosAmpErr[3][9] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_610MeVCM4 = Phip_610MeVCM4_Para->GetAsymmetry(Phip_610MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_610MeVCM4->SetName("ParaPerpAsymmPhip610MeVCM4");
  ParaPerpAsymmPhip_610MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for 600-620MeV Photon Energy (ThetaCM 90-120)");
  ParaPerpAsymmPhip_610MeVCM4->Fit("CosFit", "LL");
  pCosAmp[3][10] = CosFit->GetParameter(0);
  pCosAmpErr[3][10] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_630MeVCM4 = Phip_630MeVCM4_Para->GetAsymmetry(Phip_630MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_630MeVCM4->SetName("ParaPerpAsymmPhip630MeVCM4");
  ParaPerpAsymmPhip_630MeVCM4->SetTitle("Proton Para/Perp Phi Asymmetry for 620-640MeV Photon Energy (ThetaCM 90-120)");
  ParaPerpAsymmPhip_630MeVCM4->Fit("CosFit", "LL");
  pCosAmp[3][11] = CosFit->GetParameter(0);
  pCosAmpErr[3][11] = CosFit->GetParError(0);

  ///////////////////////////////////////////
  //////////////////CM5 p ///////////////////
  ///////////////////////////////////////////

  ParaPerpAsymmPhip_410MeVCM5 = Phip_410MeVCM5_Para->GetAsymmetry(Phip_410MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_410MeVCM5->SetName("ParaPerpAsymmPhip410MeVCM5");
  ParaPerpAsymmPhip_410MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for 400-420MeV Photon Energy (ThetaCM 120-150)");
  ParaPerpAsymmPhip_410MeVCM5->Fit("CosFit", "LL");
  pCosAmp[4][0] = CosFit->GetParameter(0);
  pCosAmpErr[4][0] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_430MeVCM5 = Phip_430MeVCM5_Para->GetAsymmetry(Phip_430MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_430MeVCM5->SetName("ParaPerpAsymmPhip430MeVCM5");
  ParaPerpAsymmPhip_430MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for 420-440MeV Photon Energy (ThetaCM 120-150)");
  ParaPerpAsymmPhip_430MeVCM5->Fit("CosFit", "LL");
  pCosAmp[4][1] = CosFit->GetParameter(0);
  pCosAmpErr[4][1] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_450MeVCM5 = Phip_450MeVCM5_Para->GetAsymmetry(Phip_450MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_450MeVCM5->SetName("ParaPerpAsymmPhip450MeVCM5");
  ParaPerpAsymmPhip_450MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for 440-460MeV Photon Energy (ThetaCM 120-150)");
  ParaPerpAsymmPhip_450MeVCM5->Fit("CosFit", "LL");
  pCosAmp[4][2] = CosFit->GetParameter(0);
  pCosAmpErr[4][2] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_470MeVCM5 = Phip_470MeVCM5_Para->GetAsymmetry(Phip_470MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_470MeVCM5->SetName("ParaPerpAsymmPhip470MeVCM5");
  ParaPerpAsymmPhip_470MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for 460-480MeV Photon Energy (ThetaCM 120-150)");
  ParaPerpAsymmPhip_470MeVCM5->Fit("CosFit", "LL");
  pCosAmp[4][3] = CosFit->GetParameter(0);
  pCosAmpErr[4][3] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_490MeVCM5 = Phip_490MeVCM5_Para->GetAsymmetry(Phip_490MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_490MeVCM5->SetName("ParaPerpAsymmPhip490MeVCM5");
  ParaPerpAsymmPhip_490MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for 480-500MeV Photon Energy (ThetaCM 120-150)");
  ParaPerpAsymmPhip_490MeVCM5->Fit("CosFit", "LL");
  pCosAmp[4][4] = CosFit->GetParameter(0);
  pCosAmpErr[4][4] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_510MeVCM5 = Phip_510MeVCM5_Para->GetAsymmetry(Phip_510MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_510MeVCM5->SetName("ParaPerpAsymmPhip510MeVCM5");
  ParaPerpAsymmPhip_510MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for 500-520MeV Photon Energy (ThetaCM 120-150)");
  ParaPerpAsymmPhip_510MeVCM5->Fit("CosFit", "LL");
  pCosAmp[4][5] = CosFit->GetParameter(0);
  pCosAmpErr[4][5] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_530MeVCM5 = Phip_530MeVCM5_Para->GetAsymmetry(Phip_530MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_530MeVCM5->SetName("ParaPerpAsymmPhip530MeVCM5");
  ParaPerpAsymmPhip_530MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for 520-540MeV Photon Energy (ThetaCM 120-150)");
  ParaPerpAsymmPhip_530MeVCM5->Fit("CosFit", "LL");
  pCosAmp[4][6] = CosFit->GetParameter(0);
  pCosAmpErr[4][6] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_550MeVCM5 = Phip_550MeVCM5_Para->GetAsymmetry(Phip_550MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_550MeVCM5->SetName("ParaPerpAsymmPhip550MeVCM5");
  ParaPerpAsymmPhip_550MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for 540-560MeV Photon Energy (ThetaCM 120-150)");
  ParaPerpAsymmPhip_550MeVCM5->Fit("CosFit", "LL");
  pCosAmp[4][7] = CosFit->GetParameter(0);
  pCosAmpErr[4][7] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_570MeVCM5 = Phip_570MeVCM5_Para->GetAsymmetry(Phip_570MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_570MeVCM5->SetName("ParaPerpAsymmPhip570MeVCM5");
  ParaPerpAsymmPhip_570MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for 560-580MeV Photon Energy (ThetaCM 120-150)");
  ParaPerpAsymmPhip_570MeVCM5->Fit("CosFit", "LL");
  pCosAmp[4][8] = CosFit->GetParameter(0);
  pCosAmpErr[4][8] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_590MeVCM5 = Phip_590MeVCM5_Para->GetAsymmetry(Phip_590MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_590MeVCM5->SetName("ParaPerpAsymmPhip590MeVCM5");
  ParaPerpAsymmPhip_590MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for 580-600MeV Photon Energy (ThetaCM 120-150)");
  ParaPerpAsymmPhip_590MeVCM5->Fit("CosFit", "LL");
  pCosAmp[4][9] = CosFit->GetParameter(0);
  pCosAmpErr[4][9] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_610MeVCM5 = Phip_610MeVCM5_Para->GetAsymmetry(Phip_610MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_610MeVCM5->SetName("ParaPerpAsymmPhip610MeVCM5");
  ParaPerpAsymmPhip_610MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for 600-620MeV Photon Energy (ThetaCM 120-150)");
  ParaPerpAsymmPhip_610MeVCM5->Fit("CosFit", "LL");
  pCosAmp[4][10] = CosFit->GetParameter(0);
  pCosAmpErr[4][10] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_630MeVCM5 = Phip_630MeVCM5_Para->GetAsymmetry(Phip_630MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_630MeVCM5->SetName("ParaPerpAsymmPhip630MeVCM5");
  ParaPerpAsymmPhip_630MeVCM5->SetTitle("Proton Para/Perp Phi Asymmetry for 620-640MeV Photon Energy (ThetaCM 120-150)");
  ParaPerpAsymmPhip_630MeVCM5->Fit("CosFit", "LL");
  pCosAmp[4][11] = CosFit->GetParameter(0);
  pCosAmpErr[4][11] = CosFit->GetParError(0);

  ///////////////////////////////////////////
  //////////////////CM6 p ///////////////////
  ///////////////////////////////////////////

  ParaPerpAsymmPhip_410MeVCM6 = Phip_410MeVCM6_Para->GetAsymmetry(Phip_410MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_410MeVCM6->SetName("ParaPerpAsymmPhip410MeVCM6");
  ParaPerpAsymmPhip_410MeVCM6->SetTitle("Proton Para/Perp Phi Asymmetry for 400-420MeV Photon Energy (ThetaCM 150-180)");
  ParaPerpAsymmPhip_410MeVCM6->Fit("CosFit", "LL");
  pCosAmp[5][0] = CosFit->GetParameter(0);
  pCosAmpErr[5][0] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_430MeVCM6 = Phip_430MeVCM6_Para->GetAsymmetry(Phip_430MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_430MeVCM6->SetName("ParaPerpAsymmPhip430MeVCM6");
  ParaPerpAsymmPhip_430MeVCM6->SetTitle("Proton Para/Perp Phi Asymmetry for 420-440MeV Photon Energy (ThetaCM 150-180)");
  ParaPerpAsymmPhip_430MeVCM6->Fit("CosFit", "LL");
  pCosAmp[5][1] = CosFit->GetParameter(0);
  pCosAmpErr[5][1] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_450MeVCM6 = Phip_450MeVCM6_Para->GetAsymmetry(Phip_450MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_450MeVCM6->SetName("ParaPerpAsymmPhip450MeVCM6");
  ParaPerpAsymmPhip_450MeVCM6->SetTitle("Proton Para/Perp Phi Asymmetry for 440-460MeV Photon Energy (ThetaCM 150-180)");
  ParaPerpAsymmPhip_450MeVCM6->Fit("CosFit", "LL");
  pCosAmp[5][2] = CosFit->GetParameter(0);
  pCosAmpErr[5][2] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_470MeVCM6 = Phip_470MeVCM6_Para->GetAsymmetry(Phip_470MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_470MeVCM6->SetName("ParaPerpAsymmPhip470MeVCM6");
  ParaPerpAsymmPhip_470MeVCM6->SetTitle("Proton Para/Perp Phi Asymmetry for 460-480MeV Photon Energy (ThetaCM 150-180)");
  ParaPerpAsymmPhip_470MeVCM6->Fit("CosFit", "LL");
  pCosAmp[5][3] = CosFit->GetParameter(0);
  pCosAmpErr[5][3] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_490MeVCM6 = Phip_490MeVCM6_Para->GetAsymmetry(Phip_490MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_490MeVCM6->SetName("ParaPerpAsymmPhip490MeVCM6");
  ParaPerpAsymmPhip_490MeVCM6->SetTitle("Proton Para/Perp Phi Asymmetry for 480-500MeV Photon Energy (ThetaCM 150-180)");
  ParaPerpAsymmPhip_490MeVCM6->Fit("CosFit", "LL");
  pCosAmp[5][4] = CosFit->GetParameter(0);
  pCosAmpErr[5][4] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_510MeVCM6 = Phip_510MeVCM6_Para->GetAsymmetry(Phip_510MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_510MeVCM6->SetName("ParaPerpAsymmPhip510MeVCM6");
  ParaPerpAsymmPhip_510MeVCM6->SetTitle("Proton Para/Perp Phi Asymmetry for 500-520MeV Photon Energy (ThetaCM 150-180)");
  ParaPerpAsymmPhip_510MeVCM6->Fit("CosFit", "LL");
  pCosAmp[5][5] = CosFit->GetParameter(0);
  pCosAmpErr[5][5] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_530MeVCM6 = Phip_530MeVCM6_Para->GetAsymmetry(Phip_530MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_530MeVCM6->SetName("ParaPerpAsymmPhip530MeVCM6");
  ParaPerpAsymmPhip_530MeVCM6->SetTitle("Proton Para/Perp Phi Asymmetry for 520-540MeV Photon Energy (ThetaCM 150-180)");
  ParaPerpAsymmPhip_530MeVCM6->Fit("CosFit", "LL");
  pCosAmp[5][6] = CosFit->GetParameter(0);
  pCosAmpErr[5][6] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_550MeVCM6 = Phip_550MeVCM6_Para->GetAsymmetry(Phip_550MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_550MeVCM6->SetName("ParaPerpAsymmPhip550MeVCM6");
  ParaPerpAsymmPhip_550MeVCM6->SetTitle("Proton Para/Perp Phi Asymmetry for 540-560MeV Photon Energy (ThetaCM 150-180)");
  ParaPerpAsymmPhip_550MeVCM6->Fit("CosFit", "LL");
  pCosAmp[5][7] = CosFit->GetParameter(0);
  pCosAmpErr[5][7] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_570MeVCM6 = Phip_570MeVCM6_Para->GetAsymmetry(Phip_570MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_570MeVCM6->SetName("ParaPerpAsymmPhip570MeVCM6");
  ParaPerpAsymmPhip_570MeVCM6->SetTitle("Proton Para/Perp Phi Asymmetry for 560-580MeV Photon Energy (ThetaCM 150-180)");
  ParaPerpAsymmPhip_570MeVCM6->Fit("CosFit", "LL");
  pCosAmp[5][8] = CosFit->GetParameter(0);
  pCosAmpErr[5][8] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_590MeVCM6 = Phip_590MeVCM6_Para->GetAsymmetry(Phip_590MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_590MeVCM6->SetName("ParaPerpAsymmPhip590MeVCM6");
  ParaPerpAsymmPhip_590MeVCM6->SetTitle("Proton Para/Perp Phi Asymmetry for 580-600MeV Photon Energy (ThetaCM 150-180)");
  ParaPerpAsymmPhip_590MeVCM6->Fit("CosFit", "LL");
  pCosAmp[5][9] = CosFit->GetParameter(0);
  pCosAmpErr[5][9] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_610MeVCM6 = Phip_610MeVCM6_Para->GetAsymmetry(Phip_610MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_610MeVCM6->SetName("ParaPerpAsymmPhip610MeVCM6");
  ParaPerpAsymmPhip_610MeVCM6->SetTitle("Proton Para/Perp Phi Asymmetry for 600-620MeV Photon Energy (ThetaCM 150-180)");
  ParaPerpAsymmPhip_610MeVCM6->Fit("CosFit", "LL");
  pCosAmp[5][10] = CosFit->GetParameter(0);
  pCosAmpErr[5][10] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_630MeVCM6 = Phip_630MeVCM6_Para->GetAsymmetry(Phip_630MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_630MeVCM6->SetName("ParaPerpAsymmPhip630MeVCM6");
  ParaPerpAsymmPhip_630MeVCM6->SetTitle("Proton Para/Perp Phi Asymmetry for 620-640MeV Photon Energy (ThetaCM 150-180)");
  ParaPerpAsymmPhip_630MeVCM6->Fit("CosFit", "LL");
  pCosAmp[5][11] = CosFit->GetParameter(0);
  pCosAmpErr[5][11] = CosFit->GetParError(0);

  /////////////////////////////////////
  ////////////// Phin /////////////////
  /////////////////////////////////////

  ParaPerpAsymmPhin_410MeVCM1 = Phin_410MeVCM1_Para->GetAsymmetry(Phin_410MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_410MeVCM1->SetName("ParaPerpAsymmPhin410MeVCM1");
  ParaPerpAsymmPhin_410MeVCM1->SetTitle("Neutron Para/Perp Phi Asymmetry for 400-420MeV Photon Energy (ThetaCM 0-30)");
  ParaPerpAsymmPhin_410MeVCM1->Fit("CosFit", "LL");
  nCosAmp[0][0] = CosFit->GetParameter(0);
  nCosAmpErr[0][0] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_430MeVCM1 = Phin_430MeVCM1_Para->GetAsymmetry(Phin_430MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_430MeVCM1->SetName("ParaPerpAsymmPhin430MeVCM1");
  ParaPerpAsymmPhin_430MeVCM1->SetTitle("Neutron Para/Perp Phi Asymmetry for 420-440MeV Photon Energy (ThetaCM 0-30)");
  ParaPerpAsymmPhin_430MeVCM1->Fit("CosFit", "LL");
  nCosAmp[0][1] = CosFit->GetParameter(0);
  nCosAmpErr[0][1] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_450MeVCM1 = Phin_450MeVCM1_Para->GetAsymmetry(Phin_450MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_450MeVCM1->SetName("ParaPerpAsymmPhin450MeVCM1");
  ParaPerpAsymmPhin_450MeVCM1->SetTitle("Neutron Para/Perp Phi Asymmetry for 440-460MeV Photon Energy (ThetaCM 0-30)");
  ParaPerpAsymmPhin_450MeVCM1->Fit("CosFit", "LL");
  nCosAmp[0][2] = CosFit->GetParameter(0);
  nCosAmpErr[0][2] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_470MeVCM1 = Phin_470MeVCM1_Para->GetAsymmetry(Phin_470MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_470MeVCM1->SetName("ParaPerpAsymmPhin470MeVCM1");
  ParaPerpAsymmPhin_470MeVCM1->SetTitle("Neutron Para/Perp Phi Asymmetry for 460-480MeV Photon Energy (ThetaCM 0-30)");
  ParaPerpAsymmPhin_470MeVCM1->Fit("CosFit", "LL");
  nCosAmp[0][3] = CosFit->GetParameter(0);
  nCosAmpErr[0][3] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_490MeVCM1 = Phin_490MeVCM1_Para->GetAsymmetry(Phin_490MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_490MeVCM1->SetName("ParaPerpAsymmPhin490MeVCM1");
  ParaPerpAsymmPhin_490MeVCM1->SetTitle("Neutron Para/Perp Phi Asymmetry for 480-500MeV Photon Energy (ThetaCM 0-30)");
  ParaPerpAsymmPhin_490MeVCM1->Fit("CosFit", "LL");
  nCosAmp[0][4] = CosFit->GetParameter(0);
  nCosAmpErr[0][4] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_510MeVCM1 = Phin_510MeVCM1_Para->GetAsymmetry(Phin_510MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_510MeVCM1->SetName("ParaPerpAsymmPhin510MeVCM1");
  ParaPerpAsymmPhin_510MeVCM1->SetTitle("Neutron Para/Perp Phi Asymmetry for 500-520MeV Photon Energy (ThetaCM 0-30)");
  ParaPerpAsymmPhin_510MeVCM1->Fit("CosFit", "LL");
  nCosAmp[0][5] = CosFit->GetParameter(0);
  nCosAmpErr[0][5] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_530MeVCM1 = Phin_530MeVCM1_Para->GetAsymmetry(Phin_530MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_530MeVCM1->SetName("ParaPerpAsymmPhin530MeVCM1");
  ParaPerpAsymmPhin_530MeVCM1->SetTitle("Neutron Para/Perp Phi Asymmetry for 520-540MeV Photon Energy (ThetaCM 0-30)");
  ParaPerpAsymmPhin_530MeVCM1->Fit("CosFit", "LL");
  nCosAmp[0][6] = CosFit->GetParameter(0);
  nCosAmpErr[0][6] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_550MeVCM1 = Phin_550MeVCM1_Para->GetAsymmetry(Phin_550MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_550MeVCM1->SetName("ParaPerpAsymmPhin550MeVCM1");
  ParaPerpAsymmPhin_550MeVCM1->SetTitle("Neutron Para/Perp Phi Asymmetry for 540-560MeV Photon Energy (ThetaCM 0-30)");
  ParaPerpAsymmPhin_550MeVCM1->Fit("CosFit", "LL");
  nCosAmp[0][7] = CosFit->GetParameter(0);
  nCosAmpErr[0][7] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_570MeVCM1 = Phin_570MeVCM1_Para->GetAsymmetry(Phin_570MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_570MeVCM1->SetName("ParaPerpAsymmPhin570MeVCM1");
  ParaPerpAsymmPhin_570MeVCM1->SetTitle("Neutron Para/Perp Phi Asymmetry for 560-580MeV Photon Energy (ThetaCM 0-30)");
  ParaPerpAsymmPhin_570MeVCM1->Fit("CosFit", "LL");
  nCosAmp[0][8] = CosFit->GetParameter(0);
  nCosAmpErr[0][8] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_590MeVCM1 = Phin_590MeVCM1_Para->GetAsymmetry(Phin_590MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_590MeVCM1->SetName("ParaPerpAsymmPhin590MeVCM1");
  ParaPerpAsymmPhin_590MeVCM1->SetTitle("Neutron Para/Perp Phi Asymmetry for 580-600MeV Photon Energy (ThetaCM 0-30)");
  ParaPerpAsymmPhin_590MeVCM1->Fit("CosFit", "LL");
  nCosAmp[0][9] = CosFit->GetParameter(0);
  nCosAmpErr[0][9] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_610MeVCM1 = Phin_610MeVCM1_Para->GetAsymmetry(Phin_610MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_610MeVCM1->SetName("ParaPerpAsymmPhin610MeVCM1");
  ParaPerpAsymmPhin_610MeVCM1->SetTitle("Neutron Para/Perp Phi Asymmetry for 600-620MeV Photon Energy (ThetaCM 0-30)");
  ParaPerpAsymmPhin_610MeVCM1->Fit("CosFit", "LL");
  nCosAmp[0][10] = CosFit->GetParameter(0);
  nCosAmpErr[0][10] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_630MeVCM1 = Phin_630MeVCM1_Para->GetAsymmetry(Phin_630MeVCM1_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_630MeVCM1->SetName("ParaPerpAsymmPhin630MeVCM1");
  ParaPerpAsymmPhin_630MeVCM1->SetTitle("Neutron Para/Perp Phi Asymmetry for 620-640MeV Photon Energy (ThetaCM 0-30)");
  ParaPerpAsymmPhin_630MeVCM1->Fit("CosFit", "LL");
  nCosAmp[0][11] = CosFit->GetParameter(0);
  nCosAmpErr[0][11] = CosFit->GetParError(0);

  ///////////////////////////////////////////
  //////////////////CM2 p ///////////////////
  ///////////////////////////////////////////

  ParaPerpAsymmPhin_410MeVCM2 = Phin_410MeVCM2_Para->GetAsymmetry(Phin_410MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_410MeVCM2->SetName("ParaPerpAsymmPhin410MeVCM2");
  ParaPerpAsymmPhin_410MeVCM2->SetTitle("Neutron Para/Perp Phi Asymmetry for 400-420MeV Photon Energy (ThetaCM 30-60)");
  ParaPerpAsymmPhin_410MeVCM2->Fit("CosFit", "LL");
  nCosAmp[1][0] = CosFit->GetParameter(0);
  nCosAmpErr[1][0] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_430MeVCM2 = Phin_430MeVCM2_Para->GetAsymmetry(Phin_430MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_430MeVCM2->SetName("ParaPerpAsymmPhin430MeVCM2");
  ParaPerpAsymmPhin_430MeVCM2->SetTitle("Neutron Para/Perp Phi Asymmetry for 420-440MeV Photon Energy (ThetaCM 30-60)");
  ParaPerpAsymmPhin_430MeVCM2->Fit("CosFit", "LL");
  nCosAmp[1][1] = CosFit->GetParameter(0);
  nCosAmpErr[1][1] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_450MeVCM2 = Phin_450MeVCM2_Para->GetAsymmetry(Phin_450MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_450MeVCM2->SetName("ParaPerpAsymmPhin450MeVCM2");
  ParaPerpAsymmPhin_450MeVCM2->SetTitle("Neutron Para/Perp Phi Asymmetry for 440-460MeV Photon Energy (ThetaCM 30-60)");
  ParaPerpAsymmPhin_450MeVCM2->Fit("CosFit", "LL");
  nCosAmp[1][2] = CosFit->GetParameter(0);
  nCosAmpErr[1][2] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_470MeVCM2 = Phin_470MeVCM2_Para->GetAsymmetry(Phin_470MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_470MeVCM2->SetName("ParaPerpAsymmPhin470MeVCM2");
  ParaPerpAsymmPhin_470MeVCM2->SetTitle("Neutron Para/Perp Phi Asymmetry for 460-480MeV Photon Energy (ThetaCM 30-60)");
  ParaPerpAsymmPhin_470MeVCM2->Fit("CosFit", "LL");
  nCosAmp[1][3] = CosFit->GetParameter(0);
  nCosAmpErr[1][3] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_490MeVCM2 = Phin_490MeVCM2_Para->GetAsymmetry(Phin_490MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_490MeVCM2->SetName("ParaPerpAsymmPhin490MeVCM2");
  ParaPerpAsymmPhin_490MeVCM2->SetTitle("Neutron Para/Perp Phi Asymmetry for 480-500MeV Photon Energy (ThetaCM 30-60)");
  ParaPerpAsymmPhin_490MeVCM2->Fit("CosFit", "LL");
  nCosAmp[1][4] = CosFit->GetParameter(0);
  nCosAmpErr[1][4] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_510MeVCM2 = Phin_510MeVCM2_Para->GetAsymmetry(Phin_510MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_510MeVCM2->SetName("ParaPerpAsymmPhin510MeVCM2");
  ParaPerpAsymmPhin_510MeVCM2->SetTitle("Neutron Para/Perp Phi Asymmetry for 500-520MeV Photon Energy (ThetaCM 30-60)");
  ParaPerpAsymmPhin_510MeVCM2->Fit("CosFit", "LL");
  nCosAmp[1][5] = CosFit->GetParameter(0);
  nCosAmpErr[1][5] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_530MeVCM2 = Phin_530MeVCM2_Para->GetAsymmetry(Phin_530MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_530MeVCM2->SetName("ParaPerpAsymmPhin530MeVCM2");
  ParaPerpAsymmPhin_530MeVCM2->SetTitle("Neutron Para/Perp Phi Asymmetry for 520-540MeV Photon Energy (ThetaCM 30-60)");
  ParaPerpAsymmPhin_530MeVCM2->Fit("CosFit", "LL");
  nCosAmp[1][6] = CosFit->GetParameter(0);
  nCosAmpErr[1][6] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_550MeVCM2 = Phin_550MeVCM2_Para->GetAsymmetry(Phin_550MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_550MeVCM2->SetName("ParaPerpAsymmPhin550MeVCM2");
  ParaPerpAsymmPhin_550MeVCM2->SetTitle("Neutron Para/Perp Phi Asymmetry for 540-560MeV Photon Energy (ThetaCM 30-60)");
  ParaPerpAsymmPhin_550MeVCM2->Fit("CosFit", "LL");
  nCosAmp[1][7] = CosFit->GetParameter(0);
  nCosAmpErr[1][7] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_570MeVCM2 = Phin_570MeVCM2_Para->GetAsymmetry(Phin_570MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_570MeVCM2->SetName("ParaPerpAsymmPhin570MeVCM2");
  ParaPerpAsymmPhin_570MeVCM2->SetTitle("Neutron Para/Perp Phi Asymmetry for 560-580MeV Photon Energy (ThetaCM 30-60)");
  ParaPerpAsymmPhin_570MeVCM2->Fit("CosFit", "LL");
  nCosAmp[1][8] = CosFit->GetParameter(0);
  nCosAmpErr[1][8] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_590MeVCM2 = Phin_590MeVCM2_Para->GetAsymmetry(Phin_590MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_590MeVCM2->SetName("ParaPerpAsymmPhin590MeVCM2");
  ParaPerpAsymmPhin_590MeVCM2->SetTitle("Neutron Para/Perp Phi Asymmetry for 580-600MeV Photon Energy (ThetaCM 30-60)");
  ParaPerpAsymmPhin_590MeVCM2->Fit("CosFit", "LL");
  nCosAmp[1][9] = CosFit->GetParameter(0);
  nCosAmpErr[1][9] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_610MeVCM2 = Phin_610MeVCM2_Para->GetAsymmetry(Phin_610MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_610MeVCM2->SetName("ParaPerpAsymmPhin610MeVCM2");
  ParaPerpAsymmPhin_610MeVCM2->SetTitle("Neutron Para/Perp Phi Asymmetry for 600-620MeV Photon Energy (ThetaCM 30-60)");
  ParaPerpAsymmPhin_610MeVCM2->Fit("CosFit", "LL");
  nCosAmp[1][10] = CosFit->GetParameter(0);
  nCosAmpErr[1][10] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_630MeVCM2 = Phin_630MeVCM2_Para->GetAsymmetry(Phin_630MeVCM2_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_630MeVCM2->SetName("ParaPerpAsymmPhin630MeVCM2");
  ParaPerpAsymmPhin_630MeVCM2->SetTitle("Neutron Para/Perp Phi Asymmetry for 620-640MeV Photon Energy (ThetaCM 30-60)");
  ParaPerpAsymmPhin_630MeVCM2->Fit("CosFit", "LL");
  nCosAmp[1][11] = CosFit->GetParameter(0);
  nCosAmpErr[1][11] = CosFit->GetParError(0);

  ///////////////////////////////////////////
  //////////////////CM3 p ///////////////////
  ///////////////////////////////////////////

  ParaPerpAsymmPhin_410MeVCM3 = Phin_410MeVCM3_Para->GetAsymmetry(Phin_410MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_410MeVCM3->SetName("ParaPerpAsymmPhin410MeVCM3");
  ParaPerpAsymmPhin_410MeVCM3->SetTitle("Neutron Para/Perp Phi Asymmetry for 400-420MeV Photon Energy (ThetaCM 60-90)");
  ParaPerpAsymmPhin_410MeVCM3->Fit("CosFit", "LL");
  nCosAmp[2][0] = CosFit->GetParameter(0);
  nCosAmpErr[2][0] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_430MeVCM3 = Phin_430MeVCM3_Para->GetAsymmetry(Phin_430MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_430MeVCM3->SetName("ParaPerpAsymmPhin430MeVCM3");
  ParaPerpAsymmPhin_430MeVCM3->SetTitle("Neutron Para/Perp Phi Asymmetry for 420-440MeV Photon Energy (ThetaCM 60-90)");
  ParaPerpAsymmPhin_430MeVCM3->Fit("CosFit", "LL");
  nCosAmp[2][1] = CosFit->GetParameter(0);
  nCosAmpErr[2][1] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_450MeVCM3 = Phin_450MeVCM3_Para->GetAsymmetry(Phin_450MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_450MeVCM3->SetName("ParaPerpAsymmPhin450MeVCM3");
  ParaPerpAsymmPhin_450MeVCM3->SetTitle("Neutron Para/Perp Phi Asymmetry for 440-460MeV Photon Energy (ThetaCM 60-90)");
  ParaPerpAsymmPhin_450MeVCM3->Fit("CosFit", "LL");
  nCosAmp[2][2] = CosFit->GetParameter(0);
  nCosAmpErr[2][2] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_470MeVCM3 = Phin_470MeVCM3_Para->GetAsymmetry(Phin_470MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_470MeVCM3->SetName("ParaPerpAsymmPhin470MeVCM3");
  ParaPerpAsymmPhin_470MeVCM3->SetTitle("Neutron Para/Perp Phi Asymmetry for 460-480MeV Photon Energy (ThetaCM 60-90)");
  ParaPerpAsymmPhin_470MeVCM3->Fit("CosFit", "LL");
  nCosAmp[2][3] = CosFit->GetParameter(0);
  nCosAmpErr[2][3] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_490MeVCM3 = Phin_490MeVCM3_Para->GetAsymmetry(Phin_490MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_490MeVCM3->SetName("ParaPerpAsymmPhin490MeVCM3");
  ParaPerpAsymmPhin_490MeVCM3->SetTitle("Neutron Para/Perp Phi Asymmetry for 480-500MeV Photon Energy (ThetaCM 60-90)");
  ParaPerpAsymmPhin_490MeVCM3->Fit("CosFit", "LL");
  nCosAmp[2][4] = CosFit->GetParameter(0);
  nCosAmpErr[2][4] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_510MeVCM3 = Phin_510MeVCM3_Para->GetAsymmetry(Phin_510MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_510MeVCM3->SetName("ParaPerpAsymmPhin510MeVCM3");
  ParaPerpAsymmPhin_510MeVCM3->SetTitle("Neutron Para/Perp Phi Asymmetry for 500-520MeV Photon Energy (ThetaCM 60-90)");
  ParaPerpAsymmPhin_510MeVCM3->Fit("CosFit", "LL");
  nCosAmp[2][5] = CosFit->GetParameter(0);
  nCosAmpErr[2][5] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_530MeVCM3 = Phin_530MeVCM3_Para->GetAsymmetry(Phin_530MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_530MeVCM3->SetName("ParaPerpAsymmPhin530MeVCM3");
  ParaPerpAsymmPhin_530MeVCM3->SetTitle("Neutron Para/Perp Phi Asymmetry for 520-540MeV Photon Energy (ThetaCM 60-90)");
  ParaPerpAsymmPhin_530MeVCM3->Fit("CosFit", "LL");
  nCosAmp[2][6] = CosFit->GetParameter(0);
  nCosAmpErr[2][6] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_550MeVCM3 = Phin_550MeVCM3_Para->GetAsymmetry(Phin_550MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_550MeVCM3->SetName("ParaPerpAsymmPhin550MeVCM3");
  ParaPerpAsymmPhin_550MeVCM3->SetTitle("Neutron Para/Perp Phi Asymmetry for 540-560MeV Photon Energy (ThetaCM 60-90)");
  ParaPerpAsymmPhin_550MeVCM3->Fit("CosFit", "LL");
  nCosAmp[2][7] = CosFit->GetParameter(0);
  nCosAmpErr[2][7] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_570MeVCM3 = Phin_570MeVCM3_Para->GetAsymmetry(Phin_570MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_570MeVCM3->SetName("ParaPerpAsymmPhin570MeVCM3");
  ParaPerpAsymmPhin_570MeVCM3->SetTitle("Neutron Para/Perp Phi Asymmetry for 560-580MeV Photon Energy (ThetaCM 60-90)");
  ParaPerpAsymmPhin_570MeVCM3->Fit("CosFit", "LL");
  nCosAmp[2][8] = CosFit->GetParameter(0);
  nCosAmpErr[2][8] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_590MeVCM3 = Phin_590MeVCM3_Para->GetAsymmetry(Phin_590MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_590MeVCM3->SetName("ParaPerpAsymmPhin590MeVCM3");
  ParaPerpAsymmPhin_590MeVCM3->SetTitle("Neutron Para/Perp Phi Asymmetry for 580-600MeV Photon Energy (ThetaCM 60-90)");
  ParaPerpAsymmPhin_590MeVCM3->Fit("CosFit", "LL");
  nCosAmp[2][9] = CosFit->GetParameter(0);
  nCosAmpErr[2][9] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_610MeVCM3 = Phin_610MeVCM3_Para->GetAsymmetry(Phin_610MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_610MeVCM3->SetName("ParaPerpAsymmPhin610MeVCM3");
  ParaPerpAsymmPhin_610MeVCM3->SetTitle("Neutron Para/Perp Phi Asymmetry for 600-620MeV Photon Energy (ThetaCM 60-90)");
  ParaPerpAsymmPhin_610MeVCM3->Fit("CosFit", "LL");
  nCosAmp[2][10] = CosFit->GetParameter(0);
  nCosAmpErr[2][10] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_630MeVCM3 = Phin_630MeVCM3_Para->GetAsymmetry(Phin_630MeVCM3_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_630MeVCM3->SetName("ParaPerpAsymmPhin630MeVCM3");
  ParaPerpAsymmPhin_630MeVCM3->SetTitle("Neutron Para/Perp Phi Asymmetry for 620-640MeV Photon Energy (ThetaCM 60-90)");
  ParaPerpAsymmPhin_630MeVCM3->Fit("CosFit", "LL");
  nCosAmp[2][11] = CosFit->GetParameter(0);
  nCosAmpErr[2][11] = CosFit->GetParError(0);

  ///////////////////////////////////////////
  //////////////////CM4 p ///////////////////
  ///////////////////////////////////////////

  ParaPerpAsymmPhin_410MeVCM4 = Phin_410MeVCM4_Para->GetAsymmetry(Phin_410MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_410MeVCM4->SetName("ParaPerpAsymmPhin410MeVCM4");
  ParaPerpAsymmPhin_410MeVCM4->SetTitle("Neutron Para/Perp Phi Asymmetry for 400-420MeV Photon Energy (ThetaCM 90-120)");
  ParaPerpAsymmPhin_410MeVCM4->Fit("CosFit", "LL");
  nCosAmp[3][0] = CosFit->GetParameter(0);
  nCosAmpErr[3][0] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_430MeVCM4 = Phin_430MeVCM4_Para->GetAsymmetry(Phin_430MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_430MeVCM4->SetName("ParaPerpAsymmPhin430MeVCM4");
  ParaPerpAsymmPhin_430MeVCM4->SetTitle("Neutron Para/Perp Phi Asymmetry for 420-440MeV Photon Energy (ThetaCM 90-120)");
  ParaPerpAsymmPhin_430MeVCM4->Fit("CosFit", "LL");
  nCosAmp[3][1] = CosFit->GetParameter(0);
  nCosAmpErr[3][1] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_450MeVCM4 = Phin_450MeVCM4_Para->GetAsymmetry(Phin_450MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_450MeVCM4->SetName("ParaPerpAsymmPhin450MeVCM4");
  ParaPerpAsymmPhin_450MeVCM4->SetTitle("Neutron Para/Perp Phi Asymmetry for 440-460MeV Photon Energy (ThetaCM 90-120)");
  ParaPerpAsymmPhin_450MeVCM4->Fit("CosFit", "LL");
  nCosAmp[3][2] = CosFit->GetParameter(0);
  nCosAmpErr[3][2] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_470MeVCM4 = Phin_470MeVCM4_Para->GetAsymmetry(Phin_470MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_470MeVCM4->SetName("ParaPerpAsymmPhin470MeVCM4");
  ParaPerpAsymmPhin_470MeVCM4->SetTitle("Neutron Para/Perp Phi Asymmetry for 460-480MeV Photon Energy (ThetaCM 90-120)");
  ParaPerpAsymmPhin_470MeVCM4->Fit("CosFit", "LL");
  nCosAmp[3][3] = CosFit->GetParameter(0);
  nCosAmpErr[3][3] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_490MeVCM4 = Phin_490MeVCM4_Para->GetAsymmetry(Phin_490MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_490MeVCM4->SetName("ParaPerpAsymmPhin490MeVCM4");
  ParaPerpAsymmPhin_490MeVCM4->SetTitle("Neutron Para/Perp Phi Asymmetry for 480-500MeV Photon Energy (ThetaCM 90-120)");
  ParaPerpAsymmPhin_490MeVCM4->Fit("CosFit", "LL");
  nCosAmp[3][4] = CosFit->GetParameter(0);
  nCosAmpErr[3][4] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_510MeVCM4 = Phin_510MeVCM4_Para->GetAsymmetry(Phin_510MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_510MeVCM4->SetName("ParaPerpAsymmPhin510MeVCM4");
  ParaPerpAsymmPhin_510MeVCM4->SetTitle("Neutron Para/Perp Phi Asymmetry for 500-520MeV Photon Energy (ThetaCM 90-120)");
  ParaPerpAsymmPhin_510MeVCM4->Fit("CosFit", "LL");
  nCosAmp[3][5] = CosFit->GetParameter(0);
  nCosAmpErr[3][5] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_530MeVCM4 = Phin_530MeVCM4_Para->GetAsymmetry(Phin_530MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_530MeVCM4->SetName("ParaPerpAsymmPhin530MeVCM4");
  ParaPerpAsymmPhin_530MeVCM4->SetTitle("Neutron Para/Perp Phi Asymmetry for 520-540MeV Photon Energy (ThetaCM 90-120)");
  ParaPerpAsymmPhin_530MeVCM4->Fit("CosFit", "LL");
  nCosAmp[3][6] = CosFit->GetParameter(0);
  nCosAmpErr[3][6] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_550MeVCM4 = Phin_550MeVCM4_Para->GetAsymmetry(Phin_550MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_550MeVCM4->SetName("ParaPerpAsymmPhin550MeVCM4");
  ParaPerpAsymmPhin_550MeVCM4->SetTitle("Neutron Para/Perp Phi Asymmetry for 540-560MeV Photon Energy (ThetaCM 90-120)");
  ParaPerpAsymmPhin_550MeVCM4->Fit("CosFit", "LL");
  nCosAmp[3][7] = CosFit->GetParameter(0);
  nCosAmpErr[3][7] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_570MeVCM4 = Phin_570MeVCM4_Para->GetAsymmetry(Phin_570MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_570MeVCM4->SetName("ParaPerpAsymmPhin570MeVCM4");
  ParaPerpAsymmPhin_570MeVCM4->SetTitle("Neutron Para/Perp Phi Asymmetry for 560-580MeV Photon Energy (ThetaCM 90-120)");
  ParaPerpAsymmPhin_570MeVCM4->Fit("CosFit", "LL");
  nCosAmp[3][8] = CosFit->GetParameter(0);
  nCosAmpErr[3][8] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_590MeVCM4 = Phin_590MeVCM4_Para->GetAsymmetry(Phin_590MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_590MeVCM4->SetName("ParaPerpAsymmPhin590MeVCM4");
  ParaPerpAsymmPhin_590MeVCM4->SetTitle("Neutron Para/Perp Phi Asymmetry for 580-600MeV Photon Energy (ThetaCM 90-120)");
  ParaPerpAsymmPhin_590MeVCM4->Fit("CosFit", "LL");
  nCosAmp[3][9] = CosFit->GetParameter(0);
  nCosAmpErr[3][9] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_610MeVCM4 = Phin_610MeVCM4_Para->GetAsymmetry(Phin_610MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_610MeVCM4->SetName("ParaPerpAsymmPhin610MeVCM4");
  ParaPerpAsymmPhin_610MeVCM4->SetTitle("Neutron Para/Perp Phi Asymmetry for 600-620MeV Photon Energy (ThetaCM 90-120)");
  ParaPerpAsymmPhin_610MeVCM4->Fit("CosFit", "LL");
  nCosAmp[3][10] = CosFit->GetParameter(0);
  nCosAmpErr[3][10] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_630MeVCM4 = Phin_630MeVCM4_Para->GetAsymmetry(Phin_630MeVCM4_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_630MeVCM4->SetName("ParaPerpAsymmPhin630MeVCM4");
  ParaPerpAsymmPhin_630MeVCM4->SetTitle("Neutron Para/Perp Phi Asymmetry for 620-640MeV Photon Energy (ThetaCM 90-120)");
  ParaPerpAsymmPhin_630MeVCM4->Fit("CosFit", "LL");
  nCosAmp[3][11] = CosFit->GetParameter(0);
  nCosAmpErr[3][11] = CosFit->GetParError(0);

  ///////////////////////////////////////////
  //////////////////CM5 p ///////////////////
  ///////////////////////////////////////////

  ParaPerpAsymmPhin_410MeVCM5 = Phin_410MeVCM5_Para->GetAsymmetry(Phin_410MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_410MeVCM5->SetName("ParaPerpAsymmPhin410MeVCM5");
  ParaPerpAsymmPhin_410MeVCM5->SetTitle("Neutron Para/Perp Phi Asymmetry for 400-420MeV Photon Energy (ThetaCM 120-150)");
  ParaPerpAsymmPhin_410MeVCM5->Fit("CosFit", "LL");
  nCosAmp[4][0] = CosFit->GetParameter(0);
  nCosAmpErr[4][0] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_430MeVCM5 = Phin_430MeVCM5_Para->GetAsymmetry(Phin_430MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_430MeVCM5->SetName("ParaPerpAsymmPhin430MeVCM5");
  ParaPerpAsymmPhin_430MeVCM5->SetTitle("Neutron Para/Perp Phi Asymmetry for 420-440MeV Photon Energy (ThetaCM 120-150)");
  ParaPerpAsymmPhin_430MeVCM5->Fit("CosFit", "LL");
  nCosAmp[4][1] = CosFit->GetParameter(0);
  nCosAmpErr[4][1] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_450MeVCM5 = Phin_450MeVCM5_Para->GetAsymmetry(Phin_450MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_450MeVCM5->SetName("ParaPerpAsymmPhin450MeVCM5");
  ParaPerpAsymmPhin_450MeVCM5->SetTitle("Neutron Para/Perp Phi Asymmetry for 440-460MeV Photon Energy (ThetaCM 120-150)");
  ParaPerpAsymmPhin_450MeVCM5->Fit("CosFit", "LL");
  nCosAmp[4][2] = CosFit->GetParameter(0);
  nCosAmpErr[4][2] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_470MeVCM5 = Phin_470MeVCM5_Para->GetAsymmetry(Phin_470MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_470MeVCM5->SetName("ParaPerpAsymmPhin470MeVCM5");
  ParaPerpAsymmPhin_470MeVCM5->SetTitle("Neutron Para/Perp Phi Asymmetry for 460-480MeV Photon Energy (ThetaCM 120-150)");
  ParaPerpAsymmPhin_470MeVCM5->Fit("CosFit", "LL");
  nCosAmp[4][3] = CosFit->GetParameter(0);
  nCosAmpErr[4][3] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_490MeVCM5 = Phin_490MeVCM5_Para->GetAsymmetry(Phin_490MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_490MeVCM5->SetName("ParaPerpAsymmPhin490MeVCM5");
  ParaPerpAsymmPhin_490MeVCM5->SetTitle("Neutron Para/Perp Phi Asymmetry for 480-500MeV Photon Energy (ThetaCM 120-150)");
  ParaPerpAsymmPhin_490MeVCM5->Fit("CosFit", "LL");
  nCosAmp[4][4] = CosFit->GetParameter(0);
  nCosAmpErr[4][4] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_510MeVCM5 = Phin_510MeVCM5_Para->GetAsymmetry(Phin_510MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_510MeVCM5->SetName("ParaPerpAsymmPhin510MeVCM5");
  ParaPerpAsymmPhin_510MeVCM5->SetTitle("Neutron Para/Perp Phi Asymmetry for 500-520MeV Photon Energy (ThetaCM 120-150)");
  ParaPerpAsymmPhin_510MeVCM5->Fit("CosFit", "LL");
  nCosAmp[4][5] = CosFit->GetParameter(0);
  nCosAmpErr[4][5] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_530MeVCM5 = Phin_530MeVCM5_Para->GetAsymmetry(Phin_530MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_530MeVCM5->SetName("ParaPerpAsymmPhin530MeVCM5");
  ParaPerpAsymmPhin_530MeVCM5->SetTitle("Neutron Para/Perp Phi Asymmetry for 520-540MeV Photon Energy (ThetaCM 120-150)");
  ParaPerpAsymmPhin_530MeVCM5->Fit("CosFit", "LL");
  nCosAmp[4][6] = CosFit->GetParameter(0);
  nCosAmpErr[4][6] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_550MeVCM5 = Phin_550MeVCM5_Para->GetAsymmetry(Phin_550MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_550MeVCM5->SetName("ParaPerpAsymmPhin550MeVCM5");
  ParaPerpAsymmPhin_550MeVCM5->SetTitle("Neutron Para/Perp Phi Asymmetry for 540-560MeV Photon Energy (ThetaCM 120-150)");
  ParaPerpAsymmPhin_550MeVCM5->Fit("CosFit", "LL");
  nCosAmp[4][7] = CosFit->GetParameter(0);
  nCosAmpErr[4][7] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_570MeVCM5 = Phin_570MeVCM5_Para->GetAsymmetry(Phin_570MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_570MeVCM5->SetName("ParaPerpAsymmPhin570MeVCM5");
  ParaPerpAsymmPhin_570MeVCM5->SetTitle("Neutron Para/Perp Phi Asymmetry for 560-580MeV Photon Energy (ThetaCM 120-150)");
  ParaPerpAsymmPhin_570MeVCM5->Fit("CosFit", "LL");
  nCosAmp[4][8] = CosFit->GetParameter(0);
  nCosAmpErr[4][8] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_590MeVCM5 = Phin_590MeVCM5_Para->GetAsymmetry(Phin_590MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_590MeVCM5->SetName("ParaPerpAsymmPhin590MeVCM5");
  ParaPerpAsymmPhin_590MeVCM5->SetTitle("Neutron Para/Perp Phi Asymmetry for 580-600MeV Photon Energy (ThetaCM 120-150)");
  ParaPerpAsymmPhin_590MeVCM5->Fit("CosFit", "LL");
  nCosAmp[4][9] = CosFit->GetParameter(0);
  nCosAmpErr[4][9] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_610MeVCM5 = Phin_610MeVCM5_Para->GetAsymmetry(Phin_610MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_610MeVCM5->SetName("ParaPerpAsymmPhin610MeVCM5");
  ParaPerpAsymmPhin_610MeVCM5->SetTitle("Neutron Para/Perp Phi Asymmetry for 600-620MeV Photon Energy (ThetaCM 120-150)");
  ParaPerpAsymmPhin_610MeVCM5->Fit("CosFit", "LL");
  nCosAmp[4][10] = CosFit->GetParameter(0);
  nCosAmpErr[4][10] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_630MeVCM5 = Phin_630MeVCM5_Para->GetAsymmetry(Phin_630MeVCM5_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_630MeVCM5->SetName("ParaPerpAsymmPhin630MeVCM5");
  ParaPerpAsymmPhin_630MeVCM5->SetTitle("Neutron Para/Perp Phi Asymmetry for 620-640MeV Photon Energy (ThetaCM 120-150)");
  ParaPerpAsymmPhin_630MeVCM5->Fit("CosFit", "LL");
  nCosAmp[4][11] = CosFit->GetParameter(0);
  nCosAmpErr[4][11] = CosFit->GetParError(0);

  ///////////////////////////////////////////
  //////////////////CM6 p ///////////////////
  ///////////////////////////////////////////

  ParaPerpAsymmPhin_410MeVCM6 = Phin_410MeVCM6_Para->GetAsymmetry(Phin_410MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_410MeVCM6->SetName("ParaPerpAsymmPhin410MeVCM6");
  ParaPerpAsymmPhin_410MeVCM6->SetTitle("Neutron Para/Perp Phi Asymmetry for 400-420MeV Photon Energy (ThetaCM 150-180)");
  ParaPerpAsymmPhin_410MeVCM6->Fit("CosFit", "LL");
  nCosAmp[5][0] = CosFit->GetParameter(0);
  nCosAmpErr[5][0] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_430MeVCM6 = Phin_430MeVCM6_Para->GetAsymmetry(Phin_430MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_430MeVCM6->SetName("ParaPerpAsymmPhin430MeVCM6");
  ParaPerpAsymmPhin_430MeVCM6->SetTitle("Neutron Para/Perp Phi Asymmetry for 420-440MeV Photon Energy (ThetaCM 150-180)");
  ParaPerpAsymmPhin_430MeVCM6->Fit("CosFit", "LL");
  nCosAmp[5][1] = CosFit->GetParameter(0);
  nCosAmpErr[5][1] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_450MeVCM6 = Phin_450MeVCM6_Para->GetAsymmetry(Phin_450MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_450MeVCM6->SetName("ParaPerpAsymmPhin450MeVCM6");
  ParaPerpAsymmPhin_450MeVCM6->SetTitle("Neutron Para/Perp Phi Asymmetry for 440-460MeV Photon Energy (ThetaCM 150-180)");
  ParaPerpAsymmPhin_450MeVCM6->Fit("CosFit", "LL");
  nCosAmp[5][2] = CosFit->GetParameter(0);
  nCosAmpErr[5][2] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_470MeVCM6 = Phin_470MeVCM6_Para->GetAsymmetry(Phin_470MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_470MeVCM6->SetName("ParaPerpAsymmPhin470MeVCM6");
  ParaPerpAsymmPhin_470MeVCM6->SetTitle("Neutron Para/Perp Phi Asymmetry for 460-480MeV Photon Energy (ThetaCM 150-180)");
  ParaPerpAsymmPhin_470MeVCM6->Fit("CosFit", "LL");
  nCosAmp[5][3] = CosFit->GetParameter(0);
  nCosAmpErr[5][3] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_490MeVCM6 = Phin_490MeVCM6_Para->GetAsymmetry(Phin_490MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_490MeVCM6->SetName("ParaPerpAsymmPhin490MeVCM6");
  ParaPerpAsymmPhin_490MeVCM6->SetTitle("Neutron Para/Perp Phi Asymmetry for 480-500MeV Photon Energy (ThetaCM 150-180)");
  ParaPerpAsymmPhin_490MeVCM6->Fit("CosFit", "LL");
  nCosAmp[5][4] = CosFit->GetParameter(0);
  nCosAmpErr[5][4] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_510MeVCM6 = Phin_510MeVCM6_Para->GetAsymmetry(Phin_510MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_510MeVCM6->SetName("ParaPerpAsymmPhin510MeVCM6");
  ParaPerpAsymmPhin_510MeVCM6->SetTitle("Neutron Para/Perp Phi Asymmetry for 500-520MeV Photon Energy (ThetaCM 150-180)");
  ParaPerpAsymmPhin_510MeVCM6->Fit("CosFit", "LL");
  nCosAmp[5][5] = CosFit->GetParameter(0);
  nCosAmpErr[5][5] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_530MeVCM6 = Phin_530MeVCM6_Para->GetAsymmetry(Phin_530MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_530MeVCM6->SetName("ParaPerpAsymmPhin530MeVCM6");
  ParaPerpAsymmPhin_530MeVCM6->SetTitle("Neutron Para/Perp Phi Asymmetry for 520-540MeV Photon Energy (ThetaCM 150-180)");
  ParaPerpAsymmPhin_530MeVCM6->Fit("CosFit", "LL");
  nCosAmp[5][6] = CosFit->GetParameter(0);
  nCosAmpErr[5][6] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_550MeVCM6 = Phin_550MeVCM6_Para->GetAsymmetry(Phin_550MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_550MeVCM6->SetName("ParaPerpAsymmPhin550MeVCM6");
  ParaPerpAsymmPhin_550MeVCM6->SetTitle("Neutron Para/Perp Phi Asymmetry for 540-560MeV Photon Energy (ThetaCM 150-180)");
  ParaPerpAsymmPhin_550MeVCM6->Fit("CosFit", "LL");
  nCosAmp[5][7] = CosFit->GetParameter(0);
  nCosAmpErr[5][7] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_570MeVCM6 = Phin_570MeVCM6_Para->GetAsymmetry(Phin_570MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_570MeVCM6->SetName("ParaPerpAsymmPhin570MeVCM6");
  ParaPerpAsymmPhin_570MeVCM6->SetTitle("Neutron Para/Perp Phi Asymmetry for 560-580MeV Photon Energy (ThetaCM 150-180)");
  ParaPerpAsymmPhin_570MeVCM6->Fit("CosFit", "LL");
  nCosAmp[5][8] = CosFit->GetParameter(0);
  nCosAmpErr[5][8] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_590MeVCM6 = Phin_590MeVCM6_Para->GetAsymmetry(Phin_590MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_590MeVCM6->SetName("ParaPerpAsymmPhin590MeVCM6");
  ParaPerpAsymmPhin_590MeVCM6->SetTitle("Neutron Para/Perp Phi Asymmetry for 580-600MeV Photon Energy (ThetaCM 150-180)");
  ParaPerpAsymmPhin_590MeVCM6->Fit("CosFit", "LL");
  nCosAmp[5][9] = CosFit->GetParameter(0);
  nCosAmpErr[5][9] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_610MeVCM6 = Phin_610MeVCM6_Para->GetAsymmetry(Phin_610MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_610MeVCM6->SetName("ParaPerpAsymmPhin610MeVCM6");
  ParaPerpAsymmPhin_610MeVCM6->SetTitle("Neutron Para/Perp Phi Asymmetry for 600-620MeV Photon Energy (ThetaCM 150-180)");
  ParaPerpAsymmPhin_610MeVCM6->Fit("CosFit", "LL");
  nCosAmp[5][10] = CosFit->GetParameter(0);
  nCosAmpErr[5][10] = CosFit->GetParError(0);

  ParaPerpAsymmPhin_630MeVCM6 = Phin_630MeVCM6_Para->GetAsymmetry(Phin_630MeVCM6_Perp, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhin_630MeVCM6->SetName("ParaPerpAsymmPhin630MeVCM6");
  ParaPerpAsymmPhin_630MeVCM6->SetTitle("Neutron Para/Perp Phi Asymmetry for 620-640MeV Photon Energy (ThetaCM 150-180)");
  ParaPerpAsymmPhin_630MeVCM6->Fit("CosFit", "LL");
  nCosAmp[5][11] = CosFit->GetParameter(0);
  nCosAmpErr[5][11] = CosFit->GetParError(0);

  TFile f1("ParaPerpAsymm_Total_9.root", "RECREATE");

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

  /////////////////////////////////////
  ////////////// Phin /////////////////
  /////////////////////////////////////

  ParaPerpAsymmPhin_410MeVCM1->Write();
  ParaPerpAsymmPhin_430MeVCM1->Write();
  ParaPerpAsymmPhin_450MeVCM1->Write();
  ParaPerpAsymmPhin_470MeVCM1->Write();
  ParaPerpAsymmPhin_490MeVCM1->Write();
  ParaPerpAsymmPhin_510MeVCM1->Write();
  ParaPerpAsymmPhin_530MeVCM1->Write();
  ParaPerpAsymmPhin_550MeVCM1->Write();
  ParaPerpAsymmPhin_570MeVCM1->Write();
  ParaPerpAsymmPhin_590MeVCM1->Write();
  ParaPerpAsymmPhin_610MeVCM1->Write();
  ParaPerpAsymmPhin_630MeVCM1->Write();
  
  ParaPerpAsymmPhin_410MeVCM2->Write();
  ParaPerpAsymmPhin_430MeVCM2->Write();
  ParaPerpAsymmPhin_450MeVCM2->Write();
  ParaPerpAsymmPhin_470MeVCM2->Write();
  ParaPerpAsymmPhin_490MeVCM2->Write();
  ParaPerpAsymmPhin_510MeVCM2->Write();
  ParaPerpAsymmPhin_530MeVCM2->Write();
  ParaPerpAsymmPhin_550MeVCM2->Write();
  ParaPerpAsymmPhin_570MeVCM2->Write();
  ParaPerpAsymmPhin_590MeVCM2->Write();
  ParaPerpAsymmPhin_610MeVCM2->Write();
  ParaPerpAsymmPhin_630MeVCM2->Write();

  ParaPerpAsymmPhin_410MeVCM3->Write();
  ParaPerpAsymmPhin_430MeVCM3->Write();
  ParaPerpAsymmPhin_450MeVCM3->Write();
  ParaPerpAsymmPhin_470MeVCM3->Write();
  ParaPerpAsymmPhin_490MeVCM3->Write();
  ParaPerpAsymmPhin_510MeVCM3->Write();
  ParaPerpAsymmPhin_530MeVCM3->Write();
  ParaPerpAsymmPhin_550MeVCM3->Write();
  ParaPerpAsymmPhin_570MeVCM3->Write();
  ParaPerpAsymmPhin_590MeVCM3->Write();
  ParaPerpAsymmPhin_610MeVCM3->Write();
  ParaPerpAsymmPhin_630MeVCM3->Write();

  ParaPerpAsymmPhin_410MeVCM4->Write();
  ParaPerpAsymmPhin_430MeVCM4->Write();
  ParaPerpAsymmPhin_450MeVCM4->Write();
  ParaPerpAsymmPhin_470MeVCM4->Write();
  ParaPerpAsymmPhin_490MeVCM4->Write();
  ParaPerpAsymmPhin_510MeVCM4->Write();
  ParaPerpAsymmPhin_530MeVCM4->Write();
  ParaPerpAsymmPhin_550MeVCM4->Write();
  ParaPerpAsymmPhin_570MeVCM4->Write();
  ParaPerpAsymmPhin_590MeVCM4->Write();
  ParaPerpAsymmPhin_610MeVCM4->Write();
  ParaPerpAsymmPhin_630MeVCM4->Write();

  ParaPerpAsymmPhin_410MeVCM5->Write();
  ParaPerpAsymmPhin_430MeVCM5->Write();
  ParaPerpAsymmPhin_450MeVCM5->Write();
  ParaPerpAsymmPhin_470MeVCM5->Write();
  ParaPerpAsymmPhin_490MeVCM5->Write();
  ParaPerpAsymmPhin_510MeVCM5->Write();
  ParaPerpAsymmPhin_530MeVCM5->Write();
  ParaPerpAsymmPhin_550MeVCM5->Write();
  ParaPerpAsymmPhin_570MeVCM5->Write();
  ParaPerpAsymmPhin_590MeVCM5->Write();
  ParaPerpAsymmPhin_610MeVCM5->Write();
  ParaPerpAsymmPhin_630MeVCM5->Write();

  ParaPerpAsymmPhin_410MeVCM6->Write();
  ParaPerpAsymmPhin_430MeVCM6->Write();
  ParaPerpAsymmPhin_450MeVCM6->Write();
  ParaPerpAsymmPhin_470MeVCM6->Write();
  ParaPerpAsymmPhin_490MeVCM6->Write();
  ParaPerpAsymmPhin_510MeVCM6->Write();
  ParaPerpAsymmPhin_530MeVCM6->Write();
  ParaPerpAsymmPhin_550MeVCM6->Write();
  ParaPerpAsymmPhin_570MeVCM6->Write();
  ParaPerpAsymmPhin_590MeVCM6->Write();
  ParaPerpAsymmPhin_610MeVCM6->Write();
  ParaPerpAsymmPhin_630MeVCM6->Write();

  //Define new tree to store parameters in
  TTree* tree = new TTree("Parameter_Values", "Tree_of_Values");
  
  // Define branches to store parameters, (Branch Name, Variable, Type of Variable)
  tree->Branch("pCosAmp", &pCosA, "pCosA/D");
  tree->Branch("pCosAmpErr", &pCosAErr, "pCosAErr/D");
  tree->Branch("nCosAmp", &nCosA, "nCosA/D");
  tree->Branch("nCosAmpErr", &nCosAErr, "nCosAErr/D");
  
  // Fill branches (and hence tree) with corresponding parameters from above
  for (Int_t i = 0; i < 6; i++){
    for (Int_t m = 0; m < 12; m++){
    
      pCosA = pCosAmp[i][m];
      pCosAErr = pCosAmpErr[i][m];
      nCosA = nCosAmp[i][m];
      nCosAErr = nCosAmpErr[i][m];
    
      tree->Fill();
    
    }
  }

  f1.Write();

}
