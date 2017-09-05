#include "./includes_SigmaAsymm_NoScatt.h"

void SigmaAsymm_NoScatt_Helicity_Test(){

  TF1 *CosFunc = new TF1("CosFit", "[0]*cos((2*x*TMath::DegToRad())+acos(0))");
  CosFunc->SetParNames("Amplitude");

  double pCosAmpNegHel[10][20]; // Format of array is Theta bin (x) by Egamma bin (y), 6 theta bins of 30, 12 20MeV Egamma bins
  double pCosAmpErrNegHel[10][20];
  double pCosA425NegHel;
  double pCosAErr425NegHel;
  double pCosA435NegHel;
  double pCosAErr435NegHel;
  double pCosA445NegHel;
  double pCosAErr445NegHel;
  double pCosA455NegHel;
  double pCosAErr455NegHel;
  double pCosA465NegHel;
  double pCosAErr465NegHel;
  double pCosA475NegHel;
  double pCosAErr475NegHel;
  double pCosA485NegHel;
  double pCosAErr485NegHel;
  double pCosA495NegHel;
  double pCosAErr495NegHel;
  double pCosA505NegHel;
  double pCosAErr505NegHel;
  double pCosA515NegHel;
  double pCosAErr515NegHel;
  double pCosA525NegHel;
  double pCosAErr525NegHel;
  double pCosA535NegHel;
  double pCosAErr535NegHel;
  double pCosA545NegHel;
  double pCosAErr545NegHel;
  double pCosA555NegHel;
  double pCosAErr555NegHel;
  double pCosA565NegHel;
  double pCosAErr565NegHel;
  double pCosA575NegHel;
  double pCosAErr575NegHel;
  double pCosA585NegHel;
  double pCosAErr585NegHel;
  double pCosA595NegHel;
  double pCosAErr595NegHel;
  double pCosA605NegHel;
  double pCosAErr605NegHel;
  double pCosA615NegHel;
  double pCosAErr615NegHel;
  double pCosAmpPosHel[10][20]; // Format of array is Theta bin (x) by Egamma bin (y), 6 theta bins of 30, 12 20MeV Egamma bins
  double pCosAmpErrPosHel[10][20];
  double pCosA425PosHel;
  double pCosAErr425PosHel;
  double pCosA435PosHel;
  double pCosAErr435PosHel;
  double pCosA445PosHel;
  double pCosAErr445PosHel;
  double pCosA455PosHel;
  double pCosAErr455PosHel;
  double pCosA465PosHel;
  double pCosAErr465PosHel;
  double pCosA475PosHel;
  double pCosAErr475PosHel;
  double pCosA485PosHel;
  double pCosAErr485PosHel;
  double pCosA495PosHel;
  double pCosAErr495PosHel;
  double pCosA505PosHel;
  double pCosAErr505PosHel;
  double pCosA515PosHel;
  double pCosAErr515PosHel;
  double pCosA525PosHel;
  double pCosAErr525PosHel;
  double pCosA535PosHel;
  double pCosAErr535PosHel;
  double pCosA545PosHel;
  double pCosAErr545PosHel;
  double pCosA555PosHel;
  double pCosAErr555PosHel;
  double pCosA565PosHel;
  double pCosAErr565PosHel;
  double pCosA575PosHel;
  double pCosAErr575PosHel;
  double pCosA585PosHel;
  double pCosAErr585PosHel;
  double pCosA595PosHel;
  double pCosAErr595PosHel;
  double pCosA605PosHel;
  double pCosAErr605PosHel;
  double pCosA615PosHel;
  double pCosAErr615PosHel;

  TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/ParaPerp_NoScatt_Total_10_Combined.root"); // Open the latest PTotal combined file to load histograms from
  NPara = Eg_Para->GetEntries();
  NPerp = Eg_Perp->GetEntries();
  ScaleFactor = NPara/NPerp;
  ScaleFactorErr = sqrt( (NPara/((TMath::Power(NPerp,2)))) + (((TMath::Power(NPara,2)))/(TMath::Power(NPerp,3))) ); // Error Propagation of above formula, see notebook 5

  ParaPerpAsymmPhip_425MeVCM1_NegHel = Phip_425MeVCM1_Para_NegHel->GetAsymmetry(Phip_425MeVCM1_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_425MeVCM1_NegHel->SetName("ParaPerpAsymmPhip425MeVCM1");
  ParaPerpAsymmPhip_425MeVCM1_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 420-430MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_425MeVCM1_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[0][0] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[0][0] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_435MeVCM1_NegHel = Phip_435MeVCM1_Para_NegHel->GetAsymmetry(Phip_435MeVCM1_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_435MeVCM1_NegHel->SetName("ParaPerpAsymmPhip435MeVCM1");
  ParaPerpAsymmPhip_435MeVCM1_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 430-440MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_435MeVCM1_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[0][1] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[0][1] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_445MeVCM1_NegHel = Phip_445MeVCM1_Para_NegHel->GetAsymmetry(Phip_445MeVCM1_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_445MeVCM1_NegHel->SetName("ParaPerpAsymmPhip445MeVCM1");
  ParaPerpAsymmPhip_445MeVCM1_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 440-450MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_445MeVCM1_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[0][2] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[0][2] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_455MeVCM1_NegHel = Phip_455MeVCM1_Para_NegHel->GetAsymmetry(Phip_455MeVCM1_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_455MeVCM1_NegHel->SetName("ParaPerpAsymmPhip455MeVCM1");
  ParaPerpAsymmPhip_455MeVCM1_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 450-460MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_455MeVCM1_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[0][3] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[0][3] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_465MeVCM1_NegHel = Phip_465MeVCM1_Para_NegHel->GetAsymmetry(Phip_465MeVCM1_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_465MeVCM1_NegHel->SetName("ParaPerpAsymmPhip465MeVCM1");
  ParaPerpAsymmPhip_465MeVCM1_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 460-470MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_465MeVCM1_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[0][4] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[0][4] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_475MeVCM1_NegHel = Phip_475MeVCM1_Para_NegHel->GetAsymmetry(Phip_475MeVCM1_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_475MeVCM1_NegHel->SetName("ParaPerpAsymmPhip475MeVCM1");
  ParaPerpAsymmPhip_475MeVCM1_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 470-480MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_475MeVCM1_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[0][5] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[0][5] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_485MeVCM1_NegHel = Phip_485MeVCM1_Para_NegHel->GetAsymmetry(Phip_485MeVCM1_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_485MeVCM1_NegHel->SetName("ParaPerpAsymmPhip485MeVCM1");
  ParaPerpAsymmPhip_485MeVCM1_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 480-490MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_485MeVCM1_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[0][6] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[0][6] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_495MeVCM1_NegHel = Phip_495MeVCM1_Para_NegHel->GetAsymmetry(Phip_495MeVCM1_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_495MeVCM1_NegHel->SetName("ParaPerpAsymmPhip495MeVCM1");
  ParaPerpAsymmPhip_495MeVCM1_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 490-500MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_495MeVCM1_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[0][7] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[0][7] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_505MeVCM1_NegHel = Phip_505MeVCM1_Para_NegHel->GetAsymmetry(Phip_505MeVCM1_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_505MeVCM1_NegHel->SetName("ParaPerpAsymmPhip505MeVCM1");
  ParaPerpAsymmPhip_505MeVCM1_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 500-510MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_505MeVCM1_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[0][8] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[0][8] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_515MeVCM1_NegHel = Phip_515MeVCM1_Para_NegHel->GetAsymmetry(Phip_515MeVCM1_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_515MeVCM1_NegHel->SetName("ParaPerpAsymmPhip515MeVCM1");
  ParaPerpAsymmPhip_515MeVCM1_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 510-520MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_515MeVCM1_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[0][9] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[0][9] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_525MeVCM1_NegHel = Phip_525MeVCM1_Para_NegHel->GetAsymmetry(Phip_525MeVCM1_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_525MeVCM1_NegHel->SetName("ParaPerpAsymmPhip525MeVCM1");
  ParaPerpAsymmPhip_525MeVCM1_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 520-530MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_525MeVCM1_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[0][10] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[0][10] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_535MeVCM1_NegHel = Phip_535MeVCM1_Para_NegHel->GetAsymmetry(Phip_535MeVCM1_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_535MeVCM1_NegHel->SetName("ParaPerpAsymmPhip535MeVCM1");
  ParaPerpAsymmPhip_535MeVCM1_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 530-540MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_535MeVCM1_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[0][11] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[0][11] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_545MeVCM1_NegHel = Phip_545MeVCM1_Para_NegHel->GetAsymmetry(Phip_545MeVCM1_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_545MeVCM1_NegHel->SetName("ParaPerpAsymmPhip545MeVCM1");
  ParaPerpAsymmPhip_545MeVCM1_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 540-550MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_545MeVCM1_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[0][12] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[0][12] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_555MeVCM1_NegHel = Phip_555MeVCM1_Para_NegHel->GetAsymmetry(Phip_555MeVCM1_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_555MeVCM1_NegHel->SetName("ParaPerpAsymmPhip555MeVCM1");
  ParaPerpAsymmPhip_555MeVCM1_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 550-560MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_555MeVCM1_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[0][13] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[0][13] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_565MeVCM1_NegHel = Phip_565MeVCM1_Para_NegHel->GetAsymmetry(Phip_565MeVCM1_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_565MeVCM1_NegHel->SetName("ParaPerpAsymmPhip565MeVCM1");
  ParaPerpAsymmPhip_565MeVCM1_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 560-570MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_565MeVCM1_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[0][14] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[0][14] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_575MeVCM1_NegHel = Phip_575MeVCM1_Para_NegHel->GetAsymmetry(Phip_575MeVCM1_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_575MeVCM1_NegHel->SetName("ParaPerpAsymmPhip575MeVCM1");
  ParaPerpAsymmPhip_575MeVCM1_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 570-580MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_575MeVCM1_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[0][15] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[0][15] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_585MeVCM1_NegHel = Phip_585MeVCM1_Para_NegHel->GetAsymmetry(Phip_585MeVCM1_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_585MeVCM1_NegHel->SetName("ParaPerpAsymmPhip585MeVCM1");
  ParaPerpAsymmPhip_585MeVCM1_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 580-590MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_585MeVCM1_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[0][16] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[0][16] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_595MeVCM1_NegHel = Phip_595MeVCM1_Para_NegHel->GetAsymmetry(Phip_595MeVCM1_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_595MeVCM1_NegHel->SetName("ParaPerpAsymmPhip595MeVCM1");
  ParaPerpAsymmPhip_595MeVCM1_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 590-600MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_595MeVCM1_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[0][17] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[0][17] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_605MeVCM1_NegHel = Phip_605MeVCM1_Para_NegHel->GetAsymmetry(Phip_605MeVCM1_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_605MeVCM1_NegHel->SetName("ParaPerpAsymmPhip605MeVCM1");
  ParaPerpAsymmPhip_605MeVCM1_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 600-610MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_605MeVCM1_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[0][18] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[0][18] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_615MeVCM1_NegHel = Phip_615MeVCM1_Para_NegHel->GetAsymmetry(Phip_615MeVCM1_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_615MeVCM1_NegHel->SetName("ParaPerpAsymmPhip615MeVCM1");
  ParaPerpAsymmPhip_615MeVCM1_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 610-620MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_615MeVCM1_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[0][19] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[0][19] = CosFit->GetParError(0);

  ///////////////////////////////////////////
  //////////////////CM2 p ///////////////////
  ///////////////////////////////////////////

  ParaPerpAsymmPhip_425MeVCM2_NegHel = Phip_425MeVCM2_Para_NegHel->GetAsymmetry(Phip_425MeVCM2_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_425MeVCM2_NegHel->SetName("ParaPerpAsymmPhip425MeVCM2");
  ParaPerpAsymmPhip_425MeVCM2_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 420-430MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_425MeVCM2_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[1][0] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[1][0] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_435MeVCM2_NegHel = Phip_435MeVCM2_Para_NegHel->GetAsymmetry(Phip_435MeVCM2_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_435MeVCM2_NegHel->SetName("ParaPerpAsymmPhip435MeVCM2");
  ParaPerpAsymmPhip_435MeVCM2_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 430-440MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_435MeVCM2_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[1][1] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[1][1] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_445MeVCM2_NegHel = Phip_445MeVCM2_Para_NegHel->GetAsymmetry(Phip_445MeVCM2_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_445MeVCM2_NegHel->SetName("ParaPerpAsymmPhip445MeVCM2");
  ParaPerpAsymmPhip_445MeVCM2_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 440-450MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_445MeVCM2_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[1][2] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[1][2] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_455MeVCM2_NegHel = Phip_455MeVCM2_Para_NegHel->GetAsymmetry(Phip_455MeVCM2_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_455MeVCM2_NegHel->SetName("ParaPerpAsymmPhip455MeVCM2");
  ParaPerpAsymmPhip_455MeVCM2_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 450-460MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_455MeVCM2_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[1][3] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[1][3] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_465MeVCM2_NegHel = Phip_465MeVCM2_Para_NegHel->GetAsymmetry(Phip_465MeVCM2_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_465MeVCM2_NegHel->SetName("ParaPerpAsymmPhip465MeVCM2");
  ParaPerpAsymmPhip_465MeVCM2_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 460-470MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_465MeVCM2_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[1][4] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[1][4] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_475MeVCM2_NegHel = Phip_475MeVCM2_Para_NegHel->GetAsymmetry(Phip_475MeVCM2_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_475MeVCM2_NegHel->SetName("ParaPerpAsymmPhip475MeVCM2");
  ParaPerpAsymmPhip_475MeVCM2_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 470-480MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_475MeVCM2_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[1][5] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[1][5] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_485MeVCM2_NegHel = Phip_485MeVCM2_Para_NegHel->GetAsymmetry(Phip_485MeVCM2_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_485MeVCM2_NegHel->SetName("ParaPerpAsymmPhip485MeVCM2");
  ParaPerpAsymmPhip_485MeVCM2_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 480-490MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_485MeVCM2_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[1][6] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[1][6] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_495MeVCM2_NegHel = Phip_495MeVCM2_Para_NegHel->GetAsymmetry(Phip_495MeVCM2_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_495MeVCM2_NegHel->SetName("ParaPerpAsymmPhip495MeVCM2");
  ParaPerpAsymmPhip_495MeVCM2_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 490-500MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_495MeVCM2_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[1][7] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[1][7] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_505MeVCM2_NegHel = Phip_505MeVCM2_Para_NegHel->GetAsymmetry(Phip_505MeVCM2_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_505MeVCM2_NegHel->SetName("ParaPerpAsymmPhip505MeVCM2");
  ParaPerpAsymmPhip_505MeVCM2_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 500-510MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_505MeVCM2_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[1][8] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[1][8] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_515MeVCM2_NegHel = Phip_515MeVCM2_Para_NegHel->GetAsymmetry(Phip_515MeVCM2_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_515MeVCM2_NegHel->SetName("ParaPerpAsymmPhip515MeVCM2");
  ParaPerpAsymmPhip_515MeVCM2_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 510-520MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_515MeVCM2_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[1][9] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[1][9] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_525MeVCM2_NegHel = Phip_525MeVCM2_Para_NegHel->GetAsymmetry(Phip_525MeVCM2_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_525MeVCM2_NegHel->SetName("ParaPerpAsymmPhip525MeVCM2");
  ParaPerpAsymmPhip_525MeVCM2_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 520-530MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_525MeVCM2_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[1][10] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[1][10] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_535MeVCM2_NegHel = Phip_535MeVCM2_Para_NegHel->GetAsymmetry(Phip_535MeVCM2_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_535MeVCM2_NegHel->SetName("ParaPerpAsymmPhip535MeVCM2");
  ParaPerpAsymmPhip_535MeVCM2_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 530-540MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_535MeVCM2_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[1][11] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[1][11] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_545MeVCM2_NegHel = Phip_545MeVCM2_Para_NegHel->GetAsymmetry(Phip_545MeVCM2_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_545MeVCM2_NegHel->SetName("ParaPerpAsymmPhip545MeVCM2");
  ParaPerpAsymmPhip_545MeVCM2_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 540-550MeV Photon Energy -ve Hel (CosTheta0.8-0.6)");
  ParaPerpAsymmPhip_545MeVCM2_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[1][12] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[1][12] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_555MeVCM2_NegHel = Phip_555MeVCM2_Para_NegHel->GetAsymmetry(Phip_555MeVCM2_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_555MeVCM2_NegHel->SetName("ParaPerpAsymmPhip555MeVCM2");
  ParaPerpAsymmPhip_555MeVCM2_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 550-560MeV Photon Energy -ve Hel (CosTheta0.8-0.6)");
  ParaPerpAsymmPhip_555MeVCM2_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[1][13] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[1][13] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_565MeVCM2_NegHel = Phip_565MeVCM2_Para_NegHel->GetAsymmetry(Phip_565MeVCM2_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_565MeVCM2_NegHel->SetName("ParaPerpAsymmPhip565MeVCM2");
  ParaPerpAsymmPhip_565MeVCM2_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 560-570MeV Photon Energy -ve Hel (CosTheta0.8-0.6)");
  ParaPerpAsymmPhip_565MeVCM2_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[1][14] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[1][14] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_575MeVCM2_NegHel = Phip_575MeVCM2_Para_NegHel->GetAsymmetry(Phip_575MeVCM2_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_575MeVCM2_NegHel->SetName("ParaPerpAsymmPhip575MeVCM2");
  ParaPerpAsymmPhip_575MeVCM2_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 570-580MeV Photon Energy -ve Hel (CosTheta0.8-0.6)");
  ParaPerpAsymmPhip_575MeVCM2_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[1][15] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[1][15] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_585MeVCM2_NegHel = Phip_585MeVCM2_Para_NegHel->GetAsymmetry(Phip_585MeVCM2_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_585MeVCM2_NegHel->SetName("ParaPerpAsymmPhip585MeVCM2");
  ParaPerpAsymmPhip_585MeVCM2_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 580-590MeV Photon Energy -ve Hel (CosTheta0.8-0.6)");
  ParaPerpAsymmPhip_585MeVCM2_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[1][16] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[1][16] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_595MeVCM2_NegHel = Phip_595MeVCM2_Para_NegHel->GetAsymmetry(Phip_595MeVCM2_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_595MeVCM2_NegHel->SetName("ParaPerpAsymmPhip595MeVCM2");
  ParaPerpAsymmPhip_595MeVCM2_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 590-600MeV Photon Energy -ve Hel (CosTheta0.8-0.6)");
  ParaPerpAsymmPhip_595MeVCM2_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[1][17] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[1][17] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_605MeVCM2_NegHel = Phip_605MeVCM2_Para_NegHel->GetAsymmetry(Phip_605MeVCM2_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_605MeVCM2_NegHel->SetName("ParaPerpAsymmPhip605MeVCM2");
  ParaPerpAsymmPhip_605MeVCM2_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 600-610MeV Photon Energy -ve Hel (CosTheta0.8-0.6)");
  ParaPerpAsymmPhip_605MeVCM2_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[1][18] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[1][18] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_615MeVCM2_NegHel = Phip_615MeVCM2_Para_NegHel->GetAsymmetry(Phip_615MeVCM2_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_615MeVCM2_NegHel->SetName("ParaPerpAsymmPhip615MeVCM2");
  ParaPerpAsymmPhip_615MeVCM2_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 610-620MeV Photon Energy -ve Hel (CosTheta0.8-0.6)");
  ParaPerpAsymmPhip_615MeVCM2_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[1][19] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[1][19] = CosFit->GetParError(0);

  ///////////////////////////////////////////
  //////////////////CM3 p ///////////////////
  ///////////////////////////////////////////

  ParaPerpAsymmPhip_425MeVCM3_NegHel = Phip_425MeVCM3_Para_NegHel->GetAsymmetry(Phip_425MeVCM3_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_425MeVCM3_NegHel->SetName("ParaPerpAsymmPhip425MeVCM3");
  ParaPerpAsymmPhip_425MeVCM3_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 420-430MeV Photon Energy -ve Hel (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_425MeVCM3_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[2][0] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[2][0] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_435MeVCM3_NegHel = Phip_435MeVCM3_Para_NegHel->GetAsymmetry(Phip_435MeVCM3_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_435MeVCM3_NegHel->SetName("ParaPerpAsymmPhip435MeVCM3");
  ParaPerpAsymmPhip_435MeVCM3_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 430-440MeV Photon Energy -ve Hel (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_435MeVCM3_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[2][1] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[2][1] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_445MeVCM3_NegHel = Phip_445MeVCM3_Para_NegHel->GetAsymmetry(Phip_445MeVCM3_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_445MeVCM3_NegHel->SetName("ParaPerpAsymmPhip445MeVCM3");
  ParaPerpAsymmPhip_445MeVCM3_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 440-450MeV Photon Energy -ve Hel (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_445MeVCM3_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[2][2] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[2][2] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_455MeVCM3_NegHel = Phip_455MeVCM3_Para_NegHel->GetAsymmetry(Phip_455MeVCM3_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_455MeVCM3_NegHel->SetName("ParaPerpAsymmPhip455MeVCM3");
  ParaPerpAsymmPhip_455MeVCM3_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 450-460MeV Photon Energy -ve Hel (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_455MeVCM3_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[2][3] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[2][3] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_465MeVCM3_NegHel = Phip_465MeVCM3_Para_NegHel->GetAsymmetry(Phip_465MeVCM3_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_465MeVCM3_NegHel->SetName("ParaPerpAsymmPhip465MeVCM3");
  ParaPerpAsymmPhip_465MeVCM3_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 460-470MeV Photon Energy -ve Hel (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_465MeVCM3_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[2][4] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[2][4] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_475MeVCM3_NegHel = Phip_475MeVCM3_Para_NegHel->GetAsymmetry(Phip_475MeVCM3_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_475MeVCM3_NegHel->SetName("ParaPerpAsymmPhip475MeVCM3");
  ParaPerpAsymmPhip_475MeVCM3_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 470-480MeV Photon Energy -ve Hel (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_475MeVCM3_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[2][5] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[2][5] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_485MeVCM3_NegHel = Phip_485MeVCM3_Para_NegHel->GetAsymmetry(Phip_485MeVCM3_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_485MeVCM3_NegHel->SetName("ParaPerpAsymmPhip485MeVCM3");
  ParaPerpAsymmPhip_485MeVCM3_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 480-490MeV Photon Energy -ve Hel (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_485MeVCM3_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[2][6] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[2][6] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_495MeVCM3_NegHel = Phip_495MeVCM3_Para_NegHel->GetAsymmetry(Phip_495MeVCM3_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_495MeVCM3_NegHel->SetName("ParaPerpAsymmPhip495MeVCM3");
  ParaPerpAsymmPhip_495MeVCM3_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 490-500MeV Photon Energy -ve Hel (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_495MeVCM3_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[2][7] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[2][7] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_505MeVCM3_NegHel = Phip_505MeVCM3_Para_NegHel->GetAsymmetry(Phip_505MeVCM3_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_505MeVCM3_NegHel->SetName("ParaPerpAsymmPhip505MeVCM3");
  ParaPerpAsymmPhip_505MeVCM3_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 500-510MeV Photon Energy -ve Hel (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_505MeVCM3_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[2][8] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[2][8] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_515MeVCM3_NegHel = Phip_515MeVCM3_Para_NegHel->GetAsymmetry(Phip_515MeVCM3_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_515MeVCM3_NegHel->SetName("ParaPerpAsymmPhip515MeVCM3");
  ParaPerpAsymmPhip_515MeVCM3_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 510-520MeV Photon Energy -ve Hel (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_515MeVCM3_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[2][9] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[2][9] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_525MeVCM3_NegHel = Phip_525MeVCM3_Para_NegHel->GetAsymmetry(Phip_525MeVCM3_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_525MeVCM3_NegHel->SetName("ParaPerpAsymmPhip525MeVCM3");
  ParaPerpAsymmPhip_525MeVCM3_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 520-530MeV Photon Energy -ve Hel (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_525MeVCM3_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[2][10] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[2][10] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_535MeVCM3_NegHel = Phip_535MeVCM3_Para_NegHel->GetAsymmetry(Phip_535MeVCM3_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_535MeVCM3_NegHel->SetName("ParaPerpAsymmPhip535MeVCM3");
  ParaPerpAsymmPhip_535MeVCM3_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 530-540MeV Photon Energy -ve Hel (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_535MeVCM3_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[2][11] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[2][11] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_545MeVCM3_NegHel = Phip_545MeVCM3_Para_NegHel->GetAsymmetry(Phip_545MeVCM3_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_545MeVCM3_NegHel->SetName("ParaPerpAsymmPhip545MeVCM3");
  ParaPerpAsymmPhip_545MeVCM3_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 540-550MeV Photon Energy -ve Hel (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_545MeVCM3_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[2][12] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[2][12] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_555MeVCM3_NegHel = Phip_555MeVCM3_Para_NegHel->GetAsymmetry(Phip_555MeVCM3_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_555MeVCM3_NegHel->SetName("ParaPerpAsymmPhip555MeVCM3");
  ParaPerpAsymmPhip_555MeVCM3_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 550-560MeV Photon Energy -ve Hel (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_555MeVCM3_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[2][13] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[2][13] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_565MeVCM3_NegHel = Phip_565MeVCM3_Para_NegHel->GetAsymmetry(Phip_565MeVCM3_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_565MeVCM3_NegHel->SetName("ParaPerpAsymmPhip565MeVCM3");
  ParaPerpAsymmPhip_565MeVCM3_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 560-570MeV Photon Energy -ve Hel (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_565MeVCM3_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[2][14] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[2][14] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_575MeVCM3_NegHel = Phip_575MeVCM3_Para_NegHel->GetAsymmetry(Phip_575MeVCM3_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_575MeVCM3_NegHel->SetName("ParaPerpAsymmPhip575MeVCM3");
  ParaPerpAsymmPhip_575MeVCM3_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 570-580MeV Photon Energy -ve Hel (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_575MeVCM3_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[2][15] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[2][15] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_585MeVCM3_NegHel = Phip_585MeVCM3_Para_NegHel->GetAsymmetry(Phip_585MeVCM3_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_585MeVCM3_NegHel->SetName("ParaPerpAsymmPhip585MeVCM3");
  ParaPerpAsymmPhip_585MeVCM3_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 580-590MeV Photon Energy -ve Hel (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_585MeVCM3_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[2][16] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[2][16] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_595MeVCM3_NegHel = Phip_595MeVCM3_Para_NegHel->GetAsymmetry(Phip_595MeVCM3_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_595MeVCM3_NegHel->SetName("ParaPerpAsymmPhip595MeVCM3");
  ParaPerpAsymmPhip_595MeVCM3_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 590-600MeV Photon Energy -ve Hel (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_595MeVCM3_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[2][17] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[2][17] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_605MeVCM3_NegHel = Phip_605MeVCM3_Para_NegHel->GetAsymmetry(Phip_605MeVCM3_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_605MeVCM3_NegHel->SetName("ParaPerpAsymmPhip605MeVCM3");
  ParaPerpAsymmPhip_605MeVCM3_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 600-610MeV Photon Energy -ve Hel (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_605MeVCM3_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[2][18] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[2][18] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_615MeVCM3_NegHel = Phip_615MeVCM3_Para_NegHel->GetAsymmetry(Phip_615MeVCM3_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_615MeVCM3_NegHel->SetName("ParaPerpAsymmPhip615MeVCM3");
  ParaPerpAsymmPhip_615MeVCM3_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 610-620MeV Photon Energy -ve Hel (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_615MeVCM3_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[2][19] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[2][19] = CosFit->GetParError(0);

  ///////////////////////////////////////////
  //////////////////CM4 p ///////////////////
  ///////////////////////////////////////////

  ParaPerpAsymmPhip_425MeVCM4_NegHel = Phip_425MeVCM4_Para_NegHel->GetAsymmetry(Phip_425MeVCM4_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_425MeVCM4_NegHel->SetName("ParaPerpAsymmPhip425MeVCM4");
  ParaPerpAsymmPhip_425MeVCM4_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 420-430MeV Photon Energy -ve Hel (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_425MeVCM4_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[3][0] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[3][0] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_435MeVCM4_NegHel = Phip_435MeVCM4_Para_NegHel->GetAsymmetry(Phip_435MeVCM4_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_435MeVCM4_NegHel->SetName("ParaPerpAsymmPhip435MeVCM4");
  ParaPerpAsymmPhip_435MeVCM4_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 430-440MeV Photon Energy -ve Hel (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_435MeVCM4_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[3][1] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[3][1] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_445MeVCM4_NegHel = Phip_445MeVCM4_Para_NegHel->GetAsymmetry(Phip_445MeVCM4_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_445MeVCM4_NegHel->SetName("ParaPerpAsymmPhip445MeVCM4");
  ParaPerpAsymmPhip_445MeVCM4_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 440-450MeV Photon Energy -ve Hel (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_445MeVCM4_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[3][2] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[3][2] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_455MeVCM4_NegHel = Phip_455MeVCM4_Para_NegHel->GetAsymmetry(Phip_455MeVCM4_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_455MeVCM4_NegHel->SetName("ParaPerpAsymmPhip455MeVCM4");
  ParaPerpAsymmPhip_455MeVCM4_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 450-460MeV Photon Energy -ve Hel (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_455MeVCM4_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[3][3] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[3][3] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_465MeVCM4_NegHel = Phip_465MeVCM4_Para_NegHel->GetAsymmetry(Phip_465MeVCM4_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_465MeVCM4_NegHel->SetName("ParaPerpAsymmPhip465MeVCM4");
  ParaPerpAsymmPhip_465MeVCM4_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 460-470MeV Photon Energy -ve Hel (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_465MeVCM4_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[3][4] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[3][4] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_475MeVCM4_NegHel = Phip_475MeVCM4_Para_NegHel->GetAsymmetry(Phip_475MeVCM4_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_475MeVCM4_NegHel->SetName("ParaPerpAsymmPhip475MeVCM4");
  ParaPerpAsymmPhip_475MeVCM4_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 470-480MeV Photon Energy -ve Hel (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_475MeVCM4_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[3][5] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[3][5] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_485MeVCM4_NegHel = Phip_485MeVCM4_Para_NegHel->GetAsymmetry(Phip_485MeVCM4_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_485MeVCM4_NegHel->SetName("ParaPerpAsymmPhip485MeVCM4");
  ParaPerpAsymmPhip_485MeVCM4_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 480-490MeV Photon Energy -ve Hel (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_485MeVCM4_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[3][6] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[3][6] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_495MeVCM4_NegHel = Phip_495MeVCM4_Para_NegHel->GetAsymmetry(Phip_495MeVCM4_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_495MeVCM4_NegHel->SetName("ParaPerpAsymmPhip495MeVCM4");
  ParaPerpAsymmPhip_495MeVCM4_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 490-500MeV Photon Energy -ve Hel (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_495MeVCM4_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[3][7] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[3][7] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_505MeVCM4_NegHel = Phip_505MeVCM4_Para_NegHel->GetAsymmetry(Phip_505MeVCM4_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_505MeVCM4_NegHel->SetName("ParaPerpAsymmPhip505MeVCM4");
  ParaPerpAsymmPhip_505MeVCM4_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 500-510MeV Photon Energy -ve Hel (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_505MeVCM4_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[3][8] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[3][8] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_515MeVCM4_NegHel = Phip_515MeVCM4_Para_NegHel->GetAsymmetry(Phip_515MeVCM4_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_515MeVCM4_NegHel->SetName("ParaPerpAsymmPhip515MeVCM4");
  ParaPerpAsymmPhip_515MeVCM4_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 510-520MeV Photon Energy -ve Hel (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_515MeVCM4_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[3][9] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[3][9] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_525MeVCM4_NegHel = Phip_525MeVCM4_Para_NegHel->GetAsymmetry(Phip_525MeVCM4_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_525MeVCM4_NegHel->SetName("ParaPerpAsymmPhip525MeVCM4");
  ParaPerpAsymmPhip_525MeVCM4_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 520-530MeV Photon Energy -ve Hel (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_525MeVCM4_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[3][10] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[3][10] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_535MeVCM4_NegHel = Phip_535MeVCM4_Para_NegHel->GetAsymmetry(Phip_535MeVCM4_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_535MeVCM4_NegHel->SetName("ParaPerpAsymmPhip535MeVCM4");
  ParaPerpAsymmPhip_535MeVCM4_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 530-540MeV Photon Energy -ve Hel (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_535MeVCM4_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[3][11] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[3][11] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_545MeVCM4_NegHel = Phip_545MeVCM4_Para_NegHel->GetAsymmetry(Phip_545MeVCM4_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_545MeVCM4_NegHel->SetName("ParaPerpAsymmPhip545MeVCM4");
  ParaPerpAsymmPhip_545MeVCM4_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 540-550MeV Photon Energy -ve Hel (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_545MeVCM4_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[3][12] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[3][12] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_555MeVCM4_NegHel = Phip_555MeVCM4_Para_NegHel->GetAsymmetry(Phip_555MeVCM4_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_555MeVCM4_NegHel->SetName("ParaPerpAsymmPhip555MeVCM4");
  ParaPerpAsymmPhip_555MeVCM4_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 550-560MeV Photon Energy -ve Hel (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_555MeVCM4_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[3][13] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[3][13] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_565MeVCM4_NegHel = Phip_565MeVCM4_Para_NegHel->GetAsymmetry(Phip_565MeVCM4_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_565MeVCM4_NegHel->SetName("ParaPerpAsymmPhip565MeVCM4");
  ParaPerpAsymmPhip_565MeVCM4_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 560-570MeV Photon Energy -ve Hel (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_565MeVCM4_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[3][14] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[3][14] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_575MeVCM4_NegHel = Phip_575MeVCM4_Para_NegHel->GetAsymmetry(Phip_575MeVCM4_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_575MeVCM4_NegHel->SetName("ParaPerpAsymmPhip575MeVCM4");
  ParaPerpAsymmPhip_575MeVCM4_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 570-580MeV Photon Energy -ve Hel (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_575MeVCM4_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[3][15] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[3][15] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_585MeVCM4_NegHel = Phip_585MeVCM4_Para_NegHel->GetAsymmetry(Phip_585MeVCM4_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_585MeVCM4_NegHel->SetName("ParaPerpAsymmPhip585MeVCM4");
  ParaPerpAsymmPhip_585MeVCM4_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 580-590MeV Photon Energy -ve Hel (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_585MeVCM4_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[3][16] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[3][16] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_595MeVCM4_NegHel = Phip_595MeVCM4_Para_NegHel->GetAsymmetry(Phip_595MeVCM4_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_595MeVCM4_NegHel->SetName("ParaPerpAsymmPhip595MeVCM4");
  ParaPerpAsymmPhip_595MeVCM4_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 590-600MeV Photon Energy -ve Hel (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_595MeVCM4_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[3][17] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[3][17] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_605MeVCM4_NegHel = Phip_605MeVCM4_Para_NegHel->GetAsymmetry(Phip_605MeVCM4_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_605MeVCM4_NegHel->SetName("ParaPerpAsymmPhip605MeVCM4");
  ParaPerpAsymmPhip_605MeVCM4_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 600-610MeV Photon Energy -ve Hel (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_605MeVCM4_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[3][18] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[3][18] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_615MeVCM4_NegHel = Phip_615MeVCM4_Para_NegHel->GetAsymmetry(Phip_615MeVCM4_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_615MeVCM4_NegHel->SetName("ParaPerpAsymmPhip615MeVCM4");
  ParaPerpAsymmPhip_615MeVCM4_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 610-620MeV Photon Energy -ve Hel (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_615MeVCM4_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[3][19] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[3][19] = CosFit->GetParError(0);

  ///////////////////////////////////////////
  //////////////////CM5 p ///////////////////
  ///////////////////////////////////////////

  ParaPerpAsymmPhip_425MeVCM5_NegHel = Phip_425MeVCM5_Para_NegHel->GetAsymmetry(Phip_425MeVCM5_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_425MeVCM5_NegHel->SetName("ParaPerpAsymmPhip425MeVCM5");
  ParaPerpAsymmPhip_425MeVCM5_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 420-430MeV Photon Energy -ve Hel (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_425MeVCM5_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[4][0] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[4][0] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_435MeVCM5_NegHel = Phip_435MeVCM5_Para_NegHel->GetAsymmetry(Phip_435MeVCM5_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_435MeVCM5_NegHel->SetName("ParaPerpAsymmPhip435MeVCM5");
  ParaPerpAsymmPhip_435MeVCM5_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 430-440MeV Photon Energy -ve Hel (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_435MeVCM5_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[4][1] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[4][1] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_445MeVCM5_NegHel = Phip_445MeVCM5_Para_NegHel->GetAsymmetry(Phip_445MeVCM5_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_445MeVCM5_NegHel->SetName("ParaPerpAsymmPhip445MeVCM5");
  ParaPerpAsymmPhip_445MeVCM5_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 440-450MeV Photon Energy -ve Hel (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_445MeVCM5_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[4][2] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[4][2] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_455MeVCM5_NegHel = Phip_455MeVCM5_Para_NegHel->GetAsymmetry(Phip_455MeVCM5_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_455MeVCM5_NegHel->SetName("ParaPerpAsymmPhip455MeVCM5");
  ParaPerpAsymmPhip_455MeVCM5_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 450-460MeV Photon Energy -ve Hel (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_455MeVCM5_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[4][3] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[4][3] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_465MeVCM5_NegHel = Phip_465MeVCM5_Para_NegHel->GetAsymmetry(Phip_465MeVCM5_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_465MeVCM5_NegHel->SetName("ParaPerpAsymmPhip465MeVCM5");
  ParaPerpAsymmPhip_465MeVCM5_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 460-470MeV Photon Energy -ve Hel (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_465MeVCM5_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[4][4] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[4][4] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_475MeVCM5_NegHel = Phip_475MeVCM5_Para_NegHel->GetAsymmetry(Phip_475MeVCM5_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_475MeVCM5_NegHel->SetName("ParaPerpAsymmPhip475MeVCM5");
  ParaPerpAsymmPhip_475MeVCM5_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 470-480MeV Photon Energy -ve Hel (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_475MeVCM5_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[4][5] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[4][5] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_485MeVCM5_NegHel = Phip_485MeVCM5_Para_NegHel->GetAsymmetry(Phip_485MeVCM5_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_485MeVCM5_NegHel->SetName("ParaPerpAsymmPhip485MeVCM5");
  ParaPerpAsymmPhip_485MeVCM5_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 480-490MeV Photon Energy -ve Hel (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_485MeVCM5_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[4][6] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[4][6] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_495MeVCM5_NegHel = Phip_495MeVCM5_Para_NegHel->GetAsymmetry(Phip_495MeVCM5_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_495MeVCM5_NegHel->SetName("ParaPerpAsymmPhip495MeVCM5");
  ParaPerpAsymmPhip_495MeVCM5_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 490-500MeV Photon Energy -ve Hel (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_495MeVCM5_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[4][7] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[4][7] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_505MeVCM5_NegHel = Phip_505MeVCM5_Para_NegHel->GetAsymmetry(Phip_505MeVCM5_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_505MeVCM5_NegHel->SetName("ParaPerpAsymmPhip505MeVCM5");
  ParaPerpAsymmPhip_505MeVCM5_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 500-510MeV Photon Energy -ve Hel (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_505MeVCM5_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[4][8] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[4][8] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_515MeVCM5_NegHel = Phip_515MeVCM5_Para_NegHel->GetAsymmetry(Phip_515MeVCM5_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_515MeVCM5_NegHel->SetName("ParaPerpAsymmPhip515MeVCM5");
  ParaPerpAsymmPhip_515MeVCM5_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 510-520MeV Photon Energy -ve Hel (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_515MeVCM5_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[4][9] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[4][9] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_525MeVCM5_NegHel = Phip_525MeVCM5_Para_NegHel->GetAsymmetry(Phip_525MeVCM5_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_525MeVCM5_NegHel->SetName("ParaPerpAsymmPhip525MeVCM5");
  ParaPerpAsymmPhip_525MeVCM5_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 520-530MeV Photon Energy -ve Hel (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_525MeVCM5_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[4][10] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[4][10] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_535MeVCM5_NegHel = Phip_535MeVCM5_Para_NegHel->GetAsymmetry(Phip_535MeVCM5_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_535MeVCM5_NegHel->SetName("ParaPerpAsymmPhip535MeVCM5");
  ParaPerpAsymmPhip_535MeVCM5_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 530-540MeV Photon Energy -ve Hel (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_535MeVCM5_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[4][11] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[4][11] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_545MeVCM5_NegHel = Phip_545MeVCM5_Para_NegHel->GetAsymmetry(Phip_545MeVCM5_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_545MeVCM5_NegHel->SetName("ParaPerpAsymmPhip545MeVCM5");
  ParaPerpAsymmPhip_545MeVCM5_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 540-550MeV Photon Energy -ve Hel (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_545MeVCM5_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[4][12] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[4][12] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_555MeVCM5_NegHel = Phip_555MeVCM5_Para_NegHel->GetAsymmetry(Phip_555MeVCM5_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_555MeVCM5_NegHel->SetName("ParaPerpAsymmPhip555MeVCM5");
  ParaPerpAsymmPhip_555MeVCM5_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 550-560MeV Photon Energy -ve Hel (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_555MeVCM5_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[4][13] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[4][13] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_565MeVCM5_NegHel = Phip_565MeVCM5_Para_NegHel->GetAsymmetry(Phip_565MeVCM5_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_565MeVCM5_NegHel->SetName("ParaPerpAsymmPhip565MeVCM5");
  ParaPerpAsymmPhip_565MeVCM5_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 560-570MeV Photon Energy -ve Hel (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_565MeVCM5_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[4][14] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[4][14] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_575MeVCM5_NegHel = Phip_575MeVCM5_Para_NegHel->GetAsymmetry(Phip_575MeVCM5_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_575MeVCM5_NegHel->SetName("ParaPerpAsymmPhip575MeVCM5");
  ParaPerpAsymmPhip_575MeVCM5_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 570-580MeV Photon Energy -ve Hel (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_575MeVCM5_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[4][15] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[4][15] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_585MeVCM5_NegHel = Phip_585MeVCM5_Para_NegHel->GetAsymmetry(Phip_585MeVCM5_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_585MeVCM5_NegHel->SetName("ParaPerpAsymmPhip585MeVCM5");
  ParaPerpAsymmPhip_585MeVCM5_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 580-590MeV Photon Energy -ve Hel (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_585MeVCM5_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[4][16] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[4][16] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_595MeVCM5_NegHel = Phip_595MeVCM5_Para_NegHel->GetAsymmetry(Phip_595MeVCM5_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_595MeVCM5_NegHel->SetName("ParaPerpAsymmPhip595MeVCM5");
  ParaPerpAsymmPhip_595MeVCM5_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 590-600MeV Photon Energy -ve Hel (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_595MeVCM5_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[4][17] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[4][17] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_605MeVCM5_NegHel = Phip_605MeVCM5_Para_NegHel->GetAsymmetry(Phip_605MeVCM5_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_605MeVCM5_NegHel->SetName("ParaPerpAsymmPhip605MeVCM5");
  ParaPerpAsymmPhip_605MeVCM5_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 600-610MeV Photon Energy -ve Hel (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_605MeVCM5_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[4][18] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[4][18] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_615MeVCM5_NegHel = Phip_615MeVCM5_Para_NegHel->GetAsymmetry(Phip_615MeVCM5_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_615MeVCM5_NegHel->SetName("ParaPerpAsymmPhip615MeVCM5");
  ParaPerpAsymmPhip_615MeVCM5_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 610-620MeV Photon Energy -ve Hel (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_615MeVCM5_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[4][19] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[4][19] = CosFit->GetParError(0);

  ///////////////////////////////////////////
  //////////////////CM6 p ///////////////////
  ///////////////////////////////////////////

  ParaPerpAsymmPhip_425MeVCM6_NegHel = Phip_425MeVCM6_Para_NegHel->GetAsymmetry(Phip_425MeVCM6_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_425MeVCM6_NegHel->SetName("ParaPerpAsymmPhip425MeVCM6");
  ParaPerpAsymmPhip_425MeVCM6_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 420-430MeV Photon Energy -ve Hel (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_425MeVCM6_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[5][0] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[5][0] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_435MeVCM6_NegHel = Phip_435MeVCM6_Para_NegHel->GetAsymmetry(Phip_435MeVCM6_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_435MeVCM6_NegHel->SetName("ParaPerpAsymmPhip435MeVCM6");
  ParaPerpAsymmPhip_435MeVCM6_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 430-440MeV Photon Energy -ve Hel (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_435MeVCM6_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[5][1] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[5][1] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_445MeVCM6_NegHel = Phip_445MeVCM6_Para_NegHel->GetAsymmetry(Phip_445MeVCM6_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_445MeVCM6_NegHel->SetName("ParaPerpAsymmPhip445MeVCM6");
  ParaPerpAsymmPhip_445MeVCM6_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 440-450MeV Photon Energy -ve Hel (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_445MeVCM6_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[5][2] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[5][2] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_455MeVCM6_NegHel = Phip_455MeVCM6_Para_NegHel->GetAsymmetry(Phip_455MeVCM6_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_455MeVCM6_NegHel->SetName("ParaPerpAsymmPhip455MeVCM6");
  ParaPerpAsymmPhip_455MeVCM6_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 450-460MeV Photon Energy -ve Hel (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_455MeVCM6_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[5][3] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[5][3] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_465MeVCM6_NegHel = Phip_465MeVCM6_Para_NegHel->GetAsymmetry(Phip_465MeVCM6_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_465MeVCM6_NegHel->SetName("ParaPerpAsymmPhip465MeVCM6");
  ParaPerpAsymmPhip_465MeVCM6_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 460-470MeV Photon Energy -ve Hel (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_465MeVCM6_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[5][4] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[5][4] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_475MeVCM6_NegHel = Phip_475MeVCM6_Para_NegHel->GetAsymmetry(Phip_475MeVCM6_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_475MeVCM6_NegHel->SetName("ParaPerpAsymmPhip475MeVCM6");
  ParaPerpAsymmPhip_475MeVCM6_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 470-480MeV Photon Energy -ve Hel (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_475MeVCM6_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[5][5] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[5][5] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_485MeVCM6_NegHel = Phip_485MeVCM6_Para_NegHel->GetAsymmetry(Phip_485MeVCM6_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_485MeVCM6_NegHel->SetName("ParaPerpAsymmPhip485MeVCM6");
  ParaPerpAsymmPhip_485MeVCM6_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 480-490MeV Photon Energy -ve Hel (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_485MeVCM6_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[5][6] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[5][6] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_495MeVCM6_NegHel = Phip_495MeVCM6_Para_NegHel->GetAsymmetry(Phip_495MeVCM6_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_495MeVCM6_NegHel->SetName("ParaPerpAsymmPhip495MeVCM6");
  ParaPerpAsymmPhip_495MeVCM6_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 490-500MeV Photon Energy -ve Hel (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_495MeVCM6_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[5][7] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[5][7] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_505MeVCM6_NegHel = Phip_505MeVCM6_Para_NegHel->GetAsymmetry(Phip_505MeVCM6_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_505MeVCM6_NegHel->SetName("ParaPerpAsymmPhip505MeVCM6");
  ParaPerpAsymmPhip_505MeVCM6_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 500-510MeV Photon Energy -ve Hel (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_505MeVCM6_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[5][8] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[5][8] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_515MeVCM6_NegHel = Phip_515MeVCM6_Para_NegHel->GetAsymmetry(Phip_515MeVCM6_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_515MeVCM6_NegHel->SetName("ParaPerpAsymmPhip515MeVCM6");
  ParaPerpAsymmPhip_515MeVCM6_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 510-520MeV Photon Energy -ve Hel (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_515MeVCM6_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[5][9] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[5][9] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_525MeVCM6_NegHel = Phip_525MeVCM6_Para_NegHel->GetAsymmetry(Phip_525MeVCM6_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_525MeVCM6_NegHel->SetName("ParaPerpAsymmPhip525MeVCM6");
  ParaPerpAsymmPhip_525MeVCM6_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 520-530MeV Photon Energy -ve Hel (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_525MeVCM6_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[5][10] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[5][10] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_535MeVCM6_NegHel = Phip_535MeVCM6_Para_NegHel->GetAsymmetry(Phip_535MeVCM6_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_535MeVCM6_NegHel->SetName("ParaPerpAsymmPhip535MeVCM6");
  ParaPerpAsymmPhip_535MeVCM6_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 530-540MeV Photon Energy -ve Hel (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_535MeVCM6_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[5][11] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[5][11] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_545MeVCM6_NegHel = Phip_545MeVCM6_Para_NegHel->GetAsymmetry(Phip_545MeVCM6_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_545MeVCM6_NegHel->SetName("ParaPerpAsymmPhip545MeVCM6");
  ParaPerpAsymmPhip_545MeVCM6_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 540-550MeV Photon Energy -ve Hel (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_545MeVCM6_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[5][12] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[5][12] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_555MeVCM6_NegHel = Phip_555MeVCM6_Para_NegHel->GetAsymmetry(Phip_555MeVCM6_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_555MeVCM6_NegHel->SetName("ParaPerpAsymmPhip555MeVCM6");
  ParaPerpAsymmPhip_555MeVCM6_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 550-560MeV Photon Energy -ve Hel (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_555MeVCM6_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[5][13] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[5][13] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_565MeVCM6_NegHel = Phip_565MeVCM6_Para_NegHel->GetAsymmetry(Phip_565MeVCM6_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_565MeVCM6_NegHel->SetName("ParaPerpAsymmPhip565MeVCM6");
  ParaPerpAsymmPhip_565MeVCM6_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 560-570MeV Photon Energy -ve Hel (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_565MeVCM6_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[5][14] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[5][14] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_575MeVCM6_NegHel = Phip_575MeVCM6_Para_NegHel->GetAsymmetry(Phip_575MeVCM6_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_575MeVCM6_NegHel->SetName("ParaPerpAsymmPhip575MeVCM6");
  ParaPerpAsymmPhip_575MeVCM6_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 570-580MeV Photon Energy -ve Hel (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_575MeVCM6_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[5][15] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[5][15] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_585MeVCM6_NegHel = Phip_585MeVCM6_Para_NegHel->GetAsymmetry(Phip_585MeVCM6_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_585MeVCM6_NegHel->SetName("ParaPerpAsymmPhip585MeVCM6");
  ParaPerpAsymmPhip_585MeVCM6_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 580-590MeV Photon Energy -ve Hel (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_585MeVCM6_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[5][16] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[5][16] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_595MeVCM6_NegHel = Phip_595MeVCM6_Para_NegHel->GetAsymmetry(Phip_595MeVCM6_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_595MeVCM6_NegHel->SetName("ParaPerpAsymmPhip595MeVCM6");
  ParaPerpAsymmPhip_595MeVCM6_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 590-600MeV Photon Energy -ve Hel (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_595MeVCM6_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[5][17] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[5][17] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_605MeVCM6_NegHel = Phip_605MeVCM6_Para_NegHel->GetAsymmetry(Phip_605MeVCM6_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_605MeVCM6_NegHel->SetName("ParaPerpAsymmPhip605MeVCM6");
  ParaPerpAsymmPhip_605MeVCM6_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 600-610MeV Photon Energy -ve Hel (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_605MeVCM6_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[5][18] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[5][18] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_615MeVCM6_NegHel = Phip_615MeVCM6_Para_NegHel->GetAsymmetry(Phip_615MeVCM6_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_615MeVCM6_NegHel->SetName("ParaPerpAsymmPhip615MeVCM6");
  ParaPerpAsymmPhip_615MeVCM6_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 610-620MeV Photon Energy -ve Hel (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_615MeVCM6_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[5][19] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[5][19] = CosFit->GetParError(0);

  ///////////////////////////////////////////
  //////////////////CM7 p ///////////////////
  ///////////////////////////////////////////

  ParaPerpAsymmPhip_425MeVCM7_NegHel = Phip_425MeVCM7_Para_NegHel->GetAsymmetry(Phip_425MeVCM7_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_425MeVCM7_NegHel->SetName("ParaPerpAsymmPhip425MeVCM7");
  ParaPerpAsymmPhip_425MeVCM7_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 420-430MeV Photon Energy -ve Hel (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_425MeVCM7_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[6][0] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[6][0] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_435MeVCM7_NegHel = Phip_435MeVCM7_Para_NegHel->GetAsymmetry(Phip_435MeVCM7_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_435MeVCM7_NegHel->SetName("ParaPerpAsymmPhip435MeVCM7");
  ParaPerpAsymmPhip_435MeVCM7_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 430-440MeV Photon Energy -ve Hel (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_435MeVCM7_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[6][1] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[6][1] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_445MeVCM7_NegHel = Phip_445MeVCM7_Para_NegHel->GetAsymmetry(Phip_445MeVCM7_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_445MeVCM7_NegHel->SetName("ParaPerpAsymmPhip445MeVCM7");
  ParaPerpAsymmPhip_445MeVCM7_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 440-450MeV Photon Energy -ve Hel (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_445MeVCM7_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[6][2] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[6][2] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_455MeVCM7_NegHel = Phip_455MeVCM7_Para_NegHel->GetAsymmetry(Phip_455MeVCM7_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_455MeVCM7_NegHel->SetName("ParaPerpAsymmPhip455MeVCM7");
  ParaPerpAsymmPhip_455MeVCM7_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 450-460MeV Photon Energy -ve Hel (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_455MeVCM7_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[6][3] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[6][3] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_465MeVCM7_NegHel = Phip_465MeVCM7_Para_NegHel->GetAsymmetry(Phip_465MeVCM7_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_465MeVCM7_NegHel->SetName("ParaPerpAsymmPhip465MeVCM7");
  ParaPerpAsymmPhip_465MeVCM7_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 460-470MeV Photon Energy -ve Hel (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_465MeVCM7_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[6][4] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[6][4] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_475MeVCM7_NegHel = Phip_475MeVCM7_Para_NegHel->GetAsymmetry(Phip_475MeVCM7_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_475MeVCM7_NegHel->SetName("ParaPerpAsymmPhip475MeVCM7");
  ParaPerpAsymmPhip_475MeVCM7_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 470-480MeV Photon Energy -ve Hel (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_475MeVCM7_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[6][5] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[6][5] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_485MeVCM7_NegHel = Phip_485MeVCM7_Para_NegHel->GetAsymmetry(Phip_485MeVCM7_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_485MeVCM7_NegHel->SetName("ParaPerpAsymmPhip485MeVCM7");
  ParaPerpAsymmPhip_485MeVCM7_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 480-490MeV Photon Energy -ve Hel (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_485MeVCM7_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[6][6] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[6][6] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_495MeVCM7_NegHel = Phip_495MeVCM7_Para_NegHel->GetAsymmetry(Phip_495MeVCM7_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_495MeVCM7_NegHel->SetName("ParaPerpAsymmPhip495MeVCM7");
  ParaPerpAsymmPhip_495MeVCM7_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 490-500MeV Photon Energy -ve Hel (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_495MeVCM7_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[6][7] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[6][7] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_505MeVCM7_NegHel = Phip_505MeVCM7_Para_NegHel->GetAsymmetry(Phip_505MeVCM7_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_505MeVCM7_NegHel->SetName("ParaPerpAsymmPhip505MeVCM7");
  ParaPerpAsymmPhip_505MeVCM7_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 500-510MeV Photon Energy -ve Hel (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_505MeVCM7_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[6][8] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[6][8] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_515MeVCM7_NegHel = Phip_515MeVCM7_Para_NegHel->GetAsymmetry(Phip_515MeVCM7_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_515MeVCM7_NegHel->SetName("ParaPerpAsymmPhip515MeVCM7");
  ParaPerpAsymmPhip_515MeVCM7_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 510-520MeV Photon Energy -ve Hel (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_515MeVCM7_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[6][9] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[6][9] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_525MeVCM7_NegHel = Phip_525MeVCM7_Para_NegHel->GetAsymmetry(Phip_525MeVCM7_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_525MeVCM7_NegHel->SetName("ParaPerpAsymmPhip525MeVCM7");
  ParaPerpAsymmPhip_525MeVCM7_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 520-530MeV Photon Energy -ve Hel (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_525MeVCM7_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[6][10] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[6][10] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_535MeVCM7_NegHel = Phip_535MeVCM7_Para_NegHel->GetAsymmetry(Phip_535MeVCM7_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_535MeVCM7_NegHel->SetName("ParaPerpAsymmPhip535MeVCM7");
  ParaPerpAsymmPhip_535MeVCM7_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 530-540MeV Photon Energy -ve Hel (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_535MeVCM7_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[6][11] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[6][11] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_545MeVCM7_NegHel = Phip_545MeVCM7_Para_NegHel->GetAsymmetry(Phip_545MeVCM7_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_545MeVCM7_NegHel->SetName("ParaPerpAsymmPhip545MeVCM7");
  ParaPerpAsymmPhip_545MeVCM7_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 540-550MeV Photon Energy -ve Hel (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_545MeVCM7_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[6][12] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[6][12] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_555MeVCM7_NegHel = Phip_555MeVCM7_Para_NegHel->GetAsymmetry(Phip_555MeVCM7_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_555MeVCM7_NegHel->SetName("ParaPerpAsymmPhip555MeVCM7");
  ParaPerpAsymmPhip_555MeVCM7_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 550-560MeV Photon Energy -ve Hel (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_555MeVCM7_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[6][13] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[6][13] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_565MeVCM7_NegHel = Phip_565MeVCM7_Para_NegHel->GetAsymmetry(Phip_565MeVCM7_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_565MeVCM7_NegHel->SetName("ParaPerpAsymmPhip565MeVCM7");
  ParaPerpAsymmPhip_565MeVCM7_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 560-570MeV Photon Energy -ve Hel (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_565MeVCM7_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[6][14] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[6][14] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_575MeVCM7_NegHel = Phip_575MeVCM7_Para_NegHel->GetAsymmetry(Phip_575MeVCM7_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_575MeVCM7_NegHel->SetName("ParaPerpAsymmPhip575MeVCM7");
  ParaPerpAsymmPhip_575MeVCM7_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 570-580MeV Photon Energy -ve Hel (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_575MeVCM7_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[6][15] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[6][15] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_585MeVCM7_NegHel = Phip_585MeVCM7_Para_NegHel->GetAsymmetry(Phip_585MeVCM7_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_585MeVCM7_NegHel->SetName("ParaPerpAsymmPhip585MeVCM7");
  ParaPerpAsymmPhip_585MeVCM7_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 580-590MeV Photon Energy -ve Hel (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_585MeVCM7_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[6][16] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[6][16] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_595MeVCM7_NegHel = Phip_595MeVCM7_Para_NegHel->GetAsymmetry(Phip_595MeVCM7_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_595MeVCM7_NegHel->SetName("ParaPerpAsymmPhip595MeVCM7");
  ParaPerpAsymmPhip_595MeVCM7_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 590-600MeV Photon Energy -ve Hel (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_595MeVCM7_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[6][17] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[6][17] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_605MeVCM7_NegHel = Phip_605MeVCM7_Para_NegHel->GetAsymmetry(Phip_605MeVCM7_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_605MeVCM7_NegHel->SetName("ParaPerpAsymmPhip605MeVCM7");
  ParaPerpAsymmPhip_605MeVCM7_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 600-610MeV Photon Energy -ve Hel (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_605MeVCM7_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[6][18] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[6][18] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_615MeVCM7_NegHel = Phip_615MeVCM7_Para_NegHel->GetAsymmetry(Phip_615MeVCM7_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_615MeVCM7_NegHel->SetName("ParaPerpAsymmPhip615MeVCM7");
  ParaPerpAsymmPhip_615MeVCM7_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 610-620MeV Photon Energy -ve Hel (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_615MeVCM7_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[6][19] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[6][19] = CosFit->GetParError(0);

  ///////////////////////////////////////////
  //////////////////CM8 p ///////////////////
  ///////////////////////////////////////////

  ParaPerpAsymmPhip_425MeVCM8_NegHel = Phip_425MeVCM8_Para_NegHel->GetAsymmetry(Phip_425MeVCM8_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_425MeVCM8_NegHel->SetName("ParaPerpAsymmPhip425MeVCM8");
  ParaPerpAsymmPhip_425MeVCM8_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 420-430MeV Photon Energy -ve Hel (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_425MeVCM8_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[7][0] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[7][0] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_435MeVCM8_NegHel = Phip_435MeVCM8_Para_NegHel->GetAsymmetry(Phip_435MeVCM8_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_435MeVCM8_NegHel->SetName("ParaPerpAsymmPhip435MeVCM8");
  ParaPerpAsymmPhip_435MeVCM8_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 430-440MeV Photon Energy -ve Hel (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_435MeVCM8_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[7][1] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[7][1] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_445MeVCM8_NegHel = Phip_445MeVCM8_Para_NegHel->GetAsymmetry(Phip_445MeVCM8_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_445MeVCM8_NegHel->SetName("ParaPerpAsymmPhip445MeVCM8");
  ParaPerpAsymmPhip_445MeVCM8_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 440-450MeV Photon Energy -ve Hel (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_445MeVCM8_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[7][2] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[7][2] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_455MeVCM8_NegHel = Phip_455MeVCM8_Para_NegHel->GetAsymmetry(Phip_455MeVCM8_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_455MeVCM8_NegHel->SetName("ParaPerpAsymmPhip455MeVCM8");
  ParaPerpAsymmPhip_455MeVCM8_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 450-460MeV Photon Energy -ve Hel (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_455MeVCM8_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[7][3] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[7][3] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_465MeVCM8_NegHel = Phip_465MeVCM8_Para_NegHel->GetAsymmetry(Phip_465MeVCM8_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_465MeVCM8_NegHel->SetName("ParaPerpAsymmPhip465MeVCM8");
  ParaPerpAsymmPhip_465MeVCM8_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 460-470MeV Photon Energy -ve Hel (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_465MeVCM8_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[7][4] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[7][4] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_475MeVCM8_NegHel = Phip_475MeVCM8_Para_NegHel->GetAsymmetry(Phip_475MeVCM8_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_475MeVCM8_NegHel->SetName("ParaPerpAsymmPhip475MeVCM8");
  ParaPerpAsymmPhip_475MeVCM8_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 470-480MeV Photon Energy -ve Hel (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_475MeVCM8_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[7][5] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[7][5] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_485MeVCM8_NegHel = Phip_485MeVCM8_Para_NegHel->GetAsymmetry(Phip_485MeVCM8_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_485MeVCM8_NegHel->SetName("ParaPerpAsymmPhip485MeVCM8");
  ParaPerpAsymmPhip_485MeVCM8_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 480-490MeV Photon Energy -ve Hel (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_485MeVCM8_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[7][6] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[7][6] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_495MeVCM8_NegHel = Phip_495MeVCM8_Para_NegHel->GetAsymmetry(Phip_495MeVCM8_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_495MeVCM8_NegHel->SetName("ParaPerpAsymmPhip495MeVCM8");
  ParaPerpAsymmPhip_495MeVCM8_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 490-500MeV Photon Energy -ve Hel (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_495MeVCM8_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[7][7] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[7][7] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_505MeVCM8_NegHel = Phip_505MeVCM8_Para_NegHel->GetAsymmetry(Phip_505MeVCM8_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_505MeVCM8_NegHel->SetName("ParaPerpAsymmPhip505MeVCM8");
  ParaPerpAsymmPhip_505MeVCM8_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 500-510MeV Photon Energy -ve Hel (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_505MeVCM8_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[7][8] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[7][8] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_515MeVCM8_NegHel = Phip_515MeVCM8_Para_NegHel->GetAsymmetry(Phip_515MeVCM8_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_515MeVCM8_NegHel->SetName("ParaPerpAsymmPhip515MeVCM8");
  ParaPerpAsymmPhip_515MeVCM8_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 510-520MeV Photon Energy -ve Hel (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_515MeVCM8_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[7][9] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[7][9] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_525MeVCM8_NegHel = Phip_525MeVCM8_Para_NegHel->GetAsymmetry(Phip_525MeVCM8_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_525MeVCM8_NegHel->SetName("ParaPerpAsymmPhip525MeVCM8");
  ParaPerpAsymmPhip_525MeVCM8_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 520-530MeV Photon Energy -ve Hel (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_525MeVCM8_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[7][10] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[7][10] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_535MeVCM8_NegHel = Phip_535MeVCM8_Para_NegHel->GetAsymmetry(Phip_535MeVCM8_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_535MeVCM8_NegHel->SetName("ParaPerpAsymmPhip535MeVCM8");
  ParaPerpAsymmPhip_535MeVCM8_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 530-540MeV Photon Energy -ve Hel (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_535MeVCM8_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[7][11] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[7][11] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_545MeVCM8_NegHel = Phip_545MeVCM8_Para_NegHel->GetAsymmetry(Phip_545MeVCM8_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_545MeVCM8_NegHel->SetName("ParaPerpAsymmPhip545MeVCM8");
  ParaPerpAsymmPhip_545MeVCM8_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 540-550MeV Photon Energy -ve Hel (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_545MeVCM8_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[7][12] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[7][12] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_555MeVCM8_NegHel = Phip_555MeVCM8_Para_NegHel->GetAsymmetry(Phip_555MeVCM8_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_555MeVCM8_NegHel->SetName("ParaPerpAsymmPhip555MeVCM8");
  ParaPerpAsymmPhip_555MeVCM8_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 550-560MeV Photon Energy -ve Hel (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_555MeVCM8_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[7][13] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[7][13] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_565MeVCM8_NegHel = Phip_565MeVCM8_Para_NegHel->GetAsymmetry(Phip_565MeVCM8_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_565MeVCM8_NegHel->SetName("ParaPerpAsymmPhip565MeVCM8");
  ParaPerpAsymmPhip_565MeVCM8_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 560-570MeV Photon Energy -ve Hel (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_565MeVCM8_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[7][14] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[7][14] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_575MeVCM8_NegHel = Phip_575MeVCM8_Para_NegHel->GetAsymmetry(Phip_575MeVCM8_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_575MeVCM8_NegHel->SetName("ParaPerpAsymmPhip575MeVCM8");
  ParaPerpAsymmPhip_575MeVCM8_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 570-580MeV Photon Energy -ve Hel (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_575MeVCM8_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[7][15] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[7][15] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_585MeVCM8_NegHel = Phip_585MeVCM8_Para_NegHel->GetAsymmetry(Phip_585MeVCM8_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_585MeVCM8_NegHel->SetName("ParaPerpAsymmPhip585MeVCM8");
  ParaPerpAsymmPhip_585MeVCM8_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 580-590MeV Photon Energy -ve Hel (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_585MeVCM8_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[7][16] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[7][16] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_595MeVCM8_NegHel = Phip_595MeVCM8_Para_NegHel->GetAsymmetry(Phip_595MeVCM8_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_595MeVCM8_NegHel->SetName("ParaPerpAsymmPhip595MeVCM8");
  ParaPerpAsymmPhip_595MeVCM8_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 590-600MeV Photon Energy -ve Hel (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_595MeVCM8_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[7][17] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[7][17] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_605MeVCM8_NegHel = Phip_605MeVCM8_Para_NegHel->GetAsymmetry(Phip_605MeVCM8_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_605MeVCM8_NegHel->SetName("ParaPerpAsymmPhip605MeVCM8");
  ParaPerpAsymmPhip_605MeVCM8_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 600-610MeV Photon Energy -ve Hel (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_605MeVCM8_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[7][18] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[7][18] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_615MeVCM8_NegHel = Phip_615MeVCM8_Para_NegHel->GetAsymmetry(Phip_615MeVCM8_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_615MeVCM8_NegHel->SetName("ParaPerpAsymmPhip615MeVCM8");
  ParaPerpAsymmPhip_615MeVCM8_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 610-620MeV Photon Energy -ve Hel (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_615MeVCM8_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[7][19] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[7][19] = CosFit->GetParError(0);

  ///////////////////////////////////////////
  //////////////////CM9 p ///////////////////
  ///////////////////////////////////////////

  ParaPerpAsymmPhip_425MeVCM9_NegHel = Phip_425MeVCM9_Para_NegHel->GetAsymmetry(Phip_425MeVCM9_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_425MeVCM9_NegHel->SetName("ParaPerpAsymmPhip425MeVCM9");
  ParaPerpAsymmPhip_425MeVCM9_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 420-430MeV Photon Energy -ve Hel (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_425MeVCM9_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[8][0] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[8][0] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_435MeVCM9_NegHel = Phip_435MeVCM9_Para_NegHel->GetAsymmetry(Phip_435MeVCM9_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_435MeVCM9_NegHel->SetName("ParaPerpAsymmPhip435MeVCM9");
  ParaPerpAsymmPhip_435MeVCM9_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 430-440MeV Photon Energy -ve Hel (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_435MeVCM9_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[8][1] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[8][1] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_445MeVCM9_NegHel = Phip_445MeVCM9_Para_NegHel->GetAsymmetry(Phip_445MeVCM9_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_445MeVCM9_NegHel->SetName("ParaPerpAsymmPhip445MeVCM9");
  ParaPerpAsymmPhip_445MeVCM9_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 440-450MeV Photon Energy -ve Hel (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_445MeVCM9_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[8][2] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[8][2] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_455MeVCM9_NegHel = Phip_455MeVCM9_Para_NegHel->GetAsymmetry(Phip_455MeVCM9_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_455MeVCM9_NegHel->SetName("ParaPerpAsymmPhip455MeVCM9");
  ParaPerpAsymmPhip_455MeVCM9_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 450-460MeV Photon Energy -ve Hel (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_455MeVCM9_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[8][3] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[8][3] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_465MeVCM9_NegHel = Phip_465MeVCM9_Para_NegHel->GetAsymmetry(Phip_465MeVCM9_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_465MeVCM9_NegHel->SetName("ParaPerpAsymmPhip465MeVCM9");
  ParaPerpAsymmPhip_465MeVCM9_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 460-470MeV Photon Energy -ve Hel (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_465MeVCM9_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[8][4] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[8][4] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_475MeVCM9_NegHel = Phip_475MeVCM9_Para_NegHel->GetAsymmetry(Phip_475MeVCM9_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_475MeVCM9_NegHel->SetName("ParaPerpAsymmPhip475MeVCM9");
  ParaPerpAsymmPhip_475MeVCM9_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 470-480MeV Photon Energy -ve Hel (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_475MeVCM9_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[8][5] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[8][5] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_485MeVCM9_NegHel = Phip_485MeVCM9_Para_NegHel->GetAsymmetry(Phip_485MeVCM9_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_485MeVCM9_NegHel->SetName("ParaPerpAsymmPhip485MeVCM9");
  ParaPerpAsymmPhip_485MeVCM9_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 480-490MeV Photon Energy -ve Hel (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_485MeVCM9_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[8][6] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[8][6] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_495MeVCM9_NegHel = Phip_495MeVCM9_Para_NegHel->GetAsymmetry(Phip_495MeVCM9_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_495MeVCM9_NegHel->SetName("ParaPerpAsymmPhip495MeVCM9");
  ParaPerpAsymmPhip_495MeVCM9_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 490-500MeV Photon Energy -ve Hel (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_495MeVCM9_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[8][7] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[8][7] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_505MeVCM9_NegHel = Phip_505MeVCM9_Para_NegHel->GetAsymmetry(Phip_505MeVCM9_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_505MeVCM9_NegHel->SetName("ParaPerpAsymmPhip505MeVCM9");
  ParaPerpAsymmPhip_505MeVCM9_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 500-510MeV Photon Energy -ve Hel (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_505MeVCM9_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[8][8] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[8][8] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_515MeVCM9_NegHel = Phip_515MeVCM9_Para_NegHel->GetAsymmetry(Phip_515MeVCM9_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_515MeVCM9_NegHel->SetName("ParaPerpAsymmPhip515MeVCM9");
  ParaPerpAsymmPhip_515MeVCM9_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 510-520MeV Photon Energy -ve Hel (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_515MeVCM9_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[8][9] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[8][9] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_525MeVCM9_NegHel = Phip_525MeVCM9_Para_NegHel->GetAsymmetry(Phip_525MeVCM9_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_525MeVCM9_NegHel->SetName("ParaPerpAsymmPhip525MeVCM9");
  ParaPerpAsymmPhip_525MeVCM9_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 520-530MeV Photon Energy -ve Hel (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_525MeVCM9_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[8][10] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[8][10] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_535MeVCM9_NegHel = Phip_535MeVCM9_Para_NegHel->GetAsymmetry(Phip_535MeVCM9_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_535MeVCM9_NegHel->SetName("ParaPerpAsymmPhip535MeVCM9");
  ParaPerpAsymmPhip_535MeVCM9_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 530-540MeV Photon Energy -ve Hel (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_535MeVCM9_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[8][11] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[8][11] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_545MeVCM9_NegHel = Phip_545MeVCM9_Para_NegHel->GetAsymmetry(Phip_545MeVCM9_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_545MeVCM9_NegHel->SetName("ParaPerpAsymmPhip545MeVCM9");
  ParaPerpAsymmPhip_545MeVCM9_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 540-550MeV Photon Energy -ve Hel (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_545MeVCM9_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[8][12] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[8][12] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_555MeVCM9_NegHel = Phip_555MeVCM9_Para_NegHel->GetAsymmetry(Phip_555MeVCM9_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_555MeVCM9_NegHel->SetName("ParaPerpAsymmPhip555MeVCM9");
  ParaPerpAsymmPhip_555MeVCM9_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 550-560MeV Photon Energy -ve Hel (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_555MeVCM9_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[8][13] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[8][13] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_565MeVCM9_NegHel = Phip_565MeVCM9_Para_NegHel->GetAsymmetry(Phip_565MeVCM9_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_565MeVCM9_NegHel->SetName("ParaPerpAsymmPhip565MeVCM9");
  ParaPerpAsymmPhip_565MeVCM9_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 560-570MeV Photon Energy -ve Hel (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_565MeVCM9_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[8][14] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[8][14] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_575MeVCM9_NegHel = Phip_575MeVCM9_Para_NegHel->GetAsymmetry(Phip_575MeVCM9_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_575MeVCM9_NegHel->SetName("ParaPerpAsymmPhip575MeVCM9");
  ParaPerpAsymmPhip_575MeVCM9_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 570-580MeV Photon Energy -ve Hel (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_575MeVCM9_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[8][15] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[8][15] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_585MeVCM9_NegHel = Phip_585MeVCM9_Para_NegHel->GetAsymmetry(Phip_585MeVCM9_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_585MeVCM9_NegHel->SetName("ParaPerpAsymmPhip585MeVCM9");
  ParaPerpAsymmPhip_585MeVCM9_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 580-590MeV Photon Energy -ve Hel (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_585MeVCM9_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[8][16] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[8][16] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_595MeVCM9_NegHel = Phip_595MeVCM9_Para_NegHel->GetAsymmetry(Phip_595MeVCM9_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_595MeVCM9_NegHel->SetName("ParaPerpAsymmPhip595MeVCM9");
  ParaPerpAsymmPhip_595MeVCM9_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 590-600MeV Photon Energy -ve Hel (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_595MeVCM9_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[8][17] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[8][17] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_605MeVCM9_NegHel = Phip_605MeVCM9_Para_NegHel->GetAsymmetry(Phip_605MeVCM9_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_605MeVCM9_NegHel->SetName("ParaPerpAsymmPhip605MeVCM9");
  ParaPerpAsymmPhip_605MeVCM9_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 600-610MeV Photon Energy -ve Hel (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_605MeVCM9_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[8][18] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[8][18] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_615MeVCM9_NegHel = Phip_615MeVCM9_Para_NegHel->GetAsymmetry(Phip_615MeVCM9_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_615MeVCM9_NegHel->SetName("ParaPerpAsymmPhip615MeVCM9");
  ParaPerpAsymmPhip_615MeVCM9_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 610-620MeV Photon Energy -ve Hel (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_615MeVCM9_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[8][19] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[8][19] = CosFit->GetParError(0);

  ///////////////////////////////////////////
  //////////////////CM10 p ///////////////////
  ///////////////////////////////////////////

  ParaPerpAsymmPhip_425MeVCM10_NegHel = Phip_425MeVCM10_Para_NegHel->GetAsymmetry(Phip_425MeVCM10_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_425MeVCM10_NegHel->SetName("ParaPerpAsymmPhip425MeVCM10");
  ParaPerpAsymmPhip_425MeVCM10_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 420-430MeV Photon Energy -ve Hel (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_425MeVCM10_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[9][0] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[9][0] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_435MeVCM10_NegHel = Phip_435MeVCM10_Para_NegHel->GetAsymmetry(Phip_435MeVCM10_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_435MeVCM10_NegHel->SetName("ParaPerpAsymmPhip435MeVCM10");
  ParaPerpAsymmPhip_435MeVCM10_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 430-440MeV Photon Energy -ve Hel (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_435MeVCM10_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[9][1] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[9][1] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_445MeVCM10_NegHel = Phip_445MeVCM10_Para_NegHel->GetAsymmetry(Phip_445MeVCM10_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_445MeVCM10_NegHel->SetName("ParaPerpAsymmPhip445MeVCM10");
  ParaPerpAsymmPhip_445MeVCM10_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 440-450MeV Photon Energy -ve Hel (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_445MeVCM10_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[9][2] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[9][2] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_455MeVCM10_NegHel = Phip_455MeVCM10_Para_NegHel->GetAsymmetry(Phip_455MeVCM10_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_455MeVCM10_NegHel->SetName("ParaPerpAsymmPhip455MeVCM10");
  ParaPerpAsymmPhip_455MeVCM10_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 450-460MeV Photon Energy -ve Hel (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_455MeVCM10_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[9][3] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[9][3] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_465MeVCM10_NegHel = Phip_465MeVCM10_Para_NegHel->GetAsymmetry(Phip_465MeVCM10_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_465MeVCM10_NegHel->SetName("ParaPerpAsymmPhip465MeVCM10");
  ParaPerpAsymmPhip_465MeVCM10_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 460-470MeV Photon Energy -ve Hel (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_465MeVCM10_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[9][4] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[9][4] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_475MeVCM10_NegHel = Phip_475MeVCM10_Para_NegHel->GetAsymmetry(Phip_475MeVCM10_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_475MeVCM10_NegHel->SetName("ParaPerpAsymmPhip475MeVCM10");
  ParaPerpAsymmPhip_475MeVCM10_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 470-480MeV Photon Energy -ve Hel (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_475MeVCM10_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[9][5] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[9][5] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_485MeVCM10_NegHel = Phip_485MeVCM10_Para_NegHel->GetAsymmetry(Phip_485MeVCM10_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_485MeVCM10_NegHel->SetName("ParaPerpAsymmPhip485MeVCM10");
  ParaPerpAsymmPhip_485MeVCM10_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 480-490MeV Photon Energy -ve Hel (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_485MeVCM10_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[9][6] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[9][6] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_495MeVCM10_NegHel = Phip_495MeVCM10_Para_NegHel->GetAsymmetry(Phip_495MeVCM10_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_495MeVCM10_NegHel->SetName("ParaPerpAsymmPhip495MeVCM10");
  ParaPerpAsymmPhip_495MeVCM10_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 490-500MeV Photon Energy -ve Hel (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_495MeVCM10_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[9][7] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[9][7] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_505MeVCM10_NegHel = Phip_505MeVCM10_Para_NegHel->GetAsymmetry(Phip_505MeVCM10_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_505MeVCM10_NegHel->SetName("ParaPerpAsymmPhip505MeVCM10");
  ParaPerpAsymmPhip_505MeVCM10_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 500-510MeV Photon Energy -ve Hel (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_505MeVCM10_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[9][8] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[9][8] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_515MeVCM10_NegHel = Phip_515MeVCM10_Para_NegHel->GetAsymmetry(Phip_515MeVCM10_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_515MeVCM10_NegHel->SetName("ParaPerpAsymmPhip515MeVCM10");
  ParaPerpAsymmPhip_515MeVCM10_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 510-520MeV Photon Energy -ve Hel (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_515MeVCM10_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[9][9] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[9][9] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_525MeVCM10_NegHel = Phip_525MeVCM10_Para_NegHel->GetAsymmetry(Phip_525MeVCM10_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_525MeVCM10_NegHel->SetName("ParaPerpAsymmPhip525MeVCM10");
  ParaPerpAsymmPhip_525MeVCM10_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 520-530MeV Photon Energy -ve Hel (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_525MeVCM10_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[9][10] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[9][10] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_535MeVCM10_NegHel = Phip_535MeVCM10_Para_NegHel->GetAsymmetry(Phip_535MeVCM10_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_535MeVCM10_NegHel->SetName("ParaPerpAsymmPhip535MeVCM10");
  ParaPerpAsymmPhip_535MeVCM10_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 530-540MeV Photon Energy -ve Hel (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_535MeVCM10_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[9][11] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[9][11] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_545MeVCM10_NegHel = Phip_545MeVCM10_Para_NegHel->GetAsymmetry(Phip_545MeVCM10_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_545MeVCM10_NegHel->SetName("ParaPerpAsymmPhip545MeVCM10");
  ParaPerpAsymmPhip_545MeVCM10_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 540-550MeV Photon Energy -ve Hel (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_545MeVCM10_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[9][12] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[9][12] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_555MeVCM10_NegHel = Phip_555MeVCM10_Para_NegHel->GetAsymmetry(Phip_555MeVCM10_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_555MeVCM10_NegHel->SetName("ParaPerpAsymmPhip555MeVCM10");
  ParaPerpAsymmPhip_555MeVCM10_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 550-560MeV Photon Energy -ve Hel (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_555MeVCM10_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[9][13] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[9][13] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_565MeVCM10_NegHel = Phip_565MeVCM10_Para_NegHel->GetAsymmetry(Phip_565MeVCM10_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_565MeVCM10_NegHel->SetName("ParaPerpAsymmPhip565MeVCM10");
  ParaPerpAsymmPhip_565MeVCM10_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 560-570MeV Photon Energy -ve Hel (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_565MeVCM10_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[9][14] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[9][14] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_575MeVCM10_NegHel = Phip_575MeVCM10_Para_NegHel->GetAsymmetry(Phip_575MeVCM10_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_575MeVCM10_NegHel->SetName("ParaPerpAsymmPhip575MeVCM10");
  ParaPerpAsymmPhip_575MeVCM10_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 570-580MeV Photon Energy -ve Hel (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_575MeVCM10_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[9][15] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[9][15] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_585MeVCM10_NegHel = Phip_585MeVCM10_Para_NegHel->GetAsymmetry(Phip_585MeVCM10_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_585MeVCM10_NegHel->SetName("ParaPerpAsymmPhip585MeVCM10");
  ParaPerpAsymmPhip_585MeVCM10_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 580-590MeV Photon Energy -ve Hel (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_585MeVCM10_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[9][16] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[9][16] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_595MeVCM10_NegHel = Phip_595MeVCM10_Para_NegHel->GetAsymmetry(Phip_595MeVCM10_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_595MeVCM10_NegHel->SetName("ParaPerpAsymmPhip595MeVCM10");
  ParaPerpAsymmPhip_595MeVCM10_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 590-600MeV Photon Energy -ve Hel (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_595MeVCM10_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[9][17] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[9][17] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_605MeVCM10_NegHel = Phip_605MeVCM10_Para_NegHel->GetAsymmetry(Phip_605MeVCM10_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_605MeVCM10_NegHel->SetName("ParaPerpAsymmPhip605MeVCM10");
  ParaPerpAsymmPhip_605MeVCM10_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 600-610MeV Photon Energy -ve Hel (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_605MeVCM10_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[9][18] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[9][18] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_615MeVCM10_NegHel = Phip_615MeVCM10_Para_NegHel->GetAsymmetry(Phip_615MeVCM10_Perp_NegHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_615MeVCM10_NegHel->SetName("ParaPerpAsymmPhip615MeVCM10");
  ParaPerpAsymmPhip_615MeVCM10_NegHel->SetTitle("Proton Para/Perp Phi Asymmetry for 610-620MeV Photon Energy -ve Hel (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_615MeVCM10_NegHel->Fit("CosFit", "Q");
  pCosAmpNegHel[9][19] = CosFit->GetParameter(0);
  pCosAmpErrNegHel[9][19] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_425MeVCM1_PosHel = Phip_425MeVCM1_Para_PosHel->GetAsymmetry(Phip_425MeVCM1_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_425MeVCM1_PosHel->SetName("ParaPerpAsymmPhip425MeVCM1");
  ParaPerpAsymmPhip_425MeVCM1_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 420-430MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_425MeVCM1_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[0][0] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[0][0] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_435MeVCM1_PosHel = Phip_435MeVCM1_Para_PosHel->GetAsymmetry(Phip_435MeVCM1_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_435MeVCM1_PosHel->SetName("ParaPerpAsymmPhip435MeVCM1");
  ParaPerpAsymmPhip_435MeVCM1_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 430-440MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_435MeVCM1_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[0][1] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[0][1] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_445MeVCM1_PosHel = Phip_445MeVCM1_Para_PosHel->GetAsymmetry(Phip_445MeVCM1_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_445MeVCM1_PosHel->SetName("ParaPerpAsymmPhip445MeVCM1");
  ParaPerpAsymmPhip_445MeVCM1_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 440-450MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_445MeVCM1_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[0][2] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[0][2] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_455MeVCM1_PosHel = Phip_455MeVCM1_Para_PosHel->GetAsymmetry(Phip_455MeVCM1_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_455MeVCM1_PosHel->SetName("ParaPerpAsymmPhip455MeVCM1");
  ParaPerpAsymmPhip_455MeVCM1_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 450-460MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_455MeVCM1_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[0][3] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[0][3] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_465MeVCM1_PosHel = Phip_465MeVCM1_Para_PosHel->GetAsymmetry(Phip_465MeVCM1_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_465MeVCM1_PosHel->SetName("ParaPerpAsymmPhip465MeVCM1");
  ParaPerpAsymmPhip_465MeVCM1_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 460-470MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_465MeVCM1_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[0][4] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[0][4] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_475MeVCM1_PosHel = Phip_475MeVCM1_Para_PosHel->GetAsymmetry(Phip_475MeVCM1_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_475MeVCM1_PosHel->SetName("ParaPerpAsymmPhip475MeVCM1");
  ParaPerpAsymmPhip_475MeVCM1_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 470-480MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_475MeVCM1_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[0][5] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[0][5] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_485MeVCM1_PosHel = Phip_485MeVCM1_Para_PosHel->GetAsymmetry(Phip_485MeVCM1_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_485MeVCM1_PosHel->SetName("ParaPerpAsymmPhip485MeVCM1");
  ParaPerpAsymmPhip_485MeVCM1_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 480-490MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_485MeVCM1_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[0][6] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[0][6] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_495MeVCM1_PosHel = Phip_495MeVCM1_Para_PosHel->GetAsymmetry(Phip_495MeVCM1_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_495MeVCM1_PosHel->SetName("ParaPerpAsymmPhip495MeVCM1");
  ParaPerpAsymmPhip_495MeVCM1_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 490-500MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_495MeVCM1_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[0][7] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[0][7] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_505MeVCM1_PosHel = Phip_505MeVCM1_Para_PosHel->GetAsymmetry(Phip_505MeVCM1_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_505MeVCM1_PosHel->SetName("ParaPerpAsymmPhip505MeVCM1");
  ParaPerpAsymmPhip_505MeVCM1_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 500-510MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_505MeVCM1_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[0][8] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[0][8] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_515MeVCM1_PosHel = Phip_515MeVCM1_Para_PosHel->GetAsymmetry(Phip_515MeVCM1_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_515MeVCM1_PosHel->SetName("ParaPerpAsymmPhip515MeVCM1");
  ParaPerpAsymmPhip_515MeVCM1_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 510-520MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_515MeVCM1_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[0][9] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[0][9] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_525MeVCM1_PosHel = Phip_525MeVCM1_Para_PosHel->GetAsymmetry(Phip_525MeVCM1_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_525MeVCM1_PosHel->SetName("ParaPerpAsymmPhip525MeVCM1");
  ParaPerpAsymmPhip_525MeVCM1_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 520-530MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_525MeVCM1_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[0][10] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[0][10] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_535MeVCM1_PosHel = Phip_535MeVCM1_Para_PosHel->GetAsymmetry(Phip_535MeVCM1_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_535MeVCM1_PosHel->SetName("ParaPerpAsymmPhip535MeVCM1");
  ParaPerpAsymmPhip_535MeVCM1_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 530-540MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_535MeVCM1_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[0][11] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[0][11] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_545MeVCM1_PosHel = Phip_545MeVCM1_Para_PosHel->GetAsymmetry(Phip_545MeVCM1_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_545MeVCM1_PosHel->SetName("ParaPerpAsymmPhip545MeVCM1");
  ParaPerpAsymmPhip_545MeVCM1_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 540-550MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_545MeVCM1_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[0][12] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[0][12] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_555MeVCM1_PosHel = Phip_555MeVCM1_Para_PosHel->GetAsymmetry(Phip_555MeVCM1_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_555MeVCM1_PosHel->SetName("ParaPerpAsymmPhip555MeVCM1");
  ParaPerpAsymmPhip_555MeVCM1_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 550-560MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_555MeVCM1_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[0][13] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[0][13] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_565MeVCM1_PosHel = Phip_565MeVCM1_Para_PosHel->GetAsymmetry(Phip_565MeVCM1_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_565MeVCM1_PosHel->SetName("ParaPerpAsymmPhip565MeVCM1");
  ParaPerpAsymmPhip_565MeVCM1_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 560-570MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_565MeVCM1_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[0][14] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[0][14] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_575MeVCM1_PosHel = Phip_575MeVCM1_Para_PosHel->GetAsymmetry(Phip_575MeVCM1_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_575MeVCM1_PosHel->SetName("ParaPerpAsymmPhip575MeVCM1");
  ParaPerpAsymmPhip_575MeVCM1_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 570-580MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_575MeVCM1_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[0][15] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[0][15] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_585MeVCM1_PosHel = Phip_585MeVCM1_Para_PosHel->GetAsymmetry(Phip_585MeVCM1_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_585MeVCM1_PosHel->SetName("ParaPerpAsymmPhip585MeVCM1");
  ParaPerpAsymmPhip_585MeVCM1_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 580-590MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_585MeVCM1_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[0][16] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[0][16] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_595MeVCM1_PosHel = Phip_595MeVCM1_Para_PosHel->GetAsymmetry(Phip_595MeVCM1_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_595MeVCM1_PosHel->SetName("ParaPerpAsymmPhip595MeVCM1");
  ParaPerpAsymmPhip_595MeVCM1_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 590-600MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_595MeVCM1_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[0][17] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[0][17] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_605MeVCM1_PosHel = Phip_605MeVCM1_Para_PosHel->GetAsymmetry(Phip_605MeVCM1_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_605MeVCM1_PosHel->SetName("ParaPerpAsymmPhip605MeVCM1");
  ParaPerpAsymmPhip_605MeVCM1_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 600-610MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_605MeVCM1_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[0][18] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[0][18] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_615MeVCM1_PosHel = Phip_615MeVCM1_Para_PosHel->GetAsymmetry(Phip_615MeVCM1_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_615MeVCM1_PosHel->SetName("ParaPerpAsymmPhip615MeVCM1");
  ParaPerpAsymmPhip_615MeVCM1_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 610-620MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_615MeVCM1_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[0][19] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[0][19] = CosFit->GetParError(0);

  ///////////////////////////////////////////
  //////////////////CM2 p ///////////////////
  ///////////////////////////////////////////

  ParaPerpAsymmPhip_425MeVCM2_PosHel = Phip_425MeVCM2_Para_PosHel->GetAsymmetry(Phip_425MeVCM2_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_425MeVCM2_PosHel->SetName("ParaPerpAsymmPhip425MeVCM2");
  ParaPerpAsymmPhip_425MeVCM2_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 420-430MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_425MeVCM2_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[1][0] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[1][0] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_435MeVCM2_PosHel = Phip_435MeVCM2_Para_PosHel->GetAsymmetry(Phip_435MeVCM2_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_435MeVCM2_PosHel->SetName("ParaPerpAsymmPhip435MeVCM2");
  ParaPerpAsymmPhip_435MeVCM2_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 430-440MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_435MeVCM2_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[1][1] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[1][1] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_445MeVCM2_PosHel = Phip_445MeVCM2_Para_PosHel->GetAsymmetry(Phip_445MeVCM2_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_445MeVCM2_PosHel->SetName("ParaPerpAsymmPhip445MeVCM2");
  ParaPerpAsymmPhip_445MeVCM2_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 440-450MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_445MeVCM2_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[1][2] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[1][2] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_455MeVCM2_PosHel = Phip_455MeVCM2_Para_PosHel->GetAsymmetry(Phip_455MeVCM2_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_455MeVCM2_PosHel->SetName("ParaPerpAsymmPhip455MeVCM2");
  ParaPerpAsymmPhip_455MeVCM2_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 450-460MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_455MeVCM2_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[1][3] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[1][3] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_465MeVCM2_PosHel = Phip_465MeVCM2_Para_PosHel->GetAsymmetry(Phip_465MeVCM2_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_465MeVCM2_PosHel->SetName("ParaPerpAsymmPhip465MeVCM2");
  ParaPerpAsymmPhip_465MeVCM2_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 460-470MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_465MeVCM2_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[1][4] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[1][4] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_475MeVCM2_PosHel = Phip_475MeVCM2_Para_PosHel->GetAsymmetry(Phip_475MeVCM2_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_475MeVCM2_PosHel->SetName("ParaPerpAsymmPhip475MeVCM2");
  ParaPerpAsymmPhip_475MeVCM2_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 470-480MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_475MeVCM2_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[1][5] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[1][5] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_485MeVCM2_PosHel = Phip_485MeVCM2_Para_PosHel->GetAsymmetry(Phip_485MeVCM2_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_485MeVCM2_PosHel->SetName("ParaPerpAsymmPhip485MeVCM2");
  ParaPerpAsymmPhip_485MeVCM2_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 480-490MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_485MeVCM2_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[1][6] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[1][6] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_495MeVCM2_PosHel = Phip_495MeVCM2_Para_PosHel->GetAsymmetry(Phip_495MeVCM2_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_495MeVCM2_PosHel->SetName("ParaPerpAsymmPhip495MeVCM2");
  ParaPerpAsymmPhip_495MeVCM2_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 490-500MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_495MeVCM2_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[1][7] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[1][7] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_505MeVCM2_PosHel = Phip_505MeVCM2_Para_PosHel->GetAsymmetry(Phip_505MeVCM2_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_505MeVCM2_PosHel->SetName("ParaPerpAsymmPhip505MeVCM2");
  ParaPerpAsymmPhip_505MeVCM2_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 500-510MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_505MeVCM2_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[1][8] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[1][8] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_515MeVCM2_PosHel = Phip_515MeVCM2_Para_PosHel->GetAsymmetry(Phip_515MeVCM2_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_515MeVCM2_PosHel->SetName("ParaPerpAsymmPhip515MeVCM2");
  ParaPerpAsymmPhip_515MeVCM2_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 510-520MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_515MeVCM2_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[1][9] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[1][9] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_525MeVCM2_PosHel = Phip_525MeVCM2_Para_PosHel->GetAsymmetry(Phip_525MeVCM2_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_525MeVCM2_PosHel->SetName("ParaPerpAsymmPhip525MeVCM2");
  ParaPerpAsymmPhip_525MeVCM2_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 520-530MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_525MeVCM2_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[1][10] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[1][10] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_535MeVCM2_PosHel = Phip_535MeVCM2_Para_PosHel->GetAsymmetry(Phip_535MeVCM2_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_535MeVCM2_PosHel->SetName("ParaPerpAsymmPhip535MeVCM2");
  ParaPerpAsymmPhip_535MeVCM2_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 530-540MeV Photon Energy -ve Hel (CosTheta1-0.8)");
  ParaPerpAsymmPhip_535MeVCM2_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[1][11] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[1][11] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_545MeVCM2_PosHel = Phip_545MeVCM2_Para_PosHel->GetAsymmetry(Phip_545MeVCM2_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_545MeVCM2_PosHel->SetName("ParaPerpAsymmPhip545MeVCM2");
  ParaPerpAsymmPhip_545MeVCM2_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 540-550MeV Photon Energy -ve Hel (CosTheta0.8-0.6)");
  ParaPerpAsymmPhip_545MeVCM2_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[1][12] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[1][12] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_555MeVCM2_PosHel = Phip_555MeVCM2_Para_PosHel->GetAsymmetry(Phip_555MeVCM2_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_555MeVCM2_PosHel->SetName("ParaPerpAsymmPhip555MeVCM2");
  ParaPerpAsymmPhip_555MeVCM2_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 550-560MeV Photon Energy -ve Hel (CosTheta0.8-0.6)");
  ParaPerpAsymmPhip_555MeVCM2_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[1][13] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[1][13] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_565MeVCM2_PosHel = Phip_565MeVCM2_Para_PosHel->GetAsymmetry(Phip_565MeVCM2_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_565MeVCM2_PosHel->SetName("ParaPerpAsymmPhip565MeVCM2");
  ParaPerpAsymmPhip_565MeVCM2_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 560-570MeV Photon Energy -ve Hel (CosTheta0.8-0.6)");
  ParaPerpAsymmPhip_565MeVCM2_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[1][14] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[1][14] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_575MeVCM2_PosHel = Phip_575MeVCM2_Para_PosHel->GetAsymmetry(Phip_575MeVCM2_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_575MeVCM2_PosHel->SetName("ParaPerpAsymmPhip575MeVCM2");
  ParaPerpAsymmPhip_575MeVCM2_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 570-580MeV Photon Energy -ve Hel (CosTheta0.8-0.6)");
  ParaPerpAsymmPhip_575MeVCM2_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[1][15] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[1][15] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_585MeVCM2_PosHel = Phip_585MeVCM2_Para_PosHel->GetAsymmetry(Phip_585MeVCM2_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_585MeVCM2_PosHel->SetName("ParaPerpAsymmPhip585MeVCM2");
  ParaPerpAsymmPhip_585MeVCM2_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 580-590MeV Photon Energy -ve Hel (CosTheta0.8-0.6)");
  ParaPerpAsymmPhip_585MeVCM2_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[1][16] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[1][16] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_595MeVCM2_PosHel = Phip_595MeVCM2_Para_PosHel->GetAsymmetry(Phip_595MeVCM2_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_595MeVCM2_PosHel->SetName("ParaPerpAsymmPhip595MeVCM2");
  ParaPerpAsymmPhip_595MeVCM2_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 590-600MeV Photon Energy -ve Hel (CosTheta0.8-0.6)");
  ParaPerpAsymmPhip_595MeVCM2_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[1][17] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[1][17] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_605MeVCM2_PosHel = Phip_605MeVCM2_Para_PosHel->GetAsymmetry(Phip_605MeVCM2_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_605MeVCM2_PosHel->SetName("ParaPerpAsymmPhip605MeVCM2");
  ParaPerpAsymmPhip_605MeVCM2_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 600-610MeV Photon Energy -ve Hel (CosTheta0.8-0.6)");
  ParaPerpAsymmPhip_605MeVCM2_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[1][18] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[1][18] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_615MeVCM2_PosHel = Phip_615MeVCM2_Para_PosHel->GetAsymmetry(Phip_615MeVCM2_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_615MeVCM2_PosHel->SetName("ParaPerpAsymmPhip615MeVCM2");
  ParaPerpAsymmPhip_615MeVCM2_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 610-620MeV Photon Energy -ve Hel (CosTheta0.8-0.6)");
  ParaPerpAsymmPhip_615MeVCM2_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[1][19] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[1][19] = CosFit->GetParError(0);

  ///////////////////////////////////////////
  //////////////////CM3 p ///////////////////
  ///////////////////////////////////////////

  ParaPerpAsymmPhip_425MeVCM3_PosHel = Phip_425MeVCM3_Para_PosHel->GetAsymmetry(Phip_425MeVCM3_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_425MeVCM3_PosHel->SetName("ParaPerpAsymmPhip425MeVCM3");
  ParaPerpAsymmPhip_425MeVCM3_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 420-430MeV Photon Energy -ve Hel (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_425MeVCM3_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[2][0] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[2][0] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_435MeVCM3_PosHel = Phip_435MeVCM3_Para_PosHel->GetAsymmetry(Phip_435MeVCM3_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_435MeVCM3_PosHel->SetName("ParaPerpAsymmPhip435MeVCM3");
  ParaPerpAsymmPhip_435MeVCM3_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 430-440MeV Photon Energy -ve Hel (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_435MeVCM3_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[2][1] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[2][1] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_445MeVCM3_PosHel = Phip_445MeVCM3_Para_PosHel->GetAsymmetry(Phip_445MeVCM3_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_445MeVCM3_PosHel->SetName("ParaPerpAsymmPhip445MeVCM3");
  ParaPerpAsymmPhip_445MeVCM3_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 440-450MeV Photon Energy -ve Hel (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_445MeVCM3_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[2][2] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[2][2] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_455MeVCM3_PosHel = Phip_455MeVCM3_Para_PosHel->GetAsymmetry(Phip_455MeVCM3_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_455MeVCM3_PosHel->SetName("ParaPerpAsymmPhip455MeVCM3");
  ParaPerpAsymmPhip_455MeVCM3_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 450-460MeV Photon Energy -ve Hel (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_455MeVCM3_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[2][3] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[2][3] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_465MeVCM3_PosHel = Phip_465MeVCM3_Para_PosHel->GetAsymmetry(Phip_465MeVCM3_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_465MeVCM3_PosHel->SetName("ParaPerpAsymmPhip465MeVCM3");
  ParaPerpAsymmPhip_465MeVCM3_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 460-470MeV Photon Energy -ve Hel (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_465MeVCM3_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[2][4] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[2][4] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_475MeVCM3_PosHel = Phip_475MeVCM3_Para_PosHel->GetAsymmetry(Phip_475MeVCM3_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_475MeVCM3_PosHel->SetName("ParaPerpAsymmPhip475MeVCM3");
  ParaPerpAsymmPhip_475MeVCM3_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 470-480MeV Photon Energy -ve Hel (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_475MeVCM3_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[2][5] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[2][5] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_485MeVCM3_PosHel = Phip_485MeVCM3_Para_PosHel->GetAsymmetry(Phip_485MeVCM3_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_485MeVCM3_PosHel->SetName("ParaPerpAsymmPhip485MeVCM3");
  ParaPerpAsymmPhip_485MeVCM3_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 480-490MeV Photon Energy -ve Hel (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_485MeVCM3_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[2][6] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[2][6] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_495MeVCM3_PosHel = Phip_495MeVCM3_Para_PosHel->GetAsymmetry(Phip_495MeVCM3_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_495MeVCM3_PosHel->SetName("ParaPerpAsymmPhip495MeVCM3");
  ParaPerpAsymmPhip_495MeVCM3_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 490-500MeV Photon Energy -ve Hel (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_495MeVCM3_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[2][7] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[2][7] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_505MeVCM3_PosHel = Phip_505MeVCM3_Para_PosHel->GetAsymmetry(Phip_505MeVCM3_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_505MeVCM3_PosHel->SetName("ParaPerpAsymmPhip505MeVCM3");
  ParaPerpAsymmPhip_505MeVCM3_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 500-510MeV Photon Energy -ve Hel (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_505MeVCM3_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[2][8] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[2][8] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_515MeVCM3_PosHel = Phip_515MeVCM3_Para_PosHel->GetAsymmetry(Phip_515MeVCM3_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_515MeVCM3_PosHel->SetName("ParaPerpAsymmPhip515MeVCM3");
  ParaPerpAsymmPhip_515MeVCM3_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 510-520MeV Photon Energy -ve Hel (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_515MeVCM3_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[2][9] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[2][9] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_525MeVCM3_PosHel = Phip_525MeVCM3_Para_PosHel->GetAsymmetry(Phip_525MeVCM3_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_525MeVCM3_PosHel->SetName("ParaPerpAsymmPhip525MeVCM3");
  ParaPerpAsymmPhip_525MeVCM3_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 520-530MeV Photon Energy -ve Hel (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_525MeVCM3_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[2][10] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[2][10] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_535MeVCM3_PosHel = Phip_535MeVCM3_Para_PosHel->GetAsymmetry(Phip_535MeVCM3_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_535MeVCM3_PosHel->SetName("ParaPerpAsymmPhip535MeVCM3");
  ParaPerpAsymmPhip_535MeVCM3_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 530-540MeV Photon Energy -ve Hel (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_535MeVCM3_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[2][11] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[2][11] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_545MeVCM3_PosHel = Phip_545MeVCM3_Para_PosHel->GetAsymmetry(Phip_545MeVCM3_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_545MeVCM3_PosHel->SetName("ParaPerpAsymmPhip545MeVCM3");
  ParaPerpAsymmPhip_545MeVCM3_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 540-550MeV Photon Energy -ve Hel (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_545MeVCM3_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[2][12] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[2][12] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_555MeVCM3_PosHel = Phip_555MeVCM3_Para_PosHel->GetAsymmetry(Phip_555MeVCM3_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_555MeVCM3_PosHel->SetName("ParaPerpAsymmPhip555MeVCM3");
  ParaPerpAsymmPhip_555MeVCM3_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 550-560MeV Photon Energy -ve Hel (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_555MeVCM3_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[2][13] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[2][13] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_565MeVCM3_PosHel = Phip_565MeVCM3_Para_PosHel->GetAsymmetry(Phip_565MeVCM3_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_565MeVCM3_PosHel->SetName("ParaPerpAsymmPhip565MeVCM3");
  ParaPerpAsymmPhip_565MeVCM3_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 560-570MeV Photon Energy -ve Hel (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_565MeVCM3_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[2][14] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[2][14] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_575MeVCM3_PosHel = Phip_575MeVCM3_Para_PosHel->GetAsymmetry(Phip_575MeVCM3_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_575MeVCM3_PosHel->SetName("ParaPerpAsymmPhip575MeVCM3");
  ParaPerpAsymmPhip_575MeVCM3_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 570-580MeV Photon Energy -ve Hel (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_575MeVCM3_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[2][15] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[2][15] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_585MeVCM3_PosHel = Phip_585MeVCM3_Para_PosHel->GetAsymmetry(Phip_585MeVCM3_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_585MeVCM3_PosHel->SetName("ParaPerpAsymmPhip585MeVCM3");
  ParaPerpAsymmPhip_585MeVCM3_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 580-590MeV Photon Energy -ve Hel (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_585MeVCM3_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[2][16] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[2][16] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_595MeVCM3_PosHel = Phip_595MeVCM3_Para_PosHel->GetAsymmetry(Phip_595MeVCM3_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_595MeVCM3_PosHel->SetName("ParaPerpAsymmPhip595MeVCM3");
  ParaPerpAsymmPhip_595MeVCM3_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 590-600MeV Photon Energy -ve Hel (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_595MeVCM3_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[2][17] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[2][17] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_605MeVCM3_PosHel = Phip_605MeVCM3_Para_PosHel->GetAsymmetry(Phip_605MeVCM3_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_605MeVCM3_PosHel->SetName("ParaPerpAsymmPhip605MeVCM3");
  ParaPerpAsymmPhip_605MeVCM3_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 600-610MeV Photon Energy -ve Hel (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_605MeVCM3_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[2][18] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[2][18] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_615MeVCM3_PosHel = Phip_615MeVCM3_Para_PosHel->GetAsymmetry(Phip_615MeVCM3_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_615MeVCM3_PosHel->SetName("ParaPerpAsymmPhip615MeVCM3");
  ParaPerpAsymmPhip_615MeVCM3_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 610-620MeV Photon Energy -ve Hel (CosTheta0.6-0.4)");
  ParaPerpAsymmPhip_615MeVCM3_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[2][19] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[2][19] = CosFit->GetParError(0);

  ///////////////////////////////////////////
  //////////////////CM4 p ///////////////////
  ///////////////////////////////////////////

  ParaPerpAsymmPhip_425MeVCM4_PosHel = Phip_425MeVCM4_Para_PosHel->GetAsymmetry(Phip_425MeVCM4_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_425MeVCM4_PosHel->SetName("ParaPerpAsymmPhip425MeVCM4");
  ParaPerpAsymmPhip_425MeVCM4_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 420-430MeV Photon Energy -ve Hel (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_425MeVCM4_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[3][0] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[3][0] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_435MeVCM4_PosHel = Phip_435MeVCM4_Para_PosHel->GetAsymmetry(Phip_435MeVCM4_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_435MeVCM4_PosHel->SetName("ParaPerpAsymmPhip435MeVCM4");
  ParaPerpAsymmPhip_435MeVCM4_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 430-440MeV Photon Energy -ve Hel (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_435MeVCM4_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[3][1] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[3][1] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_445MeVCM4_PosHel = Phip_445MeVCM4_Para_PosHel->GetAsymmetry(Phip_445MeVCM4_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_445MeVCM4_PosHel->SetName("ParaPerpAsymmPhip445MeVCM4");
  ParaPerpAsymmPhip_445MeVCM4_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 440-450MeV Photon Energy -ve Hel (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_445MeVCM4_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[3][2] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[3][2] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_455MeVCM4_PosHel = Phip_455MeVCM4_Para_PosHel->GetAsymmetry(Phip_455MeVCM4_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_455MeVCM4_PosHel->SetName("ParaPerpAsymmPhip455MeVCM4");
  ParaPerpAsymmPhip_455MeVCM4_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 450-460MeV Photon Energy -ve Hel (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_455MeVCM4_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[3][3] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[3][3] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_465MeVCM4_PosHel = Phip_465MeVCM4_Para_PosHel->GetAsymmetry(Phip_465MeVCM4_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_465MeVCM4_PosHel->SetName("ParaPerpAsymmPhip465MeVCM4");
  ParaPerpAsymmPhip_465MeVCM4_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 460-470MeV Photon Energy -ve Hel (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_465MeVCM4_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[3][4] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[3][4] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_475MeVCM4_PosHel = Phip_475MeVCM4_Para_PosHel->GetAsymmetry(Phip_475MeVCM4_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_475MeVCM4_PosHel->SetName("ParaPerpAsymmPhip475MeVCM4");
  ParaPerpAsymmPhip_475MeVCM4_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 470-480MeV Photon Energy -ve Hel (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_475MeVCM4_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[3][5] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[3][5] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_485MeVCM4_PosHel = Phip_485MeVCM4_Para_PosHel->GetAsymmetry(Phip_485MeVCM4_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_485MeVCM4_PosHel->SetName("ParaPerpAsymmPhip485MeVCM4");
  ParaPerpAsymmPhip_485MeVCM4_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 480-490MeV Photon Energy -ve Hel (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_485MeVCM4_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[3][6] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[3][6] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_495MeVCM4_PosHel = Phip_495MeVCM4_Para_PosHel->GetAsymmetry(Phip_495MeVCM4_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_495MeVCM4_PosHel->SetName("ParaPerpAsymmPhip495MeVCM4");
  ParaPerpAsymmPhip_495MeVCM4_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 490-500MeV Photon Energy -ve Hel (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_495MeVCM4_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[3][7] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[3][7] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_505MeVCM4_PosHel = Phip_505MeVCM4_Para_PosHel->GetAsymmetry(Phip_505MeVCM4_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_505MeVCM4_PosHel->SetName("ParaPerpAsymmPhip505MeVCM4");
  ParaPerpAsymmPhip_505MeVCM4_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 500-510MeV Photon Energy -ve Hel (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_505MeVCM4_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[3][8] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[3][8] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_515MeVCM4_PosHel = Phip_515MeVCM4_Para_PosHel->GetAsymmetry(Phip_515MeVCM4_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_515MeVCM4_PosHel->SetName("ParaPerpAsymmPhip515MeVCM4");
  ParaPerpAsymmPhip_515MeVCM4_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 510-520MeV Photon Energy -ve Hel (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_515MeVCM4_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[3][9] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[3][9] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_525MeVCM4_PosHel = Phip_525MeVCM4_Para_PosHel->GetAsymmetry(Phip_525MeVCM4_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_525MeVCM4_PosHel->SetName("ParaPerpAsymmPhip525MeVCM4");
  ParaPerpAsymmPhip_525MeVCM4_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 520-530MeV Photon Energy -ve Hel (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_525MeVCM4_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[3][10] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[3][10] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_535MeVCM4_PosHel = Phip_535MeVCM4_Para_PosHel->GetAsymmetry(Phip_535MeVCM4_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_535MeVCM4_PosHel->SetName("ParaPerpAsymmPhip535MeVCM4");
  ParaPerpAsymmPhip_535MeVCM4_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 530-540MeV Photon Energy -ve Hel (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_535MeVCM4_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[3][11] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[3][11] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_545MeVCM4_PosHel = Phip_545MeVCM4_Para_PosHel->GetAsymmetry(Phip_545MeVCM4_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_545MeVCM4_PosHel->SetName("ParaPerpAsymmPhip545MeVCM4");
  ParaPerpAsymmPhip_545MeVCM4_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 540-550MeV Photon Energy -ve Hel (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_545MeVCM4_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[3][12] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[3][12] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_555MeVCM4_PosHel = Phip_555MeVCM4_Para_PosHel->GetAsymmetry(Phip_555MeVCM4_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_555MeVCM4_PosHel->SetName("ParaPerpAsymmPhip555MeVCM4");
  ParaPerpAsymmPhip_555MeVCM4_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 550-560MeV Photon Energy -ve Hel (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_555MeVCM4_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[3][13] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[3][13] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_565MeVCM4_PosHel = Phip_565MeVCM4_Para_PosHel->GetAsymmetry(Phip_565MeVCM4_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_565MeVCM4_PosHel->SetName("ParaPerpAsymmPhip565MeVCM4");
  ParaPerpAsymmPhip_565MeVCM4_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 560-570MeV Photon Energy -ve Hel (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_565MeVCM4_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[3][14] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[3][14] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_575MeVCM4_PosHel = Phip_575MeVCM4_Para_PosHel->GetAsymmetry(Phip_575MeVCM4_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_575MeVCM4_PosHel->SetName("ParaPerpAsymmPhip575MeVCM4");
  ParaPerpAsymmPhip_575MeVCM4_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 570-580MeV Photon Energy -ve Hel (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_575MeVCM4_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[3][15] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[3][15] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_585MeVCM4_PosHel = Phip_585MeVCM4_Para_PosHel->GetAsymmetry(Phip_585MeVCM4_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_585MeVCM4_PosHel->SetName("ParaPerpAsymmPhip585MeVCM4");
  ParaPerpAsymmPhip_585MeVCM4_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 580-590MeV Photon Energy -ve Hel (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_585MeVCM4_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[3][16] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[3][16] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_595MeVCM4_PosHel = Phip_595MeVCM4_Para_PosHel->GetAsymmetry(Phip_595MeVCM4_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_595MeVCM4_PosHel->SetName("ParaPerpAsymmPhip595MeVCM4");
  ParaPerpAsymmPhip_595MeVCM4_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 590-600MeV Photon Energy -ve Hel (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_595MeVCM4_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[3][17] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[3][17] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_605MeVCM4_PosHel = Phip_605MeVCM4_Para_PosHel->GetAsymmetry(Phip_605MeVCM4_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_605MeVCM4_PosHel->SetName("ParaPerpAsymmPhip605MeVCM4");
  ParaPerpAsymmPhip_605MeVCM4_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 600-610MeV Photon Energy -ve Hel (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_605MeVCM4_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[3][18] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[3][18] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_615MeVCM4_PosHel = Phip_615MeVCM4_Para_PosHel->GetAsymmetry(Phip_615MeVCM4_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_615MeVCM4_PosHel->SetName("ParaPerpAsymmPhip615MeVCM4");
  ParaPerpAsymmPhip_615MeVCM4_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 610-620MeV Photon Energy -ve Hel (CosTheta0.4-0.2)");
  ParaPerpAsymmPhip_615MeVCM4_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[3][19] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[3][19] = CosFit->GetParError(0);

  ///////////////////////////////////////////
  //////////////////CM5 p ///////////////////
  ///////////////////////////////////////////

  ParaPerpAsymmPhip_425MeVCM5_PosHel = Phip_425MeVCM5_Para_PosHel->GetAsymmetry(Phip_425MeVCM5_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_425MeVCM5_PosHel->SetName("ParaPerpAsymmPhip425MeVCM5");
  ParaPerpAsymmPhip_425MeVCM5_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 420-430MeV Photon Energy -ve Hel (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_425MeVCM5_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[4][0] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[4][0] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_435MeVCM5_PosHel = Phip_435MeVCM5_Para_PosHel->GetAsymmetry(Phip_435MeVCM5_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_435MeVCM5_PosHel->SetName("ParaPerpAsymmPhip435MeVCM5");
  ParaPerpAsymmPhip_435MeVCM5_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 430-440MeV Photon Energy -ve Hel (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_435MeVCM5_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[4][1] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[4][1] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_445MeVCM5_PosHel = Phip_445MeVCM5_Para_PosHel->GetAsymmetry(Phip_445MeVCM5_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_445MeVCM5_PosHel->SetName("ParaPerpAsymmPhip445MeVCM5");
  ParaPerpAsymmPhip_445MeVCM5_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 440-450MeV Photon Energy -ve Hel (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_445MeVCM5_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[4][2] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[4][2] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_455MeVCM5_PosHel = Phip_455MeVCM5_Para_PosHel->GetAsymmetry(Phip_455MeVCM5_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_455MeVCM5_PosHel->SetName("ParaPerpAsymmPhip455MeVCM5");
  ParaPerpAsymmPhip_455MeVCM5_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 450-460MeV Photon Energy -ve Hel (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_455MeVCM5_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[4][3] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[4][3] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_465MeVCM5_PosHel = Phip_465MeVCM5_Para_PosHel->GetAsymmetry(Phip_465MeVCM5_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_465MeVCM5_PosHel->SetName("ParaPerpAsymmPhip465MeVCM5");
  ParaPerpAsymmPhip_465MeVCM5_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 460-470MeV Photon Energy -ve Hel (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_465MeVCM5_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[4][4] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[4][4] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_475MeVCM5_PosHel = Phip_475MeVCM5_Para_PosHel->GetAsymmetry(Phip_475MeVCM5_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_475MeVCM5_PosHel->SetName("ParaPerpAsymmPhip475MeVCM5");
  ParaPerpAsymmPhip_475MeVCM5_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 470-480MeV Photon Energy -ve Hel (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_475MeVCM5_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[4][5] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[4][5] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_485MeVCM5_PosHel = Phip_485MeVCM5_Para_PosHel->GetAsymmetry(Phip_485MeVCM5_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_485MeVCM5_PosHel->SetName("ParaPerpAsymmPhip485MeVCM5");
  ParaPerpAsymmPhip_485MeVCM5_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 480-490MeV Photon Energy -ve Hel (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_485MeVCM5_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[4][6] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[4][6] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_495MeVCM5_PosHel = Phip_495MeVCM5_Para_PosHel->GetAsymmetry(Phip_495MeVCM5_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_495MeVCM5_PosHel->SetName("ParaPerpAsymmPhip495MeVCM5");
  ParaPerpAsymmPhip_495MeVCM5_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 490-500MeV Photon Energy -ve Hel (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_495MeVCM5_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[4][7] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[4][7] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_505MeVCM5_PosHel = Phip_505MeVCM5_Para_PosHel->GetAsymmetry(Phip_505MeVCM5_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_505MeVCM5_PosHel->SetName("ParaPerpAsymmPhip505MeVCM5");
  ParaPerpAsymmPhip_505MeVCM5_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 500-510MeV Photon Energy -ve Hel (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_505MeVCM5_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[4][8] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[4][8] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_515MeVCM5_PosHel = Phip_515MeVCM5_Para_PosHel->GetAsymmetry(Phip_515MeVCM5_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_515MeVCM5_PosHel->SetName("ParaPerpAsymmPhip515MeVCM5");
  ParaPerpAsymmPhip_515MeVCM5_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 510-520MeV Photon Energy -ve Hel (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_515MeVCM5_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[4][9] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[4][9] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_525MeVCM5_PosHel = Phip_525MeVCM5_Para_PosHel->GetAsymmetry(Phip_525MeVCM5_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_525MeVCM5_PosHel->SetName("ParaPerpAsymmPhip525MeVCM5");
  ParaPerpAsymmPhip_525MeVCM5_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 520-530MeV Photon Energy -ve Hel (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_525MeVCM5_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[4][10] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[4][10] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_535MeVCM5_PosHel = Phip_535MeVCM5_Para_PosHel->GetAsymmetry(Phip_535MeVCM5_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_535MeVCM5_PosHel->SetName("ParaPerpAsymmPhip535MeVCM5");
  ParaPerpAsymmPhip_535MeVCM5_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 530-540MeV Photon Energy -ve Hel (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_535MeVCM5_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[4][11] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[4][11] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_545MeVCM5_PosHel = Phip_545MeVCM5_Para_PosHel->GetAsymmetry(Phip_545MeVCM5_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_545MeVCM5_PosHel->SetName("ParaPerpAsymmPhip545MeVCM5");
  ParaPerpAsymmPhip_545MeVCM5_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 540-550MeV Photon Energy -ve Hel (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_545MeVCM5_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[4][12] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[4][12] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_555MeVCM5_PosHel = Phip_555MeVCM5_Para_PosHel->GetAsymmetry(Phip_555MeVCM5_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_555MeVCM5_PosHel->SetName("ParaPerpAsymmPhip555MeVCM5");
  ParaPerpAsymmPhip_555MeVCM5_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 550-560MeV Photon Energy -ve Hel (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_555MeVCM5_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[4][13] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[4][13] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_565MeVCM5_PosHel = Phip_565MeVCM5_Para_PosHel->GetAsymmetry(Phip_565MeVCM5_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_565MeVCM5_PosHel->SetName("ParaPerpAsymmPhip565MeVCM5");
  ParaPerpAsymmPhip_565MeVCM5_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 560-570MeV Photon Energy -ve Hel (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_565MeVCM5_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[4][14] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[4][14] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_575MeVCM5_PosHel = Phip_575MeVCM5_Para_PosHel->GetAsymmetry(Phip_575MeVCM5_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_575MeVCM5_PosHel->SetName("ParaPerpAsymmPhip575MeVCM5");
  ParaPerpAsymmPhip_575MeVCM5_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 570-580MeV Photon Energy -ve Hel (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_575MeVCM5_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[4][15] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[4][15] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_585MeVCM5_PosHel = Phip_585MeVCM5_Para_PosHel->GetAsymmetry(Phip_585MeVCM5_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_585MeVCM5_PosHel->SetName("ParaPerpAsymmPhip585MeVCM5");
  ParaPerpAsymmPhip_585MeVCM5_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 580-590MeV Photon Energy -ve Hel (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_585MeVCM5_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[4][16] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[4][16] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_595MeVCM5_PosHel = Phip_595MeVCM5_Para_PosHel->GetAsymmetry(Phip_595MeVCM5_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_595MeVCM5_PosHel->SetName("ParaPerpAsymmPhip595MeVCM5");
  ParaPerpAsymmPhip_595MeVCM5_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 590-600MeV Photon Energy -ve Hel (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_595MeVCM5_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[4][17] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[4][17] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_605MeVCM5_PosHel = Phip_605MeVCM5_Para_PosHel->GetAsymmetry(Phip_605MeVCM5_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_605MeVCM5_PosHel->SetName("ParaPerpAsymmPhip605MeVCM5");
  ParaPerpAsymmPhip_605MeVCM5_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 600-610MeV Photon Energy -ve Hel (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_605MeVCM5_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[4][18] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[4][18] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_615MeVCM5_PosHel = Phip_615MeVCM5_Para_PosHel->GetAsymmetry(Phip_615MeVCM5_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_615MeVCM5_PosHel->SetName("ParaPerpAsymmPhip615MeVCM5");
  ParaPerpAsymmPhip_615MeVCM5_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 610-620MeV Photon Energy -ve Hel (CosTheta0.2-0.0)");
  ParaPerpAsymmPhip_615MeVCM5_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[4][19] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[4][19] = CosFit->GetParError(0);

  ///////////////////////////////////////////
  //////////////////CM6 p ///////////////////
  ///////////////////////////////////////////

  ParaPerpAsymmPhip_425MeVCM6_PosHel = Phip_425MeVCM6_Para_PosHel->GetAsymmetry(Phip_425MeVCM6_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_425MeVCM6_PosHel->SetName("ParaPerpAsymmPhip425MeVCM6");
  ParaPerpAsymmPhip_425MeVCM6_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 420-430MeV Photon Energy -ve Hel (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_425MeVCM6_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[5][0] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[5][0] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_435MeVCM6_PosHel = Phip_435MeVCM6_Para_PosHel->GetAsymmetry(Phip_435MeVCM6_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_435MeVCM6_PosHel->SetName("ParaPerpAsymmPhip435MeVCM6");
  ParaPerpAsymmPhip_435MeVCM6_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 430-440MeV Photon Energy -ve Hel (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_435MeVCM6_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[5][1] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[5][1] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_445MeVCM6_PosHel = Phip_445MeVCM6_Para_PosHel->GetAsymmetry(Phip_445MeVCM6_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_445MeVCM6_PosHel->SetName("ParaPerpAsymmPhip445MeVCM6");
  ParaPerpAsymmPhip_445MeVCM6_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 440-450MeV Photon Energy -ve Hel (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_445MeVCM6_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[5][2] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[5][2] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_455MeVCM6_PosHel = Phip_455MeVCM6_Para_PosHel->GetAsymmetry(Phip_455MeVCM6_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_455MeVCM6_PosHel->SetName("ParaPerpAsymmPhip455MeVCM6");
  ParaPerpAsymmPhip_455MeVCM6_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 450-460MeV Photon Energy -ve Hel (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_455MeVCM6_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[5][3] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[5][3] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_465MeVCM6_PosHel = Phip_465MeVCM6_Para_PosHel->GetAsymmetry(Phip_465MeVCM6_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_465MeVCM6_PosHel->SetName("ParaPerpAsymmPhip465MeVCM6");
  ParaPerpAsymmPhip_465MeVCM6_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 460-470MeV Photon Energy -ve Hel (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_465MeVCM6_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[5][4] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[5][4] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_475MeVCM6_PosHel = Phip_475MeVCM6_Para_PosHel->GetAsymmetry(Phip_475MeVCM6_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_475MeVCM6_PosHel->SetName("ParaPerpAsymmPhip475MeVCM6");
  ParaPerpAsymmPhip_475MeVCM6_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 470-480MeV Photon Energy -ve Hel (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_475MeVCM6_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[5][5] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[5][5] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_485MeVCM6_PosHel = Phip_485MeVCM6_Para_PosHel->GetAsymmetry(Phip_485MeVCM6_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_485MeVCM6_PosHel->SetName("ParaPerpAsymmPhip485MeVCM6");
  ParaPerpAsymmPhip_485MeVCM6_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 480-490MeV Photon Energy -ve Hel (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_485MeVCM6_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[5][6] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[5][6] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_495MeVCM6_PosHel = Phip_495MeVCM6_Para_PosHel->GetAsymmetry(Phip_495MeVCM6_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_495MeVCM6_PosHel->SetName("ParaPerpAsymmPhip495MeVCM6");
  ParaPerpAsymmPhip_495MeVCM6_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 490-500MeV Photon Energy -ve Hel (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_495MeVCM6_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[5][7] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[5][7] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_505MeVCM6_PosHel = Phip_505MeVCM6_Para_PosHel->GetAsymmetry(Phip_505MeVCM6_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_505MeVCM6_PosHel->SetName("ParaPerpAsymmPhip505MeVCM6");
  ParaPerpAsymmPhip_505MeVCM6_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 500-510MeV Photon Energy -ve Hel (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_505MeVCM6_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[5][8] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[5][8] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_515MeVCM6_PosHel = Phip_515MeVCM6_Para_PosHel->GetAsymmetry(Phip_515MeVCM6_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_515MeVCM6_PosHel->SetName("ParaPerpAsymmPhip515MeVCM6");
  ParaPerpAsymmPhip_515MeVCM6_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 510-520MeV Photon Energy -ve Hel (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_515MeVCM6_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[5][9] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[5][9] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_525MeVCM6_PosHel = Phip_525MeVCM6_Para_PosHel->GetAsymmetry(Phip_525MeVCM6_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_525MeVCM6_PosHel->SetName("ParaPerpAsymmPhip525MeVCM6");
  ParaPerpAsymmPhip_525MeVCM6_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 520-530MeV Photon Energy -ve Hel (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_525MeVCM6_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[5][10] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[5][10] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_535MeVCM6_PosHel = Phip_535MeVCM6_Para_PosHel->GetAsymmetry(Phip_535MeVCM6_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_535MeVCM6_PosHel->SetName("ParaPerpAsymmPhip535MeVCM6");
  ParaPerpAsymmPhip_535MeVCM6_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 530-540MeV Photon Energy -ve Hel (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_535MeVCM6_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[5][11] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[5][11] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_545MeVCM6_PosHel = Phip_545MeVCM6_Para_PosHel->GetAsymmetry(Phip_545MeVCM6_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_545MeVCM6_PosHel->SetName("ParaPerpAsymmPhip545MeVCM6");
  ParaPerpAsymmPhip_545MeVCM6_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 540-550MeV Photon Energy -ve Hel (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_545MeVCM6_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[5][12] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[5][12] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_555MeVCM6_PosHel = Phip_555MeVCM6_Para_PosHel->GetAsymmetry(Phip_555MeVCM6_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_555MeVCM6_PosHel->SetName("ParaPerpAsymmPhip555MeVCM6");
  ParaPerpAsymmPhip_555MeVCM6_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 550-560MeV Photon Energy -ve Hel (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_555MeVCM6_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[5][13] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[5][13] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_565MeVCM6_PosHel = Phip_565MeVCM6_Para_PosHel->GetAsymmetry(Phip_565MeVCM6_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_565MeVCM6_PosHel->SetName("ParaPerpAsymmPhip565MeVCM6");
  ParaPerpAsymmPhip_565MeVCM6_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 560-570MeV Photon Energy -ve Hel (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_565MeVCM6_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[5][14] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[5][14] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_575MeVCM6_PosHel = Phip_575MeVCM6_Para_PosHel->GetAsymmetry(Phip_575MeVCM6_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_575MeVCM6_PosHel->SetName("ParaPerpAsymmPhip575MeVCM6");
  ParaPerpAsymmPhip_575MeVCM6_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 570-580MeV Photon Energy -ve Hel (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_575MeVCM6_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[5][15] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[5][15] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_585MeVCM6_PosHel = Phip_585MeVCM6_Para_PosHel->GetAsymmetry(Phip_585MeVCM6_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_585MeVCM6_PosHel->SetName("ParaPerpAsymmPhip585MeVCM6");
  ParaPerpAsymmPhip_585MeVCM6_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 580-590MeV Photon Energy -ve Hel (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_585MeVCM6_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[5][16] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[5][16] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_595MeVCM6_PosHel = Phip_595MeVCM6_Para_PosHel->GetAsymmetry(Phip_595MeVCM6_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_595MeVCM6_PosHel->SetName("ParaPerpAsymmPhip595MeVCM6");
  ParaPerpAsymmPhip_595MeVCM6_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 590-600MeV Photon Energy -ve Hel (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_595MeVCM6_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[5][17] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[5][17] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_605MeVCM6_PosHel = Phip_605MeVCM6_Para_PosHel->GetAsymmetry(Phip_605MeVCM6_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_605MeVCM6_PosHel->SetName("ParaPerpAsymmPhip605MeVCM6");
  ParaPerpAsymmPhip_605MeVCM6_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 600-610MeV Photon Energy -ve Hel (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_605MeVCM6_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[5][18] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[5][18] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_615MeVCM6_PosHel = Phip_615MeVCM6_Para_PosHel->GetAsymmetry(Phip_615MeVCM6_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_615MeVCM6_PosHel->SetName("ParaPerpAsymmPhip615MeVCM6");
  ParaPerpAsymmPhip_615MeVCM6_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 610-620MeV Photon Energy -ve Hel (CosTheta0.0-(-0.2))");
  ParaPerpAsymmPhip_615MeVCM6_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[5][19] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[5][19] = CosFit->GetParError(0);

  ///////////////////////////////////////////
  //////////////////CM7 p ///////////////////
  ///////////////////////////////////////////

  ParaPerpAsymmPhip_425MeVCM7_PosHel = Phip_425MeVCM7_Para_PosHel->GetAsymmetry(Phip_425MeVCM7_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_425MeVCM7_PosHel->SetName("ParaPerpAsymmPhip425MeVCM7");
  ParaPerpAsymmPhip_425MeVCM7_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 420-430MeV Photon Energy -ve Hel (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_425MeVCM7_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[6][0] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[6][0] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_435MeVCM7_PosHel = Phip_435MeVCM7_Para_PosHel->GetAsymmetry(Phip_435MeVCM7_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_435MeVCM7_PosHel->SetName("ParaPerpAsymmPhip435MeVCM7");
  ParaPerpAsymmPhip_435MeVCM7_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 430-440MeV Photon Energy -ve Hel (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_435MeVCM7_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[6][1] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[6][1] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_445MeVCM7_PosHel = Phip_445MeVCM7_Para_PosHel->GetAsymmetry(Phip_445MeVCM7_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_445MeVCM7_PosHel->SetName("ParaPerpAsymmPhip445MeVCM7");
  ParaPerpAsymmPhip_445MeVCM7_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 440-450MeV Photon Energy -ve Hel (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_445MeVCM7_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[6][2] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[6][2] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_455MeVCM7_PosHel = Phip_455MeVCM7_Para_PosHel->GetAsymmetry(Phip_455MeVCM7_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_455MeVCM7_PosHel->SetName("ParaPerpAsymmPhip455MeVCM7");
  ParaPerpAsymmPhip_455MeVCM7_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 450-460MeV Photon Energy -ve Hel (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_455MeVCM7_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[6][3] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[6][3] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_465MeVCM7_PosHel = Phip_465MeVCM7_Para_PosHel->GetAsymmetry(Phip_465MeVCM7_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_465MeVCM7_PosHel->SetName("ParaPerpAsymmPhip465MeVCM7");
  ParaPerpAsymmPhip_465MeVCM7_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 460-470MeV Photon Energy -ve Hel (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_465MeVCM7_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[6][4] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[6][4] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_475MeVCM7_PosHel = Phip_475MeVCM7_Para_PosHel->GetAsymmetry(Phip_475MeVCM7_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_475MeVCM7_PosHel->SetName("ParaPerpAsymmPhip475MeVCM7");
  ParaPerpAsymmPhip_475MeVCM7_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 470-480MeV Photon Energy -ve Hel (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_475MeVCM7_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[6][5] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[6][5] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_485MeVCM7_PosHel = Phip_485MeVCM7_Para_PosHel->GetAsymmetry(Phip_485MeVCM7_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_485MeVCM7_PosHel->SetName("ParaPerpAsymmPhip485MeVCM7");
  ParaPerpAsymmPhip_485MeVCM7_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 480-490MeV Photon Energy -ve Hel (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_485MeVCM7_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[6][6] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[6][6] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_495MeVCM7_PosHel = Phip_495MeVCM7_Para_PosHel->GetAsymmetry(Phip_495MeVCM7_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_495MeVCM7_PosHel->SetName("ParaPerpAsymmPhip495MeVCM7");
  ParaPerpAsymmPhip_495MeVCM7_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 490-500MeV Photon Energy -ve Hel (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_495MeVCM7_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[6][7] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[6][7] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_505MeVCM7_PosHel = Phip_505MeVCM7_Para_PosHel->GetAsymmetry(Phip_505MeVCM7_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_505MeVCM7_PosHel->SetName("ParaPerpAsymmPhip505MeVCM7");
  ParaPerpAsymmPhip_505MeVCM7_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 500-510MeV Photon Energy -ve Hel (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_505MeVCM7_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[6][8] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[6][8] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_515MeVCM7_PosHel = Phip_515MeVCM7_Para_PosHel->GetAsymmetry(Phip_515MeVCM7_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_515MeVCM7_PosHel->SetName("ParaPerpAsymmPhip515MeVCM7");
  ParaPerpAsymmPhip_515MeVCM7_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 510-520MeV Photon Energy -ve Hel (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_515MeVCM7_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[6][9] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[6][9] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_525MeVCM7_PosHel = Phip_525MeVCM7_Para_PosHel->GetAsymmetry(Phip_525MeVCM7_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_525MeVCM7_PosHel->SetName("ParaPerpAsymmPhip525MeVCM7");
  ParaPerpAsymmPhip_525MeVCM7_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 520-530MeV Photon Energy -ve Hel (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_525MeVCM7_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[6][10] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[6][10] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_535MeVCM7_PosHel = Phip_535MeVCM7_Para_PosHel->GetAsymmetry(Phip_535MeVCM7_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_535MeVCM7_PosHel->SetName("ParaPerpAsymmPhip535MeVCM7");
  ParaPerpAsymmPhip_535MeVCM7_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 530-540MeV Photon Energy -ve Hel (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_535MeVCM7_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[6][11] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[6][11] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_545MeVCM7_PosHel = Phip_545MeVCM7_Para_PosHel->GetAsymmetry(Phip_545MeVCM7_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_545MeVCM7_PosHel->SetName("ParaPerpAsymmPhip545MeVCM7");
  ParaPerpAsymmPhip_545MeVCM7_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 540-550MeV Photon Energy -ve Hel (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_545MeVCM7_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[6][12] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[6][12] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_555MeVCM7_PosHel = Phip_555MeVCM7_Para_PosHel->GetAsymmetry(Phip_555MeVCM7_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_555MeVCM7_PosHel->SetName("ParaPerpAsymmPhip555MeVCM7");
  ParaPerpAsymmPhip_555MeVCM7_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 550-560MeV Photon Energy -ve Hel (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_555MeVCM7_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[6][13] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[6][13] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_565MeVCM7_PosHel = Phip_565MeVCM7_Para_PosHel->GetAsymmetry(Phip_565MeVCM7_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_565MeVCM7_PosHel->SetName("ParaPerpAsymmPhip565MeVCM7");
  ParaPerpAsymmPhip_565MeVCM7_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 560-570MeV Photon Energy -ve Hel (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_565MeVCM7_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[6][14] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[6][14] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_575MeVCM7_PosHel = Phip_575MeVCM7_Para_PosHel->GetAsymmetry(Phip_575MeVCM7_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_575MeVCM7_PosHel->SetName("ParaPerpAsymmPhip575MeVCM7");
  ParaPerpAsymmPhip_575MeVCM7_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 570-580MeV Photon Energy -ve Hel (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_575MeVCM7_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[6][15] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[6][15] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_585MeVCM7_PosHel = Phip_585MeVCM7_Para_PosHel->GetAsymmetry(Phip_585MeVCM7_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_585MeVCM7_PosHel->SetName("ParaPerpAsymmPhip585MeVCM7");
  ParaPerpAsymmPhip_585MeVCM7_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 580-590MeV Photon Energy -ve Hel (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_585MeVCM7_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[6][16] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[6][16] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_595MeVCM7_PosHel = Phip_595MeVCM7_Para_PosHel->GetAsymmetry(Phip_595MeVCM7_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_595MeVCM7_PosHel->SetName("ParaPerpAsymmPhip595MeVCM7");
  ParaPerpAsymmPhip_595MeVCM7_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 590-600MeV Photon Energy -ve Hel (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_595MeVCM7_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[6][17] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[6][17] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_605MeVCM7_PosHel = Phip_605MeVCM7_Para_PosHel->GetAsymmetry(Phip_605MeVCM7_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_605MeVCM7_PosHel->SetName("ParaPerpAsymmPhip605MeVCM7");
  ParaPerpAsymmPhip_605MeVCM7_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 600-610MeV Photon Energy -ve Hel (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_605MeVCM7_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[6][18] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[6][18] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_615MeVCM7_PosHel = Phip_615MeVCM7_Para_PosHel->GetAsymmetry(Phip_615MeVCM7_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_615MeVCM7_PosHel->SetName("ParaPerpAsymmPhip615MeVCM7");
  ParaPerpAsymmPhip_615MeVCM7_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 610-620MeV Photon Energy -ve Hel (CosTheta-0.2-(-0.4))");
  ParaPerpAsymmPhip_615MeVCM7_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[6][19] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[6][19] = CosFit->GetParError(0);

  ///////////////////////////////////////////
  //////////////////CM8 p ///////////////////
  ///////////////////////////////////////////

  ParaPerpAsymmPhip_425MeVCM8_PosHel = Phip_425MeVCM8_Para_PosHel->GetAsymmetry(Phip_425MeVCM8_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_425MeVCM8_PosHel->SetName("ParaPerpAsymmPhip425MeVCM8");
  ParaPerpAsymmPhip_425MeVCM8_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 420-430MeV Photon Energy -ve Hel (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_425MeVCM8_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[7][0] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[7][0] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_435MeVCM8_PosHel = Phip_435MeVCM8_Para_PosHel->GetAsymmetry(Phip_435MeVCM8_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_435MeVCM8_PosHel->SetName("ParaPerpAsymmPhip435MeVCM8");
  ParaPerpAsymmPhip_435MeVCM8_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 430-440MeV Photon Energy -ve Hel (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_435MeVCM8_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[7][1] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[7][1] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_445MeVCM8_PosHel = Phip_445MeVCM8_Para_PosHel->GetAsymmetry(Phip_445MeVCM8_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_445MeVCM8_PosHel->SetName("ParaPerpAsymmPhip445MeVCM8");
  ParaPerpAsymmPhip_445MeVCM8_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 440-450MeV Photon Energy -ve Hel (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_445MeVCM8_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[7][2] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[7][2] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_455MeVCM8_PosHel = Phip_455MeVCM8_Para_PosHel->GetAsymmetry(Phip_455MeVCM8_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_455MeVCM8_PosHel->SetName("ParaPerpAsymmPhip455MeVCM8");
  ParaPerpAsymmPhip_455MeVCM8_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 450-460MeV Photon Energy -ve Hel (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_455MeVCM8_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[7][3] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[7][3] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_465MeVCM8_PosHel = Phip_465MeVCM8_Para_PosHel->GetAsymmetry(Phip_465MeVCM8_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_465MeVCM8_PosHel->SetName("ParaPerpAsymmPhip465MeVCM8");
  ParaPerpAsymmPhip_465MeVCM8_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 460-470MeV Photon Energy -ve Hel (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_465MeVCM8_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[7][4] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[7][4] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_475MeVCM8_PosHel = Phip_475MeVCM8_Para_PosHel->GetAsymmetry(Phip_475MeVCM8_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_475MeVCM8_PosHel->SetName("ParaPerpAsymmPhip475MeVCM8");
  ParaPerpAsymmPhip_475MeVCM8_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 470-480MeV Photon Energy -ve Hel (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_475MeVCM8_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[7][5] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[7][5] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_485MeVCM8_PosHel = Phip_485MeVCM8_Para_PosHel->GetAsymmetry(Phip_485MeVCM8_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_485MeVCM8_PosHel->SetName("ParaPerpAsymmPhip485MeVCM8");
  ParaPerpAsymmPhip_485MeVCM8_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 480-490MeV Photon Energy -ve Hel (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_485MeVCM8_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[7][6] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[7][6] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_495MeVCM8_PosHel = Phip_495MeVCM8_Para_PosHel->GetAsymmetry(Phip_495MeVCM8_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_495MeVCM8_PosHel->SetName("ParaPerpAsymmPhip495MeVCM8");
  ParaPerpAsymmPhip_495MeVCM8_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 490-500MeV Photon Energy -ve Hel (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_495MeVCM8_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[7][7] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[7][7] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_505MeVCM8_PosHel = Phip_505MeVCM8_Para_PosHel->GetAsymmetry(Phip_505MeVCM8_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_505MeVCM8_PosHel->SetName("ParaPerpAsymmPhip505MeVCM8");
  ParaPerpAsymmPhip_505MeVCM8_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 500-510MeV Photon Energy -ve Hel (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_505MeVCM8_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[7][8] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[7][8] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_515MeVCM8_PosHel = Phip_515MeVCM8_Para_PosHel->GetAsymmetry(Phip_515MeVCM8_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_515MeVCM8_PosHel->SetName("ParaPerpAsymmPhip515MeVCM8");
  ParaPerpAsymmPhip_515MeVCM8_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 510-520MeV Photon Energy -ve Hel (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_515MeVCM8_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[7][9] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[7][9] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_525MeVCM8_PosHel = Phip_525MeVCM8_Para_PosHel->GetAsymmetry(Phip_525MeVCM8_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_525MeVCM8_PosHel->SetName("ParaPerpAsymmPhip525MeVCM8");
  ParaPerpAsymmPhip_525MeVCM8_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 520-530MeV Photon Energy -ve Hel (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_525MeVCM8_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[7][10] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[7][10] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_535MeVCM8_PosHel = Phip_535MeVCM8_Para_PosHel->GetAsymmetry(Phip_535MeVCM8_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_535MeVCM8_PosHel->SetName("ParaPerpAsymmPhip535MeVCM8");
  ParaPerpAsymmPhip_535MeVCM8_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 530-540MeV Photon Energy -ve Hel (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_535MeVCM8_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[7][11] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[7][11] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_545MeVCM8_PosHel = Phip_545MeVCM8_Para_PosHel->GetAsymmetry(Phip_545MeVCM8_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_545MeVCM8_PosHel->SetName("ParaPerpAsymmPhip545MeVCM8");
  ParaPerpAsymmPhip_545MeVCM8_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 540-550MeV Photon Energy -ve Hel (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_545MeVCM8_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[7][12] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[7][12] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_555MeVCM8_PosHel = Phip_555MeVCM8_Para_PosHel->GetAsymmetry(Phip_555MeVCM8_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_555MeVCM8_PosHel->SetName("ParaPerpAsymmPhip555MeVCM8");
  ParaPerpAsymmPhip_555MeVCM8_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 550-560MeV Photon Energy -ve Hel (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_555MeVCM8_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[7][13] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[7][13] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_565MeVCM8_PosHel = Phip_565MeVCM8_Para_PosHel->GetAsymmetry(Phip_565MeVCM8_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_565MeVCM8_PosHel->SetName("ParaPerpAsymmPhip565MeVCM8");
  ParaPerpAsymmPhip_565MeVCM8_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 560-570MeV Photon Energy -ve Hel (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_565MeVCM8_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[7][14] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[7][14] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_575MeVCM8_PosHel = Phip_575MeVCM8_Para_PosHel->GetAsymmetry(Phip_575MeVCM8_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_575MeVCM8_PosHel->SetName("ParaPerpAsymmPhip575MeVCM8");
  ParaPerpAsymmPhip_575MeVCM8_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 570-580MeV Photon Energy -ve Hel (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_575MeVCM8_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[7][15] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[7][15] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_585MeVCM8_PosHel = Phip_585MeVCM8_Para_PosHel->GetAsymmetry(Phip_585MeVCM8_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_585MeVCM8_PosHel->SetName("ParaPerpAsymmPhip585MeVCM8");
  ParaPerpAsymmPhip_585MeVCM8_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 580-590MeV Photon Energy -ve Hel (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_585MeVCM8_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[7][16] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[7][16] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_595MeVCM8_PosHel = Phip_595MeVCM8_Para_PosHel->GetAsymmetry(Phip_595MeVCM8_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_595MeVCM8_PosHel->SetName("ParaPerpAsymmPhip595MeVCM8");
  ParaPerpAsymmPhip_595MeVCM8_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 590-600MeV Photon Energy -ve Hel (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_595MeVCM8_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[7][17] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[7][17] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_605MeVCM8_PosHel = Phip_605MeVCM8_Para_PosHel->GetAsymmetry(Phip_605MeVCM8_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_605MeVCM8_PosHel->SetName("ParaPerpAsymmPhip605MeVCM8");
  ParaPerpAsymmPhip_605MeVCM8_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 600-610MeV Photon Energy -ve Hel (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_605MeVCM8_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[7][18] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[7][18] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_615MeVCM8_PosHel = Phip_615MeVCM8_Para_PosHel->GetAsymmetry(Phip_615MeVCM8_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_615MeVCM8_PosHel->SetName("ParaPerpAsymmPhip615MeVCM8");
  ParaPerpAsymmPhip_615MeVCM8_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 610-620MeV Photon Energy -ve Hel (CosTheta-0.4-(-0.6))");
  ParaPerpAsymmPhip_615MeVCM8_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[7][19] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[7][19] = CosFit->GetParError(0);

  ///////////////////////////////////////////
  //////////////////CM9 p ///////////////////
  ///////////////////////////////////////////

  ParaPerpAsymmPhip_425MeVCM9_PosHel = Phip_425MeVCM9_Para_PosHel->GetAsymmetry(Phip_425MeVCM9_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_425MeVCM9_PosHel->SetName("ParaPerpAsymmPhip425MeVCM9");
  ParaPerpAsymmPhip_425MeVCM9_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 420-430MeV Photon Energy -ve Hel (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_425MeVCM9_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[8][0] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[8][0] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_435MeVCM9_PosHel = Phip_435MeVCM9_Para_PosHel->GetAsymmetry(Phip_435MeVCM9_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_435MeVCM9_PosHel->SetName("ParaPerpAsymmPhip435MeVCM9");
  ParaPerpAsymmPhip_435MeVCM9_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 430-440MeV Photon Energy -ve Hel (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_435MeVCM9_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[8][1] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[8][1] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_445MeVCM9_PosHel = Phip_445MeVCM9_Para_PosHel->GetAsymmetry(Phip_445MeVCM9_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_445MeVCM9_PosHel->SetName("ParaPerpAsymmPhip445MeVCM9");
  ParaPerpAsymmPhip_445MeVCM9_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 440-450MeV Photon Energy -ve Hel (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_445MeVCM9_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[8][2] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[8][2] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_455MeVCM9_PosHel = Phip_455MeVCM9_Para_PosHel->GetAsymmetry(Phip_455MeVCM9_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_455MeVCM9_PosHel->SetName("ParaPerpAsymmPhip455MeVCM9");
  ParaPerpAsymmPhip_455MeVCM9_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 450-460MeV Photon Energy -ve Hel (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_455MeVCM9_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[8][3] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[8][3] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_465MeVCM9_PosHel = Phip_465MeVCM9_Para_PosHel->GetAsymmetry(Phip_465MeVCM9_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_465MeVCM9_PosHel->SetName("ParaPerpAsymmPhip465MeVCM9");
  ParaPerpAsymmPhip_465MeVCM9_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 460-470MeV Photon Energy -ve Hel (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_465MeVCM9_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[8][4] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[8][4] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_475MeVCM9_PosHel = Phip_475MeVCM9_Para_PosHel->GetAsymmetry(Phip_475MeVCM9_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_475MeVCM9_PosHel->SetName("ParaPerpAsymmPhip475MeVCM9");
  ParaPerpAsymmPhip_475MeVCM9_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 470-480MeV Photon Energy -ve Hel (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_475MeVCM9_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[8][5] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[8][5] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_485MeVCM9_PosHel = Phip_485MeVCM9_Para_PosHel->GetAsymmetry(Phip_485MeVCM9_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_485MeVCM9_PosHel->SetName("ParaPerpAsymmPhip485MeVCM9");
  ParaPerpAsymmPhip_485MeVCM9_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 480-490MeV Photon Energy -ve Hel (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_485MeVCM9_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[8][6] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[8][6] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_495MeVCM9_PosHel = Phip_495MeVCM9_Para_PosHel->GetAsymmetry(Phip_495MeVCM9_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_495MeVCM9_PosHel->SetName("ParaPerpAsymmPhip495MeVCM9");
  ParaPerpAsymmPhip_495MeVCM9_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 490-500MeV Photon Energy -ve Hel (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_495MeVCM9_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[8][7] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[8][7] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_505MeVCM9_PosHel = Phip_505MeVCM9_Para_PosHel->GetAsymmetry(Phip_505MeVCM9_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_505MeVCM9_PosHel->SetName("ParaPerpAsymmPhip505MeVCM9");
  ParaPerpAsymmPhip_505MeVCM9_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 500-510MeV Photon Energy -ve Hel (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_505MeVCM9_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[8][8] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[8][8] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_515MeVCM9_PosHel = Phip_515MeVCM9_Para_PosHel->GetAsymmetry(Phip_515MeVCM9_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_515MeVCM9_PosHel->SetName("ParaPerpAsymmPhip515MeVCM9");
  ParaPerpAsymmPhip_515MeVCM9_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 510-520MeV Photon Energy -ve Hel (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_515MeVCM9_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[8][9] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[8][9] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_525MeVCM9_PosHel = Phip_525MeVCM9_Para_PosHel->GetAsymmetry(Phip_525MeVCM9_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_525MeVCM9_PosHel->SetName("ParaPerpAsymmPhip525MeVCM9");
  ParaPerpAsymmPhip_525MeVCM9_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 520-530MeV Photon Energy -ve Hel (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_525MeVCM9_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[8][10] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[8][10] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_535MeVCM9_PosHel = Phip_535MeVCM9_Para_PosHel->GetAsymmetry(Phip_535MeVCM9_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_535MeVCM9_PosHel->SetName("ParaPerpAsymmPhip535MeVCM9");
  ParaPerpAsymmPhip_535MeVCM9_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 530-540MeV Photon Energy -ve Hel (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_535MeVCM9_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[8][11] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[8][11] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_545MeVCM9_PosHel = Phip_545MeVCM9_Para_PosHel->GetAsymmetry(Phip_545MeVCM9_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_545MeVCM9_PosHel->SetName("ParaPerpAsymmPhip545MeVCM9");
  ParaPerpAsymmPhip_545MeVCM9_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 540-550MeV Photon Energy -ve Hel (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_545MeVCM9_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[8][12] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[8][12] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_555MeVCM9_PosHel = Phip_555MeVCM9_Para_PosHel->GetAsymmetry(Phip_555MeVCM9_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_555MeVCM9_PosHel->SetName("ParaPerpAsymmPhip555MeVCM9");
  ParaPerpAsymmPhip_555MeVCM9_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 550-560MeV Photon Energy -ve Hel (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_555MeVCM9_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[8][13] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[8][13] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_565MeVCM9_PosHel = Phip_565MeVCM9_Para_PosHel->GetAsymmetry(Phip_565MeVCM9_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_565MeVCM9_PosHel->SetName("ParaPerpAsymmPhip565MeVCM9");
  ParaPerpAsymmPhip_565MeVCM9_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 560-570MeV Photon Energy -ve Hel (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_565MeVCM9_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[8][14] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[8][14] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_575MeVCM9_PosHel = Phip_575MeVCM9_Para_PosHel->GetAsymmetry(Phip_575MeVCM9_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_575MeVCM9_PosHel->SetName("ParaPerpAsymmPhip575MeVCM9");
  ParaPerpAsymmPhip_575MeVCM9_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 570-580MeV Photon Energy -ve Hel (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_575MeVCM9_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[8][15] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[8][15] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_585MeVCM9_PosHel = Phip_585MeVCM9_Para_PosHel->GetAsymmetry(Phip_585MeVCM9_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_585MeVCM9_PosHel->SetName("ParaPerpAsymmPhip585MeVCM9");
  ParaPerpAsymmPhip_585MeVCM9_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 580-590MeV Photon Energy -ve Hel (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_585MeVCM9_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[8][16] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[8][16] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_595MeVCM9_PosHel = Phip_595MeVCM9_Para_PosHel->GetAsymmetry(Phip_595MeVCM9_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_595MeVCM9_PosHel->SetName("ParaPerpAsymmPhip595MeVCM9");
  ParaPerpAsymmPhip_595MeVCM9_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 590-600MeV Photon Energy -ve Hel (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_595MeVCM9_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[8][17] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[8][17] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_605MeVCM9_PosHel = Phip_605MeVCM9_Para_PosHel->GetAsymmetry(Phip_605MeVCM9_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_605MeVCM9_PosHel->SetName("ParaPerpAsymmPhip605MeVCM9");
  ParaPerpAsymmPhip_605MeVCM9_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 600-610MeV Photon Energy -ve Hel (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_605MeVCM9_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[8][18] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[8][18] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_615MeVCM9_PosHel = Phip_615MeVCM9_Para_PosHel->GetAsymmetry(Phip_615MeVCM9_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_615MeVCM9_PosHel->SetName("ParaPerpAsymmPhip615MeVCM9");
  ParaPerpAsymmPhip_615MeVCM9_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 610-620MeV Photon Energy -ve Hel (CosTheta-0.6-(-0.8))");
  ParaPerpAsymmPhip_615MeVCM9_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[8][19] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[8][19] = CosFit->GetParError(0);

  ///////////////////////////////////////////
  //////////////////CM10 p ///////////////////
  ///////////////////////////////////////////

  ParaPerpAsymmPhip_425MeVCM10_PosHel = Phip_425MeVCM10_Para_PosHel->GetAsymmetry(Phip_425MeVCM10_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_425MeVCM10_PosHel->SetName("ParaPerpAsymmPhip425MeVCM10");
  ParaPerpAsymmPhip_425MeVCM10_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 420-430MeV Photon Energy -ve Hel (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_425MeVCM10_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[9][0] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[9][0] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_435MeVCM10_PosHel = Phip_435MeVCM10_Para_PosHel->GetAsymmetry(Phip_435MeVCM10_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_435MeVCM10_PosHel->SetName("ParaPerpAsymmPhip435MeVCM10");
  ParaPerpAsymmPhip_435MeVCM10_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 430-440MeV Photon Energy -ve Hel (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_435MeVCM10_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[9][1] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[9][1] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_445MeVCM10_PosHel = Phip_445MeVCM10_Para_PosHel->GetAsymmetry(Phip_445MeVCM10_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_445MeVCM10_PosHel->SetName("ParaPerpAsymmPhip445MeVCM10");
  ParaPerpAsymmPhip_445MeVCM10_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 440-450MeV Photon Energy -ve Hel (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_445MeVCM10_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[9][2] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[9][2] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_455MeVCM10_PosHel = Phip_455MeVCM10_Para_PosHel->GetAsymmetry(Phip_455MeVCM10_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_455MeVCM10_PosHel->SetName("ParaPerpAsymmPhip455MeVCM10");
  ParaPerpAsymmPhip_455MeVCM10_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 450-460MeV Photon Energy -ve Hel (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_455MeVCM10_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[9][3] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[9][3] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_465MeVCM10_PosHel = Phip_465MeVCM10_Para_PosHel->GetAsymmetry(Phip_465MeVCM10_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_465MeVCM10_PosHel->SetName("ParaPerpAsymmPhip465MeVCM10");
  ParaPerpAsymmPhip_465MeVCM10_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 460-470MeV Photon Energy -ve Hel (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_465MeVCM10_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[9][4] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[9][4] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_475MeVCM10_PosHel = Phip_475MeVCM10_Para_PosHel->GetAsymmetry(Phip_475MeVCM10_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_475MeVCM10_PosHel->SetName("ParaPerpAsymmPhip475MeVCM10");
  ParaPerpAsymmPhip_475MeVCM10_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 470-480MeV Photon Energy -ve Hel (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_475MeVCM10_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[9][5] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[9][5] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_485MeVCM10_PosHel = Phip_485MeVCM10_Para_PosHel->GetAsymmetry(Phip_485MeVCM10_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_485MeVCM10_PosHel->SetName("ParaPerpAsymmPhip485MeVCM10");
  ParaPerpAsymmPhip_485MeVCM10_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 480-490MeV Photon Energy -ve Hel (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_485MeVCM10_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[9][6] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[9][6] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_495MeVCM10_PosHel = Phip_495MeVCM10_Para_PosHel->GetAsymmetry(Phip_495MeVCM10_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_495MeVCM10_PosHel->SetName("ParaPerpAsymmPhip495MeVCM10");
  ParaPerpAsymmPhip_495MeVCM10_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 490-500MeV Photon Energy -ve Hel (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_495MeVCM10_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[9][7] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[9][7] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_505MeVCM10_PosHel = Phip_505MeVCM10_Para_PosHel->GetAsymmetry(Phip_505MeVCM10_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_505MeVCM10_PosHel->SetName("ParaPerpAsymmPhip505MeVCM10");
  ParaPerpAsymmPhip_505MeVCM10_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 500-510MeV Photon Energy -ve Hel (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_505MeVCM10_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[9][8] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[9][8] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_515MeVCM10_PosHel = Phip_515MeVCM10_Para_PosHel->GetAsymmetry(Phip_515MeVCM10_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_515MeVCM10_PosHel->SetName("ParaPerpAsymmPhip515MeVCM10");
  ParaPerpAsymmPhip_515MeVCM10_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 510-520MeV Photon Energy -ve Hel (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_515MeVCM10_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[9][9] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[9][9] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_525MeVCM10_PosHel = Phip_525MeVCM10_Para_PosHel->GetAsymmetry(Phip_525MeVCM10_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_525MeVCM10_PosHel->SetName("ParaPerpAsymmPhip525MeVCM10");
  ParaPerpAsymmPhip_525MeVCM10_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 520-530MeV Photon Energy -ve Hel (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_525MeVCM10_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[9][10] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[9][10] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_535MeVCM10_PosHel = Phip_535MeVCM10_Para_PosHel->GetAsymmetry(Phip_535MeVCM10_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_535MeVCM10_PosHel->SetName("ParaPerpAsymmPhip535MeVCM10");
  ParaPerpAsymmPhip_535MeVCM10_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 530-540MeV Photon Energy -ve Hel (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_535MeVCM10_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[9][11] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[9][11] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_545MeVCM10_PosHel = Phip_545MeVCM10_Para_PosHel->GetAsymmetry(Phip_545MeVCM10_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_545MeVCM10_PosHel->SetName("ParaPerpAsymmPhip545MeVCM10");
  ParaPerpAsymmPhip_545MeVCM10_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 540-550MeV Photon Energy -ve Hel (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_545MeVCM10_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[9][12] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[9][12] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_555MeVCM10_PosHel = Phip_555MeVCM10_Para_PosHel->GetAsymmetry(Phip_555MeVCM10_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_555MeVCM10_PosHel->SetName("ParaPerpAsymmPhip555MeVCM10");
  ParaPerpAsymmPhip_555MeVCM10_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 550-560MeV Photon Energy -ve Hel (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_555MeVCM10_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[9][13] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[9][13] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_565MeVCM10_PosHel = Phip_565MeVCM10_Para_PosHel->GetAsymmetry(Phip_565MeVCM10_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_565MeVCM10_PosHel->SetName("ParaPerpAsymmPhip565MeVCM10");
  ParaPerpAsymmPhip_565MeVCM10_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 560-570MeV Photon Energy -ve Hel (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_565MeVCM10_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[9][14] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[9][14] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_575MeVCM10_PosHel = Phip_575MeVCM10_Para_PosHel->GetAsymmetry(Phip_575MeVCM10_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_575MeVCM10_PosHel->SetName("ParaPerpAsymmPhip575MeVCM10");
  ParaPerpAsymmPhip_575MeVCM10_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 570-580MeV Photon Energy -ve Hel (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_575MeVCM10_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[9][15] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[9][15] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_585MeVCM10_PosHel = Phip_585MeVCM10_Para_PosHel->GetAsymmetry(Phip_585MeVCM10_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_585MeVCM10_PosHel->SetName("ParaPerpAsymmPhip585MeVCM10");
  ParaPerpAsymmPhip_585MeVCM10_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 580-590MeV Photon Energy -ve Hel (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_585MeVCM10_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[9][16] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[9][16] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_595MeVCM10_PosHel = Phip_595MeVCM10_Para_PosHel->GetAsymmetry(Phip_595MeVCM10_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_595MeVCM10_PosHel->SetName("ParaPerpAsymmPhip595MeVCM10");
  ParaPerpAsymmPhip_595MeVCM10_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 590-600MeV Photon Energy -ve Hel (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_595MeVCM10_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[9][17] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[9][17] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_605MeVCM10_PosHel = Phip_605MeVCM10_Para_PosHel->GetAsymmetry(Phip_605MeVCM10_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_605MeVCM10_PosHel->SetName("ParaPerpAsymmPhip605MeVCM10");
  ParaPerpAsymmPhip_605MeVCM10_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 600-610MeV Photon Energy -ve Hel (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_605MeVCM10_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[9][18] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[9][18] = CosFit->GetParError(0);

  ParaPerpAsymmPhip_615MeVCM10_PosHel = Phip_615MeVCM10_Para_PosHel->GetAsymmetry(Phip_615MeVCM10_Perp_PosHel, ScaleFactor, ScaleFactorErr);
  ParaPerpAsymmPhip_615MeVCM10_PosHel->SetName("ParaPerpAsymmPhip615MeVCM10");
  ParaPerpAsymmPhip_615MeVCM10_PosHel->SetTitle("Proton Para/Perp Phi Asymmetry for 610-620MeV Photon Energy -ve Hel (CosTheta-0.8-(-1.0))");
  ParaPerpAsymmPhip_615MeVCM10_PosHel->Fit("CosFit", "Q");
  pCosAmpPosHel[9][19] = CosFit->GetParameter(0);
  pCosAmpErrPosHel[9][19] = CosFit->GetParError(0);

  TFile f1("ParaPerpAsymm_NoScatt_Total_10.root", "RECREATE");

  ParaPerpAsymmPhip_425MeVCM1_NegHel->Write();
  ParaPerpAsymmPhip_435MeVCM1_NegHel->Write();
  ParaPerpAsymmPhip_445MeVCM1_NegHel->Write();
  ParaPerpAsymmPhip_455MeVCM1_NegHel->Write();
  ParaPerpAsymmPhip_465MeVCM1_NegHel->Write();
  ParaPerpAsymmPhip_475MeVCM1_NegHel->Write();
  ParaPerpAsymmPhip_485MeVCM1_NegHel->Write();
  ParaPerpAsymmPhip_495MeVCM1_NegHel->Write();
  ParaPerpAsymmPhip_505MeVCM1_NegHel->Write();
  ParaPerpAsymmPhip_515MeVCM1_NegHel->Write();
  ParaPerpAsymmPhip_525MeVCM1_NegHel->Write();
  ParaPerpAsymmPhip_535MeVCM1_NegHel->Write();
  ParaPerpAsymmPhip_545MeVCM1_NegHel->Write();
  ParaPerpAsymmPhip_555MeVCM1_NegHel->Write();
  ParaPerpAsymmPhip_565MeVCM1_NegHel->Write();
  ParaPerpAsymmPhip_575MeVCM1_NegHel->Write();
  ParaPerpAsymmPhip_585MeVCM1_NegHel->Write();
  ParaPerpAsymmPhip_595MeVCM1_NegHel->Write();
  ParaPerpAsymmPhip_605MeVCM1_NegHel->Write();
  ParaPerpAsymmPhip_615MeVCM1_NegHel->Write();

  ParaPerpAsymmPhip_425MeVCM2_NegHel->Write();
  ParaPerpAsymmPhip_435MeVCM2_NegHel->Write();
  ParaPerpAsymmPhip_445MeVCM2_NegHel->Write();
  ParaPerpAsymmPhip_455MeVCM2_NegHel->Write();
  ParaPerpAsymmPhip_465MeVCM2_NegHel->Write();
  ParaPerpAsymmPhip_475MeVCM2_NegHel->Write();
  ParaPerpAsymmPhip_485MeVCM2_NegHel->Write();
  ParaPerpAsymmPhip_495MeVCM2_NegHel->Write();
  ParaPerpAsymmPhip_505MeVCM2_NegHel->Write();
  ParaPerpAsymmPhip_515MeVCM2_NegHel->Write();
  ParaPerpAsymmPhip_525MeVCM2_NegHel->Write();
  ParaPerpAsymmPhip_535MeVCM2_NegHel->Write();
  ParaPerpAsymmPhip_545MeVCM2_NegHel->Write();
  ParaPerpAsymmPhip_555MeVCM2_NegHel->Write();
  ParaPerpAsymmPhip_565MeVCM2_NegHel->Write();
  ParaPerpAsymmPhip_575MeVCM2_NegHel->Write();
  ParaPerpAsymmPhip_585MeVCM2_NegHel->Write();
  ParaPerpAsymmPhip_595MeVCM2_NegHel->Write();
  ParaPerpAsymmPhip_605MeVCM2_NegHel->Write();
  ParaPerpAsymmPhip_615MeVCM2_NegHel->Write();

  ParaPerpAsymmPhip_425MeVCM3_NegHel->Write();
  ParaPerpAsymmPhip_435MeVCM3_NegHel->Write();
  ParaPerpAsymmPhip_445MeVCM3_NegHel->Write();
  ParaPerpAsymmPhip_455MeVCM3_NegHel->Write();
  ParaPerpAsymmPhip_465MeVCM3_NegHel->Write();
  ParaPerpAsymmPhip_475MeVCM3_NegHel->Write();
  ParaPerpAsymmPhip_485MeVCM3_NegHel->Write();
  ParaPerpAsymmPhip_495MeVCM3_NegHel->Write();
  ParaPerpAsymmPhip_505MeVCM3_NegHel->Write();
  ParaPerpAsymmPhip_515MeVCM3_NegHel->Write();
  ParaPerpAsymmPhip_525MeVCM3_NegHel->Write();
  ParaPerpAsymmPhip_535MeVCM3_NegHel->Write();
  ParaPerpAsymmPhip_545MeVCM3_NegHel->Write();
  ParaPerpAsymmPhip_555MeVCM3_NegHel->Write();
  ParaPerpAsymmPhip_565MeVCM3_NegHel->Write();
  ParaPerpAsymmPhip_575MeVCM3_NegHel->Write();
  ParaPerpAsymmPhip_585MeVCM3_NegHel->Write();
  ParaPerpAsymmPhip_595MeVCM3_NegHel->Write();
  ParaPerpAsymmPhip_605MeVCM3_NegHel->Write();
  ParaPerpAsymmPhip_615MeVCM3_NegHel->Write();

  ParaPerpAsymmPhip_425MeVCM4_NegHel->Write();
  ParaPerpAsymmPhip_435MeVCM4_NegHel->Write();
  ParaPerpAsymmPhip_445MeVCM4_NegHel->Write();
  ParaPerpAsymmPhip_455MeVCM4_NegHel->Write();
  ParaPerpAsymmPhip_465MeVCM4_NegHel->Write();
  ParaPerpAsymmPhip_475MeVCM4_NegHel->Write();
  ParaPerpAsymmPhip_485MeVCM4_NegHel->Write();
  ParaPerpAsymmPhip_495MeVCM4_NegHel->Write();
  ParaPerpAsymmPhip_505MeVCM4_NegHel->Write();
  ParaPerpAsymmPhip_515MeVCM4_NegHel->Write();
  ParaPerpAsymmPhip_525MeVCM4_NegHel->Write();
  ParaPerpAsymmPhip_535MeVCM4_NegHel->Write();
  ParaPerpAsymmPhip_545MeVCM4_NegHel->Write();
  ParaPerpAsymmPhip_555MeVCM4_NegHel->Write();
  ParaPerpAsymmPhip_565MeVCM4_NegHel->Write();
  ParaPerpAsymmPhip_575MeVCM4_NegHel->Write();
  ParaPerpAsymmPhip_585MeVCM4_NegHel->Write();
  ParaPerpAsymmPhip_595MeVCM4_NegHel->Write();
  ParaPerpAsymmPhip_605MeVCM4_NegHel->Write();
  ParaPerpAsymmPhip_615MeVCM4_NegHel->Write();

  ParaPerpAsymmPhip_425MeVCM5_NegHel->Write();
  ParaPerpAsymmPhip_435MeVCM5_NegHel->Write();
  ParaPerpAsymmPhip_445MeVCM5_NegHel->Write();
  ParaPerpAsymmPhip_455MeVCM5_NegHel->Write();
  ParaPerpAsymmPhip_465MeVCM5_NegHel->Write();
  ParaPerpAsymmPhip_475MeVCM5_NegHel->Write();
  ParaPerpAsymmPhip_485MeVCM5_NegHel->Write();
  ParaPerpAsymmPhip_495MeVCM5_NegHel->Write();
  ParaPerpAsymmPhip_505MeVCM5_NegHel->Write();
  ParaPerpAsymmPhip_515MeVCM5_NegHel->Write();
  ParaPerpAsymmPhip_525MeVCM5_NegHel->Write();
  ParaPerpAsymmPhip_535MeVCM5_NegHel->Write();
  ParaPerpAsymmPhip_545MeVCM5_NegHel->Write();
  ParaPerpAsymmPhip_555MeVCM5_NegHel->Write();
  ParaPerpAsymmPhip_565MeVCM5_NegHel->Write();
  ParaPerpAsymmPhip_575MeVCM5_NegHel->Write();
  ParaPerpAsymmPhip_585MeVCM5_NegHel->Write();
  ParaPerpAsymmPhip_595MeVCM5_NegHel->Write();
  ParaPerpAsymmPhip_605MeVCM5_NegHel->Write();
  ParaPerpAsymmPhip_615MeVCM5_NegHel->Write();

  ParaPerpAsymmPhip_425MeVCM6_NegHel->Write();
  ParaPerpAsymmPhip_435MeVCM6_NegHel->Write();
  ParaPerpAsymmPhip_445MeVCM6_NegHel->Write();
  ParaPerpAsymmPhip_455MeVCM6_NegHel->Write();
  ParaPerpAsymmPhip_465MeVCM6_NegHel->Write();
  ParaPerpAsymmPhip_475MeVCM6_NegHel->Write();
  ParaPerpAsymmPhip_485MeVCM6_NegHel->Write();
  ParaPerpAsymmPhip_495MeVCM6_NegHel->Write();
  ParaPerpAsymmPhip_505MeVCM6_NegHel->Write();
  ParaPerpAsymmPhip_515MeVCM6_NegHel->Write();
  ParaPerpAsymmPhip_525MeVCM6_NegHel->Write();
  ParaPerpAsymmPhip_535MeVCM6_NegHel->Write();
  ParaPerpAsymmPhip_545MeVCM6_NegHel->Write();
  ParaPerpAsymmPhip_555MeVCM6_NegHel->Write();
  ParaPerpAsymmPhip_565MeVCM6_NegHel->Write();
  ParaPerpAsymmPhip_575MeVCM6_NegHel->Write();
  ParaPerpAsymmPhip_585MeVCM6_NegHel->Write();
  ParaPerpAsymmPhip_595MeVCM6_NegHel->Write();
  ParaPerpAsymmPhip_605MeVCM6_NegHel->Write();
  ParaPerpAsymmPhip_615MeVCM6_NegHel->Write();

  ParaPerpAsymmPhip_425MeVCM7_NegHel->Write();
  ParaPerpAsymmPhip_435MeVCM7_NegHel->Write();
  ParaPerpAsymmPhip_445MeVCM7_NegHel->Write();
  ParaPerpAsymmPhip_455MeVCM7_NegHel->Write();
  ParaPerpAsymmPhip_465MeVCM7_NegHel->Write();
  ParaPerpAsymmPhip_475MeVCM7_NegHel->Write();
  ParaPerpAsymmPhip_485MeVCM7_NegHel->Write();
  ParaPerpAsymmPhip_495MeVCM7_NegHel->Write();
  ParaPerpAsymmPhip_505MeVCM7_NegHel->Write();
  ParaPerpAsymmPhip_515MeVCM7_NegHel->Write();
  ParaPerpAsymmPhip_525MeVCM7_NegHel->Write();
  ParaPerpAsymmPhip_535MeVCM7_NegHel->Write();
  ParaPerpAsymmPhip_545MeVCM7_NegHel->Write();
  ParaPerpAsymmPhip_555MeVCM7_NegHel->Write();
  ParaPerpAsymmPhip_565MeVCM7_NegHel->Write();
  ParaPerpAsymmPhip_575MeVCM7_NegHel->Write();
  ParaPerpAsymmPhip_585MeVCM7_NegHel->Write();
  ParaPerpAsymmPhip_595MeVCM7_NegHel->Write();
  ParaPerpAsymmPhip_605MeVCM7_NegHel->Write();
  ParaPerpAsymmPhip_615MeVCM7_NegHel->Write();

  ParaPerpAsymmPhip_425MeVCM8_NegHel->Write();
  ParaPerpAsymmPhip_435MeVCM8_NegHel->Write();
  ParaPerpAsymmPhip_445MeVCM8_NegHel->Write();
  ParaPerpAsymmPhip_455MeVCM8_NegHel->Write();
  ParaPerpAsymmPhip_465MeVCM8_NegHel->Write();
  ParaPerpAsymmPhip_475MeVCM8_NegHel->Write();
  ParaPerpAsymmPhip_485MeVCM8_NegHel->Write();
  ParaPerpAsymmPhip_495MeVCM8_NegHel->Write();
  ParaPerpAsymmPhip_505MeVCM8_NegHel->Write();
  ParaPerpAsymmPhip_515MeVCM8_NegHel->Write();
  ParaPerpAsymmPhip_525MeVCM8_NegHel->Write();
  ParaPerpAsymmPhip_535MeVCM8_NegHel->Write();
  ParaPerpAsymmPhip_545MeVCM8_NegHel->Write();
  ParaPerpAsymmPhip_555MeVCM8_NegHel->Write();
  ParaPerpAsymmPhip_565MeVCM8_NegHel->Write();
  ParaPerpAsymmPhip_575MeVCM8_NegHel->Write();
  ParaPerpAsymmPhip_585MeVCM8_NegHel->Write();
  ParaPerpAsymmPhip_595MeVCM8_NegHel->Write();
  ParaPerpAsymmPhip_605MeVCM8_NegHel->Write();
  ParaPerpAsymmPhip_615MeVCM8_NegHel->Write();

  ParaPerpAsymmPhip_425MeVCM9_NegHel->Write();
  ParaPerpAsymmPhip_435MeVCM9_NegHel->Write();
  ParaPerpAsymmPhip_445MeVCM9_NegHel->Write();
  ParaPerpAsymmPhip_455MeVCM9_NegHel->Write();
  ParaPerpAsymmPhip_465MeVCM9_NegHel->Write();
  ParaPerpAsymmPhip_475MeVCM9_NegHel->Write();
  ParaPerpAsymmPhip_485MeVCM9_NegHel->Write();
  ParaPerpAsymmPhip_495MeVCM9_NegHel->Write();
  ParaPerpAsymmPhip_505MeVCM9_NegHel->Write();
  ParaPerpAsymmPhip_515MeVCM9_NegHel->Write();
  ParaPerpAsymmPhip_525MeVCM9_NegHel->Write();
  ParaPerpAsymmPhip_535MeVCM9_NegHel->Write();
  ParaPerpAsymmPhip_545MeVCM9_NegHel->Write();
  ParaPerpAsymmPhip_555MeVCM9_NegHel->Write();
  ParaPerpAsymmPhip_565MeVCM9_NegHel->Write();
  ParaPerpAsymmPhip_575MeVCM9_NegHel->Write();
  ParaPerpAsymmPhip_585MeVCM9_NegHel->Write();
  ParaPerpAsymmPhip_595MeVCM9_NegHel->Write();
  ParaPerpAsymmPhip_605MeVCM9_NegHel->Write();
  ParaPerpAsymmPhip_615MeVCM9_NegHel->Write();

  ParaPerpAsymmPhip_425MeVCM10_NegHel->Write();
  ParaPerpAsymmPhip_435MeVCM10_NegHel->Write();
  ParaPerpAsymmPhip_445MeVCM10_NegHel->Write();
  ParaPerpAsymmPhip_455MeVCM10_NegHel->Write();
  ParaPerpAsymmPhip_465MeVCM10_NegHel->Write();
  ParaPerpAsymmPhip_475MeVCM10_NegHel->Write();
  ParaPerpAsymmPhip_485MeVCM10_NegHel->Write();
  ParaPerpAsymmPhip_495MeVCM10_NegHel->Write();
  ParaPerpAsymmPhip_505MeVCM10_NegHel->Write();
  ParaPerpAsymmPhip_515MeVCM10_NegHel->Write();
  ParaPerpAsymmPhip_525MeVCM10_NegHel->Write();
  ParaPerpAsymmPhip_535MeVCM10_NegHel->Write();
  ParaPerpAsymmPhip_545MeVCM10_NegHel->Write();
  ParaPerpAsymmPhip_555MeVCM10_NegHel->Write();
  ParaPerpAsymmPhip_565MeVCM10_NegHel->Write();
  ParaPerpAsymmPhip_575MeVCM10_NegHel->Write();
  ParaPerpAsymmPhip_585MeVCM10_NegHel->Write();
  ParaPerpAsymmPhip_595MeVCM10_NegHel->Write();
  ParaPerpAsymmPhip_605MeVCM10_NegHel->Write();
  ParaPerpAsymmPhip_615MeVCM10_NegHel->Write();

  ParaPerpAsymmPhip_425MeVCM1_PosHel->Write();
  ParaPerpAsymmPhip_435MeVCM1_PosHel->Write();
  ParaPerpAsymmPhip_445MeVCM1_PosHel->Write();
  ParaPerpAsymmPhip_455MeVCM1_PosHel->Write();
  ParaPerpAsymmPhip_465MeVCM1_PosHel->Write();
  ParaPerpAsymmPhip_475MeVCM1_PosHel->Write();
  ParaPerpAsymmPhip_485MeVCM1_PosHel->Write();
  ParaPerpAsymmPhip_495MeVCM1_PosHel->Write();
  ParaPerpAsymmPhip_505MeVCM1_PosHel->Write();
  ParaPerpAsymmPhip_515MeVCM1_PosHel->Write();
  ParaPerpAsymmPhip_525MeVCM1_PosHel->Write();
  ParaPerpAsymmPhip_535MeVCM1_PosHel->Write();
  ParaPerpAsymmPhip_545MeVCM1_PosHel->Write();
  ParaPerpAsymmPhip_555MeVCM1_PosHel->Write();
  ParaPerpAsymmPhip_565MeVCM1_PosHel->Write();
  ParaPerpAsymmPhip_575MeVCM1_PosHel->Write();
  ParaPerpAsymmPhip_585MeVCM1_PosHel->Write();
  ParaPerpAsymmPhip_595MeVCM1_PosHel->Write();
  ParaPerpAsymmPhip_605MeVCM1_PosHel->Write();
  ParaPerpAsymmPhip_615MeVCM1_PosHel->Write();

  ParaPerpAsymmPhip_425MeVCM2_PosHel->Write();
  ParaPerpAsymmPhip_435MeVCM2_PosHel->Write();
  ParaPerpAsymmPhip_445MeVCM2_PosHel->Write();
  ParaPerpAsymmPhip_455MeVCM2_PosHel->Write();
  ParaPerpAsymmPhip_465MeVCM2_PosHel->Write();
  ParaPerpAsymmPhip_475MeVCM2_PosHel->Write();
  ParaPerpAsymmPhip_485MeVCM2_PosHel->Write();
  ParaPerpAsymmPhip_495MeVCM2_PosHel->Write();
  ParaPerpAsymmPhip_505MeVCM2_PosHel->Write();
  ParaPerpAsymmPhip_515MeVCM2_PosHel->Write();
  ParaPerpAsymmPhip_525MeVCM2_PosHel->Write();
  ParaPerpAsymmPhip_535MeVCM2_PosHel->Write();
  ParaPerpAsymmPhip_545MeVCM2_PosHel->Write();
  ParaPerpAsymmPhip_555MeVCM2_PosHel->Write();
  ParaPerpAsymmPhip_565MeVCM2_PosHel->Write();
  ParaPerpAsymmPhip_575MeVCM2_PosHel->Write();
  ParaPerpAsymmPhip_585MeVCM2_PosHel->Write();
  ParaPerpAsymmPhip_595MeVCM2_PosHel->Write();
  ParaPerpAsymmPhip_605MeVCM2_PosHel->Write();
  ParaPerpAsymmPhip_615MeVCM2_PosHel->Write();

  ParaPerpAsymmPhip_425MeVCM3_PosHel->Write();
  ParaPerpAsymmPhip_435MeVCM3_PosHel->Write();
  ParaPerpAsymmPhip_445MeVCM3_PosHel->Write();
  ParaPerpAsymmPhip_455MeVCM3_PosHel->Write();
  ParaPerpAsymmPhip_465MeVCM3_PosHel->Write();
  ParaPerpAsymmPhip_475MeVCM3_PosHel->Write();
  ParaPerpAsymmPhip_485MeVCM3_PosHel->Write();
  ParaPerpAsymmPhip_495MeVCM3_PosHel->Write();
  ParaPerpAsymmPhip_505MeVCM3_PosHel->Write();
  ParaPerpAsymmPhip_515MeVCM3_PosHel->Write();
  ParaPerpAsymmPhip_525MeVCM3_PosHel->Write();
  ParaPerpAsymmPhip_535MeVCM3_PosHel->Write();
  ParaPerpAsymmPhip_545MeVCM3_PosHel->Write();
  ParaPerpAsymmPhip_555MeVCM3_PosHel->Write();
  ParaPerpAsymmPhip_565MeVCM3_PosHel->Write();
  ParaPerpAsymmPhip_575MeVCM3_PosHel->Write();
  ParaPerpAsymmPhip_585MeVCM3_PosHel->Write();
  ParaPerpAsymmPhip_595MeVCM3_PosHel->Write();
  ParaPerpAsymmPhip_605MeVCM3_PosHel->Write();
  ParaPerpAsymmPhip_615MeVCM3_PosHel->Write();

  ParaPerpAsymmPhip_425MeVCM4_PosHel->Write();
  ParaPerpAsymmPhip_435MeVCM4_PosHel->Write();
  ParaPerpAsymmPhip_445MeVCM4_PosHel->Write();
  ParaPerpAsymmPhip_455MeVCM4_PosHel->Write();
  ParaPerpAsymmPhip_465MeVCM4_PosHel->Write();
  ParaPerpAsymmPhip_475MeVCM4_PosHel->Write();
  ParaPerpAsymmPhip_485MeVCM4_PosHel->Write();
  ParaPerpAsymmPhip_495MeVCM4_PosHel->Write();
  ParaPerpAsymmPhip_505MeVCM4_PosHel->Write();
  ParaPerpAsymmPhip_515MeVCM4_PosHel->Write();
  ParaPerpAsymmPhip_525MeVCM4_PosHel->Write();
  ParaPerpAsymmPhip_535MeVCM4_PosHel->Write();
  ParaPerpAsymmPhip_545MeVCM4_PosHel->Write();
  ParaPerpAsymmPhip_555MeVCM4_PosHel->Write();
  ParaPerpAsymmPhip_565MeVCM4_PosHel->Write();
  ParaPerpAsymmPhip_575MeVCM4_PosHel->Write();
  ParaPerpAsymmPhip_585MeVCM4_PosHel->Write();
  ParaPerpAsymmPhip_595MeVCM4_PosHel->Write();
  ParaPerpAsymmPhip_605MeVCM4_PosHel->Write();
  ParaPerpAsymmPhip_615MeVCM4_PosHel->Write();

  ParaPerpAsymmPhip_425MeVCM5_PosHel->Write();
  ParaPerpAsymmPhip_435MeVCM5_PosHel->Write();
  ParaPerpAsymmPhip_445MeVCM5_PosHel->Write();
  ParaPerpAsymmPhip_455MeVCM5_PosHel->Write();
  ParaPerpAsymmPhip_465MeVCM5_PosHel->Write();
  ParaPerpAsymmPhip_475MeVCM5_PosHel->Write();
  ParaPerpAsymmPhip_485MeVCM5_PosHel->Write();
  ParaPerpAsymmPhip_495MeVCM5_PosHel->Write();
  ParaPerpAsymmPhip_505MeVCM5_PosHel->Write();
  ParaPerpAsymmPhip_515MeVCM5_PosHel->Write();
  ParaPerpAsymmPhip_525MeVCM5_PosHel->Write();
  ParaPerpAsymmPhip_535MeVCM5_PosHel->Write();
  ParaPerpAsymmPhip_545MeVCM5_PosHel->Write();
  ParaPerpAsymmPhip_555MeVCM5_PosHel->Write();
  ParaPerpAsymmPhip_565MeVCM5_PosHel->Write();
  ParaPerpAsymmPhip_575MeVCM5_PosHel->Write();
  ParaPerpAsymmPhip_585MeVCM5_PosHel->Write();
  ParaPerpAsymmPhip_595MeVCM5_PosHel->Write();
  ParaPerpAsymmPhip_605MeVCM5_PosHel->Write();
  ParaPerpAsymmPhip_615MeVCM5_PosHel->Write();

  ParaPerpAsymmPhip_425MeVCM6_PosHel->Write();
  ParaPerpAsymmPhip_435MeVCM6_PosHel->Write();
  ParaPerpAsymmPhip_445MeVCM6_PosHel->Write();
  ParaPerpAsymmPhip_455MeVCM6_PosHel->Write();
  ParaPerpAsymmPhip_465MeVCM6_PosHel->Write();
  ParaPerpAsymmPhip_475MeVCM6_PosHel->Write();
  ParaPerpAsymmPhip_485MeVCM6_PosHel->Write();
  ParaPerpAsymmPhip_495MeVCM6_PosHel->Write();
  ParaPerpAsymmPhip_505MeVCM6_PosHel->Write();
  ParaPerpAsymmPhip_515MeVCM6_PosHel->Write();
  ParaPerpAsymmPhip_525MeVCM6_PosHel->Write();
  ParaPerpAsymmPhip_535MeVCM6_PosHel->Write();
  ParaPerpAsymmPhip_545MeVCM6_PosHel->Write();
  ParaPerpAsymmPhip_555MeVCM6_PosHel->Write();
  ParaPerpAsymmPhip_565MeVCM6_PosHel->Write();
  ParaPerpAsymmPhip_575MeVCM6_PosHel->Write();
  ParaPerpAsymmPhip_585MeVCM6_PosHel->Write();
  ParaPerpAsymmPhip_595MeVCM6_PosHel->Write();
  ParaPerpAsymmPhip_605MeVCM6_PosHel->Write();
  ParaPerpAsymmPhip_615MeVCM6_PosHel->Write();

  ParaPerpAsymmPhip_425MeVCM7_PosHel->Write();
  ParaPerpAsymmPhip_435MeVCM7_PosHel->Write();
  ParaPerpAsymmPhip_445MeVCM7_PosHel->Write();
  ParaPerpAsymmPhip_455MeVCM7_PosHel->Write();
  ParaPerpAsymmPhip_465MeVCM7_PosHel->Write();
  ParaPerpAsymmPhip_475MeVCM7_PosHel->Write();
  ParaPerpAsymmPhip_485MeVCM7_PosHel->Write();
  ParaPerpAsymmPhip_495MeVCM7_PosHel->Write();
  ParaPerpAsymmPhip_505MeVCM7_PosHel->Write();
  ParaPerpAsymmPhip_515MeVCM7_PosHel->Write();
  ParaPerpAsymmPhip_525MeVCM7_PosHel->Write();
  ParaPerpAsymmPhip_535MeVCM7_PosHel->Write();
  ParaPerpAsymmPhip_545MeVCM7_PosHel->Write();
  ParaPerpAsymmPhip_555MeVCM7_PosHel->Write();
  ParaPerpAsymmPhip_565MeVCM7_PosHel->Write();
  ParaPerpAsymmPhip_575MeVCM7_PosHel->Write();
  ParaPerpAsymmPhip_585MeVCM7_PosHel->Write();
  ParaPerpAsymmPhip_595MeVCM7_PosHel->Write();
  ParaPerpAsymmPhip_605MeVCM7_PosHel->Write();
  ParaPerpAsymmPhip_615MeVCM7_PosHel->Write();

  ParaPerpAsymmPhip_425MeVCM8_PosHel->Write();
  ParaPerpAsymmPhip_435MeVCM8_PosHel->Write();
  ParaPerpAsymmPhip_445MeVCM8_PosHel->Write();
  ParaPerpAsymmPhip_455MeVCM8_PosHel->Write();
  ParaPerpAsymmPhip_465MeVCM8_PosHel->Write();
  ParaPerpAsymmPhip_475MeVCM8_PosHel->Write();
  ParaPerpAsymmPhip_485MeVCM8_PosHel->Write();
  ParaPerpAsymmPhip_495MeVCM8_PosHel->Write();
  ParaPerpAsymmPhip_505MeVCM8_PosHel->Write();
  ParaPerpAsymmPhip_515MeVCM8_PosHel->Write();
  ParaPerpAsymmPhip_525MeVCM8_PosHel->Write();
  ParaPerpAsymmPhip_535MeVCM8_PosHel->Write();
  ParaPerpAsymmPhip_545MeVCM8_PosHel->Write();
  ParaPerpAsymmPhip_555MeVCM8_PosHel->Write();
  ParaPerpAsymmPhip_565MeVCM8_PosHel->Write();
  ParaPerpAsymmPhip_575MeVCM8_PosHel->Write();
  ParaPerpAsymmPhip_585MeVCM8_PosHel->Write();
  ParaPerpAsymmPhip_595MeVCM8_PosHel->Write();
  ParaPerpAsymmPhip_605MeVCM8_PosHel->Write();
  ParaPerpAsymmPhip_615MeVCM8_PosHel->Write();

  ParaPerpAsymmPhip_425MeVCM9_PosHel->Write();
  ParaPerpAsymmPhip_435MeVCM9_PosHel->Write();
  ParaPerpAsymmPhip_445MeVCM9_PosHel->Write();
  ParaPerpAsymmPhip_455MeVCM9_PosHel->Write();
  ParaPerpAsymmPhip_465MeVCM9_PosHel->Write();
  ParaPerpAsymmPhip_475MeVCM9_PosHel->Write();
  ParaPerpAsymmPhip_485MeVCM9_PosHel->Write();
  ParaPerpAsymmPhip_495MeVCM9_PosHel->Write();
  ParaPerpAsymmPhip_505MeVCM9_PosHel->Write();
  ParaPerpAsymmPhip_515MeVCM9_PosHel->Write();
  ParaPerpAsymmPhip_525MeVCM9_PosHel->Write();
  ParaPerpAsymmPhip_535MeVCM9_PosHel->Write();
  ParaPerpAsymmPhip_545MeVCM9_PosHel->Write();
  ParaPerpAsymmPhip_555MeVCM9_PosHel->Write();
  ParaPerpAsymmPhip_565MeVCM9_PosHel->Write();
  ParaPerpAsymmPhip_575MeVCM9_PosHel->Write();
  ParaPerpAsymmPhip_585MeVCM9_PosHel->Write();
  ParaPerpAsymmPhip_595MeVCM9_PosHel->Write();
  ParaPerpAsymmPhip_605MeVCM9_PosHel->Write();
  ParaPerpAsymmPhip_615MeVCM9_PosHel->Write();

  ParaPerpAsymmPhip_425MeVCM10_PosHel->Write();
  ParaPerpAsymmPhip_435MeVCM10_PosHel->Write();
  ParaPerpAsymmPhip_445MeVCM10_PosHel->Write();
  ParaPerpAsymmPhip_455MeVCM10_PosHel->Write();
  ParaPerpAsymmPhip_465MeVCM10_PosHel->Write();
  ParaPerpAsymmPhip_475MeVCM10_PosHel->Write();
  ParaPerpAsymmPhip_485MeVCM10_PosHel->Write();
  ParaPerpAsymmPhip_495MeVCM10_PosHel->Write();
  ParaPerpAsymmPhip_505MeVCM10_PosHel->Write();
  ParaPerpAsymmPhip_515MeVCM10_PosHel->Write();
  ParaPerpAsymmPhip_525MeVCM10_PosHel->Write();
  ParaPerpAsymmPhip_535MeVCM10_PosHel->Write();
  ParaPerpAsymmPhip_545MeVCM10_PosHel->Write();
  ParaPerpAsymmPhip_555MeVCM10_PosHel->Write();
  ParaPerpAsymmPhip_565MeVCM10_PosHel->Write();
  ParaPerpAsymmPhip_575MeVCM10_PosHel->Write();
  ParaPerpAsymmPhip_585MeVCM10_PosHel->Write();
  ParaPerpAsymmPhip_595MeVCM10_PosHel->Write();
  ParaPerpAsymmPhip_605MeVCM10_PosHel->Write();
  ParaPerpAsymmPhip_615MeVCM10_PosHel->Write();

  //Define new tree to store parameters in
  TTree* tree = new TTree("Parameter_Values", "Tree_of_Values");

  // Define branches to store parameters, (Branch Name, Variable, Type of Variable)
  tree->Branch("pCosAmp425NegHel", &pCosA425NegHel, "pCosA425NegHel/D");
  tree->Branch("pCosAmpErr425NegHel", &pCosAErr425NegHel, "pCosAErr425NegHel/D");
  tree->Branch("pCosAmp435NegHel", &pCosA435NegHel, "pCosA435NegHel/D");
  tree->Branch("pCosAmpErr435NegHel", &pCosAErr435NegHel, "pCosAErr435NegHel/D");
  tree->Branch("pCosAmp445NegHel", &pCosA445NegHel, "pCosA445NegHel/D");
  tree->Branch("pCosAmpErr445NegHel", &pCosAErr445NegHel, "pCosAErr445NegHel/D");
  tree->Branch("pCosAmp455NegHel", &pCosA455NegHel, "pCosA455NegHel/D");
  tree->Branch("pCosAmpErr455NegHel", &pCosAErr455NegHel, "pCosAErr455NegHel/D");
  tree->Branch("pCosAmp465NegHel", &pCosA465NegHel, "pCosA465NegHel/D");
  tree->Branch("pCosAmpErr465NegHel", &pCosAErr465NegHel, "pCosAErr465NegHel/D");
  tree->Branch("pCosAmp475NegHel", &pCosA475NegHel, "pCosA475NegHel/D");
  tree->Branch("pCosAmpErr475NegHel", &pCosAErr475NegHel, "pCosAErr475NegHel/D");
  tree->Branch("pCosAmp485NegHel", &pCosA485NegHel, "pCosA485NegHel/D");
  tree->Branch("pCosAmpErr485NegHel", &pCosAErr485NegHel, "pCosAErr485NegHel/D");
  tree->Branch("pCosAmp495NegHel", &pCosA495NegHel, "pCosA495NegHel/D");
  tree->Branch("pCosAmpErr495NegHel", &pCosAErr495NegHel, "pCosAErr495NegHel/D");
  tree->Branch("pCosAmp505NegHel", &pCosA505NegHel, "pCosA505NegHel/D");
  tree->Branch("pCosAmpErr505NegHel", &pCosAErr505NegHel, "pCosAErr505NegHel/D");
  tree->Branch("pCosAmp515NegHel", &pCosA515NegHel, "pCosA515NegHel/D");
  tree->Branch("pCosAmpErr515NegHel", &pCosAErr515NegHel, "pCosAErr515NegHel/D");
  tree->Branch("pCosAmp525NegHel", &pCosA525NegHel, "pCosA525NegHel/D");
  tree->Branch("pCosAmpErr525NegHel", &pCosAErr525NegHel, "pCosAErr525NegHel/D");
  tree->Branch("pCosAmp535NegHel", &pCosA535NegHel, "pCosA535NegHel/D");
  tree->Branch("pCosAmpErr535NegHel", &pCosAErr535NegHel, "pCosAErr535NegHel/D");
  tree->Branch("pCosAmp545NegHel", &pCosA545NegHel, "pCosA545NegHel/D");
  tree->Branch("pCosAmpErr545NegHel", &pCosAErr545NegHel, "pCosAErr545NegHel/D");
  tree->Branch("pCosAmp555NegHel", &pCosA555NegHel, "pCosA555NegHel/D");
  tree->Branch("pCosAmpErr555NegHel", &pCosAErr555NegHel, "pCosAErr555NegHel/D");
  tree->Branch("pCosAmp565NegHel", &pCosA565NegHel, "pCosA565NegHel/D");
  tree->Branch("pCosAmpErr565NegHel", &pCosAErr565NegHel, "pCosAErr565NegHel/D");
  tree->Branch("pCosAmp575NegHel", &pCosA575NegHel, "pCosA575NegHel/D");
  tree->Branch("pCosAmpErr575NegHel", &pCosAErr575NegHel, "pCosAErr575NegHel/D");
  tree->Branch("pCosAmp585NegHel", &pCosA585NegHel, "pCosA585NegHel/D");
  tree->Branch("pCosAmpErr585NegHel", &pCosAErr585NegHel, "pCosAErr585NegHel/D");
  tree->Branch("pCosAmp595NegHel", &pCosA595NegHel, "pCosA595NegHel/D");
  tree->Branch("pCosAmpErr595NegHel", &pCosAErr595NegHel, "pCosAErr595NegHel/D");
  tree->Branch("pCosAmp605NegHel", &pCosA605NegHel, "pCosA605NegHel/D");
  tree->Branch("pCosAmpErr605NegHel", &pCosAErr605NegHel, "pCosAErr605NegHel/D");
  tree->Branch("pCosAmp615NegHel", &pCosA615NegHel, "pCosA615NegHel/D");
  tree->Branch("pCosAmpErr615NegHel", &pCosAErr615NegHel, "pCosAErr615NegHel/D");

  tree->Branch("pCosAmp425PosHel", &pCosA425PosHel, "pCosA425PosHel/D");
  tree->Branch("pCosAmpErr425PosHel", &pCosAErr425PosHel, "pCosAErr425PosHel/D");
  tree->Branch("pCosAmp435PosHel", &pCosA435PosHel, "pCosA435PosHel/D");
  tree->Branch("pCosAmpErr435PosHel", &pCosAErr435PosHel, "pCosAErr435PosHel/D");
  tree->Branch("pCosAmp445PosHel", &pCosA445PosHel, "pCosA445PosHel/D");
  tree->Branch("pCosAmpErr445PosHel", &pCosAErr445PosHel, "pCosAErr445PosHel/D");
  tree->Branch("pCosAmp455PosHel", &pCosA455PosHel, "pCosA455PosHel/D");
  tree->Branch("pCosAmpErr455PosHel", &pCosAErr455PosHel, "pCosAErr455PosHel/D");
  tree->Branch("pCosAmp465PosHel", &pCosA465PosHel, "pCosA465PosHel/D");
  tree->Branch("pCosAmpErr465PosHel", &pCosAErr465PosHel, "pCosAErr465PosHel/D");
  tree->Branch("pCosAmp475PosHel", &pCosA475PosHel, "pCosA475PosHel/D");
  tree->Branch("pCosAmpErr475PosHel", &pCosAErr475PosHel, "pCosAErr475PosHel/D");
  tree->Branch("pCosAmp485PosHel", &pCosA485PosHel, "pCosA485PosHel/D");
  tree->Branch("pCosAmpErr485PosHel", &pCosAErr485PosHel, "pCosAErr485PosHel/D");
  tree->Branch("pCosAmp495PosHel", &pCosA495PosHel, "pCosA495PosHel/D");
  tree->Branch("pCosAmpErr495PosHel", &pCosAErr495PosHel, "pCosAErr495PosHel/D");
  tree->Branch("pCosAmp505PosHel", &pCosA505PosHel, "pCosA505PosHel/D");
  tree->Branch("pCosAmpErr505PosHel", &pCosAErr505PosHel, "pCosAErr505PosHel/D");
  tree->Branch("pCosAmp515PosHel", &pCosA515PosHel, "pCosA515PosHel/D");
  tree->Branch("pCosAmpErr515PosHel", &pCosAErr515PosHel, "pCosAErr515PosHel/D");
  tree->Branch("pCosAmp525PosHel", &pCosA525PosHel, "pCosA525PosHel/D");
  tree->Branch("pCosAmpErr525PosHel", &pCosAErr525PosHel, "pCosAErr525PosHel/D");
  tree->Branch("pCosAmp535PosHel", &pCosA535PosHel, "pCosA535PosHel/D");
  tree->Branch("pCosAmpErr535PosHel", &pCosAErr535PosHel, "pCosAErr535PosHel/D");
  tree->Branch("pCosAmp545PosHel", &pCosA545PosHel, "pCosA545PosHel/D");
  tree->Branch("pCosAmpErr545PosHel", &pCosAErr545PosHel, "pCosAErr545PosHel/D");
  tree->Branch("pCosAmp555PosHel", &pCosA555PosHel, "pCosA555PosHel/D");
  tree->Branch("pCosAmpErr555PosHel", &pCosAErr555PosHel, "pCosAErr555PosHel/D");
  tree->Branch("pCosAmp565PosHel", &pCosA565PosHel, "pCosA565PosHel/D");
  tree->Branch("pCosAmpErr565PosHel", &pCosAErr565PosHel, "pCosAErr565PosHel/D");
  tree->Branch("pCosAmp575PosHel", &pCosA575PosHel, "pCosA575PosHel/D");
  tree->Branch("pCosAmpErr575PosHel", &pCosAErr575PosHel, "pCosAErr575PosHel/D");
  tree->Branch("pCosAmp585PosHel", &pCosA585PosHel, "pCosA585PosHel/D");
  tree->Branch("pCosAmpErr585PosHel", &pCosAErr585PosHel, "pCosAErr585PosHel/D");
  tree->Branch("pCosAmp595PosHel", &pCosA595PosHel, "pCosA595PosHel/D");
  tree->Branch("pCosAmpErr595PosHel", &pCosAErr595PosHel, "pCosAErr595PosHel/D");
  tree->Branch("pCosAmp605PosHel", &pCosA605PosHel, "pCosA605PosHel/D");
  tree->Branch("pCosAmpErr605PosHel", &pCosAErr605PosHel, "pCosAErr605PosHel/D");
  tree->Branch("pCosAmp615PosHel", &pCosA615PosHel, "pCosA615PosHel/D");
  tree->Branch("pCosAmpErr615PosHel", &pCosAErr615PosHel, "pCosAErr615PosHel/D");

  // Fill branches (and hence tree) with corresponding parameters from above
  for (Int_t m = 0; m < 10; m++){

    pCosA425NegHel = pCosAmpNegHel[m][0];
    pCosAErr425NegHel = pCosAmpErrNegHel[m][0];
    pCosA435NegHel = pCosAmpNegHel[m][1];
    pCosAErr435NegHel = pCosAmpErrNegHel[m][1];
    pCosA445NegHel = pCosAmpNegHel[m][2];
    pCosAErr445NegHel = pCosAmpErrNegHel[m][2];
    pCosA455NegHel = pCosAmpNegHel[m][3];
    pCosAErr455NegHel = pCosAmpErrNegHel[m][3];
    pCosA465NegHel = pCosAmpNegHel[m][4];
    pCosAErr465NegHel = pCosAmpErrNegHel[m][4];
    pCosA475NegHel = pCosAmpNegHel[m][5];
    pCosAErr475NegHel = pCosAmpErrNegHel[m][5];
    pCosA485NegHel = pCosAmpNegHel[m][6];
    pCosAErr485NegHel = pCosAmpErrNegHel[m][6];
    pCosA495NegHel = pCosAmpNegHel[m][7];
    pCosAErr495NegHel = pCosAmpErrNegHel[m][7];
    pCosA505NegHel = pCosAmpNegHel[m][8];
    pCosAErr505NegHel = pCosAmpErrNegHel[m][8];
    pCosA515NegHel = pCosAmpNegHel[m][9];
    pCosAErr515NegHel = pCosAmpErrNegHel[m][9];
    pCosA525NegHel = pCosAmpNegHel[m][10];
    pCosAErr525NegHel = pCosAmpErrNegHel[m][10];
    pCosA535NegHel = pCosAmpNegHel[m][11];
    pCosAErr535NegHel = pCosAmpErrNegHel[m][11];
    pCosA545NegHel = pCosAmpNegHel[m][12];
    pCosAErr545NegHel= pCosAmpErrNegHel[m][12];
    //pCosA555NegHel = pCosAmpNegHel[m][13];
    //pCosAErr555NegHel = pCosAmpErrNegHel[m][13];
    pCosA555NegHel = 0;
    pCosAErr555NegHel = 0;
    pCosA565NegHel = pCosAmpNegHel[m][14];
    pCosAErr565NegHel = pCosAmpErrNegHel[m][14];
    pCosA575NegHel = pCosAmpNegHel[m][15];
    pCosAErr575NegHel = pCosAmpErrNegHel[m][15];
    pCosA585NegHel = pCosAmpNegHel[m][16];
    pCosAErr585NegHel = pCosAmpErrNegHel[m][16];
    pCosA595NegHel = pCosAmpNegHel[m][17];
    pCosAErr595NegHel = pCosAmpErrNegHel[m][17];
    pCosA605NegHel = pCosAmpNegHel[m][18];
    pCosAErr605NegHel = pCosAmpErrNegHel[m][18];
    pCosA615NegHel = pCosAmpNegHel[m][19];
    pCosAErr615NegHel = pCosAmpErrNegHel[m][19];

    pCosA425PosHel = pCosAmpPosHel[m][0];
    pCosAErr425PosHel = pCosAmpErrPosHel[m][0];
    pCosA435PosHel = pCosAmpPosHel[m][1];
    pCosAErr435PosHel = pCosAmpErrPosHel[m][1];
    pCosA445PosHel = pCosAmpPosHel[m][2];
    pCosAErr445PosHel = pCosAmpErrPosHel[m][2];
    pCosA455PosHel = pCosAmpPosHel[m][3];
    pCosAErr455PosHel = pCosAmpErrPosHel[m][3];
    pCosA465PosHel = pCosAmpPosHel[m][4];
    pCosAErr465PosHel = pCosAmpErrPosHel[m][4];
    pCosA475PosHel = pCosAmpPosHel[m][5];
    pCosAErr475PosHel = pCosAmpErrPosHel[m][5];
    pCosA485PosHel = pCosAmpPosHel[m][6];
    pCosAErr485PosHel = pCosAmpErrPosHel[m][6];
    pCosA495PosHel = pCosAmpPosHel[m][7];
    pCosAErr495PosHel = pCosAmpErrPosHel[m][7];
    pCosA505PosHel = pCosAmpPosHel[m][8];
    pCosAErr505PosHel = pCosAmpErrPosHel[m][8];
    pCosA515PosHel = pCosAmpPosHel[m][9];
    pCosAErr515PosHel = pCosAmpErrPosHel[m][9];
    pCosA525PosHel = pCosAmpPosHel[m][10];
    pCosAErr525PosHel = pCosAmpErrPosHel[m][10];
    pCosA535PosHel = pCosAmpPosHel[m][11];
    pCosAErr535PosHel = pCosAmpErrPosHel[m][11];
    pCosA545PosHel = pCosAmpPosHel[m][12];
    pCosAErr545PosHel = pCosAmpErrPosHel[m][12];
    //pCosA555PosHel = pCosAmpPosHel[m][13];
    //pCosAErr555PosHel = pCosAmpErrPosHel[m][13];
    pCosA555PosHel = 0;
    pCosAErr555PosHel = 0;
    pCosA565PosHel = pCosAmpPosHel[m][14];
    pCosAErr565PosHel = pCosAmpErrPosHel[m][14];
    pCosA575PosHel = pCosAmpPosHel[m][15];
    pCosAErr575PosHel = pCosAmpErrPosHel[m][15];
    pCosA585PosHel = pCosAmpPosHel[m][16];
    pCosAErr585PosHel = pCosAmpErrPosHel[m][16];
    pCosA595PosHel = pCosAmpPosHel[m][17];
    pCosAErr595PosHel = pCosAmpErrPosHel[m][17];
    pCosA605PosHel = pCosAmpPosHel[m][18];
    pCosAErr605PosHel = pCosAmpErrPosHel[m][18];
    pCosA615PosHel = pCosAmpPosHel[m][19];
    pCosAErr615PosHel = pCosAmpErrPosHel[m][19];

    tree->Fill();
  }

  f1.Write();

}
