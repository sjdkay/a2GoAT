#include "./includes.h"

void PlotSaver(){

  TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/Amo/Physics_Total_12_21_2_17.root"); // Open the latest PTotal file to load histograms from
  TText *warn = new TText(0, 0 ,"PRELIMINARY"); // Preliminary warning label text

  // Missing mass (with banana cut) across various photon energy bins
  TCanvas *canvas = new TCanvas("canvas","canvas", 1920, 1080);
  MMp200300->SetLineColor(1);
  MMp300400->SetLineColor(2);   
  MMp400500->SetLineColor(3);
  MMp500600->SetLineColor(4);
  MMp600700->SetLineColor(6);
  MMp700800->SetLineColor(7);
  MMp800900->SetLineColor(807);
  MMp200300->SetTitle("Missing Mass as Seen By Proton Across EGamma Bins");
  MMp200300->SetStats(kFALSE);
  MMp200300->Draw();
  MMp300400->Draw("Same");
  MMp400500->Draw("Same");
  MMp500600->Draw("Same");
  MMp600700->Draw("Same");
  MMp700800->Draw("Same");
  MMp800900->Draw("Same");
  leg = new TLegend(0.8, 0.70, 0.9, 0.90);
  leg->AddEntry(MMp200300, "EGamma 200-300", "l");
  leg->AddEntry(MMp300400, "EGamma 300-400", "l");
  leg->AddEntry(MMp400500, "EGamma 400-500", "l");
  leg->AddEntry(MMp500600, "EGamma 500-600", "l");
  leg->AddEntry(MMp600700, "EGamma 600-700", "l");
  leg->AddEntry(MMp700800, "EGamma 700-800", "l");
  leg->AddEntry(MMp800900, "EGamma 800-900", "l");
  leg->Draw("Same");
  canvas->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/MMpComparison_EGamma_12.png");
  canvas->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/MMpComparison_EGamma_12.pdf");

  // Missing mass plot with/without banana cut
  TCanvas *canvas2 = new TCanvas("canvas2","canvas2", 1920, 1080);
  MMpEpCorrected->SetLineColor(4);
  MMpEpCorrectedCut->SetLineColor(2);
  MMpEpCorrected->SetTitle("Missing Mass as Seen by Proton");
  MMpEpCorrected->SetStats(kFALSE);
  leg = new TLegend(0.8, 0.80, 0.9, 0.90);
  leg->AddEntry(MMpEpCorrected, "MMp from EpCorrected", "l");
  leg->AddEntry(MMpEpCorrectedCut, "MMp from EpCorrected with p banana cut", "l");
  MMpEpCorrected->Draw();
  MMpEpCorrectedCut->Draw("Same");
  leg->Draw("Same");
  canvas2->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/MMpEpCorr_CutComparison_12.png");
  canvas2->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/MMpEpCorr_CutComparison_12.pdf");

  // Opening angle plot, comparison between full data and banana cut
  TCanvas *canvas4 = new TCanvas("canvas4","canvas4", 1920, 1080);
  OAngle->SetLineColor(4);
  OAngleCut->SetLineColor(2);
  OAngle->SetTitle("Opening angle between P and N Vectors");
  OAngle->SetStats(kFALSE);
  leg = new TLegend(0.75, 0.80, 0.9, 0.90);
  leg->AddEntry(OAngle, "Opening Angle", "l");
  leg->AddEntry(OAngleCut, "Opening Angle with p banana cut", "l");
  OAngle->Draw();
  OAngleCut->Draw("Same");
  leg->Draw("Same");
  canvas4->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/OAngle_CutComparison_12.png");
  canvas4->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/OAngle_CutComparison_12.pdf");

  // Banana plot
  TCanvas *canvas5 = new TCanvas("canvas5","canvas5", 1920, 1080);
  E_dE->SetXTitle("CB Energy/MeV");
  E_dE->SetYTitle("PID Energy/MeV");
  E_dE->Draw("Col");
  TFile *f2 = new TFile("/scratch/Mainz_Software/a2GoAT/configfiles/cuts/CB_DeltaE-E_Proton_7_12_16.root"); // Open the proton cut file used
  Proton->SetLineWidth(5);
  Proton->Draw("Same");
  canvas5->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/E_dE_12.png");
  canvas5->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/E_dE_12.pdf");
  TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/Amo/Physics_Total_12_21_2_17.root"); // Open the latest PTotal file to load histograms from

  // Banana plot showing cut region
  TCanvas *canvas6 = new TCanvas("canvas6","canvas6", 1920, 1080);
  E_dE_Cut->SetXTitle("CB Energy/MeV");
  E_dE_Cut->SetYTitle("PID Energy/MeV");
  E_dE_Cut->Draw("Col");
  TFile *f2 = new TFile("/scratch/Mainz_Software/a2GoAT/configfiles/cuts/CB_DeltaE-E_Proton_7_12_16.root"); // Open the proton cut file used
  Proton->SetLineWidth(5);
  Proton->Draw("Same");
  canvas6->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/E_dE_Cut_12.png");
  canvas6->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/E_dE_Cut_12.pdf");
  TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/Amo/Physics_Total_12_21_2_17.root"); // Open the latest PTotal file to load histograms from
  
  // Banana plot for KinEpdEp
  TCanvas *canvas7 = new TCanvas("canvas7","canvas7", 1920, 1080);
  KinEp_dE->SetXTitle("KinEp/MeV");
  KinEp_dE->SetYTitle("PID Energy/MeV");
  KinEp_dE->Draw("Col");
  TFile *f2 = new TFile("/scratch/Mainz_Software/a2GoAT/configfiles/cuts/CB_DeltaE-E_ProtonKinGood_03_02_17.root"); // Open the proton cut file used
  ProtonKinGood->SetLineWidth(5);
  ProtonKinGood->Draw("Same");
  TFile *f3 = new TFile("/scratch/Mainz_Software/a2GoAT/configfiles/cuts/CB_DeltaE-E_ProtonKinBad_15_12_16.root"); // Open the proton cut file used
  ProtonKinBad->SetLineWidth(5);
  ProtonKinBad->SetLineColor(2);
  ProtonKinBad->Draw("Same");
  canvas7->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/KinEp_dE_12.png");
  canvas7->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/KinEp_dE_12.pdf");
  TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/Amo/Physics_Total_12_21_2_17.root"); // Open the latest PTotal file to load histograms from

  // Banana plot showing good cut region for KinEpdEp
  TCanvas *canvas8 = new TCanvas("canvas8","canvas8", 1920, 1080);
  KinEp_dE_GoodCut->SetXTitle("KinEp/MeV");
  KinEp_dE_GoodCut->SetYTitle("PID Energy/MeV");
  KinEp_dE_GoodCut->Draw("Col");
  TFile *f2 = new TFile("/scratch/Mainz_Software/a2GoAT/configfiles/cuts/CB_DeltaE-E_ProtonKinGood_03_02_17.root"); // Open the proton cut file used
  ProtonKinGood->SetLineWidth(5);
  ProtonKinGood->Draw("Same");
  canvas8->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/KinEp_dE_GoodCut_12.png");
  canvas8->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/KinEp_dE_GoodCut_12.pdf");
  TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/Amo/Physics_Total_12_21_2_17.root"); // Open the latest PTotal file to load histograms from

  // ThetaSc Plot
  TCanvas *canvas9 = new TCanvas("canvas9","canvas9", 1920, 1080);
  Theta_Scattered->SetXTitle("Theta/Degrees");
  Theta_Scattered->SetLineColor(4);
  Theta_Scattered->SetStats(kFALSE);
  Theta_Scattered->Draw("Col");
  canvas9->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/ThetaSc_12.png");
  canvas9->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/ThetaSc_12.pdf");

  // PhiSc Plot
  TCanvas *canvas10 = new TCanvas("canvas10","canvas10", 1920, 1080);
  Phi_Scattered->SetXTitle("Phi/Degrees");
  Phi_Scattered->SetLineColor(4);
  Phi_Scattered->SetStats(kFALSE);
  Phi_Scattered->Draw("Col");
  canvas10->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/PhiSc_12.png");
  canvas10->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/PhiSc_12.pdf");

  // PhiSc Plot across Photon E bins
  TCanvas *canvas11 = new TCanvas("canvas11","canvas11", 1920, 1080);
  Phi_Scattered_375MeV->SetLineColor(417);
  Phi_Scattered_425MeV->SetLineColor(4);
  Phi_Scattered_475MeV->SetLineColor(2);
  Phi_Scattered_525MeV->SetLineColor(1);
  Phi_Scattered_375MeV->SetStats(kFALSE);
  Phi_Scattered_375MeV->SetTitle("Phi in Scattered Frame Across EGamma Bins (350-550MeV)");
  Phi_Scattered_375MeV->SetXTitle("Phi/Degrees");
  Phi_Scattered_375MeV->Draw();
  Phi_Scattered_425MeV->Draw("Same");
  Phi_Scattered_475MeV->Draw("Same");
  Phi_Scattered_525MeV->Draw("Same");
  leg = new TLegend(0.8, 0.70, 0.9, 0.90);
  leg->AddEntry(Phi_Scattered_375MeV, "EGamma 350-400", "l");
  leg->AddEntry(Phi_Scattered_425MeV, "EGamma 400-450", "l");
  leg->AddEntry(Phi_Scattered_475MeV, "EGamma 450-500", "l");
  leg->AddEntry(Phi_Scattered_525MeV, "EGamma 500-550", "l");
  leg->Draw("Same");
  canvas11->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/PhiSc_EGamma_Low_12.png");
  canvas11->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/PhiSc_EGamma_Low_12.pdf");

  // PhiSc Plot across Photon E bins
  TCanvas *canvas12 = new TCanvas("canvas12","canvas12", 1920, 1080);
  Phi_Scattered_575MeV->SetLineColor(1);
  Phi_Scattered_625MeV->SetLineColor(807);
  Phi_Scattered_675MeV->SetLineColor(3);
  Phi_Scattered_725MeV->SetLineColor(4);
  Phi_Scattered_775MeV->SetLineColor(880);
  Phi_Scattered_825MeV->SetLineColor(2);
  Phi_Scattered_875MeV->SetLineColor(30);
  Phi_Scattered_575MeV->SetStats(kFALSE);
  Phi_Scattered_575MeV->SetTitle("Phi in Scattered Frame Across EGamma Bins (550-900MeV)");
  Phi_Scattered_575MeV->SetXTitle("Phi/Degrees");
  Phi_Scattered_575MeV->Draw();
  Phi_Scattered_625MeV->Draw("Same");
  Phi_Scattered_675MeV->Draw("Same");
  Phi_Scattered_725MeV->Draw("Same");
  Phi_Scattered_775MeV->Draw("Same");
  Phi_Scattered_825MeV->Draw("Same");
  Phi_Scattered_875MeV->Draw("Same");
  leg = new TLegend(0.8, 0.70, 0.9, 0.90);
  leg->AddEntry(Phi_Scattered_575MeV, "EGamma 550-600", "l");
  leg->AddEntry(Phi_Scattered_625MeV, "EGamma 600-650", "l");
  leg->AddEntry(Phi_Scattered_675MeV, "EGamma 650-700", "l");
  leg->AddEntry(Phi_Scattered_725MeV, "EGamma 700-750", "l");
  leg->AddEntry(Phi_Scattered_775MeV, "EGamma 750-800", "l");
  leg->AddEntry(Phi_Scattered_825MeV, "EGamma 800-850", "l");
  leg->AddEntry(Phi_Scattered_875MeV, "EGamma 850-900", "l");
  leg->Draw("Same");
  canvas12->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/PhiSc_EGamma_High_12.png");
  canvas12->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/PhiSc_EGamma_High_12.pdf");

  //PhiSc as fn of ThetaSc
  TCanvas *canvas13 = new TCanvas("canvas13","canvas13", 1920, 1080);
  ThetaScPhiSc->SetXTitle("ThetaSc/Deg");
  ThetaScPhiSc->SetYTitle("PhiSc/Deg");
  ThetaScPhiSc->SetStats(kFALSE);
  ThetaScPhiSc->Draw("Col");
  canvas13->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/ThetaScPhiSc_12.png");
  canvas13->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/ThetaScPhiSc_12.pdf");

  //OAngle Egammacut
  TCanvas *canvas14 = new TCanvas("canvas14","canvas14", 1920, 1080);
  OAngleCut200400->SetXTitle("Opening Angle/Deg");
  OAngleCut200400->SetStats(kFALSE);
  OAngleCut200400->Draw();
  canvas14->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/OAngleCut200400_12.png");
  canvas14->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/OAngleCut200400_12.pdf");

  //OAnglePhiDiffLab
  TCanvas *canvas15 = new TCanvas("canvas15","canvas15", 1920, 1080);
  OAnglePhiDiffLab->SetXTitle("Opening Angle/Deg");
  OAnglePhiDiffLab->SetYTitle("PhiScLab/Deg");
  OAnglePhiDiffLab->SetStats(kFALSE);
  OAnglePhiDiffLab->Draw("Col");
  canvas15->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/OAnglePhiDiffLab_12.png");
  canvas15->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/OAnglePhiDiffLab_12.pdf");

  //ThetaDiffPhiDiffLab
  TCanvas *canvas16 = new TCanvas("canvas16","canvas16", 1920, 1080);
  ThetaDiffPhiDiffLab->SetXTitle("ThetaDiff/Deg");
  ThetaDiffPhiDiffLab->SetYTitle("PhiScLab/Deg");
  ThetaDiffPhiDiffLab->SetStats(kFALSE);
  ThetaDiffPhiDiffLab->Draw("Col");
  canvas16->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/ThetaDiffPhiDiffLab_12.png");
  canvas16->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_23_01_17/ThetaDiffPhiDiffLab_12.pdf");

}
