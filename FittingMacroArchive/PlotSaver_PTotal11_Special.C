#include "./includes.h"

void PlotSaver_PTotal11_Special(){

  TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_15_11_16/Physics_Total_11_15_12_16.root"); // Open the latest PTotal file to load histograms from
  TText *warn = new TText(0, 0 ,"PRELIMINARY"); // Preliminary warning label text

  // Missing mass (with good banana cut) across various photon energy bins
  TCanvas *canvas = new TCanvas("canvas","canvas", 1920, 1080);
  MMp200300GoodCut->SetLineColor(1);
  MMp300400GoodCut->SetLineColor(2);   
  MMp400500GoodCut->SetLineColor(3);
  MMp500600GoodCut->SetLineColor(4);
  MMp600700GoodCut->SetLineColor(6);
  MMp700800GoodCut->SetLineColor(7);
  MMp800900GoodCut->SetLineColor(807);
  MMp200300GoodCut->SetTitle("Missing Mass as Seen By Proton Across EGamma Bins (Good P Cut)");
  MMp200300GoodCut->SetStats(kFALSE);
  MMp200300GoodCut->Draw();
  MMp300400GoodCut->Draw("Same");
  MMp400500GoodCut->Draw("Same");
  MMp500600GoodCut->Draw("Same");
  MMp600700GoodCut->Draw("Same");
  MMp700800GoodCut->Draw("Same");
  MMp800900GoodCut->Draw("Same");
  leg = new TLegend(0.8, 0.70, 0.9, 0.90);
  leg->AddEntry(MMp200300GoodCut, "EGamma 200-300", "l");
  leg->AddEntry(MMp300400GoodCut, "EGamma 300-400", "l");
  leg->AddEntry(MMp400500GoodCut, "EGamma 400-500", "l");
  leg->AddEntry(MMp500600GoodCut, "EGamma 500-600", "l");
  leg->AddEntry(MMp600700GoodCut, "EGamma 600-700", "l");
  leg->AddEntry(MMp700800GoodCut, "EGamma 700-800", "l");
  leg->AddEntry(MMp800900GoodCut, "EGamma 800-900", "l");
  leg->Draw("Same");
  canvas->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_15_11_16/MMpComparison_EGamma_GoodCut_11.png");
  canvas->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_15_11_16/MMpComparison_EGamma_GoodCut_11.pdf");

  // Missing mass (with bad banana cut) across various photon energy bins
  TCanvas *canvas2 = new TCanvas("canvas2","canvas2", 1920, 1080);
  MMp300400BadCut->SetLineColor(2);   
  MMp400500BadCut->SetLineColor(3);
  MMp500600BadCut->SetLineColor(4);
  MMp600700BadCut->SetLineColor(6);
  MMp700800BadCut->SetLineColor(7);
  MMp800900BadCut->SetLineColor(807);
  MMp500600BadCut->SetTitle("Missing Mass as Seen By Proton Across EGamma Bins (Bad P Cut)");
  MMp500600BadCut->SetStats(kFALSE);
  MMp500600BadCut->Draw();
  MMp400500BadCut->Draw("Same");
  MMp300400BadCut->Draw("Same");
  MMp600700BadCut->Draw("Same");
  MMp700800BadCut->Draw("Same");
  MMp800900BadCut->Draw("Same");
  leg = new TLegend(0.8, 0.70, 0.9, 0.90);
  leg->AddEntry(MMp300400BadCut, "EGamma 300-400", "l");
  leg->AddEntry(MMp400500BadCut, "EGamma 400-500", "l");
  leg->AddEntry(MMp500600BadCut, "EGamma 500-600", "l");
  leg->AddEntry(MMp600700BadCut, "EGamma 600-700", "l");
  leg->AddEntry(MMp700800BadCut, "EGamma 700-800", "l");
  leg->AddEntry(MMp800900BadCut, "EGamma 800-900", "l");
  leg->Draw("Same");
  canvas2->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_15_11_16/MMpComparison_EGamma_BadCut_11.png");
  canvas2->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_15_11_16/MMpComparison_EGamma_BadCut_11.pdf");

  // Missing mass plot with/without good/bad banana cut
  TCanvas *canvas3 = new TCanvas("canvas3","canvas3", 1920, 1080);
  MMpEpCorrected->SetLineColor(4);
  MMpEpCorrectedGoodCut->SetLineColor(2);
  MMpEpCorrectedBadCut->SetLineColor(3);
  MMpEpCorrected->SetTitle("Missing Mass as Seen by Proton");
  MMpEpCorrected->SetStats(kFALSE);
  leg = new TLegend(0.8, 0.80, 0.9, 0.90);
  leg->AddEntry(MMpEpCorrected, "MMp from EpCorrected", "l");
  leg->AddEntry(MMpEpCorrectedGoodCut, "MMp from EpCorrected with good p banana cut", "l");
  leg->AddEntry(MMpEpCorrectedBadCut, "MMp from EpCorrected with bad p banana cut", "l");
  MMpEpCorrected->Draw();
  MMpEpCorrectedGoodCut->Draw("Same");
  MMpEpCorrectedBadCut->Draw("Same");
  leg->Draw("Same");
  canvas3->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_15_11_16/MMpEpCorr_CutComparison_11.png");
  canvas3->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_15_11_16/MMpEpCorr_CutComparison_11.pdf");

  // Diff between ThetanWC and ThetanRec across Photon energy range for good p cut
  TCanvas *canvas4 = new TCanvas("canvas4","canvas4", 1920, 1080);
  ThetanWCThetanRecDiff200300GoodCut->SetLineColor(1);
  ThetanWCThetanRecDiff300400GoodCut->SetLineColor(2);   
  ThetanWCThetanRecDiff400500GoodCut->SetLineColor(3);
  ThetanWCThetanRecDiff500600GoodCut->SetLineColor(4);
  ThetanWCThetanRecDiff600700GoodCut->SetLineColor(6);
  ThetanWCThetanRecDiff700800GoodCut->SetLineColor(7);
  ThetanWCThetanRecDiff800900GoodCut->SetLineColor(807);
  ThetanWCThetanRecDiff200300GoodCut->SetTitle("Differnce Between ThetanWC and ThetanRec Across EGamma Bins (Good p cut)");
  ThetanWCThetanRecDiff200300GoodCut->SetStats(kFALSE);
  ThetanWCThetanRecDiff200300GoodCut->Draw();
  ThetanWCThetanRecDiff300400GoodCut->Draw("Same");
  ThetanWCThetanRecDiff400500GoodCut->Draw("Same");
  ThetanWCThetanRecDiff500600GoodCut->Draw("Same");
  ThetanWCThetanRecDiff600700GoodCut->Draw("Same");
  ThetanWCThetanRecDiff700800GoodCut->Draw("Same");
  ThetanWCThetanRecDiff800900GoodCut->Draw("Same");
  leg = new TLegend(0.65, 0.70, 0.9, 0.90);
  leg->AddEntry(ThetanWCThetanRecDiff200300GoodCut, "Thetan WCvsRec Diff over EGamma 200-300", "l");
  leg->AddEntry(ThetanWCThetanRecDiff300400GoodCut, "Thetan WCvsRec Diff over EGamma 300-400", "l");
  leg->AddEntry(ThetanWCThetanRecDiff400500GoodCut, "Thetan WCvsRec Diff over EGamma 400-500", "l");
  leg->AddEntry(ThetanWCThetanRecDiff500600GoodCut, "Thetan WCvsRec Diff over EGamma 500-600", "l");
  leg->AddEntry(ThetanWCThetanRecDiff600700GoodCut, "Thetan WCvsRec Diff over EGamma 600-700", "l");
  leg->AddEntry(ThetanWCThetanRecDiff700800GoodCut, "Thetan WCvsRec Diff over EGamma 700-800", "l");
  leg->AddEntry(ThetanWCThetanRecDiff800900GoodCut, "Thetan WCvsRec Diff over EGamma 800-900", "l");
  leg->Draw("Same");
  canvas4->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_15_11_16/ThetanWCThetanRecDiff_EGamma_GoodCut_11.png");
  canvas4->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_15_11_16/ThetanWCThetanRecDiff_EGamma_GoodCut_11.pdf");

  // Diff between ThetanWC and ThetanRec across Photon energy range for bad p cut
  TCanvas *canvas5 = new TCanvas("canvas5","canvas5", 1920, 1080);
  ThetanWCThetanRecDiff300400BadCut->SetLineColor(2);   
  ThetanWCThetanRecDiff400500BadCut->SetLineColor(3);
  ThetanWCThetanRecDiff500600BadCut->SetLineColor(4);
  ThetanWCThetanRecDiff600700BadCut->SetLineColor(6);
  ThetanWCThetanRecDiff700800BadCut->SetLineColor(7);
  ThetanWCThetanRecDiff800900BadCut->SetLineColor(1);
  ThetanWCThetanRecDiff500600BadCut->SetTitle("Differnce Between ThetanWC and ThetanRec Across EGamma Bins (Bad p cut)");
  ThetanWCThetanRecDiff500600BadCut->SetStats(kFALSE);
  ThetanWCThetanRecDiff500600BadCut->Draw();
  ThetanWCThetanRecDiff400500BadCut->Draw("Same");
  ThetanWCThetanRecDiff300400BadCut->Draw("Same");
  ThetanWCThetanRecDiff600700BadCut->Draw("Same");
  ThetanWCThetanRecDiff700800BadCut->Draw("Same");
  ThetanWCThetanRecDiff800900BadCut->Draw("Same");
  leg = new TLegend(0.65, 0.70, 0.9, 0.90);
  leg->AddEntry(ThetanWCThetanRecDiff300400BadCut, "Thetan WCvsRec Diff over EGamma 300-400", "l");
  leg->AddEntry(ThetanWCThetanRecDiff400500BadCut, "Thetan WCvsRec Diff over EGamma 400-500", "l");
  leg->AddEntry(ThetanWCThetanRecDiff500600BadCut, "Thetan WCvsRec Diff over EGamma 500-600", "l");
  leg->AddEntry(ThetanWCThetanRecDiff600700BadCut, "Thetan WCvsRec Diff over EGamma 600-700", "l");
  leg->AddEntry(ThetanWCThetanRecDiff700800BadCut, "Thetan WCvsRec Diff over EGamma 700-800", "l");
  leg->AddEntry(ThetanWCThetanRecDiff800900BadCut, "Thetan WCvsRec Diff over EGamma 800-900", "l");
  leg->Draw("Same");
  canvas5->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_15_11_16/ThetanWCThetanRecDiff_EGamma_BadCut_11.png");
  canvas5->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_15_11_16/ThetanWCThetanRecDiff_EGamma_BadCut_11.pdf");

  // Opening angle plot, comparison between full data and good/bad banana cut
  TCanvas *canvas6 = new TCanvas("canvas6","canvas6", 1920, 1080);
  OAngle->SetLineColor(4);
  OAngleGoodCut->SetLineColor(2);
  OAngleBadCut->SetLineColor(3);
  OAngle->SetTitle("Opening angle between P and N Vectors");
  OAngle->SetStats(kFALSE);
  leg = new TLegend(0.75, 0.80, 0.9, 0.90);
  leg->AddEntry(OAngle, "Opening Angle", "l");
  leg->AddEntry(OAngleGoodCut, "Opening Angle with good p banana cut", "l");
  leg->AddEntry(OAngleBadCut, "Opening Angle with bad p banana cut", "l");
  OAngle->Draw();
  OAngleGoodCut->Draw("Same");
  OAngleBadCut->Draw("Same");
  leg->Draw("Same");
  canvas6->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_15_11_16/OAngle_CutComparison_11.png");
  canvas6->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_15_11_16/OAngle_CutComparison_11.pdf");

  // Banana plot
  TCanvas *canvas7 = new TCanvas("canvas7","canvas7", 1920, 1080);
  E_dE->SetXTitle("CB Energy/MeV");
  E_dE->SetYTitle("PID Energy/MeV");
  E_dE->Draw("Col");
  TFile *f2 = new TFile("/scratch/Mainz_Software/a2GoAT/configfiles/cuts/CB_DeltaE-E_Proton_7_12_16.root"); // Open the proton cut file used
  Proton->SetLineWidth(5);
  Proton->Draw("Same");
  canvas7->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_15_11_16/E_dE_11.png");
  canvas7->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_15_11_16/E_dE_11.pdf");
  TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_15_11_16/Physics_Total_11_15_12_16.root"); // Open the latest PTotal file to load histograms from

  // Banana plot showing cut region
  TCanvas *canvas8 = new TCanvas("canvas8","canvas8", 1920, 1080);
  E_dE_Cut->SetXTitle("CB Energy/MeV"); //Need to change this to E_dE_Cut for future Total files
  E_dE_Cut->SetYTitle("PID Energy/MeV");
  E_dE_Cut->Draw("Col");
  TFile *f2 = new TFile("/scratch/Mainz_Software/a2GoAT/configfiles/cuts/CB_DeltaE-E_Proton_7_12_16.root"); // Open the proton cut file used
  Proton->SetLineWidth(5);
  Proton->Draw("Same");
  canvas8->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_15_11_16/E_dE_Cut_11.png");
  canvas8->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_15_11_16/E_dE_Cut_11.pdf");
  TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_15_11_16/Physics_Total_11_15_12_16.root"); // Open the latest PTotal file to load histograms from

  // Banana plot for KinEpdEp
  TCanvas *canvas5 = new TCanvas("canvas9","canvas9", 1920, 1080);
  KinEp_dE->SetXTitle("KinEp/MeV");
  KinEp_dE->SetYTitle("PID Energy/MeV");
  KinEp_dE->Draw("Col");
  TFile *f2 = new TFile("/scratch/Mainz_Software/a2GoAT/configfiles/cuts/CB_DeltaE-E_ProtonKinGood_15_12_16.root"); // Open the proton cut file used
  ProtonKinGood->SetLineWidth(5);
  ProtonKinGood->Draw("Same");
  TFile *f3 = new TFile("/scratch/Mainz_Software/a2GoAT/configfiles/cuts/CB_DeltaE-E_ProtonKinBad_15_12_16.root"); // Open the proton cut file used
  ProtonKinBad->SetLineWidth(5);
  ProtonKinBad->SetLineColor(2);
  ProtonKinBad->Draw("Same");
  canvas9->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_15_11_16/KinEp_dE_11.png");
  canvas9->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_15_11_16/KinEp_dE_11.pdf");
  TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_15_11_16/Physics_Total_11_15_12_16.root"); // Open the latest PTotal file to load histograms from

  // Banana plot showing good cut region for KinEpdEp
  TCanvas *canvas10 = new TCanvas("canvas10","canvas10", 1920, 1080);
  KinEp_dE_GoodCut->SetXTitle("CB Energy/MeV"); //Need to change this to E_dE_Cut for future Total files
  KinEp_dE_GoodCut->SetYTitle("PID Energy/MeV");
  KinEp_dE_GoodCut->Draw("Col");
  TFile *f2 = new TFile("/scratch/Mainz_Software/a2GoAT/configfiles/cuts/CB_DeltaE-E_ProtonKinGood_15_12_16.root"); // Open the proton cut file used
  ProtonKinGood->SetLineWidth(5);
  ProtonKinGood->Draw("Same");
  canvas10->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_15_11_16/KinEp_dE_GoodCut_11.png");
  canvas10->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_15_11_16/KinEp_dE_GoodCut_11.pdf");
  TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_15_11_16/Physics_Total_11_15_12_16.root"); // Open the latest PTotal file to load histograms from
 
  // Banana plot showing bad cut region for KinEpdEp
  TCanvas *canvas11 = new TCanvas("canvas11","canvas11", 1920, 1080);
  KinEp_dE_BadCut->SetXTitle("CB Energy/MeV"); //Need to change this to E_dE_Cut for future Total files
  KinEp_dE_BadCut->SetYTitle("PID Energy/MeV");
  KinEp_dE_BadCut->Draw("Col");
  TFile *f2 = new TFile("/scratch/Mainz_Software/a2GoAT/configfiles/cuts/CB_DeltaE-E_ProtonKinBad_15_12_16.root"); // Open the proton cut file used
  ProtonKinBad->SetLineWidth(5);
  ProtonKinBad->Draw("Same");
  canvas11->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_15_11_16/KinEp_dE_BadCut_11.png");
  canvas11->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_15_11_16/KinEp_dE_BadCut_11.pdf");
  TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_15_11_16/Physics_Total_11_15_12_16.root"); // Open the latest PTotal file to load histograms from

}
