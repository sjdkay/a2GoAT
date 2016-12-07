#include "./includes.h"

void PlotSaver(){

  TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_15_11_16/Physics_Total_8_05_12_16.root"); // Open the latest PTotal file to load histograms from
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
  canvas->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_15_11_16/MMpComparison_EGamma_8.png");
  canvas->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_15_11_16/MMpComparison_EGamma_8.pdf");

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
  canvas2->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_15_11_16/MMpEpCorr_CutComparison_8.png");
  canvas2->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_15_11_16/MMpEpCorr_CutComparison_8.pdf");

  // Diff between ThetanWC and ThetanRec across Photon energy range
  TCanvas *canvas3 = new TCanvas("canvas3","canvas3", 1920, 1080);
  ThetanWCThetanRecDiff200300->SetLineColor(1);
  ThetanWCThetanRecDiff300400->SetLineColor(2);   
  ThetanWCThetanRecDiff400500->SetLineColor(3);
  ThetanWCThetanRecDiff500600->SetLineColor(4);
  ThetanWCThetanRecDiff600700->SetLineColor(6);
  ThetanWCThetanRecDiff700800->SetLineColor(7);
  ThetanWCThetanRecDiff800900->SetLineColor(807);
  ThetanWCThetanRecDiff200300->SetTitle("Differnce Between ThetanWC and ThetanRec Across EGamma Bins");
  ThetanWCThetanRecDiff200300->SetStats(kFALSE);
  ThetanWCThetanRecDiff200300->Draw();
  ThetanWCThetanRecDiff300400->Draw("Same");
  ThetanWCThetanRecDiff400500->Draw("Same");
  ThetanWCThetanRecDiff500600->Draw("Same");
  ThetanWCThetanRecDiff600700->Draw("Same");
  ThetanWCThetanRecDiff700800->Draw("Same");
  ThetanWCThetanRecDiff800900->Draw("Same");
  leg = new TLegend(0.65, 0.70, 0.9, 0.90);
  leg->AddEntry(ThetanWCThetanRecDiff200300, "Thetan WCvsRec Diff over EGamma 200-300", "l");
  leg->AddEntry(ThetanWCThetanRecDiff300400, "Thetan WCvsRec Diff over EGamma 300-400", "l");
  leg->AddEntry(ThetanWCThetanRecDiff400500, "Thetan WCvsRec Diff over EGamma 400-500", "l");
  leg->AddEntry(ThetanWCThetanRecDiff500600, "Thetan WCvsRec Diff over EGamma 500-600", "l");
  leg->AddEntry(ThetanWCThetanRecDiff600700, "Thetan WCvsRec Diff over EGamma 600-700", "l");
  leg->AddEntry(ThetanWCThetanRecDiff700800, "Thetan WCvsRec Diff over EGamma 700-800", "l");
  leg->AddEntry(ThetanWCThetanRecDiff800900, "Thetan WCvsRec Diff over EGamma 800-900", "l");
  leg->Draw("Same");
  canvas3->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_15_11_16/ThetanWCThetanRecDiff_EGamma_8.png");
  canvas3->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_15_11_16/ThetanWCThetanRecDiff_EGamma_8.pdf");

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
  canvas4->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_15_11_16/OAngle_CutComparison_8.png");
  canvas4->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_15_11_16/OAngle_CutComparison_8.pdf");

  // Banana plot
  TCanvas *canvas5 = new TCanvas("canvas5","canvas5", 1920, 1080);
  E_dE->SetXTitle("CB Energy/MeV");
  E_dE->SetYTitle("PID Energy/MeV");
  E_dE->Draw("Col");
  TFile *f2 = new TFile("/scratch/Mainz_Software/a2GoAT/configfiles/cuts/CB_DeltaE-E_Proton_7_12_16.root"); // Open the proton cut file used
  Proton->SetLineWidth(5);
  Proton->Draw("Same");
  canvas5->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_15_11_16/E_dE_8.png");
  canvas5->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_15_11_16/E_dE_8.pdf");
  TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_15_11_16/Physics_Total_8_05_12_16.root"); // Open the latest PTotal file to load histograms from

  // Banana plot showing cut region
  TCanvas *canvas6 = new TCanvas("canvas6","canvas6", 1920, 1080);
  ECB_dE_Cut->SetXTitle("CB Energy/MeV"); //Need to change this to E_dE_Cut for future Total files, GH2 was misnamed in code
  ECB_dE_Cut->SetYTitle("PID Energy/MeV");
  ECB_dE_Cut->Draw("Col");
  TFile *f2 = new TFile("/scratch/Mainz_Software/a2GoAT/configfiles/cuts/CB_DeltaE-E_Proton_7_12_16.root"); // Open the proton cut file used
  Proton->SetLineWidth(5);
  Proton->Draw("Same");
  canvas6->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_15_11_16/E_dE_Cut_8.png");
  canvas6->SaveAs("/home/s1427339/Documents/Hadron\ Physics/Aug_16_Data_Plots/GoAT_15_11_16/E_dE_Cut_8.pdf");
  TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_15_11_16/Physics_Total_8_05_12_16.root"); // Open the latest PTotal file to load histograms from
  
}
