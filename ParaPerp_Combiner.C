#include "./includes_ParaPerpCombiner.h"

void ParaPerp_Combiner(){

  TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/Para/Physics_Total_Para_2_6_3_17.root"); // Open latest Para file
  TH1D* time_Para = (TH1D*)f->Get("time")->Clone();
  time_Para->SetName("time_Para");
  TH1D* time_cut_Para = (TH1D*)f->Get("time_cut")->Clone();
  time_Para->SetName("time_cut_Para");

  TFile f2("ParaPerp_Total_2_Combined_TEST.root", "RECREATE");
  time_Para->Write();
  f2.Write();
 
}
