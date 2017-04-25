#include "./includes.h"

void ParaPerp_Combiner(){

  double NPara;
  double NPerp;
  double ScalingFactor;
  TH1D* PhiScNegHel_Combined;
  TH1D* PhiScPosHel_Combined;
  TH1D* PhiSc275NegHel_Combined;
  TH1D* PhiSc325NegHel_Combined;
  TH1D* PhiSc375NegHel_Combined;
  TH1D* PhiSc425NegHel_Combined;
  TH1D* PhiSc475NegHel_Combined;
  TH1D* PhiSc525NegHel_Combined;
  TH1D* PhiSc575NegHel_Combined;
  TH1D* PhiSc625NegHel_Combined;
  TH1D* PhiSc675NegHel_Combined;
  TH1D* PhiSc725NegHel_Combined;
  TH1D* PhiSc775NegHel_Combined;
  TH1D* PhiSc825NegHel_Combined;
  TH1D* PhiSc875NegHel_Combined;
  TH1D* PhiSc275PosHel_Combined;
  TH1D* PhiSc325PosHel_Combined;
  TH1D* PhiSc375PosHel_Combined;
  TH1D* PhiSc425PosHel_Combined;
  TH1D* PhiSc475PosHel_Combined;
  TH1D* PhiSc525PosHel_Combined;
  TH1D* PhiSc575PosHel_Combined;
  TH1D* PhiSc625PosHel_Combined;
  TH1D* PhiSc675PosHel_Combined;
  TH1D* PhiSc725PosHel_Combined;
  TH1D* PhiSc775PosHel_Combined;
  TH1D* PhiSc825PosHel_Combined;
  TH1D* PhiSc875PosHel_Combined;
  int NBins;
  double BinValue;
  double AdjBinValue;
  double TmpBinValue;

  TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/Para/Physics_Total_Para_2_6_3_17.root"); // Open latest Para file
  NPara = KinEp_dE_GoodCut->GetEntries();
  TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/Perp/Physics_Total_Perp_2_6_3_17.root"); // Open latest Perp file
  NPerp = KinEp_dE_GoodCut->GetEntries();
  ScalingFactor = NPara/NPerp;

  TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/Para/Physics_Total_Para_2_6_3_17.root"); // Open latest Para file
  
  NBins = (PhiScNegHel->GetSize())-2;
  PhiScNegHel_Combined = new TH1D("PhiScNegHel_Combined", "Scattetred Proton Phi Distribution in Rotated Frame for -ve Helicity", NBins, -180, 180);

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue = PhiScNegHel -> GetBinContent(i+1);
    PhiScNegHel_Combined -> SetBinContent(i+1, BinValue);

  }

  NBins = (PhiScPosHel->GetSize())-2;
  PhiScPosHel_Combined = new TH1D("PhiScPosHel_Combined", "Scattetred Proton Phi Distribution in Rotated Frame for +ve Helicity", NBins, -180, 180);

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue = PhiScPosHel -> GetBinContent(i+1);
    PhiScPosHel_Combined -> SetBinContent(i+1, BinValue);

  }

  NBins = (Phi_Scattered_275MeV_NegHel->GetSize())-2;
  PhiSc275NegHel_Combined = new TH1D("PhiSc275NegHel_Combined", "Scattetred Proton Phi Distribution in Rotated Frame for Photon Energies of 275pm25MeV for -ve Helicity", NBins, -180, 180);

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue = Phi_Scattered_275MeV_NegHel -> GetBinContent(i+1);
    PhiSc275NegHel_Combined -> SetBinContent(i+1, BinValue);

  }

  NBins = (Phi_Scattered_325MeV_NegHel->GetSize())-2;
  PhiSc325NegHel_Combined = new TH1D("PhiSc325NegHel_Combined", "Scattetred Proton Phi Distribution in Rotated Frame for Photon Energies of 325pm25MeV for -ve Helicity", NBins, -180, 180);

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue = Phi_Scattered_325MeV_NegHel -> GetBinContent(i+1);
    PhiSc325NegHel_Combined -> SetBinContent(i+1, BinValue);

  }

  NBins = (Phi_Scattered_375MeV_NegHel->GetSize())-2;
  PhiSc375NegHel_Combined = new TH1D("PhiSc375NegHel_Combined", "Scattetred Proton Phi Distribution in Rotated Frame for Photon Energies of 375pm25MeV for -ve Helicity", NBins, -180, 180);

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue = Phi_Scattered_375MeV_NegHel -> GetBinContent(i+1);
    PhiSc375NegHel_Combined -> SetBinContent(i+1, BinValue);

  }

  NBins = (Phi_Scattered_425MeV_NegHel->GetSize())-2;
  PhiSc425NegHel_Combined = new TH1D("PhiSc425NegHel_Combined", "Scattetred Proton Phi Distribution in Rotated Frame for Photon Energies of 425pm25MeV for -ve Helicity", NBins, -180, 180);

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue = Phi_Scattered_425MeV_NegHel -> GetBinContent(i+1);
    PhiSc425NegHel_Combined -> SetBinContent(i+1, BinValue);

  }

  NBins = (Phi_Scattered_475MeV_NegHel->GetSize())-2;
  PhiSc475NegHel_Combined = new TH1D("PhiSc475NegHel_Combined", "Scattetred Proton Phi Distribution in Rotated Frame for Photon Energies of 475pm25MeV for -ve Helicity", NBins, -180, 180);

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue = Phi_Scattered_475MeV_NegHel -> GetBinContent(i+1);
    PhiSc475NegHel_Combined -> SetBinContent(i+1, BinValue);

  }

  NBins = (Phi_Scattered_525MeV_NegHel->GetSize())-2;
  PhiSc525NegHel_Combined = new TH1D("PhiSc525NegHel_Combined", "Scattetred Proton Phi Distribution in Rotated Frame for Photon Energies of 525pm25MeV for -ve Helicity", NBins, -180, 180);

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue = Phi_Scattered_525MeV_NegHel-> GetBinContent(i+1);
    PhiSc525NegHel_Combined -> SetBinContent(i+1, BinValue);

  }

  NBins = (Phi_Scattered_575MeV_NegHel->GetSize())-2;
  PhiSc575NegHel_Combined = new TH1D("PhiSc575NegHel_Combined", "Scattetred Proton Phi Distribution in Rotated Frame for Photon Energies of 575pm25MeV for -ve Helicity", NBins, -180, 180);

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue = Phi_Scattered_575MeV_NegHel -> GetBinContent(i+1);
    PhiSc575NegHel_Combined -> SetBinContent(i+1, BinValue);

  }

  NBins = (Phi_Scattered_625MeV_NegHel->GetSize())-2;
  PhiSc625NegHel_Combined = new TH1D("PhiSc625NegHel_Combined", "Scattetred Proton Phi Distribution in Rotated Frame for Photon Energies of 625pm25MeV for -ve Helicity", NBins, -180, 180);

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue =  Phi_Scattered_625MeV_NegHel -> GetBinContent(i+1);
    PhiSc625NegHel_Combined -> SetBinContent(i+1, BinValue);

  }

  NBins = (Phi_Scattered_675MeV_NegHel->GetSize())-2;
  PhiSc675NegHel_Combined = new TH1D("PhiSc675NegHel_Combined", "Scattetred Proton Phi Distribution in Rotated Frame for Photon Energies of 675pm25MeV for -ve Helicity", NBins, -180, 180);

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue =  Phi_Scattered_675MeV_NegHel-> GetBinContent(i+1);
    PhiSc675NegHel_Combined -> SetBinContent(i+1, BinValue);

  }

  NBins = (Phi_Scattered_725MeV_NegHel->GetSize())-2;
  PhiSc725NegHel_Combined = new TH1D("PhiSc725NegHel_Combined", "Scattetred Proton Phi Distribution in Rotated Frame for Photon Energies of 725pm25MeV for -ve Helicity", NBins, -180, 180);

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue = Phi_Scattered_725MeV_NegHel-> GetBinContent(i+1);
    PhiSc725NegHel_Combined -> SetBinContent(i+1, BinValue);

  }

  NBins = (Phi_Scattered_775MeV_NegHel->GetSize())-2;
  PhiSc775NegHel_Combined = new TH1D("PhiSc775NegHel_Combined", "Scattetred Proton Phi Distribution in Rotated Frame for Photon Energies of 775pm25MeV for -ve Helicity", NBins, -180, 180);

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue = Phi_Scattered_775MeV_NegHel -> GetBinContent(i+1);
    PhiSc775NegHel_Combined -> SetBinContent(i+1, BinValue);

  }

  NBins = (Phi_Scattered_825MeV_NegHel->GetSize())-2;
  PhiSc825NegHel_Combined = new TH1D("PhiSc825NegHel_Combined", "Scattetred Proton Phi Distribution in Rotated Frame for Photon Energies of 825pm25MeV for -ve Helicity", NBins, -180, 180);

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue = Phi_Scattered_825MeV_NegHel -> GetBinContent(i+1);
    PhiSc825NegHel_Combined -> SetBinContent(i+1, BinValue);

  }

  NBins = (Phi_Scattered_875MeV_NegHel->GetSize())-2;
  PhiSc875NegHel_Combined = new TH1D("PhiSc875NegHel_Combined", "Scattetred Proton Phi Distribution in Rotated Frame for Photon Energies of 875pm25MeV for -ve Helicity", NBins, -180, 180);

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue =  Phi_Scattered_875MeV_NegHel-> GetBinContent(i+1);
    PhiSc875NegHel_Combined -> SetBinContent(i+1, BinValue);

  }

  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // ALL NEG HEL FOR PARA DONE
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------

  NBins = (Phi_Scattered_275MeV_PosHel->GetSize())-2;
  PhiSc275PosHel_Combined = new TH1D("PhiSc275PosHel_Combined", "Scattetred Proton Phi Distribution in Rotated Frame for Photon Energies of 275pm25MeV for +ve Helicity", NBins, -180, 180);

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue = Phi_Scattered_275MeV_PosHel -> GetBinContent(i+1);
    PhiSc275PosHel_Combined -> SetBinContent(i+1, BinValue);

  }

  NBins = (Phi_Scattered_325MeV_PosHel->GetSize())-2;
  PhiSc325PosHel_Combined = new TH1D("PhiSc325PosHel_Combined", "Scattetred Proton Phi Distribution in Rotated Frame for Photon Energies of 325pm25MeV for +ve Helicity", NBins, -180, 180);

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue = Phi_Scattered_325MeV_PosHel -> GetBinContent(i+1);
    PhiSc325PosHel_Combined -> SetBinContent(i+1, BinValue);

  }

  NBins = (Phi_Scattered_375MeV_PosHel->GetSize())-2;
  PhiSc375PosHel_Combined = new TH1D("PhiSc375PosHel_Combined", "Scattetred Proton Phi Distribution in Rotated Frame for Photon Energies of 375pm25MeV for +ve Helicity", NBins, -180, 180);

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue = Phi_Scattered_375MeV_PosHel -> GetBinContent(i+1);
    PhiSc375PosHel_Combined -> SetBinContent(i+1, BinValue);

  }

  NBins = (Phi_Scattered_425MeV_PosHel->GetSize())-2;
  PhiSc425PosHel_Combined = new TH1D("PhiSc425PosHel_Combined", "Scattetred Proton Phi Distribution in Rotated Frame for Photon Energies of 425pm25MeV for +ve Helicity", NBins, -180, 180);

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue = Phi_Scattered_425MeV_PosHel -> GetBinContent(i+1);
    PhiSc425PosHel_Combined -> SetBinContent(i+1, BinValue);

  }

  NBins = (Phi_Scattered_475MeV_PosHel->GetSize())-2;
  PhiSc475PosHel_Combined = new TH1D("PhiSc475PosHel_Combined", "Scattetred Proton Phi Distribution in Rotated Frame for Photon Energies of 475pm25MeV for +ve Helicity", NBins, -180, 180);

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue = Phi_Scattered_475MeV_PosHel -> GetBinContent(i+1);
    PhiSc475PosHel_Combined -> SetBinContent(i+1, BinValue);

  }

  NBins = (Phi_Scattered_525MeV_PosHel->GetSize())-2;
  PhiSc525PosHel_Combined = new TH1D("PhiSc525PosHel_Combined", "Scattetred Proton Phi Distribution in Rotated Frame for Photon Energies of 525pm25MeV for +ve Helicity", NBins, -180, 180);

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue = Phi_Scattered_525MeV_PosHel-> GetBinContent(i+1);
    PhiSc525PosHel_Combined -> SetBinContent(i+1, BinValue);

  }

  NBins = (Phi_Scattered_575MeV_PosHel->GetSize())-2;
  PhiSc575PosHel_Combined = new TH1D("PhiSc575PosHel_Combined", "Scattetred Proton Phi Distribution in Rotated Frame for Photon Energies of 575pm25MeV for +ve Helicity", NBins, -180, 180);

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue = Phi_Scattered_575MeV_PosHel -> GetBinContent(i+1);
    PhiSc575PosHel_Combined -> SetBinContent(i+1, BinValue);

  }

  NBins = (Phi_Scattered_625MeV_PosHel->GetSize())-2;
  PhiSc625PosHel_Combined = new TH1D("PhiSc625PosHel_Combined", "Scattetred Proton Phi Distribution in Rotated Frame for Photon Energies of 625pm25MeV for +ve Helicity", NBins, -180, 180);

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue =  Phi_Scattered_625MeV_PosHel -> GetBinContent(i+1);
    PhiSc625PosHel_Combined -> SetBinContent(i+1, BinValue);

  }

  NBins = (Phi_Scattered_675MeV_PosHel->GetSize())-2;
  PhiSc675PosHel_Combined = new TH1D("PhiSc675PosHel_Combined", "Scattetred Proton Phi Distribution in Rotated Frame for Photon Energies of 675pm25MeV for +ve Helicity", NBins, -180, 180);

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue =  Phi_Scattered_675MeV_PosHel -> GetBinContent(i+1);
    PhiSc675PosHel_Combined -> SetBinContent(i+1, BinValue);

  }

  NBins = (Phi_Scattered_725MeV_PosHel->GetSize())-2;
  PhiSc725PosHel_Combined = new TH1D("PhiSc725PosHel_Combined", "Scattetred Proton Phi Distribution in Rotated Frame for Photon Energies of 725pm25MeV for +ve Helicity", NBins, -180, 180);

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue = Phi_Scattered_725MeV_PosHel -> GetBinContent(i+1);
    PhiSc725PosHel_Combined -> SetBinContent(i+1, BinValue);

  }

  NBins = (Phi_Scattered_775MeV_PosHel->GetSize())-2;
  PhiSc775PosHel_Combined = new TH1D("PhiSc775PosHel_Combined", "Scattetred Proton Phi Distribution in Rotated Frame for Photon Energies of 775pm25MeV for +ve Helicity", NBins, -180, 180);

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue = Phi_Scattered_775MeV_PosHel -> GetBinContent(i+1);
    PhiSc775PosHel_Combined -> SetBinContent(i+1, BinValue);

  }

  NBins = (Phi_Scattered_825MeV_PosHel->GetSize())-2;
  PhiSc825PosHel_Combined = new TH1D("PhiSc825PosHel_Combined", "Scattetred Proton Phi Distribution in Rotated Frame for Photon Energies of 825pm25MeV for +ve Helicity", NBins, -180, 180);

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue = Phi_Scattered_825MeV_PosHel -> GetBinContent(i+1);
    PhiSc825PosHel_Combined -> SetBinContent(i+1, BinValue);

  }

  NBins = (Phi_Scattered_875MeV_PosHel->GetSize())-2;
  PhiSc875PosHel_Combined = new TH1D("PhiSc875PosHel_Combined", "Scattetred Proton Phi Distribution in Rotated Frame for Photon Energies of 875pm25MeV for +ve Helicity", NBins, -180, 180);

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue =  Phi_Scattered_875MeV_PosHel-> GetBinContent(i+1);
    PhiSc875PosHel_Combined -> SetBinContent(i+1, BinValue);

  }

  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // ALL PARA ONES FILLED NOW DO PERP
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------

  TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/Perp/Physics_Total_Perp_2_6_3_17.root"); // Open latest Perp file
  
  NBins = (PhiScNegHel->GetSize())-2;

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue = PhiScNegHel -> GetBinContent(i+1);
    AdjBinValue = BinValue*ScalingFactor;
    TmpBinValue = PhiScNegHel_Combined -> GetBinContent(i+1);
    PhiScNegHel_Combined -> SetBinContent(i+1, (AdjBinValue + TmpBinValue));

  }

  NBins = (PhiScPosHel->GetSize())-2;

  for (Int_t i = 0; i < NBins; i++){

    BinValue = PhiScPosHel -> GetBinContent(i+1);
    AdjBinValue = BinValue*ScalingFactor;
    TmpBinValue = PhiScPosHel_Combined -> GetBinContent(i+1);
    PhiScPosHel_Combined -> SetBinContent(i+1, (AdjBinValue + TmpBinValue));

  }

  NBins = (Phi_Scattered_275MeV_NegHel->GetSize())-2;

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue = Phi_Scattered_275MeV_NegHel -> GetBinContent(i+1);
    AdjBinValue = BinValue*ScalingFactor;
    TmpBinValue = PhiScNegHel_Combined -> GetBinContent(i+1);
    PhiSc275NegHel_Combined -> SetBinContent(i+1, (AdjBinValue + TmpBinValue));

  }
 
  NBins = (Phi_Scattered_325MeV_NegHel->GetSize())-2;

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue = Phi_Scattered_325MeV_NegHel -> GetBinContent(i+1);
    AdjBinValue = BinValue*ScalingFactor;
    TmpBinValue = PhiScNegHel_Combined -> GetBinContent(i+1);
    PhiSc325NegHel_Combined -> SetBinContent(i+1, (AdjBinValue + TmpBinValue));

  }

  NBins = (Phi_Scattered_375MeV_NegHel->GetSize())-2;

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue = Phi_Scattered_375MeV_NegHel -> GetBinContent(i+1);
    AdjBinValue = BinValue*ScalingFactor;
    TmpBinValue = PhiScNegHel_Combined -> GetBinContent(i+1);
    PhiSc375NegHel_Combined -> SetBinContent(i+1, (AdjBinValue + TmpBinValue));

  }

  NBins = (Phi_Scattered_425MeV_NegHel->GetSize())-2;

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue = Phi_Scattered_425MeV_NegHel -> GetBinContent(i+1);
    AdjBinValue = BinValue*ScalingFactor;
    TmpBinValue = PhiScNegHel_Combined -> GetBinContent(i+1);
    PhiSc425NegHel_Combined -> SetBinContent(i+1, (AdjBinValue + TmpBinValue));

  }

  NBins = (Phi_Scattered_475MeV_NegHel->GetSize())-2;

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue = Phi_Scattered_475MeV_NegHel -> GetBinContent(i+1);
    AdjBinValue = BinValue*ScalingFactor;
    TmpBinValue = PhiScNegHel_Combined -> GetBinContent(i+1);
    PhiSc475NegHel_Combined -> SetBinContent(i+1, (AdjBinValue + TmpBinValue));

  }

  NBins = (Phi_Scattered_525MeV_NegHel->GetSize())-2;

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue = Phi_Scattered_525MeV_NegHel -> GetBinContent(i+1);
    AdjBinValue = BinValue*ScalingFactor;
    TmpBinValue = PhiScNegHel_Combined -> GetBinContent(i+1);
    PhiSc525NegHel_Combined -> SetBinContent(i+1, (AdjBinValue + TmpBinValue));

  }

  NBins = (Phi_Scattered_575MeV_NegHel->GetSize())-2;

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue = Phi_Scattered_575MeV_NegHel -> GetBinContent(i+1);
    AdjBinValue = BinValue*ScalingFactor;
    TmpBinValue = PhiScNegHel_Combined -> GetBinContent(i+1);
    PhiSc575NegHel_Combined -> SetBinContent(i+1, (AdjBinValue + TmpBinValue));

  }

  NBins = (Phi_Scattered_625MeV_NegHel->GetSize())-2;

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue = Phi_Scattered_625MeV_NegHel -> GetBinContent(i+1);
    AdjBinValue = BinValue*ScalingFactor;
    TmpBinValue = PhiScNegHel_Combined -> GetBinContent(i+1);
    PhiSc625NegHel_Combined -> SetBinContent(i+1, (AdjBinValue + TmpBinValue));

  }

  NBins = (Phi_Scattered_675MeV_NegHel->GetSize())-2;

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue = Phi_Scattered_675MeV_NegHel -> GetBinContent(i+1);
    AdjBinValue = BinValue*ScalingFactor;
    TmpBinValue = PhiScNegHel_Combined -> GetBinContent(i+1);
    PhiSc675NegHel_Combined -> SetBinContent(i+1, (AdjBinValue + TmpBinValue));

  }

  NBins = (Phi_Scattered_725MeV_NegHel->GetSize())-2;

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue = Phi_Scattered_725MeV_NegHel -> GetBinContent(i+1);
    AdjBinValue = BinValue*ScalingFactor;
    TmpBinValue = PhiScNegHel_Combined -> GetBinContent(i+1);
    PhiSc725NegHel_Combined -> SetBinContent(i+1, (AdjBinValue + TmpBinValue));

  }

  NBins = (Phi_Scattered_775MeV_NegHel->GetSize())-2;

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue = Phi_Scattered_775MeV_NegHel -> GetBinContent(i+1);
    AdjBinValue = BinValue*ScalingFactor;
    TmpBinValue = PhiScNegHel_Combined -> GetBinContent(i+1);
    PhiSc775NegHel_Combined -> SetBinContent(i+1, (AdjBinValue + TmpBinValue));

  }

  NBins = (Phi_Scattered_825MeV_NegHel->GetSize())-2;

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue = Phi_Scattered_825MeV_NegHel -> GetBinContent(i+1);
    AdjBinValue = BinValue*ScalingFactor;
    TmpBinValue = PhiScNegHel_Combined -> GetBinContent(i+1);
    PhiSc825NegHel_Combined -> SetBinContent(i+1, (AdjBinValue + TmpBinValue));

  }

  NBins = (Phi_Scattered_875MeV_NegHel->GetSize())-2;

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue = Phi_Scattered_875MeV_NegHel -> GetBinContent(i+1);
    AdjBinValue = BinValue*ScalingFactor;
    TmpBinValue = PhiScNegHel_Combined -> GetBinContent(i+1);
    PhiSc875NegHel_Combined -> SetBinContent(i+1, (AdjBinValue + TmpBinValue));

  }

  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // ALL NEG HEL FOR PERP DONE
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------

 NBins = (Phi_Scattered_275MeV_PosHel->GetSize())-2;

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue = Phi_Scattered_275MeV_PosHel -> GetBinContent(i+1);
    AdjBinValue = BinValue*ScalingFactor;
    TmpBinValue = PhiScPosHel_Combined -> GetBinContent(i+1);
    PhiSc275PosHel_Combined -> SetBinContent(i+1, (AdjBinValue + TmpBinValue));

  }
 
  NBins = (Phi_Scattered_325MeV_PosHel->GetSize())-2;

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue = Phi_Scattered_325MeV_PosHel -> GetBinContent(i+1);
    AdjBinValue = BinValue*ScalingFactor;
    TmpBinValue = PhiScPosHel_Combined -> GetBinContent(i+1);
    PhiSc325PosHel_Combined -> SetBinContent(i+1, (AdjBinValue + TmpBinValue));

  }

  NBins = (Phi_Scattered_375MeV_PosHel->GetSize())-2;

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue = Phi_Scattered_375MeV_PosHel -> GetBinContent(i+1);
    AdjBinValue = BinValue*ScalingFactor;
    TmpBinValue = PhiScPosHel_Combined -> GetBinContent(i+1);
    PhiSc375PosHel_Combined -> SetBinContent(i+1, (AdjBinValue + TmpBinValue));

  }

  NBins = (Phi_Scattered_425MeV_PosHel->GetSize())-2;

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue = Phi_Scattered_425MeV_PosHel -> GetBinContent(i+1);
    AdjBinValue = BinValue*ScalingFactor;
    TmpBinValue = PhiScPosHel_Combined -> GetBinContent(i+1);
    PhiSc425PosHel_Combined -> SetBinContent(i+1, (AdjBinValue + TmpBinValue));

  }

  NBins = (Phi_Scattered_475MeV_PosHel->GetSize())-2;

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue = Phi_Scattered_475MeV_PosHel -> GetBinContent(i+1);
    AdjBinValue = BinValue*ScalingFactor;
    TmpBinValue = PhiScPosHel_Combined -> GetBinContent(i+1);
    PhiSc475PosHel_Combined -> SetBinContent(i+1, (AdjBinValue + TmpBinValue));

  }

  NBins = (Phi_Scattered_525MeV_PosHel->GetSize())-2;

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue = Phi_Scattered_525MeV_PosHel -> GetBinContent(i+1);
    AdjBinValue = BinValue*ScalingFactor;
    TmpBinValue = PhiScPosHel_Combined -> GetBinContent(i+1);
    PhiSc525PosHel_Combined -> SetBinContent(i+1, (AdjBinValue + TmpBinValue));

  }

  NBins = (Phi_Scattered_575MeV_PosHel->GetSize())-2;

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue = Phi_Scattered_575MeV_PosHel -> GetBinContent(i+1);
    AdjBinValue = BinValue*ScalingFactor;
    TmpBinValue = PhiScPosHel_Combined -> GetBinContent(i+1);
    PhiSc575PosHel_Combined -> SetBinContent(i+1, (AdjBinValue + TmpBinValue));

  }

  NBins = (Phi_Scattered_625MeV_PosHel->GetSize())-2;

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue = Phi_Scattered_625MeV_PosHel -> GetBinContent(i+1);
    AdjBinValue = BinValue*ScalingFactor;
    TmpBinValue = PhiScPosHel_Combined -> GetBinContent(i+1);
    PhiSc625PosHel_Combined -> SetBinContent(i+1, (AdjBinValue + TmpBinValue));

  }

  NBins = (Phi_Scattered_675MeV_PosHel->GetSize())-2;

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue = Phi_Scattered_675MeV_PosHel -> GetBinContent(i+1);
    AdjBinValue = BinValue*ScalingFactor;
    TmpBinValue = PhiScPosHel_Combined -> GetBinContent(i+1);
    PhiSc675PosHel_Combined -> SetBinContent(i+1, (AdjBinValue + TmpBinValue));

  }

  NBins = (Phi_Scattered_725MeV_PosHel->GetSize())-2;

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue = Phi_Scattered_725MeV_PosHel -> GetBinContent(i+1);
    AdjBinValue = BinValue*ScalingFactor;
    TmpBinValue = PhiScPosHel_Combined -> GetBinContent(i+1);
    PhiSc725PosHel_Combined -> SetBinContent(i+1, (AdjBinValue + TmpBinValue));

  }

  NBins = (Phi_Scattered_775MeV_PosHel->GetSize())-2;

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue = Phi_Scattered_775MeV_PosHel -> GetBinContent(i+1);
    AdjBinValue = BinValue*ScalingFactor;
    TmpBinValue = PhiScPosHel_Combined -> GetBinContent(i+1);
    PhiSc775PosHel_Combined -> SetBinContent(i+1, (AdjBinValue + TmpBinValue));

  }

  NBins = (Phi_Scattered_825MeV_PosHel->GetSize())-2;

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue = Phi_Scattered_825MeV_PosHel -> GetBinContent(i+1);
    AdjBinValue = BinValue*ScalingFactor;
    TmpBinValue = PhiScPosHel_Combined -> GetBinContent(i+1);
    PhiSc825PosHel_Combined -> SetBinContent(i+1, (AdjBinValue + TmpBinValue));

  }

  NBins = (Phi_Scattered_875MeV_PosHel->GetSize())-2;

  for (Int_t i = 0; i < NBins; i++){
    
    BinValue = Phi_Scattered_875MeV_PosHel -> GetBinContent(i+1);
    AdjBinValue = BinValue*ScalingFactor;
    TmpBinValue = PhiScPosHel_Combined -> GetBinContent(i+1);
    PhiSc875PosHel_Combined -> SetBinContent(i+1, (AdjBinValue + TmpBinValue));

  }

  // Define new file to store fit parameters
  TFile f1("ParaPerp_Total_2_Combined.root", "RECREATE");

   PhiScNegHel_Combined->Write();
   PhiScPosHel_Combined->Write();
   PhiSc275NegHel_Combined->Write();
   PhiSc325NegHel_Combined->Write();
   PhiSc375NegHel_Combined->Write();
   PhiSc425NegHel_Combined->Write();
   PhiSc475NegHel_Combined->Write();
   PhiSc525NegHel_Combined->Write();
   PhiSc575NegHel_Combined->Write();
   PhiSc625NegHel_Combined->Write();
   PhiSc675NegHel_Combined->Write();
   PhiSc725NegHel_Combined->Write();
   PhiSc775NegHel_Combined->Write();
   PhiSc825NegHel_Combined->Write();
   PhiSc875NegHel_Combined->Write();
   PhiSc275PosHel_Combined->Write();
   PhiSc325PosHel_Combined->Write();
   PhiSc375PosHel_Combined->Write();
   PhiSc425PosHel_Combined->Write();
   PhiSc475PosHel_Combined->Write();
   PhiSc525PosHel_Combined->Write();
   PhiSc575PosHel_Combined->Write();
   PhiSc625PosHel_Combined->Write();
   PhiSc675PosHel_Combined->Write();
   PhiSc725PosHel_Combined->Write();
   PhiSc775PosHel_Combined->Write();
   PhiSc825PosHel_Combined->Write();
   PhiSc875PosHel_Combined->Write();
  
  f1.Write();

}
