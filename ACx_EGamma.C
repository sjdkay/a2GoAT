#include "./includes.h"

void ACx_EGamma(){

  TFile *f1= TFile::Open("/scratch/Mainz_Software/a2GoAT/SinPhiFitValues.root");
  TTree *t1 = (TTree*)f1->Get("Parameter_Values");

  Double_t Parameters[10][6];
  Double_t Par1, Par1Err, Par2, Par2Err, Par3, Par3Err;

  t1->SetBranchAddress("Par1", &Par1);
  t1->SetBranchAddress("Par1Err", &Par1Err);
  t1->SetBranchAddress("Par2", &Par2);
  t1->SetBranchAddress("Par2Err", &Par2Err);
  t1->SetBranchAddress("Par3", &Par3);
  t1->SetBranchAddress("Par3Err", &Par3Err);

  for (Int_t k = 0; k < 10; k++){

    Parameter_Values->GetEntry(k);
    Parameters[k][0] = Par1;
    Parameters[k][1] = Par1Err;
    Parameters[k][2] = Par2;
    Parameters[k][3] = Par2Err;
    Parameters[k][4] = Par3;
    Parameters[k][5] = Par3Err;

  }
