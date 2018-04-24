#include "./includes.h"

void CxBGFitting() {

    TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/Amo/Physics_Total_109_12_4_18.root"); // Open the latest PTotal file to load histograms from

    TH1D* MM_Eg_CM_Tot[7][5];
    char Histname[100];

    for (Int_t i = 0; i < 7; i++){
        for(Int_t j = 0; j < 5; j++){
            sprintf(Histname, "MM_Eg%i_CM%i_Tot", 265+(i*70), j+1);
            MM_Eg_CM_Tot[i][j] = (TH1D*)f->Get(Histname);
        }
    }

    // Want to fit to peaks and add a background function in too
    // Looks like one peak at roughly nucleon mass and another at Nucleon + pion

    TFile *f1 = new TFile("CxBGFitting_109_V1.root", "RECREATE"); // Open the latest PTotal file to load histograms from

    for (Int_t i = 0; i < 7; i++){
        for(Int_t j = 0; j < 5; j++){
            MM_Eg_CM_Tot[i][j]->Write();
        }
    }

    f1->Write();

}
