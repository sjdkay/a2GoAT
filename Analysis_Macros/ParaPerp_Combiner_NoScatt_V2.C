#include "./includes_ParaPerpCombiner_NoScatt_V2.h"

void ParaPerp_Combiner_NoScatt_V2(){

    char HistName[60];
    char NewHistName[60];
    TH1F* PhipPara[21][20];
    TH1F* PhipPerp[21][20];

    TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/Para/NoScatt/Physics_Total_Para_NoScatt_14_9_4_18.root"); // Open latest Para file

    TH1D* Eg_Para = (TH1D*)f->Get("Eg2")->Clone();
    Eg_Para->SetName("Eg_Para");

    for(Int_t i = 0; i < 21; i++){ // Energy
        for(Int_t j = 0; j < 20; j++){ // Theta
            sprintf(HistName, "Phip_%iMeVCM%i", 415+(i*10) , j+1);
            sprintf(NewHistName, "Phip_%iMeVCM%i_Para", 415+(i*10) , j+1);
            PhipPara[i][j] = ((TH1F*)f->Get(HistName));
            PhipPara[i][j]->SetName(NewHistName);
        }
    }


    /////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////
    ////////////////////////// Para Done ////////////////////////////
    /////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////

    TFile *f1 = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/Perp/NoScatt/Physics_Total_Perp_NoScatt_14_9_4_18.root"); // Open latest Para file

    TH1D* Eg_Perp = (TH1D*)f1->Get("Eg2")->Clone();
    Eg_Perp->SetName("Eg_Perp");

    for(Int_t i = 0; i < 21; i++){ // Energy
        for(Int_t j = 0; j < 20; j++){ // Theta
            sprintf(HistName, "Phip_%iMeVCM%i", 415+(i*10) , j+1);
            sprintf(NewHistName, "Phip_%iMeVCM%i_Perp", 415+(i*10) , j+1);
            PhipPerp[i][j] = ((TH1F*)f1->Get(HistName));
            PhipPerp[i][j]->SetName(NewHistName);
        }
    }

    TFile f2("ParaPerp_NS14_Combined_V2.root", "RECREATE");

    Eg_Para->Write();
    Eg_Perp->Write();

    for(Int_t i = 0; i < 21; i++){ // Energy
        for(Int_t j = 0; j < 20; j++){ // Theta
            PhipPara[i][j]->Write();
            PhipPerp[i][j]->Write();
        }
    }

    f2.Write();

}
