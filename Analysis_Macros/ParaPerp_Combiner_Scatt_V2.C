void ParaPerp_Combiner_Scatt_V2(){

    char HistName[60];
    char NewHistName[60];
    TH1F* PhipPara[10][5];
    TH1F* PhipPerp[10][5];

    TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/Para/Physics_Total_Para_35_17_5_18.root"); // Open latest Para file

    TH1D* Eg_Para = (TH1D*)f->Get("Eg")->Clone();
    Eg_Para->SetName("Eg_Para");

    for(Int_t i = 0; i < 10; i++){ // Energy
        for(Int_t j = 0; j < 5; j++){ // Theta
            sprintf(HistName, "Phip_%iMeVCM%i", 415+(i*20) , j+1);
            sprintf(NewHistName, "Phip_%iMeVCM%i_Para", 415+(i*20) , j+1);
            PhipPara[i][j] = ((TH1F*)f->Get(HistName));
            PhipPara[i][j]->SetName(NewHistName);
        }
    }


    /////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////
    ////////////////////////// Para Done ////////////////////////////
    /////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////

    TFile *f1 = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/Perp/Physics_Total_Perp_35_17_5_18.root"); // Open latest Para file

    TH1D* Eg_Perp = (TH1D*)f1->Get("Eg")->Clone();
    Eg_Perp->SetName("Eg_Perp");

    for(Int_t i = 0; i < 10; i++){ // Energy
        for(Int_t j = 0; j < 5; j++){ // Theta
            sprintf(HistName, "Phip_%iMeVCM%i", 415+(i*20) , j+1);
            sprintf(NewHistName, "Phip_%iMeVCM%i_Perp", 415+(i*20) , j+1);
            PhipPerp[i][j] = ((TH1F*)f1->Get(HistName));
            PhipPerp[i][j]->SetName(NewHistName);
        }
    }

    TFile f2("ParaPerp_S35_Combined.root", "RECREATE");

    Eg_Para->Write();
    Eg_Perp->Write();

    for(Int_t i = 0; i < 10; i++){ // Energy
        for(Int_t j = 0; j < 5; j++){ // Theta
            PhipPara[i][j]->Write();
            PhipPerp[i][j]->Write();
        }
    }

    f2.Write();

}
