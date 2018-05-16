#include "./includes_SigmaAsymm_Scatt_V2.h"

void SigmaAsymm_Scatt_V2(){

    TF1 *CosFunc = new TF1("CosFit", "[0]*cos((2*x*TMath::DegToRad())+acos(0))");
    CosFunc->SetParNames("Amplitude");

    TH1F* AsymmHists[10][5];

    char ParaHistName[60];
    char PerpHistName[60];
    char AsymmHistName[60];

    TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/ParaPerp_S34_Combined.root"); // Open the latest PTotal combined file to load histograms from
    NPara = Eg_Para->GetEntries();
    NPerp = Eg_Perp->GetEntries();
    ScaleFactor = NPara/NPerp;
    ScaleFactorErr = sqrt( (NPara/((TMath::Power(NPerp,2)))) + (((TMath::Power(NPara,2)))/(TMath::Power(NPerp,3))) ); // Error Propagation of above formula, see notebook 5

    for(Int_t i = 0; i < 10; i++){ // Energy
        for(Int_t j = 0; j < 5; j++){ // Theta
            sprintf(ParaHistName, "Phip_%iMeVCM%i_Para", 410+(i*20) , j+1);
            sprintf(PerpHistName, "Phip_%iMeVCM%i_Perp", 410+(i*20) , j+1);
            sprintf(AsymmHistName, "Sigma_%iMeV_CM%i_Hist", 410+(i*20) , j+1);
            AsymmHists[i][j] = (TH1F*) (((TH1F*)f->Get(ParaHistName))->GetAsymmetry(((TH1F*)f->Get(PerpHistName)), ScaleFactor, ScaleFactorErr)));
            AsymmHists[i][j]->SetName(AsymmHistName);
            AsymmHists[i][j]-> Fit("CosFit", "Q");
            cout << "NDOF " << CosFit->GetNDF() << "   " << "Chi2 " << CosFit->GetChisquare() << "   " << "Chi2/DoF " << CosFit->GetChisquare()/CosFit->GetNDF() << endl;
            pCosAmp[i][j] = CosFit->GetParameter(0);
            pCosAmpErr[i][j] = CosFit->GetParError(0);
        }
    }

    TFile f1("ParaPerpAsymm_S34.root", "RECREATE");

    //Define new tree to store parameters in
    TTree* tree = new TTree("Parameter_Values", "Tree_of_Values");

    // Define branches to store parameters, (Branch Name, Variable, Type of Variable)
    tree->Branch("pCosAmp410", &pCosA410, "pCosA410/D");
    tree->Branch("pCosAmpErr410", &pCosAErr410, "pCosAErr410/D");
    tree->Branch("pCosAmp430", &pCosA430, "pCosA430/D");
    tree->Branch("pCosAmpErr430", &pCosAErr430, "pCosAErr430/D");
    tree->Branch("pCosAmp450", &pCosA450, "pCosA450/D");
    tree->Branch("pCosAmpErr450", &pCosAErr450, "pCosAErr450/D");
    tree->Branch("pCosAmp470", &pCosA470, "pCosA470/D");
    tree->Branch("pCosAmpErr470", &pCosAErr470, "pCosAErr470/D");
    tree->Branch("pCosAmp490", &pCosA490, "pCosA490/D");
    tree->Branch("pCosAmpErr490", &pCosAErr490, "pCosAErr490/D");
    tree->Branch("pCosAmp510", &pCosA510, "pCosA510/D");
    tree->Branch("pCosAmpErr510", &pCosAErr510, "pCosAErr510/D");
    tree->Branch("pCosAmp530", &pCosA530, "pCosA530/D");
    tree->Branch("pCosAmpErr530", &pCosAErr530, "pCosAErr530/D");
    tree->Branch("pCosAmp550", &pCosA550, "pCosA550/D");
    tree->Branch("pCosAmpErr550", &pCosAErr550, "pCosAErr550/D");
    tree->Branch("pCosAmp570", &pCosA570, "pCosA570/D");
    tree->Branch("pCosAmpErr570", &pCosAErr570, "pCosAErr570/D");
    tree->Branch("pCosAmp590", &pCosA590, "pCosA590/D");
    tree->Branch("pCosAmpErr590", &pCosAErr590, "pCosAErr590/D");

    // Fill branches (and hence tree) with corresponding parameters from above
    for (Int_t m = 0; m < 5; m++){
        pCosA410 = pCosAmp[0][m];
        pCosAErr410 = pCosAmpErr[0][m];
        pCosA430 = pCosAmp[1][m];
        pCosAErr430 = pCosAmpErr[1][m];
        pCosA450 = pCosAmp[2][m];
        pCosAErr450 = pCosAmpErr[2][m];
        pCosA470 = pCosAmp[3][m];
        pCosAErr470 = pCosAmpErr[3][m];
        pCosA490 = pCosAmp[4][m];
        pCosAErr490 = pCosAmpErr[4][m];
        pCosA510 = pCosAmp[5][m];
        pCosAErr510 = pCosAmpErr[5][m];
        pCosA530 = pCosAmp[6][m];
        pCosAErr530 = pCosAmpErr[6][m];
        pCosA550 = pCosAmp[7][m];
        pCosAErr550 = pCosAmpErr[7][m];
        pCosA570 = pCosAmp[8][m];
        pCosAErr570 = pCosAmpErr[8][m];
        pCosA590 = pCosAmp[9][m];
        pCosAErr590 = pCosAmpErr[9][m];

        tree->Fill();
    }

    for(Int_t i = 0; i < 10; i++){ // Energy
        for(Int_t j = 0; j < 5; j++){ // Theta
            AsymmHists[i][j]->Write();
        }
    }

    f1.Write();

}
