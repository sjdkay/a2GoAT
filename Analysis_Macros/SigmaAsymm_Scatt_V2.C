#include "./includes_SigmaAsymm_Scatt_V2.h"

void SigmaAsymm_Scatt_V2(){

    TF1 *CosFunc = new TF1("CosFit", "[0]*cos((2*x*TMath::DegToRad())+acos(0))");
    CosFunc->SetParNames("Amplitude");

    TH1F* AsymmHists[10][5];

    char ParaHistName[60];
    char PerpHistName[60];
    char AsymmHistName[60];

    TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/ParaPerp_S35_Combined.root"); // Open the latest PTotal combined file to load histograms from
    NPara = Eg_Para->GetEntries();
    NPerp = Eg_Perp->GetEntries();
    ScaleFactor = NPara/NPerp;
    ScaleFactorErr = sqrt( (NPara/((TMath::Power(NPerp,2)))) + (((TMath::Power(NPara,2)))/(TMath::Power(NPerp,3))) ); // Error Propagation of above formula, see notebook 5

    for(Int_t i = 0; i < 10; i++){ // Energy
        for(Int_t j = 0; j < 5; j++){ // Theta
            sprintf(ParaHistName, "Phip_%iMeVCM%i_Para", 415+(i*20) , j+1);
            sprintf(PerpHistName, "Phip_%iMeVCM%i_Perp", 415+(i*20) , j+1);
            sprintf(AsymmHistName, "Sigma_%iMeV_CM%i_Hist", 415+(i*20) , j+1);
            AsymmHists[i][j] = (TH1F*) (((TH1F*)f->Get(ParaHistName))->GetAsymmetry(((TH1F*)f->Get(PerpHistName)), ScaleFactor, ScaleFactorErr)));
            AsymmHists[i][j]->SetName(AsymmHistName);
            AsymmHists[i][j]-> Fit("CosFit", "Q");
            cout << "NDOF " << CosFit->GetNDF() << "   " << "Chi2 " << CosFit->GetChisquare() << "   " << "Chi2/DoF " << CosFit->GetChisquare()/CosFit->GetNDF() << endl;
            pCosAmp[i][j] = CosFit->GetParameter(0);
            pCosAmpErr[i][j] = CosFit->GetParError(0);
        }
    }

    TFile f1("ParaPerpAsymm_S35.root", "RECREATE");

    //Define new tree to store parameters in
    TTree* tree = new TTree("Parameter_Values", "Tree_of_Values");

    // Define branches to store parameters, (Branch Name, Variable, Type of Variable)
    tree->Branch("pCosAmp415", &pCosA415, "pCosA415/D");
    tree->Branch("pCosAmpErr415", &pCosAErr415, "pCosAErr415/D");
    tree->Branch("pCosAmp435", &pCosA435, "pCosA435/D");
    tree->Branch("pCosAmpErr435", &pCosAErr435, "pCosAErr435/D");
    tree->Branch("pCosAmp455", &pCosA455, "pCosA455/D");
    tree->Branch("pCosAmpErr455", &pCosAErr455, "pCosAErr455/D");
    tree->Branch("pCosAmp475", &pCosA475, "pCosA475/D");
    tree->Branch("pCosAmpErr475", &pCosAErr475, "pCosAErr475/D");
    tree->Branch("pCosAmp495", &pCosA495, "pCosA495/D");
    tree->Branch("pCosAmpErr495", &pCosAErr495, "pCosAErr495/D");
    tree->Branch("pCosAmp515", &pCosA515, "pCosA515/D");
    tree->Branch("pCosAmpErr515", &pCosAErr515, "pCosAErr515/D");
    tree->Branch("pCosAmp535", &pCosA535, "pCosA535/D");
    tree->Branch("pCosAmpErr535", &pCosAErr535, "pCosAErr535/D");
    tree->Branch("pCosAmp555", &pCosA555, "pCosA555/D");
    tree->Branch("pCosAmpErr555", &pCosAErr555, "pCosAErr555/D");
    tree->Branch("pCosAmp575", &pCosA575, "pCosA575/D");
    tree->Branch("pCosAmpErr575", &pCosAErr575, "pCosAErr575/D");
    tree->Branch("pCosAmp595", &pCosA595, "pCosA595/D");
    tree->Branch("pCosAmpErr595", &pCosAErr595, "pCosAErr595/D");

    // Fill branches (and hence tree) with corresponding parameters from above
    for (Int_t m = 0; m < 5; m++){
        pCosA415 = pCosAmp[0][m];
        pCosAErr415 = pCosAmpErr[0][m];
        pCosA435 = pCosAmp[1][m];
        pCosAErr435 = pCosAmpErr[1][m];
        pCosA455 = pCosAmp[2][m];
        pCosAErr455 = pCosAmpErr[2][m];
        pCosA475 = pCosAmp[3][m];
        pCosAErr475 = pCosAmpErr[3][m];
        pCosA495 = pCosAmp[4][m];
        pCosAErr495 = pCosAmpErr[4][m];
        pCosA515 = pCosAmp[5][m];
        pCosAErr515 = pCosAmpErr[5][m];
        pCosA535 = pCosAmp[6][m];
        pCosAErr535 = pCosAmpErr[6][m];
        pCosA555 = pCosAmp[7][m];
        pCosAErr555 = pCosAmpErr[7][m];
        pCosA575 = pCosAmp[8][m];
        pCosAErr575 = pCosAmpErr[8][m];
        pCosA595 = pCosAmp[9][m];
        pCosAErr595 = pCosAmpErr[9][m];

        tree->Fill();
    }

    for(Int_t i = 0; i < 10; i++){ // Energy
        for(Int_t j = 0; j < 5; j++){ // Theta
            AsymmHists[i][j]->Write();
        }
    }

    f1.Write();

}
