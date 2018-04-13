#include "./includes_SigmaAsymm_NoScatt_V2.h"

void SigmaAsymm_NoScatt_V2(){

    TF1 *CosFunc = new TF1("CosFit", "[0]*cos((2*x*TMath::DegToRad())+acos(0))");
    CosFunc->SetParNames("Amplitude");

    TH1F* AsymmHists[21][20];

    char ParaHistName[60];
    char PerpHistName[60];
    char AsymmHistName[60];

    TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/ParaPerp_NS14_Combined_V2.root"); // Open the latest PTotal combined file to load histograms from
    NPara = Eg_Para->GetEntries();
    NPerp = Eg_Perp->GetEntries();
    ScaleFactor = NPara/NPerp;
    ScaleFactorErr = sqrt( (NPara/((TMath::Power(NPerp,2)))) + (((TMath::Power(NPara,2)))/(TMath::Power(NPerp,3))) ); // Error Propagation of above formula, see notebook 5


    for(Int_t i = 0; i < 21; i++){ // Energy
        for(Int_t j = 0; j < 20; j++){ // Theta
            sprintf(ParaHistName, "Phip_%iMeVCM%i_Para", 415+(i*10) , j+1);
            sprintf(PerpHistName, "Phip_%iMeVCM%i_Perp", 415+(i*10) , j+1);
            sprintf(AsymmHistName, "Sigma_%iMeV_CM%i_Hist", 415+(i*10) , j+1);
            AsymmHists[i][j] = (TH1F*) (((TH1F*)f->Get(ParaHistName))->GetAsymmetry(((TH1F*)f->Get(PerpHistName)), ScaleFactor, ScaleFactorErr)));
            AsymmHists[i][j]->SetName(AsymmHistName);
            AsymmHists[i][j]-> Fit("CosFit", "Q");
            pCosAmp[i][j] = CosFit->GetParameter(0);
            pCosAmpErr[i][j] = CosFit->GetParError(0);
        }
    }

    TFile f1("ParaPerpAsymm_NS14_V2.root", "RECREATE");

    //Define new tree to store parameters in
    TTree* tree = new TTree("Parameter_Values", "Tree_of_Values");

    // Define branches to store parameters, (Branch Name, Variable, Type of Variable)
    tree->Branch("pCosAmp415", &pCosA415, "pCosA415/D");
    tree->Branch("pCosAmpErr415", &pCosAErr415, "pCosAErr415/D");
    tree->Branch("pCosAmp425", &pCosA425, "pCosA425/D");
    tree->Branch("pCosAmpErr425", &pCosAErr425, "pCosAErr425/D");
    tree->Branch("pCosAmp435", &pCosA435, "pCosA435/D");
    tree->Branch("pCosAmpErr435", &pCosAErr435, "pCosAErr435/D");
    tree->Branch("pCosAmp445", &pCosA445, "pCosA445/D");
    tree->Branch("pCosAmpErr445", &pCosAErr445, "pCosAErr445/D");
    tree->Branch("pCosAmp455", &pCosA455, "pCosA455/D");
    tree->Branch("pCosAmpErr455", &pCosAErr455, "pCosAErr455/D");
    tree->Branch("pCosAmp465", &pCosA465, "pCosA465/D");
    tree->Branch("pCosAmpErr465", &pCosAErr465, "pCosAErr465/D");
    tree->Branch("pCosAmp475", &pCosA475, "pCosA475/D");
    tree->Branch("pCosAmpErr475", &pCosAErr475, "pCosAErr475/D");
    tree->Branch("pCosAmp485", &pCosA485, "pCosA485/D");
    tree->Branch("pCosAmpErr485", &pCosAErr485, "pCosAErr485/D");
    tree->Branch("pCosAmp495", &pCosA495, "pCosA495/D");
    tree->Branch("pCosAmpErr495", &pCosAErr495, "pCosAErr495/D");
    tree->Branch("pCosAmp505", &pCosA505, "pCosA505/D");
    tree->Branch("pCosAmpErr505", &pCosAErr505, "pCosAErr505/D");
    tree->Branch("pCosAmp515", &pCosA515, "pCosA515/D");
    tree->Branch("pCosAmpErr515", &pCosAErr515, "pCosAErr515/D");
    tree->Branch("pCosAmp525", &pCosA525, "pCosA525/D");
    tree->Branch("pCosAmpErr525", &pCosAErr525, "pCosAErr525/D");
    tree->Branch("pCosAmp535", &pCosA535, "pCosA535/D");
    tree->Branch("pCosAmpErr535", &pCosAErr535, "pCosAErr535/D");
    tree->Branch("pCosAmp545", &pCosA545, "pCosA545/D");
    tree->Branch("pCosAmpErr545", &pCosAErr545, "pCosAErr545/D");
    tree->Branch("pCosAmp555", &pCosA555, "pCosA555/D");
    tree->Branch("pCosAmpErr555", &pCosAErr555, "pCosAErr555/D");
    tree->Branch("pCosAmp565", &pCosA565, "pCosA565/D");
    tree->Branch("pCosAmpErr565", &pCosAErr565, "pCosAErr565/D");
    tree->Branch("pCosAmp575", &pCosA575, "pCosA575/D");
    tree->Branch("pCosAmpErr575", &pCosAErr575, "pCosAErr575/D");
    tree->Branch("pCosAmp585", &pCosA585, "pCosA585/D");
    tree->Branch("pCosAmpErr585", &pCosAErr585, "pCosAErr585/D");
    tree->Branch("pCosAmp595", &pCosA595, "pCosA595/D");
    tree->Branch("pCosAmpErr595", &pCosAErr595, "pCosAErr595/D");
    tree->Branch("pCosAmp605", &pCosA605, "pCosA605/D");
    tree->Branch("pCosAmpErr605", &pCosAErr605, "pCosAErr605/D");
    tree->Branch("pCosAmp615", &pCosA615, "pCosA615/D");
    tree->Branch("pCosAmpErr615", &pCosAErr615, "pCosAErr615/D");

    // Fill branches (and hence tree) with corresponding parameters from above
    for (Int_t m = 0; m < 20; m++){
        pCosA415 = pCosAmp[0][m];
        pCosAErr415 = pCosAmpErr[0][m];
        pCosA425 = pCosAmp[1][m];
        pCosAErr425 = pCosAmpErr[1][m];
        pCosA435 = pCosAmp[2][m];
        pCosAErr435 = pCosAmpErr[2][m];
        pCosA445 = pCosAmp[3][m];
        pCosAErr445 = pCosAmpErr[3][m];
        pCosA455 = pCosAmp[4][m];
        pCosAErr455 = pCosAmpErr[4][m];
        pCosA465 = pCosAmp[5][m];
        pCosAErr465 = pCosAmpErr[5][m];
        pCosA475 = pCosAmp[6][m];
        pCosAErr475 = pCosAmpErr[6][m];
        pCosA485 = pCosAmp[7][m];
        pCosAErr485 = pCosAmpErr[7][m];
        pCosA495 = pCosAmp[8][m];
        pCosAErr495 = pCosAmpErr[8][m];
        pCosA505 = pCosAmp[9][m];
        pCosAErr505 = pCosAmpErr[9][m];
        pCosA515 = pCosAmp[10][m];
        pCosAErr515 = pCosAmpErr[10][m];
        pCosA525 = pCosAmp[11][m];
        pCosAErr525 = pCosAmpErr[11][m];
        pCosA535 = pCosAmp[12][m];
        pCosAErr535 = pCosAmpErr[12][m];
        pCosA545 = pCosAmp[13][m];
        pCosAErr545= pCosAmpErr[13][m];
        //pCosA555 = pCosAmp[14][m];
        //pCosAErr555 = pCosAmpErr[14][m];
        pCosA555 = 0;
        pCosAErr555 = 0;
        pCosA565 = pCosAmp[15][m];
        pCosAErr565 = pCosAmpErr[15][m];
        pCosA575 = pCosAmp[16][m];
        pCosAErr575 = pCosAmpErr[16][m];
        pCosA585 = pCosAmp[17][m];
        pCosAErr585 = pCosAmpErr[17][m];
        pCosA595 = pCosAmp[18][m];
        pCosAErr595 = pCosAmpErr[18][m];
        pCosA605 = pCosAmp[19][m];
        pCosAErr605 = pCosAmpErr[19][m];
        pCosA615 = pCosAmp[20][m];
        pCosAErr615 = pCosAmpErr[20][m];

        tree->Fill();
    }

    for(Int_t i = 0; i < 21; i++){ // Energy
        for(Int_t j = 0; j < 20; j++){ // Theta
            AsymmHists[i][j]->Write();
        }
    }

    f1.Write();

}
