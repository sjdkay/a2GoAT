#include "./includes_SigmaAsymm_NoScatt.h"

void SigmaAsymm_NoScatt_V2(){

    TF1 *CosFunc = new TF1("CosFit", "[0]*cos((2*x*TMath::DegToRad())+acos(0))");
    CosFunc->SetParNames("Amplitude");

    double pCosAmp[20][21]; // Format of array is Theta bin (x) by Egamma bin (y)
    double pCosAmpErr[20][21];
    double pCosA415;
    double pCosAErr415;
    double pCosA425;
    double pCosAErr425;
    double pCosA435;
    double pCosAErr435;
    double pCosA445;
    double pCosAErr445;
    double pCosA455;
    double pCosAErr455;
    double pCosA465;
    double pCosAErr465;
    double pCosA475;
    double pCosAErr475;
    double pCosA485;
    double pCosAErr485;
    double pCosA495;
    double pCosAErr495;
    double pCosA505;
    double pCosAErr505;
    double pCosA515;
    double pCosAErr515;
    double pCosA525;
    double pCosAErr525;
    double pCosA535;
    double pCosAErr535;
    double pCosA545;
    double pCosAErr545;
    double pCosA555;
    double pCosAErr555;
    double pCosA565;
    double pCosAErr565;
    double pCosA575;
    double pCosAErr575;
    double pCosA585;
    double pCosAErr585;
    double pCosA595;
    double pCosAErr595;
    double pCosA605;
    double pCosAErr605;
    double pCosA615;
    double pCosAErr615;

    TList *Phip415AsymmList = new TList;
    TList *Phip425AsymmList = new TList;
    TList *Phip435AsymmList = new TList;
    TList *Phip445AsymmList = new TList;
    TList *Phip455AsymmList = new TList;
    TList *Phip465AsymmList = new TList;
    TList *Phip475AsymmList = new TList;
    TList *Phip485AsymmList = new TList;
    TList *Phip495AsymmList = new TList;
    TList *Phip505AsymmList = new TList;
    TList *Phip515AsymmList = new TList;
    TList *Phip525AsymmList = new TList;
    TList *Phip535AsymmList = new TList;
    TList *Phip545AsymmList = new TList;
    TList *Phip555AsymmList = new TList;
    TList *Phip565AsymmList = new TList;
    TList *Phip575AsymmList = new TList;
    TList *Phip585AsymmList = new TList;
    TList *Phip595AsymmList = new TList;
    TList *Phip605AsymmList = new TList;
    TList *Phip615AsymmList = new TList;

    TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/ParaPerp_NS14_Combined.root"); // Open the latest PTotal combined file to load histograms from
    NPara = Eg_Para->GetEntries();
    NPerp = Eg_Perp->GetEntries();
    ScaleFactor = NPara/NPerp;
    ScaleFactorErr = sqrt( (NPara/((TMath::Power(NPerp,2)))) + (((TMath::Power(NPara,2)))/(TMath::Power(NPerp,3))) ); // Error Propagation of above formula, see notebook 5

    for(Int_t i = 0; i < 20; i++){
        Sig415Asymm[i] = (TH1D*) ((TH1D*) Phip415ParaList->At(i))->GetAsymmetry(((TH1D*) Phip415PerpList->At(i)), ScaleFactor, ScaleFactorErr);
        Sig415Asymm[i]-> Fit("CosFit", "Q");
        Phip415AsymmList->Add(Sig415Asymm[i]);
        pCosAmp[i][0] = CosFit->GetParameter(0);
        pCosAmpErr[i][0] = CosFit->GetParError(0);
    }

    for(Int_t i = 0; i < 20; i++){
        Sig425Asymm[i] = (TH1D*) ((TH1D*) Phip425ParaList->At(i))->GetAsymmetry(((TH1D*) Phip425PerpList->At(i)), ScaleFactor, ScaleFactorErr);
        Sig425Asymm[i]-> Fit("CosFit", "Q");
        Phip425AsymmList->Add(Sig425Asymm[i]);
        pCosAmp[i][1] = CosFit->GetParameter(0);
        pCosAmpErr[i][1] = CosFit->GetParError(0);
    }

    for(Int_t i = 0; i < 20; i++){
        Sig435Asymm[i] = (TH1D*) ((TH1D*) Phip435ParaList->At(i))->GetAsymmetry(((TH1D*) Phip435PerpList->At(i)), ScaleFactor, ScaleFactorErr);
        Sig435Asymm[i]-> Fit("CosFit", "Q");
        Phip435AsymmList->Add(Sig435Asymm[i]);
        pCosAmp[i][2] = CosFit->GetParameter(0);
        pCosAmpErr[i][2] = CosFit->GetParError(0);
    }

    for(Int_t i = 0; i < 20; i++){
        Sig445Asymm[i] = (TH1D*) ((TH1D*) Phip445ParaList->At(i))->GetAsymmetry(((TH1D*) Phip445PerpList->At(i)), ScaleFactor, ScaleFactorErr);
        Sig445Asymm[i]-> Fit("CosFit", "Q");
        Phip445AsymmList->Add(Sig445Asymm[i]);
        pCosAmp[i][3] = CosFit->GetParameter(0);
        pCosAmpErr[i][3] = CosFit->GetParError(0);
    }

    for(Int_t i = 0; i < 20; i++){
        Sig455Asymm[i] = (TH1D*) ((TH1D*) Phip455ParaList->At(i))->GetAsymmetry(((TH1D*) Phip455PerpList->At(i)), ScaleFactor, ScaleFactorErr);
        Sig455Asymm[i]-> Fit("CosFit", "Q");
        Phip455AsymmList->Add(Sig455Asymm[i]);
        pCosAmp[i][4] = CosFit->GetParameter(0);
        pCosAmpErr[i][4] = CosFit->GetParError(0);
    }

    for(Int_t i = 0; i < 20; i++){
        Sig465Asymm[i] = (TH1D*) ((TH1D*) Phip465ParaList->At(i))->GetAsymmetry(((TH1D*) Phip465PerpList->At(i)), ScaleFactor, ScaleFactorErr);
        Sig465Asymm[i]-> Fit("CosFit", "Q");
        Phip465AsymmList->Add(Sig465Asymm[i]);
        pCosAmp[i][5] = CosFit->GetParameter(0);
        pCosAmpErr[i][5] = CosFit->GetParError(0);
    }

    for(Int_t i = 0; i < 20; i++){
        Sig475Asymm[i] = (TH1D*) ((TH1D*) Phip475ParaList->At(i))->GetAsymmetry(((TH1D*) Phip475PerpList->At(i)), ScaleFactor, ScaleFactorErr);
        Sig475Asymm[i]-> Fit("CosFit", "Q");
        Phip475AsymmList->Add(Sig475Asymm[i]);
        pCosAmp[i][6] = CosFit->GetParameter(0);
        pCosAmpErr[i][6] = CosFit->GetParError(0);
    }

    for(Int_t i = 0; i < 20; i++){
        Sig485Asymm[i] = (TH1D*) ((TH1D*) Phip485ParaList->At(i))->GetAsymmetry(((TH1D*) Phip485PerpList->At(i)), ScaleFactor, ScaleFactorErr);
        Sig485Asymm[i]-> Fit("CosFit", "Q");
        Phip485AsymmList->Add(Sig485Asymm[i]);
        pCosAmp[i][7] = CosFit->GetParameter(0);
        pCosAmpErr[i][7] = CosFit->GetParError(0);
    }

    for(Int_t i = 0; i < 20; i++){
        Sig495Asymm[i] = (TH1D*) ((TH1D*) Phip495ParaList->At(i))->GetAsymmetry(((TH1D*) Phip495PerpList->At(i)), ScaleFactor, ScaleFactorErr);
        Sig495Asymm[i]-> Fit("CosFit", "Q");
        Phip495AsymmList->Add(Sig495Asymm[i]);
        pCosAmp[i][8] = CosFit->GetParameter(0);
        pCosAmpErr[i][8] = CosFit->GetParError(0);
    }

    for(Int_t i = 0; i < 20; i++){
        Sig505Asymm[i] = (TH1D*) ((TH1D*) Phip505ParaList->At(i))->GetAsymmetry(((TH1D*) Phip505PerpList->At(i)), ScaleFactor, ScaleFactorErr);
        Sig505Asymm[i]-> Fit("CosFit", "Q");
        Phip505AsymmList->Add(Sig505Asymm[i]);
        pCosAmp[i][9] = CosFit->GetParameter(0);
        pCosAmpErr[i][9] = CosFit->GetParError(0);
    }

    for(Int_t i = 0; i < 20; i++){
        Sig515Asymm[i] = (TH1D*) ((TH1D*) Phip515ParaList->At(i))->GetAsymmetry(((TH1D*) Phip515PerpList->At(i)), ScaleFactor, ScaleFactorErr);
        Sig515Asymm[i]-> Fit("CosFit", "Q");
        Phip515AsymmList->Add(Sig515Asymm[i]);
        pCosAmp[i][10] = CosFit->GetParameter(0);
        pCosAmpErr[i][10] = CosFit->GetParError(0);
    }

    for(Int_t i = 0; i < 20; i++){
        Sig525Asymm[i] = (TH1D*) ((TH1D*) Phip525ParaList->At(i))->GetAsymmetry(((TH1D*) Phip525PerpList->At(i)), ScaleFactor, ScaleFactorErr);
        Sig525Asymm[i]-> Fit("CosFit", "Q");
        Phip525AsymmList->Add(Sig525Asymm[i]);
        pCosAmp[i][11] = CosFit->GetParameter(0);
        pCosAmpErr[i][11] = CosFit->GetParError(0);
    }

    for(Int_t i = 0; i < 20; i++){
        Sig535Asymm[i] = (TH1D*) ((TH1D*) Phip535ParaList->At(i))->GetAsymmetry(((TH1D*) Phip535PerpList->At(i)), ScaleFactor, ScaleFactorErr);
        Sig535Asymm[i]-> Fit("CosFit", "Q");
        Phip535AsymmList->Add(Sig535Asymm[i]);
        pCosAmp[i][12] = CosFit->GetParameter(0);
        pCosAmpErr[i][12] = CosFit->GetParError(0);
    }

    for(Int_t i = 0; i < 20; i++){
        Sig545Asymm[i] = (TH1D*) ((TH1D*) Phip545ParaList->At(i))->GetAsymmetry(((TH1D*) Phip545PerpList->At(i)), ScaleFactor, ScaleFactorErr);
        Sig545Asymm[i]-> Fit("CosFit", "Q");
        Phip545AsymmList->Add(Sig545Asymm[i]);
        pCosAmp[i][13] = CosFit->GetParameter(0);
        pCosAmpErr[i][13] = CosFit->GetParError(0);
    }

    for(Int_t i = 0; i < 20; i++){
        Sig555Asymm[i] = (TH1D*) ((TH1D*) Phip555ParaList->At(i))->GetAsymmetry(((TH1D*) Phip555PerpList->At(i)), ScaleFactor, ScaleFactorErr);
        Sig555Asymm[i]-> Fit("CosFit", "Q");
        Phip555AsymmList->Add(Sig555Asymm[i]);
        pCosAmp[i][14] = CosFit->GetParameter(0);
        pCosAmpErr[i][14] = CosFit->GetParError(0);
    }

    for(Int_t i = 0; i < 20; i++){
        Sig565Asymm[i] = (TH1D*) ((TH1D*) Phip565ParaList->At(i))->GetAsymmetry(((TH1D*) Phip565PerpList->At(i)), ScaleFactor, ScaleFactorErr);
        Sig565Asymm[i]-> Fit("CosFit", "Q");
        Phip565AsymmList->Add(Sig565Asymm[i]);
        pCosAmp[i][15] = CosFit->GetParameter(0);
        pCosAmpErr[i][15] = CosFit->GetParError(0);
    }

    for(Int_t i = 0; i < 20; i++){
        Sig575Asymm[i] = (TH1D*) ((TH1D*) Phip575ParaList->At(i))->GetAsymmetry(((TH1D*) Phip575PerpList->At(i)), ScaleFactor, ScaleFactorErr);
        Sig575Asymm[i]-> Fit("CosFit", "Q");
        Phip575AsymmList->Add(Sig575Asymm[i]);
        pCosAmp[i][16] = CosFit->GetParameter(0);
        pCosAmpErr[i][16] = CosFit->GetParError(0);
    }

    for(Int_t i = 0; i < 20; i++){
        Sig585Asymm[i] = (TH1D*) ((TH1D*) Phip585ParaList->At(i))->GetAsymmetry(((TH1D*) Phip585PerpList->At(i)), ScaleFactor, ScaleFactorErr);
        Sig585Asymm[i]-> Fit("CosFit", "Q");
        Phip585AsymmList->Add(Sig585Asymm[i]);
        pCosAmp[i][17] = CosFit->GetParameter(0);
        pCosAmpErr[i][17] = CosFit->GetParError(0);
    }

    for(Int_t i = 0; i < 20; i++){
        Sig595Asymm[i] = (TH1D*) ((TH1D*) Phip595ParaList->At(i))->GetAsymmetry(((TH1D*) Phip595PerpList->At(i)), ScaleFactor, ScaleFactorErr);
        Sig595Asymm[i]-> Fit("CosFit", "Q");
        Phip595AsymmList->Add(Sig595Asymm[i]);
        pCosAmp[i][18] = CosFit->GetParameter(0);
        pCosAmpErr[i][18] = CosFit->GetParError(0);
    }

    for(Int_t i = 0; i < 20; i++){
        Sig605Asymm[i] = (TH1D*) ((TH1D*) Phip605ParaList->At(i))->GetAsymmetry(((TH1D*) Phip605PerpList->At(i)), ScaleFactor, ScaleFactorErr);
        Sig605Asymm[i]-> Fit("CosFit", "Q");
        Phip605AsymmList->Add(Sig605Asymm[i]);
        pCosAmp[i][19] = CosFit->GetParameter(0);
        pCosAmpErr[i][19] = CosFit->GetParError(0);
    }

    for(Int_t i = 0; i < 20; i++){
        Sig615Asymm[i] = (TH1D*) ((TH1D*) Phip615ParaList->At(i))->GetAsymmetry(((TH1D*) Phip615PerpList->At(i)), ScaleFactor, ScaleFactorErr);
        Sig615Asymm[i]-> Fit("CosFit", "Q");
        Phip615AsymmList->Add(Sig615Asymm[i]);
        pCosAmp[i][20] = CosFit->GetParameter(0);
        pCosAmpErr[i][20] = CosFit->GetParError(0);
    }

    TFile f1("ParaPerpAsymm_NS14.root", "RECREATE");

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
        pCosA415 = pCosAmp[m][0];
        pCosAErr415 = pCosAmpErr[m][0];
        pCosA425 = pCosAmp[m][1];
        pCosAErr425 = pCosAmpErr[m][1];
        pCosA435 = pCosAmp[m][2];
        pCosAErr435 = pCosAmpErr[m][2];
        pCosA445 = pCosAmp[m][3];
        pCosAErr445 = pCosAmpErr[m][3];
        pCosA455 = pCosAmp[m][4];
        pCosAErr455 = pCosAmpErr[m][4];
        pCosA465 = pCosAmp[m][5];
        pCosAErr465 = pCosAmpErr[m][5];
        pCosA475 = pCosAmp[m][6];
        pCosAErr475 = pCosAmpErr[m][6];
        pCosA485 = pCosAmp[m][7];
        pCosAErr485 = pCosAmpErr[m][7];
        pCosA495 = pCosAmp[m][8];
        pCosAErr495 = pCosAmpErr[m][8];
        pCosA505 = pCosAmp[m][9];
        pCosAErr505 = pCosAmpErr[m][9];
        pCosA515 = pCosAmp[m][10];
        pCosAErr515 = pCosAmpErr[m][10];
        pCosA525 = pCosAmp[m][11];
        pCosAErr525 = pCosAmpErr[m][11];
        pCosA535 = pCosAmp[m][12];
        pCosAErr535 = pCosAmpErr[m][12];
        pCosA545 = pCosAmp[m][13];
        pCosAErr545= pCosAmpErr[m][13];
        //pCosA555 = pCosAmp[m][14];
        //pCosAErr555 = pCosAmpErr[m][14];
        pCosA555 = 0;
        pCosAErr555 = 0;
        pCosA565 = pCosAmp[m][15];
        pCosAErr565 = pCosAmpErr[m][15];
        pCosA575 = pCosAmp[m][16];
        pCosAErr575 = pCosAmpErr[m][16];
        pCosA585 = pCosAmp[m][17];
        pCosAErr585 = pCosAmpErr[m][17];
        pCosA595 = pCosAmp[m][18];
        pCosAErr595 = pCosAmpErr[m][18];
        pCosA605 = pCosAmp[m][19];
        pCosAErr605 = pCosAmpErr[m][19];
        pCosA615 = pCosAmp[m][20];
        pCosAErr615 = pCosAmpErr[m][20];

        tree->Fill();
    }

    Phip415AsymmList -> Write("Phip415AsymmList", TObject::kSingleKey);
    Phip425AsymmList -> Write("Phip425AsymmList", TObject::kSingleKey);
    Phip435AsymmList -> Write("Phip435AsymmList", TObject::kSingleKey);
    Phip445AsymmList -> Write("Phip445AsymmList", TObject::kSingleKey);
    Phip455AsymmList -> Write("Phip455AsymmList", TObject::kSingleKey);
    Phip465AsymmList -> Write("Phip465AsymmList", TObject::kSingleKey);
    Phip475AsymmList -> Write("Phip475AsymmList", TObject::kSingleKey);
    Phip485AsymmList -> Write("Phip485AsymmList", TObject::kSingleKey);
    Phip495AsymmList -> Write("Phip495AsymmList", TObject::kSingleKey);
    Phip505AsymmList -> Write("Phip505AsymmList", TObject::kSingleKey);
    Phip515AsymmList -> Write("Phip515AsymmList", TObject::kSingleKey);
    Phip525AsymmList -> Write("Phip525AsymmList", TObject::kSingleKey);
    Phip535AsymmList -> Write("Phip535AsymmList", TObject::kSingleKey);
    Phip545AsymmList -> Write("Phip545AsymmList", TObject::kSingleKey);
    Phip555AsymmList -> Write("Phip555AsymmList", TObject::kSingleKey);
    Phip565AsymmList -> Write("Phip565AsymmList", TObject::kSingleKey);
    Phip575AsymmList -> Write("Phip575AsymmList", TObject::kSingleKey);
    Phip585AsymmList -> Write("Phip585AsymmList", TObject::kSingleKey);
    Phip595AsymmList -> Write("Phip595AsymmList", TObject::kSingleKey);
    Phip605AsymmList -> Write("Phip605AsymmList", TObject::kSingleKey);
    Phip615AsymmList -> Write("Phip615AsymmList", TObject::kSingleKey);

    f1.Write();

}
