#include "./includes_Sigma_Systematics.h"

void Sigma_Systematics(){

    TFile *f1= TFile::Open("/scratch/Mainz_Software/a2GoAT/Results/ParaPerpAsymm_NS18.root"); // 2 Sigma missing mass cut fil
    TTree *t1 = (TTree*)f1->Get("Parameter_Values");

    TGraphErrors* SigmaSystPlots[21];
    char name[21];
    char title[60];

    // Set branch addresses to get values from
    t1->SetBranchAddress("pCosAmp415", &pCosAmp415);
    t1->SetBranchAddress("pCosAmpErr415", &pCosAmpErr415);
    t1->SetBranchAddress("pCosAmp425", &pCosAmp425);
    t1->SetBranchAddress("pCosAmpErr425", &pCosAmpErr425);
    t1->SetBranchAddress("pCosAmp435", &pCosAmp435);
    t1->SetBranchAddress("pCosAmpErr435", &pCosAmpErr435);
    t1->SetBranchAddress("pCosAmp445", &pCosAmp445);
    t1->SetBranchAddress("pCosAmpErr445", &pCosAmpErr445);
    t1->SetBranchAddress("pCosAmp455", &pCosAmp455);
    t1->SetBranchAddress("pCosAmpErr455", &pCosAmpErr455);
    t1->SetBranchAddress("pCosAmp465", &pCosAmp465);
    t1->SetBranchAddress("pCosAmpErr465", &pCosAmpErr465);
    t1->SetBranchAddress("pCosAmp475", &pCosAmp475);
    t1->SetBranchAddress("pCosAmpErr475", &pCosAmpErr475);
    t1->SetBranchAddress("pCosAmp485", &pCosAmp485);
    t1->SetBranchAddress("pCosAmpErr485", &pCosAmpErr485);
    t1->SetBranchAddress("pCosAmp495", &pCosAmp495);
    t1->SetBranchAddress("pCosAmpErr495", &pCosAmpErr495);
    t1->SetBranchAddress("pCosAmp505", &pCosAmp505);
    t1->SetBranchAddress("pCosAmpErr505", &pCosAmpErr505);
    t1->SetBranchAddress("pCosAmp515", &pCosAmp515);
    t1->SetBranchAddress("pCosAmpErr515", &pCosAmpErr515);
    t1->SetBranchAddress("pCosAmp525", &pCosAmp525);
    t1->SetBranchAddress("pCosAmpErr525", &pCosAmpErr525);
    t1->SetBranchAddress("pCosAmp535", &pCosAmp535);
    t1->SetBranchAddress("pCosAmpErr535", &pCosAmpErr535);
    t1->SetBranchAddress("pCosAmp545", &pCosAmp545);
    t1->SetBranchAddress("pCosAmpErr545", &pCosAmpErr545);
    t1->SetBranchAddress("pCosAmp555", &pCosAmp555);
    t1->SetBranchAddress("pCosAmpErr555", &pCosAmpErr555);
    t1->SetBranchAddress("pCosAmp565", &pCosAmp565);
    t1->SetBranchAddress("pCosAmpErr565", &pCosAmpErr565);
    t1->SetBranchAddress("pCosAmp575", &pCosAmp575);
    t1->SetBranchAddress("pCosAmpErr575", &pCosAmpErr575);
    t1->SetBranchAddress("pCosAmp585", &pCosAmp585);
    t1->SetBranchAddress("pCosAmpErr585", &pCosAmpErr585);
    t1->SetBranchAddress("pCosAmp595", &pCosAmp595);
    t1->SetBranchAddress("pCosAmpErr595", &pCosAmpErr595);
    t1->SetBranchAddress("pCosAmp605", &pCosAmp605);
    t1->SetBranchAddress("pCosAmpErr605", &pCosAmpErr605);
    t1->SetBranchAddress("pCosAmp615", &pCosAmp615);
    t1->SetBranchAddress("pCosAmpErr615", &pCosAmpErr615);

    for (Int_t k = 0; k < 20; k++){
        Parameter_Values->GetEntry(k);
        pValues415[k] = pCosAmp415;
        pErrValues415[k] = pCosAmpErr415;
        pValues425[k] = pCosAmp425;
        pErrValues425[k] = pCosAmpErr425;
        pValues435[k] = pCosAmp435;
        pErrValues435[k] = pCosAmpErr435;
        pValues445[k] = pCosAmp445;
        pErrValues445[k] = pCosAmpErr445;
        pValues455[k] = pCosAmp455;
        pErrValues455[k] = pCosAmpErr455;
        pValues465[k] = pCosAmp465;
        pErrValues465[k] = pCosAmpErr465;
        pValues475[k] = pCosAmp475;
        pErrValues475[k] = pCosAmpErr475;
        pValues485[k] = pCosAmp485;
        pErrValues485[k] = pCosAmpErr485;
        pValues495[k] = pCosAmp495;
        pErrValues495[k] = pCosAmpErr495;
        pValues505[k] = pCosAmp505;
        pErrValues505[k] = pCosAmpErr505;
        pValues515[k] = pCosAmp515;
        pErrValues515[k] = pCosAmpErr515;
        pValues525[k] = pCosAmp525;
        pErrValues525[k] = pCosAmpErr525;
        pValues535[k] = pCosAmp535;
        pErrValues535[k] = pCosAmpErr535;
        pValues545[k] = pCosAmp545;
        pErrValues545[k] = pCosAmpErr545;
        pValues555[k] = pCosAmp555;
        pErrValues555[k] = pCosAmpErr555;
        pValues565[k] = pCosAmp565;
        pErrValues565[k] = pCosAmpErr565;
        pValues575[k] = pCosAmp575;
        pErrValues575[k] = pCosAmpErr575;
        pValues585[k] = pCosAmp585;
        pErrValues585[k] = pCosAmpErr585;
        pValues595[k] = pCosAmp595;
        pErrValues595[k] = pCosAmpErr595;
        pValues605[k] = pCosAmp605;
        pErrValues605[k] = pCosAmpErr605;
        pValues615[k] = pCosAmp615;
        pErrValues615[k] = pCosAmpErr615;

    }

    TFile *f3 = TFile::Open("/scratch/Mainz_Software/a2GoAT/LinPol_Aug16.root"); // Open linear polarisation plot

    for (Int_t i = 0; i < 18; i++){

        SigmaValues[0][0][i] = pValues415[i+1]/(Graph->Eval(425,0));
        SigmaErrValues[0][0][i] = pErrValues415[i+1]/(Graph->Eval(425,0));
        SigmaValues[0][1][i] = pValues425[i+1]/(Graph->Eval(425,0));
        SigmaErrValues[0][1][i] = pErrValues425[i+1]/(Graph->Eval(425,0));
        SigmaValues[0][2][i] = pValues435[i+1]/(Graph->Eval(435,0));
        SigmaErrValues[0][2][i] = pErrValues435[i+1]/(Graph->Eval(435,0));
        SigmaValues[0][3][i] = pValues445[i+1]/(Graph->Eval(445,0));
        SigmaErrValues[0][3][i] = pErrValues445[i+1]/(Graph->Eval(445,0));
        SigmaValues[0][4][i] = pValues455[i+1]/(Graph->Eval(455,0));
        SigmaErrValues[0][4][i] = pErrValues455[i+1]/(Graph->Eval(455,0));
        SigmaValues[0][5][i] = pValues465[i+1]/(Graph->Eval(465,0));
        SigmaErrValues[0][5][i] = pErrValues465[i+1]/(Graph->Eval(465,0));
        SigmaValues[0][6][i] = pValues475[i+1]/(Graph->Eval(475,0));
        SigmaErrValues[0][6][i] = pErrValues475[i+1]/(Graph->Eval(475,0));
        SigmaValues[0][7][i] = pValues485[i+1]/(Graph->Eval(485,0));
        SigmaErrValues[0][7][i] = pErrValues485[i+1]/(Graph->Eval(485,0));
        SigmaValues[0][8][i] = pValues495[i+1]/(Graph->Eval(495,0));
        SigmaErrValues[0][8][i] = pErrValues495[i+1]/(Graph->Eval(495,0));
        SigmaValues[0][9][i] = pValues505[i+1]/(Graph->Eval(505,0));
        SigmaErrValues[0][9][i] = pErrValues505[i+1]/(Graph->Eval(505,0));
        SigmaValues[0][10][i] = pValues515[i+1]/(Graph->Eval(515,0));
        SigmaErrValues[0][10][i] = pErrValues515[i+1]/(Graph->Eval(515,0));
        SigmaValues[0][11][i] = pValues525[i+1]/(Graph->Eval(525,0));
        SigmaErrValues[0][11][i] = pErrValues525[i+1]/(Graph->Eval(525,0));
        SigmaValues[0][12][i] = pValues535[i+1]/(Graph->Eval(535,0));
        SigmaErrValues[0][12][i] = pErrValues535[i+1]/(Graph->Eval(535,0));
        SigmaValues[0][13][i] = pValues545[i+1]/(Graph->Eval(545,0));
        SigmaErrValues[0][13][i] = pErrValues545[i+1]/(Graph->Eval(545,0));
        SigmaValues[0][14][i] = pValues555[i+1]/(Graph->Eval(555,0));
        SigmaErrValues[0][14][i] = pErrValues555[i+1]/(Graph->Eval(555,0));
        SigmaValues[0][15][i] = pValues565[i+1]/(Graph->Eval(565,0));
        SigmaErrValues[0][15][i] = pErrValues565[i+1]/(Graph->Eval(565,0));
        SigmaValues[0][16][i] = pValues575[i+1]/(Graph->Eval(575,0));
        SigmaErrValues[0][16][i] = pErrValues575[i+1]/(Graph->Eval(575,0));
        SigmaValues[0][17][i] = pValues585[i+1]/(Graph->Eval(585,0));
        SigmaErrValues[0][17][i] = pErrValues585[i+1]/(Graph->Eval(585,0));
        SigmaValues[0][18][i] = pValues595[i+1]/(Graph->Eval(595,0));
        SigmaErrValues[0][18][i] = pErrValues595[i+1]/(Graph->Eval(595,0));
        SigmaValues[0][19][i] = pValues605[i+1]/(Graph->Eval(605,0));
        SigmaErrValues[0][19][i] = pErrValues605[i+1]/(Graph->Eval(605,0));
        SigmaValues[0][20][i] = pValues615[i+1]/(Graph->Eval(615,0));
        SigmaErrValues[0][20][i] = pErrValues615[i+1]/(Graph->Eval(615,0));

    }

    TFile *f2= TFile::Open("/scratch/Mainz_Software/a2GoAT/Results/ParaPerpAsymm_NS17.root"); //1 Sigma missing mass cut file
    TTree *t2 = (TTree*)f2->Get("Parameter_Values");

    // Set branch addresses to get values from
    t2->SetBranchAddress("pCosAmp415", &pCosAmp415);
    t2->SetBranchAddress("pCosAmpErr415", &pCosAmpErr415);
    t2->SetBranchAddress("pCosAmp425", &pCosAmp425);
    t2->SetBranchAddress("pCosAmpErr425", &pCosAmpErr425);
    t2->SetBranchAddress("pCosAmp435", &pCosAmp435);
    t2->SetBranchAddress("pCosAmpErr435", &pCosAmpErr435);
    t2->SetBranchAddress("pCosAmp445", &pCosAmp445);
    t2->SetBranchAddress("pCosAmpErr445", &pCosAmpErr445);
    t2->SetBranchAddress("pCosAmp455", &pCosAmp455);
    t2->SetBranchAddress("pCosAmpErr455", &pCosAmpErr455);
    t2->SetBranchAddress("pCosAmp465", &pCosAmp465);
    t2->SetBranchAddress("pCosAmpErr465", &pCosAmpErr465);
    t2->SetBranchAddress("pCosAmp475", &pCosAmp475);
    t2->SetBranchAddress("pCosAmpErr475", &pCosAmpErr475);
    t2->SetBranchAddress("pCosAmp485", &pCosAmp485);
    t2->SetBranchAddress("pCosAmpErr485", &pCosAmpErr485);
    t2->SetBranchAddress("pCosAmp495", &pCosAmp495);
    t2->SetBranchAddress("pCosAmpErr495", &pCosAmpErr495);
    t2->SetBranchAddress("pCosAmp505", &pCosAmp505);
    t2->SetBranchAddress("pCosAmpErr505", &pCosAmpErr505);
    t2->SetBranchAddress("pCosAmp515", &pCosAmp515);
    t2->SetBranchAddress("pCosAmpErr515", &pCosAmpErr515);
    t2->SetBranchAddress("pCosAmp525", &pCosAmp525);
    t2->SetBranchAddress("pCosAmpErr525", &pCosAmpErr525);
    t2->SetBranchAddress("pCosAmp535", &pCosAmp535);
    t2->SetBranchAddress("pCosAmpErr535", &pCosAmpErr535);
    t2->SetBranchAddress("pCosAmp545", &pCosAmp545);
    t2->SetBranchAddress("pCosAmpErr545", &pCosAmpErr545);
    t2->SetBranchAddress("pCosAmp555", &pCosAmp555);
    t2->SetBranchAddress("pCosAmpErr555", &pCosAmpErr555);
    t2->SetBranchAddress("pCosAmp565", &pCosAmp565);
    t2->SetBranchAddress("pCosAmpErr565", &pCosAmpErr565);
    t2->SetBranchAddress("pCosAmp575", &pCosAmp575);
    t2->SetBranchAddress("pCosAmpErr575", &pCosAmpErr575);
    t2->SetBranchAddress("pCosAmp585", &pCosAmp585);
    t2->SetBranchAddress("pCosAmpErr585", &pCosAmpErr585);
    t2->SetBranchAddress("pCosAmp595", &pCosAmp595);
    t2->SetBranchAddress("pCosAmpErr595", &pCosAmpErr595);
    t2->SetBranchAddress("pCosAmp605", &pCosAmp605);
    t2->SetBranchAddress("pCosAmpErr605", &pCosAmpErr605);
    t2->SetBranchAddress("pCosAmp615", &pCosAmp615);
    t2->SetBranchAddress("pCosAmpErr615", &pCosAmpErr615);

    for (Int_t k = 0; k < 20; k++){
        Parameter_Values->GetEntry(k);
        pValues415[k] = pCosAmp415;
        pErrValues415[k] = pCosAmpErr415;
        pValues425[k] = pCosAmp425;
        pErrValues425[k] = pCosAmpErr425;
        pValues435[k] = pCosAmp435;
        pErrValues435[k] = pCosAmpErr435;
        pValues445[k] = pCosAmp445;
        pErrValues445[k] = pCosAmpErr445;
        pValues455[k] = pCosAmp455;
        pErrValues455[k] = pCosAmpErr455;
        pValues465[k] = pCosAmp465;
        pErrValues465[k] = pCosAmpErr465;
        pValues475[k] = pCosAmp475;
        pErrValues475[k] = pCosAmpErr475;
        pValues485[k] = pCosAmp485;
        pErrValues485[k] = pCosAmpErr485;
        pValues495[k] = pCosAmp495;
        pErrValues495[k] = pCosAmpErr495;
        pValues505[k] = pCosAmp505;
        pErrValues505[k] = pCosAmpErr505;
        pValues515[k] = pCosAmp515;
        pErrValues515[k] = pCosAmpErr515;
        pValues525[k] = pCosAmp525;
        pErrValues525[k] = pCosAmpErr525;
        pValues535[k] = pCosAmp535;
        pErrValues535[k] = pCosAmpErr535;
        pValues545[k] = pCosAmp545;
        pErrValues545[k] = pCosAmpErr545;
        pValues555[k] = pCosAmp555;
        pErrValues555[k] = pCosAmpErr555;
        pValues565[k] = pCosAmp565;
        pErrValues565[k] = pCosAmpErr565;
        pValues575[k] = pCosAmp575;
        pErrValues575[k] = pCosAmpErr575;
        pValues585[k] = pCosAmp585;
        pErrValues585[k] = pCosAmpErr585;
        pValues595[k] = pCosAmp595;
        pErrValues595[k] = pCosAmpErr595;
        pValues605[k] = pCosAmp605;
        pErrValues605[k] = pCosAmpErr605;
        pValues615[k] = pCosAmp615;
        pErrValues615[k] = pCosAmpErr615;

    }

    TFile *f4 = TFile::Open("/scratch/Mainz_Software/a2GoAT/LinPol_Aug16.root"); // Open linear polarisation plot

    for (Int_t i = 0; i < 18; i++){

        SigmaValues[1][0][i] = pValues415[i+1]/(Graph->Eval(425,0));
        SigmaErrValues[1][0][i] = pErrValues415[i+1]/(Graph->Eval(425,0));
        SigmaValues[1][1][i] = pValues425[i+1]/(Graph->Eval(425,0));
        SigmaErrValues[1][1][i] = pErrValues425[i+1]/(Graph->Eval(425,0));
        SigmaValues[1][2][i] = pValues435[i+1]/(Graph->Eval(435,0));
        SigmaErrValues[1][2][i] = pErrValues435[i+1]/(Graph->Eval(435,0));
        SigmaValues[1][3][i] = pValues445[i+1]/(Graph->Eval(445,0));
        SigmaErrValues[1][3][i] = pErrValues445[i+1]/(Graph->Eval(445,0));
        SigmaValues[1][4][i] = pValues455[i+1]/(Graph->Eval(455,0));
        SigmaErrValues[1][4][i] = pErrValues455[i+1]/(Graph->Eval(455,0));
        SigmaValues[1][5][i] = pValues465[i+1]/(Graph->Eval(465,0));
        SigmaErrValues[1][5][i] = pErrValues465[i+1]/(Graph->Eval(465,0));
        SigmaValues[1][6][i] = pValues475[i+1]/(Graph->Eval(475,0));
        SigmaErrValues[1][6][i] = pErrValues475[i+1]/(Graph->Eval(475,0));
        SigmaValues[1][7][i] = pValues485[i+1]/(Graph->Eval(485,0));
        SigmaErrValues[1][7][i] = pErrValues485[i+1]/(Graph->Eval(485,0));
        SigmaValues[1][8][i] = pValues495[i+1]/(Graph->Eval(495,0));
        SigmaErrValues[1][8][i] = pErrValues495[i+1]/(Graph->Eval(495,0));
        SigmaValues[1][9][i] = pValues505[i+1]/(Graph->Eval(505,0));
        SigmaErrValues[1][9][i] = pErrValues505[i+1]/(Graph->Eval(505,0));
        SigmaValues[1][10][i] = pValues515[i+1]/(Graph->Eval(515,0));
        SigmaErrValues[1][10][i] = pErrValues515[i+1]/(Graph->Eval(515,0));
        SigmaValues[1][11][i] = pValues525[i+1]/(Graph->Eval(525,0));
        SigmaErrValues[1][11][i] = pErrValues525[i+1]/(Graph->Eval(525,0));
        SigmaValues[1][12][i] = pValues535[i+1]/(Graph->Eval(535,0));
        SigmaErrValues[1][12][i] = pErrValues535[i+1]/(Graph->Eval(535,0));
        SigmaValues[1][13][i] = pValues545[i+1]/(Graph->Eval(545,0));
        SigmaErrValues[1][13][i] = pErrValues545[i+1]/(Graph->Eval(545,0));
        SigmaValues[1][14][i] = pValues555[i+1]/(Graph->Eval(555,0));
        SigmaErrValues[1][14][i] = pErrValues555[i+1]/(Graph->Eval(555,0));
        SigmaValues[1][15][i] = pValues565[i+1]/(Graph->Eval(565,0));
        SigmaErrValues[1][15][i] = pErrValues565[i+1]/(Graph->Eval(565,0));
        SigmaValues[1][16][i] = pValues575[i+1]/(Graph->Eval(575,0));
        SigmaErrValues[1][16][i] = pErrValues575[i+1]/(Graph->Eval(575,0));
        SigmaValues[1][17][i] = pValues585[i+1]/(Graph->Eval(585,0));
        SigmaErrValues[1][17][i] = pErrValues585[i+1]/(Graph->Eval(585,0));
        SigmaValues[1][18][i] = pValues595[i+1]/(Graph->Eval(595,0));
        SigmaErrValues[1][18][i] = pErrValues595[i+1]/(Graph->Eval(595,0));
        SigmaValues[1][19][i] = pValues605[i+1]/(Graph->Eval(605,0));
        SigmaErrValues[1][19][i] = pErrValues605[i+1]/(Graph->Eval(605,0));
        SigmaValues[1][20][i] = pValues615[i+1]/(Graph->Eval(615,0));
        SigmaErrValues[1][20][i] = pErrValues615[i+1]/(Graph->Eval(615,0));

    }

    TFile f5("Sigma_Systematic_18_17.root", "RECREATE");

    TF1 *Line = new TF1("Line", "[0]", -1, 1);
    double x[18] = {0.85, 0.75, 0.65, 0.55, 0.45, 0.35, 0.25, 0.15, 0.05, -0.05, -0.15, -0.25, -0.35, -0.45, -0.55, -0.65, -0.75, -0.85}; // Need to adjust
    double ex[18] = {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05}; // Need to adjust

    for(Int_t i = 0; i < 21; i ++){ // Egamma value
        for(Int_t j = 0; j < 18; j ++){ // CM value
            SystValues[i][j] = SigmaValues[0][i][j] - SigmaValues[1][i][j];
            SystErrValues[i][j] = sqrt( ((SigmaErrValues[0][i][j])**2) + ((SigmaErrValues[1][i][j])**2) );
        }
        sprintf(name, "SigmaSyst_%i", 415+(i*10));
        sprintf(title, "#Sigma_{2#sigma} - #Sigma_{1#sigma} E_{#gamma} %i #pm 10 MeV", 415+(i*10));
        SigmaSystPlots[i] = new TGraphErrors(18 , x, SystValues[i] , ex, SystErrValues[i]) ;
        SigmaSystPlots[i]->Fit("Line", "M");
        SigmaSystPlots[i]->SetMarkerColor(4);
        SigmaSystPlots[i]->SetLineColor(4);
        SigmaSystPlots[i]->SetMarkerStyle(8);
        SigmaSystPlots[i]->SetMarkerSize(1);
        SigmaSystPlots[i]->GetXaxis()->SetTitle("Cos#theta_{CM}");
        SigmaSystPlots[i]->GetXaxis()->SetRangeUser(-1, 1);
        SigmaSystPlots[i]->GetYaxis()->SetTitle("#Sigma_{2#sigma} - #Sigma_{1#sigma}");
        SigmaSystPlots[i]->SetName(name);
        SigmaSystPlots[i]->SetTitle(title);
        SigmaSystPlots[i]->Write();

        LinePar[i] = Line->GetParameter(0);
        LineParErr[i] = Line->GetParError(0);

        cout << (LinePar[i])/(LineParErr[i]) << endl;
    }

    TCanvas *canvas20 = new TCanvas("canvas20","canvas20", 1920, 1080);
    canvas20->Divide(5,4);
    for(int i = 1 ; i < 22 ; i++){
        if(i == 15) continue;
        if (i < 15) canvas20->cd(i);
        else if (i > 15) canvas20->cd(i-1);
        SigmaSystPlots[i-1]->Draw("AEP");
//        TLine *Plus1Sigma = new TLine(-1, LinePar[i-1] + LineParErr[i-1], 1, LinePar[i-1] + LineParErr[i-1]);
//        TLine *Minus1Sigma = new TLine(-1, LinePar[i-1] - LineParErr[i-1], 1, LinePar[i-1] - LineParErr[i-1]);
//        Plus1Sigma->SetLineColor(2);
//        Plus1Sigma->SetLineStyle(9);
//        Plus1Sigma->Draw("SAME");
//        Minus1Sigma->SetLineColor(2);
//        Minus1Sigma->SetLineStyle(9);
//        Minus1Sigma->Draw("SAME");
    }

    canvas20->Write();

    f5.Write();

}
