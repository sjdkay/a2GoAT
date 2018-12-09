#include "./includes_ParaPerp_Combiner_Unpolarised_V3.h"

// Version 3 utilises the 2D histograms of AEff as a fn of PhiSc for each bin
// AEff is extracted from this and averaged to calculate the value of Cx'
// Obsolete Histograms are no longer added as neccessary


void ParaPerp_Combiner_Unpolarised_V3() {

    char name[60];
    char title[60];
    char PhiScHistName[60];
    char PosHelHistName[60];
    char NegHelHistName[60];
    char PosHelHistName2[60];
    char NegHelHistName2[60];

    TFile *f = new TFile("/d4tb1/sjdkay/MainzRecoil/GoAT_Output_Files/nobackup/GoAT_Files/Para/Physics_Total_Para_54_40MeV_CentBins.root"); // Open latest Para file

    TH1D* Eg_Para = (TH1D*)f->Get("Eg")->Clone();
    Eg_Para->SetName("Eg_Para");

    for(Int_t A = 0; A < 8; A++){
        for(Int_t B = 0; B < 5; B++){
            sprintf(PosHelHistName, "PhiSc%iPosHelCM%i", 340+(A*100), B+1);
            sprintf(NegHelHistName, "PhiSc%iNegHelCM%i", 340+(A*100), B+1);
            sprintf(PosHelHistName2, "PhiScAEff%iPosHelCM%i", 340+(A*100), B+1);
            sprintf(NegHelHistName2, "PhiScAEff%iNegHelCM%i", 340+(A*100), B+1);
            PhiScPosHelPara[A][B] = (TH1D*) (((TH1D*)f->Get(PosHelHistName)));
            PhiScNegHelPara[A][B] = (TH1D*) (((TH1D*)f->Get(NegHelHistName)));
            PhiScAEffPosHelPara[A][B] = (TH2D*) (((TH2D*)f->Get(PosHelHistName2)));
            PhiScAEffNegHelPara[A][B] = (TH2D*) (((TH2D*)f->Get(NegHelHistName2)));
        }
    }

    TFile *f1 = new TFile("/d4tb1/sjdkay/MainzRecoil/GoAT_Output_Files/nobackup/GoAT_Files/Perp/Physics_Total_Perp_54_40MeV_CentBins.root"); // Open latest Perp file

    TH1D* Eg_Perp = (TH1D*)f1->Get("Eg")->Clone();
    Eg_Perp->SetName("Eg_Perp");

    for(Int_t A = 0; A < 8; A++){
        for(Int_t B = 0; B < 5; B++){
            sprintf(PosHelHistName, "PhiSc%iPosHelCM%i", 340+(A*100), B+1);
            sprintf(NegHelHistName, "PhiSc%iNegHelCM%i", 340+(A*100), B+1);
            sprintf(PosHelHistName2, "PhiScAEff%iPosHelCM%i", 340+(A*100), B+1);
            sprintf(NegHelHistName2, "PhiScAEff%iNegHelCM%i", 340+(A*100), B+1);
            PhiScPosHelPerp[A][B] = (TH1D*) (((TH1D*)f1->Get(PosHelHistName)));
            PhiScNegHelPerp[A][B] = (TH1D*) (((TH1D*)f1->Get(NegHelHistName)));
            PhiScAEffPosHelPerp[A][B] = (TH2D*) (((TH2D*)f1->Get(PosHelHistName2)));
            PhiScAEffNegHelPerp[A][B] = (TH2D*) (((TH2D*)f1->Get(NegHelHistName2)));
        }
    }

    TFile f2("ParaPerp_Total_54_Combined_Unpolarised_40MeV_Cent.root", "RECREATE");

    NPara = Eg_Para->GetEntries();
    NPerp = Eg_Perp->GetEntries();
    ScaleFactor = NPara/NPerp;
    ScaleFactorErr = sqrt( (NPara/((TMath::Power(NPerp,2)))) + (((TMath::Power(NPara,2)))/(TMath::Power(NPerp,3))) ); // Error Propagation of above formula, see notebook 5

    Eg = new TH1D( "Eg", "E_{#gamma} Distribution", 200, 100, 1600 );

    for(Int_t A = 0; A < 8; A++){
        for(Int_t B = 0; B < 5; B++){
            PhiScPosHel[A][B] = new TH1D(Form("PhiSc%iPosHelCM%i", 340+(A*100), B+1), Form("#phi_{Sc} E_{#gamma}%i #pm 50MeV CM%i for +ve Helicity", 340+(A*100), B+1), 10, -1*acos(-1), acos(-1));
            PhiScNegHel[A][B] = new TH1D(Form("PhiSc%iNegHelCM%i", 340+(A*100), B+1), Form("#phi_{Sc} E_{#gamma}%i #pm 50MeV CM%i for -ve Helicity", 340+(A*100), B+1), 10, -1*acos(-1), acos(-1));
            PhiScAEffPosHel[A][B] = new TH2D(Form("PhiScAEff%iPosHelCM%i", 340+(A*100), B+1), Form("#phi_{Sc} E_{#gamma}%i #pm 50MeV CM%i for +ve Helicity", 340+(A*100), B+1), 10, -1*acos(-1), acos(-1), 500, -1, 1);
            PhiScAEffNegHel[A][B] = new TH2D(Form("PhiScAEff%iNegHelCM%i", 340+(A*100), B+1), Form("#phi_{Sc} E_{#gamma}%i #pm 50MeV CM%i for -ve Helicity", 340+(A*100), B+1), 10, -1*acos(-1), acos(-1), 500, -1, 1);
        }
    }

    Eg->Add(Eg_Para);
    Eg_Perp->Scale(ScaleFactor);
    Eg->Add(Eg_Perp);

    Eg->Write();

    for(Int_t A = 0; A < 8; A++){
        for(Int_t B = 0; B < 5; B++){
            PhiScPosHel[A][B]->Add(PhiScPosHelPara[A][B]);
            PhiScPosHelPerp[A][B]->Scale(ScaleFactor);
            PhiScPosHel[A][B]->Add(PhiScPosHelPerp[A][B]);

            PhiScNegHel[A][B]->Add(PhiScNegHelPara[A][B]);
            PhiScNegHelPerp[A][B]->Scale(ScaleFactor);
            PhiScNegHel[A][B]->Add(PhiScNegHelPerp[A][B]);

            PhiScAEffPosHel[A][B]->Add(PhiScAEffPosHelPara[A][B]);
            PhiScAEffPosHelPerp[A][B]->Scale(ScaleFactor);
            PhiScAEffPosHel[A][B]->Add(PhiScAEffPosHelPerp[A][B]);

            PhiScAEffNegHel[A][B]->Add(PhiScAEffNegHelPara[A][B]);
            PhiScAEffNegHelPerp[A][B]->Scale(ScaleFactor);
            PhiScAEffNegHel[A][B]->Add(PhiScAEffNegHelPerp[A][B]);

            PhiScPosHel[A][B]->Write();
            PhiScNegHel[A][B]->Write();
            PhiScAEffPosHel[A][B]->Write();
            PhiScAEffNegHel[A][B]->Write();

        }
    }

    f2.Write();
}
