#include "./includes_ParaPerp_Combiner_Unpolarised_V2.h"

void ParaPerp_Combiner_Unpolarised_V2() {

    char name[60];
    char title[60];
    char PosHelHistName[60];
    char NegHelHistName[60];

    TFile *f = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/Para/Physics_Total_Para_34_14_5_18.root"); // Open latest Para file

    TH1D* Eg_Para = (TH1D*)f->Get("Eg")->Clone();
    Eg_Para->SetName("Eg_Para");

    for(Int_t A = 0; A < 6; A++){
        for(Int_t B = 0; B < 3; B++){
            sprintf(PosHelHistName, "PhiSc%iPosHelCM%i", 200+(A*100), B+1);
            sprintf(NegHelHistName, "PhiSc%iNegHelCM%i", 200+(A*100), B+1);
            PhiScPosHelPara[A][B] = (TH1D*) (((TH1D*)f->Get(PosHelHistName)));
            PhiScNegHelPara[A][B] = (TH1D*) (((TH1D*)f->Get(NegHelHistName)));
        }
    }

    TFile *f1 = new TFile("/scratch/Mainz_Software/Data/GoAT_Output/GoAT_23_01_17/Perp/Physics_Total_Perp_34_14_5_18.root"); // Open latest Perp file

    TH1D* Eg_Perp = (TH1D*)f1->Get("Eg")->Clone();
    Eg_Perp->SetName("Eg_Perp");

    for(Int_t A = 0; A < 6; A++){
        for(Int_t B = 0; B < 3; B++){
            sprintf(PosHelHistName, "PhiSc%iPosHelCM%i", 200+(A*100), B+1);
            sprintf(NegHelHistName, "PhiSc%iNegHelCM%i", 200+(A*100), B+1);
            PhiScPosHelPerp[A][B] = (TH1D*) (((TH1D*)f1->Get(PosHelHistName)));
            PhiScNegHelPerp[A][B] = (TH1D*) (((TH1D*)f1->Get(NegHelHistName)));
        }
    }

    NPara = Eg_Para->GetEntries();
    NPerp = Eg_Perp->GetEntries();
    ScaleFactor = NPara/NPerp;
    ScaleFactorErr = sqrt( (NPara/((TMath::Power(NPerp,2)))) + (((TMath::Power(NPara,2)))/(TMath::Power(NPerp,3))) ); // Error Propagation of above formula, see notebook 5

    Eg = new TH1D( "Eg", "E_{#gamma} Distribution", 200, 100, 1600 );

    for(Int_t A = 0; A < 6; A++){
        for(Int_t B = 0; B < 3; B++){
            PhiScPosHel[A][B] = new TH1D(Form("PhiSc%iPosHelCM%i", 200+(A*100), B+1), Form("#phi_{Sc} E_{#gamma}%i #pm 50MeV CM%i for +ve Helicity", 200+(A*100), B+1), 10, -4, 4);
            PhiScNegHel[A][B] = new TH1D(Form("PhiSc%iNegHelCM%i", 200+(A*100), B+1), Form("#phi_{Sc} E_{#gamma}%i #pm 50MeV CM%i for -ve Helicity", 200+(A*100), B+1), 10, -4, 4);
        }
    }

    Eg->Add(Eg_Para);
    Eg_Perp->Scale(ScaleFactor);
    Eg->Add(Eg_Perp);

    for(Int_t A = 0; A < 6; A++){
        for(Int_t B = 0; B < 3; B++){
            PhiScPosHel[A][B]->Add(PhiScPosHelPara[A][B]);
            PhiScPosHelPerp[A][B]->Scale(ScaleFactor);
            PhiScPosHel[A][B]->Add(PhiScPosHelPerp[A][B]);

            PhiScNegHel[A][B]->Add(PhiScNegHelPara[A][B]);
            PhiScNegHelPerp[A][B]->Scale(ScaleFactor);
            PhiScNegHel[A][B]->Add(PhiScNegHelPerp[A][B]);
        }
    }

    TFile f2("ParaPerp_Total_34_Combined_Unpolarised.root", "RECREATE");

    Eg->Write();

    for(Int_t A = 0; A < 6; A++){
        for(Int_t B = 0; B < 3; B++){
            PhiScPosHel[A][B]->Write();
            PhiScNegHel[A][B]->Write();
        }
    }

    f2.Write();
}
