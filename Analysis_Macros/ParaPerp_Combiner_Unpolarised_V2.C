#include "./includes_ParaPerp_Combiner_Unpolarised_V2.h"

void ParaPerp_Combiner_Unpolarised_V2() {

    char name[60];
    char title[60];
    char PhiScHistName[60];
    char PosHelHistName[60];
    char NegHelHistName[60];
    char NeutronEThetaScName[60];

    TFile *f = new TFile("/d4tb1/sjdkay/MainzRecoil/GoAT_Output_Files/nobackup/GoAT_Files/Para/Physics_Total_Para_53_1_10_18.root"); // Open latest Para file

    TH1D* Eg_Para = (TH1D*)f->Get("Eg")->Clone();
    Eg_Para->SetName("Eg_Para");

    for(Int_t A = 0; A < 8; A++){
        for(Int_t B = 0; B < 3; B++){
            sprintf(PosHelHistName, "PhiSc%iPosHelCM%i", 250+(A*100), B+1);
            sprintf(NegHelHistName, "PhiSc%iNegHelCM%i", 250+(A*100), B+1);
            sprintf(NeutronEThetaScName, "NeutronEThetaSc%iCM%i", 250+(A*100), B+1);
            PhiScPosHelPara[A][B] = (TH1D*) (((TH1D*)f->Get(PosHelHistName)));
            PhiScNegHelPara[A][B] = (TH1D*) (((TH1D*)f->Get(NegHelHistName)));
            NeutronEThetaScPara[A][B] = (TH2D*) ((TH2D*)f->Get(NeutronEThetaScName)));
        }
    }

//    for(Int_t A = 0; A < 12; A++){
//        sprintf(PhiScHistName, "PhiSc%i", 200+(A*50));
//        PhiScPara[A] = (TH1D*) (((TH1D*)f->Get(PhiScHistName)));
//    }

    TFile *f1 = new TFile("/d4tb1/sjdkay/MainzRecoil/GoAT_Output_Files/nobackup/GoAT_Files/Perp/Physics_Total_Perp_53_1_10_18.root"); // Open latest Perp file

    TH1D* Eg_Perp = (TH1D*)f1->Get("Eg")->Clone();
    Eg_Perp->SetName("Eg_Perp");

    for(Int_t A = 0; A < 8; A++){
        for(Int_t B = 0; B < 3; B++){
            sprintf(PosHelHistName, "PhiSc%iPosHelCM%i", 250+(A*100), B+1);
            sprintf(NegHelHistName, "PhiSc%iNegHelCM%i", 250+(A*100), B+1);
            sprintf(NeutronEThetaScName, "NeutronEThetaSc%iCM%i", 250+(A*100), B+1);
            PhiScPosHelPerp[A][B] = (TH1D*) (((TH1D*)f1->Get(PosHelHistName)));
            PhiScNegHelPerp[A][B] = (TH1D*) (((TH1D*)f1->Get(NegHelHistName)));
            NeutronEThetaScPerp[A][B] = (TH2D*) ((TH2D*)f1->Get(NeutronEThetaScName)));
        }
    }

//    for(Int_t A = 0; A < 12; A++){
//        sprintf(PhiScHistName, "PhiSc%i", 200+(A*50));
//        PhiScPerp[A] = (TH1D*) (((TH1D*)f1->Get(PhiScHistName)));
//    }

    TFile f2("ParaPerp_Total_53_Combined_Unpolarised.root", "RECREATE");

    NPara = Eg_Para->GetEntries();
    NPerp = Eg_Perp->GetEntries();
    ScaleFactor = NPara/NPerp;
    ScaleFactorErr = sqrt( (NPara/((TMath::Power(NPerp,2)))) + (((TMath::Power(NPara,2)))/(TMath::Power(NPerp,3))) ); // Error Propagation of above formula, see notebook 5

    Eg = new TH1D( "Eg", "E_{#gamma} Distribution", 200, 100, 1600 );

    for(Int_t A = 0; A < 8; A++){
        for(Int_t B = 0; B < 3; B++){
            PhiScPosHel[A][B] = new TH1D(Form("PhiSc%iPosHelCM%i", 250+(A*100), B+1), Form("#phi_{Sc} E_{#gamma}%i #pm 50MeV CM%i for +ve Helicity", 250+(A*100), B+1), 10, -1*acos(-1), acos(-1));
            PhiScNegHel[A][B] = new TH1D(Form("PhiSc%iNegHelCM%i", 250+(A*100), B+1), Form("#phi_{Sc} E_{#gamma}%i #pm 50MeV CM%i for -ve Helicity", 250+(A*100), B+1), 10, -1*acos(-1), acos(-1));
            NeutronEThetaSc[A][B] = new TH2D(Form ("NeutronEThetaSc%iCM%i", 250+(A*100), B+1), Form ("#theta_{Sc}(E_{n}) %iMeV CM%i", 250+(A*100), B+1), 200, 0, 800, 200, 0, 90);
        }
    }

//    for(Int_t X = 0; X < 12; X++){
//        PhiSc[X] = new TH1D(Form("PhiSc%i", 200+(X*50)), Form("#phi_{Sc} E_{#gamma}%i #pm 25MeV", 200+(X*50)), 10, -4, 4);
//    }

    Eg->Add(Eg_Para);
    Eg_Perp->Scale(ScaleFactor);
    Eg->Add(Eg_Perp);

    Eg->Write();

    for(Int_t A = 0; A < 8; A++){
        for(Int_t B = 0; B < 3; B++){
            PhiScPosHel[A][B]->Add(PhiScPosHelPara[A][B]);
            PhiScPosHelPerp[A][B]->Scale(ScaleFactor);
            PhiScPosHel[A][B]->Add(PhiScPosHelPerp[A][B]);

            PhiScNegHel[A][B]->Add(PhiScNegHelPara[A][B]);
            PhiScNegHelPerp[A][B]->Scale(ScaleFactor);
            PhiScNegHel[A][B]->Add(PhiScNegHelPerp[A][B]);

            NeutronEThetaSc[A][B]->Add(NeutronEThetaScPara[A][B]);
            NeutronEThetaScPerp[A][B]->Scale(ScaleFactor);
            NeutronEThetaSc[A][B]->Add(NeutronEThetaScPerp[A][B]);

            PhiScPosHel[A][B]->Write();
            PhiScNegHel[A][B]->Write();
            NeutronEThetaSc[A][B]->Write();
        }
    }

//    for(Int_t X = 0; X < 12; X++){
//        PhiSc[X]->Add(PhiScPara[X]);
//        PhiScPerp[X]->Scale(ScaleFactor);
//        PhiSc[X]->Add(PhiScPerp[X]);
//        PhiSc[X]->Write();
//    }

    f2.Write();
}
