#include "./includes.h"

// define a function with 4 parameters
Double_t fitf(Double_t *x,Double_t *par)
{
    Double_t fitval = 0;
    fitval =  ((par[0]*sin(x[0]))/(1 + (par[1]*cos(x[0]))));
    return fitval;
}

void CxAsymm_V2() {

    double InitialSinAmp[8][8];
    double InitialSinAmpErr[8][8];
    double SinAmp[8][8]; // Format of arrays is Theta Bin Defined by i, Energy bin defined by j
    double SinAmpErr[8][8];
    double CosAmp[8][8]; // Format of arrays is Theta Bin Defined by i, Energy bin defined by j
    double CosAmpErr[8][8];
    Int_t i;
    double ISinAm265;
    double ISinAmErr265;
    double SinAm265;
    double SinAmErr265;
    double CosAm265;
    double CosAmErr265;
    double ISinAm335;
    double ISinAmErr335;
    double SinAm335;
    double SinAmErr335;
    double CosAm335;
    double CosAmErr335;
    double ISinAm405;
    double ISinAmErr405;
    double SinAm405;
    double SinAmErr405;
    double CosAm405;
    double CosAmErr405;
    double ISinAm475;
    double ISinAmErr475;
    double SinAm475;
    double SinAmErr475;
    double CosAm475;
    double CosAmErr475;
    double ISinAm545;
    double ISinAmErr545;
    double SinAm545;
    double SinAmErr545;
    double CosAm545;
    double CosAmErr545;
    double ISinAm615;
    double ISinAmErr615;
    double SinAm615;
    double SinAmErr615;
    double CosAm615;
    double CosAmErr615;
    double ISinAm685;
    double ISinAmErr685;
    double SinAm685;
    double SinAmErr685;
    double CosAm685;
    double CosAmErr685;

    Double_t Cx265[8], Cx335[8], Cx405[8], Cx475[8], Cx545[8], Cx615[8], Cx685[8];
    Double_t CxErr265[8], CxErr335[8], CxErr405[8], CxErr475[8], CxErr545[8], CxErr615[8], CxErr685[8];

    TF1 *AsymmFunc = new TF1("AsymmFit",  fitf, -3.0, 3.0, 2); //Give a name and range to the fitting funcion
    AsymmFunc->SetParNames("SinAmp", "CosAmp"); //Name the parameters
    AsymmFunc->SetParameter(0, 0);
    TF1 *SinFunc = new TF1("SinFit", "[0]*sin(x*TMath::DegToRad())", -3, 3);
    SinFunc->SetParNames("InitialSinAmp");
    TFile *f = new TFile("/scratch/Mainz_Software/a2GoAT/AmoTotal_2_96_27.root"); // Open the latest PTotal file to load histograms from
    TF1 *Pn90CM = new TF1("Pn90CM", "1.64576-2.95484*(x/1000)+0.684577*(x/1000)**2-0.65*90**2/4/((x-560)**2+90**2/4)+(5.32305-35.3819*(x/1000)+70.145*(x/1000)**2-44.2899*(x/1000)**3)",300,700);

    TList *PhiSc265NegHelList = new TList;
    TList *PhiSc265PosHelList = new TList;
    TList *PhiSc265AsymmList = new TList;
    TList *PhiSc335NegHelList = new TList;
    TList *PhiSc335PosHelList = new TList;
    TList *PhiSc335AsymmList = new TList;
    TList *PhiSc405NegHelList = new TList;
    TList *PhiSc405PosHelList = new TList;
    TList *PhiSc405AsymmList = new TList;
    TList *PhiSc475NegHelList = new TList;
    TList *PhiSc475PosHelList = new TList;
    TList *PhiSc475AsymmList = new TList;
    TList *PhiSc545NegHelList = new TList;
    TList *PhiSc545PosHelList = new TList;
    TList *PhiSc545AsymmList = new TList;
    TList *PhiSc615NegHelList = new TList;
    TList *PhiSc615PosHelList = new TList;
    TList *PhiSc615AsymmList = new TList;
    TList *PhiSc685NegHelList = new TList;
    TList *PhiSc685PosHelList = new TList;
    TList *PhiSc685AsymmList = new TList;

    PhiSc265NegHelList->Add(Phi_Scattered_265MeV_NegHelCM1);
    PhiSc265NegHelList->Add(Phi_Scattered_265MeV_NegHelCM2);
    PhiSc265NegHelList->Add(Phi_Scattered_265MeV_NegHelCM3);
    PhiSc265NegHelList->Add(Phi_Scattered_265MeV_NegHelCM4);
    PhiSc265NegHelList->Add(Phi_Scattered_265MeV_NegHelCM5);
    PhiSc265NegHelList->Add(Phi_Scattered_265MeV_NegHelCM6);
    PhiSc265NegHelList->Add(Phi_Scattered_265MeV_NegHelCM7);
    PhiSc265NegHelList->Add(Phi_Scattered_265MeV_NegHelCM8);
    PhiSc265PosHelList->Add(Phi_Scattered_265MeV_PosHelCM1);
    PhiSc265PosHelList->Add(Phi_Scattered_265MeV_PosHelCM2);
    PhiSc265PosHelList->Add(Phi_Scattered_265MeV_PosHelCM3);
    PhiSc265PosHelList->Add(Phi_Scattered_265MeV_PosHelCM4);
    PhiSc265PosHelList->Add(Phi_Scattered_265MeV_PosHelCM5);
    PhiSc265PosHelList->Add(Phi_Scattered_265MeV_PosHelCM6);
    PhiSc265PosHelList->Add(Phi_Scattered_265MeV_PosHelCM7);
    PhiSc265PosHelList->Add(Phi_Scattered_265MeV_PosHelCM8);

    PhiSc335NegHelList->Add(Phi_Scattered_335MeV_NegHelCM1);
    PhiSc335NegHelList->Add(Phi_Scattered_335MeV_NegHelCM2);
    PhiSc335NegHelList->Add(Phi_Scattered_335MeV_NegHelCM3);
    PhiSc335NegHelList->Add(Phi_Scattered_335MeV_NegHelCM4);
    PhiSc335NegHelList->Add(Phi_Scattered_335MeV_NegHelCM5);
    PhiSc335NegHelList->Add(Phi_Scattered_335MeV_NegHelCM6);
    PhiSc335NegHelList->Add(Phi_Scattered_335MeV_NegHelCM7);
    PhiSc335NegHelList->Add(Phi_Scattered_335MeV_NegHelCM8);
    PhiSc335PosHelList->Add(Phi_Scattered_335MeV_PosHelCM1);
    PhiSc335PosHelList->Add(Phi_Scattered_335MeV_PosHelCM2);
    PhiSc335PosHelList->Add(Phi_Scattered_335MeV_PosHelCM3);
    PhiSc335PosHelList->Add(Phi_Scattered_335MeV_PosHelCM4);
    PhiSc335PosHelList->Add(Phi_Scattered_335MeV_PosHelCM5);
    PhiSc335PosHelList->Add(Phi_Scattered_335MeV_PosHelCM6);
    PhiSc335PosHelList->Add(Phi_Scattered_335MeV_PosHelCM7);
    PhiSc335PosHelList->Add(Phi_Scattered_335MeV_PosHelCM8);

    PhiSc405NegHelList->Add(Phi_Scattered_405MeV_NegHelCM1);
    PhiSc405NegHelList->Add(Phi_Scattered_405MeV_NegHelCM2);
    PhiSc405NegHelList->Add(Phi_Scattered_405MeV_NegHelCM3);
    PhiSc405NegHelList->Add(Phi_Scattered_405MeV_NegHelCM4);
    PhiSc405NegHelList->Add(Phi_Scattered_405MeV_NegHelCM5);
    PhiSc405NegHelList->Add(Phi_Scattered_405MeV_NegHelCM6);
    PhiSc405NegHelList->Add(Phi_Scattered_405MeV_NegHelCM7);
    PhiSc405NegHelList->Add(Phi_Scattered_405MeV_NegHelCM8);
    PhiSc405PosHelList->Add(Phi_Scattered_405MeV_PosHelCM1);
    PhiSc405PosHelList->Add(Phi_Scattered_405MeV_PosHelCM2);
    PhiSc405PosHelList->Add(Phi_Scattered_405MeV_PosHelCM3);
    PhiSc405PosHelList->Add(Phi_Scattered_405MeV_PosHelCM4);
    PhiSc405PosHelList->Add(Phi_Scattered_405MeV_PosHelCM5);
    PhiSc405PosHelList->Add(Phi_Scattered_405MeV_PosHelCM6);
    PhiSc405PosHelList->Add(Phi_Scattered_405MeV_PosHelCM7);
    PhiSc405PosHelList->Add(Phi_Scattered_405MeV_PosHelCM8);

    PhiSc475NegHelList->Add(Phi_Scattered_475MeV_NegHelCM1);
    PhiSc475NegHelList->Add(Phi_Scattered_475MeV_NegHelCM2);
    PhiSc475NegHelList->Add(Phi_Scattered_475MeV_NegHelCM3);
    PhiSc475NegHelList->Add(Phi_Scattered_475MeV_NegHelCM4);
    PhiSc475NegHelList->Add(Phi_Scattered_475MeV_NegHelCM5);
    PhiSc475NegHelList->Add(Phi_Scattered_475MeV_NegHelCM6);
    PhiSc475NegHelList->Add(Phi_Scattered_475MeV_NegHelCM7);
    PhiSc475NegHelList->Add(Phi_Scattered_475MeV_NegHelCM8);
    PhiSc475PosHelList->Add(Phi_Scattered_475MeV_PosHelCM1);
    PhiSc475PosHelList->Add(Phi_Scattered_475MeV_PosHelCM2);
    PhiSc475PosHelList->Add(Phi_Scattered_475MeV_PosHelCM3);
    PhiSc475PosHelList->Add(Phi_Scattered_475MeV_PosHelCM4);
    PhiSc475PosHelList->Add(Phi_Scattered_475MeV_PosHelCM5);
    PhiSc475PosHelList->Add(Phi_Scattered_475MeV_PosHelCM6);
    PhiSc475PosHelList->Add(Phi_Scattered_475MeV_PosHelCM7);
    PhiSc475PosHelList->Add(Phi_Scattered_475MeV_PosHelCM8);

    PhiSc545NegHelList->Add(Phi_Scattered_545MeV_NegHelCM1);
    PhiSc545NegHelList->Add(Phi_Scattered_545MeV_NegHelCM2);
    PhiSc545NegHelList->Add(Phi_Scattered_545MeV_NegHelCM3);
    PhiSc545NegHelList->Add(Phi_Scattered_545MeV_NegHelCM4);
    PhiSc545NegHelList->Add(Phi_Scattered_545MeV_NegHelCM5);
    PhiSc545NegHelList->Add(Phi_Scattered_545MeV_NegHelCM6);
    PhiSc545NegHelList->Add(Phi_Scattered_545MeV_NegHelCM7);
    PhiSc545NegHelList->Add(Phi_Scattered_545MeV_NegHelCM8);
    PhiSc545PosHelList->Add(Phi_Scattered_545MeV_PosHelCM1);
    PhiSc545PosHelList->Add(Phi_Scattered_545MeV_PosHelCM2);
    PhiSc545PosHelList->Add(Phi_Scattered_545MeV_PosHelCM3);
    PhiSc545PosHelList->Add(Phi_Scattered_545MeV_PosHelCM4);
    PhiSc545PosHelList->Add(Phi_Scattered_545MeV_PosHelCM5);
    PhiSc545PosHelList->Add(Phi_Scattered_545MeV_PosHelCM6);
    PhiSc545PosHelList->Add(Phi_Scattered_545MeV_PosHelCM7);
    PhiSc545PosHelList->Add(Phi_Scattered_545MeV_PosHelCM8);

    PhiSc615NegHelList->Add(Phi_Scattered_615MeV_NegHelCM1);
    PhiSc615NegHelList->Add(Phi_Scattered_615MeV_NegHelCM2);
    PhiSc615NegHelList->Add(Phi_Scattered_615MeV_NegHelCM3);
    PhiSc615NegHelList->Add(Phi_Scattered_615MeV_NegHelCM4);
    PhiSc615NegHelList->Add(Phi_Scattered_615MeV_NegHelCM5);
    PhiSc615NegHelList->Add(Phi_Scattered_615MeV_NegHelCM6);
    PhiSc615NegHelList->Add(Phi_Scattered_615MeV_NegHelCM7);
    PhiSc615NegHelList->Add(Phi_Scattered_615MeV_NegHelCM8);
    PhiSc615PosHelList->Add(Phi_Scattered_615MeV_PosHelCM1);
    PhiSc615PosHelList->Add(Phi_Scattered_615MeV_PosHelCM2);
    PhiSc615PosHelList->Add(Phi_Scattered_615MeV_PosHelCM3);
    PhiSc615PosHelList->Add(Phi_Scattered_615MeV_PosHelCM4);
    PhiSc615PosHelList->Add(Phi_Scattered_615MeV_PosHelCM5);
    PhiSc615PosHelList->Add(Phi_Scattered_615MeV_PosHelCM6);
    PhiSc615PosHelList->Add(Phi_Scattered_615MeV_PosHelCM7);
    PhiSc615PosHelList->Add(Phi_Scattered_615MeV_PosHelCM8);

    PhiSc685NegHelList->Add(Phi_Scattered_685MeV_NegHelCM1);
    PhiSc685NegHelList->Add(Phi_Scattered_685MeV_NegHelCM2);
    PhiSc685NegHelList->Add(Phi_Scattered_685MeV_NegHelCM3);
    PhiSc685NegHelList->Add(Phi_Scattered_685MeV_NegHelCM4);
    PhiSc685NegHelList->Add(Phi_Scattered_685MeV_NegHelCM5);
    PhiSc685NegHelList->Add(Phi_Scattered_685MeV_NegHelCM6);
    PhiSc685NegHelList->Add(Phi_Scattered_685MeV_NegHelCM7);
    PhiSc685NegHelList->Add(Phi_Scattered_685MeV_NegHelCM8);
    PhiSc685PosHelList->Add(Phi_Scattered_685MeV_PosHelCM1);
    PhiSc685PosHelList->Add(Phi_Scattered_685MeV_PosHelCM2);
    PhiSc685PosHelList->Add(Phi_Scattered_685MeV_PosHelCM3);
    PhiSc685PosHelList->Add(Phi_Scattered_685MeV_PosHelCM4);
    PhiSc685PosHelList->Add(Phi_Scattered_685MeV_PosHelCM5);
    PhiSc685PosHelList->Add(Phi_Scattered_685MeV_PosHelCM6);
    PhiSc685PosHelList->Add(Phi_Scattered_685MeV_PosHelCM7);
    PhiSc685PosHelList->Add(Phi_Scattered_685MeV_PosHelCM8);



}
