#include "./includes.h"

void CircPol(){

   TCanvas *canvas = new TCanvas("canvas","canvas", 1920, 1080);
   Double_t x[1257], y[1257];
   Double_t Pm = 0.7666;
   Double_t Em = 1557;
   Int_t n = 1257;
   for (Int_t i=0; i<n; i++) {
     x[i] = i+300;
     y[i] = (Pm*x[i]*(Em+((1/3)*(Pm-x[i]))))/((TMath::Power(Em,2))+(TMath::Power((Em-x[i]),2))-((2/3)*Em*(Em-x[i])));
   }
   TGraph* CircPolGamma = new TGraph(n,x,y);
   CircPolGamma->SetTitle("p^{#odot}_{#gamma}(E_{#gamma})");
   CircPolGamma->GetXaxis()->SetTitle("E_{#gamma}/MeV");
   CircPolGamma->GetYaxis()->SetTitle("p^{#odot}_{#gamma}");
   CircPolGamma->GetXaxis()->SetLimits(300, 1557);
   CircPolGamma->SetLineWidth(4);
   CircPolGamma->Draw("AC");
   CircPolGamma->SaveAs("CircPol_Aug16.root");
}
