#include "TROOT.h"
#include <TChain.h>
#include "TFile.h"

#include "TObject.h"
#include "TApplication.h"
#include "TVector3.h"
#include "TVector2.h"
#include "TString.h"
#include "TChain.h"
#include "TGraphErrors.h"
#include "TLorentzVector.h"

#include "TStyle.h"
#include "TColor.h"
#include "TLegend.h"

#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TMath.h"

#include "TFile.h"
#include "TTree.h"

#include <stdio.h>
#include <iostream>
#include <unistd.h>

#include <fstream>

#include <string.h>

using namespace std;
using namespace TMath;

    Double_t NPara;
    Double_t NPerp;
    Double_t ScaleFactor;
    Double_t ScaleFactorErr;

    double pCosAmp[21][20]; // Format of array is Egamma bin (x) Theta bin (y)
    double pCosAmpErr[21][20];
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



