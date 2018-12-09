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
    Double_t NPara2;
    Double_t NPerp2;
    Double_t ScaleFactor2;
    Double_t ScaleFactorErr2;

    double pCosAmp[21][20]; // Format of array is Egamma bin (x) Theta bin (y)
    double pCosAmpErr[21][20];
    double pCosAmp2[21][20];
    double pCosAmpErr2[21][20];

    Double_t SigmaValues[2][21][18];
    Double_t SigmaErrValues[2][21][18];
    Double_t SystValues[21][18];
    Double_t SystErrValues[21][18];
    Double_t LinePar[21];
    Double_t LineParErr[21];
