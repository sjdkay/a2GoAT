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

    double pCosAmp[10][5]; // Format of array is Egamma bin (x) Theta bin (y)
    double pCosAmpErr[10][5];
    double pCosA410;
    double pCosAErr410;
    double pCosA430;
    double pCosAErr430;
    double pCosA450;
    double pCosAErr450;
    double pCosA470;
    double pCosAErr470;
    double pCosA490;
    double pCosAErr490;
    double pCosA510;
    double pCosAErr510;
    double pCosA530;
    double pCosAErr530;
    double pCosA550;
    double pCosAErr550;
    double pCosA570;
    double pCosAErr570;
    double pCosA590;
    double pCosAErr590;

