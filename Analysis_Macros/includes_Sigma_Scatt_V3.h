#include "TROOT.h"
#include <TChain.h>
#include "TFile.h"
#include "TList.h"

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

Double_t SigmaValues[10][5];
Double_t SigmaErrValues[10][5];
