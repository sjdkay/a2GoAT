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

double pCosAmp[21][20]; // Format of array is Egamma bin (x) Theta bin (y)
double pCosAmpErr[21][20];

Double_t LegPar[8][21];
Double_t LegParErr[8][21];
Double_t SigmaValues[21][18];
Double_t SigmaErrValues[21][18];
double p0, p0Err, p1, p1Err, p2, p2Err, p3, p3Err, p4, p4Err, p5, p5Err, p6, p6Err, p7, p7Err;
