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

Double_t RebinVal;
Double_t BinWidth;
double Offset[4][13]; // Format of arrays is Theta Bin Defined by i, Energy bin defined by j
double OffsetErr[4][13];
double SinAmp[4][13]; // Format of arrays is Theta Bin Defined by i, Energy bin defined by j
double SinAmpErr[4][13];
double CosAmp[4][13]; // Format of arrays is Theta Bin Defined by i, Energy bin defined by j
double CosAmpErr[4][13];
Int_t i;
double OffsetAllTheta;
double OffsetAllThetaErr;
double OffsetTheta010;
double OffsetTheta010Err;
double OffsetTheta1020;
double OffsetTheta1020Err;
double OffsetTheta2030;
double OffsetTheta2030Err;
double SinAmpAllTheta;
double SinAmpAllThetaErr;
double SinAmpTheta010;
double SinAmpTheta010Err;
double SinAmpTheta1020;
double SinAmpTheta1020Err;
double SinAmpTheta2030;
double SinAmpTheta2030Err;
double CosAmpAllTheta;
double CosAmpAllThetaErr;
double CosAmpTheta010;
double CosAmpTheta010Err;
double CosAmpTheta1020;
double CosAmpTheta1020Err;
double CosAmpTheta2030;
double CosAmpTheta2030Err;
TH1D* histNeg;
TH1D* histPos;
Char_t Title;
Char_t GraphPDF;
Char_t GraphPNG;
