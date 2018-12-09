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

    Double_t pValues410[5], pValues430[5], pValues450[5], pValues470[5], pValues490[5], pValues510[5], pValues530[5], pValues550[5], pValues570[5], pValues590[5];
    Double_t pErrValues410[5], pErrValues430[5], pErrValues450[5], pErrValues470[5], pErrValues490[5], pErrValues510[5], pErrValues530[5], pErrValues550[5], pErrValues570[5], pErrValues590[5];
    Double_t LegPar[8][10];
    Double_t LegParErr[8][10];
    Double_t SigmaValues[10][5];
    Double_t SigmaErrValues[10][5];
    double p0, p0Err, p1, p1Err, p2, p2Err, p3, p3Err, p4, p4Err, p5, p5Err, p6, p6Err, p7, p7Err;

    Double_t pCosAmp410, pCosAmpErr410;
    Double_t pCosAmp430, pCosAmpErr430;
    Double_t pCosAmp450, pCosAmpErr450;
    Double_t pCosAmp470, pCosAmpErr470;
    Double_t pCosAmp490, pCosAmpErr490;
    Double_t pCosAmp510, pCosAmpErr510;
    Double_t pCosAmp530, pCosAmpErr530;
    Double_t pCosAmp550, pCosAmpErr550;
    Double_t pCosAmp570, pCosAmpErr570;
    Double_t pCosAmp590, pCosAmpErr590;
