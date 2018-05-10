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

TH1D* Eg;
TH1D* Eg_Para;
TH1D* Eg_Perp;

TH1D* PhiScPosHel[7][5];
TH1D* PhiScNegHel[7][5];

TH1D* PhiScPosHelPara[7][5];
TH1D* PhiScNegHelPara[7][5];

TH1D* PhiScPosHelPerp[7][5];
TH1D* PhiScNegHelPerp[7][5];

Double_t NPara;
Double_t NPerp;
Double_t ScaleFactor;
Double_t ScaleFactorErr;

