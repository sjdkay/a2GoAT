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

    TH1D* Sig415Asymm[20];
    TH1D* Sig425Asymm[20];
    TH1D* Sig435Asymm[20];
    TH1D* Sig445Asymm[20];
    TH1D* Sig455Asymm[20];
    TH1D* Sig465Asymm[20];
    TH1D* Sig475Asymm[20];
    TH1D* Sig485Asymm[20];
    TH1D* Sig495Asymm[20];
    TH1D* Sig505Asymm[20];
    TH1D* Sig515Asymm[20];
    TH1D* Sig525Asymm[20];
    TH1D* Sig535Asymm[20];
    TH1D* Sig545Asymm[20];
    TH1D* Sig555Asymm[20];
    TH1D* Sig565Asymm[20];
    TH1D* Sig575Asymm[20];
    TH1D* Sig585Asymm[20];
    TH1D* Sig595Asymm[20];
    TH1D* Sig605Asymm[20];
    TH1D* Sig615Asymm[20];
    TH1D* Sig625Asymm[20];



