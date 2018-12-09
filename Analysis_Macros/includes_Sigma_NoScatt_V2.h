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

    Double_t pValues415[20], pValues425[20], pValues435[20], pValues445[20], pValues455[20], pValues465[20], pValues475[20], pValues485[20], pValues495[20], pValues505[20], pValues515[20], pValues525[20], pValues535[20], pValues545[20], pValues555[20], pValues565[20], pValues575[20], pValues585[20], pValues595[20], pValues605[20], pValues615[20];
    Double_t pErrValues415[20], pErrValues425[20], pErrValues435[20], pErrValues445[20], pErrValues455[20], pErrValues465[20], pErrValues475[20], pErrValues485[20], pErrValues495[20], pErrValues505[20], pErrValues515[20], pErrValues525[20], pErrValues535[20], pErrValues545[20], pErrValues555[20], pErrValues565[20], pErrValues575[20], pErrValues585[20], pErrValues595[20], pErrValues605[20], pErrValues615[20];
    Double_t LegPar[8][21];
    Double_t LegParErr[8][21];
    Double_t SigmaValues[21][18];
    Double_t SigmaErrValues[21][18];
    double p0, p0Err, p1, p1Err, p2, p2Err, p3, p3Err, p4, p4Err, p5, p5Err, p6, p6Err, p7, p7Err;

    Double_t pCosAmp415, pCosAmpErr415;
    Double_t pCosAmp425, pCosAmpErr425;
    Double_t pCosAmp435, pCosAmpErr435;
    Double_t pCosAmp445, pCosAmpErr445;
    Double_t pCosAmp455, pCosAmpErr455;
    Double_t pCosAmp465, pCosAmpErr465;
    Double_t pCosAmp475, pCosAmpErr475;
    Double_t pCosAmp485, pCosAmpErr485;
    Double_t pCosAmp495, pCosAmpErr495;
    Double_t pCosAmp505, pCosAmpErr505;
    Double_t pCosAmp515, pCosAmpErr515;
    Double_t pCosAmp525, pCosAmpErr525;
    Double_t pCosAmp535, pCosAmpErr535;
    Double_t pCosAmp545, pCosAmpErr545;
    Double_t pCosAmp555, pCosAmpErr555;
    Double_t pCosAmp565, pCosAmpErr565;
    Double_t pCosAmp575, pCosAmpErr575;
    Double_t pCosAmp585, pCosAmpErr585;
    Double_t pCosAmp595, pCosAmpErr595;
    Double_t pCosAmp605, pCosAmpErr605;
    Double_t pCosAmp615, pCosAmpErr615;
