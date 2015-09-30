#ifndef __P1pi0_h__
#define __P1pi0_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
using namespace std;
#include "GTreeManager.h"
#include "PPhysics.h"
#include "TCutG.h"
#include "TObject.h"
#include "TGraph.h"
#include "TRandom2.h"
#include "TMath.h"

class	P1pi0 : public PPhysics
{
private:

  TH1*	TaggerAccScal;

  Int_t NTag;
  Int_t NTrack;
  Int_t i;

  double_t Time;
  double_t TaggerTime;
  double_t EGamma;
  double_t Mn;
  double_t Mp;
  double_t Mpi;
  double_t Theta1;
  double_t Theta2;
  double_t z1;
  double_t z2;
  double_t E1;
  double_t E2;
  double_t MWPCHits_1;
  double_t MWPCHits_2;
  double_t PIDHits1;
  double_t PIDHits2;

  TLorentzVector GV1;
  TLorentzVector GV2;
  TLorentzVector Deut;

  TRandom2 rGen;

  GH1* Ek;

  GH1*	time;
  GH1*	time_cut;

protected:
    virtual Bool_t  Start();
    virtual void    ProcessEvent();
    virtual void    ProcessScalerRead();
    virtual Bool_t  Write();

public:

    P1pi0();
    virtual ~P1pi0();
    virtual Bool_t  Init();
    Int_t GetEvent();
    TLorentzVector InitialVect();
    Double_t InitialProp();
    Double_t MCTrue();
    void FillHists();

};
#endif
