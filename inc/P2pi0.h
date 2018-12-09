#ifndef __P2pi0_h__
#define __P2pi0_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
using namespace std;
#include "GTreeManager.h"
#include "PPhysics.h"
#include "TCutG.h"#ifndef __P2pi0_h__
#define __P2pi0_h__

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

class	P2pi0 : public PPhysics
{
private:

  TH1*	TaggerAccScal;

  Int_t NTag;
  Int_t NTrack;
  Int_t NPhotons;
  Int_t i;
  Int_t k;

  double_t Time;
  double_t TaggerTime;
  double_t EGamma;
  double_t Mn;
  double_t Mp;
  double_t Mpi;
  double_t MWPCHits_1;
  double_t MWPCHits_2;
  double_t MWPCHits_3;
  double_t MWPCHits_4;
  double_t PIDHits1;
  double_t PIDHits2;
  double_t PIDHits3;
  double_t PIDHits4;

  TLorentzVector GV1;
  TLorentzVector GV2;
  TLorentzVector GV3;
  TLorentzVector GV4;
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

    P2pi0();
    virtual ~P2pi0();
    virtual Bool_t  Init();
    Int_t GetEvent();
    TLorentzVector InitialVect();
    Double_t InitialProp();
    Double_t MCTrue();
    void FillHists();

};
#endif
