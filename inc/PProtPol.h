#ifndef __PProtPol_h__
#define __PProtPol_h__

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

class	PProtPol : public PPhysics
{
private:

  TH1*	TaggerAccScal;

  Int_t NP;
  Int_t NPi;
  Int_t NRoo;
  Int_t NTag;
  Int_t NTrack;

  double_t Time;
  double_t TaggerTime;
  TLorentzVector Gamma;
  TLorentzVector Deut;

  GH1*	time;
  GH1*	time_cut;

  char cutfilename[256];
  char cutname[256];
  TFile* CutFile;
  TCutG* Cut;
  TCutG* Cut_proton;
  TCutG* Cut_proton_4Sig;
  TCutG* Cut_proton_Full;
  TCutG* Cut_pion;
  TCutG* Cut_neutron;
  TCutG* Cut_CB_proton;
  TCutG* Cut_CB_proton_4Sig;
  TCutG* Cut_CB_proton_Full;
  TCutG* Cut_CB_pion;
  TCutG* Cut_CB_neutron;

protected:
    virtual Bool_t  Start();
    virtual void    ProcessEvent();
    virtual void    ProcessScalerRead();
    virtual Bool_t  Write();

public:

    PProtPol();
    virtual ~PProtPol();
    virtual Bool_t  Init();
    TCutG* OpenCutFile(Char_t* filename, Char_t* cutname);
    Int_t GetEvent();
    TLorentzVector InitialVect();
    Double_t InitialProp();
    Double_t PNProp(Int_t ProtonParticleNumber);
    TLorentzVector PNVect(Int_t ProtonParticleNumber);
    void LabBoost();
    void LabScatter();
    void nFrameScatter();
    void FillHists();

};
#endif
