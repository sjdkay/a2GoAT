// GoAT Physics analysis to determine spin polarisation of protonsn
// from deuterium photodisintegration

#include "PProtPol.h"

PProtPol::PProtPol()
{
  time = new GH1("time", 	"time", 	1400, -700, 700);
  time_cut = new GH1("time_cut", 	"time_cut", 	1400, -700, 700);
  
}

PProtPol::~PProtPol()
{
}

Bool_t	PProtPol::Init()
{
  cout << "Initialising physics analysis..." << endl;
  cout << "--------------------------------------------------" << endl << endl;
  
  if(!InitBackgroundCuts()) return kFALSE;
  if(!InitTargetMass()) return kFALSE;
  if(!InitTaggerChannelCuts()) return kFALSE;
  if(!InitTaggerScalers()) return kFALSE;
  cout << "--------------------------------------------------" << endl;
  return kTRUE;
}

Bool_t	PProtPol::Start()
{
  if(!IsGoATFile())
    {
      cout << "ERROR: Input File is not a GoAT file." << endl;
      return kFALSE;
    }
  
  SetAsPhysicsFile();
  
  NP = 0;
  NPi = 0;
  NRoo = 0;
  Deut = TLorentzVector (0., 0., 0., 1875.613); // 4-Vector of Deuterium target
  
  Cut_CB_proton = OpenCutFile("configfiles/cuts/CB_DeltaE-E_Proton_27_02_15.root", "Proton");
  Cut_proton = Cut_CB_proton;
  Cut_CB_pion = OpenCutFile("configfiles/cuts/CB_DeltaE-E_Pion_27_02_15.root", "Pion");
  Cut_pion = Cut_CB_pion;
  cout << endl;
  
  TraverseValidEvents(); // This loops over each event as in old file and calls ProcessEvent() each loop
  
  return kTRUE;
}

void	PProtPol::ProcessEvent()
{
  
  GetEvent();
  if (NRoo != 0) return;
  InitialVect();
  InitialProp();
  
  for (Int_t i=0; i < NTrack; i++){ // Currently nothing relies upon i!

    for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++){
            
    }
  }
}

Int_t PProtPol::GetEvent(){

  NTrack = GetTracks()->GetNTracks();
  NP = GetProtons()-> GetNParticles();
  NPi = GetChargedPions()->GetNParticles();
  NRoo = GetRootinos()->GetNParticles();
  NTag = GetTagger()->GetNTagged();
  return NTrack, NP, NPi, NRoo, NTag;

}

TLorentzVector PProtPol::InitialVect()
{
}

Double_t PProtPol::InitialProp()
{
}

void	PProtPol::ProcessScalerRead()
{
	// Fill Tagger Scalers // Currently this seems to fill the file with loads of "TaggerAccScal" histograms
	//FillScalers(GetTC_scaler_min(),GetTC_scaler_max(),TaggerAccScal); // Don't know if these are needed so cut out for now
}

Bool_t	PProtPol::Write()
{
  // Write some TH1s - currently none to write so commented out
  // GTreeManager::Write(TaggerAccScal); // This spams the file with 2500+ histograms of "TaggerAccScal" so commented out
  
  // Write all GH1's easily
  GTreeManager::Write();
}

TCutG*	PProtPol::OpenCutFile(Char_t* filename, Char_t* cutname)
{
  CutFile = new TFile(filename, "READ");
  
  if( !CutFile || !CutFile->IsOpen() ) {
    cerr << "Can't open cut file: " << filename << endl;
    throw false;
  }
  
  // Try to find a TCutG with the name we want
  // GetObject checks the type to be TCutG,
    // see http://root.cern.ch/root/html534/TDirectory.html#TDirectory:GetObject
  CutFile->GetObject(cutname, Cut);
  
  if( !Cut ) {
    cerr << "Could not find a TCutG with the name " << cutname << " in " << filename << endl;
    throw false;
  }
  
  TCutG* Cut_clone = Cut;
  CutFile->Close();
  cout << "cut file " << filename <<
    " opened (Cut-name = " << cutname << ")"<< endl;
  return Cut_clone;
  
}
