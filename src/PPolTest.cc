// GoAT Physics to test Polarimeter acceptance from MC simulation
// This version does not include TAPS!

#include "PPolTest.h"

PPolTest::~PPolTest()
{
}

Bool_t	PPolTest::Init()
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

Bool_t	PPolTest::Start()
{
  if(!IsGoATFile())
    {
      cout << "ERROR: Input File is not a GoAT file." << endl;
      return kFALSE;
    }

  MCDataCheck();

  if (MCData == kFALSE){
    cout << "Error : Input file is not an MC File" << endl;
    cout << "This code is only designed to work with an MC input" << endl;
    return kFALSE;
  }

  SetAsPhysicsFile();

  i = 0; // Integer counter
  k = 0;
  d = 54.2; // Distance from centre of target to centre of PID
  NP = 0; // Set number of Protons to 0 before checking
  NPi = 0; // Set number of pions to 0 before checking
  NRoo = 0; // Set number of Rootinos to 0 before checking
  Deut = TLorentzVector (0., 0., 0., 1875.613); // 4-Vector of Deuterium target, assume at rest
  Mn = 939.565; // Mass of neutron
  Mp = 938.272; // Mass of proton

  Cut_CB_proton = OpenCutFile("configfiles/cuts/CB_DeltaE-E_Proton_29_07_15.root", "Proton"); // These will need adjusting with new Acqu files
  Cut_proton = Cut_CB_proton;
  Cut_CB_pion = OpenCutFile("configfiles/cuts/CB_DeltaE-E_Pion_29_07_15.root", "Pion");
  Cut_pion = Cut_CB_pion;
  Cut_CB_ROI = OpenCutFile("configfiles/cuts/CB_DeltaE-E_ROI_05_08_15.root", "ROI");
  Cut_ROI = Cut_CB_ROI;
  // Cut_CB_neutron = OpenCutFile("configfiles/cuts/CB_DeltaE-E_Neutron_25_02_15.root", "Neutron");
  // Cut_neutron = Cut_CB_neutron; // This only needs to be here if we get simulation to show where we expect neutrons
  cout << endl;

  MCDataCheck();
  if (MCData == kFALSE){
    cout << "Error : Input file is not an MC File" << endl;
    cout << "This code is only designed to work with an MC input" << endl;
    return kFALSE;
  }

  TraverseValidEvents(); // This loops over each event as in old file and calls ProcessEvent() each loop

  cout << k << endl;

  return kTRUE;
}

void	PPolTest::ProcessEvent()
{

  GetEvent(); // Function gets number of tracks/protons/pions e.t.c.


  MCTrueID(); // Get True ID for each particle
  ParticleAssignment(); // Assign Boolean flags to each
  GetTrueTheta(); // Get true theta for the neutron in each event
  if ((ThetanTrue == 1000 ) || (ThetapTrue == 1000)) return; // Theta is set to 1000 if the event does not feature a proton or a neutron, shouldn't occur
  if (ThetanTrue > 20) return; // We only want to select events where the neutron not going into CB
  if ( (ThetapTrue < 20) || (ThetapTrue > 160)) return; // Want proton to be going into CB
  i++; // Add one to our counter, this is the number of events where the neutron went forward and the proton went into CB

  // Next thing we want to check is how many events out of those we want actually have a proton being DETECTED in the CB
  // If our proton is going into CB expect it to produce a track
  // Do NOT expect neutron to give a track is it is going into TAPS which is switched off in the simulation currently

  GetDetectorHits(); // Get # of hits in various detectors, expect a hit in all for a proton


   for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++){

   }

  FillHists(); // Fill histograms with data generated

    //}

  // cout << endl; // Use to separate out events

  //}
}

void	PPolTest::ProcessScalerRead()
{
	// Fill Tagger Scalers // Currently this seems to fill the file with loads of "TaggerAccScal" histograms
	//FillScalers(GetTC_scaler_min(),GetTC_scaler_max(),TaggerAccScal); // Don't know if these are needed so cut out for now
}

TCutG*	PPolTest::OpenCutFile(Char_t* filename, Char_t* cutname)
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
  cout << "cut file " << filename << " opened (Cut-name = " << cutname << ")"<< endl;
  return Cut_clone;
}

Bool_t PPolTest::MCDataCheck(){

  if (GetScalers()->GetNEntries() == 0) // MC Data has no scaler entries so if 0, data gets a flag to denote it as MC
  {
      MCData = kTRUE;
      TRandom2 *rGen = new TRandom2(0); // Define new random generator for use in MC smearing fn
  }

  else if (GetScalers()->GetNEntries() != 0) // This flag is listed as false if the number of scaler entries does not equal 0

  {
      MCData = kFALSE;
  }

  return MCData;
}

Int_t PPolTest::GetEvent() // Gets basic info on particles for event
{
  NTrack = GetTracks()->GetNTracks();
  NP = GetProtons()->GetNParticles();
  NPi = GetChargedPions()->GetNParticles();
  NRoo = GetRootinos()->GetNParticles();
  NTag = GetTagger()->GetNTagged();
  return NTrack, NP, NPi, NRoo, NTag;
}

Int_t PPolTest::MCTrueID()
{

    MCTrueID1 = GetGeant()->GetTrueID(0);
    MCTrueID2 = GetGeant()->GetTrueID(1);
    return MCTrueID1, MCTrueID2;

}

Bool_t PPolTest::ParticleAssignment(){

    if (MCTrueID1 == 14){
        Proton1 == kTRUE;
        Proton2 == kFALSE;
    }

    if (MCTrueID2 == 14){
        Proton1 == kFALSE;
        Proton2 == kTRUE;
    }


    if (MCTrueID1 == 13){
        Neutron1 == kTRUE;
        Neutron2 == kFALSE;
    }

    if (MCTrueID2 == 13){
        Neutron1 == kFALSE;
        Neutron2 == kTRUE;
    }

    else{
        Proton1 == kFALSE;
        Proton2 == kFALSE;
        Neutron1 == kFALSE;
        Neutron2 == kFALSE;
    }

}

Double_t PPolTest::GetTrueTheta(){

    if (Neutron1 == kTRUE){
    NeutronTrueVect = GetGeant()->GetTrueVector(0);
    NeutronTrue3Vect = NeutronTrueVect.Vect();
    ThetanTrue = NeutronTrue3Vect.Theta() * TMath::RadToDeg();
    }

    if (Neutron2 == kTRUE){
    NeutronTrueVect = GetGeant()->GetTrueVector(1);
    NeutronTrue3Vect = NeutronTrueVect.Vect();
    ThetanTrue = NeutronTrue3Vect.Theta() * TMath::RadToDeg();
    }

    else if ((Neutron1 == kFALSE) && (Neutron2 == kFALSE)) ThetanTrue = 1000;

    if (Proton1 == kTRUE){
    ProtonTrueVect = GetGeant()->GetTrueVector(0);
    ProtonTrue3Vect = ProtonTrueVect.Vect();
    ThetapTrue = ProtonTrue3Vect.Theta() * TMath::RadToDeg();
    }

    if (Proton2 == kTRUE){
    ProtonTrueVect = GetGeant()->GetTrueVector(1);
    ProtonTrue3Vect = ProtonTrueVect.Vect();
    ThetapTrue = ProtonTrue3Vect.Theta() * TMath::RadToDeg();
    }

    else if ((Proton1 == kFALSE) && (Proton2 == kFALSE)) ThetapTrue = 1000;

    return ThetanTrue, ThetapTrue;

}

PPolTest::PPolTest() // Define a load of histograms to fill
{

  time = new GH1("time", 	"time", 	1400, -700, 700);
  time_cut = new GH1("time_cut", 	"time_cut", 	1400, -700, 700);

}

void PPolTest::FillHists()
{
}


Bool_t	PPolTest::Write(){
  // Write some TH1s - currently none to write so commented out
  // GTreeManager::Write(TaggerAccScal); // This spams the file with 2500+ histograms of "TaggerAccScal" so commented out

  // Write all GH1's easily

  GTreeManager::Write();
}
