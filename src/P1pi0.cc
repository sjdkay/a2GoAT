// GoAT Physics analysis to identify neutrons from deuterium photodisintegration
// Various properties of neutrons/protons identified plotted in histograms
// Main aim is to determine spin polarisation of neutrons

#include "P1pi0.h"

P1pi0::~P1pi0()
{
}

Bool_t	P1pi0::Init()
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

Bool_t	P1pi0::Start()
{
  if(!IsGoATFile())
    {
      cout << "ERROR: Input File is not a GoAT file." << endl;
      return kFALSE;
    }

  SetAsPhysicsFile();

  i = 0; // Integer counter
  Deut = TLorentzVector (0., 0., 0., 1875.613); // 4-Vector of Deuterium target, assume at rest
  Mn = 939.565; // Mass of neutron
  Mp = 938.272; // Mass of proton
  Mpi = 139.57; // Mass of charged pions

  cout << endl;

  TraverseValidEvents(); // This loops over each event as in old file and calls ProcessEvent() each loop

  cout << " Detected " << i << " Double Gamma Scatter events out of 10^7 Pi0's" << endl;

}

void	P1pi0::ProcessEvent()
{
  GetEvent(); // Function gets number of tracks/protons/pions e.t.c.
  if (NTrack !=2) return; // Ensures two track event
  if (NPhotons != 0) return;
  InitialVect(); // Function gets vectors of identified tracks and returns them
  InitialProp(); // Function gets initial properties (energy, vertex e.t.c.) of identified tracks
  //MCTrue(); // Get MC True data and evaluate it
  //for (Int_t i=0; i < NTrack; i++){ // Currently nothing relies upon i!

  for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++){

    // Time = ( GetTagger()->GetTaggedTime(j) - GetTracks()->GetTime(i) ); // maybe move this to AFTER the cuts once the Eg-EpSum loop has been checked?
    TaggerTime = GetTagger()->GetTaggedTime(j); // Get tagged time for event
    EGamma = (GetTagger()->GetTaggedEnergy(j)); // Get Photon energy for event

    if (PIDHits1 != 0 || PIDHits2 !=0) continue; // If either track has hit in PID continue
    if (MWPCHits_1 == 0 || MWPCHits_2 == 0) continue; // If either track does NOT have hit in MWPC continue

    i++;

    FillHists(); // Fill histograms with data generated

    //}
  }

  // cout << endl; // Use to separate out events

  //}
}

void	P1pi0::ProcessScalerRead()
{
	// Fill Tagger Scalers // Currently this seems to fill the file with loads of "TaggerAccScal" histograms
	//FillScalers(GetTC_scaler_min(),GetTC_scaler_max(),TaggerAccScal); // Don't know if these are needed so cut out for now
}


Int_t P1pi0::GetEvent() // Gets basic info on particles for event
{
  NTrack = GetTracks()->GetNTracks();
  NTag = GetTagger()->GetNTagged();
  NPhotons = GetPhotons()->GetNParticles();
  return NTrack, NTag, NPhotons;
}

TLorentzVector P1pi0::InitialVect() // Defines initial vectors
{
  GV1 = GetTracks()->GetVector(0, Mpi);
  GV2 = GetTracks()->GetVector(1, Mpi); // Set both to have proton mass for now
  return GV1, GV2; // Returns intial 4-vectors for use in later functions
}

Double_t P1pi0::InitialProp() // Defines initial particle properties
{
  Theta1 = (GV1.Theta()) * TMath::RadToDeg();
  Theta2 = (GV2.Theta()) * TMath::RadToDeg();
  z1 = GetTracks()->GetPseudoVertexZ(0);
  z2 = GetTracks()->GetPseudoVertexZ(1);
  E1 = GetTracks()->GetClusterEnergy(0);
  E2 = GetTracks()->GetClusterEnergy(1);
  MWPCHits_1 = GetDetectorHits()->GetMWPCHits(0);
  MWPCHits_2 = GetDetectorHits()->GetMWPCHits(1);
  PIDHits1 = GetDetectorHits()->GetPIDHits(0);
  PIDHits2 = GetDetectorHits()->GetPIDHits(1);

  return Theta1, Theta2, z1, z2, E1, E2, MWPCHits_1, MWPCHits_2, PIDHits1, PIDHits2; // Returns various quantities used in later functions
}

Double_t P1pi0::MCTrue()
{

    //Unsure how to actually extract MC true data from its relevant tree?
    Double_t Example = 5;
    // Get MC True info and do some stuff
    return Example;

}

P1pi0::P1pi0() // Define a load of histograms to fill
{

  Ek = new GH1( "Ek", "Particle Energy Distribution", 100, 0, 1600 );

}

void P1pi0::FillHists() // Fill histograms with stuff we've calculated
{

  Ek->Fill (E1, TaggerTime);
  Ek->Fill (E2, TaggerTime);

}

Bool_t	P1pi0::Write(){
  // Write some TH1s - currently none to write so commented out
  // GTreeManager::Write(TaggerAccScal); // This spams the file with 2500+ histograms of "TaggerAccScal" so commented out

  // Write all GH1's easily

  GTreeManager::Write();

}
