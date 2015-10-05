// GoAT Physics analysis to identify neutrons from deuterium photodisintegration
// Various properties of neutrons/protons identified plotted in histograms
// Main aim is to determine spin polarisation of neutrons

#include "P2pi0.h"

P2pi0::~P2pi0()
{
}

Bool_t	P2pi0::Init()
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

Bool_t	P2pi0::Start()
{
  if(!IsGoATFile())
    {
      cout << "ERROR: Input File is not a GoAT file." << endl;
      return kFALSE;
    }

  SetAsPhysicsFile();

  i = 0; // Integer counter
  k = 0; // Another integer counter
  Deut = TLorentzVector (0., 0., 0., 1875.613); // 4-Vector of Deuterium target, assume at rest
  Mn = 939.565; // Mass of neutron
  Mp = 938.272; // Mass of proton
  Mpi = 139.57; // Mass of charged pions

  cout << endl;

  TraverseValidEvents(); // This loops over each event as in old file and calls ProcessEvent() each loop

  cout << " Detected " << (i + (k*2)) << " Double Gamma Scatter events out of 10^7 Pi0's" << endl;
  cout << i << " Events where only one pair detected" << endl;
  cout << k << " Events where two pairs were detected" << endl;

}

void	P2pi0::ProcessEvent()
{
  GetEvent(); // Function gets number of tracks/protons/pions e.t.c.
  if (NTrack !=4) return; // Ensures two track event
  if (NPhotons != 0) { // If we don't have 0 photons, want one OR two (one of the others could scatter still)
    if(NPhotons != 2 || NPhotons != 1) return;
  }
  InitialVect(); // Function gets vectors of identified tracks and returns them
  InitialProp(); // Function gets initial properties (energy, vertex e.t.c.) of identified tracks
  //MCTrue(); // Get MC True data and evaluate it
  //for (Int_t i=0; i < NTrack; i++){ // Currently nothing relies upon i!

  for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++){

    // Time = ( GetTagger()->GetTaggedTime(j) - GetTracks()->GetTime(i) ); // maybe move this to AFTER the cuts once the Eg-EpSum loop has been checked?
    TaggerTime = GetTagger()->GetTaggedTime(j); // Get tagged time for event
    EGamma = (GetTagger()->GetTaggedEnergy(j)); // Get Photon energy for event

    // If neither particle pair has both tracks giving no hits in the PID continue
    if ((PIDHits1 != 0 || PIDHits2 !=0) && (PIDHits3 != 0 || PIDHits4 !=0)) continue;

    // If first pair has tracks giving no tracks in PID but other pair does, check MWPC for first pair
    if ((PIDHits1 == 0 && PIDHits2 == 0) && (PIDHits3 != 0 || PIDHits4 !=0)){

        if (MWPCHits_1 == 0 || MWPCHits_2 == 0) continue; // Check get hits in MWPC for first pair
        i++; // If hits in MWPC for both add to counter i

    }

    // If second pair has tracks giving no tracks in PID but other pair does, check MWPC for second pair
    if ((PIDHits1 != 0 || PIDHits2 !=0) && (PIDHits3 == 0 && PIDHits4 == 0)){

        if (MWPCHits_3 == 0 || MWPCHits_4 == 0) continue; // Check get hits in MWPC for second pair
        i++; // If hits in MWPC for both add to counter i

    }

    // If both particle pairs have no hits in PID, check MWPC for both
    if ((PIDHits1 == 0 && PIDHits2 == 0) && (PIDHits3 == 0 && PIDHits4 == 0)){

        if (MWPCHits_1 == 0 || MWPCHits_2 == 0) {
            if (MWPCHits_3 == 0 || MWPCHits_4 == 0) continue;
            i++;
        }

        if (MWPCHits_3 == 0 || MWPCHits_4 == 0) {
            if (MWPCHits_1 == 0 || MWPCHits_2 == 0) continue;
            i++;
        }

        if ((MWPCHits_1 != 0 && MWPCHits_2 != 0) && (MWPCHits_3 != 0 && MWPCHits_4 != 0)) k++;

    }

    FillHists(); // Fill histograms with data generated

    //}
  }

  // cout << endl; // Use to separate out events

  //}
}

void	P2pi0::ProcessScalerRead()
{
	// Fill Tagger Scalers // Currently this seems to fill the file with loads of "TaggerAccScal" histograms
	//FillScalers(GetTC_scaler_min(),GetTC_scaler_max(),TaggerAccScal); // Don't know if these are needed so cut out for now
}


Int_t P2pi0::GetEvent() // Gets basic info on particles for event
{
  NTrack = GetTracks()->GetNTracks();
  NTag = GetTagger()->GetNTagged();
  NPhotons = GetPhotons()->GetNParticles();
  return NTrack, NTag, NPhotons;
}

TLorentzVector P2pi0::InitialVect() // Defines initial vectors
{
  GV1 = GetTracks()->GetVector(0, Mpi);
  GV2 = GetTracks()->GetVector(1, Mpi);// Set both to have pion mass for now
  GV3 = GetTracks()->GetVector(2, Mpi);
  GV4 = GetTracks()->GetVector(3, Mpi);
  return GV1, GV2, GV3, GV4; // Returns intial 4-vectors for use in later functions
}

Double_t P2pi0::InitialProp() // Defines initial particle properties
{

  MWPCHits_1 = GetDetectorHits()->GetMWPCHits(0);
  MWPCHits_2 = GetDetectorHits()->GetMWPCHits(1);
  MWPCHits_3 = GetDetectorHits()->GetMWPCHits(2);
  MWPCHits_4 = GetDetectorHits()->GetMWPCHits(3);
  PIDHits1 = GetDetectorHits()->GetPIDHits(0);
  PIDHits2 = GetDetectorHits()->GetPIDHits(1);
  PIDHits3 = GetDetectorHits()->GetPIDHits(2);
  PIDHits4 = GetDetectorHits()->GetPIDHits(3);

  return MWPCHits_1, MWPCHits_2, MWPCHits_3, MWPCHits_4, PIDHits1, PIDHits2, PIDHits3, PIDHits4; // Returns various quantities used in later functions
}

Double_t P2pi0::MCTrue()
{

    //Unsure how to actually extract MC true data from its relevant tree?
    Double_t Example = 5;
    // Get MC True info and do some stuff
    return Example;

}

P2pi0::P2pi0() // Define a load of histograms to fill
{

  Ek = new GH1( "Ek", "Particle Energy Distribution", 100, 0, 1600 );

}

void P2pi0::FillHists() // Fill histograms with stuff we've calculated
{

}

Bool_t	P2pi0::Write(){
  // Write some TH1s - currently none to write so commented out
  // GTreeManager::Write(TaggerAccScal); // This spams the file with 2500+ histograms of "TaggerAccScal" so commented out

  // Write all GH1's easily

  GTreeManager::Write();

}
