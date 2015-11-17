// GoAT Physics to test Polarimeter acceptance from MC simulation
// This version does not include TAPS!

#include "PPiPolTest.h"

PPiPolTest::~PPiPolTest()
{
}

Bool_t	PPiPolTest::Init()
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

Bool_t	PPiPolTest::Start()
{
  if(!IsGoATFile())
    {
      cout << "ERROR: Input File is not a GoAT file." << endl;
      return kFALSE;
    }

  MCDataCheck();

  if (MCData == kFALSE){
    cout << "ERROR: Input file is not an MC File" << endl;
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

  cout << i << "  Events with forward neutron, proton into CB" << endl;
  cout << k << "  Events with forward neutron, proton into CB that is detected" << endl;

  return kTRUE;
}

void	PPiPolTest::ProcessEvent()
{

  GetEvent(); // Function gets number of tracks/protons/pions e.t.c.

  for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++){

    TaggerTime = GetTagger()->GetTaggedTime(j); // Get tagged time for event
    EGamma = (GetTagger()->GetTaggedEnergy(j)); // Get Photon energy for event

    MCTrueID(); // Get True ID for each particle
    ParticleAssignment(); // Assign Boolean flags to each
    GetTrueTheta(); // Get true theta for the neutron in each event
    if ((ThetaNTrue == 1000 ) || (ThetaPiTrue == 1000)) return; // Theta is set to 1000 if the event does not feature a proton or a neutron, shouldn't occur
    if (ThetaNTrue > 20) return; // We only want to select events where the neutron not going into CB
    if ( (ThetaPiTrue < 20) || (ThetaPiTrue > 160)) return; // Want proton to be going into CB
    i++; // Add one to our counter, this is the number of events where the neutron went forward and the proton went into CB

    TrueThetaPi->Fill(ThetaPiTrue, TaggerTime);
    TrueThetaN->Fill(ThetaNTrue, TaggerTime);
    TrueThetaPiThetaN->Fill(ThetaPiTrue, ThetaNTrue, TaggerTime);
    ThetaPiEPiTrue->Fill(ThetaPiTrue, EPiTrue, TaggerTime);
    EPiEgTrue->Fill(EPiTrue, EgTrue, TaggerTime);
    //TrueEgThetapEp->Fill(EgTrue, ThetapTrue, EpTrue, TaggerTime); // Cause a seg fault currently, too big?
    //TrueEgThetanEn->Fill(EgTrue, ThetanTrue, EnTrue, TaggerTime);

    // Next thing we want to check is how many events out of those we want actually have a proton being DETECTED in the CB
    // If our proton is going into CB expect it to produce a track
    // Do NOT expect neutron to give a track is it is going into TAPS which is switched off in the simulation currently

    if (NTrack != 1) return; // In the case we want there should only be one track, is this true here?
    // Need seperate case for Pi0's as would expect two gamma tracks which we reconstruct from

    GetDetectorHitInfo(); // Get # of hits in various detectors, expect a hit in all for a proton
    if ((CBHits == 0 || (PIDHits == 0)) || (MWPCHits == 0)) return; // If any of the detectors we expect a hit in does not detect a hit drop out
    GetTrackVector();
    GetTrackInfo(); // Get info for the track
    MCSmearing(); // Smear dE values
    k++; // At this point add to our second counter, these are the events where the proton has been detected
    // Do we want to check this track looks like a proton here first or not? (EdE)

    // At this stage we should probably add our values to plots just to see what energy protons we're left with for a given photon energy
    // need to background subtract at this stage though as we need to use the tagger info
    // We probably want to plot both the measured energy and the actual MC energy the particle had
    // Note we may need to add a tagger loop earlier to look at the energy each proton had for all the forward neutrons
    // Not just the ones we end up detecting

    // Time = ( GetTagger()->GetTaggedTime(j) - GetTracks()->GetTime(i) ); // maybe move this to AFTER the cuts once the Eg-EpSum loop has been checked?

    FillHists(); // Fill histograms with data generated

   }

    //}

  // cout << endl; // Use to separate out events

  //}
}

void	PPiPolTest::ProcessScalerRead()
{
	// Fill Tagger Scalers // Currently this seems to fill the file with loads of "TaggerAccScal" histograms
	//FillScalers(GetTC_scaler_min(),GetTC_scaler_max(),TaggerAccScal); // Don't know if these are needed so cut out for now
}

TCutG*	PPiPolTest::OpenCutFile(Char_t* filename, Char_t* cutname)
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

Bool_t PPiPolTest::MCDataCheck(){

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

Int_t PPiPolTest::GetEvent() // Gets basic info on particles for event
{
  NTrack = GetTracks()->GetNTracks();
  NP = GetProtons()->GetNParticles();
  NPi = GetChargedPions()->GetNParticles();
  NRoo = GetRootinos()->GetNParticles();
  NTag = GetTagger()->GetNTagged();

  return NTrack, NP, NPi, NRoo, NTag;
}

Int_t PPiPolTest::MCTrueID()
{

    MCTrueID1 = GetGeant()->GetTrueID(0);
    MCTrueID2 = GetGeant()->GetTrueID(1);

    return MCTrueID1, MCTrueID2;

}

Bool_t PPiPolTest::ParticleAssignment(){

    if (MCTrueID1 == 13 || MCTrueID1 == 14){
        Nucleon1 = kTRUE;
        Nucleon2 = kFALSE;
    }

    if (MCTrueID2 == 13 || MCTrueID2 == 14){
        Nucleon1 = kFALSE;
        Nucleon2 = kTRUE;
    }

    if (MCTrueID1 == 7 || MCTrueID1 == 8 || MCTrueID1 == 9){
        Pion1 = kTRUE;
        Pion2 = kFALSE;
        if(MCTrueID1 == 7) NeutralPion = kTRUE;
    }

    if (MCTrueID2 == 7 || MCTrueID2 == 8 || MCTrueID2 == 9){
        Pion1 = kFALSE;
        Pion2 = kTRUE;
        if(MCTrueID1 == 7) NeutralPion = kTRUE;
    }

    else{
        Nucleon1 = kFALSE;
        Nucleon2 = kFALSE;
        Pion1 = kFALSE;
        Pion2 = kFALSE;
        NeutralPion = kFALSE;
    }

    return Nucleon1, Nucleon2, Pion1, Pion2, NeutralPion;

}

Double_t PPiPolTest::GetTrueTheta(){

    if (Nucleon1 == kTRUE){
    NucleonTrueVect = GetGeant()->GetTrueVector(0);
    NucleonTrue3Vect = NucleonTrueVect.Vect();
    ThetaNTrue = NucleonTrue3Vect.Theta() * TMath::RadToDeg();
    ENTrue = (NucleonTrueVect(3)*1000) - Mn;
    TrueGamma = GetGeant()->GetBeam();
    EgTrue = TrueGamma(3)*1000;
    }

    if (Nucleon2 == kTRUE){
    NucleonTrueVect = GetGeant()->GetTrueVector(1);
    NucleonTrue3Vect = NucleonTrueVect.Vect();
    ThetaNTrue = NucleonTrue3Vect.Theta() * TMath::RadToDeg();
    ENTrue = (NucleonTrueVect(3)*1000) - Mn;
    TrueGamma = GetGeant()->GetBeam();
    EgTrue = TrueGamma(3)*1000;
    }

    else if ((Nucleon1 == kFALSE) && (Nucleon2 == kFALSE)) ThetaNTrue = 1000;

    if (Pion1 == kTRUE){
    PionTrueVect = GetGeant()->GetTrueVector(0);
    PionTrue3Vect = PionTrueVect.Vect();
    ThetaPiTrue = PionTrue3Vect.Theta() * TMath::RadToDeg();
    EPiTrue = (PionTrueVect(3)*1000)- Mp; // In this case the value in brackets is the component of the vector!
    TrueGamma = GetGeant()->GetBeam();
    EgTrue = TrueGamma(3)*1000;
    }

    if (Pion2 == kTRUE){
    PionTrueVect = GetGeant()->GetTrueVector(1);
    PionTrue3Vect = PionTrueVect.Vect();
    ThetaPiTrue = PionTrue3Vect.Theta() * TMath::RadToDeg();
    EPiTrue = (PionTrueVect(3)*1000)- Mp;
    TrueGamma = GetGeant()->GetBeam();
    EgTrue = TrueGamma(3)*1000;
    }

    else if ((Pion1 == kFALSE) && (Pion2 == kFALSE)) ThetaPiTrue = 1000, EPiTrue = 0, ENTrue = 0;

    return ThetaNTrue, ThetaPiTrue, EPiTrue, ENTrue, EgTrue;

}

Int_t PPiPolTest::GetDetectorHitInfo(){

    PIDHits = GetDetectorHits()->GetPIDHits(0);
    MWPCHits = GetDetectorHits()->GetMWPCHits(0);
    CBHits = GetDetectorHits()->GetNaIHits(0);

    return PIDHits, MWPCHits, CBHits;

}

TLorentzVector PPiPolTest::GetTrackVector(){

  GV1 = GetTracks()->GetVector(0, Mp);

  return GV1;

}

Double_t PPiPolTest::GetTrackInfo(){

  Theta1 = (GV1.Theta()) * TMath::RadToDeg();
  Phi1 = (GV1.Phi()) * TMath::RadToDeg();
  z1 = GetTracks()->GetPseudoVertexZ(0);
  E1 = GetTracks()->GetClusterEnergy(0);

  return Theta1, Phi1, z1, E1; // Returns various quantities used in later functions

}

PPiPolTest::PPiPolTest() // Define a load of histograms to fill
{

  time = new GH1("time", 	"time", 	1400, -700, 700);
  time_cut = new GH1("time_cut", 	"time_cut", 	1400, -700, 700);

  TrueThetaPi = new GH1 ("TrueThetaPi", "True Theta of Pion", 90, 0, 180);
  TrueThetaN = new GH1 ("TrueThetaN", "True Theta of Nucleon", 90, 0, 180);

  TrueThetaPiThetaN = new GH2("TrueThetaPiThetaN", "True Theta for Pion vs True Theta for Nucleon", 90, 0, 180, 90, 0, 180);
  ThetaPiEPiTrue = new GH2("ThetaPiEPiTrue", "True Theta vs Energy for Pions", 90, 0, 180, 200, 0, 400);
  ThetaPiEPi = new GH2("ThetaPiEPi", "Theta vs Energy for Pions", 90, 0, 180, 200, 0, 800);
  EPiEgTrue = new GH2("EPiEgTrue", "True Pion Energy vs True Photon Energy", 200, 0, 400, 400, 0, 1600);
  EPiEg = new GH2("EPiEg", "Pion Energy vs Photon Energy", 200, 0, 400, 400, 0, 1600);

  //TrueEgThetapEp = new GH3 ("TrueEgThetapEp", "Proton Energy as Fn of Eg and Theta (all TRUE values)", 400, 0, 1600, 60, 0, 180, 200, 0, 800 );
  //TrueEgThetapEp = new GH3 ("TrueEgThetapEn", "Neutron Energy as Fn of Eg and Theta (all TRUE values)", 400, 0, 1600, 60, 0, 180, 200, 0, 800 );
  EgThetaEPi = new GH3("EgThetaEPi", "Pion Energy as a Fn of Eg and Theta", 400, 0, 1600, 60, 0, 180, 200, 0, 800 );
  EgThetaEPiTrue = new GH3("EgThetaEPiTrue", "True Pion Energy as a Fn of Eg and True Theta", 400, 0, 1600, 60, 0, 180, 200, 0, 400 );

}

Double_t PPiPolTest::MCSmearing() // Smear dE values for MC data to represent Energy resolution of PID
{
  dE1 = rGen.Gaus(GetTracks()->GetVetoEnergy(0) , (0.29*(sqrt(GetTracks()->GetVetoEnergy(0)))));

  if (dE1 < 0) dE1 = 0.01;

  return dE1;
}

void PPiPolTest::FillHists()
{
    ThetaPiEPi->Fill(Theta1, E1, TaggerTime);
    EgThetaEPi->Fill(EGamma, Theta1, E1, TaggerTime);
    EgThetaEPiTrue->Fill(EGamma, ThetaPiTrue, EPiTrue, TaggerTime);
    EPiEg->Fill(E1, EGamma, TaggerTime);
}


Bool_t	PPiPolTest::Write(){
  // Write some TH1s - currently none to write so commented out
  // GTreeManager::Write(TaggerAccScal); // This spams the file with 2500+ histograms of "TaggerAccScal" so commented out

  // Write all GH1's easily

  GTreeManager::Write();
}
