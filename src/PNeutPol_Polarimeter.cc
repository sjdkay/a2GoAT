// GoAT Physics analysis to identify neutrons from deuterium photodisintegration
// Various properties of neutrons/protons identified plotted in histograms
// Main aim is to determine spin polarisation of neutrons

#include "PNeutPol_Polarimeter.h"

PNeutPol_Polarimeter::~PNeutPol_Polarimeter()
{
}

Bool_t	PNeutPol_Polarimeter::Init()
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

Bool_t	PNeutPol_Polarimeter::Start()
{
  if(!IsGoATFile())
    {
      cout << "ERROR: Input File is not a GoAT file." << endl;
      return kFALSE;
    }

  SetAsPhysicsFile();

  i = 0; // Integer counter
  k = 0;
  d = 52; // Distance from centre of target to centre of Polarimeter
  NP = 0; // Set number of Protons to 0 before checking
  NPi = 0; // Set number of pions to 0 before checking
  NRoo = 0; // Set number of Rootinos to 0 before checking
  Deut = TLorentzVector (0., 0., 0., 1875.613); // 4-Vector of Deuterium target, assume at rest
  Mn = 939.565; // Mass of neutron in MeV
  Mp = 938.272; // Mass of proton in MeV

  Cut_CB_proton = OpenCutFile("configfiles/cuts/CB_DeltaE-E_Proton_29_07_15.root", "Proton"); // These will need adjusting with new Acqu files
  Cut_proton = Cut_CB_proton;
  Cut_CB_pion = OpenCutFile("configfiles/cuts/CB_DeltaE-E_Pion_29_07_15.root", "Pion");
  Cut_pion = Cut_CB_pion;
  //Cut_CB_ROI = OpenCutFile("configfiles/cuts/CB_DeltaE-E_ROI_05_08_15.root", "ROI");
  //Cut_ROI = Cut_CB_ROI; // Uncomment these two if you want some ROI cut to be loaded and used
  //Cut_CB_neutron = OpenCutFile("configfiles/cuts/CB_DeltaE-E_Neutron_25_02_15.root", "Neutron");
  //Cut_neutron = Cut_CB_neutron; // This only needs to be here if we get simulation to show where we expect neutrons
  cout << endl;

  MCDataCheck();

  if (MCData == kTRUE) MCHists();

  TraverseValidEvents(); // This loops over each event as in old file and calls ProcessEvent() each loop

  //cout << k << endl;

  return kTRUE;
}

void	PNeutPol_Polarimeter::ProcessEvent()
{

  GetEvent(); // Function gets number of tracks/protons/pions e.t.c.
  if (NRoo !=0) return; // Goes to next event if any "rootinos" found
  if (NTrack !=2) return; // Ensures two track event
  DetectorCheck(); // Function checks detector numbers for each track

  // Currently try to accept only MWPC+PID+NaI and NaI + MWPC
  // If track 1 only gives signals in MWPC it is the neutron
  if((Detectors1 == 7) && (Detectors2 == 5))
  {
    Proton1 = kTRUE;
    Proton2 = kFALSE;
  }

  // If track 2 only gives signals in MWPC it is the neutron
  else if((Detectors1 == 5) && (Detectors2 == 7))
  {
    Proton1 = kFALSE;
    Proton2 = kTRUE;
  }

  // Drop out on ANY other condition (for now)
  else
  {
    return;
  }

  InitialVect(); // Function gets vectors of identified tracks and returns them
  InitialProp(); // Function gets initial properties (energy, vertex e.t.c.) of identified tracks

  if ( MCData == kTRUE)
  {
      MCSmearing(); // Smear dE values for MC data
      MCTrueID();
      MCTrueVectors();
      MCTheta1True = (MCTrueVect1.Theta()) * TMath::RadToDeg();
      MCTheta2True = (MCTrueVect2.Theta()) * TMath::RadToDeg();
      MCE1True = (MCTrueVect1.Energy()*1000)-Mp; // *1000 to get MeV
      MCE2True = (MCTrueVect2.Energy()*1000)-Mn;
  }

  //for (Int_t i=0; i < NTrack; i++){ // Currently nothing relies upon i!

  for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++){

    //cout << DetectorsSum << endl;
    // Time = ( GetTagger()->GetTaggedTime(j) - GetTracks()->GetTime(i) ); // maybe move this to AFTER the cuts once the Eg-EpSum loop has been checked?
    TaggerTime = GetTagger()->GetTaggedTime(j); // Get tagged time for event
    EGamma = (GetTagger()->GetTaggedEnergy(j)); // Get Photon energy for event

    if (Proton1 == kTRUE)
    {
        PNProp(1);
        PNVect(1);
    }

    else if (Proton2 == kTRUE)
    {
        PNProp(2);
        PNVect(2);
    }

    else
    {
        continue;
    }

    //if (i == 0) { // These properties get defined for each photon

    GVp3 = GVp.Vect(); // Generate some 3-vectors from the 4-vectors we have
    GVn3 = GVn.Vect();
    WC3Vectors();
    WCAngles();
    Gamma = TLorentzVector (0., 0., EGamma , EGamma); // 4-Vector of Photon beam
    Gamma3 = Gamma.Vect(); // Convert photon beam 4-vector to 3-vector
    B = (Deut + Gamma).Beta(); // Calculate Beta
    b = TVector3(0., 0., B); // Define boost vector
    LabAngles(); // Get angles in lab based on track info

    // Cut on difference between Phip and PhinRec next - If Diff =/= 180 cut
    //PhiDiff = abs (Phip - PhinRec);
    //if ((PhiDiff < 165) || (PhiDiff > 195)) continue;

    //FillTime(*GetProtons(),time); Needs to be get tracks not protons now
    //FillTimeCut(*GetProtons(),time_cut);

    if (MCData == kTRUE)
    {
        // Get some MCTrue parameters and MAYBE fill some histograms with them
        MCTheta1 = (GetTracks()->GetVector(0, Mp)).Theta() * TMath::RadToDeg();
        MCTheta2 = (GetTracks()->GetVector(1, Mp)).Theta() * TMath::RadToDeg();
        MCE1 = GetTracks()->GetClusterEnergy(0);
        MCE2 = GetTracks()->GetClusterEnergy(1);
    }

    //k++;
    FillHists(); // Fill histograms with data generated

    //}
  }

  // cout << endl; // Use to separate out events

  //}
}

void	PNeutPol_Polarimeter::ProcessScalerRead()
{
	// Fill Tagger Scalers // Currently this seems to fill the file with loads of "TaggerAccScal" histograms
	//FillScalers(GetTC_scaler_min(),GetTC_scaler_max(),TaggerAccScal); // Don't know if these are needed so cut out for now
}

TCutG*	PNeutPol_Polarimeter::OpenCutFile(Char_t* filename, Char_t* cutname)
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

Bool_t PNeutPol_Polarimeter::MCDataCheck(){

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

Int_t PNeutPol_Polarimeter::GetEvent() // Gets basic info on particles for event
{
  NTrack = GetTracks()->GetNTracks();
  NP = GetProtons()->GetNParticles();
  NPi = GetChargedPions()->GetNParticles();
  NRoo = GetRootinos()->GetNParticles();
  NTag = GetTagger()->GetNTagged();
  return NTrack, NP, NPi, NRoo, NTag;
}

TLorentzVector PNeutPol_Polarimeter::InitialVect() // Defines initial vectors
{
  GV1 = GetTracks()->GetVector(0, Mp);
  GV2 = GetTracks()->GetVector(1, Mp); // Set both to have proton mass for now
  return GV1, GV2; // Returns intial 4-vectors for use in later functions
}

Double_t PNeutPol_Polarimeter::InitialProp() // Defines initial particle properties
{
  Theta1 = (GV1.Theta()) * TMath::RadToDeg();
  Theta2 = (GV2.Theta()) * TMath::RadToDeg();
  Phi1 = (GV1.Phi()) * TMath::RadToDeg();
  Phi2 = (GV2.Phi()) * TMath::RadToDeg();
  z1 = GetTracks()->GetPseudoVertexZ(0);
  z2 = GetTracks()->GetPseudoVertexZ(1);
  E1 = GetTracks()->GetClusterEnergy(0);
  E2 = GetTracks()->GetClusterEnergy(1);
  dE1 = GetTracks()->GetVetoEnergy(0);
  dE2 = GetTracks()->GetVetoEnergy(1);
  WC1X1 = GetMWPCHitsChris()->GetMWPCChamber1X(0); // Don't need tree getters
  WC1Y1 = GetMWPCHitsChris()->GetMWPCChamber1Y(0); // Just need fn to get these values directly
  WC1Z1 = GetMWPCHitsChris()->GetMWPCChamber1Z(0);
  WC1X2 = GetMWPCHitsChris()->GetMWPCChamber1X(1);
  WC1Y2 = GetMWPCHitsChris()->GetMWPCChamber1Y(1);
  WC1Z2 = GetMWPCHitsChris()->GetMWPCChamber1Z(1);

  return Theta1, Theta2, Phi1, Phi2, z1, z2, E1, E2, dE1, dE2, WC1X1, WC1X2, WC1Y1, WC1Y2, WC1Z1, WC1Z2; // Returns various quantities used in later functions
}

Int_t PNeutPol_Polarimeter::DetectorCheck()
{
    Detectors1 = GetTracks()->GetDetectors(0); //Gets number for detectors that registered hits
    Detectors2 = GetTracks()->GetDetectors(1); // 7 = NaI + PID + MWPC, 5 = NaI + MWPC
    return Detectors1, Detectors2;
}

Double_t PNeutPol_Polarimeter::MCSmearing() // Smear dE values for MC data to represent Energy resolution of PID
{
  dE1 = rGen.Gaus(GetTracks()->GetVetoEnergy(0) , (0.29*(sqrt(GetTracks()->GetVetoEnergy(0)))));
  dE2 = rGen.Gaus(GetTracks()->GetVetoEnergy(1) , (0.29*(sqrt(GetTracks()->GetVetoEnergy(1)))));

  if (dE1 < 0) dE1 = 0.01;
  if (dE2 < 0) dE2 = 0.01;

  return dE1, dE2;
}

Int_t PNeutPol_Polarimeter::MCTrueID()
{

    MCTrueID1 = GetGeant()->GetTrueID(0);
    MCTrueID2 = GetGeant()->GetTrueID(1);
    return MCTrueID1, MCTrueID2;

}

TLorentzVector PNeutPol_Polarimeter::MCTrueVectors()
{

    MCTrueVect1 = GetGeant()->GetTrueVector(0);
    MCTrueVect2 = GetGeant()->GetTrueVector(1);

    return MCTrueVect1, MCTrueVect2;

}

Double_t PNeutPol_Polarimeter::PNProp(Int_t ProtonParticleNumber) // Define properties of proton and neutron from particles that correspond to each
{
  if(ProtonParticleNumber == 1)
    {
      Zp = z1; // First particle is proton, second neutron
      Zn = z2;
      zdiff = abs (Zp - Zn);
      Ep = E1;
      En = E2;
      dEp = dE1;
      dEn = dE2;
      WC1pX = WC1X1;
      WC1pY = WC1Y1;
      WC1pZ = WC1Z1;
      WC1nX = WC1X2;
      WC1nY = WC1Y2;
      WC1nZ = WC1Z2;
    }

  if(ProtonParticleNumber == 2)
    {
      Zp = z2; // First particle is neutron, second is proton
      Zn = z1;
      zdiff = abs (Zp - Zn);
      Ep = E2; // Therefore the quantity mmp is the amount of missing mass we see when we do a kinematics calculation USING the proton
      En = E1;
      dEp = dE2;
      dEn = dE1;
      WC1pX = WC1X2;
      WC1pY = WC1Y2;
      WC1pZ = WC1Z2;
      WC1nX = WC1X1;
      WC1nY = WC1Y1;
      WC1nZ = WC1Z1;
    }

  return Zp, Zn, zdiff, Ep, En, dEp, dEn, WC1nX, WC1nY, WC1nZ, WC1pX, WC1pY, WC1pZ;
}

TLorentzVector PNeutPol_Polarimeter::PNVect(Int_t ProtonParticleNumber) // Define vectors for p and n in similar manner to properties above
{
  if(ProtonParticleNumber == 1)
    {
      GVp = GV1;
      GV2 = GetTracks()->GetVector(1, Mn);
      GVn = GV2;
    }

  if(ProtonParticleNumber == 2)

    {
      GVp = GV2;
      GV1 = GetTracks()->GetVector(0, Mn); // Since we've decided this particle is a neutron, set its mass to Mn
      GVn = GV1; // The neutron vector as measured by the vertex information
    }

  return GVp, GVn;
}

TVector3 PNeutPol_Polarimeter::WC3Vectors()
{
    WC3Vectp.SetXYZ(WC1pX, WC1pY, WC1pZ);
    WC3Vectn.SetXYZ(WC1nX, WC1nY, WC1nZ);

    return WC3Vectp, WC3Vectn;
}

Double_t PNeutPol_Polarimeter::WCAngles()
{
  WCThetap = WC3Vectp.Theta() * TMath::RadToDeg(); // Angles from WC hit positons
  WCPhip = WC3Vectp.Phi() * TMath::RadToDeg();
  WCThetan = WC3Vectn.Theta() * TMath::RadToDeg();
  WCPhin = WC3Vectn.Phi() * TMath::RadToDeg();
  PhiWCDiff = abs (WCPhip-WCPhin);

  return WCThetap, WCThetan, WCPhip, WCPhin, PhiWCDiff;
}

Double_t PNeutPol_Polarimeter::LabAngles()
{
  Thetap = GVp3.Theta() * TMath::RadToDeg(); // Lab frame angles for proton/neutron
  Phip = GVp3.Phi() * TMath::RadToDeg();
  Thetan = GVn3.Theta() * TMath::RadToDeg();
  Phin = GVn3.Phi() * TMath::RadToDeg();

  return Thetan, Phin, Thetap, Phip;
}

PNeutPol_Polarimeter::PNeutPol_Polarimeter() // Define a load of histograms to fill
{
  time = new GH1("time", 	"time", 	1400, -700, 700);
  time_cut = new GH1("time_cut", 	"time_cut", 	1400, -700, 700);

  Zp_Vert = new GH1("Zp_Vertex", "Proton Z Vertex Distribution", 300, -150, 150 );
  Zn_Vert = new GH1("Zn_Vertex", "Neutron Z Vertex Distribution", 300, -150, 150 );
  Ekp = new GH1( "Ekp", "Proton Energy Distribution", 100, 0, 500 );
  Ekn = new GH1( "Ekn", "Neutron Energy Distribution", 100, 0, 500 );
  EkSum = new GH1( "Ek Sum", "Particle Energy Sum Distribution", 300, 0, 900 );
  Eg = new GH1( "Eg", "Photon Energy Distribution", 200, 100, 1600 );
  ThetaProt = new GH1( "ThetaProt", " Proton Theta Distribution", 180, 0, 180 );
  ThetaNeut = new GH1( "ThetaNeut", " Neutron Theta Distribution", 180, 0, 180 );
  PhiProt = new GH1( "PhiProt", " Proton Phi Distribution", 180, -180, 180 );
  PhiNeut = new GH1( "PhiNeut", " Neutron Phi Distribution", 180, -180, 180 );
  WCPhiDifference = new GH1 ("WCPhiDifference", "WC Phi Difference Between p and n", 180, 0, 360);
  WCThetaProt = new GH1 ("WCThetaProt", "WC Theta for p", 180, 0, 180);
  WCThetaNeut = new GH1 ("WCThetaNeut", "WC Theta for n", 180, 0, 180);
  WCPhiProt = new GH1 ("WCPhiProt", "WC Phi for p", 180, -180, 180);
  WCPhiNeut = new GH1 ("WCPhiNeut", "WC Phi for n", 180, -180, 180);

  E_dE = new GH2("E_dE", "EdE Plot", 125, 0, 500, 125, 0, 7);
  E_dE_p = new GH2("E_dE_p", "EdE Plot for Protons", 125, 0, 500, 125, 0, 7);
  E_dE_n = new GH2("E_dE_n", "EdE Plot for Neutrons", 125, 0, 500, 125, 0, 7);

  EkEg = new GH2("EkEg", "Ek vs Eg for all Particles", 100, 0, 500, 100, 100, 1600);
  EkEg_p = new GH2("EkEg_p", "Ek vs Eg for all Protons", 100, 0, 500, 100, 100, 1600);
  EkEg_n = new GH2("EkEg_n", "Ek vs Eg for all Neutrons", 100, 0, 500, 100, 100, 1600);

  Thetap_ThetaWCp = new GH2 ("Thetap_ThetaWCp", "Theta Track vs Theta WC for Protons", 90, 0, 180, 90, 0 180);
}

void PNeutPol_Polarimeter::FillHists()
{

  Zp_Vert->Fill(Zp, TaggerTime);
  Zn_Vert->Fill(Zn, TaggerTime);
  Ekp->Fill(Ep, TaggerTime);
  Ekn->Fill(En, TaggerTime);
  EkSum->Fill((Ep + En),TaggerTime);
  Eg->Fill(EGamma, TaggerTime);
  ThetaProt->Fill(Thetap, TaggerTime);
  ThetaNeut->Fill(Thetan, TaggerTime);
  PhiProt->Fill(Phip, TaggerTime);
  PhiNeut->Fill(Phin, TaggerTime);
  WCPhiDifference->Fill(PhiWCDiff);
  WCThetaProt->Fill(WCThetap, TaggerTime);
  WCThetaNeut->Fill(WCThetan, TaggerTime);
  WCPhiProt->Fill(WCPhip, TaggerTime);
  WCPhiNeut->Fill(WCPhin, TaggerTime);
  E_dE->Fill(Ep, dEp, TaggerTime);
  E_dE->Fill(En, dEn, TaggerTime);
  E_dE_p->Fill(Ep, dEp, TaggerTime);
  E_dE_n->Fill(En, dEn, TaggerTime);
  EkEg->Fill(Ep, EGamma, TaggerTime);
  EkEg->Fill(En, EGamma, TaggerTime);
  EkEg_p->Fill(Ep, EGamma, TaggerTime);
  EkEg_n->Fill(En, EGamma, TaggerTime);
  Thetap_ThetaWCp->Fill(Thetap, WCThetap, TaggerTime)

}

void PNeutPol_Polarimeter::MCHists()
{

}

Bool_t	PNeutPol_Polarimeter::Write(){
  // Write some TH1s - currently none to write so commented out
  // GTreeManager::Write(TaggerAccScal); // This spams the file with 2500+ histograms of "TaggerAccScal" so commented out

  // Write all GH1's easily

  GTreeManager::Write();
}
