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

  d = 52; // Distance from centre of target to centre of Polarimeter
  NP = 0; // Set number of Protons to 0 before checking
  NPi = 0; // Set number of pions to 0 before checking
  NRoo = 0; // Set number of Rootinos to 0 before checking
  Mn = 939.565; // Mass of neutron in MeV
  Mp = 938.272; // Mass of proton in MeV

  Cut_CB_proton = OpenCutFile("configfiles/cuts/CB_DeltaE-E_Proton_29_07_15.root", "Proton"); // These will need adjusting with new Acqu files
  Cut_proton = Cut_CB_proton;
  Cut_CB_pion = OpenCutFile("configfiles/cuts/CB_DeltaE-E_Pion_29_07_15.root", "Pion");
  Cut_pion = Cut_CB_pion;
  cout << endl;

  MCDataCheck();

  if (MCData == kTRUE) MCHists();

  TraverseValidEvents(); // This loops over each event as in old file and calls ProcessEvent() each loop

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
  if((Detectors1 == 7) && (Detectors2 == 1))
  {
    Proton1 = kTRUE;
    Proton2 = kFALSE;
  }

  // If track 2 only gives signals in MWPC it is the neutron
  else if((Detectors1 == 1) && (Detectors2 == 7))
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

  for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++){

    TaggerTime = GetTagger()->GetTaggedTime(j); // Get tagged time for event
    EGamma = (GetTagger()->GetTaggedEnergy(j)); // Get Photon energy for event
    if (EGamma < 400 || EGamma > 600) continue;

    if (Proton1 == kTRUE)
    {
        PNProp(1);
        PNVect(1);
        PIDEle = GetTracks()->GetCentralVeto(0);
    }

    else if (Proton2 == kTRUE)
    {
        PNProp(2);
        PNVect(2);
        PIDEle = GetTracks()->GetCentralVeto(1);
    }

    else
    {
        continue;
    }

    if( Zp > 40 || Zp < -50) continue;

    GVp3 = GVp.Vect(); // Generate some 3-vectors from the 4-vectors we have
    GVn3 = GVn.Vect();
    WC3Vectors(WC1pX, WC1pY, WC1pZ, WC1nX, WC1nY, WC1nZ);
    WCAngles(WC3Vectp, WC3Vectn);

    LabAngles(); // Get angles in lab based on track info

    KinEp = CalcKinEnergy(WCThetap, EGamma);
    KinEDiff = abs(KinEp - Ep);

    if (KinEDiff > 100) continue; // If difference between CB energy and calculated Energy for proton > 100MeV continue

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


  }
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
  WC1X1 = GetMWPCHitsChris()->GetMWPCChamber1X(0);
  WC1Y1 = GetMWPCHitsChris()->GetMWPCChamber1Y(0);
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

TVector3 PNeutPol_Polarimeter::WC3Vectors(Double_t WCpX, Double_t WCpY, Double_t WCpZ, Double_t WCnX, Double_t WCnY, Double_t WCnZ)
{
    WC3Vectp.SetXYZ(WC1pX, WC1pY, WC1pZ);
    WC3Vectn.SetXYZ(WC1nX, WC1nY, WC1nZ);

    return WC3Vectp, WC3Vectn;
}

Double_t PNeutPol_Polarimeter::WCAngles(TVector3 MWPCp3Vector, TVector3 MWPCn3Vector)
{
  WCThetap = MWPCp3Vector.Theta() * TMath::RadToDeg(); // Angles from WC hit positons
  WCPhip = MWPCp3Vector.Phi() * TMath::RadToDeg();
  WCThetan = MWPCn3Vector.Theta() * TMath::RadToDeg();
  WCPhin = MWPCn3Vector.Phi() * TMath::RadToDeg();
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
  time = new TH1D("time", 	"time", 	1400, -700, 700);
  time_cut = new TH1D("time_cut", 	"time_cut", 	1400, -700, 700);

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
  EpKin = new GH1 ("EpKin", "Ep Calculated from Ep/Thetap", 100, 0, 500);
  EpEKinDiff = new GH1 ("EpEKinDiff", "Difference Between EpKin & Ep", 100, 0, 500);
  WCXp = new GH1("WCXp", "WC X Position for Proton", 200, -100, 100);
  WCYp = new GH1("WCYp", "WC Y Position for Proton", 200, -100, 100);
  WCZp = new GH1("WCZp", "WC Z Position for Proton", 200, -500, 500);
  WCXn = new GH1("WCXn", "WC X Position for Neutron", 200, -100, 100);
  WCYn = new GH1("WCYn", "WC Y Position for Neutron", 200, -100, 100);
  WCZn = new GH1("WCZn", "WC Z Position for Neutron", 200, -500, 500);

  E_dE = new GH2("E_dE", "EdE Plot", 100, 0, 500, 100, 0, 5);
  E_dE_ThetaCut = new GH2 ("E_dE_ThetaCut", "EdE Plot with Theta Cut (35 -45)", 100, 0, 500, 100, 0, 5);
  ECB_dE = new GH2 ("ECB_dE", "CBEdE Plot", 100, 0, 500, 100, 0, 5);
  ECB_dE_ThetaCut = new GH2 ("ECB_dE_ThetaCut", "CBEdE Plot with Theta (35 -45)", 100, 0, 500, 100, 0, 5);

  EpEpKinDiff= new GH2 ("EpEpkinDiff", "Proton CB Energy as fn of EpEKinDiff", 100, 0, 500, 100, 0, 200);

  E_dE_Ele1 = new GH2("E_dE_Ele1", "EdE Plot for PID Element 1", 125, 0, 500, 125, 0, 5);
  E_dE_Ele2 = new GH2("E_dE_Ele2", "EdE Plot for PID Element 2", 125, 0, 500, 125, 0, 5);
  E_dE_Ele3 = new GH2("E_dE_Ele3", "EdE Plot for PID Element 3", 125, 0, 500, 125, 0, 5);
  E_dE_Ele4 = new GH2("E_dE_Ele4", "EdE Plot for PID Element 4", 125, 0, 500, 125, 0, 5);
  E_dE_Ele5 = new GH2("E_dE_Ele5", "EdE Plot for PID Element 5", 125, 0, 500, 125, 0, 5);
  E_dE_Ele6 = new GH2("E_dE_Ele6", "EdE Plot for PID Element 6", 125, 0, 500, 125, 0, 5);
  //E_dE_Ele7 = new GH2("E_dE_Ele7", "EdE Plot for PID Element 7", 125, 0, 500, 125, 0, 5);
  //E_dE_Ele8 = new GH2("E_dE_Ele8", "EdE Plot for PID Element 8", 125, 0, 500, 125, 0, 5);
  //E_dE_Ele9 = new GH2("E_dE_Ele9", "EdE Plot for PID Element 9", 125, 0, 500, 125, 0, 5);
  //E_dE_Ele10 = new GH2("E_dE_Ele10", "EdE Plot for PID Element 10", 125, 0, 500, 125, 0, 5);
  //E_dE_Ele11 = new GH2("E_dE_Ele11", "EdE Plot for PID Element 11", 125, 0, 500, 125, 0, 5);
  //E_dE_Ele12 = new GH2("E_dE_Ele12", "EdE Plot for PID Element 12", 125, 0, 500, 125, 0, 5);
  //E_dE_Ele13 = new GH2("E_dE_Ele13", "EdE Plot for PID Element 13", 125, 0, 500, 125, 0, 5);
  //E_dE_Ele14 = new GH2("E_dE_Ele14", "EdE Plot for PID Element 14", 125, 0, 500, 125, 0, 5);
  //E_dE_Ele15 = new GH2("E_dE_Ele15", "EdE Plot for PID Element 15", 125, 0, 500, 125, 0, 5);
  //E_dE_Ele16 = new GH2("E_dE_Ele16", "EdE Plot for PID Element 16", 125, 0, 500, 125, 0, 5);
  //E_dE_Ele17 = new GH2("E_dE_Ele17", "EdE Plot for PID Element 17", 125, 0, 500, 125, 0, 5);
  //E_dE_Ele18 = new GH2("E_dE_Ele18", "EdE Plot for PID Element 18", 125, 0, 500, 125, 0, 5);
  //E_dE_Ele19 = new GH2("E_dE_Ele19", "EdE Plot for PID Element 19", 125, 0, 500, 125, 0, 5);
  //E_dE_Ele20 = new GH2("E_dE_Ele20", "EdE Plot for PID Element 20", 125, 0, 500, 125, 0, 5);
  //E_dE_Ele21 = new GH2("E_dE_Ele21", "EdE Plot for PID Element 21", 125, 0, 500, 125, 0, 5);
  //E_dE_Ele22 = new GH2("E_dE_Ele22", "EdE Plot for PID Element 22", 125, 0, 500, 125, 0, 5);
  //E_dE_Ele23 = new GH2("E_dE_Ele23", "EdE Plot for PID Element 23", 125, 0, 500, 125, 0, 5);
  //E_dE_Ele24 = new GH2("E_dE_Ele24", "EdE Plot for PID Element 24", 125, 0, 500, 125, 0, 5);

  ThetaPidE = new GH2("ThetaPidE", "Theta as fn of PID Energy", 100, 0, 180, 100, 0, 5);
  ThetaPidE1 = new GH2("ThetaPidE1", "Theta as fn of PID Energy for Element 1", 100, 0, 180, 100, 0, 5);
  ThetaPidE2 = new GH2("ThetaPidE2", "Theta as fn of PID Energy for Element 2", 100, 0, 180, 100, 0, 5);
  ThetaPidE3 = new GH2("ThetaPidE3", "Theta as fn of PID Energy for Element 3", 100, 0, 180, 100, 0, 5);
  ThetaPidE4 = new GH2("ThetaPidE4", "Theta as fn of PID Energy for Element 4", 100, 0, 180, 100, 0, 5);
  ThetaPidE5 = new GH2("ThetaPidE5", "Theta as fn of PID Energy for Element 5", 100, 0, 180, 100, 0, 5);
  ThetaPidE6= new GH2("ThetaPidE6", "Theta as fn of PID Energy for Element 6", 100, 0, 180, 100, 0, 5);
  //ThetaPidE7 = new GH2("ThetaPidE7", "Theta as fn of PID Energy for Element 7", 100, 0, 180, 100, 0, 5);
  //ThetaPidE8 = new GH2("ThetaPidE8", "Theta as fn of PID Energy for Element 8", 100, 0, 180, 100, 0, 5);
  //ThetaPidE9 = new GH2("ThetaPidE9", "Theta as fn of PID Energy for Element 9", 100, 0, 180, 100, 0, 5);
  //ThetaPidE10 = new GH2("ThetaPidE10", "Theta as fn of PID Energy for Element 10", 100, 0, 180, 100, 0, 5);
  //ThetaPidE11 = new GH2("ThetaPidE11", "Theta as fn of PID Energy for Element 11", 100, 0, 180, 100, 0, 5);
  //ThetaPidE12 = new GH2("ThetaPidE12", "Theta as fn of PID Energy for Element 12", 100, 0, 180, 100, 0, 5);
  //ThetaPidE13 = new GH2("ThetaPidE13", "Theta as fn of PID Energy for Element 13", 100, 0, 180, 100, 0, 5);
  //ThetaPidE14 = new GH2("ThetaPidE14", "Theta as fn of PID Energy for Element 14", 100, 0, 180, 100, 0, 5);
  //ThetaPidE15 = new GH2("ThetaPidE15", "Theta as fn of PID Energy for Element 15", 100, 0, 180, 100, 0, 5);
  //ThetaPidE16 = new GH2("ThetaPidE16", "Theta as fn of PID Energy for Element 16", 100, 0, 180, 100, 0, 5);
  //ThetaPidE17 = new GH2("ThetaPidE17", "Theta as fn of PID Energy for Element 17", 100, 0, 180, 100, 0, 5);
  //ThetaPidE18 = new GH2("ThetaPidE18", "Theta as fn of PID Energy for Element 18", 100, 0, 180, 100, 0, 5);
  //ThetaPidE19 = new GH2("ThetaPidE19", "Theta as fn of PID Energy for Element 19", 100, 0, 180, 100, 0, 5);
  //ThetaPidE20 = new GH2("ThetaPidE20", "Theta as fn of PID Energy for Element 20", 100, 0, 180, 100, 0, 5);
  //ThetaPidE21 = new GH2("ThetaPidE21", "Theta as fn of PID Energy for Element 21", 100, 0, 180, 100, 0, 5);
  //ThetaPidE22 = new GH2("ThetaPidE22", "Theta as fn of PID Energy for Element 22", 100, 0, 180, 100, 0, 5);
  //ThetaPidE23 = new GH2("ThetaPidE23", "Theta as fn of PID Energy for Element 23", 100, 0, 180, 100, 0, 5);
  //ThetaPidE24 = new GH2("ThetaPidE24", "Theta as fn of PID Energy for Element 24", 100, 0, 180, 100, 0, 5);
}

void PNeutPol_Polarimeter::FillHists()
{
  time->Fill(TaggerTime);
  if (-5 < TaggerTime && TaggerTime < 20) time_cut->Fill(TaggerTime);

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
  E_dE->Fill(KinEp, dEp, TaggerTime);
  ECB_dE->Fill(Ep, dEp, TaggerTime);
  EpKin->Fill(KinEp, TaggerTime);
  EpEKinDiff->Fill(KinEDiff, TaggerTime);
  EpEpKinDiff->Fill(Ep, KinEDiff, TaggerTime);
  WCXp->Fill(WC1pX, TaggerTime);
  WCYp->Fill(WC1pY, TaggerTime);
  WCZp->Fill(WC1pZ, TaggerTime);
  WCXn->Fill(WC1nX, TaggerTime);
  WCYn->Fill(WC1nY, TaggerTime);
  WCZn->Fill(WC1nZ, TaggerTime);

  if( 35 < Thetap && Thetap < 45) E_dE_ThetaCut->Fill(KinEp, dEp, TaggerTime);
  if( 35 < Thetap && Thetap < 45) ECB_dE_ThetaCut->Fill(Ep, dEp, TaggerTime);

  if(PIDEle == 0) E_dE_Ele1->Fill(KinEp, dEp, TaggerTime);
  if(PIDEle == 1) E_dE_Ele2->Fill(KinEp, dEp, TaggerTime);
  if(PIDEle == 2) E_dE_Ele3->Fill(KinEp, dEp, TaggerTime);
  if(PIDEle == 3) E_dE_Ele4->Fill(KinEp, dEp, TaggerTime);
  if(PIDEle == 4) E_dE_Ele5->Fill(KinEp, dEp, TaggerTime);
  if(PIDEle == 5) E_dE_Ele6->Fill(KinEp, dEp, TaggerTime);
  //if(PIDEle == 6) E_dE_Ele7->Fill(KinEp, dEp, TaggerTime);
  //if(PIDEle == 7) E_dE_Ele8->Fill(KinEp, dEp, TaggerTime);
  //if(PIDEle == 8) E_dE_Ele9->Fill(KinEp, dEp, TaggerTime);
  //if(PIDEle == 9) E_dE_Ele10->Fill(KinEp, dEp, TaggerTime);
  //if(PIDEle == 10) E_dE_Ele11->Fill(KinEp, dEp, TaggerTime);
  //if(PIDEle == 11) E_dE_Ele12->Fill(KinEp, dEp, TaggerTime);
  //if(PIDEle == 12) E_dE_Ele13->Fill(KinEp, dEp, TaggerTime);
  //if(PIDEle == 13) E_dE_Ele14->Fill(KinEp, dEp, TaggerTime);
  //if(PIDEle == 14) E_dE_Ele15->Fill(KinEp, dEp, TaggerTime);
  //if(PIDEle == 15) E_dE_Ele16->Fill(KinEp, dEp, TaggerTime);
  //if(PIDEle == 16) E_dE_Ele17->Fill(KinEp, dEp, TaggerTime);
  //if(PIDEle == 17) E_dE_Ele18->Fill(KinEp, dEp, TaggerTime);
  //if(PIDEle == 18) E_dE_Ele19->Fill(KinEp, dEp, TaggerTime);
  //if(PIDEle == 19) E_dE_Ele20->Fill(KinEp, dEp, TaggerTime);
  //if(PIDEle == 20) E_dE_Ele21->Fill(KinEp, dEp, TaggerTime);
  //if(PIDEle == 21) E_dE_Ele22->Fill(KinEp, dEp, TaggerTime);
  //if(PIDEle == 22) E_dE_Ele23->Fill(KinEp, dEp, TaggerTime);
  //if(PIDEle == 23) E_dE_Ele24->Fill(KinEp, dEp, TaggerTime);

  if (200 < KinEp && KinEp < 240){

    ThetaPidE->Fill(Thetap, dEp, TaggerTime);
    if(PIDEle == 0) ThetaPidE1->Fill(Thetap, dEp, TaggerTime);
    if(PIDEle == 1) ThetaPidE2->Fill(Thetap, dEp, TaggerTime);
    if(PIDEle == 2) ThetaPidE3->Fill(Thetap, dEp, TaggerTime);
    if(PIDEle == 3) ThetaPidE4->Fill(Thetap, dEp, TaggerTime);
    if(PIDEle == 4) ThetaPidE5->Fill(Thetap, dEp, TaggerTime);
    if(PIDEle == 5) ThetaPidE6->Fill(Thetap, dEp, TaggerTime);
    //if(PIDEle == 6) ThetaPidE7->Fill(Thetap, dEp, TaggerTime);
    //if(PIDEle == 7) ThetaPidE8->Fill(Thetap, dEp, TaggerTime);
    //if(PIDEle == 8) ThetaPidE9->Fill(Thetap, dEp, TaggerTime);
    //if(PIDEle == 9) ThetaPidE10->Fill(Thetap, dEp, TaggerTime);
    //if(PIDEle == 10) ThetaPidE11->Fill(Thetap, dEp, TaggerTime);
    //if(PIDEle == 11) ThetaPidE12->Fill(Thetap, dEp, TaggerTime);
    //if(PIDEle == 12) ThetaPidE13->Fill(Thetap, dEp, TaggerTime);
    //if(PIDEle == 13) ThetaPidE14->Fill(Thetap, dEp, TaggerTime);
    //if(PIDEle == 14) ThetaPidE15->Fill(Thetap, dEp, TaggerTime);
    //if(PIDEle == 15) ThetaPidE16->Fill(Thetap, dEp, TaggerTime);
    //if(PIDEle == 16) ThetaPidE17->Fill(Thetap, dEp, TaggerTime);
    //if(PIDEle == 17) ThetaPidE18->Fill(Thetap, dEp, TaggerTime);
    //if(PIDEle == 18) ThetaPidE19->Fill(Thetap, dEp, TaggerTime);
    //if(PIDEle == 19) ThetaPidE20->Fill(Thetap, dEp, TaggerTime);
    //if(PIDEle == 20) ThetaPidE21->Fill(Thetap, dEp, TaggerTime);
    //if(PIDEle == 21) ThetaPidE22->Fill(Thetap, dEp, TaggerTime);
    //if(PIDEle == 22) ThetaPidE23->Fill(Thetap, dEp, TaggerTime);
    //if(PIDEle == 23) ThetaPidE24->Fill(Thetap, dEp, TaggerTime);
  }

}

void PNeutPol_Polarimeter::MCHists()
{

}

Bool_t	PNeutPol_Polarimeter::Write(){
  // Write all GH1's easily

  GTreeManager::Write();
}
