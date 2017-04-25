// GoAT Physics analysis to identify neutrons from deuterium photodisintegration
// Various properties of neutrons/protons identified plotted in histograms
// Main aim is to determine spin polarisation of neutrons
// For use on linearly polarised data files

#include "PNeutPol_Polarimeter_Lin.h"

PNeutPol_Polarimeter_Lin::~PNeutPol_Polarimeter_Lin()
{
}

Bool_t	PNeutPol_Polarimeter_Lin::Init()
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

Bool_t	PNeutPol_Polarimeter_Lin::Start()
{
  if(!IsGoATFile())
  {
    cout << "ERROR: Input File is not a GoAT file." << endl;
    return kFALSE;
  }

  SetAsPhysicsFile();

  EventCounter = 0;
  EventCounterTrackCut = 0;
  EventCounterZCut = 0;
  EventCounterCoplanarCut = 0;
  NP = 0; // Set number of Protons to 0 before checking
  NPi = 0; // Set number of pions to 0 before checking
  NRoo = 0; // Set number of Rootinos to 0 before checking
  Mn = 939.565; // Mass of neutron in MeV
  Mp = 938.272; // Mass of proton in MeV
  Md = 1875.613; //Mass of Deuterium in MeV
  Mpi = 139.57018; // Mass of charged pion in MeV
  Deut = TLorentzVector (0., 0., 0., 1875.613); // 4-Vector of Deuterium target, assume at rest
  Neut = TLorentzVector (0., 0., 0., 939.565); // 4-Vector of Deuterium target, assume at rest

  Cut_CB_proton = OpenCutFile("configfiles/cuts/CB_DeltaE-E_Proton_7_12_16.root", "Proton"); // These will need adjusting with new Acqu files
  Cut_proton = Cut_CB_proton;
  Cut_CB_pion = OpenCutFile("configfiles/cuts/CB_DeltaE-E_Pion_29_07_15.root", "Pion");
  Cut_pion = Cut_CB_pion;
  Cut_CB_protonKinGood = OpenCutFile("configfiles/cuts/CB_DeltaE-E_ProtonKinGood_18_03_17.root", "ProtonKinGood"); // These will need adjusting with new Acqu files
  Cut_protonKinGood = Cut_CB_protonKinGood;
  Cut_CB_protonKinBad = OpenCutFile("configfiles/cuts/CB_DeltaE-E_ProtonKinBad_15_12_16.root", "ProtonKinBad");
  Cut_protonKinBad = Cut_CB_protonKinBad;
  cout << endl;

  TraverseValidEvents(); // This loops over each event as in old file and calls ProcessEvent() each loop

  cout << EventCounter << " Events in file " << EventCounterTrackCut << " Events After Track Cut " << EventCounterZCut << " Events after Z cut " << EventCounterCoplanarCut << " Events after Coplanarity Cut" << endl;

  return kTRUE;
}

void	PNeutPol_Polarimeter_Lin::ProcessEvent()
{

  NTrack = GetTracks()->GetNTracks();
  NP = GetProtons()->GetNParticles();
  NPi = GetChargedPions()->GetNParticles();
  NRoo = GetRootinos()->GetNParticles();
  NTag = GetTagger()->GetNTagged();
  EventCounter++;
  if (NRoo !=0) return; // Goes to next event if any "rootinos" found
  if (NTrack !=2) return; // Ensures two track event
  Detectors1 = GetTracks()->GetDetectors(0); //Gets number for detectors that registered hits
  Detectors2 = GetTracks()->GetDetectors(1); // 7 = NaI + PID + MWPC, 5 = NaI + MWPC

  //Condition is to now keep any ALL + CB/MWPC events
  //If track 1 only gives signals in CB it is the neutron

  if((Detectors1 == 7) && (Detectors2 == 5))
  {
    Proton1 = kTRUE;
    Proton2 = kFALSE;
    if (GetTracks()->GetMWPC0Energy(0) == 0) return; // If no hit in first chamber for p drop out
    if (GetTracks()->GetMWPC0Energy(1) == 0) return; // If no hit in first chamber for n drop out
  }

  // If track 2 only gives signals in MWPC and CB it is the neutron
  else if((Detectors1 == 5) && (Detectors2 == 7))
  {
    Proton1 = kFALSE;
    Proton2 = kTRUE;
    if (GetTracks()->GetMWPC0Energy(1) == 0) return; // If no hit in first chamber for p drop out
    if (GetTracks()->GetMWPC0Energy(0) == 0) return; // If no hit in first chamber for n drop out
  }

  // Drop out on ANY other condition (for now)
  else
  {
    return;
  }

  EventCounterTrackCut++;

  if (Proton1 == kTRUE)
  {
    GVp = GetTracks()->GetVector(0, Mp);
    GVn = GetTracks()->GetVector(1, Mn);
    Zp = GetTracks()->GetPseudoVertexZ(0); // First particle is proton, second neutron
    Zn = GetTracks()->GetPseudoVertexZ(1);
    Ep = GetTracks()->GetClusterEnergy(0);
    En = GetTracks()->GetClusterEnergy(1);
    dEp = GetTracks()->GetVetoEnergy(0);
    dEn = GetTracks()->GetVetoEnergy(1);
    WC1pX = GetTracks()->GetMWPC0PosX(0);
    WC1pY = GetTracks()->GetMWPC0PosY(0);
    WC1pZ = GetTracks()->GetMWPC0PosZ(0);
    WC1nX = GetTracks()->GetMWPC0PosX(1);
    WC1nY = GetTracks()->GetMWPC0PosY(1);
    WC1nZ = GetTracks()->GetMWPC0PosZ(1);
  }

  else if (Proton2 == kTRUE)
  {
    GVp = GetTracks()->GetVector(1, Mp);
    GVn = GetTracks()->GetVector(0, Mn);
    Zp = GetTracks()->GetPseudoVertexZ(1); // First particle is neutron, second is proton
    Zn = GetTracks()->GetPseudoVertexZ(0);
    Ep = GetTracks()->GetClusterEnergy(1); // Therefore the quantity mmp is the amount of missing mass we see when we do a kinematics calculation USING the proton
    En = GetTracks()->GetClusterEnergy(0);
    dEp = GetTracks()->GetVetoEnergy(1);
    dEn = GetTracks()->GetVetoEnergy(0);
    WC1pX = GetTracks()->GetMWPC0PosX(1);
    WC1pY = GetTracks()->GetMWPC0PosY(1);
    WC1pZ = GetTracks()->GetMWPC0PosZ(1);
    WC1nX = GetTracks()->GetMWPC0PosX(0);
    WC1nY = GetTracks()->GetMWPC0PosY(0);
    WC1nZ = GetTracks()->GetMWPC0PosZ(0);
  }

  else
  {
    return;
  }

  GVp3 = GVp.Vect(); // Generate some 3-vectors from the 4-vectors we have
  GVn3 = GVn.Vect();
  WC3Vectp.SetXYZ(WC1pX, WC1pY, WC1pZ);
  WC3Vectn.SetXYZ(WC1nX, WC1nY, WC1nZ);
  WCThetap = WC3Vectp.Theta()*TMath::RadToDeg(); // Angles from WC hit positons
  WCThetapRad = WC3Vectp.Theta();
  WCPhip = WC3Vectp.Phi()*TMath::RadToDeg();
  WCPhipRad = WC3Vectp.Phi();
  WCThetan = WC3Vectn.Theta()*TMath::RadToDeg();
  WCPhin = WC3Vectn.Phi()*TMath::RadToDeg();
  PhiWCDiff = abs (WCPhip-WCPhin);

  //if( Zp > 60 || Zp < -60) return; // Particles selected out from other parts tend to be inside anyway, skip this?
  //if( Zp > 200 || Zp < 150) return; // Select out windows

  EventCounterZCut++;

  //if ( PhiWCDiff > 195 || PhiWCDiff < 165) return;

  EventCounterCoplanarCut++;

  for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
  {

    TaggerTime = GetTagger()->GetTaggedTime(j); // Get tagged time for event
    EGamma = (GetTagger()->GetTaggedEnergy(j)); // Get Photon energy for event
    Gamma = TLorentzVector (0., 0., EGamma , EGamma); // 4-Vector of Photon beam
    B = (Deut + Gamma).Beta(); // Calculate Lorentz Beta
    b = TVector3(0., 0., B); // Define boost vector
    GVpB = GVp;
    GVpB.Boost(b); // Boost GVp to CM frame
    ThetapCM = (GVpB.Theta())*TMath::RadToDeg(); // Get Theta of proton in CM frame

    Thetap = GVp3.Theta()*TMath::RadToDeg(); // Lab frame angles for proton/neutron
    Phip = GVp3.Phi()*TMath::RadToDeg();
    Thetan = GVn3.Theta()*TMath::RadToDeg();
    Phin = GVn3.Phi()*TMath::RadToDeg();

    EpCorr = EpPolCorrect(Ep, WCThetap);
    EpDiff = abs(EpCorr - Ep);

    // Gamma(d,p)n
    KinEp = CalcKinEnergy(WCThetap, EGamma, Md, 0., Mp, Mn); // Calculate kin E of proton assuming pn production
    RecKinProton = LProton4VectorKin(KinEp, WCThetapRad, WCPhipRad);
    RecKinNeutron = LNeutron4VectorKin(RecKinProton);
    ThetanRec = (RecKinNeutron.Theta()) * TMath::RadToDeg();
    PhinRec = (RecKinNeutron.Phi()) * TMath::RadToDeg();
    WCZnRec = 72/tan(RecKinNeutron.Theta());

    // Gamma(n,p)Pi (P detected correct)
    // Assume proton track is proton and "neutron" track is from charged pion
    KinEpPi = CalcKinEnergy(WCThetap, EGamma, Mn, 0, Mp, Mpi); // Calculate kin E of proton assuming g(n, p) pi
    RecKinProtonPi = LProton4VectorKin(KinEpPi, WCThetapRad, WCPhipRad); // Get Proton 4 vector from calculated kin E
    RecKinPion = LPion4VectorKin(RecKinProtonPi); // Get Pion 4 vector from 4 momenta conservation
    ThetaPiRec = (RecKinPion.Theta())*TMath::RadToDeg();
    PhiPiRec = (RecKinPion.Phi())*TMath::RadToDeg();
    ThetaPiRecDiff = abs(ThetaPiRec - Thetan);

    // Gamma(n,p)Pi (Pion detected correct)
    // Assume proton track is pion and "neutron" track is from proton
    KinPi = CalcKinEnergy(WCThetap, EGamma, Mn, 0, Mpi, Mp); // Calculate kin E of pion
    RecKinPionP = LProton4VectorKin(KinPi, WCThetapRad, WCPhipRad); // Get Pion 4 vector from calculated kinE
    RecKinPPi = LPion4VectorKin(RecKinPionP); // Get Proton 4 vector from 4 momenta conservation
    ThetapRec = (RecKinPPi.Theta())*TMath::RadToDeg();
    PhipRec = (RecKinPPi.Phi())*TMath::RadToDeg();
    ThetapRecDiff = abs (ThetapRec - Thetan);

    KinEDiff = KinEp - EpCorr;

    RecProtonEpCorr = LProton4VectorKin(EpCorr, WCThetapRad, WCPhipRad);
    RecNeutronEpCorr = LNeutron4VectorKin(RecProtonEpCorr);
    MMpEpCorr = RecNeutronEpCorr.M();
    RecProtonEpCorr3 = RecProtonEpCorr.Vect();
    RecNeutronEpCorr3 = RecNeutronEpCorr.Vect();

    P3Vect = RecKinProton.Vect();
    N3Vect = RecKinNeutron.Vect();
    OpeningAngle = (N3Vect.Angle(GVn3))*TMath::RadToDeg();

    ThetanDiff = abs(ThetanRec - WCThetan);
    PhinDiff = abs(PhinRec - WCPhin);

    TVector3 ScattAngles = ScatteredFrameAngles(RecNeutronEpCorr3, GVp3, GVn3, Gamma);
    ScattTheta = ScattAngles(0); // Theta is 1st component in vector fn returns above
    ScattPhi = ScattAngles(1); // Phi is 2nd component

    if(Cut_protonKinGood -> IsInside(KinEp, dEp) == kFALSE) continue; // If E loss correct proton is NOT inside p banana drop out
    if(ThetaPiRec > 20) continue;
    //if ( 850 > MMpEpCorr || 1050 < MMpEpCorr) continue;
    if (ScattTheta > 60) continue;
    //if (ScattPhi > 170) continue; // Exclude values  at edges for now
    //if (ScattPhi < -170) continue;
    //if (ScattPhi < -165 || ScattPhi > 165) continue;

    //if (abs(KinEDiff) > 100) continue; // If difference between CB energy and calculated Energy for proton > 100MeV continue

    //k++;
    FillHists(); // Fill histograms with data generated
  }
}

void	PNeutPol_Polarimeter_Lin::ProcessScalerRead()
{
	// Fill Tagger Scalers // Currently this seems to fill the file with loads of "TaggerAccScal" histograms
	//FillScalers(GetTC_scaler_min(),GetTC_scaler_max(),TaggerAccScal); // Don't know if these are needed so cut out for now
}

TCutG*	PNeutPol_Polarimeter_Lin::OpenCutFile(Char_t* filename, Char_t* cutname)
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

TLorentzVector PNeutPol_Polarimeter_Lin::LProton4VectorKin(Double_t KinE, Double_t Theta, Double_t Phi)
{
    EpTot = KinE + Mp;
    Pp = sqrt (TMath::Power(EpTot,2) - TMath::Power(Mp,2));
    Ppx = Pp*sin(Theta)*cos(Phi);
    Ppy = Pp*sin(Theta)*sin(Phi);
    Ppz = Pp*cos(Theta);
    P4Vect.SetE(EpTot);
    P4Vect.SetPx(Ppx);
    P4Vect.SetPy(Ppy);
    P4Vect.SetPz(Ppz);

    return P4Vect;
}

TLorentzVector PNeutPol_Polarimeter_Lin::LNeutron4VectorKin(TLorentzVector ProtonKinVector)
{

    N4Vect = (Gamma + Deut) - ProtonKinVector;

    return N4Vect;
}

TLorentzVector PNeutPol_Polarimeter_Lin::LPion4VectorKin(TLorentzVector ProtonKinVector)
{

    Pi4Vect = (Gamma + Neut) - ProtonKinVector;

    return Pi4Vect;
}

PNeutPol_Polarimeter_Lin::PNeutPol_Polarimeter_Lin() // Define a load of histograms to fill
{
  time = new TH1D("time", 	"time", 	1400, -700, 700);
  time_cut = new TH1D("time_cut", 	"time_cut", 	1400, -700, 700);

  Eg = new GH1( "Eg", "Photon Energy Distribution", 200, 100, 1600 );
  WCPhiDifference = new GH1 ("WCPhiDifference", "WC Phi Difference Between p and n", 180, 0, 360);
  EpKin = new GH1 ("EpKin", "Ep Calculated from Ep/Thetap", 100, 0, 500);
  EpCorrected = new GH1 ("EpCorrected", "Ep Corrected for Energy Loss in Polarimeter ", 100, 0, 500);
  OAngle = new GH1 ("OAngle", "Opening Angle between P and N Vectors", 180, 0, 180);
  WCZnRecon = new GH1 ("WCZnRecon", "WCZ Hit Position from Reconstructed n Vector", 200, 0, 400);

  ThetaSc =  new GH1( "Theta_Scattered", "Scattered Proton Theta Distribution in Rotated Frame", 180, 0, 180 );
  PhiSc = new GH1( "Phi_Scattered", "Scattered Proton Phi Distribution in Rotated Frame", 90, -180, 180 );

  EpKinEpCorrDiff = new GH1("EpKinEpCorrDiff", "Difference Between EpKin and EpCorr", 300, -300, 300);
  EpEpCorrDiff = new GH1("EpEpCorrDiff", "Difference Between Ep and EpCorr", 200, 0, 200);

  MMpEpCorrected = new GH1 ("MMpEpCorrected", "Missing mass seen by Proton (E Loss Corrected)", 400, 0, 2000);

  MMpEpCorrectedCut =  new GH1 ("MMpEpCorrectedCut", "Missing mass seen by Proton (E Loss Corrected, P Banana Cut)", 400, 0, 2000);
  OAngleCut = new GH1 ("OAngleCut", "Opening Angle between P and N Vectors (P Banana Cut)", 180, 0, 180);
  OAngleCut200400 = new GH1 ("OAngleCut200400", "Opening Angle between P and N Vectors (P Banana Cut, 200-400MeV Gamma)", 180, 0, 180);
  EgCut = new GH1( "EgCut", "Photon Energy Distribution (P Banana Cut)", 400, 100, 1600 );

  ZpDist = new GH1 ("ZpDist", "Proton Pseudo Z Vertex Distribution", 200, -400, 400);
  ZpPhiScatNeg180 = new GH1("ZpPhiScatNeg180", "Proton Pseudo Vertex Z for events with PhiSc ~ -ve180", 200, -200, 200);
  ZpPhiScat0 = new GH1("ZpPhiScat0", "Proton Pseudo Vertex Z for events with PhiSc ~ 0", 200, -200, 200);
  ZpPhiScatPos180 = new GH1("ZpPhiScatPos180", "Proton Pseudo Vertex Z for events with PhiSc ~ 180", 200, -200, 200);

  // MMp across photon E bins
  MMp200300 = new GH1("MMp200300", "Missing mass as seen by Proton (200-300MeV Photon Energy)", 400, 0, 2000);
  MMp300400 = new GH1("MMp300400", "Missing mass as seen by Proton (300-400MeV Photon Energy)", 400, 0, 2000);
  MMp400500 = new GH1("MMp400500", "Missing mass as seen by Proton (400-500MeV Photon Energy)", 400, 0, 2000);
  MMp500600 = new GH1("MMp500600", "Missing mass as seen by Proton (500-600MeV Photon Energy)", 400, 0, 2000);
  MMp600700 = new GH1("MMp600700", "Missing mass as seen by Proton (600-700MeV Photon Energy)", 400, 0, 2000);
  MMp700800 = new GH1("MMp700800", "Missing mass as seen by Proton (700-800MeV Photon Energy)", 400, 0, 2000);
  MMp800900 = new GH1("MMp800900", "Missing mass as seen by Proton (800-900MeV Photon Energy)", 400, 0, 2000);

  // Angles of neutron in scattered frame across EGamma bins
  PhiSc410 = new GH1("Phi_Scattered_410MeV", "Scattered Proton Phi Distribution in Rotated Frame for Photon Energies of 410pm10MeV", 2, -180, 180);
  PhiSc430 = new GH1("Phi_Scattered_430MeV", "Scattered Proton Phi Distribution in Rotated Frame for Photon Energies of 430pm10MeV", 2, -180, 180);
  PhiSc450 = new GH1("Phi_Scattered_450MeV", "Scattered Proton Phi Distribution in Rotated Frame for Photon Energies of 450pm10MeV", 2, -180, 180);
  PhiSc470 = new GH1("Phi_Scattered_470MeV", "Scattered Proton Phi Distribution in Rotated Frame for Photon Energies of 470pm10MeV", 2, -180, 180);
  PhiSc490 = new GH1("Phi_Scattered_490MeV", "Scattered Proton Phi Distribution in Rotated Frame for Photon Energies of 490pm10MeV", 2, -180, 180);
  PhiSc510 = new GH1("Phi_Scattered_510MeV", "Scattered Proton Phi Distribution in Rotated Frame for Photon Energies of 510pm10MeV", 2, -180, 180);
  PhiSc530 = new GH1("Phi_Scattered_530MeV", "Scattered Proton Phi Distribution in Rotated Frame for Photon Energies of 530pm10MeV", 2, -180, 180);
  PhiSc550 = new GH1("Phi_Scattered_550MeV", "Scattered Proton Phi Distribution in Rotated Frame for Photon Energies of 550pm10MeV", 2, -180, 180);
  PhiSc570 = new GH1("Phi_Scattered_570MeV", "Scattered Proton Phi Distribution in Rotated Frame for Photon Energies of 570pm10MeV", 2, -180, 180);
  PhiSc590 = new GH1("Phi_Scattered_590MeV", "Scattered Proton Phi Distribution in Rotated Frame for Photon Energies of 590pm10MeV", 2, -180, 180);
  PhiSc610 = new GH1("Phi_Scattered_610MeV", "Scattered Proton Phi Distribution in Rotated Frame for Photon Energies of 610pm10MeV", 2, -180, 180);
  PhiSc630 = new GH1("Phi_Scattered_630MeV", "Scattered Proton Phi Distribution in Rotated Frame for Photon Energies of 630pm10MeV", 2, -180, 180);

  // Proton Phi dists across EGamma bins
  Phip410CM1 = new GH1("Phip_410MeVCM1", "Proton Phi Distribution for Photon Energies of 410pm10MeV (ThetaCM0-30)", 10, -180, 180);
  Phip430CM1 = new GH1("Phip_430MeVCM1", "Proton Phi Distribution for Photon Energies of 430pm10MeV (ThetaCM0-30)", 10, -180, 180);
  Phip450CM1 = new GH1("Phip_450MeVCM1", "Proton Phi Distribution for Photon Energies of 450pm10MeV (ThetaCM0-30)", 10, -180, 180);
  Phip470CM1 = new GH1("Phip_470MeVCM1", "Proton Phi Distribution for Photon Energies of 470pm10MeV (ThetaCM0-30)", 10, -180, 180);
  Phip490CM1 = new GH1("Phip_490MeVCM1", "Proton Phi Distribution for Photon Energies of 490pm10MeV (ThetaCM0-30)", 10, -180, 180);
  Phip510CM1 = new GH1("Phip_510MeVCM1", "Proton Phi Distribution for Photon Energies of 510pm10MeV (ThetaCM0-30)", 10, -180, 180);
  Phip530CM1 = new GH1("Phip_530MeVCM1", "Proton Phi Distribution for Photon Energies of 530pm10MeV (ThetaCM0-30)", 10, -180, 180);
  Phip550CM1 = new GH1("Phip_550MeVCM1", "Proton Phi Distribution for Photon Energies of 550pm10MeV (ThetaCM0-30)", 10, -180, 180);
  Phip570CM1 = new GH1("Phip_570MeVCM1", "Proton Phi Distribution for Photon Energies of 570pm10MeV (ThetaCM0-30)", 10, -180, 180);
  Phip590CM1 = new GH1("Phip_590MeVCM1", "Proton Phi Distribution for Photon Energies of 590pm10MeV (ThetaCM0-30)", 10, -180, 180);
  Phip610CM1 = new GH1("Phip_610MeVCM1", "Proton Phi Distribution for Photon Energies of 610pm10MeV (ThetaCM0-30)", 10, -180, 180);
  Phip630CM1 = new GH1("Phip_630MeVCM1", "Proton Phi Distribution for Photon Energies of 630pm10MeV (ThetaCM0-30)", 10, -180, 180);

  // Angles of neutron in scattered frame across EGamma bins
  Phin410CM1 = new GH1("Phin_410MeVCM1", "Neutron Phi Distribution for Photon Energies of 410pm10MeV (ThetaCM0-30)", 10, -180, 180);
  Phin430CM1 = new GH1("Phin_430MeVCM1", "Neutron Phi Distribution for Photon Energies of 430pm10MeV (ThetaCM0-30)", 10, -180, 180);
  Phin450CM1 = new GH1("Phin_450MeVCM1", "Neutron Phi Distribution for Photon Energies of 450pm10MeV (ThetaCM0-30)", 10, -180, 180);
  Phin470CM1 = new GH1("Phin_470MeVCM1", "Neutron Phi Distribution for Photon Energies of 470pm10MeV (ThetaCM0-30)", 10, -180, 180);
  Phin490CM1 = new GH1("Phin_490MeVCM1", "Neutron Phi Distribution for Photon Energies of 490pm10MeV (ThetaCM0-30)", 10, -180, 180);
  Phin510CM1 = new GH1("Phin_510MeVCM1", "Neutron Phi Distribution for Photon Energies of 510pm10MeV (ThetaCM0-30)", 10, -180, 180);
  Phin530CM1 = new GH1("Phin_530MeVCM1", "Neutron Phi Distribution for Photon Energies of 530pm10MeV (ThetaCM0-30)", 10, -180, 180);
  Phin550CM1 = new GH1("Phin_550MeVCM1", "Neutron Phi Distribution for Photon Energies of 550pm10MeV (ThetaCM0-30)", 10, -180, 180);
  Phin570CM1 = new GH1("Phin_570MeVCM1", "Neutron Phi Distribution for Photon Energies of 570pm10MeV (ThetaCM0-30)", 10, -180, 180);
  Phin590CM1 = new GH1("Phin_590MeVCM1", "Neutron Phi Distribution for Photon Energies of 590pm10MeV (ThetaCM0-30)", 10, -180, 180);
  Phin610CM1 = new GH1("Phin_610MeVCM1", "Neutron Phi Distribution for Photon Energies of 610pm10MeV (ThetaCM0-30)", 10, -180, 180);
  Phin630CM1 = new GH1("Phin_630MeVCM1", "Neutron Phi Distribution for Photon Energies of 630pm10MeV (ThetaCM0-30)", 10, -180, 180);

  // Proton Phi dists across EGamma bins
  Phip410CM2 = new GH1("Phip_410MeVCM2", "Proton Phi Distribution for Photon Energies of 410pm10MeV (ThetaCM30-60)", 10, -180, 180);
  Phip430CM2 = new GH1("Phip_430MeVCM2", "Proton Phi Distribution for Photon Energies of 430pm10MeV (ThetaCM30-60)", 10, -180, 180);
  Phip450CM2 = new GH1("Phip_450MeVCM2", "Proton Phi Distribution for Photon Energies of 450pm10MeV (ThetaCM30-60)", 10, -180, 180);
  Phip470CM2 = new GH1("Phip_470MeVCM2", "Proton Phi Distribution for Photon Energies of 470pm10MeV (ThetaCM30-60)", 10, -180, 180);
  Phip490CM2 = new GH1("Phip_490MeVCM2", "Proton Phi Distribution for Photon Energies of 490pm10MeV (ThetaCM30-60)", 10, -180, 180);
  Phip510CM2 = new GH1("Phip_510MeVCM2", "Proton Phi Distribution for Photon Energies of 510pm10MeV (ThetaCM30-60)", 10, -180, 180);
  Phip530CM2 = new GH1("Phip_530MeVCM2", "Proton Phi Distribution for Photon Energies of 530pm10MeV (ThetaCM30-60)", 10, -180, 180);
  Phip550CM2 = new GH1("Phip_550MeVCM2", "Proton Phi Distribution for Photon Energies of 550pm10MeV (ThetaCM30-60)", 10, -180, 180);
  Phip570CM2 = new GH1("Phip_570MeVCM2", "Proton Phi Distribution for Photon Energies of 570pm10MeV (ThetaCM30-60)", 10, -180, 180);
  Phip590CM2 = new GH1("Phip_590MeVCM2", "Proton Phi Distribution for Photon Energies of 590pm10MeV (ThetaCM30-60)", 10, -180, 180);
  Phip610CM2 = new GH1("Phip_610MeVCM2", "Proton Phi Distribution for Photon Energies of 610pm10MeV (ThetaCM30-60)", 10, -180, 180);
  Phip630CM2 = new GH1("Phip_630MeVCM2", "Proton Phi Distribution for Photon Energies of 630pm10MeV (ThetaCM30-60)", 10, -180, 180);

  // Angles of neutron in scattered frame across EGamma bins
  Phin410CM2 = new GH1("Phin_410MeVCM2", "Neutron Phi Distribution for Photon Energies of 410pm10MeV (ThetaCM30-60)", 10, -180, 180);
  Phin430CM2 = new GH1("Phin_430MeVCM2", "Neutron Phi Distribution for Photon Energies of 430pm10MeV (ThetaCM30-60)", 10, -180, 180);
  Phin450CM2 = new GH1("Phin_450MeVCM2", "Neutron Phi Distribution for Photon Energies of 450pm10MeV (ThetaCM30-60)", 10, -180, 180);
  Phin470CM2 = new GH1("Phin_470MeVCM2", "Neutron Phi Distribution for Photon Energies of 470pm10MeV (ThetaCM30-60)", 10, -180, 180);
  Phin490CM2 = new GH1("Phin_490MeVCM2", "Neutron Phi Distribution for Photon Energies of 490pm10MeV (ThetaCM30-60)", 10, -180, 180);
  Phin510CM2 = new GH1("Phin_510MeVCM2", "Neutron Phi Distribution for Photon Energies of 510pm10MeV (ThetaCM30-60)", 10, -180, 180);
  Phin530CM2 = new GH1("Phin_530MeVCM2", "Neutron Phi Distribution for Photon Energies of 530pm10MeV (ThetaCM30-60)", 10, -180, 180);
  Phin550CM2 = new GH1("Phin_550MeVCM2", "Neutron Phi Distribution for Photon Energies of 550pm10MeV (ThetaCM30-60)", 10, -180, 180);
  Phin570CM2 = new GH1("Phin_570MeVCM2", "Neutron Phi Distribution for Photon Energies of 570pm10MeV (ThetaCM30-60)", 10, -180, 180);
  Phin590CM2 = new GH1("Phin_590MeVCM2", "Neutron Phi Distribution for Photon Energies of 590pm10MeV (ThetaCM30-60)", 10, -180, 180);
  Phin610CM2 = new GH1("Phin_610MeVCM2", "Neutron Phi Distribution for Photon Energies of 610pm10MeV (ThetaCM30-60)", 10, -180, 180);
  Phin630CM2 = new GH1("Phin_630MeVCM2", "Neutron Phi Distribution for Photon Energies of 630pm10MeV (ThetaCM30-60)", 10, -180, 180);

  // Proton Phi dists across EGamma bins
  Phip410CM3 = new GH1("Phip_410MeVCM3", "Proton Phi Distribution for Photon Energies of 410pm10MeV (ThetaCM60-90)", 10, -180, 180);
  Phip430CM3 = new GH1("Phip_430MeVCM3", "Proton Phi Distribution for Photon Energies of 430pm10MeV (ThetaCM60-90)", 10, -180, 180);
  Phip450CM3 = new GH1("Phip_450MeVCM3", "Proton Phi Distribution for Photon Energies of 450pm10MeV (ThetaCM60-90)", 10, -180, 180);
  Phip470CM3 = new GH1("Phip_470MeVCM3", "Proton Phi Distribution for Photon Energies of 470pm10MeV (ThetaCM60-90)", 10, -180, 180);
  Phip490CM3 = new GH1("Phip_490MeVCM3", "Proton Phi Distribution for Photon Energies of 490pm10MeV (ThetaCM60-90)", 10, -180, 180);
  Phip510CM3 = new GH1("Phip_510MeVCM3", "Proton Phi Distribution for Photon Energies of 510pm10MeV (ThetaCM60-90)", 10, -180, 180);
  Phip530CM3 = new GH1("Phip_530MeVCM3", "Proton Phi Distribution for Photon Energies of 530pm10MeV (ThetaCM60-90)", 10, -180, 180);
  Phip550CM3 = new GH1("Phip_550MeVCM3", "Proton Phi Distribution for Photon Energies of 550pm10MeV (ThetaCM60-90)", 10, -180, 180);
  Phip570CM3 = new GH1("Phip_570MeVCM3", "Proton Phi Distribution for Photon Energies of 570pm10MeV (ThetaCM60-90)", 10, -180, 180);
  Phip590CM3 = new GH1("Phip_590MeVCM3", "Proton Phi Distribution for Photon Energies of 590pm10MeV (ThetaCM60-90)", 10, -180, 180);
  Phip610CM3 = new GH1("Phip_610MeVCM3", "Proton Phi Distribution for Photon Energies of 610pm10MeV (ThetaCM60-90)", 10, -180, 180);
  Phip630CM3 = new GH1("Phip_630MeVCM3", "Proton Phi Distribution for Photon Energies of 630pm10MeV (ThetaCM60-90)", 10, -180, 180);

  // Angles of neutron in scattered frame across EGamma bins
  Phin410CM3 = new GH1("Phin_410MeVCM3", "Neutron Phi Distribution for Photon Energies of 410pm10MeV (ThetaCM60-90)", 10, -180, 180);
  Phin430CM3 = new GH1("Phin_430MeVCM3", "Neutron Phi Distribution for Photon Energies of 430pm10MeV (ThetaCM60-90)", 10, -180, 180);
  Phin450CM3 = new GH1("Phin_450MeVCM3", "Neutron Phi Distribution for Photon Energies of 450pm10MeV (ThetaCM60-90)", 10, -180, 180);
  Phin470CM3 = new GH1("Phin_470MeVCM3", "Neutron Phi Distribution for Photon Energies of 470pm10MeV (ThetaCM60-90)", 10, -180, 180);
  Phin490CM3 = new GH1("Phin_490MeVCM3", "Neutron Phi Distribution for Photon Energies of 490pm10MeV (ThetaCM60-90)", 10, -180, 180);
  Phin510CM3 = new GH1("Phin_510MeVCM3", "Neutron Phi Distribution for Photon Energies of 510pm10MeV (ThetaCM60-90)", 10, -180, 180);
  Phin530CM3 = new GH1("Phin_530MeVCM3", "Neutron Phi Distribution for Photon Energies of 530pm10MeV (ThetaCM60-90)", 10, -180, 180);
  Phin550CM3 = new GH1("Phin_550MeVCM3", "Neutron Phi Distribution for Photon Energies of 550pm10MeV (ThetaCM60-90)", 10, -180, 180);
  Phin570CM3 = new GH1("Phin_570MeVCM3", "Neutron Phi Distribution for Photon Energies of 570pm10MeV (ThetaCM60-90)", 10, -180, 180);
  Phin590CM3 = new GH1("Phin_590MeVCM3", "Neutron Phi Distribution for Photon Energies of 590pm10MeV (ThetaCM60-90)", 10, -180, 180);
  Phin610CM3 = new GH1("Phin_610MeVCM3", "Neutron Phi Distribution for Photon Energies of 610pm10MeV (ThetaCM60-90)", 10, -180, 180);
  Phin630CM3 = new GH1("Phin_630MeVCM3", "Neutron Phi Distribution for Photon Energies of 630pm10MeV (ThetaCM60-90)", 10, -180, 180);

  // Proton Phi dists across EGamma bins
  Phip410CM4 = new GH1("Phip_410MeVCM4", "Proton Phi Distribution for Photon Energies of 410pm10MeV (ThetaCM90-120)", 10, -180, 180);
  Phip430CM4 = new GH1("Phip_430MeVCM4", "Proton Phi Distribution for Photon Energies of 430pm10MeV (ThetaCM90-120)", 10, -180, 180);
  Phip450CM4 = new GH1("Phip_450MeVCM4", "Proton Phi Distribution for Photon Energies of 450pm10MeV (ThetaCM90-120)", 10, -180, 180);
  Phip470CM4 = new GH1("Phip_470MeVCM4", "Proton Phi Distribution for Photon Energies of 470pm10MeV (ThetaCM90-120)", 10, -180, 180);
  Phip490CM4 = new GH1("Phip_490MeVCM4", "Proton Phi Distribution for Photon Energies of 490pm10MeV (ThetaCM90-120)", 10, -180, 180);
  Phip510CM4 = new GH1("Phip_510MeVCM4", "Proton Phi Distribution for Photon Energies of 510pm10MeV (ThetaCM90-120)", 10, -180, 180);
  Phip530CM4 = new GH1("Phip_530MeVCM4", "Proton Phi Distribution for Photon Energies of 530pm10MeV (ThetaCM90-120)", 10, -180, 180);
  Phip550CM4 = new GH1("Phip_550MeVCM4", "Proton Phi Distribution for Photon Energies of 550pm10MeV (ThetaCM90-120)", 10, -180, 180);
  Phip570CM4 = new GH1("Phip_570MeVCM4", "Proton Phi Distribution for Photon Energies of 570pm10MeV (ThetaCM90-120)", 10, -180, 180);
  Phip590CM4 = new GH1("Phip_590MeVCM4", "Proton Phi Distribution for Photon Energies of 590pm10MeV (ThetaCM90-120)", 10, -180, 180);
  Phip610CM4 = new GH1("Phip_610MeVCM4", "Proton Phi Distribution for Photon Energies of 610pm10MeV (ThetaCM90-120)", 10, -180, 180);
  Phip630CM4 = new GH1("Phip_630MeVCM4", "Proton Phi Distribution for Photon Energies of 630pm10MeV (ThetaCM90-120)", 10, -180, 180);

  // Angles of neutron in scattered frame across EGamma bins
  Phin410CM4 = new GH1("Phin_410MeVCM4", "Neutron Phi Distribution for Photon Energies of 410pm10MeV (ThetaCM90-120)", 10, -180, 180);
  Phin430CM4 = new GH1("Phin_430MeVCM4", "Neutron Phi Distribution for Photon Energies of 430pm10MeV (ThetaCM90-120)", 10, -180, 180);
  Phin450CM4 = new GH1("Phin_450MeVCM4", "Neutron Phi Distribution for Photon Energies of 450pm10MeV (ThetaCM90-120)", 10, -180, 180);
  Phin470CM4 = new GH1("Phin_470MeVCM4", "Neutron Phi Distribution for Photon Energies of 470pm10MeV (ThetaCM90-120)", 10, -180, 180);
  Phin490CM4 = new GH1("Phin_490MeVCM4", "Neutron Phi Distribution for Photon Energies of 490pm10MeV (ThetaCM90-120)", 10, -180, 180);
  Phin510CM4 = new GH1("Phin_510MeVCM4", "Neutron Phi Distribution for Photon Energies of 510pm10MeV (ThetaCM90-120)", 10, -180, 180);
  Phin530CM4 = new GH1("Phin_530MeVCM4", "Neutron Phi Distribution for Photon Energies of 530pm10MeV (ThetaCM90-120)", 10, -180, 180);
  Phin550CM4 = new GH1("Phin_550MeVCM4", "Neutron Phi Distribution for Photon Energies of 550pm10MeV (ThetaCM90-120)", 10, -180, 180);
  Phin570CM4 = new GH1("Phin_570MeVCM4", "Neutron Phi Distribution for Photon Energies of 570pm10MeV (ThetaCM90-120)", 10, -180, 180);
  Phin590CM4 = new GH1("Phin_590MeVCM4", "Neutron Phi Distribution for Photon Energies of 590pm10MeV (ThetaCM90-120)", 10, -180, 180);
  Phin610CM4 = new GH1("Phin_610MeVCM4", "Neutron Phi Distribution for Photon Energies of 610pm10MeV (ThetaCM90-120)", 10, -180, 180);
  Phin630CM4 = new GH1("Phin_630MeVCM4", "Neutron Phi Distribution for Photon Energies of 630pm10MeV (ThetaCM90-120)", 10, -180, 180);

  // Proton Phi dists across EGamma bins
  Phip410CM5 = new GH1("Phip_410MeVCM5", "Proton Phi Distribution for Photon Energies of 410pm10MeV (ThetaCM120-150)", 10, -180, 180);
  Phip430CM5 = new GH1("Phip_430MeVCM5", "Proton Phi Distribution for Photon Energies of 430pm10MeV (ThetaCM120-150)", 10, -180, 180);
  Phip450CM5 = new GH1("Phip_450MeVCM5", "Proton Phi Distribution for Photon Energies of 450pm10MeV (ThetaCM120-150)", 10, -180, 180);
  Phip470CM5 = new GH1("Phip_470MeVCM5", "Proton Phi Distribution for Photon Energies of 470pm10MeV (ThetaCM120-150)", 10, -180, 180);
  Phip490CM5 = new GH1("Phip_490MeVCM5", "Proton Phi Distribution for Photon Energies of 490pm10MeV (ThetaCM120-150)", 10, -180, 180);
  Phip510CM5 = new GH1("Phip_510MeVCM5", "Proton Phi Distribution for Photon Energies of 510pm10MeV (ThetaCM120-150)", 10, -180, 180);
  Phip530CM5 = new GH1("Phip_530MeVCM5", "Proton Phi Distribution for Photon Energies of 530pm10MeV (ThetaCM120-150)", 10, -180, 180);
  Phip550CM5 = new GH1("Phip_550MeVCM5", "Proton Phi Distribution for Photon Energies of 550pm10MeV (ThetaCM120-150)", 10, -180, 180);
  Phip570CM5 = new GH1("Phip_570MeVCM5", "Proton Phi Distribution for Photon Energies of 570pm10MeV (ThetaCM120-150)", 10, -180, 180);
  Phip590CM5 = new GH1("Phip_590MeVCM5", "Proton Phi Distribution for Photon Energies of 590pm10MeV (ThetaCM120-150)", 10, -180, 180);
  Phip610CM5 = new GH1("Phip_610MeVCM5", "Proton Phi Distribution for Photon Energies of 610pm10MeV (ThetaCM120-150)", 10, -180, 180);
  Phip630CM5 = new GH1("Phip_630MeVCM5", "Proton Phi Distribution for Photon Energies of 630pm10MeV (ThetaCM120-150)", 10, -180, 180);

  // Angles of neutron in scattered frame across EGamma bins
  Phin410CM5 = new GH1("Phin_410MeVCM5", "Neutron Phi Distribution for Photon Energies of 410pm10MeV (ThetaCM120-150)", 10, -180, 180);
  Phin430CM5 = new GH1("Phin_430MeVCM5", "Neutron Phi Distribution for Photon Energies of 430pm10MeV (ThetaCM120-150)", 10, -180, 180);
  Phin450CM5 = new GH1("Phin_450MeVCM5", "Neutron Phi Distribution for Photon Energies of 450pm10MeV (ThetaCM120-150)", 10, -180, 180);
  Phin470CM5 = new GH1("Phin_470MeVCM5", "Neutron Phi Distribution for Photon Energies of 470pm10MeV (ThetaCM120-150)", 10, -180, 180);
  Phin490CM5 = new GH1("Phin_490MeVCM5", "Neutron Phi Distribution for Photon Energies of 490pm10MeV (ThetaCM120-150)", 10, -180, 180);
  Phin510CM5 = new GH1("Phin_510MeVCM5", "Neutron Phi Distribution for Photon Energies of 510pm10MeV (ThetaCM120-150)", 10, -180, 180);
  Phin530CM5 = new GH1("Phin_530MeVCM5", "Neutron Phi Distribution for Photon Energies of 530pm10MeV (ThetaCM120-150)", 10, -180, 180);
  Phin550CM5 = new GH1("Phin_550MeVCM5", "Neutron Phi Distribution for Photon Energies of 550pm10MeV (ThetaCM120-150)", 10, -180, 180);
  Phin570CM5 = new GH1("Phin_570MeVCM5", "Neutron Phi Distribution for Photon Energies of 570pm10MeV (ThetaCM120-150)", 10, -180, 180);
  Phin590CM5 = new GH1("Phin_590MeVCM5", "Neutron Phi Distribution for Photon Energies of 590pm10MeV (ThetaCM120-150)", 10, -180, 180);
  Phin610CM5 = new GH1("Phin_610MeVCM5", "Neutron Phi Distribution for Photon Energies of 610pm10MeV (ThetaCM120-150)", 10, -180, 180);
  Phin630CM5 = new GH1("Phin_630MeVCM5", "Neutron Phi Distribution for Photon Energies of 630pm10MeV (ThetaCM120-150)", 10, -180, 180);

  // Proton Phi dists across EGamma bins
  Phip410CM6 = new GH1("Phip_410MeVCM6", "Proton Phi Distribution for Photon Energies of 410pm10MeV (ThetaCM150-180)", 10, -180, 180);
  Phip430CM6 = new GH1("Phip_430MeVCM6", "Proton Phi Distribution for Photon Energies of 430pm10MeV (ThetaCM150-180)", 10, -180, 180);
  Phip450CM6 = new GH1("Phip_450MeVCM6", "Proton Phi Distribution for Photon Energies of 450pm10MeV (ThetaCM150-180)", 10, -180, 180);
  Phip470CM6 = new GH1("Phip_470MeVCM6", "Proton Phi Distribution for Photon Energies of 470pm10MeV (ThetaCM150-180)", 10, -180, 180);
  Phip490CM6 = new GH1("Phip_490MeVCM6", "Proton Phi Distribution for Photon Energies of 490pm10MeV (ThetaCM150-180)", 10, -180, 180);
  Phip510CM6 = new GH1("Phip_510MeVCM6", "Proton Phi Distribution for Photon Energies of 510pm10MeV (ThetaCM150-180)", 10, -180, 180);
  Phip530CM6 = new GH1("Phip_530MeVCM6", "Proton Phi Distribution for Photon Energies of 530pm10MeV (ThetaCM150-180)", 10, -180, 180);
  Phip550CM6 = new GH1("Phip_550MeVCM6", "Proton Phi Distribution for Photon Energies of 550pm10MeV (ThetaCM150-180)", 10, -180, 180);
  Phip570CM6 = new GH1("Phip_570MeVCM6", "Proton Phi Distribution for Photon Energies of 570pm10MeV (ThetaCM150-180)", 10, -180, 180);
  Phip590CM6 = new GH1("Phip_590MeVCM6", "Proton Phi Distribution for Photon Energies of 590pm10MeV (ThetaCM150-180)", 10, -180, 180);
  Phip610CM6 = new GH1("Phip_610MeVCM6", "Proton Phi Distribution for Photon Energies of 610pm10MeV (ThetaCM150-180)", 10, -180, 180);
  Phip630CM6 = new GH1("Phip_630MeVCM6", "Proton Phi Distribution for Photon Energies of 630pm10MeV (ThetaCM150-180)", 10, -180, 180);

  // Angles of neutron in scattered frame across EGamma bins
  Phin410CM6 = new GH1("Phin_410MeVCM6", "Neutron Phi Distribution for Photon Energies of 410pm10MeV (ThetaCM150-180)", 10, -180, 180);
  Phin430CM6 = new GH1("Phin_430MeVCM6", "Neutron Phi Distribution for Photon Energies of 430pm10MeV (ThetaCM150-180)", 10, -180, 180);
  Phin450CM6 = new GH1("Phin_450MeVCM6", "Neutron Phi Distribution for Photon Energies of 450pm10MeV (ThetaCM150-180)", 10, -180, 180);
  Phin470CM6 = new GH1("Phin_470MeVCM6", "Neutron Phi Distribution for Photon Energies of 470pm10MeV (ThetaCM150-180)", 10, -180, 180);
  Phin490CM6 = new GH1("Phin_490MeVCM6", "Neutron Phi Distribution for Photon Energies of 490pm10MeV (ThetaCM150-180)", 10, -180, 180);
  Phin510CM6 = new GH1("Phin_510MeVCM6", "Neutron Phi Distribution for Photon Energies of 510pm10MeV (ThetaCM150-180)", 10, -180, 180);
  Phin530CM6 = new GH1("Phin_530MeVCM6", "Neutron Phi Distribution for Photon Energies of 530pm10MeV (ThetaCM150-180)", 10, -180, 180);
  Phin550CM6 = new GH1("Phin_550MeVCM6", "Neutron Phi Distribution for Photon Energies of 550pm10MeV (ThetaCM150-180)", 10, -180, 180);
  Phin570CM6 = new GH1("Phin_570MeVCM6", "Neutron Phi Distribution for Photon Energies of 570pm10MeV (ThetaCM150-180)", 10, -180, 180);
  Phin590CM6 = new GH1("Phin_590MeVCM6", "Neutron Phi Distribution for Photon Energies of 590pm10MeV (ThetaCM150-180)", 10, -180, 180);
  Phin610CM6 = new GH1("Phin_610MeVCM6", "Neutron Phi Distribution for Photon Energies of 610pm10MeV (ThetaCM150-180)", 10, -180, 180);
  Phin630CM6 = new GH1("Phin_630MeVCM6", "Neutron Phi Distribution for Photon Energies of 630pm10MeV (ThetaCM150-180)", 10, -180, 180);

  ThetaRecPiDiff = new GH1 ("ThetaRecPiDiff", "Difference between ThetaPiRec and Thetan", 200, 0, 180);
  ThetanThetaRecPi = new GH2 ("ThetanThetaRecPi", "Thetan vs ThetaPiRec", 100, 0, 180, 100, 0, 180);
  ThetanThetaRecPiDiff = new GH2 ("ThetanThetaRecPiDiff", "Thetan vs (ThetaPiRec - Thetan)", 100, 0, 180, 100, 0, 180);
  ThetaRecPDiff = new GH1 ("ThetaRecPDiff", "Difference between ThetaPRec and Thetan", 200, 0, 180);
  ThetanThetaRecP = new GH2 ("ThetanThetaRecP", "Thetan vs ThetaPRec", 100, 0, 180, 100, 0, 180);
  ThetanThetaRecPDiff = new GH2 ("ThetanThetaRecPDiff", "Thetan vs (ThetaPRec - Thetan)", 100, 0, 180, 100, 0, 180);

  E_dE = new GH2 ("E_dE", "EdE Plot With E Loss Adjustment", 100, 0, 500, 100, 0, 5);
  E_dE_KinCut = new GH2 ("E_dE_KinCut", "EdE Plot (With cut on kinematic proton banana + E Loss)", 100, 0, 500, 100, 0, 5);
  KinEp_dE = new GH2 ("KinEp_dE", "KinEpdE Plot", 100, 0, 500, 100, 0, 5);
  //KinEp_dE_GoodCut = new GH2 ("KinEp_dE_GoodCut", "KinEpdE Plot With Good Proton Cut", 100, 0, 500, 100, 0, 5);
  ThetaScPhiSc = new GH2 ("ThetaScPhiSc", "Phi as a function of Theta (Both in rotated frame)", 100, 0, 180, 100, -180, 180);
  E_KinEp = new GH2 ("E_KinEp", "Kinematic Energy of Proton as a function of CB energy", 100, 0, 500, 100, 0, 500);
  PhinDiffWCZRec = new GH2 ("PhinDiffWCZRec", "Difference between WC Phi and Reconstructed Phi as a fn of WCZ Hit Position", 100, 0, 200, 100, 0, 180);
  PhinDiffWCZRec_KinCut = new GH2 ("PhinDiffWCZRec_KinCut", "Difference between WC Phi and Reconstructed Phi as a fn of WCZ Hit Position (Kin P Banana Cut)", 200, -300, 300, 200, 0, 180);

}

void PNeutPol_Polarimeter_Lin::FillHists()
{
    time->Fill(TaggerTime);
    if (-5 < TaggerTime && TaggerTime < 20) time_cut->Fill(TaggerTime);

    Eg->Fill(EGamma, TaggerTime);
    E_dE->Fill(EpCorr, dEp, TaggerTime);
    KinEp_dE->Fill(KinEp, dEp, TaggerTime);
    EpKin->Fill(KinEp, TaggerTime);
    EpCorrected->Fill(EpCorr, TaggerTime);
    EpKinEpCorrDiff->Fill(KinEDiff, TaggerTime);
    EpEpCorrDiff->Fill(EpDiff, TaggerTime);
    MMpEpCorrected->Fill(MMpEpCorr, TaggerTime);
    OAngle->Fill(OpeningAngle, TaggerTime);
    WCZnRecon->Fill(WCZnRec, TaggerTime);
    ZpDist->Fill(Zp, TaggerTime);

    WCPhiDifference->Fill(PhiWCDiff);
    E_KinEp->Fill(EpCorr, KinEp, TaggerTime);
    PhinDiffWCZRec->Fill(WCZnRec, PhinDiff, TaggerTime);
    ThetaRecPiDiff->Fill(ThetaPiRecDiff, TaggerTime);
    ThetanThetaRecPi->Fill(Thetan, ThetaPiRec, TaggerTime);
    ThetanThetaRecPiDiff->Fill(Thetan, ThetaPiRecDiff, TaggerTime);
    ThetaRecPDiff->Fill(ThetapRecDiff, TaggerTime);
    ThetanThetaRecP->Fill(Thetan, ThetapRec, TaggerTime);
    ThetanThetaRecPDiff->Fill(Thetan, ThetapRecDiff, TaggerTime);

    // Fill events inside good proton banana on KinEpdE plot
    if(Cut_protonKinGood -> IsInside(KinEp, dEp) == kTRUE)
    {
        //KinEp_dE_GoodCut->Fill(KinEp, dEp, TaggerTime);
        MMpEpCorrectedCut->Fill(MMpEpCorr, TaggerTime);
        EgCut->Fill(EGamma, TaggerTime);
        OAngleCut->Fill(OpeningAngle, TaggerTime);
        ThetaSc -> Fill(ScattTheta, TaggerTime);
        PhiSc -> Fill(ScattPhi, TaggerTime);
        ThetaScPhiSc->Fill(ScattTheta, ScattPhi, TaggerTime);
        E_dE_KinCut->Fill(EpCorr, dEp, TaggerTime);
        PhinDiffWCZRec_KinCut->Fill(WCZnRec, PhinDiff, TaggerTime);

        if(ScattPhi < -165){
            ZpPhiScatNeg180->Fill(Zp, TaggerTime);
        }

        if(ScattPhi < 15 && ScattPhi > -15){
            ZpPhiScat0->Fill(Zp, TaggerTime);
        }

        if(ScattPhi > 165){
            ZpPhiScatPos180->Fill(Zp, TaggerTime);
        }

        if(200 < EGamma && EGamma < 300){
            MMp200300->Fill(MMpEpCorr, TaggerTime);
            OAngleCut200400->Fill(OpeningAngle, TaggerTime);
        }

        else if(300 < EGamma && EGamma < 400){
            MMp300400->Fill(MMpEpCorr, TaggerTime);
            OAngleCut200400->Fill(OpeningAngle, TaggerTime);
        }

        else if(400 < EGamma && EGamma < 500){
            MMp400500->Fill(MMpEpCorr, TaggerTime);
        }

        else if(500 < EGamma && EGamma < 600){
            MMp500600->Fill(MMpEpCorr, TaggerTime);
        }

        else if(600 < EGamma && EGamma < 700){
            MMp600700->Fill(MMpEpCorr, TaggerTime);
        }

        else if(700 < EGamma && EGamma < 800){
            MMp700800->Fill(MMpEpCorr, TaggerTime);
        }

        else if(800 < EGamma && EGamma < 900){
            MMp800900->Fill(MMpEpCorr, TaggerTime);
        }

        if ( 400 < EGamma && EGamma < 420) {
            PhiSc410->Fill(ScattPhi, TaggerTime);

            if(0 < ThetapCM && ThetapCM < 30){
                Phip410CM1->Fill(WCPhip, TaggerTime);
                Phin410CM1->Fill(PhinRec, TaggerTime);
                }

            else if(30 < ThetapCM && ThetapCM < 60){
                Phip410CM2->Fill(WCPhip, TaggerTime);
                Phin410CM2->Fill(PhinRec, TaggerTime);
                }

            else if(60 < ThetapCM && ThetapCM < 90){
                Phip410CM3->Fill(WCPhip, TaggerTime);
                Phin410CM3->Fill(PhinRec, TaggerTime);
                }

            else if(90 < ThetapCM && ThetapCM < 120){
                Phip410CM4->Fill(WCPhip, TaggerTime);
                Phin410CM4->Fill(PhinRec, TaggerTime);
                }

            else if(120 < ThetapCM && ThetapCM < 150){
                Phip410CM5->Fill(WCPhip, TaggerTime);
                Phin410CM5->Fill(PhinRec, TaggerTime);
                }

            else if(150 < ThetapCM && ThetapCM < 180){
                Phip410CM6->Fill(WCPhip, TaggerTime);
                Phin410CM6->Fill(PhinRec, TaggerTime);
                }
        }

        else if ( 420 < EGamma && EGamma < 440) {
            PhiSc430->Fill(ScattPhi, TaggerTime);

            if(0 < ThetapCM && ThetapCM < 30){
                Phip430CM1->Fill(WCPhip, TaggerTime);
                Phin430CM1->Fill(PhinRec, TaggerTime);
                }

            else if(30 < ThetapCM && ThetapCM < 60){
                Phip430CM2->Fill(WCPhip, TaggerTime);
                Phin430CM2->Fill(PhinRec, TaggerTime);
                }

            else if(60 < ThetapCM && ThetapCM < 90){
                Phip430CM3->Fill(WCPhip, TaggerTime);
                Phin430CM3->Fill(PhinRec, TaggerTime);
                }

            else if(90 < ThetapCM && ThetapCM < 120){
                Phip430CM4->Fill(WCPhip, TaggerTime);
                Phin430CM4->Fill(PhinRec, TaggerTime);
                }

            else if(120 < ThetapCM && ThetapCM < 150){
                Phip430CM5->Fill(WCPhip, TaggerTime);
                Phin430CM5->Fill(PhinRec, TaggerTime);
                }

            else if(150 < ThetapCM && ThetapCM < 180){
                Phip430CM6->Fill(WCPhip, TaggerTime);
                Phin430CM6->Fill(PhinRec, TaggerTime);
                }
        }

        else if ( 440 < EGamma && EGamma < 460) {
            PhiSc450->Fill(ScattPhi, TaggerTime);

            if(0 < ThetapCM && ThetapCM < 30){
                Phip450CM1->Fill(WCPhip, TaggerTime);
                Phin450CM1->Fill(PhinRec, TaggerTime);
                }

            else if(30 < ThetapCM && ThetapCM < 60){
                Phip450CM2->Fill(WCPhip, TaggerTime);
                Phin450CM2->Fill(PhinRec, TaggerTime);
                }

            else if(60 < ThetapCM && ThetapCM < 90){
                Phip450CM3->Fill(WCPhip, TaggerTime);
                Phin450CM3->Fill(PhinRec, TaggerTime);
                }

            else if(90 < ThetapCM && ThetapCM < 120){
                Phip450CM4->Fill(WCPhip, TaggerTime);
                Phin450CM4->Fill(PhinRec, TaggerTime);
                }

            else if(120 < ThetapCM && ThetapCM < 150){
                Phip450CM5->Fill(WCPhip, TaggerTime);
                Phin450CM5->Fill(PhinRec, TaggerTime);
                }

            else if(150 < ThetapCM && ThetapCM < 180){
                Phip450CM6->Fill(WCPhip, TaggerTime);
                Phin450CM6->Fill(PhinRec, TaggerTime);
                }
        }

        else if ( 460 < EGamma && EGamma < 480) {
            PhiSc470->Fill(ScattPhi, TaggerTime);

            if(0 < ThetapCM && ThetapCM < 30){
                Phip470CM1->Fill(WCPhip, TaggerTime);
                Phin470CM1->Fill(PhinRec, TaggerTime);
                }

            else if(30 < ThetapCM && ThetapCM < 60){
                Phip470CM2->Fill(WCPhip, TaggerTime);
                Phin470CM2->Fill(PhinRec, TaggerTime);
                }

            else if(60 < ThetapCM && ThetapCM < 90){
                Phip470CM3->Fill(WCPhip, TaggerTime);
                Phin470CM3->Fill(PhinRec, TaggerTime);
                }

            else if(90 < ThetapCM && ThetapCM < 120){
                Phip470CM4->Fill(WCPhip, TaggerTime);
                Phin470CM4->Fill(PhinRec, TaggerTime);
                }

            else if(120 < ThetapCM && ThetapCM < 150){
                Phip470CM5->Fill(WCPhip, TaggerTime);
                Phin470CM5->Fill(PhinRec, TaggerTime);
                }

            else if(150 < ThetapCM && ThetapCM < 180){
                Phip470CM6->Fill(WCPhip, TaggerTime);
                Phin470CM6->Fill(PhinRec, TaggerTime);
                }
        }

        else if ( 480 < EGamma && EGamma < 500) {
            PhiSc490->Fill(ScattPhi, TaggerTime);

            if(0 < ThetapCM && ThetapCM < 30){
                Phip490CM1->Fill(WCPhip, TaggerTime);
                Phin490CM1->Fill(PhinRec, TaggerTime);
                }

            else if(30 < ThetapCM && ThetapCM < 60){
                Phip490CM2->Fill(WCPhip, TaggerTime);
                Phin490CM2->Fill(PhinRec, TaggerTime);
                }

            else if(60 < ThetapCM && ThetapCM < 90){
                Phip490CM3->Fill(WCPhip, TaggerTime);
                Phin490CM3->Fill(PhinRec, TaggerTime);
                }

            else if(90 < ThetapCM && ThetapCM < 120){
                Phip490CM4->Fill(WCPhip, TaggerTime);
                Phin490CM4->Fill(PhinRec, TaggerTime);
                }

            else if(120 < ThetapCM && ThetapCM < 150){
                Phip490CM5->Fill(WCPhip, TaggerTime);
                Phin490CM5->Fill(PhinRec, TaggerTime);
                }

            else if(150 < ThetapCM && ThetapCM < 180){
                Phip490CM6->Fill(WCPhip, TaggerTime);
                Phin490CM6->Fill(PhinRec, TaggerTime);
                }
        }

        else if ( 500 < EGamma && EGamma < 520) {
            PhiSc510->Fill(ScattPhi, TaggerTime);

            if(0 < ThetapCM && ThetapCM < 30){
                Phip510CM1->Fill(WCPhip, TaggerTime);
                Phin510CM1->Fill(PhinRec, TaggerTime);
                }

            else if(30 < ThetapCM && ThetapCM < 60){
                Phip510CM2->Fill(WCPhip, TaggerTime);
                Phin510CM2->Fill(PhinRec, TaggerTime);
                }

            else if(60 < ThetapCM && ThetapCM < 90){
                Phip510CM3->Fill(WCPhip, TaggerTime);
                Phin510CM3->Fill(PhinRec, TaggerTime);
                }

            else if(90 < ThetapCM && ThetapCM < 120){
                Phip510CM4->Fill(WCPhip, TaggerTime);
                Phin510CM4->Fill(PhinRec, TaggerTime);
                }

            else if(120 < ThetapCM && ThetapCM < 150){
                Phip510CM5->Fill(WCPhip, TaggerTime);
                Phin510CM5->Fill(PhinRec, TaggerTime);
                }

            else if(150 < ThetapCM && ThetapCM < 180){
                Phip510CM6->Fill(WCPhip, TaggerTime);
                Phin510CM6->Fill(PhinRec, TaggerTime);
                }
        }

        else if ( 520 < EGamma && EGamma < 540) {
            PhiSc530->Fill(ScattPhi, TaggerTime);

            if(0 < ThetapCM && ThetapCM < 30){
                Phip530CM1->Fill(WCPhip, TaggerTime);
                Phin530CM1->Fill(PhinRec, TaggerTime);
                }

            else if(30 < ThetapCM && ThetapCM < 60){
                Phip530CM2->Fill(WCPhip, TaggerTime);
                Phin530CM2->Fill(PhinRec, TaggerTime);
                }

            else if(60 < ThetapCM && ThetapCM < 90){
                Phip530CM3->Fill(WCPhip, TaggerTime);
                Phin530CM3->Fill(PhinRec, TaggerTime);
                }

            else if(90 < ThetapCM && ThetapCM < 120){
                Phip530CM4->Fill(WCPhip, TaggerTime);
                Phin530CM4->Fill(PhinRec, TaggerTime);
                }

            else if(120 < ThetapCM && ThetapCM < 150){
                Phip530CM5->Fill(WCPhip, TaggerTime);
                Phin530CM5->Fill(PhinRec, TaggerTime);
                }

            else if(150 < ThetapCM && ThetapCM < 180){
                Phip530CM6->Fill(WCPhip, TaggerTime);
                Phin530CM6->Fill(PhinRec, TaggerTime);
                }
        }

        else if ( 540 < EGamma && EGamma < 560) {
            PhiSc550->Fill(ScattPhi, TaggerTime);

            if(0 < ThetapCM && ThetapCM < 30){
                Phip550CM1->Fill(WCPhip, TaggerTime);
                Phin550CM1->Fill(PhinRec, TaggerTime);
                }

            else if(30 < ThetapCM && ThetapCM < 60){
                Phip550CM2->Fill(WCPhip, TaggerTime);
                Phin550CM2->Fill(PhinRec, TaggerTime);
                }

            else if(60 < ThetapCM && ThetapCM < 90){
                Phip550CM3->Fill(WCPhip, TaggerTime);
                Phin550CM3->Fill(PhinRec, TaggerTime);
                }

            else if(90 < ThetapCM && ThetapCM < 120){
                Phip550CM4->Fill(WCPhip, TaggerTime);
                Phin550CM4->Fill(PhinRec, TaggerTime);
                }

            else if(120 < ThetapCM && ThetapCM < 150){
                Phip550CM5->Fill(WCPhip, TaggerTime);
                Phin550CM5->Fill(PhinRec, TaggerTime);
                }

            else if(150 < ThetapCM && ThetapCM < 180){
                Phip550CM6->Fill(WCPhip, TaggerTime);
                Phin550CM6->Fill(PhinRec, TaggerTime);
                }
        }

        else if ( 560 < EGamma && EGamma < 580) {
            PhiSc570->Fill(ScattPhi, TaggerTime);

            if(0 < ThetapCM && ThetapCM < 30){
                Phip570CM1->Fill(WCPhip, TaggerTime);
                Phin570CM1->Fill(PhinRec, TaggerTime);
                }

            else if(30 < ThetapCM && ThetapCM < 60){
                Phip570CM2->Fill(WCPhip, TaggerTime);
                Phin570CM2->Fill(PhinRec, TaggerTime);
                }

            else if(60 < ThetapCM && ThetapCM < 90){
                Phip570CM3->Fill(WCPhip, TaggerTime);
                Phin570CM3->Fill(PhinRec, TaggerTime);
                }

            else if(90 < ThetapCM && ThetapCM < 120){
                Phip570CM4->Fill(WCPhip, TaggerTime);
                Phin570CM4->Fill(PhinRec, TaggerTime);
                }

            else if(120 < ThetapCM && ThetapCM < 150){
                Phip570CM5->Fill(WCPhip, TaggerTime);
                Phin570CM5->Fill(PhinRec, TaggerTime);
                }

            else if(150 < ThetapCM && ThetapCM < 180){
                Phip570CM6->Fill(WCPhip, TaggerTime);
                Phin570CM6->Fill(PhinRec, TaggerTime);
                }
        }

        else if ( 580 < EGamma && EGamma < 600) {
            PhiSc590->Fill(ScattPhi, TaggerTime);

            if(0 < ThetapCM && ThetapCM < 30){
                Phip590CM1->Fill(WCPhip, TaggerTime);
                Phin590CM1->Fill(PhinRec, TaggerTime);
                }

            else if(30 < ThetapCM && ThetapCM < 60){
                Phip590CM2->Fill(WCPhip, TaggerTime);
                Phin590CM2->Fill(PhinRec, TaggerTime);
                }

            else if(60 < ThetapCM && ThetapCM < 90){
                Phip590CM3->Fill(WCPhip, TaggerTime);
                Phin590CM3->Fill(PhinRec, TaggerTime);
                }

            else if(90 < ThetapCM && ThetapCM < 120){
                Phip590CM4->Fill(WCPhip, TaggerTime);
                Phin590CM4->Fill(PhinRec, TaggerTime);
                }

            else if(120 < ThetapCM && ThetapCM < 150){
                Phip590CM5->Fill(WCPhip, TaggerTime);
                Phin590CM5->Fill(PhinRec, TaggerTime);
                }

            else if(150 < ThetapCM && ThetapCM < 180){
                Phip590CM6->Fill(WCPhip, TaggerTime);
                Phin590CM6->Fill(PhinRec, TaggerTime);
                }
        }

        else if ( 600 < EGamma && EGamma < 620) {
            PhiSc610->Fill(ScattPhi, TaggerTime);

            if(0 < ThetapCM && ThetapCM < 30){
                Phip610CM1->Fill(WCPhip, TaggerTime);
                Phin610CM1->Fill(PhinRec, TaggerTime);
                }

            else if(30 < ThetapCM && ThetapCM < 60){
                Phip610CM2->Fill(WCPhip, TaggerTime);
                Phin610CM2->Fill(PhinRec, TaggerTime);
                }

            else if(60 < ThetapCM && ThetapCM < 90){
                Phip610CM3->Fill(WCPhip, TaggerTime);
                Phin610CM3->Fill(PhinRec, TaggerTime);
                }

            else if(90 < ThetapCM && ThetapCM < 120){
                Phip610CM4->Fill(WCPhip, TaggerTime);
                Phin610CM4->Fill(PhinRec, TaggerTime);
                }

            else if(120 < ThetapCM && ThetapCM < 150){
                Phip610CM5->Fill(WCPhip, TaggerTime);
                Phin610CM5->Fill(PhinRec, TaggerTime);
                }

            else if(150 < ThetapCM && ThetapCM < 180){
                Phip610CM6->Fill(WCPhip, TaggerTime);
                Phin610CM6->Fill(PhinRec, TaggerTime);
                }
        }

        else if ( 620 < EGamma && EGamma < 640) {
            PhiSc630->Fill(ScattPhi, TaggerTime);

            if(0 < ThetapCM && ThetapCM < 30){
                Phip630CM1->Fill(WCPhip, TaggerTime);
                Phin630CM1->Fill(PhinRec, TaggerTime);
                }

            else if(30 < ThetapCM && ThetapCM < 60){
                Phip630CM2->Fill(WCPhip, TaggerTime);
                Phin630CM2->Fill(PhinRec, TaggerTime);
                }

            else if(60 < ThetapCM && ThetapCM < 90){
                Phip630CM3->Fill(WCPhip, TaggerTime);
                Phin630CM3->Fill(PhinRec, TaggerTime);
                }

            else if(90 < ThetapCM && ThetapCM < 120){
                Phip630CM4->Fill(WCPhip, TaggerTime);
                Phin630CM4->Fill(PhinRec, TaggerTime);
                }

            else if(120 < ThetapCM && ThetapCM < 150){
                Phip630CM5->Fill(WCPhip, TaggerTime);
                Phin630CM5->Fill(PhinRec, TaggerTime);
                }

            else if(150 < ThetapCM && ThetapCM < 180){
                Phip630CM6->Fill(WCPhip, TaggerTime);
                Phin630CM6->Fill(PhinRec, TaggerTime);
                }
        }
    }
}

Bool_t	PNeutPol_Polarimeter_Lin::Write(){
  // Write all GH1's easily

  GTreeManager::Write();
}
