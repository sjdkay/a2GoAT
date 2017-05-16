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

  Cut_CB_proton = OpenCutFile("configfiles/cuts/CB_DeltaE-E_Proton_11_05_17.root", "Proton"); // These will need adjusting with new Acqu files
  Cut_proton = Cut_CB_proton;
  Cut_CB_pion = OpenCutFile("configfiles/cuts/CB_DeltaE-E_Pion_29_07_15.root", "Pion");
  Cut_pion = Cut_CB_pion;
  Cut_CB_protonKinGood = OpenCutFile("configfiles/cuts/CB_DeltaE-E_ProtonKinGood_11_05_17.root", "ProtonKinGood"); // These will need adjusting with new Acqu files
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
  EventNum = GetEventNumber();

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
    CosThetapCM = cos (ThetapCM * TMath::DegToRad());

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

    if(Cut_proton -> IsInside(EpCorr, dEp) == kFALSE) continue; // If E loss correct proton is NOT inside p banana drop out
    if(Cut_protonKinGood -> IsInside(KinEp, dEp) == kFALSE) continue; // If KinE proton is NOT inside p banana drop out
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

  Eg = new GH1( "Eg", "Photon Energy Distribution", 200, 100, 1600);
  WCPhiDifference = new GH1 ("WCPhiDifference", "WC Phi Difference Between p and n", 180, 0, 360);
  EpKin = new GH1 ("EpKin", "Ep Calculated from Ep/Thetap", 100, 0, 500);
  EpCorrected = new GH1 ("EpCorrected", "Ep Corrected for Energy Loss in Polarimeter ", 100, 0, 500);
  OAngle = new GH1 ("OAngle", "Opening Angle between P and N Vectors", 180, 0, 180);
  WCZnRecon = new GH1 ("WCZnRecon", "WCZ Hit Position from Reconstructed n Vector", 200, 0, 400);

  ThetaSc =  new GH1( "Theta_Scattered", "Scattered Proton Theta Distribution in Rotated Frame", 180, 0, 180);
  PhiSc = new GH1( "Phi_Scattered", "Scattered Proton Phi Distribution in Rotated Frame", 90, -180, 180 );

  EpKinEpCorrDiff = new GH1("EpKinEpCorrDiff", "Difference Between EpKin and EpCorr", 300, -300, 300);
  EpEpCorrDiff = new GH1("EpEpCorrDiff", "Difference Between Ep and EpCorr", 200, 0, 200);
  OAngle200400 = new GH1 ("OAngle200400", "Opening Angle between P and N Vectors (P Banana Cut, 200-400MeV Gamma)", 180, 0, 180);
  MMpEpCorrected = new GH1 ("MMpEpCorrected", "Missing mass seen by Proton (E Loss Corrected)", 400, 0, 2000);

  ZpDist = new GH1 ("ZpDist", "Proton Pseudo Z Vertex Distribution", 200, -400, 400);
  ZpPhiScatNeg180 = new GH1("ZpPhiScatNeg180", "Proton Pseudo Vertex Z for events with PhiSc ~ -ve180", 200, -200, 200);
  ZpPhiScat0 = new GH1("ZpPhiScat0", "Proton Pseudo Vertex Z for events with PhiSc ~ 0", 200, -200, 200);
  ZpPhiScatPos180 = new GH1("ZpPhiScatPos180", "Proton Pseudo Vertex Z for events with PhiSc ~ 180", 200, -200, 200);

  // MMp across photon E bins
  MMp200300 = new GH1("MMp200300", "Missing mass as seen by Proton (200-300MeV E_{#gamma})", 400, 0, 2000);
  MMp300400 = new GH1("MMp300400", "Missing mass as seen by Proton (300-400MeV E_{#gamma})", 400, 0, 2000);
  MMp400500 = new GH1("MMp400500", "Missing mass as seen by Proton (400-500MeV E_{#gamma})", 400, 0, 2000);
  MMp500600 = new GH1("MMp500600", "Missing mass as seen by Proton (500-600MeV E_{#gamma})", 400, 0, 2000);
  MMp600700 = new GH1("MMp600700", "Missing mass as seen by Proton (600-700MeV E_{#gamma})", 400, 0, 2000);
  MMp700800 = new GH1("MMp700800", "Missing mass as seen by Proton (700-800MeV E_{#gamma})", 400, 0, 2000);
  MMp800900 = new GH1("MMp800900", "Missing mass as seen by Proton (800-900MeV E_{#gamma})", 400, 0, 2000);

  // Angles of neutron in scattered frame across EGamma bins
  PhiSc410 = new GH1("Phi_Scattered_410MeV", "#phi_{pScat} Distribution in Rotated Frame for E_{#gamma} 410 #pm 10MeV", 2, -180, 180);
  PhiSc430 = new GH1("Phi_Scattered_430MeV", "#phi_{pScat} Distribution in Rotated Frame for E_{#gamma} 430 #pm 10MeV", 2, -180, 180);
  PhiSc450 = new GH1("Phi_Scattered_450MeV", "#phi_{pScat} Distribution in Rotated Frame for E_{#gamma} 450 #pm 10MeV", 2, -180, 180);
  PhiSc470 = new GH1("Phi_Scattered_470MeV", "#phi_{pScat} Distribution in Rotated Frame for E_{#gamma} 470 #pm 10MeV", 2, -180, 180);
  PhiSc490 = new GH1("Phi_Scattered_490MeV", "#phi_{pScat} Distribution in Rotated Frame for E_{#gamma} 490 #pm 10MeV", 2, -180, 180);
  PhiSc510 = new GH1("Phi_Scattered_510MeV", "#phi_{pScat} Distribution in Rotated Frame for E_{#gamma} 510 #pm 10MeV", 2, -180, 180);
  PhiSc530 = new GH1("Phi_Scattered_530MeV", "#phi_{pScat} Distribution in Rotated Frame for E_{#gamma} 530 #pm 10MeV", 2, -180, 180);
  PhiSc550 = new GH1("Phi_Scattered_550MeV", "#phi_{pScat} Distribution in Rotated Frame for E_{#gamma} 550 #pm 10MeV", 2, -180, 180);
  PhiSc570 = new GH1("Phi_Scattered_570MeV", "#phi_{pScat} Distribution in Rotated Frame for E_{#gamma} 570 #pm 10MeV", 2, -180, 180);
  PhiSc590 = new GH1("Phi_Scattered_590MeV", "#phi_{pScat} Distribution in Rotated Frame for E_{#gamma} 590 #pm 10MeV", 2, -180, 180);
  PhiSc610 = new GH1("Phi_Scattered_610MeV", "#phi_{pScat} Distribution in Rotated Frame for E_{#gamma} 610 #pm 10MeV", 2, -180, 180);
  PhiSc630 = new GH1("Phi_Scattered_630MeV", "#phi_{pScat} Distribution in Rotated Frame for E_{#gamma} 630 #pm 10MeV", 2, -180, 180);

  // #phi_{p} dists across EGamma bins
  Phip410CM1 = new GH1("Phip_410MeVCM1", "#phi_{p} Distribution for E_{#gamma} 410 #pm 10MeV (Cos#theta_{CM}1-0.8)", 10, -180, 180);
  Phip430CM1 = new GH1("Phip_430MeVCM1", "#phi_{p} Distribution for E_{#gamma} 430 #pm 10MeV (Cos#theta_{CM}1-0.8)", 10, -180, 180);
  Phip450CM1 = new GH1("Phip_450MeVCM1", "#phi_{p} Distribution for E_{#gamma} 450 #pm 10MeV (Cos#theta_{CM}1-0.8)", 10, -180, 180);
  Phip470CM1 = new GH1("Phip_470MeVCM1", "#phi_{p} Distribution for E_{#gamma} 470 #pm 10MeV (Cos#theta_{CM}1-0.8)", 10, -180, 180);
  Phip490CM1 = new GH1("Phip_490MeVCM1", "#phi_{p} Distribution for E_{#gamma} 490 #pm 10MeV (Cos#theta_{CM}1-0.8)", 10, -180, 180);
  Phip510CM1 = new GH1("Phip_510MeVCM1", "#phi_{p} Distribution for E_{#gamma} 510 #pm 10MeV (Cos#theta_{CM}1-0.8)", 10, -180, 180);
  Phip530CM1 = new GH1("Phip_530MeVCM1", "#phi_{p} Distribution for E_{#gamma} 530 #pm 10MeV (Cos#theta_{CM}1-0.8)", 10, -180, 180);
  Phip550CM1 = new GH1("Phip_550MeVCM1", "#phi_{p} Distribution for E_{#gamma} 550 #pm 10MeV (Cos#theta_{CM}1-0.8)", 10, -180, 180);
  Phip570CM1 = new GH1("Phip_570MeVCM1", "#phi_{p} Distribution for E_{#gamma} 570 #pm 10MeV (Cos#theta_{CM}1-0.8)", 10, -180, 180);
  Phip590CM1 = new GH1("Phip_590MeVCM1", "#phi_{p} Distribution for E_{#gamma} 590 #pm 10MeV (Cos#theta_{CM}1-0.8)", 10, -180, 180);
  Phip610CM1 = new GH1("Phip_610MeVCM1", "#phi_{p} Distribution for E_{#gamma} 610 #pm 10MeV (Cos#theta_{CM}1-0.8)", 10, -180, 180);
  Phip630CM1 = new GH1("Phip_630MeVCM1", "#phi_{p} Distribution for E_{#gamma} 630 #pm 10MeV (Cos#theta_{CM}1-0.8)", 10, -180, 180);

  // #phi_{p} dists across EGamma bins
  Phip410CM2 = new GH1("Phip_410MeVCM2", "#phi_{p} Distribution for E_{#gamma} 410 #pm 10MeV (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
  Phip430CM2 = new GH1("Phip_430MeVCM2", "#phi_{p} Distribution for E_{#gamma} 430 #pm 10MeV (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
  Phip450CM2 = new GH1("Phip_450MeVCM2", "#phi_{p} Distribution for E_{#gamma} 450 #pm 10MeV (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
  Phip470CM2 = new GH1("Phip_470MeVCM2", "#phi_{p} Distribution for E_{#gamma} 470 #pm 10MeV (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
  Phip490CM2 = new GH1("Phip_490MeVCM2", "#phi_{p} Distribution for E_{#gamma} 490 #pm 10MeV (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
  Phip510CM2 = new GH1("Phip_510MeVCM2", "#phi_{p} Distribution for E_{#gamma} 510 #pm 10MeV (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
  Phip530CM2 = new GH1("Phip_530MeVCM2", "#phi_{p} Distribution for E_{#gamma} 530 #pm 10MeV (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
  Phip550CM2 = new GH1("Phip_550MeVCM2", "#phi_{p} Distribution for E_{#gamma} 550 #pm 10MeV (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
  Phip570CM2 = new GH1("Phip_570MeVCM2", "#phi_{p} Distribution for E_{#gamma} 570 #pm 10MeV (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
  Phip590CM2 = new GH1("Phip_590MeVCM2", "#phi_{p} Distribution for E_{#gamma} 590 #pm 10MeV (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
  Phip610CM2 = new GH1("Phip_610MeVCM2", "#phi_{p} Distribution for E_{#gamma} 610 #pm 10MeV (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
  Phip630CM2 = new GH1("Phip_630MeVCM2", "#phi_{p} Distribution for E_{#gamma} 630 #pm 10MeV (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);

  // #phi_{p} dists across EGamma bins
  Phip410CM3 = new GH1("Phip_410MeVCM3", "#phi_{p} Distribution for E_{#gamma} 410 #pm 10MeV (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
  Phip430CM3 = new GH1("Phip_430MeVCM3", "#phi_{p} Distribution for E_{#gamma} 430 #pm 10MeV (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
  Phip450CM3 = new GH1("Phip_450MeVCM3", "#phi_{p} Distribution for E_{#gamma} 450 #pm 10MeV (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
  Phip470CM3 = new GH1("Phip_470MeVCM3", "#phi_{p} Distribution for E_{#gamma} 470 #pm 10MeV (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
  Phip490CM3 = new GH1("Phip_490MeVCM3", "#phi_{p} Distribution for E_{#gamma} 490 #pm 10MeV (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
  Phip510CM3 = new GH1("Phip_510MeVCM3", "#phi_{p} Distribution for E_{#gamma} 510 #pm 10MeV (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
  Phip530CM3 = new GH1("Phip_530MeVCM3", "#phi_{p} Distribution for E_{#gamma} 530 #pm 10MeV (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
  Phip550CM3 = new GH1("Phip_550MeVCM3", "#phi_{p} Distribution for E_{#gamma} 550 #pm 10MeV (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
  Phip570CM3 = new GH1("Phip_570MeVCM3", "#phi_{p} Distribution for E_{#gamma} 570 #pm 10MeV (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
  Phip590CM3 = new GH1("Phip_590MeVCM3", "#phi_{p} Distribution for E_{#gamma} 590 #pm 10MeV (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
  Phip610CM3 = new GH1("Phip_610MeVCM3", "#phi_{p} Distribution for E_{#gamma} 610 #pm 10MeV (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
  Phip630CM3 = new GH1("Phip_630MeVCM3", "#phi_{p} Distribution for E_{#gamma} 630 #pm 10MeV (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);

  // #phi_{p} dists across EGamma bins
  Phip410CM4 = new GH1("Phip_410MeVCM4", "#phi_{p} Distribution for E_{#gamma} 410 #pm 10MeV (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
  Phip430CM4 = new GH1("Phip_430MeVCM4", "#phi_{p} Distribution for E_{#gamma} 430 #pm 10MeV (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
  Phip450CM4 = new GH1("Phip_450MeVCM4", "#phi_{p} Distribution for E_{#gamma} 450 #pm 10MeV (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
  Phip470CM4 = new GH1("Phip_470MeVCM4", "#phi_{p} Distribution for E_{#gamma} 470 #pm 10MeV (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
  Phip490CM4 = new GH1("Phip_490MeVCM4", "#phi_{p} Distribution for E_{#gamma} 490 #pm 10MeV (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
  Phip510CM4 = new GH1("Phip_510MeVCM4", "#phi_{p} Distribution for E_{#gamma} 510 #pm 10MeV (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
  Phip530CM4 = new GH1("Phip_530MeVCM4", "#phi_{p} Distribution for E_{#gamma} 530 #pm 10MeV (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
  Phip550CM4 = new GH1("Phip_550MeVCM4", "#phi_{p} Distribution for E_{#gamma} 550 #pm 10MeV (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
  Phip570CM4 = new GH1("Phip_570MeVCM4", "#phi_{p} Distribution for E_{#gamma} 570 #pm 10MeV (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
  Phip590CM4 = new GH1("Phip_590MeVCM4", "#phi_{p} Distribution for E_{#gamma} 590 #pm 10MeV (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
  Phip610CM4 = new GH1("Phip_610MeVCM4", "#phi_{p} Distribution for E_{#gamma} 610 #pm 10MeV (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
  Phip630CM4 = new GH1("Phip_630MeVCM4", "#phi_{p} Distribution for E_{#gamma} 630 #pm 10MeV (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);

  // #phi_{p} dists across EGamma bins
  Phip410CM5 = new GH1("Phip_410MeVCM5", "#phi_{p} Distribution for E_{#gamma} 410 #pm 10MeV (Cos#theta_{CM}0.2-0)", 10, -180, 180);
  Phip430CM5 = new GH1("Phip_430MeVCM5", "#phi_{p} Distribution for E_{#gamma} 430 #pm 10MeV (Cos#theta_{CM}0.2-0)", 10, -180, 180);
  Phip450CM5 = new GH1("Phip_450MeVCM5", "#phi_{p} Distribution for E_{#gamma} 450 #pm 10MeV (Cos#theta_{CM}0.2-0)", 10, -180, 180);
  Phip470CM5 = new GH1("Phip_470MeVCM5", "#phi_{p} Distribution for E_{#gamma} 470 #pm 10MeV (Cos#theta_{CM}0.2-0)", 10, -180, 180);
  Phip490CM5 = new GH1("Phip_490MeVCM5", "#phi_{p} Distribution for E_{#gamma} 490 #pm 10MeV (Cos#theta_{CM}0.2-0)", 10, -180, 180);
  Phip510CM5 = new GH1("Phip_510MeVCM5", "#phi_{p} Distribution for E_{#gamma} 510 #pm 10MeV (Cos#theta_{CM}0.2-0)", 10, -180, 180);
  Phip530CM5 = new GH1("Phip_530MeVCM5", "#phi_{p} Distribution for E_{#gamma} 530 #pm 10MeV (Cos#theta_{CM}0.2-0)", 10, -180, 180);
  Phip550CM5 = new GH1("Phip_550MeVCM5", "#phi_{p} Distribution for E_{#gamma} 550 #pm 10MeV (Cos#theta_{CM}0.2-0)", 10, -180, 180);
  Phip570CM5 = new GH1("Phip_570MeVCM5", "#phi_{p} Distribution for E_{#gamma} 570 #pm 10MeV (Cos#theta_{CM}0.2-0)", 10, -180, 180);
  Phip590CM5 = new GH1("Phip_590MeVCM5", "#phi_{p} Distribution for E_{#gamma} 590 #pm 10MeV (Cos#theta_{CM}0.2-0)", 10, -180, 180);
  Phip610CM5 = new GH1("Phip_610MeVCM5", "#phi_{p} Distribution for E_{#gamma} 610 #pm 10MeV (Cos#theta_{CM}0.2-0)", 10, -180, 180);
  Phip630CM5 = new GH1("Phip_630MeVCM5", "#phi_{p} Distribution for E_{#gamma} 630 #pm 10MeV (Cos#theta_{CM}0.2-0)", 10, -180, 180);

  // #phi_{p} dists across EGamma bins
  Phip410CM6 = new GH1("Phip_410MeVCM6", "#phi_{p} Distribution for E_{#gamma} 410 #pm 10MeV (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
  Phip430CM6 = new GH1("Phip_430MeVCM6", "#phi_{p} Distribution for E_{#gamma} 430 #pm 10MeV (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
  Phip450CM6 = new GH1("Phip_450MeVCM6", "#phi_{p} Distribution for E_{#gamma} 450 #pm 10MeV (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
  Phip470CM6 = new GH1("Phip_470MeVCM6", "#phi_{p} Distribution for E_{#gamma} 470 #pm 10MeV (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
  Phip490CM6 = new GH1("Phip_490MeVCM6", "#phi_{p} Distribution for E_{#gamma} 490 #pm 10MeV (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
  Phip510CM6 = new GH1("Phip_510MeVCM6", "#phi_{p} Distribution for E_{#gamma} 510 #pm 10MeV (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
  Phip530CM6 = new GH1("Phip_530MeVCM6", "#phi_{p} Distribution for E_{#gamma} 530 #pm 10MeV (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
  Phip550CM6 = new GH1("Phip_550MeVCM6", "#phi_{p} Distribution for E_{#gamma} 550 #pm 10MeV (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
  Phip570CM6 = new GH1("Phip_570MeVCM6", "#phi_{p} Distribution for E_{#gamma} 570 #pm 10MeV (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
  Phip590CM6 = new GH1("Phip_590MeVCM6", "#phi_{p} Distribution for E_{#gamma} 590 #pm 10MeV (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
  Phip610CM6 = new GH1("Phip_610MeVCM6", "#phi_{p} Distribution for E_{#gamma} 610 #pm 10MeV (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
  Phip630CM6 = new GH1("Phip_630MeVCM6", "#phi_{p} Distribution for E_{#gamma} 630 #pm 10MeV (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);

  // #phi_{p} dists across EGamma bins
  Phip410CM7 = new GH1("Phip_410MeVCM7", "#phi_{p} Distribution for E_{#gamma} 410 #pm 10MeV (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
  Phip430CM7 = new GH1("Phip_430MeVCM7", "#phi_{p} Distribution for E_{#gamma} 430 #pm 10MeV (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
  Phip450CM7 = new GH1("Phip_450MeVCM7", "#phi_{p} Distribution for E_{#gamma} 450 #pm 10MeV (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
  Phip470CM7 = new GH1("Phip_470MeVCM7", "#phi_{p} Distribution for E_{#gamma} 470 #pm 10MeV (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
  Phip490CM7 = new GH1("Phip_490MeVCM7", "#phi_{p} Distribution for E_{#gamma} 490 #pm 10MeV (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
  Phip510CM7 = new GH1("Phip_510MeVCM7", "#phi_{p} Distribution for E_{#gamma} 510 #pm 10MeV (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
  Phip530CM7 = new GH1("Phip_530MeVCM7", "#phi_{p} Distribution for E_{#gamma} 530 #pm 10MeV (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
  Phip550CM7 = new GH1("Phip_550MeVCM7", "#phi_{p} Distribution for E_{#gamma} 550 #pm 10MeV (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
  Phip570CM7 = new GH1("Phip_570MeVCM7", "#phi_{p} Distribution for E_{#gamma} 570 #pm 10MeV (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
  Phip590CM7 = new GH1("Phip_590MeVCM7", "#phi_{p} Distribution for E_{#gamma} 590 #pm 10MeV (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
  Phip610CM7 = new GH1("Phip_610MeVCM7", "#phi_{p} Distribution for E_{#gamma} 610 #pm 10MeV (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
  Phip630CM7 = new GH1("Phip_630MeVCM7", "#phi_{p} Distribution for E_{#gamma} 630 #pm 10MeV (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);

  // #phi_{p} dists across EGamma bins
  Phip410CM8 = new GH1("Phip_410MeVCM8", "#phi_{p} Distribution for E_{#gamma} 410 #pm 10MeV (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
  Phip430CM8 = new GH1("Phip_430MeVCM8", "#phi_{p} Distribution for E_{#gamma} 430 #pm 10MeV (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
  Phip450CM8 = new GH1("Phip_450MeVCM8", "#phi_{p} Distribution for E_{#gamma} 450 #pm 10MeV (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
  Phip470CM8 = new GH1("Phip_470MeVCM8", "#phi_{p} Distribution for E_{#gamma} 470 #pm 10MeV (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
  Phip490CM8 = new GH1("Phip_490MeVCM8", "#phi_{p} Distribution for E_{#gamma} 490 #pm 10MeV (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
  Phip510CM8 = new GH1("Phip_510MeVCM8", "#phi_{p} Distribution for E_{#gamma} 510 #pm 10MeV (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
  Phip530CM8 = new GH1("Phip_530MeVCM8", "#phi_{p} Distribution for E_{#gamma} 530 #pm 10MeV (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
  Phip550CM8 = new GH1("Phip_550MeVCM8", "#phi_{p} Distribution for E_{#gamma} 550 #pm 10MeV (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
  Phip570CM8 = new GH1("Phip_570MeVCM8", "#phi_{p} Distribution for E_{#gamma} 570 #pm 10MeV (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
  Phip590CM8 = new GH1("Phip_590MeVCM8", "#phi_{p} Distribution for E_{#gamma} 590 #pm 10MeV (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
  Phip610CM8 = new GH1("Phip_610MeVCM8", "#phi_{p} Distribution for E_{#gamma} 610 #pm 10MeV (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
  Phip630CM8 = new GH1("Phip_630MeVCM8", "#phi_{p} Distribution for E_{#gamma} 630 #pm 10MeV (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);

  // #phi_{p} dists across EGamma bins
  Phip410CM9 = new GH1("Phip_410MeVCM9", "#phi_{p} Distribution for E_{#gamma} 410 #pm 10MeV (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
  Phip430CM9 = new GH1("Phip_430MeVCM9", "#phi_{p} Distribution for E_{#gamma} 430 #pm 10MeV (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
  Phip450CM9 = new GH1("Phip_450MeVCM9", "#phi_{p} Distribution for E_{#gamma} 450 #pm 10MeV (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
  Phip470CM9 = new GH1("Phip_470MeVCM9", "#phi_{p} Distribution for E_{#gamma} 470 #pm 10MeV (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
  Phip490CM9 = new GH1("Phip_490MeVCM9", "#phi_{p} Distribution for E_{#gamma} 490 #pm 10MeV (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
  Phip510CM9 = new GH1("Phip_510MeVCM9", "#phi_{p} Distribution for E_{#gamma} 510 #pm 10MeV (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
  Phip530CM9 = new GH1("Phip_530MeVCM9", "#phi_{p} Distribution for E_{#gamma} 530 #pm 10MeV (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
  Phip550CM9 = new GH1("Phip_550MeVCM9", "#phi_{p} Distribution for E_{#gamma} 550 #pm 10MeV (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
  Phip570CM9 = new GH1("Phip_570MeVCM9", "#phi_{p} Distribution for E_{#gamma} 570 #pm 10MeV (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
  Phip590CM9 = new GH1("Phip_590MeVCM9", "#phi_{p} Distribution for E_{#gamma} 590 #pm 10MeV (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
  Phip610CM9 = new GH1("Phip_610MeVCM9", "#phi_{p} Distribution for E_{#gamma} 610 #pm 10MeV (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
  Phip630CM9 = new GH1("Phip_630MeVCM9", "#phi_{p} Distribution for E_{#gamma} 630 #pm 10MeV (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);

  // #phi_{p} dists across EGamma bins
  Phip410CM10 = new GH1("Phip_410MeVCM10", "#phi_{p} Distribution for E_{#gamma} 410 #pm 10MeV (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
  Phip430CM10 = new GH1("Phip_430MeVCM10", "#phi_{p} Distribution for E_{#gamma} 430 #pm 10MeV (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
  Phip450CM10 = new GH1("Phip_450MeVCM10", "#phi_{p} Distribution for E_{#gamma} 450 #pm 10MeV (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
  Phip470CM10 = new GH1("Phip_470MeVCM10", "#phi_{p} Distribution for E_{#gamma} 470 #pm 10MeV (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
  Phip490CM10 = new GH1("Phip_490MeVCM10", "#phi_{p} Distribution for E_{#gamma} 490 #pm 10MeV (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
  Phip510CM10 = new GH1("Phip_510MeVCM10", "#phi_{p} Distribution for E_{#gamma} 510 #pm 10MeV (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
  Phip530CM10 = new GH1("Phip_530MeVCM10", "#phi_{p} Distribution for E_{#gamma} 530 #pm 10MeV (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
  Phip550CM10 = new GH1("Phip_550MeVCM10", "#phi_{p} Distribution for E_{#gamma} 550 #pm 10MeV (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
  Phip570CM10 = new GH1("Phip_570MeVCM10", "#phi_{p} Distribution for E_{#gamma} 570 #pm 10MeV (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
  Phip590CM10 = new GH1("Phip_590MeVCM10", "#phi_{p} Distribution for E_{#gamma} 590 #pm 10MeV (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
  Phip610CM10 = new GH1("Phip_610MeVCM10", "#phi_{p} Distribution for E_{#gamma} 610 #pm 10MeV (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
  Phip630CM10 = new GH1("Phip_630MeVCM10", "#phi_{p} Distribution for E_{#gamma} 630 #pm 10MeV (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);

  ThetanDist = new GH1 ("ThetanDist", "#theta_{n} Distribution", 200, 0, 180);
  ThetanRecDist = new GH1 ("ThetanRecDist", "Reconstructed #theta_{n} Distribution", 200, 0, 180);
  ThetanDiffDist = new GH1 ("ThetanDiffDist", "Difference Between #theta_{n} and  #theta_{nRec}", 200, 0, 180);
  ThetanDiffZp = new GH2 ("ThetanDiffZp", "Diff(#theta_{n} - #theta_{nRec}) as a Fn of Z_{p}", 200, 0, 180, 200, -100, 100);
  ThetaRecPiDiff = new GH1 ("ThetaRecPiDiff", "Difference between #theta_{#pi Rec} and #theta_{n}", 200, 0, 180);
  ThetanThetaRecPi = new GH2 ("ThetanThetaRecPi", "#theta_{n} vs #theta_{#pi Rec}", 100, 0, 180, 100, 0, 180);
  ThetanThetaRecPiDiff = new GH2 ("ThetanThetaRecPiDiff", "#theta_{n} vs (#theta_{#pi Rec} - #theta_{n})", 100, 0, 180, 100, 0, 180);
  ThetaRecPDiff = new GH1 ("ThetaRecPDiff", "Difference between #theta_{pRec} and #theta_{n}", 200, 0, 180);
  ThetanThetaRecP = new GH2 ("ThetanThetaRecP", "#theta_{n} vs #theta_{pRec}", 100, 0, 180, 100, 0, 180);
  ThetanThetaRecPDiff = new GH2 ("ThetanThetaRecPDiff", "#theta_{n} vs (#theta_{pRec}  - #theta_{n})", 100, 0, 180, 100, 0, 180);

  E_dE = new GH2 ("E_dE", "EdE Plot With E Loss Adjustment", 100, 0, 500, 100, 0, 5);
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

    ThetanDist->Fill(Thetan, TaggerTime);
    ThetanRecDist->Fill(ThetanRec, TaggerTime);
    ThetanDiffDist->Fill(abs(Thetan-ThetanRec), TaggerTime);
    ThetanDiffZp->Fill(abs(Thetan-ThetanRec), Zp, TaggerTime);
    PhinDiffWCZRec->Fill(WCZnRec, PhinDiff, TaggerTime);
    ThetaRecPiDiff->Fill(ThetaPiRecDiff, TaggerTime);
    ThetanThetaRecPi->Fill(Thetan, ThetaPiRec, TaggerTime);
    ThetanThetaRecPiDiff->Fill(Thetan, ThetaPiRecDiff, TaggerTime);
    ThetaRecPDiff->Fill(ThetapRecDiff, TaggerTime);
    ThetanThetaRecP->Fill(Thetan, ThetapRec, TaggerTime);
    ThetanThetaRecPDiff->Fill(Thetan, ThetapRecDiff, TaggerTime);

    ThetaSc -> Fill(ScattTheta, TaggerTime);
    PhiSc -> Fill(ScattPhi, TaggerTime);
    ThetaScPhiSc->Fill(ScattTheta, ScattPhi, TaggerTime);
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
        OAngle200400->Fill(OpeningAngle, TaggerTime);
    }

    else if(300 < EGamma && EGamma < 400){
        MMp300400->Fill(MMpEpCorr, TaggerTime);
        OAngle200400->Fill(OpeningAngle, TaggerTime);
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

        if(1 > CosThetapCM && CosThetapCM > 0.8 ){
            Phip410CM1->Fill(WCPhip, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.6){
            Phip410CM2->Fill(WCPhip, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.4){
            Phip410CM3->Fill(WCPhip, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.2){
            Phip410CM4->Fill(WCPhip, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0){
            Phip410CM5->Fill(WCPhip, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.2){
            Phip410CM6->Fill(WCPhip, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.4){
            Phip410CM7->Fill(WCPhip, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.6){
            Phip410CM8->Fill(WCPhip, TaggerTime);
        }

        else if(-0.6> CosThetapCM && CosThetapCM > -0.8){
            Phip410CM9->Fill(WCPhip, TaggerTime);
        }

        else if(-0.8> CosThetapCM && CosThetapCM > -1){
            Phip410CM10->Fill(WCPhip, TaggerTime);
        }
    }

    else if ( 420 < EGamma && EGamma < 440) {
        PhiSc430->Fill(ScattPhi, TaggerTime);

        if(1 > CosThetapCM && CosThetapCM > 0.8 ){
            Phip430CM1->Fill(WCPhip, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.6){
            Phip430CM2->Fill(WCPhip, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.4){
            Phip430CM3->Fill(WCPhip, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.2){
            Phip430CM4->Fill(WCPhip, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0){
            Phip430CM5->Fill(WCPhip, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.2){
            Phip430CM6->Fill(WCPhip, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.4){
            Phip430CM7->Fill(WCPhip, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.6){
            Phip430CM8->Fill(WCPhip, TaggerTime);
        }

        else if(-0.6> CosThetapCM && CosThetapCM > -0.8){
            Phip430CM9->Fill(WCPhip, TaggerTime);
        }

        else if(-0.8> CosThetapCM && CosThetapCM > -1){
            Phip430CM10->Fill(WCPhip, TaggerTime);
        }
    }

    else if ( 440 < EGamma && EGamma < 460) {
        PhiSc450->Fill(ScattPhi, TaggerTime);

        if(1 > CosThetapCM && CosThetapCM > 0.8 ){
            Phip450CM1->Fill(WCPhip, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.6){
            Phip450CM2->Fill(WCPhip, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.4){
            Phip450CM3->Fill(WCPhip, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.2){
            Phip450CM4->Fill(WCPhip, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0){
            Phip450CM5->Fill(WCPhip, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.2){
            Phip450CM6->Fill(WCPhip, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.4){
            Phip450CM7->Fill(WCPhip, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.6){
            Phip450CM8->Fill(WCPhip, TaggerTime);
        }

        else if(-0.6> CosThetapCM && CosThetapCM > -0.8){
            Phip450CM9->Fill(WCPhip, TaggerTime);
        }

        else if(-0.8> CosThetapCM && CosThetapCM > -1){
            Phip450CM10->Fill(WCPhip, TaggerTime);
        }
    }

    else if ( 460 < EGamma && EGamma < 480) {
        PhiSc470->Fill(ScattPhi, TaggerTime);

        if(1 > CosThetapCM && CosThetapCM > 0.8 ){
            Phip470CM1->Fill(WCPhip, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.6){
            Phip470CM2->Fill(WCPhip, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.4){
            Phip470CM3->Fill(WCPhip, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.2){
            Phip470CM4->Fill(WCPhip, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0){
            Phip470CM5->Fill(WCPhip, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.2){
            Phip470CM6->Fill(WCPhip, TaggerTime);
        }
        else if(-0.2 > CosThetapCM && CosThetapCM > -0.4){
            Phip470CM7->Fill(WCPhip, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.6){
            Phip470CM8->Fill(WCPhip, TaggerTime);
        }

        else if(-0.6> CosThetapCM && CosThetapCM > -0.8){
            Phip470CM9->Fill(WCPhip, TaggerTime);
        }

        else if(-0.8> CosThetapCM && CosThetapCM > -1){
            Phip470CM10->Fill(WCPhip, TaggerTime);
        }
    }

    else if ( 480 < EGamma && EGamma < 500) {
        PhiSc490->Fill(ScattPhi, TaggerTime);

        if(1 > CosThetapCM && CosThetapCM > 0.8 ){
            Phip490CM1->Fill(WCPhip, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.6){
            Phip490CM2->Fill(WCPhip, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.4){
            Phip490CM3->Fill(WCPhip, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.2){
            Phip490CM4->Fill(WCPhip, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0){
            Phip490CM5->Fill(WCPhip, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.2){
            Phip490CM6->Fill(WCPhip, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.4){
            Phip490CM7->Fill(WCPhip, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.6){
            Phip490CM8->Fill(WCPhip, TaggerTime);
        }

        else if(-0.6> CosThetapCM && CosThetapCM > -0.8){
            Phip490CM9->Fill(WCPhip, TaggerTime);
        }

        else if(-0.8> CosThetapCM && CosThetapCM > -1){
            Phip490CM10->Fill(WCPhip, TaggerTime);
        }
    }

    else if ( 500 < EGamma && EGamma < 520) {
        PhiSc510->Fill(ScattPhi, TaggerTime);

        if(1 > CosThetapCM && CosThetapCM > 0.8 ){
            Phip510CM1->Fill(WCPhip, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.6){
            Phip510CM2->Fill(WCPhip, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.4){
            Phip510CM3->Fill(WCPhip, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.2){
            Phip510CM4->Fill(WCPhip, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0){
            Phip510CM5->Fill(WCPhip, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.2){
            Phip510CM6->Fill(WCPhip, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.4){
            Phip510CM7->Fill(WCPhip, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.6){
            Phip510CM8->Fill(WCPhip, TaggerTime);
        }

        else if(-0.6> CosThetapCM && CosThetapCM > -0.8){
            Phip510CM9->Fill(WCPhip, TaggerTime);
        }

        else if(-0.8> CosThetapCM && CosThetapCM > -1){
            Phip510CM10->Fill(WCPhip, TaggerTime);
        }
    }

    else if ( 520 < EGamma && EGamma < 540) {
        PhiSc530->Fill(ScattPhi, TaggerTime);

        if(1 > CosThetapCM && CosThetapCM > 0.8 ){
            Phip530CM1->Fill(WCPhip, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.6){
            Phip530CM2->Fill(WCPhip, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.4){
            Phip530CM3->Fill(WCPhip, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.2){
            Phip530CM4->Fill(WCPhip, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0){
            Phip530CM5->Fill(WCPhip, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.2){
            Phip530CM6->Fill(WCPhip, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.4){
            Phip530CM7->Fill(WCPhip, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.6){
            Phip530CM8->Fill(WCPhip, TaggerTime);
        }

        else if(-0.6> CosThetapCM && CosThetapCM > -0.8){
            Phip530CM9->Fill(WCPhip, TaggerTime);
        }

        else if(-0.8> CosThetapCM && CosThetapCM > -1){
            Phip530CM10->Fill(WCPhip, TaggerTime);
        }
    }

    else if ( 540 < EGamma && EGamma < 560) {
        PhiSc550->Fill(ScattPhi, TaggerTime);

        if(1 > CosThetapCM && CosThetapCM > 0.8 ){
            Phip550CM1->Fill(WCPhip, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.6){
            Phip550CM2->Fill(WCPhip, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.4){
            Phip550CM3->Fill(WCPhip, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.2){
            Phip550CM4->Fill(WCPhip, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0){
            Phip550CM5->Fill(WCPhip, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.2){
            Phip550CM6->Fill(WCPhip, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.4){
            Phip550CM7->Fill(WCPhip, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.6){
            Phip550CM8->Fill(WCPhip, TaggerTime);
        }

        else if(-0.6> CosThetapCM && CosThetapCM > -0.8){
            Phip550CM9->Fill(WCPhip, TaggerTime);
        }

        else if(-0.8> CosThetapCM && CosThetapCM > -1){
            Phip550CM10->Fill(WCPhip, TaggerTime);
        }
    }

    else if ( 560 < EGamma && EGamma < 580) {
        PhiSc570->Fill(ScattPhi, TaggerTime);

        if(1 > CosThetapCM && CosThetapCM > 0.8 ){
            Phip570CM1->Fill(WCPhip, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.6){
            Phip570CM2->Fill(WCPhip, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.4){
            Phip570CM3->Fill(WCPhip, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.2){
            Phip570CM4->Fill(WCPhip, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0){
            Phip570CM5->Fill(WCPhip, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.2){
            Phip570CM6->Fill(WCPhip, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.4){
            Phip570CM7->Fill(WCPhip, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.6){
            Phip570CM8->Fill(WCPhip, TaggerTime);
        }

        else if(-0.6> CosThetapCM && CosThetapCM > -0.8){
            Phip570CM9->Fill(WCPhip, TaggerTime);
        }

        else if(-0.8> CosThetapCM && CosThetapCM > -1){
            Phip570CM10->Fill(WCPhip, TaggerTime);
        }
    }

    else if ( 580 < EGamma && EGamma < 600) {
        PhiSc590->Fill(ScattPhi, TaggerTime);

        if(1 > CosThetapCM && CosThetapCM > 0.8 ){
            Phip590CM1->Fill(WCPhip, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.6){
            Phip590CM2->Fill(WCPhip, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.4){
            Phip590CM3->Fill(WCPhip, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.2){
            Phip590CM4->Fill(WCPhip, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0){
            Phip590CM5->Fill(WCPhip, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.2){
            Phip590CM6->Fill(WCPhip, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.4){
            Phip590CM7->Fill(WCPhip, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.6){
            Phip590CM8->Fill(WCPhip, TaggerTime);
        }

        else if(-0.6> CosThetapCM && CosThetapCM > -0.8){
            Phip590CM9->Fill(WCPhip, TaggerTime);
        }

        else if(-0.8> CosThetapCM && CosThetapCM > -1){
            Phip590CM10->Fill(WCPhip, TaggerTime);
        }
    }

    else if ( 600 < EGamma && EGamma < 620) {
        PhiSc610->Fill(ScattPhi, TaggerTime);

        if(1 > CosThetapCM && CosThetapCM > 0.8 ){
            Phip610CM1->Fill(WCPhip, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.6){
            Phip610CM2->Fill(WCPhip, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.4){
            Phip610CM3->Fill(WCPhip, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.2){
            Phip610CM4->Fill(WCPhip, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0){
            Phip610CM5->Fill(WCPhip, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.2){
            Phip610CM6->Fill(WCPhip, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.4){
            Phip610CM7->Fill(WCPhip, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.6){
            Phip610CM8->Fill(WCPhip, TaggerTime);
        }

        else if(-0.6> CosThetapCM && CosThetapCM > -0.8){
            Phip610CM9->Fill(WCPhip, TaggerTime);
        }

        else if(-0.8> CosThetapCM && CosThetapCM > -1){
            Phip610CM10->Fill(WCPhip, TaggerTime);
        }
    }

    else if ( 620 < EGamma && EGamma < 640) {
        PhiSc630->Fill(ScattPhi, TaggerTime);

        if(1 > CosThetapCM && CosThetapCM > 0.8 ){
            Phip630CM1->Fill(WCPhip, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.6){
            Phip630CM2->Fill(WCPhip, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.4){
            Phip630CM3->Fill(WCPhip, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.2){
            Phip630CM4->Fill(WCPhip, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0){
            Phip630CM5->Fill(WCPhip, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.2){
            Phip630CM6->Fill(WCPhip, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.4){
            Phip630CM7->Fill(WCPhip, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.6){
            Phip630CM8->Fill(WCPhip, TaggerTime);
        }

        else if(-0.6> CosThetapCM && CosThetapCM > -0.8){
            Phip630CM9->Fill(WCPhip, TaggerTime);
        }

        else if(-0.8> CosThetapCM && CosThetapCM > -1){
            Phip630CM10->Fill(WCPhip, TaggerTime);
        }
    }
}

Bool_t	PNeutPol_Polarimeter_Lin::Write(){
  // Write all GH1's easily

  GTreeManager::Write();
}
