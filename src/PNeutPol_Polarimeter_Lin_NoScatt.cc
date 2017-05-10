// GoAT Physics analysis to identify neutrons from deuterium photodisintegration
// Various properties of neutrons/protons identified plotted in histograms
// Main aim is to determine spin polarisation of neutrons
// For use on linearly polarised data files

#include "PNeutPol_Polarimeter_Lin_NoScatt.h"

PNeutPol_Polarimeter_Lin_NoScatt::~PNeutPol_Polarimeter_Lin_NoScatt()
{
}

Bool_t	PNeutPol_Polarimeter_Lin_NoScatt::Init()
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

Bool_t	PNeutPol_Polarimeter_Lin_NoScatt::Start()
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

  return kTRUE;
}

void	PNeutPol_Polarimeter_Lin_NoScatt::ProcessEvent()
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

  if((Detectors1 == 7) && (Detectors2 == 1))
  {
    Proton1 = kTRUE;
    Proton2 = kFALSE;
    if (GetTracks()->GetMWPC0Energy(0) == 0) return; // If no hit in first chamber for p drop out
  }

  // If track 2 only gives signals in MWPC and CB it is the neutron
  else if((Detectors1 == 1) && (Detectors2 == 7))
  {
    Proton1 = kFALSE;
    Proton2 = kTRUE;
    if (GetTracks()->GetMWPC0Energy(1) == 0) return; // If no hit in first chamber for p drop out
  }

  // Drop out on ANY other condition (for now)
  else
  {
    return;
  }

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
  }

  else
  {
    return;
  }

  GVp3 = GVp.Vect(); // Generate some 3-vectors from the 4-vectors we have
  GVn3 = GVn.Vect();
  WC3Vectp.SetXYZ(WC1pX, WC1pY, WC1pZ);
  WCThetap = WC3Vectp.Theta()*TMath::RadToDeg(); // Angles from WC hit positons
  WCThetapRad = WC3Vectp.Theta();
  WCPhip = WC3Vectp.Phi()*TMath::RadToDeg();
  WCPhipRad = WC3Vectp.Phi();

  if( Zp > 60 || Zp < -60) return; // Particles selected out from other parts tend to be inside anyway, skip this?
  //if( Zp > 200 || Zp < 150) return; // Select out windows

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

    if(Cut_protonKinGood -> IsInside(KinEp, dEp) == kFALSE) continue; // If E loss correct proton is NOT inside p banana drop out
    if(ThetaPiRec > 20) continue;
    //if ( 850 > MMpEpCorr || 1050 < MMpEpCorr) continue;

    //if (abs(KinEDiff) > 100) continue; // If difference between CB energy and calculated Energy for proton > 100MeV continue

    //k++;
    FillHists(); // Fill histograms with data generated
  }
}

void	PNeutPol_Polarimeter_Lin_NoScatt::ProcessScalerRead()
{
	// Fill Tagger Scalers // Currently this seems to fill the file with loads of "TaggerAccScal" histograms
	//FillScalers(GetTC_scaler_min(),GetTC_scaler_max(),TaggerAccScal); // Don't know if these are needed so cut out for now
}

TCutG*	PNeutPol_Polarimeter_Lin_NoScatt::OpenCutFile(Char_t* filename, Char_t* cutname)
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

TLorentzVector PNeutPol_Polarimeter_Lin_NoScatt::LProton4VectorKin(Double_t KinE, Double_t Theta, Double_t Phi)
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

TLorentzVector PNeutPol_Polarimeter_Lin_NoScatt::LNeutron4VectorKin(TLorentzVector ProtonKinVector)
{

    N4Vect = (Gamma + Deut) - ProtonKinVector;

    return N4Vect;
}

TLorentzVector PNeutPol_Polarimeter_Lin_NoScatt::LPion4VectorKin(TLorentzVector ProtonKinVector)
{

    Pi4Vect = (Gamma + Neut) - ProtonKinVector;

    return Pi4Vect;
}

PNeutPol_Polarimeter_Lin_NoScatt::PNeutPol_Polarimeter_Lin_NoScatt() // Define a load of histograms to fill
{
  time = new TH1D("time", 	"time", 	1400, -700, 700);
  time_cut = new TH1D("time_cut", 	"time_cut", 	1400, -700, 700);

  Eg = new GH1( "Eg", "Photon Energy Distribution", 200, 100, 1600);
  EpKin = new GH1 ("EpKin", "Ep Calculated from Ep/Thetap", 100, 0, 500);
  EpCorrected = new GH1 ("EpCorrected", "Ep Corrected for Energy Loss in Polarimeter ", 100, 0, 500);
  OAngle = new GH1 ("OAngle", "Opening Angle between P and N Vectors", 180, 0, 180);
  WCZnRecon = new GH1 ("WCZnRecon", "WCZ Hit Position from Reconstructed n Vector", 200, 0, 400);

  EpKinEpCorrDiff = new GH1("EpKinEpCorrDiff", "Difference Between EpKin and EpCorr", 300, -300, 300);
  EpEpCorrDiff = new GH1("EpEpCorrDiff", "Difference Between Ep and EpCorr", 200, 0, 200);
  OAngle200400 = new GH1 ("OAngle200400", "Opening Angle between P and N Vectors (P Banana Cut, 200-400MeV Gamma)", 180, 0, 180);
  MMpEpCorrected = new GH1 ("MMpEpCorrected", "Missing mass seen by Proton (E Loss Corrected)", 400, 0, 2000);

  ZpDist = new GH1 ("ZpDist", "Proton Pseudo Z Vertex Distribution", 200, -400, 400);

  // MMp across photon E bins
  MMp200300 = new GH1("MMp200300", "Missing mass as seen by Proton (200-300MeV Photon Energy)", 400, 0, 2000);
  MMp300400 = new GH1("MMp300400", "Missing mass as seen by Proton (300-400MeV Photon Energy)", 400, 0, 2000);
  MMp400500 = new GH1("MMp400500", "Missing mass as seen by Proton (400-500MeV Photon Energy)", 400, 0, 2000);
  MMp500600 = new GH1("MMp500600", "Missing mass as seen by Proton (500-600MeV Photon Energy)", 400, 0, 2000);
  MMp600700 = new GH1("MMp600700", "Missing mass as seen by Proton (600-700MeV Photon Energy)", 400, 0, 2000);
  MMp700800 = new GH1("MMp700800", "Missing mass as seen by Proton (700-800MeV Photon Energy)", 400, 0, 2000);
  MMp800900 = new GH1("MMp800900", "Missing mass as seen by Proton (800-900MeV Photon Energy)", 400, 0, 2000);

  // Proton Phi dists across EGamma bins
  Phip425CM1 = new GH1("Phip_425MeVCM1", "Proton Phi Distribution for Photon Energies of 425pm5MeV (CosTheta1-0.8)", 10, -180, 180);
  Phip435CM1 = new GH1("Phip_435MeVCM1", "Proton Phi Distribution for Photon Energies of 435pm5MeV (CosTheta1-0.8)", 10, -180, 180);
  Phip445CM1 = new GH1("Phip_445MeVCM1", "Proton Phi Distribution for Photon Energies of 445pm5MeV (CosTheta1-0.8)", 10, -180, 180);
  Phip455CM1 = new GH1("Phip_455MeVCM1", "Proton Phi Distribution for Photon Energies of 455pm5MeV (CosTheta1-0.8)", 10, -180, 180);
  Phip465CM1 = new GH1("Phip_465MeVCM1", "Proton Phi Distribution for Photon Energies of 465pm5MeV (CosTheta1-0.8)", 10, -180, 180);
  Phip475CM1 = new GH1("Phip_475MeVCM1", "Proton Phi Distribution for Photon Energies of 475pm5MeV (CosTheta1-0.8)", 10, -180, 180);
  Phip485CM1 = new GH1("Phip_485MeVCM1", "Proton Phi Distribution for Photon Energies of 485pm5MeV (CosTheta1-0.8)", 10, -180, 180);
  Phip495CM1 = new GH1("Phip_495MeVCM1", "Proton Phi Distribution for Photon Energies of 495pm5MeV (CosTheta1-0.8)", 10, -180, 180);
  Phip505CM1 = new GH1("Phip_505MeVCM1", "Proton Phi Distribution for Photon Energies of 505pm5MeV (CosTheta1-0.8)", 10, -180, 180);
  Phip515CM1 = new GH1("Phip_515MeVCM1", "Proton Phi Distribution for Photon Energies of 515pm5MeV (CosTheta1-0.8)", 10, -180, 180);
  Phip525CM1 = new GH1("Phip_525MeVCM1", "Proton Phi Distribution for Photon Energies of 525pm5MeV (CosTheta1-0.8)", 10, -180, 180);
  Phip535CM1 = new GH1("Phip_535MeVCM1", "Proton Phi Distribution for Photon Energies of 535pm5MeV (CosTheta1-0.8)", 10, -180, 180);
  Phip545CM1 = new GH1("Phip_545MeVCM1", "Proton Phi Distribution for Photon Energies of 545pm5MeV (CosTheta1-0.8)", 10, -180, 180);
  Phip555CM1 = new GH1("Phip_555MeVCM1", "Proton Phi Distribution for Photon Energies of 555pm5MeV (CosTheta1-0.8)", 10, -180, 180);
  Phip565CM1 = new GH1("Phip_565MeVCM1", "Proton Phi Distribution for Photon Energies of 535pm5MeV (CosTheta1-0.8)", 10, -180, 180);
  Phip575CM1 = new GH1("Phip_575MeVCM1", "Proton Phi Distribution for Photon Energies of 575pm5MeV (CosTheta1-0.8)", 10, -180, 180);
  Phip585CM1 = new GH1("Phip_585MeVCM1", "Proton Phi Distribution for Photon Energies of 585pm5MeV (CosTheta1-0.8)", 10, -180, 180);
  Phip595CM1 = new GH1("Phip_595MeVCM1", "Proton Phi Distribution for Photon Energies of 595pm5MeV (CosTheta1-0.8)", 10, -180, 180);
  Phip605CM1 = new GH1("Phip_605MeVCM1", "Proton Phi Distribution for Photon Energies of 605pm5MeV (CosTheta1-0.8)", 10, -180, 180);
  Phip615CM1 = new GH1("Phip_615MeVCM1", "Proton Phi Distribution for Photon Energies of 615pm5MeV (CosTheta1-0.8)", 10, -180, 180);
  Phip625CM1 = new GH1("Phip_625MeVCM1", "Proton Phi Distribution for Photon Energies of 625pm5MeV (CosTheta1-0.8)", 10, -180, 180);

  // Proton Phi dists across EGamma bins
  Phip425CM2 = new GH1("Phip_425MeVCM2", "Proton Phi Distribution for Photon Energies of 425pm5MeV (CosTheta0.8-0.6)", 10, -180, 180);
  Phip435CM2 = new GH1("Phip_435MeVCM2", "Proton Phi Distribution for Photon Energies of 435pm5MeV (CosTheta0.8-0.6)", 10, -180, 180);
  Phip445CM2 = new GH1("Phip_445MeVCM2", "Proton Phi Distribution for Photon Energies of 445pm5MeV (CosTheta0.8-0.6)", 10, -180, 180);
  Phip455CM2 = new GH1("Phip_455MeVCM2", "Proton Phi Distribution for Photon Energies of 455pm5MeV (CosTheta0.8-0.6)", 10, -180, 180);
  Phip465CM2 = new GH1("Phip_465MeVCM2", "Proton Phi Distribution for Photon Energies of 465pm5MeV (CosTheta0.8-0.6)", 10, -180, 180);
  Phip475CM2 = new GH1("Phip_475MeVCM2", "Proton Phi Distribution for Photon Energies of 475pm5MeV (CosTheta0.8-0.6)", 10, -180, 180);
  Phip485CM2 = new GH1("Phip_485MeVCM2", "Proton Phi Distribution for Photon Energies of 485pm5MeV (CosTheta0.8-0.6)", 10, -180, 180);
  Phip495CM2 = new GH1("Phip_495MeVCM2", "Proton Phi Distribution for Photon Energies of 495pm5MeV (CosTheta0.8-0.6)", 10, -180, 180);
  Phip505CM2 = new GH1("Phip_505MeVCM2", "Proton Phi Distribution for Photon Energies of 505pm5MeV (CosTheta0.8-0.6)", 10, -180, 180);
  Phip515CM2 = new GH1("Phip_515MeVCM2", "Proton Phi Distribution for Photon Energies of 515pm5MeV (CosTheta0.8-0.6)", 10, -180, 180);
  Phip525CM2 = new GH1("Phip_525MeVCM2", "Proton Phi Distribution for Photon Energies of 525pm5MeV (CosTheta0.8-0.6)", 10, -180, 180);
  Phip535CM2 = new GH1("Phip_535MeVCM2", "Proton Phi Distribution for Photon Energies of 535pm5MeV (CosTheta0.8-0.6)", 10, -180, 180);
  Phip545CM2 = new GH1("Phip_545MeVCM2", "Proton Phi Distribution for Photon Energies of 545pm5MeV (CosTheta0.8-0.6)", 10, -180, 180);
  Phip555CM2 = new GH1("Phip_555MeVCM2", "Proton Phi Distribution for Photon Energies of 555pm5MeV (CosTheta0.8-0.6)", 10, -180, 180);
  Phip565CM2 = new GH1("Phip_565MeVCM2", "Proton Phi Distribution for Photon Energies of 535pm5MeV (CosTheta0.8-0.6)", 10, -180, 180);
  Phip575CM2 = new GH1("Phip_575MeVCM2", "Proton Phi Distribution for Photon Energies of 575pm5MeV (CosTheta0.8-0.6)", 10, -180, 180);
  Phip585CM2 = new GH1("Phip_585MeVCM2", "Proton Phi Distribution for Photon Energies of 585pm5MeV (CosTheta0.8-0.6)", 10, -180, 180);
  Phip595CM2 = new GH1("Phip_595MeVCM2", "Proton Phi Distribution for Photon Energies of 595pm5MeV (CosTheta0.8-0.6)", 10, -180, 180);
  Phip605CM2 = new GH1("Phip_605MeVCM2", "Proton Phi Distribution for Photon Energies of 605pm5MeV (CosTheta0.8-0.6)", 10, -180, 180);
  Phip615CM2 = new GH1("Phip_615MeVCM2", "Proton Phi Distribution for Photon Energies of 615pm5MeV (CosTheta0.8-0.6)", 10, -180, 180);
  Phip625CM2 = new GH1("Phip_625MeVCM2", "Proton Phi Distribution for Photon Energies of 625pm5MeV (CosTheta0.8-0.6)", 10, -180, 180);

  // Proton Phi dists across EGamma bins
  Phip425CM3 = new GH1("Phip_425MeVCM3", "Proton Phi Distribution for Photon Energies of 425pm5MeV (CosTheta0.6-0.4)", 10, -180, 180);
  Phip435CM3 = new GH1("Phip_435MeVCM3", "Proton Phi Distribution for Photon Energies of 435pm5MeV (CosTheta0.6-0.4)", 10, -180, 180);
  Phip445CM3 = new GH1("Phip_445MeVCM3", "Proton Phi Distribution for Photon Energies of 445pm5MeV (CosTheta0.6-0.4)", 10, -180, 180);
  Phip455CM3 = new GH1("Phip_455MeVCM3", "Proton Phi Distribution for Photon Energies of 455pm5MeV (CosTheta0.6-0.4)", 10, -180, 180);
  Phip465CM3 = new GH1("Phip_465MeVCM3", "Proton Phi Distribution for Photon Energies of 465pm5MeV (CosTheta0.6-0.4)", 10, -180, 180);
  Phip475CM3 = new GH1("Phip_475MeVCM3", "Proton Phi Distribution for Photon Energies of 475pm5MeV (CosTheta0.6-0.4)", 10, -180, 180);
  Phip485CM3 = new GH1("Phip_485MeVCM3", "Proton Phi Distribution for Photon Energies of 485pm5MeV (CosTheta0.6-0.4)", 10, -180, 180);
  Phip495CM3 = new GH1("Phip_495MeVCM3", "Proton Phi Distribution for Photon Energies of 495pm5MeV (CosTheta0.6-0.4)", 10, -180, 180);
  Phip505CM3 = new GH1("Phip_505MeVCM3", "Proton Phi Distribution for Photon Energies of 505pm5MeV (CosTheta0.6-0.4)", 10, -180, 180);
  Phip515CM3 = new GH1("Phip_515MeVCM3", "Proton Phi Distribution for Photon Energies of 515pm5MeV (CosTheta0.6-0.4)", 10, -180, 180);
  Phip525CM3 = new GH1("Phip_525MeVCM3", "Proton Phi Distribution for Photon Energies of 525pm5MeV (CosTheta0.6-0.4)", 10, -180, 180);
  Phip535CM3 = new GH1("Phip_535MeVCM3", "Proton Phi Distribution for Photon Energies of 535pm5MeV (CosTheta0.6-0.4)", 10, -180, 180);
  Phip545CM3 = new GH1("Phip_545MeVCM3", "Proton Phi Distribution for Photon Energies of 545pm5MeV (CosTheta0.6-0.4)", 10, -180, 180);
  Phip555CM3 = new GH1("Phip_555MeVCM3", "Proton Phi Distribution for Photon Energies of 555pm5MeV (CosTheta0.6-0.4)", 10, -180, 180);
  Phip565CM3 = new GH1("Phip_565MeVCM3", "Proton Phi Distribution for Photon Energies of 535pm5MeV (CosTheta0.6-0.4)", 10, -180, 180);
  Phip575CM3 = new GH1("Phip_575MeVCM3", "Proton Phi Distribution for Photon Energies of 575pm5MeV (CosTheta0.6-0.4)", 10, -180, 180);
  Phip585CM3 = new GH1("Phip_585MeVCM3", "Proton Phi Distribution for Photon Energies of 585pm5MeV (CosTheta0.6-0.4)", 10, -180, 180);
  Phip595CM3 = new GH1("Phip_595MeVCM3", "Proton Phi Distribution for Photon Energies of 595pm5MeV (CosTheta0.6-0.4)", 10, -180, 180);
  Phip605CM3 = new GH1("Phip_605MeVCM3", "Proton Phi Distribution for Photon Energies of 605pm5MeV (CosTheta0.6-0.4)", 10, -180, 180);
  Phip615CM3 = new GH1("Phip_615MeVCM3", "Proton Phi Distribution for Photon Energies of 615pm5MeV (CosTheta0.6-0.4)", 10, -180, 180);
  Phip625CM3 = new GH1("Phip_625MeVCM3", "Proton Phi Distribution for Photon Energies of 625pm5MeV (CosTheta0.6-0.4)", 10, -180, 180);

  // Proton Phi dists across EGamma bins
  Phip425CM4 = new GH1("Phip_425MeVCM4", "Proton Phi Distribution for Photon Energies of 425pm5MeV (CosTheta0.4-0.2)", 10, -180, 180);
  Phip435CM4 = new GH1("Phip_435MeVCM4", "Proton Phi Distribution for Photon Energies of 435pm5MeV (CosTheta0.4-0.2)", 10, -180, 180);
  Phip445CM4 = new GH1("Phip_445MeVCM4", "Proton Phi Distribution for Photon Energies of 445pm5MeV (CosTheta0.4-0.2)", 10, -180, 180);
  Phip455CM4 = new GH1("Phip_455MeVCM4", "Proton Phi Distribution for Photon Energies of 455pm5MeV (CosTheta0.4-0.2)", 10, -180, 180);
  Phip465CM4 = new GH1("Phip_465MeVCM4", "Proton Phi Distribution for Photon Energies of 465pm5MeV (CosTheta0.4-0.2)", 10, -180, 180);
  Phip475CM4 = new GH1("Phip_475MeVCM4", "Proton Phi Distribution for Photon Energies of 475pm5MeV (CosTheta0.4-0.2)", 10, -180, 180);
  Phip485CM4 = new GH1("Phip_485MeVCM4", "Proton Phi Distribution for Photon Energies of 485pm5MeV (CosTheta0.4-0.2)", 10, -180, 180);
  Phip495CM4 = new GH1("Phip_495MeVCM4", "Proton Phi Distribution for Photon Energies of 495pm5MeV (CosTheta0.4-0.2)", 10, -180, 180);
  Phip505CM4 = new GH1("Phip_505MeVCM4", "Proton Phi Distribution for Photon Energies of 505pm5MeV (CosTheta0.4-0.2)", 10, -180, 180);
  Phip515CM4 = new GH1("Phip_515MeVCM4", "Proton Phi Distribution for Photon Energies of 515pm5MeV (CosTheta0.4-0.2)", 10, -180, 180);
  Phip525CM4 = new GH1("Phip_525MeVCM4", "Proton Phi Distribution for Photon Energies of 525pm5MeV (CosTheta0.4-0.2)", 10, -180, 180);
  Phip535CM4 = new GH1("Phip_535MeVCM4", "Proton Phi Distribution for Photon Energies of 535pm5MeV (CosTheta0.4-0.2)", 10, -180, 180);
  Phip545CM4 = new GH1("Phip_545MeVCM4", "Proton Phi Distribution for Photon Energies of 545pm5MeV (CosTheta0.4-0.2)", 10, -180, 180);
  Phip555CM4 = new GH1("Phip_555MeVCM4", "Proton Phi Distribution for Photon Energies of 555pm5MeV (CosTheta0.4-0.2)", 10, -180, 180);
  Phip565CM4 = new GH1("Phip_565MeVCM4", "Proton Phi Distribution for Photon Energies of 535pm5MeV (CosTheta0.4-0.2)", 10, -180, 180);
  Phip575CM4 = new GH1("Phip_575MeVCM4", "Proton Phi Distribution for Photon Energies of 575pm5MeV (CosTheta0.4-0.2)", 10, -180, 180);
  Phip585CM4 = new GH1("Phip_585MeVCM4", "Proton Phi Distribution for Photon Energies of 585pm5MeV (CosTheta0.4-0.2)", 10, -180, 180);
  Phip595CM4 = new GH1("Phip_595MeVCM4", "Proton Phi Distribution for Photon Energies of 595pm5MeV (CosTheta0.4-0.2)", 10, -180, 180);
  Phip605CM4 = new GH1("Phip_605MeVCM4", "Proton Phi Distribution for Photon Energies of 605pm5MeV (CosTheta0.4-0.2)", 10, -180, 180);
  Phip615CM4 = new GH1("Phip_615MeVCM4", "Proton Phi Distribution for Photon Energies of 615pm5MeV (CosTheta0.4-0.2)", 10, -180, 180);
  Phip625CM4 = new GH1("Phip_625MeVCM4", "Proton Phi Distribution for Photon Energies of 625pm5MeV (CosTheta0.4-0.2)", 10, -180, 180);

  // Proton Phi dists across EGamma bins
  Phip425CM5 = new GH1("Phip_425MeVCM5", "Proton Phi Distribution for Photon Energies of 425pm5MeV (CosTheta0.2-0)", 10, -180, 180);
  Phip435CM5 = new GH1("Phip_435MeVCM5", "Proton Phi Distribution for Photon Energies of 435pm5MeV (CosTheta0.2-0)", 10, -180, 180);
  Phip445CM5 = new GH1("Phip_445MeVCM5", "Proton Phi Distribution for Photon Energies of 445pm5MeV (CosTheta0.2-0)", 10, -180, 180);
  Phip455CM5 = new GH1("Phip_455MeVCM5", "Proton Phi Distribution for Photon Energies of 455pm5MeV (CosTheta0.2-0)", 10, -180, 180);
  Phip465CM5 = new GH1("Phip_465MeVCM5", "Proton Phi Distribution for Photon Energies of 465pm5MeV (CosTheta0.2-0)", 10, -180, 180);
  Phip475CM5 = new GH1("Phip_475MeVCM5", "Proton Phi Distribution for Photon Energies of 475pm5MeV (CosTheta0.2-0)", 10, -180, 180);
  Phip485CM5 = new GH1("Phip_485MeVCM5", "Proton Phi Distribution for Photon Energies of 485pm5MeV (CosTheta0.2-0)", 10, -180, 180);
  Phip495CM5 = new GH1("Phip_495MeVCM5", "Proton Phi Distribution for Photon Energies of 495pm5MeV (CosTheta0.2-0)", 10, -180, 180);
  Phip505CM5 = new GH1("Phip_505MeVCM5", "Proton Phi Distribution for Photon Energies of 505pm5MeV (CosTheta0.2-0)", 10, -180, 180);
  Phip515CM5 = new GH1("Phip_515MeVCM5", "Proton Phi Distribution for Photon Energies of 515pm5MeV (CosTheta0.2-0)", 10, -180, 180);
  Phip525CM5 = new GH1("Phip_525MeVCM5", "Proton Phi Distribution for Photon Energies of 525pm5MeV (CosTheta0.2-0)", 10, -180, 180);
  Phip535CM5 = new GH1("Phip_535MeVCM5", "Proton Phi Distribution for Photon Energies of 535pm5MeV (CosTheta0.2-0)", 10, -180, 180);
  Phip545CM5 = new GH1("Phip_545MeVCM5", "Proton Phi Distribution for Photon Energies of 545pm5MeV (CosTheta0.2-0)", 10, -180, 180);
  Phip555CM5 = new GH1("Phip_555MeVCM5", "Proton Phi Distribution for Photon Energies of 555pm5MeV (CosTheta0.2-0)", 10, -180, 180);
  Phip565CM5 = new GH1("Phip_565MeVCM5", "Proton Phi Distribution for Photon Energies of 535pm5MeV (CosTheta0.2-0)", 10, -180, 180);
  Phip575CM5 = new GH1("Phip_575MeVCM5", "Proton Phi Distribution for Photon Energies of 575pm5MeV (CosTheta0.2-0)", 10, -180, 180);
  Phip585CM5 = new GH1("Phip_585MeVCM5", "Proton Phi Distribution for Photon Energies of 585pm5MeV (CosTheta0.2-0)", 10, -180, 180);
  Phip595CM5 = new GH1("Phip_595MeVCM5", "Proton Phi Distribution for Photon Energies of 595pm5MeV (CosTheta0.2-0)", 10, -180, 180);
  Phip605CM5 = new GH1("Phip_605MeVCM5", "Proton Phi Distribution for Photon Energies of 605pm5MeV (CosTheta0.2-0)", 10, -180, 180);
  Phip615CM5 = new GH1("Phip_615MeVCM5", "Proton Phi Distribution for Photon Energies of 615pm5MeV (CosTheta0.2-0)", 10, -180, 180);
  Phip625CM5 = new GH1("Phip_625MeVCM5", "Proton Phi Distribution for Photon Energies of 625pm5MeV (CosTheta0.2-0)", 10, -180, 180);

  // Proton Phi dists across EGamma bins
  Phip425CM6 = new GH1("Phip_425MeVCM6", "Proton Phi Distribution for Photon Energies of 425pm5MeV (CosTheta0-(-0.2))", 10, -180, 180);
  Phip435CM6 = new GH1("Phip_435MeVCM6", "Proton Phi Distribution for Photon Energies of 435pm5MeV (CosTheta0-(-0.2))", 10, -180, 180);
  Phip445CM6 = new GH1("Phip_445MeVCM6", "Proton Phi Distribution for Photon Energies of 445pm5MeV (CosTheta0-(-0.2))", 10, -180, 180);
  Phip455CM6 = new GH1("Phip_455MeVCM6", "Proton Phi Distribution for Photon Energies of 455pm5MeV (CosTheta0-(-0.2))", 10, -180, 180);
  Phip465CM6 = new GH1("Phip_465MeVCM6", "Proton Phi Distribution for Photon Energies of 465pm5MeV (CosTheta0-(-0.2))", 10, -180, 180);
  Phip475CM6 = new GH1("Phip_475MeVCM6", "Proton Phi Distribution for Photon Energies of 475pm5MeV (CosTheta0-(-0.2))", 10, -180, 180);
  Phip485CM6 = new GH1("Phip_485MeVCM6", "Proton Phi Distribution for Photon Energies of 485pm5MeV (CosTheta0-(-0.2))", 10, -180, 180);
  Phip495CM6 = new GH1("Phip_495MeVCM6", "Proton Phi Distribution for Photon Energies of 495pm5MeV (CosTheta0-(-0.2))", 10, -180, 180);
  Phip505CM6 = new GH1("Phip_505MeVCM6", "Proton Phi Distribution for Photon Energies of 505pm5MeV (CosTheta0-(-0.2))", 10, -180, 180);
  Phip515CM6 = new GH1("Phip_515MeVCM6", "Proton Phi Distribution for Photon Energies of 515pm5MeV (CosTheta0-(-0.2))", 10, -180, 180);
  Phip525CM6 = new GH1("Phip_525MeVCM6", "Proton Phi Distribution for Photon Energies of 525pm5MeV (CosTheta0-(-0.2))", 10, -180, 180);
  Phip535CM6 = new GH1("Phip_535MeVCM6", "Proton Phi Distribution for Photon Energies of 535pm5MeV (CosTheta0-(-0.2))", 10, -180, 180);
  Phip545CM6 = new GH1("Phip_545MeVCM6", "Proton Phi Distribution for Photon Energies of 545pm5MeV (CosTheta0-(-0.2))", 10, -180, 180);
  Phip555CM6 = new GH1("Phip_555MeVCM6", "Proton Phi Distribution for Photon Energies of 555pm5MeV (CosTheta0-(-0.2))", 10, -180, 180);
  Phip565CM6 = new GH1("Phip_565MeVCM6", "Proton Phi Distribution for Photon Energies of 535pm5MeV (CosTheta0-(-0.2))", 10, -180, 180);
  Phip575CM6 = new GH1("Phip_575MeVCM6", "Proton Phi Distribution for Photon Energies of 575pm5MeV (CosTheta0-(-0.2))", 10, -180, 180);
  Phip585CM6 = new GH1("Phip_585MeVCM6", "Proton Phi Distribution for Photon Energies of 585pm5MeV (CosTheta0-(-0.2))", 10, -180, 180);
  Phip595CM6 = new GH1("Phip_595MeVCM6", "Proton Phi Distribution for Photon Energies of 595pm5MeV (CosTheta0-(-0.2))", 10, -180, 180);
  Phip605CM6 = new GH1("Phip_605MeVCM6", "Proton Phi Distribution for Photon Energies of 605pm5MeV (CosTheta0-(-0.2))", 10, -180, 180);
  Phip615CM6 = new GH1("Phip_615MeVCM6", "Proton Phi Distribution for Photon Energies of 615pm5MeV (CosTheta0-(-0.2))", 10, -180, 180);
  Phip625CM6 = new GH1("Phip_625MeVCM6", "Proton Phi Distribution for Photon Energies of 625pm5MeV (CosTheta0-(-0.2))", 10, -180, 180);

  // Proton Phi dists across EGamma bins
  Phip425CM7 = new GH1("Phip_425MeVCM7", "Proton Phi Distribution for Photon Energies of 425pm5MeV (CosTheta-0.2-(-0.4))", 10, -180, 180);
  Phip435CM7 = new GH1("Phip_435MeVCM7", "Proton Phi Distribution for Photon Energies of 435pm5MeV (CosTheta-0.2-(-0.4))", 10, -180, 180);
  Phip445CM7 = new GH1("Phip_445MeVCM7", "Proton Phi Distribution for Photon Energies of 445pm5MeV (CosTheta-0.2-(-0.4))", 10, -180, 180);
  Phip455CM7 = new GH1("Phip_455MeVCM7", "Proton Phi Distribution for Photon Energies of 455pm5MeV (CosTheta-0.2-(-0.4))", 10, -180, 180);
  Phip465CM7 = new GH1("Phip_465MeVCM7", "Proton Phi Distribution for Photon Energies of 465pm5MeV (CosTheta-0.2-(-0.4))", 10, -180, 180);
  Phip475CM7 = new GH1("Phip_475MeVCM7", "Proton Phi Distribution for Photon Energies of 475pm5MeV (CosTheta-0.2-(-0.4))", 10, -180, 180);
  Phip485CM7 = new GH1("Phip_485MeVCM7", "Proton Phi Distribution for Photon Energies of 485pm5MeV (CosTheta-0.2-(-0.4))", 10, -180, 180);
  Phip495CM7 = new GH1("Phip_495MeVCM7", "Proton Phi Distribution for Photon Energies of 495pm5MeV (CosTheta-0.2-(-0.4))", 10, -180, 180);
  Phip505CM7 = new GH1("Phip_505MeVCM7", "Proton Phi Distribution for Photon Energies of 505pm5MeV (CosTheta-0.2-(-0.4))", 10, -180, 180);
  Phip515CM7 = new GH1("Phip_515MeVCM7", "Proton Phi Distribution for Photon Energies of 515pm5MeV (CosTheta-0.2-(-0.4))", 10, -180, 180);
  Phip525CM7 = new GH1("Phip_525MeVCM7", "Proton Phi Distribution for Photon Energies of 525pm5MeV (CosTheta-0.2-(-0.4))", 10, -180, 180);
  Phip535CM7 = new GH1("Phip_535MeVCM7", "Proton Phi Distribution for Photon Energies of 535pm5MeV (CosTheta-0.2-(-0.4))", 10, -180, 180);
  Phip545CM7 = new GH1("Phip_545MeVCM7", "Proton Phi Distribution for Photon Energies of 545pm5MeV (CosTheta-0.2-(-0.4))", 10, -180, 180);
  Phip555CM7 = new GH1("Phip_555MeVCM7", "Proton Phi Distribution for Photon Energies of 555pm5MeV (CosTheta-0.2-(-0.4))", 10, -180, 180);
  Phip565CM7 = new GH1("Phip_565MeVCM7", "Proton Phi Distribution for Photon Energies of 535pm5MeV (CosTheta-0.2-(-0.4))", 10, -180, 180);
  Phip575CM7 = new GH1("Phip_575MeVCM7", "Proton Phi Distribution for Photon Energies of 575pm5MeV (CosTheta-0.2-(-0.4))", 10, -180, 180);
  Phip585CM7 = new GH1("Phip_585MeVCM7", "Proton Phi Distribution for Photon Energies of 585pm5MeV (CosTheta-0.2-(-0.4))", 10, -180, 180);
  Phip595CM7 = new GH1("Phip_595MeVCM7", "Proton Phi Distribution for Photon Energies of 595pm5MeV (CosTheta-0.2-(-0.4))", 10, -180, 180);
  Phip605CM7 = new GH1("Phip_605MeVCM7", "Proton Phi Distribution for Photon Energies of 605pm5MeV (CosTheta-0.2-(-0.4))", 10, -180, 180);
  Phip615CM7 = new GH1("Phip_615MeVCM7", "Proton Phi Distribution for Photon Energies of 615pm5MeV (CosTheta-0.2-(-0.4))", 10, -180, 180);
  Phip625CM7 = new GH1("Phip_625MeVCM7", "Proton Phi Distribution for Photon Energies of 625pm5MeV (CosTheta-0.2-(-0.4))", 10, -180, 180);

  // Proton Phi dists across EGamma bins
  Phip425CM8 = new GH1("Phip_425MeVCM8", "Proton Phi Distribution for Photon Energies of 425pm5MeV (CosTheta-0.4-(-0.6))", 10, -180, 180);
  Phip435CM8 = new GH1("Phip_435MeVCM8", "Proton Phi Distribution for Photon Energies of 435pm5MeV (CosTheta-0.4-(-0.6))", 10, -180, 180);
  Phip445CM8 = new GH1("Phip_445MeVCM8", "Proton Phi Distribution for Photon Energies of 445pm5MeV (CosTheta-0.4-(-0.6))", 10, -180, 180);
  Phip455CM8 = new GH1("Phip_455MeVCM8", "Proton Phi Distribution for Photon Energies of 455pm5MeV (CosTheta-0.4-(-0.6))", 10, -180, 180);
  Phip465CM8 = new GH1("Phip_465MeVCM8", "Proton Phi Distribution for Photon Energies of 465pm5MeV (CosTheta-0.4-(-0.6))", 10, -180, 180);
  Phip475CM8 = new GH1("Phip_475MeVCM8", "Proton Phi Distribution for Photon Energies of 475pm5MeV (CosTheta-0.4-(-0.6))", 10, -180, 180);
  Phip485CM8 = new GH1("Phip_485MeVCM8", "Proton Phi Distribution for Photon Energies of 485pm5MeV (CosTheta-0.4-(-0.6))", 10, -180, 180);
  Phip495CM8 = new GH1("Phip_495MeVCM8", "Proton Phi Distribution for Photon Energies of 495pm5MeV (CosTheta-0.4-(-0.6))", 10, -180, 180);
  Phip505CM8 = new GH1("Phip_505MeVCM8", "Proton Phi Distribution for Photon Energies of 505pm5MeV (CosTheta-0.4-(-0.6))", 10, -180, 180);
  Phip515CM8 = new GH1("Phip_515MeVCM8", "Proton Phi Distribution for Photon Energies of 515pm5MeV (CosTheta-0.4-(-0.6))", 10, -180, 180);
  Phip525CM8 = new GH1("Phip_525MeVCM8", "Proton Phi Distribution for Photon Energies of 525pm5MeV (CosTheta-0.4-(-0.6))", 10, -180, 180);
  Phip535CM8 = new GH1("Phip_535MeVCM8", "Proton Phi Distribution for Photon Energies of 535pm5MeV (CosTheta-0.4-(-0.6))", 10, -180, 180);
  Phip545CM8 = new GH1("Phip_545MeVCM8", "Proton Phi Distribution for Photon Energies of 545pm5MeV (CosTheta-0.4-(-0.6))", 10, -180, 180);
  Phip555CM8 = new GH1("Phip_555MeVCM8", "Proton Phi Distribution for Photon Energies of 555pm5MeV (CosTheta-0.4-(-0.6))", 10, -180, 180);
  Phip565CM8 = new GH1("Phip_565MeVCM8", "Proton Phi Distribution for Photon Energies of 535pm5MeV (CosTheta-0.4-(-0.6))", 10, -180, 180);
  Phip575CM8 = new GH1("Phip_575MeVCM8", "Proton Phi Distribution for Photon Energies of 575pm5MeV (CosTheta-0.4-(-0.6))", 10, -180, 180);
  Phip585CM8 = new GH1("Phip_585MeVCM8", "Proton Phi Distribution for Photon Energies of 585pm5MeV (CosTheta-0.4-(-0.6))", 10, -180, 180);
  Phip595CM8 = new GH1("Phip_595MeVCM8", "Proton Phi Distribution for Photon Energies of 595pm5MeV (CosTheta-0.4-(-0.6))", 10, -180, 180);
  Phip605CM8 = new GH1("Phip_605MeVCM8", "Proton Phi Distribution for Photon Energies of 605pm5MeV (CosTheta-0.4-(-0.6))", 10, -180, 180);
  Phip615CM8 = new GH1("Phip_615MeVCM8", "Proton Phi Distribution for Photon Energies of 615pm5MeV (CosTheta-0.4-(-0.6))", 10, -180, 180);
  Phip625CM8 = new GH1("Phip_625MeVCM8", "Proton Phi Distribution for Photon Energies of 625pm5MeV (CosTheta-0.4-(-0.6))", 10, -180, 180);

  // Proton Phi dists across EGamma bins
  Phip425CM9 = new GH1("Phip_425MeVCM9", "Proton Phi Distribution for Photon Energies of 425pm5MeV (CosTheta-0.6-(-0.8))", 10, -180, 180);
  Phip435CM9 = new GH1("Phip_435MeVCM9", "Proton Phi Distribution for Photon Energies of 435pm5MeV (CosTheta-0.6-(-0.8))", 10, -180, 180);
  Phip445CM9 = new GH1("Phip_445MeVCM9", "Proton Phi Distribution for Photon Energies of 445pm5MeV (CosTheta-0.6-(-0.8))", 10, -180, 180);
  Phip455CM9 = new GH1("Phip_455MeVCM9", "Proton Phi Distribution for Photon Energies of 455pm5MeV (CosTheta-0.6-(-0.8))", 10, -180, 180);
  Phip465CM9 = new GH1("Phip_465MeVCM9", "Proton Phi Distribution for Photon Energies of 465pm5MeV (CosTheta-0.6-(-0.8))", 10, -180, 180);
  Phip475CM9 = new GH1("Phip_475MeVCM9", "Proton Phi Distribution for Photon Energies of 475pm5MeV (CosTheta-0.6-(-0.8))", 10, -180, 180);
  Phip485CM9 = new GH1("Phip_485MeVCM9", "Proton Phi Distribution for Photon Energies of 485pm5MeV (CosTheta-0.6-(-0.8))", 10, -180, 180);
  Phip495CM9 = new GH1("Phip_495MeVCM9", "Proton Phi Distribution for Photon Energies of 495pm5MeV (CosTheta-0.6-(-0.8))", 10, -180, 180);
  Phip505CM9 = new GH1("Phip_505MeVCM9", "Proton Phi Distribution for Photon Energies of 505pm5MeV (CosTheta-0.6-(-0.8))", 10, -180, 180);
  Phip515CM9 = new GH1("Phip_515MeVCM9", "Proton Phi Distribution for Photon Energies of 515pm5MeV (CosTheta-0.6-(-0.8))", 10, -180, 180);
  Phip525CM9 = new GH1("Phip_525MeVCM9", "Proton Phi Distribution for Photon Energies of 525pm5MeV (CosTheta-0.6-(-0.8))", 10, -180, 180);
  Phip535CM9 = new GH1("Phip_535MeVCM9", "Proton Phi Distribution for Photon Energies of 535pm5MeV (CosTheta-0.6-(-0.8))", 10, -180, 180);
  Phip545CM9 = new GH1("Phip_545MeVCM9", "Proton Phi Distribution for Photon Energies of 545pm5MeV (CosTheta-0.6-(-0.8))", 10, -180, 180);
  Phip555CM9 = new GH1("Phip_555MeVCM9", "Proton Phi Distribution for Photon Energies of 555pm5MeV (CosTheta-0.6-(-0.8))", 10, -180, 180);
  Phip565CM9 = new GH1("Phip_565MeVCM9", "Proton Phi Distribution for Photon Energies of 535pm5MeV (CosTheta-0.6-(-0.8))", 10, -180, 180);
  Phip575CM9 = new GH1("Phip_575MeVCM9", "Proton Phi Distribution for Photon Energies of 575pm5MeV (CosTheta-0.6-(-0.8))", 10, -180, 180);
  Phip585CM9 = new GH1("Phip_585MeVCM9", "Proton Phi Distribution for Photon Energies of 585pm5MeV (CosTheta-0.6-(-0.8))", 10, -180, 180);
  Phip595CM9 = new GH1("Phip_595MeVCM9", "Proton Phi Distribution for Photon Energies of 595pm5MeV (CosTheta-0.6-(-0.8))", 10, -180, 180);
  Phip605CM9 = new GH1("Phip_605MeVCM9", "Proton Phi Distribution for Photon Energies of 605pm5MeV (CosTheta-0.6-(-0.8))", 10, -180, 180);
  Phip615CM9 = new GH1("Phip_615MeVCM9", "Proton Phi Distribution for Photon Energies of 615pm5MeV (CosTheta-0.6-(-0.8))", 10, -180, 180);
  Phip625CM9 = new GH1("Phip_625MeVCM9", "Proton Phi Distribution for Photon Energies of 625pm5MeV (CosTheta-0.6-(-0.8))", 10, -180, 180);

  // Proton Phi dists across EGamma bins
  Phip425CM10 = new GH1("Phip_425MeVCM10", "Proton Phi Distribution for Photon Energies of 425pm5MeV (CosTheta-0.8-(-1))", 10, -180, 180);
  Phip435CM10 = new GH1("Phip_435MeVCM10", "Proton Phi Distribution for Photon Energies of 435pm5MeV (CosTheta-0.8-(-1))", 10, -180, 180);
  Phip445CM10 = new GH1("Phip_445MeVCM10", "Proton Phi Distribution for Photon Energies of 445pm5MeV (CosTheta-0.8-(-1))", 10, -180, 180);
  Phip455CM10 = new GH1("Phip_455MeVCM10", "Proton Phi Distribution for Photon Energies of 455pm5MeV (CosTheta-0.8-(-1))", 10, -180, 180);
  Phip465CM10 = new GH1("Phip_465MeVCM10", "Proton Phi Distribution for Photon Energies of 465pm5MeV (CosTheta-0.8-(-1))", 10, -180, 180);
  Phip475CM10 = new GH1("Phip_475MeVCM10", "Proton Phi Distribution for Photon Energies of 475pm5MeV (CosTheta-0.8-(-1))", 10, -180, 180);
  Phip485CM10 = new GH1("Phip_485MeVCM10", "Proton Phi Distribution for Photon Energies of 485pm5MeV (CosTheta-0.8-(-1))", 10, -180, 180);
  Phip495CM10 = new GH1("Phip_495MeVCM10", "Proton Phi Distribution for Photon Energies of 495pm5MeV (CosTheta-0.8-(-1))", 10, -180, 180);
  Phip505CM10 = new GH1("Phip_505MeVCM10", "Proton Phi Distribution for Photon Energies of 505pm5MeV (CosTheta-0.8-(-1))", 10, -180, 180);
  Phip515CM10 = new GH1("Phip_515MeVCM10", "Proton Phi Distribution for Photon Energies of 515pm5MeV (CosTheta-0.8-(-1))", 10, -180, 180);
  Phip525CM10 = new GH1("Phip_525MeVCM10", "Proton Phi Distribution for Photon Energies of 525pm5MeV (CosTheta-0.8-(-1))", 10, -180, 180);
  Phip535CM10 = new GH1("Phip_535MeVCM10", "Proton Phi Distribution for Photon Energies of 535pm5MeV (CosTheta-0.8-(-1))", 10, -180, 180);
  Phip545CM10 = new GH1("Phip_545MeVCM10", "Proton Phi Distribution for Photon Energies of 545pm5MeV (CosTheta-0.8-(-1))", 10, -180, 180);
  Phip555CM10 = new GH1("Phip_555MeVCM10", "Proton Phi Distribution for Photon Energies of 555pm5MeV (CosTheta-0.8-(-1))", 10, -180, 180);
  Phip565CM10 = new GH1("Phip_565MeVCM10", "Proton Phi Distribution for Photon Energies of 535pm5MeV (CosTheta-0.8-(-1))", 10, -180, 180);
  Phip575CM10 = new GH1("Phip_575MeVCM10", "Proton Phi Distribution for Photon Energies of 575pm5MeV (CosTheta-0.8-(-1))", 10, -180, 180);
  Phip585CM10 = new GH1("Phip_585MeVCM10", "Proton Phi Distribution for Photon Energies of 585pm5MeV (CosTheta-0.8-(-1))", 10, -180, 180);
  Phip595CM10 = new GH1("Phip_595MeVCM10", "Proton Phi Distribution for Photon Energies of 595pm5MeV (CosTheta-0.8-(-1))", 10, -180, 180);
  Phip605CM10 = new GH1("Phip_605MeVCM10", "Proton Phi Distribution for Photon Energies of 605pm5MeV (CosTheta-0.8-(-1))", 10, -180, 180);
  Phip615CM10 = new GH1("Phip_615MeVCM10", "Proton Phi Distribution for Photon Energies of 615pm5MeV (CosTheta-0.8-(-1))", 10, -180, 180);
  Phip625CM10 = new GH1("Phip_625MeVCM10", "Proton Phi Distribution for Photon Energies of 625pm5MeV (CosTheta-0.8-(-1))", 10, -180, 180);

  ThetaRecPiDiff = new GH1 ("ThetaRecPiDiff", "Difference between ThetaPiRec and Thetan", 200, 0, 180);
  ThetanThetaRecPi = new GH2 ("ThetanThetaRecPi", "Thetan vs ThetaPiRec", 100, 0, 180, 100, 0, 180);
  ThetanThetaRecPiDiff = new GH2 ("ThetanThetaRecPiDiff", "Thetan vs (ThetaPiRec - Thetan)", 100, 0, 180, 100, 0, 180);
  ThetaRecPDiff = new GH1 ("ThetaRecPDiff", "Difference between ThetaPRec and Thetan", 200, 0, 180);
  ThetanThetaRecP = new GH2 ("ThetanThetaRecP", "Thetan vs ThetaPRec", 100, 0, 180, 100, 0, 180);
  ThetanThetaRecPDiff = new GH2 ("ThetanThetaRecPDiff", "Thetan vs (ThetaPRec - Thetan)", 100, 0, 180, 100, 0, 180);

  E_dE = new GH2 ("E_dE", "EdE Plot With E Loss Adjustment", 100, 0, 500, 100, 0, 5);
  KinEp_dE = new GH2 ("KinEp_dE", "KinEpdE Plot", 100, 0, 500, 100, 0, 5);
  E_KinEp = new GH2 ("E_KinEp", "Kinematic Energy of Proton as a function of CB energy", 100, 0, 500, 100, 0, 500);

}

void PNeutPol_Polarimeter_Lin_NoScatt::FillHists()
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

    if ( 420 < EGamma && EGamma < 430) {

        if(1 > CosThetapCM && CosThetapCM > 0.8 ){
            Phip425CM1->Fill(WCPhip, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.6){
            Phip425CM2->Fill(WCPhip, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.4){
            Phip425CM3->Fill(WCPhip, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.2){
            Phip425CM4->Fill(WCPhip, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0){
            Phip425CM5->Fill(WCPhip, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.2){
            Phip425CM6->Fill(WCPhip, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.4){
            Phip425CM7->Fill(WCPhip, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.6){
            Phip425CM8->Fill(WCPhip, TaggerTime);
        }

        else if(-0.6> CosThetapCM && CosThetapCM > -0.8){
            Phip425CM9->Fill(WCPhip, TaggerTime);
        }

        else if(-0.8> CosThetapCM && CosThetapCM > -1){
            Phip425CM10->Fill(WCPhip, TaggerTime);
        }
    }

    else if ( 430 < EGamma && EGamma < 440) {

        if(1 > CosThetapCM && CosThetapCM > 0.8 ){
            Phip435CM1->Fill(WCPhip, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.6){
            Phip435CM2->Fill(WCPhip, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.4){
            Phip435CM3->Fill(WCPhip, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.2){
            Phip435CM4->Fill(WCPhip, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0){
            Phip435CM5->Fill(WCPhip, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.2){
            Phip435CM6->Fill(WCPhip, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.4){
            Phip435CM7->Fill(WCPhip, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.6){
            Phip435CM8->Fill(WCPhip, TaggerTime);
        }

        else if(-0.6> CosThetapCM && CosThetapCM > -0.8){
            Phip435CM9->Fill(WCPhip, TaggerTime);
        }

        else if(-0.8> CosThetapCM && CosThetapCM > -1){
            Phip435CM10->Fill(WCPhip, TaggerTime);
        }
    }

    else if ( 440 < EGamma && EGamma < 450) {

        if(1 > CosThetapCM && CosThetapCM > 0.8 ){
            Phip445CM1->Fill(WCPhip, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.6){
            Phip445CM2->Fill(WCPhip, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.4){
            Phip445CM3->Fill(WCPhip, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.2){
            Phip445CM4->Fill(WCPhip, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0){
            Phip445CM5->Fill(WCPhip, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.2){
            Phip445CM6->Fill(WCPhip, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.4){
            Phip445CM7->Fill(WCPhip, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.6){
            Phip445CM8->Fill(WCPhip, TaggerTime);
        }

        else if(-0.6> CosThetapCM && CosThetapCM > -0.8){
            Phip445CM9->Fill(WCPhip, TaggerTime);
        }

        else if(-0.8> CosThetapCM && CosThetapCM > -1){
            Phip445CM10->Fill(WCPhip, TaggerTime);
        }
    }

    else if ( 450 < EGamma && EGamma < 460) {

        if(1 > CosThetapCM && CosThetapCM > 0.8 ){
            Phip455CM1->Fill(WCPhip, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.6){
            Phip455CM2->Fill(WCPhip, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.4){
            Phip455CM3->Fill(WCPhip, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.2){
            Phip455CM4->Fill(WCPhip, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0){
            Phip455CM5->Fill(WCPhip, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.2){
            Phip455CM6->Fill(WCPhip, TaggerTime);
        }
        else if(-0.2 > CosThetapCM && CosThetapCM > -0.4){
            Phip455CM7->Fill(WCPhip, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.6){
            Phip455CM8->Fill(WCPhip, TaggerTime);
        }

        else if(-0.6> CosThetapCM && CosThetapCM > -0.8){
            Phip455CM9->Fill(WCPhip, TaggerTime);
        }

        else if(-0.8> CosThetapCM && CosThetapCM > -1){
            Phip455CM10->Fill(WCPhip, TaggerTime);
        }
    }

    else if ( 460 < EGamma && EGamma < 470) {

        if(1 > CosThetapCM && CosThetapCM > 0.8 ){
            Phip465CM1->Fill(WCPhip, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.6){
            Phip465CM2->Fill(WCPhip, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.4){
            Phip465CM3->Fill(WCPhip, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.2){
            Phip465CM4->Fill(WCPhip, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0){
            Phip465CM5->Fill(WCPhip, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.2){
            Phip465CM6->Fill(WCPhip, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.4){
            Phip465CM7->Fill(WCPhip, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.6){
            Phip465CM8->Fill(WCPhip, TaggerTime);
        }

        else if(-0.6> CosThetapCM && CosThetapCM > -0.8){
            Phip465CM9->Fill(WCPhip, TaggerTime);
        }

        else if(-0.8> CosThetapCM && CosThetapCM > -1){
            Phip465CM10->Fill(WCPhip, TaggerTime);
        }
    }

    else if ( 470 < EGamma && EGamma < 480) {

        if(1 > CosThetapCM && CosThetapCM > 0.8 ){
            Phip475CM1->Fill(WCPhip, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.6){
            Phip475CM2->Fill(WCPhip, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.4){
            Phip475CM3->Fill(WCPhip, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.2){
            Phip475CM4->Fill(WCPhip, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0){
            Phip475CM5->Fill(WCPhip, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.2){
            Phip475CM6->Fill(WCPhip, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.4){
            Phip475CM7->Fill(WCPhip, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.6){
            Phip475CM8->Fill(WCPhip, TaggerTime);
        }

        else if(-0.6> CosThetapCM && CosThetapCM > -0.8){
            Phip475CM9->Fill(WCPhip, TaggerTime);
        }

        else if(-0.8> CosThetapCM && CosThetapCM > -1){
            Phip475CM10->Fill(WCPhip, TaggerTime);
        }
    }

    else if ( 480 < EGamma && EGamma < 490) {

        if(1 > CosThetapCM && CosThetapCM > 0.8 ){
            Phip485CM1->Fill(WCPhip, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.6){
            Phip485CM2->Fill(WCPhip, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.4){
            Phip485CM3->Fill(WCPhip, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.2){
            Phip485CM4->Fill(WCPhip, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0){
            Phip485CM5->Fill(WCPhip, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.2){
            Phip485CM6->Fill(WCPhip, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.4){
            Phip485CM7->Fill(WCPhip, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.6){
            Phip485CM8->Fill(WCPhip, TaggerTime);
        }

        else if(-0.6> CosThetapCM && CosThetapCM > -0.8){
            Phip485CM9->Fill(WCPhip, TaggerTime);
        }

        else if(-0.8> CosThetapCM && CosThetapCM > -1){
            Phip485CM10->Fill(WCPhip, TaggerTime);
        }
    }

    else if ( 490 < EGamma && EGamma < 500) {

        if(1 > CosThetapCM && CosThetapCM > 0.8 ){
            Phip495CM1->Fill(WCPhip, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.6){
            Phip495CM2->Fill(WCPhip, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.4){
            Phip495CM3->Fill(WCPhip, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.2){
            Phip495CM4->Fill(WCPhip, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0){
            Phip495CM5->Fill(WCPhip, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.2){
            Phip495CM6->Fill(WCPhip, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.4){
            Phip495CM7->Fill(WCPhip, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.6){
            Phip495CM8->Fill(WCPhip, TaggerTime);
        }

        else if(-0.6> CosThetapCM && CosThetapCM > -0.8){
            Phip495CM9->Fill(WCPhip, TaggerTime);
        }

        else if(-0.8> CosThetapCM && CosThetapCM > -1){
            Phip495CM10->Fill(WCPhip, TaggerTime);
        }
    }

    else if ( 500 < EGamma && EGamma < 510) {

        if(1 > CosThetapCM && CosThetapCM > 0.8 ){
            Phip505CM1->Fill(WCPhip, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.6){
            Phip505CM2->Fill(WCPhip, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.4){
            Phip505CM3->Fill(WCPhip, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.2){
            Phip505CM4->Fill(WCPhip, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0){
            Phip505CM5->Fill(WCPhip, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.2){
            Phip505CM6->Fill(WCPhip, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.4){
            Phip505CM7->Fill(WCPhip, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.6){
            Phip505CM8->Fill(WCPhip, TaggerTime);
        }

        else if(-0.6> CosThetapCM && CosThetapCM > -0.8){
            Phip505CM9->Fill(WCPhip, TaggerTime);
        }

        else if(-0.8> CosThetapCM && CosThetapCM > -1){
            Phip505CM10->Fill(WCPhip, TaggerTime);
        }
    }

    else if ( 510 < EGamma && EGamma < 520) {

        if(1 > CosThetapCM && CosThetapCM > 0.8 ){
            Phip515CM1->Fill(WCPhip, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.6){
            Phip515CM2->Fill(WCPhip, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.4){
            Phip515CM3->Fill(WCPhip, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.2){
            Phip515CM4->Fill(WCPhip, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0){
            Phip515CM5->Fill(WCPhip, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.2){
            Phip515CM6->Fill(WCPhip, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.4){
            Phip515CM7->Fill(WCPhip, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.6){
            Phip515CM8->Fill(WCPhip, TaggerTime);
        }

        else if(-0.6> CosThetapCM && CosThetapCM > -0.8){
            Phip515CM9->Fill(WCPhip, TaggerTime);
        }

        else if(-0.8> CosThetapCM && CosThetapCM > -1){
            Phip515CM10->Fill(WCPhip, TaggerTime);
        }
    }

    else if ( 520 < EGamma && EGamma < 530) {

        if(1 > CosThetapCM && CosThetapCM > 0.8 ){
            Phip525CM1->Fill(WCPhip, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.6){
            Phip525CM2->Fill(WCPhip, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.4){
            Phip525CM3->Fill(WCPhip, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.2){
            Phip525CM4->Fill(WCPhip, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0){
            Phip525CM5->Fill(WCPhip, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.2){
            Phip525CM6->Fill(WCPhip, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.4){
            Phip525CM7->Fill(WCPhip, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.6){
            Phip525CM8->Fill(WCPhip, TaggerTime);
        }

        else if(-0.6> CosThetapCM && CosThetapCM > -0.8){
            Phip525CM9->Fill(WCPhip, TaggerTime);
        }

        else if(-0.8> CosThetapCM && CosThetapCM > -1){
            Phip525CM10->Fill(WCPhip, TaggerTime);
        }
    }

    else if ( 530 < EGamma && EGamma < 540) {

        if(1 > CosThetapCM && CosThetapCM > 0.8 ){
            Phip535CM1->Fill(WCPhip, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.6){
            Phip535CM2->Fill(WCPhip, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.4){
            Phip535CM3->Fill(WCPhip, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.2){
            Phip535CM4->Fill(WCPhip, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0){
            Phip535CM5->Fill(WCPhip, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.2){
            Phip535CM6->Fill(WCPhip, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.4){
            Phip535CM7->Fill(WCPhip, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.6){
            Phip535CM8->Fill(WCPhip, TaggerTime);
        }

        else if(-0.6> CosThetapCM && CosThetapCM > -0.8){
            Phip535CM9->Fill(WCPhip, TaggerTime);
        }

        else if(-0.8> CosThetapCM && CosThetapCM > -1){
            Phip535CM10->Fill(WCPhip, TaggerTime);
        }
    }

    else if ( 540 < EGamma && EGamma < 550) {

        if(1 > CosThetapCM && CosThetapCM > 0.8 ){
            Phip545CM1->Fill(WCPhip, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.6){
            Phip545CM2->Fill(WCPhip, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.4){
            Phip545CM3->Fill(WCPhip, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.2){
            Phip545CM4->Fill(WCPhip, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0){
            Phip545CM5->Fill(WCPhip, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.2){
            Phip545CM6->Fill(WCPhip, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.4){
            Phip545CM7->Fill(WCPhip, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.6){
            Phip545CM8->Fill(WCPhip, TaggerTime);
        }

        else if(-0.6> CosThetapCM && CosThetapCM > -0.8){
            Phip545CM9->Fill(WCPhip, TaggerTime);
        }

        else if(-0.8> CosThetapCM && CosThetapCM > -1){
            Phip545CM10->Fill(WCPhip, TaggerTime);
        }
    }

    else if ( 550 < EGamma && EGamma < 560) {

        if(1 > CosThetapCM && CosThetapCM > 0.8 ){
            Phip555CM1->Fill(WCPhip, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.6){
            Phip555CM2->Fill(WCPhip, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.4){
            Phip555CM3->Fill(WCPhip, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.2){
            Phip555CM4->Fill(WCPhip, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0){
            Phip555CM5->Fill(WCPhip, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.2){
            Phip555CM6->Fill(WCPhip, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.4){
            Phip555CM7->Fill(WCPhip, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.6){
            Phip555CM8->Fill(WCPhip, TaggerTime);
        }

        else if(-0.6> CosThetapCM && CosThetapCM > -0.8){
            Phip555CM9->Fill(WCPhip, TaggerTime);
        }

        else if(-0.8> CosThetapCM && CosThetapCM > -1){
            Phip555CM10->Fill(WCPhip, TaggerTime);
        }
    }

    else if ( 560 < EGamma && EGamma < 570) {

        if(1 > CosThetapCM && CosThetapCM > 0.8 ){
            Phip565CM1->Fill(WCPhip, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.6){
            Phip565CM2->Fill(WCPhip, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.4){
            Phip565CM3->Fill(WCPhip, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.2){
            Phip565CM4->Fill(WCPhip, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0){
            Phip565CM5->Fill(WCPhip, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.2){
            Phip565CM6->Fill(WCPhip, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.4){
            Phip565CM7->Fill(WCPhip, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.6){
            Phip565CM8->Fill(WCPhip, TaggerTime);
        }

        else if(-0.6> CosThetapCM && CosThetapCM > -0.8){
            Phip565CM9->Fill(WCPhip, TaggerTime);
        }

        else if(-0.8> CosThetapCM && CosThetapCM > -1){
            Phip565CM10->Fill(WCPhip, TaggerTime);
        }
    }

    else if ( 570 < EGamma && EGamma < 580) {

        if(1 > CosThetapCM && CosThetapCM > 0.8 ){
            Phip575CM1->Fill(WCPhip, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.6){
            Phip575CM2->Fill(WCPhip, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.4){
            Phip575CM3->Fill(WCPhip, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.2){
            Phip575CM4->Fill(WCPhip, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0){
            Phip575CM5->Fill(WCPhip, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.2){
            Phip575CM6->Fill(WCPhip, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.4){
            Phip575CM7->Fill(WCPhip, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.6){
            Phip575CM8->Fill(WCPhip, TaggerTime);
        }

        else if(-0.6> CosThetapCM && CosThetapCM > -0.8){
            Phip575CM9->Fill(WCPhip, TaggerTime);
        }

        else if(-0.8> CosThetapCM && CosThetapCM > -1){
            Phip575CM10->Fill(WCPhip, TaggerTime);
        }
    }

    else if ( 580 < EGamma && EGamma < 590) {

        if(1 > CosThetapCM && CosThetapCM > 0.8 ){
            Phip585CM1->Fill(WCPhip, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.6){
            Phip585CM2->Fill(WCPhip, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.4){
            Phip585CM3->Fill(WCPhip, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.2){
            Phip585CM4->Fill(WCPhip, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0){
            Phip585CM5->Fill(WCPhip, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.2){
            Phip585CM6->Fill(WCPhip, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.4){
            Phip585CM7->Fill(WCPhip, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.6){
            Phip585CM8->Fill(WCPhip, TaggerTime);
        }

        else if(-0.6> CosThetapCM && CosThetapCM > -0.8){
            Phip585CM9->Fill(WCPhip, TaggerTime);
        }

        else if(-0.8> CosThetapCM && CosThetapCM > -1){
            Phip585CM10->Fill(WCPhip, TaggerTime);
        }
    }

    else if ( 590 < EGamma && EGamma < 600) {

        if(1 > CosThetapCM && CosThetapCM > 0.8 ){
            Phip595CM1->Fill(WCPhip, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.6){
            Phip595CM2->Fill(WCPhip, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.4){
            Phip595CM3->Fill(WCPhip, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.2){
            Phip595CM4->Fill(WCPhip, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0){
            Phip595CM5->Fill(WCPhip, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.2){
            Phip595CM6->Fill(WCPhip, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.4){
            Phip595CM7->Fill(WCPhip, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.6){
            Phip595CM8->Fill(WCPhip, TaggerTime);
        }

        else if(-0.6> CosThetapCM && CosThetapCM > -0.8){
            Phip595CM9->Fill(WCPhip, TaggerTime);
        }

        else if(-0.8> CosThetapCM && CosThetapCM > -1){
            Phip595CM10->Fill(WCPhip, TaggerTime);
        }
    }

    else if ( 600 < EGamma && EGamma < 610) {

        if(1 > CosThetapCM && CosThetapCM > 0.8 ){
            Phip605CM1->Fill(WCPhip, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.6){
            Phip605CM2->Fill(WCPhip, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.4){
            Phip605CM3->Fill(WCPhip, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.2){
            Phip605CM4->Fill(WCPhip, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0){
            Phip605CM5->Fill(WCPhip, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.2){
            Phip605CM6->Fill(WCPhip, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.4){
            Phip605CM7->Fill(WCPhip, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.6){
            Phip605CM8->Fill(WCPhip, TaggerTime);
        }

        else if(-0.6> CosThetapCM && CosThetapCM > -0.8){
            Phip605CM9->Fill(WCPhip, TaggerTime);
        }

        else if(-0.8> CosThetapCM && CosThetapCM > -1){
            Phip605CM10->Fill(WCPhip, TaggerTime);
        }
    }

    else if ( 610 < EGamma && EGamma < 620) {

        if(1 > CosThetapCM && CosThetapCM > 0.8 ){
            Phip615CM1->Fill(WCPhip, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.6){
            Phip615CM2->Fill(WCPhip, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.4){
            Phip615CM3->Fill(WCPhip, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.2){
            Phip615CM4->Fill(WCPhip, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0){
            Phip615CM5->Fill(WCPhip, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.2){
            Phip615CM6->Fill(WCPhip, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.4){
            Phip615CM7->Fill(WCPhip, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.6){
            Phip615CM8->Fill(WCPhip, TaggerTime);
        }

        else if(-0.6> CosThetapCM && CosThetapCM > -0.8){
            Phip615CM9->Fill(WCPhip, TaggerTime);
        }

        else if(-0.8> CosThetapCM && CosThetapCM > -1){
            Phip615CM10->Fill(WCPhip, TaggerTime);
        }
    }
}

Bool_t	PNeutPol_Polarimeter_Lin_NoScatt::Write(){
  // Write all GH1's easily

  GTreeManager::Write();
}
