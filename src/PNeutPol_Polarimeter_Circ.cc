// GoAT Physics analysis to identify neutrons from deuterium photodisintegration
// Various properties of neutrons/protons identified plotted in histograms
// Main aim is to determine spin polarisation of neutrons
// For use on circularly polarised data files

#include "PNeutPol_Polarimeter_Circ.h"

PNeutPol_Polarimeter_Circ::~PNeutPol_Polarimeter_Circ()
{
}

Bool_t	PNeutPol_Polarimeter_Circ::Init()
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

Bool_t	PNeutPol_Polarimeter_Circ::Start()
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

void	PNeutPol_Polarimeter_Circ::ProcessEvent()
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

  //if( Zp > 60 || Zp < -60) return;
  //if( Zp > 200 || Zp < 150) return; // Select out windows

  EventCounterZCut++;

  //if ( PhiWCDiff > 195 || PhiWCDiff < 165) return;

  EventCounterCoplanarCut++;

  for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
  {

    TaggerTime = GetTagger()->GetTaggedTime(j); // Get tagged time for event
    EGamma = (GetTagger()->GetTaggedEnergy(j)); // Get Photon energy for event
    Gamma = TLorentzVector (0., 0., EGamma , EGamma); // 4-Vector of Photon beam
    BeamHelicity = GetTrigger()->GetHelicity();

    Thetap = GVp3.Theta()*TMath::RadToDeg(); // Lab frame angles for proton/neutron
    Phip = GVp3.Phi()*TMath::RadToDeg();
    Thetan = GVn3.Theta()*TMath::RadToDeg();
    Phin = GVn3.Phi()*TMath::RadToDeg();

    EpCorr = EpPolCorrect(Ep, WCThetap);
    EpDiff = abs(EpCorr - Ep);

    // Gamma(d,p)n
    KinEp = CalcKinEnergy(WCThetap, EGamma, Md, 0., Mp, Mn); // Calculate kin E of proton assuming pn production
    RecKinProton = CProton4VectorKin(KinEp, WCThetapRad, WCPhipRad);
    RecKinNeutron = CNeutron4VectorKin(RecKinProton);
    ThetanRec = (RecKinNeutron.Theta()) * TMath::RadToDeg();
    WCZnRec = 72/tan(RecKinNeutron.Theta());

    // Gamma(n,p)Pi (P detected correct)
    // Assume proton track is proton and "neutron" track is from charged pion
    KinEpPi = CalcKinEnergy(WCThetap, EGamma, Mn, 0, Mp, Mpi); // Calculate kin E of proton assuming g(n, p) pi
    RecKinProtonPi = CProton4VectorKin(KinEpPi, WCThetapRad, WCPhipRad); // Get Proton 4 vector from calculated kin E
    RecKinPion = CPion4VectorKin(RecKinProtonPi); // Get Pion 4 vector from 4 momenta conservation
    ThetaPiRec = (RecKinPion.Theta())*TMath::RadToDeg();
    PhiPiRec = (RecKinPion.Phi())*TMath::RadToDeg();
    ThetaPiRecDiff = abs(ThetaPiRec - Thetan);

    // Gamma(n,p)Pi (Pion detected correct)
    // Assume proton track is pion and "neutron" track is from proton
    KinPi = CalcKinEnergy(WCThetap, EGamma, Mn, 0, Mpi, Mp); // Calculate kin E of pion
    RecKinPionP = CProton4VectorKin(KinPi, WCThetapRad, WCPhipRad); // Get Pion 4 vector from calculated kinE
    RecKinPPi = CPion4VectorKin(RecKinPionP); // Get Proton 4 vector from 4 momenta conservation
    ThetapRec = (RecKinPPi.Theta())*TMath::RadToDeg();
    PhipRec = (RecKinPPi.Phi())*TMath::RadToDeg();
    ThetapRecDiff = abs (ThetapRec - Thetan);

    KinEDiff = KinEp - EpCorr;

    RecProtonEpCorr = CProton4VectorKin(EpCorr, WCThetapRad, WCPhipRad);
    RecNeutronEpCorr = CNeutron4VectorKin(RecProtonEpCorr);
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
    //if (ScattPhi > -165 && ScattPhi < 165) continue;
    //if (abs(KinEDiff) > 100) continue; // If difference between CB energy and calculated Energy for proton > 100MeV continue

    //k++;
    FillHists(); // Fill histograms with data generated
  }
}

void	PNeutPol_Polarimeter_Circ::ProcessScalerRead()
{
	// Fill Tagger Scalers // Currently this seems to fill the file with loads of "TaggerAccScal" histograms
	//FillScalers(GetTC_scaler_min(),GetTC_scaler_max(),TaggerAccScal); // Don't know if these are needed so cut out for now
}

TCutG*	PNeutPol_Polarimeter_Circ::OpenCutFile(Char_t* filename, Char_t* cutname)
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

TLorentzVector PNeutPol_Polarimeter_Circ::CProton4VectorKin(Double_t KinE, Double_t Theta, Double_t Phi)
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

TLorentzVector PNeutPol_Polarimeter_Circ::CNeutron4VectorKin(TLorentzVector ProtonKinVector)
{

    N4Vect = (Gamma + Deut) - ProtonKinVector;

    return N4Vect;
}

TLorentzVector PNeutPol_Polarimeter_Circ::CPion4VectorKin(TLorentzVector ProtonKinVector)
{

    Pi4Vect = (Gamma + Neut) - ProtonKinVector;

    return Pi4Vect;
}

PNeutPol_Polarimeter_Circ::PNeutPol_Polarimeter_Circ() // Define a load of histograms to fill
{
  time = new TH1D("time", 	"time", 	1400, -700, 700);
  time_cut = new TH1D("time_cut", 	"time_cut", 	1400, -700, 700);

  Eg = new GH1( "Eg", "Photon Energy Distribution", 200, 100, 1600 );
  WCPhiDifference = new GH1 ("WCPhiDifference", "WC Phi Difference Between p and n", 180, 0, 360);
  EpKin = new GH1 ("EpKin", "Ep Calculated from Ep/Thetap", 100, 0, 500);
  EpCorrected = new GH1 ("EpCorrected", "Ep Corrected for Energy Loss in Polarimeter ", 100, 0, 500);
  OAngle = new GH1 ("OAngle", "Opening Angle between P and N Vectors", 180, 0, 180);
  WCZnRecon = new GH1 ("WCZnRecon", "WCZ Hit Position from Reconstructed n Vector", 200, 0, 400);
  WCZnRecon = new GH1 ("WCZnRecon", "WCZ Hit Position from Reconstructed n Vector", 200, 0, 400);

  ThetaSc =  new GH1( "Theta_Scattered", "Scattered Proton Theta Distribution in Rotated Frame", 180, 0, 180 );
  PhiSc = new GH1( "Phi_Scattered", "Scattered Proton Phi Distribution in Rotated Frame", 90, -180, 180 );

  EpKinEpCorrDiff = new GH1("EpKinEpCorrDiff", "Difference Between EpKin and EpCorr", 300, -300, 300);
  EpEpCorrDiff = new GH1("EpEpCorrDiff", "Difference Between Ep and EpCorr", 200, 0, 200);
  MMpEpCorrected = new GH1 ("MMpEpCorrected", "Missing mass seen by Proton (E Loss Corrected)", 400, 0, 2000);

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
  PhiSc275 = new GH1( "Phi_Scattered_275MeV", "Scattered Proton Phi Distribution in Rotated Frame for Photon Energies of 275pm25MeV", 2, -180, 180);
  PhiSc325 = new GH1( "Phi_Scattered_325MeV", "Scattered Proton Phi Distribution in Rotated Frame for Photon Energies of 325pm25MeV", 2, -180, 180);
  PhiSc375 = new GH1( "Phi_Scattered_375MeV", "Scattered Proton Phi Distribution in Rotated Frame for Photon Energies of 375pm25MeV", 2, -180, 180);
  PhiSc425 = new GH1( "Phi_Scattered_425MeV", "Scattered Proton Phi Distribution in Rotated Frame for Photon Energies of 425pm25MeV", 2, -180, 180);
  PhiSc475 = new GH1( "Phi_Scattered_475MeV", "Scattered Proton Phi Distribution in Rotated Frame for Photon Energies of 475pm25MeV", 2, -180, 180);
  PhiSc525 = new GH1( "Phi_Scattered_525MeV", "Scattered Proton Phi Distribution in Rotated Frame for Photon Energies of 525pm25MeV", 2, -180, 180);
  PhiSc575 = new GH1( "Phi_Scattered_575MeV", "Scattered Proton Phi Distribution in Rotated Frame for Photon Energies of 575pm25MeV", 2, -180, 180);
  PhiSc625 = new GH1( "Phi_Scattered_625MeV", "Scattered Proton Phi Distribution in Rotated Frame for Photon Energies of 625pm25MeV", 2, -180, 180);
  PhiSc675 = new GH1( "Phi_Scattered_675MeV", "Scattered Proton Phi Distribution in Rotated Frame for Photon Energies of 675pm25MeV", 2, -180, 180);
  PhiSc725 = new GH1( "Phi_Scattered_725MeV", "Scattered Proton Phi Distribution in Rotated Frame for Photon Energies of 725pm25MeV", 2, -180, 180);
  PhiSc775 = new GH1( "Phi_Scattered_775MeV", "Scattered Proton Phi Distribution in Rotated Frame for Photon Energies of 775pm25MeV", 2, -180, 180);
  PhiSc825 = new GH1( "Phi_Scattered_825MeV", "Scattered Proton Phi Distribution in Rotated Frame for Photon Energies of 825pm25MeV", 2, -180, 180);
  PhiSc875 = new GH1( "Phi_Scattered_875MeV", "Scattered Proton Phi Distribution in Rotated Frame for Photon Energies of 875pm25MeV", 2, -180, 180);

  PhiScNegHel = new GH1("PhiScNegHel", "Scattered Proton Phi Distribution in Rotated Frame for -ve Helicity", 2, -180, 180);
  PhiScPosHel = new GH1("PhiScPosHel", "Scattered Proton Phi Distribution in Rotated Frame for +ve Helicity", 2, -180, 180);

  // Angles of neutron in scattered frame across EGamma bins for negative helicity
  PhiSc275NegHel = new GH1( "Phi_Scattered_275MeV_NegHel", "Scattered Proton Phi Distribution in Rotated Frame for Photon Energies of 275pm25MeV for -ve Helicity", 2, -180, 180);
  PhiSc325NegHel = new GH1( "Phi_Scattered_325MeV_NegHel", "Scattered Proton Phi Distribution in Rotated Frame for Photon Energies of 325pm25MeV for -ve Helicity", 2, -180, 180);
  PhiSc375NegHel = new GH1( "Phi_Scattered_375MeV_NegHel", "Scattered Proton Phi Distribution in Rotated Frame for Photon Energies of 375pm25MeV for -ve Helicity", 2, -180, 180);
  PhiSc425NegHel = new GH1( "Phi_Scattered_425MeV_NegHel", "Scattered Proton Phi Distribution in Rotated Frame for Photon Energies of 425pm25MeV for -ve Helicity", 2, -180, 180);
  PhiSc475NegHel = new GH1( "Phi_Scattered_475MeV_NegHel", "Scattered Proton Phi Distribution in Rotated Frame for Photon Energies of 475pm25MeV for -ve Helicity", 2, -180, 180);
  PhiSc525NegHel = new GH1( "Phi_Scattered_525MeV_NegHel", "Scattered Proton Phi Distribution in Rotated Frame for Photon Energies of 525pm25MeV for -ve Helicity", 2, -180, 180);
  PhiSc575NegHel = new GH1( "Phi_Scattered_575MeV_NegHel", "Scattered Proton Phi Distribution in Rotated Frame for Photon Energies of 575pm25MeV for -ve Helicity", 2, -180, 180);
  PhiSc625NegHel = new GH1( "Phi_Scattered_625MeV_NegHel", "Scattered Proton Phi Distribution in Rotated Frame for Photon Energies of 625pm25MeV for -ve Helicity", 2, -180, 180);
  PhiSc675NegHel = new GH1( "Phi_Scattered_675MeV_NegHel", "Scattered Proton Phi Distribution in Rotated Frame for Photon Energies of 675pm25MeV for -ve Helicity", 2, -180, 180);
  PhiSc725NegHel = new GH1( "Phi_Scattered_725MeV_NegHel", "Scattered Proton Phi Distribution in Rotated Frame for Photon Energies of 725pm25MeV for -ve Helicity", 2, -180, 180);
  PhiSc775NegHel = new GH1( "Phi_Scattered_775MeV_NegHel", "Scattered Proton Phi Distribution in Rotated Frame for Photon Energies of 775pm25MeV for -ve Helicity", 2, -180, 180);
  PhiSc825NegHel = new GH1( "Phi_Scattered_825MeV_NegHel", "Scattered Proton Phi Distribution in Rotated Frame for Photon Energies of 825pm25MeV for -ve Helicity", 2, -180, 180);
  PhiSc875NegHel = new GH1( "Phi_Scattered_875MeV_NegHel", "Scattered Proton Phi Distribution in Rotated Frame for Photon Energies of 875pm25MeV for -ve Helicity", 2, -180, 180);

  // Angles of neutron in scattered frame across EGamma bins for positive helicity
  PhiSc275PosHel = new GH1( "Phi_Scattered_275MeV_PosHel", "Scattered Proton Phi Distribution in Rotated Frame for Photon Energies of 275pm25MeV for +ve Helicity", 2, -180, 180);
  PhiSc325PosHel = new GH1( "Phi_Scattered_325MeV_PosHel", "Scattered Proton Phi Distribution in Rotated Frame for Photon Energies of 325pm25MeV for +ve Helicity", 2, -180, 180);
  PhiSc375PosHel = new GH1( "Phi_Scattered_375MeV_PosHel", "Scattered Proton Phi Distribution in Rotated Frame for Photon Energies of 375pm25MeV for +ve Helicity", 2, -180, 180);
  PhiSc425PosHel = new GH1( "Phi_Scattered_425MeV_PosHel", "Scattered Proton Phi Distribution in Rotated Frame for Photon Energies of 425pm25MeV for +ve Helicity", 2, -180, 180);
  PhiSc475PosHel = new GH1( "Phi_Scattered_475MeV_PosHel", "Scattered Proton Phi Distribution in Rotated Frame for Photon Energies of 475pm25MeV for +ve Helicity", 2, -180, 180);
  PhiSc525PosHel = new GH1( "Phi_Scattered_525MeV_PosHel", "Scattered Proton Phi Distribution in Rotated Frame for Photon Energies of 525pm25MeV for +ve Helicity", 2, -180, 180);
  PhiSc575PosHel = new GH1( "Phi_Scattered_575MeV_PosHel", "Scattered Proton Phi Distribution in Rotated Frame for Photon Energies of 575pm25MeV for +ve Helicity", 2, -180, 180);
  PhiSc625PosHel = new GH1( "Phi_Scattered_625MeV_PosHel", "Scattered Proton Phi Distribution in Rotated Frame for Photon Energies of 625pm25MeV for +ve Helicity", 2, -180, 180);
  PhiSc675PosHel = new GH1( "Phi_Scattered_675MeV_PosHel", "Scattered Proton Phi Distribution in Rotated Frame for Photon Energies of 675pm25MeV for +ve Helicity", 2, -180, 180);
  PhiSc725PosHel = new GH1( "Phi_Scattered_725MeV_PosHel", "Scattered Proton Phi Distribution in Rotated Frame for Photon Energies of 725pm25MeV for +ve Helicity", 2, -180, 180);
  PhiSc775PosHel = new GH1( "Phi_Scattered_775MeV_PosHel", "Scattered Proton Phi Distribution in Rotated Frame for Photon Energies of 775pm25MeV for +ve Helicity", 2, -180, 180);
  PhiSc825PosHel = new GH1( "Phi_Scattered_825MeV_PosHel", "Scattered Proton Phi Distribution in Rotated Frame for Photon Energies of 825pm25MeV for +ve Helicity", 2, -180, 180);
  PhiSc875PosHel = new GH1( "Phi_Scattered_875MeV_PosHel", "Scattered Proton Phi Distribution in Rotated Frame for Photon Energies of 875pm25MeV for +ve Helicity", 2, -180, 180);

  ThetaRecPiDiff = new GH1 ("ThetaRecPiDiff", "Difference between ThetaPiRec and Thetan", 200, 0, 180);
  ThetanThetaRecPi = new GH2 ("ThetanThetaRecPi", "Thetan vs ThetaPiRec", 100, 0, 180, 100, 0, 180);
  ThetanThetaRecPiDiff = new GH2 ("ThetanThetaRecPiDiff", "Thetan vs (ThetaPiRec - Thetan)", 100, 0, 180, 100, 0, 180);
  ThetaRecPDiff = new GH1 ("ThetaRecPDiff", "Difference between ThetaPRec and Thetan", 200, 0, 180);
  ThetanThetaRecP = new GH2 ("ThetanThetaRecP", "Thetan vs ThetaPRec", 100, 0, 180, 100, 0, 180);
  ThetanThetaRecPDiff = new GH2 ("ThetanThetaRecPDiff", "Thetan vs (ThetaPRec - Thetan)", 100, 0, 180, 100, 0, 180);

  E_dE = new GH2 ("E_dE", "EdE Plot With E Loss Adjustment", 100, 0, 500, 100, 0, 5);
  KinEp_dE = new GH2 ("KinEp_dE", "KinEpdE Plot", 100, 0, 500, 100, 0, 5);
  //KinEp_dE_GoodCut = new GH2 ("KinEp_dE_GoodCut", "KinEpdE Plot With Good Proton Cut", 100, 0, 500, 100, 0, 5);
  ThetaScPhiSc = new GH2 ("ThetaScPhiSc", "Phi as a function of Theta (Both in rotated frame)", 100, 0, 180, 100, -180, 180);
  E_KinEp = new GH2 ("E_KinEp", "Kinematic Energy of Proton as a function of CB energy", 100, 0, 500, 100, 0, 500);
  PhinDiffWCZRec = new GH2 ("PhinDiffWCZRec", "Difference between WC Phi and Reconstructed Phi as a fn of WCZ Hit Position", 100, 0, 200, 100, 0, 180);
  ThetaDiffPhiDiff = new GH2 ("ThetaDiffPhiDiff", "PhiDiff as a Fn of ThetaDiff (Detected-Rec)", 100, 0, 180, 100, 0, 180);

}

void PNeutPol_Polarimeter_Circ::FillHists()
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
    ThetaDiffPhiDiff->Fill((abs(Thetan-ThetanRec)), (abs(Phin-PhinRec)), TaggerTime );

    ThetaSc -> Fill(ScattTheta, TaggerTime);
    PhiSc -> Fill(ScattPhi, TaggerTime);
    ThetaScPhiSc->Fill(ScattTheta, ScattPhi, TaggerTime);

    if(ScattPhi < -165){
        ZpPhiScatNeg180->Fill(Zp, TaggerTime);
    }

    if(ScattPhi < 15 && ScattPhi > -15){
        ZpPhiScat0->Fill(Zp, TaggerTime);
    }

    if(ScattPhi > 165){
        ZpPhiScatPos180->Fill(Zp, TaggerTime);
    }

    if (BeamHelicity == kFALSE) PhiScNegHel->Fill(ScattPhi, TaggerTime);
    if (BeamHelicity == kTRUE) PhiScPosHel->Fill(ScattPhi, TaggerTime);

    if(200 < EGamma && EGamma < 300){
        MMp200300->Fill(MMpEpCorr, TaggerTime);
    }

    else if(300 < EGamma && EGamma < 400){
        MMp300400->Fill(MMpEpCorr, TaggerTime);
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

    if ( 250 < EGamma && EGamma < 300) {
        PhiSc275->Fill(ScattPhi, TaggerTime);
        if (BeamHelicity == kFALSE) PhiSc275NegHel->Fill(ScattPhi, TaggerTime);
        else if (BeamHelicity == kTRUE) PhiSc275PosHel->Fill(ScattPhi, TaggerTime);
    }

    else if ( 300 < EGamma && EGamma < 350) {
        PhiSc325->Fill(ScattPhi, TaggerTime);
        if (BeamHelicity == kFALSE) PhiSc325NegHel->Fill(ScattPhi, TaggerTime);
        else if (BeamHelicity == kTRUE) PhiSc325PosHel->Fill(ScattPhi, TaggerTime);
    }

    else if ( 350 < EGamma && EGamma < 400) {
        PhiSc375->Fill(ScattPhi, TaggerTime);
        if (BeamHelicity == kFALSE) PhiSc375NegHel->Fill(ScattPhi, TaggerTime);
        else if (BeamHelicity == kTRUE) PhiSc375PosHel->Fill(ScattPhi, TaggerTime);
    }

    else if ( 400 < EGamma && EGamma < 450) {
        PhiSc425->Fill(ScattPhi, TaggerTime);
        if (BeamHelicity == kFALSE) PhiSc425NegHel->Fill(ScattPhi, TaggerTime);
        else if (BeamHelicity == kTRUE) PhiSc425PosHel->Fill(ScattPhi, TaggerTime);
    }

    else if ( 450 < EGamma && EGamma < 500) {
        PhiSc475->Fill(ScattPhi, TaggerTime);
        if (BeamHelicity == kFALSE) PhiSc475NegHel->Fill(ScattPhi, TaggerTime);
        else if (BeamHelicity == kTRUE) PhiSc475PosHel->Fill(ScattPhi, TaggerTime);
    }

    else if ( 500 < EGamma && EGamma < 550) {
        PhiSc525->Fill(ScattPhi, TaggerTime);
        if (BeamHelicity == kFALSE) PhiSc525NegHel->Fill(ScattPhi, TaggerTime);
        else if (BeamHelicity == kTRUE) PhiSc525PosHel->Fill(ScattPhi, TaggerTime);
    }

    else if ( 550 < EGamma && EGamma < 600) {
        PhiSc575->Fill(ScattPhi, TaggerTime);
        if (BeamHelicity == kFALSE) PhiSc575NegHel->Fill(ScattPhi, TaggerTime);
        else if (BeamHelicity == kTRUE) PhiSc575PosHel->Fill(ScattPhi, TaggerTime);
    }

    else if ( 600 < EGamma && EGamma < 650) {
        PhiSc625->Fill(ScattPhi, TaggerTime);
        if (BeamHelicity == kFALSE) PhiSc625NegHel->Fill(ScattPhi, TaggerTime);
        else if (BeamHelicity == kTRUE) PhiSc625PosHel->Fill(ScattPhi, TaggerTime);
    }

    else if ( 650 < EGamma && EGamma < 700) {
        PhiSc675->Fill(ScattPhi, TaggerTime);
        if (BeamHelicity == kFALSE) PhiSc675NegHel->Fill(ScattPhi, TaggerTime);
        else if (BeamHelicity == kTRUE) PhiSc675PosHel->Fill(ScattPhi, TaggerTime);
    }

    else if ( 700 < EGamma && EGamma < 750) {
        PhiSc725->Fill(ScattPhi, TaggerTime);
        if (BeamHelicity == kFALSE) PhiSc725NegHel->Fill(ScattPhi, TaggerTime);
        else if (BeamHelicity == kTRUE) PhiSc725PosHel->Fill(ScattPhi, TaggerTime);
    }

    else if ( 750 < EGamma && EGamma < 800) {
        PhiSc775->Fill(ScattPhi, TaggerTime);
        if (BeamHelicity == kFALSE) PhiSc775NegHel->Fill(ScattPhi, TaggerTime);
        else if (BeamHelicity == kTRUE) PhiSc775PosHel->Fill(ScattPhi, TaggerTime);
    }

    else if ( 800 < EGamma && EGamma < 850) {
        PhiSc825->Fill(ScattPhi, TaggerTime);
        if (BeamHelicity == kFALSE) PhiSc825NegHel->Fill(ScattPhi, TaggerTime);
        else if (BeamHelicity == kTRUE) PhiSc825PosHel->Fill(ScattPhi, TaggerTime);
    }

    else if ( 850 < EGamma && EGamma < 900) {
        PhiSc875->Fill(ScattPhi, TaggerTime);
        if (BeamHelicity == kFALSE) PhiSc875NegHel->Fill(ScattPhi, TaggerTime);
        else if (BeamHelicity == kTRUE) PhiSc875PosHel->Fill(ScattPhi, TaggerTime);
    }
}

Bool_t	PNeutPol_Polarimeter_Circ::Write(){
  // Write all GH1's easily

  GTreeManager::Write();
}
