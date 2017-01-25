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
  Deut = TLorentzVector (0., 0., 0., 1875.613); // 4-Vector of Deuterium target, assume at rest

  Cut_CB_proton = OpenCutFile("configfiles/cuts/CB_DeltaE-E_Proton_7_12_16.root", "Proton"); // These will need adjusting with new Acqu files
  Cut_proton = Cut_CB_proton;
  Cut_CB_pion = OpenCutFile("configfiles/cuts/CB_DeltaE-E_Pion_29_07_15.root", "Pion");
  Cut_pion = Cut_CB_pion;
  Cut_CB_protonKinGood = OpenCutFile("configfiles/cuts/CB_DeltaE-E_ProtonKinGood_20_1_17.root", "ProtonKinGood"); // These will need adjusting with new Acqu files
  Cut_protonKinGood = Cut_CB_protonKinGood;
  Cut_CB_protonKinBad = OpenCutFile("configfiles/cuts/CB_DeltaE-E_ProtonKinBad_15_12_16.root", "ProtonKinBad");
  Cut_protonKinBad = Cut_CB_protonKinBad;
  cout << endl;

  TraverseValidEvents(); // This loops over each event as in old file and calls ProcessEvent() each loop

  cout << EventCounter << " Events in file " << EventCounterTrackCut << " Events After Track Cut " << EventCounterZCut << " Events after Z cut " << EventCounterCoplanarCut << " Events after Coplanarity Cut" << endl;

  return kTRUE;
}

void	PNeutPol_Polarimeter::ProcessEvent()
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

  zdiff = abs(Zp - Zn);
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

  if( Zp > 40 || Zp < -50) return;

  EventCounterZCut++;

  if ( PhiWCDiff > 195 || PhiWCDiff < 165) return;

  EventCounterCoplanarCut++;

  for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++){

    TaggerTime = GetTagger()->GetTaggedTime(j); // Get tagged time for event
    EGamma = (GetTagger()->GetTaggedEnergy(j)); // Get Photon energy for event
    Gamma = TLorentzVector (0., 0., EGamma , EGamma); // 4-Vector of Photon beam
    BeamHelicity = GetTrigger()->GetHelicity();
    if (BeamHelicity == kFALSE) cout << "False" << endl;
    else if (BeamHelicity == kTRUE) cout<< "True" <<endl;

    Thetap = GVp3.Theta()*TMath::RadToDeg(); // Lab frame angles for proton/neutron
    Phip = GVp3.Phi()*TMath::RadToDeg();
    Thetan = GVn3.Theta()*TMath::RadToDeg();
    Phin = GVn3.Phi()*TMath::RadToDeg();

    KinEp = CalcKinEnergy(WCThetap, EGamma, Md, 0., Mp, Mn);
    EpCorr = EpPolCorrect(Ep, WCThetap);
    KinEpMB = CalcKinEnergyMB(WCThetap, EGamma, Md, 0., Mp, Mn);
    EpDiff = abs(EpCorr - Ep);
    KinEDiff = KinEp - EpCorr;

    RecKinProton = Proton4VectorKin(KinEp, WCThetapRad, WCPhipRad);
    RecKinNeutron = Neutron4VectorKin(RecKinProton);
    ThetanRec = (RecKinNeutron.Theta()) * TMath::RadToDeg();
    MMpKin = RecKinNeutron.M();

    RecProtonEpCorr = Proton4VectorKin(EpCorr, WCThetapRad, WCPhipRad);
    RecNeutronEpCorr = Neutron4VectorKin(RecProtonEpCorr);
    MMpEpCorr = RecNeutronEpCorr.M();
    RecProtonEpCorr3 = RecProtonEpCorr.Vect();
    RecNeutronEpCorr3 = RecNeutronEpCorr.Vect();

    P3Vect = RecKinProton.Vect();
    N3Vect = RecKinNeutron.Vect();
    OpeningAngle = (N3Vect.Angle(GVn3))*TMath::RadToDeg();

    ThetanDiff = abs(ThetanRec - WCThetan);

    TVector3 ScattAngles = ScatteredFrameAngles(RecNeutronEpCorr3, GVp3, GVn3, Gamma);
    ScattTheta = ScattAngles(0);
    ScattPhi = ScattAngles(1);

    if( 850 > MMpEpCorr || 1050 < MMpEpCorr) continue;

    //if (ScattTheta > 90) continue;

    //if (abs(KinEDiff) > 100) continue; // If difference between CB energy and calculated Energy for proton > 100MeV continue

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

TLorentzVector PNeutPol_Polarimeter::Proton4VectorKin(Double_t KinE, Double_t Theta, Double_t Phi)
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

TLorentzVector PNeutPol_Polarimeter::Neutron4VectorKin(TLorentzVector ProtonKinVector)
{

    N4Vect = (Gamma + Deut) - ProtonKinVector;

    return N4Vect;
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
  EpCorrected = new GH1 ("EpCorrected", "Ep Corrected for Energy Loss in Polarimeter ", 100, 0, 500);
  OAngle = new GH1 ("OAngle", "Opening Angle between P and N Vectors", 180, 0, 180);

  ThetaSc =  new GH1( "Theta_Scattered", "Scattetred Proton Theta Distribution in Rotated Frame", 180, 0, 180 );
  PhiSc = new GH1( "Phi_Scattered", "Scattetred Proton Phi Distribution in Rotated Frame", 90, -180, 180 );

  EpKinEpCorrDiff = new GH1("EpKinEpCorrDiff", "Difference Between EpKin and EpCorr", 300, -300, 300);
  EpEpCorrDiff = new GH1("EpEpCorrDiff", "Difference Between Ep and EpCorr", 200, 0, 200);

  WCXp = new GH1("WCXp", "WC X Position for Proton", 200, -100, 100);
  WCYp = new GH1("WCYp", "WC Y Position for Proton", 200, -100, 100);
  WCZp = new GH1("WCZp", "WC Z Position for Proton", 200, -500, 500);
  WCXn = new GH1("WCXn", "WC X Position for Neutron", 200, -100, 100);
  WCYn = new GH1("WCYn", "WC Y Position for Neutron", 200, -100, 100);
  WCZn = new GH1("WCZn", "WC Z Position for Neutron", 200, -500, 500);
  MMp = new GH1 ("MMp", "Missing mass seen by Proton", 400, 800, 1000);
  MMpEpCorrected = new GH1 ("MMpEpCorrected", "Missing mass seen by Proton (E Loss Corrected)", 400, 0, 2000);

  MMpEpCorrectedCut =  new GH1 ("MMpEpCorrectedCut", "Missing mass seen by Proton (E Loss Corrected, P Banana Cut)", 400, 0, 2000);
  ThetanWCThetanRecDiff = new GH1 ("ThetanWCThetanRecDiff", "Difference between ThetaWC and ThetaRec for n", 180, 0, 180);
  OAngleCut = new GH1 ("OAngleCut", "Opening Angle between P and N Vectors (P Banana Cut)", 180, 0, 180);
  OAngleCut200400 = new GH1 ("OAngleCut200400", "Opening Angle between P and N Vectors (P Banana Cut, 200-400MeV Gamma)", 180, 0, 180);
  EgCut = new GH1( "EgCut", "Photon Energy Distribution (P Banana Cut)", 400, 100, 1600 );

  // ThetaWCThetaRecDiff across EGamma Bins
  ThetanWCThetanRecDiff200300 = new GH1 ("ThetanWCThetanRecDiff200300", "Difference between ThetaWC and ThetaRec for n (200-300MeV Photon Energy)", 180, 0, 180);
  ThetanWCThetanRecDiff300400 = new GH1 ("ThetanWCThetanRecDiff300400", "Difference between ThetaWC and ThetaRec for n (300-400MeV Photon Energy)", 180, 0, 180);
  ThetanWCThetanRecDiff400500 = new GH1 ("ThetanWCThetanRecDiff400500", "Difference between ThetaWC and ThetaRec for n (400-500MeV Photon Energy)", 180, 0, 180);
  ThetanWCThetanRecDiff500600 = new GH1 ("ThetanWCThetanRecDiff500600", "Difference between ThetaWC and ThetaRec for n (500-600MeV Photon Energy)", 180, 0, 180);
  ThetanWCThetanRecDiff600700 = new GH1 ("ThetanWCThetanRecDiff600700", "Difference between ThetaWC and ThetaRec for n (600-700MeV Photon Energy)", 180, 0, 180);
  ThetanWCThetanRecDiff700800 = new GH1 ("ThetanWCThetanRecDiff700800", "Difference between ThetaWC and ThetaRec for n (700-800MeV Photon Energy)", 180, 0, 180);
  ThetanWCThetanRecDiff800900 = new GH1 ("ThetanWCThetanRecDiff800900", "Difference between ThetaWC and ThetaRec for n (800-900MeV Photon Energy)", 180, 0, 180);

  // MMp across photon E bins
  MMp200300 = new GH1("MMp200300", "Missing mass as seen by Proton (200-300MeV Photon Energy)", 400, 0, 2000);
  MMp300400 = new GH1("MMp300400", "Missing mass as seen by Proton (300-400MeV Photon Energy)", 400, 0, 2000);
  MMp400500 = new GH1("MMp400500", "Missing mass as seen by Proton (400-500MeV Photon Energy)", 400, 0, 2000);
  MMp500600 = new GH1("MMp500600", "Missing mass as seen by Proton (500-600MeV Photon Energy)", 400, 0, 2000);
  MMp600700 = new GH1("MMp600700", "Missing mass as seen by Proton (600-700MeV Photon Energy)", 400, 0, 2000);
  MMp700800 = new GH1("MMp700800", "Missing mass as seen by Proton (700-800MeV Photon Energy)", 400, 0, 2000);
  MMp800900 = new GH1("MMp800900", "Missing mass as seen by Proton (800-900MeV Photon Energy)", 400, 0, 2000);

  // Angles of neutron in scattered frame across EGamma bins
  PhiSc275 = new GH1( "Phi_Scattered_275MeV", "Scattetred Proton Phi Distribution in Rotated Frame for Photon Energies of 275pm25MeV", 90, -180, 180);
  PhiSc325 = new GH1( "Phi_Scattered_325MeV", "Scattetred Proton Phi Distribution in Rotated Frame for Photon Energies of 325pm25MeV", 90, -180, 180);
  PhiSc375 = new GH1( "Phi_Scattered_375MeV", "Scattetred Proton Phi Distribution in Rotated Frame for Photon Energies of 375pm25MeV", 90, -180, 180);
  PhiSc425 = new GH1( "Phi_Scattered_425MeV", "Scattetred Proton Phi Distribution in Rotated Frame for Photon Energies of 425pm25MeV", 90, -180, 180);
  PhiSc475 = new GH1( "Phi_Scattered_475MeV", "Scattetred Proton Phi Distribution in Rotated Frame for Photon Energies of 475pm25MeV", 90, -180, 180);
  PhiSc525 = new GH1( "Phi_Scattered_525MeV", "Scattetred Proton Phi Distribution in Rotated Frame for Photon Energies of 525pm25MeV", 90, -180, 180);
  PhiSc575 = new GH1( "Phi_Scattered_575MeV", "Scattetred Proton Phi Distribution in Rotated Frame for Photon Energies of 575pm25MeV", 90, -180, 180);
  PhiSc625 = new GH1( "Phi_Scattered_625MeV", "Scattetred Proton Phi Distribution in Rotated Frame for Photon Energies of 625pm25MeV", 90, -180, 180);
  PhiSc675 = new GH1( "Phi_Scattered_675MeV", "Scattetred Proton Phi Distribution in Rotated Frame for Photon Energies of 675pm25MeV", 90, -180, 180);
  PhiSc725 = new GH1( "Phi_Scattered_725MeV", "Scattetred Proton Phi Distribution in Rotated Frame for Photon Energies of 725pm25MeV", 90, -180, 180);
  PhiSc775 = new GH1( "Phi_Scattered_775MeV", "Scattetred Proton Phi Distribution in Rotated Frame for Photon Energies of 775pm25MeV", 90, -180, 180);
  PhiSc825 = new GH1( "Phi_Scattered_825MeV", "Scattetred Proton Phi Distribution in Rotated Frame for Photon Energies of 825pm25MeV", 90, -180, 180);
  PhiSc875 = new GH1( "Phi_Scattered_875MeV", "Scattetred Proton Phi Distribution in Rotated Frame for Photon Energies of 875pm25MeV", 90, -180, 180);

  PhiScNegHel = new GH1("PhiScNegHel", "Scattetred Proton Phi Distribution in Rotated Frame for -ve Helicity", 90, -180, 180);
  PhiScPosHel = new GH1("PhiScPosHel", "Scattetred Proton Phi Distribution in Rotated Frame for +ve Helicity", 90, -180, 180);

  E_dE = new GH2 ("E_dE", "EdE Plot With E Loss Adjustment", 100, 0, 500, 100, 0, 5);
  E_dE_Cut = new GH2 ("E_dE_Cut", "EdE Plot (With cut on proton banana + E Loss)", 100, 0, 500, 100, 0, 5);
  KinEp_dE = new GH2 ("KinEp_dE", "KinEpdE Plot", 100, 0, 500, 100, 0, 5);
  KinEp_dE_GoodCut = new GH2 ("KinEp_dE_GoodCut", "KinEpdE Plot With Good Proton Cut", 100, 0, 500, 100, 0, 5);
  ThetaScPhiSc = new GH2 ("ThetaScPhiSc", "Phi as a function of Theta (Both in rotated frame)", 100, 0, 180, 100, -180, 180);

  //MMp as fn of ThetaP/EpKin across Photon E bins
  //MMpThetap200300 = new GH2("MMpThetap200300", "MMp as a fn of WC Thetap (200-300MeV Photon Energy)", 150, 0, 2000, 150, 0, 180);
  //MMpThetap300400 = new GH2("MMpThetap300400", "MMp as a fn of WC Thetap (300-400MeV Photon Energy)", 150, 0, 2000, 150, 0, 180);
  //MMpThetap400500 = new GH2("MMpThetap400500", "MMp as a fn of WC Thetap (400-500MeV Photon Energy)", 150, 0, 2000, 150, 0, 180);
  //MMpThetap500600 = new GH2("MMpThetap500600", "MMp as a fn of WC Thetap (500-600MeV Photon Energy)", 150, 0, 2000, 150, 0, 180);
  //MMpThetap600700 = new GH2("MMpThetap600700", "MMp as a fn of WC Thetap (600-700MeV Photon Energy)", 150, 0, 2000, 150, 0, 180);
  //MMpThetap700800 = new GH2("MMpThetap700800", "MMp as a fn of WC Thetap (700-800MeV Photon Energy)", 150, 0, 2000, 150, 0, 180);
  //MMpThetap800900 = new GH2("MMpThetap800900", "MMp as a fn of WC Thetap (800-900MeV Photon Energy)", 150, 0, 2000, 150, 0, 180);
  //MMpEpKin200300 = new GH2("MMpEpKin200300", "MMp as a fn of EpKin (200-300MeV Photon Energy)", 150, 0, 2000, 150, 80, 500);
  //MMpEpKin300400 = new GH2("MMpEpKin300400", "MMp as a fn of EpKin (300-400MeV Photon Energy)", 150, 0, 2000, 150, 80, 500);
  //MMpEpKin400500 = new GH2("MMpEpKin400500", "MMp as a fn of EpKin (400-500MeV Photon Energy)", 150, 0, 2000, 150, 80, 500);
  //MMpEpKin500600 = new GH2("MMpEpKin500600", "MMp as a fn of EpKin (500-600MeV Photon Energy)", 150, 0, 2000, 150, 80, 500);
  //MMpEpKin600700 = new GH2("MMpEpKin600700", "MMp as a fn of EpKin (600-700MeV Photon Energy)", 150, 0, 2000, 150, 80, 500);
  //MMpEpKin700800 = new GH2("MMpEpKin700800", "MMp as a fn of EpKin (700-800MeV Photon Energy)", 150, 0, 2000, 150, 80, 500);
  //MMpEpKin800900 = new GH2("MMpEpKin800900", "MMp as a fn of EpKin (800-900MeV Photon Energy)", 150, 0, 2000, 150, 80, 500);
}

void PNeutPol_Polarimeter::FillHists()
{
  time->Fill(TaggerTime);
  if (-5 < TaggerTime && TaggerTime < 20) time_cut->Fill(TaggerTime);

  Zp_Vert->Fill(Zp, TaggerTime);
  Ekp->Fill(Ep, TaggerTime);
  Ekn->Fill(En, TaggerTime);
  EkSum->Fill((Ep + En),TaggerTime);
  Eg->Fill(EGamma, TaggerTime);
  ThetaProt->Fill(Thetap, TaggerTime);
  ThetaNeut->Fill(Thetan, TaggerTime);
  PhiProt->Fill(Phip, TaggerTime);
  PhiNeut->Fill(Phin, TaggerTime);
  WCThetaProt->Fill(WCThetap, TaggerTime);
  WCPhiProt->Fill(WCPhip, TaggerTime);
  E_dE->Fill(EpCorr, dEp, TaggerTime);
  KinEp_dE->Fill(KinEp, dEp, TaggerTime);
  EpKin->Fill(KinEp, TaggerTime);
  EpCorrected->Fill(EpCorr, TaggerTime);
  EpKinEpCorrDiff->Fill(KinEDiff, TaggerTime);
  EpEpCorrDiff->Fill(EpDiff, TaggerTime);
  WCXp->Fill(WC1pX, TaggerTime);
  WCYp->Fill(WC1pY, TaggerTime);
  WCZp->Fill(WC1pZ, TaggerTime);
  MMp->Fill(MMpKin, TaggerTime);
  MMpEpCorrected->Fill(MMpEpCorr, TaggerTime);
  OAngle->Fill(OpeningAngle, TaggerTime);
  Zn_Vert->Fill(Zn, TaggerTime);
  WCThetaNeut->Fill(WCThetan, TaggerTime);
  WCPhiNeut->Fill(WCPhin, TaggerTime);
  WCXn->Fill(WC1nX, TaggerTime);
  WCYn->Fill(WC1nY, TaggerTime);
  WCZn->Fill(WC1nZ, TaggerTime);
  WCPhiDifference->Fill(PhiWCDiff);

  if(Cut_proton -> IsInside(EpCorr, dEp) == kTRUE){
    E_dE_Cut->Fill(EpCorr, dEp, TaggerTime);
  }

  // Fill events inside good proton banana on KinEpdE plot
  if(Cut_protonKinGood -> IsInside(KinEp, dEp) == kTRUE)
  {
    KinEp_dE_GoodCut->Fill(KinEp, dEp, TaggerTime);
    MMpEpCorrectedCut->Fill(MMpEpCorr, TaggerTime);
    EgCut->Fill(EGamma, TaggerTime);
    OAngleCut->Fill(OpeningAngle, TaggerTime);
    ThetanWCThetanRecDiff->Fill(ThetanDiff, TaggerTime);
    ThetaSc -> Fill(ScattTheta, TaggerTime);
    PhiSc -> Fill(ScattPhi, TaggerTime);
    ThetaScPhiSc->Fill(ScattTheta, ScattPhi, TaggerTime);

    //if (BeamHelicity == kFALSE) PhiScNegHel->Fill(ScattPhi, TaggerTime);
    //if (BeamHelicity == kTRUE) PhiScPosHel->Fill(ScattPhi, TaggerTime);

    if(200 < EGamma && EGamma < 300){
        MMp200300->Fill(MMpEpCorr, TaggerTime);
        //MMpThetap200300->Fill(MMpEpCorr, WCThetap, TaggerTime);
        //MMpEpKin200300->Fill(MMpEpCorr, KinEp, TaggerTime);
        ThetanWCThetanRecDiff200300->Fill(ThetanDiff, TaggerTime);
        OAngleCut200400->Fill(OpeningAngle, TaggerTime);
    }

    else if(300 < EGamma && EGamma < 400){
        MMp300400->Fill(MMpEpCorr, TaggerTime);
        //MMpThetap300400->Fill(MMpEpCorr, WCThetap, TaggerTime);
        //MMpEpKin300400->Fill(MMpEpCorr, KinEp, TaggerTime);
        ThetanWCThetanRecDiff300400->Fill(ThetanDiff, TaggerTime);
        OAngleCut200400->Fill(OpeningAngle, TaggerTime);
    }

    else if(400 < EGamma && EGamma < 500){
        MMp400500->Fill(MMpEpCorr, TaggerTime);
        //MMpThetap400500->Fill(MMpEpCorr, WCThetap, TaggerTime);
        //MMpEpKin400500->Fill(MMpEpCorr, KinEp, TaggerTime);
        ThetanWCThetanRecDiff400500->Fill(ThetanDiff, TaggerTime);
    }

    else if(500 < EGamma && EGamma < 600){
        MMp500600->Fill(MMpEpCorr, TaggerTime);
        //MMpThetap500600->Fill(MMpEpCorr, WCThetap, TaggerTime);
        //MMpEpKin500600->Fill(MMpEpCorr, KinEp, TaggerTime);
        ThetanWCThetanRecDiff500600->Fill(ThetanDiff, TaggerTime);
    }

    else if(600 < EGamma && EGamma < 700){
        MMp600700->Fill(MMpEpCorr, TaggerTime);
        //MMpThetap600700->Fill(MMpEpCorr, WCThetap, TaggerTime);
        //MMpEpKin600700->Fill(MMpEpCorr, KinEp, TaggerTime);
        ThetanWCThetanRecDiff600700->Fill(ThetanDiff, TaggerTime);
    }

    else if(700 < EGamma && EGamma < 800){
        MMp700800->Fill(MMpEpCorr, TaggerTime);
        //MMpThetap700800->Fill(MMpEpCorr, WCThetap, TaggerTime);
        //MMpEpKin700800->Fill(MMpEpCorr, KinEp, TaggerTime);
        ThetanWCThetanRecDiff700800->Fill(ThetanDiff, TaggerTime);
    }

    else if(800 < EGamma && EGamma < 900){
        MMp800900->Fill(MMpEpCorr, TaggerTime);
        //MMpThetap800900->Fill(MMpEpCorr, WCThetap, TaggerTime);
        //MMpEpKin800900->Fill(MMpEpCorr, KinEp, TaggerTime);
        ThetanWCThetanRecDiff800900->Fill(ThetanDiff, TaggerTime);
    }

    if ( 250 < EGamma && EGamma < 300) PhiSc275->Fill(ScattPhi, TaggerTime);
    if ( 300 < EGamma && EGamma < 350) PhiSc325->Fill(ScattPhi, TaggerTime);
    if ( 350 < EGamma && EGamma < 400) PhiSc375->Fill(ScattPhi, TaggerTime);
    if ( 400 < EGamma && EGamma < 450) PhiSc425->Fill(ScattPhi, TaggerTime);
    if ( 450 < EGamma && EGamma < 500) PhiSc475->Fill(ScattPhi, TaggerTime);
    if ( 500 < EGamma && EGamma < 550) PhiSc525->Fill(ScattPhi, TaggerTime);
    if ( 550 < EGamma && EGamma < 600) PhiSc575->Fill(ScattPhi, TaggerTime);
    if ( 600 < EGamma && EGamma < 650) PhiSc625->Fill(ScattPhi, TaggerTime);
    if ( 650 < EGamma && EGamma < 700) PhiSc675->Fill(ScattPhi, TaggerTime);
    if ( 700 < EGamma && EGamma < 750) PhiSc725->Fill(ScattPhi, TaggerTime);
    if ( 750 < EGamma && EGamma < 800) PhiSc775->Fill(ScattPhi, TaggerTime);
    if ( 800 < EGamma && EGamma < 850) PhiSc825->Fill(ScattPhi, TaggerTime);
    if ( 850 < EGamma && EGamma < 900) PhiSc875->Fill(ScattPhi, TaggerTime);

  }

}

Bool_t	PNeutPol_Polarimeter::Write(){
  // Write all GH1's easily

  GTreeManager::Write();
}
