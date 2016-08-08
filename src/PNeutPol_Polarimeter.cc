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
  InitialVect(); // Function gets vectors of identified tracks and returns them
  GV1_3 = GV1.Vect();
  GV2_3 = GV2.Vect();
  InitialProp(); // Function gets initial properties (energy, vertex e.t.c.) of identified tracks
  DetectorCheck(); // Function checks detector numbers for each track
  DetectorsSum = Detectors1 + Detectors2; // To make splitting into separate histograms less of a pain in the arse

  //cout << DetectorsSum << endl;

  // Currently try to accept only MWPC+PID+NaI and NaI + MWPC OR Only NaI events
  // If track 1 only gives signals in MWPC it is the neutron
  if((Detectors1 == 7) && (Detectors2 == 1))
  {
    Proton1 = kTRUE;
    Proton2 = kFALSE;
    //cout << Detectors1 << "   " << Detectors2 << endl;
  }

  // If track 2 only gives signals in MWPC it is the neutron
  else if((Detectors1 == 1) && (Detectors2 == 7))
  {
    Proton1 = kFALSE;
    Proton2 = kTRUE;
  }

  else if((Detectors1 == 7) && (Detectors2 == 5))
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

    //if (dE1 < 0.5 || dE2 < 0.5) continue; // Cut out low PID energy events // This won't work if PID has no hits as we desire!

    // If neither particle is in the proton banana region continue
    //if (Cut_proton -> IsInside(E1, dE1) == kFALSE && Cut_proton -> IsInside(E2, dE2) == kFALSE) continue; // if neither particle is in proton region drop out
    //if (Cut_proton -> IsInside(E1, dE1) == kFALSE && Cut_pion -> IsInside(E1, dE1) == kFALSE) continue; // If not in proton or pion region drop out
    //if (Cut_proton -> IsInside(E2, dE2) == kFALSE && Cut_pion -> IsInside(E2, dE2) == kFALSE) continue;

    //if (i == 0) { // These properties get defined for each photon

    Gamma = TLorentzVector (0., 0., EGamma , EGamma); // 4-Vector of Photon beam
    Gamma3 = Gamma.Vect(); // Convert photon beam 4-vector to 3-vector
    B = (Deut + Gamma).Beta(); // Calculte Beta
    b = TVector3(0., 0., B); // Define boost vector
    ReconstructVectors();
    GV1Rec_3 = GV1Rec.Vect();
    GV2Rec_3 = GV2Rec.Vect();
    mm1 = GV1Rec.M(); // Calculate missing mass of each particle
    mm2 = GV2Rec.M();
    mm1Diff = abs(Mn-mm1); // Look at difference of missing mass from neutron mass
    mm2Diff = abs(Mn-mm2);

    FillTime(*GetProtons(),time);
    FillTimeCut(*GetProtons(),time_cut);

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

    GVp3 = GVp.Vect(); // Generate some 3-vectors from the 4-vectors we have
    GVn3 = GVn.Vect();
    GVpCalc3 = GVpCalc.Vect();
    GVnCalc3 = GVnCalc.Vect();

    if ((mmn < 850) || (mmn > 1050)) continue; //If missing mass for particle that we think is the neutron is not correct, continue
    if (Cut_proton -> IsInside(Ep, dEp) == kFALSE) continue; // If proton not in banana drop out

    LabBoost(); // Boost particles in lab frame and return results
    LabScatter(); // Work out scattering angle in lab frame and return results
    // Drop this for now need to think about it
    //WCVertex(GVn3, GVnCalc3, Zp, Zn); // Calculate Z vertex location in WC from measured and reconstructed vectors

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

  return Theta1, Theta2, Phi1, Phi2, z1, z2, E1, E2, dE1, dE2; // Returns various quantities used in later functions
}

Int_t PNeutPol_Polarimeter::DetectorCheck()
{
    Detectors1 = GetTracks()->GetDetectors(0); //Gets number for detectors that registered hits
    Detectors2 = GetTracks()->GetDetectors(1); // 7 = NaI + PID + MWPC, 5 = NaI + MWPC
    return Detectors1, Detectors2;
}

TLorentzVector PNeutPol_Polarimeter::ReconstructVectors()
{

    GV1Rec = (Gamma + Deut) - GV2; //Assume GV2 correct, reconstruct 1st vector
    GV2Rec = (Gamma + Deut) - GV1;

    return GV1Rec, GV2Rec;

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
      mmp = mm2;
      mmn = mm1;
      Ep = E1;
      En = E2;
      dEp = dE1;
      dEn = dE2;
    }

  if(ProtonParticleNumber == 2)
    {
      Zp = z2; // First particle is neutron, second is proton
      Zn = z1;
      zdiff = abs (Zp - Zn);
      mmp = mm1; // Note this looks incorrect but MM1 is the missing mass of particle 1 as calculated using particle 2, we have determined that particle 2 is the proton here
      mmn = mm2;
      Ep = E2; // Therefore the quantity mmp is the amount of missing mass we see when we do a kinematics calculation USING the proton
      En = E1;
      dEp = dE2;
      dEn = dE1;

    }

  return Zp, Zn, zdiff, mmp, mmn, Ep, En, dEp, dEn;
}

TLorentzVector PNeutPol_Polarimeter::PNVect(Int_t ProtonParticleNumber) // Define vectors for p and n in similar manner to properties above
{
  if(ProtonParticleNumber == 1)
    {
      GVp = GV1;
      GV2 = GetTracks()->GetVector(1, Mn);
      GVn = GV2;
      GVpCalc = GV1Rec;
      GVnCalc = GV2Rec;
    }

  if(ProtonParticleNumber == 2)

    {
      GVp = GV2;
      GV1 = GetTracks()->GetVector(0, Mn); // Since we've decided this particle is a neutron, set its mass to Mn
      GVn = GV1; // The neutron vector as measured by the vertex information
      GVpCalc = GV2Rec;
      GVnCalc = GV1Rec;
    }

  return GVp, GVn, GVpCalc, GVnCalc;
}

Double_t PNeutPol_Polarimeter::LabBoost() // Boost particles from lab frame to CM
{
  GVpB = GVp; //Reset the boost vector on each photon
  GVnB = GVn;
  GVpB.Boost(b); // Do boost after p/n identification now
  GVnB.Boost(b);
  ThetapB = (GVpB.Theta()) * TMath::RadToDeg();
  ThetanB = (GVnB.Theta()) * TMath::RadToDeg();
  PhipB = (GVpB.Phi()) * TMath::RadToDeg();
  PhinB = (GVnB.Phi()) * TMath::RadToDeg();

  return ThetapB, ThetanB, PhipB, PhinB;
}

Double_t PNeutPol_Polarimeter::LabScatter()
{
  Thetap = GVp3.Theta() * TMath::RadToDeg(); // Lab frame angles for proton/neutron
  Phip = GVp3.Phi() * TMath::RadToDeg();
  Thetan = GVn3.Theta() * TMath::RadToDeg();
  Phin = GVn3.Phi() * TMath::RadToDeg();
  ThetanCalc = abs(GVnCalc3.Theta() * TMath::RadToDeg());
  PhinCalc = (GVnCalc3.Phi()) * TMath::RadToDeg();

  return Thetan, Phin, Thetap, Phip, ThetanCalc, PhinCalc;
}

//Double_t PNeutPol::WCVertex(TVector3 MeasuredVector, TVector3 ReconstructedVector, double_t ReconstructorZ, double_t MeasuredZ) // Calculate location of WC Z vertex for measured and reconstructed vector
//{
  // ReconstructorZ = Z vertex of the particle we are using to reconstruct track of other, Measured Z = Z vertex of as measured of the particle we're reconstructing
  //lrec = d/(tan(ReconstructedVector.Theta())); // Calculate how far along from interaction point in target the polarimeter interaction point was
  //zWCRec = ReconstructorZ + lrec; // Assume Zp is correct so neutron actually did come from here
  //l = d/tan(MeasuredVector.Theta());
  //zWC = MeasuredZ + l; // The location of the interaction point as determined by the detected proton

  //return zWCRec, zWC;
//}

PNeutPol_Polarimeter::PNeutPol_Polarimeter() // Define a load of histograms to fill
{
  time = new GH1("time", 	"time", 	1400, -700, 700);
  time_cut = new GH1("time_cut", 	"time_cut", 	1400, -700, 700);

  Zp_Vert = new GH1("Zp_Vertex", "Proton Z Vertex Distribution", 300, -150, 150 );
  Zn_Vert = new GH1("Zn_Vertex", "Neutron Z Vertex Distribution", 300, -150, 150 );
  Zp_Vert_DetSum8 = new GH1("Zp_Vertex_DetSum8", "Proton Z Vertex Distribution DetSum8", 300, -150, 150 );
  Zn_Vert_DetSum8 = new GH1("Zn_Vertex_DetSum8", "Neutron Z Vertex Distribution DetSum8", 300, -150, 150 );
  Zp_Vert_DetSum12 = new GH1("Zp_Vertex_DetSum12", "Proton Z Vertex Distribution DetSum12", 300, -150, 150 );
  Zn_Vert_DetSum12 = new GH1("Zn_Vertex_DetSum12", "Neutron Z Vertex Distribution DetSum12", 300, -150, 150 );
  Ekp = new GH1( "Ekp", "Proton Energy Distribution", 100, 0, 500 );
  Ekn = new GH1( "Ekn", "Neutron Energy Distribution", 100, 0, 500 );
  Ekp_DetSum8 = new GH1( "Ekp_DetSum8", "Proton Energy Distribution DetSum8", 100, 0, 500 );
  Ekn_DetSum8 = new GH1( "Ekn_DetSum8", "Neutron Energy Distribution DetSum8", 100, 0, 500 );
  Ekp_DetSum12 = new GH1( "Ekp_DetSum12", "Proton Energy Distribution DetSum12", 100, 0, 500 );
  Ekn_DetSum12 = new GH1( "Ekn_DetSum12", "Neutron Energy Distribution DetSum12", 100, 0, 500 );
  EkSum = new GH1( "Ek Sum", "Particle Energy Sum Distribution", 300, 0, 900 );
  EkSum_DetSum8 = new GH1( "Ek Sum DetSum8", "Particle Energy Sum Distribution DetSum8", 300, 0, 900 );
  EkSum_DetSum12 = new GH1( "Ek Sum DetsSum13", "Particle Energy Sum Distribution DetSum12", 300, 0, 900 );
  Eg = new GH1( "Eg", "Photon Energy Distribution", 200, 100, 1600 );
  Eg_DetSum8 = new GH1( "Eg_DetSum8", "Photon Energy Distribution DetSum8", 200, 100, 1600 );
  Eg_DetSum12 = new GH1( "Eg_DetSum12", "Photon Energy Distribution DetSum 12", 200, 100, 1600 );
  ThetaProt = new GH1( "ThetaProt", " Proton Theta Distribution", 180, 0, 180 );
  ThetaProt_DetSum8 = new GH1( "ThetaProt_DetSum8", " Proton Theta Distribution DetSum8", 180, 0, 180 );
  ThetaProt_DetSum12 = new GH1( "ThetaProt_DetSum12", " Proton Theta Distribution DetSum12", 180, 0, 180 );
  ThetaNeut = new GH1( "ThetaNeut", " Neutron Theta Distribution", 180, 0, 180 );
  ThetaNeut_DetSum8 = new GH1( "ThetaNeut_DetSum8", " Neutron Theta Distribution DetSum8", 180, 0, 180 );
  ThetaNeut_DetSum12 = new GH1( "ThetaNeut_DetSum12", " Neutron Theta Distribution DetSum12", 180, 0, 180 );
  PhiProt = new GH1( "PhiProt", " Proton Phi Distribution", 180, -180, 180 );
  PhiProt_DetSum8 = new GH1( "PhiProt_DetSum8", " Proton Phi Distribution DetSum8", 180, -180, 180 );
  PhiProt_DetSum12 = new GH1( "PhiProt_DetSum12", " Proton Phi Distribution DetSum12", 180, -180, 180 );
  PhiNeut = new GH1( "PhiNeut", " Neutron Phi Distribution", 180, -180, 180 );
  PhiNeut_DetSum8 = new GH1( "PhiNeut_DetSum8", " Neutron Phi Distribution DetSum8", 180, -180, 180 );
  PhiNeut_DetSum12 = new GH1( "PhiNeut_DetSum12", " Neutron Phi Distribution DetSum12", 180, -180, 180 );
  MM_Proton = new GH1("MM_Proton", 	"MM_Proton", 	 	300,   800, 1100);
  MM_Proton_DetSum8 = new GH1("MM_Proton_DetSum8", 	"MM_Proton DetSum8", 	 	300,   800, 1100);
  MM_Proton_DetSum12 = new GH1("MM_Proton_DetSum12", 	"MM_Proton DetSum12", 	 	300,   800, 1100);

  E_dE = new GH2("E_dE", "EdE Plot", 125, 0, 500, 125, 0, 7);
  E_dE_p = new GH2("E_dE_p", "EdE Plot for Protons", 125, 0, 500, 125, 0, 7);
  E_dE_n = new GH2("E_dE_n", "EdE Plot for Neutrons", 125, 0, 500, 125, 0, 7);
  E_dE_DetSum8 = new GH2("E_dE_DetSum8", "EdE Plot DetSum8", 125, 0, 500, 125, 0, 7);
  E_dE_p_DetSum8 = new GH2("E_dE_p_DetSum8", "EdE Plot for Protons DetSum8", 125, 0, 500, 125, 0, 7);
  E_dE_n_DetSum8 = new GH2("E_dE_n_DetSum8", "EdE Plot for Neutrons DetSum8", 125, 0, 500, 125, 0, 7);
  E_dE_DetSum12 = new GH2("E_dE_DetSum12", "EdE Plot DetSum12", 125, 0, 500, 125, 0, 7);
  E_dE_p_DetSum12 = new GH2("E_dE_p_DetSum12", "EdE Plot for Protons DetSum12", 125, 0, 500, 125, 0, 7);
  E_dE_n_DetSum12 = new GH2("E_dE_n_DetSum12", "EdE Plot for Neutrons DetSum12", 125, 0, 500, 125, 0, 7);

  EkEg = new GH2("EkEg", "Ek vs Eg for all Particles", 100, 0, 500, 100, 100, 1600);
  EkEg_p = new GH2("EkEg_p", "Ek vs Eg for all Protons", 100, 0, 500, 100, 100, 1600);
  EkEg_n = new GH2("EkEg_n", "Ek vs Eg for all Neutrons", 100, 0, 500, 100, 100, 1600);
  EkEg_DetSum8 = new GH2("EkEg_DetSum8", "Ek vs Eg for all Particles DetSum8", 100, 0, 500, 100, 100, 1600);
  EkEg_p_DetSum8 = new GH2("EkEg_p_DetSum8", "Ek vs Eg for all Protons DetSum8", 100, 0, 500, 100, 100, 1600);
  EkEg_n_DetSum8 = new GH2("EkEg_n_DetSum8", "Ek vs Eg for all Neutrons DetSum8", 100, 0, 500, 100, 100, 1600);
  EkEg_DetSum12 = new GH2("EkEg_DetSum12", "Ek vs Eg for all Particles DetSum12", 100, 0, 500, 100, 100, 1600);
  EkEg_p_DetSum12 = new GH2("EkEg_p_DetSum12", "Ek vs Eg for all Protons DetSum12", 100, 0, 500, 100, 100, 1600);
  EkEg_n_DetSum12 = new GH2("EkEg_n_DetSum12", "Ek vs Eg for all Neutrons DetSum12", 100, 0, 500, 100, 100, 1600);

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
  MM_Proton->Fill(mmp, TaggerTime);
  E_dE->Fill(Ep, dEp, TaggerTime);
  E_dE->Fill(En, dEn, TaggerTime);
  E_dE_p->Fill(Ep, dEp, TaggerTime);
  E_dE_n->Fill(En, dEn, TaggerTime);
  EkEg->Fill(Ep, EGamma, TaggerTime);
  EkEg->Fill(En, EGamma, TaggerTime);
  EkEg_p->Fill(Ep, EGamma, TaggerTime);
  EkEg_n->Fill(En, EGamma, TaggerTime);

  if (DetectorsSum == 8) {

    Zp_Vert_DetSum8->Fill(Zp, TaggerTime);
    Zn_Vert_DetSum8->Fill(Zn, TaggerTime);
    Ekp_DetSum8->Fill(Ep, TaggerTime);
    Ekn_DetSum8->Fill(En, TaggerTime);
    EkSum_DetSum8->Fill((Ep + En),TaggerTime);
    Eg_DetSum8->Fill(EGamma, TaggerTime);
    ThetaProt_DetSum8->Fill(Thetap, TaggerTime);
    ThetaNeut_DetSum8->Fill(Thetan, TaggerTime);
    PhiProt_DetSum8->Fill(Phip, TaggerTime);
    PhiNeut_DetSum8->Fill(Phin, TaggerTime);
    MM_Proton_DetSum8->Fill(mmp, TaggerTime);
    E_dE_DetSum8->Fill(Ep, dEp, TaggerTime);
    E_dE_DetSum8->Fill(En, dEn, TaggerTime);
    E_dE_p_DetSum8->Fill(Ep, dEp, TaggerTime);
    E_dE_n_DetSum8->Fill(En, dEn, TaggerTime);
    EkEg_DetSum8->Fill(Ep, EGamma, TaggerTime);
    EkEg_DetSum8->Fill(En, EGamma, TaggerTime);
    EkEg_p_DetSum8->Fill(Ep, EGamma, TaggerTime);
    EkEg_n_DetSum8->Fill(En, EGamma, TaggerTime);

  }

  else if (DetectorsSum == 12) {

    Zp_Vert_DetSum12->Fill(Zp, TaggerTime);
    Zn_Vert_DetSum12->Fill(Zn, TaggerTime);
    Ekp_DetSum12->Fill(Ep, TaggerTime);
    Ekn_DetSum12->Fill(En, TaggerTime);
    EkSum_DetSum12->Fill((Ep + En),TaggerTime);
    Eg_DetSum12->Fill(EGamma, TaggerTime);
    ThetaProt_DetSum12->Fill(Thetap, TaggerTime);
    ThetaNeut_DetSum12->Fill(Thetan, TaggerTime);
    PhiProt_DetSum12->Fill(Phip, TaggerTime);
    PhiNeut_DetSum12->Fill(Phin, TaggerTime);
    MM_Proton_DetSum12->Fill(mmp, TaggerTime);
    E_dE_DetSum12->Fill(Ep, dEp, TaggerTime);
    E_dE_DetSum12->Fill(En, dEn, TaggerTime);
    E_dE_p_DetSum12->Fill(Ep, dEp, TaggerTime);
    E_dE_n_DetSum12->Fill(En, dEn, TaggerTime);
    EkEg_DetSum12->Fill(Ep, EGamma, TaggerTime);
    EkEg_DetSum12->Fill(En, EGamma, TaggerTime);
    EkEg_p_DetSum12->Fill(Ep, EGamma, TaggerTime);
    EkEg_n_DetSum12->Fill(En, EGamma, TaggerTime);

  }

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
