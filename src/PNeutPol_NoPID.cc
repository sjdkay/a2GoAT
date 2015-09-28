// GoAT Physics analysis to identify neutrons from deuterium photodisintegration
// Various properties of neutrons/protons identified plotted in histograms
// Main aim is to determine spin polarisation of neutrons

#include "PNeutPol_NoPID.h"

PNeutPol_NoPID::~PNeutPol_NoPID()
{
}

Bool_t	PNeutPol_NoPID::Init()
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

Bool_t	PNeutPol_NoPID::Start()
{
  if(!IsGoATFile())
    {
      cout << "ERROR: Input File is not a GoAT file." << endl;
      return kFALSE;
    }

  SetAsPhysicsFile();

  i = 0; // Integer counter
  d = 54.2; // Distance from centre of target to centre of PID
  NP = 0; // Set number of Protons to 0 before checking
  NPi = 0; // Set number of pions to 0 before checking
  NRoo = 0; // Set number of Rootinos to 0 before checking
  Deut = TLorentzVector (0., 0., 0., 1875.613); // 4-Vector of Deuterium target, assume at rest
  Mn = 939.565; // Mass of neutron
  Mp = 938.272; // Mass of proton

  cout << endl;

  if (GetScalers()->GetNEntries() == 0) // MC Data has no scaler entries so if 0, data gets a flag to denote it as MC
  {
      MCData = kTRUE;
      TRandom2 *rGen = new TRandom2(0); // Define new random generator for use in MC smearing fn
  }

  else if (GetScalers()->GetNEntries() != 0) // This flag is listed as false if the number of scaler entries does not equal 0

  {
      MCData = kFALSE;
  }

  if (MCData == kTRUE) MCHists();

  TraverseValidEvents(); // This loops over each event as in old file and calls ProcessEvent() each loop

  return kTRUE;
}

void	PNeutPol_NoPID::ProcessEvent()
{

  GetEvent(); // Function gets number of tracks/protons/pions e.t.c.
  if (GetTracks()->GetMWPC0Energy(0) == 0 || GetTracks()->GetMWPC0Energy(1) == 0 ) return;
  if (GetTracks()->GetMWPC1Energy(0) == 0 || GetTracks()->GetMWPC1Energy(1) == 0 ) return;
  if (NRoo !=0) return; // Goes to next event if any "rootinos" found
  if (NTrack !=2) return; // Ensures two track event
  InitialVect(); // Function gets vectors of identified tracks and returns them
  InitialProp(); // Function gets initial properties (energy, vertex e.t.c.) of identified tracks
  if ( MCData == kTRUE)
  {
      MCSmearing(); // Smear dE values for MC data
      MCTrue(); // Get MC True data and evaluate it
      pPIDElement = GetDetectorHits()->GetPIDHits(0);
  }
  //for (Int_t i=0; i < NTrack; i++){ // Currently nothing relies upon i!

  for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++){

    // Time = ( GetTagger()->GetTaggedTime(j) - GetTracks()->GetTime(i) ); // maybe move this to AFTER the cuts once the Eg-EpSum loop has been checked?
    TaggerTime = GetTagger()->GetTaggedTime(j); // Get tagged time for event
    EGamma = (GetTagger()->GetTaggedEnergy(j)); // Get Photon energy for event

    //if (i == 0) { // These properties get defined for each photon

    Gamma = TLorentzVector (0., 0., EGamma , EGamma); // 4-Vector of Photon beam
    Gamma3 = Gamma.Vect(); // Convert photon beam 4-vector to 3-vector
    B = (Deut + Gamma).Beta(); // Calculte Beta
    b = TVector3(0., 0., B); // Define boost vector
    mm1 = ((Gamma + Deut) - GV2).M(); // Calculate missing mass of each particle
    mm2 = ((Gamma + Deut) - GV1).M();
    mm1Diff = abs(Mn-mm1); // Look at difference of missing mass from neutron mass
    mm2Diff = abs(Mn-mm2);

    //}

    // Cut on MM before looking at regions, require both MM values to look good
    if ( ((mm1 < 850 || mm1 > 1050) == kTRUE) || (((mm2 < 850 || mm2 > 1050)) == kTRUE)) continue; // If either missing mass is outside of the range, drop out

    FillTime(*GetProtons(),time);
    FillTimeCut(*GetProtons(),time_cut);

    //if (i == 0) { // First loop sees if both are inside the protons region, checks only once per event on FIRST track


    if (mm1Diff < mm2Diff) { //If mm1 closer to nucleon than mm2, particle 2 is the "Good event" i.e. the initial proton

        PNProp(2); //These functions are given the particle number of the proton (2 here)
        PNVect(2); // They return the vectors and properties of the neutron and proton from this

    }

    if (mm1Diff > mm2Diff) { // If mm2 closer to nucleon than mm1, particle 1 is the "Good event"

        PNProp(1);
        PNVect(1);

    }


    // Now do loops which check if only one of the protons is inside the proton banana


    if( Zp > 60 || Zp < -60) continue; // Cut if proton vertex not where we expect it to be

    //mmadj = mmp - Mn; // Look at difference from what the missing mass SHOULD be
    //GVpUnadjE = GVp(3); // Take the initial energy of the proton
    //GVp(3) = GVpUnadjE + mmadj; // Add on the adjustment from the missing mass to the energy of the 4-vector
    mmn = ((Gamma+Deut)-GVn).M(); // Recalculate mmn using the Neutron vector with a corrected mass
    GVnCalc = ((Gamma + Deut) - GVp); // Calculate the original neutron 4-vector from conservation of 4 momenta
    GVp3 = GVp.Vect(); // Generate some 3-vectors from the 4-vectors we have
    GVn3 = GVn.Vect();
    GVnCalc3 = GVnCalc.Vect();

    NeutronEnergy(); // Calculate the neutron energy in two ways and look at difference
    LabBoost(); // Boost particles in lab frame and return results
    LabScatter(); // Work out scattering angle in lab frame and return results
    nFrameScatter(); // Work out neutron scattering in its frame and return results
    WCVertex(); // Calculate Z vertex location in WC from measured and reconstructed vectors
    if (ScattTheta > 90) continue;
    if (zWCRec-zWC > 20 || zWCRec-zWC < -20) continue;
    //if (Cut_proton -> IsInside(En, dEn) == kFALSE) continue; //If the neutron is inside the proton region drop out
    // This has been done as the neutrons in the proton region have potentially not scattered in the PID

    FillHists(); // Fill histograms with data generated

    //}
  }

  // cout << endl; // Use to separate out events

  //}
}

void	PNeutPol_NoPID::ProcessScalerRead()
{
	// Fill Tagger Scalers // Currently this seems to fill the file with loads of "TaggerAccScal" histograms
	//FillScalers(GetTC_scaler_min(),GetTC_scaler_max(),TaggerAccScal); // Don't know if these are needed so cut out for now
}

TCutG*	PNeutPol_NoPID::OpenCutFile(Char_t* filename, Char_t* cutname)
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

Int_t PNeutPol_NoPID::GetEvent() // Gets basic info on particles for event
{
  NTrack = GetTracks()->GetNTracks();
  NP = GetProtons()->GetNParticles();
  NPi = GetChargedPions()->GetNParticles();
  NRoo = GetRootinos()->GetNParticles();
  NTag = GetTagger()->GetNTagged();
  return NTrack, NP, NPi, NRoo, NTag;
}

TLorentzVector PNeutPol_NoPID::InitialVect() // Defines initial vectors
{
  GV1 = GetTracks()->GetVector(0, Mp);
  GV2 = GetTracks()->GetVector(1, Mp); // Set both to have proton mass for now
  return GV1, GV2; // Returns intial 4-vectors for use in later functions
}

Double_t PNeutPol_NoPID::InitialProp() // Defines initial particle properties
{
  Theta1 = (GV1.Theta()) * TMath::RadToDeg();
  Theta2 = (GV2.Theta()) * TMath::RadToDeg();
  z1 = GetTracks()->GetPseudoVertexZ(0);
  z2 = GetTracks()->GetPseudoVertexZ(1);
  E1 = GetTracks()->GetClusterEnergy(0);
  E2 = GetTracks()->GetClusterEnergy(1);
  dE1 = GetTracks()->GetVetoEnergy(0);
  dE2 = GetTracks()->GetVetoEnergy(1);

  return Theta1, Theta2, z1, z2, E1, E2, dE1, dE2; // Returns various quantities used in later functions
}

Double_t PNeutPol_NoPID::MCSmearing() // Smear dE values for MC data to represent Energy resolution of PID
{
  dE1 = rGen.Gaus(GetTracks()->GetVetoEnergy(0) , (0.29*(sqrt(GetTracks()->GetVetoEnergy(0)))));
  dE2 = rGen.Gaus(GetTracks()->GetVetoEnergy(1) , (0.29*(sqrt(GetTracks()->GetVetoEnergy(1)))));

  if (dE1 < 0) dE1 = 0.01;
  if (dE2 < 0) dE2 = 0.01;

  return dE1, dE2;
}

Double_t PNeutPol_NoPID::MCTrue()
{

    //Unsure how to actually extract MC true data from its relevant tree?
    Double_t Example = 5;
    // Get MC True info and do some stuff
    return Example;

}

Double_t PNeutPol_NoPID::PNProp(Int_t ProtonParticleNumber) // Define properties of proton and neutron from particles that correspond to each
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
      pPIDElement = GetDetectorHits()->GetPIDHits(0);
      nPIDElement = GetDetectorHits()->GetPIDHits(1);
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
      pPIDElement = GetDetectorHits()->GetPIDHits(1);
      nPIDElement = GetDetectorHits()->GetPIDHits(0);
    }

  return Zp, Zn, zdiff, mmp, mmn, Ep, En, dEp, dEn, pPIDElement, nPIDElement;
}

TLorentzVector PNeutPol_NoPID::PNVect(Int_t ProtonParticleNumber) // Define vectors for p and n in similar manner to properties above
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

Double_t PNeutPol_NoPID::NeutronEnergy() // Calculate the neutron energy from the reconstructed n vector and from kinematics

{

  EnVectCalc = GVnCalc(3)-Mn;
  EnKinCalc = ((TMath::Power((Mn + Mp),2)) * En)/(2*Mn*Mp*(1-cos(GVnB.Theta())));

  return EnVectCalc;

}

Double_t PNeutPol_NoPID::LabBoost() // Boost particles from lab frame to CM
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

Double_t PNeutPol_NoPID::LabScatter()
{
  Thetap = GVp3.Theta() * TMath::RadToDeg(); // Lab frame angles for proton/neutron
  Phip = GVp3.Phi() * TMath::RadToDeg();
  Thetan = GVn3.Theta() * TMath::RadToDeg();
  Phin = GVn3.Phi() * TMath::RadToDeg();
  ThetanCalc = abs(GVnCalc3.Theta() * TMath::RadToDeg());
  PhinCalc = (GVnCalc3.Phi()) * TMath::RadToDeg();

  return Thetan, Phin, Thetap, Phip, ThetanCalc, PhinCalc;
}

Double_t PNeutPol_NoPID::WCVertex() // Calculate location of WC Z vertex for measured and reconstructed vector
{

  lrec = d/(tan(GVnCalc3.Theta())); // Calculate how far along from interaction point the PID interaction was
  zWCRec = Zp + lrec; // Assume Zp is correct so neutron actually did come from here
  l = d/tan(GVn3.Theta());
  zWC = Zn + l; // The location of the interaction point as determined by the detected proton

  return zWCRec, zWC;
}

TVector3 PNeutPol_NoPID::DefineAxes(TVector3 ProtonVector, TVector3 ReconstuctedNeutronVector)
{
  fZ = (ReconstuctedNeutronVector.Unit()); // Define axes of the plane
  fY = ((Gamma3.Cross(ProtonVector)).Unit());
  fX = ((fY.Cross(fZ)).Unit());

  return fX, fY, fZ;
}

Double_t PNeutPol_NoPID::nFrameScatter()
{
  DefineAxes(GVp3, GVnCalc3);
  ScattZ = fZ.Angle(GVn3); // Gives the angles between the axes defined above and the Scattered Proton vector
  ScattY = fY.Angle(GVn3);
  ScattX = fX.Angle(GVn3);
  ScattTheta = ScattZ * TMath::RadToDeg(); // Get the angle of the scattered particle in frame of initial particle
  ScattPhi = (atan2(cos(ScattY),cos(ScattX)))* TMath::RadToDeg();

  return ScattTheta, ScattPhi;
}

PNeutPol_NoPID::PNeutPol_NoPID() // Define a load of histograms to fill
{
  time = new GH1("time", 	"time", 	1400, -700, 700);
  time_cut = new GH1("time_cut", 	"time_cut", 	1400, -700, 700);

  MM_Proton = new GH1("MM_Proton", 	"MM_Proton", 	 	300,   800, 1100);
  MM_Neutron = new GH1("MM_Neutron", 	"MM_Neutron", 	 	300,  800, 1100);

  //TaggerAccScal = new TH1D("TaggerAccScal","TaggerAccScal",352,0,352);

  Ek = new GH1( "Ek", "Particle Energy Distribution", 100, 0, 500 );
  Eg = new GH1( "Eg", "Photon Energy Distribution", 100, 100, 900 );
  EkSum = new GH1( "Ek Sum", "Particle Energy Sum Distribution", 300, 0, 900 );
  ENeutronVectCalc = new GH1 ( "EnVC", "Neutron Energy Calculated from Reconstructed Vectors", 300, 0, 600);
  ENeutronKinCalc = new GH1 ( "EnKC", "Neutron Energy Calculated from Kinematics", 300, 0, 600);
  ENeutronDiff = new GH1 ("EnDiff", "Difference in Calculated Neutron Energies (ZKin - ZVect)", 250, -100, 900);
  ThetaCM = new GH1( "ThetaCM", "Theta (CM) Distribution", 180, 0, 180 );
  Z_Vert = new GH1("Z_Vertex", "Z Vertex Distribution", 300, -150, 150 );
  Zp_Vert = new GH1("Zp_Vertex", "Proton Z Vertex Distribution", 300, -150, 150 );
  Zn_Vert = new GH1("Zn_Vertex", "Neutron Z Vertex Distribution", 300, -150, 150 );
  Z_WireChamber = new GH1("Z_Wire_Chamber", "Wire Chamber Z Vertex Distribution", 300, -150, 150 );
  Z_WireChamberRec = new GH1("Z_Wire_Chamber_Reconstructed", "Reconstructed Wire Chamber Z Vertex Distribution", 300, -150, 150 );
  Z_WireChamberDifference = new GH1("Z_Wire_Chamber_Difference", "Wire Chamber Z Vertex Difference Distribution", 300, -150, 150 );
  ThetaCMProton = new GH1( "ThetaCMProton", " Proton Theta (CM) Distribution", 180, 0, 180 );
  ThetaCMNeutron = new GH1( "ThetaCMNeutron", "Neutron Theta (CM) Distribution", 180, 0, 180 );
  PhiCMProton = new GH1 ("PhiCMProton", "Proton Phi (CM) Distribution", 36, -180, 180);
  PhiCMNeutron = new GH1 ("PhiCMNeutron", "Neutron Phi (CM) Distribution", 36, -180, 180);

  CM150 = new GH1("CM_150MeV", "Theta (CM) Distribution for photon energies of 150pm50 MeV", 160, 0, 160);
  CM250 = new GH1("CM_250MeV", "Theta (CM) Distribution for photon energies of 250pm50 MeV", 160, 0, 160);
  CM350 = new GH1("CM_350MeV", "Theta (CM) Distribution for photon energies of 350pm50 MeV", 160, 0, 160);
  CM450 = new GH1("CM_450MeV", "Theta (CM) Distribution for photon energies of 450pm50 MeV", 80, 0, 160);
  CM550 = new GH1("CM_550MeV", "Theta (CM) Distribution for photon energies of 550pm50 MeV", 40, 0, 160);

  ThetaScLab =  new GH1( "Theta_Neutron_Lab Frame", "Theta_Neutron_Lab_Frame", 180, 0, 180);
  ThetaSc =  new GH1( "Theta_Scattered", "Scattetred Proton Theta Distribution in Rotated Frame", 180, 0, 180 );
  PhiSc = new GH1( "Phi_Scattered", "Scattetred Proton Phi Distribution in Rotated Frame", 36, -180, 180 );
  PhiScCut = new GH1( "Phi_Scattered_Cut", "Scattetred Proton Phi Distribution in Rotated Frame With Cut on Incident Theta", 36, -180, 180 );
  PhiSc125 = new GH1( "Phi_Scattered_125MeV", "Scattetred Proton Phi Distribution in Rotated Frame for Photon Energies of 125pm25MeV", 36, -180, 180);
  PhiSc175 = new GH1( "Phi_Scattered_175MeV", "Scattetred Proton Phi Distribution in Rotated Frame for Photon Energies of 175pm25MeV", 36, -180, 180);
  PhiSc225 = new GH1( "Phi_Scattered_225MeV", "Scattetred Proton Phi Distribution in Rotated Frame for Photon Energies of 225pm25MeV", 36, -180, 180);
  PhiSc275 = new GH1( "Phi_Scattered_275MeV", "Scattetred Proton Phi Distribution in Rotated Frame for Photon Energies of 275pm25MeV", 36, -180, 180);
  PhiSc325 = new GH1( "Phi_Scattered_325MeV", "Scattetred Proton Phi Distribution in Rotated Frame for Photon Energies of 325pm25MeV", 36, -180, 180);
  PhiSc375 = new GH1( "Phi_Scattered_375MeV", "Scattetred Proton Phi Distribution in Rotated Frame for Photon Energies of 375pm25MeV", 36, -180, 180);
  PhiSc425 = new GH1( "Phi_Scattered_425MeV", "Scattetred Proton Phi Distribution in Rotated Frame for Photon Energies of 425pm25MeV", 36, -180, 180);
  PhiSc475 = new GH1( "Phi_Scattered_475MeV", "Scattetred Proton Phi Distribution in Rotated Frame for Photon Energies of 475pm25MeV", 36, -180, 180);
  PhiSc525 = new GH1( "Phi_Scattered_525MeV", "Scattetred Proton Phi Distribution in Rotated Frame for Photon Energies of 525pm25MeV", 36, -180, 180);
  PhiSc575 = new GH1( "Phi_Scattered_575MeV", "Scattetred Proton Phi Distribution in Rotated Frame for Photon Energies of 575pm25MeV", 36, -180, 180);
  PhiSc125Cut = new GH1( "Phi_Scattered_125MeV_Cut", "Phi_Scattered_125MeV_Cut", 36, -180, 180);
  PhiSc175Cut = new GH1( "Phi_Scattered_175MeV_Cut", "Phi_Scattered_175MeV_Cut", 36, -180, 180);
  PhiSc225Cut = new GH1( "Phi_Scattered_225MeV_Cut", "Phi_Scattered_225MeV_Cut", 36, -180, 180);
  PhiSc275Cut = new GH1( "Phi_Scattered_275MeV_Cut", "Phi_Scattered_275MeV_Cut", 36, -180, 180);
  PhiSc325Cut = new GH1( "Phi_Scattered_325MeV_Cut", "Phi_Scattered_325MeV_Cut", 36, -180, 180);
  PhiSc375Cut = new GH1( "Phi_Scattered_375MeV_Cut", "Phi_Scattered_375MeV_Cut", 36, -180, 180);
  PhiSc425Cut = new GH1( "Phi_Scattered_425MeV_Cut", "Phi_Scattered_425MeV_Cut", 36, -180, 180);
  PhiSc475Cut = new GH1( "Phi_Scattered_475MeV_Cut", "Phi_Scattered_475MeV_Cut", 36, -180, 180);
  PhiSc525Cut = new GH1( "Phi_Scattered_525MeV_Cut", "Phi_Scattered_525MeV_Cut", 36, -180, 180);
  PhiSc575Cut = new GH1( "Phi_Scattered_575MeV_Cut", "Phi_Scattered_575MeV_Cut", 36, -180, 180);

  // PhiScCut = new GH1( "Phi Scattered Cut", "Phi Scattered Cut", 160, 0, 160 );

  E_dE = new GH2("E_dE", "EdE Plot", 150, 0, 500, 150, 0, 7);
  E_dE_Proton = new GH2("E_dE_Proton", "EdE Plot for Protons", 150, 0, 500, 150, 0, 7);
  E_dE_Neutron = new GH2("E_dE_Neutron", "EdE Plot for Neutrons", 150, 0, 500, 150, 0, 7);
  //E_dE_LowPhi = new GH2("E_dE_LowPhi", "EdE Plot for Low Phi Events", 150, 0, 500, 150, 0, 7);
  //E_dE_LowPhi_p = new GH2("E_dE_LowPhi_p", "EdE Plot for Low Phi Events for Protons", 150, 0, 500, 150, 0, 7);
  //E_dE_LowPhi_n = new GH2("E_dE_LowPhi_n", "EdE Plot for Low Phi Events for Neutrons", 150, 0, 500, 150, 0, 7);
  //PhiScLow = new GH1( "Phi_Scattered_Low", "Phi_Scattered", 36, -180, 180);
  //ThetanCM_PhinCM = new GH2 ("ThetaCM_vs_PhiCM_Neutron", "ThetaCM_vs_PhiCM_Neutron", 180, 0, 180, 180, -180, 180);
  //Thetan_Vs_Phin_Lab = new GH2("Theta vs Phi Lab" , "Theta vs Phi Lab", 180, 0, 180, 180, -180, 180);
  // Theta_Vs_Phi = new GH2("Theta vs Phi" , "Theta vs Phi", 180, 0, 50, 180, -180, 180); //These plots aren't showing what it used to, not sure if correct/still useful or not
  //PhiSc_dEn = new GH2("PhiSc_dEn" , "PhiSc_dEn", 180, -180, 180, 180, 0, 10);
  //mmE = new GH2 ("mmE", "Missing_Mass_as_a_Function_of_Energy", 200, 0, 600, 200, 850, 1050);
  //PhidEFixp = new GH2 ("PhidEFixp", "Phi_dE_Distribution", 180, -180, 180, 180, 0, 5);
  //PhidECorrFixp = new GH2 ("PhidECorrFixp", "Phi_dECorr_Distribution", 180, -180, 180, 180, 0, 5);

}

void PNeutPol_NoPID::FillHists()
{
  // Fill histograms with stuff we've calculated
  E_dE->Fill(E1, dE1, TaggerTime);
  E_dE->Fill(E2, dE2, TaggerTime);
  MM_Proton->Fill(mmp, TaggerTime);
  MM_Neutron->Fill(mmn, TaggerTime);
  E_dE_Proton->Fill(Ep, dEp, TaggerTime);
  E_dE_Neutron->Fill(En, dEn, TaggerTime);
  Z_Vert->Fill(Zp, TaggerTime);
  Z_Vert->Fill(Zn, TaggerTime);
  Zp_Vert->Fill(Zp, TaggerTime);
  Zn_Vert->Fill(Zn, TaggerTime);
  Z_WireChamber->Fill(zWC, TaggerTime);
  Z_WireChamberRec->Fill(zWCRec, TaggerTime);
  Z_WireChamberDifference->Fill((zWCRec-zWC), TaggerTime);
  ThetaCMProton->Fill(ThetapB, TaggerTime);
  ThetaCMNeutron->Fill(ThetanB, TaggerTime);
  PhiCMProton->Fill(PhipB, TaggerTime);
  PhiCMNeutron->Fill(PhinB, TaggerTime);
  ThetaCM->Fill(ThetapB, TaggerTime);
  ThetaCM->Fill(ThetanB, TaggerTime);
  EkSum->Fill ((E1+E2), TaggerTime);
  Ek->Fill (E1, TaggerTime);
  Ek->Fill (E2, TaggerTime);
  Eg->Fill((EGamma), TaggerTime);
  ENeutronVectCalc->Fill(EnVectCalc, TaggerTime);
  ENeutronKinCalc->Fill(EnKinCalc, TaggerTime);
  ENeutronDiff->Fill(EnKinCalc - EnVectCalc, TaggerTime);
  ThetaScLab -> Fill(abs(Thetan), TaggerTime);
  //Thetan_Vs_Phin_Lab -> Fill(abs(Thetan), Phin, TaggerTime); // We want Phi of scattered p and theta of scattered p in lab frame?
  //ThetanCM_PhinCM ->Fill(ThetanB, PhinB, TaggerTime);
  ThetaSc -> Fill(ScattTheta, TaggerTime);
  PhiSc -> Fill(ScattPhi, TaggerTime);
  //Theta_Vs_Phi -> Fill(ScattTheta, ScattPhi, TaggerTime);
  //PhiSc_dEn -> Fill(ScattPhi, dEn, TaggerTime);
  //mmE ->Fill(Ep, mmp, TaggerTime);
  //if (((Ep > 82 && Ep < 118) == kTRUE) && ((Thetap > 81 && Thetap < 99) == kTRUE))
  //PhidEFixp->Fill(Phip, (dEp*sin(GVp.Theta())), TaggerTime);
  //if (((Ep > 82 && Ep < 118) == kTRUE) && ((Thetap > 81 && Thetap < 99) == kTRUE))
  //PhidECorrFixp->Fill(Phip, dEp, TaggerTime);

 // if (-10 < ScattPhi && ScattPhi < 10) // Disabled when not needed, these slow down process a lot
  //{
    //  E_dE_LowPhi->Fill(Ep, dEp, TaggerTime);
    //  E_dE_LowPhi->Fill(En, dEn, TaggerTime);
    //  E_dE_LowPhi_p->Fill(Ep, dEp, TaggerTime);
    //  E_dE_LowPhi_n->Fill(En, dEn, TaggerTime);
    //  PhiScLow->Fill(ScattPhi, TaggerTime);
  //}


  if (MCData == kTRUE && nBanana == kTRUE) // Fill some histograms for MC data
    {
            EgMC_In->Fill(EGamma, TaggerTime);
            ThetapMC_In->Fill(Thetap, TaggerTime);
            ThetanMC_In->Fill(Thetan, TaggerTime);
            ThetanMC_Rec_In->Fill(ThetanCalc, TaggerTime);
            EpMC_In->Fill(Ep, TaggerTime);
            EnMC_In->Fill(En, TaggerTime);
            PID_Phi_In->Fill(nPIDElement, Phin, TaggerTime);
            PID_Theta_In->Fill(nPIDElement, Thetan, TaggerTime);
            Thetan_dE_MC_In->Fill(dEn, Thetan, TaggerTime);
            ThetanRec_dE_MC_In->Fill(dEn, ThetanCalc, TaggerTime);
            Phin_dE_MC_In->Fill(dEn, Phin, TaggerTime);
    }

  if (MCData == kTRUE && nBanana == kFALSE)
    {
            EgMC_Out->Fill(EGamma, TaggerTime);
            ThetapMC_Out->Fill(Thetap, TaggerTime);
            ThetanMC_Out->Fill(Thetan, TaggerTime);
            ThetanMC_Rec_Out->Fill(ThetanCalc, TaggerTime);
            EpMC_Out->Fill(Ep, TaggerTime);
            EnMC_Out->Fill(En, TaggerTime);
            PID_Phi_Out->Fill(nPIDElement, Phin , TaggerTime);
            PID_Theta_Out->Fill(nPIDElement, Thetan , TaggerTime);
            Thetan_dE_MC_Out->Fill(dEn, Thetan, TaggerTime);
            ThetanRec_dE_MC_Out->Fill(dEn, ThetanCalc, TaggerTime);
            Phin_dE_MC_Out->Fill(dEn, Phin, TaggerTime);

    }

  // Fill PhiScattering angles and ThetaCM angles for various energy regions
  if ( 100 < EGamma && EGamma < 200) CM150->Fill(ThetapB, TaggerTime);
  if ( 200 < EGamma && EGamma < 300) CM250->Fill(ThetapB, TaggerTime);
  if ( 300 < EGamma && EGamma < 400) CM350->Fill(ThetapB, TaggerTime);
  if ( 400 < EGamma && EGamma < 500) CM450->Fill(ThetapB, TaggerTime);
  if ( 500 < EGamma && EGamma < 600) CM550->Fill(ThetapB, TaggerTime);

  if ( 100 < EGamma && EGamma < 150) PhiSc125->Fill(ScattPhi, TaggerTime);
  if ( 150 < EGamma && EGamma < 200) PhiSc175->Fill(ScattPhi, TaggerTime);
  if ( 200 < EGamma && EGamma < 250) PhiSc225->Fill(ScattPhi, TaggerTime);
  if ( 250 < EGamma && EGamma < 300) PhiSc275->Fill(ScattPhi, TaggerTime);
  if ( 300 < EGamma && EGamma < 350) PhiSc325->Fill(ScattPhi, TaggerTime);
  if ( 350 < EGamma && EGamma < 400) PhiSc375->Fill(ScattPhi, TaggerTime);
  if ( 400 < EGamma && EGamma < 450) PhiSc425->Fill(ScattPhi, TaggerTime);
  if ( 450 < EGamma && EGamma < 500) PhiSc475->Fill(ScattPhi, TaggerTime);
  if ( 500 < EGamma && EGamma < 550) PhiSc525->Fill(ScattPhi, TaggerTime);
  if ( 550 < EGamma && EGamma < 600) PhiSc575->Fill(ScattPhi, TaggerTime);

  if(ThetanCalc > 50 && ThetanCalc < 130){

    PhiScCut -> Fill(ScattPhi, TaggerTime);

    if ( 100 < EGamma && EGamma < 150) PhiSc125Cut->Fill(ScattPhi, TaggerTime);
    if ( 150 < EGamma && EGamma < 200) PhiSc175Cut->Fill(ScattPhi, TaggerTime);
    if ( 200 < EGamma && EGamma < 250) PhiSc225Cut->Fill(ScattPhi, TaggerTime);
    if ( 250 < EGamma && EGamma < 300) PhiSc275Cut->Fill(ScattPhi, TaggerTime);
    if ( 300 < EGamma && EGamma < 350) PhiSc325Cut->Fill(ScattPhi, TaggerTime);
    if ( 350 < EGamma && EGamma < 400) PhiSc375Cut->Fill(ScattPhi, TaggerTime);
    if ( 400 < EGamma && EGamma < 450) PhiSc425Cut->Fill(ScattPhi, TaggerTime);
    if ( 450 < EGamma && EGamma < 500) PhiSc475Cut->Fill(ScattPhi, TaggerTime);
    if ( 500 < EGamma && EGamma < 550) PhiSc525Cut->Fill(ScattPhi, TaggerTime);
    if ( 550 < EGamma && EGamma < 600) PhiSc575Cut->Fill(ScattPhi, TaggerTime);
  }
}

void PNeutPol_NoPID::MCHists()
{
  EgMC_In = new GH1("EgMC_In", "Photon Energy Distribution (MC In)", 100, 100 , 900);
  ThetapMC_In = new GH1("ThetapMC_In", "Proton Theta Distribution (MC In)", 180, 0 , 180);
  ThetanMC_In = new GH1("ThetanMC_In", "Neutron Theta Distribution (MC In)", 180, 0 , 180);
  ThetanMC_Rec_In = new GH1("ThetanMC_Rec_In", "Neutron Reconstructed Theta Distribution (MC In)", 180, 0 , 180);
  EpMC_In = new GH1("EpMC_In", "Proton Energy Distribution (MC In)", 100, 0 , 500);
  EnMC_In = new GH1("EnMC_In", "Neutron Energy Distribution (MC In)", 100, 0 , 500);
  EgMC_Out = new GH1("EgMC_Out", "Photon Energy Distribution (MC Out)", 100, 100 , 900);
  ThetapMC_Out = new GH1("ThetapMC_Out", "Proton Theta Distribution (MC Out)", 180, 0 , 180);
  ThetanMC_Out = new GH1("ThetanMC_Out", "Neutron Theta Distribution (MC Out)", 180, 0 , 180);
  ThetanMC_Rec_Out = new GH1("ThetanMC_Rec_Out", "Neutron Reconstructed Theta Distribution (MC Out)", 180, 0 , 180);
  EpMC_Out = new GH1("EpMC_Out", "Proton Energy Distribution (MC Out)", 100, 0 , 500);
  EnMC_Out = new GH1("EnMC_Out", "Neutron Energy Distribution (MC Out)", 100, 0 , 500);
  PID_Phi_In = new GH2("PID_Phi_In", "PID Element vs Phi for Neutrons (MC In)", 24, 0, 24, 72, -180, 180);
  PID_Phi_Out= new GH2("PID_Phi_Out", "PID Element vs Phi for Neutrons (MC Out)", 24, 0, 24, 72, -180, 180);
  PID_Theta_In = new GH2("PID_Theta_In", "PID Element vs Theta for Neutrons (MC In)", 24, 0, 24, 36, 0, 180);
  PID_Theta_Out = new GH2("PID_Theta_Out", "PID Element vs Theta for Neutrons (MC Out)", 24, 0, 24, 36, 0, 180);
  Thetan_dE_MC_In = new GH2 ("Theta_dE_MC_In", "PID energy vs Neutron Theta Distribution (MC In)", 150, 0, 8, 36, 0, 180);
  ThetanRec_dE_MC_In = new GH2 ("ThetaRec_dE_MC_In", "PID energy vs Neutron Reconstructed Theta Distribution (MC In)", 150, 0, 8, 36, 0, 180);
  Phin_dE_MC_In  = new GH2 ("Phi_dE_MC_In", "PID energy vs Neutron Phi Distribution (MC In)", 150, 0, 8, 72, -180, 180);
  Thetan_dE_MC_Out = new GH2 ("Theta_dE_MC_Out", "PID energy vs Neutron Theta Distribution (MC Out)", 150, 0, 3, 36, 0, 180);
  ThetanRec_dE_MC_Out = new GH2 ("ThetaRec_dE_MC_Out", "PID energy vs Neutron Reconstructed Theta Distribution (MC Out)", 150, 0, 3, 36, 0, 180);
  Phin_dE_MC_Out  = new GH2 ("Phi_dE_MC_Out", "PID energy vs Neutron Phi Distribution (MC Out)", 150, 0, 3, 72, -180, 180);
}

Bool_t	PNeutPol_NoPID::Write(){
  // Write some TH1s - currently none to write so commented out
  // GTreeManager::Write(TaggerAccScal); // This spams the file with 2500+ histograms of "TaggerAccScal" so commented out

  // Write all GH1's easily

  GTreeManager::Write();
}
