// GoAT Physics analysis to identify neutrons from deuterium photodisintegration
// Various properties of neutrons/protons identified plotted in histograms
// Main aim is to determine spin polarisation of neutrons

#include "MC_Analysis_1_Part.h"

MC_Analysis_1_Part::~MC_Analysis_1_Part()
{
}

Bool_t	MC_Analysis_1_Part::Init()
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

Bool_t	MC_Analysis_1_Part::Start()
{
  if(!IsGoATFile())
    {
      cout << "ERROR: Input File is not a GoAT file." << endl;
      return kFALSE;
    }

  SetAsPhysicsFile();

  d = 54.2;
  NP = 0;
  NPi = 0;
  NRoo = 0;
  Deut = TLorentzVector (0., 0., 0., 1875.613); // 4-Vector of Deuterium target, assume at rest
  TRandom2 *rGen = new TRandom2(0);

  Cut_CB_proton = OpenCutFile("configfiles/cuts/CB_DeltaE-E_Proton_27_02_15.root", "Proton");
  Cut_proton = Cut_CB_proton;
  Cut_CB_pion = OpenCutFile("configfiles/cuts/CB_DeltaE-E_Pion_27_02_15.root", "Pion");
  Cut_pion = Cut_CB_pion;
  // Cut_CB_neutron = OpenCutFile("configfiles/cuts/CB_DeltaE-E_Neutron_25_02_15.root", "Neutron");
  // Cut_neutron = Cut_CB_neutron; // This only needs to be here if we get simulation to show where we expect neutrons
  cout << endl;

  TraverseValidEvents(); // This loops over each event as in old file and calls ProcessEvent() each loop

  return kTRUE;
}

void	MC_Analysis_1_Part::ProcessEvent()
{

  GetEvent(); // Function gets number of tracks/protons/pions e.t.c.
  if (NRoo !=0) return; // Goes to next event if any "rootinos" found
  if (NTrack !=1) return; // Ensures two track event
  InitialVect(); // Function gets vectors of identified tracks and returns them
  InitialProp(); // Function gets initial properties (energy, vertex e.t.c.) of identified tracks
  MCSmearing(); // Function smears dE, switch off for real data

  //for (Int_t i=0; i < NTrack; i++){ // Currently nothing relies upon i!

  for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++){

    // Time = ( GetTagger()->GetTaggedTime(j) - GetTracks()->GetTime(i) ); // maybe move this to AFTER the cuts once the Eg-EpSum loop has been checked?
    TaggerTime = GetTagger()->GetTaggedTime(j);
    EGamma = (GetTagger()->GetTaggedEnergy(j));

    // if(((EGamma - (E1 + E2)) > 200) || ((EGamma - (E1 + E2)) < 0)) continue; // Remove this cut?

    //if (i == 0) { // These properties get defined for each photon

    Gamma = TLorentzVector (0., 0., EGamma , EGamma); // 4-Vector of Photon beam
    Gamma3 = Gamma.Vect();
    B = (Deut + Gamma).Beta();
    b = TVector3(0., 0., B);

    //}
    LabBoost();
    FillTime(*GetProtons(),time);
    FillTimeCut(*GetProtons(),time_cut);
    FillHists(); // Fill histograms with data generated

    //}
  }

  // cout << endl; // Use to separate out events

  //}
}

void	MC_Analysis_1_Part::ProcessScalerRead()
{
	// Fill Tagger Scalers // Currently this seems to fill the file with loads of "TaggerAccScal" histograms
	//FillScalers(GetTC_scaler_min(),GetTC_scaler_max(),TaggerAccScal); // Don't know if these are needed so cut out for now
}

TCutG*	MC_Analysis_1_Part::OpenCutFile(Char_t* filename, Char_t* cutname)
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

Int_t MC_Analysis_1_Part::GetEvent()
{
  NTrack = GetTracks()->GetNTracks();
  NP = GetProtons()->GetNParticles();
  NPi = GetChargedPions()->GetNParticles();
  NRoo = GetRootinos()->GetNParticles();
  NTag = GetTagger()->GetNTagged();
  return NTrack, NP, NPi, NRoo, NTag;
}

TLorentzVector MC_Analysis_1_Part::InitialVect()
{
  GV1 = GetTracks()->GetVector(0, 938.272);
  return GV1; // Returns intial 4-vectors for use in later functions
}

Double_t MC_Analysis_1_Part::InitialProp()
{
  Theta1 = (GV1.Theta()) * TMath::RadToDeg();
  Phi1 = (GV1.Phi()) * TMath::RadToDeg();
  z1 = GetTracks()->GetPseudoVertexZ(0);
  E1 = GetTracks()->GetClusterEnergy(0);
  dE1 = GetTracks()->GetVetoEnergy(0);
  dE1Corr = dE1/(sin(GV1.Theta()));
  //PIDElement = GetPIDHitIndex(0)
  return Theta1, Phi1, z1, E1, dE1, dE1Corr, PIDElement; // Returns various quantities used in later functions
}

Double_t MC_Analysis_1_Part::MCSmearing()
{
  dE1 = rGen.Gaus(GetTracks()->GetVetoEnergy(0) , (0.63/(sqrt(GetTracks()->GetVetoEnergy(0)))));

  if (dE1 < 0) dE1 = 0.01;

  dE1Corr = dE1/(sin(GV1.Theta()));

  return dE1, dE1Corr;
}

Double_t MC_Analysis_1_Part::LabBoost()
{
  GV1B = GV1; //Reset the boost vector on each photon
  GV1B.Boost(b); // Do boost after p/n identification now
  ThetaB = (GV1B.Theta()) * TMath::RadToDeg();

  return ThetaB;
}

MC_Analysis_1_Part::MC_Analysis_1_Part()
{
  time = new GH1("time", 	"time", 	1400, -700, 700);
  time_cut = new GH1("time_cut", 	"time_cut", 	1400, -700, 700);

  // TaggerAccScal = new TH1D("TaggerAccScal","TaggerAccScal",352,0,352);

  Ek = new GH1( "Ek", "Particle Energy Distribution", 100, 0, 500 );
  Eg = new GH1( "Eg", "Photon Energy Distribution", 100, 100, 900 );
  Theta = new GH1( "Theta", "Theta (lab) Distribution", 180, 0, 180 );
  ThetaCM = new GH1( "ThetaCM", "Theta (CM) Distribution", 180, 0, 180 );
  Z_Vert = new GH1("Z_Vertex", "Z Vertex Distribution", 300, -150, 150 );

  CM150 = new GH1("CM_150MeV", "Theta (CM) Distribution for photon energies of 150pm50 MeV", 160, 0, 160);
  CM250 = new GH1("CM_250MeV", "Theta (CM) Distribution for photon energies of 250pm50 MeV", 160, 0, 160);
  CM350 = new GH1("CM_350MeV", "Theta (CM) Distribution for photon energies of 350pm50 MeV", 160, 0, 160);
  CM450 = new GH1("CM_450MeV", "Theta (CM) Distribution for photon energies of 450pm50 MeV", 80, 0, 160);
  CM550 = new GH1("CM_550MeV", "Theta (CM) Distribution for photon energies of 550pm50 MeV", 40, 0, 160);

  // PhiScCut = new GH1( "Phi Scattered Cut", "Phi Scattered Cut", 160, 0, 160 );

  E_dE = new GH2("E_dE", "E_dE", 150, 0, 500, 150, 0, 7);
  //ThetaPhi = new GH2("ThetaPhi", "ThetaPhi", 180, 0, 180, 180, -180, 180);
  //ThetaPhiIn = new GH2("ThetaPhiIn", "ThetaPhiIn", 180, 0, 180, 180, -180, 180);
  //ThetaPhiOut = new GH2("ThetaPhiOut", "ThetaPhiOut", 180, 0, 180, 180, -180, 180);
  PhidEFix = new GH2 ("PhidEFix", "Phi_dE_Distribution_at_E=100pm18_&_Theta=90pm9", 180, -180, 180, 180, 0, 5);

}

void MC_Analysis_1_Part::FillHists()
{
  // Fill histograms with stuff we've calculated
  E_dE->Fill(E1, dE1, TaggerTime);
  Z_Vert->Fill(z1, TaggerTime);
  Theta->Fill(Theta1, TaggerTime);
  ThetaCM->Fill(ThetaB, TaggerTime);
  Ek->Fill (E1, TaggerTime);
  Eg->Fill((EGamma), TaggerTime);
  if (((E1 > 82 && E1 < 118) == kTRUE) && ((Theta1 > 81 && Theta1 < 99) == kTRUE)) PhidEFix->Fill(Phi1, dE1, TaggerTime);

  // Fill PhiScattering angles and ThetaCM angles for various energy regions
  if ( 100 < EGamma && EGamma < 200) CM150->Fill(ThetaB, TaggerTime);
  if ( 200 < EGamma && EGamma < 300) CM250->Fill(ThetaB, TaggerTime);
  if ( 300 < EGamma && EGamma < 400) CM350->Fill(ThetaB, TaggerTime);
  if ( 400 < EGamma && EGamma < 500) CM450->Fill(ThetaB, TaggerTime);
  if ( 500 < EGamma && EGamma < 600) CM550->Fill(ThetaB, TaggerTime);

  //ThetaPhi->Fill(Theta1, Phi1, TaggerTime);

  //if (Cut_proton -> IsInside(E1, dE1) == kTRUE) ThetaPhiIn->Fill(Theta1, Phi1, TaggerTime);
  //if (Cut_proton -> IsInside(E1, dE1) == kFALSE) ThetaPhiOut->Fill(Theta1, Phi1, TaggerTime);

}

Bool_t	MC_Analysis_1_Part::Write(){
  // Write some TH1s - currently none to write so commented out
  // GTreeManager::Write(TaggerAccScal); // This spams the file with 2500+ histograms of "TaggerAccScal" so commented out

  // Write all GH1's easily
  GTreeManager::Write();
}
