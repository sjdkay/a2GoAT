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
    PromptLow = -5;
    PromptHigh = 20;
    RandomLow1 = -120;
    RandomHigh1 = -20;
    RandomLow2 = 35;
    RandomHigh2 = 135;
    PvRratio = (PromptHigh - PromptLow)/( (RandomHigh1 - RandomLow1) + (RandomHigh2 - RandomLow2));

    EWidth = 100; //Fix the width of the Energy and CosTheta bins used later
    CosThetaWidth = (2./3.);
    PhiScWidth = 0.8;

    NP = 0; // Set number of Protons to 0 before checking
    NPi = 0; // Set number of pions to 0 before checking
    NRoo = 0; // Set number of Rootinos to 0 before checking
    Mn = 939.565; // Mass of neutron in MeV
    Mp = 938.272; // Mass of proton in MeV
    Md = 1875.613; //Mass of Deuterium in MeV
    Mpi = 139.57018; // Mass of charged pion in MeV
    Deut = TLorentzVector (0., 0., 0., 1875.613); // 4-Vector of Deuterium target, assume at rest
    Pp1 = new TLorentzVector(0.,0.,0.,0.);
    Pp2 = new TLorentzVector(0.,0.,0.,0.);
    Pp1C = new TLorentzVector(0.,0.,0.,0.);
    Pp2C = new TLorentzVector(0.,0.,0.,0.);
    Pbeam = new TLorentzVector(0.,0.,0.,0.);
    PbeamC = new TLorentzVector(0.,0.,0.,0.);
    PtargetC = new TLorentzVector(0.,0.,0.,0.);

    Cut_CB_proton = OpenCutFile("configfiles/cuts/CB_DeltaE-E_Proton_24_4_18.root", "Proton"); // These will need adjusting with new Acqu files
    Cut_proton = Cut_CB_proton;
    Cut_CB_pion = OpenCutFile("configfiles/cuts/CB_DeltaE-E_Pion_29_07_15.root", "Pion");
    Cut_pion = Cut_CB_pion;
    Cut_CB_protonKinGood = OpenCutFile("configfiles/cuts/CB_DeltaE-E_ProtonKinGood_11_05_17.root", "ProtonKinGood"); // These will need adjusting with new Acqu files
    Cut_protonKinGood = Cut_CB_protonKinGood;
    Cut_CB_protonKinBad = OpenCutFile("configfiles/cuts/CB_DeltaE-E_ProtonKinBad_15_12_16.root", "ProtonKinBad");
    Cut_protonKinBad = Cut_CB_protonKinBad;
    cout << endl;

    MCData = MCDataCheck(); //Check scaler entries and deterimine if MC Data or not
    TraverseValidEvents(); // This loops over each event as in old file and calls ProcessEvent() each loop
    BGSub();

    return kTRUE;
}

void	PNeutPol_Polarimeter_Circ::ProcessEvent()
{
    // Set up APLCON
    APLCON::Fit_Settings_t settings = APLCON::Fit_Settings_t::Default;
    settings.MaxIterations = 10;
    settings.DebugLevel = 0;
    APLCON kinfit("EMcons", settings);

    auto EnergyMomentumBalance = [] (const vector< vector<double> >& particles) -> vector<double>
    {
        const TLorentzVector target(0,0,0, 1875.613);
        // assume first particle is beam photon
        double tmpMom;
        double tmpMass;
        tmpMass=0;
        tmpMom=sqrt((particles[0][0]+tmpMass)*(particles[0][0]+tmpMass)-tmpMass*tmpMass);

        TLorentzVector tmpBeam (tmpMom*sin(particles[0][1])*cos(particles[0][2]),tmpMom*sin(particles[0][1])*sin(particles[0][2]),tmpMom*cos(particles[0][1]),particles[0][0]+tmpMass);

        tmpMass=938.272;
        tmpMom=sqrt((particles[1][0]+tmpMass)*(particles[1][0]+tmpMass)-tmpMass*tmpMass);
        TLorentzVector tmpProton (tmpMom*sin(particles[1][1])*cos(particles[1][2]),tmpMom*sin(particles[1][1])*sin(particles[1][2]),tmpMom*cos(particles[1][1]),particles[1][0]+tmpMass);

        tmpMass=939.565;
        tmpMom=sqrt((particles[2][0]+tmpMass)*(particles[2][0]+tmpMass)-tmpMass*tmpMass);
        TLorentzVector tmpNeutron (tmpMom*sin(particles[2][1])*cos(particles[2][2]),tmpMom*sin(particles[2][1])*sin(particles[2][2]),tmpMom*cos(particles[2][1]),particles[2][0]+tmpMass);

        TLorentzVector diff = target + tmpBeam - tmpProton - tmpNeutron;
        return {diff.X(), diff.Y(), diff.Z(), diff.T()};
    };

    NTrack = GetTracks()->GetNTracks();
    NP = GetProtons()->GetNParticles();
    NPi = GetChargedPions()->GetNParticles();
    NRoo = GetRootinos()->GetNParticles();
    NTag = GetTagger()->GetNTagged();
    if (MCData == kFALSE){
        //if (NRoo !=0) return; // Goes to next event if any "rootinos" found, remove this?
    }
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
        if (GetTracks()->GetMWPC1Energy(0) == 0) return; // If no hit in second chamber for p drop out
        if (GetTracks()->GetMWPC0Energy(1) == 0) return; // If no hit in first chamber for n drop out
        if (GetTracks()->GetMWPC1Energy(1) == 0) return; // If no hit in second chamber for p drop out
    }

    // If track 2 only gives signals in MWPC and CB it is the neutron
    else if((Detectors1 == 5) && (Detectors2 == 7))
    {
        Proton1 = kFALSE;
        Proton2 = kTRUE;
        if (GetTracks()->GetMWPC0Energy(1) == 0) return; // If no hit in first chamber for p drop out
        if (GetTracks()->GetMWPC1Energy(1) == 0) return; // If no hit in second chamber for p drop out
        if (GetTracks()->GetMWPC0Energy(0) == 0) return; // If no hit in first chamber for n drop out
        if (GetTracks()->GetMWPC1Energy(0) == 0) return; // If no hit in second chamber for n drop out
    }

    // Drop out on ANY other condition (for now)
    else
    {
        return;
    }

    EventNum = GetEventNumber();

    if (Proton1 == kTRUE)
    {
        GVp = GetTracks()->GetVector(0, Mp);
        GVn = GetTracks()->GetVector(1, Mn);
        Timep = GetTracks()->GetTime(0);
        Timen = GetTracks()->GetTime(1);
        Thp = GetTracks()->GetTheta(0);
        ThpRad = GetTracks()->GetThetaRad(0);
        Thn = GetTracks()->GetTheta(1);
        ThnRad = GetTracks()->GetThetaRad(1);
        Php = GetTracks()->GetPhi(0);
        PhpRad = GetTracks()->GetPhiRad(0);
        Phn = GetTracks()->GetPhi(1);
        Ep = GetTracks()->GetClusterEnergy(0);
        En = GetTracks()->GetClusterEnergy(1);
        dEp = GetTracks()->GetVetoEnergy(0);
        dEn = GetTracks()->GetVetoEnergy(1);
        MWPC0pE = GetTracks()->GetMWPC0Energy(0);
        MWPC1pE = GetTracks()->GetMWPC1Energy(0);
        MWPC0nE = GetTracks()->GetMWPC0Energy(1);
        MWPC1nE = GetTracks()->GetMWPC1Energy(1);
        pVertex.SetXYZ(GetTracks()->GetPseudoVertexX(0), GetTracks()->GetPseudoVertexY(0), GetTracks()->GetPseudoVertexZ(0)); // First particle is proton, second is neutron
        nVertex.SetXYZ(GetTracks()->GetPseudoVertexX(1), GetTracks()->GetPseudoVertexY(1), GetTracks()->GetPseudoVertexZ(1));
        WC13Vectp.SetXYZ(GetTracks()->GetMWPC0PosX(0), GetTracks()->GetMWPC0PosY(0), GetTracks()->GetMWPC0PosZ(0));
        WC23Vectp.SetXYZ(GetTracks()->GetMWPC1PosX(0), GetTracks()->GetMWPC1PosY(0), GetTracks()->GetMWPC1PosZ(0));
        WC13Vectn.SetXYZ(GetTracks()->GetMWPC0PosX(1), GetTracks()->GetMWPC0PosY(1), GetTracks()->GetMWPC0PosZ(1));
        WC23Vectn.SetXYZ(GetTracks()->GetMWPC1PosX(1), GetTracks()->GetMWPC1PosY(1), GetTracks()->GetMWPC1PosZ(1));
    }

    else if (Proton2 == kTRUE)
    {
        GVp = GetTracks()->GetVector(1, Mp);
        GVn = GetTracks()->GetVector(0, Mn);
        Timep = GetTracks()->GetTime(1);
        Timen = GetTracks()->GetTime(0);
        Thp = GetTracks()->GetTheta(1);
        ThpRad = GetTracks()->GetThetaRad(1);
        Thn = GetTracks()->GetTheta(0);
        ThnRad = GetTracks()->GetThetaRad(0);
        Php = GetTracks()->GetPhi(1);
        PhpRad = GetTracks()->GetPhiRad(1);
        Phn = GetTracks()->GetPhi(0);
        Ep = GetTracks()->GetClusterEnergy(1); // Therefore the quantity mmp is the amount of missing mass we see when we do a kinematics calculation USING the proton
        En = GetTracks()->GetClusterEnergy(0);
        dEp = GetTracks()->GetVetoEnergy(1);
        dEn = GetTracks()->GetVetoEnergy(0);
        MWPC0pE = GetTracks()->GetMWPC0Energy(1);
        MWPC1pE = GetTracks()->GetMWPC1Energy(1);
        MWPC0nE = GetTracks()->GetMWPC0Energy(0);
        MWPC1nE = GetTracks()->GetMWPC1Energy(0);
        pVertex.SetXYZ(GetTracks()->GetPseudoVertexX(1), GetTracks()->GetPseudoVertexY(1), GetTracks()->GetPseudoVertexZ(1)); // First particle is neutron, second is proton
        nVertex.SetXYZ(GetTracks()->GetPseudoVertexX(0), GetTracks()->GetPseudoVertexY(0), GetTracks()->GetPseudoVertexZ(0));
        WC13Vectp.SetXYZ(GetTracks()->GetMWPC0PosX(1), GetTracks()->GetMWPC0PosY(1), GetTracks()->GetMWPC0PosZ(1));
        WC23Vectp.SetXYZ(GetTracks()->GetMWPC1PosX(1), GetTracks()->GetMWPC1PosY(1), GetTracks()->GetMWPC1PosZ(1));
        WC13Vectn.SetXYZ(GetTracks()->GetMWPC0PosX(0), GetTracks()->GetMWPC0PosY(0), GetTracks()->GetMWPC0PosZ(0));
        WC23Vectn.SetXYZ(GetTracks()->GetMWPC1PosX(0), GetTracks()->GetMWPC1PosY(0), GetTracks()->GetMWPC1PosZ(0));
    }

    else
    {
        return;
    }

    if( pVertex(2) > 60 || pVertex(2) < -60) return; // Particles selected out from other parts tend to be inside anyway, skip this?
    //if( nVertex(2) > 60 || nVertex(2) < -60) return;
    EpCorr = EpPolCorrect(Ep, Thp); //correct Ep for energy loss in polarimeter

    if(Cut_proton -> IsInside(EpCorr, dEp) == kFALSE) return; // If E loss correct proton is NOT inside p banana drop out
    if(MWPC0pE  < 100) return;
    if(MWPC0nE  < 100) return;
    if(MWPC1pE  < 100) return;
    if(MWPC1nE  < 100) return;
    if((MWPC0nE + MWPC1nE) > 1000) return;

    EpDiff = abs(EpCorr - Ep);

    Pp = sqrt (TMath::Power((Mp + EpCorr),2) - TMath::Power(Mp,2));
    Pn = sqrt (TMath::Power((En + Mn ),2) - TMath::Power(Mn,2));
    GVpCorr = TLorentzVector(Pp*sin(ThpRad)*cos(PhpRad), Pp*sin(ThpRad)*sin(PhpRad), Pp*cos(ThpRad), EpCorr+Mp);

    // MUST FEED IN n PHI IN RAD NOT DEG!
    GVnCorr =  CNeutron4VectorCorr(pVertex(2), GVn, En, Pn , Mn, PhnRad);
    ThetanCorr = (GVnCorr.Theta())*TMath::RadToDeg();
    PhinCorr = (GVnCorr.Phi())*TMath::RadToDeg();

    WC1Phip = (WC13Vectp.Phi());
    WC2Phip = (WC23Vectp.Phi());
    WC1Phin = (WC13Vectn.Phi());
    WC2Phin = (WC23Vectn.Phi());
    WCPhiDiffp = abs(WC1Phip - WC2Phip);
    WCPhiDiffn = abs(WC1Phin - WC2Phin);
    WCDirp = WC23Vectp - WC13Vectp;
    MWPCPLp = WCDirp.Mag();

    GVn3 = GVn.Vect();

    kinfit.LinkVariable("beamF",    beamF.Link(),       beamF.LinkSigma());
    kinfit.LinkVariable("protonF",    protonF.Link(),       protonF.LinkSigma());
    kinfit.LinkVariable("neutronF",    neutronF.Link(),       neutronF.LinkSigma());

    vector<string> all_names = {"beamF", "protonF", "neutronF"};
    kinfit.AddConstraint("EnergyMomentumBalance", all_names, EnergyMomentumBalance);

    for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
    {
        TaggerTime = GetTagger()->GetTaggedTime(j); // Get tagged time for event
        EGamma = (GetTagger()->GetTaggedEnergy(j)); // Get Photon energy for event
        BeamHelicity = GetTrigger()->GetHelicity();
        Gamma = TLorentzVector (0., 0., EGamma , EGamma); // 4-Vector of Photon beam
        B = (Deut + Gamma);
        b = -1*B.BoostVector();

        // Boost p and n to CM frame
        pKinB = pKin;
        pKinB.Boost(b);
        ThetapCM = (pKinB.Theta())*TMath::RadToDeg();
        CosThetapCM = cos (pKinB.Theta());

        RecNeutronEpCorr = CNeutron4VectorKin(GVpCorr);
        MMpEpCorr = RecNeutronEpCorr.M();

        if (((MMpEpCorr < 887.8) == kTRUE) || ((MMpEpCorr > 994.2) == kTRUE)) continue; // Force a missing mass cut

        // Gamma(d,p)n calc n from kinematics
        KinEp = CalcKinEnergy(Thp, EGamma, Md, 0., Mp, Mn); // Calculate kin E of proton assuming pn production
        pKin = CProton4VectorKin(KinEp, ThpRad, PhpRad);
        nKin = CNeutron4VectorKin(pKin);
        ThetanRec = (nKin.Theta()) * TMath::RadToDeg();
        PhinRec = (nKin.Phi()) * TMath::RadToDeg();
        pKin3 = pKin.Vect();
        nKin3 = nKin.Vect();

        nKinB = nKin;
        nKinB.Boost(b);
        ThetanCM = (nKinB.Theta())*TMath::RadToDeg();
//        if (80 > ThetanCM || ThetanCM > 100) continue;

        KinEDiff = KinEp - EpCorr;

        ThetanDiff = ThetanCorr-ThetanRec;
        PhinDiff = PhinCorr-PhinRec;

        TVector3 ScattAngles = ScatteredFrameAngles(nKin3, pKin3 , GVn3 , Gamma);
        ScattThetaRad = ScattAngles(2);
        ScattTheta = ScattThetaRad*TMath::RadToDeg();
        ScattPhiRad = ScattAngles(0);
        ScattPhi = ScattPhiRad*TMath::RadToDeg();
        if (ScattPhiRad > -0.01 && ScattPhiRad < 0.01) continue; // Kill excessively small values
        if (abs(ScattPhiRad) > 3.13) continue;
        Wgt = 1 + ((-0.25546)+(-0.870945*(EGamma/1000))+(1.4198*(TMath::Power((EGamma/1000),2))*cos(ScattPhiRad)));

        DOCAVertex1 = TVector3(0., 0., 0.);
        DOCAVertex2 = TVector3(0., 0., 0.);
        DOCA = Calc_dtfInterDOCA((nKin3.Unit()), (GVn3.Unit()), pVertex, nVertex, DOCAVertex1, DOCAVertex2);
        POCAx = DOCAVertex1.X()-(DOCAVertex1.X()-DOCAVertex2.X())/2.0;
        POCAy = DOCAVertex1.Y()-(DOCAVertex1.Y()-DOCAVertex2.Y())/2.0;
        POCAz = DOCAVertex1.Z()-(DOCAVertex1.Z()-DOCAVertex2.Z())/2.0;
        r = sqrt((TMath::Power(POCAx,2))+(TMath::Power(POCAy,2)));
        POCA = TVector3(POCAx, POCAy, POCAz);

        if( r < 35 ) return; // Ensure POCA is at polarimeter radius
        if(ScattTheta > 70) return;

        tdif = TaggerTime - Timep;

        //Calculate vertex position of interaction from MWPC info only, produces same result as pseudo vertex
        MWPCpDir0.SetXYZ(WC23Vectp(0) - WC13Vectp(0), WC23Vectp(1) - WC13Vectp(1), WC23Vectp(2) - WC13Vectp(2));
        MWPCpDir1.SetXYZ(0, 0, 1);
        MWPCpPerp.SetXYZ(-1*WC13Vectp(0), -1*WC13Vectp(1), -1*WC13Vectp(2));
        num0 = (MWPCpDir0.Dot(MWPCpPerp)) -(MWPCpDir0.Dot(MWPCpDir1))*(MWPCpDir1.Dot(MWPCpPerp));
        denum0 = (MWPCpDir0.Dot(MWPCpDir0))-(MWPCpDir0.Dot(MWPCpDir1))*(MWPCpDir0.Dot(MWPCpDir1));
        tsk = num0/denum0;
        MWPCpCalcVertex.SetXYZ(WC13Vectp(0)+tsk*MWPCpDir0.X(), WC13Vectp(1)+tsk*MWPCpDir0.Y(), WC13Vectp(2)+tsk*MWPCpDir0.Z());
        ZpMWPC = MWPCpCalcVertex.Z();

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

TLorentzVector PNeutPol_Polarimeter_Circ::CNeutron4VectorCorr(Double_t ZVert, TLorentzVector n4Vector, Double_t nE, Double_t MagP, Double_t nMass, Double_t nPhi)
{
    Ncor1 = 0.0886404-0.000555077*ZVert+0.000914921*ZVert*ZVert-7.6616e-06*ZVert*ZVert*ZVert;
    Ncor2 = 0.991;
    Ncor3 = 0.612847+0.153167*ZVert-0.00106208*ZVert*ZVert;
    NcorR = Ncor1+Ncor2*(n4Vector.Theta()*180/acos(-1))+Ncor3*sin(n4Vector.Theta());
    NcorRR = NcorR/180.0*acos(-1); //I.e. convert to Rad

    N4VectCorr =  TLorentzVector (MagP*sin(NcorRR)*cos(nPhi),MagP*sin(NcorRR)*sin(nPhi) , MagP*cos(NcorRR) , nE+nMass);

    return N4VectCorr;
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

PNeutPol_Polarimeter_Circ::PNeutPol_Polarimeter_Circ() // Define a load of histograms to fill
{
    Eg = new TH1D( "Eg", "E_{#gamma} Distribution", 200, 100, 1600 );
    EgPrompt = new TH1D( "EgPrompt", "E_{#gamma} Distribution", 200, 100, 1600 );
    EgRandom = new TH1D( "EgRandom", "E_{#gamma} Distribution", 200, 100, 1600 );

    ThetaSc = new TH1D ("ThetaSc", "#theta_{Sc} Distribution", 200, 0, 4);
    ThetaScPrompt = new TH1D ("ThetaScPrompt", "#theta_{Sc} Distribution", 200, 0, 4);
    ThetaScRandom = new TH1D ("ThetaScRandom", "#theta_{Sc} Distribution", 200, 0, 4);

    RPoca = new TH1D("RPoca", "R_{POCA} Distribution", 200, 0, 200);
    RPocaPrompt = new TH1D("RPocaPrompt", "R_{POCA} Distribution", 200, 0, 200);
    RPocaRandom = new TH1D("RPocaRandom", "R_{POCA} Distribution", 200, 0, 200);

    Thetap = new TH1D ("Thetap", "#theta_{p} Distribution", 200, 0, 4);
    ThetapPrompt = new TH1D ("ThetapPrompt", "#theta_{p} Distribution", 200, 0, 4);
    ThetapRandom = new TH1D ("ThetapRandom", "#theta_{p} Distribution", 200, 0, 4);

    Phip = new TH1D( "Phip", "#phi_{sc} Proton Distribution", 200, -4, 4);
    PhipPrompt = new TH1D( "PhipPrompt", "#phi_{sc} Proton Distribution", 200, -4, 4);
    PhipRandom = new TH1D( "PhipRandom", "#phi_{sc} Proton Distribution", 200, -4, 4);

    EProton = new TH1D( "EProton", "E_{p} Distribution", 200, 0, 400);
    EProtonPrompt = new TH1D( "EProtonPrompt", "E_{p} Distribution", 200, 0, 400);
    EProtonRandom = new TH1D( "EProtonRandom", "E_{p} Distribution", 200, 0, 400);

    MMProton = new TH1D( "MMProton", "MM_{p} Distribution", 200, 0, 2000);
    MMProtonPrompt = new TH1D( "MMProtonPrompt", "MM_{p} Distribution", 200, 0, 2000);
    MMProtonRandom = new TH1D( "MMProtonRandom", "MM_{p} Distribution", 200, 0, 2000);

    EdE = new TH2D ("EdE", "dE_{p} as a function of E_{p}", 200, 0, 1000, 200, 0, 10);
    EdEPrompt = new TH2D ("EdEPrompt", "dE_{p} as a function of E_{p}", 200, 0, 1000, 200, 0, 10);
    EdERandom = new TH2D ("EdERandom", "dE_{p} as a function of E_{p}", 200, 0, 1000, 200, 0, 10);

    ThetaScThetap = new TH2D ("ThetaScThetap", "#theta_{p}(#theta_{Sc})", 200, 0, 3, 200, 0, 3);
    ThetaScThetapPrompt = new TH2D ("ThetaScThetapPrompt", "#theta_{p}(#theta_{Sc})", 200, 0, 3, 200, 0, 3);
    ThetaScThetapRandom = new TH2D ("ThetaScThetapRandom", "#theta_{p}(#theta_{Sc})", 200, 0, 3, 200, 0, 3);

    PhiScThetap = new TH2D ("PhiScThetap", "#theta_{p}(#phi_{Sc})", 200, -4, 4, 200,  0, 3);
    PhiScThetapPrompt = new TH2D ("PhiScThetapPrompt", "#theta_{p}(#phi_{Sc})", 200, -4, 4, 200, 0, 3);
    PhiScThetapRandom = new TH2D ("PhiScThetapRandom", "#theta_{p}(#phi_{Sc})", 200, -4, 4, 200,  0, 3);

    rPocaThetap = new TH2D ("rPocaThetap", "#theta_{p}(r_{POCA})", 200, 0, 200, 200,  0, 3);
    rPocaThetapPrompt = new TH2D ("rPocaThetapPrompt", "#theta_{p}(r_{POCA})", 200, 0, 200, 200, 0, 3);
    rPocaThetapRandom = new TH2D ("rPocaThetapRandom", "#theta_{p}(r_{POCA})", 200, 0, 200, 200, 0, 3);

    ThetaDiffPhiDiff = new TH2D("ThetaDiffPhiDiff", "#phi_{nDiff}(#theta_{nDiff})", 200, -180, 180, 200, -180, 180);
    ThetaDiffPhiDiffPrompt = new TH2D("ThetaDiffPhiDiffPrompt", "#phi_{nDiff}(#theta_{nDiff})", 200, -180, 180, 200, -180, 180);
    ThetaDiffPhiDiffRandom = new TH2D("ThetaDiffPhiDiffRandom", "#phi_{nDiff}(#theta_{nDiff})", 200, -180, 180, 200, -180, 180);

    ThetapEgamma = new TH2D("ThetapEg", "#theta_{p}(E_{#gamma})", 200, 0, 1600, 200, 0, 3);
    ThetapEgammaPrompt = new TH2D("ThetapEgPrompt", "#theta_{p}(E_{#gamma})", 200, 0, 1600, 200, 0, 3);
    ThetapEgammaRandom = new TH2D("ThetapEgRandom", "#theta_{p}(E_{#gamma})", 200, 0, 1600, 200, 0, 3);

    ThetapMMp = new TH2D("ThetapMMp", "#theta_{p}(MM_{p})", 200, 800, 1200, 200, 0, 3);
    ThetapMMpPrompt = new TH2D("ThetapMMpPrompt", "#theta_{p}(MM_{p})", 200, 800, 1200, 200, 0, 3);
    ThetapMMpRandom = new TH2D("ThetapMMpRandom", "#theta_{p}(MM_{p})", 200, 800, 1200, 200, 0, 3);

    ThetanEgamma = new TH2D("ThetanEg", "#theta_{n}(E_{#gamma})", 200, 0, 1600, 200, 0, 3);
    ThetanEgammaPrompt = new TH2D("ThetanEgPrompt", "#theta_{n}(E_{#gamma})", 200, 0, 1600, 200, 0, 3);
    ThetanEgammaRandom = new TH2D("ThetanEgRandom", "#theta_{n}(E_{#gamma})", 200, 0, 1600, 200, 0, 3);

    ThetanMMp = new TH2D("ThetanMMp", "#theta_{n}(MM_{p})", 200, 800, 1200, 200, 0, 3);
    ThetanMMpPrompt = new TH2D("ThetanMMpPrompt", "#theta_{n}(MM_{p})", 200, 800, 1200, 200, 0, 3);
    ThetanMMpRandom = new TH2D("ThetanMMpRandom", "#theta_{n}(MM_{p})", 200, 800, 1200, 200, 0, 3);

    PhiScThetan = new TH2D ("PhiScThetan", "#theta_{n}(#phi_{Sc})", 200, -4, 4, 200,  0, 3);
    PhiScThetanPrompt = new TH2D ("PhiScThetanPrompt", "#theta_{n}(#phi_{Sc})", 200, -4, 4, 200, 0, 3);
    PhiScThetanRandom = new TH2D ("PhiScThetanRandom", "#theta_{n}(#phi_{Sc})", 200, -4, 4, 200,  0, 3);

    PhiScThetapCM = new TH2D ("PhiScThetapCM", "#theta_{pCM}(#phi_{Sc})", 200, -4, 4, 200,  0, 3);
    PhiScThetapCMPrompt = new TH2D ("PhiScThetapCMPrompt", "#theta_{pCM}(#phi_{Sc})", 200, -4, 4, 200, 0, 3);
    PhiScThetapCMRandom = new TH2D ("PhiScThetapCMRandom", "#theta_{pCM}(#phi_{Sc})", 200, -4, 4, 200,  0, 3);

    WeightEg = new TH2D("WeightEg", "Weight(E_{#gamma})", 200, 0, 1600, 200, -1, 1);
    WeightEgPrompt = new TH2D("WeightEgPrompt", "Weight(E_{#gamma})", 200, 0, 1600, 200, -1, 1);
    WeightEgRandom = new TH2D("WeightEgRandom", "Weight(E_{#gamma})", 200, 0, 1600, 200, -1, 1);

    WeightPhiSc = new TH2D("WeightPhiSc", "Weight(#phi_{Sc})", 200, -4, 4, 200, -1, 1);
    WeightPhiScPrompt = new TH2D("WeightPhiScPrompt", "Weight(#phi_{Sc})", 200, -4, 4, 200, -1, 1);
    WeightPhiScRandom = new TH2D("WeightPhiScRandom", "Weight(#phi_{Sc})", 200, -4, 4, 200, -1, 1);

    ThetapThetan = new TH2D ("ThetapThetan", "#theta_{p}(#theta_{n})", 200, 0, 4, 200, 0, 4);
    ThetapThetanPrompt = new TH2D ("ThetapThetanPrompt", "#theta_{p}(#theta_{n})", 200, 0, 4, 200, 0, 4);
    ThetapThetanRandom = new TH2D ("ThetapThetanRandom", "#theta_{p}(#theta_{n})", 200, 0, 4, 200, 0, 4);

    for(Int_t A = 0; A < 6; A++){
        for(Int_t B = 0; B < 3; B++){
            PhiScPosHel[A][B] = new TH1D(Form("PhiSc%iPosHelCM%i", 200+(A*100), B+1), Form("#phi_{Sc} E_{#gamma}%i #pm 50MeV CM%i for +ve Helicity", 200+(A*100), B+1), 10, -4, 4);
            PhiScPosHelPrompt[A][B] = new TH1D(Form("PhiSc%iPosHelCM%iPrompt", 200+(A*100), B+1), Form("#phi_{Sc} E_{#gamma}%i #pm 50MeV CM%i for +ve Helicity", 200+(A*100), B+1), 10, -4, 4);
            PhiScPosHelRandom[A][B] = new TH1D(Form("PhiSc%iPosHelCM%iRandom", 200+(A*100), B+1), Form("#phi_{Sc} E_{#gamma}%i #pm 50MeV CM%i for +ve Helicity", 200+(A*100), B+1), 10, -4, 4);
            PhiScNegHel[A][B] = new TH1D(Form("PhiSc%iNegHelCM%i", 200+(A*100), B+1), Form("#phi_{Sc} E_{#gamma}%i #pm 50MeV CM%i for -ve Helicity", 200+(A*100), B+1), 10, -4, 4);
            PhiScNegHelPrompt[A][B] = new TH1D(Form("PhiSc%iNegHelCM%iPrompt", 200+(A*100), B+1), Form("#phi_{Sc} E_{#gamma}%i #pm 50MeV CM%i for -ve Helicity", 200+(A*100), B+1), 10, -4, 4);
            PhiScNegHelRandom[A][B] = new TH1D(Form("PhiSc%iNegHelCM%iRandom", 200+(A*100), B+1), Form("#phi_{Sc} E_{#gamma}%i #pm 50MeV CM%i for -ve Helicity", 265+(A*70), B+1), 10, -4, 4);
        }
    }


    for(Int_t X = 0; X < 12; X++){
        MMpEgamma[X] = new TH1D(Form("MMp_%iMeV", 225+(X*50)), Form("#M_{p} %i #pm 25MeV", 225+(X*50)), 400, 0, 2000);
        MMpEgammaPrompt[X] = new TH1D(Form("MMp_%iMeVPrompt", 225+(X*50)), Form("MM_{p} %i #pm 25MeV", 225+(X*50)), 400, 0, 2000);
        MMpEgammaRandom[X] = new TH1D(Form("MMp_%iMeVRandom", 225+(X*50)), Form("MM_{p} %i #pm 25MeV", 225+(X*50)), 400, 0, 2000);
    }

}

void PNeutPol_Polarimeter_Circ::FillHists()
{
    if (GHistBGSub::IsPrompt(TaggerTime) == kTRUE){
        EgPrompt->Fill(EGamma);

        ThetaScPrompt->Fill(ScattThetaRad);
        RPocaPrompt->Fill(r);
        ThetapPrompt->Fill(ThpRad);
        PhipPrompt->Fill(PhpRad);
        EProtonPrompt->Fill(Ep);
        MMProtonPrompt->Fill(MMpEpCorr);
        EdEPrompt->Fill(EpCorr, dEp);
        ThetaScThetapPrompt->Fill(ScattThetaRad, ThpRad);
        PhiScThetapPrompt->Fill(ScattPhiRad, ThpRad);
        rPocaThetapPrompt->Fill(r, ThpRad);
        ThetaDiffPhiDiffPrompt->Fill(ThetanDiff, PhinDiff);
        ThetapEgammaPrompt->Fill(EGamma, ThpRad);
        ThetapMMpPrompt->Fill(MMpEpCorr, ThpRad);

        ThetanEgammaPrompt->Fill(EGamma, ThnRad);
        ThetanMMpPrompt->Fill(MMpEpCorr, ThnRad);
        PhiScThetanPrompt->Fill(ScattPhiRad, ThnRad);

        PhiScThetapCMPrompt->Fill(ScattPhiRad, pKinB.Theta());

        WeightEgPrompt->Fill(EGamma, Wgt);
        WeightPhiScPrompt->Fill(ScattPhiRad, Wgt);

        ThetapThetanPrompt->Fill(ThpRad, nKin.Theta());

        for(Int_t d = 0; d < 6; d++){ //Energy bins
            ELow = 200 + (d*EWidth);
            EHigh = 300 + (d*EWidth);
            if( ELow < EGamma && EGamma < EHigh){
                for(Int_t e = 0; e < 3; e++){
                    CosThetaLow = (1./3.) - (e*CosThetaWidth);
                    CosThetaHigh = 1 - (e*CosThetaWidth);
                    if(CosThetaHigh > CosThetapCM && CosThetapCM > CosThetaLow){
                        if(BeamHelicity == kTRUE) PhiScPosHelPrompt[d][e]->Fill(ScattPhiRad);
                        else if(BeamHelicity == kFALSE) PhiScNegHelPrompt[d][e]->Fill(ScattPhiRad);
                    }
                }
            }
        }

        for(Int_t Y =0; Y < 12; Y++){
            double_t EMMLOW = 200 + (Y*50);
            double_t EMMHIGH = 250 + (Y*50);
            if( EMMLOW < EGamma && EGamma < EMMHIGH) MMpEgammaPrompt[Y]->Fill(MMpEpCorr);
        }
    }

    if (GHistBGSub::IsRandom(TaggerTime) == kTRUE){
        EgRandom->Fill(EGamma);

        ThetaScRandom->Fill(ScattThetaRad);
        RPocaRandom->Fill(r);
        ThetapRandom->Fill(ThpRad);
        PhipRandom->Fill(PhpRad);
        EProtonRandom->Fill(Ep);
        MMProtonRandom->Fill(MMpEpCorr);
        EdERandom->Fill(EpCorr, dEp);
        ThetaScThetapRandom->Fill(ScattThetaRad, ThpRad);
        PhiScThetapRandom->Fill(ScattPhiRad, ThpRad);
        rPocaThetapRandom->Fill(r, ThpRad);
        ThetaDiffPhiDiffRandom->Fill(ThetanDiff, PhinDiff);
        ThetapEgammaRandom->Fill(EGamma, ThpRad);
        ThetapMMpRandom->Fill(MMpEpCorr, ThpRad);

        ThetanEgammaRandom->Fill(EGamma, ThnRad);
        ThetanMMpRandom->Fill(MMpEpCorr, ThnRad);
        PhiScThetanRandom->Fill(ScattPhiRad, ThnRad);

        PhiScThetapCMRandom->Fill(ScattPhiRad, pKinB.Theta());

        WeightEgRandom->Fill(EGamma, Wgt);
        WeightPhiScRandom->Fill(ScattPhiRad, Wgt);

        ThetapThetanRandom->Fill(ThpRad, nKin.Theta());

        for(Int_t d = 0; d < 6; d++){ //Energy bins
            ELow = 200 + (d*EWidth);
            EHigh = 300 + (d*EWidth);
            if( ELow < EGamma && EGamma < EHigh){
                for(Int_t e = 0; e < 3; e++){
                    CosThetaLow = (1./3.) - (e*CosThetaWidth);
                    CosThetaHigh = 1 - (e*CosThetaWidth);
                    if(CosThetaHigh > CosThetapCM && CosThetapCM > CosThetaLow){
                        if(BeamHelicity == kTRUE) PhiScPosHelRandom[d][e]->Fill(ScattPhiRad);
                        else if(BeamHelicity == kFALSE) PhiScNegHelRandom[d][e]->Fill(ScattPhiRad);
                    }
                }
            }
        }

        for(Int_t Y =0; Y < 12; Y++){
            double_t EMMLOW = 200 + (Y*50);
            double_t EMMHIGH = 250 + (Y*50);
            if( EMMLOW < EGamma && EGamma < EMMHIGH) MMpEgammaRandom[Y]->Fill(MMpEpCorr);
        }
    }
}

void PNeutPol_Polarimeter_Circ::BGSub(){

    Eg->Add(EgPrompt);
    Eg->Add(EgRandom, -PvRratio);
    ThetaSc->Add(ThetaScPrompt);
    ThetaSc->Add(ThetaScRandom, -PvRratio);
    RPoca->Add(RPocaPrompt);
    RPoca->Add(RPocaRandom, -PvRratio);
    Thetap->Add(ThetapPrompt);
    Thetap->Add(ThetapRandom, -PvRratio);
    Phip->Add(PhipPrompt);
    Phip->Add(PhipRandom, -PvRratio);
    EProton->Add(EProtonPrompt);
    EProton->Add(EProtonRandom, -PvRratio);
    MMProton->Add(MMProtonPrompt);
    MMProton->Add(MMProtonRandom, -PvRratio);

    EdE->Add(EdEPrompt);
    EdE->Add(EdERandom, -PvRratio);
    ThetaScThetap->Add(ThetaScThetapPrompt);
    ThetaScThetap->Add(ThetaScThetapRandom, -PvRratio);
    PhiScThetap->Add(PhiScThetapPrompt);
    PhiScThetap->Add(PhiScThetapRandom, -PvRratio);
    rPocaThetap->Add(rPocaThetapPrompt);
    rPocaThetap->Add(rPocaThetapRandom, -PvRratio);
    ThetaDiffPhiDiff->Add(ThetaDiffPhiDiffPrompt);
    ThetaDiffPhiDiff->Add(ThetaDiffPhiDiffRandom, -PvRratio);
    ThetapEgamma->Add(ThetapEgammaPrompt);
    ThetapEgamma->Add(ThetapEgammaRandom, -PvRratio);
    ThetapMMp->Add(ThetapMMpPrompt);
    ThetapMMp->Add(ThetapMMpRandom, -PvRratio);

    ThetanEgamma->Add(ThetanEgammaPrompt);
    ThetanEgamma->Add(ThetanEgammaRandom, -PvRratio);
    ThetanMMp->Add(ThetanMMpPrompt);
    ThetanMMp->Add(ThetanMMpRandom, -PvRratio);
    PhiScThetan->Add(PhiScThetanPrompt);
    PhiScThetan->Add(PhiScThetanRandom, -PvRratio);

    PhiScThetapCM->Add(PhiScThetapCMPrompt);
    PhiScThetapCM->Add(PhiScThetapCMRandom, -PvRratio);

    WeightEg->Add(WeightEgPrompt);
    WeightEg->Add(WeightEgRandom, -PvRratio);
    WeightPhiSc->Add(WeightPhiScPrompt);
    WeightPhiSc->Add(WeightPhiScRandom, -PvRratio);

    ThetapThetan->Add(ThetapThetanPrompt);
    ThetapThetan->Add(ThetapThetanRandom, -PvRratio);

    for(Int_t E = 0; E < 6; E++){
        for(Int_t F = 0; F < 3; F++){
            PhiScNegHel[E][F]->Add(PhiScNegHelPrompt[E][F]);
            PhiScNegHel[E][F]->Add(PhiScNegHelRandom[E][F], -PvRratio);
            PhiScPosHel[E][F]->Add(PhiScPosHelPrompt[E][F]);
            PhiScPosHel[E][F]->Add(PhiScPosHelRandom[E][F], -PvRratio);
        }
    }


    for(Int_t XY = 0; XY < 12; XY++){
        MMpEgamma[XY]->Add(MMpEgammaPrompt[XY]);
        MMpEgamma[XY]->Add(MMpEgammaRandom[XY], -PvRratio);
    }

}


Bool_t	PNeutPol_Polarimeter_Circ::Write(){
    GTreeManager::Write();
}
