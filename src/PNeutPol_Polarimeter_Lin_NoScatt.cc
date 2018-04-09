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

    PromptLow = -5;
    PromptHigh = 20;
    RandomLow1 = -135;
    RandomHigh1 = -35;
    RandomLow2 = 35;
    RandomHigh2 = 135;
    PvRratio = (PromptHigh - PromptLow)/( (RandomHigh1 - RandomLow1) + (RandomHigh2 - RandomLow2));

    k = 0;
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

    MCData = MCDataCheck(); //Check scaler entries and deterimine if MC Data or not
    TraverseValidEvents(); // This loops over each event as in old file and calls ProcessEvent() each loop
    BGSub();

    return kTRUE;
}

void	PNeutPol_Polarimeter_Lin_NoScatt::ProcessEvent()
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

    EventNumber = GetEventNumber();
    NTrack = GetTracks()->GetNTracks();
    NP = GetProtons()->GetNParticles();
    NPi = GetChargedPions()->GetNParticles();
    NRoo = GetRootinos()->GetNParticles();
    NTag = GetTagger()->GetNTagged();
    if (MCData == kFALSE){
        // if (NRoo !=0) return; // Goes to next event if any "rootinos" found, remove this?
    }
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
        Php = GetTracks()->GetPhi(0);
        PhpRad = GetTracks()->GetPhiRad(0);
        Phn = GetTracks()->GetPhi(1);
        Xp = GetTracks()->GetPseudoVertexX(0);
        Yp = GetTracks()->GetPseudoVertexY(0);
        Zp = GetTracks()->GetPseudoVertexZ(0); // First particle is proton, second neutron
        Zn = GetTracks()->GetPseudoVertexZ(1);
        Ep = GetTracks()->GetClusterEnergy(0);
        En = GetTracks()->GetClusterEnergy(1);
        dEp = GetTracks()->GetVetoEnergy(0);
        dEn = GetTracks()->GetVetoEnergy(1);
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
        Php = GetTracks()->GetPhi(1);
        PhpRad = GetTracks()->GetPhiRad(1);
        Phn = GetTracks()->GetPhi(0);
        Xp = GetTracks()->GetPseudoVertexX(1);
        Yp = GetTracks()->GetPseudoVertexY(1);
        Zp = GetTracks()->GetPseudoVertexZ(1); // First particle is neutron, second is proton
        Zn = GetTracks()->GetPseudoVertexZ(0);
        Ep = GetTracks()->GetClusterEnergy(1); // Therefore the quantity mmp is the amount of missing mass we see when we do a kinematics calculation USING the proton
        En = GetTracks()->GetClusterEnergy(0);
        dEp = GetTracks()->GetVetoEnergy(1);
        dEn = GetTracks()->GetVetoEnergy(0);
    }

    else
    {
        return;
    }

    PhiDiff = abs (Php-Phn);

    if( Zp > 60 || Zp < -60) return; // Particles selected out from other parts tend to be inside anyway, skip this?
    if ( PhiDiff > 195 || PhiDiff < 165) return;

    EpCorr = EpPolCorrect(Ep, Thp); //correct Ep for energy loss in polarimeter

    if(Cut_proton -> IsInside(EpCorr, dEp) == kFALSE) return; // If E loss correct proton is NOT inside p banana drop out

    EpDiff = abs(EpCorr - Ep);

    Pp = sqrt (TMath::Power((Mp + EpCorr),2) - TMath::Power(Mp,2));
    Pn = sqrt (TMath::Power((En + Mn ),2) - TMath::Power(Mn,2));
    GVpCorr = TLorentzVector(Pp*sin(Thp)*cos(Php), Pp*sin(Thp)*sin(Php), Pp*cos(Thp), EpCorr+Mp);

    GVnCorr =  LNeutron4VectorCorr(Zp, GVn, En, Pn , Mn, Phn);
    ThetanCorr = (GVnCorr.Theta())*TMath::RadToDeg();

    GVpCorr3 = GVpCorr.Vect();
    GVnCorr3 = GVnCorr.Vect();
    pVertex = TVector3(Xp, Yp, Zp);

    for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
    {
        TaggerTime = GetTagger()->GetTaggedTime(j); // Get tagged time for event
        EGamma = (GetTagger()->GetTaggedEnergy(j)); // Get Photon energy for event
        BeamHelicity = GetTrigger()->GetHelicity();
        Gamma = TLorentzVector (0., 0., EGamma , EGamma); // 4-Vector of Photon beam
        B = (Deut + Gamma);
        b = -1*B.BoostVector();

        // Gamma(d,p)n calc n from kinematics
        KinEp = CalcKinEnergy(Thp, EGamma, Md, 0., Mp, Mn); // Calculate kin E of proton assuming pn production
        pKin = LProton4VectorKin(KinEp, ThpRad, PhpRad);
        nKin = LNeutron4VectorKin(pKin);
        ThetanRec = (nKin.Theta()) * TMath::RadToDeg();
        PhinRec = (nKin.Phi()) * TMath::RadToDeg();
        pKin3 = pKin.Vect();
        nKin3 = nKin.Vect();

        // Boost p and n to CM frame
        pKinB = pKin;
        pKinB.Boost(b);
        ThetapCM = (pKinB.Theta())*TMath::RadToDeg();
        CosThetapCM = cos (pKinB.Theta());
        nKinB = nKin;
        nKinB.Boost(b);
        ThetanCM = (nKinB.Theta())*TMath::RadToDeg();
        //if (80 > ThetanCM || ThetanCM > 100) continue;

        WCZnRec = 72/tan(RecKinNeutron.Theta());

        // Gamma(n,p)Pi (P detected correct)
        // Assume proton track is proton and "neutron" track is from charged pion
        KinEpPi = CalcKinEnergy(Thp, EGamma, Mn, 0, Mp, Mpi); // Calculate kin E of proton assuming g(n, p) pi
        RecKinProtonPi = LProton4VectorKin(KinEpPi, ThpRad, PhpRad); // Get Proton 4 vector from calculated kin E
        RecKinPion = LPion4VectorKin(RecKinProtonPi); // Get Pion 4 vector from 4 momenta conservation
        ThetaPiRec = (RecKinPion.Theta())*TMath::RadToDeg();
        PhiPiRec = (RecKinPion.Phi())*TMath::RadToDeg();
        ThetaPiRecDiff = ThetaPiRec - ThetanCorr;

        // Gamma(n,p)Pi (Pion detected correct)
        // Assume proton track is pion and "neutron" track is from proton
        KinPi = CalcKinEnergy(Thp, EGamma, Mn, 0, Mpi, Mp); // Calculate kin E of pion
        RecKinPionP = LProton4VectorKin(KinPi, ThpRad, PhpRad); // Get Pion 4 vector from calculated kinE
        RecKinPPi = LPion4VectorKin(RecKinPionP); // Get Proton 4 vector from 4 momenta conservation
        ThetapRec = (RecKinPPi.Theta())*TMath::RadToDeg();
        PhipRec = (RecKinPPi.Phi())*TMath::RadToDeg();
        ThetapRecDiff = ThetapRec - ThetanCorr;


        KinEDiff = KinEp - EpCorr;

        RecProtonEpCorr = LProton4VectorKin(EpCorr, ThpRad, PhpRad);
        RecNeutronEpCorr = LNeutron4VectorKin(RecProtonEpCorr);
        MMpEpCorr = RecNeutronEpCorr.M();
        RecProtonEpCorr3 = RecProtonEpCorr.Vect();
        RecNeutronEpCorr3 = RecNeutronEpCorr.Vect();

        P3Vect = RecKinProton.Vect();
        N3Vect = RecKinNeutron.Vect();
        OpeningAngle = (N3Vect.Angle(GVnCorr3))*TMath::RadToDeg();

        ThetanDiff = abs(ThetanRec - ThetanCorr);
        PhinDiff = abs(PhinRec - Phn);

        TVector3 ScattAngles = ScatteredFrameAngles(RecNeutronEpCorr3, GVpCorr3, GVnCorr3, Gamma);
        ScattTheta = ScattAngles(0); // Theta is 1st component in vector fn returns above
        ScattPhi = ScattAngles(1); // Phi is 2nd component


        if(Cut_protonKinGood -> IsInside(KinEp, dEp) == kFALSE) continue; // If KinE proton is NOT inside p banana drop out
        if (((MMpEpCorr < 800) == kTRUE) || ((MMpEpCorr > 1100) == kTRUE)) continue; // Force a missing mass cut
        if ( (ThetanCorr-ThetanRec) < -15 || (ThetanCorr-ThetanRec) > 15) continue;

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

TLorentzVector PNeutPol_Polarimeter_Lin_NoScatt::LNeutron4VectorCorr(Double_t ZVert, TLorentzVector n4Vector, Double_t nE, Double_t MagP, Double_t nMass, Double_t nPhi)
{
    Ncor1 = 0.0886404-0.000555077*ZVert+0.000914921*ZVert*ZVert-7.6616e-06*ZVert*ZVert*ZVert;
    Ncor2=0.991;
    Ncor3= 0.612847+0.153167*ZVert-0.00106208*ZVert*ZVert;
    NcorR=Ncor1+Ncor2*(n4Vector.Theta()*180/acos(-1))+Ncor3*sin(n4Vector.Theta());
    NcorRR=NcorR/180.0*acos(-1);

    N4VectCorr =  TLorentzVector (MagP*sin(NcorRR)*cos(nPhi),MagP*sin(NcorRR)*sin(nPhi) , MagP*cos(NcorRR) , nE+nMass);

    return N4VectCorr;
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
    OAngle = new GH1 ("OAngle", "Opening Angle between P and N Vectors", 180, 0, 180);
    MMpEpCorrected = new GH1 ("MMpEpCorrected", "Missing mass seen by Proton (E Loss Corrected)", 400, 0, 2000);
    ZpDist = new GH1 ("ZpDist", "Proton Pseudo Z Vertex Distribution", 200, -400, 400);
    ThetanDist = new GH1 ("ThetanDist", "#theta_{n} Distribution", 200, 0, 180);

    E_dE = new GH2 ("E_dE", "EdE Plot With E Loss Adjustment", 100, 0, 500, 100, 0, 5);
    DeutKinPiKin = new GH2 ("DeutKinPiKin", "(#theta_{nRec} - #theta_{n}) vs (#theta_{#pi Rec} - #theta_{n})", 200, -180, 180, 200, -180, 180);

    // MMp across photon E bins
    MMp200300 = new GH1("MMp200300", "Missing mass as seen by Proton (200-300MeV Photon Energy)", 400, 0, 2000);
    MMp300400 = new GH1("MMp300400", "Missing mass as seen by Proton (300-400MeV Photon Energy)", 400, 0, 2000);
    MMp400500 = new GH1("MMp400500", "Missing mass as seen by Proton (400-500MeV Photon Energy)", 400, 0, 2000);
    MMp500600 = new GH1("MMp500600", "Missing mass as seen by Proton (500-600MeV Photon Energy)", 400, 0, 2000);
    MMp600700 = new GH1("MMp600700", "Missing mass as seen by Proton (600-700MeV Photon Energy)", 400, 0, 2000);
    MMp700800 = new GH1("MMp700800", "Missing mass as seen by Proton (700-800MeV Photon Energy)", 400, 0, 2000);
    MMp800900 = new GH1("MMp800900", "Missing mass as seen by Proton (800-900MeV Photon Energy)", 400, 0, 2000);

    // Proton Phi dists across EGamma bins
    Phip415CM1 = new GH1("Phip_415MeVCM1", "#phi_{p} Distribution for E_{#gamma} 415 #pm 5MeV (Cos#theta_{CM}1-0.9)", 10, -180, 180);
    Phip425CM1 = new GH1("Phip_425MeVCM1", "#phi_{p} Distribution for E_{#gamma} 425 #pm 5MeV (Cos#theta_{CM}1-0.9)", 10, -180, 180);
    Phip435CM1 = new GH1("Phip_435MeVCM1", "#phi_{p} Distribution for E_{#gamma} 435 #pm 5MeV (Cos#theta_{CM}1-0.9)", 10, -180, 180);
    Phip445CM1 = new GH1("Phip_445MeVCM1", "#phi_{p} Distribution for E_{#gamma} 445 #pm 5MeV (Cos#theta_{CM}1-0.9)", 10, -180, 180);
    Phip455CM1 = new GH1("Phip_455MeVCM1", "#phi_{p} Distribution for E_{#gamma} 455 #pm 5MeV (Cos#theta_{CM}1-0.9)", 10, -180, 180);
    Phip465CM1 = new GH1("Phip_465MeVCM1", "#phi_{p} Distribution for E_{#gamma} 465 #pm 5MeV (Cos#theta_{CM}1-0.9)", 10, -180, 180);
    Phip475CM1 = new GH1("Phip_475MeVCM1", "#phi_{p} Distribution for E_{#gamma} 475 #pm 5MeV (Cos#theta_{CM}1-0.9)", 10, -180, 180);
    Phip485CM1 = new GH1("Phip_485MeVCM1", "#phi_{p} Distribution for E_{#gamma} 485 #pm 5MeV (Cos#theta_{CM}1-0.9)", 10, -180, 180);
    Phip495CM1 = new GH1("Phip_495MeVCM1", "#phi_{p} Distribution for E_{#gamma} 495 #pm 5MeV (Cos#theta_{CM}1-0.9)", 10, -180, 180);
    Phip505CM1 = new GH1("Phip_505MeVCM1", "#phi_{p} Distribution for E_{#gamma} 505 #pm 5MeV (Cos#theta_{CM}1-0.9)", 10, -180, 180);
    Phip515CM1 = new GH1("Phip_515MeVCM1", "#phi_{p} Distribution for E_{#gamma} 515 #pm 5MeV (Cos#theta_{CM}1-0.9)", 10, -180, 180);
    Phip525CM1 = new GH1("Phip_525MeVCM1", "#phi_{p} Distribution for E_{#gamma} 525 #pm 5MeV (Cos#theta_{CM}1-0.9)", 10, -180, 180);
    Phip535CM1 = new GH1("Phip_535MeVCM1", "#phi_{p} Distribution for E_{#gamma} 535 #pm 5MeV (Cos#theta_{CM}1-0.9)", 10, -180, 180);
    Phip545CM1 = new GH1("Phip_545MeVCM1", "#phi_{p} Distribution for E_{#gamma} 545 #pm 5MeV (Cos#theta_{CM}1-0.9)", 10, -180, 180);
    Phip555CM1 = new GH1("Phip_555MeVCM1", "#phi_{p} Distribution for E_{#gamma} 555 #pm 5MeV (Cos#theta_{CM}1-0.9)", 10, -180, 180);
    Phip565CM1 = new GH1("Phip_565MeVCM1", "#phi_{p} Distribution for E_{#gamma} 565 #pm 5MeV (Cos#theta_{CM}1-0.9)", 10, -180, 180);
    Phip575CM1 = new GH1("Phip_575MeVCM1", "#phi_{p} Distribution for E_{#gamma} 575 #pm 5MeV (Cos#theta_{CM}1-0.9)", 10, -180, 180);
    Phip585CM1 = new GH1("Phip_585MeVCM1", "#phi_{p} Distribution for E_{#gamma} 585 #pm 5MeV (Cos#theta_{CM}1-0.9)", 10, -180, 180);
    Phip595CM1 = new GH1("Phip_595MeVCM1", "#phi_{p} Distribution for E_{#gamma} 595 #pm 5MeV (Cos#theta_{CM}1-0.9)", 10, -180, 180);
    Phip605CM1 = new GH1("Phip_605MeVCM1", "#phi_{p} Distribution for E_{#gamma} 605 #pm 5MeV (Cos#theta_{CM}1-0.9)", 10, -180, 180);
    Phip615CM1 = new GH1("Phip_615MeVCM1", "#phi_{p} Distribution for E_{#gamma} 615 #pm 5MeV (Cos#theta_{CM}1-0.9)", 10, -180, 180);

    // #phi_{p} dists across EGamma bins
    Phip415CM2 = new GH1("Phip_415MeVCM2", "#phi_{p} Distribution for E_{#gamma} 415 #pm 5MeV (Cos#theta_{CM}0.9-0.8)", 10, -180, 180);
    Phip425CM2 = new GH1("Phip_425MeVCM2", "#phi_{p} Distribution for E_{#gamma} 425 #pm 5MeV (Cos#theta_{CM}0.9-0.8)", 10, -180, 180);
    Phip435CM2 = new GH1("Phip_435MeVCM2", "#phi_{p} Distribution for E_{#gamma} 435 #pm 5MeV (Cos#theta_{CM}0.9-0.8)", 10, -180, 180);
    Phip445CM2 = new GH1("Phip_445MeVCM2", "#phi_{p} Distribution for E_{#gamma} 445 #pm 5MeV (Cos#theta_{CM}0.9-0.8)", 10, -180, 180);
    Phip455CM2 = new GH1("Phip_455MeVCM2", "#phi_{p} Distribution for E_{#gamma} 455 #pm 5MeV (Cos#theta_{CM}0.9-0.8)", 10, -180, 180);
    Phip465CM2 = new GH1("Phip_465MeVCM2", "#phi_{p} Distribution for E_{#gamma} 465 #pm 5MeV (Cos#theta_{CM}0.9-0.8)", 10, -180, 180);
    Phip475CM2 = new GH1("Phip_475MeVCM2", "#phi_{p} Distribution for E_{#gamma} 475 #pm 5MeV (Cos#theta_{CM}0.9-0.8)", 10, -180, 180);
    Phip485CM2 = new GH1("Phip_485MeVCM2", "#phi_{p} Distribution for E_{#gamma} 485 #pm 5MeV (Cos#theta_{CM}0.9-0.8)", 10, -180, 180);
    Phip495CM2 = new GH1("Phip_495MeVCM2", "#phi_{p} Distribution for E_{#gamma} 495 #pm 5MeV (Cos#theta_{CM}0.9-0.8)", 10, -180, 180);
    Phip505CM2 = new GH1("Phip_505MeVCM2", "#phi_{p} Distribution for E_{#gamma} 505 #pm 5MeV (Cos#theta_{CM}0.9-0.8)", 10, -180, 180);
    Phip515CM2 = new GH1("Phip_515MeVCM2", "#phi_{p} Distribution for E_{#gamma} 515 #pm 5MeV (Cos#theta_{CM}0.9-0.8)", 10, -180, 180);
    Phip525CM2 = new GH1("Phip_525MeVCM2", "#phi_{p} Distribution for E_{#gamma} 525 #pm 5MeV (Cos#theta_{CM}0.9-0.8)", 10, -180, 180);
    Phip535CM2 = new GH1("Phip_535MeVCM2", "#phi_{p} Distribution for E_{#gamma} 535 #pm 5MeV (Cos#theta_{CM}0.9-0.8)", 10, -180, 180);
    Phip545CM2 = new GH1("Phip_545MeVCM2", "#phi_{p} Distribution for E_{#gamma} 545 #pm 5MeV (Cos#theta_{CM}0.9-0.8)", 10, -180, 180);
    Phip555CM2 = new GH1("Phip_555MeVCM2", "#phi_{p} Distribution for E_{#gamma} 555 #pm 5MeV (Cos#theta_{CM}0.9-0.8)", 10, -180, 180);
    Phip565CM2 = new GH1("Phip_565MeVCM2", "#phi_{p} Distribution for E_{#gamma} 565 #pm 5MeV (Cos#theta_{CM}0.9-0.8)", 10, -180, 180);
    Phip575CM2 = new GH1("Phip_575MeVCM2", "#phi_{p} Distribution for E_{#gamma} 575 #pm 5MeV (Cos#theta_{CM}0.9-0.8)", 10, -180, 180);
    Phip585CM2 = new GH1("Phip_585MeVCM2", "#phi_{p} Distribution for E_{#gamma} 585 #pm 5MeV (Cos#theta_{CM}0.9-0.8)", 10, -180, 180);
    Phip595CM2 = new GH1("Phip_595MeVCM2", "#phi_{p} Distribution for E_{#gamma} 595 #pm 5MeV (Cos#theta_{CM}0.9-0.8)", 10, -180, 180);
    Phip605CM2 = new GH1("Phip_605MeVCM2", "#phi_{p} Distribution for E_{#gamma} 605 #pm 5MeV (Cos#theta_{CM}0.9-0.8)", 10, -180, 180);
    Phip615CM2 = new GH1("Phip_615MeVCM2", "#phi_{p} Distribution for E_{#gamma} 615 #pm 5MeV (Cos#theta_{CM}0.9-0.8)", 10, -180, 180);

    // #phi_{p} dists across EGamma bins
    Phip415CM3 = new GH1("Phip_415MeVCM3", "#phi_{p} Distribution for E_{#gamma} 415 #pm 5MeV (Cos#theta_{CM}0.8-0.7)", 10, -180, 180);
    Phip425CM3 = new GH1("Phip_425MeVCM3", "#phi_{p} Distribution for E_{#gamma} 425 #pm 5MeV (Cos#theta_{CM}0.8-0.7)", 10, -180, 180);
    Phip435CM3 = new GH1("Phip_435MeVCM3", "#phi_{p} Distribution for E_{#gamma} 435 #pm 5MeV (Cos#theta_{CM}0.8-0.7)", 10, -180, 180);
    Phip445CM3 = new GH1("Phip_445MeVCM3", "#phi_{p} Distribution for E_{#gamma} 445 #pm 5MeV (Cos#theta_{CM}0.8-0.7)", 10, -180, 180);
    Phip455CM3 = new GH1("Phip_455MeVCM3", "#phi_{p} Distribution for E_{#gamma} 455 #pm 5MeV (Cos#theta_{CM}0.8-0.7)", 10, -180, 180);
    Phip465CM3 = new GH1("Phip_465MeVCM3", "#phi_{p} Distribution for E_{#gamma} 465 #pm 5MeV (Cos#theta_{CM}0.8-0.7)", 10, -180, 180);
    Phip475CM3 = new GH1("Phip_475MeVCM3", "#phi_{p} Distribution for E_{#gamma} 475 #pm 5MeV (Cos#theta_{CM}0.8-0.7)", 10, -180, 180);
    Phip485CM3 = new GH1("Phip_485MeVCM3", "#phi_{p} Distribution for E_{#gamma} 485 #pm 5MeV (Cos#theta_{CM}0.8-0.7)", 10, -180, 180);
    Phip495CM3 = new GH1("Phip_495MeVCM3", "#phi_{p} Distribution for E_{#gamma} 495 #pm 5MeV (Cos#theta_{CM}0.8-0.7)", 10, -180, 180);
    Phip505CM3 = new GH1("Phip_505MeVCM3", "#phi_{p} Distribution for E_{#gamma} 505 #pm 5MeV (Cos#theta_{CM}0.8-0.7)", 10, -180, 180);
    Phip515CM3 = new GH1("Phip_515MeVCM3", "#phi_{p} Distribution for E_{#gamma} 515 #pm 5MeV (Cos#theta_{CM}0.8-0.7)", 10, -180, 180);
    Phip525CM3 = new GH1("Phip_525MeVCM3", "#phi_{p} Distribution for E_{#gamma} 525 #pm 5MeV (Cos#theta_{CM}0.8-0.7)", 10, -180, 180);
    Phip535CM3 = new GH1("Phip_535MeVCM3", "#phi_{p} Distribution for E_{#gamma} 535 #pm 5MeV (Cos#theta_{CM}0.8-0.7)", 10, -180, 180);
    Phip545CM3 = new GH1("Phip_545MeVCM3", "#phi_{p} Distribution for E_{#gamma} 545 #pm 5MeV (Cos#theta_{CM}0.8-0.7)", 10, -180, 180);
    Phip555CM3 = new GH1("Phip_555MeVCM3", "#phi_{p} Distribution for E_{#gamma} 555 #pm 5MeV (Cos#theta_{CM}0.8-0.7)", 10, -180, 180);
    Phip565CM3 = new GH1("Phip_565MeVCM3", "#phi_{p} Distribution for E_{#gamma} 565 #pm 5MeV (Cos#theta_{CM}0.8-0.7)", 10, -180, 180);
    Phip575CM3 = new GH1("Phip_575MeVCM3", "#phi_{p} Distribution for E_{#gamma} 575 #pm 5MeV (Cos#theta_{CM}0.8-0.7)", 10, -180, 180);
    Phip585CM3 = new GH1("Phip_585MeVCM3", "#phi_{p} Distribution for E_{#gamma} 585 #pm 5MeV (Cos#theta_{CM}0.8-0.7)", 10, -180, 180);
    Phip595CM3 = new GH1("Phip_595MeVCM3", "#phi_{p} Distribution for E_{#gamma} 595 #pm 5MeV (Cos#theta_{CM}0.8-0.7)", 10, -180, 180);
    Phip605CM3 = new GH1("Phip_605MeVCM3", "#phi_{p} Distribution for E_{#gamma} 605 #pm 5MeV (Cos#theta_{CM}0.8-0.7)", 10, -180, 180);
    Phip615CM3 = new GH1("Phip_615MeVCM3", "#phi_{p} Distribution for E_{#gamma} 615 #pm 5MeV (Cos#theta_{CM}0.8-0.7)", 10, -180, 180);

    // #phi_{p} dists across EGamma bins
    Phip415CM4 = new GH1("Phip_415MeVCM4", "#phi_{p} Distribution for E_{#gamma} 415 #pm 5MeV (Cos#theta_{CM}0.7-0.6)", 10, -180, 180);
    Phip425CM4 = new GH1("Phip_425MeVCM4", "#phi_{p} Distribution for E_{#gamma} 425 #pm 5MeV (Cos#theta_{CM}0.7-0.6)", 10, -180, 180);
    Phip435CM4 = new GH1("Phip_435MeVCM4", "#phi_{p} Distribution for E_{#gamma} 435 #pm 5MeV (Cos#theta_{CM}0.7-0.6)", 10, -180, 180);
    Phip445CM4 = new GH1("Phip_445MeVCM4", "#phi_{p} Distribution for E_{#gamma} 445 #pm 5MeV (Cos#theta_{CM}0.7-0.6)", 10, -180, 180);
    Phip455CM4 = new GH1("Phip_455MeVCM4", "#phi_{p} Distribution for E_{#gamma} 455 #pm 5MeV (Cos#theta_{CM}0.7-0.6)", 10, -180, 180);
    Phip465CM4 = new GH1("Phip_465MeVCM4", "#phi_{p} Distribution for E_{#gamma} 465 #pm 5MeV (Cos#theta_{CM}0.7-0.6)", 10, -180, 180);
    Phip475CM4 = new GH1("Phip_475MeVCM4", "#phi_{p} Distribution for E_{#gamma} 475 #pm 5MeV (Cos#theta_{CM}0.7-0.6)", 10, -180, 180);
    Phip485CM4 = new GH1("Phip_485MeVCM4", "#phi_{p} Distribution for E_{#gamma} 485 #pm 5MeV (Cos#theta_{CM}0.7-0.6)", 10, -180, 180);
    Phip495CM4 = new GH1("Phip_495MeVCM4", "#phi_{p} Distribution for E_{#gamma} 495 #pm 5MeV (Cos#theta_{CM}0.7-0.6)", 10, -180, 180);
    Phip505CM4 = new GH1("Phip_505MeVCM4", "#phi_{p} Distribution for E_{#gamma} 505 #pm 5MeV (Cos#theta_{CM}0.7-0.6)", 10, -180, 180);
    Phip515CM4 = new GH1("Phip_515MeVCM4", "#phi_{p} Distribution for E_{#gamma} 515 #pm 5MeV (Cos#theta_{CM}0.7-0.6)", 10, -180, 180);
    Phip525CM4 = new GH1("Phip_525MeVCM4", "#phi_{p} Distribution for E_{#gamma} 525 #pm 5MeV (Cos#theta_{CM}0.7-0.6)", 10, -180, 180);
    Phip535CM4 = new GH1("Phip_535MeVCM4", "#phi_{p} Distribution for E_{#gamma} 535 #pm 5MeV (Cos#theta_{CM}0.7-0.6)", 10, -180, 180);
    Phip545CM4 = new GH1("Phip_545MeVCM4", "#phi_{p} Distribution for E_{#gamma} 545 #pm 5MeV (Cos#theta_{CM}0.7-0.6)", 10, -180, 180);
    Phip555CM4 = new GH1("Phip_555MeVCM4", "#phi_{p} Distribution for E_{#gamma} 555 #pm 5MeV (Cos#theta_{CM}0.7-0.6)", 10, -180, 180);
    Phip565CM4 = new GH1("Phip_565MeVCM4", "#phi_{p} Distribution for E_{#gamma} 565 #pm 5MeV (Cos#theta_{CM}0.7-0.6)", 10, -180, 180);
    Phip575CM4 = new GH1("Phip_575MeVCM4", "#phi_{p} Distribution for E_{#gamma} 575 #pm 5MeV (Cos#theta_{CM}0.7-0.6)", 10, -180, 180);
    Phip585CM4 = new GH1("Phip_585MeVCM4", "#phi_{p} Distribution for E_{#gamma} 585 #pm 5MeV (Cos#theta_{CM}0.7-0.6)", 10, -180, 180);
    Phip595CM4 = new GH1("Phip_595MeVCM4", "#phi_{p} Distribution for E_{#gamma} 595 #pm 5MeV (Cos#theta_{CM}0.7-0.6)", 10, -180, 180);
    Phip605CM4 = new GH1("Phip_605MeVCM4", "#phi_{p} Distribution for E_{#gamma} 605 #pm 5MeV (Cos#theta_{CM}0.7-0.6)", 10, -180, 180);
    Phip615CM4 = new GH1("Phip_615MeVCM4", "#phi_{p} Distribution for E_{#gamma} 615 #pm 5MeV (Cos#theta_{CM}0.7-0.6)", 10, -180, 180);

    // #phi_{p} dists across EGamma bins
    Phip415CM5 = new GH1("Phip_415MeVCM5", "#phi_{p} Distribution for E_{#gamma} 415 #pm 5MeV (Cos#theta_{CM}0.6-0.5)", 10, -180, 180);
    Phip425CM5 = new GH1("Phip_425MeVCM5", "#phi_{p} Distribution for E_{#gamma} 425 #pm 5MeV (Cos#theta_{CM}0.6-0.5)", 10, -180, 180);
    Phip435CM5 = new GH1("Phip_435MeVCM5", "#phi_{p} Distribution for E_{#gamma} 435 #pm 5MeV (Cos#theta_{CM}0.6-0.5)", 10, -180, 180);
    Phip445CM5 = new GH1("Phip_445MeVCM5", "#phi_{p} Distribution for E_{#gamma} 445 #pm 5MeV (Cos#theta_{CM}0.6-0.5)", 10, -180, 180);
    Phip455CM5 = new GH1("Phip_455MeVCM5", "#phi_{p} Distribution for E_{#gamma} 455 #pm 5MeV (Cos#theta_{CM}0.6-0.5)", 10, -180, 180);
    Phip465CM5 = new GH1("Phip_465MeVCM5", "#phi_{p} Distribution for E_{#gamma} 465 #pm 5MeV (Cos#theta_{CM}0.6-0.5)", 10, -180, 180);
    Phip475CM5 = new GH1("Phip_475MeVCM5", "#phi_{p} Distribution for E_{#gamma} 475 #pm 5MeV (Cos#theta_{CM}0.6-0.5)", 10, -180, 180);
    Phip485CM5 = new GH1("Phip_485MeVCM5", "#phi_{p} Distribution for E_{#gamma} 485 #pm 5MeV (Cos#theta_{CM}0.6-0.5)", 10, -180, 180);
    Phip495CM5 = new GH1("Phip_495MeVCM5", "#phi_{p} Distribution for E_{#gamma} 495 #pm 5MeV (Cos#theta_{CM}0.6-0.5)", 10, -180, 180);
    Phip505CM5 = new GH1("Phip_505MeVCM5", "#phi_{p} Distribution for E_{#gamma} 505 #pm 5MeV (Cos#theta_{CM}0.6-0.5)", 10, -180, 180);
    Phip515CM5 = new GH1("Phip_515MeVCM5", "#phi_{p} Distribution for E_{#gamma} 515 #pm 5MeV (Cos#theta_{CM}0.6-0.5)", 10, -180, 180);
    Phip525CM5 = new GH1("Phip_525MeVCM5", "#phi_{p} Distribution for E_{#gamma} 525 #pm 5MeV (Cos#theta_{CM}0.6-0.5)", 10, -180, 180);
    Phip535CM5 = new GH1("Phip_535MeVCM5", "#phi_{p} Distribution for E_{#gamma} 535 #pm 5MeV (Cos#theta_{CM}0.6-0.5)", 10, -180, 180);
    Phip545CM5 = new GH1("Phip_545MeVCM5", "#phi_{p} Distribution for E_{#gamma} 545 #pm 5MeV (Cos#theta_{CM}0.6-0.5)", 10, -180, 180);
    Phip555CM5 = new GH1("Phip_555MeVCM5", "#phi_{p} Distribution for E_{#gamma} 555 #pm 5MeV (Cos#theta_{CM}0.6-0.5)", 10, -180, 180);
    Phip565CM5 = new GH1("Phip_565MeVCM5", "#phi_{p} Distribution for E_{#gamma} 565 #pm 5MeV (Cos#theta_{CM}0.6-0.5)", 10, -180, 180);
    Phip575CM5 = new GH1("Phip_575MeVCM5", "#phi_{p} Distribution for E_{#gamma} 575 #pm 5MeV (Cos#theta_{CM}0.6-0.5)", 10, -180, 180);
    Phip585CM5 = new GH1("Phip_585MeVCM5", "#phi_{p} Distribution for E_{#gamma} 585 #pm 5MeV (Cos#theta_{CM}0.6-0.5)", 10, -180, 180);
    Phip595CM5 = new GH1("Phip_595MeVCM5", "#phi_{p} Distribution for E_{#gamma} 595 #pm 5MeV (Cos#theta_{CM}0.6-0.5)", 10, -180, 180);
    Phip605CM5 = new GH1("Phip_605MeVCM5", "#phi_{p} Distribution for E_{#gamma} 605 #pm 5MeV (Cos#theta_{CM}0.6-0.5)", 10, -180, 180);
    Phip615CM5 = new GH1("Phip_615MeVCM5", "#phi_{p} Distribution for E_{#gamma} 615 #pm 5MeV (Cos#theta_{CM}0.6-0.5)", 10, -180, 180);

    // #phi_{p} dists across EGamma bins
    Phip415CM6 = new GH1("Phip_415MeVCM6", "#phi_{p} Distribution for E_{#gamma} 415 #pm 5MeV (Cos#theta_{CM}0.5-0.4)", 10, -180, 180);
    Phip425CM6 = new GH1("Phip_425MeVCM6", "#phi_{p} Distribution for E_{#gamma} 425 #pm 5MeV (Cos#theta_{CM}0.5-0.4)", 10, -180, 180);
    Phip435CM6 = new GH1("Phip_435MeVCM6", "#phi_{p} Distribution for E_{#gamma} 435 #pm 5MeV (Cos#theta_{CM}0.5-0.4)", 10, -180, 180);
    Phip445CM6 = new GH1("Phip_445MeVCM6", "#phi_{p} Distribution for E_{#gamma} 445 #pm 5MeV (Cos#theta_{CM}0.5-0.4)", 10, -180, 180);
    Phip455CM6 = new GH1("Phip_455MeVCM6", "#phi_{p} Distribution for E_{#gamma} 455 #pm 5MeV (Cos#theta_{CM}0.5-0.4)", 10, -180, 180);
    Phip465CM6 = new GH1("Phip_465MeVCM6", "#phi_{p} Distribution for E_{#gamma} 465 #pm 5MeV (Cos#theta_{CM}0.5-0.4)", 10, -180, 180);
    Phip475CM6 = new GH1("Phip_475MeVCM6", "#phi_{p} Distribution for E_{#gamma} 475 #pm 5MeV (Cos#theta_{CM}0.5-0.4)", 10, -180, 180);
    Phip485CM6 = new GH1("Phip_485MeVCM6", "#phi_{p} Distribution for E_{#gamma} 485 #pm 5MeV (Cos#theta_{CM}0.5-0.4)", 10, -180, 180);
    Phip495CM6 = new GH1("Phip_495MeVCM6", "#phi_{p} Distribution for E_{#gamma} 495 #pm 5MeV (Cos#theta_{CM}0.5-0.4)", 10, -180, 180);
    Phip505CM6 = new GH1("Phip_505MeVCM6", "#phi_{p} Distribution for E_{#gamma} 505 #pm 5MeV (Cos#theta_{CM}0.5-0.4)", 10, -180, 180);
    Phip515CM6 = new GH1("Phip_515MeVCM6", "#phi_{p} Distribution for E_{#gamma} 515 #pm 5MeV (Cos#theta_{CM}0.5-0.4)", 10, -180, 180);
    Phip525CM6 = new GH1("Phip_525MeVCM6", "#phi_{p} Distribution for E_{#gamma} 525 #pm 5MeV (Cos#theta_{CM}0.5-0.4)", 10, -180, 180);
    Phip535CM6 = new GH1("Phip_535MeVCM6", "#phi_{p} Distribution for E_{#gamma} 535 #pm 5MeV (Cos#theta_{CM}0.5-0.4)", 10, -180, 180);
    Phip545CM6 = new GH1("Phip_545MeVCM6", "#phi_{p} Distribution for E_{#gamma} 545 #pm 5MeV (Cos#theta_{CM}0.5-0.4)", 10, -180, 180);
    Phip555CM6 = new GH1("Phip_555MeVCM6", "#phi_{p} Distribution for E_{#gamma} 555 #pm 5MeV (Cos#theta_{CM}0.5-0.4)", 10, -180, 180);
    Phip565CM6 = new GH1("Phip_565MeVCM6", "#phi_{p} Distribution for E_{#gamma} 565 #pm 5MeV (Cos#theta_{CM}0.5-0.4)", 10, -180, 180);
    Phip575CM6 = new GH1("Phip_575MeVCM6", "#phi_{p} Distribution for E_{#gamma} 575 #pm 5MeV (Cos#theta_{CM}0.5-0.4)", 10, -180, 180);
    Phip585CM6 = new GH1("Phip_585MeVCM6", "#phi_{p} Distribution for E_{#gamma} 585 #pm 5MeV (Cos#theta_{CM}0.5-0.4)", 10, -180, 180);
    Phip595CM6 = new GH1("Phip_595MeVCM6", "#phi_{p} Distribution for E_{#gamma} 595 #pm 5MeV (Cos#theta_{CM}0.5-0.4)", 10, -180, 180);
    Phip605CM6 = new GH1("Phip_605MeVCM6", "#phi_{p} Distribution for E_{#gamma} 605 #pm 5MeV (Cos#theta_{CM}0.5-0.4)", 10, -180, 180);
    Phip615CM6 = new GH1("Phip_615MeVCM6", "#phi_{p} Distribution for E_{#gamma} 615 #pm 5MeV (Cos#theta_{CM}0.5-0.4)", 10, -180, 180);

    // #phi_{p} dists across EGamma bins
    Phip415CM7 = new GH1("Phip_415MeVCM7", "#phi_{p} Distribution for E_{#gamma} 415 #pm 5MeV (Cos#theta_{CM}0.4-0.3)", 10, -180, 180);
    Phip425CM7 = new GH1("Phip_425MeVCM7", "#phi_{p} Distribution for E_{#gamma} 425 #pm 5MeV (Cos#theta_{CM}0.4-0.3)", 10, -180, 180);
    Phip435CM7 = new GH1("Phip_435MeVCM7", "#phi_{p} Distribution for E_{#gamma} 435 #pm 5MeV (Cos#theta_{CM}0.4-0.3)", 10, -180, 180);
    Phip445CM7 = new GH1("Phip_445MeVCM7", "#phi_{p} Distribution for E_{#gamma} 445 #pm 5MeV (Cos#theta_{CM}0.4-0.3)", 10, -180, 180);
    Phip455CM7 = new GH1("Phip_455MeVCM7", "#phi_{p} Distribution for E_{#gamma} 455 #pm 5MeV (Cos#theta_{CM}0.4-0.3)", 10, -180, 180);
    Phip465CM7 = new GH1("Phip_465MeVCM7", "#phi_{p} Distribution for E_{#gamma} 465 #pm 5MeV (Cos#theta_{CM}0.4-0.3)", 10, -180, 180);
    Phip475CM7 = new GH1("Phip_475MeVCM7", "#phi_{p} Distribution for E_{#gamma} 475 #pm 5MeV (Cos#theta_{CM}0.4-0.3)", 10, -180, 180);
    Phip485CM7 = new GH1("Phip_485MeVCM7", "#phi_{p} Distribution for E_{#gamma} 485 #pm 5MeV (Cos#theta_{CM}0.4-0.3)", 10, -180, 180);
    Phip495CM7 = new GH1("Phip_495MeVCM7", "#phi_{p} Distribution for E_{#gamma} 495 #pm 5MeV (Cos#theta_{CM}0.4-0.3)", 10, -180, 180);
    Phip505CM7 = new GH1("Phip_505MeVCM7", "#phi_{p} Distribution for E_{#gamma} 505 #pm 5MeV (Cos#theta_{CM}0.4-0.3)", 10, -180, 180);
    Phip515CM7 = new GH1("Phip_515MeVCM7", "#phi_{p} Distribution for E_{#gamma} 515 #pm 5MeV (Cos#theta_{CM}0.4-0.3)", 10, -180, 180);
    Phip525CM7 = new GH1("Phip_525MeVCM7", "#phi_{p} Distribution for E_{#gamma} 525 #pm 5MeV (Cos#theta_{CM}0.4-0.3)", 10, -180, 180);
    Phip535CM7 = new GH1("Phip_535MeVCM7", "#phi_{p} Distribution for E_{#gamma} 535 #pm 5MeV (Cos#theta_{CM}0.4-0.3)", 10, -180, 180);
    Phip545CM7 = new GH1("Phip_545MeVCM7", "#phi_{p} Distribution for E_{#gamma} 545 #pm 5MeV (Cos#theta_{CM}0.4-0.3)", 10, -180, 180);
    Phip555CM7 = new GH1("Phip_555MeVCM7", "#phi_{p} Distribution for E_{#gamma} 555 #pm 5MeV (Cos#theta_{CM}0.4-0.3)", 10, -180, 180);
    Phip565CM7 = new GH1("Phip_565MeVCM7", "#phi_{p} Distribution for E_{#gamma} 565 #pm 5MeV (Cos#theta_{CM}0.4-0.3)", 10, -180, 180);
    Phip575CM7 = new GH1("Phip_575MeVCM7", "#phi_{p} Distribution for E_{#gamma} 575 #pm 5MeV (Cos#theta_{CM}0.4-0.3)", 10, -180, 180);
    Phip585CM7 = new GH1("Phip_585MeVCM7", "#phi_{p} Distribution for E_{#gamma} 585 #pm 5MeV (Cos#theta_{CM}0.4-0.3)", 10, -180, 180);
    Phip595CM7 = new GH1("Phip_595MeVCM7", "#phi_{p} Distribution for E_{#gamma} 595 #pm 5MeV (Cos#theta_{CM}0.4-0.3)", 10, -180, 180);
    Phip605CM7 = new GH1("Phip_605MeVCM7", "#phi_{p} Distribution for E_{#gamma} 605 #pm 5MeV (Cos#theta_{CM}0.4-0.3)", 10, -180, 180);
    Phip615CM7 = new GH1("Phip_615MeVCM7", "#phi_{p} Distribution for E_{#gamma} 615 #pm 5MeV (Cos#theta_{CM}0.4-0.3)", 10, -180, 180);

    // #phi_{p} dists across EGamma bins
    Phip415CM8 = new GH1("Phip_415MeVCM8", "#phi_{p} Distribution for E_{#gamma} 415 #pm 5MeV (Cos#theta_{CM}0.3-0.2)", 10, -180, 180);
    Phip425CM8 = new GH1("Phip_425MeVCM8", "#phi_{p} Distribution for E_{#gamma} 425 #pm 5MeV (Cos#theta_{CM}0.3-0.2)", 10, -180, 180);
    Phip435CM8 = new GH1("Phip_435MeVCM8", "#phi_{p} Distribution for E_{#gamma} 435 #pm 5MeV (Cos#theta_{CM}0.3-0.2)", 10, -180, 180);
    Phip445CM8 = new GH1("Phip_445MeVCM8", "#phi_{p} Distribution for E_{#gamma} 445 #pm 5MeV (Cos#theta_{CM}0.3-0.2)", 10, -180, 180);
    Phip455CM8 = new GH1("Phip_455MeVCM8", "#phi_{p} Distribution for E_{#gamma} 455 #pm 5MeV (Cos#theta_{CM}0.3-0.2)", 10, -180, 180);
    Phip465CM8 = new GH1("Phip_465MeVCM8", "#phi_{p} Distribution for E_{#gamma} 465 #pm 5MeV (Cos#theta_{CM}0.3-0.2)", 10, -180, 180);
    Phip475CM8 = new GH1("Phip_475MeVCM8", "#phi_{p} Distribution for E_{#gamma} 475 #pm 5MeV (Cos#theta_{CM}0.3-0.2)", 10, -180, 180);
    Phip485CM8 = new GH1("Phip_485MeVCM8", "#phi_{p} Distribution for E_{#gamma} 485 #pm 5MeV (Cos#theta_{CM}0.3-0.2)", 10, -180, 180);
    Phip495CM8 = new GH1("Phip_495MeVCM8", "#phi_{p} Distribution for E_{#gamma} 495 #pm 5MeV (Cos#theta_{CM}0.3-0.2)", 10, -180, 180);
    Phip505CM8 = new GH1("Phip_505MeVCM8", "#phi_{p} Distribution for E_{#gamma} 505 #pm 5MeV (Cos#theta_{CM}0.3-0.2)", 10, -180, 180);
    Phip515CM8 = new GH1("Phip_515MeVCM8", "#phi_{p} Distribution for E_{#gamma} 515 #pm 5MeV (Cos#theta_{CM}0.3-0.2)", 10, -180, 180);
    Phip525CM8 = new GH1("Phip_525MeVCM8", "#phi_{p} Distribution for E_{#gamma} 525 #pm 5MeV (Cos#theta_{CM}0.3-0.2)", 10, -180, 180);
    Phip535CM8 = new GH1("Phip_535MeVCM8", "#phi_{p} Distribution for E_{#gamma} 535 #pm 5MeV (Cos#theta_{CM}0.3-0.2)", 10, -180, 180);
    Phip545CM8 = new GH1("Phip_545MeVCM8", "#phi_{p} Distribution for E_{#gamma} 545 #pm 5MeV (Cos#theta_{CM}0.3-0.2)", 10, -180, 180);
    Phip555CM8 = new GH1("Phip_555MeVCM8", "#phi_{p} Distribution for E_{#gamma} 555 #pm 5MeV (Cos#theta_{CM}0.3-0.2)", 10, -180, 180);
    Phip565CM8 = new GH1("Phip_565MeVCM8", "#phi_{p} Distribution for E_{#gamma} 565 #pm 5MeV (Cos#theta_{CM}0.3-0.2)", 10, -180, 180);
    Phip575CM8 = new GH1("Phip_575MeVCM8", "#phi_{p} Distribution for E_{#gamma} 575 #pm 5MeV (Cos#theta_{CM}0.3-0.2)", 10, -180, 180);
    Phip585CM8 = new GH1("Phip_585MeVCM8", "#phi_{p} Distribution for E_{#gamma} 585 #pm 5MeV (Cos#theta_{CM}0.3-0.2)", 10, -180, 180);
    Phip595CM8 = new GH1("Phip_595MeVCM8", "#phi_{p} Distribution for E_{#gamma} 595 #pm 5MeV (Cos#theta_{CM}0.3-0.2)", 10, -180, 180);
    Phip605CM8 = new GH1("Phip_605MeVCM8", "#phi_{p} Distribution for E_{#gamma} 605 #pm 5MeV (Cos#theta_{CM}0.3-0.2)", 10, -180, 180);
    Phip615CM8 = new GH1("Phip_615MeVCM8", "#phi_{p} Distribution for E_{#gamma} 615 #pm 5MeV (Cos#theta_{CM}0.3-0.2)", 10, -180, 180);

    // #phi_{p} dists across EGamma bins
    Phip415CM9 = new GH1("Phip_415MeVCM9", "#phi_{p} Distribution for E_{#gamma} 415 #pm 5MeV (Cos#theta_{CM}0.2-0.1)", 10, -180, 180);
    Phip425CM9 = new GH1("Phip_425MeVCM9", "#phi_{p} Distribution for E_{#gamma} 425 #pm 5MeV (Cos#theta_{CM}0.2-0.1)", 10, -180, 180);
    Phip435CM9 = new GH1("Phip_435MeVCM9", "#phi_{p} Distribution for E_{#gamma} 435 #pm 5MeV (Cos#theta_{CM}0.2-0.1)", 10, -180, 180);
    Phip445CM9 = new GH1("Phip_445MeVCM9", "#phi_{p} Distribution for E_{#gamma} 445 #pm 5MeV (Cos#theta_{CM}0.2-0.1)", 10, -180, 180);
    Phip455CM9 = new GH1("Phip_455MeVCM9", "#phi_{p} Distribution for E_{#gamma} 455 #pm 5MeV (Cos#theta_{CM}0.2-0.1)", 10, -180, 180);
    Phip465CM9 = new GH1("Phip_465MeVCM9", "#phi_{p} Distribution for E_{#gamma} 465 #pm 5MeV (Cos#theta_{CM}0.2-0.1)", 10, -180, 180);
    Phip475CM9 = new GH1("Phip_475MeVCM9", "#phi_{p} Distribution for E_{#gamma} 475 #pm 5MeV (Cos#theta_{CM}0.2-0.1)", 10, -180, 180);
    Phip485CM9 = new GH1("Phip_485MeVCM9", "#phi_{p} Distribution for E_{#gamma} 485 #pm 5MeV (Cos#theta_{CM}0.2-0.1)", 10, -180, 180);
    Phip495CM9 = new GH1("Phip_495MeVCM9", "#phi_{p} Distribution for E_{#gamma} 495 #pm 5MeV (Cos#theta_{CM}0.2-0.1)", 10, -180, 180);
    Phip505CM9 = new GH1("Phip_505MeVCM9", "#phi_{p} Distribution for E_{#gamma} 505 #pm 5MeV (Cos#theta_{CM}0.2-0.1)", 10, -180, 180);
    Phip515CM9 = new GH1("Phip_515MeVCM9", "#phi_{p} Distribution for E_{#gamma} 515 #pm 5MeV (Cos#theta_{CM}0.2-0.1)", 10, -180, 180);
    Phip525CM9 = new GH1("Phip_525MeVCM9", "#phi_{p} Distribution for E_{#gamma} 525 #pm 5MeV (Cos#theta_{CM}0.2-0.1)", 10, -180, 180);
    Phip535CM9 = new GH1("Phip_535MeVCM9", "#phi_{p} Distribution for E_{#gamma} 535 #pm 5MeV (Cos#theta_{CM}0.2-0.1)", 10, -180, 180);
    Phip545CM9 = new GH1("Phip_545MeVCM9", "#phi_{p} Distribution for E_{#gamma} 545 #pm 5MeV (Cos#theta_{CM}0.2-0.1)", 10, -180, 180);
    Phip555CM9 = new GH1("Phip_555MeVCM9", "#phi_{p} Distribution for E_{#gamma} 555 #pm 5MeV (Cos#theta_{CM}0.2-0.1)", 10, -180, 180);
    Phip565CM9 = new GH1("Phip_565MeVCM9", "#phi_{p} Distribution for E_{#gamma} 565 #pm 5MeV (Cos#theta_{CM}0.2-0.1)", 10, -180, 180);
    Phip575CM9 = new GH1("Phip_575MeVCM9", "#phi_{p} Distribution for E_{#gamma} 575 #pm 5MeV (Cos#theta_{CM}0.2-0.1)", 10, -180, 180);
    Phip585CM9 = new GH1("Phip_585MeVCM9", "#phi_{p} Distribution for E_{#gamma} 585 #pm 5MeV (Cos#theta_{CM}0.2-0.1)", 10, -180, 180);
    Phip595CM9 = new GH1("Phip_595MeVCM9", "#phi_{p} Distribution for E_{#gamma} 595 #pm 5MeV (Cos#theta_{CM}0.2-0.1)", 10, -180, 180);
    Phip605CM9 = new GH1("Phip_605MeVCM9", "#phi_{p} Distribution for E_{#gamma} 605 #pm 5MeV (Cos#theta_{CM}0.2-0.1)", 10, -180, 180);
    Phip615CM9 = new GH1("Phip_615MeVCM9", "#phi_{p} Distribution for E_{#gamma} 615 #pm 5MeV (Cos#theta_{CM}0.2-0.1)", 10, -180, 180);

    // #phi_{p} dists across EGamma bins
    Phip415CM10 = new GH1("Phip_415MeVCM10", "#phi_{p} Distribution for E_{#gamma} 415 #pm 5MeV (Cos#theta_{CM}0.1-0)", 10, -180, 180);
    Phip425CM10 = new GH1("Phip_425MeVCM10", "#phi_{p} Distribution for E_{#gamma} 425 #pm 5MeV (Cos#theta_{CM}0.1-0)", 10, -180, 180);
    Phip435CM10 = new GH1("Phip_435MeVCM10", "#phi_{p} Distribution for E_{#gamma} 435 #pm 5MeV (Cos#theta_{CM}0.1-0)", 10, -180, 180);
    Phip445CM10 = new GH1("Phip_445MeVCM10", "#phi_{p} Distribution for E_{#gamma} 445 #pm 5MeV (Cos#theta_{CM}0.1-0)", 10, -180, 180);
    Phip455CM10 = new GH1("Phip_455MeVCM10", "#phi_{p} Distribution for E_{#gamma} 455 #pm 5MeV (Cos#theta_{CM}0.1-0)", 10, -180, 180);
    Phip465CM10 = new GH1("Phip_465MeVCM10", "#phi_{p} Distribution for E_{#gamma} 465 #pm 5MeV (Cos#theta_{CM}0.1-0)", 10, -180, 180);
    Phip475CM10 = new GH1("Phip_475MeVCM10", "#phi_{p} Distribution for E_{#gamma} 475 #pm 5MeV (Cos#theta_{CM}0.1-0)", 10, -180, 180);
    Phip485CM10 = new GH1("Phip_485MeVCM10", "#phi_{p} Distribution for E_{#gamma} 485 #pm 5MeV (Cos#theta_{CM}0.1-0)", 10, -180, 180);
    Phip495CM10 = new GH1("Phip_495MeVCM10", "#phi_{p} Distribution for E_{#gamma} 495 #pm 5MeV (Cos#theta_{CM}0.1-0)", 10, -180, 180);
    Phip505CM10 = new GH1("Phip_505MeVCM10", "#phi_{p} Distribution for E_{#gamma} 505 #pm 5MeV (Cos#theta_{CM}0.1-0)", 10, -180, 180);
    Phip515CM10 = new GH1("Phip_515MeVCM10", "#phi_{p} Distribution for E_{#gamma} 515 #pm 5MeV (Cos#theta_{CM}0.1-0)", 10, -180, 180);
    Phip525CM10 = new GH1("Phip_525MeVCM10", "#phi_{p} Distribution for E_{#gamma} 525 #pm 5MeV (Cos#theta_{CM}0.1-0)", 10, -180, 180);
    Phip535CM10 = new GH1("Phip_535MeVCM10", "#phi_{p} Distribution for E_{#gamma} 535 #pm 5MeV (Cos#theta_{CM}0.1-0)", 10, -180, 180);
    Phip545CM10 = new GH1("Phip_545MeVCM10", "#phi_{p} Distribution for E_{#gamma} 545 #pm 5MeV (Cos#theta_{CM}0.1-0)", 10, -180, 180);
    Phip555CM10 = new GH1("Phip_555MeVCM10", "#phi_{p} Distribution for E_{#gamma} 555 #pm 5MeV (Cos#theta_{CM}0.1-0)", 10, -180, 180);
    Phip565CM10 = new GH1("Phip_565MeVCM10", "#phi_{p} Distribution for E_{#gamma} 565 #pm 5MeV (Cos#theta_{CM}0.1-0)", 10, -180, 180);
    Phip575CM10 = new GH1("Phip_575MeVCM10", "#phi_{p} Distribution for E_{#gamma} 575 #pm 5MeV (Cos#theta_{CM}0.1-0)", 10, -180, 180);
    Phip585CM10 = new GH1("Phip_585MeVCM10", "#phi_{p} Distribution for E_{#gamma} 585 #pm 5MeV (Cos#theta_{CM}0.1-0)", 10, -180, 180);
    Phip595CM10 = new GH1("Phip_595MeVCM10", "#phi_{p} Distribution for E_{#gamma} 595 #pm 5MeV (Cos#theta_{CM}0.1-0)", 10, -180, 180);
    Phip605CM10 = new GH1("Phip_605MeVCM10", "#phi_{p} Distribution for E_{#gamma} 605 #pm 5MeV (Cos#theta_{CM}0.1-0)", 10, -180, 180);
    Phip615CM10 = new GH1("Phip_615MeVCM10", "#phi_{p} Distribution for E_{#gamma} 615 #pm 5MeV (Cos#theta_{CM}0.1-0)", 10, -180, 180);

    Phip415CM11 = new GH1("Phip_415MeVCM11", "#phi_{p} Distribution for E_{#gamma} 415 #pm 5MeV (Cos#theta_{CM}0-(-0.1))", 10, -180, 180);
    Phip425CM11 = new GH1("Phip_425MeVCM11", "#phi_{p} Distribution for E_{#gamma} 425 #pm 5MeV (Cos#theta_{CM}0-(-0.1))", 10, -180, 180);
    Phip435CM11 = new GH1("Phip_435MeVCM11", "#phi_{p} Distribution for E_{#gamma} 435 #pm 5MeV (Cos#theta_{CM}0-(-0.1))", 10, -180, 180);
    Phip445CM11 = new GH1("Phip_445MeVCM11", "#phi_{p} Distribution for E_{#gamma} 445 #pm 5MeV (Cos#theta_{CM}0-(-0.1))", 10, -180, 180);
    Phip455CM11 = new GH1("Phip_455MeVCM11", "#phi_{p} Distribution for E_{#gamma} 455 #pm 5MeV (Cos#theta_{CM}0-(-0.1))", 10, -180, 180);
    Phip465CM11 = new GH1("Phip_465MeVCM11", "#phi_{p} Distribution for E_{#gamma} 465 #pm 5MeV (Cos#theta_{CM}0-(-0.1))", 10, -180, 180);
    Phip475CM11 = new GH1("Phip_475MeVCM11", "#phi_{p} Distribution for E_{#gamma} 475 #pm 5MeV (Cos#theta_{CM}0-(-0.1))", 10, -180, 180);
    Phip485CM11 = new GH1("Phip_485MeVCM11", "#phi_{p} Distribution for E_{#gamma} 485 #pm 5MeV (Cos#theta_{CM}0-(-0.1))", 10, -180, 180);
    Phip495CM11 = new GH1("Phip_495MeVCM11", "#phi_{p} Distribution for E_{#gamma} 495 #pm 5MeV (Cos#theta_{CM}0-(-0.1))", 10, -180, 180);
    Phip505CM11 = new GH1("Phip_505MeVCM11", "#phi_{p} Distribution for E_{#gamma} 505 #pm 5MeV (Cos#theta_{CM}0-(-0.1))", 10, -180, 180);
    Phip515CM11 = new GH1("Phip_515MeVCM11", "#phi_{p} Distribution for E_{#gamma} 515 #pm 5MeV (Cos#theta_{CM}0-(-0.1))", 10, -180, 180);
    Phip525CM11 = new GH1("Phip_525MeVCM11", "#phi_{p} Distribution for E_{#gamma} 525 #pm 5MeV (Cos#theta_{CM}0-(-0.1))", 10, -180, 180);
    Phip535CM11 = new GH1("Phip_535MeVCM11", "#phi_{p} Distribution for E_{#gamma} 535 #pm 5MeV (Cos#theta_{CM}0-(-0.1))", 10, -180, 180);
    Phip545CM11 = new GH1("Phip_545MeVCM11", "#phi_{p} Distribution for E_{#gamma} 545 #pm 5MeV (Cos#theta_{CM}0-(-0.1))", 10, -180, 180);
    Phip555CM11 = new GH1("Phip_555MeVCM11", "#phi_{p} Distribution for E_{#gamma} 555 #pm 5MeV (Cos#theta_{CM}0-(-0.1))", 10, -180, 180);
    Phip565CM11 = new GH1("Phip_565MeVCM11", "#phi_{p} Distribution for E_{#gamma} 565 #pm 5MeV (Cos#theta_{CM}0-(-0.1))", 10, -180, 180);
    Phip575CM11 = new GH1("Phip_575MeVCM11", "#phi_{p} Distribution for E_{#gamma} 575 #pm 5MeV (Cos#theta_{CM}0-(-0.1))", 10, -180, 180);
    Phip585CM11 = new GH1("Phip_585MeVCM11", "#phi_{p} Distribution for E_{#gamma} 585 #pm 5MeV (Cos#theta_{CM}0-(-0.1))", 10, -180, 180);
    Phip595CM11 = new GH1("Phip_595MeVCM11", "#phi_{p} Distribution for E_{#gamma} 595 #pm 5MeV (Cos#theta_{CM}0-(-0.1))", 10, -180, 180);
    Phip605CM11 = new GH1("Phip_605MeVCM11", "#phi_{p} Distribution for E_{#gamma} 605 #pm 5MeV (Cos#theta_{CM}0-(-0.1))", 10, -180, 180);
    Phip615CM11 = new GH1("Phip_615MeVCM11", "#phi_{p} Distribution for E_{#gamma} 615 #pm 5MeV (Cos#theta_{CM}0-(-0.1))", 10, -180, 180);

    Phip415CM12 = new GH1("Phip_415MeVCM12", "#phi_{p} Distribution for E_{#gamma} 415 #pm 5MeV (Cos#theta_{CM}-0.1-(-0.2))", 10, -180, 180);
    Phip425CM12 = new GH1("Phip_425MeVCM12", "#phi_{p} Distribution for E_{#gamma} 425 #pm 5MeV (Cos#theta_{CM}-0.1-(-0.2))", 10, -180, 180);
    Phip435CM12 = new GH1("Phip_435MeVCM12", "#phi_{p} Distribution for E_{#gamma} 435 #pm 5MeV (Cos#theta_{CM}-0.1-(-0.2))", 10, -180, 180);
    Phip445CM12 = new GH1("Phip_445MeVCM12", "#phi_{p} Distribution for E_{#gamma} 445 #pm 5MeV (Cos#theta_{CM}-0.1-(-0.2))", 10, -180, 180);
    Phip455CM12 = new GH1("Phip_455MeVCM12", "#phi_{p} Distribution for E_{#gamma} 455 #pm 5MeV (Cos#theta_{CM}-0.1-(-0.2))", 10, -180, 180);
    Phip465CM12 = new GH1("Phip_465MeVCM12", "#phi_{p} Distribution for E_{#gamma} 465 #pm 5MeV (Cos#theta_{CM}-0.1-(-0.2))", 10, -180, 180);
    Phip475CM12 = new GH1("Phip_475MeVCM12", "#phi_{p} Distribution for E_{#gamma} 475 #pm 5MeV (Cos#theta_{CM}-0.1-(-0.2))", 10, -180, 180);
    Phip485CM12 = new GH1("Phip_485MeVCM12", "#phi_{p} Distribution for E_{#gamma} 485 #pm 5MeV (Cos#theta_{CM}-0.1-(-0.2))", 10, -180, 180);
    Phip495CM12 = new GH1("Phip_495MeVCM12", "#phi_{p} Distribution for E_{#gamma} 495 #pm 5MeV (Cos#theta_{CM}-0.1-(-0.2))", 10, -180, 180);
    Phip505CM12 = new GH1("Phip_505MeVCM12", "#phi_{p} Distribution for E_{#gamma} 505 #pm 5MeV (Cos#theta_{CM}-0.1-(-0.2))", 10, -180, 180);
    Phip515CM12 = new GH1("Phip_515MeVCM12", "#phi_{p} Distribution for E_{#gamma} 515 #pm 5MeV (Cos#theta_{CM}-0.1-(-0.2))", 10, -180, 180);
    Phip525CM12 = new GH1("Phip_525MeVCM12", "#phi_{p} Distribution for E_{#gamma} 525 #pm 5MeV (Cos#theta_{CM}-0.1-(-0.2))", 10, -180, 180);
    Phip535CM12 = new GH1("Phip_535MeVCM12", "#phi_{p} Distribution for E_{#gamma} 535 #pm 5MeV (Cos#theta_{CM}-0.1-(-0.2))", 10, -180, 180);
    Phip545CM12 = new GH1("Phip_545MeVCM12", "#phi_{p} Distribution for E_{#gamma} 545 #pm 5MeV (Cos#theta_{CM}-0.1-(-0.2))", 10, -180, 180);
    Phip555CM12 = new GH1("Phip_555MeVCM12", "#phi_{p} Distribution for E_{#gamma} 555 #pm 5MeV (Cos#theta_{CM}-0.1-(-0.2))", 10, -180, 180);
    Phip565CM12 = new GH1("Phip_565MeVCM12", "#phi_{p} Distribution for E_{#gamma} 565 #pm 5MeV (Cos#theta_{CM}-0.1-(-0.2))", 10, -180, 180);
    Phip575CM12 = new GH1("Phip_575MeVCM12", "#phi_{p} Distribution for E_{#gamma} 575 #pm 5MeV (Cos#theta_{CM}-0.1-(-0.2))", 10, -180, 180);
    Phip585CM12 = new GH1("Phip_585MeVCM12", "#phi_{p} Distribution for E_{#gamma} 585 #pm 5MeV (Cos#theta_{CM}-0.1-(-0.2))", 10, -180, 180);
    Phip595CM12 = new GH1("Phip_595MeVCM12", "#phi_{p} Distribution for E_{#gamma} 595 #pm 5MeV (Cos#theta_{CM}-0.1-(-0.2))", 10, -180, 180);
    Phip605CM12 = new GH1("Phip_605MeVCM12", "#phi_{p} Distribution for E_{#gamma} 605 #pm 5MeV (Cos#theta_{CM}-0.1-(-0.2))", 10, -180, 180);
    Phip615CM12 = new GH1("Phip_615MeVCM12", "#phi_{p} Distribution for E_{#gamma} 615 #pm 5MeV (Cos#theta_{CM}-0.1-(-0.2))", 10, -180, 180);

    Phip415CM13 = new GH1("Phip_415MeVCM13", "#phi_{p} Distribution for E_{#gamma} 415 #pm 5MeV (Cos#theta_{CM}-0.2-(-0.3))", 10, -180, 180);
    Phip425CM13 = new GH1("Phip_425MeVCM13", "#phi_{p} Distribution for E_{#gamma} 425 #pm 5MeV (Cos#theta_{CM}-0.2-(-0.3))", 10, -180, 180);
    Phip435CM13 = new GH1("Phip_435MeVCM13", "#phi_{p} Distribution for E_{#gamma} 435 #pm 5MeV (Cos#theta_{CM}-0.2-(-0.3))", 10, -180, 180);
    Phip445CM13 = new GH1("Phip_445MeVCM13", "#phi_{p} Distribution for E_{#gamma} 445 #pm 5MeV (Cos#theta_{CM}-0.2-(-0.3))", 10, -180, 180);
    Phip455CM13 = new GH1("Phip_455MeVCM13", "#phi_{p} Distribution for E_{#gamma} 455 #pm 5MeV (Cos#theta_{CM}-0.2-(-0.3))", 10, -180, 180);
    Phip465CM13 = new GH1("Phip_465MeVCM13", "#phi_{p} Distribution for E_{#gamma} 465 #pm 5MeV (Cos#theta_{CM}-0.2-(-0.3))", 10, -180, 180);
    Phip475CM13 = new GH1("Phip_475MeVCM13", "#phi_{p} Distribution for E_{#gamma} 475 #pm 5MeV (Cos#theta_{CM}-0.2-(-0.3))", 10, -180, 180);
    Phip485CM13 = new GH1("Phip_485MeVCM13", "#phi_{p} Distribution for E_{#gamma} 485 #pm 5MeV (Cos#theta_{CM}-0.2-(-0.3))", 10, -180, 180);
    Phip495CM13 = new GH1("Phip_495MeVCM13", "#phi_{p} Distribution for E_{#gamma} 495 #pm 5MeV (Cos#theta_{CM}-0.2-(-0.3))", 10, -180, 180);
    Phip505CM13 = new GH1("Phip_505MeVCM13", "#phi_{p} Distribution for E_{#gamma} 505 #pm 5MeV (Cos#theta_{CM}-0.2-(-0.3))", 10, -180, 180);
    Phip515CM13 = new GH1("Phip_515MeVCM13", "#phi_{p} Distribution for E_{#gamma} 515 #pm 5MeV (Cos#theta_{CM}-0.2-(-0.3))", 10, -180, 180);
    Phip525CM13 = new GH1("Phip_525MeVCM13", "#phi_{p} Distribution for E_{#gamma} 525 #pm 5MeV (Cos#theta_{CM}-0.2-(-0.3))", 10, -180, 180);
    Phip535CM13 = new GH1("Phip_535MeVCM13", "#phi_{p} Distribution for E_{#gamma} 535 #pm 5MeV (Cos#theta_{CM}-0.2-(-0.3))", 10, -180, 180);
    Phip545CM13 = new GH1("Phip_545MeVCM13", "#phi_{p} Distribution for E_{#gamma} 545 #pm 5MeV (Cos#theta_{CM}-0.2-(-0.3))", 10, -180, 180);
    Phip555CM13 = new GH1("Phip_555MeVCM13", "#phi_{p} Distribution for E_{#gamma} 555 #pm 5MeV (Cos#theta_{CM}-0.2-(-0.3))", 10, -180, 180);
    Phip565CM13 = new GH1("Phip_565MeVCM13", "#phi_{p} Distribution for E_{#gamma} 565 #pm 5MeV (Cos#theta_{CM}-0.2-(-0.3))", 10, -180, 180);
    Phip575CM13 = new GH1("Phip_575MeVCM13", "#phi_{p} Distribution for E_{#gamma} 575 #pm 5MeV (Cos#theta_{CM}-0.2-(-0.3))", 10, -180, 180);
    Phip585CM13 = new GH1("Phip_585MeVCM13", "#phi_{p} Distribution for E_{#gamma} 585 #pm 5MeV (Cos#theta_{CM}-0.2-(-0.3))", 10, -180, 180);
    Phip595CM13 = new GH1("Phip_595MeVCM13", "#phi_{p} Distribution for E_{#gamma} 595 #pm 5MeV (Cos#theta_{CM}-0.2-(-0.3))", 10, -180, 180);
    Phip605CM13 = new GH1("Phip_605MeVCM13", "#phi_{p} Distribution for E_{#gamma} 605 #pm 5MeV (Cos#theta_{CM}-0.2-(-0.3))", 10, -180, 180);
    Phip615CM13 = new GH1("Phip_615MeVCM13", "#phi_{p} Distribution for E_{#gamma} 615 #pm 5MeV (Cos#theta_{CM}-0.2-(-0.3))", 10, -180, 180);

    Phip415CM14 = new GH1("Phip_415MeVCM14", "#phi_{p} Distribution for E_{#gamma} 415 #pm 5MeV (Cos#theta_{CM}-0.3-(-0.4))", 10, -180, 180);
    Phip425CM14 = new GH1("Phip_425MeVCM14", "#phi_{p} Distribution for E_{#gamma} 425 #pm 5MeV (Cos#theta_{CM}-0.3-(-0.4))", 10, -180, 180);
    Phip435CM14 = new GH1("Phip_435MeVCM14", "#phi_{p} Distribution for E_{#gamma} 435 #pm 5MeV (Cos#theta_{CM}-0.3-(-0.4))", 10, -180, 180);
    Phip445CM14 = new GH1("Phip_445MeVCM14", "#phi_{p} Distribution for E_{#gamma} 445 #pm 5MeV (Cos#theta_{CM}-0.3-(-0.4))", 10, -180, 180);
    Phip455CM14 = new GH1("Phip_455MeVCM14", "#phi_{p} Distribution for E_{#gamma} 455 #pm 5MeV (Cos#theta_{CM}-0.3-(-0.4))", 10, -180, 180);
    Phip465CM14 = new GH1("Phip_465MeVCM14", "#phi_{p} Distribution for E_{#gamma} 465 #pm 5MeV (Cos#theta_{CM}-0.3-(-0.4))", 10, -180, 180);
    Phip475CM14 = new GH1("Phip_475MeVCM14", "#phi_{p} Distribution for E_{#gamma} 475 #pm 5MeV (Cos#theta_{CM}-0.3-(-0.4))", 10, -180, 180);
    Phip485CM14 = new GH1("Phip_485MeVCM14", "#phi_{p} Distribution for E_{#gamma} 485 #pm 5MeV (Cos#theta_{CM}-0.3-(-0.4))", 10, -180, 180);
    Phip495CM14 = new GH1("Phip_495MeVCM14", "#phi_{p} Distribution for E_{#gamma} 495 #pm 5MeV (Cos#theta_{CM}-0.3-(-0.4))", 10, -180, 180);
    Phip505CM14 = new GH1("Phip_505MeVCM14", "#phi_{p} Distribution for E_{#gamma} 505 #pm 5MeV (Cos#theta_{CM}-0.3-(-0.4))", 10, -180, 180);
    Phip515CM14 = new GH1("Phip_515MeVCM14", "#phi_{p} Distribution for E_{#gamma} 515 #pm 5MeV (Cos#theta_{CM}-0.3-(-0.4))", 10, -180, 180);
    Phip525CM14 = new GH1("Phip_525MeVCM14", "#phi_{p} Distribution for E_{#gamma} 525 #pm 5MeV (Cos#theta_{CM}-0.3-(-0.4))", 10, -180, 180);
    Phip535CM14 = new GH1("Phip_535MeVCM14", "#phi_{p} Distribution for E_{#gamma} 535 #pm 5MeV (Cos#theta_{CM}-0.3-(-0.4))", 10, -180, 180);
    Phip545CM14 = new GH1("Phip_545MeVCM14", "#phi_{p} Distribution for E_{#gamma} 545 #pm 5MeV (Cos#theta_{CM}-0.3-(-0.4))", 10, -180, 180);
    Phip555CM14 = new GH1("Phip_555MeVCM14", "#phi_{p} Distribution for E_{#gamma} 555 #pm 5MeV (Cos#theta_{CM}-0.3-(-0.4))", 10, -180, 180);
    Phip565CM14 = new GH1("Phip_565MeVCM14", "#phi_{p} Distribution for E_{#gamma} 565 #pm 5MeV (Cos#theta_{CM}-0.3-(-0.4))", 10, -180, 180);
    Phip575CM14 = new GH1("Phip_575MeVCM14", "#phi_{p} Distribution for E_{#gamma} 575 #pm 5MeV (Cos#theta_{CM}-0.3-(-0.4))", 10, -180, 180);
    Phip585CM14 = new GH1("Phip_585MeVCM14", "#phi_{p} Distribution for E_{#gamma} 585 #pm 5MeV (Cos#theta_{CM}-0.3-(-0.4))", 10, -180, 180);
    Phip595CM14 = new GH1("Phip_595MeVCM14", "#phi_{p} Distribution for E_{#gamma} 595 #pm 5MeV (Cos#theta_{CM}-0.3-(-0.4))", 10, -180, 180);
    Phip605CM14 = new GH1("Phip_605MeVCM14", "#phi_{p} Distribution for E_{#gamma} 605 #pm 5MeV (Cos#theta_{CM}-0.3-(-0.4))", 10, -180, 180);
    Phip615CM14 = new GH1("Phip_615MeVCM14", "#phi_{p} Distribution for E_{#gamma} 615 #pm 5MeV (Cos#theta_{CM}-0.3-(-0.4))", 10, -180, 180);

    Phip415CM15 = new GH1("Phip_415MeVCM15", "#phi_{p} Distribution for E_{#gamma} 415 #pm 5MeV (Cos#theta_{CM}-0.4-(-0.5))", 10, -180, 180);
    Phip425CM15 = new GH1("Phip_425MeVCM15", "#phi_{p} Distribution for E_{#gamma} 425 #pm 5MeV (Cos#theta_{CM}-0.4-(-0.5))", 10, -180, 180);
    Phip435CM15 = new GH1("Phip_435MeVCM15", "#phi_{p} Distribution for E_{#gamma} 435 #pm 5MeV (Cos#theta_{CM}-0.4-(-0.5))", 10, -180, 180);
    Phip445CM15 = new GH1("Phip_445MeVCM15", "#phi_{p} Distribution for E_{#gamma} 445 #pm 5MeV (Cos#theta_{CM}-0.4-(-0.5))", 10, -180, 180);
    Phip455CM15 = new GH1("Phip_455MeVCM15", "#phi_{p} Distribution for E_{#gamma} 455 #pm 5MeV (Cos#theta_{CM}-0.4-(-0.5))", 10, -180, 180);
    Phip465CM15 = new GH1("Phip_465MeVCM15", "#phi_{p} Distribution for E_{#gamma} 465 #pm 5MeV (Cos#theta_{CM}-0.4-(-0.5))", 10, -180, 180);
    Phip475CM15 = new GH1("Phip_475MeVCM15", "#phi_{p} Distribution for E_{#gamma} 475 #pm 5MeV (Cos#theta_{CM}-0.4-(-0.5))", 10, -180, 180);
    Phip485CM15 = new GH1("Phip_485MeVCM15", "#phi_{p} Distribution for E_{#gamma} 485 #pm 5MeV (Cos#theta_{CM}-0.4-(-0.5))", 10, -180, 180);
    Phip495CM15 = new GH1("Phip_495MeVCM15", "#phi_{p} Distribution for E_{#gamma} 495 #pm 5MeV (Cos#theta_{CM}-0.4-(-0.5))", 10, -180, 180);
    Phip505CM15 = new GH1("Phip_505MeVCM15", "#phi_{p} Distribution for E_{#gamma} 505 #pm 5MeV (Cos#theta_{CM}-0.4-(-0.5))", 10, -180, 180);
    Phip515CM15 = new GH1("Phip_515MeVCM15", "#phi_{p} Distribution for E_{#gamma} 515 #pm 5MeV (Cos#theta_{CM}-0.4-(-0.5))", 10, -180, 180);
    Phip525CM15 = new GH1("Phip_525MeVCM15", "#phi_{p} Distribution for E_{#gamma} 525 #pm 5MeV (Cos#theta_{CM}-0.4-(-0.5))", 10, -180, 180);
    Phip535CM15 = new GH1("Phip_535MeVCM15", "#phi_{p} Distribution for E_{#gamma} 535 #pm 5MeV (Cos#theta_{CM}-0.4-(-0.5))", 10, -180, 180);
    Phip545CM15 = new GH1("Phip_545MeVCM15", "#phi_{p} Distribution for E_{#gamma} 545 #pm 5MeV (Cos#theta_{CM}-0.4-(-0.5))", 10, -180, 180);
    Phip555CM15 = new GH1("Phip_555MeVCM15", "#phi_{p} Distribution for E_{#gamma} 555 #pm 5MeV (Cos#theta_{CM}-0.4-(-0.5))", 10, -180, 180);
    Phip565CM15 = new GH1("Phip_565MeVCM15", "#phi_{p} Distribution for E_{#gamma} 565 #pm 5MeV (Cos#theta_{CM}-0.4-(-0.5))", 10, -180, 180);
    Phip575CM15 = new GH1("Phip_575MeVCM15", "#phi_{p} Distribution for E_{#gamma} 575 #pm 5MeV (Cos#theta_{CM}-0.4-(-0.5))", 10, -180, 180);
    Phip585CM15 = new GH1("Phip_585MeVCM15", "#phi_{p} Distribution for E_{#gamma} 585 #pm 5MeV (Cos#theta_{CM}-0.4-(-0.5))", 10, -180, 180);
    Phip595CM15 = new GH1("Phip_595MeVCM15", "#phi_{p} Distribution for E_{#gamma} 595 #pm 5MeV (Cos#theta_{CM}-0.4-(-0.5))", 10, -180, 180);
    Phip605CM15 = new GH1("Phip_605MeVCM15", "#phi_{p} Distribution for E_{#gamma} 605 #pm 5MeV (Cos#theta_{CM}-0.4-(-0.5))", 10, -180, 180);
    Phip615CM15 = new GH1("Phip_615MeVCM15", "#phi_{p} Distribution for E_{#gamma} 615 #pm 5MeV (Cos#theta_{CM}-0.4-(-0.5))", 10, -180, 180);

    Phip415CM16 = new GH1("Phip_415MeVCM16", "#phi_{p} Distribution for E_{#gamma} 415 #pm 5MeV (Cos#theta_{CM}-0.5-(-0.6))", 10, -180, 180);
    Phip425CM16 = new GH1("Phip_425MeVCM16", "#phi_{p} Distribution for E_{#gamma} 425 #pm 5MeV (Cos#theta_{CM}-0.5-(-0.6))", 10, -180, 180);
    Phip435CM16 = new GH1("Phip_435MeVCM16", "#phi_{p} Distribution for E_{#gamma} 435 #pm 5MeV (Cos#theta_{CM}-0.5-(-0.6))", 10, -180, 180);
    Phip445CM16 = new GH1("Phip_445MeVCM16", "#phi_{p} Distribution for E_{#gamma} 445 #pm 5MeV (Cos#theta_{CM}-0.5-(-0.6))", 10, -180, 180);
    Phip455CM16 = new GH1("Phip_455MeVCM16", "#phi_{p} Distribution for E_{#gamma} 455 #pm 5MeV (Cos#theta_{CM}-0.5-(-0.6))", 10, -180, 180);
    Phip465CM16 = new GH1("Phip_465MeVCM16", "#phi_{p} Distribution for E_{#gamma} 465 #pm 5MeV (Cos#theta_{CM}-0.5-(-0.6))", 10, -180, 180);
    Phip475CM16 = new GH1("Phip_475MeVCM16", "#phi_{p} Distribution for E_{#gamma} 475 #pm 5MeV (Cos#theta_{CM}-0.5-(-0.6))", 10, -180, 180);
    Phip485CM16 = new GH1("Phip_485MeVCM16", "#phi_{p} Distribution for E_{#gamma} 485 #pm 5MeV (Cos#theta_{CM}-0.5-(-0.6))", 10, -180, 180);
    Phip495CM16 = new GH1("Phip_495MeVCM16", "#phi_{p} Distribution for E_{#gamma} 495 #pm 5MeV (Cos#theta_{CM}-0.5-(-0.6))", 10, -180, 180);
    Phip505CM16 = new GH1("Phip_505MeVCM16", "#phi_{p} Distribution for E_{#gamma} 505 #pm 5MeV (Cos#theta_{CM}-0.5-(-0.6))", 10, -180, 180);
    Phip515CM16 = new GH1("Phip_515MeVCM16", "#phi_{p} Distribution for E_{#gamma} 515 #pm 5MeV (Cos#theta_{CM}-0.5-(-0.6))", 10, -180, 180);
    Phip525CM16 = new GH1("Phip_525MeVCM16", "#phi_{p} Distribution for E_{#gamma} 525 #pm 5MeV (Cos#theta_{CM}-0.5-(-0.6))", 10, -180, 180);
    Phip535CM16 = new GH1("Phip_535MeVCM16", "#phi_{p} Distribution for E_{#gamma} 535 #pm 5MeV (Cos#theta_{CM}-0.5-(-0.6))", 10, -180, 180);
    Phip545CM16 = new GH1("Phip_545MeVCM16", "#phi_{p} Distribution for E_{#gamma} 545 #pm 5MeV (Cos#theta_{CM}-0.5-(-0.6))", 10, -180, 180);
    Phip555CM16 = new GH1("Phip_555MeVCM16", "#phi_{p} Distribution for E_{#gamma} 555 #pm 5MeV (Cos#theta_{CM}-0.5-(-0.6))", 10, -180, 180);
    Phip565CM16 = new GH1("Phip_565MeVCM16", "#phi_{p} Distribution for E_{#gamma} 565 #pm 5MeV (Cos#theta_{CM}-0.5-(-0.6))", 10, -180, 180);
    Phip575CM16 = new GH1("Phip_575MeVCM16", "#phi_{p} Distribution for E_{#gamma} 575 #pm 5MeV (Cos#theta_{CM}-0.5-(-0.6))", 10, -180, 180);
    Phip585CM16 = new GH1("Phip_585MeVCM16", "#phi_{p} Distribution for E_{#gamma} 585 #pm 5MeV (Cos#theta_{CM}-0.5-(-0.6))", 10, -180, 180);
    Phip595CM16 = new GH1("Phip_595MeVCM16", "#phi_{p} Distribution for E_{#gamma} 595 #pm 5MeV (Cos#theta_{CM}-0.5-(-0.6))", 10, -180, 180);
    Phip605CM16 = new GH1("Phip_605MeVCM16", "#phi_{p} Distribution for E_{#gamma} 605 #pm 5MeV (Cos#theta_{CM}-0.5-(-0.6))", 10, -180, 180);
    Phip615CM16 = new GH1("Phip_615MeVCM16", "#phi_{p} Distribution for E_{#gamma} 615 #pm 5MeV (Cos#theta_{CM}-0.5-(-0.6))", 10, -180, 180);

    Phip415CM17 = new GH1("Phip_415MeVCM17", "#phi_{p} Distribution for E_{#gamma} 415 #pm 5MeV (Cos#theta_{CM}-0.6-(-0.7))", 10, -180, 180);
    Phip425CM17 = new GH1("Phip_425MeVCM17", "#phi_{p} Distribution for E_{#gamma} 425 #pm 5MeV (Cos#theta_{CM}-0.6-(-0.7))", 10, -180, 180);
    Phip435CM17 = new GH1("Phip_435MeVCM17", "#phi_{p} Distribution for E_{#gamma} 435 #pm 5MeV (Cos#theta_{CM}-0.6-(-0.7))", 10, -180, 180);
    Phip445CM17 = new GH1("Phip_445MeVCM17", "#phi_{p} Distribution for E_{#gamma} 445 #pm 5MeV (Cos#theta_{CM}-0.6-(-0.7))", 10, -180, 180);
    Phip455CM17 = new GH1("Phip_455MeVCM17", "#phi_{p} Distribution for E_{#gamma} 455 #pm 5MeV (Cos#theta_{CM}-0.6-(-0.7))", 10, -180, 180);
    Phip465CM17 = new GH1("Phip_465MeVCM17", "#phi_{p} Distribution for E_{#gamma} 465 #pm 5MeV (Cos#theta_{CM}-0.6-(-0.7))", 10, -180, 180);
    Phip475CM17 = new GH1("Phip_475MeVCM17", "#phi_{p} Distribution for E_{#gamma} 475 #pm 5MeV (Cos#theta_{CM}-0.6-(-0.7))", 10, -180, 180);
    Phip485CM17 = new GH1("Phip_485MeVCM17", "#phi_{p} Distribution for E_{#gamma} 485 #pm 5MeV (Cos#theta_{CM}-0.6-(-0.7))", 10, -180, 180);
    Phip495CM17 = new GH1("Phip_495MeVCM17", "#phi_{p} Distribution for E_{#gamma} 495 #pm 5MeV (Cos#theta_{CM}-0.6-(-0.7))", 10, -180, 180);
    Phip505CM17 = new GH1("Phip_505MeVCM17", "#phi_{p} Distribution for E_{#gamma} 505 #pm 5MeV (Cos#theta_{CM}-0.6-(-0.7))", 10, -180, 180);
    Phip515CM17 = new GH1("Phip_515MeVCM17", "#phi_{p} Distribution for E_{#gamma} 515 #pm 5MeV (Cos#theta_{CM}-0.6-(-0.7))", 10, -180, 180);
    Phip525CM17 = new GH1("Phip_525MeVCM17", "#phi_{p} Distribution for E_{#gamma} 525 #pm 5MeV (Cos#theta_{CM}-0.6-(-0.7))", 10, -180, 180);
    Phip535CM17 = new GH1("Phip_535MeVCM17", "#phi_{p} Distribution for E_{#gamma} 535 #pm 5MeV (Cos#theta_{CM}-0.6-(-0.7))", 10, -180, 180);
    Phip545CM17 = new GH1("Phip_545MeVCM17", "#phi_{p} Distribution for E_{#gamma} 545 #pm 5MeV (Cos#theta_{CM}-0.6-(-0.7))", 10, -180, 180);
    Phip555CM17 = new GH1("Phip_555MeVCM17", "#phi_{p} Distribution for E_{#gamma} 555 #pm 5MeV (Cos#theta_{CM}-0.6-(-0.7))", 10, -180, 180);
    Phip565CM17 = new GH1("Phip_565MeVCM17", "#phi_{p} Distribution for E_{#gamma} 565 #pm 5MeV (Cos#theta_{CM}-0.6-(-0.7))", 10, -180, 180);
    Phip575CM17 = new GH1("Phip_575MeVCM17", "#phi_{p} Distribution for E_{#gamma} 575 #pm 5MeV (Cos#theta_{CM}-0.6-(-0.7))", 10, -180, 180);
    Phip585CM17 = new GH1("Phip_585MeVCM17", "#phi_{p} Distribution for E_{#gamma} 585 #pm 5MeV (Cos#theta_{CM}-0.6-(-0.7))", 10, -180, 180);
    Phip595CM17 = new GH1("Phip_595MeVCM17", "#phi_{p} Distribution for E_{#gamma} 595 #pm 5MeV (Cos#theta_{CM}-0.6-(-0.7))", 10, -180, 180);
    Phip605CM17 = new GH1("Phip_605MeVCM17", "#phi_{p} Distribution for E_{#gamma} 605 #pm 5MeV (Cos#theta_{CM}-0.6-(-0.7))", 10, -180, 180);
    Phip615CM17 = new GH1("Phip_615MeVCM17", "#phi_{p} Distribution for E_{#gamma} 615 #pm 5MeV (Cos#theta_{CM}-0.6-(-0.7))", 10, -180, 180);

    Phip415CM18 = new GH1("Phip_415MeVCM18", "#phi_{p} Distribution for E_{#gamma} 415 #pm 5MeV (Cos#theta_{CM}-0.7-(-0.8))", 10, -180, 180);
    Phip425CM18 = new GH1("Phip_425MeVCM18", "#phi_{p} Distribution for E_{#gamma} 425 #pm 5MeV (Cos#theta_{CM}-0.7-(-0.8))", 10, -180, 180);
    Phip435CM18 = new GH1("Phip_435MeVCM18", "#phi_{p} Distribution for E_{#gamma} 435 #pm 5MeV (Cos#theta_{CM}-0.7-(-0.8))", 10, -180, 180);
    Phip445CM18 = new GH1("Phip_445MeVCM18", "#phi_{p} Distribution for E_{#gamma} 445 #pm 5MeV (Cos#theta_{CM}-0.7-(-0.8))", 10, -180, 180);
    Phip455CM18 = new GH1("Phip_455MeVCM18", "#phi_{p} Distribution for E_{#gamma} 455 #pm 5MeV (Cos#theta_{CM}-0.7-(-0.8))", 10, -180, 180);
    Phip465CM18 = new GH1("Phip_465MeVCM18", "#phi_{p} Distribution for E_{#gamma} 465 #pm 5MeV (Cos#theta_{CM}-0.7-(-0.8))", 10, -180, 180);
    Phip475CM18 = new GH1("Phip_475MeVCM18", "#phi_{p} Distribution for E_{#gamma} 475 #pm 5MeV (Cos#theta_{CM}-0.7-(-0.8))", 10, -180, 180);
    Phip485CM18 = new GH1("Phip_485MeVCM18", "#phi_{p} Distribution for E_{#gamma} 485 #pm 5MeV (Cos#theta_{CM}-0.7-(-0.8))", 10, -180, 180);
    Phip495CM18 = new GH1("Phip_495MeVCM18", "#phi_{p} Distribution for E_{#gamma} 495 #pm 5MeV (Cos#theta_{CM}-0.7-(-0.8))", 10, -180, 180);
    Phip505CM18 = new GH1("Phip_505MeVCM18", "#phi_{p} Distribution for E_{#gamma} 505 #pm 5MeV (Cos#theta_{CM}-0.7-(-0.8))", 10, -180, 180);
    Phip515CM18 = new GH1("Phip_515MeVCM18", "#phi_{p} Distribution for E_{#gamma} 515 #pm 5MeV (Cos#theta_{CM}-0.7-(-0.8))", 10, -180, 180);
    Phip525CM18 = new GH1("Phip_525MeVCM18", "#phi_{p} Distribution for E_{#gamma} 525 #pm 5MeV (Cos#theta_{CM}-0.7-(-0.8))", 10, -180, 180);
    Phip535CM18 = new GH1("Phip_535MeVCM18", "#phi_{p} Distribution for E_{#gamma} 535 #pm 5MeV (Cos#theta_{CM}-0.7-(-0.8))", 10, -180, 180);
    Phip545CM18 = new GH1("Phip_545MeVCM18", "#phi_{p} Distribution for E_{#gamma} 545 #pm 5MeV (Cos#theta_{CM}-0.7-(-0.8))", 10, -180, 180);
    Phip555CM18 = new GH1("Phip_555MeVCM18", "#phi_{p} Distribution for E_{#gamma} 555 #pm 5MeV (Cos#theta_{CM}-0.7-(-0.8))", 10, -180, 180);
    Phip565CM18 = new GH1("Phip_565MeVCM18", "#phi_{p} Distribution for E_{#gamma} 565 #pm 5MeV (Cos#theta_{CM}-0.7-(-0.8))", 10, -180, 180);
    Phip575CM18 = new GH1("Phip_575MeVCM18", "#phi_{p} Distribution for E_{#gamma} 575 #pm 5MeV (Cos#theta_{CM}-0.7-(-0.8))", 10, -180, 180);
    Phip585CM18 = new GH1("Phip_585MeVCM18", "#phi_{p} Distribution for E_{#gamma} 585 #pm 5MeV (Cos#theta_{CM}-0.7-(-0.8))", 10, -180, 180);
    Phip595CM18 = new GH1("Phip_595MeVCM18", "#phi_{p} Distribution for E_{#gamma} 595 #pm 5MeV (Cos#theta_{CM}-0.7-(-0.8))", 10, -180, 180);
    Phip605CM18 = new GH1("Phip_605MeVCM18", "#phi_{p} Distribution for E_{#gamma} 605 #pm 5MeV (Cos#theta_{CM}-0.7-(-0.8))", 10, -180, 180);
    Phip615CM18 = new GH1("Phip_615MeVCM18", "#phi_{p} Distribution for E_{#gamma} 615 #pm 5MeV (Cos#theta_{CM}-0.7-(-0.8))", 10, -180, 180);

    Phip415CM19 = new GH1("Phip_415MeVCM19", "#phi_{p} Distribution for E_{#gamma} 415 #pm 5MeV (Cos#theta_{CM}-0.8-(-0.9))", 10, -180, 180);
    Phip425CM19 = new GH1("Phip_425MeVCM19", "#phi_{p} Distribution for E_{#gamma} 425 #pm 5MeV (Cos#theta_{CM}-0.8-(-0.9))", 10, -180, 180);
    Phip435CM19 = new GH1("Phip_435MeVCM19", "#phi_{p} Distribution for E_{#gamma} 435 #pm 5MeV (Cos#theta_{CM}-0.8-(-0.9))", 10, -180, 180);
    Phip445CM19 = new GH1("Phip_445MeVCM19", "#phi_{p} Distribution for E_{#gamma} 445 #pm 5MeV (Cos#theta_{CM}-0.8-(-0.9))", 10, -180, 180);
    Phip455CM19 = new GH1("Phip_455MeVCM19", "#phi_{p} Distribution for E_{#gamma} 455 #pm 5MeV (Cos#theta_{CM}-0.8-(-0.9))", 10, -180, 180);
    Phip465CM19 = new GH1("Phip_465MeVCM19", "#phi_{p} Distribution for E_{#gamma} 465 #pm 5MeV (Cos#theta_{CM}-0.8-(-0.9))", 10, -180, 180);
    Phip475CM19 = new GH1("Phip_475MeVCM19", "#phi_{p} Distribution for E_{#gamma} 475 #pm 5MeV (Cos#theta_{CM}-0.8-(-0.9))", 10, -180, 180);
    Phip485CM19 = new GH1("Phip_485MeVCM19", "#phi_{p} Distribution for E_{#gamma} 485 #pm 5MeV (Cos#theta_{CM}-0.8-(-0.9))", 10, -180, 180);
    Phip495CM19 = new GH1("Phip_495MeVCM19", "#phi_{p} Distribution for E_{#gamma} 495 #pm 5MeV (Cos#theta_{CM}-0.8-(-0.9))", 10, -180, 180);
    Phip505CM19 = new GH1("Phip_505MeVCM19", "#phi_{p} Distribution for E_{#gamma} 505 #pm 5MeV (Cos#theta_{CM}-0.8-(-0.9))", 10, -180, 180);
    Phip515CM19 = new GH1("Phip_515MeVCM19", "#phi_{p} Distribution for E_{#gamma} 515 #pm 5MeV (Cos#theta_{CM}-0.8-(-0.9))", 10, -180, 180);
    Phip525CM19 = new GH1("Phip_525MeVCM19", "#phi_{p} Distribution for E_{#gamma} 525 #pm 5MeV (Cos#theta_{CM}-0.8-(-0.9))", 10, -180, 180);
    Phip535CM19 = new GH1("Phip_535MeVCM19", "#phi_{p} Distribution for E_{#gamma} 535 #pm 5MeV (Cos#theta_{CM}-0.8-(-0.9))", 10, -180, 180);
    Phip545CM19 = new GH1("Phip_545MeVCM19", "#phi_{p} Distribution for E_{#gamma} 545 #pm 5MeV (Cos#theta_{CM}-0.8-(-0.9))", 10, -180, 180);
    Phip555CM19 = new GH1("Phip_555MeVCM19", "#phi_{p} Distribution for E_{#gamma} 555 #pm 5MeV (Cos#theta_{CM}-0.8-(-0.9))", 10, -180, 180);
    Phip565CM19 = new GH1("Phip_565MeVCM19", "#phi_{p} Distribution for E_{#gamma} 565 #pm 5MeV (Cos#theta_{CM}-0.8-(-0.9))", 10, -180, 180);
    Phip575CM19 = new GH1("Phip_575MeVCM19", "#phi_{p} Distribution for E_{#gamma} 575 #pm 5MeV (Cos#theta_{CM}-0.8-(-0.9))", 10, -180, 180);
    Phip585CM19 = new GH1("Phip_585MeVCM19", "#phi_{p} Distribution for E_{#gamma} 585 #pm 5MeV (Cos#theta_{CM}-0.8-(-0.9))", 10, -180, 180);
    Phip595CM19 = new GH1("Phip_595MeVCM19", "#phi_{p} Distribution for E_{#gamma} 595 #pm 5MeV (Cos#theta_{CM}-0.8-(-0.9))", 10, -180, 180);
    Phip605CM19 = new GH1("Phip_605MeVCM19", "#phi_{p} Distribution for E_{#gamma} 605 #pm 5MeV (Cos#theta_{CM}-0.8-(-0.9))", 10, -180, 180);
    Phip615CM19 = new GH1("Phip_615MeVCM19", "#phi_{p} Distribution for E_{#gamma} 615 #pm 5MeV (Cos#theta_{CM}-0.8-(-0.9))", 10, -180, 180);

    Phip415CM20 = new GH1("Phip_415MeVCM20", "#phi_{p} Distribution for E_{#gamma} 415 #pm 5MeV (Cos#theta_{CM}-0.9-(-1.0))", 10, -180, 180);
    Phip425CM20 = new GH1("Phip_425MeVCM20", "#phi_{p} Distribution for E_{#gamma} 425 #pm 5MeV (Cos#theta_{CM}-0.9-(-1.0))", 10, -180, 180);
    Phip435CM20 = new GH1("Phip_435MeVCM20", "#phi_{p} Distribution for E_{#gamma} 435 #pm 5MeV (Cos#theta_{CM}-0.9-(-1.0))", 10, -180, 180);
    Phip445CM20 = new GH1("Phip_445MeVCM20", "#phi_{p} Distribution for E_{#gamma} 445 #pm 5MeV (Cos#theta_{CM}-0.9-(-1.0))", 10, -180, 180);
    Phip455CM20 = new GH1("Phip_455MeVCM20", "#phi_{p} Distribution for E_{#gamma} 455 #pm 5MeV (Cos#theta_{CM}-0.9-(-1.0))", 10, -180, 180);
    Phip465CM20 = new GH1("Phip_465MeVCM20", "#phi_{p} Distribution for E_{#gamma} 465 #pm 5MeV (Cos#theta_{CM}-0.9-(-1.0))", 10, -180, 180);
    Phip475CM20 = new GH1("Phip_475MeVCM20", "#phi_{p} Distribution for E_{#gamma} 475 #pm 5MeV (Cos#theta_{CM}-0.9-(-1.0))", 10, -180, 180);
    Phip485CM20 = new GH1("Phip_485MeVCM20", "#phi_{p} Distribution for E_{#gamma} 485 #pm 5MeV (Cos#theta_{CM}-0.9-(-1.0))", 10, -180, 180);
    Phip495CM20 = new GH1("Phip_495MeVCM20", "#phi_{p} Distribution for E_{#gamma} 495 #pm 5MeV (Cos#theta_{CM}-0.9-(-1.0))", 10, -180, 180);
    Phip505CM20 = new GH1("Phip_505MeVCM20", "#phi_{p} Distribution for E_{#gamma} 505 #pm 5MeV (Cos#theta_{CM}-0.9-(-1.0))", 10, -180, 180);
    Phip515CM20 = new GH1("Phip_515MeVCM20", "#phi_{p} Distribution for E_{#gamma} 515 #pm 5MeV (Cos#theta_{CM}-0.9-(-1.0))", 10, -180, 180);
    Phip525CM20 = new GH1("Phip_525MeVCM20", "#phi_{p} Distribution for E_{#gamma} 525 #pm 5MeV (Cos#theta_{CM}-0.9-(-1.0))", 10, -180, 180);
    Phip535CM20 = new GH1("Phip_535MeVCM20", "#phi_{p} Distribution for E_{#gamma} 535 #pm 5MeV (Cos#theta_{CM}-0.9-(-1.0))", 10, -180, 180);
    Phip545CM20 = new GH1("Phip_545MeVCM20", "#phi_{p} Distribution for E_{#gamma} 545 #pm 5MeV (Cos#theta_{CM}-0.9-(-1.0))", 10, -180, 180);
    Phip555CM20 = new GH1("Phip_555MeVCM20", "#phi_{p} Distribution for E_{#gamma} 555 #pm 5MeV (Cos#theta_{CM}-0.9-(-1.0))", 10, -180, 180);
    Phip565CM20 = new GH1("Phip_565MeVCM20", "#phi_{p} Distribution for E_{#gamma} 565 #pm 5MeV (Cos#theta_{CM}-0.9-(-1.0))", 10, -180, 180);
    Phip575CM20 = new GH1("Phip_575MeVCM20", "#phi_{p} Distribution for E_{#gamma} 575 #pm 5MeV (Cos#theta_{CM}-0.9-(-1.0))", 10, -180, 180);
    Phip585CM20 = new GH1("Phip_585MeVCM20", "#phi_{p} Distribution for E_{#gamma} 585 #pm 5MeV (Cos#theta_{CM}-0.9-(-1.0))", 10, -180, 180);
    Phip595CM20 = new GH1("Phip_595MeVCM20", "#phi_{p} Distribution for E_{#gamma} 595 #pm 5MeV (Cos#theta_{CM}-0.9-(-1.0))", 10, -180, 180);
    Phip605CM20 = new GH1("Phip_605MeVCM20", "#phi_{p} Distribution for E_{#gamma} 605 #pm 5MeV (Cos#theta_{CM}-0.9-(-1.0))", 10, -180, 180);
    Phip615CM20 = new GH1("Phip_615MeVCM20", "#phi_{p} Distribution for E_{#gamma} 615 #pm 5MeV (Cos#theta_{CM}-0.9-(-1.0))", 10, -180, 180);

    Eg2 = new TH1D( "Eg2", "E_{#gamma} Distribution", 200, 100, 1600 );
    EgPrompt = new TH1D( "EgPrompt", "E_{#gamma} Distribution", 200, 100, 1600 );
    EgRandom = new TH1D( "EgRandom", "E_{#gamma} Distribution", 200, 100, 1600 );

    Phip320 = new TH1D("Phip320", "#phi_{p} (300-340MeV)", 24, -4, 4);
    Phip360 = new TH1D("Phip360", "#phi_{p} (340-380MeV)", 24, -4, 4);
    Phip400 = new TH1D("Phip400", "#phi_{p} (380-420MeV)", 24, -4, 4);
    Phip440 = new TH1D("Phip440", "#phi_{p} (420-460MeV)", 24, -4, 4);
    Phip480 = new TH1D("Phip480", "#phi_{p} (460-500MeV)", 24, -4, 4);
    Phip520 = new TH1D("Phip520", "#phi_{p} (500-540MeV)", 24, -4, 4);
    Phip560 = new TH1D("Phip560", "#phi_{p} (540-580MeV)", 24, -4, 4);
    Phip600 = new TH1D("Phip600", "#phi_{p} (580-620MeV)", 24, -4, 4);
    Phip640 = new TH1D("Phip640", "#phi_{p} (620-660MeV)", 24, -4, 4);
    Phip680 = new TH1D("Phip680", "#phi_{p} (660-700MeV)", 24, -4, 4);
    Phip320Random = new TH1D("Phip320Random", "#phi_{p} (300-340MeV)", 24, -4, 4);
    Phip360Random = new TH1D("Phip360Random", "#phi_{p} (340-380MeV)", 24, -4, 4);
    Phip400Random = new TH1D("Phip400Random", "#phi_{p} (380-420MeV)", 24, -4, 4);
    Phip440Random = new TH1D("Phip440Random", "#phi_{p} (420-460MeV)", 24, -4, 4);
    Phip480Random = new TH1D("Phip480Random", "#phi_{p} (460-500MeV)", 24, -4, 4);
    Phip520Random = new TH1D("Phip520Random", "#phi_{p} (500-540MeV)", 24, -4, 4);
    Phip560Random = new TH1D("Phip560Random", "#phi_{p} (540-580MeV)", 24, -4, 4);
    Phip600Random = new TH1D("Phip600Random", "#phi_{p} (580-620MeV)", 24, -4, 4);
    Phip640Random = new TH1D("Phip640Random", "#phi_{p} (620-660MeV)", 24, -4, 4);
    Phip680Random = new TH1D("Phip680Random", "#phi_{p} (660-700MeV)", 24, -4, 4);
    Phip320Prompt = new TH1D("Phip320Prompt", "#phi_{p} (300-340MeV)", 24, -4, 4);
    Phip360Prompt = new TH1D("Phip360Prompt", "#phi_{p} (340-380MeV)", 24, -4, 4);
    Phip400Prompt = new TH1D("Phip400Prompt", "#phi_{p} (380-420MeV)", 24, -4, 4);
    Phip440Prompt = new TH1D("Phip440Prompt", "#phi_{p} (420-460MeV)", 24, -4, 4);
    Phip480Prompt = new TH1D("Phip480Prompt", "#phi_{p} (460-500MeV)", 24, -4, 4);
    Phip520Prompt = new TH1D("Phip520Prompt", "#phi_{p} (500-540MeV)", 24, -4, 4);
    Phip560Prompt = new TH1D("Phip560Prompt", "#phi_{p} (540-580MeV)", 24, -4, 4);
    Phip600Prompt = new TH1D("Phip600Prompt", "#phi_{p} (580-620MeV)", 24, -4, 4);
    Phip640Prompt = new TH1D("Phip640Prompt", "#phi_{p} (620-660MeV)", 24, -4, 4);
    Phip680Prompt = new TH1D("Phip680Prompt", "#phi_{p} (660-700MeV)", 24, -4, 4);

    //    // Angles of neutron in scattered frame across EGamma bins for negative helicity
    Phip265NegHelCM1 = new TH1D( "Phip_265MeV_NegHelCM1", "#phi_{p} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 20, -4, 4);
    Phip335NegHelCM1 = new TH1D( "Phip_335MeV_NegHelCM1", "#phi_{p} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 20, -4, 4);
    Phip405NegHelCM1 = new TH1D( "Phip_405MeV_NegHelCM1", "#phi_{p} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 20, -4, 4);
    Phip475NegHelCM1 = new TH1D( "Phip_475MeV_NegHelCM1", "#phi_{p} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 20, -4, 4);
    Phip545NegHelCM1 = new TH1D( "Phip_545MeV_NegHelCM1", "#phi_{p} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 20, -4, 4);
    Phip615NegHelCM1 = new TH1D( "Phip_615MeV_NegHelCM1", "#phi_{p} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 20, -4, 4);
    Phip685NegHelCM1 = new TH1D( "Phip_685MeV_NegHelCM1", "#phi_{p} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 20, -4, 4);

    Phip265NegHelCM2 = new TH1D( "Phip_265MeV_NegHelCM2", "#phi_{p} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 20, -4, 4);
    Phip335NegHelCM2 = new TH1D( "Phip_335MeV_NegHelCM2", "#phi_{p} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 20, -4, 4);
    Phip405NegHelCM2 = new TH1D( "Phip_405MeV_NegHelCM2", "#phi_{p} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 20, -4, 4);
    Phip475NegHelCM2 = new TH1D( "Phip_475MeV_NegHelCM2", "#phi_{p} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 20, -4, 4);
    Phip545NegHelCM2 = new TH1D( "Phip_545MeV_NegHelCM2", "#phi_{p} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 20, -4, 4);
    Phip615NegHelCM2 = new TH1D( "Phip_615MeV_NegHelCM2", "#phi_{p} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 20, -4, 4);
    Phip685NegHelCM2 = new TH1D( "Phip_685MeV_NegHelCM2", "#phi_{p} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 20, -4, 4);

    Phip265NegHelCM3 = new TH1D( "Phip_265MeV_NegHelCM3", "#phi_{p} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 20, -4, 4);
    Phip335NegHelCM3 = new TH1D( "Phip_335MeV_NegHelCM3", "#phi_{p} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 20, -4, 4);
    Phip405NegHelCM3 = new TH1D( "Phip_405MeV_NegHelCM3", "#phi_{p} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 20, -4, 4);
    Phip475NegHelCM3 = new TH1D( "Phip_475MeV_NegHelCM3", "#phi_{p} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 20, -4, 4);
    Phip545NegHelCM3 = new TH1D( "Phip_545MeV_NegHelCM3", "#phi_{p} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 20, -4, 4);
    Phip615NegHelCM3 = new TH1D( "Phip_615MeV_NegHelCM3", "#phi_{p} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 20, -4, 4);
    Phip685NegHelCM3 = new TH1D( "Phip_685MeV_NegHelCM3", "#phi_{p} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 20, -4, 4);

    Phip265NegHelCM4 = new TH1D( "Phip_265MeV_NegHelCM4", "#phi_{p} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 20, -4, 4);
    Phip335NegHelCM4 = new TH1D( "Phip_335MeV_NegHelCM4", "#phi_{p} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 20, -4, 4);
    Phip405NegHelCM4 = new TH1D( "Phip_405MeV_NegHelCM4", "#phi_{p} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 20, -4, 4);
    Phip475NegHelCM4 = new TH1D( "Phip_475MeV_NegHelCM4", "#phi_{p} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 20, -4, 4);
    Phip545NegHelCM4 = new TH1D( "Phip_545MeV_NegHelCM4", "#phi_{p} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 20, -4, 4);
    Phip615NegHelCM4 = new TH1D( "Phip_615MeV_NegHelCM4", "#phi_{p} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 20, -4, 4);
    Phip685NegHelCM4 = new TH1D( "Phip_685MeV_NegHelCM4", "#phi_{p} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 20, -4, 4);

    Phip265NegHelCM5 = new TH1D( "Phip_265MeV_NegHelCM5", "#phi_{p} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for -ve Helicity", 20, -4, 4);
    Phip335NegHelCM5 = new TH1D( "Phip_335MeV_NegHelCM5", "#phi_{p} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for -ve Helicity", 20, -4, 4);
    Phip405NegHelCM5 = new TH1D( "Phip_405MeV_NegHelCM5", "#phi_{p} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for -ve Helicity", 20, -4, 4);
    Phip475NegHelCM5 = new TH1D( "Phip_475MeV_NegHelCM5", "#phi_{p} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for -ve Helicity", 20, -4, 4);
    Phip545NegHelCM5 = new TH1D( "Phip_545MeV_NegHelCM5", "#phi_{p} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for -ve Helicity", 20, -4, 4);
    Phip615NegHelCM5 = new TH1D( "Phip_615MeV_NegHelCM5", "#phi_{p} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for -ve Helicity", 20, -4, 4);
    Phip685NegHelCM5 = new TH1D( "Phip_685MeV_NegHelCM5", "#phi_{p} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for -ve Helicity", 20, -4, 4);

    Phip265NegHelCM6 = new TH1D( "Phip_265MeV_NegHelCM6", "#phi_{p} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 20, -4, 4);
    Phip335NegHelCM6 = new TH1D( "Phip_335MeV_NegHelCM6", "#phi_{p} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 20, -4, 4);
    Phip405NegHelCM6 = new TH1D( "Phip_405MeV_NegHelCM6", "#phi_{p} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 20, -4, 4);
    Phip475NegHelCM6 = new TH1D( "Phip_475MeV_NegHelCM6", "#phi_{p} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 20, -4, 4);
    Phip545NegHelCM6 = new TH1D( "Phip_545MeV_NegHelCM6", "#phi_{p} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 20, -4, 4);
    Phip615NegHelCM6 = new TH1D( "Phip_615MeV_NegHelCM6", "#phi_{p} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 20, -4, 4);
    Phip685NegHelCM6 = new TH1D( "Phip_685MeV_NegHelCM6", "#phi_{p} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 20, -4, 4);

    Phip265NegHelCM7 = new TH1D( "Phip_265MeV_NegHelCM7", "#phi_{p} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 20, -4, 4);
    Phip335NegHelCM7 = new TH1D( "Phip_335MeV_NegHelCM7", "#phi_{p} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 20, -4, 4);
    Phip405NegHelCM7 = new TH1D( "Phip_405MeV_NegHelCM7", "#phi_{p} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 20, -4, 4);
    Phip475NegHelCM7 = new TH1D( "Phip_475MeV_NegHelCM7", "#phi_{p} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 20, -4, 4);
    Phip545NegHelCM7 = new TH1D( "Phip_545MeV_NegHelCM7", "#phi_{p} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 20, -4, 4);
    Phip615NegHelCM7 = new TH1D( "Phip_615MeV_NegHelCM7", "#phi_{p} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 20, -4, 4);
    Phip685NegHelCM7 = new TH1D( "Phip_685MeV_NegHelCM7", "#phi_{p} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 20, -4, 4);

    Phip265NegHelCM8 = new TH1D( "Phip_265MeV_NegHelCM8", "#phi_{p} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 20, -4, 4);
    Phip335NegHelCM8 = new TH1D( "Phip_335MeV_NegHelCM8", "#phi_{p} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 20, -4, 4);
    Phip405NegHelCM8 = new TH1D( "Phip_405MeV_NegHelCM8", "#phi_{p} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 20, -4, 4);
    Phip475NegHelCM8 = new TH1D( "Phip_475MeV_NegHelCM8", "#phi_{p} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 20, -4, 4);
    Phip545NegHelCM8 = new TH1D( "Phip_545MeV_NegHelCM8", "#phi_{p} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 20, -4, 4);
    Phip615NegHelCM8 = new TH1D( "Phip_615MeV_NegHelCM8", "#phi_{p} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 20, -4, 4);
    Phip685NegHelCM8 = new TH1D( "Phip_685MeV_NegHelCM8", "#phi_{p} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 20, -4, 4);

    // Angles of neutron in scattered frame across EGamma bins for positive helicity
    Phip265PosHelCM1 = new TH1D( "Phip_265MeV_PosHelCM1", "#phi_{p} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 20, -4, 4);
    Phip335PosHelCM1 = new TH1D( "Phip_335MeV_PosHelCM1", "#phi_{p} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 20, -4, 4);
    Phip405PosHelCM1 = new TH1D( "Phip_405MeV_PosHelCM1", "#phi_{p} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 20, -4, 4);
    Phip475PosHelCM1 = new TH1D( "Phip_475MeV_PosHelCM1", "#phi_{p} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 20, -4, 4);
    Phip545PosHelCM1 = new TH1D( "Phip_545MeV_PosHelCM1", "#phi_{p} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 20, -4, 4);
    Phip615PosHelCM1 = new TH1D( "Phip_615MeV_PosHelCM1", "#phi_{p} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 20, -4, 4);
    Phip685PosHelCM1 = new TH1D( "Phip_685MeV_PosHelCM1", "#phi_{p} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 20, -4, 4);

    Phip265PosHelCM2 = new TH1D( "Phip_265MeV_PosHelCM2", "#phi_{p} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 20, -4, 4);
    Phip335PosHelCM2 = new TH1D( "Phip_335MeV_PosHelCM2", "#phi_{p} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 20, -4, 4);
    Phip405PosHelCM2 = new TH1D( "Phip_405MeV_PosHelCM2", "#phi_{p} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 20, -4, 4);
    Phip475PosHelCM2 = new TH1D( "Phip_475MeV_PosHelCM2", "#phi_{p} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 20, -4, 4);
    Phip545PosHelCM2 = new TH1D( "Phip_545MeV_PosHelCM2", "#phi_{p} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 20, -4, 4);
    Phip615PosHelCM2 = new TH1D( "Phip_615MeV_PosHelCM2", "#phi_{p} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 20, -4, 4);
    Phip685PosHelCM2 = new TH1D( "Phip_685MeV_PosHelCM2", "#phi_{p} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 20, -4, 4);

    Phip265PosHelCM3 = new TH1D( "Phip_265MeV_PosHelCM3", "#phi_{p} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 20, -4, 4);
    Phip335PosHelCM3 = new TH1D( "Phip_335MeV_PosHelCM3", "#phi_{p} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 20, -4, 4);
    Phip405PosHelCM3 = new TH1D( "Phip_405MeV_PosHelCM3", "#phi_{p} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 20, -4, 4);
    Phip475PosHelCM3 = new TH1D( "Phip_475MeV_PosHelCM3", "#phi_{p} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 20, -4, 4);
    Phip545PosHelCM3 = new TH1D( "Phip_545MeV_PosHelCM3", "#phi_{p} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 20, -4, 4);
    Phip615PosHelCM3 = new TH1D( "Phip_615MeV_PosHelCM3", "#phi_{p} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 20, -4, 4);
    Phip685PosHelCM3 = new TH1D( "Phip_685MeV_PosHelCM3", "#phi_{p} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 20, -4, 4);

    Phip265PosHelCM4 = new TH1D( "Phip_265MeV_PosHelCM4", "#phi_{p} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 20, -4, 4);
    Phip335PosHelCM4 = new TH1D( "Phip_335MeV_PosHelCM4", "#phi_{p} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 20, -4, 4);
    Phip405PosHelCM4 = new TH1D( "Phip_405MeV_PosHelCM4", "#phi_{p} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 20, -4, 4);
    Phip475PosHelCM4 = new TH1D( "Phip_475MeV_PosHelCM4", "#phi_{p} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 20, -4, 4);
    Phip545PosHelCM4 = new TH1D( "Phip_545MeV_PosHelCM4", "#phi_{p} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 20, -4, 4);
    Phip615PosHelCM4 = new TH1D( "Phip_615MeV_PosHelCM4", "#phi_{p} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 20, -4, 4);
    Phip685PosHelCM4 = new TH1D( "Phip_685MeV_PosHelCM4", "#phi_{p} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 20, -4, 4);

    Phip265PosHelCM5 = new TH1D( "Phip_265MeV_PosHelCM5", "#phi_{p} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for +ve Helicity", 20, -4, 4);
    Phip335PosHelCM5 = new TH1D( "Phip_335MeV_PosHelCM5", "#phi_{p} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for +ve Helicity", 20, -4, 4);
    Phip405PosHelCM5 = new TH1D( "Phip_405MeV_PosHelCM5", "#phi_{p} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for +ve Helicity", 20, -4, 4);
    Phip475PosHelCM5 = new TH1D( "Phip_475MeV_PosHelCM5", "#phi_{p} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for +ve Helicity", 20, -4, 4);
    Phip545PosHelCM5 = new TH1D( "Phip_545MeV_PosHelCM5", "#phi_{p} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for +ve Helicity", 20, -4, 4);
    Phip615PosHelCM5 = new TH1D( "Phip_615MeV_PosHelCM5", "#phi_{p} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for +ve Helicity", 20, -4, 4);
    Phip685PosHelCM5 = new TH1D( "Phip_685MeV_PosHelCM5", "#phi_{p} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for +ve Helicity", 20, -4, 4);

    Phip265PosHelCM6 = new TH1D( "Phip_265MeV_PosHelCM6", "#phi_{p} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 20, -4, 4);
    Phip335PosHelCM6 = new TH1D( "Phip_335MeV_PosHelCM6", "#phi_{p} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 20, -4, 4);
    Phip405PosHelCM6 = new TH1D( "Phip_405MeV_PosHelCM6", "#phi_{p} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 20, -4, 4);
    Phip475PosHelCM6 = new TH1D( "Phip_475MeV_PosHelCM6", "#phi_{p} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 20, -4, 4);
    Phip545PosHelCM6 = new TH1D( "Phip_545MeV_PosHelCM6", "#phi_{p} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 20, -4, 4);
    Phip615PosHelCM6 = new TH1D( "Phip_615MeV_PosHelCM6", "#phi_{p} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 20, -4, 4);
    Phip685PosHelCM6 = new TH1D( "Phip_685MeV_PosHelCM6", "#phi_{p} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 20, -4, 4);

    Phip265PosHelCM7 = new TH1D( "Phip_265MeV_PosHelCM7", "#phi_{p} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 20, -4, 4);
    Phip335PosHelCM7 = new TH1D( "Phip_335MeV_PosHelCM7", "#phi_{p} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 20, -4, 4);
    Phip405PosHelCM7 = new TH1D( "Phip_405MeV_PosHelCM7", "#phi_{p} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 20, -4, 4);
    Phip475PosHelCM7 = new TH1D( "Phip_475MeV_PosHelCM7", "#phi_{p} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 20, -4, 4);
    Phip545PosHelCM7 = new TH1D( "Phip_545MeV_PosHelCM7", "#phi_{p} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 20, -4, 4);
    Phip615PosHelCM7 = new TH1D( "Phip_615MeV_PosHelCM7", "#phi_{p} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 20, -4, 4);
    Phip685PosHelCM7 = new TH1D( "Phip_685MeV_PosHelCM7", "#phi_{p} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 20, -4, 4);

    Phip265PosHelCM8 = new TH1D( "Phip_265MeV_PosHelCM8", "#phi_{p} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 20, -4, 4);
    Phip335PosHelCM8 = new TH1D( "Phip_335MeV_PosHelCM8", "#phi_{p} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 20, -4, 4);
    Phip405PosHelCM8 = new TH1D( "Phip_405MeV_PosHelCM8", "#phi_{p} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 20, -4, 4);
    Phip475PosHelCM8 = new TH1D( "Phip_475MeV_PosHelCM8", "#phi_{p} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 20, -4, 4);
    Phip545PosHelCM8 = new TH1D( "Phip_545MeV_PosHelCM8", "#phi_{p} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 20, -4, 4);
    Phip615PosHelCM8 = new TH1D( "Phip_615MeV_PosHelCM8", "#phi_{p} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 20, -4, 4);
    Phip685PosHelCM8 = new TH1D( "Phip_685MeV_PosHelCM8", "#phi_{p} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 20, -4, 4);

    Phip265NegHelCM1Prompt = new TH1D( "Phip_265MeV_NegHelCM1Prompt", "#phi_{p} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 20, -4, 4);
    Phip335NegHelCM1Prompt = new TH1D( "Phip_335MeV_NegHelCM1Prompt", "#phi_{p} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 20, -4, 4);
    Phip405NegHelCM1Prompt = new TH1D( "Phip_405MeV_NegHelCM1Prompt", "#phi_{p} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 20, -4, 4);
    Phip475NegHelCM1Prompt = new TH1D( "Phip_475MeV_NegHelCM1Prompt", "#phi_{p} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 20, -4, 4);
    Phip545NegHelCM1Prompt = new TH1D( "Phip_545MeV_NegHelCM1Prompt", "#phi_{p} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 20, -4, 4);
    Phip615NegHelCM1Prompt = new TH1D( "Phip_615MeV_NegHelCM1Prompt", "#phi_{p} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 20, -4, 4);
    Phip685NegHelCM1Prompt = new TH1D( "Phip_685MeV_NegHelCM1Prompt", "#phi_{p} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 20, -4, 4);

    Phip265NegHelCM2Prompt = new TH1D( "Phip_265MeV_NegHelCM2Prompt", "#phi_{p} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 20, -4, 4);
    Phip335NegHelCM2Prompt = new TH1D( "Phip_335MeV_NegHelCM2Prompt", "#phi_{p} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 20, -4, 4);
    Phip405NegHelCM2Prompt = new TH1D( "Phip_405MeV_NegHelCM2Prompt", "#phi_{p} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 20, -4, 4);
    Phip475NegHelCM2Prompt = new TH1D( "Phip_475MeV_NegHelCM2Prompt", "#phi_{p} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 20, -4, 4);
    Phip545NegHelCM2Prompt = new TH1D( "Phip_545MeV_NegHelCM2Prompt", "#phi_{p} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 20, -4, 4);
    Phip615NegHelCM2Prompt = new TH1D( "Phip_615MeV_NegHelCM2Prompt", "#phi_{p} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 20, -4, 4);
    Phip685NegHelCM2Prompt = new TH1D( "Phip_685MeV_NegHelCM2Prompt", "#phi_{p} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 20, -4, 4);

    Phip265NegHelCM3Prompt = new TH1D( "Phip_265MeV_NegHelCM3Prompt", "#phi_{p} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 20, -4, 4);
    Phip335NegHelCM3Prompt = new TH1D( "Phip_335MeV_NegHelCM3Prompt", "#phi_{p} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 20, -4, 4);
    Phip405NegHelCM3Prompt = new TH1D( "Phip_405MeV_NegHelCM3Prompt", "#phi_{p} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 20, -4, 4);
    Phip475NegHelCM3Prompt = new TH1D( "Phip_475MeV_NegHelCM3Prompt", "#phi_{p} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 20, -4, 4);
    Phip545NegHelCM3Prompt = new TH1D( "Phip_545MeV_NegHelCM3Prompt", "#phi_{p} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 20, -4, 4);
    Phip615NegHelCM3Prompt = new TH1D( "Phip_615MeV_NegHelCM3Prompt", "#phi_{p} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 20, -4, 4);
    Phip685NegHelCM3Prompt = new TH1D( "Phip_685MeV_NegHelCM3Prompt", "#phi_{p} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 20, -4, 4);

    Phip265NegHelCM4Prompt = new TH1D( "Phip_265MeV_NegHelCM4Prompt", "#phi_{p} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 20, -4, 4);
    Phip335NegHelCM4Prompt = new TH1D( "Phip_335MeV_NegHelCM4Prompt", "#phi_{p} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 20, -4, 4);
    Phip405NegHelCM4Prompt = new TH1D( "Phip_405MeV_NegHelCM4Prompt", "#phi_{p} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 20, -4, 4);
    Phip475NegHelCM4Prompt = new TH1D( "Phip_475MeV_NegHelCM4Prompt", "#phi_{p} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 20, -4, 4);
    Phip545NegHelCM4Prompt = new TH1D( "Phip_545MeV_NegHelCM4Prompt", "#phi_{p} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 20, -4, 4);
    Phip615NegHelCM4Prompt = new TH1D( "Phip_615MeV_NegHelCM4Prompt", "#phi_{p} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 20, -4, 4);
    Phip685NegHelCM4Prompt = new TH1D( "Phip_685MeV_NegHelCM4Prompt", "#phi_{p} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 20, -4, 4);

    Phip265NegHelCM5Prompt = new TH1D( "Phip_265MeV_NegHelCM5Prompt", "#phi_{p} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for -ve Helicity", 20, -4, 4);
    Phip335NegHelCM5Prompt = new TH1D( "Phip_335MeV_NegHelCM5Prompt", "#phi_{p} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for -ve Helicity", 20, -4, 4);
    Phip405NegHelCM5Prompt = new TH1D( "Phip_405MeV_NegHelCM5Prompt", "#phi_{p} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for -ve Helicity", 20, -4, 4);
    Phip475NegHelCM5Prompt = new TH1D( "Phip_475MeV_NegHelCM5Prompt", "#phi_{p} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for -ve Helicity", 20, -4, 4);
    Phip545NegHelCM5Prompt = new TH1D( "Phip_545MeV_NegHelCM5Prompt", "#phi_{p} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for -ve Helicity", 20, -4, 4);
    Phip615NegHelCM5Prompt = new TH1D( "Phip_615MeV_NegHelCM5Prompt", "#phi_{p} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for -ve Helicity", 20, -4, 4);
    Phip685NegHelCM5Prompt = new TH1D( "Phip_685MeV_NegHelCM5Prompt", "#phi_{p} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for -ve Helicity", 20, -4, 4);

    Phip265NegHelCM6Prompt = new TH1D( "Phip_265MeV_NegHelCM6Prompt", "#phi_{p} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 20, -4, 4);
    Phip335NegHelCM6Prompt = new TH1D( "Phip_335MeV_NegHelCM6Prompt", "#phi_{p} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 20, -4, 4);
    Phip405NegHelCM6Prompt = new TH1D( "Phip_405MeV_NegHelCM6Prompt", "#phi_{p} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 20, -4, 4);
    Phip475NegHelCM6Prompt = new TH1D( "Phip_475MeV_NegHelCM6Prompt", "#phi_{p} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 20, -4, 4);
    Phip545NegHelCM6Prompt = new TH1D( "Phip_545MeV_NegHelCM6Prompt", "#phi_{p} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 20, -4, 4);
    Phip615NegHelCM6Prompt = new TH1D( "Phip_615MeV_NegHelCM6Prompt", "#phi_{p} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 20, -4, 4);
    Phip685NegHelCM6Prompt = new TH1D( "Phip_685MeV_NegHelCM6Prompt", "#phi_{p} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 20, -4, 4);

    Phip265NegHelCM7Prompt = new TH1D( "Phip_265MeV_NegHelCM7Prompt", "#phi_{p} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 20, -4, 4);
    Phip335NegHelCM7Prompt = new TH1D( "Phip_335MeV_NegHelCM7Prompt", "#phi_{p} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 20, -4, 4);
    Phip405NegHelCM7Prompt = new TH1D( "Phip_405MeV_NegHelCM7Prompt", "#phi_{p} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 20, -4, 4);
    Phip475NegHelCM7Prompt = new TH1D( "Phip_475MeV_NegHelCM7Prompt", "#phi_{p} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 20, -4, 4);
    Phip545NegHelCM7Prompt = new TH1D( "Phip_545MeV_NegHelCM7Prompt", "#phi_{p} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 20, -4, 4);
    Phip615NegHelCM7Prompt = new TH1D( "Phip_615MeV_NegHelCM7Prompt", "#phi_{p} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 20, -4, 4);
    Phip685NegHelCM7Prompt = new TH1D( "Phip_685MeV_NegHelCM7Prompt", "#phi_{p} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 20, -4, 4);

    Phip265NegHelCM8Prompt = new TH1D( "Phip_265MeV_NegHelCM8Prompt", "#phi_{p} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 20, -4, 4);
    Phip335NegHelCM8Prompt = new TH1D( "Phip_335MeV_NegHelCM8Prompt", "#phi_{p} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 20, -4, 4);
    Phip405NegHelCM8Prompt = new TH1D( "Phip_405MeV_NegHelCM8Prompt", "#phi_{p} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 20, -4, 4);
    Phip475NegHelCM8Prompt = new TH1D( "Phip_475MeV_NegHelCM8Prompt", "#phi_{p} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 20, -4, 4);
    Phip545NegHelCM8Prompt = new TH1D( "Phip_545MeV_NegHelCM8Prompt", "#phi_{p} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 20, -4, 4);
    Phip615NegHelCM8Prompt = new TH1D( "Phip_615MeV_NegHelCM8Prompt", "#phi_{p} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 20, -4, 4);
    Phip685NegHelCM8Prompt = new TH1D( "Phip_685MeV_NegHelCM8Prompt", "#phi_{p} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 20, -4, 4);

    // Angles of neutron in scattered frame across EGamma bins for positive helicity
    Phip265PosHelCM1Prompt = new TH1D( "Phip_265MeV_PosHelCM1Prompt", "#phi_{p} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 20, -4, 4);
    Phip335PosHelCM1Prompt = new TH1D( "Phip_335MeV_PosHelCM1Prompt", "#phi_{p} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 20, -4, 4);
    Phip405PosHelCM1Prompt = new TH1D( "Phip_405MeV_PosHelCM1Prompt", "#phi_{p} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 20, -4, 4);
    Phip475PosHelCM1Prompt = new TH1D( "Phip_475MeV_PosHelCM1Prompt", "#phi_{p} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 20, -4, 4);
    Phip545PosHelCM1Prompt = new TH1D( "Phip_545MeV_PosHelCM1Prompt", "#phi_{p} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 20, -4, 4);
    Phip615PosHelCM1Prompt = new TH1D( "Phip_615MeV_PosHelCM1Prompt", "#phi_{p} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 20, -4, 4);
    Phip685PosHelCM1Prompt = new TH1D( "Phip_685MeV_PosHelCM1Prompt", "#phi_{p} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 20, -4, 4);

    Phip265PosHelCM2Prompt = new TH1D( "Phip_265MeV_PosHelCM2Prompt", "#phi_{p} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 20, -4, 4);
    Phip335PosHelCM2Prompt = new TH1D( "Phip_335MeV_PosHelCM2Prompt", "#phi_{p} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 20, -4, 4);
    Phip405PosHelCM2Prompt = new TH1D( "Phip_405MeV_PosHelCM2Prompt", "#phi_{p} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 20, -4, 4);
    Phip475PosHelCM2Prompt = new TH1D( "Phip_475MeV_PosHelCM2Prompt", "#phi_{p} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 20, -4, 4);
    Phip545PosHelCM2Prompt = new TH1D( "Phip_545MeV_PosHelCM2Prompt", "#phi_{p} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 20, -4, 4);
    Phip615PosHelCM2Prompt = new TH1D( "Phip_615MeV_PosHelCM2Prompt", "#phi_{p} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 20, -4, 4);
    Phip685PosHelCM2Prompt = new TH1D( "Phip_685MeV_PosHelCM2Prompt", "#phi_{p} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 20, -4, 4);

    Phip265PosHelCM3Prompt = new TH1D( "Phip_265MeV_PosHelCM3Prompt", "#phi_{p} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 20, -4, 4);
    Phip335PosHelCM3Prompt = new TH1D( "Phip_335MeV_PosHelCM3Prompt", "#phi_{p} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 20, -4, 4);
    Phip405PosHelCM3Prompt = new TH1D( "Phip_405MeV_PosHelCM3Prompt", "#phi_{p} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 20, -4, 4);
    Phip475PosHelCM3Prompt = new TH1D( "Phip_475MeV_PosHelCM3Prompt", "#phi_{p} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 20, -4, 4);
    Phip545PosHelCM3Prompt = new TH1D( "Phip_545MeV_PosHelCM3Prompt", "#phi_{p} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 20, -4, 4);
    Phip615PosHelCM3Prompt = new TH1D( "Phip_615MeV_PosHelCM3Prompt", "#phi_{p} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 20, -4, 4);
    Phip685PosHelCM3Prompt = new TH1D( "Phip_685MeV_PosHelCM3Prompt", "#phi_{p} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 20, -4, 4);

    Phip265PosHelCM4Prompt = new TH1D( "Phip_265MeV_PosHelCM4Prompt", "#phi_{p} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 20, -4, 4);
    Phip335PosHelCM4Prompt = new TH1D( "Phip_335MeV_PosHelCM4Prompt", "#phi_{p} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 20, -4, 4);
    Phip405PosHelCM4Prompt = new TH1D( "Phip_405MeV_PosHelCM4Prompt", "#phi_{p} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 20, -4, 4);
    Phip475PosHelCM4Prompt = new TH1D( "Phip_475MeV_PosHelCM4Prompt", "#phi_{p} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 20, -4, 4);
    Phip545PosHelCM4Prompt = new TH1D( "Phip_545MeV_PosHelCM4Prompt", "#phi_{p} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 20, -4, 4);
    Phip615PosHelCM4Prompt = new TH1D( "Phip_615MeV_PosHelCM4Prompt", "#phi_{p} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 20, -4, 4);
    Phip685PosHelCM4Prompt = new TH1D( "Phip_685MeV_PosHelCM4Prompt", "#phi_{p} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 20, -4, 4);

    Phip265PosHelCM5Prompt = new TH1D( "Phip_265MeV_PosHelCM5Prompt", "#phi_{p} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for +ve Helicity", 20, -4, 4);
    Phip335PosHelCM5Prompt = new TH1D( "Phip_335MeV_PosHelCM5Prompt", "#phi_{p} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for +ve Helicity", 20, -4, 4);
    Phip405PosHelCM5Prompt = new TH1D( "Phip_405MeV_PosHelCM5Prompt", "#phi_{p} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for +ve Helicity", 20, -4, 4);
    Phip475PosHelCM5Prompt = new TH1D( "Phip_475MeV_PosHelCM5Prompt", "#phi_{p} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for +ve Helicity", 20, -4, 4);
    Phip545PosHelCM5Prompt = new TH1D( "Phip_545MeV_PosHelCM5Prompt", "#phi_{p} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for +ve Helicity", 20, -4, 4);
    Phip615PosHelCM5Prompt = new TH1D( "Phip_615MeV_PosHelCM5Prompt", "#phi_{p} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for +ve Helicity", 20, -4, 4);
    Phip685PosHelCM5Prompt = new TH1D( "Phip_685MeV_PosHelCM5Prompt", "#phi_{p} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for +ve Helicity", 20, -4, 4);

    Phip265PosHelCM6Prompt = new TH1D( "Phip_265MeV_PosHelCM6Prompt", "#phi_{p} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 20, -4, 4);
    Phip335PosHelCM6Prompt = new TH1D( "Phip_335MeV_PosHelCM6Prompt", "#phi_{p} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 20, -4, 4);
    Phip405PosHelCM6Prompt = new TH1D( "Phip_405MeV_PosHelCM6Prompt", "#phi_{p} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 20, -4, 4);
    Phip475PosHelCM6Prompt = new TH1D( "Phip_475MeV_PosHelCM6Prompt", "#phi_{p} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 20, -4, 4);
    Phip545PosHelCM6Prompt = new TH1D( "Phip_545MeV_PosHelCM6Prompt", "#phi_{p} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 20, -4, 4);
    Phip615PosHelCM6Prompt = new TH1D( "Phip_615MeV_PosHelCM6Prompt", "#phi_{p} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 20, -4, 4);
    Phip685PosHelCM6Prompt = new TH1D( "Phip_685MeV_PosHelCM6Prompt", "#phi_{p} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 20, -4, 4);

    Phip265PosHelCM7Prompt = new TH1D( "Phip_265MeV_PosHelCM7Prompt", "#phi_{p} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 20, -4, 4);
    Phip335PosHelCM7Prompt = new TH1D( "Phip_335MeV_PosHelCM7Prompt", "#phi_{p} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 20, -4, 4);
    Phip405PosHelCM7Prompt = new TH1D( "Phip_405MeV_PosHelCM7Prompt", "#phi_{p} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 20, -4, 4);
    Phip475PosHelCM7Prompt = new TH1D( "Phip_475MeV_PosHelCM7Prompt", "#phi_{p} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 20, -4, 4);
    Phip545PosHelCM7Prompt = new TH1D( "Phip_545MeV_PosHelCM7Prompt", "#phi_{p} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 20, -4, 4);
    Phip615PosHelCM7Prompt = new TH1D( "Phip_615MeV_PosHelCM7Prompt", "#phi_{p} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 20, -4, 4);
    Phip685PosHelCM7Prompt = new TH1D( "Phip_685MeV_PosHelCM7Prompt", "#phi_{p} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 20, -4, 4);

    Phip265PosHelCM8Prompt = new TH1D( "Phip_265MeV_PosHelCM8Prompt", "#phi_{p} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 20, -4, 4);
    Phip335PosHelCM8Prompt = new TH1D( "Phip_335MeV_PosHelCM8Prompt", "#phi_{p} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 20, -4, 4);
    Phip405PosHelCM8Prompt = new TH1D( "Phip_405MeV_PosHelCM8Prompt", "#phi_{p} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 20, -4, 4);
    Phip475PosHelCM8Prompt = new TH1D( "Phip_475MeV_PosHelCM8Prompt", "#phi_{p} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 20, -4, 4);
    Phip545PosHelCM8Prompt = new TH1D( "Phip_545MeV_PosHelCM8Prompt", "#phi_{p} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 20, -4, 4);
    Phip615PosHelCM8Prompt = new TH1D( "Phip_615MeV_PosHelCM8Prompt", "#phi_{p} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 20, -4, 4);
    Phip685PosHelCM8Prompt = new TH1D( "Phip_685MeV_PosHelCM8Prompt", "#phi_{p} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 20, -4, 4);

    Phip265NegHelCM1Random = new TH1D( "Phip_265MeV_NegHelCM1Random", "#phi_{p} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 20, -4, 4);
    Phip335NegHelCM1Random = new TH1D( "Phip_335MeV_NegHelCM1Random", "#phi_{p} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 20, -4, 4);
    Phip405NegHelCM1Random = new TH1D( "Phip_405MeV_NegHelCM1Random", "#phi_{p} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 20, -4, 4);
    Phip475NegHelCM1Random = new TH1D( "Phip_475MeV_NegHelCM1Random", "#phi_{p} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 20, -4, 4);
    Phip545NegHelCM1Random = new TH1D( "Phip_545MeV_NegHelCM1Random", "#phi_{p} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 20, -4, 4);
    Phip615NegHelCM1Random = new TH1D( "Phip_615MeV_NegHelCM1Random", "#phi_{p} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 20, -4, 4);
    Phip685NegHelCM1Random = new TH1D( "Phip_685MeV_NegHelCM1Random", "#phi_{p} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 20, -4, 4);

    Phip265NegHelCM2Random = new TH1D( "Phip_265MeV_NegHelCM2Random", "#phi_{p} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 20, -4, 4);
    Phip335NegHelCM2Random = new TH1D( "Phip_335MeV_NegHelCM2Random", "#phi_{p} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 20, -4, 4);
    Phip405NegHelCM2Random = new TH1D( "Phip_405MeV_NegHelCM2Random", "#phi_{p} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 20, -4, 4);
    Phip475NegHelCM2Random = new TH1D( "Phip_475MeV_NegHelCM2Random", "#phi_{p} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 20, -4, 4);
    Phip545NegHelCM2Random = new TH1D( "Phip_545MeV_NegHelCM2Random", "#phi_{p} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 20, -4, 4);
    Phip615NegHelCM2Random = new TH1D( "Phip_615MeV_NegHelCM2Random", "#phi_{p} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 20, -4, 4);
    Phip685NegHelCM2Random = new TH1D( "Phip_685MeV_NegHelCM2Random", "#phi_{p} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 20, -4, 4);

    Phip265NegHelCM3Random = new TH1D( "Phip_265MeV_NegHelCM3Random", "#phi_{p} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 20, -4, 4);
    Phip335NegHelCM3Random = new TH1D( "Phip_335MeV_NegHelCM3Random", "#phi_{p} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 20, -4, 4);
    Phip405NegHelCM3Random = new TH1D( "Phip_405MeV_NegHelCM3Random", "#phi_{p} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 20, -4, 4);
    Phip475NegHelCM3Random = new TH1D( "Phip_475MeV_NegHelCM3Random", "#phi_{p} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 20, -4, 4);
    Phip545NegHelCM3Random = new TH1D( "Phip_545MeV_NegHelCM3Random", "#phi_{p} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 20, -4, 4);
    Phip615NegHelCM3Random = new TH1D( "Phip_615MeV_NegHelCM3Random", "#phi_{p} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 20, -4, 4);
    Phip685NegHelCM3Random = new TH1D( "Phip_685MeV_NegHelCM3Random", "#phi_{p} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 20, -4, 4);

    Phip265NegHelCM4Random = new TH1D( "Phip_265MeV_NegHelCM4Random", "#phi_{p} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 20, -4, 4);
    Phip335NegHelCM4Random = new TH1D( "Phip_335MeV_NegHelCM4Random", "#phi_{p} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 20, -4, 4);
    Phip405NegHelCM4Random = new TH1D( "Phip_405MeV_NegHelCM4Random", "#phi_{p} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 20, -4, 4);
    Phip475NegHelCM4Random = new TH1D( "Phip_475MeV_NegHelCM4Random", "#phi_{p} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 20, -4, 4);
    Phip545NegHelCM4Random = new TH1D( "Phip_545MeV_NegHelCM4Random", "#phi_{p} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 20, -4, 4);
    Phip615NegHelCM4Random = new TH1D( "Phip_615MeV_NegHelCM4Random", "#phi_{p} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 20, -4, 4);
    Phip685NegHelCM4Random = new TH1D( "Phip_685MeV_NegHelCM4Random", "#phi_{p} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 20, -4, 4);

    Phip265NegHelCM5Random = new TH1D( "Phip_265MeV_NegHelCM5Random", "#phi_{p} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for -ve Helicity", 20, -4, 4);
    Phip335NegHelCM5Random = new TH1D( "Phip_335MeV_NegHelCM5Random", "#phi_{p} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for -ve Helicity", 20, -4, 4);
    Phip405NegHelCM5Random = new TH1D( "Phip_405MeV_NegHelCM5Random", "#phi_{p} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for -ve Helicity", 20, -4, 4);
    Phip475NegHelCM5Random = new TH1D( "Phip_475MeV_NegHelCM5Random", "#phi_{p} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for -ve Helicity", 20, -4, 4);
    Phip545NegHelCM5Random = new TH1D( "Phip_545MeV_NegHelCM5Random", "#phi_{p} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for -ve Helicity", 20, -4, 4);
    Phip615NegHelCM5Random = new TH1D( "Phip_615MeV_NegHelCM5Random", "#phi_{p} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for -ve Helicity", 20, -4, 4);
    Phip685NegHelCM5Random = new TH1D( "Phip_685MeV_NegHelCM5Random", "#phi_{p} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for -ve Helicity", 20, -4, 4);

    Phip265NegHelCM6Random = new TH1D( "Phip_265MeV_NegHelCM6Random", "#phi_{p} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 20, -4, 4);
    Phip335NegHelCM6Random = new TH1D( "Phip_335MeV_NegHelCM6Random", "#phi_{p} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 20, -4, 4);
    Phip405NegHelCM6Random = new TH1D( "Phip_405MeV_NegHelCM6Random", "#phi_{p} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 20, -4, 4);
    Phip475NegHelCM6Random = new TH1D( "Phip_475MeV_NegHelCM6Random", "#phi_{p} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 20, -4, 4);
    Phip545NegHelCM6Random = new TH1D( "Phip_545MeV_NegHelCM6Random", "#phi_{p} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 20, -4, 4);
    Phip615NegHelCM6Random = new TH1D( "Phip_615MeV_NegHelCM6Random", "#phi_{p} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 20, -4, 4);
    Phip685NegHelCM6Random = new TH1D( "Phip_685MeV_NegHelCM6Random", "#phi_{p} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 20, -4, 4);

    Phip265NegHelCM7Random = new TH1D( "Phip_265MeV_NegHelCM7Random", "#phi_{p} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 20, -4, 4);
    Phip335NegHelCM7Random = new TH1D( "Phip_335MeV_NegHelCM7Random", "#phi_{p} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 20, -4, 4);
    Phip405NegHelCM7Random = new TH1D( "Phip_405MeV_NegHelCM7Random", "#phi_{p} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 20, -4, 4);
    Phip475NegHelCM7Random = new TH1D( "Phip_475MeV_NegHelCM7Random", "#phi_{p} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 20, -4, 4);
    Phip545NegHelCM7Random = new TH1D( "Phip_545MeV_NegHelCM7Random", "#phi_{p} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 20, -4, 4);
    Phip615NegHelCM7Random = new TH1D( "Phip_615MeV_NegHelCM7Random", "#phi_{p} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 20, -4, 4);
    Phip685NegHelCM7Random = new TH1D( "Phip_685MeV_NegHelCM7Random", "#phi_{p} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 20, -4, 4);

    Phip265NegHelCM8Random = new TH1D( "Phip_265MeV_NegHelCM8Random", "#phi_{p} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 20, -4, 4);
    Phip335NegHelCM8Random = new TH1D( "Phip_335MeV_NegHelCM8Random", "#phi_{p} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 20, -4, 4);
    Phip405NegHelCM8Random = new TH1D( "Phip_405MeV_NegHelCM8Random", "#phi_{p} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 20, -4, 4);
    Phip475NegHelCM8Random = new TH1D( "Phip_475MeV_NegHelCM8Random", "#phi_{p} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 20, -4, 4);
    Phip545NegHelCM8Random = new TH1D( "Phip_545MeV_NegHelCM8Random", "#phi_{p} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 20, -4, 4);
    Phip615NegHelCM8Random = new TH1D( "Phip_615MeV_NegHelCM8Random", "#phi_{p} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 20, -4, 4);
    Phip685NegHelCM8Random = new TH1D( "Phip_685MeV_NegHelCM8Random", "#phi_{p} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 20, -4, 4);

    // Angles of neutron in scattered frame across EGamma bins for positive helicity
    Phip265PosHelCM1Random = new TH1D( "Phip_265MeV_PosHelCM1Random", "#phi_{p} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 20, -4, 4);
    Phip335PosHelCM1Random = new TH1D( "Phip_335MeV_PosHelCM1Random", "#phi_{p} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 20, -4, 4);
    Phip405PosHelCM1Random = new TH1D( "Phip_405MeV_PosHelCM1Random", "#phi_{p} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 20, -4, 4);
    Phip475PosHelCM1Random = new TH1D( "Phip_475MeV_PosHelCM1Random", "#phi_{p} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 20, -4, 4);
    Phip545PosHelCM1Random = new TH1D( "Phip_545MeV_PosHelCM1Random", "#phi_{p} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 20, -4, 4);
    Phip615PosHelCM1Random = new TH1D( "Phip_615MeV_PosHelCM1Random", "#phi_{p} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 20, -4, 4);
    Phip685PosHelCM1Random = new TH1D( "Phip_685MeV_PosHelCM1Random", "#phi_{p} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 20, -4, 4);

    Phip265PosHelCM2Random = new TH1D( "Phip_265MeV_PosHelCM2Random", "#phi_{p} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 20, -4, 4);
    Phip335PosHelCM2Random = new TH1D( "Phip_335MeV_PosHelCM2Random", "#phi_{p} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 20, -4, 4);
    Phip405PosHelCM2Random = new TH1D( "Phip_405MeV_PosHelCM2Random", "#phi_{p} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 20, -4, 4);
    Phip475PosHelCM2Random = new TH1D( "Phip_475MeV_PosHelCM2Random", "#phi_{p} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 20, -4, 4);
    Phip545PosHelCM2Random = new TH1D( "Phip_545MeV_PosHelCM2Random", "#phi_{p} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 20, -4, 4);
    Phip615PosHelCM2Random = new TH1D( "Phip_615MeV_PosHelCM2Random", "#phi_{p} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 20, -4, 4);
    Phip685PosHelCM2Random = new TH1D( "Phip_685MeV_PosHelCM2Random", "#phi_{p} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 20, -4, 4);

    Phip265PosHelCM3Random = new TH1D( "Phip_265MeV_PosHelCM3Random", "#phi_{p} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 20, -4, 4);
    Phip335PosHelCM3Random = new TH1D( "Phip_335MeV_PosHelCM3Random", "#phi_{p} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 20, -4, 4);
    Phip405PosHelCM3Random = new TH1D( "Phip_405MeV_PosHelCM3Random", "#phi_{p} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 20, -4, 4);
    Phip475PosHelCM3Random = new TH1D( "Phip_475MeV_PosHelCM3Random", "#phi_{p} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 20, -4, 4);
    Phip545PosHelCM3Random = new TH1D( "Phip_545MeV_PosHelCM3Random", "#phi_{p} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 20, -4, 4);
    Phip615PosHelCM3Random = new TH1D( "Phip_615MeV_PosHelCM3Random", "#phi_{p} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 20, -4, 4);
    Phip685PosHelCM3Random = new TH1D( "Phip_685MeV_PosHelCM3Random", "#phi_{p} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 20, -4, 4);

    Phip265PosHelCM4Random = new TH1D( "Phip_265MeV_PosHelCM4Random", "#phi_{p} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 20, -4, 4);
    Phip335PosHelCM4Random = new TH1D( "Phip_335MeV_PosHelCM4Random", "#phi_{p} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 20, -4, 4);
    Phip405PosHelCM4Random = new TH1D( "Phip_405MeV_PosHelCM4Random", "#phi_{p} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 20, -4, 4);
    Phip475PosHelCM4Random = new TH1D( "Phip_475MeV_PosHelCM4Random", "#phi_{p} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 20, -4, 4);
    Phip545PosHelCM4Random = new TH1D( "Phip_545MeV_PosHelCM4Random", "#phi_{p} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 20, -4, 4);
    Phip615PosHelCM4Random = new TH1D( "Phip_615MeV_PosHelCM4Random", "#phi_{p} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 20, -4, 4);
    Phip685PosHelCM4Random = new TH1D( "Phip_685MeV_PosHelCM4Random", "#phi_{p} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 20, -4, 4);

    Phip265PosHelCM5Random = new TH1D( "Phip_265MeV_PosHelCM5Random", "#phi_{p} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for +ve Helicity", 20, -4, 4);
    Phip335PosHelCM5Random = new TH1D( "Phip_335MeV_PosHelCM5Random", "#phi_{p} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for +ve Helicity", 20, -4, 4);
    Phip405PosHelCM5Random = new TH1D( "Phip_405MeV_PosHelCM5Random", "#phi_{p} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for +ve Helicity", 20, -4, 4);
    Phip475PosHelCM5Random = new TH1D( "Phip_475MeV_PosHelCM5Random", "#phi_{p} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for +ve Helicity", 20, -4, 4);
    Phip545PosHelCM5Random = new TH1D( "Phip_545MeV_PosHelCM5Random", "#phi_{p} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for +ve Helicity", 20, -4, 4);
    Phip615PosHelCM5Random = new TH1D( "Phip_615MeV_PosHelCM5Random", "#phi_{p} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for +ve Helicity", 20, -4, 4);
    Phip685PosHelCM5Random = new TH1D( "Phip_685MeV_PosHelCM5Random", "#phi_{p} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for +ve Helicity", 20, -4, 4);

    Phip265PosHelCM6Random = new TH1D( "Phip_265MeV_PosHelCM6Random", "#phi_{p} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 20, -4, 4);
    Phip335PosHelCM6Random = new TH1D( "Phip_335MeV_PosHelCM6Random", "#phi_{p} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 20, -4, 4);
    Phip405PosHelCM6Random = new TH1D( "Phip_405MeV_PosHelCM6Random", "#phi_{p} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 20, -4, 4);
    Phip475PosHelCM6Random = new TH1D( "Phip_475MeV_PosHelCM6Random", "#phi_{p} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 20, -4, 4);
    Phip545PosHelCM6Random = new TH1D( "Phip_545MeV_PosHelCM6Random", "#phi_{p} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 20, -4, 4);
    Phip615PosHelCM6Random = new TH1D( "Phip_615MeV_PosHelCM6Random", "#phi_{p} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 20, -4, 4);
    Phip685PosHelCM6Random = new TH1D( "Phip_685MeV_PosHelCM6Random", "#phi_{p} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 20, -4, 4);

    Phip265PosHelCM7Random = new TH1D( "Phip_265MeV_PosHelCM7Random", "#phi_{p} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 20, -4, 4);
    Phip335PosHelCM7Random = new TH1D( "Phip_335MeV_PosHelCM7Random", "#phi_{p} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 20, -4, 4);
    Phip405PosHelCM7Random = new TH1D( "Phip_405MeV_PosHelCM7Random", "#phi_{p} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 20, -4, 4);
    Phip475PosHelCM7Random = new TH1D( "Phip_475MeV_PosHelCM7Random", "#phi_{p} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 20, -4, 4);
    Phip545PosHelCM7Random = new TH1D( "Phip_545MeV_PosHelCM7Random", "#phi_{p} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 20, -4, 4);
    Phip615PosHelCM7Random = new TH1D( "Phip_615MeV_PosHelCM7Random", "#phi_{p} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 20, -4, 4);
    Phip685PosHelCM7Random = new TH1D( "Phip_685MeV_PosHelCM7Random", "#phi_{p} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 20, -4, 4);

    Phip265PosHelCM8Random = new TH1D( "Phip_265MeV_PosHelCM8Random", "#phi_{p} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 20, -4, 4);
    Phip335PosHelCM8Random = new TH1D( "Phip_335MeV_PosHelCM8Random", "#phi_{p} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 20, -4, 4);
    Phip405PosHelCM8Random = new TH1D( "Phip_405MeV_PosHelCM8Random", "#phi_{p} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 20, -4, 4);
    Phip475PosHelCM8Random = new TH1D( "Phip_475MeV_PosHelCM8Random", "#phi_{p} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 20, -4, 4);
    Phip545PosHelCM8Random = new TH1D( "Phip_545MeV_PosHelCM8Random", "#phi_{p} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 20, -4, 4);
    Phip615PosHelCM8Random = new TH1D( "Phip_615MeV_PosHelCM8Random", "#phi_{p} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 20, -4, 4);
    Phip685PosHelCM8Random = new TH1D( "Phip_685MeV_PosHelCM8Random", "#phi_{p} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 20, -4, 4);
}

void PNeutPol_Polarimeter_Lin_NoScatt::FillHists()
{
    time->Fill(TaggerTime);
    if (-5 < TaggerTime && TaggerTime < 20) time_cut->Fill(TaggerTime);

    Eg->Fill(EGamma, TaggerTime);
    OAngle->Fill(OpeningAngle, TaggerTime);
    MMpEpCorrected->Fill(MMpEpCorr, TaggerTime);
    ZpDist->Fill(Zp, TaggerTime);
    ThetanDist->Fill(Thn, TaggerTime);

    E_dE->Fill(EpCorr, dEp, TaggerTime);
    DeutKinPiKin->Fill(ThetanRec-Thn, ThetaPiRecDiff, TaggerTime);

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

    if ( 410 < EGamma && EGamma < 420) {

        if(1 > CosThetapCM && CosThetapCM > 0.9 ){
            Phip415CM1->Fill(Php, TaggerTime);
        }

        else if(0.9 > CosThetapCM && CosThetapCM > 0.8){
            Phip415CM2->Fill(Php, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.7){
            Phip415CM3->Fill(Php, TaggerTime);
        }

        else if(0.7 > CosThetapCM && CosThetapCM > 0.6){
            Phip415CM4->Fill(Php, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.5){
            Phip415CM5->Fill(Php, TaggerTime);
        }

        else if(0.5 > CosThetapCM && CosThetapCM > 0.4){
            Phip415CM6->Fill(Php, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.3){
            Phip415CM7->Fill(Php, TaggerTime);
        }

        else if(0.3 > CosThetapCM && CosThetapCM > 0.2){
            Phip415CM8->Fill(Php, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0.1){
            Phip415CM9->Fill(Php, TaggerTime);
        }

        else if(0.1 > CosThetapCM && CosThetapCM > 0){
            Phip415CM10->Fill(Php, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.1){
            Phip415CM11->Fill(Php, TaggerTime);
        }

        else if(-0.1 > CosThetapCM && CosThetapCM > -0.2){
            Phip415CM12->Fill(Php, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.3){
            Phip415CM13->Fill(Php, TaggerTime);
        }

        else if(-0.3 > CosThetapCM && CosThetapCM > -0.4){
            Phip415CM14->Fill(Php, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.5){
            Phip415CM15->Fill(Php, TaggerTime);
        }

        else if(-0.5 > CosThetapCM && CosThetapCM > -0.6){
            Phip415CM16->Fill(Php, TaggerTime);
        }

        else if(-0.6 > CosThetapCM && CosThetapCM > -0.7){
            Phip415CM17->Fill(Php, TaggerTime);
        }

        else if(-0.7 > CosThetapCM && CosThetapCM > -0.8){
            Phip415CM18->Fill(Php, TaggerTime);
        }

        else if(-0.8 > CosThetapCM && CosThetapCM > -0.9){
            Phip415CM19->Fill(Php, TaggerTime);
        }

        else if(-0.9 > CosThetapCM && CosThetapCM > -1){
            Phip415CM20->Fill(Php, TaggerTime);
        }
    }


    if ( 420 < EGamma && EGamma < 430) {

        if(1 > CosThetapCM && CosThetapCM > 0.9 ){
            Phip425CM1->Fill(Php, TaggerTime);
        }

        else if(0.9 > CosThetapCM && CosThetapCM > 0.8){
            Phip425CM2->Fill(Php, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.7){
            Phip425CM3->Fill(Php, TaggerTime);
        }

        else if(0.7 > CosThetapCM && CosThetapCM > 0.6){
            Phip425CM4->Fill(Php, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.5){
            Phip425CM5->Fill(Php, TaggerTime);
        }

        else if(0.5 > CosThetapCM && CosThetapCM > 0.4){
            Phip425CM6->Fill(Php, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.3){
            Phip425CM7->Fill(Php, TaggerTime);
        }

        else if(0.3 > CosThetapCM && CosThetapCM > 0.2){
            Phip425CM8->Fill(Php, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0.1){
            Phip425CM9->Fill(Php, TaggerTime);
        }

        else if(0.1 > CosThetapCM && CosThetapCM > 0){
            Phip425CM10->Fill(Php, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.1){
            Phip425CM11->Fill(Php, TaggerTime);
        }

        else if(-0.1 > CosThetapCM && CosThetapCM > -0.2){
            Phip425CM12->Fill(Php, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.3){
            Phip425CM13->Fill(Php, TaggerTime);
        }

        else if(-0.3 > CosThetapCM && CosThetapCM > -0.4){
            Phip425CM14->Fill(Php, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.5){
            Phip425CM15->Fill(Php, TaggerTime);
        }

        else if(-0.5 > CosThetapCM && CosThetapCM > -0.6){
            Phip425CM16->Fill(Php, TaggerTime);
        }

        else if(-0.6 > CosThetapCM && CosThetapCM > -0.7){
            Phip425CM17->Fill(Php, TaggerTime);
        }

        else if(-0.7 > CosThetapCM && CosThetapCM > -0.8){
            Phip425CM18->Fill(Php, TaggerTime);
        }

        else if(-0.8 > CosThetapCM && CosThetapCM > -0.9){
            Phip425CM19->Fill(Php, TaggerTime);
        }

        else if(-0.9 > CosThetapCM && CosThetapCM > -1){
            Phip425CM20->Fill(Php, TaggerTime);
        }
    }

    else if ( 430 < EGamma && EGamma < 440) {

        if(1 > CosThetapCM && CosThetapCM > 0.9 ){
            Phip435CM1->Fill(Php, TaggerTime);
        }

        else if(0.9 > CosThetapCM && CosThetapCM > 0.8){
            Phip435CM2->Fill(Php, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.7){
            Phip435CM3->Fill(Php, TaggerTime);
        }

        else if(0.7 > CosThetapCM && CosThetapCM > 0.6){
            Phip435CM4->Fill(Php, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.5){
            Phip435CM5->Fill(Php, TaggerTime);
        }

        else if(0.5 > CosThetapCM && CosThetapCM > 0.4){
            Phip435CM6->Fill(Php, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.3){
            Phip435CM7->Fill(Php, TaggerTime);
        }

        else if(0.3 > CosThetapCM && CosThetapCM > 0.2){
            Phip435CM8->Fill(Php, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0.1){
            Phip435CM9->Fill(Php, TaggerTime);
        }

        else if(0.1 > CosThetapCM && CosThetapCM > 0){
            Phip435CM10->Fill(Php, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.1){
            Phip435CM11->Fill(Php, TaggerTime);
        }

        else if(-0.1 > CosThetapCM && CosThetapCM > -0.2){
            Phip435CM12->Fill(Php, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.3){
            Phip435CM13->Fill(Php, TaggerTime);
        }

        else if(-0.3 > CosThetapCM && CosThetapCM > -0.4){
            Phip435CM14->Fill(Php, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.5){
            Phip435CM15->Fill(Php, TaggerTime);
        }

        else if(-0.5 > CosThetapCM && CosThetapCM > -0.6){
            Phip435CM16->Fill(Php, TaggerTime);
        }

        else if(-0.6 > CosThetapCM && CosThetapCM > -0.7){
            Phip435CM17->Fill(Php, TaggerTime);
        }

        else if(-0.7 > CosThetapCM && CosThetapCM > -0.8){
            Phip435CM18->Fill(Php, TaggerTime);
        }

        else if(-0.8 > CosThetapCM && CosThetapCM > -0.9){
            Phip435CM19->Fill(Php, TaggerTime);
        }

        else if(-0.9 > CosThetapCM && CosThetapCM > -1){
            Phip435CM20->Fill(Php, TaggerTime);
        }
    }

    else if ( 440 < EGamma && EGamma < 450) {

        if(1 > CosThetapCM && CosThetapCM > 0.9 ){
            Phip445CM1->Fill(Php, TaggerTime);
        }

        else if(0.9 > CosThetapCM && CosThetapCM > 0.8){
            Phip445CM2->Fill(Php, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.7){
            Phip445CM3->Fill(Php, TaggerTime);
        }

        else if(0.7 > CosThetapCM && CosThetapCM > 0.6){
            Phip445CM4->Fill(Php, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.5){
            Phip445CM5->Fill(Php, TaggerTime);
        }

        else if(0.5 > CosThetapCM && CosThetapCM > 0.4){
            Phip445CM6->Fill(Php, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.3){
            Phip445CM7->Fill(Php, TaggerTime);
        }

        else if(0.3 > CosThetapCM && CosThetapCM > 0.2){
            Phip445CM8->Fill(Php, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0.1){
            Phip445CM9->Fill(Php, TaggerTime);
        }

        else if(0.1 > CosThetapCM && CosThetapCM > 0){
            Phip445CM10->Fill(Php, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.1){
            Phip445CM11->Fill(Php, TaggerTime);
        }

        else if(-0.1 > CosThetapCM && CosThetapCM > -0.2){
            Phip445CM12->Fill(Php, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.3){
            Phip445CM13->Fill(Php, TaggerTime);
        }

        else if(-0.3 > CosThetapCM && CosThetapCM > -0.4){
            Phip445CM14->Fill(Php, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.5){
            Phip445CM15->Fill(Php, TaggerTime);
        }

        else if(-0.5 > CosThetapCM && CosThetapCM > -0.6){
            Phip445CM16->Fill(Php, TaggerTime);
        }

        else if(-0.6 > CosThetapCM && CosThetapCM > -0.7){
            Phip445CM17->Fill(Php, TaggerTime);
        }

        else if(-0.7 > CosThetapCM && CosThetapCM > -0.8){
            Phip445CM18->Fill(Php, TaggerTime);
        }

        else if(-0.8 > CosThetapCM && CosThetapCM > -0.9){
            Phip445CM19->Fill(Php, TaggerTime);
        }

        else if(-0.9 > CosThetapCM && CosThetapCM > -1){
            Phip445CM20->Fill(Php, TaggerTime);
        }
    }

    else if ( 450 < EGamma && EGamma < 460) {

        if(1 > CosThetapCM && CosThetapCM > 0.9 ){
            Phip455CM1->Fill(Php, TaggerTime);
        }

        else if(0.9 > CosThetapCM && CosThetapCM > 0.8){
            Phip455CM2->Fill(Php, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.7){
            Phip455CM3->Fill(Php, TaggerTime);
        }

        else if(0.7 > CosThetapCM && CosThetapCM > 0.6){
            Phip455CM4->Fill(Php, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.5){
            Phip455CM5->Fill(Php, TaggerTime);
        }

        else if(0.5 > CosThetapCM && CosThetapCM > 0.4){
            Phip455CM6->Fill(Php, TaggerTime);
        }
        else if(0.4 > CosThetapCM && CosThetapCM > 0.3){
            Phip455CM7->Fill(Php, TaggerTime);
        }

        else if(0.3 > CosThetapCM && CosThetapCM > 0.2){
            Phip455CM8->Fill(Php, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0.1){
            Phip455CM9->Fill(Php, TaggerTime);
        }

        else if(0.1 > CosThetapCM && CosThetapCM > 0){
            Phip455CM10->Fill(Php, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.1){
            Phip455CM11->Fill(Php, TaggerTime);
        }

        else if(-0.1 > CosThetapCM && CosThetapCM > -0.2){
            Phip455CM12->Fill(Php, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.3){
            Phip455CM13->Fill(Php, TaggerTime);
        }

        else if(-0.3 > CosThetapCM && CosThetapCM > -0.4){
            Phip455CM14->Fill(Php, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.5){
            Phip455CM15->Fill(Php, TaggerTime);
        }

        else if(-0.5 > CosThetapCM && CosThetapCM > -0.6){
            Phip455CM16->Fill(Php, TaggerTime);
        }

        else if(-0.6 > CosThetapCM && CosThetapCM > -0.7){
            Phip455CM17->Fill(Php, TaggerTime);
        }

        else if(-0.7 > CosThetapCM && CosThetapCM > -0.8){
            Phip455CM18->Fill(Php, TaggerTime);
        }

        else if(-0.8 > CosThetapCM && CosThetapCM > -0.9){
            Phip455CM19->Fill(Php, TaggerTime);
        }

        else if(-0.9 > CosThetapCM && CosThetapCM > -1){
            Phip455CM20->Fill(Php, TaggerTime);
        }
    }

    else if ( 460 < EGamma && EGamma < 470) {

        if(1 > CosThetapCM && CosThetapCM > 0.9 ){
            Phip465CM1->Fill(Php, TaggerTime);
        }

        else if(0.9 > CosThetapCM && CosThetapCM > 0.8){
            Phip465CM2->Fill(Php, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.7){
            Phip465CM3->Fill(Php, TaggerTime);
        }

        else if(0.7 > CosThetapCM && CosThetapCM > 0.6){
            Phip465CM4->Fill(Php, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.5){
            Phip465CM5->Fill(Php, TaggerTime);
        }

        else if(0.5 > CosThetapCM && CosThetapCM > 0.4){
            Phip465CM6->Fill(Php, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.3){
            Phip465CM7->Fill(Php, TaggerTime);
        }

        else if(0.3 > CosThetapCM && CosThetapCM > 0.2){
            Phip465CM8->Fill(Php, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0.1){
            Phip465CM9->Fill(Php, TaggerTime);
        }

        else if(0.1 > CosThetapCM && CosThetapCM > 0){
            Phip465CM10->Fill(Php, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.1){
            Phip465CM11->Fill(Php, TaggerTime);
        }

        else if(-0.1 > CosThetapCM && CosThetapCM > -0.2){
            Phip465CM12->Fill(Php, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.3){
            Phip465CM13->Fill(Php, TaggerTime);
        }

        else if(-0.3 > CosThetapCM && CosThetapCM > -0.4){
            Phip465CM14->Fill(Php, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.5){
            Phip465CM15->Fill(Php, TaggerTime);
        }

        else if(-0.5 > CosThetapCM && CosThetapCM > -0.6){
            Phip465CM16->Fill(Php, TaggerTime);
        }

        else if(-0.6 > CosThetapCM && CosThetapCM > -0.7){
            Phip465CM17->Fill(Php, TaggerTime);
        }

        else if(-0.7 > CosThetapCM && CosThetapCM > -0.8){
            Phip465CM18->Fill(Php, TaggerTime);
        }

        else if(-0.8 > CosThetapCM && CosThetapCM > -0.9){
            Phip465CM19->Fill(Php, TaggerTime);
        }

        else if(-0.9 > CosThetapCM && CosThetapCM > -1){
            Phip465CM20->Fill(Php, TaggerTime);
        }
    }

    else if ( 470 < EGamma && EGamma < 480) {

        if(1 > CosThetapCM && CosThetapCM > 0.9 ){
            Phip475CM1->Fill(Php, TaggerTime);
        }

        else if(0.9 > CosThetapCM && CosThetapCM > 0.8){
            Phip475CM2->Fill(Php, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.7){
            Phip475CM3->Fill(Php, TaggerTime);
        }

        else if(0.7 > CosThetapCM && CosThetapCM > 0.6){
            Phip475CM4->Fill(Php, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.5){
            Phip475CM5->Fill(Php, TaggerTime);
        }

        else if(0.5 > CosThetapCM && CosThetapCM > 0.4){
            Phip475CM6->Fill(Php, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.3){
            Phip475CM7->Fill(Php, TaggerTime);
        }

        else if(0.3 > CosThetapCM && CosThetapCM > 0.2){
            Phip475CM8->Fill(Php, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0.1){
            Phip475CM9->Fill(Php, TaggerTime);
        }

        else if(0.1 > CosThetapCM && CosThetapCM > 0){
            Phip475CM10->Fill(Php, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.1){
            Phip475CM11->Fill(Php, TaggerTime);
        }

        else if(-0.1 > CosThetapCM && CosThetapCM > -0.2){
            Phip475CM12->Fill(Php, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.3){
            Phip475CM13->Fill(Php, TaggerTime);
        }

        else if(-0.3 > CosThetapCM && CosThetapCM > -0.4){
            Phip475CM14->Fill(Php, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.5){
            Phip475CM15->Fill(Php, TaggerTime);
        }

        else if(-0.5 > CosThetapCM && CosThetapCM > -0.6){
            Phip475CM16->Fill(Php, TaggerTime);
        }

        else if(-0.6 > CosThetapCM && CosThetapCM > -0.7){
            Phip475CM17->Fill(Php, TaggerTime);
        }

        else if(-0.7 > CosThetapCM && CosThetapCM > -0.8){
            Phip475CM18->Fill(Php, TaggerTime);
        }

        else if(-0.8 > CosThetapCM && CosThetapCM > -0.9){
            Phip475CM19->Fill(Php, TaggerTime);
        }

        else if(-0.9 > CosThetapCM && CosThetapCM > -1){
            Phip475CM20->Fill(Php, TaggerTime);
        }
    }

    else if ( 480 < EGamma && EGamma < 490) {

        if(1 > CosThetapCM && CosThetapCM > 0.9 ){
            Phip485CM1->Fill(Php, TaggerTime);
        }

        else if(0.9 > CosThetapCM && CosThetapCM > 0.8){
            Phip485CM2->Fill(Php, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.7){
            Phip485CM3->Fill(Php, TaggerTime);
        }

        else if(0.7 > CosThetapCM && CosThetapCM > 0.6){
            Phip485CM4->Fill(Php, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.5){
            Phip485CM5->Fill(Php, TaggerTime);
        }

        else if(0.5 > CosThetapCM && CosThetapCM > 0.4){
            Phip485CM6->Fill(Php, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.3){
            Phip485CM7->Fill(Php, TaggerTime);
        }

        else if(0.3 > CosThetapCM && CosThetapCM > 0.2){
            Phip485CM8->Fill(Php, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0.1){
            Phip485CM9->Fill(Php, TaggerTime);
        }

        else if(0.1 > CosThetapCM && CosThetapCM > 0){
            Phip485CM10->Fill(Php, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.1){
            Phip485CM11->Fill(Php, TaggerTime);
        }

        else if(-0.1 > CosThetapCM && CosThetapCM > -0.2){
            Phip485CM12->Fill(Php, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.3){
            Phip485CM13->Fill(Php, TaggerTime);
        }

        else if(-0.3 > CosThetapCM && CosThetapCM > -0.4){
            Phip485CM14->Fill(Php, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.5){
            Phip485CM15->Fill(Php, TaggerTime);
        }

        else if(-0.5 > CosThetapCM && CosThetapCM > -0.6){
            Phip485CM16->Fill(Php, TaggerTime);
        }

        else if(-0.6 > CosThetapCM && CosThetapCM > -0.7){
            Phip485CM17->Fill(Php, TaggerTime);
        }

        else if(-0.7 > CosThetapCM && CosThetapCM > -0.8){
            Phip485CM18->Fill(Php, TaggerTime);
        }

        else if(-0.8 > CosThetapCM && CosThetapCM > -0.9){
            Phip485CM19->Fill(Php, TaggerTime);
        }

        else if(-0.9 > CosThetapCM && CosThetapCM > -1){
            Phip485CM20->Fill(Php, TaggerTime);
        }
    }

    else if ( 490 < EGamma && EGamma < 500) {

        if(1 > CosThetapCM && CosThetapCM > 0.9 ){
            Phip495CM1->Fill(Php, TaggerTime);
        }

        else if(0.9 > CosThetapCM && CosThetapCM > 0.8){
            Phip495CM2->Fill(Php, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.7){
            Phip495CM3->Fill(Php, TaggerTime);
        }

        else if(0.7 > CosThetapCM && CosThetapCM > 0.6){
            Phip495CM4->Fill(Php, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.5){
            Phip495CM5->Fill(Php, TaggerTime);
        }

        else if(0.5 > CosThetapCM && CosThetapCM > 0.4){
            Phip495CM6->Fill(Php, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.3){
            Phip495CM7->Fill(Php, TaggerTime);
        }

        else if(0.3 > CosThetapCM && CosThetapCM > 0.2){
            Phip495CM8->Fill(Php, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0.1){
            Phip495CM9->Fill(Php, TaggerTime);
        }

        else if(0.1 > CosThetapCM && CosThetapCM > 0){
            Phip495CM10->Fill(Php, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.1){
            Phip495CM11->Fill(Php, TaggerTime);
        }

        else if(-0.1 > CosThetapCM && CosThetapCM > -0.2){
            Phip495CM12->Fill(Php, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.3){
            Phip495CM13->Fill(Php, TaggerTime);
        }

        else if(-0.3 > CosThetapCM && CosThetapCM > -0.4){
            Phip495CM14->Fill(Php, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.5){
            Phip495CM15->Fill(Php, TaggerTime);
        }

        else if(-0.5 > CosThetapCM && CosThetapCM > -0.6){
            Phip495CM16->Fill(Php, TaggerTime);
        }

        else if(-0.6 > CosThetapCM && CosThetapCM > -0.7){
            Phip495CM17->Fill(Php, TaggerTime);
        }

        else if(-0.7 > CosThetapCM && CosThetapCM > -0.8){
            Phip495CM18->Fill(Php, TaggerTime);
        }

        else if(-0.8 > CosThetapCM && CosThetapCM > -0.9){
            Phip495CM19->Fill(Php, TaggerTime);
        }

        else if(-0.9 > CosThetapCM && CosThetapCM > -1){
            Phip495CM20->Fill(Php, TaggerTime);
        }
    }

    else if ( 500 < EGamma && EGamma < 510) {

        if(1 > CosThetapCM && CosThetapCM > 0.9 ){
            Phip505CM1->Fill(Php, TaggerTime);
        }

        else if(0.9 > CosThetapCM && CosThetapCM > 0.8){
            Phip505CM2->Fill(Php, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.7){
            Phip505CM3->Fill(Php, TaggerTime);
        }

        else if(0.7 > CosThetapCM && CosThetapCM > 0.6){
            Phip505CM4->Fill(Php, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.5){
            Phip505CM5->Fill(Php, TaggerTime);
        }

        else if(0.5 > CosThetapCM && CosThetapCM > 0.4){
            Phip505CM6->Fill(Php, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.3){
            Phip505CM7->Fill(Php, TaggerTime);
        }

        else if(0.3 > CosThetapCM && CosThetapCM > 0.2){
            Phip505CM8->Fill(Php, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0.1){
            Phip505CM9->Fill(Php, TaggerTime);
        }

        else if(0.1 > CosThetapCM && CosThetapCM > 0){
            Phip505CM10->Fill(Php, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.1){
            Phip505CM11->Fill(Php, TaggerTime);
        }

        else if(-0.1 > CosThetapCM && CosThetapCM > -0.2){
            Phip505CM12->Fill(Php, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.3){
            Phip505CM13->Fill(Php, TaggerTime);
        }

        else if(-0.3 > CosThetapCM && CosThetapCM > -0.4){
            Phip505CM14->Fill(Php, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.5){
            Phip505CM15->Fill(Php, TaggerTime);
        }

        else if(-0.5 > CosThetapCM && CosThetapCM > -0.6){
            Phip505CM16->Fill(Php, TaggerTime);
        }

        else if(-0.6 > CosThetapCM && CosThetapCM > -0.7){
            Phip505CM17->Fill(Php, TaggerTime);
        }

        else if(-0.7 > CosThetapCM && CosThetapCM > -0.8){
            Phip505CM18->Fill(Php, TaggerTime);
        }

        else if(-0.8 > CosThetapCM && CosThetapCM > -0.9){
            Phip505CM19->Fill(Php, TaggerTime);
        }

        else if(-0.9 > CosThetapCM && CosThetapCM > -1){
            Phip505CM20->Fill(Php, TaggerTime);
        }
    }

    else if ( 510 < EGamma && EGamma < 520) {

        if(1 > CosThetapCM && CosThetapCM > 0.9 ){
            Phip515CM1->Fill(Php, TaggerTime);
        }

        else if(0.9 > CosThetapCM && CosThetapCM > 0.8){
            Phip515CM2->Fill(Php, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.7){
            Phip515CM3->Fill(Php, TaggerTime);
        }

        else if(0.7 > CosThetapCM && CosThetapCM > 0.6){
            Phip515CM4->Fill(Php, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.5){
            Phip515CM5->Fill(Php, TaggerTime);
        }

        else if(0.5 > CosThetapCM && CosThetapCM > 0.4){
            Phip515CM6->Fill(Php, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.3){
            Phip515CM7->Fill(Php, TaggerTime);
        }

        else if(0.3 > CosThetapCM && CosThetapCM > 0.2){
            Phip515CM8->Fill(Php, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0.1){
            Phip515CM9->Fill(Php, TaggerTime);
        }

        else if(0.1 > CosThetapCM && CosThetapCM > 0){
            Phip515CM10->Fill(Php, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.1){
            Phip515CM11->Fill(Php, TaggerTime);
        }

        else if(-0.1 > CosThetapCM && CosThetapCM > -0.2){
            Phip515CM12->Fill(Php, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.3){
            Phip515CM13->Fill(Php, TaggerTime);
        }

        else if(-0.3 > CosThetapCM && CosThetapCM > -0.4){
            Phip515CM14->Fill(Php, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.5){
            Phip515CM15->Fill(Php, TaggerTime);
        }

        else if(-0.5 > CosThetapCM && CosThetapCM > -0.6){
            Phip515CM16->Fill(Php, TaggerTime);
        }

        else if(-0.6 > CosThetapCM && CosThetapCM > -0.7){
            Phip515CM17->Fill(Php, TaggerTime);
        }

        else if(-0.7 > CosThetapCM && CosThetapCM > -0.8){
            Phip515CM18->Fill(Php, TaggerTime);
        }

        else if(-0.8 > CosThetapCM && CosThetapCM > -0.9){
            Phip515CM19->Fill(Php, TaggerTime);
        }

        else if(-0.9 > CosThetapCM && CosThetapCM > -1){
            Phip515CM20->Fill(Php, TaggerTime);
        }
    }

    else if ( 520 < EGamma && EGamma < 530) {

        if(1 > CosThetapCM && CosThetapCM > 0.9 ){
            Phip525CM1->Fill(Php, TaggerTime);
        }

        else if(0.9 > CosThetapCM && CosThetapCM > 0.8){
            Phip525CM2->Fill(Php, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.7){
            Phip525CM3->Fill(Php, TaggerTime);
        }

        else if(0.7 > CosThetapCM && CosThetapCM > 0.6){
            Phip525CM4->Fill(Php, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.5){
            Phip525CM5->Fill(Php, TaggerTime);
        }

        else if(0.5 > CosThetapCM && CosThetapCM > 0.4){
            Phip525CM6->Fill(Php, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.3){
            Phip525CM7->Fill(Php, TaggerTime);
        }

        else if(0.3 > CosThetapCM && CosThetapCM > 0.2){
            Phip525CM8->Fill(Php, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0.1){
            Phip525CM9->Fill(Php, TaggerTime);
        }

        else if(0.1 > CosThetapCM && CosThetapCM > 0){
            Phip525CM10->Fill(Php, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.1){
            Phip525CM11->Fill(Php, TaggerTime);
        }

        else if(-0.1 > CosThetapCM && CosThetapCM > -0.2){
            Phip525CM12->Fill(Php, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.3){
            Phip525CM13->Fill(Php, TaggerTime);
        }

        else if(-0.3 > CosThetapCM && CosThetapCM > -0.4){
            Phip525CM14->Fill(Php, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.5){
            Phip525CM15->Fill(Php, TaggerTime);
        }

        else if(-0.5 > CosThetapCM && CosThetapCM > -0.6){
            Phip525CM16->Fill(Php, TaggerTime);
        }

        else if(-0.6 > CosThetapCM && CosThetapCM > -0.7){
            Phip525CM17->Fill(Php, TaggerTime);
        }

        else if(-0.7 > CosThetapCM && CosThetapCM > -0.8){
            Phip525CM18->Fill(Php, TaggerTime);
        }

        else if(-0.8 > CosThetapCM && CosThetapCM > -0.9){
            Phip525CM19->Fill(Php, TaggerTime);
        }

        else if(-0.9 > CosThetapCM && CosThetapCM > -1){
            Phip525CM20->Fill(Php, TaggerTime);
        }
    }

    else if ( 530 < EGamma && EGamma < 540) {

        if(1 > CosThetapCM && CosThetapCM > 0.9 ){
            Phip535CM1->Fill(Php, TaggerTime);
        }

        else if(0.9 > CosThetapCM && CosThetapCM > 0.8){
            Phip535CM2->Fill(Php, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.7){
            Phip535CM3->Fill(Php, TaggerTime);
        }

        else if(0.7 > CosThetapCM && CosThetapCM > 0.6){
            Phip535CM4->Fill(Php, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.5){
            Phip535CM5->Fill(Php, TaggerTime);
        }

        else if(0.5 > CosThetapCM && CosThetapCM > 0.4){
            Phip535CM6->Fill(Php, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.3){
            Phip535CM7->Fill(Php, TaggerTime);
        }

        else if(0.3 > CosThetapCM && CosThetapCM > 0.2){
            Phip535CM8->Fill(Php, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0.1){
            Phip535CM9->Fill(Php, TaggerTime);
        }

        else if(0.1 > CosThetapCM && CosThetapCM > 0){
            Phip535CM10->Fill(Php, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.1){
            Phip535CM11->Fill(Php, TaggerTime);
        }

        else if(-0.1 > CosThetapCM && CosThetapCM > -0.2){
            Phip535CM12->Fill(Php, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.3){
            Phip535CM13->Fill(Php, TaggerTime);
        }

        else if(-0.3 > CosThetapCM && CosThetapCM > -0.4){
            Phip535CM14->Fill(Php, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.5){
            Phip535CM15->Fill(Php, TaggerTime);
        }

        else if(-0.5 > CosThetapCM && CosThetapCM > -0.6){
            Phip535CM16->Fill(Php, TaggerTime);
        }

        else if(-0.6 > CosThetapCM && CosThetapCM > -0.7){
            Phip535CM17->Fill(Php, TaggerTime);
        }

        else if(-0.7 > CosThetapCM && CosThetapCM > -0.8){
            Phip535CM18->Fill(Php, TaggerTime);
        }

        else if(-0.8 > CosThetapCM && CosThetapCM > -0.9){
            Phip535CM19->Fill(Php, TaggerTime);
        }

        else if(-0.9 > CosThetapCM && CosThetapCM > -1){
            Phip535CM20->Fill(Php, TaggerTime);
        }
    }

    else if ( 540 < EGamma && EGamma < 550) {

        if(1 > CosThetapCM && CosThetapCM > 0.9 ){
            Phip545CM1->Fill(Php, TaggerTime);
        }

        else if(0.9 > CosThetapCM && CosThetapCM > 0.8){
            Phip545CM2->Fill(Php, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.7){
            Phip545CM3->Fill(Php, TaggerTime);
        }

        else if(0.7 > CosThetapCM && CosThetapCM > 0.6){
            Phip545CM4->Fill(Php, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.5){
            Phip545CM5->Fill(Php, TaggerTime);
        }

        else if(0.5 > CosThetapCM && CosThetapCM > 0.4){
            Phip545CM6->Fill(Php, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.3){
            Phip545CM7->Fill(Php, TaggerTime);
        }

        else if(0.3 > CosThetapCM && CosThetapCM > 0.2){
            Phip545CM8->Fill(Php, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0.1){
            Phip545CM9->Fill(Php, TaggerTime);
        }

        else if(0.1 > CosThetapCM && CosThetapCM > 0){
            Phip545CM10->Fill(Php, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.1){
            Phip545CM11->Fill(Php, TaggerTime);
        }

        else if(-0.1 > CosThetapCM && CosThetapCM > -0.2){
            Phip545CM12->Fill(Php, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.3){
            Phip545CM13->Fill(Php, TaggerTime);
        }

        else if(-0.3 > CosThetapCM && CosThetapCM > -0.4){
            Phip545CM14->Fill(Php, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.5){
            Phip545CM15->Fill(Php, TaggerTime);
        }

        else if(-0.5 > CosThetapCM && CosThetapCM > -0.6){
            Phip545CM16->Fill(Php, TaggerTime);
        }

        else if(-0.6 > CosThetapCM && CosThetapCM > -0.7){
            Phip545CM17->Fill(Php, TaggerTime);
        }

        else if(-0.7 > CosThetapCM && CosThetapCM > -0.8){
            Phip545CM18->Fill(Php, TaggerTime);
        }

        else if(-0.8 > CosThetapCM && CosThetapCM > -0.9){
            Phip545CM19->Fill(Php, TaggerTime);
        }

        else if(-0.9 > CosThetapCM && CosThetapCM > -1){
            Phip545CM20->Fill(Php, TaggerTime);
        }
    }

    else if ( 550 < EGamma && EGamma < 560) {

        if(1 > CosThetapCM && CosThetapCM > 0.9 ){
            Phip555CM1->Fill(Php, TaggerTime);
        }

        else if(0.9 > CosThetapCM && CosThetapCM > 0.8){
            Phip555CM2->Fill(Php, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.7){
            Phip555CM3->Fill(Php, TaggerTime);
        }

        else if(0.7 > CosThetapCM && CosThetapCM > 0.6){
            Phip555CM4->Fill(Php, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.5){
            Phip555CM5->Fill(Php, TaggerTime);
        }

        else if(0.5 > CosThetapCM && CosThetapCM > 0.4){
            Phip555CM6->Fill(Php, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.3){
            Phip555CM7->Fill(Php, TaggerTime);
        }

        else if(0.3 > CosThetapCM && CosThetapCM > 0.2){
            Phip555CM8->Fill(Php, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0.1){
            Phip555CM9->Fill(Php, TaggerTime);
        }

        else if(0.1 > CosThetapCM && CosThetapCM > 0){
            Phip555CM10->Fill(Php, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.1){
            Phip555CM11->Fill(Php, TaggerTime);
        }

        else if(-0.1 > CosThetapCM && CosThetapCM > -0.2){
            Phip555CM12->Fill(Php, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.3){
            Phip555CM13->Fill(Php, TaggerTime);
        }

        else if(-0.3 > CosThetapCM && CosThetapCM > -0.4){
            Phip555CM14->Fill(Php, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.5){
            Phip555CM15->Fill(Php, TaggerTime);
        }

        else if(-0.5 > CosThetapCM && CosThetapCM > -0.6){
            Phip555CM16->Fill(Php, TaggerTime);
        }

        else if(-0.6 > CosThetapCM && CosThetapCM > -0.7){
            Phip555CM17->Fill(Php, TaggerTime);
        }

        else if(-0.7 > CosThetapCM && CosThetapCM > -0.8){
            Phip555CM18->Fill(Php, TaggerTime);
        }

        else if(-0.8 > CosThetapCM && CosThetapCM > -0.9){
            Phip555CM19->Fill(Php, TaggerTime);
        }

        else if(-0.9 > CosThetapCM && CosThetapCM > -1){
            Phip555CM20->Fill(Php, TaggerTime);
        }
    }

    else if ( 560 < EGamma && EGamma < 570) {

        if(1 > CosThetapCM && CosThetapCM > 0.9 ){
            Phip565CM1->Fill(Php, TaggerTime);
        }

        else if(0.9 > CosThetapCM && CosThetapCM > 0.8){
            Phip565CM2->Fill(Php, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.7){
            Phip565CM3->Fill(Php, TaggerTime);
        }

        else if(0.7 > CosThetapCM && CosThetapCM > 0.6){
            Phip565CM4->Fill(Php, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.5){
            Phip565CM5->Fill(Php, TaggerTime);
        }

        else if(0.5 > CosThetapCM && CosThetapCM > 0.4){
            Phip565CM6->Fill(Php, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.3){
            Phip565CM7->Fill(Php, TaggerTime);
        }

        else if(0.3 > CosThetapCM && CosThetapCM > 0.2){
            Phip565CM8->Fill(Php, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0.1){
            Phip565CM9->Fill(Php, TaggerTime);
        }

        else if(0.1 > CosThetapCM && CosThetapCM > 0){
            Phip565CM10->Fill(Php, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.1){
            Phip565CM11->Fill(Php, TaggerTime);
        }

        else if(-0.1 > CosThetapCM && CosThetapCM > -0.2){
            Phip565CM12->Fill(Php, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.3){
            Phip565CM13->Fill(Php, TaggerTime);
        }

        else if(-0.3 > CosThetapCM && CosThetapCM > -0.4){
            Phip565CM14->Fill(Php, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.5){
            Phip565CM15->Fill(Php, TaggerTime);
        }

        else if(-0.5 > CosThetapCM && CosThetapCM > -0.6){
            Phip565CM16->Fill(Php, TaggerTime);
        }

        else if(-0.6 > CosThetapCM && CosThetapCM > -0.7){
            Phip565CM17->Fill(Php, TaggerTime);
        }

        else if(-0.7 > CosThetapCM && CosThetapCM > -0.8){
            Phip565CM18->Fill(Php, TaggerTime);
        }

        else if(-0.8 > CosThetapCM && CosThetapCM > -0.9){
            Phip565CM19->Fill(Php, TaggerTime);
        }

        else if(-0.9 > CosThetapCM && CosThetapCM > -1){
            Phip565CM20->Fill(Php, TaggerTime);
        }
    }

    else if ( 570 < EGamma && EGamma < 580) {

        if(1 > CosThetapCM && CosThetapCM > 0.9 ){
            Phip575CM1->Fill(Php, TaggerTime);
        }

        else if(0.9 > CosThetapCM && CosThetapCM > 0.8){
            Phip575CM2->Fill(Php, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.7){
            Phip575CM3->Fill(Php, TaggerTime);
        }

        else if(0.7 > CosThetapCM && CosThetapCM > 0.6){
            Phip575CM4->Fill(Php, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.5){
            Phip575CM5->Fill(Php, TaggerTime);
        }

        else if(0.5 > CosThetapCM && CosThetapCM > 0.4){
            Phip575CM6->Fill(Php, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.3){
            Phip575CM7->Fill(Php, TaggerTime);
        }

        else if(0.3 > CosThetapCM && CosThetapCM > 0.2){
            Phip575CM8->Fill(Php, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0.1){
            Phip575CM9->Fill(Php, TaggerTime);
        }

        else if(0.1 > CosThetapCM && CosThetapCM > 0){
            Phip575CM10->Fill(Php, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.1){
            Phip575CM11->Fill(Php, TaggerTime);
        }

        else if(-0.1 > CosThetapCM && CosThetapCM > -0.2){
            Phip575CM12->Fill(Php, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.3){
            Phip575CM13->Fill(Php, TaggerTime);
        }

        else if(-0.3 > CosThetapCM && CosThetapCM > -0.4){
            Phip575CM14->Fill(Php, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.5){
            Phip575CM15->Fill(Php, TaggerTime);
        }

        else if(-0.5 > CosThetapCM && CosThetapCM > -0.6){
            Phip575CM16->Fill(Php, TaggerTime);
        }

        else if(-0.6 > CosThetapCM && CosThetapCM > -0.7){
            Phip575CM17->Fill(Php, TaggerTime);
        }

        else if(-0.7 > CosThetapCM && CosThetapCM > -0.8){
            Phip575CM18->Fill(Php, TaggerTime);
        }

        else if(-0.8 > CosThetapCM && CosThetapCM > -0.9){
            Phip575CM19->Fill(Php, TaggerTime);
        }

        else if(-0.9 > CosThetapCM && CosThetapCM > -1){
            Phip575CM20->Fill(Php, TaggerTime);
        }
    }

    else if ( 580 < EGamma && EGamma < 590) {

        if(1 > CosThetapCM && CosThetapCM > 0.9 ){
            Phip585CM1->Fill(Php, TaggerTime);
        }

        else if(0.9 > CosThetapCM && CosThetapCM > 0.8){
            Phip585CM2->Fill(Php, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.7){
            Phip585CM3->Fill(Php, TaggerTime);
        }

        else if(0.7 > CosThetapCM && CosThetapCM > 0.6){
            Phip585CM4->Fill(Php, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.5){
            Phip585CM5->Fill(Php, TaggerTime);
        }

        else if(0.5 > CosThetapCM && CosThetapCM > 0.4){
            Phip585CM6->Fill(Php, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.3){
            Phip585CM7->Fill(Php, TaggerTime);
        }

        else if(0.3 > CosThetapCM && CosThetapCM > 0.2){
            Phip585CM8->Fill(Php, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0.1){
            Phip585CM9->Fill(Php, TaggerTime);
        }

        else if(0.1 > CosThetapCM && CosThetapCM > 0){
            Phip585CM10->Fill(Php, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.1){
            Phip585CM11->Fill(Php, TaggerTime);
        }

        else if(-0.1 > CosThetapCM && CosThetapCM > -0.2){
            Phip585CM12->Fill(Php, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.3){
            Phip585CM13->Fill(Php, TaggerTime);
        }

        else if(-0.3 > CosThetapCM && CosThetapCM > -0.4){
            Phip585CM14->Fill(Php, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.5){
            Phip585CM15->Fill(Php, TaggerTime);
        }

        else if(-0.5 > CosThetapCM && CosThetapCM > -0.6){
            Phip585CM16->Fill(Php, TaggerTime);
        }

        else if(-0.6 > CosThetapCM && CosThetapCM > -0.7){
            Phip585CM17->Fill(Php, TaggerTime);
        }

        else if(-0.7 > CosThetapCM && CosThetapCM > -0.8){
            Phip585CM18->Fill(Php, TaggerTime);
        }

        else if(-0.8 > CosThetapCM && CosThetapCM > -0.9){
            Phip585CM19->Fill(Php, TaggerTime);
        }

        else if(-0.9 > CosThetapCM && CosThetapCM > -1){
            Phip585CM20->Fill(Php, TaggerTime);
        }
    }

    else if ( 590 < EGamma && EGamma < 600) {

        if(1 > CosThetapCM && CosThetapCM > 0.9 ){
            Phip595CM1->Fill(Php, TaggerTime);
        }

        else if(0.9 > CosThetapCM && CosThetapCM > 0.8){
            Phip595CM2->Fill(Php, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.7){
            Phip595CM3->Fill(Php, TaggerTime);
        }

        else if(0.7 > CosThetapCM && CosThetapCM > 0.6){
            Phip595CM4->Fill(Php, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.5){
            Phip595CM5->Fill(Php, TaggerTime);
        }

        else if(0.5 > CosThetapCM && CosThetapCM > 0.4){
            Phip595CM6->Fill(Php, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.3){
            Phip595CM7->Fill(Php, TaggerTime);
        }

        else if(0.3 > CosThetapCM && CosThetapCM > 0.2){
            Phip595CM8->Fill(Php, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0.1){
            Phip595CM9->Fill(Php, TaggerTime);
        }

        else if(0.1 > CosThetapCM && CosThetapCM > 0){
            Phip595CM10->Fill(Php, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.1){
            Phip595CM11->Fill(Php, TaggerTime);
        }

        else if(-0.1 > CosThetapCM && CosThetapCM > -0.2){
            Phip595CM12->Fill(Php, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.3){
            Phip595CM13->Fill(Php, TaggerTime);
        }

        else if(-0.3 > CosThetapCM && CosThetapCM > -0.4){
            Phip595CM14->Fill(Php, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.5){
            Phip595CM15->Fill(Php, TaggerTime);
        }

        else if(-0.5 > CosThetapCM && CosThetapCM > -0.6){
            Phip595CM16->Fill(Php, TaggerTime);
        }

        else if(-0.6 > CosThetapCM && CosThetapCM > -0.7){
            Phip595CM17->Fill(Php, TaggerTime);
        }

        else if(-0.7 > CosThetapCM && CosThetapCM > -0.8){
            Phip595CM18->Fill(Php, TaggerTime);
        }

        else if(-0.8 > CosThetapCM && CosThetapCM > -0.9){
            Phip595CM19->Fill(Php, TaggerTime);
        }

        else if(-0.9 > CosThetapCM && CosThetapCM > -1){
            Phip595CM20->Fill(Php, TaggerTime);
        }
    }

    else if ( 600 < EGamma && EGamma < 610) {

        if(1 > CosThetapCM && CosThetapCM > 0.9 ){
            Phip605CM1->Fill(Php, TaggerTime);
        }

        else if(0.9 > CosThetapCM && CosThetapCM > 0.8){
            Phip605CM2->Fill(Php, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.7){
            Phip605CM3->Fill(Php, TaggerTime);
        }

        else if(0.7 > CosThetapCM && CosThetapCM > 0.6){
            Phip605CM4->Fill(Php, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.5){
            Phip605CM5->Fill(Php, TaggerTime);
        }

        else if(0.5 > CosThetapCM && CosThetapCM > 0.4){
            Phip605CM6->Fill(Php, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.3){
            Phip605CM7->Fill(Php, TaggerTime);
        }

        else if(0.3 > CosThetapCM && CosThetapCM > 0.2){
            Phip605CM8->Fill(Php, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0.1){
            Phip605CM9->Fill(Php, TaggerTime);
        }

        else if(0.1 > CosThetapCM && CosThetapCM > 0){
            Phip605CM10->Fill(Php, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.1){
            Phip605CM11->Fill(Php, TaggerTime);
        }

        else if(-0.1 > CosThetapCM && CosThetapCM > -0.2){
            Phip605CM12->Fill(Php, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.3){
            Phip605CM13->Fill(Php, TaggerTime);
        }

        else if(-0.3 > CosThetapCM && CosThetapCM > -0.4){
            Phip605CM14->Fill(Php, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.5){
            Phip605CM15->Fill(Php, TaggerTime);
        }

        else if(-0.5 > CosThetapCM && CosThetapCM > -0.6){
            Phip605CM16->Fill(Php, TaggerTime);
        }

        else if(-0.6 > CosThetapCM && CosThetapCM > -0.7){
            Phip605CM17->Fill(Php, TaggerTime);
        }

        else if(-0.7 > CosThetapCM && CosThetapCM > -0.8){
            Phip605CM18->Fill(Php, TaggerTime);
        }
        else if(-0.8 > CosThetapCM && CosThetapCM > -0.9){
            Phip605CM19->Fill(Php, TaggerTime);
        }
        else if(-0.9 > CosThetapCM && CosThetapCM > -1){
            Phip605CM20->Fill(Php, TaggerTime);
        }
    }

    else if ( 610 < EGamma && EGamma < 620) {

        if(1 > CosThetapCM && CosThetapCM > 0.9 ){
            Phip615CM1->Fill(Php, TaggerTime);
        }

        else if(0.9 > CosThetapCM && CosThetapCM > 0.8){
            Phip615CM2->Fill(Php, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.7){
            Phip615CM3->Fill(Php, TaggerTime);
        }

        else if(0.7 > CosThetapCM && CosThetapCM > 0.6){
            Phip615CM4->Fill(Php, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.5){
            Phip615CM5->Fill(Php, TaggerTime);
        }

        else if(0.5 > CosThetapCM && CosThetapCM > 0.4){
            Phip615CM6->Fill(Php, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.3){
            Phip615CM7->Fill(Php, TaggerTime);
        }

        else if(0.3 > CosThetapCM && CosThetapCM > 0.2){
            Phip615CM8->Fill(Php, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0.1){
            Phip615CM9->Fill(Php, TaggerTime);
        }

        else if(0.1 > CosThetapCM && CosThetapCM > 0){
            Phip615CM10->Fill(Php, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.1){
            Phip615CM11->Fill(Php, TaggerTime);
        }

        else if(-0.1 > CosThetapCM && CosThetapCM > -0.2){
            Phip615CM12->Fill(Php, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.3){
            Phip615CM13->Fill(Php, TaggerTime);
        }

        else if(-0.3 > CosThetapCM && CosThetapCM > -0.4){
            Phip615CM14->Fill(Php, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.5){
            Phip615CM15->Fill(Php, TaggerTime);
        }

        else if(-0.5 > CosThetapCM && CosThetapCM > -0.6){
            Phip615CM16->Fill(Php, TaggerTime);
        }

        else if(-0.6 > CosThetapCM && CosThetapCM > -0.7){
            Phip615CM17->Fill(Php, TaggerTime);
        }

        else if(-0.7 > CosThetapCM && CosThetapCM > -0.8){
            Phip615CM18->Fill(Php, TaggerTime);
        }

        else if(-0.8 > CosThetapCM && CosThetapCM > -0.9){
            Phip615CM19->Fill(Php, TaggerTime);
        }

        else if(-0.9 > CosThetapCM && CosThetapCM > -1){
            Phip615CM20->Fill(Php, TaggerTime);
        }
    }

    if (GHistBGSub::IsPrompt(TaggerTime) == kTRUE){
        EgPrompt->Fill(EGamma);
        if(ThetanCM  > 80 && ThetanCM < 100){
            if ( 300 < EGamma && EGamma < 340) Phip320Prompt->Fill(PhpRad);
            if ( 340 < EGamma && EGamma < 380) Phip360Prompt->Fill(PhpRad);
            if ( 380 < EGamma && EGamma < 420) Phip400Prompt->Fill(PhpRad);
            if ( 420 < EGamma && EGamma < 460) Phip440Prompt->Fill(PhpRad);
            if ( 460 < EGamma && EGamma < 500) Phip480Prompt->Fill(PhpRad);
            if ( 500 < EGamma && EGamma < 540) Phip520Prompt->Fill(PhpRad);
            if ( 540 < EGamma && EGamma < 580) Phip560Prompt->Fill(PhpRad);
            if ( 580 < EGamma && EGamma < 620) Phip600Prompt->Fill(PhpRad);
            if ( 620 < EGamma && EGamma < 660) Phip640Prompt->Fill(PhpRad);
            if ( 660 < EGamma && EGamma < 700) Phip680Prompt->Fill(PhpRad);

            if ( 230 < EGamma && EGamma < 300) {
                if(1 > CosThetapCM && CosThetapCM > 0.75 ){
                    if (BeamHelicity == kFALSE) Phip265NegHelCM1Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip265PosHelCM1Prompt->Fill(PhpRad);
                }

                else if(0.75 > CosThetapCM && CosThetapCM > 0.5){
                    if (BeamHelicity == kFALSE) Phip265NegHelCM2Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip265PosHelCM2Prompt->Fill(PhpRad);
                }

                else if(0.5 > CosThetapCM && CosThetapCM > 0.25){
                    if (BeamHelicity == kFALSE) Phip265NegHelCM3Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip265PosHelCM3Prompt->Fill(PhpRad);
                }

                else if(0.25 > CosThetapCM && CosThetapCM > 0.0){
                    if (BeamHelicity == kFALSE) Phip265NegHelCM4Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip265PosHelCM4Prompt->Fill(PhpRad);
                }

                else if(0> CosThetapCM && CosThetapCM > -0.25){
                    if (BeamHelicity == kFALSE) Phip265NegHelCM5Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip265PosHelCM5Prompt->Fill(PhpRad);
                }

                else if(-0.25 > CosThetapCM && CosThetapCM > -0.5){
                    if (BeamHelicity == kFALSE) Phip265NegHelCM6Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip265PosHelCM6Prompt->Fill(PhpRad);
                }

                else if(-0.5 > CosThetapCM && CosThetapCM > -0.75){
                    if (BeamHelicity == kFALSE) Phip265NegHelCM7Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip265PosHelCM7Prompt->Fill(PhpRad);
                }

                else if(-0.75 > CosThetapCM && CosThetapCM > -1){
                    if (BeamHelicity == kFALSE) Phip265NegHelCM8Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip265PosHelCM8Prompt->Fill(PhpRad);
                }
            }


            if ( 300 < EGamma && EGamma < 370) {
                if(1 > CosThetapCM && CosThetapCM > 0.75 ){
                    if (BeamHelicity == kFALSE) Phip335NegHelCM1Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip335PosHelCM1Prompt->Fill(PhpRad);
                }

                else if(0.75 > CosThetapCM && CosThetapCM > 0.5){
                    if (BeamHelicity == kFALSE) Phip335NegHelCM2Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip335PosHelCM2Prompt->Fill(PhpRad);
                }

                else if(0.5 > CosThetapCM && CosThetapCM > 0.25){
                    if (BeamHelicity == kFALSE) Phip335NegHelCM3Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip335PosHelCM3Prompt->Fill(PhpRad);
                }

                else if(0.25 > CosThetapCM && CosThetapCM > 0.0){
                    if (BeamHelicity == kFALSE) Phip335NegHelCM4Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip335PosHelCM4Prompt->Fill(PhpRad);
                }

                else if(0> CosThetapCM && CosThetapCM > -0.25){
                    if (BeamHelicity == kFALSE) Phip335NegHelCM5Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip335PosHelCM5Prompt->Fill(PhpRad);
                }

                else if(-0.25 > CosThetapCM && CosThetapCM > -0.5){
                    if (BeamHelicity == kFALSE) Phip335NegHelCM6Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip335PosHelCM6Prompt->Fill(PhpRad);
                }

                else if(-0.5 > CosThetapCM && CosThetapCM > -0.75){
                    if (BeamHelicity == kFALSE) Phip335NegHelCM7Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip335PosHelCM7Prompt->Fill(PhpRad);
                }

                else if(-0.75 > CosThetapCM && CosThetapCM > -1){
                    if (BeamHelicity == kFALSE) Phip335NegHelCM8Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip335PosHelCM8Prompt->Fill(PhpRad);
                }
            }

            if ( 370 < EGamma && EGamma < 440) {
                if(1 > CosThetapCM && CosThetapCM > 0.75 ){
                    if (BeamHelicity == kFALSE) Phip405NegHelCM1Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip405PosHelCM1Prompt->Fill(PhpRad);
                }

                else if(0.75 > CosThetapCM && CosThetapCM > 0.5){
                    if (BeamHelicity == kFALSE) Phip405NegHelCM2Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip405PosHelCM2Prompt->Fill(PhpRad);
                }

                else if(0.5 > CosThetapCM && CosThetapCM > 0.25){
                    if (BeamHelicity == kFALSE) Phip405NegHelCM3Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip405PosHelCM3Prompt->Fill(PhpRad);
                }

                else if(0.25 > CosThetapCM && CosThetapCM > 0.0){
                    if (BeamHelicity == kFALSE) Phip405NegHelCM4Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip405PosHelCM4Prompt->Fill(PhpRad);
                }

                else if(0> CosThetapCM && CosThetapCM > -0.25){
                    if (BeamHelicity == kFALSE) Phip405NegHelCM5Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip405PosHelCM5Prompt->Fill(PhpRad);
                }

                else if(-0.25 > CosThetapCM && CosThetapCM > -0.5){
                    if (BeamHelicity == kFALSE) Phip405NegHelCM6Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip405PosHelCM6Prompt->Fill(PhpRad);
                }

                else if(-0.5 > CosThetapCM && CosThetapCM > -0.75){
                    if (BeamHelicity == kFALSE) Phip405NegHelCM7Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip405PosHelCM7Prompt->Fill(PhpRad);
                }

                else if(-0.75 > CosThetapCM && CosThetapCM > -1){
                    if (BeamHelicity == kFALSE) Phip405NegHelCM8Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip405PosHelCM8Prompt->Fill(PhpRad);
                }
            }

            if ( 440 < EGamma && EGamma < 510) {
                if(1 > CosThetapCM && CosThetapCM > 0.75 ){
                    if (BeamHelicity == kFALSE) Phip475NegHelCM1Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip475PosHelCM1Prompt->Fill(PhpRad);
                }

                else if(0.75 > CosThetapCM && CosThetapCM > 0.5){
                    if (BeamHelicity == kFALSE) Phip475NegHelCM2Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip475PosHelCM2Prompt->Fill(PhpRad);
                }

                else if(0.5 > CosThetapCM && CosThetapCM > 0.25){
                    if (BeamHelicity == kFALSE) Phip475NegHelCM3Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip475PosHelCM3Prompt->Fill(PhpRad);
                }

                else if(0.25 > CosThetapCM && CosThetapCM > 0.0){
                    if (BeamHelicity == kFALSE) Phip475NegHelCM4Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip475PosHelCM4Prompt->Fill(PhpRad);
                }

                else if(0> CosThetapCM && CosThetapCM > -0.25){
                    if (BeamHelicity == kFALSE) Phip475NegHelCM5Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip475PosHelCM5Prompt->Fill(PhpRad);
                }

                else if(-0.25 > CosThetapCM && CosThetapCM > -0.5){
                    if (BeamHelicity == kFALSE) Phip475NegHelCM6Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip475PosHelCM6Prompt->Fill(PhpRad);
                }

                else if(-0.5 > CosThetapCM && CosThetapCM > -0.75){
                    if (BeamHelicity == kFALSE) Phip475NegHelCM7Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip475PosHelCM7Prompt->Fill(PhpRad);
                }

                else if(-0.75 > CosThetapCM && CosThetapCM > -1){
                    if (BeamHelicity == kFALSE) Phip475NegHelCM8Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip475PosHelCM8Prompt->Fill(PhpRad);
                }
            }

            if ( 510 < EGamma && EGamma < 580) {
                if(1 > CosThetapCM && CosThetapCM > 0.75 ){
                    if (BeamHelicity == kFALSE) Phip545NegHelCM1Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip545PosHelCM1Prompt->Fill(PhpRad);
                }

                else if(0.75 > CosThetapCM && CosThetapCM > 0.5){
                    if (BeamHelicity == kFALSE) Phip545NegHelCM2Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip545PosHelCM2Prompt->Fill(PhpRad);
                }

                else if(0.5 > CosThetapCM && CosThetapCM > 0.25){
                    if (BeamHelicity == kFALSE) Phip545NegHelCM3Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip545PosHelCM3Prompt->Fill(PhpRad);
                }

                else if(0.25 > CosThetapCM && CosThetapCM > 0.0){
                    if (BeamHelicity == kFALSE) Phip545NegHelCM4Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip545PosHelCM4Prompt->Fill(PhpRad);
                }

                else if(0> CosThetapCM && CosThetapCM > -0.25){
                    if (BeamHelicity == kFALSE) Phip545NegHelCM5Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip545PosHelCM5Prompt->Fill(PhpRad);
                }

                else if(-0.25 > CosThetapCM && CosThetapCM > -0.5){
                    if (BeamHelicity == kFALSE) Phip545NegHelCM6Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip545PosHelCM6Prompt->Fill(PhpRad);
                }

                else if(-0.5 > CosThetapCM && CosThetapCM > -0.75){
                    if (BeamHelicity == kFALSE) Phip545NegHelCM7Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip545PosHelCM7Prompt->Fill(PhpRad);
                }

                else if(-0.75 > CosThetapCM && CosThetapCM > -1){
                    if (BeamHelicity == kFALSE) Phip545NegHelCM8Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip545PosHelCM8Prompt->Fill(PhpRad);
                }
            }

            if ( 580 < EGamma && EGamma < 650) {
                if(1 > CosThetapCM && CosThetapCM > 0.75 ){
                    if (BeamHelicity == kFALSE) Phip615NegHelCM1Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip615PosHelCM1Prompt->Fill(PhpRad);
                }

                else if(0.75 > CosThetapCM && CosThetapCM > 0.5){
                    if (BeamHelicity == kFALSE) Phip615NegHelCM2Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip615PosHelCM2Prompt->Fill(PhpRad);
                }

                else if(0.5 > CosThetapCM && CosThetapCM > 0.25){
                    if (BeamHelicity == kFALSE) Phip615NegHelCM3Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip615PosHelCM3Prompt->Fill(PhpRad);
                }

                else if(0.25 > CosThetapCM && CosThetapCM > 0.0){
                    if (BeamHelicity == kFALSE) Phip615NegHelCM4Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip615PosHelCM4Prompt->Fill(PhpRad);
                }

                else if(0> CosThetapCM && CosThetapCM > -0.25){
                    if (BeamHelicity == kFALSE) Phip615NegHelCM5Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip615PosHelCM5Prompt->Fill(PhpRad);
                }

                else if(-0.25 > CosThetapCM && CosThetapCM > -0.5){
                    if (BeamHelicity == kFALSE) Phip615NegHelCM6Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip615PosHelCM6Prompt->Fill(PhpRad);
                }

                else if(-0.5 > CosThetapCM && CosThetapCM > -0.75){
                    if (BeamHelicity == kFALSE) Phip615NegHelCM7Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip615PosHelCM7Prompt->Fill(PhpRad);
                }

                else if(-0.75 > CosThetapCM && CosThetapCM > -1){
                    if (BeamHelicity == kFALSE) Phip615NegHelCM8Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip615PosHelCM8Prompt->Fill(PhpRad);
                }
            }

            if ( 650 < EGamma && EGamma < 720) {
                if(1 > CosThetapCM && CosThetapCM > 0.75 ){
                    if (BeamHelicity == kFALSE) Phip685NegHelCM1Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip685PosHelCM1Prompt->Fill(PhpRad);
                }

                else if(0.75 > CosThetapCM && CosThetapCM > 0.5){
                    if (BeamHelicity == kFALSE) Phip685NegHelCM2Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip685PosHelCM2Prompt->Fill(PhpRad);
                }

                else if(0.5 > CosThetapCM && CosThetapCM > 0.25){
                    if (BeamHelicity == kFALSE) Phip685NegHelCM3Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip685PosHelCM3Prompt->Fill(PhpRad);
                }

                else if(0.25 > CosThetapCM && CosThetapCM > 0.0){
                    if (BeamHelicity == kFALSE) Phip685NegHelCM4Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip685PosHelCM4Prompt->Fill(PhpRad);
                }

                else if(0> CosThetapCM && CosThetapCM > -0.25){
                    if (BeamHelicity == kFALSE) Phip685NegHelCM5Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip685PosHelCM5Prompt->Fill(PhpRad);
                }

                else if(-0.25 > CosThetapCM && CosThetapCM > -0.5){
                    if (BeamHelicity == kFALSE) Phip685NegHelCM6Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip685PosHelCM6Prompt->Fill(PhpRad);
                }

                else if(-0.5 > CosThetapCM && CosThetapCM > -0.75){
                    if (BeamHelicity == kFALSE) Phip685NegHelCM7Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip685PosHelCM7Prompt->Fill(PhpRad);
                }

                else if(-0.75 > CosThetapCM && CosThetapCM > -1){
                    if (BeamHelicity == kFALSE) Phip685NegHelCM8Prompt->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip685PosHelCM8Prompt->Fill(PhpRad);
                }
            }
        }
    }

    if (GHistBGSub::IsRandom(TaggerTime) == kTRUE){
        EgRandom->Fill(EGamma);
        if(ThetanCM  > 80 && ThetanCM < 100){
            if ( 300 < EGamma && EGamma < 340) Phip320Random->Fill(PhpRad);
            if ( 340 < EGamma && EGamma < 380) Phip360Random->Fill(PhpRad);
            if ( 380 < EGamma && EGamma < 420) Phip400Random->Fill(PhpRad);
            if ( 420 < EGamma && EGamma < 460) Phip440Random->Fill(PhpRad);
            if ( 460 < EGamma && EGamma < 500) Phip480Random->Fill(PhpRad);
            if ( 500 < EGamma && EGamma < 540) Phip520Random->Fill(PhpRad);
            if ( 540 < EGamma && EGamma < 580) Phip560Random->Fill(PhpRad);
            if ( 580 < EGamma && EGamma < 620) Phip600Random->Fill(PhpRad);
            if ( 620 < EGamma && EGamma < 660) Phip640Random->Fill(PhpRad);
            if ( 660 < EGamma && EGamma < 700) Phip680Random->Fill(PhpRad);

            if ( 230 < EGamma && EGamma < 300) {
                if(1 > CosThetapCM && CosThetapCM > 0.75 ){
                    if (BeamHelicity == kFALSE) Phip265NegHelCM1Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip265PosHelCM1Random->Fill(PhpRad);
                }

                else if(0.75 > CosThetapCM && CosThetapCM > 0.5){
                    if (BeamHelicity == kFALSE) Phip265NegHelCM2Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip265PosHelCM2Random->Fill(PhpRad);
                }

                else if(0.5 > CosThetapCM && CosThetapCM > 0.25){
                    if (BeamHelicity == kFALSE) Phip265NegHelCM3Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip265PosHelCM3Random->Fill(PhpRad);
                }

                else if(0.25 > CosThetapCM && CosThetapCM > 0.0){
                    if (BeamHelicity == kFALSE) Phip265NegHelCM4Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip265PosHelCM4Random->Fill(PhpRad);
                }

                else if(0> CosThetapCM && CosThetapCM > -0.25){
                    if (BeamHelicity == kFALSE) Phip265NegHelCM5Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip265PosHelCM5Random->Fill(PhpRad);
                }

                else if(-0.25 > CosThetapCM && CosThetapCM > -0.5){
                    if (BeamHelicity == kFALSE) Phip265NegHelCM6Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip265PosHelCM6Random->Fill(PhpRad);
                }

                else if(-0.5 > CosThetapCM && CosThetapCM > -0.75){
                    if (BeamHelicity == kFALSE) Phip265NegHelCM7Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip265PosHelCM7Random->Fill(PhpRad);
                }

                else if(-0.75 > CosThetapCM && CosThetapCM > -1){
                    if (BeamHelicity == kFALSE) Phip265NegHelCM8Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip265PosHelCM8Random->Fill(PhpRad);
                }
            }


            if ( 300 < EGamma && EGamma < 370) {
                if(1 > CosThetapCM && CosThetapCM > 0.75 ){
                    if (BeamHelicity == kFALSE) Phip335NegHelCM1Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip335PosHelCM1Random->Fill(PhpRad);
                }

                else if(0.75 > CosThetapCM && CosThetapCM > 0.5){
                    if (BeamHelicity == kFALSE) Phip335NegHelCM2Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip335PosHelCM2Random->Fill(PhpRad);
                }

                else if(0.5 > CosThetapCM && CosThetapCM > 0.25){
                    if (BeamHelicity == kFALSE) Phip335NegHelCM3Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip335PosHelCM3Random->Fill(PhpRad);
                }

                else if(0.25 > CosThetapCM && CosThetapCM > 0.0){
                    if (BeamHelicity == kFALSE) Phip335NegHelCM4Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip335PosHelCM4Random->Fill(PhpRad);
                }

                else if(0> CosThetapCM && CosThetapCM > -0.25){
                    if (BeamHelicity == kFALSE) Phip335NegHelCM5Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip335PosHelCM5Random->Fill(PhpRad);
                }

                else if(-0.25 > CosThetapCM && CosThetapCM > -0.5){
                    if (BeamHelicity == kFALSE) Phip335NegHelCM6Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip335PosHelCM6Random->Fill(PhpRad);
                }

                else if(-0.5 > CosThetapCM && CosThetapCM > -0.75){
                    if (BeamHelicity == kFALSE) Phip335NegHelCM7Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip335PosHelCM7Random->Fill(PhpRad);
                }

                else if(-0.75 > CosThetapCM && CosThetapCM > -1){
                    if (BeamHelicity == kFALSE) Phip335NegHelCM8Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip335PosHelCM8Random->Fill(PhpRad);
                }
            }

            if ( 370 < EGamma && EGamma < 440) {
                if(1 > CosThetapCM && CosThetapCM > 0.75 ){
                    if (BeamHelicity == kFALSE) Phip405NegHelCM1Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip405PosHelCM1Random->Fill(PhpRad);
                }

                else if(0.75 > CosThetapCM && CosThetapCM > 0.5){
                    if (BeamHelicity == kFALSE) Phip405NegHelCM2Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip405PosHelCM2Random->Fill(PhpRad);
                }

                else if(0.5 > CosThetapCM && CosThetapCM > 0.25){
                    if (BeamHelicity == kFALSE) Phip405NegHelCM3Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip405PosHelCM3Random->Fill(PhpRad);
                }

                else if(0.25 > CosThetapCM && CosThetapCM > 0.0){
                    if (BeamHelicity == kFALSE) Phip405NegHelCM4Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip405PosHelCM4Random->Fill(PhpRad);
                }

                else if(0> CosThetapCM && CosThetapCM > -0.25){
                    if (BeamHelicity == kFALSE) Phip405NegHelCM5Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip405PosHelCM5Random->Fill(PhpRad);
                }

                else if(-0.25 > CosThetapCM && CosThetapCM > -0.5){
                    if (BeamHelicity == kFALSE) Phip405NegHelCM6Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip405PosHelCM6Random->Fill(PhpRad);
                }

                else if(-0.5 > CosThetapCM && CosThetapCM > -0.75){
                    if (BeamHelicity == kFALSE) Phip405NegHelCM7Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip405PosHelCM7Random->Fill(PhpRad);
                }

                else if(-0.75 > CosThetapCM && CosThetapCM > -1){
                    if (BeamHelicity == kFALSE) Phip405NegHelCM8Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip405PosHelCM8Random->Fill(PhpRad);
                }
            }

            if ( 440 < EGamma && EGamma < 510) {
                if(1 > CosThetapCM && CosThetapCM > 0.75 ){
                    if (BeamHelicity == kFALSE) Phip475NegHelCM1Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip475PosHelCM1Random->Fill(PhpRad);
                }

                else if(0.75 > CosThetapCM && CosThetapCM > 0.5){
                    if (BeamHelicity == kFALSE) Phip475NegHelCM2Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip475PosHelCM2Random->Fill(PhpRad);
                }

                else if(0.5 > CosThetapCM && CosThetapCM > 0.25){
                    if (BeamHelicity == kFALSE) Phip475NegHelCM3Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip475PosHelCM3Random->Fill(PhpRad);
                }

                else if(0.25 > CosThetapCM && CosThetapCM > 0.0){
                    if (BeamHelicity == kFALSE) Phip475NegHelCM4Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip475PosHelCM4Random->Fill(PhpRad);
                }

                else if(0> CosThetapCM && CosThetapCM > -0.25){
                    if (BeamHelicity == kFALSE) Phip475NegHelCM5Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip475PosHelCM5Random->Fill(PhpRad);
                }

                else if(-0.25 > CosThetapCM && CosThetapCM > -0.5){
                    if (BeamHelicity == kFALSE) Phip475NegHelCM6Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip475PosHelCM6Random->Fill(PhpRad);
                }

                else if(-0.5 > CosThetapCM && CosThetapCM > -0.75){
                    if (BeamHelicity == kFALSE) Phip475NegHelCM7Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip475PosHelCM7Random->Fill(PhpRad);
                }

                else if(-0.75 > CosThetapCM && CosThetapCM > -1){
                    if (BeamHelicity == kFALSE) Phip475NegHelCM8Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip475PosHelCM8Random->Fill(PhpRad);
                }
            }

            if ( 510 < EGamma && EGamma < 580) {
                if(1 > CosThetapCM && CosThetapCM > 0.75 ){
                    if (BeamHelicity == kFALSE) Phip545NegHelCM1Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip545PosHelCM1Random->Fill(PhpRad);
                }

                else if(0.75 > CosThetapCM && CosThetapCM > 0.5){
                    if (BeamHelicity == kFALSE) Phip545NegHelCM2Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip545PosHelCM2Random->Fill(PhpRad);
                }

                else if(0.5 > CosThetapCM && CosThetapCM > 0.25){
                    if (BeamHelicity == kFALSE) Phip545NegHelCM3Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip545PosHelCM3Random->Fill(PhpRad);
                }

                else if(0.25 > CosThetapCM && CosThetapCM > 0.0){
                    if (BeamHelicity == kFALSE) Phip545NegHelCM4Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip545PosHelCM4Random->Fill(PhpRad);
                }

                else if(0> CosThetapCM && CosThetapCM > -0.25){
                    if (BeamHelicity == kFALSE) Phip545NegHelCM5Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip545PosHelCM5Random->Fill(PhpRad);
                }

                else if(-0.25 > CosThetapCM && CosThetapCM > -0.5){
                    if (BeamHelicity == kFALSE) Phip545NegHelCM6Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip545PosHelCM6Random->Fill(PhpRad);
                }

                else if(-0.5 > CosThetapCM && CosThetapCM > -0.75){
                    if (BeamHelicity == kFALSE) Phip545NegHelCM7Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip545PosHelCM7Random->Fill(PhpRad);
                }

                else if(-0.75 > CosThetapCM && CosThetapCM > -1){
                    if (BeamHelicity == kFALSE) Phip545NegHelCM8Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip545PosHelCM8Random->Fill(PhpRad);
                }
            }

            if ( 580 < EGamma && EGamma < 650) {
                if(1 > CosThetapCM && CosThetapCM > 0.75 ){
                    if (BeamHelicity == kFALSE) Phip615NegHelCM1Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip615PosHelCM1Random->Fill(PhpRad);
                }

                else if(0.75 > CosThetapCM && CosThetapCM > 0.5){
                    if (BeamHelicity == kFALSE) Phip615NegHelCM2Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip615PosHelCM2Random->Fill(PhpRad);
                }

                else if(0.5 > CosThetapCM && CosThetapCM > 0.25){
                    if (BeamHelicity == kFALSE) Phip615NegHelCM3Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip615PosHelCM3Random->Fill(PhpRad);
                }

                else if(0.25 > CosThetapCM && CosThetapCM > 0.0){
                    if (BeamHelicity == kFALSE) Phip615NegHelCM4Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip615PosHelCM4Random->Fill(PhpRad);
                }

                else if(0> CosThetapCM && CosThetapCM > -0.25){
                    if (BeamHelicity == kFALSE) Phip615NegHelCM5Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip615PosHelCM5Random->Fill(PhpRad);
                }

                else if(-0.25 > CosThetapCM && CosThetapCM > -0.5){
                    if (BeamHelicity == kFALSE) Phip615NegHelCM6Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip615PosHelCM6Random->Fill(PhpRad);
                }

                else if(-0.5 > CosThetapCM && CosThetapCM > -0.75){
                    if (BeamHelicity == kFALSE) Phip615NegHelCM7Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip615PosHelCM7Random->Fill(PhpRad);
                }

                else if(-0.75 > CosThetapCM && CosThetapCM > -1){
                    if (BeamHelicity == kFALSE) Phip615NegHelCM8Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip615PosHelCM8Random->Fill(PhpRad);
                }
            }

            if ( 650 < EGamma && EGamma < 720) {
                if(1 > CosThetapCM && CosThetapCM > 0.75 ){
                    if (BeamHelicity == kFALSE) Phip685NegHelCM1Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip685PosHelCM1Random->Fill(PhpRad);
                }

                else if(0.75 > CosThetapCM && CosThetapCM > 0.5){
                    if (BeamHelicity == kFALSE) Phip685NegHelCM2Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip685PosHelCM2Random->Fill(PhpRad);
                }

                else if(0.5 > CosThetapCM && CosThetapCM > 0.25){
                    if (BeamHelicity == kFALSE) Phip685NegHelCM3Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip685PosHelCM3Random->Fill(PhpRad);
                }

                else if(0.25 > CosThetapCM && CosThetapCM > 0.0){
                    if (BeamHelicity == kFALSE) Phip685NegHelCM4Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip685PosHelCM4Random->Fill(PhpRad);
                }

                else if(0> CosThetapCM && CosThetapCM > -0.25){
                    if (BeamHelicity == kFALSE) Phip685NegHelCM5Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip685PosHelCM5Random->Fill(PhpRad);
                }

                else if(-0.25 > CosThetapCM && CosThetapCM > -0.5){
                    if (BeamHelicity == kFALSE) Phip685NegHelCM6Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip685PosHelCM6Random->Fill(PhpRad);
                }

                else if(-0.5 > CosThetapCM && CosThetapCM > -0.75){
                    if (BeamHelicity == kFALSE) Phip685NegHelCM7Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip685PosHelCM7Random->Fill(PhpRad);
                }

                else if(-0.75 > CosThetapCM && CosThetapCM > -1){
                    if (BeamHelicity == kFALSE) Phip685NegHelCM8Random->Fill(PhpRad);
                    else if (BeamHelicity == kTRUE) Phip685PosHelCM8Random->Fill(PhpRad);
                }
            }
        }
    }
}

void PNeutPol_Polarimeter_Lin_NoScatt::BGSub(){

    Eg2->Add(EgPrompt);
    Eg2->Add(EgRandom, -PvRratio);
    Phip320->Add(Phip320Prompt);
    Phip320->Add(Phip320Random, -PvRratio);
    Phip360->Add(Phip360Prompt);
    Phip360->Add(Phip360Random, -PvRratio);
    Phip400->Add(Phip400Prompt);
    Phip400->Add(Phip400Random, -PvRratio);
    Phip440->Add(Phip440Prompt);
    Phip440->Add(Phip440Random, -PvRratio);
    Phip480->Add(Phip480Prompt);
    Phip480->Add(Phip480Random, -PvRratio);
    Phip520->Add(Phip520Prompt);
    Phip520->Add(Phip520Random, -PvRratio);
    Phip560->Add(Phip560Prompt);
    Phip560->Add(Phip560Random, -PvRratio);
    Phip600->Add(Phip600Prompt);
    Phip600->Add(Phip600Random, -PvRratio);
    Phip640->Add(Phip640Prompt);
    Phip640->Add(Phip640Random, -PvRratio);
    Phip680->Add(Phip680Prompt);
    Phip680->Add(Phip680Random, -PvRratio);

    Phip265NegHelCM1->Add(Phip265NegHelCM1Prompt);
    Phip335NegHelCM1->Add(Phip335NegHelCM1Prompt);
    Phip405NegHelCM1->Add(Phip405NegHelCM1Prompt);
    Phip475NegHelCM1->Add(Phip475NegHelCM1Prompt);
    Phip545NegHelCM1->Add(Phip545NegHelCM1Prompt);
    Phip615NegHelCM1->Add(Phip615NegHelCM1Prompt);
    Phip685NegHelCM1->Add(Phip685NegHelCM1Prompt);
    Phip265NegHelCM2->Add(Phip265NegHelCM2Prompt);
    Phip335NegHelCM2->Add(Phip335NegHelCM2Prompt);
    Phip405NegHelCM2->Add(Phip405NegHelCM2Prompt);
    Phip475NegHelCM2->Add(Phip475NegHelCM2Prompt);
    Phip545NegHelCM2->Add(Phip545NegHelCM2Prompt);
    Phip615NegHelCM2->Add(Phip615NegHelCM2Prompt);
    Phip685NegHelCM2->Add(Phip685NegHelCM2Prompt);
    Phip265NegHelCM3->Add(Phip265NegHelCM3Prompt);
    Phip335NegHelCM3->Add(Phip335NegHelCM3Prompt);
    Phip405NegHelCM3->Add(Phip405NegHelCM3Prompt);
    Phip475NegHelCM3->Add(Phip475NegHelCM3Prompt);
    Phip545NegHelCM3->Add(Phip545NegHelCM3Prompt);
    Phip615NegHelCM3->Add(Phip615NegHelCM3Prompt);
    Phip685NegHelCM3->Add(Phip685NegHelCM3Prompt);
    Phip265NegHelCM4->Add(Phip265NegHelCM4Prompt);
    Phip335NegHelCM4->Add(Phip335NegHelCM4Prompt);
    Phip405NegHelCM4->Add(Phip405NegHelCM4Prompt);
    Phip475NegHelCM4->Add(Phip475NegHelCM4Prompt);
    Phip545NegHelCM4->Add(Phip545NegHelCM4Prompt);
    Phip615NegHelCM4->Add(Phip615NegHelCM4Prompt);
    Phip685NegHelCM4->Add(Phip685NegHelCM4Prompt);
    Phip265NegHelCM5->Add(Phip265NegHelCM5Prompt);
    Phip335NegHelCM5->Add(Phip335NegHelCM5Prompt);
    Phip405NegHelCM5->Add(Phip405NegHelCM5Prompt);
    Phip475NegHelCM5->Add(Phip475NegHelCM5Prompt);
    Phip545NegHelCM5->Add(Phip545NegHelCM5Prompt);
    Phip615NegHelCM5->Add(Phip615NegHelCM5Prompt);
    Phip685NegHelCM5->Add(Phip685NegHelCM5Prompt);
    Phip265NegHelCM6->Add(Phip265NegHelCM6Prompt);
    Phip335NegHelCM6->Add(Phip335NegHelCM6Prompt);
    Phip405NegHelCM6->Add(Phip405NegHelCM6Prompt);
    Phip475NegHelCM6->Add(Phip475NegHelCM6Prompt);
    Phip545NegHelCM6->Add(Phip545NegHelCM6Prompt);
    Phip615NegHelCM6->Add(Phip615NegHelCM6Prompt);
    Phip685NegHelCM6->Add(Phip685NegHelCM6Prompt);
    Phip265NegHelCM7->Add(Phip265NegHelCM7Prompt);
    Phip335NegHelCM7->Add(Phip335NegHelCM7Prompt);
    Phip405NegHelCM7->Add(Phip405NegHelCM7Prompt);
    Phip475NegHelCM7->Add(Phip475NegHelCM7Prompt);
    Phip545NegHelCM7->Add(Phip545NegHelCM7Prompt);
    Phip615NegHelCM7->Add(Phip615NegHelCM7Prompt);
    Phip685NegHelCM7->Add(Phip685NegHelCM7Prompt);
    Phip265NegHelCM8->Add(Phip265NegHelCM8Prompt);
    Phip335NegHelCM8->Add(Phip335NegHelCM8Prompt);
    Phip405NegHelCM8->Add(Phip405NegHelCM8Prompt);
    Phip475NegHelCM8->Add(Phip475NegHelCM8Prompt);
    Phip545NegHelCM8->Add(Phip545NegHelCM8Prompt);
    Phip615NegHelCM8->Add(Phip615NegHelCM8Prompt);
    Phip685NegHelCM8->Add(Phip685NegHelCM8Prompt);
    Phip265PosHelCM1->Add(Phip265PosHelCM1Prompt);
    Phip335PosHelCM1->Add(Phip335PosHelCM1Prompt);
    Phip405PosHelCM1->Add(Phip405PosHelCM1Prompt);
    Phip475PosHelCM1->Add(Phip475PosHelCM1Prompt);
    Phip545PosHelCM1->Add(Phip545PosHelCM1Prompt);
    Phip615PosHelCM1->Add(Phip615PosHelCM1Prompt);
    Phip685PosHelCM1->Add(Phip685PosHelCM1Prompt);
    Phip265PosHelCM2->Add(Phip265PosHelCM2Prompt);
    Phip335PosHelCM2->Add(Phip335PosHelCM2Prompt);
    Phip405PosHelCM2->Add(Phip405PosHelCM2Prompt);
    Phip475PosHelCM2->Add(Phip475PosHelCM2Prompt);
    Phip545PosHelCM2->Add(Phip545PosHelCM2Prompt);
    Phip615PosHelCM2->Add(Phip615PosHelCM2Prompt);
    Phip685PosHelCM2->Add(Phip685PosHelCM2Prompt);
    Phip265PosHelCM3->Add(Phip265PosHelCM3Prompt);
    Phip335PosHelCM3->Add(Phip335PosHelCM3Prompt);
    Phip405PosHelCM3->Add(Phip405PosHelCM3Prompt);
    Phip475PosHelCM3->Add(Phip475PosHelCM3Prompt);
    Phip545PosHelCM3->Add(Phip545PosHelCM3Prompt);
    Phip615PosHelCM3->Add(Phip615PosHelCM3Prompt);
    Phip685PosHelCM3->Add(Phip685PosHelCM3Prompt);
    Phip265PosHelCM4->Add(Phip265PosHelCM4Prompt);
    Phip335PosHelCM4->Add(Phip335PosHelCM4Prompt);
    Phip405PosHelCM4->Add(Phip405PosHelCM4Prompt);
    Phip475PosHelCM4->Add(Phip475PosHelCM4Prompt);
    Phip545PosHelCM4->Add(Phip545PosHelCM4Prompt);
    Phip615PosHelCM4->Add(Phip615PosHelCM4Prompt);
    Phip685PosHelCM4->Add(Phip685PosHelCM4Prompt);
    Phip265PosHelCM5->Add(Phip265PosHelCM5Prompt);
    Phip335PosHelCM5->Add(Phip335PosHelCM5Prompt);
    Phip405PosHelCM5->Add(Phip405PosHelCM5Prompt);
    Phip475PosHelCM5->Add(Phip475PosHelCM5Prompt);
    Phip545PosHelCM5->Add(Phip545PosHelCM5Prompt);
    Phip615PosHelCM5->Add(Phip615PosHelCM5Prompt);
    Phip685PosHelCM5->Add(Phip685PosHelCM5Prompt);
    Phip265PosHelCM6->Add(Phip265PosHelCM6Prompt);
    Phip335PosHelCM6->Add(Phip335PosHelCM6Prompt);
    Phip405PosHelCM6->Add(Phip405PosHelCM6Prompt);
    Phip475PosHelCM6->Add(Phip475PosHelCM6Prompt);
    Phip545PosHelCM6->Add(Phip545PosHelCM6Prompt);
    Phip615PosHelCM6->Add(Phip615PosHelCM6Prompt);
    Phip685PosHelCM6->Add(Phip685PosHelCM6Prompt);
    Phip265PosHelCM7->Add(Phip265PosHelCM7Prompt);
    Phip335PosHelCM7->Add(Phip335PosHelCM7Prompt);
    Phip405PosHelCM7->Add(Phip405PosHelCM7Prompt);
    Phip475PosHelCM7->Add(Phip475PosHelCM7Prompt);
    Phip545PosHelCM7->Add(Phip545PosHelCM7Prompt);
    Phip615PosHelCM7->Add(Phip615PosHelCM7Prompt);
    Phip685PosHelCM7->Add(Phip685PosHelCM7Prompt);
    Phip265PosHelCM8->Add(Phip265PosHelCM8Prompt);
    Phip335PosHelCM8->Add(Phip335PosHelCM8Prompt);
    Phip405PosHelCM8->Add(Phip405PosHelCM8Prompt);
    Phip475PosHelCM8->Add(Phip475PosHelCM8Prompt);
    Phip545PosHelCM8->Add(Phip545PosHelCM8Prompt);
    Phip615PosHelCM8->Add(Phip615PosHelCM8Prompt);
    Phip685PosHelCM8->Add(Phip685PosHelCM8Prompt);

    Phip265NegHelCM1->Add(Phip265NegHelCM1Random, -PvRratio);
    Phip335NegHelCM1->Add(Phip335NegHelCM1Random, -PvRratio);
    Phip405NegHelCM1->Add(Phip405NegHelCM1Random, -PvRratio);
    Phip475NegHelCM1->Add(Phip475NegHelCM1Random, -PvRratio);
    Phip545NegHelCM1->Add(Phip545NegHelCM1Random, -PvRratio);
    Phip615NegHelCM1->Add(Phip615NegHelCM1Random, -PvRratio);
    Phip685NegHelCM1->Add(Phip685NegHelCM1Random, -PvRratio);
    Phip265NegHelCM2->Add(Phip265NegHelCM2Random, -PvRratio);
    Phip335NegHelCM2->Add(Phip335NegHelCM2Random, -PvRratio);
    Phip405NegHelCM2->Add(Phip405NegHelCM2Random, -PvRratio);
    Phip475NegHelCM2->Add(Phip475NegHelCM2Random, -PvRratio);
    Phip545NegHelCM2->Add(Phip545NegHelCM2Random, -PvRratio);
    Phip615NegHelCM2->Add(Phip615NegHelCM2Random, -PvRratio);
    Phip685NegHelCM2->Add(Phip685NegHelCM2Random, -PvRratio);
    Phip265NegHelCM3->Add(Phip265NegHelCM3Random, -PvRratio);
    Phip335NegHelCM3->Add(Phip335NegHelCM3Random, -PvRratio);
    Phip405NegHelCM3->Add(Phip405NegHelCM3Random, -PvRratio);
    Phip475NegHelCM3->Add(Phip475NegHelCM3Random, -PvRratio);
    Phip545NegHelCM3->Add(Phip545NegHelCM3Random, -PvRratio);
    Phip615NegHelCM3->Add(Phip615NegHelCM3Random, -PvRratio);
    Phip685NegHelCM3->Add(Phip685NegHelCM3Random, -PvRratio);
    Phip265NegHelCM4->Add(Phip265NegHelCM4Random, -PvRratio);
    Phip335NegHelCM4->Add(Phip335NegHelCM4Random, -PvRratio);
    Phip405NegHelCM4->Add(Phip405NegHelCM4Random, -PvRratio);
    Phip475NegHelCM4->Add(Phip475NegHelCM4Random, -PvRratio);
    Phip545NegHelCM4->Add(Phip545NegHelCM4Random, -PvRratio);
    Phip615NegHelCM4->Add(Phip615NegHelCM4Random, -PvRratio);
    Phip685NegHelCM4->Add(Phip685NegHelCM4Random, -PvRratio);
    Phip265NegHelCM5->Add(Phip265NegHelCM5Random, -PvRratio);
    Phip335NegHelCM5->Add(Phip335NegHelCM5Random, -PvRratio);
    Phip405NegHelCM5->Add(Phip405NegHelCM5Random, -PvRratio);
    Phip475NegHelCM5->Add(Phip475NegHelCM5Random, -PvRratio);
    Phip545NegHelCM5->Add(Phip545NegHelCM5Random, -PvRratio);
    Phip615NegHelCM5->Add(Phip615NegHelCM5Random, -PvRratio);
    Phip685NegHelCM5->Add(Phip685NegHelCM5Random, -PvRratio);
    Phip265NegHelCM6->Add(Phip265NegHelCM6Random, -PvRratio);
    Phip335NegHelCM6->Add(Phip335NegHelCM6Random, -PvRratio);
    Phip405NegHelCM6->Add(Phip405NegHelCM6Random, -PvRratio);
    Phip475NegHelCM6->Add(Phip475NegHelCM6Random, -PvRratio);
    Phip545NegHelCM6->Add(Phip545NegHelCM6Random, -PvRratio);
    Phip615NegHelCM6->Add(Phip615NegHelCM6Random, -PvRratio);
    Phip685NegHelCM6->Add(Phip685NegHelCM6Random, -PvRratio);
    Phip265NegHelCM7->Add(Phip265NegHelCM7Random, -PvRratio);
    Phip335NegHelCM7->Add(Phip335NegHelCM7Random, -PvRratio);
    Phip405NegHelCM7->Add(Phip405NegHelCM7Random, -PvRratio);
    Phip475NegHelCM7->Add(Phip475NegHelCM7Random, -PvRratio);
    Phip545NegHelCM7->Add(Phip545NegHelCM7Random, -PvRratio);
    Phip615NegHelCM7->Add(Phip615NegHelCM7Random, -PvRratio);
    Phip685NegHelCM7->Add(Phip685NegHelCM7Random, -PvRratio);
    Phip265NegHelCM8->Add(Phip265NegHelCM8Random, -PvRratio);
    Phip335NegHelCM8->Add(Phip335NegHelCM8Random, -PvRratio);
    Phip405NegHelCM8->Add(Phip405NegHelCM8Random, -PvRratio);
    Phip475NegHelCM8->Add(Phip475NegHelCM8Random, -PvRratio);
    Phip545NegHelCM8->Add(Phip545NegHelCM8Random, -PvRratio);
    Phip615NegHelCM8->Add(Phip615NegHelCM8Random, -PvRratio);
    Phip685NegHelCM8->Add(Phip685NegHelCM8Random, -PvRratio);
    Phip265PosHelCM1->Add(Phip265PosHelCM1Random, -PvRratio);
    Phip335PosHelCM1->Add(Phip335PosHelCM1Random, -PvRratio);
    Phip405PosHelCM1->Add(Phip405PosHelCM1Random, -PvRratio);
    Phip475PosHelCM1->Add(Phip475PosHelCM1Random, -PvRratio);
    Phip545PosHelCM1->Add(Phip545PosHelCM1Random, -PvRratio);
    Phip615PosHelCM1->Add(Phip615PosHelCM1Random, -PvRratio);
    Phip685PosHelCM1->Add(Phip685PosHelCM1Random, -PvRratio);
    Phip265PosHelCM2->Add(Phip265PosHelCM2Random, -PvRratio);
    Phip335PosHelCM2->Add(Phip335PosHelCM2Random, -PvRratio);
    Phip405PosHelCM2->Add(Phip405PosHelCM2Random, -PvRratio);
    Phip475PosHelCM2->Add(Phip475PosHelCM2Random, -PvRratio);
    Phip545PosHelCM2->Add(Phip545PosHelCM2Random, -PvRratio);
    Phip615PosHelCM2->Add(Phip615PosHelCM2Random, -PvRratio);
    Phip685PosHelCM2->Add(Phip685PosHelCM2Random, -PvRratio);
    Phip265PosHelCM3->Add(Phip265PosHelCM3Random, -PvRratio);
    Phip335PosHelCM3->Add(Phip335PosHelCM3Random, -PvRratio);
    Phip405PosHelCM3->Add(Phip405PosHelCM3Random, -PvRratio);
    Phip475PosHelCM3->Add(Phip475PosHelCM3Random, -PvRratio);
    Phip545PosHelCM3->Add(Phip545PosHelCM3Random, -PvRratio);
    Phip615PosHelCM3->Add(Phip615PosHelCM3Random, -PvRratio);
    Phip685PosHelCM3->Add(Phip685PosHelCM3Random, -PvRratio);
    Phip265PosHelCM4->Add(Phip265PosHelCM4Random, -PvRratio);
    Phip335PosHelCM4->Add(Phip335PosHelCM4Random, -PvRratio);
    Phip405PosHelCM4->Add(Phip405PosHelCM4Random, -PvRratio);
    Phip475PosHelCM4->Add(Phip475PosHelCM4Random, -PvRratio);
    Phip545PosHelCM4->Add(Phip545PosHelCM4Random, -PvRratio);
    Phip615PosHelCM4->Add(Phip615PosHelCM4Random, -PvRratio);
    Phip685PosHelCM4->Add(Phip685PosHelCM4Random, -PvRratio);
    Phip265PosHelCM5->Add(Phip265PosHelCM5Random, -PvRratio);
    Phip335PosHelCM5->Add(Phip335PosHelCM5Random, -PvRratio);
    Phip405PosHelCM5->Add(Phip405PosHelCM5Random, -PvRratio);
    Phip475PosHelCM5->Add(Phip475PosHelCM5Random, -PvRratio);
    Phip545PosHelCM5->Add(Phip545PosHelCM5Random, -PvRratio);
    Phip615PosHelCM5->Add(Phip615PosHelCM5Random, -PvRratio);
    Phip685PosHelCM5->Add(Phip685PosHelCM5Random, -PvRratio);
    Phip265PosHelCM6->Add(Phip265PosHelCM6Random, -PvRratio);
    Phip335PosHelCM6->Add(Phip335PosHelCM6Random, -PvRratio);
    Phip405PosHelCM6->Add(Phip405PosHelCM6Random, -PvRratio);
    Phip475PosHelCM6->Add(Phip475PosHelCM6Random, -PvRratio);
    Phip545PosHelCM6->Add(Phip545PosHelCM6Random, -PvRratio);
    Phip615PosHelCM6->Add(Phip615PosHelCM6Random, -PvRratio);
    Phip685PosHelCM6->Add(Phip685PosHelCM6Random, -PvRratio);
    Phip265PosHelCM7->Add(Phip265PosHelCM7Random, -PvRratio);
    Phip335PosHelCM7->Add(Phip335PosHelCM7Random, -PvRratio);
    Phip405PosHelCM7->Add(Phip405PosHelCM7Random, -PvRratio);
    Phip475PosHelCM7->Add(Phip475PosHelCM7Random, -PvRratio);
    Phip545PosHelCM7->Add(Phip545PosHelCM7Random, -PvRratio);
    Phip615PosHelCM7->Add(Phip615PosHelCM7Random, -PvRratio);
    Phip685PosHelCM7->Add(Phip685PosHelCM7Random, -PvRratio);
    Phip265PosHelCM8->Add(Phip265PosHelCM8Random, -PvRratio);
    Phip335PosHelCM8->Add(Phip335PosHelCM8Random, -PvRratio);
    Phip405PosHelCM8->Add(Phip405PosHelCM8Random, -PvRratio);
    Phip475PosHelCM8->Add(Phip475PosHelCM8Random, -PvRratio);
    Phip545PosHelCM8->Add(Phip545PosHelCM8Random, -PvRratio);
    Phip615PosHelCM8->Add(Phip615PosHelCM8Random, -PvRratio);
    Phip685PosHelCM8->Add(Phip685PosHelCM8Random, -PvRratio);

}

Bool_t	PNeutPol_Polarimeter_Lin_NoScatt::Write(){
  // Write all GH1's easily

  GTreeManager::Write();
}
