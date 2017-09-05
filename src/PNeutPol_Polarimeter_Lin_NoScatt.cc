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
        GVpCorrB = GVp;
        GVpCorrB.Boost(b); // Boost GVp to CM frame
        ThetapCM = (GVpCorrB.Theta())*TMath::RadToDeg(); // Get Theta of proton in CM frame
        CosThetapCM = cos (GVpCorrB.Theta());

        // Gamma(d,p)n
        KinEp = CalcKinEnergy(Thp, EGamma, Md, 0., Mp, Mn); // Calculate kin E of proton assuming pn production
        RecKinProton = LProton4VectorKin(KinEp, ThpRad, PhpRad);
        RecKinNeutron = LNeutron4VectorKin(RecKinProton);
        ThetanRec = (RecKinNeutron.Theta()) * TMath::RadToDeg();
        PhinRec = (RecKinNeutron.Phi()) * TMath::RadToDeg();
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

        beamF.SetFromVector(Gamma); // Set Lorentz vectors for use in APLCON
        protonF.SetFromVector(GVpCorr);
        neutronF.SetFromVector(GVnCorr);

        neutronF.Ek_Sigma=0; // Set errors on Ek, Theta and Phi for n/p/Photon
        neutronF.Theta_Sigma=0.0474839+0.00425626*GVnCorr.Theta();
        neutronF.Phi_Sigma=0.112339-0.0761341*GVnCorr.Theta()+0.0244866*GVnCorr.Theta()*GVnCorr.Theta();

        protonF.Ek_Sigma=(0.045+0.043*(GVpCorr.E()-GVpCorr.M()))*4;
        protonF.Theta_Sigma=0.00920133-0.00511389*GVpCorr.Theta()+0.00307915*GVpCorr.Theta()*GVpCorr.Theta();
        protonF.Phi_Sigma=0.00974036+0.00411955*GVpCorr.Theta()-0.0096472*GVpCorr.Theta()*GVpCorr.Theta()+0.00414428*GVpCorr.Theta()*GVpCorr.Theta()*GVpCorr.Theta();
        if(EpCorr>410)protonF.Ek_Sigma=0;

        beamF.Ek_Sigma=0.87-0.0000011*(EGamma-177);
        beamF.Theta_Sigma=1e-3;
        beamF.Phi_Sigma=1e-3;

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

//    // Proton Phi dists across EGamma bins
//    Phip425CM1 = new GH1("Phip_425MeVCM1", "#phi_{p} Distribution for E_{#gamma} 425 #pm 5MeV (Cos#theta_{CM}1-0.8)", 10, -180, 180);
//    Phip435CM1 = new GH1("Phip_435MeVCM1", "#phi_{p} Distribution for E_{#gamma} 435 #pm 5MeV (Cos#theta_{CM}1-0.8)", 10, -180, 180);
//    Phip445CM1 = new GH1("Phip_445MeVCM1", "#phi_{p} Distribution for E_{#gamma} 445 #pm 5MeV (Cos#theta_{CM}1-0.8)", 10, -180, 180);
//    Phip455CM1 = new GH1("Phip_455MeVCM1", "#phi_{p} Distribution for E_{#gamma} 455 #pm 5MeV (Cos#theta_{CM}1-0.8)", 10, -180, 180);
//    Phip465CM1 = new GH1("Phip_465MeVCM1", "#phi_{p} Distribution for E_{#gamma} 465 #pm 5MeV (Cos#theta_{CM}1-0.8)", 10, -180, 180);
//    Phip475CM1 = new GH1("Phip_475MeVCM1", "#phi_{p} Distribution for E_{#gamma} 475 #pm 5MeV (Cos#theta_{CM}1-0.8)", 10, -180, 180);
//    Phip485CM1 = new GH1("Phip_485MeVCM1", "#phi_{p} Distribution for E_{#gamma} 485 #pm 5MeV (Cos#theta_{CM}1-0.8)", 10, -180, 180);
//    Phip495CM1 = new GH1("Phip_495MeVCM1", "#phi_{p} Distribution for E_{#gamma} 495 #pm 5MeV (Cos#theta_{CM}1-0.8)", 10, -180, 180);
//    Phip505CM1 = new GH1("Phip_505MeVCM1", "#phi_{p} Distribution for E_{#gamma} 505 #pm 5MeV (Cos#theta_{CM}1-0.8)", 10, -180, 180);
//    Phip515CM1 = new GH1("Phip_515MeVCM1", "#phi_{p} Distribution for E_{#gamma} 515 #pm 5MeV (Cos#theta_{CM}1-0.8)", 10, -180, 180);
//    Phip525CM1 = new GH1("Phip_525MeVCM1", "#phi_{p} Distribution for E_{#gamma} 525 #pm 5MeV (Cos#theta_{CM}1-0.8)", 10, -180, 180);
//    Phip535CM1 = new GH1("Phip_535MeVCM1", "#phi_{p} Distribution for E_{#gamma} 535 #pm 5MeV (Cos#theta_{CM}1-0.8)", 10, -180, 180);
//    Phip545CM1 = new GH1("Phip_545MeVCM1", "#phi_{p} Distribution for E_{#gamma} 545 #pm 5MeV (Cos#theta_{CM}1-0.8)", 10, -180, 180);
//    Phip555CM1 = new GH1("Phip_555MeVCM1", "#phi_{p} Distribution for E_{#gamma} 555 #pm 5MeV (Cos#theta_{CM}1-0.8)", 10, -180, 180);
//    Phip565CM1 = new GH1("Phip_565MeVCM1", "#phi_{p} Distribution for E_{#gamma} 565 #pm 5MeV (Cos#theta_{CM}1-0.8)", 10, -180, 180);
//    Phip575CM1 = new GH1("Phip_575MeVCM1", "#phi_{p} Distribution for E_{#gamma} 575 #pm 5MeV (Cos#theta_{CM}1-0.8)", 10, -180, 180);
//    Phip585CM1 = new GH1("Phip_585MeVCM1", "#phi_{p} Distribution for E_{#gamma} 585 #pm 5MeV (Cos#theta_{CM}1-0.8)", 10, -180, 180);
//    Phip595CM1 = new GH1("Phip_595MeVCM1", "#phi_{p} Distribution for E_{#gamma} 595 #pm 5MeV (Cos#theta_{CM}1-0.8)", 10, -180, 180);
//    Phip605CM1 = new GH1("Phip_605MeVCM1", "#phi_{p} Distribution for E_{#gamma} 605 #pm 5MeV (Cos#theta_{CM}1-0.8)", 10, -180, 180);
//    Phip615CM1 = new GH1("Phip_615MeVCM1", "#phi_{p} Distribution for E_{#gamma} 615 #pm 5MeV (Cos#theta_{CM}1-0.8)", 10, -180, 180);
//
//    // #phi_{p} dists across EGamma bins
//    Phip425CM2 = new GH1("Phip_425MeVCM2", "#phi_{p} Distribution for E_{#gamma} 425 #pm 5MeV (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
//    Phip435CM2 = new GH1("Phip_435MeVCM2", "#phi_{p} Distribution for E_{#gamma} 435 #pm 5MeV (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
//    Phip445CM2 = new GH1("Phip_445MeVCM2", "#phi_{p} Distribution for E_{#gamma} 445 #pm 5MeV (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
//    Phip455CM2 = new GH1("Phip_455MeVCM2", "#phi_{p} Distribution for E_{#gamma} 455 #pm 5MeV (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
//    Phip465CM2 = new GH1("Phip_465MeVCM2", "#phi_{p} Distribution for E_{#gamma} 465 #pm 5MeV (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
//    Phip475CM2 = new GH1("Phip_475MeVCM2", "#phi_{p} Distribution for E_{#gamma} 475 #pm 5MeV (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
//    Phip485CM2 = new GH1("Phip_485MeVCM2", "#phi_{p} Distribution for E_{#gamma} 485 #pm 5MeV (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
//    Phip495CM2 = new GH1("Phip_495MeVCM2", "#phi_{p} Distribution for E_{#gamma} 495 #pm 5MeV (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
//    Phip505CM2 = new GH1("Phip_505MeVCM2", "#phi_{p} Distribution for E_{#gamma} 505 #pm 5MeV (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
//    Phip515CM2 = new GH1("Phip_515MeVCM2", "#phi_{p} Distribution for E_{#gamma} 515 #pm 5MeV (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
//    Phip525CM2 = new GH1("Phip_525MeVCM2", "#phi_{p} Distribution for E_{#gamma} 525 #pm 5MeV (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
//    Phip535CM2 = new GH1("Phip_535MeVCM2", "#phi_{p} Distribution for E_{#gamma} 535 #pm 5MeV (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
//    Phip545CM2 = new GH1("Phip_545MeVCM2", "#phi_{p} Distribution for E_{#gamma} 545 #pm 5MeV (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
//    Phip555CM2 = new GH1("Phip_555MeVCM2", "#phi_{p} Distribution for E_{#gamma} 555 #pm 5MeV (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
//    Phip565CM2 = new GH1("Phip_565MeVCM2", "#phi_{p} Distribution for E_{#gamma} 565 #pm 5MeV (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
//    Phip575CM2 = new GH1("Phip_575MeVCM2", "#phi_{p} Distribution for E_{#gamma} 575 #pm 5MeV (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
//    Phip585CM2 = new GH1("Phip_585MeVCM2", "#phi_{p} Distribution for E_{#gamma} 585 #pm 5MeV (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
//    Phip595CM2 = new GH1("Phip_595MeVCM2", "#phi_{p} Distribution for E_{#gamma} 595 #pm 5MeV (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
//    Phip605CM2 = new GH1("Phip_605MeVCM2", "#phi_{p} Distribution for E_{#gamma} 605 #pm 5MeV (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
//    Phip615CM2 = new GH1("Phip_615MeVCM2", "#phi_{p} Distribution for E_{#gamma} 615 #pm 5MeV (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
//
//    // #phi_{p} dists across EGamma bins
//    Phip425CM3 = new GH1("Phip_425MeVCM3", "#phi_{p} Distribution for E_{#gamma} 425 #pm 5MeV (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
//    Phip435CM3 = new GH1("Phip_435MeVCM3", "#phi_{p} Distribution for E_{#gamma} 435 #pm 5MeV (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
//    Phip445CM3 = new GH1("Phip_445MeVCM3", "#phi_{p} Distribution for E_{#gamma} 445 #pm 5MeV (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
//    Phip455CM3 = new GH1("Phip_455MeVCM3", "#phi_{p} Distribution for E_{#gamma} 455 #pm 5MeV (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
//    Phip465CM3 = new GH1("Phip_465MeVCM3", "#phi_{p} Distribution for E_{#gamma} 465 #pm 5MeV (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
//    Phip475CM3 = new GH1("Phip_475MeVCM3", "#phi_{p} Distribution for E_{#gamma} 475 #pm 5MeV (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
//    Phip485CM3 = new GH1("Phip_485MeVCM3", "#phi_{p} Distribution for E_{#gamma} 485 #pm 5MeV (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
//    Phip495CM3 = new GH1("Phip_495MeVCM3", "#phi_{p} Distribution for E_{#gamma} 495 #pm 5MeV (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
//    Phip505CM3 = new GH1("Phip_505MeVCM3", "#phi_{p} Distribution for E_{#gamma} 505 #pm 5MeV (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
//    Phip515CM3 = new GH1("Phip_515MeVCM3", "#phi_{p} Distribution for E_{#gamma} 515 #pm 5MeV (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
//    Phip525CM3 = new GH1("Phip_525MeVCM3", "#phi_{p} Distribution for E_{#gamma} 525 #pm 5MeV (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
//    Phip535CM3 = new GH1("Phip_535MeVCM3", "#phi_{p} Distribution for E_{#gamma} 535 #pm 5MeV (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
//    Phip545CM3 = new GH1("Phip_545MeVCM3", "#phi_{p} Distribution for E_{#gamma} 545 #pm 5MeV (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
//    Phip555CM3 = new GH1("Phip_555MeVCM3", "#phi_{p} Distribution for E_{#gamma} 555 #pm 5MeV (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
//    Phip565CM3 = new GH1("Phip_565MeVCM3", "#phi_{p} Distribution for E_{#gamma} 565 #pm 5MeV (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
//    Phip575CM3 = new GH1("Phip_575MeVCM3", "#phi_{p} Distribution for E_{#gamma} 575 #pm 5MeV (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
//    Phip585CM3 = new GH1("Phip_585MeVCM3", "#phi_{p} Distribution for E_{#gamma} 585 #pm 5MeV (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
//    Phip595CM3 = new GH1("Phip_595MeVCM3", "#phi_{p} Distribution for E_{#gamma} 595 #pm 5MeV (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
//    Phip605CM3 = new GH1("Phip_605MeVCM3", "#phi_{p} Distribution for E_{#gamma} 605 #pm 5MeV (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
//    Phip615CM3 = new GH1("Phip_615MeVCM3", "#phi_{p} Distribution for E_{#gamma} 615 #pm 5MeV (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
//
//    // #phi_{p} dists across EGamma bins
//    Phip425CM4 = new GH1("Phip_425MeVCM4", "#phi_{p} Distribution for E_{#gamma} 425 #pm 5MeV (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
//    Phip435CM4 = new GH1("Phip_435MeVCM4", "#phi_{p} Distribution for E_{#gamma} 435 #pm 5MeV (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
//    Phip445CM4 = new GH1("Phip_445MeVCM4", "#phi_{p} Distribution for E_{#gamma} 445 #pm 5MeV (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
//    Phip455CM4 = new GH1("Phip_455MeVCM4", "#phi_{p} Distribution for E_{#gamma} 455 #pm 5MeV (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
//    Phip465CM4 = new GH1("Phip_465MeVCM4", "#phi_{p} Distribution for E_{#gamma} 465 #pm 5MeV (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
//    Phip475CM4 = new GH1("Phip_475MeVCM4", "#phi_{p} Distribution for E_{#gamma} 475 #pm 5MeV (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
//    Phip485CM4 = new GH1("Phip_485MeVCM4", "#phi_{p} Distribution for E_{#gamma} 485 #pm 5MeV (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
//    Phip495CM4 = new GH1("Phip_495MeVCM4", "#phi_{p} Distribution for E_{#gamma} 495 #pm 5MeV (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
//    Phip505CM4 = new GH1("Phip_505MeVCM4", "#phi_{p} Distribution for E_{#gamma} 505 #pm 5MeV (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
//    Phip515CM4 = new GH1("Phip_515MeVCM4", "#phi_{p} Distribution for E_{#gamma} 515 #pm 5MeV (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
//    Phip525CM4 = new GH1("Phip_525MeVCM4", "#phi_{p} Distribution for E_{#gamma} 525 #pm 5MeV (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
//    Phip535CM4 = new GH1("Phip_535MeVCM4", "#phi_{p} Distribution for E_{#gamma} 535 #pm 5MeV (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
//    Phip545CM4 = new GH1("Phip_545MeVCM4", "#phi_{p} Distribution for E_{#gamma} 545 #pm 5MeV (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
//    Phip555CM4 = new GH1("Phip_555MeVCM4", "#phi_{p} Distribution for E_{#gamma} 555 #pm 5MeV (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
//    Phip565CM4 = new GH1("Phip_565MeVCM4", "#phi_{p} Distribution for E_{#gamma} 565 #pm 5MeV (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
//    Phip575CM4 = new GH1("Phip_575MeVCM4", "#phi_{p} Distribution for E_{#gamma} 575 #pm 5MeV (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
//    Phip585CM4 = new GH1("Phip_585MeVCM4", "#phi_{p} Distribution for E_{#gamma} 585 #pm 5MeV (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
//    Phip595CM4 = new GH1("Phip_595MeVCM4", "#phi_{p} Distribution for E_{#gamma} 595 #pm 5MeV (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
//    Phip605CM4 = new GH1("Phip_605MeVCM4", "#phi_{p} Distribution for E_{#gamma} 605 #pm 5MeV (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
//    Phip615CM4 = new GH1("Phip_615MeVCM4", "#phi_{p} Distribution for E_{#gamma} 615 #pm 5MeV (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
//
//    // #phi_{p} dists across EGamma bins
//    Phip425CM5 = new GH1("Phip_425MeVCM5", "#phi_{p} Distribution for E_{#gamma} 425 #pm 5MeV (Cos#theta_{CM}0.2-0)", 10, -180, 180);
//    Phip435CM5 = new GH1("Phip_435MeVCM5", "#phi_{p} Distribution for E_{#gamma} 435 #pm 5MeV (Cos#theta_{CM}0.2-0)", 10, -180, 180);
//    Phip445CM5 = new GH1("Phip_445MeVCM5", "#phi_{p} Distribution for E_{#gamma} 445 #pm 5MeV (Cos#theta_{CM}0.2-0)", 10, -180, 180);
//    Phip455CM5 = new GH1("Phip_455MeVCM5", "#phi_{p} Distribution for E_{#gamma} 455 #pm 5MeV (Cos#theta_{CM}0.2-0)", 10, -180, 180);
//    Phip465CM5 = new GH1("Phip_465MeVCM5", "#phi_{p} Distribution for E_{#gamma} 465 #pm 5MeV (Cos#theta_{CM}0.2-0)", 10, -180, 180);
//    Phip475CM5 = new GH1("Phip_475MeVCM5", "#phi_{p} Distribution for E_{#gamma} 475 #pm 5MeV (Cos#theta_{CM}0.2-0)", 10, -180, 180);
//    Phip485CM5 = new GH1("Phip_485MeVCM5", "#phi_{p} Distribution for E_{#gamma} 485 #pm 5MeV (Cos#theta_{CM}0.2-0)", 10, -180, 180);
//    Phip495CM5 = new GH1("Phip_495MeVCM5", "#phi_{p} Distribution for E_{#gamma} 495 #pm 5MeV (Cos#theta_{CM}0.2-0)", 10, -180, 180);
//    Phip505CM5 = new GH1("Phip_505MeVCM5", "#phi_{p} Distribution for E_{#gamma} 505 #pm 5MeV (Cos#theta_{CM}0.2-0)", 10, -180, 180);
//    Phip515CM5 = new GH1("Phip_515MeVCM5", "#phi_{p} Distribution for E_{#gamma} 515 #pm 5MeV (Cos#theta_{CM}0.2-0)", 10, -180, 180);
//    Phip525CM5 = new GH1("Phip_525MeVCM5", "#phi_{p} Distribution for E_{#gamma} 525 #pm 5MeV (Cos#theta_{CM}0.2-0)", 10, -180, 180);
//    Phip535CM5 = new GH1("Phip_535MeVCM5", "#phi_{p} Distribution for E_{#gamma} 535 #pm 5MeV (Cos#theta_{CM}0.2-0)", 10, -180, 180);
//    Phip545CM5 = new GH1("Phip_545MeVCM5", "#phi_{p} Distribution for E_{#gamma} 545 #pm 5MeV (Cos#theta_{CM}0.2-0)", 10, -180, 180);
//    Phip555CM5 = new GH1("Phip_555MeVCM5", "#phi_{p} Distribution for E_{#gamma} 555 #pm 5MeV (Cos#theta_{CM}0.2-0)", 10, -180, 180);
//    Phip565CM5 = new GH1("Phip_565MeVCM5", "#phi_{p} Distribution for E_{#gamma} 565 #pm 5MeV (Cos#theta_{CM}0.2-0)", 10, -180, 180);
//    Phip575CM5 = new GH1("Phip_575MeVCM5", "#phi_{p} Distribution for E_{#gamma} 575 #pm 5MeV (Cos#theta_{CM}0.2-0)", 10, -180, 180);
//    Phip585CM5 = new GH1("Phip_585MeVCM5", "#phi_{p} Distribution for E_{#gamma} 585 #pm 5MeV (Cos#theta_{CM}0.2-0)", 10, -180, 180);
//    Phip595CM5 = new GH1("Phip_595MeVCM5", "#phi_{p} Distribution for E_{#gamma} 595 #pm 5MeV (Cos#theta_{CM}0.2-0)", 10, -180, 180);
//    Phip605CM5 = new GH1("Phip_605MeVCM5", "#phi_{p} Distribution for E_{#gamma} 605 #pm 5MeV (Cos#theta_{CM}0.2-0)", 10, -180, 180);
//    Phip615CM5 = new GH1("Phip_615MeVCM5", "#phi_{p} Distribution for E_{#gamma} 615 #pm 5MeV (Cos#theta_{CM}0.2-0)", 10, -180, 180);
//
//    // #phi_{p} dists across EGamma bins
//    Phip425CM6 = new GH1("Phip_425MeVCM6", "#phi_{p} Distribution for E_{#gamma} 425 #pm 5MeV (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
//    Phip435CM6 = new GH1("Phip_435MeVCM6", "#phi_{p} Distribution for E_{#gamma} 435 #pm 5MeV (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
//    Phip445CM6 = new GH1("Phip_445MeVCM6", "#phi_{p} Distribution for E_{#gamma} 445 #pm 5MeV (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
//    Phip455CM6 = new GH1("Phip_455MeVCM6", "#phi_{p} Distribution for E_{#gamma} 455 #pm 5MeV (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
//    Phip465CM6 = new GH1("Phip_465MeVCM6", "#phi_{p} Distribution for E_{#gamma} 465 #pm 5MeV (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
//    Phip475CM6 = new GH1("Phip_475MeVCM6", "#phi_{p} Distribution for E_{#gamma} 475 #pm 5MeV (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
//    Phip485CM6 = new GH1("Phip_485MeVCM6", "#phi_{p} Distribution for E_{#gamma} 485 #pm 5MeV (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
//    Phip495CM6 = new GH1("Phip_495MeVCM6", "#phi_{p} Distribution for E_{#gamma} 495 #pm 5MeV (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
//    Phip505CM6 = new GH1("Phip_505MeVCM6", "#phi_{p} Distribution for E_{#gamma} 505 #pm 5MeV (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
//    Phip515CM6 = new GH1("Phip_515MeVCM6", "#phi_{p} Distribution for E_{#gamma} 515 #pm 5MeV (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
//    Phip525CM6 = new GH1("Phip_525MeVCM6", "#phi_{p} Distribution for E_{#gamma} 525 #pm 5MeV (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
//    Phip535CM6 = new GH1("Phip_535MeVCM6", "#phi_{p} Distribution for E_{#gamma} 535 #pm 5MeV (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
//    Phip545CM6 = new GH1("Phip_545MeVCM6", "#phi_{p} Distribution for E_{#gamma} 545 #pm 5MeV (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
//    Phip555CM6 = new GH1("Phip_555MeVCM6", "#phi_{p} Distribution for E_{#gamma} 555 #pm 5MeV (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
//    Phip565CM6 = new GH1("Phip_565MeVCM6", "#phi_{p} Distribution for E_{#gamma} 565 #pm 5MeV (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
//    Phip575CM6 = new GH1("Phip_575MeVCM6", "#phi_{p} Distribution for E_{#gamma} 575 #pm 5MeV (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
//    Phip585CM6 = new GH1("Phip_585MeVCM6", "#phi_{p} Distribution for E_{#gamma} 585 #pm 5MeV (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
//    Phip595CM6 = new GH1("Phip_595MeVCM6", "#phi_{p} Distribution for E_{#gamma} 595 #pm 5MeV (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
//    Phip605CM6 = new GH1("Phip_605MeVCM6", "#phi_{p} Distribution for E_{#gamma} 605 #pm 5MeV (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
//    Phip615CM6 = new GH1("Phip_615MeVCM6", "#phi_{p} Distribution for E_{#gamma} 615 #pm 5MeV (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
//
//    // #phi_{p} dists across EGamma bins
//    Phip425CM7 = new GH1("Phip_425MeVCM7", "#phi_{p} Distribution for E_{#gamma} 425 #pm 5MeV (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
//    Phip435CM7 = new GH1("Phip_435MeVCM7", "#phi_{p} Distribution for E_{#gamma} 435 #pm 5MeV (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
//    Phip445CM7 = new GH1("Phip_445MeVCM7", "#phi_{p} Distribution for E_{#gamma} 445 #pm 5MeV (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
//    Phip455CM7 = new GH1("Phip_455MeVCM7", "#phi_{p} Distribution for E_{#gamma} 455 #pm 5MeV (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
//    Phip465CM7 = new GH1("Phip_465MeVCM7", "#phi_{p} Distribution for E_{#gamma} 465 #pm 5MeV (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
//    Phip475CM7 = new GH1("Phip_475MeVCM7", "#phi_{p} Distribution for E_{#gamma} 475 #pm 5MeV (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
//    Phip485CM7 = new GH1("Phip_485MeVCM7", "#phi_{p} Distribution for E_{#gamma} 485 #pm 5MeV (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
//    Phip495CM7 = new GH1("Phip_495MeVCM7", "#phi_{p} Distribution for E_{#gamma} 495 #pm 5MeV (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
//    Phip505CM7 = new GH1("Phip_505MeVCM7", "#phi_{p} Distribution for E_{#gamma} 505 #pm 5MeV (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
//    Phip515CM7 = new GH1("Phip_515MeVCM7", "#phi_{p} Distribution for E_{#gamma} 515 #pm 5MeV (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
//    Phip525CM7 = new GH1("Phip_525MeVCM7", "#phi_{p} Distribution for E_{#gamma} 525 #pm 5MeV (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
//    Phip535CM7 = new GH1("Phip_535MeVCM7", "#phi_{p} Distribution for E_{#gamma} 535 #pm 5MeV (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
//    Phip545CM7 = new GH1("Phip_545MeVCM7", "#phi_{p} Distribution for E_{#gamma} 545 #pm 5MeV (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
//    Phip555CM7 = new GH1("Phip_555MeVCM7", "#phi_{p} Distribution for E_{#gamma} 555 #pm 5MeV (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
//    Phip565CM7 = new GH1("Phip_565MeVCM7", "#phi_{p} Distribution for E_{#gamma} 565 #pm 5MeV (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
//    Phip575CM7 = new GH1("Phip_575MeVCM7", "#phi_{p} Distribution for E_{#gamma} 575 #pm 5MeV (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
//    Phip585CM7 = new GH1("Phip_585MeVCM7", "#phi_{p} Distribution for E_{#gamma} 585 #pm 5MeV (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
//    Phip595CM7 = new GH1("Phip_595MeVCM7", "#phi_{p} Distribution for E_{#gamma} 595 #pm 5MeV (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
//    Phip605CM7 = new GH1("Phip_605MeVCM7", "#phi_{p} Distribution for E_{#gamma} 605 #pm 5MeV (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
//    Phip615CM7 = new GH1("Phip_615MeVCM7", "#phi_{p} Distribution for E_{#gamma} 615 #pm 5MeV (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
//
//    // #phi_{p} dists across EGamma bins
//    Phip425CM8 = new GH1("Phip_425MeVCM8", "#phi_{p} Distribution for E_{#gamma} 425 #pm 5MeV (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
//    Phip435CM8 = new GH1("Phip_435MeVCM8", "#phi_{p} Distribution for E_{#gamma} 435 #pm 5MeV (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
//    Phip445CM8 = new GH1("Phip_445MeVCM8", "#phi_{p} Distribution for E_{#gamma} 445 #pm 5MeV (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
//    Phip455CM8 = new GH1("Phip_455MeVCM8", "#phi_{p} Distribution for E_{#gamma} 455 #pm 5MeV (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
//    Phip465CM8 = new GH1("Phip_465MeVCM8", "#phi_{p} Distribution for E_{#gamma} 465 #pm 5MeV (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
//    Phip475CM8 = new GH1("Phip_475MeVCM8", "#phi_{p} Distribution for E_{#gamma} 475 #pm 5MeV (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
//    Phip485CM8 = new GH1("Phip_485MeVCM8", "#phi_{p} Distribution for E_{#gamma} 485 #pm 5MeV (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
//    Phip495CM8 = new GH1("Phip_495MeVCM8", "#phi_{p} Distribution for E_{#gamma} 495 #pm 5MeV (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
//    Phip505CM8 = new GH1("Phip_505MeVCM8", "#phi_{p} Distribution for E_{#gamma} 505 #pm 5MeV (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
//    Phip515CM8 = new GH1("Phip_515MeVCM8", "#phi_{p} Distribution for E_{#gamma} 515 #pm 5MeV (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
//    Phip525CM8 = new GH1("Phip_525MeVCM8", "#phi_{p} Distribution for E_{#gamma} 525 #pm 5MeV (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
//    Phip535CM8 = new GH1("Phip_535MeVCM8", "#phi_{p} Distribution for E_{#gamma} 535 #pm 5MeV (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
//    Phip545CM8 = new GH1("Phip_545MeVCM8", "#phi_{p} Distribution for E_{#gamma} 545 #pm 5MeV (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
//    Phip555CM8 = new GH1("Phip_555MeVCM8", "#phi_{p} Distribution for E_{#gamma} 555 #pm 5MeV (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
//    Phip565CM8 = new GH1("Phip_565MeVCM8", "#phi_{p} Distribution for E_{#gamma} 565 #pm 5MeV (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
//    Phip575CM8 = new GH1("Phip_575MeVCM8", "#phi_{p} Distribution for E_{#gamma} 575 #pm 5MeV (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
//    Phip585CM8 = new GH1("Phip_585MeVCM8", "#phi_{p} Distribution for E_{#gamma} 585 #pm 5MeV (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
//    Phip595CM8 = new GH1("Phip_595MeVCM8", "#phi_{p} Distribution for E_{#gamma} 595 #pm 5MeV (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
//    Phip605CM8 = new GH1("Phip_605MeVCM8", "#phi_{p} Distribution for E_{#gamma} 605 #pm 5MeV (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
//    Phip615CM8 = new GH1("Phip_615MeVCM8", "#phi_{p} Distribution for E_{#gamma} 615 #pm 5MeV (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
//
//    // #phi_{p} dists across EGamma bins
//    Phip425CM9 = new GH1("Phip_425MeVCM9", "#phi_{p} Distribution for E_{#gamma} 425 #pm 5MeV (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
//    Phip435CM9 = new GH1("Phip_435MeVCM9", "#phi_{p} Distribution for E_{#gamma} 435 #pm 5MeV (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
//    Phip445CM9 = new GH1("Phip_445MeVCM9", "#phi_{p} Distribution for E_{#gamma} 445 #pm 5MeV (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
//    Phip455CM9 = new GH1("Phip_455MeVCM9", "#phi_{p} Distribution for E_{#gamma} 455 #pm 5MeV (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
//    Phip465CM9 = new GH1("Phip_465MeVCM9", "#phi_{p} Distribution for E_{#gamma} 465 #pm 5MeV (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
//    Phip475CM9 = new GH1("Phip_475MeVCM9", "#phi_{p} Distribution for E_{#gamma} 475 #pm 5MeV (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
//    Phip485CM9 = new GH1("Phip_485MeVCM9", "#phi_{p} Distribution for E_{#gamma} 485 #pm 5MeV (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
//    Phip495CM9 = new GH1("Phip_495MeVCM9", "#phi_{p} Distribution for E_{#gamma} 495 #pm 5MeV (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
//    Phip505CM9 = new GH1("Phip_505MeVCM9", "#phi_{p} Distribution for E_{#gamma} 505 #pm 5MeV (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
//    Phip515CM9 = new GH1("Phip_515MeVCM9", "#phi_{p} Distribution for E_{#gamma} 515 #pm 5MeV (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
//    Phip525CM9 = new GH1("Phip_525MeVCM9", "#phi_{p} Distribution for E_{#gamma} 525 #pm 5MeV (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
//    Phip535CM9 = new GH1("Phip_535MeVCM9", "#phi_{p} Distribution for E_{#gamma} 535 #pm 5MeV (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
//    Phip545CM9 = new GH1("Phip_545MeVCM9", "#phi_{p} Distribution for E_{#gamma} 545 #pm 5MeV (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
//    Phip555CM9 = new GH1("Phip_555MeVCM9", "#phi_{p} Distribution for E_{#gamma} 555 #pm 5MeV (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
//    Phip565CM9 = new GH1("Phip_565MeVCM9", "#phi_{p} Distribution for E_{#gamma} 565 #pm 5MeV (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
//    Phip575CM9 = new GH1("Phip_575MeVCM9", "#phi_{p} Distribution for E_{#gamma} 575 #pm 5MeV (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
//    Phip585CM9 = new GH1("Phip_585MeVCM9", "#phi_{p} Distribution for E_{#gamma} 585 #pm 5MeV (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
//    Phip595CM9 = new GH1("Phip_595MeVCM9", "#phi_{p} Distribution for E_{#gamma} 595 #pm 5MeV (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
//    Phip605CM9 = new GH1("Phip_605MeVCM9", "#phi_{p} Distribution for E_{#gamma} 605 #pm 5MeV (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
//    Phip615CM9 = new GH1("Phip_615MeVCM9", "#phi_{p} Distribution for E_{#gamma} 615 #pm 5MeV (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
//
//    // #phi_{p} dists across EGamma bins
//    Phip425CM10 = new GH1("Phip_425MeVCM10", "#phi_{p} Distribution for E_{#gamma} 425 #pm 5MeV (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
//    Phip435CM10 = new GH1("Phip_435MeVCM10", "#phi_{p} Distribution for E_{#gamma} 435 #pm 5MeV (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
//    Phip445CM10 = new GH1("Phip_445MeVCM10", "#phi_{p} Distribution for E_{#gamma} 445 #pm 5MeV (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
//    Phip455CM10 = new GH1("Phip_455MeVCM10", "#phi_{p} Distribution for E_{#gamma} 455 #pm 5MeV (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
//    Phip465CM10 = new GH1("Phip_465MeVCM10", "#phi_{p} Distribution for E_{#gamma} 465 #pm 5MeV (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
//    Phip475CM10 = new GH1("Phip_475MeVCM10", "#phi_{p} Distribution for E_{#gamma} 475 #pm 5MeV (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
//    Phip485CM10 = new GH1("Phip_485MeVCM10", "#phi_{p} Distribution for E_{#gamma} 485 #pm 5MeV (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
//    Phip495CM10 = new GH1("Phip_495MeVCM10", "#phi_{p} Distribution for E_{#gamma} 495 #pm 5MeV (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
//    Phip505CM10 = new GH1("Phip_505MeVCM10", "#phi_{p} Distribution for E_{#gamma} 505 #pm 5MeV (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
//    Phip515CM10 = new GH1("Phip_515MeVCM10", "#phi_{p} Distribution for E_{#gamma} 515 #pm 5MeV (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
//    Phip525CM10 = new GH1("Phip_525MeVCM10", "#phi_{p} Distribution for E_{#gamma} 525 #pm 5MeV (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
//    Phip535CM10 = new GH1("Phip_535MeVCM10", "#phi_{p} Distribution for E_{#gamma} 535 #pm 5MeV (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
//    Phip545CM10 = new GH1("Phip_545MeVCM10", "#phi_{p} Distribution for E_{#gamma} 545 #pm 5MeV (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
//    Phip555CM10 = new GH1("Phip_555MeVCM10", "#phi_{p} Distribution for E_{#gamma} 555 #pm 5MeV (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
//    Phip565CM10 = new GH1("Phip_565MeVCM10", "#phi_{p} Distribution for E_{#gamma} 565 #pm 5MeV (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
//    Phip575CM10 = new GH1("Phip_575MeVCM10", "#phi_{p} Distribution for E_{#gamma} 575 #pm 5MeV (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
//    Phip585CM10 = new GH1("Phip_585MeVCM10", "#phi_{p} Distribution for E_{#gamma} 585 #pm 5MeV (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
//    Phip595CM10 = new GH1("Phip_595MeVCM10", "#phi_{p} Distribution for E_{#gamma} 595 #pm 5MeV (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
//    Phip605CM10 = new GH1("Phip_605MeVCM10", "#phi_{p} Distribution for E_{#gamma} 605 #pm 5MeV (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
//    Phip615CM10 = new GH1("Phip_615MeVCM10", "#phi_{p} Distribution for E_{#gamma} 615 #pm 5MeV (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
// Proton Phi dists across EGamma bins
    Phip425CM1PosHel = new GH1("Phip_425MeVCM1PosHel", "#phi_{p} Distribution for E_{#gamma} 425 #pm 5MeV +ve Hel (Cos#theta_{CM}1-0.8)", 10, -180, 180);
    Phip435CM1PosHel = new GH1("Phip_435MeVCM1PosHel", "#phi_{p} Distribution for E_{#gamma} 435 #pm 5MeV +ve Hel (Cos#theta_{CM}1-0.8)", 10, -180, 180);
    Phip445CM1PosHel = new GH1("Phip_445MeVCM1PosHel", "#phi_{p} Distribution for E_{#gamma} 445 #pm 5MeV +ve Hel (Cos#theta_{CM}1-0.8)", 10, -180, 180);
    Phip455CM1PosHel = new GH1("Phip_455MeVCM1PosHel", "#phi_{p} Distribution for E_{#gamma} 455 #pm 5MeV +ve Hel (Cos#theta_{CM}1-0.8)", 10, -180, 180);
    Phip465CM1PosHel = new GH1("Phip_465MeVCM1PosHel", "#phi_{p} Distribution for E_{#gamma} 465 #pm 5MeV +ve Hel (Cos#theta_{CM}1-0.8)", 10, -180, 180);
    Phip475CM1PosHel = new GH1("Phip_475MeVCM1PosHel", "#phi_{p} Distribution for E_{#gamma} 475 #pm 5MeV +ve Hel (Cos#theta_{CM}1-0.8)", 10, -180, 180);
    Phip485CM1PosHel = new GH1("Phip_485MeVCM1PosHel", "#phi_{p} Distribution for E_{#gamma} 485 #pm 5MeV +ve Hel (Cos#theta_{CM}1-0.8)", 10, -180, 180);
    Phip495CM1PosHel = new GH1("Phip_495MeVCM1PosHel", "#phi_{p} Distribution for E_{#gamma} 495 #pm 5MeV +ve Hel (Cos#theta_{CM}1-0.8)", 10, -180, 180);
    Phip505CM1PosHel = new GH1("Phip_505MeVCM1PosHel", "#phi_{p} Distribution for E_{#gamma} 505 #pm 5MeV +ve Hel (Cos#theta_{CM}1-0.8)", 10, -180, 180);
    Phip515CM1PosHel = new GH1("Phip_515MeVCM1PosHel", "#phi_{p} Distribution for E_{#gamma} 515 #pm 5MeV +ve Hel (Cos#theta_{CM}1-0.8)", 10, -180, 180);
    Phip525CM1PosHel = new GH1("Phip_525MeVCM1PosHel", "#phi_{p} Distribution for E_{#gamma} 525 #pm 5MeV +ve Hel (Cos#theta_{CM}1-0.8)", 10, -180, 180);
    Phip535CM1PosHel = new GH1("Phip_535MeVCM1PosHel", "#phi_{p} Distribution for E_{#gamma} 535 #pm 5MeV +ve Hel (Cos#theta_{CM}1-0.8)", 10, -180, 180);
    Phip545CM1PosHel = new GH1("Phip_545MeVCM1PosHel", "#phi_{p} Distribution for E_{#gamma} 545 #pm 5MeV +ve Hel (Cos#theta_{CM}1-0.8)", 10, -180, 180);
    Phip555CM1PosHel = new GH1("Phip_555MeVCM1PosHel", "#phi_{p} Distribution for E_{#gamma} 555 #pm 5MeV +ve Hel (Cos#theta_{CM}1-0.8)", 10, -180, 180);
    Phip565CM1PosHel = new GH1("Phip_565MeVCM1PosHel", "#phi_{p} Distribution for E_{#gamma} 565 #pm 5MeV +ve Hel (Cos#theta_{CM}1-0.8)", 10, -180, 180);
    Phip575CM1PosHel = new GH1("Phip_575MeVCM1PosHel", "#phi_{p} Distribution for E_{#gamma} 575 #pm 5MeV +ve Hel (Cos#theta_{CM}1-0.8)", 10, -180, 180);
    Phip585CM1PosHel = new GH1("Phip_585MeVCM1PosHel", "#phi_{p} Distribution for E_{#gamma} 585 #pm 5MeV +ve Hel (Cos#theta_{CM}1-0.8)", 10, -180, 180);
    Phip595CM1PosHel = new GH1("Phip_595MeVCM1PosHel", "#phi_{p} Distribution for E_{#gamma} 595 #pm 5MeV +ve Hel (Cos#theta_{CM}1-0.8)", 10, -180, 180);
    Phip605CM1PosHel = new GH1("Phip_605MeVCM1PosHel", "#phi_{p} Distribution for E_{#gamma} 605 #pm 5MeV +ve Hel (Cos#theta_{CM}1-0.8)", 10, -180, 180);
    Phip615CM1PosHel = new GH1("Phip_615MeVCM1PosHel", "#phi_{p} Distribution for E_{#gamma} 615 #pm 5MeV +ve Hel (Cos#theta_{CM}1-0.8)", 10, -180, 180);

    // #phi_{p} dists across EGamma bins
    Phip425CM2PosHel = new GH1("Phip_425MeVCM2PosHel", "#phi_{p} Distribution for E_{#gamma} 425 #pm 5MeV +ve Hel (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
    Phip435CM2PosHel = new GH1("Phip_435MeVCM2PosHel", "#phi_{p} Distribution for E_{#gamma} 435 #pm 5MeV +ve Hel (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
    Phip445CM2PosHel = new GH1("Phip_445MeVCM2PosHel", "#phi_{p} Distribution for E_{#gamma} 445 #pm 5MeV +ve Hel (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
    Phip455CM2PosHel = new GH1("Phip_455MeVCM2PosHel", "#phi_{p} Distribution for E_{#gamma} 455 #pm 5MeV +ve Hel (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
    Phip465CM2PosHel = new GH1("Phip_465MeVCM2PosHel", "#phi_{p} Distribution for E_{#gamma} 465 #pm 5MeV +ve Hel (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
    Phip475CM2PosHel = new GH1("Phip_475MeVCM2PosHel", "#phi_{p} Distribution for E_{#gamma} 475 #pm 5MeV +ve Hel (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
    Phip485CM2PosHel = new GH1("Phip_485MeVCM2PosHel", "#phi_{p} Distribution for E_{#gamma} 485 #pm 5MeV +ve Hel (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
    Phip495CM2PosHel = new GH1("Phip_495MeVCM2PosHel", "#phi_{p} Distribution for E_{#gamma} 495 #pm 5MeV +ve Hel (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
    Phip505CM2PosHel = new GH1("Phip_505MeVCM2PosHel", "#phi_{p} Distribution for E_{#gamma} 505 #pm 5MeV +ve Hel (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
    Phip515CM2PosHel = new GH1("Phip_515MeVCM2PosHel", "#phi_{p} Distribution for E_{#gamma} 515 #pm 5MeV +ve Hel (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
    Phip525CM2PosHel = new GH1("Phip_525MeVCM2PosHel", "#phi_{p} Distribution for E_{#gamma} 525 #pm 5MeV +ve Hel (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
    Phip535CM2PosHel = new GH1("Phip_535MeVCM2PosHel", "#phi_{p} Distribution for E_{#gamma} 535 #pm 5MeV +ve Hel (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
    Phip545CM2PosHel = new GH1("Phip_545MeVCM2PosHel", "#phi_{p} Distribution for E_{#gamma} 545 #pm 5MeV +ve Hel (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
    Phip555CM2PosHel = new GH1("Phip_555MeVCM2PosHel", "#phi_{p} Distribution for E_{#gamma} 555 #pm 5MeV +ve Hel (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
    Phip565CM2PosHel = new GH1("Phip_565MeVCM2PosHel", "#phi_{p} Distribution for E_{#gamma} 565 #pm 5MeV +ve Hel (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
    Phip575CM2PosHel = new GH1("Phip_575MeVCM2PosHel", "#phi_{p} Distribution for E_{#gamma} 575 #pm 5MeV +ve Hel (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
    Phip585CM2PosHel = new GH1("Phip_585MeVCM2PosHel", "#phi_{p} Distribution for E_{#gamma} 585 #pm 5MeV +ve Hel (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
    Phip595CM2PosHel = new GH1("Phip_595MeVCM2PosHel", "#phi_{p} Distribution for E_{#gamma} 595 #pm 5MeV +ve Hel (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
    Phip605CM2PosHel = new GH1("Phip_605MeVCM2PosHel", "#phi_{p} Distribution for E_{#gamma} 605 #pm 5MeV +ve Hel (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
    Phip615CM2PosHel = new GH1("Phip_615MeVCM2PosHel", "#phi_{p} Distribution for E_{#gamma} 615 #pm 5MeV +ve Hel (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);

    // #phi_{p} dists across EGamma bins
    Phip425CM3PosHel = new GH1("Phip_425MeVCM3PosHel", "#phi_{p} Distribution for E_{#gamma} 425 #pm 5MeV +ve Hel (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
    Phip435CM3PosHel = new GH1("Phip_435MeVCM3PosHel", "#phi_{p} Distribution for E_{#gamma} 435 #pm 5MeV +ve Hel (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
    Phip445CM3PosHel = new GH1("Phip_445MeVCM3PosHel", "#phi_{p} Distribution for E_{#gamma} 445 #pm 5MeV +ve Hel (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
    Phip455CM3PosHel = new GH1("Phip_455MeVCM3PosHel", "#phi_{p} Distribution for E_{#gamma} 455 #pm 5MeV +ve Hel (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
    Phip465CM3PosHel = new GH1("Phip_465MeVCM3PosHel", "#phi_{p} Distribution for E_{#gamma} 465 #pm 5MeV +ve Hel (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
    Phip475CM3PosHel = new GH1("Phip_475MeVCM3PosHel", "#phi_{p} Distribution for E_{#gamma} 475 #pm 5MeV +ve Hel (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
    Phip485CM3PosHel = new GH1("Phip_485MeVCM3PosHel", "#phi_{p} Distribution for E_{#gamma} 485 #pm 5MeV +ve Hel (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
    Phip495CM3PosHel = new GH1("Phip_495MeVCM3PosHel", "#phi_{p} Distribution for E_{#gamma} 495 #pm 5MeV +ve Hel (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
    Phip505CM3PosHel = new GH1("Phip_505MeVCM3PosHel", "#phi_{p} Distribution for E_{#gamma} 505 #pm 5MeV +ve Hel (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
    Phip515CM3PosHel = new GH1("Phip_515MeVCM3PosHel", "#phi_{p} Distribution for E_{#gamma} 515 #pm 5MeV +ve Hel (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
    Phip525CM3PosHel = new GH1("Phip_525MeVCM3PosHel", "#phi_{p} Distribution for E_{#gamma} 525 #pm 5MeV +ve Hel (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
    Phip535CM3PosHel = new GH1("Phip_535MeVCM3PosHel", "#phi_{p} Distribution for E_{#gamma} 535 #pm 5MeV +ve Hel (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
    Phip545CM3PosHel = new GH1("Phip_545MeVCM3PosHel", "#phi_{p} Distribution for E_{#gamma} 545 #pm 5MeV +ve Hel (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
    Phip555CM3PosHel = new GH1("Phip_555MeVCM3PosHel", "#phi_{p} Distribution for E_{#gamma} 555 #pm 5MeV +ve Hel (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
    Phip565CM3PosHel = new GH1("Phip_565MeVCM3PosHel", "#phi_{p} Distribution for E_{#gamma} 565 #pm 5MeV +ve Hel (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
    Phip575CM3PosHel = new GH1("Phip_575MeVCM3PosHel", "#phi_{p} Distribution for E_{#gamma} 575 #pm 5MeV +ve Hel (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
    Phip585CM3PosHel = new GH1("Phip_585MeVCM3PosHel", "#phi_{p} Distribution for E_{#gamma} 585 #pm 5MeV +ve Hel (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
    Phip595CM3PosHel = new GH1("Phip_595MeVCM3PosHel", "#phi_{p} Distribution for E_{#gamma} 595 #pm 5MeV +ve Hel (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
    Phip605CM3PosHel = new GH1("Phip_605MeVCM3PosHel", "#phi_{p} Distribution for E_{#gamma} 605 #pm 5MeV +ve Hel (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
    Phip615CM3PosHel = new GH1("Phip_615MeVCM3PosHel", "#phi_{p} Distribution for E_{#gamma} 615 #pm 5MeV +ve Hel (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);

    // #phi_{p} dists across EGamma bins
    Phip425CM4PosHel = new GH1("Phip_425MeVCM4PosHel", "#phi_{p} Distribution for E_{#gamma} 425 #pm 5MeV +ve Hel (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
    Phip435CM4PosHel = new GH1("Phip_435MeVCM4PosHel", "#phi_{p} Distribution for E_{#gamma} 435 #pm 5MeV +ve Hel (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
    Phip445CM4PosHel = new GH1("Phip_445MeVCM4PosHel", "#phi_{p} Distribution for E_{#gamma} 445 #pm 5MeV +ve Hel (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
    Phip455CM4PosHel = new GH1("Phip_455MeVCM4PosHel", "#phi_{p} Distribution for E_{#gamma} 455 #pm 5MeV +ve Hel (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
    Phip465CM4PosHel = new GH1("Phip_465MeVCM4PosHel", "#phi_{p} Distribution for E_{#gamma} 465 #pm 5MeV +ve Hel (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
    Phip475CM4PosHel = new GH1("Phip_475MeVCM4PosHel", "#phi_{p} Distribution for E_{#gamma} 475 #pm 5MeV +ve Hel (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
    Phip485CM4PosHel = new GH1("Phip_485MeVCM4PosHel", "#phi_{p} Distribution for E_{#gamma} 485 #pm 5MeV +ve Hel (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
    Phip495CM4PosHel = new GH1("Phip_495MeVCM4PosHel", "#phi_{p} Distribution for E_{#gamma} 495 #pm 5MeV +ve Hel (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
    Phip505CM4PosHel = new GH1("Phip_505MeVCM4PosHel", "#phi_{p} Distribution for E_{#gamma} 505 #pm 5MeV +ve Hel (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
    Phip515CM4PosHel = new GH1("Phip_515MeVCM4PosHel", "#phi_{p} Distribution for E_{#gamma} 515 #pm 5MeV +ve Hel (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
    Phip525CM4PosHel = new GH1("Phip_525MeVCM4PosHel", "#phi_{p} Distribution for E_{#gamma} 525 #pm 5MeV +ve Hel (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
    Phip535CM4PosHel = new GH1("Phip_535MeVCM4PosHel", "#phi_{p} Distribution for E_{#gamma} 535 #pm 5MeV +ve Hel (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
    Phip545CM4PosHel = new GH1("Phip_545MeVCM4PosHel", "#phi_{p} Distribution for E_{#gamma} 545 #pm 5MeV +ve Hel (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
    Phip555CM4PosHel = new GH1("Phip_555MeVCM4PosHel", "#phi_{p} Distribution for E_{#gamma} 555 #pm 5MeV +ve Hel (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
    Phip565CM4PosHel = new GH1("Phip_565MeVCM4PosHel", "#phi_{p} Distribution for E_{#gamma} 565 #pm 5MeV +ve Hel (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
    Phip575CM4PosHel = new GH1("Phip_575MeVCM4PosHel", "#phi_{p} Distribution for E_{#gamma} 575 #pm 5MeV +ve Hel (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
    Phip585CM4PosHel = new GH1("Phip_585MeVCM4PosHel", "#phi_{p} Distribution for E_{#gamma} 585 #pm 5MeV +ve Hel (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
    Phip595CM4PosHel = new GH1("Phip_595MeVCM4PosHel", "#phi_{p} Distribution for E_{#gamma} 595 #pm 5MeV +ve Hel (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
    Phip605CM4PosHel = new GH1("Phip_605MeVCM4PosHel", "#phi_{p} Distribution for E_{#gamma} 605 #pm 5MeV +ve Hel (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
    Phip615CM4PosHel = new GH1("Phip_615MeVCM4PosHel", "#phi_{p} Distribution for E_{#gamma} 615 #pm 5MeV +ve Hel (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);

    // #phi_{p} dists across EGamma bins
    Phip425CM5PosHel = new GH1("Phip_425MeVCM5PosHel", "#phi_{p} Distribution for E_{#gamma} 425 #pm 5MeV +ve Hel (Cos#theta_{CM}0.2-0)", 10, -180, 180);
    Phip435CM5PosHel = new GH1("Phip_435MeVCM5PosHel", "#phi_{p} Distribution for E_{#gamma} 435 #pm 5MeV +ve Hel (Cos#theta_{CM}0.2-0)", 10, -180, 180);
    Phip445CM5PosHel = new GH1("Phip_445MeVCM5PosHel", "#phi_{p} Distribution for E_{#gamma} 445 #pm 5MeV +ve Hel (Cos#theta_{CM}0.2-0)", 10, -180, 180);
    Phip455CM5PosHel = new GH1("Phip_455MeVCM5PosHel", "#phi_{p} Distribution for E_{#gamma} 455 #pm 5MeV +ve Hel (Cos#theta_{CM}0.2-0)", 10, -180, 180);
    Phip465CM5PosHel = new GH1("Phip_465MeVCM5PosHel", "#phi_{p} Distribution for E_{#gamma} 465 #pm 5MeV +ve Hel (Cos#theta_{CM}0.2-0)", 10, -180, 180);
    Phip475CM5PosHel = new GH1("Phip_475MeVCM5PosHel", "#phi_{p} Distribution for E_{#gamma} 475 #pm 5MeV +ve Hel (Cos#theta_{CM}0.2-0)", 10, -180, 180);
    Phip485CM5PosHel = new GH1("Phip_485MeVCM5PosHel", "#phi_{p} Distribution for E_{#gamma} 485 #pm 5MeV +ve Hel (Cos#theta_{CM}0.2-0)", 10, -180, 180);
    Phip495CM5PosHel = new GH1("Phip_495MeVCM5PosHel", "#phi_{p} Distribution for E_{#gamma} 495 #pm 5MeV +ve Hel (Cos#theta_{CM}0.2-0)", 10, -180, 180);
    Phip505CM5PosHel = new GH1("Phip_505MeVCM5PosHel", "#phi_{p} Distribution for E_{#gamma} 505 #pm 5MeV +ve Hel (Cos#theta_{CM}0.2-0)", 10, -180, 180);
    Phip515CM5PosHel = new GH1("Phip_515MeVCM5PosHel", "#phi_{p} Distribution for E_{#gamma} 515 #pm 5MeV +ve Hel (Cos#theta_{CM}0.2-0)", 10, -180, 180);
    Phip525CM5PosHel = new GH1("Phip_525MeVCM5PosHel", "#phi_{p} Distribution for E_{#gamma} 525 #pm 5MeV +ve Hel (Cos#theta_{CM}0.2-0)", 10, -180, 180);
    Phip535CM5PosHel = new GH1("Phip_535MeVCM5PosHel", "#phi_{p} Distribution for E_{#gamma} 535 #pm 5MeV +ve Hel (Cos#theta_{CM}0.2-0)", 10, -180, 180);
    Phip545CM5PosHel = new GH1("Phip_545MeVCM5PosHel", "#phi_{p} Distribution for E_{#gamma} 545 #pm 5MeV +ve Hel (Cos#theta_{CM}0.2-0)", 10, -180, 180);
    Phip555CM5PosHel = new GH1("Phip_555MeVCM5PosHel", "#phi_{p} Distribution for E_{#gamma} 555 #pm 5MeV +ve Hel (Cos#theta_{CM}0.2-0)", 10, -180, 180);
    Phip565CM5PosHel = new GH1("Phip_565MeVCM5PosHel", "#phi_{p} Distribution for E_{#gamma} 565 #pm 5MeV +ve Hel (Cos#theta_{CM}0.2-0)", 10, -180, 180);
    Phip575CM5PosHel = new GH1("Phip_575MeVCM5PosHel", "#phi_{p} Distribution for E_{#gamma} 575 #pm 5MeV +ve Hel (Cos#theta_{CM}0.2-0)", 10, -180, 180);
    Phip585CM5PosHel = new GH1("Phip_585MeVCM5PosHel", "#phi_{p} Distribution for E_{#gamma} 585 #pm 5MeV +ve Hel (Cos#theta_{CM}0.2-0)", 10, -180, 180);
    Phip595CM5PosHel = new GH1("Phip_595MeVCM5PosHel", "#phi_{p} Distribution for E_{#gamma} 595 #pm 5MeV +ve Hel (Cos#theta_{CM}0.2-0)", 10, -180, 180);
    Phip605CM5PosHel = new GH1("Phip_605MeVCM5PosHel", "#phi_{p} Distribution for E_{#gamma} 605 #pm 5MeV +ve Hel (Cos#theta_{CM}0.2-0)", 10, -180, 180);
    Phip615CM5PosHel = new GH1("Phip_615MeVCM5PosHel", "#phi_{p} Distribution for E_{#gamma} 615 #pm 5MeV +ve Hel (Cos#theta_{CM}0.2-0)", 10, -180, 180);

    // #phi_{p} dists across EGamma bins
    Phip425CM6PosHel = new GH1("Phip_425MeVCM6PosHel", "#phi_{p} Distribution for E_{#gamma} 425 #pm 5MeV +ve Hel (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
    Phip435CM6PosHel = new GH1("Phip_435MeVCM6PosHel", "#phi_{p} Distribution for E_{#gamma} 435 #pm 5MeV +ve Hel (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
    Phip445CM6PosHel = new GH1("Phip_445MeVCM6PosHel", "#phi_{p} Distribution for E_{#gamma} 445 #pm 5MeV +ve Hel (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
    Phip455CM6PosHel = new GH1("Phip_455MeVCM6PosHel", "#phi_{p} Distribution for E_{#gamma} 455 #pm 5MeV +ve Hel (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
    Phip465CM6PosHel = new GH1("Phip_465MeVCM6PosHel", "#phi_{p} Distribution for E_{#gamma} 465 #pm 5MeV +ve Hel (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
    Phip475CM6PosHel = new GH1("Phip_475MeVCM6PosHel", "#phi_{p} Distribution for E_{#gamma} 475 #pm 5MeV +ve Hel (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
    Phip485CM6PosHel = new GH1("Phip_485MeVCM6PosHel", "#phi_{p} Distribution for E_{#gamma} 485 #pm 5MeV +ve Hel (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
    Phip495CM6PosHel = new GH1("Phip_495MeVCM6PosHel", "#phi_{p} Distribution for E_{#gamma} 495 #pm 5MeV +ve Hel (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
    Phip505CM6PosHel = new GH1("Phip_505MeVCM6PosHel", "#phi_{p} Distribution for E_{#gamma} 505 #pm 5MeV +ve Hel (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
    Phip515CM6PosHel = new GH1("Phip_515MeVCM6PosHel", "#phi_{p} Distribution for E_{#gamma} 515 #pm 5MeV +ve Hel (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
    Phip525CM6PosHel = new GH1("Phip_525MeVCM6PosHel", "#phi_{p} Distribution for E_{#gamma} 525 #pm 5MeV +ve Hel (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
    Phip535CM6PosHel = new GH1("Phip_535MeVCM6PosHel", "#phi_{p} Distribution for E_{#gamma} 535 #pm 5MeV +ve Hel (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
    Phip545CM6PosHel = new GH1("Phip_545MeVCM6PosHel", "#phi_{p} Distribution for E_{#gamma} 545 #pm 5MeV +ve Hel (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
    Phip555CM6PosHel = new GH1("Phip_555MeVCM6PosHel", "#phi_{p} Distribution for E_{#gamma} 555 #pm 5MeV +ve Hel (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
    Phip565CM6PosHel = new GH1("Phip_565MeVCM6PosHel", "#phi_{p} Distribution for E_{#gamma} 565 #pm 5MeV +ve Hel (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
    Phip575CM6PosHel = new GH1("Phip_575MeVCM6PosHel", "#phi_{p} Distribution for E_{#gamma} 575 #pm 5MeV +ve Hel (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
    Phip585CM6PosHel = new GH1("Phip_585MeVCM6PosHel", "#phi_{p} Distribution for E_{#gamma} 585 #pm 5MeV +ve Hel (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
    Phip595CM6PosHel = new GH1("Phip_595MeVCM6PosHel", "#phi_{p} Distribution for E_{#gamma} 595 #pm 5MeV +ve Hel (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
    Phip605CM6PosHel = new GH1("Phip_605MeVCM6PosHel", "#phi_{p} Distribution for E_{#gamma} 605 #pm 5MeV +ve Hel (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
    Phip615CM6PosHel = new GH1("Phip_615MeVCM6PosHel", "#phi_{p} Distribution for E_{#gamma} 615 #pm 5MeV +ve Hel (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);

    // #phi_{p} dists across EGamma bins
    Phip425CM7PosHel = new GH1("Phip_425MeVCM7PosHel", "#phi_{p} Distribution for E_{#gamma} 425 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
    Phip435CM7PosHel = new GH1("Phip_435MeVCM7PosHel", "#phi_{p} Distribution for E_{#gamma} 435 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
    Phip445CM7PosHel = new GH1("Phip_445MeVCM7PosHel", "#phi_{p} Distribution for E_{#gamma} 445 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
    Phip455CM7PosHel = new GH1("Phip_455MeVCM7PosHel", "#phi_{p} Distribution for E_{#gamma} 455 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
    Phip465CM7PosHel = new GH1("Phip_465MeVCM7PosHel", "#phi_{p} Distribution for E_{#gamma} 465 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
    Phip475CM7PosHel = new GH1("Phip_475MeVCM7PosHel", "#phi_{p} Distribution for E_{#gamma} 475 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
    Phip485CM7PosHel = new GH1("Phip_485MeVCM7PosHel", "#phi_{p} Distribution for E_{#gamma} 485 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
    Phip495CM7PosHel = new GH1("Phip_495MeVCM7PosHel", "#phi_{p} Distribution for E_{#gamma} 495 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
    Phip505CM7PosHel = new GH1("Phip_505MeVCM7PosHel", "#phi_{p} Distribution for E_{#gamma} 505 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
    Phip515CM7PosHel = new GH1("Phip_515MeVCM7PosHel", "#phi_{p} Distribution for E_{#gamma} 515 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
    Phip525CM7PosHel = new GH1("Phip_525MeVCM7PosHel", "#phi_{p} Distribution for E_{#gamma} 525 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
    Phip535CM7PosHel = new GH1("Phip_535MeVCM7PosHel", "#phi_{p} Distribution for E_{#gamma} 535 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
    Phip545CM7PosHel = new GH1("Phip_545MeVCM7PosHel", "#phi_{p} Distribution for E_{#gamma} 545 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
    Phip555CM7PosHel = new GH1("Phip_555MeVCM7PosHel", "#phi_{p} Distribution for E_{#gamma} 555 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
    Phip565CM7PosHel = new GH1("Phip_565MeVCM7PosHel", "#phi_{p} Distribution for E_{#gamma} 565 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
    Phip575CM7PosHel = new GH1("Phip_575MeVCM7PosHel", "#phi_{p} Distribution for E_{#gamma} 575 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
    Phip585CM7PosHel = new GH1("Phip_585MeVCM7PosHel", "#phi_{p} Distribution for E_{#gamma} 585 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
    Phip595CM7PosHel = new GH1("Phip_595MeVCM7PosHel", "#phi_{p} Distribution for E_{#gamma} 595 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
    Phip605CM7PosHel = new GH1("Phip_605MeVCM7PosHel", "#phi_{p} Distribution for E_{#gamma} 605 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
    Phip615CM7PosHel = new GH1("Phip_615MeVCM7PosHel", "#phi_{p} Distribution for E_{#gamma} 615 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);

    // #phi_{p} dists across EGamma bins
    Phip425CM8PosHel = new GH1("Phip_425MeVCM8PosHel", "#phi_{p} Distribution for E_{#gamma} 425 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
    Phip435CM8PosHel = new GH1("Phip_435MeVCM8PosHel", "#phi_{p} Distribution for E_{#gamma} 435 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
    Phip445CM8PosHel = new GH1("Phip_445MeVCM8PosHel", "#phi_{p} Distribution for E_{#gamma} 445 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
    Phip455CM8PosHel = new GH1("Phip_455MeVCM8PosHel", "#phi_{p} Distribution for E_{#gamma} 455 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
    Phip465CM8PosHel = new GH1("Phip_465MeVCM8PosHel", "#phi_{p} Distribution for E_{#gamma} 465 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
    Phip475CM8PosHel = new GH1("Phip_475MeVCM8PosHel", "#phi_{p} Distribution for E_{#gamma} 475 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
    Phip485CM8PosHel = new GH1("Phip_485MeVCM8PosHel", "#phi_{p} Distribution for E_{#gamma} 485 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
    Phip495CM8PosHel = new GH1("Phip_495MeVCM8PosHel", "#phi_{p} Distribution for E_{#gamma} 495 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
    Phip505CM8PosHel = new GH1("Phip_505MeVCM8PosHel", "#phi_{p} Distribution for E_{#gamma} 505 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
    Phip515CM8PosHel = new GH1("Phip_515MeVCM8PosHel", "#phi_{p} Distribution for E_{#gamma} 515 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
    Phip525CM8PosHel = new GH1("Phip_525MeVCM8PosHel", "#phi_{p} Distribution for E_{#gamma} 525 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
    Phip535CM8PosHel = new GH1("Phip_535MeVCM8PosHel", "#phi_{p} Distribution for E_{#gamma} 535 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
    Phip545CM8PosHel = new GH1("Phip_545MeVCM8PosHel", "#phi_{p} Distribution for E_{#gamma} 545 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
    Phip555CM8PosHel = new GH1("Phip_555MeVCM8PosHel", "#phi_{p} Distribution for E_{#gamma} 555 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
    Phip565CM8PosHel = new GH1("Phip_565MeVCM8PosHel", "#phi_{p} Distribution for E_{#gamma} 565 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
    Phip575CM8PosHel = new GH1("Phip_575MeVCM8PosHel", "#phi_{p} Distribution for E_{#gamma} 575 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
    Phip585CM8PosHel = new GH1("Phip_585MeVCM8PosHel", "#phi_{p} Distribution for E_{#gamma} 585 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
    Phip595CM8PosHel = new GH1("Phip_595MeVCM8PosHel", "#phi_{p} Distribution for E_{#gamma} 595 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
    Phip605CM8PosHel = new GH1("Phip_605MeVCM8PosHel", "#phi_{p} Distribution for E_{#gamma} 605 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
    Phip615CM8PosHel = new GH1("Phip_615MeVCM8PosHel", "#phi_{p} Distribution for E_{#gamma} 615 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);

    // #phi_{p} dists across EGamma bins
    Phip425CM9PosHel = new GH1("Phip_425MeVCM9PosHel", "#phi_{p} Distribution for E_{#gamma} 425 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
    Phip435CM9PosHel = new GH1("Phip_435MeVCM9PosHel", "#phi_{p} Distribution for E_{#gamma} 435 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
    Phip445CM9PosHel = new GH1("Phip_445MeVCM9PosHel", "#phi_{p} Distribution for E_{#gamma} 445 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
    Phip455CM9PosHel = new GH1("Phip_455MeVCM9PosHel", "#phi_{p} Distribution for E_{#gamma} 455 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
    Phip465CM9PosHel = new GH1("Phip_465MeVCM9PosHel", "#phi_{p} Distribution for E_{#gamma} 465 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
    Phip475CM9PosHel = new GH1("Phip_475MeVCM9PosHel", "#phi_{p} Distribution for E_{#gamma} 475 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
    Phip485CM9PosHel = new GH1("Phip_485MeVCM9PosHel", "#phi_{p} Distribution for E_{#gamma} 485 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
    Phip495CM9PosHel = new GH1("Phip_495MeVCM9PosHel", "#phi_{p} Distribution for E_{#gamma} 495 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
    Phip505CM9PosHel = new GH1("Phip_505MeVCM9PosHel", "#phi_{p} Distribution for E_{#gamma} 505 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
    Phip515CM9PosHel = new GH1("Phip_515MeVCM9PosHel", "#phi_{p} Distribution for E_{#gamma} 515 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
    Phip525CM9PosHel = new GH1("Phip_525MeVCM9PosHel", "#phi_{p} Distribution for E_{#gamma} 525 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
    Phip535CM9PosHel = new GH1("Phip_535MeVCM9PosHel", "#phi_{p} Distribution for E_{#gamma} 535 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
    Phip545CM9PosHel = new GH1("Phip_545MeVCM9PosHel", "#phi_{p} Distribution for E_{#gamma} 545 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
    Phip555CM9PosHel = new GH1("Phip_555MeVCM9PosHel", "#phi_{p} Distribution for E_{#gamma} 555 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
    Phip565CM9PosHel = new GH1("Phip_565MeVCM9PosHel", "#phi_{p} Distribution for E_{#gamma} 565 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
    Phip575CM9PosHel = new GH1("Phip_575MeVCM9PosHel", "#phi_{p} Distribution for E_{#gamma} 575 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
    Phip585CM9PosHel = new GH1("Phip_585MeVCM9PosHel", "#phi_{p} Distribution for E_{#gamma} 585 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
    Phip595CM9PosHel = new GH1("Phip_595MeVCM9PosHel", "#phi_{p} Distribution for E_{#gamma} 595 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
    Phip605CM9PosHel = new GH1("Phip_605MeVCM9PosHel", "#phi_{p} Distribution for E_{#gamma} 605 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
    Phip615CM9PosHel = new GH1("Phip_615MeVCM9PosHel", "#phi_{p} Distribution for E_{#gamma} 615 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);

    // #phi_{p} dists across EGamma bins
    Phip425CM10PosHel = new GH1("Phip_425MeVCM10PosHel", "#phi_{p} Distribution for E_{#gamma} 425 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
    Phip435CM10PosHel = new GH1("Phip_435MeVCM10PosHel", "#phi_{p} Distribution for E_{#gamma} 435 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
    Phip445CM10PosHel = new GH1("Phip_445MeVCM10PosHel", "#phi_{p} Distribution for E_{#gamma} 445 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
    Phip455CM10PosHel = new GH1("Phip_455MeVCM10PosHel", "#phi_{p} Distribution for E_{#gamma} 455 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
    Phip465CM10PosHel = new GH1("Phip_465MeVCM10PosHel", "#phi_{p} Distribution for E_{#gamma} 465 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
    Phip475CM10PosHel = new GH1("Phip_475MeVCM10PosHel", "#phi_{p} Distribution for E_{#gamma} 475 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
    Phip485CM10PosHel = new GH1("Phip_485MeVCM10PosHel", "#phi_{p} Distribution for E_{#gamma} 485 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
    Phip495CM10PosHel = new GH1("Phip_495MeVCM10PosHel", "#phi_{p} Distribution for E_{#gamma} 495 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
    Phip505CM10PosHel = new GH1("Phip_505MeVCM10PosHel", "#phi_{p} Distribution for E_{#gamma} 505 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
    Phip515CM10PosHel = new GH1("Phip_515MeVCM10PosHel", "#phi_{p} Distribution for E_{#gamma} 515 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
    Phip525CM10PosHel = new GH1("Phip_525MeVCM10PosHel", "#phi_{p} Distribution for E_{#gamma} 525 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
    Phip535CM10PosHel = new GH1("Phip_535MeVCM10PosHel", "#phi_{p} Distribution for E_{#gamma} 535 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
    Phip545CM10PosHel = new GH1("Phip_545MeVCM10PosHel", "#phi_{p} Distribution for E_{#gamma} 545 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
    Phip555CM10PosHel = new GH1("Phip_555MeVCM10PosHel", "#phi_{p} Distribution for E_{#gamma} 555 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
    Phip565CM10PosHel = new GH1("Phip_565MeVCM10PosHel", "#phi_{p} Distribution for E_{#gamma} 565 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
    Phip575CM10PosHel = new GH1("Phip_575MeVCM10PosHel", "#phi_{p} Distribution for E_{#gamma} 575 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
    Phip585CM10PosHel = new GH1("Phip_585MeVCM10PosHel", "#phi_{p} Distribution for E_{#gamma} 585 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
    Phip595CM10PosHel = new GH1("Phip_595MeVCM10PosHel", "#phi_{p} Distribution for E_{#gamma} 595 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
    Phip605CM10PosHel = new GH1("Phip_605MeVCM10PosHel", "#phi_{p} Distribution for E_{#gamma} 605 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
    Phip615CM10PosHel = new GH1("Phip_615MeVCM10PosHel", "#phi_{p} Distribution for E_{#gamma} 615 #pm 5MeV +ve Hel (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);

    Phip425CM1NegHel = new GH1("Phip_425MeVCM1NegHel", "#phi_{p} Distribution for E_{#gamma} 425 #pm 5MeV -ve Hel (Cos#theta_{CM}1-0.8)", 10, -180, 180);
    Phip435CM1NegHel = new GH1("Phip_435MeVCM1NegHel", "#phi_{p} Distribution for E_{#gamma} 435 #pm 5MeV -ve Hel (Cos#theta_{CM}1-0.8)", 10, -180, 180);
    Phip445CM1NegHel = new GH1("Phip_445MeVCM1NegHel", "#phi_{p} Distribution for E_{#gamma} 445 #pm 5MeV -ve Hel (Cos#theta_{CM}1-0.8)", 10, -180, 180);
    Phip455CM1NegHel = new GH1("Phip_455MeVCM1NegHel", "#phi_{p} Distribution for E_{#gamma} 455 #pm 5MeV -ve Hel (Cos#theta_{CM}1-0.8)", 10, -180, 180);
    Phip465CM1NegHel = new GH1("Phip_465MeVCM1NegHel", "#phi_{p} Distribution for E_{#gamma} 465 #pm 5MeV -ve Hel (Cos#theta_{CM}1-0.8)", 10, -180, 180);
    Phip475CM1NegHel = new GH1("Phip_475MeVCM1NegHel", "#phi_{p} Distribution for E_{#gamma} 475 #pm 5MeV -ve Hel (Cos#theta_{CM}1-0.8)", 10, -180, 180);
    Phip485CM1NegHel = new GH1("Phip_485MeVCM1NegHel", "#phi_{p} Distribution for E_{#gamma} 485 #pm 5MeV -ve Hel (Cos#theta_{CM}1-0.8)", 10, -180, 180);
    Phip495CM1NegHel = new GH1("Phip_495MeVCM1NegHel", "#phi_{p} Distribution for E_{#gamma} 495 #pm 5MeV -ve Hel (Cos#theta_{CM}1-0.8)", 10, -180, 180);
    Phip505CM1NegHel = new GH1("Phip_505MeVCM1NegHel", "#phi_{p} Distribution for E_{#gamma} 505 #pm 5MeV -ve Hel (Cos#theta_{CM}1-0.8)", 10, -180, 180);
    Phip515CM1NegHel = new GH1("Phip_515MeVCM1NegHel", "#phi_{p} Distribution for E_{#gamma} 515 #pm 5MeV -ve Hel (Cos#theta_{CM}1-0.8)", 10, -180, 180);
    Phip525CM1NegHel = new GH1("Phip_525MeVCM1NegHel", "#phi_{p} Distribution for E_{#gamma} 525 #pm 5MeV -ve Hel (Cos#theta_{CM}1-0.8)", 10, -180, 180);
    Phip535CM1NegHel = new GH1("Phip_535MeVCM1NegHel", "#phi_{p} Distribution for E_{#gamma} 535 #pm 5MeV -ve Hel (Cos#theta_{CM}1-0.8)", 10, -180, 180);
    Phip545CM1NegHel = new GH1("Phip_545MeVCM1NegHel", "#phi_{p} Distribution for E_{#gamma} 545 #pm 5MeV -ve Hel (Cos#theta_{CM}1-0.8)", 10, -180, 180);
    Phip555CM1NegHel = new GH1("Phip_555MeVCM1NegHel", "#phi_{p} Distribution for E_{#gamma} 555 #pm 5MeV -ve Hel (Cos#theta_{CM}1-0.8)", 10, -180, 180);
    Phip565CM1NegHel = new GH1("Phip_565MeVCM1NegHel", "#phi_{p} Distribution for E_{#gamma} 565 #pm 5MeV -ve Hel (Cos#theta_{CM}1-0.8)", 10, -180, 180);
    Phip575CM1NegHel = new GH1("Phip_575MeVCM1NegHel", "#phi_{p} Distribution for E_{#gamma} 575 #pm 5MeV -ve Hel (Cos#theta_{CM}1-0.8)", 10, -180, 180);
    Phip585CM1NegHel = new GH1("Phip_585MeVCM1NegHel", "#phi_{p} Distribution for E_{#gamma} 585 #pm 5MeV -ve Hel (Cos#theta_{CM}1-0.8)", 10, -180, 180);
    Phip595CM1NegHel = new GH1("Phip_595MeVCM1NegHel", "#phi_{p} Distribution for E_{#gamma} 595 #pm 5MeV -ve Hel (Cos#theta_{CM}1-0.8)", 10, -180, 180);
    Phip605CM1NegHel = new GH1("Phip_605MeVCM1NegHel", "#phi_{p} Distribution for E_{#gamma} 605 #pm 5MeV -ve Hel (Cos#theta_{CM}1-0.8)", 10, -180, 180);
    Phip615CM1NegHel = new GH1("Phip_615MeVCM1NegHel", "#phi_{p} Distribution for E_{#gamma} 615 #pm 5MeV -ve Hel (Cos#theta_{CM}1-0.8)", 10, -180, 180);

    // #phi_{p} dists across EGamma bins
    Phip425CM2NegHel = new GH1("Phip_425MeVCM2NegHel", "#phi_{p} Distribution for E_{#gamma} 425 #pm 5MeV -ve Hel (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
    Phip435CM2NegHel = new GH1("Phip_435MeVCM2NegHel", "#phi_{p} Distribution for E_{#gamma} 435 #pm 5MeV -ve Hel (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
    Phip445CM2NegHel = new GH1("Phip_445MeVCM2NegHel", "#phi_{p} Distribution for E_{#gamma} 445 #pm 5MeV -ve Hel (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
    Phip455CM2NegHel = new GH1("Phip_455MeVCM2NegHel", "#phi_{p} Distribution for E_{#gamma} 455 #pm 5MeV -ve Hel (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
    Phip465CM2NegHel = new GH1("Phip_465MeVCM2NegHel", "#phi_{p} Distribution for E_{#gamma} 465 #pm 5MeV -ve Hel (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
    Phip475CM2NegHel = new GH1("Phip_475MeVCM2NegHel", "#phi_{p} Distribution for E_{#gamma} 475 #pm 5MeV -ve Hel (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
    Phip485CM2NegHel = new GH1("Phip_485MeVCM2NegHel", "#phi_{p} Distribution for E_{#gamma} 485 #pm 5MeV -ve Hel (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
    Phip495CM2NegHel = new GH1("Phip_495MeVCM2NegHel", "#phi_{p} Distribution for E_{#gamma} 495 #pm 5MeV -ve Hel (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
    Phip505CM2NegHel = new GH1("Phip_505MeVCM2NegHel", "#phi_{p} Distribution for E_{#gamma} 505 #pm 5MeV -ve Hel (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
    Phip515CM2NegHel = new GH1("Phip_515MeVCM2NegHel", "#phi_{p} Distribution for E_{#gamma} 515 #pm 5MeV -ve Hel (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
    Phip525CM2NegHel = new GH1("Phip_525MeVCM2NegHel", "#phi_{p} Distribution for E_{#gamma} 525 #pm 5MeV -ve Hel (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
    Phip535CM2NegHel = new GH1("Phip_535MeVCM2NegHel", "#phi_{p} Distribution for E_{#gamma} 535 #pm 5MeV -ve Hel (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
    Phip545CM2NegHel = new GH1("Phip_545MeVCM2NegHel", "#phi_{p} Distribution for E_{#gamma} 545 #pm 5MeV -ve Hel (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
    Phip555CM2NegHel = new GH1("Phip_555MeVCM2NegHel", "#phi_{p} Distribution for E_{#gamma} 555 #pm 5MeV -ve Hel (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
    Phip565CM2NegHel = new GH1("Phip_565MeVCM2NegHel", "#phi_{p} Distribution for E_{#gamma} 565 #pm 5MeV -ve Hel (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
    Phip575CM2NegHel = new GH1("Phip_575MeVCM2NegHel", "#phi_{p} Distribution for E_{#gamma} 575 #pm 5MeV -ve Hel (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
    Phip585CM2NegHel = new GH1("Phip_585MeVCM2NegHel", "#phi_{p} Distribution for E_{#gamma} 585 #pm 5MeV -ve Hel (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
    Phip595CM2NegHel = new GH1("Phip_595MeVCM2NegHel", "#phi_{p} Distribution for E_{#gamma} 595 #pm 5MeV -ve Hel (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
    Phip605CM2NegHel = new GH1("Phip_605MeVCM2NegHel", "#phi_{p} Distribution for E_{#gamma} 605 #pm 5MeV -ve Hel (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
    Phip615CM2NegHel = new GH1("Phip_615MeVCM2NegHel", "#phi_{p} Distribution for E_{#gamma} 615 #pm 5MeV -ve Hel (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);

    // #phi_{p} dists across EGamma bins
    Phip425CM3NegHel = new GH1("Phip_425MeVCM3NegHel", "#phi_{p} Distribution for E_{#gamma} 425 #pm 5MeV -ve Hel (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
    Phip435CM3NegHel = new GH1("Phip_435MeVCM3NegHel", "#phi_{p} Distribution for E_{#gamma} 435 #pm 5MeV -ve Hel (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
    Phip445CM3NegHel = new GH1("Phip_445MeVCM3NegHel", "#phi_{p} Distribution for E_{#gamma} 445 #pm 5MeV -ve Hel (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
    Phip455CM3NegHel = new GH1("Phip_455MeVCM3NegHel", "#phi_{p} Distribution for E_{#gamma} 455 #pm 5MeV -ve Hel (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
    Phip465CM3NegHel = new GH1("Phip_465MeVCM3NegHel", "#phi_{p} Distribution for E_{#gamma} 465 #pm 5MeV -ve Hel (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
    Phip475CM3NegHel = new GH1("Phip_475MeVCM3NegHel", "#phi_{p} Distribution for E_{#gamma} 475 #pm 5MeV -ve Hel (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
    Phip485CM3NegHel = new GH1("Phip_485MeVCM3NegHel", "#phi_{p} Distribution for E_{#gamma} 485 #pm 5MeV -ve Hel (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
    Phip495CM3NegHel = new GH1("Phip_495MeVCM3NegHel", "#phi_{p} Distribution for E_{#gamma} 495 #pm 5MeV -ve Hel (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
    Phip505CM3NegHel = new GH1("Phip_505MeVCM3NegHel", "#phi_{p} Distribution for E_{#gamma} 505 #pm 5MeV -ve Hel (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
    Phip515CM3NegHel = new GH1("Phip_515MeVCM3NegHel", "#phi_{p} Distribution for E_{#gamma} 515 #pm 5MeV -ve Hel (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
    Phip525CM3NegHel = new GH1("Phip_525MeVCM3NegHel", "#phi_{p} Distribution for E_{#gamma} 525 #pm 5MeV -ve Hel (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
    Phip535CM3NegHel = new GH1("Phip_535MeVCM3NegHel", "#phi_{p} Distribution for E_{#gamma} 535 #pm 5MeV -ve Hel (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
    Phip545CM3NegHel = new GH1("Phip_545MeVCM3NegHel", "#phi_{p} Distribution for E_{#gamma} 545 #pm 5MeV -ve Hel (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
    Phip555CM3NegHel = new GH1("Phip_555MeVCM3NegHel", "#phi_{p} Distribution for E_{#gamma} 555 #pm 5MeV -ve Hel (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
    Phip565CM3NegHel = new GH1("Phip_565MeVCM3NegHel", "#phi_{p} Distribution for E_{#gamma} 565 #pm 5MeV -ve Hel (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
    Phip575CM3NegHel = new GH1("Phip_575MeVCM3NegHel", "#phi_{p} Distribution for E_{#gamma} 575 #pm 5MeV -ve Hel (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
    Phip585CM3NegHel = new GH1("Phip_585MeVCM3NegHel", "#phi_{p} Distribution for E_{#gamma} 585 #pm 5MeV -ve Hel (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
    Phip595CM3NegHel = new GH1("Phip_595MeVCM3NegHel", "#phi_{p} Distribution for E_{#gamma} 595 #pm 5MeV -ve Hel (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
    Phip605CM3NegHel = new GH1("Phip_605MeVCM3NegHel", "#phi_{p} Distribution for E_{#gamma} 605 #pm 5MeV -ve Hel (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
    Phip615CM3NegHel = new GH1("Phip_615MeVCM3NegHel", "#phi_{p} Distribution for E_{#gamma} 615 #pm 5MeV -ve Hel (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);

    // #phi_{p} dists across EGamma bins
    Phip425CM4NegHel = new GH1("Phip_425MeVCM4NegHel", "#phi_{p} Distribution for E_{#gamma} 425 #pm 5MeV -ve Hel (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
    Phip435CM4NegHel = new GH1("Phip_435MeVCM4NegHel", "#phi_{p} Distribution for E_{#gamma} 435 #pm 5MeV -ve Hel (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
    Phip445CM4NegHel = new GH1("Phip_445MeVCM4NegHel", "#phi_{p} Distribution for E_{#gamma} 445 #pm 5MeV -ve Hel (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
    Phip455CM4NegHel = new GH1("Phip_455MeVCM4NegHel", "#phi_{p} Distribution for E_{#gamma} 455 #pm 5MeV -ve Hel (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
    Phip465CM4NegHel = new GH1("Phip_465MeVCM4NegHel", "#phi_{p} Distribution for E_{#gamma} 465 #pm 5MeV -ve Hel (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
    Phip475CM4NegHel = new GH1("Phip_475MeVCM4NegHel", "#phi_{p} Distribution for E_{#gamma} 475 #pm 5MeV -ve Hel (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
    Phip485CM4NegHel = new GH1("Phip_485MeVCM4NegHel", "#phi_{p} Distribution for E_{#gamma} 485 #pm 5MeV -ve Hel (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
    Phip495CM4NegHel = new GH1("Phip_495MeVCM4NegHel", "#phi_{p} Distribution for E_{#gamma} 495 #pm 5MeV -ve Hel (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
    Phip505CM4NegHel = new GH1("Phip_505MeVCM4NegHel", "#phi_{p} Distribution for E_{#gamma} 505 #pm 5MeV -ve Hel (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
    Phip515CM4NegHel = new GH1("Phip_515MeVCM4NegHel", "#phi_{p} Distribution for E_{#gamma} 515 #pm 5MeV -ve Hel (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
    Phip525CM4NegHel = new GH1("Phip_525MeVCM4NegHel", "#phi_{p} Distribution for E_{#gamma} 525 #pm 5MeV -ve Hel (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
    Phip535CM4NegHel = new GH1("Phip_535MeVCM4NegHel", "#phi_{p} Distribution for E_{#gamma} 535 #pm 5MeV -ve Hel (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
    Phip545CM4NegHel = new GH1("Phip_545MeVCM4NegHel", "#phi_{p} Distribution for E_{#gamma} 545 #pm 5MeV -ve Hel (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
    Phip555CM4NegHel = new GH1("Phip_555MeVCM4NegHel", "#phi_{p} Distribution for E_{#gamma} 555 #pm 5MeV -ve Hel (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
    Phip565CM4NegHel = new GH1("Phip_565MeVCM4NegHel", "#phi_{p} Distribution for E_{#gamma} 565 #pm 5MeV -ve Hel (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
    Phip575CM4NegHel = new GH1("Phip_575MeVCM4NegHel", "#phi_{p} Distribution for E_{#gamma} 575 #pm 5MeV -ve Hel (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
    Phip585CM4NegHel = new GH1("Phip_585MeVCM4NegHel", "#phi_{p} Distribution for E_{#gamma} 585 #pm 5MeV -ve Hel (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
    Phip595CM4NegHel = new GH1("Phip_595MeVCM4NegHel", "#phi_{p} Distribution for E_{#gamma} 595 #pm 5MeV -ve Hel (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
    Phip605CM4NegHel = new GH1("Phip_605MeVCM4NegHel", "#phi_{p} Distribution for E_{#gamma} 605 #pm 5MeV -ve Hel (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
    Phip615CM4NegHel = new GH1("Phip_615MeVCM4NegHel", "#phi_{p} Distribution for E_{#gamma} 615 #pm 5MeV -ve Hel (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);

    // #phi_{p} dists across EGamma bins
    Phip425CM5NegHel = new GH1("Phip_425MeVCM5NegHel", "#phi_{p} Distribution for E_{#gamma} 425 #pm 5MeV -ve Hel (Cos#theta_{CM}0.2-0)", 10, -180, 180);
    Phip435CM5NegHel = new GH1("Phip_435MeVCM5NegHel", "#phi_{p} Distribution for E_{#gamma} 435 #pm 5MeV -ve Hel (Cos#theta_{CM}0.2-0)", 10, -180, 180);
    Phip445CM5NegHel = new GH1("Phip_445MeVCM5NegHel", "#phi_{p} Distribution for E_{#gamma} 445 #pm 5MeV -ve Hel (Cos#theta_{CM}0.2-0)", 10, -180, 180);
    Phip455CM5NegHel = new GH1("Phip_455MeVCM5NegHel", "#phi_{p} Distribution for E_{#gamma} 455 #pm 5MeV -ve Hel (Cos#theta_{CM}0.2-0)", 10, -180, 180);
    Phip465CM5NegHel = new GH1("Phip_465MeVCM5NegHel", "#phi_{p} Distribution for E_{#gamma} 465 #pm 5MeV -ve Hel (Cos#theta_{CM}0.2-0)", 10, -180, 180);
    Phip475CM5NegHel = new GH1("Phip_475MeVCM5NegHel", "#phi_{p} Distribution for E_{#gamma} 475 #pm 5MeV -ve Hel (Cos#theta_{CM}0.2-0)", 10, -180, 180);
    Phip485CM5NegHel = new GH1("Phip_485MeVCM5NegHel", "#phi_{p} Distribution for E_{#gamma} 485 #pm 5MeV -ve Hel (Cos#theta_{CM}0.2-0)", 10, -180, 180);
    Phip495CM5NegHel = new GH1("Phip_495MeVCM5NegHel", "#phi_{p} Distribution for E_{#gamma} 495 #pm 5MeV -ve Hel (Cos#theta_{CM}0.2-0)", 10, -180, 180);
    Phip505CM5NegHel = new GH1("Phip_505MeVCM5NegHel", "#phi_{p} Distribution for E_{#gamma} 505 #pm 5MeV -ve Hel (Cos#theta_{CM}0.2-0)", 10, -180, 180);
    Phip515CM5NegHel = new GH1("Phip_515MeVCM5NegHel", "#phi_{p} Distribution for E_{#gamma} 515 #pm 5MeV -ve Hel (Cos#theta_{CM}0.2-0)", 10, -180, 180);
    Phip525CM5NegHel = new GH1("Phip_525MeVCM5NegHel", "#phi_{p} Distribution for E_{#gamma} 525 #pm 5MeV -ve Hel (Cos#theta_{CM}0.2-0)", 10, -180, 180);
    Phip535CM5NegHel = new GH1("Phip_535MeVCM5NegHel", "#phi_{p} Distribution for E_{#gamma} 535 #pm 5MeV -ve Hel (Cos#theta_{CM}0.2-0)", 10, -180, 180);
    Phip545CM5NegHel = new GH1("Phip_545MeVCM5NegHel", "#phi_{p} Distribution for E_{#gamma} 545 #pm 5MeV -ve Hel (Cos#theta_{CM}0.2-0)", 10, -180, 180);
    Phip555CM5NegHel = new GH1("Phip_555MeVCM5NegHel", "#phi_{p} Distribution for E_{#gamma} 555 #pm 5MeV -ve Hel (Cos#theta_{CM}0.2-0)", 10, -180, 180);
    Phip565CM5NegHel = new GH1("Phip_565MeVCM5NegHel", "#phi_{p} Distribution for E_{#gamma} 565 #pm 5MeV -ve Hel (Cos#theta_{CM}0.2-0)", 10, -180, 180);
    Phip575CM5NegHel = new GH1("Phip_575MeVCM5NegHel", "#phi_{p} Distribution for E_{#gamma} 575 #pm 5MeV -ve Hel (Cos#theta_{CM}0.2-0)", 10, -180, 180);
    Phip585CM5NegHel = new GH1("Phip_585MeVCM5NegHel", "#phi_{p} Distribution for E_{#gamma} 585 #pm 5MeV -ve Hel (Cos#theta_{CM}0.2-0)", 10, -180, 180);
    Phip595CM5NegHel = new GH1("Phip_595MeVCM5NegHel", "#phi_{p} Distribution for E_{#gamma} 595 #pm 5MeV -ve Hel (Cos#theta_{CM}0.2-0)", 10, -180, 180);
    Phip605CM5NegHel = new GH1("Phip_605MeVCM5NegHel", "#phi_{p} Distribution for E_{#gamma} 605 #pm 5MeV -ve Hel (Cos#theta_{CM}0.2-0)", 10, -180, 180);
    Phip615CM5NegHel = new GH1("Phip_615MeVCM5NegHel", "#phi_{p} Distribution for E_{#gamma} 615 #pm 5MeV -ve Hel (Cos#theta_{CM}0.2-0)", 10, -180, 180);

    // #phi_{p} dists across EGamma bins
    Phip425CM6NegHel = new GH1("Phip_425MeVCM6NegHel", "#phi_{p} Distribution for E_{#gamma} 425 #pm 5MeV -ve Hel (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
    Phip435CM6NegHel = new GH1("Phip_435MeVCM6NegHel", "#phi_{p} Distribution for E_{#gamma} 435 #pm 5MeV -ve Hel (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
    Phip445CM6NegHel = new GH1("Phip_445MeVCM6NegHel", "#phi_{p} Distribution for E_{#gamma} 445 #pm 5MeV -ve Hel (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
    Phip455CM6NegHel = new GH1("Phip_455MeVCM6NegHel", "#phi_{p} Distribution for E_{#gamma} 455 #pm 5MeV -ve Hel (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
    Phip465CM6NegHel = new GH1("Phip_465MeVCM6NegHel", "#phi_{p} Distribution for E_{#gamma} 465 #pm 5MeV -ve Hel (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
    Phip475CM6NegHel = new GH1("Phip_475MeVCM6NegHel", "#phi_{p} Distribution for E_{#gamma} 475 #pm 5MeV -ve Hel (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
    Phip485CM6NegHel = new GH1("Phip_485MeVCM6NegHel", "#phi_{p} Distribution for E_{#gamma} 485 #pm 5MeV -ve Hel (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
    Phip495CM6NegHel = new GH1("Phip_495MeVCM6NegHel", "#phi_{p} Distribution for E_{#gamma} 495 #pm 5MeV -ve Hel (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
    Phip505CM6NegHel = new GH1("Phip_505MeVCM6NegHel", "#phi_{p} Distribution for E_{#gamma} 505 #pm 5MeV -ve Hel (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
    Phip515CM6NegHel = new GH1("Phip_515MeVCM6NegHel", "#phi_{p} Distribution for E_{#gamma} 515 #pm 5MeV -ve Hel (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
    Phip525CM6NegHel = new GH1("Phip_525MeVCM6NegHel", "#phi_{p} Distribution for E_{#gamma} 525 #pm 5MeV -ve Hel (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
    Phip535CM6NegHel = new GH1("Phip_535MeVCM6NegHel", "#phi_{p} Distribution for E_{#gamma} 535 #pm 5MeV -ve Hel (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
    Phip545CM6NegHel = new GH1("Phip_545MeVCM6NegHel", "#phi_{p} Distribution for E_{#gamma} 545 #pm 5MeV -ve Hel (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
    Phip555CM6NegHel = new GH1("Phip_555MeVCM6NegHel", "#phi_{p} Distribution for E_{#gamma} 555 #pm 5MeV -ve Hel (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
    Phip565CM6NegHel = new GH1("Phip_565MeVCM6NegHel", "#phi_{p} Distribution for E_{#gamma} 565 #pm 5MeV -ve Hel (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
    Phip575CM6NegHel = new GH1("Phip_575MeVCM6NegHel", "#phi_{p} Distribution for E_{#gamma} 575 #pm 5MeV -ve Hel (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
    Phip585CM6NegHel = new GH1("Phip_585MeVCM6NegHel", "#phi_{p} Distribution for E_{#gamma} 585 #pm 5MeV -ve Hel (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
    Phip595CM6NegHel = new GH1("Phip_595MeVCM6NegHel", "#phi_{p} Distribution for E_{#gamma} 595 #pm 5MeV -ve Hel (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
    Phip605CM6NegHel = new GH1("Phip_605MeVCM6NegHel", "#phi_{p} Distribution for E_{#gamma} 605 #pm 5MeV -ve Hel (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
    Phip615CM6NegHel = new GH1("Phip_615MeVCM6NegHel", "#phi_{p} Distribution for E_{#gamma} 615 #pm 5MeV -ve Hel (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);

    // #phi_{p} dists across EGamma bins
    Phip425CM7NegHel = new GH1("Phip_425MeVCM7NegHel", "#phi_{p} Distribution for E_{#gamma} 425 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
    Phip435CM7NegHel = new GH1("Phip_435MeVCM7NegHel", "#phi_{p} Distribution for E_{#gamma} 435 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
    Phip445CM7NegHel = new GH1("Phip_445MeVCM7NegHel", "#phi_{p} Distribution for E_{#gamma} 445 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
    Phip455CM7NegHel = new GH1("Phip_455MeVCM7NegHel", "#phi_{p} Distribution for E_{#gamma} 455 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
    Phip465CM7NegHel = new GH1("Phip_465MeVCM7NegHel", "#phi_{p} Distribution for E_{#gamma} 465 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
    Phip475CM7NegHel = new GH1("Phip_475MeVCM7NegHel", "#phi_{p} Distribution for E_{#gamma} 475 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
    Phip485CM7NegHel = new GH1("Phip_485MeVCM7NegHel", "#phi_{p} Distribution for E_{#gamma} 485 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
    Phip495CM7NegHel = new GH1("Phip_495MeVCM7NegHel", "#phi_{p} Distribution for E_{#gamma} 495 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
    Phip505CM7NegHel = new GH1("Phip_505MeVCM7NegHel", "#phi_{p} Distribution for E_{#gamma} 505 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
    Phip515CM7NegHel = new GH1("Phip_515MeVCM7NegHel", "#phi_{p} Distribution for E_{#gamma} 515 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
    Phip525CM7NegHel = new GH1("Phip_525MeVCM7NegHel", "#phi_{p} Distribution for E_{#gamma} 525 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
    Phip535CM7NegHel = new GH1("Phip_535MeVCM7NegHel", "#phi_{p} Distribution for E_{#gamma} 535 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
    Phip545CM7NegHel = new GH1("Phip_545MeVCM7NegHel", "#phi_{p} Distribution for E_{#gamma} 545 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
    Phip555CM7NegHel = new GH1("Phip_555MeVCM7NegHel", "#phi_{p} Distribution for E_{#gamma} 555 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
    Phip565CM7NegHel = new GH1("Phip_565MeVCM7NegHel", "#phi_{p} Distribution for E_{#gamma} 565 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
    Phip575CM7NegHel = new GH1("Phip_575MeVCM7NegHel", "#phi_{p} Distribution for E_{#gamma} 575 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
    Phip585CM7NegHel = new GH1("Phip_585MeVCM7NegHel", "#phi_{p} Distribution for E_{#gamma} 585 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
    Phip595CM7NegHel = new GH1("Phip_595MeVCM7NegHel", "#phi_{p} Distribution for E_{#gamma} 595 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
    Phip605CM7NegHel = new GH1("Phip_605MeVCM7NegHel", "#phi_{p} Distribution for E_{#gamma} 605 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
    Phip615CM7NegHel = new GH1("Phip_615MeVCM7NegHel", "#phi_{p} Distribution for E_{#gamma} 615 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);

    // #phi_{p} dists across EGamma bins
    Phip425CM8NegHel = new GH1("Phip_425MeVCM8NegHel", "#phi_{p} Distribution for E_{#gamma} 425 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
    Phip435CM8NegHel = new GH1("Phip_435MeVCM8NegHel", "#phi_{p} Distribution for E_{#gamma} 435 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
    Phip445CM8NegHel = new GH1("Phip_445MeVCM8NegHel", "#phi_{p} Distribution for E_{#gamma} 445 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
    Phip455CM8NegHel = new GH1("Phip_455MeVCM8NegHel", "#phi_{p} Distribution for E_{#gamma} 455 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
    Phip465CM8NegHel = new GH1("Phip_465MeVCM8NegHel", "#phi_{p} Distribution for E_{#gamma} 465 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
    Phip475CM8NegHel = new GH1("Phip_475MeVCM8NegHel", "#phi_{p} Distribution for E_{#gamma} 475 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
    Phip485CM8NegHel = new GH1("Phip_485MeVCM8NegHel", "#phi_{p} Distribution for E_{#gamma} 485 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
    Phip495CM8NegHel = new GH1("Phip_495MeVCM8NegHel", "#phi_{p} Distribution for E_{#gamma} 495 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
    Phip505CM8NegHel = new GH1("Phip_505MeVCM8NegHel", "#phi_{p} Distribution for E_{#gamma} 505 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
    Phip515CM8NegHel = new GH1("Phip_515MeVCM8NegHel", "#phi_{p} Distribution for E_{#gamma} 515 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
    Phip525CM8NegHel = new GH1("Phip_525MeVCM8NegHel", "#phi_{p} Distribution for E_{#gamma} 525 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
    Phip535CM8NegHel = new GH1("Phip_535MeVCM8NegHel", "#phi_{p} Distribution for E_{#gamma} 535 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
    Phip545CM8NegHel = new GH1("Phip_545MeVCM8NegHel", "#phi_{p} Distribution for E_{#gamma} 545 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
    Phip555CM8NegHel = new GH1("Phip_555MeVCM8NegHel", "#phi_{p} Distribution for E_{#gamma} 555 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
    Phip565CM8NegHel = new GH1("Phip_565MeVCM8NegHel", "#phi_{p} Distribution for E_{#gamma} 565 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
    Phip575CM8NegHel = new GH1("Phip_575MeVCM8NegHel", "#phi_{p} Distribution for E_{#gamma} 575 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
    Phip585CM8NegHel = new GH1("Phip_585MeVCM8NegHel", "#phi_{p} Distribution for E_{#gamma} 585 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
    Phip595CM8NegHel = new GH1("Phip_595MeVCM8NegHel", "#phi_{p} Distribution for E_{#gamma} 595 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
    Phip605CM8NegHel = new GH1("Phip_605MeVCM8NegHel", "#phi_{p} Distribution for E_{#gamma} 605 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
    Phip615CM8NegHel = new GH1("Phip_615MeVCM8NegHel", "#phi_{p} Distribution for E_{#gamma} 615 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);

    // #phi_{p} dists across EGamma bins
    Phip425CM9NegHel = new GH1("Phip_425MeVCM9NegHel", "#phi_{p} Distribution for E_{#gamma} 425 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
    Phip435CM9NegHel = new GH1("Phip_435MeVCM9NegHel", "#phi_{p} Distribution for E_{#gamma} 435 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
    Phip445CM9NegHel = new GH1("Phip_445MeVCM9NegHel", "#phi_{p} Distribution for E_{#gamma} 445 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
    Phip455CM9NegHel = new GH1("Phip_455MeVCM9NegHel", "#phi_{p} Distribution for E_{#gamma} 455 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
    Phip465CM9NegHel = new GH1("Phip_465MeVCM9NegHel", "#phi_{p} Distribution for E_{#gamma} 465 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
    Phip475CM9NegHel = new GH1("Phip_475MeVCM9NegHel", "#phi_{p} Distribution for E_{#gamma} 475 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
    Phip485CM9NegHel = new GH1("Phip_485MeVCM9NegHel", "#phi_{p} Distribution for E_{#gamma} 485 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
    Phip495CM9NegHel = new GH1("Phip_495MeVCM9NegHel", "#phi_{p} Distribution for E_{#gamma} 495 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
    Phip505CM9NegHel = new GH1("Phip_505MeVCM9NegHel", "#phi_{p} Distribution for E_{#gamma} 505 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
    Phip515CM9NegHel = new GH1("Phip_515MeVCM9NegHel", "#phi_{p} Distribution for E_{#gamma} 515 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
    Phip525CM9NegHel = new GH1("Phip_525MeVCM9NegHel", "#phi_{p} Distribution for E_{#gamma} 525 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
    Phip535CM9NegHel = new GH1("Phip_535MeVCM9NegHel", "#phi_{p} Distribution for E_{#gamma} 535 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
    Phip545CM9NegHel = new GH1("Phip_545MeVCM9NegHel", "#phi_{p} Distribution for E_{#gamma} 545 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
    Phip555CM9NegHel = new GH1("Phip_555MeVCM9NegHel", "#phi_{p} Distribution for E_{#gamma} 555 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
    Phip565CM9NegHel = new GH1("Phip_565MeVCM9NegHel", "#phi_{p} Distribution for E_{#gamma} 565 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
    Phip575CM9NegHel = new GH1("Phip_575MeVCM9NegHel", "#phi_{p} Distribution for E_{#gamma} 575 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
    Phip585CM9NegHel = new GH1("Phip_585MeVCM9NegHel", "#phi_{p} Distribution for E_{#gamma} 585 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
    Phip595CM9NegHel = new GH1("Phip_595MeVCM9NegHel", "#phi_{p} Distribution for E_{#gamma} 595 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
    Phip605CM9NegHel = new GH1("Phip_605MeVCM9NegHel", "#phi_{p} Distribution for E_{#gamma} 605 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
    Phip615CM9NegHel = new GH1("Phip_615MeVCM9NegHel", "#phi_{p} Distribution for E_{#gamma} 615 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);

    // #phi_{p} dists across EGamma bins
    Phip425CM10NegHel = new GH1("Phip_425MeVCM10NegHel", "#phi_{p} Distribution for E_{#gamma} 425 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
    Phip435CM10NegHel = new GH1("Phip_435MeVCM10NegHel", "#phi_{p} Distribution for E_{#gamma} 435 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
    Phip445CM10NegHel = new GH1("Phip_445MeVCM10NegHel", "#phi_{p} Distribution for E_{#gamma} 445 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
    Phip455CM10NegHel = new GH1("Phip_455MeVCM10NegHel", "#phi_{p} Distribution for E_{#gamma} 455 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
    Phip465CM10NegHel = new GH1("Phip_465MeVCM10NegHel", "#phi_{p} Distribution for E_{#gamma} 465 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
    Phip475CM10NegHel = new GH1("Phip_475MeVCM10NegHel", "#phi_{p} Distribution for E_{#gamma} 475 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
    Phip485CM10NegHel = new GH1("Phip_485MeVCM10NegHel", "#phi_{p} Distribution for E_{#gamma} 485 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
    Phip495CM10NegHel = new GH1("Phip_495MeVCM10NegHel", "#phi_{p} Distribution for E_{#gamma} 495 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
    Phip505CM10NegHel = new GH1("Phip_505MeVCM10NegHel", "#phi_{p} Distribution for E_{#gamma} 505 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
    Phip515CM10NegHel = new GH1("Phip_515MeVCM10NegHel", "#phi_{p} Distribution for E_{#gamma} 515 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
    Phip525CM10NegHel = new GH1("Phip_525MeVCM10NegHel", "#phi_{p} Distribution for E_{#gamma} 525 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
    Phip535CM10NegHel = new GH1("Phip_535MeVCM10NegHel", "#phi_{p} Distribution for E_{#gamma} 535 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
    Phip545CM10NegHel = new GH1("Phip_545MeVCM10NegHel", "#phi_{p} Distribution for E_{#gamma} 545 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
    Phip555CM10NegHel = new GH1("Phip_555MeVCM10NegHel", "#phi_{p} Distribution for E_{#gamma} 555 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
    Phip565CM10NegHel = new GH1("Phip_565MeVCM10NegHel", "#phi_{p} Distribution for E_{#gamma} 565 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
    Phip575CM10NegHel = new GH1("Phip_575MeVCM10NegHel", "#phi_{p} Distribution for E_{#gamma} 575 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
    Phip585CM10NegHel = new GH1("Phip_585MeVCM10NegHel", "#phi_{p} Distribution for E_{#gamma} 585 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
    Phip595CM10NegHel = new GH1("Phip_595MeVCM10NegHel", "#phi_{p} Distribution for E_{#gamma} 595 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
    Phip605CM10NegHel = new GH1("Phip_605MeVCM10NegHel", "#phi_{p} Distribution for E_{#gamma} 605 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
    Phip615CM10NegHel = new GH1("Phip_615MeVCM10NegHel", "#phi_{p} Distribution for E_{#gamma} 615 #pm 5MeV -ve Hel (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);


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


    if ( 420 < EGamma && EGamma < 430) {

        if(1 > CosThetapCM && CosThetapCM > 0.8 ){
            Phip425CM1->Fill(Php, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.6){
            Phip425CM2->Fill(Php, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.4){
            Phip425CM3->Fill(Php, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.2){
            Phip425CM4->Fill(Php, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0){
            Phip425CM5->Fill(Php, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.2){
            Phip425CM6->Fill(Php, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.4){
            Phip425CM7->Fill(Php, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.6){
            Phip425CM8->Fill(Php, TaggerTime);
        }

        else if(-0.6> CosThetapCM && CosThetapCM > -0.8){
            Phip425CM9->Fill(Php, TaggerTime);
        }

        else if(-0.8> CosThetapCM && CosThetapCM > -1){
            Phip425CM10->Fill(Php, TaggerTime);
        }
    }

    else if ( 430 < EGamma && EGamma < 440) {

        if(1 > CosThetapCM && CosThetapCM > 0.8 ){
            Phip435CM1->Fill(Php, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.6){
            Phip435CM2->Fill(Php, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.4){
            Phip435CM3->Fill(Php, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.2){
            Phip435CM4->Fill(Php, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0){
            Phip435CM5->Fill(Php, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.2){
            Phip435CM6->Fill(Php, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.4){
            Phip435CM7->Fill(Php, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.6){
            Phip435CM8->Fill(Php, TaggerTime);
        }

        else if(-0.6> CosThetapCM && CosThetapCM > -0.8){
            Phip435CM9->Fill(Php, TaggerTime);
        }

        else if(-0.8> CosThetapCM && CosThetapCM > -1){
            Phip435CM10->Fill(Php, TaggerTime);
        }
    }

    else if ( 440 < EGamma && EGamma < 450) {

        if(1 > CosThetapCM && CosThetapCM > 0.8 ){
            Phip445CM1->Fill(Php, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.6){
            Phip445CM2->Fill(Php, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.4){
            Phip445CM3->Fill(Php, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.2){
            Phip445CM4->Fill(Php, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0){
            Phip445CM5->Fill(Php, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.2){
            Phip445CM6->Fill(Php, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.4){
            Phip445CM7->Fill(Php, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.6){
            Phip445CM8->Fill(Php, TaggerTime);
        }

        else if(-0.6> CosThetapCM && CosThetapCM > -0.8){
            Phip445CM9->Fill(Php, TaggerTime);
        }

        else if(-0.8> CosThetapCM && CosThetapCM > -1){
            Phip445CM10->Fill(Php, TaggerTime);
        }
    }

    else if ( 450 < EGamma && EGamma < 460) {

        if(1 > CosThetapCM && CosThetapCM > 0.8 ){
            Phip455CM1->Fill(Php, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.6){
            Phip455CM2->Fill(Php, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.4){
            Phip455CM3->Fill(Php, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.2){
            Phip455CM4->Fill(Php, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0){
            Phip455CM5->Fill(Php, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.2){
            Phip455CM6->Fill(Php, TaggerTime);
        }
        else if(-0.2 > CosThetapCM && CosThetapCM > -0.4){
            Phip455CM7->Fill(Php, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.6){
            Phip455CM8->Fill(Php, TaggerTime);
        }

        else if(-0.6> CosThetapCM && CosThetapCM > -0.8){
            Phip455CM9->Fill(Php, TaggerTime);
        }

        else if(-0.8> CosThetapCM && CosThetapCM > -1){
            Phip455CM10->Fill(Php, TaggerTime);
        }
    }

    else if ( 460 < EGamma && EGamma < 470) {

        if(1 > CosThetapCM && CosThetapCM > 0.8 ){
            Phip465CM1->Fill(Php, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.6){
            Phip465CM2->Fill(Php, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.4){
            Phip465CM3->Fill(Php, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.2){
            Phip465CM4->Fill(Php, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0){
            Phip465CM5->Fill(Php, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.2){
            Phip465CM6->Fill(Php, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.4){
            Phip465CM7->Fill(Php, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.6){
            Phip465CM8->Fill(Php, TaggerTime);
        }

        else if(-0.6> CosThetapCM && CosThetapCM > -0.8){
            Phip465CM9->Fill(Php, TaggerTime);
        }

        else if(-0.8> CosThetapCM && CosThetapCM > -1){
            Phip465CM10->Fill(Php, TaggerTime);
        }
    }

    else if ( 470 < EGamma && EGamma < 480) {

        if(1 > CosThetapCM && CosThetapCM > 0.8 ){
            Phip475CM1->Fill(Php, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.6){
            Phip475CM2->Fill(Php, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.4){
            Phip475CM3->Fill(Php, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.2){
            Phip475CM4->Fill(Php, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0){
            Phip475CM5->Fill(Php, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.2){
            Phip475CM6->Fill(Php, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.4){
            Phip475CM7->Fill(Php, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.6){
            Phip475CM8->Fill(Php, TaggerTime);
        }

        else if(-0.6> CosThetapCM && CosThetapCM > -0.8){
            Phip475CM9->Fill(Php, TaggerTime);
        }

        else if(-0.8> CosThetapCM && CosThetapCM > -1){
            Phip475CM10->Fill(Php, TaggerTime);
        }
    }

    else if ( 480 < EGamma && EGamma < 490) {

        if(1 > CosThetapCM && CosThetapCM > 0.8 ){
            Phip485CM1->Fill(Php, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.6){
            Phip485CM2->Fill(Php, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.4){
            Phip485CM3->Fill(Php, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.2){
            Phip485CM4->Fill(Php, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0){
            Phip485CM5->Fill(Php, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.2){
            Phip485CM6->Fill(Php, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.4){
            Phip485CM7->Fill(Php, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.6){
            Phip485CM8->Fill(Php, TaggerTime);
        }

        else if(-0.6> CosThetapCM && CosThetapCM > -0.8){
            Phip485CM9->Fill(Php, TaggerTime);
        }

        else if(-0.8> CosThetapCM && CosThetapCM > -1){
            Phip485CM10->Fill(Php, TaggerTime);
        }
    }

    else if ( 490 < EGamma && EGamma < 500) {

        if(1 > CosThetapCM && CosThetapCM > 0.8 ){
            Phip495CM1->Fill(Php, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.6){
            Phip495CM2->Fill(Php, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.4){
            Phip495CM3->Fill(Php, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.2){
            Phip495CM4->Fill(Php, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0){
            Phip495CM5->Fill(Php, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.2){
            Phip495CM6->Fill(Php, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.4){
            Phip495CM7->Fill(Php, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.6){
            Phip495CM8->Fill(Php, TaggerTime);
        }

        else if(-0.6> CosThetapCM && CosThetapCM > -0.8){
            Phip495CM9->Fill(Php, TaggerTime);
        }

        else if(-0.8> CosThetapCM && CosThetapCM > -1){
            Phip495CM10->Fill(Php, TaggerTime);
        }
    }

    else if ( 500 < EGamma && EGamma < 510) {

        if(1 > CosThetapCM && CosThetapCM > 0.8 ){
            Phip505CM1->Fill(Php, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.6){
            Phip505CM2->Fill(Php, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.4){
            Phip505CM3->Fill(Php, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.2){
            Phip505CM4->Fill(Php, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0){
            Phip505CM5->Fill(Php, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.2){
            Phip505CM6->Fill(Php, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.4){
            Phip505CM7->Fill(Php, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.6){
            Phip505CM8->Fill(Php, TaggerTime);
        }

        else if(-0.6> CosThetapCM && CosThetapCM > -0.8){
            Phip505CM9->Fill(Php, TaggerTime);
        }

        else if(-0.8> CosThetapCM && CosThetapCM > -1){
            Phip505CM10->Fill(Php, TaggerTime);
        }
    }

    else if ( 510 < EGamma && EGamma < 520) {

        if(1 > CosThetapCM && CosThetapCM > 0.8 ){
            Phip515CM1->Fill(Php, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.6){
            Phip515CM2->Fill(Php, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.4){
            Phip515CM3->Fill(Php, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.2){
            Phip515CM4->Fill(Php, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0){
            Phip515CM5->Fill(Php, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.2){
            Phip515CM6->Fill(Php, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.4){
            Phip515CM7->Fill(Php, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.6){
            Phip515CM8->Fill(Php, TaggerTime);
        }

        else if(-0.6> CosThetapCM && CosThetapCM > -0.8){
            Phip515CM9->Fill(Php, TaggerTime);
        }

        else if(-0.8> CosThetapCM && CosThetapCM > -1){
            Phip515CM10->Fill(Php, TaggerTime);
        }
    }

    else if ( 520 < EGamma && EGamma < 530) {

        if(1 > CosThetapCM && CosThetapCM > 0.8 ){
            Phip525CM1->Fill(Php, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.6){
            Phip525CM2->Fill(Php, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.4){
            Phip525CM3->Fill(Php, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.2){
            Phip525CM4->Fill(Php, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0){
            Phip525CM5->Fill(Php, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.2){
            Phip525CM6->Fill(Php, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.4){
            Phip525CM7->Fill(Php, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.6){
            Phip525CM8->Fill(Php, TaggerTime);
        }

        else if(-0.6> CosThetapCM && CosThetapCM > -0.8){
            Phip525CM9->Fill(Php, TaggerTime);
        }

        else if(-0.8> CosThetapCM && CosThetapCM > -1){
            Phip525CM10->Fill(Php, TaggerTime);
        }
    }

    else if ( 530 < EGamma && EGamma < 540) {

        if(1 > CosThetapCM && CosThetapCM > 0.8 ){
            Phip535CM1->Fill(Php, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.6){
            Phip535CM2->Fill(Php, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.4){
            Phip535CM3->Fill(Php, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.2){
            Phip535CM4->Fill(Php, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0){
            Phip535CM5->Fill(Php, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.2){
            Phip535CM6->Fill(Php, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.4){
            Phip535CM7->Fill(Php, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.6){
            Phip535CM8->Fill(Php, TaggerTime);
        }

        else if(-0.6> CosThetapCM && CosThetapCM > -0.8){
            Phip535CM9->Fill(Php, TaggerTime);
        }

        else if(-0.8> CosThetapCM && CosThetapCM > -1){
            Phip535CM10->Fill(Php, TaggerTime);
        }
    }

    else if ( 540 < EGamma && EGamma < 550) {

        if(1 > CosThetapCM && CosThetapCM > 0.8 ){
            Phip545CM1->Fill(Php, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.6){
            Phip545CM2->Fill(Php, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.4){
            Phip545CM3->Fill(Php, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.2){
            Phip545CM4->Fill(Php, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0){
            Phip545CM5->Fill(Php, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.2){
            Phip545CM6->Fill(Php, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.4){
            Phip545CM7->Fill(Php, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.6){
            Phip545CM8->Fill(Php, TaggerTime);
        }

        else if(-0.6> CosThetapCM && CosThetapCM > -0.8){
            Phip545CM9->Fill(Php, TaggerTime);
        }

        else if(-0.8> CosThetapCM && CosThetapCM > -1){
            Phip545CM10->Fill(Php, TaggerTime);
        }
    }

    else if ( 550 < EGamma && EGamma < 560) {

        if(1 > CosThetapCM && CosThetapCM > 0.8 ){
            Phip555CM1->Fill(Php, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.6){
            Phip555CM2->Fill(Php, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.4){
            Phip555CM3->Fill(Php, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.2){
            Phip555CM4->Fill(Php, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0){
            Phip555CM5->Fill(Php, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.2){
            Phip555CM6->Fill(Php, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.4){
            Phip555CM7->Fill(Php, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.6){
            Phip555CM8->Fill(Php, TaggerTime);
        }

        else if(-0.6> CosThetapCM && CosThetapCM > -0.8){
            Phip555CM9->Fill(Php, TaggerTime);
        }

        else if(-0.8> CosThetapCM && CosThetapCM > -1){
            Phip555CM10->Fill(Php, TaggerTime);
        }
    }

    else if ( 560 < EGamma && EGamma < 570) {

        if(1 > CosThetapCM && CosThetapCM > 0.8 ){
            Phip565CM1->Fill(Php, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.6){
            Phip565CM2->Fill(Php, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.4){
            Phip565CM3->Fill(Php, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.2){
            Phip565CM4->Fill(Php, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0){
            Phip565CM5->Fill(Php, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.2){
            Phip565CM6->Fill(Php, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.4){
            Phip565CM7->Fill(Php, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.6){
            Phip565CM8->Fill(Php, TaggerTime);
        }

        else if(-0.6> CosThetapCM && CosThetapCM > -0.8){
            Phip565CM9->Fill(Php, TaggerTime);
        }

        else if(-0.8> CosThetapCM && CosThetapCM > -1){
            Phip565CM10->Fill(Php, TaggerTime);
        }
    }

    else if ( 570 < EGamma && EGamma < 580) {

        if(1 > CosThetapCM && CosThetapCM > 0.8 ){
            Phip575CM1->Fill(Php, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.6){
            Phip575CM2->Fill(Php, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.4){
            Phip575CM3->Fill(Php, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.2){
            Phip575CM4->Fill(Php, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0){
            Phip575CM5->Fill(Php, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.2){
            Phip575CM6->Fill(Php, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.4){
            Phip575CM7->Fill(Php, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.6){
            Phip575CM8->Fill(Php, TaggerTime);
        }

        else if(-0.6> CosThetapCM && CosThetapCM > -0.8){
            Phip575CM9->Fill(Php, TaggerTime);
        }

        else if(-0.8> CosThetapCM && CosThetapCM > -1){
            Phip575CM10->Fill(Php, TaggerTime);
        }
    }

    else if ( 580 < EGamma && EGamma < 590) {

        if(1 > CosThetapCM && CosThetapCM > 0.8 ){
            Phip585CM1->Fill(Php, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.6){
            Phip585CM2->Fill(Php, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.4){
            Phip585CM3->Fill(Php, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.2){
            Phip585CM4->Fill(Php, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0){
            Phip585CM5->Fill(Php, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.2){
            Phip585CM6->Fill(Php, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.4){
            Phip585CM7->Fill(Php, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.6){
            Phip585CM8->Fill(Php, TaggerTime);
        }

        else if(-0.6> CosThetapCM && CosThetapCM > -0.8){
            Phip585CM9->Fill(Php, TaggerTime);
        }

        else if(-0.8> CosThetapCM && CosThetapCM > -1){
            Phip585CM10->Fill(Php, TaggerTime);
        }
    }

    else if ( 590 < EGamma && EGamma < 600) {

        if(1 > CosThetapCM && CosThetapCM > 0.8 ){
            Phip595CM1->Fill(Php, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.6){
            Phip595CM2->Fill(Php, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.4){
            Phip595CM3->Fill(Php, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.2){
            Phip595CM4->Fill(Php, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0){
            Phip595CM5->Fill(Php, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.2){
            Phip595CM6->Fill(Php, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.4){
            Phip595CM7->Fill(Php, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.6){
            Phip595CM8->Fill(Php, TaggerTime);
        }

        else if(-0.6> CosThetapCM && CosThetapCM > -0.8){
            Phip595CM9->Fill(Php, TaggerTime);
        }

        else if(-0.8> CosThetapCM && CosThetapCM > -1){
            Phip595CM10->Fill(Php, TaggerTime);
        }
    }

    else if ( 600 < EGamma && EGamma < 610) {

        if(1 > CosThetapCM && CosThetapCM > 0.8 ){
            Phip605CM1->Fill(Php, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.6){
            Phip605CM2->Fill(Php, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.4){
            Phip605CM3->Fill(Php, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.2){
            Phip605CM4->Fill(Php, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0){
            Phip605CM5->Fill(Php, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.2){
            Phip605CM6->Fill(Php, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.4){
            Phip605CM7->Fill(Php, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.6){
            Phip605CM8->Fill(Php, TaggerTime);
        }

        else if(-0.6> CosThetapCM && CosThetapCM > -0.8){
            Phip605CM9->Fill(Php, TaggerTime);
        }

        else if(-0.8> CosThetapCM && CosThetapCM > -1){
            Phip605CM10->Fill(Php, TaggerTime);
        }
    }

    else if ( 610 < EGamma && EGamma < 620) {

        if(1 > CosThetapCM && CosThetapCM > 0.8 ){
            Phip615CM1->Fill(Php, TaggerTime);
        }

        else if(0.8 > CosThetapCM && CosThetapCM > 0.6){
            Phip615CM2->Fill(Php, TaggerTime);
        }

        else if(0.6 > CosThetapCM && CosThetapCM > 0.4){
            Phip615CM3->Fill(Php, TaggerTime);
        }

        else if(0.4 > CosThetapCM && CosThetapCM > 0.2){
            Phip615CM4->Fill(Php, TaggerTime);
        }

        else if(0.2 > CosThetapCM && CosThetapCM > 0){
            Phip615CM5->Fill(Php, TaggerTime);
        }

        else if(0 > CosThetapCM && CosThetapCM > -0.2){
            Phip615CM6->Fill(Php, TaggerTime);
        }

        else if(-0.2 > CosThetapCM && CosThetapCM > -0.4){
            Phip615CM7->Fill(Php, TaggerTime);
        }

        else if(-0.4 > CosThetapCM && CosThetapCM > -0.6){
            Phip615CM8->Fill(Php, TaggerTime);
        }

        else if(-0.6> CosThetapCM && CosThetapCM > -0.8){
            Phip615CM9->Fill(Php, TaggerTime);
        }

        else if(-0.8> CosThetapCM && CosThetapCM > -1){
            Phip615CM10->Fill(Php, TaggerTime);
        }
    }
}

Bool_t	PNeutPol_Polarimeter_Lin_NoScatt::Write(){
  // Write all GH1's easily

  GTreeManager::Write();
}
