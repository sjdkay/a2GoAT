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

void	PNeutPol_Polarimeter_Lin::ProcessEvent()
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
        Xn = GetTracks()->GetPseudoVertexX(1);
        Yn = GetTracks()->GetPseudoVertexY(1);
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
        WC2nX = GetTracks()->GetMWPC1PosX(1);
        WC2nY = GetTracks()->GetMWPC1PosY(1);
        WC2nZ = GetTracks()->GetMWPC1PosZ(1);
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
        Xn = GetTracks()->GetPseudoVertexX(0);
        Yn = GetTracks()->GetPseudoVertexY(0);
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
        WC2nX = GetTracks()->GetMWPC1PosX(0);
        WC2nY = GetTracks()->GetMWPC1PosY(0);
        WC2nZ = GetTracks()->GetMWPC1PosZ(0);
    }

    else
    {
        return;
    }


    if( Zp > 60 || Zp < -60) return; // Particles selected out from other parts tend to be inside anyway, skip this?

    EpCorr = EpPolCorrect(Ep, Thp); //correct Ep for energy loss in polarimeter

    if(Cut_proton -> IsInside(EpCorr, dEp) == kFALSE) return; // If E loss correct proton is NOT inside p banana drop out

    EpDiff = abs(EpCorr - Ep);

    Pp = sqrt (TMath::Power((Mp + EpCorr),2) - TMath::Power(Mp,2));
    Pn = sqrt (TMath::Power((En + Mn ),2) - TMath::Power(Mn,2));
    GVpCorr = TLorentzVector(Pp*sin(Thp)*cos(Php), Pp*sin(Thp)*sin(Php), Pp*cos(Thp), EpCorr+Mp);

    GVnCorr =  LNeutron4VectorCorr(Zp, GVn, En, Pn , Mn, Phn);
    ThetanCorr = (GVnCorr.Theta())*TMath::RadToDeg();

    WC13Vectn.SetXYZ(WC1nX, WC1nY, WC1nZ);
    WC23Vectn.SetXYZ(WC2nX, WC2nY, WC2nZ);

    GVpCorr3 = GVpCorr.Vect();
    GVnCorr3 = GVnCorr.Vect();
    GVn3Unit = (GVn.Vect()).Unit();
    GVnCorr3Unit = GVnCorr3.Unit();
    pVertex = TVector3(Xp, Yp, Zp);
    nVertex = TVector3(Xn, Yn, Zn);

    kinfit.LinkVariable("beamF",    beamF.Link(),       beamF.LinkSigma());
    kinfit.LinkVariable("protonF",    protonF.Link(),       protonF.LinkSigma());
    kinfit.LinkVariable("neutronF",    neutronF.Link(),       neutronF.LinkSigma());

    vector<string> all_names = {"beamF", "protonF", "neutronF"};
    kinfit.AddConstraint("EnergyMomentumBalance", all_names, EnergyMomentumBalance);

    for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
    {
        TaggerTime = GetTagger()->GetTaggedTime(j); // Get tagged time for event
        EGamma = (GetTagger()->GetTaggedEnergy(j)); // Get Photon energy for event
        //BeamHelicity = GetTrigger()->GetHelicity();
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

        if(Cut_protonKinGood -> IsInside(KinEp, dEp) == kFALSE) continue; // If KinE proton is NOT inside p banana drop out

        RecProtonEpCorr = LProton4VectorKin(EpCorr, ThpRad, PhpRad);
        RecNeutronEpCorr = LNeutron4VectorKin(RecProtonEpCorr);
        MMpEpCorr = RecNeutronEpCorr.M();

        if (((MMpEpCorr < 800) == kTRUE) || ((MMpEpCorr > 1100) == kTRUE)) continue; // Force a missing mass cut

        PhiDiff = abs (Php-Phn); // This will always be 180?
        // if ( PhiDiff > 195 || PhiDiff < 165) return;

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

        RecProtonEpCorr3 = RecProtonEpCorr.Vect();
        RecNeutronEpCorr3 = RecNeutronEpCorr.Vect();

        P3Vect = RecKinProton.Vect();
        N3Vect = RecKinNeutron.Vect();
        N3VectUnit = N3Vect.Unit();
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

        DOCAVertex1 = TVector3(0., 0., 0.);
        DOCAVertex2 = TVector3(0., 0., 0.);
        DOCA = Calc_dtfInterDOCA(N3VectUnit, GVn3Unit, pVertex, nVertex, DOCAVertex1, DOCAVertex2);
        POCAx = DOCAVertex1.X()-(DOCAVertex1.X()-DOCAVertex2.X())/2.0;
        POCAy = DOCAVertex1.Y()-(DOCAVertex1.Y()-DOCAVertex2.Y())/2.0;
        POCAz = DOCAVertex1.Z()-(DOCAVertex1.Z()-DOCAVertex2.Z())/2.0;
        r = sqrt((TMath::Power(POCAx,2))+(TMath::Power(POCAy,2)));
        POCA = TVector3(POCAx, POCAy, POCAz);
        //if (ScattTheta > 60) continue;

        if( r > 75 || r < 35 ) return; // Ensure POCA is at polarimeter radius

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

TLorentzVector PNeutPol_Polarimeter_Lin::LNeutron4VectorCorr(Double_t ZVert, TLorentzVector n4Vector, Double_t nE, Double_t MagP, Double_t nMass, Double_t nPhi)
{
    Ncor1 = 0.0886404-0.000555077*ZVert+0.000914921*ZVert*ZVert-7.6616e-06*ZVert*ZVert*ZVert;
    Ncor2=0.991;
    Ncor3= 0.612847+0.153167*ZVert-0.00106208*ZVert*ZVert;
    NcorR=Ncor1+Ncor2*(n4Vector.Theta()*180/acos(-1))+Ncor3*sin(n4Vector.Theta());
    NcorRR=NcorR/180.0*acos(-1);

    N4VectCorr =  TLorentzVector (MagP*sin(NcorRR)*cos(nPhi),MagP*sin(NcorRR)*sin(nPhi) , MagP*cos(NcorRR) , nE+nMass);

    return N4VectCorr;
}

TLorentzVector PNeutPol_Polarimeter_Lin::LProton4VectorKin(Double_t KinE, Double_t Theta, Double_t Phi)
{
    EpTot = KinE + Mp;
    PpKin = sqrt (TMath::Power(EpTot,2) - TMath::Power(Mp,2));
    Ppx = PpKin*sin(Theta)*cos(Phi);
    Ppy = PpKin*sin(Theta)*sin(Phi);
    Ppz = PpKin*cos(Theta);
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
    PhiDifference = new GH1 ("PhiDifference", "#phi_{Diff} Between p and n", 180, 0, 360);
    EpKin = new GH1 ("EpKin", "Ep Calculated from Ep/Thetap", 100, 0, 500);
    EpCorrected = new GH1 ("EpCorrected", "Ep Corrected for Energy Loss in Polarimeter ", 100, 0, 500);
    OAngle = new GH1 ("OAngle", "Opening Angle between P and N Vectors", 180, 0, 180);
    WCZnRecon = new GH1 ("WCZnRecon", "WCZ Hit Position from Reconstructed n Vector", 200, -400, 400);

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

    // Proton Phi dists across EGamma bins
    Phip435CM1 = new GH1("Phip_435MeVCM1", "#phi_{p} Distribution for E_{#gamma} 435 #pm 10MeV (Cos#theta_{CM}1-0.8)", 10, -180, 180);
    Phip455CM1 = new GH1("Phip_455MeVCM1", "#phi_{p} Distribution for E_{#gamma} 455 #pm 10MeV (Cos#theta_{CM}1-0.8)", 10, -180, 180);
    Phip475CM1 = new GH1("Phip_475MeVCM1", "#phi_{p} Distribution for E_{#gamma} 475 #pm 10MeV (Cos#theta_{CM}1-0.8)", 10, -180, 180);
    Phip495CM1 = new GH1("Phip_495MeVCM1", "#phi_{p} Distribution for E_{#gamma} 495 #pm 10MeV (Cos#theta_{CM}1-0.8)", 10, -180, 180);
    Phip515CM1 = new GH1("Phip_515MeVCM1", "#phi_{p} Distribution for E_{#gamma} 515 #pm 10MeV (Cos#theta_{CM}1-0.8)", 10, -180, 180);
    Phip535CM1 = new GH1("Phip_535MeVCM1", "#phi_{p} Distribution for E_{#gamma} 535 #pm 10MeV (Cos#theta_{CM}1-0.8)", 10, -180, 180);
    Phip555CM1 = new GH1("Phip_555MeVCM1", "#phi_{p} Distribution for E_{#gamma} 555 #pm 10MeV (Cos#theta_{CM}1-0.8)", 10, -180, 180);
    Phip575CM1 = new GH1("Phip_575MeVCM1", "#phi_{p} Distribution for E_{#gamma} 575 #pm 10MeV (Cos#theta_{CM}1-0.8)", 10, -180, 180);
    Phip595CM1 = new GH1("Phip_595MeVCM1", "#phi_{p} Distribution for E_{#gamma} 595 #pm 10MeV (Cos#theta_{CM}1-0.8)", 10, -180, 180);
    Phip615CM1 = new GH1("Phip_615MeVCM1", "#phi_{p} Distribution for E_{#gamma} 615 #pm 10MeV (Cos#theta_{CM}1-0.8)", 10, -180, 180);

    // #phi_{p} dists across EGamma bins
    Phip435CM2 = new GH1("Phip_435MeVCM2", "#phi_{p} Distribution for E_{#gamma} 435 #pm 10MeV (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
    Phip455CM2 = new GH1("Phip_455MeVCM2", "#phi_{p} Distribution for E_{#gamma} 455 #pm 10MeV (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
    Phip475CM2 = new GH1("Phip_475MeVCM2", "#phi_{p} Distribution for E_{#gamma} 475 #pm 10MeV (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
    Phip495CM2 = new GH1("Phip_495MeVCM2", "#phi_{p} Distribution for E_{#gamma} 495 #pm 10MeV (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
    Phip515CM2 = new GH1("Phip_515MeVCM2", "#phi_{p} Distribution for E_{#gamma} 515 #pm 10MeV (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
    Phip535CM2 = new GH1("Phip_535MeVCM2", "#phi_{p} Distribution for E_{#gamma} 535 #pm 10MeV (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
    Phip555CM2 = new GH1("Phip_555MeVCM2", "#phi_{p} Distribution for E_{#gamma} 555 #pm 10MeV (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
    Phip575CM2 = new GH1("Phip_575MeVCM2", "#phi_{p} Distribution for E_{#gamma} 575 #pm 10MeV (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
    Phip595CM2 = new GH1("Phip_595MeVCM2", "#phi_{p} Distribution for E_{#gamma} 595 #pm 10MeV (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);
    Phip615CM2 = new GH1("Phip_615MeVCM2", "#phi_{p} Distribution for E_{#gamma} 615 #pm 10MeV (Cos#theta_{CM}0.8-0.6)", 10, -180, 180);

    // #phi_{p} dists across EGamma bins
    Phip435CM3 = new GH1("Phip_435MeVCM3", "#phi_{p} Distribution for E_{#gamma} 435 #pm 10MeV (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
    Phip455CM3 = new GH1("Phip_455MeVCM3", "#phi_{p} Distribution for E_{#gamma} 455 #pm 10MeV (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
    Phip475CM3 = new GH1("Phip_475MeVCM3", "#phi_{p} Distribution for E_{#gamma} 475 #pm 10MeV (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
    Phip495CM3 = new GH1("Phip_495MeVCM3", "#phi_{p} Distribution for E_{#gamma} 495 #pm 10MeV (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
    Phip515CM3 = new GH1("Phip_515MeVCM3", "#phi_{p} Distribution for E_{#gamma} 515 #pm 10MeV (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
    Phip535CM3 = new GH1("Phip_535MeVCM3", "#phi_{p} Distribution for E_{#gamma} 535 #pm 10MeV (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
    Phip555CM3 = new GH1("Phip_555MeVCM3", "#phi_{p} Distribution for E_{#gamma} 555 #pm 10MeV (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
    Phip575CM3 = new GH1("Phip_575MeVCM3", "#phi_{p} Distribution for E_{#gamma} 575 #pm 10MeV (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
    Phip595CM3 = new GH1("Phip_595MeVCM3", "#phi_{p} Distribution for E_{#gamma} 595 #pm 10MeV (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);
    Phip615CM3 = new GH1("Phip_615MeVCM3", "#phi_{p} Distribution for E_{#gamma} 615 #pm 10MeV (Cos#theta_{CM}0.6-0.4)", 10, -180, 180);

    // #phi_{p} dists across EGamma bins
    Phip435CM4 = new GH1("Phip_435MeVCM4", "#phi_{p} Distribution for E_{#gamma} 435 #pm 10MeV (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
    Phip455CM4 = new GH1("Phip_455MeVCM4", "#phi_{p} Distribution for E_{#gamma} 455 #pm 10MeV (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
    Phip475CM4 = new GH1("Phip_475MeVCM4", "#phi_{p} Distribution for E_{#gamma} 475 #pm 10MeV (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
    Phip495CM4 = new GH1("Phip_495MeVCM4", "#phi_{p} Distribution for E_{#gamma} 495 #pm 10MeV (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
    Phip515CM4 = new GH1("Phip_515MeVCM4", "#phi_{p} Distribution for E_{#gamma} 515 #pm 10MeV (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
    Phip535CM4 = new GH1("Phip_535MeVCM4", "#phi_{p} Distribution for E_{#gamma} 535 #pm 10MeV (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
    Phip555CM4 = new GH1("Phip_555MeVCM4", "#phi_{p} Distribution for E_{#gamma} 555 #pm 10MeV (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
    Phip575CM4 = new GH1("Phip_575MeVCM4", "#phi_{p} Distribution for E_{#gamma} 575 #pm 10MeV (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
    Phip595CM4 = new GH1("Phip_595MeVCM4", "#phi_{p} Distribution for E_{#gamma} 595 #pm 10MeV (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);
    Phip615CM4 = new GH1("Phip_615MeVCM4", "#phi_{p} Distribution for E_{#gamma} 615 #pm 10MeV (Cos#theta_{CM}0.4-0.2)", 10, -180, 180);

    // #phi_{p} dists across EGamma bins
    Phip435CM5 = new GH1("Phip_435MeVCM5", "#phi_{p} Distribution for E_{#gamma} 435 #pm 10MeV (Cos#theta_{CM}0.2-0)", 10, -180, 180);
    Phip455CM5 = new GH1("Phip_455MeVCM5", "#phi_{p} Distribution for E_{#gamma} 455 #pm 10MeV (Cos#theta_{CM}0.2-0)", 10, -180, 180);
    Phip475CM5 = new GH1("Phip_475MeVCM5", "#phi_{p} Distribution for E_{#gamma} 475 #pm 10MeV (Cos#theta_{CM}0.2-0)", 10, -180, 180);
    Phip495CM5 = new GH1("Phip_495MeVCM5", "#phi_{p} Distribution for E_{#gamma} 495 #pm 10MeV (Cos#theta_{CM}0.2-0)", 10, -180, 180);
    Phip515CM5 = new GH1("Phip_515MeVCM5", "#phi_{p} Distribution for E_{#gamma} 515 #pm 10MeV (Cos#theta_{CM}0.2-0)", 10, -180, 180);
    Phip535CM5 = new GH1("Phip_535MeVCM5", "#phi_{p} Distribution for E_{#gamma} 535 #pm 10MeV (Cos#theta_{CM}0.2-0)", 10, -180, 180);
    Phip555CM5 = new GH1("Phip_555MeVCM5", "#phi_{p} Distribution for E_{#gamma} 555 #pm 10MeV (Cos#theta_{CM}0.2-0)", 10, -180, 180);
    Phip575CM5 = new GH1("Phip_575MeVCM5", "#phi_{p} Distribution for E_{#gamma} 575 #pm 10MeV (Cos#theta_{CM}0.2-0)", 10, -180, 180);
    Phip595CM5 = new GH1("Phip_595MeVCM5", "#phi_{p} Distribution for E_{#gamma} 595 #pm 10MeV (Cos#theta_{CM}0.2-0)", 10, -180, 180);
    Phip615CM5 = new GH1("Phip_615MeVCM5", "#phi_{p} Distribution for E_{#gamma} 615 #pm 10MeV (Cos#theta_{CM}0.2-0)", 10, -180, 180);

    // #phi_{p} dists across EGamma bins
    Phip435CM6 = new GH1("Phip_435MeVCM6", "#phi_{p} Distribution for E_{#gamma} 435 #pm 10MeV (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
    Phip455CM6 = new GH1("Phip_455MeVCM6", "#phi_{p} Distribution for E_{#gamma} 455 #pm 10MeV (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
    Phip475CM6 = new GH1("Phip_475MeVCM6", "#phi_{p} Distribution for E_{#gamma} 475 #pm 10MeV (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
    Phip495CM6 = new GH1("Phip_495MeVCM6", "#phi_{p} Distribution for E_{#gamma} 495 #pm 10MeV (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
    Phip515CM6 = new GH1("Phip_515MeVCM6", "#phi_{p} Distribution for E_{#gamma} 515 #pm 10MeV (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
    Phip535CM6 = new GH1("Phip_535MeVCM6", "#phi_{p} Distribution for E_{#gamma} 535 #pm 10MeV (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
    Phip555CM6 = new GH1("Phip_555MeVCM6", "#phi_{p} Distribution for E_{#gamma} 555 #pm 10MeV (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
    Phip575CM6 = new GH1("Phip_575MeVCM6", "#phi_{p} Distribution for E_{#gamma} 575 #pm 10MeV (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
    Phip595CM6 = new GH1("Phip_595MeVCM6", "#phi_{p} Distribution for E_{#gamma} 595 #pm 10MeV (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);
    Phip615CM6 = new GH1("Phip_615MeVCM6", "#phi_{p} Distribution for E_{#gamma} 615 #pm 10MeV (Cos#theta_{CM}0-(-0.2))", 10, -180, 180);

    // #phi_{p} dists across EGamma bins
    Phip435CM7 = new GH1("Phip_435MeVCM7", "#phi_{p} Distribution for E_{#gamma} 435 #pm 10MeV (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
    Phip455CM7 = new GH1("Phip_455MeVCM7", "#phi_{p} Distribution for E_{#gamma} 455 #pm 10MeV (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
    Phip475CM7 = new GH1("Phip_475MeVCM7", "#phi_{p} Distribution for E_{#gamma} 475 #pm 10MeV (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
    Phip495CM7 = new GH1("Phip_495MeVCM7", "#phi_{p} Distribution for E_{#gamma} 495 #pm 10MeV (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
    Phip515CM7 = new GH1("Phip_515MeVCM7", "#phi_{p} Distribution for E_{#gamma} 515 #pm 10MeV (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
    Phip535CM7 = new GH1("Phip_535MeVCM7", "#phi_{p} Distribution for E_{#gamma} 535 #pm 10MeV (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
    Phip555CM7 = new GH1("Phip_555MeVCM7", "#phi_{p} Distribution for E_{#gamma} 555 #pm 10MeV (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
    Phip575CM7 = new GH1("Phip_575MeVCM7", "#phi_{p} Distribution for E_{#gamma} 575 #pm 10MeV (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
    Phip595CM7 = new GH1("Phip_595MeVCM7", "#phi_{p} Distribution for E_{#gamma} 595 #pm 10MeV (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);
    Phip615CM7 = new GH1("Phip_615MeVCM7", "#phi_{p} Distribution for E_{#gamma} 615 #pm 10MeV (Cos#theta_{CM}-0.2-(-0.4))", 10, -180, 180);

    // #phi_{p} dists across EGamma bins
    Phip435CM8 = new GH1("Phip_435MeVCM8", "#phi_{p} Distribution for E_{#gamma} 435 #pm 10MeV (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
    Phip455CM8 = new GH1("Phip_455MeVCM8", "#phi_{p} Distribution for E_{#gamma} 455 #pm 10MeV (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
    Phip475CM8 = new GH1("Phip_475MeVCM8", "#phi_{p} Distribution for E_{#gamma} 475 #pm 10MeV (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
    Phip495CM8 = new GH1("Phip_495MeVCM8", "#phi_{p} Distribution for E_{#gamma} 495 #pm 10MeV (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
    Phip515CM8 = new GH1("Phip_515MeVCM8", "#phi_{p} Distribution for E_{#gamma} 515 #pm 10MeV (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
    Phip535CM8 = new GH1("Phip_535MeVCM8", "#phi_{p} Distribution for E_{#gamma} 535 #pm 10MeV (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
    Phip555CM8 = new GH1("Phip_555MeVCM8", "#phi_{p} Distribution for E_{#gamma} 555 #pm 10MeV (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
    Phip575CM8 = new GH1("Phip_575MeVCM8", "#phi_{p} Distribution for E_{#gamma} 575 #pm 10MeV (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
    Phip595CM8 = new GH1("Phip_595MeVCM8", "#phi_{p} Distribution for E_{#gamma} 595 #pm 10MeV (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);
    Phip615CM8 = new GH1("Phip_615MeVCM8", "#phi_{p} Distribution for E_{#gamma} 615 #pm 10MeV (Cos#theta_{CM}-0.4-(-0.6))", 10, -180, 180);

    // #phi_{p} dists across EGamma bins
    Phip435CM9 = new GH1("Phip_435MeVCM9", "#phi_{p} Distribution for E_{#gamma} 435 #pm 10MeV (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
    Phip455CM9 = new GH1("Phip_455MeVCM9", "#phi_{p} Distribution for E_{#gamma} 455 #pm 10MeV (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
    Phip475CM9 = new GH1("Phip_475MeVCM9", "#phi_{p} Distribution for E_{#gamma} 475 #pm 10MeV (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
    Phip495CM9 = new GH1("Phip_495MeVCM9", "#phi_{p} Distribution for E_{#gamma} 495 #pm 10MeV (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
    Phip515CM9 = new GH1("Phip_515MeVCM9", "#phi_{p} Distribution for E_{#gamma} 515 #pm 10MeV (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
    Phip535CM9 = new GH1("Phip_535MeVCM9", "#phi_{p} Distribution for E_{#gamma} 535 #pm 10MeV (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
    Phip555CM9 = new GH1("Phip_555MeVCM9", "#phi_{p} Distribution for E_{#gamma} 555 #pm 10MeV (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
    Phip575CM9 = new GH1("Phip_575MeVCM9", "#phi_{p} Distribution for E_{#gamma} 575 #pm 10MeV (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
    Phip595CM9 = new GH1("Phip_595MeVCM9", "#phi_{p} Distribution for E_{#gamma} 595 #pm 10MeV (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);
    Phip615CM9 = new GH1("Phip_615MeVCM9", "#phi_{p} Distribution for E_{#gamma} 615 #pm 10MeV (Cos#theta_{CM}-0.6-(-0.8))", 10, -180, 180);

    // #phi_{p} dists across EGamma bins
    Phip435CM10 = new GH1("Phip_435MeVCM10", "#phi_{p} Distribution for E_{#gamma} 435 #pm 10MeV (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
    Phip455CM10 = new GH1("Phip_455MeVCM10", "#phi_{p} Distribution for E_{#gamma} 455 #pm 10MeV (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
    Phip475CM10 = new GH1("Phip_475MeVCM10", "#phi_{p} Distribution for E_{#gamma} 475 #pm 10MeV (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
    Phip495CM10 = new GH1("Phip_495MeVCM10", "#phi_{p} Distribution for E_{#gamma} 495 #pm 10MeV (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
    Phip515CM10 = new GH1("Phip_515MeVCM10", "#phi_{p} Distribution for E_{#gamma} 515 #pm 10MeV (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
    Phip535CM10 = new GH1("Phip_535MeVCM10", "#phi_{p} Distribution for E_{#gamma} 535 #pm 10MeV (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
    Phip555CM10 = new GH1("Phip_555MeVCM10", "#phi_{p} Distribution for E_{#gamma} 555 #pm 10MeV (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
    Phip575CM10 = new GH1("Phip_575MeVCM10", "#phi_{p} Distribution for E_{#gamma} 575 #pm 10MeV (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
    Phip595CM10 = new GH1("Phip_595MeVCM10", "#phi_{p} Distribution for E_{#gamma} 595 #pm 10MeV (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);
    Phip615CM10 = new GH1("Phip_615MeVCM10", "#phi_{p} Distribution for E_{#gamma} 615 #pm 10MeV (Cos#theta_{CM}-0.8-(-1))", 10, -180, 180);

    ThetanDist = new GH1 ("ThetanDist", "#theta_{n} Distribution", 200, 0, 180);
    ThetanRecDist = new GH1 ("ThetanRecDist", "Reconstructed #theta_{n} Distribution", 200, 0, 180);
    ThetanDiffDist = new GH1 ("ThetanDiffDist", "Difference Between #theta_{n} and  #theta_{nRec}", 200, -90, 90);
    ThetanDiffZp = new GH2 ("ThetanDiffZp", "Diff(#theta_{n} - #theta_{nRec}) as a Fn of Z_{p}", 200, -90, 90, 200, -100, 100);

    ThetanCorrDist = new GH1 ("ThetanCorrDist", "#theta_{nCorr} Distribution", 200, 0, 180);
    ThetanCorrDiffDist = new GH1 ("ThetanCorrDiffDist", "Difference Between #theta_{n} and #theta_{nCorr} Distribution", 200, -90, 90);
    ThetanCorrRecDiffDist = new GH1 ("ThetanCorrRecDiffDist", "Difference Between #theta_{nCorr} and  #theta_{nRec}", 200, -90, 90);
    ThetanCorrDiffZp = new GH2 ("ThetanCorrDiffZp", "Diff(#theta_{nCorr} - #theta_{nRec}) as a Fn of Z_{p}", 200, -90, 90, 200, -100, 100);

    ThetaRecPiDiff = new GH1 ("ThetaRecPiDiff", "Difference between #theta_{#pi Rec} and #theta_{n}", 200, -90, 90);
    ThetanThetaRecPi = new GH2 ("ThetanThetaRecPi", "#theta_{n} vs #theta_{#pi rec}", 100, 0, 180, 100, 0, 180);
    ThetanThetaRecPiDiff = new GH2 ("ThetanThetaRecPiDiff", "#theta_{n} vs (#theta_{#pi Rec} - #theta_{n})", 100, 0, 180, 100, -90, 90);

    ThetaRecPDiff = new GH1 ("ThetaRecPDiff", "Difference between #theta_{pRec} and #theta_{n}", 200, -90, 90);
    ThetanThetaRecP = new GH2 ("ThetanThetaRecP", "#theta_{n} vs #theta_{pRec}", 100, 0, 180, 100, 0, 180);
    ThetanThetaRecPDiff = new GH2 ("ThetanThetaRecPDiff", "#theta_{n} vs (#theta_{pRec} - #theta_{n})", 100, 0, 180, 100, -90, 90);

    DeutKinPiKin = new GH2 ("DeutKinPiKin", "(#theta_{nRec} - #theta_{n}) vs (#theta_{#pi Rec} - #theta_{n})", 200, -180, 180, 200, -180, 180);

    E_dE = new GH2 ("E_dE", "EdE Plot With E Loss Adjustment", 100, 0, 500, 100, 0, 5);
    KinEp_dE = new GH2 ("KinEp_dE", "KinEpdE Plot", 100, 0, 500, 100, 0, 5);
    //KinEp_dE_GoodCut = new GH2 ("KinEp_dE_GoodCut", "KinEpdE Plot With Good Proton Cut", 100, 0, 500, 100, 0, 5);
    ThetaScPhiSc = new GH2 ("ThetaScPhiSc", "Phi as a function of Theta (Both in rotated frame)", 100, 0, 180, 100, -180, 180);
    E_KinEp = new GH2 ("E_KinEp", "Kinematic Energy of Proton as a function of CB energy", 100, 0, 500, 100, 0, 500);
    PhinDiffWCZRec = new GH2 ("PhinDiffWCZRec", "Difference between WC Phi and Reconstructed Phi as a fn of WCZ Hit Position", 100, 0, 200, 100, 0, 180);
    PhinDiffWCZRec_KinCut = new GH2 ("PhinDiffWCZRec_KinCut", "Difference between WC Phi and Reconstructed Phi as a fn of WCZ Hit Position (Kin P Banana Cut)", 200, -300, 300, 200, 0, 180);

    ClosestApproach = new GH1("ClosestApproach", "DOCA of n and p' vectors", 200, -200, 200);
    POCAr = new GH1("POCAr", "r_{POCA}", 200, 0, 300);
    ScatterVertexZ = new GH1("ScatterVertexZ", "Z_{POCA}", 200, -200, 200);
    ScatterVertexZr = new GH2("ScatterVertexZr", "Z_{POCA} vs r_{POCA}", 200, -200, 200, 200, 0, 200);
    ScatterVertexXY = new GH2("ScatterVertexXY", "XY Vertex Point of Scatter from DOCA Method", 100, -80, 80, 100, -80, 80);
    ScatterVertex = new GH3("ScatterVertex", "Vertex Point of Scatter from DOCA Method", 100, -80, 80, 100, -80, 80, 100, -200, 200);

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

    PhiDifference->Fill(PhiDiff);
    E_KinEp->Fill(EpCorr, KinEp, TaggerTime);

    ThetanDist->Fill(Thn, TaggerTime);
    ThetanRecDist->Fill(ThetanRec, TaggerTime);
    ThetanDiffDist->Fill(Thetan-ThetanRec, TaggerTime);
    ThetanDiffZp->Fill(Thn-ThetanRec, Zp, TaggerTime);

    ThetanCorrDist->Fill(ThetanCorr, TaggerTime);
    ThetanCorrDiffDist->Fill(Thn-ThetanCorr, TaggerTime);
    ThetanCorrRecDiffDist->Fill(ThetanCorr-ThetanRec, TaggerTime);
    ThetanCorrDiffZp ->Fill(ThetanCorr-ThetanRec, Zp, TaggerTime);

    ThetaRecPiDiff->Fill(ThetaPiRecDiff, TaggerTime);
    ThetanThetaRecPi->Fill(Thn, ThetaPiRec, TaggerTime);
    ThetanThetaRecPiDiff->Fill(Thn, ThetaPiRecDiff, TaggerTime);

    ThetaRecPDiff->Fill(ThetapRecDiff, TaggerTime);
    ThetanThetaRecP->Fill(Thn, ThetapRec, TaggerTime);
    ThetanThetaRecPDiff->Fill(Thn, ThetapRecDiff, TaggerTime);

    DeutKinPiKin->Fill(ThetanRec-Thn, ThetaPiRecDiff, TaggerTime);

    ClosestApproach->Fill(DOCA, TaggerTime);
    POCAr->Fill(r, TaggerTime);
    ScatterVertexZ->Fill(POCAz, TaggerTime);
    ScatterVertexZr->Fill(POCAz, r, TaggerTime);
    ScatterVertexXY->Fill(POCAx, POCAy, TaggerTime);
    ScatterVertex->Fill(POCAx, POCAy, POCAz, TaggerTime);

    PhinDiffWCZRec->Fill(WCZnRec, PhinDiff, TaggerTime);

    ThetaRecPDiff->Fill(ThetapRecDiff, TaggerTime);
    ThetanThetaRecP->Fill(Thn, ThetapRec, TaggerTime);
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

    if ( 425 < EGamma && EGamma < 445) {

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


    else if ( 445 < EGamma && EGamma < 465) {

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

    else if ( 465 < EGamma && EGamma < 485) {

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

    else if ( 485 < EGamma && EGamma < 505) {

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

    else if ( 505 < EGamma && EGamma < 525) {

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

    else if ( 525 < EGamma && EGamma < 545) {

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

    else if ( 545 < EGamma && EGamma < 565) {

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


    else if ( 565 < EGamma && EGamma < 585) {

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

    else if ( 585 < EGamma && EGamma < 605) {

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

    else if ( 605 < EGamma && EGamma < 625) {

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

Bool_t	PNeutPol_Polarimeter_Lin::Write(){
  // Write all GH1's easily

  GTreeManager::Write();
}
