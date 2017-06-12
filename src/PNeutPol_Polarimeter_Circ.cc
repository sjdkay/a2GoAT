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

    NP = 0; // Set number of Protons to 0 before checking
    NPi = 0; // Set number of pions to 0 before checking
    NRoo = 0; // Set number of Rootinos to 0 before checking
    Mn = 939.565; // Mass of neutron in MeV
    Mp = 938.272; // Mass of proton in MeV
    Md = 1875.613; //Mass of Deuterium in MeV
    Mpi = 139.57018; // Mass of charged pion in MeV
    Deut = TLorentzVector (0., 0., 0., 1875.613); // 4-Vector of Deuterium target, assume at rest
    Neut = TLorentzVector (0., 0., 0., 939.565); // 4-Vector of Neutron target, assume at rest

    Cut_CB_proton = OpenCutFile("configfiles/cuts/CB_DeltaE-E_Proton_11_05_17.root", "Proton"); // These will need adjusting with new Acqu files
    Cut_proton = Cut_CB_proton;
    Cut_CB_pion = OpenCutFile("configfiles/cuts/CB_DeltaE-E_Pion_29_07_15.root", "Pion");
    Cut_pion = Cut_CB_pion;
    Cut_CB_protonKinGood = OpenCutFile("configfiles/cuts/CB_DeltaE-E_ProtonKinGood_11_05_17.root", "ProtonKinGood"); // These will need adjusting with new Acqu files
    Cut_protonKinGood = Cut_CB_protonKinGood;
    Cut_CB_protonKinBad = OpenCutFile("configfiles/cuts/CB_DeltaE-E_ProtonKinBad_15_12_16.root", "ProtonKinBad");
    Cut_protonKinBad = Cut_CB_protonKinBad;
    cout << endl;

    TraverseValidEvents(); // This loops over each event as in old file and calls ProcessEvent() each loop

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

    GVnCorr =  CNeutron4VectorCorr(Zp, GVn, En, Pn , Mn, Phn);
    ThetanCorr = (GVnCorr.Theta())*TMath::RadToDeg();

    WC3Vectp.SetXYZ(WC1pX, WC1pY, WC1pZ);
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
        RecKinProton = CProton4VectorKin(KinEp, ThpRad, PhpRad);
        RecKinNeutron = CNeutron4VectorKin(RecKinProton);
        ThetanRec = (RecKinNeutron.Theta()) * TMath::RadToDeg();
        PhinRec = (RecKinNeutron.Phi()) * TMath::RadToDeg();
        WCZnRec = 72/tan(RecKinNeutron.Theta());

        if(Cut_protonKinGood -> IsInside(KinEp, dEp) == kFALSE) continue; // If KinE proton is NOT inside p banana drop out

        RecProtonEpCorr = CProton4VectorKin(EpCorr, ThpRad, PhpRad);
        RecNeutronEpCorr = CNeutron4VectorKin(RecProtonEpCorr);
        MMpEpCorr = RecNeutronEpCorr.M();

        if (((MMpEpCorr < 800) == kTRUE) || ((MMpEpCorr > 1100) == kTRUE)) continue; // Force a missing mass cut

        PhiDiff = abs (Php-Phn); // This will always be 180?
        // if ( PhiDiff > 195 || PhiDiff < 165) return;

        // Gamma(n,p)Pi (P detected correct)
        // Assume proton track is proton and "neutron" track is from charged pion
        KinEpPi = CalcKinEnergy(Thp, EGamma, Mn, 0, Mp, Mpi); // Calculate kin E of proton assuming g(n, p) pi
        RecKinProtonPi = CProton4VectorKin(KinEpPi, ThpRad, PhpRad); // Get Proton 4 vector from calculated kin E
        RecKinPion = CPion4VectorKin(RecKinProtonPi); // Get Pion 4 vector from 4 momenta conservation
        ThetaPiRec = (RecKinPion.Theta())*TMath::RadToDeg();
        PhiPiRec = (RecKinPion.Phi())*TMath::RadToDeg();
        ThetaPiRecDiff = ThetaPiRec - ThetanCorr;

        // Gamma(n,p)Pi (Pion detected correct)
        // Assume proton track is pion and "neutron" track is from proton
        KinPi = CalcKinEnergy(Thp, EGamma, Mn, 0, Mpi, Mp); // Calculate kin E of pion
        RecKinPionP = CProton4VectorKin(KinPi, ThpRad, PhpRad); // Get Pion 4 vector from calculated kinE
        RecKinPPi = CPion4VectorKin(RecKinPionP); // Get Proton 4 vector from 4 momenta conservation
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
    Ncor2=0.991;
    Ncor3= 0.612847+0.153167*ZVert-0.00106208*ZVert*ZVert;
    NcorR=Ncor1+Ncor2*(n4Vector.Theta()*180/acos(-1))+Ncor3*sin(n4Vector.Theta());
    NcorRR=NcorR/180.0*acos(-1);

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

TLorentzVector PNeutPol_Polarimeter_Circ::CPion4VectorKin(TLorentzVector ProtonKinVector)
{

    Pi4Vect = (Gamma + Neut) - ProtonKinVector;

    return Pi4Vect;
}

PNeutPol_Polarimeter_Circ::PNeutPol_Polarimeter_Circ() // Define a load of histograms to fill
{
    time = new TH1D("time", 	"time", 	1400, -700, 700);
    time_cut = new TH1D("time_cut", 	"time_cut", 	1400, -700, 700);

    Eg = new GH1( "Eg", "E_{#gamma} Distribution", 200, 100, 1600 );
    PhiDifference = new GH1 ("PhiDifference", "#phi_{Diff} Between p and n", 180, 0, 360);
    EpKin = new GH1 ("EpKin", "E_{p} Calculated from E_{p}/#theta_{p}", 100, 0, 500);
    EpCorrected = new GH1 ("EpCorrected", "E_{pCorr}", 100, 0, 500);
    OAngle = new GH1 ("OAngle", "Opening Angle between P and N Vectors", 180, 0, 180);
    WCZnRecon = new GH1 ("WCZnRecon", "WCZ Hit Position from Reconstructed n Vector", 400, -400, 400);

    ThetaSc =  new GH1( "Theta_Scattered", "Scattered Proton Theta Distribution in Rotated Frame", 180, 0, 180 );
    PhiSc = new GH1( "Phi_Scattered", "Scattered Proton Phi Distribution in Rotated Frame", 90, -180, 180 );

    EpKinEpCorrDiff = new GH1("EpKinEpCorrDiff", "Difference Between E_{pKin} and E_{pCorr}", 300, -300, 300);
    EpEpCorrDiff = new GH1("EpEpCorrDiff", "Difference Between E_{p} and E_{pCorr}", 200, 0, 200);
    MMpEpCorrected = new GH1 ("MMpEpCorrected", "Missing mass seen by Proton (E Loss Corrected)", 400, 0, 2000);

    ZpDist = new GH1 ("ZpDist", "Proton Pseudo Z Vertex Distribution", 200, -400, 400);
    ZpPhiScatNeg180 = new GH1("ZpPhiScatNeg180", "Proton Pseudo Vertex Z for events with #phi_{Sc} ~ -ve180", 200, -200, 200);
    ZpPhiScat0 = new GH1("ZpPhiScat0", "Proton Pseudo Vertex Z for events with #phi_{Sc} ~ 0", 200, -200, 200);
    ZpPhiScatPos180 = new GH1("ZpPhiScatPos180", "Proton Pseudo Vertex Z for events with #phi_{Sc} ~ 180", 200, -200, 200);

    // MMp across photon E bins
    MMp200300 = new GH1("MMp200300", "Missing mass as seen by Proton (200-300MeV E_{#gamma})", 400, 0, 2000);
    MMp300400 = new GH1("MMp300400", "Missing mass as seen by Proton (300-400MeV E_{#gamma})", 400, 0, 2000);
    MMp400500 = new GH1("MMp400500", "Missing mass as seen by Proton (400-500MeV E_{#gamma})", 400, 0, 2000);
    MMp500600 = new GH1("MMp500600", "Missing mass as seen by Proton (500-600MeV E_{#gamma})", 400, 0, 2000);
    MMp600700 = new GH1("MMp600700", "Missing mass as seen by Proton (600-700MeV E_{#gamma})", 400, 0, 2000);
    MMp700800 = new GH1("MMp700800", "Missing mass as seen by Proton (700-800MeV E_{#gamma})", 400, 0, 2000);
    MMp800900 = new GH1("MMp800900", "Missing mass as seen by Proton (800-900MeV E_{#gamma})", 400, 0, 2000);

    // Angles of neutron in scattered frame across EGamma bins for negative helicity
    PhiSc315NegHelCM1 = new GH1( "Phi_Scattered_315MeV_NegHelCM1", "#phi_{Sc} E_{#gamma}315 #pm 15MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 10, -180, 180);
    PhiSc345NegHelCM1 = new GH1( "Phi_Scattered_345MeV_NegHelCM1", "#phi_{Sc} E_{#gamma}345 #pm 15MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 10, -180, 180);
    PhiSc375NegHelCM1 = new GH1( "Phi_Scattered_375MeV_NegHelCM1", "#phi_{Sc} E_{#gamma}375 #pm 15MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 10, -180, 180);
    PhiSc405NegHelCM1 = new GH1( "Phi_Scattered_405MeV_NegHelCM1", "#phi_{Sc} E_{#gamma}405 #pm 15MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 10, -180, 180);
    PhiSc435NegHelCM1 = new GH1( "Phi_Scattered_435MeV_NegHelCM1", "#phi_{Sc} E_{#gamma}435 #pm 15MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 10, -180, 180);
    PhiSc465NegHelCM1 = new GH1( "Phi_Scattered_465MeV_NegHelCM1", "#phi_{Sc} E_{#gamma}465 #pm 15MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 10, -180, 180);
    PhiSc495NegHelCM1 = new GH1( "Phi_Scattered_495MeV_NegHelCM1", "#phi_{Sc} E_{#gamma}495 #pm 15MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 10, -180, 180);
    PhiSc525NegHelCM1 = new GH1( "Phi_Scattered_525MeV_NegHelCM1", "#phi_{Sc} E_{#gamma}525 #pm 15MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 10, -180, 180);
    PhiSc555NegHelCM1 = new GH1( "Phi_Scattered_555MeV_NegHelCM1", "#phi_{Sc} E_{#gamma}555 #pm 15MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 10, -180, 180);
    PhiSc585NegHelCM1 = new GH1( "Phi_Scattered_585MeV_NegHelCM1", "#phi_{Sc} E_{#gamma}585 #pm 15MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 10, -180, 180);
    PhiSc615NegHelCM1 = new GH1( "Phi_Scattered_615MeV_NegHelCM1", "#phi_{Sc} E_{#gamma}615 #pm 15MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 10, -180, 180);
    PhiSc645NegHelCM1 = new GH1( "Phi_Scattered_645MeV_NegHelCM1", "#phi_{Sc} E_{#gamma}645 #pm 15MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 10, -180, 180);
    PhiSc675NegHelCM1 = new GH1( "Phi_Scattered_675MeV_NegHelCM1", "#phi_{Sc} E_{#gamma}675 #pm 15MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 10, -180, 180);

    PhiSc315NegHelCM2 = new GH1( "Phi_Scattered_315MeV_NegHelCM2", "#phi_{Sc} E_{#gamma}315 #pm 15MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 10, -180, 180);
    PhiSc345NegHelCM2 = new GH1( "Phi_Scattered_345MeV_NegHelCM2", "#phi_{Sc} E_{#gamma}345 #pm 15MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 10, -180, 180);
    PhiSc375NegHelCM2 = new GH1( "Phi_Scattered_375MeV_NegHelCM2", "#phi_{Sc} E_{#gamma}375 #pm 15MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 10, -180, 180);
    PhiSc405NegHelCM2 = new GH1( "Phi_Scattered_405MeV_NegHelCM2", "#phi_{Sc} E_{#gamma}405 #pm 15MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 10, -180, 180);
    PhiSc435NegHelCM2 = new GH1( "Phi_Scattered_435MeV_NegHelCM2", "#phi_{Sc} E_{#gamma}435 #pm 15MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 10, -180, 180);
    PhiSc465NegHelCM2 = new GH1( "Phi_Scattered_465MeV_NegHelCM2", "#phi_{Sc} E_{#gamma}465 #pm 15MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 10, -180, 180);
    PhiSc495NegHelCM2 = new GH1( "Phi_Scattered_495MeV_NegHelCM2", "#phi_{Sc} E_{#gamma}495 #pm 15MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 10, -180, 180);
    PhiSc525NegHelCM2 = new GH1( "Phi_Scattered_525MeV_NegHelCM2", "#phi_{Sc} E_{#gamma}525 #pm 15MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 10, -180, 180);
    PhiSc555NegHelCM2 = new GH1( "Phi_Scattered_555MeV_NegHelCM2", "#phi_{Sc} E_{#gamma}555 #pm 15MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 10, -180, 180);
    PhiSc585NegHelCM2 = new GH1( "Phi_Scattered_585MeV_NegHelCM2", "#phi_{Sc} E_{#gamma}585 #pm 15MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 10, -180, 180);
    PhiSc615NegHelCM2 = new GH1( "Phi_Scattered_615MeV_NegHelCM2", "#phi_{Sc} E_{#gamma}615 #pm 15MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 10, -180, 180);
    PhiSc645NegHelCM2 = new GH1( "Phi_Scattered_645MeV_NegHelCM2", "#phi_{Sc} E_{#gamma}645 #pm 15MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 10, -180, 180);
    PhiSc675NegHelCM2 = new GH1( "Phi_Scattered_675MeV_NegHelCM2", "#phi_{Sc} E_{#gamma}675 #pm 15MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 10, -180, 180);

    PhiSc315NegHelCM3 = new GH1( "Phi_Scattered_315MeV_NegHelCM3", "#phi_{Sc} E_{#gamma}315 #pm 15MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 10, -180, 180);
    PhiSc345NegHelCM3 = new GH1( "Phi_Scattered_345MeV_NegHelCM3", "#phi_{Sc} E_{#gamma}345 #pm 15MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 10, -180, 180);
    PhiSc375NegHelCM3 = new GH1( "Phi_Scattered_375MeV_NegHelCM3", "#phi_{Sc} E_{#gamma}375 #pm 15MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 10, -180, 180);
    PhiSc405NegHelCM3 = new GH1( "Phi_Scattered_405MeV_NegHelCM3", "#phi_{Sc} E_{#gamma}405 #pm 15MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 10, -180, 180);
    PhiSc435NegHelCM3 = new GH1( "Phi_Scattered_435MeV_NegHelCM3", "#phi_{Sc} E_{#gamma}435 #pm 15MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 10, -180, 180);
    PhiSc465NegHelCM3 = new GH1( "Phi_Scattered_465MeV_NegHelCM3", "#phi_{Sc} E_{#gamma}465 #pm 15MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 10, -180, 180);
    PhiSc495NegHelCM3 = new GH1( "Phi_Scattered_495MeV_NegHelCM3", "#phi_{Sc} E_{#gamma}495 #pm 15MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 10, -180, 180);
    PhiSc525NegHelCM3 = new GH1( "Phi_Scattered_525MeV_NegHelCM3", "#phi_{Sc} E_{#gamma}525 #pm 15MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 10, -180, 180);
    PhiSc555NegHelCM3 = new GH1( "Phi_Scattered_555MeV_NegHelCM3", "#phi_{Sc} E_{#gamma}555 #pm 15MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 10, -180, 180);
    PhiSc585NegHelCM3 = new GH1( "Phi_Scattered_585MeV_NegHelCM3", "#phi_{Sc} E_{#gamma}585 #pm 15MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 10, -180, 180);
    PhiSc615NegHelCM3 = new GH1( "Phi_Scattered_615MeV_NegHelCM3", "#phi_{Sc} E_{#gamma}615 #pm 15MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 10, -180, 180);
    PhiSc645NegHelCM3 = new GH1( "Phi_Scattered_645MeV_NegHelCM3", "#phi_{Sc} E_{#gamma}645 #pm 15MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 10, -180, 180);
    PhiSc675NegHelCM3 = new GH1( "Phi_Scattered_675MeV_NegHelCM3", "#phi_{Sc} E_{#gamma}675 #pm 15MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 10, -180, 180);

    PhiSc315NegHelCM4 = new GH1( "Phi_Scattered_315MeV_NegHelCM4", "#phi_{Sc} E_{#gamma}315 #pm 15MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 10, -180, 180);
    PhiSc345NegHelCM4 = new GH1( "Phi_Scattered_345MeV_NegHelCM4", "#phi_{Sc} E_{#gamma}345 #pm 15MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 10, -180, 180);
    PhiSc375NegHelCM4 = new GH1( "Phi_Scattered_375MeV_NegHelCM4", "#phi_{Sc} E_{#gamma}375 #pm 15MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 10, -180, 180);
    PhiSc405NegHelCM4 = new GH1( "Phi_Scattered_405MeV_NegHelCM4", "#phi_{Sc} E_{#gamma}405 #pm 15MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 10, -180, 180);
    PhiSc435NegHelCM4 = new GH1( "Phi_Scattered_435MeV_NegHelCM4", "#phi_{Sc} E_{#gamma}435 #pm 15MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 10, -180, 180);
    PhiSc465NegHelCM4 = new GH1( "Phi_Scattered_465MeV_NegHelCM4", "#phi_{Sc} E_{#gamma}465 #pm 15MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 10, -180, 180);
    PhiSc495NegHelCM4 = new GH1( "Phi_Scattered_495MeV_NegHelCM4", "#phi_{Sc} E_{#gamma}495 #pm 15MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 10, -180, 180);
    PhiSc525NegHelCM4 = new GH1( "Phi_Scattered_525MeV_NegHelCM4", "#phi_{Sc} E_{#gamma}525 #pm 15MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 10, -180, 180);
    PhiSc555NegHelCM4 = new GH1( "Phi_Scattered_555MeV_NegHelCM4", "#phi_{Sc} E_{#gamma}555 #pm 15MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 10, -180, 180);
    PhiSc585NegHelCM4 = new GH1( "Phi_Scattered_585MeV_NegHelCM4", "#phi_{Sc} E_{#gamma}585 #pm 15MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 10, -180, 180);
    PhiSc615NegHelCM4 = new GH1( "Phi_Scattered_615MeV_NegHelCM4", "#phi_{Sc} E_{#gamma}615 #pm 15MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 10, -180, 180);
    PhiSc645NegHelCM4 = new GH1( "Phi_Scattered_645MeV_NegHelCM4", "#phi_{Sc} E_{#gamma}645 #pm 15MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 10, -180, 180);
    PhiSc675NegHelCM4 = new GH1( "Phi_Scattered_675MeV_NegHelCM4", "#phi_{Sc} E_{#gamma}675 #pm 15MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 10, -180, 180);

    PhiSc315NegHelCM5 = new GH1( "Phi_Scattered_315MeV_NegHelCM5", "#phi_{Sc} E_{#gamma}315 #pm 15MeV (Cos#theta_{CM}0-(-0.25))) for -ve Helicity", 10, -180, 180);
    PhiSc345NegHelCM5 = new GH1( "Phi_Scattered_345MeV_NegHelCM5", "#phi_{Sc} E_{#gamma}345 #pm 15MeV (Cos#theta_{CM}0-(-0.25))) for -ve Helicity", 10, -180, 180);
    PhiSc375NegHelCM5 = new GH1( "Phi_Scattered_375MeV_NegHelCM5", "#phi_{Sc} E_{#gamma}375 #pm 15MeV (Cos#theta_{CM}0-(-0.25))) for -ve Helicity", 10, -180, 180);
    PhiSc405NegHelCM5 = new GH1( "Phi_Scattered_405MeV_NegHelCM5", "#phi_{Sc} E_{#gamma}405 #pm 15MeV (Cos#theta_{CM}0-(-0.25))) for -ve Helicity", 10, -180, 180);
    PhiSc435NegHelCM5 = new GH1( "Phi_Scattered_435MeV_NegHelCM5", "#phi_{Sc} E_{#gamma}435 #pm 15MeV (Cos#theta_{CM}0-(-0.25))) for -ve Helicity", 10, -180, 180);
    PhiSc465NegHelCM5 = new GH1( "Phi_Scattered_465MeV_NegHelCM5", "#phi_{Sc} E_{#gamma}465 #pm 15MeV (Cos#theta_{CM}0-(-0.25))) for -ve Helicity", 10, -180, 180);
    PhiSc495NegHelCM5 = new GH1( "Phi_Scattered_495MeV_NegHelCM5", "#phi_{Sc} E_{#gamma}495 #pm 15MeV (Cos#theta_{CM}0-(-0.25))) for -ve Helicity", 10, -180, 180);
    PhiSc525NegHelCM5 = new GH1( "Phi_Scattered_525MeV_NegHelCM5", "#phi_{Sc} E_{#gamma}525 #pm 15MeV (Cos#theta_{CM}0-(-0.25))) for -ve Helicity", 10, -180, 180);
    PhiSc555NegHelCM5 = new GH1( "Phi_Scattered_555MeV_NegHelCM5", "#phi_{Sc} E_{#gamma}555 #pm 15MeV (Cos#theta_{CM}0-(-0.25))) for -ve Helicity", 10, -180, 180);
    PhiSc585NegHelCM5 = new GH1( "Phi_Scattered_585MeV_NegHelCM5", "#phi_{Sc} E_{#gamma}585 #pm 15MeV (Cos#theta_{CM}0-(-0.25))) for -ve Helicity", 10, -180, 180);
    PhiSc615NegHelCM5 = new GH1( "Phi_Scattered_615MeV_NegHelCM5", "#phi_{Sc} E_{#gamma}615 #pm 15MeV (Cos#theta_{CM}0-(-0.25))) for -ve Helicity", 10, -180, 180);
    PhiSc645NegHelCM5 = new GH1( "Phi_Scattered_645MeV_NegHelCM5", "#phi_{Sc} E_{#gamma}645 #pm 15MeV (Cos#theta_{CM}0-(-0.25))) for -ve Helicity", 10, -180, 180);
    PhiSc675NegHelCM5 = new GH1( "Phi_Scattered_675MeV_NegHelCM5", "#phi_{Sc} E_{#gamma}675 #pm 15MeV (Cos#theta_{CM}0-(-0.25))) for -ve Helicity", 10, -180, 180);

    PhiSc315NegHelCM6 = new GH1( "Phi_Scattered_315MeV_NegHelCM6", "#phi_{Sc} E_{#gamma}315 #pm 15MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 10, -180, 180);
    PhiSc345NegHelCM6 = new GH1( "Phi_Scattered_345MeV_NegHelCM6", "#phi_{Sc} E_{#gamma}345 #pm 15MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 10, -180, 180);
    PhiSc375NegHelCM6 = new GH1( "Phi_Scattered_375MeV_NegHelCM6", "#phi_{Sc} E_{#gamma}375 #pm 15MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 10, -180, 180);
    PhiSc405NegHelCM6 = new GH1( "Phi_Scattered_405MeV_NegHelCM6", "#phi_{Sc} E_{#gamma}405 #pm 15MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 10, -180, 180);
    PhiSc435NegHelCM6 = new GH1( "Phi_Scattered_435MeV_NegHelCM6", "#phi_{Sc} E_{#gamma}435 #pm 15MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 10, -180, 180);
    PhiSc465NegHelCM6 = new GH1( "Phi_Scattered_465MeV_NegHelCM6", "#phi_{Sc} E_{#gamma}465 #pm 15MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 10, -180, 180);
    PhiSc495NegHelCM6 = new GH1( "Phi_Scattered_495MeV_NegHelCM6", "#phi_{Sc} E_{#gamma}495 #pm 15MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 10, -180, 180);
    PhiSc525NegHelCM6 = new GH1( "Phi_Scattered_525MeV_NegHelCM6", "#phi_{Sc} E_{#gamma}525 #pm 15MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 10, -180, 180);
    PhiSc555NegHelCM6 = new GH1( "Phi_Scattered_555MeV_NegHelCM6", "#phi_{Sc} E_{#gamma}555 #pm 15MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 10, -180, 180);
    PhiSc585NegHelCM6 = new GH1( "Phi_Scattered_585MeV_NegHelCM6", "#phi_{Sc} E_{#gamma}585 #pm 15MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 10, -180, 180);
    PhiSc615NegHelCM6 = new GH1( "Phi_Scattered_615MeV_NegHelCM6", "#phi_{Sc} E_{#gamma}615 #pm 15MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 10, -180, 180);
    PhiSc645NegHelCM6 = new GH1( "Phi_Scattered_645MeV_NegHelCM6", "#phi_{Sc} E_{#gamma}645 #pm 15MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 10, -180, 180);
    PhiSc675NegHelCM6 = new GH1( "Phi_Scattered_675MeV_NegHelCM6", "#phi_{Sc} E_{#gamma}675 #pm 15MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 10, -180, 180);

    PhiSc315NegHelCM7 = new GH1( "Phi_Scattered_315MeV_NegHelCM7", "#phi_{Sc} E_{#gamma}315 #pm 15MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 10, -180, 180);
    PhiSc345NegHelCM7 = new GH1( "Phi_Scattered_345MeV_NegHelCM7", "#phi_{Sc} E_{#gamma}345 #pm 15MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 10, -180, 180);
    PhiSc375NegHelCM7 = new GH1( "Phi_Scattered_375MeV_NegHelCM7", "#phi_{Sc} E_{#gamma}375 #pm 15MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 10, -180, 180);
    PhiSc405NegHelCM7 = new GH1( "Phi_Scattered_405MeV_NegHelCM7", "#phi_{Sc} E_{#gamma}405 #pm 15MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 10, -180, 180);
    PhiSc435NegHelCM7 = new GH1( "Phi_Scattered_435MeV_NegHelCM7", "#phi_{Sc} E_{#gamma}435 #pm 15MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 10, -180, 180);
    PhiSc465NegHelCM7 = new GH1( "Phi_Scattered_465MeV_NegHelCM7", "#phi_{Sc} E_{#gamma}465 #pm 15MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 10, -180, 180);
    PhiSc495NegHelCM7 = new GH1( "Phi_Scattered_495MeV_NegHelCM7", "#phi_{Sc} E_{#gamma}495 #pm 15MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 10, -180, 180);
    PhiSc525NegHelCM7 = new GH1( "Phi_Scattered_525MeV_NegHelCM7", "#phi_{Sc} E_{#gamma}525 #pm 15MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 10, -180, 180);
    PhiSc555NegHelCM7 = new GH1( "Phi_Scattered_555MeV_NegHelCM7", "#phi_{Sc} E_{#gamma}555 #pm 15MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 10, -180, 180);
    PhiSc585NegHelCM7 = new GH1( "Phi_Scattered_585MeV_NegHelCM7", "#phi_{Sc} E_{#gamma}585 #pm 15MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 10, -180, 180);
    PhiSc615NegHelCM7 = new GH1( "Phi_Scattered_615MeV_NegHelCM7", "#phi_{Sc} E_{#gamma}615 #pm 15MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 10, -180, 180);
    PhiSc645NegHelCM7 = new GH1( "Phi_Scattered_645MeV_NegHelCM7", "#phi_{Sc} E_{#gamma}645 #pm 15MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 10, -180, 180);
    PhiSc675NegHelCM7 = new GH1( "Phi_Scattered_675MeV_NegHelCM7", "#phi_{Sc} E_{#gamma}675 #pm 15MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 10, -180, 180);

    PhiSc315NegHelCM8 = new GH1( "Phi_Scattered_315MeV_NegHelCM8", "#phi_{Sc} E_{#gamma}315 #pm 15MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 10, -180, 180);
    PhiSc345NegHelCM8 = new GH1( "Phi_Scattered_345MeV_NegHelCM8", "#phi_{Sc} E_{#gamma}345 #pm 15MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 10, -180, 180);
    PhiSc375NegHelCM8 = new GH1( "Phi_Scattered_375MeV_NegHelCM8", "#phi_{Sc} E_{#gamma}375 #pm 15MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 10, -180, 180);
    PhiSc405NegHelCM8 = new GH1( "Phi_Scattered_405MeV_NegHelCM8", "#phi_{Sc} E_{#gamma}405 #pm 15MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 10, -180, 180);
    PhiSc435NegHelCM8 = new GH1( "Phi_Scattered_435MeV_NegHelCM8", "#phi_{Sc} E_{#gamma}435 #pm 15MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 10, -180, 180);
    PhiSc465NegHelCM8 = new GH1( "Phi_Scattered_465MeV_NegHelCM8", "#phi_{Sc} E_{#gamma}465 #pm 15MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 10, -180, 180);
    PhiSc495NegHelCM8 = new GH1( "Phi_Scattered_495MeV_NegHelCM8", "#phi_{Sc} E_{#gamma}495 #pm 15MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 10, -180, 180);
    PhiSc525NegHelCM8 = new GH1( "Phi_Scattered_525MeV_NegHelCM8", "#phi_{Sc} E_{#gamma}525 #pm 15MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 10, -180, 180);
    PhiSc555NegHelCM8 = new GH1( "Phi_Scattered_555MeV_NegHelCM8", "#phi_{Sc} E_{#gamma}555 #pm 15MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 10, -180, 180);
    PhiSc585NegHelCM8 = new GH1( "Phi_Scattered_585MeV_NegHelCM8", "#phi_{Sc} E_{#gamma}585 #pm 15MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 10, -180, 180);
    PhiSc615NegHelCM8 = new GH1( "Phi_Scattered_615MeV_NegHelCM8", "#phi_{Sc} E_{#gamma}615 #pm 15MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 10, -180, 180);
    PhiSc645NegHelCM8 = new GH1( "Phi_Scattered_645MeV_NegHelCM8", "#phi_{Sc} E_{#gamma}645 #pm 15MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 10, -180, 180);
    PhiSc675NegHelCM8 = new GH1( "Phi_Scattered_675MeV_NegHelCM8", "#phi_{Sc} E_{#gamma}675 #pm 15MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 10, -180, 180);

    // Angles of neutron in scattered frame across EGamma bins for positive helicity
    PhiSc315PosHelCM1 = new GH1( "Phi_Scattered_315MeV_PosHelCM1", "#phi_{Sc} E_{#gamma}315 #pm 15MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 10, -180, 180);
    PhiSc345PosHelCM1 = new GH1( "Phi_Scattered_345MeV_PosHelCM1", "#phi_{Sc} E_{#gamma}345 #pm 15MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 10, -180, 180);
    PhiSc375PosHelCM1 = new GH1( "Phi_Scattered_375MeV_PosHelCM1", "#phi_{Sc} E_{#gamma}375 #pm 15MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 10, -180, 180);
    PhiSc405PosHelCM1 = new GH1( "Phi_Scattered_405MeV_PosHelCM1", "#phi_{Sc} E_{#gamma}405 #pm 15MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 10, -180, 180);
    PhiSc435PosHelCM1 = new GH1( "Phi_Scattered_435MeV_PosHelCM1", "#phi_{Sc} E_{#gamma}435 #pm 15MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 10, -180, 180);
    PhiSc465PosHelCM1 = new GH1( "Phi_Scattered_465MeV_PosHelCM1", "#phi_{Sc} E_{#gamma}465 #pm 15MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 10, -180, 180);
    PhiSc495PosHelCM1 = new GH1( "Phi_Scattered_495MeV_PosHelCM1", "#phi_{Sc} E_{#gamma}495 #pm 15MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 10, -180, 180);
    PhiSc525PosHelCM1 = new GH1( "Phi_Scattered_525MeV_PosHelCM1", "#phi_{Sc} E_{#gamma}525 #pm 15MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 10, -180, 180);
    PhiSc555PosHelCM1 = new GH1( "Phi_Scattered_555MeV_PosHelCM1", "#phi_{Sc} E_{#gamma}555 #pm 15MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 10, -180, 180);
    PhiSc585PosHelCM1 = new GH1( "Phi_Scattered_585MeV_PosHelCM1", "#phi_{Sc} E_{#gamma}585 #pm 15MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 10, -180, 180);
    PhiSc615PosHelCM1 = new GH1( "Phi_Scattered_615MeV_PosHelCM1", "#phi_{Sc} E_{#gamma}615 #pm 15MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 10, -180, 180);
    PhiSc645PosHelCM1 = new GH1( "Phi_Scattered_645MeV_PosHelCM1", "#phi_{Sc} E_{#gamma}645 #pm 15MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 10, -180, 180);
    PhiSc675PosHelCM1 = new GH1( "Phi_Scattered_675MeV_PosHelCM1", "#phi_{Sc} E_{#gamma}675 #pm 15MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 10, -180, 180);

    PhiSc315PosHelCM2 = new GH1( "Phi_Scattered_315MeV_PosHelCM2", "#phi_{Sc} E_{#gamma}315 #pm 15MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 10, -180, 180);
    PhiSc345PosHelCM2 = new GH1( "Phi_Scattered_345MeV_PosHelCM2", "#phi_{Sc} E_{#gamma}345 #pm 15MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 10, -180, 180);
    PhiSc375PosHelCM2 = new GH1( "Phi_Scattered_375MeV_PosHelCM2", "#phi_{Sc} E_{#gamma}375 #pm 15MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 10, -180, 180);
    PhiSc405PosHelCM2 = new GH1( "Phi_Scattered_405MeV_PosHelCM2", "#phi_{Sc} E_{#gamma}405 #pm 15MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 10, -180, 180);
    PhiSc435PosHelCM2 = new GH1( "Phi_Scattered_435MeV_PosHelCM2", "#phi_{Sc} E_{#gamma}435 #pm 15MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 10, -180, 180);
    PhiSc465PosHelCM2 = new GH1( "Phi_Scattered_465MeV_PosHelCM2", "#phi_{Sc} E_{#gamma}465 #pm 15MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 10, -180, 180);
    PhiSc495PosHelCM2 = new GH1( "Phi_Scattered_495MeV_PosHelCM2", "#phi_{Sc} E_{#gamma}495 #pm 15MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 10, -180, 180);
    PhiSc525PosHelCM2 = new GH1( "Phi_Scattered_525MeV_PosHelCM2", "#phi_{Sc} E_{#gamma}525 #pm 15MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 10, -180, 180);
    PhiSc555PosHelCM2 = new GH1( "Phi_Scattered_555MeV_PosHelCM2", "#phi_{Sc} E_{#gamma}555 #pm 15MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 10, -180, 180);
    PhiSc585PosHelCM2 = new GH1( "Phi_Scattered_585MeV_PosHelCM2", "#phi_{Sc} E_{#gamma}585 #pm 15MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 10, -180, 180);
    PhiSc615PosHelCM2 = new GH1( "Phi_Scattered_615MeV_PosHelCM2", "#phi_{Sc} E_{#gamma}615 #pm 15MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 10, -180, 180);
    PhiSc645PosHelCM2 = new GH1( "Phi_Scattered_645MeV_PosHelCM2", "#phi_{Sc} E_{#gamma}645 #pm 15MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 10, -180, 180);
    PhiSc675PosHelCM2 = new GH1( "Phi_Scattered_675MeV_PosHelCM2", "#phi_{Sc} E_{#gamma}675 #pm 15MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 10, -180, 180);

    PhiSc315PosHelCM3 = new GH1( "Phi_Scattered_315MeV_PosHelCM3", "#phi_{Sc} E_{#gamma}315 #pm 15MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 10, -180, 180);
    PhiSc345PosHelCM3 = new GH1( "Phi_Scattered_345MeV_PosHelCM3", "#phi_{Sc} E_{#gamma}345 #pm 15MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 10, -180, 180);
    PhiSc375PosHelCM3 = new GH1( "Phi_Scattered_375MeV_PosHelCM3", "#phi_{Sc} E_{#gamma}375 #pm 15MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 10, -180, 180);
    PhiSc405PosHelCM3 = new GH1( "Phi_Scattered_405MeV_PosHelCM3", "#phi_{Sc} E_{#gamma}405 #pm 15MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 10, -180, 180);
    PhiSc435PosHelCM3 = new GH1( "Phi_Scattered_435MeV_PosHelCM3", "#phi_{Sc} E_{#gamma}435 #pm 15MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 10, -180, 180);
    PhiSc465PosHelCM3 = new GH1( "Phi_Scattered_465MeV_PosHelCM3", "#phi_{Sc} E_{#gamma}465 #pm 15MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 10, -180, 180);
    PhiSc495PosHelCM3 = new GH1( "Phi_Scattered_495MeV_PosHelCM3", "#phi_{Sc} E_{#gamma}495 #pm 15MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 10, -180, 180);
    PhiSc525PosHelCM3 = new GH1( "Phi_Scattered_525MeV_PosHelCM3", "#phi_{Sc} E_{#gamma}525 #pm 15MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 10, -180, 180);
    PhiSc555PosHelCM3 = new GH1( "Phi_Scattered_555MeV_PosHelCM3", "#phi_{Sc} E_{#gamma}555 #pm 15MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 10, -180, 180);
    PhiSc585PosHelCM3 = new GH1( "Phi_Scattered_585MeV_PosHelCM3", "#phi_{Sc} E_{#gamma}585 #pm 15MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 10, -180, 180);
    PhiSc615PosHelCM3 = new GH1( "Phi_Scattered_615MeV_PosHelCM3", "#phi_{Sc} E_{#gamma}615 #pm 15MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 10, -180, 180);
    PhiSc645PosHelCM3 = new GH1( "Phi_Scattered_645MeV_PosHelCM3", "#phi_{Sc} E_{#gamma}645 #pm 15MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 10, -180, 180);
    PhiSc675PosHelCM3 = new GH1( "Phi_Scattered_675MeV_PosHelCM3", "#phi_{Sc} E_{#gamma}675 #pm 15MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 10, -180, 180);

    PhiSc315PosHelCM4 = new GH1( "Phi_Scattered_315MeV_PosHelCM4", "#phi_{Sc} E_{#gamma}315 #pm 15MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 10, -180, 180);
    PhiSc345PosHelCM4 = new GH1( "Phi_Scattered_345MeV_PosHelCM4", "#phi_{Sc} E_{#gamma}345 #pm 15MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 10, -180, 180);
    PhiSc375PosHelCM4 = new GH1( "Phi_Scattered_375MeV_PosHelCM4", "#phi_{Sc} E_{#gamma}375 #pm 15MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 10, -180, 180);
    PhiSc405PosHelCM4 = new GH1( "Phi_Scattered_405MeV_PosHelCM4", "#phi_{Sc} E_{#gamma}405 #pm 15MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 10, -180, 180);
    PhiSc435PosHelCM4 = new GH1( "Phi_Scattered_435MeV_PosHelCM4", "#phi_{Sc} E_{#gamma}435 #pm 15MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 10, -180, 180);
    PhiSc465PosHelCM4 = new GH1( "Phi_Scattered_465MeV_PosHelCM4", "#phi_{Sc} E_{#gamma}465 #pm 15MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 10, -180, 180);
    PhiSc495PosHelCM4 = new GH1( "Phi_Scattered_495MeV_PosHelCM4", "#phi_{Sc} E_{#gamma}495 #pm 15MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 10, -180, 180);
    PhiSc525PosHelCM4 = new GH1( "Phi_Scattered_525MeV_PosHelCM4", "#phi_{Sc} E_{#gamma}525 #pm 15MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 10, -180, 180);
    PhiSc555PosHelCM4 = new GH1( "Phi_Scattered_555MeV_PosHelCM4", "#phi_{Sc} E_{#gamma}555 #pm 15MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 10, -180, 180);
    PhiSc585PosHelCM4 = new GH1( "Phi_Scattered_585MeV_PosHelCM4", "#phi_{Sc} E_{#gamma}585 #pm 15MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 10, -180, 180);
    PhiSc615PosHelCM4 = new GH1( "Phi_Scattered_615MeV_PosHelCM4", "#phi_{Sc} E_{#gamma}615 #pm 15MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 10, -180, 180);
    PhiSc645PosHelCM4 = new GH1( "Phi_Scattered_645MeV_PosHelCM4", "#phi_{Sc} E_{#gamma}645 #pm 15MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 10, -180, 180);
    PhiSc675PosHelCM4 = new GH1( "Phi_Scattered_675MeV_PosHelCM4", "#phi_{Sc} E_{#gamma}675 #pm 15MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 10, -180, 180);

    PhiSc315PosHelCM5 = new GH1( "Phi_Scattered_315MeV_PosHelCM5", "#phi_{Sc} E_{#gamma}315 #pm 15MeV (Cos#theta_{CM}0-(-0.25))) for +ve Helicity", 10, -180, 180);
    PhiSc345PosHelCM5 = new GH1( "Phi_Scattered_345MeV_PosHelCM5", "#phi_{Sc} E_{#gamma}345 #pm 15MeV (Cos#theta_{CM}0-(-0.25))) for +ve Helicity", 10, -180, 180);
    PhiSc375PosHelCM5 = new GH1( "Phi_Scattered_375MeV_PosHelCM5", "#phi_{Sc} E_{#gamma}375 #pm 15MeV (Cos#theta_{CM}0-(-0.25))) for +ve Helicity", 10, -180, 180);
    PhiSc405PosHelCM5 = new GH1( "Phi_Scattered_405MeV_PosHelCM5", "#phi_{Sc} E_{#gamma}405 #pm 15MeV (Cos#theta_{CM}0-(-0.25))) for +ve Helicity", 10, -180, 180);
    PhiSc435PosHelCM5 = new GH1( "Phi_Scattered_435MeV_PosHelCM5", "#phi_{Sc} E_{#gamma}435 #pm 15MeV (Cos#theta_{CM}0-(-0.25))) for +ve Helicity", 10, -180, 180);
    PhiSc465PosHelCM5 = new GH1( "Phi_Scattered_465MeV_PosHelCM5", "#phi_{Sc} E_{#gamma}465 #pm 15MeV (Cos#theta_{CM}0-(-0.25))) for +ve Helicity", 10, -180, 180);
    PhiSc495PosHelCM5 = new GH1( "Phi_Scattered_495MeV_PosHelCM5", "#phi_{Sc} E_{#gamma}495 #pm 15MeV (Cos#theta_{CM}0-(-0.25))) for +ve Helicity", 10, -180, 180);
    PhiSc525PosHelCM5 = new GH1( "Phi_Scattered_525MeV_PosHelCM5", "#phi_{Sc} E_{#gamma}525 #pm 15MeV (Cos#theta_{CM}0-(-0.25))) for +ve Helicity", 10, -180, 180);
    PhiSc555PosHelCM5 = new GH1( "Phi_Scattered_555MeV_PosHelCM5", "#phi_{Sc} E_{#gamma}555 #pm 15MeV (Cos#theta_{CM}0-(-0.25))) for +ve Helicity", 10, -180, 180);
    PhiSc585PosHelCM5 = new GH1( "Phi_Scattered_585MeV_PosHelCM5", "#phi_{Sc} E_{#gamma}585 #pm 15MeV (Cos#theta_{CM}0-(-0.25))) for +ve Helicity", 10, -180, 180);
    PhiSc615PosHelCM5 = new GH1( "Phi_Scattered_615MeV_PosHelCM5", "#phi_{Sc} E_{#gamma}615 #pm 15MeV (Cos#theta_{CM}0-(-0.25))) for +ve Helicity", 10, -180, 180);
    PhiSc645PosHelCM5 = new GH1( "Phi_Scattered_645MeV_PosHelCM5", "#phi_{Sc} E_{#gamma}645 #pm 15MeV (Cos#theta_{CM}0-(-0.25))) for +ve Helicity", 10, -180, 180);
    PhiSc675PosHelCM5 = new GH1( "Phi_Scattered_675MeV_PosHelCM5", "#phi_{Sc} E_{#gamma}675 #pm 15MeV (Cos#theta_{CM}0-(-0.25))) for +ve Helicity", 10, -180, 180);

    PhiSc315PosHelCM6 = new GH1( "Phi_Scattered_315MeV_PosHelCM6", "#phi_{Sc} E_{#gamma}315 #pm 15MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 10, -180, 180);
    PhiSc345PosHelCM6 = new GH1( "Phi_Scattered_345MeV_PosHelCM6", "#phi_{Sc} E_{#gamma}345 #pm 15MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 10, -180, 180);
    PhiSc375PosHelCM6 = new GH1( "Phi_Scattered_375MeV_PosHelCM6", "#phi_{Sc} E_{#gamma}375 #pm 15MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 10, -180, 180);
    PhiSc405PosHelCM6 = new GH1( "Phi_Scattered_405MeV_PosHelCM6", "#phi_{Sc} E_{#gamma}405 #pm 15MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 10, -180, 180);
    PhiSc435PosHelCM6 = new GH1( "Phi_Scattered_435MeV_PosHelCM6", "#phi_{Sc} E_{#gamma}435 #pm 15MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 10, -180, 180);
    PhiSc465PosHelCM6 = new GH1( "Phi_Scattered_465MeV_PosHelCM6", "#phi_{Sc} E_{#gamma}465 #pm 15MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 10, -180, 180);
    PhiSc495PosHelCM6 = new GH1( "Phi_Scattered_495MeV_PosHelCM6", "#phi_{Sc} E_{#gamma}495 #pm 15MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 10, -180, 180);
    PhiSc525PosHelCM6 = new GH1( "Phi_Scattered_525MeV_PosHelCM6", "#phi_{Sc} E_{#gamma}525 #pm 15MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 10, -180, 180);
    PhiSc555PosHelCM6 = new GH1( "Phi_Scattered_555MeV_PosHelCM6", "#phi_{Sc} E_{#gamma}555 #pm 15MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 10, -180, 180);
    PhiSc585PosHelCM6 = new GH1( "Phi_Scattered_585MeV_PosHelCM6", "#phi_{Sc} E_{#gamma}585 #pm 15MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 10, -180, 180);
    PhiSc615PosHelCM6 = new GH1( "Phi_Scattered_615MeV_PosHelCM6", "#phi_{Sc} E_{#gamma}615 #pm 15MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 10, -180, 180);
    PhiSc645PosHelCM6 = new GH1( "Phi_Scattered_645MeV_PosHelCM6", "#phi_{Sc} E_{#gamma}645 #pm 15MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 10, -180, 180);
    PhiSc675PosHelCM6 = new GH1( "Phi_Scattered_675MeV_PosHelCM6", "#phi_{Sc} E_{#gamma}675 #pm 15MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 10, -180, 180);

    PhiSc315PosHelCM7 = new GH1( "Phi_Scattered_315MeV_PosHelCM7", "#phi_{Sc} E_{#gamma}315 #pm 15MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 10, -180, 180);
    PhiSc345PosHelCM7 = new GH1( "Phi_Scattered_345MeV_PosHelCM7", "#phi_{Sc} E_{#gamma}345 #pm 15MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 10, -180, 180);
    PhiSc375PosHelCM7 = new GH1( "Phi_Scattered_375MeV_PosHelCM7", "#phi_{Sc} E_{#gamma}375 #pm 15MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 10, -180, 180);
    PhiSc405PosHelCM7 = new GH1( "Phi_Scattered_405MeV_PosHelCM7", "#phi_{Sc} E_{#gamma}405 #pm 15MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 10, -180, 180);
    PhiSc435PosHelCM7 = new GH1( "Phi_Scattered_435MeV_PosHelCM7", "#phi_{Sc} E_{#gamma}435 #pm 15MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 10, -180, 180);
    PhiSc465PosHelCM7 = new GH1( "Phi_Scattered_465MeV_PosHelCM7", "#phi_{Sc} E_{#gamma}465 #pm 15MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 10, -180, 180);
    PhiSc495PosHelCM7 = new GH1( "Phi_Scattered_495MeV_PosHelCM7", "#phi_{Sc} E_{#gamma}495 #pm 15MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 10, -180, 180);
    PhiSc525PosHelCM7 = new GH1( "Phi_Scattered_525MeV_PosHelCM7", "#phi_{Sc} E_{#gamma}525 #pm 15MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 10, -180, 180);
    PhiSc555PosHelCM7 = new GH1( "Phi_Scattered_555MeV_PosHelCM7", "#phi_{Sc} E_{#gamma}555 #pm 15MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 10, -180, 180);
    PhiSc585PosHelCM7 = new GH1( "Phi_Scattered_585MeV_PosHelCM7", "#phi_{Sc} E_{#gamma}585 #pm 15MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 10, -180, 180);
    PhiSc615PosHelCM7 = new GH1( "Phi_Scattered_615MeV_PosHelCM7", "#phi_{Sc} E_{#gamma}615 #pm 15MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 10, -180, 180);
    PhiSc645PosHelCM7 = new GH1( "Phi_Scattered_645MeV_PosHelCM7", "#phi_{Sc} E_{#gamma}645 #pm 15MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 10, -180, 180);
    PhiSc675PosHelCM7 = new GH1( "Phi_Scattered_675MeV_PosHelCM7", "#phi_{Sc} E_{#gamma}675 #pm 15MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 10, -180, 180);

    PhiSc315PosHelCM8 = new GH1( "Phi_Scattered_315MeV_PosHelCM8", "#phi_{Sc} E_{#gamma}315 #pm 15MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 10, -180, 180);
    PhiSc345PosHelCM8 = new GH1( "Phi_Scattered_345MeV_PosHelCM8", "#phi_{Sc} E_{#gamma}345 #pm 15MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 10, -180, 180);
    PhiSc375PosHelCM8 = new GH1( "Phi_Scattered_375MeV_PosHelCM8", "#phi_{Sc} E_{#gamma}375 #pm 15MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 10, -180, 180);
    PhiSc405PosHelCM8 = new GH1( "Phi_Scattered_405MeV_PosHelCM8", "#phi_{Sc} E_{#gamma}405 #pm 15MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 10, -180, 180);
    PhiSc435PosHelCM8 = new GH1( "Phi_Scattered_435MeV_PosHelCM8", "#phi_{Sc} E_{#gamma}435 #pm 15MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 10, -180, 180);
    PhiSc465PosHelCM8 = new GH1( "Phi_Scattered_465MeV_PosHelCM8", "#phi_{Sc} E_{#gamma}465 #pm 15MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 10, -180, 180);
    PhiSc495PosHelCM8 = new GH1( "Phi_Scattered_495MeV_PosHelCM8", "#phi_{Sc} E_{#gamma}495 #pm 15MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 10, -180, 180);
    PhiSc525PosHelCM8 = new GH1( "Phi_Scattered_525MeV_PosHelCM8", "#phi_{Sc} E_{#gamma}525 #pm 15MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 10, -180, 180);
    PhiSc555PosHelCM8 = new GH1( "Phi_Scattered_555MeV_PosHelCM8", "#phi_{Sc} E_{#gamma}555 #pm 15MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 10, -180, 180);
    PhiSc585PosHelCM8 = new GH1( "Phi_Scattered_585MeV_PosHelCM8", "#phi_{Sc} E_{#gamma}585 #pm 15MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 10, -180, 180);
    PhiSc615PosHelCM8 = new GH1( "Phi_Scattered_615MeV_PosHelCM8", "#phi_{Sc} E_{#gamma}615 #pm 15MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 10, -180, 180);
    PhiSc645PosHelCM8 = new GH1( "Phi_Scattered_645MeV_PosHelCM8", "#phi_{Sc} E_{#gamma}645 #pm 15MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 10, -180, 180);
    PhiSc675PosHelCM8 = new GH1( "Phi_Scattered_675MeV_PosHelCM8", "#phi_{Sc} E_{#gamma}675 #pm 15MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 10, -180, 180);

    ThetanDist = new GH1 ("ThetanDist", "#theta_{n} Distribution", 200, 0, 180);
    ThetanRecDist = new GH1 ("ThetanRecDist", "Reconstructed #theta_{n} Distribution", 200, 0, 180);
    ThetanDiffDist = new GH1 ("ThetanDiffDist", "Difference Between #theta_{n} and  #theta_{nRec}", 200, -90, 90);
    ThetanDiffZp = new GH2 ("ThetanDiffZp", "Diff(#theta_{n} - #theta_{nRec}) as a Fn of Z_{p}", 200, -90, 90, 200, -100, 100);

    ThetanCorrDist = new GH1 ("ThetanCorrDist", "#theta_{nCorr} Distribution", 200, 0, 180);
    ThetanCorrDiffDist = new GH1 ("ThetanCorrDiffDist", "Difference Between #theta_{n} and #theta_{nCorr} Distribution", 200, -90, 90);
    ThetanCorrRecDiffDist = new GH1 ("ThetanCorrRecDiffDist", "Difference Between #theta_{nCorr} and  #theta_{nRec}", 200, -90, 90);
    ThetanCorrDiffZp = new GH2 ("ThetanCorrDiffZp", "Diff(#theta_{nCorr} - #theta_{nRec}) as a Fn of Z_{p}", 200, -90, 90, 200, -100, 100);

    ThetaRecPiDiff = new GH1 ("ThetaRecPiDiff", "Difference between #theta_{#pi Rec} and #theta_{nCorr}", 200, -90, 90);
    ThetanThetaRecPi = new GH2 ("ThetanThetaRecPi", "#theta_{nCorr} vs #theta_{#pi rec}", 100, 0, 180, 100, 0, 180);
    ThetanThetaRecPiDiff = new GH2 ("ThetanThetaRecPiDiff", "#theta_{nCorr} vs (#theta_{#pi Rec} - #theta_{nCorr})", 100, 0, 180, 100, -90, 90);

    ThetaRecPDiff = new GH1 ("ThetaRecPDiff", "Difference between #theta_{pRec} and #theta_{nCorr}", 200, -90, 90);
    ThetanThetaRecP = new GH2 ("ThetanThetaRecP", "#theta_{nCorr} vs #theta_{pRec}", 100, 0, 180, 100, 0, 180);
    ThetanThetaRecPDiff = new GH2 ("ThetanThetaRecPDiff", "#theta_{nCorr} vs (#theta_{pRec} - #theta_{nCorr})", 100, 0, 180, 100, -90, 90);

    DeutKinPiKin = new GH2 ("DeutKinPiKin", "(#theta_{nRec} - #theta_{nCorr}) vs (#theta_{#pi Rec} - #theta_{nCorr})", 200, -180, 180, 200, -180, 180);

    E_dE = new GH2 ("E_dE", "EdE Plot With E Loss Adjustment", 100, 0, 500, 100, 0, 5);
    KinEp_dE = new GH2 ("KinEp_dE", "KinEpdE Plot", 100, 0, 500, 100, 0, 5);
    PhiDiffThetaSc = new GH2("PhiDiffThetaSc", "#theta_{Sc} as a fn of #phi_{Diff}", 100, 0, 360, 100, 0, 180);
    //KinEp_dE_GoodCut = new GH2 ("KinEp_dE_GoodCut", "KinEpdE Plot With Good Proton Cut", 100, 0, 500, 100, 0, 5);
    ThetaScPhiSc = new GH2 ("ThetaScPhiSc", "#Phi_{Sc} as a function of #theta_{Sc}", 100, 0, 180, 100, -180, 180);
    EpCorr_KinEp = new GH2 ("EpCorr_KinEp", "E_{pKin} as a fn of E_{pCorr}", 100, 0, 500, 100, 0, 500);
    PhinDiffWCZRec = new GH2 ("PhinDiffWCZRec", "Difference between WC Phi and Reconstructed Phi as a fn of WCZ Hit Position", 100, 0, 200, 100, 0, 180);
    ThetaDiffPhiDiff = new GH2 ("ThetaDiffPhiDiff", "PhiDiff as a Fn of ThetaDiff (Detected-Rec)", 100, -180, 180, 100, -180, 180);

    ClosestApproach = new GH1("ClosestApproach", "DOCA of n and p' vectors", 200, -200, 200);
    POCAr = new GH1("POCAr", "r_{POCA}", 200, 0, 300);
    ScatterVertexZ = new GH1("ScatterVertexZ", "Z_{POCA}", 200, -200, 200);
    ScatterVertexZr = new GH2("ScatterVertexZr", "Z_{POCA} vs r_{POCA}", 200, -200, 200, 200, 0, 200);
    ScatterVertexXY = new GH2("ScatterVertexXY", "XY Vertex Point of Scatter from DOCA Method", 100, -80, 80, 100, -80, 80);
    ScatterVertex = new GH3("ScatterVertex", "Vertex Point of Scatter from DOCA Method", 100, -80, 80, 100, -80, 80, 100, -200, 200);
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
    PhiDiffThetaSc->Fill(PhiDiff, ScattTheta, TaggerTime);

    PhiDifference->Fill(PhiDiff);
    EpCorr_KinEp->Fill(EpCorr, KinEp, TaggerTime);
    PhinDiffWCZRec->Fill(WCZnRec, PhinDiff, TaggerTime);

    ThetanDist->Fill(Thn, TaggerTime);
    ThetanRecDist->Fill(ThetanRec, TaggerTime);
    ThetanDiffDist->Fill(Thetan-ThetanRec, TaggerTime);
    ThetanDiffZp->Fill(Thn-ThetanRec, Zp, TaggerTime);

    ThetanCorrDist->Fill(ThetanCorr, TaggerTime);
    ThetanCorrDiffDist->Fill(Thn-ThetanCorr, TaggerTime);
    ThetanCorrRecDiffDist->Fill(ThetanCorr-ThetanRec, TaggerTime);
    ThetanCorrDiffZp ->Fill(ThetanCorr-ThetanRec, Zp, TaggerTime);

    ThetaRecPiDiff->Fill(ThetaPiRecDiff, TaggerTime);
    ThetanThetaRecPi->Fill(ThetanCorr, ThetaPiRec, TaggerTime);
    ThetanThetaRecPiDiff->Fill(ThetanCorr, ThetaPiRecDiff, TaggerTime);

    ThetaRecPDiff->Fill(ThetapRecDiff, TaggerTime);
    ThetanThetaRecP->Fill(ThetanCorr, ThetapRec, TaggerTime);
    ThetanThetaRecPDiff->Fill(ThetanCorr, ThetapRecDiff, TaggerTime);

    ThetaDiffPhiDiff->Fill((ThetanCorr-ThetanRec), (Phn-PhinRec), TaggerTime);

    DeutKinPiKin->Fill(ThetanRec-ThetanCorr, ThetaPiRecDiff, TaggerTime);

    ClosestApproach->Fill(DOCA, TaggerTime);
    POCAr->Fill(r, TaggerTime);
    ScatterVertexZ->Fill(POCAz, TaggerTime);
    ScatterVertexZr->Fill(POCAz, r, TaggerTime);
    ScatterVertexXY->Fill(POCAx, POCAy, TaggerTime);
    ScatterVertex->Fill(POCAx, POCAy, POCAz, TaggerTime);

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


    if ( 300 < EGamma && EGamma < 330) {

        if(1 > CosThetapCM && CosThetapCM > 0.75 ){
            if (BeamHelicity == kFALSE) PhiSc315NegHelCM1->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc315PosHelCM1->Fill(ScattPhi, TaggerTime);
        }

        else if(0.75 > CosThetapCM && CosThetapCM > 0.5){
            if (BeamHelicity == kFALSE) PhiSc315NegHelCM2->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc315PosHelCM2->Fill(ScattPhi, TaggerTime);
        }

        else if(0.5 > CosThetapCM && CosThetapCM > 0.25){
            if (BeamHelicity == kFALSE) PhiSc315NegHelCM3->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc315PosHelCM3->Fill(ScattPhi, TaggerTime);
        }

        else if(0.5 > CosThetapCM && CosThetapCM > 0.0){
            if (BeamHelicity == kFALSE) PhiSc315NegHelCM4->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc315PosHelCM4->Fill(ScattPhi, TaggerTime);
        }

        else if(0.0 > CosThetapCM && CosThetapCM > -0.25){
            if (BeamHelicity == kFALSE) PhiSc315NegHelCM5->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc315PosHelCM5->Fill(ScattPhi, TaggerTime);
        }

        else if(-0.25 > CosThetapCM && CosThetapCM > -0.5){
            if (BeamHelicity == kFALSE) PhiSc315NegHelCM6->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc315PosHelCM6->Fill(ScattPhi, TaggerTime);
        }

        else if(-0.5 > CosThetapCM && CosThetapCM > -0.75){
            if (BeamHelicity == kFALSE) PhiSc315NegHelCM7->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc315PosHelCM7->Fill(ScattPhi, TaggerTime);
        }

        else if(-0.75 > CosThetapCM && CosThetapCM > -1){
            if (BeamHelicity == kFALSE) PhiSc315NegHelCM8->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc315PosHelCM8->Fill(ScattPhi, TaggerTime);
        }
    }

    if ( 330 < EGamma && EGamma < 360) {

        if(1 > CosThetapCM && CosThetapCM > 0.75 ){
            if (BeamHelicity == kFALSE) PhiSc345NegHelCM1->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc345PosHelCM1->Fill(ScattPhi, TaggerTime);
        }

        else if(0.75 > CosThetapCM && CosThetapCM > 0.5){
            if (BeamHelicity == kFALSE) PhiSc345NegHelCM2->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc345PosHelCM2->Fill(ScattPhi, TaggerTime);
        }

        else if(0.5 > CosThetapCM && CosThetapCM > 0.25){
            if (BeamHelicity == kFALSE) PhiSc345NegHelCM3->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc345PosHelCM3->Fill(ScattPhi, TaggerTime);
        }

        else if(0.5 > CosThetapCM && CosThetapCM > 0.0){
            if (BeamHelicity == kFALSE) PhiSc345NegHelCM4->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc345PosHelCM4->Fill(ScattPhi, TaggerTime);
        }

        else if(0.0 > CosThetapCM && CosThetapCM > -0.25){
            if (BeamHelicity == kFALSE) PhiSc345NegHelCM5->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc345PosHelCM5->Fill(ScattPhi, TaggerTime);
        }

        else if(-0.25 > CosThetapCM && CosThetapCM > -0.5){
            if (BeamHelicity == kFALSE) PhiSc345NegHelCM6->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc345PosHelCM6->Fill(ScattPhi, TaggerTime);
        }

        else if(-0.5 > CosThetapCM && CosThetapCM > -0.75){
            if (BeamHelicity == kFALSE) PhiSc345NegHelCM7->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc345PosHelCM7->Fill(ScattPhi, TaggerTime);
        }

        else if(-0.75 > CosThetapCM && CosThetapCM > -1){
            if (BeamHelicity == kFALSE) PhiSc345NegHelCM8->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc345PosHelCM8->Fill(ScattPhi, TaggerTime);
        }
    }

    if ( 360 < EGamma && EGamma < 390) {

        if(1 > CosThetapCM && CosThetapCM > 0.75 ){
            if (BeamHelicity == kFALSE) PhiSc375NegHelCM1->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc375PosHelCM1->Fill(ScattPhi, TaggerTime);
        }

        else if(0.75 > CosThetapCM && CosThetapCM > 0.5){
            if (BeamHelicity == kFALSE) PhiSc375NegHelCM2->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc375PosHelCM2->Fill(ScattPhi, TaggerTime);
        }

        else if(0.5 > CosThetapCM && CosThetapCM > 0.25){
            if (BeamHelicity == kFALSE) PhiSc375NegHelCM3->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc375PosHelCM3->Fill(ScattPhi, TaggerTime);
        }

        else if(0.5 > CosThetapCM && CosThetapCM > 0.0){
            if (BeamHelicity == kFALSE) PhiSc375NegHelCM4->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc375PosHelCM4->Fill(ScattPhi, TaggerTime);
        }

        else if(0.0 > CosThetapCM && CosThetapCM > -0.25){
            if (BeamHelicity == kFALSE) PhiSc375NegHelCM5->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc375PosHelCM5->Fill(ScattPhi, TaggerTime);
        }

        else if(-0.25 > CosThetapCM && CosThetapCM > -0.5){
            if (BeamHelicity == kFALSE) PhiSc375NegHelCM6->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc375PosHelCM6->Fill(ScattPhi, TaggerTime);
        }

        else if(-0.5 > CosThetapCM && CosThetapCM > -0.75){
            if (BeamHelicity == kFALSE) PhiSc375NegHelCM7->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc375PosHelCM7->Fill(ScattPhi, TaggerTime);
        }

        else if(-0.75 > CosThetapCM && CosThetapCM > -1){
            if (BeamHelicity == kFALSE) PhiSc375NegHelCM8->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc375PosHelCM8->Fill(ScattPhi, TaggerTime);
        }
    }

    if ( 390 < EGamma && EGamma < 420) {

        if(1 > CosThetapCM && CosThetapCM > 0.75 ){
            if (BeamHelicity == kFALSE) PhiSc405NegHelCM1->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc405PosHelCM1->Fill(ScattPhi, TaggerTime);
        }

        else if(0.75 > CosThetapCM && CosThetapCM > 0.5){
            if (BeamHelicity == kFALSE) PhiSc405NegHelCM2->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc405PosHelCM2->Fill(ScattPhi, TaggerTime);
        }

        else if(0.5 > CosThetapCM && CosThetapCM > 0.25){
            if (BeamHelicity == kFALSE) PhiSc405NegHelCM3->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc405PosHelCM3->Fill(ScattPhi, TaggerTime);
        }

        else if(0.5 > CosThetapCM && CosThetapCM > 0.0){
            if (BeamHelicity == kFALSE) PhiSc405NegHelCM4->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc405PosHelCM4->Fill(ScattPhi, TaggerTime);
        }

        else if(0.0 > CosThetapCM && CosThetapCM > -0.25){
            if (BeamHelicity == kFALSE) PhiSc405NegHelCM5->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc405PosHelCM5->Fill(ScattPhi, TaggerTime);
        }

        else if(-0.25 > CosThetapCM && CosThetapCM > -0.5){
            if (BeamHelicity == kFALSE) PhiSc405NegHelCM6->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc405PosHelCM6->Fill(ScattPhi, TaggerTime);
        }

        else if(-0.5 > CosThetapCM && CosThetapCM > -0.75){
            if (BeamHelicity == kFALSE) PhiSc405NegHelCM7->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc405PosHelCM7->Fill(ScattPhi, TaggerTime);
        }

        else if(-0.75 > CosThetapCM && CosThetapCM > -1){
            if (BeamHelicity == kFALSE) PhiSc405NegHelCM8->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc405PosHelCM8->Fill(ScattPhi, TaggerTime);
        }
    }

    if ( 420 < EGamma && EGamma < 450) {

        if(1 > CosThetapCM && CosThetapCM > 0.75 ){
            if (BeamHelicity == kFALSE) PhiSc435NegHelCM1->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc435PosHelCM1->Fill(ScattPhi, TaggerTime);
        }

        else if(0.75 > CosThetapCM && CosThetapCM > 0.5){
            if (BeamHelicity == kFALSE) PhiSc435NegHelCM2->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc435PosHelCM2->Fill(ScattPhi, TaggerTime);
        }

        else if(0.5 > CosThetapCM && CosThetapCM > 0.25){
            if (BeamHelicity == kFALSE) PhiSc435NegHelCM3->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc435PosHelCM3->Fill(ScattPhi, TaggerTime);
        }

        else if(0.5 > CosThetapCM && CosThetapCM > 0.0){
            if (BeamHelicity == kFALSE) PhiSc435NegHelCM4->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc435PosHelCM4->Fill(ScattPhi, TaggerTime);
        }

        else if(0.0 > CosThetapCM && CosThetapCM > -0.25){
            if (BeamHelicity == kFALSE) PhiSc435NegHelCM5->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc435PosHelCM5->Fill(ScattPhi, TaggerTime);
        }

        else if(-0.25 > CosThetapCM && CosThetapCM > -0.5){
            if (BeamHelicity == kFALSE) PhiSc435NegHelCM6->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc435PosHelCM6->Fill(ScattPhi, TaggerTime);
        }

        else if(-0.5 > CosThetapCM && CosThetapCM > -0.75){
            if (BeamHelicity == kFALSE) PhiSc435NegHelCM7->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc435PosHelCM7->Fill(ScattPhi, TaggerTime);
        }

        else if(-0.75 > CosThetapCM && CosThetapCM > -1){
            if (BeamHelicity == kFALSE) PhiSc435NegHelCM8->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc435PosHelCM8->Fill(ScattPhi, TaggerTime);
        }
    }

    if ( 450 < EGamma && EGamma < 480) {

        if(1 > CosThetapCM && CosThetapCM > 0.75 ){
            if (BeamHelicity == kFALSE) PhiSc465NegHelCM1->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc465PosHelCM1->Fill(ScattPhi, TaggerTime);
        }

        else if(0.75 > CosThetapCM && CosThetapCM > 0.5){
            if (BeamHelicity == kFALSE) PhiSc465NegHelCM2->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc465PosHelCM2->Fill(ScattPhi, TaggerTime);
        }

        else if(0.5 > CosThetapCM && CosThetapCM > 0.25){
            if (BeamHelicity == kFALSE) PhiSc465NegHelCM3->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc465PosHelCM3->Fill(ScattPhi, TaggerTime);
        }

        else if(0.5 > CosThetapCM && CosThetapCM > 0.0){
            if (BeamHelicity == kFALSE) PhiSc465NegHelCM4->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc465PosHelCM4->Fill(ScattPhi, TaggerTime);
        }

        else if(0.0 > CosThetapCM && CosThetapCM > -0.25){
            if (BeamHelicity == kFALSE) PhiSc465NegHelCM5->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc465PosHelCM5->Fill(ScattPhi, TaggerTime);
        }

        else if(-0.25 > CosThetapCM && CosThetapCM > -0.5){
            if (BeamHelicity == kFALSE) PhiSc465NegHelCM6->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc465PosHelCM6->Fill(ScattPhi, TaggerTime);
        }

        else if(-0.5 > CosThetapCM && CosThetapCM > -0.75){
            if (BeamHelicity == kFALSE) PhiSc465NegHelCM7->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc465PosHelCM7->Fill(ScattPhi, TaggerTime);
        }

        else if(-0.75 > CosThetapCM && CosThetapCM > -1){
            if (BeamHelicity == kFALSE) PhiSc465NegHelCM8->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc465PosHelCM8->Fill(ScattPhi, TaggerTime);
        }
    }

    if ( 480 < EGamma && EGamma < 510) {

        if(1 > CosThetapCM && CosThetapCM > 0.75 ){
            if (BeamHelicity == kFALSE) PhiSc495NegHelCM1->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc495PosHelCM1->Fill(ScattPhi, TaggerTime);
        }

        else if(0.75 > CosThetapCM && CosThetapCM > 0.5){
            if (BeamHelicity == kFALSE) PhiSc495NegHelCM2->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc495PosHelCM2->Fill(ScattPhi, TaggerTime);
        }

        else if(0.5 > CosThetapCM && CosThetapCM > 0.25){
            if (BeamHelicity == kFALSE) PhiSc495NegHelCM3->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc495PosHelCM3->Fill(ScattPhi, TaggerTime);
        }

        else if(0.5 > CosThetapCM && CosThetapCM > 0.0){
            if (BeamHelicity == kFALSE) PhiSc495NegHelCM4->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc495PosHelCM4->Fill(ScattPhi, TaggerTime);
        }

        else if(0.0 > CosThetapCM && CosThetapCM > -0.25){
            if (BeamHelicity == kFALSE) PhiSc495NegHelCM5->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc495PosHelCM5->Fill(ScattPhi, TaggerTime);
        }

        else if(-0.25 > CosThetapCM && CosThetapCM > -0.5){
            if (BeamHelicity == kFALSE) PhiSc495NegHelCM6->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc495PosHelCM6->Fill(ScattPhi, TaggerTime);
        }

        else if(-0.5 > CosThetapCM && CosThetapCM > -0.75){
            if (BeamHelicity == kFALSE) PhiSc495NegHelCM7->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc495PosHelCM7->Fill(ScattPhi, TaggerTime);
        }

        else if(-0.75 > CosThetapCM && CosThetapCM > -1){
            if (BeamHelicity == kFALSE) PhiSc495NegHelCM8->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc495PosHelCM8->Fill(ScattPhi, TaggerTime);
        }
    }

    if ( 510 < EGamma && EGamma < 540) {

        if(1 > CosThetapCM && CosThetapCM > 0.75 ){
            if (BeamHelicity == kFALSE) PhiSc525NegHelCM1->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc525PosHelCM1->Fill(ScattPhi, TaggerTime);
        }

        else if(0.75 > CosThetapCM && CosThetapCM > 0.5){
            if (BeamHelicity == kFALSE) PhiSc525NegHelCM2->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc525PosHelCM2->Fill(ScattPhi, TaggerTime);
        }

        else if(0.5 > CosThetapCM && CosThetapCM > 0.25){
            if (BeamHelicity == kFALSE) PhiSc525NegHelCM3->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc525PosHelCM3->Fill(ScattPhi, TaggerTime);
        }

        else if(0.5 > CosThetapCM && CosThetapCM > 0.0){
            if (BeamHelicity == kFALSE) PhiSc525NegHelCM4->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc525PosHelCM4->Fill(ScattPhi, TaggerTime);
        }

        else if(0.0 > CosThetapCM && CosThetapCM > -0.25){
            if (BeamHelicity == kFALSE) PhiSc525NegHelCM5->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc525PosHelCM5->Fill(ScattPhi, TaggerTime);
        }

        else if(-0.25 > CosThetapCM && CosThetapCM > -0.5){
            if (BeamHelicity == kFALSE) PhiSc525NegHelCM6->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc525PosHelCM6->Fill(ScattPhi, TaggerTime);
        }

        else if(-0.5 > CosThetapCM && CosThetapCM > -0.75){
            if (BeamHelicity == kFALSE) PhiSc525NegHelCM7->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc525PosHelCM7->Fill(ScattPhi, TaggerTime);
        }

        else if(-0.75 > CosThetapCM && CosThetapCM > -1){
            if (BeamHelicity == kFALSE) PhiSc525NegHelCM8->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc525PosHelCM8->Fill(ScattPhi, TaggerTime);
        }
    }

    if ( 540 < EGamma && EGamma < 570) {

        if(1 > CosThetapCM && CosThetapCM > 0.75 ){
            if (BeamHelicity == kFALSE) PhiSc555NegHelCM1->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc555PosHelCM1->Fill(ScattPhi, TaggerTime);
        }

        else if(0.75 > CosThetapCM && CosThetapCM > 0.5){
            if (BeamHelicity == kFALSE) PhiSc555NegHelCM2->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc555PosHelCM2->Fill(ScattPhi, TaggerTime);
        }

        else if(0.5 > CosThetapCM && CosThetapCM > 0.25){
            if (BeamHelicity == kFALSE) PhiSc555NegHelCM3->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc555PosHelCM3->Fill(ScattPhi, TaggerTime);
        }

        else if(0.5 > CosThetapCM && CosThetapCM > 0.0){
            if (BeamHelicity == kFALSE) PhiSc555NegHelCM4->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc555PosHelCM4->Fill(ScattPhi, TaggerTime);
        }

        else if(0.0 > CosThetapCM && CosThetapCM > -0.25){
            if (BeamHelicity == kFALSE) PhiSc555NegHelCM5->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc555PosHelCM5->Fill(ScattPhi, TaggerTime);
        }

        else if(-0.25 > CosThetapCM && CosThetapCM > -0.5){
            if (BeamHelicity == kFALSE) PhiSc555NegHelCM6->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc555PosHelCM6->Fill(ScattPhi, TaggerTime);
        }

        else if(-0.5 > CosThetapCM && CosThetapCM > -0.75){
            if (BeamHelicity == kFALSE) PhiSc555NegHelCM7->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc555PosHelCM7->Fill(ScattPhi, TaggerTime);
        }

        else if(-0.75 > CosThetapCM && CosThetapCM > -1){
            if (BeamHelicity == kFALSE) PhiSc555NegHelCM8->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc555PosHelCM8->Fill(ScattPhi, TaggerTime);
        }
    }

    if ( 570 < EGamma && EGamma < 600) {

        if(1 > CosThetapCM && CosThetapCM > 0.75 ){
            if (BeamHelicity == kFALSE) PhiSc585NegHelCM1->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc585PosHelCM1->Fill(ScattPhi, TaggerTime);
        }

        else if(0.75 > CosThetapCM && CosThetapCM > 0.5){
            if (BeamHelicity == kFALSE) PhiSc585NegHelCM2->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc585PosHelCM2->Fill(ScattPhi, TaggerTime);
        }

        else if(0.5 > CosThetapCM && CosThetapCM > 0.25){
            if (BeamHelicity == kFALSE) PhiSc585NegHelCM3->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc585PosHelCM3->Fill(ScattPhi, TaggerTime);
        }

        else if(0.5 > CosThetapCM && CosThetapCM > 0.0){
            if (BeamHelicity == kFALSE) PhiSc585NegHelCM4->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc585PosHelCM4->Fill(ScattPhi, TaggerTime);
        }

        else if(0.0 > CosThetapCM && CosThetapCM > -0.25){
            if (BeamHelicity == kFALSE) PhiSc585NegHelCM5->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc585PosHelCM5->Fill(ScattPhi, TaggerTime);
        }

        else if(-0.25 > CosThetapCM && CosThetapCM > -0.5){
            if (BeamHelicity == kFALSE) PhiSc585NegHelCM6->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc585PosHelCM6->Fill(ScattPhi, TaggerTime);
        }

        else if(-0.5 > CosThetapCM && CosThetapCM > -0.75){
            if (BeamHelicity == kFALSE) PhiSc585NegHelCM7->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc585PosHelCM7->Fill(ScattPhi, TaggerTime);
        }

        else if(-0.75 > CosThetapCM && CosThetapCM > -1){
            if (BeamHelicity == kFALSE) PhiSc585NegHelCM8->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc585PosHelCM8->Fill(ScattPhi, TaggerTime);
        }
    }

    if ( 600 < EGamma && EGamma < 630) {

        if(1 > CosThetapCM && CosThetapCM > 0.75 ){
            if (BeamHelicity == kFALSE) PhiSc615NegHelCM1->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc615PosHelCM1->Fill(ScattPhi, TaggerTime);
        }

        else if(0.75 > CosThetapCM && CosThetapCM > 0.5){
            if (BeamHelicity == kFALSE) PhiSc615NegHelCM2->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc615PosHelCM2->Fill(ScattPhi, TaggerTime);
        }

        else if(0.5 > CosThetapCM && CosThetapCM > 0.25){
            if (BeamHelicity == kFALSE) PhiSc615NegHelCM3->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc615PosHelCM3->Fill(ScattPhi, TaggerTime);
        }

        else if(0.5 > CosThetapCM && CosThetapCM > 0.0){
            if (BeamHelicity == kFALSE) PhiSc615NegHelCM4->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc615PosHelCM4->Fill(ScattPhi, TaggerTime);
        }

        else if(0.0 > CosThetapCM && CosThetapCM > -0.25){
            if (BeamHelicity == kFALSE) PhiSc615NegHelCM5->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc615PosHelCM5->Fill(ScattPhi, TaggerTime);
        }

        else if(-0.25 > CosThetapCM && CosThetapCM > -0.5){
            if (BeamHelicity == kFALSE) PhiSc615NegHelCM6->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc615PosHelCM6->Fill(ScattPhi, TaggerTime);
        }

        else if(-0.5 > CosThetapCM && CosThetapCM > -0.75){
            if (BeamHelicity == kFALSE) PhiSc615NegHelCM7->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc615PosHelCM7->Fill(ScattPhi, TaggerTime);
        }

        else if(-0.75 > CosThetapCM && CosThetapCM > -1){
            if (BeamHelicity == kFALSE) PhiSc615NegHelCM8->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc615PosHelCM8->Fill(ScattPhi, TaggerTime);
        }
    }

    if ( 630 < EGamma && EGamma < 660) {

        if(1 > CosThetapCM && CosThetapCM > 0.75 ){
            if (BeamHelicity == kFALSE) PhiSc645NegHelCM1->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc645PosHelCM1->Fill(ScattPhi, TaggerTime);
        }

        else if(0.75 > CosThetapCM && CosThetapCM > 0.5){
            if (BeamHelicity == kFALSE) PhiSc645NegHelCM2->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc645PosHelCM2->Fill(ScattPhi, TaggerTime);
        }

        else if(0.5 > CosThetapCM && CosThetapCM > 0.25){
            if (BeamHelicity == kFALSE) PhiSc645NegHelCM3->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc645PosHelCM3->Fill(ScattPhi, TaggerTime);
        }

        else if(0.5 > CosThetapCM && CosThetapCM > 0.0){
            if (BeamHelicity == kFALSE) PhiSc645NegHelCM4->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc645PosHelCM4->Fill(ScattPhi, TaggerTime);
        }

        else if(0.0 > CosThetapCM && CosThetapCM > -0.25){
            if (BeamHelicity == kFALSE) PhiSc645NegHelCM5->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc645PosHelCM5->Fill(ScattPhi, TaggerTime);
        }

        else if(-0.25 > CosThetapCM && CosThetapCM > -0.5){
            if (BeamHelicity == kFALSE) PhiSc645NegHelCM6->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc645PosHelCM6->Fill(ScattPhi, TaggerTime);
        }

        else if(-0.5 > CosThetapCM && CosThetapCM > -0.75){
            if (BeamHelicity == kFALSE) PhiSc645NegHelCM7->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc645PosHelCM7->Fill(ScattPhi, TaggerTime);
        }

        else if(-0.75 > CosThetapCM && CosThetapCM > -1){
            if (BeamHelicity == kFALSE) PhiSc645NegHelCM8->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc645PosHelCM8->Fill(ScattPhi, TaggerTime);
        }
    }

    if ( 660 < EGamma && EGamma < 690) {

        if(1 > CosThetapCM && CosThetapCM > 0.75 ){
            if (BeamHelicity == kFALSE) PhiSc675NegHelCM1->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc675PosHelCM1->Fill(ScattPhi, TaggerTime);
        }

        else if(0.75 > CosThetapCM && CosThetapCM > 0.5){
            if (BeamHelicity == kFALSE) PhiSc675NegHelCM2->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc675PosHelCM2->Fill(ScattPhi, TaggerTime);
        }

        else if(0.5 > CosThetapCM && CosThetapCM > 0.25){
            if (BeamHelicity == kFALSE) PhiSc675NegHelCM3->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc675PosHelCM3->Fill(ScattPhi, TaggerTime);
        }

        else if(0.5 > CosThetapCM && CosThetapCM > 0.0){
            if (BeamHelicity == kFALSE) PhiSc675NegHelCM4->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc675PosHelCM4->Fill(ScattPhi, TaggerTime);
        }

        else if(0.0 > CosThetapCM && CosThetapCM > -0.25){
            if (BeamHelicity == kFALSE) PhiSc675NegHelCM5->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc675PosHelCM5->Fill(ScattPhi, TaggerTime);
        }

        else if(-0.25 > CosThetapCM && CosThetapCM > -0.5){
            if (BeamHelicity == kFALSE) PhiSc675NegHelCM6->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc675PosHelCM6->Fill(ScattPhi, TaggerTime);
        }

        else if(-0.5 > CosThetapCM && CosThetapCM > -0.75){
            if (BeamHelicity == kFALSE) PhiSc675NegHelCM7->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc675PosHelCM7->Fill(ScattPhi, TaggerTime);
        }

        else if(-0.75 > CosThetapCM && CosThetapCM > -1){
            if (BeamHelicity == kFALSE) PhiSc675NegHelCM8->Fill(ScattPhi, TaggerTime);
            else if (BeamHelicity == kTRUE) PhiSc675PosHelCM8->Fill(ScattPhi, TaggerTime);
        }
    }
}

Bool_t	PNeutPol_Polarimeter_Circ::Write(){
  // Write all GH1's easily

  GTreeManager::Write();
}
