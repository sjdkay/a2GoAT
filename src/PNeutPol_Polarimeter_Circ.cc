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
    RandomLow1 = -135;
    RandomHigh1 = -35;
    RandomLow2 = 35;
    RandomHigh2 = 135;
    PvRratio = (PromptHigh - PromptLow)/( (RandomHigh1 - RandomLow1) + (RandomHigh2 - RandomLow2));

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
    if( nVertex(2) > 60 || nVertex(2) < -60) return;
    EpCorr = EpPolCorrect(Ep, Thp); //correct Ep for energy loss in polarimeter

    if(Cut_proton -> IsInside(EpCorr, dEp) == kFALSE) return; // If E loss correct proton is NOT inside p banana drop out
    if(MWPC0pE  < 150) return;
    if(MWPC0nE  < 150) return;
    if(MWPC1pE  < 150) return;
    if(MWPC1nE  < 150) return;
    if((MWPC0nE + MWPC1nE) > 700) return;

    EpDiff = abs(EpCorr - Ep);

    Pp = sqrt (TMath::Power((Mp + EpCorr),2) - TMath::Power(Mp,2));
    Pn = sqrt (TMath::Power((En + Mn ),2) - TMath::Power(Mn,2));
    GVpCorr = TLorentzVector(Pp*sin(Thp)*cos(Php), Pp*sin(Thp)*sin(Php), Pp*cos(Thp), EpCorr+Mp);

    // MUST FEED IN n PHI IN RAD NOT DEG!
    GVnCorr =  CNeutron4VectorCorr(pVertex(2), GVn, En, Pn , Mn, PhnRad);
    ThetanCorr = (GVnCorr.Theta())*TMath::RadToDeg();

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

        // Gamma(d,p)n calc n from kinematics
        KinEp = CalcKinEnergy(Thp, EGamma, Md, 0., Mp, Mn); // Calculate kin E of proton assuming pn production
        pKin = CProton4VectorKin(KinEp, ThpRad, PhpRad);
        nKin = CNeutron4VectorKin(pKin);
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
        if (80 > ThetanCM || ThetanCM > 100) continue;

        RecNeutronEpCorr = CNeutron4VectorKin(GVpCorr);
        MMpEpCorr = RecNeutronEpCorr.M();

        if (((MMpEpCorr < 800) == kTRUE) || ((MMpEpCorr > 1100) == kTRUE)) continue; // Force a missing mass cut

        KinEDiff = KinEp - EpCorr;

        ThetanDiff = ThetanCorr-ThetanRec;
        PhinDiff = Phn-PhinRec;

        TVector3 ScattAngles = ScatteredFrameAngles(nKin3, pKin3 , GVn3 , Gamma);
        ScattThetaRad = ScattAngles(2);
        ScattTheta = ScattThetaRad*TMath::RadToDeg();
        ScattPhiRad = ScattAngles(0);
        ScattPhi = ScattPhiRad*TMath::RadToDeg();
        if (ScattPhiRad > -0.01 && ScattPhiRad < 0.01) continue; // Kill excessively small values
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
        if(ScattTheta > 45) return;

        //APLCON
        beamF.SetFromVector(Gamma); // Set Lorentz vectors for use in APLCON
        protonF.SetFromVector(GVpCorr);
        neutronF.SetFromVector(GVn);

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

    PhiSc = new TH1D( "Phi_Scattered", "#phi_{sc} Proton Distribution", 200, -4, 4);
    PhiScPrompt = new TH1D( "Phi_Scattered", "#phi_{sc} Proton Distribution", 200, -4, 4);
    PhiScRandom = new TH1D( "Phi_Scattered", "#phi_{sc} Proton Distribution", 200, -4, 4);

    PhiSc320 = new TH1D("PhiSc320", "#phi_{Sc} (300-340MeV)", 24, -4, 4);
    PhiSc360 = new TH1D("PhiSc360", "#phi_{Sc} (340-380MeV)", 24, -4, 4);
    PhiSc400 = new TH1D("PhiSc400", "#phi_{Sc} (380-420MeV)", 24, -4, 4);
    PhiSc440 = new TH1D("PhiSc440", "#phi_{Sc} (420-460MeV)", 24, -4, 4);
    PhiSc480 = new TH1D("PhiSc480", "#phi_{Sc} (460-500MeV)", 24, -4, 4);
    PhiSc520 = new TH1D("PhiSc520", "#phi_{Sc} (500-540MeV)", 24, -4, 4);
    PhiSc560 = new TH1D("PhiSc560", "#phi_{Sc} (540-580MeV)", 24, -4, 4);
    PhiSc600 = new TH1D("PhiSc600", "#phi_{Sc} (580-620MeV)", 24, -4, 4);
    PhiSc640 = new TH1D("PhiSc640", "#phi_{Sc} (620-660MeV)", 24, -4, 4);
    PhiSc680 = new TH1D("PhiSc680", "#phi_{Sc} (660-700MeV)", 24, -4, 4);
    PhiSc320Random = new TH1D("PhiSc320Random", "#phi_{Sc} (300-340MeV)", 24, -4, 4);
    PhiSc360Random = new TH1D("PhiSc360Random", "#phi_{Sc} (340-380MeV)", 24, -4, 4);
    PhiSc400Random = new TH1D("PhiSc400Random", "#phi_{Sc} (380-420MeV)", 24, -4, 4);
    PhiSc440Random = new TH1D("PhiSc440Random", "#phi_{Sc} (420-460MeV)", 24, -4, 4);
    PhiSc480Random = new TH1D("PhiSc480Random", "#phi_{Sc} (460-500MeV)", 24, -4, 4);
    PhiSc520Random = new TH1D("PhiSc520Random", "#phi_{Sc} (500-540MeV)", 24, -4, 4);
    PhiSc560Random = new TH1D("PhiSc560Random", "#phi_{Sc} (540-580MeV)", 24, -4, 4);
    PhiSc600Random = new TH1D("PhiSc600Random", "#phi_{Sc} (580-620MeV)", 24, -4, 4);
    PhiSc640Random = new TH1D("PhiSc640Random", "#phi_{Sc} (620-660MeV)", 24, -4, 4);
    PhiSc680Random = new TH1D("PhiSc680Random", "#phi_{Sc} (660-700MeV)", 24, -4, 4);
    PhiSc320Prompt = new TH1D("PhiSc320Prompt", "#phi_{Sc} (300-340MeV)", 24, -4, 4);
    PhiSc360Prompt = new TH1D("PhiSc360Prompt", "#phi_{Sc} (340-380MeV)", 24, -4, 4);
    PhiSc400Prompt = new TH1D("PhiSc400Prompt", "#phi_{Sc} (380-420MeV)", 24, -4, 4);
    PhiSc440Prompt = new TH1D("PhiSc440Prompt", "#phi_{Sc} (420-460MeV)", 24, -4, 4);
    PhiSc480Prompt = new TH1D("PhiSc480Prompt", "#phi_{Sc} (460-500MeV)", 24, -4, 4);
    PhiSc520Prompt = new TH1D("PhiSc520Prompt", "#phi_{Sc} (500-540MeV)", 24, -4, 4);
    PhiSc560Prompt = new TH1D("PhiSc560Prompt", "#phi_{Sc} (540-580MeV)", 24, -4, 4);
    PhiSc600Prompt = new TH1D("PhiSc600Prompt", "#phi_{Sc} (580-620MeV)", 24, -4, 4);
    PhiSc640Prompt = new TH1D("PhiSc640Prompt", "#phi_{Sc} (620-660MeV)", 24, -4, 4);
    PhiSc680Prompt = new TH1D("PhiSc680Prompt", "#phi_{Sc} (660-700MeV)", 24, -4, 4);

    //    // Angles of neutron in scattered frame across EGamma bins for negative helicity
    PhiSc265NegHelCM1 = new TH1D( "Phi_Scattered_265MeV_NegHelCM1", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 20, -4, 4);
    PhiSc335NegHelCM1 = new TH1D( "Phi_Scattered_335MeV_NegHelCM1", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 20, -4, 4);
    PhiSc405NegHelCM1 = new TH1D( "Phi_Scattered_405MeV_NegHelCM1", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 20, -4, 4);
    PhiSc475NegHelCM1 = new TH1D( "Phi_Scattered_475MeV_NegHelCM1", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 20, -4, 4);
    PhiSc545NegHelCM1 = new TH1D( "Phi_Scattered_545MeV_NegHelCM1", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 20, -4, 4);
    PhiSc615NegHelCM1 = new TH1D( "Phi_Scattered_615MeV_NegHelCM1", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 20, -4, 4);
    PhiSc685NegHelCM1 = new TH1D( "Phi_Scattered_685MeV_NegHelCM1", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 20, -4, 4);

    PhiSc265NegHelCM2 = new TH1D( "Phi_Scattered_265MeV_NegHelCM2", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 20, -4, 4);
    PhiSc335NegHelCM2 = new TH1D( "Phi_Scattered_335MeV_NegHelCM2", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 20, -4, 4);
    PhiSc405NegHelCM2 = new TH1D( "Phi_Scattered_405MeV_NegHelCM2", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 20, -4, 4);
    PhiSc475NegHelCM2 = new TH1D( "Phi_Scattered_475MeV_NegHelCM2", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 20, -4, 4);
    PhiSc545NegHelCM2 = new TH1D( "Phi_Scattered_545MeV_NegHelCM2", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 20, -4, 4);
    PhiSc615NegHelCM2 = new TH1D( "Phi_Scattered_615MeV_NegHelCM2", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 20, -4, 4);
    PhiSc685NegHelCM2 = new TH1D( "Phi_Scattered_685MeV_NegHelCM2", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 20, -4, 4);

    PhiSc265NegHelCM3 = new TH1D( "Phi_Scattered_265MeV_NegHelCM3", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 20, -4, 4);
    PhiSc335NegHelCM3 = new TH1D( "Phi_Scattered_335MeV_NegHelCM3", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 20, -4, 4);
    PhiSc405NegHelCM3 = new TH1D( "Phi_Scattered_405MeV_NegHelCM3", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 20, -4, 4);
    PhiSc475NegHelCM3 = new TH1D( "Phi_Scattered_475MeV_NegHelCM3", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 20, -4, 4);
    PhiSc545NegHelCM3 = new TH1D( "Phi_Scattered_545MeV_NegHelCM3", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 20, -4, 4);
    PhiSc615NegHelCM3 = new TH1D( "Phi_Scattered_615MeV_NegHelCM3", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 20, -4, 4);
    PhiSc685NegHelCM3 = new TH1D( "Phi_Scattered_685MeV_NegHelCM3", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 20, -4, 4);

    PhiSc265NegHelCM4 = new TH1D( "Phi_Scattered_265MeV_NegHelCM4", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 20, -4, 4);
    PhiSc335NegHelCM4 = new TH1D( "Phi_Scattered_335MeV_NegHelCM4", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 20, -4, 4);
    PhiSc405NegHelCM4 = new TH1D( "Phi_Scattered_405MeV_NegHelCM4", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 20, -4, 4);
    PhiSc475NegHelCM4 = new TH1D( "Phi_Scattered_475MeV_NegHelCM4", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 20, -4, 4);
    PhiSc545NegHelCM4 = new TH1D( "Phi_Scattered_545MeV_NegHelCM4", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 20, -4, 4);
    PhiSc615NegHelCM4 = new TH1D( "Phi_Scattered_615MeV_NegHelCM4", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 20, -4, 4);
    PhiSc685NegHelCM4 = new TH1D( "Phi_Scattered_685MeV_NegHelCM4", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 20, -4, 4);

    PhiSc265NegHelCM5 = new TH1D( "Phi_Scattered_265MeV_NegHelCM5", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for -ve Helicity", 20, -4, 4);
    PhiSc335NegHelCM5 = new TH1D( "Phi_Scattered_335MeV_NegHelCM5", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for -ve Helicity", 20, -4, 4);
    PhiSc405NegHelCM5 = new TH1D( "Phi_Scattered_405MeV_NegHelCM5", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for -ve Helicity", 20, -4, 4);
    PhiSc475NegHelCM5 = new TH1D( "Phi_Scattered_475MeV_NegHelCM5", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for -ve Helicity", 20, -4, 4);
    PhiSc545NegHelCM5 = new TH1D( "Phi_Scattered_545MeV_NegHelCM5", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for -ve Helicity", 20, -4, 4);
    PhiSc615NegHelCM5 = new TH1D( "Phi_Scattered_615MeV_NegHelCM5", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for -ve Helicity", 20, -4, 4);
    PhiSc685NegHelCM5 = new TH1D( "Phi_Scattered_685MeV_NegHelCM5", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for -ve Helicity", 20, -4, 4);

    PhiSc265NegHelCM6 = new TH1D( "Phi_Scattered_265MeV_NegHelCM6", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 20, -4, 4);
    PhiSc335NegHelCM6 = new TH1D( "Phi_Scattered_335MeV_NegHelCM6", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 20, -4, 4);
    PhiSc405NegHelCM6 = new TH1D( "Phi_Scattered_405MeV_NegHelCM6", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 20, -4, 4);
    PhiSc475NegHelCM6 = new TH1D( "Phi_Scattered_475MeV_NegHelCM6", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 20, -4, 4);
    PhiSc545NegHelCM6 = new TH1D( "Phi_Scattered_545MeV_NegHelCM6", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 20, -4, 4);
    PhiSc615NegHelCM6 = new TH1D( "Phi_Scattered_615MeV_NegHelCM6", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 20, -4, 4);
    PhiSc685NegHelCM6 = new TH1D( "Phi_Scattered_685MeV_NegHelCM6", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 20, -4, 4);

    PhiSc265NegHelCM7 = new TH1D( "Phi_Scattered_265MeV_NegHelCM7", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 20, -4, 4);
    PhiSc335NegHelCM7 = new TH1D( "Phi_Scattered_335MeV_NegHelCM7", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 20, -4, 4);
    PhiSc405NegHelCM7 = new TH1D( "Phi_Scattered_405MeV_NegHelCM7", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 20, -4, 4);
    PhiSc475NegHelCM7 = new TH1D( "Phi_Scattered_475MeV_NegHelCM7", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 20, -4, 4);
    PhiSc545NegHelCM7 = new TH1D( "Phi_Scattered_545MeV_NegHelCM7", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 20, -4, 4);
    PhiSc615NegHelCM7 = new TH1D( "Phi_Scattered_615MeV_NegHelCM7", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 20, -4, 4);
    PhiSc685NegHelCM7 = new TH1D( "Phi_Scattered_685MeV_NegHelCM7", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 20, -4, 4);

    PhiSc265NegHelCM8 = new TH1D( "Phi_Scattered_265MeV_NegHelCM8", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 20, -4, 4);
    PhiSc335NegHelCM8 = new TH1D( "Phi_Scattered_335MeV_NegHelCM8", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 20, -4, 4);
    PhiSc405NegHelCM8 = new TH1D( "Phi_Scattered_405MeV_NegHelCM8", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 20, -4, 4);
    PhiSc475NegHelCM8 = new TH1D( "Phi_Scattered_475MeV_NegHelCM8", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 20, -4, 4);
    PhiSc545NegHelCM8 = new TH1D( "Phi_Scattered_545MeV_NegHelCM8", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 20, -4, 4);
    PhiSc615NegHelCM8 = new TH1D( "Phi_Scattered_615MeV_NegHelCM8", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 20, -4, 4);
    PhiSc685NegHelCM8 = new TH1D( "Phi_Scattered_685MeV_NegHelCM8", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 20, -4, 4);

    // Angles of neutron in scattered frame across EGamma bins for positive helicity
    PhiSc265PosHelCM1 = new TH1D( "Phi_Scattered_265MeV_PosHelCM1", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 20, -4, 4);
    PhiSc335PosHelCM1 = new TH1D( "Phi_Scattered_335MeV_PosHelCM1", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 20, -4, 4);
    PhiSc405PosHelCM1 = new TH1D( "Phi_Scattered_405MeV_PosHelCM1", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 20, -4, 4);
    PhiSc475PosHelCM1 = new TH1D( "Phi_Scattered_475MeV_PosHelCM1", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 20, -4, 4);
    PhiSc545PosHelCM1 = new TH1D( "Phi_Scattered_545MeV_PosHelCM1", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 20, -4, 4);
    PhiSc615PosHelCM1 = new TH1D( "Phi_Scattered_615MeV_PosHelCM1", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 20, -4, 4);
    PhiSc685PosHelCM1 = new TH1D( "Phi_Scattered_685MeV_PosHelCM1", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 20, -4, 4);

    PhiSc265PosHelCM2 = new TH1D( "Phi_Scattered_265MeV_PosHelCM2", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 20, -4, 4);
    PhiSc335PosHelCM2 = new TH1D( "Phi_Scattered_335MeV_PosHelCM2", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 20, -4, 4);
    PhiSc405PosHelCM2 = new TH1D( "Phi_Scattered_405MeV_PosHelCM2", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 20, -4, 4);
    PhiSc475PosHelCM2 = new TH1D( "Phi_Scattered_475MeV_PosHelCM2", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 20, -4, 4);
    PhiSc545PosHelCM2 = new TH1D( "Phi_Scattered_545MeV_PosHelCM2", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 20, -4, 4);
    PhiSc615PosHelCM2 = new TH1D( "Phi_Scattered_615MeV_PosHelCM2", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 20, -4, 4);
    PhiSc685PosHelCM2 = new TH1D( "Phi_Scattered_685MeV_PosHelCM2", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 20, -4, 4);

    PhiSc265PosHelCM3 = new TH1D( "Phi_Scattered_265MeV_PosHelCM3", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 20, -4, 4);
    PhiSc335PosHelCM3 = new TH1D( "Phi_Scattered_335MeV_PosHelCM3", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 20, -4, 4);
    PhiSc405PosHelCM3 = new TH1D( "Phi_Scattered_405MeV_PosHelCM3", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 20, -4, 4);
    PhiSc475PosHelCM3 = new TH1D( "Phi_Scattered_475MeV_PosHelCM3", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 20, -4, 4);
    PhiSc545PosHelCM3 = new TH1D( "Phi_Scattered_545MeV_PosHelCM3", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 20, -4, 4);
    PhiSc615PosHelCM3 = new TH1D( "Phi_Scattered_615MeV_PosHelCM3", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 20, -4, 4);
    PhiSc685PosHelCM3 = new TH1D( "Phi_Scattered_685MeV_PosHelCM3", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 20, -4, 4);

    PhiSc265PosHelCM4 = new TH1D( "Phi_Scattered_265MeV_PosHelCM4", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 20, -4, 4);
    PhiSc335PosHelCM4 = new TH1D( "Phi_Scattered_335MeV_PosHelCM4", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 20, -4, 4);
    PhiSc405PosHelCM4 = new TH1D( "Phi_Scattered_405MeV_PosHelCM4", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 20, -4, 4);
    PhiSc475PosHelCM4 = new TH1D( "Phi_Scattered_475MeV_PosHelCM4", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 20, -4, 4);
    PhiSc545PosHelCM4 = new TH1D( "Phi_Scattered_545MeV_PosHelCM4", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 20, -4, 4);
    PhiSc615PosHelCM4 = new TH1D( "Phi_Scattered_615MeV_PosHelCM4", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 20, -4, 4);
    PhiSc685PosHelCM4 = new TH1D( "Phi_Scattered_685MeV_PosHelCM4", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 20, -4, 4);

    PhiSc265PosHelCM5 = new TH1D( "Phi_Scattered_265MeV_PosHelCM5", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for +ve Helicity", 20, -4, 4);
    PhiSc335PosHelCM5 = new TH1D( "Phi_Scattered_335MeV_PosHelCM5", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for +ve Helicity", 20, -4, 4);
    PhiSc405PosHelCM5 = new TH1D( "Phi_Scattered_405MeV_PosHelCM5", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for +ve Helicity", 20, -4, 4);
    PhiSc475PosHelCM5 = new TH1D( "Phi_Scattered_475MeV_PosHelCM5", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for +ve Helicity", 20, -4, 4);
    PhiSc545PosHelCM5 = new TH1D( "Phi_Scattered_545MeV_PosHelCM5", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for +ve Helicity", 20, -4, 4);
    PhiSc615PosHelCM5 = new TH1D( "Phi_Scattered_615MeV_PosHelCM5", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for +ve Helicity", 20, -4, 4);
    PhiSc685PosHelCM5 = new TH1D( "Phi_Scattered_685MeV_PosHelCM5", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for +ve Helicity", 20, -4, 4);

    PhiSc265PosHelCM6 = new TH1D( "Phi_Scattered_265MeV_PosHelCM6", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 20, -4, 4);
    PhiSc335PosHelCM6 = new TH1D( "Phi_Scattered_335MeV_PosHelCM6", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 20, -4, 4);
    PhiSc405PosHelCM6 = new TH1D( "Phi_Scattered_405MeV_PosHelCM6", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 20, -4, 4);
    PhiSc475PosHelCM6 = new TH1D( "Phi_Scattered_475MeV_PosHelCM6", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 20, -4, 4);
    PhiSc545PosHelCM6 = new TH1D( "Phi_Scattered_545MeV_PosHelCM6", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 20, -4, 4);
    PhiSc615PosHelCM6 = new TH1D( "Phi_Scattered_615MeV_PosHelCM6", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 20, -4, 4);
    PhiSc685PosHelCM6 = new TH1D( "Phi_Scattered_685MeV_PosHelCM6", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 20, -4, 4);

    PhiSc265PosHelCM7 = new TH1D( "Phi_Scattered_265MeV_PosHelCM7", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 20, -4, 4);
    PhiSc335PosHelCM7 = new TH1D( "Phi_Scattered_335MeV_PosHelCM7", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 20, -4, 4);
    PhiSc405PosHelCM7 = new TH1D( "Phi_Scattered_405MeV_PosHelCM7", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 20, -4, 4);
    PhiSc475PosHelCM7 = new TH1D( "Phi_Scattered_475MeV_PosHelCM7", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 20, -4, 4);
    PhiSc545PosHelCM7 = new TH1D( "Phi_Scattered_545MeV_PosHelCM7", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 20, -4, 4);
    PhiSc615PosHelCM7 = new TH1D( "Phi_Scattered_615MeV_PosHelCM7", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 20, -4, 4);
    PhiSc685PosHelCM7 = new TH1D( "Phi_Scattered_685MeV_PosHelCM7", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 20, -4, 4);

    PhiSc265PosHelCM8 = new TH1D( "Phi_Scattered_265MeV_PosHelCM8", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 20, -4, 4);
    PhiSc335PosHelCM8 = new TH1D( "Phi_Scattered_335MeV_PosHelCM8", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 20, -4, 4);
    PhiSc405PosHelCM8 = new TH1D( "Phi_Scattered_405MeV_PosHelCM8", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 20, -4, 4);
    PhiSc475PosHelCM8 = new TH1D( "Phi_Scattered_475MeV_PosHelCM8", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 20, -4, 4);
    PhiSc545PosHelCM8 = new TH1D( "Phi_Scattered_545MeV_PosHelCM8", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 20, -4, 4);
    PhiSc615PosHelCM8 = new TH1D( "Phi_Scattered_615MeV_PosHelCM8", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 20, -4, 4);
    PhiSc685PosHelCM8 = new TH1D( "Phi_Scattered_685MeV_PosHelCM8", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 20, -4, 4);

    PhiSc265NegHelCM1Prompt = new TH1D( "Phi_Scattered_265MeV_NegHelCM1Prompt", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 20, -4, 4);
    PhiSc335NegHelCM1Prompt = new TH1D( "Phi_Scattered_335MeV_NegHelCM1Prompt", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 20, -4, 4);
    PhiSc405NegHelCM1Prompt = new TH1D( "Phi_Scattered_405MeV_NegHelCM1Prompt", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 20, -4, 4);
    PhiSc475NegHelCM1Prompt = new TH1D( "Phi_Scattered_475MeV_NegHelCM1Prompt", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 20, -4, 4);
    PhiSc545NegHelCM1Prompt = new TH1D( "Phi_Scattered_545MeV_NegHelCM1Prompt", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 20, -4, 4);
    PhiSc615NegHelCM1Prompt = new TH1D( "Phi_Scattered_615MeV_NegHelCM1Prompt", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 20, -4, 4);
    PhiSc685NegHelCM1Prompt = new TH1D( "Phi_Scattered_685MeV_NegHelCM1Prompt", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 20, -4, 4);

    PhiSc265NegHelCM2Prompt = new TH1D( "Phi_Scattered_265MeV_NegHelCM2Prompt", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 20, -4, 4);
    PhiSc335NegHelCM2Prompt = new TH1D( "Phi_Scattered_335MeV_NegHelCM2Prompt", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 20, -4, 4);
    PhiSc405NegHelCM2Prompt = new TH1D( "Phi_Scattered_405MeV_NegHelCM2Prompt", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 20, -4, 4);
    PhiSc475NegHelCM2Prompt = new TH1D( "Phi_Scattered_475MeV_NegHelCM2Prompt", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 20, -4, 4);
    PhiSc545NegHelCM2Prompt = new TH1D( "Phi_Scattered_545MeV_NegHelCM2Prompt", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 20, -4, 4);
    PhiSc615NegHelCM2Prompt = new TH1D( "Phi_Scattered_615MeV_NegHelCM2Prompt", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 20, -4, 4);
    PhiSc685NegHelCM2Prompt = new TH1D( "Phi_Scattered_685MeV_NegHelCM2Prompt", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 20, -4, 4);

    PhiSc265NegHelCM3Prompt = new TH1D( "Phi_Scattered_265MeV_NegHelCM3Prompt", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 20, -4, 4);
    PhiSc335NegHelCM3Prompt = new TH1D( "Phi_Scattered_335MeV_NegHelCM3Prompt", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 20, -4, 4);
    PhiSc405NegHelCM3Prompt = new TH1D( "Phi_Scattered_405MeV_NegHelCM3Prompt", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 20, -4, 4);
    PhiSc475NegHelCM3Prompt = new TH1D( "Phi_Scattered_475MeV_NegHelCM3Prompt", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 20, -4, 4);
    PhiSc545NegHelCM3Prompt = new TH1D( "Phi_Scattered_545MeV_NegHelCM3Prompt", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 20, -4, 4);
    PhiSc615NegHelCM3Prompt = new TH1D( "Phi_Scattered_615MeV_NegHelCM3Prompt", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 20, -4, 4);
    PhiSc685NegHelCM3Prompt = new TH1D( "Phi_Scattered_685MeV_NegHelCM3Prompt", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 20, -4, 4);

    PhiSc265NegHelCM4Prompt = new TH1D( "Phi_Scattered_265MeV_NegHelCM4Prompt", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 20, -4, 4);
    PhiSc335NegHelCM4Prompt = new TH1D( "Phi_Scattered_335MeV_NegHelCM4Prompt", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 20, -4, 4);
    PhiSc405NegHelCM4Prompt = new TH1D( "Phi_Scattered_405MeV_NegHelCM4Prompt", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 20, -4, 4);
    PhiSc475NegHelCM4Prompt = new TH1D( "Phi_Scattered_475MeV_NegHelCM4Prompt", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 20, -4, 4);
    PhiSc545NegHelCM4Prompt = new TH1D( "Phi_Scattered_545MeV_NegHelCM4Prompt", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 20, -4, 4);
    PhiSc615NegHelCM4Prompt = new TH1D( "Phi_Scattered_615MeV_NegHelCM4Prompt", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 20, -4, 4);
    PhiSc685NegHelCM4Prompt = new TH1D( "Phi_Scattered_685MeV_NegHelCM4Prompt", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 20, -4, 4);

    PhiSc265NegHelCM5Prompt = new TH1D( "Phi_Scattered_265MeV_NegHelCM5Prompt", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for -ve Helicity", 20, -4, 4);
    PhiSc335NegHelCM5Prompt = new TH1D( "Phi_Scattered_335MeV_NegHelCM5Prompt", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for -ve Helicity", 20, -4, 4);
    PhiSc405NegHelCM5Prompt = new TH1D( "Phi_Scattered_405MeV_NegHelCM5Prompt", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for -ve Helicity", 20, -4, 4);
    PhiSc475NegHelCM5Prompt = new TH1D( "Phi_Scattered_475MeV_NegHelCM5Prompt", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for -ve Helicity", 20, -4, 4);
    PhiSc545NegHelCM5Prompt = new TH1D( "Phi_Scattered_545MeV_NegHelCM5Prompt", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for -ve Helicity", 20, -4, 4);
    PhiSc615NegHelCM5Prompt = new TH1D( "Phi_Scattered_615MeV_NegHelCM5Prompt", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for -ve Helicity", 20, -4, 4);
    PhiSc685NegHelCM5Prompt = new TH1D( "Phi_Scattered_685MeV_NegHelCM5Prompt", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for -ve Helicity", 20, -4, 4);

    PhiSc265NegHelCM6Prompt = new TH1D( "Phi_Scattered_265MeV_NegHelCM6Prompt", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 20, -4, 4);
    PhiSc335NegHelCM6Prompt = new TH1D( "Phi_Scattered_335MeV_NegHelCM6Prompt", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 20, -4, 4);
    PhiSc405NegHelCM6Prompt = new TH1D( "Phi_Scattered_405MeV_NegHelCM6Prompt", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 20, -4, 4);
    PhiSc475NegHelCM6Prompt = new TH1D( "Phi_Scattered_475MeV_NegHelCM6Prompt", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 20, -4, 4);
    PhiSc545NegHelCM6Prompt = new TH1D( "Phi_Scattered_545MeV_NegHelCM6Prompt", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 20, -4, 4);
    PhiSc615NegHelCM6Prompt = new TH1D( "Phi_Scattered_615MeV_NegHelCM6Prompt", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 20, -4, 4);
    PhiSc685NegHelCM6Prompt = new TH1D( "Phi_Scattered_685MeV_NegHelCM6Prompt", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 20, -4, 4);

    PhiSc265NegHelCM7Prompt = new TH1D( "Phi_Scattered_265MeV_NegHelCM7Prompt", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 20, -4, 4);
    PhiSc335NegHelCM7Prompt = new TH1D( "Phi_Scattered_335MeV_NegHelCM7Prompt", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 20, -4, 4);
    PhiSc405NegHelCM7Prompt = new TH1D( "Phi_Scattered_405MeV_NegHelCM7Prompt", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 20, -4, 4);
    PhiSc475NegHelCM7Prompt = new TH1D( "Phi_Scattered_475MeV_NegHelCM7Prompt", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 20, -4, 4);
    PhiSc545NegHelCM7Prompt = new TH1D( "Phi_Scattered_545MeV_NegHelCM7Prompt", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 20, -4, 4);
    PhiSc615NegHelCM7Prompt = new TH1D( "Phi_Scattered_615MeV_NegHelCM7Prompt", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 20, -4, 4);
    PhiSc685NegHelCM7Prompt = new TH1D( "Phi_Scattered_685MeV_NegHelCM7Prompt", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 20, -4, 4);

    PhiSc265NegHelCM8Prompt = new TH1D( "Phi_Scattered_265MeV_NegHelCM8Prompt", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 20, -4, 4);
    PhiSc335NegHelCM8Prompt = new TH1D( "Phi_Scattered_335MeV_NegHelCM8Prompt", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 20, -4, 4);
    PhiSc405NegHelCM8Prompt = new TH1D( "Phi_Scattered_405MeV_NegHelCM8Prompt", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 20, -4, 4);
    PhiSc475NegHelCM8Prompt = new TH1D( "Phi_Scattered_475MeV_NegHelCM8Prompt", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 20, -4, 4);
    PhiSc545NegHelCM8Prompt = new TH1D( "Phi_Scattered_545MeV_NegHelCM8Prompt", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 20, -4, 4);
    PhiSc615NegHelCM8Prompt = new TH1D( "Phi_Scattered_615MeV_NegHelCM8Prompt", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 20, -4, 4);
    PhiSc685NegHelCM8Prompt = new TH1D( "Phi_Scattered_685MeV_NegHelCM8Prompt", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 20, -4, 4);

    // Angles of neutron in scattered frame across EGamma bins for positive helicity
    PhiSc265PosHelCM1Prompt = new TH1D( "Phi_Scattered_265MeV_PosHelCM1Prompt", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 20, -4, 4);
    PhiSc335PosHelCM1Prompt = new TH1D( "Phi_Scattered_335MeV_PosHelCM1Prompt", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 20, -4, 4);
    PhiSc405PosHelCM1Prompt = new TH1D( "Phi_Scattered_405MeV_PosHelCM1Prompt", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 20, -4, 4);
    PhiSc475PosHelCM1Prompt = new TH1D( "Phi_Scattered_475MeV_PosHelCM1Prompt", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 20, -4, 4);
    PhiSc545PosHelCM1Prompt = new TH1D( "Phi_Scattered_545MeV_PosHelCM1Prompt", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 20, -4, 4);
    PhiSc615PosHelCM1Prompt = new TH1D( "Phi_Scattered_615MeV_PosHelCM1Prompt", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 20, -4, 4);
    PhiSc685PosHelCM1Prompt = new TH1D( "Phi_Scattered_685MeV_PosHelCM1Prompt", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 20, -4, 4);

    PhiSc265PosHelCM2Prompt = new TH1D( "Phi_Scattered_265MeV_PosHelCM2Prompt", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 20, -4, 4);
    PhiSc335PosHelCM2Prompt = new TH1D( "Phi_Scattered_335MeV_PosHelCM2Prompt", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 20, -4, 4);
    PhiSc405PosHelCM2Prompt = new TH1D( "Phi_Scattered_405MeV_PosHelCM2Prompt", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 20, -4, 4);
    PhiSc475PosHelCM2Prompt = new TH1D( "Phi_Scattered_475MeV_PosHelCM2Prompt", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 20, -4, 4);
    PhiSc545PosHelCM2Prompt = new TH1D( "Phi_Scattered_545MeV_PosHelCM2Prompt", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 20, -4, 4);
    PhiSc615PosHelCM2Prompt = new TH1D( "Phi_Scattered_615MeV_PosHelCM2Prompt", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 20, -4, 4);
    PhiSc685PosHelCM2Prompt = new TH1D( "Phi_Scattered_685MeV_PosHelCM2Prompt", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 20, -4, 4);

    PhiSc265PosHelCM3Prompt = new TH1D( "Phi_Scattered_265MeV_PosHelCM3Prompt", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 20, -4, 4);
    PhiSc335PosHelCM3Prompt = new TH1D( "Phi_Scattered_335MeV_PosHelCM3Prompt", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 20, -4, 4);
    PhiSc405PosHelCM3Prompt = new TH1D( "Phi_Scattered_405MeV_PosHelCM3Prompt", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 20, -4, 4);
    PhiSc475PosHelCM3Prompt = new TH1D( "Phi_Scattered_475MeV_PosHelCM3Prompt", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 20, -4, 4);
    PhiSc545PosHelCM3Prompt = new TH1D( "Phi_Scattered_545MeV_PosHelCM3Prompt", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 20, -4, 4);
    PhiSc615PosHelCM3Prompt = new TH1D( "Phi_Scattered_615MeV_PosHelCM3Prompt", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 20, -4, 4);
    PhiSc685PosHelCM3Prompt = new TH1D( "Phi_Scattered_685MeV_PosHelCM3Prompt", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 20, -4, 4);

    PhiSc265PosHelCM4Prompt = new TH1D( "Phi_Scattered_265MeV_PosHelCM4Prompt", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 20, -4, 4);
    PhiSc335PosHelCM4Prompt = new TH1D( "Phi_Scattered_335MeV_PosHelCM4Prompt", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 20, -4, 4);
    PhiSc405PosHelCM4Prompt = new TH1D( "Phi_Scattered_405MeV_PosHelCM4Prompt", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 20, -4, 4);
    PhiSc475PosHelCM4Prompt = new TH1D( "Phi_Scattered_475MeV_PosHelCM4Prompt", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 20, -4, 4);
    PhiSc545PosHelCM4Prompt = new TH1D( "Phi_Scattered_545MeV_PosHelCM4Prompt", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 20, -4, 4);
    PhiSc615PosHelCM4Prompt = new TH1D( "Phi_Scattered_615MeV_PosHelCM4Prompt", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 20, -4, 4);
    PhiSc685PosHelCM4Prompt = new TH1D( "Phi_Scattered_685MeV_PosHelCM4Prompt", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 20, -4, 4);

    PhiSc265PosHelCM5Prompt = new TH1D( "Phi_Scattered_265MeV_PosHelCM5Prompt", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for +ve Helicity", 20, -4, 4);
    PhiSc335PosHelCM5Prompt = new TH1D( "Phi_Scattered_335MeV_PosHelCM5Prompt", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for +ve Helicity", 20, -4, 4);
    PhiSc405PosHelCM5Prompt = new TH1D( "Phi_Scattered_405MeV_PosHelCM5Prompt", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for +ve Helicity", 20, -4, 4);
    PhiSc475PosHelCM5Prompt = new TH1D( "Phi_Scattered_475MeV_PosHelCM5Prompt", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for +ve Helicity", 20, -4, 4);
    PhiSc545PosHelCM5Prompt = new TH1D( "Phi_Scattered_545MeV_PosHelCM5Prompt", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for +ve Helicity", 20, -4, 4);
    PhiSc615PosHelCM5Prompt = new TH1D( "Phi_Scattered_615MeV_PosHelCM5Prompt", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for +ve Helicity", 20, -4, 4);
    PhiSc685PosHelCM5Prompt = new TH1D( "Phi_Scattered_685MeV_PosHelCM5Prompt", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for +ve Helicity", 20, -4, 4);

    PhiSc265PosHelCM6Prompt = new TH1D( "Phi_Scattered_265MeV_PosHelCM6Prompt", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 20, -4, 4);
    PhiSc335PosHelCM6Prompt = new TH1D( "Phi_Scattered_335MeV_PosHelCM6Prompt", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 20, -4, 4);
    PhiSc405PosHelCM6Prompt = new TH1D( "Phi_Scattered_405MeV_PosHelCM6Prompt", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 20, -4, 4);
    PhiSc475PosHelCM6Prompt = new TH1D( "Phi_Scattered_475MeV_PosHelCM6Prompt", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 20, -4, 4);
    PhiSc545PosHelCM6Prompt = new TH1D( "Phi_Scattered_545MeV_PosHelCM6Prompt", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 20, -4, 4);
    PhiSc615PosHelCM6Prompt = new TH1D( "Phi_Scattered_615MeV_PosHelCM6Prompt", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 20, -4, 4);
    PhiSc685PosHelCM6Prompt = new TH1D( "Phi_Scattered_685MeV_PosHelCM6Prompt", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 20, -4, 4);

    PhiSc265PosHelCM7Prompt = new TH1D( "Phi_Scattered_265MeV_PosHelCM7Prompt", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 20, -4, 4);
    PhiSc335PosHelCM7Prompt = new TH1D( "Phi_Scattered_335MeV_PosHelCM7Prompt", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 20, -4, 4);
    PhiSc405PosHelCM7Prompt = new TH1D( "Phi_Scattered_405MeV_PosHelCM7Prompt", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 20, -4, 4);
    PhiSc475PosHelCM7Prompt = new TH1D( "Phi_Scattered_475MeV_PosHelCM7Prompt", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 20, -4, 4);
    PhiSc545PosHelCM7Prompt = new TH1D( "Phi_Scattered_545MeV_PosHelCM7Prompt", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 20, -4, 4);
    PhiSc615PosHelCM7Prompt = new TH1D( "Phi_Scattered_615MeV_PosHelCM7Prompt", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 20, -4, 4);
    PhiSc685PosHelCM7Prompt = new TH1D( "Phi_Scattered_685MeV_PosHelCM7Prompt", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 20, -4, 4);

    PhiSc265PosHelCM8Prompt = new TH1D( "Phi_Scattered_265MeV_PosHelCM8Prompt", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 20, -4, 4);
    PhiSc335PosHelCM8Prompt = new TH1D( "Phi_Scattered_335MeV_PosHelCM8Prompt", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 20, -4, 4);
    PhiSc405PosHelCM8Prompt = new TH1D( "Phi_Scattered_405MeV_PosHelCM8Prompt", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 20, -4, 4);
    PhiSc475PosHelCM8Prompt = new TH1D( "Phi_Scattered_475MeV_PosHelCM8Prompt", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 20, -4, 4);
    PhiSc545PosHelCM8Prompt = new TH1D( "Phi_Scattered_545MeV_PosHelCM8Prompt", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 20, -4, 4);
    PhiSc615PosHelCM8Prompt = new TH1D( "Phi_Scattered_615MeV_PosHelCM8Prompt", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 20, -4, 4);
    PhiSc685PosHelCM8Prompt = new TH1D( "Phi_Scattered_685MeV_PosHelCM8Prompt", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 20, -4, 4);

    PhiSc265NegHelCM1Random = new TH1D( "Phi_Scattered_265MeV_NegHelCM1Random", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 20, -4, 4);
    PhiSc335NegHelCM1Random = new TH1D( "Phi_Scattered_335MeV_NegHelCM1Random", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 20, -4, 4);
    PhiSc405NegHelCM1Random = new TH1D( "Phi_Scattered_405MeV_NegHelCM1Random", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 20, -4, 4);
    PhiSc475NegHelCM1Random = new TH1D( "Phi_Scattered_475MeV_NegHelCM1Random", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 20, -4, 4);
    PhiSc545NegHelCM1Random = new TH1D( "Phi_Scattered_545MeV_NegHelCM1Random", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 20, -4, 4);
    PhiSc615NegHelCM1Random = new TH1D( "Phi_Scattered_615MeV_NegHelCM1Random", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 20, -4, 4);
    PhiSc685NegHelCM1Random = new TH1D( "Phi_Scattered_685MeV_NegHelCM1Random", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}1-0.75)) for -ve Helicity", 20, -4, 4);

    PhiSc265NegHelCM2Random = new TH1D( "Phi_Scattered_265MeV_NegHelCM2Random", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 20, -4, 4);
    PhiSc335NegHelCM2Random = new TH1D( "Phi_Scattered_335MeV_NegHelCM2Random", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 20, -4, 4);
    PhiSc405NegHelCM2Random = new TH1D( "Phi_Scattered_405MeV_NegHelCM2Random", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 20, -4, 4);
    PhiSc475NegHelCM2Random = new TH1D( "Phi_Scattered_475MeV_NegHelCM2Random", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 20, -4, 4);
    PhiSc545NegHelCM2Random = new TH1D( "Phi_Scattered_545MeV_NegHelCM2Random", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 20, -4, 4);
    PhiSc615NegHelCM2Random = new TH1D( "Phi_Scattered_615MeV_NegHelCM2Random", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 20, -4, 4);
    PhiSc685NegHelCM2Random = new TH1D( "Phi_Scattered_685MeV_NegHelCM2Random", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for -ve Helicity", 20, -4, 4);

    PhiSc265NegHelCM3Random = new TH1D( "Phi_Scattered_265MeV_NegHelCM3Random", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 20, -4, 4);
    PhiSc335NegHelCM3Random = new TH1D( "Phi_Scattered_335MeV_NegHelCM3Random", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 20, -4, 4);
    PhiSc405NegHelCM3Random = new TH1D( "Phi_Scattered_405MeV_NegHelCM3Random", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 20, -4, 4);
    PhiSc475NegHelCM3Random = new TH1D( "Phi_Scattered_475MeV_NegHelCM3Random", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 20, -4, 4);
    PhiSc545NegHelCM3Random = new TH1D( "Phi_Scattered_545MeV_NegHelCM3Random", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 20, -4, 4);
    PhiSc615NegHelCM3Random = new TH1D( "Phi_Scattered_615MeV_NegHelCM3Random", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 20, -4, 4);
    PhiSc685NegHelCM3Random = new TH1D( "Phi_Scattered_685MeV_NegHelCM3Random", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for -ve Helicity", 20, -4, 4);

    PhiSc265NegHelCM4Random = new TH1D( "Phi_Scattered_265MeV_NegHelCM4Random", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 20, -4, 4);
    PhiSc335NegHelCM4Random = new TH1D( "Phi_Scattered_335MeV_NegHelCM4Random", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 20, -4, 4);
    PhiSc405NegHelCM4Random = new TH1D( "Phi_Scattered_405MeV_NegHelCM4Random", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 20, -4, 4);
    PhiSc475NegHelCM4Random = new TH1D( "Phi_Scattered_475MeV_NegHelCM4Random", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 20, -4, 4);
    PhiSc545NegHelCM4Random = new TH1D( "Phi_Scattered_545MeV_NegHelCM4Random", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 20, -4, 4);
    PhiSc615NegHelCM4Random = new TH1D( "Phi_Scattered_615MeV_NegHelCM4Random", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 20, -4, 4);
    PhiSc685NegHelCM4Random = new TH1D( "Phi_Scattered_685MeV_NegHelCM4Random", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for -ve Helicity", 20, -4, 4);

    PhiSc265NegHelCM5Random = new TH1D( "Phi_Scattered_265MeV_NegHelCM5Random", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for -ve Helicity", 20, -4, 4);
    PhiSc335NegHelCM5Random = new TH1D( "Phi_Scattered_335MeV_NegHelCM5Random", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for -ve Helicity", 20, -4, 4);
    PhiSc405NegHelCM5Random = new TH1D( "Phi_Scattered_405MeV_NegHelCM5Random", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for -ve Helicity", 20, -4, 4);
    PhiSc475NegHelCM5Random = new TH1D( "Phi_Scattered_475MeV_NegHelCM5Random", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for -ve Helicity", 20, -4, 4);
    PhiSc545NegHelCM5Random = new TH1D( "Phi_Scattered_545MeV_NegHelCM5Random", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for -ve Helicity", 20, -4, 4);
    PhiSc615NegHelCM5Random = new TH1D( "Phi_Scattered_615MeV_NegHelCM5Random", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for -ve Helicity", 20, -4, 4);
    PhiSc685NegHelCM5Random = new TH1D( "Phi_Scattered_685MeV_NegHelCM5Random", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for -ve Helicity", 20, -4, 4);

    PhiSc265NegHelCM6Random = new TH1D( "Phi_Scattered_265MeV_NegHelCM6Random", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 20, -4, 4);
    PhiSc335NegHelCM6Random = new TH1D( "Phi_Scattered_335MeV_NegHelCM6Random", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 20, -4, 4);
    PhiSc405NegHelCM6Random = new TH1D( "Phi_Scattered_405MeV_NegHelCM6Random", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 20, -4, 4);
    PhiSc475NegHelCM6Random = new TH1D( "Phi_Scattered_475MeV_NegHelCM6Random", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 20, -4, 4);
    PhiSc545NegHelCM6Random = new TH1D( "Phi_Scattered_545MeV_NegHelCM6Random", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 20, -4, 4);
    PhiSc615NegHelCM6Random = new TH1D( "Phi_Scattered_615MeV_NegHelCM6Random", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 20, -4, 4);
    PhiSc685NegHelCM6Random = new TH1D( "Phi_Scattered_685MeV_NegHelCM6Random", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for -ve Helicity", 20, -4, 4);

    PhiSc265NegHelCM7Random = new TH1D( "Phi_Scattered_265MeV_NegHelCM7Random", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 20, -4, 4);
    PhiSc335NegHelCM7Random = new TH1D( "Phi_Scattered_335MeV_NegHelCM7Random", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 20, -4, 4);
    PhiSc405NegHelCM7Random = new TH1D( "Phi_Scattered_405MeV_NegHelCM7Random", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 20, -4, 4);
    PhiSc475NegHelCM7Random = new TH1D( "Phi_Scattered_475MeV_NegHelCM7Random", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 20, -4, 4);
    PhiSc545NegHelCM7Random = new TH1D( "Phi_Scattered_545MeV_NegHelCM7Random", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 20, -4, 4);
    PhiSc615NegHelCM7Random = new TH1D( "Phi_Scattered_615MeV_NegHelCM7Random", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 20, -4, 4);
    PhiSc685NegHelCM7Random = new TH1D( "Phi_Scattered_685MeV_NegHelCM7Random", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for -ve Helicity", 20, -4, 4);

    PhiSc265NegHelCM8Random = new TH1D( "Phi_Scattered_265MeV_NegHelCM8Random", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 20, -4, 4);
    PhiSc335NegHelCM8Random = new TH1D( "Phi_Scattered_335MeV_NegHelCM8Random", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 20, -4, 4);
    PhiSc405NegHelCM8Random = new TH1D( "Phi_Scattered_405MeV_NegHelCM8Random", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 20, -4, 4);
    PhiSc475NegHelCM8Random = new TH1D( "Phi_Scattered_475MeV_NegHelCM8Random", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 20, -4, 4);
    PhiSc545NegHelCM8Random = new TH1D( "Phi_Scattered_545MeV_NegHelCM8Random", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 20, -4, 4);
    PhiSc615NegHelCM8Random = new TH1D( "Phi_Scattered_615MeV_NegHelCM8Random", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 20, -4, 4);
    PhiSc685NegHelCM8Random = new TH1D( "Phi_Scattered_685MeV_NegHelCM8Random", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for -ve Helicity", 20, -4, 4);

    // Angles of neutron in scattered frame across EGamma bins for positive helicity
    PhiSc265PosHelCM1Random = new TH1D( "Phi_Scattered_265MeV_PosHelCM1Random", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 20, -4, 4);
    PhiSc335PosHelCM1Random = new TH1D( "Phi_Scattered_335MeV_PosHelCM1Random", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 20, -4, 4);
    PhiSc405PosHelCM1Random = new TH1D( "Phi_Scattered_405MeV_PosHelCM1Random", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 20, -4, 4);
    PhiSc475PosHelCM1Random = new TH1D( "Phi_Scattered_475MeV_PosHelCM1Random", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 20, -4, 4);
    PhiSc545PosHelCM1Random = new TH1D( "Phi_Scattered_545MeV_PosHelCM1Random", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 20, -4, 4);
    PhiSc615PosHelCM1Random = new TH1D( "Phi_Scattered_615MeV_PosHelCM1Random", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 20, -4, 4);
    PhiSc685PosHelCM1Random = new TH1D( "Phi_Scattered_685MeV_PosHelCM1Random", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}1-0.75)) for +ve Helicity", 20, -4, 4);

    PhiSc265PosHelCM2Random = new TH1D( "Phi_Scattered_265MeV_PosHelCM2Random", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 20, -4, 4);
    PhiSc335PosHelCM2Random = new TH1D( "Phi_Scattered_335MeV_PosHelCM2Random", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 20, -4, 4);
    PhiSc405PosHelCM2Random = new TH1D( "Phi_Scattered_405MeV_PosHelCM2Random", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 20, -4, 4);
    PhiSc475PosHelCM2Random = new TH1D( "Phi_Scattered_475MeV_PosHelCM2Random", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 20, -4, 4);
    PhiSc545PosHelCM2Random = new TH1D( "Phi_Scattered_545MeV_PosHelCM2Random", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 20, -4, 4);
    PhiSc615PosHelCM2Random = new TH1D( "Phi_Scattered_615MeV_PosHelCM2Random", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 20, -4, 4);
    PhiSc685PosHelCM2Random = new TH1D( "Phi_Scattered_685MeV_PosHelCM2Random", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.75-0.5)) for +ve Helicity", 20, -4, 4);

    PhiSc265PosHelCM3Random = new TH1D( "Phi_Scattered_265MeV_PosHelCM3Random", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 20, -4, 4);
    PhiSc335PosHelCM3Random = new TH1D( "Phi_Scattered_335MeV_PosHelCM3Random", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 20, -4, 4);
    PhiSc405PosHelCM3Random = new TH1D( "Phi_Scattered_405MeV_PosHelCM3Random", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 20, -4, 4);
    PhiSc475PosHelCM3Random = new TH1D( "Phi_Scattered_475MeV_PosHelCM3Random", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 20, -4, 4);
    PhiSc545PosHelCM3Random = new TH1D( "Phi_Scattered_545MeV_PosHelCM3Random", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 20, -4, 4);
    PhiSc615PosHelCM3Random = new TH1D( "Phi_Scattered_615MeV_PosHelCM3Random", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 20, -4, 4);
    PhiSc685PosHelCM3Random = new TH1D( "Phi_Scattered_685MeV_PosHelCM3Random", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.5-0.25)) for +ve Helicity", 20, -4, 4);

    PhiSc265PosHelCM4Random = new TH1D( "Phi_Scattered_265MeV_PosHelCM4Random", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 20, -4, 4);
    PhiSc335PosHelCM4Random = new TH1D( "Phi_Scattered_335MeV_PosHelCM4Random", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 20, -4, 4);
    PhiSc405PosHelCM4Random = new TH1D( "Phi_Scattered_405MeV_PosHelCM4Random", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 20, -4, 4);
    PhiSc475PosHelCM4Random = new TH1D( "Phi_Scattered_475MeV_PosHelCM4Random", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 20, -4, 4);
    PhiSc545PosHelCM4Random = new TH1D( "Phi_Scattered_545MeV_PosHelCM4Random", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 20, -4, 4);
    PhiSc615PosHelCM4Random = new TH1D( "Phi_Scattered_615MeV_PosHelCM4Random", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 20, -4, 4);
    PhiSc685PosHelCM4Random = new TH1D( "Phi_Scattered_685MeV_PosHelCM4Random", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.25-0.0)) for +ve Helicity", 20, -4, 4);

    PhiSc265PosHelCM5Random = new TH1D( "Phi_Scattered_265MeV_PosHelCM5Random", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for +ve Helicity", 20, -4, 4);
    PhiSc335PosHelCM5Random = new TH1D( "Phi_Scattered_335MeV_PosHelCM5Random", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for +ve Helicity", 20, -4, 4);
    PhiSc405PosHelCM5Random = new TH1D( "Phi_Scattered_405MeV_PosHelCM5Random", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for +ve Helicity", 20, -4, 4);
    PhiSc475PosHelCM5Random = new TH1D( "Phi_Scattered_475MeV_PosHelCM5Random", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for +ve Helicity", 20, -4, 4);
    PhiSc545PosHelCM5Random = new TH1D( "Phi_Scattered_545MeV_PosHelCM5Random", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for +ve Helicity", 20, -4, 4);
    PhiSc615PosHelCM5Random = new TH1D( "Phi_Scattered_615MeV_PosHelCM5Random", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for +ve Helicity", 20, -4, 4);
    PhiSc685PosHelCM5Random = new TH1D( "Phi_Scattered_685MeV_PosHelCM5Random", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}0.0-(-0.25))) for +ve Helicity", 20, -4, 4);

    PhiSc265PosHelCM6Random = new TH1D( "Phi_Scattered_265MeV_PosHelCM6Random", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 20, -4, 4);
    PhiSc335PosHelCM6Random = new TH1D( "Phi_Scattered_335MeV_PosHelCM6Random", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 20, -4, 4);
    PhiSc405PosHelCM6Random = new TH1D( "Phi_Scattered_405MeV_PosHelCM6Random", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 20, -4, 4);
    PhiSc475PosHelCM6Random = new TH1D( "Phi_Scattered_475MeV_PosHelCM6Random", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 20, -4, 4);
    PhiSc545PosHelCM6Random = new TH1D( "Phi_Scattered_545MeV_PosHelCM6Random", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 20, -4, 4);
    PhiSc615PosHelCM6Random = new TH1D( "Phi_Scattered_615MeV_PosHelCM6Random", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 20, -4, 4);
    PhiSc685PosHelCM6Random = new TH1D( "Phi_Scattered_685MeV_PosHelCM6Random", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}-0.25-(-0.5))) for +ve Helicity", 20, -4, 4);

    PhiSc265PosHelCM7Random = new TH1D( "Phi_Scattered_265MeV_PosHelCM7Random", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 20, -4, 4);
    PhiSc335PosHelCM7Random = new TH1D( "Phi_Scattered_335MeV_PosHelCM7Random", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 20, -4, 4);
    PhiSc405PosHelCM7Random = new TH1D( "Phi_Scattered_405MeV_PosHelCM7Random", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 20, -4, 4);
    PhiSc475PosHelCM7Random = new TH1D( "Phi_Scattered_475MeV_PosHelCM7Random", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 20, -4, 4);
    PhiSc545PosHelCM7Random = new TH1D( "Phi_Scattered_545MeV_PosHelCM7Random", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 20, -4, 4);
    PhiSc615PosHelCM7Random = new TH1D( "Phi_Scattered_615MeV_PosHelCM7Random", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 20, -4, 4);
    PhiSc685PosHelCM7Random = new TH1D( "Phi_Scattered_685MeV_PosHelCM7Random", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}-0.5-(-0.75))) for +ve Helicity", 20, -4, 4);

    PhiSc265PosHelCM8Random = new TH1D( "Phi_Scattered_265MeV_PosHelCM8Random", "#phi_{Sc} E_{#gamma}265 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 20, -4, 4);
    PhiSc335PosHelCM8Random = new TH1D( "Phi_Scattered_335MeV_PosHelCM8Random", "#phi_{Sc} E_{#gamma}335 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 20, -4, 4);
    PhiSc405PosHelCM8Random = new TH1D( "Phi_Scattered_405MeV_PosHelCM8Random", "#phi_{Sc} E_{#gamma}405 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 20, -4, 4);
    PhiSc475PosHelCM8Random = new TH1D( "Phi_Scattered_475MeV_PosHelCM8Random", "#phi_{Sc} E_{#gamma}475 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 20, -4, 4);
    PhiSc545PosHelCM8Random = new TH1D( "Phi_Scattered_545MeV_PosHelCM8Random", "#phi_{Sc} E_{#gamma}545 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 20, -4, 4);
    PhiSc615PosHelCM8Random = new TH1D( "Phi_Scattered_615MeV_PosHelCM8Random", "#phi_{Sc} E_{#gamma}615 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 20, -4, 4);
    PhiSc685PosHelCM8Random = new TH1D( "Phi_Scattered_685MeV_PosHelCM8Random", "#phi_{Sc} E_{#gamma}685 #pm 35MeV (Cos#theta_{CM}-0.75-(-1.0))) for +ve Helicity", 20, -4, 4);
}

void PNeutPol_Polarimeter_Circ::FillHists()
{
    if (GHistBGSub::IsPrompt(TaggerTime) == kTRUE){
        EgPrompt->Fill(EGamma);
        PhiScPrompt->Fill(ScattPhiRad);
        //if(ThetanCM  > 80 && ThetanCM < 100){
            if ( 300 < EGamma && EGamma < 340) PhiSc320Prompt->Fill(ScattPhiRad);
            if ( 340 < EGamma && EGamma < 380) PhiSc360Prompt->Fill(ScattPhiRad);
            if ( 380 < EGamma && EGamma < 420) PhiSc400Prompt->Fill(ScattPhiRad);
            if ( 420 < EGamma && EGamma < 460) PhiSc440Prompt->Fill(ScattPhiRad);
            if ( 460 < EGamma && EGamma < 500) PhiSc480Prompt->Fill(ScattPhiRad);
            if ( 500 < EGamma && EGamma < 540) PhiSc520Prompt->Fill(ScattPhiRad);
            if ( 540 < EGamma && EGamma < 580) PhiSc560Prompt->Fill(ScattPhiRad);
            if ( 580 < EGamma && EGamma < 620) PhiSc600Prompt->Fill(ScattPhiRad);
            if ( 620 < EGamma && EGamma < 660) PhiSc640Prompt->Fill(ScattPhiRad);
            if ( 660 < EGamma && EGamma < 700) PhiSc680Prompt->Fill(ScattPhiRad);

            if ( 230 < EGamma && EGamma < 300) {
                if(1 > CosThetapCM && CosThetapCM > 0.75 ){
                    if (BeamHelicity == kFALSE) PhiSc265NegHelCM1Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc265PosHelCM1Prompt->Fill(ScattPhiRad);
                }

                else if(0.75 > CosThetapCM && CosThetapCM > 0.5){
                    if (BeamHelicity == kFALSE) PhiSc265NegHelCM2Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc265PosHelCM2Prompt->Fill(ScattPhiRad);
                }

                else if(0.5 > CosThetapCM && CosThetapCM > 0.25){
                    if (BeamHelicity == kFALSE) PhiSc265NegHelCM3Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc265PosHelCM3Prompt->Fill(ScattPhiRad);
                }

                else if(0.5 > CosThetapCM && CosThetapCM > 0.0){
                    if (BeamHelicity == kFALSE) PhiSc265NegHelCM4Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc265PosHelCM4Prompt->Fill(ScattPhiRad);
                }

                else if(0.0 > CosThetapCM && CosThetapCM > -0.25){
                    if (BeamHelicity == kFALSE) PhiSc265NegHelCM5Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc265PosHelCM5Prompt->Fill(ScattPhiRad);
                }

                else if(-0.25 > CosThetapCM && CosThetapCM > -0.5){
                    if (BeamHelicity == kFALSE) PhiSc265NegHelCM6Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc265PosHelCM6Prompt->Fill(ScattPhiRad);
                }

                else if(-0.5 > CosThetapCM && CosThetapCM > -0.75){
                    if (BeamHelicity == kFALSE) PhiSc265NegHelCM7Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc265PosHelCM7Prompt->Fill(ScattPhiRad);
                }

                else if(-0.75 > CosThetapCM && CosThetapCM > -1){
                    if (BeamHelicity == kFALSE) PhiSc265NegHelCM8Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc265PosHelCM8Prompt->Fill(ScattPhiRad);
                }
            }


            if ( 300 < EGamma && EGamma < 370) {
                if(1 > CosThetapCM && CosThetapCM > 0.75 ){
                    if (BeamHelicity == kFALSE) PhiSc335NegHelCM1Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc335PosHelCM1Prompt->Fill(ScattPhiRad);
                }

                else if(0.75 > CosThetapCM && CosThetapCM > 0.5){
                    if (BeamHelicity == kFALSE) PhiSc335NegHelCM2Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc335PosHelCM2Prompt->Fill(ScattPhiRad);
                }

                else if(0.5 > CosThetapCM && CosThetapCM > 0.25){
                    if (BeamHelicity == kFALSE) PhiSc335NegHelCM3Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc335PosHelCM3Prompt->Fill(ScattPhiRad);
                }

                else if(0.5 > CosThetapCM && CosThetapCM > 0.0){
                    if (BeamHelicity == kFALSE) PhiSc335NegHelCM4Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc335PosHelCM4Prompt->Fill(ScattPhiRad);
                }

                else if(0.0 > CosThetapCM && CosThetapCM > -0.25){
                    if (BeamHelicity == kFALSE) PhiSc335NegHelCM5Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc335PosHelCM5Prompt->Fill(ScattPhiRad);
                }

                else if(-0.25 > CosThetapCM && CosThetapCM > -0.5){
                    if (BeamHelicity == kFALSE) PhiSc335NegHelCM6Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc335PosHelCM6Prompt->Fill(ScattPhiRad);
                }

                else if(-0.5 > CosThetapCM && CosThetapCM > -0.75){
                    if (BeamHelicity == kFALSE) PhiSc335NegHelCM7Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc335PosHelCM7Prompt->Fill(ScattPhiRad);
                }

                else if(-0.75 > CosThetapCM && CosThetapCM > -1){
                    if (BeamHelicity == kFALSE) PhiSc335NegHelCM8Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc335PosHelCM8Prompt->Fill(ScattPhiRad);
                }
            }

            if ( 370 < EGamma && EGamma < 440) {
                if(1 > CosThetapCM && CosThetapCM > 0.75 ){
                    if (BeamHelicity == kFALSE) PhiSc405NegHelCM1Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc405PosHelCM1Prompt->Fill(ScattPhiRad);
                }

                else if(0.75 > CosThetapCM && CosThetapCM > 0.5){
                    if (BeamHelicity == kFALSE) PhiSc405NegHelCM2Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc405PosHelCM2Prompt->Fill(ScattPhiRad);
                }

                else if(0.5 > CosThetapCM && CosThetapCM > 0.25){
                    if (BeamHelicity == kFALSE) PhiSc405NegHelCM3Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc405PosHelCM3Prompt->Fill(ScattPhiRad);
                }

                else if(0.5 > CosThetapCM && CosThetapCM > 0.0){
                    if (BeamHelicity == kFALSE) PhiSc405NegHelCM4Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc405PosHelCM4Prompt->Fill(ScattPhiRad);
                }

                else if(0.0 > CosThetapCM && CosThetapCM > -0.25){
                    if (BeamHelicity == kFALSE) PhiSc405NegHelCM5Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc405PosHelCM5Prompt->Fill(ScattPhiRad);
                }

                else if(-0.25 > CosThetapCM && CosThetapCM > -0.5){
                    if (BeamHelicity == kFALSE) PhiSc405NegHelCM6Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc405PosHelCM6Prompt->Fill(ScattPhiRad);
                }

                else if(-0.5 > CosThetapCM && CosThetapCM > -0.75){
                    if (BeamHelicity == kFALSE) PhiSc405NegHelCM7Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc405PosHelCM7Prompt->Fill(ScattPhiRad);
                }

                else if(-0.75 > CosThetapCM && CosThetapCM > -1){
                    if (BeamHelicity == kFALSE) PhiSc405NegHelCM8Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc405PosHelCM8Prompt->Fill(ScattPhiRad);
                }
            }

            if ( 440 < EGamma && EGamma < 510) {
                if(1 > CosThetapCM && CosThetapCM > 0.75 ){
                    if (BeamHelicity == kFALSE) PhiSc475NegHelCM1Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc475PosHelCM1Prompt->Fill(ScattPhiRad);
                }

                else if(0.75 > CosThetapCM && CosThetapCM > 0.5){
                    if (BeamHelicity == kFALSE) PhiSc475NegHelCM2Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc475PosHelCM2Prompt->Fill(ScattPhiRad);
                }

                else if(0.5 > CosThetapCM && CosThetapCM > 0.25){
                    if (BeamHelicity == kFALSE) PhiSc475NegHelCM3Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc475PosHelCM3Prompt->Fill(ScattPhiRad);
                }

                else if(0.5 > CosThetapCM && CosThetapCM > 0.0){
                    if (BeamHelicity == kFALSE) PhiSc475NegHelCM4Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc475PosHelCM4Prompt->Fill(ScattPhiRad);
                }

                else if(0.0 > CosThetapCM && CosThetapCM > -0.25){
                    if (BeamHelicity == kFALSE) PhiSc475NegHelCM5Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc475PosHelCM5Prompt->Fill(ScattPhiRad);
                }

                else if(-0.25 > CosThetapCM && CosThetapCM > -0.5){
                    if (BeamHelicity == kFALSE) PhiSc475NegHelCM6Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc475PosHelCM6Prompt->Fill(ScattPhiRad);
                }

                else if(-0.5 > CosThetapCM && CosThetapCM > -0.75){
                    if (BeamHelicity == kFALSE) PhiSc475NegHelCM7Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc475PosHelCM7Prompt->Fill(ScattPhiRad);
                }

                else if(-0.75 > CosThetapCM && CosThetapCM > -1){
                    if (BeamHelicity == kFALSE) PhiSc475NegHelCM8Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc475PosHelCM8Prompt->Fill(ScattPhiRad);
                }
            }

            if ( 510 < EGamma && EGamma < 580) {
                if(1 > CosThetapCM && CosThetapCM > 0.75 ){
                    if (BeamHelicity == kFALSE) PhiSc545NegHelCM1Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc545PosHelCM1Prompt->Fill(ScattPhiRad);
                }

                else if(0.75 > CosThetapCM && CosThetapCM > 0.5){
                    if (BeamHelicity == kFALSE) PhiSc545NegHelCM2Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc545PosHelCM2Prompt->Fill(ScattPhiRad);
                }

                else if(0.5 > CosThetapCM && CosThetapCM > 0.25){
                    if (BeamHelicity == kFALSE) PhiSc545NegHelCM3Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc545PosHelCM3Prompt->Fill(ScattPhiRad);
                }

                else if(0.5 > CosThetapCM && CosThetapCM > 0.0){
                    if (BeamHelicity == kFALSE) PhiSc545NegHelCM4Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc545PosHelCM4Prompt->Fill(ScattPhiRad);
                }

                else if(0.0 > CosThetapCM && CosThetapCM > -0.25){
                    if (BeamHelicity == kFALSE) PhiSc545NegHelCM5Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc545PosHelCM5Prompt->Fill(ScattPhiRad);
                }

                else if(-0.25 > CosThetapCM && CosThetapCM > -0.5){
                    if (BeamHelicity == kFALSE) PhiSc545NegHelCM6Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc545PosHelCM6Prompt->Fill(ScattPhiRad);
                }

                else if(-0.5 > CosThetapCM && CosThetapCM > -0.75){
                    if (BeamHelicity == kFALSE) PhiSc545NegHelCM7Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc545PosHelCM7Prompt->Fill(ScattPhiRad);
                }

                else if(-0.75 > CosThetapCM && CosThetapCM > -1){
                    if (BeamHelicity == kFALSE) PhiSc545NegHelCM8Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc545PosHelCM8Prompt->Fill(ScattPhiRad);
                }
            }

            if ( 580 < EGamma && EGamma < 650) {
                if(1 > CosThetapCM && CosThetapCM > 0.75 ){
                    if (BeamHelicity == kFALSE) PhiSc615NegHelCM1Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc615PosHelCM1Prompt->Fill(ScattPhiRad);
                }

                else if(0.75 > CosThetapCM && CosThetapCM > 0.5){
                    if (BeamHelicity == kFALSE) PhiSc615NegHelCM2Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc615PosHelCM2Prompt->Fill(ScattPhiRad);
                }

                else if(0.5 > CosThetapCM && CosThetapCM > 0.25){
                    if (BeamHelicity == kFALSE) PhiSc615NegHelCM3Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc615PosHelCM3Prompt->Fill(ScattPhiRad);
                }

                else if(0.5 > CosThetapCM && CosThetapCM > 0.0){
                    if (BeamHelicity == kFALSE) PhiSc615NegHelCM4Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc615PosHelCM4Prompt->Fill(ScattPhiRad);
                }

                else if(0.0 > CosThetapCM && CosThetapCM > -0.25){
                    if (BeamHelicity == kFALSE) PhiSc615NegHelCM5Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc615PosHelCM5Prompt->Fill(ScattPhiRad);
                }

                else if(-0.25 > CosThetapCM && CosThetapCM > -0.5){
                    if (BeamHelicity == kFALSE) PhiSc615NegHelCM6Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc615PosHelCM6Prompt->Fill(ScattPhiRad);
                }

                else if(-0.5 > CosThetapCM && CosThetapCM > -0.75){
                    if (BeamHelicity == kFALSE) PhiSc615NegHelCM7Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc615PosHelCM7Prompt->Fill(ScattPhiRad);
                }

                else if(-0.75 > CosThetapCM && CosThetapCM > -1){
                    if (BeamHelicity == kFALSE) PhiSc615NegHelCM8Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc615PosHelCM8Prompt->Fill(ScattPhiRad);
                }
            }

            if ( 650 < EGamma && EGamma < 720) {
                if(1 > CosThetapCM && CosThetapCM > 0.75 ){
                    if (BeamHelicity == kFALSE) PhiSc685NegHelCM1Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc685PosHelCM1Prompt->Fill(ScattPhiRad);
                }

                else if(0.75 > CosThetapCM && CosThetapCM > 0.5){
                    if (BeamHelicity == kFALSE) PhiSc685NegHelCM2Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc685PosHelCM2Prompt->Fill(ScattPhiRad);
                }

                else if(0.5 > CosThetapCM && CosThetapCM > 0.25){
                    if (BeamHelicity == kFALSE) PhiSc685NegHelCM3Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc685PosHelCM3Prompt->Fill(ScattPhiRad);
                }

                else if(0.5 > CosThetapCM && CosThetapCM > 0.0){
                    if (BeamHelicity == kFALSE) PhiSc685NegHelCM4Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc685PosHelCM4Prompt->Fill(ScattPhiRad);
                }

                else if(0.0 > CosThetapCM && CosThetapCM > -0.25){
                    if (BeamHelicity == kFALSE) PhiSc685NegHelCM5Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc685PosHelCM5Prompt->Fill(ScattPhiRad);
                }

                else if(-0.25 > CosThetapCM && CosThetapCM > -0.5){
                    if (BeamHelicity == kFALSE) PhiSc685NegHelCM6Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc685PosHelCM6Prompt->Fill(ScattPhiRad);
                }

                else if(-0.5 > CosThetapCM && CosThetapCM > -0.75){
                    if (BeamHelicity == kFALSE) PhiSc685NegHelCM7Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc685PosHelCM7Prompt->Fill(ScattPhiRad);
                }

                else if(-0.75 > CosThetapCM && CosThetapCM > -1){
                    if (BeamHelicity == kFALSE) PhiSc685NegHelCM8Prompt->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc685PosHelCM8Prompt->Fill(ScattPhiRad);
                }
            }
        //}
    }

    if (GHistBGSub::IsRandom(TaggerTime) == kTRUE){
        EgRandom->Fill(EGamma);
        PhiScRandom->Fill(ScattPhiRad);
        //if(ThetanCM  > 80 && ThetanCM < 100){
            if ( 300 < EGamma && EGamma < 340) PhiSc320Random->Fill(ScattPhiRad);
            if ( 340 < EGamma && EGamma < 380) PhiSc360Random->Fill(ScattPhiRad);
            if ( 380 < EGamma && EGamma < 420) PhiSc400Random->Fill(ScattPhiRad);
            if ( 420 < EGamma && EGamma < 460) PhiSc440Random->Fill(ScattPhiRad);
            if ( 460 < EGamma && EGamma < 500) PhiSc480Random->Fill(ScattPhiRad);
            if ( 500 < EGamma && EGamma < 540) PhiSc520Random->Fill(ScattPhiRad);
            if ( 540 < EGamma && EGamma < 580) PhiSc560Random->Fill(ScattPhiRad);
            if ( 580 < EGamma && EGamma < 620) PhiSc600Random->Fill(ScattPhiRad);
            if ( 620 < EGamma && EGamma < 660) PhiSc640Random->Fill(ScattPhiRad);
            if ( 660 < EGamma && EGamma < 700) PhiSc680Random->Fill(ScattPhiRad);

            if ( 230 < EGamma && EGamma < 300) {
                if(1 > CosThetapCM && CosThetapCM > 0.75 ){
                    if (BeamHelicity == kFALSE) PhiSc265NegHelCM1Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc265PosHelCM1Random->Fill(ScattPhiRad);
                }

                else if(0.75 > CosThetapCM && CosThetapCM > 0.5){
                    if (BeamHelicity == kFALSE) PhiSc265NegHelCM2Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc265PosHelCM2Random->Fill(ScattPhiRad);
                }

                else if(0.5 > CosThetapCM && CosThetapCM > 0.25){
                    if (BeamHelicity == kFALSE) PhiSc265NegHelCM3Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc265PosHelCM3Random->Fill(ScattPhiRad);
                }

                else if(0.5 > CosThetapCM && CosThetapCM > 0.0){
                    if (BeamHelicity == kFALSE) PhiSc265NegHelCM4Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc265PosHelCM4Random->Fill(ScattPhiRad);
                }

                else if(0.0 > CosThetapCM && CosThetapCM > -0.25){
                    if (BeamHelicity == kFALSE) PhiSc265NegHelCM5Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc265PosHelCM5Random->Fill(ScattPhiRad);
                }

                else if(-0.25 > CosThetapCM && CosThetapCM > -0.5){
                    if (BeamHelicity == kFALSE) PhiSc265NegHelCM6Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc265PosHelCM6Random->Fill(ScattPhiRad);
                }

                else if(-0.5 > CosThetapCM && CosThetapCM > -0.75){
                    if (BeamHelicity == kFALSE) PhiSc265NegHelCM7Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc265PosHelCM7Random->Fill(ScattPhiRad);
                }

                else if(-0.75 > CosThetapCM && CosThetapCM > -1){
                    if (BeamHelicity == kFALSE) PhiSc265NegHelCM8Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc265PosHelCM8Random->Fill(ScattPhiRad);
                }
            }


            if ( 300 < EGamma && EGamma < 370) {
                if(1 > CosThetapCM && CosThetapCM > 0.75 ){
                    if (BeamHelicity == kFALSE) PhiSc335NegHelCM1Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc335PosHelCM1Random->Fill(ScattPhiRad);
                }

                else if(0.75 > CosThetapCM && CosThetapCM > 0.5){
                    if (BeamHelicity == kFALSE) PhiSc335NegHelCM2Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc335PosHelCM2Random->Fill(ScattPhiRad);
                }

                else if(0.5 > CosThetapCM && CosThetapCM > 0.25){
                    if (BeamHelicity == kFALSE) PhiSc335NegHelCM3Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc335PosHelCM3Random->Fill(ScattPhiRad);
                }

                else if(0.5 > CosThetapCM && CosThetapCM > 0.0){
                    if (BeamHelicity == kFALSE) PhiSc335NegHelCM4Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc335PosHelCM4Random->Fill(ScattPhiRad);
                }

                else if(0.0 > CosThetapCM && CosThetapCM > -0.25){
                    if (BeamHelicity == kFALSE) PhiSc335NegHelCM5Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc335PosHelCM5Random->Fill(ScattPhiRad);
                }

                else if(-0.25 > CosThetapCM && CosThetapCM > -0.5){
                    if (BeamHelicity == kFALSE) PhiSc335NegHelCM6Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc335PosHelCM6Random->Fill(ScattPhiRad);
                }

                else if(-0.5 > CosThetapCM && CosThetapCM > -0.75){
                    if (BeamHelicity == kFALSE) PhiSc335NegHelCM7Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc335PosHelCM7Random->Fill(ScattPhiRad);
                }

                else if(-0.75 > CosThetapCM && CosThetapCM > -1){
                    if (BeamHelicity == kFALSE) PhiSc335NegHelCM8Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc335PosHelCM8Random->Fill(ScattPhiRad);
                }
            }

            if ( 370 < EGamma && EGamma < 440) {
                if(1 > CosThetapCM && CosThetapCM > 0.75 ){
                    if (BeamHelicity == kFALSE) PhiSc405NegHelCM1Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc405PosHelCM1Random->Fill(ScattPhiRad);
                }

                else if(0.75 > CosThetapCM && CosThetapCM > 0.5){
                    if (BeamHelicity == kFALSE) PhiSc405NegHelCM2Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc405PosHelCM2Random->Fill(ScattPhiRad);
                }

                else if(0.5 > CosThetapCM && CosThetapCM > 0.25){
                    if (BeamHelicity == kFALSE) PhiSc405NegHelCM3Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc405PosHelCM3Random->Fill(ScattPhiRad);
                }

                else if(0.5 > CosThetapCM && CosThetapCM > 0.0){
                    if (BeamHelicity == kFALSE) PhiSc405NegHelCM4Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc405PosHelCM4Random->Fill(ScattPhiRad);
                }

                else if(0.0 > CosThetapCM && CosThetapCM > -0.25){
                    if (BeamHelicity == kFALSE) PhiSc405NegHelCM5Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc405PosHelCM5Random->Fill(ScattPhiRad);
                }

                else if(-0.25 > CosThetapCM && CosThetapCM > -0.5){
                    if (BeamHelicity == kFALSE) PhiSc405NegHelCM6Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc405PosHelCM6Random->Fill(ScattPhiRad);
                }

                else if(-0.5 > CosThetapCM && CosThetapCM > -0.75){
                    if (BeamHelicity == kFALSE) PhiSc405NegHelCM7Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc405PosHelCM7Random->Fill(ScattPhiRad);
                }

                else if(-0.75 > CosThetapCM && CosThetapCM > -1){
                    if (BeamHelicity == kFALSE) PhiSc405NegHelCM8Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc405PosHelCM8Random->Fill(ScattPhiRad);
                }
            }

            if ( 440 < EGamma && EGamma < 510) {
                if(1 > CosThetapCM && CosThetapCM > 0.75 ){
                    if (BeamHelicity == kFALSE) PhiSc475NegHelCM1Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc475PosHelCM1Random->Fill(ScattPhiRad);
                }

                else if(0.75 > CosThetapCM && CosThetapCM > 0.5){
                    if (BeamHelicity == kFALSE) PhiSc475NegHelCM2Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc475PosHelCM2Random->Fill(ScattPhiRad);
                }

                else if(0.5 > CosThetapCM && CosThetapCM > 0.25){
                    if (BeamHelicity == kFALSE) PhiSc475NegHelCM3Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc475PosHelCM3Random->Fill(ScattPhiRad);
                }

                else if(0.5 > CosThetapCM && CosThetapCM > 0.0){
                    if (BeamHelicity == kFALSE) PhiSc475NegHelCM4Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc475PosHelCM4Random->Fill(ScattPhiRad);
                }

                else if(0.0 > CosThetapCM && CosThetapCM > -0.25){
                    if (BeamHelicity == kFALSE) PhiSc475NegHelCM5Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc475PosHelCM5Random->Fill(ScattPhiRad);
                }

                else if(-0.25 > CosThetapCM && CosThetapCM > -0.5){
                    if (BeamHelicity == kFALSE) PhiSc475NegHelCM6Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc475PosHelCM6Random->Fill(ScattPhiRad);
                }

                else if(-0.5 > CosThetapCM && CosThetapCM > -0.75){
                    if (BeamHelicity == kFALSE) PhiSc475NegHelCM7Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc475PosHelCM7Random->Fill(ScattPhiRad);
                }

                else if(-0.75 > CosThetapCM && CosThetapCM > -1){
                    if (BeamHelicity == kFALSE) PhiSc475NegHelCM8Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc475PosHelCM8Random->Fill(ScattPhiRad);
                }
            }

            if ( 510 < EGamma && EGamma < 580) {
                if(1 > CosThetapCM && CosThetapCM > 0.75 ){
                    if (BeamHelicity == kFALSE) PhiSc545NegHelCM1Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc545PosHelCM1Random->Fill(ScattPhiRad);
                }

                else if(0.75 > CosThetapCM && CosThetapCM > 0.5){
                    if (BeamHelicity == kFALSE) PhiSc545NegHelCM2Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc545PosHelCM2Random->Fill(ScattPhiRad);
                }

                else if(0.5 > CosThetapCM && CosThetapCM > 0.25){
                    if (BeamHelicity == kFALSE) PhiSc545NegHelCM3Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc545PosHelCM3Random->Fill(ScattPhiRad);
                }

                else if(0.5 > CosThetapCM && CosThetapCM > 0.0){
                    if (BeamHelicity == kFALSE) PhiSc545NegHelCM4Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc545PosHelCM4Random->Fill(ScattPhiRad);
                }

                else if(0.0 > CosThetapCM && CosThetapCM > -0.25){
                    if (BeamHelicity == kFALSE) PhiSc545NegHelCM5Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc545PosHelCM5Random->Fill(ScattPhiRad);
                }

                else if(-0.25 > CosThetapCM && CosThetapCM > -0.5){
                    if (BeamHelicity == kFALSE) PhiSc545NegHelCM6Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc545PosHelCM6Random->Fill(ScattPhiRad);
                }

                else if(-0.5 > CosThetapCM && CosThetapCM > -0.75){
                    if (BeamHelicity == kFALSE) PhiSc545NegHelCM7Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc545PosHelCM7Random->Fill(ScattPhiRad);
                }

                else if(-0.75 > CosThetapCM && CosThetapCM > -1){
                    if (BeamHelicity == kFALSE) PhiSc545NegHelCM8Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc545PosHelCM8Random->Fill(ScattPhiRad);
                }
            }

            if ( 580 < EGamma && EGamma < 650) {
                if(1 > CosThetapCM && CosThetapCM > 0.75 ){
                    if (BeamHelicity == kFALSE) PhiSc615NegHelCM1Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc615PosHelCM1Random->Fill(ScattPhiRad);
                }

                else if(0.75 > CosThetapCM && CosThetapCM > 0.5){
                    if (BeamHelicity == kFALSE) PhiSc615NegHelCM2Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc615PosHelCM2Random->Fill(ScattPhiRad);
                }

                else if(0.5 > CosThetapCM && CosThetapCM > 0.25){
                    if (BeamHelicity == kFALSE) PhiSc615NegHelCM3Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc615PosHelCM3Random->Fill(ScattPhiRad);
                }

                else if(0.5 > CosThetapCM && CosThetapCM > 0.0){
                    if (BeamHelicity == kFALSE) PhiSc615NegHelCM4Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc615PosHelCM4Random->Fill(ScattPhiRad);
                }

                else if(0.0 > CosThetapCM && CosThetapCM > -0.25){
                    if (BeamHelicity == kFALSE) PhiSc615NegHelCM5Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc615PosHelCM5Random->Fill(ScattPhiRad);
                }

                else if(-0.25 > CosThetapCM && CosThetapCM > -0.5){
                    if (BeamHelicity == kFALSE) PhiSc615NegHelCM6Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc615PosHelCM6Random->Fill(ScattPhiRad);
                }

                else if(-0.5 > CosThetapCM && CosThetapCM > -0.75){
                    if (BeamHelicity == kFALSE) PhiSc615NegHelCM7Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc615PosHelCM7Random->Fill(ScattPhiRad);
                }

                else if(-0.75 > CosThetapCM && CosThetapCM > -1){
                    if (BeamHelicity == kFALSE) PhiSc615NegHelCM8Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc615PosHelCM8Random->Fill(ScattPhiRad);
                }
            }

            if ( 650 < EGamma && EGamma < 720) {
                if(1 > CosThetapCM && CosThetapCM > 0.75 ){
                    if (BeamHelicity == kFALSE) PhiSc685NegHelCM1Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc685PosHelCM1Random->Fill(ScattPhiRad);
                }

                else if(0.75 > CosThetapCM && CosThetapCM > 0.5){
                    if (BeamHelicity == kFALSE) PhiSc685NegHelCM2Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc685PosHelCM2Random->Fill(ScattPhiRad);
                }

                else if(0.5 > CosThetapCM && CosThetapCM > 0.25){
                    if (BeamHelicity == kFALSE) PhiSc685NegHelCM3Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc685PosHelCM3Random->Fill(ScattPhiRad);
                }

                else if(0.5 > CosThetapCM && CosThetapCM > 0.0){
                    if (BeamHelicity == kFALSE) PhiSc685NegHelCM4Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc685PosHelCM4Random->Fill(ScattPhiRad);
                }

                else if(0.0 > CosThetapCM && CosThetapCM > -0.25){
                    if (BeamHelicity == kFALSE) PhiSc685NegHelCM5Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc685PosHelCM5Random->Fill(ScattPhiRad);
                }

                else if(-0.25 > CosThetapCM && CosThetapCM > -0.5){
                    if (BeamHelicity == kFALSE) PhiSc685NegHelCM6Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc685PosHelCM6Random->Fill(ScattPhiRad);
                }

                else if(-0.5 > CosThetapCM && CosThetapCM > -0.75){
                    if (BeamHelicity == kFALSE) PhiSc685NegHelCM7Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc685PosHelCM7Random->Fill(ScattPhiRad);
                }

                else if(-0.75 > CosThetapCM && CosThetapCM > -1){
                    if (BeamHelicity == kFALSE) PhiSc685NegHelCM8Random->Fill(ScattPhiRad);
                    else if (BeamHelicity == kTRUE) PhiSc685PosHelCM8Random->Fill(ScattPhiRad);
                }
            }
        //}
    }
}

void PNeutPol_Polarimeter_Circ::BGSub(){

    Eg->Add(EgPrompt);
    Eg->Add(EgRandom, -PvRratio);
    PhiSc->Add(PhiScPrompt);
    PhiSc->Add(PhiScRandom, -PvRratio);
    PhiSc320->Add(PhiSc320Prompt);
    PhiSc320->Add(PhiSc320Random, -PvRratio);
    PhiSc360->Add(PhiSc360Prompt);
    PhiSc360->Add(PhiSc360Random, -PvRratio);
    PhiSc400->Add(PhiSc400Prompt);
    PhiSc400->Add(PhiSc400Random, -PvRratio);
    PhiSc440->Add(PhiSc440Prompt);
    PhiSc440->Add(PhiSc440Random, -PvRratio);
    PhiSc480->Add(PhiSc480Prompt);
    PhiSc480->Add(PhiSc480Random, -PvRratio);
    PhiSc520->Add(PhiSc520Prompt);
    PhiSc520->Add(PhiSc520Random, -PvRratio);
    PhiSc560->Add(PhiSc560Prompt);
    PhiSc560->Add(PhiSc560Random, -PvRratio);
    PhiSc600->Add(PhiSc600Prompt);
    PhiSc600->Add(PhiSc600Random, -PvRratio);
    PhiSc640->Add(PhiSc640Prompt);
    PhiSc640->Add(PhiSc640Random, -PvRratio);
    PhiSc680->Add(PhiSc680Prompt);
    PhiSc680->Add(PhiSc680Random, -PvRratio);

    PhiSc265NegHelCM1->Add(PhiSc265NegHelCM1Prompt);
    PhiSc335NegHelCM1->Add(PhiSc335NegHelCM1Prompt);
    PhiSc405NegHelCM1->Add(PhiSc405NegHelCM1Prompt);
    PhiSc475NegHelCM1->Add(PhiSc475NegHelCM1Prompt);
    PhiSc545NegHelCM1->Add(PhiSc545NegHelCM1Prompt);
    PhiSc615NegHelCM1->Add(PhiSc615NegHelCM1Prompt);
    PhiSc685NegHelCM1->Add(PhiSc685NegHelCM1Prompt);
    PhiSc265NegHelCM2->Add(PhiSc265NegHelCM2Prompt);
    PhiSc335NegHelCM2->Add(PhiSc335NegHelCM2Prompt);
    PhiSc405NegHelCM2->Add(PhiSc405NegHelCM2Prompt);
    PhiSc475NegHelCM2->Add(PhiSc475NegHelCM2Prompt);
    PhiSc545NegHelCM2->Add(PhiSc545NegHelCM2Prompt);
    PhiSc615NegHelCM2->Add(PhiSc615NegHelCM2Prompt);
    PhiSc685NegHelCM2->Add(PhiSc685NegHelCM2Prompt);
    PhiSc265NegHelCM3->Add(PhiSc265NegHelCM3Prompt);
    PhiSc335NegHelCM3->Add(PhiSc335NegHelCM3Prompt);
    PhiSc405NegHelCM3->Add(PhiSc405NegHelCM3Prompt);
    PhiSc475NegHelCM3->Add(PhiSc475NegHelCM3Prompt);
    PhiSc545NegHelCM3->Add(PhiSc545NegHelCM3Prompt);
    PhiSc615NegHelCM3->Add(PhiSc615NegHelCM3Prompt);
    PhiSc685NegHelCM3->Add(PhiSc685NegHelCM3Prompt);
    PhiSc265NegHelCM4->Add(PhiSc265NegHelCM4Prompt);
    PhiSc335NegHelCM4->Add(PhiSc335NegHelCM4Prompt);
    PhiSc405NegHelCM4->Add(PhiSc405NegHelCM4Prompt);
    PhiSc475NegHelCM4->Add(PhiSc475NegHelCM4Prompt);
    PhiSc545NegHelCM4->Add(PhiSc545NegHelCM4Prompt);
    PhiSc615NegHelCM4->Add(PhiSc615NegHelCM4Prompt);
    PhiSc685NegHelCM4->Add(PhiSc685NegHelCM4Prompt);
    PhiSc265NegHelCM5->Add(PhiSc265NegHelCM5Prompt);
    PhiSc335NegHelCM5->Add(PhiSc335NegHelCM5Prompt);
    PhiSc405NegHelCM5->Add(PhiSc405NegHelCM5Prompt);
    PhiSc475NegHelCM5->Add(PhiSc475NegHelCM5Prompt);
    PhiSc545NegHelCM5->Add(PhiSc545NegHelCM5Prompt);
    PhiSc615NegHelCM5->Add(PhiSc615NegHelCM5Prompt);
    PhiSc685NegHelCM5->Add(PhiSc685NegHelCM5Prompt);
    PhiSc265NegHelCM6->Add(PhiSc265NegHelCM6Prompt);
    PhiSc335NegHelCM6->Add(PhiSc335NegHelCM6Prompt);
    PhiSc405NegHelCM6->Add(PhiSc405NegHelCM6Prompt);
    PhiSc475NegHelCM6->Add(PhiSc475NegHelCM6Prompt);
    PhiSc545NegHelCM6->Add(PhiSc545NegHelCM6Prompt);
    PhiSc615NegHelCM6->Add(PhiSc615NegHelCM6Prompt);
    PhiSc685NegHelCM6->Add(PhiSc685NegHelCM6Prompt);
    PhiSc265NegHelCM7->Add(PhiSc265NegHelCM7Prompt);
    PhiSc335NegHelCM7->Add(PhiSc335NegHelCM7Prompt);
    PhiSc405NegHelCM7->Add(PhiSc405NegHelCM7Prompt);
    PhiSc475NegHelCM7->Add(PhiSc475NegHelCM7Prompt);
    PhiSc545NegHelCM7->Add(PhiSc545NegHelCM7Prompt);
    PhiSc615NegHelCM7->Add(PhiSc615NegHelCM7Prompt);
    PhiSc685NegHelCM7->Add(PhiSc685NegHelCM7Prompt);
    PhiSc265NegHelCM8->Add(PhiSc265NegHelCM8Prompt);
    PhiSc335NegHelCM8->Add(PhiSc335NegHelCM8Prompt);
    PhiSc405NegHelCM8->Add(PhiSc405NegHelCM8Prompt);
    PhiSc475NegHelCM8->Add(PhiSc475NegHelCM8Prompt);
    PhiSc545NegHelCM8->Add(PhiSc545NegHelCM8Prompt);
    PhiSc615NegHelCM8->Add(PhiSc615NegHelCM8Prompt);
    PhiSc685NegHelCM8->Add(PhiSc685NegHelCM8Prompt);
    PhiSc265PosHelCM1->Add(PhiSc265PosHelCM1Prompt);
    PhiSc335PosHelCM1->Add(PhiSc335PosHelCM1Prompt);
    PhiSc405PosHelCM1->Add(PhiSc405PosHelCM1Prompt);
    PhiSc475PosHelCM1->Add(PhiSc475PosHelCM1Prompt);
    PhiSc545PosHelCM1->Add(PhiSc545PosHelCM1Prompt);
    PhiSc615PosHelCM1->Add(PhiSc615PosHelCM1Prompt);
    PhiSc685PosHelCM1->Add(PhiSc685PosHelCM1Prompt);
    PhiSc265PosHelCM2->Add(PhiSc265PosHelCM2Prompt);
    PhiSc335PosHelCM2->Add(PhiSc335PosHelCM2Prompt);
    PhiSc405PosHelCM2->Add(PhiSc405PosHelCM2Prompt);
    PhiSc475PosHelCM2->Add(PhiSc475PosHelCM2Prompt);
    PhiSc545PosHelCM2->Add(PhiSc545PosHelCM2Prompt);
    PhiSc615PosHelCM2->Add(PhiSc615PosHelCM2Prompt);
    PhiSc685PosHelCM2->Add(PhiSc685PosHelCM2Prompt);
    PhiSc265PosHelCM3->Add(PhiSc265PosHelCM3Prompt);
    PhiSc335PosHelCM3->Add(PhiSc335PosHelCM3Prompt);
    PhiSc405PosHelCM3->Add(PhiSc405PosHelCM3Prompt);
    PhiSc475PosHelCM3->Add(PhiSc475PosHelCM3Prompt);
    PhiSc545PosHelCM3->Add(PhiSc545PosHelCM3Prompt);
    PhiSc615PosHelCM3->Add(PhiSc615PosHelCM3Prompt);
    PhiSc685PosHelCM3->Add(PhiSc685PosHelCM3Prompt);
    PhiSc265PosHelCM4->Add(PhiSc265PosHelCM4Prompt);
    PhiSc335PosHelCM4->Add(PhiSc335PosHelCM4Prompt);
    PhiSc405PosHelCM4->Add(PhiSc405PosHelCM4Prompt);
    PhiSc475PosHelCM4->Add(PhiSc475PosHelCM4Prompt);
    PhiSc545PosHelCM4->Add(PhiSc545PosHelCM4Prompt);
    PhiSc615PosHelCM4->Add(PhiSc615PosHelCM4Prompt);
    PhiSc685PosHelCM4->Add(PhiSc685PosHelCM4Prompt);
    PhiSc265PosHelCM5->Add(PhiSc265PosHelCM5Prompt);
    PhiSc335PosHelCM5->Add(PhiSc335PosHelCM5Prompt);
    PhiSc405PosHelCM5->Add(PhiSc405PosHelCM5Prompt);
    PhiSc475PosHelCM5->Add(PhiSc475PosHelCM5Prompt);
    PhiSc545PosHelCM5->Add(PhiSc545PosHelCM5Prompt);
    PhiSc615PosHelCM5->Add(PhiSc615PosHelCM5Prompt);
    PhiSc685PosHelCM5->Add(PhiSc685PosHelCM5Prompt);
    PhiSc265PosHelCM6->Add(PhiSc265PosHelCM6Prompt);
    PhiSc335PosHelCM6->Add(PhiSc335PosHelCM6Prompt);
    PhiSc405PosHelCM6->Add(PhiSc405PosHelCM6Prompt);
    PhiSc475PosHelCM6->Add(PhiSc475PosHelCM6Prompt);
    PhiSc545PosHelCM6->Add(PhiSc545PosHelCM6Prompt);
    PhiSc615PosHelCM6->Add(PhiSc615PosHelCM6Prompt);
    PhiSc685PosHelCM6->Add(PhiSc685PosHelCM6Prompt);
    PhiSc265PosHelCM7->Add(PhiSc265PosHelCM7Prompt);
    PhiSc335PosHelCM7->Add(PhiSc335PosHelCM7Prompt);
    PhiSc405PosHelCM7->Add(PhiSc405PosHelCM7Prompt);
    PhiSc475PosHelCM7->Add(PhiSc475PosHelCM7Prompt);
    PhiSc545PosHelCM7->Add(PhiSc545PosHelCM7Prompt);
    PhiSc615PosHelCM7->Add(PhiSc615PosHelCM7Prompt);
    PhiSc685PosHelCM7->Add(PhiSc685PosHelCM7Prompt);
    PhiSc265PosHelCM8->Add(PhiSc265PosHelCM8Prompt);
    PhiSc335PosHelCM8->Add(PhiSc335PosHelCM8Prompt);
    PhiSc405PosHelCM8->Add(PhiSc405PosHelCM8Prompt);
    PhiSc475PosHelCM8->Add(PhiSc475PosHelCM8Prompt);
    PhiSc545PosHelCM8->Add(PhiSc545PosHelCM8Prompt);
    PhiSc615PosHelCM8->Add(PhiSc615PosHelCM8Prompt);
    PhiSc685PosHelCM8->Add(PhiSc685PosHelCM8Prompt);

    PhiSc265NegHelCM1->Add(PhiSc265NegHelCM1Random, -PvRratio);
    PhiSc335NegHelCM1->Add(PhiSc335NegHelCM1Random, -PvRratio);
    PhiSc405NegHelCM1->Add(PhiSc405NegHelCM1Random, -PvRratio);
    PhiSc475NegHelCM1->Add(PhiSc475NegHelCM1Random, -PvRratio);
    PhiSc545NegHelCM1->Add(PhiSc545NegHelCM1Random, -PvRratio);
    PhiSc615NegHelCM1->Add(PhiSc615NegHelCM1Random, -PvRratio);
    PhiSc685NegHelCM1->Add(PhiSc685NegHelCM1Random, -PvRratio);
    PhiSc265NegHelCM2->Add(PhiSc265NegHelCM2Random, -PvRratio);
    PhiSc335NegHelCM2->Add(PhiSc335NegHelCM2Random, -PvRratio);
    PhiSc405NegHelCM2->Add(PhiSc405NegHelCM2Random, -PvRratio);
    PhiSc475NegHelCM2->Add(PhiSc475NegHelCM2Random, -PvRratio);
    PhiSc545NegHelCM2->Add(PhiSc545NegHelCM2Random, -PvRratio);
    PhiSc615NegHelCM2->Add(PhiSc615NegHelCM2Random, -PvRratio);
    PhiSc685NegHelCM2->Add(PhiSc685NegHelCM2Random, -PvRratio);
    PhiSc265NegHelCM3->Add(PhiSc265NegHelCM3Random, -PvRratio);
    PhiSc335NegHelCM3->Add(PhiSc335NegHelCM3Random, -PvRratio);
    PhiSc405NegHelCM3->Add(PhiSc405NegHelCM3Random, -PvRratio);
    PhiSc475NegHelCM3->Add(PhiSc475NegHelCM3Random, -PvRratio);
    PhiSc545NegHelCM3->Add(PhiSc545NegHelCM3Random, -PvRratio);
    PhiSc615NegHelCM3->Add(PhiSc615NegHelCM3Random, -PvRratio);
    PhiSc685NegHelCM3->Add(PhiSc685NegHelCM3Random, -PvRratio);
    PhiSc265NegHelCM4->Add(PhiSc265NegHelCM4Random, -PvRratio);
    PhiSc335NegHelCM4->Add(PhiSc335NegHelCM4Random, -PvRratio);
    PhiSc405NegHelCM4->Add(PhiSc405NegHelCM4Random, -PvRratio);
    PhiSc475NegHelCM4->Add(PhiSc475NegHelCM4Random, -PvRratio);
    PhiSc545NegHelCM4->Add(PhiSc545NegHelCM4Random, -PvRratio);
    PhiSc615NegHelCM4->Add(PhiSc615NegHelCM4Random, -PvRratio);
    PhiSc685NegHelCM4->Add(PhiSc685NegHelCM4Random, -PvRratio);
    PhiSc265NegHelCM5->Add(PhiSc265NegHelCM5Random, -PvRratio);
    PhiSc335NegHelCM5->Add(PhiSc335NegHelCM5Random, -PvRratio);
    PhiSc405NegHelCM5->Add(PhiSc405NegHelCM5Random, -PvRratio);
    PhiSc475NegHelCM5->Add(PhiSc475NegHelCM5Random, -PvRratio);
    PhiSc545NegHelCM5->Add(PhiSc545NegHelCM5Random, -PvRratio);
    PhiSc615NegHelCM5->Add(PhiSc615NegHelCM5Random, -PvRratio);
    PhiSc685NegHelCM5->Add(PhiSc685NegHelCM5Random, -PvRratio);
    PhiSc265NegHelCM6->Add(PhiSc265NegHelCM6Random, -PvRratio);
    PhiSc335NegHelCM6->Add(PhiSc335NegHelCM6Random, -PvRratio);
    PhiSc405NegHelCM6->Add(PhiSc405NegHelCM6Random, -PvRratio);
    PhiSc475NegHelCM6->Add(PhiSc475NegHelCM6Random, -PvRratio);
    PhiSc545NegHelCM6->Add(PhiSc545NegHelCM6Random, -PvRratio);
    PhiSc615NegHelCM6->Add(PhiSc615NegHelCM6Random, -PvRratio);
    PhiSc685NegHelCM6->Add(PhiSc685NegHelCM6Random, -PvRratio);
    PhiSc265NegHelCM7->Add(PhiSc265NegHelCM7Random, -PvRratio);
    PhiSc335NegHelCM7->Add(PhiSc335NegHelCM7Random, -PvRratio);
    PhiSc405NegHelCM7->Add(PhiSc405NegHelCM7Random, -PvRratio);
    PhiSc475NegHelCM7->Add(PhiSc475NegHelCM7Random, -PvRratio);
    PhiSc545NegHelCM7->Add(PhiSc545NegHelCM7Random, -PvRratio);
    PhiSc615NegHelCM7->Add(PhiSc615NegHelCM7Random, -PvRratio);
    PhiSc685NegHelCM7->Add(PhiSc685NegHelCM7Random, -PvRratio);
    PhiSc265NegHelCM8->Add(PhiSc265NegHelCM8Random, -PvRratio);
    PhiSc335NegHelCM8->Add(PhiSc335NegHelCM8Random, -PvRratio);
    PhiSc405NegHelCM8->Add(PhiSc405NegHelCM8Random, -PvRratio);
    PhiSc475NegHelCM8->Add(PhiSc475NegHelCM8Random, -PvRratio);
    PhiSc545NegHelCM8->Add(PhiSc545NegHelCM8Random, -PvRratio);
    PhiSc615NegHelCM8->Add(PhiSc615NegHelCM8Random, -PvRratio);
    PhiSc685NegHelCM8->Add(PhiSc685NegHelCM8Random, -PvRratio);
    PhiSc265PosHelCM1->Add(PhiSc265PosHelCM1Random, -PvRratio);
    PhiSc335PosHelCM1->Add(PhiSc335PosHelCM1Random, -PvRratio);
    PhiSc405PosHelCM1->Add(PhiSc405PosHelCM1Random, -PvRratio);
    PhiSc475PosHelCM1->Add(PhiSc475PosHelCM1Random, -PvRratio);
    PhiSc545PosHelCM1->Add(PhiSc545PosHelCM1Random, -PvRratio);
    PhiSc615PosHelCM1->Add(PhiSc615PosHelCM1Random, -PvRratio);
    PhiSc685PosHelCM1->Add(PhiSc685PosHelCM1Random, -PvRratio);
    PhiSc265PosHelCM2->Add(PhiSc265PosHelCM2Random, -PvRratio);
    PhiSc335PosHelCM2->Add(PhiSc335PosHelCM2Random, -PvRratio);
    PhiSc405PosHelCM2->Add(PhiSc405PosHelCM2Random, -PvRratio);
    PhiSc475PosHelCM2->Add(PhiSc475PosHelCM2Random, -PvRratio);
    PhiSc545PosHelCM2->Add(PhiSc545PosHelCM2Random, -PvRratio);
    PhiSc615PosHelCM2->Add(PhiSc615PosHelCM2Random, -PvRratio);
    PhiSc685PosHelCM2->Add(PhiSc685PosHelCM2Random, -PvRratio);
    PhiSc265PosHelCM3->Add(PhiSc265PosHelCM3Random, -PvRratio);
    PhiSc335PosHelCM3->Add(PhiSc335PosHelCM3Random, -PvRratio);
    PhiSc405PosHelCM3->Add(PhiSc405PosHelCM3Random, -PvRratio);
    PhiSc475PosHelCM3->Add(PhiSc475PosHelCM3Random, -PvRratio);
    PhiSc545PosHelCM3->Add(PhiSc545PosHelCM3Random, -PvRratio);
    PhiSc615PosHelCM3->Add(PhiSc615PosHelCM3Random, -PvRratio);
    PhiSc685PosHelCM3->Add(PhiSc685PosHelCM3Random, -PvRratio);
    PhiSc265PosHelCM4->Add(PhiSc265PosHelCM4Random, -PvRratio);
    PhiSc335PosHelCM4->Add(PhiSc335PosHelCM4Random, -PvRratio);
    PhiSc405PosHelCM4->Add(PhiSc405PosHelCM4Random, -PvRratio);
    PhiSc475PosHelCM4->Add(PhiSc475PosHelCM4Random, -PvRratio);
    PhiSc545PosHelCM4->Add(PhiSc545PosHelCM4Random, -PvRratio);
    PhiSc615PosHelCM4->Add(PhiSc615PosHelCM4Random, -PvRratio);
    PhiSc685PosHelCM4->Add(PhiSc685PosHelCM4Random, -PvRratio);
    PhiSc265PosHelCM5->Add(PhiSc265PosHelCM5Random, -PvRratio);
    PhiSc335PosHelCM5->Add(PhiSc335PosHelCM5Random, -PvRratio);
    PhiSc405PosHelCM5->Add(PhiSc405PosHelCM5Random, -PvRratio);
    PhiSc475PosHelCM5->Add(PhiSc475PosHelCM5Random, -PvRratio);
    PhiSc545PosHelCM5->Add(PhiSc545PosHelCM5Random, -PvRratio);
    PhiSc615PosHelCM5->Add(PhiSc615PosHelCM5Random, -PvRratio);
    PhiSc685PosHelCM5->Add(PhiSc685PosHelCM5Random, -PvRratio);
    PhiSc265PosHelCM6->Add(PhiSc265PosHelCM6Random, -PvRratio);
    PhiSc335PosHelCM6->Add(PhiSc335PosHelCM6Random, -PvRratio);
    PhiSc405PosHelCM6->Add(PhiSc405PosHelCM6Random, -PvRratio);
    PhiSc475PosHelCM6->Add(PhiSc475PosHelCM6Random, -PvRratio);
    PhiSc545PosHelCM6->Add(PhiSc545PosHelCM6Random, -PvRratio);
    PhiSc615PosHelCM6->Add(PhiSc615PosHelCM6Random, -PvRratio);
    PhiSc685PosHelCM6->Add(PhiSc685PosHelCM6Random, -PvRratio);
    PhiSc265PosHelCM7->Add(PhiSc265PosHelCM7Random, -PvRratio);
    PhiSc335PosHelCM7->Add(PhiSc335PosHelCM7Random, -PvRratio);
    PhiSc405PosHelCM7->Add(PhiSc405PosHelCM7Random, -PvRratio);
    PhiSc475PosHelCM7->Add(PhiSc475PosHelCM7Random, -PvRratio);
    PhiSc545PosHelCM7->Add(PhiSc545PosHelCM7Random, -PvRratio);
    PhiSc615PosHelCM7->Add(PhiSc615PosHelCM7Random, -PvRratio);
    PhiSc685PosHelCM7->Add(PhiSc685PosHelCM7Random, -PvRratio);
    PhiSc265PosHelCM8->Add(PhiSc265PosHelCM8Random, -PvRratio);
    PhiSc335PosHelCM8->Add(PhiSc335PosHelCM8Random, -PvRratio);
    PhiSc405PosHelCM8->Add(PhiSc405PosHelCM8Random, -PvRratio);
    PhiSc475PosHelCM8->Add(PhiSc475PosHelCM8Random, -PvRratio);
    PhiSc545PosHelCM8->Add(PhiSc545PosHelCM8Random, -PvRratio);
    PhiSc615PosHelCM8->Add(PhiSc615PosHelCM8Random, -PvRratio);
    PhiSc685PosHelCM8->Add(PhiSc685PosHelCM8Random, -PvRratio);

}


Bool_t	PNeutPol_Polarimeter_Circ::Write(){
    GTreeManager::Write();
}
