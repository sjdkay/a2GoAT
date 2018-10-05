// GoAT Physics analysis to identify neutrons from deuterium photodisintegration
// Various properties of neutrons/protons identified plotted in histograms
// Main aim is to determine spin polarisation of neutrons
// For use on circularly polarised data OR Linearly Polarised files if to be combined files

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
    CosThetaWidth = (0.4);
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

    //Cut_CB_proton = OpenCutFile("configfiles/cuts/CB_DeltaE-E_Proton_24_4_18.root", "Proton");
    Cut_CB_proton = OpenCutFile("configfiles/cuts/CB_DeltaE-E_Proton_25_5_18_V2.root", "Proton");
    Cut_proton = Cut_CB_proton;
    Cut_CB_pion = OpenCutFile("configfiles/cuts/CB_DeltaE-E_Pion_29_07_15.root", "Pion");
    Cut_pion = Cut_CB_pion;
    Cut_CB_protonKinGood = OpenCutFile("configfiles/cuts/CB_DeltaE-E_ProtonKinGood_11_05_17.root", "ProtonKinGood");
    Cut_protonKinGood = Cut_CB_protonKinGood;
    Cut_CB_protonKinBad = OpenCutFile("configfiles/cuts/CB_DeltaE-E_ProtonKinBad_15_12_16.root", "ProtonKinBad");
    Cut_protonKinBad = Cut_CB_protonKinBad;

    fAy = new TFile("/home/sjdkay/work/A2/a2GoAT/npAy.root");

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
    }

    // If track 2 only gives signals in MWPC and CB it is the neutron
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

    EpDiff = abs(EpCorr - Ep);

    Pp = sqrt (TMath::Power((Mp + EpCorr),2) - TMath::Power(Mp,2));
    Pn = sqrt (TMath::Power((En + Mn ),2) - TMath::Power(Mn,2));
    GVpCorr = TLorentzVector(Pp*sin(ThpRad)*cos(PhpRad), Pp*sin(ThpRad)*sin(PhpRad), Pp*cos(ThpRad), EpCorr+Mp);
    MWPC0pEVeto = MWPC0pE*sin(GVpCorr.Theta());
    MWPC1pEVeto = MWPC1pE*sin(GVpCorr.Theta());
    MWPC0nEVeto = MWPC0pE*sin(GVpCorr.Theta());
    MWPC1nEVeto = MWPC1pE*sin(GVpCorr.Theta());

    if((MWPC0pE != 0) && (MWPC1pE != 0)){ // If there is a non zero energy deposit in each MWPC for proton set one MWPC energy sum
        MWPCpETot2 = MWPC0pE + MWPC1pE;
        MWPCpEVetoTot2 = MWPC0pEVeto + MWPC1pEVeto;
        MWPCpEVetoTot1 = 0;
    }
    else if (((MWPC0pE != 0) && (MWPC1pE == 0)) == kTRUE || ((MWPC0pE == 0) && (MWPC1pE != 0)) == kTRUE){ // If only detect proton in one chamber set different ESum
        MWPCpETot1 = MWPC0pE + MWPC1pE;
        MWPCpEVetoTot1 = MWPC0pEVeto + MWPC1pEVeto;
        MWPCpEVetoTot2 = 0;
    }
    if((MWPC0nE != 0) && (MWPC1nE != 0)){ // If there is a non zero energy deposit in each MWPC for neutron set one MWPC energy sum
        MWPCnETot2 = MWPC0nE + MWPC1nE;
        MWPCnEVetoTot2 = MWPC0nEVeto + MWPC1nEVeto;
        MWPCnEVetoTot1 = 0;
    }
    else if (((MWPC0nE != 0) && (MWPC1nE == 0)) == kTRUE || ((MWPC0nE == 0) && (MWPC1nE != 0)) == kTRUE){ // If only detect neutron in one chamber set different ESum
        MWPCnETot1 = MWPC0nE + MWPC1nE;
        MWPCnEVetoTot1 = MWPC0nEVeto + MWPC1nEVeto;
        MWPCnEVetoTot2 = 0;
    }

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
        pSc = GVnCorr;
        pScB = pSc; //Boost the scattered proton track to CM frame
        pScB.Boost(b);
        CosThetapScCM = cos (pScB.Theta());

        RecNeutronEpCorr = CNeutron4VectorKin(GVpCorr);
        MMpEpCorr = RecNeutronEpCorr.M();

        if (((MMpEpCorr < 800) == kTRUE) || ((MMpEpCorr > 1300) == kTRUE)) continue; // Force a missing mass cut

        // Gamma(d,p)n calc n from kinematics
        KinEp = CalcKinEnergy(Thp, EGamma, Md, 0., Mp, Mn); // Calculate kin E of proton assuming pn production
        pKin = CProton4VectorKin(KinEp, ThpRad, PhpRad);
        nKin = CNeutron4VectorKin(pKin);
        ThetanRec = (nKin.Theta()) * TMath::RadToDeg();
        PhinRec = (nKin.Phi()) * TMath::RadToDeg();
        pKin3 = pKin.Vect();
        nKin3 = nKin.Vect();
        nKinE = nKin.E() - Mn;

        nKinB = nKin;
        nKinB.Boost(b);
        ThetanCM = (nKinB.Theta())*TMath::RadToDeg();
        //if (80 > ThetanCM || ThetanCM > 100) continue;

        KinEDiff = KinEp - EpCorr;

        ThetanDiff = ThetanCorr-ThetanRec;
        PhinDiff = PhinCorr-PhinRec;

        TVector3 ScattAngles = ScatteredFrameAngles(nKin3, pKinB.Vect(), GVn3 , Gamma);
        ScattThetaRad = ScattAngles(2);
        ScattTheta = ScattThetaRad*TMath::RadToDeg();
        ScattPhiRad = ScattAngles(0);
        ScattPhi = ScattPhiRad*TMath::RadToDeg();
        //if (ScattPhiRad > -0.01 && ScattPhiRad < 0.01) continue; // Kill excessively small values
        //if (abs(ScattPhiRad) > 3.13) continue;
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
        if((MWPC0nE != 0) && (MWPC1nE != 0)){
            if(MWPCnEVetoTot2 < 89.326) return; // If Total MWPC veto energy for scattered proton track less than X return
        }
        else if (((MWPC0nE != 0) && (MWPC1nE == 0)) == kTRUE || ((MWPC0nE == 0) && (MWPC1nE != 0)) == kTRUE){ // If only detect neutron in one chamber set different ESum
            if(MWPCnEVetoTot1 < 76.7) return;
        }

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

        if( ((nKinE < 1050) == kTRUE) && ((ScattTheta < 90) == kTRUE) ){
            AEff = ((TH2F*)fAy->Get("nppnAy"))->Interpolate(nKinE, ScattTheta);
        }

        else if (((nKinE > 1050) == kTRUE) || ((ScattTheta > 90) == kTRUE)) return;

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

    ThetaSc = new TH1D ("ThetaSc", "#theta_{Sc} Distribution", 200, 0, acos(-1));
    ThetaScPrompt = new TH1D ("ThetaScPrompt", "#theta_{Sc} Distribution", 200, 0, acos(-1));
    ThetaScRandom = new TH1D ("ThetaScRandom", "#theta_{Sc} Distribution", 200, 0, acos(-1));

    ZpDist = new GH1 ("ZpDist", "Proton Pseudo Z Vertex Distribution", 200, -400, 400);

    RPoca = new TH1D("RPoca", "R_{POCA} Distribution", 200, 0, 200);
    RPocaPrompt = new TH1D("RPocaPrompt", "R_{POCA} Distribution", 200, 0, 200);
    RPocaRandom = new TH1D("RPocaRandom", "R_{POCA} Distribution", 200, 0, 200);

    Thetap = new TH1D ("Thetap", "#theta_{p} Distribution", 200, 0, acos(-1));
    ThetapPrompt = new TH1D ("ThetapPrompt", "#theta_{p} Distribution", 200, 0, acos(-1));
    ThetapRandom = new TH1D ("ThetapRandom", "#theta_{p} Distribution", 200, 0, acos(-1));

    Phip = new TH1D( "Phip", "#phi_{p} Proton Distribution", 200, -1*acos(-1), acos(-1));
    PhipPrompt = new TH1D( "PhipPrompt", "#phi_{p} Proton Distribution", 200, -1*acos(-1), acos(-1));
    PhipRandom = new TH1D( "PhipRandom", "#phi_{p} Proton Distribution", 200, -1*acos(-1), acos(-1));

    EProton = new TH1D( "EProton", "E_{p} Distribution", 200, 0, 400);
    EProtonPrompt = new TH1D( "EProtonPrompt", "E_{p} Distribution", 200, 0, 400);
    EProtonRandom = new TH1D( "EProtonRandom", "E_{p} Distribution", 200, 0, 400);

    MMProton = new TH1D( "MMProton", "MM_{p} Distribution", 200, 0, 2000);
    MMProtonPrompt = new TH1D( "MMProtonPrompt", "MM_{p} Distribution", 200, 0, 2000);
    MMProtonRandom = new TH1D( "MMProtonRandom", "MM_{p} Distribution", 200, 0, 2000);

    EdE = new TH2D ("EdE", "dE_{p} as a function of E_{p}", 200, 0, 1000, 200, 0, 10);
    EdEPrompt = new TH2D ("EdEPrompt", "dE_{p} as a function of E_{p}", 200, 0, 1000, 200, 0, 10);
    EdERandom = new TH2D ("EdERandom", "dE_{p} as a function of E_{p}", 200, 0, 1000, 200, 0, 10);

    ThetaScThetap = new TH2D ("ThetaScThetap", "#theta_{p}(#theta_{Sc})", 200, 0, acos(-1), 200, 0, acos(-1));
    ThetaScThetapPrompt = new TH2D ("ThetaScThetapPrompt", "#theta_{p}(#theta_{Sc})", 200, 0, acos(-1), 200, 0, acos(-1));
    ThetaScThetapRandom = new TH2D ("ThetaScThetapRandom", "#theta_{p}(#theta_{Sc})", 200, 0, acos(-1), 200, 0, acos(-1));

    PhiScThetap = new TH2D ("PhiScThetap", "#theta_{p}(#phi_{Sc})", 200, -1*acos(-1), acos(-1), 200,  0, acos(-1));
    PhiScThetapPrompt = new TH2D ("PhiScThetapPrompt", "#theta_{p}(#phi_{Sc})", 200, -1*acos(-1), acos(-1), 200, 0, acos(-1));
    PhiScThetapRandom = new TH2D ("PhiScThetapRandom", "#theta_{p}(#phi_{Sc})", 200, -1*acos(-1), acos(-1), 200,  0, acos(-1));

    rPocaThetap = new TH2D ("rPocaThetap", "#theta_{p}(r_{POCA})", 200, 0, 200, 200,  0, acos(-1));
    rPocaThetapPrompt = new TH2D ("rPocaThetapPrompt", "#theta_{p}(r_{POCA})", 200, 0, 200, 200, 0, acos(-1));
    rPocaThetapRandom = new TH2D ("rPocaThetapRandom", "#theta_{p}(r_{POCA})", 200, 0, 200, 200, 0, acos(-1));

    ThetaDiffPhiDiff = new TH2D("ThetaDiffPhiDiff", "#phi_{nDiff}(#theta_{nDiff})", 200, -180, 180, 200, -180, 180);
    ThetaDiffPhiDiffPrompt = new TH2D("ThetaDiffPhiDiffPrompt", "#phi_{nDiff}(#theta_{nDiff})", 200, -180, 180, 200, -180, 180);
    ThetaDiffPhiDiffRandom = new TH2D("ThetaDiffPhiDiffRandom", "#phi_{nDiff}(#theta_{nDiff})", 200, -180, 180, 200, -180, 180);

    ThetapEgamma = new TH2D("ThetapEg", "#theta_{p}(E_{#gamma})", 200, 0, 1600, 200, 0, acos(-1));
    ThetapEgammaPrompt = new TH2D("ThetapEgPrompt", "#theta_{p}(E_{#gamma})", 200, 0, 1600, 200, 0, acos(-1));
    ThetapEgammaRandom = new TH2D("ThetapEgRandom", "#theta_{p}(E_{#gamma})", 200, 0, 1600, 200, 0, acos(-1));

    ThetapMMp = new TH2D("ThetapMMp", "#theta_{p}(MM_{p})", 200, 800, 1200, 200, 0, acos(-1));
    ThetapMMpPrompt = new TH2D("ThetapMMpPrompt", "#theta_{p}(MM_{p})", 200, 800, 1200, 200, 0, acos(-1));
    ThetapMMpRandom = new TH2D("ThetapMMpRandom", "#theta_{p}(MM_{p})", 200, 800, 1200, 200, 0, acos(-1));

    ThetanEgamma = new TH2D("ThetanEg", "#theta_{n}(E_{#gamma})", 200, 0, 1600, 200, 0, acos(-1));
    ThetanEgammaPrompt = new TH2D("ThetanEgPrompt", "#theta_{n}(E_{#gamma})", 200, 0, 1600, 200, 0, acos(-1));
    ThetanEgammaRandom = new TH2D("ThetanEgRandom", "#theta_{n}(E_{#gamma})", 200, 0, 1600, 200, 0, acos(-1));

    ThetanMMp = new TH2D("ThetanMMp", "#theta_{n}(MM_{p})", 200, 800, 1200, 200, 0, acos(-1));
    ThetanMMpPrompt = new TH2D("ThetanMMpPrompt", "#theta_{n}(MM_{p})", 200, 800, 1200, 200, 0, acos(-1));
    ThetanMMpRandom = new TH2D("ThetanMMpRandom", "#theta_{n}(MM_{p})", 200, 800, 1200, 200, 0, acos(-1));

    PhiScThetan = new TH2D ("PhiScThetan", "#theta_{n}(#phi_{Sc})", 200, -1*acos(-1), acos(-1), 200,  0, acos(-1));
    PhiScThetanPrompt = new TH2D ("PhiScThetanPrompt", "#theta_{n}(#phi_{Sc})", 200, -1*acos(-1), acos(-1), 200, 0, acos(-1));
    PhiScThetanRandom = new TH2D ("PhiScThetanRandom", "#theta_{n}(#phi_{Sc})", 200, -1*acos(-1), acos(-1), 200,  0, acos(-1));

    PhiScThetapCM = new TH2D ("PhiScThetapCM", "#theta_{pCM}(#phi_{Sc})", 200, -1*acos(-1), acos(-1), 200,  0, acos(-1));
    PhiScThetapCMPrompt = new TH2D ("PhiScThetapCMPrompt", "#theta_{pCM}(#phi_{Sc})", 200, -1*acos(-1), acos(-1), 200, 0, acos(-1));
    PhiScThetapCMRandom = new TH2D ("PhiScThetapCMRandom", "#theta_{pCM}(#phi_{Sc})", 200, -1* acos(-1), acos(-1), 200,  0, acos(-1));

    WeightEg = new TH2D("WeightEg", "Weight(E_{#gamma})", 200, 0, 1600, 200, -1, 1);
    WeightEgPrompt = new TH2D("WeightEgPrompt", "Weight(E_{#gamma})", 200, 0, 1600, 200, -1, 1);
    WeightEgRandom = new TH2D("WeightEgRandom", "Weight(E_{#gamma})", 200, 0, 1600, 200, -1, 1);

    WeightPhiSc = new TH2D("WeightPhiSc", "Weight(#phi_{Sc})", 200, -1*acos(-1), acos(-1), 200, -1, 1);
    WeightPhiScPrompt = new TH2D("WeightPhiScPrompt", "Weight(#phi_{Sc})", 200, -1*acos(-1), acos(-1), 200, -1, 1);
    WeightPhiScRandom = new TH2D("WeightPhiScRandom", "Weight(#phi_{Sc})", 200, -1*acos(-1), acos(-1), 200, -1, 1);

    ThetapThetan = new TH2D ("ThetapThetan", "#theta_{p}(#theta_{n})", 200, 0, acos(-1), 200, 0, acos(-1));
    ThetapThetanPrompt = new TH2D ("ThetapThetanPrompt", "#theta_{p}(#theta_{n})", 200, 0, acos(-1), 200, 0, acos(-1));
    ThetapThetanRandom = new TH2D ("ThetapThetanRandom", "#theta_{p}(#theta_{n})", 200, 0, acos(-1), 200, 0, acos(-1));

    NeutronE = new TH1D("NeutronE", "E_{n}", 400, 0, 500);
    NeutronEPrompt = new TH1D("NeutronEPrompt", "E_{n}", 400, 0, 500);
    NeutronERandom = new TH1D("NeutronERandom", "E_{n}", 400, 0, 500);

    NeutronEThetaScFull = new TH2D("NeutronEThetaScFull", "#theta_{Sc}(E_{n})", 200, 0, 500, 200, 0, 90);
    NeutronEThetaScFullPrompt = new TH2D("NeutronEThetaScFullPrompt", "#theta_{Sc}(E_{n})", 200, 0, 500, 200, 0, 90);
    NeutronEThetaScFullRandom = new TH2D("NeutronEThetaScFullRandom", "#theta_{Sc}(E_{n})", 200, 0, 500, 200, 0, 90);

    NeutronEEg = new TH2D("NeutronEEg", "E_{#gamma}(E_{n})", 200, 0, 500, 200, 0, 1600);
    NeutronEEgPrompt = new TH2D("NeutronEEgPrompt", "E_{#gamma}(E_{n})", 200, 0, 500, 200, 0, 1600);
    NeutronEEgRandom = new TH2D("NeutronEEgRandom", "E_{#gamma}(E_{n})", 200, 0, 500, 200, 0, 1600);

    ThetaScEg = new TH2D("ThetaScEg", "E_{#gamma}(#theta_{Sc})", 200, 0, 90, 200, 0, 1600);
    ThetaScEgPrompt = new TH2D("ThetaScEgPrompt", "E_{#gamma}(#theta_{Sc})", 200, 0, 90, 200, 0, 1600);
    ThetaScEgRandom = new TH2D("ThetaScEgRandom", "E_{#gamma}(#theta_{Sc})", 200, 0, 90, 200, 0, 1600);

    EpMWPCEpVetoTot1 = new TH2D("EpMWPCEpVetoTot1", "EVeto^{Tot1}_MWPCp(E_{p})", 200, 0, 500, 200, 0, 2000);
    EpMWPCEpVetoTot1Prompt = new TH2D("EpMWPCEpVetoTot1Prompt", "EVeto^{Tot1}_MWPCp(E_{p})", 200, 0, 500, 200, 0, 2000);
    EpMWPCEpVetoTot1Random = new TH2D("EpMWPCEpVetoTot1Random", "EVeto^{Tot1}_MWPCp(E_{p})", 200, 0, 500, 200, 0, 2000);

    EnMWPCEnVetoTot1 = new TH2D("EnMWPCEnVetoTot1", "EVeto^{Tot1}_MWPCn(E_{n})", 200, 0, 500, 200, 0, 400);
    EnMWPCEnVetoTot1Prompt = new TH2D("EnMWPCEnVetoTot1Prompt", "E^{Tot1}_MWPCn(E_{n})", 200, 0, 500, 200, 0, 400);
    EnMWPCEnVetoTot1Random = new TH2D("EnMWPCEnVetoTot1Random", "E^{Tot1}_MWPCn(E_{n})", 200, 0, 500, 200, 0, 400);

    EpMWPCEpVetoTot2 = new TH2D("EpMWPCEpVetoTot2", "EVeto^{Tot1}_MWPCp(E_{p})", 200, 0, 500, 200, 0, 2000);
    EpMWPCEpVetoTot2Prompt = new TH2D("EpMWPCEpVetoTot2Prompt", "EVeto^{Tot2}_MWPCp(E_{p})", 200, 0, 500, 200, 0, 2000);
    EpMWPCEpVetoTot2Random = new TH2D("EpMWPCEpVetoTot2Random", "EVeto^{Tot2}_MWPCp(E_{p})", 200, 0, 500, 200, 0, 2000);

    EnMWPCEnVetoTot2 = new TH2D("EnMWPCEnVetoTot2", "EVeto^{Tot2}_MWPCn(E_{n})", 200, 0, 500, 200, 0, 1600);
    EnMWPCEnVetoTot2Prompt = new TH2D("EnMWPCEnVetoTot2Prompt", "E^{Tot2}_MWPCn(E_{n})", 200, 0, 500, 200, 0, 1600);
    EnMWPCEnVetoTot2Random = new TH2D("EnMWPCEnVetoTot2Random", "E^{Tot2}_MWPCn(E_{n})", 200, 0, 500, 200, 0, 1600);

    ThetapMWPCEpTot1 = new TH2D("ThetapMWPCEpTot1", "E^{Tot1}_MWPCp(#theta_{p})", 200, 0, 180, 200, 0, 400);
    ThetapMWPCEpTot1Prompt = new TH2D("ThetapMWPCEpTot1Prompt", "E^{Tot1}_MWPCp(#theta_{p})", 200, 0, 180, 200, 0, 400);
    ThetapMWPCEpTot1Random = new TH2D("ThetapMWPCEpTot1Random", "E^{Tot1}_MWPCp(#theta_{p})", 200, 0, 180, 200, 0, 400);

    ThetanMWPCEnTot1 = new TH2D("ThetanMWPCEnTot1", "E^{Tot1}_MWPCp(#theta_{n})", 200, 0, 180, 200, 0, 400);
    ThetanMWPCEnTot1Prompt = new TH2D("ThetanMWPCEnTot1Prompt", "E^{Tot1}_MWPCn(#theta_{n})", 200, 0, 180, 200, 0, 400);
    ThetanMWPCEnTot1Random = new TH2D("ThetanMWPCEnTot1Random", "E^{Tot1}_MWPCn(#theta_{n})", 200, 0, 180, 200, 0, 400);

    ThetapMWPCEpTot2 = new TH2D("ThetapMWPCEpTot2", "E^{Tot2}_MWPCp(#theta_{p})", 200, 0, 180, 200, 0, 1600);
    ThetapMWPCEpTot2Prompt = new TH2D("ThetapMWPCEpTot2Prompt", "E^{Tot2}_MWPCp(#theta_{p})", 200, 0, 180, 200, 0, 1600);
    ThetapMWPCEpTot2Random = new TH2D("ThetapMWPCEpTot2Random", "E^{Tot2}_MWPCp(#theta_{p})", 200, 0, 180, 200, 0, 1600);

    ThetanMWPCEnTot2 = new TH2D("ThetanMWPCEnTot2", "E^{Tot2}_MWPCp(#theta_{n})", 200, 0, 180, 200, 0, 1600);
    ThetanMWPCEnTot2Prompt = new TH2D("ThetanMWPCEnTot2Prompt", "E^{Tot2}_MWPCn(#theta_{n})", 200, 0, 180, 200, 0, 1600);
    ThetanMWPCEnTot2Random = new TH2D("ThetanMWPCEnTot2Random", "E^{Tot2}_MWPCn(#theta_{n})", 200, 0, 180, 200, 0, 1600);

    DOCADist = new TH1D("DOCADist", "DOCA", 500, 0, 1000);
    DOCADistPrompt = new TH1D("DOCADistPrompt", "DOCA", 500, 0, 1000);
    DOCADistRandom = new TH1D("DOCADistRandom", "DOCA", 500, 0, 1000);

    DOCArPOCA = new TH2D("DOCArPOCA", "DOCA(r_{POCA})", 200, 0, 200, 200, 0, 1000);
    DOCArPOCAPrompt = new TH2D("DOCArPOCAPrompt", "DOCA(r_{POCA})", 200, 0, 200, 200, 0, 1000);
    DOCArPOCARandom = new TH2D("DOCArPOCARandom", "DOCA(r_{POCA})", 200, 0, 200, 200, 0, 1000);

    ThetaScPhiSc = new TH2D("ThetaScPhiSc", "#phi_{Sc}(#theta_{Sc})", 200, 0 , 1, 200, -1*acos(-1) , acos(-1));
    ThetaScPhiScPrompt = new TH2D("ThetaScPhiScPrompt", "#phi_{Sc}(#theta_{Sc})", 200, 0 , 1, 200, -1*acos(-1) , acos(-1));
    ThetaScPhiScRandom = new TH2D("ThetaScPhiScRandom", "#phi_{Sc}(#theta_{Sc})", 200, 0 , 1, 200, -1*acos(-1) , acos(-1));

    PhiScPosHelFull = new TH1D("PhiScPosHelFull", "#phi_{Sc} +ve Helicity", 200, -1*acos(-1), acos(-1));
    PhiScPosHelFullPrompt = new TH1D("PhiScPosHelFullPrompt", "#phi_{Sc} +ve Helicity", 200, -1*acos(-1), acos(-1));
    PhiScPosHelFullRandom = new TH1D("PhiScPosHelFullRandom", "#phi_{Sc} +ve Helicity", 200, -1*acos(-1), acos(-1));

    PhiScNegHelFull = new TH1D("PhiScNegHelFull", "#phi_{Sc} +ve Helicity", 200, -1*acos(-1), acos(-1));
    PhiScNegHelFullPrompt = new TH1D("PhiScNegHelFullPrompt", "#phi_{Sc} +ve Helicity", 200, -1*acos(-1), acos(-1));
    PhiScNegHelFullRandom = new TH1D("PhiScNegHelFullRandom", "#phi_{Sc} +ve Helicity", 200, -1*acos(-1), acos(-1));

    MMEg = new TH2D("MMEg", "MM_{p}(E_{#gamma})", 200, 0, 1000, 200, 0, 2000);
    MMEgPrompt = new TH2D("MMEgPrompt", "MM_{p}(E_{#gamma})", 200, 0, 1000, 200, 0, 2000);
    MMEgRandom = new TH2D("MMEgRadnom", "MM_{p}(E_{#gamma})", 200, 0, 1000, 200, 0, 2000);

    PhiScEg = new TH2D ("PhiScEg", "#phi_{Sc}(E_{#gamma})", 200, 0 ,1000, 200, -1*acos(-1), acos(-1));
    PhiScEgPrompt = new TH2D ("PhiScEgPrompt", "#phi_{Sc}(E_{#gamma})", 200, 0, 1000, 200, -1*acos(-1), acos(-1));
    PhiScEgRandom = new TH2D ("PhiScEgRandom", "#phi_{Sc}(E_{#gamma})", 200, 0, 1000, 200, -1*acos(-1), acos(-1));

    PhiScZp = new TH2D ("PhiScZp", "#phi_{Sc}(Z_{p})", 200, -100 ,100, 200, -1*acos(-1), acos(-1));
    PhiScZpPrompt = new TH2D ("PhiScZpPrompt", "#phi_{Sc}(Z_{p})", 200, -100, 100, 200, -1*acos(-1), acos(-1));
    PhiScZpRandom = new TH2D ("PhiScZpRandom", "#phi_{Sc}(Z_{p})", 200, -100, 100, 200, -1*acos(-1), acos(-1));

    for(Int_t A = 0; A < 8; A++){
        for(Int_t B = 0; B < 5; B++){
            PhiScPosHel[A][B] = new TH1D(Form("PhiSc%iPosHelCM%i", 340+(A*100), B+1), Form("#phi_{Sc} E_{#gamma}%i #pm 50MeV CM%i for +ve Helicity", 340+(A*100), B+1), 10, -1*acos(-1) , acos(-1));
            PhiScPosHelPrompt[A][B] = new TH1D(Form("PhiSc%iPosHelCM%iPrompt", 340+(A*100), B+1), Form("#phi_{Sc} E_{#gamma}%i #pm 50MeV CM%i for +ve Helicity", 340+(A*100), B+1), 10, -1*acos(-1) , acos(-1));
            PhiScPosHelRandom[A][B] = new TH1D(Form("PhiSc%iPosHelCM%iRandom", 340+(A*100), B+1), Form("#phi_{Sc} E_{#gamma}%i #pm 50MeV CM%i for +ve Helicity", 340+(A*100), B+1), 10, -1*acos(-1) , acos(-1));
            PhiScNegHel[A][B] = new TH1D(Form("PhiSc%iNegHelCM%i", 340+(A*100), B+1), Form("#phi_{Sc} E_{#gamma}%i #pm 50MeV CM%i for -ve Helicity", 340+(A*100), B+1), 10, -1*acos(-1) , acos(-1));
            PhiScNegHelPrompt[A][B] = new TH1D(Form("PhiSc%iNegHelCM%iPrompt", 340+(A*100), B+1), Form("#phi_{Sc} E_{#gamma}%i #pm 50MeV CM%i for -ve Helicity", 340+(A*100), B+1), 10, -1*acos(-1) , acos(-1));
            PhiScNegHelRandom[A][B] = new TH1D(Form("PhiSc%iNegHelCM%iRandom", 340+(A*100), B+1), Form("#phi_{Sc} E_{#gamma}%i #pm 50MeV CM%i for -ve Helicity", 340+(A*100), B+1), 10, -1*acos(-1) , acos(-1));
            PhiScAEffPosHel[A][B] = new TH2D(Form("PhiScAEff%iPosHelCM%i", 340+(A*100), B+1), Form("#phi_{Sc} E_{#gamma}%i #pm 50MeV CM%i for +ve Helicity", 340+(A*100), B+1), 10, -1*acos(-1), acos(-1), 500, -1, 1);
            PhiScAEffPosHelPrompt[A][B] = new TH2D(Form("PhiScAEff%iPosHelCM%iPrompt", 340+(A*100), B+1), Form("#phi_{Sc} E_{#gamma}%i #pm 50MeV CM%i for +ve Helicity", 340+(A*100), B+1), 10, -1*acos(-1), acos(-1), 500, -1, 1);
            PhiScAEffPosHelRandom[A][B] = new TH2D(Form("PhiScAEff%iPosHelCM%iRandom", 340+(A*100), B+1), Form("#phi_{Sc} E_{#gamma}%i #pm 50MeV CM%i for +ve Helicity", 340+(A*100), B+1), 10, -1*acos(-1), acos(-1), 500, -1, 1);
            PhiScAEffNegHel[A][B] = new TH2D(Form("PhiScAEff%iNegHelCM%i", 340+(A*100), B+1), Form("#phi_{Sc} E_{#gamma}%i #pm 50MeV CM%i for -ve Helicity", 340+(A*100), B+1), 10, -1*acos(-1), acos(-1), 500, -1, 1);
            PhiScAEffNegHelPrompt[A][B] = new TH2D(Form("PhiScAEff%iNegHelCM%iPrompt", 340+(A*100), B+1), Form("#phi_{Sc} E_{#gamma}%i #pm 50MeV CM%i for -ve Helicity", 340+(A*100), B+1), 10, -1*acos(-1), acos(-1), 500, -1, 1);
            PhiScAEffNegHelRandom[A][B] = new TH2D(Form("PhiScAEff%iNegHelCM%iRandom", 340+(A*100), B+1), Form("#phi_{Sc} E_{#gamma}%i #pm 50MeV CM%i for -ve Helicity", 340+(A*100), B+1), 10, -1*acos(-1), acos(-1), 500, -1, 1);

//            NeutronEThetaSc[A][B] = new TH2D(Form ("NeutronEThetaSc%iCM%i", 250+(A*100), B+1), Form ("#theta_{Sc}(E_{n}) %iMeV CM%i", 250+(A*100), B+1), 200, 0, 800, 200, 0, 90);
//            NeutronEThetaScPrompt[A][B] = new TH2D(Form ("NeutronEThetaSc%iCM%iPrompt", 250+(A*100), B+1), Form ("#theta_{Sc}(E_{n}) %iMeV CM%i", 250+(A*100), B+1), 200, 0, 800, 200, 0, 90);
//            NeutronEThetaScRandom[A][B] = new TH2D(Form ("NeutronEThetaSc%iCM%iRandom", 250+(A*100), B+1), Form ("#theta_{Sc}(E_{n}) %iMeV CM%i", 250+(A*100), B+1), 200, 0, 800, 200, 0, 90);
        }
    }

    for(Int_t X = 0; X < 12; X++){
        MMpEgamma[X] = new TH1D(Form("MMp_%iMeV", 200+(X*50)), Form("MM_{p} %i #pm 25MeV", 200+(X*50)), 400, 0, 2000);
        MMpEgammaPrompt[X] = new TH1D(Form("MMp_%iMeVPrompt", 200+(X*50)), Form("MM_{p} %i #pm 25MeV", 200+(X*50)), 400, 0, 2000);
        MMpEgammaRandom[X] = new TH1D(Form("MMp_%iMeVRandom", 200+(X*50)), Form("MM_{p} %i #pm 25MeV", 200+(X*50)), 400, 0, 2000);
        PhiSc[X] = new TH1D(Form("PhiSc%i", 200+(X*50)), Form("#phi_{Sc} E_{#gamma}%i #pm 25MeV", 200+(X*50)), 10, -1*acos(-1) , acos(-1));
        PhiScPrompt[X] = new TH1D(Form("PhiSc%iPrompt", 200+(X*50)), Form("#phi_{Sc} E_{#gamma}%i #pm 25MeV", 200+(X*50)), 10, -1*acos(-1) , acos(-1));
        PhiScRandom[X] = new TH1D(Form("PhiSc%iRandom", 200+(X*50)), Form("#phi_{Sc} E_{#gamma}%i #pm 25MeV", 200+(X*50)), 10, -1*acos(-1) , acos(-1));
    }

    for(Int_t C = 0; C < 10; C++){
        for(Int_t D = 0; D < 5; D++){
            PhipSet[C][D] = new TH1D(Form("Phip_%iMeVCM%i", 415+(C*20), D+1), Form("#phi_{p} %i #pm 10MeV CM%i", 415+(C*20), D+1), 10, -180, 180);
            PhipSetPrompt[C][D] = new TH1D(Form("Phip_%iMeVCM%iPrompt", 415+(C*20), D+1), Form("#phi_{p} %i #pm 10MeV CM%i", 415+(C*20), D+1), 10, -180, 180);
            PhipSetRandom[C][D] = new TH1D(Form("Phip_%iMeVCM%iRandom", 415+(C*20), D+1), Form("#phi_{p} %i #pm 10MeV CM%i", 415+(C*20), D+1), 10, -180, 180);
        }
    }

}

void PNeutPol_Polarimeter_Circ::FillHists()
{
    ZpDist->Fill(pVertex(2), TaggerTime);

    // Cutting on MM as fn of energy HERE to avoid looping issues above
    double MMCutArray[6][2] = {{939.7, 14.1}, {943.3, 20.4}, {945.7, 28.5}, {954.6, 47.1}, {948.1, 41.9}, {948.7, 50.4}};

    for(Int_t M = 0; M <6; M++){ // For EGamma between 200-500 cut on MM differently for each energy
        double_t EgMMLow = 200 + (M*50);
        double_t EgMMHigh = 250 + (M*50);
        if( EgMMLow < EGamma && EGamma < EgMMHigh){ // Cut on missing mass depending upon enegry
            if (((MMpEpCorr < (MMCutArray[M][0] - (2*MMCutArray[M][1]))) == kTRUE) || (MMpEpCorr > (MMCutArray[M][0] + (2*MMCutArray[M][1]))) == kTRUE) return;
        }
    }

    if (200 > EGamma || EGamma > 500){ // If outside 200-700 MeV cut on 2Sigma of main MM peak
        if (((MMpEpCorr < 884.4) == kTRUE) || ((MMpEpCorr > 1003.6) == kTRUE)) return;
    }

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

        NeutronEPrompt->Fill(nKinE);
        NeutronEThetaScFullPrompt->Fill(nKinE, ScattTheta);
        NeutronEEgPrompt->Fill(nKinE, EGamma);
        ThetaScEgPrompt->Fill(ScattTheta, EGamma);

        ThetapThetanPrompt->Fill(ThpRad, nKin.Theta());

        DOCADistPrompt->Fill(DOCA);
        DOCArPOCAPrompt->Fill(r, DOCA);
        ThetaScPhiScPrompt->Fill(ScattThetaRad, ScattPhiRad);
        MMEgPrompt->Fill(EGamma, MMpEpCorr);
        PhiScEgPrompt->Fill(EGamma, ScattPhiRad);
        PhiScZpPrompt->Fill(pVertex(2), ScattPhiRad);

        if(BeamHelicity == kTRUE) PhiScPosHelFullPrompt->Fill(ScattPhiRad);
        else if(BeamHelicity == kFALSE) PhiScNegHelFullPrompt->Fill(ScattPhiRad);

        if((MWPC0pE != 0) && (MWPC1pE != 0)){ // If there is a non zero energy deposit in each MWPC for proton set one MWPC energy sum
            EpMWPCEpVetoTot2Prompt->Fill(EpCorr, MWPCpEVetoTot2);
            ThetapMWPCEpTot2Prompt->Fill(Thp, MWPCpETot2);
        }
        else if (((MWPC0pE != 0) && (MWPC1pE == 0)) == kTRUE || ((MWPC0pE == 0) && (MWPC1pE != 0)) == kTRUE){ // If only detect proton in one chamber set different ESum
            EpMWPCEpVetoTot1Prompt->Fill(EpCorr, MWPCpEVetoTot1);
            ThetapMWPCEpTot1Prompt->Fill(Thp, MWPCpETot1);
        }
        if((MWPC0nE != 0) && (MWPC1nE != 0)){ // If there is a non zero energy deposit in each MWPC for neutron set one MWPC energy sum
            EnMWPCEnVetoTot2Prompt->Fill(GVnCorr.E() - Mn, MWPCnEVetoTot2);
            ThetanMWPCEnTot2Prompt->Fill(Thn, MWPCnETot2);
        }
        else if (((MWPC0nE != 0) && (MWPC1nE == 0)) == kTRUE || ((MWPC0nE == 0) && (MWPC1nE != 0)) == kTRUE){ // If only detect neutron in one chamber set different ESum
            EnMWPCEnVetoTot1Prompt->Fill(GVnCorr.E() - Mn, MWPCnEVetoTot1);
            ThetanMWPCEnTot1Prompt->Fill(Thn, MWPCnETot1);
        }

        for(Int_t d = 0; d < 8; d++){ //Energy bins
            ELow = 290 + (d*EWidth);
            EHigh = 390 + (d*EWidth);
            if( ELow < EGamma && EGamma < EHigh){
                for(Int_t e = 0; e < 5; e++){
                    CosThetaLow = (0.6) - (e*CosThetaWidth);
                    CosThetaHigh = 1 - (e*CosThetaWidth);
                    if(CosThetaHigh > CosThetapCM && CosThetapCM > CosThetaLow){
//                        NeutronEThetaScPrompt[d][e]->Fill(nKinE, ScattTheta);
                        if(BeamHelicity == kTRUE) {
                            PhiScPosHelPrompt[d][e]->Fill(ScattPhiRad);
                            PhiScAEffPosHelPrompt[d][e]->Fill(ScattPhiRad, AEff);
                        }
                        else if(BeamHelicity == kFALSE){
                            PhiScNegHelPrompt[d][e]->Fill(ScattPhiRad);
                            PhiScAEffNegHelPrompt[d][e]->Fill(ScattPhiRad, AEff);
                        }
                    }
                }
            }
        }

        for(Int_t f = 0; f < 10; f++){ //Energy bins
            ELow2 = 415 + (f*20);
            EHigh2 = 435 + (f*20);
            if( ELow2 < EGamma && EGamma < EHigh2){
                for(Int_t g = 0; g < 5; g++){ //CM bins
                    CosThetaLow2 =  0.6 - (g*0.4);
                    CosThetaHigh2 = 1 - (g*0.4);
                    if(CosThetaHigh2 > CosThetapCM && CosThetapCM > CosThetaLow2){
                        PhipSetPrompt[f][g]->Fill(Php);
                    }
                }
            }
        }

        for(Int_t Y =0; Y < 12; Y++){
            double_t EMMLOW = 175 + (Y*50);
            double_t EMMHIGH = 225 + (Y*50);
            if( EMMLOW < EGamma && EGamma < EMMHIGH){
                MMpEgammaPrompt[Y]->Fill(MMpEpCorr);
                PhiScPrompt[Y]->Fill(ScattPhiRad, Wgt);
            }
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
        NeutronERandom->Fill(nKinE);
        NeutronEThetaScFullRandom->Fill(nKinE, ScattTheta);
        NeutronEEgRandom->Fill(nKinE, EGamma);
        ThetaScEgRandom->Fill(ScattTheta, EGamma);

        DOCADistRandom->Fill(DOCA);
        DOCArPOCARandom->Fill(r, DOCA);
        ThetaScPhiScRandom->Fill(ScattThetaRad, ScattPhiRad);
        MMEgRandom->Fill(EGamma, MMpEpCorr);

        PhiScEgRandom->Fill(EGamma, ScattPhiRad);
        PhiScZpRandom->Fill(pVertex(2), ScattPhiRad);

        if(BeamHelicity == kTRUE) PhiScPosHelFullRandom->Fill(ScattPhiRad);
        else if(BeamHelicity == kFALSE) PhiScNegHelFullRandom->Fill(ScattPhiRad);

        if((MWPC0pE != 0) && (MWPC1pE != 0)){ // If there is a non zero energy deposit in each MWPC for proton set one MWPC energy sum
            EpMWPCEpVetoTot2Random->Fill(EpCorr, MWPCpEVetoTot2);
            ThetapMWPCEpTot2Random->Fill(Thp, MWPCpETot2);
        }
        else if (((MWPC0pE != 0) && (MWPC1pE == 0)) == kTRUE || ((MWPC0pE == 0) && (MWPC1pE != 0)) == kTRUE){ // If only detect proton in one chamber set different ESum
            EpMWPCEpVetoTot1Random->Fill(EpCorr, MWPCpEVetoTot1);
            ThetapMWPCEpTot1Random->Fill(Thp, MWPCpETot1);
        }
        if((MWPC0nE != 0) && (MWPC1nE != 0)){ // If there is a non zero energy deposit in each MWPC for neutron set one MWPC energy sum
            EnMWPCEnVetoTot2Random->Fill(GVnCorr.E() - Mn, MWPCnEVetoTot2);
            ThetanMWPCEnTot2Random->Fill(Thn, MWPCnETot2);
        }
        else if (((MWPC0nE != 0) && (MWPC1nE == 0)) == kTRUE || ((MWPC0nE == 0) && (MWPC1nE != 0)) == kTRUE){ // If only detect neutron in one chamber set different ESum
            EnMWPCEnVetoTot1Random->Fill(GVnCorr.E() - Mn, MWPCnEVetoTot1);
            ThetanMWPCEnTot1Random->Fill(Thn, MWPCnETot1);
        }

        ThetapThetanRandom->Fill(ThpRad, nKin.Theta());

        for(Int_t d = 0; d < 8; d++){ //Energy bins
            ELow = 290 + (d*EWidth);
            EHigh = 390 + (d*EWidth);
            if( ELow < EGamma && EGamma < EHigh){
                for(Int_t e = 0; e < 5; e++){
                    CosThetaLow = (0.6) - (e*CosThetaWidth);
                    CosThetaHigh = 1 - (e*CosThetaWidth);
                    if(CosThetaHigh > CosThetapCM && CosThetapCM > CosThetaLow){
//                        NeutronEThetaScRandom[d][e]->Fill(nKinE, ScattTheta);
                        if(BeamHelicity == kTRUE) {
                            PhiScPosHelRandom[d][e]->Fill(ScattPhiRad);
                            PhiScAEffPosHelRandom[d][e]->Fill(ScattPhiRad, AEff);
                        }
                        else if(BeamHelicity == kFALSE){
                            PhiScNegHelRandom[d][e]->Fill(ScattPhiRad);
                            PhiScAEffNegHelRandom[d][e]->Fill(ScattPhiRad, AEff);
                        }
                    }
                }
            }
        }

        for(Int_t f = 0; f < 10; f++){ //Energy bins
            ELow2 = 415 + (f*20);
            EHigh2 = 435 + (f*20);
            if( ELow2 < EGamma && EGamma < EHigh2){
                for(Int_t g = 0; g < 5; g++){ //CM bins
                    CosThetaLow2 = 0.6 - (g*0.4);
                    CosThetaHigh2 = 1 - (g*0.4);
                    if(CosThetaHigh2 > CosThetapCM && CosThetapCM > CosThetaLow2){
                        PhipSetRandom[f][g]->Fill(Php);
                    }
                }
            }
        }


        for(Int_t Y =0; Y < 12; Y++){
            double_t EMMLOW = 175 + (Y*50);
            double_t EMMHIGH = 225 + (Y*50);
            if( EMMLOW < EGamma && EGamma < EMMHIGH){
                MMpEgammaRandom[Y]->Fill(MMpEpCorr);
                PhiScRandom[Y]->Fill(ScattPhiRad, Wgt);
            }
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

    NeutronE->Add(NeutronEPrompt);
    NeutronE->Add(NeutronERandom, -PvRratio);

    NeutronEThetaScFull->Add(NeutronEThetaScFullPrompt);
    NeutronEThetaScFull->Add(NeutronEThetaScFullRandom, -PvRratio);

    NeutronEEg->Add(NeutronEEgPrompt);
    NeutronEEg->Add(NeutronEEgRandom, -PvRratio);

    ThetaScEg->Add(ThetaScEgPrompt);
    ThetaScEg->Add(ThetaScEgRandom, -PvRratio);

    EpMWPCEpVetoTot1->Add(EpMWPCEpVetoTot1Prompt);
    EpMWPCEpVetoTot1->Add(EpMWPCEpVetoTot1Random, -PvRratio);
    EpMWPCEpVetoTot2->Add(EpMWPCEpVetoTot2Prompt);
    EpMWPCEpVetoTot2->Add(EpMWPCEpVetoTot2Random, -PvRratio);
    EnMWPCEnVetoTot1->Add(EnMWPCEnVetoTot1Prompt);
    EnMWPCEnVetoTot1->Add(EnMWPCEnVetoTot1Random, -PvRratio);
    EnMWPCEnVetoTot2->Add(EnMWPCEnVetoTot2Prompt);
    EnMWPCEnVetoTot2->Add(EnMWPCEnVetoTot2Random, -PvRratio);

    ThetapMWPCEpTot1->Add(ThetapMWPCEpTot1Prompt);
    ThetapMWPCEpTot1->Add(ThetapMWPCEpTot1Random, -PvRratio);
    ThetapMWPCEpTot2->Add(ThetapMWPCEpTot2Prompt);
    ThetapMWPCEpTot2->Add(ThetapMWPCEpTot2Random, -PvRratio);
    ThetanMWPCEnTot1->Add(ThetanMWPCEnTot1Prompt);
    ThetanMWPCEnTot1->Add(ThetanMWPCEnTot1Random, -PvRratio);
    ThetanMWPCEnTot2->Add(ThetanMWPCEnTot2Prompt);
    ThetanMWPCEnTot2->Add(ThetanMWPCEnTot2Random, -PvRratio);

    DOCADist->Add(DOCADistPrompt);
    DOCADist->Add(DOCADistRandom, -PvRratio);
    DOCArPOCA->Add(DOCArPOCAPrompt);
    DOCArPOCA->Add(DOCArPOCARandom, -PvRratio);
    ThetaScPhiSc->Add(ThetaScPhiScPrompt);
    ThetaScPhiSc->Add(ThetaScPhiScRandom, -PvRratio);
    PhiScPosHelFull->Add(PhiScPosHelFullPrompt);
    PhiScPosHelFull->Add(PhiScPosHelFullRandom, -PvRratio);
    PhiScNegHelFull->Add(PhiScNegHelFullPrompt);
    PhiScNegHelFull->Add(PhiScNegHelFullRandom, -PvRratio);
    MMEg->Add(MMEgPrompt);
    MMEg->Add(MMEgRandom, -PvRratio);
    PhiScEg->Add(PhiScEgPrompt);
    PhiScEg->Add(PhiScEgRandom, -PvRratio);
    PhiScZp->Add(PhiScZpPrompt);
    PhiScZp->Add(PhiScZpRandom, -PvRratio);

    for(Int_t E = 0; E < 8; E++){
        for(Int_t F = 0; F < 5; F++){
            PhiScNegHel[E][F]->Add(PhiScNegHelPrompt[E][F]);
            PhiScNegHel[E][F]->Add(PhiScNegHelRandom[E][F], -PvRratio);
            PhiScPosHel[E][F]->Add(PhiScPosHelPrompt[E][F]);
            PhiScPosHel[E][F]->Add(PhiScPosHelRandom[E][F], -PvRratio);
            PhiScAEffNegHel[E][F]->Add(PhiScAEffNegHelPrompt[E][F]);
            PhiScAEffNegHel[E][F]->Add(PhiScAEffNegHelRandom[E][F], -PvRratio);
            PhiScAEffPosHel[E][F]->Add(PhiScAEffPosHelPrompt[E][F]);
            PhiScAEffPosHel[E][F]->Add(PhiScAEffPosHelRandom[E][F], -PvRratio);
//            NeutronEThetaSc[E][F]->Add(NeutronEThetaScPrompt[E][F]);
//            NeutronEThetaSc[E][F]->Add(NeutronEThetaScRandom[E][F], -PvRratio);
        }
    }


    for(Int_t XY = 0; XY < 12; XY++){
        MMpEgamma[XY]->Add(MMpEgammaPrompt[XY]);
        MMpEgamma[XY]->Add(MMpEgammaRandom[XY], -PvRratio);
        PhiSc[XY]->Add(PhiScPrompt[XY]);
        PhiSc[XY]->Add(PhiScRandom[XY], -PvRratio);
    }

    for(Int_t G = 0; G < 10; G++){
        for(Int_t H = 0; H < 5; H++){
            PhipSet[G][H]->Add(PhipSetPrompt[G][H]);
            PhipSet[G][H]->Add(PhipSetRandom[G][H], -PvRratio);
        }
    }

}


Bool_t	PNeutPol_Polarimeter_Circ::Write(){
    GTreeManager::Write();
}
