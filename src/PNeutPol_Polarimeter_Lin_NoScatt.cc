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

    EWidth = 70; //Fix the width of the Energy and CosTheta bins used later, for Helicity plots
    EWidth2 = 10; // For Sigma plots
    CosThetaWidth = 0.4;
    CosThetaWidth2 = 0.1;

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

    Cut_CB_proton = OpenCutFile("configfiles/cuts/CB_DeltaE-E_Proton_24_4_18.root", "Proton"); // These will need adjusting with new Acqu files
    Cut_proton = Cut_CB_proton;
    Cut_CB_pion = OpenCutFile("configfiles/cuts/CB_DeltaE-E_Pion_29_07_15.root", "Pion");
    Cut_pion = Cut_CB_pion;
    Cut_CB_protonKinGood = OpenCutFile("configfiles/cuts/CB_DeltaE-E_ProtonKinGood_26_4_18.root", "ProtonKinGood"); // These will need adjusting with new Acqu files
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

    //if(Cut_proton -> IsInside(EpCorr, dEp) == kFALSE) return; // If E loss correct proton is NOT inside p banana drop out

    EpDiff = abs(EpCorr - Ep);

    Pp = sqrt (TMath::Power((Mp + EpCorr),2) - TMath::Power(Mp,2));
    Pn = sqrt (TMath::Power((En + Mn ),2) - TMath::Power(Mn,2));
    GVpCorr = TLorentzVector(Pp*sin(ThpRad)*cos(PhpRad), Pp*sin(ThpRad)*sin(PhpRad), Pp*cos(ThpRad), EpCorr+Mp);

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
        //pKinB = GVpCorr;
        pKinB.Boost(b);
        ThetapCM = (pKinB.Theta())*TMath::RadToDeg();
        CosThetapCM = cos (pKinB.Theta());
        nKinB = nKin;
        nKinB.Boost(b);
        ThetanCM = (nKinB.Theta())*TMath::RadToDeg();

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

        \\if (((MMpEpCorr < 800) == kTRUE) || ((MMpEpCorr > 1300) == kTRUE)) continue; // Very rough cut on MM

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
    PhiDiffDist = new GH1 ("PhiDiffDist", "#phi_{p} - #phi_{n}", 360, 0, 360);
    OAngle = new GH1 ("OAngle", "Opening Angle between P and N Vectors", 180, 0, 180);
    MMpEpCorrected = new GH1 ("MMpEpCorrected", "Missing mass seen by Proton (E Loss Corrected)", 400, 0, 2000);
    ZpDist = new GH1 ("ZpDist", "Proton Pseudo Z Vertex Distribution", 200, -400, 400);
    ThetanDist = new GH1 ("ThetanDist", "#theta_{n} Distribution", 200, 0, 180);
    ThetanDiffDist = new GH1("ThetanDiffDist", "#theta_{nMeas} - #theta_{nRec}", 360, -180, 180);

    Thetap = new TH1D ("Thetap", "#theta_{p} Distribution", 200, 0, 4);
    ThetapPrompt = new TH1D ("ThetapPrompt", "#theta_{p} Distribution", 200, 0, 4);
    ThetapRandom = new TH1D ("ThetapRandom", "#theta_{p} Distribution", 200, 0, 4);

    EdE = new TH2D ("EdE", "dE_{p} as a function of E_{p}", 200, 0, 1000, 200, 0, 10);
    EdEPrompt = new TH2D ("EdEPrompt", "dE_{p} as a function of E_{p}", 200, 0, 1000, 200, 0, 10);
    EdERandom = new TH2D ("EdERandom", "dE_{p} as a function of E_{p}", 200, 0, 1000, 200, 0, 10);
    E_dEKin = new GH2 ("E_dEKin", "E_{pKin}-dE", 100, 0, 500, 100, 0, 10);
    DeutKinPiKin = new GH2 ("DeutKinPiKin", "(#theta_{nRec} - #theta_{n}) vs (#theta_{#pi Rec} - #theta_{n})", 200, -180, 180, 200, -180, 180);

    Eg2 = new TH1D( "Eg2", "E_{#gamma} Distribution", 200, 100, 1600 );
    EgPrompt = new TH1D( "EgPrompt", "E_{#gamma} Distribution", 200, 100, 1600 );
    EgRandom = new TH1D( "EgRandom", "E_{#gamma} Distribution", 200, 100, 1600 );

    for(Int_t A = 0; A < 7; A++){
        for(Int_t B = 0; B < 5; B++){
            PhipPosHel[A][B] = new TH1D(Form("Phip%iPosHelCM%i", 265+(A*70), B+1), Form("#phi_{p} E_{#gamma}%i #pm 35MeV CM%i for +ve Helicity", 265+(A*70), B+1), 10, -4, 4);
            PhipPosHelPrompt[A][B] = new TH1D(Form("Phip%iPosHelCM%iPrompt", 265+(A*70), B+1), Form("#phi_{p} E_{#gamma}%i #pm 35MeV CM%i for +ve Helicity", 265+(A*70), B+1), 10, -4, 4);
            PhipPosHelRandom[A][B] = new TH1D(Form("Phip%iPosHelCM%iRandom", 265+(A*70), B+1), Form("#phi_{p} E_{#gamma}%i #pm 35MeV CM%i for +ve Helicity", 265+(A*70), B+1), 10, -4, 4);
            PhipNegHel[A][B] = new TH1D(Form("Phip%iNegHelCM%i", 265+(A*70), B+1), Form("#phi_{p} E_{#gamma}%i #pm 35MeV CM%i for -ve Helicity", 265+(A*70), B+1), 10, -4, 4);
            PhipNegHelPrompt[A][B] = new TH1D(Form("Phip%iNegHelCM%iPrompt", 265+(A*70), B+1), Form("#phi_{p} E_{#gamma}%i #pm 35MeV CM%i for -ve Helicity", 265+(A*70), B+1), 10, -4, 4);
            PhipNegHelRandom[A][B] = new TH1D(Form("Phip%iNegHelCM%iRandom", 265+(A*70), B+1), Form("#phi_{p} E_{#gamma}%i #pm 35MeV CM%i for -ve Helicity", 265+(A*70), B+1), 10, -4, 4);
        }
    }

    for(Int_t C = 0; C < 21; C++){
        for(Int_t D = 0; D < 20; D++){
            PhipSet[C][D] = new TH1D(Form("Phip_%iMeVCM%i", 410+(C*10), D+1), Form("#phi_{p} %i #pm 5MeV CM%i", 410+(C*10), D+1), 10, -180, 180);
            PhipSetPrompt[C][D] = new TH1D(Form("Phip_%iMeVCM%iPrompt", 410+(C*10), D+1), Form("#phi_{p} %i #pm 5MeV CM%i", 410+(C*10), D+1), 10, -180, 180);
            PhipSetRandom[C][D] = new TH1D(Form("Phip_%iMeVCM%iRandom", 410+(C*10), D+1), Form("#phi_{p} %i #pm 5MeV CM%i", 410+(C*10), D+1), 10, -180, 180);
        }
    }

    for(Int_t X = 0; X < 12; X++){
        MMpEgamma[X] = new TH1D(Form("MMp_%iMeV", 225+(X*50)), Form("#MM_{p} %i #pm 25MeV", 225+(X*50)), 400, 0, 2000);
        MMpEgammaPrompt[X] = new TH1D(Form("MMp_%iMeVPrompt", 225+(X*50)), Form("#MM_{p} %i #pm 25MeV", 225+(X*50)), 400, 0, 2000);
        MMpEgammaRandom[X] = new TH1D(Form("MMp_%iMeVRandom", 225+(X*50)), Form("#MM_{p} %i #pm 25MeV", 225+(X*50)), 400, 0, 2000);
    }

}

void PNeutPol_Polarimeter_Lin_NoScatt::FillHists()
{
    // Cutting on MM as fn of energy HERE to avoid looping issues above
    double MMCutArray[10][2] = {{937, 14.4}, {939.5, 20.4}, {942.2, 29.1}, {944.1, 36.5}, {946, 49.7}, {952.2, 59.3}, {951, 62.7}, {956.4, 69.9}, {959.0, 79.9}, {981.3, 109.2}};

    for(Int_t M = 0; M <10; M++){ // For EGamma between 200-700 cut on MM differently for each energy
        double_t EgMMLow = 200 + (M*50);
        double_t EgMMHigh = 250 + (M*50);
        if( EgMMLow < EGamma && EGamma < EgMMHigh){ // Cut on missing mass depending upon enegry
            if (((MMpEpCorr < (MMCutArray[M][0] - (2*MMCutArray[M][1]))) == kTRUE) || (MMpEpCorr > (MMCutArray[M][0] + (2*MMCutArray[M][1]))) == kTRUE) return;
        }
    }

    if (200 > EGamma || EGamma > 700){ // If outside 200-700 MeV cut on 2Sigma of main MM peak
        if (((MMpEpCorr < 872.1) == kTRUE) || ((MMpEpCorr > 950) == kTRUE)) return;
    }

    time->Fill(TaggerTime);
    if (-5 < TaggerTime && TaggerTime < 20) time_cut->Fill(TaggerTime);

    Eg->Fill(EGamma, TaggerTime);
    PhiDiffDist->Fill(PhiDiff, TaggerTime);
    OAngle->Fill(OpeningAngle, TaggerTime);
    MMpEpCorrected->Fill(MMpEpCorr, TaggerTime);
    ZpDist->Fill(Zp, TaggerTime);
    ThetanDist->Fill(Thn, TaggerTime);
    ThetanDiffDist->Fill(ThetanCorr-ThetanRec, TaggerTime);

    E_dEKin->Fill(KinEp, dEp, TaggerTime);
    DeutKinPiKin->Fill(ThetanRec-Thn, ThetaPiRecDiff, TaggerTime);

    if (GHistBGSub::IsPrompt(TaggerTime) == kTRUE){
        EgPrompt->Fill(EGamma);
        EdEPrompt->Fill(EpCorr, dEp);
        ThetapPrompt->Fill(ThpRad);

        for(Int_t d = 0; d < 7; d++){ //Energy bins
            ELow = 230 + (d*EWidth);
            EHigh = 300 + (d*EWidth);
            if( ELow < EGamma && EGamma < EHigh){
                for(Int_t e = 0; e < 5; e++){
                    CosThetaLow = 0.6 - (e*CosThetaWidth);
                    CosThetaHigh = 1 - (e*CosThetaWidth);
                    if(CosThetaHigh > CosThetapCM && CosThetapCM > CosThetaLow){
                        if(BeamHelicity == kTRUE) PhipPosHelPrompt[d][e]->Fill(PhpRad);
                        else if(BeamHelicity == kFALSE) PhipNegHelPrompt[d][e]->Fill(PhpRad);
                    }
                }
            }
        }

        for(Int_t f = 0; f < 21; f++){ //Energy bins
            ELow2 = 410 + (f*EWidth2);
            EHigh2 = 420 + (f*EWidth2);
            if( ELow2 < EGamma && EGamma < EHigh2){
                for(Int_t g = 0; g < 20; g++){ //CM bins
                    CosThetaLow2 = 0.9 - (g*CosThetaWidth2);
                    CosThetaHigh2 = 1 - (g*CosThetaWidth2);
                    if(CosThetaHigh2 > CosThetapCM && CosThetapCM > CosThetaLow2){
                        PhipSetPrompt[f][g]->Fill(Php);
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
        EdERandom->Fill(EpCorr, dEp);
        ThetapRandom->Fill(ThpRad);

        for(Int_t d = 0; d < 7; d++){ //Energy bins
            ELow = 230 + (d*EWidth);
            EHigh = 300 + (d*EWidth);
            if( ELow < EGamma && EGamma < EHigh){
                for(Int_t e = 0; e < 5; e++){
                    CosThetaLow = 0.6 - (e*CosThetaWidth);
                    CosThetaHigh = 1 - (e*CosThetaWidth);
                    if(CosThetaHigh > CosThetapCM && CosThetapCM > CosThetaLow){
                        if(BeamHelicity == kTRUE) PhipPosHelRandom[d][e]->Fill(PhpRad);
                        else if(BeamHelicity == kFALSE) PhipNegHelRandom[d][e]->Fill(PhpRad);
                    }
                }
            }
        }

        for(Int_t f = 0; f < 21; f++){ //Energy bins
            ELow2 = 410 + (f*EWidth2);
            EHigh2 = 420 + (f*EWidth2);
            if( ELow2 < EGamma && EGamma < EHigh2){
                for(Int_t g = 0; g < 20; g++){ //CM bins
                    CosThetaLow2 = 0.9 - (g*CosThetaWidth2);
                    CosThetaHigh2 = 1 - (g*CosThetaWidth2);
                    if(CosThetaHigh2 > CosThetapCM && CosThetapCM > CosThetaLow2){
                        PhipSetRandom[f][g]->Fill(Php);
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

void PNeutPol_Polarimeter_Lin_NoScatt::BGSub(){

    Eg2->Add(EgPrompt);
    Eg2->Add(EgRandom, -PvRratio);
    EdE->Add(EdEPrompt);
    EdE->Add(EdERandom, -PvRratio);
    Thetap->Add(ThetapPrompt);
    Thetap->Add(ThetapRandom, -PvRratio);

    for(Int_t E = 0; E < 7; E++){
        for(Int_t F = 0; F < 5; F++){
            PhipNegHel[E][F]->Add(PhipNegHelPrompt[E][F]);
            PhipNegHel[E][F]->Add(PhipNegHelRandom[E][F], -PvRratio);
            PhipNegHel[E][F]->Add(PhipNegHelPrompt[E][F]);
            PhipNegHel[E][F]->Add(PhipNegHelRandom[E][F], -PvRratio);
        }
    }

    for(Int_t G = 0; G < 21; G++){
        for(Int_t H = 0; H < 20; H++){
            PhipSet[G][H]->Add(PhipSetPrompt[G][H]);
            PhipSet[G][H]->Add(PhipSetRandom[G][H], -PvRratio);
        }
    }

    for(Int_t XY = 0; XY < 12; XY++){
        MMpEgamma[XY]->Add(MMpEgammaPrompt[XY]);
        MMpEgamma[XY]->Add(MMpEgammaRandom[XY], -PvRratio);
    }

}

Bool_t	PNeutPol_Polarimeter_Lin_NoScatt::Write(){
  // Write all GH1's easily

  GTreeManager::Write();
}
