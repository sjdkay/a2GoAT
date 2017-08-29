#ifndef __CINT__

#include "PPhysics.h"

PPhysics::PPhysics()
{
	TC_cut_min = 0;
	TC_cut_max = 352;
}

PPhysics::~PPhysics()
{
}

Bool_t	PPhysics::Init()
{
	return kTRUE;
}

void	PPhysics::Reconstruct()
{
}

// ----------------------------------------------------------------------------------------
// TH1 routines
// ----------------------------------------------------------------------------------------
void PPhysics::FillScalers(Int_t low_scaler_number, Int_t high_scaler_number, TH1* hist)
{
	Int_t nFillScalers = high_scaler_number - low_scaler_number + 1;

	if( nFillScalers < hist->GetNbinsX())
	{
	    cout << "Error: FillScalers - histogram has insufficient bins for range" << endl;
	    return;
	}

	// To properly accumulate, create a histogram for this scaler read
	// cloning input histogram means the axis will be equivalent
	TH1* hist_current_SR = (TH1D*) hist->Clone();
	hist_current_SR->Reset();

    // Loop over scaler range, don't pull anything higher than the real # GetScalers()
	if (low_scaler_number  < 0)
	{
		cout << "FillScalers given scaler number outside range: " << low_scaler_number << endl;
		cout << "Setting lower limit to zero and continuing" << endl;
		low_scaler_number = 0;
	}
    if (high_scaler_number > GetScalers()->GetNScalers())
	{
		cout << "FillScalers given scaler number outside range: " << high_scaler_number << endl;
		cout << "Setting upper limit to "<< high_scaler_number << " and continuing" << endl;
        high_scaler_number = GetScalers()->GetNScalers();
	}

    for (Int_t i = low_scaler_number; i <= high_scaler_number; i++)
	{
		Int_t bin = i - low_scaler_number;
        hist_current_SR->SetBinContent(bin,GetScalers()->GetScaler(i));
	}

	// Add to accumulated
	hist->Add(hist_current_SR);
}

void PPhysics::FillMissingMass(const GTreeParticle& tree, TH1* Hprompt, TH1* Hrandom)
{
    for (Int_t i = 0; i < tree.GetNParticles(); i++)
    {
    for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
	{
        	FillMissingMass(tree, i, j, Hprompt, Hrandom);
	}
    }
}

void PPhysics::FillMissingMass(const GTreeParticle& tree, Int_t particle_index, TH1* Hprompt, TH1* Hrandom)
{
    for (Int_t i = 0; i < GetTagger()->GetNTagged(); i++)
	{
        	FillMissingMass(tree, particle_index, i, Hprompt, Hrandom);
	}
}

void PPhysics::FillMissingMass(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index, TH1* Hprompt, TH1* Hrandom)
{
    time = GetTagger()->GetTaggedTime(tagger_index) - tree.GetTime(particle_index);
    missingp4 = CalcMissingP4(tree, particle_index,tagger_index);

	if (GHistBGSub::IsPrompt(time)) Hprompt->Fill(missingp4.M());
	if (GHistBGSub::IsRandom(time)) Hrandom->Fill(missingp4.M());
}

void PPhysics::FillTime(const GTreeParticle& tree, TH1* Hist)
{
    for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
        for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
		{
            // Is tagger channel rejected by user?
            if(GetTagger()->GetTaggedChannel(j) < TC_cut_min) continue;
            if(GetTagger()->GetTaggedChannel(j) > TC_cut_max) continue;

                time = GetTagger()->GetTaggedTime(j) - tree.GetTime(i);
			Hist->Fill(time);
		}
	}
}

void PPhysics::FillTime(const GTreeParticle& tree, Int_t particle_index, TH1* Hist)
{
    for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
	{
    // Is tagger channel rejected by user?
        if(GetTagger()->GetTaggedChannel(j) < TC_cut_min) continue;
        if(GetTagger()->GetTaggedChannel(j) > TC_cut_max) continue;

        time = GetTagger()->GetTaggedTime(j) - tree.GetTime(particle_index);
	Hist->Fill(time);
	}
}

void PPhysics::FillTimeCut(const GTreeParticle& tree, TH1* Hist)
{
    for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
        for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
		{
            // Is tagger channel rejected by user?
            if(GetTagger()->GetTaggedChannel(j) < TC_cut_min) continue;
            if(GetTagger()->GetTaggedChannel(j) > TC_cut_max) continue;

                    time = GetTagger()->GetTaggedTime(j) - tree.GetTime(i);
			if((GHistBGSub::IsPrompt(time)) || (GHistBGSub::IsRandom(time))) Hist->Fill(time);
		}
	}
}

void PPhysics::FillTimeCut(const GTreeParticle& tree, Int_t particle_index, TH1* Hist)
{
    for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
	{
        // Is tagger channel rejected by user?
        if(GetTagger()->GetTaggedChannel(j) < TC_cut_min) continue;
        if(GetTagger()->GetTaggedChannel(j) > TC_cut_max) continue;

        time = GetTagger()->GetTaggedTime(j) - tree.GetTime(particle_index);
		if((GHistBGSub::IsPrompt(time)) || (GHistBGSub::IsRandom(time))) Hist->Fill(time);
	}
}

void PPhysics::FillMass(const GTreeParticle& tree, TH1* Hist)
{
    for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
        Hist->Fill(tree.GetMass(i));
	}
}

void PPhysics::FillMass(const GTreeParticle& tree, Int_t particle_index, TH1* Hist)
{
    Hist->Fill(tree.GetMass(particle_index));
}

void PPhysics::FillBeamAsymmetry(const GTreeParticle& tree, Int_t particle_index, TH1* Hprompt, TH1* Hrandom, Double_t MM_min, Double_t MM_max)
{
    for (Int_t i = 0; i < GetTagger()->GetNTagged(); i++)
    {
        FillBeamAsymmetry(tree, particle_index, i, Hprompt, Hrandom, MM_min, MM_max);
    }

}

void PPhysics::FillBeamAsymmetry(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index, TH1* Hprompt, TH1* Hrandom, Double_t MM_min, Double_t MM_max)
{
    // Is tagger channel rejected by user?
//    cout << tagger->GetTaggedChannel(tagger_index) << " " << TC_cut_min << " " << TC_cut_max << endl;
    if(GetTagger()->GetTaggedChannel(tagger_index) < TC_cut_min) return;
    if(GetTagger()->GetTaggedChannel(tagger_index) > TC_cut_max) return;

    // calc particle time diff
    time = GetTagger()->GetTaggedTime(tagger_index) - tree.GetTime(particle_index);
//    cout << "time " << time << endl;
    // calc missing p4
    missingp4 = CalcMissingP4(tree, particle_index,tagger_index);
 //   cout << "MM " << missingp4.M() << endl;
    if((missingp4.M() < MM_min) || (missingp4.M() > MM_max)) return;

   if (GHistBGSub::IsPrompt(time)) Hprompt->Fill(tree.GetPhi(particle_index)); //cout << "prompt" << endl;}
   if (GHistBGSub::IsRandom(time)) Hrandom->Fill(tree.GetPhi(particle_index));	//cout << "random" << endl;}
}

Double_t PPhysics::CalcCoplanarity(const GTreeParticle& tree1, Int_t particle_index1, const GTreeParticle& tree2, Int_t particle_index2)
{
   Double_t phi1 = tree1.GetPhi(particle_index1);
   Double_t phi2 = tree2.GetPhi(particle_index2);
   Double_t phidiff = TMath::Abs(phi1 - phi2);

   return phidiff;
}

Double_t PPhysics::CalcKinEnergy(Double_t PrimaryTheta, Double_t BeamEnergy, Double_t TargetMass, Double_t BeamMass, Double_t PrimaryMass, Double_t SecondaryMass) // Not working correctly?
{
    // Primary Theta/Mass is the angle of the particle we want to calculate the energy for
    // Secondary mass is the mass of the other particle in the reaction
    // Adapted from fortran fn, function takes Proton theta and the beam energy
    // to calculate the initial energy of the proton

    // Fortran (Dan) Version

    Double_t PrimaryThetaRad = PrimaryTheta*TMath::DegToRad(); // Convert input theta to radians
    Double_t Beta = cos(PrimaryThetaRad);
    Double_t P0 = BeamEnergy + BeamMass;
    Double_t E0 = (BeamEnergy + TargetMass); //Caution! Ensure both in same units! 1875.613 is Deuterium mass in MeV/C^2
    Double_t M2 = (((TMath::Power(E0,2))-(TMath::Power(P0,2))) + (TMath::Power(PrimaryMass,2)) - (TMath::Power(SecondaryMass,2))); //M2a + proton mass squared - Neutron mass squared

    Double_t a = 4*((TMath::Power((Beta*P0),2)) - (TMath::Power(E0,2)));
    Double_t b = 4*Beta*M2*P0;
    Double_t c = (TMath::Power(M2,2)) - (4*((TMath::Power((E0*PrimaryMass),2))));

    Double_t d = (TMath::Power(b,2)) - (4*a*c);
    Double_t P = (-b - (sqrt(d)))/(2*a); // Magnitude of 4-momentum for primary
    Double_t P_Energy_a = sqrt((TMath::Power(P,2)) + (TMath::Power(PrimaryMass,2))); // Energy of primary
    Double_t P_Energy_b = P_Energy_a - PrimaryMass; // Kinetic energy of primary

    return P_Energy_b;
}

// Calculate kinetic energy of proton from CB energy and proton angle with polarimeter in place
// Uses parameterisation worked out by Mikhail Bashkanov
Double_t PPhysics::EpPolCorrect(Double_t ProtE, Double_t ProtTheta){

    A = CoeffA(ProtTheta);
    B = CoeffB(ProtTheta);
    C = CoeffC(ProtTheta);

    Double_t EKinMB = (A*exp(B*ProtE)) + C + ProtE;

    return EKinMB;
}

Double_t PPhysics::CalcKinEnergyMB(Double_t PrimaryTheta, Double_t BeamEnergy, Double_t TargetMass, Double_t BeamMass, Double_t PrimaryMass, Double_t SecondaryMass)
{
    // Mikhail version
    // Working for some values but others return a NAN error!

    Double_t PrimaryThetaRad = PrimaryTheta*TMath::DegToRad(); // Convert input theta to radians
    Double_t Beta = cos(PrimaryThetaRad);

    Double_t E1MB = BeamEnergy + BeamMass;
    Double_t M1MB = BeamMass;
    Double_t M2MB = TargetMass;
    Double_t P1MB = sqrt((TMath::Power(E1MB,2)) - (TMath::Power(M1MB,2)));
    Double_t P2MB = 0;
    Double_t A1MB = (2*M2MB) + (2*E1MB);
    Double_t A2MB = ((TMath::Power(SecondaryMass,2)) - (TMath::Power(M1MB,2)) - ((TMath::Power(PrimaryMass,2)) - (2*E1MB*M2MB) - (TMath::Power(M2MB,2))));

    Double_t V1MB = 2*A1MB*A2MB;
    Double_t V2MB = 0.;
    Double_t V3MB = (-4*(TMath::Power(P1MB,2)))*(((TMath::Power(PrimaryMass,2))*(TMath::Power(A1MB,2)))-(TMath::Power(A2MB,2)));
    Double_t V4MB = 16*(TMath::Power(PrimaryMass,2)*(TMath::Power(P1MB,4)));
    Double_t V5MB = TMath::Power(A1MB,2);
    Double_t V6MB = -4*TMath::Power(P1MB,2);

    Double_t V2MB2 =0.;
    Double_t V3MB2 =4*V3MB;
    Double_t V4MB2 =4*V4MB;

    Double_t e5MB = (V1MB + sqrt(V2MB2 + (V3MB2*TMath::Power(Beta,2))+(V4MB2*TMath::Power(Beta,4))))/2/(V5MB+(V6MB*TMath::Power(Beta,2)))-PrimaryMass;
    Double_t e6MB = (V1MB - sqrt(V2MB2 + (V3MB2*TMath::Power(Beta,2))+(V4MB2*TMath::Power(Beta,4))))/2/(V5MB+(V6MB*TMath::Power(Beta,2)))-PrimaryMass;;

    if((PrimaryThetaRad < 0.5*acos(-1)) == kTRUE) e4MB = e5MB;
    else if ((PrimaryThetaRad < 0.5*acos(-1)) == kFALSE) e4MB = e6MB;

    return e4MB;
}

Double_t PPhysics::CoeffA(Double_t ProtTheta){ // Calculate Coefficient A for MB parameterisation

    Double_t Theta = ProtTheta * TMath::DegToRad();
    Double_t CoeA = 201.915 - (57.9314*sin(Theta));

    return CoeA;
}

Double_t PPhysics::CoeffB(Double_t ProtTheta){ // Calculate Coefficient B for MB parameterisation

    Double_t Theta = ProtTheta * TMath::DegToRad();
    Double_t CoeB = -0.000800067 - (0.00451967*sin(Theta));

    return CoeB;
}

Double_t PPhysics::CoeffC(Double_t ProtTheta){ // Calculate Coefficient C for MB parameterisation

    Double_t Theta = ProtTheta * TMath::DegToRad();
    Double_t CoeC = -82.3023 + (23.2409*sin(Theta));

    return CoeC;
}

TVector3 PPhysics::ScatteredFrameAngles(TVector3 InitialVect, TVector3 RealPVect, TVector3 ScattVector, TLorentzVector GammaVect)
{
    TVector3 ValueHolder;

    TVector3 RecNeut = (InitialVect.Unit()); //neutron angle, z-axis
    //cout << "Rec Neut" << "   " << RecNeut(0) << "   " << RecNeut(1) << "   " << RecNeut(2) << endl;
    TVector3 ZAxis = (InitialVect.Unit()); //neutron angle, z-axis
    //cout << "Z Axis" << "   " << ZAxis(0) << "   " << ZAxis(1) << "   " << ZAxis(2) << endl;

    TVector3 ScattNeut = (ScattVector.Unit()); //recoil proton angle
    double_t TT1 = RecNeut.Angle(ScattNeut);
    double_t tmpPh = TMath::ATan2((RecNeut.Py()), (RecNeut.Px()))-0.5*acos(-1);
    //cout << "Scatt Neut" << "   " << ScattNeut(0) << "   " << ScattNeut(1) << "   " << ScattNeut(2) << endl;

    TVector3 BeamVect = (GammaVect.Vect()).Unit();
    TVector3 ProtVect = (RealPVect.Unit());
    TVector3 YAxis = BeamVect.Cross(ProtVect); // Y-Axis
    TVector3 XAxis = YAxis.Cross(ZAxis); // X-Axis

    //cout << "Y Axis" << "   " << YAxis(0) << "   " << YAxis(1) << "   " << YAxis(2) << endl;
    //cout << "X Axis" << "   " << XAxis(0) << "   " << XAxis(1) << "   " << XAxis(2) << endl;

    double_t tmpL = cos(TT1);
    //cout << "tmpL" << "   " << tmpL << endl;
    TVector3 V1(RecNeut.X()*tmpL, RecNeut.Y()*tmpL, RecNeut.Z()*tmpL); // Projection of lengths of n vector in scattered frame
    TVector3 V2 = ScattNeut-V1; // XY projection of recoil vector in recoil frame
    //cout << "V1" << "   " << V1(0) << "   " << V1(1) << "   " << V1(2) << endl;
    //cout << "V2" << "   " << V2(0) << "   " << V2(1) << "   " << V2(2) << endl;
    //cout << YAxis.Angle(V2) << "   " << XAxis.Angle(V2) << endl;

    double_t Phi = TMath::ATan2(cos(YAxis.Angle(V2)), cos(XAxis.Angle(V2)));

    ValueHolder.SetXYZ(Phi*TMath::RadToDeg(), tmpPh, TT1*TMath::RadToDeg());

    return ValueHolder;
}

Double_t PPhysics::Calc_dtfInterDOCA(const TVector3 &locUnitDir1, const TVector3 &locUnitDir2, const TVector3 &locVertex1, const TVector3 &locVertex2, TVector3 &locInterDOCA1, TVector3 &locInterDOCA2){
    //originated from code by Jörn Langheinrich
    //you can use this function to find the DOCA to a fixed point by calling this function with locUnitDir1 and 2 parallel, and the fixed vertex as locVertex2
    Double_t locUnitDot = locUnitDir1*locUnitDir2;
    Double_t locDenominator = locUnitDot*locUnitDot - 1.0; /// scalar product of directions
    Double_t locDistVertToInterDOCA1 = 0.0, locDistVertToInterDOCA2 = 0.0; //distance from vertex to DOCA point

    if(fabs(locDenominator) < 1.0e-15) //parallel
        locDistVertToInterDOCA1 = (locVertex2 - locVertex1)*locUnitDir2/locUnitDot; //the opposite
    else{
        Double_t locA = (locVertex1 - locVertex2)*locUnitDir1;
        Double_t locB = (locVertex1 - locVertex2)*locUnitDir2;
        locDistVertToInterDOCA1 = (locA - locUnitDot*locB)/locDenominator;
        locDistVertToInterDOCA2 = (locUnitDot*locA - locB)/locDenominator;
    }

    locInterDOCA1 = locVertex1 + locDistVertToInterDOCA1*locUnitDir1; //intersection point of DOCA line and track 1
    locInterDOCA2 = locVertex2 + locDistVertToInterDOCA2*locUnitDir2; //intersection point of DOCA line and track 2
    Double_t locDOCA = (locInterDOCA1 - locInterDOCA2).Mag();
    return ((locVertex2.Z() > locVertex1.Z()) ? locDOCA : -1.0*locDOCA); // Returns DOCA, also want POCA
}

Bool_t PPhysics::MCDataCheck(){

    Bool_t MC;
    if (GetScalers()->GetNEntries() == 0) // MC Data has no scaler entries so if 0, data gets a flag to denote it as MC
    {
        MC = kTRUE;
    }

    else if (GetScalers()->GetNEntries() != 0) // This flag is listed as false if the number of scaler entries does not equal 0
    {
        MC = kFALSE;
    }

    return MC;
}

// ----------------------------------------------------------------------------------------
// GH1 routines
// ----------------------------------------------------------------------------------------

void PPhysics::FillMissingMass(const GTreeParticle& tree, GH1* gHist, Bool_t TaggerBinning)
{
	for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
        for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
		{
			FillMissingMass(tree, i, j, gHist, TaggerBinning);
		}
	}
}

void PPhysics::FillMissingMass(const GTreeParticle& tree, Int_t particle_index, GH1* gHist, Bool_t TaggerBinning)
{
    for (Int_t i = 0; i < GetTagger()->GetNTagged(); i++)
	{
        FillMissingMass(tree, particle_index, i, gHist, TaggerBinning);
	}
}

void PPhysics::FillMissingMass(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index, GH1* gHist, Bool_t TaggerBinning)
{
    // Is tagger channel rejected by user?
    if(GetTagger()->GetTaggedChannel(tagger_index) < TC_cut_min) return;
    if(GetTagger()->GetTaggedChannel(tagger_index) > TC_cut_max) return;

    // calc particle time diff
    time = GetTagger()->GetTaggedTime(tagger_index) - tree.GetTime(particle_index);

    // calc missing p4
    missingp4 = CalcMissingP4(tree, particle_index,tagger_index);

   // Fill GH1
   if(TaggerBinning)   gHist->Fill(missingp4.M(),time, GetTagger()->GetTaggedChannel(tagger_index));
   else gHist->Fill(missingp4.M(),time);

}

Double_t PPhysics::CalcMissingMass(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index)
{
    missingp4 	= CalcMissingP4(tree, particle_index, tagger_index);

	return missingp4.M();
}

Double_t PPhysics::CalcMissingEnergy(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index)
{
    missingp4 	= CalcMissingP4(tree,particle_index, tagger_index);

	return missingp4.T();
}

TLorentzVector PPhysics::CalcMissingP4(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index)
{
    particle	= tree.Particle(particle_index);
    beam 		= TLorentzVector(0.,0.,GetTagger()->GetTaggedEnergy(tagger_index),GetTagger()->GetTaggedEnergy(tagger_index));
	missingp4 	= beam + target - particle;

	return missingp4;
}

void PPhysics::FillBeamAsymmetry(const GTreeParticle& tree, Int_t particle_index, GH1* gHist, Bool_t TaggerBinning, Double_t MM_min, Double_t MM_max)
{
    for (Int_t i = 0; i < GetTagger()->GetNTagged(); i++)
    {
        FillBeamAsymmetry(tree, particle_index, i, gHist, TaggerBinning);
    }

}

void PPhysics::FillBeamAsymmetry(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index, GH1* gHist, Bool_t TaggerBinning, Double_t MM_min, Double_t MM_max)
{
    // Is tagger channel rejected by user?
    if(GetTagger()->GetTaggedChannel(tagger_index) < TC_cut_min) return;
    if(GetTagger()->GetTaggedChannel(tagger_index) > TC_cut_max) return;

    // calc particle time diff
    time = GetTagger()->GetTaggedTime(tagger_index) - tree.GetTime(particle_index);

    // calc missing p4
    missingp4 = CalcMissingP4(tree, particle_index,tagger_index);
    if((missingp4.M() < MM_min) || (missingp4.M() > MM_max)) return;

   if(TaggerBinning)   gHist->Fill(tree.GetPhi(particle_index),time,GetTagger()->GetTaggedChannel(tagger_index));
   else gHist->Fill(tree.GetPhi(particle_index),time);

}

void PPhysics::FillTime(const GTreeParticle& tree, GH1* gHist)
{
    for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
        for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
		{
            // Is tagger channel rejected by user?
            if(GetTagger()->GetTaggedChannel(j) < TC_cut_min) continue;
            if(GetTagger()->GetTaggedChannel(j) > TC_cut_max) continue;

                time = GetTagger()->GetTaggedTime(j) - tree.GetTime(i);
			gHist->Fill(time);
		}
	}
}

void PPhysics::FillTime(const GTreeParticle& tree, Int_t particle_index, GH1* gHist)
{
    for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
	{
        // Is tagger channel rejected by user?
        if(GetTagger()->GetTaggedChannel(j) < TC_cut_min) continue;
        if(GetTagger()->GetTaggedChannel(j) > TC_cut_max) continue;

        time = GetTagger()->GetTaggedTime(j) - tree.GetTime(particle_index);
		gHist->Fill(time);
	}
}

void PPhysics::FillTimeCut(const GTreeParticle& tree, GH1* gHist)
{
    for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
        for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
		{
            // Is tagger channel rejected by user?
            if(GetTagger()->GetTaggedChannel(j) < TC_cut_min) continue;
            if(GetTagger()->GetTaggedChannel(j) > TC_cut_max) continue;

                    time = GetTagger()->GetTaggedTime(j) - tree.GetTime(i);
			if((GHistBGSub::IsPrompt(time)) || (GHistBGSub::IsRandom(time))) gHist->Fill(time);
		}
	}
}

void PPhysics::FillTimeCut(const GTreeParticle& tree, Int_t particle_index, GH1* gHist)
{
    for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
	{
        // Is tagger channel rejected by user?
        if(GetTagger()->GetTaggedChannel(j) < TC_cut_min) continue;
        if(GetTagger()->GetTaggedChannel(j) > TC_cut_max) continue;

        time = GetTagger()->GetTaggedTime(j) - tree.GetTime(particle_index);
		if((GHistBGSub::IsPrompt(time)) || (GHistBGSub::IsRandom(time))) gHist->Fill(time);
	}
}

void PPhysics::FillMass(const GTreeParticle& tree, GH1* gHist)
{
    for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
        gHist->Fill(tree.GetMass(i));
	}
}

void PPhysics::FillMass(const GTreeParticle& tree, Int_t particle_index, GH1* gHist)
{
    gHist->Fill(tree.GetMass(particle_index));
}

Bool_t 	PPhysics::Write()
{
	return kTRUE;
}

// Some common initialisation stuff
Bool_t 	PPhysics::InitBackgroundCuts()
{
	// Set background cuts
	Double_t p1, p2, r1, r2;
	string config = ReadConfig("Set-Prompt-Cut");
	if(strcmp(config.c_str(), "nokey") == 0)
		cout << "No BG subtraction - At least 1 prompt and random cut required" << endl;
	else if(sscanf( config.c_str(), "%lf %lf\n", &p1, &p2) == 2)
	{
	   config = ReadConfig("Add-Random-Cut",0);
	   if(strcmp(config.c_str(), "nokey") == 0)
	   	cout << "No BG subtraction - At least 1 random cut required" << endl;
	   else if(sscanf( config.c_str(), "%lf %lf\n", &r1, &r2) == 2)
	   {
		cout << "Init BG cuts:" << endl;
		cout << "prompt(" << p1 << "," << p2 << ") ";
		cout << "random(" << r1 << "," << r2 << ") " << endl;

		GHistBGSub::InitCuts(p1,p2,r1,r2);

		// Look for additional random windows
		Int_t instance = 1;
		do
		{
			config = ReadConfig("Add-Random-Cut",instance);
			if(sscanf( config.c_str(), "%lf %lf\n", &r1, &r2) == 2)
			{
				cout << "Adding random cuts: ";
				cout << "random(" << r1 << "," << r2 << ") " << endl;

				GHistBGSub::AddRandCut(r1,r2);
			}
			instance++;
		} while (strcmp(config.c_str(), "nokey") != 0);
	   }
	   else {cout << "Random window not set correctly" << endl; return kFALSE;}
	}
	else {cout << "Prompt window not set correctly" << endl; return kFALSE;}

	cout << endl;
	return kTRUE;

}

Bool_t 	PPhysics::InitTargetMass()
{
	Double_t mass;
	string config = ReadConfig("Target-Mass");
	if(strcmp(config.c_str(), "nokey") == 0)
	{
		cout << "Target mass unknown!" << endl;
	}
	else if(sscanf( config.c_str(), "%lf\n", &mass) == 1)
	{
		cout << "Setting Target mass: " << mass << " MeV" << endl;
		SetTarget(mass);
	}
	else
	{
		cout << "Target Mass not set correctly" << endl;
		return kFALSE;
	}

	cout << endl;
	return kTRUE;

}

Bool_t 	PPhysics::InitTaggerChannelCuts()
{
	Double_t tc1, tc2;
	string config = ReadConfig("Tagger-Channel-Cut");
	if(sscanf( config.c_str(), "%lf %lf\n", &tc1, &tc2) == 2)
	{
		if ((tc1 < 0) || (tc1 > 352))
		{
           cout << "Invalid tagger channel cut: " << tc1 << endl;
		   return kFALSE;
		}
		else if ((tc2 < 0) || (tc2 > 352))
		{
           cout << "Invalid tagger channel cut: " << tc2 << endl;
		   return kFALSE;
		}

        cout << "Setting cut on tagger channels: " << tc1 << " to " << tc2 << endl;
		SetTC_cut(tc1,tc2);
	}
	else if(strcmp(config.c_str(), "nokey") != 0)
	{
		cout << "Tagger Channel cut not set correctly" << endl;
		return kFALSE;
	}

	cout << endl;
	return kTRUE;

}

Bool_t 	PPhysics::InitTaggerScalers()
{
	Int_t sc1, sc2;
	string config = ReadConfig("Tagger-Scalers");
	if(sscanf( config.c_str(), "%d %d\n", &sc1, &sc2) == 2)
	{
		cout << "Setting Tagger scaler channels: " << sc1 << " to " << sc2 << endl;
        SetTC_scalers(sc1,sc2);
	}
	else if(strcmp(config.c_str(), "nokey") != 0)
	{
        cout << "Tagger Channel GetScalers() not set correctly" << endl;
		return kFALSE;
	}

	cout << endl;
	return kTRUE;

}
#endif
