
#include "GSort.h"


using namespace std;



GSort::GSort() :
     SP_n(0),
     SP_type(0),
     SP_condition(0),
     SP_theta_min(0),
     SP_theta_max(0),
     SN_n(0),
     SN_type(0),
     SN_condition(0),
     SN_theta_min(0),
     SN_theta_max(0)     
{
	SP_n 		= new Int_t[50];
    SP_type		= new GTreeParticle*[50];
    SP_condition= new Sort_Condition[50];
	SP_theta_min= new Double_t[50];
	SP_theta_max= new Double_t[50];

	SN_n 		= new Int_t[50];
    SN_type		= new Int_t[50];
    SN_condition= new Sort_Condition[50];
	SN_theta_min= new Double_t[50];
	SN_theta_max= new Double_t[50];
}

GSort::~GSort()
{
}

Bool_t	GSort::Init()
{
	// Cut on CB Energy Sum (raw data cut)
    string config = ReadConfig("SortRaw-CBEnergySum");
	if (strcmp(config.c_str(), "nokey") == 0) SortRawCBESum = 0;
	else if( sscanf(config.c_str(),"%lf %s\n",&SR_CBESum, string_in1) == 2 )
	{
		CheckConfigCondition(string_in1, &SR_CBESum_condition, string_out1);
		if(SR_CBESum_condition == -1) SortRawCBESum = 0;
		else SortRawCBESum = 1;
		
		if(SortRawCBESum == 1)
			cout << "Sort: Crystal Ball Energy Sum " 
				 << SR_CBESum << " "<< string_out1 << endl;
		else 
		{
			cout << "SortRaw-CBEnergySum cut set improperly" <<endl;
			return kFALSE;
		 }
		 cout << endl;
	}	
	else {
		SortRawCBESum = 0;
	}	

	// Cut on raw number of particle tracks (raw data cut)
    config = ReadConfig("SortRaw-NTracks");
    if (strcmp(config.c_str(), "nokey") == 0) SortNTracks = 0;
	else if( sscanf( config.c_str(), "%d %s %d %s %d %s\n", 
        &SR_nTracks_total, 	string_in1,
        &SR_nTracks_CB,   	string_in2,
        &SR_nTracks_TAPS,  	string_in3) == 6 )
	{		
        SortNTracks = 1;

		// Total Number of particle tracks
        CheckConfigCondition(string_in1, &SR_nTracks_total_condition, string_out1);
        if(SR_nTracks_total_condition == -1) SortNTracks = 0;
		
		// Number of particle tracks in CB
        CheckConfigCondition(string_in2, &SR_nTracks_CB_condition, string_out2);
        if(SR_nTracks_CB_condition == -1) SortNTracks = 0;

		// Number of particle tracks in TAPS
        CheckConfigCondition(string_in3, &SR_nTracks_TAPS_condition, string_out3);
        if(SR_nTracks_TAPS_condition == -1) SortNTracks = 0;

        if(SortNTracks == 1)
		{
            cout << "Sort: Total # of Tracks before reconstruction   "
                 << SR_nTracks_total << " "<< string_out1 << endl;
            cout << "Sort: # of Tracks before reconstruction in CB   "
                 << SR_nTracks_CB << " "<< string_out2 << endl;
            cout << "Sort: # of Tracks before reconstruction in TAPS "
                 << SR_nTracks_TAPS << " "<< string_out3 << endl;
		 }
		 else
		 {
            cout << "SortRaw-NTracks cut set improperly" <<endl;
			return kFALSE;
		 }
		 cout << endl;
	}
	else {
        SortNTracks = 0;
	}
	
	// Cut on reconstructed number of particles
	config = ReadConfig("Sort-NParticles");
	if (strcmp(config.c_str(), "nokey") == 0) SortNParticles = 0;
	else if( sscanf( config.c_str(), "%d %s\n", 
		&S_nParticles, string_in1) == 2 )
	{

		CheckConfigCondition(string_in1, &S_nParticles_condition, string_out1);
		if(S_nParticles_condition == -1) SortNParticles = 0;
		else SortNParticles = 1;
		
		if(SortNParticles == 1)
			cout << "Sort: Total # of Particles after reconstruction    " 
				 << S_nParticles << " "<< string_out1 << endl;
		else 
		{
			cout << "Sort-NParticles cut set improperly" <<endl;
			cout << "Continuing with this cut turned off!" << endl;
		 }
	}	
	else {
		SortNParticles = 0;
	}	
	
	// Cut on individual particle properties
	Int_t instance = 0;
	n_cut_SP = 0;
	n_cut_SN = 0;
	do
	{
        char            particleName[256], cond[256];
        Int_t           num;
        Sort_Condition  condition;
        Double_t        th_min, th_max;

		config = ReadConfig("Sort-Particle",instance);
        if( sscanf( config.c_str(), "%s %d %s %lf %lf\n", particleName, &num, cond, &th_min, &th_max) == 5 )
		{			
			CheckConfigCondition(cond, &condition, string_out1);
			
            if((strcmp(particleName,"charged") == 0) || (strcmp(particleName,"neutral") == 0))
			{
				SN_n[n_cut_SN] = num;
                if (strcmp(particleName,"charged") == 0) SN_type[n_cut_SN] = 1;
				else SN_type[n_cut_SN] = 0;
				SN_condition[n_cut_SN] = condition;
				SN_theta_min[n_cut_SN] = th_min;
				SN_theta_max[n_cut_SN] = th_max;
				n_cut_SN++;
				
                cout << "Sort: " << particleName << " (" << num << cond <<","<<th_min<<"," << th_max<<")" << endl;
			}	
            else// if(pdgDB->GetParticle(particleName) != 0x0) // Check for particle in pdg database
            {
                Int_t     i       =0;
                while(i<GetTreeList().GetEntries())
                {
                    if(strcmp(GetTreeList()[i]->GetName(),particleName) == 0)
                    {
                        SP_type[n_cut_SP] 	   = (GTreeParticle*)GetTreeList()[i];
                        break;
                    }
                    i++;
                }
                if(i<GetTreeList().GetEntries())
                {
                    SP_n[n_cut_SP] 		   = num;
                    SP_condition[n_cut_SP] = condition;
                    SP_theta_min[n_cut_SP] = th_min;
                    SP_theta_max[n_cut_SP] = th_max;
                    n_cut_SP++;
                    cout << "Sort: " << particleName << " (" << num << cond <<","<<th_min<<"," << th_max<<")" << endl;
                }
                else
                {
                    cout << "tree for particle type "  << particleName  << " does not exist." << endl;
                    return kFALSE;
                }
			}
            /*
            else
            {
                cout << endl << "ERROR unknown particle type ("  << particleName
					 << ") set by Sort-Particle." << endl << endl;
				return kFALSE;
			}
            */
		}	
		else if (strcmp(config.c_str(), "nokey") != 0)
		{
			cout << "Sort-Particle cut was set improperly" << endl;
			return kFALSE;		
		}	
		instance++;
	} while (strcmp(config.c_str(), "nokey") != 0);

	cout << endl;
	return kTRUE;
}

void	GSort::Reconstruct()
{
}

Bool_t GSort::SortAnalyseEvent()
{

	// Sort on raw variables before analysis (increases speed)
    if(SortNTracks == 1)
	{
        switch (SR_nTracks_total_condition) 	// Total number of tracks
		{
            case Condion_EqualOrMore:
                if (GetTracks()->GetNTracks() < SR_nTracks_total) 	return kFALSE;
				break;
            case Condion_EqualOrLess:
                if (GetTracks()->GetNTracks() > SR_nTracks_total) 	return kFALSE;
				break;
            case Condion_Equal:
                if (GetTracks()->GetNTracks() != SR_nTracks_total) 	return kFALSE;
				break;
            case Condion_NONE:
                return kFALSE;
		}
		
        switch (SR_nTracks_CB_condition) 	// Number of tracks in CB
		{
            case Condion_EqualOrMore:
                if (GetTracks()->GetNCB() < SR_nTracks_CB) 			return kFALSE;
				break;
            case Condion_EqualOrLess:
                if (GetTracks()->GetNCB() > SR_nTracks_CB) 			return kFALSE;
				break;
            case Condion_Equal:
                if (GetTracks()->GetNCB() != SR_nTracks_CB) 			return kFALSE;
				break;
            case Condion_NONE:
                return kFALSE;
		}
		
        switch (SR_nTracks_TAPS_condition) 	// Number of tracks in TAPS
		{
            case Condion_EqualOrMore:
                if (GetTracks()->GetNTAPS() < SR_nTracks_TAPS) 		return kFALSE;
				break;
            case Condion_EqualOrLess:
                if (GetTracks()->GetNTAPS() > SR_nTracks_TAPS) 		return kFALSE;
				break;
            case Condion_Equal:
                if (GetTracks()->GetNTAPS() != SR_nTracks_TAPS)		return kFALSE;
				break;
            case Condion_NONE:
                return kFALSE;
		}
	}
	
	if(SortRawCBESum == 1)
	{
		switch (SR_CBESum_condition) 	// Crystal Ball Energy Sum
		{
            case Condion_EqualOrMore:
                if (GetTrigger()->GetEnergySum() < SR_CBESum) 				return kFALSE;
				break;
            case Condion_EqualOrLess:
                if (GetTrigger()->GetEnergySum() > SR_CBESum) 				return kFALSE;
				break;
            case Condion_Equal:
                if (GetTrigger()->GetEnergySum() != SR_CBESum)				return kFALSE;
				break;
            case Condion_NONE:
                return kFALSE;
		}
	}	
	
	return kTRUE;
}


Bool_t GSort::SortFillEvent()
{
    if(SortNParticles == 1)
	{
		switch (S_nParticles_condition) // Number of reconstructed part
		{
            case Condion_EqualOrMore:
                if(GetNReconstructed() < S_nParticles)  return kFALSE;
				break;
            case Condion_EqualOrLess:
                if(GetNReconstructed() > S_nParticles)  return kFALSE;
				break;
            case Condion_Equal:
                if(GetNReconstructed() != S_nParticles) return kFALSE;
				break;
            case Condion_NONE:
                return kFALSE;
		}
	}

    for(Int_t i = 0; i < n_cut_SN; i++)
	{
		if(!SortOnNeutrality(SN_type[i],	
							SN_n[i], 
							SN_condition[i], 
							SN_theta_min[i], 
							SN_theta_max[i]))		return kFALSE;
    }

	for(Int_t i = 0; i < n_cut_SP; i++)
	{
        if(!SortOnParticle( *SP_type[i],
							SP_n[i], 
							SP_condition[i], 
							SP_theta_min[i], 
							SP_theta_max[i]))		return kFALSE;
	}		
			
    // No cut failed, so return TRUE
	return kTRUE;

}


Bool_t	GSort::SortOnParticle(const GTreeParticle &tree, Int_t Num, Int_t cond, Double_t ThetaMin, Double_t ThetaMax)
{
	Int_t NumberFound = 0;
	
    for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
        // Check theta limits
        if ((tree.GetTheta(i) < ThetaMin) ||
            (tree.GetTheta(i) >= ThetaMax))
			return kFALSE;
        NumberFound++;
    }

    switch (cond)
	{
        case Condion_EqualOrMore:
			if (NumberFound < Num) 	return kFALSE;
			break;
        case Condion_EqualOrLess:
			if (NumberFound > Num) 	return kFALSE;
			break;
        case Condion_Equal:
			if (NumberFound != Num)	return kFALSE;
			break;
        case Condion_NONE:
            return kFALSE;
	}

	return kTRUE;
}

Bool_t	GSort::SortOnNeutrality(Bool_t charge, Int_t Num, Sort_Condition cond, Double_t ThetaMin, Double_t ThetaMax)
{
	Int_t NumberFound = 0;
	
    if(charge)
	{
        for (Int_t i = 0; i <GetElectrons()->GetNParticles(); i++)
        {
            //Check theta limits
            if ((GetElectrons()->GetTheta(i) < ThetaMin) ||
                (GetElectrons()->GetTheta(i) >= ThetaMax))
                continue;
            NumberFound++;
        }
        for (Int_t i = 0; i <GetChargedPions()->GetNParticles(); i++)
        {
            //Check theta limits
            if ((GetChargedPions()->GetTheta(i) < ThetaMin) ||
                (GetChargedPions()->GetTheta(i) >= ThetaMax))
                continue;
            NumberFound++;
        }
        for (Int_t i = 0; i <GetProtons()->GetNParticles(); i++)
        {
            //Check theta limits
            if ((GetProtons()->GetTheta(i) < ThetaMin) ||
                (GetProtons()->GetTheta(i) >= ThetaMax))
                continue;
            NumberFound++;
        }
        for (Int_t i = 0; i <GetRootinos()->GetNParticles(); i++)
        {
            //Check theta limits
            if ((GetRootinos()->GetTheta(i) < ThetaMin) ||
                (GetRootinos()->GetTheta(i) >= ThetaMax))
                continue;
            NumberFound++;
        }
    }
    else
    {
        for (Int_t i = 0; i < GetPhotons()->GetNParticles(); i++)
        {
            //Check theta limits
            if ((GetPhotons()->GetTheta(i) < ThetaMin) ||
                (GetPhotons()->GetTheta(i) >= ThetaMax))
                continue;
            NumberFound++;
        }
        for (Int_t i = 0; i < GetNeutrons()->GetNParticles(); i++)
        {
            //Check theta limits
            if ((GetNeutrons()->GetTheta(i) < ThetaMin) ||
                (GetNeutrons()->GetTheta(i) >= ThetaMax))
                continue;
            NumberFound++;
        }
	}	

	switch (cond)
	{
        case Condion_EqualOrMore:
			if (NumberFound < Num) 	return kFALSE;
			break;
        case Condion_EqualOrLess:
			if (NumberFound > Num) 	return kFALSE;
			break;
        case Condion_Equal:
			if (NumberFound != Num)	return kFALSE;
			break;
        case Condion_NONE:
            return kFALSE;
	}

	return kTRUE;
}

void GSort::CheckConfigCondition(char string[], Sort_Condition *condition, std::string& string_out)
{
		if(strcmp(string,"+") == 0) 
		{
			string_out = "or more ";
            (*condition) = Condion_EqualOrMore;
		}
		else if (strcmp(string,"-") == 0) 
		{
			string_out = "or less ";
            (*condition) = Condion_EqualOrLess;
		}		
		else if (strcmp(string,"=") == 0) 
		{
			string_out = "exactly ";
            (*condition) = Condion_Equal;
		}
		else
		{
			string_out = "error ";
            (*condition) = Condion_NONE;
		}

}
