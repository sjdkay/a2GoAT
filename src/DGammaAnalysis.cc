#ifndef __CINT__

#include "DGammaAnalysis.h"

// Heavily altered version of PPi0Example.cc

int main(int argc, char *argv[])
{

	clock_t start, end;
	start = clock();

	// Initialise strings
	std::string configfile = "";
	std::string serverfile = "";
	std::string dir_in = "";
	std::string dir_out = "";
	std::string file_in = "";
	std::string file_out = "";
	std::string pre_in = "";
	std::string pre_out = "";

	Int_t length;
	std::string flag;

	if(argc == 1)
	{
		cout << "Please provide a config file" << endl;
		return 0;
	}
	else if(argc == 2) configfile = argv[1];
	else
	{
	    for(int i=1; i<argc; i++)
		{
			flag = argv[i];
			if(flag.find_first_of("-") == 0)
			{
				i++;
				flag.erase(0,1);
				if(strcmp(flag.c_str(), "s") == 0) serverfile = argv[i];
				else if(strcmp(flag.c_str(), "d") == 0) dir_in = argv[i];
				else if(strcmp(flag.c_str(), "D") == 0) dir_out = argv[i];
				else if(strcmp(flag.c_str(), "f") == 0) file_in = argv[i];
				else if(strcmp(flag.c_str(), "F") == 0) file_out = argv[i];
				else if(strcmp(flag.c_str(), "p") == 0) pre_in = argv[i];
				else if(strcmp(flag.c_str(), "P") == 0) pre_out = argv[i];
				else
				{
					cout << "Unknown flag " << flag << endl;
					return 0;
				}
			}
			else configfile = argv[i];
		}
	}

	// Check that config file exists:
	ifstream cfile(configfile.c_str());
	if(!cfile)
	{
		cout << "Config file '" << configfile << "' could not be found." << endl;
		return 0;
	}

	// If server file is specified, check that it exists
	if(serverfile.length() > 0)
	{
		// Check that file exists:
		ifstream sfile(serverfile.c_str());
		if(!sfile)
		{
			cout << "Server file '" << serverfile << "' could not be found" << endl;
			return 0;
		}
	}
	// If no server file is specified, allow for checking in the config file
	else serverfile = configfile;

	DGammaAnalysis* DGA = new DGammaAnalysis;

	// If unset, scan server or config file for file settings
	if(dir_in.length() == 0)
	{
		flag = DGA->ReadConfig("Input-Directory",0,(Char_t*)serverfile.c_str());	
		flag.erase(0,flag.find_first_not_of(" "));
		if(strcmp(flag.c_str(),"nokey") != 0) dir_in = flag;
	}
	
	if(dir_out.length() == 0)
	{	
		flag = DGA->ReadConfig("Output-Directory",0,(Char_t*)serverfile.c_str());	
		flag.erase(0,flag.find_first_not_of(" "));
		if(strcmp(flag.c_str(),"nokey") != 0) dir_out = flag;
	}
	
	if(file_in.length() == 0)
	{	
		flag = DGA->ReadConfig("Input-File",0,(Char_t*)serverfile.c_str());	
		flag.erase(0,flag.find_first_not_of(" "));
		if(strcmp(flag.c_str(),"nokey") != 0) file_in = flag;
	}
	
	if(file_out.length() == 0)
	{	
		flag = DGA->ReadConfig("Output-File",0,(Char_t*)serverfile.c_str());	
		flag.erase(0,flag.find_first_not_of(" "));
		if(strcmp(flag.c_str(),"nokey") != 0) file_out = flag;
	}
	
	if(pre_in.length() == 0)
	{	
		flag = DGA->ReadConfig("Input-Prefix",0,(Char_t*)serverfile.c_str());	
		flag.erase(0,flag.find_first_not_of(" "));
		if(strcmp(flag.c_str(),"nokey") != 0) pre_in = flag;
	}
	
	if(pre_out.length() == 0)
	{	
		flag = DGA->ReadConfig("Output-Prefix",0,(Char_t*)serverfile.c_str());	
		flag.erase(0,flag.find_first_not_of(" "));
		if(strcmp(flag.c_str(),"nokey") != 0) pre_out = flag;
	}
	// Finished scanning for file settings
	
	// Fix directories to include final slash if not there
	if(dir_in.find_last_of("/") != (dir_in.length()-1)) dir_in += "/";
	if(dir_out.find_last_of("/") != (dir_out.length()-1)) dir_out += "/";

	// Output user settings (Set to defaults if still unspecified)
	cout << endl << "User inputs" << endl;
	cout << "Config file:      '" << configfile << "' chosen" << endl;
	if(dir_in.length() != 0)  	cout << "Input directory:  '" << dir_in << "' chosen" << endl;
	if(dir_out.length() != 0)  	cout << "Output directory: '" << dir_out << "' chosen" << endl;
	if(file_in.length() != 0)  	cout << "Input file:       '" << file_in << "' chosen" << endl;
	if(file_out.length() != 0) 	cout << "Output file:      '" << file_out << "' chosen" << endl;
	if(pre_in.length() != 0)  	cout << "Input prefix:     '" << pre_in << "' chosen" << endl;
	else { pre_in = "GoAT"; 	cout << "Input prefix:     '" << pre_in << "' chosen by default" << endl; }
	if(pre_out.length() != 0)  	cout << "Output prefix:    '" << pre_out << "' chosen" << endl;	
	else { pre_out = "Physics"; 	cout << "Output prefix:    '" << pre_out << "' chosen by default" << endl; }
	cout << endl;
	
	// Perform full initialisation 
	if(!DGA->Init(configfile.c_str()))
	{
		cout << "ERROR: DGammaAnalysis Init failed!" << endl;
		return 0;
	}

	std::string file;
	std::string prefix;
	std::string suffix;

	Int_t files_found = 0;
	// If input file is specified, use it
	if(file_in.length() > 0)
	{
		cout << "Searching for input file(s)" << endl;
		file = file_in;
		length = file.length();
		// File should at least have '.root' at the end
		if(length >= 5)
		{
			// Add input directory to it
			file_in = dir_in+file_in;
			cout << "Input file  '" << file_in << "' chosen" << endl;

			// If output file is specified, use it
			if(file_out.length() > 0) file_out = dir_out+file_out;
			// If output file is not specified, build it
			else
			{
				// If output directory is not specified, build it
				if(dir_out.length() == 0)
				{
					prefix = file.substr(0,file.find_last_of("/")+1);
					dir_out = dir_in+prefix;
				}
				// If input prefix doesn't match, simply prepend output prefix to the file name
				if(file.find(pre_in)>file.length()) suffix = ("_"+file.substr(file.find_last_of("/")+1,length-(file.find_last_of("/")+1)));
				// If input prefix does match, switch prefixes
				else suffix = file.substr(file.find_last_of("/")+1+pre_in.length(),length-(file.find_last_of("/")+1+pre_in.length()));
				// Build output file name
				file_out = dir_out+pre_out+suffix;
			}
			
			cout << "Output file '" << file_out << "' chosen" << endl << endl;
			if(!DGA->File(file_in.c_str(), file_out.c_str())) cout << "ERROR: DGammaAnalysis failed on file " << file_in << "!" << endl;
			files_found++;
		}
	}
	// Otherwise scan input directory for matching files
	else
	{
		cout << "Searching input directory for files matching input prefix" << endl;
		cout << "Input prefix  '" << pre_in << "' chosen" << endl;
		cout << "Output prefix '" << pre_out << "' chosen" << endl;
		
		// If output directory is not specified, use the input directory
		if(dir_in.length()  == 0) dir_in = "./";
		if(dir_out.length() == 0) dir_out = dir_in;

		// Create list of files in input directory
		TSystemFile *sys_file;
		TSystemDirectory *sys_dir = new TSystemDirectory("files",dir_in.c_str());
		TList *file_list = sys_dir->GetListOfFiles();
		file_list->Sort();
		TIter file_iter(file_list);

		// Iterate over files
		while((sys_file=(TSystemFile*)file_iter()))
		{
			file = sys_file->GetName();
			length = file.length();
			// File should at least have '.root' at the end
			if(length >= 5)
			{
				//Check that prefixes and suffixes match
				prefix = file.substr(0,pre_in.length());
				suffix = file.substr(length-5,5);
				if(((strcmp(prefix.c_str(),pre_in.c_str()) == 0)) && (strcmp(suffix.c_str(),".root") == 0))
				{
					// Build input file name
					file_in = dir_in+file;
					// Build output file name
					suffix = file.substr(pre_in.length(),length-pre_in.length());
					file_out = dir_out+pre_out+suffix;					

					files_found++;
					// Run DGammaAnalysis
					if(!DGA->File(file_in.c_str(), file_out.c_str())) 
						cout << "ERROR: DGammaAnalysis failed on file " << file_in << "!" << endl;

				}
			}
		}
	}
	if (files_found == 0)
	{
		cout << "ERROR: No GoAT files found!" << endl;
		return 0;
	}

	end = clock();
	cout << "Time required for execution: "
	<< (double)(end-start)/CLOCKS_PER_SEC
	<< " seconds." << "\n\n";

	return 0;
}

DGammaAnalysis::DGammaAnalysis() 
{ 
}

DGammaAnalysis::~DGammaAnalysis()
{
}

Bool_t DGammaAnalysis::Init(const char* configfile)
{
	
  // Initialise shared pdg database
  pdgDB = TDatabasePDG::Instance();

  // Set by user in the future...
  SetTarget(1875.6);
	
  Double_t Prompt_low 	=  -15;
  Double_t Prompt_high 	=   15;
  Double_t Random_low1 	= -100;
  Double_t Random_high1 	=  -40;
  Double_t Random_low2 	=   35;
  Double_t Random_high2 	=   95;
  
  SetPromptWindow(Prompt_low, Prompt_high);
  SetRandomWindow1(Random_low1, Random_high1);
  SetRandomWindow2(Random_low2, Random_high2);
  SetPvRratio();

  return kTRUE;

}

Bool_t DGammaAnalysis::File(const char* file_in, const char* file_out)
{

  OpenGoATFile(file_in, "READ");
  OpenHistFile(file_out);
  DefineHistograms();

  cout << "Setting up tree files:" << endl;
  if(!OpenTreeParticles(GoATFile)) 	return kFALSE;
  if(!OpenTreeTagger(GoATFile))		return kFALSE;
  cout << endl;
  
  cout << "Determining valid for analysis:" << endl;	
  if(!FindValidGoATEvents())			return kFALSE;	
  cout << endl;
  
  Analyse();	
  return kTRUE;

}

void DGammaAnalysis::Analyse()
{

  Cut_CB_proton = OpenCutFile("configfiles/cuts/CB_DeltaE-E_Proton.root", "Proton");
  Cut_proton = Cut_CB_proton;
  Cut_CB_proton_4Sig = OpenCutFile("configfiles/cuts/CB_DeltaE-E_Proton_4Sigma.root", "Proton");
  Cut_proton_4Sig = Cut_CB_proton_4Sig;
  Cut_CB_proton_Full = OpenCutFile("configfiles/cuts/CB_DeltaE-E_ProtonFull.root", "Proton");
  Cut_proton_Full = Cut_CB_proton_Full;
  Cut_CB_pion = OpenCutFile("configfiles/cuts/CB_DeltaE-E_Pion.root", "Pion");
  Cut_pion = Cut_CB_pion;
  // Cut_CB_neutron = OpenCutFile("configfiles/cuts/CB_DeltaE-E_Neutron.root", "Neutron");
  // Cut_neutron = Cut_CB_neutron; // This only needs to be here if we get simulation to show where we expect neutrons
  cout << endl;

  N_P = 0, N_P2 = 0, N_Part = 0, k = 0;
  TraverseGoATEntries(); //This functions calls Reconstruct() inside it, looks at all GoAT entries
  cout << "Total particles found: " << N_Part << "  Total Protons found: " << N_P <<" ... after cuts: "<< (N_P2) << endl;
  Diff = ((N_P)-(N_P2));
  DiffD = double (Diff);
  N_PD = double (N_P);
  cout << "Percentage of protons lost due to cuts: " << (DiffD/N_PD)*100 << " %" << endl;
  Diff = ((N_Part) - (N_P2));
  DiffD = double (Diff);
  N_PD = double (N_Part);
  cout << "Percentage of particles lost due to cuts: " << (DiffD/N_PD)*100 << " %" << endl;
      
  PostReconstruction();		
  WriteHistograms();
  CloseHistFile();	

}

void DGammaAnalysis::Reconstruct() // Starts at event 0 so 0 - X events hence extra two protons
{
  
  // Potentially move several of the cuts (Theta, dE, Phi) out of the tagger and proton loop and to this point??
  // The two for loops below are themselves within a loop in TraverseGoATEntries so this should be possible, shouldn't make a difference though
  // As a sanity check move them here and check the output is the same?

	// This for loop examines the number of tagged photons for each event and loops over them all
	
	// k = 0; // Integer which can be used for counting something in the loop if desired
	
	for (Int_t i = 0; i < GoATTree_GetNParticles(); i++) {
							    	  
	  // This for loop loops over all the particles in an event, i.e. the # of protons

	  for (Int_t j = 0; j < GetNTagged(); j++) {    	   	    
	    
	    if (j == 0) N_Part++;
	    
	    // Check PDG: if both particles not a proton, continue
		    
	    if ( GoATTree_GetPDG(0) != pdgDB->GetParticle("proton")->PdgCode() ) continue; // Be careful with this cut with new GoAT files
	    if ( GoATTree_GetPDG(1) != pdgDB->GetParticle("proton")->PdgCode() ) continue; 
					
	    if (j == 0) N_P++;
				     
	    if ( ( GetPhotonBeam_E(j) )- (GoATTree_GetEk(0) ) - ( GoATTree_GetEk(1) ) > 80) continue; // This is the cut that now removes the vast majority of events
	    if ( ( GetPhotonBeam_E(j) )- (GoATTree_GetEk(0) ) - ( GoATTree_GetEk(1) ) < -40 ) continue;  
	    
	    // Cut if dE for each proton equal, with the Eg-Esum and Phi cuts this cut actually removes very few events

	    if ( GoATTree_Get_dE(0) == GoATTree_Get_dE(1) ) continue;
	    
	    // Remove events that are not approx back to back
	    
	    if ( ( abs(GoATTree_GetPhi(0) - GoATTree_GetPhi(1))) > 200 ) continue;
	    if ( ( abs(GoATTree_GetPhi(0) - GoATTree_GetPhi(1))) < 160 ) continue;
	    
	    // This cut checks that both of the particles in the event are random or they are both prompt

	    if ( IsPrompt(GetTagged_t(j) - GoATTree_GetTime(0)) && !IsPrompt(GetTagged_t(j) - GoATTree_GetTime(1))  ) continue;
	    if ( IsPrompt(GetTagged_t(j) - GoATTree_GetTime(1)) && !IsPrompt(GetTagged_t(j) - GoATTree_GetTime(0))  ) continue;
	    
	    // Do Boost for each proton 4-vector, since this is before filling histograms this COULD be used to do a cut
	    // Also if moved to be above theta could cut on CM theta instead of lab theta
	    
	    GV1 = GetGoATVector(0); 
	    GV2 = GetGoATVector(1); 
	    Gamma = TLorentzVector(0., 0., GetPhotonBeam_E(j), GetPhotonBeam_E(j)); // 4-Vector of Photon beam
	    Deut = TLorentzVector (0., 0., 0., 1875.613); // 4-Vector of Deuterium target
	    Theta1 = (GV1.Theta()) * TMath::RadToDeg();
	    Theta2 = (GV2.Theta()) * TMath::RadToDeg();
	    B = (Deut + Gamma).Beta();
	    b = TVector3(0., 0., B);
	    GV1.Boost(b);
	    GV2.Boost(b);
	    Theta1B = (GV1.Theta()) * TMath::RadToDeg();
	    Theta2B = (GV2.Theta()) * TMath::RadToDeg();
	    
	    // Fill timing histogram (all PDG matching proton)
	    
	    FillTimePDG( pdgDB->GetParticle("proton")->PdgCode(), time_proton );
	    
	    FillMissingMassPair( i, j, prompt_proton, random_proton );
	    
	    // Look at prompt photons for each proton
	    
	    if ( IsPrompt( GetTagged_t(j) - GoATTree_GetTime(i) ) ) {
	      
	      PEp -> Fill( GoATTree_GetEk(i) );
	      PTheta -> Fill( GoATTree_GetTheta(i) ); 
	      EpdE -> Fill(GoATTree_GetEk(i), GoATTree_Get_dE(i));
	      N_P2++;
	      
	      if(Cut_proton->IsInside(GoATTree_GetEk(i), GoATTree_Get_dE(i))) { // Use loops like this to check which cut particle is present in
		k++; // Idea is that one proton is in our more accurate banana and the other is somewhere else
	      }
	      
	      if ((i % 2) != 0) {
		
		P1Calc = (Gamma + Deut) - GetGoATVector(i);
		P1CalcTheta = (P2Calc.Theta()) * TMath::RadToDeg();
		P1ThetaDiff = abs( Theta1 - P1CalcTheta );
		P2Calc = (Gamma + Deut) - GetGoATVector(i-1); // Calculate 4-vector of second particle using first
		P2CalcTheta = (P2Calc.Theta()) * TMath::RadToDeg(); // Conservation of 4-momentum, if other particle really is a proton we should get same angle
		P2ThetaDiff = abs( Theta2 - P2CalcTheta );
		
		// l1 = 56.2/(tan(GoATTree_GetTheta(0) * TMath::DegToRad()));
		// ThetaNReal1 = (atan(56.2/( l1 + GoATTree_GetWC_Vertex_Z(0)))) *  TMath::RadToDeg();
		// l2 = 56.2/(tan(GoATTree_GetTheta(1) * TMath::DegToRad()));
		// ThetaNReal2 = (atan(56.2/( l2 + GoATTree_GetWC_Vertex_Z(1)))) *  TMath::RadToDeg();
		// cout<< l1 << "    " << ThetaNReal1 << "     " << Theta1 << "    " << "    " << l2 << "    " << ThetaNReal2 << "    " << Theta2  << endl;
		// cout << j << "    " << Theta1 << "    " << Theta2 << "    " << GetPhotonBeam_E(j) << "    " << GetTagged_t(j) << endl;
		
		PEgPrompt -> Fill( GetPhotonBeam_E(j) );
		Eg_EpsumPrompt -> Fill( ( GetPhotonBeam_E(j) )- (GoATTree_GetEk(0) ) - ( GoATTree_GetEk(1) ) );
		PEpTot -> Fill( ( GoATTree_GetEk(0) + GoATTree_GetEk(1) ) ); 
		PThetaCMPrompt -> Fill(Theta1B);
		PThetaCMPrompt -> Fill(Theta2B);
		P2CDiff -> Fill(P2ThetaDiff);
		
	      }
	    }
	    
	    // Look at random photons for each proton
	    
	    else if ( IsRandom( GetTagged_t(j) - GoATTree_GetTime(i) ) ) {			   			   
	      
	      if ((i % 2) != 0) {
		
		Eg_EpsumRandom -> Fill( (GetPhotonBeam_E(j))- (GoATTree_GetEk(0)) - (GoATTree_GetEk(1)));
		PEgRandom -> Fill( GetPhotonBeam_E(j) );
		PThetaCMRandom -> Fill(Theta1B);
		PThetaCMRandom -> Fill(Theta2B);
		
	      }											     
	    }			 
	  }  
	}

        // cout << endl; // Insert this to separate each event out when printing some info from each event
	if (GetGoATEvent() != 0) {
	  if(GetGoATEvent() % 1000 == 0) cout << "Event: "<< GetGoATEvent() << " Total Protons found: " << N_P << " ... after cuts: "<<N_P2 << endl;
	}
}

void DGammaAnalysis::PostReconstruction()
{
 
     cout << "Performing post reconstruction." << endl;

     RandomSubtraction( prompt_proton, random_proton, proton );
     RandomSubtraction( Eg_EpsumPrompt, Eg_EpsumRandom, Eg_Epsum );
     RandomSubtraction( PEgPrompt, PEgRandom, PEg );
     RandomSubtraction( PThetaCMPrompt, PThetaCMRandom, PThetaCM );
		
     ShowTimeCuts( time_proton, time_proton_cuts );

}

void DGammaAnalysis::DefineHistograms()
{
 	gROOT->cd();

	time_proton = new TH1D( "time_proton", "time_proton", 1000, -500, 500 );
	time_proton_cuts = new TH1D( "time_proton_cuts", "time_proton_cuts", 1000, -500, 500 );
	
	Eg_Epsum = new TH1D( "Eg - Epsum", "Eg - Epsum", 100, -60, 100);
	PEp = new TH1D( "P_Ep", "P_Ep", 100, 0, 500 );
	PEg = new TH1D( "photonbeam_E", "photonbeam_E", 100, 100, 900 );
	PEpTot = new TH1D( "P_Ep_Total", "P_Ep_Total", 300, 0, 900 );	
	PTheta = new TH1D( "P_Theta", "P_Theta", 150, 0, 160 );
	PThetaCM = new TH1D( "P_ThetaCM", "P_ThetaCM", 150, 0, 160 );
	EpdE = new TH2D("E_dE", "E_dE", 150, 0, 500, 150, 0, 8);
	// EpdEProton = new TH2D("E_dE_Proton", "E_dE_Proton", 150, 0, 500, 150, 0, 8);
	// EpdEPion = new TH2D("E_dE_Pion", "E_dE_Pion", 150, 0, 500, 150, 0, 8);
	proton = new TH1D( "proton", "proton", 1500, 0, 1500 );	
	P2CDiff = new TH1D( "P2CDiff", "P2CDiff", 100, 0, 130 );

	Eg_EpsumPrompt = new TH1D( "Eg - Epsum_Prompt", "Eg - Epsum_Prompt", 100, -60, 100 );
	Eg_EpsumRandom = new TH1D( "Eg - Epsum_Random", "Eg - Epsum_Random", 100, -60, 100 );
	PEgPrompt = new TH1D( "photonbeam_E_Prompt", "photonbeam_E_Prompt", 100, 100, 900 );
	PEgRandom = new TH1D( "photonbeam_E_Random", "photonbeam_E_Random", 100, 100, 900 );	
	PThetaCMPrompt = new TH1D( "P_ThetaCM_Prompt", "P_ThetaCM_Prompt", 150, 0, 160 );
	PThetaCMRandom = new TH1D( "P_ThetaCM_Random", "P_ThetaCM_Random", 150, 0, 160 );
	prompt_proton = new TH1D( "prompt_proton", "prompt_proton", 1500, 0, 1500 );
	random_proton = new TH1D( "random_proton", "random_proton",1500, 0, 1500 );
		
}

Bool_t DGammaAnalysis::WriteHistograms(TFile* pfile)
{
	
  cout << "Writing histograms." << endl;
		
	if(!pfile) return kFALSE;
	pfile->cd();

	gROOT->GetList()->Write();
	gROOT->GetList()->Delete();
		
	return kTRUE;
	
}

TCutG*	DGammaAnalysis::OpenCutFile(Char_t* filename, Char_t* cutname)
{
	CutFile = new TFile(filename, "READ");

    if( !CutFile || !CutFile->IsOpen() ) {
        cerr << "Can't open cut file: " << filename << endl;
        throw false;
    }


    // Try to find a TCutG with the name we want
    // GetObject checks the type to be TCutG,
    // see http://root.cern.ch/root/html534/TDirectory.html#TDirectory:GetObject
    CutFile->GetObject(cutname, Cut);

    if( !Cut ) {
        cerr << "Could not find a TCutG with the name " << cutname << " in " << filename << endl;
        throw false;
    }
	
	TCutG* Cut_clone = Cut;
	CutFile->Close();

	cout << "cut file " << filename << 
	 " opened (Cut-name = " << cutname << ")"<< endl;
	return Cut_clone;
	
}


#endif

