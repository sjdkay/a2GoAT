#ifndef __CINT__

#include "DGammaAnalysis.h"

// Heavily altered version of PPi0Example.cc, unsure if "Reconstruction" and "Post Reconstruction" steps needed in this case
// Also unsure as to what exactly the "Analyse" step is doing

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
	else { pre_out = "Proton"; 	cout << "Output prefix:    '" << pre_out << "' chosen by default" << endl; }
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

Bool_t	DGammaAnalysis::Init(const char* configfile)
{
	// Initialise shared pdg database
	pdgDB = TDatabasePDG::Instance();

	// Set by user in the future...
	SetTarget(1875.6);
	
	Double_t Prompt_low 	=  -20;
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

Bool_t	DGammaAnalysis::File(const char* file_in, const char* file_out)
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

void	DGammaAnalysis::Analyse()
{

  TraverseGoATEntries(); //This just looks through all GoAT events, how does it know they're Pi0's/protons?
  	cout << "Total Protons found: " << N_P << endl << endl;
	
  	PostReconstruction();		
  	WriteHistograms();
  	CloseHistFile();	

}

//Change this to look for protons from PDG? Don't want time though want energy/angle. Also don't need to reconstruct protons? Is this phase needed?
void	DGammaAnalysis::Reconstruct()
{
  	if(GetGoATEvent() == 0) N_P = 0;
		else if(GetGoATEvent() % 1000 == 0) cout << "Event: "<< GetGoATEvent() << " Total Protons found: " << N_P << endl;

		for (Int_t i = 0; i < GoATTree_GetNParticles(); i++)
	{
		// Count total protons
	         	if(GoATTree_GetPDG(i) == pdgDB->GetParticle("proton")->PdgCode()) 	N_P++;
		
		//Check PDG: Not proton, continue
			if (GoATTree_GetPDG(i) != pdgDB->GetParticle("proton")->PdgCode()) continue; 
		
		// Check Charge: Not +1, continue
			if (GoATTree_GetCharge(i) != 1) continue;
			
			//Fill some histograms with events
			
			cout<<GoATTree_GetEk(i)<<"   "<<GoATTree_GetTheta(i)<<endl; //Use these to get parameters we want
			// Now need to fill the correct histograms using these values I.e. fill PEk e.t.c. defined below
		
		// Fill missing mass for particle i
		//	FillMissingMass(i, MM_prompt_pi0_n_2g, MM_prompt_pi0_n_2g);
	}
}

void  DGammaAnalysis::PostReconstruction()
{
  //    cout << "Performing post reconstruction." << endl;

  //	RandomSubtraction(MM_prompt_pi0,MM_random_pi0, MM_pi0);		
  //	RandomSubtraction(MM_prompt_pi0_n_2g,MM_random_pi0_n_2g, MM_pi0_n_2g);	
		
  //	ShowTimeCuts(time_pi0, time_pi0_cuts);

}

void	DGammaAnalysis::DefineHistograms()
{
 	gROOT->cd();
       
	PEk = new TH1D("P_Ek", "P_Ek", 10000, 0, 1000);
	PTheta = new TH1D("P_Theta", "P_Theta", 10000, -180, 180);
	PEg = new TH1D("photonbeam_E", "photonbeam_E", 100000, 0, 1500);
	EpTp = new TH2D("EpTp", "EpTp", 2000 ,0, 1000, 500, -180, 180);
   	
}

Bool_t 	DGammaAnalysis::WriteHistograms(TFile* pfile)
{
	cout << "Writing histograms." << endl;
		
	if(!pfile) return kFALSE;
	pfile->cd();

	gROOT->GetList()->Write();
	gROOT->GetList()->Delete();
		
	return kTRUE;
}


#endif

//So far should open and check files and then prepare output files to put stuff in, NOT certain that PvRratio needs to be done.
//Now need to try and get script to do what we want - e.g. Get a plot of Ek vs photonbeam_E, Ek vs Theta, need to define classes correctly for this
//Need to change define histograms to the ones we actually want
//Then need to write them too
