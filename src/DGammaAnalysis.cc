#ifndef __CINT__

#include "DGammaAnalysis.h"

// Heavily altered version of PPi0Example.cc
// Unsure as to what exactly the "Analyse" step is doing

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

  TraverseGoATEntries(); //This functions calls Reconstruct() inside it, looks at all GoAT entries
  cout << "Total Protons found: " << N_P <<" ... after cuts: "<< (N_P2) << endl;
   Diff = ((N_P)-(N_P2));
   DiffD = double (Diff);
   N_PD = double (N_P);
   cout << "Percentage of events lost after cuts: " << (DiffD/N_PD)*100 << " %" << endl;
	
  	PostReconstruction();		
  	WriteHistograms();
  	CloseHistFile();	

}

void	DGammaAnalysis::Reconstruct() // Starts at event 0 so 0 - X events hence extra two protons
{
  if(GetGoATEvent() == 0) N_P = 0, N_P2 = 0;
  else if(GetGoATEvent() % 1000 == 0) cout << "Event: "<< GetGoATEvent() << " Total Protons found: " << N_P << " ... after cuts: "<<N_P2 << endl;

  	// Fill timing histogram (all PDG matching proton)
	FillTimePDG(pdgDB->GetParticle("proton")->PdgCode(), time_proton);

	// Fill missing mass (all PDG matching proton)
	MissingMassPDG(pdgDB->GetParticle("proton")->PdgCode(), prompt_proton, random_proton);
   
	// This for loop examines the number of tagged photons for each event and loops over them all
	
	// k = 0;
	
	for (Int_t j = 0; j < GetNTagged(); j++) {
	  	  
	  // if (k > 1) continue;
	  // if (IsPrompt(GetTagged_t(j), -20, 15) == kTRUE) k++;
							    	  
	  // This for loop loops over all the particles in an event, i.e. the # of protons

	    for (Int_t i = 0; i < GoATTree_GetNParticles(); i++)
		  { 
   
		      // Check PDG: Not proton, continue
		    
		      if ( GoATTree_GetPDG(i) != pdgDB->GetParticle("proton")->PdgCode() ) continue; 
			
		      // Check Charge: Not +1, continue
		      
			if ( GoATTree_GetCharge(i) != 1 ) continue;
			
			if (j == 0) N_P++;
			
			// Cut if dE for each proton equal
			
			if ( GoATTree_Get_dE(0) == GoATTree_Get_dE(1) ) continue;
			
			// Cuts to remove sections of theta
			
			if ( GoATTree_GetTheta(0) < 90 ){
			  if ( GoATTree_GetTheta(1) < 90 ) continue;
			}
			
			 // if (GoATTree_GetTheta(0) < 90){
			 // if (GoATTree_GetTheta(1) > 90) continue;
			 // }
			
			 // if (GoATTree_GetTheta(0) > 90){
			 // if (GoATTree_GetTheta(1) < 90) continue;
			 // }
			
			if (GoATTree_GetTheta(0) > 90){
			   if (GoATTree_GetTheta(1) > 90) continue;
			 }
			
			// Remove events that are not approx back to back
			
			if ( (abs(GoATTree_GetPhi(0) - GoATTree_GetPhi(1))) > 200 ) continue;
			if ( (abs(GoATTree_GetPhi(0) - GoATTree_GetPhi(1))) < 160 ) continue;
			
			// Remove events if Eg - Epsum out of a certain window

			// if ((EgPrompt - (GoATTree_GetEk(0)) - (GoATTree_GetEk(1))) > 200) continue;
			// if ((EgPrompt - (GoATTree_GetEk(0)) - (GoATTree_GetEk(1))) < -100) continue;

			 // Look at prompt photons for each proton

			 if ( IsPrompt(GetTagged_t(j), -20, 15) == kTRUE ){
			   
			   PEpPrompt->Fill( GoATTree_GetEk(i) );
			   N_P2++;

			   if ((i % 2) != 0) {
			     PEgPrompt->Fill( GetPhotonBeam_E(j) );
			     Eg_EpsumPrompt->Fill( ( GetPhotonBeam_E(j) )- (GoATTree_GetEk(0) ) - ( GoATTree_GetEk(1) ) );
			     PEpTotPrompt->Fill( (GoATTree_GetEk(0) + GoATTree_GetEk(1) ) );
			   }
			   //cout << GetPhotonBeam_E(j)<<"   "<< j <<endl;

			 }

			 // look at random photons for each proton

			 else if (IsRandom(GetTagged_t(j), -100, -40, 35, 95) == kTRUE){			   
			   
			   PEpRandom->Fill( GoATTree_GetEk(i) );
			   
			   if ((i % 2) != 0) {
			     Eg_EpsumRandom->Fill( (GetPhotonBeam_E(j))- (GoATTree_GetEk(0)) - (GoATTree_GetEk(1)));
			     PEgRandom->Fill( GetPhotonBeam_E(j) );
			     PEpTotRandom->Fill( (GoATTree_GetEk(0) + GoATTree_GetEk(1) ) );
			   }
											     
			 }

			 // PEp->Fill(GoATTree_GetEk(i));
			 // PTheta->Fill(GoATTree_GetTheta(i));
			 // PPhi->Fill(GoATTree_GetPhi(i));
			 // EpEg->Fill(GoATTree_GetEk(i),GetPhotonBeam_E(i));
			 // EpTp->Fill(GoATTree_GetEk(i), GoATTree_GetTheta(i));
			 // PVX->Fill(GoATTree_GetWC_Vertex_X(i));
			 // PVY->Fill(GoATTree_GetWC_Vertex_Y(i));
			 // PVZ->Fill(GoATTree_GetWC_Vertex_Z(i));
			 // nTAPS->Fill(GetBaF2_PbWO4_Hits(i));
			 // EpdE->Fill(GoATTree_GetEk(i), GoATTree_Get_dE(i));
			 // TpPp->Fill(GoATTree_GetTheta(i), GoATTree_GetPhi(i)); 
			 // TpdE->Fill(GoATTree_GetTheta(i), GoATTree_Get_dE(i));
			 // EpPVZ->Fill(GoATTree_GetEk(i), GoATTree_GetWC_Vertex_Z(i));
			 
			 // if(GetPhotonBeam_E(i) > 450) EpdE_HighEg->Fill(GoATTree_GetEk(i), GoATTree_Get_dE(i)); // 7/8 events high Eg
			 // else if (GetPhotonBeam_E(i) < 450) EpdE_LowEg->Fill(GoATTree_GetEk(i), GoATTree_Get_dE(i));
			 
			 // If loop only selects odd numbers, i.e the second proton, first is particle 0, second is particle 1
			 // This loop is the one used to make comparison hisotgrams
			 
			  if ((i % 2) != 0){
			   
			   // PEg->Fill(GetPhotonBeam_E(j));
			   
			   // PEpTot->Fill(((GoATTree_GetEk(i)) + (GoATTree_GetEk(i-1))));
			   // PPhiDiff->Fill(abs(GoATTree_GetPhi(i) - GoATTree_GetPhi(i-1)));
			   // dE1_dE2->Fill(GoATTree_Get_dE(i-1), GoATTree_Get_dE(i));
			   // Ep1_Ep2->Fill(GoATTree_GetEk(i-1), GoATTree_GetEk(i));
			   // PVZ1_PVZ2->Fill(GoATTree_GetWC_Vertex_Z(i-1), GoATTree_GetWC_Vertex_Z(i));
			   // PTheta1_PTheta2->Fill(GoATTree_GetTheta(i-1), GoATTree_GetTheta(i));
			   // Ep1dE1->Fill(GoATTree_GetEk(i-1), GoATTree_Get_dE(i-1));			  
			   // Ep2dE2->Fill(GoATTree_GetEk(i), GoATTree_Get_dE(i));
			   
			   // cout<< EgPrompt << "    " << GoATTree_GetEk(0) << "    " << GoATTree_GetEk(1) << "    " << (EgPrompt - (GoATTree_GetEk(i)) - (GoATTree_GetEk(i-1))) << endl;
			   
			    GV1 = GetGoATVector(i-1);
			    GV2 = GetGoATVector(i);
			    sum = GV1 + GV2;
			   
			   // cout<<G1(0)<<"   "<<G2(0)<<"    "<<sum(0)<<endl;
			   
			   // cout<<GoATTree_GetEk(i-1)<<"    "<<GoATTree_GetEk(i)<<"   "<<GoATTree_Get_dE(i-1)<<"    "<<GoATTree_Get_dE(i)<<"    "<<GoATTree_GetWC_Vertex_Z(i-1)<<"    "<<GoATTree_GetWC_Vertex_Z(i)<<endl;
			   
			  }
			 
		  }  
	}	    
}

void  DGammaAnalysis::PostReconstruction()
{
 
     cout << "Performing post reconstruction." << endl;

     RandomSubtraction( prompt_proton, random_proton, proton );
     RandomSubtraction( Eg_EpsumPrompt, Eg_EpsumRandom, Eg_Epsum );
     RandomSubtraction( PEpPrompt, PEpRandom, PEp );
     RandomSubtraction( PEgPrompt, PEgRandom, PEg );
     RandomSubtraction( PEpTotPrompt, PEpTotRandom, PEpTot );
		
     ShowTimeCuts( time_proton, time_proton_cuts );

}



void	DGammaAnalysis::DefineHistograms()
{
 	gROOT->cd();

	time_proton = new TH1D( "time_proton", "time_proton", 1000, -500, 500 );
	time_proton_cuts = new TH1D( "time_proton_cuts", "time_proton_cuts", 1000, -500, 500 );
	prompt_proton 	= new TH1D( "prompt_proton", "prompt_proton", 1500, 0, 1500 );
	random_proton 	= new TH1D( "random_proton", "random_proton",1500, 0, 1500 );
	proton	      	= new TH1D( "proton", "proton", 1500, 0, 1500 );
	
	Eg_Epsum = new TH1D( "Eg - Epsum", "Eg - Epsum", 100, -500, 800);
	Eg_EpsumPrompt = new TH1D( "Eg - Epsum_Prompt", "Eg - Epsum_Prompt", 100, -500, 800 );
	Eg_EpsumRandom = new TH1D( "Eg - Epsum_Random", "Eg - Epsum_Random", 100, -500, 800 );
	PEp = new TH1D( "P_Ep", "P_Ep", 100, 0, 500 );
	PEpPrompt = new TH1D( "P_Ep_Prompt", "P_Ep_Prompt", 100, 0, 500 );
	PEpRandom = new TH1D( "P_Ep_Random", "P_Ep_Random", 100, 0, 500 );
	PEg = new TH1D( "photonbeam_E", "photonbeam_E", 100, 100, 900 );
	PEgPrompt = new TH1D( "photonbeam_E_Prompt", "photonbeam_E_Prompt", 100, 100, 900 );
	PEgRandom = new TH1D( "photonbeam_E_Random", "photonbeam_E_Random", 100, 100, 900 );	
	PEpTot = new TH1D( "P_Ep_Total", "P_Ep_Total", 300, 0, 900 );
	PEpTotPrompt = new TH1D( "P_Ep_Total_Prompt", "P_Ep_Total_Prompt", 300, 0, 900 );
	PEpTotRandom = new TH1D( "P_Ep_Total_Random", "P_Ep_Total_Random", 300, 0, 900 );

	// PTheta = new TH1D("P_Theta", "P_Theta", 150, 0, 180);
	// PPhi = new TH1D("PPhi", "PPhi", 500, -180, 180);
	// PEpTot = new TH1D("P_Ep_Total", "P_Ep_Total", 300, 0, 900);
	// PPhiDiff = new TH1D("Phi1-Phi2", "Phi1-Phi2", 300, 0, 360);
	// PVX = new TH1D("X_Vertex", "X_Vertex", 200, -50, 50);
	// PVY = new TH1D("Y_Vertex", "Y_Vertex", 200, -60 , 60);
	// PVZ = new TH1D("Z_Vertex", "Z_Vertex", 200, -300, 300);	
	// nTAPS = new TH1D("TAPS_Hits", "TAPS_Hits", 5, 0, 10);

	// EpdE = new TH2D("E_dE", "E_dE", 150, 0, 500, 150, 0, 8); 
	// EpEg = new TH2D("EpEg", "EpEg", 150, 0, 500, 200, 100, 900);
	// EpTp = new TH2D("EpTp", "EpTp", 150 ,0, 500, 150, 0, 180);
	// TpPp = new TH2D("Theta_Phi", "Theta_Phi", 150, 0, 180, 500, -180, 180);
	// TpdE = new TH2D("Theta_dE", "Theta_dE", 150, 0, 180, 150, 0, 8);
	// EpPVZ = new TH2D("Ep_Z_Vertex", "Ep_Z_Vertex", 150, 0, 500, 200, -300, 300);

	// dE1_dE2 = new TH2D("dE1_dE2", "dE1_dE2", 300, 0, 8, 300, 0, 8);
	// Ep1_Ep2 = new TH2D("Ep1_Ep2", "Ep1_Ep2", 150, 0, 500, 150, 0, 500);
	// PVZ1_PVZ2 = new TH2D("Z_Vertex1_Z_Vertex2","Z_Vertex1_Z_Vertex2", 400, -300, 300, 400, -300, 300);
	// PTheta1_PTheta2 = new TH2D("Theta1_Theta2", "Theta1_Theta2", 150, 0, 180, 150, 0, 180);
	// Ep1dE1 = new TH2D("E1_dE1", "E1_dE1", 150, 0, 500, 150, 0, 8);
	// Ep2dE2 = new TH2D("E2_dE2", "E2_dE2", 150, 0, 500, 150, 0, 8);
	
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

