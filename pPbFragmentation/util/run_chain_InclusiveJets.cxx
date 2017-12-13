#include "xAODRootAccess/Init.h"
#include "SampleHandler/SampleHandler.h"
#include "SampleHandler/ToolsDiscovery.h"
#include "EventLoop/Job.h"
#include "EventLoop/DirectDriver.h"
#include "EventLoopGrid/PrunDriver.h"
#include "EventLoopGrid/GridDriver.h"
#include "SampleHandler/DiskListLocal.h"
#include <EventLoopAlgs/NTupleSvc.h>
#include <EventLoop/OutputStream.h>
#include <TSystem.h>

#include "SampleHandler/ScanDir.h"

#include <string>
#include <vector>
#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <exception>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <boost/shared_ptr.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include "pPbFragmentation/InclusiveJetsEventLoop.h"
#include "pPbFragmentation/BaseClass.h"

int main(int argc, char *argv[])
{
	
	//configuration
	int data_switch = 0; //
	int dataset = 0; // 0 - pPb, 1 - 2p76TeV pp, 2 - 7TeV pp , 3 - 5p02TeV pp, 4 - 5p02TeV PbPb
	int isMB=0;
	int isHerwig=0;
	int num_evt;
	int trk_selection=1;
	int jet_radius = 4;
	int run=218391;
	int doClusters=0;
	float trkptBkgrThreshold=4;
	int truth_only=0;
	string reco_jet_collection;
	string test_reco_jet_collection;
	string truth_jet_collection;
	string grl;
	std::string input_directory;
	bool isGridJob;
	std::string submitDir;
	std::string InDS;
	std::string OutDS;
	float dR_truth_matching;
	int centrality_scheme;
	int nFilesPerJob;
	std::string output_file_name="Ntuple";
	bool doForward=false;
	std::string grid_configuration="";
	
	//Boost configuration
	//1) command line only: options can only be given on command line, not in config file
	boost::program_options::options_description cmd_only_options("command line only options");
	std::string config_file_name;

	cmd_only_options.add_options() //note unusual syntax when adding options!
		("help,h","produce help message")
		("config,c",boost::program_options::value<std::string>(&config_file_name),"name of configuration file");
	
	//2) main options: most likely to be set by user, can be specified both via command line or config
	//explaination is included in help message
	boost::program_options::options_description main_options("main options");

	main_options.add_options()
		("output_file",boost::program_options::value<std::string>(&output_file_name)->default_value("ntuple"),"name of output root file")
		("grl,g",boost::program_options::value<std::string>(&grl)->default_value("GRL_pPb_5p02TeV_2013.xml"),"grl file name")
		("reco_jet_collection",boost::program_options::value<std::string>(&reco_jet_collection)->default_value("antikt4HIItrEM"),"Jet collection")
		("test_reco_jet_collection",boost::program_options::value<std::string>(&test_reco_jet_collection)->default_value("none"),"Test Jet collection")
		("truth_jet_collection",boost::program_options::value<std::string>(&truth_jet_collection)->default_value("antikt4Truth"),"Truth jet collection")
		("dataset",boost::program_options::value<int>(&dataset)->default_value(-1),"Type of input data")
		("isMB",boost::program_options::value<int>(&isMB)->default_value(0),"MB or HP")
		("isHerwig",boost::program_options::value<int>(&isHerwig)->default_value(0),"Pythia or Herwig")
		("data_switch",boost::program_options::value<int>(&data_switch)->default_value(0),"MC or data mode")
		("num_evt,n",boost::program_options::value<int>(&num_evt)->default_value(-1),"number of events, -1 runs all events")
		("doClusters",boost::program_options::value<int>(&doClusters)->default_value(0),"do HI Clusters study")
		("track_cuts,t",boost::program_options::value<int>(&trk_selection)->default_value(0),"apply different tracking cut")
		("run,r",boost::program_options::value<int>(&run)->default_value(0),"run number")
		("UE_threshold",boost::program_options::value<float>(&trkptBkgrThreshold)->default_value(4),"UE threshold")
		("truth_only",boost::program_options::value<int>(&truth_only)->default_value(0),"Construct only truth tree")
		("isGridJob",boost::program_options::value<bool>(&isGridJob)->default_value(0),"is it grid job?")
		("input_directory",boost::program_options::value<std::string>(&input_directory)->default_value(""),"name of input directory containing all files")
		("submit_directory",boost::program_options::value<std::string>(&submitDir)->default_value("submitDir"),"name of output directory")
		("InDS,i",boost::program_options::value<std::string>(&InDS)->default_value(""),"InDS for grid job")
		("OutDS,o",boost::program_options::value<std::string>(&OutDS)->default_value(""),"OutDS for grid job")
		("dR_truth_matching",boost::program_options::value<float>(&dR_truth_matching)->default_value(0.3),"dR truth matching parameter")
		("nFilesPerJob",boost::program_options::value<int>(&nFilesPerJob)->default_value(1),"Number of files per grid job")
		("jet_radius",boost::program_options::value<int>(&jet_radius)->default_value(4),"Jet radius")
		("centrality_scheme,s",boost::program_options::value<int>(&centrality_scheme)->default_value(1),"Centrality scheme")
		("grid_configuration",boost::program_options::value<std::string>(&grid_configuration)->default_value(""),"Settings for grid configuration")
		;
	 
	if (1<=dataset && dataset<=3) centrality_scheme = 1;
	//combine options types for parsing
	//all options may be specified on command line
	boost::program_options::options_description cmdline_options; 
	cmdline_options.add(cmd_only_options).add(main_options);
	
	//all options except command line only may be specified in config file
	boost::program_options::options_description config_options; 
	config_options.add(main_options);
	
	boost::program_options::variables_map vm;

	//first parse command line
	try
	{
		boost::program_options::store(boost::program_options::parse_command_line(argc, argv, cmdline_options), vm);
		boost::program_options::notify(vm);
	}
	catch(std::exception& e)
	{
		std::cerr << "Bad command line argument" << std::endl;
		std::cerr << e.what() << std::endl;
		return 1;
	}

	//if config was specified, also parse config file
	if(vm.count("config"))
	{
		ifstream config_stream(config_file_name.c_str());
		try
		{
			boost::program_options::store(boost::program_options::parse_config_file(config_stream,cmdline_options), vm);
			boost::program_options::notify(vm);
		}
		catch(std::exception& e)
		{
			std::cerr << "Bad config file argument" << std::endl;
			std::cerr << e.what() << std::endl;
			return 1;
		}
	}
	
	//Need to be in MC mode for performance

	cout << endl << endl << "*********Configuration**********" << endl;
	if (data_switch==0) cout << "Run Number:" << run <<endl;
	if (dataset==0) {cout << "Using p+Pb 5.02 TeV setup" << endl;}
	if (dataset==1) {cout << "Using p+p 2.76 TeV setup" << endl;}
	if (dataset==2 && !isHerwig) {cout << "Using p+p 7 TeV setup" << endl;}
	if (dataset==2 && isHerwig) {cout << "Using 5.02 TeV Herwig setup" << endl;}
	if (truth_only && isHerwig) {cout << "Using 5.02 TeV Herwig setup" << endl;}
	
	if (dataset==3) {cout << "Using p+p 5.02 TeV setup" << endl;}
	if (dataset==4) {cout << "Using Pb+Pb 5.02 TeV setup" << endl;}
	
	if (isMB==0) {cout << "Using HP data" << endl;}
	if (isMB==1) {cout << "Using MB data" << endl;} 
	if (isHerwig==1) {cout << "Using Herwig" << endl;}
	
	if (data_switch==0) cout << "Running in data mode" << endl;
	if (data_switch==1) cout << "Running in MC mode" << endl;
	
	cout << "Input directory:  " << input_directory << endl;
	cout << "Output directory: " << submitDir << endl;
	cout << "Output ntuple: " << output_file_name << endl;

	cout << "Specific cut selection: " << trk_selection << endl;
	cout << "dR truth matching parameter: " << dR_truth_matching << endl; 
	cout << "Jet collection: " << reco_jet_collection << endl;
	if (strcmp (test_reco_jet_collection.c_str(),"none") != 0) cout << "Using test jet collection: " << test_reco_jet_collection << endl;
	if (data_switch==1) cout << "Truth jet collection: " << truth_jet_collection << endl;
	cout << "grl: " << grl << endl;
	cout << "UE threshold: " << trkptBkgrThreshold << " GeV" << endl;
	if (truth_only==1) cout << "Running in truth mode" << endl;
	cout << "Centrality scheme: "<< centrality_scheme << endl;
	
	if (strcmp (grid_configuration.c_str(),"") != 0)  cout << "Additional grid configuration: " << grid_configuration.c_str() << endl;
	
	cout << "********************************" << endl << endl << endl;
	
	// Set up the job for xAOD access:
	xAOD::Init().ignore();

	// Construct the samples to run on:
	SH::SampleHandler sh;

	// Get input file (! be careful about path -- not to include last folder !)
	if (!isGridJob){
		SH::ScanDir().filePattern("*").scan(sh, input_directory);
	}
	else {
		SH::scanDQ2 (sh, InDS.c_str());
		sh.setMetaString( "nc_grid_filter", "*AOD*");
	}
	sh.setMetaString( "nc_tree", "CollectionTree" );

	// Print what we found:
	sh.print();

	// Create an EventLoop job:
	cout << "Creating EventLoop job" << endl;
	EL::Job job;

	//Set outputFile
	EL::OutputStream output(output_file_name.c_str());
	job.outputAdd (output);
	EL::NTupleSvc *ntuple = new EL::NTupleSvc(output_file_name.c_str());
	job.algsAdd (ntuple);

	job.sampleHandler( sh );
	 
	cout << "Seting maximum events to " << num_evt << endl;
	job.options()->setDouble (EL::Job::optMaxEvents, num_evt);
			
	// To automatically delete submitDir
	job.options()->setDouble(EL::Job::optRemoveSubmitDir, 1);


	// Add our analysis to the job:
	cout << "Add our analysis to the job" << endl;
	BaseClass* alg = new InclusiveJetsEventLoop();

//     InclusiveJetsEventLoop *alg = new InclusiveJetsEventLoop();
        
        
//         printf("HH 1 \n");
	 

	 
        
	//Set parameters
	alg->_data_switch = data_switch;
	alg->_dataset = dataset;
	alg->_GRL = grl;
	alg->_reco_jet_collection=reco_jet_collection.c_str();
	alg->_test_reco_jet_collection=test_reco_jet_collection.c_str();
	alg->_truth_jet_collection=truth_jet_collection.c_str();
	alg->_centrality_scheme = centrality_scheme;
	alg->_jet_radius = jet_radius;
	alg->_truth_only = truth_only;
	alg->_isMB = isMB;
	alg->_isHerwig = isHerwig;
	alg->_Run_Number= run;
	alg->_trkptBkgrThreshold=trkptBkgrThreshold; 
	alg->_dR_truth_matching = dR_truth_matching;
	alg->_trk_selection=trk_selection;
	alg->_doForward=doForward;
	
	//Initialzie trigger

    
    //<RS>
//     alg->SetTrigger_chains();
    //</RS>

	job.algsAdd( alg );
//         printf("HH 4 \n");
	alg->_outputName = output_file_name.c_str(); // give the name of the output to our algorithm
	// Run the job using the local/direct driver:
	cout << "Run the job" << endl;
	//Split level protection
	//@MS 20160511 -- does not compile, please fix:
//         printf("HH 5 \n");
	    job.options()->setString(EL::Job::optXaodAccessMode, EL::Job::optXaodAccessMode_athena);
//         printf("HH 6 \n");
	
	if(!isGridJob){
		 EL::DirectDriver driver;
// 	         printf("HH 7 \n");
         driver.submit( job, submitDir );
	}
	else {
		EL::PrunDriver driver;
		driver.options()->setString("nc_outputSampleName",OutDS.c_str());
		if (nFilesPerJob > 0) driver.options()->setDouble("nc_nFilesPerJob", nFilesPerJob);
		driver.options()->setString("nc_cmtConfig", "x86_64-slc6-gcc49-opt");
		driver.options()->setDouble(EL::Job::optGridMergeOutput, 1); //run merging jobs for all samples before downloading (recommended)
		if (strcmp (grid_configuration.c_str(),"") != 0) job.options()->setString (EL::Job::optSubmitFlags, grid_configuration.c_str()); //allow task duplication
		driver.submitOnly( job, submitDir );
		 
	}
	cout << "We are done!" << endl;
	
	std::cout << "done looping" << std::endl;
}
