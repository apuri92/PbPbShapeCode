#include "xAODRootAccess/Init.h"
#include "SampleHandler/SampleHandler.h"
#include "SampleHandler/ToolsDiscovery.h"
#include "EventLoop/Job.h"
#include "EventLoop/DirectDriver.h"
#include "EventLoop/CondorDriver.h"
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
#include "pPbFragmentation/pPbFragmentation.h"
#include "pPbFragmentation/PbPbFragmentation.h"
#include "pPbFragmentation/PbPbFFShape.h"
#include "pPbFragmentation/MBUEEstimator.h"
#include "pPbFragmentation/Performance.h"
#include "pPbFragmentation/TrackingPerformance.h"
#include "pPbFragmentation/BaseClass.h"
#include "pPbFragmentation/JetPerformance.h"

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
	int performance_mode;
	int jet_performance_mode;
	bool MBUE_mode;
	int run=218391;
	int doClusters=0;
	float trkptBkgrThreshold=6;
	float jetptBkgrThreshold;
	int truth_only=0;
	int DxAODMode;
	int ff_shape_mode;
	int trk_perf_mode;
	string reco_jet_collection;
	string test_reco_jet_collection;
	string truth_jet_collection;
	string grl;
	string cut_level;
	float dR_max;
	std::string input_directory;
	bool isGridJob;
	bool isCondor;
	std::string submitDir;
	std::string InDS;
	std::string OutDS;
	float dR_truth_matching;
	int centrality_scheme;
	int nFilesPerJob;
	float pTtrkCut, pTjetCut, truthpTjetCut;
	float mcProbCut; 
	std::string output_file_name="ntuple";
	bool doForward=false;
	int doSlimTree;
	std::string grid_configuration="";
	string weight_file;
	string centrality_weight;
	double pt_iso;
	double JERBalancecut;
	bool applyReweighting;
	int uncert_index;
	bool doPileupRejection;
	int useCharge;
	bool correctTrackpT;
	bool eff_jety;
	bool UseAltzDef;
	bool doFJR;
	float maxjetdeltaR;
	bool doJPRCorrection;
	float jet_y_cut;
	bool doCoarsTrackpT;
	int PythiaPowheg = 0;

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
		("cut_level,C",boost::program_options::value<std::string>(&cut_level)->default_value("ppTight"),"Trk cut level")
		("reco_jet_collection",boost::program_options::value<std::string>(&reco_jet_collection)->default_value("antikt4HIItrEM"),"Jet collection")
		("test_reco_jet_collection",boost::program_options::value<std::string>(&test_reco_jet_collection)->default_value("none"),"Test Jet collection")
		("truth_jet_collection",boost::program_options::value<std::string>(&truth_jet_collection)->default_value("antikt4Truth"),"Truth jet collection")
		("dataset",boost::program_options::value<int>(&dataset)->default_value(-1),"Type of input data")
		("isMB",boost::program_options::value<int>(&isMB)->default_value(0),"MB or HP")
		("isHerwig",boost::program_options::value<int>(&isHerwig)->default_value(0),"Pythia or Herwig")
		("data_switch",boost::program_options::value<int>(&data_switch)->default_value(0),"MC or data mode")
		("num_evt,n",boost::program_options::value<int>(&num_evt)->default_value(-1),"number of events, -1 runs all events")
		("performance_mode",boost::program_options::value<int>(&performance_mode)->default_value(0),"performance mode 0<->FF, 1<->Perf")
		("MBUE_mode",boost::program_options::value<bool>(&MBUE_mode)->default_value(0),"Run the UE estimate from MB events")
		("DxAODMode",boost::program_options::value<int>(&DxAODMode)->default_value(0),"PbPb mode 0<->Off, 1<->On (runs on DxAODs")
		("ff_shape_mode",boost::program_options::value<int>(&ff_shape_mode)->default_value(0),"FF and Jet shape mode 0<->Off, 1<->On (runs on DxAODs)")
		("trk_perf_mode",boost::program_options::value<int>(&trk_perf_mode)->default_value(0),"Track Performance mode 0<->Off, 1<->On (runs on DxAODs")
		("jet_performance_mode",boost::program_options::value<int>(&jet_performance_mode)->default_value(0),"jet performance mode 0<->FF, 1<->Perf")
		("doClusters",boost::program_options::value<int>(&doClusters)->default_value(0),"do HI Clusters study")
		("track_cuts,t",boost::program_options::value<int>(&trk_selection)->default_value(0),"apply different tracking cut")
		("eff_jety",boost::program_options::value<bool>(&eff_jety)->default_value(1),"1<-> do track efficiency in jet y bins, 0<-> do track efficiency in track eta bins")
		("run,r",boost::program_options::value<int>(&run)->default_value(0),"run number")
		("UE_trk_threshold",boost::program_options::value<float>(&trkptBkgrThreshold)->default_value(10),"UE threshold")
		("UE_jet_threshold",boost::program_options::value<float>(&jetptBkgrThreshold)->default_value(90),"UE threshold")
		("truth_only",boost::program_options::value<int>(&truth_only)->default_value(0),"Construct only truth tree")
		("isGridJob",boost::program_options::value<bool>(&isGridJob)->default_value(0),"is it grid job?")
		("isCondor",boost::program_options::value<bool>(&isCondor)->default_value(0),"is it running on condor?")
		("input_directory",boost::program_options::value<std::string>(&input_directory)->default_value("/afs/cern.ch/work/m/mrybar/xAOD/"),"name of input directory containing all files")
		("submit_directory",boost::program_options::value<std::string>(&submitDir)->default_value("submitDir"),"name of output directory")
		("InDS,i",boost::program_options::value<std::string>(&InDS)->default_value(""),"InDS for grid job")
		("OutDS,o",boost::program_options::value<std::string>(&OutDS)->default_value(""),"OutDS for grid job")
		("dR_truth_matching",boost::program_options::value<float>(&dR_truth_matching)->default_value(0.2),"dR truth matching parameter")
		("nFilesPerJob",boost::program_options::value<int>(&nFilesPerJob)->default_value(1),"Number of files per grid job")
		("jet_radius",boost::program_options::value<int>(&jet_radius)->default_value(4),"Jet radius")
		("dR_max",boost::program_options::value<float>(&dR_max)->default_value(-1.),"Jet Shape Radius")
		("centrality_scheme,s",boost::program_options::value<int>(&centrality_scheme)->default_value(1),"Centrality scheme")
		("track_pT_cut",boost::program_options::value<float>(&pTtrkCut)->default_value(0.5),"Track pT cut")
		("jet_pT_cut",boost::program_options::value<float>(&pTjetCut)->default_value(10),"Jet pT cut")
		("truth_jet_pT_cut",boost::program_options::value<float>(&truthpTjetCut)->default_value(10),"Truth jet pT cut")
		("mc_prob_cut",boost::program_options::value<float>(&mcProbCut)->default_value(0.2),"MC probability cut")
		("doSlimTree",boost::program_options::value<int>(&doSlimTree)->default_value(0),"Produce smaller tree")
		("doForward",boost::program_options::value<bool>(&doForward)->default_value(0),"Use forward jet triggers")
		("grid_configuration",boost::program_options::value<std::string>(&grid_configuration)->default_value(""),"Settings for grid configuration")
		("pt_iso",boost::program_options::value<double>(&pt_iso)->default_value(-1),"Jet pT isolation requirement, -1 <=> 100p of jet pT, 0-1: isolating pT requirement defined in percenatges of jet pT")
		("JERBalancecut",boost::program_options::value<double>(&JERBalancecut)->default_value(999.),"Trut-to-reco jet balance cut in sigmas")
		("applyReweighting",boost::program_options::value<bool>(&applyReweighting)->default_value(0),"apply reweighting to match shape between data and MC?")
		("uncertainty,u",boost::program_options::value<int>(&uncert_index)->default_value(0),"Systematic uncertainty")
		("doPileupRejection",boost::program_options::value<bool>(&doPileupRejection)->default_value(0),"Reject pileup")
		("useCharge",boost::program_options::value<int>(&useCharge)->default_value(0),"To be set to 0 (no selection), or tp -+1 to select only pozitive or negative charge tracks")
		("correctTrackpT",boost::program_options::value<bool>(&correctTrackpT)->default_value(1),"Correct the track momentum for alignemnt, now default;")
		("UseAltzDef",boost::program_options::value<bool>(&UseAltzDef)->default_value(0),"Alternative z definiton (from 7 TeV pp paper)")
		("doFJR",boost::program_options::value<bool>(&doFJR)->default_value(0),"do fake jet rejection")
		("UEmaxjetdeltaR",boost::program_options::value<float>(&maxjetdeltaR)->default_value(0.8),"maximum distance of another jet to the random cone")
		("doJPRCorrection",boost::program_options::value<bool>(&doJPRCorrection)->default_value(1),"Improve the jet position resolution by R=0.2 jets")
		("jet_y_cut",boost::program_options::value<float>(&jet_y_cut)->default_value(2.1),"max jet rapidity")
		("doCoarsTrackpT",boost::program_options::value<bool>(&doCoarsTrackpT)->default_value(0),"Use coars binning from the shape analysis")
		("PythiaPowheg,u",boost::program_options::value<int>(&PythiaPowheg)->default_value(0),"PythiaPowheg or Pythia")
		;
	 
	if (1<=dataset && dataset<=3) centrality_scheme = 1;
	if (!jet_performance_mode) doForward=false;
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
	if (performance_mode==1) data_switch = 1;

	cout << endl << endl << "*********Configuration**********" << endl;
	if (data_switch==0) cout << "Run Number:" << run <<endl;
	if (dataset==0) {cout << "Using p+Pb 5.02 TeV setup" << endl;}
	if (dataset==1) {cout << "Using p+p 2.76 TeV setup" << endl;}
	if (dataset==2 && !isHerwig) {cout << "Using p+p 7 TeV setup" << endl;}
	if (dataset==2 && isHerwig) {cout << "Using 5.02 TeV Herwig setup" << endl;}
	if (truth_only && isHerwig) {cout << "Using 5.02 TeV Herwig setup" << endl;}
	
	if (dataset==3 && PythiaPowheg) {cout << "Using p+p 5.02 TeV Pythia PowHeg" << endl;}
	if (dataset==3 && !PythiaPowheg) {cout << "Using p+p 5.02 TeV Pythia" << endl;}
	if (dataset==4) {cout << "Using Pb+Pb 5.02 TeV setup" << endl;}
	
	if (isMB==0) {cout << "Using HP data" << endl;}
	if (isMB==1) {cout << "Using MB data" << endl;} 
	if (isHerwig==1) {cout << "Using Herwig" << endl;}
	
	if (data_switch==0) cout << "Running in data mode" << endl;
	if (data_switch==1) cout << "Running in MC mode" << endl;
	
	if (performance_mode==0) cout << "Running in FF mode" << endl;
	if (performance_mode==1) cout << "Running in track performance mode" << endl;
	if (jet_performance_mode==1) cout << "Running in jet performance mode" << endl;
	
	cout << "Input directory:  " << input_directory << endl;
	cout << "Output directory: " << submitDir << endl;

	if (DxAODMode || ff_shape_mode || trk_perf_mode) cout << "Specific cut selection: " << cut_level << endl;
	else cout << "Specific cut selection: " << trk_selection << endl;

	cout << "dR truth matching parameter: " << dR_truth_matching << endl; 
	cout << "Jet collection: " << reco_jet_collection << endl;
	if (strcmp (test_reco_jet_collection.c_str(),"none") != 0) cout << "Using test jet collection: " << test_reco_jet_collection << endl;
	if (data_switch==1) cout << "Truth jet collection: " << truth_jet_collection << endl;
	cout << "grl: " << grl << endl;
	cout << "UE track threshold: " << trkptBkgrThreshold << " GeV" << endl;
	cout << "UE jet threshold: " << jetptBkgrThreshold << " GeV" << endl;
	if (truth_only==1) cout << "Running in truth mode" << endl;
	cout << "Centrality scheme: "<< centrality_scheme << endl;
	cout << "mc probability cut: "<< mcProbCut << endl;
	cout << "track pt cut: "<< pTtrkCut << endl;
	cout << "jet pt cut: "<< pTjetCut << endl;
	cout << "jet isolation pT cut: "<< pt_iso << endl;
	cout << "JER balance cut: "<< JERBalancecut << endl;
	cout << "dR shape range: " << dR_max << endl;
	if (UseAltzDef) cout << "Using alternative z definiton" << endl;
	if (doFJR) cout << "Using fake jet rejection" << endl;

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
	BaseClass* alg;

	if (DxAODMode + ff_shape_mode + trk_perf_mode + jet_performance_mode + performance_mode > 1) { cout << "*** SELECT PROPER RUNNING MODE. QUITTING. ****" << endl; return 0; }
	else if (DxAODMode) alg = new PbPbFragmentation();
	else if (ff_shape_mode) alg = new PbPbFFShape();
	else if (trk_perf_mode) alg = new TrackingPerformance();
	else if (jet_performance_mode) alg = new JetPerformance();
	else if (performance_mode) alg = new Performance();
	else if (MBUE_mode) alg = new MBUEEstimator();
	else alg = new pPbFragmentation();

	//TODO write a copy const of base class
	 
	//Set parameters
	alg->_data_switch = data_switch;
	alg->_dataset = dataset;
	alg->_GRL = grl;
	alg->_cut_level = cut_level;
	alg->_reco_jet_collection=reco_jet_collection.c_str();
	alg->_test_reco_jet_collection=test_reco_jet_collection.c_str();
	alg->_truth_jet_collection=truth_jet_collection.c_str();
	alg->_centrality_scheme = centrality_scheme;
	alg->_jet_radius = jet_radius;
	alg->_dR_max = dR_max;
	alg->_truth_only = truth_only;
	alg->_isMB = isMB;
	alg->_isHerwig = isHerwig;
	alg->_Run_Number= run;
	alg->_trkptBkgrThreshold=trkptBkgrThreshold;
	alg->_jetptBkgrThreshold=jetptBkgrThreshold;
	alg->_dR_truth_matching = dR_truth_matching;
	alg->_trk_selection=trk_selection;
	alg->_eff_jety=eff_jety;
	alg->_pTtrkCut=pTtrkCut;
	alg->_pTjetCut=pTjetCut;
	alg->_truthpTjetCut=truthpTjetCut;
	alg->_mcProbCut=mcProbCut;
	alg->_doSlimTree=doSlimTree;
	alg->_doClusters=doClusters;
	alg->_doForward=doForward;
	alg->_nTrkSelTools=4; // configure TrackSelectionTools in PbFragmentation::initialize
	alg->_pt_iso=pt_iso;
	alg->_JERBalancecut=JERBalancecut;
	alg->_applyReweighting=applyReweighting;
	alg->_uncert_index=uncert_index;
	alg->_doPileupRejection=doPileupRejection;
	alg->_useCharge=useCharge;
	alg->_correctTrackpT=correctTrackpT;
	alg->_UseAltzDef=UseAltzDef;
	alg->_doFJR=doFJR;
	alg->_maxjetdeltaR=maxjetdeltaR;
	alg->_doJPRCorrection=doJPRCorrection;
	alg->_jet_y_cut=jet_y_cut;
	alg->_doCoarsTrackpT=doCoarsTrackpT;
	alg->_PythiaPowheg=PythiaPowheg;

	//Initialzie trigger
	alg->SetTrigger_chains();
	job.algsAdd( alg );
	alg->_outputName = output_file_name.c_str(); // give the name of the output to our algorithm
	// Run the job using the local/direct driver:
	cout << "Run the job" << endl;
	//Split level protection
    job.options()->setString(EL::Job::optXaodAccessMode, EL::Job::optXaodAccessMode_athena);
	
	if(!isGridJob){
		 EL::DirectDriver driver;
		 driver.submit( job, submitDir );
	}
	else if(!isGridJob && isCondor){
		EL::CondorDriver driver;
		job.options()->setString(EL::Job::optCondorConf, "getenv = true\naccounting_group = group_atlas.uillu");
		driver.submit( job, submitDir );
	}else {
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
