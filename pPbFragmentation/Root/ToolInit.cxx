#define ToolInit_cxx
#include "pPbFragmentation/PbPbFragmentation.h"
#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/tools/Message.h"
#include <string>
#include <iostream>

#define EL_RETURN_CHECK( CONTEXT, EXP )			\
do {							\
if( ! EXP.isSuccess() ) {				\
Error( CONTEXT,					\
XAOD_MESSAGE( "Failed to execute: %s" ),	\
#EXP );					\
return EL::StatusCode::FAILURE;			\
}							\
} while( false )

using namespace std;
using namespace MTCorrector;

EL::StatusCode PbPbFragmentation :: setupJob (EL::Job& job)
{
	// let's initialize the algorithm to use the xAODRootAccess package
	job.useXAOD ();

	EL_RETURN_CHECK( "setupJob()", xAOD::Init() ); // call before opening first file
	cout << " Job setup done!" << endl;    
    _dR_max = (float) _jet_radius/10.;
	std::cout << "[PbPbFragmentation() : initialized with jet radius = " << _jet_radius << " dR max: " << _dR_max << std::endl;	
	
	return EL::StatusCode::SUCCESS;
}

EL::StatusCode PbPbFragmentation :: initialize ()
{
//	count number of events
	cout << " Starting initialization" << endl;
	m_eventCounter = 0;

	xAOD::TEvent* event = wk()->xaodEvent();

	// as a check, let's see the number of events in our xAOD
	Info("initialize()", "Number of events = %lli", event->getEntries() ); // print long long int

	// check if the event is MC
	const xAOD::EventInfo* eventInfo = 0;
	EL_RETURN_CHECK("execute",event->retrieve( eventInfo, "EventInfo"));
	bool isMC = false;

	if(eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) ) isMC = true;
	if (isMC != _data_switch)
	{
		cout << "MC/DATA MISMATCH. PLEASE CHECK IF THIS IS OVERLAY" << endl;
		isMC = 1; //to account for data overlay
	}

	if (isMC) { cout << "******** IS SIMULATION *********" << endl;}
	if (!isMC) {cout << "******** IS DATA *********" << endl;}
	if (!_isMB) { cout << "******** HP MODE (MB = 0) *********" << endl;}


	// Initialize and configure trigger tools
	if (_data_switch==0)
	{
		m_trigConfigTool = new TrigConf::xAODConfigTool("xAODConfigTool"); // gives us access to the meta-data
		m_trigConfigTool->msg().setLevel( MSG::ERROR );
		EL_RETURN_CHECK("initialize()",m_trigConfigTool->initialize());
		ToolHandle< TrigConf::ITrigConfigTool > trigConfigHandle( m_trigConfigTool );

		m_trigDecisionTool = new Trig::TrigDecisionTool("TrigDecisionTool");
		m_trigDecisionTool->msg().setLevel( MSG::ERROR );
		EL_RETURN_CHECK("initialize()",m_trigDecisionTool->setProperty( "ConfigTool", trigConfigHandle ));
		EL_RETURN_CHECK("initialize()",m_trigDecisionTool->setProperty( "TrigDecisionKey", "xTrigDecision"));
		EL_RETURN_CHECK("initialize()", m_trigDecisionTool->initialize() );

		cout << "Adding following " << _nTriggers << " triggers: ";
		for (int i=0;i<_nTriggers;i++){
			cout << trigger_chains.at(i) << ", ";
			_chainGroup.push_back(m_trigDecisionTool->getChainGroup(trigger_chains.at(i)));
		}
		
		TString xfn = gSystem->GetFromPipe("echo $ROOTCOREBIN");
		f_trigger_RunNumber_prescale = new TFile(xfn + "/../pPbFragmentation/data/TriggerPrescales.root","READ");
		h2_trigger_RunNumber_prescale = (TH2F*)f_trigger_RunNumber_prescale->Get("h2_Trig_RunNumber_prescale");
		
		cout << endl << "Initialize triggers finished" << endl;
	}

	//Track Selector Tool
	m_trackSelectorTool = new InDet::InDetTrackSelectionTool("InDetTrackSelectorTool");
	TrackHelperTools::SetCutLevel(m_trackSelectorTool, _cut_level.c_str());
	EL_RETURN_CHECK("initialize()",m_trackSelectorTool->initialize());
	
	//Selector for FJR
	m_trackSelectorTool_FJR = new InDet::InDetTrackSelectionTool("InDetTrackSelectorTool_FJR");
	TrackHelperTools::SetCutLevel(m_trackSelectorTool_FJR, "FJR");
	EL_RETURN_CHECK("initialize()",m_trackSelectorTool_FJR->initialize());

	//Setup the electron calibration tool
	m_EgammaCalibrationAndSmearingTool  = new CP::EgammaCalibrationAndSmearingTool("EgammaCalibrationAndSmearingTool");
	EL_RETURN_CHECK("initialize()", m_EgammaCalibrationAndSmearingTool->setProperty("ESModel", "es2015PRE") );
	EL_RETURN_CHECK("initialize()", m_EgammaCalibrationAndSmearingTool->setProperty("ResolutionType", "SigmaEff90") );
	EL_RETURN_CHECK("initialize()", m_EgammaCalibrationAndSmearingTool->setProperty("decorrelationModel", "FULL_v1") );
	EL_RETURN_CHECK("initialize()", m_EgammaCalibrationAndSmearingTool->initialize() );

	m_LHToolTight2015= new AsgElectronLikelihoodTool ("m_LHToolTight2015");

	// initialize the primary vertex container for the tool to have access
	// to the number of vertices used to adapt cuts based on the pileup
	EL_RETURN_CHECK("initialize()", m_LHToolTight2015->setProperty("primaryVertexContainer","PrimaryVertices") );

	// define the config files
	std::string confDir = "ElectronPhotonSelectorTools/offline/mc15_20150712/";
	EL_RETURN_CHECK("initialize()", m_LHToolTight2015->setProperty("ConfigFile",confDir+"ElectronLikelihoodTightOfflineConfig2015.conf" ) );
	EL_RETURN_CHECK("initialize()", m_LHToolTight2015->initialize() );

	// GRL
	TString xfn = gSystem->GetFromPipe("echo $ROOTCOREBIN");
	TString xmlfile = xfn + "/../pPbFragmentation/data/"+ _GRL;

	m_grl = new GoodRunsListSelectionTool("GoodRunsListSelectionTool");
	std::vector<std::string> vecStringGRL;
	vecStringGRL.push_back(xmlfile.Data());
	EL_RETURN_CHECK("initialize()",m_grl->setProperty( "GoodRunsListVec", vecStringGRL));
	EL_RETURN_CHECK("initialize()",m_grl->setProperty("PassThrough", false)); // if true (default) will ignore result of GRL and will just pass all events
	EL_RETURN_CHECK("initialize()",m_grl->initialize());

	//Calibration
	const std::string name = "PbPbFragmentation"; //string describing the current thread, for logging
	TString jetAlgo = "AntiKt4HI"; //String describing your jet collection, for example AntiKt4EMTopo or AntiKt4LCTopo (see below)
	TString config = "JES_MC15CHI_060316.config"; //Path to global config used to initialize the tool (see below)
	TString calibSeq = "EtaJES_DEV"; //String describing the calibration sequence to apply (see below)
									 //bool isData = false; -- not used //bool describing if the events are data or from simulation

	//insitu calibration
	TString jetAlgo_insitu = "AntiKt4EMTopo"; //String describing your jet collection, for example AntiKt4EMTopo or AntiKt4LCTopo (see below)
	TString config_insitu = "JES_2015dataset_recommendation_Feb2016.config"; //Path to global config used to initialize the tool (see below)
	const std::string name_insitu = "insitu"; //string describing the current thread, for logging
	TString calibSeq_insitu = "Insitu_DEV"; //String describing the calibration sequence to apply (see below)

	//Call the constructor. The default constructor can also be used if the arguments are set with python configuration instead
	m_jetCalibration = new JetCalibrationTool(name, jetAlgo, config, calibSeq, true);
	m_jetCalibration_insitu = new JetCalibrationTool(name_insitu, jetAlgo_insitu, config_insitu, calibSeq_insitu, true);

	//Initialize the tool
	EL_RETURN_CHECK("initialize()",m_jetCalibration->initializeTool(name));
	EL_RETURN_CHECK("initialize()",m_jetCalibration_insitu->initializeTool(name_insitu));

	//Jet Cleaning
	// initialize and configure the jet cleaning tool
	m_jetCleaning = new JetCleaningTool("JetCleaning");
	m_jetCleaning->msg().setLevel( MSG::DEBUG );
	EL_RETURN_CHECK("initialize()",m_jetCleaning->setProperty( "CutLevel", "LooseBad"));
	EL_RETURN_CHECK("initialize()",m_jetCleaning->setProperty("DoUgly", false));
	EL_RETURN_CHECK("initialize()",m_jetCleaning->initialize());

	//Uncertainty
	//JES
    //jesProv = new JetUncertaintiesTool("JESProvider");
    //EL_RETURN_CHECK("initialize()",jesProv->setProperty("JetDefinition","AntiKt4EMTopo"));
    //EL_RETURN_CHECK("initialize()",jesProv->setProperty("MCType","MC15"));
    //EL_RETURN_CHECK("initialize()",jesProv->setProperty("ConfigFile","JES_2015/ICHEP2016/JES2015_19NP.config"));
	
	//JER
	jerTool = new JERTool("JERTool");
	smearTool = new JERSmearingTool("JERSmearingTool");
	//EL_RETURN_CHECK("initialize()",jerTool->setProperty("PlotFileName", "JetResolution/Prerec2015_xCalib_2012JER_ReducedTo9NP_Plots_v2.root") );
	//EL_RETURN_CHECK("initialize()",jerTool->setProperty("CollectionName", "AntiKt4EMTopoJets") );

	// Configure the JERSmearingTool
	//smearTool->msg().setLevel(MSG::DEBUG);
	//ToolHandle<IJERTool> jerHandle(jerTool->name());
	//EL_RETURN_CHECK("initialize()",smearTool->setProperty("JERTool", jerHandle) );
	//EL_RETURN_CHECK("initialize()",smearTool->setProperty("ApplyNominalSmearing", false) );
	//EL_RETURN_CHECK("initialize()",smearTool->setProperty("isMC", true) );
	//EL_RETURN_CHECK("initialize()",smearTool->setProperty("SystematicMode", "Full") );

	// Initialize the tools
	//EL_RETURN_CHECK("initialize()",jesProv->initialize());
	
	//EL_RETURN_CHECK("initialize()",jerTool->initialize() );
	//EL_RETURN_CHECK("initialize()",smearTool->initialize() );

	//d0 function
	f_d0_cut = new TF1("f1", "[0]*exp([1]*x)+[2]*exp([3]*x)", 0.4, 500);
	f_d0_cut->SetParameters(0.472367, -0.149934, 0.193095, 0.000337765);
		
	//UEEstimator
	uee = new UEEstimator();
	uee->ptBkgrThreshold = _trkptBkgrThreshold;
	uee->jetptBkgrThreshold = _jetptBkgrThreshold;
	uee->m_maxjetdeltaR = _maxjetdeltaR;
	uee->initShapeUE(isMC);
	
	//Uncert tool
	uncertprovider = new UncertProvider(_uncert_index,_mcProbCut,_cut_level.c_str(), GetCentralityNBins(31)-1, _eff_jety);
	_mcProbCut = uncertprovider->GetMCProb();
	cout << "mc prob: " << _mcProbCut <<endl;
	//uncertprovider->JES_tool=jesProv;
	
	//Pileup tool
	
	// ZDCAnalysisTool
	m_zdcTools = new ZDC::ZdcAnalysisTool("ZdcAnalysisTool");
	// HIPileupTool
	m_hiPileup = new HI::HIPileupTool("PileupTool");
		
	if (_doPileupRejection) {
		EL_RETURN_CHECK("initialize()",m_hiPileup->initialize());
		EL_RETURN_CHECK("initialize()",m_zdcTools->initializeTool());
	}

	
	cout << " Initialization done" << endl;
	return EL::StatusCode::SUCCESS;
}

