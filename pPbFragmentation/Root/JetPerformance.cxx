#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include "pPbFragmentation/JetPerformance.h"
#include "pPbFragmentation/SEBCorrectorTool.h"
#include "pPbCentrality/pPbMinBiasUtil.h"
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/tools/Message.h"
#include "xAODJet/JetContainer.h"
#include "xAODMuon/MuonContainer.h"
#include "xAODTrigger/JetRoIContainer.h"
#include "xAODTrigger/JetRoIAuxContainer.h"
#include <TFile.h>
#include <TSystem.h>

using namespace std;
using namespace JetHelperTools;
using namespace MTCorrector;

ClassImp(JetPerformance)

#define EL_RETURN_CHECK( CONTEXT, EXP )			\
do {							\
  if( ! EXP.isSuccess() ) {				\
    Error( CONTEXT,					\
    XAOD_MESSAGE( "Failed to execute: %s" ),	\
    #EXP );					\
    return EL::StatusCode::FAILURE;			\
    }							\
    } while( false )


JetPerformance :: JetPerformance ()
{
}

EL::StatusCode JetPerformance :: setupJob (EL::Job& job)
{
	// let's initialize the algorithm to use the xAODRootAccess package
	job.useXAOD ();
	
	EL_RETURN_CHECK( "setupJob()", xAOD::Init() ); // call before opening first file
	cout << " Job setup done!" << endl;
	
	switch (_jet_radius){
		case 2:
		_dR_max = 0.2;
		break;
		case 3:
			_dR_max = 0.3;
		break;
		case 4:
			_dR_max = 0.4;
		break;
		case 5:
			_dR_max = 0.5;
		break;
		case 6:
			_dR_max = 0.6;
		break;  
	}
  
	std::cout << "[JetPerformance() : initialized with jet radius = " << _jet_radius << std::endl;
		  
	return EL::StatusCode::SUCCESS;
}

EL::StatusCode JetPerformance :: histInitialize ()
{   
	cout << " Setting  histograms" << endl;
	jet_tree = false;
	
	h_FCal_Et = new TH1D("h_FCal_Et",";FCal E_{T};N",100,0,5);
	h_RejectionHisto = new TH1D("RejectionHisto","RejectionHisto",9,0,9);
	SetRejectionHistogram(h_RejectionHisto);
	h_DAQErrors= new TH1D("DAQErrors","DAQErrors",4,0,4);
	
	Double_t ptTrkBins[100], etaTrkBins[100], phiTrkBins[100],finehitsBins[1000],d0z0Bins[1000],ptJetBins[100],mJetBins[100];
	Int_t ptTrkBinsN = 35, etaTrkBinsN = 50, phiTrkBinsN = 100, finehitsBinsN = 100, d0z0BinsN=600, ptJetBinsN, mJetBinsN;
	Double_t PVBins[3]={0,1,2};
	int PVBinsN=2;
   
    SetupBinning(0, "pt-jet-PbPb", ptJetBins, ptJetBinsN);
    SetupBinning(0, "m-jet", mJetBins, mJetBinsN);
	SetupBinning(0, "pt-trk", ptTrkBins, ptTrkBinsN);
	SetupBinning(0, "eta-trk", etaTrkBins, etaTrkBinsN);
	SetupBinning(0, "phi-trk", phiTrkBins, phiTrkBinsN);
	SetupBinning(0, "hits_fine", finehitsBins, finehitsBinsN);
	SetupBinning(0, "d0z0", d0z0Bins, d0z0BinsN);
   
	//Basic histograms
	hET_ETsub = new TH3D("hET_ETsub","hET_ETsub",100, 0, 500,100,-5,5,200,-50,50);
	hET_ETsub_v_dEta = new TH3D("hET_ETsub_v_dEta","hET_ETsub_v_dEta",100, 0, 500,100,-5,5,200,-50,50);
	h3_Jet_EtaPt_EtsubPerNconst = new TH3D("h3_Jet_EtaPt_EtsubPerNconst", ";p_{T,jet} [GeV]; #eta; Number of constituents", 100, 0, 500, 100, -5, 5, 200, -1, 1);
	h3_Jet_EtaPt_Nconst = new TH3D("h2_Jet_EtaPtNconst", ";p_{T,jet} [GeV]; #eta; Number of constituents", 100, 0, 500, 100, -5, 5, 100, 0, 100);
	htrig_reco_pt = new TH2D("htrig_reco_pt","htrig_reco_pt;reco jet p_{T} [GeV]; trigger jet p_{T} [GeV]",100,0,500,100,0,500);
	
	// calo jets versus clusters
    h2_JetCl_DiffEtRaw = new TH2D("h2_JetCl_DiffEtRaw", ";E_{T,jet} [GeV];E_{T,jet} - #Sigma E_{T,cluster}^{Raw} [GeV];", 100, 0, 500, 400, -100, 100); // jet pt [GeV]
    h2_JetCl_DiffEtAlt = new TH2D("h2_JetCl_DiffEtAlt", ";E_{T,jet} [GeV];E_{T,jet} - #Sigma E_{T,cluster}^{Alt} [GeV];", 100, 0, 500, 400, -100, 100); // jet pt [GeV]
    h2_JetCl_PtNconst1 = new TH2D("h2_JetCl_PtNconst1", ";p_{T,jet} [GeV]; N_{const} - N_{cl+}", 100, 0, 500, 200, -100, 100);
    h2_JetCl_PtNconst2 = new TH2D("h2_JetCl_PtNconst2", ";p_{T,jet} [GeV]; N_{const} + N_{cl-}", 100, 0, 500, 100, 0, 100);
	h3_JetCl_NegEt = new TH3D("h3_JetCl_NegEt", ";p_{T,jet} [GeV]; subtracted p_{T} [GeV]; #Sigma E_{T,cluster}^{Alt, negative} [GeV] ", 100, 0, 500, 200,-50,50, 300, -20, 10);
	
	h_triggercounter = new TH2D("h_triggercounter","h_triggercounter",_nTriggers,0,_nTriggers,2,-0.5,1.5);
	SetTrigger_hist(h_triggercounter);
    
    h3_HLT_jet_spect = new TH3D("h3_HLT_jet_spect", ";p_{T,jet} [GeV]; #eta_{jet}; #phi_{jet}", 400,0,400,100, -5., 5., 64, -TMath::Pi(), TMath::Pi()); // jet pt [GeV]
     
	TH3D* temphist_3D = nullptr;
	TH2D* temphist_2D = nullptr;
	TH1D* temphist_1D = nullptr;
	for (int i=0;i<_nTriggers;i++){
		temphist_1D = new TH1D(Form("hjetpt_%.0f",trigger_thresholds.at(i)),Form("hjetpt_%.0f",trigger_thresholds.at(i)),500,0,500);
		hjetpt_trig.push_back(temphist_1D);
	}
		
	 
	wk()->addOutput (h_RejectionHisto);
	wk()->addOutput (h_DAQErrors);
	wk()->addOutput (h_FCal_Et);
	
	wk()->addOutput (hET_ETsub);
	wk()->addOutput (hET_ETsub_v_dEta);
	wk()->addOutput (h3_Jet_EtaPt_EtsubPerNconst);
	wk()->addOutput (h3_Jet_EtaPt_Nconst);
	wk()->addOutput (htrig_reco_pt);
	wk()->addOutput (h_triggercounter);
	
	/*
	wk()->addOutput (h2_JetCl_DiffEtRaw);
    wk()->addOutput (h2_JetCl_DiffEtAlt);
    wk()->addOutput (h2_JetCl_PtNconst1);
    wk()->addOutput (h2_JetCl_PtNconst2);
    wk()->addOutput (h3_JetCl_NegEt);
    */
    
    wk()->addOutput (h3_HLT_jet_spect);
	
	/*
	for (int i=0;i<_nTriggers;i++){
		wk()->addOutput (hjetpt_trig.at(i));
	}
	*/
	for (int i=0;i<GetCentralityNBins(_centrality_scheme);i++)
	{	
		temphist_2D = new TH2D(Form("h_jet_v_mass_cent%i",i),Form("h_jet_v_mass_cent%i",i),ptJetBinsN, ptJetBins,mJetBinsN, mJetBins);
		h_jet_v_mass.push_back(temphist_2D);
		h_jet_v_mass.at(i)->Sumw2();
		wk()->addOutput (h_jet_v_mass.at(i));
		//in MC only
		if (_data_switch==1){
			temphist_2D = new TH2D(Form("h_truth_jet_v_mass_cent%i",i),Form("h_truth_jet_v_mass_cent%i",i),ptJetBinsN, ptJetBins,mJetBinsN, mJetBins);
			h_truth_jet_v_mass.push_back(temphist_2D);
			h_truth_jet_v_mass.at(i)->Sumw2();
			wk()->addOutput (h_truth_jet_v_mass.at(i));
		}
		
	}
	
	
	
	cout << " Histograms  ready, now setting tree" << endl;
	//TFile *outputFile = wk()->getOutputFile (_outputName);
	tree_performance = new TTree ("tree_performance", "Performance tree");
	//tree_performance->SetDirectory (outputFile);
	wk()->addOutput(tree_performance);

	tree_performance->Branch("event_n",&event_n,"event_n/I");
	tree_performance->Branch("run_n",&run_n,"run_n/I");
	tree_performance->Branch("lbn_n",&lbn_n,"lbn_n/I");
	tree_performance->Branch("FCalEt",&FCal_Et,"FCalEt/F");
	tree_performance->Branch("weight",&weight,"weight/F");
	tree_performance->Branch("centrality",&jet_centrality,"centrality/I");
	for (int i=0;i<_nTriggers;i++){
		tree_performance->Branch(Form("event_isTriggered_%s",trigger_chains.at(i).c_str()),&event_isTriggered[i],Form("event_isTriggered_%i/O",i));
		tree_performance->Branch(Form("trigger_prescale_%s",trigger_chains.at(i).c_str()),&trig_prescale[i],Form("trigger_prescale_%i/F",i));
	}
	
	if(jet_tree){
		

		tree_performance->Branch("jet_pt_EM",&jet_pt_EM,"jet_pt_EM/F");
		tree_performance->Branch("jet_m_EM",&jet_m_EM,"jet_m_EM/F");
		tree_performance->Branch("jet_pt_xcalib",&jet_pt_xcalib,"jet_pt_xcalib/F");
		//tree_performance->Branch("jet_pt_prexcalib",&jet_pt_prexcalib,"jet_pt_prexcalib/F");
		//tree_performance->Branch("jet_pt_seb",&jet_pt_seb,"jet_pt_seb/F");
		tree_performance->Branch("jet_pt_unsubtracted",&jet_pt_unsubtracted,"jet_pt_unsubtracted/F");
		//tree_performance->Branch("jet_pt",&jet_pt,"jet_pt/F");
		tree_performance->Branch("jet_phi",&jet_phi,"jet_phi/F");
		tree_performance->Branch("jet_eta",&jet_eta,"jet_eta/F");
		tree_performance->Branch("jet_m",&jet_m,"jet_m/F");
		tree_performance->Branch("jet_nConst",&jet_nConst,"jet_nConst/I");
		tree_performance->Branch("jet_isGood",&jet_isGood,"jet_isGood/I");
		//tree_performance->Branch("jet_NBJ_pT",&jet_NBJ_pT,"jet_NBJ_pT/F");
		//tree_performance->Branch("jet_isMuonIsolated",&jet_isMuonIsolated,"jet_isMuonIsolated/I");
	
		/*
		tree_performance->Branch("jet_jetQuality",&jet_jetQuality,"jet_jetQuality/F");
		tree_performance->Branch("jet_jetTime",&jet_jetTime,"jet_jetTime/F");
		tree_performance->Branch("jet_jethecq",&jet_jethecq,"jet_jethecq/F");
		tree_performance->Branch("jet_jetnegE",&jet_jetnegE,"jet_jetnegE/F");
		tree_performance->Branch("jet_jetemf",&jet_jetemf,"jet_jetemf/F");
		tree_performance->Branch("jet_jethecf",&jet_jethecf,"jet_jethecf/F");
		tree_performance->Branch("jet_jetfracSamplingMax",&jet_jetfracSamplingMax,"jet_jetfracSamplingMax/F");
		tree_performance->Branch("jet_jetchf",&jet_jetchf,"jet_jetchf/F");
		tree_performance->Branch("jet_jetBchCorrCell",&jet_jetBchCorrCell,"jet_jetBchCorrCell/F");
		tree_performance->Branch("jet_jetLArBadHVEnergyFrac",&jet_jetLArBadHVEnergyFrac,"jet_jetLArBadHVEnergyFrac/F");
		tree_performance->Branch("jet_jetLArBadHVNCell",&jet_jetLArBadHVNCell,"jet_jetLArBadHVNCell/I");
	    */
		if(strcmp (_test_reco_jet_collection.c_str(),"none") != 0){
			tree_performance->Branch("test_jet_pt",&test_jet_pt,"test_jet_pt/F");
			tree_performance->Branch("test_jet_pt_EM",&test_jet_pt_EM,"test_jet_pt_EM/F");
			tree_performance->Branch("test_jet_m_EM",&test_jet_m_EM,"test_jet_m_EM/F");
			tree_performance->Branch("test_jet_phi",&test_jet_phi,"test_jet_phi/F");
			tree_performance->Branch("test_jet_eta",&test_jet_eta,"test_jet_eta/F");
			tree_performance->Branch("test_jet_m",&test_jet_m,"test_jet_m/F");
			tree_performance->Branch("test_jet_isGood",&test_jet_isGood,"test_jet_isGood/I");
		}
	    /*
		tree_performance->Branch("jet_ClSub_et",&jet_ClSub_et,"jet_ClSub_et/F");
	 	tree_performance->Branch("jet_ClUnsub_et",&jet_ClUnsub_et,"jet_ClUnsub_et/F");
		tree_performance->Branch("jet_Clneg_et",&jet_Clneg_et,"jet_Clneg_et/F");
		tree_performance->Branch("jet_Clpost_et",&jet_Clpost_et,"jet_Clpost_et/F");
		tree_performance->Branch("ClUnsub_et",&ClUnsub_et);
		tree_performance->Branch("ClUnsub_eta",&ClUnsub_eta);
		tree_performance->Branch("ClUnsub_phi",&ClUnsub_phi);
		*/
		if(_data_switch==1){
			tree_performance->Branch("truth_jet_pt",&truth_jet_pt,"truth_jet_pt/F");
			tree_performance->Branch("truth_jet_phi",&truth_jet_phi,"truth_jet_phi/F");
			tree_performance->Branch("truth_jet_eta",&truth_jet_eta,"truth_jet_eta/F");
			tree_performance->Branch("truth_jet_m",&truth_jet_m,"truth_jet_m/F");
			tree_performance->Branch("truth_reco_jet_dR",&truth_reco_jet_dR,"truth_reco_jet_dR/F");
			//tree_performance->Branch("truth_jet_NBJ_pT",&truth_jet_NBJ_pT,"truth_jet_NBJ_pT/F");
		}
	
		truth_tree_performance = new TTree ("truth_tree_performance", "Performance tree");
		//truth_tree_performance->SetDirectory (outputFile);
		wk()->addOutput(truth_tree_performance);

		if(_data_switch==1){
			truth_tree_performance->Branch("event_n",&event_n,"event_n/I");
			truth_tree_performance->Branch("run_n",&run_n,"run_n/I");
			truth_tree_performance->Branch("FCalEt",&FCal_Et,"FCalEt/F");
			truth_tree_performance->Branch("truth_jet_pt",&truth_jet_pt,"truth_jet_pt/F");
			truth_tree_performance->Branch("truth_jet_phi",&truth_jet_phi,"truth_jet_phi/F");
			truth_tree_performance->Branch("truth_jet_eta",&truth_jet_eta,"truth_jet_eta/F");
			truth_tree_performance->Branch("truth_jet_m",&truth_jet_m,"truth_jet_m/F");
		}
	}
	//Event tree
	else {

		tree_performance->Branch("jet_pt_EM",&jet_pt_EM_vector);
		tree_performance->Branch("jet_pt_xcalib",&jet_pt_xcalib_vector);
		//tree_performance->Branch("jet_pt_unsubtracted",&jet_pt_unsubtracted_vector,"jet_pt_unsubtracted/F");
		tree_performance->Branch("jet_phi",&jet_phi_vector);
		tree_performance->Branch("jet_eta",&jet_eta_vector);
		tree_performance->Branch("jet_m",&jet_m_vector);
		tree_performance->Branch("jet_nConst",&jet_nConst_vector);
		tree_performance->Branch("jet_isGood",&Is_jet_Good);
		tree_performance->Branch("jet_trigger_prescale",&jet_TrigPresc_vector);
		
		tree_performance->Branch("muon_pt",&muon_pt_vector);
		tree_performance->Branch("muon_phi",&muon_phi_vector);
		tree_performance->Branch("muon_eta",&muon_eta_vector);
		tree_performance->Branch("muon_charge",&muon_charge_vector);
		tree_performance->Branch("muon_quality",&muon_quality_vector);
	
		
		if(_data_switch==1){
			tree_performance->Branch("truth_jet_pt",&truth_jet_pt_vector);
			tree_performance->Branch("truth_jet_phi",&truth_jet_phi_vector);
			tree_performance->Branch("truth_jet_eta",&truth_jet_eta_vector);
			tree_performance->Branch("truth_muon_pt",&truth_muon_pt_vector);
			tree_performance->Branch("truth_muon_phi",&truth_muon_phi_vector);
			tree_performance->Branch("truth_muon_eta",&truth_muon_eta_vector);
			tree_performance->Branch("truth_muon_charge",&truth_muon_charge_vector);
			tree_performance->Branch("truth_jet_m",&truth_jet_m_vector);
			tree_performance->Branch("truth_reco_jet_dR",&truth_reco_jet_dR_vector);
			tree_performance->Branch("truth_matched_index",&truth_matched_index);
		}	
	}	
	cout << " tree  ready" << endl; 
	
	return EL::StatusCode::SUCCESS;
	
	return EL::StatusCode::SUCCESS;
}

EL::StatusCode JetPerformance :: fileExecute ()
{
	// Here you do everything that needs to be done exactly once for every
	// single file, e.g. collect a list of all lumi-blocks processed
	return EL::StatusCode::SUCCESS;
}

EL::StatusCode JetPerformance :: changeInput (bool firstFile)
{
	return EL::StatusCode::SUCCESS;
}

EL::StatusCode JetPerformance :: initialize ()
{
	// count number of events
	
	cout << " Starting initialization" << endl;
	m_eventCounter = 0;
	
	xAOD::TEvent* event = wk()->xaodEvent();
	
	// as a check, let's see the number of events in our xAOD
	Info("initialize()", "Number of events = %lli", event->getEntries() ); // print long long int
	
	// Initialize and configure trigger tools  
	if (_data_switch==0){
		m_trigConfigTool = new TrigConf::xAODConfigTool("xAODConfigTool"); // gives us access to the meta-data
		m_trigConfigTool->msg().setLevel( MSG::ERROR ); 
		m_trigConfigTool->initialize();
		ToolHandle< TrigConf::ITrigConfigTool > trigConfigHandle( m_trigConfigTool );
		m_trigDecisionTool = new Trig::TrigDecisionTool("TrigDecisionTool");
		m_trigDecisionTool->msg().setLevel( MSG::ERROR ); 
		m_trigDecisionTool->setProperty( "ConfigTool", trigConfigHandle ); // connect the TrigDecisionTool to the ConfigTool
		m_trigDecisionTool->setProperty( "TrigDecisionKey", "xTrigDecision");
		m_trigDecisionTool->initialize();
		
		cout << "Adding following " << _nTriggers << " triggers: ";
		for (int i=0;i<_nTriggers;i++){
			cout << trigger_chains.at(i) << ", ";
			_chainGroup.push_back(m_trigDecisionTool->getChainGroup(trigger_chains.at(i)));
		}	
		cout << endl << "Initialize triggers finished" << endl;
	}
	
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
	
	//Test collection calibration
	const std::string name_test = "TestCollection"; //string describing the current thread, for logging
	TString jetAlgo_test = "AntiKt4EMTopo"; //String describing your jet collection, for example AntiKt4EMTopo or AntiKt4LCTopo (see below)
	TString config_test = "JES_MC15cRecommendation_May2016.config"; //Path to global config used to initialize the tool (see below)
	TString calibSeq_test = "JetArea_Residual_Origin_DEV"; //String describing the calibration sequence to apply (see below)
	//TString calibSeq_test = "JetArea_Residual_Origin_EtaJES_GSC_DEV"; //String describing the calibration sequence to apply (see below)
									 //bool isData = false; -- not used //bool describing if the events are data or from simulatio
	//if (_data_switch==0) calibSeq_test = "JetArea_Residual_Origin_EtaJES_GSC_Insitu_DEV";
	if (_data_switch==0) calibSeq_test = "JetArea_Residual_Origin_DEV";
	//Call the constructor. The default constructor can also be used if the arguments are set with python configuration instead
	m_jetCalibration_test = new JetCalibrationTool(name_test, jetAlgo_test, config_test, calibSeq_test, true);
	
	//Initialize the tool
	EL_RETURN_CHECK("initialize()",m_jetCalibration_test->initializeTool(name_test));
	
	//Jet Cleaning
	// initialize and configure the jet cleaning tool
	m_jetCleaning = new JetCleaningTool("JetCleaning");
	m_jetCleaning->msg().setLevel( MSG::DEBUG ); 
	EL_RETURN_CHECK("initialize()",m_jetCleaning->setProperty( "CutLevel", "LooseBad"));
	EL_RETURN_CHECK("initialize()",m_jetCleaning->setProperty("DoUgly", false));
	EL_RETURN_CHECK("initialize()",m_jetCleaning->initialize());
	
	//JetCorrector
	jetcorr = new JetCorrector();
	
	//MuonSelectionTool 
	m_muonSelection = new CP::MuonSelectionTool("MuonSelection");
	m_muonSelection->msg().setLevel( MSG::ERROR ); 
	m_muonSelection->setProperty( "MaxEta", 2.5 ); 
	m_muonSelection->setProperty( "MuQuality", 2); // 0-tight, 1-medium, 2-loose, 3-very loose
	m_muonSelection->setProperty("TrtCutOff",true);
	EL_RETURN_CHECK("initialize()", m_muonSelection->initialize() );
	 
	//Muon corection tool (MC only) 
	m_muonCalibrationAndSmearingTool = new CP::MuonCalibrationAndSmearingTool( "MuonCorrectionTool" );
	m_muonCalibrationAndSmearingTool->setProperty("Year", "Data15");
	m_muonCalibrationAndSmearingTool->setProperty("Release", "Recs2016_15_07");
	m_muonCalibrationAndSmearingTool->setProperty("SagittaCorr", false);
	m_muonCalibrationAndSmearingTool->setProperty("doSagittaMCDistortion", false); 
	m_muonCalibrationAndSmearingTool->setProperty("StatComb", true); 
	EL_RETURN_CHECK("initialize()", m_muonCalibrationAndSmearingTool->initialize() );

	cout << " Initialization done" << endl;
	return EL::StatusCode::SUCCESS;
}

//Loop over events
EL::StatusCode JetPerformance :: execute (){

	xAOD::TEvent* event = wk()->xaodEvent();
	
	// Event counter
	int statSize=1;
	if(m_eventCounter!=0)
	{
		double power=std::floor(log10(m_eventCounter));
		statSize=(int)std::pow(10.,power);
	}
	if(m_eventCounter%statSize==0) std::cout << "Event: " << m_eventCounter << std::endl;
	m_eventCounter++;
	

	//All events
	h_RejectionHisto->Fill(0.5);
		
	//---------------------------
	//     Event information
	//--------------------------- 
		
	const xAOD::EventInfo* eventInfo = 0;
	EL_RETURN_CHECK("execute",event->retrieve( eventInfo, "EventInfo"));
		
	event_n = eventInfo->eventNumber();
	run_n = eventInfo->runNumber();
	lbn_n = eventInfo->lumiBlock();
	
	if (eventInfo->actualInteractionsPerCrossing() > 0) cout << "Pileup" << eventInfo->actualInteractionsPerCrossing() << endl;
		
	//if(m_eventCounter%statSize==0) cout << "EventNumber " << event_n << endl;
		
	FCal_Et = 0;
	int cent_bin=0;
	//TODO FCAl weights if needed
	if (_centrality_scheme>1) {
		//Centrality
		const xAOD::HIEventShapeContainer* calos = 0;
		EL_RETURN_CHECK("execute()",event->retrieve( calos, "CaloSums" ));
		FCal_Et = 0;
		int x = 0;
		xAOD::HIEventShapeContainer::const_iterator calo_itr = calos->begin();
		xAOD::HIEventShapeContainer::const_iterator calo_end = calos->end();
		for( ; calo_itr != calo_end; ++calo_itr ) {
			if (x == 5) {
				FCal_Et = ((*calo_itr)->et() * 0.001 * 0.001 );
			}
			x++;
		}
		cent_bin = GetCentralityBin(_centrality_scheme, FCal_Et);
		h_FCal_Et->Fill(FCal_Et);
	}  
	if (cent_bin < 0) {
		h_RejectionHisto->Fill(1.5);
		return EL::StatusCode::SUCCESS;
	}
	
		
	// check if the event is data or MC
	bool isMC = false;
	bool isHIJING = false; //For centrality
	// check if the event is MC
	if(eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) )
	{
		isMC = true;
		isHIJING = true;
		_data_switch=1;
	}
	else
	{
		const xAOD::TruthParticleContainer * particles = 0;
		if( event->xAOD::TVirtualEvent::retrieve(particles, "TruthParticles", true) )
		{
			// this is overlay
			isMC = true;
			isHIJING = false;
			_data_switch=1;
		}
		else
		{
			_data_switch=0;
		}
	} 
	  
	// GRL
	if(!isMC){ 
		if(!m_grl->passRunLB(*eventInfo)){
			h_RejectionHisto->Fill(2.5);
			return EL::StatusCode::SUCCESS; // go to next event
		}
	}
	
	//Vertex requirement
	const xAOD::VertexContainer * vertices = 0;
	if ( !event->retrieve( vertices, "PrimaryVertices" ).isSuccess() ){
		Error("execute()", "Failed to retrieve VertexContainer container. Exiting." );
		return EL::StatusCode::FAILURE;
	}
		
	if(vertices->size()<2) {
	 	h_RejectionHisto->Fill(3.5);
	 	return EL::StatusCode::SUCCESS;
	}

	xAOD::VertexContainer::const_iterator vtx_itr = vertices->begin();
	xAOD::VertexContainer::const_iterator vtx_end = vertices->end();
	// find primary vertex
	const xAOD::Vertex* primaryVertex = 0;
	for(;vtx_itr!=vtx_end;++vtx_itr)
	{
		if((*vtx_itr)->vertexType()==xAOD::VxType::PriVtx) {
			primaryVertex = (*vtx_itr);
			break;
		}
	}
	

	
	//DAQ errors
	if(!isMC){
		if(   (eventInfo->errorState(xAOD::EventInfo::LAr)==xAOD::EventInfo::Error ) || (eventInfo->errorState(xAOD::EventInfo::Tile)==xAOD::EventInfo::Error ) || (eventInfo->errorState(xAOD::EventInfo::SCT)==xAOD::EventInfo::Error ) || (eventInfo->isEventFlagBitSet(xAOD::EventInfo::Core, 18) ) ){
			h_RejectionHisto->Fill(4.5);
			return EL::StatusCode::SUCCESS; // go to next event
		}
	}
 	 
	h_RejectionHisto->Fill(5.5);
	
	// trigger
	if (_data_switch==0)
	{
		int event_passed_trigger=0;

		for (int i=0;i<_nTriggers;i++){

			event_isTriggered[i] = false;
			event_isTriggered[i] =  _chainGroup.at(i)->isPassed();
			trig_prescale[i] =  _chainGroup.at(i)->getPrescale();
			h_triggercounter->Fill(i, (Double_t) event_isTriggered[i]);
			if(event_isTriggered[i]) event_passed_trigger=1;
		}

		if(!event_passed_trigger) return EL::StatusCode::SUCCESS; // go to next event
		else h_RejectionHisto->Fill(8.5);
	}
	
	//Tracks
	const xAOD::TrackParticleContainer* recoTracks = 0;
	if ( !event->retrieve( recoTracks, "InDetTrackParticles" ).isSuccess() ){
		Error("execute()", "Failed to retrieve Reconstructed Track container. Exiting." );
		return EL::StatusCode::FAILURE;
	}
	
	
	float event_weight = 1;
	double max_pt = 1;
	bool keep_event = false;
	TLorentzVector jet4vector;
	
	//----------------------------
    // Container with clusters
    //--------------------------- 
    const xAOD::CaloClusterContainer* cls = 0;
    if (_doClusters) EL_RETURN_CHECK("execute()",event->retrieve( cls, "HIClusters" ));
	
	//Jet vectors
	vector<float> jet_pt_SEB_vector,jet_pt_prexcalib_vector,test_jet_phi_vector,test_jet_eta_vector,test_jet_m_vector,test_jet_pt_EM_vector,test_jet_m_EM_vector,test_jet_pt_vector;
	vector<int> Is_test_jet_Good, Is_dummyJet,jet_IsTrig_vector;
	vector<float> truth_jet_y_vector;
	vector<int> truth_jet_indices, truth_jet_isDummy_vector;
	vector<int> hasTruth,jet_jetLArBadHVNCell_vector;
	vector<int> isTriggered[_nTriggers];
	vector<float> jet_NBJ_pT_vector,truth_jet_NBJ_pT_vector;
	vector<float> trig_EF_jet_pt, trig_EF_jet_phi, trig_EF_jet_eta;
	vector<float> antikt2_pt,antikt2_phi,antikt2_eta;
	vector<float> jet_jetQuality_vector,jet_jetTime_vector,jet_jethecq_vector,jet_jetnegE_vector,jet_jetemf_vector,jet_jethecf_vector,jet_jetchf_vector,jet_jetfracSamplingMax_vector,jet_jetLArBadHVEnergyFrac_vector,jet_jetBchCorrCell_vector;
	
	// Clear vectors
	hasTruth.clear();
	
	truth_jet_indices.clear();
	truth_jet_pt_vector.clear();
	truth_jet_eta_vector.clear();
	truth_jet_phi_vector.clear();
	truth_jet_m_vector.clear();
	truth_jet_y_vector.clear();
	
	jet_pt_EM_vector.clear();
	jet_m_EM_vector.clear();
	jet_pt_unsubtracted_vector.clear();
	jet_pt_SEB_vector.clear();
	jet_pt_prexcalib_vector.clear();
	jet_pt_xcalib_vector.clear();
	jet_phi_vector.clear();
	jet_eta_vector.clear();
	jet_m_vector.clear();
	jet_nConst_vector.clear();
	
	jet_jetQuality_vector.clear();
	jet_jetTime_vector.clear();
	jet_jethecq_vector.clear();
	jet_jetnegE_vector.clear();
	jet_jetemf_vector.clear();
	jet_jethecf_vector.clear();
	jet_jetfracSamplingMax_vector.clear();
	jet_jetchf_vector.clear();
	jet_jetBchCorrCell_vector.clear();
	jet_jetLArBadHVEnergyFrac_vector.clear();
	jet_jetLArBadHVNCell_vector.clear();
	Is_jet_Good.clear();
	
	jet_IsTrig_vector.clear();
	jet_TrigPresc_vector.clear();
	
	muon_eta_vector.clear();
	muon_phi_vector.clear();
	muon_pt_vector.clear();
	muon_charge_vector.clear();
	muon_quality_vector.clear();
	
	truth_muon_eta_vector.clear();
	truth_muon_phi_vector.clear();
	truth_muon_pt_vector.clear();
	truth_muon_charge_vector.clear();
		
	for (int j=0;j<_nTriggers;j++){
		isTriggered[j].clear();
	}
		
	//***** TRUTH JETS *****
	 
	// ---- GETTING TRUTH JETS ----
	std::vector<float> JetTruthPt, JetTruthPhi, JetTruthEta;
	const xAOD::JetContainer * jet_truth = 0;
	if(_data_switch == 1)
	{
		EL_RETURN_CHECK("execute()",event->retrieve( jet_truth, _truth_jet_collection.c_str() ));

		truth_jet_pt_vector.clear(); truth_jet_phi_vector.clear(); truth_jet_eta_vector.clear();;

		xAOD::JetContainer::const_iterator jet_itr = jet_truth->begin();
		xAOD::JetContainer::const_iterator jet_end = jet_truth->end();
		for( ; jet_itr != jet_end; ++jet_itr )
		{
			xAOD::JetFourMom_t jet_truth_4mom = (*jet_itr)->jetP4();

			truth_jet_pt    = (jet_truth_4mom.pt() * 0.001 );
			truth_jet_eta    = (jet_truth_4mom.eta());
			truth_jet_phi    = (jet_truth_4mom.phi());
			truth_jet_m   = (jet_truth_4mom.M()*0.001);

			if (truth_jet_pt>max_pt) 	//event weight from leading truth jet
			{
				event_weight = jetcorr->GetJetWeight(truth_jet_pt, truth_jet_eta, truth_jet_phi);
				max_pt = truth_jet_pt;
			}
		
			//if (pt < _truthpTjetCut) continue;
			//if (fabs(eta)>2.1) continue;

			//filling truth pt/eta/phi vectors
			truth_jet_pt_vector.push_back(truth_jet_pt);
			truth_jet_m_vector.push_back(truth_jet_m);
			truth_jet_phi_vector.push_back(truth_jet_phi);
			truth_jet_eta_vector.push_back(truth_jet_eta);
			////Get rapidity
			jet4vector.SetPtEtaPhiM(truth_jet_pt, truth_jet_eta, truth_jet_phi, truth_jet_m);
			truth_jet_y_vector.push_back(jet4vector.Rapidity());
			
			if (fabs(truth_jet_eta)>2.1) continue;
			h_truth_jet_v_mass.at(cent_bin)->Fill(truth_jet_pt,truth_jet_m,event_weight);
			if (truth_jet_pt < _pTjetCut) continue;
			truth_tree_performance->Fill();
		}
	}
	weight=event_weight;
	
	//***** Reco jets*****
		

	xAOD::TStore store; //For calibration
	const xAOD::JetContainer* jets = 0;
	EL_RETURN_CHECK("execute()",event->retrieve( jets, _reco_jet_collection.c_str() ));
	
	xAOD::JetContainer::const_iterator jet_itr = jets->begin();
	xAOD::JetContainer::const_iterator jet_end = jets->end();
	
	xAOD::JetContainer* updatedjets = new xAOD::JetContainer();
	xAOD::AuxContainerBase* updatedjetsAux = new xAOD::AuxContainerBase();
	updatedjets->setStore( updatedjetsAux );
		
	int jet_counter=0;
	for( ; jet_itr != jet_end; ++jet_itr ) {
		
		xAOD::Jet newjet;// = new xAOD::Jet();
		newjet.makePrivateStore( **jet_itr );

		const xAOD::JetFourMom_t jet_4mom_def = newjet.jetP4();
		float def_jet_pt  = (jet_4mom_def.pt() * 0.001);
		
		//cout << " Def:  " << def_jet_pt << endl;

		xAOD::JetFourMom_t jet_4mom = newjet.jetP4("JetSubtractedScaleMomentum"); //getting SubtractedScale instead of EMScale because EMScale is not in DFAntiKt4HI
		float uncalib_jet_pt  = (jet_4mom.pt() * 0.001);
		float uncalib_jet_m  = (jet_4mom.M() * 0.001);

		const xAOD::JetFourMom_t jet_4mom_unsubtracted = newjet.jetP4("JetUnsubtractedScaleMomentum");
		float unsubtracted_jet_pt  = (jet_4mom_unsubtracted.pt() * 0.001);

		hET_ETsub->Fill(def_jet_pt,jet_4mom_def.eta(),unsubtracted_jet_pt-uncalib_jet_pt);

		//newjet.setJetP4("JetPileupScaleMomentum",jet_4mom); //Setting PileupScale and ConstitScale because they are not in DFAntiKt4HI
		newjet.setJetP4("JetConstitScaleMomentum",jet_4mom_unsubtracted); //Required
		//EL_RETURN_CHECK("execute()", m_jetCalibration->applyCalibration( newjet ) );
		
		//cout << " Calib " << (newjet.pt() * 0.001);
		
		//For non DF collections
		if (_dataset == 3 && _reco_jet_collection.find("DF") == std::string::npos) // pp: calibrate with sequence set in ShapeToolInit
		{
			newjet.setJetP4("JetEMScaleMomentum",jet_4mom); //Required
			EL_RETURN_CHECK("execute()", m_jetCalibration->applyCalibration( newjet ) ); //calibrates with sequence EtaJes_Insitu for data, EtaJes for MC
		}
		
		const xAOD::JetFourMom_t jet_4mom_xcalib = newjet.jetP4();
		newjet.setJetP4("JetGSCScaleMomentum", jet_4mom_xcalib);

		//Cross-calibration
		if (_data_switch==0) EL_RETURN_CHECK("execute()", m_jetCalibration_insitu->applyCalibration( newjet ) );
		
		//Jet quality moment
		if(_dataset==3 && !m_jetCleaning->accept( **jet_itr )) Is_jet_Good.push_back(0);
		else Is_jet_Good.push_back(1);
		
		jet_nConst = (*jet_itr)->numConstituents();
		jet_pt  = (newjet.pt() * 0.001);
		jet_eta = newjet.eta();
		jet_phi = newjet.phi();
		jet_m = newjet.m()*0.001;
		
		jet_pt_EM_vector.push_back(uncalib_jet_pt);
		jet_m_EM_vector.push_back(uncalib_jet_m);
		jet_pt_unsubtracted_vector.push_back(unsubtracted_jet_pt);
		jet_pt_xcalib_vector.push_back(jet_pt);
		
		jet_phi_vector.push_back(jet_phi);
		jet_eta_vector.push_back(jet_eta);
		jet_m_vector.push_back(jet_m);
		jet_nConst_vector.push_back(jet_nConst);
		
		//jet quality

        //@MS 20160511: failing on compile => commented out, please fix
		float jetQuality      = 0; //(*jet_itr)->getAttribute<float>(xAOD::JetAttribute::LArQuality);
        float jetTime         = (*jet_itr)->getAttribute<float>(xAOD::JetAttribute::Timing);
        float hecq            = (*jet_itr)->getAttribute<float>(xAOD::JetAttribute::HECQuality);
        float negE            = (*jet_itr)->getAttribute<float>(xAOD::JetAttribute::NegativeE);
        float emf             = (*jet_itr)->getAttribute<float>(xAOD::JetAttribute::EMFrac);
        float hecf            = (*jet_itr)->getAttribute<float>(xAOD::JetAttribute::HECFrac); 
		float fracSamplingMax = (*jet_itr)->getAttribute<float>(xAOD::JetAttribute::FracSamplingMax);
		float BchCorrCell             = (*jet_itr)->getAttribute<float>(xAOD::JetAttribute::BchCorrCell);

        //@MS 20160511: failing on compile => commented out, please fix
        float LArBadHVEnergyFrac       = 0;  //(*jet_itr)->getAttribute<float>(xAOD::JetAttribute::LArBadHVEnergyFrac); 
		float LArBadHVNCell = 0; //(*jet_itr)->getAttribute<int>(xAOD::JetAttribute::LArBadHVNCell);
		std::vector<float> SumPtTrkPt1000;
		(*jet_itr)->getAttribute(xAOD::JetAttribute::SumPtTrkPt1000,SumPtTrkPt1000);
		float chf             = SumPtTrkPt1000.size() > 0 ? SumPtTrkPt1000.at(0)/(*jet_itr)->pt() : -1;
		
		jet_jetQuality_vector.push_back(jetQuality);
		jet_jetTime_vector.push_back(jetTime);
		jet_jethecq_vector.push_back(hecq);
		jet_jetnegE_vector.push_back(negE);
		jet_jetemf_vector.push_back(emf);
		jet_jethecf_vector.push_back(hecf);
		jet_jetfracSamplingMax_vector.push_back(fracSamplingMax);
		jet_jetchf_vector.push_back(chf);
		jet_jetBchCorrCell_vector.push_back(BchCorrCell);
		jet_jetLArBadHVEnergyFrac_vector.push_back(LArBadHVEnergyFrac);
		jet_jetLArBadHVNCell_vector.push_back(LArBadHVNCell);
		
		if (_data_switch==0) //is data, so need prescales and trigger decisions (at jet level)
		{
			bool is_trig = false;
			double presc = -1;

			for (int k=0;k<_nTriggers;k++)
			{
				if(event_isTriggered[k] && (jet_pt > jet_pt_trig[k][0] && jet_pt <= jet_pt_trig[k][1]))
				{
					is_trig = true;
					presc =  trig_prescale[k];
					break;
				}
			}
			jet_TrigPresc_vector.push_back(presc);
			jet_IsTrig_vector.push_back(is_trig);		
		}
		
		//Get mass
		/*
		const xAOD::JetConstituentVector constituents = (*jet_itr)->getConstituents();
		int iterator=0;
		for (xAOD::JetConstituentVector::iterator itr = constituents.begin(); itr != constituents.end(); ++itr){ 
			iterator++; 
			const xAOD::CaloCluster* cl=static_cast<const xAOD::CaloCluster*>(itr->rawConstituent());
			float sumET= cl->rawE() * 0.001 / std::cosh( cl->rawEta() ); 
		}
		*/
		
	}
	
	
	//Test jet loop
	if (strcmp (_test_reco_jet_collection.c_str(),"none") != 0){	
		const xAOD::JetContainer* test_jets = 0;
		EL_RETURN_CHECK("execute()",event->retrieve( test_jets, _test_reco_jet_collection.c_str() ));
		
		xAOD::JetContainer::const_iterator test_jet_itr = test_jets->begin();
		xAOD::JetContainer::const_iterator test_jet_end = test_jets->end();
		
		xAOD::JetContainer* updatedtestjets = new xAOD::JetContainer();
		xAOD::AuxContainerBase* updatedtestjetsAux = new xAOD::AuxContainerBase();
		updatedtestjets->setStore( updatedtestjetsAux );

		
		for( ; test_jet_itr != test_jet_end; ++test_jet_itr ) {
		
			xAOD::Jet newjet_test;
			newjet_test.makePrivateStore( **test_jet_itr );
			const xAOD::JetFourMom_t jet_4mom = ( *test_jet_itr )->jetP4("JetPileupScaleMomentum");
			float uncalib_jet_pt  = (jet_4mom.pt() * 0.001);
			float uncalib_jet_m  = (jet_4mom.M() * 0.001);
			EL_RETURN_CHECK("execute()", m_jetCalibration_test->applyCalibration( newjet_test ) );
						
			jet_pt = newjet_test.pt()*0.001;
			jet_eta = newjet_test.eta();
			jet_phi = newjet_test.phi();
			jet_m   = newjet_test.m()*0.001;
				
			//Jet quality moment
			if( !m_jetCleaning->accept( **test_jet_itr )) Is_test_jet_Good.push_back(0);
			else Is_test_jet_Good.push_back(1);
				
			test_jet_pt_EM_vector.push_back(uncalib_jet_pt);
			test_jet_m_EM_vector.push_back(uncalib_jet_m);
			test_jet_pt_vector.push_back(jet_pt);		
			test_jet_phi_vector.push_back(jet_phi);
			test_jet_eta_vector.push_back(jet_eta);
			test_jet_m_vector.push_back(jet_m);
		}
	}
	
	//Truh matching, MC only	  
	if(_data_switch == 1){ 
		for(unsigned int i=0; i<jet_pt_xcalib_vector.size(); i++){
			jet_pt = jet_pt_xcalib_vector.at(i) / 1.0;
			//if(jet_pt < 20.0)continue;
			jet_eta = jet_eta_vector.at(i);
			jet_phi = jet_phi_vector.at(i);
			bool hasTruthJet=false;
			int TM_index=-1;
			float dR_min=999.;
			for(unsigned int j=0; j<truth_jet_pt_vector.size(); j++){
				truth_jet_pt = truth_jet_pt_vector.at(j);
				truth_jet_m = truth_jet_m_vector.at(j);
				truth_jet_eta = truth_jet_eta_vector.at(j);
				truth_jet_phi = truth_jet_phi_vector.at(j);
				float R = DeltaR(truth_jet_phi,truth_jet_eta,jet_phi,jet_eta);
				if(R < 0.2 && R < dR_min){
					hasTruthJet = true;
					dR_min=R;
					TM_index = j;
				}
			}
			if (hasTruthJet) hasTruth.push_back(1);   	
			else hasTruth.push_back(0);
			truth_jet_indices.push_back(TM_index);
		}
	}
	
	//Muons
	const xAOD::MuonContainer * muons = 0;
	EL_RETURN_CHECK("execute()",event->retrieve(muons,"Muons"));
	
	// create a shallow copy of the muons container 
	std::pair< xAOD::MuonContainer*, xAOD::ShallowAuxContainer* > muons_shallowCopy = xAOD::shallowCopyContainer( *muons ); 
	    
	xAOD::MuonContainer::iterator muonSC_itr = (muons_shallowCopy.first)->begin(); 
	xAOD::MuonContainer::iterator muonSC_end = (muons_shallowCopy.first)->end(); 
	for( ; muonSC_itr != muonSC_end; ++muonSC_itr ) { 
		if(m_muonCalibrationAndSmearingTool->applyCorrection(**muonSC_itr) == CP::CorrectionCode::Error){ // apply correction and check return code 
		Error("execute()", "MuonCalibrationAndSmearingTool returns Error CorrectionCode"); 
		}
		muon_pt_vector.push_back((*muonSC_itr)->pt() * 0.001);
		muon_eta_vector.push_back((*muonSC_itr)->eta());
		muon_phi_vector.push_back((*muonSC_itr)->phi());
		muon_charge_vector.push_back((*muonSC_itr)->charge());
		muon_quality_vector.push_back(m_muonSelection->getQuality(**muonSC_itr));
		if ((*muonSC_itr)->pt() * 0.001 > 5.) keep_event = true;
	}
	//Truth muons
	if(_data_switch == 1){
	
		const xAOD::TruthParticleContainer * particles = 0;
		if ( !event->retrieve( particles, "TruthParticles" ).isSuccess() ){
			Error("execute()", "Failed to retrieve TruthParticle container. Exiting." );
			return EL::StatusCode::FAILURE;
		}
		xAOD::TruthParticleContainer::const_iterator truth_itr = particles->begin();
		xAOD::TruthParticleContainer::const_iterator truth_end = particles->end();
		for( ; truth_itr!=truth_end; ++truth_itr){

			int ty=TrackHelperTools::getTypeTruth((*truth_itr)->barcode(),(*truth_itr)->pdgId(), (*truth_itr)->status(), (*truth_itr)->charge());				
			if(ty!=1 && ty!=5) continue;				
			//get the tracks....
			float pt = (*truth_itr)->pt()/ 1000.0;				
			float eta = (*truth_itr)->eta();				
			float phi = (*truth_itr)->phi();
			float charge = (*truth_itr)->charge();

			if (fabs((*truth_itr)->pdgId())==13) {truth_muon_pt_vector.push_back(pt);truth_muon_eta_vector.push_back(eta);truth_muon_phi_vector.push_back(phi);truth_muon_charge_vector.push_back(charge);}
		} 
	} 


	
	//Here is the analysis code
	for(unsigned int i=0; i<jet_pt_xcalib_vector.size(); i++){
			
		float jet_weight =1.;		
		jet_weight *=event_weight;
		if (_data_switch==0)
		{
			if (!jet_IsTrig_vector.at(i)) continue;
			jet_weight = jet_TrigPresc_vector.at(i);
		}
		
		jet_pt = jet_pt_xcalib_vector.at(i);
		jet_eta = jet_eta_vector.at(i); 
		jet_phi = jet_phi_vector.at(i);
		jet_m = jet_m_vector.at(i);
		jet_Et = sqrt(pow(jet_pt,2)+pow(jet_m,2));
		jet_pt_EM = jet_pt_EM_vector.at(i);
		jet_m_EM = jet_m_EM_vector.at(i);
		jet_pt_unsubtracted = jet_pt_unsubtracted_vector.at(i);
		//jet_pt_seb = jet_pt_SEB_vector.at(i);
		//jet_pt_prexcalib = jet_pt_prexcalib_vector.at(i);
		jet_pt_xcalib = jet_pt_xcalib_vector.at(i);
		jet_nConst = jet_nConst_vector.at(i);
		jet_isGood = Is_jet_Good.at(i);
		jet_isMuonIsolated = true;
		jet_centrality=cent_bin;
		
		jet_jetQuality = jet_jetQuality_vector.at(i);
		jet_jetTime = jet_jetTime_vector.at(i);
		jet_jethecq = jet_jethecq_vector.at(i);
		jet_jetnegE = jet_jetnegE_vector.at(i);
		jet_jetemf = jet_jetemf_vector.at(i);
		jet_jethecf = jet_jethecf_vector.at(i);
		jet_jetfracSamplingMax = jet_jetfracSamplingMax_vector.at(i);
		jet_jetchf = jet_jetchf_vector.at(i);
		jet_jetBchCorrCell = jet_jetBchCorrCell_vector.at(i);
		jet_jetLArBadHVEnergyFrac = jet_jetLArBadHVEnergyFrac_vector.at(i);
		jet_jetLArBadHVNCell = jet_jetLArBadHVNCell_vector.at(i);
		//weight=jet_weight;		 
		
		//10 GeV reco cut
		if (jet_pt_xcalib<10.) continue;
		
		if(_data_switch == 1){
			//cout << hasTruth.size() << endl;
			jet_hasTruth = hasTruth.at(i);
			
			if (jet_hasTruth){
				truth_jet_eta = truth_jet_eta_vector.at(truth_jet_indices[i]);
				truth_jet_phi = truth_jet_phi_vector.at(truth_jet_indices[i]);
				truth_jet_pt = truth_jet_pt_vector.at(truth_jet_indices[i]);
				truth_jet_m = truth_jet_m_vector.at(truth_jet_indices[i]);
				truth_reco_jet_dR = DeltaR(truth_jet_phi,truth_jet_eta,jet_phi,jet_eta);
				truth_reco_jet_dR_vector.push_back(truth_reco_jet_dR);
				truth_matched_index.push_back(truth_jet_indices[i]);
			}
			else{
				truth_jet_eta = -100.;
				truth_jet_phi = -100.;
				truth_jet_pt = 0;
				truth_reco_jet_dR_vector.push_back(-1);
				truth_matched_index.push_back(-1);
			}
			
		}
		
		//test jets
		if(strcmp (_test_reco_jet_collection.c_str(),"none") != 0){
			float dR_min = 999.;
			int index=-1;
			for(unsigned int j=0; j<test_jet_pt_EM_vector.size(); j++){			
				float R = DeltaR(test_jet_phi_vector.at(j),test_jet_eta_vector.at(j),jet_phi,jet_eta);
				if (R<0.4 && R<dR_min){
					index=j;
				}
			}
			if (index>=0){
				test_jet_eta = test_jet_eta_vector.at(index);		 
				test_jet_phi = test_jet_phi_vector.at(index);
				test_jet_m = test_jet_m_vector.at(index);
				test_jet_pt_EM = test_jet_pt_EM_vector.at(index);
				test_jet_m_EM = test_jet_m_EM_vector.at(index);
				test_jet_pt = test_jet_pt_vector.at(index);
				test_jet_isGood = Is_test_jet_Good.at(index);
			}	
		}
		//Only good jets to histograms
		if (jet_isGood) {
						
			hET_ETsub->Fill(jet_pt_xcalib,jet_eta,jet_pt_unsubtracted-jet_pt_EM);
			h3_Jet_EtaPt_EtsubPerNconst->Fill(jet_pt_xcalib,jet_eta,(jet_pt_unsubtracted-jet_pt_EM)/jet_nConst);
			h3_Jet_EtaPt_Nconst->Fill(jet_pt_xcalib,jet_eta,jet_nConst);
		
			if (jet_pt>50.){
				for(unsigned int j=0; j<jet_pt_xcalib_vector.size(); j++){
					if (i==j) continue;
					hET_ETsub_v_dEta->Fill(jet_pt_xcalib_vector.at(j),jet_eta-jet_eta_vector.at(j),jet_pt_unsubtracted-jet_pt_EM);
				}	
			}
			if (fabs(jet_eta)<=2.1) h_jet_v_mass.at(cent_bin)->Fill(jet_pt,jet_m,jet_weight);
		}
		if (_doClusters) {
		   
		   ClUnsub_et.clear();
		   ClUnsub_eta.clear();
		   ClUnsub_phi.clear();
		
		   xAOD::CaloClusterContainer::const_iterator cl_itr = cls->begin();
		   xAOD::CaloClusterContainer::const_iterator cl_end = cls->end();

		   Int_t n_pos=0, n_neg=0;
		   Float_t etRaw = 0, etAlt = 0, etCl = 0 , et_pos = 0, et_neg = 0;
		   for( ; cl_itr != cl_end; ++cl_itr ) {
		      Float_t deltaR = DeltaR(jet_phi, jet_eta, (*cl_itr)->phi(), (*cl_itr)->eta() );
		      if (deltaR > 0.4) continue;
		      etRaw += (*cl_itr)->rawE() / cosh( (*cl_itr)->rawEta() ) * 0.001;
		      etAlt += (*cl_itr)->altE() / cosh( (*cl_itr)->altEta() ) * 0.001;

		      if ((*cl_itr)->altE() > 0) {n_pos++;et_pos+=(*cl_itr)->altE() / cosh( (*cl_itr)->altEta() ) * 0.001;}
		      if ((*cl_itr)->altE() < 0) {n_neg++;et_neg+=(*cl_itr)->altE() / cosh( (*cl_itr)->altEta() ) * 0.001;}
		      ClUnsub_et.push_back((*cl_itr)->rawE() / cosh( (*cl_itr)->rawEta() ) * 0.001);
		      ClUnsub_eta.push_back((*cl_itr)->rawEta() );
		      ClUnsub_phi.push_back((*cl_itr)->rawPhi() );
		   }     
		   
		   if (jet_isGood){
			   h2_JetCl_DiffEtRaw->Fill( jet_Et, etRaw - jet_Et );
			   h2_JetCl_DiffEtAlt->Fill( jet_Et, etAlt - jet_Et );
			   
			   h3_JetCl_NegEt->Fill( jet_Et, jet_pt_unsubtracted-jet_pt_EM, et_neg );

			   h2_JetCl_PtNconst1->Fill( jet_pt, jet_nConst - n_pos );
			   h2_JetCl_PtNconst2->Fill( jet_pt, jet_nConst + n_neg );
		   }
		   
		   jet_ClSub_et = etRaw;
		   jet_ClUnsub_et = etAlt;
		   jet_Clneg_et = et_neg;
		   jet_Clpost_et = et_pos;
		   
		
       } // doClusters		
      if (jet_pt < _pTjetCut) continue;
      keep_event = true;
      if (jet_tree) tree_performance->Fill();
		
	}
	if (!jet_tree && keep_event) tree_performance->Fill();
	//delete shallow copy of muon container 
	delete muons_shallowCopy.first; 
	delete muons_shallowCopy.second; 
	return EL::StatusCode::SUCCESS;
}

EL::StatusCode JetPerformance :: postExecute (){
	return EL::StatusCode::SUCCESS;
}

EL::StatusCode JetPerformance :: finalize (){
	//xAOD::TEvent* event = wk()->xaodEvent();
	
	//cleaning muons
	if(m_muonCalibrationAndSmearingTool){ 
		delete m_muonCalibrationAndSmearingTool; 
		m_muonCalibrationAndSmearingTool = 0; 
	 } 
	 if(m_muonSelection){ 
		delete m_muonSelection; 
		m_muonSelection = 0; 
	 }
	
	// cleaning up trigger tools
	if (_data_switch==0){	  
		if( m_trigConfigTool ) {
			delete m_trigConfigTool;
			m_trigConfigTool = 0;
		}
		
		if( m_trigDecisionTool ) {
			delete m_trigDecisionTool;
			m_trigDecisionTool = 0;
		}
	}
	
	// cleaning GRL
	if (m_grl) {
		delete m_grl;
		m_grl = 0;
	}
	
	//cleaning cleaning :)
	if( m_jetCleaning ) {
		delete m_jetCleaning;
		m_jetCleaning = 0;
	}
	
	return EL::StatusCode::SUCCESS;
}
	
EL::StatusCode JetPerformance :: histFinalize (){  
	cout<<"Events = "<< m_eventCounter<<endl;
	return EL::StatusCode::SUCCESS;
}
