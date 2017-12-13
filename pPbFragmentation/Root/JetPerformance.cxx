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
	
	h_FCal_Et = new TH1D("h_FCal_Et",";FCal E_{T};N",100,0,5);
	h_RejectionHisto = new TH1D("RejectionHisto","RejectionHisto",7,0,7);
	SetRejectionHistogram(h_RejectionHisto);
	h_DAQErrors= new TH1D("DAQErrors","DAQErrors",4,0,4);
	
	Double_t ptTrkBins[1000], etaTrkBins[1000], phiTrkBins[1000],finehitsBins[1000],d0z0Bins[1000];
	Int_t ptTrkBinsN = 35, etaTrkBinsN = 50, phiTrkBinsN = 100, finehitsBinsN = 100, d0z0BinsN=600;
	Double_t PVBins[3]={0,1,2};
	int PVBinsN=2;
   
   
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
	
	wk()->addOutput (h2_JetCl_DiffEtRaw);
    wk()->addOutput (h2_JetCl_DiffEtAlt);
    wk()->addOutput (h2_JetCl_PtNconst1);
    wk()->addOutput (h2_JetCl_PtNconst2);
    wk()->addOutput (h3_JetCl_NegEt);
    
    wk()->addOutput (h3_HLT_jet_spect);
	
	for (int i=0;i<_nTriggers;i++){
		wk()->addOutput (hjetpt_trig.at(i));
	}
	
	cout << " Histograms  ready, now setting tree" << endl;
	TFile *outputFile = wk()->getOutputFile (_outputName);
	tree_performance = new TTree ("tree_performance", "Performance tree"); // TODO - track charge
	tree_performance->SetDirectory (outputFile);

	tree_performance->Branch("event_n",&event_n,"event_n/I");
	tree_performance->Branch("run_n",&run_n,"run_n/I");
	tree_performance->Branch("lbn_n",&lbn_n,"lbn_n/I");
	tree_performance->Branch("FCalEt",&FCalEt,"FCalEt/F");
	
	for (int i=0;i<_nTriggers;i++){
		tree_performance->Branch(Form("event_isTriggered_%i",i),&event_isTriggered[i],Form("event_isTriggered_%i/O",i));
		tree_performance->Branch(Form("trigger_prescale_%i",i),&trig_prescale[i],Form("trigger_prescale_%i/F",i));
		tree_performance->Branch(Form("jet_isTriggered_%i",i),&jet_isTriggered[i],Form("jet_isTriggered_%i/I",i));
	}

	tree_performance->Branch("jet_pt_EM",&jet_pt_EM,"jet_pt_EM/F");
	tree_performance->Branch("jet_pt_xcalib",&jet_pt_xcalib,"jet_pt_xcalib/F");
	tree_performance->Branch("jet_pt_prexcalib",&jet_pt_prexcalib,"jet_pt_prexcalib/F");
	tree_performance->Branch("jet_pt_seb",&jet_pt_seb,"jet_pt_seb/F");
	tree_performance->Branch("jet_pt_unsubtracted",&jet_pt_unsubtracted,"jet_pt_unsubtracted/F");
	tree_performance->Branch("jet_pt",&jet_pt,"jet_pt/F");
	tree_performance->Branch("jet_phi",&jet_phi,"jet_phi/F");
	tree_performance->Branch("jet_eta",&jet_eta,"jet_eta/F");
	tree_performance->Branch("jet_m",&jet_m,"jet_m/F");
	tree_performance->Branch("jet_nConst",&jet_nConst,"jet_nConst/I");
	tree_performance->Branch("jet_isGood",&jet_isGood,"jet_isGood/I");
	tree_performance->Branch("jet_NBJ_pT",&jet_NBJ_pT,"jet_NBJ_pT/F");
	tree_performance->Branch("jet_isMuonIsolated",&jet_isMuonIsolated,"jet_isMuonIsolated/I");
	tree_performance->Branch("jet_centrality",&jet_centrality,"jet_centrality/I");
	
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
	
	if(strcmp (_test_reco_jet_collection.c_str(),"none") != 0){
		tree_performance->Branch("test_jet_pt_EM",&test_jet_pt_EM,"test_jet_pt_EM/F");
		tree_performance->Branch("test_jet_phi",&test_jet_phi,"test_jet_phi/F");
		tree_performance->Branch("test_jet_eta",&test_jet_eta,"test_jet_eta/F");
		tree_performance->Branch("test_jet_m",&test_jet_m,"test_jet_m/F");
		tree_performance->Branch("test_jet_isGood",&test_jet_isGood,"test_jet_isGood/I");
	}

	tree_performance->Branch("jet_ClSub_et",&jet_ClSub_et,"jet_ClSub_et/F");
 	tree_performance->Branch("jet_ClUnsub_et",&jet_ClUnsub_et,"jet_ClUnsub_et/F");
	tree_performance->Branch("jet_Clneg_et",&jet_Clneg_et,"jet_Clneg_et/F");
	tree_performance->Branch("jet_Clpost_et",&jet_Clpost_et,"jet_Clpost_et/F");
	tree_performance->Branch("ClUnsub_et",&ClUnsub_et);
	tree_performance->Branch("ClUnsub_eta",&ClUnsub_eta);
	tree_performance->Branch("ClUnsub_phi",&ClUnsub_phi);

	if(_data_switch==1){
		tree_performance->Branch("truth_jet_pt",&truth_jet_pt,"truth_jet_pt/F");
		tree_performance->Branch("truth_jet_phi",&truth_jet_phi,"truth_jet_phi/F");
		tree_performance->Branch("truth_jet_eta",&truth_jet_eta,"truth_jet_eta/F");
		tree_performance->Branch("truth_jet_m",&truth_jet_m,"truth_jet_m/F");
		tree_performance->Branch("truth_reco_jet_dR",&truth_reco_jet_dR,"truth_reco_jet_dR/F");
		tree_performance->Branch("truth_jet_NBJ_pT",&truth_jet_NBJ_pT,"truth_jet_NBJ_pT/F");
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
    const std::string name = "pPbFragmentation"; //string describing the current thread, for logging
    TString jetAlgo = "AntiKt4HI"; //String describing your jet collection, for example AntiKt4EMTopo or AntiKt4LCTopo (see below)
    TString config = "JES_MC15CHI_042316.config"; //Path to global config used to initialize the tool (see below)
    TString calibSeq = "EtaJES_DEV"; //String describing the calibration sequence to apply (see below)
    bool isData = false; //bool describing if the events are data or from simulation
	
	//Call the constructor. The default constructor can also be used if the arguments are set with python configuration instead
	m_jetCalibration = new JetCalibrationTool(name, jetAlgo, config, calibSeq, isData);
	
	//Initialize the tool
	EL_RETURN_CHECK("initialize()",m_jetCalibration->initializeTool(name));
	
	//Jet Cleaning
	// initialize and configure the jet cleaning tool
	m_jetCleaning = new JetCleaningTool("JetCleaning");
	m_jetCleaning->msg().setLevel( MSG::DEBUG ); 
	EL_RETURN_CHECK("initialize()",m_jetCleaning->setProperty( "CutLevel", "LooseBad"));
	EL_RETURN_CHECK("initialize()",m_jetCleaning->setProperty("DoUgly", false));
	EL_RETURN_CHECK("initialize()",m_jetCleaning->initialize());

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
	
	if(m_eventCounter%statSize==0) cout << "EventNumber " << event_n << endl;
		
	double FCal_Et = 0;
	int cent_bin=0;
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
	// check if the event is MC
	if(eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) ){
		isMC = true;
		_data_switch=1; 
	}
	else{
		_data_switch=0;
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
	
	//One evnt only
	//if (1585081075!=event_n) return EL::StatusCode::SUCCESS;
	
	//TODO vertex position?
	
	//TODO rapidity gap?
	
	//DAQ errors 
        /*
        @MS 20160511: failing on compile => commented out, please fix
	if(!isMC){
		if(   (eventInfo->errorState(xAOD::EventInfo::LAr)==xAOD::EventInfo::Error ) || (eventInfo->errorState(xAOD::EventInfo::Tile)==xAOD::EventInfo::Error ) || (eventInfo->errorState(xAOD::EventInfo::SCT)==xAOD::EventInfo::Error ) || (eventInfo->isEventFlagBitSet(xAOD::EventInfo::Core, 18) ) ){
			h_RejectionHisto->Fill(4.5);
			if(   eventInfo->errorState(xAOD::EventInfo::LAr)==xAOD::EventInfo::Error ) h_DAQErrors->Fill(0.5); if(   eventInfo->errorState(xAOD::EventInfo::Tile)==xAOD::EventInfo::Error ) h_DAQErrors->Fill(1.5); if(   eventInfo->errorState(xAOD::EventInfo::SCT)==xAOD::EventInfo::Error ) h_DAQErrors->Fill(2.5); if(   eventInfo->isEventFlagBitSet(xAOD::EventInfo::Core, 18) ) h_DAQErrors->Fill(3.5);
			//return EL::StatusCode::SUCCESS; // go to the next event
		}
	}
        */
 	 
	h_RejectionHisto->Fill(5.5);
	
	// trigger
	
	
	if (_truth_only==0) {
		
		for (int i=0;i<_nTriggers;i++){
			trig_prescale[i] = 1.0;
			event_isTriggered[i] = false;
			jet_isTriggered[i] = 0;
		}

		if (_data_switch==0){
			
			int event_passed_trigger=0;
			
			for (int i=0;i<_nTriggers;i++){
				event_isTriggered[i] =  _chainGroup.at(i)->isPassed();
				trig_prescale[i] =  _chainGroup.at(i)->getPrescale();
				h_triggercounter->Fill(i, (Double_t) event_isTriggered[i]);
				if(event_isTriggered[i]) event_passed_trigger=1;
			}
			
			if(!event_passed_trigger) return EL::StatusCode::SUCCESS; // go to next event
			else h_RejectionHisto->Fill(6.5);
		}
	}  
	
	
	//Tracks
	const xAOD::TrackParticleContainer* recoTracks = 0;
	if ( !event->retrieve( recoTracks, "InDetTrackParticles" ).isSuccess() ){
		Error("execute()", "Failed to retrieve Reconstructed Track container. Exiting." );
		return EL::StatusCode::FAILURE;
	}
	
	
	//----------------------------
    // Container with clusters
    //--------------------------- 
    const xAOD::CaloClusterContainer* cls = 0;
    if (_doClusters) EL_RETURN_CHECK("execute()",event->retrieve( cls, "HIClusters" ));
	
	//Jet vectors
	vector<float> jet_pt_EM_vector, test_jet_pt_EM_vector, jet_pt_unsubtracted_vector, jet_pt_SEB_vector,jet_pt_prexcalib_vector,jet_pt_xcalib_vector,jet_phi_vector,jet_eta_vector,jet_m_vector,test_jet_phi_vector,test_jet_eta_vector,test_jet_m_vector;
	vector<int> Is_jet_Good, Is_test_jet_Good, Is_dummyJet;
	vector<float> truth_jet_eta_vector,truth_jet_m_vector,truth_jet_phi_vector,truth_jet_pt_vector;
	vector<int> truth_jet_indices, truth_jet_isDummy_vector,jet_nConst_vector;
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
	
	jet_pt_EM_vector.clear();
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
		
	for (int j=0;j<_nTriggers;j++){
		isTriggered[j].clear();
	}
	
	
		
	//***** HLT Jets *****
	
	
	if(_data_switch==0 && !_isMB) {   
		const xAOD::JetContainer * hlt_jet = 0;
		string trigger_container="HLT_xAOD__JetContainer_" + _trigger_collection;
		  
		if( ! event->retrieve( hlt_jet, trigger_container.c_str()).isSuccess() ) {
			Error("execute()", Form("failed to retrieve %s", trigger_container.c_str()));
		}
		
		const xAOD::Jet *hlt_jet_obj = 0;
		//cout << "n HLT jets: " << hlt_jet->size() << endl;
		for(unsigned int HLTjet_itr = 0; HLTjet_itr < hlt_jet->size(); ++HLTjet_itr) {
			hlt_jet_obj = hlt_jet->at(HLTjet_itr);
			trig_EF_jet_pt.push_back(  hlt_jet_obj->pt() );
			trig_EF_jet_phi.push_back( hlt_jet_obj->phi() );
			trig_EF_jet_eta.push_back( hlt_jet_obj->eta() );
			h3_HLT_jet_spect->Fill(hlt_jet_obj->pt() * 0.001  , hlt_jet_obj->eta(), hlt_jet_obj->phi());   
		}
	}
	
	
	//***** TRUTH JETS *****
	 
	std::vector<float> JetTruthPt, JetTruthPhi, JetTruthEta;
	const xAOD::JetContainer * jet_truth = 0;
	if(isMC) { 
		if( ! event->retrieve( jet_truth, _truth_jet_collection.c_str()).isSuccess() ) {
			Error("execute()", Form("failed to retrieve %s", _truth_jet_collection.c_str()));
		}
		
		truth_jet_pt_vector.clear(); truth_jet_phi_vector.clear(); truth_jet_eta_vector.clear();truth_jet_m_vector.clear();
		
		const xAOD::Jet *truth_jet_obj = 0;
		for(unsigned int jet_itr = 0; jet_itr < jet_truth->size(); ++jet_itr) {
			truth_jet_obj = jet_truth->at(jet_itr);
		
			truth_jet_pt_vector.push_back(  truth_jet_obj->pt() * 0.001 );
			truth_jet_phi_vector.push_back( truth_jet_obj->phi() );
			truth_jet_eta_vector.push_back( truth_jet_obj->eta() );
			truth_jet_m_vector.push_back( truth_jet_obj->m() );  
			truth_jet_isDummy_vector.push_back(0);
		}
		truth_jet_NBJ_pT_vector = MTCorrector::GetIsolation(truth_jet_pt_vector,truth_jet_eta_vector,truth_jet_phi_vector,_jet_radius);
	}
	
	
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
		jet_counter++;
		xAOD::Jet* newjet = new xAOD::Jet();
		newjet->makePrivateStore( **jet_itr );
		updatedjets->push_back( newjet );	
		
		const xAOD::JetFourMom_t jet_4mom_def = newjet->jetP4();
		const xAOD::JetFourMom_t jet_4mom = newjet->jetP4("JetEMScaleMomentum");
		float def_jet_pt  = (jet_4mom_def.pt() * 0.001);
		float unsubtracted_jet_pt  = (jet_4mom.pt() * 0.001);  
		float uncalib_jet_pt  = (jet_4mom.pt() * 0.001);
		
		if (_reco_jet_collection.find("HI") != std::string::npos) {
			const xAOD::JetFourMom_t jet_4mom_unsubtracted = newjet->jetP4("JetUnsubtractedScaleMomentum");
			unsubtracted_jet_pt  = (jet_4mom_unsubtracted.pt() * 0.001);   
		}
		
		newjet->setJetP4("JetPileupScaleMomentum", jet_4mom);	
		
		EL_RETURN_CHECK("execute()", m_jetCalibration->applyCalibration( *newjet ) );
		
		if (_reco_jet_collection.find("HI") != std::string::npos) {
			jet_pt  = (newjet->pt() * 0.001);
			jet_eta = newjet->eta();
			jet_phi = newjet->phi();
			jet_m   = newjet->m()*0.001;
		}
		else{
			jet_pt  = (jet_4mom_def.pt() * 0.001);
			jet_eta = jet_4mom_def.eta();
			jet_phi = jet_4mom_def.phi();
			jet_m   = jet_4mom_def.M()*0.001;
		}
		
		//Jet quality moment
		if( !m_jetCleaning->accept( **jet_itr )) Is_jet_Good.push_back(0);
		else Is_jet_Good.push_back(1);
		
		jet_nConst = (*jet_itr)->numConstituents();
		
		jet_pt_EM_vector.push_back(uncalib_jet_pt);
		jet_pt_unsubtracted_vector.push_back(unsubtracted_jet_pt);
		jet_pt_SEB_vector.push_back(jet_pt);
		jet_pt_prexcalib_vector.push_back(jet_4mom_def.pt() * 0.001); //TODO right now is default
		jet_pt_xcalib_vector.push_back(jet_pt);
		
		jet_phi_vector.push_back(jet_phi);
		jet_eta_vector.push_back(jet_eta);
		jet_m_vector.push_back(jet_m);
		jet_nConst_vector.push_back(jet_nConst);
		
		//cout << "EM " << uncalib_jet_pt << " def " << jet_4mom_def.pt() * 0.001 << " calib " << jet_pt << " eta " << jet_eta << " phi " << jet_phi << endl; 
		
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

	}
	
	//Test jet loop
	if (strcmp (_test_reco_jet_collection.c_str(),"none") != 0){	
		const xAOD::JetContainer* test_jets = 0;
		EL_RETURN_CHECK("execute()",event->retrieve( test_jets, _test_reco_jet_collection.c_str() ));
	
		xAOD::JetContainer::const_iterator test_jet_itr = test_jets->begin();
		xAOD::JetContainer::const_iterator test_jet_end = test_jets->end();
		for( ; test_jet_itr != test_jet_end; ++test_jet_itr ) {
		
			//TODO calibration
			const xAOD::JetFourMom_t jet_4mom = ( *test_jet_itr )->jetP4("JetPileupScaleMomentum");
			float uncalib_jet_pt  = (jet_4mom.pt() * 0.001);
						
			jet_eta = jet_4mom.eta();
			jet_phi = jet_4mom.phi();
			jet_m   = jet_4mom.M()*0.001;
				
			//Jet quality moment
			if( !m_jetCleaning->accept( **test_jet_itr )) Is_test_jet_Good.push_back(0);
			else Is_test_jet_Good.push_back(1);
				
			test_jet_pt_EM_vector.push_back(uncalib_jet_pt);		
			test_jet_phi_vector.push_back(jet_phi);
			test_jet_eta_vector.push_back(jet_eta);
			test_jet_m_vector.push_back(jet_m);
		}
	}
	//TODO isolation after calibration
	//test_jet_NBJ_pT_vector = MTCorrector::GetIsolation(jet_pt_xcalib_vector,jet_eta_vector,jet_phi_vector,_jet_radius);
	
	//Trigger handeling, Data only
	
	if (_data_switch==0){	
		
		for(unsigned int i=0; i<jet_pt_xcalib_vector.size(); i++){
			
			jet_pt = jet_pt_xcalib_vector.at(i);			 
			jet_eta = jet_eta_vector.at(i);
			jet_phi = jet_phi_vector.at(i);
			
			for (int j=0;j<_nTriggers;j++){
				trigger[j] = false;
			}
			//Only for HP
			if (_isMB==0){
				//matching to trigger jets
				
				for(unsigned int j=0; j<trig_EF_jet_pt.size(); j++){	 
					float trig_jet_eta = trig_EF_jet_eta.at(j);
					float trig_jet_phi = trig_EF_jet_phi.at(j);
					float trig_jet_pt = trig_EF_jet_pt.at(j) / 1000.0;
					
					for (int k=0;k<_nTriggers;k++){			 
						if(event_isTriggered[k] && trig_jet_pt > trigger_thresholds.at(k)){
							if(jet_pt > jet_pt_trig[k][0] && jet_pt < jet_pt_trig[k][1]) {				 
								float R = DeltaR(trig_jet_phi,trig_jet_eta,jet_phi,jet_eta);
								if(R < 0.4) {trigger[k] = true; htrig_reco_pt->Fill(jet_pt,trig_jet_pt);}	
							}					
						}
					}	
				}
				
				//No matching tro trigger jets
				for (int k=0;k<_nTriggers;k++){			 
					if(event_isTriggered[k]){
						if(jet_pt > jet_pt_trig[k][0] && jet_pt < jet_pt_trig[k][1]) //comment to match Yakov analysis 
							{				 
							trigger[k] = true; 	
						}					
					}
				}
			}
			for (int j=0;j<_nTriggers;j++){			 	
				if (trigger[j] && _isMB==0){
					isTriggered[j].push_back(1);
					hjetpt_trig[j]->Fill(jet_pt,_chainGroup[j]->getPrescale());		 	
				}
				else
				{
					if(_isMB==1 && event_isTriggered[j]) {
						isTriggered[j].push_back(1);
						if(!Is_dummyJet.at(i)) hjetpt_trig[j]->Fill(jet_pt,_chainGroup[j]->getPrescale());
					}
					else isTriggered[j].push_back(0);
				}
			}	
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
	
	//Here is the analysis code
	
	for(unsigned int i=0; i<jet_pt_xcalib_vector.size(); i++){
			
		int triggered_at_least_once=0;
		if (_data_switch==0){
			for (int j=0;j<_nTriggers;j++){
				jet_isTriggered[j] = isTriggered[j].at(i);
				if(jet_isTriggered[j]){
					triggered_at_least_once=1;
				}
			}			
			//skip ther rest
			if(!triggered_at_least_once) continue;
		}
		
		jet_pt = jet_pt_xcalib_vector.at(i);
		jet_eta = jet_eta_vector.at(i); 
		jet_phi = jet_phi_vector.at(i);
		jet_m = jet_m_vector.at(i);
		jet_Et = sqrt(pow(jet_pt,2)+pow(jet_m,2));
		jet_pt_EM = jet_pt_EM_vector.at(i);
		jet_pt_unsubtracted = jet_pt_unsubtracted_vector.at(i);
		jet_pt_seb = jet_pt_SEB_vector.at(i);
		jet_pt_prexcalib = jet_pt_prexcalib_vector.at(i);
		jet_pt_xcalib = jet_pt_xcalib_vector.at(i);
		jet_nConst = jet_nConst_vector.at(i);
		jet_isGood = Is_jet_Good.at(i);
		//jet_NBJ_pT = jet_NBJ_pT_vector.at(i);
		//jet_isMuonIsolated = Is_muon_Isolated.at(i); //TODO if needed
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
			
		//10 GeV reco cut
		
		if (jet_pt_xcalib<10.) continue;
		
		if(_data_switch == 1){
			jet_hasTruth = hasTruth.at(i);
			if (jet_hasTruth){
				truth_jet_eta = truth_jet_eta_vector.at(truth_jet_indices[i]);
				truth_jet_phi = truth_jet_phi_vector.at(truth_jet_indices[i]);
				truth_jet_pt = truth_jet_pt_vector.at(truth_jet_indices[i]);
				truth_jet_m = truth_jet_m_vector.at(truth_jet_indices[i]);
				truth_jet_NBJ_pT = truth_jet_NBJ_pT_vector.at(truth_jet_indices[i]);
			}
			else{
				truth_jet_eta = -100.;
				truth_jet_phi = -100.;
				truth_jet_pt = 0;
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
       tree_performance->Fill();
		
	}
	
	return EL::StatusCode::SUCCESS;
}

EL::StatusCode JetPerformance :: postExecute (){
	return EL::StatusCode::SUCCESS;
}

EL::StatusCode JetPerformance :: finalize (){
	//xAOD::TEvent* event = wk()->xaodEvent();
	
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
