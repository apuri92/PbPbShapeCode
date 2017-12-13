#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include "pPbFragmentation/pPbFragmentation.h"
#include "pPbFragmentation/SEBCorrectorTool.h"
#include "pPbCentrality/pPbMinBiasUtil.h"
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/tools/Message.h"
#include "xAODJet/JetContainer.h"
#include "xAODTrigger/JetRoIContainer.h"
#include "xAODTrigger/JetRoIAuxContainer.h"
#include <TFile.h>
#include <TSystem.h>

using namespace std;
using namespace JetHelperTools;
using namespace TrackHelperTools;
using namespace MTCorrector;

ClassImp(pPbFragmentation)

#define EL_RETURN_CHECK( CONTEXT, EXP )			\
do {							\
  if( ! EXP.isSuccess() ) {				\
    Error( CONTEXT,					\
    XAOD_MESSAGE( "Failed to execute: %s" ),	\
    #EXP );					\
    return EL::StatusCode::FAILURE;			\
    }							\
    } while( false )


pPbFragmentation :: pPbFragmentation ()
{
}

EL::StatusCode pPbFragmentation :: setupJob (EL::Job& job)
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

	//_barcodeMin=0;
	//_barcodeMax=9999;
	//if (_isHerwig) _barcodeMax=1e5;
	  
	std::cout << "[pPbFragmentation() : initialized with jet radius = " << _jet_radius << std::endl;
	//std::cout << "[pPbFragmentation() : initialized barcodes in range " << _barcodeMin << " to " << _barcodeMax << std::endl;
	
	//UEEstimator
	if(!_doSlimTree)
	{
		uee = new UEEstimator();
		uee->ptBkgrThreshold = _trkptBkgrThreshold;
		uee->jetptBkgrThreshold = 20.;
	}
	else
	{
		uee=0;
	}
	
	return EL::StatusCode::SUCCESS;
}

EL::StatusCode pPbFragmentation :: histInitialize ()
{   
	cout << " Setting  histograms" << endl;
	
	h_FCal_Et = new TH1D("h_FCal_Et",";FCal E_{T};N",100,0,10);
	h_RejectionHisto = new TH1D("RejectionHisto","RejectionHisto",9,0,9);
	SetRejectionHistogram(h_RejectionHisto);
	
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
	if(!_doSlimTree)
	{
		hET_ETsub = new TH3D("hET_ETsub","hET_ETsub",200,0,200,100,-5,5,200,-50,50);
		//htrig_reco_pt = new TH2D("htrig_reco_pt","htrig_reco_pt;reco jet p_{T} [GeV]; trigger jet p_{T} [GeV]",100,0,200,100,0,200);
	}
	
	h_triggercounter = new TH2D("h_triggercounter","h_triggercounter",_nTriggers,0,_nTriggers,2,-0.5,1.5);
	SetTrigger_hist(h_triggercounter);
   
	TH3D* temphist_3D = nullptr;
	TH2D* temphist_2D = nullptr;
	TH1D* temphist_1D = nullptr;
	for (int i=0;i<_nTriggers;i++){
		temphist_1D = new TH1D(Form("hjetpt_%.0f",trigger_thresholds.at(i)),Form("hjetpt_%.0f",trigger_thresholds.at(i)),500,0,500);
		hjetpt_trig.push_back(temphist_1D);
		if(!_doSlimTree)
		{
			temphist_2D = new TH2D(Form("hjetz_%.0f",trigger_thresholds.at(i)),Form("hjetz_%.0f",trigger_thresholds.at(i)),300,0,1.5,100,0,500);
			hjetz_trig.push_back(temphist_2D);
			temphist_2D = new TH2D(Form("hjetz_UE_%.0f",trigger_thresholds.at(i)),Form("hjetz_UE_%.0f",trigger_thresholds.at(i)),300,0,1.5,100,0,500);
			hjetz_UE_trig.push_back(temphist_2D);
		}
	}

	if(!_doSlimTree)
	{
		hz = new TH1D("hz","hz",300,0,1.5);
		//hz_UE = new TH1D("hz_UE","hz_UE",300,0,1.5);   
		truth_hjetz = new TH2D("truth_hjetz","truth_hjetz",300,0,1.5,100,0,500);
	
		h_reco_trk_map = new TH3D("h_reco_trk_map","h_reco_trk_map;p_{T};#eta;#phi",ptTrkBinsN, ptTrkBins,etaTrkBinsN,etaTrkBins,phiTrkBinsN,phiTrkBins);
		h_reco_trk_map_nocuts = new TH3D("h_reco_trk_map_nocuts","h_reco_trk_map_nocuts;p_{T};#eta;#phi",ptTrkBinsN, ptTrkBins,etaTrkBinsN,etaTrkBins,phiTrkBinsN,phiTrkBins);
	}
      
	if(!_doSlimTree)
	{
		for (int i=0;i<GetCentralityNBins(_centrality_scheme);i++){
			string i_cent="";
			if (_centrality_scheme>1) i_cent=Form("_%i",i);
			temphist_3D = new TH3D(Form("h_PixHits%s",i_cent.c_str()),Form("h_PixHits%s",i_cent.c_str()),ptTrkBinsN, ptTrkBins,etaTrkBinsN, etaTrkBins,finehitsBinsN,finehitsBins);
			h_PixHits.push_back(temphist_3D); 
			temphist_3D = new TH3D(Form("h_SCTHits%s",i_cent.c_str()),Form("h_SCTHits%s",i_cent.c_str()),ptTrkBinsN, ptTrkBins,etaTrkBinsN, etaTrkBins,finehitsBinsN,finehitsBins);
			h_SCTHits.push_back(temphist_3D); 
			temphist_2D = new TH2D(Form("h_d0%s",i_cent.c_str()),Form("h_d0%s",i_cent.c_str()),ptTrkBinsN,ptTrkBins,d0z0BinsN,d0z0Bins);
			h_d0.push_back(temphist_2D);
			temphist_2D = new TH2D(Form("h_z0sintheta%s",i_cent.c_str()),Form("h_z0sintheta%s",i_cent.c_str()),ptTrkBinsN, ptTrkBins,d0z0BinsN,d0z0Bins);
			h_z0sintheta.push_back(temphist_2D); 
			temphist_3D = new TH3D(Form("h_vx%s",i_cent.c_str()),"h_vx;v_{x};v_{y};v_{z}",200,-2,2,200,-2,2,200,-200,200);
			h_vx.push_back(temphist_3D);
			temphist_3D = new TH3D(Form("h_BL%s",i_cent.c_str()),"h_BL;BL_{x};BL_{y};BL_{z}",200,-1,1,200,-1,1,200,-10,10);
			h_BL.push_back(temphist_3D); 
			wk()->addOutput (h_PixHits.at(i));
			wk()->addOutput (h_SCTHits.at(i));
			wk()->addOutput (h_d0.at(i));
			wk()->addOutput (h_z0sintheta.at(i));
			wk()->addOutput (h_vx.at(i));
			wk()->addOutput (h_BL.at(i));
		}
		
		const int NCuts=8;
		string CutsOn[NCuts]={"Reco","d0","z0sintheta","SI_hits","Pix_holes","SI_holes","Shared_SI_hits","Gen"};
		for (int i=0;i<NCuts;i++){
			temphist_3D = new TH3D(Form("h_cut_flow_%s",CutsOn[i].c_str()),Form("h_cut_flow_%s",CutsOn[i].c_str()),ptTrkBinsN, ptTrkBins,etaTrkBinsN,etaTrkBins,PVBinsN,PVBins);
			h_cut_flow.push_back(temphist_3D);  
			wk()->addOutput (h_cut_flow.at(i));
		}
	}
		
	wk()->addOutput (h_RejectionHisto);
	wk()->addOutput (h_FCal_Et);
	wk()->addOutput (h_triggercounter);
	
	if(!_doSlimTree)
	{
		wk()->addOutput (hET_ETsub);
		//wk()->addOutput (htrig_reco_pt);
		wk()->addOutput (hz);
		//wk()->addOutput (hz_UE);
		wk()->addOutput (truth_hjetz);
		wk()->addOutput (h_reco_trk_map);
		wk()->addOutput (h_reco_trk_map_nocuts);
	}
	
	
	
	
	
	for (int i=0;i<_nTriggers;i++)
	{
		wk()->addOutput (hjetpt_trig.at(i));
		if(!_doSlimTree)
		{
			wk()->addOutput (hjetz_trig.at(i));
			wk()->addOutput (hjetz_UE_trig.at(i));
		}
	}
	
	
	cout << " Histograms  ready, now setting tree" << endl;
	TFile *outputFile = wk()->getOutputFile(_outputName);
	tree_ff = new TTree ("tree_ff", "Reco Jets for FF");
	tree_ff->SetDirectory (outputFile);
	if (_truth_only==0){
		if(!_doSlimTree) tree_ff->Branch("event_n",&event_n,"event_n/I");
		tree_ff->Branch("run_n",&run_n,"run_n/I");
		tree_ff->Branch("lbn_n",&lbn_n,"lbn_n/I");
		tree_ff->Branch("FCalEt",&FCalEt,"FCalEt/F");
		tree_ff->Branch("multiplicity",&multiplicity,"multiplicity/I");
	
		for (int i=0;i<_nTriggers;i++){
			tree_ff->Branch(Form("event_isTriggered_%i",i),&event_isTriggered[i],Form("event_isTriggered_%i/O",i));
			tree_ff->Branch(Form("trigger_prescale_%i",i),&trig_prescale[i],Form("trigger_prescale_%i/F",i));
			tree_ff->Branch(Form("jet_isTriggered_%i",i),&jet_isTriggered[i],Form("jet_isTriggered_%i/I",i));
		}

		//tree_ff->Branch("jet_pt_EM",&jet_pt_EM,"jet_pt_EM/F");
		tree_ff->Branch("jet_pt_xcalib",&jet_pt_xcalib,"jet_pt_xcalib/F");
		if(!_doSlimTree)
		{
			tree_ff->Branch("jet_pt_prexcalib",&jet_pt_prexcalib,"jet_pt_prexcalib/F");
			//tree_ff->Branch("jet_pt_seb",&jet_pt_seb,"jet_pt_seb/F");
			//tree_ff->Branch("jet_pt",&jet_pt,"jet_pt/F");
			if(_data_switch==1)
			{
				tree_ff->Branch("jet_uJER",&jet_uJER,"jet_uJER/F");
				tree_ff->Branch("jet_uJES",&jet_uJES);
			}
		}
		tree_ff->Branch("jet_phi",&jet_phi,"jet_phi/F");
		tree_ff->Branch("jet_eta",&jet_eta,"jet_eta/F");
		tree_ff->Branch("jet_isGood",&jet_isGood,"jet_isGood/I");
		if(!_doSlimTree)
		{
			tree_ff->Branch("jet_m",&jet_m,"jet_m/F");
			tree_ff->Branch("jet_a2_pt",&jet_a2_pt,"jet_a2_pt/F");
			tree_ff->Branch("jet_a2_phi",&jet_a2_phi,"jet_a2_phi/F");
			tree_ff->Branch("jet_a2_eta",&jet_a2_eta,"jet_a2_eta/F");
			tree_ff->Branch("jet_NBJ_pT",&jet_NBJ_pT,"jet_NBJ_pT/F");
			tree_ff->Branch("jet_Muon_pT",&jet_Muon_pT,"jet_Muon_pT/F");
			tree_ff->Branch("jet_electron_pT",&jet_electron_pT,"jet_electron_pT/F"); 
			tree_ff->Branch("jet_centrality",&jet_centrality,"jet_centrality/I");
		}
		tree_ff->Branch("jet_isDummy",&jet_isDummy,"jet_isDummy/I");
		
		tree_ff->Branch("track_pt",&track_pt);
		if(!_doSlimTree) tree_ff->Branch("track_vertex_type",&track_vertex_type);
		tree_ff->Branch("track_phi",&track_phi);
		tree_ff->Branch("track_eta",&track_eta);
		tree_ff->Branch("track_isMuon",&track_isMuon);
		tree_ff->Branch("track_multiJetMatch",&track_multiJetMatch);
		tree_ff->Branch("track_multiJetMatch_B",&track_multiJetMatch_B);
		tree_ff->Branch("track_charge",&track_charge);
		
		for(int n=0;n<_nTrkSelTools;++n)
		{
			tree_ff->Branch(Form("track_passedTrkSelTool_%d",n),&track_passedTrkSelTool[n]);
		}
		
		tree_ff->Branch("track_nPixHits",&track_nPixHits);
		if(!_doSlimTree || _centrality_scheme<30) tree_ff->Branch("track_nPixHoles",&track_nPixHoles);
		//tree_ff->Branch("track_nPixDeadS",&track_nPixDeadS); 
		tree_ff->Branch("track_nShPixH",&track_nShPixH);
		tree_ff->Branch("track_nSCTHits",&track_nSCTHits);
		if(!_doSlimTree || _centrality_scheme<30) tree_ff->Branch("track_nSCTHoles",&track_nSCTHoles);
		//tree_ff->Branch("track_nSCTDeadS",&track_nSCTDeadS);
		tree_ff->Branch("track_nShSCTH",&track_nShSCTH);
		tree_ff->Branch("track_nTRTHits",&track_nTRTHits);
		tree_ff->Branch("track_nIBLHits",&track_nIBLHits);
		tree_ff->Branch("track_expIBLHits",&track_expIBLHits);
		tree_ff->Branch("track_nBLHits",&track_nBLHits);
		tree_ff->Branch("track_expBLHits",&track_expBLHits);
		
		tree_ff->Branch("track_d0",&track_d0);
		tree_ff->Branch("track_z0sintheta",&track_z0sintheta);
		tree_ff->Branch("track_ed0",&track_ed0);
		tree_ff->Branch("track_ez0sintheta",&track_ez0sintheta);
		if(_data_switch==1){
			tree_ff->Branch("track_mc_pt",&track_mc_pt);
			tree_ff->Branch("track_mc_eta",&track_mc_eta);
			tree_ff->Branch("track_mc_phi",&track_mc_phi);
			tree_ff->Branch("track_mc_probability",&track_mc_probability);
			tree_ff->Branch("track_mc_pdg",&track_mc_pdg);
			tree_ff->Branch("track_mc_barcode",&track_mc_barcode);
			tree_ff->Branch("track_mc_type",&track_mc_type);
			tree_ff->Branch("track_mc_charge",&track_mc_charge);
		}
		if(!_doSlimTree) tree_ff->Branch("track_z",&track_z);

		if(_data_switch==1){
			tree_ff->Branch("truth_jet_pt",&truth_jet_pt,"truth_jet_pt/F");
			tree_ff->Branch("truth_jet_phi",&truth_jet_phi,"truth_jet_phi/F");
			tree_ff->Branch("truth_jet_eta",&truth_jet_eta,"truth_jet_eta/F");
			if(!_doSlimTree)
			{
				tree_ff->Branch("truth_jet_m",&truth_jet_m,"truth_jet_m/F");
				tree_ff->Branch("truth_reco_jet_dR",&truth_reco_jet_dR,"truth_reco_jet_dR/F");
			}
		}
		
		//TRUTH tree
		truth_tree_ff = new TTree ("truth_tree_ff", "Reco Jets for FF");
  		truth_tree_ff->SetDirectory (outputFile);
		if(_data_switch==1){
			if(!_doSlimTree) truth_tree_ff->Branch("event_n",&event_n,"event_n/I");
			truth_tree_ff->Branch("run_n",&run_n,"run_n/I");
			truth_tree_ff->Branch("truth_jet_pt",&truth_jet_pt,"truth_jet_pt/F");
			truth_tree_ff->Branch("truth_jet_phi",&truth_jet_phi,"truth_jet_phi/F");
			truth_tree_ff->Branch("truth_jet_eta",&truth_jet_eta,"truth_jet_eta/F");
			truth_tree_ff->Branch("truth_jet_isDummy",&truth_jet_isDummy,"truth_jet_isDummy/I");	
			if(!_doSlimTree) 
			{
				truth_tree_ff->Branch("truth_jet_m",&truth_jet_m,"truth_jet_m/F");
				truth_tree_ff->Branch("truth_jet_Muon_pT",&truth_jet_Muon_pT,"truth_jet_Muon_pT/F");
				truth_tree_ff->Branch("truth_jet_electron_pT",&truth_jet_electron_pT,"truth_jet_electron_pT/F");
				truth_tree_ff->Branch("truth_jet_NBJ_pT",&truth_jet_NBJ_pT,"truth_jet_NBJ_pT/F");
			}
			
			truth_tree_ff->Branch("truth_track_pt",&truth_track_pt);
			truth_tree_ff->Branch("truth_track_phi",&truth_track_phi);
			truth_tree_ff->Branch("truth_track_eta",&truth_track_eta);
			if(!_doSlimTree) truth_tree_ff->Branch("truth_track_z",&truth_track_z);
			truth_tree_ff->Branch("truth_track_pdg",&truth_track_pdg);
			truth_tree_ff->Branch("truth_track_barcode",&truth_track_barcode);
			truth_tree_ff->Branch("truth_track_charge",&truth_track_charge);
			truth_tree_ff->Branch("truth_track_multiJetMatch",&truth_track_multiJetMatch);
			truth_tree_ff->Branch("truth_track_multiJetMatch_B",&truth_track_multiJetMatch_B);

			if(!_doSlimTree) truth_tree_ff->Branch("jet_centrality",&jet_centrality,"jet_centrality/I");
			truth_tree_ff->Branch("FCalEt",&FCalEt,"FCalEt/F");
			truth_tree_ff->Branch("multiplicity",&multiplicity,"multiplicity/I");
		}
	}
	
	cout << " tree  ready" << endl; 
	
	return EL::StatusCode::SUCCESS;
}

EL::StatusCode pPbFragmentation :: fileExecute ()
{
	// Here you do everything that needs to be done exactly once for every
	// single file, e.g. collect a list of all lumi-blocks processed
	return EL::StatusCode::SUCCESS;
}

EL::StatusCode pPbFragmentation :: changeInput (bool firstFile)
{
	return EL::StatusCode::SUCCESS;
}

EL::StatusCode pPbFragmentation :: initialize ()
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
		EL_RETURN_CHECK("initialize()", m_trigDecisionTool->initialize() );
		
		cout << "Adding following " << _nTriggers << " triggers: ";
		for (int i=0;i<_nTriggers;i++){
			cout << trigger_chains.at(i) << ", ";
			_chainGroup.push_back(m_trigDecisionTool->getChainGroup(trigger_chains.at(i)));
		}	
		cout << endl << "Initialize triggers finished" << endl;
	}
	
	//MuonSelectionTool 
	m_muonSelection = new CP::MuonSelectionTool("MuonSelection");
	m_muonSelection->msg().setLevel( MSG::ERROR ); 
	m_muonSelection->setProperty( "MaxEta", 2.5 ); 
	m_muonSelection->setProperty( "MuQuality", 2); // 0-tight, 1-medium, 2-loose, 3-very loose
	EL_RETURN_CHECK("initialize()", m_muonSelection->initialize() );
	 
	//Muon corection tool (MC only) 
	m_muonCalibrationAndSmearingTool = new CP::MuonCalibrationAndSmearingTool( "MuonCorrectionTool" ); 
	EL_RETURN_CHECK("initialize()", m_muonCalibrationAndSmearingTool->initialize() ); 

	//Track Selection Tool
	//_nTrkSelTools=2; - set in run_chain_pPbFragmentation.cxx
	m_trackSelectionTool[0] = new InDet::InDetTrackSelectionTool("InDetTrackSelectionTool_0");
	EL_RETURN_CHECK("initialize()",m_trackSelectionTool[0]->setProperty("CutLevel","HITight"));

	m_trackSelectionTool[1] = new InDet::InDetTrackSelectionTool("InDetTrackSelectionTool_1");
	// tohle je skoro TightPrimary, ale mame jiny cut na IBL+BL hity (to je to useMinBiasInnermostLayersCut) a jsou ignorovany shared hity (to nejak neslo vypnout, proto jsem to tu zkoporoval)
	EL_RETURN_CHECK("initialize()",m_trackSelectionTool[1]->setProperty("CutLevel","NoCut"));
	EL_RETURN_CHECK("initialize()",m_trackSelectionTool[1]->setProperty("maxAbsEta",2.5));
	EL_RETURN_CHECK("initialize()",m_trackSelectionTool[1]->setProperty("minNSiHits", 9));
	EL_RETURN_CHECK("initialize()",m_trackSelectionTool[1]->setProperty("maxNSiHoles", 2));
	EL_RETURN_CHECK("initialize()",m_trackSelectionTool[1]->setProperty("maxNPixelHoles", 0));
	EL_RETURN_CHECK("initialize()",m_trackSelectionTool[1]->setProperty("minEtaForStrictNSiHitsCut", 1.65));
	EL_RETURN_CHECK("initialize()",m_trackSelectionTool[1]->setProperty("minNSiHitsAboveEtaCutoff", 11));
	EL_RETURN_CHECK("initialize()",m_trackSelectionTool[1]->setProperty("useMinBiasInnermostLayersCut",1)); // s timhle asi projdou tracky, ktere nemaji hity v IBL ani v BL a ktere nemaji expected hity v BL ani v IBL; podle me takove tracky nechceme (protoze ma platit nIBL+nBL>0), ale jine cuty by se o to mohly postarat...
	
	for(int n=0;n<_nTrkSelTools;++n)
	{
		EL_RETURN_CHECK("initialize()",m_trackSelectionTool[n]->initialize());
	}
	
	
	//Corrections for track pT sagitta bias
	//m_trkBiasingTool = new InDet::InDetTrackBiasingTool("InDetTrackBiasingTool");
	//EL_RETURN_CHECK("initialize()", m_trkBiasingTool->initialize() );

	if(!_doSlimTree)
	{	
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
			 
		// initialize the tools 
		EL_RETURN_CHECK("initialize()", m_LHToolTight2015->initialize() ); 	
	}
	else
	{
		m_EgammaCalibrationAndSmearingTool=0;
		m_LHToolTight2015=0;
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

	if(!_doSlimTree)
	{	
		//Unceratinty
		//JES
		jesProv = new JetUncertaintiesTool("JESProvider");
		EL_RETURN_CHECK("initialize()",jesProv->setProperty("JetDefinition","AntiKt4EMTopo"));
		EL_RETURN_CHECK("initialize()",jesProv->setProperty("MCType","MC15"));
		EL_RETURN_CHECK("initialize()",jesProv->setProperty("ConfigFile","JES_2015/ICHEP2016/JES2015_AllNuisanceParameters.config"));
	
		//JER
		jerTool = new JERTool("JERTool");
		smearTool = new JERSmearingTool("JERSmearingTool");
		EL_RETURN_CHECK("initialize()",jerTool->setProperty("PlotFileName", "JetResolution/Prerec2015_xCalib_2012JER_ReducedTo9NP_Plots_v2.root") );
		EL_RETURN_CHECK("initialize()",jerTool->setProperty("CollectionName", "AntiKt4EMTopoJets") );

		// Configure the JERSmearingTool
		smearTool->msg().setLevel(MSG::DEBUG);
		ToolHandle<IJERTool> jerHandle(jerTool->name());
		EL_RETURN_CHECK("initialize()",smearTool->setProperty("JERTool", jerHandle) );
		EL_RETURN_CHECK("initialize()",smearTool->setProperty("ApplyNominalSmearing", false) );
		EL_RETURN_CHECK("initialize()",smearTool->setProperty("isMC", true) );
		EL_RETURN_CHECK("initialize()",smearTool->setProperty("SystematicMode", "Full") );

		// Initialize the tools
		EL_RETURN_CHECK("initialize()",jerTool->initialize() );
		EL_RETURN_CHECK("initialize()",smearTool->initialize() );
	
		// Initialize the tool
		EL_RETURN_CHECK("initialize()",jesProv->initialize());
	}
	else
	{
		jesProv=0;
		jerTool=0;
		smearTool=0;
	}
	
	cout << " Initialization done" << endl;
	return EL::StatusCode::SUCCESS;
}

//Loop over events
EL::StatusCode pPbFragmentation :: execute (){

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
	bool keep = true;
	h_RejectionHisto->Fill(0.5);
		
	//---------------------------
	//     Event information
	//--------------------------- 
		
	const xAOD::EventInfo* eventInfo = 0;
	EL_RETURN_CHECK("execute",event->retrieve( eventInfo, "EventInfo"));
		
	event_n = eventInfo->eventNumber();
	run_n = eventInfo->runNumber();
	lbn_n = eventInfo->lumiBlock();
		
	// check if the event is data or MC
	bool isMC = false;
	if(eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) )
	{
		isMC = true;
		_data_switch=1; 
	}
	else
	{
		const xAOD::TruthParticleContainer * particles = 0;
		if( event->xAOD::TVirtualEvent::retrieve(particles, "TruthParticles", true) )
		{
			// this is overlay
			isMC = true;
			_data_switch=1; 
		}
		else
		{
			_data_switch=0;
		}
	} 
	
	
	FCalEt = 0;
	int cent_bin=0;
	if (_centrality_scheme>1) {
		//Centrality
		const xAOD::HIEventShapeContainer* calos = 0;
		EL_RETURN_CHECK("execute()",event->retrieve( calos, "CaloSums" ));
		FCalEt = 0;
		int x = 0;
		xAOD::HIEventShapeContainer::const_iterator calo_itr = calos->begin();
		xAOD::HIEventShapeContainer::const_iterator calo_end = calos->end();
		for( ; calo_itr != calo_end; ++calo_itr ) {
			if (x == 5) {
				FCalEt = ((*calo_itr)->et() * 0.001 * 0.001 );
			}
			x++;
		}
		cent_bin = GetCentralityBin(_centrality_scheme, FCalEt, isMC);
		h_FCal_Et->Fill(FCalEt);
	}  
	if (cent_bin < 0) {
		if(cent_bin==-2) Error("execute()", "Unknown centrality scheme" );
		h_RejectionHisto->Fill(1.5);
		keep = false;
	}
	
	// GRL
	if(!isMC){ 
		if(!m_grl->passRunLB(*eventInfo)){
			h_RejectionHisto->Fill(2.5);
			keep = false;
		}
	}
	
	//Vertex requirement
	const xAOD::VertexContainer * vertices = 0;
	if ( !event->retrieve( vertices, "PrimaryVertices" ).isSuccess() ){
		Error("execute()", "Failed to retrieve VertexContainer container. Exiting." );
		keep = false;
	}
		
	if(vertices->size()<2) {
		h_RejectionHisto->Fill(3.5);
		keep = false;
	}	
	xAOD::VertexContainer::const_iterator vtx_itr = vertices->begin();
	xAOD::VertexContainer::const_iterator vtx_end = vertices->end();
	// find primary vertex
	const xAOD::Vertex* primaryVertex = 0;
	for(;vtx_itr!=vtx_end;++vtx_itr)
	{
		if((*vtx_itr)->vertexType()==xAOD::VxType::PriVtx) {
			primaryVertex = (*vtx_itr);
			//if(cent_bin>=0 && !_doSlimTree) h_vx.at(cent_bin)->Fill(primaryVertex->x(),primaryVertex->y(),primaryVertex->z());
			if(cent_bin>=0 && !_doSlimTree) h_vx.at(cent_bin)->Fill(0.0,0.0,primaryVertex->z());
			break;
		}
	}

	if(primaryVertex)
	{
		if (fabs(primaryVertex->z())>150.){
			h_RejectionHisto->Fill(5.5);
			keep = false;
		}
	}
	else
	{
		h_RejectionHisto->Fill(5.5);
		keep = false;
	}
	
	//TODO rapidity gap?
	
	//DAQ errors 
	if(!isMC){
		if(   (eventInfo->errorState(xAOD::EventInfo::LAr)==xAOD::EventInfo::Error ) || (eventInfo->errorState(xAOD::EventInfo::Tile)==xAOD::EventInfo::Error ) || (eventInfo->errorState(xAOD::EventInfo::SCT)==xAOD::EventInfo::Error ) || (eventInfo->isEventFlagBitSet(xAOD::EventInfo::Core, 18) ) ){
			h_RejectionHisto->Fill(4.5);			
			keep = false;
		}
	}
	
	if (!keep) return EL::StatusCode::SUCCESS; // go to the next event 
	 
	h_RejectionHisto->Fill(7.5);
	
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
			else h_RejectionHisto->Fill(8.5);
		}
	}  
	
	
	//Tracks
	const xAOD::TrackParticleContainer* recoTracks = 0;
	if ( !event->retrieve( recoTracks, "InDetTrackParticles" ).isSuccess() ){
		Error("execute()", "Failed to retrieve Reconstructed Track container. Exiting." );
		return EL::StatusCode::FAILURE;
	}
	multiplicity=0;
	for (const auto& trk : *recoTracks)
	{
		if(trk->pt()>1000.) ++multiplicity;
	}
			
	
	//Muons
	const xAOD::MuonContainer * muons = 0;
	EL_RETURN_CHECK("execute()",event->retrieve(muons,"Muons"));
	
	// create a shallow copy of the muons container 
	std::pair< xAOD::MuonContainer*, xAOD::ShallowAuxContainer* > muons_shallowCopy = xAOD::shallowCopyContainer( *muons ); 
	
	//if (isMC){      
		xAOD::MuonContainer::iterator muonSC_itr = (muons_shallowCopy.first)->begin(); 
		xAOD::MuonContainer::iterator muonSC_end = (muons_shallowCopy.first)->end(); 
		for( ; muonSC_itr != muonSC_end; ++muonSC_itr ) { 
			if(m_muonCalibrationAndSmearingTool->applyCorrection(**muonSC_itr) == CP::CorrectionCode::Error){ // apply correction and check return code 
			Error("execute()", "MuonCalibrationAndSmearingTool returns Error CorrectionCode"); 
			}  
		} 
	//}
	
	//electrons
	const xAOD::ElectronContainer* electrons = 0; 
	EL_RETURN_CHECK("execute()",event->retrieve( electrons, "Electrons" )); 
	
	std::pair< xAOD::ElectronContainer*, xAOD::ShallowAuxContainer* > electrons_shallowCopy; 
	electrons_shallowCopy = xAOD::shallowCopyContainer( *electrons );
	
	//Jet vectors
	vector<float> jet_pt_EM_vector, jet_pt_SEB_vector,jet_pt_prexcalib_vector,jet_pt_xcalib_vector,jet_phi_vector,jet_eta_vector,jet_m_vector,jet_uJER_vector;
	vector<vector<float>> jet_uJES_vector;
	vector<int> Is_jet_Good, Is_dummyJet;
	vector<float> truth_jet_eta_vector,truth_jet_m_vector,truth_jet_phi_vector,truth_jet_pt_vector;
	vector<int> truth_jet_indices, truth_jet_isDummy_vector;
	vector<int> hasTruth;
	vector<int> isTriggered[_nTriggers];
	vector<float> jet_NBJ_pT_vector,truth_jet_NBJ_pT_vector;
	vector<float> trig_EF_jet_pt, trig_EF_jet_phi, trig_EF_jet_eta;
	vector<float> antikt2_pt,antikt2_phi,antikt2_eta;
	
	// Clear vectors
	hasTruth.clear();
	
	truth_jet_indices.clear();
	truth_jet_pt_vector.clear();
	truth_jet_eta_vector.clear();
	truth_jet_phi_vector.clear();
	truth_jet_m_vector.clear();
	jet_uJER_vector.clear();
	jet_uJES_vector.clear();
	truth_jet_isDummy_vector.clear();
	
	jet_pt_EM_vector.clear();
	jet_pt_SEB_vector.clear();
	jet_pt_prexcalib_vector.clear();
	jet_pt_xcalib_vector.clear();
	jet_phi_vector.clear();
	jet_eta_vector.clear();
	jet_m_vector.clear();
	
	Is_jet_Isolated.clear();
	jet_NBJ_pT_vector.clear();
	truth_jet_NBJ_pT_vector.clear();
	InJet_muon_pT.clear(); 
	InJet_electron_pT.clear();
	Is_jet_Good.clear();
	
	for (int j=0;j<_nTriggers;j++){
		isTriggered[j].clear();
	}
	
	truth_reco_jet_dR=0.0;
	
	//***** HLT Jets *****
	if(_data_switch==0 && !_isMB) {   
		const xAOD::JetContainer * hlt_jet = 0;
		string trigger_container="HLT_xAOD__JetContainer_" + _trigger_collection;
		  
		if( ! event->retrieve( hlt_jet, trigger_container.c_str()).isSuccess() ) {
			Error("execute()", Form("failed to retrieve %s", trigger_container.c_str()));
		}
		
		const xAOD::Jet *hlt_jet_obj = 0;
		for(unsigned int HLTjet_itr = 0; HLTjet_itr < hlt_jet->size(); ++HLTjet_itr) {
			hlt_jet_obj = hlt_jet->at(HLTjet_itr);
			trig_EF_jet_pt.push_back(  hlt_jet_obj->pt() );
			trig_EF_jet_phi.push_back( hlt_jet_obj->phi() );
			trig_EF_jet_eta.push_back( hlt_jet_obj->eta() );   
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
			truth_jet_m_vector.push_back( truth_jet_obj->m() * 0.001 );  
			truth_jet_isDummy_vector.push_back(0);
		}
		if(!_doSlimTree) truth_jet_NBJ_pT_vector = MTCorrector::GetIsolation(truth_jet_pt_vector,truth_jet_eta_vector,truth_jet_phi_vector,_jet_radius);
	}
	
	//***** Reco jets*****
	if (_truth_only==0) {
		xAOD::TStore store; //For calibration
		const xAOD::JetContainer* jets = 0;
		EL_RETURN_CHECK("execute()",event->retrieve( jets, _reco_jet_collection.c_str() ));
		
		xAOD::JetContainer::const_iterator jet_itr = jets->begin();
		xAOD::JetContainer::const_iterator jet_end = jets->end();
		
	/*	xAOD::JetContainer* updatedjets = new xAOD::JetContainer();
		xAOD::AuxContainerBase* updatedjetsAux = new xAOD::AuxContainerBase();
		updatedjets->setStore( updatedjetsAux );
	*/	
		int jet_counter=0;
		for( ; jet_itr != jet_end; ++jet_itr ) {
			jet_counter++;

			xAOD::Jet* newjet = new xAOD::Jet();
			newjet->makePrivateStore( **jet_itr );
			//updatedjets->push_back( newjet );
			
			const xAOD::JetFourMom_t jet_4mom_def = newjet->jetP4();
			xAOD::JetFourMom_t jet_4mom;
			if (_centrality_scheme==30) jet_4mom = newjet->jetP4("JetSubtractedScaleMomentum");
			else jet_4mom = newjet->jetP4("JetEMScaleMomentum");
			float def_jet_pt  = (jet_4mom_def.pt() * 0.001);
			float unsubtracted_jet_pt  = (jet_4mom.pt() * 0.001);  
			float uncalib_jet_pt  = (jet_4mom.pt() * 0.001);
			
			if (_reco_jet_collection.find("HI") != std::string::npos) {
				const xAOD::JetFourMom_t jet_4mom_unsubtracted = newjet->jetP4("JetUnsubtractedScaleMomentum");
				unsubtracted_jet_pt  = (jet_4mom_unsubtracted.pt() * 0.001);   
			}
			
			if(!_doSlimTree) hET_ETsub->Fill(def_jet_pt,jet_4mom_def.eta(),unsubtracted_jet_pt-uncalib_jet_pt);
			
			if (_centrality_scheme==30) newjet->setJetP4("JetConstitScaleMomentum",jet_4mom);
			newjet->setJetP4("JetPileupScaleMomentum", jet_4mom);	
			
			EL_RETURN_CHECK("execute()", m_jetCalibration->applyCalibration( *newjet ) );
			float jet_pt_prexcalib = (newjet->pt() * 0.001);
			
			//Insitu corrections, need to be done at GSW scale
			const xAOD::JetFourMom_t jet_4mom_xcalib = newjet->jetP4();
			newjet->setJetP4("JetGSCScaleMomentum", jet_4mom_xcalib);
                               
			if (_data_switch==0) EL_RETURN_CHECK("execute()", m_jetCalibration_insitu->applyCalibration( *newjet ) );
			
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
			
			if (fabs(jet_eta)>3.0) continue;
			
			//Jet quality moment
			if( !m_jetCleaning->accept( **jet_itr )) Is_jet_Good.push_back(0);
			else Is_jet_Good.push_back(1);
			
			jet_pt_EM_vector.push_back(uncalib_jet_pt);
			jet_pt_SEB_vector.push_back(jet_pt);
			jet_pt_prexcalib_vector.push_back(jet_pt_prexcalib);
			jet_pt_xcalib_vector.push_back(jet_pt);
			
			jet_phi_vector.push_back(jet_phi);
			jet_eta_vector.push_back(jet_eta);
			jet_m_vector.push_back(jet_m); 
			Is_dummyJet.push_back(0);
			
			//Uncertainties
			if (isMC)
			{
				// Get the resolution in MC and data
				//double resMC = jerTool->getRelResolutionMC(newjet); -- not used
				// Get the resolution uncertainty
				JER::Uncert errType=JER::JER_CROSS_CALIB_ALL; // JER probably no needed at all
				double uncert = -1.0;
				if(!_doSlimTree) uncert = jerTool->getUncertainty(newjet, errType);
				jet_uJER_vector.push_back(uncert);
								
				vector<float> uJES;
				if(!_doSlimTree)
				{
					for (size_t iComp = 0; iComp < jesProv->getNumComponents(); ++iComp)
					{
						double unc = -1.0;
						unc = jesProv->getUncertainty(iComp,*newjet);
						uJES.push_back(unc); // 1+unc for an upward shift, use 1-unc instead for a downward shift
					}
				}
				else
				{
					uJES.push_back(-1.0);
				}
				jet_uJES_vector.push_back(uJES);	
			}

			delete newjet;
		}
		
		if(!_doSlimTree)
		{
			jet_NBJ_pT_vector = MTCorrector::GetIsolation(jet_pt_xcalib_vector,jet_eta_vector,jet_phi_vector,_jet_radius);
			InJet_muon_pT = FilterJetsByMuons(jet_eta_vector,jet_phi_vector,muons_shallowCopy,_dR_max); 
			InJet_electron_pT = FilterJetsByElectrons(jet_eta_vector,jet_phi_vector,electrons_shallowCopy,_dR_max);
		}
		
		// add dummy jet
		if(_isMB)
		{
			Is_jet_Good.push_back(0);
			jet_pt_EM_vector.push_back(-1.);
			jet_pt_SEB_vector.push_back(-1.);
			jet_pt_prexcalib_vector.push_back(-1.);
			jet_pt_xcalib_vector.push_back(-1.);
			 
			jet_phi_vector.push_back(0.);
			jet_eta_vector.push_back(-99.);
			jet_m_vector.push_back(-1.); 
			Is_dummyJet.push_back(1);
			
			jet_NBJ_pT_vector.push_back(-1.);
			InJet_muon_pT.push_back(0);
			InJet_electron_pT.push_back(0);
			
			if (isMC)
			{
				jet_uJER_vector.push_back(-1);
				vector<float> uJES;
				if(!_doSlimTree)
				{
					for (size_t iComp = 0; iComp < jesProv->getNumComponents(); ++iComp)
					{
					   uJES.push_back(-1);
					}
				}
				else
				{
					uJES.push_back(-1);
				}
				jet_uJES_vector.push_back(uJES);
			}
			
		}
		
		//Get the R=0.2 jets
		if (_reco_jet_collection.find("HI") != std::string::npos) {
			const xAOD::JetContainer* a2_jets = 0;
			EL_RETURN_CHECK("execute()",event->retrieve( a2_jets, "AntiKt2HIJets" ));

			// loop over the jets in the container
			xAOD::JetContainer::const_iterator a2_jet_itr = a2_jets->begin();
			xAOD::JetContainer::const_iterator a2_jet_end = a2_jets->end();
			
			for( ; a2_jet_itr != a2_jet_end; ++a2_jet_itr ) {
				antikt2_pt.push_back((*a2_jet_itr)->pt());
				antikt2_eta.push_back((*a2_jet_itr)->eta());
				antikt2_phi.push_back((*a2_jet_itr)->phi());
			} 
		}
		
		//Data only
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
					/*
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
					*/
					//No matching to trigger jets
					for (int k=0;k<_nTriggers;k++){			 
						if(event_isTriggered[k]){
							if(jet_pt > jet_pt_trig[k][0] && jet_pt < jet_pt_trig[k][1]) {				 
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
		
		//MC only 	  
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
		
		for(unsigned int i=0; i<jet_pt_xcalib_vector.size(); i++){
			
			int triggered_at_least_once=0;
			if (_data_switch==0){
				for (int j=0;j<_nTriggers;j++){
					jet_isTriggered[j] = isTriggered[j].at(i);
					if(jet_isTriggered[j]){
						triggered_at_least_once=1;
						//break;
					}
				}
				
				//skip ther rest
				if(!triggered_at_least_once) continue;
			}
			
			jet_pt = jet_pt_xcalib_vector.at(i);
			jet_eta = jet_eta_vector.at(i);		 
			jet_phi = jet_phi_vector.at(i);
			jet_m = jet_m_vector.at(i);
			jet_pt_EM = jet_pt_EM_vector.at(i);
			jet_pt_seb = jet_pt_SEB_vector.at(i);
			jet_pt_prexcalib = jet_pt_prexcalib_vector.at(i);
			jet_pt_xcalib = jet_pt_xcalib_vector.at(i);
			jet_isGood = Is_jet_Good.at(i);
			if(!_doSlimTree) 
			{
				jet_NBJ_pT = jet_NBJ_pT_vector.at(i);
				jet_Muon_pT = InJet_muon_pT.at(i);
				jet_electron_pT = InJet_electron_pT.at(i);
			}
			jet_centrality=cent_bin;
			jet_isDummy = Is_dummyJet.at(i);
			if (isMC && !_doSlimTree)
			{
				jet_uJER= jet_uJER_vector.at(i);
				jet_uJES= jet_uJES_vector.at(i);
			}			 
			
			
			//10 GeV reco cut
			if (jet_pt_xcalib<_pTjetCut && !jet_isDummy) continue;
			
			if(_data_switch == 1){
				jet_hasTruth = hasTruth.at(i);
				if (jet_hasTruth){
					truth_jet_eta = truth_jet_eta_vector.at(truth_jet_indices[i]);
					truth_jet_phi = truth_jet_phi_vector.at(truth_jet_indices[i]);
					truth_jet_pt = truth_jet_pt_vector.at(truth_jet_indices[i]);
					truth_jet_m = truth_jet_m_vector.at(truth_jet_indices[i]);
				}
				else{
					truth_jet_eta = -100.;
					truth_jet_phi = -100.;
					truth_jet_pt = 0;
				}
			}
			
			//Matching to a2 jest
			jet_a2_pt = -999; jet_a2_eta = -999; jet_a2_phi = -999;
			float dR_min=999;
			for(unsigned int j=0; j<antikt2_pt.size(); j++){
				float a2_jet_pt = antikt2_pt.at(j)/1000.;
				float a2_jet_eta = antikt2_eta.at(j);
				float a2_jet_phi = antikt2_phi.at(j);
				float R = DeltaR(a2_jet_phi,a2_jet_eta,jet_phi,jet_eta);
				if(R < 0.3 && R<dR_min){
					jet_a2_pt=a2_jet_pt;
					jet_a2_eta=a2_jet_eta;
					jet_a2_phi=a2_jet_phi;
					dR_min = R;
				}
			}
			
			
			track_eta.clear();
			track_phi.clear();
			track_pt.clear();
			track_vertex_type.clear();
			track_isMuon.clear();
			track_multiJetMatch.clear();
			track_multiJetMatch_B.clear();
			track_charge.clear();
			for(int n=0;n<_nTrkSelTools;++n)
			{
				track_passedTrkSelTool[n].clear();
			}
			
			track_nPixHits.clear();
			track_nPixHoles.clear();
			//track_nPixDeadS.clear();
			track_nShPixH.clear();
			track_nSCTHits.clear();
			track_nSCTHoles.clear();
			//track_nSCTDeadS.clear();
			track_nShSCTH.clear();
			track_nTRTHits.clear();
			track_nIBLHits.clear();
			track_expIBLHits.clear();
			track_nBLHits.clear();
			track_expBLHits.clear();
			
			track_d0.clear();
			track_z0sintheta.clear();
			track_ed0.clear();
			track_ez0sintheta.clear();
			
			track_z.clear();
			track_mc_pt.clear();
			track_mc_eta.clear();
			track_mc_phi.clear();
			track_mc_probability.clear();
			track_mc_pdg.clear();
			track_mc_barcode.clear();
			track_mc_type.clear();
			track_mc_charge.clear();
			
			for (const auto& trk : *recoTracks) {
				//get the tracks....
				
				/*if(trk->pt() > 100. && trk->eta()>0 && trk->eta()<0.7 && trk->phi()<0 && trk->phi()>-1.5)
				{
					std::cout<<"run "<<run_n<<std::endl;
					xAOD::TrackParticle * newTrk = nullptr;
					m_trkBiasingTool->correctedCopy( *trk, newTrk );
					std::cout<<" old pT = "<<trk->pt()/1000.<<" eta "<<trk->eta()<<" phi "<<trk->phi()<<" ..  new pT "<<newTrk->pt()/1000.<<"  "<<trk->pt()/1000./(1+trk->pt()/1000000.*0.922689*trk->charge())<<std::endl;
				}*/

				float pt = trk->pt()/1000.;
				float eta = trk->eta();
				float phi = trk->phi();
				
				if (fabs(eta) > 2.5) continue;
				if (pt < _pTtrkCut) continue;
				
				// reject TRT-seeded tracks
				if(_centrality_scheme<30)
				{
					std::bitset < xAOD::NumberOfTrackRecoInfo > author = trk->patternRecoInfo();
					if(author[xAOD::TRTSeededTrackFinder] || author[xAOD::TRTStandalone]) continue;
				}
				
				//track parameters
				int nPixHits = trk->auxdata< unsigned char >("numberOfPixelHits") + trk->auxdata< unsigned char >("numberOfPixelDeadSensors");
				int nPixHoles = trk->auxdata< unsigned char >("numberOfPixelHoles");
				//int nPixDeadS = trk->auxdata< unsigned char >("numberOfPixelDeadSensors");
				int nShPixH= trk->auxdata< unsigned char >("numberOfPixelSharedHits");
				int nSCTHits = trk->auxdata< unsigned char >("numberOfSCTHits") + trk->auxdata< unsigned char >("numberOfSCTDeadSensors");
				int nSCTHoles = trk->auxdata< unsigned char >("numberOfSCTHoles");
				//int nSCTDeadS = trk->auxdata< unsigned char >("numberOfSCTDeadSensors");
				int nShSCTH = trk->auxdata< unsigned char >("numberOfSCTSharedHits");
				int nTRTHits = trk->auxdata< unsigned char >("numberOfTRTHits");  
				int nIBLHits = trk->auxdata< unsigned char >("numberOfInnermostPixelLayerHits");
				int nBLHits = trk->auxdata< unsigned char >("numberOfNextToInnermostPixelLayerHits");
				int expIBLHits = trk->auxdata< unsigned char >("expectInnermostPixelLayerHit");
				int expBLHits = trk->auxdata< unsigned char >("expectNextToInnermostPixelLayerHit");
				
				double d0 = trk->d0();
				double theta = trk->theta();
				double z0pv=(trk->z0()+trk->vz()-(*vtx_itr)->z())*sin(theta);	// pp: trk->z0() - w.r.t. BS
				
				float ed0=sqrt(trk->definingParametersCovMatrix()(0,0));
				float ez0pv=sqrt( (trk->definingParametersCovMatrix()(1,1))*pow(sin(theta),2) +
				                  (trk->definingParametersCovMatrix()(3,3))*pow(z0pv*cos(theta),2)+
				                2*(trk->definingParametersCovMatrix()(1,3))*fabs(sin(theta)*z0pv*cos(theta)) );
				if(ed0==0.0) ed0=1e-10;
				if(ez0pv==0.0) ez0pv=1e-10;
							
				// check if it's muon
				
				float isMuon=-1;
				//xAOD::MuonContainer::const_iterator muon_itr = muons->begin();
				//xAOD::MuonContainer::const_iterator muon_end = muons->end();
				
				xAOD::MuonContainer::const_iterator muon_itr = ( muons_shallowCopy.first)->begin();
				xAOD::MuonContainer::const_iterator muon_end = ( muons_shallowCopy.first)->end();

				for( ; muon_itr != muon_end; ++muon_itr )
				{
					ElementLink< xAOD::TrackParticleContainer > trkmulink = (*muon_itr)->auxdata<ElementLink< xAOD::TrackParticleContainer > >("inDetTrackParticleLink");
					if(trkmulink.isValid() && m_muonSelection->accept(**muon_itr) )
					{
						if(TMath::Abs((*trkmulink)->pt() - trk->pt())<1e-6 && TMath::Abs((*trkmulink)->eta() - trk->eta())<1e-6 && TMath::Abs((*trkmulink)->phi() - trk->phi())<1e-6)
						{
							isMuon=(*muon_itr)->pt()/1000.;
							break;
						}
					}
				}
				
				//Check the  track-vertex association
				/*				
				NoVtx   = 0, ///< Dummy vertex. TrackParticle was not used in vertex fit
				PriVtx  = 1, ///< Primary vertex
				SecVtx  = 2, ///< Secondary vertex
				PileUp  = 3, ///< Pile-up vertex
				ConvVtx = 4, ///< Conversion vertex
				V0Vtx   = 5, ///< Vertex from V0 decay
				KinkVtx = 6, ///< Kink vertex
				NotSpecified = -99 ///< Default value, no explicit type set
				*/
				ElementLink< xAOD::VertexContainer > vtxLink = trk->auxdata<ElementLink< xAOD::VertexContainer> >("vertexLink");
				int vertex_type=0;
				if(vtxLink.isValid()) vertex_type =  (*vtxLink)->vertexType();
				//tracking cuts
				//Only tracks associated with a jet
				float R = DeltaR(phi,eta,jet_phi,jet_eta);
				float jetRmin=999., jetPTmax=0.;
				unsigned  int jetMatchCount=0, jetRminIndex=-1, jetPTmaxIndex=-1;
				for(unsigned int ii=0; ii<jet_pt_xcalib_vector.size(); ii++){
					if(jet_pt_xcalib_vector.at(ii)<_pTjetCut || jet_isDummy) continue;
					float d=DeltaR(phi,eta,jet_phi_vector.at(ii),jet_eta_vector.at(ii));
					if(d < _dR_max)
					{
						++jetMatchCount;
						if(d<jetRmin){
							jetRmin=d;
							jetRminIndex=ii;
						}
						
						if(jet_pt_xcalib_vector.at(ii)>jetPTmax)
						{
							jetPTmax=jet_pt_xcalib_vector.at(ii);
							jetPTmaxIndex=ii;
						}
					}
				}
				
				if (jet_pt<110.  && jet_pt>80. && R <= _dR_max) {
				//if (i==0) {
					if(!_doSlimTree)
					{
						h_PixHits.at(cent_bin)->Fill(pt,eta,nPixHits);
						h_SCTHits.at(cent_bin)->Fill(pt,eta,nSCTHits);
						h_d0.at(cent_bin)->Fill(pt,d0);
						h_z0sintheta.at(cent_bin)->Fill(pt,z0pv);
					
						//all tracks
						h_reco_trk_map_nocuts->Fill(pt,eta,phi);
		            }
				}
				bool isTruthMatched=false;
				if(_data_switch==1){ //look for truth tracks in MC
					ElementLink< xAOD::TruthParticleContainer > truthLink = trk->auxdata<ElementLink< xAOD::TruthParticleContainer > >("truthParticleLink");				   
					if(truthLink.isValid() && trk->auxdata<float>("truthMatchProbability") > _mcProbCut) isTruthMatched = true;	   
				}
				if(!_doSlimTree) h_cut_flow[0]->Fill(pt,eta,(int)isTruthMatched+0.5);
				
				//Impact parameter cuts
				bool pased=true;
				//if(fabs(d0) > 2) pased=false; //(2,1.5,1)
				//else h_cut_flow[1]->Fill(pt,eta,(int)isTruthMatched+0.5);
				//if(fabs(z0pv) > 2) pased=false; //(2,1.5,1)
				//else h_cut_flow[2]->Fill(pt,eta,(int)isTruthMatched+0.5);
				
				if(!_doSlimTree) 
				{
					if(fabs(d0) < 2) h_cut_flow[1]->Fill(pt,eta,(int)isTruthMatched+0.5);
					if(fabs(z0pv) < 2) h_cut_flow[2]->Fill(pt,eta,(int)isTruthMatched+0.5);
				}
				//Silicon cuts
				if(nSCTHits + nPixHits < 7) pased=false;
				else if(!_doSlimTree) h_cut_flow[3]->Fill(pt,eta,(int)isTruthMatched+0.5);
				if(nPixHoles > 1) pased=false;
				else if(!_doSlimTree) h_cut_flow[4]->Fill(pt,eta,(int)isTruthMatched+0.5);
				if(nPixHoles + nSCTHoles > 2) pased=false;
				else if(!_doSlimTree) h_cut_flow[5]->Fill(pt,eta,(int)isTruthMatched+0.5);
				
				if (!pased) continue;
				
				//track spectra after loose cuts
				
				if(fabs(d0)<2 && fabs(z0pv)<2 && !_doSlimTree)
				{
					if(i==0) h_reco_trk_map->Fill(pt,eta,phi);
					//h_BL.at(cent_bin)->Fill(trk->vx(),trk->vy(),trk->vz());
					h_BL.at(cent_bin)->Fill(0.0,0.0,trk->vz());
					if((nShPixH + nShSCTH/2.)<=1) h_cut_flow[6]->Fill(pt,eta,(int)isTruthMatched+0.5);
				}
				
				//Filling the FF tree
				if(R <= _dR_max || jet_isDummy) {
					//get the tracks....
								
					track_eta.push_back(eta);
					track_phi.push_back(phi);
					track_pt.push_back(pt);
					track_vertex_type.push_back(vertex_type);
					track_isMuon.push_back(isMuon);
					if(jetMatchCount>=2)
					{
						if(i==jetRminIndex) track_multiJetMatch.push_back(1);
						else track_multiJetMatch.push_back(2);
						
						if(i==jetPTmaxIndex) track_multiJetMatch_B.push_back(1);
						else track_multiJetMatch_B.push_back(2);
					}
					else
					{
						track_multiJetMatch.push_back(0);
						track_multiJetMatch_B.push_back(0);
					}
					
					if(trk->charge()>0) track_charge.push_back(1);
					else track_charge.push_back(-1);
					
					for(int n=0;n<_nTrkSelTools;++n)
					{
						if(m_trackSelectionTool[n]->accept(*trk,(*vtx_itr))) track_passedTrkSelTool[n].push_back(1);
						else track_passedTrkSelTool[n].push_back(0);
					}
					
					track_nPixHits.push_back(nPixHits);
					track_nPixHoles.push_back(nPixHoles);
					//track_nPixDeadS.push_back(nPixDeadS);
					track_nShPixH.push_back(nShPixH);
					track_nSCTHits.push_back(nSCTHits);
					track_nSCTHoles.push_back(nSCTHoles);
					//track_nSCTDeadS.push_back(nSCTDeadS);
					track_nShSCTH.push_back(nShSCTH);
					track_nTRTHits.push_back(nTRTHits);
					track_nIBLHits.push_back(nIBLHits);
					track_expIBLHits.push_back(expIBLHits);
					track_nBLHits.push_back(nBLHits);
					track_expBLHits.push_back(expBLHits);

					track_d0.push_back(d0);
					track_z0sintheta.push_back(z0pv);
					track_ed0.push_back(ed0);
					track_ez0sintheta.push_back(ez0pv);
										//Get z			
					track_z.push_back(cos(R)*pt / jet_pt);
					//cout << "z: " << cos(R)*pt / jet_pt << " jet_pt: " << jet_pt << " trk pt: " << pt << endl; 
					
					if(_data_switch==1){ //look for truth tracks in MC
						ElementLink< xAOD::TruthParticleContainer > truthLink = trk->auxdata<ElementLink< xAOD::TruthParticleContainer > >("truthParticleLink");				   
						if(truthLink.isValid()){
						
							int trktype = getTypeReco((*truthLink)->barcode(),(*truthLink)->pdgId(),(*truthLink)->status(),(*truthLink)->charge(),trk->auxdata<float>("truthMatchProbability"),_mcProbCut);
						
							track_mc_pt.push_back((*truthLink)->pt()*0.001); //bring into GeV
							track_mc_eta.push_back((*truthLink)->eta());
							track_mc_phi.push_back((*truthLink)->phi());
							track_mc_probability.push_back(trk->auxdata<float>("truthMatchProbability"));
							track_mc_pdg.push_back((*truthLink)->pdgId());
							track_mc_barcode.push_back((*truthLink)->barcode());
							track_mc_type.push_back(trktype);
							track_mc_charge.push_back((*truthLink)->threeCharge()/3);
						}
						else
						{
							track_mc_pt.push_back(-999.);
							track_mc_eta.push_back(-999.);
							track_mc_phi.push_back(-999.);
							track_mc_probability.push_back(-999.);
							track_mc_pdg.push_back(-999);
							track_mc_barcode.push_back(-999);
							track_mc_type.push_back(0); // fake (no associated truth particle)
							track_mc_charge.push_back(0);
						}
					}
					
					if(!_doSlimTree)
					{
						if(fabs(d0)<2 && fabs(z0pv)<2)
						{
							if(_data_switch == 1) {
								hjetz_trig[0]->Fill(track_z.back(),jet_pt);
								if (jet_pt > 30.){
									hz->Fill(track_z.back());
								}	
							}	
							else{	
								for (int k=0;k<_nTriggers;k++){
									if(jet_isTriggered[k]) {
										hjetz_trig[k]->Fill(track_z.back(),jet_pt);
										hz->Fill(track_z.back());
									}	
								}	
							}
						}
					}
				}
				
			} // end reco track loop
			
			
			tree_ff->Fill();
			
			
		}// end reco jet loop 
	}
	

	//Loop over truth jets
	if(_data_switch == 1){
	
		const xAOD::TruthParticleContainer * particles = 0;
		if ( !event->retrieve( particles, "TruthParticles" ).isSuccess() ){
			Error("execute()", "Failed to retrieve TruthParticle container. Exiting." );
			return EL::StatusCode::FAILURE;
		}
		
		// add dummy jet for MB
		if(_isMB)
		{
			truth_jet_pt_vector.push_back(-1); 
			truth_jet_phi_vector.push_back(0);
			truth_jet_eta_vector.push_back(-99.);
			truth_jet_m_vector.push_back(-1);
			truth_jet_isDummy_vector.push_back(1); 
			truth_jet_NBJ_pT_vector.push_back(-1);
		}
		
		
		
		for(unsigned int i=0; i<truth_jet_pt_vector.size(); i++){
			truth_jet_pt = truth_jet_pt_vector.at(i);
			if (truth_jet_pt<_pTjetCut && !truth_jet_isDummy_vector.at(i)) continue;
			truth_jet_eta = truth_jet_eta_vector.at(i);
			truth_jet_phi = truth_jet_phi_vector.at(i);
			truth_jet_m = truth_jet_m_vector.at(i);
			if(!_doSlimTree) truth_jet_NBJ_pT = truth_jet_NBJ_pT_vector.at(i);
			truth_jet_isDummy = truth_jet_isDummy_vector.at(i);
						 
			//truth_jet_flavour = GetJetFlavour(truth_jet_phi,truth_jet_eta,mc_unstable_gen_pt,mc_unstable_gen_eta,mc_unstable_gen_phi,mc_unstable_pdg);
			float max_mu_pT = 0.; 
			float max_el_pT = 0.;
			
			xAOD::TruthParticleContainer::const_iterator truth_itr = particles->begin();
			xAOD::TruthParticleContainer::const_iterator truth_end = particles->end();
			
			truth_track_eta.clear();
			truth_track_phi.clear();
			truth_track_pt.clear();
			truth_track_z.clear();  
			truth_track_pdg.clear();
			truth_track_barcode.clear();
			truth_track_charge.clear();
			truth_track_multiJetMatch.clear();
			truth_track_multiJetMatch_B.clear();
			
			jet_centrality=cent_bin;
			
			for( ; truth_itr!=truth_end; ++truth_itr){
			
				int ty=getTypeTruth((*truth_itr)->barcode(),(*truth_itr)->pdgId(),(*truth_itr)->status(),(*truth_itr)->charge());
				
				if(ty!=1 && ty!=5) continue;
				
				//get the tracks....
				float pt = (*truth_itr)->pt()/ 1000.0;
				if (fabs(pt)<_pTtrkCut) continue;
				
				float eta = (*truth_itr)->eta();
				if (fabs(eta)>2.5) continue;
				
				float phi = (*truth_itr)->phi();
				
				if (i==0 && !_doSlimTree) h_cut_flow[7]->Fill(pt,eta,1.5);
				
				//Only tracks associated with a jet		
				float R = DeltaR(phi,eta,truth_jet_phi,truth_jet_eta);
				if(R > _dR_max && !truth_jet_isDummy) continue;
				
				float jetRmin=999., jetPTmax=0. ;
				unsigned int jetMatchCount=0, jetRminIndex=-1, jetPTmaxIndex=-1;
				
				for(unsigned int ii=0; ii<truth_jet_pt_vector.size(); ii++){
					if(truth_jet_pt_vector.at(ii)<_pTjetCut || truth_jet_isDummy) continue;
					float d=DeltaR(phi,eta,truth_jet_phi_vector.at(ii),truth_jet_eta_vector.at(ii));
					if(d < _dR_max)
					{
						++jetMatchCount;
						if(d<jetRmin){
							jetRmin=d;
							jetRminIndex=ii;
						}
						if(truth_jet_pt_vector.at(ii)>jetPTmax)
						{
							jetPTmax=truth_jet_pt_vector.at(ii);
							jetPTmaxIndex=ii;
						}
					}
				}
								
				
				//get the tracks....
				truth_track_eta.push_back(eta);
				truth_track_phi.push_back(phi);
				truth_track_pt.push_back(pt);
				truth_track_pdg.push_back((*truth_itr)->pdgId());
				truth_track_barcode.push_back((*truth_itr)->barcode());
				truth_track_charge.push_back((*truth_itr)->threeCharge()/3);
				if(jetMatchCount>=2)
				{
					if(i==jetRminIndex) truth_track_multiJetMatch.push_back(1);
					else truth_track_multiJetMatch.push_back(2);
					
					if(i==jetPTmaxIndex) truth_track_multiJetMatch_B.push_back(1);
					else truth_track_multiJetMatch_B.push_back(2);
				}
				else
				{
					truth_track_multiJetMatch.push_back(0);
					truth_track_multiJetMatch_B.push_back(0);
				}
				
				//Get z			
				truth_track_z.push_back(cos(R)*pt / truth_jet_pt);		
				if(!truth_jet_isDummy && !_doSlimTree) truth_hjetz->Fill(truth_track_z.back(),truth_jet_pt);
				
				//test heavy flavours
				if (fabs((*truth_itr)->pdgId())==13 && pt > max_mu_pT) max_mu_pT = pt; 
 	         	//test electrons  
 	 	    	if (fabs((*truth_itr)->pdgId())==11 && pt > max_el_pT) max_el_pT = pt;
			}
			truth_jet_Muon_pT=max_mu_pT;
			truth_jet_electron_pT=max_el_pT;
			
			truth_tree_ff->Fill();
			
		} // end loop over truth jets
	}

	//delete shallow copy of muon container 
	delete muons_shallowCopy.first; 
	delete muons_shallowCopy.second; 
	
	delete electrons_shallowCopy.first; 
	delete electrons_shallowCopy.second;
	
	return EL::StatusCode::SUCCESS;
}

EL::StatusCode pPbFragmentation :: postExecute (){
	return EL::StatusCode::SUCCESS;
}

EL::StatusCode pPbFragmentation :: finalize (){
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
	
	//cleaning muons
	if(m_muonCalibrationAndSmearingTool){ 
		delete m_muonCalibrationAndSmearingTool; 
		m_muonCalibrationAndSmearingTool = 0; 
	 } 
	 if(m_muonSelection){ 
		delete m_muonSelection; 
		m_muonSelection = 0; 
	 }
	 
	 //cleaning electrons
	 if(m_EgammaCalibrationAndSmearingTool)
	 {
		delete m_EgammaCalibrationAndSmearingTool;
		m_EgammaCalibrationAndSmearingTool=0;
	 }
	 if(m_LHToolTight2015)
	 {
		delete m_LHToolTight2015;
		m_LHToolTight2015=0;
	 }
	 
	 //cleaning tracks
	 for(int n=0;n<_nTrkSelTools;++n)
	 {
		if(m_trackSelectionTool[n])
		{
			delete m_trackSelectionTool[n];
			m_trackSelectionTool[n]=0;
		}
	 }
	
	//cleaning jets
	if(jesProv)
	{
		delete jesProv;
		jesProv=0;
	}
	if(jerTool)
	{
		delete jerTool;
		jerTool=0;
	}
	if(smearTool)
	{
		delete smearTool;
		smearTool=0;
	}
	
	return EL::StatusCode::SUCCESS;
}
	
EL::StatusCode pPbFragmentation :: histFinalize ()
{  
	cout<<"Events = "<< m_eventCounter<<endl;
	return EL::StatusCode::SUCCESS;
}
