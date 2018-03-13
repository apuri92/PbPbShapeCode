#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include "pPbFragmentation/MBUEEstimator.h"
#include "pPbCentrality/pPbMinBiasUtil.h"
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/tools/Message.h"
#include "xAODJet/JetContainer.h"
#include "xAODTrigger/JetRoIContainer.h"
#include "xAODTrigger/JetRoIAuxContainer.h"
#include <TFile.h>
#include <TSystem.h>

#include "xAODHIEvent/HIEventShapeContainer.h"
#include "xAODForward/ZdcModuleContainer.h"

using namespace std;
using namespace JetHelperTools;
using namespace TrackHelperTools;
using namespace MTCorrector;

ClassImp(MBUEEstimator)

#define EL_RETURN_CHECK( CONTEXT, EXP )			\
do {							\
if( ! EXP.isSuccess() ) {				\
Error( CONTEXT,					\
XAOD_MESSAGE( "Failed to execute: %s" ),	\
#EXP );					\
return EL::StatusCode::FAILURE;			\
}							\
} while( false )


MBUEEstimator :: MBUEEstimator ()
{
}

EL::StatusCode MBUEEstimator :: fileExecute ()
{
	// Here you do everything that needs to be done exactly once for every
	// single file, e.g. collect a list of all lumi-blocks processed
	return EL::StatusCode::SUCCESS;
}

EL::StatusCode MBUEEstimator :: changeInput (bool firstFile)
{
	return EL::StatusCode::SUCCESS;
}

EL::StatusCode MBUEEstimator :: setupJob (EL::Job& job)
{
	// let's initialize the algorithm to use the xAODRootAccess package
	job.useXAOD ();

	EL_RETURN_CHECK( "setupJob()", xAOD::Init() ); // call before opening first file
	cout << " Job setup done!" << endl;    
	
	return EL::StatusCode::SUCCESS;
}

EL::StatusCode MBUEEstimator :: initialize ()
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
		m_trigDecisionTool->msg().setLevel( MSG::ERROR);
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
	
	
	// GRL
	TString xfn = gSystem->GetFromPipe("echo $ROOTCOREBIN");
	TString xmlfile = xfn + "/../pPbFragmentation/data/"+ _GRL;

	m_grl = new GoodRunsListSelectionTool("GoodRunsListSelectionTool");
	std::vector<std::string> vecStringGRL;
	vecStringGRL.push_back(xmlfile.Data());
	EL_RETURN_CHECK("initialize()",m_grl->setProperty( "GoodRunsListVec", vecStringGRL));
	EL_RETURN_CHECK("initialize()",m_grl->setProperty("PassThrough", false)); // if true (default) will ignore result of GRL and will just pass all events
	EL_RETURN_CHECK("initialize()",m_grl->initialize());

	//d0 function
	f_d0_cut = new TF1("f1", "[0]*exp([1]*x)+[2]*exp([3]*x)", 0.4, 500);
	f_d0_cut->SetParameters(0.472367, -0.149934, 0.193095, 0.000337765);
		

	//Uncert tool
	uncertprovider = new UncertProvider(_uncert_index,_mcProbCut,_cut_level.c_str(), GetCentralityNBins(31)-1, _eff_jety);
	_mcProbCut = uncertprovider->GetMCProb();
	cout << "mc prob: " << _mcProbCut <<endl;
	
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

EL::StatusCode MBUEEstimator :: histInitialize ()
{
	
	FCal_only = false;
	cout << " Setting  histograms" << endl;
	
	int ptTrkBinsN, etaTrkBinsN, phiTrkBinsN;
	double ptTrkBins[1000],  etaTrkBins[1000], phiTrkBins[1000];

	SetupBinning(0, "pt-trk-shape", ptTrkBins, ptTrkBinsN);
	SetupBinning(0, "eta-trk-coars", etaTrkBins, etaTrkBinsN);
	SetupBinning(0, "phi-trk-coars", phiTrkBins, phiTrkBinsN);
	
	//Track corrector
	int _nCentbins = GetCentralityNBins(_centrality_scheme);
	trkcorr = new TrackCorrector(_cut_level.c_str(),GetCentralityNBins(31)-1,_eff_jety);
	for (int i = 0; i < ptTrkBinsN-1; i++){ 
		if (ptTrkBins[i]<_pTtrkCut && _pTtrkCut < ptTrkBins[i+1]) {trkcorr->trkpTThreshold = ptTrkBins[i]; break;}
	}	
	
	//Jet corrector
	jetcorr = new JetCorrector();
	
	//Basic histograms

	h_FCal_Et = new TH1D("h_FCal_Et",";FCal E_{T};N",1200,0.,6.);
	h_FCal_Et->Sumw2();
	h_FCal_Et_HP = new TH1D("h_FCal_Et_HP",";FCal E_{T};N",1200,0.,6.);
	h_FCal_Et_HP->Sumw2();
	h_FCal_Et_MC = new TH1D("h_FCal_Et_MC",";FCal E_{T};N",1200,0.,6.);
	h_FCal_Et_MC->Sumw2();
	
	h_dPsi = new TH1D("h_dPsi","h_dPsi",32,-TMath::Pi(),TMath::Pi());
	h_Psi = new TH1D("h_Psi","h_Psi",32,-TMath::Pi(),TMath::Pi());
	h_jet_phi = new TH1D("h_jet_phi","h_jet_hi",32,-TMath::Pi(),TMath::Pi());
	h_jet_phi_v_Psi = new TH2D("h_jet_phi_v_Psi","h_jet_v_Psi",32,-TMath::Pi(),TMath::Pi(),32,-TMath::Pi(),TMath::Pi());
	h_jet_phi_v_dPsi = new TH2D("h_jet_phi_v_dPsi","h_jet_v_dPsi",32,-TMath::Pi(),TMath::Pi(),32,-TMath::Pi(),TMath::Pi());
	
	wk()->addOutput (h_dPsi);
	wk()->addOutput (h_Psi);
	wk()->addOutput (h_jet_phi);
	wk()->addOutput (h_jet_phi_v_Psi);
	wk()->addOutput (h_jet_phi_v_dPsi);
	
	h_RejectionHisto = new TH1D("RejectionHisto","RejectionHisto",9,0,9);
	SetRejectionHistogram(h_RejectionHisto);
	
	h_centrality_HP = new TH1D("Centrality_HP","Centrality",10,0,10);
	h_centrality_HP->Sumw2();
	h_centrality_MC = new TH1D("Centrality_MC","Centrality",10,0,10);
	h_centrality_MC->Sumw2();
	
	h_triggercounter = new TH2D("h_triggercounter","h_triggercounter",_nTriggers,0,_nTriggers,2,-0.5,1.5);
	SetTrigger_hist(h_triggercounter);
	
	wk()->addOutput (h_FCal_Et);
	wk()->addOutput (h_FCal_Et_HP);
	wk()->addOutput (h_FCal_Et_MC);
	wk()->addOutput (h_RejectionHisto);
	wk()->addOutput (h_triggercounter);
	wk()->addOutput (h_centrality_HP);
	wk()->addOutput (h_centrality_MC);
	
	if (FCal_only) return EL::StatusCode::SUCCESS;
	
	TH3D* temphist_3D = nullptr;
	TH1D* temphist_1D = nullptr;
	
	int ndPsibins = 16;
	int ndRbins = 13;
	//h_trk_dNdEtadPhidpT =  vector<vector<TH3D*> > (ndPsibins, vector<TH3D*>(_nCentbins));
	h_UE_dNdEtadPhidpT_HP =  vector<vector<vector<TH3D*>>> (ndPsibins, vector<vector<TH3D*>> (_nCentbins,vector<TH3D*>(ndRbins)));
	h_jet_v_Psi_HP =  vector<TH1D*>(_nCentbins);
	h_UE_dNdEtadPhidpT_MC =  vector<vector<vector<TH3D*>>> (ndPsibins, vector<vector<TH3D*>> (_nCentbins,vector<TH3D*>(ndRbins)));
	h_jet_v_Psi_MC =  vector<TH1D*>(_nCentbins);
	for (int j=0;j<_nCentbins;j++){
		temphist_1D = new TH1D(Form("h_jet_v_Psi_HP_cent%i",j),Form("h_jet_v_Psi_cent%i",j),16,0,16);
		h_jet_v_Psi_HP.at(j) = temphist_1D;
		h_jet_v_Psi_HP.at(j)->Sumw2();
		wk()->addOutput (h_jet_v_Psi_HP.at(j));
		temphist_1D = new TH1D(Form("h_jet_v_Psi_MC_cent%i",j),Form("h_jet_v_Psi_cent%i",j),16,0,16);
		h_jet_v_Psi_MC.at(j) = temphist_1D;
		h_jet_v_Psi_MC.at(j)->Sumw2();
		wk()->addOutput (h_jet_v_Psi_MC.at(j));
		for (int i=0;i<ndPsibins;i++){
			//cout << "  Bins " << ptTrkBinsN << " " << etaTrkBinsN << " " << phiTrkBinsN << endl;
			//temphist_3D = new TH3D(Form("h_trk_dNdEtadPhidpT_dPsi%i_cent%i",i,j),Form("h_trk_dNdEtadPhidpT_dPsi%i_cent%i",i,j),ptTrkBinsN, ptTrkBins,etaTrkBinsN,etaTrkBins,phiTrkBinsN,phiTrkBins);
			//h_trk_dNdEtadPhidpT.at(i).at(j) = temphist_3D;			        
            //h_trk_dNdEtadPhidpT.at(i).at(j)->Sumw2();
            //wk()->addOutput (h_trk_dNdEtadPhidpT.at(i).at(j));
            
            for (int k=0;k<ndRbins;k++){
		        temphist_3D = new TH3D(Form("h_UE_dNdEtadPhidpT_HP_dPsi%i_cent%i_dR%i",i,j,k),Form("h_UE_dNdEtadPhidpT_dPsi%i_cent%i_dR%i",i,j,k),ptTrkBinsN, ptTrkBins,etaTrkBinsN,etaTrkBins,phiTrkBinsN,phiTrkBins);
				h_UE_dNdEtadPhidpT_HP.at(i).at(j).at(k) = temphist_3D;			        
		        h_UE_dNdEtadPhidpT_HP.at(i).at(j).at(k)->Sumw2();
		        wk()->addOutput (h_UE_dNdEtadPhidpT_HP.at(i).at(j).at(k));
		        
		        temphist_3D = new TH3D(Form("h_UE_dNdEtadPhidpT_MC_dPsi%i_cent%i_dR%i",i,j,k),Form("h_UE_dNdEtadPhidpT_dPsi%i_cent%i_dR%i",i,j,k),ptTrkBinsN, ptTrkBins,etaTrkBinsN,etaTrkBins,phiTrkBinsN,phiTrkBins);
				h_UE_dNdEtadPhidpT_MC.at(i).at(j).at(k) = temphist_3D;			        
		        h_UE_dNdEtadPhidpT_MC.at(i).at(j).at(k)->Sumw2();
		        wk()->addOutput (h_UE_dNdEtadPhidpT_MC.at(i).at(j).at(k));
		    }    
		}	
	}
	cout << " Histograms ready" << endl;
	return EL::StatusCode::SUCCESS;
}

//Loop over events
EL::StatusCode MBUEEstimator :: execute (){

	xAOD::TEvent* event = wk()->xaodEvent();
	
	//cout << "test" << endl;
	
	float dRbins[14]={0. , 0.05 , 0.1 , 0.15 , 0.2 , 0.25 , 0.3 , 0.4 , 0.5 , 0.6 , 0.7 , 0.8 , 1.0 , 1.2 };
	// Event counter
	int statSize=1;
	
	if(m_eventCounter!=0)
	{
		double power=std::floor(log10(m_eventCounter));
		statSize=(int)std::pow(10.,power);
	}
	if(m_eventCounter%statSize==0) std::cout << "Event: " << m_eventCounter << std::endl;  
	m_eventCounter++;
	//if (m_eventCounter < 1916) return EL::StatusCode::SUCCESS; 
 	//std::cout << "Event: " << m_eventCounter << std::endl;

	//All events
	bool keep = true;
	h_RejectionHisto->Fill(0.5);

	//---------------------------
	//     Event information 
	//---------------------------

	const xAOD::EventInfo* eventInfo = 0;
	EL_RETURN_CHECK("execute",event->retrieve( eventInfo, "EventInfo"));

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

	//Get centrality bin and centile. Centile used for MB weighting (from MB_FCal_Normalization.txt)
	FCalEt = 0;
	int cent_bin = 0;
	int cent_bin_corse = 0;
	double event_weight_fcal_2 = 1.;
	double event_weight_fcal_3 = 1.;
	float Psi=-999;
	int dPsi_bin;
	//Centrality
	const xAOD::HIEventShapeContainer* calos=0;
	EL_RETURN_CHECK("execute",event->retrieve( calos, "CaloSums"));
	if (_centrality_scheme>1)
	{	
		FCalEt=calos->at(5)->et()*1e-6;
		cent_bin = GetCentralityBin(_centrality_scheme, FCalEt,  isHIJING );
		cent_bin_corse = GetCentralityBin(31, FCalEt,  isHIJING ); //Need for some tools
		event_weight_fcal_2 = jetcorr->GetFCalWeight(FCalEt,2);
		event_weight_fcal_3 = jetcorr->GetFCalWeight(FCalEt,3);
		Psi = GetEventPlane(calos);
	}

	if (cent_bin < 0) {
		if(cent_bin==-2) Error("execute()", "Unknown centrality scheme" );
		h_RejectionHisto->Fill(1.5);
		keep = false;
	}

	// GRL
	if(!isMC)
	{
		if(!m_grl->passRunLB(*eventInfo))
		{
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


	// find primary vertex
	xAOD::VertexContainer::const_iterator vtx_itr = vertices->begin();
	xAOD::VertexContainer::const_iterator vtx_end = vertices->end();
	const xAOD::Vertex* primaryVertex = 0;
	for(;vtx_itr!=vtx_end;++vtx_itr)
	{
		if((*vtx_itr)->vertexType()==xAOD::VxType::PriVtx) {
			primaryVertex = (*vtx_itr);
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

	//DAQ errors
	if(!isMC){
		if(   (eventInfo->errorState(xAOD::EventInfo::LAr)==xAOD::EventInfo::Error ) || (eventInfo->errorState(xAOD::EventInfo::Tile)==xAOD::EventInfo::Error ) || (eventInfo->errorState(xAOD::EventInfo::SCT)==xAOD::EventInfo::Error ) || (eventInfo->isEventFlagBitSet(xAOD::EventInfo::Core, 18) ) ){
			h_RejectionHisto->Fill(4.5);
			keep = false;
		}
	}

	//Pileup	
	if (_doPileupRejection){		
		bool m_is_pileup = false;
		if (!isMC) {	
			const xAOD::ZdcModuleContainer* zdcMod = 0;
			EL_RETURN_CHECK("execute",event->retrieve( zdcMod, "ZdcModules")); 
			// ZDC
			m_zdcTools->reprocessZdc();
			// is Pileup
			m_is_pileup = m_hiPileup->is_pileup( *calos, *zdcMod); // SAVE pileup Decision HERE 0 = NO pileup, 1 = pileup
		}
		else m_is_pileup = (FCalEt > 4.8); //Remove pileup in MC
				
		if (m_is_pileup){
			h_RejectionHisto->Fill(6.5);
			keep = false;
		}
	}
	

	if (!keep) return EL::StatusCode::SUCCESS; // go to the next event

	h_RejectionHisto->Fill(7.5);

	// trigger
	double trigger_prescale = 1.;
	if (_data_switch==0)
	{
		int event_passed_trigger=0;

		for (int i=0;i<_nTriggers;i++){

			event_isTriggered[i] = false;

			event_isTriggered[i] =  _chainGroup.at(i)->isPassed();
			h_triggercounter->Fill(i, (Double_t) event_isTriggered[i]);
			if(event_isTriggered[i]) {
				event_passed_trigger=1;
				trigger_prescale=trigger_PS.at(i); 
			}	
		}

		if(!event_passed_trigger) return EL::StatusCode::SUCCESS; // go to next event
		else h_RejectionHisto->Fill(8.5);
	}

	h_FCal_Et_HP->Fill(FCalEt, event_weight_fcal_2*trigger_prescale);
	h_FCal_Et_MC->Fill(FCalEt, event_weight_fcal_3*trigger_prescale);
	h_FCal_Et->Fill(FCalEt,trigger_prescale);
	h_centrality_HP->Fill(cent_bin,event_weight_fcal_2*trigger_prescale);
	h_centrality_MC->Fill(cent_bin,event_weight_fcal_3*trigger_prescale);
	
	if (FCal_only) return EL::StatusCode::SUCCESS;

	//Tracks
	const xAOD::TrackParticleContainer* recoTracks = 0;
	EL_RETURN_CHECK("execute",event->retrieve( recoTracks, "InDetTrackParticles"));
		
	float event_weight = 1;
	double max_pt = 1;
	
	//cout << "Psi " << Psi << endl;
	
	h_Psi->Fill(Psi);
	
	for (int dRbin = 0; dRbin<13;dRbin++){		
		for (int phibin = 1; phibin <=h_UE_dNdEtadPhidpT_HP.at(0).at(0).at(0)->GetZaxis()->GetNbins();phibin++){	
			float jet_phi = h_UE_dNdEtadPhidpT_HP.at(0).at(0).at(0)->GetZaxis()->GetBinCenter(phibin);
			float jet_dPsi_bin = GetPsiBin(DeltaPsi(jet_phi,Psi));
			h_dPsi->Fill(DeltaPsi(jet_phi,Psi));
			h_jet_phi->Fill(jet_phi);
			h_jet_phi_v_Psi->Fill(jet_phi,Psi);
			h_jet_phi_v_dPsi->Fill(jet_phi,DeltaPsi(jet_phi,Psi));
			for (int etabin = 1; etabin <=h_UE_dNdEtadPhidpT_HP.at(0).at(0).at(0)->GetYaxis()->GetNbins() ;etabin++){
				float jet_eta = h_UE_dNdEtadPhidpT_HP.at(0).at(0).at(0)->GetYaxis()->GetBinCenter(etabin);
				if (fabs(jet_eta)>1.3) continue;
				for (const auto& trk : *recoTracks) {
					float eta = trk->eta();
					float phi = trk->phi();
					float dR = DeltaR(jet_phi, jet_eta, phi, eta);
					if (dR>=dRbins[dRbin+1] || dR<dRbins[dRbin]) continue;
					float pt = trk->pt()/1000.;
					//if (fabs(jet_eta-eta)>dRbins[dRbin+1]) continue;
						
					//correct the alignement
					bool isMatchedToTruthParticle = false;
					if (isMC) {
						ElementLink< xAOD::TruthParticleContainer > truthLink = trk->auxdata<ElementLink< xAOD::TruthParticleContainer > >("truthParticleLink");				   
						float mcprob =trk->auxdata<float>("truthMatchProbability");
						if(truthLink.isValid() && mcprob > _mcProbCut) isMatchedToTruthParticle = true;
					}
					if (_correctTrackpT && (!isMC || !isMatchedToTruthParticle)) trkcorr->correctChTrackpT(pt, eta, phi, trk->charge());
					if (fabs(eta) > 2.5) continue;
					//Charge cut
					if (_useCharge!=0 && ((int)trk->charge())!=_useCharge) continue; 
					if (isMC && uncertprovider->uncert_class==5) uncertprovider->UncerTrackMomentum(pt, eta, phi, trk->charge() );
					if (pt < _pTtrkCut) continue; //min pT cut
					double d0 = trk->d0();
					double d0_cut = f_d0_cut->Eval(pt);
					if(fabs(d0) > d0_cut) continue; //pT dependant d0 cut
					if(!m_trackSelectorTool->accept(*trk)) continue; //track selector tool
			
					
					float eff_uncertainty = 0;
					if (_uncert_index > 0 && uncertprovider->uncert_class==4) eff_uncertainty = uncertprovider->CorrectTrackEff(pt,eta, 1.0, cent_bin_corse); 
					//float eff_weight = trkcorr->get_effcorr(pt, eta, cent_bin_corse, eff_uncertainty, 10, 0., 1.0); 
					float eff_weight = trkcorr->get_effcorr(pt, eta, cent_bin_corse, 0, _dataset);
			
					//h_trk_dNdEtadPhidpT.at(dPsi_bin).at(cent_bin)->Fill(pt,eta,phi,event_weight_fcal*eff_weight);
					//cout << " jet_phi " << jet_phi << " Psi " << Psi << "dPsi " << DeltaPsi(jet_phi,Psi) << " bin " <<  jet_dPsi_bin << endl;
					h_UE_dNdEtadPhidpT_HP.at(jet_dPsi_bin).at(cent_bin).at(dRbin)->Fill(pt,jet_eta,jet_phi,event_weight_fcal_2*eff_weight*trigger_prescale);
					h_UE_dNdEtadPhidpT_MC.at(jet_dPsi_bin).at(cent_bin).at(dRbin)->Fill(pt,jet_eta,jet_phi,event_weight_fcal_3*eff_weight*trigger_prescale);
				}
			}
		}
			
	}
	
	//Get normalization
	for (int phibin = 1; phibin <=h_UE_dNdEtadPhidpT_HP.at(0).at(0).at(0)->GetZaxis()->GetNbins();phibin++){
		float jet_phi = h_UE_dNdEtadPhidpT_HP.at(0).at(0).at(0)->GetZaxis()->GetBinCenter(phibin);
		float jet_dPsi_bin = GetPsiBin(DeltaPsi(jet_phi,Psi));
		h_jet_v_Psi_HP.at(cent_bin)->Fill(jet_dPsi_bin,event_weight_fcal_2/h_UE_dNdEtadPhidpT_HP.at(0).at(0).at(0)->GetZaxis()->GetNbins()*trigger_prescale);
		h_jet_v_Psi_MC.at(cent_bin)->Fill(jet_dPsi_bin,event_weight_fcal_3/h_UE_dNdEtadPhidpT_MC.at(0).at(0).at(0)->GetZaxis()->GetNbins()*trigger_prescale);
	}

		
	return EL::StatusCode::SUCCESS;
}

EL::StatusCode MBUEEstimator :: postExecute (){
	return EL::StatusCode::SUCCESS;
}

EL::StatusCode MBUEEstimator :: finalize (){
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

	if(m_trackSelectorTool)
	{
		delete m_trackSelectorTool;
		m_trackSelectorTool=0;
	}

	return EL::StatusCode::SUCCESS;
}

EL::StatusCode MBUEEstimator :: histFinalize (){
	cout<<"Events = "<< m_eventCounter<<endl;
	return EL::StatusCode::SUCCESS;
}
