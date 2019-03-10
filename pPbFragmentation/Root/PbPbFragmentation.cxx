#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include "pPbFragmentation/PbPbFragmentation.h"
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

#include "xAODHIEvent/HIEventShapeContainer.h"
#include "xAODForward/ZdcModuleContainer.h"

using namespace std;
using namespace JetHelperTools;
using namespace TrackHelperTools;
using namespace MTCorrector;

ClassImp(PbPbFragmentation)

#define EL_RETURN_CHECK( CONTEXT, EXP )			\
do {							\
if( ! EXP.isSuccess() ) {				\
Error( CONTEXT,					\
XAOD_MESSAGE( "Failed to execute: %s" ),	\
#EXP );					\
return EL::StatusCode::FAILURE;			\
}							\
} while( false )


PbPbFragmentation :: PbPbFragmentation ()
{
}

EL::StatusCode PbPbFragmentation :: fileExecute ()
{
	// Here you do everything that needs to be done exactly once for every
	// single file, e.g. collect a list of all lumi-blocks processed
	return EL::StatusCode::SUCCESS;
}

EL::StatusCode PbPbFragmentation :: changeInput (bool firstFile)
{
	return EL::StatusCode::SUCCESS;
}

//Loop over events
EL::StatusCode PbPbFragmentation :: execute (){

	xAOD::TEvent* event = wk()->xaodEvent();
	
	// Event counter
	int statSize=1;
	//For testing
	
	if(m_eventCounter!=0)
	{
		double power=std::floor(log10(m_eventCounter));
		statSize=(int)std::pow(10.,power);
	}
	if(m_eventCounter%statSize==0) std::cout << "Event: " << m_eventCounter << std::endl;
	m_eventCounter++;
	//if (m_eventCounter < 7000) return EL::StatusCode::SUCCESS;
	//std::cout << "Event: " << m_eventCounter << std::endl;

	//All events
	bool keep = true;
	h_RejectionHisto->Fill(0.5);

	//---------------------------
	//     Event information - TODO write an event selector
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
	int cent_bin_fine = 0;
	int cent_bin_corse = 0;
	double event_weight_fcal = 1;
	//Centrality
	const xAOD::HIEventShapeContainer* calos=0;
	EL_RETURN_CHECK("execute",event->retrieve( calos, "CaloSums"));
	if (_centrality_scheme>1)
	{	
		FCalEt=calos->at(5)->et()*1e-6;
		cent_bin = GetCentralityBin(_centrality_scheme, FCalEt,  isHIJING );
		cent_bin_fine = GetCentralityBin(30, FCalEt,  isHIJING ); //Need for some tools
		cent_bin_corse = GetCentralityBin(31, FCalEt,  isHIJING ); //Need for some tools
		if (isMC) event_weight_fcal = jetcorr->GetFCalWeight(FCalEt,1);
		h_centrality->Fill(cent_bin,event_weight_fcal);
		
		//Get HI clusters for flow
		//const xAOD::CaloClusterContainer *hiclus(0);
		//EL_RETURN_CHECK("execute",event->retrieve(hiclus,"HIClusters") );
		//cout << "Psi 1: " << GetEventPlane(hiclus) << " Psi 2: " << GetEventPlane(calos) << endl;
		uee->Psi = GetEventPlane(calos);
		uee->Psi3 = GetEventPlane(calos,3); 
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
			//No vx and vy info so removing line below
//			if(cent_bin>=0) h_vx.at(cent_bin)->Fill(primaryVertex->x(),primaryVertex->y(),primaryVertex->z());
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
			//const xAOD::HIEventShapeContainer* hiev = 0;
			//EL_RETURN_CHECK("execute",event->retrieve( hiev, "HIEventShape"));
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

	//TODO MBTS timing cut?

	h_RejectionHisto->Fill(7.5);

	// trigger
	if (_data_switch==0)
	{
		int event_passed_trigger=0;

		for (int i=0;i<_nTriggers;i++){

			event_isTriggered[i] = false;

			event_isTriggered[i] =  _chainGroup.at(i)->isPassed();
			h_triggercounter->Fill(i, (Double_t) event_isTriggered[i]);
			if(event_isTriggered[i]) event_passed_trigger=1;
		}

		if(!event_passed_trigger) return EL::StatusCode::SUCCESS; // go to next event
		else h_RejectionHisto->Fill(8.5);
	}

	h_FCal_Et->Fill(FCalEt, event_weight_fcal); //filled here to get proper event weight

	//Tracks
	const xAOD::TrackParticleContainer* recoTracks = 0;
	EL_RETURN_CHECK("execute",event->retrieve( recoTracks, "InDetTrackParticles"));
		
	//TODO	electrons
//	const xAOD::ElectronContainer* electrons = 0;
//	EL_RETURN_CHECK("execute()",event->retrieve( electrons, "Electrons" ));
//
//	std::pair< xAOD::ElectronContainer*, xAOD::ShallowAuxContainer* > electrons_shallowCopy;
//	electrons_shallowCopy = xAOD::shallowCopyContainer( *electrons );

	//Jet vectors
	vector<float> jet_pt_xcalib_vector,jet_phi_vector,jet_eta_vector, jet_y_vector, jet_TrigPresc_vector;
	vector<bool> jet_isolated_vector, jet_IsTrig_vector, isFake_vector;
	vector<float> truth_jet_eta_vector,truth_jet_phi_vector,truth_jet_pt_vector, truth_jet_y_vector;
	vector<bool> truth_jet_isolated_vector;
	vector<int> truth_jet_indices;
	vector<int> TruthJetIndex;
	vector<int> isTriggered[_nTriggers];
	vector<float> antikt2_pt,antikt2_phi,antikt2_eta;

	vector<float> jet_uJER_vector;
	vector<vector<float>> jet_uJES_vector;

	float event_weight = 1;
	double max_pt = 1;

	TLorentzVector jet4vector;

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

			double pt    = (jet_truth_4mom.pt() * 0.001 );
			double eta    = (jet_truth_4mom.eta());
			double phi    = (jet_truth_4mom.phi());
			double m   = (jet_truth_4mom.M()*0.001);

			if (pt>max_pt) 	//event weight from leading truth jet
			{
				event_weight = jetcorr->GetJetWeight(pt, eta, phi);
				max_pt = pt;
			}

			//if (pt < _truthpTjetCut) continue;
			//if (fabs(eta)>2.1) continue;

			//filling truth pt/eta/phi vectors
			truth_jet_pt_vector.push_back(pt);
			truth_jet_phi_vector.push_back(phi);
			truth_jet_eta_vector.push_back(eta);
			////Get rapidity
			jet4vector.SetPtEtaPhiM(pt, eta, phi, m);
			truth_jet_y_vector.push_back(jet4vector.Rapidity());
		}
		event_weight = event_weight*event_weight_fcal; //event weight is only set if MC. Otherwise default is 1.

		//getting truth isolation vector (vector<bool>)
		truth_jet_isolated_vector = MTCorrector::GetIsolation(truth_jet_pt_vector,truth_jet_eta_vector,truth_jet_phi_vector, 1.0, _pt_iso); // -1 is for pT_iso == jet pT
	}

	// ---- GETTING RECO JETS ----
	xAOD::TStore *store = new xAOD::TStore; //For calibration
	const xAOD::JetContainer* jets = 0;
	EL_RETURN_CHECK("execute()",event->retrieve( jets, _reco_jet_collection.c_str() ));
	
	xAOD::JetContainer* updatedjets = new xAOD::JetContainer();
	xAOD::AuxContainerBase* updatedjetsAux = new xAOD::AuxContainerBase();
	updatedjets->setStore( updatedjetsAux );
	
	store->record(updatedjets,"updatedjets");
	store->record(updatedjetsAux,"updatedjetsAux");

	xAOD::JetContainer::const_iterator jet_itr = jets->begin();
	xAOD::JetContainer::const_iterator jet_end = jets->end();
	
	const xAOD::JetContainer * track_jets = 0;
	if (_doFJR) {	
		EL_RETURN_CHECK("execute()",event->retrieve( track_jets, "AntiKt4HITrackJets" ));
	}
	
	for( ; jet_itr != jet_end; ++jet_itr )
	{
		xAOD::Jet newjet;// = new xAOD::Jet();
		newjet.makePrivateStore( **jet_itr );

		const xAOD::JetFourMom_t jet_4mom_def = newjet.jetP4();
		float def_jet_pt  = (jet_4mom_def.pt() * 0.001);
		
		//cout << " Def:  " << def_jet_pt << endl;

		xAOD::JetFourMom_t jet_4mom = newjet.jetP4("JetSubtractedScaleMomentum"); //getting SubtractedScale instead of EMScale because EMScale is not in DFAntiKt4HI
		float uncalib_jet_pt  = (jet_4mom.pt() * 0.001);

		const xAOD::JetFourMom_t jet_4mom_unsubtracted = newjet.jetP4("JetUnsubtractedScaleMomentum");
		float unsubtracted_jet_pt  = (jet_4mom_unsubtracted.pt() * 0.001);

		hET_ETsub->Fill(def_jet_pt,jet_4mom_def.eta(),unsubtracted_jet_pt-uncalib_jet_pt);

		//newjet.setJetP4("JetPileupScaleMomentum",jet_4mom); //Setting PileupScale and ConstitScale because they are not in DFAntiKt4HI
		newjet.setJetP4("JetConstitScaleMomentum",jet_4mom_unsubtracted); //Required
		//EL_RETURN_CHECK("execute()", m_jetCalibration->applyCalibration( newjet ) );
		
		//cout << " Calib " << (newjet.pt() * 0.001);
		
		const xAOD::JetFourMom_t jet_4mom_xcalib = newjet.jetP4();
		newjet.setJetP4("JetGSCScaleMomentum", jet_4mom_xcalib);

		//Cross-calibration
		if (_data_switch==0) EL_RETURN_CHECK("execute()", m_jetCalibration_insitu->applyCalibration( newjet ) );
		
		//cout << " Cross-Calib " << (newjet.pt() * 0.001) << endl;
		
		//Uncertainties
		if (isMC && (_uncert_index > 0) && ((uncertprovider->uncert_class<4) || (uncertprovider->uncert_class==6)) ) {
			xAOD::JetContainer::const_iterator truth_jet_itr = jet_truth->begin();
			xAOD::JetContainer::const_iterator truth_jet_end = jet_truth->end();
			xAOD::Jet *truthMjet= new xAOD::Jet();
			float dRmin=999;
			for( ; truth_jet_itr != truth_jet_end; ++truth_jet_itr )
			{
				xAOD::JetFourMom_t jet_truth_4mom = (*truth_jet_itr)->jetP4();
				double teta    = (jet_truth_4mom.eta());
				double tphi    = (jet_truth_4mom.phi());
				float R = DeltaR(newjet.phi(),newjet.eta(),tphi,teta);
				if (R>_dR_truth_matching || R>dRmin) continue;
				//truthMjet.makePrivateStore();
				//truthMjet.setJetP4( jet_truth_4mom );
				truthMjet = (xAOD::Jet *)(*truth_jet_itr);
				dRmin=R;
			}
			if (dRmin<_dR_truth_matching) {
				uncertprovider->CorrectJet(eventInfo, &newjet,truthMjet, cent_bin_fine,FCalEt);
				//truthMjet.releasePrivateStore();
			}
			//if (truthMjet) delete truthMjet;				
		}
		jet_pt  = (newjet.pt() * 0.001);
		jet_eta = newjet.eta();
		jet_phi = newjet.phi();
		jet_m = newjet.m()*0.001;
				
		//TODO: jet cleaning
		//		//Jet quality moment
		//		if( !m_jetCleaning->accept( **jet_itr )) Is_jet_Good.push_back(0);
		//		else Is_jet_Good.push_back(1);

		
		jet_pt_xcalib_vector.push_back(jet_pt);
		jet_phi_vector.push_back(jet_phi);
		jet_eta_vector.push_back(jet_eta);
		
		////Get rapidity
		jet4vector.SetPtEtaPhiM(jet_pt, jet_eta, jet_phi, jet_m);
		jet_y_vector.push_back(jet4vector.Rapidity());

		
		if (!isMC) //is data, so need prescales and trigger decisions (at jet level)
		{
			bool is_trig = false;
			double presc = -1;

			for (int k=_first_trigger;k<_nTriggers;k++)
			{
				if(event_isTriggered[k] && (jet_pt > jet_pt_trig[k][0] && jet_pt <= jet_pt_trig[k][1]))
				{
					is_trig = true;
					presc = h2_trigger_RunNumber_prescale->GetBinContent ( h2_trigger_RunNumber_prescale->GetXaxis()->FindBin ( eventInfo->runNumber() ),  k + 2 ); //j40 is bin 3
					//TODO temporary workaround for j40 trigger:
					if (k==_first_trigger) presc = 523.453;
					//cout << "trigger_chains: " << trigger_chains.at(k) << " bin: " << k +2 << endl; 
					break;
				}
			}
			jet_TrigPresc_vector.push_back(presc);
			jet_IsTrig_vector.push_back(is_trig);
			
		}
		//Fake jet rejection
		if (_doFJR) {
			//isFake_vector.push_back(GetFJR(jet_eta, jet_phi, track_jets));
			isFake_vector.push_back(MTCorrector::GetFJR(&newjet,m_trackSelectorTool_FJR));
		}
		
		//delete newjet;
	}

	store->clear();
	
	//delete store;

	//isolate reco jet
	jet_isolated_vector = MTCorrector::GetIsolation(jet_pt_xcalib_vector,jet_eta_vector,jet_phi_vector, 1.0, _pt_iso); // -1 is for pT_iso == jet pT
	//TODO electrons
	//	InJet_electron_pT = FilterJetsByElectrons(jet_eta_vector,jet_phi_vector,electrons_shallowCopy,_dR_max);

	//if MC, match to truth
	if (_data_switch) TruthJetIndex = TruthMatching(jet_pt_xcalib_vector,jet_eta_vector,jet_phi_vector,
											   truth_jet_pt_vector,truth_jet_eta_vector,truth_jet_phi_vector,
											   _dR_truth_matching);

	//Get the R=0.2 jets
	if (_reco_jet_collection.find("HI") != std::string::npos)
	{
		const xAOD::JetContainer* a2_jets = 0;
		EL_RETURN_CHECK("execute()",event->retrieve( a2_jets, "DFAntiKt2HIJets" ));

		// loop over the jets in the container
		xAOD::JetContainer::const_iterator a2_jet_itr = a2_jets->begin();
		xAOD::JetContainer::const_iterator a2_jet_end = a2_jets->end();

		for( ; a2_jet_itr != a2_jet_end; ++a2_jet_itr )
		{
			antikt2_pt.push_back((*a2_jet_itr)->pt());
			antikt2_eta.push_back((*a2_jet_itr)->eta());
			antikt2_phi.push_back((*a2_jet_itr)->phi());
		}
	}

	std::vector<float> trk_good_eta;
	std::vector<float> trk_good_phi;
	std::vector<float> trk_good_pt;
	
	//Loop over tracks to exclude UE cones
	for (const auto& trk : *recoTracks) {
			//get the tracks....
			float pt = trk->pt()/1000.;
			//cout << "pt " << pt << endl;
			float eta = trk->eta();
			float phi = trk->phi();
			//correct the alignement
			bool isMatchedToTruthParticle = false;
			if (isMC) {
				ElementLink< xAOD::TruthParticleContainer > truthLink = trk->auxdata<ElementLink< xAOD::TruthParticleContainer > >("truthParticleLink");				   
				float mcprob =trk->auxdata<float>("truthMatchProbability");
				if(truthLink.isValid() && mcprob > _mcProbCut) isMatchedToTruthParticle = true;
			}
			if (_correctTrackpT && (!isMC || !isMatchedToTruthParticle)) trkcorr->correctChTrackpT(pt, eta, phi, trk->charge());
			if (fabs(eta) > 2.5) continue;
			//deriv_val->Fill(pt,eventInfo->runNumber(),eventInfo->lumiBlock());
			//Charge cut
			if (_useCharge!=0 && ((int)trk->charge())!=_useCharge) continue; 
			if (isMC && uncertprovider->uncert_class==5) uncertprovider->UncerTrackMomentum(pt, eta, phi, trk->charge() );
			if (pt < trkcorr->trkpTThreshold) continue; //min pT cut
			double d0 = trk->d0();
			double d0_cut = f_d0_cut->Eval(pt);
			if(fabs(d0) > d0_cut) continue; //pT dependant d0 cut
			if(!m_trackSelectorTool->accept(*trk)) continue; //track selector tool
			
			//Additional track selection
			//if (!trkcorr->PassTracktoJetBalance(pt, jet_pt, eta, jet_eta,cent_bin_fine)) continue;	
			
			trk_good_eta.push_back(eta);
			trk_good_phi.push_back(phi);
			trk_good_pt.push_back(pt);
	}
	
	//Construct UE disitrbutions
	uee->ExcludeConesByJetandTrack(trk_good_pt,trk_good_eta,trk_good_phi,jet_pt_xcalib_vector,jet_eta_vector,jet_phi_vector);
	//uee->ExcludeConesByTrk(trk_good_pt,trk_good_eta,trk_good_phi);
	trk_good_eta.clear();
	trk_good_phi.clear();
	trk_good_pt.clear();
	
	for(unsigned int i=0; i<jet_pt_xcalib_vector.size(); i++)
	{

		//Reco jets entereing the analysis are
		//in MC: isolated and matched to isolated truth
		//in Data: isolated, and have passed the trigger requirement
		
		float jet_weight =1.;		
		jet_weight *=event_weight;
		
		jet_pt = jet_pt_xcalib_vector.at(i);
		jet_eta = jet_eta_vector.at(i);
		jet_y = jet_y_vector.at(i);
		jet_phi = jet_phi_vector.at(i);
		
		if (fabs(jet_eta)>_jet_y_cut) continue;
		if (jet_pt < _pTjetCut) continue;
		
		//if (fabs(jet_eta) > 2.) cout << "diff " << jet_y-jet_eta << " jet_y " << jet_y << " jet_eta " << jet_eta << " pt " << jet_pt << endl;
		
		if (fabs(jet_y)>_jet_y_cut) continue; //cut on rapidity (simultaniously with 2.1 on pseudorapidity)
		if (_doFJR && isFake_vector.at(i)) continue;
				
		if (_data_switch==0)
		{
			if (!jet_IsTrig_vector.at(i)) continue;
			jet_weight = jet_TrigPresc_vector.at(i);
		}
		else if (_data_switch==1)
		{
			int truthindex=TruthJetIndex.at(i);
			if (truthindex<0) continue; //Matching to truh jets
			if (!jetcorr->MCJetJERClean(truth_jet_pt_vector.at(truthindex),jet_pt,truth_jet_eta_vector.at(truthindex),cent_bin_fine) ) continue; //cut on JER balance
			//Reweighting
			if (_applyReweighting) jet_weight*=jetcorr->GetJetReweightingFactor(truth_jet_pt_vector.at(truthindex),truth_jet_eta_vector.at(truthindex),cent_bin);
		}
		
		//Run-by-run stability
		if (jet_pt > 126 && jet_pt < 158){
			h_Njets_v_Run->Fill(GetRunNumberBin(eventInfo->runNumber())+0.5,jet_weight);	
		}
		
		for(unsigned int j=0; j<jet_pt_xcalib_vector.size(); j++)
		{
			if (j==i) continue;
			if (jet_pt < jet_pt_xcalib_vector.at(j)) continue;
			h_RdR.at(cent_bin)->Fill(jet_pt,jet_pt_xcalib_vector.at(j),DeltaR(jet_phi_vector.at(j),jet_eta_vector.at(j),jet_phi,jet_eta),jet_weight);
		}
		if(!jet_isolated_vector.at(i)) continue;
		
		int y_bin = jetcorr->GetJetYBin(jet_y);

		//		jet_electron_pT = InJet_electron_pT.at(i);
		//		if (isMC)
		//		{
		//			jet_uJER= jet_uJER_vector.at(i);
		//			jet_uJES= jet_uJES_vector.at(i);
		//		}

		if (_doJPRCorrection){	
			float dR_R2R4=R2Matching(jet_eta, jet_phi, antikt2_eta, antikt2_phi, 0.3); //Matching to R = 0.2 jets and adjusting jet eta/phi if matched
			h_R2vR4.at(cent_bin)->Fill(jet_pt, dR_R2R4,jet_weight);
		}
		//fill ff normalization histogram
		h_reco_jet_spectrum.at(y_bin).at(cent_bin)->Fill(jet_pt, jet_weight);
		h_reco_jet_spectrum.at(jetcorr->nJetYBins - 1).at(cent_bin)->Fill(jet_pt, jet_weight);
		h_reco_jet_spectrum_fine.at(cent_bin)->Fill(jet_pt, jet_weight);
		
		h_reco_jet_spectrum_counts.at(y_bin).at(cent_bin)->Fill(jet_pt);
		h_reco_jet_spectrum_counts.at(jetcorr->nJetYBins - 1).at(cent_bin)->Fill(jet_pt);
		
		//Reponses
		if (_data_switch==1) {
			bool passed_truth_cut = true;
			if (truth_jet_pt_vector.at(TruthJetIndex.at(i)) < _truthpTjetCut) passed_truth_cut = false;
			if (fabs(truth_jet_eta_vector.at(TruthJetIndex.at(i)))>_jet_y_cut) passed_truth_cut = false;
			if (passed_truth_cut)  {
				ff_jetResponse.at(y_bin).at(cent_bin)->Fill(jet_pt, truth_jet_pt_vector.at(TruthJetIndex.at(i)), jet_weight );
				ff_jetResponse.at(jetcorr->nJetYBins - 1).at(cent_bin)->Fill(jet_pt, truth_jet_pt_vector.at(TruthJetIndex.at(i)), jet_weight );
				response_jet.at(y_bin).at(cent_bin)->Fill(jet_pt, truth_jet_pt_vector.at(TruthJetIndex.at(i)), jet_weight );
				response_jet.at(jetcorr->nJetYBins - 1).at(cent_bin)->Fill(jet_pt, truth_jet_pt_vector.at(TruthJetIndex.at(i)), jet_weight );
				
				h_reco_jet_spectrum_counts_TMR.at(y_bin).at(cent_bin)->Fill(jet_pt);
				h_reco_jet_spectrum_counts_TMR.at(jetcorr->nJetYBins - 1).at(cent_bin)->Fill(jet_pt);
			}	
		}
		
		float z_max=0;
		float pT_max=0;
		int trk_multiplicity[10]; for (int nMultThreshold=0; nMultThreshold<trkcorr->nMultThresholds; nMultThreshold++) trk_multiplicity[nMultThreshold]=0;
		float UE_sum=0;
		float altUE_sum=0;
		for (const auto& trk : *recoTracks) {
			//get the tracks....
			float pt = trk->pt()/1000.;
			float eta = trk->eta();
			float phi = trk->phi();
			bool isMatchedToTruthParticle = false;
			if (isMC) {
				ElementLink< xAOD::TruthParticleContainer > truthLink = trk->auxdata<ElementLink< xAOD::TruthParticleContainer > >("truthParticleLink");				   
				float mcprob =trk->auxdata<float>("truthMatchProbability");
				if(truthLink.isValid() && mcprob > _mcProbCut) isMatchedToTruthParticle = true;
			}
			if (_correctTrackpT && (!isMC || !isMatchedToTruthParticle)) trkcorr->correctChTrackpT(pt, eta, phi, trk->charge());
			if (fabs(eta) > 2.5) continue;
			//Moving track pt in systematic variation  if needed
			if (_useCharge!=0 && ((int)trk->charge())!=_useCharge) continue;
			if (isMC && uncertprovider->uncert_class==5) uncertprovider->UncerTrackMomentum(pt, eta, phi, trk->charge() );
			if (pt < trkcorr->trkpTThreshold) continue; //min pT cut
			
			float R = DeltaR(phi,eta,jet_phi,jet_eta);
			
			//Tracking validation histograms
			if (jet_pt>251. && jet_pt<316. && R < _dR_max) {
				//track parameters
				int nPixHits = trk->auxdata< unsigned char >("numberOfPixelHits") + trk->auxdata< unsigned char >("numberOfPixelDeadSensors");
				int nSCTHits = trk->auxdata< unsigned char >("numberOfSCTHits") + trk->auxdata< unsigned char >("numberOfSCTDeadSensors");
				double d0 = trk->d0();
				double theta = trk->theta();
				double z0pv=(trk->z0()+trk->vz()-(*vtx_itr)->z())*sin(theta);	// pp: trk->z0() - w.r.t. BS
			
				h_PixHits.at(cent_bin)->Fill(pt,eta,nPixHits);
				h_SCTHits.at(cent_bin)->Fill(pt,eta,nSCTHits);
				h_d0.at(cent_bin)->Fill(pt,d0);
				h_z0sintheta.at(cent_bin)->Fill(pt,z0pv);
				h_reco_trk_map_nocuts->Fill(pt,eta,phi);
			}

			double d0 = trk->d0();
			double d0_cut = f_d0_cut->Eval(pt);
			if(fabs(d0) > d0_cut) continue; //pT dependant d0 cut
			if(!m_trackSelectorTool->accept(*trk)) continue; //track selector tool
			
			//Additional track selection
			if (!trkcorr->PassTracktoJetBalance(pt, jet_pt, eta, jet_eta,cent_bin_corse)) continue;
			
			//Efficiency correction;
			float eff_uncertainty = 0;
			if (_uncert_index > 0 && uncertprovider->uncert_class==4) eff_uncertainty = uncertprovider->CorrectTrackEff(jet_pt, jet_y, pt,eta, R, cent_bin_corse); 
			float eff_weight = trkcorr->get_effcorr(pt, eta, cent_bin_corse, eff_uncertainty, jet_pt, jet_y, R); //TODO efficiency down to 1 GeV
			
			double z = GetZ(pt,eta,phi,jet_pt,jet_eta,jet_phi,_UseAltzDef);
			
			//required to be within jet, need to be separated for UEEstimator
			if (R < _dR_max) {				
				ff_raw.at(y_bin).at(cent_bin)->Fill(z,jet_pt, jet_weight*eff_weight);
				ff_raw.at(jetcorr->nJetYBins - 1).at(cent_bin)->Fill(z,jet_pt, jet_weight*eff_weight);
				ChPS_raw.at(y_bin).at(cent_bin)->Fill(pt,jet_pt, jet_weight*eff_weight);
				ChPS_raw.at(jetcorr->nJetYBins - 1).at(cent_bin)->Fill(pt,jet_pt, jet_weight*eff_weight);
				
				//Run-by-run stability
				if ((jet_pt > 126 && jet_pt < 158) && (z > 0.1 && z<0.2 )){
					h_FF_v_Run->Fill(GetRunNumberBin(eventInfo->runNumber())+0.5,jet_weight*eff_weight);	
				}

				ff_raw_fine.at(y_bin).at(cent_bin)->Fill(z,jet_pt, jet_weight*eff_weight);
				ff_raw_fine.at(jetcorr->nJetYBins - 1).at(cent_bin)->Fill(z,jet_pt, jet_weight*eff_weight);
				ChPS_raw_fine.at(y_bin).at(cent_bin)->Fill(pt,jet_pt, jet_weight*eff_weight);
				ChPS_raw_fine.at(jetcorr->nJetYBins - 1).at(cent_bin)->Fill(pt,jet_pt, jet_weight*eff_weight);
			
				//Trackibg validation histograms
				if (jet_pt>80. && jet_pt<110.) h_reco_trk_map->Fill(pt,eta,phi);
				
				for (int nMultThreshold=0; nMultThreshold<trkcorr->nMultThresholds; nMultThreshold++) {if (pt > trkcorr->MultThresholds[nMultThreshold]) {trk_multiplicity[nMultThreshold]++;}}
			}
			
			//UE distributions
			uee->FindCone(pt,eta,phi);
			
		    float deltaRBkgr = uee->GetDeltaRToConeAxis();
		    float maxpT = uee->GetMaxConepT();
		    int ConeIndex = uee->GetMaxConeIndex();
		    if (deltaRBkgr <= _dR_max){ 
		           
		           //cout << " max pt " << maxpT <<  " deltaRBkgr " << deltaRBkgr << " ConeIndex " << ConeIndex << endl;
		           
		           float w_eta  = uee->CalculateEtaWeight(pt,eta,jet_eta,cent_bin_fine);
				   float w_ncones = uee->GetNConesWeight();
				   float w_flow=1;
				   if (_dataset==4) {
				   	w_flow = uee->CalculateFlowWeight( pt, eta, phi, jet_phi,  FCalEt );
				   	//w_flow = w_flow * uee->CalculateV3Weight( pt, eta, phi, jet_phi,  FCalEt );
				   }
				   float w_bkgr = w_eta * w_ncones * w_flow;
		    		
		           float EtaBkgr = uee->GetetaOfConeAxis();
		           float PhiBkgr = uee->GetphiOfConeAxis();
		           float z_UE = GetZ(pt,eta,phi,jet_pt,EtaBkgr,PhiBkgr,_UseAltzDef);

		           ff_raw_UE.at(y_bin).at(cent_bin)->Fill(z_UE,jet_pt, jet_weight*eff_weight*w_bkgr);
				   ff_raw_UE.at(jetcorr->nJetYBins - 1).at(cent_bin)->Fill(z_UE,jet_pt, jet_weight*eff_weight*w_bkgr);
				   ChPS_raw_UE.at(y_bin).at(cent_bin)->Fill(pt,jet_pt, jet_weight*eff_weight*w_bkgr);
				   ChPS_raw_UE.at(jetcorr->nJetYBins - 1).at(cent_bin)->Fill(pt,jet_pt, jet_weight*eff_weight*w_bkgr);
					
				   ff_raw_UE_fine.at(y_bin).at(cent_bin)->Fill(z_UE,jet_pt, jet_weight*eff_weight*w_bkgr);
				   ff_raw_UE_fine.at(jetcorr->nJetYBins - 1).at(cent_bin)->Fill(z_UE,jet_pt, jet_weight*eff_weight*w_bkgr);
				   ChPS_raw_UE_fine.at(y_bin).at(cent_bin)->Fill(pt,jet_pt, jet_weight*eff_weight*w_bkgr);
				   ChPS_raw_UE_fine.at(jetcorr->nJetYBins - 1).at(cent_bin)->Fill(pt,jet_pt, jet_weight*eff_weight*w_bkgr);
				   
				   //UE versus jet response
				   if (_data_switch==1){
				   		if (pt > 1.0 && pt < 3.98) UE_sum+=(pt*eff_weight*w_bkgr);
				   }
				   
			}
			

			if(_data_switch==1 && R < _dR_max)  //look for truth tracks in MC, used for response matrices
			{
				//Only truth jets > 40 GeV and <2.1 in responses 
				bool passed_truth_cut = true;
				if (truth_jet_pt_vector.at(TruthJetIndex.at(i)) < _truthpTjetCut) passed_truth_cut = false;
				if (fabs(truth_jet_eta_vector.at(TruthJetIndex.at(i)))>_jet_y_cut) passed_truth_cut = false;
				if (!passed_truth_cut)  continue;
				
				bool isFake=true;
				ElementLink< xAOD::TruthParticleContainer > truthLink = trk->auxdata<ElementLink< xAOD::TruthParticleContainer > >("truthParticleLink");				   
				//cout << "truthLink.isValid() "  << truthLink.isValid() << endl;
				float mcprob =trk->auxdata<float>("truthMatchProbability");
				if(truthLink.isValid() && mcprob > _mcProbCut){	
					
					int trktype = getTypeReco((*truthLink)->barcode(),(*truthLink)->pdgId(),(*truthLink)->status(),(*truthLink)->charge(),mcprob,_mcProbCut);
					//cout << "trktype "  << trktype << endl;
					if ( trktype == 1  || trktype == 5 )  {					
						float track_mc_pt=(*truthLink)->pt()*0.001; //bring into GeV
						float track_mc_eta=(*truthLink)->eta();
						float track_mc_phi=(*truthLink)->phi();
						int track_mc_pdg=(*truthLink)->pdgId();
						int track_mc_barcode=(*truthLink)->barcode();
						float track_mc_charge= (*truthLink)->threeCharge()/3;
						
						float R_truth = DeltaR(track_mc_phi,track_mc_eta,truth_jet_phi_vector.at(TruthJetIndex.at(i)),truth_jet_eta_vector.at(TruthJetIndex.at(i)) );
						double z_truth = GetZ(track_mc_pt, track_mc_eta,track_mc_phi, truth_jet_pt_vector.at(TruthJetIndex.at(i)),truth_jet_eta_vector.at(TruthJetIndex.at(i)),truth_jet_phi_vector.at(TruthJetIndex.at(i)),_UseAltzDef);
						
						float dpT_weight =1.;
						float z_weight =1.;
						float dpT_weight_fine =1.;
						float z_weight_fine =1.;
						if (_applyReweighting) {	
							dpT_weight = jetcorr->GetCHPSReweightingFactor(track_mc_pt, truth_jet_pt_vector.at(TruthJetIndex.at(i)), truth_jet_eta_vector.at(TruthJetIndex.at(i)), cent_bin,0); 
							z_weight = jetcorr->GetFFReweightingFactor(z_truth, truth_jet_pt_vector.at(TruthJetIndex.at(i)), truth_jet_eta_vector.at(TruthJetIndex.at(i)), cent_bin,0);
							dpT_weight_fine = jetcorr->GetCHPSReweightingFactor(track_mc_pt, truth_jet_pt_vector.at(TruthJetIndex.at(i)), truth_jet_eta_vector.at(TruthJetIndex.at(i)), cent_bin,1); 
							z_weight_fine = jetcorr->GetFFReweightingFactor(z_truth, truth_jet_pt_vector.at(TruthJetIndex.at(i)), truth_jet_eta_vector.at(TruthJetIndex.at(i)), cent_bin,1);
							if ((isnan(z_weight) || isnan(dpT_weight) ) || (z_weight!=z_weight) || (z_weight< 0.4 || z_weight > 2.)) cout << endl << endl << "WARNING: reweighting factor is NaN or strange size of " << z_weight << endl << endl;
							if ((isnan((float)z) || isnan((float)z_truth) ) || ((z!=z) || (z_truth!=z_truth)) || (z< 0. || z > 2.) || (z_truth< 0. || z_truth > 2.)) cout << endl << endl << "WARNING: potentional problem with z "  << endl << endl;
						}
						
						//cout << "jet " << truth_jet_pt_vector.at(TruthJetIndex.at(i)) << " track " << track_mc_pt <<" dpT_weight " << dpT_weight << " z_weight " << z_weight << " eff_weight " << eff_weight << endl;
						
						ff_trackpTResponse.at(y_bin).at(cent_bin)->Fill(pt, track_mc_pt, truth_jet_pt_vector.at(TruthJetIndex.at(i)), jet_weight*eff_weight*dpT_weight );
						ff_trackpTResponse.at(jetcorr->nJetYBins - 1).at(cent_bin)->Fill(pt, track_mc_pt, truth_jet_pt_vector.at(TruthJetIndex.at(i)), jet_weight*eff_weight*dpT_weight );
						
						ff_trackzResponse.at(y_bin).at(cent_bin)->Fill(z, z_truth, truth_jet_pt_vector.at(TruthJetIndex.at(i)), jet_weight*eff_weight *z_weight );
						ff_trackzResponse.at(jetcorr->nJetYBins - 1).at(cent_bin)->Fill(z, z_truth, truth_jet_pt_vector.at(TruthJetIndex.at(i)), jet_weight*eff_weight*z_weight );
						
						ff_trackzResponse_counts.at(y_bin).at(cent_bin)->Fill(z, z_truth, truth_jet_pt_vector.at(TruthJetIndex.at(i)));
						ff_trackzResponse_counts.at(jetcorr->nJetYBins - 1).at(cent_bin)->Fill(z, z_truth, truth_jet_pt_vector.at(TruthJetIndex.at(i)));
						
						if (z_truth > 1.) {ff_zg1_v_PDG.at(cent_bin)->Fill(z_truth,fabs(track_mc_pdg)); if (fabs(track_mc_pdg>6000)) cout << "particle " << track_mc_pdg << " with truth z: " << z_truth << endl;}
						
						response_ChPS.at(y_bin).at(cent_bin)->Fill(pt, jet_pt, track_mc_pt, truth_jet_pt_vector.at(TruthJetIndex.at(i)), jet_weight*eff_weight*dpT_weight );
						response_ChPS.at(jetcorr->nJetYBins - 1).at(cent_bin)->Fill(pt, jet_pt, track_mc_pt , truth_jet_pt_vector.at(TruthJetIndex.at(i)), jet_weight*eff_weight*dpT_weight );
						
						response_FF.at(y_bin).at(cent_bin)->Fill(z, jet_pt , z_truth, truth_jet_pt_vector.at(TruthJetIndex.at(i)), jet_weight*eff_weight*z_weight );
						response_FF.at(jetcorr->nJetYBins - 1).at(cent_bin)->Fill(z, jet_pt, z_truth, truth_jet_pt_vector.at(TruthJetIndex.at(i)), jet_weight*eff_weight*z_weight );
						
						response_ChPS_fine.at(y_bin).at(cent_bin)->Fill(pt, jet_pt, track_mc_pt, truth_jet_pt_vector.at(TruthJetIndex.at(i)), jet_weight*eff_weight*dpT_weight_fine );
						response_ChPS_fine.at(jetcorr->nJetYBins - 1).at(cent_bin)->Fill(pt, jet_pt, track_mc_pt , truth_jet_pt_vector.at(TruthJetIndex.at(i)), jet_weight*eff_weight*dpT_weight_fine );
						
						response_FF_fine.at(y_bin).at(cent_bin)->Fill(z, jet_pt , z_truth, truth_jet_pt_vector.at(TruthJetIndex.at(i)), jet_weight*eff_weight*z_weight_fine );
						response_FF_fine.at(jetcorr->nJetYBins - 1).at(cent_bin)->Fill(z, jet_pt, z_truth, truth_jet_pt_vector.at(TruthJetIndex.at(i)), jet_weight*eff_weight*z_weight_fine );
												
						//Counts for truncation
						response_ChPS.at(y_bin).at(cent_bin)->Fill(pt, jet_pt, track_mc_pt, truth_jet_pt_vector.at(TruthJetIndex.at(i)));
						response_ChPS.at(jetcorr->nJetYBins - 1).at(cent_bin)->Fill(pt, jet_pt, track_mc_pt , truth_jet_pt_vector.at(TruthJetIndex.at(i)));
						
						response_FF.at(y_bin).at(cent_bin)->Fill(z, jet_pt , z_truth, truth_jet_pt_vector.at(TruthJetIndex.at(i)));
						response_FF.at(jetcorr->nJetYBins - 1).at(cent_bin)->Fill(z, jet_pt, z_truth, truth_jet_pt_vector.at(TruthJetIndex.at(i)));
						
						response_ChPS_fine.at(y_bin).at(cent_bin)->Fill(pt, jet_pt, track_mc_pt, truth_jet_pt_vector.at(TruthJetIndex.at(i)));
						response_ChPS_fine.at(jetcorr->nJetYBins - 1).at(cent_bin)->Fill(pt, jet_pt, track_mc_pt , truth_jet_pt_vector.at(TruthJetIndex.at(i)));
						
						response_FF_fine.at(y_bin).at(cent_bin)->Fill(z, jet_pt , z_truth, truth_jet_pt_vector.at(TruthJetIndex.at(i)));
						response_FF_fine.at(jetcorr->nJetYBins - 1).at(cent_bin)->Fill(z, jet_pt, z_truth, truth_jet_pt_vector.at(TruthJetIndex.at(i)));
						
						
						isFake=false;
						
						//Only for truth matched tracks
						if (z_max < z) z_max = z;
						if (pT_max < pt) pT_max = pt;
					}
				 }
				 //Fake/UE tracks
				 if (isFake){
				 	//TODO check weighting of UE
				 	ff_UE_z.at(y_bin).at(cent_bin)->Fill(z,jet_pt, jet_weight*eff_weight);
					ff_UE_z.at(jetcorr->nJetYBins - 1).at(cent_bin)->Fill(z,jet_pt, jet_weight*eff_weight);
					ff_UE_pT.at(y_bin).at(cent_bin)->Fill(pt,jet_pt, jet_weight*eff_weight);
					ff_UE_pT.at(jetcorr->nJetYBins - 1).at(cent_bin)->Fill(pt,jet_pt, jet_weight*eff_weight);
					
					ff_UE_z_fine.at(y_bin).at(cent_bin)->Fill(z,jet_pt, jet_weight*eff_weight);
					ff_UE_z_fine.at(jetcorr->nJetYBins - 1).at(cent_bin)->Fill(z,jet_pt, jet_weight*eff_weight);
					ff_UE_pT_fine.at(y_bin).at(cent_bin)->Fill(pt,jet_pt, jet_weight*eff_weight);
					ff_UE_pT_fine.at(jetcorr->nJetYBins - 1).at(cent_bin)->Fill(pt,jet_pt, jet_weight*eff_weight);
					
					if (pt > 1.0 && pt < 3.98) altUE_sum+=pt*eff_weight;
				 }

			}

		} // end reco track loop
		//Alternative UE
		int i_dPsi = GetPsiBin(DeltaPsi(jet_phi, uee->Psi));
		for (int i_dR = 0; i_dR < 1; i_dR++) //Up to 0.4
		{
			for (int i_pt = 0; i_pt < 7; i_pt++)
			{
				double UE_err = -1;
				double UE_val = uee->getShapeUE(i_dR, i_dPsi, i_pt, cent_bin, jet_eta, jet_phi, UE_err);
				double trk_bin_center = ChPS_MB_raw_UE.at(0).at(0)->GetXaxis()->GetBinCenter(i_pt+1);
				//cout << " c" << cent_bin <<" dr " << i_dR << " i_dPsi " << i_dPsi << " i_pt " << i_pt << " trk_bin_center " << trk_bin_center << " jet_eta " << jet_eta << " jet_phi " << jet_phi << " UE " << UE_val << endl;
				ChPS_MB_raw_UE.at(y_bin).at(cent_bin)->Fill(trk_bin_center, jet_pt, UE_val*jet_weight);
				ChPS_MB_raw_UE.at(jetcorr->nJetYBins - 1).at(cent_bin)->Fill(trk_bin_center, jet_pt, UE_val*jet_weight);
				ChPS_MB_raw_UE_err.at(y_bin).at(cent_bin)->Fill(trk_bin_center, jet_pt, UE_err*jet_weight);
				ChPS_MB_raw_UE_err.at(jetcorr->nJetYBins - 1).at(cent_bin)->Fill(trk_bin_center, jet_pt, UE_err*jet_weight);
			}
		}
		
		
		//JES plots
		if(_data_switch==1) {
			ff_UE_pT_fine_response.at(y_bin).at(cent_bin)->Fill(altUE_sum,jet_pt, (jet_pt - truth_jet_pt_vector.at(TruthJetIndex.at(i))) /truth_jet_pt_vector.at(TruthJetIndex.at(i)), jet_weight);
			ff_UE_pT_fine_response.at(jetcorr->nJetYBins - 1).at(cent_bin)->Fill(altUE_sum,jet_pt,(jet_pt - truth_jet_pt_vector.at(TruthJetIndex.at(i)))/truth_jet_pt_vector.at(TruthJetIndex.at(i)), jet_weight);				 
			ff_UE_pT_fine_response.at(y_bin).at(GetCentralityNBins(_centrality_scheme)-1)->Fill(altUE_sum,jet_pt, (jet_pt - truth_jet_pt_vector.at(TruthJetIndex.at(i))) /truth_jet_pt_vector.at(TruthJetIndex.at(i)), jet_weight);
			ff_UE_pT_fine_response.at(jetcorr->nJetYBins - 1).at(GetCentralityNBins(_centrality_scheme)-1)->Fill(altUE_sum,jet_pt,(jet_pt - truth_jet_pt_vector.at(TruthJetIndex.at(i)))/truth_jet_pt_vector.at(TruthJetIndex.at(i)), jet_weight);				 
			
			
			ChPS_raw_UE_fine_response.at(y_bin).at(cent_bin)->Fill(UE_sum,jet_pt, (jet_pt- truth_jet_pt_vector.at(TruthJetIndex.at(i)))/truth_jet_pt_vector.at(TruthJetIndex.at(i)), jet_weight);
			ChPS_raw_UE_fine_response.at(jetcorr->nJetYBins - 1).at(cent_bin)->Fill(UE_sum,jet_pt, (jet_pt- truth_jet_pt_vector.at(TruthJetIndex.at(i)))/truth_jet_pt_vector.at(TruthJetIndex.at(i)), jet_weight);			   	
			ChPS_raw_UE_fine_response.at(y_bin).at(GetCentralityNBins(_centrality_scheme)-1)->Fill(UE_sum,jet_pt, (jet_pt- truth_jet_pt_vector.at(TruthJetIndex.at(i)))/truth_jet_pt_vector.at(TruthJetIndex.at(i)), jet_weight);
			ChPS_raw_UE_fine_response.at(jetcorr->nJetYBins - 1).at(GetCentralityNBins(_centrality_scheme)-1)->Fill(UE_sum,jet_pt, (jet_pt- truth_jet_pt_vector.at(TruthJetIndex.at(i)))/truth_jet_pt_vector.at(TruthJetIndex.at(i)), jet_weight);			   	
			
			
			UE_pT_Mcorrelation_v_response_fine.at(y_bin).at(cent_bin)->Fill((altUE_sum-UE_sum)/UE_sum,jet_pt, (jet_pt - truth_jet_pt_vector.at(TruthJetIndex.at(i)))/truth_jet_pt_vector.at(TruthJetIndex.at(i)) , jet_weight);
			UE_pT_Mcorrelation_v_response_fine.at(jetcorr->nJetYBins - 1).at(cent_bin)->Fill((altUE_sum-UE_sum)/UE_sum,jet_pt,(jet_pt- truth_jet_pt_vector.at(TruthJetIndex.at(i)))/truth_jet_pt_vector.at(TruthJetIndex.at(i)), jet_weight);
			UE_pT_Mcorrelation_v_response_fine.at(y_bin).at(GetCentralityNBins(_centrality_scheme)-1)->Fill((altUE_sum-UE_sum)/UE_sum,jet_pt, (jet_pt - truth_jet_pt_vector.at(TruthJetIndex.at(i)))/truth_jet_pt_vector.at(TruthJetIndex.at(i)) , jet_weight);
			UE_pT_Mcorrelation_v_response_fine.at(jetcorr->nJetYBins - 1).at(GetCentralityNBins(_centrality_scheme)-1)->Fill((altUE_sum-UE_sum)/UE_sum,jet_pt,(jet_pt- truth_jet_pt_vector.at(TruthJetIndex.at(i)))/truth_jet_pt_vector.at(TruthJetIndex.at(i)), jet_weight);
			
			
			JES_v_max_pT.at(cent_bin)->Fill(pT_max,(jet_pt - truth_jet_pt_vector.at(TruthJetIndex.at(i)) )/truth_jet_pt_vector.at(TruthJetIndex.at(i)),truth_jet_pt_vector.at(TruthJetIndex.at(i)),jet_weight);
			JES_v_max_z.at(cent_bin)->Fill(z_max,(jet_pt - truth_jet_pt_vector.at(TruthJetIndex.at(i)) )/truth_jet_pt_vector.at(TruthJetIndex.at(i)),truth_jet_pt_vector.at(TruthJetIndex.at(i)),jet_weight);
		}
		for (int nMultThreshold=0;nMultThreshold<trkcorr->nMultThresholds;nMultThreshold++) {h_jetpT_v_multiplicity.at(cent_bin)->Fill(jet_pt,nMultThreshold,trk_multiplicity[nMultThreshold]);}		
	}// end reco jet loop
	
	

	//Loop over truth jets
	if(_data_switch == 1)
	{
		const xAOD::TruthParticleContainer * particles = 0;
		EL_RETURN_CHECK("execute",event->retrieve( particles, "TruthParticles"));


		//Count truth jets for alternative FF evaluation
		/*
		for(unsigned int i=0; i<truth_jet_pt_vector.size(); i++)
		{
			float jet_weight = 1.;		
		    jet_weight *=event_weight;

			truth_jet_pt = truth_jet_pt_vector.at(i);
			truth_jet_y = truth_jet_y_vector.at(i);
			truth_jet_eta = truth_jet_eta_vector.at(i);
			
			if (truth_jet_pt< _truthpTjetCut) continue;
			if (fabs(truth_jet_eta)>2.1) continue;
			if (fabs(truth_jet_y)>2.1) continue; //cut on rapidity (simultaniously with 2.1 on pseudorapidity)
			if(!truth_jet_isolated_vector.at(i)) continue;
			h_event_jet_counts->Fill(truth_jet_pt, jet_weight);		
		}	
		*/
		for(unsigned int i=0; i<truth_jet_pt_vector.size(); i++)
		{

			float jet_weight = 1.;		
		    jet_weight *=event_weight;

			truth_jet_pt = truth_jet_pt_vector.at(i);
			truth_jet_eta = truth_jet_eta_vector.at(i);
			truth_jet_y = truth_jet_y_vector.at(i);
			truth_jet_phi = truth_jet_phi_vector.at(i);			
			int y_bin = jetcorr->GetJetYBin(truth_jet_y);
			
			if (truth_jet_pt< _truthpTjetCut) continue;
			//if (fabs(truth_jet_eta)>_jet_y_cut) continue;
			if (fabs(truth_jet_y)>_jet_y_cut) continue; //cut on rapidity (simultaniously with 2.1 on pseudorapidity)
			//if(!truth_jet_isolated_vector.at(i)) continue;
			
			//TODO do we want to reweigth also truth spectrum?
			//if (_applyReweighting) jet_weight*=jetcorr->GetJetReweightingFactor(truth_jet_pt,truth_jet_eta,cent_bin);
			
			h_true_jet_spectrum.at(y_bin).at(cent_bin)->Fill(truth_jet_pt, jet_weight);
			h_true_jet_spectrum.at(jetcorr->nJetYBins - 1).at(cent_bin)->Fill(truth_jet_pt, jet_weight);
			h_true_jet_spectrum.at(y_bin).at(GetCentralityNBins(_centrality_scheme)-1)->Fill(truth_jet_pt, jet_weight);
			h_true_jet_spectrum.at(jetcorr->nJetYBins - 1).at(GetCentralityNBins(_centrality_scheme)-1)->Fill(truth_jet_pt, jet_weight);			
			h_true_jet_spectrum_fine.at(cent_bin)->Fill(truth_jet_pt, jet_weight);
			h_true_jet_spectrum_fine.at(GetCentralityNBins(_centrality_scheme)-1)->Fill(truth_jet_pt, jet_weight);

			xAOD::TruthParticleContainer::const_iterator truth_itr = particles->begin();
			xAOD::TruthParticleContainer::const_iterator truth_end = particles->end();

			for( ; truth_itr!=truth_end; ++truth_itr)
			{
				if (_useCharge!=0 && ((int)(*truth_itr)->charge())!=_useCharge) continue;
				int ty=getTypeTruth((*truth_itr)->barcode(),(*truth_itr)->pdgId(),(*truth_itr)->status(),(*truth_itr)->charge());
				if(ty!=1 && ty!=5) continue;

				//get the tracks....
				float eta = (*truth_itr)->eta();
				float phi = (*truth_itr)->phi();
				float pt = (*truth_itr)->pt()/ 1000.0;

				if (fabs(eta)>2.5) continue;
				if (fabs(pt)<trkcorr->trkpTThreshold) continue;

				//Only tracks associated with a jet
				float R = DeltaR(phi,eta,truth_jet_phi,truth_jet_eta);
				if(R > _dR_max) continue;

				//Get z
				double truth_z = GetZ(pt,eta,phi,truth_jet_pt,truth_jet_eta,truth_jet_phi,_UseAltzDef);				
				
				ff_truth.at(y_bin).at(cent_bin)->Fill(truth_z,truth_jet_pt, jet_weight);
				ff_truth.at(jetcorr->nJetYBins - 1).at(cent_bin)->Fill(truth_z,truth_jet_pt, jet_weight);
			    ChPS_truth.at(y_bin).at(cent_bin)->Fill(pt,truth_jet_pt, jet_weight);
			    ChPS_truth.at(jetcorr->nJetYBins - 1).at(cent_bin)->Fill(pt,truth_jet_pt, jet_weight);
			    
			    ff_truth_fine.at(y_bin).at(cent_bin)->Fill(truth_z,truth_jet_pt, jet_weight);
				ff_truth_fine.at(jetcorr->nJetYBins - 1).at(cent_bin)->Fill(truth_z,truth_jet_pt, jet_weight);
				//ff_truth_alt_fine.at(jetcorr->nJetYBins - 1).at(cent_bin)->Fill(truth_z,truth_jet_pt, jet_weight/h_event_jet_counts->GetBinContent(h_event_jet_counts->FindBin(truth_jet_pt)));
			    ChPS_truth_fine.at(y_bin).at(cent_bin)->Fill(pt,truth_jet_pt, jet_weight);
			    ChPS_truth_fine.at(jetcorr->nJetYBins - 1).at(cent_bin)->Fill(pt,truth_jet_pt, jet_weight);
			    
			    //Centrality inclusive
			    ff_truth.at(y_bin).at(GetCentralityNBins(_centrality_scheme)-1)->Fill(truth_z,truth_jet_pt, jet_weight);
				ff_truth.at(jetcorr->nJetYBins - 1).at(GetCentralityNBins(_centrality_scheme)-1)->Fill(truth_z,truth_jet_pt, jet_weight);
			    ChPS_truth.at(y_bin).at(GetCentralityNBins(_centrality_scheme)-1)->Fill(pt,truth_jet_pt, jet_weight);
			    ChPS_truth.at(jetcorr->nJetYBins - 1).at(GetCentralityNBins(_centrality_scheme)-1)->Fill(pt,truth_jet_pt, jet_weight);
			    
			    ff_truth_fine.at(y_bin).at(GetCentralityNBins(_centrality_scheme)-1)->Fill(truth_z,truth_jet_pt, jet_weight);
				ff_truth_fine.at(jetcorr->nJetYBins - 1).at(GetCentralityNBins(_centrality_scheme)-1)->Fill(truth_z,truth_jet_pt, jet_weight);
			    ChPS_truth_fine.at(y_bin).at(GetCentralityNBins(_centrality_scheme)-1)->Fill(pt,truth_jet_pt, jet_weight);
			    ChPS_truth_fine.at(jetcorr->nJetYBins - 1).at(GetCentralityNBins(_centrality_scheme)-1)->Fill(pt,truth_jet_pt, jet_weight);
			}
		}
	}


//	delete electrons_shallowCopy.first;
//	delete electrons_shallowCopy.second;

	jet_pt_xcalib_vector.clear();
	jet_phi_vector.clear();
	jet_eta_vector.clear();
	jet_y_vector.clear();

	truth_jet_indices.clear();
	truth_jet_pt_vector.clear();
	truth_jet_eta_vector.clear();
	truth_jet_phi_vector.clear();
	TruthJetIndex.clear();

	for (int j=0;j<_nTriggers;j++)
	{
		isTriggered[j].clear();
	}

	antikt2_pt.clear();
	antikt2_eta.clear();
	antikt2_phi.clear();

	jet_pt_xcalib_vector.shrink_to_fit();
	jet_phi_vector.shrink_to_fit();
	jet_eta_vector.shrink_to_fit();

	truth_jet_indices.shrink_to_fit();
	truth_jet_pt_vector.shrink_to_fit();
	truth_jet_eta_vector.shrink_to_fit();
	truth_jet_y_vector.shrink_to_fit();
	truth_jet_phi_vector.shrink_to_fit();
	TruthJetIndex.shrink_to_fit();

	for (int j=0;j<_nTriggers;j++){
		isTriggered[j].shrink_to_fit();
	}

	antikt2_pt.shrink_to_fit();
	antikt2_eta.shrink_to_fit();
	antikt2_phi.shrink_to_fit();
	
	//h_event_jet_counts->Reset();
	
	return EL::StatusCode::SUCCESS;
}

EL::StatusCode PbPbFragmentation :: postExecute (){
	return EL::StatusCode::SUCCESS;
}

EL::StatusCode PbPbFragmentation :: finalize (){
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


	if(m_trackSelectorTool)
	{
		delete m_trackSelectorTool;
		m_trackSelectorTool=0;
	}

	//cleaning jets
	/*
	if(jesProv)
	{
		delete jesProv;
		jesProv=0;
	}
	*/
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

EL::StatusCode PbPbFragmentation :: histFinalize (){
	cout<<"Events = "<< m_eventCounter<<endl;
	return EL::StatusCode::SUCCESS;
}
