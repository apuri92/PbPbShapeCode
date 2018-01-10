#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include "pPbFragmentation/PbPbFFShape.h"
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

ClassImp(PbPbFFShape)

#define EL_RETURN_CHECK( CONTEXT, EXP )			\
do {							\
if( ! EXP.isSuccess() ) {				\
Error( CONTEXT,					\
XAOD_MESSAGE( "Failed to execute: %s" ),	\
#EXP );					\
return EL::StatusCode::FAILURE;			\
}							\
} while( false )


PbPbFFShape :: PbPbFFShape ()
{
}

EL::StatusCode PbPbFFShape :: fileExecute ()
{
	// Here you do everything that needs to be done exactly once for every
	// single file, e.g. collect a list of all lumi-blocks processed
	return EL::StatusCode::SUCCESS;
}

EL::StatusCode PbPbFFShape :: changeInput (bool firstFile)
{
	return EL::StatusCode::SUCCESS;
}

//Loop over events
EL::StatusCode PbPbFFShape :: execute (){

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
//	if (m_eventCounter != 246) return EL::StatusCode::SUCCESS;
//	cout << "Event: " << m_eventCounter << endl;
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
	double event_weight_fcal = 1;
    int n_cent_bins = GetCentralityNBins(_centrality_scheme);

    if (_centrality_scheme>1)
	{
		//Centrality
		const xAOD::HIEventShapeContainer* calos=0;
		EL_RETURN_CHECK("execute",event->retrieve( calos, "CaloSums"));
		FCalEt=calos->at(5)->et()*1e-6;
		cent_bin = GetCentralityBin(_centrality_scheme, FCalEt,  isHIJING );
		cent_bin_fine = GetCentralityBin(_centrality_scheme, FCalEt,  isHIJING ); //Need for some tools
		if (isMC) event_weight_fcal = jetcorr->GetFCalWeight(FCalEt);
		h_centrality->Fill(cent_bin,event_weight_fcal);
		h_centrality->Fill(n_cent_bins-1,event_weight_fcal);

		//Get HI clusters for flow
		//const xAOD::CaloClusterContainer *hiclus(0);
		//EL_RETURN_CHECK("execute",event->retrieve(hiclus,"HIClusters") );
		//cout << "Psi 1: " << GetEventPlane(hiclus) << " Psi 2: " << GetEventPlane(calos) << endl;
		uee->Psi = GetEventPlane(calos);
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
			const xAOD::HIEventShapeContainer* hiev = 0;
			EL_RETURN_CHECK("execute",event->retrieve( hiev, "HIEventShape"));
			// ZDC
			m_zdcTools->reprocessZdc();

			// is Pileup
			m_is_pileup = m_hiPileup->is_pileup( *hiev, *zdcMod); // SAVE pileup Decision HERE 0 = NO pileup, 1 = pileup
		}
		//bool m_is_pileup = (FCalEt > 5.);

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
	vector<bool> jet_isolated_vector, jet_IsTrig_vector;//, isFake_vector;
	vector<float> truth_jet_eta_vector,truth_jet_phi_vector,truth_jet_pt_vector, truth_jet_y_vector;
	vector<bool> truth_jet_isolated_vector;
	vector<int> truth_jet_indices;
	vector<int> TruthJetIndex;
	vector<int> isTriggered[_nTriggers];
	vector<float> jet_uJER_vector;
	vector<vector<float>> jet_uJES_vector;

	float event_weight = 1;
	double max_pt = 1;

	TLorentzVector jet4vector;

	// ---- GETTING TRUTH JETS ----
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

			if (pt < _truthpTjetCut) continue;
//			if (fabs(eta)>2.5 - _dR_max) continue;

			//filling truth pt/eta/phi vectors
			truth_jet_pt_vector.push_back(pt);
			truth_jet_phi_vector.push_back(phi);
			truth_jet_eta_vector.push_back(eta);
			//Get rapidity
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

	//In PbPb MC official derivations, the default jet kinematics in DFAntiKt4HI are already calibrated. Nothing more is needed
	//In PbPb data official derivations, the default jet kinematics in DFAntiKt4HI are already calibrated upto EtaJes. Set Constit to Unsubtracted, and GSC to default (this is calibrated upto EtaJES). Apply the insitu calibration to default.

	//in pp MC, there is no DF container. The AntiKt4HIJets container has unsubtracted, subtracted = EMScale, default (incorrectly calibrated). To calibrate these, start from the Subtracted or EMScale. Set Consitutuent scale to usubtracted.
	//In pp data, get jets calibrated upto EtaJes as above. Set Constit to unsubtratced, GSC to the etaJes calibrated scale, and then apply insitu to the EtaJES calibrated scale

	for( ; jet_itr != jet_end; ++jet_itr )
	{
		xAOD::Jet newjet;
		newjet.makePrivateStore( **jet_itr );

		const xAOD::JetFourMom_t jet_4mom_def = newjet.jetP4();
		const xAOD::JetFourMom_t jet_4mom_subtr = newjet.jetP4("JetSubtractedScaleMomentum");
		const xAOD::JetFourMom_t jet_4mom_unsubtr = newjet.jetP4("JetUnsubtractedScaleMomentum");

		float jet_pt_subtr  = (jet_4mom_subtr.pt() * 0.001);
		float jet_pt_unsubtr  = (jet_4mom_unsubtr.pt() * 0.001);
		float jet_pt_def  = (jet_4mom_def.pt() * 0.001);
		hET_ETsub->Fill(jet_pt_def,jet_4mom_def.eta(),jet_pt_unsubtr-jet_pt_subtr);

		//PbPb Data
		if (_data_switch == 0 && _dataset == 4)
		{
			newjet.setJetP4("JetConstitScaleMomentum",jet_4mom_unsubtr); //Required
			newjet.setJetP4("JetEtaJESScaleMomentum",jet_4mom_def); //Required
			newjet.setJetP4("JetGSCScaleMomentum", jet_4mom_def); //Already etajes calibrated
			EL_RETURN_CHECK("execute()", m_jetCalibration->applyCalibration( newjet ) );
		}

		//PbPb MC
		if (_data_switch == 1 && _dataset == 4)
		{
			//do nothing
		}

		//pp
		if (_dataset == 3) //calibrate with sequence set in ShapeToolInit
		{
			newjet.setJetP4("JetConstitScaleMomentum",jet_4mom_unsubtr); //Required
			EL_RETURN_CHECK("execute()", m_jetCalibration->applyCalibration( newjet ) ); //calibrates with sequence EtaJes_Insitu for data, EtaJes for MC
		}

		float jet_pt  = (newjet.pt() * 0.001);
		float jet_eta = newjet.eta();
		float jet_phi = newjet.phi();
		float jet_m = newjet.m()*0.001;


		//Uncertainties

		if (isMC && _uncert_index > 0 && uncertprovider->uncert_class<4 || uncertprovider->uncert_class==6) {
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
				uncertprovider->CorrectJet(&newjet,truthMjet,cent_bin_fine,FCalEt);
				//truthMjet.releasePrivateStore();
			}
			//if (truthMjet) delete truthMjet;
		}

		if (fabs(jet_eta)> (2.5 - _dR_max)) continue;
		if (jet_pt <= _pTjetCut) continue;

		//TODO: jet cleaning
		//		//Jet quality moment
		//		if( !m_jetCleaning->accept( **jet_itr )) Is_jet_Good.push_back(0);
		//		else Is_jet_Good.push_back(1);

		jet_pt_xcalib_vector.push_back(jet_pt);
		jet_phi_vector.push_back(jet_phi);
		jet_eta_vector.push_back(jet_eta);

		//Get rapidity
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
		//delete newjet;
	}

	store->clear();

	//isolate reco jet
	jet_isolated_vector = MTCorrector::GetIsolation(jet_pt_xcalib_vector, jet_eta_vector,jet_phi_vector, 1.0, _pt_iso); // -1 is for pT_iso == jet pT

	//if MC, match to truth
	if (_data_switch) TruthJetIndex = TruthMatching(jet_pt_xcalib_vector,jet_eta_vector,jet_phi_vector,
													truth_jet_pt_vector,truth_jet_eta_vector,truth_jet_phi_vector,
													_dR_truth_matching);


	for(unsigned int i=0; i<jet_pt_xcalib_vector.size(); i++)
	{
		//Reco jets entering the analysis are
		//in MC: isolated and matched to isolated truth
		//in Data: isolated, and have passed the trigger requirement

		float jet_weight = 1.;
		jet_weight *=event_weight;

		jet_pt = jet_pt_xcalib_vector.at(i);
		jet_eta = jet_eta_vector.at(i);
		jet_y = jet_y_vector.at(i);
		jet_phi = jet_phi_vector.at(i);

		if (fabs(jet_y) > (2.5 - _dR_max)) continue; //cut on rapidity (simultaniously with 2.5-drmax on pseudorapidity)
		if(!jet_isolated_vector.at(i)) continue;
		if (jet_pt < _pTjetCut) continue;
		int y_bin = jetcorr->GetJetYBin(jet_y);

		if (_data_switch==0)
		{
			if (!jet_IsTrig_vector.at(i)) continue;
			jet_weight = jet_TrigPresc_vector.at(i);
		}
		else if (_data_switch==1)
		{
			h_reco_pre_truth_match.at(cent_bin)->Fill(jet_pt, jet_eta, jet_y, jet_weight);
            h_reco_pre_truth_match.at(n_cent_bins-1)->Fill(jet_pt, jet_eta, jet_y, jet_weight);

			int truthindex=TruthJetIndex.at(i);
			if (truthindex<0) continue; //Matching to truth jets

			h_reco_post_truth_match.at(cent_bin)->Fill(jet_pt, jet_eta, jet_y, jet_weight);
            h_reco_post_truth_match.at(n_cent_bins-1)->Fill(jet_pt, jet_eta, jet_y, jet_weight);

			float delta_r = DeltaR(jet_phi,jet_eta,truth_jet_phi_vector.at(truthindex),truth_jet_eta_vector.at(truthindex));
			h_reco_truth_matched.at(cent_bin)->Fill(jet_pt, truth_jet_pt_vector.at(truthindex), delta_r, jet_weight);
            h_reco_truth_matched.at(n_cent_bins-1)->Fill(jet_pt, truth_jet_pt_vector.at(truthindex), delta_r, jet_weight);
			h_reco_truth_comparison.at(cent_bin)->Fill(truth_jet_pt_vector.at(truthindex), (jet_pt/truth_jet_pt_vector.at(truthindex) - 1), truth_jet_eta_vector.at(truthindex), jet_weight);
            h_reco_truth_comparison.at(n_cent_bins-1)->Fill(truth_jet_pt_vector.at(truthindex), (jet_pt/truth_jet_pt_vector.at(truthindex) - 1), truth_jet_eta_vector.at(truthindex), jet_weight);

			if (!jetcorr->MCJetJERClean(truth_jet_pt_vector.at(truthindex),jet_pt,truth_jet_eta_vector.at(truthindex),cent_bin_fine) ) continue; //cut on JER balance
			//Reweighting
			if (_applyReweighting) jet_weight*=jetcorr->GetJetReweightingFactor(truth_jet_pt_vector.at(truthindex),truth_jet_eta_vector.at(truthindex),cent_bin_fine); //TODO cent_bin_fine -> cent_bin when available
		}

		//fill ff normalization histogram
		h_reco_jet_spectrum.at(y_bin).at(cent_bin)->Fill(jet_pt, jet_weight);
		h_reco_jet_spectrum.at(jetcorr->nJetYBins - 1).at(cent_bin)->Fill(jet_pt, jet_weight);
        h_reco_jet_spectrum.at(y_bin).at(n_cent_bins-1)->Fill(jet_pt, jet_weight);
        h_reco_jet_spectrum.at(jetcorr->nJetYBins - 1).at(n_cent_bins-1)->Fill(jet_pt, jet_weight);

		//Reponses
		if (_data_switch==1)
		{
			bool passed_truth_cuts = true;
			if (truth_jet_pt_vector.at(TruthJetIndex.at(i)) < _truthpTjetCut) passed_truth_cuts = false;
			if (!truth_jet_isolated_vector.at(TruthJetIndex.at(i))) passed_truth_cuts = false;
			if (truth_jet_y_vector.at(TruthJetIndex.at(i)) > (2.5 - _dR_max))  passed_truth_cuts = false;
			if (!passed_truth_cuts) continue;

			ff_jetResponse.at(y_bin).at(cent_bin)->Fill(jet_pt, truth_jet_pt_vector.at(TruthJetIndex.at(i)), jet_weight );
			ff_jetResponse.at(jetcorr->nJetYBins - 1).at(cent_bin)->Fill(jet_pt, truth_jet_pt_vector.at(TruthJetIndex.at(i)), jet_weight );
            ff_jetResponse.at(y_bin).at(n_cent_bins-1)->Fill(jet_pt, truth_jet_pt_vector.at(TruthJetIndex.at(i)), jet_weight );
            ff_jetResponse.at(jetcorr->nJetYBins - 1).at(n_cent_bins-1)->Fill(jet_pt, truth_jet_pt_vector.at(TruthJetIndex.at(i)), jet_weight );

			response_jet.at(y_bin).at(cent_bin)->Fill(jet_pt, truth_jet_pt_vector.at(TruthJetIndex.at(i)), jet_weight );
			response_jet.at(jetcorr->nJetYBins - 1).at(cent_bin)->Fill(jet_pt, truth_jet_pt_vector.at(TruthJetIndex.at(i)), jet_weight );
            response_jet.at(y_bin).at(n_cent_bins-1)->Fill(jet_pt, truth_jet_pt_vector.at(TruthJetIndex.at(i)), jet_weight );
            response_jet.at(jetcorr->nJetYBins - 1).at(n_cent_bins-1)->Fill(jet_pt, truth_jet_pt_vector.at(TruthJetIndex.at(i)), jet_weight );

            truth_jet_y = truth_jet_y_vector.at(TruthJetIndex.at(i));
			int y_truth_bin = jetcorr->GetJetYBin(truth_jet_y);

            h_true_jet_spectrum_matched.at(y_truth_bin).at(cent_bin)->Fill(truth_jet_pt_vector.at(TruthJetIndex.at(i)),jet_weight);
			h_true_jet_spectrum_matched.at(jetcorr->nJetYBins - 1).at(cent_bin)->Fill(truth_jet_pt_vector.at(TruthJetIndex.at(i)), jet_weight);
			h_true_jet_spectrum_matched.at(y_truth_bin).at(n_cent_bins-1)->Fill(truth_jet_pt_vector.at(TruthJetIndex.at(i)), jet_weight);
			h_true_jet_spectrum_matched.at(jetcorr->nJetYBins - 1).at(n_cent_bins-1)->Fill(truth_jet_pt_vector.at(TruthJetIndex.at(i)), jet_weight);

			h_reco_jet_spectrum_matched.at(y_bin).at(cent_bin)->Fill(jet_pt,jet_weight);
			h_reco_jet_spectrum_matched.at(jetcorr->nJetYBins - 1).at(cent_bin)->Fill(jet_pt,jet_weight);
			h_reco_jet_spectrum_matched.at(y_bin).at(n_cent_bins-1)->Fill(jet_pt,jet_weight);
			h_reco_jet_spectrum_matched.at(jetcorr->nJetYBins - 1).at(n_cent_bins-1)->Fill(jet_pt,jet_weight);

		}



		float z_max=0;
		float pT_max=0;
		int trk_multiplicity[10]; for (int nMultThreshold=0; nMultThreshold<trkcorr->nMultThresholds; nMultThreshold++) trk_multiplicity[nMultThreshold]=0;
		for (const auto& trk : *recoTracks)
		{
			//get the tracks....
			float pt = trk->pt()/1000.;
			float eta = trk->eta();
			float phi = trk->phi();
			if (_correctTrackpT && !isMC) trkcorr->correctChTrackpT(pt, eta, phi, trk->charge());
			if (fabs(eta) >= 2.5) continue;
			//Moving track pt in systematic variation  if needed
			if (_useCharge!=0 && ((int)trk->charge())!=_useCharge) continue;
			if (isMC && uncertprovider->uncert_class==5) uncertprovider->UncerTrackMomentum(pt, eta, phi, trk->charge() );
			if (pt < _pTtrkCut) continue; //min pT cut

			float R = DeltaR(phi,eta,jet_phi,jet_eta);
			if (R > 2.0) continue; //broad cut on dR to remove unecessary effcorr warnings for tracks too far from the jet anwyay

			//Tracking validation histograms
			if (jet_pt>80. && jet_pt<110. && R < _dR_max) {
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

			//Efficiency correction;
			float eff_uncertainty = 0;
			if (_uncert_index > 0 && uncertprovider->uncert_class==4) eff_uncertainty = uncertprovider->CorrectTrackEff(pt,eta, R, cent_bin);
			float eff_weight = trkcorr->get_effcorr(pt, eta, cent_bin, 0, _dataset);
			double z = cos(R)*pt / jet_pt;

			//required to be within jet, need to be separated for UEEstimator
			int dr_bin = trkcorr->GetdRBin(R);
			if (R < _dR_max)
			{
				ff_raw.at(dr_bin).at(cent_bin)->Fill(z,jet_pt, jet_weight*eff_weight);
                ff_raw.at(dr_bin).at(n_cent_bins-1)->Fill(z,jet_pt, jet_weight*eff_weight);

                ChPS_raw.at(dr_bin).at(cent_bin)->Fill(pt,jet_pt, jet_weight*eff_weight);
                ChPS_raw.at(dr_bin).at(n_cent_bins-1)->Fill(pt,jet_pt, jet_weight*eff_weight);
				//Tracking validation histograms
				if (jet_pt>80. && jet_pt<110.) h_reco_trk_map->Fill(pt,eta,phi);

				for (int nMultThreshold=0; nMultThreshold<trkcorr->nMultThresholds; nMultThreshold++) {if (pt > trkcorr->MultThresholds[nMultThreshold]) {trk_multiplicity[nMultThreshold]++;}}
			}



			if(_data_switch==1)
			{
				//look for truth tracks in MC, used for response matrices
				//Only truth jets > 40 GeV and <2.1 in responses

				float matched_truth_jet_pt = truth_jet_pt_vector.at(TruthJetIndex.at(i));


				bool isFake=true;
				ElementLink< xAOD::TruthParticleContainer > truthLink = trk->auxdata<ElementLink< xAOD::TruthParticleContainer > >("truthParticleLink");
				float mcprob =trk->auxdata<float>("truthMatchProbability");
				if(truthLink.isValid() && mcprob > _mcProbCut)
				{
					int trktype = getTypeReco((*truthLink)->barcode(),(*truthLink)->pdgId(),(*truthLink)->status(),(*truthLink)->charge(),mcprob,_mcProbCut);
					if ( trktype == 1  || trktype == 5 )
					{
						float track_mc_pt=(*truthLink)->pt()*0.001; //bring into GeV
						float track_mc_eta=(*truthLink)->eta();
						float track_mc_phi=(*truthLink)->phi();

						if (fabs(track_mc_eta) >= 2.5) continue;

						int track_mc_pdg=(*truthLink)->pdgId();
						int track_mc_barcode=(*truthLink)->barcode();
						float track_mc_charge= (*truthLink)->threeCharge()/3;

						float R_truth = DeltaR(track_mc_phi,track_mc_eta,truth_jet_phi_vector.at(TruthJetIndex.at(i)),truth_jet_eta_vector.at(TruthJetIndex.at(i)) );
						double z_truth = cos(R_truth)*track_mc_pt / matched_truth_jet_pt;

						float dpT_weight =1.;
						float z_weight =1.;
						if (_applyReweighting)
						{
							dpT_weight = jetcorr->GetCHPSReweightingFactor(track_mc_pt, matched_truth_jet_pt, truth_jet_eta_vector.at(TruthJetIndex.at(i)), cent_bin_fine); //TODO cent_bin_fine -> cent_bin when available
							z_weight = jetcorr->GetFFReweightingFactor(z_truth, matched_truth_jet_pt, truth_jet_eta_vector.at(TruthJetIndex.at(i)), cent_bin_fine); //TODO cent_bin_fine -> cent_bin when available
						}

						//R_trk_jet
						float R_reco_reco = DeltaR(phi,eta,jet_phi,jet_eta); //(reco_reco and truth_reco) and (reco_truth and truth_truth) are same
						float R_reco_truth = DeltaR(phi,eta,truth_jet_phi_vector.at(TruthJetIndex.at(i)),truth_jet_eta_vector.at(TruthJetIndex.at(i)) );
						float R_truth_reco = DeltaR(track_mc_phi,track_mc_eta,jet_phi,jet_eta);
						float R_truth_truth = DeltaR(track_mc_phi,track_mc_eta,truth_jet_phi_vector.at(TruthJetIndex.at(i)),truth_jet_eta_vector.at(TruthJetIndex.at(i)) );

						int dr_bin_reco_reco = trkcorr->GetdRBin(R_reco_reco);
						int dr_bin_reco_truth = trkcorr->GetdRBin(R_reco_truth);
						int dr_bin_truth_reco = trkcorr->GetdRBin(R_truth_reco);
						int dr_bin_truth_truth = trkcorr->GetdRBin(R_truth_truth);

						float eff_weight = trkcorr->get_effcorr(pt, eta, cent_bin, 0, _dataset);

						if (R_truth_truth < _dR_max)
						{
							int deta_bin = trkcorr->GetdRBin(fabs(DeltaEta(track_mc_eta,truth_jet_eta_vector.at(TruthJetIndex.at(i)))));
							ChPS_raw_tt_deta.at(deta_bin).at(cent_bin)->Fill(track_mc_pt,matched_truth_jet_pt, jet_weight);
                            ChPS_raw_tt_deta.at(deta_bin).at(n_cent_bins-1)->Fill(track_mc_pt,matched_truth_jet_pt, jet_weight);

							int dphi_bin = trkcorr->GetdRBin(fabs(DeltaPhi(track_mc_phi,truth_jet_phi_vector.at(TruthJetIndex.at(i)))));
							ChPS_raw_tt_dphi.at(dphi_bin).at(cent_bin)->Fill(track_mc_pt,matched_truth_jet_pt, jet_weight);
                            ChPS_raw_tt_dphi.at(dphi_bin).at(n_cent_bins-1)->Fill(track_mc_pt,matched_truth_jet_pt, jet_weight);

							ChPS_raw_tt.at(dr_bin_truth_truth).at(cent_bin)->Fill(track_mc_pt, matched_truth_jet_pt, jet_weight);
                            ChPS_raw_tt.at(dr_bin_truth_truth).at(n_cent_bins-1)->Fill(track_mc_pt, matched_truth_jet_pt, jet_weight);

						}
						if (R_reco_truth < _dR_max) ChPS_raw_rt.at(dr_bin_reco_truth).at(cent_bin)->Fill(pt, matched_truth_jet_pt, jet_weight*eff_weight);
						if (R_truth_reco < _dR_max) ChPS_raw_tr.at(dr_bin_truth_reco).at(cent_bin)->Fill(track_mc_pt, jet_pt, jet_weight);
						if (R_reco_reco < _dR_max)
                        {
                            ChPS_raw_rr.at(dr_bin_reco_reco).at(cent_bin)->Fill(pt, jet_pt, jet_weight*eff_weight);
                            ChPS_raw_rr.at(dr_bin_reco_reco).at(n_cent_bins-1)->Fill(pt, jet_pt, jet_weight*eff_weight);
                        }

						if (R_truth_truth < _dR_max) ChPS_raw_tt_mod.at(dr_bin_truth_truth).at(cent_bin)->Fill(track_mc_pt, matched_truth_jet_pt, jet_weight);

						if (R < _dR_max)
						{
							ff_trackpTResponse.at(dr_bin).at(cent_bin)->Fill(pt, track_mc_pt, matched_truth_jet_pt, jet_weight*eff_weight*dpT_weight );
							ff_trackzResponse.at(dr_bin).at(cent_bin)->Fill(z, z_truth, matched_truth_jet_pt, jet_weight*eff_weight *z_weight );

                            ff_trackpTResponse.at(dr_bin).at(n_cent_bins-1)->Fill(pt, track_mc_pt, matched_truth_jet_pt, jet_weight*eff_weight*dpT_weight );
                            ff_trackzResponse.at(dr_bin).at(n_cent_bins-1)->Fill(z, z_truth, matched_truth_jet_pt, jet_weight*eff_weight *z_weight );

							response_ChPS.at(dr_bin).at(cent_bin)->Fill(pt, jet_pt, track_mc_pt, matched_truth_jet_pt, jet_weight*eff_weight*dpT_weight );
							response_FF.at(dr_bin).at(cent_bin)->Fill(z, jet_pt , z_truth, matched_truth_jet_pt, jet_weight*eff_weight*z_weight );

                            response_ChPS.at(dr_bin).at(n_cent_bins-1)->Fill(pt, jet_pt, track_mc_pt, matched_truth_jet_pt, jet_weight*eff_weight*dpT_weight );
                            response_FF.at(dr_bin).at(n_cent_bins-1)->Fill(z, jet_pt , z_truth, matched_truth_jet_pt, jet_weight*eff_weight*z_weight );

							double truth_pt = matched_truth_jet_pt;
							int jetpt_bin = jetcorr->GetJetpTBin(truth_pt, (TAxis*)reco_posRes_ChPS.at(0).at(0)->GetYaxis());
							h_dR_change.at(jetpt_bin).at(cent_bin)->Fill(R_truth_truth, R_reco_reco, track_mc_pt, jet_weight);
                            h_dR_change.at(jetpt_bin).at(n_cent_bins-1)->Fill(R_truth_truth, R_reco_reco, track_mc_pt, jet_weight);

							if (R_reco_reco < _dR_max) reco_posRes_ChPS.at(dr_bin_reco_reco).at(cent_bin)->Fill(track_mc_pt, matched_truth_jet_pt, jet_weight);
							if (R_truth_truth < _dR_max) truth_posRes_ChPS.at(dr_bin_truth_truth).at(cent_bin)->Fill(track_mc_pt, matched_truth_jet_pt, jet_weight);
							if (R_reco_truth < _dR_max) rt_posRes_ChPS.at(dr_bin_reco_truth).at(cent_bin)->Fill(track_mc_pt, matched_truth_jet_pt, jet_weight);
							if (R_truth_reco < _dR_max) tr_posRes_ChPS.at(dr_bin_truth_reco).at(cent_bin)->Fill(track_mc_pt, matched_truth_jet_pt, jet_weight);

							//Only for truth matched tracks
							if (z > z_max) z_max = z;
							if (pt > pT_max) pT_max = pt;

						}
						isFake=false;
						h_dR_binning->Fill(R_reco_reco);
					}

				} //end truthlinked trk loop
				  //Fake/UE tracks
				if (isFake)
				{
					//TODO check weighting of UE
					if (R < _dR_max)
					{
						ff_UE_z.at(dr_bin).at(cent_bin)->Fill(z,jet_pt, jet_weight*eff_weight);
						ff_UE_pT.at(dr_bin).at(cent_bin)->Fill(pt,jet_pt, jet_weight*eff_weight);

                        ff_UE_z.at(dr_bin).at(n_cent_bins-1)->Fill(z,jet_pt, jet_weight*eff_weight);
                        ff_UE_pT.at(dr_bin).at(n_cent_bins-1)->Fill(pt,jet_pt, jet_weight*eff_weight);
                        
						int jetpt_bin = jetcorr->GetJetpTBin(jet_pt, (TAxis*)reco_posRes_ChPS.at(0).at(0)->GetYaxis());
						UE_distr.at(jetpt_bin).at(cent_bin)->Fill(R, pt, eta, jet_weight*eff_weight);
                        UE_distr.at(jetpt_bin).at(n_cent_bins-1)->Fill(R, pt, eta, jet_weight*eff_weight);
					}
				}
			}
		} // end reco track loop


		//JES plots
		if(_data_switch==1)
		{
			JES_v_max_pT.at(cent_bin)->Fill(pT_max,(jet_pt - truth_jet_pt_vector.at(TruthJetIndex.at(i)) )/truth_jet_pt_vector.at(TruthJetIndex.at(i)),truth_jet_pt_vector.at(TruthJetIndex.at(i)),jet_weight);
			JES_v_max_z.at(cent_bin)->Fill(z_max,(jet_pt - truth_jet_pt_vector.at(TruthJetIndex.at(i)) )/truth_jet_pt_vector.at(TruthJetIndex.at(i)),truth_jet_pt_vector.at(TruthJetIndex.at(i)),jet_weight);

			float truth_reco_R = DeltaR(jet_phi,jet_eta,truth_jet_phi_vector.at(TruthJetIndex.at(i)),truth_jet_eta_vector.at(TruthJetIndex.at(i)));

			h_jet_pos_v_zmax.at(cent_bin)->Fill(z_max,truth_reco_R,truth_jet_pt_vector.at(TruthJetIndex.at(i)),jet_weight);
			h_jet_pos_v_ptmax.at(cent_bin)->Fill(pT_max,truth_reco_R,truth_jet_pt_vector.at(TruthJetIndex.at(i)),jet_weight);

		}
		for (int nMultThreshold=0;nMultThreshold<trkcorr->nMultThresholds;nMultThreshold++) {h_jetpT_v_multiplicity.at(cent_bin)->Fill(jet_pt,nMultThreshold,trk_multiplicity[nMultThreshold]);}
	}// end reco jet loop



	std::vector<int> RecoJetIndex;

	//Loop over truth jets
	if(_data_switch == 1)
	{
		const xAOD::TruthParticleContainer * particles = 0;
		EL_RETURN_CHECK("execute",event->retrieve( particles, "TruthParticles"));

		RecoJetIndex = TruthMatching(truth_jet_pt_vector,truth_jet_eta_vector,truth_jet_phi_vector,
									 jet_pt_xcalib_vector,jet_eta_vector,jet_phi_vector,
									 _dR_truth_matching);

		for(unsigned int i=0; i<truth_jet_pt_vector.size(); i++)
		{
			float dR_reco_truth = -1.;

			if (RecoJetIndex.at(i) >= 0)
			{
				dR_reco_truth = DeltaR(truth_jet_phi_vector.at(i),truth_jet_eta_vector.at(i),jet_phi_vector.at(RecoJetIndex.at(i)),jet_eta_vector.at(RecoJetIndex.at(i)) );
			}



			float jet_weight = 1.;
			jet_weight *=event_weight;

			truth_jet_pt = truth_jet_pt_vector.at(i);
			truth_jet_eta = truth_jet_eta_vector.at(i);
			truth_jet_y = truth_jet_y_vector.at(i);
			truth_jet_phi = truth_jet_phi_vector.at(i);
			int y_bin = jetcorr->GetJetYBin(truth_jet_y);

			if(!truth_jet_isolated_vector.at(i)) continue;

			if (RecoJetIndex.at(i) >= 0)
			{
				h_jet_for_eff.at(cent_bin)->Fill(truth_jet_pt, fabs(truth_jet_eta), truth_jet_phi, jet_weight);
                h_jet_for_eff.at(n_cent_bins-1)->Fill(truth_jet_pt, fabs(truth_jet_eta), truth_jet_phi, jet_weight);
			}

			h_jet_for_eff_full.at(cent_bin)->Fill(truth_jet_pt, fabs(truth_jet_eta), truth_jet_phi, jet_weight);
            h_jet_for_eff_full.at(n_cent_bins-1)->Fill(truth_jet_pt, fabs(truth_jet_eta), truth_jet_phi, jet_weight);

			if (truth_jet_pt < _truthpTjetCut) continue;
			if (fabs(truth_jet_y)>(2.5 - _dR_max)) continue; //cut on rapidity (simultaniously with 2.5-drmax on pseudorapidity)

			//TODO do we want to reweigth also truth spectrum?
			//if (_applyReweighting) jet_weight*=jetcorr->GetJetReweightingFactor(truth_jet_pt,truth_jet_eta,cent_bin);

			h_true_jet_spectrum.at(y_bin).at(cent_bin)->Fill(truth_jet_pt, jet_weight);
			h_true_jet_spectrum.at(jetcorr->nJetYBins - 1).at(cent_bin)->Fill(truth_jet_pt, jet_weight);
			h_true_jet_spectrum.at(y_bin).at(n_cent_bins-1)->Fill(truth_jet_pt, jet_weight);
			h_true_jet_spectrum.at(jetcorr->nJetYBins - 1).at(n_cent_bins-1)->Fill(truth_jet_pt, jet_weight);

			xAOD::TruthParticleContainer::const_iterator truth_itr = particles->begin();
			xAOD::TruthParticleContainer::const_iterator truth_end = particles->end();

			double truth_z_max = 0.;
			double truth_pt_max = 0.;

			for( ; truth_itr!=truth_end; ++truth_itr)
			{
				if (_useCharge!=0 && ((int)(*truth_itr)->charge())!=_useCharge) continue;
				int ty=getTypeTruth((*truth_itr)->barcode(),(*truth_itr)->pdgId(),(*truth_itr)->status(),(*truth_itr)->charge());
				if(ty!=1 && ty!=5) continue;

				//get the tracks....
				float eta = (*truth_itr)->eta();
				float phi = (*truth_itr)->phi();
				float pt = (*truth_itr)->pt()/ 1000.0;

				if (fabs(eta)>=2.5) continue;
				if (pt<_pTtrkCut) continue;

				//Only tracks associated with a jet
				float R = DeltaR(phi,eta,truth_jet_phi,truth_jet_eta);
				if(R > _dR_max) continue;

				//Get z
				double truth_z = cos(R)*pt / truth_jet_pt;

				if (truth_z > truth_z_max) truth_z_max = truth_z;
				if (pt >  truth_pt_max ) truth_pt_max = pt;

				if (R < 0.4 && dR_reco_truth > 0.)
				{
					h_jet_pos_v_truth_pt.at(cent_bin)->Fill(pt,dR_reco_truth,truth_jet_pt_vector.at(i),jet_weight);
				}

				//alternative
				//TVector3 jet3vector; jet3vector.SetPtEtaPhi(truth_jet_pt,truth_jet_eta,truth_jet_phi);
				//TVector3 particle3vector; particle3vector.SetPtEtaPhi(pt,eta,phi);
				//double truth_z =  (jet3vector * particle3vector)/ pow(truth_jet_pt*cosh(truth_jet_eta),2);

				int dr_bin = trkcorr->GetdRBin(R);

				ff_truth.at(dr_bin).at(cent_bin)->Fill(truth_z,truth_jet_pt, jet_weight);
				ChPS_truth.at(dr_bin).at(cent_bin)->Fill(pt,truth_jet_pt, jet_weight);

				int deta_bin = trkcorr->GetdRBin(fabs(DeltaEta(eta, truth_jet_eta)));
				ChPS_truth_deta.at(deta_bin).at(cent_bin)->Fill(pt,truth_jet_pt, jet_weight);

				int dphi_bin = trkcorr->GetdRBin(fabs(DeltaPhi(phi,truth_jet_phi)));
				ChPS_truth_dphi.at(dphi_bin).at(cent_bin)->Fill(pt,truth_jet_pt, jet_weight);

				//Centrality inclusive
				ff_truth.at(dr_bin).at(n_cent_bins-1)->Fill(truth_z,truth_jet_pt, jet_weight);
				ChPS_truth.at(dr_bin).at(n_cent_bins-1)->Fill(pt,truth_jet_pt, jet_weight);
			}


			if (dR_reco_truth > 0.)
			{
				h_jet_pos_v_truth_zmax.at(cent_bin)->Fill(truth_z_max,dR_reco_truth,truth_jet_pt_vector.at(i),jet_weight);
				h_jet_pos_v_truth_ptmax.at(cent_bin)->Fill(truth_pt_max,dR_reco_truth,truth_jet_pt_vector.at(i),jet_weight);
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
	RecoJetIndex.clear();

	for (int j=0;j<_nTriggers;j++)
	{
		isTriggered[j].clear();
	}

	jet_pt_xcalib_vector.shrink_to_fit();
	jet_phi_vector.shrink_to_fit();
	jet_eta_vector.shrink_to_fit();

	truth_jet_indices.shrink_to_fit();
	truth_jet_pt_vector.shrink_to_fit();
	truth_jet_eta_vector.shrink_to_fit();
	truth_jet_y_vector.shrink_to_fit();
	truth_jet_phi_vector.shrink_to_fit();
	TruthJetIndex.shrink_to_fit();
	RecoJetIndex.shrink_to_fit();

	for (int j=0;j<_nTriggers;j++){
		isTriggered[j].shrink_to_fit();
	}

	return EL::StatusCode::SUCCESS;
}

EL::StatusCode PbPbFFShape :: postExecute (){
	return EL::StatusCode::SUCCESS;
}

EL::StatusCode PbPbFFShape :: finalize (){
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

EL::StatusCode PbPbFFShape :: histFinalize (){
	cout<<"Events = "<< m_eventCounter<<endl;
	return EL::StatusCode::SUCCESS;
}
