#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include "pPbFragmentation/TrackingPerformance.h"
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

ClassImp(TrackingPerformance)

#define EL_RETURN_CHECK( CONTEXT, EXP )			\
do {							\
if( ! EXP.isSuccess() ) {				\
Error( CONTEXT,					\
XAOD_MESSAGE( "Failed to execute: %s" ),	\
#EXP );					\
return EL::StatusCode::FAILURE;			\
}							\
} while( false )


TrackingPerformance :: TrackingPerformance ()
{
}

EL::StatusCode TrackingPerformance :: setupJob (EL::Job& job)
{
	// let's initialize the algorithm to use the xAODRootAccess package
	job.useXAOD ();

	EL_RETURN_CHECK( "setupJob()", xAOD::Init() ); // call before opening first file
	cout << " Job setup done!" << endl;

	nCentBins = GetCentralityNBins(_centrality_scheme);
	cout << "Number of centrality bins: " << nCentBins <<endl; //all bins + 1 inclusive

	std::cout << "[TrackingPerformance() : initialized with jet radius = " << _jet_radius << std::endl;

	return EL::StatusCode::SUCCESS;
}

EL::StatusCode TrackingPerformance :: histInitialize ()
{

	trkcorr = new TrackCorrector(_cut_level.c_str(),nCentBins-1,_eff_jety);
	trkcorr->drmax = _dR_max;
	trkcorr->InitdRBinRange();
	int _ndRBins = trkcorr->ndRBins + 1;


	cout << " Setting  histograms" << endl;


	int ptJetBinsN, etaJetBinsN, phiJetBinsN, ptTrkBinsN, etaTrkBinsN, phiTrkBinsN, d0z0BinsN, finehitsBinsN, MCProbBinsN=10, trk_resN=2000,  etaFineBinsN, respBinsN, fine_respBinsN;
	double ptJetBins[1000], etaJetBins[1000], phiJetBins[1000], ptTrkBins[1000], etaTrkBins[1000], phiTrkBins[1000], d0z0Bins[1000], finehitsBins[1000], MCProbBins[1000], trk_res[1000], etaFineBins[1000], respBins[1000], fine_respBins[1000];

	SetupBinning(0, "pt-jet-PbPb", ptJetBins, ptJetBinsN);
	SetupBinning(0, "eta-jet", etaJetBins, etaJetBinsN);
	SetupBinning(0, "phi-trk", phiJetBins, phiJetBinsN);
	SetupBinning(0, "pt-trk", ptTrkBins, ptTrkBinsN);
	SetupBinning(0, "eta-trk", etaTrkBins, etaTrkBinsN);
	SetupBinning(0, "phi-trk", phiTrkBins, phiTrkBinsN);
	SetupBinning(0, "d0z0", d0z0Bins, d0z0BinsN);
	SetupBinning(0, "hits_fine", finehitsBins, finehitsBinsN);
	SetupBinning(0, "MC-prob", MCProbBins, MCProbBinsN);
	SetupBinning(0, "trk_res", trk_res, trk_resN);
	SetupBinning(0, "eta-fine", etaFineBins, etaFineBinsN);
	SetupBinning(0, "resp", respBins, respBinsN);
	SetupBinning(0, "resp-fine", fine_respBins, fine_respBinsN);
	Double_t PVBins[3]={0,1,2};
	int PVBinsN=2;

	TH2D* temphist_2D = nullptr;
	TH3D* temphist_3D = nullptr;

	h_FCal_Et = new TH1D("h_FCal_Et",";FCal E_{T};N",100,0,5);
	h_FCal_Et->Sumw2();

	h_RejectionHisto = new TH1D("RejectionHisto","RejectionHisto",9,0,9);
	SetRejectionHistogram(h_RejectionHisto);

	h_centrality = new TH1D("Centrality","Centrality",10,0,10);
	h_centrality->Sumw2();

	h_trk_resolution = new TH3D("h_trk_resolution","h_trk_resolution",ptTrkBinsN, ptTrkBins,etaTrkBinsN,etaTrkBins,trk_resN,trk_res);
	h_trk_resolution->Sumw2();

	h_reco_trk_map = new TH3D("h_reco_trk_map","h_reco_trk_map;p_{T};#eta;#phi",ptTrkBinsN, ptTrkBins,etaTrkBinsN,etaTrkBins,phiTrkBinsN,phiTrkBins);
	h_reco_trk_map->Sumw2();

	h_reco_trk_map_nocuts = new TH3D("h_reco_trk_map_nocuts","h_reco_trk_map_nocuts;p_{T};#eta;#phi",ptTrkBinsN, ptTrkBins,etaTrkBinsN,etaTrkBins,phiTrkBinsN,phiTrkBins);
	h_reco_trk_map_nocuts->Sumw2();

	h_truth_trk_map = new TH3D("h_truth_trk_map","h_truth_trk_map;p_{T};#eta;#phi",ptTrkBinsN, ptTrkBins,etaTrkBinsN,etaTrkBins,phiTrkBinsN,phiTrkBins);
	h_truth_trk_map->Sumw2();

	h_dR_axis = new TH1D(Form("h_dR_axis"),Form("h_dR_axis"),_ndRBins-1,trkcorr->dRrange);


	wk()->addOutput (h_FCal_Et);
	wk()->addOutput (h_centrality);
	wk()->addOutput (h_RejectionHisto);
	wk()->addOutput (h_trk_resolution);
	wk()->addOutput (h_reco_trk_map);
	wk()->addOutput (h_reco_trk_map_nocuts);
	wk()->addOutput (h_truth_trk_map);
	wk()->addOutput (h_dR_axis);

	CutsOn.push_back("Reco"); //Must be first
	CutsOn.push_back("d0");
	CutsOn.push_back("z0sintheta");
	CutsOn.push_back("d0Sign");
	CutsOn.push_back("z0sinthetaSign");
	CutsOn.push_back("Z0");
	CutsOn.push_back("Z0SinTheta");
	CutsOn.push_back("InnermostLayersHits");
	CutsOn.push_back("SiHits");
	CutsOn.push_back("PixelHits");
	CutsOn.push_back("SctHits");
	CutsOn.push_back("TrtHits");
	CutsOn.push_back("FitQuality");
	CutsOn.push_back("AllCuts"); //Must be last


	NCuts = CutsOn.size();

	h_eff_dR_matched =  vector<vector<TH3D*> > (_ndRBins, vector<TH3D*>(nCentBins));
	h_eff_dR =  vector<vector<TH3D*> > (_ndRBins, vector<TH3D*>(nCentBins));
	h_eff_deta_matched =  vector<vector<TH3D*> > (_ndRBins, vector<TH3D*>(nCentBins));
	h_eff_deta =  vector<vector<TH3D*> > (_ndRBins, vector<TH3D*>(nCentBins));
	h_eff_dphi_matched =  vector<vector<TH3D*> > (_ndRBins, vector<TH3D*>(nCentBins));
	h_eff_dphi =  vector<vector<TH3D*> > (_ndRBins, vector<TH3D*>(nCentBins));


	for (int i=0;i<nCentBins;i++)
	{
		temphist_3D = new TH3D(Form("h_reco_Injet_matched_cent%i",i),Form("h_reco_Injet_matched_cent%i",i),ptJetBinsN,ptJetBins,ptTrkBinsN, ptTrkBins,etaFineBinsN,etaFineBins);
		h_reco_Injet_matched.push_back(temphist_3D);


		temphist_3D = new TH3D(Form("h_reco_jet_map_cent%i",i),Form("h_reco_Injet_matched_cent%i",i),ptJetBinsN, ptJetBins, etaJetBinsN, etaJetBins, phiJetBinsN, phiJetBins);
		h_reco_jet_map.push_back(temphist_3D);

		temphist_3D = new TH3D(Form("h_eff_matched_cent%i",i),Form("h_eff_matched_cent%i",i),ptJetBinsN,ptJetBins,ptTrkBinsN, ptTrkBins,etaFineBinsN,etaFineBins);
		h_eff_matched.push_back(temphist_3D);

		temphist_3D = new TH3D(Form("h_eff_total_cent%i",i),Form("h_eff_total_cent%i",i),ptJetBinsN,ptJetBins,ptTrkBinsN, ptTrkBins,etaFineBinsN,etaFineBins);
		h_eff_total.push_back(temphist_3D);


		temphist_3D = new TH3D(Form("h_eff_Injet_matched_cent%i",i),Form("h_eff_Injet_matched_cent%i",i),ptJetBinsN,ptJetBins,ptTrkBinsN, ptTrkBins,etaFineBinsN,etaFineBins);
		h_eff_Injet_matched.push_back(temphist_3D);

		temphist_3D = new TH3D(Form("h_eff_Injet_cent%i",i),Form("h_eff_Injet_cent%i",i),ptJetBinsN,ptJetBins,ptTrkBinsN, ptTrkBins,etaFineBinsN,etaFineBins);
		h_eff_Injet.push_back(temphist_3D);

		temphist_3D = new TH3D(Form("h_eff_Injet_entries_cent%i",i),Form("h_eff_Injet_entries_cent%i",i),ptJetBinsN,ptJetBins,ptTrkBinsN, ptTrkBins,etaFineBinsN,etaFineBins);
		h_eff_Injet_entries.push_back(temphist_3D);

		temphist_3D = new TH3D(Form("h_eff_entries_cent%i",i),Form("h_eff_entries_cent%i",i),ptJetBinsN,ptJetBins,ptTrkBinsN, ptTrkBins,etaFineBinsN,etaFineBins);
		h_eff_entries.push_back(temphist_3D);

		temphist_3D = new TH3D(Form("h_trk_foreff_full_cent%i",i),Form("h_trk_foreff_full_cent%i",i), phiTrkBinsN, phiTrkBins, ptTrkBinsN, ptTrkBins, etaTrkBinsN, etaTrkBins);
		h_trk_foreff_full.push_back(temphist_3D);

		temphist_3D = new TH3D(Form("h_trk_foreff_matched_cent%i",i),Form("h_trk_foreff_matched_cent%i",i), phiTrkBinsN, phiTrkBins, ptTrkBinsN, ptTrkBins, etaTrkBinsN, etaTrkBins);
		h_trk_foreff_matched.push_back(temphist_3D);

        temphist_3D = new TH3D(Form("h_trk_foreff_entries_cent%i",i),Form("h_trk_foreff_entries_cent%i",i), phiTrkBinsN, phiTrkBins, ptTrkBinsN, ptTrkBins, etaTrkBinsN, etaTrkBins);
        h_trk_foreff_entries.push_back(temphist_3D);
        
		temphist_3D = new TH3D(Form("h_fake_v_jet_cent%i",i),Form("h_fake_v_jet_cent%i",i),ptJetBinsN,ptJetBins,ptTrkBinsN, ptTrkBins,etaFineBinsN,etaFineBins);
		h_fake_v_jet.push_back(temphist_3D);

		temphist_3D = new TH3D(Form("h_fake_v_jet_PV_cent%i",i),Form("h_fake_v_jet_PV_cent%i",i),ptJetBinsN,ptJetBins,ptTrkBinsN, ptTrkBins,etaFineBinsN,etaFineBins);
		h_fake_v_jet_PV.push_back(temphist_3D);

		temphist_3D = new TH3D(Form("h_trk_scale_cent%i",i),Form("h_trk_scale_cent%i",i),ptTrkBinsN, ptTrkBins,respBinsN,respBins,etaTrkBinsN, etaTrkBins);
		h_trk_scale.push_back(temphist_3D);
		h_trk_scale.at(i)->Sumw2();

		temphist_3D = new TH3D(Form("h_trk_eta_scale_cent%i",i),Form("h_trk_eta_scale_cent%i",i),ptTrkBinsN, ptTrkBins,fine_respBinsN,fine_respBins,etaTrkBinsN, etaTrkBins);
		h_trk_eta_scale.push_back(temphist_3D);
		h_trk_eta_scale.at(i)->Sumw2();

		temphist_3D = new TH3D(Form("h_trk_phi_scale_cent%i",i),Form("h_trk_phi_scale_cent%i",i),ptTrkBinsN, ptTrkBins,fine_respBinsN,fine_respBins,etaTrkBinsN, etaTrkBins);
		h_trk_phi_scale.push_back(temphist_3D);
		h_trk_phi_scale.at(i)->Sumw2();

		temphist_3D = new TH3D(Form("h_trk_R_scale_cent%i",i),Form("h_trk_R_scale_cent%i",i),ptTrkBinsN, ptTrkBins,fine_respBinsN,fine_respBins,etaTrkBinsN, etaTrkBins);
		h_trk_R_scale.push_back(temphist_3D);
		h_trk_R_scale.at(i)->Sumw2();

		temphist_3D = new TH3D(Form("h_truth_outside_cent%i",i),Form("h_truth_outside_cent%i",i), ptTrkBinsN, ptTrkBins, etaTrkBinsN, etaTrkBins, phiTrkBinsN, phiTrkBins);
		h_truth_outside.push_back(temphist_3D);
		h_truth_outside.at(i)->Sumw2();

		temphist_3D = new TH3D(Form("h_reco_outside_cent%i",i),Form("h_reco_outside_cent%i",i), ptTrkBinsN, ptTrkBins, etaTrkBinsN, etaTrkBins, phiTrkBinsN, phiTrkBins);
		h_reco_outside.push_back(temphist_3D);
		h_reco_outside.at(i)->Sumw2();

		temphist_3D = new TH3D(Form("h_trk_eff_matched_map_cent%i",i),Form("h_trk_eff_matched_map_cent%i",i), ptTrkBinsN, ptTrkBins, etaTrkBinsN, etaTrkBins, phiTrkBinsN, phiTrkBins);
		h_trk_eff_matched_map.push_back(temphist_3D);
		h_trk_eff_matched_map.at(i)->Sumw2();

		temphist_3D = new TH3D(Form("h_trk_eff_map_cent%i",i),Form("h_trk_eff_map_cent%i",i), ptTrkBinsN, ptTrkBins, etaTrkBinsN, etaTrkBins, phiTrkBinsN, phiTrkBins);
		h_trk_eff_map.push_back(temphist_3D);
		h_trk_eff_map.at(i)->Sumw2();

		for (int j = 0; j < _ndRBins; j++)
		{
			temphist_3D = new TH3D(Form("h_eff_matched_dR%i_cent%i",j,i),Form("h_eff_matched_dR%i_cent%i",j,i),ptJetBinsN,ptJetBins,ptTrkBinsN, ptTrkBins,etaFineBinsN,etaFineBins);
			h_eff_dR_matched.at(j).at(i) = temphist_3D;

			temphist_3D = new TH3D(Form("h_eff_dR%i_cent%i",j,i),Form("h_eff_dR%i_cent%i",j,i),ptJetBinsN,ptJetBins,ptTrkBinsN, ptTrkBins,etaFineBinsN,etaFineBins);
			h_eff_dR.at(j).at(i) = temphist_3D;

			temphist_3D = new TH3D(Form("h_eff_matched_deta%i_cent%i",j,i),Form("h_eff_matched_deta%i_cent%i",j,i),ptJetBinsN,ptJetBins,ptTrkBinsN, ptTrkBins,etaFineBinsN,etaFineBins);
			h_eff_deta_matched.at(j).at(i) = temphist_3D;

			temphist_3D = new TH3D(Form("h_eff_deta%i_cent%i",j,i),Form("h_eff_deta%i_cent%i",j,i),ptJetBinsN,ptJetBins,ptTrkBinsN, ptTrkBins,etaFineBinsN,etaFineBins);
			h_eff_deta.at(j).at(i) = temphist_3D;

			temphist_3D = new TH3D(Form("h_eff_matched_dphi%i_cent%i",j,i),Form("h_eff_matched_dphi%i_cent%i",j,i),ptJetBinsN,ptJetBins,ptTrkBinsN, ptTrkBins,etaFineBinsN,etaFineBins);
			h_eff_dphi_matched.at(j).at(i) = temphist_3D;

			temphist_3D = new TH3D(Form("h_eff_dphi%i_cent%i",j,i),Form("h_eff_dphi%i_cent%i",j,i),ptJetBinsN,ptJetBins,ptTrkBinsN, ptTrkBins,etaFineBinsN,etaFineBins);
			h_eff_dphi.at(j).at(i) = temphist_3D;

			h_eff_dR_matched.at(j).at(i)->Sumw2();
			h_eff_dR.at(j).at(i)->Sumw2();
			h_eff_deta_matched.at(j).at(i)->Sumw2();
			h_eff_deta.at(j).at(i)->Sumw2();
			h_eff_dphi_matched.at(j).at(i)->Sumw2();
			h_eff_dphi.at(j).at(i)->Sumw2();

			wk()->addOutput (h_eff_dR_matched.at(j).at(i));
			wk()->addOutput (h_eff_dR.at(j).at(i));
			wk()->addOutput (h_eff_deta_matched.at(j).at(i));
			wk()->addOutput (h_eff_deta.at(j).at(i));
			wk()->addOutput (h_eff_dphi_matched.at(j).at(i));
			wk()->addOutput (h_eff_dphi.at(j).at(i));
		}



		for (int j = 0 ; j < NCuts ; j++)
		{
			temphist_3D = new TH3D(Form("h_cut_flow_%s_cent%i",CutsOn.at(j).c_str(),i),Form("h_cut_flow_%s_c%i",CutsOn.at(j).c_str(),i), ptTrkBinsN, ptTrkBins, etaTrkBinsN, etaTrkBins, phiTrkBinsN, phiJetBins);
			h_tmp.push_back(temphist_3D);
		}

		h_cut_flow_cent.push_back(h_tmp);
		h_tmp.clear();


		wk()->addOutput (h_reco_Injet_matched.at(i));
		wk()->addOutput (h_reco_jet_map.at(i));
		wk()->addOutput (h_eff_matched.at(i));
		wk()->addOutput (h_eff_total.at(i));

		wk()->addOutput (h_eff_Injet_matched.at(i));
		wk()->addOutput (h_eff_Injet.at(i));
		wk()->addOutput (h_trk_eff_matched_map.at(i));
		wk()->addOutput (h_trk_eff_map.at(i));
		wk()->addOutput (h_eff_Injet_entries.at(i));
		wk()->addOutput (h_eff_entries.at(i));
        wk()->addOutput (h_trk_foreff_entries.at(i));
		wk()->addOutput (h_trk_foreff_matched.at(i));
		wk()->addOutput (h_trk_foreff_full.at(i));
		wk()->addOutput (h_fake_v_jet.at(i));
		wk()->addOutput (h_fake_v_jet_PV.at(i));
		wk()->addOutput (h_trk_scale.at(i));
		wk()->addOutput (h_trk_eta_scale.at(i));
		wk()->addOutput (h_trk_phi_scale.at(i));
		wk()->addOutput (h_trk_R_scale.at(i));
		wk()->addOutput (h_truth_outside.at(i));
		wk()->addOutput (h_reco_outside.at(i));

		for (int j = 0 ; j < NCuts ; j++)
		{
			wk()->addOutput (h_cut_flow_cent.at(i).at(j));
		}


	}


	h_d0_nocut = new TH2D("h_d0_nocut","h_d0_nocut", ptTrkBinsN, ptTrkBins, d0z0BinsN, d0z0Bins);
	h_d0_cut = new TH2D("h_d0_cut","h_d0_cut", ptTrkBinsN, ptTrkBins, d0z0BinsN, d0z0Bins);

	h_d0sign_wrt_truth = new TH3D("h_d0sign_wrt_truth","h_d0sign_wrt_truth", ptTrkBinsN, ptTrkBins, etaTrkBinsN, etaTrkBins, phiTrkBinsN, phiJetBins);
	h_z0sign_wrt_truth = new TH3D("h_z0sign_wrt_truth","h_z0sign_wrt_truth", ptTrkBinsN, ptTrkBins, etaTrkBinsN, etaTrkBins, phiTrkBinsN, phiJetBins);

	h_d0sign_subset = new TH3D("h_d0sign_subset","h_d0sign_subset", ptTrkBinsN, ptTrkBins, etaTrkBinsN, etaTrkBins, phiTrkBinsN, phiJetBins);
	h_z0sign_subset = new TH3D("h_z0sign_subset","h_z0sign_subset", ptTrkBinsN, ptTrkBins, etaTrkBinsN, etaTrkBins, phiTrkBinsN, phiJetBins);

	wk()->addOutput (h_d0sign_wrt_truth);
	wk()->addOutput (h_z0sign_wrt_truth);
	wk()->addOutput (h_d0sign_subset);
	wk()->addOutput (h_z0sign_subset);

	wk()->addOutput (h_d0_cut);
	wk()->addOutput (h_d0_nocut);


	cout << " Histograms  ready" << endl;
	cout << " tree  ready" << endl;

	cout << "Parametrization of d0 cut" << endl;
	f_d0_cut = new TF1("f1", "[0]*exp([1]*x)+[2]*exp([3]*x)", 0.4, 500);
	f_d0_cut->SetParameters(0.472367, -0.149934, 0.193095, 0.000337765);

	return EL::StatusCode::SUCCESS;
}

EL::StatusCode TrackingPerformance :: fileExecute ()
{
	// Here you do everything that needs to be done exactly once for every
	// single file, e.g. collect a list of all lumi-blocks processed
	return EL::StatusCode::SUCCESS;
}



EL::StatusCode TrackingPerformance :: changeInput (bool firstFile)
{
	return EL::StatusCode::SUCCESS;
}



EL::StatusCode TrackingPerformance :: initialize ()
{
	// count number of events
	cout << " Starting initialization" << endl;
	m_eventCounter = 0;

	xAOD::TEvent* event = wk()->xaodEvent();

	// as a check, let's see the number of events in our xAOD
	Info("initialize()", "Number of events = %lli", event->getEntries() ); // print long long int

	//Track Selector Tool
	m_trackSelectorTool = new InDet::InDetTrackSelectionTool("InDetTrackSelectorTool");
	SetCutLevel(m_trackSelectorTool, _cut_level.c_str());
	EL_RETURN_CHECK("initialize()",m_trackSelectorTool->initialize());

	//Calibration
	const std::string name = "TrackingPerformance"; //string describing the current thread, for logging
	TString jetAlgo = "AntiKt4HI"; //String describing your jet collection, for example AntiKt4EMTopo or AntiKt4LCTopo (see below)
	TString config = "JES_MC15c_HI_Nov2016.config"; //Path to global config used to initialize the tool (see below)
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

	//Jet corrector
	jetcorr = new JetCorrector();
	jetcorr->is_pp = (_dataset==3); //dataset 3 is 5 TeV pp MC
	cout << " Initialization done" << endl;
	return EL::StatusCode::SUCCESS;
}

//Loop over events
EL::StatusCode TrackingPerformance :: execute (){

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


//	if (m_eventCounter != 591) return EL::StatusCode::SUCCESS;
	//cout << "Event: " << m_eventCounter << endl;


	//All events
	h_RejectionHisto->Fill(0.5);

	//---------------------------
	//     Event information
	//---------------------------
/*
	const xAOD::EventInfo* ei = nullptr;
	EL_RETURN_CHECK("execute",event->retrieve( ei, "EventInfo"));

	//create pointer to an EventInfo object which may or may not be the “main” one
	const xAOD::EventInfo* ei_mc=nullptr;
	if(ei->eventType( xAOD::EventInfo::IS_SIMULATION )) ei_mc=ei; //If the EventInfo knows its simulation, then it should already have weights etc.
	else
	{
		//loop on the subEvents and find the first one that is simulation and grab a pointer to it
		for(const auto& subEvent : ei->subEvents())
		{
			const xAOD::EventInfo* se = subEvent.ptr();
			if( ! se ) continue;
			if(se->eventType( xAOD::EventInfo::IS_SIMULATION ))
			{
				ei_mc=se;
			 break;
			}
		}
	}

	cout << Form("evntwegight new: %f", ei_mc->mcEventWeights().at(0)) << endl;
*/

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
	double FCalEt = 0;
	int cent_bin = 0;
	int cent_bin_fine = 0;

	double event_weight_fcal = 1;
	if (_centrality_scheme>1)
	{
		//Centrality
		const xAOD::HIEventShapeContainer* calos=0;
		EL_RETURN_CHECK("execute",event->retrieve( calos, "CaloSums"));
		FCalEt=calos->at(5)->et()*1e-6;
		cent_bin = GetCentralityBin(_centrality_scheme, FCalEt, isHIJING);
		cent_bin_fine = GetCentralityBin(_centrality_scheme, FCalEt,  isHIJING ); //Need for some tools
//		event_weight_fcal = jetcorr->GetFCalWeight(FCalEt);

//		if (isMC && isHIJING) event_weight_fcal = 1;
//		cout << event_weight_fcal << endl;
		h_centrality->Fill(cent_bin,event_weight_fcal);
		h_centrality->Fill(nCentBins-1,event_weight_fcal);
}

	h_FCal_Et->Fill(FCalEt, event_weight_fcal); //filled here to get proper event weight

	if (cent_bin < 0)
	{
		h_RejectionHisto->Fill(1.5);
		return EL::StatusCode::SUCCESS;
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
		if (fabs(primaryVertex->z())>150.)
		{
			h_RejectionHisto->Fill(5.5);
			return EL::StatusCode::SUCCESS;
		}
	}
	else
	{
		h_RejectionHisto->Fill(5.5);
		return EL::StatusCode::SUCCESS;
	}

	//Tracks
	const xAOD::TrackParticleContainer* recoTracks = 0;
	EL_RETURN_CHECK("execute",event->retrieve( recoTracks, "InDetTrackParticles"));


	//Jet vectors
	vector<float> jet_pt_xcalib_vector,jet_phi_vector,jet_eta_vector, jet_m_vector;
	vector<bool> jet_isolated_vector, truth_jet_isolated_vector;
	vector<float> truth_jet_eta_vector,truth_jet_phi_vector,truth_jet_pt_vector;
	vector<int> TruthJetIndex;

	TLorentzVector jet4vector;

	float event_weight = 1;
	double max_pt = 1;

	//***** Truth jets *****
	const xAOD::JetContainer * truth_jet = 0;
	EL_RETURN_CHECK("execute()",event->retrieve( truth_jet, _truth_jet_collection.c_str() ));
	xAOD::JetContainer::const_iterator truth_jet_itr = truth_jet->begin();
	xAOD::JetContainer::const_iterator truth_jet_end = truth_jet->end();
	for( ; truth_jet_itr != truth_jet_end; ++truth_jet_itr )
	{
		xAOD::JetFourMom_t jet_truth_4mom = (*truth_jet_itr)->jetP4();

		double pt    = (jet_truth_4mom.pt() * 0.001 );
		double eta    = (jet_truth_4mom.eta());
		double phi    = (jet_truth_4mom.phi());

		if (pt>max_pt) 	//event weight from leading truth jet
		{
			event_weight = jetcorr->GetJetWeight(pt, eta, phi);
			max_pt = pt;
			if (isMC && isHIJING) event_weight = 1;
		}

		if (pt < _truthpTjetCut) continue;
		if (fabs(eta)>(2.5 - _dR_max)) continue;

		//filling truth pt/eta/phi vectors
		truth_jet_pt_vector.push_back(pt);
		truth_jet_phi_vector.push_back(phi);
		truth_jet_eta_vector.push_back(eta);
	}
	truth_jet_isolated_vector = MTCorrector::GetIsolation(truth_jet_pt_vector,truth_jet_eta_vector,truth_jet_phi_vector, 1.0, _pt_iso); // -1 is for pT_iso == jet pT

	event_weight = event_weight*event_weight_fcal; //event weight is only set if MC. Otherwise default is 1.

	//***** Reco jets*****
	xAOD::TStore *store = new xAOD::TStore; //For calibration
	const xAOD::JetContainer* jets = 0;
	EL_RETURN_CHECK("execute()",event->retrieve( jets, _reco_jet_collection.c_str() ));

	xAOD::JetContainer* updatedjets = new xAOD::JetContainer();
	xAOD::AuxContainerBase* updatedjetsAux = new xAOD::AuxContainerBase();
	updatedjets->setStore( updatedjetsAux );

	xAOD::JetContainer::const_iterator jet_itr = jets->begin();
	xAOD::JetContainer::const_iterator jet_end = jets->end();

	EL_RETURN_CHECK("execute()",store->record( updatedjets, "updatedjets" ));
	EL_RETURN_CHECK("execute()",store->record( updatedjetsAux, "updatedjetsAux" ));

	//In PbPb MC official derivations, the default jet kinematics in DFAntiKt4HI are already calibrated. Nothing more is needed
	//in pp MC, there is no DF container. The AntiKt4HIJets container has unsubtracted, subtracted = EMScale, and default (incorrectly calibrated). To calibrate these, Set Consitutuent scale to usubtracted and calibrate


	for( ; jet_itr != jet_end; ++jet_itr )
	{
		xAOD::Jet newjet;// = new xAOD::Jet();
		newjet.makePrivateStore( **jet_itr );

		const xAOD::JetFourMom_t jet_4mom_def = newjet.jetP4();
		const xAOD::JetFourMom_t jet_4mom_subtr = newjet.jetP4("JetSubtractedScaleMomentum");
		const xAOD::JetFourMom_t jet_4mom_unsubtr = newjet.jetP4("JetUnsubtractedScaleMomentum");


		//if PbPb MC, no etajes calibration needed use defaults
		if (_dataset == 4) //do nothing

		//if pp MC, default is incorrect. Set constit to unsubtracted scale, and calibrate
		if (_dataset == 3)
		{
			newjet.setJetP4("JetConstitScaleMomentum",jet_4mom_unsubtr);
			EL_RETURN_CHECK("execute()", m_jetCalibration->applyCalibration( newjet ) );
		}

		xAOD::JetFourMom_t jet_4mom_calib = newjet.jetP4();
		float jet_pt  = (jet_4mom_calib.pt() * 0.001);
		float jet_eta  = (jet_4mom_calib.eta());
		float jet_phi  = (jet_4mom_calib.phi());
		float jet_m   = (newjet.m()*0.001);

		if (fabs(jet_eta)>2.5-_dR_max) continue;
		if (jet_pt < _pTjetCut) continue;

		jet_pt_xcalib_vector.push_back(jet_pt);
		jet_phi_vector.push_back(jet_phi);
		jet_eta_vector.push_back(jet_eta);
		jet_m_vector.push_back(jet_m);
	}


	store->clear();

	delete store;

	jet_isolated_vector = MTCorrector::GetIsolation(jet_pt_xcalib_vector, jet_eta_vector,jet_phi_vector, 1.0, _pt_iso); // -1 is for pT_iso == jet pT

	TruthJetIndex = TruthMatching(jet_pt_xcalib_vector,jet_eta_vector,jet_phi_vector,
								  truth_jet_pt_vector,truth_jet_eta_vector,truth_jet_phi_vector,
								  _dR_truth_matching);


	//Get Truth tracks
	const xAOD::TruthParticleContainer * particles = 0;
	EL_RETURN_CHECK("execute()",event->retrieve( particles, "TruthParticles" ));

	xAOD::TruthParticleContainer::const_iterator truth_itr = particles->begin();
	xAOD::TruthParticleContainer::const_iterator truth_end = particles->end();

	//Track vectors
	vector<float> trk_good_pt, trk_good_eta, trk_good_phi;
	vector<int>   trk_good_isMatched, trk_good_vertextype;
	vector<float> trk_good_matched_eta, trk_good_matched_phi, trk_good_matched_pt;
	vector<float> trk_good_matched_eta_reco, trk_good_matched_phi_reco, trk_good_matched_pt_reco;

	//Loop over reconstructed tracks
	for (const auto& trk : *recoTracks)
    {
		if (_useCharge!=0 && ((int)trk->charge())!=_useCharge) continue;

		//get the tracks....
		float pt = trk->pt()/1000.;
		float eta = trk->eta();
		float phi = trk->phi();
		if (fabs(eta)>=2.5) continue;
		if (pt < _pTtrkCut) continue;

		//Tracks have to be in jet
		bool good_jet = false;

		for(unsigned int i=0; i<jet_pt_xcalib_vector.size(); i++)
		{
			jet_eta = jet_eta_vector.at(i);
			jet_phi = jet_phi_vector.at(i);
			jet_pt = jet_pt_xcalib_vector.at(i);
			jet_m   = jet_m_vector.at(i);

			float R = DeltaR(phi,eta,jet_phi,jet_eta);

			int truthindex=TruthJetIndex.at(i);
			if (truthindex<0) continue; //Reco jet must be matched to truth jet
			if (fabs(jet_y) > (2.5 - _dR_max)) continue; //cut on rapidity (simultaniously with 2.5-drmax on pseudorapidity)
			if (R > _dR_max) continue; //track must be in jet
			
			good_jet = true;
		}

		double d0 = trk->d0();
		double theta = trk->theta();
		double z0pv=(trk->z0()+trk->vz()-(*vtx_itr)->z())*sin(theta);	// pp: trk->z0() - w.r.t. BS

		float ed0=sqrt(trk->definingParametersCovMatrix()(0,0));
		float ez0pv=sqrt( (trk->definingParametersCovMatrix()(1,1))*pow(sin(theta),2) +
						 (trk->definingParametersCovMatrix()(3,3))*pow(z0pv*cos(theta),2)+
						 2*(trk->definingParametersCovMatrix()(1,3))*fabs(sin(theta)*z0pv*cos(theta)) );
		if(ed0==0.0) ed0=1e-10;
		if(ez0pv==0.0) ez0pv=1e-10;

		int nShPixH= trk->auxdata< unsigned char >("numberOfPixelSharedHits");
		int nShSCTH = trk->auxdata< unsigned char >("numberOfSCTSharedHits");

		float sigd0 = fabs(d0/ed0);
		float sigz0 = fabs(z0pv/ez0pv);


		//all tracks spectra
		if (good_jet) h_reco_trk_map_nocuts->Fill(pt,eta,phi, event_weight);

		//Before cuts
		if (good_jet) h_cut_flow_cent.at(cent_bin).at(0)->Fill(pt, eta, phi, event_weight);
		if (good_jet) h_cut_flow_cent.at(nCentBins-1).at(0)->Fill(pt, eta, phi, event_weight);

		//Vertex association
		ElementLink< xAOD::VertexContainer > vtxLink = trk->auxdata<ElementLink< xAOD::VertexContainer> >("vertexLink");
		int vertex_type=0;
		if(vtxLink.isValid()) vertex_type =  (*vtxLink)->vertexType();

		//Truth matching
		bool isTruthMatched=false;
		ElementLink< xAOD::TruthParticleContainer > truthLink = trk->auxdata<ElementLink< xAOD::TruthParticleContainer > >("truthParticleLink");
		float mcprob = trk->auxdata<float>("truthMatchProbability");

		if(truthLink.isValid() && mcprob > _mcProbCut)
		{
			int trktype = getTypeReco((*truthLink)->barcode(),(*truthLink)->pdgId(),(*truthLink)->status(),(*truthLink)->charge(),mcprob,_mcProbCut);
			if (( trktype==1 || trktype==5 )) isTruthMatched=true;
		}

		//0		CutsOn.push_back("Reco"); //Must be first
		//1		CutsOn.push_back("d0");
		//2		CutsOn.push_back("z0sintheta");
		//3		CutsOn.push_back("d0Sign");
		//4		CutsOn.push_back("z0sinthetaSign");
		//5		CutsOn.push_back("Z0");
		//6		CutsOn.push_back("Z0SinTheta");
		//7		CutsOn.push_back("InnermostLayersHits");
		//8		CutsOn.push_back("SiHits");
		//9		CutsOn.push_back("PixelHits");
		//10	CutsOn.push_back("SctHits");
		//11	CutsOn.push_back("TrtHits");
		//12	CutsOn.push_back("FitQuality");
		//13	CutsOn.push_back("AllCuts"); //Must be last

		if (good_jet) h_d0_nocut->Fill(pt, d0);

		//d0 (from total)
		double d0_cut;
		if (!_cut_level.std::string::compare("NoCuts")) d0_cut = 999.0;
		else d0_cut = f_d0_cut->Eval(pt);

		if (fabs(d0) < d0_cut && good_jet)
		{
			h_cut_flow_cent.at(cent_bin).at(1)->Fill(pt, eta, phi, event_weight);
			h_cut_flow_cent.at(nCentBins-1).at(1)->Fill(pt, eta, phi, event_weight);
			h_d0_cut->Fill(pt, d0);
		}

		//z0sintheta (from total)
		double z0sintheta_cut = 1.0;
		if (fabs(z0pv) < z0sintheta_cut && good_jet)
		{
			h_cut_flow_cent.at(cent_bin).at(2)->Fill(pt, eta, phi, event_weight);
			h_cut_flow_cent.at(nCentBins-1).at(2)->Fill(pt, eta, phi, event_weight);
			//			h_z0pvs_tmp->Fill(z0pv);
		}

		// d0/z0sintheta sign.
		double sig_cut = 3.;
		if (sigd0 < sig_cut && good_jet)
		{
			h_cut_flow_cent.at(cent_bin).at(3)->Fill(pt, eta, phi, event_weight);
			h_cut_flow_cent.at(nCentBins-1).at(3)->Fill(pt, eta, phi, event_weight);
		}
		if (sigz0 < sig_cut && good_jet)
		{
			h_cut_flow_cent.at(cent_bin).at(4)->Fill(pt, eta, phi, event_weight);
			h_cut_flow_cent.at(nCentBins-1).at(4)->Fill(pt, eta, phi, event_weight);
		}

		if (sigd0 < sig_cut && isTruthMatched && good_jet) h_d0sign_wrt_truth->Fill(pt, eta, phi, event_weight);
		if (sigz0 < sig_cut && isTruthMatched && good_jet) h_z0sign_wrt_truth->Fill(pt, eta, phi, event_weight);

		//Track selector tool cut
		bool skip_track = false;
		if(!m_trackSelectorTool->accept(*trk, *vtx_itr )) skip_track = true;

		// Z0
		if (good_jet)
		{
			if ( m_trackSelectorTool->getTAccept().getCutResult("Z0") )
			{
				h_cut_flow_cent.at(cent_bin).at(5)->Fill(pt, eta, phi, event_weight);
				h_cut_flow_cent.at(nCentBins-1).at(5)->Fill(pt, eta, phi, event_weight);
			}

			// Z0SinTheta
			if ( m_trackSelectorTool->getTAccept().getCutResult("Z0SinTheta") )
			{
				h_cut_flow_cent.at(cent_bin).at(6)->Fill(pt, eta, phi, event_weight);
				h_cut_flow_cent.at(nCentBins-1).at(6)->Fill(pt, eta, phi, event_weight);
			}

			// InnermostLayersHits
			if ( m_trackSelectorTool->getTAccept().getCutResult("InnermostLayersHits") )
			{
				h_cut_flow_cent.at(cent_bin).at(7)->Fill(pt, eta, phi, event_weight);
				h_cut_flow_cent.at(nCentBins-1).at(7)->Fill(pt, eta, phi, event_weight);
			}

			// Silicon hits
			if ( m_trackSelectorTool->getTAccept().getCutResult("SiHits") )
			{
				h_cut_flow_cent.at(cent_bin).at(8)->Fill(pt, eta, phi, event_weight);
				h_cut_flow_cent.at(nCentBins-1).at(8)->Fill(pt, eta, phi, event_weight);
			}

			// Pixel hits
			if ( m_trackSelectorTool->getTAccept().getCutResult("PixelHits") )
			{
				h_cut_flow_cent.at(cent_bin).at(9)->Fill(pt, eta, phi, event_weight);
				h_cut_flow_cent.at(nCentBins-1).at(9)->Fill(pt, eta, phi, event_weight);
			}

			// SCT hits
			if ( m_trackSelectorTool->getTAccept().getCutResult("SctHits") )
			{
				h_cut_flow_cent.at(cent_bin).at(10)->Fill(pt, eta, phi, event_weight);
				h_cut_flow_cent.at(nCentBins-1).at(10)->Fill(pt, eta, phi, event_weight);
			}

			// TRT hits
			if ( m_trackSelectorTool->getTAccept().getCutResult("TRTHits") )
			{
				h_cut_flow_cent.at(cent_bin).at(11)->Fill(pt, eta, phi, event_weight);
				h_cut_flow_cent.at(nCentBins-1).at(11)->Fill(pt, eta, phi, event_weight);
			}

			// Fit quality
			if ( m_trackSelectorTool->getTAccept().getCutResult("FitQuality") )
			{
				h_cut_flow_cent.at(cent_bin).at(12)->Fill(pt, eta, phi, event_weight);
				h_cut_flow_cent.at(nCentBins-1).at(12)->Fill(pt, eta, phi, event_weight);
			}
			
		}

		//d0 cut
		if (skip_track) continue;
		if(fabs(d0) > d0_cut) continue; //pT dependent d0

		if (sigd0 < sig_cut && good_jet) h_d0sign_subset->Fill(pt, eta, phi, event_weight);
		if (sigz0 < sig_cut && good_jet) h_z0sign_subset->Fill(pt, eta, phi, event_weight);

		//determining track resolution
		if(isTruthMatched)
		{
			float reco_pt = pt;
			float reco_eta = eta;
			float reco_phi = phi;

			float truth_eta = ((*truthLink)->eta());
			float truth_pt = ((*truthLink)->pt())/1000.;
			float truth_phi = ((*truthLink)->phi());

			float dR = DeltaR(reco_phi,reco_eta,truth_phi,truth_eta);

			h_trk_scale.at(cent_bin)->Fill(truth_pt,(truth_pt/reco_pt) - 1,truth_eta, event_weight);
			h_trk_scale.at(nCentBins-1)->Fill(truth_pt,(truth_pt/reco_pt) - 1,truth_eta, event_weight);

			h_trk_eta_scale.at(cent_bin)->Fill(truth_pt,(reco_eta-truth_eta),truth_eta, event_weight);
			h_trk_eta_scale.at(nCentBins-1)->Fill(truth_pt,(reco_eta-truth_eta),truth_eta, event_weight);

			h_trk_phi_scale.at(cent_bin)->Fill(truth_pt,(reco_phi-truth_phi),truth_eta, event_weight);
			h_trk_phi_scale.at(nCentBins-1)->Fill(truth_pt,(reco_phi-truth_phi),truth_eta, event_weight);

			h_trk_R_scale.at(cent_bin)->Fill(truth_pt, dR, truth_eta, event_weight);
			h_trk_R_scale.at(nCentBins-1)->Fill(truth_pt, dR, truth_eta, event_weight);
		}


		if (good_jet)
		{
			h_reco_trk_map->Fill(pt,eta,phi, event_weight);
			h_cut_flow_cent.at(cent_bin).at(NCuts-1)->Fill(pt, eta, phi, event_weight);
			h_cut_flow_cent.at(nCentBins-1).at(NCuts-1)->Fill(pt, eta, phi, event_weight);
		}

		trk_good_eta.push_back( eta );
		trk_good_phi.push_back( phi );
		trk_good_pt.push_back( pt );
		trk_good_vertextype.push_back(vertex_type);

		if(isTruthMatched)
		{
			trk_good_matched_eta.push_back((*truthLink)->eta());
			trk_good_matched_pt.push_back((*truthLink)->pt()*0.001);
			trk_good_matched_phi.push_back((*truthLink)->phi());

			trk_good_matched_pt_reco.push_back(pt);
			trk_good_matched_eta_reco.push_back(eta);
			trk_good_matched_phi_reco.push_back(phi);

			trk_good_isMatched.push_back(1);
//            cout << Form("event:%i	pt:	%f->%f	,eta: %f->%f	, phi:%f->%f", m_eventCounter, pt, (*truthLink)->pt()*0.001, eta, (*truthLink)->eta(), phi, (*truthLink)->phi()) << endl;
            
            
			float MPS = (pt - (*truthLink)->pt()*0.001)/((*truthLink)->pt()*0.001);
			h_trk_resolution->Fill((*truthLink)->pt()*0.001,(*truthLink)->eta(),MPS, event_weight);
		}
		else
		{
			trk_good_isMatched.push_back(0);
		}
	}

	//Reco jets
	bool isFirstPass = true;

	for(unsigned int i=0; i<jet_pt_xcalib_vector.size(); i++)
	{
		jet_pt = jet_pt_xcalib_vector.at(i);
		jet_eta = jet_eta_vector.at(i);
		jet_phi = jet_phi_vector.at(i);
		jet_m   = jet_m_vector.at(i);
		if (jet_pt < 63.1)  continue;
		jet4vector.SetPtEtaPhiM(jet_pt, jet_eta, jet_phi, jet_m);
		double rapidity = jet4vector.Rapidity();
		double jet_y = rapidity;

		int truthindex=TruthJetIndex.at(i);
		if (truthindex<0) continue; //Reco jet must be matched to truth jet
		if (jet_pt < _pTjetCut)  continue;
		if (fabs(jet_y) > (2.5 - _dR_max)) continue; //cut on rapidity (simultaniously with 2.5-drmax on pseudorapidity) //removed 12/8 for tracking efficienciy

		h_reco_jet_map.at(cent_bin)->Fill(jet_pt, jet_eta, jet_phi, event_weight);
        h_reco_jet_map.at(nCentBins-1)->Fill(jet_pt, jet_eta, jet_phi, event_weight);

		for(unsigned int j = 0; j< trk_good_pt.size(); j++)
		{
			int tru_itr = 0;
			float eta = trk_good_eta.at(j);
			float phi = trk_good_phi.at(j);
			float pt = trk_good_pt.at(j);

			float R = DeltaR(phi,eta,jet_phi,jet_eta);
			if (R > _dR_max) continue;

			if (trk_good_isMatched.at(j)==1)
			{
				float eta_truth = trk_good_matched_eta.at(tru_itr);
				float phi_truth = trk_good_matched_phi.at(tru_itr);
				float pt_truth = trk_good_matched_pt.at(tru_itr);

				float R_truth = DeltaR(phi_truth,eta_truth,jet_phi,jet_eta);
				if (R_truth > _dR_max)
                {
                    h_truth_outside.at(cent_bin)->Fill(pt_truth, eta_truth, phi_truth, event_weight);
                    h_truth_outside.at(nCentBins-1)->Fill(pt_truth, eta_truth, phi_truth, event_weight);
                }

				tru_itr++;
			}
		}

		for(unsigned int j = 0; j< trk_good_matched_pt.size(); j++)
		{
			int tru_itr = 0;
			float pt = trk_good_matched_pt.at(j);
			float phi = trk_good_matched_phi.at(j);
			float eta = trk_good_matched_eta.at(j);

			float R = DeltaR(phi,eta,jet_phi,jet_eta);
			if (R > _dR_max) continue;

			float eta_reco = trk_good_matched_eta_reco.at(j);
			float phi_reco = trk_good_matched_phi_reco.at(j);
			float pt_reco = trk_good_matched_pt_reco.at(j);

			float R_reco = DeltaR(phi_reco,eta_reco,jet_phi,jet_eta);
			if (R_reco > _dR_max)
            {
                h_reco_outside.at(cent_bin)->Fill(pt_reco, eta_reco, phi_reco, event_weight);
                h_reco_outside.at(nCentBins-1)->Fill(pt_reco, eta_reco, phi_reco, event_weight);
            }
		}

		//All good tracks
		for(unsigned int j = 0; j< trk_good_pt.size(); j++)
		{
			if (trk_good_isMatched.at(j)==1) continue;
			float eta = trk_good_eta.at(j);
			float phi = trk_good_phi.at(j);
			float pt = trk_good_pt.at(j);

			//track-to-jet balance cut 3 sigma cut
			float R = DeltaR(phi,eta,jet_phi,jet_eta);
			if (R > _dR_max) continue;

			h_fake_v_jet[cent_bin]->Fill(jet_pt,pt,eta, event_weight);
			h_fake_v_jet[nCentBins-1]->Fill(jet_pt,pt,eta, event_weight);
			if (trk_good_vertextype.at(j)==1)
			{
				h_fake_v_jet_PV[cent_bin]->Fill(jet_pt,pt,eta, event_weight);
				h_fake_v_jet_PV[nCentBins-1]->Fill(jet_pt,pt,eta, event_weight);
			}
		}

		//Only those truth with matched reco tracks
		for(unsigned int j = 0; j< trk_good_matched_pt.size(); j++)
		{
			float eta = trk_good_matched_eta.at(j);
			float phi = trk_good_matched_phi.at(j);
			float pt = trk_good_matched_pt.at(j);
			float pt_reco = trk_good_matched_pt_reco.at(j);

			//Only tracks associated with a jet
			float R = DeltaR(phi,eta,jet_phi,jet_eta);
			if(R > _dR_max) continue;

			int deta_bin = trkcorr->GetdRBin(fabs(DeltaEta(eta,jet_eta)));
			h_eff_deta_matched.at(deta_bin).at(cent_bin)->Fill(jet_pt,pt,fabs(rapidity), event_weight);
            h_eff_deta_matched.at(deta_bin).at(nCentBins-1)->Fill(jet_pt,pt,fabs(rapidity), event_weight);

			int dphi_bin = trkcorr->GetdRBin(fabs(DeltaPhi(phi,jet_phi)));
			h_eff_dphi_matched.at(dphi_bin).at(cent_bin)->Fill(jet_pt,pt,fabs(rapidity), event_weight);
            h_eff_dphi_matched.at(dphi_bin).at(nCentBins-1)->Fill(jet_pt,pt,fabs(rapidity), event_weight);

			int dR_bin = trkcorr->GetdRBin(R);
			h_eff_dR_matched.at(dR_bin).at(cent_bin)->Fill(jet_pt,pt,fabs(rapidity), event_weight);
            h_eff_dR_matched.at(dR_bin).at(nCentBins-1)->Fill(jet_pt,pt,fabs(rapidity), event_weight);

			h_eff_matched.at(cent_bin)->Fill(jet_pt, pt, eta, event_weight);
            h_eff_matched.at(nCentBins-1)->Fill(jet_pt, pt, eta, event_weight);

            h_trk_eff_matched_map.at(cent_bin)->Fill(pt, eta, phi, event_weight);
            h_trk_eff_matched_map.at(nCentBins-1)->Fill(pt, eta, phi, event_weight);

			h_reco_Injet_matched[cent_bin]->Fill(jet_pt, pt_reco, eta, event_weight);
			h_reco_Injet_matched[nCentBins-1]->Fill(jet_pt,pt_reco,eta, event_weight);

			h_eff_Injet_matched[cent_bin]->Fill(jet_pt, pt, fabs(rapidity), event_weight);
			h_eff_Injet_matched[nCentBins-1]->Fill(jet_pt, pt, fabs(rapidity), event_weight);

			float eff_weight = 1;
			if (isFirstPass)
			{
//                eff_weight = trkcorr->get_effcorr(pt, eta, cent_bin, 0, _dataset);
                h_trk_foreff_matched[cent_bin]->Fill(phi, pt, eta, event_weight*eff_weight);
				h_trk_foreff_matched[nCentBins-1]->Fill(phi, pt, eta, event_weight*eff_weight);
			}

		}

		//get the MC tracks....
		truth_itr = particles->begin();
		for( ; truth_itr!=truth_end; ++truth_itr)
		{
			if( fabs((*truth_itr)->charge())<0.5 ) continue;
			if( ((*truth_itr)->status())!=1) continue;

			int ty=getTypeTruth((*truth_itr)->barcode(),(*truth_itr)->pdgId(),(*truth_itr)->status(),(*truth_itr)->charge());
			if(ty!=1 && ty!=5) continue;

            
			//get the tracks....
			float eta = (*truth_itr)->eta();
			float phi = (*truth_itr)->phi();
			float pt = (*truth_itr)->pt()/ 1000.0;

			float R = DeltaR(phi,eta,jet_phi,jet_eta);
			if (fabs(eta)>=2.5) continue;
			if (R > _dR_max) continue;

			int deta_bin = trkcorr->GetdRBin(fabs(DeltaEta(eta,jet_eta)));
			h_eff_deta.at(deta_bin).at(cent_bin)->Fill(jet_pt,pt,fabs(rapidity), event_weight);
            h_eff_deta.at(deta_bin).at(nCentBins-1)->Fill(jet_pt,pt,fabs(rapidity), event_weight);

			int dphi_bin = trkcorr->GetdRBin(fabs(DeltaPhi(phi,jet_phi)));
			h_eff_dphi.at(dphi_bin).at(cent_bin)->Fill(jet_pt,pt,fabs(rapidity), event_weight);
            h_eff_dphi.at(dphi_bin).at(nCentBins-1)->Fill(jet_pt,pt,fabs(rapidity), event_weight);

			int dR_bin = trkcorr->GetdRBin(R);
			h_eff_dR.at(dR_bin).at(cent_bin)->Fill(jet_pt,pt,fabs(rapidity), event_weight);
            h_eff_dR.at(dR_bin).at(nCentBins-1)->Fill(jet_pt,pt,fabs(rapidity), event_weight);

			h_eff_total.at(cent_bin)->Fill(jet_pt, pt, eta, event_weight);
            h_eff_total.at(nCentBins-1)->Fill(jet_pt, pt, eta, event_weight);
            
            h_trk_eff_map.at(cent_bin)->Fill(pt, eta, phi, event_weight);
            h_trk_eff_map.at(nCentBins-1)->Fill(pt, eta, phi, event_weight);

            h_eff_entries.at(cent_bin)->Fill(jet_pt, pt, eta);
            h_eff_entries.at(nCentBins-1)->Fill(jet_pt, pt, eta);

			h_eff_Injet[cent_bin]->Fill(jet_pt,pt,fabs(rapidity), event_weight);
			h_eff_Injet[nCentBins-1]->Fill(jet_pt,pt, fabs(rapidity), event_weight);
			
            h_eff_Injet_entries[cent_bin]->Fill(jet_pt,pt,fabs(rapidity));
			h_eff_Injet_entries[nCentBins-1]->Fill(jet_pt,pt,fabs(rapidity));

			if (isFirstPass)
			{
				h_trk_foreff_full[cent_bin]->Fill(phi, pt, eta, event_weight);
				h_trk_foreff_full[nCentBins-1]->Fill(phi, pt, eta, event_weight);
                
                h_trk_foreff_entries[cent_bin]->Fill(phi, pt, eta);
                h_trk_foreff_entries[nCentBins-1]->Fill(phi, pt, eta);
                h_truth_trk_map->Fill(pt,eta,phi, event_weight);
			}
		}

		isFirstPass = false;
	}

	// Clear vectors

	trk_good_matched_eta.clear();
	trk_good_matched_phi.clear();
	trk_good_matched_pt.clear();
	trk_good_matched_pt_reco.clear();
	trk_good_eta.clear();
	trk_good_phi.clear();
	trk_good_pt.clear();
	trk_good_isMatched.clear();
	trk_good_vertextype.clear();

	jet_pt_xcalib_vector.clear();
	jet_phi_vector.clear();
	jet_eta_vector.clear();
	truth_jet_pt_vector.clear();
	truth_jet_eta_vector.clear();
	truth_jet_phi_vector.clear();
	TruthJetIndex.clear();


	trk_good_matched_eta.shrink_to_fit();
	trk_good_matched_phi.shrink_to_fit();
	trk_good_matched_pt.shrink_to_fit();
	trk_good_matched_pt_reco.shrink_to_fit();
	trk_good_eta.shrink_to_fit();
	trk_good_phi.shrink_to_fit();
	trk_good_pt.shrink_to_fit();
	trk_good_isMatched.shrink_to_fit();
	trk_good_vertextype.shrink_to_fit();

	jet_pt_xcalib_vector.shrink_to_fit();
	jet_phi_vector.shrink_to_fit();
	jet_eta_vector.shrink_to_fit();
	truth_jet_pt_vector.shrink_to_fit();
	truth_jet_eta_vector.shrink_to_fit();
	truth_jet_phi_vector.shrink_to_fit();
	TruthJetIndex.shrink_to_fit();

	return EL::StatusCode::SUCCESS;
}

EL::StatusCode TrackingPerformance :: postExecute (){
	return EL::StatusCode::SUCCESS;
}

EL::StatusCode TrackingPerformance :: finalize (){
	xAOD::TEvent* event = wk()->xaodEvent();
		EL_RETURN_CHECK( "Finalize", m_trackSelectorTool->finalize() );

	//cleaning cleaning :)
	if( m_jetCleaning ) delete m_jetCleaning; m_jetCleaning = 0;
	if( m_jetCalibration ) delete m_jetCalibration; m_jetCalibration = 0;
	if( m_trackSelectorTool ) delete m_trackSelectorTool; m_trackSelectorTool=0;
	
	return EL::StatusCode::SUCCESS;
}

EL::StatusCode TrackingPerformance :: histFinalize (){
	cout<<"Events = "<< m_eventCounter<<endl;
	return EL::StatusCode::SUCCESS;
}
