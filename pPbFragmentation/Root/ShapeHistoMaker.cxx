#define ShapeToolInit_cxx
#include "pPbFragmentation/PbPbFFShape.h"
#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>

using namespace std;
using namespace MTCorrector;

EL::StatusCode PbPbFFShape :: histInitialize ()
{
	cout << " Setting  histograms" << endl;
	TH3D* temphist_3D = nullptr;
	TH2D* temphist_2D = nullptr;
	TH1D* temphist_1D = nullptr;


	//Track corrector
	int _nCentbins = GetCentralityNBins(_centrality_scheme);
	trkcorr = new TrackCorrector(_cut_level.c_str(), GetCentralityNBins(31)-1,_eff_jety);
	trkcorr->drmax = _dR_max;
	trkcorr->InitdRBinRange();

	//Jet corrector
	jetcorr = new JetCorrector();
	int _nJetYBins = jetcorr->nJetYBins;
	int _ndRBins = trkcorr->ndRBins + 1;

	jetcorr->JERcut = _JERBalancecut;
	jetcorr->min_jet_pt = _pTjetCut+1;
	jetcorr->max_jet_pt = 500.; //Maximum jet pt for reweighting
	jetcorr->m_isMB = _isMB;
    jetcorr->is_pp = (_dataset == 3);


	int ptJetBinsN, etaJetBinsN, phiJetBinsN, ptTrkBinsN, etaTrkBinsN, phiTrkBinsN, zBinsN, zBinsFineN, d0z0BinsN, respBinsN, finehitsBinsN, dR_resBinsN;
	double ptJetBins[1000], etaJetBins[1000], phiJetBins[1000], ptTrkBins[1000], etaTrkBins[1000], phiTrkBins[1000], zBins[1000], zBinsFine[1000],d0z0Bins[1000], respBins[1000], finehitsBins[1000], dR_resBins[1000];

	SetupBinning(0, "pt-jet-PbPb", ptJetBins, ptJetBinsN);
	SetupBinning(0, "eta-jet", etaJetBins, etaJetBinsN);
	SetupBinning(0, "phi-trk", phiJetBins, phiJetBinsN);
	SetupBinning(0, "pt-trk-shape", ptTrkBins, ptTrkBinsN);
	SetupBinning(0, "eta-trk", etaTrkBins, etaTrkBinsN);
	SetupBinning(0, "phi-trk", phiTrkBins, phiTrkBinsN);
	SetupBinning(0, "z-rebin", zBins, zBinsN);
	SetupBinning(0, "z", zBinsFine, zBinsFineN);
	SetupBinning(0, "d0z0", d0z0Bins, d0z0BinsN);
	SetupBinning(0, "resp", respBins, respBinsN);
	SetupBinning(0, "hits_fine", finehitsBins, finehitsBinsN);
	SetupBinning(0, "dr_fine", dR_resBins, dR_resBinsN);

	Double_t PVBins[3]={0,1,2};
	int PVBinsN=2;


	//Axis histograms
	h_jet_pt_eta_phi = new TH3D("h_jet_pt_eta_phi","h_jet_pt_eta_phi",ptJetBinsN, ptJetBins, etaJetBinsN, etaJetBins, phiJetBinsN, phiJetBins);
	h_trk_pt_eta_phi = new TH3D("h_trk_pt_eta_phi","h_trk_pt_eta_phi",ptTrkBinsN, ptTrkBins, etaTrkBinsN, etaTrkBins, phiTrkBinsN, phiTrkBins);
	h_z_zfine_d0z0 = new TH3D("h_z_zfine_d0z0","h_z_zfine_d0z0",zBinsN, zBins, zBinsFineN, zBinsFine, d0z0BinsN, d0z0Bins);
	h_resp_hitsfine_drfine = new TH3D("h_resp_hitsfine_drfine","h_resp_hitsfine_drfine",respBinsN, respBins, finehitsBinsN, finehitsBins, dR_resBinsN, dR_resBins);
	wk()->addOutput (h_jet_pt_eta_phi);
	wk()->addOutput (h_trk_pt_eta_phi);
	wk()->addOutput (h_z_zfine_d0z0);
	wk()->addOutput (h_resp_hitsfine_drfine);

	//Debugging histograms
	for (int i=0;i<GetCentralityNBins(_centrality_scheme);i++)
	{
		temphist_3D = new TH3D(Form("h_reco_pre_truth_match_cent%i",i),Form("h_reco_pre_truth_match_cent%i",i),ptJetBinsN, ptJetBins, etaJetBinsN, etaJetBins, etaJetBinsN, etaJetBins);
		temphist_3D->Sumw2();
		h_reco_pre_truth_match.push_back(temphist_3D);

		temphist_3D = new TH3D(Form("h_reco_post_truth_match_cent%i",i),Form("h_reco_post_truth_match_cent%i",i),ptJetBinsN, ptJetBins, etaJetBinsN, etaJetBins, etaJetBinsN, etaJetBins);
		temphist_3D->Sumw2();
		h_reco_post_truth_match.push_back(temphist_3D);

		temphist_3D = new TH3D(Form("h_reco_truth_matched_cent%i",i),Form("h_reco_truth_matched_cent%i",i),ptJetBinsN, ptJetBins, ptJetBinsN, ptJetBins, dR_resBinsN, dR_resBins);
		temphist_3D->Sumw2();
		h_reco_truth_matched.push_back(temphist_3D);

		temphist_3D = new TH3D(Form("h_reco_truth_comparison_cent%i",i),Form("h_reco_truth_comparison_cent%i",i), ptJetBinsN, ptJetBins, respBinsN, respBins, etaJetBinsN, etaJetBins);
		temphist_3D->Sumw2();
		h_reco_truth_comparison.push_back(temphist_3D);

		temphist_3D = new TH3D(Form("h_jet_for_eff_cent%i",i),Form("h_jet_for_eff_cent%i",i), ptJetBinsN, ptJetBins, etaJetBinsN, etaJetBins, phiJetBinsN, phiJetBins);
		temphist_3D->Sumw2();
		h_jet_for_eff.push_back(temphist_3D);

		temphist_3D = new TH3D(Form("h_jet_for_eff_full_cent%i",i),Form("h_jet_for_eff_full_cent%i",i), ptJetBinsN, ptJetBins, etaJetBinsN, etaJetBins, phiJetBinsN, phiJetBins);
		temphist_3D->Sumw2();
		h_jet_for_eff_full.push_back(temphist_3D);

		wk()->addOutput (h_reco_pre_truth_match.at(i));
		wk()->addOutput (h_reco_post_truth_match.at(i));
		wk()->addOutput (h_reco_truth_matched.at(i));
		wk()->addOutput (h_reco_truth_comparison.at(i));
		wk()->addOutput (h_jet_for_eff.at(i));
		wk()->addOutput (h_jet_for_eff_full.at(i));

	}

	//Basic histograms
	h_FCal_Et = new TH1D("h_FCal_Et",";FCal E_{T};N",140,0,7);
	h_FCal_Et->Sumw2();

	h_RejectionHisto = new TH1D("RejectionHisto","RejectionHisto",9,0,9);
	SetRejectionHistogram(h_RejectionHisto);

	h_centrality = new TH1D("Centrality","Centrality",10,0,10);
	h_centrality->Sumw2();

	hET_ETsub = new TH3D("hET_ETsub","hET_ETsub",200,0,200,100,-5,5,200,-50,50);
	hET_ETsub->Sumw2();

	h_dR_binning = new TH1D("h_dR_binning","h_dR_binning",_ndRBins-1,trkcorr->dRrange);
	h_dR_binning->Sumw2();

	h_triggercounter = new TH2D("h_triggercounter","h_triggercounter",_nTriggers,0,_nTriggers,2,-0.5,1.5);
	h_triggercounter->Sumw2();
	SetTrigger_hist(h_triggercounter);

	h_reco_trk_map = new TH3D("h_reco_trk_map","h_reco_trk_map;p_{T};#eta;#phi",ptTrkBinsN, ptTrkBins,etaTrkBinsN,etaTrkBins,phiTrkBinsN,phiTrkBins);
	h_reco_trk_map->Sumw2();

	h_reco_trk_map_nocuts = new TH3D("h_reco_trk_map_nocuts","h_reco_trk_map_nocuts;p_{T};#eta;#phi",ptTrkBinsN, ptTrkBins,etaTrkBinsN,etaTrkBins,phiTrkBinsN,phiTrkBins);
	h_reco_trk_map_nocuts->Sumw2();

	wk()->addOutput (h_FCal_Et);
	wk()->addOutput (h_RejectionHisto);
	wk()->addOutput (hET_ETsub);
	wk()->addOutput (h_triggercounter);
	wk()->addOutput (h_reco_trk_map);
	wk()->addOutput (h_reco_trk_map_nocuts);
	wk()->addOutput (h_centrality);
	wk()->addOutput (h_dR_binning);
	
	//wk()->addOutput (deriv_val);


	const int NCuts=8;
	string CutsOn[NCuts]={"Reco","d0","z0sintheta","SI_hits","Pix_holes","SI_holes","Shared_SI_hits","Gen"};

	/*
	 for (int i=0;i<NCuts;i++)
	 {
		temphist_3D = new TH3D(Form("h_cut_flow_%s",CutsOn[i].c_str()),Form("h_cut_flow_%s",CutsOn[i].c_str()),ptTrkBinsN, ptTrkBins,etaTrkBinsN,etaTrkBins,PVBinsN,PVBins);
		temphist_3D->Sumw2();
		h_cut_flow.push_back(temphist_3D);
		wk()->addOutput (h_cut_flow.at(i));
	 }
	 */
	for (int i=0;i<GetCentralityNBins(_centrality_scheme);i++)
	{
		temphist_3D = new TH3D(Form("h_PixHits_cent%i",i),Form("h_PixHits_cent%i",i),ptTrkBinsN, ptTrkBins,etaTrkBinsN, etaTrkBins,finehitsBinsN,finehitsBins);
		temphist_3D->Sumw2();
		h_PixHits.push_back(temphist_3D);

		temphist_3D = new TH3D(Form("h_SCTHits_cent%i",i),Form("h_SCTHits_cent%i",i),ptTrkBinsN, ptTrkBins,etaTrkBinsN, etaTrkBins,finehitsBinsN,finehitsBins);
		temphist_3D->Sumw2();
		h_SCTHits.push_back(temphist_3D);

		temphist_2D = new TH2D(Form("h_d0_cent%i",i),Form("h_d0_cent%i",i),ptTrkBinsN,ptTrkBins,d0z0BinsN,d0z0Bins);
		temphist_2D->Sumw2();
		h_d0.push_back(temphist_2D);

		temphist_2D = new TH2D(Form("h_z0sintheta_cent%i",i),Form("h_z0sintheta_cent%i",i),ptTrkBinsN, ptTrkBins,d0z0BinsN,d0z0Bins);
		temphist_2D->Sumw2();
		h_z0sintheta.push_back(temphist_2D);

		temphist_3D = new TH3D(Form("h_vx_cent%i",i),"h_vx;v_{x};v_{y};v_{z}",200,-2,2,200,-2,2,200,-200,200);
		temphist_3D->Sumw2();
		h_vx.push_back(temphist_3D);

		temphist_3D = new TH3D(Form("h_jetpT_v_multiplicity_cent%i",i),"h_jetpT_v_multiplicity;jet p_{T} GeV;p_{T}^{trk,min} [GeV];multiplicity",20,0,1000,6,0,6,50,0,50);
		h_jetpT_v_multiplicity.push_back(temphist_3D);
		h_jetpT_v_multiplicity.at(i)->Sumw2();

		//temphist_3D = new TH3D(Form("h_BL_cent%i",i),"h_BL;BL_{x};BL_{y};BL_{z}",200,-1,1,200,-1,1,200,-10,10);
		//h_BL.push_back(temphist_3D);

		temphist_2D = new TH2D(Form("h_R2vR4_cent%i",i),Form("h_R2vR4_cent%i",i),40,100,500,80,-0.4,0.4);
		temphist_2D->Sumw2();
		h_R2vR4.push_back(temphist_2D);

		if (_data_switch==1){
			temphist_3D = new TH3D(Form("JES_v_max_z_cent%i",i),Form("JES_v_max_z_cent%i",i),zBinsN, zBins,respBinsN,respBins, ptJetBinsN, ptJetBins);
			temphist_3D->Sumw2();
			JES_v_max_z.push_back(temphist_3D);
			temphist_3D = new TH3D(Form("JES_v_max_pT_cent%i",i),Form("JES_v_max_pT_cent%i",i),ptTrkBinsN, ptTrkBins, respBinsN,respBins, ptJetBinsN, ptJetBins);
			temphist_3D->Sumw2();
			JES_v_max_pT.push_back(temphist_3D);
			temphist_3D = new TH3D(Form("h_jet_pos_v_zmax_cent%i",i),Form("h_jet_pos_v_zmax_cent%i",i),zBinsN, zBins, dR_resBinsN,dR_resBins, ptJetBinsN, ptJetBins);
			temphist_3D->Sumw2();
			h_jet_pos_v_zmax.push_back(temphist_3D);
			temphist_3D = new TH3D(Form("h_jet_pos_v_truth_zmax_cent%i",i),Form("h_jet_pos_v_truth_zmax_cent%i",i),zBinsN, zBins, dR_resBinsN,dR_resBins, ptJetBinsN, ptJetBins);
			temphist_3D->Sumw2();
			h_jet_pos_v_truth_zmax.push_back(temphist_3D);
			temphist_3D = new TH3D(Form("h_jet_pos_v_ptmax_cent%i",i),Form("h_jet_pos_v_ptmax_cent%i",i),ptTrkBinsN, ptTrkBins, dR_resBinsN,dR_resBins, ptJetBinsN, ptJetBins);
			temphist_3D->Sumw2();
			h_jet_pos_v_ptmax.push_back(temphist_3D);
			temphist_3D = new TH3D(Form("h_jet_pos_v_truth_ptmax_cent%i",i),Form("h_jet_pos_v_truth_ptmax_cent%i",i),ptTrkBinsN, ptTrkBins, dR_resBinsN,dR_resBins, ptJetBinsN, ptJetBins);
			temphist_3D->Sumw2();
			h_jet_pos_v_truth_ptmax.push_back(temphist_3D);

			temphist_3D = new TH3D(Form("h_jet_pos_v_pt_cent%i",i),Form("h_jet_pos_v_pt_cent%i",i),ptTrkBinsN, ptTrkBins, dR_resBinsN,dR_resBins, ptJetBinsN, ptJetBins);
			temphist_3D->Sumw2();
			h_jet_pos_v_pt.push_back(temphist_3D);
			temphist_3D = new TH3D(Form("h_jet_pos_v_truth_pt_cent%i",i),Form("h_jet_pos_v_truth_pt_cent%i",i),ptTrkBinsN, ptTrkBins, dR_resBinsN,dR_resBins, ptJetBinsN, ptJetBins);
			temphist_3D->Sumw2();
			h_jet_pos_v_truth_pt.push_back(temphist_3D);

		}

		wk()->addOutput (h_PixHits.at(i));
		wk()->addOutput (h_SCTHits.at(i));
		wk()->addOutput (h_d0.at(i));
		wk()->addOutput (h_z0sintheta.at(i));
		wk()->addOutput (h_vx.at(i));

		wk()->addOutput (h_R2vR4.at(i));

		wk()->addOutput (h_jetpT_v_multiplicity.at(i));

		if (_data_switch==1){
			wk()->addOutput (JES_v_max_pT.at(i));
			wk()->addOutput (JES_v_max_z.at(i));
			wk()->addOutput (h_jet_pos_v_zmax.at(i));
			wk()->addOutput (h_jet_pos_v_truth_zmax.at(i));
			wk()->addOutput (h_jet_pos_v_ptmax.at(i));
			wk()->addOutput (h_jet_pos_v_truth_ptmax.at(i));

		}
		//wk()->addOutput (h_BL.at(i));

	}



	ff_raw =  vector<vector<TH2D*> > (_ndRBins, vector<TH2D*>(_nCentbins));
	ff_raw_UE =  vector<vector<TH2D*> > (_ndRBins, vector<TH2D*>(_nCentbins));
	ChPS_raw =  vector<vector<TH2D*> > (_ndRBins, vector<TH2D*>(_nCentbins));
	ChPS_raw_rr =  vector<vector<TH2D*> > (_ndRBins, vector<TH2D*>(_nCentbins));
	ChPS_raw_rt =  vector<vector<TH2D*> > (_ndRBins, vector<TH2D*>(_nCentbins));
	ChPS_raw_tr =  vector<vector<TH2D*> > (_ndRBins, vector<TH2D*>(_nCentbins));
	ChPS_raw_tt =  vector<vector<TH2D*> > (_ndRBins, vector<TH2D*>(_nCentbins));
	ChPS_raw_tt_deta =  vector<vector<TH2D*> > (_ndRBins, vector<TH2D*>(_nCentbins));
	ChPS_raw_tt_dphi =  vector<vector<TH2D*> > (_ndRBins, vector<TH2D*>(_nCentbins));
	ChPS_raw_tt_mod =  vector<vector<TH2D*> > (_ndRBins, vector<TH2D*>(_nCentbins));
	ChPS_raw_UE =  vector<vector<TH2D*> > (_ndRBins, vector<TH2D*>(_nCentbins));
	h_dR_change =  vector<vector<TH3D*> > (ptJetBinsN, vector<TH3D*>(_nCentbins));
	UE_distr =  vector<vector<TH3D*> > (ptJetBinsN, vector<TH3D*>(_nCentbins));


	h_reco_jet_spectrum =  vector<vector<TH1D*> > (_nJetYBins, vector<TH1D*>(_nCentbins));
	//in MC only
	if (_data_switch==1)
	{
		ff_UE_z =  vector<vector<TH2D*> > (_ndRBins, vector<TH2D*>(_nCentbins));
		ff_UE_pT =  vector<vector<TH2D*> > (_ndRBins, vector<TH2D*>(_nCentbins));
		ff_truth =  vector<vector<TH2D*> > (_ndRBins, vector<TH2D*>(_nCentbins));
		ChPS_truth =  vector<vector<TH2D*> > (_ndRBins, vector<TH2D*>(_nCentbins));
		ChPS_truth_deta =  vector<vector<TH2D*> > (_ndRBins, vector<TH2D*>(_nCentbins));
		ChPS_truth_dphi =  vector<vector<TH2D*> > (_ndRBins, vector<TH2D*>(_nCentbins));

		reco_posRes_ChPS =  vector<vector<TH2D*> > (_ndRBins, vector<TH2D*>(_nCentbins));
		rt_posRes_ChPS =  vector<vector<TH2D*> > (_ndRBins, vector<TH2D*>(_nCentbins));
		tr_posRes_ChPS =  vector<vector<TH2D*> > (_ndRBins, vector<TH2D*>(_nCentbins));
		truth_posRes_ChPS =  vector<vector<TH2D*> > (_ndRBins, vector<TH2D*>(_nCentbins));

		ff_trackpTResponse =  vector<vector<TH3D*> > (_ndRBins, vector<TH3D*>(_nCentbins));
		ff_trackzResponse =  vector<vector<TH3D*> > (_ndRBins, vector<TH3D*>(_nCentbins));
		response_FF =  vector<vector<RooUnfoldResponse*> > (_ndRBins, vector<RooUnfoldResponse*>(_nCentbins));
		response_ChPS =  vector<vector<RooUnfoldResponse*> > (_ndRBins, vector<RooUnfoldResponse*>(_nCentbins));

		h_true_jet_spectrum =  vector<vector<TH1D*> > (_nJetYBins, vector<TH1D*>(_nCentbins));
		h_true_jet_spectrum_matched =  vector<vector<TH1D*> > (_nJetYBins, vector<TH1D*>(_nCentbins));
		h_reco_jet_spectrum_matched =  vector<vector<TH1D*> > (_nJetYBins, vector<TH1D*>(_nCentbins));
		ff_jetResponse =  vector<vector<TH2D*> > (_nJetYBins, vector<TH2D*>(_nCentbins));
		response_jet =  vector<vector<RooUnfoldResponse*> > (_nJetYBins, vector<RooUnfoldResponse*>(_nCentbins));
	}


	for (int j=0;j<_nCentbins;j++)
	{

		for (int i=0;i<ptJetBinsN;i++)
		{
			temphist_3D = new TH3D(Form("h_dR_change_jetpt%i_cent%i",i,j),Form("h_dR_change_jetpt%i_cent%i",i,j),_ndRBins-1,trkcorr->dRrange,_ndRBins-1,trkcorr->dRrange,ptTrkBinsN, ptTrkBins);
			h_dR_change.at(i).at(j) = temphist_3D;
			h_dR_change.at(i).at(j)->Sumw2();
			wk()->addOutput (h_dR_change.at(i).at(j));

			temphist_3D = new TH3D(Form("h_UE_distr_jetpt%i_cent%i",i,j),Form("h_UE_distr_jetpt%i_cent%i",i,j),_ndRBins-1,trkcorr->dRrange, ptTrkBinsN, ptTrkBins, etaTrkBinsN, etaTrkBins);
			UE_distr.at(i).at(j) = temphist_3D;
			UE_distr.at(i).at(j)->Sumw2();
			wk()->addOutput (UE_distr.at(i).at(j));


		}

		for (int i=0;i<_ndRBins;i++)
		{
			//D(z)
			temphist_2D = new TH2D(Form("ff_raw_0_dR%i_cent%i",i,j),Form("ff_raw_0_dR%i_cent%i",i,j),zBinsN, zBins, ptJetBinsN, ptJetBins);
			ff_raw.at(i).at(j) = temphist_2D;

			temphist_2D = new TH2D(Form("ff_raw_1_dR%i_cent%i",i,j),Form("ff_raw_1_dR%i_cent%i",i,j),zBinsN, zBins, ptJetBinsN, ptJetBins);
			ff_raw_UE.at(i).at(j) = temphist_2D;

			//D(Pt)
			temphist_2D = new TH2D(Form("ChPS_raw_0_dR%i_cent%i",i,j),Form("ChPS_raw_0_dR%i_cent%i",i,j),ptTrkBinsN, ptTrkBins, ptJetBinsN, ptJetBins);
			ChPS_raw.at(i).at(j) = temphist_2D;

			temphist_2D = new TH2D(Form("ChPS_raw_rr_dR%i_cent%i",i,j),Form("ChPS_raw_rr_dR%i_cent%i",i,j),ptTrkBinsN, ptTrkBins, ptJetBinsN, ptJetBins);
			ChPS_raw_rr.at(i).at(j) = temphist_2D;

			temphist_2D = new TH2D(Form("ChPS_raw_rt_dR%i_cent%i",i,j),Form("ChPS_raw_rt_dR%i_cent%i",i,j),ptTrkBinsN, ptTrkBins, ptJetBinsN, ptJetBins);
			ChPS_raw_rt.at(i).at(j) = temphist_2D;

			temphist_2D = new TH2D(Form("ChPS_raw_tr_dR%i_cent%i",i,j),Form("ChPS_raw_tr_dR%i_cent%i",i,j),ptTrkBinsN, ptTrkBins, ptJetBinsN, ptJetBins);
			ChPS_raw_tr.at(i).at(j) = temphist_2D;

			temphist_2D = new TH2D(Form("ChPS_raw_tt_dR%i_cent%i",i,j),Form("ChPS_raw_tt_dR%i_cent%i",i,j),ptTrkBinsN, ptTrkBins, ptJetBinsN, ptJetBins);
			ChPS_raw_tt.at(i).at(j) = temphist_2D;

			temphist_2D = new TH2D(Form("ChPS_raw_tt_deta%i_cent%i",i,j),Form("ChPS_raw_tt_deta%i_cent%i",i,j),ptTrkBinsN, ptTrkBins, ptJetBinsN, ptJetBins);
			ChPS_raw_tt_deta.at(i).at(j) = temphist_2D;

			temphist_2D = new TH2D(Form("ChPS_raw_tt_dphi%i_cent%i",i,j),Form("ChPS_raw_tt_dphi%i_cent%i",i,j),ptTrkBinsN, ptTrkBins, ptJetBinsN, ptJetBins);
			ChPS_raw_tt_dphi.at(i).at(j) = temphist_2D;

			temphist_2D = new TH2D(Form("ChPS_raw_1_dR%i_cent%i",i,j),Form("ChPS_raw_1_dR%i_cent%i",i,j),ptTrkBinsN, ptTrkBins, ptJetBinsN, ptJetBins);
			ChPS_raw_UE.at(i).at(j) = temphist_2D;

			temphist_2D = new TH2D(Form("ChPS_raw_tt_mod_dR%i_cent%i",i,j),Form("ChPS_raw_tt_mod_dR%i_cent%i",i,j),ptTrkBinsN, ptTrkBins, ptJetBinsN, ptJetBins);
			ChPS_raw_tt_mod.at(i).at(j) = temphist_2D;

			ff_raw.at(i).at(j)->Sumw2();
			ff_raw_UE.at(i).at(j)->Sumw2();

			ChPS_raw.at(i).at(j)->Sumw2();
			ChPS_raw_rr.at(i).at(j)->Sumw2();
			ChPS_raw_rt.at(i).at(j)->Sumw2();
			ChPS_raw_tr.at(i).at(j)->Sumw2();
			ChPS_raw_tt.at(i).at(j)->Sumw2();
			ChPS_raw_tt_deta.at(i).at(j)->Sumw2();
			ChPS_raw_tt_dphi.at(i).at(j)->Sumw2();
			ChPS_raw_tt_mod.at(i).at(j)->Sumw2();
			ChPS_raw_UE.at(i).at(j)->Sumw2();

			wk()->addOutput (ff_raw.at(i).at(j));
			wk()->addOutput (ff_raw_UE.at(i).at(j));
			wk()->addOutput (ChPS_raw.at(i).at(j));
			wk()->addOutput (ChPS_raw_rr.at(i).at(j));
			wk()->addOutput (ChPS_raw_tr.at(i).at(j));
			wk()->addOutput (ChPS_raw_rt.at(i).at(j));
			wk()->addOutput (ChPS_raw_tt.at(i).at(j));
			wk()->addOutput (ChPS_raw_tt_deta.at(i).at(j));
			wk()->addOutput (ChPS_raw_tt_dphi.at(i).at(j));
			wk()->addOutput (ChPS_raw_tt_mod.at(i).at(j));
			wk()->addOutput (ChPS_raw_UE.at(i).at(j));

			if (_data_switch==1)
			{
				temphist_2D = new TH2D(Form("ff_UE_z_dR%i_cent%i",i,j),Form("ff_UE_z_dR%i_cent%i",i,j),zBinsN, zBins, ptJetBinsN, ptJetBins);
				ff_UE_z.at(i).at(j) = temphist_2D;
				temphist_2D = new TH2D(Form("ff_UE_pT_dR%i_cent%i",i,j),Form("ff_UE_pT_dR%i_cent%i",i,j),ptTrkBinsN, ptTrkBins, ptJetBinsN, ptJetBins);
				ff_UE_pT.at(i).at(j) = temphist_2D;

				//Truth
				temphist_2D = new TH2D(Form("ff_truth_dR%i_cent%i",i,j),Form("ff_truth_dR%i_cent%i",i,j),zBinsN, zBins, ptJetBinsN, ptJetBins);
				ff_truth.at(i).at(j) = temphist_2D;
				temphist_2D = new TH2D(Form("ChPS_truth_dR%i_cent%i",i,j),Form("ChPS_truth_dR%i_cent%i",i,j),ptTrkBinsN, ptTrkBins, ptJetBinsN, ptJetBins);
				ChPS_truth.at(i).at(j) = temphist_2D;
				temphist_2D = new TH2D(Form("ChPS_truth_deta%i_cent%i",i,j),Form("ChPS_truth_deta%i_cent%i",i,j),ptTrkBinsN, ptTrkBins, ptJetBinsN, ptJetBins);
				ChPS_truth_deta.at(i).at(j) = temphist_2D;
				temphist_2D = new TH2D(Form("ChPS_truth_dphi%i_cent%i",i,j),Form("ChPS_truth_dphi%i_cent%i",i,j),ptTrkBinsN, ptTrkBins, ptJetBinsN, ptJetBins);
				ChPS_truth_dphi.at(i).at(j) = temphist_2D;


				temphist_3D = new TH3D(Form("ff_trackpTResponse_dR%i_cent%i",i,j),Form("ff_trackpTResponse_dR%i_cent%i",i,j),ptTrkBinsN, ptTrkBins,ptTrkBinsN, ptTrkBins, ptJetBinsN, ptJetBins);
				ff_trackpTResponse.at(i).at(j) = temphist_3D;
				temphist_3D = new TH3D(Form("ff_trackzResponse_dR%i_cent%i",i,j),Form("ff_trackzResponse_dR%i_cent%i",i,j),zBinsN, zBins, zBinsN, zBins, ptJetBinsN, ptJetBins);
				ff_trackzResponse.at(i).at(j) = temphist_3D;

				//Responses
				response_FF.at(i).at(j) = new RooUnfoldResponse(ff_raw.at(i).at(j),ff_truth.at(i).at(j));
				response_ChPS.at(i).at(j) = new RooUnfoldResponse(ChPS_raw.at(i).at(j),ChPS_truth.at(i).at(j));


				//Truth
				temphist_2D = new TH2D(Form("reco_posRes_ChPS_dR%i_cent%i",i,j),Form("reco_posRes_ChPS_dR%i_cent%i",i,j),ptTrkBinsN, ptTrkBins, ptJetBinsN, ptJetBins);
				reco_posRes_ChPS.at(i).at(j) = temphist_2D;

				temphist_2D = new TH2D(Form("rt_posRes_ChPS_dR%i_cent%i",i,j),Form("rt_posRes_ChPS_dR%i_cent%i",i,j),ptTrkBinsN, ptTrkBins, ptJetBinsN, ptJetBins);
				rt_posRes_ChPS.at(i).at(j) = temphist_2D;

				temphist_2D = new TH2D(Form("tr_posRes_ChPS_dR%i_cent%i",i,j),Form("tr_posRes_ChPS_dR%i_cent%i",i,j),ptTrkBinsN, ptTrkBins, ptJetBinsN, ptJetBins);
				tr_posRes_ChPS.at(i).at(j) = temphist_2D;

				temphist_2D = new TH2D(Form("truth_posRes_ChPS_dR%i_cent%i",i,j),Form("truth_posRes_ChPS_dR%i_cent%i",i,j),ptTrkBinsN, ptTrkBins, ptJetBinsN, ptJetBins);
				truth_posRes_ChPS.at(i).at(j) = temphist_2D;

				ff_UE_z.at(i).at(j)->Sumw2();
				ff_UE_pT.at(i).at(j)->Sumw2();
				ff_truth.at(i).at(j)->Sumw2();
				ChPS_truth.at(i).at(j)->Sumw2();
				ChPS_truth_deta.at(i).at(j)->Sumw2();
				ChPS_truth_dphi.at(i).at(j)->Sumw2();

				ff_trackpTResponse.at(i).at(j)->Sumw2();
				ff_trackzResponse.at(i).at(j)->Sumw2();

				reco_posRes_ChPS.at(i).at(j)->Sumw2();
				rt_posRes_ChPS.at(i).at(j)->Sumw2();
				tr_posRes_ChPS.at(i).at(j)->Sumw2();
				truth_posRes_ChPS.at(i).at(j)->Sumw2();

				wk()->addOutput (ff_UE_z.at(i).at(j));
				wk()->addOutput (ff_UE_pT.at(i).at(j));
				wk()->addOutput (ff_truth.at(i).at(j));
				wk()->addOutput (ChPS_truth.at(i).at(j));
				wk()->addOutput (ChPS_truth_deta.at(i).at(j));
				wk()->addOutput (ChPS_truth_dphi.at(i).at(j));

				wk()->addOutput (ff_trackpTResponse.at(i).at(j));
				wk()->addOutput (ff_trackzResponse.at(i).at(j));

				wk()->addOutput (response_FF.at(i).at(j));
				wk()->addOutput (response_ChPS.at(i).at(j));

				wk()->addOutput (reco_posRes_ChPS.at(i).at(j));
				wk()->addOutput (rt_posRes_ChPS.at(i).at(j));
				wk()->addOutput (tr_posRes_ChPS.at(i).at(j));
				wk()->addOutput (truth_posRes_ChPS.at(i).at(j));



			}

		}

		for (int i=0;i<_nJetYBins;i++)
		{

			//temphist_2D = new TH2D(Form("ff_truth_matched_y%i_cent%i",i,j),Form("ff_truth_matched_y%i_cent%i",i,j),zBinsN, zBins, ptJetBinsN, ptJetBins);
			//ff_truth_matched.at(i).at(j) = temphist_2D;
			//temphist_2D = new TH2D(Form("ChPS_matched_truth_y%i_cent%i",i,j),Form("ChPS_matched_truth_y%i_cent%i",i,j),ptTrkBinsN, ptTrkBins, ptJetBinsN, ptJetBins);
			//ChPS_truth_matched.at(i).at(j) = temphist_2D;

			//Jet spectra
			temphist_1D = new TH1D(Form("h_reco_jet_spectrum_y%i_cent%i",i,j),Form("h_reco_jet_spectrum_y%i_cent%i",i,j),ptJetBinsN, ptJetBins);
			h_reco_jet_spectrum.at(i).at(j) = temphist_1D;

			//TODO to be enabled when needed
			//temphist_1D = new TH1D(Form("h_reco_jet_spectrum_weighted_y%i_cent%i",i,j),Form("h_reco_jet_spectrum_weighted_y%i_cent%i",i,j),ptJetBinsN, ptJetBins);
			//h_reco_jet_spectrum_weighted.at(i).at(j) = temphist_1D;
			//temphist_1D = new TH1D(Form("h_true_jet_spectrum_weighted_y%i_cent%i",i,j),Form("h_true_jet_spectrum_weighted_y%i_cent%i",i,j),ptJetBinsN, ptJetBins);
			//h_true_jet_spectrum_weighted.at(i).at(j) = temphist_1D;

			h_reco_jet_spectrum.at(i).at(j)->Sumw2();
			//h_reco_jet_spectrum_weighted.at(i).at(j)->Sumw2();
			//h_true_jet_spectrum_weighted.at(i).at(j)->Sumw2();


			wk()->addOutput (h_reco_jet_spectrum.at(i).at(j));
			//wk()->addOutput (h_reco_jet_spectrum_weighted.at(i).at(j));
			//wk()->addOutput (h_true_jet_spectrum_weighted.at(i).at(j));


			//MC only
			if (_data_switch==1)
			{

				//Truth
				temphist_1D = new TH1D(Form("h_true_jet_spectrum_y%i_cent%i",i,j),Form("h_true_jet_spectrum_y%i_cent%i",i,j),ptJetBinsN, ptJetBins);
				h_true_jet_spectrum.at(i).at(j) = temphist_1D;

				temphist_1D = new TH1D(Form("h_reco_jet_spectrum_matched_y%i_cent%i",i,j),Form("h_reco_jet_spectrum_matched_y%i_cent%i",i,j),ptJetBinsN, ptJetBins);
				h_reco_jet_spectrum_matched.at(i).at(j) = temphist_1D;

				temphist_1D = new TH1D(Form("h_true_jet_spectrum_matched_y%i_cent%i",i,j),Form("h_true_jet_spectrum_matched_y%i_cent%i",i,j),ptJetBinsN, ptJetBins);
				h_true_jet_spectrum_matched.at(i).at(j) = temphist_1D;

				temphist_2D = new TH2D(Form("ff_jetResponse_y%i_cent%i",i,j),Form("ff_jetResponse_y%i_cent%i",i,j),ptJetBinsN, ptJetBins, ptJetBinsN, ptJetBins);
				ff_jetResponse.at(i).at(j) = temphist_2D;

				response_jet.at(i).at(j) = new RooUnfoldResponse(h_reco_jet_spectrum.at(i).at(j),h_true_jet_spectrum.at(i).at(j));


				//Responses
				h_true_jet_spectrum.at(i).at(j)->Sumw2();
				h_reco_jet_spectrum_matched.at(i).at(j)->Sumw2();
				h_true_jet_spectrum_matched.at(i).at(j)->Sumw2();

				ff_jetResponse.at(i).at(j)->Sumw2();

				wk()->addOutput (h_true_jet_spectrum.at(i).at(j));
				wk()->addOutput (h_true_jet_spectrum_matched.at(i).at(j));
				wk()->addOutput (h_reco_jet_spectrum_matched.at(i).at(j));
				wk()->addOutput (ff_jetResponse.at(i).at(j));
				wk()->addOutput (response_jet.at(i).at(j));
			}

		}
		
//		temphist_1D = new TH1D(Form("h_reco_jet_spectrum_fine_cent%i",j),Form("h_reco_jet_spectrum_fine_cent%i",j),50,40,140);
//		h_reco_jet_spectrum_fine.push_back(temphist_1D);
//		wk()->addOutput (h_reco_jet_spectrum_fine.at(j));

//		if (_data_switch==1){
//			temphist_1D = new TH1D(Form("h_true_jet_spectrum_fine_cent%i",j),Form("h_true_jet_spectrum_fine_cent%i",j),50,40,140);
//			h_true_jet_spectrum_fine.push_back(temphist_1D);
//			wk()->addOutput (h_true_jet_spectrum_fine.at(j));
//		}
	}

	cout << " Histograms ready" << endl;


	return EL::StatusCode::SUCCESS;
}
