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


	int ptJetBinsN, etaJetBinsN, phiJetBinsN, ptTrkBinsN, etaTrkBinsN, phiTrkBinsN, zBinsN, zBinsFineN, d0z0BinsN, respBinsN, finehitsBinsN, dR_resBinsN, PsiBinsN;
	double ptJetBins[1000], etaJetBins[1000], phiJetBins[1000], ptTrkBins[1000], etaTrkBins[1000], phiTrkBins[1000], zBins[1000], zBinsFine[1000],d0z0Bins[1000], respBins[1000], finehitsBins[1000], dR_resBins[1000], PsiBins[1000];

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
	SetupBinning(0, "PsiBins", PsiBins, PsiBinsN);

	Double_t PVBins[3]={0,1,2};
	int PVBinsN=2;

	h_fcal_change = new TH2D("h_fcal_change","h_fcal_change",10,0,10,10,0,10);
	wk()->addOutput (h_fcal_change);

	h_fcal_diff = new TH1D("h_fcal_diff","h_fcal_diff",400,-0.1,0.9);
	wk()->addOutput (h_fcal_diff);
	
	//Axis histograms
	h_jet_pt_eta_phi = new TH3D("h_jet_pt_eta_phi","h_jet_pt_eta_phi",ptJetBinsN, ptJetBins, etaJetBinsN, etaJetBins, phiJetBinsN, phiJetBins);
	h_trk_pt_eta_phi = new TH3D("h_trk_pt_eta_phi","h_trk_pt_eta_phi",ptTrkBinsN, ptTrkBins, etaTrkBinsN, etaTrkBins, phiTrkBinsN, phiTrkBins);
	wk()->addOutput (h_jet_pt_eta_phi);
	wk()->addOutput (h_trk_pt_eta_phi);


	double tmp_finerBinsN = 300, tmp_rBinsN = 13, tmp_coneBinsN = 51;
	double tmp_finerBins[400], tmp_rBins[200], tmp_coneBins[50];

	tmp_rBins[0] = 0.; tmp_rBins[1] = 0.05; tmp_rBins[2] = 0.1; tmp_rBins[3] = 0.15; tmp_rBins[4] = 0.2; tmp_rBins[5] = 0.25; tmp_rBins[6] = 0.3; tmp_rBins[7] = 0.4; tmp_rBins[8] = 0.5; tmp_rBins[9] = 0.6; tmp_rBins[10] = 0.7; tmp_rBins[11] = 0.8; tmp_rBins[12] = 1.0; tmp_rBins[13] = 1.2;

	for (int i = 0; i <= tmp_finerBinsN; i++)
	{
		tmp_finerBins[i] = -1.5+i*0.01;
		cout << tmp_finerBins[i] << ", ";
	}

	cout << "cone bins" << endl;
	for (int i = 0; i <= tmp_coneBinsN; i++)
	{
		tmp_coneBins[i] = -0.5+i;
		cout << tmp_coneBins[i] << ", ";
	}

	h_tmp_trk = new TH3D("h_tmp_trk","h_tmp_trk",ptTrkBinsN, ptTrkBins, etaTrkBinsN, etaTrkBins, phiTrkBinsN, phiTrkBins);
	wk()->addOutput(h_tmp_trk);
	h_tmp_dR = new TH1D("h_tmp_dR","h_tmp_dR",120,0,1.2);
	wk()->addOutput(h_tmp_dR);
	h_tmp_coneIndex = new TH1D("h_tmp_coneIndex","h_tmp_coneIndex",100,0,100);
	wk()->addOutput(h_tmp_coneIndex);
	h_tmp_dRBin = new TH1D("h_tmp_dRBin","h_tmp_dRBin",20,0,20);
	wk()->addOutput(h_tmp_dRBin);
	h_tmp_rdEtadPhi = new TH3D("h_tmp_rdEtadPhi","h_tmp_rdEtadPhi;r;#delta#eta;#delta#phi",tmp_rBinsN, tmp_rBins, tmp_finerBinsN, tmp_finerBins, tmp_finerBinsN, tmp_finerBins);
	wk()->addOutput(h_tmp_rdEtadPhi);
	h_tmp_cone_stats = new TH1D("h_tmp_cone_stats","h_tmp_cone_stats;# Valid Cone ; Events;",51,-0.5,50.5);
	wk()->addOutput(h_tmp_cone_stats);
//	h_tmp_cone_stats = new TH3D("h_tmp_cone_stats","h_tmp_cone_stats# Valid Cone;p_T^{Trk};p_T^{jet}",tmp_coneBinsN, tmp_coneBins, ptTrkBinsN, ptTrkBins, ptJetBinsN, ptJetBins);
//	wk()->addOutput(h_tmp_cone_stats);


	h_cone_map = new TH2D("h_cone_map","h_cone_map;#eta_{cone};#phi_{cone}",100, -2.5, 2.5, 140, -3.5, 3.5);
	wk()->addOutput(h_cone_map);

	//Debugging histograms
	for (int i=0;i<GetCentralityNBins(_centrality_scheme);i++)
	{
		temphist_3D = new TH3D(Form("h_reco_truth_matched_cent%i",i),Form("h_reco_truth_matched_cent%i",i),ptJetBinsN, ptJetBins, ptJetBinsN, ptJetBins, dR_resBinsN, dR_resBins);
		temphist_3D->Sumw2();
		h_reco_truth_matched.push_back(temphist_3D);

		temphist_3D = new TH3D(Form("h_jet_for_eff_cent%i",i),Form("h_jet_for_eff_cent%i",i), ptJetBinsN, ptJetBins, etaJetBinsN, etaJetBins, phiJetBinsN, phiJetBins);
		temphist_3D->Sumw2();
		h_jet_for_eff.push_back(temphist_3D);

		temphist_3D = new TH3D(Form("h_jet_for_eff_full_cent%i",i),Form("h_jet_for_eff_full_cent%i",i), ptJetBinsN, ptJetBins, etaJetBinsN, etaJetBins, phiJetBinsN, phiJetBins);
		temphist_3D->Sumw2();
		h_jet_for_eff_full.push_back(temphist_3D);

		temphist_3D = new TH3D(Form("h_jet_psi3_cent%i",i),Form("h_jet_psi3_cent%i",i), ptJetBinsN, ptJetBins, etaJetBinsN, etaJetBins, phiJetBinsN, phiJetBins);
		temphist_3D->Sumw2();
		h_jet_psi3.push_back(temphist_3D);

		temphist_3D = new TH3D(Form("h_reco_jets_cent%i",i),Form("h_reco_jets_cent%i",i), ptJetBinsN, ptJetBins, etaTrkBinsN, etaTrkBins, phiJetBinsN, phiJetBins);
		temphist_3D->Sumw2();
		h_reco_jets.push_back(temphist_3D);

		wk()->addOutput (h_reco_truth_matched.at(i));
		wk()->addOutput (h_jet_for_eff.at(i));
		wk()->addOutput (h_jet_for_eff_full.at(i));
		wk()->addOutput (h_jet_psi3.at(i));
		wk()->addOutput (h_reco_jets.at(i));
	}

	//Basic histograms
	h_FCal_Et = new TH1D("h_FCal_Et",";FCal E_{T};N",1200,0.,6.);
	h_FCal_Et->Sumw2();

	h_fcal_mc = new TH1D("h_fcal_mc",";FCal E_{T};N",1200,0.,6.);
	h_fcal_mc->Sumw2();

	h_fcal_mbov = new TH1D("h_fcal_mbov",";FCal E_{T};N",1200,0.,6.);
	h_fcal_mbov->Sumw2();

	h_FCal_Et_unw = new TH1D("h_FCal_Et_unw",";FCal E_{T};N",1200,0.,6.);
	h_FCal_Et_unw->Sumw2();

	h_FCal_Et_restr = new TH1D("h_FCal_Et_restr",";FCal E_{T};N",1200,0.,6.);
	h_FCal_Et_restr->Sumw2();

	h_RejectionHisto = new TH1D("RejectionHisto","RejectionHisto",9,0,9);
	SetRejectionHistogram(h_RejectionHisto);

	h_centrality = new TH1D("Centrality","Centrality",10,0,10);
	h_centrality->Sumw2();

	h_cent_mc = new TH1D("h_cent_mc","h_cent_mc",10,0,10);
	h_cent_mc->Sumw2();
	
	h_cent_mbov = new TH1D("h_cent_mbov","h_cent_mbov",10,0,10);
	h_cent_mbov->Sumw2();

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
	wk()->addOutput (h_fcal_mc);
	wk()->addOutput (h_fcal_mbov);
	wk()->addOutput (h_cent_mc);
	wk()->addOutput (h_cent_mbov);

	wk()->addOutput (h_FCal_Et_unw);
	wk()->addOutput (h_FCal_Et_restr);
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

//		temphist_3D = new TH3D(Form("h_jetpT_v_multiplicity_cent%i",i),"h_jetpT_v_multiplicity;jet p_{T} GeV;p_{T}^{trk,min} [GeV];multiplicity",20,0,1000,6,0,6,50,0,50);
//		h_jetpT_v_multiplicity.push_back(temphist_3D);
//		h_jetpT_v_multiplicity.at(i)->Sumw2();


		wk()->addOutput (h_PixHits.at(i));
		wk()->addOutput (h_SCTHits.at(i));
		wk()->addOutput (h_d0.at(i));
		wk()->addOutput (h_z0sintheta.at(i));
//		wk()->addOutput (h_jetpT_v_multiplicity.at(i));

	}


	ChPS_raw =  vector<vector<TH2D*> > (_ndRBins, vector<TH2D*>(_nCentbins));

	ChPS_cone_UE =  vector<vector<TH2D*> > (_ndRBins, vector<TH2D*>(_nCentbins));
	ChPS_MB_UE =  vector<vector<TH2D*> > (_ndRBins, vector<TH2D*>(_nCentbins));
	ChPS_MB_UE_err =  vector<vector<TH2D*> > (_ndRBins, vector<TH2D*>(_nCentbins));
	ChPS_FS_UE =  vector<vector<TH2D*> > (_ndRBins, vector<TH2D*>(_nCentbins));
	ChPS_FNS_UE =  vector<vector<TH2D*> > (_ndRBins, vector<TH2D*>(_nCentbins));
	h_dR_change =  vector<vector<TH3D*> > (ptJetBinsN, vector<TH3D*>(_nCentbins));

	h_reco_jet_spectrum =  vector<vector<TH1D*> > (_nJetYBins, vector<TH1D*>(_nCentbins));

	//in MC only
	if (_data_switch==1)
	{
		ChPS_TM_UE =  vector<vector<TH2D*> > (_ndRBins, vector<TH2D*>(_nCentbins));
		ChPS_truth =  vector<vector<TH2D*> > (_ndRBins, vector<TH2D*>(_nCentbins));

		ff_trackpTResponse =  vector<vector<TH3D*> > (_ndRBins, vector<TH3D*>(_nCentbins));
		response_ChPS =  vector<vector<RooUnfoldResponse*> > (_ndRBins, vector<RooUnfoldResponse*>(_nCentbins));

		h_true_jet_spectrum =  vector<vector<TH1D*> > (_nJetYBins, vector<TH1D*>(_nCentbins));
		ff_jetResponse =  vector<vector<TH2D*> > (_nJetYBins, vector<TH2D*>(_nCentbins));
		response_jet =  vector<vector<RooUnfoldResponse*> > (_nJetYBins, vector<RooUnfoldResponse*>(_nCentbins));
	}

	//test UE histograms
	cone_norm_jet = vector<TH1D*> (_nCentbins);
	MB_norm_jet = vector<TH1D*> (_nCentbins);
	TM_norm_jet = vector<TH1D*> (_nCentbins);
	FS_norm_jet = vector<TH1D*> (_nCentbins);

	h_UE_dNdEtadPhidpT =  vector<vector<vector<vector<TH3D*>>>> (ptJetBinsN, vector<vector<vector<TH3D*>>> (PsiBinsN, vector<vector<TH3D*>> (_nCentbins, vector<TH3D*>(_ndRBins))));
	h_jet_v_Psi =  vector<vector<TH3D*>> (ptJetBinsN, vector<TH3D*>(_nCentbins));

	for (int j=0;j<_nCentbins;j++)
	{
		temphist_1D = new TH1D(Form("cone_norm_jet_cent%i",j),Form("cone_norm_jet_cent%i",j),ptJetBinsN, ptJetBins);
		cone_norm_jet.at(j) = temphist_1D;
		wk()->addOutput (cone_norm_jet.at(j));

		temphist_1D = new TH1D(Form("MB_norm_jet_cent%i",j),Form("MB_norm_jet_cent%i",j),ptJetBinsN, ptJetBins);
		MB_norm_jet.at(j) = temphist_1D;
		wk()->addOutput (MB_norm_jet.at(j));

		temphist_1D = new TH1D(Form("TM_norm_jet_cent%i",j),Form("TM_norm_jet_cent%i",j),ptJetBinsN, ptJetBins);
		TM_norm_jet.at(j) = temphist_1D;
		wk()->addOutput (TM_norm_jet.at(j));

		temphist_1D = new TH1D(Form("FS_norm_jet_cent%i",j),Form("FS_norm_jet_cent%i",j),ptJetBinsN, ptJetBins);
		FS_norm_jet.at(j) = temphist_1D;
		wk()->addOutput (FS_norm_jet.at(j));


		for (int i=0;i<ptJetBinsN;i++)
		{
			if (j == _nCentbins - 1) continue;

			temphist_3D = new TH3D(Form("h_jet_v_Psi_cent%i_jetpt%i",j, i),Form("h_jet_v_Psi_cent%i_jetpt%i",j, i),PsiBinsN,PsiBins,etaTrkBinsN,etaTrkBins,phiTrkBinsN,phiTrkBins);
			h_jet_v_Psi.at(i).at(j) = temphist_3D;
			h_jet_v_Psi.at(i).at(j)->Sumw2();
			wk()->addOutput (h_jet_v_Psi.at(i).at(j));

			if (derive_UE_mode)
			{
				//restrict to trk < 10 GeV, jet > 100, jet < 400, cent != 6
				if (i >= lo_jetpt_bin && i <= hi_jetpt_bin)
				{
					for (int m=0;m<10;m++) //psibins
					{
						for (int k=0;k<_ndRBins;k++)
						{
							temphist_3D = new TH3D(Form("h_UE_jetpt%i_dPsi%i_cent%i_dR%i",i,m,j,k),Form("h_UE_jetpt%i_dPsi%i_cent%i_dR%i",i,m,j,k),ptTrkBinsN, ptTrkBins,etaTrkBinsN,etaTrkBins,phiTrkBinsN,phiTrkBins);
							h_UE_dNdEtadPhidpT.at(i).at(m).at(j).at(k) = temphist_3D;
							h_UE_dNdEtadPhidpT.at(i).at(m).at(j).at(k)->Sumw2();
							wk()->addOutput (h_UE_dNdEtadPhidpT.at(i).at(m).at(j).at(k));
						}
					}
				}
			}
		}


	}

	for (int j=0;j<_nCentbins;j++)
	{

		for (int i=0;i<ptJetBinsN;i++)
		{
			temphist_3D = new TH3D(Form("h_dR_change_jetpt%i_cent%i",i,j),Form("h_dR_change_jetpt%i_cent%i",i,j),_ndRBins-1,trkcorr->dRrange,_ndRBins-1,trkcorr->dRrange,ptTrkBinsN, ptTrkBins);
			h_dR_change.at(i).at(j) = temphist_3D;
			h_dR_change.at(i).at(j)->Sumw2();
			wk()->addOutput (h_dR_change.at(i).at(j));
		}


		for (int i=0;i<_ndRBins;i++)
		{
			//D(Pt)
			temphist_2D = new TH2D(Form("ChPS_raw_0_dR%i_cent%i",i,j),Form("ChPS_raw_0_dR%i_cent%i",i,j),ptTrkBinsN, ptTrkBins, ptJetBinsN, ptJetBins);
			ChPS_raw.at(i).at(j) = temphist_2D;

			temphist_2D = new TH2D(Form("ChPS_cone_UE_dR%i_cent%i",i,j),Form("ChPS_cone_UE_dR%i_cent%i",i,j),ptTrkBinsN, ptTrkBins, ptJetBinsN, ptJetBins);
			ChPS_cone_UE.at(i).at(j) = temphist_2D;

			temphist_2D = new TH2D(Form("ChPS_MB_UE_dR%i_cent%i",i,j),Form("ChPS_MB_UE_dR%i_cent%i",i,j),ptTrkBinsN, ptTrkBins, ptJetBinsN, ptJetBins);
			ChPS_MB_UE.at(i).at(j) = temphist_2D;

			temphist_2D = new TH2D(Form("ChPS_MB_UE_err_dR%i_cent%i",i,j),Form("ChPS_MB_UE_err_dR%i_cent%i",i,j),ptTrkBinsN, ptTrkBins, ptJetBinsN, ptJetBins);
			ChPS_MB_UE_err.at(i).at(j) = temphist_2D;

			temphist_2D = new TH2D(Form("ChPS_FS_UE_dR%i_cent%i",i,j),Form("ChPS_FS_UE_dR%i_cent%i",i,j),ptTrkBinsN, ptTrkBins, ptJetBinsN, ptJetBins);
			ChPS_FS_UE.at(i).at(j) = temphist_2D;

			temphist_2D = new TH2D(Form("ChPS_FNS_UE_dR%i_cent%i",i,j),Form("ChPS_FNS_UE_dR%i_cent%i",i,j),ptTrkBinsN, ptTrkBins, ptJetBinsN, ptJetBins);
			ChPS_FNS_UE.at(i).at(j) = temphist_2D;


			ChPS_raw.at(i).at(j)->Sumw2();
			ChPS_cone_UE.at(i).at(j)->Sumw2();
			ChPS_MB_UE.at(i).at(j)->Sumw2();
			ChPS_MB_UE_err.at(i).at(j)->Sumw2();
			ChPS_FS_UE.at(i).at(j)->Sumw2();
			ChPS_FNS_UE.at(i).at(j)->Sumw2();

			wk()->addOutput (ChPS_raw.at(i).at(j));
			wk()->addOutput (ChPS_cone_UE.at(i).at(j));
			wk()->addOutput (ChPS_MB_UE.at(i).at(j));
			wk()->addOutput (ChPS_MB_UE_err.at(i).at(j));
			wk()->addOutput (ChPS_FS_UE.at(i).at(j));
			wk()->addOutput (ChPS_FNS_UE.at(i).at(j));

			if (_data_switch==1)
			{
				temphist_2D = new TH2D(Form("ChPS_TM_UE_dR%i_cent%i",i,j),Form("ChPS_TM_UE_dR%i_cent%i",i,j),ptTrkBinsN, ptTrkBins, ptJetBinsN, ptJetBins);
				ChPS_TM_UE.at(i).at(j) = temphist_2D;

				temphist_2D = new TH2D(Form("ChPS_truth_dR%i_cent%i",i,j),Form("ChPS_truth_dR%i_cent%i",i,j),ptTrkBinsN, ptTrkBins, ptJetBinsN, ptJetBins);
				ChPS_truth.at(i).at(j) = temphist_2D;

				//Truth
				temphist_3D = new TH3D(Form("ff_trackpTResponse_dR%i_cent%i",i,j),Form("ff_trackpTResponse_dR%i_cent%i",i,j),ptTrkBinsN, ptTrkBins,ptTrkBinsN, ptTrkBins, ptJetBinsN, ptJetBins);
				ff_trackpTResponse.at(i).at(j) = temphist_3D;

				//Responses
				response_ChPS.at(i).at(j) = new RooUnfoldResponse(ChPS_raw.at(i).at(j),ChPS_truth.at(i).at(j));

				//Truth
				ChPS_TM_UE.at(i).at(j)->Sumw2();
				ChPS_truth.at(i).at(j)->Sumw2();
				ff_trackpTResponse.at(i).at(j)->Sumw2();

				wk()->addOutput (ChPS_TM_UE.at(i).at(j));
				wk()->addOutput (ChPS_truth.at(i).at(j));
				wk()->addOutput (ff_trackpTResponse.at(i).at(j));
				wk()->addOutput (response_ChPS.at(i).at(j));
			}

		}


		for (int i=0;i<_nJetYBins;i++)
		{

			//Jet spectra
			temphist_1D = new TH1D(Form("h_reco_jet_spectrum_y%i_cent%i",i,j),Form("h_reco_jet_spectrum_y%i_cent%i",i,j),ptJetBinsN, ptJetBins);
			h_reco_jet_spectrum.at(i).at(j) = temphist_1D;
			h_reco_jet_spectrum.at(i).at(j)->Sumw2();
			wk()->addOutput (h_reco_jet_spectrum.at(i).at(j));

			//MC only
			if (_data_switch==1)
			{
				//Truth
				temphist_1D = new TH1D(Form("h_true_jet_spectrum_y%i_cent%i",i,j),Form("h_true_jet_spectrum_y%i_cent%i",i,j),ptJetBinsN, ptJetBins);
				h_true_jet_spectrum.at(i).at(j) = temphist_1D;

				temphist_2D = new TH2D(Form("ff_jetResponse_y%i_cent%i",i,j),Form("ff_jetResponse_y%i_cent%i",i,j),ptJetBinsN, ptJetBins, ptJetBinsN, ptJetBins);
				ff_jetResponse.at(i).at(j) = temphist_2D;
				response_jet.at(i).at(j) = new RooUnfoldResponse(h_reco_jet_spectrum.at(i).at(j),h_true_jet_spectrum.at(i).at(j));

				//Responses
				h_true_jet_spectrum.at(i).at(j)->Sumw2();
				ff_jetResponse.at(i).at(j)->Sumw2();

				wk()->addOutput (h_true_jet_spectrum.at(i).at(j));
				wk()->addOutput (ff_jetResponse.at(i).at(j));
				wk()->addOutput (response_jet.at(i).at(j));
			}
		}

	}

	cout << " Histograms ready" << endl;


	return EL::StatusCode::SUCCESS;
}
