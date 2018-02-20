#define HistoMaker_cxx
#include "pPbFragmentation/PbPbFragmentation.h"
#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>

using namespace std;
using namespace MTCorrector;

EL::StatusCode PbPbFragmentation :: histInitialize ()
{
	cout << " Setting  histograms" << endl;
	
	int ptJetBinsN, etaJetBinsN, phiJetBinsN, ptTrkBinsN, ptTrkBinsFineN, ptTrkBinsSumN , etaTrkBinsN, phiTrkBinsN, zBinsN, zBinsFineN, d0z0BinsN, respBinsN, finehitsBinsN, RdRBinsN, RunBinsN;
	double ptJetBins[1000], etaJetBins[1000], phiJetBins[1000], ptTrkBins[1000], ptTrkBinsFine[1000], ptTrkBinsSum[1000], RdRBins[20], etaTrkBins[1000], phiTrkBins[1000], zBins[1000], zBinsFine[1000],d0z0Bins[1000], respBins[1000], RunBins[50], finehitsBins[1000];

	SetupBinning(0, "pt-jet-PbPb", ptJetBins, ptJetBinsN);
	SetupBinning(0, "eta-jet", etaJetBins, etaJetBinsN);
	SetupBinning(0, "phi-trk", phiJetBins, phiJetBinsN);
	SetupBinning(0, "pt-trk-rebin", ptTrkBins, ptTrkBinsN);
	SetupBinning(0, "pt-trk", ptTrkBinsFine, ptTrkBinsFineN);
	SetupBinning(0, "pt-trk-sum", ptTrkBinsSum, ptTrkBinsSumN);
	SetupBinning(0, "eta-trk", etaTrkBins, etaTrkBinsN);
	SetupBinning(0, "phi-trk", phiTrkBins, phiTrkBinsN);
	SetupBinning(0, "z-rebin", zBins, zBinsN);
	SetupBinning(0, "z", zBinsFine, zBinsFineN);
	SetupBinning(0, "d0z0", d0z0Bins, d0z0BinsN);
	SetupBinning(0, "resp", respBins, respBinsN);
	SetupBinning(0, "hits_fine", finehitsBins, finehitsBinsN);
	SetupBinning(0, "dR-RdR", RdRBins, RdRBinsN);
	SetupBinning(0, "PbPb_runs", RunBins, RunBinsN);
	
	//Track corrector
	int _nCentbins = GetCentralityNBins(_centrality_scheme);
	trkcorr = new TrackCorrector(_cut_level.c_str(),GetCentralityNBins(31)-1,_eff_jety);	for (int i = 0; i < ptTrkBinsN-1; i++){ 
		if (ptTrkBins[i]<_pTtrkCut && _pTtrkCut < ptTrkBins[i+1]) {trkcorr->trkpTThreshold = ptTrkBins[i]; break;}
	}
	
	
	//Jet corrector
	jetcorr = new JetCorrector();	
	int _nJetYBins = jetcorr->nJetYBins;
	
	jetcorr->JERcut = _JERBalancecut;
	//jetcorr->min_jet_pt = _pTjetCut+1;
	jetcorr->min_jet_pt = 100.; //TODO
	jetcorr->max_jet_pt = 630.944; //Maximum jet pt for reweighting 
	jetcorr->m_isMB = _isMB;
	

	

	Double_t PVBins[3]={0,1,2};
	int PVBinsN=2;

	//Basic histograms

	h_FCal_Et = new TH1D("h_FCal_Et",";FCal E_{T};N",140,0,7);
	h_FCal_Et->Sumw2();

	h_RejectionHisto = new TH1D("RejectionHisto","RejectionHisto",9,0,9);
	SetRejectionHistogram(h_RejectionHisto);
	
	h_centrality = new TH1D("Centrality","Centrality",10,0,10);
	h_centrality->Sumw2();
	
	h_Njets_v_Run = new TH1D("Njets_v_Run","Njets_v_Run",RunBinsN,0,RunBinsN);
	h_Njets_v_Run->Sumw2();
	
	h_FF_v_Run = new TH1D("FF_v_Run","FF_v_Run",RunBinsN,0,RunBinsN);
	h_FF_v_Run->Sumw2();

	for (int i = 0; i<RunBinsN; i++){
		h_Njets_v_Run->GetXaxis()->SetBinLabel(i+1,Form("%i",(int)RunBins[i]));
		h_FF_v_Run->GetXaxis()->SetBinLabel(i+1,Form("%i",(int)RunBins[i]));
	}

	hET_ETsub = new TH3D("hET_ETsub","hET_ETsub",200,0,200,100,-5,5,200,-50,50);
	hET_ETsub->Sumw2();
	
	
	//deriv_val = new TH3F("deriv_val","deriv_val;p_{T} [GeV];RunNumber;lbn",10,0,10,1225,286710,287935,200,0,1000);
	//deriv_val->Sumw2();

	h_triggercounter = new TH2D("h_triggercounter","h_triggercounter",_nTriggers,0,_nTriggers,2,-0.5,1.5);
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
	wk()->addOutput (h_Njets_v_Run);
	wk()->addOutput (h_FF_v_Run);
	//wk()->addOutput (deriv_val);
	
	//Jet helper histogram
	h_event_jet_counts = new TH1D("h_event_jet_counts","h_event_jet_counts",ptJetBinsN, ptJetBins);
	h_event_jet_counts->Sumw2();
	
	TH3D* temphist_3D = nullptr;
	TH2D* temphist_2D = nullptr;
	TH1D* temphist_1D = nullptr;
	

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
		h_PixHits.push_back(temphist_3D);

		temphist_3D = new TH3D(Form("h_SCTHits_cent%i",i),Form("h_SCTHits_cent%i",i),ptTrkBinsN, ptTrkBins,etaTrkBinsN, etaTrkBins,finehitsBinsN,finehitsBins);
		h_SCTHits.push_back(temphist_3D);

		temphist_2D = new TH2D(Form("h_d0_cent%i",i),Form("h_d0_cent%i",i),ptTrkBinsN,ptTrkBins,d0z0BinsN,d0z0Bins);
		h_d0.push_back(temphist_2D);

		temphist_2D = new TH2D(Form("h_z0sintheta_cent%i",i),Form("h_z0sintheta_cent%i",i),ptTrkBinsN, ptTrkBins,d0z0BinsN,d0z0Bins);
		h_z0sintheta.push_back(temphist_2D);

		//temphist_3D = new TH3D(Form("h_vx_cent%i",i),"h_vx;v_{x};v_{y};v_{z}",200,-2,2,200,-2,2,200,-200,200);
		//h_vx.push_back(temphist_3D);
		
		temphist_3D = new TH3D(Form("h_jetpT_v_multiplicity_cent%i",i),"h_jetpT_v_multiplicity;jet p_{T} GeV;p_{T}^{trk,min} [GeV];multiplicity",20,0,1000,6,0,6,50,0,50);
		h_jetpT_v_multiplicity.push_back(temphist_3D);
		h_jetpT_v_multiplicity.at(i)->Sumw2();

		//temphist_3D = new TH3D(Form("h_BL_cent%i",i),"h_BL;BL_{x};BL_{y};BL_{z}",200,-1,1,200,-1,1,200,-10,10);
		//h_BL.push_back(temphist_3D);
		
		temphist_2D = new TH2D(Form("h_R2vR4_cent%i",i),Form("h_R2vR4_cent%i",i),40,100,500,80,-0.4,0.4);
		h_R2vR4.push_back(temphist_2D);
		
		temphist_3D = new TH3D(Form("h_RdR_cent%i",i),Form("h_RdR_cent%i",i),ptJetBinsN, ptJetBins,ptJetBinsN, ptJetBins,RdRBinsN, RdRBins);
		h_RdR.push_back(temphist_3D);
		
		temphist_2D = new TH2D(Form("ff_zg1_v_PDG_cent%i",i),Form("ff_zg1_v_PDG_cent%i",i),20,1.,2.,6000,0,6000);
		ff_zg1_v_PDG.push_back(temphist_2D);
		
		if (_data_switch==1){
			temphist_3D = new TH3D(Form("JES_v_max_z_cent%i",i),Form("JES_v_max_z_cent%i",i),zBinsN, zBins,respBinsN,respBins, ptJetBinsN, ptJetBins);
			JES_v_max_z.push_back(temphist_3D);
			temphist_3D = new TH3D(Form("JES_v_max_pT_cent%i",i),Form("JES_v_max_pT_cent%i",i),ptTrkBinsN, ptTrkBins, respBinsN,respBins, ptJetBinsN, ptJetBins);
			JES_v_max_pT.push_back(temphist_3D);
		}
		
		h_RdR.at(i)->Sumw2();
		
		wk()->addOutput (h_PixHits.at(i));
		wk()->addOutput (h_SCTHits.at(i));
		wk()->addOutput (h_d0.at(i));
		wk()->addOutput (h_z0sintheta.at(i));
		//wk()->addOutput (h_vx.at(i));
		
		wk()->addOutput (h_R2vR4.at(i));
		wk()->addOutput (h_RdR.at(i));
		
		wk()->addOutput (h_jetpT_v_multiplicity.at(i));
		
		if (_data_switch==1){
			wk()->addOutput (JES_v_max_pT.at(i));
			wk()->addOutput (JES_v_max_z.at(i));
			wk()->addOutput (ff_zg1_v_PDG.at(i));
		}	
		//wk()->addOutput (h_BL.at(i));

	}
	
	ff_raw =  vector<vector<TH2D*> > (_nJetYBins, vector<TH2D*>(_nCentbins));
	ff_raw_fine =  vector<vector<TH2D*> > (_nJetYBins, vector<TH2D*>(_nCentbins));
    ff_raw_UE =  vector<vector<TH2D*> > (_nJetYBins, vector<TH2D*>(_nCentbins));
    ff_raw_UE_fine =  vector<vector<TH2D*> > (_nJetYBins, vector<TH2D*>(_nCentbins));    
    ChPS_raw =  vector<vector<TH2D*> > (_nJetYBins, vector<TH2D*>(_nCentbins));
    ChPS_raw_fine =  vector<vector<TH2D*> > (_nJetYBins, vector<TH2D*>(_nCentbins));
    ChPS_raw_UE =  vector<vector<TH2D*> > (_nJetYBins, vector<TH2D*>(_nCentbins));
    ChPS_raw_UE_fine =  vector<vector<TH2D*> > (_nJetYBins, vector<TH2D*>(_nCentbins));
    ChPS_raw_UE_fine_response =  vector<vector<TH3D*> > (_nJetYBins, vector<TH3D*>(_nCentbins));
    h_reco_jet_spectrum =  vector<vector<TH1D*> > (_nJetYBins, vector<TH1D*>(_nCentbins));
    h_reco_jet_spectrum_counts =  vector<vector<TH1D*> > (_nJetYBins, vector<TH1D*>(_nCentbins));
    h_reco_jet_spectrum_counts_TMR =  vector<vector<TH1D*> > (_nJetYBins, vector<TH1D*>(_nCentbins));
    
    //in MC only
    if (_data_switch==1){
    	ff_UE_z =  vector<vector<TH2D*> > (_nJetYBins, vector<TH2D*>(_nCentbins));
    	ff_UE_pT =  vector<vector<TH2D*> > (_nJetYBins, vector<TH2D*>(_nCentbins));
    	ff_truth =  vector<vector<TH2D*> > (_nJetYBins, vector<TH2D*>(_nCentbins));
    	ChPS_truth =  vector<vector<TH2D*> > (_nJetYBins, vector<TH2D*>(_nCentbins));
        UE_pT_Mcorrelation_v_response_fine =  vector<vector<TH3D*> > (_nJetYBins, vector<TH3D*>(_nCentbins));
    	
    	ff_UE_z_fine =  vector<vector<TH2D*> > (_nJetYBins, vector<TH2D*>(_nCentbins));
    	ff_UE_pT_fine =  vector<vector<TH2D*> > (_nJetYBins, vector<TH2D*>(_nCentbins));
    	ff_UE_pT_fine_response =  vector<vector<TH3D*> > (_nJetYBins, vector<TH3D*>(_nCentbins));
    	ff_truth_fine =  vector<vector<TH2D*> > (_nJetYBins, vector<TH2D*>(_nCentbins));
    	//ff_truth_alt_fine =  vector<vector<TH2D*> > (_nJetYBins, vector<TH2D*>(_nCentbins));
    	ChPS_truth_fine =  vector<vector<TH2D*> > (_nJetYBins, vector<TH2D*>(_nCentbins));
    	
    	
    	h_true_jet_spectrum =  vector<vector<TH1D*> > (_nJetYBins, vector<TH1D*>(_nCentbins));
    	
    	ff_jetResponse =  vector<vector<TH2D*> > (_nJetYBins, vector<TH2D*>(_nCentbins));
    	ff_trackpTResponse =  vector<vector<TH3D*> > (_nJetYBins, vector<TH3D*>(_nCentbins));
    	ff_trackzResponse =  vector<vector<TH3D*> > (_nJetYBins, vector<TH3D*>(_nCentbins));
    	ff_trackzResponse_counts =  vector<vector<TH3D*> > (_nJetYBins, vector<TH3D*>(_nCentbins));
    	response_FF =  vector<vector<RooUnfoldResponse*> > (_nJetYBins, vector<RooUnfoldResponse*>(_nCentbins));
    	response_ChPS =  vector<vector<RooUnfoldResponse*> > (_nJetYBins, vector<RooUnfoldResponse*>(_nCentbins));
    	response_FF_fine =  vector<vector<RooUnfoldResponse*> > (_nJetYBins, vector<RooUnfoldResponse*>(_nCentbins));
    	response_ChPS_fine =  vector<vector<RooUnfoldResponse*> > (_nJetYBins, vector<RooUnfoldResponse*>(_nCentbins));
    	
    	response_jet =  vector<vector<RooUnfoldResponse*> > (_nJetYBins, vector<RooUnfoldResponse*>(_nCentbins));
    	
    	response_FF_counts =  vector<vector<RooUnfoldResponse*> > (_nJetYBins, vector<RooUnfoldResponse*>(_nCentbins));
    	response_ChPS_counts =  vector<vector<RooUnfoldResponse*> > (_nJetYBins, vector<RooUnfoldResponse*>(_nCentbins));
    	response_FF_fine_counts =  vector<vector<RooUnfoldResponse*> > (_nJetYBins, vector<RooUnfoldResponse*>(_nCentbins));
    	response_ChPS_fine_counts =  vector<vector<RooUnfoldResponse*> > (_nJetYBins, vector<RooUnfoldResponse*>(_nCentbins));
    	
    }
	
	
	for (int j=0;j<_nCentbins;j++){
		for (int i=0;i<_nJetYBins;i++){
		
			//D(z)
			temphist_2D = new TH2D(Form("ff_raw_0_%i_cent%i",i,j),Form("ff_raw_0_%i_cent%i",i,j),zBinsN, zBins, ptJetBinsN, ptJetBins);
			ff_raw.at(i).at(j) = temphist_2D;			
			temphist_2D = new TH2D(Form("ff_raw_fine_0_%i_cent%i",i,j),Form("ff_raw_fine_0_%i_cent%i",i,j),zBinsFineN, zBinsFine, ptJetBinsN, ptJetBins);
			ff_raw_fine.at(i).at(j) = temphist_2D;
			temphist_2D = new TH2D(Form("ff_raw_1_%i_cent%i",i,j),Form("ff_raw_1_%i_cent%i",i,j),zBinsN, zBins, ptJetBinsN, ptJetBins);
			ff_raw_UE.at(i).at(j) = temphist_2D;
			temphist_2D = new TH2D(Form("ff_raw_fine_1_%i_cent%i",i,j),Form("ff_raw_fine_1_%i_cent%i",i,j),zBinsFineN, zBinsFine, ptJetBinsN, ptJetBins);
			ff_raw_UE_fine.at(i).at(j) = temphist_2D;
			
			//D(Pt)
			temphist_2D = new TH2D(Form("ChPS_raw_0_%i_cent%i",i,j),Form("ChPS_raw_0_%i_cent%i",i,j),ptTrkBinsN, ptTrkBins, ptJetBinsN, ptJetBins);
			ChPS_raw.at(i).at(j) = temphist_2D;			
			temphist_2D = new TH2D(Form("ChPS_raw_fine_0_%i_cent%i",i,j),Form("ChPS_raw_fine_0_%i_cent%i",i,j),ptTrkBinsFineN, ptTrkBinsFine, ptJetBinsN, ptJetBins);
			ChPS_raw_fine.at(i).at(j) = temphist_2D;	
			temphist_2D = new TH2D(Form("ChPS_raw_1_%i_cent%i",i,j),Form("ChPS_raw_1_%i_cent%i",i,j),ptTrkBinsN, ptTrkBins, ptJetBinsN, ptJetBins);
			ChPS_raw_UE.at(i).at(j) = temphist_2D;	
			temphist_2D = new TH2D(Form("ChPS_raw_fine_1_%i_cent%i",i,j),Form("ChPS_raw_fine_1_%i_cent%i",i,j),ptTrkBinsFineN, ptTrkBinsFine, ptJetBinsN, ptJetBins);
			ChPS_raw_UE_fine.at(i).at(j) = temphist_2D;
			temphist_3D = new TH3D(Form("ChPS_raw_fine_1_%i_cent%i_response",i,j),Form("ChPS_raw_fine_1_%i_cent%i_response",i,j),ptTrkBinsSumN, ptTrkBinsSum, ptJetBinsN, ptJetBins,respBinsN,respBins);
			ChPS_raw_UE_fine_response.at(i).at(j) = temphist_3D;	
			
			//temphist_2D = new TH2D(Form("ff_truth_matched_%i_cent%i",i,j),Form("ff_truth_matched_%i_cent%i",i,j),zBinsN, zBins, ptJetBinsN, ptJetBins);
			//ff_truth_matched.at(i).at(j) = temphist_2D;
			//temphist_2D = new TH2D(Form("ChPS_matched_truth_%i_cent%i",i,j),Form("ChPS_matched_truth_%i_cent%i",i,j),ptTrkBinsN, ptTrkBins, ptJetBinsN, ptJetBins);
			//ChPS_truth_matched.at(i).at(j) = temphist_2D;
			
			//Jet spectra
			temphist_1D = new TH1D(Form("h_reco_jet_spectrum_%i_cent%i",i,j),Form("h_reco_jet_spectrum_%i_cent%i",i,j),ptJetBinsN, ptJetBins);
			h_reco_jet_spectrum.at(i).at(j) = temphist_1D;
			temphist_1D = new TH1D(Form("h_reco_jet_spectrum_counts_%i_cent%i",i,j),Form("h_reco_jet_spectrum_counts_%i_cent%i",i,j),ptJetBinsN, ptJetBins);
			h_reco_jet_spectrum_counts.at(i).at(j) = temphist_1D;
			temphist_1D = new TH1D(Form("h_reco_jet_spectrum_counts_TMR_%i_cent%i",i,j),Form("h_reco_jet_spectrum_counts_TMR_%i_cent%i",i,j),ptJetBinsN, ptJetBins);
			h_reco_jet_spectrum_counts_TMR.at(i).at(j) = temphist_1D;
			//temphist_1D = new TH1D(Form("h_true_jet_spectrum_matched_%i_cent%i",i,j),Form("h_true_jet_spectrum_matched_%i_cent%i",i,j),ptJetBinsN, ptJetBins);
			//h_true_jet_spectrum_matched.at(i).at(j) = temphist_1D;
			//TODO to be enabled when needed
			//temphist_1D = new TH1D(Form("h_reco_jet_spectrum_weighted_%i_cent%i",i,j),Form("h_reco_jet_spectrum_weighted_%i_cent%i",i,j),ptJetBinsN, ptJetBins);
			//h_reco_jet_spectrum_weighted.at(i).at(j) = temphist_1D;
			//temphist_1D = new TH1D(Form("h_true_jet_spectrum_weighted_%i_cent%i",i,j),Form("h_true_jet_spectrum_weighted_%i_cent%i",i,j),ptJetBinsN, ptJetBins);
			//h_true_jet_spectrum_weighted.at(i).at(j) = temphist_1D;
			
            
            ff_raw.at(i).at(j)->Sumw2();
            ff_raw_fine.at(i).at(j)->Sumw2();
            ff_raw_UE.at(i).at(j)->Sumw2();
            ff_raw_UE_fine.at(i).at(j)->Sumw2();
            ChPS_raw.at(i).at(j)->Sumw2();
            ChPS_raw_fine.at(i).at(j)->Sumw2();
            ChPS_raw_UE.at(i).at(j)->Sumw2();
            ChPS_raw_UE_fine.at(i).at(j)->Sumw2();
            ChPS_raw_UE_fine_response.at(i).at(j)->Sumw2();
            h_reco_jet_spectrum.at(i).at(j)->Sumw2();
            h_reco_jet_spectrum_counts.at(i).at(j)->Sumw2();
            h_reco_jet_spectrum_counts_TMR.at(i).at(j)->Sumw2();
            //h_reco_jet_spectrum_weighted.at(i).at(j)->Sumw2();
            //h_true_jet_spectrum_weighted.at(i).at(j)->Sumw2();
           
            
            wk()->addOutput (ff_raw.at(i).at(j));
            wk()->addOutput (ff_raw_fine.at(i).at(j));
            wk()->addOutput (ff_raw_UE.at(i).at(j));
            wk()->addOutput (ff_raw_UE_fine.at(i).at(j));
            wk()->addOutput (ChPS_raw.at(i).at(j));
            wk()->addOutput (ChPS_raw_fine.at(i).at(j));
            wk()->addOutput (ChPS_raw_UE.at(i).at(j));
            wk()->addOutput (ChPS_raw_UE_fine.at(i).at(j));
            wk()->addOutput (ChPS_raw_UE_fine_response.at(i).at(j));
            
            wk()->addOutput (h_reco_jet_spectrum.at(i).at(j));
            wk()->addOutput (h_reco_jet_spectrum_counts.at(i).at(j));
            wk()->addOutput (h_reco_jet_spectrum_counts_TMR.at(i).at(j));
            //wk()->addOutput (h_reco_jet_spectrum_weighted.at(i).at(j));
            //wk()->addOutput (h_true_jet_spectrum_weighted.at(i).at(j));
            
            
            //MC only
            if (_data_switch==1){
		        temphist_2D = new TH2D(Form("ff_UE_z_%i_cent%i",i,j),Form("ff_UE_z_%i_cent%i",i,j),zBinsN, zBins, ptJetBinsN, ptJetBins);
				ff_UE_z.at(i).at(j) = temphist_2D;
				temphist_2D = new TH2D(Form("ff_UE_pT_%i_cent%i",i,j),Form("ff_UE_pT_%i_cent%i",i,j),ptTrkBinsN, ptTrkBins, ptJetBinsN, ptJetBins);
				ff_UE_pT.at(i).at(j) = temphist_2D;				
				
				temphist_2D = new TH2D(Form("ff_UE_z_fine_%i_cent%i",i,j),Form("ff_UE_z_fine_%i_cent%i",i,j),zBinsFineN, zBinsFine, ptJetBinsN, ptJetBins);
				ff_UE_z_fine.at(i).at(j) = temphist_2D;
				temphist_2D = new TH2D(Form("ff_UE_pT_fine_%i_cent%i",i,j),Form("ff_UE_pT_fine_%i_cent%i",i,j),ptTrkBinsFineN, ptTrkBinsFine, ptJetBinsN, ptJetBins);
				ff_UE_pT_fine.at(i).at(j) = temphist_2D;
				temphist_3D = new TH3D(Form("ff_UE_pT_fine_%i_cent%i_response",i,j),Form("ff_UE_pT_fine_%i_cent%i_response",i,j),ptTrkBinsSumN, ptTrkBinsSum, ptJetBinsN, ptJetBins,respBinsN,respBins);
				ff_UE_pT_fine_response.at(i).at(j) = temphist_3D;
				temphist_3D = new TH3D(Form("UE_pT_Mcorrelation_v_response_%i_cent%i_response_fine",i,j),Form("UE_pT_Mcorrelation_v_response_%i_cent%i_response_fine",i,j),respBinsN,respBins, ptJetBinsN, ptJetBins,respBinsN,respBins);
				UE_pT_Mcorrelation_v_response_fine.at(i).at(j) = temphist_3D;
				
			
				//Truth
				temphist_2D = new TH2D(Form("ff_truth_%i_cent%i",i,j),Form("ff_truth_%i_cent%i",i,j),zBinsN, zBins, ptJetBinsN, ptJetBins);
				ff_truth.at(i).at(j) = temphist_2D;
				temphist_2D = new TH2D(Form("ChPS_truth_%i_cent%i",i,j),Form("ChPS_truth_%i_cent%i",i,j),ptTrkBinsN, ptTrkBins, ptJetBinsN, ptJetBins);
				ChPS_truth.at(i).at(j) = temphist_2D;
				
				
				temphist_2D = new TH2D(Form("ff_truth_fine_%i_cent%i",i,j),Form("ff_truth_fine_%i_cent%i",i,j),zBinsFineN, zBinsFine, ptJetBinsN, ptJetBins);
				ff_truth_fine.at(i).at(j) = temphist_2D;
				//temphist_2D = new TH2D(Form("ff_truth_alt_fine_%i_cent%i",i,j),Form("ff_truth_alt_fine_%i_cent%i",i,j),zBinsFineN, zBinsFine, ptJetBinsN, ptJetBins);
				//ff_truth_alt_fine.at(i).at(j) = temphist_2D;
				temphist_2D = new TH2D(Form("ChPS_truth_fine_%i_cent%i",i,j),Form("ChPS_truth_fine_%i_cent%i",i,j),ptTrkBinsFineN, ptTrkBinsFine, ptJetBinsN, ptJetBins);
				ChPS_truth_fine.at(i).at(j) = temphist_2D;
			
				temphist_1D = new TH1D(Form("h_true_jet_spectrum_%i_cent%i",i,j),Form("h_true_jet_spectrum_%i_cent%i",i,j),ptJetBinsN, ptJetBins);
				h_true_jet_spectrum.at(i).at(j) = temphist_1D;
			    
			    
				//Responses
				temphist_2D = new TH2D(Form("ff_jetResponse_%i_cent%i",i,j),Form("ff_jetResponse_%i_cent%i",i,j),ptJetBinsN, ptJetBins, ptJetBinsN, ptJetBins);
				ff_jetResponse.at(i).at(j) = temphist_2D;			
				temphist_3D = new TH3D(Form("ff_trackpTResponse_%i_cent%i",i,j),Form("ff_trackpTResponse_%i_cent%i",i,j),ptTrkBinsFineN, ptTrkBinsFine,ptTrkBinsFineN, ptTrkBinsFine, ptJetBinsN, ptJetBins);
				ff_trackpTResponse.at(i).at(j) = temphist_3D;
				temphist_3D = new TH3D(Form("ff_trackzResponse_%i_cent%i",i,j),Form("ff_trackzResponse_%i_cent%i",i,j),zBinsFineN, zBinsFine, zBinsFineN, zBinsFine, ptJetBinsN, ptJetBins);
				ff_trackzResponse.at(i).at(j) = temphist_3D;
				temphist_3D = new TH3D(Form("ff_trackzResponse_counts_%i_cent%i",i,j),Form("ff_trackzResponse_counts_%i_cent%i",i,j),zBinsFineN, zBinsFine, zBinsFineN, zBinsFine, ptJetBinsN, ptJetBins);
				ff_trackzResponse_counts.at(i).at(j) = temphist_3D;
			
				response_FF.at(i).at(j) = new RooUnfoldResponse(ff_raw.at(i).at(j),ff_truth.at(i).at(j));
		        response_ChPS.at(i).at(j) = new RooUnfoldResponse(ChPS_raw.at(i).at(j),ChPS_truth.at(i).at(j));
		        response_jet.at(i).at(j) = new RooUnfoldResponse(h_reco_jet_spectrum.at(i).at(j),h_true_jet_spectrum.at(i).at(j));
		        
		        
		        response_FF_fine.at(i).at(j) = new RooUnfoldResponse(ff_raw_fine.at(i).at(j),ff_truth_fine.at(i).at(j));
		        response_ChPS_fine.at(i).at(j) = new RooUnfoldResponse(ChPS_raw_fine.at(i).at(j),ChPS_truth_fine.at(i).at(j));
		        
		        response_FF_counts.at(i).at(j) = new RooUnfoldResponse(ff_raw.at(i).at(j),ff_truth.at(i).at(j),Form("4D_trackzResponse_%i_cent%i_counts",i,j));
		        response_ChPS_counts.at(i).at(j) = new RooUnfoldResponse(ChPS_raw.at(i).at(j),ChPS_truth.at(i).at(j),Form("4D_trackpTResponse_%i_cent%i_counts",i,j));
		        
		        
		        response_FF_fine_counts.at(i).at(j) = new RooUnfoldResponse(ff_raw_fine.at(i).at(j),ff_truth_fine.at(i).at(j),Form("4D_trackzResponse_%i_cent%i_fine_counts",i,j));
		        response_ChPS_fine_counts.at(i).at(j) = new RooUnfoldResponse(ChPS_raw_fine.at(i).at(j),ChPS_truth_fine.at(i).at(j),Form("4D_trackpTResponse_%i_cent%i_fine_counts",i,j));
		    	
		        ff_UE_z.at(i).at(j)->Sumw2();
		        ff_UE_pT.at(i).at(j)->Sumw2();
		        ff_truth.at(i).at(j)->Sumw2();
		        ChPS_truth.at(i).at(j)->Sumw2();
		        
		        ff_UE_z_fine.at(i).at(j)->Sumw2();
		        ff_UE_pT_fine.at(i).at(j)->Sumw2();
		        ff_UE_pT_fine_response.at(i).at(j)->Sumw2();
		        UE_pT_Mcorrelation_v_response_fine.at(i).at(j)->Sumw2();
		        ff_truth_fine.at(i).at(j)->Sumw2();
		        //ff_truth_alt_fine.at(i).at(j)->Sumw2();
		        ChPS_truth_fine.at(i).at(j)->Sumw2();
		        
		        h_true_jet_spectrum.at(i).at(j)->Sumw2();
		        
		        ff_jetResponse.at(i).at(j)->Sumw2();
		        ff_trackpTResponse.at(i).at(j)->Sumw2();
		        ff_trackzResponse.at(i).at(j)->Sumw2();
		        ff_trackzResponse_counts.at(i).at(j)->Sumw2();
            	
		    	    
		        wk()->addOutput (ff_UE_z.at(i).at(j));
		        wk()->addOutput (ff_UE_pT.at(i).at(j));
		        wk()->addOutput (ff_truth.at(i).at(j));
		        wk()->addOutput (ChPS_truth.at(i).at(j));
		        wk()->addOutput (ff_UE_z_fine.at(i).at(j));
		        wk()->addOutput (ff_UE_pT_fine.at(i).at(j));
		        wk()->addOutput (ff_UE_pT_fine_response.at(i).at(j));
		        wk()->addOutput (UE_pT_Mcorrelation_v_response_fine.at(i).at(j));
		        wk()->addOutput (ff_truth_fine.at(i).at(j));
		        //wk()->addOutput (ff_truth_alt_fine.at(i).at(j));
		        wk()->addOutput (ChPS_truth_fine.at(i).at(j));
		        
		        wk()->addOutput (h_true_jet_spectrum.at(i).at(j));
		        
		        wk()->addOutput (ff_jetResponse.at(i).at(j));
		        wk()->addOutput (ff_trackpTResponse.at(i).at(j));
		        wk()->addOutput (ff_trackzResponse.at(i).at(j));
		        wk()->addOutput (ff_trackzResponse_counts.at(i).at(j));
		        wk()->addOutput (response_FF.at(i).at(j));
		        wk()->addOutput (response_ChPS.at(i).at(j));
		        wk()->addOutput (response_FF_fine.at(i).at(j));
		        wk()->addOutput (response_ChPS_fine.at(i).at(j));
		        wk()->addOutput (response_jet.at(i).at(j));
		        
		        //wk()->addOutput (response_FF_counts.at(i).at(j));
		        //wk()->addOutput (response_ChPS_counts.at(i).at(j));
		        //wk()->addOutput (response_FF_fine_counts.at(i).at(j));
		        //wk()->addOutput (response_ChPS_fine_counts.at(i).at(j));
            }	
            			
		}
		temphist_1D = new TH1D(Form("h_reco_jet_spectrum_fine_cent%i",j),Form("h_reco_jet_spectrum_fine_cent%i",j),50,40,140);
		h_reco_jet_spectrum_fine.push_back(temphist_1D);
		wk()->addOutput (h_reco_jet_spectrum_fine.at(j));
		
		if (_data_switch==1){
			temphist_1D = new TH1D(Form("h_true_jet_spectrum_fine_cent%i",j),Form("h_true_jet_spectrum_fine_cent%i",j),50,40,140);
			h_true_jet_spectrum_fine.push_back(temphist_1D);
			wk()->addOutput (h_true_jet_spectrum_fine.at(j));
		}	
	}
	
	cout << " Histograms ready" << endl;


	return EL::StatusCode::SUCCESS;
}
