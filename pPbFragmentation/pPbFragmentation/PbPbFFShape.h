
#ifndef PbPbFFShape_h
#define PbPbFFShape_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TMath.h>
#include <TRandom3.h>
//#include <THnSparse.h>

#include "BaseClass.h"

#include <EventLoop/Algorithm.h>
#include "xAODEventInfo/EventInfo.h"
#include "TrigDecisionTool/TrigDecisionTool.h"
#include "TrigConfxAOD/xAODConfigTool.h"
#include "GoodRunsLists/GoodRunsListSelectionTool.h"
#include "JetCalibTools/JetCalibrationTool.h"
#include <boost/regex.hpp>
#include "JetSelectorTools/JetCleaningTool.h"
#include "xAODTracking/TrackParticleContainer.h"
#include "xAODTracking/VertexContainer.h"

#define private public
#include "xAODHIEvent/HIEventShapeAuxContainer.h"
#undef private
#include "xAODHIEvent/HIEventShapeContainer.h"

#include "pPbFragmentation/JetHelperTools.h"
#include "pPbFragmentation/FragmentationHelperTools.h"
#include "pPbFragmentation/UEEstimator.h"
#include "pPbFragmentation/GlobalHelper.h"
#include "pPbFragmentation/TrackHelperTools.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthVertexContainer.h"
#include "pPbFragmentation/TrackCorrector.h"
#include "pPbFragmentation/JetCorrector.h"
#include "pPbFragmentation/UncertProvider.h"

#include "Asg_RooUnfold/RooUnfoldResponse.h"

#include "HIEventUtils/HIPileupTool.h"
#include "ZdcAnalysis/ZdcAnalysisTool.h"

#include <iostream>

// Fixed size dimensions of array or collections stored in the TTree if any.

using namespace std;

class PbPbFFShape : public BaseClass{
	
	public :


	float _b_track_pt[100];
	float _b_track_eta[100];
	float _b_track_phi[100];
	
	
	//Internal configuration	
	UEEstimator* uee; //!
	
	UncertProvider* uncertprovider; //!
	
	TrackCorrector* trkcorr; //!
	JetCorrector* jetcorr; //!

	//Axis histograms
	TH3D *h_jet_pt_eta_phi; //!
	TH3D *h_trk_pt_eta_phi; //!


	//Basic histograms
	TH3D *hET_ETsub; //!
	TH1D *h_FCal_Et_restr; //!
	TH2D *h_fcal_change; //!
	TH1D *h_fcal_diff; //!


	TH3D* h_tmp_trk; //!
	TH1D* h_tmp_dR; //!
	TH1D* h_tmp_coneIndex; //!
	TH1D* h_tmp_dRBin; //!
	TH3D* h_tmp_rdEtadPhi; //!
	TH2D* h_cone_map; //!
	TH1D* h_tmp_cone_stats; //!


	TH1 *h_event_rN; //!
	TH1 *h_fcal_mc; //!
	TH1 *h_fcal_mbov; //!
	TH1 *h_cent_mc; //!
	TH1 *h_cent_mbov; //!

	TH1D * h_centrality; //!
	TH1D *h_dR_binning; //!

	HI::HIPileupTool *m_hiPileup;    //!
	ZDC::ZdcAnalysisTool *m_zdcTools; //!

	vector<vector<TH3D*>> h_dR_change; //!

	bool derive_UE_mode = 0;
//	int lo_jetpt_bin = 6;
//	int hi_jetpt_bin = 8;
	int lo_jetpt_bin = 9;
	int hi_jetpt_bin = 11;

	//debugging histograms
	vector<TH3D*> h_reco_truth_matched; //!
	vector<TH3D*> h_jet_for_eff; //!
	vector<TH3D*> h_jet_for_eff_full; //!
	vector<TH3D*> h_jet_psi3; //!
	vector<TH3D*> h_reco_jets; //!


	//Raw disitrbutions
	vector<vector<TH2D*>> ChPS_raw; //!

	//UE
	vector<vector<vector<vector<TH3D*>>>> h_UE_dNdEtadPhidpT; //!
	vector<vector<TH3D*>> h_jet_v_Psi; //!
//	THnSparse* h_UE_dNdEtadPhidpT; //!
//	double fill_variable[7]; //!


	//test UE
	vector<TH1D*> cone_norm_jet;//!
	vector<TH1D*> MB_norm_jet;//!
	vector<TH2D*> MB_norm_jet_rN;//!
	vector<TH1D*> TM_norm_jet;//!
	vector<TH2D*> TM_norm_jet_rN;//!
	vector<TH1D*> FS_norm_jet;//!

	//Raw disitrbutions from random cones
	vector<vector<TH2D*>> ChPS_cone_UE; //!
	vector<vector<TH2D*>> ChPS_MB_UE; //!
	vector<vector<TH3D*>> ChPS_MB_UE_rN; //!
	vector<vector<TH2D*>> ChPS_MB_UE_err; //!
	vector<vector<TH2D*>> ChPS_FS_UE; //!
	vector<vector<TH2D*>> ChPS_FNS_UE; //!
	vector<vector<TH2D*>> ChPS_circ_UE; //!

	//Truth distributions
	vector<vector<TH2D*>> ChPS_truth; //!
	vector<vector<TH2D*>> ChPS_truth_deta; //!
	vector<vector<TH2D*>> ChPS_truth_dphi; //!

	//Truth matched distributions
	//vector<vector<TH2D*>> ChPS_truth_matched; //!
	
	//UE unmatched disitrbution
	vector<vector<TH2D*>> ChPS_TM_UE; //!
	vector<vector<TH3D*>> ChPS_TM_UE_rN; //!

	vector<vector<TH1D*>> h_reco_jet_spectrum; //!
    vector<vector<TH1D*>> h_true_jet_spectrum; //!
	vector<vector<TH1D*>> h_true_jet_spectrum_UE_norm; //!


   	//Responses
   	//4D
    vector<vector<RooUnfoldResponse *>> response_ChPS; //!

	//2D
    vector<vector<RooUnfoldResponse *>> response_jet; //!


	//2D for quick check
	vector<vector<TH2D*>> ff_jetResponse; //!
	vector<vector<TH3D*>> ff_trackpTResponse; //!

//	vector<TH3D*> h_jetpT_v_multiplicity; //!

	TF1 *f_d0_cut; //!
	
	bool event_isTriggered[10]; //!
	bool trigger[10]; //!
	int jet_isTriggered[10]; //!
	float trig_prescale[10]; //!
	
	float FCalEt;

	float truth_jet_pt, truth_jet_phi, truth_jet_eta, truth_jet_y, truth_jet_m, truth_reco_jet_dR,truth_jet_Muon_pT,truth_jet_electron_pT; //!
	int truth_jet_isDummy,truth_jet_centrality,truth_jet_flavour; //!


	//float jet_uJER;
	//vector<float> jet_uJES; //!
	
	int jet_isGood,jet_isDummy,jet_hasTruth,jet_centrality,jet_flavour; //!
	float jet_NBJ_pT,truth_jet_NBJ_pT, jet_Muon_pT, jet_electron_pT; //!
	vector<int> Is_jet_Isolated; //!

	vector<float> track_pt, track_eta, track_phi, track_z,InJet_muon_pT,InJet_electron_pT; //!
	vector<int> track_charge; //!
	vector<float> track_UE_charge,track_UE_z,track_UE_weight; //!
	vector<float> track_d0,track_z0sintheta,track_ed0,track_ez0sintheta; //!
	vector<int> track_nPixHits,track_nPixHoles,track_nPixDeadS,track_nShPixH, track_nSCTHits,track_nSCTHoles,track_nSCTDeadS,track_nShSCTH, track_nTRTHits; //!
	vector<int> track_nIBLHits,track_expIBLHits, track_nBLHits,track_expBLHits, track_vertex_type, track_multiJetMatch, track_multiJetMatch_B; //!
	vector<int> track_passedTrkSelTool; //!
	
	vector<float> track_mc_phi, track_mc_eta, track_mc_pt, track_mc_probability, track_isMuon; //!
	vector<int> track_mc_pdg, track_mc_barcode, track_mc_type, track_mc_charge; //!
	
	vector<float> truth_track_pt, truth_track_eta, truth_track_phi, truth_track_z; //!
	vector<int> truth_track_pdg, truth_track_barcode, truth_track_charge, truth_track_multiJetMatch, truth_track_multiJetMatch_B; //!
	
	// this is a standard constructor
	PbPbFFShape ();
 
	// these are the functions inherited from Algorithm
	virtual EL::StatusCode setupJob (EL::Job& job);
	virtual EL::StatusCode fileExecute ();
	virtual EL::StatusCode histInitialize ();
	virtual EL::StatusCode changeInput (bool firstFile);
	virtual EL::StatusCode initialize ();
	virtual EL::StatusCode execute ();
	virtual EL::StatusCode postExecute ();
	virtual EL::StatusCode finalize ();
	virtual EL::StatusCode histFinalize ();	
		
	// this is needed to distribute the algorithm to the workers
	ClassDef(PbPbFFShape, 1);
	
};

#endif
