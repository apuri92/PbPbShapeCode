
#ifndef TrackingPerformance_h
#define TrackingPerformance_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TMath.h>

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
#include <iostream>

// Fixed size dimensions of array or collections stored in the TTree if any.

using namespace std;
class TrackingPerformance : public BaseClass{
	public :



	//Internal configuration

	int NCuts; //!
	vector<string> CutsOn; //!
	int nCentBins;


	vector<TH3D*> h_reco_Injet_matched; //!
	vector<TH3D*> h_eff_Injet_matched; //!
	vector<vector<TH3D*>> h_eff_dR_matched; //!
	vector<vector<TH3D*>> h_eff_dR; //!
	vector<vector<TH3D*>> h_eff_deta_matched; //!
	vector<vector<TH3D*>> h_eff_deta; //!
	vector<vector<TH3D*>> h_eff_dphi_matched; //!
	vector<vector<TH3D*>> h_eff_dphi; //!

	vector<TH3D*> h_eff_matched; //!
	vector<TH3D*> h_eff_total; //!

	vector<TH3D*> h_eff_Injet; //!
	vector<TH3D*> h_eff_entries; //!
	vector<TH3D*> h_eff_Injet_entries; //!
	vector<TH3D*> h_trk_foreff_matched; //!
	vector<TH3D*> h_trk_foreff_full; //!
    vector<TH3D*> h_trk_foreff_entries; //!
    vector<TH3D*> h_fake_v_jet; //!
	vector<TH3D*> h_fake_v_jet_PV; //!
	vector<TH3D*> h_trk_scale; //!
	vector<TH3D*> h_trk_eta_scale; //!
	vector<TH3D*> h_trk_phi_scale; //!
	vector<TH3D*> h_trk_R_scale; //!
	vector<TH3D*> h_truth_outside; //!
	vector<TH3D*> h_reco_outside; //!
	vector<vector<TH3D*> > h_cut_flow_cent; //!
	vector<TH3D*> h_tmp; //!
	vector<TH3D*> h_trk_eff_matched_map; //!
	vector<TH3D*> h_trk_eff_map; //!
	vector<TH3D*> h_reco_jet_map; //!

	TH3D* h_trk_resolution; //!
	TH3D* h_truth_trk_map; //!

	TrackCorrector* trkcorr; //!
	JetCorrector* jetcorr; //!


//	TH1 *h_z0pvs_before; //!
//	TH1 *h_z0pvs_after; //!
//	TH1 *h_z0pvs_applied; //!
//	TH1 *h_z0pvs_tmp; //!

	TH3 *h_d0sign_wrt_truth; //!
	TH3 *h_z0sign_wrt_truth; //!
	TH3 *h_d0sign_subset; //!
	TH3 *h_z0sign_subset; //!

	TH2 *h_d0_cut; //!
	TH2 *h_d0_nocut; //!
	TH1* h_dR_axis; //!
	// this is a standard constructor
	TrackingPerformance ();

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
	ClassDef(TrackingPerformance, 1);

};

#endif
