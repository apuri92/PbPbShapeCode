
#ifndef pPbFragmentation_h
#define pPbFragmentation_h

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

#include <iostream>

// Fixed size dimensions of array or collections stored in the TTree if any.

using namespace std;

class pPbFragmentation : public BaseClass{
	
	public :
	 
	float _b_track_pt[100];
	float _b_track_eta[100];
	float _b_track_phi[100];
	
	
	//Internal configuration
	
	float _dR_max;
	
	UEEstimator* uee; //!
				
	//Basic histograms
	TH3D *hET_ETsub; //!
	TH2D *htrig_reco_pt; //!
	
	vector<TH1D*> hjetpt_trig;
	vector<TH2D*> hjetz_trig;
	vector<TH2D*> hjetz_UE_trig;
	
	TH1D *hz; //!
	TH1D *hz_UE; //!
	TH2D *truth_hjetz; //!

	//trigger
	/*
	string _trigger_collection;
	string trigger_chains[10];
	float trigger_thresholds[10];
	float jet_pt_trig[10][2];
	int _nTriggers;
	*/
	
	bool event_isTriggered[10]; //!
	bool trigger[10]; //!
	int jet_isTriggered[10]; //!
	float trig_prescale[10]; //!
	
	//Output tree
	TTree *tree_ff; //!
	TTree *truth_tree_ff; //!
		
	float FCalEt; //!
	int multiplicity; //!
	
	float truth_jet_pt, truth_jet_phi, truth_jet_eta, truth_jet_m, truth_reco_jet_dR,truth_jet_Muon_pT,truth_jet_electron_pT; //!
	int truth_jet_isDummy,truth_jet_flavour; //!
	
	float jet_uJER;//!
	vector<float> jet_uJES; //!
	
	int jet_isGood,jet_isDummy,jet_hasTruth,jet_centrality,jet_flavour; //!
	float jet_NBJ_pT,truth_jet_NBJ_pT, jet_Muon_pT, jet_electron_pT; //!
	vector<int> Is_jet_Isolated; //!
	
	vector<float> track_pt, track_eta, track_phi, track_z,InJet_muon_pT,InJet_electron_pT; //!
	vector<int> track_charge; //!
	vector<float> track_UE_charge,track_UE_z,track_UE_weight; //!
	vector<float> track_d0,track_z0sintheta,track_ed0,track_ez0sintheta; //!
	vector<int> track_nPixHits,track_nPixHoles,track_nShPixH, track_nSCTHits,track_nSCTHoles,track_nShSCTH, track_nTRTHits; //!
	//vector<int> track_nPixDeadS,track_nSCTDeadS; //!
	vector<int> track_nIBLHits,track_expIBLHits, track_nBLHits,track_expBLHits, track_vertex_type, track_multiJetMatch, track_multiJetMatch_B; //!
	vector<int> track_passedTrkSelTool[10]; //!
	
	vector<float> track_mc_phi, track_mc_eta, track_mc_pt, track_mc_probability, track_isMuon; //!
	vector<int> track_mc_pdg, track_mc_barcode, track_mc_type, track_mc_charge; //!
	
	vector<float> truth_track_pt, truth_track_eta, truth_track_phi, truth_track_z; //!
	vector<int> truth_track_pdg, truth_track_barcode, truth_track_charge, truth_track_multiJetMatch, truth_track_multiJetMatch_B; //!
	
	// this is a standard constructor
	pPbFragmentation ();
 
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
	ClassDef(pPbFragmentation, 1);
	
};

#endif
