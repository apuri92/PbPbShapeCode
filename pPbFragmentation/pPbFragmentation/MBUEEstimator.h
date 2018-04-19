
#ifndef MBUEEstimator_h
#define MBUEEstimator_h

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

#include "BaseClass.h"

#include <EventLoop/Algorithm.h>
#include "xAODEventInfo/EventInfo.h"
#include "TrigDecisionTool/TrigDecisionTool.h"
#include "TrigConfxAOD/xAODConfigTool.h"
#include "GoodRunsLists/GoodRunsListSelectionTool.h"
#include <boost/regex.hpp>
#include "xAODTracking/TrackParticleContainer.h"
#include "xAODTracking/VertexContainer.h"

#define private public
#include "xAODHIEvent/HIEventShapeAuxContainer.h"
#undef private
#include "xAODHIEvent/HIEventShapeContainer.h"

#include "pPbFragmentation/FragmentationHelperTools.h"
#include "pPbFragmentation/GlobalHelper.h"
#include "pPbFragmentation/TrackHelperTools.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthVertexContainer.h"
#include "pPbFragmentation/TrackCorrector.h"
#include "pPbFragmentation/UncertProvider.h"
#include "pPbFragmentation/JetCorrector.h"

#include "HIEventUtils/HIPileupTool.h"
#include "ZdcAnalysis/ZdcAnalysisTool.h"

#include <iostream>

// Fixed size dimensions of array or collections stored in the TTree if any.

using namespace std;

class MBUEEstimator : public BaseClass{
	
	public :

	
	//Internal configuration

	float _dR_max;
		
	UncertProvider* uncertprovider; //!
	
	TrackCorrector* trkcorr; //!
	
	JetCorrector* jetcorr; //!
				
	bool m_shapeFF;
	
	
	
	//vector<vector<TH3D*>> h_trk_dNdEtadPhidpT; //!
	vector<vector<vector<vector<TH3D*>>>> h_UE_dNdEtadPhidpT_HP; //!
	vector<vector<vector<vector<TH3D*>>>> h_UE_dNdEtadPhidpT_HP_fine; //!
	vector<vector<TH3D*>> h_jet_v_Psi_HP; //!
	//vector<vector<vector<TH3D*>>> h_UE_dNdEtadPhidpT_MC; //!
	//vector<vector<vector<TH3D*>>> h_UE_dNdEtadPhidpT_MC_fine; //!
	//vector<TH3D*> h_jet_v_Psi_MC; //!
	
	TH2D * map_excluded_jets;//!
	
	TH1D * h_centrality_HP; //!
	TH1D * h_FCal_Et_HP; //!
	TH1D * h_centrality_MC; //!
	TH1D * h_FCal_Et_MC; //!
	TH1D * h_dPsi; //!
	TH1D * h_Psi; //!
	TH1D * h_jet_phi; //!
	TH2D * h_jet_phi_v_dPsi; //!
	TH2D * h_jet_phi_v_Psi; //!
	
	HI::HIPileupTool                      *m_hiPileup;    //!
	ZDC::ZdcAnalysisTool *m_zdcTools; //!
		
	TF1 *f_d0_cut; //!
	
	bool event_isTriggered[10]; //!
	bool trigger[10]; //!
	int jet_isTriggered[10]; //!
	float trig_prescale[10]; //!
	
	float FCalEt;
	
	bool FCal_only;

	
	
	vector<float> track_pt, track_eta, track_phi, track_z,InJet_muon_pT,InJet_electron_pT; //!
	vector<int> track_charge; //!
	vector<float> track_d0,track_z0sintheta,track_ed0,track_ez0sintheta; //!
	vector<int> track_nPixHits,track_nPixHoles,track_nPixDeadS,track_nShPixH, track_nSCTHits,track_nSCTHoles,track_nSCTDeadS,track_nShSCTH, track_nTRTHits; //!
	vector<int> track_nIBLHits,track_expIBLHits, track_nBLHits,track_expBLHits, track_vertex_type, track_multiJetMatch, track_multiJetMatch_B; //!
	vector<int> track_passedTrkSelTool; //!
	
	vector<float> track_mc_phi, track_mc_eta, track_mc_pt, track_mc_probability, track_isMuon; //!
	vector<int> track_mc_pdg, track_mc_barcode, track_mc_type, track_mc_charge; //!
	
	vector<float> truth_track_pt, truth_track_eta, truth_track_phi, truth_track_z; //!
	vector<int> truth_track_pdg, truth_track_barcode, truth_track_charge, truth_track_multiJetMatch, truth_track_multiJetMatch_B; //!
	
	// this is a standard constructor
	MBUEEstimator ();
 
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
	ClassDef(MBUEEstimator, 1);
	
};

#endif
