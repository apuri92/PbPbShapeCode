
#ifndef JetPerformance_h
#define JetPerformance_h

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
#include "pPbFragmentation/JetCorrector.h"

#define private public
#include "xAODHIEvent/HIEventShapeAuxContainer.h"
#undef private
#include "xAODHIEvent/HIEventShapeContainer.h"

#include "pPbFragmentation/JetHelperTools.h"
#include "pPbFragmentation/GlobalHelper.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthVertexContainer.h"
#include <iostream>

// Fixed size dimensions of array or collections stored in the TTree if any.

using namespace std;

class JetPerformance : public BaseClass{
	
	public :
	 

	
	//Internal configuration
	
	float _dR_max;
	
	//Output tree
	TTree *tree_performance; //!
	TTree *truth_tree_performance; //!
	
				
	//Basic histograms
	TH3D *hET_ETsub; //!
	TH3D *hET_ETsub_v_dEta; //!
	TH2D *htrig_reco_pt; //!
	TH3D *h3_Jet_EtaPt_Nconst; //!
	TH3D *h3_Jet_EtaPt_EtsubPerNconst; //!
	
	TH2D * h2_JetCl_DiffEtRaw; //!
    TH2D * h2_JetCl_DiffEtAlt; //!
    TH2D * h2_JetCl_DiffEtCal; //!
    TH2D * h2_JetCl_PtNconst1; //!
    TH2D * h2_JetCl_PtNconst2; //!
	TH3D * h3_JetCl_NegEt; //!
	TH3D * h3_HLT_jet_spect; //!

	TH1D * h_DAQErrors; //!
	
	vector<TH1D*> hjetpt_trig; //!
	vector<TH2D*> hjetz_trig; //!
	vector<TH2D*> hjetz_UE_trig; //!
	
	vector<TH2D*> h_jet_v_mass; //!
	vector<TH2D*> h_truth_jet_v_mass; //!
	
	bool event_isTriggered[10]; //!
	bool trigger[10]; //!
	int jet_isTriggered[10]; //!
	float trig_prescale[10]; //!
	
	float jet_ClSub_et, jet_ClUnsub_et, jet_Clneg_et, jet_Clpost_et; //!
			
	float FCalEt; //!
	
	float jet_NBJ_pT;
	float truth_jet_pt, truth_jet_phi, truth_jet_eta, truth_jet_m, truth_reco_jet_dR, truth_jet_NBJ_pT; //!
	int truth_jet_isMuonIsolated,truth_jet_isDummy,truth_jet_centrality,truth_jet_flavour; //!
	
	int jet_isGood,test_jet_isGood,jet_isMuonIsolated,jet_isDummy,jet_hasTruth,jet_centrality,jet_flavour,jet_jetLArBadHVNCell; //!
	float test_jet_eta, test_jet_phi,test_jet_m, test_jet_pt_EM; //!
	float jet_jetQuality,jet_jetTime,jet_jethecq,jet_jetnegE,jet_jetemf,jet_jethecf,jet_jetfracSamplingMax,jet_jetchf,jet_jetBchCorrCell,jet_jetLArBadHVEnergyFrac,weight; //!
	
	vector<float> ClUnsub_et, ClUnsub_eta, ClUnsub_phi; //!
	
	JetCorrector* jetcorr; //!
	
	// this is a standard constructor
	JetPerformance ();
 
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
	ClassDef(JetPerformance, 1);
	
};

#endif
