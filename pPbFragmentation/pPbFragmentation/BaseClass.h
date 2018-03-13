#ifndef BaseClass_h
#define BaseClass_h

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
#include <xAODEgamma/ElectronContainer.h>
#include "xAODCore/ShallowAuxContainer.h"
#include "xAODCore/ShallowCopy.h"

#include "MuonSelectorTools/MuonSelectionTool.h"
#include "MuonMomentumCorrections/MuonCalibrationAndSmearingTool.h"

#include "ElectronPhotonFourMomentumCorrection/EgammaCalibrationAndSmearingTool.h"
#include "ElectronPhotonSelectorTools/AsgElectronLikelihoodTool.h"
#include "ElectronPhotonSelectorTools/AsgForwardElectronIsEMSelector.h"

#include "InDetTrackSelectionTool/InDetTrackSelectionTool.h"
//#include "InDetTrackSystematicsTools/InDetTrackBiasingTool.h"

#define private public
#include "xAODHIEvent/HIEventShapeAuxContainer.h"
#undef private
#include "xAODHIEvent/HIEventShapeContainer.h"

#include "pPbFragmentation/FragmentationHelperTools.h"
#include "pPbFragmentation/JetHelperTools.h"
#include "pPbFragmentation/UEEstimator.h"
#include "pPbFragmentation/GlobalHelper.h"
#include "pPbFragmentation/TrackHelperTools.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthVertexContainer.h"
#include "xAODMuon/MuonContainer.h"
#include "JetUncertainties/JetUncertaintiesTool.h"
#include "JetResolution/JERTool.h"
#include "JetResolution/JERSmearingTool.h"
#include <iostream>

using namespace std;

class BaseClass : public EL::Algorithm
{

	public:

	//Configuration from main macro
	int _data_switch; //
	int _dataset; //
	int _isMB; //
	int _isHerwig; //
	int _jet_radius; //
	float _trkptBkgrThreshold; //
	string _reco_jet_collection; //
	string _test_reco_jet_collection; //
	string _truth_jet_collection; //
	string _GRL; //
	string _cut_level; //
	float _dR_max; //
	int _centrality_scheme; //
	int _Run_Number; //
	int _truth_only; //
	float _dR_truth_matching; //
	int _trk_selection; //
	float _pTtrkCut; //
	float _pTjetCut; //
	float _truthpTjetCut; //
	float _mcProbCut; //
	int _doClusters; //
	int _doSlimTree;
	bool _doForward; //
	std::string _outputName; //
	double _pt_iso; //
	double _JERBalancecut; //
	//trigger
	string _trigger_collection; //
	vector <string> trigger_chains; //TODO
	vector <float> trigger_PS; //TODO
	vector <float> trigger_thresholds; //TODO
	vector <vector<float>> jet_pt_trig; //TODO
	int _nTriggers; //
	bool _applyReweighting; //
	float _jetptBkgrThreshold; //
	int _uncert_index; //
	bool _doPileupRejection;
	int _useCharge;
	bool _correctTrackpT;
	bool _eff_jety;
	bool _UseAltzDef;
	bool _doFJR;
	float _maxjetdeltaR;
	bool _doJPRCorrection;
	float _jet_y_cut;
	bool _doCoarsTrackpT;
	int _PythiaPowheg;

	bool event_isTriggered[10]; //!	//TODO
	bool trigger[10]; //!	//TODO
	int jet_isTriggered[10]; //! //TODO
	float trig_prescale[10]; //! //TODO
	
	//Evnets 
	int m_eventCounter; //!
	int event_n, run_n, lbn_n; //!
	
	//Jets
	float jet_pt, jet_Et, jet_eta, jet_y, jet_phi,jet_m, jet_a2_pt, jet_a2_eta, jet_a2_phi, jet_pt_seb, jet_pt_prexcalib, jet_pt_xcalib, jet_pt_unsubtracted, jet_pt_EM; //!
	int jet_nConst; //!
	//MC	
	//int _barcodeMin; //!
	//int _barcodeMax; //!
	
	//Histograms

	TH1D *h_RejectionHisto; //!
	TH2D *h_triggercounter; //!
	TH1D *h_FCal_Et; //!
	TH1D *h_FCal_Et_unw; //!

	TH3D* h_reco_trk_map; //!
	TH3D* h_reco_trk_map_nocuts; //!
	
	vector<TH3D*> h_PixHits; //!
	vector<TH3D*> h_SCTHits; //!
	vector<TH2D*> h_d0; //!
	vector<TH2D*> h_z0sintheta; //!
	vector<TH3D*> h_vx; //!
	vector<TH3D*> h_BL; //!
	
	vector<TH3D*> h_cut_flow; //!
	vector<TH3D*> h_cut_flow_PV; //!
	vector<TH3D*> h_cut_flow_PV_jet; //!
	vector<TH3D*> h_cut_flow_jet; //!
	
	TF1* f_d0_cut; //!
	
	//Vertex parameters
	int vx_n; //!
	vector<float> vx_sumPt; //!
	vector<int> vx_nTracks; //!
	
	//Timing
	float mbtime_timeA; //!
	float mbtime_timeC; //!
	
		
	//Jet cleaning
	JetCleaningTool *m_jetCleaning; //!
	//JER tool
	JERTool *jerTool; //!
	JERSmearingTool *smearTool; //!

	//Uncertainties tool
	JetUncertaintiesTool *jesProv; //!	
		
	//GRL
	GoodRunsListSelectionTool *m_grl; //!
	
	//Muon selection and correction
	CP::MuonSelectionTool *m_muonSelection; //!
	CP::MuonCalibrationAndSmearingTool *m_muonCalibrationAndSmearingTool; //!
	
	//Electron selection and correction
	CP::EgammaCalibrationAndSmearingTool* m_EgammaCalibrationAndSmearingTool;  //!
	AsgElectronLikelihoodTool*      m_LHToolTight2015;  //!
	
	//Trigger tools member variables
	Trig::TrigDecisionTool *m_trigDecisionTool; //!
	TrigConf::xAODConfigTool *m_trigConfigTool; //!
	JetCalibrationTool * m_jetCalibration; //!
	JetCalibrationTool * m_jetCalibration_insitu; //!
	vector<const Trig::ChainGroup*> _chainGroup; //!
	vector<const Trig::ChainGroup*> _referenceChainGroup; //!
	int _first_trigger;
	
	//Prescale sets
	TFile * f_trigger_RunNumber_prescale; //!
	TH2F * h2_trigger_RunNumber_prescale; //!

	//Track Selection Tool
	int _nTrkSelTools;
	InDet::InDetTrackSelectionTool * m_trackSelectionTool[10]; //!
	InDet::InDetTrackSelectionTool * m_trackSelectorTool; //!
	InDet::InDetTrackSelectionTool * m_trackSelectorTool_FJR; //!

	//Corrections for track pT sagitta bias
	//InDet::InDetTrackBiasingTool * m_trkBiasingTool; //!

	//Bool_t IsGoodMuon(int index);
	void SetTrigger_chains();
	void SetTrigger_hist(TH2D* h);
	
	//virtual float dR(float, float, float, float);
	virtual bool HasVertex() const;
	virtual bool HasGoodTiming() const;
	virtual bool HasPileUp() const;
	virtual std::vector<float> FilterJetsByMuons(std::vector<float>& jetEta, std::vector<float>& jetPhi, std::pair< xAOD::MuonContainer*, xAOD::ShallowAuxContainer* >& muon_container, Float_t deltaRMin); 
    virtual std::vector<float> FilterJetsByElectrons(std::vector<float>& jetEta, std::vector<float>& jetPhi, std::pair< xAOD::ElectronContainer*, xAOD::ShallowAuxContainer* >& electron_container, Float_t deltaRMin);
	// this is a standard constructor
	BaseClass ();
	BaseClass(const BaseClass& base);
	
	// this is needed to distribute the algorithm to the workers
	ClassDef(BaseClass, 1);
};

#endif
