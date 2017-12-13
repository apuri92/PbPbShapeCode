
#ifndef InclusiveJetsEventLoop_h
#define InclusiveJetsEventLoop_h

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
#include "pPbFragmentation/SEBCorrectorTool.h"
#include "pPbCentrality/pPbMinBiasUtil.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthVertexContainer.h"
#include <iostream>

// Fixed size dimensions of array or collections stored in the TTree if any.

using namespace std;

class InclusiveJetsEventLoop : public BaseClass{
  
public :
  
  // basic diagnostic plots for global event
  TH1D* h_FCal_Et_withJet;	//!
  
  // basic diagnostic plots for clusters
  TH1F* h1_FCalEtClPt_1;	//!
  TH1F* h1_FCalEtCl_1;	//!
  TH1F* h1_FCalEtClPt_2;	//!
  TH1F* h1_FCalEtCl_2;	//!
  TH2F* h2_cl_EtaPhi1[7]; 	//!
  TH2F* h2_cl_EtaPhi2[7]; 	//!
  TH2F* h2_cl_EtaPhiPt1[7]; 	//!
  TH2F* h2_cl_EtaPhiPt2[7]; 	//!
  TH3F* h3_cl_EtaPhiFCalEtPt_1;	//!
  TH3F* h3_cl_EtaPhiFCalEt_1;	//!			
  TH3F* h3_cl_EtaPhiFCalEtPt_2;	//!			
  TH3F* h3_cl_EtaPhiFCalEt_2;	//!			
  
  // basic diagnostic plots for all three jet parameters
  TH2F* h2_jet42_Pt;	//!
  TH2F* h2_jet32_Pt;	//!
  TH1F* h1_jet4_Pt;	//!
  TH1F* h1_jet3_Pt;	//!
  TH1F* h1_jet2_Pt;	//!
  //<RS>
  
  TH1F* h1_stat[7];  //!
  TH1F* h1_NJet[7]; //!
 
  TH1F* h1_pt_spect[7][10];   //!
  TH1F* h1_pt_spect_trig[7][10];  //!
  
  TH1F* h1_pt_spect_truth[7];   //!
  TH1F* h1_pt_spect_truth2[7];   //!
  TH1F* h1_pt_spect_truth_match[7];  //!
  TH1F* h1_pt_spect_truth_match2[7];  //!
  
  TH2F* h2_response[7]; //!
  TH2F* h2_response2[7]; //!
  TH2F* h2_response3[7]; //!
  
  TH2F* h2_jet4_JES[7];   //!
  TH2F* h2_jet3_JES[7];   //!
  TH2F* h2_jet2_JES[7];   //!
  
  TH3F* h3_jet4_JES_eta[7]; //!
  TH3F* h3_jet4_JES_phi[7]; //!
  TH2F* h2_jet4_JES_inPlane[7]; //!
  TH2F* h2_jet4_JES_outPlane[7]; //!
  TH2F* h2_jet4_JES_dpT_inPlane[7]; //!
  TH2F* h2_jet4_JES_dpT_outPlane[7]; //!
  TH2F* h2_jet4_JES_eta1[7]; //!
  TH2F* h2_jet4_JES_eta2[7]; //!
  
  TH2F* h2_jet4_dR[7];   //!
  TH2F* h2_jet4_deta[7];   //!
  TH2F* h2_jet4_dphi[7];   //!
  
  TH2F* h2_jet4_JES_truth[7];   //!
  TH3F* h3_jet4_JES_truth_eta[7]; //!
  TH3F* h3_jet4_JES_truth_phi[7]; //!
  TH2F* h2_jet4_JES_truth_inPlane[7]; //!
  TH2F* h2_jet4_JES_truth_outPlane[7]; //!
  TH2F* h2_jet4_JES_truth_dpT_inPlane[7]; //!
  TH2F* h2_jet4_JES_truth_dpT_outPlane[7]; //!
  TH2F* h2_jet4_JES_truth_eta1[7]; //!
  TH2F* h2_jet4_JES_truth_eta2[7]; //!
  
  TH2F* h2_jet4_truth_dR[7];   //!
  TH2F* h2_jet4_truth_deta[7];   //!
  TH2F* h2_jet4_truth_dphi[7];   //!
  
  TH2F* h2_jet4_truth_etaPhi[7]; //!
  TH2F* h2_jet4_truth_etaPhi_weightedJES[7]; //!
  
  Float_t vN_fcal; //!
  
  TH2F* h2_jet4_etaPhi[7]; //!
  TH2F* h2_jet4_etaPhi_weightedJES[7]; //!
  
  TH2F* h2_jet4FJR_PtNConst_binning2[7]; //!
  
  TH2F* h2_jet4FJR_FCalEt_subPt;  //!
  
  //</RS>
  
  TH2F* h2_jet4_PtEta[7];   //!
  TH2F* h2_jet4FJR_PtEta[7];   //!
  TH2F* h2_jet4FJR_PtSubtrPt[7];   //!
  TH2F* h2_jet4FJR_PtNConst[7];   //!
  
  TH2F* h2_jet3_PtEta[7];   //!
  TH2F* h2_jet3FJR_PtEta[7];   //!
  TH2F* h2_jet3FJR_PtSubtrPt[7];   //!
  TH2F* h2_jet3FJR_PtNConst[7];   //!
  
  TH2F* h2_jet2_PtEta[7];   //!
  TH2F* h2_jet2FJR_PtEta[7];   //!
  TH2F* h2_jet2FJR_PtSubtrPt[7];   //!
  TH2F* h2_jet2FJR_PtNConst[7];   //!
  
  //-- More diag plots for R=0.4 jets
  
  TH2F* h2_jet4FJR_EtaPhi[7];   //!
  TH2F* h2_jet4FJR_EtaPhiAvgNconst[7];   //!
  TH2F* h2_jet4FJR_EtaPhiAvgSubtrPt[7];   //!
  TH2F* h2_jet4FJR_EtaPhiDevSubtrPt[7];   //!
  
  // Tree
  TTree *tree; //!
  
  Bool_t event_isTriggered[10]; //!
  Bool_t trigger[10]; //!
  Int_t jet_isTriggered[10]; //!
  Float_t trig_prescale[10]; //!
  Float_t FCalEt; //!
  //UInt_t lbn_n, event_n, run_n; ... these are members of the algorithm
  
  Int_t jet4_n; //!
  vector<float> *jet4_pt, *jet4_eta, *jet4_phi, *jet4_m;	//!
  Int_t jet4Truth_n; //!
  vector<float> *jet4Truth_pt, *jet4Truth_eta, *jet4Truth_phi, *jet4Truth_m;  //!
  
  // this is a standard constructor
  InclusiveJetsEventLoop ();
  
  
  void setZero(int pos, vector<float>* v);
  
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
  ClassDef(InclusiveJetsEventLoop, 1);
  
};

#endif
