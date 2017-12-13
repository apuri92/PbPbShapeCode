
#ifndef Performance_h
#define Performance_h

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
class Performance : public BaseClass{
public :
    
   
   
   //Internal configuration
   
   float _dR_max;
    
   bool event_isTriggered[10]; //!
   bool trigger[10]; //!
   int jet_isTriggered[10]; //!
   float trig_prescale[10]; //!
   
   int nSCTHits_cut;
   int nPixHits_cut;
   int nBLHits_cut;
   int nIBLHits_cut;
   int nSIHits_cut; 
   int nSIHoles_cut;
   int nPixHoles_cut;
   int nSCTHoles_cut;
   int nSISHits_cut;
   bool Dod0Param;
   int nTRTHits_cut;
   float sig_cut;
   
   int NCuts;
   
   int nCentBins;
   
   double d0_cut;
   double z0sintheta_cut;
   double trk_mc_probability_cut;
   	      
   vector<TH3D*> h_Injet_nPixSharedHits;
   vector<TH3D*> h_Injet_nSCTSharedHits;
   vector<TH3D*> h_Injet_nBLSharedHits;   
   
   vector<TH3D*> h_reco_Injet_matched;
   vector<TH3D*> h_eff_Injet_matched;
   vector<TH3D*> h_eff_matched;
   vector<TH3D*> h_eff_Injet;
   vector<TH3D*> h_eff;
   vector<TH3D*> h_fake_v_jet;
   vector<TH3D*> h_fake_v_jet_PV;
   
   TH3D* h_trk_resolution; //!
   TH3D* h_truth_trk_map; //!
   TH3D* h_mc_prob; //!
   TH3D* h_mc_prob_v_dR; //!
   TH3D* h_mc_prob_v_dR_injet; //!
   TH3D* h_mc_prob_v_dpt; //!
   TH3D* h_mc_prob_v_dpt_injet; //!

   
   
   // this is a standard constructor
  Performance ();
 
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
  ClassDef(Performance, 1);
   
};

#endif
