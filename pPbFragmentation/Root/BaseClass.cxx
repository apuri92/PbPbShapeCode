#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include "pPbFragmentation/pPbFragmentation.h"
#include "pPbFragmentation/PbPbFragmentation.h"
#include "pPbFragmentation/PbPbFFShape.h"
#include "pPbFragmentation/SEBCorrectorTool.h"
#include "pPbCentrality/pPbMinBiasUtil.h"
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/tools/Message.h"
#include "xAODJet/JetContainer.h"
#include "xAODTrigger/JetRoIContainer.h"
#include "xAODTrigger/JetRoIAuxContainer.h"
#include <TFile.h>
#include <TSystem.h>


// this is needed to distribute the algorithm to the workers
ClassImp(BaseClass)


/// Helper macro for checking xAOD::TReturnCode return values
#define EL_RETURN_CHECK( CONTEXT, EXP )			\
  do {							\
    if( ! EXP.isSuccess() ) {				\
      Error( CONTEXT,					\
	     XAOD_MESSAGE( "Failed to execute: %s" ),	\
	     #EXP );					\
      return EL::StatusCode::FAILURE;			\
    }							\
  } while( false )


BaseClass :: BaseClass ()
{
  
  
  
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  Note that you should only put
  // the most basic initialization here, since this method will be
  // called on both the submission and the worker node.  Most of your
  // initialization code will go into histInitialize() and
  // initialize(). 
}

BaseClass :: BaseClass (const BaseClass& base)  {
	
	_data_switch = base._data_switch;
	_dataset = base._dataset;
	_isMB = base._isMB;
	_isHerwig = base._isHerwig;
	_jet_radius = base._jet_radius;
	_trkptBkgrThreshold = base._trkptBkgrThreshold;
	_reco_jet_collection = base._reco_jet_collection;
	_test_reco_jet_collection = base._test_reco_jet_collection;
	_truth_jet_collection=base._truth_jet_collection;
	_GRL=base._GRL;
	_cut_level=base._cut_level;
	_dR_max = base._dR_max;
	_centrality_scheme=base._centrality_scheme;
	_Run_Number=base._Run_Number;
	_truth_only = base._truth_only;
	_dR_truth_matching = base._dR_truth_matching;
	_trk_selection = base._trk_selection;
	_outputName=base._outputName;
	_pTtrkCut=base._pTtrkCut;
	_pTjetCut=base._pTjetCut;
	_truthpTjetCut=base._truthpTjetCut;
	_mcProbCut=base._mcProbCut;
	_doClusters=base._doClusters;
	_doSlimTree=base._doSlimTree;
	_doForward=base._doForward;
	_pt_iso=base._pt_iso;
	_JERBalancecut=base._JERBalancecut;
	_applyReweighting=base._applyReweighting;
	_jetptBkgrThreshold=base._jetptBkgrThreshold;
	_uncert_index=base._uncert_index;
	_useCharge=base._useCharge;
	_doPileupRejection=base._doPileupRejection;
	_correctTrackpT=base._correctTrackpT;
	_eff_jety=base._eff_jety;
	_UseAltzDef=base._UseAltzDef;
	_doFJR=base._doFJR;
	_maxjetdeltaR=base._maxjetdeltaR;
	_doJPRCorrection=base._doJPRCorrection;
	_PythiaPowheg = base._PythiaPowheg;
}













