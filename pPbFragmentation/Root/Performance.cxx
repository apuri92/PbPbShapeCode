#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include "pPbFragmentation/Performance.h"
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

using namespace std;
using namespace JetHelperTools;
using namespace TrackHelperTools;
using namespace MTCorrector;

ClassImp(Performance)

#define EL_RETURN_CHECK( CONTEXT, EXP )			\
do {							\
  if( ! EXP.isSuccess() ) {				\
    Error( CONTEXT,					\
    XAOD_MESSAGE( "Failed to execute: %s" ),	\
    #EXP );					\
    return EL::StatusCode::FAILURE;			\
    }							\
    } while( false )


	Performance :: Performance ()
    {
    }
    
    EL::StatusCode Performance :: setupJob (EL::Job& job)
    {    
      // let's initialize the algorithm to use the xAODRootAccess package
      job.useXAOD ();
      
      EL_RETURN_CHECK( "setupJob()", xAOD::Init() ); // call before opening first file
      cout << " Job setup done!" << endl;
            
      switch (_jet_radius){
	  	case 2:
	  		_dR_max = 0.2;
	  	break;
	  	case 3:
	  		_dR_max = 0.3;
	  	break;
	  	case 4:
	  		_dR_max = 0.4;
	  	break;
	  	case 5:
	  		_dR_max = 0.5;
	  	break;
	  	case 6:
	  		_dR_max = 0.6;
	  	break;  
	  }

	  //_barcodeMin=0;
	  //_barcodeMax=9999;
	  //if (_isHerwig) _barcodeMax=1e5;
	  
	  nCentBins = GetCentralityNBins(_centrality_scheme);
	  cout << "Number of centrality bins: " << nCentBins <<endl; 
	  
	  std::cout << "[Performance() : initialized with jet radius = " << _jet_radius << std::endl;
	  //std::cout << "[Performance() : initialized barcodes in range " << _barcodeMin << " to " << _barcodeMax << std::endl;
      
      //Tracking cuts
      SetTrackCuts(_trk_selection,NCuts,Dod0Param,nSISHits_cut, nSIHits_cut,nSIHoles_cut,nPixHoles_cut,nSCTHits_cut,nPixHits_cut,nBLHits_cut, nIBLHits_cut,d0_cut,z0sintheta_cut,nTRTHits_cut,sig_cut);		  
      return EL::StatusCode::SUCCESS;
    }
    
    EL::StatusCode Performance :: histInitialize ()
    {   
      cout << " Setting  histograms" << endl;
      
      h_FCal_Et = new TH1D("h_FCal_Et",";FCal E_{T};N",100,0,5);
	  h_RejectionHisto = new TH1D("RejectionHisto","RejectionHisto",7,0,7);
      SetRejectionHistogram(h_RejectionHisto);
      
      Double_t zBins[1000], ptTrkBins[1000], ptJetBins[1000], phiTrkBins[1000] ,etaTrkBins[1000],MCProbBins[1000],trk_res[1000];
      Int_t ptJetBinsN = 14, ptTrkBinsN = 35, etaTrkBinsN = 50,phiTrkBinsN = 100,MCProbBinsN=10,trk_resN=2000;
      Double_t PVBins[3]={0,1,2};
      int PVBinsN=2;
      SetupBinning(0, "pt-trk", ptTrkBins, ptTrkBinsN);
      SetupBinning(0, "pt-jet", ptJetBins, ptJetBinsN);
      SetupBinning(0, "eta-trk", etaTrkBins, etaTrkBinsN);
      SetupBinning(0, "phi-trk", phiTrkBins, phiTrkBinsN);
      SetupBinning(0, "MC-prob", MCProbBins, MCProbBinsN);
      SetupBinning(0, "trk_res", trk_res, trk_resN);
      
      h_triggercounter = new TH2D("h_triggercounter","h_triggercounter",_nTriggers,0,_nTriggers,2,0,2);
      SetTrigger_hist(h_triggercounter);
      h_trk_resolution = new TH3D("h_trk_resolution","h_trk_resolution",ptTrkBinsN, ptTrkBins,etaTrkBinsN,etaTrkBins,trk_resN,trk_res);
      
      h_reco_trk_map = new TH3D("h_reco_trk_map","h_reco_trk_map;p_{T};#eta;#phi",ptTrkBinsN, ptTrkBins,etaTrkBinsN,etaTrkBins,phiTrkBinsN,phiTrkBins);
      h_reco_trk_map_nocuts = new TH3D("h_reco_trk_map_nocuts","h_reco_trk_map_nocuts;p_{T};#eta;#phi",ptTrkBinsN, ptTrkBins,etaTrkBinsN,etaTrkBins,phiTrkBinsN,phiTrkBins);
      h_truth_trk_map = new TH3D("h_truth_trk_map","h_truth_trk_map;p_{T};#eta;#phi",ptTrkBinsN, ptTrkBins,etaTrkBinsN,etaTrkBins,phiTrkBinsN,phiTrkBins);
      	   
	  TH3D* temphist_3D = nullptr;
	  
	  for (int i=0;i<GetCentralityNBins(_centrality_scheme);i++){
	   		temphist_3D = new TH3D(Form("h_reco_Injet_matched_%i",i),Form("h_reco_Injet_matched_%i",i),ptJetBinsN,ptJetBins,ptTrkBinsN, ptTrkBins,etaTrkBinsN,etaTrkBins);
	   		h_reco_Injet_matched.push_back(temphist_3D);
	   		temphist_3D = new TH3D(Form("h_eff_Injet_matched_%i",i),Form("h_eff_Injet_matched_%i",i),ptJetBinsN,ptJetBins,ptTrkBinsN, ptTrkBins,etaTrkBinsN,etaTrkBins);
	   		h_eff_Injet_matched.push_back(temphist_3D);
	   		temphist_3D = new TH3D(Form("h_eff_matched_%i",i),Form("h_eff_matched_%i",i),ptJetBinsN,ptJetBins,ptTrkBinsN, ptTrkBins,etaTrkBinsN,etaTrkBins);
	   		h_eff_matched.push_back(temphist_3D);
	   		temphist_3D = new TH3D(Form("h_eff_Injet_%i",i),Form("h_eff_Injet_%i",i),ptJetBinsN,ptJetBins,ptTrkBinsN, ptTrkBins,etaTrkBinsN,etaTrkBins);
	   		h_eff_Injet.push_back(temphist_3D);
	   		temphist_3D = new TH3D(Form("h_eff_%i",i),Form("h_eff_%i",i),ptJetBinsN,ptJetBins,ptTrkBinsN, ptTrkBins,etaTrkBinsN,etaTrkBins);
	   		h_eff.push_back(temphist_3D);
	   		temphist_3D = new TH3D(Form("h_fake_v_jet_%i",i),Form("h_fake_v_jet_%i",i),ptJetBinsN,ptJetBins,ptTrkBinsN, ptTrkBins,etaTrkBinsN,etaTrkBins);
	   	    h_fake_v_jet.push_back(temphist_3D);
	   	    temphist_3D = new TH3D(Form("h_fake_v_jet_PV_%i",i),Form("h_fake_v_jet_PV_%i",i),ptJetBinsN,ptJetBins,ptTrkBinsN, ptTrkBins,etaTrkBinsN,etaTrkBins);
	   	    h_fake_v_jet_PV.push_back(temphist_3D);
	   	    
	   	    wk()->addOutput (h_reco_Injet_matched.at(i));
	   	    wk()->addOutput (h_eff_Injet_matched.at(i));
	   	    wk()->addOutput (h_eff_matched.at(i));
	   	    wk()->addOutput (h_eff_Injet.at(i));
	   	    wk()->addOutput (h_eff.at(i));
	   	    wk()->addOutput (h_fake_v_jet.at(i));
	   	    wk()->addOutput (h_fake_v_jet_PV.at(i));
	   }
	   
	   TH2D* temphist_2D = nullptr;
	   const int n = 13;
	   std::array<string, n> CutsOn;
	   if (_trk_selection==2 || _trk_selection==6 || _trk_selection==7 ) CutsOn={"Reco","d0","z0sintheta","SI_hits","Pix_holes","SI_holes", "Shared_SI_hits","TRT_hits","Sign","Pix_hits","PixExpected","Allcuts","Gen"};
	   else if (_trk_selection<20) CutsOn={"Reco","d0","z0sintheta","SI_hits","Pix_holes","SI_holes", "Shared_SI_hits","TRT_hits","Sign","Allcuts","Gen","",""};
	   else CutsOn={"Reco","d0","z0sintheta","SCT_hits","IBL_hit","SI_holes","Pix_hits","TRT_hits","Sign","Allcuts","Gen","",""};
	   
	   for (int i=0;i<NCuts;i++){
	   	    temphist_3D = new TH3D(Form("h_cut_flow_%s",CutsOn[i].c_str()),Form("h_cut_flow_%s",CutsOn[i].c_str()),ptTrkBinsN, ptTrkBins,etaTrkBinsN,etaTrkBins,PVBinsN,PVBins);
	   		h_cut_flow.push_back(temphist_3D);  
	   		wk()->addOutput (h_cut_flow.at(i));
	   		temphist_3D = new TH3D(Form("h_cut_flow_PV_%s",CutsOn[i].c_str()),Form("h_cut_flow_PV_%s",CutsOn[i].c_str()),ptTrkBinsN, ptTrkBins,etaTrkBinsN,etaTrkBins,PVBinsN,PVBins);
	   		h_cut_flow_PV.push_back(temphist_3D);  
	   		wk()->addOutput (h_cut_flow_PV.at(i));
	   		temphist_3D = new TH3D(Form("h_cut_flow_PV_jet_%s",CutsOn[i].c_str()),Form("h_cut_flow_PV_jet_%s",CutsOn[i].c_str()),ptTrkBinsN, ptTrkBins,etaTrkBinsN,etaTrkBins,PVBinsN,PVBins);
	   		h_cut_flow_PV_jet.push_back(temphist_3D);  
	   		wk()->addOutput (h_cut_flow_PV_jet.at(i));
	   		temphist_3D = new TH3D(Form("h_cut_flow_jet_%s",CutsOn[i].c_str()),Form("h_cut_flow_jet_%s",CutsOn[i].c_str()),ptTrkBinsN, ptTrkBins,etaTrkBinsN,etaTrkBins,PVBinsN,PVBins);
	   		h_cut_flow_jet.push_back(temphist_3D);  
	   		wk()->addOutput (h_cut_flow_jet.at(i));
	  }
	  
	  Double_t finehitsBins[1000],d0z0Bins[1000];
      Int_t finehitsBinsN = 100, d0z0BinsN=600;
         
      SetupBinning(0, "hits_fine", finehitsBins, finehitsBinsN);
      SetupBinning(0, "d0z0", d0z0Bins, d0z0BinsN);
	   
	  for (int i=0;i<GetCentralityNBins(_centrality_scheme);i++){
		  string i_cent="";
		  if (_centrality_scheme>1) i_cent=Form("_%i",i);
		  temphist_3D = new TH3D(Form("h_PixHits%s",i_cent.c_str()),Form("h_PixHits%s",i_cent.c_str()),ptTrkBinsN, ptTrkBins,etaTrkBinsN, etaTrkBins,finehitsBinsN,finehitsBins);
		  h_PixHits.push_back(temphist_3D);
		  temphist_3D = new TH3D(Form("h_SCTHits%s",i_cent.c_str()),Form("h_SCTHits%s",i_cent.c_str()),ptTrkBinsN, ptTrkBins,etaTrkBinsN, etaTrkBins,finehitsBinsN,finehitsBins);
		  h_SCTHits.push_back(temphist_3D); 
		  temphist_2D = new TH2D(Form("h_d0%s",i_cent.c_str()),Form("h_d0%s",i_cent.c_str()),ptTrkBinsN,ptTrkBins,d0z0BinsN,d0z0Bins);
		  h_d0.push_back(temphist_2D);
		  temphist_2D = new TH2D(Form("h_z0sintheta%s",i_cent.c_str()),Form("h_z0sintheta%s",i_cent.c_str()),ptTrkBinsN, ptTrkBins,d0z0BinsN,d0z0Bins);
		  h_z0sintheta.push_back(temphist_2D); 
		  temphist_3D = new TH3D(Form("h_vx%s",i_cent.c_str()),"h_vx;v_{x};v_{y};v_{z}",100,-0.6,-0.4,100,-0.6,-0.4,200,-200,200);
		  h_vx.push_back(temphist_3D); 
		  wk()->addOutput (h_PixHits.at(i));
		  wk()->addOutput (h_SCTHits.at(i));
		  wk()->addOutput (h_d0.at(i));
		  wk()->addOutput (h_z0sintheta.at(i));
		  wk()->addOutput (h_vx.at(i));
      }
       
      h_mc_prob = new TH3D("h_mc_prob","h_mc_prob",ptTrkBinsN, ptTrkBins,etaTrkBinsN, etaTrkBins,d0z0BinsN,d0z0Bins);
      h_mc_prob_v_dR = new TH3D("h_mc_prob_v_dR","h_mc_prob_v_dR",ptTrkBinsN, ptTrkBins,d0z0BinsN,d0z0Bins,d0z0BinsN,d0z0Bins);
      h_mc_prob_v_dR_injet = new TH3D("h_mc_prob_v_dR_injet","h_mc_prob_v_dR_injet",ptTrkBinsN, ptTrkBins,d0z0BinsN,d0z0Bins,d0z0BinsN,d0z0Bins);
      h_mc_prob_v_dpt = new TH3D("h_mc_prob_v_dpt","h_mc_prob_v_dpt",ptTrkBinsN, ptTrkBins,d0z0BinsN,d0z0Bins,trk_resN,trk_res);
      h_mc_prob_v_dpt_injet = new TH3D("h_mc_prob_v_dpt_injet","h_mc_prob_v_dpt_injet",ptTrkBinsN, ptTrkBins,d0z0BinsN,d0z0Bins,trk_resN,trk_res);
	   

      wk()->addOutput (h_RejectionHisto);
      wk()->addOutput (h_FCal_Et);
      
      wk()->addOutput (h_triggercounter);
      wk()->addOutput (h_trk_resolution);
      wk()->addOutput (h_reco_trk_map);
      wk()->addOutput (h_reco_trk_map_nocuts);
      wk()->addOutput (h_truth_trk_map);
      
      wk()->addOutput (h_mc_prob);
      wk()->addOutput (h_mc_prob_v_dR);
      wk()->addOutput (h_mc_prob_v_dR_injet);
      wk()->addOutput (h_mc_prob_v_dpt);
      wk()->addOutput (h_mc_prob_v_dpt_injet);
      
      
      cout << " Histograms  ready" << endl;    
      cout << " tree  ready" << endl; 
      
      cout << "Parametrization of d0 cut" << endl;
      f_d0_cut = new TF1("f1", "[0]*exp([1]*x)+[2]*exp([3]*x)", 0.4, 500);
      f_d0_cut->SetParameters(0.472367, -0.149934, 0.193095, 0.000337765);
       
      return EL::StatusCode::SUCCESS;
    }
       
    EL::StatusCode Performance :: fileExecute ()
    {
      // Here you do everything that needs to be done exactly once for every
      // single file, e.g. collect a list of all lumi-blocks processed
      return EL::StatusCode::SUCCESS;
    }
    
    
    
    EL::StatusCode Performance :: changeInput (bool firstFile)
    {
      return EL::StatusCode::SUCCESS;
    }
    
    
    
    EL::StatusCode Performance :: initialize ()
    {
      // count number of events
      
      cout << " Starting initialization" << endl;
      m_eventCounter = 0;
      
      xAOD::TEvent* event = wk()->xaodEvent();
      
      // as a check, let's see the number of events in our xAOD
      Info("initialize()", "Number of events = %lli", event->getEntries() ); // print long long int
      
      // Initialize and configure trigger tools  
	  if (_data_switch==0){		  
		  m_trigConfigTool = new TrigConf::xAODConfigTool("xAODConfigTool"); // gives us access to the meta-data
		  m_trigConfigTool->msg().setLevel( MSG::ERROR ); 
		  m_trigConfigTool->initialize();
		  ToolHandle< TrigConf::ITrigConfigTool > trigConfigHandle( m_trigConfigTool );
		  m_trigDecisionTool = new Trig::TrigDecisionTool("TrigDecisionTool");
		  m_trigDecisionTool->msg().setLevel( MSG::ERROR ); 
		  m_trigDecisionTool->setProperty( "ConfigTool", trigConfigHandle ); // connect the TrigDecisionTool to the ConfigTool
		  m_trigDecisionTool->setProperty( "TrigDecisionKey", "xTrigDecision");
		  m_trigDecisionTool->initialize();
		  
		  cout << "Adding following " << _nTriggers << " triggers: ";
		  for (int i=0;i<_nTriggers;i++){
		     cout << trigger_chains.at(i) << ", ";
		    _chainGroup.push_back(m_trigDecisionTool->getChainGroup(trigger_chains.at(i)));
		  }	
		  cout << endl << "Initialize triggers finished" << endl;
      }
      
      // GRL
      TString xfn = gSystem->GetFromPipe("echo $ROOTCOREBIN");
      TString xmlfile = xfn + "/../pPbFragmentation/data/"+ _GRL;
      
      m_grl = new GoodRunsListSelectionTool("GoodRunsListSelectionTool");  
      std::vector<std::string> vecStringGRL;
      vecStringGRL.push_back(xmlfile.Data());
      EL_RETURN_CHECK("initialize()",m_grl->setProperty( "GoodRunsListVec", vecStringGRL));
      EL_RETURN_CHECK("initialize()",m_grl->setProperty("PassThrough", false)); // if true (default) will ignore result of GRL and will just pass all events
      EL_RETURN_CHECK("initialize()",m_grl->initialize());
      
      //Calibration  
      const std::string name = "pPbFragmentation"; //string describing the current thread, for logging
      TString jetAlgo = "AntiKt4HI"; //String describing your jet collection, for example AntiKt4EMTopo or AntiKt4LCTopo (see below)
      TString config = "JES_MC15CHI_042316.config"; //Path to global config used to initialize the tool (see below)
      TString calibSeq = "EtaJES_DEV"; //String describing the calibration sequence to apply (see below)
      bool isData = false; //bool describing if the events are data or from simulation
      
      //Call the constructor. The default constructor can also be used if the arguments are set with python configuration instead
      m_jetCalibration = new JetCalibrationTool(name, jetAlgo, config, calibSeq, isData);
      
      //Initialize the tool
      EL_RETURN_CHECK("initialize()",m_jetCalibration->initializeTool(name));
      
      //Jet Cleaning
      // initialize and configure the jet cleaning tool
      m_jetCleaning = new JetCleaningTool("JetCleaning");
      m_jetCleaning->msg().setLevel( MSG::DEBUG ); 
      EL_RETURN_CHECK("initialize()",m_jetCleaning->setProperty( "CutLevel", "LooseBad"));
      EL_RETURN_CHECK("initialize()",m_jetCleaning->setProperty("DoUgly", false));
      EL_RETURN_CHECK("initialize()",m_jetCleaning->initialize());

      cout << " Initialization done" << endl;
      return EL::StatusCode::SUCCESS;
    }

//Loop over events
EL::StatusCode Performance :: execute (){
   
   xAOD::TEvent* event = wk()->xaodEvent();
      
   // Event counter
   int statSize=1;
   if(m_eventCounter!=0)
   {
      double power=std::floor(log10(m_eventCounter));
      statSize=(int)std::pow(10.,power);
   }
   if(m_eventCounter%statSize==0) {
   	std::cout << "Event: " << m_eventCounter << std::endl;
   	cout << " _mcProbCut " << _mcProbCut << endl;
   }	
   m_eventCounter++;
   
   
   //All events
   h_RejectionHisto->Fill(0.5);
      
   //---------------------------
   //     Event information
   //--------------------------- 
      
   const xAOD::EventInfo* eventInfo = 0;
   EL_RETURN_CHECK("execute",event->retrieve( eventInfo, "EventInfo"));
      
   event_n = eventInfo->eventNumber();
   run_n = eventInfo->runNumber();
   int lbn = eventInfo->lumiBlock();
      
   double FCal_Et = 0;
   int cent_bin=0;
   if (_centrality_scheme>1) {
        //Centrality
        const xAOD::HIEventShapeContainer* calos = 0;
        EL_RETURN_CHECK("execute()",event->retrieve( calos, "CaloSums" ));
        FCal_Et = 0;
        int x = 0;
        xAOD::HIEventShapeContainer::const_iterator calo_itr = calos->begin();
        xAOD::HIEventShapeContainer::const_iterator calo_end = calos->end();
        for( ; calo_itr != calo_end; ++calo_itr ) {
          if (x == 5) {
            FCal_Et = ((*calo_itr)->et() * 0.001 * 0.001 );
          }
          x++;
        }
        cent_bin = GetCentralityBin(_centrality_scheme, FCal_Et);
        h_FCal_Et->Fill(FCal_Et);
   }     
   if (cent_bin < 0) {
     h_RejectionHisto->Fill(1.5);
     return EL::StatusCode::SUCCESS;
   }
      
   // check if the event is data or MC
   bool isMC = false;
   // check if the event is MC
   if(eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) ){
      isMC = true; 
    } 
      
   // GRL
   if(!isMC){ 
      if(!m_grl->passRunLB(*eventInfo)){
        h_RejectionHisto->Fill(2.5);
        return EL::StatusCode::SUCCESS; // go to next event
      }
    }
      
    //Vertex requirement
   
    const xAOD::VertexContainer * vertices = 0;
    if ( !event->retrieve( vertices, "PrimaryVertices" ).isSuccess() ){
        Error("execute()", "Failed to retrieve VertexContainer container. Exiting." );
        return EL::StatusCode::FAILURE;
    }
        
    if(vertices->size()<2) {
     	h_RejectionHisto->Fill(3.5);
     	return EL::StatusCode::SUCCESS;
    }	
    xAOD::VertexContainer::const_iterator vtx_itr = vertices->begin();
    
    xAOD::VertexContainer::const_iterator vxIter    = vertices->begin();
    xAOD::VertexContainer::const_iterator vxIterEnd = vertices->end();
    const xAOD::Vertex* primaryVertex = 0;
    for ( size_t ivtx = 0; vxIter != vxIterEnd; ++vxIter, ++ivtx ){
        // the first and only primary vertex candidate is picked
        if ( (*vxIter)->vertexType() ==  xAOD::VxType::PriVtx){
            primaryVertex = (*vxIter);
            h_vx.at(cent_bin)->Fill(primaryVertex->x(),primaryVertex->y(),primaryVertex->z());
            break;
        }
    }
    if (fabs(primaryVertex->z())>150.){
		 return EL::StatusCode::SUCCESS;
	}
    //TODO rapidity gap?
    
    //DAQ errors 
    if(!isMC){
		if(   (eventInfo->errorState(xAOD::EventInfo::LAr)==xAOD::EventInfo::Error ) || (eventInfo->errorState(xAOD::EventInfo::Tile)==xAOD::EventInfo::Error ) || (eventInfo->errorState(xAOD::EventInfo::SCT)==xAOD::EventInfo::Error ) || (eventInfo->isEventFlagBitSet(xAOD::EventInfo::Core, 18) ) ){
		  h_RejectionHisto->Fill(4.5);
		  return EL::StatusCode::SUCCESS; // go to the next event
		}
	}
    
     
     h_RejectionHisto->Fill(5.5);
    
    //Tracks    
    const xAOD::TrackParticleContainer* recoTracks = 0;
	if ( !event->retrieve( recoTracks, "InDetTrackParticles" ).isSuccess() ){
	   Error("execute()", "Failed to retrieve Reconstructed Track container. Exiting." );
	   return EL::StatusCode::FAILURE;
	}
    
    //Jet vectors
    vector<float> jet_pt_EM_vector, jet_pt_SEB_vector,jet_pt_prexcalib_vector,jet_pt_xcalib_vector,jet_phi_vector,jet_eta_vector,jet_m_vector;
	vector<int> Is_jet_Good;
	vector<float> truth_jet_eta_vector,truth_jet_m_vector,truth_jet_phi_vector,truth_jet_pt_vector;
	vector<int> truth_jet_indices;
	vector<int> hasTruth;
	vector<int> isTriggered[_nTriggers];
	vector<float> jet_NBJ_pT_vector,truth_jet_NBJ_pT_vector;
	       
    //***** Reco jets*****
    
	  xAOD::TStore store; //For calibration
	  const xAOD::JetContainer* jets = 0;
	  EL_RETURN_CHECK("execute()",event->retrieve( jets, _reco_jet_collection.c_str() ));
	      
     
	  xAOD::JetContainer::const_iterator jet_itr = jets->begin();
	  xAOD::JetContainer::const_iterator jet_end = jets->end();
	  
	  xAOD::JetContainer* updatedjets = new xAOD::JetContainer();
	  xAOD::AuxContainerBase* updatedjetsAux = new xAOD::AuxContainerBase();
	  updatedjets->setStore( updatedjetsAux );
	  
	  int jet_counter=0;
	  for( ; jet_itr != jet_end; ++jet_itr ) {
	    jet_counter++;
	    xAOD::Jet* newjet = new xAOD::Jet();
	    newjet->makePrivateStore( **jet_itr );
	    updatedjets->push_back( newjet );    
	    
	    const xAOD::JetFourMom_t jet_4mom_def = newjet->jetP4();        
	    const xAOD::JetFourMom_t jet_4mom = newjet->jetP4("JetPileupScaleMomentum");
	    
	    float uncalib_jet_pt  = (jet_4mom.pt() * 0.001);
	               
	    newjet->setJetP4("JetPileupScaleMomentum", jet_4mom);    
	    
	    EL_RETURN_CHECK("execute()", m_jetCalibration->applyCalibration( *newjet ) );
	    
	    if (_reco_jet_collection.find("HI") != std::string::npos) {
			jet_pt  = (newjet->pt() * 0.001);
			jet_eta  = newjet->eta();
			jet_phi  = newjet->phi();
			jet_m  = newjet->m()*0.001;
		}
		else{
			jet_pt  = (jet_4mom_def.pt() * 0.001);
			jet_eta  = jet_4mom_def.eta();
			jet_phi  = jet_4mom_def.phi();
			jet_m  = jet_4mom_def.M()*0.001;
		}	
	    
	    if (fabs(jet_eta)>2.1) continue;
	    
	    //Jet quality moment
		if( !m_jetCleaning->accept( **jet_itr )) continue;
			
		jet_pt_EM_vector.push_back(uncalib_jet_pt);
		jet_pt_SEB_vector.push_back(jet_pt);
		jet_pt_prexcalib_vector.push_back(jet_pt);
		jet_pt_xcalib_vector.push_back(jet_pt);
		
		jet_phi_vector.push_back(jet_phi);
		jet_eta_vector.push_back(jet_eta);
		jet_m_vector.push_back(jet_m); 
	  }
	 
	 
	 //Get Truth tracks
	 const xAOD::TruthParticleContainer * particles = 0;
	      if ( !event->retrieve( particles, "TruthParticles" ).isSuccess() ){
	              Error("execute()", "Failed to retrieve TruthParticle container. Exiting." );
	              return EL::StatusCode::FAILURE;
	      }
	 xAOD::TruthParticleContainer::const_iterator truth_itr = particles->begin();
	 xAOD::TruthParticleContainer::const_iterator truth_beggining = particles->begin();
     xAOD::TruthParticleContainer::const_iterator truth_end = particles->end();
	 
	  
	 //Get track isolation
	 vector<float> *truth_trk_dRmin = new std::vector<float>();
	 
	 //TODO Some problem is here...
 /*    for( ; truth_itr!=truth_end; ++truth_itr){								
		//get the tracks....
		float eta = (*truth_itr)->eta();
		float phi = (*truth_itr)->phi();
		float _dR_min = 999.;
		xAOD::TruthParticleContainer::const_iterator truth_itr_in = particles->begin();
        xAOD::TruthParticleContainer::const_iterator truth_end_in = particles->end();
		for( ; truth_itr_in!=truth_end_in; ++truth_itr_in){						
	         if (truth_itr==truth_itr_in) continue;
	         if( fabs((*truth_itr_in)->charge())<0.5 ) continue;
	         if( ((*truth_itr_in)->status())!=1) continue;
	            
	         int ty=getTypeTruth((*truth_itr_in)->barcode(),(*truth_itr_in)->pdgId(),(*truth_itr_in)->status(),(*truth_itr_in)->charge());
				
		     if(ty!=1) continue;
				
		     //get the tracks....
		     float eta_close = (*truth_itr_in)->eta();
			 float phi_close = (*truth_itr_in)->phi();

			 float R = DeltaR(phi,eta,phi_close,eta_close);
			 if (R<_dR_min) _dR_min = R;  
		}		
		truth_trk_dRmin->push_back(_dR_min);
     }*/
      
     	       
     //Loop over reconstructed tracks	 
	 vector<float> *trk_good_eta = new std::vector<float>();
	 vector<float> *trk_good_phi = new std::vector<float>();
	 vector<float> *trk_good_pt = new std::vector<float>();
	 vector<int> *trk_good_isMatched = new std::vector<int>();
	 vector<int> *trk_good_vertextype = new std::vector<int>();
	 vector<float> *trk_good_nHits = new std::vector<float>();
	 vector<float> *trk_nHits = new std::vector<float>();
	 vector<float> *trk_good_matched_eta = new std::vector<float>();
	 vector<float> *trk_good_matched_phi = new std::vector<float>();
	 vector<float> *trk_good_matched_pt = new std::vector<float>();
	 vector<float> *trk_good_matched_pt_reco = new std::vector<float>();
	 vector<int> *trk_good_matched_index = new std::vector<int>();
	 vector<int> *trk_good_matched_truth_index = new std::vector<int>();
	 float mcprob;
			 
	 //Get only good 	 
	 for (const auto& trk : *recoTracks) {
						
		//get the tracks....
	    float pt = trk->pt()/1000.;
		float eta = trk->eta();
		float phi = trk->phi();
			    			
		//track parameters
		int nPixHits = trk->auxdata< unsigned char >("numberOfPixelHits");
		int nPixHoles = trk->auxdata< unsigned char >("numberOfPixelHoles");
		int nPixDeadS = trk->auxdata< unsigned char >("numberOfPixelDeadSensors");
		int nShPixH= trk->auxdata< unsigned char >("numberOfPixelSharedHits");
		int nSCTHits = trk->auxdata< unsigned char >("numberOfSCTHits");				
		int nSCTHoles = trk->auxdata< unsigned char >("numberOfSCTHoles");
		int nSCTDeadS = trk->auxdata< unsigned char >("numberOfSCTDeadSensors");
		int nShSCTH = trk->auxdata< unsigned char >("numberOfSCTSharedHits");
		int nTRTHits = trk->auxdata< unsigned char >("numberOfTRTHits");  
		int nIBLHits = trk->auxdata< unsigned char >("numberOfInnermostPixelLayerHits");
		int nBLHits = trk->auxdata< unsigned char >("numberOfNextToInnermostPixelLayerHits");
		int expIBLHits = trk->auxdata< unsigned char >("expectInnermostPixelLayerHit");
		int expBLHits = trk->auxdata< unsigned char >("expectNextToInnermostPixelLayerHit");
		int nPixSh= trk->auxdata< unsigned char >("numberOfPixelSharedHits");			
		int nSCTSh = trk->auxdata< unsigned char >("numberOfSCTSharedHits");			
				
		double d0 = trk->d0();
		double theta = trk->theta();
		double z0pv=(trk->z0()+trk->vz()-(*vtx_itr)->z())*sin(theta);	// pp: trk->z0() - w.r.t. BS
				
		float ed0=sqrt(trk->definingParametersCovMatrix()(0,0));
		float ez0pv=sqrt( (trk->definingParametersCovMatrix()(1,1))*pow(sin(theta),2) +
		                 (trk->definingParametersCovMatrix()(3,3))*pow(z0pv*cos(theta),2)+
		                 2*(trk->definingParametersCovMatrix()(1,3))*fabs(sin(theta)*z0pv*cos(theta)) );
		if(ed0==0.0) ed0=1e-10;
		if(ez0pv==0.0) ez0pv=1e-10;
		
		float sigd0 = fabs(d0/ed0);
        float sigz0 = fabs(z0pv/ez0pv);
		

        // reject TRT-seeded tracks
		std::bitset < xAOD::NumberOfTrackRecoInfo > author = trk->patternRecoInfo();
		if(author[xAOD::TRTSeededTrackFinder] || author[xAOD::TRTStandalone]) continue;
		

		
		//Vertex association
		ElementLink< xAOD::VertexContainer > vtxLink = trk->auxdata<ElementLink< xAOD::VertexContainer> >("vertexLink");
		int vertex_type=0;
		if(vtxLink.isValid()) vertex_type =  (*vtxLink)->vertexType();
		
		if (vertex_type==1) {
			h_PixHits.at(cent_bin)->Fill(pt,eta,nPixHits);
			h_SCTHits.at(cent_bin)->Fill(pt,eta,nSCTHits);
			h_d0.at(cent_bin)->Fill(pt,d0);
			h_z0sintheta.at(cent_bin)->Fill(pt,z0pv);
		}
		
		//Jet association
		float dR_min = 999.;
		float jet_pt_min = 999.;
		for(unsigned int i=0; i<jet_pt_xcalib_vector.size(); i++){
			jet_pt = jet_pt_xcalib_vector.at(i);
			if (jet_pt < 45.) continue;
			jet_eta = jet_eta_vector.at(i);
			jet_phi = jet_phi_vector.at(i);
			if (fabs(jet_eta)>2.1) continue;
			float R = DeltaR(phi,eta,jet_phi,jet_eta);
			if ( (R > _dR_max) || (pt > 1.3 * jet_pt) || dR_min < R ) continue;
			dR_min = R;	jet_pt_min = jet_pt;	
		}	
	    	   
		//Truth matching
		bool isTruthMatched=false; 
		ElementLink< xAOD::TruthParticleContainer > truthLink = trk->auxdata<ElementLink< xAOD::TruthParticleContainer > >("truthParticleLink");
		mcprob=-1.;
		if(truthLink.isValid()){
			int trktype = getTypeReco((*truthLink)->barcode(),(*truthLink)->pdgId(),(*truthLink)->status(),(*truthLink)->charge(),trk->auxdata<float>("truthMatchProbability"),_mcProbCut);
			mcprob=trk->auxdata<float>("truthMatchProbability");  					   
			if ((trktype==1 || trktype==5) && mcprob > _mcProbCut){	
				isTruthMatched=true;
			}		
		}
		
		//all tracks spectra
		h_reco_trk_map_nocuts->Fill(pt,eta,phi);
			
		//tracking cuts
		h_cut_flow[0]->Fill(pt,eta,(int)isTruthMatched+0.5);
		if (vertex_type==1) h_cut_flow_PV[0]->Fill(pt,eta,(int)isTruthMatched+0.5);		   	
		if (dR_min < _dR_max) h_cut_flow_jet[0]->Fill(pt,eta,(int)isTruthMatched+0.5);
		if (vertex_type==1 && dR_min < _dR_max) h_cut_flow_PV_jet[0]->Fill(pt,eta,(int)isTruthMatched+0.5);
		   
		//cuts
		bool pased=true;
		
		//d0 cut
		if (Dod0Param) d0_cut = f_d0_cut->Eval(pt);
		if(fabs(d0) > d0_cut) pased=false; //(2,1.5,1)
		else {
			h_cut_flow[1]->Fill(pt,eta,(int)isTruthMatched+0.5);
			if (vertex_type==1) h_cut_flow_PV[1]->Fill(pt,eta,(int)isTruthMatched+0.5);	
			if (dR_min < _dR_max) h_cut_flow_jet[1]->Fill(pt,eta,(int)isTruthMatched+0.5);
			if (vertex_type==1 && dR_min < _dR_max) h_cut_flow_PV_jet[1]->Fill(pt,eta,(int)isTruthMatched+0.5);
		}	
				
		//z0 cut
		if(fabs(z0pv) > z0sintheta_cut ) pased=false; //(2,1.5,1)
		else {
			h_cut_flow[2]->Fill(pt,eta,(int)isTruthMatched+0.5);
			if (vertex_type==1) h_cut_flow_PV[2]->Fill(pt,eta,(int)isTruthMatched+0.5);	
			if (dR_min < _dR_max) h_cut_flow_jet[2]->Fill(pt,eta,(int)isTruthMatched+0.5);
			if (vertex_type==1 && dR_min < _dR_max) h_cut_flow_PV_jet[2]->Fill(pt,eta,(int)isTruthMatched+0.5);
		}
		
		//Silicon cuts
		if(nSCTHits + nPixHits < nSIHits_cut) pased=false;
		else {
			h_cut_flow[3]->Fill(pt,eta,(int)isTruthMatched+0.5);
			if (vertex_type==1) h_cut_flow_PV[3]->Fill(pt,eta,(int)isTruthMatched+0.5);	
			if (dR_min < _dR_max) h_cut_flow_jet[3]->Fill(pt,eta,(int)isTruthMatched+0.5);
			if (vertex_type==1 && dR_min < _dR_max) h_cut_flow_PV_jet[3]->Fill(pt,eta,(int)isTruthMatched+0.5);
		}	
		
		//Pixel holes
		if(nPixHoles > nPixHoles_cut) pased=false;
		else {
			h_cut_flow[4]->Fill(pt,eta,(int)isTruthMatched+0.5);
			if (vertex_type==1) h_cut_flow_PV[4]->Fill(pt,eta,(int)isTruthMatched+0.5);	
			if (dR_min < _dR_max) h_cut_flow_jet[4]->Fill(pt,eta,(int)isTruthMatched+0.5);
			if (vertex_type==1 && dR_min < _dR_max) h_cut_flow_PV_jet[4]->Fill(pt,eta,(int)isTruthMatched+0.5);
		}	
		
		//SI holes
		if(nPixHoles + nSCTHoles > nSIHoles_cut) pased=false;
		else {
			h_cut_flow[5]->Fill(pt,eta,(int)isTruthMatched+0.5);
			if (vertex_type==1) h_cut_flow_PV[5]->Fill(pt,eta,(int)isTruthMatched+0.5);	
			if (dR_min < _dR_max) h_cut_flow_jet[5]->Fill(pt,eta,(int)isTruthMatched+0.5);
			if (vertex_type==1 && dR_min < _dR_max) h_cut_flow_PV_jet[5]->Fill(pt,eta,(int)isTruthMatched+0.5);
		}
		
		//Shared hits cut
		if ((nPixSh + nSCTSh/2.)>nSISHits_cut) ;//continue;
		else {
			h_cut_flow[6]->Fill(pt,eta,(int)isTruthMatched+0.5);
			if (vertex_type==1) h_cut_flow_PV[6]->Fill(pt,eta,(int)isTruthMatched+0.5);	
			if (dR_min < _dR_max) h_cut_flow_jet[6]->Fill(pt,eta,(int)isTruthMatched+0.5);
			if (vertex_type==1 && dR_min < _dR_max) h_cut_flow_PV_jet[6]->Fill(pt,eta,(int)isTruthMatched+0.5);
		}
		
		//TRT hits
		if (nTRTHits < nTRTHits_cut) pased=false;
		else {
			h_cut_flow[7]->Fill(pt,eta,(int)isTruthMatched+0.5);
			if (vertex_type==1) h_cut_flow_PV[7]->Fill(pt,eta,(int)isTruthMatched+0.5);	
			if (dR_min < _dR_max) h_cut_flow_jet[7]->Fill(pt,eta,(int)isTruthMatched+0.5);
			if (vertex_type==1 && dR_min < _dR_max) h_cut_flow_PV_jet[7]->Fill(pt,eta,(int)isTruthMatched+0.5);
		}
		
		//Significant cuts
		if (sigd0 > sig_cut || sigz0 > sig_cut) pased=false;
		else {
			h_cut_flow[8]->Fill(pt,eta,(int)isTruthMatched+0.5);
			if (vertex_type==1) h_cut_flow_PV[8]->Fill(pt,eta,(int)isTruthMatched+0.5);	
			if (dR_min < _dR_max) h_cut_flow_jet[8]->Fill(pt,eta,(int)isTruthMatched+0.5);
			if (vertex_type==1 && dR_min < _dR_max) h_cut_flow_PV_jet[8]->Fill(pt,eta,(int)isTruthMatched+0.5);
		}
													
		if (_trk_selection==2 || (_trk_selection==6 || _trk_selection==7)) { //tight cuts
			if (fabs(eta)>1.65 && ((nSCTHits + nPixHits) < 11) ) pased=false;
			if (fabs(eta)<=1.65 && ((nSCTHits + nPixHits) < 9) ) pased=false;
			if (nIBLHits+nBLHits < 1) pased=false;
			else {
				h_cut_flow[9]->Fill(pt,eta,(int)isTruthMatched+0.5);
				if (vertex_type==1) h_cut_flow_PV[9]->Fill(pt,eta,(int)isTruthMatched+0.5);	
				if (dR_min < _dR_max) h_cut_flow_jet[9]->Fill(pt,eta,(int)isTruthMatched+0.5);
				if (vertex_type==1 && dR_min < _dR_max) h_cut_flow_PV_jet[9]->Fill(pt,eta,(int)isTruthMatched+0.5);
			}
			if (expIBLHits && nIBLHits==0) pased=false;
			else if (!expIBLHits && (expBLHits && nBLHits==0)) pased=false;
			else {
				h_cut_flow[10]->Fill(pt,eta,(int)isTruthMatched+0.5);
				if (vertex_type==1) h_cut_flow_PV[10]->Fill(pt,eta,(int)isTruthMatched+0.5);
				if (dR_min < _dR_max) h_cut_flow_jet[10]->Fill(pt,eta,(int)isTruthMatched+0.5);
				if (vertex_type==1 && dR_min < _dR_max) h_cut_flow_PV_jet[10]->Fill(pt,eta,(int)isTruthMatched+0.5);
			}	
		}
		
		if (pased) {
		    h_cut_flow[NCuts-2]->Fill(pt,eta,(int)isTruthMatched+0.5);
			if (vertex_type==1) h_cut_flow_PV[NCuts-2]->Fill(pt,eta,(int)isTruthMatched+0.5);
		   	if (dR_min < _dR_max) h_cut_flow_jet[NCuts-2]->Fill(pt,eta,(int)isTruthMatched+0.5);   	
		   	if (vertex_type==1 && dR_min < _dR_max) h_cut_flow_PV_jet[NCuts-2]->Fill(pt,eta,(int)isTruthMatched+0.5);
		   	h_reco_trk_map->Fill(pt,eta,phi);
		}
		
		//Fake study
		if(truthLink.isValid() && pased){
		    float MPS = (pt - (*truthLink)->pt()*0.001)/((*truthLink)->pt()*0.001);	
		    h_mc_prob->Fill(pt,eta,mcprob);
			h_mc_prob_v_dR->Fill(pt,mcprob,DeltaR(phi,eta,(*truthLink)->phi(),(*truthLink)->eta()));
			h_mc_prob_v_dpt->Fill(pt,mcprob,MPS);
			if (dR_min < _dR_max) {
				h_mc_prob_v_dR_injet->Fill(pt,mcprob,DeltaR(phi,eta,(*truthLink)->phi(),(*truthLink)->eta()));
				h_mc_prob_v_dpt_injet->Fill(pt,mcprob,MPS);
			}				
		}
							
		if (pt > _pTtrkCut && pased) {
			trk_good_eta->push_back( eta );
			trk_good_phi->push_back( phi );
			trk_good_pt->push_back( pt );
			trk_good_nHits->push_back(nIBLHits+nBLHits+nPixHits+nSCTHits);		
			trk_good_vertextype->push_back(vertex_type);
			
			if(isTruthMatched){	
					trk_good_matched_eta->push_back((*truthLink)->eta());
					trk_good_matched_pt->push_back((*truthLink)->pt());
					trk_good_matched_pt_reco->push_back(pt);
					trk_good_matched_phi->push_back((*truthLink)->phi());
					trk_good_isMatched->push_back(1);
					float MPS = (pt - (*truthLink)->pt()*0.001)/((*truthLink)->pt()*0.001);
					h_trk_resolution->Fill((*truthLink)->pt()*0.001,(*truthLink)->eta(),MPS);
			}
			else  trk_good_isMatched->push_back(0);				
		}	
	 }      
     //Reco jets
     bool isFirstPass = true; 
     for(unsigned int i=0; i<jet_pt_xcalib_vector.size(); i++){
		jet_pt = jet_pt_xcalib_vector.at(i);
		jet_eta = jet_eta_vector.at(i);
		jet_phi = jet_phi_vector.at(i);
	    	    
	    if (fabs(jet_eta)>2.1) continue;
		//All good tracks 	    
		for(unsigned int j = 0; j< trk_good_pt->size(); j++){
			if (trk_good_isMatched->at(j)==1) continue;
			float eta = trk_good_eta->at(j);
			float phi = trk_good_phi->at(j);
			float pt = trk_good_pt->at(j);
			//track-to-jet balance cut 3 sigma cut
			if (pt > 1.3 * jet_pt) continue;				
			float R = DeltaR(phi,eta,jet_phi,jet_eta);
			if (R > _dR_max) continue;
			h_fake_v_jet[cent_bin]->Fill(jet_pt,pt,eta);
			h_fake_v_jet[nCentBins-1]->Fill(jet_pt,pt,eta);
			if (trk_good_vertextype->at(j)==1){
				h_fake_v_jet_PV[cent_bin]->Fill(jet_pt,pt,eta);
				h_fake_v_jet_PV[nCentBins-1]->Fill(jet_pt,pt,eta);
			}	 
		}			    
	    //get the MC tracks....
        
        truth_itr = particles->begin();     	     
		for( ; truth_itr!=truth_end; ++truth_itr){						
	       if( fabs((*truth_itr)->charge())<0.5 ) continue;
	       if( ((*truth_itr)->status())!=1) continue;
	            
	       int ty=getTypeTruth((*truth_itr)->barcode(),(*truth_itr)->pdgId(),(*truth_itr)->status(),(*truth_itr)->charge());	
		   if(ty!=1 && ty!=5) continue;
				
		   //get the tracks....
		   float eta = (*truth_itr)->eta();
		   float phi = (*truth_itr)->phi();
		   float pt = (*truth_itr)->pt()/ 1000.0;
		  			
		   if (fabs(eta)>2.5) continue;
			
		   //Only tracks associated with a jet		
		   float R = DeltaR(phi,eta,jet_phi,jet_eta);	
		   if(R < _dR_max ) {
		    	h_eff_Injet[cent_bin]->Fill(jet_pt,pt,eta);
				h_eff_Injet[nCentBins-1]->Fill(jet_pt,pt,eta);	 	
			}
			if (isFirstPass){
				h_eff[cent_bin]->Fill(jet_pt,pt,eta);
				h_eff[nCentBins-1]->Fill(jet_pt,pt,eta);
				h_cut_flow[NCuts-1]->Fill(pt,eta,1.5);
				h_cut_flow_PV[NCuts-1]->Fill(pt,eta,1.5);
				h_truth_trk_map->Fill(pt,eta,phi);	
			}
			if (R < _dR_max && jet_pt > 35) {
			   	h_cut_flow_jet[NCuts-1]->Fill(pt,eta,1.5);	
				h_cut_flow_PV_jet[NCuts-1]->Fill(pt,eta,1.5);
			}
							
		 }
		 
		 //Only reco track
		 //cout <<  "test " << trk_good_matched_pt->size() << endl;
		 for(unsigned int j = 0; j< trk_good_matched_pt->size(); j++){
		 	float eta = trk_good_matched_eta->at(j);
			float phi = trk_good_matched_phi->at(j);
			float pt = trk_good_matched_pt->at(j) / 1000.0;
			float pt_reco = trk_good_matched_pt_reco->at(j);
			if (pt_reco > 1.3 * jet_pt) continue;															
		    //Only tracks associated with a jet	
			float R = DeltaR(phi,eta,jet_phi,jet_eta);	
			if(R < _dR_max) { 				
				h_eff_Injet_matched[cent_bin]->Fill(jet_pt,pt,eta);
				h_eff_Injet_matched[nCentBins-1]->Fill(jet_pt,pt,eta);
				h_reco_Injet_matched[cent_bin]->Fill(jet_pt,pt_reco,eta);
				h_reco_Injet_matched[nCentBins-1]->Fill(jet_pt,pt_reco,eta);
			}
			if (isFirstPass){
				h_eff_matched[cent_bin]->Fill(jet_pt,pt,eta);
				h_eff_matched[nCentBins-1]->Fill(jet_pt,pt,eta);
			}		 
		 }
		 isFirstPass=false;
	  }
	  delete trk_good_matched_eta;
      delete trk_good_matched_phi;
      delete trk_good_matched_pt;
      delete trk_good_matched_pt_reco;
      delete trk_good_nHits;
      delete trk_nHits;
      delete truth_trk_dRmin;
      
      delete trk_good_eta;
      delete trk_good_phi;
      delete trk_good_pt;
      delete trk_good_isMatched;
      delete trk_good_vertextype;	  		 
	 
      
      // Clear vectors	  
      
	  jet_pt_EM_vector.clear();
      jet_pt_SEB_vector.clear();
      jet_pt_prexcalib_vector.clear();
      jet_pt_xcalib_vector.clear();
      jet_phi_vector.clear();
      jet_eta_vector.clear();
      jet_m_vector.clear();
      

    
    return EL::StatusCode::SUCCESS;
}

EL::StatusCode Performance :: postExecute (){
      return EL::StatusCode::SUCCESS;
}
     
EL::StatusCode Performance :: finalize (){
      xAOD::TEvent* event = wk()->xaodEvent();
      
      // cleaning up trigger tools
	  if (_data_switch==0){	  
		  if( m_trigConfigTool ) {
		    delete m_trigConfigTool;
		    m_trigConfigTool = 0;
		  }
		  if( m_trigDecisionTool ) {
		    delete m_trigDecisionTool;
		    m_trigDecisionTool = 0;
		  }
      }
      // cleaning GRL
      if (m_grl) {
        delete m_grl;
        m_grl = 0;
      }
      
      //cleaning cleaning :)
      if( m_jetCleaning ) {
	    delete m_jetCleaning;
	    m_jetCleaning = 0;
	  }
 
      return EL::StatusCode::SUCCESS;
}
    
EL::StatusCode Performance :: histFinalize (){  
      cout<<"Events = "<< m_eventCounter<<endl;
      return EL::StatusCode::SUCCESS;
}
