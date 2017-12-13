#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include "pPbFragmentation/InclusiveJetsEventLoop.h"
#include "pPbFragmentation/InclusiveJetsHelper.h"
#include "pPbFragmentation/GlobalHelper.h"
#include "pPbCentrality/pPbMinBiasUtil.h"
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/tools/Message.h"
#include "xAODJet/JetContainer.h"
#include "xAODMuon/MuonContainer.h"
#include "xAODTrigger/JetRoIContainer.h"
#include "xAODTrigger/JetRoIAuxContainer.h"
#include "xAODCaloEvent/CaloClusterContainer.h"
#include <TFile.h>
#include <TSystem.h>

using namespace std;
using namespace JetHelperTools;

ClassImp (InclusiveJetsEventLoop)

#define EL_RETURN_CHECK( CONTEXT, EXP )     \
  do {              \
    if( ! EXP.isSuccess() ) {       \
      Error( CONTEXT,         \
             XAOD_MESSAGE( "Failed to execute: %s" ),  \
             #EXP );         \
      return EL::StatusCode::FAILURE;     \
      }             \
    } while( false )


InclusiveJetsEventLoop :: InclusiveJetsEventLoop () {
  }

EL::StatusCode InclusiveJetsEventLoop :: setupJob (EL::Job& job) {
  // let's initialize the algorithm to use the xAODRootAccess package
  job.useXAOD ();

  EL_RETURN_CHECK ("setupJob()", xAOD::Init());  // call before opening first file
  cout << " Job setup done!" << endl;

  std::cout << "[InclusiveJetsEventLoop() : initialized with some setting " << std::endl;

  return EL::StatusCode::SUCCESS;
  }

EL::StatusCode InclusiveJetsEventLoop :: histInitialize () {
  cout << " Setting  histograms" << endl;

  h_FCal_Et  = new TH1D ("h_FCal_Et", ";FCal E_{T};N", 200, 0, 8);
  h_FCal_Et_withJet = new TH1D ("h_FCal_Et_withJet", ";FCal E_{T};N", 200, 0, 8);
  h_RejectionHisto = new TH1D ("RejectionHisto", "RejectionHisto", 7, 0, 7);
  SetRejectionHistogram (h_RejectionHisto);

  //<RS>
  
  h2_jet4FJR_FCalEt_subPt = new TH2F( "h2_jet4FJR_FCalEt_subPt", ";FCal E_{T};p^{unsub}_{T}-p^{sub}_{T} [GeV]", 200, 0, 8, 400, 0, 200);
  wk()->addOutput (h2_jet4FJR_FCalEt_subPt);
  h2_jet4FJR_FCalEt_subPt->Sumw2();
  
  double ptJetBins_All[] = {21.544, 23.263, 25.119, 27.122, 29.286, 31.622, 34.145, 36.869, 39.8110,  42.9869,  46.4162, 50.1190,  54.1172, 58.4344, 63.096,  68.129,  73.564,  79.433,  85.770,  92.612,  100.000,  107.977,  116.591,  125.892,  135.935,  146.779,  158.489,  171.132, 184.784, 199.525,  215.442,  232.628,  251.186,  271.224,  292.861,  316.224,  341.450,  368.689,  398.101,  429.859,  464.151,  501.178, 541.159, 584.328, 630.942, 681.275, 735.623, 794.306, 857.671, 926.091, 999.968, 1079.739, 1165.873, 1258.879, 1359.303, 1467.738 };

  double ptTruthJetBins_All[] = { 5.011, 6.309, 7.943, 9.999, 12.589, 15.848, 19.953, 25.119, 31.623, 39.811, 50.119, 63.096, 79.433, 100.000, 125.892,  158.488,  199.525,  251.186,  316.224,  398.101,  501.178,  630.944, 794.308, 999.970, 1258.883, 1584.833 };

  int ptJetBins_N = sizeof (ptJetBins_All) / sizeof (double) - 1;
  int ptTruthJetBins_N = sizeof (ptTruthJetBins_All) / sizeof (double) - 1;
  //</RS>

  
  
  h2_jet42_Pt = new TH2F ("h2_jet42_Pt", ";p_{T,jet}^{R=0.4} [GeV];p_{T,jet}^{R=0.2} [GeV];", ptJetBins_N, ptJetBins_All, ptJetBins_N, ptJetBins_All); // jet pt [GeV]
  h2_jet32_Pt = new TH2F ("h2_jet32_Pt", ";p_{T,jet}^{R=0.3} [GeV];p_{T,jet}^{R=0.2} [GeV];", ptJetBins_N, ptJetBins_All, ptJetBins_N, ptJetBins_All); // jet pt [GeV]
  h1_jet4_Pt = new TH1F ("h1_jet4_Pt", ";p_{T,jet}^{R=0.4} [GeV];", ptJetBins_N, ptJetBins_All);
  h1_jet3_Pt = new TH1F ("h1_jet3_Pt", ";p_{T,jet}^{R=0.3} [GeV];", ptJetBins_N, ptJetBins_All);
  h1_jet2_Pt = new TH1F ("h1_jet2_Pt", ";p_{T,jet}^{R=0.2} [GeV];", ptJetBins_N, ptJetBins_All);

  h1_FCalEtClPt_1 = new TH1F ("h1_FCalEtClPt_1", ";FCal E_{T} [TeV]; p_{T,cluster} [GeV];", 128, -0.4, 6);
  h1_FCalEtCl_1 = new TH1F ("h1_FCalEtCl_1", ";FCal E_{T} [TeV];", 128, -0.4, 6);
  h1_FCalEtClPt_2 = new TH1F ("h1_FCalEtClPt_2", ";FCal E_{T} [TeV]; p_{T,cluster} [GeV];", 128, -0.4, 6);
  h1_FCalEtCl_2 = new TH1F ("h1_FCalEtCl_2", ";FCal E_{T} [TeV];", 128, -0.4, 6);

  h3_cl_EtaPhiFCalEtPt_1 = new TH3F ("h3_cl_EtaPhiFCalEtPt_1", "#eta;#phi;FCal E_{T};", 21, -2.1, 2.1, 16, -TMath::Pi(), TMath::Pi(), 128, -0.4, 6);
  h3_cl_EtaPhiFCalEt_1 = new TH3F ("h3_cl_EtaPhiFCalEt_1", "#eta;#phi;FCal E_{T};", 21, -2.1, 2.1, 16, -TMath::Pi(), TMath::Pi(), 128, -0.4, 6);
  h3_cl_EtaPhiFCalEtPt_2 = new TH3F ("h3_cl_EtaPhiFCalEtPt_2", "#eta;#phi;FCal E_{T};", 21, -2.1, 2.1, 16, -TMath::Pi(), TMath::Pi(), 128, -0.4, 6);
  h3_cl_EtaPhiFCalEt_2 = new TH3F ("h3_cl_EtaPhiFCalEt_2", "#eta;#phi;FCal E_{T};", 21, -2.1, 2.1, 16, -TMath::Pi(), TMath::Pi(), 128, -0.4, 6);


  wk()->addOutput (h_FCal_Et);
  wk()->addOutput (h_FCal_Et_withJet);
  wk()->addOutput (h_RejectionHisto);
  wk()->addOutput (h2_jet42_Pt);
  wk()->addOutput (h2_jet32_Pt);
  wk()->addOutput (h1_jet4_Pt);
  wk()->addOutput (h1_jet3_Pt);
  wk()->addOutput (h1_jet2_Pt);

  h1_FCalEtClPt_1->Sumw2();
  h1_FCalEtCl_1->Sumw2();
  h1_FCalEtClPt_2->Sumw2();
  h1_FCalEtCl_2->Sumw2();

  wk()->addOutput (h1_FCalEtClPt_1);
  wk()->addOutput (h1_FCalEtCl_1);
  wk()->addOutput (h1_FCalEtClPt_2);
  wk()->addOutput (h1_FCalEtCl_2);

  h3_cl_EtaPhiFCalEtPt_1->Sumw2();
  h3_cl_EtaPhiFCalEt_1->Sumw2();
  h3_cl_EtaPhiFCalEtPt_2->Sumw2();
  h3_cl_EtaPhiFCalEt_2->Sumw2();

  wk()->addOutput (h3_cl_EtaPhiFCalEtPt_1);
  wk()->addOutput (h3_cl_EtaPhiFCalEt_1);
  wk()->addOutput (h3_cl_EtaPhiFCalEtPt_2);
  wk()->addOutput (h3_cl_EtaPhiFCalEt_2);



  int Ncent_bin = GetCentralityNBins (_centrality_scheme);

  for (int i = 0; i < Ncent_bin; i++)  {

    // basic diagnostic plots for clusters

    h2_cl_EtaPhi1[i] = new TH2F (Form ("h2_cl_EtaPhi1_c%i", i), ";#eta;#phi;", 42, -2.1, 2.1, 64, -TMath::Pi(), TMath::Pi());
    h2_cl_EtaPhi2[i] = new TH2F (Form ("h2_cl_EtaPhi2_c%i", i), ";#eta;#phi;", 42, -2.1, 2.1, 64, -TMath::Pi(), TMath::Pi());
    h2_cl_EtaPhiPt1[i] = new TH2F (Form ("h2_cl_EtaPhiPt1_c%i", i), ";#eta;#phi;", 42, -2.1, 2.1, 64, -TMath::Pi(), TMath::Pi());
    h2_cl_EtaPhiPt2[i] = new TH2F (Form ("h2_cl_EtaPhiPt2_c%i", i), ";#eta;#phi;", 42, -2.1, 2.1, 64, -TMath::Pi(), TMath::Pi());

    wk()->addOutput (h2_cl_EtaPhi1[i]); //!
    wk()->addOutput (h2_cl_EtaPhi2[i]); //!
    wk()->addOutput (h2_cl_EtaPhiPt1[i]); //!
    wk()->addOutput (h2_cl_EtaPhiPt2[i]); //!

    // basic diagnostic plots for all three jet parameters

    h2_jet4_PtEta[i] = new TH2F (Form ("h2_jet4_PtEta_c%i", i), ";p_{T,jet}^{R=0.4} [GeV];y;", ptJetBins_N, ptJetBins_All, 42, -2.1, 2.1); // jet pt [GeV]
    h2_jet4FJR_PtEta[i] = new TH2F (Form ("h2_jet4FJR_PtEta_c%i", i), ";p_{T,jet}^{R=0.4} [GeV];y;", ptJetBins_N, ptJetBins_All, 42, -2.1, 2.1); // jet pt [GeV]
    h2_jet4FJR_PtSubtrPt[i] = new TH2F (Form ("h2_jet4FJR_PtSubtrPt_c%i", i), ";p_{T,jet}^{R=0.4} [GeV];subtracted p_{T};", ptJetBins_N, ptJetBins_All, 200, 0, 200); // jet pt [GeV]
    h2_jet4FJR_PtNConst[i] = new TH2F (Form ("h2_jet4FJR_PtNConst_c%i", i), ";p_{T,jet}^{R=0.4} [GeV];N_{const};", ptJetBins_N, ptJetBins_All, 100, 0, 100); // jet pt [GeV]
    h2_jet4FJR_PtNConst_binning2[i] = new TH2F (Form ("h2_jet4FJR_PtNConst_binning2_c%i", i), ";p_{T,jet}^{R=0.4} [GeV];N_{const};", 500, 0, 1000,  100, 0, 100); // jet pt [GeV]

    h2_jet3_PtEta[i] = new TH2F (Form ("h2_jet3_PtEta_c%i", i), ";p_{T,jet}^{R=0.3} [GeV];y;", ptJetBins_N, ptJetBins_All, 42, -2.1, 2.1); // jet pt [GeV]
    h2_jet3FJR_PtEta[i] = new TH2F (Form ("h2_jet3FJR_PtEta_c%i", i), ";p_{T,jet}^{R=0.3} [GeV];y;", ptJetBins_N, ptJetBins_All, 42, -2.1, 2.1); // jet pt [GeV]
    h2_jet3FJR_PtSubtrPt[i] = new TH2F (Form ("h2_jet3FJR_PtSubtrPt_c%i", i), ";p_{T,jet}^{R=0.3} [GeV];subtracted p_{T};", ptJetBins_N, ptJetBins_All, 200, 0, 200); // jet pt [GeV]
    h2_jet3FJR_PtNConst[i] = new TH2F (Form ("h2_jet3FJR_PtNConst_c%i", i), ";p_{T,jet}^{R=0.3} [GeV];N_{const};", ptJetBins_N, ptJetBins_All, 100, 0, 100); // jet pt [GeV]

    h2_jet2_PtEta[i] = new TH2F (Form ("h2_jet2_PtEta_c%i", i), ";p_{T,jet}^{R=0.2} [GeV];y;", ptJetBins_N, ptJetBins_All, 42, -2.1, 2.1); // jet pt [GeV]
    h2_jet2FJR_PtEta[i] = new TH2F (Form ("h2_jet2FJR_PtEta_c%i", i), ";p_{T,jet}^{R=0.2} [GeV];y;", ptJetBins_N, ptJetBins_All, 42, -2.1, 2.1); // jet pt [GeV]
    h2_jet2FJR_PtSubtrPt[i] = new TH2F (Form ("h2_jet2FJR_PtSubtrPt_c%i", i), ";p_{T,jet}^{R=0.2} [GeV];subtracted p_{T};", ptJetBins_N, ptJetBins_All, 200, 0, 200); // jet pt [GeV]
    h2_jet2FJR_PtNConst[i] = new TH2F (Form ("h2_jet2FJR_PtNConst_c%i", i), ";p_{T,jet}^{R=0.2} [GeV];N_{const};", ptJetBins_N, ptJetBins_All, 100, 0, 100); // jet pt [GeV]

    h2_jet4_PtEta[i]->Sumw2();
    h2_jet4FJR_PtEta[i]->Sumw2();
    h2_jet4FJR_PtSubtrPt[i]->Sumw2();
    h2_jet4FJR_PtNConst[i]->Sumw2();
    h2_jet4FJR_PtNConst_binning2[i]->Sumw2();
    wk()->addOutput (h2_jet4_PtEta[i]);
    wk()->addOutput (h2_jet4FJR_PtEta[i]);
    wk()->addOutput (h2_jet4FJR_PtSubtrPt[i]);
    wk()->addOutput (h2_jet4FJR_PtNConst[i]);
    wk()->addOutput (h2_jet4FJR_PtNConst_binning2[i]);

    h2_jet3_PtEta[i]->Sumw2();
    h2_jet3FJR_PtEta[i]->Sumw2();
    h2_jet3FJR_PtSubtrPt[i]->Sumw2();
    h2_jet3FJR_PtNConst[i]->Sumw2();
    wk()->addOutput (h2_jet3_PtEta[i]);
    wk()->addOutput (h2_jet3FJR_PtEta[i]);
    wk()->addOutput (h2_jet3FJR_PtSubtrPt[i]);
    wk()->addOutput (h2_jet3FJR_PtNConst[i]);

    h2_jet2_PtEta[i]->Sumw2();
    h2_jet2FJR_PtEta[i]->Sumw2();
    h2_jet2FJR_PtSubtrPt[i]->Sumw2();
    h2_jet2FJR_PtNConst[i]->Sumw2();
    wk()->addOutput (h2_jet2_PtEta[i]);
    wk()->addOutput (h2_jet2FJR_PtEta[i]);
    wk()->addOutput (h2_jet2FJR_PtSubtrPt[i]);
    wk()->addOutput (h2_jet2FJR_PtNConst[i]);

    //-- More diag plots for R=0.4 jets
    h2_jet4FJR_EtaPhi[i] = new TH2F (Form ("h2_jet4FJR_EtaPhi_c%i", i), ";y;#phi;", 42, -2.1, 2.1, 64, -TMath::Pi(), TMath::Pi());
    h2_jet4FJR_EtaPhiAvgNconst[i] = new TH2F (Form ("h2_jet4FJR_AvgEtaPhiNConst_c%i", i), ";y;#phi;N_{const}", 42, -2.1, 2.1, 64, -TMath::Pi(), TMath::Pi());
    h2_jet4FJR_EtaPhiAvgSubtrPt[i] = new TH2F (Form ("h2_jet4FJR_AvgEtaPhiSubtrPt_c%i", i), ";y;#phi;<p_{T,subtr}>", 42, -2.1, 2.1, 64, -TMath::Pi(), TMath::Pi());
    h2_jet4FJR_EtaPhiDevSubtrPt[i] = new TH2F (Form ("h2_jet4FJR_DevEtaPhiSubtrPt_c%i", i), ";y;#phi;#sigma(p_{T,subtr})", 42, -2.1, 2.1, 64, -TMath::Pi(), TMath::Pi());

    h2_jet4FJR_EtaPhi[i]->Sumw2();
    h2_jet4FJR_EtaPhiAvgNconst[i]->Sumw2();
    h2_jet4FJR_EtaPhiAvgSubtrPt[i]->Sumw2();
    h2_jet4FJR_EtaPhiDevSubtrPt[i]->Sumw2();

    wk()->addOutput (h2_jet4FJR_EtaPhi[i]);
    wk()->addOutput (h2_jet4FJR_EtaPhiAvgNconst[i]);
    wk()->addOutput (h2_jet4FJR_EtaPhiAvgSubtrPt[i]);
    wk()->addOutput (h2_jet4FJR_EtaPhiDevSubtrPt[i]);

    //<RS>

    h1_NJet[i] = new TH1F (Form ("h1_NJet_c%i", i), ";jet p_{T} [GeV];NJet", ptJetBins_N, ptJetBins_All);
    h1_NJet[i]->Sumw2();
    wk()->addOutput (h1_NJet[i]);

    h1_stat[i] = new TH1F (Form ("Statistics_c%i", i), "", 10, 0, 10);
    h1_stat[i]->Sumw2();
    wk()->addOutput (h1_stat[i]);



    if (_dR_truth_matching > 0) {

      h1_pt_spect_truth[i] = new TH1F (Form ("h1_pt_spect_truth_c%i", i), ";jet p_{T}^{truth} [GeV];Effi;", ptTruthJetBins_N, ptTruthJetBins_All);
      h1_pt_spect_truth2[i] = new TH1F (Form ("h1_pt_spect_truth2_c%i", i), ";jet p_{T}^{truth} [GeV];Effi;", ptTruthJetBins_N, ptTruthJetBins_All);
      h1_pt_spect_truth_match[i] =  new TH1F (Form ("h1_pt_spect_truth_match_c%i", i), ";jet p_{T}^{truth} [GeV];Effi;", ptTruthJetBins_N, ptTruthJetBins_All);
      h1_pt_spect_truth_match2[i] = new TH1F (Form ("h1_pt_spect_truth_match2_c%i", i), ";jet p_{T}^{truth} [GeV];Effi;", ptTruthJetBins_N, ptTruthJetBins_All);

      h1_pt_spect_truth[i]->Sumw2();
      h1_pt_spect_truth2[i]->Sumw2();
      h1_pt_spect_truth_match[i]->Sumw2();
      h1_pt_spect_truth_match2[i]->Sumw2();
      wk()->addOutput (h1_pt_spect_truth[i]);
      wk()->addOutput (h1_pt_spect_truth2[i]);
      wk()->addOutput (h1_pt_spect_truth_match[i]);
      wk()->addOutput (h1_pt_spect_truth_match2[i]);


      h2_response[i] = new TH2F (Form ("h2_response_c%i", i), ";jet p_{T}^{reco} [GeV];jet p_{T}^{truth} [GeV]", ptJetBins_N, ptJetBins_All,  ptTruthJetBins_N, ptTruthJetBins_All);
      h2_response[i]->Sumw2();
      wk()->addOutput (h2_response[i]);

      h2_response2[i] = new TH2F (Form ("h2_response2_c%i", i), ";jet p_{T}^{reco} [GeV];jet p_{T}^{truth} [GeV]", ptJetBins_N, ptJetBins_All, ptTruthJetBins_N, ptTruthJetBins_All);
      h2_response2[i]->Sumw2();
      wk()->addOutput (h2_response2[i]);

      h2_response3[i] = new TH2F (Form ("h2_response3_c%i", i), ";jet p_{T}^{reco} [GeV];jet p_{T}^{truth} [GeV]", ptJetBins_N, ptJetBins_All,  ptTruthJetBins_N, ptTruthJetBins_All);
      h2_response3[i]->Sumw2();
      wk()->addOutput (h2_response3[i]);



      h2_jet4_JES[i] = new TH2F (Form ("h2_jet4_JES_c%i", i) , ";jet p_{T} [GeV];(p_{T}^{reco}-p_{T}^{truth})/p_{T}^{truth}", ptTruthJetBins_N, ptTruthJetBins_All,  1000, -2, 2);
      h2_jet3_JES[i] = new TH2F (Form ("h2_jet3_JES_c%i", i) , ";jet p_{T} [GeV];(p_{T}^{reco}-p_{T}^{truth})/p_{T}^{truth}", ptTruthJetBins_N, ptTruthJetBins_All,  1000, -2, 2);
      h2_jet2_JES[i] = new TH2F (Form ("h2_jet2_JES_c%i", i) , ";jet p_{T} [GeV];(p_{T}^{reco}-p_{T}^{truth})/p_{T}^{truth}", ptTruthJetBins_N, ptTruthJetBins_All,  1000, -2, 2);
      h2_jet4_JES_eta1[i] = new TH2F (Form ("h2_jet4_JES_eta1_c%i", i) , ";jet p_{T} [GeV];(p_{T}^{reco}-p_{T}^{truth})/p_{T}^{truth}", ptTruthJetBins_N, ptTruthJetBins_All,  1000, -2, 2);
      h2_jet4_JES_eta2[i] = new TH2F (Form ("h2_jet4_JES_eta2_c%i", i) , ";jet p_{T} [GeV];(p_{T}^{reco}-p_{T}^{truth})/p_{T}^{truth}", ptTruthJetBins_N, ptTruthJetBins_All,  1000, -2, 2);

      h2_jet4_dR[i] = new TH2F (Form ("h2_jet4_dR_c%i", i), ";jet p_{T};#Delta R",  ptTruthJetBins_N, ptTruthJetBins_All,  200, 0, 1);
      h2_jet4_deta[i] = new TH2F (Form ("h2_jet4_deta_c%i", i), ";jet p_{T};#eta",  ptTruthJetBins_N, ptTruthJetBins_All,  400, -1, 1);
      h2_jet4_dphi[i] = new TH2F (Form ("h2_jet4_dphi_c%i", i), ";jet p_{T};#phi",  ptTruthJetBins_N, ptTruthJetBins_All,  400, -1, 1);

      h2_jet4_JES[i]->Sumw2();
      h2_jet3_JES[i]->Sumw2();
      h2_jet2_JES[i]->Sumw2();
      h2_jet4_dR[i]->Sumw2();
      h2_jet4_JES_eta1[i]->Sumw2();
      h2_jet4_JES_eta2[i]->Sumw2();

      h2_jet4_deta[i]->Sumw2();
      h2_jet4_dphi[i]->Sumw2();

      wk()->addOutput (h2_jet4_JES[i]);
      wk()->addOutput (h2_jet3_JES[i]);
      wk()->addOutput (h2_jet2_JES[i]);
      wk()->addOutput (h2_jet4_JES_eta1[i]);
      wk()->addOutput (h2_jet4_JES_eta2[i]);
      wk()->addOutput (h2_jet4_dR[i]);
      wk()->addOutput (h2_jet4_deta[i]);
      wk()->addOutput (h2_jet4_dphi[i]);



//           h3_jet4_JES_eta[i] = new TH3F(Form("h3_jet4_JES_eta_c%i", i) ,";jet p_{T} [GeV];(E_{T}^{truth}-E_{T}^{reco})/E_{T}^{truth};#eta",ptJetBins_N, ptJetBins_All,  1000,-2,2,  42, -2.1, 2.1);
//           h3_jet4_JES_phi[i] = new TH3F(Form("h3_jet4_JES_phi_c%i", i) ,";jet p_{T} [GeV];(E_{T}^{truth}-E_{T}^{reco})/E_{T}^{truth};#phi",ptJetBins_N, ptJetBins_All,  1000,-2,2,  64, -TMath::Pi(), TMath::Pi());

      h3_jet4_JES_eta[i] = new TH3F (Form ("h3_jet4_JES_eta_c%i", i) , ";jet p_{T} [GeV];(p_{T}^{reco}-p_{T}^{truth})/p_{T}^{truth};#eta^{truth}", 500, 0, 1000,  1000, -2, 2,  42, -2.1, 2.1);
      h3_jet4_JES_eta[i]->GetXaxis()->Set (ptTruthJetBins_N, ptTruthJetBins_All);
      h3_jet4_JES_phi[i] = new TH3F (Form ("h3_jet4_JES_phi_c%i", i) , ";jet p_{T} [GeV];(p_{T}^{reco}-p_{T}^{truth})/p_{T}^{truth};#phi^{truth}", 500, 0, 1000,  1000, -2, 2,  64, -TMath::Pi(), TMath::Pi());
      h3_jet4_JES_phi[i]->GetXaxis()->Set (ptTruthJetBins_N, ptTruthJetBins_All);

      h2_jet4_JES_inPlane[i] = new TH2F (Form ("h2_jet4_JES_inPlane_c%i", i) , ";jet p_{T} [GeV];(p_{T}^{reco}-p_{T}^{truth})/p_{T}^{truth}", ptTruthJetBins_N, ptTruthJetBins_All,  1000, -2, 2);
      h2_jet4_JES_outPlane[i] = new TH2F (Form ("h2_jet4_JES_outPlane_c%i", i) , ";jet p_{T} [GeV];(p_{T}^{reco}-p_{T}^{truth})/p_{T}^{truth}", ptTruthJetBins_N, ptTruthJetBins_All,  1000, -2, 2);
      h2_jet4_JES_dpT_inPlane[i] = new TH2F (Form ("h2_jet4_JES_dpT_inPlane_c%i", i) , ";jet p_{T} [GeV];p_{T}^{reco}-p_{T}^{truth} [GeV]", ptTruthJetBins_N, ptTruthJetBins_All,  2000, -500, 500);
      h2_jet4_JES_dpT_outPlane[i] = new TH2F (Form ("h2_jet4_JES_dpT_outPlane_c%i", i) , ";jet p_{T} [GeV];p_{T}^{reco}-p_{T}^{truth} [GeV]", ptTruthJetBins_N, ptTruthJetBins_All,  2000, -500, 500);

      h3_jet4_JES_eta[i]->Sumw2();
      h3_jet4_JES_phi[i]->Sumw2();
      h2_jet4_JES_inPlane[i]->Sumw2();
      h2_jet4_JES_outPlane[i]->Sumw2();
      h2_jet4_JES_dpT_inPlane[i]->Sumw2();
      h2_jet4_JES_dpT_outPlane[i]->Sumw2();

      wk()->addOutput (h3_jet4_JES_eta[i]);
      wk()->addOutput (h3_jet4_JES_phi[i]);
      wk()->addOutput (h2_jet4_JES_inPlane[i]);
      wk()->addOutput (h2_jet4_JES_outPlane[i]);
      wk()->addOutput (h2_jet4_JES_dpT_inPlane[i]);
      wk()->addOutput (h2_jet4_JES_dpT_outPlane[i]);


      h2_jet4_etaPhi[i] = new TH2F (Form ("h2_jet4_etaPhi_c%i", i) , "#eta-#phi map;#eta;#phi",  42, -2.1, 2.1,  64, -TMath::Pi(), TMath::Pi());
      h2_jet4_etaPhi_weightedJES[i] = new TH2F (Form ("h2_jet4_etaPhi_weightedJES_c%i", i) , "#eta-#phi map weighted with JES;#eta;#phi",  42, -2.1, 2.1,  64, -TMath::Pi(), TMath::Pi());

      h2_jet4_etaPhi[i]->Sumw2();
      h2_jet4_etaPhi_weightedJES[i]->Sumw2();
      wk()->addOutput (h2_jet4_etaPhi[i]);
      wk()->addOutput (h2_jet4_etaPhi_weightedJES[i]);




      ////////////////////////////////////////////////////////



      h2_jet4_JES_truth[i] = new TH2F (Form ("h2_jet4_JES_truth_c%i", i) , ";jet p_{T} [GeV];(p_{T}^{reco}-p_{T}^{truth})/p_{T}^{truth}", ptTruthJetBins_N, ptTruthJetBins_All,  1000, -2, 2);
      h2_jet4_JES_truth_eta1[i] = new TH2F (Form ("h2_jet4_JES_truth_eta1_c%i", i) , ";jet p_{T} [GeV];(p_{T}^{reco}-p_{T}^{truth})/p_{T}^{truth}", ptTruthJetBins_N, ptTruthJetBins_All,  1000, -2, 2);
      h2_jet4_JES_truth_eta2[i] = new TH2F (Form ("h2_jet4_JES_truth_eta2_c%i", i) , ";jet p_{T} [GeV];(p_{T}^{reco}-p_{T}^{truth})/p_{T}^{truth}", ptTruthJetBins_N, ptTruthJetBins_All,  1000, -2, 2);

      h2_jet4_truth_dR[i] = new TH2F (Form ("h2_jet4_truth_dR_c%i", i), ";jet p_{T};#Delta R",  ptTruthJetBins_N, ptTruthJetBins_All,  200, 0, 1);
      h2_jet4_truth_deta[i] = new TH2F (Form ("h2_jet4_truth_deta_c%i", i), ";jet p_{T};#eta",  ptTruthJetBins_N, ptTruthJetBins_All,  400, -1, 1);
      h2_jet4_truth_dphi[i] = new TH2F (Form ("h2_jet4_truth_dphi_c%i", i), ";jet p_{T};#phi",  ptTruthJetBins_N, ptTruthJetBins_All,  400, -1, 1);

      h2_jet4_JES_truth[i]->Sumw2();
      h2_jet4_JES_truth_eta1[i]->Sumw2();
      h2_jet4_JES_truth_eta2[i]->Sumw2();
      h2_jet4_truth_dR[i]->Sumw2();
      h2_jet4_truth_deta[i]->Sumw2();
      h2_jet4_truth_dphi[i]->Sumw2();

      wk()->addOutput (h2_jet4_JES_truth[i]);
      wk()->addOutput (h2_jet4_JES_truth_eta1[i]);
      wk()->addOutput (h2_jet4_JES_truth_eta2[i]);
      wk()->addOutput (h2_jet4_truth_dR[i]);
      wk()->addOutput (h2_jet4_truth_deta[i]);
      wk()->addOutput (h2_jet4_truth_dphi[i]);


      ///


      h3_jet4_JES_truth_eta[i] = new TH3F (Form ("h3_jet4_JES_truth_eta_c%i", i) , ";jet p_{T} [GeV];(p_{T}^{reco}-p_{T}^{truth})/p_{T}^{truth};#eta^{truth}", 500, 0, 1000,  1000, -2, 2,  42, -2.1, 2.1);
      h3_jet4_JES_truth_eta[i]->GetXaxis()->Set (ptTruthJetBins_N, ptTruthJetBins_All);
      h3_jet4_JES_truth_phi[i] = new TH3F (Form ("h3_jet4_JES_truth_phi_c%i", i) , ";jet p_{T} [GeV];(p_{T}^{reco}-p_{T}^{truth})/p_{T}^{truth};#phi^{truth}", 500, 0, 1000,  1000, -2, 2,  64, -TMath::Pi(), TMath::Pi());
      h3_jet4_JES_truth_phi[i]->GetXaxis()->Set (ptTruthJetBins_N, ptTruthJetBins_All);

      h2_jet4_JES_truth_inPlane[i] = new TH2F (Form ("h2_jet4_JES_truth_inPlane_c%i", i) , ";jet p_{T} [GeV];(p_{T}^{reco}-p_{T}^{truth})/p_{T}^{truth}", ptTruthJetBins_N, ptTruthJetBins_All,  1000, -2, 2);
      h2_jet4_JES_truth_outPlane[i] = new TH2F (Form ("h2_jet4_JES_truth_outPlane_c%i", i) , ";jet p_{T} [GeV];(p_{T}^{reco}-p_{T}^{truth})/p_{T}^{truth}", ptTruthJetBins_N, ptTruthJetBins_All,  1000, -2, 2);
      h2_jet4_JES_truth_dpT_inPlane[i] = new TH2F (Form ("h2_jet4_JES_truth_dpT_inPlane_c%i", i) , ";jet p_{T} [GeV];p_{T}^{reco}-p_{T}^{truth} [GeV]", ptTruthJetBins_N, ptTruthJetBins_All,  2000, -500, 500);
      h2_jet4_JES_truth_dpT_outPlane[i] = new TH2F (Form ("h2_jet4_JES_truth_dpT_outPlane_c%i", i) , ";jet p_{T} [GeV];p_{T}^{reco}-p_{T}^{truth} [GeV]", ptTruthJetBins_N, ptTruthJetBins_All,  2000, -500, 500);

      h3_jet4_JES_truth_eta[i]->Sumw2();
      h3_jet4_JES_truth_phi[i]->Sumw2();
      h2_jet4_JES_truth_inPlane[i]->Sumw2();
      h2_jet4_JES_truth_outPlane[i]->Sumw2();
      h2_jet4_JES_truth_dpT_inPlane[i]->Sumw2();
      h2_jet4_JES_truth_dpT_outPlane[i]->Sumw2();

      wk()->addOutput (h3_jet4_JES_truth_eta[i]);
      wk()->addOutput (h3_jet4_JES_truth_phi[i]);
      wk()->addOutput (h2_jet4_JES_truth_inPlane[i]);
      wk()->addOutput (h2_jet4_JES_truth_outPlane[i]);
      wk()->addOutput (h2_jet4_JES_truth_dpT_inPlane[i]);
      wk()->addOutput (h2_jet4_JES_truth_dpT_outPlane[i]);


      h2_jet4_truth_etaPhi[i] = new TH2F (Form ("h2_jet4_truth_etaPhi_c%i", i) , "#eta-#phi map;#eta;#phi",  42, -2.1, 2.1,  64, -TMath::Pi(), TMath::Pi());
      h2_jet4_truth_etaPhi_weightedJES[i] = new TH2F (Form ("h2_jet4_truth_etaPhi_weightedJES_c%i", i) , "#eta-#phi map weighted with JES;#eta;#phi",  42, -2.1, 2.1,  64, -TMath::Pi(), TMath::Pi());

      h2_jet4_truth_etaPhi[i]->Sumw2();
      h2_jet4_truth_etaPhi_weightedJES[i]->Sumw2();
      wk()->addOutput (h2_jet4_truth_etaPhi[i]);
      wk()->addOutput (h2_jet4_truth_etaPhi_weightedJES[i]);













      //////////////////////////////////////////////////////

      }


    for (int n_trig = 0; n_trig < _nTriggers; ++n_trig)  {


      h1_pt_spect[i][n_trig] = new TH1F (Form ("h_pt_spect_%s_c%i", trigger_chains.at (n_trig).c_str(), i), Form ("h_pt_spect_%s", trigger_chains.at (n_trig).c_str()), 500, 0, 1000);
      h1_pt_spect_trig[i][n_trig] = new TH1F (Form ("h_pt_spect_trig_%s_c%i", trigger_chains.at (n_trig).c_str(), i), Form ("h_pt_spect_trig_%s", trigger_chains.at (n_trig).c_str()), 500, 0, 1000);


      h1_pt_spect[i][n_trig]->Sumw2();
      h1_pt_spect_trig[i][n_trig]->Sumw2();
      wk()->addOutput (h1_pt_spect[i][n_trig]);
      wk()->addOutput (h1_pt_spect_trig[i][n_trig]);
      }


    //</RS>
    }

  cout << " Histograms ready, now setting tree" << endl;
  TFile* outputFile = wk()->getOutputFile (_outputName);
  tree = new TTree ("tree", "Reco Jets for FF");
  tree->SetDirectory (outputFile);

  tree->Branch ("event_n", &event_n, "event_n/I");
  tree->Branch ("run_n", &run_n, "run_n/I");
  tree->Branch ("lbn_n", &lbn_n, "lbn_n/I");
  tree->Branch ("FCalEt", &FCalEt, "FCalEt/F");

  for (int i = 0; i < _nTriggers; i++) {
    tree->Branch (Form ("event_isTriggered_%i", i), &event_isTriggered[i], Form ("event_isTriggered_%i/O", i));
    tree->Branch (Form ("trigger_prescale_%i", i), &trig_prescale[i], Form ("trigger_prescale_%i/F", i));
    tree->Branch (Form ("jet_isTriggered_%i", i), &jet_isTriggered[i], Form ("jet_isTriggered_%i/I", i));
    }

  jet4_pt = new vector<float>;
  jet4_eta = new vector<float>;
  jet4_phi = new vector<float>;
  jet4_m = new vector<float>;

  tree->Branch ("AntiKt4HIJets_n", &jet4_n, "jet4_n/I");
  tree->Branch ("AntiKt4HIJets_pt", &jet4_pt);
  tree->Branch ("AntiKt4HIJets_eta", &jet4_eta);
  tree->Branch ("AntiKt4HIJets_phi", &jet4_phi);
  tree->Branch ("AntiKt4HIJets_m", &jet4_m);
  tree->Branch ("v2_fcal", &vN_fcal, "vN_fcal/F");



  if (_dR_truth_matching > 0) {  // here will be branches with truth jets

    jet4Truth_pt = new vector<float>;
    jet4Truth_eta = new vector<float>;
    jet4Truth_phi = new vector<float>;
    jet4Truth_m = new vector<float>;

    tree->Branch ("AntiKt4TruthJets_n", &jet4Truth_n, "jet4Truth_n/I");
    tree->Branch ("AntiKt4TruthJets_pt", &jet4Truth_pt);
    tree->Branch ("AntiKt4TruthJets_eta", &jet4Truth_eta);
    tree->Branch ("AntiKt4TruthJets_phi", &jet4Truth_phi);
    tree->Branch ("AntiKt4TruthJets_m", &jet4Truth_m);


    }

  return EL::StatusCode::SUCCESS;
  }

EL::StatusCode InclusiveJetsEventLoop :: fileExecute () {
  // Here you do everything that needs to be done exactly once for every
  // single file, e.g. collect a list of all lumi-blocks processed
  return EL::StatusCode::SUCCESS;
  }

EL::StatusCode InclusiveJetsEventLoop :: changeInput (bool firstFile) {
  return EL::StatusCode::SUCCESS;
  }

EL::StatusCode InclusiveJetsEventLoop :: initialize () {
  // count number of events

  cout << " Starting initialization" << endl;
  m_eventCounter = 0;

  xAOD::TEvent* event = wk()->xaodEvent();

  // as a check, let's see the number of events in our xAOD
  Info ("initialize()", "Number of events = %lli", event->getEntries()); // print long long int

  // Initialize and configure trigger tools
  if (_data_switch == 0) {
    m_trigConfigTool = new TrigConf::xAODConfigTool ("xAODConfigTool"); // gives us access to the meta-data
    m_trigConfigTool->msg().setLevel (MSG::ERROR);
    m_trigConfigTool->initialize();
    ToolHandle< TrigConf::ITrigConfigTool > trigConfigHandle (m_trigConfigTool);
    m_trigDecisionTool = new Trig::TrigDecisionTool ("TrigDecisionTool");
    m_trigDecisionTool->msg().setLevel (MSG::ERROR);
    m_trigDecisionTool->setProperty ("ConfigTool", trigConfigHandle);  // connect the TrigDecisionTool to the ConfigTool
    m_trigDecisionTool->setProperty ("TrigDecisionKey", "xTrigDecision");
    m_trigDecisionTool->initialize();

    cout << "Adding following " << _nTriggers << " triggers: ";

    for (int i = 0; i < _nTriggers; i++) {
      cout << trigger_chains.at (i) << ", ";
      _chainGroup.push_back (m_trigDecisionTool->getChainGroup (trigger_chains.at (i)));
      }

    cout << endl << "Initialize triggers finished" << endl;
    }

  printf ("HERE GOOD 1\n");

  // GRL
  TString xfn = gSystem->GetFromPipe ("echo $ROOTCOREBIN");
  TString xmlfile = xfn + "/../pPbFragmentation/data/" + _GRL;

  m_grl = new GoodRunsListSelectionTool ("GoodRunsListSelectionTool");
  std::vector<std::string> vecStringGRL;
  vecStringGRL.push_back (xmlfile.Data());

  EL_RETURN_CHECK ("initialize()", m_grl->setProperty ("GoodRunsListVec", vecStringGRL));
  EL_RETURN_CHECK ("initialize()", m_grl->setProperty ("PassThrough", false)); // if true (default) will ignore result of GRL and will just pass all events
  EL_RETURN_CHECK ("initialize()", m_grl->initialize());



  //Calibration
  const std::string name = "InclusiveJetsEventLoop"; //string describing the current thread, for logging
  TString jetAlgo = "AntiKt4HI"; //String describing your jet collection, for example AntiKt4EMTopo or AntiKt4LCTopo (see below)
  //mine
  //   TString config = "JES_MC15CHI_042316.config"; //Path to global config used to initialize the tool (see below)
  //Dennis
  TString config = "JES_MC15CHI_060316.config"; //Path to global config used to initialize the tool (see below)
  TString calibSeq = "EtaJES_DEV"; //String describing the calibration sequence to apply (see below)

  //TODO fix
  bool isData = true; //bool describing if the events are data or from simulation

  printf ("HERE GOOD 2\n");

  //insitu calibration
  TString jetAlgo_insitu = "AntiKt4EMTopo"; //String describing your jet collection, for example AntiKt4EMTopo or AntiKt4LCTopo (see below)
  TString config_insitu = "JES_2015dataset_recommendation_Feb2016.config"; //Path to global config used to initialize the tool (see below)
  const std::string name_insitu = "insitu"; //string describing the current thread, for logging
  TString calibSeq_insitu = "Insitu_DEV"; //String describing the calibration sequence to apply (see below)

  //Call the constructor. The default constructor can also be used if the arguments are set with python configuration instead
  m_jetCalibration = new JetCalibrationTool (name, jetAlgo, config, calibSeq, false);
  m_jetCalibration_insitu = new JetCalibrationTool (name_insitu, jetAlgo_insitu, config_insitu, calibSeq_insitu, isData);

  //Initialize the tool
  EL_RETURN_CHECK ("initialize()", m_jetCalibration->initializeTool (name));
  EL_RETURN_CHECK ("initialize()", m_jetCalibration_insitu->initializeTool (name_insitu));

  //Jet Cleaning
  // initialize and configure the jet cleaning tool
  m_jetCleaning = new JetCleaningTool ("JetCleaning");
  m_jetCleaning->msg().setLevel (MSG::DEBUG);
  EL_RETURN_CHECK ("initialize()", m_jetCleaning->setProperty ("CutLevel", "LooseBad"));
  EL_RETURN_CHECK ("initialize()", m_jetCleaning->setProperty ("DoUgly", false));
  EL_RETURN_CHECK ("initialize()", m_jetCleaning->initialize());

  printf ("HERE GOOD 3\n");

  cout << " Initialization done" << endl;
  return EL::StatusCode::SUCCESS;
  }

//Loop over events
EL::StatusCode InclusiveJetsEventLoop :: execute () {

  xAOD::TEvent* event = wk()->xaodEvent();

  // Event counter
  int statSize = 1;

  if (m_eventCounter != 0) {
    double power = std::floor (log10 (m_eventCounter));
    statSize = (int) std::pow (10., power);
    }

  if (m_eventCounter % statSize == 0)
    std::cout << "Event: " << m_eventCounter << std::endl;

  m_eventCounter++;

  h_RejectionHisto->Fill (0.5);

  //-------------------------------------
  //     Event cleaning -- begining
  //-------------------------------------

  const xAOD::EventInfo* eventInfo = 0;
  EL_RETURN_CHECK ("execute", event->retrieve (eventInfo, "EventInfo"));


  event_n = eventInfo->eventNumber();
  run_n = eventInfo->runNumber();
  lbn_n = eventInfo->lumiBlock();
  Int_t bcid = eventInfo->bcid();


  FCalEt = 0;
  int cent_bin = -1;

  if (_centrality_scheme > 1) {
    //Centrality
    const xAOD::HIEventShapeContainer* calos = 0;
    EL_RETURN_CHECK ("execute()", event->retrieve (calos, "CaloSums"));
    int x = 0;
    xAOD::HIEventShapeContainer::const_iterator calo_itr = calos->begin();
    xAOD::HIEventShapeContainer::const_iterator calo_end = calos->end();

    for (; calo_itr != calo_end; ++calo_itr) {
      if (x == 5) {
        FCalEt = ( (*calo_itr)->et() * 0.001 * 0.001);
        }

      x++;
      }

    cent_bin = GetCentralityBin (_centrality_scheme, FCalEt);
//           cout<<"cent = "<<cent_bin<<endl;
//           cout<<"fcal = "<<FCalEt<<endl;
    h_FCal_Et->Fill (FCalEt);
    }


  if (cent_bin < 0) {
    h_RejectionHisto->Fill (1.5);
    //                 cout<<"!!. Bad centrality = "<<FCalEt<<endl;
    return EL::StatusCode::SUCCESS;
    }

  //         printf("... cent_bin: %i , FCalEt: %.2f \n", cent_bin, FCalEt);




  //////////////////////////////////
  //  reaction plane
  //////////
  vN_fcal = 0;

  if (_centrality_scheme > 1) {

    // const xAOD::HIEventShapeContainer* calos = 0; initialized above but it was crashing down
    const xAOD::HIEventShapeContainer* calos2 = 0;
    EL_RETURN_CHECK ("execute()", event->retrieve (calos2, "CaloSums"));
//           float psiN_FCal=0;
    float FCalEtHelp = 0;
    unsigned int vn = 2;
    unsigned int harmonic = vn - 1; // n-1

    float harmonic_float = harmonic;

    for (const xAOD::HIEventShape * sh : *calos2)  {

      std::string summary;

      if (sh->isAvailable<std::string> ("Summary"))
        summary = sh->auxdata<std::string> ("Summary");

      if (summary.compare ("FCal") == 0)  {

        FCalEtHelp = sh->et(); //another way how to get FCal_Et
        float qx = sh->etCos().at (harmonic);
        float qy = sh->etSin().at (harmonic);
        //               psiN_FCal=std::atan2(qy,qx)/harmonic_float;
        vN_fcal = std::sqrt (qx + qx + qy * qy) / FCalEtHelp;
//               cout<<"Fcal_Et = "<<FCalEtHelp<<endl;
//               cout<<"qx = "<<qx<<endl;
//               cout<<"qy = "<<qy<<endl;
//               cout<<"phiN_FCal = "<<psiN_FCal<<endl;
//               cout<<"vN_fcal = "<<vN_fcal<<endl<<endl;

        break;
        }
      }

    }


  // check if the event is data or MC
  bool isMC = false;

  //<RS>
  if (_data_switch == 1)
    isMC = true;

  /*
   *    //BUG next line doesn't work for MC PbPb JZ2 sample
   *    //probably only for MC with data overlay
   *    if(eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) ){
   *
   *      isMC = true;
   *      _data_switch=1;
  }
  else{
  _data_switch=0;
  }
  */
  //</RS>

  // GRL
  if (!isMC) {
    if (!m_grl->passRunLB (*eventInfo)) {
      h_RejectionHisto->Fill (2.5);
      //   return EL::StatusCode::SUCCESS; // go to next event
      }
    }



  //Vertex requirement
  const xAOD::VertexContainer* vertices = 0;

  if (!event->retrieve (vertices, "PrimaryVertices").isSuccess()) {
    Error ("execute()", "Failed to retrieve VertexContainer container. Exiting.");
    return EL::StatusCode::FAILURE;
    }


  if (vertices->size() < 2) {
    //   printf("size of : %i \n", (int)vertices->size());
    h_RejectionHisto->Fill (3.5);
    return EL::StatusCode::SUCCESS;
    }

  xAOD::VertexContainer::const_iterator vtx_itr = vertices->begin();
  xAOD::VertexContainer::const_iterator vtx_end = vertices->end();
  // find primary vertex
  const xAOD::Vertex* primaryVertex = 0;

  for (; vtx_itr != vtx_end; ++vtx_itr) {
    if ( (*vtx_itr)->vertexType() == xAOD::VxType::PriVtx) {
      primaryVertex = (*vtx_itr);
      //h_vx.at(cent_bin)->Fill(primaryVertex->x(),primaryVertex->y(),primaryVertex->z());
      break;
      }
    }

  //TODO vertex position?

  //TODO rapidity gap?

  //DAQ errors
  if (!isMC) {
    if ( (eventInfo->errorState (xAOD::EventInfo::LAr) == xAOD::EventInfo::Error) || (eventInfo->errorState (xAOD::EventInfo::Tile) == xAOD::EventInfo::Error) || (eventInfo->errorState (xAOD::EventInfo::SCT) == xAOD::EventInfo::Error) || (eventInfo->isEventFlagBitSet (xAOD::EventInfo::Core, 18))) {
      h_RejectionHisto->Fill (4.5);
      return EL::StatusCode::SUCCESS; // go to the next event
      }
    }


  h_RejectionHisto->Fill (5.5);

  h1_stat[cent_bin]->Fill (1.5);

  //stupid but it doesn't work in Finalize()
  for (int i = 0; i < GetCentralityNBins (_centrality_scheme); ++i)  {
    h1_stat[i]->SetBinContent (3, GetCentralityNBins (_centrality_scheme));
    }

  /*
   *  // trigger
   *  if (_truth_only==0) {
   *
   *    for (int i=0;i<_nTriggers;i++){
   *      trig_prescale[i] = 1.0;
   *      event_isTriggered[i] = false;
   *      jet_isTriggered[i] = 0;
  }
  */


  //<RS>

  //----------------------------
  //Triggers
  //----------------------------

  if (_data_switch == 0) {

    int event_passed_trigger = 0;

    for (int i = 0; i < _nTriggers; i++) {
      event_isTriggered[i] =  _chainGroup.at (i)->isPassed();
      trig_prescale[i] =  _chainGroup.at (i)->getPrescale();

      //h_triggercounter->Fill(i, (Double_t) event_isTriggered[i]);
      if (event_isTriggered[i])
        event_passed_trigger = 1;
      }

//           if(!event_passed_trigger) return EL::StatusCode::SUCCESS; // go to next event
//           else h_RejectionHisto->Fill(6.5);
    }



  //</RS>



  //      Event cleaning -- end




  

   //  MC event weight for Powheg
   //TODO implement to *.cfg file
//     const xAOD::EventInfo* ei = 0;
//     EL_RETURN_CHECK ("execute()", event->retrieve (ei, "EventInfo"));
// 
//     float MCEventWeight = 1;
// 
//     if (ei->mcEventWeights().size() > 0)
//       MCEventWeight = ei->mcEventWeights() [0];

//  different approach
//     float weight = 1;
//     const std::vector< float > weights = eventInfo->mcEventWeights();
//     if( weights.size() > 0 ) weight = weights[0];

//     cout << "weight = " << weight << endl;














  //-------------------------------------
  //  Filling reco jets
  //-------------------------------------
  const xAOD::JetContainer* jets4 = 0;
  EL_RETURN_CHECK ("execute()", event->retrieve (jets4, "AntiKt4HIJets"));       //

  const xAOD::JetContainer* jets3 = 0;
  EL_RETURN_CHECK ("execute()", event->retrieve (jets3, "AntiKt3HIJets"));    //

  const xAOD::JetContainer* jets2 = 0;
  EL_RETURN_CHECK ("execute()", event->retrieve (jets2, "AntiKt2HIJets"));    //






  //<RS>
  //for calibration
  xAOD::TStore store;
  xAOD::JetContainer* updatedjets = new xAOD::JetContainer();
  xAOD::AuxContainerBase* updatedjetsAux = new xAOD::AuxContainerBase();
  updatedjets->setStore (updatedjetsAux);
  //</RS>

  //truth Jets
  const xAOD::JetContainer* jets4Truth = 0;
  const xAOD::JetContainer* jets3Truth = 0;
  const xAOD::JetContainer* jets2Truth = 0;

  if (_dR_truth_matching > 0)  {
    EL_RETURN_CHECK ("execute()", event->retrieve (jets4Truth, "AntiKt4TruthJets"));       //
    EL_RETURN_CHECK ("execute()", event->retrieve (jets3Truth, "AntiKt3TruthJets"));    //
    EL_RETURN_CHECK ("execute()", event->retrieve (jets2Truth, "AntiKt2TruthJets"));    //
    }

  // trkJets
  const xAOD::JetContainer* trkjets4 = 0;
  EL_RETURN_CHECK ("execute()", event->retrieve (trkjets4, "AntiKt4HITrackJets"));    //

  xAOD::JetContainer::const_iterator jet4_itr = jets4->begin();
  xAOD::JetContainer::const_iterator jet4_end = jets4->end();

  xAOD::JetContainer::const_iterator jet3_itr = jets3->begin();
  xAOD::JetContainer::const_iterator jet3_end = jets3->end();

  xAOD::JetContainer::const_iterator jet2_itr = jets2->begin();
  xAOD::JetContainer::const_iterator jet2_end = jets2->end();


  xAOD::JetContainer::const_iterator jet4Truth_itr ;
  xAOD::JetContainer::const_iterator jet4Truth_end ;

  xAOD::JetContainer::const_iterator jet3Truth_itr ;
  xAOD::JetContainer::const_iterator jet3Truth_end ;
  xAOD::JetContainer::const_iterator jet2Truth_itr ;
  xAOD::JetContainer::const_iterator jet2Truth_end ;









  if (_dR_truth_matching > 0) {
    //truth itr
    jet4Truth_itr = jets4Truth->begin();
    jet4Truth_end = jets4Truth->end();

    jet3Truth_itr = jets3Truth->begin();
    jet3Truth_end = jets3Truth->end();

    jet2Truth_itr = jets2Truth->begin();
    jet2Truth_end = jets2Truth->end();
    }

  xAOD::JetContainer::const_iterator trkjet4_itr = trkjets4->begin();
  xAOD::JetContainer::const_iterator trkjet4_end = trkjets4->end();

  jet4_pt->clear();
  jet4_eta->clear();
  jet4_phi->clear();
  jet4_m->clear();
  jet4_n = 0;

  if (_dR_truth_matching > 0)  {
    jet4Truth_pt->clear();
    jet4Truth_eta->clear();
    jet4Truth_phi->clear();
    jet4Truth_m->clear();
    jet4Truth_n = 0;
    }

  //         cout<<"swith 1= "<<_data_switch<<endl;
  Bool_t isEventWithHighPtJet = false;



  //   TEST
//   cout<<"for i"<<endl;
//   const xAOD::Jet *jet4_obj = 0;
//   for(int i = 0; i < jets4->size(); ++i)  {
//     jet4_obj = jets4->at(i);
//     cout<<jet4_obj->pt() * 0.001<<endl;
//
//   }
//
//   cout<<"ITER"<<endl;
//   for (; jet4_itr != jet4_end; ++jet4_itr) {
//     xAOD::Jet* uncalibjet = new xAOD::Jet();
//     uncalibjet->makePrivateStore (**jet4_itr);
//     const xAOD::JetFourMom_t jet_4mom = uncalibjet->jetP4 ();
//     uncalibjet->setJetP4 (jet_4mom);
// //     uncalibjet->setJetP4("JetConstitScaleMomentum",jet_4mom);
//
//     cout<< (uncalibjet->pt() * 0.001)<<endl;;
//   }
//   cout<<"END"<<endl<<endl;




  //**********  RECO JETS  *********//
  //JetPileupScaleMomentum
  for (; jet4_itr != jet4_end; ++jet4_itr) {

    

    xAOD::Jet* uncalibjet = new xAOD::Jet();
    uncalibjet->makePrivateStore (**jet4_itr);
    const xAOD::JetFourMom_t jet_4mom = uncalibjet->jetP4 ("JetSubtractedScaleMomentum");


    //OLD
    uncalibjet->setJetP4 ("JetConstitScaleMomentum", jet_4mom);
    uncalibjet->setJetP4 ("JetPileupScaleMomentum", jet_4mom);

    //DF
//     uncalibjet->setJetP4 ("JetSubtractedScaleMomentum", jet_4mom);

    EL_RETURN_CHECK ("execute()", m_jetCalibration->applyCalibration (*uncalibjet));



    Float_t jetPtSubtracted   = (*jet4_itr)->auxdata<float> ("JetSubtractedScaleMomentum_pt");
    Float_t jetPtUnsubtracted = (*jet4_itr)->auxdata<float> ("JetUnsubtractedScaleMomentum_pt");
    Float_t subtrPt = (jetPtUnsubtracted - jetPtSubtracted) * 0.001;
    Int_t nconst = (*jet4_itr)->numConstituents();

    Float_t jetPt  = (uncalibjet->pt() * 0.001);
    Float_t jetEta  = uncalibjet->eta();
    Float_t jetPhi  = uncalibjet->phi();
    Float_t jetM = uncalibjet->m() * 0.001;

    delete uncalibjet;
//           cout<<"After Calib ******* "<<endl;
//           cout<<"pt = "<<jetPt<<endl;
//           cout<<"eta = "<<jetEta<<endl;
//           cout<<"phi = "<<jetPhi<<endl;
//           cout<<"m = "<<jetM<<endl<<endl;


    jet4_n++;
    jet4_pt->push_back (jetPt);
    jet4_eta->push_back (jetEta);
    jet4_phi->push_back (jetPhi);
    jet4_m->push_back (jetM);




    //end of calibration

    h1_NJet[cent_bin]->Fill (jetPt);


    //          trigger effi
    for (int n_trig = 0; n_trig < _nTriggers; ++ n_trig)  {

//             cout<<"prescale "<<n_trig<<" = "<<trig_prescale[n_trig]<<endl;
      h1_pt_spect[cent_bin][n_trig]->Fill (jetPt, trig_prescale[n_trig]);

      if (event_isTriggered[n_trig])  {
        h1_pt_spect_trig[cent_bin][n_trig]->Fill (jetPt, trig_prescale[n_trig]);
        }
      }

    //</RS>



    h1_jet4_Pt->Fill (jetPt);
    h2_jet4_PtEta[cent_bin]->Fill (jetPt, jetEta);

    jet2_itr = jets2->begin();

    for (; jet2_itr != jet2_end; ++jet2_itr) {
      if (DeltaR (jetPhi, jetEta, (*jet2_itr)->phi(), (*jet2_itr)->eta()) > 0.3)
        continue;

      h2_jet42_Pt->Fill (jetPt, (*jet2_itr)->pt() * 0.001);
      }

    // FJR using track jets
    Bool_t isGoodJet = false;

    for (; trkjet4_itr != trkjet4_end; ++trkjet4_itr) {
      if ( (*trkjet4_itr)->pt() * 0.001 < 4)
        continue;  // for the first run it is 7 GeV

      if (DeltaR (jetPhi, jetEta, (*trkjet4_itr)->phi(), (*trkjet4_itr)->eta()) > 0.3)
        continue;

      isGoodJet = true;
      h2_jet4FJR_PtEta[cent_bin]->Fill (jetPt, jetEta);
      h2_jet4FJR_PtSubtrPt[cent_bin]->Fill (jetPt, subtrPt);
      
      if( (fabs(jetEta) < 0.3) && (jetPt > 60.))
        h2_jet4FJR_FCalEt_subPt->Fill(FCalEt, subtrPt);
      
      h2_jet4FJR_PtNConst[cent_bin]->Fill (jetPt, nconst);
      h2_jet4FJR_PtNConst_binning2[cent_bin]->Fill (jetPt, nconst);

      if (jetPt > 100) {
        isEventWithHighPtJet = true;
        h2_jet4FJR_EtaPhi[cent_bin]->Fill (jetEta, jetPhi);
        h2_jet4FJR_EtaPhiAvgNconst[cent_bin]->Fill (jetEta, jetPhi, nconst);
        h2_jet4FJR_EtaPhiAvgSubtrPt[cent_bin]->Fill (jetEta, jetPhi, subtrPt);
        h2_jet4FJR_EtaPhiDevSubtrPt[cent_bin]->Fill (jetEta, jetPhi, subtrPt * subtrPt);
        }

      }

    //BEGIN JES
    if (_dR_truth_matching > 0)  {   //only for MC

      Float_t dR_min = _dR_truth_matching;
      Float_t pT_min = -999.;
      Float_t eta_min = -999.;
      Float_t phi_min = -999.;
      Float_t deta_min = -999.;
      Float_t dphi_min = -999.;

//       cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;
      for (; jet4Truth_itr != jet4Truth_end; ++jet4Truth_itr)  {




        Float_t jetPtTruth = (*jet4Truth_itr)->pt() * 0.001; // GeV
        Float_t jetEtaTruth = (*jet4Truth_itr)->eta();
        Float_t jetPhiTruth = (*jet4Truth_itr)->phi();
        Float_t jetMTruth = (*jet4Truth_itr)->m() * 0.001;


        h1_pt_spect_truth2[cent_bin]->Fill (jetPtTruth);

        Float_t dR = JetHelperTools::DeltaR (jetPhi, jetEta, jetPhiTruth, jetEtaTruth);
        //                     cout<<"dR = "<<dR<<endl;

        h2_response3[cent_bin]->Fill (jetPt, jetPtTruth);

        if (dR < _dR_truth_matching)
          h2_response[cent_bin]->Fill (jetPt, jetPtTruth);

//         cout << "dR_ALL = " << dR << endl;
        if (dR < dR_min)    {
          dR_min = dR;
          deta_min = jetEta - jetEtaTruth;
          dphi_min = JetHelperTools::DeltaPhi (jetPhi, jetPhiTruth);
          pT_min = jetPtTruth;
          eta_min = jetEtaTruth;
          phi_min = jetPhiTruth;

          h1_pt_spect_truth_match2[cent_bin]->Fill (jetPtTruth);


//           cout << "dR_min = " << dR << endl;
//           cout << "JES = " << (jetPt - pT_min) / pT_min << endl;
//           cout << "pT_min = " << pT_min << endl << endl;
//           cout << " *************************" << endl;
          }
        } //end for loop truth jets



      if (dR_min < _dR_truth_matching) {



//         cout << "cent = " << cent_bin << endl;
//         cout << "jetPt = " << jetPt << endl;
//         cout << "truthPt = " << pT_min << endl;
//
//         cout << "jetEta = " << jetEta << endl;
//         cout << "truthEta = " << eta_min << endl;
//         cout << "jetPhi = " << jetPhi << endl;
//         cout << "truthPhi = " << phi_min << endl;
//         cout << "dRmin = " << dR_min << endl;
//         cout << "JES = " << (jetPt - pT_min) / pT_min << endl << endl;



        h2_jet4_dR[cent_bin]->Fill (pT_min, dR_min);
        h2_jet4_deta[cent_bin]->Fill (pT_min, deta_min);
        h2_jet4_dphi[cent_bin]->Fill (pT_min, dphi_min);
        h2_jet4_JES[cent_bin]->Fill (pT_min, (jetPt - pT_min) / pT_min);

        if (fabs (eta_min) < 1.0)
          h2_jet4_JES_eta1[cent_bin]->Fill (pT_min, (jetPt - pT_min) / pT_min);
        else
          h2_jet4_JES_eta2[cent_bin]->Fill (pT_min, (jetPt - pT_min) / pT_min);



        h1_pt_spect_truth_match[cent_bin]->Fill (pT_min);
//           h2_response2[cent_bin]->Fill (jetPt, jetPtTruth);
        h2_response2[cent_bin]->Fill (jetPt, pT_min);

        h3_jet4_JES_eta[cent_bin]->Fill (pT_min, (jetPt - pT_min) / pT_min, eta_min);
        h3_jet4_JES_phi[cent_bin]->Fill (pT_min, (jetPt - pT_min) / pT_min, phi_min);

        if (pT_min > 50.)  {
          h2_jet4_etaPhi[cent_bin]->Fill (jetEta, jetPhi, 1.);
          h2_jet4_etaPhi_weightedJES[cent_bin]->Fill (jetEta, jetPhi, (jetPt - pT_min) / pT_min);
          }

        double deltaPsi = JetHelperTools::DeltaPsi (jetPhi, vN_fcal);

        if ( (deltaPsi > 0) && (deltaPsi < TMath::Pi() / 4))         {
          h2_jet4_JES_inPlane[cent_bin]->Fill (pT_min, (jetPt - pT_min) / pT_min);
          h2_jet4_JES_dpT_inPlane[cent_bin]->Fill (pT_min, (jetPt - pT_min));
          }
        else
          if ( (deltaPsi > TMath::Pi() / 4) && (deltaPsi < TMath::Pi() / 2))  {
            h2_jet4_JES_outPlane[cent_bin]->Fill (pT_min, (jetPt - pT_min) / pT_min);
            h2_jet4_JES_dpT_outPlane[cent_bin]->Fill (pT_min, (jetPt - pT_min));

            }
        }

      }

    }

  //END for loop over jets


  //BEGIN truthJets

  if (_dR_truth_matching > 0) {
    //truth itr
    jet4Truth_itr = jets4Truth->begin();
    jet4Truth_end = jets4Truth->end();

    for (; jet4Truth_itr != jet4Truth_end; ++jet4Truth_itr)  {

      Float_t jetPtTruth = (*jet4Truth_itr)->pt() * 0.001; // GeV
      Float_t jetEtaTruth = (*jet4Truth_itr)->eta();
      Float_t jetPhiTruth = (*jet4Truth_itr)->phi();
      Float_t jetMTruth = (*jet4Truth_itr)->m() * 0.001;

      h1_pt_spect_truth[cent_bin]->Fill (jetPtTruth);

      //WARNING check if it is on the right place
      jet4Truth_n++;
      jet4Truth_pt->push_back (jetPtTruth);
      jet4Truth_eta->push_back (jetEtaTruth);
      jet4Truth_phi->push_back (jetPhiTruth);
      jet4Truth_m->push_back (jetMTruth);

      }
    }








  /*************************************************************************************************
   *
   *                   Jet R = 0.3
   *
   *
   */



  jet3_itr = jets3->begin();
  trkjet4_itr = trkjets4->begin();

  for (; jet3_itr != jet3_end; ++jet3_itr) {

    Float_t jetPt = (*jet3_itr)->pt() * 0.001; // GeV
    Float_t jetEta = (*jet3_itr)->eta();
    Float_t jetPhi = (*jet3_itr)->phi();

    Float_t jetPtSubtracted   = (*jet3_itr)->auxdata<float> ("JetSubtractedScaleMomentum_pt");
    Float_t jetPtUnsubtracted = (*jet3_itr)->auxdata<float> ("JetUnsubtractedScaleMomentum_pt");
    Float_t subtrPt = (jetPtUnsubtracted - jetPtSubtracted) * 0.001;
    Int_t nconst = (*jet3_itr)->numConstituents();


    h1_jet3_Pt->Fill (jetPt);
    h2_jet3_PtEta[cent_bin]->Fill (jetPt, jetEta);

    jet2_itr = jets2->begin();

    for (; jet2_itr != jet2_end; ++jet2_itr) {
      if (DeltaR (jetPhi, jetEta, (*jet2_itr)->phi(), (*jet2_itr)->eta()) > 0.3)
        continue;

      h2_jet32_Pt->Fill (jetPt, (*jet2_itr)->pt() * 0.001);
      }

    // FJR using track jets
    for (; trkjet4_itr != trkjet4_end; ++trkjet4_itr) {
      if ( (*trkjet4_itr)->pt() * 0.001 < 7)
        continue;

      if (DeltaR (jetPhi, jetEta, (*trkjet4_itr)->phi(), (*trkjet4_itr)->eta()) > 0.3)
        continue;

      h2_jet3FJR_PtEta[cent_bin]->Fill (jetPt, jetEta);
      h2_jet3FJR_PtSubtrPt[cent_bin]->Fill (jetPt, subtrPt);
      h2_jet3FJR_PtNConst[cent_bin]->Fill (jetPt, nconst);
      }

    } // end for loop over jets

    
    

  jet2_itr = jets2->begin();
  trkjet4_itr = trkjets4->begin();

  for (; jet2_itr != jet2_end; ++jet2_itr) {

    Float_t jetPt = (*jet2_itr)->pt() * 0.001; // GeV
    Float_t jetEta = (*jet2_itr)->eta();
    Float_t jetPhi = (*jet2_itr)->phi();

    Float_t jetPtSubtracted   = (*jet2_itr)->auxdata<float> ("JetSubtractedScaleMomentum_pt");
    Float_t jetPtUnsubtracted = (*jet2_itr)->auxdata<float> ("JetUnsubtractedScaleMomentum_pt");
    Float_t subtrPt = (jetPtUnsubtracted - jetPtSubtracted) * 0.001;
    Int_t nconst = (*jet2_itr)->numConstituents();

    h1_jet2_Pt->Fill (jetPt);
    h2_jet2_PtEta[cent_bin]->Fill (jetPt, jetEta);

    // FJR using track jets
    for (; trkjet4_itr != trkjet4_end; ++trkjet4_itr) {
      if ( (*trkjet4_itr)->pt() * 0.001 < 7)
        continue;

      if (DeltaR (jetPhi, jetEta, (*trkjet4_itr)->phi(), (*trkjet4_itr)->eta()) > 0.3)
        continue;

      h2_jet2FJR_PtEta[cent_bin]->Fill (jetPt, jetEta);
      h2_jet2FJR_PtSubtrPt[cent_bin]->Fill (jetPt, subtrPt);
      h2_jet2FJR_PtNConst[cent_bin]->Fill (jetPt, nconst);
      }

    }

  // Fill the FCal Et if there was a high-pt jet
  if (isEventWithHighPtJet) {
    //    printf("Now filling the histograms \n");
    //    cout << FCalEt << endl;
    h_FCal_Et_withJet->Fill (FCalEt);
    }

  //-------------------------------------
  //  Filling clusters
  //-------------------------------------

  /*
    const xAOD::CaloClusterContainer* cls = 0;

    EL_RETURN_CHECK ("execute()", event->retrieve (cls, "HIClusters"));

    xAOD::CaloClusterContainer::const_iterator cl_itr = cls->begin();
    xAOD::CaloClusterContainer::const_iterator cl_end = cls->end();

    for (; cl_itr != cl_end; ++cl_itr) {
      float clPt = (*cl_itr)->pt() * 0.001;
      float clEta = (*cl_itr)->rawEta();
      float clPhi = (*cl_itr)->rawPhi();

      if ( (isEventWithHighPtJet) && ( (run_n == 286717) || (run_n == 286711) || (run_n == 287560))) {
        if ( (bcid == 1) || (bcid == 454)) {
          h2_cl_EtaPhi1[cent_bin]->Fill (clEta, clPhi);
          h2_cl_EtaPhiPt1[cent_bin]->Fill (clEta, clPhi, clPt);
          h1_FCalEtClPt_1->Fill (FCalEt, clPt);
          h1_FCalEtCl_1->Fill (FCalEt);
          h3_cl_EtaPhiFCalEtPt_1->Fill (clEta, clPhi, FCalEt, clPt);
          h3_cl_EtaPhiFCalEt_1->Fill (clEta, clPhi, FCalEt);
          }
        else {
          h2_cl_EtaPhi2[cent_bin]->Fill (clEta, clPhi);
          h2_cl_EtaPhiPt2[cent_bin]->Fill (clEta, clPhi, clPt);
          h1_FCalEtClPt_2->Fill (FCalEt, clPt);
          h1_FCalEtCl_2->Fill (FCalEt);
          h3_cl_EtaPhiFCalEtPt_2->Fill (clEta, clPhi, FCalEt, clPt);
          h3_cl_EtaPhiFCalEt_2->Fill (clEta, clPhi, FCalEt);
          }
        }
      }

  */

//   cout << "Fcal = " << FCalEt << endl;
  tree->Fill();


  //
  //     iterating over truth jets, finding the closest reco jet and book JES
  //     isolation to truth has to by applied
  //

  //    reco stored in vector jet4_pt, etc., number of jets in jet4_n
  //    truth in vector jetTruth_pt, ...   ,   --||--          jet4_Truth_n
  //    one cannot change this vector before  <tree->Fill();>


  //BEGIN of isolation

  // for truth jest only, in there is a jet in dR < 0.6 with pT > 7 GeV then reject both of them
  bool appyIsolation = true;

  if (appyIsolation)  {



    Int_t  njets = 0.;

//     Float_t deltaRMin = 2.5 * 0.4;

    Float_t deltaRMin = 0.6;
    Float_t pT_minCut = 7;
    njets = jet4Truth_n;
    //  vector for booking jets that are isolated, one jet can be store more then once
    //  this will cause a problem when removing jets for e.g. jet4Truth_pt vector
    std::vector<int> Isolated_pos;
    Isolated_pos.clear();



    for (int i = 0; i < njets - 1; i++)  {

      if (jet4Truth_pt->at (i) < 0)  {
//       cout<<"recon->g_pt("<<i<<") < 0  "<<recon->g_pt(i)<<endl;
        continue;
        }

      for (int j = i + 1; j < njets; j++)  {

        if (jet4Truth_pt->at (i) < 0) {
//          cout<<"recon->g_pt("<<j<<") < 0  "<<recon->g_pt(j)<<endl;
          continue;
          }

        float R = JetHelperTools::DeltaR (jet4Truth_phi->at (i), jet4Truth_eta->at (i), jet4Truth_phi->at (j), jet4Truth_eta->at (j));



        if (R < deltaRMin && jet4Truth_pt->at (j) > pT_minCut)  {
//           cout << "i = " << i << endl;
//           cout << "j = " << j << endl;
//           cout << "R = " << R << endl;
//           cout << "jet4Truth_pt[j] = " << jet4Truth_pt->at (j) << endl;


          Isolated_pos.push_back (j);
          Isolated_pos.push_back (i);


          }

        }
      }

//     cout << "size = " << Isolated_pos.size() << endl;
//
//     for (auto & i : Isolated_pos)  {
//       cout << "pos = " << i << endl;
//       }
//
//
//
//     cout << "   Before" << endl;
//
//     for (unsigned int i = 0; i < jet4Truth_pt->size(); ++i)  {
//       cout << jet4Truth_pt->at (i) << endl;
//       }

    for (unsigned int i = 0; i < Isolated_pos.size(); ++i)  {


      setZero (Isolated_pos.at (i), jet4Truth_pt);
      setZero (Isolated_pos.at (i), jet4Truth_eta);
      setZero (Isolated_pos.at (i), jet4Truth_phi);
      setZero (Isolated_pos.at (i), jet4Truth_m);

      // we do not modify jet4Truth_n
      }


//     cout << "   After" << endl;
//
//     for (unsigned int i = 0; i < jet4Truth_pt->size(); ++i)  {
//       cout << jet4Truth_pt->at (i) << endl;
//       }
//
//     cout << "***************************************" << endl;



    }



  //BEGIN of JES, for truth jet find the closest reco jet and book JES

  for (int i_truth = 0; i_truth < jet4Truth_n; ++i_truth)  {

    if (jet4Truth_pt < 0)
      continue;

//     cout<<endl<<"NEW "<<endl<<endl;
    
    
    float dR_min = _dR_truth_matching;
    //for reco jets
    float pT_min = -999;
    float deta_min = -999;
    float dphi_min = -999;
    float eta_min = -999;
    float phi_min = -999;

    float jetTruthPt = jet4Truth_pt->at (i_truth);
    float jetTruthEta = jet4Truth_eta->at (i_truth);
    float jetTruthPhi = jet4Truth_phi->at (i_truth);




    for (int i_reco = 0; i_reco < jet4_n; ++i_reco)  {

      float dR = JetHelperTools::DeltaR (jet4Truth_phi->at (i_truth), jet4Truth_eta->at (i_truth), jet4_phi->at (i_reco), jet4_eta->at (i_reco));

//       cout << "i truth = " << i_truth << endl;
//       cout << "i reco = " << i_reco << endl;
//       cout << "dR = " << dR << endl;
//       cout << "dR_min = " << dR << endl;

      if (dR < dR_min)  {





        dR_min = dR;
        pT_min = jet4_pt->at (i_reco);

        //

        deta_min = jet4_eta->at (i_reco)  - jet4Truth_eta->at (i_truth);
        dphi_min = JetHelperTools::DeltaPhi (jet4_phi->at (i_reco), jet4Truth_phi->at (i_truth));

        eta_min = jet4_eta->at (i_reco);
        phi_min = jet4_phi->at (i_reco);
//         cout << endl << " *********************************" << endl;
//         cout << "dR < dR_min" << endl;
//         cout << "dR_min = " << dR << endl;
//         cout << "JES = " << (jetTruthPt - pT_min) / jetTruthPt << endl;
//         cout << "pT_min = " << pT_min << endl;
//         cout << "pT_truth = " << jetTruthPt << endl;
//         cout << " *********************************" << endl << endl;

        }

      }

    //END of loop ever reco jets

    //fill JES
    if (dR_min < _dR_truth_matching) {
      float JES = (jetTruthPt - pT_min) / jetTruthPt;

      h2_jet4_JES_truth[cent_bin]->Fill (jetTruthPt, JES);

      h2_jet4_truth_dR[cent_bin]->Fill (jetTruthPt, dR_min);
      h2_jet4_truth_deta[cent_bin]->Fill (jetTruthPt, deta_min);
      h2_jet4_truth_dphi[cent_bin]->Fill (jetTruthPt, dphi_min);


      if (fabs (jetTruthEta) < 1.0)
        h2_jet4_JES_truth_eta1[cent_bin]->Fill (jetTruthPt, JES);
      else
        h2_jet4_JES_truth_eta2[cent_bin]->Fill (jetTruthPt, JES);

      
//         h1_pt_spect_truth_match[cent_bin]->Fill (pT_min);
//         h2_response2[cent_bin]->Fill (jetPt, jetPtTruth);
//         h2_response2[cent_bin]->Fill (jetTruthPt, pT_min);

      h3_jet4_JES_truth_eta[cent_bin]->Fill (jetTruthPt, JES, jetTruthEta);
      h3_jet4_JES_truth_phi[cent_bin]->Fill (jetTruthPt, JES, jetTruthPhi);

      if (jetTruthPt > 50.)  {
        h2_jet4_truth_etaPhi[cent_bin]->Fill (jetTruthEta, jetTruthPhi, 1.);
        h2_jet4_truth_etaPhi_weightedJES[cent_bin]->Fill (jetTruthEta, jetTruthPhi, JES);
        }

      double deltaPsi = JetHelperTools::DeltaPsi (jetTruthPhi, vN_fcal);

      if ( (deltaPsi > 0) && (deltaPsi < TMath::Pi() / 4))         {
        h2_jet4_JES_truth_inPlane[cent_bin]->Fill (jetTruthPt, JES);
        h2_jet4_JES_truth_dpT_inPlane[cent_bin]->Fill (jetTruthPt, (jetTruthPt - pT_min));
        }
      else
        if ( (deltaPsi > TMath::Pi() / 4) && (deltaPsi < TMath::Pi() / 2))  {
          h2_jet4_JES_truth_outPlane[cent_bin]->Fill (jetTruthPt, JES);
          h2_jet4_JES_truth_dpT_outPlane[cent_bin]->Fill (jetTruthPt, (jetTruthPt - pT_min));

          }

      }

    }

  //END of loop ever truth jets















  return EL::StatusCode::SUCCESS;
  }

EL::StatusCode InclusiveJetsEventLoop :: postExecute () {
  return EL::StatusCode::SUCCESS;
  }

EL::StatusCode InclusiveJetsEventLoop :: finalize () {
  //xAOD::TEvent* event = wk()->xaodEvent();

  // cleaning up trigger tools
  if (_data_switch == 0) {
    if (m_trigConfigTool) {
      delete m_trigConfigTool;
      m_trigConfigTool = 0;
      }

    if (m_trigDecisionTool) {
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
  if (m_jetCleaning) {
    delete m_jetCleaning;
    m_jetCleaning = 0;
    }

  /*
    // for some reason this does not work
    for (int cent_bin = 0; cent_bin < 7; cent_bin++) {
      h2_jet4FJR_EtaPhiAvgNconst[cent_bin]->Divide (h2_jet4FJR_EtaPhi[cent_bin]);
      //h2_jet4FJR_EtaPhiAvgSubtrPt[cent_bin]->Divide( h2_jet4FJR_EtaPhi[cent_bin] );
      //h2_jet4FJR_EtaPhiDevSubtrPt[cent_bin]->Divide( h2_jet4FJR_EtaPhi[cent_bin] );

      for (int i = 1; i <= h2_jet4FJR_EtaPhiAvgSubtrPt[cent_bin]->GetXaxis()->GetNbins(); i++) {
        for (int j = 1; j <= h2_jet4FJR_EtaPhiAvgSubtrPt[cent_bin]->GetYaxis()->GetNbins(); j++) {

          Float_t avg = h2_jet4FJR_EtaPhiAvgSubtrPt[cent_bin]->GetBinContent (i, j);
          Float_t dev = h2_jet4FJR_EtaPhiDevSubtrPt[cent_bin]->GetBinContent (i, j);
          Float_t num = h2_jet4FJR_EtaPhi[cent_bin]->GetBinContent (i, j);
          dev = (num > 0) ? dev / num : dev;
          avg = (num > 0) ? avg / num : avg;

          if ( (dev - avg * avg) < 0)
            printf ("!!. stdev negativ \n");

          Float_t content = sqrt (fabs (dev - avg * avg));
          h2_jet4FJR_EtaPhiDevSubtrPt[cent_bin]->SetBinContent (i, j, content);
          h2_jet4FJR_EtaPhiAvgSubtrPt[cent_bin]->SetBinContent (i, j, avg);

          if (cent_bin == 0)
            printf ("%i, %i | %.2f, %.2f | %f  -- %f , %f \n", i, j, h2_jet4FJR_EtaPhiAvgSubtrPt[cent_bin]->GetXaxis()->GetBinCenter (i), h2_jet4FJR_EtaPhiAvgSubtrPt[cent_bin]->GetYaxis()->GetBinCenter (j), content, dev, avg);
          }
        }
      }
  */
  return EL::StatusCode::SUCCESS;
  }

EL::StatusCode InclusiveJetsEventLoop :: histFinalize () {


  cout << "Events = " << m_eventCounter << endl;
  return EL::StatusCode::SUCCESS;
  }






/**
 * @brief in vector @param v sets -1 to position @param pos
 *        if size of the vector @param v is smaller then position to be set to -1 write a message and return
 *
 * @param pos ...
 * @param v ...
 * @return void
 */
void InclusiveJetsEventLoop::setZero (int pos, vector< float >* v) {


  if (v->size() < (unsigned int) pos)  {
    cout << "!!! v->size() < pos" << endl;
    return;
    }
  else {

    v->at (pos) = -1;


    }




  }
