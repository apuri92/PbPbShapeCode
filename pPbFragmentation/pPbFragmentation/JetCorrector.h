#ifndef JetCorrector_H
#define JetCorrector_H

#include "TFile.h"
#include "TH2F.h"
#include "TF1.h"
#include "pPbFragmentation/JetHelperTools.h"
#include <fstream>
#include <TSystem.h>
#include <TAxis.h>

using namespace std;
using namespace JetHelperTools;

class JetCorrector
{
  private:
    // meaning of these quantities is explained in the constructor where the default values are assigned
    TFile * weight_file;
    
    TH3D *event_weight_histo;
	TH1D *cent_weight_histo;
	TH1D* FCal_HP_v_MB_weights_histo;
	TH1D* FCal_HP_v_MBOV_weights_histo;
	TH1D* FCal_HP_v_MCOV_weights_histo;
	TH1D* FCal_MCOV_v_MB_weights_histo;
	TH1D* FCal_MCOV_v_MBOV_weights_histo;
	vector<double> range_lo;
	vector<double> range_hi;
	vector<int> centiles;
	
	TFile * _f_JER;
    TH1F* _h_JER[8][8];
    Int_t _nEta_JER;
    Int_t _nCent_JER;
    Int_t _nCent_reweighting;
    
    TFile * _f_reweighting;
    TAxis * jet_pt_binning;
    TF1 *jet_spectra_weight[8][8];
	TH1D *FF_weight[8][8][20];
	TH1D *CHPS_weight[8][8][20];
	TH1D *FF_weight_fine[8][8][20];
	TH1D *CHPS_weight_fine[8][8][20];
	int etabin;
	TFile * HP_v_MB_FCAl_weight_file;

	TFile * eta_factors;
	TH1D * eta_weight[6][20];

  public:
    int nJetYBins;
    float JERcut;
    std::string _weight_file; 
	std::string _centrality_weight_file;
	std:: string _HP_v_MB_FCAl_weight_file;
	std:: string _eta_weight_file;

	float min_jet_pt;
	float max_jet_pt;
	bool m_isMB;
	bool is_pp;

    JetCorrector() 
     {

         m_isMB = false;
         nJetYBins = 5;			// Default number of y bins
         JERcut = 999; //Cut in sigma of JER
		 is_pp = false; //implemented for running over pp MC where no reweighting is required
         //Track-to-jet balance
		 _nEta_JER = 8;			// number of eta bins 
		 _nCent_JER = 8;			// number of centrality bins      
		 TString xfn = gSystem->GetFromPipe("echo $ROOTCOREBIN");  
		 _f_JER = new TFile(xfn + "/../pPbFragmentation/data/jet_perf_histos.root","read");
		 //Shifted by unity in Laura's file 
		 for(int i=0;i<_nEta_JER;i++){
			for(int j=0;j<_nCent_JER;j++){
				_h_JER[i][j]=(TH1F*)_f_JER->Get(Form("h1_JER_pT_y%i_c%i",i,j+1));
			}				
		 }
         
         //Reweighting
         _nCent_reweighting = 7;			// number of centrality bins      
		 _f_reweighting = new TFile(xfn + "/../pPbFragmentation/data/spectra_weights_PbPb.root","read");
		 jet_pt_binning = (TAxis*) ((TH1F*)_f_reweighting->Get("h_reco_jet_spectrum_4_cent0_system_0_PbPb"))->GetXaxis(); 
		 		 
		 //Now inclusive in eta 
		 //for(int etabin=0;etabin<nJetYBins;etabin++){
		 int etabin=4;
			for(int j=0;j<_nCent_reweighting;j++){
				//Spectra weights
				jet_spectra_weight[etabin][j]=(TF1*)_f_reweighting->Get(Form("jet_spectra_weight_%i_cent%i_PbPb",etabin,j));
				for(int k=1;k<=jet_pt_binning->GetNbins();k++){
					//FF weights
					FF_weight[etabin][j][k] = (TH1D*)_f_reweighting->Get(Form("ff_weight_%i_cent%i_system_0_PbPb_jet_pt%i",etabin,j,k));
					CHPS_weight[etabin][j][k] = (TH1D*)_f_reweighting->Get(Form("CHPS_weight_%i_cent%i_system_0_PbPb_jet_pt%i",etabin,j,k));
					FF_weight_fine[etabin][j][k] = (TH1D*)_f_reweighting->Get(Form("ff_weight_fine_%i_cent%i_system_0_PbPb_jet_pt%i",etabin,j,k));
					CHPS_weight_fine[etabin][j][k] = (TH1D*)_f_reweighting->Get(Form("CHPS_weight_fine_%i_cent%i_system_0_PbPb_jet_pt%i",etabin,j,k));
				}
			}				
		 //}
         
         //Event weights
         _HP_v_MB_FCAl_weight_file="fcal_weights.root";
         _weight_file="Powheg.reweight.root";
	     _centrality_weight_file="MB_FCal_Normalization.txt";
		 _eta_weight_file="eta_factors.root";

		TString HP_v_MB_FCAl_weight = xfn + "/../pPbFragmentation/data/"+ _HP_v_MB_FCAl_weight_file;
		TString event_weight_file = xfn + "/../pPbFragmentation/data/"+ _weight_file;
		TString fcal_weight_file = xfn +"/../pPbFragmentation/data/"+ _centrality_weight_file;
		TString eta_weight_file = xfn +"/../pPbFragmentation/data/"+ _eta_weight_file;


		HP_v_MB_FCAl_weight_file = new TFile(HP_v_MB_FCAl_weight.Data());
		weight_file = new TFile(event_weight_file.Data());
		event_weight_histo = (TH3D*)weight_file->Get("h3_pT_y_phi_rw");
		cent_weight_histo = (TH1D*)weight_file->Get("h1_cent_rw");
		FCal_HP_v_MB_weights_histo = (TH1D*)HP_v_MB_FCAl_weight_file->Get("weight_MB_to_HP");
		FCal_HP_v_MBOV_weights_histo = (TH1D*)HP_v_MB_FCAl_weight_file->Get("weight_MBov_to_HP");
		FCal_HP_v_MCOV_weights_histo = (TH1D*)HP_v_MB_FCAl_weight_file->Get("weight_MC_to_HP");
		FCal_MCOV_v_MB_weights_histo = (TH1D*)HP_v_MB_FCAl_weight_file->Get("weight_MB_to_MC");
		FCal_MCOV_v_MBOV_weights_histo = (TH1D*)HP_v_MB_FCAl_weight_file->Get("weight_MBov_to_MC");

		std::ifstream ifs (fcal_weight_file.Data(), std::ifstream::in);
		if(!ifs)
		{
			cout << endl << "Failed to open file " << fcal_weight_file.Data();
		}

		int tmp1; double tmp2, tmp3;
		while (ifs >> tmp1 >> tmp2 >> tmp3)
		{
			range_lo.push_back(tmp2);
			range_hi.push_back(tmp3);
			centiles.push_back(tmp1);
			if(!ifs.good())break;
		}    

		 eta_factors = new TFile(eta_weight_file.Data());
		 int n_pt_bins = ((TAxis*)eta_factors->Get("jetpT_binning"))->GetNbins();
		 for (int i_cent = 0; i_cent < 6; i_cent++)
		 {
			 for (int i_jet = 0; i_jet < n_pt_bins; i_jet++)
			 {
				 std::string name = Form("eta_factor_cent%i_jet%i", i_cent, i_jet);
				 eta_weight[i_cent][i_jet] = (TH1D*)eta_factors->Get(name.c_str());
			 }
		 }


     }// end of constructor
		
	float GetFCalWeight(float FCalEt);
	float GetFCalWeight(float FCalEt, int sample=1);
	//float GetFCalHPWeight(float FCalEt);
	float GetJetWeight(double pt, double eta, double phi);
	int GetJetYBin(float y);
	int GetJetpTBin(float pT, TAxis* pt_bins);
	Int_t GetJetEtaBin(float eta);
	bool MCJetJERClean(float truth_jet_pt,float reco_jet_pt, float truth_jet_eta, int cent);
	float GetJER(float truth_jet_pt, float truth_jet_eta, int cent);
	float GetJetReweightingFactor(double pt, double eta, int cent);
	float GetFFReweightingFactor(double z, double jet_pt, double jet_eta, int cent, bool isFine);
	float GetCHPSReweightingFactor(double pt, double jet_pt, double jet_eta, int cent, bool isFine);
	float GetEtaReweightingFactor(double jet_pt, double jet_eta, int cent);
    ~JetCorrector() {}
};


#endif
