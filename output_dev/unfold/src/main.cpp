
#include <iostream>
using std::cout;
using std::endl;
#include <vector>

#include "TRandom.h"
#include <TEnv.h>
#include "/Applications/root/RooUnfold/src/RooUnfoldResponse.h"
#include "/Applications/root/RooUnfold/src/RooUnfold.h"
#include "/Applications/root/RooUnfold/src/RooUnfoldBayes.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TSystem.h"
#include "TAxis.h"
#include "TError.h"
#include "TCanvas.h"
#include "TLine.h"


int main()
{
	gErrorIgnoreLevel = 3001;
	cout << "Unfolding..." << endl;

	//	##############	Reading config	##############"
	TEnv *m_config = new TEnv();
	m_config->ReadFile("ff_config.cfg", EEnvLevel(1));
	m_config->Print();

	std::string dataset_type = "PbPb"; dataset_type = m_config->GetValue("dataset_type", dataset_type.c_str());
	std::string tracking_cut = "ppTight"; tracking_cut = m_config->GetValue("tracking_cut", tracking_cut.c_str());
	int centrality_scheme = 31; centrality_scheme = m_config->GetValue("centrality_scheme", centrality_scheme);
	int isMC = 1; isMC = m_config->GetValue("isMC", isMC);

	std::string did = "data";
	if (isMC) did = "MC_JZ_comb";

	//	##############	Config done	##############"




	int is_data = 0;
	std::string ChPS_raw_type = "ChPS_raw_0";
	int n_unfold = 4;
	bool do_unfolding = 1;
	bool do_UE_subtr = 1;
	bool do_BbB = 1;
    std::string dataset_type = "PbPb";
    
	TEnv *m_config = new TEnv();
	m_config->ReadFile("config.cfg", EEnvLevel(1));
	m_config->Print();

	is_data = m_config->GetValue("is_data", is_data);
	ChPS_raw_type = m_config->GetValue("ChPS_raw_type", ChPS_raw_type.c_str());
	n_unfold = m_config->GetValue("n_unfold_iterations", n_unfold);
	do_UE_subtr = m_config->GetValue("do_UE_subtr", do_UE_subtr);
	do_unfolding = m_config->GetValue("do_unfolding", do_unfolding);
	do_BbB = m_config->GetValue("do_BbB", do_BbB);
    dataset_type = m_config->GetValue("dataset_type", dataset_type.c_str());


	TFile *f_mc = new TFile(Form("../FF_MC_JZ_comb_out_histo_%s_5p02_r001.root", dataset_type.c_str()));
    cout << f_mc->GetName() << endl;
	TFile *f_data; //using MC for closure test, data for actual data
	if (is_data) f_data = new TFile(Form("../FF_data_out_histo_%s_5p02_r001.root", dataset_type.c_str()));
	else f_data = new TFile(Form("../FF_MC_JZ_comb_out_histo_%s_5p02_r001.root", dataset_type.c_str()));

	std::string dataset = "data";
	if (!is_data) dataset = "mc";
	TFile *dr_factors = new TFile("dR_factors.root"); //getting bin by bin correction factors

	std::string evol = "noSubtr_noUnf_noBbB";
	if (do_UE_subtr) evol = "Subtr_noUnf_noBbB";
	if (do_unfolding) evol = "Subtr_Unf_noBbB";
	if (do_BbB) evol = "Subtr_Unf_BbB";

	TFile *f_output = new TFile(Form("unfolded_%s_%s_%s_%s.root",ChPS_raw_type.c_str(),dataset.c_str(), dataset_type.c_str(), evol.c_str() ),"recreate");
	std::string name;

	int N_CENT = 6;
	int N_Y = 5;


	TAxis* dR_binning = (TAxis*)((TH3*)f_mc->Get("h_dR_change_jetpt0_cent0"))->GetXaxis();
	TAxis* jetpT_binning = (TAxis*)((TH3*)f_mc->Get("ChPS_raw_0_dR0_cent0"))->GetYaxis();
	TAxis* trkpT_binning = (TAxis*)((TH3*)f_mc->Get("ChPS_raw_0_dR0_cent0"))->GetXaxis();

	int N_dR = dR_binning->GetNbins();
	int N_jetpt = jetpT_binning->GetNbins();
	int N_trkpt = trkpT_binning->GetNbins();


	f_output->cd();
	dR_binning->Write("dR_binning");
	jetpT_binning->Write("jetpT_binning");
	trkpT_binning->Write("trkpT_binning");

	TCanvas *c = new TCanvas("c","c",800,600);
	c->Divide(4,2);
	TLine *line = new TLine();
	line->SetLineStyle(3);
	line->SetLineColor(kBlack);

//	c->Print("inclusive_jet_spect_closure.pdf(","Title: Start");
	
	for (int i_cent = 0; i_cent < N_CENT; i_cent++)
	{
        
		TH1* h_unfolded_jet_spect_y4;
		TH1* h_truth_jet_spect_y4;
		TH1* h_raw_jet_spect_y4;

		for (int i_y = 0; i_y < N_Y; i_y++)
		{
			std::string jet_spect_hist;
			if (ChPS_raw_type == "ChPS_raw_0") jet_spect_hist = Form("h_reco_jet_spectrum_y%i_cent%i",i_y, i_cent);
			if (ChPS_raw_type == "ChPS_raw_rr") jet_spect_hist = Form("h_reco_jet_spectrum_matched_y%i_cent%i",i_y, i_cent);
			if (ChPS_raw_type == "ChPS_raw_tr") jet_spect_hist = Form("h_reco_jet_spectrum_matched_y%i_cent%i",i_y, i_cent);
			if (ChPS_raw_type == "ChPS_raw_rt") jet_spect_hist = Form("h_true_jet_spectrum_matched_y%i_cent%i",i_y, i_cent);
			if (ChPS_raw_type == "ChPS_raw_tt") jet_spect_hist = Form("h_true_jet_spectrum_matched_y%i_cent%i",i_y, i_cent);
			if (ChPS_raw_type == "ChPS_raw_tt_mod") jet_spect_hist = Form("h_true_jet_spectrum_matched_y%i_cent%i",i_y, i_cent);

            cout << jet_spect_hist << endl;
			TH1* h_raw = (TH1*)((TH1*)f_data->Get(jet_spect_hist.c_str()))->Clone(jet_spect_hist.c_str());
			h_raw->SetName(jet_spect_hist.c_str());
			TH1* h_truth = (TH1*)((TH1*)f_mc->Get(Form("h_true_jet_spectrum_y%i_cent%i",i_y, i_cent)))->Clone(Form("true_jet_spectrum_y%i_cent%i",i_y, i_cent));
			TH1* h_unfolded = (TH1*)h_truth->Clone(Form("c_unf_y%i_c%i", i_y, i_cent));
            h_unfolded->Reset();
			h_raw->Sumw2();
			h_truth->Sumw2();
//            h_raw->Print("all");
//            h_truth->Print("all");
            
//            if ((dataset_type == "pp" && i_cent == 5) ||
//                dataset_type == "PbPb")
            {
                //unfold only if using proper reco ChPS (not if using *reco matched* truth)
                if ((ChPS_raw_type == "ChPS_raw_0" || ChPS_raw_type == "ChPS_raw_rr") && do_unfolding)
                {
                    RooUnfoldResponse* r_response = (RooUnfoldResponse*)f_mc->Get(Form("h_reco_jet_spectrum_y%i_cent%i_h_true_jet_spectrum_y%i_cent%i",i_y, i_cent, i_y, i_cent));
                    r_response->Print("all");
                    RooUnfoldBayes unfold(r_response,h_raw,4);
                    unfold.SetVerbose(0);
                    h_unfolded = (TH1D*)unfold.Hreco();
                    h_unfolded->Print("all");
                    name = Form("Unfolded_jetSpect_c%i_y%i",i_cent, i_y);
                    h_unfolded->SetName(name.c_str());
                }
                else
                {
                    name = Form("nonUnfolded_jetSpect_c%i_y%i",i_cent, i_y);
                    h_unfolded = (TH1*)h_raw->Clone(name.c_str());
                    h_unfolded->SetName(name.c_str());
                }
                
                h_unfolded->Sumw2();
                
                if (i_y == 4)
                {
                    TH1* h_ratio = (TH1*)h_unfolded->Clone(Form("ratio_c%i_y%i", i_cent, i_y));
                    h_ratio->Divide(h_truth);
                    
                    TH1* h_raw_ratio = (TH1*)h_raw->Clone(Form("raw_c%i_y%i", i_cent, i_y));
                    h_raw_ratio->Divide(h_truth);
                    
                    c->cd(i_cent+1);
                    h_ratio->GetYaxis()->SetRangeUser(0.,2.);
                    h_ratio->GetXaxis()->SetRangeUser(80,1000);
                    
                    h_ratio->SetMarkerStyle(20);
                    h_raw_ratio->SetMarkerStyle(20);
                    
                    h_ratio->SetMarkerColor(kBlack);
                    h_raw_ratio->SetMarkerColor(kRed);
                    
                    h_ratio->Draw();
                    h_raw_ratio->Draw("same");
                    
                    line->DrawLine(80,1,1000,1);
                    
                    gPad->SetLogx();
                    h_unfolded->Print("all");
                    c->Print("inclusive_jet_spect_closure.pdf");//,Form("Title: ratio_c%i_y%i", i_cent, i_y));
                }
                
                //Use n_jets from spectra inclusive in y
                if (i_y == 4)
                {
                    h_raw_jet_spect_y4 = (TH1*)h_raw->Clone(Form("h_n_raw_jet_spect_y%i_cent%i",i_y, i_cent));
                    h_unfolded_jet_spect_y4 = (TH1*)h_unfolded->Clone(Form("h_n_unf_jet_spect_y%i_cent%i",i_y, i_cent));
                    h_truth_jet_spect_y4 = (TH1*)h_truth->Clone(Form("h_n_tru_jet_spect_y%i_cent%i",i_y, i_cent));;
                }
                
            }
            f_output->cd();
            h_raw->Write(Form("h_raw_jet_spect_cent%i_y%i",i_cent, i_y));
            h_truth->Write(Form("h_true_jet_spect_cent%i_y%i",i_cent, i_y));
            h_unfolded->Write(Form("h_unfolded_jet_spect_cent%i_y%i",i_cent, i_y));

            cout << Form("h_unfolded_jet_spect_cent%i_y%i",i_cent, i_y) << endl;
        }

		TH2* h_raw_injet;
		TH2* h_UE_injet;


		for (int i_dR = 0; i_dR < N_dR; i_dR++)
		{
            cout << __LINE__ << " " << Form("%s_dR%i_cent%i",ChPS_raw_type.c_str(), i_dR, i_cent) << endl;

			TH2* h_raw = (TH2*)f_data->Get(Form("%s_dR%i_cent%i",ChPS_raw_type.c_str(), i_dR, i_cent)); //ChPS_raw_0 or ChPS_raw_rr or ChPS_raw_rt or ChPS_raw_tr or ChPS_raw_tt
            cout << __LINE__ << " " << Form("%s_dR%i_cent%i",ChPS_raw_type.c_str(), i_dR, i_cent) << endl;
            TH2* h_truth = (TH2*)f_mc->Get(Form("ChPS_truth_dR%i_cent%i",i_dR, i_cent));
			TH2* h_unfolded = (TH2*) h_raw->Clone(Form("nonUnfolded_ChPS_c%i_dR%i",i_cent, i_dR));
            h_unfolded->Reset();
			TH2* h_UE = (TH2*)f_mc->Get(Form("ff_UE_pT_dR%i_cent%i",i_dR, i_cent));
            cout << __LINE__ << " " << Form("%s_dR%i_cent%i",ChPS_raw_type.c_str(), i_dR, i_cent) << endl;
            
            
            h_raw->Sumw2();
            h_truth->Sumw2();
            cout << __LINE__ << " " << Form("%s_dR%i_cent%i",ChPS_raw_type.c_str(), i_dR, i_cent) << endl;
            
//            if ((dataset_type == "pp" && i_cent == 5) ||
//                dataset_type == "PbPb")
            {
                
                cout << Form("%s_dR%i_cent%i",ChPS_raw_type.c_str(), i_dR, i_cent) << endl;
                //UE needs to take into account different number of reco jets in MC and data
                TH1* h_reco_mc_jet_spect = (TH1*)((TH1*)f_mc->Get(Form("h_reco_jet_spectrum_y%i_cent%i", 4, i_cent)))->Clone(Form("norm_reco_mc_jet_y%i_cent%i", 4, i_cent));
                cout << __LINE__ << endl;
                h_reco_mc_jet_spect->SetName(Form("norm_reco_mc_jet_y4_cent%i", i_cent));
                TH1* h_reco_data_jet_spect = (TH1*)((TH1*)f_data->Get(Form("h_reco_jet_spectrum_y%i_cent%i", 4, i_cent)))->Clone(Form("norm_reco_data_jet_y%i_cent%i", 4, i_cent));
                h_reco_data_jet_spect->SetName(Form("norm_reco_data_jet_y4_cent%i", i_cent));
                
                h_reco_mc_jet_spect->Sumw2();
                h_reco_data_jet_spect->Sumw2();
                
                for (int i_jet_bin = 1; i_jet_bin <= N_jetpt; i_jet_bin++)
                {
                    double n_jets_mc = h_reco_mc_jet_spect->GetBinContent(i_jet_bin);
                    double n_jets_data = h_reco_data_jet_spect->GetBinContent(i_jet_bin);
                    
                    if (n_jets_mc == 0) continue;
                    
                    for (int i_trk_bin = 1; i_trk_bin <= N_jetpt; i_trk_bin++)
                    {
                        double updated_UE = h_UE->GetBinContent(i_trk_bin, i_jet_bin) * n_jets_data / n_jets_mc;
                        double updated_UE_err = h_UE->GetBinError(i_trk_bin, i_jet_bin) * n_jets_data / n_jets_mc;
                        
                        h_UE->SetBinContent(i_trk_bin, i_jet_bin, updated_UE);
                        h_UE->SetBinError(i_trk_bin, i_jet_bin, updated_UE_err); //scaled errors
                    }
                }
                h_UE->Sumw2();
                
                //Do UE subtraction and unfolding and bin by bin correction if using proper reco ChPS
                if (ChPS_raw_type == "ChPS_raw_0" && do_UE_subtr) h_raw->Add(h_UE,-1); //End UE Subtraction
                
                //Unfolding
                if ((ChPS_raw_type == "ChPS_raw_0" || ChPS_raw_type == "ChPS_raw_rr") && do_unfolding)
                {
                    RooUnfoldResponse* r_response = (RooUnfoldResponse*)f_mc->Get(Form("ChPS_raw_0_dR%i_cent%i_ChPS_truth_dR%i_cent%i", i_dR, i_cent, i_dR, i_cent));
                    RooUnfoldBayes unfold(r_response, h_raw, n_unfold);
                    unfold.SetVerbose(0);
                    
                    h_unfolded = (TH2D*)unfold.Hreco(); //errors handled internally
                    name = Form("Unfolded_ChPS_c%i_dR%i",i_cent, i_dR);
                    h_unfolded->SetName(name.c_str());
                    
                } //End unfolding
                else
                {
                    name = Form("nonUnfolded_ChPS_c%i_dR%i",i_cent, i_dR);
                    h_unfolded = (TH2*)h_raw->Clone(name.c_str());
                    h_unfolded->SetName(name.c_str());
                }
                
                h_unfolded->Sumw2();
                
                
                
                //Bin by bin correction factors
                if ((ChPS_raw_type == "ChPS_raw_0" || ChPS_raw_type == "ChPS_raw_rr") && do_BbB)
                {
                    for (int i_jet_bin = 0; i_jet_bin < N_jetpt; i_jet_bin++)
                    {
                        name = Form("h_bin_by_bin_cent%i_jetpt%i_dR%i",i_cent, i_jet_bin, i_dR);
                        TH1* h_bin_by_bin = (TH1*)dr_factors->Get(name.c_str());
                        h_bin_by_bin->SetName(name.c_str());
                        
                        cout << name << endl;
                        for (int i_trk_bin = 0; i_trk_bin < N_trkpt; i_trk_bin++)
                        {
                            double original = h_unfolded->GetBinContent(i_trk_bin+1, i_jet_bin+1);
                            double bin_by_bin = h_bin_by_bin->GetBinContent(i_trk_bin+1);
                            double original_err = h_unfolded->GetBinError(i_trk_bin+1, i_jet_bin+1);
                            double bin_by_bin_err = h_bin_by_bin->GetBinError(i_trk_bin+1);
                            
                            if (original == 0 || bin_by_bin == 0) continue;
                            
                            double corrected = original * bin_by_bin;
                            double corrected_err = fabs(corrected) * sqrt(pow(original_err/original,2) + pow(bin_by_bin_err/bin_by_bin,2));
                            
                            h_unfolded->SetBinContent(i_trk_bin+1, i_jet_bin+1, corrected);
                            h_unfolded->SetBinError(i_trk_bin+1, i_jet_bin+1, corrected_err); //errors propagated
                            
                        } //end trk bin loop
                    } //end jet bin loop
                }//End bin by bin correction (done after unfolding)
                
                
                //doing per jet normalization
                for (int i_jet_bin = 0; i_jet_bin < N_jetpt; i_jet_bin++)
                {
                    double n_jets_raw = h_raw_jet_spect_y4->GetBinContent(i_jet_bin+1);
                    double n_jets_unf = h_unfolded_jet_spect_y4->GetBinContent(i_jet_bin+1);
                    double n_jets_tru = h_truth_jet_spect_y4->GetBinContent(i_jet_bin+1);
                    
                    if (n_jets_raw == 0 || n_jets_unf == 0 || n_jets_tru == 0) continue;
                    
                    for (int i_trk_bin = 0; i_trk_bin < N_trkpt; i_trk_bin++)
                    {
                        double updated_raw = h_raw->GetBinContent(i_trk_bin+1,i_jet_bin+1) / n_jets_raw;
                        double updated_reco = h_unfolded->GetBinContent(i_trk_bin+1,i_jet_bin+1) / n_jets_unf;
                        double updated_truth = h_truth->GetBinContent(i_trk_bin+1,i_jet_bin+1) / n_jets_tru;
                        double updated_UE = h_UE->GetBinContent(i_trk_bin+1,i_jet_bin+1) / n_jets_raw;
                        
                        double updated_raw_err =  h_raw->GetBinError(i_trk_bin+1,i_jet_bin+1) / n_jets_raw;
                        double updated_reco_err = h_unfolded->GetBinError(i_trk_bin+1,i_jet_bin+1) / n_jets_unf;
                        double updated_truth_err = h_truth->GetBinError(i_trk_bin+1,i_jet_bin+1) / n_jets_tru;
                        double updated_UE_err = h_UE->GetBinError(i_trk_bin+1,i_jet_bin+1) / n_jets_raw;
                        
                        h_raw->SetBinContent(i_trk_bin+1, i_jet_bin+1, updated_raw);
                        h_unfolded->SetBinContent(i_trk_bin+1, i_jet_bin+1, updated_reco);
                        h_truth->SetBinContent(i_trk_bin+1, i_jet_bin+1, updated_truth);
                        h_UE->SetBinContent(i_trk_bin+1, i_jet_bin+1, updated_UE);
                        
                        h_raw->SetBinError(i_trk_bin+1, i_jet_bin+1, updated_raw_err); //scaled errors
                        h_unfolded->SetBinError(i_trk_bin+1, i_jet_bin+1, updated_reco_err); //scaled errors
                        h_truth->SetBinError(i_trk_bin+1, i_jet_bin+1, updated_truth_err); //scaled errors
                        h_UE->SetBinError(i_trk_bin+1, i_jet_bin+1, updated_UE_err); //scaled errors
                    }
                }
                
                if (i_dR == 0)
                {
                    h_UE_injet = (TH2*)h_UE->Clone(Form("UE_inJet_cent%i", i_cent));;
                    h_raw_injet = (TH2*)h_raw->Clone(Form("UE_inJet_cent%i", i_cent));;
                }
                
                if (i_dR > 0 && i_dR < 7)
                {
                    h_UE_injet->Add(h_UE);
                    h_raw_injet->Add(h_raw);
                }
            }
            cout << __LINE__ << " " << Form("%s_dR%i_cent%i",ChPS_raw_type.c_str(), i_dR, i_cent) << endl;

			f_output->cd();

			h_raw->Write(Form("h_%s_cent%i_dR%i",ChPS_raw_type.c_str(),i_cent, i_dR));
            cout << __LINE__ << " " << Form("%s_dR%i_cent%i",ChPS_raw_type.c_str(), i_dR, i_cent) << endl;
			h_truth->Write(Form("h_true_ChPS_cent%i_dR%i",i_cent, i_dR));
            cout << __LINE__ << " " << Form("%s_dR%i_cent%i",ChPS_raw_type.c_str(), i_dR, i_cent) << endl;
			h_unfolded->Write(Form("h_unfolded_%s_cent%i_dR%i",ChPS_raw_type.c_str(), i_cent, i_dR));
            cout << __LINE__ << " " << Form("%s_dR%i_cent%i",ChPS_raw_type.c_str(), i_dR, i_cent) << endl;
			h_UE->Write(Form("h_UE_cent%i_dR%i", i_cent, i_dR));
            cout << __LINE__ << " " << Form("%s_dR%i_cent%i",ChPS_raw_type.c_str(), i_dR, i_cent) << endl;

		}

		f_output->cd();
		h_raw_injet->Write(Form("h_%s_injet_cent%i",ChPS_raw_type.c_str(),i_cent));
		h_UE_injet->Write(Form("h_UE_injet_cent%i", i_cent));
	}

//	((TH3*)f_mc->Get("h_dR_change_jetpt0_cent0"))->Write("h_dR_change_jetpt0_cent0");
//	f_output->Close();
	


	cout << "Done unfolding" << endl;
	return 0;
}
