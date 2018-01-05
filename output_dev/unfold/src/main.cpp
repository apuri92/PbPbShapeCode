#include <iostream>
using std::cout;
using std::endl;
#include <vector>

#include "TRandom.h"
#include <TEnv.h>
#include "/Applications/root_v6.10.08/RooUnfold/src/RooUnfoldResponse.h"
#include "/Applications/root_v6.10.08/RooUnfold/src/RooUnfold.h"
#include "/Applications/root_v6.10.08/RooUnfold/src/RooUnfoldBayes.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TSystem.h"
#include "TAxis.h"
#include "TError.h"
#include "TCanvas.h"
#include "TLine.h"
#include "../../functions/global_variables.h"

int main()
{
	gErrorIgnoreLevel = 3001;
	cout << "Unfolding..." << endl;
//	gROOT->LoadMacro("AtlasStyle.C");
//	SetAtlasStyle();
	//	##############	Reading config	##############"
	TEnv *m_config = new TEnv();
	m_config->ReadFile("ff_config.cfg", EEnvLevel(1));
	m_config->Print();

	std::string dataset_type = "PbPb"; dataset_type = m_config->GetValue("dataset_type", dataset_type.c_str());
	std::string tracking_cut = "ppTight"; tracking_cut = m_config->GetValue("tracking_cut", tracking_cut.c_str());
	int centrality_scheme = 31; centrality_scheme = m_config->GetValue("centrality_scheme", centrality_scheme);
	int isMC = 1; isMC = m_config->GetValue("isMC", isMC);
//	int n_unfold = 4; n_unfold = m_config->GetValue("n_unfold", n_unfold);
	int n_unfold = 4;
	std::string did = "data";
	if (isMC) did = "MC_JZ_comb";

	//	##############	Config done	##############"

	TFile *f_mc = new TFile(Form("../raw_results/FF_MC_JZ_comb_out_histo_%s_5p02_r001.root", dataset_type.c_str()));
	TFile *f_data = new TFile(Form("../raw_results/FF_%s_out_histo_%s_5p02_r001.root", did.c_str(), dataset_type.c_str()));
	TFile *dr_factors = new TFile(Form("posCorr_factors_%s.root", dataset_type.c_str()));
	TFile *f_output = new TFile(Form("unfolded_%s_%s.root",did.c_str(), dataset_type.c_str()),"recreate");
	std::string name;

//	int N_CENT = 6;
	//remove when rerun is complete
	int N_Y = 5;
	n_cent_cuts = 6;

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

	TCanvas *c = new TCanvas("c","c",1200,600);
	c->Divide(4,2);
	TLine *line = new TLine();
	line->SetLineStyle(3);
	line->SetLineColor(kBlack);

//	c->Print("inclusive_jet_spect_closure.pdf(","Title: Start");

	for (int i_cent = 0; i_cent < n_cent_cuts; i_cent++)
	{
		cout << Form("Done cent%i", i_cent) << endl;
		//ChPS: raw_0, rr(truth matched), truth, raw_unfolded, raw_rr_unfo
		//jet spect: reco, reco_matched, truth_matched, truth

		TH1 *h_reco_jet, *h_reco_jet_matched, *h_truth_jet, *h_truth_jet_matched, *h_reco_unfolded, *h_reco_matched_unfolded;


		for (int i_y = 0; i_y < N_Y; i_y++)
		{

			name = Form("h_reco_jet_spectrum_y%i_cent%i",i_y, i_cent);
			TH1* h_reco_jet_y_c = (TH1*)((TH1*)f_data->Get(name.c_str()))->Clone(Form("reco_jet_y%i_c%i",i_y, i_cent));
			h_reco_jet_y_c->Sumw2();

			name = Form("h_reco_jet_spectrum_matched_y%i_cent%i",i_y, i_cent);
			TH1* h_reco_jet_matched_y_c = (TH1*)((TH1*)f_mc->Get(name.c_str()))->Clone(Form("reco_jet_matched_y%i_c%i",i_y, i_cent));
			h_reco_jet_matched_y_c->Sumw2();

			name = Form("h_true_jet_spectrum_y%i_cent%i",i_y, i_cent);
			TH1* h_true_jet_y_c = (TH1*)((TH1*)f_mc->Get(name.c_str()))->Clone(Form("true_jet_y%i_c%i",i_y, i_cent));
			h_true_jet_y_c->Sumw2();

			name = Form("h_true_jet_spectrum_matched_y%i_cent%i",i_y, i_cent);
			TH1* h_true_jet_matched_y_c = (TH1*)((TH1*)f_mc->Get(name.c_str()))->Clone(Form("true_jet_matched_y%i_c%i",i_y, i_cent));
			h_true_jet_matched_y_c->Sumw2();

			TH1* h_unfolded_jet_y_c = (TH1*)h_true_jet_y_c->Clone(Form("h_raw_unfolded_y%i_c%i", i_y, i_cent));
            h_unfolded_jet_y_c->Reset();

			TH1* h_matched_unfolded_jet_y_c = (TH1*)h_true_jet_matched_y_c->Clone(Form("h_raw_matched_unfolded_y%i_c%i", i_y, i_cent));
			h_matched_unfolded_jet_y_c->Reset();

			//unfold only if using proper reco ChPS (not if using *reco matched* truth)
			name = Form("h_reco_jet_spectrum_y%i_cent%i_h_true_jet_spectrum_y%i_cent%i",i_y, i_cent, i_y, i_cent);
			RooUnfoldResponse* r_response = (RooUnfoldResponse*)f_mc->Get(name.c_str());

			RooUnfoldBayes unfold(r_response, h_reco_jet_y_c,4);
			unfold.SetVerbose(0);
			if (h_reco_jet_y_c->GetEntries() != 0)
			{
				h_unfolded_jet_y_c = (TH1D*)unfold.Hreco();
				name = Form("unfolded_y%i_c%i", i_y, i_cent);
				h_unfolded_jet_y_c->SetName(name.c_str());
				h_unfolded_jet_y_c->Sumw2();
			}

			RooUnfoldBayes matched_unfold(r_response, h_reco_jet_matched_y_c,4);
			matched_unfold.SetVerbose(0);
			if (h_reco_jet_matched_y_c->GetEntries() != 0)
			{
				h_matched_unfolded_jet_y_c = (TH1D*)matched_unfold.Hreco();
				name = Form("matched_unfolded_y%i_c%i", i_y, i_cent);
				h_matched_unfolded_jet_y_c->SetName(name.c_str());
				h_matched_unfolded_jet_y_c->Sumw2();
			}


			if (i_y == 4)
			{
				h_reco_jet = (TH1*)h_reco_jet_y_c->Clone(Form("h_inc_reco_jet_y%i_c%i",i_y, i_cent));
				h_reco_jet->Sumw2();

				h_reco_jet_matched = (TH1*)h_reco_jet_matched_y_c->Clone(Form("h_inc_reco_matched_jet_y%i_c%i",i_y, i_cent));
				h_reco_jet_matched->Sumw2();

				h_truth_jet = (TH1*)h_true_jet_y_c->Clone(Form("h_inc_truth_matched_jet_y%i_c%i",i_y, i_cent));
				h_truth_jet->Sumw2();

				h_truth_jet_matched = (TH1*)h_true_jet_matched_y_c->Clone(Form("h_inc_truth_jet_y%i_c%i",i_y, i_cent));
				h_truth_jet_matched->Sumw2();

				h_reco_unfolded = (TH1*)h_unfolded_jet_y_c->Clone(Form("h_inc_reco_unf_jet_y%i_c%i",i_y, i_cent));
				h_reco_unfolded->Sumw2();

				h_reco_matched_unfolded = (TH1*)h_matched_unfolded_jet_y_c->Clone(Form("h_inc_reco_matched_unf_jet_y%i_c%i",i_y, i_cent));
				h_reco_matched_unfolded->Sumw2();
//
//				//draw raw/unfolded/true jet spectra
//				TH1* h_ratio_unf_raw = (TH1*)h_reco_unfolded->Clone(Form("ratio_unf_raw_inc_y%i_c%i", i_y, i_cent));
//				h_ratio_unf_raw->Divide(h_reco_jet);
//
//				TH1* h_ratio_unf_truth = (TH1*)h_reco_unfolded->Clone(Form("ratio_unf_truth_inc_y%i_c%i", i_y, i_cent));
//				h_ratio_unf_truth->Divide(h_truth_jet);
//
//				TH1* h_ratio_raw_truth = (TH1*)h_reco_jet->Clone(Form("ratio_raw_truth_inc_y%i_c%i", i_y, i_cent));
//				h_ratio_raw_truth->Divide(h_truth_jet);
//
//				c->cd(i_cent+1);
//				h_ratio_raw_truth->GetYaxis()->SetRangeUser(0.,2.);
//				h_ratio_unf_raw->GetXaxis()->SetRangeUser(80,1000);
//				h_ratio_unf_truth->GetXaxis()->SetRangeUser(80,1000);
//				h_ratio_raw_truth->GetXaxis()->SetRangeUser(80,1000);
//
////				SetHStyle(h_ratio_unf_raw, 0);
////				SetHStyle(h_ratio_unf_truth, 1);
////				SetHStyle(h_ratio_raw_truth, 2);
//
//				h_ratio_unf_raw->Draw();
//				h_ratio_unf_truth->Draw("same");
//				h_ratio_raw_truth->Draw("same");
//
//				line->DrawLine(80,1,1000,1);
//
//				gPad->SetLogx();
//				if (i_cent == n_cent_cuts-1)
//				{
//					//					c->Print("inclusive_jet_spect_closure.pdf");//,Form("Title: ratio_c%i_y%i", i_cent, i_y));
//				}
//remove drawing from here. put in in post unfoldin code so you can use hstyle etc. will need to be done inclusive anyway

            }



            f_output->cd();
			name = Form("h_reco_jet_y%i_c%i", i_y, i_cent);
			h_reco_jet_y_c->SetTitle(name.c_str());
			h_reco_jet_y_c->Write(name.c_str());

			name = Form("h_true_jet_y%i_c%i", i_y, i_cent);
			h_true_jet_y_c->SetTitle(name.c_str());
			h_true_jet_y_c->Write(name.c_str());

			name = Form("h_unfolded_jet_y%i_c%i", i_y, i_cent);
			h_unfolded_jet_y_c->SetTitle(name.c_str());
			h_unfolded_jet_y_c->Write(name.c_str());

			name = Form("h_reco_jet_matched_y%i_c%i", i_y, i_cent);
			h_reco_jet_matched_y_c->SetTitle(name.c_str());
			h_reco_jet_matched_y_c->Write(name.c_str());

			name = Form("h_true_jet_matched_y%i_c%i", i_y, i_cent);
			h_true_jet_matched_y_c->SetTitle(name.c_str());
			h_true_jet_matched_y_c->Write(name.c_str());
		}


		TH2* h_raw_injet;
		TH2* h_raw_subtr_injet;
		TH2* h_raw_subtr_unf_injet;
		TH2* h_raw_subtr_unf_bbb_injet;

		TH2* h_raw_rr_injet;
		TH2* h_raw_rr_unf_injet;
		TH2* h_raw_rr_unf_bbb_injet;

		TH2* h_UE_injet;

		for (int i_dR = 0; i_dR < N_dR; i_dR++)
		{
			TH2* h_raw = (TH2*)f_data->Get(Form("ChPS_raw_0_dR%i_cent%i", i_dR, i_cent));
			h_raw->Sumw2();

			TH2* h_raw_rr = (TH2*)f_mc->Get(Form("ChPS_raw_rr_dR%i_cent%i", i_dR, i_cent));
			h_raw_rr->Sumw2();

			TH2* h_truth = (TH2*)f_mc->Get(Form("ChPS_truth_dR%i_cent%i", i_dR, i_cent));
			h_truth->Sumw2();

			TH2* h_UE = (TH2*)f_mc->Get(Form("ff_UE_pT_dR%i_cent%i",i_dR, i_cent));
			h_UE->Sumw2();

			//normalize UE by # of jets in data/MC
			TH1* h_reco_mc_jet_spect = (TH1*)((TH1*)f_mc->Get(Form("h_reco_jet_spectrum_y%i_cent%i", 4, i_cent)))->Clone(Form("norm_reco_mc_jet_y%i_cent%i", 4, i_cent));
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

			//UE subtraction
			TH2* h_raw_subtr = (TH2*)h_raw->Clone(Form("h_raw_subtr_dR%i_c%i", i_dR, i_cent));
			h_raw_subtr->Add(h_UE, -1);

			//Unfolding
			TH2* h_raw_subtr_unf = (TH2*)h_raw_subtr->Clone(Form("h_raw_subtr_unf_dR%i_c%i", i_dR, i_cent));
			TH2* h_raw_rr_unf = (TH2*)h_raw_rr->Clone(Form("h_raw_subtr_unf_dR%i_c%i", i_dR, i_cent));

			RooUnfoldResponse* r_response = (RooUnfoldResponse*)f_mc->Get(Form("ChPS_raw_0_dR%i_cent%i_ChPS_truth_dR%i_cent%i", i_dR, i_cent, i_dR, i_cent));

			RooUnfoldBayes unfold_raw(r_response, h_raw_subtr_unf, n_unfold);
			unfold_raw.SetVerbose(0);
			h_raw_subtr_unf = (TH2D*)unfold_raw.Hreco(); //errors handled internally
			name = Form("h_raw_subtr_unf_dR%i_c%i", i_dR, i_cent);
			h_raw_subtr_unf->SetName(name.c_str());
			h_raw_subtr_unf->Sumw2();


			RooUnfoldBayes unfold_raw_rr(r_response, h_raw_rr_unf, n_unfold);
			unfold_raw_rr.SetVerbose(0);
			h_raw_rr_unf = (TH2D*)unfold_raw_rr.Hreco(); //errors handled internally
			name = Form("h_raw_rr_unf_dR%i_c%i", i_dR, i_cent);
			h_raw_rr_unf->SetName(name.c_str());
			h_raw_rr_unf->Sumw2();

			//Bin by bin correction factors
			TH2* h_raw_subtr_unf_bbb = (TH2*)h_raw_subtr_unf->Clone(Form("h_raw_subtr_unf_bbb_dR%i_cent%i", i_dR, i_cent));
			TH2* h_raw_rr_unf_bbb = (TH2*)h_raw_rr_unf->Clone(Form("h_raw_subtr_unf_bbb_dR%i_cent%i", i_dR, i_cent));

			for (int i_jet_bin = 0; i_jet_bin < N_jetpt; i_jet_bin++)
			{
				name = Form("h_bin_by_bin_cent%i_jetpt%i_dR%i", i_cent, i_jet_bin, i_dR);
				TH1* h_bin_by_bin = (TH1*)dr_factors->Get(name.c_str());
				h_bin_by_bin->SetName(name.c_str());

				for (int i_trk_bin = 0; i_trk_bin < N_trkpt; i_trk_bin++)
				{
					double bin_by_bin = h_bin_by_bin->GetBinContent(i_trk_bin+1);
					double bin_by_bin_err = h_bin_by_bin->GetBinError(i_trk_bin+1);

					if (bin_by_bin == 0) continue;
					double original, original_err, corrected, corrected_err;

					original = h_raw_subtr_unf_bbb->GetBinContent(i_trk_bin+1, i_jet_bin+1);
					original_err = h_raw_subtr_unf_bbb->GetBinError(i_trk_bin+1, i_jet_bin+1);
					if (original != 0)
					{
						corrected = original * bin_by_bin;
						corrected_err = fabs(corrected) * sqrt(pow(original_err/original,2) + pow(bin_by_bin_err/bin_by_bin,2));
						h_raw_subtr_unf_bbb->SetBinContent(i_trk_bin+1, i_jet_bin+1, corrected);
						h_raw_subtr_unf_bbb->SetBinError(i_trk_bin+1, i_jet_bin+1, corrected_err); //errors propagated
					}

					original = h_raw_rr_unf_bbb->GetBinContent(i_trk_bin+1, i_jet_bin+1);
					original_err = h_raw_rr_unf_bbb->GetBinError(i_trk_bin+1, i_jet_bin+1);
					if (original != 0)
					{
						corrected = original * bin_by_bin;
						corrected_err = fabs(corrected) * sqrt(pow(original_err/original,2) + pow(bin_by_bin_err/bin_by_bin,2));
						h_raw_rr_unf_bbb->SetBinContent(i_trk_bin+1, i_jet_bin+1, corrected);
						h_raw_rr_unf_bbb->SetBinError(i_trk_bin+1, i_jet_bin+1, corrected_err); //errors propagated
					}

				} //end trk bin loop
			} //end jet bin loop

			//doing per jet normalization
			for (int i_jet_bin = 0; i_jet_bin < N_jetpt; i_jet_bin++)
			{
				double n_jets_raw = h_reco_jet->GetBinContent(i_jet_bin+1);
				double n_jets_unf = h_reco_unfolded->GetBinContent(i_jet_bin+1);

				double n_jets_raw_rr = h_reco_jet_matched->GetBinContent(i_jet_bin+1);
				double n_jets_unf_rr = h_reco_matched_unfolded->GetBinContent(i_jet_bin+1);

				double n_jets_tru = h_truth_jet->GetBinContent(i_jet_bin+1);

				
				if (n_jets_raw == 0 || n_jets_unf == 0 ||
					n_jets_raw_rr == 0 || n_jets_unf_rr == 0 ||
					n_jets_tru == 0) continue;

				for (int i_trk_bin = 0; i_trk_bin < N_trkpt; i_trk_bin++)
				{

					//ChPS_raw
					double updated_raw = h_raw->GetBinContent(i_trk_bin+1,i_jet_bin+1) / n_jets_raw;
					double updated_raw_subtr = h_raw_subtr->GetBinContent(i_trk_bin+1,i_jet_bin+1) / n_jets_raw;
					double updated_raw_subtr_unf = h_raw_subtr_unf->GetBinContent(i_trk_bin+1,i_jet_bin+1) / n_jets_unf;
					double updated_raw_subtr_unf_bbb = h_raw_subtr_unf_bbb->GetBinContent(i_trk_bin+1,i_jet_bin+1) / n_jets_unf;

					double updated_raw_err = h_raw->GetBinError(i_trk_bin+1,i_jet_bin+1) / n_jets_raw;
					double updated_raw_subtr_err = h_raw_subtr->GetBinError(i_trk_bin+1,i_jet_bin+1) / n_jets_raw;
					double updated_raw_subtr_unf_err = h_raw_subtr_unf->GetBinError(i_trk_bin+1,i_jet_bin+1) / n_jets_unf;
					double updated_raw_subtr_unf_bbb_err = h_raw_subtr_unf_bbb->GetBinError(i_trk_bin+1,i_jet_bin+1) / n_jets_unf;

					h_raw->SetBinContent(i_trk_bin+1,i_jet_bin+1, updated_raw);
					h_raw_subtr->SetBinContent(i_trk_bin+1,i_jet_bin+1, updated_raw_subtr);
					h_raw_subtr_unf->SetBinContent(i_trk_bin+1,i_jet_bin+1, updated_raw_subtr_unf);
					h_raw_subtr_unf_bbb->SetBinContent(i_trk_bin+1,i_jet_bin+1, updated_raw_subtr_unf_bbb);

					h_raw->SetBinError(i_trk_bin+1,i_jet_bin+1, updated_raw_err);
					h_raw_subtr->SetBinError(i_trk_bin+1,i_jet_bin+1, updated_raw_subtr_err);
					h_raw_subtr_unf->SetBinError(i_trk_bin+1,i_jet_bin+1, updated_raw_subtr_unf_err);
					h_raw_subtr_unf_bbb->SetBinError(i_trk_bin+1,i_jet_bin+1, updated_raw_subtr_unf_bbb_err);

					//ChPS_raw_rr
					double updated_raw_rr = h_raw_rr->GetBinContent(i_trk_bin+1,i_jet_bin+1) / n_jets_raw_rr;
					double updated_raw_rr_unf = h_raw_rr_unf->GetBinContent(i_trk_bin+1,i_jet_bin+1) / n_jets_unf_rr;
					double updated_raw_rr_unf_bbb = h_raw_rr_unf_bbb->GetBinContent(i_trk_bin+1,i_jet_bin+1) / n_jets_unf_rr;

					double updated_raw_rr_err = h_raw_rr->GetBinError(i_trk_bin+1,i_jet_bin+1) / n_jets_raw_rr;
					double updated_raw_rr_unf_err = h_raw_rr_unf->GetBinError(i_trk_bin+1,i_jet_bin+1) / n_jets_unf_rr;
					double updated_raw_rr_unf_bbb_err = h_raw_rr_unf_bbb->GetBinError(i_trk_bin+1,i_jet_bin+1) / n_jets_unf_rr;

					h_raw_rr->SetBinContent(i_trk_bin+1,i_jet_bin+1, updated_raw_rr);
					h_raw_rr_unf->SetBinContent(i_trk_bin+1,i_jet_bin+1, updated_raw_rr_unf);
					h_raw_rr_unf_bbb->SetBinContent(i_trk_bin+1,i_jet_bin+1, updated_raw_rr_unf_bbb);

					h_raw_rr->SetBinError(i_trk_bin+1,i_jet_bin+1, updated_raw_rr_err);
					h_raw_rr_unf->SetBinError(i_trk_bin+1,i_jet_bin+1, updated_raw_rr_unf_err);
					h_raw_rr_unf_bbb->SetBinError(i_trk_bin+1,i_jet_bin+1, updated_raw_rr_unf_bbb_err);


					
					//ChPS_truth
					double updated_truth = h_raw_rr->GetBinContent(i_trk_bin+1,i_jet_bin+1) / n_jets_tru;

					double updated_truth_err = h_raw_rr->GetBinError(i_trk_bin+1,i_jet_bin+1) / n_jets_tru;

					h_truth->SetBinContent(i_trk_bin+1,i_jet_bin+1, updated_truth);
					h_truth->SetBinError(i_trk_bin+1,i_jet_bin+1, updated_truth_err);

					//UE
					double updated_UE = h_UE->GetBinContent(i_trk_bin+1,i_jet_bin+1) / n_jets_raw;

					double updated_UE_err = h_UE->GetBinError(i_trk_bin+1,i_jet_bin+1) / n_jets_raw;

					h_UE->SetBinContent(i_trk_bin+1,i_jet_bin+1, updated_UE);
					h_UE->SetBinError(i_trk_bin+1,i_jet_bin+1, updated_UE_err);

				}
			}

			if (i_dR == 0)
			{

				h_raw_injet = (TH2*)h_raw->Clone(Form("raw_injet_c%i", i_cent));
				h_raw_subtr_injet = (TH2*)h_raw_subtr->Clone(Form("raw_subtr_injet_c%i", i_cent));
				h_raw_subtr_unf_injet = (TH2*)h_raw_subtr_unf->Clone(Form("raw_subtr_unf_injet_c%i", i_cent));
				h_raw_subtr_unf_bbb_injet = (TH2*)h_raw_subtr_unf_bbb->Clone(Form("raw_subtr_unf_bbb_injet_c%i", i_cent));

				h_raw_rr_injet = (TH2*)h_raw_rr->Clone(Form("raw_rr_injet_c%i", i_cent));
				h_raw_rr_unf_injet = (TH2*)h_raw_rr_unf->Clone(Form("raw_rr_unf_injet_c%i", i_cent));
				h_raw_rr_unf_bbb_injet = (TH2*)h_raw_rr_unf_bbb->Clone(Form("raw_rr_unf_bbb_injet_c%i", i_cent));

				h_UE_injet = (TH2*)h_UE->Clone(Form("UE_injet_c%i", i_cent));
			}
			if (i_dR > 0 && i_dR < 7)
			{
				h_raw_injet->Add(h_raw);
				h_raw_subtr_injet->Add(h_raw_subtr);
				h_raw_subtr_unf_injet->Add(h_raw_subtr_unf);
				h_raw_subtr_unf_bbb_injet->Add(h_raw_subtr_unf_bbb);

				h_raw_rr_injet->Add(h_raw_rr);
				h_raw_rr_unf_injet->Add(h_raw_rr_unf);
				h_raw_rr_unf_bbb_injet->Add(h_raw_rr_unf_bbb);

				h_UE_injet->Add(h_UE);
			}

			f_output->cd();

			name = Form("raw_dR%i_cent%i", i_dR, i_cent);
			h_raw->SetTitle(name.c_str());
			h_raw->Write(name.c_str());

			name = Form("raw_subtr_dR%i_cent%i", i_dR, i_cent);
			h_raw_subtr->SetTitle(name.c_str());
			h_raw_subtr->Write(name.c_str());

			name = Form("raw_subtr_unf_dR%i_cent%i", i_dR, i_cent);
			h_raw_subtr_unf->SetTitle(name.c_str());
			h_raw_subtr_unf->Write(name.c_str());

			name = Form("raw_subtr_unf_bbb_dR%i_cent%i", i_dR, i_cent);
			h_raw_subtr_unf_bbb->SetTitle(name.c_str());
			h_raw_subtr_unf_bbb->Write(name.c_str());

			name = Form("raw_rr_dR%i_cent%i", i_dR, i_cent);
			h_raw_rr->SetTitle(name.c_str());
			h_raw_rr->Write(name.c_str());

			name = Form("raw_rr_unf_dR%i_cent%i", i_dR, i_cent);
			h_raw_rr_unf->SetTitle(name.c_str());
			h_raw_rr_unf->Write(name.c_str());

			name = Form("raw_rr_unf_bbb_dR%i_cent%i", i_dR, i_cent);
			h_raw_rr_unf_bbb->SetTitle(name.c_str());
			h_raw_rr_unf_bbb->Write(name.c_str());

			name = Form("h_UE_dR%i_cent%i", i_dR, i_cent);
			h_UE->SetTitle(name.c_str());
			h_UE->Write(name.c_str());
		}

		f_output->cd();

		name = Form("raw_injet_cent%i", i_cent);
		h_raw_injet->SetTitle(name.c_str());
		h_raw_injet->Write(name.c_str());

		name = Form("raw_subtr_injet_cent%i", i_cent);
		h_raw_subtr_injet->SetTitle(name.c_str());
		h_raw_subtr_injet->Write(name.c_str());

		name = Form("raw_subtr_unf_injet_cent%i", i_cent);
		h_raw_subtr_unf_injet->SetTitle(name.c_str());
		h_raw_subtr_unf_injet->Write(name.c_str());

		name = Form("raw_subtr_unf_bbb_injet_cent%i", i_cent);
		h_raw_subtr_unf_bbb_injet->SetTitle(name.c_str());
		h_raw_subtr_unf_bbb_injet->Write(name.c_str());

		name = Form("raw_rr_injet_cent%i", i_cent);
		h_raw_rr_injet->SetTitle(name.c_str());
		h_raw_rr_injet->Write(name.c_str());

		name = Form("raw_rr_unf_injet_cent%i", i_cent);
		h_raw_rr_unf_injet->SetTitle(name.c_str());
		h_raw_rr_unf_injet->Write(name.c_str());

		name = Form("raw_rr_unf_bbb_injet_cent%i", i_cent);
		h_raw_rr_unf_bbb_injet->SetTitle(name.c_str());
		h_raw_rr_unf_bbb_injet->Write(name.c_str());

		name = Form("h_UE_injet_cent%i", i_cent);
		h_UE_injet->SetTitle(name.c_str());
		h_UE_injet->Write(name.c_str());
	}

	cout << "Done unfolding" << endl;
	return 0;
}
