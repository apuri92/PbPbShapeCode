#include "../functions/global_variables.h"
#include "TEnv.h"
#include "TGaxis.h"

void draw_ChPS(string config_file = "ff_config.cfg")
{
	SetAtlasStyle();
	gErrorIgnoreLevel = 3001;

	cout << "Drawing ChPS..." << endl;
	//	##############	Reading config	##############"
	TEnv *m_config = new TEnv();
	m_config->ReadFile(config_file.c_str(), EEnvLevel(1));
	m_config->Print();

	std::string dataset_type = "PbPb"; dataset_type = m_config->GetValue("dataset_type", dataset_type.c_str());
	std::string tracking_cut = "ppTight"; tracking_cut = m_config->GetValue("tracking_cut", tracking_cut.c_str());
	int centrality_scheme = 31; centrality_scheme = m_config->GetValue("centrality_scheme", centrality_scheme);
	int isMC = 1; isMC = m_config->GetValue("isMC", isMC);
	int n_unfold = 4; n_unfold = m_config->GetValue("n_unfold", n_unfold);
	bool diagnostic = false; diagnostic = m_config->GetValue("diagnostic_mode", diagnostic);

	std::string did = "data";
	if (isMC) did = "MC";
	//	##############	Config done	##############"

	TFile *f_input = new TFile(Form("unfolded_%s_%s.root",did.c_str(), dataset_type.c_str()));
	TFile *f_output = new TFile(Form("final_ChPS_%s_%s.root",did.c_str(), dataset_type.c_str()), "recreate");

	TAxis* dR_binning = (TAxis*)f_input->Get("dR_binning");
	TAxis* jetpT_binning = (TAxis*)f_input->Get("jetpT_binning");
	TAxis* trkpT_binning = (TAxis*)f_input->Get("trkpT_binning");

	f_output->cd();
	dR_binning->Write("dR_binning");
	jetpT_binning->Write("jetpT_binning");
	trkpT_binning->Write("trkpT_binning");


	int N_dR = dR_binning->GetNbins();
	int N_jetpt = jetpT_binning->GetNbins();
	int N_trkpt = trkpT_binning->GetNbins();

	//raw_0
	vector<vector<vector<TH1*>>> h_ChPS_raw (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_raw_subtr (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_raw_subtr_unf (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_raw_subtr_unf_bbb (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));

	vector<vector<vector<TH1*>>> h_ChPS_final_dR (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_noUnf_dR (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_UE_dR (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));

	vector<vector<TH1*>> h_ChPS_truth_inJet (n_cent_cuts, vector<TH1*> (N_jetpt));
	vector<vector<TH1*>> h_ChPS_final_inJet (n_cent_cuts, vector<TH1*> (N_jetpt));
	vector<vector<TH1*>> h_ChPS_ratio_final_truth_inJet (n_cent_cuts, vector<TH1*> (N_jetpt));
	vector<vector<TH1*>> h_ChPS_UE_inJet (n_cent_cuts, vector<TH1*> (N_jetpt));

	vector<vector<vector<TH1*>>> h_ChPS_ratio_subtr_raw (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_ratio_unf_subtr (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_ratio_closure (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));

	vector<vector<vector<TH1*>>> h_ChPS_diff_subtr_raw (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_diff_unf_subtr (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_diff_closure (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));

	//raw_rr
	vector<vector<vector<TH1*>>> h_ChPS_raw_rr (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_raw_rr_unf (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_raw_rr_unf_bbb (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));

	vector<vector<vector<TH1*>>> h_ChPS_ratio_unf_raw_rr (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_ratio_closure_rr (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));

	vector<vector<vector<TH1*>>> h_ChPS_diff_unf_raw_rr (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_diff_closure_rr (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));

	//truth
	vector<vector<vector<TH1*>>> h_ChPS_truth (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));

	//UE
	vector<vector<vector<TH1*>>> h_ChPS_UE (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));


	/**/

	//	vector<vector<TH1*>> h_ChPS_raw_injet (n_cent_cuts, vector<TH1*> (N_jetpt));
	//	vector<vector<TH1*>> h_ChPS_UE_injet (n_cent_cuts, vector<TH1*> (N_jetpt));
	//	vector<vector<TH1*>> h_ChPS_UE_SB_injet (n_cent_cuts, vector<TH1*> (N_jetpt));

	string name;
	string pdf_label;

	TLine *line = new TLine();
	line->SetLineColor(kBlack);

	TLatex *ltx = new TLatex();
	ltx->SetTextFont(43);
	ltx->SetTextSize(12);

	double trk_pt_lo = 1.;
	double trk_pt_hi = 150.;
	double ratio_lo = 0.;
	double ratio_hi = 2.;

	int jet_pt_start = 7;
	int jet_pt_end = 11;

	bool doSmall;
	if (dataset_type == "pp") doSmall = false;
	if (dataset_type == "PbPb") doSmall = true;


	double array_dr_bins[N_dR+1];
	for (int i_dR = 0; i_dR <= N_dR; i_dR++) array_dr_bins[i_dR] = dR_binning->GetBinLowEdge(i_dR+1);

	for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
	{
		for (int i_cent = 0; i_cent < n_cent_cuts; i_cent++)
		{
			if (dataset_type == "PbPb" && i_cent == 6) continue;
			if (dataset_type == "pp" && i_cent < 6) continue;

			for (int i_trk = 0; i_trk < N_trkpt; i_trk++)
			{
				name = Form("h_ChPS_final_dR_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_ChPS_final_dR.at(i_trk).at(i_cent).at(i_jet) = new TH1D(name.c_str(), name.c_str(), N_dR, array_dr_bins);

				name = Form("h_ChPS_noUnf_dR_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_ChPS_noUnf_dR.at(i_trk).at(i_cent).at(i_jet) = new TH1D(name.c_str(), name.c_str(), N_dR, array_dr_bins);

				if (dataset_type == "PbPb")
				{
					name = Form("h_ChPS_UE_dR_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
					h_ChPS_UE_dR.at(i_trk).at(i_cent).at(i_jet) = new TH1D(name.c_str(), name.c_str(), N_dR, array_dr_bins);
				}
			}

			for (int i_dR = 0; i_dR < N_dR; i_dR++)
			{
				if (dataset_type == "PbPb")
				{
					//UE
					name = Form("h_ChPS_UE_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet);
					h_ChPS_UE.at(i_dR).at(i_cent).at(i_jet) = (TH1*)f_input->Get(name.c_str());
					h_ChPS_UE.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]");
					h_ChPS_UE.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetTitle("#frac{1}{Area N_{Jets}} #frac{dN}{dp_{T}} (UE)");
					h_ChPS_UE.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(trk_pt_lo, trk_pt_hi);
					h_ChPS_UE.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetNdivisions(504);
					f_output->cd();
					h_ChPS_UE.at(i_dR).at(i_cent).at(i_jet)->Write(name.c_str());
				}
				//truth
				name = Form("h_ChPS_truth_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet);
				h_ChPS_truth.at(i_dR).at(i_cent).at(i_jet) = (TH1*)f_input->Get(name.c_str());
				h_ChPS_truth.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]");
				h_ChPS_truth.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetTitle("#frac{1}{Area N_{Jets}} #frac{dN}{dp_{T}} (Truth matched)");
				h_ChPS_truth.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(trk_pt_lo, trk_pt_hi);
				h_ChPS_truth.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetNdivisions(504);
				f_output->cd();
				h_ChPS_truth.at(i_dR).at(i_cent).at(i_jet)->Write(name.c_str());

				//raw
				name = Form("h_ChPS_raw_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet);
				h_ChPS_raw.at(i_dR).at(i_cent).at(i_jet) = (TH1*)f_input->Get(name.c_str());
				h_ChPS_raw.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]");
				h_ChPS_raw.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetTitle("#frac{1}{Area N_{Jets}} #frac{dN}{dp_{T}}");
				h_ChPS_raw.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(trk_pt_lo, trk_pt_hi);
				h_ChPS_raw.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetNdivisions(504);

				name = Form("h_ChPS_raw_subtr_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet);
				h_ChPS_raw_subtr.at(i_dR).at(i_cent).at(i_jet) = (TH1*)f_input->Get(name.c_str());
				h_ChPS_raw_subtr.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]");
				h_ChPS_raw_subtr.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetTitle("#frac{1}{Area N_{Jets}} #frac{dN}{dp_{T}}");
				h_ChPS_raw_subtr.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(trk_pt_lo, trk_pt_hi);
				h_ChPS_raw_subtr.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetNdivisions(504);
				name = Form("h_ChPS_noUnf_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet);
				h_ChPS_raw_subtr.at(i_dR).at(i_cent).at(i_jet)->Write(name.c_str());

				name = Form("h_ChPS_raw_subtr_unf_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet);
				h_ChPS_raw_subtr_unf.at(i_dR).at(i_cent).at(i_jet) = (TH1*)f_input->Get(name.c_str());
				h_ChPS_raw_subtr_unf.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]");
				h_ChPS_raw_subtr_unf.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetTitle("#frac{1}{Area N_{Jets}} #frac{dN}{dp_{T}}");
				h_ChPS_raw_subtr_unf.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(trk_pt_lo, trk_pt_hi);
				h_ChPS_raw_subtr_unf.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetNdivisions(504);

				name = Form("h_ChPS_raw_subtr_unf_bbb_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet);
				h_ChPS_raw_subtr_unf_bbb.at(i_dR).at(i_cent).at(i_jet) = (TH1*)f_input->Get(name.c_str());
				h_ChPS_raw_subtr_unf_bbb.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]");
				h_ChPS_raw_subtr_unf_bbb.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetTitle("#frac{1}{Area N_{Jets}} #frac{dN}{dp_{T}}");
				h_ChPS_raw_subtr_unf_bbb.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(trk_pt_lo, trk_pt_hi);
				h_ChPS_raw_subtr_unf_bbb.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetNdivisions(504);
				f_output->cd();
				name = Form("h_ChPS_final_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet);
				h_ChPS_raw_subtr_unf_bbb.at(i_dR).at(i_cent).at(i_jet)->Write(name.c_str());

				for (int i_trk = 0; i_trk < N_trkpt; i_trk++)
				{
					h_ChPS_final_dR.at(i_trk).at(i_cent).at(i_jet)->SetBinContent(i_dR+1, h_ChPS_raw_subtr_unf_bbb.at(i_dR).at(i_cent).at(i_jet)->GetBinContent(i_trk+1));
					h_ChPS_final_dR.at(i_trk).at(i_cent).at(i_jet)->SetBinError(i_dR+1, h_ChPS_raw_subtr_unf_bbb.at(i_dR).at(i_cent).at(i_jet)->GetBinError(i_trk+1));
					h_ChPS_final_dR.at(i_trk).at(i_cent).at(i_jet)->GetYaxis()->SetTitle(h_ChPS_raw_subtr_unf_bbb.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->GetTitle());
					h_ChPS_final_dR.at(i_trk).at(i_cent).at(i_jet)->GetXaxis()->SetTitle("r");

					h_ChPS_noUnf_dR.at(i_trk).at(i_cent).at(i_jet)->SetBinContent(i_dR+1, h_ChPS_raw_subtr.at(i_dR).at(i_cent).at(i_jet)->GetBinContent(i_trk+1));
					h_ChPS_noUnf_dR.at(i_trk).at(i_cent).at(i_jet)->SetBinError(i_dR+1, h_ChPS_raw_subtr.at(i_dR).at(i_cent).at(i_jet)->GetBinError(i_trk+1));
					h_ChPS_noUnf_dR.at(i_trk).at(i_cent).at(i_jet)->GetYaxis()->SetTitle(h_ChPS_raw_subtr.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->GetTitle());
					h_ChPS_noUnf_dR.at(i_trk).at(i_cent).at(i_jet)->GetXaxis()->SetTitle("r");

					if (dataset_type == "PbPb")
					{
						h_ChPS_UE_dR.at(i_trk).at(i_cent).at(i_jet)->SetBinContent(i_dR+1, h_ChPS_UE.at(i_dR).at(i_cent).at(i_jet)->GetBinContent(i_trk+1));
						h_ChPS_UE_dR.at(i_trk).at(i_cent).at(i_jet)->SetBinError(i_dR+1, h_ChPS_UE.at(i_dR).at(i_cent).at(i_jet)->GetBinError(i_trk+1));
						h_ChPS_UE_dR.at(i_trk).at(i_cent).at(i_jet)->GetYaxis()->SetTitle(h_ChPS_UE.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->GetTitle());
						h_ChPS_UE_dR.at(i_trk).at(i_cent).at(i_jet)->GetXaxis()->SetTitle("r");
					}
				}

				//raw_ratios
				name = Form("h_ChPS_ratio_subtr_raw_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet);
				h_ChPS_ratio_subtr_raw.at(i_dR).at(i_cent).at(i_jet) = (TH1*)h_ChPS_raw_subtr.at(i_dR).at(i_cent).at(i_jet)->Clone(name.c_str());
				h_ChPS_ratio_subtr_raw.at(i_dR).at(i_cent).at(i_jet)->Divide(h_ChPS_raw.at(i_dR).at(i_cent).at(i_jet));
				h_ChPS_ratio_subtr_raw.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]");
				h_ChPS_ratio_subtr_raw.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetTitle("Ratio");
				h_ChPS_ratio_subtr_raw.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(trk_pt_lo, trk_pt_hi);
				h_ChPS_ratio_subtr_raw.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetRangeUser(ratio_lo, ratio_hi);
				h_ChPS_ratio_subtr_raw.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetNdivisions(504);

				name = Form("h_ChPS_ratio_unf_subtr_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet);
				h_ChPS_ratio_unf_subtr.at(i_dR).at(i_cent).at(i_jet) = (TH1*)h_ChPS_raw_subtr_unf.at(i_dR).at(i_cent).at(i_jet)->Clone(name.c_str());
				h_ChPS_ratio_unf_subtr.at(i_dR).at(i_cent).at(i_jet)->Divide(h_ChPS_raw_subtr.at(i_dR).at(i_cent).at(i_jet));
				h_ChPS_ratio_unf_subtr.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]");
				h_ChPS_ratio_unf_subtr.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetTitle("Ratio");
				h_ChPS_ratio_unf_subtr.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(trk_pt_lo, trk_pt_hi);
				h_ChPS_ratio_unf_subtr.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetRangeUser(ratio_lo, ratio_hi);
				h_ChPS_ratio_unf_subtr.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetNdivisions(504);

				name = Form("h_ChPS_ratio_closure_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet);
				h_ChPS_ratio_closure.at(i_dR).at(i_cent).at(i_jet) = (TH1*)h_ChPS_raw_subtr_unf_bbb.at(i_dR).at(i_cent).at(i_jet)->Clone(name.c_str());
				h_ChPS_ratio_closure.at(i_dR).at(i_cent).at(i_jet)->Divide(h_ChPS_truth.at(i_dR).at(i_cent).at(i_jet));
				h_ChPS_ratio_closure.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]");
				h_ChPS_ratio_closure.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetTitle("Unfolded/Truth");
				h_ChPS_ratio_closure.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(trk_pt_lo, trk_pt_hi);
				h_ChPS_ratio_closure.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetRangeUser(ratio_lo, ratio_hi);
				h_ChPS_ratio_closure.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetNdivisions(504);
				f_output->cd();
				name = Form("h_ChPS_ratio_final_truth_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet);
				h_ChPS_ratio_closure.at(i_dR).at(i_cent).at(i_jet)->Write(name.c_str());

				//do injet stuff
				if (i_dR < 7)
				{
					if (i_dR == 0)
					{
						name = Form("h_ChPS_truth_inJet_cent%i_jetpt%i", i_cent, i_jet);
						h_ChPS_truth_inJet.at(i_cent).at(i_jet) = (TH1*)h_ChPS_truth.at(i_dR).at(i_cent).at(i_jet)->Clone(name.c_str());

						name = Form("h_ChPS_final_inJet_cent%i_jetpt%i", i_cent, i_jet);
						h_ChPS_final_inJet.at(i_cent).at(i_jet) = (TH1*)h_ChPS_raw_subtr_unf_bbb.at(i_dR).at(i_cent).at(i_jet)->Clone(name.c_str());

						if (dataset_type == "PbPb")
						{
							name = Form("h_ChPS_UE_inJet_cent%i_jetpt%i", i_cent, i_jet);
							h_ChPS_UE_inJet.at(i_cent).at(i_jet) = (TH1*)h_ChPS_UE.at(i_dR).at(i_cent).at(i_jet)->Clone(name.c_str());
						}
					}
					else
					{
						h_ChPS_truth_inJet.at(i_cent).at(i_jet)->Add(h_ChPS_truth.at(i_dR).at(i_cent).at(i_jet));
						h_ChPS_final_inJet.at(i_cent).at(i_jet)->Add(h_ChPS_raw_subtr_unf_bbb.at(i_dR).at(i_cent).at(i_jet));
						if (dataset_type == "PbPb") h_ChPS_UE_inJet.at(i_cent).at(i_jet)->Add(h_ChPS_UE.at(i_dR).at(i_cent).at(i_jet));
					}

				}


			}

			name = Form("h_ChPS_ratio_final_truth_inJet_cent%i_jetpt%i", i_cent, i_jet);
			h_ChPS_ratio_final_truth_inJet.at(i_cent).at(i_jet) = (TH1*)h_ChPS_final_inJet.at(i_cent).at(i_jet)->Clone(name.c_str());
			h_ChPS_ratio_final_truth_inJet.at(i_cent).at(i_jet)->Divide(h_ChPS_truth_inJet.at(i_cent).at(i_jet));
			h_ChPS_ratio_final_truth_inJet.at(i_cent).at(i_jet)->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]");
			h_ChPS_ratio_final_truth_inJet.at(i_cent).at(i_jet)->GetYaxis()->SetTitle("Unfolded/Truth");
			h_ChPS_ratio_final_truth_inJet.at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(trk_pt_lo, trk_pt_hi);
			h_ChPS_ratio_final_truth_inJet.at(i_cent).at(i_jet)->GetYaxis()->SetRangeUser(ratio_lo, ratio_hi);
			h_ChPS_ratio_final_truth_inJet.at(i_cent).at(i_jet)->GetYaxis()->SetNdivisions(504);

			f_output->cd();
			name = Form("h_ChPS_final_injet_cent%i_jetpt%i", i_cent, i_jet);
			h_ChPS_final_inJet.at(i_cent).at(i_jet)->Write(name.c_str());
			name = Form("h_ChPS_truth_injet_cent%i_jetpt%i", i_cent, i_jet);
			h_ChPS_truth_inJet.at(i_cent).at(i_jet)->Write(name.c_str());
			name = Form("h_ChPS_ratio_final_truth_injet_cent%i_jetpt%i", i_cent, i_jet);
			h_ChPS_ratio_final_truth_inJet.at(i_cent).at(i_jet)->Write(name.c_str());



			for (int i_trk = 0; i_trk < N_trkpt; i_trk++)
			{
				name = Form("h_ChPS_final_dR_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_ChPS_final_dR.at(i_trk).at(i_cent).at(i_jet)->Write(name.c_str());

				name = Form("h_ChPS_noUnf_dR_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_ChPS_noUnf_dR.at(i_trk).at(i_cent).at(i_jet)->Write(name.c_str());

			}




		}
	}

	if (diagnostic)
	{
		cout << "Doing Evolution plots" << endl;
		TCanvas *c_evol = new TCanvas("c_evol","c_evol",900,600);
		if (dataset_type == "pp") c_evol->SetCanvasSize(600,600);
		TLegend *legend_evol = new TLegend(0.18, 0.41, 0.45, 0.68, "","brNDC");
		legend_evol->SetTextFont(43);
		legend_evol->SetBorderSize(0);
		if (dataset_type == "pp") legend_evol->SetTextSize(12);
		if (dataset_type == "PbPb") legend_evol->SetTextSize(8);


		//Draw evolution plots
		for (int i_dR = 0; i_dR < N_dR; i_dR++)
		{
			string dr_label = Form("%1.2f < dR < %1.2f", dR_binning->GetBinLowEdge(i_dR+1), dR_binning->GetBinUpEdge(i_dR+1));
//			cout << Form("Done %s", dr_label.c_str()) << endl;

			int jet_itr = 0;
			for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
			{
				string jet_label = Form("%1.0f < p_{T}^{Jet} < %1.0f", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));

				c_evol->cd();
				c_evol->Clear();
				if (dataset_type == "PbPb") c_evol->Divide(3,2);
				bool first_pass = true;

				for (int i_cent = 0; i_cent < n_cent_cuts; i_cent++)
				{
					string centrality = num_to_cent(centrality_scheme,i_cent);

					if (dataset_type == "PbPb" && i_cent == 6) continue;
					if (dataset_type == "pp" && i_cent < 6) continue;

					SetHStyle_smallify(h_ChPS_truth.at(i_dR).at(i_cent).at(i_jet),0, doSmall);

					SetHStyle_smallify(h_ChPS_raw.at(i_dR).at(i_cent).at(i_jet),1, doSmall);
					SetHStyle_smallify(h_ChPS_raw_subtr.at(i_dR).at(i_cent).at(i_jet),2, doSmall);
					SetHStyle_smallify(h_ChPS_raw_subtr_unf.at(i_dR).at(i_cent).at(i_jet),3, doSmall);
					SetHStyle_smallify(h_ChPS_raw_subtr_unf_bbb.at(i_dR).at(i_cent).at(i_jet),4, doSmall);

					SetHStyle_smallify(h_ChPS_ratio_subtr_raw.at(i_dR).at(i_cent).at(i_jet),5, doSmall);
					SetHStyle_smallify(h_ChPS_ratio_unf_subtr.at(i_dR).at(i_cent).at(i_jet),6, doSmall);
					SetHStyle_smallify(h_ChPS_ratio_closure.at(i_dR).at(i_cent).at(i_jet),7, doSmall);

					if (i_jet == jet_pt_start && i_dR == 0 && first_pass)
					{
						legend_evol->AddEntry(h_ChPS_truth.at(i_dR).at(i_cent).at(i_jet),"Truth","lp");

						legend_evol->AddEntry(h_ChPS_raw.at(i_dR).at(i_cent).at(i_jet),"Raw","lp");
						legend_evol->AddEntry(h_ChPS_raw_subtr.at(i_dR).at(i_cent).at(i_jet),"Raw+Subtr","lp");
						legend_evol->AddEntry(h_ChPS_raw_subtr_unf.at(i_dR).at(i_cent).at(i_jet),"Raw+Subtr+Unf","lp");
						legend_evol->AddEntry(h_ChPS_raw_subtr_unf_bbb.at(i_dR).at(i_cent).at(i_jet),"Raw+Subtr+Unf+BbB","lp");

						legend_evol->AddEntry(h_ChPS_ratio_subtr_raw.at(i_dR).at(i_cent).at(i_jet),"Subtr/Raw", "lp");
						legend_evol->AddEntry(h_ChPS_ratio_unf_subtr.at(i_dR).at(i_cent).at(i_jet),"Unf/Subtr", "lp");
						legend_evol->AddEntry(h_ChPS_ratio_closure.at(i_dR).at(i_cent).at(i_jet),"Unf/Truth", "lp");
					}


					if (dataset_type == "pp") c_evol->Divide(1,2);
					if (dataset_type == "PbPb") c_evol->cd(i_cent+1)->Divide(1,2);

					if (dataset_type == "pp") c_evol->cd()->cd(1);
					if (dataset_type == "PbPb") c_evol->cd(i_cent+1)->cd(1);
					gPad->SetPad(0,0.40,0.95,0.95);
					gPad->SetTopMargin(0.05);
					gPad->SetBottomMargin(0);
					gPad->SetRightMargin(0);
					h_ChPS_raw.at(i_dR).at(i_cent).at(i_jet)->Draw("");
					h_ChPS_raw_subtr.at(i_dR).at(i_cent).at(i_jet)->Draw("same");
					h_ChPS_raw_subtr_unf.at(i_dR).at(i_cent).at(i_jet)->Draw("same");
					h_ChPS_raw_subtr_unf_bbb.at(i_dR).at(i_cent).at(i_jet)->Draw("same");
					h_ChPS_truth.at(i_dR).at(i_cent).at(i_jet)->Draw("same");
					gPad->SetLogx();
					gPad->SetLogy();

					if (dataset_type == "pp") c_evol->cd()->cd(2);
					if (dataset_type == "PbPb") c_evol->cd(i_cent+1)->cd(2);
					gPad->SetPad(0,0.0,0.95,0.40);
					gPad->SetTopMargin(0);
					gPad->SetBottomMargin(0.30);
					gPad->SetRightMargin(0);
					h_ChPS_ratio_subtr_raw.at(i_dR).at(i_cent).at(i_jet)->Draw("");
					h_ChPS_ratio_unf_subtr.at(i_dR).at(i_cent).at(i_jet)->Draw("same");
					h_ChPS_ratio_closure.at(i_dR).at(i_cent).at(i_jet)->Draw("same");
					if (dataset_type == "pp") h_ChPS_ratio_closure.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetTitleOffset(3.2);
					if (dataset_type == "PbPb") h_ChPS_ratio_closure.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetTitleOffset(5);

					gPad->SetLogx();
					line->DrawLine(trk_pt_lo, 1, trk_pt_hi, 1);

					if (dataset_type == "pp") c_evol->cd();
					if (dataset_type == "PbPb") c_evol->cd(i_cent+1);
					ltx->SetTextAlign(32);
					ltx->SetTextSize(12);
					if (dataset_type == "pp") ltx->SetTextSize(16);
					ltx->DrawLatexNDC(0.93, 0.90, Form("%s", dr_label.c_str()));
					ltx->DrawLatexNDC(0.93, 0.85, Form("%s", jet_label.c_str()));
					ltx->DrawLatexNDC(0.93, 0.80, Form("%s", centrality.c_str()));
					ltx->DrawLatexNDC(0.93, 0.98, Form("%s %s", dataset_type.c_str(), did.c_str()));
					legend_evol->Draw();

					first_pass = false;
				} // end cent loop

				pdf_label = "";
				if (i_dR == 0 && i_jet == jet_pt_start) pdf_label = "(";
				if (i_dR == N_dR-1 && i_jet == jet_pt_end-1) pdf_label = ")";
				c_evol->Print(Form("evol_%s_%s.pdf%s", dataset_type.c_str(), did.c_str(), pdf_label.c_str()), Form("Title:dR%i_jetpt%i", i_dR, i_jet));

				jet_itr++;
			} //end jet loop
		} //end dR loop
	} //end diagnostic


	//Draw Final ChPS plots
	{
		cout << "Doing Final ChPS plots in jet pT" << endl;

		TCanvas *c_ChPS_final = new TCanvas("c_ChPS_final","c_ChPS_final",900,600);
		if (dataset_type == "pp") c_ChPS_final->SetCanvasSize(600,600);
		TLegend *legend_ChPS_final = new TLegend(0.20, 0.43, 0.40, 0.60, "","brNDC");
		legend_ChPS_final->SetTextFont(43);
		legend_ChPS_final->SetBorderSize(0);
		if (dataset_type == "pp") legend_ChPS_final->SetTextSize(12);
		if (dataset_type == "PbPb") legend_ChPS_final->SetTextSize(8);

		for (int i_dR = 0; i_dR < N_dR; i_dR++)
		{
			string dr_label = Form("%1.2f < dR < %1.2f", dR_binning->GetBinLowEdge(i_dR+1), dR_binning->GetBinUpEdge(i_dR+1));
//			cout << Form("Done %s", dr_label.c_str()) << endl;

			c_ChPS_final->cd();
			c_ChPS_final->Clear();
			if (dataset_type == "PbPb") c_ChPS_final->Divide(3,2);

			bool first_pass_cent = true;
			for (int i_cent = 0; i_cent < n_cent_cuts; i_cent++)
			{
				string centrality = num_to_cent(centrality_scheme,i_cent);
				if (dataset_type == "pp" && i_cent < 6) continue;
				else if (dataset_type == "PbPb" && i_cent == 6) continue;

				if (dataset_type == "pp") c_ChPS_final->Divide(1,2);
				if (dataset_type == "PbPb") c_ChPS_final->cd(i_cent+1)->Divide(1,2);

				if (dataset_type == "pp") c_ChPS_final->cd(1);
				if (dataset_type == "PbPb") c_ChPS_final->cd(i_cent+1)->cd(1);
				gPad->SetPad(0,0.40,0.95,0.95);
				gPad->SetTopMargin(0.05);
				gPad->SetBottomMargin(0);
				gPad->SetRightMargin(0);

				if (dataset_type == "pp") c_ChPS_final->cd(2);
				if (dataset_type == "PbPb") c_ChPS_final->cd(i_cent+1)->cd(2);
				gPad->SetPad(0,0.0,0.95,0.40);
				gPad->SetTopMargin(0);
				gPad->SetBottomMargin(0.35);
				gPad->SetRightMargin(0);


				int jet_itr = 0;
				for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
				{
					string jet_label = Form("%1.0f < p_{T}^{Jet} < %1.0f", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));


					SetHStyle_smallify(h_ChPS_raw_subtr_unf_bbb.at(i_dR).at(i_cent).at(i_jet), jet_itr, doSmall);
					SetHStyle_smallify(h_ChPS_ratio_closure.at(i_dR).at(i_cent).at(i_jet), jet_itr, doSmall);

					if (i_dR == 0 && first_pass_cent) legend_ChPS_final->AddEntry(h_ChPS_raw_subtr_unf_bbb.at(i_dR).at(i_cent).at(i_jet),jet_label.c_str(),"lp");

					if (dataset_type == "pp") c_ChPS_final->cd(1);
					if (dataset_type == "PbPb") c_ChPS_final->cd(i_cent+1)->cd(1);
					if (i_jet == jet_pt_start) h_ChPS_raw_subtr_unf_bbb.at(i_dR).at(i_cent).at(i_jet)->Draw("");
					else h_ChPS_raw_subtr_unf_bbb.at(i_dR).at(i_cent).at(i_jet)->Draw("same");
					gPad->SetLogx();
					gPad->SetLogy();

					if (dataset_type == "pp") c_ChPS_final->cd(2);
					if (dataset_type == "PbPb") c_ChPS_final->cd(i_cent+1)->cd(2);
					if (i_jet == jet_pt_start) h_ChPS_ratio_closure.at(i_dR).at(i_cent).at(i_jet)->Draw("");
					else h_ChPS_ratio_closure.at(i_dR).at(i_cent).at(i_jet)->Draw("same");
					if (dataset_type == "pp") h_ChPS_ratio_closure.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetTitleOffset(3.2);
					if (dataset_type == "PbPb") h_ChPS_ratio_closure.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetTitleOffset(5);

					line->DrawLine(trk_pt_lo, 1, trk_pt_hi, 1);
					gPad->SetLogx();

					jet_itr++;
				} // end jet loop

				if (dataset_type == "pp") c_ChPS_final->cd();
				if (dataset_type == "PbPb") c_ChPS_final->cd(i_cent+1);
				ltx->SetTextAlign(32);
				ltx->SetTextSize(12);
				if (dataset_type == "pp") ltx->SetTextSize(16);
				ltx->DrawLatexNDC(0.93, 0.90, Form("%s", dr_label.c_str()));
				ltx->DrawLatexNDC(0.93, 0.85, Form("%s", centrality.c_str()));
				ltx->DrawLatexNDC(0.93, 0.98, Form("%s %s", dataset_type.c_str(), did.c_str()));

				legend_ChPS_final->Draw();
				first_pass_cent = false;
			} //end cent loop

			pdf_label = "";
			if (i_dR == 0) pdf_label = "(";
			if (i_dR == N_dR-1) pdf_label = ")";
			c_ChPS_final->Print(Form("ChPS_final_%s_%s.pdf%s", dataset_type.c_str(), did.c_str(), pdf_label.c_str()), Form("Title:dR%i", i_dR));

		} //end dR loop
	}


	//Draw Final ChPS plots in jet
	{
		cout << "Doing Final ChPS plots in jet pT for R < 0.4" << endl;

		TCanvas *c_ChPS_final_inJet = new TCanvas("c_ChPS_final_inJet","c_ChPS_final_inJet",900,600);
		if (dataset_type == "pp") c_ChPS_final_inJet->SetCanvasSize(600,600);
		TLegend *legend_ChPS_final_inJet = new TLegend(0.20, 0.43, 0.40, 0.60, "","brNDC");
		legend_ChPS_final_inJet->SetTextFont(43);
		legend_ChPS_final_inJet->SetBorderSize(0);
		if (dataset_type == "pp") legend_ChPS_final_inJet->SetTextSize(12);
		if (dataset_type == "PbPb") legend_ChPS_final_inJet->SetTextSize(8);


		c_ChPS_final_inJet->cd();
		c_ChPS_final_inJet->Clear();
		if (dataset_type == "PbPb") c_ChPS_final_inJet->Divide(3,2);

		bool first_pass_cent = true;
		for (int i_cent = 0; i_cent < n_cent_cuts; i_cent++)
		{
			string centrality = num_to_cent(centrality_scheme,i_cent);
			if (dataset_type == "pp" && i_cent < 6) continue;
			else if (dataset_type == "PbPb" && i_cent == 6) continue;

			if (dataset_type == "pp") c_ChPS_final_inJet->Divide(1,2);
			if (dataset_type == "PbPb") c_ChPS_final_inJet->cd(i_cent+1)->Divide(1,2);

			if (dataset_type == "pp") c_ChPS_final_inJet->cd(1);
			if (dataset_type == "PbPb") c_ChPS_final_inJet->cd(i_cent+1)->cd(1);
			gPad->SetPad(0,0.40,0.95,0.95);
			gPad->SetTopMargin(0.05);
			gPad->SetBottomMargin(0);
			gPad->SetRightMargin(0);

			if (dataset_type == "pp") c_ChPS_final_inJet->cd(2);
			if (dataset_type == "PbPb") c_ChPS_final_inJet->cd(i_cent+1)->cd(2);
			gPad->SetPad(0,0.0,0.95,0.40);
			gPad->SetTopMargin(0);
			gPad->SetBottomMargin(0.35);
			gPad->SetRightMargin(0);

			int jet_itr = 0;
			for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
			{
				string jet_label = Form("%1.0f < p_{T}^{Jet} < %1.0f", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));

				SetHStyle_smallify(h_ChPS_final_inJet.at(i_cent).at(i_jet), jet_itr, doSmall);
				SetHStyle_smallify(h_ChPS_ratio_final_truth_inJet.at(i_cent).at(i_jet), jet_itr, doSmall);

				if (first_pass_cent) legend_ChPS_final_inJet->AddEntry(h_ChPS_final_inJet.at(i_cent).at(i_jet),jet_label.c_str(),"lp");

				if (dataset_type == "pp") c_ChPS_final_inJet->cd(1);
				if (dataset_type == "PbPb") c_ChPS_final_inJet->cd(i_cent+1)->cd(1);
				if (i_jet == jet_pt_start) h_ChPS_final_inJet.at(i_cent).at(i_jet)->Draw("");
				else h_ChPS_final_inJet.at(i_cent).at(i_jet)->Draw("same");
				gPad->SetLogx();
				gPad->SetLogy();

				if (dataset_type == "pp") c_ChPS_final_inJet->cd(2);
				if (dataset_type == "PbPb") c_ChPS_final_inJet->cd(i_cent+1)->cd(2);
				if (i_jet == jet_pt_start) h_ChPS_ratio_final_truth_inJet.at(i_cent).at(i_jet)->Draw("");
				else h_ChPS_ratio_final_truth_inJet.at(i_cent).at(i_jet)->Draw("same");
				if (dataset_type == "pp") h_ChPS_ratio_final_truth_inJet.at(i_cent).at(i_jet)->GetXaxis()->SetTitleOffset(3.2);
				if (dataset_type == "PbPb") h_ChPS_ratio_final_truth_inJet.at(i_cent).at(i_jet)->GetXaxis()->SetTitleOffset(5);

				line->DrawLine(trk_pt_lo, 1, trk_pt_hi, 1);
				gPad->SetLogx();

				jet_itr++;
			} // end jet loop

			if (dataset_type == "pp") c_ChPS_final_inJet->cd();
			if (dataset_type == "PbPb") c_ChPS_final_inJet->cd(i_cent+1);
			ltx->SetTextAlign(32);
			ltx->SetTextSize(12);
			if (dataset_type == "pp") ltx->SetTextSize(16);
			ltx->DrawLatexNDC(0.93, 0.98, Form("%s %s", dataset_type.c_str(), did.c_str()));
			ltx->DrawLatexNDC(0.93, 0.90, Form("r < R = 0.4"));
			ltx->DrawLatexNDC(0.93, 0.85, Form("%s", centrality.c_str()));

			legend_ChPS_final_inJet->Draw();
			first_pass_cent = false;
		} //end cent loop

		pdf_label = "";
		c_ChPS_final_inJet->Print(Form("ChPS_final_inJet_%s_%s.pdf", dataset_type.c_str(), did.c_str()), Form("Title:InJet"));

	}





	//Draw Final UE plots
	if (dataset_type == "PbPb")
	{
		cout << "Doing Final UE plots in jet pT" << endl;

		TCanvas *c_ChPS_UE = new TCanvas("c_ChPS_UE","c_ChPS_UE",900,600);
		if (dataset_type == "pp") c_ChPS_UE->SetCanvasSize(600,600);
		TLegend *legend_ChPS_UE = new TLegend(0.40, 0.33, 0.60, 0.70, "","brNDC");
		legend_ChPS_UE->SetTextFont(43);
		legend_ChPS_UE->SetBorderSize(0);
		if (dataset_type == "pp") legend_ChPS_UE->SetTextSize(12);
		if (dataset_type == "PbPb") legend_ChPS_UE->SetTextSize(14);

		for (int i_dR = 0; i_dR < N_dR; i_dR++)
		{
			string dr_label = Form("%1.2f < dR < %1.2f", dR_binning->GetBinLowEdge(i_dR+1), dR_binning->GetBinUpEdge(i_dR+1));
//			cout << Form("Done %s", dr_label.c_str()) << endl;

			c_ChPS_UE->cd();
			c_ChPS_UE->Clear();
			if (dataset_type == "PbPb") c_ChPS_UE->Divide(3,2);

			bool first_pass_cent = true;
			for (int i_cent = 0; i_cent < n_cent_cuts; i_cent++)
			{
				string centrality = num_to_cent(centrality_scheme,i_cent);
				if (dataset_type == "pp" && i_cent < 6) continue;
				else if (dataset_type == "PbPb" && i_cent == 6) continue;

				int jet_itr = 0;
				for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
				{
					string jet_label = Form("%1.0f < p_{T}^{Jet} < %1.0f", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));

					SetHStyle_smallify(h_ChPS_UE.at(i_dR).at(i_cent).at(i_jet), jet_itr, doSmall);
					if (i_dR == 0 && first_pass_cent) legend_ChPS_UE->AddEntry(h_ChPS_UE.at(i_dR).at(i_cent).at(i_jet),jet_label.c_str(),"lp");

					double max = h_ChPS_UE.at(i_dR).at(i_cent).at(i_jet)->GetMaximum();
					double min = h_ChPS_UE.at(i_dR).at(i_cent).at(i_jet)->GetMinimum();
					h_ChPS_UE.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetRangeUser(0.2*min,1.05*max);
					h_ChPS_UE.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetNdivisions(504);

					if (dataset_type == "pp") c_ChPS_UE->cd();
					if (dataset_type == "PbPb") c_ChPS_UE->cd(i_cent+1);
					if (i_jet == jet_pt_start) h_ChPS_UE.at(i_dR).at(i_cent).at(i_jet)->Draw("");
					else h_ChPS_UE.at(i_dR).at(i_cent).at(i_jet)->Draw("same");

					gPad->SetLogx();
					gPad->SetLogy(0);

					jet_itr++;
				} // end jet loop

				if (dataset_type == "pp") c_ChPS_UE->cd();
				if (dataset_type == "PbPb") c_ChPS_UE->cd(i_cent+1);
				ltx->SetTextAlign(32);
				ltx->SetTextSize(12);
				ltx->DrawLatexNDC(0.93, 0.90, Form("%s", dr_label.c_str()));
				ltx->DrawLatexNDC(0.93, 0.85, Form("%s", centrality.c_str()));
				ltx->DrawLatexNDC(0.93, 0.98, Form("%s %s", dataset_type.c_str(), did.c_str()));

				legend_ChPS_UE->Draw();
				first_pass_cent = false;
			} //end cent loop

			pdf_label = "";
			if (i_dR == 0) pdf_label = "(";
			if (i_dR == N_dR-1) pdf_label = ")";
			c_ChPS_UE->Print(Form("ChPS_UE_%s_%s.pdf%s", dataset_type.c_str(), did.c_str(), pdf_label.c_str()), Form("Title:dR%i", i_dR));

		} //end dR loop
	}

	//draw as function of dR
	{
		cout << "Doing Final ChPS plots in dR" << endl;

		TCanvas *c_ChPS_dR = new TCanvas("c_ChPS_dR","c_ChPS_dR",900,600);
		if (dataset_type == "pp") c_ChPS_dR->SetCanvasSize(600,600);
		TLegend *legend_ChPS_dR = new TLegend(0.19, 0.80, 0.40, 0.92, "","brNDC");
		legend_ChPS_dR->SetTextFont(43);
		legend_ChPS_dR->SetBorderSize(0);
		if (dataset_type == "pp") legend_ChPS_dR->SetTextSize(12);
		if (dataset_type == "PbPb") legend_ChPS_dR->SetTextSize(10);

		int jet_itr = 0;
		for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
		{
			string jet_label = Form("%1.0f < p_{T}^{Jet} < %1.0f", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));
			c_ChPS_dR->cd();
			c_ChPS_dR->Clear();
			if (dataset_type == "PbPb") c_ChPS_dR->Divide(3,2);

			bool first_pass_cent = true;
			for (int i_cent = 0; i_cent < n_cent_cuts; i_cent++)
			{
				string centrality = num_to_cent(centrality_scheme,i_cent);
				if (dataset_type == "pp" && i_cent < 6) continue;
				else if (dataset_type == "PbPb" && i_cent == 6) continue;

				if (dataset_type == "pp") c_ChPS_dR->cd();
				if (dataset_type == "PbPb") c_ChPS_dR->cd(i_cent+1);

				int trk_itr = 0;
				double max = 0., min = 999., tmp;
				for (int i_trk = 0; i_trk < N_trkpt; i_trk++)
				{
					if (i_trk != 2 &&
						i_trk != 5 &&
						i_trk != 7 &&
						i_trk != 8) continue;
//					1, 5, 10, 40

					string trk_label = Form("%1.2f < p_{T}^{Trk} < %1.2f", trkpT_binning->GetBinLowEdge(i_trk+1), trkpT_binning->GetBinUpEdge(i_trk+1));
//					cout << Form("Done %i %s", i_trk, trk_label.c_str()) << endl;
					SetHStyle_smallify(h_ChPS_final_dR.at(i_trk).at(i_cent).at(i_jet), trk_itr, doSmall);
					if (jet_itr == 0 && first_pass_cent) legend_ChPS_dR->AddEntry(h_ChPS_final_dR.at(i_trk).at(i_cent).at(i_jet),trk_label.c_str(),"lp");

					h_ChPS_final_dR.at(i_trk).at(i_cent).at(i_jet)->GetYaxis()->SetRangeUser(1E-5, 1E3);
					h_ChPS_final_dR.at(i_trk).at(i_cent).at(i_jet)->GetYaxis()->SetNdivisions(504);

					if (trk_itr == 0) h_ChPS_final_dR.at(i_trk).at(i_cent).at(i_jet)->Draw("");
					else h_ChPS_final_dR.at(i_trk).at(i_cent).at(i_jet)->Draw("same");

					gPad->SetLogx(0);
					gPad->SetLogy();

					trk_itr++;

				} // end trk loop

				if (dataset_type == "pp") c_ChPS_dR->cd(1);
				if (dataset_type == "PbPb") c_ChPS_dR->cd(i_cent+1);
				ltx->SetTextAlign(32);
				ltx->SetTextSize(12);
				ltx->DrawLatexNDC(0.93, 0.90, Form("%s", jet_label.c_str()));
				ltx->DrawLatexNDC(0.93, 0.85, Form("%s", centrality.c_str()));
				ltx->DrawLatexNDC(0.93, 0.98, Form("%s %s", dataset_type.c_str(), did.c_str()));

				legend_ChPS_dR->Draw();
				first_pass_cent = false;
			} //end cent loop

			pdf_label = "";
			if (i_jet == jet_pt_start) pdf_label = "(";
			if (i_jet == jet_pt_end-1) pdf_label = ")";
			c_ChPS_dR->Print(Form("ChPS_dR_%s_%s.pdf%s", dataset_type.c_str(), did.c_str(), pdf_label.c_str()), Form("Title:jetpt%i", i_jet));
//
			jet_itr++;

		} //end jet loop
	}


	//draw UE as function of dR
	if (dataset_type == "PbPb")
	{
		cout << "Doing Final UE plots in dR" << endl;

		TCanvas *c_ChPS_dR_UE = new TCanvas("c_ChPS_dR_UE","c_ChPS_dR_UE",900,600);
		if (dataset_type == "pp") c_ChPS_dR_UE->SetCanvasSize(600,600);
		TLegend *legend_ChPS_dR_UE = new TLegend(0.19, 0.80, 0.40, 0.92, "","brNDC");
		legend_ChPS_dR_UE->SetTextFont(43);
		legend_ChPS_dR_UE->SetBorderSize(0);
		if (dataset_type == "pp") legend_ChPS_dR_UE->SetTextSize(12);
		if (dataset_type == "PbPb") legend_ChPS_dR_UE->SetTextSize(10);

		int jet_itr = 0;
		for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
		{
			string jet_label = Form("%1.0f < p_{T}^{Jet} < %1.0f", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));
			c_ChPS_dR_UE->cd();
			c_ChPS_dR_UE->Clear();
			if (dataset_type == "PbPb") c_ChPS_dR_UE->Divide(3,2);

			bool first_pass_cent = true;
			for (int i_cent = 0; i_cent < n_cent_cuts; i_cent++)
			{
				string centrality = num_to_cent(centrality_scheme,i_cent);
				if (dataset_type == "pp" && i_cent < 6) continue;
				else if (dataset_type == "PbPb" && i_cent == 6) continue;

				if (dataset_type == "pp") c_ChPS_dR_UE->cd();
				if (dataset_type == "PbPb") c_ChPS_dR_UE->cd(i_cent+1);

				int trk_itr = 0;
				double max = 0., min = 999., tmp;
				for (int i_trk = 0; i_trk < N_trkpt; i_trk++)
				{
					if (i_trk != 2 &&
						i_trk != 5 &&
						i_trk != 7 &&
						i_trk != 8) continue;
					//					1, 5, 10, 40

					string trk_label = Form("%1.2f < p_{T}^{Trk} < %1.2f", trkpT_binning->GetBinLowEdge(i_trk+1), trkpT_binning->GetBinUpEdge(i_trk+1));
					SetHStyle_smallify(h_ChPS_UE_dR.at(i_trk).at(i_cent).at(i_jet), trk_itr, doSmall);
					if (jet_itr == 0 && first_pass_cent) legend_ChPS_dR_UE->AddEntry(h_ChPS_UE_dR.at(i_trk).at(i_cent).at(i_jet),trk_label.c_str(),"lp");

					h_ChPS_UE_dR.at(i_trk).at(i_cent).at(i_jet)->GetYaxis()->SetRangeUser(1E-5, 1E3);
					h_ChPS_UE_dR.at(i_trk).at(i_cent).at(i_jet)->GetYaxis()->SetNdivisions(504);

					if (trk_itr == 0) h_ChPS_UE_dR.at(i_trk).at(i_cent).at(i_jet)->Draw("");
					else h_ChPS_UE_dR.at(i_trk).at(i_cent).at(i_jet)->Draw("same");

					gPad->SetLogx(0);
					gPad->SetLogy();

					trk_itr++;

				} // end trk loop

				if (dataset_type == "pp") c_ChPS_dR_UE->cd(1);
				if (dataset_type == "PbPb") c_ChPS_dR_UE->cd(i_cent+1);
				ltx->SetTextAlign(32);
				ltx->SetTextSize(12);
				ltx->DrawLatexNDC(0.93, 0.90, Form("%s", jet_label.c_str()));
				ltx->DrawLatexNDC(0.93, 0.85, Form("%s", centrality.c_str()));
				ltx->DrawLatexNDC(0.93, 0.98, Form("%s %s", dataset_type.c_str(), did.c_str()));

				legend_ChPS_dR_UE->Draw();
				first_pass_cent = false;
			} //end cent loop

			pdf_label = "";
			if (i_jet == jet_pt_start) pdf_label = "(";
			if (i_jet == jet_pt_end-1) pdf_label = ")";
			c_ChPS_dR_UE->Print(Form("ChPS_dR_UE_%s_%s.pdf%s", dataset_type.c_str(), did.c_str(), pdf_label.c_str()), Form("Title:jetpt%i", i_jet));
			//
			jet_itr++;

		} //end jet loop
	}
	
	//draw injet stuff, jet spectra closure


}
