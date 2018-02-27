#include "../functions/global_variables.h"
#include "TEnv.h"

void comp_ChPS(bool isMC = 0)
{

	SetAtlasStyle();
	gErrorIgnoreLevel = 3001;
	std::string did = "data";
	if (isMC) did = "MC";

	cout << "Drawing ChPS..." << endl;

	TFile *f_PbPb = new TFile(Form("final_ChPS_%s_PbPb.root", did.c_str()));
	TFile *f_pp = new TFile(Form("final_ChPS_%s_pp.root", did.c_str()));

	TFile *f_FF_PbPb = new TFile("Uncertainties_eta_4_dpt_PbPb.root");
	TFile *f_FF_pp = new TFile("Uncertainties_eta_4_dpt_pp.root");
	TFile *f_FF_ratio = new TFile("Uncertainties_eta_4_dpt_ratio.root");

	TAxis* dR_binning = (TAxis*)f_PbPb->Get("dR_binning");
	TAxis* jetpT_binning = (TAxis*)f_PbPb->Get("jetpT_binning");
	TAxis* trkpT_binning = (TAxis*)f_PbPb->Get("trkpT_binning");

	int N_dR = dR_binning->GetNbins();
	int N_jetpt = jetpT_binning->GetNbins();
	int N_trkpt = trkpT_binning->GetNbins();

	double array_dr_bins[N_dR+1];
	for (int i_dR = 0; i_dR <= N_dR; i_dR++) array_dr_bins[i_dR] = dR_binning->GetBinLowEdge(i_dR+1);
	//indR
	vector<vector<vector<TH1*>>> h_ChPS_raw_subtr_PbPb_indR (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_raw_subtr_pp_indR (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_raw_subtr_ratio_indR (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));

	vector<vector<vector<TH1*>>> h_ChPS_raw_subtr_unf_PbPb_indR (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_raw_subtr_unf_pp_indR (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_raw_subtr_unf_ratio_indR (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));

	vector<vector<vector<TH1*>>> h_ChPS_final_PbPb_indR (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_final_pp_indR (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_final_ratio_indR (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));

	//injet
	vector<vector<TH1*>> h_ChPS_final_PbPb_injet (n_cent_cuts, vector<TH1*> (N_jetpt));
	vector<vector<TH1*>> h_ChPS_final_pp_injet (n_cent_cuts, vector<TH1*> (N_jetpt));
	vector<vector<TH1*>> h_ChPS_final_ratio_injet (n_cent_cuts, vector<TH1*> (N_jetpt));

	vector<vector<TGraphErrors*>> g_ChPS_final_PbPb_injet (n_cent_cuts, vector<TGraphErrors*> (N_jetpt));
	vector<vector<TGraphErrors*>> g_ChPS_final_pp_injet (n_cent_cuts, vector<TGraphErrors*> (N_jetpt));
	vector<vector<TGraphErrors*>> g_ChPS_final_ratio_injet (n_cent_cuts, vector<TGraphErrors*> (N_jetpt));


	vector<vector<TH1*>> h_FF_final_PbPb_injet (n_cent_cuts, vector<TH1*> (N_jetpt));
	vector<vector<TH1*>> h_FF_final_pp_injet (n_cent_cuts, vector<TH1*> (N_jetpt));
	vector<vector<TH1*>> h_FF_final_ratio_injet (n_cent_cuts, vector<TH1*> (N_jetpt));

	vector<vector<TGraphErrors*>> g_FF_final_PbPb_injet (n_cent_cuts, vector<TGraphErrors*> (N_jetpt));
	vector<vector<TGraphErrors*>> g_FF_final_pp_injet (n_cent_cuts, vector<TGraphErrors*> (N_jetpt));
	vector<vector<TGraphErrors*>> g_FF_final_ratio_injet (n_cent_cuts, vector<TGraphErrors*> (N_jetpt));

	vector<vector<TGraphAsymmErrors*>> g_FF_final_PbPb_sys_injet (n_cent_cuts, vector<TGraphAsymmErrors*> (N_jetpt));
	vector<vector<TGraphAsymmErrors*>> g_FF_final_pp_sys_injet (n_cent_cuts, vector<TGraphAsymmErrors*> (N_jetpt));
	vector<vector<TGraphAsymmErrors*>> g_FF_final_ratio_sys_injet (n_cent_cuts, vector<TGraphAsymmErrors*> (N_jetpt));


	//in track pT
	vector<vector<vector<TH1*>>> h_ChPS_raw_subtr_PbPb (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_raw_subtr_pp (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_raw_subtr_ratio (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));

	vector<vector<vector<TH1*>>> h_ChPS_raw_subtr_unf_PbPb (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_raw_subtr_unf_pp (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_raw_subtr_unf_ratio (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));

	vector<vector<vector<TH1*>>> h_ChPS_final_PbPb (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_final_pp (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_final_ratio (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));

	vector<vector<vector<TH1*>>> h_ChPS_truth_PbPb (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_truth_pp (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_truth_ratio (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));



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
	double ratio_hi = 5.;

	int jet_pt_start = 7;
	int jet_pt_end = 11;

	for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
	{
		for (int i_cent = 0; i_cent < 6; i_cent++)
		{

			for (int i_trk = 0; i_trk < N_trkpt; i_trk++)
			{
				if (i_cent == 0) //get the pp plots from cent6
				{
					name = Form("h_ChPS_raw_subtr_indR_trk%i_cent%i_jetpt%i", i_trk, 6, i_jet);
					h_ChPS_raw_subtr_pp_indR.at(i_trk).at(6).at(i_jet) = (TH1*)f_pp->Get(name.c_str())->Clone(Form("pp_%s",name.c_str()));

					name = Form("h_ChPS_raw_subtr_unf_indR_trk%i_cent%i_jetpt%i", i_trk, 6, i_jet);
					h_ChPS_raw_subtr_unf_pp_indR.at(i_trk).at(6).at(i_jet) = (TH1*)f_pp->Get(name.c_str())->Clone(Form("pp_%s",name.c_str()));

					name = Form("h_ChPS_raw_subtr_unf_bbb_indR_trk%i_cent%i_jetpt%i", i_trk, 6, i_jet);
					h_ChPS_final_pp_indR.at(i_trk).at(6).at(i_jet) = (TH1*)f_pp->Get(name.c_str())->Clone(Form("pp_%s",name.c_str()));


				}

				name = Form("h_ChPS_raw_subtr_indR_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_ChPS_raw_subtr_PbPb_indR.at(i_trk).at(i_cent).at(i_jet) = (TH1*)f_PbPb->Get(name.c_str())->Clone(Form("PbPb_%s",name.c_str()));

				name = Form("h_ChPS_raw_subtr_unf_indR_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_ChPS_raw_subtr_unf_PbPb_indR.at(i_trk).at(i_cent).at(i_jet) = (TH1*)f_PbPb->Get(name.c_str())->Clone(Form("PbPb_%s",name.c_str()));

				name = Form("h_ChPS_raw_subtr_unf_bbb_indR_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_ChPS_final_PbPb_indR.at(i_trk).at(i_cent).at(i_jet) = (TH1*)f_PbPb->Get(name.c_str())->Clone(Form("PbPb_%s",name.c_str()));


				//Ratios
				name = Form("h_ChPS_raw_subtr_indR_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_ChPS_raw_subtr_ratio_indR.at(i_trk).at(i_cent).at(i_jet) = (TH1*)h_ChPS_raw_subtr_PbPb_indR.at(i_trk).at(i_cent).at(i_jet)->Clone(name.c_str());
				h_ChPS_raw_subtr_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->Divide(h_ChPS_raw_subtr_pp_indR.at(i_trk).at(6).at(i_jet));
				h_ChPS_raw_subtr_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->GetYaxis()->SetTitle("Raw+Subtr #frac{PbPb}{pp} _{(#frac{1}{Area N_{Jets}} #frac{dN}{dp_{T}})}");

				name = Form("h_ChPS_raw_subtr_unf_indR_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_ChPS_raw_subtr_unf_ratio_indR.at(i_trk).at(i_cent).at(i_jet) = (TH1*)h_ChPS_raw_subtr_unf_PbPb_indR.at(i_trk).at(i_cent).at(i_jet)->Clone(name.c_str());
				h_ChPS_raw_subtr_unf_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->Divide(h_ChPS_raw_subtr_unf_pp_indR.at(i_trk).at(6).at(i_jet));
				h_ChPS_raw_subtr_unf_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->GetYaxis()->SetTitle("Raw+Subtr+Unf #frac{PbPb}{pp} _{(#frac{1}{Area N_{Jets}} #frac{dN}{dp_{T}})}");

				name = Form("h_ChPS_final_ratio_PbPb_pp_indR_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_ChPS_final_ratio_indR.at(i_trk).at(i_cent).at(i_jet) = (TH1*)h_ChPS_final_PbPb_indR.at(i_trk).at(i_cent).at(i_jet)->Clone(name.c_str());
				h_ChPS_final_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->Divide(h_ChPS_final_pp_indR.at(i_trk).at(6).at(i_jet));
				h_ChPS_final_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->GetYaxis()->SetTitle("R_{D (p_{T}, r)}");
			}

			//injet
			if (i_cent == 0)
			{
				name = Form("h_dpt_pp_%i", i_jet);
				h_FF_final_pp_injet.at(6).at(i_jet) = (TH1*)f_FF_pp->Get(name.c_str())->Clone(Form("pp_%s",name.c_str()));
				g_FF_final_pp_injet.at(6).at(i_jet) = new TGraphErrors(h_FF_final_pp_injet.at(6).at(i_jet));
				name = Form("g_dpt_pp_uncert_%i", i_jet);
				g_FF_final_pp_sys_injet.at(6).at(i_jet) = (TGraphAsymmErrors*)f_FF_pp->Get(name.c_str())->Clone(Form("pp_%s",name.c_str()));

				name = Form("h_ChPS_final_injet_cent%i_jetpt%i", 6, i_jet);
				h_ChPS_final_pp_injet.at(6).at(i_jet) = (TH1*)f_pp->Get(name.c_str())->Clone(Form("pp_%s",name.c_str()));
				g_ChPS_final_pp_injet.at(6).at(i_jet) = new TGraphErrors(h_ChPS_final_pp_injet.at(6).at(i_jet));

			}

			name = Form("h_ChPS_final_injet_cent%i_jetpt%i", i_cent, i_jet);
			h_ChPS_final_PbPb_injet.at(i_cent).at(i_jet) = (TH1*)f_PbPb->Get(name.c_str())->Clone(Form("PbPb_%s",name.c_str()));
			g_ChPS_final_PbPb_injet.at(i_cent).at(i_jet) = new TGraphErrors(h_ChPS_final_PbPb_injet.at(i_cent).at(i_jet));

			name = Form("h_ChPS_final_ratio_PbPb_pp_cent%i_jetpt%i", i_cent, i_jet);
			h_ChPS_final_ratio_injet.at(i_cent).at(i_jet) = (TH1*)h_ChPS_final_PbPb_injet.at(i_cent).at(i_jet)->Clone(name.c_str());
			h_ChPS_final_ratio_injet.at(i_cent).at(i_jet)->Divide(h_ChPS_final_pp_injet.at(6).at(i_jet));
			h_ChPS_final_ratio_injet.at(i_cent).at(i_jet)->GetYaxis()->SetTitle("R_{D (p_{T},r)} (r < 0.4)");
			g_ChPS_final_ratio_injet.at(i_cent).at(i_jet) = new TGraphErrors(h_ChPS_final_ratio_injet.at(i_cent).at(i_jet));

			h_FF_final_PbPb_injet.at(i_cent).at(i_jet) = (TH1*)f_FF_PbPb->Get(Form("h_dpt_PbPb_%i_cent%i", i_jet, i_cent));
			g_FF_final_PbPb_injet.at(i_cent).at(i_jet) = new TGraphErrors(h_FF_final_PbPb_injet.at(i_cent).at(i_jet));
			g_FF_final_PbPb_sys_injet.at(i_cent).at(i_jet) = (TGraphAsymmErrors*)f_FF_PbPb->Get(Form("g_dpt_PbPb_uncert_%i_cent%i", i_jet, i_cent));

			h_FF_final_ratio_injet.at(i_cent).at(i_jet) = (TH1*)f_FF_ratio->Get(Form("h_dpt_ratio_%i_cent%i", i_jet, i_cent));
			g_FF_final_ratio_injet.at(i_cent).at(i_jet) = new TGraphErrors(h_FF_final_ratio_injet.at(i_cent).at(i_jet));
			g_FF_final_ratio_sys_injet.at(i_cent).at(i_jet) = (TGraphAsymmErrors*)f_FF_ratio->Get(Form("g_dpt_ratio_uncert_%i_cent%i", i_jet, i_cent));




			//add as function of jet pt later
			for (int i_dR = 0; i_dR < N_dR; i_dR++)
			{
				if (i_cent == 0) //get the pp plots from cent6
				{
					name = Form("h_ChPS_raw_subtr_dR%i_cent%i_jetpt%i", i_dR, 6, i_jet);
					h_ChPS_raw_subtr_pp.at(i_dR).at(6).at(i_jet) = (TH1*)f_pp->Get(name.c_str())->Clone(Form("pp_%s",name.c_str()));

					name = Form("h_ChPS_raw_subtr_unf_dR%i_cent%i_jetpt%i", i_dR, 6, i_jet);
					h_ChPS_raw_subtr_unf_pp.at(i_dR).at(6).at(i_jet) = (TH1*)f_pp->Get(name.c_str())->Clone(Form("pp_%s",name.c_str()));

					name = Form("h_ChPS_raw_subtr_unf_bbb_dR%i_cent%i_jetpt%i", i_dR, 6, i_jet);
					h_ChPS_final_pp.at(i_dR).at(6).at(i_jet) = (TH1*)f_pp->Get(name.c_str())->Clone(Form("pp_%s",name.c_str()));

					name = Form("h_ChPS_truth_dR%i_cent%i_jetpt%i", i_dR, 6, i_jet);
					h_ChPS_truth_pp.at(i_dR).at(6).at(i_jet) = (TH1*)f_pp->Get(name.c_str())->Clone(Form("pp_%s",name.c_str()));

				}

				name = Form("h_ChPS_raw_subtr_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet);
				h_ChPS_raw_subtr_PbPb.at(i_dR).at(i_cent).at(i_jet) = (TH1*)f_PbPb->Get(name.c_str())->Clone(Form("PbPb_%s",name.c_str()));

				name = Form("h_ChPS_raw_subtr_unf_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet);
				h_ChPS_raw_subtr_unf_PbPb.at(i_dR).at(i_cent).at(i_jet) = (TH1*)f_PbPb->Get(name.c_str())->Clone(Form("PbPb_%s",name.c_str()));

				name = Form("h_ChPS_raw_subtr_unf_bbb_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet);
				h_ChPS_final_PbPb.at(i_dR).at(i_cent).at(i_jet) = (TH1*)f_PbPb->Get(name.c_str())->Clone(Form("PbPb_%s",name.c_str()));

				name = Form("h_ChPS_truth_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet);
				h_ChPS_truth_PbPb.at(i_dR).at(i_cent).at(i_jet) = (TH1*)f_PbPb->Get(name.c_str())->Clone(Form("PbPb_%s",name.c_str()));


				//Ratios
				name = Form("h_ChPS_raw_subtr_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet);
				h_ChPS_raw_subtr_ratio.at(i_dR).at(i_cent).at(i_jet) = (TH1*)h_ChPS_raw_subtr_PbPb.at(i_dR).at(i_cent).at(i_jet)->Clone(name.c_str());
				h_ChPS_raw_subtr_ratio.at(i_dR).at(i_cent).at(i_jet)->Divide(h_ChPS_raw_subtr_pp.at(i_dR).at(6).at(i_jet));
				h_ChPS_raw_subtr_ratio.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetTitle("Raw+Subtr #frac{PbPb}{pp} _{(#frac{1}{Area N_{Jets}} #frac{dN}{dp_{T}})}");

				name = Form("h_ChPS_raw_subtr_unf_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet);
				h_ChPS_raw_subtr_unf_ratio.at(i_dR).at(i_cent).at(i_jet) = (TH1*)h_ChPS_raw_subtr_unf_PbPb.at(i_dR).at(i_cent).at(i_jet)->Clone(name.c_str());
				h_ChPS_raw_subtr_unf_ratio.at(i_dR).at(i_cent).at(i_jet)->Divide(h_ChPS_raw_subtr_unf_pp.at(i_dR).at(6).at(i_jet));
				h_ChPS_raw_subtr_unf_ratio.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetTitle("Raw+Subtr+Unf #frac{PbPb}{pp} _{(#frac{1}{Area N_{Jets}} #frac{dN}{dp_{T}})}");

				name = Form("h_ChPS_final_ratio_PbPb_pp_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet);
				h_ChPS_final_ratio.at(i_dR).at(i_cent).at(i_jet) = (TH1*)h_ChPS_final_PbPb.at(i_dR).at(i_cent).at(i_jet)->Clone(name.c_str());
				h_ChPS_final_ratio.at(i_dR).at(i_cent).at(i_jet)->Divide(h_ChPS_final_pp.at(i_dR).at(6).at(i_jet));
				h_ChPS_final_ratio.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetTitle("R_{D (p_{T}, r)}");

				name = Form("h_ChPS_truth_ratio_PbPb_pp_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet);
				h_ChPS_truth_ratio.at(i_dR).at(i_cent).at(i_jet) = (TH1*)h_ChPS_final_PbPb.at(i_dR).at(i_cent).at(i_jet)->Clone(name.c_str());
				h_ChPS_truth_ratio.at(i_dR).at(i_cent).at(i_jet)->Divide(h_ChPS_final_pp.at(i_dR).at(6).at(i_jet));
				h_ChPS_truth_ratio.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetTitle("Truth R_{D (p_{T}, r)}");

			}
		}
	}


	{
		cout << "Doing Final ChPS ratio (PbPb/pp) plots in dR" << endl;

		TCanvas *c_ChPS_raw_subtr_indR = new TCanvas("c_ChPS_raw_subtr_indR","c_ChPS_raw_subtr_indR",900,600);
		TLegend *legend_ChPS_raw_subtr_indR = new TLegend(0.19, 0.70, 0.40, 0.92, "","brNDC");
		legend_ChPS_raw_subtr_indR->SetTextFont(43);
		legend_ChPS_raw_subtr_indR->SetBorderSize(0);
		legend_ChPS_raw_subtr_indR->SetTextSize(10);

		TCanvas *c_ChPS_raw_subtr_unf_indR = new TCanvas("c_ChPS_raw_subtr_unf_indR","c_ChPS_raw_subtr_unf_indR",900,600);
		TLegend *legend_ChPS_raw_subtr_unf_indR = new TLegend(0.19, 0.70, 0.40, 0.92, "","brNDC");
		legend_ChPS_raw_subtr_unf_indR->SetTextFont(43);
		legend_ChPS_raw_subtr_unf_indR->SetBorderSize(0);
		legend_ChPS_raw_subtr_unf_indR->SetTextSize(10);

		TCanvas *c_ChPS_final_indR = new TCanvas("c_ChPS_final_indR","c_ChPS_final_indR",900,600);
		TLegend *legend_ChPS_final_indR = new TLegend(0.19, 0.60, 0.40, 0.92, "","brNDC");
		legend_ChPS_final_indR->SetTextFont(43);
		legend_ChPS_final_indR->SetBorderSize(0);
		legend_ChPS_final_indR->SetTextSize(12);

		int jet_itr = 0;
		for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
		{
			string jet_label = Form("%1.0f < p_{T}^{Jet} < %1.0f", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));

			c_ChPS_raw_subtr_indR->cd();
			c_ChPS_raw_subtr_indR->Clear();
			c_ChPS_raw_subtr_indR->Divide(3,2);

			c_ChPS_raw_subtr_unf_indR->cd();
			c_ChPS_raw_subtr_unf_indR->Clear();
			c_ChPS_raw_subtr_unf_indR->Divide(3,2);

			c_ChPS_final_indR->cd();
			c_ChPS_final_indR->Clear();
			c_ChPS_final_indR->Divide(3,2);

			bool first_pass_cent = true;
			for (int i_cent = 0; i_cent < 6; i_cent++)
			{
				string centrality = num_to_cent(31,i_cent);


				int trk_itr = 0;
				for (int i_trk = 0; i_trk < N_trkpt; i_trk++)
				{
//					if (i_trk != 2 &&
//						i_trk != 5 &&
//						i_trk != 7 &&
//						i_trk != 8) continue;
					//1, 5, 10, 40

					if (i_trk < 2 || i_trk > 8) continue;

					string trk_label = Form("%1.2f < p_{T}^{Trk} < %1.2f", trkpT_binning->GetBinLowEdge(i_trk+1), trkpT_binning->GetBinUpEdge(i_trk+1));

					SetHStyle_smallify(h_ChPS_raw_subtr_ratio_indR.at(i_trk).at(i_cent).at(i_jet), trk_itr, 1);
					SetHStyle_smallify(h_ChPS_raw_subtr_unf_ratio_indR.at(i_trk).at(i_cent).at(i_jet), trk_itr, 1);
					SetHStyle_smallify(h_ChPS_final_ratio_indR.at(i_trk).at(i_cent).at(i_jet), trk_itr, 1);


					if (jet_itr == 0 && first_pass_cent) legend_ChPS_raw_subtr_indR->AddEntry(h_ChPS_raw_subtr_ratio_indR.at(i_trk).at(i_cent).at(i_jet),trk_label.c_str(),"lp");
					if (jet_itr == 0 && first_pass_cent) legend_ChPS_raw_subtr_unf_indR->AddEntry(h_ChPS_raw_subtr_unf_ratio_indR.at(i_trk).at(i_cent).at(i_jet),trk_label.c_str(),"lp");
					if (jet_itr == 0 && first_pass_cent) legend_ChPS_final_indR->AddEntry(h_ChPS_final_ratio_indR.at(i_trk).at(i_cent).at(i_jet),trk_label.c_str(),"lp");

					h_ChPS_raw_subtr_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(0, 1.2);
					h_ChPS_raw_subtr_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->GetYaxis()->SetRangeUser(ratio_lo, ratio_hi);
					h_ChPS_raw_subtr_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->GetYaxis()->SetNdivisions(504);

					h_ChPS_raw_subtr_unf_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(0, 1.2);
					h_ChPS_raw_subtr_unf_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->GetYaxis()->SetRangeUser(ratio_lo, ratio_hi);
					h_ChPS_raw_subtr_unf_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->GetYaxis()->SetNdivisions(504);

					h_ChPS_final_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(0, 1.2);
					h_ChPS_final_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->GetYaxis()->SetRangeUser(ratio_lo, ratio_hi);
					h_ChPS_final_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->GetYaxis()->SetNdivisions(504);


					c_ChPS_raw_subtr_indR->cd(i_cent+1);
					if (trk_itr == 0) h_ChPS_raw_subtr_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->Draw("");
					else h_ChPS_raw_subtr_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->Draw("same");
					gPad->SetLogx(0);
					gPad->SetLogy(0);

					c_ChPS_raw_subtr_unf_indR->cd(i_cent+1);
					if (trk_itr == 0) h_ChPS_raw_subtr_unf_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->Draw("");
					else h_ChPS_raw_subtr_unf_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->Draw("same");
					gPad->SetLogx(0);
					gPad->SetLogy(0);

					c_ChPS_final_indR->cd(i_cent+1);
					if (trk_itr == 0) h_ChPS_final_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->Draw("");
					else h_ChPS_final_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->Draw("same");
					gPad->SetLogx(0);
					gPad->SetLogy(0);

					trk_itr++;

				} // end trk loop

				c_ChPS_raw_subtr_indR->cd(i_cent+1);
				ltx->SetTextAlign(32);
				ltx->SetTextSize(12);
				ltx->DrawLatexNDC(0.93, 0.90, Form("%s", jet_label.c_str()));
				ltx->DrawLatexNDC(0.93, 0.85, Form("%s", centrality.c_str()));
				line->DrawLine(0, 1, 1.2, 1);
				legend_ChPS_raw_subtr_indR->Draw();

				c_ChPS_raw_subtr_unf_indR->cd(i_cent+1);
				ltx->SetTextAlign(32);
				ltx->SetTextSize(12);
				ltx->DrawLatexNDC(0.93, 0.90, Form("%s", jet_label.c_str()));
				ltx->DrawLatexNDC(0.93, 0.85, Form("%s", centrality.c_str()));
				line->DrawLine(0, 1, 1.2, 1);
				legend_ChPS_raw_subtr_unf_indR->Draw();

				c_ChPS_final_indR->cd(i_cent+1);
				ltx->SetTextAlign(32);
				ltx->SetTextSize(12);
				ltx->DrawLatexNDC(0.93, 0.90, Form("%s", jet_label.c_str()));
				ltx->DrawLatexNDC(0.93, 0.85, Form("%s", centrality.c_str()));
				line->DrawLine(0, 1, 1.2, 1);
				legend_ChPS_final_indR->Draw();

				first_pass_cent = false;
			} //end cent loop

			pdf_label = "";
			if (i_jet == jet_pt_start) pdf_label = "(";
			if (i_jet == jet_pt_end-1) pdf_label = ")";
			c_ChPS_raw_subtr_indR->Print(Form("ChPS_raw_subtr_ratio_dR_%s.pdf%s", did.c_str(), pdf_label.c_str()), Form("Title:jetpt%i", i_jet));
			c_ChPS_raw_subtr_unf_indR->Print(Form("ChPS_raw_subtr_unf_ratio_dR_%s.pdf%s", did.c_str(), pdf_label.c_str()), Form("Title:jetpt%i", i_jet));
			c_ChPS_final_indR->Print(Form("ChPS_final_ratio_dR_%s.pdf%s", did.c_str(), pdf_label.c_str()), Form("Title:jetpt%i", i_jet));

			jet_itr++;
		} //end jet loop
	}

	{
		cout << "Doing in jet R < 0.4" << endl;


		TCanvas *c_ChPS_final_injet = new TCanvas("c_ChPS_final_injet","c_ChPS_final_injet",900,600);
		TLegend *legend_ChPS_final_injet = new TLegend(0.19, 0.20, 0.40, 0.30, "","brNDC");
		legend_ChPS_final_injet->SetTextFont(43);
		legend_ChPS_final_injet->SetBorderSize(0);
		legend_ChPS_final_injet->SetTextSize(12);

		c_ChPS_final_injet->cd();
		c_ChPS_final_injet->Clear();
		c_ChPS_final_injet->Divide(3,2);

		TCanvas *c_ChPS_final_PbPb_injet = new TCanvas("c_ChPS_final_PbPb_injet","c_ChPS_final_PbPb_injet",900,600);
		TLegend *legend_ChPS_final_PbPb_injet = new TLegend(0.19, 0.20, 0.40, 0.30, "","brNDC");
		legend_ChPS_final_PbPb_injet->SetTextFont(43);
		legend_ChPS_final_PbPb_injet->SetBorderSize(0);
		legend_ChPS_final_PbPb_injet->SetTextSize(12);

		c_ChPS_final_PbPb_injet->cd();
		c_ChPS_final_PbPb_injet->Clear();
		c_ChPS_final_PbPb_injet->Divide(3,2);

		TCanvas *c_ChPS_final_pp_injet = new TCanvas("c_ChPS_final_pp_injet","c_ChPS_final_pp_injet",900,600);
		TLegend *legend_ChPS_final_pp_injet = new TLegend(0.19, 0.20, 0.40, 0.30, "","brNDC");
		legend_ChPS_final_pp_injet->SetTextFont(43);
		legend_ChPS_final_pp_injet->SetBorderSize(0);
		legend_ChPS_final_pp_injet->SetTextSize(12);

		c_ChPS_final_pp_injet->cd();
		c_ChPS_final_pp_injet->Clear();


		int jet_itr = 0;
		for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
		{
			string jet_label = Form("%1.0f < p_{T}^{Jet} < %1.0f", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));

			bool first_pass_cent = true;
			for (int i_cent = 0; i_cent < 6; i_cent++)
			{
				string centrality = num_to_cent(31,i_cent);

				SetHStyle_graph_smallify(g_ChPS_final_ratio_injet.at(i_cent).at(i_jet), 0, 1);
				SetHStyle_graph_smallify(g_FF_final_ratio_injet.at(i_cent).at(i_jet), 1, 1);
				SetHStyle_graph_smallify(g_FF_final_ratio_sys_injet.at(i_cent).at(i_jet), 1, 2);

				SetHStyle_graph_smallify(g_ChPS_final_PbPb_injet.at(i_cent).at(i_jet), 0, 1);
				SetHStyle_graph_smallify(g_FF_final_PbPb_injet.at(i_cent).at(i_jet), 1, 1);
				SetHStyle_graph_smallify(g_FF_final_PbPb_sys_injet.at(i_cent).at(i_jet), 1, 2);

				SetHStyle_graph_smallify(g_ChPS_final_pp_injet.at(6).at(i_jet), 0, 0);
				SetHStyle_graph_smallify(g_FF_final_pp_injet.at(6).at(i_jet), 1, 0);
				SetHStyle_graph_smallify(g_FF_final_pp_sys_injet.at(6).at(i_jet), 1, 0);


				if (first_pass_cent && jet_itr == 0)
				{
					legend_ChPS_final_injet->AddEntry(g_ChPS_final_ratio_injet.at(i_cent).at(i_jet),"Trk jet correlations","lp");
					legend_ChPS_final_injet->AddEntry(g_FF_final_ratio_injet.at(i_cent).at(i_jet),"FF Analysis","lp");

					legend_ChPS_final_PbPb_injet->AddEntry(g_ChPS_final_PbPb_injet.at(i_cent).at(i_jet),"Trk jet correlations","lp");
					legend_ChPS_final_PbPb_injet->AddEntry(g_FF_final_PbPb_injet.at(i_cent).at(i_jet),"FF Analysis","lp");

					legend_ChPS_final_pp_injet->AddEntry(g_ChPS_final_PbPb_injet.at(i_cent).at(i_jet),"Trk jet correlations","lp");
					legend_ChPS_final_pp_injet->AddEntry(g_FF_final_PbPb_injet.at(i_cent).at(i_jet),"FF Analysis","lp");
				}


				g_ChPS_final_PbPb_injet.at(i_cent).at(i_jet)->GetYaxis()->SetTitle("D (p_{T},r #leq 0.4) _{PbPb}");
				g_ChPS_final_ratio_injet.at(i_cent).at(i_jet)->GetYaxis()->SetTitle("R_{D (p_{T},r #leq 0.4)}");
				g_FF_final_PbPb_injet.at(i_cent).at(i_jet)->GetYaxis()->SetTitle("D (p_{T}) _{PbPb}");
				g_FF_final_ratio_injet.at(i_cent).at(i_jet)->GetYaxis()->SetTitle("R_{D (p_{T})}");
				g_FF_final_ratio_sys_injet.at(i_cent).at(i_jet)->GetYaxis()->SetTitle("R_{D (p_{T},r #leq 0.4)}");

				g_FF_final_ratio_sys_injet.at(i_cent).at(i_jet)->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]");
				g_ChPS_final_ratio_injet.at(i_cent).at(i_jet)->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]");
				g_FF_final_ratio_injet.at(i_cent).at(i_jet)->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]");
				g_ChPS_final_PbPb_injet.at(i_cent).at(i_jet)->GetYaxis()->SetTitle("p_{T}^{Trk} [GeV]");
				g_FF_final_PbPb_injet.at(i_cent).at(i_jet)->GetYaxis()->SetTitle("p_{T}^{Trk} [GeV]");


				//first draw FF_uncert, then FF, then ChPS

				//ratio
				g_FF_final_ratio_sys_injet.at(i_cent).at(i_jet)->GetXaxis()->SetLimits(1,400);
				g_FF_final_ratio_sys_injet.at(i_cent).at(i_jet)->GetYaxis()->SetNdivisions(504);
				g_FF_final_ratio_sys_injet.at(i_cent).at(i_jet)->GetYaxis()->SetRangeUser(0.4,2.1);

				c_ChPS_final_injet->cd(i_cent+1);
				g_FF_final_ratio_sys_injet.at(i_cent).at(i_jet)->Draw("a PE2");
				g_FF_final_ratio_injet.at(i_cent).at(i_jet)->Draw("p same");
				g_ChPS_final_ratio_injet.at(i_cent).at(i_jet)->Draw("p same");
				gPad->SetLogx();
				gPad->SetLogy(0);
				ltx->SetTextAlign(32);
				ltx->SetTextSize(14);
				ltx->DrawLatexNDC(0.93, 0.90, Form("%s", centrality.c_str()));
				ltx->DrawLatexNDC(0.93, 0.82, Form("%s", jet_label.c_str()));
				line->DrawLine(1, 1, 400, 1);
				if (i_cent == 0) legend_ChPS_final_injet->Draw();


				//PbPb
				g_FF_final_PbPb_sys_injet.at(i_cent).at(i_jet)->GetYaxis()->SetTitle("D (p_{T},r #leq 0.4) _{PbPb}");
				g_FF_final_PbPb_sys_injet.at(i_cent).at(i_jet)->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]");
				g_FF_final_PbPb_sys_injet.at(i_cent).at(i_jet)->GetXaxis()->SetLimits(1,400);
				g_FF_final_PbPb_sys_injet.at(i_cent).at(i_jet)->GetYaxis()->SetNdivisions(504);
				g_FF_final_PbPb_sys_injet.at(i_cent).at(i_jet)->GetYaxis()->SetRangeUser(1E-7,1E2);

				c_ChPS_final_PbPb_injet->cd(i_cent+1);
				g_FF_final_PbPb_sys_injet.at(i_cent).at(i_jet)->Draw("a PE2");
				g_FF_final_PbPb_injet.at(i_cent).at(i_jet)->Draw("p same");
				g_ChPS_final_PbPb_injet.at(i_cent).at(i_jet)->Draw("p same");
				gPad->SetLogx();
				gPad->SetLogy();

				ltx->SetTextAlign(32);
				ltx->SetTextSize(14);
				ltx->DrawLatexNDC(0.93, 0.90, Form("%s", centrality.c_str()));
				ltx->DrawLatexNDC(0.93, 0.82, Form("%s", jet_label.c_str()));
				if (i_cent == 0) legend_ChPS_final_PbPb_injet->Draw();

				//pp
				if (i_cent == 0)
				{
					g_FF_final_pp_sys_injet.at(6).at(i_jet)->GetYaxis()->SetTitle("D (p_{T},r #leq 0.4) _{pp}");
					g_FF_final_pp_sys_injet.at(6).at(i_jet)->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]");
					g_FF_final_pp_sys_injet.at(6).at(i_jet)->GetXaxis()->SetLimits(1,400);
					g_FF_final_pp_sys_injet.at(6).at(i_jet)->GetYaxis()->SetNdivisions(504);
					g_FF_final_pp_sys_injet.at(6).at(i_jet)->GetYaxis()->SetRangeUser(1E-7,1E2);

					c_ChPS_final_pp_injet->cd();
					g_FF_final_pp_sys_injet.at(6).at(i_jet)->Draw("a PE2");
					g_FF_final_pp_injet.at(6).at(i_jet)->Draw("p same");
					g_ChPS_final_pp_injet.at(6).at(i_jet)->Draw("p same");

					gPad->SetLogx();
					gPad->SetLogy();

					ltx->SetTextAlign(32);
					ltx->SetTextSize(14);
					ltx->DrawLatexNDC(0.93, 0.90, "Inclusive");
					ltx->DrawLatexNDC(0.93, 0.82, Form("%s", jet_label.c_str()));
					legend_ChPS_final_pp_injet->Draw();
				}


				first_pass_cent = false;

			} // end cent loop

//			c_ChPS_final_injet->cd(i_cent+1);
//			ltx->SetTextAlign(32);
//			ltx->SetTextSize(12);
//			ltx->DrawLatexNDC(0.93, 0.85, Form("%s", centrality.c_str()));
//			line->DrawLine(1, 1, 400, 1);
//			legend_ChPS_final_injet->Draw();


			pdf_label = "";
			if (i_jet == jet_pt_start) pdf_label = "(";
			if (i_jet == jet_pt_end-1) pdf_label = ")";

			c_ChPS_final_injet->Print(Form("ChPS_final_injet_ratio_%s.pdf%s", did.c_str(), pdf_label.c_str()));
			c_ChPS_final_PbPb_injet->Print(Form("ChPS_final_injet_PbPb_%s.pdf%s", did.c_str(), pdf_label.c_str()));
			c_ChPS_final_pp_injet->Print(Form("ChPS_final_injet_pp_%s.pdf%s", did.c_str(), pdf_label.c_str()));
			jet_itr++;

		} //end jet loop

	}



	{
		cout << "Doing Final ChPS ratio (PbPb/pp) plots in jetpt" << endl;

		TCanvas *c_ChPS_raw_subtr = new TCanvas("c_ChPS_raw_subtr","c_ChPS_raw_subtr",900,600);
		TLegend *legend_ChPS_raw_subtr = new TLegend(0.19, 0.80, 0.40, 0.92, "","brNDC");
		legend_ChPS_raw_subtr->SetTextFont(43);
		legend_ChPS_raw_subtr->SetBorderSize(0);
		legend_ChPS_raw_subtr->SetTextSize(10);

		TCanvas *c_ChPS_raw_subtr_unf = new TCanvas("c_ChPS_raw_subtr_unf","c_ChPS_raw_subtr_unf",900,600);
		TLegend *legend_ChPS_raw_subtr_unf = new TLegend(0.19, 0.80, 0.40, 0.92, "","brNDC");
		legend_ChPS_raw_subtr_unf->SetTextFont(43);
		legend_ChPS_raw_subtr_unf->SetBorderSize(0);
		legend_ChPS_raw_subtr_unf->SetTextSize(10);

		TCanvas *c_ChPS_final = new TCanvas("c_ChPS_final","c_ChPS_final",900,600);
		TLegend *legend_ChPS_final = new TLegend(0.19, 0.80, 0.40, 0.92, "","brNDC");
		legend_ChPS_final->SetTextFont(43);
		legend_ChPS_final->SetBorderSize(0);
		legend_ChPS_final->SetTextSize(10);




		for (int i_dR = 0; i_dR < N_dR; i_dR++)
		{
			string dr_label = Form("%1.2f < r < %1.2f", dR_binning->GetBinLowEdge(i_dR+1), dR_binning->GetBinUpEdge(i_dR+1));

			c_ChPS_raw_subtr->cd();
			c_ChPS_raw_subtr->Clear();
			c_ChPS_raw_subtr->Divide(3,2);

			c_ChPS_raw_subtr_unf->cd();
			c_ChPS_raw_subtr_unf->Clear();
			c_ChPS_raw_subtr_unf->Divide(3,2);

			c_ChPS_final->cd();
			c_ChPS_final->Clear();
			c_ChPS_final->Divide(3,2);

			bool first_pass_cent = true;
			for (int i_cent = 0; i_cent < 6; i_cent++)
			{
				string centrality = num_to_cent(31,i_cent);


				int jet_itr = 0;
				for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
				{
					string jet_label = Form("%1.0f < p_{T}^{Jet} < %1.0f", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));

					SetHStyle_smallify(h_ChPS_raw_subtr_ratio.at(i_dR).at(i_cent).at(i_jet), jet_itr, 1);
					SetHStyle_smallify(h_ChPS_raw_subtr_unf_ratio.at(i_dR).at(i_cent).at(i_jet), jet_itr, 1);
					SetHStyle_smallify(h_ChPS_final_ratio.at(i_dR).at(i_cent).at(i_jet), jet_itr, 1);


					if (i_dR == 0 && first_pass_cent) legend_ChPS_raw_subtr->AddEntry(h_ChPS_raw_subtr_ratio.at(i_dR).at(i_cent).at(i_jet),jet_label.c_str(),"lp");
					if (i_dR == 0 && first_pass_cent) legend_ChPS_raw_subtr_unf->AddEntry(h_ChPS_raw_subtr_unf_ratio.at(i_dR).at(i_cent).at(i_jet),jet_label.c_str(),"lp");
					if (i_dR == 0 && first_pass_cent) legend_ChPS_final->AddEntry(h_ChPS_final_ratio.at(i_dR).at(i_cent).at(i_jet),jet_label.c_str(),"lp");


//					h_ChPS_raw_subtr_ratio.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(0, 1.2);
					h_ChPS_raw_subtr_ratio.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetRangeUser(ratio_lo, ratio_hi);
					h_ChPS_raw_subtr_ratio.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetNdivisions(504);

//					h_ChPS_raw_subtr_unf_ratio.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(0, 1.2);
					h_ChPS_raw_subtr_unf_ratio.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetRangeUser(ratio_lo, ratio_hi);
					h_ChPS_raw_subtr_unf_ratio.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetNdivisions(504);

//					h_ChPS_final_ratio.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(0, 1.2);
					h_ChPS_final_ratio.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetRangeUser(0.4, 2.1);
					h_ChPS_final_ratio.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetNdivisions(504);


					c_ChPS_raw_subtr->cd(i_cent+1);
					if (i_jet == jet_pt_start) h_ChPS_raw_subtr_ratio.at(i_dR).at(i_cent).at(i_jet)->Draw("");
					else h_ChPS_raw_subtr_ratio.at(i_dR).at(i_cent).at(i_jet)->Draw("same");
					gPad->SetLogx();
					gPad->SetLogy(0);

					c_ChPS_raw_subtr_unf->cd(i_cent+1);
					if (i_jet == jet_pt_start) h_ChPS_raw_subtr_unf_ratio.at(i_dR).at(i_cent).at(i_jet)->Draw("");
					else h_ChPS_raw_subtr_unf_ratio.at(i_dR).at(i_cent).at(i_jet)->Draw("same");
					gPad->SetLogx();
					gPad->SetLogy(0);

					c_ChPS_final->cd(i_cent+1);
					if (i_jet == jet_pt_start) h_ChPS_final_ratio.at(i_dR).at(i_cent).at(i_jet)->Draw("");
					else h_ChPS_final_ratio.at(i_dR).at(i_cent).at(i_jet)->Draw("same");
					gPad->SetLogx();
					gPad->SetLogy(0);

					jet_itr++;

				} // end trk loop



				c_ChPS_raw_subtr->cd(i_cent+1);
				ltx->SetTextAlign(32);
				ltx->SetTextSize(12);
				ltx->DrawLatexNDC(0.93, 0.90, Form("%s", dr_label.c_str()));
				ltx->DrawLatexNDC(0.93, 0.85, Form("%s", centrality.c_str()));
				line->DrawLine(trk_pt_lo, 1, trk_pt_hi, 1);
				legend_ChPS_raw_subtr->Draw();

				c_ChPS_raw_subtr_unf->cd(i_cent+1);
				ltx->SetTextAlign(32);
				ltx->SetTextSize(12);
				ltx->DrawLatexNDC(0.93, 0.90, Form("%s", dr_label.c_str()));
				ltx->DrawLatexNDC(0.93, 0.85, Form("%s", centrality.c_str()));
				line->DrawLine(trk_pt_lo, 1, trk_pt_hi, 1);
				legend_ChPS_raw_subtr_unf->Draw();

				c_ChPS_final->cd(i_cent+1);
				ltx->SetTextAlign(32);
				ltx->SetTextSize(12);
				ltx->DrawLatexNDC(0.93, 0.90, Form("%s", dr_label.c_str()));
				ltx->DrawLatexNDC(0.93, 0.85, Form("%s", centrality.c_str()));
				line->DrawLine(trk_pt_lo, 1, trk_pt_hi, 1);
				legend_ChPS_final->Draw();

				first_pass_cent = false;
			} //end cent loop

			pdf_label = "";
			if (i_dR == 0) pdf_label = "(";
			if (i_dR == N_dR-1) pdf_label = ")";
			c_ChPS_raw_subtr->Print(Form("ChPS_raw_subtr_ratio_%s.pdf%s", did.c_str(), pdf_label.c_str()), Form("Title:dR%i", i_dR));
			c_ChPS_raw_subtr_unf->Print(Form("ChPS_raw_subtr_unf_ratio_%s.pdf%s", did.c_str(), pdf_label.c_str()), Form("Title:dR%i", i_dR));
			c_ChPS_final->Print(Form("ChPS_final_ratio_%s.pdf%s", did.c_str(), pdf_label.c_str()), Form("Title:dR%i", i_dR));
		} //end dr loop
	}


}

