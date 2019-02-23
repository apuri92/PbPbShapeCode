#include "functions/global_variables.h"


void eff_tmp()
{
	gErrorIgnoreLevel = 3001;

	int jet_pt_start = 7;
	int jet_pt_end = 11;

	int eta_start = 10;
	int eta_end = 11;

	int n_cent = 6;

	TFile *input_file = new TFile("raw_results/Perf_MC_out_histo_PbPb_5p02_r001.root");

	vector <double> eta_range;
	double x = -2.5;
	while (x <= 2.5)
	{
		eta_range.push_back(x);
		if (x < -0.9 || x >= 0.89) x = x+double(0.2);
		else x = x+double(0.3);
		if (fabs(x) < 0.0001) x = 0;
	}

	double eta_binning_new[100];
	for (int i = 0; i < eta_range.size(); i++) eta_binning_new[i] = eta_range.at(i);
	int n_eta_new = eta_range.size()-1;


	TAxis* dR_binning = (TAxis*)((TH3*)input_file->Get("h_trk_foreff_r_matched_cent0"))->GetXaxis();
	int n_dR = dR_binning->GetNbins();

	TAxis* trk_eta_binning_new = new TAxis(n_eta_new, eta_binning_new);
	int n_trk_eta_bins_new = trk_eta_binning_new->GetNbins();

	TAxis* jet_pt_binning = (TAxis*)((TH3*)input_file->Get("h_eff_matched_dR0_cent0"))->GetXaxis();
	int n_jet_pt_bins = jet_pt_binning->GetNbins();


	vector<vector<vector<vector<TH1*>>>> h_eff(n_cent, vector<vector<vector<TH1*>>> (n_dR, vector<vector<TH1*>> (n_trk_eta_bins_new, vector<TH1*>(n_jet_pt_bins))));
	TCanvas *c_eff_jetpt = new TCanvas("c_eff_jetpt","c_eff_jetpt",800,600);
	TLatex *ltx = new TLatex();
	ltx->SetTextFont(43);
	ltx->SetTextSize(12);

	string name;
	for (int i_cent = 0; i_cent < n_cent; i_cent++)
	{
		for (int i_dR = 0 ; i_dR < n_dR; i_dR++)
		{
			TH3* h_matched_3d = (TH3*)input_file->Get(Form("h_eff_matched_dR%i_cent%i", i_dR, i_cent));
			TH3* h_full_3d = (TH3*)input_file->Get(Form("h_eff_dR%i_cent%i", i_dR, i_cent));

			for (int i_eta = eta_start; i_eta < eta_end; i_eta++)
			{
				double eta_lo = trk_eta_binning_new->GetBinLowEdge(i_eta+1);
				double eta_hi = trk_eta_binning_new->GetBinUpEdge(i_eta+1);

				int eta_lo_bin = h_matched_3d->GetZaxis()->FindBin(eta_lo+0.001);
				int eta_hi_bin = h_matched_3d->GetZaxis()->FindBin(eta_hi-0.001);

//				if (eta_lo < -1.1 || eta_hi > 1.1) continue;

				int jet_itr = 0;
				for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
				{
					h_matched_3d->GetXaxis()->SetRange(i_jet+1, i_jet+1);
					h_matched_3d->GetZaxis()->SetRange(eta_hi_bin, eta_hi_bin);
					TH1* h_matched_1d = (TH1*)h_matched_3d->Project3D("y");

					h_full_3d->GetXaxis()->SetRange(i_jet+1, i_jet+1);
					h_full_3d->GetZaxis()->SetRange(eta_hi_bin, eta_hi_bin);
					TH1* h_full_1d = (TH1*)h_full_3d->Project3D("y");

//					cout << Form("ratio_dR%i_cent%i_eta%i_jet%i",i_dR, i_cent, i_eta, i_jet) << endl;

					TH1* h_ratio_1d = (TH1*)h_matched_1d->Clone(Form("ratio_dR%i_cent%i_eta%i_jet%i",i_dR, i_cent, i_eta, i_jet));
					h_ratio_1d->Divide(h_full_1d);

					h_eff[i_cent][i_dR][i_eta][i_jet] = (TH1*)h_ratio_1d->Clone(Form("h_eff_dR%i_cent%i_eta%i_jet%i", i_dR, i_cent, i_eta, i_jet));

					h_eff[i_cent][i_dR][i_eta][i_jet]->GetYaxis()->SetRangeUser(0,1.2);
					h_eff[i_cent][i_dR][i_eta][i_jet]->GetXaxis()->SetRangeUser(1,1E2);
					c_eff_jetpt->cd();
					SetHStyle_smallify(h_eff[i_cent][i_dR][i_eta][i_jet], jet_itr, 1);
					if (jet_itr == 0) h_eff[i_cent][i_dR][i_eta][i_jet]->Draw();
					else  h_eff[i_cent][i_dR][i_eta][i_jet]->Draw("same");
					gPad->SetLogx();


					jet_itr++;
				}

				string eta_label = Form("%4.2f < #eta < %4.2f", eta_lo, eta_hi);
				string centrality = num_to_cent(31, i_cent);
				string dr_label = Form("%1.2f < r < %1.2f", dR_binning->GetBinLowEdge(i_dR+1), dR_binning->GetBinUpEdge(i_dR+1));

				ltx->DrawLatexNDC(0.193, 0.90, Form("%s", dr_label.c_str()));
				ltx->DrawLatexNDC(0.193, 0.85, Form("%s", centrality.c_str()));
				ltx->DrawLatexNDC(0.193, 0.80, Form("%s", eta_label.c_str()));

				string pdf = "";
				if (i_dR == 0 && i_cent == 0 && i_eta == eta_start) pdf = "(";
				else if (i_dR == n_dR-1 && i_cent == n_cent-1 && i_eta == eta_end-1) pdf = ")";
				c_eff_jetpt->Print(Form("eff_tmp.pdf%s",pdf.c_str()), Form("Title: dR%i_cent%i_eta%i", i_dR, i_cent, i_eta));
			}
		}
	}
}
