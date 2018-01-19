#include "../functions/global_variables.h"
#include "TEnv.h"
#include "TGaxis.h"

void draw_ChPS()
{
	SetAtlasStyle();
	gErrorIgnoreLevel = 3001;

	cout << "Drawing ChPS..." << endl;
	//	##############	Reading config	##############"
	TEnv *m_config = new TEnv();
	m_config->ReadFile("ff_config.cfg", EEnvLevel(1));
	m_config->Print();

	std::string dataset_type = "PbPb"; dataset_type = m_config->GetValue("dataset_type", dataset_type.c_str());
	std::string tracking_cut = "ppTight"; tracking_cut = m_config->GetValue("tracking_cut", tracking_cut.c_str());
	int centrality_scheme = 31; centrality_scheme = m_config->GetValue("centrality_scheme", centrality_scheme);
	int isMC = 1; isMC = m_config->GetValue("isMC", isMC);
	int n_unfold = 4; n_unfold = m_config->GetValue("n_unfold", n_unfold);
	bool diagnostic = true; diagnostic = m_config->GetValue("diagnostic_mode", diagnostic);

	std::string did = "data";
	if (isMC) did = "MC_JZ_comb";
	//	##############	Config done	##############"

	TFile *f_input = new TFile(Form("unfolded_%s_%s.root",did.c_str(), dataset_type.c_str()));

	TAxis* dR_binning = (TAxis*)f_input->Get("dR_binning");
	TAxis* jetpT_binning = (TAxis*)f_input->Get("jetpT_binning");
	TAxis* trkpT_binning = (TAxis*)f_input->Get("trkpT_binning");

	int N_dR = dR_binning->GetNbins();
	int N_jetpt = jetpT_binning->GetNbins();
	int N_trkpt = trkpT_binning->GetNbins();

	//raw_0
	vector<vector<vector<TH1*>>> h_ChPS_raw (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_raw_subtr (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_raw_subtr_unf (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_raw_subtr_unf_bbb (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));

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
	double ratio_lo = 0;
	double ratio_hi = 2;

	int jet_pt_start = 7;
	int jet_pt_end = 11;

	for (int i_dR = 0; i_dR < N_dR; i_dR++)
	{
		for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
		{
			for (int i_cent = 0; i_cent < n_cent_cuts; i_cent++)
			{
				//UE
				name = Form("h_ChPS_UE_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet);
				h_ChPS_UE.at(i_dR).at(i_cent).at(i_jet) = (TH1*)f_input->Get(name.c_str());
				h_ChPS_UE.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]");
				h_ChPS_UE.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetTitle("#frac{1}{Area N_{Jets}} #frac{dN}{dp_{T}} (UE)");
				h_ChPS_UE.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(trk_pt_lo, trk_pt_hi);
				h_ChPS_UE.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetNdivisions(504);

				//truth
				name = Form("h_ChPS_truth_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet);
				h_ChPS_truth.at(i_dR).at(i_cent).at(i_jet) = (TH1*)f_input->Get(name.c_str());
				h_ChPS_truth.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]");
				h_ChPS_truth.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetTitle("#frac{1}{Area N_{Jets}} #frac{dN}{dp_{T}} (Truth matched)");
				h_ChPS_truth.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(trk_pt_lo, trk_pt_hi);
				h_ChPS_truth.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetNdivisions(504);

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

				//raw_ratios
				name = Form("h_ChPS_ratio_subtr_raw_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet);
				h_ChPS_ratio_subtr_raw.at(i_dR).at(i_cent).at(i_jet) = (TH1*)h_ChPS_raw.at(i_dR).at(i_cent).at(i_jet)->Clone(name.c_str());
				h_ChPS_ratio_subtr_raw.at(i_dR).at(i_cent).at(i_jet)->Divide(h_ChPS_raw_subtr.at(i_dR).at(i_cent).at(i_jet));
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
				h_ChPS_ratio_closure.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetTitle("Ratio");
				h_ChPS_ratio_closure.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(trk_pt_lo, trk_pt_hi);
				h_ChPS_ratio_closure.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetRangeUser(ratio_lo, ratio_hi);
				h_ChPS_ratio_closure.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetNdivisions(504);

			}
		}
	}

	if (diagnostic)
	{
		TCanvas *c_evol = new TCanvas("c_evol","c_evol",800,600);
		TLegend *legend_evol = new TLegend(0.20, 0.43, 0.60, 0.73, "","brNDC");
		legend_evol->SetTextFont(43);
		legend_evol->SetBorderSize(0);
		if (dataset_type == "pp") legend_evol->SetTextSize(12);
		if (dataset_type == "PbPb") legend_evol->SetTextSize(8);


		//Draw evolution plots
		for (int i_dR = 0; i_dR < N_dR; i_dR++)
		{
			string dr_label = Form("%1.2f < dR < %1.2f", dR_binning->GetBinLowEdge(i_dR+1), dR_binning->GetBinUpEdge(i_dR+1));
			cout << Form("Done %s", dr_label.c_str()) << endl;

			int jet_itr = 0;
			for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
			{
				string jet_label = Form("%1.2f < p_{T}^{Jet} < %1.2f", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));

				c_evol->cd();
				c_evol->Clear();
				if (dataset_type == "PbPb") c_evol->Divide(3,2);
				bool first_pass = true;

				for (int i_cent = 0; i_cent < n_cent_cuts; i_cent++)
				{
					string centrality = num_to_cent(centrality_scheme,i_cent);

					if (dataset_type == "PbPb" && i_cent == 6) continue;
					if (dataset_type == "pp" && i_cent < 6) continue;

					SetHStyle_smallify(h_ChPS_truth.at(i_dR).at(i_cent).at(i_jet),0);

					SetHStyle_smallify(h_ChPS_raw.at(i_dR).at(i_cent).at(i_jet),1);
					SetHStyle_smallify(h_ChPS_raw_subtr.at(i_dR).at(i_cent).at(i_jet),2);
					SetHStyle_smallify(h_ChPS_raw_subtr_unf.at(i_dR).at(i_cent).at(i_jet),3);
					SetHStyle_smallify(h_ChPS_raw_subtr_unf_bbb.at(i_dR).at(i_cent).at(i_jet),4);

					SetHStyle_smallify(h_ChPS_ratio_subtr_raw.at(i_dR).at(i_cent).at(i_jet),5);
					SetHStyle_smallify(h_ChPS_ratio_unf_subtr.at(i_dR).at(i_cent).at(i_jet),6);
					SetHStyle_smallify(h_ChPS_ratio_closure.at(i_dR).at(i_cent).at(i_jet),7);

					if (i_jet == jet_pt_start && i_dR == 0 && first_pass)
					{
						legend_evol->AddEntry(h_ChPS_truth.at(i_dR).at(i_cent).at(i_jet),"Truth","lp");

						legend_evol->AddEntry(h_ChPS_raw.at(i_dR).at(i_cent).at(i_jet),"Raw","lp");
						legend_evol->AddEntry(h_ChPS_raw_subtr.at(i_dR).at(i_cent).at(i_jet),"Raw+Subtr","lp");
						legend_evol->AddEntry(h_ChPS_raw_subtr_unf.at(i_dR).at(i_cent).at(i_jet),"Raw+Subtr+Unf","lp");
						legend_evol->AddEntry(h_ChPS_raw_subtr_unf_bbb.at(i_dR).at(i_cent).at(i_jet),"Raw+Subtr+Unf+BbB","lp");

						legend_evol->AddEntry(h_ChPS_ratio_subtr_raw.at(i_dR).at(i_cent).at(i_jet),"Raw/Subtr", "lp");
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
					gPad->SetBottomMargin(0.3);
					gPad->SetRightMargin(0);
					h_ChPS_ratio_subtr_raw.at(i_dR).at(i_cent).at(i_jet)->Draw("");
					h_ChPS_ratio_unf_subtr.at(i_dR).at(i_cent).at(i_jet)->Draw("same");
					h_ChPS_ratio_closure.at(i_dR).at(i_cent).at(i_jet)->Draw("same");
					gPad->SetLogx();
					line->DrawLine(trk_pt_lo, 1, trk_pt_hi, 1);

					if (dataset_type == "pp") c_evol->cd();
					if (dataset_type == "PbPb") c_evol->cd(i_cent+1);
					ltx->SetTextAlign(32);
					ltx->SetTextSize(12);
					ltx->DrawLatexNDC(0.93, 0.90, Form("%s", dr_label.c_str()));
					ltx->DrawLatexNDC(0.93, 0.85, Form("%s", jet_label.c_str()));
					ltx->DrawLatexNDC(0.93, 0.80, Form("%s", centrality.c_str()));
					legend_evol->Draw();

					first_pass = false;
				} // end cent loop

				pdf_label = "";
				if (i_dR == 0 && i_jet == jet_pt_start) pdf_label = "(";
				if (i_dR == N_dR-1 && i_jet == jet_pt_end-1) pdf_label = ")";
				c_evol->Print(Form("evol_%s.pdf%s", dataset_type.c_str(), pdf_label.c_str()), Form("Title:dR%i_jet%i", i_dR, i_jet));

				jet_itr++;
			} //end jet loop
		} //end dR loop
	} //end diagnostic


	//Draw Final ChPS plots
	{
		TCanvas *c_ChPS_final = new TCanvas("c_ChPS_final","c_ChPS_final",900,600);
		TLegend *legend_ChPS_final = new TLegend(0.20, 0.43, 0.40, 0.60, "","brNDC");
		legend_ChPS_final->SetTextFont(43);
		legend_ChPS_final->SetBorderSize(0);
		if (dataset_type == "pp") legend_ChPS_final->SetTextSize(12);
		if (dataset_type == "PbPb") legend_ChPS_final->SetTextSize(8);

		for (int i_dR = 0; i_dR < N_dR; i_dR++)
		{
			string dr_label = Form("%1.2f < dR < %1.2f", dR_binning->GetBinLowEdge(i_dR+1), dR_binning->GetBinUpEdge(i_dR+1));
			cout << Form("Done %s", dr_label.c_str()) << endl;

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
				gPad->SetBottomMargin(0.3);
				gPad->SetRightMargin(0);


				int jet_itr = 0;
				for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
				{
					string jet_label = Form("%1.2f < p_{T}^{Jet} < %1.2f", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));


					SetHStyle_smallify(h_ChPS_raw_subtr_unf_bbb.at(i_dR).at(i_cent).at(i_jet), jet_itr);
					SetHStyle_smallify(h_ChPS_ratio_closure.at(i_dR).at(i_cent).at(i_jet), jet_itr);
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
					gPad->SetLogx();

					jet_itr++;
				} // end jet loop

				if (dataset_type == "pp") c_ChPS_final->cd();
				if (dataset_type == "PbPb") c_ChPS_final->cd(i_cent+1);
				ltx->SetTextAlign(32);
				ltx->SetTextSize(12);
				ltx->DrawLatexNDC(0.93, 0.90, Form("%s", dr_label.c_str()));
				ltx->DrawLatexNDC(0.93, 0.85, Form("%s", centrality.c_str()));

				legend_ChPS_final->Draw();
				first_pass_cent = false;
			} //end cent loop

			pdf_label = "";
			if (i_dR == 0) pdf_label = "(";
			if (i_dR == N_dR-1) pdf_label = ")";
			c_ChPS_final->Print(Form("ChPS_final_%s.pdf%s", dataset_type.c_str(), pdf_label.c_str()), Form("Title:dR%i", i_dR));

		} //end dR loop
	}








//	for (int i_dR = 0; i_dR < N_dR; i_dR++)
//	{
//
//		double dR_lo = dR_binning->GetBinLowEdge(i_dR+1);
//		double dR_hi = dR_binning->GetBinUpEdge(i_dR+1);
//		double area = TMath::Pi() * ((dR_hi*dR_hi) - (dR_lo*dR_lo));
//		double area_injet = TMath::Pi() * (0.4*0.4);
//		string dr_label = Form("%1.2f < dR < %1.2f", dR_lo, dR_hi);
//		cout << Form("Done %s", dr_label.c_str()) << endl;
//
//		int jet_itr = 0;
//
//		c_ChPS_final->Clear();
//		c_ChPS_final->cd();
//		c_ChPS_final->Divide(3,2);
//
//		for (int i_jet = 0; i_jet < N_jetpt; i_jet++)
//		{
//			double jet_pt_lo = jetpT_binning->GetBinLowEdge(i_jet+1);
//			double jet_pt_hi = jetpT_binning->GetBinUpEdge(i_jet+1);
//			if (jet_pt_lo < 125. || jet_pt_hi > 350) continue;
//			string jet_label = Form("%1.2f < p_{T}^{Jet} < %1.2f", jet_pt_lo, jet_pt_hi);
//
//			int style = 0;
//
//			c_evol->cd();
//			c_evol->Clear();
//			c_evol->Divide(3,2);
//
//			for (int i_cent = 0; i_cent < n_cent_cuts; i_cent++)
//			{
//				string centrality = num_to_cent(centrality_scheme,i_cent);
//
//				if (dataset_type == "pp" && i_cent < 5)
//				{
//					c_evol->cd();
//					c_evol->Clear();
//					continue;
//				}
//				if (dataset_type == "PbPb" && i_cent == 6) continue;
//
//
//
//				double trk_pt_lo = 1.;
//				double trk_pt_hi = 300.;
//				double ratio_lo = 0;
//				double ratio_hi = 2;
//
//				//UE
//				name = Form("h_ChPS_UE_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet);
//				h_ChPS_UE.at(i_dR).at(i_cent).at(i_jet) = (TH1*)f_input->Get(name.c_str());
//				h_ChPS_UE.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]");
//				h_ChPS_UE.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetTitle("#frac{1}{Area N_{Jets}} #frac{dN}{dp_{T}} (UE)");
//				h_ChPS_UE.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(trk_pt_lo, trk_pt_hi);
//
//
//				//truth
//				name = Form("h_ChPS_truth_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet);
//				h_ChPS_truth.at(i_dR).at(i_cent).at(i_jet) = (TH1*)f_input->Get(name.c_str());
//				h_ChPS_truth.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]");
//				h_ChPS_truth.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetTitle("#frac{1}{Area N_{Jets}} #frac{dN}{dp_{T}} (Truth matched)");
//				h_ChPS_truth.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(trk_pt_lo, trk_pt_hi);
//
//				//raw
//				name = Form("h_ChPS_raw_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet);
//				h_ChPS_raw.at(i_dR).at(i_cent).at(i_jet) = (TH1*)f_input->Get(name.c_str());
//				h_ChPS_raw.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]");
//				h_ChPS_raw.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetTitle("#frac{1}{Area N_{Jets}} #frac{dN}{dp_{T}}");
//				h_ChPS_raw.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(trk_pt_lo, trk_pt_hi);
//
//				name = Form("h_ChPS_raw_subtr_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet);
//				h_ChPS_raw_subtr.at(i_dR).at(i_cent).at(i_jet) = (TH1*)f_input->Get(name.c_str());
//				h_ChPS_raw_subtr.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]");
//				h_ChPS_raw_subtr.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetTitle("#frac{1}{Area N_{Jets}} #frac{dN}{dp_{T}}");
//				h_ChPS_raw_subtr.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(trk_pt_lo, trk_pt_hi);
//
//				name = Form("h_ChPS_raw_subtr_unf_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet);
//				h_ChPS_raw_subtr_unf.at(i_dR).at(i_cent).at(i_jet) = (TH1*)f_input->Get(name.c_str());
//				h_ChPS_raw_subtr_unf.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]");
//				h_ChPS_raw_subtr_unf.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetTitle("#frac{1}{Area N_{Jets}} #frac{dN}{dp_{T}}");
//				h_ChPS_raw_subtr_unf.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(trk_pt_lo, trk_pt_hi);
//
//				name = Form("h_ChPS_raw_subtr_unf_bbb_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet);
//				h_ChPS_raw_subtr_unf_bbb.at(i_dR).at(i_cent).at(i_jet) = (TH1*)f_input->Get(name.c_str());
//				h_ChPS_raw_subtr_unf_bbb.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]");
//				h_ChPS_raw_subtr_unf_bbb.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetTitle("#frac{1}{Area N_{Jets}} #frac{dN}{dp_{T}}");
//				h_ChPS_raw_subtr_unf_bbb.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(trk_pt_lo, trk_pt_hi);
//
//				//raw_ratios
//				name = Form("h_ChPS_ratio_subtr_raw_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet);
//				h_ChPS_ratio_subtr_raw.at(i_dR).at(i_cent).at(i_jet) = (TH1*)h_ChPS_raw.at(i_dR).at(i_cent).at(i_jet)->Clone(name.c_str());
//				h_ChPS_ratio_subtr_raw.at(i_dR).at(i_cent).at(i_jet)->Divide(h_ChPS_raw_subtr.at(i_dR).at(i_cent).at(i_jet));
//				h_ChPS_ratio_subtr_raw.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]");
//				h_ChPS_ratio_subtr_raw.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetTitle("Ratio");
//				h_ChPS_ratio_subtr_raw.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(trk_pt_lo, trk_pt_hi);
//				h_ChPS_ratio_subtr_raw.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetRangeUser(ratio_lo, ratio_hi);
//
//				name = Form("h_ChPS_ratio_unf_subtr_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet);
//				h_ChPS_ratio_unf_subtr.at(i_dR).at(i_cent).at(i_jet) = (TH1*)h_ChPS_raw_subtr_unf.at(i_dR).at(i_cent).at(i_jet)->Clone(name.c_str());
//				h_ChPS_ratio_unf_subtr.at(i_dR).at(i_cent).at(i_jet)->Divide(h_ChPS_raw_subtr.at(i_dR).at(i_cent).at(i_jet));
//				h_ChPS_ratio_unf_subtr.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]");
//				h_ChPS_ratio_unf_subtr.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetTitle("Ratio");
//				h_ChPS_ratio_unf_subtr.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(trk_pt_lo, trk_pt_hi);
//				h_ChPS_ratio_unf_subtr.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetRangeUser(ratio_lo, ratio_hi);
//
//				name = Form("h_ChPS_ratio_closure_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet);
//				h_ChPS_ratio_closure.at(i_dR).at(i_cent).at(i_jet) = (TH1*)h_ChPS_raw_subtr_unf_bbb.at(i_dR).at(i_cent).at(i_jet)->Clone(name.c_str());
//				h_ChPS_ratio_closure.at(i_dR).at(i_cent).at(i_jet)->Divide(h_ChPS_truth.at(i_dR).at(i_cent).at(i_jet));
//				h_ChPS_ratio_closure.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]");
//				h_ChPS_ratio_closure.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetTitle("Ratio");
//				h_ChPS_ratio_closure.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(trk_pt_lo, trk_pt_hi);
//				h_ChPS_ratio_closure.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetRangeUser(ratio_lo, ratio_hi);
//
//				//Evol canvas
//				SetHStyle_smallify(h_ChPS_truth.at(i_dR).at(i_cent).at(i_jet),0);
//
//				SetHStyle_smallify(h_ChPS_raw.at(i_dR).at(i_cent).at(i_jet),1);
//				SetHStyle_smallify(h_ChPS_raw_subtr.at(i_dR).at(i_cent).at(i_jet),2);
//				SetHStyle_smallify(h_ChPS_raw_subtr_unf.at(i_dR).at(i_cent).at(i_jet),3);
//				SetHStyle_smallify(h_ChPS_raw_subtr_unf_bbb.at(i_dR).at(i_cent).at(i_jet),4);
//
//				SetHStyle_smallify(h_ChPS_ratio_subtr_raw.at(i_dR).at(i_cent).at(i_jet),5);
//				SetHStyle_smallify(h_ChPS_ratio_unf_subtr.at(i_dR).at(i_cent).at(i_jet),6);
//				SetHStyle_smallify(h_ChPS_ratio_closure.at(i_dR).at(i_cent).at(i_jet),7);
//
//				if (jet_itr == 0 && i_cent == 0 && i_dR == 0 )
//				{
//					legend_evol->AddEntry(h_ChPS_truth.at(i_dR).at(i_cent).at(i_jet),"Truth","lp");
//
//					legend_evol->AddEntry(h_ChPS_raw.at(i_dR).at(i_cent).at(i_jet),"Raw","lp");
//					legend_evol->AddEntry(h_ChPS_raw_subtr.at(i_dR).at(i_cent).at(i_jet),"Raw+Subtr","lp");
//					legend_evol->AddEntry(h_ChPS_raw_subtr_unf.at(i_dR).at(i_cent).at(i_jet),"Raw+Subtr+Unf","lp");
//					legend_evol->AddEntry(h_ChPS_raw_subtr_unf_bbb.at(i_dR).at(i_cent).at(i_jet),"Raw+Subtr+Unf+BbB","lp");
//
//					legend_evol->AddEntry(h_ChPS_ratio_subtr_raw.at(i_dR).at(i_cent).at(i_jet),"Raw/Subtr", "lp");
//					legend_evol->AddEntry(h_ChPS_ratio_unf_subtr.at(i_dR).at(i_cent).at(i_jet),"Unf/Subtr", "lp");
//					legend_evol->AddEntry(h_ChPS_ratio_closure.at(i_dR).at(i_cent).at(i_jet),"Unf/Truth", "lp");
//				}
//
//				c_evol->cd(i_cent+1);
//				gPad->Divide(1,2);
//
//				c_evol->cd(i_cent+1)->cd(1);
//				gPad->SetPad(0,0.40,1.0,1.0);
//				gPad->SetTopMargin(0.05);
//				gPad->SetBottomMargin(0);
//				gPad->SetRightMargin(0);
//				h_ChPS_raw.at(i_dR).at(i_cent).at(i_jet)->Draw("");
//				h_ChPS_raw_subtr.at(i_dR).at(i_cent).at(i_jet)->Draw("same");
//				h_ChPS_raw_subtr_unf.at(i_dR).at(i_cent).at(i_jet)->Draw("same");
//				h_ChPS_raw_subtr_unf_bbb.at(i_dR).at(i_cent).at(i_jet)->Draw("same");
//				h_ChPS_truth.at(i_dR).at(i_cent).at(i_jet)->Draw("same");
//				gPad->SetLogx();
//				gPad->SetLogy();
//
//				c_evol->cd(i_cent+1)->cd(2);
//				gPad->SetPad(0,0.0,1,0.40);
//				gPad->SetTopMargin(0);
//				gPad->SetBottomMargin(0.3);
//				gPad->SetRightMargin(0);
//				h_ChPS_ratio_subtr_raw.at(i_dR).at(i_cent).at(i_jet)->Draw("");
//				h_ChPS_ratio_unf_subtr.at(i_dR).at(i_cent).at(i_jet)->Draw("same");
//				h_ChPS_ratio_closure.at(i_dR).at(i_cent).at(i_jet)->Draw("same");
//				gPad->SetLogx();
//				line->DrawLine(trk_pt_lo, 1, trk_pt_hi, 1);
//				c_evol->cd(i_cent+1);
//				ltx->SetTextAlign(32);
//				ltx->SetTextSize(12);
//				ltx->DrawLatexNDC(0.95, 0.92, Form("%s", dr_label.c_str()));
//				ltx->DrawLatexNDC(0.95, 0.85, Form("%s", jet_label.c_str()));
//				ltx->DrawLatexNDC(0.95, 0.79, Form("%s", centrality.c_str()));
//				legend_evol->Draw();
//
//
//
//				//final_ChPS
//				SetHStyle_smallify(h_ChPS_raw_subtr_unf_bbb.at(i_dR).at(i_cent).at(i_jet),jet_itr);
////				cout << h_ChPS_ratio_subtr_raw.at(i_dR).at(i_cent).at(i_jet)->GetName() << " " << jet_itr << endl;
//				c_ChPS_final->cd(i_cent+1);
//				if (jet_itr == 0) h_ChPS_raw_subtr_unf_bbb.at(i_dR).at(i_cent).at(i_jet)->Draw("");
//				else h_ChPS_raw_subtr_unf_bbb.at(i_dR).at(i_cent).at(i_jet)->Draw("same");
//				gPad->SetLogx();
//				gPad->SetLogy();
//				ltx->DrawLatexNDC(0.95, 0.92, Form("%s", dr_label.c_str()));
//				ltx->DrawLatexNDC(0.95, 0.85, Form("%s", centrality.c_str()));
//
//
//
//			} // end cent loop
//
//			pdf_label = "";
//			if (i_dR == 0 && jet_itr == 0) pdf_label = "(";
//			c_evol->Print(Form("evol_%s.pdf%s", dataset_type.c_str(), pdf_label.c_str()), Form("Title:dR%i_jet%i", i_dR, i_jet));
//
//			jet_itr++;
//		} //end jet loop
//
//		pdf_label = "";
//		if (i_dR == 0 ) pdf_label = "(";
//		c_ChPS_final->Print(Form("ChPS_final_%s.pdf%s", dataset_type.c_str(), pdf_label.c_str()), Form("Title:dR%i", i_dR));
//
//
//
//	} //end dR loop
//
//	c_evol->Clear();
//	c_evol->Print(Form("evol_%s.pdf)", dataset_type.c_str()), Form("Title: End"));
//
//	c_ChPS_final->Clear();
//	c_ChPS_final->Print(Form("ChPS_final_%s.pdf)", dataset_type.c_str()), Form("Title: End"));










	/*

	//	double max = 0, min = 0;
	//
	//	bool doing_injet = true;
	//	for (int i_dR = 0; i_dR < N_dR; i_dR++)
	//	{
	//		double dR_lo = dR_binning->GetBinLowEdge(i_dR+1);
	//		double dR_hi = dR_binning->GetBinUpEdge(i_dR+1);
	//		double area = TMath::Pi() * ((dR_hi*dR_hi) - (dR_lo*dR_lo));
	//		double area_injet = TMath::Pi() * (0.4*0.4);
	//
	//		c1->Divide(3,2);
	////		c8->Divide(3,2);
	////
	////		c2->Divide(3,2);
	////		c7->Divide(3,2);
	////		c3->Divide(3,2);
	////		if (doing_injet) c4->Divide(3,2);
	////		if (doing_injet) c6->Divide(3,2);
	////		c5->Divide(3,2);
	//
	//		for (int i_cent = 0; i_cent < n_cent_cuts; i_cent++)
	//		{
	//
	//
	//			string centrality = num_to_cent(centrality_scheme,i_cent);
	//			name = Form("h_%s_cent%i_dR%i",ChPS_raw_type.c_str(), i_cent, i_dR);
	//            cout << name << endl;
	//            h_raw.at(i_cent).at(i_dR) = (TH2*)f_input->Get(name.c_str())->Clone(Form("c_%s",name.c_str()));
	//			h_raw.at(i_cent).at(i_dR)->SetName(name.c_str());
	//			h_raw.at(i_cent).at(i_dR)->Sumw2();
	//
	//
	//			name = Form("h_true_ChPS_cent%i_dR%i",i_cent, i_dR);
	//			h_true.at(i_cent).at(i_dR) = (TH2*)f_input->Get(name.c_str())->Clone(Form("c_%s",name.c_str()));
	//			h_true.at(i_cent).at(i_dR)->SetName(name.c_str());
	//			h_true.at(i_cent).at(i_dR)->Sumw2();
	//
	//			name = Form("h_unfolded_%s_cent%i_dR%i", ChPS_raw_type.c_str(), i_cent, i_dR);
	//			h_unfolded.at(i_cent).at(i_dR) = (TH2*)f_input->Get(name.c_str())->Clone(Form("c_%s",name.c_str()));
	//			h_unfolded.at(i_cent).at(i_dR)->SetName(name.c_str());
	//			h_unfolded.at(i_cent).at(i_dR)->Sumw2();
	//
	//			name = Form("h_UE_cent%i_dR%i", i_cent, i_dR);
	//			h_UE.at(i_cent).at(i_dR) = (TH2*)f_input->Get(name.c_str())->Clone(Form("c_%s",name.c_str()));
	//			h_UE.at(i_cent).at(i_dR)->SetName(name.c_str());
	//			h_UE.at(i_cent).at(i_dR)->Sumw2();
	//
	//			if (doing_injet)
	//			{
	//				name = Form("h_%s_injet_cent%i",ChPS_raw_type.c_str(), i_cent);
	//				h_raw_injet.at(i_cent) = (TH2*)f_input->Get(name.c_str())->Clone(Form("c_%s",name.c_str()));
	//				h_raw_injet.at(i_cent)->SetName(name.c_str());
	//				h_raw_injet.at(i_cent)->Sumw2();
	//
	//				name = Form("h_UE_injet_cent%i", i_cent);
	//				h_UE_injet.at(i_cent) = (TH2*)f_input->Get(name.c_str())->Clone(Form("c_%s",name.c_str()));
	//				h_UE_injet.at(i_cent)->SetName(name.c_str());
	//				h_UE_injet.at(i_cent)->Sumw2();
	//			}
	//
	//			int jet_pt_start = 7;
	//			int jet_pt_end = 10;
	//			int style = 0;
	//
	//			c1->cd(i_cent+1);
	//			gPad->Divide(1,2);
	//
	//			c2->cd(i_cent+1);
	//			gPad->Divide(1,2);
	//
	//			c8->cd(i_cent+1);
	//
	//			for (int i_jet = 0; i_jet < N_jetpt; i_jet++)
	//			{
	//
	//
	//
	//				name = Form("h_ChPS_raw_cent%i_dR%i_jetpt%i",i_cent, i_dR, i_jet);
	//				h_ChPS_raw.at(i_dR).at(i_cent).at(i_jet) = (TH1*)h_raw.at(i_cent).at(i_dR)->ProjectionX(name.c_str(), i_jet+1,i_jet+1);
	//				h_ChPS_raw.at(i_dR).at(i_cent).at(i_jet)->Scale(1.,"width");
	//				h_ChPS_raw.at(i_dR).at(i_cent).at(i_jet)->Scale(1./area);
	//				h_ChPS_raw.at(i_dR).at(i_cent).at(i_jet)->Sumw2();
	//
	//				name = Form("h_ChPS_unf_cent%i_dR%i_jetpt%i",i_cent, i_dR, i_jet);
	//				h_ChPS_unf.at(i_dR).at(i_cent).at(i_jet) = (TH1*)h_unfolded.at(i_cent).at(i_dR)->ProjectionX(name.c_str(), i_jet+1,i_jet+1);
	//				h_ChPS_unf.at(i_dR).at(i_cent).at(i_jet)->Scale(1.,"width");
	//				h_ChPS_unf.at(i_dR).at(i_cent).at(i_jet)->Scale(1./area);
	//				h_ChPS_unf.at(i_dR).at(i_cent).at(i_jet)->Sumw2();
	//
	//				name = Form("h_ChPS_tru_cent%i_dR%i_jetpt%i",i_cent, i_dR, i_jet);
	//				h_ChPS_tru.at(i_dR).at(i_cent).at(i_jet) = (TH1*)h_true.at(i_cent).at(i_dR)->ProjectionX(name.c_str(), i_jet+1, i_jet+1);
	//				h_ChPS_tru.at(i_dR).at(i_cent).at(i_jet)->Scale(1.,"width");
	//				h_ChPS_tru.at(i_dR).at(i_cent).at(i_jet)->Scale(1./area);
	//				h_ChPS_tru.at(i_dR).at(i_cent).at(i_jet)->Sumw2();
	//
	//				name = Form("h_ChPS_UE_cent%i_dR%i_jetpt%i",i_cent, i_dR, i_jet);
	//				h_ChPS_UE.at(i_dR).at(i_cent).at(i_jet) = (TH1*)h_UE.at(i_cent).at(i_dR)->ProjectionX(name.c_str(), i_jet+1, i_jet+1);
	//				h_ChPS_UE.at(i_dR).at(i_cent).at(i_jet)->Scale(1.,"width");
	//				h_ChPS_UE.at(i_dR).at(i_cent).at(i_jet)->Scale(1./area);
	//				h_ChPS_UE.at(i_dR).at(i_cent).at(i_jet)->Sumw2();
	//
	//				name = Form("h_ChPS_closure_cent%i_dR%i_jetpt%i",i_cent, i_dR, i_jet);
	//				h_ChPS_closure.at(i_dR).at(i_cent).at(i_jet) = (TH1*)h_ChPS_unf.at(i_dR).at(i_cent).at(i_jet)->Clone(name.c_str());
	//				h_ChPS_closure.at(i_dR).at(i_cent).at(i_jet)->Divide(h_ChPS_tru.at(i_dR).at(i_cent).at(i_jet));
	//
	//				name = Form("h_ChPS_closure_diff_cent%i_dR%i_jetpt%i",i_cent, i_dR, i_jet);
	//				h_ChPS_closure_diff.at(i_dR).at(i_cent).at(i_jet) = (TH1*)h_ChPS_unf.at(i_dR).at(i_cent).at(i_jet)->Clone(name.c_str());
	//				h_ChPS_closure_diff.at(i_dR).at(i_cent).at(i_jet)->Add(h_ChPS_tru.at(i_dR).at(i_cent).at(i_jet), -1);
	//
	//				name = Form("h_ChPS_ratio_cent%i_dR%i_jetpt%i",i_cent, i_dR, i_jet);
	//				h_ChPS_ratio.at(i_dR).at(i_cent).at(i_jet) = (TH1*)h_ChPS_unf.at(i_dR).at(i_cent).at(i_jet)->Clone(name.c_str());
	//				h_ChPS_ratio.at(i_dR).at(i_cent).at(i_jet)->Divide(h_ChPS_raw.at(i_dR).at(i_cent).at(i_jet));
	//
	//				name = Form("h_ChPS_UE_SB_cent%i_dR%i_jetpt%i",i_cent, i_dR, i_jet);
	//				h_ChPS_UE_SB.at(i_dR).at(i_cent).at(i_jet) = (TH1*)h_ChPS_UE.at(i_dR).at(i_cent).at(i_jet)->Clone(name.c_str());
	//				h_ChPS_UE_SB.at(i_dR).at(i_cent).at(i_jet)->Divide(h_ChPS_raw.at(i_dR).at(i_cent).at(i_jet));
	//
	//				if (doing_injet)
	//				{
	//					name = Form("h_ChPS_raw_injet_cent%i_jetpt%i",i_cent, i_jet);
	//					h_ChPS_raw_injet.at(i_cent).at(i_jet) = (TH1*)h_raw_injet.at(i_cent)->ProjectionX(name.c_str(), i_jet+1,i_jet+1);
	//					h_ChPS_raw_injet.at(i_cent).at(i_jet)->Scale(1.,"width");
	//					h_ChPS_raw_injet.at(i_cent).at(i_jet)->Scale(1./area_injet);
	//					h_ChPS_raw_injet.at(i_cent).at(i_jet)->Sumw2();
	//
	//					name = Form("h_ChPS_UE_injet_cent%i_jetpt%i",i_cent, i_jet);
	//					h_ChPS_UE_injet.at(i_cent).at(i_jet) = (TH1*)h_UE_injet.at(i_cent)->ProjectionX(name.c_str(), i_jet+1, i_jet+1);
	//					h_ChPS_UE_injet.at(i_cent).at(i_jet)->Scale(1.,"width");
	//					h_ChPS_UE_injet.at(i_cent).at(i_jet)->Scale(1./area_injet);
	//					h_ChPS_UE_injet.at(i_cent).at(i_jet)->Sumw2();
	//
	//
	//					name = Form("h_ChPS_UE_SB_injet_cent%i_jetpt%i",i_cent, i_jet);
	//					h_ChPS_UE_SB_injet.at(i_cent).at(i_jet) = (TH1*)h_ChPS_UE_injet.at(i_cent).at(i_jet)->Clone(name.c_str());
	//					h_ChPS_UE_SB_injet.at(i_cent).at(i_jet)->Divide(h_ChPS_raw_injet.at(i_cent).at(i_jet));
	//				}
	//
	//				f_output->cd();
	//
	//				name = Form("h_ChPS_raw_cent%i_dR%i_jetpt%i",i_cent, i_dR, i_jet);
	//				h_ChPS_raw.at(i_dR).at(i_cent).at(i_jet)->SetName(name.c_str());
	//				h_ChPS_raw.at(i_dR).at(i_cent).at(i_jet)->Write(name.c_str());
	//
	//				name = Form("h_ChPS_tru_cent%i_dR%i_jetpt%i",i_cent, i_dR, i_jet);
	//				h_ChPS_tru.at(i_dR).at(i_cent).at(i_jet)->SetName(name.c_str());
	//				h_ChPS_tru.at(i_dR).at(i_cent).at(i_jet)->Write(name.c_str());
	//
	//				name = Form("h_ChPS_unf_cent%i_dR%i_jetpt%i",i_cent, i_dR, i_jet);
	//				h_ChPS_unf.at(i_dR).at(i_cent).at(i_jet)->SetName(name.c_str());
	//				h_ChPS_unf.at(i_dR).at(i_cent).at(i_jet)->Write(name.c_str());
	//
	//				name = Form("h_ChPS_closure_cent%i_dR%i_jetpt%i",i_cent, i_dR, i_jet);
	//				h_ChPS_closure.at(i_dR).at(i_cent).at(i_jet)->SetName(name.c_str());
	//				h_ChPS_closure.at(i_dR).at(i_cent).at(i_jet)->Write(name.c_str());
	//
	//				name = Form("h_ChPS_closure_diff_cent%i_dR%i_jetpt%i",i_cent, i_dR, i_jet);
	//				h_ChPS_closure_diff.at(i_dR).at(i_cent).at(i_jet)->SetName(name.c_str());
	//				h_ChPS_closure_diff.at(i_dR).at(i_cent).at(i_jet)->Write(name.c_str());
	//
	//				name = Form("h_ChPS_ratio_cent%i_dR%i_jetpt%i",i_cent, i_dR, i_jet);
	//				h_ChPS_ratio.at(i_dR).at(i_cent).at(i_jet)->SetName(name.c_str());
	//				h_ChPS_ratio.at(i_dR).at(i_cent).at(i_jet)->Write(name.c_str());
	//
	//				name = Form("h_ChPS_UE_cent%i_dR%i_jetpt%i",i_cent, i_dR, i_jet);
	//				h_ChPS_UE.at(i_dR).at(i_cent).at(i_jet)->SetName(name.c_str());
	//				h_ChPS_UE.at(i_dR).at(i_cent).at(i_jet)->Write(name.c_str());
	//
	//				name = Form("h_ChPS_UE_SB_cent%i_dR%i_jetpt%i",i_cent, i_dR, i_jet);
	//				h_ChPS_UE_SB.at(i_dR).at(i_cent).at(i_jet)->SetName(name.c_str());
	//				h_ChPS_UE_SB.at(i_dR).at(i_cent).at(i_jet)->Write(name.c_str());
	//
	//				if (doing_injet)
	//				{
	//					name = Form("h_ChPS_UE_SB_injet_cent%i_dR%i_jetpt%i",i_cent, i_dR, i_jet);
	//					h_ChPS_UE_injet.at(i_cent).at(i_jet)->SetName(name.c_str());
	//					h_ChPS_UE_injet.at(i_cent).at(i_jet)->Write(name.c_str());
	//
	//					name = Form("h_ChPS_UE_SB_injet_cent%i_dR%i_jetpt%i",i_cent, i_dR, i_jet);
	//					h_ChPS_UE_SB_injet.at(i_cent).at(i_jet)->SetName(name.c_str());
	//					h_ChPS_UE_SB_injet.at(i_cent).at(i_jet)->Write(name.c_str());
	//				}
	//
	//				if (i_jet < jet_pt_start || i_jet > jet_pt_end) continue;
	//
	//				SetHStyle(h_ChPS_raw.at(i_dR).at(i_cent).at(i_jet),style);
	//				SetHStyle(h_ChPS_unf.at(i_dR).at(i_cent).at(i_jet),style);
	//				SetHStyle(h_ChPS_tru.at(i_dR).at(i_cent).at(i_jet),style);
	//				SetHStyle(h_ChPS_closure.at(i_dR).at(i_cent).at(i_jet),style);
	//				SetHStyle(h_ChPS_closure_diff.at(i_dR).at(i_cent).at(i_jet),style);
	//				SetHStyle(h_ChPS_ratio.at(i_dR).at(i_cent).at(i_jet),style);
	//				SetHStyle(h_ChPS_UE.at(i_dR).at(i_cent).at(i_jet),style);
	//				SetHStyle(h_ChPS_UE_SB.at(i_dR).at(i_cent).at(i_jet),style);
	//				if (doing_injet) SetHStyle(h_ChPS_UE_injet.at(i_cent).at(i_jet),style);
	//				if (doing_injet) SetHStyle(h_ChPS_UE_SB_injet.at(i_cent).at(i_jet),style);
	//
	//				smallify(h_ChPS_raw.at(i_dR).at(i_cent).at(i_jet));
	//				smallify(h_ChPS_unf.at(i_dR).at(i_cent).at(i_jet));
	//				smallify(h_ChPS_tru.at(i_dR).at(i_cent).at(i_jet));
	//				smallify(h_ChPS_closure.at(i_dR).at(i_cent).at(i_jet));
	//				smallify(h_ChPS_closure_diff.at(i_dR).at(i_cent).at(i_jet));
	//				smallify(h_ChPS_ratio.at(i_dR).at(i_cent).at(i_jet));
	//				smallify(h_ChPS_UE.at(i_dR).at(i_cent).at(i_jet));
	//				smallify(h_ChPS_UE_SB.at(i_dR).at(i_cent).at(i_jet));
	//				if (doing_injet) smallify(h_ChPS_UE_injet.at(i_cent).at(i_jet));
	//				if (doing_injet) smallify(h_ChPS_UE_SB_injet.at(i_cent).at(i_jet));
	//
	//
	//				double range_pt_lo = 1;
	//				double range_pt_hi = 158.;
	//				if (dR_lo >= 0.5) range_pt_hi = 60.;
	//				if (dR_lo >= 0.7) range_pt_hi = 20.;
	//
	//				h_ChPS_raw.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(range_pt_lo, range_pt_hi);
	//				h_ChPS_unf.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(range_pt_lo, range_pt_hi);
	//				h_ChPS_tru.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(range_pt_lo, range_pt_hi);
	//				h_ChPS_ratio.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(range_pt_lo, range_pt_hi);
	//				h_ChPS_closure.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(range_pt_lo, range_pt_hi);
	//				h_ChPS_closure_diff.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(range_pt_lo, range_pt_hi);
	//				h_ChPS_UE.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(range_pt_lo, range_pt_hi);
	//				h_ChPS_UE_injet.at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(range_pt_lo, range_pt_hi);
	//				h_ChPS_UE_SB.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(range_pt_lo, range_pt_hi);
	//				h_ChPS_UE_SB_injet.at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(range_pt_lo, range_pt_hi);
	//
	//				// *** C1 *** //
	//				h_ChPS_raw.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]");
	//				h_ChPS_raw.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetTitle("#frac{1}{Area N_{Jets}} #frac{dN}{dp_{T}}");
	//				h_ChPS_tru.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]");
	//				h_ChPS_tru.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetTitle("#frac{1}{Area N_{Jets}} #frac{dN}{dp_{T}}");
	//
	//
	//				h_ChPS_raw.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetNdivisions(504);
	//				h_ChPS_unf.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetNdivisions(504);
	//				h_ChPS_tru.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetNdivisions(504);
	//				h_ChPS_closure.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetNdivisions(504);
	//				h_ChPS_closure_diff.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetNdivisions(504);
	//
	//				c1->cd(i_cent+1)->cd(1);
	//				gPad->SetPad(0,0.51,1,.95);
	//				gPad->SetBottomMargin(0);
	//				gPad->SetRightMargin(0);
	//				if (i_jet == jet_pt_start) h_ChPS_raw.at(i_dR).at(i_cent).at(i_jet)->Draw();
	//				else h_ChPS_raw.at(i_dR).at(i_cent).at(i_jet)->Draw("same");
	//				gPad->SetLogx();
	//				gPad->SetLogy();
	//
	//				c1->cd(i_cent+1)->cd(2);
	//				gPad->SetPad(0, 0.0, 1, 0.51);
	//				gPad->SetTopMargin(0);
	//				gPad->SetRightMargin(0);
	//				if (i_jet == jet_pt_start) h_ChPS_tru.at(i_dR).at(i_cent).at(i_jet)->Draw();
	//				else h_ChPS_tru.at(i_dR).at(i_cent).at(i_jet)->Draw("same");
	//				gPad->SetLogx();
	//				gPad->SetLogy();
	//				// ***
	//
	//				// *** C8 *** //
	//				h_ChPS_unf.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]");
	//				h_ChPS_unf.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetTitle("#frac{1}{Area N_{Jets}} #frac{dN}{dp_{T}}");
	//
	//				h_ChPS_unf.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetNdivisions(504);
	//
	//				c8->cd(i_cent+1);
	//				if (i_jet == jet_pt_start) h_ChPS_unf.at(i_dR).at(i_cent).at(i_jet)->Draw();
	//				else h_ChPS_unf.at(i_dR).at(i_cent).at(i_jet)->Draw("same");
	//				gPad->SetLogx();
	//				gPad->SetLogy();
	//				// ***
	//
	//
	//				// *** C2 *** //
	//				h_ChPS_closure.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]");
	//				h_ChPS_closure.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetRangeUser(0.45,1.55);
	//
	//				h_ChPS_closure_diff.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]");
	//				h_ChPS_closure_diff.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetRangeUser(-0.15,0.15);
	//
	//				name = "Closure"; //default title
	//				if (ChPS_raw_type == "ChPS_raw_tt") name = "Eff. Corr. Closure";
	//				if (ChPS_raw_type == "ChPS_raw_0") name = "Unf/Truth";
	//				h_ChPS_closure.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetTitle(name.c_str());
	//				h_ChPS_closure_diff.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetTitle("Unf - Truth");
	//
	//				c2->cd(i_cent+1);
	//				c2->cd(i_cent+1)->cd(1);
	//				gPad->SetPad(0,0.51,1,.95);
	//				gPad->SetBottomMargin(0);
	//				gPad->SetRightMargin(0);
	//				if (i_jet == jet_pt_start) h_ChPS_closure.at(i_dR).at(i_cent).at(i_jet)->Draw();
	//				else h_ChPS_closure.at(i_dR).at(i_cent).at(i_jet)->Draw("same");
	//				line->SetLineStyle(3);
	//				line->DrawLine(1,1.02, 150, 1.02);
	//				line->DrawLine(1,0.98, 150, 0.98);
	//
	//				line->SetLineStyle(2);
	//				line->DrawLine(1,1, 150, 1);
	//				gPad->SetLogx();
	//				gPad->SetLogy(0);
	//
	//
	//
	//				c2->cd(i_cent+1)->cd(2);
	//				gPad->SetPad(0, 0.0, 1, 0.51);
	//				gPad->SetTopMargin(0);
	//				gPad->SetRightMargin(0);
	//				if (i_jet == jet_pt_start) h_ChPS_closure_diff.at(i_dR).at(i_cent).at(i_jet)->Draw();
	//				else h_ChPS_closure_diff.at(i_dR).at(i_cent).at(i_jet)->Draw("same");
	//				line->DrawLine(1,0,158,0);
	//				gPad->SetLogx();
	//				gPad->SetLogy(0);
	//				// ***
	//
	//
	//				// *** C7 *** //
	//				h_ChPS_ratio.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]");
	//				h_ChPS_ratio.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetRangeUser(0,2);
	//				name = "(Unf+BbB)/(Raw-UE)"; //default title
	//				if (ChPS_raw_type == "ChPS_raw_tt") name = "Eff. Corr. Closure";
	//				if (ChPS_raw_type == "ChPS_raw_0" && do_UE_subtr && do_unfolding && !do_BbB) name = "Unf/(Raw-UE)";
	//				h_ChPS_ratio.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetTitle(name.c_str());
	//
	//				c7->cd(i_cent+1);
	//				if (i_jet == jet_pt_start) h_ChPS_ratio.at(i_dR).at(i_cent).at(i_jet)->Draw();
	//				else h_ChPS_ratio.at(i_dR).at(i_cent).at(i_jet)->Draw("same");
	//				gPad->SetLogx();
	//				gPad->SetLogy(0);
	//				// ***
	//
	//				// *** C3 *** //
	//				h_ChPS_UE.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]");
	//				h_ChPS_UE.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetTitle("#frac{1}{N_{jets}^{raw} Area} #frac{dN_{UE}}{dp_{T}}");
	//				h_ChPS_UE.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetTitleFont(43);
	//				h_ChPS_UE.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetTitleSize(12);
	//				h_ChPS_UE.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetTitleOffset(3);
	//
	//				max = h_ChPS_UE.at(i_dR).at(i_cent).at(i_jet)->GetMaximum();
	//				min = h_ChPS_UE.at(i_dR).at(i_cent).at(i_jet)->GetMinimum();
	//				h_ChPS_UE.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetRangeUser(0.5*min,1.5*max);
	//				h_ChPS_UE.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetNdivisions(504);
	//
	//				c3->cd(i_cent+1);
	//				if (i_jet == jet_pt_start) h_ChPS_UE.at(i_dR).at(i_cent).at(i_jet)->Draw();
	//				else h_ChPS_UE.at(i_dR).at(i_cent).at(i_jet)->Draw("same");
	//				gPad->SetLogx();
	//				gPad->SetLogy(0);
	//				// ***
	//
	//				// *** C5 *** //
	//				h_ChPS_UE_SB.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]");
	//				h_ChPS_UE_SB.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetTitle("#frac{dN_{UE}}{dN_{ChPS}^{Subtr}}");
	//				h_ChPS_UE_SB.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetTitleFont(43);
	//				h_ChPS_UE_SB.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetTitleSize(12);
	//				h_ChPS_UE_SB.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetTitleOffset(3);
	//				h_ChPS_UE_SB.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetNdivisions(504);
	//
	//				c5->cd(i_cent+1);
	//				if (i_jet == jet_pt_start) h_ChPS_UE_SB.at(i_dR).at(i_cent).at(i_jet)->Draw();
	//				else h_ChPS_UE_SB.at(i_dR).at(i_cent).at(i_jet)->Draw("same");
	//				gPad->SetLogx();
	//				gPad->SetLogy(0);
	//				// ***
	//
	//				// *** C4 *** //
	//				if (doing_injet)
	//				{
	//					h_ChPS_UE_injet.at(i_cent).at(i_jet)->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]");
	//					h_ChPS_UE_injet.at(i_cent).at(i_jet)->GetYaxis()->SetTitle("#frac{1}{N_{jets}^{raw} Area^{Jet}} #frac{dN_{UE}}{dp_{T}} (In jet)");
	//					h_ChPS_UE_injet.at(i_cent).at(i_jet)->GetYaxis()->SetTitleFont(43);
	//					h_ChPS_UE_injet.at(i_cent).at(i_jet)->GetYaxis()->SetTitleSize(12);
	//					h_ChPS_UE_injet.at(i_cent).at(i_jet)->GetYaxis()->SetTitleOffset(3);
	//					max = h_ChPS_UE_injet.at(i_cent).at(i_jet)->GetMaximum();
	//					min = h_ChPS_UE_injet.at(i_cent).at(i_jet)->GetMinimum();
	//					h_ChPS_UE_injet.at(i_cent).at(i_jet)->GetYaxis()->SetRangeUser(0,1.5*max);
	//					h_ChPS_UE_injet.at(i_cent).at(i_jet)->GetYaxis()->SetNdivisions(504);
	//
	//					h_ChPS_UE_SB_injet.at(i_cent).at(i_jet)->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]");
	//					h_ChPS_UE_SB_injet.at(i_cent).at(i_jet)->GetYaxis()->SetTitle("#frac{dN_{UE}}{dN_{ChPS}^{Subtr}} (In jet)");
	//					h_ChPS_UE_SB_injet.at(i_cent).at(i_jet)->GetYaxis()->SetTitleFont(43);
	//					h_ChPS_UE_SB_injet.at(i_cent).at(i_jet)->GetYaxis()->SetTitleSize(12);
	//					h_ChPS_UE_SB_injet.at(i_cent).at(i_jet)->GetYaxis()->SetTitleOffset(3);
	//					max = h_ChPS_UE_SB_injet.at(i_cent).at(i_jet)->GetMaximum();
	//					min = h_ChPS_UE_SB_injet.at(i_cent).at(i_jet)->GetMinimum();
	//					h_ChPS_UE_SB_injet.at(i_cent).at(i_jet)->GetYaxis()->SetRangeUser(0,1.5*max);
	//					h_ChPS_UE_SB_injet.at(i_cent).at(i_jet)->GetYaxis()->SetNdivisions(504);
	//
	//
	//					c4->cd(i_cent+1);
	//					if (i_jet == jet_pt_start) h_ChPS_UE_injet.at(i_cent).at(i_jet)->Draw();
	//					else h_ChPS_UE_injet.at(i_cent).at(i_jet)->Draw("same");
	//					gPad->SetLogx();
	//					gPad->SetLogy(0);
	//
	//					c6->cd(i_cent+1);
	//					if (i_jet == jet_pt_start) h_ChPS_UE_SB_injet.at(i_cent).at(i_jet)->Draw();
	//					else h_ChPS_UE_SB_injet.at(i_cent).at(i_jet)->Draw("same");
	//					gPad->SetLogx();
	//					gPad->SetLogy(0);
	//
	//				}
	//				// ***
	//
	//
	//
	//				if (i_cent == 0)
	//				{
	//					double pt_lo = jetpT_binning->GetBinLowEdge(i_jet+1);
	//					double pt_hi = jetpT_binning->GetBinLowEdge(i_jet+2);
	//
	//					name = Form("%2.0f < p_{T}^{Jet} < %2.0f GeV",pt_lo, pt_hi);
	//					legend->AddEntry(h_ChPS_tru.at(i_dR).at(i_cent).at(i_jet),name.c_str(),"lp");
	//				}
	//
	//				style++;
	//			}
	//
	//			// *****
	//			c1->cd(i_cent+1);
	//			ltx->SetTextAlign(32);
	//			ltx->DrawLatexNDC(1,0.98,num_to_cent(31, i_cent).c_str());
	//			ltx->SetTextAlign(32);
	//			name = "Raw (No Subtr)";
	//			if (do_UE_subtr) name = "Raw+Subtr";
	//			if (do_UE_subtr && do_unfolding) name = "Raw+Subtr+Unf";
	//			if (do_UE_subtr && do_unfolding && do_BbB) name = "Raw+Subtr+Unf+BbB";
	//			ltx->DrawLatexNDC(0.93,0.89,name.c_str());
	//			ltx->DrawLatexNDC(0.93,0.48,"Truth");
	//			ltx->SetTextAlign(12);
	//			ltx->DrawLatexNDC(0.19,0.98,Form("%4.2f < dR < %4.2f", dR_lo, dR_hi));
	//			// *****
	//
	//			// *****
	//			c8->cd(i_cent+1);
	//			ltx->SetTextAlign(32);
	//			ltx->DrawLatexNDC(1,0.98,num_to_cent(31, i_cent).c_str());
	//			ltx->SetTextAlign(32);
	//
	//			name = "Raw (No Subtr)";
	//			if (do_UE_subtr) name = "Raw+Subtr";
	//			if (do_UE_subtr && do_unfolding) name = "Raw+Subtr+Unf";
	//			if (do_UE_subtr && do_unfolding && do_BbB) name = "Raw+Subtr+Unf+BbB";
	//			ltx->DrawLatexNDC(0.93,0.89,name.c_str());
	//			ltx->SetTextAlign(12);
	//			ltx->DrawLatexNDC(0.19,0.98,Form("%4.2f < dR < %4.2f", dR_lo, dR_hi));
	//			// *****
	//
	//
	//			// *****
	//			c2->cd(i_cent+1);
	//			ltx->SetTextAlign(32);
	//			ltx->DrawLatexNDC(0.92,0.98,num_to_cent(31, i_cent).c_str());
	//
	//			ltx->SetTextAlign(12);
	//			ltx->DrawLatexNDC(0.19,0.98,Form("%4.2f < dR < %4.2f", dR_lo, dR_hi));
	//			// *****
	//
	//			// *****
	//			c7->cd(i_cent+1);
	//			ltx->SetTextAlign(32);
	//			ltx->DrawLatexNDC(0.92,0.98,num_to_cent(31, i_cent).c_str());
	//
	//			line->SetLineStyle(3);
	//			line->DrawLine(1,1.02, 150, 1.02);
	//			line->DrawLine(1,0.98, 150, 0.98);
	//
	//			line->SetLineStyle(2);
	//			line->DrawLine(1,1, 150, 1);
	//
	//			ltx->SetTextAlign(12);
	//			ltx->DrawLatexNDC(0.19,0.98,Form("%4.2f < dR < %4.2f", dR_lo, dR_hi));
	//			// *****
	//
	//
	//			// *****
	//			c3->cd(i_cent+1);
	//			ltx->SetTextAlign(32);
	//			ltx->DrawLatexNDC(0.92,0.98,num_to_cent(31, i_cent).c_str());
	//
	//			ltx->SetTextAlign(12);
	//			ltx->DrawLatexNDC(0.19,0.92,Form("%4.2f < dR < %4.2f", dR_lo, dR_hi));
	//			// *****
	//
	//
	//			// *****
	//			c5->cd(i_cent+1);
	//			ltx->SetTextAlign(32);
	//			ltx->DrawLatexNDC(0.92,0.98,num_to_cent(31, i_cent).c_str());
	//
	//			ltx->SetTextAlign(12);
	//			ltx->DrawLatexNDC(0.19,0.92,Form("%4.2f < dR < %4.2f", dR_lo, dR_hi));
	//			// *****
	//
	//
	//			// *****
	//			if (doing_injet)
	//			{
	//				c4->cd(i_cent+1);
	//				ltx->SetTextAlign(32);
	//				ltx->DrawLatexNDC(0.92,0.98,num_to_cent(31, i_cent).c_str());
	//				ltx->SetTextAlign(12);
	//				ltx->DrawLatexNDC(0.19,0.98,Form("injet"));
	//
	//				c6->cd(i_cent+1);
	//				ltx->SetTextAlign(32);
	//				ltx->DrawLatexNDC(0.92,0.98,num_to_cent(31, i_cent).c_str());
	//				ltx->SetTextAlign(12);
	//				ltx->DrawLatexNDC(0.19,0.98,Form("injet"));
	//			}
	//			// *****
	//
	//
	//		} //end cent loop
	//
	//
	//		c1->cd(1);
	//		legend->SetX1NDC(0.20);
	//		legend->SetY1NDC(0.20);
	//		legend->SetX2NDC(0.48);
	//		legend->SetY2NDC(0.40);
	//		legend->Draw();
	//		if (i_dR == 0) name = "(";
	//		else if (i_dR == N_dR - 1) name = ")";
	//		else name = "";
	//		c1->Print(Form("%s_spectra_%s_%s_%s.pdf%s",ChPS_raw_type.c_str(), dataset.c_str(), dataset_type.c_str(), evol.c_str(), name.c_str()),Form("Title: dR%i",i_dR));
	//		c1->Clear();
	//
	//
	//		c8->cd(1);
	//		legend->SetX1NDC(0.20);
	//		legend->SetY1NDC(0.20);
	//		legend->SetX2NDC(0.48);
	//		legend->SetY2NDC(0.40);
	//		legend->Draw();
	//		if (i_dR == 0) name = "(";
	//		else if (i_dR == N_dR - 1) name = ")";
	//		else name = "";
	//		c8->Print(Form("final_ChPS_%s_%s_%s.pdf%s",dataset.c_str(), dataset_type.c_str(), evol.c_str(), name.c_str()),Form("Title: dR%i",i_dR));
	//		c8->Clear();
	//
	//
	//		c2->cd(1)->cd(1);
	//		legend->SetX1NDC(0.50);
	//		legend->SetX2NDC(0.80);
	//		legend->SetY1NDC(0.50);
	//		legend->SetY2NDC(0.85);
	//		legend->Draw();
	//		c2->cd(1);
	//		ltx->SetTextAlign(12);
	//		if (i_dR == 0) name = "(";
	//		else if (i_dR == N_dR - 1) name = ")";
	//		else name = "";
	//		c2->Print(Form("%s_closure_%s_%s_%s.pdf%s",ChPS_raw_type.c_str(),dataset.c_str(), dataset_type.c_str(), evol.c_str(), name.c_str()),Form("Title: dR%i",i_dR));
	//		c2->Print("tmp.root");
	//		c2->Clear();
	//
	//		c7->cd(1);
	//		legend->SetX1NDC(0.20);
	//		legend->SetY1NDC(0.20);
	//		legend->SetX2NDC(0.40);
	//		legend->SetY2NDC(0.40);
	//		legend->Draw();
	//		c7->cd(1);
	//		ltx->SetTextAlign(12);
	//		if (i_dR == 0) name = "(";
	//		else if (i_dR == N_dR - 1) name = ")";
	//		else name = "";
	//		c7->Print(Form("%s_ratio_%s_%s_%s.pdf%s",ChPS_raw_type.c_str(), dataset.c_str(), dataset_type.c_str(), evol.c_str(),name.c_str()),Form("Title: dR%i",i_dR));
	//		c7->Clear();
	//
	//
	//		c3->cd(1);
	//		legend->SetX1NDC(0.48);
	//		legend->SetX2NDC(0.80);
	//		legend->SetY1NDC(0.60);
	//		legend->SetY2NDC(0.80);
	//		legend->Draw();
	//		c3->cd(1);
	//		ltx->SetTextAlign(12);
	//		if (i_dR == 0) name = "(";
	//		else if (i_dR == N_dR - 1) name = ")";
	//		else name = "";
	//		c3->Print(Form("UE_%s_%s_%s.pdf%s", dataset.c_str(), dataset_type.c_str(), evol.c_str(),name.c_str()),Form("Title: dR%i",i_dR));
	//		c3->Clear();
	//
	//
	//		c5->cd(1);
	//		legend->SetX1NDC(0.48);
	//		legend->SetX2NDC(0.80);
	//		legend->SetY1NDC(0.60);
	//		legend->SetY2NDC(0.80);
	//		legend->Draw();
	//		c5->cd(1);
	//		ltx->SetTextAlign(12);
	//		if (i_dR == 0) name = "(";
	//		else if (i_dR == N_dR - 1) name = ")";
	//		else name = "";
	//		c5->Print(Form("UE_SB_%s_%s_%s.pdf%s", dataset.c_str(), dataset_type.c_str(), evol.c_str(),name.c_str()),Form("Title: dR%i",i_dR));
	//		c5->Clear();
	//
	//		if (doing_injet)
	//		{
	//			c4->cd(1);
	//			legend->SetX1NDC(0.48);
	//			legend->SetX2NDC(0.80);
	//			legend->SetY1NDC(0.60);
	//			legend->SetY2NDC(0.80);
	//			legend->Draw();
	//			c4->cd(1);
	//			ltx->SetTextAlign(12);
	//
	//			name = Form("UE_injet.pdf");
	//			c4->Print(name.c_str(),Form("Title: inJet"));
	//			c4->Clear();
	//
	//			c6->cd(1);
	//			legend->SetX1NDC(0.48);
	//			legend->SetX2NDC(0.80);
	//			legend->SetY1NDC(0.60);
	//			legend->SetY2NDC(0.80);
	//			legend->Draw();
	//			c6->cd(1);
	//			ltx->SetTextAlign(12);
	//
	//			name = Form("UE_SB_injet.pdf");
	//			c6->Print(name.c_str(),Form("Title: inJet"));
	//			c6->Clear();
	//		}
	//
	//
	//
	//		doing_injet = false;
	//		legend->Clear();
	//
	//	} //end DR loop
*/


}
