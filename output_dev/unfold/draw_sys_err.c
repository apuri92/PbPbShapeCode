#include "../functions/global_variables.h"
#include "draw_functions.c"
#include "TEnv.h"

void draw_sys_err(string config_file = "sys_config.cfg")
{
	cout << "######### DOING draw Systematics #########" << endl;
	SetAtlasStyle();
	gErrorIgnoreLevel = 3001;
	string name = "";

	//	##############	Reading config	##############"
	TEnv *m_config = new TEnv();
	m_config->ReadFile(config_file.c_str(), EEnvLevel(1));

	std::string mode = "RDpT"; mode = m_config->GetValue("mode", mode.c_str());
	int isMC = 1; isMC = m_config->GetValue("isMC", isMC);

	std::string dataset_type = "PbPb"; dataset_type = m_config->GetValue("dataset_type", dataset_type.c_str());
	int centrality_scheme = 31; centrality_scheme = m_config->GetValue("centrality_scheme", centrality_scheme);
	int verbose = 0; verbose = m_config->GetValue("verbose", verbose);

	std::string did = "data";
	if (isMC) did = "MC";

	if (verbose) m_config->Print();
	//	##############	Config done	##############"

	bool doSmall = true;
	doSmall = false;

	if (mode == "RDpT") dataset_type = "";
	else dataset_type = Form("_%s", dataset_type.c_str());

	TFile* sys_file = new TFile(Form("output_pdf_nominal/root/final_%s_sys_%s%s.root", mode.c_str(), did.c_str(), dataset_type.c_str()));

	cout << "Using files:" << endl;
	cout << sys_file->GetName() << endl;

	TAxis* dR_binning = (TAxis*)sys_file->Get("dR_binning");
	TAxis* jetpT_binning = (TAxis*)sys_file->Get("jetpT_binning");
	TAxis* trkpT_binning = (TAxis*)sys_file->Get("trkpT_binning");

	int N_dR = dR_binning->GetNbins();
	int N_jetpt = jetpT_binning->GetNbins();
	int N_trkpt = trkpT_binning->GetNbins();


	vector<TFile*> sys_files;
	vector<string> sys_names;
	vector<string> combined_sys_names;

	combined_sys_names.push_back("JER");
	combined_sys_names.push_back("JES");
	combined_sys_names.push_back("UE");
	combined_sys_names.push_back("Tracking");
	combined_sys_names.push_back("Unfolding");
	combined_sys_names.push_back("MCNonClosure");

	vector<vector<vector<TH1*>>> h_total_sys_p (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_total_sys_n (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));

	vector<vector<vector<vector<TH1*>>>> h_comb_sys_p (combined_sys_names.size(), vector<vector<vector<TH1*>>> (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt))));
	vector<vector<vector<vector<TH1*>>>> h_comb_sys_n (combined_sys_names.size(), vector<vector<vector<TH1*>>> (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt))));


	int jet_pt_start = 7;
	int jet_pt_end = 11;
	int trk_pt_start = 3;
	int trk_pt_end = 9;

	string rdptr_label = "#it{R}_{ #it{D} (#it{p}_{T}, #it{r})}";
	if (mode == "ChPS") rdptr_label = "#it{D} (#it{p}_{T}, #it{r})";
	string r_label = "#it{r}";

	TCanvas *c_sys = new TCanvas("c_sys","c_sys", 800, 600);
	TLegend *legend_sys = new TLegend(0.18, 0.18, 0.40, 0.40, "","brNDC");
	legend_sys->SetTextFont(43);
	legend_sys->SetBorderSize(0);
	legend_sys->SetTextSize(20);

	TLatex *ltx = new TLatex();
	ltx->SetTextFont(43);
	ltx->SetTextSize(12);

	TLine *line = new TLine();
	line->SetLineColor(kBlack);

	int trk_select1 = 3;
	int trk_select2 = 6;

	for (int i_cent = 0; i_cent < n_cent_cuts; i_cent++)
	{
		if ((dataset_type == "_PbPb" || mode == "RDpT") && i_cent == 6) continue;
		if (dataset_type == "_pp" && i_cent < 6) continue;
		string centrality = num_to_cent(31,i_cent);

		for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
		{
			for (int i_trk = trk_pt_start; i_trk < 11; i_trk++)
			{

				if ((i_jet != 7 && i_jet != 8 && i_jet != 9 & i_jet != 10) || (i_trk != trk_select1 && i_trk != trk_select2)) continue;

				if ((dataset_type == "_PbPb" || mode == "RDpT") && (i_cent != 0 && i_cent != 5)) continue;

				string jet_label = Form("%1.0f < #it{p}_{T}^{jet} < %1.0f GeV", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));
				string trk_label = Form("%1.1f < #it{p}_{T} < %1.1f GeV", trkpT_binning->GetBinLowEdge(i_trk+1), trkpT_binning->GetBinUpEdge(i_trk+1));

				//setup canvas
				c_sys->cd();
				c_sys->Clear();
				legend_sys->Clear();

				name = Form("h_%s_sys_trk%i_cent%i_jetpt%i_total_p",mode.c_str(), i_trk, i_cent, i_jet);
				h_total_sys_p[i_trk][i_cent][i_jet] = (TH1*)sys_file->Get(name.c_str());

				name = Form("h_%s_sys_trk%i_cent%i_jetpt%i_total_n",mode.c_str(), i_trk, i_cent, i_jet);
				h_total_sys_n[i_trk][i_cent][i_jet] = (TH1*)sys_file->Get(name.c_str());

				//set bin error so bar through point is non-zero
				for (int i = 0 ; i < N_dR; i++)
				{
					h_total_sys_p[i_trk][i_cent][i_jet]->SetBinError(i+1,0.00000001);
					h_total_sys_n[i_trk][i_cent][i_jet]->SetBinError(i+1,0.00000001);
				}

				SetHStyle_smallify(h_total_sys_p[i_trk][i_cent][i_jet],0, doSmall);
				SetHStyle_smallify(h_total_sys_n[i_trk][i_cent][i_jet],0, doSmall);
//				h_total_sys_p[i_trk][i_cent][i_jet]->SetMarkerStyle(20);
//				h_total_sys_n[i_trk][i_cent][i_jet]->SetMarkerStyle(20);
//				h_total_sys_n[i_trk][i_cent][i_jet]->SetMarkerSize(1);
//				h_total_sys_n[i_trk][i_cent][i_jet]->SetMarkerSize(1);

				h_total_sys_p[i_trk][i_cent][i_jet]->GetYaxis()->SetRangeUser(-0.34,0.34);
				h_total_sys_n[i_trk][i_cent][i_jet]->GetYaxis()->SetRangeUser(-0.34,0.34);
				if (dataset_type == "_pp")
				{
					h_total_sys_p[i_trk][i_cent][i_jet]->GetYaxis()->SetRangeUser(-0.2,0.2);
					h_total_sys_n[i_trk][i_cent][i_jet]->GetYaxis()->SetRangeUser(-0.2,0.2);
				}
				h_total_sys_p[i_trk][i_cent][i_jet]->GetXaxis()->SetRangeUser(0, 0.6);
				h_total_sys_n[i_trk][i_cent][i_jet]->GetXaxis()->SetRangeUser(0, 0.6);

				h_total_sys_p[i_trk][i_cent][i_jet]->SetLineWidth(2);
				h_total_sys_n[i_trk][i_cent][i_jet]->SetLineWidth(2);


				h_total_sys_p[i_trk][i_cent][i_jet]->GetYaxis()->SetTitle(Form("#delta%s [%%]",rdptr_label.c_str()));
				h_total_sys_n[i_trk][i_cent][i_jet]->GetYaxis()->SetTitle(Form("#delta%s [%%]",rdptr_label.c_str()));
				h_total_sys_p[i_trk][i_cent][i_jet]->GetXaxis()->SetTitle(r_label.c_str());
				h_total_sys_n[i_trk][i_cent][i_jet]->GetXaxis()->SetTitle(r_label.c_str());

				h_total_sys_p[i_trk][i_cent][i_jet]->GetYaxis()->SetTitleFont(43);
				h_total_sys_p[i_trk][i_cent][i_jet]->GetYaxis()->SetTitleSize(30);
				h_total_sys_p[i_trk][i_cent][i_jet]->GetXaxis()->SetTitleFont(43);
				h_total_sys_p[i_trk][i_cent][i_jet]->GetXaxis()->SetTitleSize(30);

				c_sys->cd();
				h_total_sys_p[i_trk][i_cent][i_jet]->Draw("lp ");
				h_total_sys_n[i_trk][i_cent][i_jet]->Draw("lp same");

				for (int i_comb_sys = 0; i_comb_sys < combined_sys_names.size(); i_comb_sys++)
				{

					name = Form("h_%s_sys_trk%i_cent%i_jetpt%i_%s_p",mode.c_str(), i_trk, i_cent, i_jet,combined_sys_names[i_comb_sys].c_str());
					h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet] = (TH1*)sys_file->Get(name.c_str());

					name = Form("h_%s_sys_trk%i_cent%i_jetpt%i_%s_n",mode.c_str(), i_trk, i_cent, i_jet,combined_sys_names[i_comb_sys].c_str());
					h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet] = (TH1*)sys_file->Get(name.c_str());

					//set bin error so bar through point is non-zero
					for (int i = 0 ; i < N_dR; i++)
					{
						h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->SetBinError(i+1,0.00000001);
						h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->SetBinError(i+1,0.00000001);
					}

					SetHStyle_smallify(h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet],i_comb_sys+1, doSmall);
					SetHStyle_smallify(h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet],i_comb_sys+1, doSmall);
//					h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->SetMarkerStyle(20);
//					h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->SetMarkerStyle(20);
//					h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->SetMarkerSize(1);
//					h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->SetMarkerSize(1);

					h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->SetLineWidth(2);
					h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->SetLineWidth(2);

					bool empty_hist_p = false;
					bool empty_hist_n = false;
					if (h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->GetMaximum() == 0 && h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->GetMinimum() == 0) empty_hist_p = true;
					if (h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->GetMaximum() == 0 && h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->GetMinimum() == 0) empty_hist_n = true;

					h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->SetLineStyle(i_comb_sys+2);
					h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->SetLineStyle(i_comb_sys+2);

					if (i_comb_sys == 3)
					{
						h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->SetLineStyle(8);
						h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->SetLineStyle(8);

					}
					if (i_comb_sys == 4)
					{
						h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->SetLineStyle(7);
						h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->SetLineStyle(7);
					}

					if (combined_sys_names[i_comb_sys] == "MCNonClosure")
					{
						legend_sys->AddEntry(h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet],"MC non-closure","lp");
					}
					else
					{
						legend_sys->AddEntry(h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet],combined_sys_names[i_comb_sys].c_str(),"lp");
					}
					h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->GetYaxis()->SetRangeUser(-0.3,0.3);
					h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->GetYaxis()->SetRangeUser(-0.3,0.3);
					h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->GetXaxis()->SetRangeUser(0, 0.6);
					h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->GetXaxis()->SetRangeUser(0, 0.6);
					h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->GetYaxis()->SetTitle(Form("#delta%s",rdptr_label.c_str()));
					h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->GetYaxis()->SetTitle(Form("#delta%s",rdptr_label.c_str()));
					h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->GetXaxis()->SetTitle(r_label.c_str());
					h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->GetXaxis()->SetTitle(r_label.c_str());

					c_sys->cd();
					if (!empty_hist_p) h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->Draw("lp same");
					if (!empty_hist_n) h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->Draw("lp same");
				}

				c_sys->cd();
				if (dataset_type == "_pp") c_sys->cd();
				legend_sys->AddEntry(h_total_sys_p[i_trk][i_cent][i_jet],"Total","lp");

				legend_sys->Draw();
				ltx->SetTextAlign(32);
				ltx->SetTextSize(18);

				ltx->SetTextSize(25);
				ltx->DrawLatexNDC(0.93, 0.90, Form("%s", trk_label.c_str()));
				ltx->DrawLatexNDC(0.93, 0.84, Form("%s", jet_label.c_str()));
				if (dataset_type != "_pp") ltx->DrawLatexNDC(0.93, 0.79, Form("%s", centrality.c_str()));

				line->SetLineStyle(1);
//				line->SetLineWidth(1);
				line->DrawLine(0, 0, 0.6, 0);

				ltx->SetTextAlign(12);
				c_sys->cd();
				ATLASLabel(0.19, 0.88, "Preliminary", "", kBlack);

				ltx->SetTextAlign(12);
				ltx->SetTextSize(23);

				if (mode == "RDpT")
				{
					ltx->DrawLatexNDC(0.19, 0.84, "Pb+Pb #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}");
					ltx->DrawLatexNDC(0.19, 0.79, "#it{pp} #sqrt{#font[12]{s}} = 5.02 TeV, 25 pb^{-1}");
				}
				else
				{
					if (dataset_type == "_PbPb") ltx->DrawLatexNDC(0.19, 0.84, "Pb+Pb #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}");
					if (dataset_type == "_pp") ltx->DrawLatexNDC(0.19, 0.84, "#it{pp} #sqrt{#font[12]{s}} = 5.02 TeV, 25 pb^{-1}");
				}

				c_sys->Print(Form("output_pdf_nominal/systematics/%s_dR_sys%s_error_trk%i_jet%i_cent%i.pdf",mode.c_str(), dataset_type.c_str(), i_trk, i_jet, i_cent));

			}
		}
	}

	cout << "######### Done Draw Systematics #########" << endl;

}




