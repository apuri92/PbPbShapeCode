#include "../functions/global_variables.h"

void draw_conf_plots(bool isMC = 0)
{
	cout << "######### DOING COMP_ChPS #########" << endl;

	SetAtlasStyle();
	gErrorIgnoreLevel = 3001;
	std::string did = "data";
	if (isMC) did = "MC";

	cout << Form("Doing in %s mode", did.c_str()) << endl;

	TFile *f_RDpT = new TFile(Form("output_pdf_nominal/root/final_RDpT_%s.root", did.c_str()));
	TFile *f_sys = new TFile(Form("output_pdf_nominal/root/final_RDpT_sys_%s.root", did.c_str()));


	TAxis* dR_binning = (TAxis*)f_RDpT->Get("dR_binning");
	TAxis* jetpT_binning = (TAxis*)f_RDpT->Get("jetpT_binning");
	TAxis* trkpT_binning = (TAxis*)f_RDpT->Get("trkpT_binning");

	int N_dR = dR_binning->GetNbins();
	int N_jetpt = jetpT_binning->GetNbins();
	int N_trkpt = trkpT_binning->GetNbins();

	//indR
	vector<vector<vector<TH1*>>> h_ChPS_final_ratio_indR (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_final_sys_Totalpos_indR (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_final_sys_Totalneg_indR (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TGraphAsymmErrors*>>> g_ChPS_final_sys_final_indR (N_trkpt, vector<vector<TGraphAsymmErrors*>> (n_cent_cuts, vector<TGraphAsymmErrors*> (N_jetpt)));

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
	double ratio_hi = 2.5;

	int jet_pt_start = 7;
	int jet_pt_end = 11;


	for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
	{
		for (int i_cent = 0; i_cent < 6; i_cent++)
		{
			for (int i_trk = 0; i_trk < N_trkpt; i_trk++)
			{
				name = Form("h_ChPS_RDpT_indR_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_ChPS_final_ratio_indR[i_trk][i_cent][i_jet] = (TH1*)f_RDpT->Get(name.c_str());

				name = Form("h_ChPS_RDpT_indR_trk%i_cent%i_jetpt%i_sys_Totalpos", i_trk, i_cent, i_jet);
				h_ChPS_final_sys_Totalpos_indR[i_trk][i_cent][i_jet] = (TH1*)f_sys->Get(name.c_str());

				name = Form("h_ChPS_RDpT_indR_trk%i_cent%i_jetpt%i_sys_Totalneg", i_trk, i_cent, i_jet);
				h_ChPS_final_sys_Totalneg_indR[i_trk][i_cent][i_jet] = (TH1*)f_sys->Get(name.c_str());

				g_ChPS_final_sys_final_indR[i_trk][i_cent][i_jet] = new TGraphAsymmErrors(h_ChPS_final_ratio_indR[i_trk][i_cent][i_jet]);
				g_ChPS_final_sys_final_indR[i_trk][i_cent][i_jet]->GetYaxis()->SetTitle(h_ChPS_final_ratio_indR[i_trk][i_cent][i_jet]->GetYaxis()->GetTitle());
				g_ChPS_final_sys_final_indR[i_trk][i_cent][i_jet]->GetXaxis()->SetTitle(h_ChPS_final_ratio_indR[i_trk][i_cent][i_jet]->GetXaxis()->GetTitle());

				double nom, err_hi, err_lo, r_position, r_width;
				for (int i_dR = 0 ; i_dR < N_dR; i_dR++)
				{
					r_position = h_ChPS_final_ratio_indR[i_trk][i_cent][i_jet]->GetBinCenter(i_dR+1);
					r_width = h_ChPS_final_ratio_indR[i_trk][i_cent][i_jet]->GetBinWidth(i_dR+1);

					nom = h_ChPS_final_ratio_indR[i_trk][i_cent][i_jet]->GetBinContent(i_dR+1);
					err_hi = nom * (h_ChPS_final_sys_Totalpos_indR[i_trk][i_cent][i_jet]->GetBinContent(i_dR+1));
					err_lo = fabs(nom * (h_ChPS_final_sys_Totalneg_indR[i_trk][i_cent][i_jet]->GetBinContent(i_dR+1)));

					g_ChPS_final_sys_final_indR[i_trk][i_cent][i_jet]->SetPoint(i_dR, r_position, nom );
					g_ChPS_final_sys_final_indR[i_trk][i_cent][i_jet]->SetPointError(i_dR, r_width/2, r_width/2, err_lo, err_hi);


				}

			}
		}
	}


	// drawing
	{
		cout << "Doing Final ChPS ratio (PbPb/pp) plots in dR" << endl;

		TCanvas *canvas = new TCanvas("canvas","canvas",800,600);
		TLegend *legend = new TLegend(0.19, 0.68, 0.40, 0.82, "","brNDC");
		legend->SetTextFont(43);
		legend->SetBorderSize(0);
		legend->SetTextSize(13);


		int jet_itr = 0;
		for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
		{
			string jet_label = Form("%1.0f < p_{T}^{Jet} < %1.0f", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));

			bool first_pass_cent = true;
			for (int i_cent = 0; i_cent < 6; i_cent++)
			{
				string centrality = num_to_cent(31,i_cent);

				int trk_itr = 0;
				canvas->cd();
				canvas->Clear();

				for (int i_trk = 0; i_trk < N_trkpt; i_trk++)
				{

					if (i_trk != 2 && i_trk != 4 && i_trk != 6 && i_trk != 8) continue;


					string trk_label = Form("%1.2f < p_{T}^{Trk} < %1.2f GeV", trkpT_binning->GetBinLowEdge(i_trk+1), trkpT_binning->GetBinUpEdge(i_trk+1));

					SetHStyle_smallify(h_ChPS_final_ratio_indR[i_trk][i_cent][i_jet], trk_itr, 0);
					SetHStyle_graph_smallify(g_ChPS_final_sys_final_indR[i_trk][i_cent][i_jet], trk_itr, 0);

					if (jet_itr == 0 && first_pass_cent) legend->AddEntry(h_ChPS_final_ratio_indR[i_trk][i_cent][i_jet],trk_label.c_str(),"lp");

					g_ChPS_final_sys_final_indR[i_trk][i_cent][i_jet]->GetXaxis()->SetRangeUser(0, 0.6);
					g_ChPS_final_sys_final_indR[i_trk][i_cent][i_jet]->GetYaxis()->SetRangeUser(0, 2.5);
					g_ChPS_final_sys_final_indR[i_trk][i_cent][i_jet]->GetYaxis()->SetNdivisions(504);

					h_ChPS_final_ratio_indR[i_trk][i_cent][i_jet]->GetXaxis()->SetRangeUser(0, 0.6);
					h_ChPS_final_ratio_indR[i_trk][i_cent][i_jet]->GetYaxis()->SetRangeUser(0, 2.5);
					h_ChPS_final_ratio_indR[i_trk][i_cent][i_jet]->GetYaxis()->SetNdivisions(504);


					if (trk_itr == 0) g_ChPS_final_sys_final_indR[i_trk][i_cent][i_jet]->Draw("a 2");
					else g_ChPS_final_sys_final_indR[i_trk][i_cent][i_jet]->Draw("same 2");
					h_ChPS_final_ratio_indR[i_trk][i_cent][i_jet]->Draw("same");

					trk_itr++;

				} // end trk loop

				canvas->cd();
				ltx->SetTextAlign(32);
				ltx->SetTextSize(16);
				ltx->DrawLatexNDC(0.93, 0.90, Form("%s", jet_label.c_str()));
				ltx->DrawLatexNDC(0.93, 0.85, Form("%s", centrality.c_str()));
				line->DrawLine(0, 1, 0.6, 1);
				legend->Draw();
				ATLASLabel(0.19, 0.88, "Internal", "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);
				first_pass_cent = false;

				if ((i_jet != 7 && i_jet != 9) || (i_cent != 0 && i_cent != 5) ) continue;
				canvas->Print(Form("output_pdf_nominal/conf/ChPS_final_ratio_dR_CONF_%s_jet%i_cent%i.pdf", did.c_str(), i_jet, i_cent));


			} //end cent loop

			jet_itr++;
		} //end jet loop
	}

	{
		cout << "Doing Final RDpT for different jets for one track pT" << endl;

		TCanvas *canvas = new TCanvas("canvas","canvas",800,600);
		TLegend *legend = new TLegend(0.2619048,0.5965217,0.3521303,0.746087,NULL,"brNDC");
		legend->SetTextFont(43);
		legend->SetBorderSize(0);
		legend->SetTextSize(13);
		TLegend *legend_open = new TLegend(0.197995,0.5965217,0.4160401,0.746087,NULL,"brNDC");
		legend_open->SetTextFont(43);
		legend_open->SetBorderSize(0);
		legend_open->SetTextSize(13);

		bool first_pass_cent = true;
		for (int i_cent = 0; i_cent < 6; i_cent++)
		{
			string centrality = num_to_cent(31,i_cent);


			canvas->cd();
			canvas->Clear();

			int trk_itr = 0;
			for (int i_trk = 0; i_trk < N_trkpt; i_trk++)
			{
				string trk_label = Form("%1.2f < p_{T}^{Trk} < %1.2f GeV", trkpT_binning->GetBinLowEdge(i_trk+1), trkpT_binning->GetBinUpEdge(i_trk+1));

				if ((i_trk != 2 && i_trk != 6)) continue;


				int jet_itr = 0;
				for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
				{
					string jet_label = Form("%1.0f < p_{T}^{Jet} < %1.0f", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));



					if (i_trk == 2)
					{
						SetHStyle(h_ChPS_final_ratio_indR[i_trk][i_cent][i_jet], jet_itr);
						SetHStyle_graph(g_ChPS_final_sys_final_indR[i_trk][i_cent][i_jet], jet_itr);
						if (first_pass_cent) legend->AddEntry(h_ChPS_final_ratio_indR[i_trk][i_cent][i_jet],jet_label.c_str(),"lp");

					}
					if (i_trk == 6)
					{
						SetHStyle_open(h_ChPS_final_ratio_indR[i_trk][i_cent][i_jet], jet_itr);
						SetHStyle_graph_open(g_ChPS_final_sys_final_indR[i_trk][i_cent][i_jet], jet_itr);
						if (first_pass_cent) legend_open->AddEntry(h_ChPS_final_ratio_indR[i_trk][i_cent][i_jet]," ","lp");

					}

					g_ChPS_final_sys_final_indR[i_trk][i_cent][i_jet]->GetXaxis()->SetRangeUser(0, 0.6);
					g_ChPS_final_sys_final_indR[i_trk][i_cent][i_jet]->GetYaxis()->SetRangeUser(0, 4.5);
					g_ChPS_final_sys_final_indR[i_trk][i_cent][i_jet]->GetYaxis()->SetNdivisions(504);

					h_ChPS_final_ratio_indR[i_trk][i_cent][i_jet]->GetXaxis()->SetRangeUser(0, 0.6);
					h_ChPS_final_ratio_indR[i_trk][i_cent][i_jet]->GetYaxis()->SetRangeUser(0, 4.5);
					h_ChPS_final_ratio_indR[i_trk][i_cent][i_jet]->GetYaxis()->SetNdivisions(504);


					if (jet_itr == 0 && trk_itr == 0) g_ChPS_final_sys_final_indR[i_trk][i_cent][i_jet]->Draw("a2");
					else g_ChPS_final_sys_final_indR[i_trk][i_cent][i_jet]->Draw("same 2");
					h_ChPS_final_ratio_indR[i_trk][i_cent][i_jet]->Draw("same");

					jet_itr++;

				} // end jet loop

				canvas->cd();
				ltx->SetTextAlign(32);
				ltx->SetTextSize(16);
				ltx->DrawLatexNDC(0.93, 0.90, Form("%s", centrality.c_str()));
				line->DrawLine(0, 1, 0.6, 1);
				legend_open->Draw();
				legend->Draw();
				ATLASLabel(0.19, 0.88, "Internal", "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}", kBlack);

				ltx->SetTextAlign(21);
				ltx->SetTextSize(12);


				if (i_trk == 2) ltx->DrawLatexNDC(0.2944862,0.7704348, Form("#splitline{%1.1f < p_{T}^{Trk} }{ < %1.1f GeV}", trkpT_binning->GetBinLowEdge(i_trk+1), trkpT_binning->GetBinUpEdge(i_trk+1) ));
				if (i_trk == 6) ltx->DrawLatexNDC(0.2192982,0.7704348, Form("#splitline{%1.1f < p_{T}^{Trk} }{ < %1.1f GeV}", trkpT_binning->GetBinLowEdge(i_trk+1), trkpT_binning->GetBinUpEdge(i_trk+1) ));


				trk_itr++;

			} //end cent loop

			if (i_cent != 0 && i_cent != 5) continue;

			canvas->Print(Form("output_pdf_nominal/conf/ChPS_final_ratio_dR_CONF_%s_trk_cent%i.pdf", did.c_str(), i_cent));

			first_pass_cent = false;
		} //end jet loop
	}

	cout << "######### DONE CONF PLOTS #########" << endl << endl;;


}

