#include "../functions/global_variables.h"

void draw_conf_plots(bool isMC = 0)
{
	cout << "######### DOING COMP_ChPS #########" << endl;

	SetAtlasStyle();
	gErrorIgnoreLevel = 3001;
	std::string did = "data";
	if (isMC) did = "MC";

	cout << Form("Doing in %s mode", did.c_str()) << endl;

	TFile *f_RDpT = new TFile(Form("output_pdf/root/final_RDpT_%s.root", did.c_str()));
	TFile *f_sys = new TFile(Form("output_pdf/root/final_RDpT_sys_%s.root", did.c_str()));


	TAxis* dR_binning = (TAxis*)f_RDpT->Get("dR_binning");
	TAxis* jetpT_binning = (TAxis*)f_RDpT->Get("jetpT_binning");
	TAxis* trkpT_binning = (TAxis*)f_RDpT->Get("trkpT_binning");

	int N_dR = dR_binning->GetNbins();
	int N_jetpt = jetpT_binning->GetNbins();
	int N_trkpt = trkpT_binning->GetNbins();

	//indR
	vector<vector<vector<TH1*>>> h_ChPS_final_ratio_indR (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_final_sys_indR (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));

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
				h_ChPS_final_ratio_indR.at(i_trk).at(i_cent).at(i_jet) = (TH1*)f_RDpT->Get(name.c_str());
			}
		}
	}


	// drawing
	{
		cout << "Doing Final ChPS ratio (PbPb/pp) plots in dR" << endl;

		TCanvas *c_ChPS_final_CONF = new TCanvas("c_ChPS_final_CONF","c_ChPS_final_CONF",800,600);
		TLegend *legend_ChPS_final_CONF = new TLegend(0.19, 0.19, 0.40, 0.35, "","brNDC");
		legend_ChPS_final_CONF->SetTextFont(43);
		legend_ChPS_final_CONF->SetBorderSize(0);
		legend_ChPS_final_CONF->SetTextSize(16);


		int jet_itr = 0;
		for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
		{
			string jet_label = Form("%1.0f < p_{T}^{Jet} < %1.0f", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));

			bool first_pass_cent = true;
			for (int i_cent = 0; i_cent < 6; i_cent++)
			{
				string centrality = num_to_cent(31,i_cent);

				int trk_itr = 0;
				c_ChPS_final_CONF->cd();
				c_ChPS_final_CONF->Clear();

				for (int i_trk = 0; i_trk < N_trkpt; i_trk++)
				{

					if (i_trk != 2 && i_trk != 4 && i_trk != 6) continue;

					string trk_label = Form("%1.2f < p_{T}^{Trk} < %1.2f", trkpT_binning->GetBinLowEdge(i_trk+1), trkpT_binning->GetBinUpEdge(i_trk+1));

					SetHStyle_smallify(h_ChPS_final_ratio_indR.at(i_trk).at(i_cent).at(i_jet), trk_itr, 0);

					if (jet_itr == 0 && first_pass_cent) legend_ChPS_final_CONF->AddEntry(h_ChPS_final_ratio_indR.at(i_trk).at(i_cent).at(i_jet),trk_label.c_str(),"lp");

					h_ChPS_final_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(0, 0.6);
					h_ChPS_final_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->GetYaxis()->SetRangeUser(0, 2.5);
					h_ChPS_final_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->GetYaxis()->SetNdivisions(504);
					if (trk_itr == 0) h_ChPS_final_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->Draw("");
					else h_ChPS_final_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->Draw("same");

					trk_itr++;

				} // end trk loop

				c_ChPS_final_CONF->cd();
				ltx->SetTextAlign(32);
				ltx->SetTextSize(18);
				ltx->DrawLatexNDC(0.93, 0.90, Form("%s", jet_label.c_str()));
				ltx->DrawLatexNDC(0.93, 0.85, Form("%s", centrality.c_str()));
				line->DrawLine(0, 1, 1.2, 1);
				legend_ChPS_final_CONF->Draw();
				ATLASLabel(0.19, 0.85, "Internal, 0.49 nb^{-1}", "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV", kBlack);
				first_pass_cent = false;

				pdf_label = "";
				if (i_jet == jet_pt_start && i_cent == 0) pdf_label = "(";
				if (i_jet == jet_pt_end-1 && i_cent == 5) pdf_label = ")";
				c_ChPS_final_CONF->Print(Form("output_pdf/ChPS_final_ratio_dR_CONF_%s.pdf%s", did.c_str(), pdf_label.c_str()), Form("Title:jetpt%i_cent%i", i_jet, i_cent));

			} //end cent loop

			jet_itr++;
		} //end jet loop
	}

	cout << "######### DONE CONF PLOTS #########" << endl << endl;;


}

