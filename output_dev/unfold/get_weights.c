#include "../functions/global_variables.h"
#include "draw_functions.c"

void get_weights(string config_file = "ff_config.cfg")
{
	SetAtlasStyle();
	gErrorIgnoreLevel = 3001;

	//	##############	Reading config	##############"
	TEnv *m_config = new TEnv();
	m_config->ReadFile(config_file.c_str(), EEnvLevel(1));

	std::string dataset_type = "PbPb"; dataset_type = m_config->GetValue("dataset_type", dataset_type.c_str());

	int centrality_scheme = 31; centrality_scheme = m_config->GetValue("centrality_scheme", centrality_scheme);
	int verbose = 0; verbose = m_config->GetValue("verbose", verbose);
	int draw_mode = 1; draw_mode = m_config->GetValue("draw_mode", draw_mode);

	if (verbose) m_config->Print();
	//	##############	Config done	##############"


	TFile *f_mc = new TFile(Form("output_pdf_nominal/root/final_ChPS_MC_%s.root", dataset_type.c_str()));
	TFile *f_data = new TFile(Form("output_pdf_nominal/root/final_ChPS_data_%s.root", dataset_type.c_str()));
	TFile *f_weights = new TFile(Form("output_pdf_nominal/root/shape_spectra_weights_%s.root", dataset_type.c_str()), "recreate");

	cout << "Using: " << endl;
	cout << f_mc->GetName() << endl;
	cout << f_data->GetName() << endl;

	TAxis* dR_binning = (TAxis*)f_mc->Get("dR_binning");
	TAxis* jetpT_binning = (TAxis*)f_mc->Get("jetpT_binning");
	TAxis* trkpT_binning = (TAxis*)f_mc->Get("trkpT_binning");

	int N_dR = dR_binning->GetNbins();
	int N_jetpt = jetpT_binning->GetNbins();
	int N_trkpt = trkpT_binning->GetNbins();

	double trk_pt_lo = 1.;
	double trk_pt_hi = 150.;

	double ratio_lo = 0;
	double ratio_hi = 2;

	int jet_pt_start = 7;
	int jet_pt_end = 11;

	int trk_pt_start = 2;
	int trk_pt_end = 5;

	bool doSmall;
	if (dataset_type == "pp") doSmall = false;
	if (dataset_type == "PbPb") doSmall = true;

	TCanvas *c = new TCanvas("c","c",800,600);
	string name;

	TLine *line = new TLine();
	line->SetLineColor(kBlack);

	TLatex *ltx = new TLatex();
	ltx->SetTextFont(43);
	ltx->SetTextSize(16);
	ltx->SetTextAlign(12);

	vector<vector<vector<TH1*>>> h_ChPS_raw_MC (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_raw_data (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));

	vector<vector<vector<TH1*>>> h_ChPS_ratio (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));

	for (int i_cent = 0; i_cent < n_cent_cuts; i_cent++)
	{
		for (int i_dR = 0; i_dR < N_dR; i_dR++)
		{
			for (int i_jet = 0; i_jet < N_jetpt; i_jet++)
			{
				if (dataset_type == "PbPb" && i_cent == 6) continue;
				if (dataset_type == "pp" && i_cent < 6) continue;

				name = Form("h_ChPS_raw_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet);
				h_ChPS_raw_MC[i_dR][i_cent][i_jet] = (TH1*)((TH1*)f_mc->Get(name.c_str()))->Clone(Form("%s_MC", name.c_str()));
				h_ChPS_raw_data[i_dR][i_cent][i_jet] = (TH1*)((TH1*)f_data->Get(name.c_str()))->Clone(Form("%s_data", name.c_str()));

				h_ChPS_ratio[i_dR][i_cent][i_jet] = (TH1*)h_ChPS_raw_data[i_dR][i_cent][i_jet]->Clone(Form("%s_ratio", name.c_str()));
				h_ChPS_ratio[i_dR][i_cent][i_jet]->Divide(h_ChPS_raw_MC[i_dR][i_cent][i_jet]);

				f_weights->cd();
				name = Form("CHPS_weight_%s_dR%i_cent%i_jet%i", dataset_type.c_str(), i_dR, i_cent, i_jet);
				h_ChPS_ratio[i_dR][i_cent][i_jet]->Write(name.c_str());
			}
		}
	}

	//compare raw DpT in data and MC
	{
		cout << "Comparing DpT in data and MC (as function of r)" << endl;

		TCanvas *c_comp_unsub = new TCanvas("c_comp_unsub","c_comp_unsub",900,600);
		TLegend *legend_comp_unsub = new TLegend(0.18, 0.18, 0.30, 0.38, "","brNDC");
		legend_comp_unsub->SetTextFont(43);
		legend_comp_unsub->SetBorderSize(0);
		legend_comp_unsub->SetTextSize(8);


		int jet_itr = 0;
		for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
		{
			string jet_label = Form("%1.0f < #it{p}_{T}^{jet} < %1.0f", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));

			c_comp_unsub->cd();
			c_comp_unsub->Clear();
			if (dataset_type == "PbPb") c_comp_unsub->Divide(3,2);

			bool first_pass = true;

			for (int i_cent = 0; i_cent < n_cent_cuts; i_cent++)
			{
				if (dataset_type == "PbPb" && i_cent == 6) continue;
				if (dataset_type == "pp" && i_cent < 6) continue;

				string centrality = num_to_cent(31,i_cent);


				int dR_itr = 0;
				if (dataset_type == "PbPb") c_comp_unsub->cd(i_cent+1);
				else c_comp_unsub->cd();
				
				for (int i_dR = 0; i_dR < N_dR; i_dR++)
				{
					string dR_label = Form("%1.2f < #it{r} < %1.2f", dR_binning->GetBinLowEdge(i_dR+1), dR_binning->GetBinUpEdge(i_dR+1));

					SetHStyle_smallify(h_ChPS_ratio.at(i_dR).at(i_cent).at(i_jet),i_dR, doSmall);

					if (i_jet == jet_pt_start && first_pass)
					{
						legend_comp_unsub->AddEntry(h_ChPS_ratio.at(i_dR).at(i_cent).at(i_jet),dR_label.c_str(),"lp");
					}

					h_ChPS_ratio.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetRangeUser(0,2);
					h_ChPS_ratio.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetNdivisions(504);

					if (dR_itr == 0) h_ChPS_ratio.at(i_dR).at(i_cent).at(i_jet)->Draw();
					else h_ChPS_ratio.at(i_dR).at(i_cent).at(i_jet)->Draw("same");

					gPad->SetLogx();
					line->DrawLine(1, 1, 1E2, 1);

					dR_itr++;

				} // end trk loop

				if (dataset_type == "PbPb") c_comp_unsub->cd(i_cent+1);
				else c_comp_unsub->cd();
				ltx->SetTextAlign(32);
				ltx->SetTextSize(12);
				ltx->DrawLatexNDC(0.93, 0.90, Form("Raw Unsubtracted Data/MC"));
				ltx->DrawLatexNDC(0.93, 0.85, Form("%s", jet_label.c_str()));
				ltx->DrawLatexNDC(0.93, 0.80, Form("%s", centrality.c_str()));
				legend_comp_unsub->Draw();

				first_pass = false;

			} //end cent loop

			
			string pdf_label = "";
			if (i_jet == jet_pt_start) pdf_label = "(";
			if (i_jet == jet_pt_end-1) pdf_label = ")";
			c_comp_unsub->Print(Form("output_pdf_nominal/%s/ChPS_weights_%s.pdf%s", dataset_type.c_str(), dataset_type.c_str(), pdf_label.c_str()), Form("Title:jetpt%i", i_jet));


			jet_itr++;
		} //end jet loop
	}

	//testing
	vector<vector<vector<TH1*>>> h_ChPS_new (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));


	{
		cout << "testing weights" << endl;

		TCanvas *c_comp_unsub = new TCanvas("c_comp_unsub","c_comp_unsub",900,600);
		TLegend *legend_comp_unsub = new TLegend(0.18, 0.18, 0.30, 0.38, "","brNDC");
		legend_comp_unsub->SetTextFont(43);
		legend_comp_unsub->SetBorderSize(0);
		legend_comp_unsub->SetTextSize(8);


		int jet_itr = 0;
		for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
		{
			string jet_label = Form("%1.0f < #it{p}_{T}^{jet} < %1.0f", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));

			for (int i_cent = 0; i_cent < n_cent_cuts; i_cent++)
			{
				if (dataset_type == "PbPb" && i_cent == 6) continue;
				if (dataset_type == "pp" && i_cent < 6) continue;
				string centrality = num_to_cent(31,i_cent);

				for (int i_dR = 0; i_dR < N_dR; i_dR++)
				{
					string dR_label = Form("%1.2f < #it{r} < %1.2f", dR_binning->GetBinLowEdge(i_dR+1), dR_binning->GetBinUpEdge(i_dR+1));


					name = Form("new_jet%i_cent%i_dR%i", i_jet, i_cent, i_dR);
					h_ChPS_new[i_dR][i_cent][i_jet] = (TH1*)h_ChPS_raw_MC[i_dR][i_cent][i_jet]->Clone(name.c_str());

					for (int i = 0 ; i < N_trkpt; i++)
					{
						double factor = h_ChPS_ratio[i_dR][i_cent][i_jet]->GetBinContent(i+1);
						double orig = h_ChPS_raw_MC[i_dR][i_cent][i_jet]->GetBinContent(i+1);
						h_ChPS_new[i_dR][i_cent][i_jet]->SetBinContent(i+1, orig*factor);
					}

					SetHStyle_smallify(h_ChPS_raw_MC[i_dR][i_cent][i_jet],0, doSmall);
					SetHStyle_smallify(h_ChPS_raw_data[i_dR][i_cent][i_jet],1, doSmall);
					SetHStyle_smallify(h_ChPS_new[i_dR][i_cent][i_jet],2, doSmall);

					legend_comp_unsub->Clear();
					legend_comp_unsub->AddEntry(h_ChPS_raw_MC[i_dR][i_cent][i_jet],"mc","lp");
					legend_comp_unsub->AddEntry(h_ChPS_raw_data[i_dR][i_cent][i_jet],"data","lp");
					legend_comp_unsub->AddEntry(h_ChPS_new[i_dR][i_cent][i_jet],"ReW","lp");

					h_ChPS_raw_MC[i_dR][i_cent][i_jet]->Draw();
					h_ChPS_raw_data[i_dR][i_cent][i_jet]->Draw("same");
					h_ChPS_new[i_dR][i_cent][i_jet]->Draw("same");

					gPad->SetLogx();
					gPad->SetLogy();


					c_comp_unsub->cd();
					ltx->SetTextAlign(32);
					ltx->SetTextSize(12);
					ltx->DrawLatexNDC(0.93, 0.90, Form("Raw Unsubtracted Data/MC"));
					ltx->DrawLatexNDC(0.93, 0.85, Form("%s", jet_label.c_str()));
					ltx->DrawLatexNDC(0.93, 0.80, Form("%s", centrality.c_str()));
					ltx->DrawLatexNDC(0.93, 0.75, Form("%s", dR_label.c_str()));
					legend_comp_unsub->Draw();

					string pdf_label = "";
					if (i_jet == jet_pt_start && i_cent == 0 && i_dR == 0) pdf_label = "(";
					else if (i_jet == jet_pt_end-1 && i_cent == 5 && i_dR == N_dR-1) pdf_label = ")";
					c_comp_unsub->Print(Form("output_pdf_nominal/%s/ChPS_weights_test_%s.pdf%s", dataset_type.c_str(), dataset_type.c_str(), pdf_label.c_str()), Form("Title:jet%i_cent%i_dR%i", i_jet, i_cent, i_dR));

				} // end trk loop
			} //end cent loop
		} //end jet loop
	}

}
