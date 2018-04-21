#include "../functions/global_variables.h"

static const int N_JET_Y = 4;
double jet_y_binning[N_JET_Y+1] = {0, 0.3, 0.8, 1.2, 1.3};
void draw_spectra(string config_file = "ff_config.cfg")
{
	cout << "######### RUNNING DRAW_Spectra #########" << endl;

	SetAtlasStyle();
	gErrorIgnoreLevel = 3001;

	//	##############	Reading config	##############"
	TEnv *m_config = new TEnv();
	m_config->ReadFile(config_file.c_str(), EEnvLevel(1));

	std::string dataset_type = "PbPb"; dataset_type = m_config->GetValue("dataset_type", dataset_type.c_str());
	int isMC = 1; isMC = m_config->GetValue("isMC", isMC);
	int sys_mode = -1; sys_mode = m_config->GetValue("sys_mode", sys_mode);

	int centrality_scheme = 31; centrality_scheme = m_config->GetValue("centrality_scheme", centrality_scheme);
	int verbose = 0; verbose = m_config->GetValue("verbose", verbose);

	std::string did = "data";
	if (isMC) did = "MC";

	if (verbose) m_config->Print();
	//	##############	Config done	##############"
	std::string sys_path = "";
	if (sys_mode == 0) sys_path = Form("_nominal");
	if (sys_mode > 0) sys_path = Form("_sys%i", sys_mode);


	TFile *f_input = new TFile(Form("output_pdf%s/root/raw_unfolded_%s_%s.root", sys_path.c_str(), did.c_str(), dataset_type.c_str()));
	cout << "Using files:" << endl;
	cout << f_input->GetName() << endl;


	vector<vector<TH1*>> h_raw (n_cent_cuts, vector<TH1*> (N_JET_Y+1));
	vector<vector<TH1*>> h_unfolded (n_cent_cuts, vector<TH1*> (N_JET_Y+1));
	vector<vector<TH1*>> h_true (n_cent_cuts, vector<TH1*> (N_JET_Y+1));
	vector<vector<TH1*>> h_closure (n_cent_cuts, vector<TH1*> (N_JET_Y+1));
	vector<vector<TH1*>> h_response_matrix (n_cent_cuts, vector<TH1*> (N_JET_Y+1));



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

	TLegend *legend = new TLegend(0.19, 0.44, 0.5, 0.7);
	legend->SetTextFont(43);
	legend->SetTextSize(12);
	legend->SetBorderSize(0);


	for (int i_y = 0; i_y < N_JET_Y+1; i_y++)
	{
		for (int i_cent = 0; i_cent < n_cent_cuts; i_cent++)
		{
			if (dataset_type == "PbPb" && i_cent == 6) continue;
			if (dataset_type == "pp" && i_cent < 6) continue;


			name = Form("h_reco_jet_y%i_c%i", i_y, i_cent);
			h_raw.at(i_cent).at(i_y) = (TH1*)f_input->Get(name.c_str());
			h_raw.at(i_cent).at(i_y)->Scale(1.,"width");

			name = Form("h_unfolded_jet_y%i_c%i", i_y, i_cent);
			h_unfolded.at(i_cent).at(i_y) = (TH1*)f_input->Get(name.c_str());
			h_unfolded.at(i_cent).at(i_y)->Scale(1.,"width");

			name = Form("h_true_jet_y%i_c%i", i_y, i_cent);
			h_true.at(i_cent).at(i_y) = (TH1*)f_input->Get(name.c_str());
			h_true.at(i_cent).at(i_y)->Scale(1.,"width");

			name = Form("h_unfolded_jet_y%i_c%i", i_y, i_cent);
			h_closure.at(i_cent).at(i_y) = (TH1*)h_unfolded.at(i_cent).at(i_y)->Clone(name.c_str());
			h_closure.at(i_cent).at(i_y)->Divide(h_true.at(i_cent).at(i_y));

			name = Form("h_response_matrix_jet_y%i_c%i", i_y, i_cent);
			h_response_matrix.at(i_cent).at(i_y) = (TH2*)f_input->Get(name.c_str());


		}

	}


	{
		cout << "Doing Final jet spectra plots" << endl;

		TCanvas *c_spect_closure = new TCanvas("c_spect_closure","c_spect_closure",800,600);
		if (dataset_type == "pp") c_spect_closure->SetCanvasSize(600,600);
		TLegend *legend_spect_closure = new TLegend(0.20, 0.43, 0.40, 0.60, "","brNDC");
		legend_spect_closure->SetTextFont(43);
		legend_spect_closure->SetBorderSize(0);
		if (dataset_type == "pp") legend_spect_closure->SetTextSize(12);
		if (dataset_type == "PbPb") legend_spect_closure->SetTextSize(8);

		TCanvas *c_resp_matrix = new TCanvas("c_resp_matrix","c_resp_matrix",800,600);
		if (dataset_type == "pp") c_resp_matrix->SetCanvasSize(600,600);

		for (int i_y = 0; i_y < N_JET_Y+1; i_y++)
		{
			string y_label = Form("%1.2f < |y^{jet}| < %1.2f", jet_y_binning[i_y], jet_y_binning[i_y+1]);
			if (i_y == N_JET_Y) y_label = Form("%1.2f < |y^{jet}| < %1.2f", jet_y_binning[0], jet_y_binning[N_JET_Y]);

			c_spect_closure->cd();
			c_spect_closure->Clear();
			if (dataset_type == "PbPb") c_spect_closure->Divide(3,2);

			c_resp_matrix->cd();
			c_resp_matrix->Clear();
			if (dataset_type == "PbPb") c_resp_matrix->Divide(3,2);

			bool first_pass_cent = true;
			for (int i_cent = 0; i_cent < n_cent_cuts; i_cent++)
			{
				string centrality = num_to_cent(centrality_scheme,i_cent);
				if (dataset_type == "pp" && i_cent < 6) continue;
				else if (dataset_type == "PbPb" && i_cent == 6) continue;

				if (dataset_type == "pp") c_spect_closure->Divide(1,2);
				if (dataset_type == "PbPb") c_spect_closure->cd(i_cent+1)->Divide(1,2);

				if (isMC)
				{
					if (dataset_type == "pp") c_spect_closure->cd(1);
					if (dataset_type == "PbPb") c_spect_closure->cd(i_cent+1)->cd(1);
					gPad->SetPad(0,0.40,0.95,0.95);
					gPad->SetTopMargin(0.05);
					gPad->SetBottomMargin(0);
					gPad->SetRightMargin(0);

					if (dataset_type == "pp") c_spect_closure->cd(2);
					if (dataset_type == "PbPb") c_spect_closure->cd(i_cent+1)->cd(2);
					gPad->SetPad(0,0.0,0.95,0.40);
					gPad->SetTopMargin(0);
					gPad->SetBottomMargin(0.35);
					gPad->SetRightMargin(0);
				}

				SetHStyle_smallify(h_raw.at(i_cent).at(i_y), 0, doSmall);
				SetHStyle_smallify(h_unfolded.at(i_cent).at(i_y), 1, doSmall);
				SetHStyle_smallify(h_true.at(i_cent).at(i_y), 2, doSmall);
				SetHStyle_smallify(h_closure.at(i_cent).at(i_y), 6, doSmall);

				double jet_pt_lo = 126;
				double jet_pt_hi = 620;
				h_raw.at(i_cent).at(i_y)->GetXaxis()->SetRangeUser(jet_pt_lo, jet_pt_hi);
				h_unfolded.at(i_cent).at(i_y)->GetXaxis()->SetRangeUser(jet_pt_lo, jet_pt_hi);
				h_true.at(i_cent).at(i_y)->GetXaxis()->SetRangeUser(jet_pt_lo, jet_pt_hi);
				h_closure.at(i_cent).at(i_y)->GetXaxis()->SetRangeUser(jet_pt_lo, jet_pt_hi);

//				h_true.at(i_cent).at(i_y)->GetXaxis()->SetNdivisions(504);
//				h_closure.at(i_cent).at(i_y)->GetXaxis()->SetNdivisions(504);

				h_true.at(i_cent).at(i_y)->GetYaxis()->SetTitle("Normalized counts [GeV^{-1}]");
				h_closure.at(i_cent).at(i_y)->GetYaxis()->SetTitle("MC Closure");

				if (doSmall)
				{
					h_true.at(i_cent).at(i_y)->GetYaxis()->SetTitleSize(12);
					h_closure.at(i_cent).at(i_y)->GetYaxis()->SetTitleSize(12);
				}

				h_closure.at(i_cent).at(i_y)->GetXaxis()->SetTitle("p_{T}^{Jet} [GeV]");

				h_closure.at(i_cent).at(i_y)->GetYaxis()->SetRangeUser(0.84, 1.16);
//				h_closure.at(i_cent).at(i_y)->GetYaxis()->SetNdivisions(504);


				if (i_y == 0 && first_pass_cent)
				{
					legend_spect_closure->AddEntry(h_raw.at(i_cent).at(i_y),"Raw","lp");
					legend_spect_closure->AddEntry(h_true.at(i_cent).at(i_y),"True","lp");
					legend_spect_closure->AddEntry(h_unfolded.at(i_cent).at(i_y),"Unfolded","lp");
					legend_spect_closure->AddEntry(h_closure.at(i_cent).at(i_y),"Closure","lp");
				}

				if (isMC)
				{
					if (dataset_type == "pp") c_spect_closure->cd(1);
					if (dataset_type == "PbPb") c_spect_closure->cd(i_cent+1)->cd(1);
				}
				else
				{
					if (dataset_type == "pp") c_spect_closure->cd();
					if (dataset_type == "PbPb") c_spect_closure->cd(i_cent+1);
				}

				h_true.at(i_cent).at(i_y)->Draw("");
				h_raw.at(i_cent).at(i_y)->Draw("same");
				h_unfolded.at(i_cent).at(i_y)->Draw("same");
				gPad->SetLogx(0);
				gPad->SetLogy();

				if (isMC)
				{
					if (dataset_type == "pp") c_spect_closure->cd(2);
					if (dataset_type == "PbPb") c_spect_closure->cd(i_cent+1)->cd(2);
					h_closure.at(i_cent).at(i_y)->GetYaxis()->SetNdivisions(504);
					h_closure.at(i_cent).at(i_y)->Draw("");
					if (dataset_type == "pp") h_closure.at(i_cent).at(i_y)->GetXaxis()->SetTitleOffset(3.2);
					if (dataset_type == "PbPb")
					{
						h_closure.at(i_cent).at(i_y)->GetXaxis()->SetTitleOffset(5);
						h_closure.at(i_cent).at(i_y)->GetYaxis()->SetTitleOffset(4.);
						h_true.at(i_cent).at(i_y)->GetYaxis()->SetTitleOffset(4.);

					}
					line->DrawLine(jet_pt_lo, 1, jet_pt_hi, 1);
					gPad->SetLogx(0);
					gPad->SetLogy(0);

				}


				if (dataset_type == "pp") c_spect_closure->cd();
				if (dataset_type == "PbPb") c_spect_closure->cd(i_cent+1);

				ltx->SetTextAlign(32);
				ltx->SetTextSize(12);
				if (dataset_type == "pp") ltx->SetTextSize(16);
				ltx->DrawLatexNDC(0.93, 0.88, Form("%s", y_label.c_str()));
				ltx->DrawLatexNDC(0.93, 0.82, Form("%s", centrality.c_str()));
				ltx->DrawLatexNDC(0.93, 0.95, Form("%s %s", dataset_type.c_str(), did.c_str()));

				legend_spect_closure->Draw();




				if (dataset_type == "pp") c_resp_matrix->cd();
				if (dataset_type == "PbPb") c_resp_matrix->cd(i_cent+1);
				h_response_matrix.at(i_cent).at(i_y)->GetXaxis()->SetTitle("p_{T}^{reco jet}");
				h_response_matrix.at(i_cent).at(i_y)->GetYaxis()->SetTitle("p_{T}^{truth jet}");
//				h_response_matrix.at(i_cent).at(i_y)->GetZaxis()->SetRangeUser(1e-10, 1e1);

				gPad->SetRightMargin(0.15);

				h_response_matrix.at(i_cent).at(i_y)->Draw("colz");
				gPad->SetLogx();
				gPad->SetLogy();
				gPad->SetLogz();
				ltx->SetTextAlign(12);
				ltx->SetTextSize(12);
				if (dataset_type == "pp") ltx->SetTextSize(16);
				ltx->DrawLatexNDC(0.19, 0.90, Form("%s", y_label.c_str()));
				ltx->DrawLatexNDC(0.19, 0.85, Form("%s", centrality.c_str()));
				ltx->DrawLatexNDC(0.19, 0.98, Form("%s %s", dataset_type.c_str(), did.c_str()));

				first_pass_cent = false;
			} //end cent loop

			pdf_label = "";
			if (i_y == 0) pdf_label = "(";
			if (i_y == N_JET_Y) pdf_label = ")";
			c_spect_closure->Print(Form("output_pdf%s/%s/spect_closure_%s_%s.pdf%s", sys_path.c_str(), dataset_type.c_str(), dataset_type.c_str(), did.c_str(), pdf_label.c_str()), Form("Title:y%i", i_y));
			c_resp_matrix->Print(Form("output_pdf%s/%s/resp_matrix_jet_%s_%s.pdf%s", sys_path.c_str(), dataset_type.c_str(), dataset_type.c_str(), did.c_str(), pdf_label.c_str()), Form("Title:y%i", i_y));

			if (i_y == 0)
			{
				c_resp_matrix->Print(Form("output_pdf%s/%s/resp_matrix_jet_%s_%s.root", sys_path.c_str(), dataset_type.c_str(), dataset_type.c_str(), did.c_str()));

			}
		} //end dR loop
	}

	cout << "######### DONE DRAW_Spectra #########" << endl << endl;;


}

