#include "../functions/global_variables.h"

void UE_factors(string config_file = "ff_config.cfg")
{
	cout << "######### GETTING UE FACTORS #########" << endl;

	SetAtlasStyle();
	gErrorIgnoreLevel = 3001;

	//	##############	Reading config	##############"
	TEnv *m_config = new TEnv();
	m_config->ReadFile(config_file.c_str(), EEnvLevel(1));
	int sys_mode = -1; sys_mode = m_config->GetValue("sys_mode", sys_mode);
	int verbose = 0; verbose = m_config->GetValue("verbose", verbose);

	if (verbose) m_config->Print();
	//	##############	Config done	##############"

	std::string sys_path = "";
	if (sys_mode == 0) sys_path = Form("nominal");
	else if (sys_mode > 50) sys_path = Form("c%i", sys_mode);
	else sys_path = Form("sys%i", sys_mode);
	TFile *input_file = new TFile(Form("../raw_results/%s/FF_MC_out_histo_PbPb_5p02_r001.root", sys_path.c_str()));

	cout << "Using file:" << endl;
	cout << input_file->GetName() << endl;

	if (sys_mode >= 0) sys_path = Form("_%s", sys_path.c_str());
	TFile *UE_factors = new TFile(Form("output_pdf%s/root/UE_factors.root", sys_path.c_str()), "recreate");


	TAxis* dR_binning = (TAxis*)((TH3*)input_file->Get("h_dR_change_jetpt0_cent0"))->GetXaxis();
	TAxis* jetpT_binning = (TAxis*)((TH3*)input_file->Get("ChPS_raw_0_dR0_cent0"))->GetYaxis();
	TAxis* trkpT_binning = (TAxis*)((TH3*)input_file->Get("ChPS_raw_0_dR0_cent0"))->GetXaxis();

	double r_max_range = 1.2;

	string name;
	TCanvas *c1 = new TCanvas("c1","c1",800,400);
	TCanvas *c2 = new TCanvas("c2","c2",800,400);
	TLatex *ltx = new TLatex();
	ltx->SetTextFont(43);
	ltx->SetTextSize(11);
	ltx->SetTextAlign(12);

	int N_dR = dR_binning->GetNbins();
	int N_jetpt = jetpT_binning->GetNbins();
	int N_trkpt = trkpT_binning->GetNbins();

	vector<vector<TH2*>> h_ratio = vector<vector<TH2*>> (13, vector<TH2*> (6));
	vector<vector<vector<TH1*>>> h_ratio_1d = vector<vector<vector<TH1*>>> (13, vector<vector<TH1*>> (6, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ratio_1d_r = vector<vector<vector<TH1*>>> (N_trkpt, vector<vector<TH1*>> (6, vector<TH1*> (N_jetpt)));

	double array_dr_bins[N_dR+1];
	for (int i_dR = 0; i_dR <= N_dR; i_dR++) array_dr_bins[i_dR] = dR_binning->GetBinLowEdge(i_dR+1);
	for (int i_cent = 0; i_cent < 6; i_cent++)
	{
		for (int i_jet = 0; i_jet < N_jetpt; i_jet++)
		{
			for (int i_trk = 0; i_trk < N_trkpt; i_trk++)
			{
				name = Form("h_ratio_1d_r_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_ratio_1d_r[i_trk][i_cent][i_jet] = new TH1D(name.c_str(), name.c_str(), N_dR, array_dr_bins);
			}
		}
	}


	for (int i_dR = 0; i_dR < 13; i_dR++)
	{
		string dr_label = Form("%1.2f < r < %1.2f", dR_binning->GetBinLowEdge(i_dR+1), dR_binning->GetBinUpEdge(i_dR+1));

		c1->Clear();
		c1->Divide(3,2);

		for (int i_cent = 0; i_cent < 6; i_cent++)
		{

			double renorm = 1;

			name = Form("ChPS_MB_UE_dR%i_cent%i", i_dR, i_cent);
			TH2* h_MB_method = (TH2*)input_file->Get(name.c_str());

			name = Form("ChPS_cone_UE_dR%i_cent%i", i_dR, i_cent);
			TH2* h_cone_method = (TH2*)input_file->Get(name.c_str());

			name = Form("ChPS_TM_UE_dR%i_cent%i", i_dR, i_cent);
			TH2* h_TM_method = (TH2*)input_file->Get(name.c_str());

			name = Form("ChPS_FNS_UE_dR%i_cent%i", i_dR, i_cent);
			TH2* h_FNS_method = (TH2*)input_file->Get(name.c_str());

			name = Form("MB_norm_jet_cent%i", i_cent);
			TH1* h_jet_spectra = (TH1*)((TH1*)input_file->Get(name.c_str()))->Clone(Form("reco_jet_y4_c%i", i_cent));
			h_jet_spectra->Sumw2();

			name = Form("cone_norm_jet_cent%i", i_cent);
			TH1* h_cone_jet_spectra = (TH1*)((TH1*)input_file->Get(name.c_str()))->Clone(Form("cone_UE_norm_y4_c%i", i_cent));
			h_cone_jet_spectra->SetName(Form("%s_mc",name.c_str()));
			h_cone_jet_spectra->Sumw2();

			for (int i_jet_bin = 1; i_jet_bin <= N_jetpt; i_jet_bin++)
			{
				double n_jets = h_jet_spectra->GetBinContent(i_jet_bin);
				double n_jets_cone = h_cone_jet_spectra->GetBinContent(i_jet_bin);

				if (n_jets == 0 || n_jets_cone == 0) continue;

				for (int i_trk_bin = 1; i_trk_bin <= N_jetpt; i_trk_bin++)
				{
					double updated_UE_MB = h_MB_method->GetBinContent(i_trk_bin, i_jet_bin) / n_jets;
					double updated_UE_TM = h_TM_method->GetBinContent(i_trk_bin, i_jet_bin) / n_jets;
					double updated_UE_cone = h_cone_method->GetBinContent(i_trk_bin, i_jet_bin) / n_jets_cone;

					double updated_UE_MB_err = h_MB_method->GetBinError(i_trk_bin, i_jet_bin) / n_jets;
					double updated_UE_TM_err = h_TM_method->GetBinError(i_trk_bin, i_jet_bin) / n_jets;
					double updated_UE_cone_err = h_cone_method->GetBinError(i_trk_bin, i_jet_bin) / n_jets_cone;

					h_MB_method->SetBinContent(i_trk_bin, i_jet_bin, updated_UE_MB);
					h_TM_method->SetBinContent(i_trk_bin, i_jet_bin, updated_UE_TM);
					h_cone_method->SetBinContent(i_trk_bin, i_jet_bin, updated_UE_cone);

					h_MB_method->SetBinError(i_trk_bin, i_jet_bin, updated_UE_MB_err);
					h_TM_method->SetBinError(i_trk_bin, i_jet_bin, updated_UE_TM_err);
					h_cone_method->SetBinError(i_trk_bin, i_jet_bin, updated_UE_cone_err);
				}
			}

			name = Form("ratio_dR%i_cent%i", i_dR, i_cent);
			h_ratio[i_dR][i_cent] = (TH2*)h_TM_method->Clone(name.c_str());
			h_ratio[i_dR][i_cent]->Divide(h_cone_method); //errors propagated assuming they are uncorrelated

			UE_factors->cd();
			name = Form("UE_ratio_dR%i_cent%i", i_dR, i_cent);
			h_ratio[i_dR][i_cent]->Write(name.c_str());

			for (int i_jet_bin = 0; i_jet_bin < N_jetpt; i_jet_bin++)
			{
				name = Form("ratio_1d_dR%i_cent%i_jet%i", i_dR, i_cent, i_jet_bin);
				h_ratio_1d[i_dR][i_cent][i_jet_bin] = (TH1*)h_ratio[i_dR][i_cent]->ProjectionX(name.c_str(), i_jet_bin+1, i_jet_bin+1);
				h_ratio_1d[i_dR][i_cent][i_jet_bin]->GetYaxis()->SetTitle("Correction factors");
				h_ratio_1d[i_dR][i_cent][i_jet_bin]->GetXaxis()->SetTitle("#it{p}_{T}^{trk}");

				for (int i_trk = 0; i_trk < N_trkpt; i_trk++)
				{
					h_ratio_1d_r[i_trk][i_cent][i_jet_bin]->SetBinContent(i_dR+1, h_ratio_1d[i_dR][i_cent][i_jet_bin]->GetBinContent(i_trk+1));
				}

			}

			c1->cd(i_cent+1);
			h_ratio[i_dR][i_cent]->SetMarkerSize(2);
			h_ratio[i_dR][i_cent]->Draw("colz text");
			h_ratio[i_dR][i_cent]->GetYaxis()->SetTitle("#it{p}_{T}^{jet}");
			h_ratio[i_dR][i_cent]->GetXaxis()->SetTitle("#it{p}_{T}^{trk}");

			h_ratio[i_dR][i_cent]->GetXaxis()->SetRangeUser(1,10);
			h_ratio[i_dR][i_cent]->GetYaxis()->SetRangeUser(90,500);
			h_ratio[i_dR][i_cent]->GetZaxis()->SetRangeUser(0,2);
			ltx->SetTextAlign(12);
			ltx->DrawLatexNDC(0.19,0.98,dr_label.c_str());
			ltx->SetTextAlign(32);
			ltx->DrawLatexNDC(0.90,0.98,num_to_cent(31,i_cent).c_str());
			gPad->SetLogx();
			gPad->SetLogy();

		}
		if (i_dR == 0) name = "(";
		else if (i_dR == 12) name = ")";
		else name = "";
		c1->Print(Form("output_pdf%s/PbPb/UE_factors.pdf%s",sys_path.c_str(), name.c_str()), Form("Title: dR%i - %s", i_dR, dr_label.c_str()));


	}


	int jet_pt_start = 7;
	int jet_pt_end = 11;
	TLine *line = new TLine();
	line->SetLineColor(kBlack);

	{
		//just the factors as function of r
		cout << "posres as Function of R" << endl;
		TCanvas *c_pos_res = new TCanvas("c_pos_res","c_pos_res",900,600);

		TLegend *legend_pos_res = new TLegend(0.55, 0.18, 0.75, 0.40, "","brNDC");
		legend_pos_res->SetTextFont(43);
		legend_pos_res->SetBorderSize(0);
		legend_pos_res->SetTextSize(10);



		int jet_itr = 0;
		for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
		{
			string jet_label = Form("%1.0f < #it{p}_{T}^{jet} < %1.0f GeV", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));
			c_pos_res->Clear();
			c_pos_res->Divide(3,2);

			for (int i_cent = 0; i_cent < 6; i_cent++)
			{
				int trk_itr = 0;
				for (int i_trk = 0; i_trk < N_trkpt; i_trk++)
				{
					if (i_trk < 2 || i_trk > 6) continue;

					string trk_label = Form("%1.1f < #it{p}_{T}^{trk} < %1.1f GeV", trkpT_binning->GetBinLowEdge(i_trk+1), trkpT_binning->GetBinUpEdge(i_trk+1));

					SetHStyle_smallify(h_ratio_1d_r[i_trk][i_cent][i_jet], trk_itr, 1);

					if (jet_itr == 0 && i_cent == 0) legend_pos_res->AddEntry(h_ratio_1d_r[i_trk][i_cent][i_jet],trk_label.c_str(),"lp");

					h_ratio_1d_r[i_trk][i_cent][i_jet]->GetYaxis()->SetRangeUser(0.9, 1.1);
					h_ratio_1d_r[i_trk][i_cent][i_jet]->GetXaxis()->SetRangeUser(0., r_max_range);
					h_ratio_1d_r[i_trk][i_cent][i_jet]->GetYaxis()->SetNdivisions(504);
					h_ratio_1d_r[i_trk][i_cent][i_jet]->GetXaxis()->SetTitle("#it{r}");
					h_ratio_1d_r[i_trk][i_cent][i_jet]->GetYaxis()->SetTitle("Correction Factors");

					c_pos_res->cd(i_cent+1);
					if (trk_itr == 0) h_ratio_1d_r[i_trk][i_cent][i_jet]->Draw("hist text");
//					else h_ratio_1d_r[i_trk][i_cent][i_jet]->Draw("same p");
					line->DrawLine(0, 1, r_max_range, 1);

					trk_itr++;
				}

				c_pos_res->cd(i_cent+1);
				ltx->SetTextAlign(32);
				ltx->DrawLatexNDC(0.94,0.92,num_to_cent(31,i_cent).c_str());
				ltx->DrawLatexNDC(0.94,0.85,jet_label.c_str());
				legend_pos_res->Draw();

			}

			if (i_jet == jet_pt_start) name = "(";
			else if (i_jet == jet_pt_end - 1) name = ")";
			else name = "";
			c_pos_res->Print(Form("output_pdf%s/PbPb/UE_factors_r.pdf%s",sys_path.c_str(), name.c_str()), Form("Title: jet%i", i_jet));
			jet_itr++;

		}
	}

	cout << "######### DONE GETTING UE FACTORS #########" << endl << endl;;

}
