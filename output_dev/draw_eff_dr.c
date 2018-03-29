#include "functions/global_variables.h"
void draw_eff_dr()
{
	gErrorIgnoreLevel = 3001;
	SetAtlasStyle();

	string name;

	//	##############	Reading config	##############"
	TEnv *m_config = new TEnv();
	m_config->ReadFile("perf_config.cfg", EEnvLevel(1));
	m_config->Print();
	
	std::string dataset_type = "PbPb"; dataset_type = m_config->GetValue("dataset_type", dataset_type.c_str());
	std::string tracking_cut = "ppTight"; tracking_cut = m_config->GetValue("tracking_cut", tracking_cut.c_str());
	int centrality_scheme = 31; centrality_scheme = m_config->GetValue("centrality_scheme", centrality_scheme);
	//	##############	Config done	##############"


	//reading from raw efficiency file
	name = Form("mc_eff_%s_trketa_dr_jetptinc_%s.root",dataset_type.c_str(), tracking_cut.c_str());
	TFile *input_file = new TFile(name.c_str());

	TAxis* trk_pt_binning = (TAxis*)input_file->Get("trk_pt_binning");
	TAxis* trk_eta_binning_new = (TAxis*)input_file->Get("new_trk_eta_binning");
	TAxis* dR_binning = (TAxis*)input_file->Get("dR_binning");

	int n_trk_pt_bins = trk_pt_binning->GetNbins();
	int n_trk_eta_bins_new = trk_eta_binning_new->GetNbins();
	int n_dR = dR_binning->GetNbins();


	TCanvas *canvas1 = new TCanvas("C1", "C1",900,600);

	TLine *line = new TLine();
	TLatex *ltx = new TLatex();
	ltx->SetNDC();
	ltx->SetTextFont(43);
	ltx->SetTextSize(15);
	ltx->SetTextAlign(11);

	TLegend *legend = new TLegend(0.20,0.17,0.80,0.40,NULL,"brNDC");
	legend->SetBorderSize(0);
	legend->SetNColumns(2);
	legend->SetTextFont(43);
	legend->SetTextSize(12);
	if (dataset_type == "PbPb") legend->SetTextSize(10);

	vector<vector<vector<TH1*>>> h_efficiency(n_cent_cuts, vector<vector<TH1*>> (n_trk_eta_bins_new, vector<TH1*> (n_dR)));
	vector<vector<vector<TGraphAsymmErrors*>>> g_efficiency(n_cent_cuts, vector<vector<TGraphAsymmErrors*>> (n_trk_eta_bins_new, vector<TGraphAsymmErrors*> (n_dR)));


	bool doSmall;
	if (dataset_type == "pp") doSmall = false;
	if (dataset_type == "PbPb") doSmall = true;

	bool first_pass = 1;

	for (int i_eta_cuts = 0; i_eta_cuts < n_trk_eta_bins_new; i_eta_cuts++)
	{
		double eta_lo = trk_eta_binning_new->GetBinLowEdge(i_eta_cuts+1);
		double eta_hi = trk_eta_binning_new->GetBinUpEdge(i_eta_cuts+1);

//		if ((eta_hi < -1. || eta_lo > 1.)) continue;

		canvas1->cd();
		canvas1->Clear();
		if (dataset_type == "PbPb") canvas1->Divide(3,2);

		for (int i_cent_cuts = 0; i_cent_cuts < n_cent_cuts; i_cent_cuts++)
		{
			if (dataset_type == "PbPb" && i_cent_cuts == 6) continue;
			if (dataset_type == "pp" && i_cent_cuts < 6) continue;

			string centrality = num_to_cent(centrality_scheme,i_cent_cuts);


			int style = 0;
			for (int i_dR = 0; i_dR < n_dR; i_dR++)
			{
				double dr_lo = dR_binning->GetBinLowEdge(i_dR+1);
				double dr_hi = dR_binning->GetBinUpEdge(i_dR+1);

				if (i_dR % 4 !=0) continue;

				name = Form("hist_eff_eta%i_cent%i_dr%i",i_eta_cuts,i_cent_cuts, i_dR);
				h_efficiency.at(i_cent_cuts).at(i_eta_cuts).at(i_dR) = (TH1*)input_file->Get(name.c_str());

				SetHStyle_smallify(h_efficiency.at(i_cent_cuts).at(i_eta_cuts).at(i_dR),style++, doSmall);

				h_efficiency.at(i_cent_cuts).at(i_eta_cuts).at(i_dR)->GetYaxis()->SetRangeUser(0.5,1.1);
				h_efficiency.at(i_cent_cuts).at(i_eta_cuts).at(i_dR)->GetXaxis()->SetRangeUser(1,200);
				h_efficiency.at(i_cent_cuts).at(i_eta_cuts).at(i_dR)->GetYaxis()->SetTitle("Efficiency");
				h_efficiency.at(i_cent_cuts).at(i_eta_cuts).at(i_dR)->GetXaxis()->SetTitle("#it{p}_{T}^{truth} [GeV]");
				h_efficiency.at(i_cent_cuts).at(i_eta_cuts).at(i_dR)->SetTitle(Form("Efficiency: %s, %4.2f < #eta < %4.2f, %2.2f < r < %2.2f",centrality.c_str(), eta_lo, eta_hi, dr_lo, dr_hi));

				if (first_pass) legend->AddEntry(h_efficiency.at(i_cent_cuts).at(i_eta_cuts).at(i_dR),Form("%4.2f < r < %4.2f",dr_lo, dr_hi),"lp");

				
				if (dataset_type == "PbPb") canvas1->cd(i_cent_cuts+1);
				else canvas1->cd();

				if (style == 0) h_efficiency.at(i_cent_cuts).at(i_eta_cuts).at(i_dR)->Draw("a p");
				else h_efficiency.at(i_cent_cuts).at(i_eta_cuts).at(i_dR)->Draw("same p");
				gPad->SetLogx();


			}
			ltx->DrawLatex(0.19,0.88,centrality.c_str());
			ltx->DrawLatex(0.19,0.80,Form("%4.2f < #eta < %4.2f", eta_lo, eta_hi));

//			if (dataset_type == "PbPb") canvas1->cd(1);
//			else canvas1->cd();
			legend->Draw();

			first_pass = 0;
		}

		string pdf = "";
//		if (i_eta_cuts == 7) pdf = "(";
//		else if (i_eta_cuts == 14) pdf = ")";
		if (i_eta_cuts == 0) pdf = "(";
		else if (i_eta_cuts == n_trk_eta_bins_new-1) pdf = ")";
		canvas1->Print(Form("output_pdf/eff_cent_trketa_dr_%s_%s.pdf%s", dataset_type.c_str(), tracking_cut.c_str(), pdf.c_str()),Form("Title: eta%i", i_eta_cuts));
	}

}

