#include "functions/global_variables.h"
void draw_eff_trketa()
{
	gErrorIgnoreLevel = 3001;
	SetAtlasStyle();

	string name;

	//	##############	Reading config	##############"
	TEnv *m_config = new TEnv();
	m_config->ReadFile("perf_config.cfg", EEnvLevel(1));
	m_config->Print();
	
	int config = 0; config = m_config->GetValue("config", config);
	std::string dataset_type = "PbPb"; dataset_type = m_config->GetValue("dataset_type", dataset_type.c_str());
	std::string tracking_cut = "ppTight"; tracking_cut = m_config->GetValue("tracking_cut", tracking_cut.c_str());
	int centrality_scheme = 31; centrality_scheme = m_config->GetValue("centrality_scheme", centrality_scheme);
	//	##############	Config done	##############"


	//reading from raw efficiency file
	name = Form("perf_output_c%i/mc_eff_%s_trketa_jetptinc_%s.root",config, dataset_type.c_str(), tracking_cut.c_str());
	TFile *input_file = new TFile(name.c_str());

//	TAxis* trk_pt_binning = (TAxis*)input_file->Get("trk_pt_binning");
	TAxis* trk_eta_binning_new = (TAxis*)input_file->Get("new_trk_eta_binning");

//    int n_trk_pt_bins = trk_pt_binning->GetNbins();
    int n_trk_eta_bins_new = trk_eta_binning_new->GetNbins();


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

	vector<vector<TH1*>> h_efficiency(n_cent_cuts, vector<TH1*> (n_trk_eta_bins_new));
	vector<vector<TGraphAsymmErrors*>> g_efficiency(n_cent_cuts, vector<TGraphAsymmErrors*> (n_trk_eta_bins_new));

	canvas1->cd();
	canvas1->Clear();
	if (dataset_type == "PbPb") canvas1->Divide(3,2);

	bool doSmall;
	if (dataset_type == "pp") doSmall = false;
	if (dataset_type == "PbPb") doSmall = true;

	bool first_pass = 1;
	for (int i_cent_cuts = 0; i_cent_cuts < n_cent_cuts; i_cent_cuts++)
	{
		if (dataset_type == "PbPb" && i_cent_cuts == 6) continue;
		if (dataset_type == "pp" && i_cent_cuts < 6) continue;

		string centrality = num_to_cent(centrality_scheme,i_cent_cuts);

		int style = 0;
		for (int i_eta_cuts = 0; i_eta_cuts < n_trk_eta_bins_new; i_eta_cuts++)
		{
			double eta_lo = trk_eta_binning_new->GetBinLowEdge(i_eta_cuts+1);
			double eta_hi = trk_eta_binning_new->GetBinUpEdge(i_eta_cuts+1);

			if ((eta_hi < -1. || eta_lo > 1.)) continue;

			name = Form("hist_eff_eta%i_cent%i", i_eta_cuts, i_cent_cuts);
			h_efficiency.at(i_cent_cuts).at(i_eta_cuts) = (TH1*)input_file->Get(name.c_str());


            SetHStyle_smallify(h_efficiency.at(i_cent_cuts).at(i_eta_cuts),style++, doSmall);

            h_efficiency.at(i_cent_cuts).at(i_eta_cuts)->GetYaxis()->SetRangeUser(0.5,1);
            h_efficiency.at(i_cent_cuts).at(i_eta_cuts)->GetXaxis()->SetRangeUser(1,200);
            h_efficiency.at(i_cent_cuts).at(i_eta_cuts)->GetYaxis()->SetTitle("Efficiency");
            h_efficiency.at(i_cent_cuts).at(i_eta_cuts)->GetXaxis()->SetTitle("#it{p}_{T}^{truth} [GeV]");
            h_efficiency.at(i_cent_cuts).at(i_eta_cuts)->SetTitle(Form("Efficiency: %s, %4.2f < #eta < %4.2f",centrality.c_str(), eta_lo, eta_hi));

			if (first_pass) legend->AddEntry(h_efficiency.at(i_cent_cuts).at(i_eta_cuts),Form("%4.2f < #eta < %4.2f",eta_lo, eta_hi),"lp");

			if (dataset_type == "PbPb") canvas1->cd(i_cent_cuts+1);
			else canvas1->cd();

			if (style == 0) h_efficiency.at(i_cent_cuts).at(i_eta_cuts)->Draw("a p");
			else h_efficiency.at(i_cent_cuts).at(i_eta_cuts)->Draw("same p");
			gPad->SetLogx();

		}
		name = Form("%s",centrality.c_str());
		ltx->SetTextAlign(11);
		ltx->DrawLatex(0.19,0.88,name.c_str());

        
//        if (i_cent_cuts == 9) name = "(";
//        else if (i_cent_cuts == 9 ) name = ")";
//        else name = "";
		first_pass = 0;

	}
	if (dataset_type == "PbPb") canvas1->cd(1);
	else canvas1->cd();
	legend->Draw();

	name = Form("perf_output_c%i/eff_cent_trketa_%s_%s.pdf",config, dataset_type.c_str(), tracking_cut.c_str());
	canvas1->Print(name.c_str());


}

