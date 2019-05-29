#include "output_dev/functions/global_variables.h"

void runStats(int config = 0)
{
	SetAtlasStyle();
	gErrorIgnoreLevel = 3001;
	double r_max_range = 0.8;
	string sys_path = "";
	if (config == 0) sys_path = Form("nominal");
	if (config > 0 && config < 100) sys_path = Form("sys%i", config);
	if (config > 100) sys_path = Form("c%i", config);

	TFile *input_file = new TFile(Form("output_dev/raw_results/%s/FF_MC_out_histo_PbPb_5p02_r001.root", sys_path.c_str()));
	TFile *input_file_data = new TFile(Form("output_dev/raw_results/%s/FF_data_out_histo_PbPb_5p02_r001.root", sys_path.c_str()));

	TAxis* dR_binning = (TAxis*)((TH3*)input_file->Get("h_dR_change_jetpt0_cent0"))->GetXaxis();
	TAxis* jetpT_binning = (TAxis*)((TH3*)input_file->Get("ChPS_raw_0_dR0_cent0"))->GetYaxis();
	TAxis* trkpT_binning = (TAxis*)((TH3*)input_file->Get("ChPS_raw_0_dR0_cent0"))->GetXaxis();
	TAxis* run_binning = (TAxis*)((TH3*)input_file->Get("h_event_rN"))->GetXaxis();

	int N_dR = dR_binning->GetNbins();
	int N_jetpt = jetpT_binning->GetNbins();
	int N_trkpt = trkpT_binning->GetNbins();
	int N_runs = run_binning->GetNbins();

//	int jet_pt_start = jetpT_binning->FindBin(127);
//	int jet_pt_end = jetpT_binning->FindBin(130);
//	int trk_pt_start = trkpT_binning->FindBin(1);
//	int trk_pt_end = trkpT_binning->FindBin(1.5);
	int run_start = 1;
	int run_end = N_runs;

	map<int, double> luminosity = {{ 286665 , 0.02071 }, { 286711 , 0.419829 }, { 286717 , 0.590276 }, { 286748 , 4.24438 }, { 286767 , 5.799 }, { 286834 , 13.1028 }, { 286854 , 13.595759 }, { 286908 , 11.00192 }, { 286990 , 10.27951 }, { 287038 , 15.9825 }, { 287044 , 23.3479 }, { 287068 , 6.83118 }, { 287222 , 1.33018 }, { 287224 , 1.932658 }, { 287259 , 17.4035 }, { 287270 , 22.0425 }, { 287281 , 24.1107 }, { 287321 , 5.49172 }, { 287330 , 21.620117 }, { 287334 , 16.3135 }, { 287378 , 16.81257 }, { 287380 , 0.319641 }, { 287382 , 17.651474 }, { 287560 , 0.572368 }, { 287594 , 12.6278 }, { 287632 , 18.9915 }, { 287706 , 26.6057 }, { 287728 , 25.93 }, { 287827 , 24.1147 }, { 287843 , 25.0225 }, { 287866 , 42.1866 }, { 287924 , 22.5426 }, { 287931 , 37.2019 }};

	TH1* h_tmp_mc = (TH1*)input_file->Get("h_event_rN");
	h_tmp_mc->SetName("h_event_rN_mc");

	TH1* h_tmp_data = (TH1*)input_file_data->Get("h_event_rN");
	h_tmp_data->SetName("h_event_rN_data");

	TH1* h_eventPercentage_mc = new TH1D("h_eventPercentage_mc","h_eventPercentage_mc",N_runs,0,N_runs);
	h_eventPercentage_mc->GetXaxis()->SetTitle("Run Number");
	h_eventPercentage_mc->GetYaxis()->SetTitle("Event Fraction");
	h_eventPercentage_mc->GetYaxis()->SetNdivisions(504);

	TH1* h_eventPercentage_data = new TH1D("h_eventPercentage_data","h_eventPercentage_data",N_runs,0,N_runs);
	h_eventPercentage_data->GetXaxis()->SetTitle("Run Number");
	h_eventPercentage_data->GetYaxis()->SetTitle("Event Fraction");
	h_eventPercentage_data->GetYaxis()->SetNdivisions(504);

	for (int i = 1; i <= h_eventPercentage_data->GetNbinsX(); i=i+1)
	{
		h_eventPercentage_mc->GetXaxis()->SetBinLabel(i,Form("%1.0f",h_tmp_mc->GetBinLowEdge(i)));
		h_eventPercentage_data->GetXaxis()->SetBinLabel(i,Form("%1.0f",h_tmp_data->GetBinLowEdge(i)));
	}

	for (int i = 1; i <= h_eventPercentage_mc->GetNbinsX(); i++)
	{
		h_eventPercentage_mc->SetBinContent(i,h_tmp_mc->GetBinContent(i));///luminosity[h_tmp_mc->GetBinLowEdge(i)]);
		cout << luminosity[h_tmp_mc->GetBinLowEdge(i)] << " " << h_tmp_mc->GetBinLowEdge(i) << endl;
		h_eventPercentage_data->SetBinContent(i,h_tmp_data->GetBinContent(i));///luminosity[h_tmp_mc->GetBinLowEdge(i)]);

	}

	h_eventPercentage_mc->Scale(1./h_eventPercentage_mc->Integral());
	h_eventPercentage_data->Scale(1./h_eventPercentage_data->Integral());

	TCanvas *c_x = new TCanvas("c_x","c_x",800,600);
	TLegend *legend_x = new TLegend(0.20, 0.6, 0.30, 0.7, "","brNDC");
	legend_x->SetTextFont(43);
	legend_x->SetBorderSize(0);
	legend_x->SetTextSize(24);
	legend_x->SetNColumns(1);
	TLatex *ltx = new TLatex();
	ltx->SetTextFont(43);
	ltx->SetTextSize(12);
	ltx->SetTextAlign(32);

	TLine *line = new TLine();

	c_x->cd();
	SetHStyle_smallify(h_eventPercentage_data, 0, 0);
	SetHStyle_smallify(h_eventPercentage_mc, 1, 0);
	h_eventPercentage_mc->GetXaxis()->SetLabelSize(14);
	h_eventPercentage_mc->LabelsOption("v");
	h_eventPercentage_data->LabelsOption("v");
//	h_eventPercentage_mc->GetYaxis()->SetTitleOffset(0.8);
	h_eventPercentage_mc->GetXaxis()->SetTitleOffset(1.4);
	h_eventPercentage_mc->GetYaxis()->SetRangeUser(0,0.15);
	h_eventPercentage_mc->Draw("hist");
	h_eventPercentage_data->Draw("hist same");
	legend_x->AddEntry(h_eventPercentage_data,"Data","lp");
	legend_x->AddEntry(h_eventPercentage_mc,"MC","lp");
	legend_x->Draw();

	c_x->Print(Form("misc_plots/EventPercentages_c%i.pdf",config));

}
