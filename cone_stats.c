#include "output_dev/functions/global_variables.h"
void cone_stats()
{
	SetAtlasStyle();
	//c992 contains cone stats for each centrality
	TFile *f_mc = new TFile(Form("output_dev/raw_results/c992/FF_MC_out_histo_PbPb_5p02_r001.root"));
	TFile *f_data = new TFile(Form("output_dev/raw_results/c992/FF_data_out_histo_PbPb_5p02_r001.root"));

	TCanvas *c = new TCanvas("c","c",900,600);
	TLegend *legend = new TLegend(0.19, 0.6, 0.4, 0.7, "","brNDC");
	legend->SetTextFont(43);
	legend->SetBorderSize(0);
	legend->SetTextSize(15);
	legend->SetNColumns(1);
	TLatex *ltx = new TLatex();
	ltx->SetTextFont(43);
	ltx->SetTextSize(18);
	ltx->SetTextAlign(12);
	c->Divide(3,2);

	for (int i = 0; i < 6; i++)
	{
		TH1* h_mc = (TH1*)f_mc->Get(Form("h_tmp_cone_stats_cent%i",i));
		h_mc->SetName(Form("MC_ConeStats_c%i",i));

		TH1* h_data = (TH1*)f_data->Get(Form("h_tmp_cone_stats_cent%i",i));
		h_data->SetName(Form("data_ConeStats_c%i",i));

		SetHStyle_open_smallify(h_mc,0,1);
		SetHStyle_smallify(h_data,0,1);
		if (i == 0)
		{
			legend->AddEntry(h_data,"Data","p");
			legend->AddEntry(h_mc,"MC","p");
		}
		h_mc->Scale(1./h_mc->Integral());
		h_data->Scale(1./h_data->Integral());

		c->cd(i+1);
		h_mc->GetXaxis()->SetRangeUser(0,10);
		h_mc->GetYaxis()->SetRangeUser(0,0.5);
		h_mc->DrawCopy(" hist text");
		h_data->DrawCopy("same");
		legend->Draw();
		ltx->DrawLatexNDC(0.22,0.87,Form("%s",num_to_cent(31,i).c_str()));

	}
	c->Print(Form("misc_conf_plots/cone_stats.pdf"));

}













