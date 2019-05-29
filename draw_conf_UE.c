#include "output_dev/functions/global_variables.h"

void draw_maps(int dR = 0, int i_trk = 2)
{
	TFile *file = new TFile("pPbFragmentation/data/UE_MC_comb.root");
	gErrorIgnoreLevel = 3001;
	TCanvas *c = new TCanvas("c","c",900,600);
	c->Divide(3,2);

	TFile *axis_files = new TFile(Form("output_dev/raw_results/nominal/FF_data_out_histo_PbPb_5p02_r001.root"));
	TAxis* dR_binning = (TAxis*)((TH3*)axis_files->Get("h_dR_change_jetpt0_cent0"))->GetXaxis();
	TAxis* jetpT_binning = (TAxis*)((TH3*)axis_files->Get("ChPS_raw_0_dR0_cent0"))->GetYaxis();
	TAxis* trkpT_binning = (TAxis*)((TH3*)axis_files->Get("ChPS_raw_0_dR0_cent0"))->GetXaxis();

	string name = "h_UE_new_MC_dR0_dPsi0_pt2_cent0_jet7";
	TLatex *ltx = new TLatex();
	ltx->SetTextFont(43);
	ltx->SetTextSize(13);
	ltx->SetTextAlign(12);

	for (int i = 0; i < 6; i++)
	{
		name = Form("h_UE_new_MC_dR%i_dPsi4_pt%i_cent%i_jet7",dR, i_trk, i);
		TH1* h0 = (TH1*)file->Get(name.c_str());

		h0->GetXaxis()->SetRangeUser(-1.69,1.7);
		h0->GetYaxis()->SetRangeUser(-3.14,3.14);
		h0->GetXaxis()->SetTitle("#eta^{Jet}");
		h0->GetYaxis()->SetTitle("#phi^{Jet}");
		c->cd(i+1);
		gPad->SetRightMargin(0.14);
		gPad->SetTopMargin(0.06);
		h0->DrawCopy("colz");
		string dRlabel = Form("%1.2f < r < %1.2f", dR_binning->GetBinLowEdge(dR+1),dR_binning->GetBinUpEdge(dR+1));
		string trklabel = Form("%1.0f < p_{T}^{trk} < %1.0f GeV", trkpT_binning->GetBinLowEdge(i_trk),trkpT_binning->GetBinUpEdge(i_trk));
		ltx->DrawLatexNDC(0.15,0.97,Form("%s, %s, %s",num_to_cent(31,i).c_str(), trklabel.c_str(), dRlabel.c_str()));

		delete h0;
	}

	c->Print(Form("misc_plots/eta_phi_map_trk%i_dR%i.pdf",i_trk, dR));
	delete c;
}

void draw_map_stat()
{
	string name = Form("./output_dev/unfold/UE_MapSystematic.root");
	TFile *file = new TFile(name.c_str());
	TCanvas *c = new TCanvas("c","c",800,600);

	TFile *axis_files = new TFile(Form("output_dev/raw_results/nominal/FF_data_out_histo_PbPb_5p02_r001.root"));
	TAxis* dR_binning = (TAxis*)((TH3*)axis_files->Get("h_dR_change_jetpt0_cent0"))->GetXaxis();
	TAxis* jetpT_binning = (TAxis*)((TH3*)axis_files->Get("ChPS_raw_0_dR0_cent0"))->GetYaxis();
	TAxis* trkpT_binning = (TAxis*)((TH3*)axis_files->Get("ChPS_raw_0_dR0_cent0"))->GetXaxis();

	TLatex *ltx = new TLatex();
	ltx->SetTextFont(43);
	ltx->SetTextSize(13);
	ltx->SetTextAlign(12);

	int trk_arr[4] = {3,5,6,7};
	int i_dR = 0, i_cent = 0, i_jet = 8;
	TLegend *legend_y = new TLegend(0.20, 0.7, 0.5, 0.9, "","brNDC");
	legend_y->SetTextFont(43);
	legend_y->SetBorderSize(0);
	legend_y->SetTextSize(14);
	legend_y->SetNColumns(1);

	for (int i = 0; i<4; i++)
	{
		name = Form("ChPS_MB_UE_dR%i_cent%i_bins_trk%i_jet%i", i_dR, i_cent, trk_arr[i], i_jet);
		TH1* h0 = (TH1*)file->Get(name.c_str());
		h0->GetXaxis()->SetRangeUser(-0.06,0.06);
		h0->GetYaxis()->SetRangeUser(0,80);
		h0->GetXaxis()->SetTitle("Relative Error");
		h0->GetYaxis()->SetTitle("Counts");
		SetHStyle_smallify(h0, i, 0);
		TF1* f0 = (TF1*)h0->GetFunction("fit");
		f0->SetLineColor(h0->GetMarkerColor());
		legend_y->AddEntry(h0,h0->GetTitle(),"l");
		for (int j = 0; j < h0->GetNbinsX(); j++)
		{
			h0->SetBinError(j+1,0.00001);
		}

		if (i == 0) h0->DrawCopy("func p");
		else h0->DrawCopy("same func p");
//		delete h0;
	}
	legend_y->Draw();
	c->Print("misc_plots/map_stat_gaus.pdf");

	int dR = 3;

	c->cd();
	c->Clear();
	name = Form("ChPS_MB_UE_dR%i_cent0_statSys", dR);
	TH1* h0 = (TH1*)file->Get(name.c_str());
	h0->Draw("colz text");

	gPad->SetLogx();
	gPad->SetLogy();
	h0->GetXaxis()->SetRangeUser(1,10);
	h0->GetYaxis()->SetRangeUser(100,390);
	h0->GetZaxis()->SetRangeUser(0.95,1.05);
	h0->GetXaxis()->SetTitle("#it{p}_{T}^{trk} [GeV]");
	h0->GetYaxis()->SetTitle("#it{p}_{T}^{Jet} [GeV]");
	h0->GetYaxis()->SetMoreLogLabels(kTRUE);
	h0->GetYaxis()->SetLabelSize(0.04);
	string dRlabel = Form("%1.2f < r < %1.2f", dR_binning->GetBinLowEdge(dR+1),dR_binning->GetBinUpEdge(dR+1));
	ltx->DrawLatexNDC(0.15,0.97,Form("%s, %s",num_to_cent(31,0).c_str(), dRlabel.c_str()));
	c->Print("misc_plots/map_stat_size.pdf");


	delete c;
}

void draw_conf_UE()
{
	SetAtlasStyle();
	int dR = 5;
	int i_trk = 2;
//	draw_maps(0, 6); draw_maps(5, 6); draw_maps(9, 6);
//	draw_maps(0, 2); draw_maps(5, 2); draw_maps(9, 2);
	draw_map_stat();

}
