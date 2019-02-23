#include "functions/global_variables.h"

void comp_jet()
{
	SetAtlasStyle();
	TFile *f_nom = new TFile("raw_results/nominal/FF_data_out_histo_PbPb_5p02_r001.root");
	TFile *f_no_iso = new TFile("raw_results/config-1/FF_data_out_histo_PbPb_5p02_r001.root");
	TFile *f_iso_12 = new TFile("raw_results/config-2/FF_data_out_histo_PbPb_5p02_r001.root");



	TCanvas *c = new TCanvas();
	TLegend *legend = new TLegend();
	legend->SetBorderSize(0);
	legend->SetTextFont(43);
	legend->SetTextSize(16);

	TLatex *ltx = new TLatex();
	ltx->SetTextFont(43);
	ltx->SetTextSize(18);
	c->Divide(3,2);

	for (int i = 0; i < 6; i++)
	{
		string name = Form("h_reco_jet_spectrum_y4_cent%i", i);


		TH1 *h_nom = (TH1*)(TH1*)f_nom->Get(name.c_str())->Clone(Form("%s_nom", name.c_str()));
		h_nom->SetName("h_nom");

		TH1 *h_no_iso = (TH1*)(TH1*)f_no_iso->Get(name.c_str())->Clone(Form("%s_no_iso", name.c_str()));
		h_no_iso->SetName("h_no_iso");

		TH1 *h_iso_12 = (TH1*)(TH1*)f_iso_12->Get(name.c_str())->Clone(Form("%s_iso_12", name.c_str()));
		h_iso_12->SetName("h_iso_12");


		h_no_iso->Divide(h_nom);
		h_iso_12->Divide(h_nom);

		h_no_iso->Scale(1E0);
		h_iso_12->Scale(2);

		if (i== 0)
		{
//			legend->AddEntry(h_iso_12,"iso 1.2 (#times 2)", "l");
			legend->AddEntry(h_no_iso,"no isolation", "l");
		}

		SetHStyle_smallify(h_no_iso, 0, 1);
		SetHStyle_smallify(h_iso_12, 1, 1);

		c->cd(i+1);


		h_no_iso->SetMarkerSize(3.);
		h_iso_12->SetMarkerSize(3.);

		h_no_iso->GetYaxis()->SetTitle("No isolation/Isolated");
		h_no_iso->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");


		h_no_iso->GetYaxis()->SetRangeUser(0,2);
		h_no_iso->GetYaxis()->SetNdivisions(505);
		h_no_iso->GetXaxis()->SetRangeUser(126, 316);
		h_no_iso->GetXaxis()->SetNdivisions(505);

		h_iso_12->GetYaxis()->SetRangeUser(0,3);
		h_iso_12->GetYaxis()->SetNdivisions(505);
		h_iso_12->GetXaxis()->SetRangeUser(126, 316);
		h_iso_12->GetXaxis()->SetNdivisions(505);

		h_no_iso->DrawCopy("hist text");
//		h_iso_12->DrawCopy("hist text same");

//		h_iso_12->DrawCopy("hist text");
//		gPad->SetLogy();

		ltx->DrawLatexNDC(0.19, 0.88, num_to_cent(31,i).c_str());
		legend->Draw();
	}


	c->Print("tmp.pdf");

}



