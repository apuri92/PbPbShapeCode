#include "output_dev/functions/global_variables.h"

void draw_integrated_UE(int config)
{
	SetAtlasStyle();
	gStyle->SetPalette(1);
	TFile *file = new TFile(Form("UE_c%i.root",config));

	//integrated over jet cone
	TCanvas *c1 = new TCanvas("c1","c1",700,800);
	c1->Divide(2,3);

	TCanvas *c2 = new TCanvas("c2","c2",900,800);
	c2->Divide(2,3);

	string name;
	TLatex *ltx = new TLatex();
	ltx->SetTextFont(43);
	ltx->SetTextSize(12);
	ltx->SetTextAlign(12);

	
	for (int i_cent = 0; i_cent < 6; i_cent++)
	{
		name = Form("Cone_Cent%i",i_cent);
		TH2* h_cone_MC = (TH2*)file->Get(name.c_str());
		h_cone_MC->SetName(Form("%s_mc",name.c_str()));

		name = Form("Cone_data_Cent%i",i_cent);
		TH2* h_cone_data = (TH2*)file->Get(name.c_str());
		h_cone_data->SetName(Form("%s_data",name.c_str()));

		name = Form("MC_data_Cent%i",i_cent);
		TH2* h_MB_data = (TH2*)file->Get(name.c_str());
		h_MB_data->SetName(Form("%s_data",name.c_str()));

		name = Form("MC_Cent%i",i_cent);
		TH2* h_MB_mc = (TH2*)file->Get(name.c_str());
		h_MB_mc->SetName(Form("%s_mc",name.c_str()));

		name = Form("TM_Cent%i",i_cent);
		TH2* h_TM = (TH2*)file->Get(name.c_str());
		h_TM->SetName(Form("%s_mc",name.c_str()));

		TH2* h_ratio_cone_data_mb_data = (TH2*)h_cone_data->Clone(Form("ratio_cone_data_mb_data_c%i",i_cent));
		h_ratio_cone_data_mb_data->Divide(h_MB_data);

		TH2* h_ratio_cone_mc_tm = (TH2*)h_cone_MC->Clone(Form("ratio_cone_mc_tm_c%i",i_cent));
		h_ratio_cone_mc_tm->Divide(h_TM);

		h_ratio_cone_data_mb_data->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]");
		h_ratio_cone_data_mb_data->GetYaxis()->SetTitle("p_{T}^{Jet} [GeV]");
		h_ratio_cone_data_mb_data->GetZaxis()->SetTitle("UE_{Cone}^{Data} / UE_{MC}^{Data}");

		h_ratio_cone_mc_tm->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]");
		h_ratio_cone_mc_tm->GetYaxis()->SetTitle("p_{T}^{Jet} [GeV]");
		h_ratio_cone_mc_tm->GetZaxis()->SetTitle("UE_{Cone}^{MC} / UE_{TM}^{MC}");

		h_ratio_cone_data_mb_data->GetXaxis()->SetMoreLogLabels(kTRUE);
		h_ratio_cone_data_mb_data->GetXaxis()->SetRangeUser(1,10);
		h_ratio_cone_data_mb_data->GetYaxis()->SetRangeUser(128,350);
		h_ratio_cone_data_mb_data->SetMarkerSize(1.5);

		h_ratio_cone_mc_tm->GetXaxis()->SetMoreLogLabels(kTRUE);
		h_ratio_cone_mc_tm->GetXaxis()->SetRangeUser(1,10);
		h_ratio_cone_mc_tm->GetYaxis()->SetRangeUser(128,350);
		h_ratio_cone_mc_tm->SetMarkerSize(1.5);

//		h_ratio->GetZaxis()->SetNdivisions(504);

		double lo_range = 0, hi_range = 2.;
		if (i_cent == 0){ lo_range = 0.95 ; hi_range = 1.01;}
		if (i_cent == 1){ lo_range = 0.96 ; hi_range = 1.03;}
		if (i_cent == 2){ lo_range = 0.975 ; hi_range = 1.04;}
		if (i_cent == 3){ lo_range = 0.99 ; hi_range = 1.15;}
		if (i_cent == 4){ lo_range = 0.95 ; hi_range = 1.1;}
		if (i_cent == 5){ lo_range = 1. ; hi_range = 1.36;}

//		h_ratio->GetZaxis()->SetRangeUser(lo_range,hi_range);


		c1->cd(i_cent+1);
		h_ratio_cone_data_mb_data->DrawCopy("colz text");
		gPad->SetLogx();
		gPad->SetRightMargin(0.18);
		ltx->SetTextAlign(32);
		ltx->DrawLatexNDC(0.90,0.98,num_to_cent(31,i_cent).c_str());

		c2->cd(i_cent+1);
		h_ratio_cone_mc_tm->DrawCopy("colz text");
		gPad->SetLogx();
		gPad->SetRightMargin(0.18);
		ltx->SetTextAlign(32);
		ltx->DrawLatexNDC(0.90,0.98,num_to_cent(31,i_cent).c_str());

	}

	c1->Print(Form("UE_cone_mb_comparison_c%i.pdf",config));
	c2->Print(Form("UE_cone_tm_comparison_c%i.pdf",config));




}
