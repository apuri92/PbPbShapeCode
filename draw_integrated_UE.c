#include "output_dev/functions/global_variables.h"
#include "get_madness_factor.c"

void draw_integrated_UE(int config = 24)
{
	SetAtlasStyle();
	gStyle->SetPalette(1);
	TFile *file = new TFile(Form("UE_c%i.root",config));

	//integrated over jet cone
	TCanvas *c1 = new TCanvas("c1","c1",700,900);
	c1->Divide(2,3);

	TCanvas *c2 = new TCanvas("c2","c2",900,600);
	c2->Divide(3,2);

	TCanvas *c3 = new TCanvas("c3","c3",600,900);
	c3->Divide(2,3);

	string name;
	TLatex *ltx = new TLatex();
	ltx->SetTextFont(43);
	ltx->SetTextSize(12);
	ltx->SetTextAlign(12);


	TLegend *legend_x_jet = new TLegend(0.20, 0.2, 0.5, 0.4, "","brNDC");
	legend_x_jet->SetTextFont(43);
	legend_x_jet->SetBorderSize(0);
	legend_x_jet->SetTextSize(12);
	legend_x_jet->SetNColumns(1);

	for (int i_cent = 0; i_cent < 6; i_cent++)
	{
		name = Form("Cone_Cent%i",i_cent);
		TH2* h_cone_MC = (TH2*)file->Get(name.c_str());
		h_cone_MC->SetName(Form("%s_mc",name.c_str()));

		name = Form("Cone_data_Cent%i",i_cent); //this is not corrected
		TH2* h_cone_data = (TH2*)file->Get(name.c_str());
		h_cone_data->SetName(Form("%s_data",name.c_str()));

		name = Form("MC_Cent%i",i_cent);
		TH2* h_MB_mc = (TH2*)file->Get(name.c_str());
		h_MB_mc->SetName(Form("%s_mc",name.c_str()));

		name = Form("MC_data_Cent%i",i_cent);
		TH2* h_MB_data = (TH2*)file->Get(name.c_str());
		h_MB_data->SetName(Form("%s_data",name.c_str()));

		name = Form("TM_Cent%i",i_cent);
		TH2* h_TM = (TH2*)file->Get(name.c_str());
		h_TM->SetName(Form("%s_mc",name.c_str()));

		double lo_range = 0.5, hi_range = 1.5;


		//correct cone method in data:
		h_cone_data->Multiply(h_TM);
		h_cone_data->Divide(h_cone_MC);

		//ratios
		TH2* h_ratio_cone_data_mb_data = (TH2*)h_cone_data->Clone(Form("ratio_cone_data_mb_data_c%i",i_cent));
		h_ratio_cone_data_mb_data->Divide(h_MB_data);

		h_ratio_cone_data_mb_data->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]");
		h_ratio_cone_data_mb_data->GetYaxis()->SetTitle("p_{T}^{Jet} [GeV]");
		h_ratio_cone_data_mb_data->GetZaxis()->SetTitle("Cone Method / MC Method");

//		h_ratio_cone_data_mb_data->GetXaxis()->SetMoreLogLabels(kTRUE);
		h_ratio_cone_data_mb_data->GetXaxis()->SetRangeUser(1,10);
		h_ratio_cone_data_mb_data->GetYaxis()->SetRangeUser(128,350);
		h_ratio_cone_data_mb_data->SetMarkerSize(1.5);

		if (i_cent == 0){ lo_range = 0.95 ; hi_range = 1.01;}
		if (i_cent == 1){ lo_range = 0.96 ; hi_range = 1.03;}
		if (i_cent == 2){ lo_range = 0.975 ; hi_range = 1.04;}
		if (i_cent == 3){ lo_range = 0.99 ; hi_range = 1.15;}
		if (i_cent == 4){ lo_range = 0.95 ; hi_range = 1.1;}
		if (i_cent == 5){ lo_range = 1. ; hi_range = 1.39;}
//		h_ratio_cone_data_mb_data->GetZaxis()->SetRangeUser(lo_range,hi_range);

		c1->cd(i_cent+1);
		h_ratio_cone_data_mb_data->DrawCopy("colz text");
		gPad->SetLogx();
		gPad->SetRightMargin(0.18);
		ltx->SetTextAlign(32);
		ltx->DrawLatexNDC(0.90,0.98,num_to_cent(31,i_cent).c_str());


		//ratios
		TH2* h_ratio_cone_mc_tm = (TH2*)h_cone_MC->Clone(Form("ratio_cone_mc_tm_c%i",i_cent));
		h_ratio_cone_mc_tm->Divide(h_TM);

		h_ratio_cone_mc_tm->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]");
		h_ratio_cone_mc_tm->GetYaxis()->SetTitle("p_{T}^{Jet} [GeV]");
		h_ratio_cone_mc_tm->GetZaxis()->SetTitle("UE_{Cone}^{MC} / UE_{TM}^{MC}");

//		h_ratio_cone_mc_tm->GetXaxis()->SetMoreLogLabels(kTRUE);
		h_ratio_cone_mc_tm->GetXaxis()->SetRangeUser(1,10);
		h_ratio_cone_mc_tm->GetYaxis()->SetRangeUser(128,350);
		h_ratio_cone_mc_tm->SetMarkerSize(1.5);



		lo_range = 0.2 ; hi_range = 1.2;
		h_ratio_cone_mc_tm->GetZaxis()->SetRangeUser(lo_range,hi_range);

		c2->cd(i_cent+1);
		h_ratio_cone_mc_tm->DrawCopy("col text");
		gPad->SetLogx();
//		gPad->SetRightMargin(0.18);
		ltx->SetTextAlign(32);
		ltx->DrawLatexNDC(0.90,0.98,num_to_cent(31,i_cent).c_str());



		TH2* h_maddness = (TH2*)h_ratio_cone_data_mb_data->Clone("h_maddness");

		for (int i_trk = 3; i_trk < 8; i_trk++)
		{
			for (int i_jet = 8; i_jet < 13; i_jet++)
			{
				double ff_number = get_madness_factor(i_cent, i_trk, i_jet);
				double orig_number = h_maddness->GetBinContent(i_trk, i_jet);
				h_maddness->SetBinContent(i_trk, i_jet, orig_number/ff_number);
			}
		}

		
//		c3->cd(i_cent+1);
//		h_maddness->DrawCopy("colz text");
//		gPad->SetLogx();
//		gPad->SetRightMargin(0.18);
//		ltx->SetTextAlign(32);
//		ltx->DrawLatexNDC(0.90,0.98,num_to_cent(31,i_cent).c_str());

		c3->cd(i_cent+1);
		for (int i_jet = 8; i_jet < 12; i_jet++)
		{
			TH1* h_maddness_proj = (TH1*)h_maddness->ProjectionX(Form("jet%i",i_jet), i_jet,i_jet);
			SetHStyle_smallify(h_maddness_proj,i_jet-8,1);
			h_maddness_proj->GetYaxis()->SetRangeUser(0.9,1.1);
			if (i_jet == 8) h_maddness_proj->DrawCopy("");
			else h_maddness_proj->DrawCopy("same");
			if (i_cent == 0)
			{
				legend_x_jet->AddEntry(h_maddness_proj,Form("%1.0f < p_{T}^{Jet} < %1.0f",h_maddness->GetYaxis()->GetBinLowEdge(i_jet),h_maddness->GetYaxis()->GetBinUpEdge(i_jet) ));
			}
		}
		legend_x_jet->Draw();
		gPad->SetLogx();
//		gPad->SetRightMargin(0.18);
		ltx->SetTextAlign(32);
		ltx->DrawLatexNDC(0.90,0.98,num_to_cent(31,i_cent).c_str());



	}

	c1->Print(Form("UE_cone_mb_comparison_c%i.pdf",config));
	c2->Print(Form("UE_cone_tm_comparison_c%i.pdf",config));
	c3->Print(Form("UE_cone_madness_c%i.pdf",config));




	{

		TFile *f_cone04 = new TFile("UE_c24.root");
		TFile *f_cone08 = new TFile("UE_c25.root");
		TCanvas *c4 = new TCanvas("c4","c4",1000,1200);

		c4->Divide(1,2);
		c4->cd(1)->Divide(3,2);
		c4->cd(2)->Divide(3,2);

		for (int i_cent = 0; i_cent < 6; i_cent++)
		{

			name = Form("Cone_data_Cent%i",i_cent); cout << "*******make sure this is uncorrected UE******" << endl;
			TH2* h_cone_data_08 = (TH2*)f_cone08->Get(name.c_str());
			h_cone_data_08->SetName(Form("%s_data_08",name.c_str()));

			TH2* h_cone_data_04 = (TH2*)f_cone04->Get(name.c_str());
			h_cone_data_04->SetName(Form("%s_data_04",name.c_str()));


			name = Form("Cone_Cent%i",i_cent);
			TH2* h_cone_MC_08 = (TH2*)f_cone08->Get(name.c_str());
			h_cone_MC_08->SetName(Form("%s_mc_08",name.c_str()));

			TH2* h_cone_MC_04 = (TH2*)f_cone04->Get(name.c_str());
			h_cone_MC_04->SetName(Form("%s_mc_04",name.c_str()));

			TH2* h_TM_08 = (TH2*)f_cone08->Get(Form("TM_Cent%i",i_cent));
			h_cone_data_08->Multiply(h_TM_08);
			h_cone_data_08->Divide(h_cone_MC_08);

			TH2* h_TM_04 = (TH2*)f_cone04->Get(Form("TM_Cent%i",i_cent));
			h_cone_data_04->Multiply(h_TM_04);
			h_cone_data_04->Divide(h_cone_MC_04);




			TH2* h_ratio_data_08_04 = (TH2*)h_cone_data_08->Clone(Form("ratio_cone_data_08_04_c%i",i_cent));
			h_ratio_data_08_04->Divide(h_cone_data_04);

			TH2* h_ratio_mc_08_04 = (TH2*)h_cone_MC_08->Clone(Form("ratio_cone_mc_08_04_c%i",i_cent));
			h_ratio_mc_08_04->Divide(h_cone_MC_04);


			double lo_range = 0, hi_range = 2.;

			h_ratio_data_08_04->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]");
			h_ratio_data_08_04->GetYaxis()->SetTitle("p_{T}^{Jet} [GeV]");
			h_ratio_data_08_04->GetXaxis()->SetMoreLogLabels(kTRUE);
			h_ratio_data_08_04->GetXaxis()->SetRangeUser(1,10);
			h_ratio_data_08_04->GetYaxis()->SetRangeUser(128,350);
			h_ratio_data_08_04->GetZaxis()->SetRangeUser(0.9,1.1);
			h_ratio_data_08_04->SetMarkerSize(1.5);

			c4->cd(1)->cd(i_cent+1);
			h_ratio_data_08_04->DrawCopy("colz text");
			gPad->SetLogx();
			gPad->SetRightMargin(0.18);
			gPad->SetTopMargin(0.07);
			ltx->SetTextAlign(32);
			ltx->DrawLatexNDC(0.90,0.98,num_to_cent(31,i_cent).c_str());
			ltx->SetTextAlign(11);
			ltx->DrawLatexNDC(0.19,0.96,"UE_{ConeR = 0.8} / UE_{ConeR = 0.4} (Data)");

			///////


			h_ratio_mc_08_04->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]");
			h_ratio_mc_08_04->GetYaxis()->SetTitle("p_{T}^{Jet} [GeV]");
			h_ratio_mc_08_04->GetXaxis()->SetMoreLogLabels(kTRUE);
			h_ratio_mc_08_04->GetXaxis()->SetRangeUser(1,10);
			h_ratio_mc_08_04->GetYaxis()->SetRangeUser(128,350);
			h_ratio_mc_08_04->GetZaxis()->SetRangeUser(0.9,1.1);
			h_ratio_mc_08_04->SetMarkerSize(1.5);

			c4->cd(2)->cd(i_cent+1);
			h_ratio_mc_08_04->DrawCopy("colz text");
			gPad->SetLogx();
			gPad->SetRightMargin(0.18);
			gPad->SetTopMargin(0.07);
			ltx->SetTextAlign(32);
			ltx->DrawLatexNDC(0.90,0.96,num_to_cent(31,i_cent).c_str());
			ltx->SetTextAlign(11);
			ltx->DrawLatexNDC(0.19,0.96,"UE_{ConeR = 0.8} / UE_{ConeR = 0.4} (MC)");

		}


		c4->Print("tmp.pdf");
	}

}
