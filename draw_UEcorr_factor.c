#include "output_dev/functions/global_variables.h"

//DUPLICATE OF output_dev/unfold/UE_FACTORS.C
void draw_UEcorr_factor()
{
	TFile *file = new TFile("UE_nominal.root");
	TAxis* dR_binning = (TAxis*)file->Get("dR_binning");
	TAxis* jetpT_binning = (TAxis*)file->Get("jetpT_binning");
	TAxis* trkpT_binning = (TAxis*)file->Get("trkpT_binning");
	SetAtlasStyle();
	gStyle->SetPalette(kRainBow);
	int N_dR = dR_binning->GetNbins();
	int N_jetpt = jetpT_binning->GetNbins();
	int N_trkpt = trkpT_binning->GetNbins();

	TCanvas *c = new TCanvas("c","c",800,400);
	c->Divide(3,2);
	TLatex *ltx = new TLatex();
	ltx->SetTextFont(43);
	ltx->SetTextSize(10);
	ltx->SetTextAlign(12);

	TH2* h_correction_inJet;
	for (int i_dR = 0; i_dR < N_dR; i_dR++)
	{

		string dr_label = Form("%1.2f < r < %1.2f", dR_binning->GetBinLowEdge(i_dR+1), dR_binning->GetBinUpEdge(i_dR+1));

		for (int i_cent = 0; i_cent < 6; i_cent++)
		{
			TH2* h_corr = (TH2*)file->Get(Form("h_cone_correction_dR%i_cent%i",i_dR, i_cent));

			c->cd(i_cent+1);
			h_corr->GetXaxis()->SetRangeUser(1,10);
			h_corr->GetYaxis()->SetRangeUser(90,500);
			h_corr->GetYaxis()->SetMoreLogLabels(kTRUE);
			h_corr->GetYaxis()->SetNoExponent(kTRUE);
			h_corr->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]");
			h_corr->GetYaxis()->SetTitle("p_{T}^{Jet} [GeV]");
			h_corr->SetMarkerSize(2);

			//		h_corr->GetZaxis()->SetRangeUser(0,1.2);

			h_corr->DrawCopy("colz text");
			gPad->SetLogx();
			gPad->SetLogy();
			ltx->SetTextAlign(12);
			ltx->DrawLatexNDC(0.19,0.985,dr_label.c_str());
			ltx->SetTextAlign(32);
			ltx->DrawLatexNDC(0.99,0.985,num_to_cent(31, i_cent).c_str());


		}
//		c->cd(1);

		string name = "";
		if (i_dR == 0) name = "(";
		else if (i_dR == N_dR - 1) name = ")";

		c->Print(Form("misc_plots/CorrFactorComparison.pdf%s",name.c_str()),Form("Title: r%i, %s",i_dR, dr_label.c_str()));

	}

}
