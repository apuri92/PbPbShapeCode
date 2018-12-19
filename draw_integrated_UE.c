#include "output_dev/functions/global_variables.h"

void UE_integrated()
{
	SetAtlasStyle();
	TFile *file = new TFile("UE_c21.root");

	//integrated over jet cone
	TCanvas *c = new TCanvas("c","c",700,900);
	c->Divide(2,3);
	string name;
	TLatex *ltx = new TLatex();
	ltx->SetTextFont(43);
	ltx->SetTextSize(12);
	ltx->SetTextAlign(12);

	for (int i_cent = 0; i_cent < 6; i_cent++)
	{
		name = Form("ChPS_cone_UE_dR0_cent%i_0_injet",i_cent);
		TH2* h_cone = (TH2*)file->Get(name.c_str());

		name = Form("ChPS_TM_UE_dR0_cent%i_0_injet",i_cent);
		TH2* h_mb = (TH2*)file->Get(name.c_str());

		TH2* h_ratio = (TH2*)h_cone->Clone(Form("ratio_c%i",i_cent));
		h_ratio->Divide(h_mb);

		h_ratio->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]");
		h_ratio->GetYaxis()->SetTitle("p_{T}^{Jet} [GeV]");
		h_ratio->GetZaxis()->SetTitle("UE_{Cone} / UE_{TM}");

		h_ratio->GetXaxis()->SetRangeUser(1,10);
		h_ratio->GetYaxis()->SetRangeUser(128,350);
//		h_ratio->GetZaxis()->SetRangeUser(0,2);
		h_ratio->GetZaxis()->SetNdivisions(504);

		h_ratio->SetMarkerSize(1.5);
		c->cd(i_cent+1);
		h_ratio->Draw("colz text");
		gPad->SetLogx();
		gPad->SetRightMargin(0.18);
//		gPad->SetLeftMargin(0.2);
//		gPad->SetTopMargin(0.2);
//		gPad->SetBottomMargin(0.2);
		ltx->SetTextAlign(32);
		ltx->DrawLatexNDC(0.90,0.98,num_to_cent(31,i_cent).c_str());

	}

	c->Print("UE_integrated_comparison.pdf");
}
