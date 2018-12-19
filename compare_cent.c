#include "output_dev/functions/global_variables.h"

void compare_cent()
{
	SetAtlasStyle();
	gErrorIgnoreLevel = 3001;

	TFile *f_MB = new TFile(Form("output_dev/raw_results/c40/FF_MC_out_histo_PbPb_5p02_r001.root"));
	TFile *f_MC = new TFile(Form("output_dev/raw_results/c41/FF_MC_out_histo_PbPb_5p02_r001.root"));

	TH1 *h_cent_mc = (TH1*)f_MC->Get("Centrality");
	h_cent_mc->SetName("h_cent_mc");

	TH1 *h_cent_mbov = (TH1*)f_MB->Get("Centrality");
	h_cent_mbov->SetName("h_cent_mbov");

	TH1 *h_cent_ratio = (TH1*)h_cent_mbov->Clone("h_fcal_ratio");
	h_cent_ratio->Divide(h_cent_mc);


	TH1 *h_fcal_mc = (TH1*)f_MC->Get("h_FCal_Et");
	h_fcal_mc->SetName("h_fcal_mc");

	TH1 *h_fcal_mbov = (TH1*)f_MB->Get("h_FCal_Et");
	h_fcal_mbov->SetName("h_fcal_mbov");

	TH1 *h_fcal_ratio = (TH1*)h_fcal_mbov->Clone("h_fcal_ratio");
	h_fcal_ratio->Divide(h_fcal_mc);


	SetHStyle(h_fcal_mc, 0);
	SetHStyle(h_fcal_mbov, 1);
	SetHStyle(h_fcal_ratio, 2);
	SetHStyle(h_cent_mc, 0);
	SetHStyle(h_cent_mbov, 1);
	SetHStyle(h_cent_ratio, 2);



	TLegend *legend = new TLegend(0.19, 0.25, 0.80, 0.42, "","brNDC");
	legend->SetTextFont(43);
	legend->SetBorderSize(0);
	legend->SetTextSize(24);

	legend->AddEntry(h_fcal_mc,"Using MC FCal E_{T}","lpf");
	legend->AddEntry(h_fcal_mbov,"Using MinBias FCal E_{T}","lpf");
	legend->AddEntry(h_fcal_ratio,"MinBias/MC","lpf");

	TLatex *ltx = new TLatex();
	ltx->SetTextFont(43);
	ltx->SetTextSize(18);
	ltx->SetTextAlign(22);



	TCanvas *c1 = new TCanvas("c","c",1200,1200);
	TLine * line = new TLine();
	c1->Divide(2,2);

	c1->cd(1);
	h_cent_mc->GetXaxis()->SetRangeUser(0,6);
	h_cent_mc->GetYaxis()->SetRangeUser(1E-9, 1E-6);
	h_cent_mc->GetXaxis()->SetTitle("CentBin");
	h_cent_mbov->GetXaxis()->SetTitle("CentBin");
	h_cent_mc->Draw("p");
	h_cent_mbov->Draw("same p");
	gPad->SetLogy();
	legend->Draw();

	ltx->SetTextAngle(0);
	ltx->DrawLatex(5.5, h_cent_mc->GetBinContent(6)*0.2E1, "60-80%");
	ltx->DrawLatex(4.5, h_cent_mc->GetBinContent(5)*0.2E1, "40-60%");
	ltx->DrawLatex(3.5, h_cent_mc->GetBinContent(4)*0.2E1, "30-40%");
	ltx->DrawLatex(2.5, h_cent_mc->GetBinContent(3)*0.2E1, "20-30%");
	ltx->DrawLatex(1.5, h_cent_mc->GetBinContent(2)*0.2E1, "10-20%");
	ltx->DrawLatex(0.5, h_cent_mc->GetBinContent(1)*0.2E1, " 0-10%");



	c1->cd(2);
	h_cent_ratio->GetXaxis()->SetRangeUser(0,6);
	h_cent_ratio->GetYaxis()->SetRangeUser(0.95,1.05);
	h_cent_ratio->GetXaxis()->SetTitle("CentBin");
	h_cent_ratio->Draw("hist text");
	line->DrawLine(0,1,6,1);

	ltx->SetTextAngle(0);
	ltx->DrawLatex(5.5, 0.985, "60-80%");
	ltx->DrawLatex(4.5, 0.98, "40-60%");
	ltx->DrawLatex(3.5, 0.985, "30-40%");
	ltx->DrawLatex(2.5, 0.98, "20-30%");
	ltx->DrawLatex(1.5, 0.985, "10-20%");
	ltx->DrawLatex(0.5, 0.98, " 0-10%");


	c1->cd(3);
//	h_fcal_mc->GetXaxis()->SetRangeUser(0,6);
	h_fcal_mc->GetYaxis()->SetRangeUser(1E-10,0.5E-8);
	h_fcal_mc->GetXaxis()->SetRangeUser(1E-2,6);
	gPad->SetLogx();
	gPad->SetLogy();
	h_fcal_mc->Draw("p");
	h_fcal_mbov->Draw("same p");

	ltx->SetTextAngle(90);
	ltx->DrawLatex((0.289595+0.063719)/2, 1.8E-9, "60-80%");
	ltx->DrawLatex((0.87541+0.289595)/2, 2.2E-9, "40-60%");
	ltx->DrawLatex((1.36875+0.87541)/2, 1.8E-9, "30-40%");
	ltx->DrawLatex((2.04651+1.36875)/2, 2.2E-9, "20-30%");
	ltx->DrawLatex((2.98931+2.04651)/2, 1.8E-9, "10-20%");
	ltx->DrawLatex((6.0+2.98931)/2, 2.2E-9, " 0-10%");

	line->DrawLine(0.063719,2E-10,0.063719,2E-9);
	line->DrawLine(0.289595,2E-10,0.289595,2E-9);
	line->DrawLine(0.87541,2E-10,0.87541,2E-9);
	line->DrawLine(1.36875,2E-10,1.36875,2E-9);
	line->DrawLine(2.04651,2E-10,2.04651,2E-9);
	line->DrawLine(2.98931,2E-10,2.98931,2E-9);
	line->DrawLine(6.00,2E-10,6.00,1E-9);


	c1->cd(4);
	h_fcal_ratio->GetXaxis()->SetRangeUser(0,6);
	h_fcal_ratio->GetYaxis()->SetRangeUser(0,2);
	h_fcal_ratio->Draw("hist");
	
	line->DrawLine(0,1,6,1);
	ltx->SetTextAngle(90);
	ltx->DrawLatex((0.289595+0.063719)/2, 1.4, "60-80%");
	ltx->DrawLatex((0.87541+0.289595)/2, 1.6, "40-60%");
	ltx->DrawLatex((1.36875+0.87541)/2, 1.4, "30-40%");
	ltx->DrawLatex((2.04651+1.36875)/2, 1.6, "20-30%");
	ltx->DrawLatex((2.98931+2.04651)/2, 1.4, "10-20%");
	ltx->DrawLatex((6.0+2.98931)/2, 1.6, " 0-10%");

	line->DrawLine(0.063719,0.5,0.063719,1.5);
	line->DrawLine(0.289595,0.5,0.289595,1.5);
	line->DrawLine(0.87541,0.5,0.87541,1.5);
	line->DrawLine(1.36875,0.5,1.36875,1.5);
	line->DrawLine(2.04651,0.5,2.04651,1.5);
	line->DrawLine(2.98931,0.5,2.98931,1.5);
	line->DrawLine(6.00,0.5,6.00,1.5);


//	line->DrawLine(2.98931,0.5,2.98931,1.5);
//	line->DrawLine(2.04651,0.5,2.04651,1.5);
//	line->DrawLine(1.36875,0.5,1.36875,1.5);
//	line->DrawLine(0.87541,0.5,0.87541,1.5);
//	line->DrawLine(0.525092,0.5,0.525092,1.5);

//	if ( 2.98931 	<=centrality && centrality< 6.00  ) return 0;		// 0-10%
//	if ( 2.04651	<=centrality && centrality< 2.98931  ) return 1;	// 10-20%
//	if ( 1.36875	<=centrality && centrality< 2.04651  ) return 2;	// 20-30%
//	if ( 0.87541	<=centrality && centrality< 1.36875  ) return 3;	// 30-40%
//	if ( 0.525092	<=centrality && centrality< 0.87541  ) return 4;	// 40-50%
//	if ( 0.289595	<=centrality && centrality< 0.525092 ) return 4;	// 50-60%
//	if ( 0.063719	<=centrality && centrality< 0.289595 ) return 5;	// 60-80%



	c1->Print("cent.pdf");
}


//
//	TFile *f1 = new TFile(Form("output_dev/raw_results/c40/FF_MC_out_histo_PbPb_5p02_r001.root"));
//
//	TH2* h_fcal_change = (TH2*)f1->Get("h_fcal_change");
//	TH1* h_cent_mc = (TH1*)h_fcal_change->ProjectionX("h_cent_mc");
//	TH1* h_cent_mbov = (TH1*)h_fcal_change->ProjectionY("h_cent_mbov");
//	TH1 *h_cent_ratio = (TH1*)h_cent_mbov->Clone("h_fcal_ratio");
//	h_cent_ratio->Divide(h_cent_mc);
//
//
//	TH1 *h_fcal_mc = (TH1*)f1->Get("h_fcal_mc");
//	TH1 *h_fcal_mbov = (TH1*)f1->Get("h_fcal_mbov");
//	TH1 *h_fcal_ratio = (TH1*)h_fcal_mbov->Clone("h_fcal_ratio");
//	h_fcal_ratio->Divide(h_fcal_mc);
//
