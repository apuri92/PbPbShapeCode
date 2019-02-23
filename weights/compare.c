#include "../output_dev/functions/global_variables.h"


void compare()
{
	SetAtlasStyle();
	TFile *new_f = new TFile("fcal_weights.root");
	TFile *old_f = new TFile("../pPbFragmentation/data/FCal_HP_v_MB_weights.root");

	TH1* new_mb_to_hp = (TH1*)new_f->Get("weight_MB_to_HP");
	TH1* new_mc_to_hp = (TH1*)new_f->Get("weight_MC_to_HP");
	TH1* new_mb_to_mc = (TH1*)new_f->Get("weight_MB_to_MC");
	TH1* new_mbov_to_mc = (TH1*)new_f->Get("weight_MBov_to_MC");
	TH1* new_mbov_to_hp = (TH1*)new_f->Get("weight_MBov_to_HP");

	TH1* old_mb_to_hp = (TH1*)old_f->Get("FCal_HP_v_MB_v2_weights");
	TH1* old_mc_to_hp = (TH1*)old_f->Get("FCal_HP_v_MCOV_weights");
	TH1* old_mb_to_mc = (TH1*)old_f->Get("FCal_MCOV_v_MB_weights");
	TH1* old_mbov_to_mc = (TH1*)old_f->Get("FCal_MCOV_v_MBOV_weights");
	TH1* old_mbov_to_hp = (TH1*)old_f->Get("FCal_HP_v_MBOV_weights");

//	new_mb_to_hp->Divide(old_mb_to_hp);
//	new_mc_to_hp->Divide(old_mc_to_hp);
//	new_mb_to_mc->Divide(old_mb_to_mc);
//	new_mbov_to_mc->Divide(old_mbov_to_mc);
//	new_mbov_to_hp->Divide(old_mbov_to_hp);


	SetHStyle_smallify(new_mb_to_hp,0, 1);
	SetHStyle_smallify(new_mc_to_hp,0, 1);
	SetHStyle_smallify(new_mb_to_mc,0, 1);
	SetHStyle_smallify(new_mbov_to_mc,0, 1);
	SetHStyle_smallify(new_mbov_to_hp,0, 1);

	SetHStyle_smallify(old_mb_to_hp, 1, 1);
	SetHStyle_smallify(old_mc_to_hp, 1, 1);
	SetHStyle_smallify(old_mb_to_mc, 1, 1);
	SetHStyle_smallify(old_mbov_to_mc, 1, 1);
	SetHStyle_smallify(old_mbov_to_hp, 1, 1);

	gStyle->SetMarkerSize(1.2);


	TCanvas *c1 = new TCanvas();
	TLatex *ltx = new TLatex();
	ltx->SetTextFont(43);
	ltx->SetTextSize(16);
	ltx->SetTextAlign(12);

	c1->Divide(1,5);

	c1->cd(1);

	new_mb_to_hp->Draw();
	old_mb_to_hp->Draw("same");
	ltx->DrawLatexNDC(0.5,0.9, "mb to hp");
	gPad->SetLogx(0);

	c1->cd(2);
	new_mc_to_hp->Draw();
	old_mc_to_hp->Draw("same");
	ltx->DrawLatexNDC(0.5,0.9, "mc to hp");
	gPad->SetLogx(0);

	c1->cd(3);
	new_mb_to_mc->Draw();
	old_mb_to_mc->Draw("same");
	ltx->DrawLatexNDC(0.5,0.9, "mb to mc");
	gPad->SetLogx(0);

	c1->cd(4);
	new_mbov_to_mc->Draw();
	old_mbov_to_mc->Draw("same");
	ltx->DrawLatexNDC(0.5,0.9, "mbov to mc");
	gPad->SetLogx(0);

	c1->cd(5);
	new_mbov_to_hp->Draw();
	old_mbov_to_hp->Draw("same");
	ltx->DrawLatexNDC(0.5,0.9, "mbov to hp");
	gPad->SetLogx(0);


	int bin = new_mb_to_hp->GetXaxis()->FindBin(1.2);
	cout << old_mb_to_hp->GetBinContent(bin) << endl << endl;

	cout << old_mb_to_mc->GetBinContent(bin) << endl;
	cout << old_mc_to_hp->GetBinContent(bin) << endl << endl;;

	cout << old_mc_to_hp->GetBinContent(bin) * old_mb_to_mc->GetBinContent(bin)<< endl;


	c1->Print("comparison.pdf");


}

