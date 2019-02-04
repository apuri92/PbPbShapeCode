void cone_stats(int config)
{
	SetAtlasStyle();
//	TFile *f_mc = new TFile("hist-local_mc_1.root");
//	TFile *f_data = new TFile("hist-local_mc_2.root");
	TFile *f_mc = new TFile(Form("output_dev/raw_results/c%i/FF_MC_out_histo_PbPb_5p02_r001.root",config));
	TFile *f_data = new TFile(Form("output_dev/raw_results/c%i/FF_data_out_histo_PbPb_5p02_r001.root",config));


	TH1* h_mc = (TH1*)f_mc->Get("h_tmp_cone_stats");
	h_mc->SetName("MC_ConeStats");

	TH1* h_data = (TH1*)f_data->Get("h_tmp_cone_stats");
	h_data->SetName("Data_ConeStats");

	TLatex *ltx = new TLatex();
	ltx->SetTextFont(43);
	ltx->SetTextSize(15);
	ltx->SetTextAlign(12);
	ltx->SetTextAngle(90);

	TLegend *legend = new TLegend(0.22, 0.80, 0.30, 0.9, "","brNDC");
	legend->SetTextFont(43);
	legend->SetBorderSize(0);
	legend->SetTextSize(28);


	TCanvas *c = new TCanvas("c","c",600,700);

	SetHStyle(h_mc,0);
	SetHStyle(h_data,1);

	h_mc->Scale(1./h_mc->Integral());
	h_data->Scale(1./h_data->Integral());

	c->cd();
//	h_mc->GetXaxis()->SetRangeUser(0.2,10);
	if (config == 24) h_mc->GetXaxis()->SetRangeUser(25,45);
	if (config == 25) h_mc->GetXaxis()->SetRangeUser(0,10);

	h_mc->GetYaxis()->SetRangeUser(0,0.3);

	h_mc->GetXaxis()->SetNdivisions(020,kTRUE);
	h_mc->GetXaxis()->SetTitle("Cones used");
	h_mc->GetYaxis()->SetTitle("% Events");
	legend->AddEntry(h_mc,"MC");
	legend->AddEntry(h_data,"Data");
	gPad->SetGridx();

	h_mc->Draw("hist");
	h_data->Draw("same hist");
	legend->Draw();

	string name;

	c->Print(Form("cone_stats_c%i.pdf",config));

}













