void EventAccept()
{
	SetAtlasStyle();

	TFile *f_PbPb = new TFile("../raw_results/nominal/FF_data_out_histo_PbPb_5p02_r001.root");
	TFile *f_pp = new TFile("../raw_results/nominal/FF_data_out_histo_pp_5p02_r001.root");


	TCanvas *c_PbPb = new TCanvas("c_PbPb","c_PbPb", 800,600);
	TCanvas *c_pp = new TCanvas("c_pp","c_pp", 800,600);

	string name = "RejectionHisto";
	TH1* h_accept_PbPb = (TH1*)(TH1*)f_PbPb->Get(name.c_str())->Clone(Form("PbPb_%s", name.c_str()));
	TH1* h_accept_pp = (TH1*)(TH1*)f_pp->Get(name.c_str())->Clone(Form("pp_%s", name.c_str()));

	h_accept_PbPb->GetYaxis()->SetTitle("#Events");
	h_accept_pp->GetYaxis()->SetTitle("#Events");

	SetHStyle(h_accept_PbPb, 1);
	SetHStyle(h_accept_pp, 1);

	TLatex *ltx = new TLatex();
	ltx->SetTextFont(43);
	ltx->SetTextSize(18);
	ltx->SetTextAlign(12);
	ltx->SetTextAngle(90);


	c_PbPb->cd();
	h_accept_PbPb->Draw("hist");
	h_accept_PbPb->GetXaxis()->SetLabelFont(43);
	h_accept_PbPb->GetXaxis()->SetLabelSize(12);

	gPad->SetLogy();

	for (int i = 0; i < h_accept_PbPb->GetXaxis()->GetNbins() ; i++)
	{
		if (i == 0) name = "Total number of events";
		else if (i == 7) name = "Accepted events";
		else if (i == 8) name = "Passed trigger selection Events";
		else name = Form("Rejected by %s", h_accept_PbPb->GetXaxis()->GetBinLabel(i+1));

		ltx->DrawLatex(i+0.5, 1E4, name.c_str());
	}





	c_pp->cd();
	h_accept_pp->Draw("hist");
	h_accept_pp->GetXaxis()->SetLabelFont(43);
	h_accept_pp->GetXaxis()->SetLabelSize(12);

	gPad->SetLogy();

	for (int i = 0; i < h_accept_PbPb->GetXaxis()->GetNbins() ; i++)
	{
		if (i == 0) name = "Total number of events";
		else if (i == 7) name = "Accepted events";
		else if (i == 8) name = "Passed trigger selection Events";
		else name = Form("Rejected by %s", h_accept_PbPb->GetXaxis()->GetBinLabel(i+1));

		ltx->DrawLatex(i+0.5, 2.5E3, name.c_str());
	}



	c_PbPb->Print("output_pdf_nominal/EventAccept_PbPb.pdf");
	c_pp->Print("output_pdf_nominal/EventAccept_pp.pdf");



}
