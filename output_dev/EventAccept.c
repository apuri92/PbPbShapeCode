void EventAccept()
{
	SetAtlasStyle();

	TFile *f_PbPb = new TFile("raw_results/FF_data_out_histo_PbPb_5p02_r001.root");
	TFile *f_pp = new TFile("raw_results/FF_data_out_histo_pp_5p02_r001.root");

	TCanvas *c_PbPb = new TCanvas("c_PbPb","c_PbPb", 800,600);
	TCanvas *c_pp = new TCanvas("c_pp","c_pp", 800,600);

	string name = "RejectionHisto";
	TH1* h_accept_PbPb = (TH1*)(TH1*)f_PbPb->Get(name.c_str())->Clone(Form("PbPb_%s", name.c_str()));
	TH1* h_accept_pp = (TH1*)(TH1*)f_pp->Get(name.c_str())->Clone(Form("pp_%s", name.c_str()));

	h_accept_PbPb->GetYaxis()->SetTitle("#Events");
	h_accept_pp->GetYaxis()->SetTitle("#Events");

	SetHStyle(h_accept_PbPb, 1);
	SetHStyle(h_accept_pp, 1);

	c_PbPb->cd();
	h_accept_PbPb->Draw("hist");
	h_accept_PbPb->GetXaxis()->LabelsOption("u");
	gPad->SetLogy();

	c_pp->cd();
	h_accept_pp->Draw("hist");
	h_accept_pp->GetXaxis()->LabelsOption("u");
	gPad->SetLogy();

	c_PbPb->Print("output_pdf/EventAccept_PbPb.pdf");
	c_pp->Print("output_pdf/EventAccept_pp.pdf");



}
