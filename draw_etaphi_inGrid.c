void draw_etaphi_inGrid()
{
	gStyle->SetPalette(kBird);
	gStyle->SetOptStat(00);
	TFile *f0 = new TFile("output_dev/raw_results/c97/FF_MC_out_histo_PbPb_5p02_r001.root");
	TFile *f1 = new TFile("output_dev/raw_results/c98/FF_MC_out_histo_PbPb_5p02_r001.root");
//	TFile *f0 = new TFile("hist-local_mc.root");


	TH3 *h0 = (TH3*)f0->Get("h_tmp_rdEtadPhi");
	h0->SetName("h0");
	h0->GetXaxis()->SetTitle("r");
	h0->GetYaxis()->SetTitle("#delta#eta");
	h0->GetZaxis()->SetTitle("#delta#phi");

	TH3 *h1 = (TH3*)f1->Get("h_tmp_rdEtadPhi");
	h1->SetName("h1");
	h1->GetXaxis()->SetTitle("r");
	h1->GetYaxis()->SetTitle("#delta#eta");
	h1->GetZaxis()->SetTitle("#delta#phi");

	TH2* h_tmp0;
	TH2* h_tmp1;
	TH2* h_ratio;
	TCanvas *c = new TCanvas("c","c",900,300);

	int first_bin = h0->GetXaxis()->FindBin(0.0);
	int last_bin = h0->GetXaxis()->FindBin(0.6);
	for (int i = first_bin; i < last_bin; i++)
	{
		h0->GetXaxis()->SetRange(i,i);
		h1->GetXaxis()->SetRange(i,i);
		c->Clear();
		c->cd();
		c->Divide(3,1);
		string name = Form("r: %1.2f - %1.2f", h0->GetXaxis()->GetBinLowEdge(i), h0->GetXaxis()->GetBinUpEdge(i));

		h_tmp0 = (TH2*)h0->Project3D("yz");
		h_tmp1 = (TH2*)h1->Project3D("yz");

		h_tmp0->SetTitle(Form("%s_0",name.c_str()));
		h_tmp1->SetTitle(Form("%s_1",name.c_str()));

		h_tmp0->GetZaxis()->SetRangeUser(0.,1E-8);
		h_tmp1->GetZaxis()->SetRangeUser(0.,1E-8);

		h_ratio = (TH2*)h_tmp0->Clone(Form("ratio_%i",i));
		h_ratio->Divide(h_tmp1);

		c->cd(1);
		h_tmp0->DrawCopy("colz");

		c->cd(2);
		h_tmp1->DrawCopy("colz");

		c->cd(3);
		h_ratio->GetZaxis()->SetRangeUser(0.9,1.1);
		h_ratio->DrawCopy("lego");

//		h_tmp->DrawCopy("colz");
		if (i == first_bin) c->Print("proj_tmp.pdf(",Form("Title: %i",i));
//		else if (i== last_bin-1) c->Print("proj_tmp.pdf)",Form("Title: %i",i));
		else c->Print("proj_tmp.pdf",Form("Title: %i",i));
	}

	h0->GetXaxis()->SetRangeUser(0.,0.6);
	h1->GetXaxis()->SetRangeUser(0.,0.6);
	h_tmp0 = (TH2*)h0->Project3D("yz");
	h_tmp1 = (TH2*)h1->Project3D("yz");
	h_ratio = (TH2*)h_tmp0->Clone(Form("ratio_x"));
	h_ratio->Divide(h_tmp1);

	c->cd(1);
	h_tmp0->DrawCopy("colz");

	c->cd(2);
	h_tmp1->DrawCopy("colz");

	c->cd(3);
	h_ratio->GetZaxis()->SetRangeUser(0.9,1.1);
	h_ratio->DrawCopy("colz");
	gPad->SetRightMargin(0.2);
	c->Print("proj_tmp.root");

	c->Print("proj_tmp.pdf)",Form("Title: full"));
}
