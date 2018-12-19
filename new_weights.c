void new_weights()
{
	TFile *f = new TFile("pPbFragmentation/data/Powheg.reweight.root");
	TFile *f_out = new TFile("new_Powheg.reweight.root","recreate");

	TH3* h_old = (TH3*)f->Get("h3_pT_y_phi_rw");
	TH1* h_x_old = (TH1*)h_old->Project3D("x");

	TH3* h_new = (TH3*)h_old->Clone("h_new");
	h_new->Reset();
	TH1* h_x_new = (TH1*)h_new->Project3D("x");

//	TF1 * func = new TF1("func","1E8 * exp(-(log(0.002*pow(x,5))))");
	TF1 * func = new TF1("func","1/(pow(x,5))");

	for (int i_x = 0; i_x < h_old->GetXaxis()->GetNbins(); i_x++)
	{
		double x_val = h_x_old->GetBinContent(i_x+1);
//		h_x_new->SetBinContent(i_x, x_val);

		for (int i_y = 0; i_y < h_old->GetYaxis()->GetNbins(); i_y++)
		{
			for (int i_z = 0; i_z < h_old->GetZaxis()->GetNbins(); i_z++)
			{
				if (h_old->GetXaxis()->GetBinCenter(i_x) < 20) continue;
				h_new->SetBinContent(i_x+1, i_y+1, i_z+1, func->Eval(h_x_new->GetXaxis()->GetBinCenter(i_x+1)));
			}
		}
	}

	TCanvas *c = new TCanvas();

	TH1* h1 = (TH1*)h_old->Project3D("x");
	TH1* h2 = (TH1*)h_new->Project3D("x");
	h1->Scale(1./h1->Integral());
	h2->Scale(1./h2->Integral());

	h1->Divide(h2);
	h1->Draw("hist");
//	h2->Draw("same p");
	gPad->SetLogx();
//	gPad->SetLogy();

	c->Print("pow.pdf");
	f_out->cd();
	h_old->Write("h_old");
	h_new->Write("h_new");



}
