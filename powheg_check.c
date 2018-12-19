void powheg_check()
{
	TFile *file = new TFile("pPbFragmentation/data/Powheg.reweight.root");

	double pt_bin[5] = {126, 158, 200, 251, 316};

	TCanvas *c = new TCanvas();

	for (int i = 0; i < 4 ; i++)
	{
		TH3* h_3d = (TH3*)file->Get("h3_pT_y_phi_rw");
		h_3d->GetXaxis()->SetRangeUser(pt_bin[i],pt_bin[i+1]);
		TH1* h_1d = (TH1*)h_3d->Project3D("y");
		h_1d->Scale(1./h_1d->Integral());
		if (i==0) h_1d->DrawCopy();
		else h_1d->DrawCopy("same");	

	}


}
