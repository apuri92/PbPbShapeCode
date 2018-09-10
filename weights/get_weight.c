	#include "../output_dev/functions/global_variables.h"

void NormInCentrality(TH1* h_in, TH1* h_out, double *fcal_values)
{
	h_in->Scale(1./h_in->Integral());
	h_out->Reset();

	for (int i = 1 ; i <= 6; i++)
	{
		int bin_adjust = 0;
		if (i>1) bin_adjust = 1;

		int lo_bin = h_in->FindBin(fcal_values[i-1])+bin_adjust;
		int hi_bin = h_in->FindBin(fcal_values[i]);

		cout << lo_bin << " " << hi_bin << endl;
		double integral = h_in->Integral(lo_bin, hi_bin);

		for (int j = lo_bin ; j <= hi_bin ; j++)
		{
			double new_fcal = h_in->GetBinContent(j)/integral;
			double new_fcal_err = h_in->GetBinError(j)/integral;

			if (h_in->GetBinCenter(j) > 4.6)
			{
				new_fcal = 1.;
				new_fcal_err = 0.;
			}

			h_out->SetBinContent(j, new_fcal);
			h_out->SetBinError(j, new_fcal_err);
		}
	}
}

void get_weight()
{
	SetAtlasStyle();

	TFile *f_data = new TFile("Data.root");
	TFile *f_mc = new TFile("MC.root");
	TFile *f_mb = new TFile("MB.root");
	TFile *f_mbov = new TFile("MBov.root");

	TFile *output = new TFile("fcal_weights.root","recreate");

	string name = "h_FCal_Et";
	TH1* h_hp_fcal = (TH1*)f_data->Get(name.c_str());
	h_hp_fcal->SetName(Form("%s_hp", name.c_str()));

	TH1* h_mc_fcal = (TH1*)f_mc->Get(Form("%s_unw",name.c_str()));
	h_mc_fcal->SetName(Form("%s_mc", name.c_str()));

	TH1* h_mb_fcal = (TH1*)f_mb->Get(name.c_str());
	h_mb_fcal->SetName(Form("%s_mb", name.c_str()));

	TH1* h_mbov_fcal = (TH1*)f_mbov->Get(name.c_str());
	h_mbov_fcal->SetName(Form("%s_mb", name.c_str()));

	name = "h_FCal_Et_new";
	TH1* h_hp_fcal_new = (TH1*)h_hp_fcal->Clone(name.c_str());
	TH1* h_mc_fcal_new = (TH1*)h_mc_fcal->Clone(name.c_str());
	TH1* h_mb_fcal_new = (TH1*)h_mb_fcal->Clone(name.c_str());
	TH1* h_mbov_fcal_new = (TH1*)h_mbov_fcal->Clone(name.c_str());

	double fcal_values[7] = {0, 0.289595, 0.87541, 1.36875, 2.04651, 2.98931, 6.00};

	//normalize to unity in each centrality bin
	NormInCentrality(h_hp_fcal, h_hp_fcal_new, fcal_values);
	NormInCentrality(h_mc_fcal, h_mc_fcal_new, fcal_values);
	NormInCentrality(h_mb_fcal, h_mb_fcal_new, fcal_values);
	NormInCentrality(h_mbov_fcal, h_mbov_fcal_new, fcal_values);

	//get weighing factors
	TH1* h_hp_mc = (TH1*)h_hp_fcal_new->Clone("MC_to_HP");
	h_hp_mc->Divide(h_mc_fcal_new);

	TH1* h_hp_mb = (TH1*)h_hp_fcal_new->Clone("MB_to_HP");
	h_hp_mb->Divide(h_mb_fcal_new);

	TH1* h_mc_mb = (TH1*)h_mc_fcal_new->Clone("MB_to_MC");
	h_mc_mb->Divide(h_mb_fcal_new);

	TH1* h_mc_mbov = (TH1*)h_mc_fcal_new->Clone("MBov_to_MC");
	h_mc_mbov->Divide(h_mbov_fcal_new);

	TH1* h_hp_mbov = (TH1*)h_hp_fcal_new->Clone("MBov_to_HP");
	h_hp_mbov->Divide(h_mbov_fcal_new);


	SetHStyle_smallify(h_hp_fcal, 0, 0);
	SetHStyle_smallify(h_mc_fcal, 1, 0);
	SetHStyle_smallify(h_mb_fcal, 2, 0);

	SetHStyle_smallify(h_hp_fcal_new, 0, 0);
	SetHStyle_smallify(h_mc_fcal_new, 1, 0);
	SetHStyle_smallify(h_mb_fcal_new, 2, 0);

	SetHStyle_smallify(h_hp_mc, 3, 0);
	SetHStyle_smallify(h_hp_mb, 4, 0);
	SetHStyle_smallify(h_mc_mb, 5, 0);



	TLegend *legend1 = new TLegend();
	legend1->AddEntry(h_hp_fcal, "Data", "lp");
	legend1->AddEntry(h_mc_fcal, "MC", "lp");
	legend1->AddEntry(h_mb_fcal, "MB", "lp");


	TLegend *legend2 = new TLegend();
	legend2->AddEntry(h_hp_mc, "Data/MC", "lp");
	legend2->AddEntry(h_hp_mb, "Data/MB", "lp");
	legend2->AddEntry(h_mc_mb, "MC/MB", "lp");


	TCanvas *c1 = new TCanvas("c1","c1",1200,1600);
	c1->cd();
	c1->Divide(1,3);
	c1->cd(1);
	h_mb_fcal->Draw("");
	h_hp_fcal->Draw("same");
	h_mc_fcal->Draw("same");
	gPad->SetLogy();
	legend1->Draw();

	c1->cd(2);
	h_mc_fcal_new->Draw("");
	h_hp_fcal_new->Draw("same");
	h_mb_fcal_new->Draw("same");
	gPad->SetLogy();


	c1->cd(3);
	h_hp_mb->Draw("");
	h_hp_mc->Draw("same");
	h_mc_mb->Draw("same");
	legend2->Draw();

	c1->Print("weights.pdf");

	output->cd();

	name = "h_FCal_Et_new";
	h_hp_fcal_new->Write(Form("%s_data", name.c_str()));
	h_mc_fcal_new->Write(Form("%s_mc", name.c_str()));
	h_mb_fcal_new->Write(Form("%s_mb", name.c_str()));

	name = "h_FCal_Et_old";
	h_hp_fcal->Write(Form("%s_data", name.c_str()));
	h_mc_fcal->Write(Form("%s_mc", name.c_str()));
	h_mb_fcal->Write(Form("%s_mb", name.c_str()));

	h_hp_mb->Write("weight_MB_to_HP");
	h_hp_mc->Write("weight_MC_to_HP");
	h_mc_mb->Write("weight_MB_to_MC");
	h_mc_mbov->Write("weight_MBov_to_MC");
	h_hp_mbov->Write("weight_MBov_to_HP");



	output->Close();



}

