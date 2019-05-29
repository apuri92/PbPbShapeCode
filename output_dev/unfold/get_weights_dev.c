#include "../functions/global_variables.h"
#include "draw_functions.c"

void get_ChPS_weights(string dataset_type, TFile *f_weights)
{
	gErrorIgnoreLevel = 3001;

	TFile *f_mc = new TFile(Form("output_pdf_nominal/root/final_ChPS_MC_%s.root", dataset_type.c_str()));
	TFile *f_data = new TFile(Form("output_pdf_nominal/root/final_ChPS_data_%s.root", dataset_type.c_str()));
//	TFile *f_weights = new TFile(Form("output_pdf_nominal/root/shape_spectra_weights_%s.root", dataset_type.c_str()), "recreate");

	cout << "Using: " << endl;
	cout << f_mc->GetName() << endl;
	cout << f_data->GetName() << endl;

	TAxis* dR_binning = (TAxis*)f_mc->Get("dR_binning");
	TAxis* jetpT_binning = (TAxis*)f_mc->Get("jetpT_binning");
	TAxis* trkpT_binning = (TAxis*)f_mc->Get("trkpT_binning");

	int N_dR = dR_binning->GetNbins();
	int N_jetpt = jetpT_binning->GetNbins();
	int N_trkpt = trkpT_binning->GetNbins();

	double trk_pt_lo = 1.;
	double trk_pt_hi = 150.;

	double ratio_lo = 0;
	double ratio_hi = 2;

	int jet_pt_start = 7;
	int jet_pt_end = 11;

	int trk_pt_start = 2;
	int trk_pt_end = 5;

	bool doSmall;
	if (dataset_type == "pp") doSmall = false;
	if (dataset_type == "PbPb") doSmall = true;

	TCanvas *c = new TCanvas("c","c",800,600);
	string name;

	TLine *line = new TLine();
	line->SetLineColor(kBlack);

	TLatex *ltx = new TLatex();
	ltx->SetTextFont(43);
	ltx->SetTextSize(16);
	ltx->SetTextAlign(12);

	vector<vector<vector<TH1*>>> h_ChPS_raw_MC (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_raw_data (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_ratio (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));

	vector<vector<vector<TH1*>>> h_ChPS_res_RW_MC (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_res_RW_data (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_res_RW_ratio (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));


	vector<vector<TH1*>> h_jet_MC (n_cent_cuts, vector<TH1*> (N_JET_Y+1));
	vector<vector<TH1*>> h_jet_data (n_cent_cuts, vector<TH1*> (N_JET_Y+1));
	vector<vector<TH1*>> h_jet_ratio (n_cent_cuts, vector<TH1*> (N_JET_Y+1));

	c->Print(Form("tmpfit_%s.pdf(", dataset_type.c_str()),"Title: start");
	for (int i_cent = 0; i_cent < n_cent_cuts; i_cent++)
	{
		if (dataset_type == "PbPb" && i_cent == 6) continue;
		if (dataset_type == "pp" && i_cent < 6) continue;

		for (int i_dR = 0; i_dR < N_dR; i_dR++)
		{
			for (int i_jet = 0; i_jet < N_jetpt; i_jet++)
			{

				name = Form("h_ChPS_raw_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet);
				h_ChPS_raw_MC[i_dR][i_cent][i_jet] = (TH1*)((TH1*)f_mc->Get(name.c_str()))->Clone(Form("%s_MC", name.c_str()));
				h_ChPS_raw_data[i_dR][i_cent][i_jet] = (TH1*)((TH1*)f_data->Get(name.c_str()))->Clone(Form("%s_data", name.c_str()));

				h_ChPS_ratio[i_dR][i_cent][i_jet] = (TH1*)h_ChPS_raw_data[i_dR][i_cent][i_jet]->Clone(Form("%s_ratio", name.c_str()));
				h_ChPS_ratio[i_dR][i_cent][i_jet]->Divide(h_ChPS_raw_MC[i_dR][i_cent][i_jet]);

				name = Form("h_ChPS_raw_subtr_unf_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet);
				h_ChPS_res_RW_MC[i_dR][i_cent][i_jet] = (TH1*)((TH1*)f_mc->Get(name.c_str()))->Clone(Form("%s_MC", name.c_str()));
				h_ChPS_res_RW_data[i_dR][i_cent][i_jet] = (TH1*)((TH1*)f_data->Get(name.c_str()))->Clone(Form("%s_data", name.c_str()));

				h_ChPS_res_RW_ratio[i_dR][i_cent][i_jet] = (TH1*)h_ChPS_res_RW_data[i_dR][i_cent][i_jet]->Clone(Form("%s_ratio", name.c_str()));
				h_ChPS_res_RW_ratio[i_dR][i_cent][i_jet]->Divide(h_ChPS_res_RW_MC[i_dR][i_cent][i_jet]);


				f_weights->cd();
				name = Form("CHPS_weight_%s_dR%i_cent%i_jet%i", dataset_type.c_str(), i_dR, i_cent, i_jet);
				h_ChPS_ratio[i_dR][i_cent][i_jet]->SetName(name.c_str());
				h_ChPS_ratio[i_dR][i_cent][i_jet]->SetTitle(name.c_str());
				h_ChPS_ratio[i_dR][i_cent][i_jet]->Write(name.c_str());

				name = Form("CHPS_resRW_weight_%s_dR%i_cent%i_jet%i", dataset_type.c_str(), i_dR, i_cent, i_jet);
				h_ChPS_res_RW_ratio[i_dR][i_cent][i_jet]->SetName(name.c_str());
				h_ChPS_res_RW_ratio[i_dR][i_cent][i_jet]->SetTitle(name.c_str());
				h_ChPS_res_RW_ratio[i_dR][i_cent][i_jet]->Write(name.c_str());
			}
		}

		double norm = 1;
		for (int i_y = 0; i_y < N_JET_Y+1; i_y++)
		{
			name = Form("h_reco_jet_y%i_c%i", i_y, i_cent);
			h_jet_MC[i_cent][i_y] = (TH1*)(TH1*)f_mc->Get(name.c_str())->Clone(Form("%s_MC", name.c_str()));
			norm = h_jet_MC[i_cent][i_y]->Integral();
			h_jet_MC[i_cent][i_y]->Scale(1./norm);

			h_jet_data[i_cent][i_y] = (TH1*)(TH1*)f_data->Get(name.c_str())->Clone(Form("%s_data", name.c_str()));
			norm = h_jet_data[i_cent][i_y]->Integral();
			h_jet_data[i_cent][i_y]->Scale(1./norm);


			h_jet_ratio[i_cent][i_y] = (TH1*)h_jet_data[i_cent][i_y]->Clone(Form("%s_ratio", name.c_str()));
			h_jet_ratio[i_cent][i_y]->Divide(h_jet_MC[i_cent][i_y]);

			TF1* f_fit = new TF1("fit","[0] + [1]*log(x) + [2]*pow(log(x),2)", 100,400);
			h_jet_ratio[i_cent][i_y]->Fit(f_fit,"qR","");

			c->cd();
			h_jet_ratio[i_cent][i_y]->GetYaxis()->SetRangeUser(0,2);
			h_jet_ratio[i_cent][i_y]->Draw("hist text");
			gPad->SetLogx();
			c->Print(Form("tmpfit_%s.pdf", dataset_type.c_str()),Form("Title: %s", name.c_str()));

			f_weights->cd();
			name = Form("jet_weight_%s_y%i_c%i",dataset_type.c_str(), i_y, i_cent);
			h_jet_ratio[i_cent][i_y]->SetName(name.c_str());
			h_jet_ratio[i_cent][i_y]->SetTitle(name.c_str());
			h_jet_ratio[i_cent][i_y]->Write(name.c_str());
			f_fit->Write(Form("fit_%s",name.c_str()));
		}

	}

	c->Print(Form("tmpfit_%s.pdf)", dataset_type.c_str()),"Title: start");

}

void get_weights_dev()
{
	TFile *f_weights = new TFile(Form("output_pdf_nominal/root/shape_spectra_weights.root"), "recreate");

	get_ChPS_weights("PbPb", f_weights);
	get_ChPS_weights("pp", f_weights);

}
