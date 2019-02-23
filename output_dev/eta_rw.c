#include "functions/global_variables.h"

void eta_rw()
{
	SetAtlasStyle();
	gErrorIgnoreLevel = 3001;


	TFile *f_mc = new TFile("raw_results/c20/FF_MC_out_histo_PbPb_5p02_r001.root");
	TFile* f_data = new TFile("raw_results/c20/FF_data_out_histo_PbPb_5p02_r001.root");
	TFile *eta_factors = new TFile("eta_factors.root","recreate");

	TAxis* dR_binning = (TAxis*)((TH3*)f_mc->Get("h_dR_change_jetpt0_cent0"))->GetXaxis();
	TAxis* jetpT_binning = (TAxis*)((TH3*)f_mc->Get("ChPS_raw_0_dR0_cent0"))->GetYaxis();
	TAxis* trkpT_binning = (TAxis*)((TH3*)f_mc->Get("ChPS_raw_0_dR0_cent0"))->GetXaxis();

	int jet_pt_start = 7;
	int jet_pt_end = 11;

	int N_dR = dR_binning->GetNbins();
	int N_jetpt = jetpT_binning->GetNbins();
	int N_trkpt = trkpT_binning->GetNbins();

	eta_factors->cd();
	jetpT_binning->Write("jetpT_binning");
	string name;

	vector<vector<TH1*>> h_ratio = vector<vector<TH1*>> (6, vector<TH1*> (N_jetpt));


	for (int i_cent = 0; i_cent < 6; i_cent++)
	{
		name = Form("h_reco_jets_cent%i",i_cent);

		TH3* jet_mc = (TH3*)f_mc->Get(name.c_str());
		jet_mc->SetName(Form("%s_mc", name.c_str()));

		TH3* jet_data = (TH3*)f_data->Get(name.c_str());
		jet_data->SetName(Form("%s_data", name.c_str()));

		for (int i_jet = 0; i_jet < N_jetpt; i_jet++)
		{

			name = Form("mc_eta_jet%i",i_jet);
			jet_mc->GetXaxis()->SetRange(i_jet+1, i_jet+1);
			TH1* eta_mc = jet_mc->Project3D("y");

			name = Form("data_eta_jet%i",i_jet);
			jet_data->GetXaxis()->SetRange(i_jet+1, i_jet+1);
			TH1* eta_data = jet_data->Project3D("y");

			eta_mc->Scale(1./eta_mc->Integral());
			eta_data->Scale(1./eta_data->Integral());

			string jet_label = Form("%1.0f < p_{T}^{Jet} < %1.0f", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));

			name = Form("eta_factor_cent%i_jet%i", i_cent, i_jet);
			h_ratio[i_cent][i_jet] = (TH1*)eta_data->Clone(name.c_str());


			if (eta_data->GetEntries() == 0)
			{
				h_ratio[i_cent][i_jet]->Reset();
				for (int i = 0; i <h_ratio[i_cent][i_jet]->GetXaxis()->GetNbins() ; i++)
				{
					h_ratio[i_cent][i_jet]->SetBinContent(i+1,1);
				}
			}
			else h_ratio[i_cent][i_jet]->Divide(eta_mc);

//			cout << Form("%s, %s, %f, %f",name.c_str(), jet_label.c_str(), eta_data->GetEntries(), h_ratio[i_cent][i_jet]->GetEntries()) << endl;

			eta_factors->cd();
			h_ratio[i_cent][i_jet]->SetName(name.c_str());
			h_ratio[i_cent][i_jet]->SetTitle(name.c_str());
			h_ratio[i_cent][i_jet]->Write(name.c_str());
		}
	}

	TCanvas *c1 = new TCanvas();
	c1->Divide(3,2);

	TLegend *legend = new TLegend(0.19, 0.70, 0.40, 0.92, "","brNDC");
	legend->SetTextFont(43);
	legend->SetBorderSize(0);
	legend->SetTextSize(10);

	TLatex *ltx = new TLatex();
	ltx->SetTextFont(43);
	ltx->SetTextSize(12);
	ltx->SetTextAlign(12);

	for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
	{
		string jet_label = Form("%1.0f < p_{T}^{Jet} < %1.0f", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));

		c1->Clear();
		c1->Divide(3,2);

		for (int i_cent = 0; i_cent < 6; i_cent++)
		{
			c1->cd(i_cent+1);
			SetHStyle_smallify(h_ratio[i_cent][i_jet], i_cent, 1);

			h_ratio[i_cent][i_jet]->GetYaxis()->SetRangeUser(0.7,1.3);
			h_ratio[i_cent][i_jet]->GetXaxis()->SetRangeUser(-0.5,0.5);

			h_ratio[i_cent][i_jet]->GetXaxis()->SetTitle("#eta^{jet}");
			h_ratio[i_cent][i_jet]->GetYaxis()->SetTitle("Factor");

			if (i_cent == 0) h_ratio[i_cent][i_jet]->Draw();
			else h_ratio[i_cent][i_jet]->Draw("same");

			ltx->SetTextAlign(12);
			ltx->DrawLatexNDC(0.19,0.90,jet_label.c_str());
			ltx->SetTextAlign(32);
			ltx->DrawLatexNDC(0.94,0.90,num_to_cent(31,i_cent).c_str());

		}

		
		string pdf_label = "";
		if (i_jet == jet_pt_start) pdf_label = "(";
		else if (i_jet == jet_pt_end-1) pdf_label = ")";
		c1->Print(Form("eta_fac.pdf%s", pdf_label.c_str()), Form("Title: jet%i", i_jet));

		//		legend->Draw();
	}
//	c1->Print("eta_fac.pdf");



//	c1->cd();
//	SetHStyle(eta_mc, 1);
//	SetHStyle(eta_data, 0);
//
//	eta_data->Draw();
//	eta_mc->Draw("same");
//
//
//	string pdf_label = "";
//	if (i_cent == 0 && i_jet == jet_pt_start) pdf_label = "(";
//	else if (i_cent == 5 && i_jet == jet_pt_end-1) pdf_label = ")";
//
//	c1->Print(Form("eta_factors.pdf%s", pdf_label.c_str()), Form("Title: jet%i_cent%i", i_jet, i_cent));


}
