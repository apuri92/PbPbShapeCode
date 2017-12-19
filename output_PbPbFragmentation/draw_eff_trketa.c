//#include "combine_eff_dev.c"
#include "extras/global_variables.h"
#include "TVirtualFitter.h"

static const int n_coarse_eta = 6;
double coarse_eta[n_coarse_eta] = {-2.5, -2, -1, 1, 2, 2.5};

void draw_eff_trketa(string cut)
{
	gStyle->SetOptTitle(0);
	SetAtlasStyle();

	string name;

	name = Form("mc_efficiency_trketa_%s.root",cut.c_str());
	TFile *input_file = new TFile(name.c_str());

	name = Form("mc_eff_fits_trketa_%s.root",cut.c_str());
	TFile *output_file = new TFile(name.c_str(),"recreate");

//	TAxis *axis_trk_eta_new = (TAxis*)input_file->Get("new_trk_eta_binning");

	TCanvas *canvas1 = new TCanvas("C1", "C1",900,600);
	TLine *line = new TLine();
	TLatex *ltx = new TLatex();
	ltx->SetNDC();
	ltx->SetTextFont(43);
	ltx->SetTextSize(15);
	ltx->SetTextAlign(11);
	TLegend *legend = new TLegend(0.21,0.20,0.70,0.40,NULL,"brNDC");
	legend->SetBorderSize(0);
	legend->SetNColumns(1);
	legend->SetTextFont(43);
	legend->SetTextSize(12);

	vector<vector<TH1*>> h_efficiency(n_cent_cuts, vector<TH1*> (n_coarse_eta));
	vector<vector<TGraphAsymmErrors*>> g_efficiency(n_cent_cuts, vector<TGraphAsymmErrors*> (n_coarse_eta));

	canvas1->Divide(3,2);

	for (int i_cent_cuts = 0; i_cent_cuts < n_cent_cuts; i_cent_cuts++)
	{
		canvas1->cd(i_cent_cuts+1);

		string centrality = num_to_cent(centrality_scheme,i_cent_cuts);
		for (int i_eta_cuts = 0; i_eta_cuts < n_coarse_eta-1; i_eta_cuts++)
		{
			double eta_lo = coarse_eta[i_eta_cuts];
			double eta_hi = coarse_eta[i_eta_cuts+1];

			name = Form("graph_eff_eta%i_cent%i", i_eta_cuts, i_cent_cuts);
			g_efficiency.at(i_cent_cuts).at(i_eta_cuts) = (TGraphAsymmErrors*)input_file->Get(name.c_str());

			name = Form("histo_eff_eta%i_cent%i", i_eta_cuts, i_cent_cuts);
			h_efficiency.at(i_cent_cuts).at(i_eta_cuts) = (TH1*)input_file->Get(name.c_str());


			g_efficiency.at(i_cent_cuts).at(i_eta_cuts)->SetMinimum(0);
			g_efficiency.at(i_cent_cuts).at(i_eta_cuts)->SetMaximum(1.2);
			g_efficiency.at(i_cent_cuts).at(i_eta_cuts)->GetXaxis()->SetLimits(1,300);
			g_efficiency.at(i_cent_cuts).at(i_eta_cuts)->GetYaxis()->SetTitle("Efficiency");
			g_efficiency.at(i_cent_cuts).at(i_eta_cuts)->GetXaxis()->SetTitle("#it{p}_{T}^{truth} [GeV]");
			g_efficiency.at(i_cent_cuts).at(i_eta_cuts)->GetXaxis()->SetLabelSize(0.04);
			g_efficiency.at(i_cent_cuts).at(i_eta_cuts)->GetXaxis()->SetTitleOffset(0.75);
			g_efficiency.at(i_cent_cuts).at(i_eta_cuts)->GetXaxis()->SetTitleSize(0.06);
			g_efficiency.at(i_cent_cuts).at(i_eta_cuts)->GetYaxis()->SetLabelSize(0.04);
			g_efficiency.at(i_cent_cuts).at(i_eta_cuts)->GetYaxis()->SetTitleOffset(1.15);
			g_efficiency.at(i_cent_cuts).at(i_eta_cuts)->GetYaxis()->SetTitleSize(0.065);
			g_efficiency.at(i_cent_cuts).at(i_eta_cuts)->SetTitle(Form("Efficiency: %s, %4.2f < #eta < %4.2f",centrality.c_str(), eta_lo, eta_hi));

			SetHStyle_graph(g_efficiency.at(i_cent_cuts).at(i_eta_cuts),i_eta_cuts);
			smallify(g_efficiency.at(i_cent_cuts).at(i_eta_cuts));

//			canvas1->cd(i_eta_cuts);
			if (i_cent_cuts == 0) legend->AddEntry(g_efficiency.at(i_cent_cuts).at(i_eta_cuts),Form("%4.2f < #eta < %4.2f",eta_lo, eta_hi),"lp");
			if (i_eta_cuts == 0) g_efficiency.at(i_cent_cuts).at(i_eta_cuts)->Draw("a p");
			else g_efficiency.at(i_cent_cuts).at(i_eta_cuts)->Draw("same p");
			gPad->SetLogx();

			output_file->cd();
			name = Form("graph_eff_eta%i_cent%i",i_eta_cuts,i_cent_cuts);
			g_efficiency.at(i_cent_cuts).at(i_eta_cuts)->SetName(name.c_str());
			g_efficiency.at(i_cent_cuts).at(i_eta_cuts)->Write(name.c_str());

			name = Form("hist_eff_eta%i_cent%i",i_eta_cuts,i_cent_cuts);
			h_efficiency.at(i_cent_cuts).at(i_eta_cuts)->SetName(name.c_str());
			h_efficiency.at(i_cent_cuts).at(i_eta_cuts)->Write(name.c_str());


		}
		name = Form("%s",centrality.c_str());
		ltx->SetTextAlign(11);
		ltx->DrawLatex(0.19,0.88,name.c_str());

	}
	canvas1->cd(1);
	legend->Draw();

	name = "eff_trketa.pdf";
	canvas1->Print(name.c_str());

}

