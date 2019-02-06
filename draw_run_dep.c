#include "output_dev/functions/global_variables.h"

void draw_run_dep(int config = 37)
{
	SetAtlasStyle();
	gErrorIgnoreLevel = 3001;
	double r_max_range = 0.8;

	string name = Form("./run_dep/run_dep_UE_c%i.root", config);
	TFile *input_file = new TFile(name.c_str());

	TAxis* run_binning = (TAxis*)input_file->Get("run_binning");
	TAxis* dR_binning = (TAxis*)input_file->Get("dR_binning");
	TAxis* jetpT_binning = (TAxis*)input_file->Get("jetpT_binning");
	TAxis* trkpT_binning = (TAxis*)input_file->Get("trkpT_binning");

	int N_dR = dR_binning->GetNbins();
	int N_jetpt = jetpT_binning->GetNbins();
	int N_trkpt = trkpT_binning->GetNbins();
	int N_runs = run_binning->GetNbins();

	int jet_pt_start = jetpT_binning->FindBin(127);
	int jet_pt_end = jetpT_binning->FindBin(315);
	int trk_pt_start = trkpT_binning->FindBin(1);
	int trk_pt_end = trkpT_binning->FindBin(9.99);
	int run_start = 15;
	int run_end = 23;


	bool doRuns[30];
	for (int i = 0; i < 30; i++) { doRuns[i] = 1;}
	double eventThreshold = 0.03;

	TH1* h_eventPercentage = (TH1*)input_file->Get("EventPercentages");
	for (int i = 1; i < h_eventPercentage->GetNbinsX(); i++)
	{
		if (h_eventPercentage->GetBinContent(i) < eventThreshold) doRuns[i] = 0;
	}


	TCanvas *c_x = new TCanvas("c_x","c_x",900,600);

	TLegend *legend_x = new TLegend(0.60, 0.45, 0.90, 0.75, "","brNDC");
	legend_x->SetTextFont(43);
	legend_x->SetBorderSize(0);
	legend_x->SetTextSize(10);
	legend_x->SetNColumns(1);
	TLatex *ltx = new TLatex();
	ltx->SetTextFont(43);
	ltx->SetTextSize(12);
	ltx->SetTextAlign(32);



	vector<TH1*> h_UE_run_dep = vector<TH1*> (N_runs+1); //histos from run_dep.c are such that h_XXX_run(i) where i is bin number (so starting at 1)
	vector<TH1*> h_UE_TM = vector<TH1*> (N_runs+1); //histos from run_dep.c are such that h_XXX_run(i) where i is bin number (so starting at 1)
	vector<TH1*> h_ratio = vector<TH1*> (N_runs+1); //histos from run_dep.c are such that h_XXX_run(i) where i is bin number (so starting at 1)

	


	for (int i_jet = jet_pt_start-1; i_jet < jet_pt_end; i_jet++)
	{
		for (int i_trk = trk_pt_start-1; i_trk < trk_pt_end; i_trk++)
		{
			string trk_label = Form("%1.1f < p_{T}^{Trk} < %1.1f GeV", trkpT_binning->GetBinLowEdge(i_trk+1), trkpT_binning->GetBinUpEdge(i_trk+1));
			string jet_label = Form("%1.0f < p_{T}^{Jet} < %1.0f GeV", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));

			c_x->cd();
			c_x->Clear();
			c_x->Divide(3,2);


			for (int i_cent = 0; i_cent < 6; i_cent++)
			{
				string cent_label = num_to_cent(31,i_cent).c_str();

				int run_itr = 1;


				for (int i_run = run_start; i_run <= run_end; i_run++)
				{
					name = Form("jet%i_trk%i_cent%i", i_jet, i_trk, i_cent);
					name = Form("UE_TM_indR_%s_run%i", name.c_str(), i_run);
					h_UE_TM[i_run] = (TH1*)input_file->Get(name.c_str())->Clone(Form("%s_clone",name.c_str()));

					name = Form("jet%i_trk%i_cent%i", i_jet, i_trk, i_cent);
					name = Form("UE_TM_rN_indR_%s_run%i", name.c_str(), i_run);
					h_UE_run_dep[i_run] = (TH1*)input_file->Get(name.c_str())->Clone(Form("%s_clone",name.c_str()));

					h_ratio[i_run] = (TH1*)h_UE_run_dep[i_run]->Clone(Form("ratio_run%i",i_run));
					h_ratio[i_run]->Divide(h_UE_TM[i_run]);

					if (!doRuns[i_run]) continue;

					SetHStyle_smallify(h_UE_TM[i_run], 0, 1);
					SetHStyle_smallify(h_UE_run_dep[i_run], run_itr, 1);
					SetHStyle_smallify(h_ratio[i_run], run_itr, 1);

					if (i_trk == trk_pt_start-1 && i_jet == jet_pt_start-1 && i_cent == 0)
					{
						if (run_itr==1) legend_x->AddEntry(h_UE_TM[i_run], "Combined", "lp");
						legend_x->AddEntry(h_UE_run_dep[i_run], Form("run %1.0f (%i)",run_binning->GetBinLowEdge(i_run), i_run), "lp");
					}

					double low_range = h_UE_TM[i_run]->GetMinimum() * 0.75;
					double hi_range = h_UE_TM[i_run]->GetMaximum() * 1.4;


					h_UE_TM[i_run]->GetXaxis()->SetRangeUser(0,r_max_range);
					h_UE_TM[i_run]->GetYaxis()->SetNdivisions(504);

					h_ratio[i_run]->GetXaxis()->SetRangeUser(0,r_max_range);
					h_ratio[i_run]->GetYaxis()->SetRangeUser(0.8,1.2);
					h_ratio[i_run]->GetYaxis()->SetNdivisions(504);

					if (run_itr == 1)
					{
						c_x->cd(i_cent+1);
						gPad->Divide(1,2);

						c_x->cd(i_cent+1)->cd(1);
						gPad->SetPad(0,0.40,0.99,0.99);
						gPad->SetTopMargin(0.05);
						gPad->SetBottomMargin(0.0);
						gPad->SetRightMargin(0);

						c_x->cd(i_cent+1)->cd(2);
						gPad->SetPad(0,0.0,0.99,0.40);
						gPad->SetTopMargin(0.0);
						gPad->SetBottomMargin(0.30);
						gPad->SetRightMargin(0);
					}

					c_x->cd(i_cent+1)->cd(1);

					if (run_itr == 1) h_UE_TM[i_run]->DrawCopy();
					h_UE_run_dep[i_run]->DrawCopy("same");

					c_x->cd(i_cent+1)->cd(2);
					if (run_itr == 1) h_ratio[i_run]->DrawCopy();
					else h_ratio[i_run]->DrawCopy("same");

					run_itr++;
				}

				c_x->cd(i_cent+1);
				ltx->DrawLatexNDC(0.98,0.94,cent_label.c_str());
				ltx->DrawLatexNDC(0.98,0.87,trk_label.c_str());
				ltx->DrawLatexNDC(0.98,0.80,jet_label.c_str());
			}
			legend_x->Draw();

			name = "";
			if (i_trk == trk_pt_start-1 && i_jet == jet_pt_start-1) name = "(";
			if (i_trk == trk_pt_end-1 && i_jet == jet_pt_end-1) name = ")";
			c_x->Print(Form("run_dep/chps_UE_run_dep.pdf%s",name.c_str()),Form("Title: jet%i_trk%i", i_jet, i_trk));
		}
	}



}
