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
	int jet_pt_end = jetpT_binning->FindBin(130);
	int trk_pt_start = trkpT_binning->FindBin(1);
	int trk_pt_end = trkpT_binning->FindBin(1.5);
	int run_start = 1;
	int run_end = N_runs;


	TH1* h_eventPercentage = new TH1D("h_eventPercentage","h_eventPercentage",N_runs,0,N_runs);
	h_eventPercentage->GetXaxis()->SetTitle("Run Number");
	h_eventPercentage->GetYaxis()->SetNdivisions(504);

	TH1* h_tmp = (TH1*)input_file->Get("EventPercentages");
	for (int i = 1; i <= h_eventPercentage->GetNbinsX(); i=i+2) h_eventPercentage->GetXaxis()->SetBinLabel(i,Form("%1.0f",h_tmp->GetBinLowEdge(i)));

	for (int i = 1; i <= h_eventPercentage->GetNbinsX(); i++)
	{
		h_eventPercentage->SetBinContent(i,h_tmp->GetBinContent(i));
	}


	TCanvas *c_x = new TCanvas("c_x","c_x",900,600);
	TCanvas *c_y = new TCanvas("c_y","c_y",900,600);

	h_eventPercentage->LabelsOption("v");
	h_eventPercentage->Draw("hist text");
	c_x->Print("run_dep/EventPercentages.pdf");

	TLegend *legend_x = new TLegend(0.20, 0.6, 0.40, 0.7, "","brNDC");
	legend_x->SetTextFont(43);
	legend_x->SetBorderSize(0);
	legend_x->SetTextSize(10);
	legend_x->SetNColumns(1);
	TLatex *ltx = new TLatex();
	ltx->SetTextFont(43);
	ltx->SetTextSize(12);
	ltx->SetTextAlign(32);



	vector<vector<vector<vector<TH1*>>>> h_UE_run_dep = vector<vector<vector<vector<TH1*>>>> (N_jetpt, vector<vector<vector<TH1*>>> (N_trkpt, vector<vector<TH1*>> (N_dR, vector<TH1*> (6)))); //histos from run_dep.c are such that h_XXX_run(i) where i is bin number (so starting at 1)
	vector<vector<TH1*>> h_jet_run_dep = vector<vector<TH1*>> (N_jetpt, vector<TH1*> (6)); //histos from run_dep.c are such that h_XXX_run(i) where i is bin number (so starting at 1)
	vector<vector<TH1*>> h_jet_run_dep_mc = vector<vector<TH1*>> (N_jetpt, vector<TH1*> (6)); //histos from run_dep.c are such that h_XXX_run(i) where i is bin number (so starting at 1)



	//initiaialize
	for (int i_jet = jet_pt_start-1; i_jet < jet_pt_end; i_jet++)
	{
		for (int i_cent = 0; i_cent < 6; i_cent++)
		{

			h_jet_run_dep[i_jet][i_cent] = (TH1*)h_eventPercentage->Clone(Form("jet_run_dep_jet%i_cent%i", i_jet, i_cent));
			h_jet_run_dep[i_jet][i_cent]->Reset();
			h_jet_run_dep[i_jet][i_cent]->GetYaxis()->SetTitle("dN/dp_{T} [Data]");

			h_jet_run_dep_mc[i_jet][i_cent] = (TH1*)h_eventPercentage->Clone(Form("jet_run_dep_mc_jet%i_cent%i", i_jet, i_cent));
			h_jet_run_dep_mc[i_jet][i_cent]->Reset();
			h_jet_run_dep_mc[i_jet][i_cent]->GetYaxis()->SetTitle("dN/dp_{T} [MC]");

			for (int i_trk = trk_pt_start-1; i_trk < trk_pt_end; i_trk++)
			{
				for (int i_dR = 0; i_dR < N_dR; i_dR++)
				{
					h_UE_run_dep[i_jet][i_trk][i_dR][i_cent] = (TH1*)h_eventPercentage->Clone(Form("UE_run_dep_jet%i_trk%i_dR%i_cent%i", i_jet, i_trk, i_dR, i_cent));
					h_UE_run_dep[i_jet][i_trk][i_dR][i_cent]->Reset();
					h_UE_run_dep[i_jet][i_trk][i_dR][i_cent]->GetYaxis()->SetTitle("UE_{TM Method}");
				}
			}

		}
	}

	double val, err;

	for (int i_jet = jet_pt_start-1; i_jet < jet_pt_end; i_jet++)
	{
		string jet_label = Form("%1.0f < p_{T}^{Jet} < %1.0f GeV", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));

		c_x->Clear();
		c_x->Divide(3,2);

		for (int i_cent = 0; i_cent < 6; i_cent++)
		{
			string cent_label = num_to_cent(31,i_cent).c_str();
			for (int i_run = 1; i_run <= N_runs; i_run++)
			{
				name = Form("data_jet_cent%i_run%i", i_cent,i_run);
				val = ((TH1*)input_file->Get(name.c_str()))->GetBinContent(i_jet+1);
				err = ((TH1*)input_file->Get(name.c_str()))->GetBinError(i_jet+1);
				h_jet_run_dep[i_jet][i_cent]->SetBinContent(i_run, val);
				h_jet_run_dep[i_jet][i_cent]->SetBinError(i_run, err);

				name = Form("TM_rN_norm_jet_cent%i_run%i", i_cent,i_run);
				val = ((TH1*)input_file->Get(name.c_str()))->GetBinContent(i_jet+1);
				err = ((TH1*)input_file->Get(name.c_str()))->GetBinError(i_jet+1);
				h_jet_run_dep_mc[i_jet][i_cent]->SetBinContent(i_run, val);
				h_jet_run_dep_mc[i_jet][i_cent]->SetBinError(i_run, err);
			}

			h_jet_run_dep[i_jet][i_cent]->Scale(1./h_jet_run_dep[i_jet][i_cent]->Integral());
			h_jet_run_dep_mc[i_jet][i_cent]->Scale(1./h_jet_run_dep_mc[i_jet][i_cent]->Integral());

			c_x->cd(i_cent+1);
			SetHStyle_smallify(h_jet_run_dep[i_jet][i_cent],0,1);
			SetHStyle_smallify(h_jet_run_dep_mc[i_jet][i_cent],1,1);


			h_jet_run_dep[i_jet][i_cent]->Draw("hist");
//			h_jet_run_dep_mc[i_jet][i_cent]->Draw("hist same");
//
//			if (i_jet == jet_pt_start-1 && i_cent == 0)
//			{
//				legend_x->AddEntry(h_jet_run_dep[i_jet][i_cent],"Data","lp");
//				legend_x->AddEntry(h_jet_run_dep_mc[i_jet][i_cent],"MC","lp");
//			}
//			legend_x->Draw();

			c_x->cd(i_cent+1);
			double x_position = 0.2;
			double y_position = 0.89;
			ltx->SetTextAlign(12);
			ltx->DrawLatexNDC(x_position,y_position,jet_label.c_str());
			ltx->DrawLatexNDC(x_position,y_position-0.05,cent_label.c_str());
		}

		name = "";
		if (i_jet == jet_pt_start-1) name = "(";
		if (i_jet == jet_pt_end-1) name = ")";

		c_x->Print(Form("run_dep/jet_run_dep.pdf%s",name.c_str()),Form("Title: jet%i", i_jet));

	}




	for (int i_jet = jet_pt_start-1; i_jet < jet_pt_end; i_jet++)
	{
		string jet_label = Form("%1.0f < p_{T}^{Jet} < %1.0f GeV", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));

		for (int i_trk = trk_pt_start-1; i_trk < trk_pt_end; i_trk++)
		{
			string trk_label = Form("%1.1f < p_{T}^{Trk} < %1.1f GeV", trkpT_binning->GetBinLowEdge(i_trk+1), trkpT_binning->GetBinUpEdge(i_trk+1));

			for (int i_dR = 0; i_dR < N_dR-2; i_dR++)
			{
				string dR_label = Form("%1.2f < r < %1.2f ", dR_binning->GetBinLowEdge(i_dR+1), dR_binning->GetBinUpEdge(i_dR+1));

				c_x->Clear();
				c_x->Divide(3,2);

//				c_y->Clear();
				for (int i_cent = 0; i_cent < 6; i_cent++)
				{
					string cent_label = num_to_cent(31,i_cent).c_str();


					//fill run dep histos here
					for (int i_run = 1; i_run <= N_runs; i_run++)
					{
						name = Form("jet%i_trk%i_cent%i", i_jet, i_trk, i_cent);
						name = Form("UE_TM_rN_indR_%s_run%i", name.c_str(), i_run);
						val = ((TH1*)input_file->Get(name.c_str()))->GetBinContent(i_dR+1);
						err = ((TH1*)input_file->Get(name.c_str()))->GetBinError(i_dR+1);

						h_UE_run_dep[i_jet][i_trk][i_dR][i_cent]->SetBinContent(i_run, val);
						h_UE_run_dep[i_jet][i_trk][i_dR][i_cent]->SetBinError(i_run, err);
					}

					c_x->cd(i_cent+1);
					h_UE_run_dep[i_jet][i_trk][i_dR][i_cent]->Draw("");

					double y_position = 0.24;
					double x_position = 0.92;

					ltx->SetTextAlign(32);
					ltx->DrawLatexNDC(x_position,y_position+0.21,cent_label.c_str());
					ltx->DrawLatexNDC(x_position,y_position+0.14,dR_label.c_str());
					ltx->DrawLatexNDC(x_position,y_position+0.07,trk_label.c_str());
					ltx->DrawLatexNDC(x_position,y_position,jet_label.c_str());


					if (i_cent == 0)
					{
//						c_y->cd();
//						h_UE_run_dep[i_jet][i_trk][i_dR][i_cent]->Draw("");

//						double y_position = 0.24;
//						double x_position = 0.92;
//
//						ltx->SetTextAlign(32);
//						ltx->DrawLatexNDC(x_position,y_position+0.21,cent_label.c_str());
//						ltx->DrawLatexNDC(x_position,y_position+0.14,dR_label.c_str());
//						ltx->DrawLatexNDC(x_position,y_position+0.07,trk_label.c_str());
//						ltx->DrawLatexNDC(x_position,y_position,jet_label.c_str());
//
//						name = "";
//						if (i_trk == trk_pt_start-1 && i_jet == jet_pt_start-1 && i_dR == 0) name = "(";
//						if (i_trk == trk_pt_end-1 && i_jet == jet_pt_end-1 && i_dR == N_dR-3) name = ")";
//						c_y->Print(Form("run_dep/chps_UE_run_dep_cent0.pdf%s",name.c_str()),Form("Title: jet%i_trk%i_dR%i", i_jet, i_trk, i_dR));

					}
				}

				name = "";
				if (i_trk == trk_pt_start-1 && i_jet == jet_pt_start-1 && i_dR == 0) name = "(";
				if (i_trk == trk_pt_end-1 && i_jet == jet_pt_end-1 && i_dR == N_dR-3) name = ")";
				c_x->Print(Form("run_dep/chps_UE_run_dep.pdf%s",name.c_str()),Form("Title: jet%i_trk%i_dR%i", i_jet, i_trk, i_dR));

			}
		}
	}




/*
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

			}
			legend_x->Draw();

			name = "";
			if (i_trk == trk_pt_start-1 && i_jet == jet_pt_start-1) name = "(";
			if (i_trk == trk_pt_end-1 && i_jet == jet_pt_end-1) name = ")";
			c_x->Print(Form("run_dep/chps_UE_run_dep.pdf%s",name.c_str()),Form("Title: jet%i_trk%i", i_jet, i_trk));
		}
	}
*/


}