#include "output_dev/functions/global_variables.h"

void draw_run_dep(int config = 38)
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

	map<int, double> luminosity = {{ 286665 , 0.02071 }, { 286711 , 0.419829 }, { 286717 , 0.590276 }, { 286748 , 4.24438 }, { 286767 , 5.799 }, { 286834 , 13.1028 }, { 286854 , 13.595759 }, { 286908 , 11.00192 }, { 286990 , 10.27951 }, { 287038 , 15.9825 }, { 287044 , 23.3479 }, { 287068 , 6.83118 }, { 287222 , 1.33018 }, { 287224 , 1.932658 }, { 287259 , 17.4035 }, { 287270 , 22.0425 }, { 287281 , 24.1107 }, { 287321 , 5.49172 }, { 287330 , 21.620117 }, { 287334 , 16.3135 }, { 287378 , 16.81257 }, { 287380 , 0.319641 }, { 287382 , 17.651474 }, { 287560 , 0.572368 }, { 287594 , 12.6278 }, { 287632 , 18.9915 }, { 287706 , 26.6057 }, { 287728 , 25.93 }, { 287827 , 24.1147 }, { 287843 , 25.0225 }, { 287866 , 42.1866 }, { 287924 , 22.5426 }, { 287931 , 37.2019 }};


	TH1* h_eventPercentage_mc = new TH1D("h_eventPercentage_mc","h_eventPercentage_mc",N_runs,0,N_runs);
	h_eventPercentage_mc->GetXaxis()->SetTitle("Run Number");
	h_eventPercentage_mc->GetYaxis()->SetTitle("Event Fraction");
	h_eventPercentage_mc->GetYaxis()->SetNdivisions(504);

	TH1* h_eventPercentage_data = new TH1D("h_eventPercentage_data","h_eventPercentage_data",N_runs,0,N_runs);
	h_eventPercentage_data->GetXaxis()->SetTitle("Run Number");
	h_eventPercentage_data->GetYaxis()->SetTitle("Event Fraction");
	h_eventPercentage_data->GetYaxis()->SetNdivisions(504);

	TH1* h_tmp_mc = (TH1*)input_file->Get("EventPercentages_mc");
	TH1* h_tmp_data = (TH1*)input_file->Get("EventPercentages_data");
	for (int i = 1; i <= h_eventPercentage_data->GetNbinsX(); i=i+1)
	{
		h_eventPercentage_mc->GetXaxis()->SetBinLabel(i,Form("%1.0f",h_tmp_mc->GetBinLowEdge(i)));
		h_eventPercentage_data->GetXaxis()->SetBinLabel(i,Form("%1.0f",h_tmp_data->GetBinLowEdge(i)));
	}

	for (int i = 1; i <= h_eventPercentage_mc->GetNbinsX(); i++)
	{
		h_eventPercentage_mc->SetBinContent(i,h_tmp_mc->GetBinContent(i));
		h_eventPercentage_data->SetBinContent(i,h_tmp_data->GetBinContent(i));
	}

	double j = 0;

	for (int i = 1; i <= h_eventPercentage_data->GetXaxis()->GetNbins(); i++)
	{
		j = j + h_eventPercentage_data->GetBinContent(i);
	}
	cout << j <<  " <-------------- 1" << endl;


	h_eventPercentage_mc->Scale(1./h_eventPercentage_mc->Integral());
	h_eventPercentage_data->Scale(1./h_eventPercentage_data->Integral());

	j = 0;

	for (int i = 1; i <= h_eventPercentage_data->GetXaxis()->GetNbins(); i++)
	{
		j = j + h_eventPercentage_data->GetBinContent(i);
	}
	cout << j <<  " <-------------- 2" << endl;



	TCanvas *c_x = new TCanvas("c_x","c_x",1200,600);
	TCanvas *c_y = new TCanvas("c_y","c_y",900,600);
	TLegend *legend_x = new TLegend(0.20, 0.6, 0.40, 0.7, "","brNDC");
	legend_x->SetTextFont(43);
	legend_x->SetBorderSize(0);
	legend_x->SetTextSize(10);
	legend_x->SetNColumns(1);
	TLatex *ltx = new TLatex();
	ltx->SetTextFont(43);
	ltx->SetTextSize(12);
	ltx->SetTextAlign(32);


	c_x->cd();
	SetHStyle_smallify(h_eventPercentage_data, 0, 1);
	SetHStyle_smallify(h_eventPercentage_mc, 1, 1);
	h_eventPercentage_mc->LabelsOption("v");
	h_eventPercentage_data->LabelsOption("v");
	h_eventPercentage_mc->Draw("hist");
	h_eventPercentage_data->Draw("hist same");
	legend_x->AddEntry(h_eventPercentage_data,"Data","lp");
	legend_x->AddEntry(h_eventPercentage_mc,"MC","lp");
	legend_x->Draw();


	c_x->Print("run_dep/EventPercentages.pdf");

	vector<vector<vector<vector<TH1*>>>> h_UE_TM_run_dep = vector<vector<vector<vector<TH1*>>>> (N_jetpt, vector<vector<vector<TH1*>>> (N_trkpt, vector<vector<TH1*>> (N_dR, vector<TH1*> (6)))); //histos from run_dep.c are such that h_XXX_run(i) where i is bin number (so starting at 1)
	vector<vector<vector<vector<TF1*>>>> func_UE_TM_run_dep = vector<vector<vector<vector<TF1*>>>> (N_jetpt, vector<vector<vector<TF1*>>> (N_trkpt, vector<vector<TF1*>> (N_dR, vector<TF1*> (6)))); //histos from run_dep.c are such that h_XXX_run(i) where i is bin number (so starting at 1)

	vector<vector<vector<vector<TH1*>>>> h_UE_MB_run_dep = vector<vector<vector<vector<TH1*>>>> (N_jetpt, vector<vector<vector<TH1*>>> (N_trkpt, vector<vector<TH1*>> (N_dR, vector<TH1*> (6))));
	vector<vector<vector<vector<TF1*>>>> func_UE_MB_run_dep = vector<vector<vector<vector<TF1*>>>> (N_jetpt, vector<vector<vector<TF1*>>> (N_trkpt, vector<vector<TF1*>> (N_dR, vector<TF1*> (6))));

	vector<vector<TH1*>> h_jet_run_dep_data = vector<vector<TH1*>> (N_jetpt, vector<TH1*> (6));
	vector<vector<TH1*>> h_jet_run_dep_mc = vector<vector<TH1*>> (N_jetpt, vector<TH1*> (6));


	vector<vector<vector<TH1*>>> jet_weighted_UE_MB = vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_trkpt, vector<TH1*> (6)));
	vector<vector<vector<TH1*>>> evt_weighted_UE_MB = vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_trkpt, vector<TH1*> (6)));


	name = Form("jet%i_trk%i_cent%i", 8, 4, 0);
	name = Form("UE_MB_indR_%s_run%i", name.c_str(), 1);
	TH1* h_tmp = (TH1*)input_file->Get(name.c_str())->Clone("h_tmp");


	//initiaialize
	for (int i_jet = jet_pt_start-1; i_jet < jet_pt_end; i_jet++)
	{
		for (int i_cent = 0; i_cent < 6; i_cent++)
		{
			h_jet_run_dep_data[i_jet][i_cent] = (TH1*)h_eventPercentage_data->Clone(Form("jet_run_dep_data_jet%i_cent%i", i_jet, i_cent));
			h_jet_run_dep_data[i_jet][i_cent]->Reset();
			h_jet_run_dep_data[i_jet][i_cent]->GetYaxis()->SetTitle("dN/dp_{T} [Data]");

			h_jet_run_dep_mc[i_jet][i_cent] = (TH1*)h_eventPercentage_mc->Clone(Form("jet_run_dep_mc_jet%i_cent%i", i_jet, i_cent));
			h_jet_run_dep_mc[i_jet][i_cent]->Reset();
			h_jet_run_dep_mc[i_jet][i_cent]->GetYaxis()->SetTitle("dN/dp_{T} [MC]");

			
			for (int i_trk = trk_pt_start-1; i_trk < trk_pt_end; i_trk++)
			{
				for (int i_dR = 0; i_dR < N_dR; i_dR++)
				{
					h_UE_TM_run_dep[i_jet][i_trk][i_dR][i_cent] = (TH1*)h_eventPercentage_mc->Clone(Form("UE_TM_run_dep_jet%i_trk%i_dR%i_cent%i", i_jet, i_trk, i_dR, i_cent));
					h_UE_TM_run_dep[i_jet][i_trk][i_dR][i_cent]->Reset();
					h_UE_TM_run_dep[i_jet][i_trk][i_dR][i_cent]->GetYaxis()->SetTitle("UE_TM_{TM Method}");

					h_UE_MB_run_dep[i_jet][i_trk][i_dR][i_cent] = (TH1*)h_eventPercentage_mc->Clone(Form("UE_MB_run_dep_jet%i_trk%i_dR%i_cent%i", i_jet, i_trk, i_dR, i_cent));
					h_UE_MB_run_dep[i_jet][i_trk][i_dR][i_cent]->Reset();
					h_UE_MB_run_dep[i_jet][i_trk][i_dR][i_cent]->GetYaxis()->SetTitle("UE_{MB Method}");
				}

				jet_weighted_UE_MB[i_jet][i_trk][i_cent] = (TH1*)h_tmp->Clone(Form("jet_w_UE_MB_jet%i_trk%i_cent%i", i_jet, i_trk, i_cent));
				jet_weighted_UE_MB[i_jet][i_trk][i_cent]->Reset();
				jet_weighted_UE_MB[i_jet][i_trk][i_cent]->GetYaxis()->SetTitle("UE_TM_{MB Method}^{JetW}");

				evt_weighted_UE_MB[i_jet][i_trk][i_cent]= (TH1*)h_tmp->Clone(Form("evt_w_UE_MB_jet%i_trk%i_cent%i", i_jet, i_trk, i_cent));
				evt_weighted_UE_MB[i_jet][i_trk][i_cent]->Reset();
				evt_weighted_UE_MB[i_jet][i_trk][i_cent]->GetYaxis()->SetTitle("UE_{MB Method}^{EvtW}");


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
				int run_number = h_tmp_mc->GetBinLowEdge(i_run);
				double lumi = 1;//luminosity[run_number];
				if (lumi < 0.0001) continue; //this ensures run is in GRL and has a finite non zero luminosity

				name = Form("MB_rN_norm_jet_cent%i_run%i", i_cent,i_run);
				val = ((TH1*)input_file->Get(name.c_str()))->GetBinContent(i_jet+1) / lumi;
				err = ((TH1*)input_file->Get(name.c_str()))->GetBinError(i_jet+1) / lumi;
				h_jet_run_dep_data[i_jet][i_cent]->SetBinContent(i_run, val);
				h_jet_run_dep_data[i_jet][i_cent]->SetBinError(i_run, err);

				name = Form("TM_rN_norm_jet_cent%i_run%i", i_cent,i_run);
				val = ((TH1*)input_file->Get(name.c_str()))->GetBinContent(i_jet+1) / lumi;
				err = ((TH1*)input_file->Get(name.c_str()))->GetBinError(i_jet+1) / lumi;
				h_jet_run_dep_mc[i_jet][i_cent]->SetBinContent(i_run, val);
				h_jet_run_dep_mc[i_jet][i_cent]->SetBinError(i_run, err);
			}

			h_jet_run_dep_data[i_jet][i_cent]->Scale(1./h_jet_run_dep_data[i_jet][i_cent]->Integral());
			h_jet_run_dep_mc[i_jet][i_cent]->Scale(1./h_jet_run_dep_mc[i_jet][i_cent]->Integral());

			c_x->cd(i_cent+1);
			SetHStyle_smallify(h_jet_run_dep_data[i_jet][i_cent],0,1);
			SetHStyle_smallify(h_jet_run_dep_mc[i_jet][i_cent],1,1);

			h_jet_run_dep_mc[i_jet][i_cent]->Draw("hist");
			h_jet_run_dep_data[i_jet][i_cent]->Draw("hist same");

			if (i_jet == jet_pt_start-1 && i_cent == 0)
			{
				legend_x->Clear();
				legend_x->AddEntry(h_jet_run_dep_data[i_jet][i_cent],"Data","lp");
				legend_x->AddEntry(h_jet_run_dep_mc[i_jet][i_cent],"MC","lp");
			}
			legend_x->Draw();

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

				c_y->Clear();
				c_y->Divide(3,2);
				for (int i_cent = 0; i_cent < 6; i_cent++)
				{
					string cent_label = num_to_cent(31,i_cent).c_str();

					double jet_w_UE = 0, sum_jet_w = 0;
					double evt_w_UE = 0, sum_evt_w = 0;

					//fill run dep histos here
					for (int i_run = 1; i_run <= N_runs; i_run++)
					{
						int run_number = h_tmp_mc->GetBinLowEdge(i_run);
						double lumi = luminosity[run_number];
						double jet_weights = h_jet_run_dep_data[i_jet][i_cent]->GetBinContent(i_run);
						double evt_weights = h_eventPercentage_data->GetBinContent(i_run);

						if (lumi < 0.001)  continue;

						name = Form("jet%i_trk%i_cent%i", i_jet, i_trk, i_cent);
						name = Form("UE_TM_rN_indR_%s_run%i", name.c_str(), i_run);
						val = ((TH1*)input_file->Get(name.c_str()))->GetBinContent(i_dR+1);
						err = ((TH1*)input_file->Get(name.c_str()))->GetBinError(i_dR+1);
						h_UE_TM_run_dep[i_jet][i_trk][i_dR][i_cent]->SetBinContent(i_run, val);
						h_UE_TM_run_dep[i_jet][i_trk][i_dR][i_cent]->SetBinError(i_run, err);

						name = Form("jet%i_trk%i_cent%i", i_jet, i_trk, i_cent);
						name = Form("UE_MB_rN_indR_%s_run%i", name.c_str(), i_run);
						val = ((TH1*)input_file->Get(name.c_str()))->GetBinContent(i_dR+1);
						err = ((TH1*)input_file->Get(name.c_str()))->GetBinError(i_dR+1);
						h_UE_MB_run_dep[i_jet][i_trk][i_dR][i_cent]->SetBinContent(i_run, val);
						h_UE_MB_run_dep[i_jet][i_trk][i_dR][i_cent]->SetBinError(i_run, err);

						jet_w_UE = jet_w_UE + jet_weights*val;
						evt_w_UE = evt_w_UE + evt_weights*val;
						sum_jet_w = sum_jet_w+jet_weights;
						sum_evt_w = sum_evt_w+evt_weights;

//						cout << Form("%1.2f * %1.5f = %1.5f, Running Total: %1.5f, Running SumW: %1.4f", val, jet_weights, jet_weights*val, jet_w_UE, sum_jet_w) << endl;
//						cout << Form("%1.2f * %1.5f = %1.5f, Running Total: %1.5f, Running SumW: %1.4f", val, evt_weights, evt_weights*val, evt_w_UE, sum_evt_w) << endl;

					}
//					cout <<  "--------------------" << endl;
//					h_eventPercentage_data->Print("all");
//					h_jet_run_dep_data[i_jet][i_cent]->Print("all");
//					h_UE_MB_run_dep[i_jet][i_trk][i_dR][i_cent]->Print("all");
//					cout <<  "--------------------" << endl;




					jet_weighted_UE_MB[i_jet][i_trk][i_cent]->SetBinContent(i_dR+1, jet_w_UE/sum_jet_w);
					evt_weighted_UE_MB[i_jet][i_trk][i_cent]->SetBinContent(i_dR+1, evt_w_UE/sum_evt_w);
					jet_weighted_UE_MB[i_jet][i_trk][i_cent]->Sumw2();//(i_dR+1, jet_w_UE/sum_jet_w);
					evt_weighted_UE_MB[i_jet][i_trk][i_cent]->Sumw2();//(i_dR+1, jet_w_UE/sum_jet_w);


					c_x->cd(i_cent+1);
					SetHStyle_smallify(h_UE_TM_run_dep[i_jet][i_trk][i_dR][i_cent],0,1);
					SetHStyle_smallify(h_UE_MB_run_dep[i_jet][i_trk][i_dR][i_cent],1,1);
					h_UE_TM_run_dep[i_jet][i_trk][i_dR][i_cent]->Draw("");
					h_UE_MB_run_dep[i_jet][i_trk][i_dR][i_cent]->Draw("same hist");

					double y_position = 0.24;
					double x_position = 0.92;

					ltx->SetTextAlign(32);
					ltx->SetTextSize(12);
					ltx->DrawLatexNDC(x_position,y_position+0.21,cent_label.c_str());
					ltx->DrawLatexNDC(x_position,y_position+0.14,dR_label.c_str());
					ltx->DrawLatexNDC(x_position,y_position+0.07,trk_label.c_str());
					ltx->DrawLatexNDC(x_position,y_position,jet_label.c_str());

					if (i_dR == N_dR - 3)
					{
						c_y->cd(i_cent+1);
						SetHStyle_smallify(jet_weighted_UE_MB[i_jet][i_trk][i_cent],2,1);
						SetHStyle_smallify(evt_weighted_UE_MB[i_jet][i_trk][i_cent],1,1);
						cout << "############" << endl;
						name = Form("jet%i_trk%i_cent%i", i_jet, i_trk, i_cent);
						name = Form("UE_TM_indR_%s_run%i", name.c_str(), 1);
						((TH1*)input_file->Get(name.c_str()))->Print("all");
						jet_weighted_UE_MB[i_jet][i_trk][i_cent]->Print("all");
						evt_weighted_UE_MB[i_jet][i_trk][i_cent]->Print("all");
						cout << "############" << endl;

						jet_weighted_UE_MB[i_jet][i_trk][i_cent]->Divide(((TH1*)input_file->Get(name.c_str())));
//						evt_weighted_UE_MB[i_jet][i_trk][i_cent]->Divide(((TH1*)input_file->Get(name.c_str())));
//						evt_weighted_UE_MB[i_jet][i_trk][i_cent]->Divide(jet_weighted_UE_MB[i_jet][i_trk][i_cent]);

//						((TH1*)input_file->Get(name.c_str()))->Draw("");
//						evt_weighted_UE_MB[i_jet][i_trk][i_cent]->Draw("hist");
						jet_weighted_UE_MB[i_jet][i_trk][i_cent]->Draw("p");
						jet_weighted_UE_MB[i_jet][i_trk][i_cent]->GetYaxis()->SetRangeUser(0.95,1.05);

						y_position = 0.24;
						x_position = 0.92;

						ltx->SetTextAlign(32);
						ltx->SetTextSize(12);
						ltx->DrawLatexNDC(x_position,y_position+0.21,cent_label.c_str());
						ltx->DrawLatexNDC(x_position,y_position+0.14,dR_label.c_str());
						ltx->DrawLatexNDC(x_position,y_position+0.07,trk_label.c_str());
						ltx->DrawLatexNDC(x_position,y_position,jet_label.c_str());
					}


//					if (i_cent == 0)
//					{
//						c_y->cd();
//						h_UE_run_dep[i_jet][i_trk][i_dR][i_cent]->GetYaxis()->SetRangeUser(0,100);
//						h_UE_run_dep[i_jet][i_trk][i_dR][i_cent]->Draw("");
//
//						func_UE_run_dep[i_jet][i_trk][i_dR][i_cent] = new TF1(Form("f_c%i",i_cent), "pol1");// h_UE_run_dep[i_jet][i_trk][i_dR][i_cent]->GetFunction();
//						func_UE_run_dep[i_jet][i_trk][i_dR][i_cent]->SetLineColor(kBlue);
//						h_UE_run_dep[i_jet][i_trk][i_dR][i_cent]->Fit(Form("f_c%i",i_cent),"RQ","",4,36);
//
//						double chi2 = func_UE_run_dep[i_jet][i_trk][i_dR][i_cent]->GetChisquare();
//						double NDF = func_UE_run_dep[i_jet][i_trk][i_dR][i_cent]->GetNDF();
//						double par0 = func_UE_run_dep[i_jet][i_trk][i_dR][i_cent]->GetParameter(0);
//						double err0 = func_UE_run_dep[i_jet][i_trk][i_dR][i_cent]->GetParError(0);
//						double par1 = func_UE_run_dep[i_jet][i_trk][i_dR][i_cent]->GetParameter(1);
//						double err1 = func_UE_run_dep[i_jet][i_trk][i_dR][i_cent]->GetParError(1);
//
//						string equation = Form("Fit: %1.2f(#pm %1.4f) + %1.4f(#pm %1.4f)x",par0, err0, par1, err1);
//						string chi2ndf = Form("#chi^{2}/NDF : %1.4f",chi2/NDF);
//
//
//						ltx->SetTextSize(16);
//						double y_position = 0.24;
//						double x_position = 0.92;
//
//						ltx->SetTextAlign(12);
//						ltx->SetTextColor(func_UE_run_dep[i_jet][i_trk][i_dR][i_cent]->GetLineColor());
//						ltx->DrawLatexNDC(0.2,y_position+0.07,equation.c_str());
//						ltx->DrawLatexNDC(0.2,y_position,chi2ndf.c_str());
//
//						ltx->SetTextAlign(32);
//						ltx->SetTextColor(kBlack);
//						ltx->DrawLatexNDC(x_position,y_position+0.21,cent_label.c_str());
//						ltx->DrawLatexNDC(x_position,y_position+0.14,dR_label.c_str());
//						ltx->DrawLatexNDC(x_position,y_position+0.07,trk_label.c_str());
//						ltx->DrawLatexNDC(x_position,y_position,jet_label.c_str());
//
//						name = "";
//						if (i_trk == trk_pt_start-1 && i_jet == jet_pt_start-1 && i_dR == 0) name = "(";
//						if (i_trk == trk_pt_end-1 && i_jet == jet_pt_end-1 && i_dR == N_dR-3) name = ")";
//						c_y->Print(Form("run_dep/chps_UE_run_dep_cent0.pdf%s",name.c_str()),Form("Title: jet%i_trk%i_dR%i", i_jet, i_trk, i_dR));
//
//					}
				}

				name = "";
				if (i_trk == trk_pt_start-1 && i_jet == jet_pt_start-1 && i_dR == 0) name = "(";
				if (i_trk == trk_pt_end-1 && i_jet == jet_pt_end-1 && i_dR == N_dR-3) name = ")";
				c_x->Print(Form("run_dep/chps_UE_run_dep.pdf%s",name.c_str()),Form("Title: jet%i_trk%i_dR%i", i_jet, i_trk, i_dR));
				c_y->Print(Form("run_dep/chps_w_UE.pdf%s",name.c_str()),Form("Title: jet%i_trk%i_dR%i", i_jet, i_trk, i_dR));


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
