#include "output_dev/functions/global_variables.h"

void run_dep(int sys_mode = 37, bool subtract = 0)
{
	SetAtlasStyle();
	gErrorIgnoreLevel = 3001;
	string name;

	string sys_path;
	double r_max_range = 0.8;
	if (sys_mode == 0) sys_path = Form("nominal");
	if (sys_mode > 0 && sys_mode < 100) sys_path = Form("c%i", sys_mode);
	if (sys_mode > 100) sys_path = Form("sys%i", sys_mode);

//	TFile *input_file = new TFile(Form("hist-local_mc.root", sys_mode));
//	TFile *input_file_data = new TFile(Form("hist-local_mc.root", sys_mode));
	TFile *input_file = new TFile(Form("output_dev/raw_results/%s/FF_MC_out_histo_PbPb_5p02_r001.root", sys_path.c_str()));

	TH1* h_tmp = (TH1*)input_file->Get("h_event_rN");
	h_tmp->Scale(1./h_tmp->Integral());

	double sum = 0;
	for (int i = 1; i <= h_tmp->GetXaxis()->GetNbins() ;i++)
	{
		if (h_tmp->GetBinContent(i) < 0.03) continue;
//		cout << h_tmp->GetBinLowEdge(i) << " " << i << " " << h_tmp->GetBinContent(i) << endl;
		sum += h_tmp->GetBinContent(i);
	}

	name = Form("./run_dep/run_dep_UE_%s.root", sys_path.c_str());
	TFile *output_file = new TFile(name.c_str(), "recreate");

	cout << "Using file " << input_file->GetName() << endl;
	output_file->cd();
	h_tmp->Write("EventPercentages");

	TAxis* dR_binning = (TAxis*)((TH3*)input_file->Get("h_dR_change_jetpt0_cent0"))->GetXaxis();
	TAxis* jetpT_binning = (TAxis*)((TH3*)input_file->Get("ChPS_raw_0_dR0_cent0"))->GetYaxis();
	TAxis* trkpT_binning = (TAxis*)((TH3*)input_file->Get("ChPS_raw_0_dR0_cent0"))->GetXaxis();
	TAxis* run_binning = (TAxis*)((TH3*)input_file->Get("h_event_rN"))->GetXaxis();

	output_file->cd();
	dR_binning->Write("dR_binning");
	jetpT_binning->Write("jetpT_binning");
	trkpT_binning->Write("trkpT_binning");
	run_binning->Write("run_binning");

	int N_dR = dR_binning->GetNbins();
	int N_jetpt = jetpT_binning->GetNbins();
	int N_trkpt = trkpT_binning->GetNbins();
	int N_runs = run_binning->GetNbins();

	int jet_pt_start = jetpT_binning->FindBin(127);
	int jet_pt_end = jetpT_binning->FindBin(315);
	int trk_pt_start = trkpT_binning->FindBin(1);
	int trk_pt_end = trkpT_binning->FindBin(9.99);
	int run_start = 1;
	int run_end = N_runs;

	//2D UE
	vector<vector<TH2*>> h_TM_2D = vector<vector<TH2*>> (N_dR, vector<TH2*> (n_cent_cuts));
	vector<vector<TH2*>> h_TM_rN_2D = vector<vector<TH2*>> (N_dR, vector<TH2*> (n_cent_cuts));

	//1D UE
	vector<vector<vector<TH1*>>> h_TM_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_dR, vector<TH1*> (n_cent_cuts)));
	vector<vector<vector<TH1*>>> h_TM_rN_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_dR, vector<TH1*> (n_cent_cuts)));

	//2D UE injet
	vector<TH2*> h_TM_injet_2D = vector<TH2*> (n_cent_cuts);
	vector<TH2*> h_TM_rN_injet_2D = vector<TH2*> (n_cent_cuts);

	//1D UE u nr
	vector<vector<vector<TH1*>>> h_TM_r_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_trkpt, vector<TH1*> (n_cent_cuts)));
	vector<vector<vector<TH1*>>> h_TM_rN_r_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_trkpt, vector<TH1*> (n_cent_cuts)));

	//2D ratio
	vector<vector<TH2*>> h_TM_rN_TM_2D = vector<vector<TH2*>> (N_dR, vector<TH2*> (n_cent_cuts));

	//1D ratio
	vector<vector<vector<TH1*>>> h_TM_rN_TM_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_dR, vector<TH1*> (n_cent_cuts)));

	//1D ratio in r
	vector<vector<vector<TH1*>>> h_TM_rN_TM_r_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_trkpt, vector<TH1*> (n_cent_cuts)));



	TCanvas *c_run_dep_jets = new TCanvas("c_run_dep_jets","c_run_dep_jets",900,600);
	TLegend *legend_jet_run_dep = new TLegend();
	legend_jet_run_dep->SetTextFont(43);
	legend_jet_run_dep->SetBorderSize(0);
	legend_jet_run_dep->SetTextSize(10);
	legend_jet_run_dep->SetNColumns(1);

	TLine *line = new TLine();
	TLatex *ltx = new TLatex();
	ltx->SetTextFont(43);
	ltx->SetTextSize(12);
	ltx->SetTextAlign(32);


	double array_dr_bins[N_dR+1];
	for (int i_dR = 0; i_dR <= N_dR; i_dR++) array_dr_bins[i_dR] = dR_binning->GetBinLowEdge(i_dR+1);

	for (int i_cent = 0; i_cent < n_cent_cuts; i_cent++)
	{
		for (int i_jet = 0; i_jet < N_jetpt; i_jet++)
		{
			for (int i_trk = 0; i_trk < N_trkpt; i_trk++)
			{
				name = Form("h_TM_r_1D_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_TM_r_1D.at(i_jet).at(i_trk).at(i_cent) = new TH1D(name.c_str(), name.c_str(), N_dR, array_dr_bins);

				name = Form("h_TM_rN_r_1D_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_TM_rN_r_1D.at(i_jet).at(i_trk).at(i_cent) = new TH1D(name.c_str(), name.c_str(), N_dR, array_dr_bins);

				name = Form("h_TM_rN_TM_r_1D_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_TM_rN_TM_r_1D.at(i_jet).at(i_trk).at(i_cent) = new TH1D(name.c_str(), name.c_str(), N_dR, array_dr_bins);
			}
		}
	}


	for (int i_run = run_start; i_run <= run_end; i_run++)
	{
		cout << Form("Run %i: %1.0f", i_run, run_binning->GetBinLowEdge(i_run)) << endl;
		c_run_dep_jets->cd();
		c_run_dep_jets->Clear();
		c_run_dep_jets->Divide(3,2);

		int run_itr = 0;
		for (int i_cent = 0; i_cent < 6; i_cent++)
		{
			string cent_label = num_to_cent(31,i_cent).c_str();

			name = Form("MB_norm_jet_cent%i", i_cent);
			TH1* TM_norm_jet = (TH1*)((TH1*)input_file->Get(name.c_str()))->Clone(Form("TM_norm_jet_c%i", i_cent));
			TM_norm_jet->SetName(Form("%s_mc",name.c_str()));
			TM_norm_jet->Sumw2();

			name = Form("TM_norm_jet_rN_cent%i", i_cent);
			TH2* TM_rN_norm_jet_2D = (TH2*)((TH2*)input_file->Get(name.c_str()))->Clone(Form("TM_rN_norm_jet_c%i", i_cent));
			TM_rN_norm_jet_2D->SetName(Form("%s_mc",name.c_str()));
			TM_rN_norm_jet_2D->Sumw2();
			TH1* TM_rN_norm_jet = (TH1*)TM_rN_norm_jet_2D->ProjectionY(Form("%s_mc_rN%i",name.c_str(), i_run), i_run,i_run);


			if (i_cent == 0 && run_itr == 0)
			{
				legend_jet_run_dep->AddEntry(TM_norm_jet, "Combined", "lp");
				legend_jet_run_dep->AddEntry(TM_rN_norm_jet, Form("run %i: %1.0f", i_run, run_binning->GetBinLowEdge(i_run)), "lp");
			}
			//drawing jets
			c_run_dep_jets->cd(i_cent+1);
			SetHStyle_smallify(TM_norm_jet,0,1);
			SetHStyle_smallify(TM_rN_norm_jet,1,1);
			TM_rN_norm_jet->GetYaxis()->SetRangeUser(1E-9,1);
			TM_rN_norm_jet->DrawCopy("");
			TM_norm_jet->DrawCopy("same");
			gPad->SetLogx();
			gPad->SetLogy();

//			if (i_run == run_start)
			{
				c_run_dep_jets->cd(i_cent+1);
				ltx->DrawLatexNDC(0.92,0.90,cent_label.c_str());
				legend_jet_run_dep->Draw();

			}
			for (int i_dR = 0; i_dR < N_dR; i_dR++)
			{
				name = Form("ChPS_TM_UE_dR%i_cent%i", i_dR, i_cent);
				h_TM_2D[i_dR][i_cent] = (TH2*)input_file->Get(name.c_str())->Clone(Form("%s_run%i",name.c_str(), i_run));

				name = Form("ChPS_TM_UE_rN_dR%i_cent%i", i_dR, i_cent);
				TH3* h_TM_rN_3D = (TH3*)input_file->Get(name.c_str())->Clone(Form("%s_run%i",name.c_str(), i_run));;
				h_TM_rN_3D->GetXaxis()->SetRange(i_run,i_run);
				h_TM_rN_2D[i_dR][i_cent] = (TH2*)h_TM_rN_3D->Project3D("zy");
				h_TM_rN_2D[i_dR][i_cent]->SetName(Form("%s_run%i",name.c_str(), i_run));

				for (int i_jet = jet_pt_start; i_jet <= jet_pt_end; i_jet++)
				{
					double n_jets_tm = TM_norm_jet->GetBinContent(i_jet);
					double n_jets_tm_rN = TM_rN_norm_jet->GetBinContent(i_jet);

					for (int i_trk = trk_pt_start; i_trk <= trk_pt_end; i_trk++)
					{
						if (n_jets_tm != 0)
						{
							double updated_UE_TM = h_TM_2D[i_dR][i_cent]->GetBinContent(i_trk, i_jet) / n_jets_tm;
							double updated_UE_TM_err = h_TM_2D[i_dR][i_cent]->GetBinError(i_trk, i_jet) / n_jets_tm;
							h_TM_2D[i_dR][i_cent]->SetBinContent(i_trk, i_jet, updated_UE_TM);
							h_TM_2D[i_dR][i_cent]->SetBinError(i_trk, i_jet, updated_UE_TM_err);
						}
						else h_TM_2D[i_dR][i_cent]->Reset();

						if (n_jets_tm_rN != 0)
						{
							double updated_UE_TM_rN = h_TM_rN_2D[i_dR][i_cent]->GetBinContent(i_trk, i_jet) / n_jets_tm_rN;
							double updated_UE_TM_rN_err = h_TM_rN_2D[i_dR][i_cent]->GetBinError(i_trk, i_jet) / n_jets_tm_rN;
							h_TM_rN_2D[i_dR][i_cent]->SetBinContent(i_trk, i_jet, updated_UE_TM_rN);
							h_TM_rN_2D[i_dR][i_cent]->SetBinError(i_trk, i_jet, updated_UE_TM_rN_err);
						}
						else h_TM_rN_2D[i_dR][i_cent]->Reset();
					}
				}


				//has to be integrated over jet cone before UE-JER correction. the correction will then be based on the integrated TM/Cone_MC
				if (dR_binning->GetBinLowEdge(i_dR+1) < 0.4)
				{

					if (i_dR == 0)
					{
						h_TM_injet_2D[i_cent] = (TH2*)h_TM_2D[i_dR][i_cent]->Clone(Form("%s_injet",h_TM_2D[i_dR][i_cent]->GetName()));
						h_TM_rN_injet_2D[i_cent] = (TH2*)h_TM_rN_2D[i_dR][i_cent]->Clone(Form("%s_injet",h_TM_rN_2D[i_dR][i_cent]->GetName()));
					}
					else
					{
						h_TM_injet_2D[i_cent]->Add(h_TM_2D[i_dR][i_cent]);
						h_TM_rN_injet_2D[i_cent]->Add(h_TM_rN_2D[i_dR][i_cent]);
					}
				}

				//compare to TM
				if (subtract)
				{
					name = Form("sub_TM_rN_TM_dR%i_cent%i", i_dR, i_cent);
					h_TM_rN_TM_2D[i_dR][i_cent] = (TH2*)h_TM_rN_2D[i_dR][i_cent]->Clone(name.c_str());
					h_TM_rN_TM_2D[i_dR][i_cent]->Add(h_TM_2D[i_dR][i_cent],-1);
				}

				else
				{
					name = Form("ratio_TM_rN_TM_dR%i_cent%i", i_dR, i_cent);
					h_TM_rN_TM_2D[i_dR][i_cent] = (TH2*)h_TM_rN_2D[i_dR][i_cent]->Clone(name.c_str());
					h_TM_rN_TM_2D[i_dR][i_cent]->Divide(h_TM_2D[i_dR][i_cent]);
				}


				double dR_lo = dR_binning->GetBinLowEdge(i_dR+1);
				double dR_hi = dR_binning->GetBinUpEdge(i_dR+1);
				double area = TMath::Pi() * ((dR_hi*dR_hi) - (dR_lo*dR_lo));

				for (int i_jet = jet_pt_start-1; i_jet < jet_pt_end; i_jet++)
				{
					h_TM_1D[i_jet][i_dR][i_cent] = (TH1*)h_TM_2D[i_dR][i_cent]->ProjectionX(Form("TM_jet%i_dr%i_cent%i", i_jet, i_dR, i_cent), i_jet+1, i_jet+1);
					h_TM_rN_1D[i_jet][i_dR][i_cent] = (TH1*)h_TM_rN_2D[i_dR][i_cent]->ProjectionX(Form("TM_rN_jet%i_dr%i_cent%i", i_jet, i_dR, i_cent), i_jet+1, i_jet+1);

					h_TM_1D[i_jet][i_dR][i_cent]->Scale(1.,"width");
					h_TM_rN_1D[i_jet][i_dR][i_cent]->Scale(1.,"width");

					h_TM_1D[i_jet][i_dR][i_cent]->Scale(1./area);
					h_TM_rN_1D[i_jet][i_dR][i_cent]->Scale(1./area);

					//dont scale ratios if dividing
					h_TM_rN_TM_1D[i_jet][i_dR][i_cent] = (TH1*)h_TM_rN_TM_2D[i_dR][i_cent]->ProjectionX(Form("TM_rN_TM_jet%i_dr%i_cent%i", i_jet, i_dR, i_cent), i_jet+1, i_jet+1);

					if (subtract)
					{
						h_TM_rN_TM_1D[i_jet][i_dR][i_cent]->Scale(1.,"width");
						h_TM_rN_TM_1D[i_jet][i_dR][i_cent]->Scale(1./area);
					}

				}


//				output_file->cd();
//				name = Form("ChPS_TM_UE_dR%i_cent%i_run%i", i_dR, i_cent,i_run);
//				h_TM_2D[i_dR][i_cent]->Write(name.c_str());
//
				name = Form("ChPS_TM_rN_UE_dR%i_cent%i_run%i", i_dR, i_cent,i_run);
				h_TM_rN_2D[i_dR][i_cent]->Write(name.c_str());


				delete h_TM_rN_3D;
			}

			//recast in terms of r
			for (int i_jet = jet_pt_start-1; i_jet < jet_pt_end; i_jet++)
			{
				for (int i_trk = trk_pt_start-1; i_trk < trk_pt_end; i_trk++)
				{
					for (int i_dR = 0; i_dR < N_dR; i_dR++)
					{
						h_TM_r_1D[i_jet][i_trk][i_cent]->SetBinContent(i_dR+1, h_TM_1D[i_jet][i_dR][i_cent]->GetBinContent(i_trk+1));
						h_TM_rN_r_1D[i_jet][i_trk][i_cent]->SetBinContent(i_dR+1, h_TM_rN_1D[i_jet][i_dR][i_cent]->GetBinContent(i_trk+1));

						h_TM_r_1D[i_jet][i_trk][i_cent]->SetBinError(i_dR+1, h_TM_1D[i_jet][i_dR][i_cent]->GetBinError(i_trk+1));
						h_TM_rN_r_1D[i_jet][i_trk][i_cent]->SetBinError(i_dR+1, h_TM_rN_1D[i_jet][i_dR][i_cent]->GetBinError(i_trk+1));

						h_TM_r_1D[i_jet][i_trk][i_cent]->GetXaxis()->SetTitle("r");
						h_TM_rN_r_1D[i_jet][i_trk][i_cent]->GetXaxis()->SetTitle("r");

						h_TM_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetTitle("D_{p_{T}} [UE]");
						h_TM_rN_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetTitle("D_{p_{T}} [UE]");


						//ratios
						h_TM_rN_TM_r_1D[i_jet][i_trk][i_cent]->SetBinContent(i_dR+1, h_TM_rN_TM_1D[i_jet][i_dR][i_cent]->GetBinContent(i_trk+1));
						h_TM_rN_TM_r_1D[i_jet][i_trk][i_cent]->SetBinError(i_dR+1, h_TM_rN_TM_1D[i_jet][i_dR][i_cent]->GetBinError(i_trk+1));
						h_TM_rN_TM_r_1D[i_jet][i_trk][i_cent]->GetXaxis()->SetTitle("r");
						h_TM_rN_TM_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetTitle("Ratio");

					}
				}
			}


			output_file->cd();
			name = Form("TM_norm_jet_cent%i_run%i", i_cent, i_run);
			TM_norm_jet->Write(name.c_str());

			name = Form("TM_rN_norm_jet_cent%i_run%i", i_cent,i_run);
			TM_rN_norm_jet->Write(name.c_str());
//
//			for (int i_dR = 0; i_dR < N_dR; i_dR++)
//			{
//				for (int i_jet = jet_pt_start-1; i_jet < jet_pt_end; i_jet++)
//				{
//					name = Form("h_ChPS_UE_TM_dR%i_cent%i_jetpt%i_run%i", i_dR, i_cent, i_jet, i_run);
//					h_TM_1D[i_jet][i_dR][i_cent]->Write(name.c_str());
//
//					name = Form("h_ChPS_UE_TM_rN_dR%i_cent%i_jetpt%i_run%i", i_dR, i_cent, i_jet, i_run);
//					h_TM_rN_1D[i_jet][i_dR][i_cent]->Write(name.c_str());
//				}
//			}
//
//			h_TM_injet_2D[i_cent]->Write(Form("MC_Cent%i_run%i",i_cent,i_run));
//			h_TM_rN_injet_2D[i_cent]->Write(Form("TM_rN_Cent%i_run%i",i_cent,i_run));

			for (int i_jet = jet_pt_start-1; i_jet < jet_pt_end; i_jet++)
			{
				for (int i_trk = trk_pt_start-1; i_trk < trk_pt_end; i_trk++)
				{
					output_file->cd();

					name = Form("jet%i_trk%i_cent%i", i_jet, i_trk, i_cent);
					h_TM_r_1D[i_jet][i_trk][i_cent]->SetName(Form("UE_TM_indR_%s_run%i", name.c_str(), i_run));
					h_TM_rN_r_1D[i_jet][i_trk][i_cent]->SetName(Form("UE_TM_rN_indR_%s_run%i", name.c_str(), i_run));
					h_TM_rN_TM_r_1D[i_jet][i_trk][i_cent]->SetName(Form("UE_TM_rN_TM_indR_%s_run%i", name.c_str(), i_run));

					h_TM_r_1D[i_jet][i_trk][i_cent]->Write(Form("UE_TM_indR_%s_run%i", name.c_str(), i_run));
					h_TM_rN_r_1D[i_jet][i_trk][i_cent]->Write(Form("UE_TM_rN_indR_%s_run%i", name.c_str(), i_run));
//					h_TM_rN_TM_r_1D[i_jet][i_trk][i_cent]->Write(Form("UE_TM_rN_TM_indR_%s_run%i", name.c_str(), i_run));
				}
			}

			delete TM_norm_jet;
			delete TM_rN_norm_jet_2D;

		}


//		{
//			//as function of r
////			cout << "Function of R" << endl;
//			TCanvas *c_x = new TCanvas("c_x","c_x",900,600);
//
//			TLegend *legend_x = new TLegend(0.20, 0.45, 0.40, 0.75, "","brNDC");
//			legend_x->SetTextFont(43);
//			legend_x->SetBorderSize(0);
//			legend_x->SetTextSize(10);
//			legend_x->SetNColumns(1);
//
//			TLegend *legend_y = new TLegend(0.50, 0.4, 0.65, 0.450, "","brNDC");
//			legend_y->SetTextFont(43);
//			legend_y->SetBorderSize(0);
//			legend_y->SetTextSize(10);
//			legend_y->SetNColumns(1);
//
//
//			for (int i_trk = trk_pt_start-1; i_trk < trk_pt_end; i_trk++)
//			{
//				string trk_label = Form("%1.1f < p_{T}^{Trk} < %1.1f GeV", trkpT_binning->GetBinLowEdge(i_trk+1), trkpT_binning->GetBinUpEdge(i_trk+1));
//
//				int jet_itr = 0;
//				for (int i_jet = jet_pt_start-1; i_jet < jet_pt_end; i_jet++)
//				{
//					string jet_label = Form("%1.0f < p_{T}^{Jet} < %1.0f GeV", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));
//					c_x->Clear();
//					c_x->Divide(3,2);
//
//					for (int i_cent = 0; i_cent < 6; i_cent++)
//					{
//						SetHStyle_smallify(h_TM_r_1D[i_jet][i_trk][i_cent], 0, 1);
//						SetHStyle_smallify(h_TM_rN_r_1D[i_jet][i_trk][i_cent], 1, 1);
//
//						SetHStyle_smallify(h_TM_rN_TM_r_1D[i_jet][i_trk][i_cent], 1, 1);
//
//						if (jet_itr == 0 && i_trk == 2 && i_cent == 0)
//						{
//							legend_x->AddEntry(h_TM_r_1D[i_jet][i_trk][i_cent],"TM","lp");
//							legend_x->AddEntry(h_TM_rN_r_1D[i_jet][i_trk][i_cent],"TM_rN","lp");
//						}
//
//						double low_range, hi_range;
//						low_range = h_TM_r_1D[i_jet][i_trk][i_cent]->GetMinimum() * 0.75;
//						hi_range = h_TM_r_1D[i_jet][i_trk][i_cent]->GetMaximum() * 1.4;
//
//						low_range = h_TM_r_1D[i_jet][i_trk][i_cent]->GetBinContent(8) * 0.75;
//						hi_range = h_TM_r_1D[i_jet][i_trk][i_cent]->GetBinContent(1) * 1.4;
//						h_TM_r_1D[i_jet][i_trk][i_cent]->GetXaxis()->SetRangeUser(0, r_max_range);
//						h_TM_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetRangeUser(low_range, hi_range);
//
//						if (i_cent == 0) {low_range = 60; hi_range = 100;}
//						if (i_cent == 1) {low_range = 40; hi_range = 60;}
//						if (i_cent == 2) {low_range = 20; hi_range = 40;}
//						if (i_cent == 3) {low_range = 15; hi_range = 25;}
//						if (i_cent == 4) {low_range = 5; hi_range = 15;}
//						if (i_cent == 5) {low_range = 1.5; hi_range = 4;}
////						h_TM_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetRangeUser(low_range, hi_range);
//
//						if (subtract) h_TM_rN_TM_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetRangeUser(-0.5,0.5);
//						h_TM_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetNdivisions(504);
//						h_TM_rN_TM_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetNdivisions(505);
//						h_TM_rN_TM_r_1D[i_jet][i_trk][i_cent]->GetXaxis()->SetRangeUser(0, r_max_range);
//						h_TM_rN_TM_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetRangeUser(0.9,1.1);
//
//						c_x->cd(i_cent+1);
//						gPad->Divide(1,2);
//
//						c_x->cd(i_cent+1)->cd(1);
//						gPad->SetPad(0,0.40,0.99,0.99);
//						gPad->SetTopMargin(0.05);
//						gPad->SetBottomMargin(0.0);
//						gPad->SetRightMargin(0);
//
//						h_TM_r_1D[i_jet][i_trk][i_cent]->Draw("");
//						h_TM_rN_r_1D[i_jet][i_trk][i_cent]->Draw("same");
//						gPad->SetLogy(0);
//
//						c_x->cd(i_cent+1)->cd(2);
//						gPad->SetPad(0,0.0,0.99,0.40);
//						gPad->SetTopMargin(0.0);
//						gPad->SetBottomMargin(0.30);
//						gPad->SetRightMargin(0);
//
//
//
//						if (subtract) h_TM_rN_TM_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetTitle("x - TM");
//						else  h_TM_rN_TM_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetTitle("x / MC");
//						h_TM_rN_TM_r_1D[i_jet][i_trk][i_cent]->Draw("");
//						if (subtract) line->DrawLine(0, 0, r_max_range, 0);
//						else line->DrawLine(0, 1, r_max_range, 1);
//
//						c_x->cd(i_cent+1);
//						ltx->SetTextAlign(32);
//						ltx->DrawLatexNDC(0.94,0.94,num_to_cent(31,i_cent).c_str());
//						ltx->DrawLatexNDC(0.94,0.88,jet_label.c_str());
//						ltx->DrawLatexNDC(0.94,0.82,trk_label.c_str());
//
//					}
//					c_x->cd(1)->cd(1);
//					legend_x->Draw();
//					c_x->cd(1)->cd(2);
//					legend_y->Draw();
//
//					jet_itr++;
//
//					if (i_trk == trk_pt_start-1 && i_jet == jet_pt_start-1) name = "(";
//					else if (i_trk == trk_pt_end-1 && i_jet == jet_pt_end - 1) name = ")";
//					else name = "";
//					if (subtract) c_x->Print(Form("run_dep/UE_x_subtr_run%i_c%i.pdf%s",i_run, sys_mode, name.c_str()), Form("Title: trk%i_jet%i", i_trk, i_jet));
//					else c_x->Print(Form("run_dep/UE_x_ratio_run%i_c%i.pdf%s", i_run, sys_mode, name.c_str()), Form("Title: trk%i_jet%i", i_trk, i_jet));
//
//				}
//
//			}
//
//			delete c_x;
//			delete legend_x;
//			delete legend_y;
//		}



		name = "";
		if (i_run == run_start) name = "(";
		else if (i_run == run_end) name = ")";
		c_run_dep_jets->Print(Form("run_dep/jet_pt_run_c%i.pdf%s", sys_mode, name.c_str()));
		legend_jet_run_dep->Clear();

	}








}
