#include "output_dev/functions/global_variables.h"

void run_dep(int sys_mode = 38)
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
	TFile *input_file_data = new TFile(Form("output_dev/raw_results/%s/FF_data_out_histo_PbPb_5p02_r001.root", sys_path.c_str()));

	TH1* h_tmp_mc = (TH1*)input_file->Get("h_event_rN");
	h_tmp_mc->SetName("h_event_rN_mc");
	//	h_tmp->Scale(1./h_tmp->Integral());

	TH1* h_tmp_data = (TH1*)input_file_data->Get("h_event_rN");
	h_tmp_data->SetName("h_event_rN_data");
//	h_tmp_data->Scale(1./h_tmp_data->Integral());
	double j = 0;

	for (int i = 1; i <= h_tmp_data->GetXaxis()->GetNbins(); i++)
	{
		j = j + h_tmp_data->GetBinContent(i);
	}
	cout << j <<  "<--------------" << endl;

	name = Form("./run_dep/run_dep_UE_%s.root", sys_path.c_str());
	TFile *output_file = new TFile(name.c_str(), "recreate");

	cout << "Using file " << input_file->GetName() << endl;
	output_file->cd();
	h_tmp_mc->Write("EventPercentages_mc");
	h_tmp_data->Write("EventPercentages_data");


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


	//1D UE
	vector<vector<vector<TH1*>>> h_TM_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_dR, vector<TH1*> (n_cent_cuts)));
	vector<vector<vector<TH1*>>> h_TM_rN_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_dR, vector<TH1*> (n_cent_cuts)));
	vector<vector<vector<TH1*>>> h_MB_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_dR, vector<TH1*> (n_cent_cuts)));
	vector<vector<vector<TH1*>>> h_MB_rN_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_dR, vector<TH1*> (n_cent_cuts)));


	//1D UE in r
	vector<vector<vector<TH1*>>> h_TM_r_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_trkpt, vector<TH1*> (n_cent_cuts)));
	vector<vector<vector<TH1*>>> h_TM_rN_r_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_trkpt, vector<TH1*> (n_cent_cuts)));
	vector<vector<vector<TH1*>>> h_MB_r_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_trkpt, vector<TH1*> (n_cent_cuts)));
	vector<vector<vector<TH1*>>> h_MB_rN_r_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_trkpt, vector<TH1*> (n_cent_cuts)));

	//2D ratio
	vector<vector<TH2*>> h_TM_rN_TM_2D = vector<vector<TH2*>> (N_dR, vector<TH2*> (n_cent_cuts));

	//1D ratio
	vector<vector<vector<TH1*>>> h_TM_rN_TM_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_dR, vector<TH1*> (n_cent_cuts)));

	//1D ratio in r
	vector<vector<vector<TH1*>>> h_TM_rN_TM_r_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_trkpt, vector<TH1*> (n_cent_cuts)));



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

				name = Form("h_MB_r_1D_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_MB_r_1D.at(i_jet).at(i_trk).at(i_cent) = new TH1D(name.c_str(), name.c_str(), N_dR, array_dr_bins);

				name = Form("h_MB_rN_r_1D_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_MB_rN_r_1D.at(i_jet).at(i_trk).at(i_cent) = new TH1D(name.c_str(), name.c_str(), N_dR, array_dr_bins);

				name = Form("h_TM_rN_TM_r_1D_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_TM_rN_TM_r_1D.at(i_jet).at(i_trk).at(i_cent) = new TH1D(name.c_str(), name.c_str(), N_dR, array_dr_bins);
			}
		}
	}


	for (int i_run = run_start; i_run <= run_end; i_run++)
	{
		cout << Form("Run %i: %1.0f", i_run, run_binning->GetBinLowEdge(i_run)) << endl;

		int run_itr = 0;
		for (int i_cent = 0; i_cent < 6; i_cent++)
		{
			string cent_label = num_to_cent(31,i_cent).c_str();

			name = Form("TM_norm_jet_cent%i", i_cent);
			TH1* TM_norm_jet = (TH1*)((TH1*)input_file->Get(name.c_str()))->Clone(Form("TM_norm_jet_c%i", i_cent));
			TM_norm_jet->SetName(Form("%s_mc",name.c_str()));
			TM_norm_jet->Sumw2();

			name = Form("TM_norm_jet_rN_cent%i", i_cent);
			TH2* TM_rN_norm_jet_2D = (TH2*)((TH2*)input_file->Get(name.c_str()))->Clone(Form("TM_rN_norm_jet_c%i", i_cent));
			TM_rN_norm_jet_2D->SetName(Form("%s_mc",name.c_str()));
			TM_rN_norm_jet_2D->Sumw2();
			TH1* TM_rN_norm_jet = (TH1*)TM_rN_norm_jet_2D->ProjectionY(Form("%s_mc_rN%i",name.c_str(), i_run), i_run,i_run);

			name = Form("MB_norm_jet_cent%i", i_cent);
			TH1* MB_norm_jet = (TH1*)((TH1*)input_file_data->Get(name.c_str()))->Clone(Form("MB_norm_jet_c%i", i_cent));
			MB_norm_jet->SetName(Form("%s_data",name.c_str()));
			MB_norm_jet->Sumw2();

			name = Form("TM_norm_jet_rN_cent%i", i_cent); //TM = MB for jet spectra, done because I forgot to fill MB_jet_norm_rN on condor
			TH2* MB_rN_norm_jet_2D = (TH2*)((TH2*)input_file_data->Get(name.c_str()))->Clone(Form("MB_rN_norm_jet_c%i", i_cent));
			MB_rN_norm_jet_2D->SetName(Form("%s_data",name.c_str()));
			MB_rN_norm_jet_2D->Sumw2();
			TH1* MB_rN_norm_jet = (TH1*)MB_rN_norm_jet_2D->ProjectionY(Form("%s_data_rN%i",name.c_str(), i_run), i_run,i_run);

			for (int i_dR = 0; i_dR < N_dR; i_dR++)
			{
				name = Form("ChPS_TM_UE_dR%i_cent%i", i_dR, i_cent);
				TH2* h_TM_2D = (TH2*)input_file->Get(name.c_str())->Clone(Form("%s_run%i",name.c_str(), i_run));

				name = Form("ChPS_TM_UE_rN_dR%i_cent%i", i_dR, i_cent);
				TH3* h_TM_rN_3D = (TH3*)input_file->Get(name.c_str())->Clone(Form("%s_run%i",name.c_str(), i_run));;
				h_TM_rN_3D->GetXaxis()->SetRange(i_run,i_run);
				TH2* h_TM_rN_2D = (TH2*)h_TM_rN_3D->Project3D("zy");
				h_TM_rN_2D->SetName(Form("%s_run%i",name.c_str(), i_run));

				name = Form("ChPS_MB_UE_dR%i_cent%i", i_dR, i_cent);
				TH2* h_MB_2D = (TH2*)input_file_data->Get(name.c_str())->Clone(Form("%s_run%i",name.c_str(), i_run));

				name = Form("ChPS_MB_UE_rN_dR%i_cent%i", i_dR, i_cent);
				TH3* h_MB_rN_3D = (TH3*)input_file_data->Get(name.c_str())->Clone(Form("%s_run%i",name.c_str(), i_run));;
				h_MB_rN_3D->GetXaxis()->SetRange(i_run,i_run);
				TH2* h_MB_rN_2D = (TH2*)h_MB_rN_3D->Project3D("zy");
				h_MB_rN_2D->SetName(Form("%s_run%i",name.c_str(), i_run));

				for (int i_jet = jet_pt_start; i_jet <= jet_pt_end; i_jet++)
				{
					double n_jets_TM = TM_norm_jet->GetBinContent(i_jet);
					double n_jets_TM_rN = TM_rN_norm_jet->GetBinContent(i_jet);

					double n_jets_MB = MB_norm_jet->GetBinContent(i_jet);
					double n_jets_MB_rN = MB_rN_norm_jet->GetBinContent(i_jet);

					for (int i_trk = trk_pt_start; i_trk <= trk_pt_end; i_trk++)
					{
						if (n_jets_TM != 0)
						{
							double updated_UE_TM = h_TM_2D->GetBinContent(i_trk, i_jet) / n_jets_TM;
							double updated_UE_TM_err = h_TM_2D->GetBinError(i_trk, i_jet) / n_jets_TM;
							h_TM_2D->SetBinContent(i_trk, i_jet, updated_UE_TM);
							h_TM_2D->SetBinError(i_trk, i_jet, updated_UE_TM_err);
						}
						else
						{
							h_TM_2D->SetBinContent(i_trk, i_jet, 0);
							h_TM_2D->SetBinError(i_trk, i_jet, 0);
						}

						if (n_jets_TM_rN != 0)
						{
							double updated_UE_TM_rN = h_TM_rN_2D->GetBinContent(i_trk, i_jet) / n_jets_TM_rN;
							double updated_UE_TM_rN_err = h_TM_rN_2D->GetBinError(i_trk, i_jet) / n_jets_TM_rN;
							h_TM_rN_2D->SetBinContent(i_trk, i_jet, updated_UE_TM_rN);
							h_TM_rN_2D->SetBinError(i_trk, i_jet, updated_UE_TM_rN_err);
						}
						else
						{
							h_TM_rN_2D->SetBinContent(i_trk, i_jet, 0);
							h_TM_rN_2D->SetBinError(i_trk, i_jet, 0);
						}


						if (n_jets_MB != 0)
						{
							double updated_UE_MB = h_MB_2D->GetBinContent(i_trk, i_jet) / n_jets_MB;
							double updated_UE_MB_err = h_MB_2D->GetBinError(i_trk, i_jet) / n_jets_MB;
							h_MB_2D->SetBinContent(i_trk, i_jet, updated_UE_MB);
							h_MB_2D->SetBinError(i_trk, i_jet, updated_UE_MB_err);
						}
						else
						{
							h_MB_2D->SetBinContent(i_trk, i_jet, 0);
							h_MB_2D->SetBinError(i_trk, i_jet, 0);
						}

						if (n_jets_MB_rN != 0)
						{
							double updated_UE_MB_rN = h_MB_rN_2D->GetBinContent(i_trk, i_jet) / n_jets_MB_rN;
							double updated_UE_MB_rN_err = h_MB_rN_2D->GetBinError(i_trk, i_jet) / n_jets_MB_rN;
							h_MB_rN_2D->SetBinContent(i_trk, i_jet, updated_UE_MB_rN);
							h_MB_rN_2D->SetBinError(i_trk, i_jet, updated_UE_MB_rN_err);
						}
						else
						{
							h_MB_rN_2D->SetBinContent(i_trk, i_jet, 0);
							h_MB_rN_2D->SetBinError(i_trk, i_jet, 0);
						}



					}
				}

				double dR_lo = dR_binning->GetBinLowEdge(i_dR+1);
				double dR_hi = dR_binning->GetBinUpEdge(i_dR+1);
				double area = TMath::Pi() * ((dR_hi*dR_hi) - (dR_lo*dR_lo));

				for (int i_jet = jet_pt_start-1; i_jet < jet_pt_end; i_jet++)
				{
					h_TM_1D[i_jet][i_dR][i_cent] = (TH1*)h_TM_2D->ProjectionX(Form("TM_jet%i_dr%i_cent%i", i_jet, i_dR, i_cent), i_jet+1, i_jet+1);
					h_TM_rN_1D[i_jet][i_dR][i_cent] = (TH1*)h_TM_rN_2D->ProjectionX(Form("TM_rN_jet%i_dr%i_cent%i", i_jet, i_dR, i_cent), i_jet+1, i_jet+1);

					h_TM_1D[i_jet][i_dR][i_cent]->Scale(1.,"width");
					h_TM_rN_1D[i_jet][i_dR][i_cent]->Scale(1.,"width");

					h_TM_1D[i_jet][i_dR][i_cent]->Scale(1./area);
					h_TM_rN_1D[i_jet][i_dR][i_cent]->Scale(1./area);


					h_MB_1D[i_jet][i_dR][i_cent] = (TH1*)h_MB_2D->ProjectionX(Form("MB_jet%i_dr%i_cent%i", i_jet, i_dR, i_cent), i_jet+1, i_jet+1);
					h_MB_rN_1D[i_jet][i_dR][i_cent] = (TH1*)h_MB_rN_2D->ProjectionX(Form("MB_rN_jet%i_dr%i_cent%i", i_jet, i_dR, i_cent), i_jet+1, i_jet+1);

					h_MB_1D[i_jet][i_dR][i_cent]->Scale(1.,"width");
					h_MB_rN_1D[i_jet][i_dR][i_cent]->Scale(1.,"width");

					h_MB_1D[i_jet][i_dR][i_cent]->Scale(1./area);
					h_MB_rN_1D[i_jet][i_dR][i_cent]->Scale(1./area);

				}

				delete h_TM_rN_3D;
				delete h_MB_rN_3D;
				delete h_TM_2D;
				delete h_MB_2D;
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

						h_TM_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetTitle("D_{p_{T}} [UE_{TM}]");
						h_TM_rN_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetTitle("D_{p_{T}} [UE_{TM}^{Run}]");

						h_MB_r_1D[i_jet][i_trk][i_cent]->SetBinContent(i_dR+1, h_MB_1D[i_jet][i_dR][i_cent]->GetBinContent(i_trk+1));
						h_MB_rN_r_1D[i_jet][i_trk][i_cent]->SetBinContent(i_dR+1, h_MB_rN_1D[i_jet][i_dR][i_cent]->GetBinContent(i_trk+1));

						h_MB_r_1D[i_jet][i_trk][i_cent]->SetBinError(i_dR+1, h_MB_1D[i_jet][i_dR][i_cent]->GetBinError(i_trk+1));
						h_MB_rN_r_1D[i_jet][i_trk][i_cent]->SetBinError(i_dR+1, h_MB_rN_1D[i_jet][i_dR][i_cent]->GetBinError(i_trk+1));

						h_MB_r_1D[i_jet][i_trk][i_cent]->GetXaxis()->SetTitle("r");
						h_MB_rN_r_1D[i_jet][i_trk][i_cent]->GetXaxis()->SetTitle("r");

						h_MB_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetTitle("D_{p_{T}} [UE_{MB}]");
						h_MB_rN_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetTitle("D_{p_{T}} [UE_{MB}^{Run}]");

					}
				}
			}


			output_file->cd();
			name = Form("TM_norm_jet_cent%i_run%i", i_cent, i_run);
			TM_norm_jet->SetName(name.c_str());
			TM_norm_jet->SetTitle(name.c_str());
			TM_norm_jet->Write(name.c_str());

			name = Form("TM_rN_norm_jet_cent%i_run%i", i_cent,i_run);
			TM_rN_norm_jet->SetName(name.c_str());
			TM_rN_norm_jet->SetTitle(name.c_str());
			TM_rN_norm_jet->Write(name.c_str());

			name = Form("MB_norm_jet_cent%i_run%i", i_cent, i_run);
			MB_norm_jet->SetName(name.c_str());
			MB_norm_jet->SetTitle(name.c_str());
			MB_norm_jet->Write(name.c_str());

			name = Form("MB_rN_norm_jet_cent%i_run%i", i_cent,i_run);
			MB_rN_norm_jet->SetName(name.c_str());
			MB_rN_norm_jet->SetTitle(name.c_str());
			MB_rN_norm_jet->Write(name.c_str());

			for (int i_jet = jet_pt_start-1; i_jet < jet_pt_end; i_jet++)
			{
				for (int i_trk = trk_pt_start-1; i_trk < trk_pt_end; i_trk++)
				{
					output_file->cd();

					name = Form("jet%i_trk%i_cent%i", i_jet, i_trk, i_cent);
					h_TM_r_1D[i_jet][i_trk][i_cent]->SetName(Form("UE_TM_indR_%s_run%i", name.c_str(), i_run));
					h_TM_rN_r_1D[i_jet][i_trk][i_cent]->SetName(Form("UE_TM_rN_indR_%s_run%i", name.c_str(), i_run));
					h_TM_r_1D[i_jet][i_trk][i_cent]->Write(Form("UE_TM_indR_%s_run%i", name.c_str(), i_run));
					h_TM_rN_r_1D[i_jet][i_trk][i_cent]->Write(Form("UE_TM_rN_indR_%s_run%i", name.c_str(), i_run));

					h_MB_r_1D[i_jet][i_trk][i_cent]->SetName(Form("UE_MB_indR_%s_run%i", name.c_str(), i_run));
					h_MB_rN_r_1D[i_jet][i_trk][i_cent]->SetName(Form("UE_MB_rN_indR_%s_run%i", name.c_str(), i_run));
					h_MB_r_1D[i_jet][i_trk][i_cent]->Write(Form("UE_MB_indR_%s_run%i", name.c_str(), i_run));
					h_MB_rN_r_1D[i_jet][i_trk][i_cent]->Write(Form("UE_MB_rN_indR_%s_run%i", name.c_str(), i_run));
				}
			}


			delete TM_norm_jet;
			delete TM_rN_norm_jet_2D;
			delete TM_rN_norm_jet;
			delete MB_norm_jet;
			delete MB_rN_norm_jet_2D;
			delete MB_rN_norm_jet;


		}

	}

}
