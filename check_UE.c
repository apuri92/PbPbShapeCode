#include "output_dev/functions/global_variables.h"

void check_UE(int sys_mode = 47, bool subtract = 0)
{
	SetAtlasStyle();
	gErrorIgnoreLevel = 3001;
	string name;

	string sys_path;
	double r_max_range = 0.8;
	if (sys_mode == 0) sys_path = Form("nominal");
	if (sys_mode > 0 && sys_mode < 100) sys_path = Form("c%i", sys_mode);
	if (sys_mode > 100) sys_path = Form("sys%i", sys_mode);
//	TFile *input_file = new TFile(Form("./hist-local_mc.root"));
//	TFile *input_file_data = new TFile(Form("./hist-local_mc.root"));


//	TFile *input_file = new TFile(Form("hist-local_mc.root", sys_mode));
//	TFile *input_file_data = new TFile(Form("hist-local_mc.root", sys_mode));
	TFile *input_file = new TFile(Form("output_dev/raw_results/%s/FF_MC_out_histo_PbPb_5p02_r001.root", sys_path.c_str()));
	TFile *input_file_data = new TFile(Form("output_dev/raw_results/%s/FF_data_out_histo_PbPb_5p02_r001.root", sys_path.c_str()));
	name = Form("./UE_%s.root", sys_path.c_str());
	if (subtract) name = Form("./UE_subtr_%s.root", sys_path.c_str());

	TFile *output_file = new TFile(name.c_str(), "recreate");

	cout << "Using file " << input_file->GetName() << endl;
	cout << "Using file " << input_file_data->GetName() << endl;

	TAxis* dR_binning = (TAxis*)((TH3*)input_file->Get("h_dR_change_jetpt0_cent0"))->GetXaxis();
	TAxis* jetpT_binning = (TAxis*)((TH3*)input_file->Get("ChPS_raw_0_dR0_cent0"))->GetYaxis();
	TAxis* trkpT_binning = (TAxis*)((TH3*)input_file->Get("ChPS_raw_0_dR0_cent0"))->GetXaxis();
	TFile *UE_factors = new TFile(Form("./output_dev/unfold/output_pdf_c%i/root/UE_factors.root", sys_mode));
	bool apply_correctionFactors = false;

	output_file->cd();
	dR_binning->Write("dR_binning");
	jetpT_binning->Write("jetpT_binning");
	trkpT_binning->Write("trkpT_binning");

	TH1* h_events = (TH1*)input_file_data->Get("Centrality");
	output_file->cd();
	h_events->Write("Centrality");

	int jet_pt_start = 6;
	int jet_pt_end = 11;

	int N_dR = dR_binning->GetNbins();
	int N_jetpt = jetpT_binning->GetNbins();
	int N_trkpt = trkpT_binning->GetNbins();

	//2D UE
	vector<vector<TH2*>> h_cone_2D = vector<vector<TH2*>> (N_dR, vector<TH2*> (n_cent_cuts));
	vector<vector<TH2*>> h_cone_data_2D = vector<vector<TH2*>> (N_dR, vector<TH2*> (n_cent_cuts));
	vector<vector<TH2*>> h_MB_2D = vector<vector<TH2*>> (N_dR, vector<TH2*> (n_cent_cuts));
	vector<vector<TH2*>> h_MB_data_2D = vector<vector<TH2*>> (N_dR, vector<TH2*> (n_cent_cuts));
	vector<vector<TH2*>> h_TM_2D = vector<vector<TH2*>> (N_dR, vector<TH2*> (n_cent_cuts));
	vector<vector<TH2*>> h_FS_2D = vector<vector<TH2*>> (N_dR, vector<TH2*> (n_cent_cuts));
	vector<vector<TH2*>> h_FNS_2D = vector<vector<TH2*>> (N_dR, vector<TH2*> (n_cent_cuts));

	//1D UE
	vector<vector<vector<TH1*>>> h_cone_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_dR, vector<TH1*> (n_cent_cuts)));
	vector<vector<vector<TH1*>>> h_cone_data_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_dR, vector<TH1*> (n_cent_cuts)));
	vector<vector<vector<TH1*>>> h_MB_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_dR, vector<TH1*> (n_cent_cuts)));
	vector<vector<vector<TH1*>>> h_MB_data_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_dR, vector<TH1*> (n_cent_cuts)));
	vector<vector<vector<TH1*>>> h_TM_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_dR, vector<TH1*> (n_cent_cuts)));
	vector<vector<vector<TH1*>>> h_FS_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_dR, vector<TH1*> (n_cent_cuts)));
	vector<vector<vector<TH1*>>> h_FNS_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_dR, vector<TH1*> (n_cent_cuts)));

	//2D UE injet
	vector<TH2*> h_cone_injet_2D =  vector<TH2*> (n_cent_cuts);
	vector<TH2*> h_cone_data_injet_2D = vector<TH2*> (n_cent_cuts);
	vector<TH2*> h_MB_injet_2D = vector<TH2*> (n_cent_cuts);
	vector<TH2*> h_MB_data_injet_2D = vector<TH2*> (n_cent_cuts);
	vector<TH2*> h_TM_injet_2D = vector<TH2*> (n_cent_cuts);


	//1D UE u nr
	vector<vector<vector<TH1*>>> h_cone_r_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_trkpt, vector<TH1*> (n_cent_cuts)));
	vector<vector<vector<TH1*>>> h_cone_data_r_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_trkpt, vector<TH1*> (n_cent_cuts)));
	vector<vector<vector<TH1*>>> h_MB_r_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_trkpt, vector<TH1*> (n_cent_cuts)));
	vector<vector<vector<TH1*>>> h_MB_data_r_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_trkpt, vector<TH1*> (n_cent_cuts)));
	vector<vector<vector<TH1*>>> h_TM_r_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_trkpt, vector<TH1*> (n_cent_cuts)));
	vector<vector<vector<TH1*>>> h_FS_r_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_trkpt, vector<TH1*> (n_cent_cuts)));
	vector<vector<vector<TH1*>>> h_FNS_r_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_trkpt, vector<TH1*> (n_cent_cuts)));

	//2D ratio
	vector<vector<TH2*>> h_cone_MB_2D = vector<vector<TH2*>> (N_dR, vector<TH2*> (n_cent_cuts));
	vector<vector<TH2*>> h_cone_data_MB_2D = vector<vector<TH2*>> (N_dR, vector<TH2*> (n_cent_cuts));
	vector<vector<TH2*>> h_MB_data_MB_2D = vector<vector<TH2*>> (N_dR, vector<TH2*> (n_cent_cuts));
	vector<vector<TH2*>> h_cone_data_MB_data_2D = vector<vector<TH2*>> (N_dR, vector<TH2*> (n_cent_cuts));
	vector<vector<TH2*>> h_TM_MB_2D = vector<vector<TH2*>> (N_dR, vector<TH2*> (n_cent_cuts));
	vector<vector<TH2*>> h_FS_MB_2D = vector<vector<TH2*>> (N_dR, vector<TH2*> (n_cent_cuts));
	vector<vector<TH2*>> h_FNS_MB_2D = vector<vector<TH2*>> (N_dR, vector<TH2*> (n_cent_cuts));

	//1D ratio
	vector<vector<vector<TH1*>>> h_cone_MB_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_dR, vector<TH1*> (n_cent_cuts)));
	vector<vector<vector<TH1*>>> h_cone_data_MB_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_dR, vector<TH1*> (n_cent_cuts)));
	vector<vector<vector<TH1*>>> h_MB_data_MB_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_dR, vector<TH1*> (n_cent_cuts)));
	vector<vector<vector<TH1*>>> h_cone_data_MB_data_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_dR, vector<TH1*> (n_cent_cuts)));
	vector<vector<vector<TH1*>>> h_TM_MB_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_dR, vector<TH1*> (n_cent_cuts)));
	vector<vector<vector<TH1*>>> h_FS_MB_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_dR, vector<TH1*> (n_cent_cuts)));
	vector<vector<vector<TH1*>>> h_FNS_MB_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_dR, vector<TH1*> (n_cent_cuts)));

	//1D ratio in r
	vector<vector<vector<TH1*>>> h_cone_MB_r_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_trkpt, vector<TH1*> (n_cent_cuts)));
	vector<vector<vector<TH1*>>> h_cone_data_MB_r_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_trkpt, vector<TH1*> (n_cent_cuts)));
	vector<vector<vector<TH1*>>> h_TM_MB_r_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_trkpt, vector<TH1*> (n_cent_cuts)));
	vector<vector<vector<TH1*>>> h_MB_data_MB_r_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_trkpt, vector<TH1*> (n_cent_cuts)));
	vector<vector<vector<TH1*>>> h_cone_data_MB_data_r_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_trkpt, vector<TH1*> (n_cent_cuts)));
	vector<vector<vector<TH1*>>> h_FS_MB_r_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_trkpt, vector<TH1*> (n_cent_cuts)));
	vector<vector<vector<TH1*>>> h_FNS_MB_r_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_trkpt, vector<TH1*> (n_cent_cuts)));

	TLine *line = new TLine();
	line->SetLineColor(kBlack);



	double array_dr_bins[N_dR+1];
	for (int i_dR = 0; i_dR <= N_dR; i_dR++) array_dr_bins[i_dR] = dR_binning->GetBinLowEdge(i_dR+1);

	for (int i_cent = 0; i_cent < n_cent_cuts; i_cent++)
	{
		for (int i_jet = 0; i_jet < N_jetpt; i_jet++)
		{
			for (int i_trk = 0; i_trk < N_trkpt; i_trk++)
			{
				name = Form("h_MB_r_1D_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_MB_r_1D.at(i_jet).at(i_trk).at(i_cent) = new TH1D(name.c_str(), name.c_str(), N_dR, array_dr_bins);

				name = Form("h_MB_data_r_1D_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_MB_data_r_1D.at(i_jet).at(i_trk).at(i_cent) = new TH1D(name.c_str(), name.c_str(), N_dR, array_dr_bins);

				name = Form("h_cone_r_1D_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_cone_r_1D.at(i_jet).at(i_trk).at(i_cent) = new TH1D(name.c_str(), name.c_str(), N_dR, array_dr_bins);

				name = Form("h_cone_data_r_1D_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_cone_data_r_1D.at(i_jet).at(i_trk).at(i_cent) = new TH1D(name.c_str(), name.c_str(), N_dR, array_dr_bins);

				name = Form("h_TM_r_1D_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_TM_r_1D.at(i_jet).at(i_trk).at(i_cent) = new TH1D(name.c_str(), name.c_str(), N_dR, array_dr_bins);

				name = Form("h_FS_r_1D_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_FS_r_1D.at(i_jet).at(i_trk).at(i_cent) = new TH1D(name.c_str(), name.c_str(), N_dR, array_dr_bins);

				name = Form("h_FNS_r_1D_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_FNS_r_1D.at(i_jet).at(i_trk).at(i_cent) = new TH1D(name.c_str(), name.c_str(), N_dR, array_dr_bins);



				name = Form("h_MB_data_MB_r_1D_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_MB_data_MB_r_1D.at(i_jet).at(i_trk).at(i_cent) = new TH1D(name.c_str(), name.c_str(), N_dR, array_dr_bins);

				name = Form("h_cone_data_MB_data_r_1D_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_cone_data_MB_data_r_1D.at(i_jet).at(i_trk).at(i_cent) = new TH1D(name.c_str(), name.c_str(), N_dR, array_dr_bins);


				name = Form("h_cone_MB_r_1D_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_cone_MB_r_1D.at(i_jet).at(i_trk).at(i_cent) = new TH1D(name.c_str(), name.c_str(), N_dR, array_dr_bins);

				name = Form("h_cone_data_MB_r_1D_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_cone_data_MB_r_1D.at(i_jet).at(i_trk).at(i_cent) = new TH1D(name.c_str(), name.c_str(), N_dR, array_dr_bins);

				name = Form("h_TM_MB_r_1D_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_TM_MB_r_1D.at(i_jet).at(i_trk).at(i_cent) = new TH1D(name.c_str(), name.c_str(), N_dR, array_dr_bins);

				name = Form("h_FS_MB_r_1D_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_FS_MB_r_1D.at(i_jet).at(i_trk).at(i_cent) = new TH1D(name.c_str(), name.c_str(), N_dR, array_dr_bins);

				name = Form("h_FNS_MB_r_1D_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_FNS_MB_r_1D.at(i_jet).at(i_trk).at(i_cent) = new TH1D(name.c_str(), name.c_str(), N_dR, array_dr_bins);
			}
		}
	}

	for (int i_cent = 0; i_cent < 6; i_cent++)
	{
		name = Form("MB_norm_jet_cent%i", i_cent);
		TH1* h_jet_spectra_data = (TH1*)((TH1*)input_file_data->Get(name.c_str()))->Clone(Form("reco_jet_y4_c%i", i_cent));
		h_jet_spectra_data->SetName(Form("%s_data",name.c_str()));
		h_jet_spectra_data->Sumw2();

		name = Form("MB_norm_jet_cent%i", i_cent);
		TH1* h_jet_spectra = (TH1*)((TH1*)input_file->Get(name.c_str()))->Clone(Form("reco_jet_UE_norm_y4_c%i", i_cent));
		h_jet_spectra->SetName(Form("%s_mc",name.c_str()));
		h_jet_spectra->Sumw2();

		name = Form("cone_norm_jet_cent%i", i_cent);
		TH1* h_cone_jet_spectra = (TH1*)((TH1*)input_file->Get(name.c_str()))->Clone(Form("cone_UE_norm_y4_c%i", i_cent));
		h_cone_jet_spectra->SetName(Form("%s_mc",name.c_str()));
		h_cone_jet_spectra->Sumw2();

		name = Form("cone_norm_jet_cent%i", i_cent);
		TH1* h_cone_jet_spectra_data = (TH1*)((TH1*)input_file_data->Get(name.c_str()))->Clone(Form("cone_data_UE_norm_y4_c%i", i_cent));
		h_cone_jet_spectra_data->SetName(Form("%s_mc",name.c_str()));
		h_cone_jet_spectra_data->Sumw2();


		name = Form("FS_norm_jet_cent%i", i_cent);
		TH1* h_true_jet_spectra = (TH1*)((TH1*)input_file->Get(name.c_str()))->Clone(Form("true_jet_UE_norm_y4_c%i", i_cent));
		h_true_jet_spectra->Sumw2();


		for (int i_dR = 0; i_dR < N_dR; i_dR++)
		{
			name = Form("ChPS_MB_UE_dR%i_cent%i", i_dR, i_cent);
			h_MB_2D[i_dR][i_cent] = (TH2*)input_file->Get(name.c_str());
			h_MB_2D[i_dR][i_cent]->SetName(Form("ChPS_MB_MC_data_UE_dR%i_cent%i", i_dR, i_cent));

			name = Form("ChPS_MB_UE_dR%i_cent%i", i_dR, i_cent);
			h_MB_data_2D[i_dR][i_cent] = (TH2*)(TH2*)input_file_data->Get(name.c_str());
			h_MB_data_2D[i_dR][i_cent]->SetName(Form("ChPS_MB_data_data_UE_dR%i_cent%i", i_dR, i_cent));

			name = Form("ChPS_cone_UE_dR%i_cent%i", i_dR, i_cent);
			h_cone_2D[i_dR][i_cent] = (TH2*)input_file->Get(name.c_str());

			name = Form("ChPS_cone_UE_dR%i_cent%i", i_dR, i_cent);
			h_cone_data_2D[i_dR][i_cent] = (TH2*)input_file_data->Get(name.c_str());
			h_cone_data_2D[i_dR][i_cent]->SetName(Form("ChPS_cone_data_UE_dR%i_cent%i", i_dR, i_cent));

			name = Form("ChPS_TM_UE_dR%i_cent%i", i_dR, i_cent);
			h_TM_2D[i_dR][i_cent] = (TH2*)input_file->Get(name.c_str());

			name = Form("ChPS_FS_UE_dR%i_cent%i",i_dR,i_cent);
			h_FS_2D[i_dR][i_cent] = (TH2*)input_file->Get(name.c_str());

			name = Form("ChPS_FNS_UE_dR%i_cent%i",i_dR,i_cent);
			h_FNS_2D[i_dR][i_cent] = (TH2*)input_file->Get(name.c_str());

			for (int i_jet = 1; i_jet <= N_jetpt; i_jet++)
			{

				double n_jets_data = h_jet_spectra_data->GetBinContent(i_jet);
				double n_jets_cone = h_cone_jet_spectra->GetBinContent(i_jet);
				double n_jets_cone_data = h_cone_jet_spectra_data->GetBinContent(i_jet);
				double n_jets = h_jet_spectra->GetBinContent(i_jet);
				double n_jets_true = h_true_jet_spectra->GetBinContent(i_jet);



				for (int i_trk = 1; i_trk <= N_trkpt; i_trk++)
				{

					if (n_jets != 0)
					{
						double updated_UE_MB = h_MB_2D[i_dR][i_cent]->GetBinContent(i_trk, i_jet) / n_jets;
						double updated_UE_MB_err = h_MB_2D[i_dR][i_cent]->GetBinError(i_trk, i_jet) / n_jets;
						h_MB_2D[i_dR][i_cent]->SetBinContent(i_trk, i_jet, updated_UE_MB);
						h_MB_2D[i_dR][i_cent]->SetBinError(i_trk, i_jet, updated_UE_MB_err);
					}
					else
					{
						h_MB_2D[i_dR][i_cent]->SetBinContent(i_trk, i_jet, 0);
						h_MB_2D[i_dR][i_cent]->SetBinError(i_trk, i_jet, 0);
					}

					if (n_jets_data != 0)
					{

						double updated_UE_MB_data = h_MB_data_2D[i_dR][i_cent]->GetBinContent(i_trk, i_jet) / n_jets_data;
						double updated_UE_MB_data_err = h_MB_data_2D[i_dR][i_cent]->GetBinError(i_trk, i_jet) / n_jets_data;
						h_MB_data_2D[i_dR][i_cent]->SetBinContent(i_trk, i_jet, updated_UE_MB_data);
						h_MB_data_2D[i_dR][i_cent]->SetBinError(i_trk, i_jet, updated_UE_MB_data_err);
					}
					else
					{
						h_MB_data_2D[i_dR][i_cent]->SetBinContent(i_trk, i_jet, 0);
						h_MB_data_2D[i_dR][i_cent]->SetBinError(i_trk, i_jet, 0);
					}

					if (n_jets_cone != 0)
					{
						double updated_UE_cone = h_cone_2D[i_dR][i_cent]->GetBinContent(i_trk, i_jet) / n_jets_cone;
						double updated_UE_cone_err = h_cone_2D[i_dR][i_cent]->GetBinError(i_trk, i_jet) / n_jets_cone;
						h_cone_2D[i_dR][i_cent]->SetBinContent(i_trk, i_jet, updated_UE_cone);
						h_cone_2D[i_dR][i_cent]->SetBinError(i_trk, i_jet, updated_UE_cone_err);
					}
					else
					{
						h_cone_2D[i_dR][i_cent]->SetBinContent(i_trk, i_jet, 0);
						h_cone_2D[i_dR][i_cent]->SetBinError(i_trk, i_jet, 0);
					}


					if (n_jets_cone_data != 0)
					{
						double updated_UE_cone_data = h_cone_data_2D[i_dR][i_cent]->GetBinContent(i_trk, i_jet) / n_jets_cone_data;
						double updated_UE_cone_data_err = h_cone_data_2D[i_dR][i_cent]->GetBinError(i_trk, i_jet) / n_jets_cone_data;
						h_cone_data_2D[i_dR][i_cent]->SetBinContent(i_trk, i_jet, updated_UE_cone_data);
						h_cone_data_2D[i_dR][i_cent]->SetBinError(i_trk, i_jet, updated_UE_cone_data_err);
					}
					else
					{
						h_cone_data_2D[i_dR][i_cent]->SetBinContent(i_trk, i_jet, 0);
						h_cone_data_2D[i_dR][i_cent]->SetBinError(i_trk, i_jet, 0);
					}

					if (n_jets != 0)
					{
						double updated_UE_TM = h_TM_2D[i_dR][i_cent]->GetBinContent(i_trk, i_jet) / n_jets;
						double updated_UE_TM_err = h_TM_2D[i_dR][i_cent]->GetBinError(i_trk, i_jet) / n_jets;
						h_TM_2D[i_dR][i_cent]->SetBinContent(i_trk, i_jet, updated_UE_TM);
						h_TM_2D[i_dR][i_cent]->SetBinError(i_trk, i_jet, updated_UE_TM_err);
					}
					else
					{
						h_TM_2D[i_dR][i_cent]->SetBinContent(i_trk, i_jet, 0);
						h_TM_2D[i_dR][i_cent]->SetBinError(i_trk, i_jet, 0);
					}

					if (n_jets_true != 0)
					{
						double updated_UE_FS = h_FS_2D[i_dR][i_cent]->GetBinContent(i_trk, i_jet) / n_jets_true;
						double updated_UE_FS_err = h_FS_2D[i_dR][i_cent]->GetBinError(i_trk, i_jet) / n_jets_true;
						h_FS_2D[i_dR][i_cent]->SetBinContent(i_trk, i_jet, updated_UE_FS);
						h_FS_2D[i_dR][i_cent]->SetBinError(i_trk, i_jet, updated_UE_FS_err);
					}
					else
					{
						h_FS_2D[i_dR][i_cent]->SetBinContent(i_trk, i_jet, 0);
						h_FS_2D[i_dR][i_cent]->SetBinError(i_trk, i_jet, 0);
					}

					if (n_jets_true != 0)
					{
						double updated_UE_FNS = h_FNS_2D[i_dR][i_cent]->GetBinContent(i_trk, i_jet) / n_jets_true;
						double updated_UE_FNS_err = h_FNS_2D[i_dR][i_cent]->GetBinError(i_trk, i_jet) / n_jets_true;
						h_FNS_2D[i_dR][i_cent]->SetBinContent(i_trk, i_jet, updated_UE_FNS);
						h_FNS_2D[i_dR][i_cent]->SetBinError(i_trk, i_jet, updated_UE_FNS_err);
					}
					else
					{
						h_FNS_2D[i_dR][i_cent]->SetBinContent(i_trk, i_jet, 0);
						h_FNS_2D[i_dR][i_cent]->SetBinError(i_trk, i_jet, 0);
					}

				}
			}


			//has to be integrated over jet cone before UE-JER correction. the correction will then be based on the integrated TM/Cone_MC
			if (dR_binning->GetBinLowEdge(i_dR+1) < 0.4)
			{

				if (i_dR == 0)
				{
					h_cone_injet_2D[i_cent] = (TH2*)h_cone_2D[i_dR][i_cent]->Clone(Form("%s_injet",h_cone_2D[i_dR][i_cent]->GetName()));
					h_cone_data_injet_2D[i_cent] = (TH2*)h_cone_data_2D[i_dR][i_cent]->Clone(Form("%s_injet",h_cone_data_2D[i_dR][i_cent]->GetName()));
					h_MB_injet_2D[i_cent] = (TH2*)h_MB_2D[i_dR][i_cent]->Clone(Form("%s_injet",h_MB_2D[i_dR][i_cent]->GetName()));
					h_MB_data_injet_2D[i_cent] = (TH2*)h_MB_data_2D[i_dR][i_cent]->Clone(Form("%s_injet",h_MB_data_2D[i_dR][i_cent]->GetName()));
					h_TM_injet_2D[i_cent] = (TH2*)h_TM_2D[i_dR][i_cent]->Clone(Form("%s_injet",h_TM_2D[i_dR][i_cent]->GetName()));
				}
				else
				{
					h_cone_injet_2D[i_cent]->Add(h_cone_2D[i_dR][i_cent]);
					h_cone_data_injet_2D[i_cent]->Add(h_cone_data_2D[i_dR][i_cent]);
					h_MB_injet_2D[i_cent]->Add(h_MB_2D[i_dR][i_cent]);
					h_MB_data_injet_2D[i_cent]->Add(h_MB_data_2D[i_dR][i_cent]);
					h_TM_injet_2D[i_cent]->Add(h_TM_2D[i_dR][i_cent]);
				}
			}
			//applying TM/x correction for x_data = cone method in data, x_data = MC method in data
			{
				TH2* h_cone_correction = (TH2*)h_TM_2D[i_dR][i_cent]->Clone(Form("h_cone_correction_dR%i_cent%i",i_dR, i_cent));
				h_cone_correction->Divide(h_cone_2D[i_dR][i_cent]);
				h_cone_data_2D[i_dR][i_cent]->Multiply(h_cone_correction);

				//MB correction not required since TM/MB_MC = 1
//				TH2* h_MB_correction = (TH2*)h_TM_2D[i_dR][i_cent]->Clone(Form("h_MB_correction_dR%i_cent%i",i_dR, i_cent));
//				h_MB_correction->Divide(h_MB_2D[i_dR][i_cent]);
//				h_MB_data_2D[i_dR][i_cent]->Multiply(h_MB_correction);

				output_file->cd();
				h_cone_correction->Write(Form("h_cone_correction_dR%i_cent%i",i_dR, i_cent));
//				h_MB_correction->Write(Form("h_MB_correction_dR%i_cent%i",i_dR, i_cent));

			}



			//compare to MB
			if (subtract)
			{
				name = Form("sub_cone_data_MB_dR%i_cent%i", i_dR, i_cent);
				h_cone_data_MB_2D[i_dR][i_cent] = (TH2*)h_cone_data_2D[i_dR][i_cent]->Clone(name.c_str());
				h_cone_data_MB_2D[i_dR][i_cent]->Add(h_MB_2D[i_dR][i_cent],-1);

				name = Form("sub_cone_MB_dR%i_cent%i", i_dR, i_cent);
				h_cone_MB_2D[i_dR][i_cent] = (TH2*)h_cone_2D[i_dR][i_cent]->Clone(name.c_str());
				h_cone_MB_2D[i_dR][i_cent]->Add(h_MB_2D[i_dR][i_cent],-1);

				name = Form("sub_MB_data_MB_dR%i_cent%i", i_dR, i_cent);
				h_MB_data_MB_2D[i_dR][i_cent] = (TH2*)h_MB_data_2D[i_dR][i_cent]->Clone(name.c_str());
				h_MB_data_MB_2D[i_dR][i_cent]->Add(h_MB_2D[i_dR][i_cent],-1);

				name = Form("sub_cone_data_MB_data_dR%i_cent%i", i_dR, i_cent);
				h_cone_data_MB_data_2D[i_dR][i_cent] = (TH2*)h_cone_data_2D[i_dR][i_cent]->Clone(name.c_str());
				h_cone_data_MB_data_2D[i_dR][i_cent]->Add(h_MB_data_2D[i_dR][i_cent],-1);

				name = Form("sub_TM_MB_dR%i_cent%i", i_dR, i_cent);
				h_TM_MB_2D[i_dR][i_cent] = (TH2*)h_TM_2D[i_dR][i_cent]->Clone(name.c_str());
				h_TM_MB_2D[i_dR][i_cent]->Add(h_MB_2D[i_dR][i_cent],-1);

				name = Form("sub_FS_MB_dR%i_cent%i", i_dR, i_cent);
				h_FS_MB_2D[i_dR][i_cent] = (TH2*)h_FS_2D[i_dR][i_cent]->Clone(name.c_str());
				h_FS_MB_2D[i_dR][i_cent]->Add(h_MB_2D[i_dR][i_cent],-1);

				name = Form("sub_FNS_MB_dR%i_cent%i", i_dR, i_cent);
				h_FNS_MB_2D[i_dR][i_cent] = (TH2*)h_FNS_2D[i_dR][i_cent]->Clone(name.c_str());
				h_FNS_MB_2D[i_dR][i_cent]->Add(h_MB_2D[i_dR][i_cent],-1);
			}

			else
			{
				name = Form("ratio_cone_data_MB_dR%i_cent%i", i_dR, i_cent);
				h_cone_data_MB_2D[i_dR][i_cent] = (TH2*)h_cone_data_2D[i_dR][i_cent]->Clone(name.c_str());
				h_cone_data_MB_2D[i_dR][i_cent]->Divide(h_MB_2D[i_dR][i_cent]);

				name = Form("ratio_cone_MB_dR%i_cent%i", i_dR, i_cent);
				h_cone_MB_2D[i_dR][i_cent] = (TH2*)h_cone_2D[i_dR][i_cent]->Clone(name.c_str());
				h_cone_MB_2D[i_dR][i_cent]->Divide(h_MB_2D[i_dR][i_cent]);

				name = Form("ratio_MB_data_MB_dR%i_cent%i", i_dR, i_cent);
				h_MB_data_MB_2D[i_dR][i_cent] = (TH2*)h_MB_data_2D[i_dR][i_cent]->Clone(name.c_str());
				h_MB_data_MB_2D[i_dR][i_cent]->Divide(h_MB_2D[i_dR][i_cent]);

				name = Form("ratio_cone_data_MB_data_MB_dR%i_cent%i", i_dR, i_cent);
				h_cone_data_MB_data_2D[i_dR][i_cent] = (TH2*)h_cone_data_2D[i_dR][i_cent]->Clone(name.c_str());
				h_cone_data_MB_data_2D[i_dR][i_cent]->Divide(h_MB_data_2D[i_dR][i_cent]);

				name = Form("ratio_TM_MB_dR%i_cent%i", i_dR, i_cent);
				h_TM_MB_2D[i_dR][i_cent] = (TH2*)h_TM_2D[i_dR][i_cent]->Clone(name.c_str());
				h_TM_MB_2D[i_dR][i_cent]->Divide(h_MB_2D[i_dR][i_cent]);

				name = Form("ratio_FS_MB_dR%i_cent%i", i_dR, i_cent);
				h_FS_MB_2D[i_dR][i_cent] = (TH2*)h_FS_2D[i_dR][i_cent]->Clone(name.c_str());
				h_FS_MB_2D[i_dR][i_cent]->Divide(h_MB_2D[i_dR][i_cent]);

				name = Form("ratio_FNS_MB_dR%i_cent%i", i_dR, i_cent);
				h_FNS_MB_2D[i_dR][i_cent] = (TH2*)h_FNS_2D[i_dR][i_cent]->Clone(name.c_str());
				h_FNS_MB_2D[i_dR][i_cent]->Divide(h_MB_2D[i_dR][i_cent]);
			}


			double dR_lo = dR_binning->GetBinLowEdge(i_dR+1);
			double dR_hi = dR_binning->GetBinUpEdge(i_dR+1);
			double area = TMath::Pi() * ((dR_hi*dR_hi) - (dR_lo*dR_lo));

			for (int i_jet = 0; i_jet < N_jetpt; i_jet++)
			{
				h_cone_1D[i_jet][i_dR][i_cent] = (TH1*)h_cone_2D[i_dR][i_cent]->ProjectionX(Form("cone_jet%i_dr%i_cent%i", i_jet, i_dR, i_cent), i_jet+1, i_jet+1);
				h_cone_data_1D[i_jet][i_dR][i_cent] = (TH1*)h_cone_data_2D[i_dR][i_cent]->ProjectionX(Form("cone_data_jet%i_dr%i_cent%i", i_jet, i_dR, i_cent), i_jet+1, i_jet+1);
				h_MB_1D[i_jet][i_dR][i_cent] = (TH1*)h_MB_2D[i_dR][i_cent]->ProjectionX(Form("MB_jet%i_dr%i_cent%i", i_jet, i_dR, i_cent), i_jet+1, i_jet+1);
				h_MB_data_1D[i_jet][i_dR][i_cent] = (TH1*)h_MB_data_2D[i_dR][i_cent]->ProjectionX(Form("MB_data_jet%i_dr%i_cent%i", i_jet, i_dR, i_cent), i_jet+1, i_jet+1);
				h_TM_1D[i_jet][i_dR][i_cent] = (TH1*)h_TM_2D[i_dR][i_cent]->ProjectionX(Form("TM_jet%i_dr%i_cent%i", i_jet, i_dR, i_cent), i_jet+1, i_jet+1);
				h_FS_1D[i_jet][i_dR][i_cent] = (TH1*)h_FS_2D[i_dR][i_cent]->ProjectionX(Form("FS_jet%i_dr%i_cent%i", i_jet, i_dR, i_cent), i_jet+1, i_jet+1);
				h_FNS_1D[i_jet][i_dR][i_cent] = (TH1*)h_FNS_2D[i_dR][i_cent]->ProjectionX(Form("FNS_jet%i_dr%i_cent%i", i_jet, i_dR, i_cent), i_jet+1, i_jet+1);

				h_cone_1D[i_jet][i_dR][i_cent]->Scale(1.,"width");
				h_cone_data_1D[i_jet][i_dR][i_cent]->Scale(1.,"width");
				h_MB_1D[i_jet][i_dR][i_cent]->Scale(1.,"width");
				h_MB_data_1D[i_jet][i_dR][i_cent]->Scale(1.,"width");
				h_TM_1D[i_jet][i_dR][i_cent]->Scale(1.,"width");
				h_FS_1D[i_jet][i_dR][i_cent]->Scale(1.,"width");
				h_FNS_1D[i_jet][i_dR][i_cent]->Scale(1.,"width");

				h_cone_1D[i_jet][i_dR][i_cent]->Scale(1./area);
				h_cone_data_1D[i_jet][i_dR][i_cent]->Scale(1./area);
				h_MB_1D[i_jet][i_dR][i_cent]->Scale(1./area);
				h_MB_data_1D[i_jet][i_dR][i_cent]->Scale(1./area);
				h_TM_1D[i_jet][i_dR][i_cent]->Scale(1./area);
				h_FS_1D[i_jet][i_dR][i_cent]->Scale(1./area);
				h_FNS_1D[i_jet][i_dR][i_cent]->Scale(1./area);

				//dont scale ratios if dividing
				h_cone_MB_1D[i_jet][i_dR][i_cent] = (TH1*)h_cone_MB_2D[i_dR][i_cent]->ProjectionX(Form("cone_MB_jet%i_dr%i_cent%i", i_jet, i_dR, i_cent), i_jet+1, i_jet+1);
				h_cone_data_MB_1D[i_jet][i_dR][i_cent] = (TH1*)h_cone_data_MB_2D[i_dR][i_cent]->ProjectionX(Form("cone_data_MB_jet%i_dr%i_cent%i", i_jet, i_dR, i_cent), i_jet+1, i_jet+1);
				h_MB_data_MB_1D[i_jet][i_dR][i_cent] = (TH1*)h_MB_data_MB_2D[i_dR][i_cent]->ProjectionX(Form("MB_data_MB_jet%i_dr%i_cent%i", i_jet, i_dR, i_cent), i_jet+1, i_jet+1);
				h_cone_data_MB_data_1D[i_jet][i_dR][i_cent] = (TH1*)h_cone_data_MB_data_2D[i_dR][i_cent]->ProjectionX(Form("cone_data_MB_data_jet%i_dr%i_cent%i", i_jet, i_dR, i_cent), i_jet+1, i_jet+1);
				h_TM_MB_1D[i_jet][i_dR][i_cent] = (TH1*)h_TM_MB_2D[i_dR][i_cent]->ProjectionX(Form("TM_MB_jet%i_dr%i_cent%i", i_jet, i_dR, i_cent), i_jet+1, i_jet+1);
				h_FS_MB_1D[i_jet][i_dR][i_cent] = (TH1*)h_FS_MB_2D[i_dR][i_cent]->ProjectionX(Form("FS_MB_jet%i_dr%i_cent%i", i_jet, i_dR, i_cent), i_jet+1, i_jet+1);
				h_FNS_MB_1D[i_jet][i_dR][i_cent] = (TH1*)h_FNS_MB_2D[i_dR][i_cent]->ProjectionX(Form("FNS_MB_jet%i_dr%i_cent%i", i_jet, i_dR, i_cent), i_jet+1, i_jet+1);

				if (subtract)
				{
					h_cone_data_MB_1D[i_jet][i_dR][i_cent]->Scale(1.,"width");
					h_cone_MB_1D[i_jet][i_dR][i_cent]->Scale(1.,"width");
					h_MB_data_MB_1D[i_jet][i_dR][i_cent]->Scale(1.,"width");
					h_cone_data_MB_data_1D[i_jet][i_dR][i_cent]->Scale(1.,"width");
					h_TM_MB_1D[i_jet][i_dR][i_cent]->Scale(1.,"width");
					h_FS_MB_1D[i_jet][i_dR][i_cent]->Scale(1.,"width");
					h_FNS_MB_1D[i_jet][i_dR][i_cent]->Scale(1.,"width");

					h_cone_data_MB_1D[i_jet][i_dR][i_cent]->Scale(1./area);
					h_cone_MB_1D[i_jet][i_dR][i_cent]->Scale(1./area);
					h_MB_data_MB_1D[i_jet][i_dR][i_cent]->Scale(1./area);
					h_cone_data_MB_data_1D[i_jet][i_dR][i_cent]->Scale(1./area);
					h_TM_MB_1D[i_jet][i_dR][i_cent]->Scale(1./area);
					h_FS_MB_1D[i_jet][i_dR][i_cent]->Scale(1./area);
					h_FNS_MB_1D[i_jet][i_dR][i_cent]->Scale(1./area);
				}

			}
		}

		//recast in terms of r
		for (int i_jet = 0; i_jet < N_jetpt; i_jet++)
		{
			for (int i_trk = 0; i_trk < N_trkpt; i_trk++)
			{
				for (int i_dR = 0; i_dR < N_dR; i_dR++)
				{
					h_cone_r_1D[i_jet][i_trk][i_cent]->SetBinContent(i_dR+1, h_cone_1D[i_jet][i_dR][i_cent]->GetBinContent(i_trk+1));
					h_cone_data_r_1D[i_jet][i_trk][i_cent]->SetBinContent(i_dR+1, h_cone_data_1D[i_jet][i_dR][i_cent]->GetBinContent(i_trk+1));
					h_MB_r_1D[i_jet][i_trk][i_cent]->SetBinContent(i_dR+1, h_MB_1D[i_jet][i_dR][i_cent]->GetBinContent(i_trk+1));
					h_MB_data_r_1D[i_jet][i_trk][i_cent]->SetBinContent(i_dR+1, h_MB_data_1D[i_jet][i_dR][i_cent]->GetBinContent(i_trk+1));
					h_TM_r_1D[i_jet][i_trk][i_cent]->SetBinContent(i_dR+1, h_TM_1D[i_jet][i_dR][i_cent]->GetBinContent(i_trk+1));
					h_FS_r_1D[i_jet][i_trk][i_cent]->SetBinContent(i_dR+1, h_FS_1D[i_jet][i_dR][i_cent]->GetBinContent(i_trk+1));
					h_FNS_r_1D[i_jet][i_trk][i_cent]->SetBinContent(i_dR+1, h_FNS_1D[i_jet][i_dR][i_cent]->GetBinContent(i_trk+1));

					h_cone_r_1D[i_jet][i_trk][i_cent]->SetBinError(i_dR+1, h_cone_1D[i_jet][i_dR][i_cent]->GetBinError(i_trk+1));
					h_cone_data_r_1D[i_jet][i_trk][i_cent]->SetBinError(i_dR+1, h_cone_data_1D[i_jet][i_dR][i_cent]->GetBinError(i_trk+1));
					h_MB_r_1D[i_jet][i_trk][i_cent]->SetBinError(i_dR+1, h_MB_1D[i_jet][i_dR][i_cent]->GetBinError(i_trk+1));
					h_MB_data_r_1D[i_jet][i_trk][i_cent]->SetBinError(i_dR+1, h_MB_data_1D[i_jet][i_dR][i_cent]->GetBinError(i_trk+1));
					h_TM_r_1D[i_jet][i_trk][i_cent]->SetBinError(i_dR+1, h_TM_1D[i_jet][i_dR][i_cent]->GetBinError(i_trk+1));
					h_FS_r_1D[i_jet][i_trk][i_cent]->SetBinError(i_dR+1, h_FS_1D[i_jet][i_dR][i_cent]->GetBinError(i_trk+1));
					h_FNS_r_1D[i_jet][i_trk][i_cent]->SetBinError(i_dR+1, h_FNS_1D[i_jet][i_dR][i_cent]->GetBinError(i_trk+1));

					h_cone_r_1D[i_jet][i_trk][i_cent]->GetXaxis()->SetTitle("r");
					h_cone_data_r_1D[i_jet][i_trk][i_cent]->GetXaxis()->SetTitle("r");
					h_MB_r_1D[i_jet][i_trk][i_cent]->GetXaxis()->SetTitle("r");
					h_MB_data_r_1D[i_jet][i_trk][i_cent]->GetXaxis()->SetTitle("r");
					h_TM_r_1D[i_jet][i_trk][i_cent]->GetXaxis()->SetTitle("r");
					h_FS_r_1D[i_jet][i_trk][i_cent]->GetXaxis()->SetTitle("r");
					h_FNS_r_1D[i_jet][i_trk][i_cent]->GetXaxis()->SetTitle("r");

					h_cone_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetTitle("D_{p_{T}} [UE]");
					h_cone_data_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetTitle("D_{p_{T}} [UE]");
					h_MB_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetTitle("D_{p_{T}} [UE]");
					h_MB_data_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetTitle("D_{p_{T}} [UE]");
					h_TM_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetTitle("D_{p_{T}} [UE]");
					h_FS_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetTitle("D_{p_{T}} [UE]");
					h_FNS_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetTitle("D_{p_{T}} [UE]");


					//ratios
					h_cone_data_MB_r_1D[i_jet][i_trk][i_cent]->SetBinContent(i_dR+1, h_cone_data_MB_1D[i_jet][i_dR][i_cent]->GetBinContent(i_trk+1));
					h_cone_MB_r_1D[i_jet][i_trk][i_cent]->SetBinContent(i_dR+1, h_cone_MB_1D[i_jet][i_dR][i_cent]->GetBinContent(i_trk+1));
					h_MB_data_MB_r_1D[i_jet][i_trk][i_cent]->SetBinContent(i_dR+1, h_MB_data_MB_1D[i_jet][i_dR][i_cent]->GetBinContent(i_trk+1));
					h_cone_data_MB_data_r_1D[i_jet][i_trk][i_cent]->SetBinContent(i_dR+1, h_cone_data_MB_data_1D[i_jet][i_dR][i_cent]->GetBinContent(i_trk+1));
					h_TM_MB_r_1D[i_jet][i_trk][i_cent]->SetBinContent(i_dR+1, h_TM_MB_1D[i_jet][i_dR][i_cent]->GetBinContent(i_trk+1));
					h_FS_MB_r_1D[i_jet][i_trk][i_cent]->SetBinContent(i_dR+1, h_FS_MB_1D[i_jet][i_dR][i_cent]->GetBinContent(i_trk+1));
					h_FNS_MB_r_1D[i_jet][i_trk][i_cent]->SetBinContent(i_dR+1, h_FNS_MB_1D[i_jet][i_dR][i_cent]->GetBinContent(i_trk+1));

					h_cone_data_MB_r_1D[i_jet][i_trk][i_cent]->SetBinError(i_dR+1, h_cone_data_MB_1D[i_jet][i_dR][i_cent]->GetBinError(i_trk+1));
					h_cone_MB_r_1D[i_jet][i_trk][i_cent]->SetBinError(i_dR+1, h_cone_MB_1D[i_jet][i_dR][i_cent]->GetBinError(i_trk+1));
					h_MB_data_MB_r_1D[i_jet][i_trk][i_cent]->SetBinError(i_dR+1, h_MB_data_MB_1D[i_jet][i_dR][i_cent]->GetBinError(i_trk+1));
					h_cone_data_MB_data_r_1D[i_jet][i_trk][i_cent]->SetBinError(i_dR+1, h_cone_data_MB_data_1D[i_jet][i_dR][i_cent]->GetBinError(i_trk+1));
					h_TM_MB_r_1D[i_jet][i_trk][i_cent]->SetBinError(i_dR+1, h_TM_MB_1D[i_jet][i_dR][i_cent]->GetBinError(i_trk+1));
					h_FS_MB_r_1D[i_jet][i_trk][i_cent]->SetBinError(i_dR+1, h_FS_MB_1D[i_jet][i_dR][i_cent]->GetBinError(i_trk+1));
					h_FNS_MB_r_1D[i_jet][i_trk][i_cent]->SetBinError(i_dR+1, h_FNS_MB_1D[i_jet][i_dR][i_cent]->GetBinError(i_trk+1));

					h_cone_data_MB_r_1D[i_jet][i_trk][i_cent]->GetXaxis()->SetTitle("r");
					h_cone_MB_r_1D[i_jet][i_trk][i_cent]->GetXaxis()->SetTitle("r");
					h_MB_data_MB_r_1D[i_jet][i_trk][i_cent]->GetXaxis()->SetTitle("r");
					h_cone_data_MB_data_r_1D[i_jet][i_trk][i_cent]->GetXaxis()->SetTitle("r");
					h_TM_MB_r_1D[i_jet][i_trk][i_cent]->GetXaxis()->SetTitle("r");
					h_FS_MB_r_1D[i_jet][i_trk][i_cent]->GetXaxis()->SetTitle("r");
					h_FNS_MB_r_1D[i_jet][i_trk][i_cent]->GetXaxis()->SetTitle("r");

					h_cone_data_MB_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetTitle("Ratio");
					h_cone_MB_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetTitle("Ratio");
					h_MB_data_MB_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetTitle("Ratio");
					h_cone_data_MB_data_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetTitle("Ratio");
					h_TM_MB_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetTitle("Ratio");
					h_FS_MB_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetTitle("Ratio");
					h_FNS_MB_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetTitle("Ratio");

				}
			}
		}


		output_file->cd();

		name = Form("h_reco_jet_spectrum_y4_cent%i", i_cent);
		h_jet_spectra_data->Write(name.c_str());

		name = Form("h_reco_jet_spectrum_UE_norm_y4_cent%i", i_cent);
		h_jet_spectra->Write(name.c_str());

		name = Form("h_cone_jet_spectrum_UE_norm_y4_cent%i", i_cent);
		h_cone_jet_spectra->Write(name.c_str());

		name = Form("h_cone_data_jet_spectrum_UE_norm_y4_cent%i", i_cent);
		h_cone_jet_spectra_data->Write(name.c_str());

		name = Form("h_true_jet_spectrum_UE_norm_y4_cent%i", i_cent);
		h_true_jet_spectra->Write(name.c_str());

		for (int i_dR = 0; i_dR < N_dR; i_dR++)
		{
			h_cone_2D[i_dR][i_cent]->Write(Form("cone_2d_dR%i_cent%i",i_dR, i_cent ));
			h_cone_data_2D[i_dR][i_cent]->Write(Form("cone_data_2d_dR%i_cent%i",i_dR, i_cent ));
			h_MB_data_2D[i_dR][i_cent]->Write(Form("MB_data_2d_dR%i_cent%i",i_dR, i_cent ));

			h_TM_2D[i_dR][i_cent]->Write(Form("TM_2d_dR%i_cent%i",i_dR, i_cent ));

			for (int i_jet = 0; i_jet < N_jetpt; i_jet++)
			{
				name = Form("h_ChPS_UE_cone_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet);
				h_cone_data_1D[i_jet][i_dR][i_cent]->Write(name.c_str());

				name = Form("h_ChPS_UE_MB_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet);
				h_MB_data_1D[i_jet][i_dR][i_cent]->Write(name.c_str());
			}

		}

		h_cone_injet_2D[i_cent]->Write(Form("Cone_Cent%i",i_cent));
		h_cone_data_injet_2D[i_cent]->Write(Form("Cone_data_Cent%i",i_cent));
		h_MB_injet_2D[i_cent]->Write(Form("MC_Cent%i",i_cent));
		h_MB_data_injet_2D[i_cent]->Write(Form("MC_data_Cent%i",i_cent));
		h_TM_injet_2D[i_cent]->Write(Form("TM_Cent%i",i_cent));


		for (int i_jet = 0; i_jet < N_jetpt; i_jet++)
		{
			for (int i_trk = 0; i_trk < N_trkpt; i_trk++)
			{
				output_file->cd();

				name = Form("jet%i_trk%i_cent%i", i_jet, i_trk, i_cent);
				h_cone_r_1D[i_jet][i_trk][i_cent]->SetName(Form("cone_MB_indR_%s", name.c_str()));
				h_cone_data_r_1D[i_jet][i_trk][i_cent]->SetName(Form("UE_cone_data_indR_%s", name.c_str()));
				h_MB_r_1D[i_jet][i_trk][i_cent]->SetName(Form("UE_MB_indR_%s", name.c_str()));
				h_MB_data_r_1D[i_jet][i_trk][i_cent]->SetName(Form("UE_MB_data_indR_%s", name.c_str()));
				h_TM_r_1D[i_jet][i_trk][i_cent]->SetName(Form("UE_TM_indR_%s", name.c_str()));
				h_FS_r_1D[i_jet][i_trk][i_cent]->SetName(Form("UE_FS_indR_%s", name.c_str()));
				h_FNS_r_1D[i_jet][i_trk][i_cent]->SetName(Form("UE_FNS_indR_%s", name.c_str()));

				h_cone_data_MB_r_1D[i_jet][i_trk][i_cent]->SetName(Form("UE_cone_data_MB_indR_%s", name.c_str()));
				h_cone_MB_r_1D[i_jet][i_trk][i_cent]->SetName(Form("UE_cone_MB_indR_%s", name.c_str()));
				h_MB_data_MB_r_1D[i_jet][i_trk][i_cent]->SetName(Form("UE_MB_data_MB_indR_%s", name.c_str()));
				h_cone_data_MB_data_r_1D[i_jet][i_trk][i_cent]->SetName(Form("UE_cone_data_MB_data_indR_%s", name.c_str()));
				h_TM_MB_r_1D[i_jet][i_trk][i_cent]->SetName(Form("UE_TM_MB_indR_%s", name.c_str()));
				h_FS_MB_r_1D[i_jet][i_trk][i_cent]->SetName(Form("UE_FS_MB_indR_%s", name.c_str()));
				h_FNS_MB_r_1D[i_jet][i_trk][i_cent]->SetName(Form("UE_FNS_MB_indR_%s", name.c_str()));


				h_cone_r_1D[i_jet][i_trk][i_cent]->Write(Form("UE_cone_indR_%s", name.c_str()));
				h_cone_data_r_1D[i_jet][i_trk][i_cent]->Write(Form("UE_cone_data_indR_%s", name.c_str()));
				h_MB_r_1D[i_jet][i_trk][i_cent]->Write(Form("UE_MB_indR_%s", name.c_str()));
				h_MB_data_r_1D[i_jet][i_trk][i_cent]->Write(Form("UE_MB_data_indR_%s", name.c_str()));
				h_TM_r_1D[i_jet][i_trk][i_cent]->Write(Form("UE_TM_indR_%s", name.c_str()));
//				h_FS_r_1D[i_jet][i_trk][i_cent]->Write(Form("UE_FS_indR_%s", name.c_str()));
//				h_FNS_r_1D[i_jet][i_trk][i_cent]->Write(Form("UE_FNS_indR_%s", name.c_str()));
//
//				h_cone_data_MB_r_1D[i_jet][i_trk][i_cent]->Write(Form("UE_cone_data_MB_indR_%s", name.c_str()));
//				h_cone_MB_r_1D[i_jet][i_trk][i_cent]->Write(Form("UE_cone_MB_indR_%s", name.c_str()));
//				h_MB_data_MB_r_1D[i_jet][i_trk][i_cent]->Write(Form("UE_MB_data_MB_indR_%s", name.c_str()));
//				h_cone_data_MB_data_r_1D[i_jet][i_trk][i_cent]->Write(Form("UE_cone_data_MB_data_indR_%s", name.c_str()));
//				h_TM_MB_r_1D[i_jet][i_trk][i_cent]->Write(Form("UE_TM_MB_indR_%s", name.c_str()));
//				h_FS_MB_r_1D[i_jet][i_trk][i_cent]->Write(Form("UE_FS_MB_indR_%s", name.c_str()));
//				h_FNS_MB_r_1D[i_jet][i_trk][i_cent]->Write(Form("UE_FNS_MB_indR_%s", name.c_str()));



			}
		}

	}



	TLatex *ltx = new TLatex();
	ltx->SetTextFont(43);
	ltx->SetTextSize(12);
	ltx->SetTextAlign(12);


	{
		//as function of r
		cout << "Function of R" << endl;
		TCanvas *c_x = new TCanvas("c_x","c_x",900,600);

		TLegend *legend_x = new TLegend(0.20, 0.45, 0.40, 0.75, "","brNDC");
		legend_x->SetTextFont(43);
		legend_x->SetBorderSize(0);
		legend_x->SetTextSize(10);
		legend_x->SetNColumns(1);

		TLegend *legend_y = new TLegend(0.50, 0.4, 0.65, 0.450, "","brNDC");
		legend_y->SetTextFont(43);
		legend_y->SetBorderSize(0);
		legend_y->SetTextSize(10);
		legend_y->SetNColumns(1);


		for (int i_trk = 0; i_trk < N_trkpt; i_trk++)
		{
			if (i_trk < 2 || i_trk > 6) continue;

			string trk_label = Form("%1.1f < p_{T}^{Trk} < %1.1f GeV", trkpT_binning->GetBinLowEdge(i_trk+1), trkpT_binning->GetBinUpEdge(i_trk+1));

			int jet_itr = 0;
			for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
			{
				string jet_label = Form("%1.0f < p_{T}^{Jet} < %1.0f GeV", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));
				c_x->Clear();
				c_x->Divide(3,2);

				for (int i_cent = 0; i_cent < 6; i_cent++)
				{

					SetHStyle_smallify(h_MB_r_1D[i_jet][i_trk][i_cent], 0, 1);
					SetHStyle_smallify(h_TM_r_1D[i_jet][i_trk][i_cent], 1, 1);
					SetHStyle_smallify(h_FS_r_1D[i_jet][i_trk][i_cent], 2, 1);
					SetHStyle_smallify(h_FNS_r_1D[i_jet][i_trk][i_cent], 3, 1);
					SetHStyle_smallify(h_MB_data_r_1D[i_jet][i_trk][i_cent], 3, 1);
					SetHStyle_smallify(h_cone_r_1D[i_jet][i_trk][i_cent], 4, 1);
					SetHStyle_smallify(h_cone_data_r_1D[i_jet][i_trk][i_cent], 5, 1);

					SetHStyle_smallify(h_TM_MB_r_1D[i_jet][i_trk][i_cent], 1, 1);
					SetHStyle_smallify(h_FS_MB_r_1D[i_jet][i_trk][i_cent], 2, 1);
					SetHStyle_smallify(h_FNS_MB_r_1D[i_jet][i_trk][i_cent], 3, 1);
					SetHStyle_smallify(h_MB_data_MB_r_1D[i_jet][i_trk][i_cent], 3, 1);
					SetHStyle_smallify(h_cone_MB_r_1D[i_jet][i_trk][i_cent], 4, 1);
					SetHStyle_smallify(h_cone_data_MB_r_1D[i_jet][i_trk][i_cent], 5, 1);
					SetHStyle_smallify(h_cone_data_MB_data_r_1D[i_jet][i_trk][i_cent], 6, 1);

					if (jet_itr == 0 && i_trk == 2 && i_cent == 0)
					{
//						legend_x->AddEntry(h_MB_r_1D[i_jet][i_trk][i_cent],"MC","lp");
						legend_x->AddEntry(h_MB_data_r_1D[i_jet][i_trk][i_cent],"MC_{data}","lp");
//						legend_x->AddEntry(h_TM_r_1D[i_jet][i_trk][i_cent],"TM","lp");
//						legend_x->AddEntry(h_cone_r_1D[i_jet][i_trk][i_cent],"cone","lp");
						legend_x->AddEntry(h_cone_data_r_1D[i_jet][i_trk][i_cent],"cone_{data}","lp");
						legend_y->AddEntry(h_cone_data_MB_data_r_1D[i_jet][i_trk][i_cent],"cone_{data}/MC_{data}","lp");
//						legend_x->AddEntry(h_FS_r_1D[i_jet][i_trk][i_cent],"FS","lp");
//						legend_x->AddEntry(h_FNS_r_1D[i_jet][i_trk][i_cent],"FNS","lp");
					}

					h_TM_r_1D[i_jet][i_trk][i_cent]->GetXaxis()->SetRangeUser(0, r_max_range);
					h_MB_r_1D[i_jet][i_trk][i_cent]->GetXaxis()->SetRangeUser(0, r_max_range);
					h_MB_data_r_1D[i_jet][i_trk][i_cent]->GetXaxis()->SetRangeUser(0, r_max_range);
					h_cone_r_1D[i_jet][i_trk][i_cent]->GetXaxis()->SetRangeUser(0, r_max_range);
					h_cone_data_r_1D[i_jet][i_trk][i_cent]->GetXaxis()->SetRangeUser(0, r_max_range);
					h_FS_r_1D[i_jet][i_trk][i_cent]->GetXaxis()->SetRangeUser(0, r_max_range);
					h_FNS_r_1D[i_jet][i_trk][i_cent]->GetXaxis()->SetRangeUser(0, r_max_range);
					h_TM_MB_r_1D[i_jet][i_trk][i_cent]->GetXaxis()->SetRangeUser(0, r_max_range);
					h_cone_data_MB_data_r_1D[i_jet][i_trk][i_cent]->GetXaxis()->SetRangeUser(0, r_max_range);
					h_cone_MB_r_1D[i_jet][i_trk][i_cent]->GetXaxis()->SetRangeUser(0, r_max_range);

					double avg =(h_cone_data_MB_data_r_1D[i_jet][i_trk][i_cent]->GetBinContent(1) + h_cone_data_MB_data_r_1D[i_jet][i_trk][i_cent]->GetBinContent(11))/2;
					if (subtract) h_TM_MB_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetRangeUser(-0.5,0.5);
//					else h_cone_data_MB_data_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetRangeUser(avg*0.98, 1.02*avg);
					h_TM_MB_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetRangeUser(0.99,1.01);
					h_TM_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetNdivisions(504);
					h_TM_MB_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetNdivisions(505);
					h_cone_data_MB_data_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetNdivisions(505);
//					h_cone_data_MB_data_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetRangeUser(0.88,1.12);


					h_cone_MB_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetRangeUser(0.85, 1.15);
					h_cone_data_MB_data_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetNdivisions(505);


					c_x->cd(i_cent+1);
					gPad->Divide(1,2);

					c_x->cd(i_cent+1)->cd(1);
					gPad->SetPad(0,0.40,0.99,0.99);
					gPad->SetTopMargin(0.05);
					gPad->SetBottomMargin(0.0);
					gPad->SetRightMargin(0);

					double low_range, hi_range;
					low_range = h_TM_r_1D[i_jet][i_trk][i_cent]->GetMinimum() * 0.75;
					hi_range = h_TM_r_1D[i_jet][i_trk][i_cent]->GetMaximum() * 1.4;

					h_TM_r_1D[i_jet][i_trk][i_cent]->GetXaxis()->SetRangeUser(0, r_max_range);
					h_MB_data_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetRangeUser(low_range, hi_range);

					if (i_cent == 0) {low_range = 30; hi_range = 100;}
					if (i_cent == 1) {low_range = 20; hi_range = 60;}
					if (i_cent == 2) {low_range = 10; hi_range = 40;}
					if (i_cent == 3) {low_range = 5; hi_range = 25;}
					if (i_cent == 4) {low_range = 5; hi_range = 15;}
					if (i_cent == 5) {low_range = 1.5; hi_range = 4;}

//					if (i_jet == jet_pt_start && i_trk == 2)
					h_TM_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetRangeUser(low_range, hi_range);


//					h_MB_data_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetRangeUser(50,100);
					h_TM_r_1D[i_jet][i_trk][i_cent]->Draw("");
//					h_MB_r_1D[i_jet][i_trk][i_cent]->Draw("hist same ");
					h_MB_data_r_1D[i_jet][i_trk][i_cent]->Draw("same");
//					h_cone_r_1D[i_jet][i_trk][i_cent]->Draw("same");
					h_cone_data_r_1D[i_jet][i_trk][i_cent]->Draw("same");
//					h_FS_r_1D[i_jet][i_trk][i_cent]->Draw("same");
//					h_FNS_r_1D[i_jet][i_trk][i_cent]->Draw("same");
					gPad->SetLogy(0);

					c_x->cd(i_cent+1)->cd(2);
					gPad->SetPad(0,0.0,0.99,0.40);
					gPad->SetTopMargin(0.0);
					gPad->SetBottomMargin(0.30);
					gPad->SetRightMargin(0);



					if (subtract) h_TM_MB_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetTitle("x - MB");
					else  h_TM_MB_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetTitle("x / MC");
					h_TM_MB_r_1D[i_jet][i_trk][i_cent]->Draw("hist text");
//					h_MB_data_MB_r_1D[i_jet][i_trk][i_cent]->Draw("same hist text");
					h_cone_data_MB_data_r_1D[i_jet][i_trk][i_cent]->Draw("hist same text");
//					h_cone_data_MB_data_r_1D[i_jet][i_trk][i_cent]->SetMarkerSize(5);
//					h_cone_MB_r_1D[i_jet][i_trk][i_cent]->Draw("hist");
//					h_cone_data_MB_r_1D[i_jet][i_trk][i_cent]->Draw("same");
//					h_FS_MB_r_1D[i_jet][i_trk][i_cent]->Draw("same");
//					h_FNS_MB_r_1D[i_jet][i_trk][i_cent]->Draw("same");
					if (subtract) line->DrawLine(0, 0, r_max_range, 0);
					else line->DrawLine(0, 1, r_max_range, 1);

					c_x->cd(i_cent+1);
					ltx->SetTextAlign(32);
					ltx->DrawLatexNDC(0.94,0.94,num_to_cent(31,i_cent).c_str());
					ltx->DrawLatexNDC(0.94,0.88,jet_label.c_str());
					ltx->DrawLatexNDC(0.94,0.82,trk_label.c_str());

				}

				c_x->cd(1)->cd(1);
				legend_x->Draw();
				c_x->cd(1)->cd(2);
				legend_y->Draw();


				jet_itr++;

				if (i_trk == 2 && i_jet == jet_pt_start) name = "(";
				else if (i_trk == 6 && i_jet == jet_pt_end - 1) name = ")";
				else name = "";
				if (subtract) c_x->Print(Form("UE_x_subtr_c%i.pdf%s", sys_mode, name.c_str()), Form("Title: trk%i_jet%i", i_trk, i_jet));
				else c_x->Print(Form("UE_x_ratio_c%i.pdf%s", sys_mode, name.c_str()), Form("Title: trk%i_jet%i", i_trk, i_jet));

			}

		}
	}


}
