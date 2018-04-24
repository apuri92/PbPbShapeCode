#include "output_dev/functions/global_variables.h"

void check_UE()
{
	SetAtlasStyle();
	gErrorIgnoreLevel = 3001;

	TFile *input_file = new TFile("output_dev/raw_results/nominal/FF_MC_out_histo_PbPb_5p02_r001.root");
	TFile *input_file_data = new TFile("output_dev/raw_results/nominal/FF_data_out_histo_PbPb_5p02_r001.root");
	TAxis* dR_binning = (TAxis*)((TH3*)input_file->Get("h_dR_change_jetpt0_cent0"))->GetXaxis();
	TAxis* jetpT_binning = (TAxis*)((TH3*)input_file->Get("ChPS_raw_0_dR0_cent0"))->GetYaxis();
	TAxis* trkpT_binning = (TAxis*)((TH3*)input_file->Get("ChPS_raw_0_dR0_cent0"))->GetXaxis();
	TFile *UE_factors = new TFile(Form("output_dev/unfold/output_pdf_nominal/root/UE_factors.root"));
	bool apply_correctionFactors = false;

	string name;

	int jet_pt_start = 7;
	int jet_pt_end = 11;

	int N_dR = dR_binning->GetNbins();
	int N_jetpt = jetpT_binning->GetNbins();
	int N_trkpt = trkpT_binning->GetNbins();

	vector<vector<TH2*>> h_MB_2D = vector<vector<TH2*>> (N_dR, vector<TH2*> (n_cent_cuts));
	vector<vector<TH2*>> h_MB_data_2D = vector<vector<TH2*>> (N_dR, vector<TH2*> (n_cent_cuts));
	vector<vector<TH2*>> h_MB_tj_2D = vector<vector<TH2*>> (N_dR, vector<TH2*> (n_cent_cuts));
	vector<vector<TH2*>> h_TM_2D = vector<vector<TH2*>> (N_dR, vector<TH2*> (n_cent_cuts));
	vector<vector<TH2*>> h_FS_2D = vector<vector<TH2*>> (N_dR, vector<TH2*> (n_cent_cuts));
	vector<vector<TH2*>> h_FNS_2D = vector<vector<TH2*>> (N_dR, vector<TH2*> (n_cent_cuts));

	vector<vector<vector<TH1*>>> h_MB_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_dR, vector<TH1*> (n_cent_cuts)));
	vector<vector<vector<TH1*>>> h_MB_data_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_dR, vector<TH1*> (n_cent_cuts)));
	vector<vector<vector<TH1*>>> h_MB_tj_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_dR, vector<TH1*> (n_cent_cuts)));
	vector<vector<vector<TH1*>>> h_TM_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_dR, vector<TH1*> (n_cent_cuts)));
	vector<vector<vector<TH1*>>> h_FS_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_dR, vector<TH1*> (n_cent_cuts)));
	vector<vector<vector<TH1*>>> h_FNS_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_dR, vector<TH1*> (n_cent_cuts)));

	vector<vector<vector<TH1*>>> h_MB_r_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_trkpt, vector<TH1*> (n_cent_cuts)));
	vector<vector<vector<TH1*>>> h_MB_data_r_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_trkpt, vector<TH1*> (n_cent_cuts)));
	vector<vector<vector<TH1*>>> h_MB_tj_r_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_trkpt, vector<TH1*> (n_cent_cuts)));
	vector<vector<vector<TH1*>>> h_TM_r_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_trkpt, vector<TH1*> (n_cent_cuts)));
	vector<vector<vector<TH1*>>> h_FS_r_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_trkpt, vector<TH1*> (n_cent_cuts)));
	vector<vector<vector<TH1*>>> h_FNS_r_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_trkpt, vector<TH1*> (n_cent_cuts)));


	vector<vector<TH2*>> h_MB_tj_MB_2D = vector<vector<TH2*>> (N_dR, vector<TH2*> (n_cent_cuts));
	vector<vector<TH2*>> h_MB_data_MB_2D = vector<vector<TH2*>> (N_dR, vector<TH2*> (n_cent_cuts));
	vector<vector<TH2*>> h_TM_MB_2D = vector<vector<TH2*>> (N_dR, vector<TH2*> (n_cent_cuts));
	vector<vector<TH2*>> h_FS_MB_2D = vector<vector<TH2*>> (N_dR, vector<TH2*> (n_cent_cuts));
	vector<vector<TH2*>> h_FNS_MB_2D = vector<vector<TH2*>> (N_dR, vector<TH2*> (n_cent_cuts));

	vector<vector<vector<TH1*>>> h_MB_tj_MB_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_dR, vector<TH1*> (n_cent_cuts)));
	vector<vector<vector<TH1*>>> h_MB_data_MB_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_dR, vector<TH1*> (n_cent_cuts)));
	vector<vector<vector<TH1*>>> h_TM_MB_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_dR, vector<TH1*> (n_cent_cuts)));
	vector<vector<vector<TH1*>>> h_FS_MB_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_dR, vector<TH1*> (n_cent_cuts)));
	vector<vector<vector<TH1*>>> h_FNS_MB_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_dR, vector<TH1*> (n_cent_cuts)));

	vector<vector<vector<TH1*>>> h_TM_MB_r_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_trkpt, vector<TH1*> (n_cent_cuts)));
	vector<vector<vector<TH1*>>> h_MB_data_MB_r_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_trkpt, vector<TH1*> (n_cent_cuts)));
	vector<vector<vector<TH1*>>> h_MB_tj_MB_r_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_trkpt, vector<TH1*> (n_cent_cuts)));
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

				name = Form("h_MB_tj_r_1D_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_MB_tj_r_1D.at(i_jet).at(i_trk).at(i_cent) = new TH1D(name.c_str(), name.c_str(), N_dR, array_dr_bins);

				name = Form("h_TM_r_1D_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_TM_r_1D.at(i_jet).at(i_trk).at(i_cent) = new TH1D(name.c_str(), name.c_str(), N_dR, array_dr_bins);

				name = Form("h_FS_r_1D_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_FS_r_1D.at(i_jet).at(i_trk).at(i_cent) = new TH1D(name.c_str(), name.c_str(), N_dR, array_dr_bins);

				name = Form("h_FNS_r_1D_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_FNS_r_1D.at(i_jet).at(i_trk).at(i_cent) = new TH1D(name.c_str(), name.c_str(), N_dR, array_dr_bins);



				name = Form("h_MB_tj_MB_r_1D_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_MB_tj_MB_r_1D.at(i_jet).at(i_trk).at(i_cent) = new TH1D(name.c_str(), name.c_str(), N_dR, array_dr_bins);

				name = Form("h_MB_data_MB_r_1D_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_MB_data_MB_r_1D.at(i_jet).at(i_trk).at(i_cent) = new TH1D(name.c_str(), name.c_str(), N_dR, array_dr_bins);

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
		for (int i_dR = 0; i_dR < N_dR; i_dR++)
		{

			name = Form("h_reco_jet_spectrum_y4_cent%i", i_cent);
			TH1* h_jet_spectra_data = (TH1*)((TH1*)input_file_data->Get(name.c_str()))->Clone(Form("reco_jet_y4_c%i", i_cent));
			h_jet_spectra_data->SetName(Form("%s_data",name.c_str()));
			h_jet_spectra_data->Sumw2();

			name = Form("h_reco_jet_spectrum_y4_cent%i", i_cent);
			TH1* h_jet_spectra = (TH1*)((TH1*)input_file->Get(name.c_str()))->Clone(Form("reco_jet_y4_c%i", i_cent));
			h_jet_spectra->Sumw2();

			name = Form("h_true_jet_spectrum_y4_cent%i", i_cent);
			TH1* h_true_jet_spectra = (TH1*)((TH1*)input_file->Get(name.c_str()))->Clone(Form("true_jet_y4_c%i", i_cent));
			h_true_jet_spectra->Sumw2();

			name = Form("ChPS_MB_UE_dR%i_cent%i", i_dR, i_cent);
			h_MB_2D[i_dR][i_cent] = (TH2*)input_file->Get(name.c_str());

			name = Form("ChPS_MB_UE_dR%i_cent%i", i_dR, i_cent);
			h_MB_data_2D[i_dR][i_cent] = (TH2*)(TH2*)input_file_data->Get(name.c_str())->Clone(Form("ChPS_MB_data_UE_dR%i_cent%i", i_dR, i_cent));

			name = Form("ChPS_MB_UE_truthjet_dR%i_cent%i", i_dR, i_cent);
			h_MB_tj_2D[i_dR][i_cent] = (TH2*)input_file->Get(name.c_str());

			name = Form("ChPS_TM_UE_dR%i_cent%i", i_dR, i_cent);
			h_TM_2D[i_dR][i_cent] = (TH2*)input_file->Get(name.c_str());

			name = Form("ChPS_FS_UE_dR%i_cent%i",i_dR,i_cent);
			h_FS_2D[i_dR][i_cent] = (TH2*)input_file->Get(name.c_str());

			name = Form("ChPS_FNS_UE_dR%i_cent%i",i_dR,i_cent);
			h_FNS_2D[i_dR][i_cent] = (TH2*)input_file->Get(name.c_str());


			if (apply_correctionFactors)
			{
				TH2* h_UE_corr_factors;
				h_UE_corr_factors = (TH2*)UE_factors->Get(Form("UE_ratio_dR%i_cent%i",i_dR, i_cent));;

				for (int i_jet = 1; i_jet <= N_jetpt; i_jet++)
				{
					for (int i_trk = 1; i_trk <= N_trkpt; i_trk++)
					{
						double orig = h_MB_2D[i_dR][i_cent]->GetBinContent(i_trk, i_jet);
						double correction = h_UE_corr_factors->GetBinContent(i_trk, i_jet); //correct number of jets, correct for UE JER correlation

						h_MB_2D[i_dR][i_cent]->SetBinContent(i_trk, i_jet, orig*correction);
					}
				}
			}
			for (int i_jet = 1; i_jet <= N_jetpt; i_jet++)
			{

				double n_jets_data = h_jet_spectra_data->GetBinContent(i_jet);
				double n_jets = h_jet_spectra->GetBinContent(i_jet);
				double n_jets_true = h_true_jet_spectra->GetBinContent(i_jet);

				if (n_jets == 0) continue;
				

				for (int i_trk = 1; i_trk <= N_trkpt; i_trk++)
				{
					double updated_UE_MB = h_MB_2D[i_dR][i_cent]->GetBinContent(i_trk, i_jet) / n_jets;
					double updated_UE_MB_data = h_MB_data_2D[i_dR][i_cent]->GetBinContent(i_trk, i_jet) / n_jets_data;
					double updated_UE_MB_tj = h_MB_tj_2D[i_dR][i_cent]->GetBinContent(i_trk, i_jet) / n_jets_true;
					double updated_UE_TM = h_TM_2D[i_dR][i_cent]->GetBinContent(i_trk, i_jet) / n_jets;
					double updated_UE_FS = h_FS_2D[i_dR][i_cent]->GetBinContent(i_trk, i_jet) / n_jets_true;
					double updated_UE_FNS = h_FNS_2D[i_dR][i_cent]->GetBinContent(i_trk, i_jet) / n_jets_true;

					double updated_UE_MB_err = h_MB_2D[i_dR][i_cent]->GetBinError(i_trk, i_jet) / n_jets;
					double updated_UE_MB_data_err = h_MB_data_2D[i_dR][i_cent]->GetBinError(i_trk, i_jet) / n_jets_data;
					double updated_UE_MB_tj_err = h_MB_tj_2D[i_dR][i_cent]->GetBinError(i_trk, i_jet) / n_jets_true;
					double updated_UE_TM_err = h_TM_2D[i_dR][i_cent]->GetBinError(i_trk, i_jet) / n_jets;
					double updated_UE_FS_err = h_FS_2D[i_dR][i_cent]->GetBinError(i_trk, i_jet) / n_jets_true;
					double updated_UE_FNS_err = h_FNS_2D[i_dR][i_cent]->GetBinError(i_trk, i_jet) / n_jets_true;

					h_MB_2D[i_dR][i_cent]->SetBinContent(i_trk, i_jet, updated_UE_MB);
					h_MB_data_2D[i_dR][i_cent]->SetBinContent(i_trk, i_jet, updated_UE_MB_data);
					h_MB_tj_2D[i_dR][i_cent]->SetBinContent(i_trk, i_jet, updated_UE_MB_tj);
					h_TM_2D[i_dR][i_cent]->SetBinContent(i_trk, i_jet, updated_UE_TM);
					h_FS_2D[i_dR][i_cent]->SetBinContent(i_trk, i_jet, updated_UE_FS);
					h_FNS_2D[i_dR][i_cent]->SetBinContent(i_trk, i_jet, updated_UE_FNS);

					h_MB_2D[i_dR][i_cent]->SetBinError(i_trk, i_jet, updated_UE_MB_err);
					h_MB_data_2D[i_dR][i_cent]->SetBinError(i_trk, i_jet, updated_UE_MB_data_err);
					h_MB_tj_2D[i_dR][i_cent]->SetBinError(i_trk, i_jet, updated_UE_MB_tj_err);
					h_TM_2D[i_dR][i_cent]->SetBinError(i_trk, i_jet, updated_UE_TM_err);
					h_FS_2D[i_dR][i_cent]->SetBinError(i_trk, i_jet, updated_UE_FS_err);
					h_FNS_2D[i_dR][i_cent]->SetBinError(i_trk, i_jet, updated_UE_FNS_err);
				}
			}


			//compare to MB
			name = Form("ratio_MB_data_MB_dR%i_cent%i", i_dR, i_cent);
			h_MB_data_MB_2D[i_dR][i_cent] = (TH2*)h_MB_data_2D[i_dR][i_cent]->Clone(name.c_str());
			h_MB_data_MB_2D[i_dR][i_cent]->Divide(h_MB_2D[i_dR][i_cent]);

			name = Form("ratio_MB_tj_MB_dR%i_cent%i", i_dR, i_cent);
			h_MB_tj_MB_2D[i_dR][i_cent] = (TH2*)h_MB_tj_2D[i_dR][i_cent]->Clone(name.c_str());
			h_MB_tj_MB_2D[i_dR][i_cent]->Divide(h_MB_2D[i_dR][i_cent]);

			name = Form("ratio_TM_MB_dR%i_cent%i", i_dR, i_cent);
			h_TM_MB_2D[i_dR][i_cent] = (TH2*)h_TM_2D[i_dR][i_cent]->Clone(name.c_str());
			h_TM_MB_2D[i_dR][i_cent]->Divide(h_MB_2D[i_dR][i_cent]);

			name = Form("ratio_FS_MB_dR%i_cent%i", i_dR, i_cent);
			h_FS_MB_2D[i_dR][i_cent] = (TH2*)h_FS_2D[i_dR][i_cent]->Clone(name.c_str());
			h_FS_MB_2D[i_dR][i_cent]->Divide(h_MB_2D[i_dR][i_cent]);

			name = Form("ratio_FNS_MB_dR%i_cent%i", i_dR, i_cent);
			h_FNS_MB_2D[i_dR][i_cent] = (TH2*)h_FNS_2D[i_dR][i_cent]->Clone(name.c_str());
			h_FNS_MB_2D[i_dR][i_cent]->Divide(h_MB_2D[i_dR][i_cent]);


			double dR_lo = dR_binning->GetBinLowEdge(i_dR+1);
			double dR_hi = dR_binning->GetBinUpEdge(i_dR+1);
			double area = TMath::Pi() * ((dR_hi*dR_hi) - (dR_lo*dR_lo));


			for (int i_jet = 0; i_jet < N_jetpt; i_jet++)
			{
				h_MB_1D[i_jet][i_dR][i_cent] = (TH1*)h_MB_2D[i_dR][i_cent]->ProjectionX(Form("MB_%i_%i_%i", i_jet, i_dR, i_cent), i_jet+1, i_jet+1);
				h_MB_data_1D[i_jet][i_dR][i_cent] = (TH1*)h_MB_data_2D[i_dR][i_cent]->ProjectionX(Form("MB_data_%i_%i_%i", i_jet, i_dR, i_cent), i_jet+1, i_jet+1);
				h_MB_tj_1D[i_jet][i_dR][i_cent] = (TH1*)h_MB_tj_2D[i_dR][i_cent]->ProjectionX(Form("MB_tj_%i_%i_%i", i_jet, i_dR, i_cent), i_jet+1, i_jet+1);
				h_TM_1D[i_jet][i_dR][i_cent] = (TH1*)h_TM_2D[i_dR][i_cent]->ProjectionX(Form("TM_%i_%i_%i", i_jet, i_dR, i_cent), i_jet+1, i_jet+1);
				h_FS_1D[i_jet][i_dR][i_cent] = (TH1*)h_FS_2D[i_dR][i_cent]->ProjectionX(Form("FS_%i_%i_%i", i_jet, i_dR, i_cent), i_jet+1, i_jet+1);
				h_FNS_1D[i_jet][i_dR][i_cent] = (TH1*)h_FNS_2D[i_dR][i_cent]->ProjectionX(Form("FNS_%i_%i_%i", i_jet, i_dR, i_cent), i_jet+1, i_jet+1);

				h_MB_1D[i_jet][i_dR][i_cent]->Scale(1.,"width");
				h_MB_data_1D[i_jet][i_dR][i_cent]->Scale(1.,"width");
				h_MB_tj_1D[i_jet][i_dR][i_cent]->Scale(1.,"width");
				h_TM_1D[i_jet][i_dR][i_cent]->Scale(1.,"width");
				h_FS_1D[i_jet][i_dR][i_cent]->Scale(1.,"width");
				h_FNS_1D[i_jet][i_dR][i_cent]->Scale(1.,"width");

				h_MB_1D[i_jet][i_dR][i_cent]->Scale(1./area);
				h_MB_data_1D[i_jet][i_dR][i_cent]->Scale(1./area);
				h_MB_tj_1D[i_jet][i_dR][i_cent]->Scale(1./area);
				h_TM_1D[i_jet][i_dR][i_cent]->Scale(1./area);
				h_FS_1D[i_jet][i_dR][i_cent]->Scale(1./area);
				h_FNS_1D[i_jet][i_dR][i_cent]->Scale(1./area);


				//dont scale ratios
				h_MB_data_MB_1D[i_jet][i_dR][i_cent] = (TH1*)h_MB_data_MB_2D[i_dR][i_cent]->ProjectionX(Form("MB_data_MB_%i_%i_%i", i_jet, i_dR, i_cent), i_jet+1, i_jet+1);
				h_MB_tj_MB_1D[i_jet][i_dR][i_cent] = (TH1*)h_MB_tj_MB_2D[i_dR][i_cent]->ProjectionX(Form("MB_tj_MB_%i_%i_%i", i_jet, i_dR, i_cent), i_jet+1, i_jet+1);
				h_TM_MB_1D[i_jet][i_dR][i_cent] = (TH1*)h_TM_MB_2D[i_dR][i_cent]->ProjectionX(Form("TM_MB_%i_%i_%i", i_jet, i_dR, i_cent), i_jet+1, i_jet+1);
				h_FS_MB_1D[i_jet][i_dR][i_cent] = (TH1*)h_FS_MB_2D[i_dR][i_cent]->ProjectionX(Form("TFS_MB_%i_%i_%i", i_jet, i_dR, i_cent), i_jet+1, i_jet+1);
				h_FNS_MB_1D[i_jet][i_dR][i_cent] = (TH1*)h_FNS_MB_2D[i_dR][i_cent]->ProjectionX(Form("FNS_MB_%i_%i_%i", i_jet, i_dR, i_cent), i_jet+1, i_jet+1);
			}

		}

		//recast in terms of r
		for (int i_jet = 0; i_jet < N_jetpt; i_jet++)
		{
			for (int i_trk = 0; i_trk < N_trkpt; i_trk++)
			{
				for (int i_dR = 0; i_dR < N_dR; i_dR++)
				{
					h_MB_r_1D[i_jet][i_trk][i_cent]->SetBinContent(i_dR+1, h_MB_1D[i_jet][i_dR][i_cent]->GetBinContent(i_trk+1));
					h_MB_data_r_1D[i_jet][i_trk][i_cent]->SetBinContent(i_dR+1, h_MB_data_1D[i_jet][i_dR][i_cent]->GetBinContent(i_trk+1));
					h_MB_tj_r_1D[i_jet][i_trk][i_cent]->SetBinContent(i_dR+1, h_MB_tj_1D[i_jet][i_dR][i_cent]->GetBinContent(i_trk+1));
					h_TM_r_1D[i_jet][i_trk][i_cent]->SetBinContent(i_dR+1, h_TM_1D[i_jet][i_dR][i_cent]->GetBinContent(i_trk+1));
					h_FS_r_1D[i_jet][i_trk][i_cent]->SetBinContent(i_dR+1, h_FS_1D[i_jet][i_dR][i_cent]->GetBinContent(i_trk+1));
					h_FNS_r_1D[i_jet][i_trk][i_cent]->SetBinContent(i_dR+1, h_FNS_1D[i_jet][i_dR][i_cent]->GetBinContent(i_trk+1));

					h_MB_r_1D[i_jet][i_trk][i_cent]->SetBinError(i_dR+1, h_MB_1D[i_jet][i_dR][i_cent]->GetBinError(i_trk+1));
					h_MB_data_r_1D[i_jet][i_trk][i_cent]->SetBinError(i_dR+1, h_MB_data_1D[i_jet][i_dR][i_cent]->GetBinError(i_trk+1));
					h_MB_tj_r_1D[i_jet][i_trk][i_cent]->SetBinError(i_dR+1, h_MB_tj_1D[i_jet][i_dR][i_cent]->GetBinError(i_trk+1));
					h_TM_r_1D[i_jet][i_trk][i_cent]->SetBinError(i_dR+1, h_TM_1D[i_jet][i_dR][i_cent]->GetBinError(i_trk+1));
					h_FS_r_1D[i_jet][i_trk][i_cent]->SetBinError(i_dR+1, h_FS_1D[i_jet][i_dR][i_cent]->GetBinError(i_trk+1));
					h_FNS_r_1D[i_jet][i_trk][i_cent]->SetBinError(i_dR+1, h_FNS_1D[i_jet][i_dR][i_cent]->GetBinError(i_trk+1));

					h_MB_r_1D[i_jet][i_trk][i_cent]->GetXaxis()->SetTitle("r");
					h_MB_data_r_1D[i_jet][i_trk][i_cent]->GetXaxis()->SetTitle("r");
					h_MB_tj_r_1D[i_jet][i_trk][i_cent]->GetXaxis()->SetTitle("r");
					h_TM_r_1D[i_jet][i_trk][i_cent]->GetXaxis()->SetTitle("r");
					h_FS_r_1D[i_jet][i_trk][i_cent]->GetXaxis()->SetTitle("r");
					h_FNS_r_1D[i_jet][i_trk][i_cent]->GetXaxis()->SetTitle("r");

					h_MB_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetTitle("D_{p_{T}} [UE]");
					h_MB_data_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetTitle("D_{p_{T}} [UE]");
					h_MB_tj_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetTitle("D_{p_{T}} [UE]");
					h_TM_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetTitle("D_{p_{T}} [UE]");
					h_FS_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetTitle("D_{p_{T}} [UE]");
					h_FNS_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetTitle("D_{p_{T}} [UE]");


					//ratios
					h_MB_data_MB_r_1D[i_jet][i_trk][i_cent]->SetBinContent(i_dR+1, h_MB_data_MB_1D[i_jet][i_dR][i_cent]->GetBinContent(i_trk+1));
					h_MB_tj_MB_r_1D[i_jet][i_trk][i_cent]->SetBinContent(i_dR+1, h_MB_tj_MB_1D[i_jet][i_dR][i_cent]->GetBinContent(i_trk+1));
					h_TM_MB_r_1D[i_jet][i_trk][i_cent]->SetBinContent(i_dR+1, h_TM_MB_1D[i_jet][i_dR][i_cent]->GetBinContent(i_trk+1));
					h_FS_MB_r_1D[i_jet][i_trk][i_cent]->SetBinContent(i_dR+1, h_FS_MB_1D[i_jet][i_dR][i_cent]->GetBinContent(i_trk+1));
					h_FNS_MB_r_1D[i_jet][i_trk][i_cent]->SetBinContent(i_dR+1, h_FNS_MB_1D[i_jet][i_dR][i_cent]->GetBinContent(i_trk+1));

					h_MB_data_MB_r_1D[i_jet][i_trk][i_cent]->SetBinError(i_dR+1, h_MB_data_MB_1D[i_jet][i_dR][i_cent]->GetBinError(i_trk+1));
					h_MB_tj_MB_r_1D[i_jet][i_trk][i_cent]->SetBinError(i_dR+1, h_MB_tj_MB_1D[i_jet][i_dR][i_cent]->GetBinError(i_trk+1));
					h_TM_MB_r_1D[i_jet][i_trk][i_cent]->SetBinError(i_dR+1, h_TM_MB_1D[i_jet][i_dR][i_cent]->GetBinError(i_trk+1));
					h_FS_MB_r_1D[i_jet][i_trk][i_cent]->SetBinError(i_dR+1, h_FS_MB_1D[i_jet][i_dR][i_cent]->GetBinError(i_trk+1));
					h_FNS_MB_r_1D[i_jet][i_trk][i_cent]->SetBinError(i_dR+1, h_FNS_MB_1D[i_jet][i_dR][i_cent]->GetBinError(i_trk+1));

					h_MB_data_MB_r_1D[i_jet][i_trk][i_cent]->GetXaxis()->SetTitle("r");
					h_MB_tj_MB_r_1D[i_jet][i_trk][i_cent]->GetXaxis()->SetTitle("r");
					h_TM_MB_r_1D[i_jet][i_trk][i_cent]->GetXaxis()->SetTitle("r");
					h_FS_MB_r_1D[i_jet][i_trk][i_cent]->GetXaxis()->SetTitle("r");
					h_FNS_MB_r_1D[i_jet][i_trk][i_cent]->GetXaxis()->SetTitle("r");

					h_MB_data_MB_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetTitle("Ratio");
					h_MB_tj_MB_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetTitle("Ratio");
					h_TM_MB_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetTitle("Ratio");
					h_FS_MB_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetTitle("Ratio");
					h_FNS_MB_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetTitle("Ratio");

				}
			}
		}
	}


	TLatex *ltx = new TLatex();
	ltx->SetTextFont(43);
	ltx->SetTextSize(12);
	ltx->SetTextAlign(12);


	//drawing 2D histo
	{
		cout << "drawing 2D histo" << endl;
		TCanvas *c_TM_MB_2D = new TCanvas("c_TM_MB_2D","c_TM_MB_2D",900,600);
		TCanvas *c_FS_MB_2D = new TCanvas("c_FS_MB_2D","c_FS_MB_2D",900,600);
		TCanvas *c_FNS_MB_2D = new TCanvas("c_FNS_MB_2D","c_FNS_MB_2D",900,600);

		for (int i_dR = 0; i_dR < N_dR; i_dR++)
		{
			string dr_label = Form("%1.2f < r < %1.2f", dR_binning->GetBinLowEdge(i_dR+1), dR_binning->GetBinUpEdge(i_dR+1));

			c_TM_MB_2D->Clear(); c_TM_MB_2D->Divide(3,2);
			c_FS_MB_2D->Clear(); c_FS_MB_2D->Divide(3,2);
			c_FNS_MB_2D->Clear(); c_FNS_MB_2D->Divide(3,2);

			for (int i_cent = 0; i_cent < 6; i_cent++)
			{
				h_TM_MB_2D[i_dR][i_cent]->GetYaxis()->SetTitle("p_{T}^{Jet}");
				h_FS_MB_2D[i_dR][i_cent]->GetYaxis()->SetTitle("p_{T}^{Jet}");
				h_FNS_MB_2D[i_dR][i_cent]->GetYaxis()->SetTitle("p_{T}^{Jet}");

				h_TM_MB_2D[i_dR][i_cent]->GetXaxis()->SetTitle("p_{T}^{Trk}");
				h_FS_MB_2D[i_dR][i_cent]->GetXaxis()->SetTitle("p_{T}^{Trk}");
				h_FNS_MB_2D[i_dR][i_cent]->GetXaxis()->SetTitle("p_{T}^{Trk}");

				h_TM_MB_2D[i_dR][i_cent]->GetXaxis()->SetRangeUser(1,10);
				h_FS_MB_2D[i_dR][i_cent]->GetXaxis()->SetRangeUser(1,10);
				h_FNS_MB_2D[i_dR][i_cent]->GetXaxis()->SetRangeUser(1,10);

				double low_range = 126;
				double hi_range = 316;

				h_TM_MB_2D[i_dR][i_cent]->GetYaxis()->SetRangeUser(90, 500);
				h_FS_MB_2D[i_dR][i_cent]->GetYaxis()->SetRangeUser(low_range, hi_range);
				h_FNS_MB_2D[i_dR][i_cent]->GetYaxis()->SetRangeUser(low_range, hi_range);

				h_TM_MB_2D[i_dR][i_cent]->GetZaxis()->SetRangeUser(0,2);
				h_FS_MB_2D[i_dR][i_cent]->GetZaxis()->SetRangeUser(0,2);
				h_FNS_MB_2D[i_dR][i_cent]->GetZaxis()->SetRangeUser(0,2);


				c_TM_MB_2D->cd(i_cent+1);
				h_TM_MB_2D[i_dR][i_cent]->Draw("colz text");

				ltx->SetTextAlign(12);
				ltx->DrawLatexNDC(0.19,0.98,dr_label.c_str());
				ltx->SetTextAlign(32);
				ltx->DrawLatexNDC(0.90,0.98,num_to_cent(31,i_cent).c_str());
				gPad->SetLogx();
				gPad->SetLogy();

				c_FS_MB_2D->cd(i_cent+1);
				h_FS_MB_2D[i_dR][i_cent]->Draw("colz text");
				ltx->SetTextAlign(12);
				ltx->DrawLatexNDC(0.19,0.98,dr_label.c_str());
				ltx->SetTextAlign(32);
				ltx->DrawLatexNDC(0.90,0.98,num_to_cent(31,i_cent).c_str());
				gPad->SetLogx();
				gPad->SetLogy(0);

				c_FNS_MB_2D->cd(i_cent+1);
				h_FNS_MB_2D[i_dR][i_cent]->Draw("colz text");
				ltx->SetTextAlign(12);
				ltx->DrawLatexNDC(0.19,0.98,dr_label.c_str());
				ltx->SetTextAlign(32);
				ltx->DrawLatexNDC(0.90,0.98,num_to_cent(31,i_cent).c_str());
				gPad->SetLogx();
				gPad->SetLogy(0);
			}

			if (i_dR == 0) name = "(";
			else if (i_dR == 12) name = ")";
			else name = "";
			c_TM_MB_2D->Print(Form("UE_TM_MB_2D.pdf%s", name.c_str()), Form("Title: dR%i - %s", i_dR, dr_label.c_str()));
//			c_FS_MB_2D->Print(Form("UE_FS_MB_2D.pdf%s", name.c_str()), Form("Title: dR%i - %s", i_dR, dr_label.c_str()));
//			c_FNS_MB_2D->Print(Form("UE_FNS_MB_2D.pdf%s", name.c_str()), Form("Title: dR%i - %s", i_dR, dr_label.c_str()));
		}
	}

	/*
	{
		//projecting over track pT
		cout << "projecting over track pT" << endl;

		TCanvas *c_TM_MB_1D = new TCanvas("c_TM_MB_1D","c_TM_MB_1D",900,600);
		TCanvas *c_FS_MB_1D = new TCanvas("c_FS_MB_1D","c_FS_MB_1D",900,600);
		TCanvas *c_FNS_MB_1D = new TCanvas("c_FNS_MB_1D","c_FNS_MB_1D",900,600);

		TLegend *legend_TM_MB = new TLegend(0.19, 0.70, 0.40, 0.92, "","brNDC");
		legend_TM_MB->SetTextFont(43);
		legend_TM_MB->SetBorderSize(0);
		legend_TM_MB->SetTextSize(10);

		TLegend *legend_FS_MB = new TLegend(0.19, 0.70, 0.40, 0.92, "","brNDC");
		legend_FS_MB->SetTextFont(43);
		legend_FS_MB->SetBorderSize(0);
		legend_FS_MB->SetTextSize(10);

		TLegend *legend_FNS_MB = new TLegend(0.19, 0.70, 0.40, 0.92, "","brNDC");
		legend_FNS_MB->SetTextFont(43);
		legend_FNS_MB->SetBorderSize(0);
		legend_FNS_MB->SetTextSize(10);


		for (int i_dR = 0; i_dR < N_dR; i_dR++)
		{
			string dr_label = Form("%1.2f < r < %1.2f", dR_binning->GetBinLowEdge(i_dR+1), dR_binning->GetBinUpEdge(i_dR+1));

			c_TM_MB_1D->Clear();
			c_TM_MB_1D->Divide(3,2);

			c_FS_MB_1D->Clear();
			c_FS_MB_1D->Divide(3,2);

			c_FNS_MB_1D->Clear();
			c_FNS_MB_1D->Divide(3,2);

			for (int i_cent = 0; i_cent < 6; i_cent++)
			{

				int jet_itr = 0;

				for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
				{
					string jet_label = Form("%1.0f < p_{T}^{Jet} < %1.0f", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));


					SetHStyle_smallify(h_TM_MB_1D[i_jet][i_dR][i_cent], jet_itr, 1);
					SetHStyle_smallify(h_FS_MB_1D[i_jet][i_dR][i_cent], jet_itr, 1);
					SetHStyle_smallify(h_FNS_MB_1D[i_jet][i_dR][i_cent], jet_itr, 1);


					if (i_cent == 0 && i_dR == 0) legend_TM_MB->AddEntry(h_TM_MB_1D[i_jet][i_dR][i_cent],jet_label.c_str(),"lp");
					if (i_cent == 0 && i_dR == 0) legend_FS_MB->AddEntry(h_FS_MB_1D[i_jet][i_dR][i_cent],jet_label.c_str(),"lp");
					if (i_cent == 0 && i_dR == 0) legend_FNS_MB->AddEntry(h_FNS_MB_1D[i_jet][i_dR][i_cent],jet_label.c_str(),"lp");


					h_TM_MB_1D[i_jet][i_dR][i_cent]->GetXaxis()->SetRangeUser(1,10);
					h_FS_MB_1D[i_jet][i_dR][i_cent]->GetXaxis()->SetRangeUser(1,10);
					h_FNS_MB_1D[i_jet][i_dR][i_cent]->GetXaxis()->SetRangeUser(1,10);

					double low_range = 0.5;
					double hi_range = 1.5;
					h_TM_MB_1D[i_jet][i_dR][i_cent]->GetYaxis()->SetRangeUser(low_range, hi_range);
					h_FS_MB_1D[i_jet][i_dR][i_cent]->GetYaxis()->SetRangeUser(low_range, hi_range);
					h_FNS_MB_1D[i_jet][i_dR][i_cent]->GetYaxis()->SetRangeUser(low_range, hi_range);


					c_TM_MB_1D->cd(i_cent+1);
					if (i_jet == 7) h_TM_MB_1D[i_jet][i_dR][i_cent]->Draw("");
					else h_TM_MB_1D[i_jet][i_dR][i_cent]->Draw("same");
					gPad->SetLogx();

					c_FS_MB_1D->cd(i_cent+1);
					if (i_jet == 7) h_FS_MB_1D[i_jet][i_dR][i_cent]->Draw("");
					else h_FS_MB_1D[i_jet][i_dR][i_cent]->Draw("same");
					gPad->SetLogx();

					c_FNS_MB_1D->cd(i_cent+1);
					if (i_jet == 7) h_FNS_MB_1D[i_jet][i_dR][i_cent]->Draw("");
					else h_FNS_MB_1D[i_jet][i_dR][i_cent]->Draw("same");
					gPad->SetLogx();

					jet_itr++;
				}

				c_TM_MB_1D->cd(i_cent+1);
				ltx->SetTextAlign(12);
				ltx->DrawLatexNDC(0.19,0.98,dr_label.c_str());
				ltx->SetTextAlign(32);
				ltx->DrawLatexNDC(0.90,0.98,num_to_cent(31,i_cent).c_str());
				legend_TM_MB->Draw();

				c_FS_MB_1D->cd(i_cent+1);
				ltx->SetTextAlign(12);
				ltx->DrawLatexNDC(0.19,0.98,dr_label.c_str());
				ltx->SetTextAlign(32);
				ltx->DrawLatexNDC(0.90,0.98,num_to_cent(31,i_cent).c_str());
				legend_FS_MB->Draw();

				c_FNS_MB_1D->cd(i_cent+1);
				ltx->SetTextAlign(12);
				ltx->DrawLatexNDC(0.19,0.98,dr_label.c_str());
				ltx->SetTextAlign(32);
				ltx->DrawLatexNDC(0.90,0.98,num_to_cent(31,i_cent).c_str());
				legend_FNS_MB->Draw();
			}
			if (i_dR == 0) name = "(";
			else if (i_dR == 12) name = ")";
			else name = "";
			c_TM_MB_1D->Print(Form("UE_TM_MB_1D.pdf%s", name.c_str()), Form("Title: dR%i - %s", i_dR, dr_label.c_str()));
			c_FS_MB_1D->Print(Form("UE_FS_MB_1D.pdf%s", name.c_str()), Form("Title: dR%i - %s", i_dR, dr_label.c_str()));
			c_FNS_MB_1D->Print(Form("UE_FNS_MB_1D.pdf%s", name.c_str()), Form("Title: dR%i - %s", i_dR, dr_label.c_str()));
		}
	}
*/
	/*
	{
		//as function of r
		cout << "Function of R" << endl;
		TCanvas *c_MB_r_1D = new TCanvas("c_MB_r_1D","c_MB_r_1D",900,600);
		TCanvas *c_TM_r_1D = new TCanvas("c_TM_r_1D","c_TM_r_1D",900,600);
		TCanvas *c_FS_r_1D = new TCanvas("c_FS_r_1D","c_FS_r_1D",900,600);
		TCanvas *c_FNS_r_1D = new TCanvas("c_FNS_r_1D","c_FNS_r_1D",900,600);

		TLegend *legend_MB_r = new TLegend(0.19, 0.70, 0.40, 0.92, "","brNDC");
		legend_MB_r->SetTextFont(43);
		legend_MB_r->SetBorderSize(0);
		legend_MB_r->SetTextSize(10);

		TLegend *legend_TM_r = new TLegend(0.19, 0.70, 0.40, 0.92, "","brNDC");
		legend_TM_r->SetTextFont(43);
		legend_TM_r->SetBorderSize(0);
		legend_TM_r->SetTextSize(10);

		TLegend *legend_FS_r = new TLegend(0.19, 0.70, 0.40, 0.92, "","brNDC");
		legend_FS_r->SetTextFont(43);
		legend_FS_r->SetBorderSize(0);
		legend_FS_r->SetTextSize(10);

		TLegend *legend_FNS_r = new TLegend(0.19, 0.70, 0.40, 0.92, "","brNDC");
		legend_FNS_r->SetTextFont(43);
		legend_FNS_r->SetBorderSize(0);
		legend_FNS_r->SetTextSize(10);

		for (int i_trk = 0; i_trk < N_trkpt; i_trk++)
		{
			if (i_trk < 2 || i_trk > 6) continue;

			string trk_label = Form("%1.2f < p_{T}^{Trk} < %1.2f", trkpT_binning->GetBinLowEdge(i_trk+1), trkpT_binning->GetBinUpEdge(i_trk+1));

			c_MB_r_1D->Clear();
			c_MB_r_1D->Divide(3,2);

			c_TM_r_1D->Clear();
			c_TM_r_1D->Divide(3,2);

			c_FS_r_1D->Clear();
			c_FS_r_1D->Divide(3,2);

			c_FNS_r_1D->Clear();
			c_FNS_r_1D->Divide(3,2);

			for (int i_cent = 0; i_cent < 6; i_cent++)
			{
				int jet_itr = 0;

				for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
				{
					string jet_label = Form("%1.0f < p_{T}^{Jet} < %1.0f", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));

					SetHStyle_smallify(h_MB_r_1D[i_jet][i_trk][i_cent], jet_itr, 1);
					SetHStyle_smallify(h_TM_r_1D[i_jet][i_trk][i_cent], jet_itr, 1);
					SetHStyle_smallify(h_FS_r_1D[i_jet][i_trk][i_cent], jet_itr, 1);
					SetHStyle_smallify(h_FNS_r_1D[i_jet][i_trk][i_cent], jet_itr, 1);

					if (i_cent == 0 && i_trk == 2) legend_MB_r->AddEntry(h_MB_r_1D[i_jet][i_trk][i_cent],jet_label.c_str(),"lp");
					if (i_cent == 0 && i_trk == 2) legend_TM_r->AddEntry(h_TM_r_1D[i_jet][i_trk][i_cent],jet_label.c_str(),"lp");
					if (i_cent == 0 && i_trk == 2) legend_FS_r->AddEntry(h_FS_r_1D[i_jet][i_trk][i_cent],jet_label.c_str(),"lp");
					if (i_cent == 0 && i_trk == 2) legend_FNS_r->AddEntry(h_FNS_r_1D[i_jet][i_trk][i_cent],jet_label.c_str(),"lp");

					double avg, low_range, hi_range;

					avg = (h_MB_r_1D[i_jet][i_trk][i_cent]->GetMaximum() + h_MB_r_1D[i_jet][i_trk][i_cent]->GetMinimum())/2;
					low_range = avg * 0.90; hi_range = avg * 1.10;

					h_MB_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetRangeUser(low_range, hi_range);
					h_TM_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetRangeUser(low_range, hi_range);
					h_FS_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetRangeUser(low_range, hi_range);
					h_FNS_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetRangeUser(low_range, hi_range);

					h_MB_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetNdivisions(504);
					h_TM_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetNdivisions(504);
					h_FS_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetNdivisions(504);
					h_FNS_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetNdivisions(504);

					c_MB_r_1D->cd(i_cent+1);
					if (jet_itr == 0) h_MB_r_1D[i_jet][i_trk][i_cent]->Draw("");
					else h_MB_r_1D[i_jet][i_trk][i_cent]->Draw("same");
//					gPad->SetLogy();

					c_TM_r_1D->cd(i_cent+1);
					if (jet_itr == 0) h_TM_r_1D[i_jet][i_trk][i_cent]->Draw("");
					else h_TM_r_1D[i_jet][i_trk][i_cent]->Draw("same");
//					gPad->SetLogy();

					c_FS_r_1D->cd(i_cent+1);
					if (jet_itr == 0) h_FS_r_1D[i_jet][i_trk][i_cent]->Draw("");
					else h_FS_r_1D[i_jet][i_trk][i_cent]->Draw("same");
//					gPad->SetLogy();

					c_FNS_r_1D->cd(i_cent+1);
					if (jet_itr == 0) h_FNS_r_1D[i_jet][i_trk][i_cent]->Draw("");
					else h_FNS_r_1D[i_jet][i_trk][i_cent]->Draw("same");
//					gPad->SetLogy();

					jet_itr++;
				}

				c_MB_r_1D->cd(i_cent+1);
				ltx->SetTextAlign(12);
				ltx->DrawLatexNDC(0.19,0.98,trk_label.c_str());
				ltx->SetTextAlign(32);
				ltx->DrawLatexNDC(0.90,0.98,num_to_cent(31,i_cent).c_str());
				legend_MB_r->Draw();

				c_TM_r_1D->cd(i_cent+1);
				ltx->SetTextAlign(12);
				ltx->DrawLatexNDC(0.19,0.98,trk_label.c_str());
				ltx->SetTextAlign(32);
				ltx->DrawLatexNDC(0.90,0.98,num_to_cent(31,i_cent).c_str());
				legend_TM_r->Draw();

				c_FS_r_1D->cd(i_cent+1);
				ltx->SetTextAlign(12);
				ltx->DrawLatexNDC(0.19,0.98,trk_label.c_str());
				ltx->SetTextAlign(32);
				ltx->DrawLatexNDC(0.90,0.98,num_to_cent(31,i_cent).c_str());
				legend_FS_r->Draw();

				c_FNS_r_1D->cd(i_cent+1);
				ltx->SetTextAlign(12);
				ltx->DrawLatexNDC(0.19,0.98,trk_label.c_str());
				ltx->SetTextAlign(32);
				ltx->DrawLatexNDC(0.90,0.98,num_to_cent(31,i_cent).c_str());
				legend_FNS_r->Draw();
			}
			if (i_trk == 2) name = "(";
			else if (i_trk == 6) name = ")";
			else name = "";
			c_MB_r_1D->Print(Form("UE_MB_r_1D.pdf%s", name.c_str()), Form("Title: trk%i - %s", i_trk, trk_label.c_str()));
			c_TM_r_1D->Print(Form("UE_TM_r_1D.pdf%s", name.c_str()), Form("Title: trk%i - %s", i_trk, trk_label.c_str()));
			c_FS_r_1D->Print(Form("UE_FS_r_1D.pdf%s", name.c_str()), Form("Title: trk%i - %s", i_trk, trk_label.c_str()));
			c_FNS_r_1D->Print(Form("UE_FNS_r_1D.pdf%s", name.c_str()), Form("Title: trk%i - %s", i_trk, trk_label.c_str()));
		}
	}
*/

	{
		//as function of r
		cout << "Function of R" << endl;
		TCanvas *c_x = new TCanvas("c_x","c_x",900,600);

		TLegend *legend_x = new TLegend(0.20, 0.10, 0.50, 0.40, "","brNDC");
		legend_x->SetTextFont(43);
		legend_x->SetBorderSize(0);
		legend_x->SetTextSize(10);


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
					SetHStyle_smallify(h_MB_tj_r_1D[i_jet][i_trk][i_cent], 4, 1);
					SetHStyle_smallify(h_MB_data_r_1D[i_jet][i_trk][i_cent], 5, 1);

					SetHStyle_smallify(h_TM_MB_r_1D[i_jet][i_trk][i_cent], 1, 1);
					SetHStyle_smallify(h_FS_MB_r_1D[i_jet][i_trk][i_cent], 2, 1);
					SetHStyle_smallify(h_FNS_MB_r_1D[i_jet][i_trk][i_cent], 3, 1);
					SetHStyle_smallify(h_MB_tj_MB_r_1D[i_jet][i_trk][i_cent], 4, 1);
					SetHStyle_smallify(h_MB_data_MB_r_1D[i_jet][i_trk][i_cent], 5, 1);

//					h_MB_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetRangeUser(1E-4,1E2);

					if (jet_itr == 0 && i_trk == 2 && i_cent == 0)
					{
						legend_x->AddEntry(h_MB_r_1D[i_jet][i_trk][i_cent],"MB","lp");
//						legend_x->AddEntry(h_MB_data_r_1D[i_jet][i_trk][i_cent],"MB_data","lp");
//						legend_x->AddEntry(h_MB_tj_r_1D[i_jet][i_trk][i_cent],"MB_{TJ}","lp");
						legend_x->AddEntry(h_TM_r_1D[i_jet][i_trk][i_cent],"TM","lp");
						legend_x->AddEntry(h_FS_r_1D[i_jet][i_trk][i_cent],"FS","lp");
						legend_x->AddEntry(h_FNS_r_1D[i_jet][i_trk][i_cent],"FNS","lp");
					}

					double avg, low_range, hi_range;
					avg = (h_MB_r_1D[i_jet][i_trk][i_cent]->GetMaximum() + h_MB_r_1D[i_jet][i_trk][i_cent]->GetMinimum())/2;
					low_range = avg * 0.50; hi_range = avg * 1.50;

//					h_MB_r_1D[i_jet][i_trk][i_cent]->GetXaxis()->SetRangeUser(0, 0.6);
					h_MB_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetRangeUser(low_range, hi_range);
					h_MB_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetNdivisions(504);


					c_x->cd(i_cent+1);
					gPad->Divide(1,2);

					c_x->cd(i_cent+1)->cd(1);
					gPad->SetPad(0,0.40,0.99,0.99);
					gPad->SetTopMargin(0.05);
					gPad->SetBottomMargin(0);
					gPad->SetRightMargin(0);
					h_MB_r_1D[i_jet][i_trk][i_cent]->Draw();
//					h_MB_data_r_1D[i_jet][i_trk][i_cent]->Draw("same");
					h_TM_r_1D[i_jet][i_trk][i_cent]->Draw("same");
					h_FS_r_1D[i_jet][i_trk][i_cent]->Draw("same");
					h_FNS_r_1D[i_jet][i_trk][i_cent]->Draw("same");
					legend_x->Draw();



					c_x->cd(i_cent+1)->cd(2);
					gPad->SetPad(0,0.0,0.99,0.40);
					gPad->SetTopMargin(0);
					gPad->SetBottomMargin(0.30);
					gPad->SetRightMargin(0);
					h_TM_MB_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetRangeUser(0.9, 1.1);
					h_TM_MB_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetNdivisions(504);
					h_MB_data_MB_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetRangeUser(0.99, 1.01);
					h_MB_data_MB_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetNdivisions(504);
					h_TM_MB_r_1D[i_jet][i_trk][i_cent]->Draw();
//					h_MB_data_MB_r_1D[i_jet][i_trk][i_cent]->Draw("");
					h_FS_MB_r_1D[i_jet][i_trk][i_cent]->Draw("same");
					h_FNS_MB_r_1D[i_jet][i_trk][i_cent]->Draw("same");
					line->DrawLine(0, 1, 0.6, 1);


					c_x->cd(i_cent+1);
					ltx->SetTextAlign(32);
					ltx->DrawLatexNDC(0.94,0.94,num_to_cent(31,i_cent).c_str());
					ltx->DrawLatexNDC(0.94,0.88,jet_label.c_str());
					ltx->DrawLatexNDC(0.94,0.82,trk_label.c_str());

				}

				jet_itr++;

				if (i_trk == 2 && i_jet == jet_pt_start) name = "(";
				else if (i_trk == 6 && i_jet == jet_pt_end - 1) name = ")";
				else name = "";
				c_x->Print(Form("UE_x.pdf%s", name.c_str()), Form("Title: trk%i_jet%i", i_trk, i_jet));

			}

		}
	}

	{
		//just the factors as function of r
		cout << "posres as Function of R" << endl;
		TCanvas *c_pos_res = new TCanvas("c_pos_res","c_pos_res",900,600);

		TLegend *legend_pos_res = new TLegend(0.20, 0.10, 0.50, 0.40, "","brNDC");
		legend_pos_res->SetTextFont(43);
		legend_pos_res->SetBorderSize(0);
		legend_pos_res->SetTextSize(10);



		int jet_itr = 0;
		for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
		{
			string jet_label = Form("%1.0f < p_{T}^{Jet} < %1.0f", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));
			c_pos_res->Clear();
			c_pos_res->Divide(3,2);

			for (int i_cent = 0; i_cent < 6; i_cent++)
			{
				int trk_itr = 0;
				for (int i_trk = 0; i_trk < N_trkpt; i_trk++)
				{
					if (i_trk < 2 || i_trk > 6) continue;

					string trk_label = Form("%1.2f < p_{T}^{Trk} < %1.2f", trkpT_binning->GetBinLowEdge(i_trk+1), trkpT_binning->GetBinUpEdge(i_trk+1));

					SetHStyle_smallify(h_TM_MB_r_1D[i_jet][i_trk][i_cent], trk_itr, 1);

					if (jet_itr == 0 && i_cent == 0)
					{
						legend_pos_res->AddEntry(h_TM_MB_r_1D[i_jet][i_trk][i_cent],trk_label.c_str(),"lp");
					}

//					double avg, low_range, hi_range;
//					avg = (h_MB_r_1D[i_jet][i_trk][i_cent]->GetMaximum() + h_MB_r_1D[i_jet][i_trk][i_cent]->GetMinimum())/2;
//					low_range = avg * 0.50; hi_range = avg * 1.50;

//					h_TM_MB_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetRangeUser(low_range, hi_range);
					h_TM_MB_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetRangeUser(0, 2);
					h_TM_MB_r_1D[i_jet][i_trk][i_cent]->GetYaxis()->SetNdivisions(504);


					c_pos_res->cd(i_cent+1);
					if (trk_itr == 0) h_TM_MB_r_1D[i_jet][i_trk][i_cent]->Draw();
					else h_TM_MB_r_1D[i_jet][i_trk][i_cent]->Draw("same");
					line->DrawLine(0, 1, 1.2, 1);


					trk_itr++;
				}

				c_pos_res->cd(i_cent+1);
				ltx->SetTextAlign(32);
				ltx->DrawLatexNDC(0.94,0.94,num_to_cent(31,i_cent).c_str());
				ltx->DrawLatexNDC(0.94,0.88,jet_label.c_str());
				legend_pos_res->Draw();

			}

			if (i_jet == jet_pt_start) name = "(";
			else if (i_jet == jet_pt_end - 1) name = ")";
			else name = "";
			cout << name << endl;
			c_pos_res->Print(Form("UE_pos_res.pdf%s", name.c_str()), Form("Title: jet%i", i_jet));
			jet_itr++;

		}
	}

}
