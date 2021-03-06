#include <iostream>
using std::cout;
using std::endl;
#include <vector>

#include "TRandom.h"
#include <TEnv.h>
#include "/Applications/root_v6.10.08/RooUnfold/src/RooUnfoldResponse.h"
#include "/Applications/root_v6.10.08/RooUnfold/src/RooUnfold.h"
#include "/Applications/root_v6.10.08/RooUnfold/src/RooUnfoldBayes.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TSystem.h"
#include "TAxis.h"
#include "TError.h"
#include "TCanvas.h"
#include "TLine.h"
#include "../../functions/global_variables.h"
#include "TMath.h"

int main(int argc, char ** argv)
{

	cout << "######### RUNNING UNFOLDING #########" << endl;
	std::string config_file = "ff_config.cfg";
	if (argc == 2) config_file = argv[1];
	gErrorIgnoreLevel = 3001;
	std::string name;

	//	##############	Reading config	##############"
	TEnv *m_config = new TEnv();
	m_config->ReadFile(config_file.c_str(), EEnvLevel(1));

	std::string dataset_type = "PbPb"; dataset_type = m_config->GetValue("dataset_type", dataset_type.c_str());
	int isMC = 1; isMC = m_config->GetValue("isMC", isMC);
	int sys_mode = -1; sys_mode = m_config->GetValue("sys_mode", sys_mode);

	int centrality_scheme = 31; centrality_scheme = m_config->GetValue("centrality_scheme", centrality_scheme);
	int n_unfold = 4; n_unfold = m_config->GetValue("n_unfold", n_unfold);
	int verbose = 0; verbose = m_config->GetValue("verbose", verbose);

	std::string did = "data";
	if (isMC) did = "MC";

	int apply_MC_nonClos = 0;
	if (sys_mode == 4) apply_MC_nonClos = 1;

	int apply_fake_uncert = 0;
	if (sys_mode == 45) apply_fake_uncert = 1;

	int apply_UE_mapStatSys = 0;
	if (sys_mode == 42) apply_UE_mapStatSys = 1;

	int apply_UE_RunDep = 0;
	if (sys_mode == 43) apply_UE_RunDep = 1;

	int apply_UE_ConeMethod = 0;
	if (sys_mode == 44) apply_UE_ConeMethod = 1;

	if (verbose) m_config->Print();
	//	##############	Config done	##############"



	cout << "isMC: " << isMC << endl;
	cout << "sys_mode: " << sys_mode << endl;
	cout << "centrality_scheme: " << centrality_scheme << endl;
	cout << "n_unfold: " << n_unfold << endl;
	cout << "verbose: " << verbose << endl;




	std::string sys_path = "";
	if (sys_mode == 0) sys_path = Form("nominal");
	else if (sys_mode > 50 || sys_mode < 0) sys_path = Form("c%i", sys_mode);
	else sys_path = Form("sys%i", sys_mode);

	TFile *f_mc = new TFile(Form("../raw_results/%s/FF_MC_out_histo_%s_5p02_r001.root", sys_path.c_str(), dataset_type.c_str()));
	TFile *f_data = new TFile(Form("../raw_results/%s/FF_%s_out_histo_%s_5p02_r001.root",sys_path.c_str(), did.c_str(), dataset_type.c_str()));

	cout << "Using files:" << endl;
	cout << f_mc->GetName() << endl;
	cout << f_data->GetName() << endl;	

//	if (sys_mode >= 0) sys_path = Form("_%s", sys_path.c_str());
	TFile *dr_factors = new TFile(Form("output_pdf_%s/root/posCorr_factors_%s.root", sys_path.c_str(), dataset_type.c_str()));
	TFile *UE_uncert, *fake_uncert;
	if (apply_UE_RunDep) UE_uncert = new TFile(Form("UE_RunDependentSys.root"));
	if (apply_UE_mapStatSys) UE_uncert = new TFile(Form("UE_MapSystematic.root"));
	if (apply_fake_uncert) fake_uncert = new TFile(Form("../raw_results/nominal/FF_MC_out_histo_pp_5p02_r001.root"));
	TFile *f_output = new TFile(Form("output_pdf_%s/root/raw_unfolded_%s_%s.root", sys_path.c_str(), did.c_str(), dataset_type.c_str()),"recreate");

//	TFile *f_output = new TFile(Form("iteration_test/raw_unfolded_iter%i_%s_%s.root",n_unfold, did.c_str(), dataset_type.c_str()),"recreate");
//	cout << f_output->GetName() << endl;
	int N_Y = 5;

	TAxis* dR_binning = (TAxis*)((TH3*)f_mc->Get("h_dR_change_jetpt0_cent0"))->GetXaxis();
	TAxis* jetpT_binning = (TAxis*)((TH3*)f_mc->Get("ChPS_raw_0_dR0_cent0"))->GetYaxis();
	TAxis* trkpT_binning = (TAxis*)((TH3*)f_mc->Get("ChPS_raw_0_dR0_cent0"))->GetXaxis();

	int N_dR = dR_binning->GetNbins();
	int N_jetpt = jetpT_binning->GetNbins();
	int N_trkpt = trkpT_binning->GetNbins();

	f_output->cd();
	dR_binning->Write("dR_binning");
	jetpT_binning->Write("jetpT_binning");
	trkpT_binning->Write("trkpT_binning");

	TCanvas *c1 = new TCanvas("c","c",1200,600);
	delete c1; //i dont know why this is necessary. crashes without creating and deleting TCanvas

	for (int i_cent = 0; i_cent < n_cent_cuts; i_cent++)
	{
		cout << Form("Starting cent%i...", i_cent) << endl;

		TH1 *h_reco_jet, *h_truth_jet, *h_reco_unfolded;
		if (dataset_type == "PbPb" && i_cent == 6) continue;
		if (dataset_type == "pp" && i_cent < 6) continue;

		for (int i_y = 0; i_y < N_Y; i_y++)
		{
			name = Form("h_reco_jet_spectrum_y%i_cent%i",i_y, i_cent);
			TH1* h_reco_jet_y_c = (TH1*)((TH1*)f_data->Get(name.c_str()))->Clone(Form("reco_jet_y%i_c%i",i_y, i_cent));
			h_reco_jet_y_c->Sumw2();

			name = Form("h_true_jet_spectrum_y%i_cent%i",i_y, i_cent);
			TH1* h_true_jet_y_c = (TH1*)((TH1*)f_mc->Get(name.c_str()))->Clone(Form("true_jet_y%i_c%i",i_y, i_cent));
			h_true_jet_y_c->Sumw2();

			TH1* h_unfolded_jet_y_c = (TH1*)h_true_jet_y_c->Clone(Form("h_raw_unfolded_y%i_c%i", i_y, i_cent));
			h_unfolded_jet_y_c->Reset();

			//unfold only if using proper reco ChPS (not if using *reco matched* truth)
			name = Form("h_reco_jet_spectrum_y%i_cent%i_h_true_jet_spectrum_y%i_cent%i",i_y, i_cent, i_y, i_cent);
			RooUnfoldResponse* r_response = (RooUnfoldResponse*)f_mc->Get(name.c_str());

			TH2* h_response_matrix = r_response->Hresponse();
			name = Form("h_response_matrix_jet_y%i_c%i", i_y, i_cent);
			h_response_matrix->SetTitle(name.c_str());
			h_response_matrix->Write(name.c_str());

			RooUnfoldBayes unfold(r_response, h_reco_jet_y_c, 4);
			unfold.SetVerbose(0);
			if (h_reco_jet_y_c->GetEntries() != 0)
			{
				h_unfolded_jet_y_c = (TH1D*)unfold.Hreco();
				name = Form("unfolded_y%i_c%i", i_y, i_cent);
				h_unfolded_jet_y_c->SetName(name.c_str());
				h_unfolded_jet_y_c->Sumw2();
			}

			delete r_response;

			if (i_y == 4)
			{
				h_reco_jet = (TH1*)h_reco_jet_y_c->Clone(Form("h_inc_reco_jet_y%i_c%i",i_y, i_cent));
				h_reco_jet->Sumw2();

				h_truth_jet = (TH1*)h_true_jet_y_c->Clone(Form("h_inc_truth_matched_jet_y%i_c%i",i_y, i_cent));
				h_truth_jet->Sumw2();

				h_reco_unfolded = (TH1*)h_unfolded_jet_y_c->Clone(Form("h_inc_reco_unf_jet_y%i_c%i",i_y, i_cent));
				h_reco_unfolded->Sumw2();
			}

			f_output->cd();
			name = Form("h_reco_jet_y%i_c%i", i_y, i_cent);
			h_reco_jet_y_c->SetTitle(name.c_str());
			h_reco_jet_y_c->Write(name.c_str());
			delete h_reco_jet_y_c;

			name = Form("h_true_jet_y%i_c%i", i_y, i_cent);
			h_true_jet_y_c->SetTitle(name.c_str());
			h_true_jet_y_c->Write(name.c_str());
			delete h_true_jet_y_c;

			name = Form("h_unfolded_jet_y%i_c%i", i_y, i_cent);
			h_unfolded_jet_y_c->SetTitle(name.c_str());
			h_unfolded_jet_y_c->Write(name.c_str());
			delete h_unfolded_jet_y_c;

		}


		TH2* h_raw_injet;
		TH2* h_raw_subtr_injet;
		TH2* h_raw_subtr_unf_injet;
		TH2* h_raw_subtr_unf_bbb_injet;

		TH2* h_UE_injet;

		TH2* h_truth_injet;


		for (int i_dR = 0; i_dR < N_dR; i_dR++)
		{
			TH2* h_raw = (TH2*)f_data->Get(Form("ChPS_raw_0_dR%i_cent%i", i_dR, i_cent));
			h_raw->Sumw2();

			TH2* h_truth = (TH2*)f_mc->Get(Form("ChPS_truth_dR%i_cent%i", i_dR, i_cent));
			h_truth->Sumw2();

			//setup UE/fakes
			TH2* h_UE_MB = (TH2*)f_data->Get(Form("ChPS_MB_UE_dR%i_cent%i",i_dR, i_cent));;
			TH2* h_UE_TM = (TH2*)f_mc->Get(Form("ChPS_TM_UE_dR%i_cent%i",i_dR, i_cent));;

			TH1* h_reco_jet_spect_mc = (TH1*)((TH1*)f_mc->Get(Form("h_reco_jet_spectrum_y4_cent%i", i_cent)))->Clone(Form("norm_reco_mc_jet_y4_cent%i", i_cent));
			h_reco_jet_spect_mc->SetName(Form("norm_reco_mc_jet_y4_cent%i", i_cent));

			TH1* h_reco_jet_spect_data = (TH1*)((TH1*)f_data->Get(Form("h_reco_jet_spectrum_y4_cent%i", i_cent)))->Clone(Form("norm_reco_data_jet_y4_cent%i", i_cent));
			h_reco_jet_spect_data->SetName(Form("norm_reco_data_jet_y4_cent%i", i_cent));


			TH2* h_UE_corrected_cone, *h_map_stat, *h_fake_rate;
			TH1* h_pp_MC_jets;
			if (dataset_type == "PbPb" && apply_UE_ConeMethod)
			{
				TH2* h_UE_cone_MC = (TH2*)f_mc->Get(Form("ChPS_cone_UE_dR%i_cent%i",i_dR, i_cent));;
				h_UE_cone_MC->SetName((Form("ChPS_cone_UE_MC_dR%i_cent%i",i_dR, i_cent)));

				TH2* h_UE_cone_data = (TH2*)f_data->Get(Form("ChPS_cone_UE_dR%i_cent%i",i_dR, i_cent));;
				h_UE_cone_data->SetName(Form("ChPS_cone_UE_data_dR%i_cent%i",i_dR, i_cent));

				h_UE_corrected_cone = (TH2*)h_UE_cone_data->Clone(Form("ChPS_cone_UE_data_corrected_dR%i_cent%i",i_dR, i_cent));
				h_UE_corrected_cone->Multiply(h_UE_TM);
				h_UE_corrected_cone->Divide(h_UE_cone_MC);
				delete h_UE_cone_MC;
				delete h_UE_cone_data;
			}

			if (dataset_type == "PbPb" && apply_UE_mapStatSys)
			{
				h_map_stat = (TH2*)UE_uncert->Get(Form("ChPS_MB_UE_dR%i_cent%i_statSys", i_dR, i_cent));
			}

			if (apply_fake_uncert)
			{
				h_fake_rate = (TH2*)(fake_uncert->Get(Form("ChPS_TM_UE_dR%i_cent6", i_dR))->Clone(Form("fake_dR%i_cent%i", i_dR, i_cent)));
				h_fake_rate->SetName(Form("fake_dR%i_cent%i", i_dR, i_cent));
				h_pp_MC_jets = (TH1*)fake_uncert->Get(Form("h_reco_jet_spectrum_y4_cent6"))->Clone(Form("fake_norm_%i",i_cent));
				h_pp_MC_jets->SetName(Form("fake_norm_%i",i_cent));
			}

			TH2* h_final_UE = (TH2*)h_UE_MB->Clone(Form("final_UE_%i_cent%i",i_dR, i_cent));;
			h_final_UE->Reset();
			h_final_UE->Sumw2();

			TH2* h_final_fake = (TH2*)h_UE_MB->Clone(Form("final_fake_%i_cent%i",i_dR, i_cent));;
			h_final_fake->Reset();
			h_final_fake->Sumw2();

			TH2* h_raw_subtr = (TH2*)h_raw->Clone(Form("h_raw_subtr_dR%i_c%i", i_dR, i_cent));
			h_raw_subtr->Reset();
			h_raw_subtr->Sumw2();

			//USING i_jet_bin AS BIN NUMBER
			for (int i_jet_bin = 1; i_jet_bin <= N_jetpt; i_jet_bin++)
			{
				if (i_jet_bin < h_reco_jet_spect_mc->FindBin(100) || i_jet_bin > h_reco_jet_spect_mc->FindBin(390)) continue;
				double n_jets_mc = h_reco_jet_spect_mc->GetBinContent(i_jet_bin);
				double n_jets_data = h_reco_jet_spect_data->GetBinContent(i_jet_bin);
				double n_jets_pp_mc = 1;
				if (apply_fake_uncert) n_jets_pp_mc = h_pp_MC_jets->GetBinContent(i_jet_bin);
				//cone method uses same number of jets as n_jets_data

				if (n_jets_mc == 0) continue;

				//if pp. Subtract only fakes (UE_TM|MC) for all pT. Need to normalize n_jets_data / n_jets_mc. This correction is 1 if running on MC since UE_TM|MC has same # of jets as in MC
				//if PbPb
				//	if < 10 GeV: Subtract (UE_MB|data). No Need to normalize n_jets_w / n_jets_unw.
				//	if >= 10 GeV: Subtract (UE_TM|MC). Need to normalize n_jets_data / n_jets_mv. This correction is 1 if running on MC since UE_TM|MC has same # of jets as in MC

				//USING i_trk_bin AS BIN NUMBER
				for (int i_trk_bin = 1; i_trk_bin <= N_trkpt; i_trk_bin++)
				{

					double UE = 0, UE_err = 0;
					double fake = 0, fake_err = 0;
					double raw = 0, raw_err = 0;
					double subtr = 0, subtr_err = 0;

					raw = h_raw->GetBinContent(i_trk_bin, i_jet_bin);
					raw_err = h_raw->GetBinError(i_trk_bin, i_jet_bin);
					if (raw == 0)
					{
						h_raw_subtr->SetBinContent(i_trk_bin, i_jet_bin, 0);
						h_raw_subtr->SetBinError(i_trk_bin, i_jet_bin, 0);
						continue;
					}

					//UE and subtr done for each individual case below
					double fake_uncert_addition = 0, fake_uncert_addition_err = 0;
					if (apply_fake_uncert)
					{
						fake_uncert_addition = h_fake_rate->GetBinContent(i_trk_bin, i_jet_bin) * n_jets_data/n_jets_pp_mc * 0.3;
						fake_uncert_addition_err = h_fake_rate->GetBinError(i_trk_bin, i_jet_bin) * n_jets_data/n_jets_pp_mc * 0.3;
					}

					if (dataset_type == "pp")
					{
						UE = 0;
						UE_err = 0;

						//scale if necessary (in MC, n_jets factor = 1)
						fake = h_UE_TM->GetBinContent(i_trk_bin, i_jet_bin) * n_jets_data/n_jets_mc;
						fake_err = h_UE_TM->GetBinError(i_trk_bin, i_jet_bin) * n_jets_data/n_jets_mc;

						if (apply_fake_uncert)
						{
							fake = 1.3 * fake;
							fake_err = 1.3 * fake_err;
						}

						subtr = raw - fake;
						if (isMC) subtr_err = sqrt(fabs( pow(raw_err,2) - pow(fake_err,2) )); //errors fully correlated
						else subtr_err = sqrt(pow(raw_err,2) + pow(fake_err,2)); //errors uncorrelated
					}

					else if (dataset_type == "PbPb")
					{
						if (i_trk_bin >= trkpT_binning->FindBin(0.9) && i_trk_bin < trkpT_binning->FindBin(10.))
						{

							if (apply_UE_ConeMethod)
							{
								//cone method is scaled correctly in both data and mc
								//UE is uncorrelated to raw in both data and mc so add in quadrature

								UE = h_UE_corrected_cone->GetBinContent(i_trk_bin, i_jet_bin);
								UE_err = h_UE_corrected_cone->GetBinError(i_trk_bin, i_jet_bin); //corrected for UE-JER above

								subtr = raw - UE;
								subtr_err = sqrt(pow(raw_err,2) + pow(UE_err,2)); //errors uncorrelated for both data and MC
							}
							else if (apply_UE_mapStatSys)
							{
								if (!isMC)
								{
									//for mapStats systematic in data, use MB method, no jet scaling required, scale by rel. err from map_stat
									//UE is uncorrelated so errors are scaled by map stat factor and added in quadrature
									UE = h_UE_MB->GetBinContent(i_trk_bin, i_jet_bin) * h_map_stat->GetBinContent(i_trk_bin, i_jet_bin);
									UE_err = h_UE_MB->GetBinError(i_trk_bin, i_jet_bin) * h_map_stat->GetBinContent(i_trk_bin, i_jet_bin);

									subtr = raw - UE;
									subtr_err = sqrt(pow(raw_err,2) + pow(UE_err,2)); //errors uncorrelated for both data and MC
								}
								else if (isMC)
								{
									//for mapStats systematic in MC, use TM method, no jet scaling required, no map_stat required
									//UE is fully correlated so errors are subtracted in quadrature
									UE = h_UE_TM->GetBinContent(i_trk_bin, i_jet_bin);
									UE_err = h_UE_TM->GetBinError(i_trk_bin, i_jet_bin);

									subtr = raw - UE;
									subtr_err = sqrt(fabs(pow(raw_err,2) - pow(UE_err,2))); //errors uncorrelated for both data and MC
								}
							}
							else if (apply_fake_uncert)
							{
								if (!isMC)
								{
									//for fake systematic in data, use MB method and add fakes from pp scaled by n_jet * 0.3
									//UE is uncorrelated so UE+fake err is just errUE and errFake added in quadrature. errFake also needs to be scaled by n_jets and 0.3
									UE = h_UE_MB->GetBinContent(i_trk_bin, i_jet_bin) + fake_uncert_addition;
									UE_err = sqrt( pow(h_UE_MB->GetBinError(i_trk_bin, i_jet_bin),2) + pow(fake_uncert_addition_err,2) );

									subtr = raw - UE;
									subtr_err = sqrt(pow(raw_err,2) + pow(UE_err,2)); //errors uncorrelated for both data and MC
								}
								else if (isMC)
								{
									//for fake systematic in MC, use TM method and add fakes from pp scaled by n_jet * 0.3
									//Do error in two steps: subtr = raw - UE_err is corr, subtr - fake is uncorr

									UE = h_UE_TM->GetBinContent(i_trk_bin, i_jet_bin) + fake_uncert_addition;
									UE_err = h_UE_TM->GetBinError(i_trk_bin, i_jet_bin); //doing fake err separately

									subtr = raw - UE;
									subtr_err = sqrt(fabs(pow(raw_err,2) - pow(UE_err,2)) + pow(fake_uncert_addition_err,2)); //raw_err and UE_err correlated, fake_uncert_err uncorrelated
								}

							}
							else //no systematic applied, this is nominal behavior
							{
								if (!isMC)
								{
									//for nominal data, use MB method, no jet scaling required,
									//UE is uncorrelated so errors added in quadrature
									UE = h_UE_MB->GetBinContent(i_trk_bin, i_jet_bin);
									UE_err = h_UE_MB->GetBinError(i_trk_bin, i_jet_bin);

									subtr = raw - UE;
									subtr_err = sqrt(pow(raw_err,2) + pow(UE_err,2)); //errors uncorrelated for both data and MC
								}
								else if (isMC)
								{
									//in MC, use TM method, no jet scaling required
									//UE is fully correlated so errors are subtracted in quadrature
									UE = h_UE_TM->GetBinContent(i_trk_bin, i_jet_bin);
									UE_err = h_UE_TM->GetBinError(i_trk_bin, i_jet_bin);

									subtr = raw - UE;
									subtr_err = sqrt(fabs(pow(raw_err,2) - pow(UE_err,2))); //errors uncorrelated for both data and MC
								}
							}
						}
						else //above 10 GeV
						{
							raw = h_raw->GetBinContent(i_trk_bin, i_jet_bin);
							raw_err = h_raw->GetBinError(i_trk_bin, i_jet_bin);

							if (apply_fake_uncert)
							{
								if (!isMC)
								{
									//for fake systematic in data above 10 GeV, use TM method and add fakes from pp scaled by n_jet * 0.3
									//UE is uncorrelated so UE+fake err is just errUE and errFake added in quadrature. errFake also needs to be scaled by n_jets and 0.3

//									UE = (h_UE_TM->GetBinContent(i_trk_bin, i_jet_bin)*n_jets_data/n_jets_mc) + fake_uncert_addition;
//									UE_err = sqrt( pow(h_UE_TM->GetBinError(i_trk_bin, i_jet_bin) * n_jets_data/n_jets_mc,2) + pow(fake_uncert_addition_err,2) );

									UE = h_UE_MB->GetBinContent(i_trk_bin, i_jet_bin);
									UE_err = h_UE_MB->GetBinError(i_trk_bin, i_jet_bin);


									subtr = raw - UE;
									subtr_err = sqrt(pow(raw_err,2) + pow(UE_err,2)); //errors uncorrelated for both data and MC
								}
								else if (isMC)
								{
									//for fake systematic in MC, use TM method and add fakes from pp scaled by n_jet * 0.3
									//Do error in two steps: subtr = raw - UE_err is corr, subtr - fake is uncorr

									UE = h_UE_TM->GetBinContent(i_trk_bin, i_jet_bin) + fake_uncert_addition;
									UE_err = h_UE_TM->GetBinError(i_trk_bin, i_jet_bin); //doing fake err separately

									subtr = raw - UE;
									subtr_err = sqrt(fabs(pow(raw_err,2) - pow(UE_err,2)) + pow(fake_uncert_addition_err,2)); //raw_err and UE_err correlated, fake_uncert_err uncorrelated
								}

							}
							else //no systematic applied, this is nominal behavior
							{
								if (!isMC)
								{
									//for nominal data above 10 GeV, use TM method, jet scaling required
									//UE is uncorrelated so errors are scaled and added in quadrature
									UE = h_UE_TM->GetBinContent(i_trk_bin, i_jet_bin) * n_jets_data/n_jets_mc;
									UE_err = h_UE_TM->GetBinError(i_trk_bin, i_jet_bin) * n_jets_data/n_jets_mc;

									subtr = raw - UE;
									subtr_err = sqrt(pow(raw_err,2) + pow(UE_err,2)); //errors uncorrelated for both data and MC
								}
								else if (isMC)
								{
									//in MC, use TM method, no jet scaling required
									//UE is fully correlated so errors are subtracted in quadrature
									UE = h_UE_TM->GetBinContent(i_trk_bin, i_jet_bin);
									UE_err = h_UE_TM->GetBinError(i_trk_bin, i_jet_bin);

									subtr = raw - UE;
									subtr_err = sqrt(fabs(pow(raw_err,2) - pow(UE_err,2))); //errors uncorrelated for both data and MC
								}
							}
						}
					}


					if (verbose && subtr < 0 && (i_jet_bin > 7 && i_jet_bin < 12) && (i_trk_bin > 2 && i_trk_bin < 10))
					{
						std::string trk_label = Form("%1.1f < trk < %1.1f", trkpT_binning->GetBinLowEdge(i_trk_bin), trkpT_binning->GetBinUpEdge(i_trk_bin));
						std::string jet_label = Form("%1.0f < jet < %1.0f", jetpT_binning->GetBinLowEdge(i_jet_bin), jetpT_binning->GetBinUpEdge(i_jet_bin));
						std::string dr_label = Form("%1.2f < r < %1.2f", dR_binning->GetBinLowEdge(i_dR+1), dR_binning->GetBinUpEdge(i_dR+1));

						cout << Form("WARNING NEGATIVE CHPS at trk%i_cent%i_jet%i dR%i, %s, %s, %s, %s", i_trk_bin-1, i_cent, i_jet_bin-1, i_dR, trk_label.c_str(), jet_label.c_str(), dr_label.c_str(), num_to_cent(31,i_cent).c_str()) << endl;

					}

					h_final_UE->SetBinContent(i_trk_bin, i_jet_bin, UE);
					h_final_UE->SetBinError(i_trk_bin, i_jet_bin, UE_err);

					h_final_fake->SetBinContent(i_trk_bin, i_jet_bin, fake);
					h_final_fake->SetBinError(i_trk_bin, i_jet_bin, fake_err);

					h_raw_subtr->SetBinContent(i_trk_bin, i_jet_bin, subtr);
					h_raw_subtr->SetBinError(i_trk_bin, i_jet_bin, subtr_err);
				}
			}


			delete h_reco_jet_spect_mc;
			delete h_reco_jet_spect_data;


			//Unfolding
			TH2* h_raw_subtr_unf = (TH2*)h_raw_subtr->Clone(Form("h_raw_subtr_unf_dR%i_c%i", i_dR, i_cent));

			RooUnfoldResponse* r_response = (RooUnfoldResponse*)f_mc->Get(Form("ChPS_raw_0_dR%i_cent%i_ChPS_truth_dR%i_cent%i", i_dR, i_cent, i_dR, i_cent));

			TH2* h_response_matrix = r_response->Hresponse();
			name = Form("h_ChPS_response_matrix_jet_dR%i_cent%i", i_dR, i_cent);
			h_response_matrix->SetTitle(name.c_str());
			h_response_matrix->Write(name.c_str());

			if (h_raw_subtr_unf->GetEntries() != 0)
			{
				RooUnfoldBayes unfold_raw(r_response, h_raw_subtr_unf, n_unfold);
				unfold_raw.SetVerbose(0);
				h_raw_subtr_unf = (TH2D*)unfold_raw.Hreco(); //errors handled internally
				name = Form("h_raw_subtr_unf_dR%i_c%i", i_dR, i_cent);
				h_raw_subtr_unf->SetName(name.c_str());
				h_raw_subtr_unf->Sumw2();
			}


			delete r_response;

			//Bin by bin correction factors
			TH2* h_raw_subtr_unf_bbb = (TH2*)h_raw_subtr_unf->Clone(Form("h_raw_subtr_unf_bbb_dR%i_cent%i", i_dR, i_cent));

			for (int i_jet_bin = 0; i_jet_bin < N_jetpt; i_jet_bin++)
			{
				name = Form("h_bin_by_bin_cent%i_jetpt%i_dR%i", i_cent, i_jet_bin, i_dR);
				TH1* h_bin_by_bin = (TH1*)dr_factors->Get(name.c_str());
				h_bin_by_bin->SetName(name.c_str());

				for (int i_trk_bin = 0; i_trk_bin < N_trkpt; i_trk_bin++)
				{
					double bin_by_bin = h_bin_by_bin->GetBinContent(i_trk_bin+1);
					double bin_by_bin_err = h_bin_by_bin->GetBinError(i_trk_bin+1);

					if (bin_by_bin == 0) continue;
					double original, original_err, corrected, corrected_err;

					original = h_raw_subtr_unf_bbb->GetBinContent(i_trk_bin+1, i_jet_bin+1);
					original_err = h_raw_subtr_unf_bbb->GetBinError(i_trk_bin+1, i_jet_bin+1);
					if (original != 0)
					{
						corrected = original * bin_by_bin;
						corrected_err = fabs(corrected) * sqrt(pow(original_err/original,2) + pow(bin_by_bin_err/bin_by_bin,2));
						h_raw_subtr_unf_bbb->SetBinContent(i_trk_bin+1, i_jet_bin+1, corrected);
						h_raw_subtr_unf_bbb->SetBinError(i_trk_bin+1, i_jet_bin+1, corrected_err); //errors propagated
					}


				} //end trk bin loop
				delete h_bin_by_bin;
			} //end jet bin loop

			//doing per jet normalization
			for (int i_jet_bin = 0; i_jet_bin < N_jetpt; i_jet_bin++)
			{
				double n_jets_raw = h_reco_jet->GetBinContent(i_jet_bin+1);
				double n_jets_unf = h_reco_unfolded->GetBinContent(i_jet_bin+1);
				double n_jets_tru = h_truth_jet->GetBinContent(i_jet_bin+1);

				for (int i_trk_bin = 0; i_trk_bin < N_trkpt; i_trk_bin++)
				{
					double updated_raw = 0, updated_raw_err = 0, updated_raw_subtr = 0, updated_raw_subtr_err = 0, updated_UE = 0, updated_UE_err = 0, updated_fake = 0, updated_fake_err = 0;
					double updated_raw_subtr_unf = 0, updated_raw_subtr_unf_err = 0, updated_raw_subtr_unf_bbb = 0, updated_raw_subtr_unf_bbb_err = 0;
					double updated_truth = 0, updated_truth_err = 0;

					if (n_jets_raw != 0)
					{
						updated_raw = h_raw->GetBinContent(i_trk_bin+1,i_jet_bin+1) / n_jets_raw;
						updated_raw_err = h_raw->GetBinError(i_trk_bin+1,i_jet_bin+1) / n_jets_raw;
						updated_raw_subtr = h_raw_subtr->GetBinContent(i_trk_bin+1,i_jet_bin+1) / n_jets_raw;
						updated_raw_subtr_err = h_raw_subtr->GetBinError(i_trk_bin+1,i_jet_bin+1) / n_jets_raw;
						updated_UE = h_final_UE->GetBinContent(i_trk_bin+1,i_jet_bin+1) / n_jets_raw;
						updated_UE_err = h_final_UE->GetBinError(i_trk_bin+1,i_jet_bin+1) / n_jets_raw;
						updated_fake = h_final_fake->GetBinContent(i_trk_bin+1,i_jet_bin+1) / n_jets_raw;
						updated_fake_err = h_final_fake->GetBinError(i_trk_bin+1,i_jet_bin+1) / n_jets_raw;
					}

					//ChPS_final - norm by raw
					h_raw->SetBinContent(i_trk_bin+1,i_jet_bin+1, updated_raw);
					h_raw->SetBinError(i_trk_bin+1,i_jet_bin+1, updated_raw_err);

					//ChPS_subtr - norm by raw
					h_raw_subtr->SetBinContent(i_trk_bin+1,i_jet_bin+1, updated_raw_subtr);
					h_raw_subtr->SetBinError(i_trk_bin+1,i_jet_bin+1, updated_raw_subtr_err);

					//UE
					h_final_UE->SetBinContent(i_trk_bin+1,i_jet_bin+1, updated_UE);
					h_final_UE->SetBinError(i_trk_bin+1,i_jet_bin+1, updated_UE_err);

					//Fakes
					h_final_fake->SetBinContent(i_trk_bin+1,i_jet_bin+1, updated_fake);
					h_final_fake->SetBinError(i_trk_bin+1,i_jet_bin+1, updated_fake_err);


					if (n_jets_unf != 0)
					{
						updated_raw_subtr_unf = h_raw_subtr_unf->GetBinContent(i_trk_bin+1,i_jet_bin+1) / n_jets_unf;
						updated_raw_subtr_unf_err = h_raw_subtr_unf->GetBinError(i_trk_bin+1,i_jet_bin+1) / n_jets_unf;
						updated_raw_subtr_unf_bbb = h_raw_subtr_unf_bbb->GetBinContent(i_trk_bin+1,i_jet_bin+1) / n_jets_unf;
						updated_raw_subtr_unf_bbb_err = h_raw_subtr_unf_bbb->GetBinError(i_trk_bin+1,i_jet_bin+1) / n_jets_unf;
					}

					//ChPS_unf - norm by unfolded
					h_raw_subtr_unf->SetBinContent(i_trk_bin+1,i_jet_bin+1, updated_raw_subtr_unf);
					h_raw_subtr_unf->SetBinError(i_trk_bin+1,i_jet_bin+1, updated_raw_subtr_unf_err);

					//ChPS_final - norm by unfolded
					h_raw_subtr_unf_bbb->SetBinContent(i_trk_bin+1,i_jet_bin+1, updated_raw_subtr_unf_bbb);
					h_raw_subtr_unf_bbb->SetBinError(i_trk_bin+1,i_jet_bin+1, updated_raw_subtr_unf_bbb_err);


					if (n_jets_tru != 0)
					{
						updated_truth = h_truth->GetBinContent(i_trk_bin+1,i_jet_bin+1) / n_jets_tru;
						updated_truth_err = h_truth->GetBinError(i_trk_bin+1,i_jet_bin+1) / n_jets_tru;
					}

					//ChPS_truth
					h_truth->SetBinContent(i_trk_bin+1,i_jet_bin+1, updated_truth);
					h_truth->SetBinError(i_trk_bin+1,i_jet_bin+1, updated_truth_err);
				}
			}


			if (apply_MC_nonClos && !isMC)
			{
				TFile * f_MC_nonClosure = new TFile(Form("output_pdf_sys3/root/final_ChPS_MC_%s.root", dataset_type.c_str())); //mc nonclosure factors will come from non-reweighed response matrices


				for (int i_jet_bin = 0; i_jet_bin < N_jetpt; i_jet_bin++)
				{
					name = Form("h_ChPS_ratio_final_truth_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet_bin);
					TH1* h_MC_nonClosure = (TH1*)f_MC_nonClosure->Get(name.c_str());
					for (int i_trk_bin = 0; i_trk_bin < N_trkpt; i_trk_bin++)
					{
						double correction = h_MC_nonClosure->GetBinContent(i_trk_bin+1);
						double orig = h_raw_subtr_unf_bbb->GetBinContent(i_trk_bin+1,i_jet_bin+1);
						double orig_err = h_raw_subtr_unf_bbb->GetBinError(i_trk_bin+1,i_jet_bin+1);

						h_raw_subtr_unf_bbb->SetBinContent(i_trk_bin+1,i_jet_bin+1, orig/correction);
						h_raw_subtr_unf_bbb->SetBinError(i_trk_bin+1,i_jet_bin+1, orig_err/correction);

					}
				}
			}


			if (i_dR == 0)
			{
				h_raw_injet = (TH2*)h_raw->Clone(Form("raw_injet_c%i", i_cent));
				h_raw_subtr_injet = (TH2*)h_raw_subtr->Clone(Form("raw_subtr_injet_c%i", i_cent));
				h_raw_subtr_unf_injet = (TH2*)h_raw_subtr_unf->Clone(Form("raw_subtr_unf_injet_c%i", i_cent));
				h_raw_subtr_unf_bbb_injet = (TH2*)h_raw_subtr_unf_bbb->Clone(Form("raw_subtr_unf_bbb_injet_c%i", i_cent));

				h_UE_injet = (TH2*)h_final_UE->Clone(Form("UE_injet_c%i", i_cent));

				h_truth_injet = (TH2*)h_truth->Clone(Form("truth_injet_c%i", i_cent));
			}
			if (i_dR > 0 && i_dR < 7)
			{
				h_raw_injet->Add(h_raw);
				h_raw_subtr_injet->Add(h_raw_subtr);
				h_raw_subtr_unf_injet->Add(h_raw_subtr_unf);
				h_raw_subtr_unf_bbb_injet->Add(h_raw_subtr_unf_bbb);

				h_UE_injet->Add(h_final_UE);

				h_truth_injet->Add(h_truth);
			}


			double dR_lo = dR_binning->GetBinLowEdge(i_dR+1);
			double dR_hi = dR_binning->GetBinUpEdge(i_dR+1);
			double area = TMath::Pi() * ((dR_hi*dR_hi) - (dR_lo*dR_lo));

			f_output->cd();

			//get 1d ChPS
			for (int i_jet_bin = 0; i_jet_bin < N_jetpt; i_jet_bin++)
			{
				//raw_0
				name = Form("h_ChPS_raw_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet_bin);
				TH1* h_ChPS_raw = (TH1*)h_raw->ProjectionX(name.c_str(), i_jet_bin+1, i_jet_bin+1);
				h_ChPS_raw->SetTitle(name.c_str());
				h_ChPS_raw->Scale(1.,"width");
				h_ChPS_raw->Scale(1./area);
				h_ChPS_raw->Write(name.c_str());

				name = Form("h_ChPS_raw_subtr_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet_bin);
				TH1* h_ChPS_raw_subtr = (TH1*)h_raw_subtr->ProjectionX(name.c_str(), i_jet_bin+1, i_jet_bin+1);
				h_ChPS_raw_subtr->SetTitle(name.c_str());
				h_ChPS_raw_subtr->Scale(1.,"width");
				h_ChPS_raw_subtr->Scale(1./area);
				h_ChPS_raw_subtr->Write(name.c_str());

				name = Form("h_ChPS_raw_subtr_unf_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet_bin);
				TH1* h_ChPS_raw_subtr_unf = (TH1*)h_raw_subtr_unf->ProjectionX(name.c_str(), i_jet_bin+1, i_jet_bin+1);
				h_ChPS_raw_subtr_unf->SetTitle(name.c_str());
				h_ChPS_raw_subtr_unf->Scale(1.,"width");
				h_ChPS_raw_subtr_unf->Scale(1./area);
				h_ChPS_raw_subtr_unf->Write(name.c_str());

				name = Form("h_ChPS_raw_subtr_unf_bbb_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet_bin);
				TH1* h_ChPS_raw_subtr_unf_bbb = (TH1*)h_raw_subtr_unf_bbb->ProjectionX(name.c_str(), i_jet_bin+1, i_jet_bin+1);
				h_ChPS_raw_subtr_unf_bbb->SetTitle(name.c_str());
				h_ChPS_raw_subtr_unf_bbb->Scale(1.,"width");
				h_ChPS_raw_subtr_unf_bbb->Scale(1./area);
				h_ChPS_raw_subtr_unf_bbb->Write(name.c_str());

				//UE
				name = Form("h_ChPS_UE_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet_bin);
				TH1* h_ChPS_UE = (TH1*)h_final_UE->ProjectionX(name.c_str(), i_jet_bin+1, i_jet_bin+1);
				h_ChPS_UE->SetTitle(name.c_str());
				h_ChPS_UE->Scale(1.,"width");
				h_ChPS_UE->Scale(1./area);
				h_ChPS_UE->Write(name.c_str());

				//Fake
				name = Form("h_ChPS_fake_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet_bin);
				TH1* h_ChPS_fake = (TH1*)h_final_fake->ProjectionX(name.c_str(), i_jet_bin+1, i_jet_bin+1);
				h_ChPS_fake->SetTitle(name.c_str());
				h_ChPS_fake->Scale(1.,"width");
				h_ChPS_fake->Scale(1./area);
				h_ChPS_fake->Write(name.c_str());


				delete h_ChPS_raw;
				delete h_ChPS_raw_subtr;
				delete h_ChPS_raw_subtr_unf;
				delete h_ChPS_raw_subtr_unf_bbb;
				delete h_ChPS_UE;
				delete h_ChPS_fake;

				//truth
				name = Form("h_ChPS_truth_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet_bin);
				TH1* h_ChPS_truth = (TH1*)h_truth->ProjectionX(name.c_str(), i_jet_bin+1, i_jet_bin+1);
				h_ChPS_truth->SetTitle(name.c_str());
				h_ChPS_truth->Scale(1.,"width");
				h_ChPS_truth->Scale(1./area);
				h_ChPS_truth->Write(name.c_str());
				delete h_ChPS_truth;

			}

			delete h_raw;
			delete h_raw_subtr;
			delete h_raw_subtr_unf;
			delete h_raw_subtr_unf_bbb;

			delete h_truth;
			delete h_final_UE;
			delete h_final_fake;

		}

		//injet ChPS
		f_output->cd();
		double area_injet = TMath::Pi() * (0.4*0.4);
		for (int i_jet_bin = 0; i_jet_bin < N_jetpt; i_jet_bin++)
		{
			//raw_0
			name = Form("h_ChPS_raw_injet_cent%i_jetpt%i", i_cent, i_jet_bin);
			TH1* h_ChPS_raw_injet = (TH1*)h_raw_injet->ProjectionX(name.c_str(), i_jet_bin+1, i_jet_bin+1);
			h_ChPS_raw_injet->SetTitle(name.c_str());
			h_ChPS_raw_injet->Scale(1.,"width");
			h_ChPS_raw_injet->Scale(1./area_injet);
			h_ChPS_raw_injet->Write(name.c_str());
			delete h_ChPS_raw_injet;

			name = Form("h_ChPS_raw_subtr_injet_cent%i_jetpt%i", i_cent, i_jet_bin);
			TH1* h_ChPS_raw_subtr_injet = (TH1*)h_raw_subtr_injet->ProjectionX(name.c_str(), i_jet_bin+1, i_jet_bin+1);
			h_ChPS_raw_subtr_injet->SetTitle(name.c_str());
			h_ChPS_raw_subtr_injet->Scale(1.,"width");
			h_ChPS_raw_subtr_injet->Scale(1./area_injet);
			h_ChPS_raw_subtr_injet->Write(name.c_str());
			delete h_ChPS_raw_subtr_injet;

			name = Form("h_ChPS_raw_subtr_unf_injet_cent%i_jetpt%i", i_cent, i_jet_bin);
			TH1* h_ChPS_raw_subtr_unf_injet = (TH1*)h_raw_subtr_unf_injet->ProjectionX(name.c_str(), i_jet_bin+1, i_jet_bin+1);
			h_ChPS_raw_subtr_unf_injet->SetTitle(name.c_str());
			h_ChPS_raw_subtr_unf_injet->Scale(1.,"width");
			h_ChPS_raw_subtr_unf_injet->Scale(1./area_injet);
			h_ChPS_raw_subtr_unf_injet->Write(name.c_str());
			delete h_ChPS_raw_subtr_unf_injet;

			name = Form("h_ChPS_raw_subtr_unf_bbb_injet_cent%i_jetpt%i", i_cent, i_jet_bin);
			TH1* h_ChPS_raw_subtr_unf_bbb_injet = (TH1*)h_raw_subtr_unf_bbb_injet->ProjectionX(name.c_str(), i_jet_bin+1, i_jet_bin+1);
			h_ChPS_raw_subtr_unf_bbb_injet->SetTitle(name.c_str());
			h_ChPS_raw_subtr_unf_bbb_injet->Scale(1.,"width");
			h_ChPS_raw_subtr_unf_bbb_injet->Scale(1./area_injet);
			h_ChPS_raw_subtr_unf_bbb_injet->Write(name.c_str());
			delete h_ChPS_raw_subtr_unf_bbb_injet;

			//truth
			name = Form("h_ChPS_truth_injet_cent%i_jetpt%i", i_cent, i_jet_bin);
			TH1* h_ChPS_truth_injet = (TH1*)h_truth_injet->ProjectionX(name.c_str(), i_jet_bin+1, i_jet_bin+1);
			h_ChPS_truth_injet->SetTitle(name.c_str());
			h_ChPS_truth_injet->Scale(1.,"width");
			h_ChPS_truth_injet->Scale(1./area_injet);
			h_ChPS_truth_injet->Write(name.c_str());
			delete h_ChPS_truth_injet;

			//UE
			name = Form("h_ChPS_UE_injet_cent%i_jetpt%i", i_cent, i_jet_bin);
			TH1* h_ChPS_UE_injet = (TH1*)h_UE_injet->ProjectionX(name.c_str(), i_jet_bin+1, i_jet_bin+1);
			h_ChPS_UE_injet->SetTitle(name.c_str());
			h_ChPS_UE_injet->Scale(1.,"width");
			h_ChPS_UE_injet->Scale(1./area_injet);
			h_ChPS_UE_injet->Write(name.c_str());
			delete h_ChPS_UE_injet;
		}

		delete h_raw_injet;
		delete h_raw_subtr_injet;
		delete h_raw_subtr_unf_injet;
		delete h_raw_subtr_unf_bbb_injet;

		delete h_UE_injet;
		delete h_truth_injet;

		delete h_reco_jet;
		delete h_truth_jet;
		delete h_reco_unfolded;
	}

	cout << "######### DONE UNFOLDING #########" << endl << endl;

	return 0;
}
