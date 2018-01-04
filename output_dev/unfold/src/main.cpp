
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

int main()
{
	gErrorIgnoreLevel = 3001;
	cout << "Unfolding..." << endl;

	//	##############	Reading config	##############"
	TEnv *m_config = new TEnv();
	m_config->ReadFile("ff_config.cfg", EEnvLevel(1));
	m_config->Print();

	std::string dataset_type = "PbPb"; dataset_type = m_config->GetValue("dataset_type", dataset_type.c_str());
	std::string tracking_cut = "ppTight"; tracking_cut = m_config->GetValue("tracking_cut", tracking_cut.c_str());
	int centrality_scheme = 31; centrality_scheme = m_config->GetValue("centrality_scheme", centrality_scheme);
	int isMC = 1; isMC = m_config->GetValue("isMC", isMC);

	std::string did = "data";
	if (isMC) did = "MC_JZ_comb";

	//	##############	Config done	##############"

	TFile *f_mc = new TFile(Form("../raw_results/FF_MC_JZ_comb_out_histo_%s_5p02_r001.root", dataset_type.c_str()));
	TFile *f_data = new TFile(Form("../raw_results/FF_%s_out_histo_%s_5p02_r001.root", did.c_str(), dataset_type.c_str()));
	TFile *dr_factors = new TFile(Form("posCorr_factors_%s.root", dataset_type.c_str()));
	TFile *f_output = new TFile(Form("unfolded_%s_%s.root",did.c_str(), dataset_type.c_str()),"recreate");
	std::string name;

//	int N_CENT = 6;
	//remove when rerun is complete
	int N_Y = 4;
	n_cent_cuts = 6;

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

	TCanvas *c = new TCanvas("c","c",800,600);
	c->Divide(4,2);
	TLine *line = new TLine();
	line->SetLineStyle(3);
	line->SetLineColor(kBlack);

//	c->Print("inclusive_jet_spect_closure.pdf(","Title: Start");

	for (int i_cent = 0; i_cent < n_cent_cuts; i_cent++)
	{
		//ChPS: raw_0, rr(truth matched), truth, raw_unfolded, raw_rr_unfo
		//jet spect: reco, reco_matched, truth_matched, truth

		TH1 *h_reco_jet, *h_reco_jet_matched, *h_truth_jet, *h_truth_jet_matched, *h_reco_unfolded, *h_reco_matched_unfolded;


		for (int i_y = 0; i_y < N_Y; i_y++)
		{

			name = Form("h_reco_jet_spectrum_y%i_cent%i",i_y, i_cent);
			TH1* h_reco_jet_y_c = (TH1*)((TH1*)f_data->Get(name.c_str()))->Clone(Form("reco_jet_y%i_c%i",i_y, i_cent));
			h_reco_jet_y_c->Sumw2();

			name = Form("h_reco_jet_spectrum_matched_y%i_cent%i",i_y, i_cent);
			TH1* h_reco_jet_matched_y_c = (TH1*)((TH1*)f_mc->Get(name.c_str()))->Clone(Form("reco_jet_matched_y%i_c%i",i_y, i_cent));
			h_reco_jet_matched_y_c->Sumw2();

			name = Form("h_true_jet_spectrum_y%i_cent%i",i_y, i_cent);
			TH1* h_true_jet_y_c = (TH1*)((TH1*)f_mc->Get(name.c_str()))->Clone(Form("true_jet_y%i_c%i",i_y, i_cent));
			h_true_jet_y_c->Sumw2();

			name = Form("h_true_jet_spectrum_matched_y%i_cent%i",i_y, i_cent);
			TH1* h_true_jet_matched_y_c = (TH1*)((TH1*)f_mc->Get(name.c_str()))->Clone(Form("true_jet_matched_y%i_c%i",i_y, i_cent));
			h_true_jet_matched_y_c->Sumw2();

			TH1* h_unfolded_jet_y_c = (TH1*)h_true_jet_y_c->Clone(Form("h_raw_unfolded_y%i_c%i", i_y, i_cent));
            h_unfolded_jet_y_c->Reset();

			TH1* h_matched_unfolded_jet_y_c = (TH1*)h_true_jet_matched_y_c->Clone(Form("h_raw_matched_unfolded_y%i_c%i", i_y, i_cent));
			h_matched_unfolded_jet_y_c->Reset();

			//unfold only if using proper reco ChPS (not if using *reco matched* truth)
			name = Form("h_reco_jet_spectrum_y%i_cent%i_h_true_jet_spectrum_y%i_cent%i",i_y, i_cent, i_y, i_cent);
			RooUnfoldResponse* r_response = (RooUnfoldResponse*)f_mc->Get(name.c_str());

			RooUnfoldBayes unfold(r_response, h_reco_jet_y_c,4);
			unfold.SetVerbose(0);
			if (h_reco_jet_y_c->GetEntries() != 0)
			{
				h_unfolded_jet_y_c = (TH1D*)unfold.Hreco();
				name = Form("unfolded_y%i_c%i", i_y, i_cent);
				h_unfolded_jet_y_c->SetName(name.c_str());
				h_unfolded_jet_y_c->Sumw2();
			}

			RooUnfoldBayes matched_unfold(r_response, h_reco_jet_matched_y_c,4);
			matched_unfold.SetVerbose(0);
			if (h_reco_jet_matched_y_c->GetEntries() != 0)
			{
				h_matched_unfolded_jet_y_c = (TH1D*)matched_unfold.Hreco();
				name = Form("matched_unfolded_y%i_c%i", i_y, i_cent);
				h_matched_unfolded_jet_y_c->SetName(name.c_str());
				h_matched_unfolded_jet_y_c->Sumw2();
			}


			if (i_y == 4)
			{
				h_reco_jet = (TH1*)h_reco_jet_y_c->Clone(Form("h_inc_reco_jet_y%i_c%i",i_y, i_cent));
				h_reco_jet->Sumw2();

				h_reco_jet_matched = (TH1*)h_reco_jet_matched_y_c->Clone(Form("h_inc_reco_matched_jet_y%i_c%i",i_y, i_cent));
				h_reco_jet_matched->Sumw2();

				h_truth_jet = (TH1*)h_true_jet_y_c->Clone(Form("h_inc_truth_matched_jet_y%i_c%i",i_y, i_cent));
				h_truth_jet->Sumw2();

				h_truth_jet_matched = (TH1*)h_true_jet_matched_y_c->Clone(Form("h_inc_truth_jet_y%i_c%i",i_y, i_cent));
				h_truth_jet_matched->Sumw2();

				h_reco_unfolded = (TH1*)h_unfolded_jet_y_c->Clone(Form("h_inc_reco_unf_jet_y%i_c%i",i_y, i_cent));
				h_reco_unfolded->Sumw2();

				h_reco_matched_unfolded = (TH1*)h_matched_unfolded_jet_y_c->Clone(Form("h_inc_reco_matched_unf_jet_y%i_c%i",i_y, i_cent));
				h_reco_matched_unfolded->Sumw2();
			}

		}
	}

//			if (i_y == 4)
//			{
//				TH1* h_ratio_unf_raw = (TH1*)h_unfolded->Clone(Form("ratio_c%i_y%i", i_cent, i_y));
//				h_ratio->Divide(h_truth);
//
//
//				TH1* h_ratio = (TH1*)h_unfolded->Clone(Form("ratio_c%i_y%i", i_cent, i_y));
//				h_ratio->Divide(h_truth);
//
//				TH1* h_raw_ratio = (TH1*)h_raw->Clone(Form("raw_c%i_y%i", i_cent, i_y));
//				h_raw_ratio->Divide(h_truth);
//
//				c->cd(i_cent+1);
//				h_ratio->GetYaxis()->SetRangeUser(0.,2.);
//				h_ratio->GetXaxis()->SetRangeUser(80,1000);
//
//				h_ratio->SetMarkerStyle(20);
//				h_raw_ratio->SetMarkerStyle(20);
//
//				h_ratio->SetMarkerColor(kBlack);
//				h_raw_ratio->SetMarkerColor(kRed);
//
//				h_ratio->Draw();
//				h_raw_ratio->Draw("same");
//
//				line->DrawLine(80,1,1000,1);
//
//				gPad->SetLogx();
//				h_unfolded->Print("all");
//				c->Print("inclusive_jet_spect_closure.pdf");//,Form("Title: ratio_c%i_y%i", i_cent, i_y));
//			}
//
//                //Use n_jets from spectra inclusive in y
//                if (i_y == 4)
//                {
//                    h_raw_jet_spect_y4 = (TH1*)h_raw->Clone(Form("h_n_raw_jet_spect_y%i_cent%i",i_y, i_cent));
//                    h_unfolded_jet_spect_y4 = (TH1*)h_unfolded->Clone(Form("h_n_unf_jet_spect_y%i_cent%i",i_y, i_cent));
//                    h_truth_jet_spect_y4 = (TH1*)h_truth->Clone(Form("h_n_tru_jet_spect_y%i_cent%i",i_y, i_cent));;
//                }
//
//            }
//            f_output->cd();
//            h_raw->Write(Form("h_raw_jet_spect_cent%i_y%i",i_cent, i_y));
//            h_truth->Write(Form("h_true_jet_spect_cent%i_y%i",i_cent, i_y));
//            h_unfolded->Write(Form("h_unfolded_jet_spect_cent%i_y%i",i_cent, i_y));
//
//            cout << Form("h_unfolded_jet_spect_cent%i_y%i",i_cent, i_y) << endl;
//        }
//
//		TH2* h_raw_injet;
//		TH2* h_UE_injet;
//
//
//		for (int i_dR = 0; i_dR < N_dR; i_dR++)
//		{
//            cout << __LINE__ << " " << Form("%s_dR%i_cent%i",ChPS_raw_type.c_str(), i_dR, i_cent) << endl;
//
//			TH2* h_raw = (TH2*)f_data->Get(Form("%s_dR%i_cent%i",ChPS_raw_type.c_str(), i_dR, i_cent)); //ChPS_raw_0 or ChPS_raw_rr or ChPS_raw_rt or ChPS_raw_tr or ChPS_raw_tt
//            cout << __LINE__ << " " << Form("%s_dR%i_cent%i",ChPS_raw_type.c_str(), i_dR, i_cent) << endl;
//            TH2* h_truth = (TH2*)f_mc->Get(Form("ChPS_truth_dR%i_cent%i",i_dR, i_cent));
//			TH2* h_unfolded = (TH2*) h_raw->Clone(Form("nonUnfolded_ChPS_c%i_dR%i",i_cent, i_dR));
//            h_unfolded->Reset();
//			TH2* h_UE = (TH2*)f_mc->Get(Form("ff_UE_pT_dR%i_cent%i",i_dR, i_cent));
//            cout << __LINE__ << " " << Form("%s_dR%i_cent%i",ChPS_raw_type.c_str(), i_dR, i_cent) << endl;
//
//
//            h_raw->Sumw2();
//            h_truth->Sumw2();
//            cout << __LINE__ << " " << Form("%s_dR%i_cent%i",ChPS_raw_type.c_str(), i_dR, i_cent) << endl;
//
////            if ((dataset_type == "pp" && i_cent == 5) ||
////                dataset_type == "PbPb")
//            {
//
//                cout << Form("%s_dR%i_cent%i",ChPS_raw_type.c_str(), i_dR, i_cent) << endl;
//                //UE needs to take into account different number of reco jets in MC and data
//                TH1* h_reco_mc_jet_spect = (TH1*)((TH1*)f_mc->Get(Form("h_reco_jet_spectrum_y%i_cent%i", 4, i_cent)))->Clone(Form("norm_reco_mc_jet_y%i_cent%i", 4, i_cent));
//
//                h_reco_mc_jet_spect->SetName(Form("norm_reco_mc_jet_y4_cent%i", i_cent));
//                TH1* h_reco_data_jet_spect = (TH1*)((TH1*)f_data->Get(Form("h_reco_jet_spectrum_y%i_cent%i", 4, i_cent)))->Clone(Form("norm_reco_data_jet_y%i_cent%i", 4, i_cent));
//                h_reco_data_jet_spect->SetName(Form("norm_reco_data_jet_y4_cent%i", i_cent));
//
//                h_reco_mc_jet_spect->Sumw2();
//                h_reco_data_jet_spect->Sumw2();
//
//                for (int i_jet_bin = 1; i_jet_bin <= N_jetpt; i_jet_bin++)
//                {
//                    double n_jets_mc = h_reco_mc_jet_spect->GetBinContent(i_jet_bin);
//                    double n_jets_data = h_reco_data_jet_spect->GetBinContent(i_jet_bin);
//
//                    if (n_jets_mc == 0) continue;
//
//                    for (int i_trk_bin = 1; i_trk_bin <= N_jetpt; i_trk_bin++)
//                    {
//                        double updated_UE = h_UE->GetBinContent(i_trk_bin, i_jet_bin) * n_jets_data / n_jets_mc;
//                        double updated_UE_err = h_UE->GetBinError(i_trk_bin, i_jet_bin) * n_jets_data / n_jets_mc;
//
//                        h_UE->SetBinContent(i_trk_bin, i_jet_bin, updated_UE);
//                        h_UE->SetBinError(i_trk_bin, i_jet_bin, updated_UE_err); //scaled errors
//                    }
//                }
//                h_UE->Sumw2();
//
//                //Do UE subtraction and unfolding and bin by bin correction if using proper reco ChPS
//                if (ChPS_raw_type == "ChPS_raw_0" && do_UE_subtr) h_raw->Add(h_UE,-1); //End UE Subtraction
//
//                //Unfolding
//                if ((ChPS_raw_type == "ChPS_raw_0" || ChPS_raw_type == "ChPS_raw_rr") && do_unfolding)
//                {
//                    RooUnfoldResponse* r_response = (RooUnfoldResponse*)f_mc->Get(Form("ChPS_raw_0_dR%i_cent%i_ChPS_truth_dR%i_cent%i", i_dR, i_cent, i_dR, i_cent));
//                    RooUnfoldBayes unfold(r_response, h_raw, n_unfold);
//                    unfold.SetVerbose(0);
//
//                    h_unfolded = (TH2D*)unfold.Hreco(); //errors handled internally
//                    name = Form("Unfolded_ChPS_c%i_dR%i",i_cent, i_dR);
//                    h_unfolded->SetName(name.c_str());
//
//                } //End unfolding
//                else
//                {
//                    name = Form("nonUnfolded_ChPS_c%i_dR%i",i_cent, i_dR);
//                    h_unfolded = (TH2*)h_raw->Clone(name.c_str());
//                    h_unfolded->SetName(name.c_str());
//                }
//
//                h_unfolded->Sumw2();
//
//
//
//                //Bin by bin correction factors
//                if ((ChPS_raw_type == "ChPS_raw_0" || ChPS_raw_type == "ChPS_raw_rr") && do_BbB)
//                {
//                    for (int i_jet_bin = 0; i_jet_bin < N_jetpt; i_jet_bin++)
//                    {
//                        name = Form("h_bin_by_bin_cent%i_jetpt%i_dR%i",i_cent, i_jet_bin, i_dR);
//                        TH1* h_bin_by_bin = (TH1*)dr_factors->Get(name.c_str());
//                        h_bin_by_bin->SetName(name.c_str());
//
//                        cout << name << endl;
//                        for (int i_trk_bin = 0; i_trk_bin < N_trkpt; i_trk_bin++)
//                        {
//                            double original = h_unfolded->GetBinContent(i_trk_bin+1, i_jet_bin+1);
//                            double bin_by_bin = h_bin_by_bin->GetBinContent(i_trk_bin+1);
//                            double original_err = h_unfolded->GetBinError(i_trk_bin+1, i_jet_bin+1);
//                            double bin_by_bin_err = h_bin_by_bin->GetBinError(i_trk_bin+1);
//
//                            if (original == 0 || bin_by_bin == 0) continue;
//
//                            double corrected = original * bin_by_bin;
//                            double corrected_err = fabs(corrected) * sqrt(pow(original_err/original,2) + pow(bin_by_bin_err/bin_by_bin,2));
//
//                            h_unfolded->SetBinContent(i_trk_bin+1, i_jet_bin+1, corrected);
//                            h_unfolded->SetBinError(i_trk_bin+1, i_jet_bin+1, corrected_err); //errors propagated
//
//                        } //end trk bin loop
//                    } //end jet bin loop
//                }//End bin by bin correction (done after unfolding)
//
//
//                //doing per jet normalization
//                for (int i_jet_bin = 0; i_jet_bin < N_jetpt; i_jet_bin++)
//                {
//                    double n_jets_raw = h_raw_jet_spect_y4->GetBinContent(i_jet_bin+1);
//                    double n_jets_unf = h_unfolded_jet_spect_y4->GetBinContent(i_jet_bin+1);
//                    double n_jets_tru = h_truth_jet_spect_y4->GetBinContent(i_jet_bin+1);
//
//                    if (n_jets_raw == 0 || n_jets_unf == 0 || n_jets_tru == 0) continue;
//
//                    for (int i_trk_bin = 0; i_trk_bin < N_trkpt; i_trk_bin++)
//                    {
//                        double updated_raw = h_raw->GetBinContent(i_trk_bin+1,i_jet_bin+1) / n_jets_raw;
//                        double updated_reco = h_unfolded->GetBinContent(i_trk_bin+1,i_jet_bin+1) / n_jets_unf;
//                        double updated_truth = h_truth->GetBinContent(i_trk_bin+1,i_jet_bin+1) / n_jets_tru;
//                        double updated_UE = h_UE->GetBinContent(i_trk_bin+1,i_jet_bin+1) / n_jets_raw;
//
//                        double updated_raw_err =  h_raw->GetBinError(i_trk_bin+1,i_jet_bin+1) / n_jets_raw;
//                        double updated_reco_err = h_unfolded->GetBinError(i_trk_bin+1,i_jet_bin+1) / n_jets_unf;
//                        double updated_truth_err = h_truth->GetBinError(i_trk_bin+1,i_jet_bin+1) / n_jets_tru;
//                        double updated_UE_err = h_UE->GetBinError(i_trk_bin+1,i_jet_bin+1) / n_jets_raw;
//
//                        h_raw->SetBinContent(i_trk_bin+1, i_jet_bin+1, updated_raw);
//                        h_unfolded->SetBinContent(i_trk_bin+1, i_jet_bin+1, updated_reco);
//                        h_truth->SetBinContent(i_trk_bin+1, i_jet_bin+1, updated_truth);
//                        h_UE->SetBinContent(i_trk_bin+1, i_jet_bin+1, updated_UE);
//
//                        h_raw->SetBinError(i_trk_bin+1, i_jet_bin+1, updated_raw_err); //scaled errors
//                        h_unfolded->SetBinError(i_trk_bin+1, i_jet_bin+1, updated_reco_err); //scaled errors
//                        h_truth->SetBinError(i_trk_bin+1, i_jet_bin+1, updated_truth_err); //scaled errors
//                        h_UE->SetBinError(i_trk_bin+1, i_jet_bin+1, updated_UE_err); //scaled errors
//                    }
//                }
//
//                if (i_dR == 0)
//                {
//                    h_UE_injet = (TH2*)h_UE->Clone(Form("UE_inJet_cent%i", i_cent));;
//                    h_raw_injet = (TH2*)h_raw->Clone(Form("UE_inJet_cent%i", i_cent));;
//                }
//
//                if (i_dR > 0 && i_dR < 7)
//                {
//                    h_UE_injet->Add(h_UE);
//                    h_raw_injet->Add(h_raw);
//                }
//            }
//            cout << __LINE__ << " " << Form("%s_dR%i_cent%i",ChPS_raw_type.c_str(), i_dR, i_cent) << endl;
//
//			f_output->cd();
//
//			h_raw->Write(Form("h_%s_cent%i_dR%i",ChPS_raw_type.c_str(),i_cent, i_dR));
//            cout << __LINE__ << " " << Form("%s_dR%i_cent%i",ChPS_raw_type.c_str(), i_dR, i_cent) << endl;
//			h_truth->Write(Form("h_true_ChPS_cent%i_dR%i",i_cent, i_dR));
//            cout << __LINE__ << " " << Form("%s_dR%i_cent%i",ChPS_raw_type.c_str(), i_dR, i_cent) << endl;
//			h_unfolded->Write(Form("h_unfolded_%s_cent%i_dR%i",ChPS_raw_type.c_str(), i_cent, i_dR));
//            cout << __LINE__ << " " << Form("%s_dR%i_cent%i",ChPS_raw_type.c_str(), i_dR, i_cent) << endl;
//			h_UE->Write(Form("h_UE_cent%i_dR%i", i_cent, i_dR));
//            cout << __LINE__ << " " << Form("%s_dR%i_cent%i",ChPS_raw_type.c_str(), i_dR, i_cent) << endl;
//
//		}
//
//		f_output->cd();
//		h_raw_injet->Write(Form("h_%s_injet_cent%i",ChPS_raw_type.c_str(),i_cent));
//		h_UE_injet->Write(Form("h_UE_injet_cent%i", i_cent));
//	}
//
////	((TH3*)f_mc->Get("h_dR_change_jetpt0_cent0"))->Write("h_dR_change_jetpt0_cent0");
////	f_output->Close();
//
//

	cout << "Done unfolding" << endl;
	return 0;
}
