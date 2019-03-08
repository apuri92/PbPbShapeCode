#include "output_dev/functions/global_variables.h"

void mapSystematic(int sys_mode = 87)
{
	SetAtlasStyle();
	gErrorIgnoreLevel = 3001;
	string name;

	string sys_path;
	double r_max_range = 0.8;
	if (sys_mode == 0) sys_path = Form("nominal");
	if (sys_mode > 0 && sys_mode < 100) sys_path = Form("c%i", sys_mode);
	if (sys_mode > 100) sys_path = Form("sys%i", sys_mode);

	TFile *input_file_data = new TFile(Form("output_dev/raw_results/%s/FF_data_out_histo_PbPb_5p02_r001.root", sys_path.c_str()));
	name = Form("./UE_MapSystematic_%s.root", sys_path.c_str());
	TFile *output_file = new TFile(name.c_str(), "recreate");

	cout << "Using file " << input_file_data->GetName() << endl;

	TAxis* dR_binning = (TAxis*)((TH3*)input_file_data->Get("h_dR_change_jetpt0_cent0"))->GetXaxis();
	TAxis* jetpT_binning = (TAxis*)((TH3*)input_file_data->Get("ChPS_raw_0_dR0_cent0"))->GetYaxis();
	TAxis* trkpT_binning = (TAxis*)((TH3*)input_file_data->Get("ChPS_raw_0_dR0_cent0"))->GetXaxis();

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

	int N_sys = 100;
	for (int i_cent = 0; i_cent < 6; i_cent++)
	{
		string cent_label = num_to_cent(31,i_cent);

		for (int i_dR = 0; i_dR < N_dR; i_dR++)
		{
			string dr_label = Form("%1.2f < r < %1.2f", dR_binning->GetBinLowEdge(i_dR+1), dR_binning->GetBinUpEdge(i_dR+1));

			name = Form("ChPS_MB_UE_dR%i_cent%i", i_dR, i_cent);
			TH2* h_nom = (TH2*)input_file_data->Get(name.c_str());
			TH2* h_sys_multiplier = (TH2*)h_nom->Clone(Form("%s_statSys",name.c_str()));
			h_sys_multiplier->Reset();
			for (int i_trk = 2; i_trk <= 8; i_trk++)
			{
				string trk_label = Form("%1.1f < p_{T}^{Trk} < %1.1f GeV", trkpT_binning->GetBinLowEdge(i_trk), trkpT_binning->GetBinUpEdge(i_trk));

				for (int i_jet = 7; i_jet <= 12; i_jet++)
				{
					string jet_label = Form("%1.0f < p_{T}^{Jet} < %1.0f GeV", jetpT_binning->GetBinLowEdge(i_jet), jetpT_binning->GetBinUpEdge(i_jet));

					double nominal = h_nom->GetBinContent(i_trk, i_jet);
					if (nominal == 0) continue;

					TH1* h_gaus = new TH1D("h_gaus","h_gaus",500,-2,2);

					//get gaussian for each bin using systematic variations
					for (int i_sys = 0; i_sys < N_sys; i_sys++)
					{
						name = Form("ChPS_MB_UE_sys%i_dR%i_cent%i",i_sys, i_dR, i_cent);
						TH2* h_sys = (TH2*)input_file_data->Get(name.c_str());
						double sys = (h_sys->GetBinContent(i_trk, i_jet) - nominal)/nominal;
						h_gaus->Fill(sys);
						delete h_sys;
					}

					//fit gaussian and get mean, put mean into TH2 histogram that will be lookup table for systematic uncertainty
					double mean = h_gaus->GetMean();
					double rms = h_gaus->GetRMS();

					TF1* f_gaus = new TF1("fit","gaus");
					h_gaus->Fit(f_gaus,"Q","");

					h_sys_multiplier->SetBinContent(i_trk, i_jet,f_gaus->GetParameter(2)+1); //+1 to get relative error using width in terms of percentage
					output_file->cd();
					name = Form("ChPS_MB_UE_dR%i_cent%i_bins_trk%i_jet%i", i_dR, i_cent, i_trk, i_jet);
					h_gaus->SetTitle(Form("%s %s %s %s", cent_label.c_str(), dr_label.c_str(), trk_label.c_str(), jet_label.c_str() ));
					h_gaus->Write(name.c_str());
					delete f_gaus;
					delete h_gaus;

				}

				cout << "Done " << Form("ChPS_MB_UE_dR%i_cent%i", i_dR, i_cent) << endl;
			}

			output_file->cd();
			name = Form("ChPS_MB_UE_dR%i_cent%i_statSys", i_dR, i_cent);
			h_sys_multiplier->Write(name.c_str());


			delete h_sys_multiplier;
			delete h_nom;
		}
	}


}
