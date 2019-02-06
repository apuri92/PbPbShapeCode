#include "output_dev/functions/global_variables.h"

void get_UE(int sys_mode = 34)
{
	gErrorIgnoreLevel = 3001;

	string sys_path;
	sys_path = Form("c%i", sys_mode);

	TFile *input_file;

	//c18 and 19 are the files that have the UE raw histos to be made into etaphi maps
	input_file = new TFile(Form("output_dev/raw_results/%s/FF_MC_out_histo_PbPb_5p02_r001.root", sys_path.c_str()));

	TFile *output_file = new TFile(Form("./UE_MC_maps_%s.root", sys_path.c_str()), "recreate");

	cout << "Using file " << input_file->GetName() << endl;

	TAxis* dR_binning = (TAxis*)((TH3*)input_file->Get("h_dR_change_jetpt0_cent0"))->GetXaxis();
	TAxis* jetpT_binning = (TAxis*)((TH3*)input_file->Get("ChPS_raw_0_dR0_cent0"))->GetYaxis();
	TAxis* trkpT_binning = (TAxis*)((TH3*)input_file->Get("ChPS_raw_0_dR0_cent0"))->GetXaxis();
	TAxis* psi_binning = (TAxis*)((TH3*)input_file->Get("h_jet_v_Psi_cent0_jetpt9"))->GetXaxis();
	TAxis* eta_binning = (TAxis*)((TH3*)input_file->Get("h_jet_v_Psi_cent0_jetpt9"))->GetYaxis();
	TAxis* phi_binning = (TAxis*)((TH3*)input_file->Get("h_jet_v_Psi_cent0_jetpt9"))->GetZaxis();

	int N_dR = dR_binning->GetNbins();
	int N_jetpt = jetpT_binning->GetNbins();
	int N_trkpt = trkpT_binning->GetNbins();
	int N_Psi = psi_binning->GetNbins();
	int N_eta = eta_binning->GetNbins();
	int N_phi = phi_binning->GetNbins();

	output_file->cd();
	dR_binning->Write("dR_binning");
	jetpT_binning->Write("jetpT_binning");
	trkpT_binning->Write("trkpT_binning");

	string name;

	vector<vector<TH3*>> h_jet_v_Psi = vector<vector<TH3*>> (N_jetpt, vector<TH3*> (n_cent_cuts));
	vector<vector<vector<vector<TH3*>>>> h_UE_dNdEtadPhidpT = vector<vector<vector<vector<TH3*>>>> (N_jetpt, vector<vector<vector<TH3*>>> (N_Psi, vector<vector<TH3*>> (n_cent_cuts, vector<TH3*> (N_dR))));
	vector<vector<vector<vector<TH3*>>>> h_UE_dNdEtadPhidpT_norm = vector<vector<vector<vector<TH3*>>>> (N_jetpt, vector<vector<vector<TH3*>>> (N_Psi, vector<vector<TH3*>> (n_cent_cuts, vector<TH3*> (N_dR))));

	vector<vector<vector<vector<vector<TH2*>>>>> h_UE_eta_phi_maps =vector<vector<vector<vector<vector<TH2*>>>>> (N_jetpt, vector<vector<vector<vector<TH2*>>>> (N_trkpt, vector<vector<vector<TH2*>>> (N_Psi, vector<vector<TH2*>> (n_cent_cuts, vector<TH2*> (N_dR)))));

	vector<vector<double>> n_jets = vector<vector<double>> (N_eta, vector<double> (N_phi));

	int jet_start = 6;
	int jet_end = 12;


	if (sys_mode == 34 )
	{
		jet_start = 6;
		jet_end = 9;
	}
	if (sys_mode == 35 )
	{
		jet_start = 9;
		jet_end = 12;
	}

	for (int i_cent = 0; i_cent < n_cent_cuts-1; i_cent++)
	{
		for (int i_jet = jet_start; i_jet < jet_end; i_jet++)
		{

			name = Form("h_jet_v_Psi_cent%i_jetpt%i",i_cent, i_jet);
			h_jet_v_Psi[i_jet][i_cent] = (TH3*)input_file->Get(name.c_str());

			for (int i_psi = 0; i_psi < 10; i_psi++)
			{

				for (int i_eta = 0; i_eta < N_eta; i_eta++)
				{
					for (int i_phi = 0; i_phi < N_phi; i_phi++)
					{
						n_jets[i_eta][i_phi] = h_jet_v_Psi[i_jet][i_cent]->GetBinContent(i_psi+1, i_eta+1, i_phi+1);
					}
				}


				for (int i_dR = 0; i_dR < N_dR; i_dR++)
				{
					name = Form("h_UE_jetpt%i_dPsi%i_cent%i_dR%i",i_jet, i_psi, i_cent, i_dR);
					h_UE_dNdEtadPhidpT[i_jet][i_psi][i_cent][i_dR] = (TH3*)input_file->Get(name.c_str());
					h_UE_dNdEtadPhidpT_norm[i_jet][i_psi][i_cent][i_dR] = (TH3*)h_UE_dNdEtadPhidpT[i_jet][i_psi][i_cent][i_dR]->Clone(Form("%s_norm",name.c_str()));
					h_UE_dNdEtadPhidpT_norm[i_jet][i_psi][i_cent][i_dR]->Reset();

					for (int i_eta = 0; i_eta < N_eta; i_eta++)
					{
						for (int i_phi = 0; i_phi < N_phi; i_phi++)
						{
							for (int i_trk = 0; i_trk < N_trkpt; i_trk++)
							{

								if (trkpT_binning->GetBinLowEdge(i_trk+1) >= 10.) continue;
								if (n_jets[i_eta][i_phi] == 0)
								{
									h_UE_dNdEtadPhidpT_norm[i_jet][i_psi][i_cent][i_dR]->SetBinContent(i_trk+1, i_eta+1, i_phi+1, 0);
								}
								else
								{
									double orig = h_UE_dNdEtadPhidpT[i_jet][i_psi][i_cent][i_dR]->GetBinContent(i_trk+1, i_eta+1, i_phi+1);
									double normalized = orig/n_jets[i_eta][i_phi];
									h_UE_dNdEtadPhidpT_norm[i_jet][i_psi][i_cent][i_dR]->SetBinContent(i_trk+1, i_eta+1, i_phi+1, normalized);
								}
							}

						}
					}

					for (int i_trk = 0; i_trk < N_trkpt; i_trk++)
					{
						if (trkpT_binning->GetBinLowEdge(i_trk+1) >= 10.) continue;
						
						h_UE_dNdEtadPhidpT_norm[i_jet][i_psi][i_cent][i_dR]->GetXaxis()->SetRange(i_trk+1,i_trk+1);
						h_UE_eta_phi_maps[i_jet][i_trk][i_psi][i_cent][i_dR] = (TH2*)h_UE_dNdEtadPhidpT_norm[i_jet][i_psi][i_cent][i_dR]->Project3D("zy");

						output_file->cd();
						string dr_label = Form("%1.2f < r < %1.2f", dR_binning->GetBinLowEdge(i_dR+1), dR_binning->GetBinUpEdge(i_dR+1));
						string psi_label = Form("%1.2f < psi < %1.2f", psi_binning->GetBinLowEdge(i_psi+1), psi_binning->GetBinUpEdge(i_psi+1));
						string trk_label = Form("%1.2f < p_{T}^{Trk} < %1.2f", trkpT_binning->GetBinLowEdge(i_trk+1), trkpT_binning->GetBinUpEdge(i_trk+1));
						string cent_label = num_to_cent(31,i_cent);

						h_UE_eta_phi_maps[i_jet][i_trk][i_psi][i_cent][i_dR]->SetTitle(Form("UE: %s, %s, %s, %s",dr_label.c_str(), psi_label.c_str(), trk_label.c_str(), cent_label.c_str()));

						name = Form("h_UE_new_MC_dR%i_dPsi%i_pt%i_cent%i_jet%i",i_dR, i_psi, i_trk+1, i_cent, i_jet);
						h_UE_eta_phi_maps[i_jet][i_trk][i_psi][i_cent][i_dR]->SetName(Form("%s",name.c_str()));
						h_UE_eta_phi_maps[i_jet][i_trk][i_psi][i_cent][i_dR]->Write(Form("%s",name.c_str()));

					}


				}

			}
		}
		cout << Form("Done cent%i", i_cent) << endl;
	}

}










