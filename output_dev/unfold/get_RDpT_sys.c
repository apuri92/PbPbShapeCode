#include "../functions/global_variables.h"
#include "draw_functions.c"
#include "TEnv.h"

void get_RDpT_sys()
{

	string name;
	TFile* pp_sys = new TFile(Form("output_pdf_nominal/root/final_ChPS_sys_data_ppx.root"));
	TFile* PbPb_sys = new TFile(Form("output_pdf_nominal/root/final_ChPS_sys_data_PbPbx.root"));
	TFile* RDpT_sys = new TFile(Form("output_pdf_nominal/root/final_RDpT_sys_data.root"));
	TFile* RDpT_sysx = new TFile(Form("output_pdf_nominal/root/final_RDpT_sys_datax.root"), "recreate");

	TAxis* dR_binning = (TAxis*)pp_sys->Get("dR_binning");
	TAxis* jetpT_binning = (TAxis*)pp_sys->Get("jetpT_binning");
	TAxis* trkpT_binning = (TAxis*)pp_sys->Get("trkpT_binning");

	int N_dR = dR_binning->GetNbins();
	int N_jetpt = jetpT_binning->GetNbins();
	int N_trkpt = trkpT_binning->GetNbins();

	RDpT_sysx->cd();
	dR_binning->Write("dR_binning");
	jetpT_binning->Write("jetpT_binning");
	trkpT_binning->Write("trkpT_binning");

	vector<string> combined_sys_names; vector<bool> sys_correlated;
	combined_sys_names.push_back("JER"); sys_correlated.push_back(1);
	combined_sys_names.push_back("JES"); sys_correlated.push_back(1);
	combined_sys_names.push_back("UE"); sys_correlated.push_back(1);
	combined_sys_names.push_back("Tracking"); sys_correlated.push_back(1);
	combined_sys_names.push_back("Unfolding"); sys_correlated.push_back(0);
	combined_sys_names.push_back("MCNonClosure"); sys_correlated.push_back(1);


	vector<vector<vector<vector<TH1*>>>> h_PbPb_comb_sys_p (combined_sys_names.size(), vector<vector<vector<TH1*>>> (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt))));
	vector<vector<vector<vector<TH1*>>>> h_PbPb_comb_sys_n (combined_sys_names.size(), vector<vector<vector<TH1*>>> (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt))));

	vector<vector<vector<vector<TH1*>>>> h_pp_comb_sys_p (combined_sys_names.size(), vector<vector<vector<TH1*>>> (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt))));
	vector<vector<vector<vector<TH1*>>>> h_pp_comb_sys_n (combined_sys_names.size(), vector<vector<vector<TH1*>>> (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt))));

	vector<vector<vector<vector<TH1*>>>> x_RDpT_comb_sys_p (combined_sys_names.size(), vector<vector<vector<TH1*>>> (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt))));
	vector<vector<vector<vector<TH1*>>>> x_RDpT_comb_sys_n (combined_sys_names.size(), vector<vector<vector<TH1*>>> (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt))));

	vector<vector<vector<vector<TH1*>>>> h_RDpT_comb_sys_p (combined_sys_names.size(), vector<vector<vector<TH1*>>> (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt))));
	vector<vector<vector<vector<TH1*>>>> h_RDpT_comb_sys_n (combined_sys_names.size(), vector<vector<vector<TH1*>>> (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt))));

	vector<vector<vector<TH1*>>> h_RDpT_total_p (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_RDpT_total_n (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));


	int jet_pt_start = 7;
	int jet_pt_end = 11;
	int trk_pt_start = 2;
	int trk_pt_end = 9;


	for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
	{
		cout << Form("Done jet%i", i_jet) << endl;
		for (int i_cent = 0; i_cent < 6; i_cent++)
		{
			for (int i_trk = 0; i_trk < N_trkpt; i_trk++)
			{
				for (int i_comb_sys = 0; i_comb_sys < combined_sys_names.size(); i_comb_sys++)
				{

					//Get pp Sys
					if (i_cent == 0)
					{
						name = Form("h_ChPS_sys_trk%i_cent%i_jetpt%i_%s_p", i_trk, 6, i_jet, combined_sys_names[i_comb_sys].c_str());
						h_pp_comb_sys_p[i_comb_sys][i_trk][6][i_jet] = (TH1*)(TH1*)pp_sys->Get(name.c_str())->Clone(Form("%s_pp", name.c_str()));

						name = Form("h_ChPS_sys_trk%i_cent%i_jetpt%i_%s_n", i_trk, 6, i_jet,combined_sys_names[i_comb_sys].c_str());
						h_pp_comb_sys_n[i_comb_sys][i_trk][6][i_jet] = (TH1*)(TH1*)pp_sys->Get(name.c_str())->Clone(Form("%s_pp", name.c_str()));
					}

					//Get PbPb Sys
					name = Form("h_ChPS_sys_trk%i_cent%i_jetpt%i_%s_p", i_trk, i_cent, i_jet,combined_sys_names[i_comb_sys].c_str());
					h_PbPb_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet] = (TH1*)(TH1*)PbPb_sys->Get(name.c_str())->Clone(Form("%s_PbPb", name.c_str()));

					name = Form("h_ChPS_sys_trk%i_cent%i_jetpt%i_%s_n", i_trk, i_cent, i_jet,combined_sys_names[i_comb_sys].c_str());
					h_PbPb_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet] = (TH1*)(TH1*)PbPb_sys->Get(name.c_str())->Clone(Form("%s_PbPb", name.c_str()));

//					//Get PbPb Sys
//					name = Form("h_RDpT_sys_trk%i_cent%i_jetpt%i_%s_p", i_trk, i_cent, i_jet,combined_sys_names[i_comb_sys].c_str());
//					x_RDpT_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet] = (TH1*)(TH1*)RDpT_sysx->Get(name.c_str())->Clone(Form("%s_RDpT", name.c_str()));
//
//					name = Form("h_RDpT_sys_trk%i_cent%i_jetpt%i_%s_n", i_trk, i_cent, i_jet,combined_sys_names[i_comb_sys].c_str());
//					x_RDpT_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet] = (TH1*)(TH1*)RDpT_sysx->Get(name.c_str())->Clone(Form("%s_RDpT", name.c_str()));
//


					//setup RDpT Sys
					name = Form("h_RDpT_sys_trk%i_cent%i_jetpt%i_%s_p", i_trk, i_cent, i_jet,combined_sys_names[i_comb_sys].c_str());
					h_RDpT_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet] = (TH1*)h_PbPb_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->Clone(name.c_str());
					h_RDpT_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->Reset();

					name = Form("h_RDpT_sys_trk%i_cent%i_jetpt%i_%s_n", i_trk, i_cent, i_jet,combined_sys_names[i_comb_sys].c_str());
					h_RDpT_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet] = (TH1*)h_PbPb_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->Clone(name.c_str());
					h_RDpT_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->Reset();

					//setup RDpT Total
					if (i_comb_sys == 0)
					{
						name = Form("h_RDpT_sys_trk%i_cent%i_jetpt%i_total_p", i_trk, i_cent, i_jet);
						h_RDpT_total_p[i_trk][i_cent][i_jet] = (TH1*)h_PbPb_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->Clone(name.c_str());
						h_RDpT_total_p[i_trk][i_cent][i_jet]->Reset();

						name = Form("h_RDpT_sys_trk%i_cent%i_jetpt%i_total_n", i_trk, i_cent, i_jet);
						h_RDpT_total_n[i_trk][i_cent][i_jet] = (TH1*)h_PbPb_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->Clone(name.c_str());
						h_RDpT_total_n[i_trk][i_cent][i_jet]->Reset();
					}

					if (i_trk == 3 && i_cent == 0 && i_jet == 7 && combined_sys_names[i_comb_sys] == "JER")
					{
//						h_PbPb_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->Print("all");
//						h_PbPb_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->Print("all");
//
//						h_pp_comb_sys_p[i_comb_sys][i_trk][6][i_jet]->Print("all");
//						h_pp_comb_sys_n[i_comb_sys][i_trk][6][i_jet]->Print("all");
//
//						x_RDpT_comb_sys_p[i_comb_sys][i_trk][6][i_jet]->Print("all");
//						x_RDpT_comb_sys_n[i_comb_sys][i_trk][6][i_jet]->Print("all");
					}



					double pbpb_err_p, pp_err_p, RDpT_err_p, pbpb_err_n, pp_err_n, RDpT_err_n;
					double pbpb_err, pp_err, RDpT_err;
					double x_RDpT_p, x_RDpT_n;
					for (int i_dR = 1; i_dR <= N_dR; i_dR++)
					{
						if (sys_correlated[i_comb_sys])
						{
							pbpb_err_p = h_PbPb_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->GetBinContent(i_dR);
							pp_err_p = h_pp_comb_sys_p[i_comb_sys][i_trk][6][i_jet]->GetBinContent(i_dR);

							pbpb_err_n = h_PbPb_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->GetBinContent(i_dR);
							pp_err_n = h_pp_comb_sys_n[i_comb_sys][i_trk][6][i_jet]->GetBinContent(i_dR);

							RDpT_err_p = (1.+pbpb_err_p)/(1.+pp_err_p) - 1;
							RDpT_err_n = (1.+pbpb_err_n)/(1.+pp_err_n) - 1;

							x_RDpT_p = x_RDpT_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->GetBinContent(i_dR);
							x_RDpT_n = x_RDpT_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->GetBinContent(i_dR);

							if (i_trk == 2 && i_jet == 7 && i_cent == 0)
							{
								cout << "******************" << endl;
								cout << x_RDpT_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->GetName() << endl;
								cout << Form("%f, %f -> %f [%f]",pbpb_err_p, pp_err_p, RDpT_err_p, x_RDpT_p ) << endl;
								cout << "------" << endl;
								cout << x_RDpT_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->GetName() << endl;
								cout << Form("%f, %f -> %f [%f]",pbpb_err_p, pp_err_p, RDpT_err_p, x_RDpT_p ) << endl;
							}


//							if (RDpT_err_p > 0 && RDpT_err_n > 0)
//							{
//								RDpT_err_p = sqrt(pow(RDpT_err_p,2) + pow(RDpT_err_n,2));
//								RDpT_err_n = 0;
//							}
//
//							if (RDpT_err_p < 0 && RDpT_err_n < 0)
//							{
//								RDpT_err_n = -sqrt(pow(RDpT_err_p,2) + pow(RDpT_err_n,2));
//								RDpT_err_p = 0;
//							}


							h_RDpT_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->SetBinContent(i_dR, RDpT_err_p);
							h_RDpT_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->SetBinContent(i_dR, RDpT_err_n);
						}
						else
						{
							pbpb_err = h_PbPb_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->GetBinContent(i_dR);
							pp_err = h_pp_comb_sys_p[i_comb_sys][i_trk][6][i_jet]->GetBinContent(i_dR);
							RDpT_err = sqrt(pow((pbpb_err),2) + pow((pp_err),2) );
							h_RDpT_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->SetBinContent(i_dR, RDpT_err);

							pbpb_err = h_PbPb_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->GetBinContent(i_dR);
							pp_err = h_pp_comb_sys_n[i_comb_sys][i_trk][6][i_jet]->GetBinContent(i_dR);
							RDpT_err = sqrt(pow((pbpb_err),2) + pow((pp_err),2) );
							h_RDpT_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->SetBinContent(i_dR, -RDpT_err);
						}


					}

					RDpT_sysx->cd();
					name = Form("h_RDpT_sys_trk%i_cent%i_jetpt%i_%s_p", i_trk, i_cent, i_jet,combined_sys_names[i_comb_sys].c_str());
					h_RDpT_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->Write(name.c_str());

					name = Form("h_RDpT_sys_trk%i_cent%i_jetpt%i_%s_n", i_trk, i_cent, i_jet,combined_sys_names[i_comb_sys].c_str());
					h_RDpT_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->Write(name.c_str());

//					if (i_trk == 3 && i_cent == 0 && i_jet == 7 && combined_sys_names[i_comb_sys] == "JER")
//					{
//						h_RDpT_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->Print("all");
//						h_RDpT_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->Print("all");
//					}

				}


				for (int i_dR = 1; i_dR <= N_dR; i_dR++)
				{
					double tmp;
					for (int i_comb_sys = 0; i_comb_sys < combined_sys_names.size(); i_comb_sys++)
					{

						if (i_comb_sys == 0)
						{
							tmp = pow(h_RDpT_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->GetBinContent(i_dR),2);
							h_RDpT_total_p[i_trk][i_cent][i_jet]->SetBinContent(i_dR, sqrt(tmp) );

							tmp = pow(h_RDpT_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->GetBinContent(i_dR),2);
							h_RDpT_total_n[i_trk][i_cent][i_jet]->SetBinContent(i_dR, -sqrt(tmp) );
						}
						else
						{
							tmp = pow(h_RDpT_total_p[i_trk][i_cent][i_jet]->GetBinContent(i_dR),2);
							tmp = tmp + pow(h_RDpT_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->GetBinContent(i_dR),2);
							h_RDpT_total_p[i_trk][i_cent][i_jet]->SetBinContent(i_dR, sqrt(tmp) );

							tmp = pow(h_RDpT_total_n[i_trk][i_cent][i_jet]->GetBinContent(i_dR),2);
							tmp = tmp + pow(h_RDpT_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->GetBinContent(i_dR),2);
							h_RDpT_total_n[i_trk][i_cent][i_jet]->SetBinContent(i_dR, -sqrt(tmp) );
						}

					}
				}

				name = Form("h_RDpT_sys_trk%i_cent%i_jetpt%i_total_p", i_trk, i_cent, i_jet);
				h_RDpT_total_p[i_trk][i_cent][i_jet]->Write(name.c_str());

				name = Form("h_RDpT_sys_trk%i_cent%i_jetpt%i_total_n", i_trk, i_cent, i_jet);
				h_RDpT_total_n[i_trk][i_cent][i_jet]->Write(name.c_str());


			}
		}
	}


}

