#include "../functions/global_variables.h"
#include "draw_functions.c"
#include "TEnv.h"

void systematics(string config_file = "sys_config.cfg")
{
	cout << "######### DOING Systematics #########" << endl;
	SetAtlasStyle();
	gErrorIgnoreLevel = 3001;
	string name = "";

	//	##############	Reading config	##############"
	TEnv *m_config = new TEnv();
	m_config->ReadFile(config_file.c_str(), EEnvLevel(1));

	std::string mode = "RDpT"; mode = m_config->GetValue("mode", mode.c_str());
	int isMC = 1; isMC = m_config->GetValue("isMC", isMC);

	std::string dataset_type = "PbPb"; dataset_type = m_config->GetValue("dataset_type", dataset_type.c_str());
	int centrality_scheme = 31; centrality_scheme = m_config->GetValue("centrality_scheme", centrality_scheme);
	int verbose = 0; verbose = m_config->GetValue("verbose", verbose);

	std::string did = "data";
	if (isMC) did = "MC";

	if (verbose) m_config->Print();
	//	##############	Config done	##############"


	if (mode == "RDpT") dataset_type = "";
	else dataset_type = Form("_%s", dataset_type.c_str());

	TFile* nom_file = new TFile(Form("output_pdf_nominal/root/final_%s_%s%s.root", mode.c_str(), did.c_str(), dataset_type.c_str()));
	TFile* output_file = new TFile(Form("output_pdf_nominal/root/final_%s_sys_%s%s.root", mode.c_str(), did.c_str(), dataset_type.c_str()), "recreate");

	cout << "Using files:" << endl;
	cout << nom_file->GetName() << endl;

	bool doSmall = true;
	if (dataset_type == "_pp") doSmall = false;
	TAxis* dR_binning = (TAxis*)nom_file->Get("dR_binning");
	TAxis* jetpT_binning = (TAxis*)nom_file->Get("jetpT_binning");
	TAxis* trkpT_binning = (TAxis*)nom_file->Get("trkpT_binning");

	int N_dR = dR_binning->GetNbins();
	int N_jetpt = jetpT_binning->GetNbins();
	int N_trkpt = trkpT_binning->GetNbins();

	output_file->cd();
	dR_binning->Write("dR_binning");
	jetpT_binning->Write("jetpT_binning");
	trkpT_binning->Write("trkpT_binning");

	vector<TFile*> sys_files;
	vector<string> sys_names;
	vector<string> combined_sys_names;

	sys_names.push_back("sys101"); //JER
	sys_names.push_back("sys102"); //Sign
	sys_names.push_back("sys103"); //Unfolding
	sys_names.push_back("sys105"); //MCProbCut
	sys_names.push_back("sys106"); //HIJES_1_P
	sys_names.push_back("sys107"); //HIJES_2_P
	sys_names.push_back("sys108"); //HIJES_1_N
	sys_names.push_back("sys109"); //HIJES_2_N
	sys_names.push_back("sys110"); //Material_P
	sys_names.push_back("sys111"); //Material_N
	sys_names.push_back("sys114"); //Tracking
	sys_names.push_back("sys115"); //TrackingRes
	sys_names.push_back("sys116"); //CentHIJES_P
	sys_names.push_back("sys117"); //CentHIJES_N
	sys_names.push_back("sys118"); //ppJES_P
	sys_names.push_back("sys119"); //ppJES_N
	sys_names.push_back("sys200"); //UE
//	sys_names.push_back("sys201"); //Fakes


	combined_sys_names.push_back("JER");
	combined_sys_names.push_back("JES");
	combined_sys_names.push_back("UE");
	combined_sys_names.push_back("Tracking");
	combined_sys_names.push_back("Unfolding");


	vector<vector<vector<TH1*>>> h_nom (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_total_sys_p (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_total_sys_n (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));


	//individual systematics
	vector<vector<vector<vector<TH1*>>>> h_sys (sys_names.size(), vector<vector<vector<TH1*>>> (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt))));
	vector<vector<vector<vector<TH1*>>>> h_sys_p (sys_names.size(), vector<vector<vector<TH1*>>> (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt))));
	vector<vector<vector<vector<TH1*>>>> h_sys_n (sys_names.size(), vector<vector<vector<TH1*>>> (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt))));

	//combined systematics
	vector<vector<vector<vector<TH1*>>>> h_comb_sys_p (combined_sys_names.size(), vector<vector<vector<TH1*>>> (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt))));
	vector<vector<vector<vector<TH1*>>>> h_comb_sys_n (combined_sys_names.size(), vector<vector<vector<TH1*>>> (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt))));


	int jet_pt_start = 7;
	int jet_pt_end = 11;
	int trk_pt_start = 2;
	int trk_pt_end = 9;

	cout << "Systematics from:" << endl;

	for (int i_sys = 0; i_sys < sys_names.size(); i_sys++)
	{
		name = Form("output_pdf_%s/root/final_%s_%s%s.root", sys_names[i_sys].c_str(), mode.c_str(), did.c_str(), dataset_type.c_str());
		sys_files.push_back( new TFile( name.c_str() ) );
		cout << sys_files[i_sys]->GetName() << endl;

		for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
		{
			for (int i_cent = 0; i_cent < n_cent_cuts; i_cent++)
			{
				if ((dataset_type == "_PbPb" || mode == "RDpT") && i_cent == 6) continue;
				if (dataset_type == "_pp" && i_cent < 6) continue;

				for (int i_trk = 0; i_trk < N_trkpt; i_trk++)
				{
					
					name = Form("h_%s_final_indR_trk%i_cent%i_jetpt%i",mode.c_str(), i_trk, i_cent, i_jet);
					if (i_sys == 0) h_nom[i_trk][i_cent][i_jet] = (TH1*)((TH1*)nom_file->Get(name.c_str()))->Clone(Form("%s_nom", name.c_str()));
					h_sys[i_sys][i_trk][i_cent][i_jet] = (TH1*)((TH1*)sys_files[i_sys]->Get(name.c_str()))->Clone(Form("%s_%s", name.c_str(), sys_names[i_sys].c_str()));
					h_sys[i_sys][i_trk][i_cent][i_jet]->Add(h_nom[i_trk][i_cent][i_jet], -1);
					h_sys[i_sys][i_trk][i_cent][i_jet]->Divide(h_nom[i_trk][i_cent][i_jet]);

//					cout << Form("%i_%i_%i %i ", i_trk, i_cent, i_jet, i_sys);
//					cout << Form("%s --> x", sys_names[i_sys].c_str()) << endl;

					name = Form("h_%s_final_indR_trk%i_cent%i_jetpt%i",mode.c_str(), i_trk, i_cent, i_jet);
					h_sys_p[i_sys][i_trk][i_cent][i_jet] = (TH1*)h_sys[i_sys][i_trk][i_cent][i_jet]->Clone(Form("%s_%s_p", name.c_str(), sys_names[i_sys].c_str()) );
					h_sys_p[i_sys][i_trk][i_cent][i_jet]->Reset();

					h_sys_n[i_sys][i_trk][i_cent][i_jet] = (TH1*)h_sys[i_sys][i_trk][i_cent][i_jet]->Clone(Form("%s_%s_n", name.c_str(), sys_names[i_sys].c_str()) );
					h_sys_n[i_sys][i_trk][i_cent][i_jet]->Reset();

					//setup combined and total histo
					if (i_sys == 0)
					{
						name = Form("h_%s_final_indR_trk%i_cent%i_jetpt%i_total_sys_pos", mode.c_str(), i_trk, i_cent, i_jet);
						h_total_sys_p[i_trk][i_cent][i_jet] = (TH1*)h_sys[i_sys][i_trk][i_cent][i_jet]->Clone(name.c_str());
						h_total_sys_p[i_trk][i_cent][i_jet]->Reset();

						name = Form("h_%s_final_indR_trk%i_cent%i_jetpt%i_total_sys_neg", mode.c_str(), i_trk, i_cent, i_jet);
						h_total_sys_n[i_trk][i_cent][i_jet] = (TH1*)h_sys[i_sys][i_trk][i_cent][i_jet]->Clone(name.c_str());
						h_total_sys_n[i_trk][i_cent][i_jet]->Reset();

						for (int i_comb_sys = 0; i_comb_sys < combined_sys_names.size(); i_comb_sys++)
						{
							name = Form("h_%s_final_indR_trk%i_cent%i_jetpt%i_comb_sys_pos",mode.c_str(), i_trk, i_cent, i_jet);
							h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet] = (TH1*)h_sys[i_sys][i_trk][i_cent][i_jet]->Clone(name.c_str());
							h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->Reset();

							name = Form("h_%s_final_indR_trk%i_cent%i_jetpt%i_comb_sys_neg",mode.c_str(), i_trk, i_cent, i_jet);
							h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet] = (TH1*)h_sys[i_sys][i_trk][i_cent][i_jet]->Clone(name.c_str());
							h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->Reset();
						}
					}


					for (int i_dR = 1; i_dR <= N_dR; i_dR++)
					{
						double tmp = h_sys[i_sys][i_trk][i_cent][i_jet]->GetBinContent(i_dR);

						//special rules
						//symmetrize for JER, UE
						if (sys_names[i_sys] == "sys101" || sys_names[i_sys] == "sys200" || sys_names[i_sys] == "sys201")
						{
							h_sys_p[i_sys][i_trk][i_cent][i_jet]->SetBinContent(i_dR, fabs(tmp));
							h_sys_n[i_sys][i_trk][i_cent][i_jet]->SetBinContent(i_dR, -fabs(tmp));
						}
						else //regular rules, setup pos and negative sides
						{
							if (tmp >= 0. ) h_sys_p[i_sys][i_trk][i_cent][i_jet]->SetBinContent(i_dR, tmp);
							else h_sys_n[i_sys][i_trk][i_cent][i_jet]->SetBinContent(i_dR, tmp);
						}

					}
					output_file->cd();

					name = Form("h_%s_sys_trk%i_cent%i_jetpt%i_%s", mode.c_str(), i_trk, i_cent, i_jet, sys_names[i_sys].c_str());
					h_sys[i_sys][i_trk][i_cent][i_jet]->Write(name.c_str());

					name = Form("h_%s_sys_trk%i_cent%i_jetpt%i_%s_p", mode.c_str(), i_trk, i_cent, i_jet, sys_names[i_sys].c_str());
					h_sys_p[i_sys][i_trk][i_cent][i_jet]->Write(name.c_str());

					name = Form("h_%s_sys_trk%i_cent%i_jetpt%i_%s_n",mode.c_str(), i_trk, i_cent, i_jet, sys_names[i_sys].c_str());
					h_sys_n[i_sys][i_trk][i_cent][i_jet]->Write(name.c_str());
				}
			}
		}
	}


	for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
	{
		for (int i_cent = 0; i_cent < n_cent_cuts; i_cent++)
		{
			if ((dataset_type == "_PbPb" || mode == "RDpT") && i_cent == 6) continue;
			if (dataset_type == "_pp" && i_cent < 6) continue;

			for (int i_trk = 0; i_trk < N_trkpt; i_trk++)
			{


				for (int i_comb_sys = 0; i_comb_sys < combined_sys_names.size(); i_comb_sys++)
				{
					for (int i_dR = 1; i_dR <= N_dR; i_dR++)
					{
						for (int i_sys = 0; i_sys < sys_names.size(); i_sys++)
						{

//							cout << Form("%i_%i_%i_%i %i_%i ", i_trk, i_cent, i_jet, i_dR, i_sys, i_comb_sys);
//							cout << Form("%s --> %s", sys_names[i_sys].c_str(), combined_sys_names[i_comb_sys].c_str()) << endl;

							double tmp;

							//each systematic is by itself
							if (

								(sys_names[i_sys] == "sys101" && combined_sys_names[i_comb_sys] == "JER") ||
								(sys_names[i_sys] == "sys103" && combined_sys_names[i_comb_sys] == "Unfolding") ||
								(sys_names[i_sys] == "sys200" && combined_sys_names[i_comb_sys] == "UE")
								)
							{
//								cout << Form("%i_%i_%i_%i %i_%i ", i_trk, i_cent, i_jet, i_dR, i_sys, i_comb_sys);
//								cout << Form("%s --> %s", sys_names[i_sys].c_str(), combined_sys_names[i_comb_sys].c_str()) << endl;

								tmp = h_sys_p[i_sys][i_trk][i_cent][i_jet]->GetBinContent(i_dR);
								h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->SetBinContent(i_dR, tmp);

								tmp = h_sys_n[i_sys][i_trk][i_cent][i_jet]->GetBinContent(i_dR);
								h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->SetBinContent(i_dR, tmp);
							}

							//each systematic needs to be combined
							if (

								(
								 (sys_names[i_sys] == "sys106" ||
								  sys_names[i_sys] == "sys107" ||
								  sys_names[i_sys] == "sys108" ||
								  sys_names[i_sys] == "sys109" ||
								  sys_names[i_sys] == "sys116" ||
								  sys_names[i_sys] == "sys117" ||
								  sys_names[i_sys] == "sys118" ||
								  sys_names[i_sys] == "sys119" ) &&
								 combined_sys_names[i_comb_sys] == "JES"
								 )

								||

								(
								 (sys_names[i_sys] == "sys102" ||
								  sys_names[i_sys] == "sys105" ||
								  sys_names[i_sys] == "sys110" ||
								  sys_names[i_sys] == "sys111" ||
								  sys_names[i_sys] == "sys114" ||
								  sys_names[i_sys] == "sys115" ||
								  sys_names[i_sys] == "sys201" ) &&
								 combined_sys_names[i_comb_sys] == "Tracking"
								 )

								)
							{
//								cout << Form("%i_%i_%i_%i %i_%i ", i_trk, i_cent, i_jet, i_dR, i_sys, i_comb_sys);
//								cout << Form("%s --> %s", sys_names[i_sys].c_str(), combined_sys_names[i_comb_sys].c_str()) << endl;

								if (sys_names[i_sys] == "sys106" ||
									sys_names[i_sys] == "sys102")
								{
									tmp = pow(h_sys_p[i_sys][i_trk][i_cent][i_jet]->GetBinContent(i_dR),2);
									h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->SetBinContent(i_dR, sqrt(tmp) );

									tmp = pow(h_sys_n[i_sys][i_trk][i_cent][i_jet]->GetBinContent(i_dR),2);
									h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->SetBinContent(i_dR, -sqrt(tmp) );
								}
								else
								{
									tmp = pow(h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->GetBinContent(i_dR),2);
									tmp = tmp + pow(h_sys_p[i_sys][i_trk][i_cent][i_jet]->GetBinContent(i_dR),2);
									h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->SetBinContent(i_dR, sqrt(tmp) );

									tmp = pow(h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->GetBinContent(i_dR),2);
									tmp = tmp + pow(h_sys_n[i_sys][i_trk][i_cent][i_jet]->GetBinContent(i_dR),2);
									h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->SetBinContent(i_dR, -sqrt(tmp) );
								}

							}

						}
					}
					name = Form("h_%s_sys_trk%i_cent%i_jetpt%i_%s_p",mode.c_str(), i_trk, i_cent, i_jet,combined_sys_names[i_comb_sys].c_str());
					h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->Write(name.c_str());

					name = Form("h_%s_sys_trk%i_cent%i_jetpt%i_%s_n",mode.c_str(), i_trk, i_cent, i_jet,combined_sys_names[i_comb_sys].c_str());
					h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->Write(name.c_str());

				}
			}
		}
	}


	for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
	{
		for (int i_cent = 0; i_cent < n_cent_cuts; i_cent++)
		{
			if ((dataset_type == "_PbPb" || mode == "RDpT") && i_cent == 6) continue;
			if (dataset_type == "_pp" && i_cent < 6) continue;

			for (int i_trk = 0; i_trk < N_trkpt; i_trk++)
			{
				for (int i_dR = 1; i_dR <= N_dR; i_dR++)
				{
					double tmp;
					for (int i_comb_sys = 0; i_comb_sys < combined_sys_names.size(); i_comb_sys++)
					{

						if (i_comb_sys == 0)
						{
							tmp = pow(h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->GetBinContent(i_dR),2);
							h_total_sys_p[i_trk][i_cent][i_jet]->SetBinContent(i_dR, sqrt(tmp) );

							tmp = pow(h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->GetBinContent(i_dR),2);
							h_total_sys_n[i_trk][i_cent][i_jet]->SetBinContent(i_dR, -sqrt(tmp) );
						}
						else
						{
							tmp = pow(h_total_sys_p[i_trk][i_cent][i_jet]->GetBinContent(i_dR),2);
							tmp = tmp + pow(h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->GetBinContent(i_dR),2);
							h_total_sys_p[i_trk][i_cent][i_jet]->SetBinContent(i_dR, sqrt(tmp) );

							tmp = pow(h_total_sys_n[i_trk][i_cent][i_jet]->GetBinContent(i_dR),2);
							tmp = tmp + pow(h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->GetBinContent(i_dR),2);
							h_total_sys_n[i_trk][i_cent][i_jet]->SetBinContent(i_dR, -sqrt(tmp) );
						}


					}
				}

				name = Form("h_%s_sys_trk%i_cent%i_jetpt%i_total_p",mode.c_str(), i_trk, i_cent, i_jet);
				h_total_sys_p[i_trk][i_cent][i_jet]->Write(name.c_str());

				name = Form("h_%s_sys_trk%i_cent%i_jetpt%i_total_n",mode.c_str(), i_trk, i_cent, i_jet);
				h_total_sys_n[i_trk][i_cent][i_jet]->Write(name.c_str());

			}
		}
	}


	//drawing

	string rdptr_label = "#it{R}_{ #it{D} (#it{p}_{T}, #it{r})}";
	if (mode == "ChPS") rdptr_label = "#it{D} (#it{p}_{T}, #it{r})";
	string r_label = "#it{r}";

	TCanvas *c_sys = new TCanvas("c_sys","c_sys", 1200, 600);
	TLegend *legend_sys = new TLegend(0.18, 0.18, 0.30, 0.45, "","brNDC");
	legend_sys->SetTextFont(43);
	legend_sys->SetBorderSize(0);
	legend_sys->SetTextSize(14);

	TLatex *ltx = new TLatex();
	ltx->SetTextFont(43);
	ltx->SetTextSize(12);

	for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
	{
		string jet_label = Form("%1.0f < #it{p}_{T}^{jet} < %1.0f GeV", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));

		for (int i_trk = trk_pt_start; i_trk < trk_pt_end; i_trk++)
		{
			string trk_label = Form("%1.1f < #it{p}_{T}^{trk} < %1.1f GeV", trkpT_binning->GetBinLowEdge(i_trk+1), trkpT_binning->GetBinUpEdge(i_trk+1));

			legend_sys->Clear();
			c_sys->cd();
			c_sys->Clear();
			c_sys->Divide(3,2);

			bool cent_first_pass = true;
			for (int i_cent = 0; i_cent < n_cent_cuts; i_cent++)
			{
				if ((dataset_type == "_PbPb" || mode == "RDpT") && i_cent == 6) continue;
				if (dataset_type == "_pp" && i_cent < 6) continue;

				string centrality = num_to_cent(31,i_cent);

				SetHStyle_smallify(h_total_sys_p[i_trk][i_cent][i_jet],0, doSmall);
				SetHStyle_smallify(h_total_sys_n[i_trk][i_cent][i_jet],0, doSmall);
				h_total_sys_p[i_trk][i_cent][i_jet]->SetLineWidth(1);
				h_total_sys_n[i_trk][i_cent][i_jet]->SetLineWidth(1);


				h_total_sys_p[i_trk][i_cent][i_jet]->GetYaxis()->SetRangeUser(-0.3,0.3);
				h_total_sys_n[i_trk][i_cent][i_jet]->GetYaxis()->SetRangeUser(-0.3,0.3);
				h_total_sys_p[i_trk][i_cent][i_jet]->GetXaxis()->SetRangeUser(0, 0.6);
				h_total_sys_n[i_trk][i_cent][i_jet]->GetXaxis()->SetRangeUser(0, 0.6);
				h_total_sys_p[i_trk][i_cent][i_jet]->GetYaxis()->SetTitle(Form("#delta%s",rdptr_label.c_str()));
				h_total_sys_n[i_trk][i_cent][i_jet]->GetYaxis()->SetTitle(Form("#delta%s",rdptr_label.c_str()));
				h_total_sys_p[i_trk][i_cent][i_jet]->GetXaxis()->SetTitle(r_label.c_str());
				h_total_sys_n[i_trk][i_cent][i_jet]->GetXaxis()->SetTitle(r_label.c_str());


				c_sys->cd(i_cent+1);
				h_total_sys_p[i_trk][i_cent][i_jet]->Draw("l");
				h_total_sys_n[i_trk][i_cent][i_jet]->Draw("l same");


				for (int i_comb_sys = 0; i_comb_sys < combined_sys_names.size(); i_comb_sys++)
				{
					SetHStyle_smallify(h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet],i_comb_sys+1, doSmall);
					SetHStyle_smallify(h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet],i_comb_sys+1, doSmall);
					h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->SetLineWidth(1);
					h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->SetLineWidth(1);

					bool empty_hist_p = false;
					bool empty_hist_n = false;
					if (h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->GetMaximum() == 0 && h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->GetMinimum() == 0) empty_hist_p = true;
					if (h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->GetMaximum() == 0 && h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->GetMinimum() == 0) empty_hist_n = true;

					if (cent_first_pass) legend_sys->AddEntry(h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet],combined_sys_names[i_comb_sys].c_str(),"l");

					h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->GetYaxis()->SetRangeUser(-0.3,0.3);
					h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->GetYaxis()->SetRangeUser(-0.3,0.3);
					h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->GetXaxis()->SetRangeUser(0, 0.6);
					h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->GetXaxis()->SetRangeUser(0, 0.6);
					h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->GetYaxis()->SetTitle(Form("#delta%s",rdptr_label.c_str()));
					h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->GetYaxis()->SetTitle(Form("#delta%s",rdptr_label.c_str()));
					h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->GetXaxis()->SetTitle(r_label.c_str());
					h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->GetXaxis()->SetTitle(r_label.c_str());

					c_sys->cd(i_cent+1);
					if (!empty_hist_p) h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->Draw("l same");
					if (!empty_hist_n) h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->Draw("l same");
				}



				c_sys->cd(i_cent+1);
				if (dataset_type == "_pp") c_sys->cd();
				if (cent_first_pass) legend_sys->AddEntry(h_total_sys_p[i_trk][i_cent][i_jet],"Total","l");

				legend_sys->Draw();
				ltx->SetTextAlign(32);
				ltx->SetTextSize(12);

				ltx->DrawLatexNDC(0.93, 0.90, Form("%s", trk_label.c_str()));
				ltx->DrawLatexNDC(0.93, 0.85, Form("%s", jet_label.c_str()));
				ltx->DrawLatexNDC(0.93, 0.80, Form("%s", centrality.c_str()));

				cent_first_pass = false;


			}


			ltx->SetTextAlign(12);


			c_sys->cd(1);
			ATLASLabel(0.19, 0.88, "     Internal", "", kBlack);
			ltx->SetTextAlign(12);
			ltx->SetTextSize(10);

			if (mode == "RDpT")
			{
				ltx->DrawLatexNDC(0.2, 0.8, "#splitline{Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}}{#it{pp}  #sqrt{#font[12]{s}} = 5.02 TeV, 25 pb^{-1}}");
			}
			else
			{
				if (dataset_type == "_PbPb") ltx->DrawLatexNDC(0.2, 0.75, "Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}");
				if (dataset_type == "_pp") ltx->DrawLatexNDC(0.2, 0.8, "#it{pp}  #sqrt{#font[12]{s}} = 5.02 TeV, 25 pb^{-1}");
			}

			string pdf_label = "";
			if (i_trk == trk_pt_start && i_jet == jet_pt_start) pdf_label = "(";
			if (i_trk == trk_pt_end-1 && i_jet == jet_pt_end-1) pdf_label = ")";
			c_sys->Print(Form("output_pdf_nominal/systematics/%s_dR_sys%s_error.pdf%s",mode.c_str(), dataset_type.c_str(), pdf_label.c_str()), Form("Title:trk%i_jetpt%i", i_trk, i_jet));


		}
	}



	cout << "######### Done Systematics #########" << endl;

}




