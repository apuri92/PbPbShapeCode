#include "../functions/global_variables.h"
#include "draw_functions.c"
#include "TEnv.h"

void systematics(string config_file = "sys_config.cfg")
{
	cout << "######### DOING Systematics #########" << endl;
	SetAtlasStyle();
	gErrorIgnoreLevel = 3001;

	TEnv *m_config = new TEnv();
	m_config->ReadFile(config_file.c_str(), EEnvLevel(1));

	std::string dataset_type = "PbPb"; dataset_type = m_config->GetValue("dataset_type", dataset_type.c_str());
	std::string mode = "RDpT"; mode = m_config->GetValue("mode", mode.c_str());
	int centrality_scheme = 31; centrality_scheme = m_config->GetValue("centrality_scheme", centrality_scheme);
	int isMC = 1; isMC = m_config->GetValue("isMC", isMC);
	int verbose = 0; verbose = m_config->GetValue("verbose", verbose);

	std::string did = "data";
	if (isMC) did = "MC";



	string name = "";

	TFile* nom_file = new TFile(Form("output_pdf_nominal/root/final_%s_%s.root", mode.c_str(), did.c_str()));
	TFile* output_file = new TFile(Form("output_pdf_nominal/root/final_%s_sys_%s.root", mode.c_str(), did.c_str()), "recreate");

	TAxis* dR_binning = (TAxis*)nom_file->Get("dR_binning");
	TAxis* jetpT_binning = (TAxis*)nom_file->Get("jetpT_binning");
	TAxis* trkpT_binning = (TAxis*)nom_file->Get("trkpT_binning");

	int N_dR = dR_binning->GetNbins();
	int N_jetpt = jetpT_binning->GetNbins();
	int N_trkpt = trkpT_binning->GetNbins();


	vector<TFile*> sys_files;
	vector<string> sys_names;
	vector<string> sys_names;

	sys_names.push_back("sys101"); //JER
	sys_names.push_back("sys102"); //Sign
//	sys_names.push_back("sys105"); //MCProbCut
	sys_names.push_back("sys106"); //HIJES_1_P
	sys_names.push_back("sys107"); //HIJES_2_P
	sys_names.push_back("sys108"); //HIJES_1_N
//	sys_names.push_back("sys109"); //HIJES_2_N
//	sys_names.push_back("sys110"); //Material_P
//	sys_names.push_back("sys111"); //Material_N
//	sys_names.push_back("sys114"); //Tracking
//	sys_names.push_back("sys115"); //TrackingRes
//	sys_names.push_back("sys116"); //CentHIJES_P
//	sys_names.push_back("sys117"); //CentHIJES_N
//	sys_names.push_back("sys118"); //ppJES_P
//	sys_names.push_back("sys199"); //ppJES_N
	sys_names.push_back("sys200"); //UE


	combined_sys_names.push_back("JER");
	combined_sys_names.push_back("Sign");
//	combined_sys_names.push_back("MCProb");
	combined_sys_names.push_back("HIJES");
//	combined_sys_names.push_back("Material");
//	combined_sys_names.push_back("Tracking");
//	combined_sys_names.push_back("TrackingResolution");
//	combined_sys_names.push_back("CentHIJES");
//	combined_sys_names.push_back("ppJES");
	combined_sys_names.push_back("UE");


	//individual systematics
	vector<vector<vector<vector<TH1*>>>> h_sys (sys_names.size(), vector<vector<vector<TH1*>>> (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt))));

	//combined totals
	vector<vector<vector<vector<TH1*>>>> h_comb_sys (combined_sys_names.size(), vector<vector<vector<TH1*>>> (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt))));
	vector<vector<vector<vector<vector<TH1*>>>>> h_comb_sys_total (2, vector<vector<vector<vector<TH1*>>>> (combined_sys_names.size(), vector<vector<vector<TH1*>>> (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)))));

	//individual totals
	vector<vector<vector<vector<TH1*>>>> h_sys_total_p (sys_names.size(), vector<vector<vector<TH1*>>> (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt))));
	vector<vector<vector<vector<TH1*>>>> h_sys_total_n (sys_names.size(), vector<vector<vector<TH1*>>> (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt))));
	vector<vector<vector<vector<TH1*>>>> h_sys_total (sys_names.size(), vector<vector<vector<TH1*>>> (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt))));



	vector<vector<vector<TH1*>>> h_sys_JERpos (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_sys_JERneg (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));

	vector<vector<vector<TH1*>>> h_sys_Signpos (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_sys_Signneg (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));

	vector<vector<vector<TH1*>>> h_sys_JESpos (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_sys_JESneg (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));

	vector<vector<vector<TH1*>>> h_sys_UEpos (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_sys_UEneg (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));


	vector<vector<vector<TH1*>>> h_sys_Totalpos (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_sys_Totalneg (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));


	int jet_pt_start = 7;
	int jet_pt_end = 11;
	int trk_pt_start = 2;
	int trk_pt_end = 9;

	for (int i = 0; i < sys_names.size(); i++)
	{
		name = Form("output_pdf_%s/root/final_RDpT_%s.root", sys_names[i].c_str(), did.c_str());
		sys_files.push_back( new TFile( name.c_str() ) );

		for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
		{
			for (int i_cent = 0; i_cent < 6; i_cent++)
			{
				for (int i_trk = 0; i_trk < N_trkpt; i_trk++)
				{

					name = Form("h_ChPS_RDpT_indR_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);

					TH1* h_nom = (TH1*)nom_file->Get(name.c_str())->Clone(Form("%s_nom", name.c_str() ));
					TH1* h_sys = (TH1*)sys_files[i]->Get(name.c_str())->Clone(Form("%s_%s", name.c_str(), sys_names[i].c_str() ));


					h_sys[i][i_trk][i_cent][i_jet] = (TH1*)h_sys->Clone(Form("%s_err", name.c_str()));
					h_sys[i][i_trk][i_cent][i_jet]->Add(h_nom, -1);
					h_sys[i][i_trk][i_cent][i_jet]->Divide(h_nom);

					h_sys[i][i_trk][i_cent][i_jet]->GetYaxis()->SetRangeUser(-1,1);


					if (i == 0)
					{

						//JER
						name = Form("h_ChPS_RDpT_indR_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
						h_sys_JERpos[i_trk][i_cent][i_jet] = (TH1*)h_sys[i][i_trk][i_cent][i_jet]->Clone(Form("%s_sys_JERpos", name.c_str()));
						h_sys_JERpos[i_trk][i_cent][i_jet]->Reset();

						h_sys_JERneg[i_trk][i_cent][i_jet] = (TH1*)h_sys[i][i_trk][i_cent][i_jet]->Clone(Form("%s_sys_JERneg", name.c_str()));
						h_sys_JERneg[i_trk][i_cent][i_jet]->Reset();


						//Sign
						name = Form("h_ChPS_RDpT_indR_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
						h_sys_Signpos[i_trk][i_cent][i_jet] = (TH1*)h_sys[i][i_trk][i_cent][i_jet]->Clone(Form("%s_sys_Signpos", name.c_str()));
						h_sys_Signpos[i_trk][i_cent][i_jet]->Reset();

						h_sys_Signneg[i_trk][i_cent][i_jet] = (TH1*)h_sys[i][i_trk][i_cent][i_jet]->Clone(Form("%s_sys_Signneg", name.c_str()));
						h_sys_Signneg[i_trk][i_cent][i_jet]->Reset();

						//JES
						name = Form("h_ChPS_RDpT_indR_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
						h_sys_JESpos[i_trk][i_cent][i_jet] = (TH1*)h_sys[i][i_trk][i_cent][i_jet]->Clone(Form("%s_sys_JESpos", name.c_str()));
						h_sys_JESpos[i_trk][i_cent][i_jet]->Reset();

						h_sys_JESneg[i_trk][i_cent][i_jet] = (TH1*)h_sys[i][i_trk][i_cent][i_jet]->Clone(Form("%s_sys_JESneg", name.c_str()));
						h_sys_JESneg[i_trk][i_cent][i_jet]->Reset();

						//UE
						name = Form("h_ChPS_RDpT_indR_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
						h_sys_UEpos[i_trk][i_cent][i_jet] = (TH1*)h_sys[i][i_trk][i_cent][i_jet]->Clone(Form("%s_sys_UE", name.c_str()));
						h_sys_UEpos[i_trk][i_cent][i_jet]->Reset();
						h_sys_UEneg[i_trk][i_cent][i_jet] = (TH1*)h_sys[i][i_trk][i_cent][i_jet]->Clone(Form("%s_sys_UE", name.c_str()));
						h_sys_UEneg[i_trk][i_cent][i_jet]->Reset();


						//Total
						name = Form("h_ChPS_RDpT_indR_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
						h_sys_Totalpos[i_trk][i_cent][i_jet] = (TH1*)h_sys[i][i_trk][i_cent][i_jet]->Clone(Form("%s_sys_Totalpos", name.c_str()));
						h_sys_Totalpos[i_trk][i_cent][i_jet]->Reset();

						h_sys_Totalneg[i_trk][i_cent][i_jet] = (TH1*)h_sys[i][i_trk][i_cent][i_jet]->Clone(Form("%s_sys_Totalneg", name.c_str()));
						h_sys_Totalneg[i_trk][i_cent][i_jet]->Reset();
					}

					output_file->cd();
					name = Form("h_ChPS_RDpT_indR_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
					h_sys[i][i_trk][i_cent][i_jet]->SetName(Form("%s_%s", name.c_str(), sys_names[i].c_str()));
					h_sys[i][i_trk][i_cent][i_jet]->SetTitle(Form("%s_%s", name.c_str(), sys_names[i].c_str()));
					h_sys[i][i_trk][i_cent][i_jet]->Write(Form("%s_%s", name.c_str(), sys_names[i].c_str()));
				}
			}
		}

	}

	double pos_JES = 0, neg_JES = 0, tmp = 0;

	for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
	{
		for (int i_cent = 0; i_cent < 6; i_cent++)
		{
			for (int i_trk = 0; i_trk < N_trkpt; i_trk++)
			{
				for (int i_dR = 1; i_dR <= N_dR; i_dR++)
				{
					for (int i = 0; i < sys_names.size(); i++)
					{


						//JER symmetrized
						if (sys_names[i] == "sys101")
						{
							tmp = h_sys[i][i_trk][i_cent][i_jet]->GetBinContent(i_dR);

							h_sys_JERpos[i_trk][i_cent][i_jet]->SetBinContent(i_dR, tmp);
							h_sys_JERneg[i_trk][i_cent][i_jet]->SetBinContent(i_dR, -tmp);
						}

						if (sys_names[i] == "sys102")
						{
							tmp = h_sys[i][i_trk][i_cent][i_jet]->GetBinContent(i_dR);

							if (tmp >= 0) h_sys_Signpos[i_trk][i_cent][i_jet]->SetBinContent(i_dR, tmp);
							else h_sys_Signneg[i_trk][i_cent][i_jet]->SetBinContent(i_dR, tmp);
						}

						//JES added up
						if (sys_names[i] == "sys106" ||
							sys_names[i] == "sys107" ||
							sys_names[i] == "sys108")

						{
							if (sys_names[i] == "sys106")
							{
								pos_JES = 0; neg_JES = 0;
							}
							tmp = pow(h_sys[i][i_trk][i_cent][i_jet]->GetBinContent(i_dR),2);

							if (tmp >= 0) pos_JES = pos_JES + tmp;
							else neg_JES = neg_JES - tmp;
						}

						//UE symmetrized
						if (sys_names[i] == "sys200")
						{
							tmp = h_sys[i][i_trk][i_cent][i_jet]->GetBinContent(i_dR);

							h_sys_UEpos[i_trk][i_cent][i_jet]->SetBinContent(i_dR, tmp);
							h_sys_UEneg[i_trk][i_cent][i_jet]->SetBinContent(i_dR, -tmp);
						}

					}

					//JES
					h_sys_JESpos[i_trk][i_cent][i_jet]->SetBinContent(i_dR, sqrt(pos_JES));
					h_sys_JESneg[i_trk][i_cent][i_jet]->SetBinContent(i_dR, -sqrt(fabs(neg_JES)));




					//totals pos
					tmp = 0;
					tmp = tmp + pow(h_sys_JESpos[i_trk][i_cent][i_jet]->GetBinContent(i_dR),2);
					tmp = tmp + pow(h_sys_Signpos[i_trk][i_cent][i_jet]->GetBinContent(i_dR),2);
					tmp = tmp + pow(h_sys_JERpos[i_trk][i_cent][i_jet]->GetBinContent(i_dR),2);
					tmp = tmp + pow(h_sys_UEpos[i_trk][i_cent][i_jet]->GetBinContent(i_dR),2);
					//add more here
					h_sys_Totalpos[i_trk][i_cent][i_jet]->SetBinContent(i_dR, sqrt(tmp));



					//total neg
					tmp = 0;
					tmp = tmp + pow(h_sys_JESneg[i_trk][i_cent][i_jet]->GetBinContent(i_dR),2);
					tmp = tmp + pow(h_sys_Signneg[i_trk][i_cent][i_jet]->GetBinContent(i_dR),2);
					tmp = tmp + pow(h_sys_JERneg[i_trk][i_cent][i_jet]->GetBinContent(i_dR),2);
					tmp = tmp + pow(h_sys_UEneg[i_trk][i_cent][i_jet]->GetBinContent(i_dR),2);
					//add more here
					h_sys_Totalneg[i_trk][i_cent][i_jet]->SetBinContent(i_dR, -sqrt(tmp));


				}

				output_file->cd();
				//JER
				name = Form("h_ChPS_RDpT_indR_trk%i_cent%i_jetpt%i_sys_JERpos", i_trk, i_cent, i_jet);
				h_sys_JERpos[i_trk][i_cent][i_jet]->SetName(name.c_str());
				h_sys_JERpos[i_trk][i_cent][i_jet]->SetTitle(name.c_str());
				h_sys_JERpos[i_trk][i_cent][i_jet]->Write(name.c_str());

				name = Form("h_ChPS_RDpT_indR_trk%i_cent%i_jetpt%i_sys_JERneg", i_trk, i_cent, i_jet);
				h_sys_JERneg[i_trk][i_cent][i_jet]->SetName(name.c_str());
				h_sys_JERneg[i_trk][i_cent][i_jet]->SetTitle(name.c_str());
				h_sys_JERneg[i_trk][i_cent][i_jet]->Write(name.c_str());

				//Sign
				name = Form("h_ChPS_RDpT_indR_trk%i_cent%i_jetpt%i_sys_Signpos", i_trk, i_cent, i_jet);
				h_sys_Signpos[i_trk][i_cent][i_jet]->SetName(name.c_str());
				h_sys_Signpos[i_trk][i_cent][i_jet]->SetTitle(name.c_str());
				h_sys_Signpos[i_trk][i_cent][i_jet]->Write(name.c_str());

				name = Form("h_ChPS_RDpT_indR_trk%i_cent%i_jetpt%i_sys_Signneg", i_trk, i_cent, i_jet);
				h_sys_Signneg[i_trk][i_cent][i_jet]->SetName(name.c_str());
				h_sys_Signneg[i_trk][i_cent][i_jet]->SetTitle(name.c_str());
				h_sys_Signneg[i_trk][i_cent][i_jet]->Write(name.c_str());

				//JES
				name = Form("h_ChPS_RDpT_indR_trk%i_cent%i_jetpt%i_sys_JESpos", i_trk, i_cent, i_jet);
				h_sys_JESpos[i_trk][i_cent][i_jet]->SetName(name.c_str());
				h_sys_JESpos[i_trk][i_cent][i_jet]->SetTitle(name.c_str());
				h_sys_JESpos[i_trk][i_cent][i_jet]->Write(name.c_str());

				name = Form("h_ChPS_RDpT_indR_trk%i_cent%i_jetpt%i_sys_JESneg", i_trk, i_cent, i_jet);
				h_sys_JESneg[i_trk][i_cent][i_jet]->SetName(name.c_str());
				h_sys_JESneg[i_trk][i_cent][i_jet]->SetTitle(name.c_str());
				h_sys_JESneg[i_trk][i_cent][i_jet]->Write(name.c_str());

				//UE
				name = Form("h_ChPS_RDpT_indR_trk%i_cent%i_jetpt%i_sys_UEpos", i_trk, i_cent, i_jet);
				h_sys_UEpos[i_trk][i_cent][i_jet]->SetName(name.c_str());
				h_sys_UEpos[i_trk][i_cent][i_jet]->SetTitle(name.c_str());
				h_sys_UEpos[i_trk][i_cent][i_jet]->Write(name.c_str());

				name = Form("h_ChPS_RDpT_indR_trk%i_cent%i_jetpt%i_sys_UEneg", i_trk, i_cent, i_jet);
				h_sys_UEneg[i_trk][i_cent][i_jet]->SetName(name.c_str());
				h_sys_UEneg[i_trk][i_cent][i_jet]->SetTitle(name.c_str());
				h_sys_UEneg[i_trk][i_cent][i_jet]->Write(name.c_str());

				//Total
				name = Form("h_ChPS_RDpT_indR_trk%i_cent%i_jetpt%i_sys_Totalpos", i_trk, i_cent, i_jet);
				h_sys_Totalpos[i_trk][i_cent][i_jet]->SetName(name.c_str());
				h_sys_Totalpos[i_trk][i_cent][i_jet]->SetTitle(name.c_str());
				h_sys_Totalpos[i_trk][i_cent][i_jet]->Write(name.c_str());

				name = Form("h_ChPS_RDpT_indR_trk%i_cent%i_jetpt%i_sys_Totalneg", i_trk, i_cent, i_jet);
				h_sys_Totalneg[i_trk][i_cent][i_jet]->SetName(name.c_str());
				h_sys_Totalneg[i_trk][i_cent][i_jet]->SetTitle(name.c_str());
				h_sys_Totalneg[i_trk][i_cent][i_jet]->Write(name.c_str());


			}
		}
	}

	//drawing
	TCanvas *c_sys = new TCanvas("c_sys","c_sys", 1200, 600);
	TLegend *legend_sys = new TLegend(0.18, 0.18, 0.30, 0.38, "","brNDC");
	legend_sys->SetTextFont(43);
	legend_sys->SetBorderSize(0);

	TLatex *ltx = new TLatex();
	ltx->SetTextFont(43);
	ltx->SetTextSize(12);

	for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
	{
		string jet_label = Form("%1.0f < p_{T}^{Jet} < %1.0f", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));

		for (int i_trk = trk_pt_start; i_trk < trk_pt_end; i_trk++)
		{
			string trk_label = Form("%1.2f < p_{T}^{Trk} < %1.2f", trkpT_binning->GetBinLowEdge(i_trk+1), trkpT_binning->GetBinUpEdge(i_trk+1));

			legend_sys->Clear();
			c_sys->cd();
			c_sys->Clear();
			c_sys->Divide(3,2);

			for (int i_cent = 0; i_cent < 6; i_cent++)
			{
				string centrality = num_to_cent(31,i_cent);

				SetHStyle_smallify(h_sys_JERpos[i_trk][i_cent][i_jet],1, 1);
				SetHStyle_smallify(h_sys_JERneg[i_trk][i_cent][i_jet],1, 1);
				SetHStyle_smallify(h_sys_Signpos[i_trk][i_cent][i_jet],2, 1);
				SetHStyle_smallify(h_sys_Signneg[i_trk][i_cent][i_jet],2, 1);
				SetHStyle_smallify(h_sys_JESpos[i_trk][i_cent][i_jet],3, 1);
				SetHStyle_smallify(h_sys_JESneg[i_trk][i_cent][i_jet],3, 1);
				SetHStyle_smallify(h_sys_UEpos[i_trk][i_cent][i_jet],4, 1);
				SetHStyle_smallify(h_sys_UEneg[i_trk][i_cent][i_jet],4, 1);
				SetHStyle_smallify(h_sys_Totalpos[i_trk][i_cent][i_jet],0, 1);
				SetHStyle_smallify(h_sys_Totalneg[i_trk][i_cent][i_jet],0, 1);

				if (i_cent == 0)
				{
					legend_sys->AddEntry(h_sys_JERpos[i_trk][i_cent][i_jet],"JER","lp");
					legend_sys->AddEntry(h_sys_JESpos[i_trk][i_cent][i_jet],"HIJES","lp");
					legend_sys->AddEntry(h_sys_Signpos[i_trk][i_cent][i_jet],"Sign.","lp");
					legend_sys->AddEntry(h_sys_UEpos[i_trk][i_cent][i_jet],"UE","lp");
					legend_sys->AddEntry(h_sys_Totalpos[i_trk][i_cent][i_jet],"Total","lp");
				}


				h_sys_JERpos[i_trk][i_cent][i_jet]->GetYaxis()->SetRangeUser(-0.3,0.3);
				h_sys_JERpos[i_trk][i_cent][i_jet]->GetXaxis()->SetRangeUser(0,0.6);

				h_sys_JERpos[i_trk][i_cent][i_jet]->GetYaxis()->SetTitle("#delta R_{D (p_{T})}");
				h_sys_JERneg[i_trk][i_cent][i_jet]->GetYaxis()->SetTitle("#delta R_{D (p_{T})}");
				h_sys_Signpos[i_trk][i_cent][i_jet]->GetYaxis()->SetTitle("#delta R_{D (p_{T})}");
				h_sys_Signneg[i_trk][i_cent][i_jet]->GetYaxis()->SetTitle("#delta R_{D (p_{T})}");
				h_sys_JESpos[i_trk][i_cent][i_jet]->GetYaxis()->SetTitle("#delta R_{D (p_{T})}");
				h_sys_JESneg[i_trk][i_cent][i_jet]->GetYaxis()->SetTitle("#delta R_{D (p_{T})}");
				h_sys_UEpos[i_trk][i_cent][i_jet]->GetYaxis()->SetTitle("#delta R_{D (p_{T})}");
				h_sys_UEneg[i_trk][i_cent][i_jet]->GetYaxis()->SetTitle("#delta R_{D (p_{T})}");

				h_sys_Totalpos[i_trk][i_cent][i_jet]->GetYaxis()->SetTitle("#delta R_{D (p_{T})}");
				h_sys_Totalneg[i_trk][i_cent][i_jet]->GetYaxis()->SetTitle("#delta R_{D (p_{T})}");



				c_sys->cd(i_cent+1);
				h_sys_JERpos[i_trk][i_cent][i_jet]->Draw("l");
				h_sys_JERneg[i_trk][i_cent][i_jet]->Draw("l same");
				h_sys_Signpos[i_trk][i_cent][i_jet]->Draw("l same");
				h_sys_Signneg[i_trk][i_cent][i_jet]->Draw("l same");
				h_sys_JESpos[i_trk][i_cent][i_jet]->Draw("l same");
				h_sys_JESneg[i_trk][i_cent][i_jet]->Draw("l same");
				h_sys_UEpos[i_trk][i_cent][i_jet]->Draw("l same");
				h_sys_UEneg[i_trk][i_cent][i_jet]->Draw("l same");

				h_sys_Totalpos[i_trk][i_cent][i_jet]->Draw("l same");
				h_sys_Totalneg[i_trk][i_cent][i_jet]->Draw("l same");


//				h_sys_JERpos[i_trk][i_cent][i_jet]->Draw("ep");
//				h_sys_JERneg[i_trk][i_cent][i_jet]->Draw("ep same");
//				h_sys_Signpos[i_trk][i_cent][i_jet]->Draw("ep same");
//				h_sys_Signneg[i_trk][i_cent][i_jet]->Draw("ep same");
//				h_sys_JESpos[i_trk][i_cent][i_jet]->Draw("ep same");
//				h_sys_JESneg[i_trk][i_cent][i_jet]->Draw("ep same");
//				h_sys_UEpos[i_trk][i_cent][i_jet]->Draw("ep same");
//				h_sys_UEneg[i_trk][i_cent][i_jet]->Draw("ep same");
//
//				h_sys_Totalpos[i_trk][i_cent][i_jet]->Draw("ep same");
//				h_sys_Totalneg[i_trk][i_cent][i_jet]->Draw("ep same");


				legend_sys->Draw();
				ltx->SetTextAlign(32);
				ltx->SetTextSize(12);

				ltx->DrawLatexNDC(0.93, 0.90, Form("%s", trk_label.c_str()));
				ltx->DrawLatexNDC(0.93, 0.85, Form("%s", jet_label.c_str()));
				ltx->DrawLatexNDC(0.93, 0.80, Form("%s", centrality.c_str()));

			}


			string pdf_label = "";
			if (i_trk == trk_pt_start && i_jet == jet_pt_start) pdf_label = "(";
			if (i_trk == trk_pt_end-1 && i_jet == jet_pt_end-1) pdf_label = ")";
			c_sys->Print(Form("output_pdf_nominal/RDpT_dR_sys_error.pdf%s", pdf_label.c_str()), Form("Title:trk%i_jetpt%i", i_trk, i_jet));


		}
	}



	cout << "######### Done Systematics #########" << endl;



}
