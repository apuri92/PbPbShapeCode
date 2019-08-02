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

	double r_max_range = 0.8;
	if (mode.compare("ChPS") !=0 ) dataset_type = "";
	else dataset_type = Form("_%s", dataset_type.c_str());

	TFile* nom_file = new TFile(Form("output_pdf_nominal/root/final_%s_%s%s.root", mode.c_str(), did.c_str(), dataset_type.c_str()));
	if (mode == "RDpT" || mode == "DeltaDpT") nom_file = new TFile(Form("output_pdf_nominal/root/final_%s_%s%s.root", "RDpT", did.c_str(), dataset_type.c_str())); //get R file because difference is saved there
	TFile* output_file = new TFile(Form("output_pdf_nominal/root/final_%s_sys_%s%s.root", mode.c_str(), did.c_str(), dataset_type.c_str()), "recreate");

	//For compiledCode
	//	TFile* nom_file = new TFile(Form("compiledCode/root/output_nominal_data.root"));
	//	TFile* output_file = new TFile(Form("compiledCode/root/final_%s_sys_%s%s.root", mode.c_str(), did.c_str(), dataset_type.c_str()), "recreate");

	TFile *f_pbpb, *f_pp;

	if (dataset_type.compare("") == 0)
	{
		f_pbpb = new TFile(Form("output_pdf_nominal/root/final_ChPS_sys_data_PbPb.root"));
		f_pp = new TFile(Form("output_pdf_nominal/root/final_ChPS_sys_data_pp.root"));
		//For compiledCode
//		f_pbpb = new TFile(Form("compiledCode/root/final_ChPS_sys_data_PbPb.root"));
//		f_pp = new TFile(Form("compiledCode/root/final_ChPS_sys_data_pp.root"));
	}


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

	sys_names.push_back("sys1"); //JER
	sys_names.push_back("sys2"); //Sign
	sys_names.push_back("sys3"); //Unfolding
	sys_names.push_back("sys4"); //MC_NonClosure
	sys_names.push_back("sys5"); //MCProbCut
	sys_names.push_back("sys6"); //HIJES_1_P
	sys_names.push_back("sys7"); //HIJES_1_N
	sys_names.push_back("sys8"); //HIJES_2_P
	sys_names.push_back("sys9"); //HIJES_2_N
	sys_names.push_back("sys10"); //Material_P
	sys_names.push_back("sys11"); //Material_N
	sys_names.push_back("sys14"); //DenseEnv
	sys_names.push_back("sys15"); //TrackingRes
	sys_names.push_back("sys16"); //CentHIJES_P
	sys_names.push_back("sys17"); //CentHIJES_N
	sys_names.push_back("sys18"); //ppJES_P_101
	sys_names.push_back("sys19"); //ppJES_N_101
	sys_names.push_back("sys20"); //ppJES_P_102
	sys_names.push_back("sys21"); //ppJES_N_102
	sys_names.push_back("sys22"); //ppJES_P_103
	sys_names.push_back("sys23"); //ppJES_N_103
	sys_names.push_back("sys42"); //UE Map stat
	sys_names.push_back("sys44"); //UE Cone method
	sys_names.push_back("sys45"); //Fakes

	map<std::string, std::string> sys_i_name = {{ "sys1" , "JER" }, { "sys2" , "Significance" }, { "sys3" , "Unfolding" }, { "sys4" , "MC NonClosure" }, { "sys5" , "MCProbCut" }, { "sys6" , "HIJES_1_P" }, { "sys7" , "HIJES_1_N" }, { "sys8" , "HIJES_2_P" }, { "sys9" , "HIJES_2_N" }, { "sys10" , "Material_P" }, { "sys11" , "Material_N" }, { "sys14" , "DenseEnv" }, { "sys15" , "TrackingRes" }, { "sys16" , "CentHIJES_P" }, { "sys17" , "CentHIJES_N" }, { "sys18" , "ppJES_P_101" }, { "sys19" , "ppJES_N_101" }, { "sys20" , "ppJES_P_102" },{ "sys21" , "ppJES_N_102" }, { "sys22" , "ppJES_P_103" }, { "sys23" , "ppJES_P_103" }, { "sys42" , "UE_{MapStat}" }, { "sys43" , "UE_{RunDep}" }, { "sys44" , "UE_{Cone}" }, { "sys45" , "Fakes" }};


	combined_sys_names.push_back("JES");
	combined_sys_names.push_back("JER");
	combined_sys_names.push_back("Tracking");
	combined_sys_names.push_back("Unfolding");
	combined_sys_names.push_back("MCNonClosure");
	combined_sys_names.push_back("UE");


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
		if (mode == "RDpT" || mode == "DeltaDpT") name = Form("output_pdf_%s/root/final_%s_%s%s.root", sys_names[i_sys].c_str(), "RDpT", did.c_str(), dataset_type.c_str()); //get R file because difference is saved there

		//For compiledCode
//		name = Form("compiledCode/root/output_%s_data.root", sys_names[i_sys].c_str());
		sys_files.push_back( new TFile( name.c_str() ) );
		cout << sys_files[i_sys]->GetName() << endl;

		for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
		{
			for (int i_cent = 0; i_cent < n_cent_cuts; i_cent++)
			{
				if ((dataset_type.compare("_pp") !=0 ) && i_cent == 6) continue;
				else if (dataset_type == "_pp" && i_cent < 6) continue;

				for (int i_trk = 0; i_trk < N_trkpt; i_trk++)
				{

					//get relative errors for all modes (pbpb, pp, RDpT, DeltaDpT)
					name = Form("h_%s_final_indR_trk%i_cent%i_jetpt%i",mode.c_str(), i_trk, i_cent, i_jet);
					//For compiledCode
//					name = Form("h%s_%s_final_indR_trk%i_cent%i_jetpt%i",dataset_type.c_str(), mode.c_str(), i_trk, i_cent, i_jet);
					if (i_sys == 0)
					{
						h_nom[i_trk][i_cent][i_jet] = (TH1*)nom_file->Get(name.c_str());
						h_nom[i_trk][i_cent][i_jet]->SetName(Form("%s_nom", name.c_str()));
					}
					h_sys[i_sys][i_trk][i_cent][i_jet] = (TH1*)sys_files[i_sys]->Get(name.c_str());
					h_sys[i_sys][i_trk][i_cent][i_jet]->SetName(Form("%s_%s", name.c_str(), sys_names[i_sys].c_str()));
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
						//symmetrize for JER (1), fakes (45), UE_ConeMethod (44), UE map stat (42), unfolding (3), nonclosure (4), significance (2)
						if (sys_names[i_sys] == "sys1" ||
							sys_names[i_sys] == "sys2" ||
							sys_names[i_sys] == "sys3" ||
							sys_names[i_sys] == "sys4" ||
							sys_names[i_sys] == "sys42" ||
							sys_names[i_sys] == "sys44" ||
							sys_names[i_sys] == "sys45" )
						{
							h_sys_p[i_sys][i_trk][i_cent][i_jet]->SetBinContent(i_dR, fabs(tmp));
							h_sys_n[i_sys][i_trk][i_cent][i_jet]->SetBinContent(i_dR, -fabs(tmp));
							h_sys_p[i_sys][i_trk][i_cent][i_jet]->SetBinError(i_dR, 0.0000001);
							h_sys_n[i_sys][i_trk][i_cent][i_jet]->SetBinError(i_dR, 0.0000001);
						}
						else //regular rules, setup pos and negative sides
						{
							if (tmp >= 0. ) h_sys_p[i_sys][i_trk][i_cent][i_jet]->SetBinContent(i_dR, tmp);
							else h_sys_n[i_sys][i_trk][i_cent][i_jet]->SetBinContent(i_dR, tmp);
							h_sys_p[i_sys][i_trk][i_cent][i_jet]->SetBinError(i_dR, 0.0000001);
							h_sys_n[i_sys][i_trk][i_cent][i_jet]->SetBinError(i_dR, 0.0000001);
						}

					}
					output_file->cd();

					name = Form("h_%s_sys_trk%i_cent%i_jetpt%i_%s", mode.c_str(), i_trk, i_cent, i_jet, sys_names[i_sys].c_str());
					h_sys[i_sys][i_trk][i_cent][i_jet]->Write(name.c_str());

					name = Form("h_%s_sys_trk%i_cent%i_jetpt%i_%s_p", mode.c_str(), i_trk, i_cent, i_jet, sys_names[i_sys].c_str());
					h_sys_p[i_sys][i_trk][i_cent][i_jet]->Write(name.c_str());

					name = Form("h_%s_sys_trk%i_cent%i_jetpt%i_%s_n",mode.c_str(), i_trk, i_cent, i_jet, sys_names[i_sys].c_str());
					h_sys_n[i_sys][i_trk][i_cent][i_jet]->Write(name.c_str());

					if (i_sys == 0)
					{
						name = Form("h_%s_final_indR_trk%i_cent%i_jetpt%i",mode.c_str(), i_trk, i_cent, i_jet);
						h_nom[i_trk][i_cent][i_jet]->Write(name.c_str());
					}
				}
			}
		}
	}

	cout << "Grouping Systematics..." << endl;
	for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
	{
		for (int i_cent = 0; i_cent < n_cent_cuts; i_cent++)
		{
			if ((dataset_type.compare("_pp") !=0 ) && i_cent == 6) continue;
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
								(sys_names[i_sys] == "sys1" && combined_sys_names[i_comb_sys] == "JER") ||
								(sys_names[i_sys] == "sys3" && combined_sys_names[i_comb_sys] == "Unfolding") ||
								(sys_names[i_sys] == "sys4" && combined_sys_names[i_comb_sys] == "MCNonClosure")
								)
							{
								//if running in combined mode, need to redo all uncorrelated systematics present in both pbpb and pp: unfolding, mc-nonclosure, done below

//								cout << Form("%i_%i_%i_%i %i_%i ", i_trk, i_cent, i_jet, i_dR, i_sys, i_comb_sys);
//								cout << Form("%s --> %s", sys_names[i_sys].c_str(), combined_sys_names[i_comb_sys].c_str()) << endl;

								tmp = h_sys_p[i_sys][i_trk][i_cent][i_jet]->GetBinContent(i_dR);
								h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->SetBinContent(i_dR, tmp);
								h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->SetBinError(i_dR, 0.0000001);

								tmp = h_sys_n[i_sys][i_trk][i_cent][i_jet]->GetBinContent(i_dR);
								h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->SetBinContent(i_dR, tmp);
								h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->SetBinError(i_dR, 0.0000001);
							}

							//each systematic needs to be combined. Most of the ones listed below are correlated. The ones that are uncorrelated (UE, CENTJES {mcnonclosure and unfolding are individual anyway}) are only different in the pbpb, so err_pp_sysi = 0, and it can be treated as uncorrelated
							if (

								(
								 (sys_names[i_sys] == "sys6" ||
								  sys_names[i_sys] == "sys7" ||
								  sys_names[i_sys] == "sys8" ||
								  sys_names[i_sys] == "sys9" ||
								  sys_names[i_sys] == "sys16" ||
								  sys_names[i_sys] == "sys17" ||
								  sys_names[i_sys] == "sys18" ||
								  sys_names[i_sys] == "sys19" ||
								  sys_names[i_sys] == "sys20" ||
								  sys_names[i_sys] == "sys21" ||
								  sys_names[i_sys] == "sys22" ||
								  sys_names[i_sys] == "sys23" ) &&
								 combined_sys_names[i_comb_sys] == "JES"
								 )

								||

								(
								 (sys_names[i_sys] == "sys2" ||
								  sys_names[i_sys] == "sys5" ||
								  sys_names[i_sys] == "sys10" ||
								  sys_names[i_sys] == "sys11" ||
								  sys_names[i_sys] == "sys14" ||
								  sys_names[i_sys] == "sys15" ||
								  sys_names[i_sys] == "sys45" ) &&
								 combined_sys_names[i_comb_sys] == "Tracking"
								 )

								||

								(
								 (sys_names[i_sys] == "sys42" ||
								  sys_names[i_sys] == "sys44" ) &&
								 combined_sys_names[i_comb_sys] == "UE"
								 )
								)
							{
//								cout << Form("%i_%i_%i_%i %i_%i ", i_trk, i_cent, i_jet, i_dR, i_sys, i_comb_sys);
//								cout << Form("%s --> %s", sys_names[i_sys].c_str(), combined_sys_names[i_comb_sys].c_str()) << endl;

								tmp = pow(h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->GetBinContent(i_dR),2);
								tmp = tmp + pow(h_sys_p[i_sys][i_trk][i_cent][i_jet]->GetBinContent(i_dR),2);
								h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->SetBinContent(i_dR, sqrt(tmp) );
								h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->SetBinError(i_dR, 0.0000001);

								tmp = pow(h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->GetBinContent(i_dR),2);
								tmp = tmp + pow(h_sys_n[i_sys][i_trk][i_cent][i_jet]->GetBinContent(i_dR),2);
								h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->SetBinContent(i_dR, -sqrt(tmp) );
								h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->SetBinError(i_dR, 0.0000001);

							}

						}
					}

					//uncorrelated uncertainties do separately so has to be done separately; MC nonclosure and unfolding are by themselves anyway

					//uncorrelated so has to be done separately
					if ( (combined_sys_names[i_comb_sys] == "MCNonClosure" ||
						  combined_sys_names[i_comb_sys] == "Unfolding" )
						&&
						(dataset_type.compare("") == 0) )


					{
						//get nominal values and relative uncertainties from pbpb and pp sys files
						name = Form("h_ChPS_final_indR_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
						TH1* h_pbpb_nom = (TH1*)f_pbpb->Get(name.c_str());
						h_pbpb_nom->SetName(Form("%s_pbpb", name.c_str()));

						name = Form("h_ChPS_final_indR_trk%i_cent%i_jetpt%i", i_trk, 6, i_jet);
						TH1* h_pp_nom = (TH1*)f_pp->Get(name.c_str());
						h_pp_nom->SetName(Form("%s_pp", name.c_str()));

						name = Form("h_ChPS_sys_trk%i_cent%i_jetpt%i_%s_p", i_trk, i_cent, i_jet,combined_sys_names[i_comb_sys].c_str());
						TH1* h_pbpb_p = (TH1*)f_pbpb->Get(name.c_str());
						h_pbpb_p->SetName(Form("%s_pbpb", name.c_str()));

						name = Form("h_ChPS_sys_trk%i_cent%i_jetpt%i_%s_p", i_trk, 6, i_jet,combined_sys_names[i_comb_sys].c_str());
						TH1* h_pp_p = (TH1*)f_pp->Get(name.c_str());
						h_pp_p->SetName(Form("%s_pp", name.c_str()));

						name = Form("h_ChPS_sys_trk%i_cent%i_jetpt%i_%s_n", i_trk, i_cent, i_jet,combined_sys_names[i_comb_sys].c_str());
						TH1* h_pbpb_n = (TH1*)f_pbpb->Get(name.c_str());
						h_pbpb_n->SetName(Form("%s_pbpb", name.c_str()));

						name = Form("h_ChPS_sys_trk%i_cent%i_jetpt%i_%s_n", i_trk, 6, i_jet,combined_sys_names[i_comb_sys].c_str());
						TH1* h_pp_n = (TH1*)f_pp->Get(name.c_str());
						h_pp_n->SetName(Form("%s_pp", name.c_str()));

						if (mode == "RDpT")
						{
							//add relative uncertainties in quadrature to get rel uncert on (grouped systematic) i
							//deltaR/R = sqrt( (deltaA/A)^2 + (deltaB/B)^2 )
							double dA_overA, dB_overB, dR_overR;
							for (int i_dR = 1; i_dR <= N_dR; i_dR++)
							{
								dA_overA = h_pbpb_p->GetBinContent(i_dR);
								dB_overB = h_pp_p->GetBinContent(i_dR);
								dR_overR = sqrt( pow(dA_overA,2) + pow(dB_overB,2) );
								h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->SetBinContent(i_dR, dR_overR);
								h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->SetBinError(i_dR, 0.000001);

								dA_overA = h_pbpb_n->GetBinContent(i_dR);
								dB_overB = h_pp_n->GetBinContent(i_dR);
								dR_overR = sqrt( pow(dA_overA,2) + pow(dB_overB,2) );
								h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->SetBinContent(i_dR, -dR_overR);
								h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->SetBinError(i_dR, 0.000001);
							}
						}
						else if (mode == "DeltaDpT")
						{
							//add absolute uncertainties in quadrature to get absolute uncert on systematic i
							//deltaR = sqrt(deltaA^2 + deltaB^2). make sure to save relative uncert
							double dA, dB, dR, dR_overR;
							for (int i_dR = 1; i_dR <= N_dR; i_dR++)
							{
								dA = h_pbpb_p->GetBinContent(i_dR) * h_pbpb_nom->GetBinContent(i_dR);
								dB = h_pp_p->GetBinContent(i_dR) * h_pp_nom->GetBinContent(i_dR);
								dR = sqrt( pow(dA,2) + pow(dB,2) );
								dR_overR = fabs(dR/(h_pbpb_nom->GetBinContent(i_dR) - h_pp_nom->GetBinContent(i_dR)));
								h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->SetBinContent(i_dR, dR_overR);
								h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->SetBinError(i_dR, 0.000001);

								dA = h_pbpb_n->GetBinContent(i_dR) * h_pbpb_nom->GetBinContent(i_dR);
								dB = h_pp_n->GetBinContent(i_dR) * h_pp_nom->GetBinContent(i_dR);
								dR = sqrt( pow(dA,2) + pow(dB,2) );
								dR_overR = fabs(dR/(h_pbpb_nom->GetBinContent(i_dR) - h_pp_nom->GetBinContent(i_dR)));
								h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->SetBinContent(i_dR, -dR_overR);
								h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->SetBinError(i_dR, 0.000001);
							}


						}

						delete h_pbpb_nom;
						delete h_pp_nom;
						delete h_pbpb_p;
						delete h_pp_p;
						delete h_pbpb_n;
						delete h_pp_n;

					}

					//remove fluctuations in JES uncert
//					if (combined_sys_names[i_comb_sys] == "JES")
//					{
//						h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->Smooth(1,"");
//						h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->Smooth(1,"");
//					}

					name = Form("h_%s_sys_trk%i_cent%i_jetpt%i_%s_p",mode.c_str(), i_trk, i_cent, i_jet,combined_sys_names[i_comb_sys].c_str());
					h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->SetTitle(name.c_str());
					h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->Write(name.c_str());

					name = Form("h_%s_sys_trk%i_cent%i_jetpt%i_%s_n",mode.c_str(), i_trk, i_cent, i_jet,combined_sys_names[i_comb_sys].c_str());
					h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->SetTitle(name.c_str());
					h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->Write(name.c_str());

				}
			}
		}
	}

	cout << "Done grouping... getting total..." << endl;
	for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
	{
		for (int i_cent = 0; i_cent < n_cent_cuts; i_cent++)
		{
			if ((dataset_type.compare("_pp") !=0 ) && i_cent == 6) continue;
			else if (dataset_type == "_pp" && i_cent < 6) continue;

			for (int i_trk = 0; i_trk < N_trkpt; i_trk++)
			{
				for (int i_dR = 1; i_dR <= N_dR; i_dR++)
				{
					double tmp_p = 0, tmp_n = 0;
					for (int i_comb_sys = 0; i_comb_sys < combined_sys_names.size(); i_comb_sys++)
					{
						tmp_p = tmp_p+pow(h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->GetBinContent(i_dR),2);
						tmp_n = tmp_n+pow(h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->GetBinContent(i_dR),2);
					}
					h_total_sys_p[i_trk][i_cent][i_jet]->SetBinContent(i_dR, sqrt(tmp_p) );
					h_total_sys_p[i_trk][i_cent][i_jet]->SetBinError(i_dR, 0.0000001 );

					h_total_sys_n[i_trk][i_cent][i_jet]->SetBinContent(i_dR, -sqrt(tmp_n) );
					h_total_sys_n[i_trk][i_cent][i_jet]->SetBinError(i_dR, 0.0000001 );
				}

				name = Form("h_%s_sys_trk%i_cent%i_jetpt%i_total_p",mode.c_str(), i_trk, i_cent, i_jet);
				h_total_sys_p[i_trk][i_cent][i_jet]->Write(name.c_str());

				name = Form("h_%s_sys_trk%i_cent%i_jetpt%i_total_n",mode.c_str(), i_trk, i_cent, i_jet);
				h_total_sys_n[i_trk][i_cent][i_jet]->Write(name.c_str());


			}
		}
	}


	//Cast in terms of track pT if in RDpT mode
	if (mode == "RDpT")
	{
		vector<vector<vector<TH1*>>> h_total_sys_p_inTrk (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
		vector<vector<vector<TH1*>>> h_total_sys_n_inTrk (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));

		for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
		{
			for (int i_cent = 0; i_cent < n_cent_cuts; i_cent++)
			{
				if ((dataset_type == "_PbPb" || mode == "RDpT" || mode == "DeltaDpT") && i_cent == 6) continue;
				else if (dataset_type == "_pp" && i_cent < 6) continue;

				for (int i_dR = 0; i_dR < N_dR; i_dR++)
				{
					name = Form("h_%s_sys_dR%i_cent%i_jetpt%i_total_p",mode.c_str(), i_dR, i_cent, i_jet);
					h_total_sys_p_inTrk[i_dR][i_cent][i_jet] = (TH1*)nom_file->Get(Form("h_%s_final_dR0_cent%i_jetpt8",mode.c_str(), i_cent))->Clone(name.c_str());
					h_total_sys_p_inTrk[i_dR][i_cent][i_jet]->Reset();

					name = Form("h_%s_sys_dR%i_cent%i_jetpt%i_total_n",mode.c_str(), i_dR, i_cent, i_jet);
					h_total_sys_n_inTrk[i_dR][i_cent][i_jet] = (TH1*)nom_file->Get(Form("h_%s_final_dR0_cent%i_jetpt8",mode.c_str(), i_cent))->Clone(name.c_str());
					h_total_sys_n_inTrk[i_dR][i_cent][i_jet]->Reset();

					double val = 0, val_err = 0;
					for (int i_trk = 0; i_trk < N_trkpt; i_trk++)
					{
						val = h_total_sys_p[i_trk][i_cent][i_jet]->GetBinContent(i_dR+1);
						val_err = h_total_sys_p[i_trk][i_cent][i_jet]->GetBinError(i_dR+1);
						h_total_sys_p_inTrk[i_dR][i_cent][i_jet]->SetBinContent(i_trk+1, val);
						h_total_sys_p_inTrk[i_dR][i_cent][i_jet]->SetBinError(i_trk+1, val_err);

						val = h_total_sys_n[i_trk][i_cent][i_jet]->GetBinContent(i_dR+1);
						val_err = h_total_sys_n[i_trk][i_cent][i_jet]->GetBinError(i_dR+1);
						h_total_sys_n_inTrk[i_dR][i_cent][i_jet]->SetBinContent(i_trk+1, val);
						h_total_sys_n_inTrk[i_dR][i_cent][i_jet]->SetBinError(i_trk+1, val_err);
					}

					name = Form("h_%s_sys_dR%i_cent%i_jetpt%i_total_p",mode.c_str(), i_dR, i_cent, i_jet);
					h_total_sys_p_inTrk[i_dR][i_cent][i_jet]->Write(name.c_str());

					name = Form("h_%s_sys_dR%i_cent%i_jetpt%i_total_n",mode.c_str(), i_dR, i_cent, i_jet);
					h_total_sys_n_inTrk[i_dR][i_cent][i_jet]->Write(name.c_str());

				}
			}
		}
	}



	cout << "Drawing... " << endl;
	//drawing

	string y_label = "";
	if (mode == "RDpT") y_label = "Rel. Unc. #it{R}_{ #it{D} (#it{p}_{T}, #it{r})}";
	if (mode == "DeltaDpT") y_label = "Rel. Unc. #delta #it{#Delta}_{ #it{D} (#it{p}_{T}, #it{r})}";
	if (mode == "ChPS") y_label = "Rel. Unc. #it{D} (#it{p}_{T}, #it{r}) [GeV^{-1}]";
	string r_label = "#it{r}";

	TCanvas *c_sys = new TCanvas("c_sys","c_sys", 1200, 600);
	if (dataset_type == "_pp") c_sys->SetCanvasSize(800,600);
	TLegend *legend_sys = new TLegend(0.19, 0.18, 0.70, 0.45, "","brNDC");
	legend_sys->SetTextFont(43);
	legend_sys->SetBorderSize(0);
	legend_sys->SetTextSize(14);
	legend_sys->SetNColumns(2);

	TLatex *ltx = new TLatex();
	ltx->SetTextFont(43);
	ltx->SetTextSize(12);
	TLine *line = new TLine();

	for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
	{
		string jet_label = Form("%1.0f < #it{p}_{T}^{jet} < %1.0f GeV", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));

		for (int i_trk = trk_pt_start; i_trk < trk_pt_end; i_trk++)
		{
			string trk_label = Form("%1.1f < #it{p}_{T} < %1.1f GeV", trkpT_binning->GetBinLowEdge(i_trk+1), trkpT_binning->GetBinUpEdge(i_trk+1));

			legend_sys->Clear();
			c_sys->cd();
			c_sys->Clear();
			c_sys->Divide(3,2);

			bool cent_first_pass = true;
			for (int i_cent = 0; i_cent < n_cent_cuts; i_cent++)
			{
				if ((dataset_type.compare("_pp") !=0 ) && i_cent == 6) continue;
				else if (dataset_type == "_pp" && i_cent < 6) continue;

				string centrality = num_to_cent(31,i_cent);

				SetHStyle_smallify(h_total_sys_p[i_trk][i_cent][i_jet],0, doSmall);
				SetHStyle_smallify(h_total_sys_n[i_trk][i_cent][i_jet],0, doSmall);
				h_total_sys_p[i_trk][i_cent][i_jet]->SetLineWidth(1);
				h_total_sys_n[i_trk][i_cent][i_jet]->SetLineWidth(1);
				h_total_sys_p[i_trk][i_cent][i_jet]->SetMarkerStyle(20);
				h_total_sys_n[i_trk][i_cent][i_jet]->SetMarkerStyle(20);

				h_total_sys_p[i_trk][i_cent][i_jet]->GetYaxis()->SetRangeUser(-0.8,0.8);
				h_total_sys_n[i_trk][i_cent][i_jet]->GetYaxis()->SetRangeUser(-0.8,0.8);
				h_total_sys_p[i_trk][i_cent][i_jet]->GetXaxis()->SetRangeUser(0, r_max_range);
				h_total_sys_n[i_trk][i_cent][i_jet]->GetXaxis()->SetRangeUser(0, r_max_range);
				h_total_sys_p[i_trk][i_cent][i_jet]->GetYaxis()->SetTitle(y_label.c_str());
				h_total_sys_n[i_trk][i_cent][i_jet]->GetYaxis()->SetTitle(y_label.c_str());
				h_total_sys_p[i_trk][i_cent][i_jet]->GetXaxis()->SetTitle(r_label.c_str());
				h_total_sys_n[i_trk][i_cent][i_jet]->GetXaxis()->SetTitle(r_label.c_str());


				c_sys->cd(i_cent+1);
//				h_total_sys_p[i_trk][i_cent][i_jet]->GetYaxis()->SetRangeUser(-0.2, 0.2);
				h_total_sys_p[i_trk][i_cent][i_jet]->Draw("p");
				h_total_sys_n[i_trk][i_cent][i_jet]->Draw("p same");

//				for (int i_sys = 0; i_sys < sys_names.size(); i_sys++)
//				{
//					SetHStyle_smallify(h_sys_p[i_sys][i_trk][i_cent][i_jet],i_sys+1,1);
//					SetHStyle_smallify(h_sys_n[i_sys][i_trk][i_cent][i_jet],i_sys+1,1);
//					h_sys_p[i_sys][i_trk][i_cent][i_jet]->SetLineStyle(1);
//					h_sys_n[i_sys][i_trk][i_cent][i_jet]->SetLineStyle(1);
//					h_sys_p[i_sys][i_trk][i_cent][i_jet]->SetLineWidth(1);
//					h_sys_n[i_sys][i_trk][i_cent][i_jet]->SetLineWidth(1);
//
//					if (cent_first_pass) legend_sys->AddEntry(h_sys_p[i_sys][i_trk][i_cent][i_jet],sys_i_name[sys_names[i_sys]].c_str(),"lp");
//					h_sys_p[i_sys][i_trk][i_cent][i_jet]->Draw("p same");
//					h_sys_n[i_sys][i_trk][i_cent][i_jet]->Draw("p same");
//
//				}
				for (int i_comb_sys = 0; i_comb_sys < combined_sys_names.size(); i_comb_sys++)
				{
					SetHStyle_smallify(h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet],i_comb_sys+1, doSmall);
					SetHStyle_smallify(h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet],i_comb_sys+1, doSmall);
					h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->SetLineWidth(1);
					h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->SetLineWidth(1);
					h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->SetLineStyle(2);
					h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->SetLineStyle(2);

					bool empty_hist_p = false;
					bool empty_hist_n = false;
					if (h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->GetMaximum() == 0 && h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->GetMinimum() == 0) empty_hist_p = true;
					if (h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->GetMaximum() == 0 && h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->GetMinimum() == 0) empty_hist_n = true;

					if (cent_first_pass) legend_sys->AddEntry(h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet],combined_sys_names[i_comb_sys].c_str(),"lp");

					//DO NOT SET RANGES ON COMB PLOTS. THIS MESSES UP THE DRAW SAME OPTION
					h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->GetXaxis()->SetRangeUser(0, r_max_range);
					h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->GetXaxis()->SetRangeUser(0, r_max_range);
					h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->GetYaxis()->SetTitle(y_label.c_str());
					h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->GetYaxis()->SetTitle(y_label.c_str());
					h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->GetXaxis()->SetTitle(r_label.c_str());
					h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->GetXaxis()->SetTitle(r_label.c_str());

					c_sys->cd(i_cent+1);
					if (!empty_hist_p) h_comb_sys_p[i_comb_sys][i_trk][i_cent][i_jet]->Draw("p same");
					if (!empty_hist_n) h_comb_sys_n[i_comb_sys][i_trk][i_cent][i_jet]->Draw("p same");
				}



				c_sys->cd(i_cent+1);
				if (dataset_type == "_pp") c_sys->cd();
				if (cent_first_pass) legend_sys->AddEntry(h_total_sys_p[i_trk][i_cent][i_jet],"Total","pl");

				legend_sys->Draw();
				ltx->SetTextAlign(31);
				ltx->SetTextSize(12);
				if (dataset_type == "_pp") ltx->SetTextSize(24);
				if (dataset_type == "_pp") legend_sys->SetTextSize(24);

				double x_left = 0.19, x_right = 0.93, y = 0.88, y_diff = 0.045;
				ltx->DrawLatexNDC(x_right, y, Form("%s", trk_label.c_str()));
				ltx->DrawLatexNDC(x_right, y=y-y_diff, Form("%s", jet_label.c_str()));
				ltx->DrawLatexNDC(x_right, y=y-y_diff, Form("%s", centrality.c_str()));

				cent_first_pass = false;


			}

			if (dataset_type == "_pp") c_sys->cd();
			else c_sys->cd(1);

			ltx->SetTextAlign(11);
			double x_left = 0.19, x_right = 0.93, y = 0.88, y_diff = 0.045;
			ltx->DrawLatexNDC(x_left, y, "#scale[1.5]{#font[72]{ATLAS} Internal}");
			if (mode == "RDpT" || mode == "DeltaDpT")
			{
				ltx->DrawLatexNDC(x_left, y=y-y_diff, "Pb+Pb #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}");
				ltx->DrawLatexNDC(x_left, y=y-y_diff, "#it{pp} #sqrt{#font[12]{s}} = 5.02 TeV, 25 pb^{-1}");
				ltx->DrawLatexNDC(x_left, y=y-y_diff, Form("anti-#font[12]{k}_{#font[12]{t}} R=0.4"));
			}
			else
			{
				if (dataset_type == "_PbPb") ltx->DrawLatexNDC(x_left, y=y-y_diff, "Pb+Pb #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}");
				if (dataset_type == "_pp") ltx->DrawLatexNDC(x_left, y=y-y_diff, "#it{pp} #sqrt{#font[12]{s}} = 5.02 TeV, 25 pb^{-1}");
				ltx->DrawLatexNDC(x_left, y=y-y_diff, Form("anti-#font[12]{k}_{#font[12]{t}} R=0.4"));
			}

			string pdf_label = "";
			if (i_trk == trk_pt_start && i_jet == jet_pt_start) pdf_label = "(";
			if (i_trk == trk_pt_end-1 && i_jet == jet_pt_end-1) pdf_label = ")";
			c_sys->Print(Form("output_pdf_nominal/systematics/Summary_%s_dR_sys%s_error.pdf%s",mode.c_str(), dataset_type.c_str(), pdf_label.c_str()), Form("Title:trk%i_jetpt%i", i_trk, i_jet));

			//For compiledCode
//			c_sys->Print(Form("compiledCode/systematics/Summary_%s_dR_sys%s_error.pdf%s",mode.c_str(), dataset_type.c_str(), pdf_label.c_str()), Form("Title:trk%i_jetpt%i", i_trk, i_jet));



		}
	}


	delete c_sys;
	delete legend_sys;
	delete ltx;
	delete line;
	delete m_config;
	delete nom_file;
	delete output_file;

	cout << "######### Done Systematics #########" << endl;

}




