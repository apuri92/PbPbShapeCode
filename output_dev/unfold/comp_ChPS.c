#include "../functions/global_variables.h"
#include "TEnv.h"

void comp_ChPS(string config_file = "ff_config.cfg")
{
	cout << "######### DOING COMP_ChPS #########" << endl;

	SetAtlasStyle();
	gErrorIgnoreLevel = 3001;

	//	##############	Reading config	##############"
	TEnv *m_config = new TEnv();
	m_config->ReadFile(config_file.c_str(), EEnvLevel(1));

	int isMC = 1; isMC = m_config->GetValue("isMC", isMC);
	int sys_mode = -1; sys_mode = m_config->GetValue("sys_mode", sys_mode);

	int centrality_scheme = 31; centrality_scheme = m_config->GetValue("centrality_scheme", centrality_scheme);
	int verbose = 0; verbose = m_config->GetValue("verbose", verbose);
	int draw_mode = 1; draw_mode = m_config->GetValue("draw_mode", draw_mode);

	double r_max_range = 0.8;
	std::string did = "data";
	if (isMC) did = "MC";

	if (verbose) m_config->Print();
	//	##############	Config done	##############"
	std::string sys_path = "";
	if (sys_mode == 0) sys_path = Form("nominal");
	else if (sys_mode > 50 || sys_mode < 0) sys_path = Form("c%i", sys_mode);
	else sys_path = Form("sys%i", sys_mode);

	string rdptr_label = "#it{R}_{#it{D} (p_{#it{T}}, #it{r})}";
	string dptr_label = "#it{D} (p_{#it{T}}, #it{r}) [GeV^{-1}]";
	string delta_dptr_label = "#Delta #it{D} (p_{#it{T}}, #it{r}) [GeV^{-1}]";

	cout << Form("Doing in %s mode", did.c_str()) << endl;

	TFile *f_PbPb = new TFile(Form("output_pdf_%s/root/final_ChPS_%s_PbPb.root", sys_path.c_str(), did.c_str()));
	TFile *f_pp = new TFile(Form("output_pdf_%s/root/final_ChPS_%s_pp.root", sys_path.c_str(), did.c_str()));
	TFile *f_RDpT = new TFile(Form("output_pdf_%s/root/final_RDpT_%s.root", sys_path.c_str(), did.c_str()), "recreate");

	cout << "Using files:" << endl;
	cout << f_PbPb->GetName() << endl;
	cout << f_pp->GetName() << endl;



	TFile *f_FF_PbPb = new TFile("FF_files/Uncertainties_eta_4_dpt_PbPb.root");
	TFile *f_FF_pp = new TFile("FF_files/Uncertainties_eta_4_dpt_pp.root");
	TFile *f_FF_ratio = new TFile("FF_files/Uncertainties_eta_4_dpt_ratio.root");

	vector<TFile*> FF_PbPb;
	for (int i = 0; i < 6; i++)
	{
		string name = Form("FF_files/dpt_PbPb_data_removeElectrons_significancecuts_FS_itr_4_itr1d_4_eta_4_cent%i.root", i);
		FF_PbPb.push_back(new TFile( name.c_str() ) );
	}

	TAxis* dR_binning = (TAxis*)f_PbPb->Get("dR_binning");
	TAxis* jetpT_binning = (TAxis*)f_PbPb->Get("jetpT_binning");
	TAxis* trkpT_binning = (TAxis*)f_PbPb->Get("trkpT_binning");

	int N_dR = dR_binning->GetNbins();
	int N_jetpt = jetpT_binning->GetNbins();
	int N_trkpt = trkpT_binning->GetNbins();


	f_RDpT->cd();
	dR_binning->Write("dR_binning");
	jetpT_binning->Write("jetpT_binning");
	trkpT_binning->Write("trkpT_binning");

	double array_dr_bins[N_dR+1];
	for (int i_dR = 0; i_dR <= N_dR; i_dR++) array_dr_bins[i_dR] = dR_binning->GetBinLowEdge(i_dR+1);
	//indR
	vector<vector<vector<TH1*>>> h_ChPS_raw_PbPb_indR (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_raw_pp_indR (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_raw_ratio_indR (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));


	vector<vector<vector<TH1*>>> h_ChPS_raw_subtr_PbPb_indR (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_raw_subtr_pp_indR (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_raw_subtr_ratio_indR (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));

	vector<vector<vector<TH1*>>> h_ChPS_raw_subtr_unf_PbPb_indR (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_raw_subtr_unf_pp_indR (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_raw_subtr_unf_ratio_indR (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));

	vector<vector<vector<TH1*>>> h_ChPS_final_PbPb_indR (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_final_pp_indR (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_final_ratio_indR (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_final_diff_indR (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));

	//injet
	vector<vector<TH1*>> h_ChPS_final_PbPb_injet (n_cent_cuts, vector<TH1*> (N_jetpt));

	vector<vector<TH1*>> h_ChPS_final_pp_injet (n_cent_cuts, vector<TH1*> (N_jetpt));
	vector<vector<TH1*>> h_ChPS_final_ratio_injet (n_cent_cuts, vector<TH1*> (N_jetpt));

	vector<vector<TGraphErrors*>> g_ChPS_final_PbPb_injet (n_cent_cuts, vector<TGraphErrors*> (N_jetpt));
	vector<vector<TGraphErrors*>> g_ChPS_final_pp_injet (n_cent_cuts, vector<TGraphErrors*> (N_jetpt));
	vector<vector<TGraphErrors*>> g_ChPS_final_ratio_injet (n_cent_cuts, vector<TGraphErrors*> (N_jetpt));


	vector<vector<TH1*>> h_FF_final_PbPb_injet (n_cent_cuts, vector<TH1*> (N_jetpt));
	vector<vector<TH1*>> h_FF_final_pp_injet (n_cent_cuts, vector<TH1*> (N_jetpt));
	vector<vector<TH1*>> h_FF_final_ratio_injet (n_cent_cuts, vector<TH1*> (N_jetpt));

	vector<vector<TGraphErrors*>> g_FF_final_PbPb_injet (n_cent_cuts, vector<TGraphErrors*> (N_jetpt));
	vector<vector<TGraphErrors*>> g_FF_final_pp_injet (n_cent_cuts, vector<TGraphErrors*> (N_jetpt));
	vector<vector<TGraphErrors*>> g_FF_final_ratio_injet (n_cent_cuts, vector<TGraphErrors*> (N_jetpt));

	vector<vector<TGraphErrors*>> fit_quality_ChPS_PbPb_injet (n_cent_cuts, vector<TGraphErrors*> (N_jetpt));
	vector<vector<TGraphErrors*>> fit_quality_ChPS_pp_injet (n_cent_cuts, vector<TGraphErrors*> (N_jetpt));
	vector<vector<TGraphErrors*>> fit_quality_FF_PbPb_injet (n_cent_cuts, vector<TGraphErrors*> (N_jetpt));
	vector<vector<TGraphErrors*>> fit_quality_FF_pp_injet (n_cent_cuts, vector<TGraphErrors*> (N_jetpt));

	vector<vector<TF1*>> fit_FF_final_PbPb_injet (n_cent_cuts, vector<TF1*> (N_jetpt));
	vector<vector<TF1*>> fit_FF_final_pp_injet (n_cent_cuts, vector<TF1*> (N_jetpt));

	vector<vector<TGraphAsymmErrors*>> g_FF_final_PbPb_sys_injet (n_cent_cuts, vector<TGraphAsymmErrors*> (N_jetpt));
	vector<vector<TGraphAsymmErrors*>> g_FF_final_pp_sys_injet (n_cent_cuts, vector<TGraphAsymmErrors*> (N_jetpt));
	vector<vector<TGraphAsymmErrors*>> g_FF_final_ratio_sys_injet (n_cent_cuts, vector<TGraphAsymmErrors*> (N_jetpt));


	//in track pT
	vector<vector<vector<TH1*>>> h_ChPS_raw_subtr_PbPb (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_raw_subtr_pp (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_raw_subtr_ratio (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));

	vector<vector<vector<TH1*>>> h_ChPS_raw_subtr_unf_PbPb (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_raw_subtr_unf_pp (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_raw_subtr_unf_ratio (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));

	vector<vector<vector<TH1*>>> h_ChPS_final_PbPb (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_final_pp (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_final_ratio (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));

	vector<vector<vector<TH1*>>> h_ChPS_truth_PbPb (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_truth_pp (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_truth_ratio (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));



	string name;
	string pdf_label;

	TLine *line = new TLine();
	line->SetLineColor(kBlack);
	line->SetLineStyle(3);

	TLatex *ltx = new TLatex();
	ltx->SetTextFont(43);
	ltx->SetTextSize(12);

	double trk_pt_lo = 1.;
	double trk_pt_hi = 63;
	double ratio_lo = 0.;
	double ratio_hi = 5;

	int jet_pt_start = 7;
	int jet_pt_end = 11;

	double area_injet = TMath::Pi() * (0.4*0.4);

	for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
	{
		for (int i_cent = 0; i_cent < 6; i_cent++)
		{

			for (int i_trk = 0; i_trk < N_trkpt; i_trk++)
			{
				if (i_cent == 0) //get the pp plots from cent6
				{
					name = Form("h_ChPS_raw_indR_trk%i_cent%i_jetpt%i", i_trk, 6, i_jet);
					h_ChPS_raw_pp_indR.at(i_trk).at(6).at(i_jet) = (TH1*)f_pp->Get(name.c_str())->Clone(Form("pp_%s",name.c_str()));

					name = Form("h_ChPS_raw_subtr_indR_trk%i_cent%i_jetpt%i", i_trk, 6, i_jet);
					h_ChPS_raw_subtr_pp_indR.at(i_trk).at(6).at(i_jet) = (TH1*)f_pp->Get(name.c_str())->Clone(Form("pp_%s",name.c_str()));

					name = Form("h_ChPS_raw_subtr_unf_indR_trk%i_cent%i_jetpt%i", i_trk, 6, i_jet);
					h_ChPS_raw_subtr_unf_pp_indR.at(i_trk).at(6).at(i_jet) = (TH1*)f_pp->Get(name.c_str())->Clone(Form("pp_%s",name.c_str()));

					name = Form("h_ChPS_final_indR_trk%i_cent%i_jetpt%i", i_trk, 6, i_jet);
					h_ChPS_final_pp_indR.at(i_trk).at(6).at(i_jet) = (TH1*)f_pp->Get(name.c_str())->Clone(Form("pp_%s",name.c_str()));


				}

				name = Form("h_ChPS_raw_indR_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_ChPS_raw_PbPb_indR.at(i_trk).at(i_cent).at(i_jet) = (TH1*)f_PbPb->Get(name.c_str())->Clone(Form("PbPb_%s",name.c_str()));

				name = Form("h_ChPS_raw_subtr_indR_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_ChPS_raw_subtr_PbPb_indR.at(i_trk).at(i_cent).at(i_jet) = (TH1*)f_PbPb->Get(name.c_str())->Clone(Form("PbPb_%s",name.c_str()));

				name = Form("h_ChPS_raw_subtr_unf_indR_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_ChPS_raw_subtr_unf_PbPb_indR.at(i_trk).at(i_cent).at(i_jet) = (TH1*)f_PbPb->Get(name.c_str())->Clone(Form("PbPb_%s",name.c_str()));

				name = Form("h_ChPS_final_indR_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_ChPS_final_PbPb_indR.at(i_trk).at(i_cent).at(i_jet) = (TH1*)f_PbPb->Get(name.c_str())->Clone(Form("PbPb_%s",name.c_str()));


				//Ratios
				name = Form("h_ChPS_raw_indR_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_ChPS_raw_ratio_indR.at(i_trk).at(i_cent).at(i_jet) = (TH1*)h_ChPS_raw_PbPb_indR.at(i_trk).at(i_cent).at(i_jet)->Clone(name.c_str());
				h_ChPS_raw_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->Divide(h_ChPS_raw_pp_indR.at(i_trk).at(6).at(i_jet));
				h_ChPS_raw_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->GetYaxis()->SetTitle("Raw #frac{PbPb}{pp} _{(#frac{1}{Area N_{Jets}} #frac{dN}{dp_{T}})}");

				name = Form("h_ChPS_raw_subtr_indR_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_ChPS_raw_subtr_ratio_indR.at(i_trk).at(i_cent).at(i_jet) = (TH1*)h_ChPS_raw_subtr_PbPb_indR.at(i_trk).at(i_cent).at(i_jet)->Clone(name.c_str());
				h_ChPS_raw_subtr_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->Divide(h_ChPS_raw_subtr_pp_indR.at(i_trk).at(6).at(i_jet));
				h_ChPS_raw_subtr_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->GetYaxis()->SetTitle("Raw+Subtr #frac{PbPb}{pp} _{(#frac{1}{Area N_{Jets}} #frac{dN}{dp_{T}})}");

				name = Form("h_ChPS_raw_subtr_unf_indR_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_ChPS_raw_subtr_unf_ratio_indR.at(i_trk).at(i_cent).at(i_jet) = (TH1*)h_ChPS_raw_subtr_unf_PbPb_indR.at(i_trk).at(i_cent).at(i_jet)->Clone(name.c_str());
				h_ChPS_raw_subtr_unf_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->Divide(h_ChPS_raw_subtr_unf_pp_indR.at(i_trk).at(6).at(i_jet));
				h_ChPS_raw_subtr_unf_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->GetYaxis()->SetTitle("Raw+Subtr+Unf #frac{PbPb}{pp} _{(#frac{1}{Area N_{Jets}} #frac{dN}{dp_{T}})}");

				name = Form("h_ChPS_final_ratio_PbPb_pp_indR_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_ChPS_final_ratio_indR.at(i_trk).at(i_cent).at(i_jet) = (TH1*)h_ChPS_final_PbPb_indR.at(i_trk).at(i_cent).at(i_jet)->Clone(name.c_str());
				h_ChPS_final_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->Divide(h_ChPS_final_pp_indR.at(i_trk).at(6).at(i_jet));
				h_ChPS_final_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->GetYaxis()->SetTitle(rdptr_label.c_str());

				name = Form("h_ChPS_final_diff_PbPb_pp_indR_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_ChPS_final_diff_indR.at(i_trk).at(i_cent).at(i_jet) = (TH1*)h_ChPS_final_PbPb_indR.at(i_trk).at(i_cent).at(i_jet)->Clone(name.c_str());
				h_ChPS_final_diff_indR.at(i_trk).at(i_cent).at(i_jet)->Add(h_ChPS_final_pp_indR.at(i_trk).at(6).at(i_jet),-1);
				h_ChPS_final_diff_indR.at(i_trk).at(i_cent).at(i_jet)->GetYaxis()->SetTitle(delta_dptr_label.c_str());

				f_RDpT->cd();
				name = Form("h_RDpT_final_indR_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_ChPS_final_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->Write(name.c_str());

				name = Form("h_DeltaDpT_final_indR_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_ChPS_final_diff_indR.at(i_trk).at(i_cent).at(i_jet)->Write(name.c_str());

			}

			//injet
			if (i_cent == 0)
			{
				name = Form("h_dpt_pp_%i", i_jet);
				h_FF_final_pp_injet.at(6).at(i_jet) = (TH1*)f_FF_pp->Get(name.c_str())->Clone(Form("pp_%s",name.c_str()));
				g_FF_final_pp_injet.at(6).at(i_jet) = new TGraphErrors(h_FF_final_pp_injet.at(6).at(i_jet));
				name = Form("g_dpt_pp_uncert_%i", i_jet);
				g_FF_final_pp_sys_injet.at(6).at(i_jet) = (TGraphAsymmErrors*)f_FF_pp->Get(name.c_str())->Clone(Form("pp_%s",name.c_str()));

				name = Form("h_ChPS_final_injet_cent%i_jetpt%i", 6, i_jet);
				h_ChPS_final_pp_injet.at(6).at(i_jet) = (TH1*)f_pp->Get(name.c_str())->Clone(Form("pp_%s",name.c_str()));
				h_ChPS_final_pp_injet.at(6).at(i_jet)->Scale(area_injet);
				g_ChPS_final_pp_injet.at(6).at(i_jet) = new TGraphErrors(h_ChPS_final_pp_injet.at(6).at(i_jet));

			}

			name = Form("h_ChPS_final_injet_cent%i_jetpt%i", i_cent, i_jet);
			h_ChPS_final_PbPb_injet.at(i_cent).at(i_jet) = (TH1*)f_PbPb->Get(name.c_str())->Clone(Form("PbPb_%s",name.c_str()));
			h_ChPS_final_PbPb_injet.at(i_cent).at(i_jet)->Scale(area_injet);
			g_ChPS_final_PbPb_injet.at(i_cent).at(i_jet) = new TGraphErrors(h_ChPS_final_PbPb_injet.at(i_cent).at(i_jet));

			name = Form("h_ChPS_final_ratio_PbPb_pp_cent%i_jetpt%i", i_cent, i_jet);
			h_ChPS_final_ratio_injet.at(i_cent).at(i_jet) = (TH1*)h_ChPS_final_PbPb_injet.at(i_cent).at(i_jet)->Clone(name.c_str());
			h_ChPS_final_ratio_injet.at(i_cent).at(i_jet)->Divide(h_ChPS_final_pp_injet.at(6).at(i_jet));
			h_ChPS_final_ratio_injet.at(i_cent).at(i_jet)->GetYaxis()->SetTitle(Form("%s (r < 0.4)", rdptr_label.c_str()));
			g_ChPS_final_ratio_injet.at(i_cent).at(i_jet) = new TGraphErrors(h_ChPS_final_ratio_injet.at(i_cent).at(i_jet));

//			h_FF_final_PbPb_injet.at(i_cent).at(i_jet) = (TH1*)f_FF_PbPb->Get(Form("h_dpt_PbPb_%i_cent%i", i_jet, i_cent));
//			g_FF_final_PbPb_injet.at(i_cent).at(i_jet) = new TGraphErrors(h_FF_final_PbPb_injet.at(i_cent).at(i_jet));
//			g_FF_final_PbPb_sys_injet.at(i_cent).at(i_jet) = (TGraphAsymmErrors*)f_FF_PbPb->Get(Form("g_dpt_PbPb_uncert_%i_cent%i", i_jet, i_cent));

			h_FF_final_PbPb_injet.at(i_cent).at(i_jet) = (TH1*)FF_PbPb[i_cent]->Get(Form("dpt_%i", i_jet));
			name = Form("FF_PbPb_dpt_cent%i_jet%i", i_cent, i_jet);
			h_FF_final_PbPb_injet.at(i_cent).at(i_jet)->SetName(name.c_str());
			g_FF_final_PbPb_injet.at(i_cent).at(i_jet) = new TGraphErrors(h_FF_final_PbPb_injet.at(i_cent).at(i_jet));
			g_FF_final_PbPb_sys_injet.at(i_cent).at(i_jet) = (TGraphAsymmErrors*)g_FF_final_PbPb_injet.at(i_cent).at(i_jet)->Clone(Form("g_Clone_%s", name.c_str()));

			h_FF_final_ratio_injet.at(i_cent).at(i_jet) = (TH1*)f_FF_ratio->Get(Form("h_dpt_ratio_%i_cent%i", i_jet, i_cent));
			g_FF_final_ratio_injet.at(i_cent).at(i_jet) = new TGraphErrors(h_FF_final_ratio_injet.at(i_cent).at(i_jet));
			g_FF_final_ratio_sys_injet.at(i_cent).at(i_jet) = (TGraphAsymmErrors*)f_FF_ratio->Get(Form("g_dpt_ratio_uncert_%i_cent%i", i_jet, i_cent));




			//add as function of jet pt later
			for (int i_dR = 0; i_dR < N_dR; i_dR++)
			{
				if (i_cent == 0) //get the pp plots from cent6
				{
					name = Form("h_ChPS_raw_subtr_dR%i_cent%i_jetpt%i", i_dR, 6, i_jet);
					h_ChPS_raw_subtr_pp.at(i_dR).at(6).at(i_jet) = (TH1*)f_pp->Get(name.c_str())->Clone(Form("pp_%s",name.c_str()));

					name = Form("h_ChPS_raw_subtr_unf_dR%i_cent%i_jetpt%i", i_dR, 6, i_jet);
					h_ChPS_raw_subtr_unf_pp.at(i_dR).at(6).at(i_jet) = (TH1*)f_pp->Get(name.c_str())->Clone(Form("pp_%s",name.c_str()));

					name = Form("h_ChPS_final_dR%i_cent%i_jetpt%i", i_dR, 6, i_jet);
					h_ChPS_final_pp.at(i_dR).at(6).at(i_jet) = (TH1*)f_pp->Get(name.c_str())->Clone(Form("pp_%s",name.c_str()));

					name = Form("h_ChPS_truth_dR%i_cent%i_jetpt%i", i_dR, 6, i_jet);
					h_ChPS_truth_pp.at(i_dR).at(6).at(i_jet) = (TH1*)f_pp->Get(name.c_str())->Clone(Form("pp_%s",name.c_str()));

				}

				name = Form("h_ChPS_raw_subtr_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet);
				h_ChPS_raw_subtr_PbPb.at(i_dR).at(i_cent).at(i_jet) = (TH1*)f_PbPb->Get(name.c_str())->Clone(Form("PbPb_%s",name.c_str()));

				name = Form("h_ChPS_raw_subtr_unf_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet);
				h_ChPS_raw_subtr_unf_PbPb.at(i_dR).at(i_cent).at(i_jet) = (TH1*)f_PbPb->Get(name.c_str())->Clone(Form("PbPb_%s",name.c_str()));

				name = Form("h_ChPS_final_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet);
				h_ChPS_final_PbPb.at(i_dR).at(i_cent).at(i_jet) = (TH1*)f_PbPb->Get(name.c_str())->Clone(Form("PbPb_%s",name.c_str()));

				name = Form("h_ChPS_truth_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet);
				h_ChPS_truth_PbPb.at(i_dR).at(i_cent).at(i_jet) = (TH1*)f_PbPb->Get(name.c_str())->Clone(Form("PbPb_%s",name.c_str()));


				//Ratios
				name = Form("h_ChPS_raw_subtr_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet);
				h_ChPS_raw_subtr_ratio.at(i_dR).at(i_cent).at(i_jet) = (TH1*)h_ChPS_raw_subtr_PbPb.at(i_dR).at(i_cent).at(i_jet)->Clone(name.c_str());
				h_ChPS_raw_subtr_ratio.at(i_dR).at(i_cent).at(i_jet)->Divide(h_ChPS_raw_subtr_pp.at(i_dR).at(6).at(i_jet));
				h_ChPS_raw_subtr_ratio.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetTitle("Raw+Subtr #frac{PbPb}{pp} _{(#frac{1}{Area N_{Jets}} #frac{dN}{dp_{T}})}");
				h_ChPS_raw_subtr_ratio.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetRangeUser(trk_pt_lo, trk_pt_hi);

				name = Form("h_ChPS_raw_subtr_unf_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet);
				h_ChPS_raw_subtr_unf_ratio.at(i_dR).at(i_cent).at(i_jet) = (TH1*)h_ChPS_raw_subtr_unf_PbPb.at(i_dR).at(i_cent).at(i_jet)->Clone(name.c_str());
				h_ChPS_raw_subtr_unf_ratio.at(i_dR).at(i_cent).at(i_jet)->Divide(h_ChPS_raw_subtr_unf_pp.at(i_dR).at(6).at(i_jet));
				h_ChPS_raw_subtr_unf_ratio.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetTitle("Raw+Subtr+Unf #frac{PbPb}{pp} _{(#frac{1}{Area N_{Jets}} #frac{dN}{dp_{T}})}");
				h_ChPS_raw_subtr_unf_ratio.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetRangeUser(trk_pt_lo, trk_pt_hi);

				name = Form("h_ChPS_final_ratio_PbPb_pp_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet);
				h_ChPS_final_ratio.at(i_dR).at(i_cent).at(i_jet) = (TH1*)h_ChPS_final_PbPb.at(i_dR).at(i_cent).at(i_jet)->Clone(name.c_str());
				h_ChPS_final_ratio.at(i_dR).at(i_cent).at(i_jet)->Divide(h_ChPS_final_pp.at(i_dR).at(6).at(i_jet));
				h_ChPS_final_ratio.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetTitle(rdptr_label.c_str());
				h_ChPS_final_ratio.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetRangeUser(trk_pt_lo, trk_pt_hi);

				name = Form("h_ChPS_truth_ratio_PbPb_pp_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet);
				h_ChPS_truth_ratio.at(i_dR).at(i_cent).at(i_jet) = (TH1*)h_ChPS_final_PbPb.at(i_dR).at(i_cent).at(i_jet)->Clone(name.c_str());
				h_ChPS_truth_ratio.at(i_dR).at(i_cent).at(i_jet)->Divide(h_ChPS_final_pp.at(i_dR).at(6).at(i_jet));
				h_ChPS_truth_ratio.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetTitle(Form("Truth %s", rdptr_label.c_str()));
				h_ChPS_truth_ratio.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetRangeUser(trk_pt_lo, trk_pt_hi);

				f_RDpT->cd();
				name = Form("h_RDpT_final_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet);
				h_ChPS_final_ratio.at(i_dR).at(i_cent).at(i_jet)->Write(name.c_str());


			}
		}
	}

	if (draw_mode)
	{
		

		// drawing
		{
			cout << "Doing Final ChPS ratio (PbPb/pp) plots in dR" << endl;
			
			TCanvas *c_ChPS_raw_indR = new TCanvas("c_ChPS_raw_indR","c_ChPS_raw_indR",900,600);
			TLegend *legend_ChPS_raw_indR = new TLegend(0.19, 0.70, 0.40, 0.92, "","brNDC");
			legend_ChPS_raw_indR->SetTextFont(43);
			legend_ChPS_raw_indR->SetBorderSize(0);
			legend_ChPS_raw_indR->SetTextSize(10);

			TCanvas *c_ChPS_raw_subtr_indR = new TCanvas("c_ChPS_raw_subtr_indR","c_ChPS_raw_subtr_indR",900,600);
			TLegend *legend_ChPS_raw_subtr_indR = new TLegend(0.19, 0.70, 0.40, 0.92, "","brNDC");
			legend_ChPS_raw_subtr_indR->SetTextFont(43);
			legend_ChPS_raw_subtr_indR->SetBorderSize(0);
			legend_ChPS_raw_subtr_indR->SetTextSize(10);
			
			TCanvas *c_ChPS_raw_subtr_unf_indR = new TCanvas("c_ChPS_raw_subtr_unf_indR","c_ChPS_raw_subtr_unf_indR",900,600);
			TLegend *legend_ChPS_raw_subtr_unf_indR = new TLegend(0.19, 0.70, 0.40, 0.92, "","brNDC");
			legend_ChPS_raw_subtr_unf_indR->SetTextFont(43);
			legend_ChPS_raw_subtr_unf_indR->SetBorderSize(0);
			legend_ChPS_raw_subtr_unf_indR->SetTextSize(10);
			
			TCanvas *c_ChPS_final_indR = new TCanvas("c_ChPS_final_indR","c_ChPS_final_indR",900,600);
			TLegend *legend_ChPS_final_indR = new TLegend(0.19, 0.60, 0.40, 0.92, "","brNDC");
			legend_ChPS_final_indR->SetTextFont(43);
			legend_ChPS_final_indR->SetBorderSize(0);
			legend_ChPS_final_indR->SetTextSize(12);
			
			TCanvas *c_ChPS_final_diff_indR = new TCanvas("c_ChPS_final_diff_indR","c_ChPS_final_diff_indR",900,600);
			TLegend *legend_ChPS_final_diff_indR = new TLegend(0.19, 0.60, 0.40, 0.92, "","brNDC");
			legend_ChPS_final_diff_indR->SetTextFont(43);
			legend_ChPS_final_diff_indR->SetBorderSize(0);
			legend_ChPS_final_diff_indR->SetTextSize(12);


			int jet_itr = 0;
			for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
			{
				string jet_label = Form("%1.0f < #it{p}_{T}^{jet} < %1.0f", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));
				
				c_ChPS_raw_indR->cd();
				c_ChPS_raw_indR->Clear();
				c_ChPS_raw_indR->Divide(3,2);

				c_ChPS_raw_subtr_indR->cd();
				c_ChPS_raw_subtr_indR->Clear();
				c_ChPS_raw_subtr_indR->Divide(3,2);
				
				c_ChPS_raw_subtr_unf_indR->cd();
				c_ChPS_raw_subtr_unf_indR->Clear();
				c_ChPS_raw_subtr_unf_indR->Divide(3,2);
				
				c_ChPS_final_indR->cd();
				c_ChPS_final_indR->Clear();
				c_ChPS_final_indR->Divide(3,2);

				c_ChPS_final_diff_indR->cd();
				c_ChPS_final_diff_indR->Clear();
				c_ChPS_final_diff_indR->Divide(3,2);

				bool first_pass_cent = true;
				for (int i_cent = 0; i_cent < 6; i_cent++)
				{
					string centrality = num_to_cent(centrality_scheme,i_cent);
					
					
					int trk_itr = 0;
					for (int i_trk = 0; i_trk < N_trkpt; i_trk++)
					{
						//					if (i_trk != 2 &&
						//						i_trk != 5 &&
						//						i_trk != 7 &&
						//						i_trk != 8) continue;
						//1, 5, 10, 40
						
						if (i_trk < 2 || i_trk > 8) continue;
						
						string trk_label = Form("%1.1f < #it{p}_{T}^{trk} < %1.1f", trkpT_binning->GetBinLowEdge(i_trk+1), trkpT_binning->GetBinUpEdge(i_trk+1));
						
						SetHStyle_smallify(h_ChPS_raw_ratio_indR.at(i_trk).at(i_cent).at(i_jet), trk_itr, 1);
						SetHStyle_smallify(h_ChPS_raw_subtr_ratio_indR.at(i_trk).at(i_cent).at(i_jet), trk_itr, 1);
						SetHStyle_smallify(h_ChPS_raw_subtr_unf_ratio_indR.at(i_trk).at(i_cent).at(i_jet), trk_itr, 1);
						SetHStyle_smallify(h_ChPS_final_ratio_indR.at(i_trk).at(i_cent).at(i_jet), trk_itr, 1);
						SetHStyle_smallify(h_ChPS_final_diff_indR.at(i_trk).at(i_cent).at(i_jet), trk_itr, 1);

						
						if (jet_itr == 0 && first_pass_cent) legend_ChPS_raw_indR->AddEntry(h_ChPS_raw_ratio_indR.at(i_trk).at(i_cent).at(i_jet),trk_label.c_str(),"lp");
						if (jet_itr == 0 && first_pass_cent) legend_ChPS_raw_subtr_indR->AddEntry(h_ChPS_raw_subtr_ratio_indR.at(i_trk).at(i_cent).at(i_jet),trk_label.c_str(),"lp");
						if (jet_itr == 0 && first_pass_cent) legend_ChPS_raw_subtr_unf_indR->AddEntry(h_ChPS_raw_subtr_unf_ratio_indR.at(i_trk).at(i_cent).at(i_jet),trk_label.c_str(),"lp");
						if (jet_itr == 0 && first_pass_cent) legend_ChPS_final_indR->AddEntry(h_ChPS_final_ratio_indR.at(i_trk).at(i_cent).at(i_jet),trk_label.c_str(),"lp");
						if (jet_itr == 0 && first_pass_cent) legend_ChPS_final_diff_indR->AddEntry(h_ChPS_final_diff_indR.at(i_trk).at(i_cent).at(i_jet),trk_label.c_str(),"lp");

						h_ChPS_raw_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(0, r_max_range);
						h_ChPS_raw_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->GetYaxis()->SetRangeUser(ratio_lo, ratio_hi);
						h_ChPS_raw_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->GetYaxis()->SetNdivisions(504);

						h_ChPS_raw_subtr_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(0, r_max_range);
						h_ChPS_raw_subtr_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->GetYaxis()->SetRangeUser(ratio_lo, ratio_hi);
						h_ChPS_raw_subtr_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->GetYaxis()->SetNdivisions(504);
						
						h_ChPS_raw_subtr_unf_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(0, r_max_range);
						h_ChPS_raw_subtr_unf_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->GetYaxis()->SetRangeUser(ratio_lo, ratio_hi);
						h_ChPS_raw_subtr_unf_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->GetYaxis()->SetNdivisions(504);
						
						h_ChPS_final_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(0, r_max_range);
						h_ChPS_final_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->GetYaxis()->SetRangeUser(ratio_lo, ratio_hi);
						h_ChPS_final_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->GetYaxis()->SetNdivisions(504);
						
						h_ChPS_final_diff_indR.at(i_trk).at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(0, r_max_range);
						h_ChPS_final_diff_indR.at(i_trk).at(i_cent).at(i_jet)->GetYaxis()->SetRangeUser(-1,9);
						h_ChPS_final_diff_indR.at(i_trk).at(i_cent).at(i_jet)->GetYaxis()->SetNdivisions(504);


						c_ChPS_raw_indR->cd(i_cent+1);
						if (trk_itr == 0) h_ChPS_raw_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->Draw("");
						else h_ChPS_raw_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->Draw("same");
						gPad->SetLogx(0);
						gPad->SetLogy(0);

						c_ChPS_raw_subtr_indR->cd(i_cent+1);
						if (trk_itr == 0) h_ChPS_raw_subtr_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->Draw("");
						else h_ChPS_raw_subtr_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->Draw("same");
						gPad->SetLogx(0);
						gPad->SetLogy(0);
						
						c_ChPS_raw_subtr_unf_indR->cd(i_cent+1);
						if (trk_itr == 0) h_ChPS_raw_subtr_unf_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->Draw("");
						else h_ChPS_raw_subtr_unf_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->Draw("same");
						gPad->SetLogx(0);
						gPad->SetLogy(0);
						
						c_ChPS_final_indR->cd(i_cent+1);
						if (trk_itr == 0) h_ChPS_final_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->Draw("");
						else h_ChPS_final_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->Draw("same");
						gPad->SetLogx(0);
						gPad->SetLogy(0);

						c_ChPS_final_diff_indR->cd(i_cent+1);
						if (trk_itr == 0) h_ChPS_final_diff_indR.at(i_trk).at(i_cent).at(i_jet)->Draw("");
						else h_ChPS_final_diff_indR.at(i_trk).at(i_cent).at(i_jet)->Draw("same");
						gPad->SetLogx(0);
						gPad->SetLogy(0);

						trk_itr++;
						
					} // end trk loop
					
					c_ChPS_raw_indR->cd(i_cent+1);
					ltx->SetTextAlign(32);
					ltx->SetTextSize(12);
					ltx->DrawLatexNDC(0.93, 0.90, Form("%s", jet_label.c_str()));
					ltx->DrawLatexNDC(0.93, 0.85, Form("%s", centrality.c_str()));
					line->DrawLine(0, 1, r_max_range, 1);
					legend_ChPS_raw_indR->Draw();

					c_ChPS_raw_subtr_indR->cd(i_cent+1);
					ltx->SetTextAlign(32);
					ltx->SetTextSize(12);
					ltx->DrawLatexNDC(0.93, 0.90, Form("%s", jet_label.c_str()));
					ltx->DrawLatexNDC(0.93, 0.85, Form("%s", centrality.c_str()));
					line->DrawLine(0, 1, r_max_range, 1);
					legend_ChPS_raw_subtr_indR->Draw();
					
					c_ChPS_raw_subtr_unf_indR->cd(i_cent+1);
					ltx->SetTextAlign(32);
					ltx->SetTextSize(12);
					ltx->DrawLatexNDC(0.93, 0.90, Form("%s", jet_label.c_str()));
					ltx->DrawLatexNDC(0.93, 0.85, Form("%s", centrality.c_str()));
					line->DrawLine(0, 1, r_max_range, 1);
					legend_ChPS_raw_subtr_unf_indR->Draw();
					
					c_ChPS_final_indR->cd(i_cent+1);
					ltx->SetTextAlign(32);
					ltx->SetTextSize(12);
					ltx->DrawLatexNDC(0.93, 0.90, Form("%s", jet_label.c_str()));
					ltx->DrawLatexNDC(0.93, 0.85, Form("%s", centrality.c_str()));
					line->DrawLine(0, 1, r_max_range, 1);
					legend_ChPS_final_indR->Draw();

					c_ChPS_final_diff_indR->cd(i_cent+1);
					ltx->SetTextAlign(32);
					ltx->SetTextSize(12);
					ltx->DrawLatexNDC(0.93, 0.90, Form("%s", jet_label.c_str()));
					ltx->DrawLatexNDC(0.93, 0.85, Form("%s", centrality.c_str()));
					line->DrawLine(0, 0, r_max_range, 0);
					legend_ChPS_final_diff_indR->Draw();

					first_pass_cent = false;
				} //end cent loop
				
				pdf_label = "";
				if (i_jet == jet_pt_start) pdf_label = "(";
				if (i_jet == jet_pt_end-1) pdf_label = ")";
				c_ChPS_raw_indR->Print(Form("output_pdf_%s/ChPS_raw_ratio_dR_%s.pdf%s", sys_path.c_str(), did.c_str(), pdf_label.c_str()), Form("Title:jetpt%i", i_jet));
				c_ChPS_raw_subtr_indR->Print(Form("output_pdf_%s/ChPS_raw_subtr_ratio_dR_%s.pdf%s", sys_path.c_str(), did.c_str(), pdf_label.c_str()), Form("Title:jetpt%i", i_jet));
				c_ChPS_raw_subtr_unf_indR->Print(Form("output_pdf_%s/ChPS_raw_subtr_unf_ratio_dR_%s.pdf%s", sys_path.c_str(), did.c_str(), pdf_label.c_str()), Form("Title:jetpt%i", i_jet));
				c_ChPS_final_indR->Print(Form("output_pdf_%s/ChPS_final_ratio_dR_%s.pdf%s", sys_path.c_str(), did.c_str(), pdf_label.c_str()), Form("Title:jetpt%i", i_jet));
				c_ChPS_final_diff_indR->Print(Form("output_pdf_%s/ChPS_final_diff_dR_%s.pdf%s", sys_path.c_str(), did.c_str(), pdf_label.c_str()), Form("Title:jetpt%i", i_jet));

				jet_itr++;
			} //end jet loop
		}




		{
			cout << "Doing pre-post unfolding RDpt in dR" << endl;
			
			TCanvas *c_ChPS_fol_unf_indR = new TCanvas("c_ChPS_fol_unf_indR","c_ChPS_fol_unf_indR",900,600);
			TLegend *legend_ChPS_fol_unf_indR = new TLegend(0.19, 0.60, 0.40, 0.92, "","brNDC");
			legend_ChPS_fol_unf_indR->SetTextFont(43);
			legend_ChPS_fol_unf_indR->SetBorderSize(0);
			legend_ChPS_fol_unf_indR->SetTextSize(12);
			
			int jet_itr = 0;
			for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
			{
				string jet_label = Form("%1.0f < #it{p}_{T}^{jet} < %1.0f", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));
				
				
				int trk_itr = 0;
				for (int i_trk = 0; i_trk < N_trkpt; i_trk++)
				{
					if (i_trk < 2 || i_trk > 8) continue;
					string trk_label = Form("%1.1f < #it{p}_{T}^{trk} < %1.1f", trkpT_binning->GetBinLowEdge(i_trk+1), trkpT_binning->GetBinUpEdge(i_trk+1));
					
					//				if (i_trk != 2 &&
					//					i_trk != 5 &&
					//					i_trk != 7 &&
					//					i_trk != 8) continue;
					//1, 5, 10, 40
					
					c_ChPS_fol_unf_indR->cd();
					c_ChPS_fol_unf_indR->Clear();
					c_ChPS_fol_unf_indR->Divide(3,2);
					
					bool first_pass_cent = true;
					for (int i_cent = 0; i_cent < 6; i_cent++)
					{
						string centrality = num_to_cent(centrality_scheme,i_cent);
						
						
						SetHStyle_smallify(h_ChPS_raw_subtr_ratio_indR.at(i_trk).at(i_cent).at(i_jet), 0, 1);
						SetHStyle_smallify(h_ChPS_raw_subtr_unf_ratio_indR.at(i_trk).at(i_cent).at(i_jet), 1, 1);
						SetHStyle_smallify(h_ChPS_final_ratio_indR.at(i_trk).at(i_cent).at(i_jet), 2, 1);
						
						
						if (jet_itr == 0 && trk_itr == 0 && first_pass_cent) legend_ChPS_fol_unf_indR->AddEntry(h_ChPS_raw_subtr_ratio_indR.at(i_trk).at(i_cent).at(i_jet),"Raw+Subtr","lp");
						if (jet_itr == 0 && trk_itr == 0 && first_pass_cent) legend_ChPS_fol_unf_indR->AddEntry(h_ChPS_raw_subtr_unf_ratio_indR.at(i_trk).at(i_cent).at(i_jet),"Raw+Subtr+Unf","lp");
						if (jet_itr == 0 && trk_itr == 0 && first_pass_cent) legend_ChPS_fol_unf_indR->AddEntry(h_ChPS_final_ratio_indR.at(i_trk).at(i_cent).at(i_jet),"Raw+Subtr+Unf+BbB","lp");
						
						h_ChPS_raw_subtr_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(0, r_max_range);
						h_ChPS_raw_subtr_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->GetYaxis()->SetRangeUser(ratio_lo, ratio_hi);
						h_ChPS_raw_subtr_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->GetYaxis()->SetNdivisions(504);
						
						h_ChPS_raw_subtr_unf_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(0, r_max_range);
						h_ChPS_raw_subtr_unf_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->GetYaxis()->SetRangeUser(ratio_lo, ratio_hi);
						h_ChPS_raw_subtr_unf_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->GetYaxis()->SetNdivisions(504);
						
						h_ChPS_final_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(0, r_max_range);
						h_ChPS_final_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->GetYaxis()->SetRangeUser(ratio_lo, ratio_hi);
						h_ChPS_final_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->GetYaxis()->SetNdivisions(504);
						
						
						c_ChPS_fol_unf_indR->cd(i_cent+1);
						h_ChPS_raw_subtr_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->GetYaxis()->SetTitle(rdptr_label.c_str());
						h_ChPS_raw_subtr_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->Draw("");
						h_ChPS_raw_subtr_unf_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->Draw("same");
						h_ChPS_final_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->Draw("same");
						gPad->SetLogx(0);
						gPad->SetLogy(0);
						
						first_pass_cent = false;
						
						c_ChPS_fol_unf_indR->cd(i_cent+1);
						ltx->SetTextAlign(32);
						ltx->SetTextSize(12);
						ltx->DrawLatexNDC(0.93, 0.95, Form("%s", centrality.c_str()));
						ltx->DrawLatexNDC(0.93, 0.88, Form("%s", trk_label.c_str()));
						ltx->DrawLatexNDC(0.93, 0.80, Form("%s", jet_label.c_str()));
						line->DrawLine(0, 1, r_max_range, 1);
						legend_ChPS_fol_unf_indR->Draw();
						
					} // end cent loop
					
					
					pdf_label = "";
					if (i_jet == jet_pt_start && i_trk == 2) pdf_label = "(";
					if (i_jet == jet_pt_end-1 && i_trk == 8) pdf_label = ")";
					c_ChPS_fol_unf_indR->Print(Form("output_pdf_%s/ChPS_fol_unf_ratio_dR_%s.pdf%s", sys_path.c_str(), did.c_str(), pdf_label.c_str()), Form("Title:jet%i_trk%i", i_jet, i_trk));
					
					trk_itr++;
				} //end trk loop
				
				jet_itr++;
			} //end jet loop
		}
		
		{
			cout << "Doing in jet R < 0.4" << endl;
			
			
			TCanvas *c_ChPS_final_injet = new TCanvas("c_ChPS_final_injet","c_ChPS_final_injet",900,600);
			TLegend *legend_ChPS_final_injet = new TLegend(0.19, 0.20, 0.40, 0.30, "","brNDC");
			legend_ChPS_final_injet->SetTextFont(43);
			legend_ChPS_final_injet->SetBorderSize(0);
			legend_ChPS_final_injet->SetTextSize(12);
			
			TCanvas *c_ChPS_final_PbPb_injet = new TCanvas("c_ChPS_final_PbPb_injet","c_ChPS_final_PbPb_injet",900,600);
			TLegend *legend_ChPS_final_PbPb_injet = new TLegend(0.19, 0.20, 0.40, 0.30, "","brNDC");
			legend_ChPS_final_PbPb_injet->SetTextFont(43);
			legend_ChPS_final_PbPb_injet->SetBorderSize(0);
			legend_ChPS_final_PbPb_injet->SetTextSize(12);
			
			TCanvas *c_ChPS_final_pp_injet = new TCanvas("c_ChPS_final_pp_injet","c_ChPS_final_pp_injet",550,550);
			TLegend *legend_ChPS_final_pp_injet = new TLegend(0.19, 0.20, 0.40, 0.30, "","brNDC");
			legend_ChPS_final_pp_injet->SetTextFont(43);
			legend_ChPS_final_pp_injet->SetBorderSize(0);
			legend_ChPS_final_pp_injet->SetTextSize(12);
			
			
			int jet_itr = 0;
			for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
			{
				string jet_label = Form("%1.0f < #it{p}_{T}^{jet} < %1.0f GeV", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));
				
				bool first_pass_cent = true;
				
				c_ChPS_final_injet->cd();
				c_ChPS_final_injet->Clear();
				c_ChPS_final_injet->Divide(3,2);
				
				c_ChPS_final_PbPb_injet->cd();
				c_ChPS_final_PbPb_injet->Clear();
				c_ChPS_final_PbPb_injet->Divide(3,2);
				
				c_ChPS_final_pp_injet->Clear();
				c_ChPS_final_pp_injet->Divide(1,2);
				
				for (int i_cent = 0; i_cent < 6; i_cent++)
				{
					string centrality = num_to_cent(centrality_scheme,i_cent);
					
					if (first_pass_cent && jet_itr == 0)
					{
						legend_ChPS_final_injet->AddEntry(g_ChPS_final_ratio_injet.at(i_cent).at(i_jet),"Trk jet correlations","lp");
						legend_ChPS_final_injet->AddEntry(g_FF_final_ratio_injet.at(i_cent).at(i_jet),"FF Analysis","lp");
						
						legend_ChPS_final_PbPb_injet->AddEntry(g_ChPS_final_PbPb_injet.at(i_cent).at(i_jet),"Trk jet correlations","lp");
						legend_ChPS_final_PbPb_injet->AddEntry(g_FF_final_PbPb_injet.at(i_cent).at(i_jet),"FF Analysis","lp");
					}
					
					
					g_ChPS_final_PbPb_injet.at(i_cent).at(i_jet)->GetYaxis()->SetTitle(Form("%s (r #leq 0.4)", rdptr_label.c_str()));
					g_ChPS_final_ratio_injet.at(i_cent).at(i_jet)->GetYaxis()->SetTitle(Form("%s (r #leq 0.4)", rdptr_label.c_str()));
					g_FF_final_PbPb_injet.at(i_cent).at(i_jet)->GetYaxis()->SetTitle(dptr_label.c_str());
					g_FF_final_ratio_injet.at(i_cent).at(i_jet)->GetYaxis()->SetTitle(rdptr_label.c_str());
					g_FF_final_ratio_sys_injet.at(i_cent).at(i_jet)->GetYaxis()->SetTitle(Form("%s (r #leq 0.4)", rdptr_label.c_str()));
					
					g_FF_final_ratio_sys_injet.at(i_cent).at(i_jet)->GetXaxis()->SetTitle("#it{p}_{T}^{trk} [GeV]");
					g_ChPS_final_ratio_injet.at(i_cent).at(i_jet)->GetXaxis()->SetTitle("#it{p}_{T}^{trk} [GeV]");
					g_FF_final_ratio_injet.at(i_cent).at(i_jet)->GetXaxis()->SetTitle("#it{p}_{T}^{trk} [GeV]");
					g_ChPS_final_PbPb_injet.at(i_cent).at(i_jet)->GetYaxis()->SetTitle("#it{p}_{T}^{trk} [GeV]");
					g_FF_final_PbPb_injet.at(i_cent).at(i_jet)->GetYaxis()->SetTitle("#it{p}_{T}^{trk} [GeV]");
					
					
					//ratio
					g_FF_final_ratio_sys_injet.at(i_cent).at(i_jet)->GetXaxis()->SetLimits(1,400);
					g_FF_final_ratio_sys_injet.at(i_cent).at(i_jet)->GetYaxis()->SetNdivisions(504);
					g_FF_final_ratio_sys_injet.at(i_cent).at(i_jet)->GetYaxis()->SetRangeUser(0.4,2.1);
					
					
					c_ChPS_final_injet->cd(i_cent+1);
					g_FF_final_ratio_sys_injet.at(i_cent).at(i_jet)->Draw("a PE2");
					g_FF_final_ratio_injet.at(i_cent).at(i_jet)->Draw("p same");
					g_ChPS_final_ratio_injet.at(i_cent).at(i_jet)->Draw("p same");
					gPad->SetLogx();
					gPad->SetLogy(0);
					ltx->SetTextAlign(32);
					ltx->SetTextSize(14);
					ltx->DrawLatexNDC(0.93, 0.90, Form("%s", centrality.c_str()));
					ltx->DrawLatexNDC(0.93, 0.82, Form("%s", jet_label.c_str()));
					line->DrawLine(1, 1, 400, 1);
					if (i_cent == 0) legend_ChPS_final_injet->Draw();
					
					
					
					
					//fits
					double par[6];
					par[0] = 10.807328; par[1] = 28.553203; par[2] = 1.702105; par[3] = -0.002884; par[4] = 0.830702; par[5] = -0.075662;
					
					name = Form("fit_FF_final_PbPb_injet_cent%i_jet%i", i_cent, i_jet);
					fit_FF_final_PbPb_injet.at(i_cent).at(i_jet) = new TF1(name.c_str(),"[0] * (pow((1+[3]*x),[1]) / pow((1+[4]*x),[2])) * exp(-[5]*x)");
					fit_FF_final_PbPb_injet.at(i_cent).at(i_jet)->SetRange(1,200);
					fit_FF_final_PbPb_injet.at(i_cent).at(i_jet)->SetParameters(par);
					g_FF_final_PbPb_injet.at(i_cent).at(i_jet)->Fit(fit_FF_final_PbPb_injet.at(i_cent).at(i_jet),"RQ0","");
					fit_FF_final_PbPb_injet.at(i_cent).at(i_jet)->SetLineColor(g_FF_final_PbPb_injet.at(i_cent).at(i_jet)->GetMarkerColor());
					fit_FF_final_PbPb_injet.at(i_cent).at(i_jet)->SetLineWidth(1);
					
					
					
					fit_quality_FF_PbPb_injet.at(i_cent).at(i_jet) = (TGraphErrors*)g_FF_final_PbPb_injet.at(i_cent).at(i_jet)->Clone(Form("fit_qual_FF_PbPb_injet_cent%i_jet%i", i_cent, i_jet));
					for (int i = 0; i < g_FF_final_PbPb_injet.at(i_cent).at(i_jet)->GetN(); i++)
					{
						double x, y, fit;
						g_FF_final_PbPb_injet.at(i_cent).at(i_jet)->GetPoint(i, x, y);
						fit = fit_FF_final_PbPb_injet.at(i_cent).at(i_jet)->Eval(x);
						fit_quality_FF_PbPb_injet.at(i_cent).at(i_jet)->SetPoint(i, x, y/fit);
					}
					
					fit_quality_ChPS_PbPb_injet.at(i_cent).at(i_jet) = (TGraphErrors*)g_ChPS_final_PbPb_injet.at(i_cent).at(i_jet)->Clone(Form("fit_qual_ChPS_PbPb_injet_cent%i_jet%i", i_cent, i_jet));
					for (int i = 0; i < g_ChPS_final_PbPb_injet.at(i_cent).at(i_jet)->GetN(); i++)
					{
						double x, y, fit;
						g_ChPS_final_PbPb_injet.at(i_cent).at(i_jet)->GetPoint(i, x, y);
						fit = fit_FF_final_PbPb_injet.at(i_cent).at(i_jet)->Eval(x);
						fit_quality_ChPS_PbPb_injet.at(i_cent).at(i_jet)->SetPoint(i, x, y/fit);
					}
					
					
					//PbPb
					g_FF_final_PbPb_sys_injet.at(i_cent).at(i_jet)->GetYaxis()->SetTitle(Form("%s (r #leq 0.4)", dptr_label.c_str()));
					fit_quality_FF_PbPb_injet.at(i_cent).at(i_jet)->GetYaxis()->SetTitle("Fit Quality");
					g_FF_final_PbPb_sys_injet.at(i_cent).at(i_jet)->GetXaxis()->SetTitle("#it{p}_{T}^{trk} [GeV]");
					fit_quality_FF_PbPb_injet.at(i_cent).at(i_jet)->GetXaxis()->SetTitle("#it{p}_{T}^{trk} [GeV]");
					g_FF_final_PbPb_sys_injet.at(i_cent).at(i_jet)->GetXaxis()->SetLimits(1,400);
					fit_quality_FF_PbPb_injet.at(i_cent).at(i_jet)->GetXaxis()->SetLimits(1,400);
					g_FF_final_PbPb_sys_injet.at(i_cent).at(i_jet)->GetYaxis()->SetNdivisions(504);
					fit_quality_FF_PbPb_injet.at(i_cent).at(i_jet)->GetYaxis()->SetNdivisions(504);
					g_FF_final_PbPb_sys_injet.at(i_cent).at(i_jet)->GetYaxis()->SetRangeUser(1E-7,1E2);
					fit_quality_FF_PbPb_injet.at(i_cent).at(i_jet)->GetYaxis()->SetRangeUser(0.65,1.25);
					
					
					SetHStyle_graph_smallify(g_ChPS_final_ratio_injet.at(i_cent).at(i_jet), 0, 1);
					SetHStyle_graph_smallify(g_FF_final_ratio_injet.at(i_cent).at(i_jet), 1, 1);
					SetHStyle_graph_smallify(g_FF_final_ratio_sys_injet.at(i_cent).at(i_jet), 1, 1);
					
					SetHStyle_graph_smallify(g_ChPS_final_PbPb_injet.at(i_cent).at(i_jet), 0, 1);
					SetHStyle_graph_smallify(g_FF_final_PbPb_injet.at(i_cent).at(i_jet), 1, 1);
					SetHStyle_graph_smallify(g_FF_final_PbPb_sys_injet.at(i_cent).at(i_jet), 1, 1);
					
					SetHStyle_graph_smallify(fit_quality_ChPS_PbPb_injet.at(i_cent).at(i_jet), 0, 1);
					SetHStyle_graph_smallify(fit_quality_FF_PbPb_injet.at(i_cent).at(i_jet), 1, 1);
					
					c_ChPS_final_PbPb_injet->cd(i_cent+1);
					gPad->Divide(1,2);
					gPad->cd(1);
					gPad->SetPad(0,0.40,0.95,0.95);
					gPad->SetTopMargin(0.05);
					gPad->SetBottomMargin(0);
					gPad->SetRightMargin(0);
					g_FF_final_PbPb_sys_injet.at(i_cent).at(i_jet)->Draw("a PE2");
					g_FF_final_PbPb_injet.at(i_cent).at(i_jet)->Draw("p same");
					fit_FF_final_PbPb_injet.at(i_cent).at(i_jet)->Draw("same");
					g_ChPS_final_PbPb_injet.at(i_cent).at(i_jet)->Draw("p same");
					gPad->SetLogx();
					gPad->SetLogy();
					
					ltx->SetTextAlign(32);
					ltx->SetTextSize(14);
					ltx->DrawLatexNDC(0.93, 0.90, Form("%s", centrality.c_str()));
					ltx->DrawLatexNDC(0.93, 0.82, Form("%s", jet_label.c_str()));
					if (i_cent == 0) legend_ChPS_final_PbPb_injet->Draw();
					
					c_ChPS_final_PbPb_injet->cd(i_cent+1);
					gPad->cd(2);
					gPad->SetPad(0,0.0,0.95,0.40);
					gPad->SetTopMargin(0);
					gPad->SetBottomMargin(0.30);
					gPad->SetRightMargin(0);
					fit_quality_FF_PbPb_injet.at(i_cent).at(i_jet)->Draw("ap");
					fit_quality_ChPS_PbPb_injet.at(i_cent).at(i_jet)->Draw("p same");
					gPad->SetLogx();
					line->DrawLine(1, 1, 400, 1);
					
					//pp
					if (i_cent == 0)
					{
						
						SetHStyle_graph_smallify(g_ChPS_final_pp_injet.at(6).at(i_jet), 0, 1);
						SetHStyle_graph_smallify(g_FF_final_pp_injet.at(6).at(i_jet), 1, 1);
						SetHStyle_graph_smallify(g_FF_final_pp_sys_injet.at(6).at(i_jet), 1, 1);
						
						if (first_pass_cent && jet_itr == 0)
						{
							legend_ChPS_final_pp_injet->AddEntry(g_ChPS_final_pp_injet.at(6).at(i_jet),"Trk jet correlations","lp");
							legend_ChPS_final_pp_injet->AddEntry(g_FF_final_pp_injet.at(6).at(i_jet),"FF Analysis","lp");
						}
						par[0] = 3.394202; par[1] = 97.809290; par[2] = 1.702168; par[3] = -0.001209; par[4] = 0.295487; par[5] = -0.095548;
						
						name = Form("fit_FF_final_pp_injet_cent%i_jet%i", 6, i_jet);
						fit_FF_final_pp_injet.at(6).at(i_jet) = new TF1(name.c_str(),"[0] * (pow((1+[3]*x),[1]) / pow((1+[4]*x),[2])) * exp(-[5]*x)");
						fit_FF_final_pp_injet.at(6).at(i_jet)->SetRange(1,200);
						fit_FF_final_pp_injet.at(6).at(i_jet)->SetParameters(par);
						g_FF_final_pp_injet.at(6).at(i_jet)->Fit(fit_FF_final_pp_injet.at(6).at(i_jet),"RQ0","");
						fit_FF_final_pp_injet.at(6).at(i_jet)->SetLineColor(g_FF_final_pp_injet.at(6).at(i_jet)->GetMarkerColor());
						
						fit_quality_FF_pp_injet.at(6).at(i_jet) = (TGraphErrors*)g_FF_final_pp_injet.at(6).at(i_jet)->Clone(Form("fit_qual_FF_pp_injet_cent%i_jet%i", 6, i_jet));
						for (int i = 0; i < g_FF_final_pp_injet.at(6).at(i_jet)->GetN(); i++)
						{
							double x, y, fit;
							g_FF_final_pp_injet.at(6).at(i_jet)->GetPoint(i, x, y);
							fit = fit_FF_final_pp_injet.at(6).at(i_jet)->Eval(x);
							fit_quality_FF_pp_injet.at(6).at(i_jet)->SetPoint(i, x, y/fit);
						}
						
						fit_quality_ChPS_pp_injet.at(6).at(i_jet) = (TGraphErrors*)g_ChPS_final_pp_injet.at(6).at(i_jet)->Clone(Form("fit_qual_ChPS_pp_injet_cent%i_jet%i", 6, i_jet));
						for (int i = 0; i < g_ChPS_final_pp_injet.at(6).at(i_jet)->GetN(); i++)
						{
							double x, y, fit;
							g_ChPS_final_pp_injet.at(6).at(i_jet)->GetPoint(i, x, y);
							fit = fit_FF_final_pp_injet.at(6).at(i_jet)->Eval(x);
							fit_quality_ChPS_pp_injet.at(6).at(i_jet)->SetPoint(i, x, y/fit);
						}
						
						g_FF_final_pp_sys_injet.at(6).at(i_jet)->GetYaxis()->SetTitle(Form("%s (#it{r} #leq 0.4)", dptr_label.c_str()));
						fit_quality_FF_pp_injet.at(6).at(i_jet)->GetYaxis()->SetTitle("Fit Quality");
						g_FF_final_pp_sys_injet.at(6).at(i_jet)->GetXaxis()->SetTitle("#it{p}_{T}^{trk} [GeV]");
						fit_quality_FF_pp_injet.at(6).at(i_jet)->GetXaxis()->SetTitle("#it{p}_{T}^{trk} [GeV]");
						
						
						g_FF_final_pp_sys_injet.at(6).at(i_jet)->GetXaxis()->SetLimits(1,400);
						fit_quality_FF_pp_injet.at(6).at(i_jet)->GetXaxis()->SetLimits(1,400);
						g_FF_final_pp_sys_injet.at(6).at(i_jet)->GetYaxis()->SetNdivisions(504);
						fit_quality_FF_pp_injet.at(6).at(i_jet)->GetYaxis()->SetNdivisions(504);
						g_FF_final_pp_sys_injet.at(6).at(i_jet)->GetYaxis()->SetRangeUser(1E-7,1E2);
						fit_quality_FF_pp_injet.at(6).at(i_jet)->GetYaxis()->SetRangeUser(0.45,1.55);
						
						
						
						c_ChPS_final_pp_injet->cd()->cd(1);
						gPad->SetPad(0,0.40,0.95,0.95);
						gPad->SetTopMargin(0.05);
						gPad->SetBottomMargin(0);
						gPad->SetRightMargin(0);
						g_FF_final_pp_sys_injet.at(6).at(i_jet)->Draw("a PE2");
						g_FF_final_pp_injet.at(6).at(i_jet)->Draw("p same");
						fit_FF_final_pp_injet.at(6).at(i_jet)->Draw("p same");
						g_ChPS_final_pp_injet.at(6).at(i_jet)->Draw("p same");
						gPad->SetLogx();
						gPad->SetLogy();
						
						ltx->SetTextAlign(32);
						ltx->SetTextSize(12);
						ltx->DrawLatexNDC(0.93, 0.90, "Inclusive");
						ltx->DrawLatexNDC(0.93, 0.82, Form("%s", jet_label.c_str()));
						legend_ChPS_final_pp_injet->Draw();
						
						c_ChPS_final_pp_injet->cd()->cd(2);
						gPad->SetPad(0,0.0,0.95,0.40);
						gPad->SetTopMargin(0);
						gPad->SetBottomMargin(0.30);
						gPad->SetRightMargin(0);
						fit_quality_FF_pp_injet.at(6).at(i_jet)->Draw("ap");
						fit_quality_ChPS_pp_injet.at(6).at(i_jet)->Draw("p same");
						gPad->SetLogx();
						line->DrawLine(1, 1, 400, 1);
						
					}
					
					first_pass_cent = false;
					
				} // end cent loop
				
				pdf_label = "";
				if (i_jet == jet_pt_start) pdf_label = "(";
				if (i_jet == jet_pt_end-1) pdf_label = ")";
				
				c_ChPS_final_injet->Print(Form("output_pdf_%s/ChPS_FF_final_injet_ratio_%s.pdf%s", sys_path.c_str(), did.c_str(), pdf_label.c_str()));
				c_ChPS_final_PbPb_injet->Print(Form("output_pdf_%s/PbPb/ChPS_FF_final_injet_PbPb_%s.pdf%s", sys_path.c_str(), did.c_str(), pdf_label.c_str()));
				c_ChPS_final_pp_injet->Print(Form("output_pdf_%s/pp/ChPS_FF_final_injet_pp_%s.pdf%s", sys_path.c_str(), did.c_str(), pdf_label.c_str()));
				jet_itr++;
				
			} //end jet loop
			
		}
		
		
		
		{
			cout << "Doing Final ChPS ratio (PbPb/pp) plots in jetpt" << endl;
			
			TCanvas *c_ChPS_raw_subtr = new TCanvas("c_ChPS_raw_subtr","c_ChPS_raw_subtr",900,600);
			TLegend *legend_ChPS_raw_subtr = new TLegend(0.19, 0.80, 0.40, 0.92, "","brNDC");
			legend_ChPS_raw_subtr->SetTextFont(43);
			legend_ChPS_raw_subtr->SetBorderSize(0);
			legend_ChPS_raw_subtr->SetTextSize(10);
			
			TCanvas *c_ChPS_raw_subtr_unf = new TCanvas("c_ChPS_raw_subtr_unf","c_ChPS_raw_subtr_unf",900,600);
			TLegend *legend_ChPS_raw_subtr_unf = new TLegend(0.19, 0.80, 0.40, 0.92, "","brNDC");
			legend_ChPS_raw_subtr_unf->SetTextFont(43);
			legend_ChPS_raw_subtr_unf->SetBorderSize(0);
			legend_ChPS_raw_subtr_unf->SetTextSize(10);
			
			TCanvas *c_ChPS_final = new TCanvas("c_ChPS_final","c_ChPS_final",900,600);
			TLegend *legend_ChPS_final = new TLegend(0.19, 0.80, 0.40, 0.92, "","brNDC");
			legend_ChPS_final->SetTextFont(43);
			legend_ChPS_final->SetBorderSize(0);
			legend_ChPS_final->SetTextSize(10);
			
			
			for (int i_dR = 0; i_dR < N_dR; i_dR++)
			{
				string dr_label = Form("%1.2f < #it{r} < %1.2f", dR_binning->GetBinLowEdge(i_dR+1), dR_binning->GetBinUpEdge(i_dR+1));
				
				c_ChPS_raw_subtr->cd();
				c_ChPS_raw_subtr->Clear();
				c_ChPS_raw_subtr->Divide(3,2);
				
				c_ChPS_raw_subtr_unf->cd();
				c_ChPS_raw_subtr_unf->Clear();
				c_ChPS_raw_subtr_unf->Divide(3,2);
				
				c_ChPS_final->cd();
				c_ChPS_final->Clear();
				c_ChPS_final->Divide(3,2);
				
				bool first_pass_cent = true;
				for (int i_cent = 0; i_cent < 6; i_cent++)
				{
					string centrality = num_to_cent(centrality_scheme,i_cent);
					
					
					int jet_itr = 0;
					for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
					{
						string jet_label = Form("%1.0f < #it{p}_{T}^{jet} < %1.0f GeV", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));
						
						SetHStyle_smallify(h_ChPS_raw_subtr_ratio.at(i_dR).at(i_cent).at(i_jet), jet_itr, 1);
						SetHStyle_smallify(h_ChPS_raw_subtr_unf_ratio.at(i_dR).at(i_cent).at(i_jet), jet_itr, 1);
						SetHStyle_smallify(h_ChPS_final_ratio.at(i_dR).at(i_cent).at(i_jet), jet_itr, 1);
						
						
						if (i_dR == 0 && first_pass_cent) legend_ChPS_raw_subtr->AddEntry(h_ChPS_raw_subtr_ratio.at(i_dR).at(i_cent).at(i_jet),jet_label.c_str(),"lp");
						if (i_dR == 0 && first_pass_cent) legend_ChPS_raw_subtr_unf->AddEntry(h_ChPS_raw_subtr_unf_ratio.at(i_dR).at(i_cent).at(i_jet),jet_label.c_str(),"lp");
						if (i_dR == 0 && first_pass_cent) legend_ChPS_final->AddEntry(h_ChPS_final_ratio.at(i_dR).at(i_cent).at(i_jet),jet_label.c_str(),"lp");
						
						
						h_ChPS_raw_subtr_ratio.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetRangeUser(ratio_lo, ratio_hi);
						h_ChPS_raw_subtr_ratio.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(trk_pt_lo, trk_pt_hi);
						h_ChPS_raw_subtr_ratio.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetNdivisions(504);
						
						h_ChPS_raw_subtr_unf_ratio.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetRangeUser(ratio_lo, ratio_hi);
						h_ChPS_raw_subtr_unf_ratio.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(trk_pt_lo, trk_pt_hi);
						h_ChPS_raw_subtr_unf_ratio.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetNdivisions(504);
						
						h_ChPS_final_ratio.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetRangeUser(ratio_lo, ratio_hi);
						h_ChPS_final_ratio.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(trk_pt_lo, trk_pt_hi);
						h_ChPS_final_ratio.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetNdivisions(504);
						
						
						c_ChPS_raw_subtr->cd(i_cent+1);
						if (i_jet == jet_pt_start) h_ChPS_raw_subtr_ratio.at(i_dR).at(i_cent).at(i_jet)->Draw("");
						else h_ChPS_raw_subtr_ratio.at(i_dR).at(i_cent).at(i_jet)->Draw("same");
						gPad->SetLogx();
						gPad->SetLogy(0);
						
						c_ChPS_raw_subtr_unf->cd(i_cent+1);
						if (i_jet == jet_pt_start) h_ChPS_raw_subtr_unf_ratio.at(i_dR).at(i_cent).at(i_jet)->Draw("");
						else h_ChPS_raw_subtr_unf_ratio.at(i_dR).at(i_cent).at(i_jet)->Draw("same");
						gPad->SetLogx();
						gPad->SetLogy(0);
						
						c_ChPS_final->cd(i_cent+1);
						if (i_jet == jet_pt_start) h_ChPS_final_ratio.at(i_dR).at(i_cent).at(i_jet)->Draw("");
						else h_ChPS_final_ratio.at(i_dR).at(i_cent).at(i_jet)->Draw("same");
						gPad->SetLogx();
						gPad->SetLogy(0);
						
						jet_itr++;
						
					} // end trk loop
					
					
					
					c_ChPS_raw_subtr->cd(i_cent+1);
					ltx->SetTextAlign(32);
					ltx->SetTextSize(12);
					ltx->DrawLatexNDC(0.93, 0.90, Form("%s", dr_label.c_str()));
					ltx->DrawLatexNDC(0.93, 0.85, Form("%s", centrality.c_str()));
					line->DrawLine(trk_pt_lo, 1, trk_pt_hi, 1);
					legend_ChPS_raw_subtr->Draw();
					
					c_ChPS_raw_subtr_unf->cd(i_cent+1);
					ltx->SetTextAlign(32);
					ltx->SetTextSize(12);
					ltx->DrawLatexNDC(0.93, 0.90, Form("%s", dr_label.c_str()));
					ltx->DrawLatexNDC(0.93, 0.85, Form("%s", centrality.c_str()));
					line->DrawLine(trk_pt_lo, 1, trk_pt_hi, 1);
					legend_ChPS_raw_subtr_unf->Draw();
					
					c_ChPS_final->cd(i_cent+1);
					ltx->SetTextAlign(32);
					ltx->SetTextSize(12);
					ltx->DrawLatexNDC(0.93, 0.90, Form("%s", dr_label.c_str()));
					ltx->DrawLatexNDC(0.93, 0.85, Form("%s", centrality.c_str()));
					line->DrawLine(trk_pt_lo, 1, trk_pt_hi, 1);
					legend_ChPS_final->Draw();
					
					first_pass_cent = false;
				} //end cent loop
				
				pdf_label = "";
				if (i_dR == 0) pdf_label = "(";
				if (i_dR == N_dR-1) pdf_label = ")";
				//			c_ChPS_raw_subtr->Print(Form("output_pdf_%s/ChPS_raw_subtr_ratio_%s.pdf%s", sys_path.c_str(), did.c_str(), pdf_label.c_str()), Form("Title:dR%i", i_dR));
				//			c_ChPS_raw_subtr_unf->Print(Form("output_pdf_%s/ChPS_raw_subtr_unf_ratio_%s.pdf%s", sys_path.c_str(), did.c_str(), pdf_label.c_str()), Form("Title:dR%i", i_dR));
				c_ChPS_final->Print(Form("output_pdf_%s/ChPS_final_ratio_%s.pdf%s", sys_path.c_str(), did.c_str(), pdf_label.c_str()), Form("Title:dR%i", i_dR));
			} //end dr loop



			c_ChPS_final->Clear();
			c_ChPS_final->Divide(3,2);

			bool first_pass_cent = true;

			legend_ChPS_final->Clear();
			legend_ChPS_final->SetX1NDC(0.2);
			legend_ChPS_final->SetY1NDC(0.20);
			legend_ChPS_final->SetX2NDC(0.4);
			legend_ChPS_final->SetY2NDC(0.40);

			for (int i_cent = 0; i_cent < 6; i_cent++)
			{
				c_ChPS_final->cd(i_cent+1);

				int dr_itr = 0;
				for (int i_dR = 0; i_dR < N_dR; i_dR++)
				{
					if (i_dR != 0 && i_dR != 3 && i_dR != 7 && i_dR != 1) continue;
					string dr_label = Form("%1.2f < #it{r} < %1.2f", dR_binning->GetBinLowEdge(i_dR+1), dR_binning->GetBinUpEdge(i_dR+1));

					if (i_cent == 0) legend_ChPS_final->AddEntry(h_ChPS_final_ratio.at(i_dR).at(i_cent).at(7), dr_label.c_str(), "lp");
					SetHStyle_smallify(h_ChPS_final_ratio.at(i_dR).at(i_cent).at(7), dr_itr, 1);

					h_ChPS_final_ratio.at(i_dR).at(i_cent).at(7)->GetYaxis()->SetRangeUser(0,2.3);
					h_ChPS_final_ratio.at(i_dR).at(i_cent).at(7)->Draw("same");
					line->DrawLine(trk_pt_lo,1.,trk_pt_hi,1.);
					gPad->SetLogx();
					dr_itr++;
				}

				legend_ChPS_final->Draw();
				ltx->SetTextAlign(32);
				ltx->DrawLatexNDC(0.93, 0.90, Form("%s", num_to_cent(31, i_cent).c_str()));
			}

			c_ChPS_final->cd(1);
			ltx->DrawLatexNDC(0.93, 0.84, "#font[72]{ATLAS} Internal");
			ltx->DrawLatexNDC(0.93, 0.79, "Pb+Pb #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}");
			ltx->DrawLatexNDC(0.93, 0.74, "#it{pp} #sqrt{#font[12]{s}} = 5.02 TeV, 25 pb^{-1}");
			ltx->DrawLatexNDC(0.93, 0.68, Form("%1.2f < p_{T}^{Jet} < %1.2f [GeV]", jetpT_binning->GetBinLowEdge(8), jetpT_binning->GetBinUpEdge(8)));
			ltx->DrawLatexNDC(0.93, 0.63, Form("anti-#font[12]{k}_{#font[12]{t}} R=0.4"));


			c_ChPS_final->Print("tmp.pdf");

		}

		{
			cout << "Doing Final ChPS ratio (PbPb/pp) plots as function of centrality" << endl;

			TCanvas *c_ChPS_final_indR_cent = new TCanvas("c_ChPS_final_indR_cent","c_ChPS_final_indR_cent",900,600);
			TLegend *legend_ChPS_final_indR_cent = new TLegend(0.19, 0.60, 0.40, 0.92, "","brNDC");
			legend_ChPS_final_indR_cent->SetTextFont(43);
			legend_ChPS_final_indR_cent->SetBorderSize(0);
			legend_ChPS_final_indR_cent->SetTextSize(12);

			int jet_itr = 0;
			for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
			{
				string jet_label = Form("%1.0f < #it{p}_{T}^{jet} < %1.0f", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));

				c_ChPS_final_indR_cent->cd();
				c_ChPS_final_indR_cent->Clear();
				c_ChPS_final_indR_cent->Divide(3,2);

				int trk_itr = 0;
				for (int i_trk = 0; i_trk < N_trkpt; i_trk++)
				{

					if (i_trk < 2 || i_trk > 7) continue;
					string trk_label = Form("%1.1f < #it{p}_{T}^{trk} < %1.1f", trkpT_binning->GetBinLowEdge(i_trk+1), trkpT_binning->GetBinUpEdge(i_trk+1));

					for (int i_cent = 0; i_cent < 6; i_cent++)
					{
						string centrality = num_to_cent(centrality_scheme,i_cent);


						SetHStyle_smallify(h_ChPS_final_ratio_indR.at(i_trk).at(i_cent).at(i_jet), i_cent, 1);


						if (jet_itr == 0 && trk_itr == 0) legend_ChPS_final_indR_cent->AddEntry(h_ChPS_final_ratio_indR.at(i_trk).at(i_cent).at(i_jet),centrality.c_str(),"lp");

						h_ChPS_final_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(0, r_max_range);
						h_ChPS_final_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->GetYaxis()->SetRangeUser(ratio_lo, ratio_hi);
						h_ChPS_final_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->GetYaxis()->SetNdivisions(504);

						c_ChPS_final_indR_cent->cd(trk_itr+1);
						if (i_cent == 0) h_ChPS_final_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->Draw("");
						else h_ChPS_final_ratio_indR.at(i_trk).at(i_cent).at(i_jet)->Draw("same");
						gPad->SetLogx(0);
						gPad->SetLogy(0);


					} // end trk loop

					c_ChPS_final_indR_cent->cd(trk_itr+1);
					ltx->SetTextAlign(32);
					ltx->SetTextSize(12);
					ltx->DrawLatexNDC(0.93, 0.90, Form("%s", jet_label.c_str()));
					ltx->DrawLatexNDC(0.93, 0.85, Form("%s", trk_label.c_str()));
					line->DrawLine(0, 1, r_max_range, 1);
					legend_ChPS_final_indR_cent->Draw();
					trk_itr++;

				} //end cent loop

				pdf_label = "";
				if (i_jet == jet_pt_start) pdf_label = "(";
				if (i_jet == jet_pt_end-1) pdf_label = ")";
				c_ChPS_final_indR_cent->Print(Form("output_pdf_%s/ChPS_final_ratio_cent_dR_%s.pdf%s", sys_path.c_str(), did.c_str(), pdf_label.c_str()), Form("Title:jetpt%i", i_jet));

				jet_itr++;
			} //end jet loop
		}


		{
			cout << "Doing Final ChPS plots" << endl;

			TCanvas *c_ChPS_final_indR = new TCanvas("c_ChPS_final_indR","c_ChPS_final_indR",1200,600);
			TLegend *legend_ChPS_final_indR = new TLegend(0.19, 0.20, 0.40, 0.35, "","brNDC");
			legend_ChPS_final_indR->SetTextFont(43);
			legend_ChPS_final_indR->SetBorderSize(0);
			legend_ChPS_final_indR->SetTextSize(12);

			int jet_itr = 0;
			for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
			{
				string jet_label = Form("%1.0f < #it{p}_{T}^{jet} < %1.0f", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));

				c_ChPS_final_indR->cd();
				c_ChPS_final_indR->Clear();
				c_ChPS_final_indR->Divide(3,2);


				for (int i_cent = 0; i_cent < 6; i_cent++)
				{
					string centrality = num_to_cent(centrality_scheme,i_cent);

					int trk_itr = 0;

					for (int i_trk = 0; i_trk < N_trkpt; i_trk++)
					{
						if (i_trk < 2 || i_trk > 5) continue;
						string trk_label = Form("%1.1f < #it{p}_{T}^{trk} < %1.1f", trkpT_binning->GetBinLowEdge(i_trk+1), trkpT_binning->GetBinUpEdge(i_trk+1));


						SetHStyle_smallify(h_ChPS_final_PbPb_indR.at(i_trk).at(i_cent).at(i_jet), trk_itr, 1);
						SetHStyle_open_smallify(h_ChPS_final_pp_indR.at(i_trk).at(6).at(i_jet), trk_itr, 1);


						if (jet_itr == 0 && i_cent == 0) legend_ChPS_final_indR->AddEntry(h_ChPS_final_PbPb_indR.at(i_trk).at(i_cent).at(i_jet),trk_label.c_str(),"lp");

						h_ChPS_final_PbPb_indR.at(i_trk).at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(0, r_max_range);
						h_ChPS_final_PbPb_indR.at(i_trk).at(i_cent).at(i_jet)->GetYaxis()->SetRangeUser(1E-3, 1E2);
						h_ChPS_final_PbPb_indR.at(i_trk).at(i_cent).at(i_jet)->GetYaxis()->SetNdivisions(504);

						c_ChPS_final_indR->cd(i_cent+1);
						if (i_trk == 0) h_ChPS_final_PbPb_indR.at(i_trk).at(i_cent).at(i_jet)->Draw("");
						else h_ChPS_final_PbPb_indR.at(i_trk).at(i_cent).at(i_jet)->Draw("same");
						h_ChPS_final_pp_indR.at(i_trk).at(6).at(i_jet)->Draw("same");
						gPad->SetLogx(0);
						gPad->SetLogy();

						trk_itr++;

					} // end trk loop

					c_ChPS_final_indR->cd(i_cent+1);
					ltx->SetTextAlign(32);
					ltx->SetTextSize(12);
					ltx->DrawLatexNDC(0.93, 0.90, Form("%s", jet_label.c_str()));
					ltx->DrawLatexNDC(0.93, 0.85, Form("%s", centrality.c_str()));
//					line->DrawLine(0, 1, r_max_range, 1);
					legend_ChPS_final_indR->Draw();

				} //end cent loop

				pdf_label = "";
				if (i_jet == jet_pt_start) pdf_label = "(";
				if (i_jet == jet_pt_end-1) pdf_label = ")";
				c_ChPS_final_indR->Print(Form("output_pdf_%s/ChPS_final_dR_%s.pdf%s", sys_path.c_str(), did.c_str(), pdf_label.c_str()), Form("Title:jetpt%i", i_jet));

				jet_itr++;
			} //end jet loop
		}



	}
	cout << "######### DONE COMP_ChPS #########" << endl << endl;;


}

