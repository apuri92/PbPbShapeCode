#include "../functions/global_variables.h"

TGraphAsymmErrors* shift(TGraphAsymmErrors* g, int variable, double shift_size = 0.0025)
{
	for (int i = 0; i < g->GetN(); i++)
	{
		double x, y, err_x;
		g->GetPoint(i, x, y);
		err_x = g->GetErrorX(i);
		double orig_pos = x+variable*shift_size;// * err_x/2;
		g->SetPoint(i, orig_pos, y);
	}

	return g;
}

TGraphAsymmErrors* shift(TGraphAsymmErrors* g, TAxis* axis, int variable, double shift_size = 0.0025)
{
	for (int i = 0; i < g->GetN(); i++)
	{
		double x, y;
		g->GetPoint(i, x, y);
		double binW = axis->GetBinWidth(axis->FindBin(x));
		double orig_pos = x+variable*shift_size * binW;
		g->SetPoint(i, orig_pos, y);
	}

	return g;
}



void draw_conf_plots()
{
	cout << "######### DOING CONF_Plots #########" << endl;

	SetAtlasStyle();
	gErrorIgnoreLevel = 3001;
//	gStyle->SetErrorX(0);
	gStyle->SetEndErrorSize(3);

	string name;
	string rdptr_label = "#it{R}_{#it{D} (#it{p}_{T}, #it{r})}";
	string deltadptr_label = "#it{#Delta} #it{D} (#it{p}_{T}, #it{r}) [GeV^{-1}]";
	string dptr_label = "#it{D} (#it{p}_{T}, #it{r}) [GeV^{-1}]";
	string r_label = "#it{r}";
	string trk_label = "#it{p}_{T}^{trk} [GeV]";
	string jet_label = "#it{p}_{T}^{jet} [GeV]";

	int isMC = 0; //will never run on MC for these plots
	std::string did = "data";

	TFile *f_RDpT = new TFile(Form("output_pdf_nominal/root/final_RDpT_%s.root", did.c_str()));
	TFile *f_RDpT_sys = new TFile(Form("output_pdf_nominal/root/final_RDpT_sys_%s.root", did.c_str()));
	TFile *f_DeltaDpT_sys = new TFile(Form("output_pdf_nominal/root/final_DeltaDpT_sys_%s.root", did.c_str()));
	TFile *f_ChPS_PbPb = new TFile(Form("output_pdf_nominal/root/final_ChPS_%s_PbPb.root", did.c_str()));
	TFile *f_ChPS_PbPb_sys = new TFile(Form("output_pdf_nominal/root/final_ChPS_sys_%s_PbPb.root", did.c_str()));
	TFile *f_ChPS_pp = new TFile(Form("output_pdf_nominal/root/final_ChPS_%s_pp.root", did.c_str()));
	TFile *f_ChPS_pp_sys = new TFile(Form("output_pdf_nominal/root/final_ChPS_sys_%s_pp.root", did.c_str()));


	TFile *f_output = new TFile(Form("output_pdf_nominal/hepdata/values.root"), "recreate");
	//Use this for new version of drawChPS (in compiledCode)
//	TFile *f_RDpT = new TFile(Form("compiledCode/root/output_nominal_data.root"));
//	TFile *f_RDpT_sys = new TFile(Form("compiledCode/root/final_RDpT_sys_%s.root", did.c_str()));
//	TFile *f_DeltaDpT_sys = new TFile(Form("compiledCode/root/final_DeltaDpT_sys_%s.root", did.c_str()));
//
//	TFile *f_ChPS_PbPb = new TFile(Form("compiledCode/root/output_nominal_data.root"));
//	TFile *f_ChPS_PbPb_sys = new TFile(Form("compiledCode/root/final_ChPS_sys_%s_PbPb.root", did.c_str()));
//	TFile *f_ChPS_pp = new TFile(Form("compiledCode/root/output_nominal_data.root"));
//	TFile *f_ChPS_pp_sys = new TFile(Form("compiledCode/root/final_ChPS_sys_%s_pp.root", did.c_str()));

	cout << "Using files:" << endl;
	cout << f_RDpT->GetName() << endl;
	cout << f_RDpT_sys->GetName() << endl;
	cout << f_DeltaDpT_sys->GetName() << endl;
	cout << f_ChPS_PbPb->GetName() << endl;
	cout << f_ChPS_PbPb_sys->GetName() << endl;
	cout << f_ChPS_pp->GetName() << endl;
	cout << f_ChPS_pp_sys->GetName() << endl;


	TAxis* dR_binning = (TAxis*)f_RDpT->Get("dR_binning");
	TAxis* jetpT_binning = (TAxis*)f_RDpT->Get("jetpT_binning");
	TAxis* trkpT_binning = (TAxis*)f_RDpT->Get("trkpT_binning");

	int N_dR = dR_binning->GetNbins();
	int N_jetpt = jetpT_binning->GetNbins();
	int N_trkpt = trkpT_binning->GetNbins();

	f_output->cd();
	dR_binning->Write("dR_binning");
	jetpT_binning->Write("jetpT_binning");
	trkpT_binning->Write("trkpT_binning");


	double r_max_range = 0.8;
	//indR
	vector<vector<vector<TH1*>>> h_RDpT_final_ratio_indR (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_RDpT_final_sys_Totalpos_indR (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_RDpT_final_sys_Totalneg_indR (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TGraphAsymmErrors*>>> g_RDpT_final_sys_indR (N_trkpt, vector<vector<TGraphAsymmErrors*>> (n_cent_cuts, vector<TGraphAsymmErrors*> (N_jetpt)));
	vector<vector<vector<TGraphAsymmErrors*>>> g_RDpT_final_stat_indR (N_trkpt, vector<vector<TGraphAsymmErrors*>> (n_cent_cuts, vector<TGraphAsymmErrors*> (N_jetpt)));

	vector<vector<vector<TH1*>>> h_DeltaDpT_final_indR (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_DeltaDpT_final_sys_Totalpos_indR (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_DeltaDpT_final_sys_Totalneg_indR (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TGraphAsymmErrors*>>> g_DeltaDpT_final_sys_indR (N_trkpt, vector<vector<TGraphAsymmErrors*>> (n_cent_cuts, vector<TGraphAsymmErrors*> (N_jetpt)));
	vector<vector<vector<TGraphAsymmErrors*>>> g_DeltaDpT_final_stat_indR (N_trkpt, vector<vector<TGraphAsymmErrors*>> (n_cent_cuts, vector<TGraphAsymmErrors*> (N_jetpt)));

	vector<vector<vector<TH1*>>> h_ChPS_PbPb_final_ratio_indR (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_PbPb_final_sys_Totalpos_indR (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_PbPb_final_sys_Totalneg_indR (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TGraphAsymmErrors*>>> g_ChPS_PbPb_final_sys_indR (N_trkpt, vector<vector<TGraphAsymmErrors*>> (n_cent_cuts, vector<TGraphAsymmErrors*> (N_jetpt)));
	vector<vector<vector<TGraphAsymmErrors*>>> g_ChPS_PbPb_final_stat_indR (N_trkpt, vector<vector<TGraphAsymmErrors*>> (n_cent_cuts, vector<TGraphAsymmErrors*> (N_jetpt)));

	vector<vector<vector<TH1*>>> h_ChPS_pp_final_ratio_indR (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_pp_final_sys_Totalpos_indR (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_pp_final_sys_Totalneg_indR (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TGraphAsymmErrors*>>> g_ChPS_pp_final_sys_indR (N_trkpt, vector<vector<TGraphAsymmErrors*>> (n_cent_cuts, vector<TGraphAsymmErrors*> (N_jetpt)));
	vector<vector<vector<TGraphAsymmErrors*>>> g_ChPS_pp_final_stat_indR (N_trkpt, vector<vector<TGraphAsymmErrors*>> (n_cent_cuts, vector<TGraphAsymmErrors*> (N_jetpt)));

	vector<vector<vector<TH1*>>> h_RDpT_final_ratio_inTrk (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_RDpT_final_sys_Totalpos_inTrk (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_RDpT_final_sys_Totalneg_inTrk (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TGraphAsymmErrors*>>> g_RDpT_final_sys_inTrk (N_dR, vector<vector<TGraphAsymmErrors*>> (n_cent_cuts, vector<TGraphAsymmErrors*> (N_jetpt)));
	vector<vector<vector<TGraphAsymmErrors*>>> g_RDpT_final_stat_inTrk (N_dR, vector<vector<TGraphAsymmErrors*>> (n_cent_cuts, vector<TGraphAsymmErrors*> (N_jetpt)));

	vector<vector<vector<TH1*>>> h_RDpT_final_ratio_inJet (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_trkpt)));
	vector<vector<vector<TH1*>>> h_RDpT_final_sys_Totalpos_inJet (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_trkpt)));
	vector<vector<vector<TH1*>>> h_RDpT_final_sys_Totalneg_inJet (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_trkpt)));
	vector<vector<vector<TGraphAsymmErrors*>>> g_RDpT_final_sys_inJet (N_dR, vector<vector<TGraphAsymmErrors*>> (n_cent_cuts, vector<TGraphAsymmErrors*> (N_trkpt)));
	vector<vector<vector<TGraphAsymmErrors*>>> g_RDpT_final_stat_inJet (N_dR, vector<vector<TGraphAsymmErrors*>> (n_cent_cuts, vector<TGraphAsymmErrors*> (N_trkpt)));

	string pdf_label;

	TLine *line = new TLine();
	line->SetLineColor(kBlack);
	line->SetLineStyle(3);
	TLatex *ltx = new TLatex();
	ltx->SetTextFont(43);
	ltx->SetTextSize(22);

	double trk_pt_lo = 1.;
	double trk_pt_hi = 63.;

	double jet_pt_lo = 126.;
	double jet_pt_hi = 316;

	double ratio_lo = 0.;
	double ratio_hi = 2.5;

	int jet_pt_start = 7;
	int jet_pt_end = 11;

	int trk_pt_start = 2;
	int trk_pt_end = 9;

	double opacity = 0.7;

	int trk_select1 = 2;
	int trk_select2 = 6;

	//get all rdpt, dpt (pbpb and pp) distributions and uncertainty (done as a function of r)
	for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
	{
		for (int i_cent = 0; i_cent < 6; i_cent++)
		{
			int trk_itr = 0;
			for (int i_trk = 0; i_trk < N_trkpt; i_trk++)
			{
				if (i_trk < 2 || i_trk > 9) continue;
				//RDpT
				name = Form("h_RDpT_final_indR_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_RDpT_final_ratio_indR[i_trk][i_cent][i_jet] = (TH1*)f_RDpT->Get(name.c_str());

				name = Form("h_RDpT_sys_trk%i_cent%i_jetpt%i_total_p", i_trk, i_cent, i_jet);
				h_RDpT_final_sys_Totalpos_indR[i_trk][i_cent][i_jet] = (TH1*)f_RDpT_sys->Get(name.c_str());

				name = Form("h_RDpT_sys_trk%i_cent%i_jetpt%i_total_n", i_trk, i_cent, i_jet);
				h_RDpT_final_sys_Totalneg_indR[i_trk][i_cent][i_jet] = (TH1*)f_RDpT_sys->Get(name.c_str());

				g_RDpT_final_sys_indR[i_trk][i_cent][i_jet] = new TGraphAsymmErrors(h_RDpT_final_ratio_indR[i_trk][i_cent][i_jet]);
				g_RDpT_final_sys_indR[i_trk][i_cent][i_jet]->GetYaxis()->SetTitle(rdptr_label.c_str());
				g_RDpT_final_sys_indR[i_trk][i_cent][i_jet]->GetXaxis()->SetTitle(r_label.c_str());

				g_RDpT_final_stat_indR[i_trk][i_cent][i_jet] = new TGraphAsymmErrors(h_RDpT_final_ratio_indR[i_trk][i_cent][i_jet]);
				g_RDpT_final_stat_indR[i_trk][i_cent][i_jet]->GetYaxis()->SetTitle(rdptr_label.c_str());
				g_RDpT_final_stat_indR[i_trk][i_cent][i_jet]->GetXaxis()->SetTitle(r_label.c_str());
	
				//DeltaDpT
				name = Form("h_DeltaDpT_final_indR_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_DeltaDpT_final_indR[i_trk][i_cent][i_jet] = (TH1*)f_RDpT->Get(name.c_str());

				name = Form("h_DeltaDpT_sys_trk%i_cent%i_jetpt%i_total_p", i_trk, i_cent, i_jet);
				h_DeltaDpT_final_sys_Totalpos_indR[i_trk][i_cent][i_jet] = (TH1*)f_DeltaDpT_sys->Get(name.c_str());

				name = Form("h_DeltaDpT_sys_trk%i_cent%i_jetpt%i_total_n", i_trk, i_cent, i_jet);
				h_DeltaDpT_final_sys_Totalneg_indR[i_trk][i_cent][i_jet] = (TH1*)f_DeltaDpT_sys->Get(name.c_str());

				g_DeltaDpT_final_sys_indR[i_trk][i_cent][i_jet] = new TGraphAsymmErrors(h_DeltaDpT_final_indR[i_trk][i_cent][i_jet]);
				g_DeltaDpT_final_sys_indR[i_trk][i_cent][i_jet]->GetYaxis()->SetTitle(deltadptr_label.c_str());
				g_DeltaDpT_final_sys_indR[i_trk][i_cent][i_jet]->GetXaxis()->SetTitle(r_label.c_str());

				g_DeltaDpT_final_stat_indR[i_trk][i_cent][i_jet] = new TGraphAsymmErrors(h_DeltaDpT_final_indR[i_trk][i_cent][i_jet]);
				g_DeltaDpT_final_stat_indR[i_trk][i_cent][i_jet]->GetYaxis()->SetTitle(deltadptr_label.c_str());
				g_DeltaDpT_final_stat_indR[i_trk][i_cent][i_jet]->GetXaxis()->SetTitle(r_label.c_str());


				//DpT PbPb
				//For compiledCode
//				name = Form("h_PbPb_ChPS_final_indR_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				name = Form("h_ChPS_final_indR_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_ChPS_PbPb_final_ratio_indR[i_trk][i_cent][i_jet] = (TH1*)f_ChPS_PbPb->Get(name.c_str());
				h_ChPS_PbPb_final_ratio_indR[i_trk][i_cent][i_jet]->SetName(Form("%s_PbPb", name.c_str()));

				name = Form("h_ChPS_sys_trk%i_cent%i_jetpt%i_total_p", i_trk, i_cent, i_jet);
				h_ChPS_PbPb_final_sys_Totalpos_indR[i_trk][i_cent][i_jet] = (TH1*)f_ChPS_PbPb_sys->Get(name.c_str());
				h_ChPS_PbPb_final_sys_Totalpos_indR[i_trk][i_cent][i_jet]->SetName(Form("%s_PbPb", name.c_str()));

				name = Form("h_ChPS_sys_trk%i_cent%i_jetpt%i_total_n", i_trk, i_cent, i_jet);
				h_ChPS_PbPb_final_sys_Totalneg_indR[i_trk][i_cent][i_jet] = (TH1*)f_ChPS_PbPb_sys->Get(name.c_str());
				h_ChPS_PbPb_final_sys_Totalneg_indR[i_trk][i_cent][i_jet]->SetName(Form("%s_PbPb", name.c_str()));

				g_ChPS_PbPb_final_sys_indR[i_trk][i_cent][i_jet] = new TGraphAsymmErrors(h_ChPS_PbPb_final_ratio_indR[i_trk][i_cent][i_jet]);
				g_ChPS_PbPb_final_sys_indR[i_trk][i_cent][i_jet]->GetYaxis()->SetTitle(dptr_label.c_str());
				g_ChPS_PbPb_final_sys_indR[i_trk][i_cent][i_jet]->GetXaxis()->SetTitle(r_label.c_str());

				g_ChPS_PbPb_final_stat_indR[i_trk][i_cent][i_jet] = new TGraphAsymmErrors(h_ChPS_PbPb_final_ratio_indR[i_trk][i_cent][i_jet]);
				g_ChPS_PbPb_final_stat_indR[i_trk][i_cent][i_jet]->GetYaxis()->SetTitle(dptr_label.c_str());
				g_ChPS_PbPb_final_stat_indR[i_trk][i_cent][i_jet]->GetXaxis()->SetTitle(r_label.c_str());


				//DpT pp
				if (i_cent == 0)
				{
					//For compiledCode
//					name = Form("h_pp_ChPS_final_indR_trk%i_cent%i_jetpt%i", i_trk, 6, i_jet);
					name = Form("h_ChPS_final_indR_trk%i_cent%i_jetpt%i", i_trk, 6, i_jet);
					h_ChPS_pp_final_ratio_indR[i_trk][6][i_jet] = (TH1*)f_ChPS_pp->Get(name.c_str());
					h_ChPS_pp_final_ratio_indR[i_trk][6][i_jet]->SetName(Form("%s_pp", name.c_str()));

					name = Form("h_ChPS_sys_trk%i_cent%i_jetpt%i_total_p", i_trk, 6, i_jet);
					h_ChPS_pp_final_sys_Totalpos_indR[i_trk][6][i_jet] = (TH1*)f_ChPS_pp_sys->Get(name.c_str());
					h_ChPS_pp_final_sys_Totalpos_indR[i_trk][6][i_jet]->SetName(Form("%s_pp", name.c_str()));

					name = Form("h_ChPS_sys_trk%i_cent%i_jetpt%i_total_n", i_trk, 6, i_jet);
					h_ChPS_pp_final_sys_Totalneg_indR[i_trk][6][i_jet] = (TH1*)f_ChPS_pp_sys->Get(name.c_str());
					h_ChPS_pp_final_sys_Totalneg_indR[i_trk][6][i_jet]->SetName(Form("%s_pp", name.c_str()));

					g_ChPS_pp_final_sys_indR[i_trk][6][i_jet] = new TGraphAsymmErrors(h_ChPS_pp_final_ratio_indR[i_trk][6][i_jet]);
					g_ChPS_pp_final_sys_indR[i_trk][6][i_jet]->GetYaxis()->SetTitle(dptr_label.c_str());
					g_ChPS_pp_final_sys_indR[i_trk][6][i_jet]->GetXaxis()->SetTitle(r_label.c_str());

					g_ChPS_pp_final_stat_indR[i_trk][6][i_jet] = new TGraphAsymmErrors(h_ChPS_pp_final_ratio_indR[i_trk][6][i_jet]);
					g_ChPS_pp_final_stat_indR[i_trk][6][i_jet]->GetYaxis()->SetTitle(dptr_label.c_str());
					g_ChPS_pp_final_stat_indR[i_trk][6][i_jet]->GetXaxis()->SetTitle(r_label.c_str());
				}


				double tmp = 0;
				double nom, sys_hi, sys_lo, r_position, r_width, stat_hi, stat_lo;
				for (int i_dR = 0 ; i_dR < N_dR; i_dR++)
				{
					if (i_dR > 10) continue;
					//RDpT
					r_position = h_RDpT_final_ratio_indR[i_trk][i_cent][i_jet]->GetBinCenter(i_dR+1);
					r_width = 0.020; //h_RDpT_final_ratio_indR[i_trk][i_cent][i_jet]->GetBinWidth(i_dR+1)/4;

					nom = h_RDpT_final_ratio_indR[i_trk][i_cent][i_jet]->GetBinContent(i_dR+1);
					sys_hi = nom * (h_RDpT_final_sys_Totalpos_indR[i_trk][i_cent][i_jet]->GetBinContent(i_dR+1));
					sys_lo = fabs(nom * (h_RDpT_final_sys_Totalneg_indR[i_trk][i_cent][i_jet]->GetBinContent(i_dR+1)));
					stat_hi = h_RDpT_final_ratio_indR[i_trk][i_cent][i_jet]->GetBinError(i_dR+1);
					stat_lo = h_RDpT_final_ratio_indR[i_trk][i_cent][i_jet]->GetBinError(i_dR+1);

					//almost 0 statistic of 25 GeV tracks in 200-251 gev jets in 60-80% collisions above 0.3
					if (i_trk == 8 && i_jet == 9 && i_dR >= 6 && i_cent == 5)
					{
						cout << "VERY SPECIFIC CUT: " << nom << endl;
						nom = -1;
					}
					// removing points with poor closure (this happens for high z particles at the jet edge
//					Jet pT [GeV]	Track pT [GeV]	R
//					126-158 GeV		6-10			r > 0.3 (jet7_trk6_dR6)
//									10-25 GeV		r > 0.3 (jet7_trk7_dR6)
//									25-63 GeV		r > 0.2 (jet7_trk8_dR4)
//					158-200 GeV		10-25 GeV		r > 0.4 (jet8_trk7_dR7)
//									25-63 GeV		r > 0.3 (jet8_trk8_dR6)
//					200-251 GeV		25-63 GeV		r > 0.3 (jet9_trk8_dR6)
//					251-316 GeV		x				x

					if ((i_jet == 7 && i_trk == 6 && i_dR >= 6)
						||
						(i_jet == 7 && i_trk == 7 && i_dR >= 6)
						||
						(i_jet == 7 && i_trk == 8 && i_dR >= 4)
						||
						(i_jet == 8 && i_trk == 7 && i_dR >= 7)
						||
						(i_jet == 8 && i_trk == 8 && i_dR >= 6)
						||
						(i_jet == 9 && i_trk == 8 && i_dR >= 6)
						)
					{
						nom = -1;
					}


					g_RDpT_final_stat_indR[i_trk][i_cent][i_jet]->SetPoint(i_dR, r_position, nom );
					g_RDpT_final_stat_indR[i_trk][i_cent][i_jet]->SetPointError(i_dR, 0, 0, stat_lo, stat_hi);

					g_RDpT_final_sys_indR[i_trk][i_cent][i_jet]->SetPoint(i_dR, r_position, nom );
					g_RDpT_final_sys_indR[i_trk][i_cent][i_jet]->SetPointError(i_dR, r_width/2, r_width/2, sys_lo, sys_hi);


					//DeltaDpT
					r_position = h_DeltaDpT_final_indR[i_trk][i_cent][i_jet]->GetBinCenter(i_dR+1);
					r_width = 0.020; //h_RDpT_final_ratio_indR[i_trk][i_cent][i_jet]->GetBinWidth(i_dR+1)/4;

					nom = h_DeltaDpT_final_indR[i_trk][i_cent][i_jet]->GetBinContent(i_dR+1);
					sys_hi = fabs(nom * (h_DeltaDpT_final_sys_Totalpos_indR[i_trk][i_cent][i_jet]->GetBinContent(i_dR+1)));
					sys_lo = fabs(nom * (h_DeltaDpT_final_sys_Totalneg_indR[i_trk][i_cent][i_jet]->GetBinContent(i_dR+1)));
					stat_hi = h_DeltaDpT_final_indR[i_trk][i_cent][i_jet]->GetBinError(i_dR+1);
					stat_lo = h_DeltaDpT_final_indR[i_trk][i_cent][i_jet]->GetBinError(i_dR+1);

					if ((i_jet == 7 && i_trk == 6 && i_dR >= 6)
						||
						(i_jet == 7 && i_trk == 7 && i_dR >= 6)
						||
						(i_jet == 7 && i_trk == 8 && i_dR >= 4)
						||
						(i_jet == 8 && i_trk == 7 && i_dR >= 7)
						||
						(i_jet == 8 && i_trk == 8 && i_dR >= 6)
						||
						(i_jet == 9 && i_trk == 8 && i_dR >= 6)
						)
					{
						nom = -5;
					}
					g_DeltaDpT_final_stat_indR[i_trk][i_cent][i_jet]->SetPoint(i_dR, r_position, nom );
					g_DeltaDpT_final_stat_indR[i_trk][i_cent][i_jet]->SetPointError(i_dR, 0, 0, stat_lo, stat_hi);

					g_DeltaDpT_final_sys_indR[i_trk][i_cent][i_jet]->SetPoint(i_dR, r_position, nom );
					g_DeltaDpT_final_sys_indR[i_trk][i_cent][i_jet]->SetPointError(i_dR, r_width/2, r_width/2, sys_lo, sys_hi);



					//DpT_PbPb
					r_position = h_ChPS_PbPb_final_ratio_indR[i_trk][i_cent][i_jet]->GetBinCenter(i_dR+1);
					r_width = 0.020; //h_ChPS_PbPb_final_ratio_indR[i_trk][i_cent][i_jet]->GetBinWidth(i_dR+1);

					nom = h_ChPS_PbPb_final_ratio_indR[i_trk][i_cent][i_jet]->GetBinContent(i_dR+1);

					if ((i_jet == 7 && i_trk == 6 && i_dR >= 6)
						||
						(i_jet == 7 && i_trk == 7 && i_dR >= 6)
						||
						(i_jet == 7 && i_trk == 8 && i_dR >= 4)
						||
						(i_jet == 8 && i_trk == 7 && i_dR >= 7)
						||
						(i_jet == 8 && i_trk == 8 && i_dR >= 6)
						||
						(i_jet == 9 && i_trk == 8 && i_dR >= 6)
						)
					{
						nom = 10E-9;
					}

					sys_hi = nom * (h_ChPS_PbPb_final_sys_Totalpos_indR[i_trk][i_cent][i_jet]->GetBinContent(i_dR+1));
					sys_lo = fabs(nom * (h_ChPS_PbPb_final_sys_Totalneg_indR[i_trk][i_cent][i_jet]->GetBinContent(i_dR+1)));
					stat_hi = h_ChPS_PbPb_final_ratio_indR[i_trk][i_cent][i_jet]->GetBinError(i_dR+1);
					stat_lo = h_ChPS_PbPb_final_ratio_indR[i_trk][i_cent][i_jet]->GetBinError(i_dR+1);

					g_ChPS_PbPb_final_stat_indR[i_trk][i_cent][i_jet]->SetPoint(i_dR, r_position, nom );
					g_ChPS_PbPb_final_stat_indR[i_trk][i_cent][i_jet]->SetPointError(i_dR, 0, 0, stat_lo, stat_hi);

					g_ChPS_PbPb_final_sys_indR[i_trk][i_cent][i_jet]->SetPoint(i_dR, r_position, nom );
					g_ChPS_PbPb_final_sys_indR[i_trk][i_cent][i_jet]->SetPointError(i_dR, r_width/2, r_width/2, sys_lo, sys_hi);

					
					//DpT_pp
					if (i_cent == 0)
					{
						r_position = h_ChPS_pp_final_ratio_indR[i_trk][6][i_jet]->GetBinCenter(i_dR+1);
						r_width = 0.02; //h_ChPS_pp_final_ratio_indR[i_trk][6][i_jet]->GetBinWidth(i_dR+1);

						nom = h_ChPS_pp_final_ratio_indR[i_trk][6][i_jet]->GetBinContent(i_dR+1);
						sys_hi = nom * (h_ChPS_pp_final_sys_Totalpos_indR[i_trk][6][i_jet]->GetBinContent(i_dR+1));
						sys_lo = fabs(nom * (h_ChPS_pp_final_sys_Totalneg_indR[i_trk][6][i_jet]->GetBinContent(i_dR+1)));
						stat_hi = h_ChPS_pp_final_ratio_indR[i_trk][6][i_jet]->GetBinError(i_dR+1);
						stat_lo = h_ChPS_pp_final_ratio_indR[i_trk][6][i_jet]->GetBinError(i_dR+1);

						if ((i_jet == 7 && i_trk == 6 && i_dR >= 6)
							||
							(i_jet == 7 && i_trk == 7 && i_dR >= 6)
							||
							(i_jet == 7 && i_trk == 8 && i_dR >= 4)
							||
							(i_jet == 8 && i_trk == 7 && i_dR >= 7)
							||
							(i_jet == 8 && i_trk == 8 && i_dR >= 6)
							||
							(i_jet == 9 && i_trk == 8 && i_dR >= 6)
							)
						{
							nom = 10E-9;
						}
						
						g_ChPS_pp_final_sys_indR[i_trk][6][i_jet]->SetPoint(i_dR, r_position, nom );
						g_ChPS_pp_final_sys_indR[i_trk][6][i_jet]->SetPointError(i_dR, r_width/2, r_width/2, sys_lo, sys_hi);

						g_ChPS_pp_final_stat_indR[i_trk][6][i_jet]->SetPoint(i_dR, r_position, nom );
						g_ChPS_pp_final_stat_indR[i_trk][6][i_jet]->SetPointError(i_dR, 0, 0, stat_lo, stat_hi);
					}

				}

//				if (i_jet == 7 && i_cent == 0 && i_trk == 4)
				{
					f_output->cd();



					//RDpT
					h_RDpT_final_ratio_indR[i_trk][i_cent][i_jet]->SetTitle(Form("RDpT Nominal - pTjet: %1.0f-%1.0f GeV, pTtrk: %1.2f-%1.2f GeV, centrality: %s", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1), trkpT_binning->GetBinLowEdge(i_trk+1), trkpT_binning->GetBinUpEdge(i_trk+1), num_to_cent(31, i_cent).c_str()));
					h_RDpT_final_sys_Totalpos_indR[i_trk][i_cent][i_jet]->SetTitle(Form("RDpT SysPos - pTjet: %1.0f-%1.0f GeV, pTtrk: %1.2f-%1.2f GeV, centrality: %s", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1), trkpT_binning->GetBinLowEdge(i_trk+1), trkpT_binning->GetBinUpEdge(i_trk+1), num_to_cent(31, i_cent).c_str()));
					h_RDpT_final_sys_Totalneg_indR[i_trk][i_cent][i_jet]->SetTitle(Form("RDpT SysNeg - pTjet: %1.0f-%1.0f GeV, pTtrk: %1.2f-%1.2f GeV, centrality: %s", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1), trkpT_binning->GetBinLowEdge(i_trk+1), trkpT_binning->GetBinUpEdge(i_trk+1), num_to_cent(31, i_cent).c_str()));

					h_RDpT_final_ratio_indR[i_trk][i_cent][i_jet]->SetName(Form("h_RDpT_trk%i_cent%i_jet%i", i_trk, i_cent, i_jet));
					h_RDpT_final_sys_Totalpos_indR[i_trk][i_cent][i_jet]->SetName(Form("h_RDpT_trk%i_cent%i_jet%i_sysP", i_trk, i_cent, i_jet));
					h_RDpT_final_sys_Totalneg_indR[i_trk][i_cent][i_jet]->SetName(Form("h_RDpT_trk%i_cent%i_jet%i_sysN", i_trk, i_cent, i_jet));

					h_RDpT_final_ratio_indR[i_trk][i_cent][i_jet]->Write(Form("h_RDpT_trk%i_cent%i_jet%i", i_trk, i_cent, i_jet));
					h_RDpT_final_sys_Totalpos_indR[i_trk][i_cent][i_jet]->Write(Form("h_RDpT_trk%i_cent%i_jet%i_sysP", i_trk, i_cent, i_jet));
					h_RDpT_final_sys_Totalneg_indR[i_trk][i_cent][i_jet]->Write(Form("h_RDpT_trk%i_cent%i_jet%i_sysN", i_trk, i_cent, i_jet));


					//DeltaDpT
					h_DeltaDpT_final_indR[i_trk][i_cent][i_jet]->SetTitle(Form("DeltaDpT Nominal - pTjet: %1.0f-%1.0f GeV, pTtrk: %1.2f-%1.2f GeV, centrality: %s", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1), trkpT_binning->GetBinLowEdge(i_trk+1), trkpT_binning->GetBinUpEdge(i_trk+1), num_to_cent(31, i_cent).c_str()));
					h_DeltaDpT_final_sys_Totalpos_indR[i_trk][i_cent][i_jet]->SetTitle(Form("DeltaDpT SysPos - pTjet: %1.0f-%1.0f GeV, pTtrk: %1.2f-%1.2f GeV, centrality: %s", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1), trkpT_binning->GetBinLowEdge(i_trk+1), trkpT_binning->GetBinUpEdge(i_trk+1), num_to_cent(31, i_cent).c_str()));
					h_DeltaDpT_final_sys_Totalneg_indR[i_trk][i_cent][i_jet]->SetTitle(Form("DeltaDpT SysNeg - pTjet: %1.0f-%1.0f GeV, pTtrk: %1.2f-%1.2f GeV, centrality: %s", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1), trkpT_binning->GetBinLowEdge(i_trk+1), trkpT_binning->GetBinUpEdge(i_trk+1), num_to_cent(31, i_cent).c_str()));

					h_DeltaDpT_final_indR[i_trk][i_cent][i_jet]->SetName(Form("h_DeltaDpT_trk%i_cent%i_jet%i", i_trk, i_cent, i_jet));
					h_DeltaDpT_final_sys_Totalpos_indR[i_trk][i_cent][i_jet]->SetName(Form("h_DeltaDpT_trk%i_cent%i_jet%i_sysP", i_trk, i_cent, i_jet));
					h_DeltaDpT_final_sys_Totalneg_indR[i_trk][i_cent][i_jet]->SetName(Form("h_DeltaDpT_trk%i_cent%i_jet%i_sysN", i_trk, i_cent, i_jet));

					h_DeltaDpT_final_indR[i_trk][i_cent][i_jet]->Write(Form("h_DeltaDpT_trk%i_cent%i_jet%i", i_trk, i_cent, i_jet));
					h_DeltaDpT_final_sys_Totalpos_indR[i_trk][i_cent][i_jet]->Write(Form("h_DeltaDpT_trk%i_cent%i_jet%i_sysP", i_trk, i_cent, i_jet));
					h_DeltaDpT_final_sys_Totalneg_indR[i_trk][i_cent][i_jet]->Write(Form("h_DeltaDpT_trk%i_cent%i_jet%i_sysN", i_trk, i_cent, i_jet));



					h_ChPS_PbPb_final_ratio_indR[i_trk][i_cent][i_jet]->SetTitle(Form("DpT PbPb Nominal - pTjet: %1.0f-%1.0f GeV, pTtrk: %1.2f-%1.2f GeV, centrality: %s", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1), trkpT_binning->GetBinLowEdge(i_trk+1), trkpT_binning->GetBinUpEdge(i_trk+1), num_to_cent(31, i_cent).c_str()));
					h_ChPS_PbPb_final_sys_Totalpos_indR[i_trk][i_cent][i_jet]->SetTitle(Form("DpT PbPb SysPos - pTjet: %1.0f-%1.0f GeV, pTtrk: %1.2f-%1.2f GeV, centrality: %s", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1), trkpT_binning->GetBinLowEdge(i_trk+1), trkpT_binning->GetBinUpEdge(i_trk+1), num_to_cent(31, i_cent).c_str()));
					h_ChPS_PbPb_final_sys_Totalneg_indR[i_trk][i_cent][i_jet]->SetTitle(Form("DpT PbPb SysNeg - pTjet: %1.0f-%1.0f GeV, pTtrk: %1.2f-%1.2f GeV, centrality: %s", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1), trkpT_binning->GetBinLowEdge(i_trk+1), trkpT_binning->GetBinUpEdge(i_trk+1), num_to_cent(31, i_cent).c_str()));

					h_ChPS_PbPb_final_ratio_indR[i_trk][i_cent][i_jet]->SetName(Form("h_DpT_PbPb_trk%i_cent%i_jet%i", i_trk, i_cent, i_jet));
					h_ChPS_PbPb_final_sys_Totalpos_indR[i_trk][i_cent][i_jet]->SetName(Form("h_DpT_PbPb_trk%i_cent%i_jet%i_sysP", i_trk, i_cent, i_jet));
					h_ChPS_PbPb_final_sys_Totalneg_indR[i_trk][i_cent][i_jet]->SetName(Form("h_DpT_PbPb_trk%i_cent%i_jet%i_sysN", i_trk, i_cent, i_jet));

					h_ChPS_PbPb_final_ratio_indR[i_trk][i_cent][i_jet]->Write(Form("h_DpT_PbPb_trk%i_cent%i_jet%i", i_trk, i_cent, i_jet));
					h_ChPS_PbPb_final_sys_Totalpos_indR[i_trk][i_cent][i_jet]->Write(Form("h_DpT_PbPb_trk%i_cent%i_jet%i_sysP", i_trk, i_cent, i_jet));
					h_ChPS_PbPb_final_sys_Totalneg_indR[i_trk][i_cent][i_jet]->Write(Form("h_DpT_PbPb_trk%i_cent%i_jet%i_sysN", i_trk, i_cent, i_jet));



					if (i_cent == 0)
					{
						h_ChPS_pp_final_ratio_indR[i_trk][6][i_jet]->SetTitle(Form("DpT pp Nominal - pTjet: %1.0f-%1.0f GeV, pTtrk: %1.2f-%1.2f GeV", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1), trkpT_binning->GetBinLowEdge(i_trk+1), trkpT_binning->GetBinUpEdge(i_trk+1)));
						h_ChPS_pp_final_sys_Totalpos_indR[i_trk][6][i_jet]->SetTitle(Form("DpT pp SysPos -  pTjet: %1.0f-%1.0f GeV, pTtrk: %1.2f-%1.2f GeV", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1), trkpT_binning->GetBinLowEdge(i_trk+1), trkpT_binning->GetBinUpEdge(i_trk+1)));
						h_ChPS_pp_final_sys_Totalneg_indR[i_trk][6][i_jet]->SetTitle(Form("DpT pp SysNeg - pTjet: %1.0f-%1.0f GeV, pTtrk: %1.2f-%1.2f GeV", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1), trkpT_binning->GetBinLowEdge(i_trk+1), trkpT_binning->GetBinUpEdge(i_trk+1)));

						h_ChPS_pp_final_ratio_indR[i_trk][6][i_jet]->SetName(Form("h_DpT_pp_trk%i_cent%i_jet%i", i_trk, 6, i_jet));
						h_ChPS_pp_final_sys_Totalpos_indR[i_trk][6][i_jet]->SetName(Form("h_DpT_pp_trk%i_cent%i_jet%i_sysP", i_trk, 6, i_jet));
						h_ChPS_pp_final_sys_Totalneg_indR[i_trk][6][i_jet]->SetName(Form("h_DpT_pp_trk%i_cent%i_jet%i_sysN", i_trk, 6, i_jet));

						h_ChPS_pp_final_ratio_indR[i_trk][6][i_jet]->Write(Form("h_DpT_pp_trk%i_cent%i_jet%i", i_trk, 6, i_jet));
						h_ChPS_pp_final_sys_Totalpos_indR[i_trk][6][i_jet]->Write(Form("h_DpT_pp_trk%i_cent%i_jet%i_sysP", i_trk, 6, i_jet));
						h_ChPS_pp_final_sys_Totalneg_indR[i_trk][6][i_jet]->Write(Form("h_DpT_pp_trk%i_cent%i_jet%i_sysN", i_trk, 6, i_jet));
					}
				}


				trk_itr++;
			}
		}
	}


	//getting rdpt with uncertainties as a function of track pt
	for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
	{
		for (int i_cent = 0; i_cent < 6; i_cent++)
		{
			int dR_itr = 0;
			for (int i_dR = 0; i_dR < N_dR; i_dR++)
			{
				if (i_dR > 10) continue;

				//RDpT
				name = Form("h_RDpT_final_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet);
				h_RDpT_final_ratio_inTrk[i_dR][i_cent][i_jet] = (TH1*)f_RDpT->Get(name.c_str());

				name = Form("h_RDpT_sys_dR%i_cent%i_jetpt%i_total_p", i_dR, i_cent, i_jet);
				h_RDpT_final_sys_Totalpos_inTrk[i_dR][i_cent][i_jet] = (TH1*)f_RDpT_sys->Get(name.c_str());

				name = Form("h_RDpT_sys_dR%i_cent%i_jetpt%i_total_n", i_dR, i_cent, i_jet);
				h_RDpT_final_sys_Totalneg_inTrk[i_dR][i_cent][i_jet] = (TH1*)f_RDpT_sys->Get(name.c_str());

				g_RDpT_final_sys_inTrk[i_dR][i_cent][i_jet] = new TGraphAsymmErrors(h_RDpT_final_ratio_inTrk[i_dR][i_cent][i_jet]);
				g_RDpT_final_sys_inTrk[i_dR][i_cent][i_jet]->GetYaxis()->SetTitle(rdptr_label.c_str());
				g_RDpT_final_sys_inTrk[i_dR][i_cent][i_jet]->GetXaxis()->SetTitle(trk_label.c_str());

				g_RDpT_final_stat_inTrk[i_dR][i_cent][i_jet] = new TGraphAsymmErrors(h_RDpT_final_ratio_inTrk[i_dR][i_cent][i_jet]);
				g_RDpT_final_stat_inTrk[i_dR][i_cent][i_jet]->GetYaxis()->SetTitle(rdptr_label.c_str());
				g_RDpT_final_stat_inTrk[i_dR][i_cent][i_jet]->GetXaxis()->SetTitle(trk_label.c_str());

				double nom, sys_hi, sys_lo, r_position, r_width, stat_hi, stat_lo;
				for (int i_trk = 0 ; i_trk < N_trkpt; i_trk++)
				{
					if (i_trk < 2 || i_trk > 9) continue;

					//RDpT
					r_position = h_RDpT_final_ratio_inTrk[i_dR][i_cent][i_jet]->GetBinCenter(i_trk+1);
					r_width = h_RDpT_final_ratio_inTrk[i_dR][i_cent][i_jet]->GetBinWidth(i_trk+1)*0.30;
					if (i_trk > 4 && i_trk < 7) r_width = h_RDpT_final_ratio_inTrk[i_dR][i_cent][i_jet]->GetBinWidth(i_trk+1)*0.40;
					if (i_trk >=7) r_width = h_RDpT_final_ratio_inTrk[i_dR][i_cent][i_jet]->GetBinWidth(i_trk+1)*0.20;

					nom = h_RDpT_final_ratio_inTrk[i_dR][i_cent][i_jet]->GetBinContent(i_trk+1);
					sys_hi = nom * (h_RDpT_final_sys_Totalpos_inTrk[i_dR][i_cent][i_jet]->GetBinContent(i_trk+1));
					sys_lo = fabs(nom * (h_RDpT_final_sys_Totalneg_inTrk[i_dR][i_cent][i_jet]->GetBinContent(i_trk+1)));
					stat_hi = h_RDpT_final_ratio_inTrk[i_dR][i_cent][i_jet]->GetBinError(i_trk+1);
					stat_lo = h_RDpT_final_ratio_inTrk[i_dR][i_cent][i_jet]->GetBinError(i_trk+1);

					if ((i_jet == 7 && i_trk == 6 && i_dR >= 6)
						||
						(i_jet == 7 && i_trk == 7 && i_dR >= 6)
						||
						(i_jet == 7 && i_trk == 8 && i_dR >= 4)
						||
						(i_jet == 8 && i_trk == 7 && i_dR >= 7)
						||
						(i_jet == 8 && i_trk == 8 && i_dR >= 6)
						||
						(i_jet == 9 && i_trk == 8 && i_dR >= 6)
						)
					{
						nom = -1;
					}
					g_RDpT_final_stat_inTrk[i_dR][i_cent][i_jet]->SetPoint(i_trk, r_position, nom );
					g_RDpT_final_stat_inTrk[i_dR][i_cent][i_jet]->SetPointError(i_trk, 0, 0, stat_lo, stat_hi);

					g_RDpT_final_sys_inTrk[i_dR][i_cent][i_jet]->SetPoint(i_trk, r_position, nom );
					g_RDpT_final_sys_inTrk[i_dR][i_cent][i_jet]->SetPointError(i_trk, r_width/2, r_width/2, sys_lo, sys_hi);

				}

				dR_itr++;
			}
		}
	}




	double *X = ((TArrayD*)jetpT_binning->GetXbins())->GetArray();
	//getting rdpt with uncertainties as a function of jet pt
	for (int i_trk = 0 ; i_trk < N_trkpt; i_trk++)
	{
		if (i_trk < 2 || i_trk > 9) continue;

		for (int i_cent = 0; i_cent < 6; i_cent++)
		{
			int dR_itr = 0;
			for (int i_dR = 0; i_dR < N_dR; i_dR++)
			{
				if (i_dR > 10) continue;

				//RDpT
				name = Form("h_RDpT_final_dR%i_cent%i_trk%i", i_dR, i_cent, i_trk);
				h_RDpT_final_ratio_inJet[i_dR][i_cent][i_trk] = new TH1D(name.c_str(), name.c_str(), jetpT_binning->GetNbins(), X);

				name = Form("h_RDpT_sys_dR%i_cent%i_trk%i_total_p", i_dR, i_cent, i_trk);
				h_RDpT_final_sys_Totalpos_inJet[i_dR][i_cent][i_trk] = new TH1D(name.c_str(), name.c_str(), jetpT_binning->GetNbins(), X);

				name = Form("h_RDpT_sys_dR%i_cent%i_trk%i_total_n", i_dR, i_cent, i_trk);
				h_RDpT_final_sys_Totalneg_inJet[i_dR][i_cent][i_trk] = new TH1D(name.c_str(), name.c_str(), jetpT_binning->GetNbins(), X);

				for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
				{
					h_RDpT_final_ratio_inJet[i_dR][i_cent][i_trk]->SetBinContent(i_jet+1, h_RDpT_final_ratio_inTrk[i_dR][i_cent][i_jet]->GetBinContent(i_trk+1));
					h_RDpT_final_ratio_inJet[i_dR][i_cent][i_trk]->SetBinError(i_jet+1, h_RDpT_final_ratio_inTrk[i_dR][i_cent][i_jet]->GetBinError(i_trk+1));

					h_RDpT_final_sys_Totalpos_inJet[i_dR][i_cent][i_trk]->SetBinContent(i_jet+1, h_RDpT_final_sys_Totalpos_inTrk[i_dR][i_cent][i_jet]->GetBinContent(i_trk+1));
					h_RDpT_final_sys_Totalpos_inJet[i_dR][i_cent][i_trk]->SetBinError(i_jet+1, h_RDpT_final_sys_Totalpos_inTrk[i_dR][i_cent][i_jet]->GetBinError(i_trk+1));

					h_RDpT_final_sys_Totalneg_inJet[i_dR][i_cent][i_trk]->SetBinContent(i_jet+1, h_RDpT_final_sys_Totalneg_inTrk[i_dR][i_cent][i_jet]->GetBinContent(i_trk+1));
					h_RDpT_final_sys_Totalneg_inJet[i_dR][i_cent][i_trk]->SetBinError(i_jet+1, h_RDpT_final_sys_Totalneg_inTrk[i_dR][i_cent][i_jet]->GetBinError(i_trk+1));
				}

				g_RDpT_final_sys_inJet[i_dR][i_cent][i_trk] = new TGraphAsymmErrors(h_RDpT_final_ratio_inJet[i_dR][i_cent][i_trk]);
				g_RDpT_final_sys_inJet[i_dR][i_cent][i_trk]->GetYaxis()->SetTitle(rdptr_label.c_str());
				g_RDpT_final_sys_inJet[i_dR][i_cent][i_trk]->GetXaxis()->SetTitle(jet_label.c_str());

				g_RDpT_final_stat_inJet[i_dR][i_cent][i_trk] = new TGraphAsymmErrors(h_RDpT_final_ratio_inJet[i_dR][i_cent][i_trk]);
				g_RDpT_final_stat_inJet[i_dR][i_cent][i_trk]->GetYaxis()->SetTitle(rdptr_label.c_str());
				g_RDpT_final_stat_inJet[i_dR][i_cent][i_trk]->GetXaxis()->SetTitle(jet_label.c_str());

				double nom, sys_hi, sys_lo, r_position, r_width, stat_hi, stat_lo;
				for (int i_jet = 0; i_jet < N_jetpt; i_jet++)
				{
					if (i_jet < 7 || i_jet > 10) continue;

					//RDpT
					r_position = h_RDpT_final_ratio_inJet[i_dR][i_cent][i_trk]->GetBinCenter(i_jet+1);
					r_width = h_RDpT_final_ratio_inJet[i_dR][i_cent][i_trk]->GetBinWidth(i_jet+1)*0.30;

					nom = h_RDpT_final_ratio_inJet[i_dR][i_cent][i_trk]->GetBinContent(i_jet+1);
					sys_hi = nom * (h_RDpT_final_sys_Totalpos_inJet[i_dR][i_cent][i_trk]->GetBinContent(i_jet+1));
					sys_lo = fabs(nom * (h_RDpT_final_sys_Totalneg_inJet[i_dR][i_cent][i_trk]->GetBinContent(i_jet+1)));
					stat_hi = h_RDpT_final_ratio_inJet[i_dR][i_cent][i_trk]->GetBinError(i_jet+1);
					stat_lo = h_RDpT_final_ratio_inJet[i_dR][i_cent][i_trk]->GetBinError(i_jet+1);

					if ((i_jet == 7 && i_trk == 6 && i_dR >= 6)
						||
						(i_jet == 7 && i_trk == 7 && i_dR >= 6)
						||
						(i_jet == 7 && i_trk == 8 && i_dR >= 4)
						||
						(i_jet == 8 && i_trk == 7 && i_dR >= 7)
						||
						(i_jet == 8 && i_trk == 8 && i_dR >= 6)
						||
						(i_jet == 9 && i_trk == 8 && i_dR >= 6)
						)
					{
						nom = -1;
					}

					g_RDpT_final_stat_inJet[i_dR][i_cent][i_trk]->SetPoint(i_jet, r_position, nom );
					g_RDpT_final_stat_inJet[i_dR][i_cent][i_trk]->SetPointError(i_jet, 0, 0, stat_lo, stat_hi);

					g_RDpT_final_sys_inJet[i_dR][i_cent][i_trk]->SetPoint(i_jet, r_position, nom );
					g_RDpT_final_sys_inJet[i_dR][i_cent][i_trk]->SetPointError(i_jet, r_width/2, r_width/2, sys_lo, sys_hi);

				}

				dR_itr++;
			}
		}
	}



	double y_range_lo = 0;
	double y_range_hi = 4;



	// DRAWING BEGINS HERRE

/*
	{
		cout << "Doing Final RDpT as a function of jet pT for different dR" << endl;

		TCanvas *canvas = new TCanvas("canvas","canvas",800,600);
		TLegend *legend = new TLegend(0.19,0.63,0.35,0.80,NULL,"brNDC");
		legend->SetTextFont(43);
		legend->SetTextSize(23);
		legend->SetBorderSize(0);
		legend->SetFillStyle(0);

		bool first_pass_cent = true;
		for (int i_cent = 0; i_cent < 6; i_cent++)
		{
			string centrality = num_to_cent(31,i_cent);

			int trk_itr = 0;
			for (int i_trk = 2; i_trk < 9; i_trk++)
			{
				string trk_label = Form("%1.1f < #it{p}_{T}^{trk} < %1.1f GeV", trkpT_binning->GetBinLowEdge(i_trk+1), trkpT_binning->GetBinUpEdge(i_trk+1));

				canvas->cd();
				canvas->Clear();

				int dR_itr = 0;

				for (int i_dR = 0; i_dR < N_dR; i_dR++)
				{
					string dR_label = Form("%1.2f < #it{r} < %1.2f ", dR_binning->GetBinLowEdge(i_dR+1), dR_binning->GetBinUpEdge(i_dR+1));
					//					if (i_dR > 9) continue;
					if (i_dR != 0 && i_dR != 3 && i_dR != 6 && i_dR != 8) continue;

					SetHStyle(h_RDpT_final_ratio_inJet[i_dR][i_cent][i_trk], dR_itr);
					SetHStyle_graph(g_RDpT_final_sys_inJet[i_dR][i_cent][i_trk], dR_itr);
					SetHStyle_graph(g_RDpT_final_stat_inJet[i_dR][i_cent][i_trk], dR_itr);
					g_RDpT_final_sys_inJet[i_dR][i_cent][i_trk]->SetFillColorAlpha(g_RDpT_final_sys_inJet[i_dR][i_cent][i_trk]->GetFillColor(),opacity);

					if (first_pass_cent && trk_itr == 0) legend->AddEntry(g_RDpT_final_stat_inJet[i_dR][i_cent][i_trk],dR_label.c_str(),"p");

					g_RDpT_final_sys_inJet[i_dR][i_cent][i_trk]->SetLineWidth(1.);
					g_RDpT_final_stat_inJet[i_dR][i_cent][i_trk]->SetLineWidth(1.);

					y_range_lo = 0.;
					y_range_hi = 4.2;

					TGraphAsymmErrors* g_sys = (TGraphAsymmErrors*)g_RDpT_final_sys_inJet[i_dR][i_cent][i_trk]->Clone("g_sys");
					TGraphAsymmErrors* g_stat = (TGraphAsymmErrors*)g_RDpT_final_stat_inJet[i_dR][i_cent][i_trk]->Clone("g_stat");

					double shift_size = 2;
					g_sys = shift(g_sys, dR_itr, shift_size);
					g_stat = shift(g_stat, dR_itr, shift_size);

					g_sys->GetXaxis()->SetLimits(126, 310);
					g_sys->GetYaxis()->SetRangeUser(y_range_lo,y_range_hi);
					g_sys->GetYaxis()->SetNdivisions(505);
					g_sys->GetXaxis()->SetMoreLogLabels(1);

					if (dR_itr == 0) g_sys->Draw("a E2");
					else g_sys->Draw("E2 same");
					g_stat->Draw("P E1");
					gPad->SetLogx(0);

					dR_itr++;
				} // end dR loop

				trk_itr++;

				canvas->cd();
				line->DrawLine(jet_pt_lo, 1, jet_pt_hi, 1);
				legend->Draw();

				double x_left = 0.19, x_right = 0.91, y = 0.88, y_diff = 0.045;
				ltx->SetTextAlign(11);
				ltx->DrawLatexNDC(x_left, y, trk_label.c_str());
				ltx->DrawLatexNDC(x_left, y-y_diff, Form("%s", centrality.c_str()));
				ltx->SetTextAlign(31);
				ltx->DrawLatexNDC(x_right, y, "#scale[1.5]{#font[72]{ATLAS}}");
				ltx->DrawLatexNDC(x_right, y=y-y_diff, "Pb+Pb #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}");
				ltx->DrawLatexNDC(x_right, y=y-y_diff, "#it{pp} #sqrt{#font[12]{s}} = 5.02 TeV, 25 pb^{-1}");
				ltx->DrawLatexNDC(x_right, y=y-y_diff, Form("anti-#font[12]{k}_{#font[12]{t}} R=0.4"));

				if (!((i_trk == 2 || i_trk == 5) && (i_cent == 0 || i_cent == 5))) continue;
				canvas->Print(Form("output_pdf_nominal/conf/RDpT_jetpt_trk%i_cent%i.pdf", i_trk, i_cent));

			} //end jet loop

			first_pass_cent = false;
		} //end cent loop
	}

 {

		cout << "Doing Final RDpT as a function of track pT for different dR" << endl;

		TCanvas *canvas = new TCanvas("canvas","canvas",800,600);
		TLegend *legend = new TLegend(0.45,0.55,0.55,0.72,NULL,"brNDC");
		legend->SetTextFont(43);
		legend->SetTextSize(25);
		legend->SetBorderSize(0);
		legend->SetFillStyle(0);

		bool first_pass_cent = true;
		for (int i_cent = 0; i_cent < 6; i_cent++)
		{
			string centrality = num_to_cent(31,i_cent);

			int jet_itr = 0;
			for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
			{
				string jet_label = Form("%1.0f < #it{p}_{T}^{jet} < %1.0f GeV", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));

				canvas->cd();
				canvas->Clear();

				int dR_itr = 0;

				for (int i_dR = 0; i_dR < N_dR; i_dR++)
				{
					string dR_label = Form("%1.2f < #it{r} < %1.2f ", dR_binning->GetBinLowEdge(i_dR+1), dR_binning->GetBinUpEdge(i_dR+1));
					if (i_dR != 0 && i_dR != 3 && i_dR != 6 && i_dR != 8) continue;

					SetHStyle(h_RDpT_final_ratio_inTrk[i_dR][i_cent][i_jet], dR_itr);
					SetHStyle_graph(g_RDpT_final_sys_inTrk[i_dR][i_cent][i_jet], dR_itr);
					SetHStyle_graph(g_RDpT_final_stat_inTrk[i_dR][i_cent][i_jet], dR_itr);
					g_RDpT_final_sys_inTrk[i_dR][i_cent][i_jet]->SetFillColorAlpha(g_RDpT_final_sys_inTrk[i_dR][i_cent][i_jet]->GetFillColor(),opacity);

					if (first_pass_cent && jet_itr == 0) legend->AddEntry(g_RDpT_final_stat_inTrk[i_dR][i_cent][i_jet],dR_label.c_str(),"p");

					g_RDpT_final_sys_inTrk[i_dR][i_cent][i_jet]->SetLineWidth(1.);
					g_RDpT_final_stat_inTrk[i_dR][i_cent][i_jet]->SetLineWidth(1.);

					y_range_lo = 0.;
					y_range_hi = 3.4;

					TGraphAsymmErrors* g_sys = (TGraphAsymmErrors*)g_RDpT_final_sys_inTrk[i_dR][i_cent][i_jet]->Clone("g_sys");
					TGraphAsymmErrors* g_stat = (TGraphAsymmErrors*)g_RDpT_final_stat_inTrk[i_dR][i_cent][i_jet]->Clone("g_stat");

					double shift_size = 0.045;
					g_sys = shift(g_sys, trkpT_binning, dR_itr, shift_size);
					g_stat = shift(g_stat, trkpT_binning, dR_itr, shift_size);

					g_sys->GetXaxis()->SetLimits(trk_pt_lo, trk_pt_hi);
					g_sys->GetYaxis()->SetRangeUser(0,y_range_hi);
					g_sys->GetYaxis()->SetNdivisions(505);

					if (dR_itr == 0) g_sys->Draw("a E2");
					else g_sys->Draw("E2 same");
					g_stat->Draw("P E1");
					gPad->SetLogx();

					dR_itr++;
				} // end dR loop

				jet_itr++;

				canvas->cd();
				line->DrawLine(trk_pt_lo, 1, trk_pt_hi, 1);
				legend->Draw();

				double x_left = 0.19, x_right = 0.92, y = 0.88, y_diff = 0.058;
				ltx->SetTextAlign(11);
				ltx->SetTextSize(28);
				ltx->DrawLatexNDC(x_left, y, jet_label.c_str());
				ltx->DrawLatexNDC(x_left, y-y_diff, Form("%s", centrality.c_str()));
				ltx->SetTextAlign(31);
				ltx->DrawLatexNDC(x_right, y, "#scale[1.5]{#font[72]{ATLAS}}");
				ltx->DrawLatexNDC(x_right, y=y-y_diff, "Pb+Pb #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}");
				ltx->DrawLatexNDC(x_right, y=y-y_diff, "#it{pp} #sqrt{#font[12]{s}} = 5.02 TeV, 25 pb^{-1}");
				ltx->DrawLatexNDC(x_right, y=y-y_diff, Form("anti-#font[12]{k}_{#font[12]{t}} R=0.4"));

				if (!((i_jet == 7 || i_jet == 9) && (i_cent == 0 || i_cent == 3 || i_cent == 5))) continue;
				canvas->Print(Form("output_pdf_nominal/conf/RDpT_trkpt_jet%i_cent%i.pdf", i_jet, i_cent));

			} //end jet loop

			first_pass_cent = false;
		} //end cent loop
	}

 {
		cout << "Doing Final RDpT as a function of track pT for different centrality" << endl;

		TCanvas *canvas = new TCanvas("canvas","canvas",800,600);
		TLegend *legend = new TLegend(0.7,0.6,0.85,0.72,NULL,"brNDC");
		legend->SetTextFont(43);
		legend->SetTextSize(23);
		legend->SetBorderSize(0);
		legend->SetNColumns(1);
		legend->SetFillStyle(0);

		int dR_itr = 0;
		for (int i_dR = 0; i_dR < N_dR; i_dR++)
		{
			string dR_label = Form("%1.2f < #it{r} < %1.2f ", dR_binning->GetBinLowEdge(i_dR+1), dR_binning->GetBinUpEdge(i_dR+1));
			if (i_dR != 0 && i_dR != 3 && i_dR != 10) continue;
			y_range_lo = 0.; y_range_hi = 4;
			if (i_dR == 0) {y_range_lo = 0.7; y_range_hi = 1.8;}
			if (i_dR == 3) {y_range_lo = 0.0; y_range_hi = 2.52;}

			int jet_itr = 0;
			for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
			{
				string jet_label = Form("%1.0f < #it{p}_{T}^{jet} < %1.0f GeV", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));

				canvas->cd();
				canvas->Clear();


				int cent_itr = 0;
				for (int i_cent = 0; i_cent < 6; i_cent++)
				{
					string centrality = num_to_cent(31,i_cent);
					if (i_cent != 0 && i_cent != 3 && i_cent != 5 ) continue;

					SetHStyle(h_RDpT_final_ratio_inTrk[i_dR][i_cent][i_jet], cent_itr);
					SetHStyle_graph(g_RDpT_final_sys_inTrk[i_dR][i_cent][i_jet], cent_itr);
					SetHStyle_graph(g_RDpT_final_stat_inTrk[i_dR][i_cent][i_jet], cent_itr);
					g_RDpT_final_sys_inTrk[i_dR][i_cent][i_jet]->SetFillColorAlpha(g_RDpT_final_sys_inTrk[i_dR][i_cent][i_jet]->GetFillColor(),opacity);

					if (dR_itr == 0 && jet_itr == 0) legend->AddEntry(g_RDpT_final_stat_inTrk[i_dR][i_cent][i_jet],centrality.c_str(),"p");

					g_RDpT_final_sys_inTrk[i_dR][i_cent][i_jet]->SetLineWidth(1.);
					g_RDpT_final_stat_inTrk[i_dR][i_cent][i_jet]->SetLineWidth(1.);

					TGraphAsymmErrors* g_sys = (TGraphAsymmErrors*)g_RDpT_final_sys_inTrk[i_dR][i_cent][i_jet]->Clone("g_sys");
					TGraphAsymmErrors* g_stat = (TGraphAsymmErrors*)g_RDpT_final_stat_inTrk[i_dR][i_cent][i_jet]->Clone("g_stat");

					double shift_size = 0.1;
					g_sys = shift(g_sys, trkpT_binning, cent_itr, shift_size);
					g_stat = shift(g_stat, trkpT_binning, cent_itr, shift_size);

					g_sys->GetXaxis()->SetLimits(trk_pt_lo, trk_pt_hi);
					g_sys->GetYaxis()->SetRangeUser(y_range_lo, y_range_hi);
					g_sys->GetYaxis()->SetNdivisions(505);

					if (i_cent == 0) g_sys->Draw("a E2");
					else g_sys->Draw("E2 same");
					g_stat->Draw("P E1");

					gPad->SetLogx();

					cent_itr++;

				} // end cent loop

				jet_itr++;

				canvas->cd();
				line->DrawLine(trk_pt_lo, 1, trk_pt_hi, 1);
				legend->Draw();

				double x_left = 0.19, x_right = 0.91, y = 0.88, y_diff = 0.045;
				ltx->SetTextAlign(11);
				ltx->DrawLatexNDC(x_left, y, jet_label.c_str());
				ltx->DrawLatexNDC(x_left, y-y_diff, Form("%s", dR_label.c_str()));
				ltx->SetTextAlign(31);
				ltx->DrawLatexNDC(x_right, y, "#scale[1.5]{#font[72]{ATLAS}}");
				ltx->DrawLatexNDC(x_right, y=y-y_diff, "Pb+Pb #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}");
				ltx->DrawLatexNDC(x_right, y=y-y_diff, "#it{pp} #sqrt{#font[12]{s}} = 5.02 TeV, 25 pb^{-1}");
				ltx->DrawLatexNDC(x_right, y=y-y_diff, Form("anti-#font[12]{k}_{#font[12]{t}} R=0.4"));

				canvas->Print(Form("output_pdf_nominal/conf/RDpT_trkpt_jet%i_dR%i.pdf", i_jet, i_dR));

				jet_itr++;
			} //end jet loop

			dR_itr++;
		} //end cent loop
	}

 {
		cout << "Doing Final RDpT ratio (PbPb/pp) plots in dR" << endl;

		TCanvas *canvas = new TCanvas("canvas","canvas",800,600);
		TLegend *legend = new TLegend(0.17, 0.55, 0.760, 0.74, "","brNDC");
		legend->SetTextFont(43);
		legend->SetBorderSize(0);
		legend->SetFillStyle(0);
		legend->SetTextSize(24);
		legend->SetNColumns(2);
		cout << legend->GetMargin() << endl;
		legend->SetMargin(0.15);
		int jet_itr = 0;
		for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
		{
			string jet_label = Form("%1.0f < #it{p}_{T}^{jet} < %1.0f GeV", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));

			bool first_pass_cent = true;
			for (int i_cent = 0; i_cent < 6; i_cent++)
			{

				string centrality = num_to_cent(31,i_cent);

				int trk_itr = 0;
				canvas->cd();
				canvas->Clear();

				for (int i_trk = trk_select1; i_trk < 9; i_trk++)
				{
					string trk_label = Form("%1.1f < #it{p}_{T} < %1.1f GeV", trkpT_binning->GetBinLowEdge(i_trk+1), trkpT_binning->GetBinUpEdge(i_trk+1));

					SetHStyle_smallify(h_RDpT_final_ratio_indR[i_trk][i_cent][i_jet], trk_itr, 0);
					SetHStyle_graph(g_RDpT_final_sys_indR[i_trk][i_cent][i_jet], trk_itr);
					SetHStyle_graph(g_RDpT_final_stat_indR[i_trk][i_cent][i_jet], trk_itr);

					if (jet_itr == 0 && first_pass_cent) legend->AddEntry(g_RDpT_final_sys_indR[i_trk][i_cent][i_jet],trk_label.c_str(),"p");

					y_range_lo = 0.;
					y_range_hi = 5;

					g_RDpT_final_sys_indR[i_trk][i_cent][i_jet]->SetLineColor(h_RDpT_final_ratio_indR[i_trk][i_cent][i_jet]->GetMarkerColor());
					g_RDpT_final_sys_indR[i_trk][i_cent][i_jet]->SetFillColorAlpha(g_RDpT_final_sys_indR[i_trk][i_cent][i_jet]->GetFillColor(),opacity);
					g_RDpT_final_sys_indR[i_trk][i_cent][i_jet]->SetLineWidth(1.);
					g_RDpT_final_stat_indR[i_trk][i_cent][i_jet]->SetLineWidth(1.);


					TGraphAsymmErrors* g_sys = (TGraphAsymmErrors*)g_RDpT_final_sys_indR[i_trk][i_cent][i_jet]->Clone("g_sys");
					TGraphAsymmErrors* g_stat = (TGraphAsymmErrors*)g_RDpT_final_stat_indR[i_trk][i_cent][i_jet]->Clone("g_stat");

					double shift_size = 0.0025;
					g_sys = shift(g_sys, trk_itr, shift_size);
					g_stat = shift(g_stat, trk_itr, shift_size);

					g_sys->GetXaxis()->SetLimits(0, 0.8);
					g_sys->GetYaxis()->SetRangeUser(0,5.2);
					g_sys->GetYaxis()->SetNdivisions(505);

					if (trk_itr == 0) g_sys->Draw("a E2 text");
					else g_sys->Draw("E2 same");
					g_stat->Draw("P E1");


					trk_itr++;

				} // end trk loop

				canvas->cd();
				line->DrawLine(0, 1, r_max_range, 1);
				ltx->SetTextSize(28);
				legend->Draw();

				double x_left = 0.19, x_right = 0.92, y = 0.88, y_diff = 0.058;
				ltx->SetTextAlign(31);
				ltx->DrawLatexNDC(x_right, y, Form("%s", jet_label.c_str()));
				ltx->DrawLatexNDC(x_right, y-y_diff, Form("%s", centrality.c_str()));
				ltx->SetTextAlign(11);
				ltx->DrawLatexNDC(x_left, y, "#scale[1.5]{#font[72]{ATLAS}}");
				ltx->DrawLatexNDC(x_left, y=y-y_diff, "Pb+Pb #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}");
				ltx->DrawLatexNDC(x_left, y=y-y_diff, "#it{pp} #sqrt{#font[12]{s}} = 5.02 TeV, 25 pb^{-1}");
				ltx->DrawLatexNDC(x_left+0.38, y, Form("anti-#font[12]{k}_{#font[12]{t}} R=0.4"));
				first_pass_cent = false;

				canvas->Print(Form("output_pdf_nominal/conf/RDpT_dR_jet%i_cent%i.pdf", i_jet, i_cent));

			} //end cent loop

			jet_itr++;
		} //end jet loop
	}


	{
		cout << "Doing Final DeltaDpT (PbPb - pp) plots in dR" << endl;

		TCanvas *canvas = new TCanvas("canvas","canvas",800,600);
		TLegend *legend = new TLegend(0.55, 0.40, 0.90, 0.80, "","brNDC");
		legend->SetTextFont(43);
		legend->SetBorderSize(0);
		legend->SetTextSize(27);
		legend->SetNColumns(1);
		legend->SetFillStyle(0);

		int jet_itr = 0;
		for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
		{
			string jet_label = Form("%1.0f < #it{p}_{T}^{jet} < %1.0f GeV", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));
			if (!(i_jet == 7 || i_jet == 8 || i_jet == 9 || i_jet == 10)) continue;

			bool first_pass_cent = true;
			for (int i_cent = 0; i_cent < 6; i_cent++)
			{
				if (i_cent != 0) continue;

				string centrality = num_to_cent(31,i_cent);

				int trk_itr = 0;
				canvas->cd();
				canvas->Clear();

				for (int i_trk = trk_select1; i_trk < 9; i_trk++)
				{
					string trk_label = Form("%1.1f < #it{p}_{T} < %1.1f GeV", trkpT_binning->GetBinLowEdge(i_trk+1), trkpT_binning->GetBinUpEdge(i_trk+1));

					SetHStyle_smallify(h_DeltaDpT_final_indR[i_trk][i_cent][i_jet], trk_itr, 0);
					SetHStyle_graph(g_DeltaDpT_final_sys_indR[i_trk][i_cent][i_jet], trk_itr);
					SetHStyle_graph(g_DeltaDpT_final_stat_indR[i_trk][i_cent][i_jet], trk_itr);

					if (jet_itr == 0 && first_pass_cent) legend->AddEntry(g_DeltaDpT_final_sys_indR[i_trk][i_cent][i_jet],trk_label.c_str(),"p");

					y_range_lo = -2.;
					y_range_hi = 16;

					g_DeltaDpT_final_sys_indR[i_trk][i_cent][i_jet]->SetLineColor(h_DeltaDpT_final_indR[i_trk][i_cent][i_jet]->GetMarkerColor());
					g_DeltaDpT_final_sys_indR[i_trk][i_cent][i_jet]->SetFillColorAlpha(g_DeltaDpT_final_sys_indR[i_trk][i_cent][i_jet]->GetFillColor(),opacity);
					g_DeltaDpT_final_sys_indR[i_trk][i_cent][i_jet]->SetLineWidth(1.);
					g_DeltaDpT_final_stat_indR[i_trk][i_cent][i_jet]->SetLineWidth(1.);

					TGraphAsymmErrors* g_sys = (TGraphAsymmErrors*)g_DeltaDpT_final_sys_indR[i_trk][i_cent][i_jet]->Clone("g_sys");
					TGraphAsymmErrors* g_stat = (TGraphAsymmErrors*)g_DeltaDpT_final_stat_indR[i_trk][i_cent][i_jet]->Clone("g_stat");

					double shift_size = 0.0040;
					g_sys = shift(g_sys, trk_itr, shift_size);
					g_stat = shift(g_stat, trk_itr, shift_size);

					g_sys->GetXaxis()->SetLimits(0, r_max_range);
					g_sys->GetYaxis()->SetRangeUser(y_range_lo, y_range_hi);

//					if (i_jet < 9) g_sys->GetYaxis()->SetRangeUser(y_range_lo, y_range_hi);
//					else if (i_jet >= 9) g_sys->GetYaxis()->SetRangeUser(y_range_lo, y_range_hi+7);
//					if (i_jet == 8) g_sys->GetYaxis()->SetRangeUser(y_range_lo-1., y_range_hi+2);
//					if (i_jet == 9) g_sys->GetYaxis()->SetRangeUser(y_range_lo-1, y_range_hi+5);
//					if (i_jet == 10) g_sys->GetYaxis()->SetRangeUser(y_range_lo-2., y_range_hi+11);
					g_sys->GetYaxis()->SetNdivisions(505);

					if (trk_itr == 0) g_sys->Draw("a E2");
					else g_sys->Draw("E2 same");
					g_stat->Draw("P E1");


					trk_itr++;

				} // end trk loop

				canvas->cd();
				line->DrawLine(0, 0, r_max_range, 0);
				legend->Draw();

				double x_left = 0.19, x_right = 0.92, y = 0.88, y_diff = 0.058;
				ltx->SetTextSize(28);
				ltx->SetTextAlign(31);
				ltx->DrawLatexNDC(x_right, y, jet_label.c_str());
				ltx->DrawLatexNDC(x_right, y-y_diff, Form("%s", centrality.c_str()));
				ltx->SetTextAlign(11);
				ltx->DrawLatexNDC(x_left, y, "#scale[1.5]{#font[72]{ATLAS}}");
				ltx->DrawLatexNDC(x_left, y=y-y_diff, "Pb+Pb #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}");
				ltx->DrawLatexNDC(x_left, y=y-y_diff, "#it{pp} #sqrt{#font[12]{s}} = 5.02 TeV, 25 pb^{-1}");
				ltx->DrawLatexNDC(x_left, y=y-y_diff, Form("anti-#font[12]{k}_{#font[12]{t}} R=0.4"));

				first_pass_cent = false;

				canvas->Print(Form("output_pdf_nominal/conf/DeltaDpT_dR_jet%i_cent%i.pdf", i_jet, i_cent));

			} //end cent loop

			jet_itr++;
		} //end jet loop
	}

	{
		trk_select1 = 3;
		cout << "Doing Final RDpT for different jets for one track pT" << endl;

		TCanvas *canvas = new TCanvas("canvas","canvas",800,600);

		TLegend *legend = new TLegend(0.6,0.55,0.90,0.82,NULL,"brNDC");
		legend->SetTextFont(43);
		legend->SetBorderSize(0);
		legend->SetTextSize(21);
		legend->SetNColumns(2);
		legend->SetFillStyle(0);

		bool first_pass_cent = true;
		for (int i_cent = 0; i_cent < 6; i_cent++)
		{
			string centrality = num_to_cent(31,i_cent);

			canvas->cd();
			canvas->Clear();


			int jet_itr = 0;
			for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
			{
				string jet_label = Form("%1.0f < #it{p}_{T}^{jet} < %1.0f GeV", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));

				int trk_itr = 0;
				for (int i_trk = 0; i_trk < N_trkpt; i_trk++)
				{
					string trk_label = Form("%1.1f < #it{p}_{T} < %1.1f GeV", trkpT_binning->GetBinLowEdge(i_trk+1), trkpT_binning->GetBinUpEdge(i_trk+1));

					if ((i_trk != trk_select1 && i_trk != trk_select2)) continue;

					if (i_trk == trk_select1)
					{
						
						SetHStyle(h_RDpT_final_ratio_indR[i_trk][i_cent][i_jet], jet_itr);
						SetHStyle_graph(g_RDpT_final_sys_indR[i_trk][i_cent][i_jet], jet_itr);
						SetHStyle_graph(g_RDpT_final_stat_indR[i_trk][i_cent][i_jet], jet_itr);
						g_RDpT_final_sys_indR[i_trk][i_cent][i_jet]->SetFillColorAlpha(g_RDpT_final_sys_indR[i_trk][i_cent][i_jet]->GetFillColor(),opacity);

						if (first_pass_cent) legend->AddEntry(g_RDpT_final_stat_indR[i_trk][i_cent][i_jet],"   ","p");

					}
					if (i_trk == trk_select2)
					{
						SetHStyle_open(h_RDpT_final_ratio_indR[i_trk][i_cent][i_jet], jet_itr);
						SetHStyle_graph(g_RDpT_final_sys_indR[i_trk][i_cent][i_jet], jet_itr);
						SetHStyle_graph_open(g_RDpT_final_stat_indR[i_trk][i_cent][i_jet], jet_itr);
						g_RDpT_final_sys_indR[i_trk][i_cent][i_jet]->SetFillColorAlpha(g_RDpT_final_sys_indR[i_trk][i_cent][i_jet]->GetFillColor(),opacity);
						if (first_pass_cent) legend->AddEntry(g_RDpT_final_stat_indR[i_trk][i_cent][i_jet],jet_label.c_str(),"p");

					}

					g_RDpT_final_sys_indR[i_trk][i_cent][i_jet]->SetLineWidth(1.);
					g_RDpT_final_stat_indR[i_trk][i_cent][i_jet]->SetLineWidth(1.);

					y_range_lo = 0.;
					y_range_hi = 6;

					TGraphAsymmErrors* g_sys = (TGraphAsymmErrors*)g_RDpT_final_sys_indR[i_trk][i_cent][i_jet]->Clone("g_sys");
					TGraphAsymmErrors* g_stat = (TGraphAsymmErrors*)g_RDpT_final_stat_indR[i_trk][i_cent][i_jet]->Clone("g_stat");

					double shift_size = 0.008;
					g_sys = shift(g_sys, jet_itr, shift_size);
					g_stat = shift(g_stat, jet_itr, shift_size);

					g_sys->GetXaxis()->SetLimits(0, r_max_range);
					g_sys->GetYaxis()->SetRangeUser(y_range_lo, y_range_hi);
					g_sys->GetYaxis()->SetNdivisions(505);

					if (trk_itr == 0 && jet_itr == 0) g_sys->Draw("a E2");
					else g_sys->Draw("E2 same");
					g_stat->Draw("P E1");

					trk_itr++;


				} // end jet loop


				jet_itr++;

			} //end cent loop

			ltx->SetTextAlign(21);
			double tmp_height = 0.87;
			double left_diff = 0;

			canvas->cd();
			line->DrawLine(0, 1, r_max_range, 1);
			legend->SetHeader("#it{p}_{T,1}   #it{p}_{T,2}");
			legend->Draw();

			double x_left = 0.19, x_right = 0.92, y = 0.88, y_diff = 0.058;
			ltx->SetTextSize(28);
			ltx->SetTextAlign(31);
			ltx->DrawLatexNDC(x_right, y=y, Form("%s", centrality.c_str()));
			ltx->SetTextAlign(11);
			ltx->DrawLatexNDC(x_left, y, "#scale[1.5]{#font[72]{ATLAS}}");
			ltx->DrawLatexNDC(x_left, y=y-y_diff, "Pb+Pb #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}");
			ltx->DrawLatexNDC(x_left, y=y-y_diff, "#it{pp} #sqrt{#font[12]{s}} = 5.02 TeV, 25 pb^{-1}");
			ltx->DrawLatexNDC(x_left, y=y-y_diff, Form("anti-#font[12]{k}_{#font[12]{t}} R=0.4"));
			ltx->DrawLatexNDC(x_left, y=y-y_diff, Form("%1.1f < #it{p}_{T,1} < %1.1f GeV",trkpT_binning->GetBinLowEdge(trk_select1+1), trkpT_binning->GetBinUpEdge(trk_select1+1)));
			ltx->DrawLatexNDC(x_left, y=y-y_diff, Form("%1.1f < #it{p}_{T,2} < %1.1f GeV",trkpT_binning->GetBinLowEdge(trk_select2+1), trkpT_binning->GetBinUpEdge(trk_select2+1)));

//			if (i_cent != 0 && i_cent != 5) continue;

			canvas->Print(Form("output_pdf_nominal/conf/RDpT_dR_trk%i_trk%i_cent%i.pdf", trk_select1, trk_select2, i_cent));

			first_pass_cent = false;
		} //end jet loop
		trk_select1 = 2;
	}
*/
	{
		cout << "Doing Final DpT spectra for PbPb and pp in dR" << endl;

		TCanvas *canvas = new TCanvas("canvas","canvas",800,900);
		double y_diff = 0.6, x_diff=0.24;
		TLegend *legend = new TLegend(0.20+x_diff, 0.20, 0.75+x_diff, 0.47,NULL,"brNDC");
		legend->SetTextFont(43);
		legend->SetBorderSize(0);
		legend->SetTextSize(17);
		legend->SetNColumns(2);
		legend->SetFillStyle(0);
		legend->SetHeader("#it{pp}   PbPb");

		TLegend *legend_open = new TLegend(0.10+x_diff, 0.26, 0.65+x_diff, 0.47, NULL,"brNDC");
		legend_open->SetTextFont(43);
		legend_open->SetBorderSize(0);
		legend_open->SetTextSize(17);
		legend_open->SetNColumns(2);
		legend_open->SetFillStyle(0);
		legend_open->SetHeader("#it{pp}   PbPb");

		int jet_itr = 0;
		for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
		{
			string jet_label = Form("%1.0f < #it{p}_{T}^{jet} < %1.0f GeV", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));

			bool first_pass_cent = true;
			canvas->cd(0);
			gPad->Clear();
			canvas->Clear();
			canvas->Divide(2,3, 0., 0., 0);
			for (int i_cent = 0; i_cent < 6; i_cent++)
			{
				string centrality = num_to_cent(31,i_cent);

				int trk_itr = 0;
				canvas->cd(i_cent+1);

				double x1, x2, y1, y2;
				double x1_m, x2_m, y1_m, y2_m;
				gPad->GetPadPar(x1, y1, x2, y2);
				x1_m = gPad->GetLeftMargin();
				x2_m = gPad->GetRightMargin();
				y1_m = gPad->GetBottomMargin();
				y2_m = gPad->GetTopMargin();

				if (i_cent%2 == 1)
				{
					x1 = 0.55;
					x1_m = 0.02;
					x2 = 1;
					x2_m = 0.04;

				}

				if (i_cent%2 == 0)
				{
					x2 = 0.54;
					x2_m = 0.05;
				}

				y1+=0.01;
				gPad->SetMargin(x1_m, x2_m, y1_m, y2_m);
				gPad->SetPad(x1, y1, x2, y2);


				for (int i_trk = trk_select1; i_trk < 9; i_trk++)
				{
					string trk_label = Form("%1.1f < #it{p}_{T} < %1.1f GeV", trkpT_binning->GetBinLowEdge(i_trk+1), trkpT_binning->GetBinUpEdge(i_trk+1));
					if (i_trk == 6) {
						trk_label = Form("  %1.1f < #it{p}_{T} < %1.1f GeV", trkpT_binning->GetBinLowEdge(i_trk+1), trkpT_binning->GetBinUpEdge(i_trk+1));
					}

					SetHStyle_smallify(h_ChPS_PbPb_final_ratio_indR[i_trk][i_cent][i_jet], trk_itr, 1);
					SetHStyle_graph_smallify(g_ChPS_PbPb_final_sys_indR[i_trk][i_cent][i_jet], trk_itr, 1);
					SetHStyle_graph_smallify(g_ChPS_PbPb_final_stat_indR[i_trk][i_cent][i_jet], trk_itr, 1);
					SetHStyle_graph_smallify(g_ChPS_pp_final_sys_indR[i_trk][6][i_jet], trk_itr, 1);

					SetHStyle_open_smallify(h_ChPS_pp_final_ratio_indR[i_trk][6][i_jet], trk_itr, 1);
					SetHStyle_graph_open(g_ChPS_pp_final_stat_indR[i_trk][6][i_jet], trk_itr);
					g_ChPS_pp_final_stat_indR[i_trk][6][i_jet]->SetMarkerSize(h_ChPS_pp_final_ratio_indR[i_trk][6][i_jet]->GetMarkerSize());
					g_ChPS_PbPb_final_stat_indR[i_trk][i_cent][i_jet]->SetLineWidth(1.);
					g_ChPS_pp_final_stat_indR[i_trk][6][i_jet]->SetLineWidth(1.);

					if (jet_itr == 0 && first_pass_cent)
					{
						if (trk_itr < 4)
						{
							legend->AddEntry(g_ChPS_pp_final_stat_indR[i_trk][6][i_jet]," ","p");
							legend->AddEntry(g_ChPS_PbPb_final_stat_indR[i_trk][i_cent][i_jet],trk_label.c_str(),"p");
						}
						else
						{
							legend_open->AddEntry(g_ChPS_pp_final_stat_indR[i_trk][6][i_jet]," ","p");
							legend_open->AddEntry(g_ChPS_PbPb_final_stat_indR[i_trk][i_cent][i_jet],trk_label.c_str(),"p");
						}

					}

					y_range_hi = 8E1;
					y_range_lo = 2E-5;

					if (i_cent > 1) {
						y_range_lo = 5E-4;
					}
					if (i_cent > 3) {
						y_range_lo = 1E-5;
					}
					g_ChPS_PbPb_final_sys_indR[i_trk][i_cent][i_jet]->SetFillColorAlpha(g_ChPS_PbPb_final_sys_indR[i_trk][i_cent][i_jet]->GetFillColor(),opacity);
					g_ChPS_pp_final_sys_indR[i_trk][6][i_jet]->SetFillColorAlpha(g_ChPS_pp_final_sys_indR[i_trk][6][i_jet]->GetFillColor(),opacity);

					TGraphAsymmErrors* g_sys_pbpb = (TGraphAsymmErrors*)g_ChPS_PbPb_final_sys_indR[i_trk][i_cent][i_jet]->Clone("g_sys_PbPb");
					TGraphAsymmErrors* g_stat_pbpb = (TGraphAsymmErrors*)g_ChPS_PbPb_final_stat_indR[i_trk][i_cent][i_jet]->Clone("g_stat_PbPb");
					TGraphAsymmErrors* g_sys_pp = (TGraphAsymmErrors*)g_ChPS_pp_final_sys_indR[i_trk][6][i_jet]->Clone("g_sys_pp");
					TGraphAsymmErrors* g_stat_pp = (TGraphAsymmErrors*)g_ChPS_pp_final_stat_indR[i_trk][6][i_jet]->Clone("g_stat_pp");

					double shift_size = 0.0025;
					g_sys_pbpb = shift(g_sys_pbpb, trk_itr, shift_size);
					g_stat_pbpb = shift(g_stat_pbpb, trk_itr, shift_size);

					g_sys_pp = shift(g_sys_pp, trk_itr, shift_size);
					g_stat_pp = shift(g_stat_pp, trk_itr, shift_size);


					g_sys_pbpb->GetXaxis()->SetRangeUser(0, r_max_range);
					g_sys_pbpb->GetYaxis()->SetRangeUser(y_range_lo, y_range_hi);
					g_sys_pbpb->GetYaxis()->SetNdivisions(505);
					g_sys_pbpb->GetYaxis()->SetTitleOffset(3.);
					if (i_cent%2 == 1) g_sys_pbpb->GetYaxis()->SetLabelOffset(4);

					g_sys_pbpb->GetYaxis()->SetTitleFont(43);
					g_sys_pbpb->GetYaxis()->SetTitleSize(23);
					g_sys_pbpb->GetYaxis()->SetLabelSize(19);

					g_sys_pbpb->GetXaxis()->SetTitleOffset(2.4);
					g_sys_pbpb->GetXaxis()->SetTitleFont(43);
					g_sys_pbpb->GetXaxis()->SetTitleSize(23);
					g_sys_pbpb->GetXaxis()->SetLabelSize(19);
					g_sys_pbpb->GetXaxis()->SetNdivisions(505);
					gPad->SetFillColor(0);
					if (trk_itr == 0) g_sys_pbpb->Draw("a E2");
					else g_sys_pbpb->Draw("E2 same");
					g_stat_pbpb->Draw("P E1");
					g_sys_pp->Draw("E2");
					g_stat_pp->Draw("P E1");


					gPad->SetLogy();

					trk_itr++;
				} // end trk loop

				first_pass_cent = false;
				double x_left = 0.19, x_right = 0.91, y = 0.88, y_diff = 0.045;
				ltx->SetTextAlign(31);
				ltx->SetTextSize(22);
				ltx->DrawLatexNDC(x_right, y, Form("%s", centrality.c_str()));

			} //end cent loop
			canvas->cd(5);
			legend->Draw();
			canvas->cd(6);
			legend_open->Draw();

			canvas->cd(1);
			double x_left = 0.22, x_right = 0.92, y = 0.87, y_diff = 0.1;
			ltx->SetTextAlign(11);
			ltx->SetTextSize(22);
			ltx->DrawLatexNDC(x_left, 0.07, "#scale[1.5]{#font[72]{ATLAS}}");
			ltx->SetTextAlign(31);
			ltx->DrawLatexNDC(x_right, y=0.17, Form("anti-#font[12]{k}_{#font[12]{t}} R=0.4"));
			ltx->DrawLatexNDC(x_right, y=y-y_diff, Form("%s", jet_label.c_str()));
			canvas->cd(2);
			ltx->DrawLatexNDC(x_right, y=0.17, "#it{pp} #sqrt{#font[12]{s}} = 5.02 TeV, 25 pb^{-1}");
			ltx->DrawLatexNDC(x_right, y=y-y_diff, "Pb+Pb #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}");


			canvas->Print(Form("output_pdf_nominal/conf/DpT_dR_jet%i.pdf", i_jet));

			jet_itr++;
		} //end jet loop
	}
/*
	{
		cout << "Doing Final RDpT ratio (PbPb/pp) plots in dR for all jet/cent" << endl;

		TCanvas *canvas = new TCanvas("canvas","canvas",800,900);
		TLegend *legend = new TLegend(0.2,0.21,0.85,0.42,NULL,"brNDC");
		legend->SetTextFont(43);
		legend->SetBorderSize(0);
		legend->SetTextSize(12);
		legend->SetFillStyle(0);
		legend->SetNColumns(2);

		int jet_itr = 0;
		for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
		{
			string jet_label = Form("%1.0f < #it{p}_{T}^{jet} < %1.0f GeV", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));

			bool first_pass_cent = true;
			canvas->cd();
			canvas->Clear();
			canvas->Divide(2,3, 0., 0., 0);

			for (int i_cent = 0; i_cent < 6; i_cent++)
			{

				string centrality = num_to_cent(31,i_cent);

				int trk_itr = 0;
				canvas->cd(i_cent+1);


				double x1, x2, y1, y2;
				double x1_m, x2_m, y1_m, y2_m;
				gPad->GetPadPar(x1, y1, x2, y2);
				x1_m = gPad->GetLeftMargin();
				x2_m = gPad->GetRightMargin();
				y1_m = gPad->GetBottomMargin();
				y2_m = gPad->GetTopMargin();

				if (i_cent%2 == 1)
				{
					x1 = 0.55;
					x1_m = 0.02;
					x2 = 1;
					x2_m = 0.04;

				}

				if (i_cent%2 == 0)
				{
					x2 = 0.54;
					x2_m = 0.05;
				}

				y1+=0.01;
				gPad->SetMargin(x1_m, x2_m, y1_m, y2_m);
				gPad->SetPad(x1, y1, x2, y2);

//
//				double xlo, ylo, xhi, yhi;
//				if (i_cent == 0) xlo = 0, xhi = 0.54, ylo = 0.68, yhi = 0.98;
//				if (i_cent == 1) xlo = 0.55, xhi = 1.00, ylo = 0.68, yhi = 0.98;
//				if (i_cent == 2) xlo = 0, xhi = 0.54, ylo-=0.3, yhi-=0.3;
//				if (i_cent == 3) xlo = 0.55, xhi = 1.00, ylo, yhi;
//				if (i_cent == 4) xlo = 0, xhi = 0.54, ylo-=0.37, yhi-=0.3;
//				if (i_cent == 5) xlo = 0.55, xhi = 1.00, ylo, yhi;
//				gPad->SetPad(xlo, ylo, xhi, yhi);
//
//				if (i_cent %2 == 1)gPad->SetLeftMargin(0.010);
//				else gPad->SetRightMargin(0.014);
//
//
//				if (i_cent == 0 || i_cent == 1) {gPad->SetBottomMargin(0.01);}
//				if (i_cent == 2 || i_cent == 3) {gPad->SetBottomMargin(0.01);}


				for (int i_trk = trk_select1; i_trk < 9; i_trk++)
				{
					string trk_label = Form("%1.1f < #it{p}_{T} < %1.1f GeV", trkpT_binning->GetBinLowEdge(i_trk+1), trkpT_binning->GetBinUpEdge(i_trk+1));
					if (i_trk<7) trk_label = Form("  %1.1f < #it{p}_{T} < %1.1f GeV", trkpT_binning->GetBinLowEdge(i_trk+1), trkpT_binning->GetBinUpEdge(i_trk+1));
					SetHStyle_smallify(h_RDpT_final_ratio_indR[i_trk][i_cent][i_jet], trk_itr, 1);
					SetHStyle_graph_smallify(g_RDpT_final_sys_indR[i_trk][i_cent][i_jet], trk_itr, 1);
					SetHStyle_graph_smallify(g_RDpT_final_stat_indR[i_trk][i_cent][i_jet], trk_itr, 1);

					if (jet_itr == 0 && first_pass_cent) legend->AddEntry(g_RDpT_final_sys_indR[i_trk][i_cent][i_jet],trk_label.c_str(),"p");

					y_range_lo = -0.2;
//					y_range_hi = 3.2;
					if (i_cent < 6) {y_range_hi = 2.6; y_range_lo = -0.3;}
					if (i_cent < 4) {y_range_hi = 3.3;}
					if (i_cent < 2) {y_range_hi = 3.1;}
					g_RDpT_final_sys_indR[i_trk][i_cent][i_jet]->SetLineColor(h_RDpT_final_ratio_indR[i_trk][i_cent][i_jet]->GetMarkerColor());
					g_RDpT_final_sys_indR[i_trk][i_cent][i_jet]->SetFillColorAlpha(g_RDpT_final_sys_indR[i_trk][i_cent][i_jet]->GetFillColor(),opacity);
					g_RDpT_final_sys_indR[i_trk][i_cent][i_jet]->SetLineWidth(1.);
					g_RDpT_final_stat_indR[i_trk][i_cent][i_jet]->SetLineWidth(1.);

					TGraphAsymmErrors* g_sys = (TGraphAsymmErrors*)g_RDpT_final_sys_indR[i_trk][i_cent][i_jet]->Clone("g_sys");
					TGraphAsymmErrors* g_stat = (TGraphAsymmErrors*)g_RDpT_final_stat_indR[i_trk][i_cent][i_jet]->Clone("g_stat");

					double shift_size = 0.0025;
					g_sys = shift(g_sys, trk_itr, shift_size);
					g_stat = shift(g_stat, trk_itr, shift_size);

					g_sys->GetXaxis()->SetLimits(0, r_max_range);
					g_sys->GetYaxis()->SetRangeUser(y_range_lo, y_range_hi);

					g_sys->GetYaxis()->SetTitleFont(43);
					g_sys->GetYaxis()->SetTitleSize(24);
					g_sys->GetYaxis()->SetLabelSize(19);
					g_sys->GetYaxis()->SetNdivisions(504);

					g_sys->GetXaxis()->SetTitleOffset(2.4);
					g_sys->GetXaxis()->SetTitleFont(43);
					g_sys->GetXaxis()->SetTitleSize(23);
					g_sys->GetXaxis()->SetLabelSize(19);
					g_sys->GetXaxis()->SetNdivisions(505);
					if (i_cent%2 == 1) g_sys->GetYaxis()->SetLabelOffset(4);

					if (trk_itr == 0) g_sys->Draw("a E2");
					else g_sys->Draw("E2 same");
					g_stat->Draw("P E1");


					trk_itr++;

				} // end trk loop
				first_pass_cent = false;
				double x_left = 0.19, x_right = 0.91, y = 0.88, y_diff = 0.045;
				ltx->SetTextAlign(11);
				ltx->SetTextSize(22);
				if (i_cent%2 == 1)x_left = 0.07;
				else x_left = 0.22;
				if (i_cent > 3) y+=0.02;

				ltx->DrawLatexNDC(x_left, y, Form("%s", centrality.c_str()));
				line->DrawLine(0, 1, r_max_range, 1);

			}
			canvas->cd(5);
			legend->Draw();

			canvas->cd(1);
			double x_left = 0.21, x_right = 0.92, y = 0.91, y_diff = 0.08;
			ltx->SetTextAlign(32);
			ltx->DrawLatexNDC(x_right, y, "#scale[1.5]{#font[72]{ATLAS}}");
			canvas->cd(6);

			ltx->SetTextAlign(32);
			ltx->DrawLatexNDC(x_right, y=y, Form("anti-#font[12]{k}_{#font[12]{t}} R=0.4"));
			ltx->DrawLatexNDC(x_right, y=y-y_diff, Form("%s", jet_label.c_str()));
			ltx->DrawLatexNDC(x_right, y=y-(6.2*y_diff), "#it{pp} #sqrt{#font[12]{s}} = 5.02 TeV, 25 pb^{-1}");
			ltx->DrawLatexNDC(x_right, y=y-y_diff, "Pb+Pb #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}");

			canvas->Print(Form("output_pdf_nominal/conf/RDpT_dR_jet%i.pdf", i_jet));

			jet_itr++;
		} //end jet loop
	}
/*
	{
		cout << "Doing Final DpT spectra for PbPb and pp in dR" << endl;

		TCanvas *canvas = new TCanvas("canvas","canvas",800,600);
		TLegend *legend = new TLegend(0.49,0.19,0.81,0.46,NULL,"brNDC");
		legend->SetTextFont(43);
		legend->SetBorderSize(0);
		legend->SetTextSize(23);
		legend->SetNColumns(2);
		legend->SetFillStyle(0);

		int jet_itr = 0;
		for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
		{
			string jet_label = Form("%1.0f < #it{p}_{T}^{jet} < %1.0f GeV", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));

			bool first_pass_cent = true;
			for (int i_cent = 0; i_cent < 6; i_cent++)
			{
				string centrality = num_to_cent(31,i_cent);

				int trk_itr = 0;
				canvas->cd();
				canvas->Clear();

				for (int i_trk = trk_select1; i_trk < 9; i_trk++)
				{
					string trk_label = Form("%1.1f < #it{p}_{T} < %1.1f GeV", trkpT_binning->GetBinLowEdge(i_trk+1), trkpT_binning->GetBinUpEdge(i_trk+1));

					SetHStyle(h_ChPS_PbPb_final_ratio_indR[i_trk][i_cent][i_jet], trk_itr);
					SetHStyle_graph(g_ChPS_PbPb_final_sys_indR[i_trk][i_cent][i_jet], trk_itr);
					SetHStyle_graph(g_ChPS_PbPb_final_stat_indR[i_trk][i_cent][i_jet], trk_itr);
					SetHStyle_graph(g_ChPS_pp_final_sys_indR[i_trk][6][i_jet], trk_itr);

					SetHStyle_open(h_ChPS_pp_final_ratio_indR[i_trk][6][i_jet], trk_itr);
					SetHStyle_graph_open(g_ChPS_pp_final_stat_indR[i_trk][6][i_jet], trk_itr);

 					g_ChPS_PbPb_final_stat_indR[i_trk][i_cent][i_jet]->SetLineWidth(1.);
 					g_ChPS_pp_final_stat_indR[i_trk][6][i_jet]->SetLineWidth(1.);

					if (jet_itr == 0 && first_pass_cent) legend->AddEntry(g_ChPS_pp_final_stat_indR[i_trk][6][i_jet]," ","p");
					if (jet_itr == 0 && first_pass_cent) legend->AddEntry(g_ChPS_PbPb_final_stat_indR[i_trk][i_cent][i_jet],trk_label.c_str(),"p");

					y_range_lo = 5E-5;
					y_range_hi = 5E2;

					if (i_jet != 7)
					{
						y_range_lo = 8E-8;
						y_range_hi = 5E3;
					}
					if (i_jet == 10)
					{
						y_range_lo = 8E-9;
						y_range_hi = 5E3;
					}

					if (i_jet == 10)
					{
						double x1_coord = legend->GetX1NDC(), x2_coord = legend->GetX2NDC(), y1_coord = legend->GetY1NDC(), y2_coord = legend->GetY2NDC();
						legend->SetX1NDC(0.2);
						legend->SetX2NDC(0.2+x2_coord-x1_coord);
						if (i_cent == 5)
						{
							legend->SetX1NDC(0.55);
							legend->SetX2NDC(0.55+x2_coord-x1_coord);
						}
					}



					g_ChPS_PbPb_final_sys_indR[i_trk][i_cent][i_jet]->GetYaxis()->SetNdivisions(505);
					g_ChPS_PbPb_final_sys_indR[i_trk][i_cent][i_jet]->SetFillColorAlpha(g_ChPS_PbPb_final_sys_indR[i_trk][i_cent][i_jet]->GetFillColor(),opacity);

					g_ChPS_pp_final_sys_indR[i_trk][6][i_jet]->GetYaxis()->SetNdivisions(505);
					g_ChPS_pp_final_sys_indR[i_trk][6][i_jet]->SetFillColorAlpha(g_ChPS_pp_final_sys_indR[i_trk][6][i_jet]->GetFillColor(),opacity);

					TGraphAsymmErrors* g_sys_pbpb = (TGraphAsymmErrors*)g_ChPS_PbPb_final_sys_indR[i_trk][i_cent][i_jet]->Clone("g_sys_PbPb");
					TGraphAsymmErrors* g_stat_pbpb = (TGraphAsymmErrors*)g_ChPS_PbPb_final_stat_indR[i_trk][i_cent][i_jet]->Clone("g_stat_PbPb");
					TGraphAsymmErrors* g_sys_pp = (TGraphAsymmErrors*)g_ChPS_pp_final_sys_indR[i_trk][6][i_jet]->Clone("g_sys_pp");
					TGraphAsymmErrors* g_stat_pp = (TGraphAsymmErrors*)g_ChPS_pp_final_stat_indR[i_trk][6][i_jet]->Clone("g_stat_pp");

					double shift_size = 0.0025;
					g_sys_pbpb = shift(g_sys_pbpb, trk_itr, shift_size);
					g_stat_pbpb = shift(g_stat_pbpb, trk_itr, shift_size);

					g_sys_pp = shift(g_sys_pp, trk_itr, shift_size);
					g_stat_pp = shift(g_stat_pp, trk_itr, shift_size);

					g_sys_pp->GetXaxis()->SetRangeUser(0, r_max_range);
					g_sys_pp->GetYaxis()->SetRangeUser(y_range_lo, y_range_hi);
					g_sys_pp->GetYaxis()->SetNdivisions(505);

					if (trk_itr == 0) g_sys_pp->Draw("a E2");
					else g_sys_pp->Draw("E2 same");
					g_stat_pp->Draw("P E1");
					g_sys_pbpb->Draw("E2");
					g_stat_pbpb->Draw("P E1");

					gPad->SetLogy();

					trk_itr++;
				} // end trk loop

				canvas->cd();
				ltx->SetTextSize(28);
				ltx->SetTextAlign(21);
				legend->SetHeader(" #scale[0.8]{#it{pp}}  #scale[0.8]{Pb+Pb}");
				legend->Draw();

				double x_left = 0.19, x_right = 0.92, y = 0.88, y_diff = 0.058;
				ltx->SetTextAlign(11);
				ltx->DrawLatexNDC(x_left, y, Form("%s", jet_label.c_str()));
				ltx->DrawLatexNDC(x_left, y-y_diff, Form("%s", centrality.c_str()));
				ltx->SetTextAlign(31);
				ltx->DrawLatexNDC(x_right, y, "#scale[1.5]{#font[72]{ATLAS}}");
				ltx->DrawLatexNDC(x_right, y=y-y_diff, "Pb+Pb #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}");
				ltx->DrawLatexNDC(x_right, y=y-y_diff, "#it{pp} #sqrt{#font[12]{s}} = 5.02 TeV, 25 pb^{-1}");
				ltx->DrawLatexNDC(x_right, y=y-y_diff, Form("anti-#font[12]{k}_{#font[12]{t}} R=0.4"));


				first_pass_cent = false;

				canvas->Print(Form("output_pdf_nominal/conf/DpT_dR_jet%i_cent%i.pdf", i_jet, i_cent));

			} //end cent loop

			jet_itr++;
		} //end jet loop
	}
*/

	{
		trk_select1 = 3;
		cout << "Doing Final RDpT spectra as function of centrality" << endl;

		TCanvas *canvas = new TCanvas("canvas","canvas",800,600);
		TLegend *legend = new TLegend(0.70,0.59,0.99,0.95,NULL,"brNDC");
		legend->SetTextFont(43);
		legend->SetBorderSize(0);
		legend->SetNColumns(2);
		legend->SetTextSize(24);
		legend->SetFillStyle(0);

		int jet_itr = 0;
		for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
		{
			string jet_label = Form("%1.0f < #it{p}_{T}^{jet} < %1.0f GeV", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));
//			if (!(i_jet == 7)) continue;

			canvas->cd();
			canvas->Clear();

			int trk_itr = 0;

			bool first_pass_cent = true;
			string centrality;

			for (int i_cent = 0; i_cent < 6; i_cent++)
			{
				centrality = num_to_cent(31,i_cent);


//			for (int i_trk = trk_select1; i_trk < 7; i_trk++)
//			{
//				if ((i_trk != trk_select1 && i_trk != trk_select2)) continue;
//				string trk_label = Form("%1.1f < #it{p}_{T} < %1.1f GeV", trkpT_binning->GetBinLowEdge(i_trk+1), trkpT_binning->GetBinUpEdge(i_trk+1));

//				bool first_pass_cent = true;
//				string centrality;

//				for (int i_cent = 0; i_cent < 6; i_cent++)
//				{
//					centrality = num_to_cent(31,i_cent);

				for (int i_trk = trk_select1; i_trk < 7; i_trk++)
				{
					if ((i_trk != trk_select1 && i_trk != trk_select2)) continue;
					string trk_label = Form("%1.1f < #it{p}_{T} < %1.1f GeV", trkpT_binning->GetBinLowEdge(i_trk+1), trkpT_binning->GetBinUpEdge(i_trk+1));

					SetHStyle(h_RDpT_final_ratio_indR[i_trk][i_cent][i_jet], i_cent);
					SetHStyle_graph(g_RDpT_final_sys_indR[i_trk][i_cent][i_jet], i_cent);

					if (i_trk == trk_select1) SetHStyle_graph(g_RDpT_final_stat_indR[i_trk][i_cent][i_jet], i_cent);
					if (i_trk == trk_select2) SetHStyle_graph_open(g_RDpT_final_stat_indR[i_trk][i_cent][i_jet], i_cent);

					if (jet_itr == 0)
					{
						name = Form("%1.1f-%1.1f",
									trkpT_binning->GetBinLowEdge(i_trk+1), trkpT_binning->GetBinUpEdge(i_trk+1));

//						if (i_trk == trk_select1) legend->AddEntry(g_RDpT_final_stat_indR[i_trk][i_cent][i_jet],name.c_str(), "p");
//						if (i_trk == trk_select2) legend->AddEntry(g_RDpT_final_stat_indR[i_trk][i_cent][i_jet],Form("%s [GeV], %s", name.c_str(), centrality.c_str()),"p");
						if (i_trk == trk_select1) legend->AddEntry(g_RDpT_final_stat_indR[i_trk][i_cent][i_jet],"  ", "p");
						if (i_trk == trk_select2) legend->AddEntry(g_RDpT_final_stat_indR[i_trk][i_cent][i_jet],Form("  %s",centrality.c_str()), "p");
					}


					g_RDpT_final_sys_indR[i_trk][i_cent][i_jet]->SetLineWidth(1.);
					g_RDpT_final_stat_indR[i_trk][i_cent][i_jet]->SetLineWidth(1.);

					y_range_lo = 0.;
					y_range_hi = 5.5;

					g_RDpT_final_sys_indR[i_trk][i_cent][i_jet]->SetLineColor(h_RDpT_final_ratio_indR[i_trk][i_cent][i_jet]->GetMarkerColor());
					g_RDpT_final_sys_indR[i_trk][i_cent][i_jet]->SetFillColorAlpha(g_RDpT_final_sys_indR[i_trk][i_cent][i_jet]->GetFillColor(),opacity);

					TGraphAsymmErrors* g_sys = (TGraphAsymmErrors*)g_RDpT_final_sys_indR[i_trk][i_cent][i_jet]->Clone("g_sys");
					TGraphAsymmErrors* g_stat = (TGraphAsymmErrors*)g_RDpT_final_stat_indR[i_trk][i_cent][i_jet]->Clone("g_stat");

					double shift_size = 0.0045;
					g_sys = shift(g_sys, i_cent, shift_size);
					g_stat = shift(g_stat, i_cent, shift_size);

					g_sys->GetXaxis()->SetLimits(0, r_max_range);
					g_sys->GetYaxis()->SetRangeUser(y_range_lo, y_range_hi);
					g_sys->GetYaxis()->SetNdivisions(505);

					if (trk_itr == 0 && i_cent == 0) g_sys->Draw("a E2");
					else g_sys->Draw("E2 same");
					g_stat->Draw("P E1");

					gPad->SetLogy(0);

					trk_itr++;

				} // end trk loop

				first_pass_cent = false;

			} //end cent loop

			canvas->cd();
			legend->SetHeader("#it{p}_{T,1}   #it{p}_{T,2}");
			legend->Draw();

			canvas->cd();
			line->DrawLine(0, 1, r_max_range, 1);

			double x_left = 0.19, x_right = 0.92, y = 0.88, y_diff = 0.058;
			ltx->SetTextAlign(11);
			ltx->SetTextSize(28);
			ltx->DrawLatexNDC(x_left, y, "#scale[1.5]{#font[72]{ATLAS}}");
			ltx->DrawLatexNDC(x_left, y=y-y_diff, "Pb+Pb #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}");
			ltx->DrawLatexNDC(x_left, y=y-y_diff, "#it{pp} #sqrt{#font[12]{s}} = 5.02 TeV, 25 pb^{-1}");
			ltx->DrawLatexNDC(x_left, y=y-y_diff, Form("%s", jet_label.c_str()));
			ltx->DrawLatexNDC(x_left, y=y-y_diff, Form("anti-#font[12]{k}_{#font[12]{t}} R=0.4"));
			ltx->DrawLatexNDC(x_left, y=y-y_diff, Form("%1.1f < #it{p}_{T,1} < %1.1f GeV",trkpT_binning->GetBinLowEdge(trk_select1+1), trkpT_binning->GetBinUpEdge(trk_select1+1)));
			ltx->DrawLatexNDC(x_left, y=y-y_diff, Form("%1.1f < #it{p}_{T,2} < %1.1f GeV",trkpT_binning->GetBinLowEdge(trk_select2+1), trkpT_binning->GetBinUpEdge(trk_select2+1)));

			canvas->Print(Form("output_pdf_nominal/conf/RDpT_dR_trk%i_trk%i_jet%i.pdf", trk_select1, trk_select2, i_jet));

			jet_itr++;
		} //end jet loop
		trk_select1 = 2;
	}

	cout << "######### DONE CONF PLOTS #########" << endl << endl;;


}

