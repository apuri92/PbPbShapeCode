#include "../functions/global_variables.h"

void get_posCorr(string config_file = "ff_config.cfg")
{
	cout << "######### GETTING POS_CORR FACTORS #########" << endl;

	SetAtlasStyle();
	string name;
	gErrorIgnoreLevel = 3001;

	//	##############	Reading config	##############"
	TEnv *m_config = new TEnv();
	m_config->ReadFile(config_file.c_str(), EEnvLevel(1));

	std::string dataset_type = "PbPb"; dataset_type = m_config->GetValue("dataset_type", dataset_type.c_str());
	int isMC = 1; isMC = m_config->GetValue("isMC", isMC);
	int sys_mode = -1; sys_mode = m_config->GetValue("sys_mode", sys_mode);

	int centrality_scheme = 31; centrality_scheme = m_config->GetValue("centrality_scheme", centrality_scheme);
	int verbose = 0; verbose = m_config->GetValue("verbose", verbose);

	std::string did = "data";
	if (isMC) did = "MC";

	if (verbose) m_config->Print();

	//	##############	Config done	##############"

	double r_max_range = 0.8;
	std::string sys_path = "";
	if (sys_mode == 0) sys_path = Form("nominal");
	else if (sys_mode > 50 || sys_mode < 0) sys_path = Form("c%i", sys_mode);
	else sys_path = Form("sys%i", sys_mode);

	name = Form("../raw_results/%s/FF_MC_out_histo_%s_5p02_r001.root", sys_path.c_str(), dataset_type.c_str());
	TFile *file = new TFile(name.c_str());

	cout << "Using files:" << endl;
	cout << file->GetName() << endl;

//	sys_path = Form("_%s", sys_path.c_str());
	TFile *output = new TFile(Form("output_pdf_%s/root/posCorr_factors_%s.root", sys_path.c_str(), dataset_type.c_str()),"recreate");

	TAxis* jet_pt_binning = (TAxis*)((TH3*)file->Get("h_reco_jet_spectrum_y0_cent0"))->GetXaxis();
	TAxis* trk_pt_binning = (TAxis*)((TH3*)file->Get("h_dR_change_jetpt0_cent0"))->GetZaxis();
	TAxis* dR_binning = (TAxis*)((TH3*)file->Get("h_dR_change_jetpt0_cent0"))->GetXaxis();

	bool doSmall;
	if (dataset_type == "pp") doSmall = false;
	if (dataset_type == "PbPb") doSmall = true;

	int ptJetBinsN = jet_pt_binning->GetNbins();
	int ptTrkBinsN = trk_pt_binning->GetNbins();
	int dRBinsN = dR_binning->GetNbins();

	TH3* h_dR_change;
	vector<vector<vector<TH2*>>> h_reco_truth =  vector<vector<vector<TH2*>>> (ptJetBinsN, vector<vector<TH2*>> (n_cent_cuts, vector<TH2*> (ptTrkBinsN)));
	vector<vector<vector<TH1*>>> h_reco_dR =  vector<vector<vector<TH1*>>> (ptJetBinsN, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (ptTrkBinsN)));
	vector<vector<vector<TH1*>>> h_truth_dR =  vector<vector<vector<TH1*>>> (ptJetBinsN, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (ptTrkBinsN)));
	vector<vector<vector<TH1*>>> h_ratio_dR =  vector<vector<vector<TH1*>>> (ptJetBinsN, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (ptTrkBinsN)));
	vector<vector<vector<TH1*>>> h_efficiency_dR =  vector<vector<vector<TH1*>>> (ptJetBinsN, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (ptTrkBinsN)));
	vector<vector<vector<TH1*>>> h_purity_dR =  vector<vector<vector<TH1*>>> (ptJetBinsN, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (ptTrkBinsN)));

	vector<vector<vector<TH1*>>> h_check =  vector<vector<vector<TH1*>>> (ptJetBinsN, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (dRBinsN)));

	for (int i_cent = 0; i_cent < n_cent_cuts; i_cent++)
	{
		for (int i_jetpt = 0; i_jetpt < ptJetBinsN; i_jetpt++)
		{
			h_dR_change = (TH3*)file->Get(Form("h_dR_change_jetpt%i_cent%i",i_jetpt,i_cent));
			h_dR_change->Sumw2();

			//Initialize check histograms
			for (int i_dR = 0; i_dR < dRBinsN; i_dR++)
			{
				h_check.at(i_jetpt).at(i_cent).at(i_dR) = (TH1*)h_dR_change->Project3D("z")->Clone(Form("c%i_j%i_r%i",i_cent, i_jetpt, i_dR));
//				h_check.at(i_jetpt).at(i_cent).at(i_dR)->SetName(Form("c%i_j%i_r%i",i_cent, i_jetpt, i_dR));
				h_check.at(i_jetpt).at(i_cent).at(i_dR)->Reset();
			}

			double pt_jet_lo = jet_pt_binning->GetBinLowEdge(i_jetpt+1);
			double pt_jet_hi = jet_pt_binning->GetBinUpEdge(i_jetpt+1);

			for (int i_pt_bin = 0; i_pt_bin < ptTrkBinsN; i_pt_bin++)
			{
				double pt_lo = h_dR_change->GetZaxis()->GetBinLowEdge(i_pt_bin+1);
				double pt_hi = h_dR_change->GetZaxis()->GetBinUpEdge(i_pt_bin+1);

				bool write_to_screen = false;

				h_dR_change->GetZaxis()->SetRange(i_pt_bin+1,i_pt_bin+1);

				//2D Response
				h_reco_truth.at(i_jetpt).at(i_cent).at(i_pt_bin) = (TH2*)h_dR_change->Project3D("yx");
				h_reco_truth.at(i_jetpt).at(i_cent).at(i_pt_bin)->SetName(Form("2D_c%i_j%i_trk%i", i_cent, i_jetpt, i_pt_bin));
				h_reco_truth.at(i_jetpt).at(i_cent).at(i_pt_bin)->Sumw2();
				name = Form("2D Resp. Cent: %s, %4.0f < #it{p}_{T}^{jet} < %4.0f GeV, %4.1f < #it{p}_{T}^{trk} < %4.1f", num_to_cent(31,i_cent).c_str(), pt_jet_lo, pt_jet_hi, pt_lo, pt_hi);
				h_reco_truth.at(i_jetpt).at(i_cent).at(i_pt_bin)->SetTitle(name.c_str());

				output->cd();
				h_reco_truth.at(i_jetpt).at(i_cent).at(i_pt_bin)->Write(Form("h_reco_truth_c%i_j%i_trk%i", i_cent, i_jetpt, i_pt_bin));

				//Reco projection
				h_reco_dR.at(i_jetpt).at(i_cent).at(i_pt_bin) = (TH1*)h_reco_truth.at(i_jetpt).at(i_cent).at(i_pt_bin)->ProjectionY(Form("reco_c%i_j%i_trk%i", i_cent, i_jetpt, i_pt_bin));
				h_reco_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->Sumw2();
				name = Form("Reco Proj. Cent: %s, %4.0f < #it{p}_{T}^{jet} < %4.0f GeV, %4.1f < #it{p}_{T}^{trk} < %4.1f", num_to_cent(31,i_cent).c_str(), pt_jet_lo, pt_jet_hi, pt_lo, pt_hi);
				h_reco_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->SetTitle(name.c_str());
				name = Form("Reco_dR_c%i_j%i_trk%i",i_cent, i_jetpt, i_pt_bin);
				h_reco_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->SetName(name.c_str());
				output->cd();
				h_reco_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->Write(Form("h_reco_c%i_j%i_trk%i", i_cent, i_jetpt, i_pt_bin));

				//Truth projection
				h_truth_dR.at(i_jetpt).at(i_cent).at(i_pt_bin) = (TH1*)h_reco_truth.at(i_jetpt).at(i_cent).at(i_pt_bin)->ProjectionX(Form("truth_c%i_j%i_trk%i", i_cent, i_jetpt, i_pt_bin));
				h_truth_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->Sumw2();
				name = Form("Truth Proj. Cent: %s, %4.0f < #it{p}_{T}^{jet} < %4.0f GeV, %4.1f < #it{p}_{T}^{trk} < %4.1f", num_to_cent(31,i_cent).c_str(), pt_jet_lo, pt_jet_hi, pt_lo, pt_hi);
				h_truth_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->SetTitle(name.c_str());
				name = Form("Truth_dR_c%i_j%i_trk%i",i_cent, i_jetpt, i_pt_bin);
				h_truth_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->SetName(name.c_str());
				output->cd();
				h_truth_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->Write(Form("h_truth_c%i_j%i_trk%i", i_cent, i_jetpt, i_pt_bin));

				//scale by annulus area
				double orig, orig_err;

				for (int i_dR = 1; i_dR <= h_truth_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->GetNbinsX(); i_dR++)
				{
					double lo = h_truth_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->GetXaxis()->GetBinLowEdge(i_dR);
					double hi = h_truth_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->GetXaxis()->GetBinUpEdge(i_dR);
					double area = TMath::Pi() * (pow(hi,2) - pow(lo,2));

					orig = h_reco_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->GetBinContent(i_dR);
					orig_err = h_reco_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->GetBinError(i_dR);
					h_reco_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->SetBinContent(i_dR, orig/area);
					h_reco_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->SetBinError(i_dR, orig_err/area);

					orig = h_truth_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->GetBinContent(i_dR);
					orig_err = h_truth_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->GetBinError(i_dR);
					h_truth_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->SetBinContent(i_dR, orig/area);
					h_truth_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->SetBinError(i_dR, orig_err/area);
				}


				//Ratios
				h_ratio_dR.at(i_jetpt).at(i_cent).at(i_pt_bin) = (TH1*)h_truth_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->Clone(Form("ratio_c%i_j%i_trk%i", i_cent, i_jetpt, i_pt_bin)); //doing truth/reco. Need to multiply this with the ChPS
				h_ratio_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->Divide(h_reco_dR.at(i_jetpt).at(i_cent).at(i_pt_bin));

				//error calculation: sqrt(T^2/R^3 ( 1 - M^2/(T*R) ))
//				for (int i = 0; i < h_ratio_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->GetXaxis()->GetNbins(); i++)
//				{
//					double T = h_truth_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->GetBinContent(i+1);
//					double R = h_reco_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->GetBinContent(i+1);
//					double M = h_reco_truth.at(i_jetpt).at(i_cent).at(i_pt_bin)->GetBinContent(i+1,i+1);
//					if (R == 0 || T == 0) continue;
//					double error = sqrt(pow(T,2)/pow(R,3));
//					double error = sqrt(((pow(M,2)/(T*R))));
//					double error = sqrt(pow(T,2)/pow(R,3) * (1-(pow(M,2)/(T*R) ) ) );
//					cout << error << endl;
//					h_ratio_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->SetBinError(i+1,error);
//				}

				name = Form("Ratio Cent: %s, %4.0f < #it{p}_{T}^{jet} < %4.0f GeV, %4.1f < #it{p}_{T}^{trk} < %4.1f", num_to_cent(31,i_cent).c_str(), pt_jet_lo, pt_jet_hi, pt_lo, pt_hi);
				h_ratio_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->SetTitle(name.c_str());

				name = Form("Ratio_dR_c%i_j%i_trk%i",i_cent, i_jetpt, i_pt_bin);
				output->cd();
				h_ratio_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->SetName(Form("ratio_c%i_j%i_trk%i", i_cent, i_jetpt, i_pt_bin));
				h_ratio_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->Write();

				//Efficiency and purity
				h_efficiency_dR.at(i_jetpt).at(i_cent).at(i_pt_bin) = (TH1*)h_ratio_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->Clone(Form("eff_c%i_j%i_trk%i", i_cent, i_jetpt, i_pt_bin));
				h_efficiency_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->Reset();
				h_purity_dR.at(i_jetpt).at(i_cent).at(i_pt_bin) = (TH1*)h_reco_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->Clone(Form("pur_c%i_j%i_trk%i", i_cent, i_jetpt, i_pt_bin));
				h_purity_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->Reset();

				TH1* h_tmp = (TH1*)h_purity_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->Clone("h_tmp");
				h_tmp->Reset();
//				cout << Form("cent%i_jet%i_trk%i",i_cent, i_jetpt, i_pt_bin) << endl;
				for (int i_dR = 1; i_dR <= dRBinsN; i_dR++)
				{
					float num = h_reco_truth.at(i_jetpt).at(i_cent).at(i_pt_bin)->GetBinContent(i_dR,i_dR);
					float den_p = 0.;
					float den_e = 0.;

					float M = h_reco_truth.at(i_jetpt).at(i_cent).at(i_pt_bin)->GetBinContent(i_dR,i_dR);

					h_tmp->SetBinContent(i_dR,M);

					for (int j = 1; j <= dRBinsN; j++)
					{
						den_p = den_p + h_reco_truth.at(i_jetpt).at(i_cent).at(i_pt_bin)->GetBinContent(i_dR, j);
						den_e = den_e + h_reco_truth.at(i_jetpt).at(i_cent).at(i_pt_bin)->GetBinContent(j, i_dR);
					}

					if (den_e == 0 || den_p == 0)  continue;
					h_purity_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->SetBinContent(i_dR,num/den_p);
					h_efficiency_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->SetBinContent(i_dR,num/den_e);

					h_check.at(i_jetpt).at(i_cent).at(i_dR-1)->SetBinContent(i_pt_bin+1,h_ratio_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->GetBinContent(i_dR));
					h_check.at(i_jetpt).at(i_cent).at(i_dR-1)->SetBinError(i_pt_bin+1,h_ratio_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->GetBinError(i_dR));

				} //End dR loop
				output->cd();
				h_tmp->Write(Form("cent%i_jet%i_trk%i",i_cent, i_jetpt, i_pt_bin));

			} //End trk pT loop
		} // End jet pT loop
	} //End centrality loop


	for (int i_cent = 0; i_cent < n_cent_cuts; i_cent++)
	{
		string cent = num_to_cent(31,i_cent);
		for (int i_jetpt = 0; i_jetpt < ptJetBinsN; i_jetpt++)
		{
			double pt_lo = jet_pt_binning->GetBinLowEdge(i_jetpt+1);
			double pt_hi = jet_pt_binning->GetBinUpEdge(i_jetpt+1);

			for (int i_dR = 0; i_dR < dRBinsN; i_dR++)
			{
				double dR_lo = dR_binning->GetBinLowEdge(i_dR+1);
				double dR_hi = dR_binning->GetBinUpEdge(i_dR+1);

				output->cd();

				name = Form("Check: %s, %4.1f < #it{p}_{T}^{jet} < %4.2f, %4.2f < dR < %4.2f", cent.c_str(), pt_lo, pt_hi, dR_lo, dR_hi);
				h_check.at(i_jetpt).at(i_cent).at(i_dR)->SetTitle(name.c_str());

				name = Form("h_bin_by_bin_cent%i_jetpt%i_dR%i",i_cent, i_jetpt, i_dR);
				h_check.at(i_jetpt).at(i_cent).at(i_dR)->Write(name.c_str());

			}
		}
	}





	TLegend *legend1 = new TLegend();
	legend1->SetBorderSize(0);
	legend1->SetTextFont(43);
	legend1->SetTextSize(12);

	TLegend *legend2 = new TLegend();
	legend2->SetBorderSize(0);
	legend2->SetTextFont(43);
	legend2->SetTextSize(12);

	TLegend *legend3 = new TLegend();
	legend3->SetBorderSize(0);
	legend3->SetTextFont(43);
	legend3->SetTextSize(12);

	TLegend *legend4 = new TLegend(0.19, 0.20, 0.70, 0.45);
	legend4->SetBorderSize(0);
	legend4->SetTextFont(43);
	legend4->SetTextSize(12);

	TLegend *legend5 = new TLegend(0.19, 0.20, 0.70, 0.45);
	legend5->SetBorderSize(0);
	legend5->SetTextFont(43);
	legend5->SetTextSize(12);

	TLatex *ltx = new TLatex();
	ltx->SetTextFont(43);
	ltx->SetTextSize(18);
	ltx->SetNDC();

	TLine *line = new TLine();
	line->SetLineStyle(2);

	TCanvas *c0 = new TCanvas("c0","c0",1200,600);
	TCanvas *c1 = new TCanvas("c1","c1",1200,600);
	TCanvas *c2 = new TCanvas("c2","c2",1200,600);
	TCanvas *c3 = new TCanvas("c3","c3",1200,600);
	TCanvas *c4 = new TCanvas("c4","c4",1200,600);
	TCanvas *c5 = new TCanvas("c5","c5",1200,600);

	bool first_cent_pass = true;
	for (int i_cent = 0; i_cent < n_cent_cuts; i_cent++)
	{
		if (dataset_type == "PbPb" && i_cent == 6) continue;
		if (dataset_type == "pp" && i_cent < 6) continue;

		c1->cd();
		c1->Clear();
		c1->Divide(4,2);

		c2->cd();
		c2->Clear();
		c2->Divide(4,2);

		c3->cd();
		c3->Clear();
		c3->Divide(4,2);

		c4->cd();
		c4->Clear();
		c4->Divide(4,2);

		c5->cd();
		c5->Clear();
		c5->Divide(4,2);

		int jet_iter = 0;

		for (int i_jetpt = 0; i_jetpt < ptJetBinsN; i_jetpt++)
		{
			double pt_jet_lo = jet_pt_binning->GetBinLowEdge(i_jetpt+1);
			double pt_jet_hi = jet_pt_binning->GetBinUpEdge(i_jetpt+1);

			if (pt_jet_hi < 130. ||
				pt_jet_lo > 300) continue;

			c0->cd();
			c0->Clear();
			c0->Divide(4,2);

			int trk_iter = 0;

			for (int i_pt_bin = 0; i_pt_bin < ptTrkBinsN; i_pt_bin++)
			{
				double pt_lo = trk_pt_binning->GetBinLowEdge(i_pt_bin+1);
				double pt_hi = trk_pt_binning->GetBinUpEdge(i_pt_bin+1);

				if (pt_lo < 1.||
					pt_hi> 80.) continue;

				//Drawing 2D response plots, each page is a centrality and jet pt bin, with each pad being a trk bin
				{
					c0->cd(trk_iter+1);
					h_reco_truth.at(i_jetpt).at(i_cent).at(i_pt_bin)->GetXaxis()->SetTitle("Truth r");
					h_reco_truth.at(i_jetpt).at(i_cent).at(i_pt_bin)->GetYaxis()->SetTitle("Reco r");
					h_reco_truth.at(i_jetpt).at(i_cent).at(i_pt_bin)->GetXaxis()->SetRangeUser(0,r_max_range);
					h_reco_truth.at(i_jetpt).at(i_cent).at(i_pt_bin)->GetYaxis()->SetRangeUser(0,r_max_range);
					h_reco_truth.at(i_jetpt).at(i_cent).at(i_pt_bin)->GetZaxis()->SetRangeUser(1E-8,1E0);

//					h_reco_truth.at(i_jetpt).at(i_cent).at(i_pt_bin)->GetZaxis()->SetRangeUser(1e-9,1);
					gPad->SetRightMargin(0.15);
					h_reco_truth.at(i_jetpt).at(i_cent).at(i_pt_bin)->Draw("colz");
					ltx->SetTextSize(12);
					ltx->SetTextAlign(12);
					ltx->DrawLatexNDC(0.19,0.88,Form("%4.1f < #it{p}_{T}^{trk} < %4.1f",pt_lo, pt_hi));
					ltx->DrawLatexNDC(0.19,0.965,Form("%s: %4.0f < #it{p}_{T}^{jet} < %4.0f GeV",num_to_cent(31,i_cent).c_str(), pt_jet_lo, pt_jet_hi));
					line->DrawLine(0,0,r_max_range,r_max_range);
					gPad->SetLogz();
				}

				//Drawing 1D reco response projections, each page is a centrality, each pad is a trk pt bin, each pad has different curves for jet pt
				{
					c1->cd(trk_iter+1);
					h_reco_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->GetXaxis()->SetTitle("Reco dR");
					SetHStyle(h_reco_dR.at(i_jetpt).at(i_cent).at(i_pt_bin),jet_iter);
					smallify(h_reco_dR.at(i_jetpt).at(i_cent).at(i_pt_bin));
					h_reco_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->GetYaxis()->SetRangeUser(1e-7, 1);
					h_reco_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->GetXaxis()->SetRangeUser(0,r_max_range);
					if (trk_iter == 0) legend1->AddEntry(h_reco_dR.at(i_jetpt).at(i_cent).at(i_pt_bin), Form("%4.0f < #it{p}_{T}^{jet} < %4.0f GeV", pt_jet_lo, pt_jet_hi));
					if (jet_iter == 0) h_reco_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->Draw("p");
					else h_reco_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->Draw("same p");
					if (jet_iter == 0) ltx->DrawLatexNDC(0.20,0.88,Form("%4.1f < #it{p}_{T}^{trk} < %4.1f",pt_lo, pt_hi));
					gPad->SetLogy();
				}

				//Drawing 1D reco response projections, each page is a centrality, each pad is a trk pt bin, each pad has different curves for jet pt
				{
					c2->cd(trk_iter+1);
					h_truth_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->GetXaxis()->SetTitle("Truth dR");
					SetHStyle(h_truth_dR.at(i_jetpt).at(i_cent).at(i_pt_bin),jet_iter);
					smallify(h_truth_dR.at(i_jetpt).at(i_cent).at(i_pt_bin));
					h_truth_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->GetYaxis()->SetRangeUser(1e-7, 1);
					h_truth_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->GetXaxis()->SetRangeUser(0,r_max_range);
					if (trk_iter == 0) legend2->AddEntry(h_truth_dR.at(i_jetpt).at(i_cent).at(i_pt_bin), Form("%4.0f < #it{p}_{T}^{jet} < %4.0f GeV", pt_jet_lo, pt_jet_hi));
					if (jet_iter == 0) h_truth_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->Draw("p");
					else h_truth_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->Draw("same p");
					if (jet_iter == 0) ltx->DrawLatexNDC(0.20,0.88,Form("%4.1f < #it{p}_{T}^{trk} < %4.1f",pt_lo, pt_hi));
					gPad->SetLogy();
				}

				//Drawing 1D ratio of response projections, each page is a centrality, each pad is a trk pt bin, each pad has different curves for jet pt
				{
					c3->cd(trk_iter+1);
					h_ratio_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->GetXaxis()->SetTitle("#it{r}");
					h_ratio_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->GetYaxis()->SetTitle("Correction Factor");
					SetHStyle(h_ratio_dR.at(i_jetpt).at(i_cent).at(i_pt_bin),jet_iter);
					smallify(h_ratio_dR.at(i_jetpt).at(i_cent).at(i_pt_bin));
					double y_max = h_ratio_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->GetBinContent(h_ratio_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->GetMaximumBin());
					double y_min = h_ratio_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->GetBinContent(h_ratio_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->GetMinimumBin());
					double diff = -1;
					if (fabs(y_max - 1)> fabs(1 - y_min)) diff = y_max - 1 + 0.07;
					else diff = 1 - y_min + 0.07;
					if (diff> 1) diff = 0.5;
					h_ratio_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->GetYaxis()->SetRangeUser(1 - diff, 1 + diff);
//					h_ratio_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->GetYaxis()->SetRangeUser(-4,5);
					h_ratio_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->GetXaxis()->SetRangeUser(0,r_max_range);
					if (trk_iter == 0) legend3->AddEntry(h_ratio_dR.at(i_jetpt).at(i_cent).at(i_pt_bin), Form("%4.0f < #it{p}_{T}^{jet} < %4.0f GeV", pt_jet_lo, pt_jet_hi));

					// removing points with poor closure (this happens for high z particles at the jet edge
					//					Jet pT [GeV]	Track pT [GeV]	R
					//					126-158 GeV	6-10			r > 0.3 (jet7_trk6_dR6)
					//								10-25 GeV	r > 0.3 (jet7_trk7_dR6)
					//								25-63 GeV	r > 0.2 (jet7_trk8_dR4)
					//					158-200 GeV	10-25 GeV	r > 0.4 (jet8_trk7_dR7)
					//								25-63 GeV	r > 0.3 (jet8_trk8_dR6)
					//					200-251 GeV	25-63 GeV	r > 0.3 (jet9_trk8_dR6)
					//					251-316 GeV	x	x

					for (int i_dR = 0; i_dR < h_ratio_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->GetXaxis()->GetNbins(); i_dR++)
					{
						if ((i_jetpt == 7 && i_pt_bin == 6 && i_dR >= 6)
							||
							(i_jetpt == 7 && i_pt_bin == 7 && i_dR >= 6)
							||
							(i_jetpt == 7 && i_pt_bin == 8 && i_dR >= 4)
							||
							(i_jetpt == 8 && i_pt_bin == 7 && i_dR >= 7)
							||
							(i_jetpt == 8 && i_pt_bin == 8 && i_dR >= 6)
							||
							(i_jetpt == 9 && i_pt_bin == 8 && i_dR >= 6)
							)
						{
							h_ratio_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->SetBinContent(i_dR+1, -1);
							h_purity_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->SetBinContent(i_dR+1, -1);
							h_efficiency_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->SetBinContent(i_dR+1, -1);
						}

					}


					if (jet_iter == 0) h_ratio_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->Draw("pe");
					else h_ratio_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->Draw("same pe");
					if (jet_iter == 0)
					{
						ltx->SetTextAlign(12);
						ltx->DrawLatexNDC(0.20,0.90,Form("%4.1f < #it{p}_{T}^{trk} < %4.1f GeV",pt_lo, pt_hi));
					}
					if (jet_iter == 0) line->DrawLine(0,1,r_max_range,1);
					ltx->SetTextAlign(32);
					if (jet_iter == 0) ltx->DrawLatexNDC(0.90,0.90,num_to_cent(31,i_cent).c_str());

				}

				//Drawing purity of response, each page is a centrality, each pad is a trk pt bin, each pad has different curves for jet pt
				{
					c4->cd(trk_iter+1);
					h_purity_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->SetTitle(0);
					h_purity_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->GetXaxis()->SetTitle(0);
					h_purity_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->GetYaxis()->SetTitle(0);
					h_purity_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->GetXaxis()->SetTitle("Reco dR");
					h_purity_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->GetYaxis()->SetTitle("Purity");
					SetHStyle(h_purity_dR.at(i_jetpt).at(i_cent).at(i_pt_bin),jet_iter);
					smallify(h_purity_dR.at(i_jetpt).at(i_cent).at(i_pt_bin));
					h_purity_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->GetYaxis()->SetRangeUser(0,1.2);
					h_purity_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->GetXaxis()->SetRangeUser(0,r_max_range);
					if (trk_iter == 0) legend4->AddEntry(h_purity_dR.at(i_jetpt).at(i_cent).at(i_pt_bin), Form("%4.0f < #it{p}_{T}^{jet} < %4.0f GeV", pt_jet_lo, pt_jet_hi));
					if (jet_iter == 0) h_purity_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->Draw("p");
					else h_purity_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->Draw("same p");
					if (jet_iter == 0) line->DrawLine(0,0.5,1.2,0.5);
					if (jet_iter == 0) line->DrawLine(0,1,r_max_range,1);
					if (jet_iter == 0)
					{
						ltx->SetTextAlign(32);
						ltx->DrawLatexNDC(0.95,0.90,num_to_cent(31,i_cent).c_str());
						ltx->SetTextAlign(12);
						ltx->DrawLatexNDC(0.17,0.90,Form("%4.1f < #it{p}_{T}^{trk} < %4.1f",pt_lo, pt_hi));
					}
				}

				//Drawing efficienct of response projections, each page is a centrality, each pad is a trk pt bin, each pad has different curves for jet pt
				{
					c5->cd(trk_iter+1);
					h_efficiency_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->GetXaxis()->SetTitle("Truth dR");
					h_efficiency_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->GetYaxis()->SetTitle("Response Efficiency");
					SetHStyle(h_efficiency_dR.at(i_jetpt).at(i_cent).at(i_pt_bin),jet_iter);
					smallify(h_efficiency_dR.at(i_jetpt).at(i_cent).at(i_pt_bin));
					h_efficiency_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->GetYaxis()->SetRangeUser(0,1.2);
					h_efficiency_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->GetXaxis()->SetRangeUser(0,r_max_range);
					if (trk_iter == 0) legend5->AddEntry(h_efficiency_dR.at(i_jetpt).at(i_cent).at(i_pt_bin), Form("%4.0f < #it{p}_{T}^{jet} < %4.0f GeV", pt_jet_lo, pt_jet_hi));
					if (jet_iter == 0) h_efficiency_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->Draw("p");
					else h_efficiency_dR.at(i_jetpt).at(i_cent).at(i_pt_bin)->Draw("same p");
					if (jet_iter == 0) line->DrawLine(0,0.5,r_max_range,0.5);
					if (jet_iter == 0) line->DrawLine(0,1,r_max_range,1);
					if (jet_iter == 0)
					{
						ltx->SetTextAlign(32);
						ltx->DrawLatexNDC(0.95,0.90,num_to_cent(31,i_cent).c_str());
						ltx->SetTextAlign(12);
						ltx->DrawLatexNDC(0.17,0.90,Form("%4.1f < #it{p}_{T}^{trk} < %4.1f",pt_lo, pt_hi));
					}

				}


				trk_iter++;
			} //End track pt loop


			//2D Response
			{
				c0->cd(1);
				ltx->SetTextAlign(12);
				if (first_cent_pass && jet_iter == 0) name = "(";
				else name = "";
				c0->Print(Form("output_pdf_%s/%s/ShapeResponse2D_%s.pdf%s",sys_path.c_str(), dataset_type.c_str(), dataset_type.c_str(), name.c_str()),Form("Title: c%i_j%i", i_cent, i_jetpt));
			}


			jet_iter++;
		} //End jet pt loop


		// Reco projection
		{
			c1->cd();
			ltx->DrawLatexNDC(0.055,0.975,num_to_cent(31,i_cent).c_str());

			c1->cd(1);
			legend1->SetX1NDC(0.2);
			legend1->SetX2NDC(0.7);
			legend1->SetY1NDC(0.2);
			legend1->SetY2NDC(0.4);
			legend1->Draw();
			if (first_cent_pass) name = "(";
			else name = "";
			c1->Print(Form("output_pdf_%s/%s/RecoProj_%s.pdf%s", sys_path.c_str(), dataset_type.c_str(), dataset_type.c_str(), name.c_str()),Form("Title: c%i", i_cent));
			legend1->Clear();
		}
		// Truth projection
		{
			c2->cd();
			ltx->DrawLatexNDC(0.055,0.975,num_to_cent(31,i_cent).c_str());

			c2->cd(1);
			legend2->SetX1NDC(0.2);
			legend2->SetX2NDC(0.7);
			legend2->SetY1NDC(0.2);
			legend2->SetY2NDC(0.4);
			legend2->Draw();
			name = Form("c%i", i_cent);
			if (first_cent_pass) name = "(";
			else name = "";
			c2->Print(Form("output_pdf_%s/%s/TruthProj_%s.pdf%s", sys_path.c_str(), dataset_type.c_str(), dataset_type.c_str(), name.c_str()),Form("Title: c%i", i_cent));
			legend2->Clear();
		}

		//Ratio plots
		{
			c3->cd();
//			ltx->DrawLatexNDC(0.055,0.975,num_to_cent(31,i_cent).c_str());

			c3->cd(1);
			legend3->SetX1NDC(0.2);
			legend3->SetX2NDC(0.7);
			legend3->SetY1NDC(0.2);
			legend3->SetY2NDC(0.4);
			legend3->Draw();
			if (first_cent_pass) name = "(";
			else name = "";
			c3->Print(Form("output_pdf_%s/%s/RatioProj_%s.pdf%s", sys_path.c_str(), dataset_type.c_str(), dataset_type.c_str(), name.c_str()),Form("Title: c%i", i_cent));
			legend3->Clear();
		}

		{
			c4->cd();
			ltx->SetTextAlign(32);

			c4->cd(1);
			legend4->Draw();
			if (first_cent_pass) name = "(";
			else name = "";
			c4->Print(Form("output_pdf_%s/%s/RespPurity_%s.pdf%s", sys_path.c_str(), dataset_type.c_str(), dataset_type.c_str(), name.c_str()),Form("Title: c%i", i_cent));
			legend4->Clear();
		}

		{
			c5->cd(1);
			legend5->Draw();
			if (first_cent_pass) name = "(";
			else name = "";
			c5->Print(Form("output_pdf_%s/%s/RespEfficiency_%s.pdf%s", sys_path.c_str(), dataset_type.c_str(), dataset_type.c_str(), name.c_str()),Form("Title: c%i", i_cent));
			legend5->Clear();
		}

		first_cent_pass = false;
	} //End centrality loop

	c0->Clear();
	c1->Clear();
	c2->Clear();
	c3->Clear();
	c4->Clear();
	c5->Clear();

	c0->Print(Form("output_pdf_%s/%s/ShapeResponse2D_%s.pdf)", sys_path.c_str(), dataset_type.c_str(), dataset_type.c_str()),"Title: End");
	c1->Print(Form("output_pdf_%s/%s/RecoProj_%s.pdf)", sys_path.c_str(), dataset_type.c_str(), dataset_type.c_str()),"Title: End");
	c2->Print(Form("output_pdf_%s/%s/TruthProj_%s.pdf)", sys_path.c_str(), dataset_type.c_str(), dataset_type.c_str()),"Title: End");
	c3->Print(Form("output_pdf_%s/%s/RatioProj_%s.pdf)", sys_path.c_str(), dataset_type.c_str(), dataset_type.c_str()),"Title: End");
	c4->Print(Form("output_pdf_%s/%s/RespPurity_%s.pdf)", sys_path.c_str(), dataset_type.c_str(), dataset_type.c_str()),"Title: End");
	c5->Print(Form("output_pdf_%s/%s/RespEfficiency_%s.pdf)", sys_path.c_str(), dataset_type.c_str(), dataset_type.c_str()),"Title: End");



	TCanvas *c6 = new TCanvas("c6","c6",900,600);
	if (dataset_type == "pp") c6->SetCanvasSize(800,600);

	c6->cd(); ltx->DrawLatexNDC(0.5,0.5,"Test");

	TLegend *legend6 = new TLegend();

//	legend6->SetNColumns(2);
	legend6->SetBorderSize(0);
	legend6->SetTextFont(43);
	legend6->SetTextSize(9);
	ltx->SetTextSize(12);
	if (dataset_type == "pp") ltx->SetTextSize(16);

	if (dataset_type == "pp") legend6->SetTextSize(16);
	if (dataset_type == "PbPb") legend6->SetTextSize(8);


	for (int i_dR = 0; i_dR < dRBinsN; i_dR++)
	{
		c6->cd();
		c6->Clear();
		if (dataset_type == "PbPb") c6->Divide(3,2);

		double dR_lo = dR_binning->GetBinLowEdge(i_dR+1);
		double dR_hi = dR_binning->GetBinUpEdge(i_dR+1);

		bool first_cent_pass = true;
		for (int i_cent = 0; i_cent < n_cent_cuts; i_cent++)
		{
			if (dataset_type == "pp") c6->cd();
			if (dataset_type == "PbPb") c6->cd(i_cent+1);
//			c6->cd(i_cent+1);
//			gPad->SetTopMargin(0.07);
//			gPad->SetRightMargin(0);

			int jet_iter = 0;
			for (int i_jetpt = 0; i_jetpt < ptJetBinsN; i_jetpt++)
			{
				double pt_jet_lo = jet_pt_binning->GetBinLowEdge(i_jetpt+1);
				double pt_jet_hi = jet_pt_binning->GetBinUpEdge(i_jetpt+1);

				if (pt_jet_hi < 127. ||
					pt_jet_lo> 315.) continue;

				if (first_cent_pass) legend6->AddEntry(h_check.at(i_jetpt).at(i_cent).at(i_dR), Form("%4.0f < #it{p}_{T}^{jet} < %4.0f GeV", pt_jet_lo, pt_jet_hi));
				h_check.at(i_jetpt).at(i_cent).at(i_dR)->GetYaxis()->SetRangeUser(0,2);
				h_check.at(i_jetpt).at(i_cent).at(i_dR)->GetXaxis()->SetRangeUser(1,100);

				h_check.at(i_jetpt).at(i_cent).at(i_dR)->GetYaxis()->SetTitle("Correction Factors");
				h_check.at(i_jetpt).at(i_cent).at(i_dR)->GetXaxis()->SetTitle("#it{p}_{T}^{trk} [GeV]");
				SetHStyle_smallify(h_check.at(i_jetpt).at(i_cent).at(i_dR),jet_iter, doSmall);
//				smallify(h_check.at(i_jetpt).at(i_cent).at(i_dR));
				if (jet_iter == 0) h_check.at(i_jetpt).at(i_cent).at(i_dR)->Draw("p");
				else h_check.at(i_jetpt).at(i_cent).at(i_dR)->Draw("same p");
				gPad->SetLogx();

				jet_iter++;
			} //End jet loop

			c6->cd(i_cent+1);
			ltx->SetTextAlign(32);
			ltx->DrawLatexNDC(0.95,0.97,num_to_cent(31, i_cent).c_str());
			first_cent_pass = false;

			if (dataset_type == "pp") c6->cd();
			if (dataset_type == "PbPb") c6->cd(i_cent+1);
			ltx->SetTextAlign(12);
			ltx->DrawLatexNDC(0.19,0.97,Form("%1.2f < r < %1.2f", dR_lo, dR_hi));


		} //End cent loop


		c6->cd(1);
		legend6->SetX1NDC(0.18);
		legend6->SetX2NDC(0.60);
		legend6->SetY1NDC(0.62);
		legend6->SetY2NDC(0.90);
		legend6->Draw();


		if (i_dR == 0) name = "(";
		else name = "";
		c6->Print(Form("output_pdf_%s/%s/pos_corr_factors_%s.pdf%s", sys_path.c_str(), dataset_type.c_str(), dataset_type.c_str(), name.c_str()),Form("Title: dR%i",i_dR));


		legend6->Clear();
	} //End dr loop

	c6->Clear();
	name = ")";
	c6->Print(Form("output_pdf_%s/%s/pos_corr_factors_%s.pdf%s", sys_path.c_str(), dataset_type.c_str(), dataset_type.c_str(), name.c_str()),Form("Title: End"));

	cout << "######### DONE POS_CORR FACTORS #########" << endl << endl;


}

