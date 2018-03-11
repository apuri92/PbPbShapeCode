#include "functions/global_variables.h"

void combine_eff_trketa(vector<TFile*>& theFiles, vector<double> w)
{
	gErrorIgnoreLevel = 3001;
	SetAtlasStyle();
	std::string name;
	int nFiles = theFiles.size();
	
	//	##############	Reading config	##############"
	TEnv *m_config = new TEnv();
	m_config->ReadFile("perf_config.cfg", EEnvLevel(1));

	std::string dataset_type = "PbPb"; dataset_type = m_config->GetValue("dataset_type", dataset_type.c_str());
	std::string tracking_cut = "ppTight"; tracking_cut = m_config->GetValue("tracking_cut", tracking_cut.c_str());
	int centrality_scheme = 31; centrality_scheme = m_config->GetValue("centrality_scheme", centrality_scheme);
	//	##############	Config done	##############"

	//output root file
	name = Form("mc_eff_%s_trketa_jetptinc_%s.root",dataset_type.c_str(), tracking_cut.c_str());
	TFile *output_file = new TFile(Form("%s",name.c_str()),"recreate");
	cout << Form("Efficiency root file: %s", output_file->GetName()) << endl;

    //setting up fine eta for eff_trketa corrections
    vector<double> eta_range;
    double x = -2.5;
    while (x <= 2.5)
    {
        eta_range.push_back(x);
        if (x < -0.9 || x >= 0.89) x = x+double(0.2);
        else x = x+double(0.3);
        if (fabs(x) < 0.0001) x = 0;
    }

    double eta_binning_new[100];
    for (int i = 0; i < eta_range.size(); i++) eta_binning_new[i] = eta_range.at(i);
    int n_eta_new = eta_range.size()-1;

	//getting binning for trk pt, new fine eta
	TAxis* trk_pt_binning = (TAxis*)((TH3*)theFiles[0]->Get("h_eff_Injet_cent0"))->GetYaxis();
	int n_trk_pt_bins = trk_pt_binning->GetNbins();

	TAxis* trk_eta_binning_new = new TAxis(n_eta_new, eta_binning_new);
	int n_trk_eta_bins_new = trk_eta_binning_new->GetNbins();


	output_file->cd();
//	trk_pt_binning->Write("trk_pt_binning");
    trk_eta_binning_new->Write("new_trk_eta_binning");

	TCanvas *canvas1 = new TCanvas("C1", "C1",0.,0.,800,400);
	TLine *line = new TLine();
	TLatex *ltx = new TLatex();
	ltx->SetNDC();
	ltx->SetTextFont(43);
	ltx->SetTextSize(15);
	ltx->SetTextAlign(22);

	TLegend *legend = new TLegend(0.35,0.15,0.55,0.54,NULL,"brNDC");
	legend->SetBorderSize(0);
	legend->SetNColumns(1);
	legend->SetTextFont(43);
	legend->SetTextSize(17);

	cout << Form("Doing in trk eta slices:") << endl;
	for (int i = 1; i <= n_trk_eta_bins_new; i++) cout << Form("%2.2f, ", trk_eta_binning_new->GetBinLowEdge(i));
	cout << Form("%2.2f", trk_eta_binning_new->GetBinUpEdge(n_trk_eta_bins_new)) << endl;

    vector<vector<TH1*>> h_efficiency(n_cent_cuts, vector<TH1*> (n_trk_eta_bins_new));
	vector<vector<TGraphAsymmErrors*>> g_efficiency(n_cent_cuts, vector<TGraphAsymmErrors*> (n_trk_eta_bins_new));

	cout << Form("Looping of centrality, eta cuts, and files") << endl;
	for (int i_cent_cuts = 0; i_cent_cuts < n_cent_cuts; i_cent_cuts++)
	{
		string centrality = num_to_cent(centrality_scheme,i_cent_cuts);
		cout << centrality << endl;

		for (int i_eta_cuts = 0; i_eta_cuts < n_trk_eta_bins_new; i_eta_cuts++)
		{
			double eta_lo = trk_eta_binning_new->GetBinLowEdge(i_eta_cuts+1);
			double eta_hi = trk_eta_binning_new->GetBinUpEdge(i_eta_cuts+1);

			vector<vector<double> > vec_e(nFiles, vector<double> (n_trk_pt_bins,-1));
			vector<vector<double> > vec_g(nFiles, vector<double> (n_trk_pt_bins,-1));
			vector<vector<double> > vec_lo_e(nFiles, vector<double> (n_trk_pt_bins,-1));
			vector<vector<double> > vec_hi_e(nFiles, vector<double> (n_trk_pt_bins,-1));

			for (int i_files = 0; i_files < nFiles; i_files++)
			{
				name = Form("total_3d_c%i_e%i_f%i", i_cent_cuts, i_eta_cuts, i_files);
				TH3* h_total_3d = (TH3D*)theFiles.at(i_files)->Get(Form("h_trk_foreff_full_cent%i",i_cent_cuts))->Clone(name.c_str());
				h_total_3d->SetName(name.c_str());

				name = Form("entries_3d_c%i_e%i_f%i",i_cent_cuts, i_eta_cuts, i_files);
                TH3* h_entries_3d = (TH3D*)theFiles.at(i_files)->Get(Form("h_trk_foreff_entries_cent%i",i_cent_cuts))->Clone(name.c_str());
				h_entries_3d->SetName(name.c_str());

				name = Form("matched_3d_c%i_e%i_f%i", i_cent_cuts, i_eta_cuts, i_files);
				TH3* h_matched_3d = (TH3D*)theFiles.at(i_files)->Get(Form("h_trk_foreff_matched_cent%i",i_cent_cuts))->Clone(name.c_str());
				h_matched_3d->SetName(name.c_str());

				int eta_lo_bin = h_total_3d->GetZaxis()->FindBin(eta_lo+0.001);
				int eta_hi_bin = h_total_3d->GetZaxis()->FindBin(eta_hi-0.001); //to make sure I get the bin which has upper edge as eta hi

				h_total_3d->GetZaxis()->SetRange(eta_lo_bin, eta_hi_bin);
				h_matched_3d->GetZaxis()->SetRange(eta_lo_bin, eta_hi_bin);
				h_entries_3d->GetZaxis()->SetRange(eta_lo_bin, eta_hi_bin);

				name = Form("total_1d_c%i_e%i_f%i",i_cent_cuts, i_eta_cuts, i_files);
				TH1* h_total_1d = (TH1*)h_total_3d->Project3D("y")->Clone(name.c_str());
				h_total_1d->Sumw2();
				h_total_1d->SetName(name.c_str());

				name = Form("matched_1d_c%i_e%i_f%i",i_cent_cuts, i_eta_cuts, i_files);
				TH1* h_matched_1d = (TH1*)h_matched_3d->Project3D("y")->Clone(name.c_str());
				h_matched_1d->Sumw2();
				h_matched_1d->SetName(name.c_str());

                name = Form("entries_1d_c%i_e%i_f%i",i_cent_cuts, i_eta_cuts, i_files);
                TH1* h_entries_1d = (TH1*)h_entries_3d->Project3D("y")->Clone(name.c_str());
                h_entries_1d->SetName(name.c_str());

                if (i_files == 0)
				{
					name = Form("h_efficiency_c%i_e%i",i_cent_cuts, i_eta_cuts);
					h_efficiency.at(i_cent_cuts).at(i_eta_cuts) = (TH1*)h_matched_1d->Clone(name.c_str());
					h_efficiency.at(i_cent_cuts).at(i_eta_cuts)->Reset();
				}

				name = Form("ratio_1d_c%i_e%i_f%i",i_cent_cuts, i_eta_cuts, i_files);
				TH1* h_ratio = (TH1*)h_matched_1d->Clone(name.c_str());
				h_ratio->SetName(name.c_str());
				h_ratio->Reset();
				h_ratio->Divide(h_matched_1d,h_total_1d,1,1,"B");

				if (1)
				{
					canvas1->Clear();

					canvas1->Divide(2,1);
					canvas1->cd(1);

					gPad->SetLogx();
					gPad->SetLogy();

					h_total_1d->SetLineColor(kBlack);
					h_matched_1d->SetLineColor(kRed);
					h_entries_1d->SetLineColor(kBlue);

					h_total_1d->GetXaxis()->SetRangeUser(1,200);
					h_matched_1d->GetXaxis()->SetRangeUser(1,200);
					h_entries_1d->GetXaxis()->SetRangeUser(1,200);
					h_ratio->GetXaxis()->SetRangeUser(1,200);

					h_total_1d->Draw("hist");
					h_matched_1d->Draw("same hist");
                    h_entries_1d->Draw("same hist text");

					name = Form("%s, %4.2f < eta < %4.2f, JZ%i", num_to_cent(centrality_scheme, i_cent_cuts).c_str(), eta_lo, eta_hi ,i_files+1);
					ltx->SetTextAlign(11);
					ltx->DrawLatex(0.19,0.955,name.c_str());

					canvas1->cd(2);
					gPad->SetLogx();

					h_ratio->SetLineColor(kRed);
					h_ratio->GetYaxis()->SetRangeUser(0,1.2);

					h_ratio->SetMarkerSize(1);
					h_ratio->Draw("p e");
					name = Form("cent%i_eta%i_file%i", i_cent_cuts,i_eta_cuts,i_files);
					ltx->SetTextAlign(11);
					ltx->DrawLatex(0.19,0.955,name.c_str());

					line->SetLineColor(kRed);
					line->SetLineStyle(3);
					line->DrawLine(1,1,200,1);

					if (i_eta_cuts == 0 && i_cent_cuts == 0 && i_files == 0) name = "(";
                    else if (i_eta_cuts == n_trk_eta_bins_new - 1  && i_cent_cuts == n_cent_cuts - 1 && i_files == nFiles - 1) name = ")";
					else name = "";
					canvas1->Print(Form("output_pdf/raw_eff_cent_trketa_%s_%s.pdf%s", dataset_type.c_str(), tracking_cut.c_str(), name.c_str()),Form("Title: eff_c%i_e%i_f%i", i_cent_cuts,i_eta_cuts,i_files));
                    canvas1->Clear();
				}

				for (int i_bin = 0; i_bin < h_matched_1d->GetNbinsX(); i_bin++)
				{
					double tmp_a = h_matched_1d->GetBinContent(i_bin+1);
					double tmp_b = h_total_1d->GetBinContent(i_bin+1);
                    int entries = h_entries_1d->GetBinContent(i_bin+1);

					if (entries < 50) continue;
//					if (tmp_a > tmp_b) continue;

					vec_e.at(i_files).at(i_bin) = h_ratio->GetBinContent(i_bin+1);
					vec_g.at(i_files).at(i_bin) = h_total_1d->GetBinContent(i_bin+1);
					vec_lo_e.at(i_files).at(i_bin) = h_ratio->GetBinError(i_bin+1);
					vec_hi_e.at(i_files).at(i_bin) = h_ratio->GetBinError(i_bin+1);
				}

				h_matched_3d->Reset();
				h_matched_1d->Reset();

				h_total_3d->Reset();
				h_total_1d->Reset();

				h_ratio->Reset();
			}

			for (int i_bin = 0; i_bin < h_efficiency.at(i_cent_cuts).at(i_eta_cuts)->GetNbinsX(); i_bin++)
			{
				double eff_numerator = 0;
				double eff_denominator = 0;

				double err_lo_numerator = 0;
				double err_hi_numerator = 0;

				double tmp = 0;

				bool no_files_used = true;
				for (int i_files = 0; i_files < nFiles; i_files++)
				{
					double e = vec_e.at(i_files).at(i_bin);
					double g = vec_g.at(i_files).at(i_bin);
					double lo_e = vec_lo_e.at(i_files).at(i_bin);
					double hi_e = vec_hi_e.at(i_files).at(i_bin);
					double weight = w.at(i_files);

					if (g <= 0) continue;

					tmp = weight * g * e;
					eff_numerator = eff_numerator + tmp;

					tmp = weight * g;
					eff_denominator = eff_denominator + tmp;

					tmp = pow(weight * g * lo_e, 2);
					err_lo_numerator = err_lo_numerator + tmp;

					tmp = pow(weight * g * hi_e, 2);
					err_hi_numerator = err_hi_numerator + tmp;

					no_files_used = false;
				}

				if (no_files_used) continue;

				double efficiency = eff_numerator/eff_denominator;
				double efficiency_lo = sqrt(err_lo_numerator)/eff_denominator;
				double efficiency_hi = sqrt(err_hi_numerator)/eff_denominator;

				if (efficiency > 1.)
				{
					cout << Form("e>1 e:%f - cent:%i	eta:%i	bin:%f - %f",efficiency, i_cent_cuts, i_eta_cuts, h_efficiency.at(i_cent_cuts).at(i_eta_cuts)->GetXaxis()->GetBinLowEdge(i_bin+1), h_efficiency.at(i_cent_cuts).at(i_eta_cuts)->GetXaxis()->GetBinUpEdge(i_bin+1)) << endl;
				}

				h_efficiency.at(i_cent_cuts).at(i_eta_cuts)->SetBinContent(i_bin+1, efficiency);
				h_efficiency.at(i_cent_cuts).at(i_eta_cuts)->SetBinError(i_bin+1, efficiency_lo);
			}

			vec_e.clear();
			vec_g.clear();
			vec_lo_e.clear();
			vec_hi_e.clear();

			output_file->cd();

			g_efficiency.at(i_cent_cuts).at(i_eta_cuts) = new TGraphAsymmErrors(h_efficiency.at(i_cent_cuts).at(i_eta_cuts));
			name = Form("Efficiency: %s, %4.2f < #eta < %4.2f",centrality.c_str(), eta_lo, eta_hi);
			g_efficiency.at(i_cent_cuts).at(i_eta_cuts)->SetTitle(name.c_str());

			name = Form("graph_eff_eta%i_cent%i",i_eta_cuts,i_cent_cuts);
			g_efficiency.at(i_cent_cuts).at(i_eta_cuts)->SetName(name.c_str());
			g_efficiency.at(i_cent_cuts).at(i_eta_cuts)->Write(name.c_str());

			h_efficiency.at(i_cent_cuts).at(i_eta_cuts)->SetTitle(Form("Efficiency: %s, %4.2f < #eta < %4.2f",centrality.c_str(), eta_lo, eta_hi));
			name = Form("hist_eff_eta%i_cent%i", i_eta_cuts, i_cent_cuts);
			h_efficiency.at(i_cent_cuts).at(i_eta_cuts)->SetName(name.c_str());
			h_efficiency.at(i_cent_cuts).at(i_eta_cuts)->Write(name.c_str());
		}
	}


	output_file->Close();
}
