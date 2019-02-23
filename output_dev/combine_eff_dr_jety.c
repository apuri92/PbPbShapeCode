#include "functions/global_variables.h"

void combine_eff_dr_jety(vector<TFile*>& theFiles, vector<double> w)
{
	gErrorIgnoreLevel = 3001;
//	SetAtlasStyle();
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
	name = Form("new_mc_eff_%s_trketa_dr_jetptinc_%s.root",dataset_type.c_str(), tracking_cut.c_str());
	TFile *output_file = new TFile(Form("%s",name.c_str()),"recreate");
	cout << Form("Efficiency root file: %s", output_file->GetName()) << endl;

	//setting up jet y bins
	vector<double> jety_range;
	jety_range.push_back(0.); jety_range.push_back(0.3); jety_range.push_back(0.8); jety_range.push_back(1.2); jety_range.push_back(1.7);

	int jety_BinsN = jety_range.size()-1;
	double jety_Bins[10];
	for (int i = 0; i < jety_BinsN+1; i++)
	{
		jety_Bins[i] = jety_range[i];
		cout << jety_Bins[i] << endl;
	}

	TAxis* jety_binning = new TAxis(jety_BinsN, jety_Bins);
	int n_jety_bins = jety_binning->GetNbins();

	TAxis* trk_pt_binning = (TAxis*)((TH3*)theFiles[0]->Get("h_eff_Injet_cent0"))->GetYaxis();
	int n_trk_pt_bins = trk_pt_binning->GetNbins();
	TAxis* dR_binning = (TAxis*)((TH3*)theFiles[0]->Get("h_trk_foreff_r_matched_cent0"))->GetXaxis();
	int n_dR = dR_binning->GetNbins();

	TAxis* jet_pt_binning = (TAxis*)((TH3*)theFiles[0]->Get("h_eff_Injet_cent0"))->GetXaxis();
	int n_jetpt_bins = jet_pt_binning->GetNbins();


	for (int i = 0; i < jety_binning->GetNbins(); i++)
	{
		cout << jety_binning->GetBinLowEdge(i+1) << " " << jety_binning->GetBinUpEdge(i+1) << endl;
	}

//	double x = -2.5;
//	while (x <= 2.5)
//	{
//		eta_range.push_back(x);
//		if (x < -0.9 || x >= 0.89) x = x+double(0.2);
//		else x = x+double(0.3);
//		if (fabs(x) < 0.0001) x = 0;
//	}
//
//	double eta_binning_new[100];
//	for (int i = 0; i < eta_range.size(); i++) eta_binning_new[i] = eta_range.at(i);
//	int n_eta_new = eta_range.size()-1;
//
//	//getting binning for trk pt, new fine eta
//	TAxis* trk_eta_binning_new = new TAxis(n_eta_new, eta_binning_new);
//	int n_trk_eta_bins_new = trk_eta_binning_new->GetNbins();


	output_file->cd();
	trk_pt_binning->Write("trk_pt_binning");
	dR_binning->Write("dR_binning");

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

	vector<vector<vector<vector<TH1*>>>> h_efficiency(n_cent_cuts, vector<vector<vector<TH1*>>> (n_dR, vector<vector<TH1*>> (n_jetpt_bins, vector<TH1*> (n_jety_bins))));
	vector<vector<vector<vector<TGraphAsymmErrors*>>>> g_efficiency(n_cent_cuts, vector<vector<vector<TGraphAsymmErrors*>>> (n_dR, vector<vector<TGraphAsymmErrors*>> (n_jetpt_bins, vector<TGraphAsymmErrors*> (n_jety_bins))));

	cout << Form("Looping of centrality, jety cuts, jetpt, dr, and files") << endl;


	for (int i_cent_cuts = 0; i_cent_cuts < 6; i_cent_cuts++)
	{
		string centrality = num_to_cent(centrality_scheme,i_cent_cuts);
		cout << centrality << endl;

		for (int i_dR = 0; i_dR < n_dR; i_dR++)
		{
			double dr_lo = dR_binning->GetBinLowEdge(i_dR+1);
			double dr_hi = dR_binning->GetBinUpEdge(i_dR+1);

			if (dr_lo > 0.8) continue;

			for (int i_jetpt_cuts = 0; i_jetpt_cuts < n_jetpt_bins; i_jetpt_cuts++)
			{
				double jetpt_lo = jet_pt_binning->GetBinLowEdge(i_jetpt_cuts+1);
				double jetpt_hi = jet_pt_binning->GetBinUpEdge(i_jetpt_cuts+1);

				if (jetpt_lo < 100 || jetpt_hi > 320) continue;


				for (int i_jety_cuts = 0; i_jety_cuts < n_jety_bins; i_jety_cuts++)
				{
					double jety_lo = jety_binning->GetBinLowEdge(i_jety_cuts+1);
					double jety_hi = jety_binning->GetBinUpEdge(i_jety_cuts+1);


					string full_label = Form("%s, %1.2f < r < %1.2f, %1.2f < jetpt < %1.2f, %1.2f < jety < %1.2f", num_to_cent(centrality_scheme, i_cent_cuts).c_str(), dr_lo, dr_hi, jetpt_lo, jetpt_hi, jety_lo, jety_hi);
					cout << full_label << endl;
					vector<vector<double> > vec_e(nFiles, vector<double> (n_trk_pt_bins,-1));
					vector<vector<double> > vec_g(nFiles, vector<double> (n_trk_pt_bins,-1));
					vector<vector<double> > vec_lo_e(nFiles, vector<double> (n_trk_pt_bins,-1));
					vector<vector<double> > vec_hi_e(nFiles, vector<double> (n_trk_pt_bins,-1));

					for (int i_files = 0; i_files < nFiles; i_files++)
					{
						name = Form("total_3d_c%i_r%i_jetpt%i_jety%i_f%i", i_cent_cuts, i_dR, i_jetpt_cuts, i_jety_cuts, i_files);
						TH3* h_total_3d = (TH3D*)theFiles.at(i_files)->Get(Form("h_eff_dR%i_cent%i",i_dR, i_cent_cuts))->Clone(name.c_str());
						h_total_3d->SetName(name.c_str());

						name = Form("entries_3d_c%i_r%i_jetpt%i_jety%i_f%i", i_cent_cuts, i_dR, i_jetpt_cuts, i_jety_cuts, i_files);
						TH3* h_entries_3d = (TH3D*)theFiles.at(i_files)->Get(Form("h_eff_dR%i_cent%i",i_dR, i_cent_cuts))->Clone(name.c_str());
						h_entries_3d->SetName(name.c_str());

						name = Form("matched_3d_c%i_r%i_jetpt%i_jety%i_f%i", i_cent_cuts, i_dR, i_jetpt_cuts, i_jety_cuts, i_files);
						TH3* h_matched_3d = (TH3D*)theFiles.at(i_files)->Get(Form("h_eff_matched_dR%i_cent%i",i_dR, i_cent_cuts))->Clone(name.c_str());
						h_matched_3d->SetName(name.c_str());

						int jety_lo_bin = h_total_3d->GetZaxis()->FindBin(jety_lo+0.001);
						int jety_hi_bin = h_total_3d->GetZaxis()->FindBin(jety_hi-0.001); //to make sure I get the bin which has upper edge as eta hi

						int jetpt_lo_bin = h_total_3d->GetXaxis()->FindBin(jetpt_lo+0.001);
						int jetpt_hi_bin = h_total_3d->GetXaxis()->FindBin(jetpt_hi-0.001); //to make sure I get the bin which has upper edge as eta hi

						h_total_3d->GetZaxis()->SetRange(jety_lo_bin, jety_hi_bin);
						h_matched_3d->GetZaxis()->SetRange(jety_lo_bin, jety_hi_bin);
						h_entries_3d->GetZaxis()->SetRange(jety_lo_bin, jety_hi_bin);

						h_total_3d->GetXaxis()->SetRange(jetpt_lo_bin, jetpt_hi_bin);
						h_matched_3d->GetXaxis()->SetRange(jetpt_lo_bin, jetpt_hi_bin);
						h_entries_3d->GetXaxis()->SetRange(jetpt_lo_bin, jetpt_hi_bin);

						name = Form("total_1d_c%i_r%i_jetpt%i_jety%i_f%i", i_cent_cuts, i_dR, i_jetpt_cuts, i_jety_cuts, i_files);
						TH1* h_total_1d = (TH1*)h_total_3d->Project3D("y")->Clone(name.c_str());
						h_total_1d->Sumw2();
						h_total_1d->SetName(name.c_str());

						name = Form("matched_1d_c%i_r%i_jetpt%i_jety%i_f%i", i_cent_cuts, i_dR, i_jetpt_cuts, i_jety_cuts, i_files);
						TH1* h_matched_1d = (TH1*)h_matched_3d->Project3D("y")->Clone(name.c_str());
						h_matched_1d->Sumw2();
						h_matched_1d->SetName(name.c_str());

						name = Form("entries_1d_c%i_r%i_jetpt%i_jety%i_f%i", i_cent_cuts, i_dR, i_jetpt_cuts, i_jety_cuts, i_files);
						TH1* h_entries_1d = (TH1*)h_entries_3d->Project3D("y")->Clone(name.c_str());
						h_entries_1d->SetName(name.c_str());

						if (i_files == 0)
						{
							name = Form("h_efficiency_c%i_r%i_jetpt%i_jety%i_f%i", i_cent_cuts, i_dR, i_jetpt_cuts, i_jety_cuts, i_files);
							h_efficiency[i_cent_cuts][i_dR][i_jetpt_cuts][i_jety_cuts] = (TH1*)h_matched_1d->Clone(name.c_str());
							h_efficiency[i_cent_cuts][i_dR][i_jetpt_cuts][i_jety_cuts]->Reset();
						}

						name = Form("ratio_1d_c%i_r%i_jetpt%i_jety%i_f%i", i_cent_cuts, i_dR, i_jetpt_cuts, i_jety_cuts, i_files);
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

							name = Form("%s JZ%i", full_label.c_str(), i_files+1);
							ltx->SetTextAlign(11);
							ltx->DrawLatex(0.19,0.955,name.c_str());

							canvas1->cd(2);
							gPad->SetLogx();

							h_ratio->SetLineColor(kRed);
							h_ratio->GetYaxis()->SetRangeUser(0,1.2);

							h_ratio->SetMarkerSize(1);
							h_ratio->Draw("p e");
							name = Form("c%i_r%i_jetpt%i_jety%i_f%i", i_cent_cuts, i_dR, i_jetpt_cuts, i_jety_cuts, i_files);
							ltx->SetTextAlign(11);
							ltx->DrawLatex(0.19,0.955,name.c_str());

							line->SetLineColor(kRed);
							line->SetLineStyle(3);
							line->DrawLine(1,1,200,1);

							if (i_jety_cuts == 0 && i_jetpt_cuts == 0 && i_cent_cuts == 0 && i_dR == 0 && i_files == 0) name = "(";
							else if (i_jety_cuts == n_jety_bins - 1 && i_jetpt_cuts == n_jetpt_bins - 1 && i_cent_cuts == n_cent_cuts - 1 && i_dR == n_dR-1, i_files == nFiles - 1) name = ")";
							else name = "";
							canvas1->Print(Form("output_pdf/new_raw_eff_cent_trketa_dr_%s_%s.pdf%s", dataset_type.c_str(), tracking_cut.c_str(), name.c_str()),Form("Title: eff_c%i_e%i_dr%i_f%i", i_cent_cuts,i_jety_cuts, i_dR, i_files));
							canvas1->Clear();
						}

						for (int i_bin = 0; i_bin < h_matched_1d->GetNbinsX(); i_bin++)
						{
							double tmp_a = h_matched_1d->GetBinContent(i_bin+1);
							double tmp_b = h_total_1d->GetBinContent(i_bin+1);
							int entries = h_entries_1d->GetBinContent(i_bin+1);

							if (entries < 5) continue;

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

					for (int i_bin = 0; i_bin < h_efficiency[i_cent_cuts][i_dR][i_jetpt_cuts][i_jety_cuts]->GetNbinsX(); i_bin++)
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
							cout << Form("e: %1.4f	%s",efficiency, full_label.c_str()) << endl;
						}

						h_efficiency[i_cent_cuts][i_dR][i_jetpt_cuts][i_jety_cuts]->SetBinContent(i_bin+1, efficiency);
						h_efficiency[i_cent_cuts][i_dR][i_jetpt_cuts][i_jety_cuts]->SetBinError(i_bin+1, efficiency_lo);
					}

					vec_e.clear();
					vec_g.clear();
					vec_lo_e.clear();
					vec_hi_e.clear();

					vec_e.shrink_to_fit();
					vec_g.shrink_to_fit();
					vec_lo_e.shrink_to_fit();
					vec_hi_e.shrink_to_fit();

					output_file->cd();

					g_efficiency[i_cent_cuts][i_dR][i_jetpt_cuts][i_jety_cuts] = new TGraphAsymmErrors(h_efficiency[i_cent_cuts][i_dR][i_jetpt_cuts][i_jety_cuts]);
					g_efficiency[i_cent_cuts][i_dR][i_jetpt_cuts][i_jety_cuts]->SetTitle(full_label.c_str());

					name = Form("graph_eff_cent%i_dr%i_jety%i_jetpt%i",i_cent_cuts, i_dR, i_jety_cuts, i_jetpt_cuts);
					g_efficiency[i_cent_cuts][i_dR][i_jetpt_cuts][i_jety_cuts]->SetName(name.c_str());
					g_efficiency[i_cent_cuts][i_dR][i_jetpt_cuts][i_jety_cuts]->Write(name.c_str());

					h_efficiency[i_cent_cuts][i_dR][i_jetpt_cuts][i_jety_cuts]->SetTitle(full_label.c_str());
					name = Form("hist_eff_cent%i_dr%i_jety%i_jetpt%i",i_cent_cuts, i_dR, i_jety_cuts, i_jetpt_cuts);
					h_efficiency[i_cent_cuts][i_dR][i_jetpt_cuts][i_jety_cuts]->SetName(name.c_str());
					h_efficiency[i_cent_cuts][i_dR][i_jetpt_cuts][i_jety_cuts]->Write(name.c_str());

					delete g_efficiency[i_cent_cuts][i_dR][i_jetpt_cuts][i_jety_cuts];
				}
			}
		}

	}

	output_file->Close();
}
