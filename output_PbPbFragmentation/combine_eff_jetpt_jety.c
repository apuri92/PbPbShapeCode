#include "extras/global_variables.h"

void combine_eff_jetpt_jety(vector<TFile*>& theFiles, double w[], string cut = "pptight")
{
	gErrorIgnoreLevel = 3001;
	cout << "*********Begin EFFICIENCY*********" << endl;

	bool draw_option = 1;
	double jet_pt_cut = 100;

	TCanvas *canvas1 = new TCanvas("C1", "C1",0.,0.,800,600);
	TCanvas *canvas2 = new TCanvas("C2", "C2",0.,0.,800,1400);
	TLine *line = new TLine();
	TLatex *ltx = new TLatex();
	ltx->SetNDC();
	ltx->SetTextFont(43);
	ltx->SetTextSize(12);
	ltx->SetTextAlign(22);

	TLegend *legend = new TLegend(0.35,0.15,0.55,0.54,NULL,"brNDC");
	legend->SetBorderSize(0);
	legend->SetNColumns(1);
	legend->SetTextFont(43);
	legend->SetTextSize(17);

	string name;
	double x;
	double y;

	cout << Form("Creating output root files") << endl;
	name = Form("mc_efficiency_jetpt_jety_%s.root", cut.c_str());
	TFile *f_eff = new TFile(name.c_str(),"recreate");

	TAxis* jet_pt_binning = (TAxis*)((TH3*)theFiles[0]->Get("h_eff_Injet_cent0"))->GetXaxis();
	int n_jetpt_cuts = jet_pt_binning->GetNbins();

	TAxis* trk_pt_binning = (TAxis*)((TH3*)theFiles[0]->Get("h_eff_Injet_cent0"))->GetYaxis();
	int n_trk_pt_bins = trk_pt_binning->GetNbins();

	f_eff->cd();
	jet_pt_binning->Write("jet_pt_binning");
	trk_pt_binning->Write("trk_pt_binning");

	cout << Form("Creating histograms") << endl;
	TH1 *h_final_efficiency_injet[n_jetpt_cuts][n_cent_cuts][n_eta_cuts];
	TGraphAsymmErrors* g_efficiency_injet[n_jetpt_cuts][n_cent_cuts][n_eta_cuts];


	cout << Form("Doing in eta slices:") << endl;
	for (int i = 0; i < n_eta_cuts; i++) cout << Form("%2.2f - %2.2f", eta_Slices[i], eta_Slices[i+1]) << endl;

	cout << Form("Looping of centrality, eta cuts, and files") << endl;

	for (int i_jetpt_cuts = 0; i_jetpt_cuts < n_jetpt_cuts; i_jetpt_cuts++)
	{
		double jet_pT_lo = jet_pt_binning->GetBinLowEdge(i_jetpt_cuts+1);
		double jet_pT_hi = jet_pt_binning->GetBinUpEdge(i_jetpt_cuts+1);

		for (int i_cent_cuts = 0; i_cent_cuts < n_cent_cuts; i_cent_cuts++)
		{
			for (int i_eta_cuts = 0; i_eta_cuts < n_eta_cuts; i_eta_cuts++)
			{
				double eta_lo = eta_Slices[i_eta_cuts];
				double eta_hi = eta_Slices[i_eta_cuts+1];

				vector<vector<double> > vec_e(nFiles, vector<double> (n_trk_pt_bins,-1));
				vector<vector<double> > vec_g(nFiles, vector<double> (n_trk_pt_bins,-1));
				vector<vector<double> > vec_lo_e(nFiles, vector<double> (n_trk_pt_bins,-1));
				vector<vector<double> > vec_hi_e(nFiles, vector<double> (n_trk_pt_bins,-1));
				vector<vector<int> > vec_entries(nFiles, vector<int> (n_trk_pt_bins,-1));

				for (int i_files = 0; i_files < nFiles; i_files++)
				{
					name = Form("h_eff_injet_tmp_c%i_j%i_e%i_f%i",i_cent_cuts, i_jetpt_cuts, i_eta_cuts, i_files);
					TH3* h_eff_injet_tmp = (TH3*)theFiles[i_files]->Get(Form("h_eff_Injet_cent%i",i_cent_cuts))->Clone(name.c_str());
					h_eff_injet_tmp->SetName(name.c_str());

					name = Form("h_eff_injet_matched_tmp_c%i_j%i_e%i_f%i",i_cent_cuts, i_jetpt_cuts, i_eta_cuts, i_files);
					TH3* h_eff_injet_matched_tmp = (TH3*)theFiles[i_files]->Get(Form("h_eff_Injet_matched_cent%i",i_cent_cuts))->Clone(name.c_str());
					h_eff_injet_matched_tmp->SetName(name.c_str());

					name = Form("h_eff_injet_entries_tmp_c%i_j%i_e%i_f%i",i_cent_cuts, i_jetpt_cuts, i_eta_cuts, i_files);
					TH3* h_eff_injet_entries_tmp = (TH3*)theFiles[i_files]->Get(Form("h_eff_Injet_entries_cent%i",i_cent_cuts))->Clone(name.c_str());
					h_eff_injet_entries_tmp->SetName(name.c_str());

					int eta_lo_bin = h_eff_injet_tmp->GetZaxis()->FindBin(eta_lo+0.0001);
					int eta_hi_bin = h_eff_injet_tmp->GetZaxis()->FindBin(eta_hi-0.0001); //to make sure I get the bin which has upper edge as eta hi

					h_eff_injet_tmp->GetZaxis()->SetRange(eta_lo_bin, eta_hi_bin);
					h_eff_injet_tmp->GetXaxis()->SetRange(i_jetpt_cuts+1, i_jetpt_cuts+1); //cut on jet pT
					name = Form("h1_eff_injet_tmp_c%i_j%i_e%i_f%i",i_cent_cuts, i_jetpt_cuts, i_eta_cuts, i_files);
					TH1* h1_eff_injet_tmp = (TH1*)h_eff_injet_tmp->Project3D("y")->Clone(name.c_str());
					h1_eff_injet_tmp->SetName(name.c_str());

					h_eff_injet_matched_tmp->GetZaxis()->SetRange(eta_lo_bin, eta_hi_bin);
					h_eff_injet_matched_tmp->GetXaxis()->SetRange(i_jetpt_cuts+1, i_jetpt_cuts+1); //cut on jet pT
					name = Form("h1_eff_injet_matched_tmp_c%i_j%i_e%i_f%i",i_cent_cuts, i_jetpt_cuts, i_eta_cuts, i_files);
					TH1* h1_eff_injet_matched_tmp = (TH1*)h_eff_injet_matched_tmp->Project3D("y")->Clone(name.c_str());
					h1_eff_injet_matched_tmp->SetName(name.c_str());

					h_eff_injet_entries_tmp->GetZaxis()->SetRange(eta_lo_bin, eta_hi_bin);
					h_eff_injet_entries_tmp->GetXaxis()->SetRange(i_jetpt_cuts+1, i_jetpt_cuts+1); //cut on jet pT
					name = Form("h1_eff_injet_entries_tmp_c%i_j%i_e%i_f%i",i_cent_cuts, i_jetpt_cuts, i_eta_cuts, i_files);
					TH1* h1_eff_injet_entries_tmp = (TH1*)h_eff_injet_entries_tmp->Project3D("y")->Clone(name.c_str());
					h1_eff_injet_entries_tmp->SetName(name.c_str());

					//used later to get binning
					if (i_files==0)
					{
						h_final_efficiency_injet[i_jetpt_cuts][i_cent_cuts][i_eta_cuts] = (TH1*)h1_eff_injet_tmp->Clone(Form("h_final_eff_injet_c%i_j%i_e%i_f%i",i_cent_cuts, i_jetpt_cuts, i_eta_cuts, i_files));
						h_final_efficiency_injet[i_jetpt_cuts][i_cent_cuts][i_eta_cuts]->Reset();
					}

					TH1* h_efficiency_injet = (TH1*)h1_eff_injet_matched_tmp->Clone(Form("h_efficiency_injet_c%i_j%i_e%i_f%i",i_cent_cuts, i_jetpt_cuts, i_eta_cuts, i_files));
					h_efficiency_injet->Reset();

					h_efficiency_injet->Divide(h1_eff_injet_matched_tmp,h1_eff_injet_tmp,1,1,"B");

					if (1)
					{
						canvas1->Clear();

						canvas1->Divide(1,2);
						canvas1->cd(1);

						gPad->SetLogx();
						gPad->SetLogy();

						h1_eff_injet_matched_tmp->SetLineColor(kRed);
						h1_eff_injet_tmp->SetLineColor(kBlack);

						double x = 1;h1_eff_injet_entries_tmp->GetMinimum();
						double y = h1_eff_injet_tmp->GetMaximum() * 1E1;
						h1_eff_injet_tmp->GetYaxis()->SetRangeUser(x,y);
						h1_eff_injet_tmp->SetTitle(Form("All vs Matched: jetpt%i_cent%i_eta%i_file%i",i_jetpt_cuts, i_cent_cuts,i_eta_cuts,i_files));
						h1_eff_injet_tmp->Draw();
						h1_eff_injet_matched_tmp->Draw("same");
						h1_eff_injet_entries_tmp->Draw("same hist text");

						canvas1->cd(2);
						gPad->SetLogx();

						h_efficiency_injet->SetLineColor(kRed);
						h_efficiency_injet->GetYaxis()->SetRangeUser(0,1.2);
						h_efficiency_injet->SetTitle(Form("Ratio Matched/All: jetpt%i_cent%i_eta%i_file%i",i_jetpt_cuts, i_cent_cuts,i_eta_cuts,i_files));

						h_efficiency_injet->Draw("hist text");

						line->SetLineColor(kRed);
						line->SetLineStyle(3);
						line->DrawLine(1E-1,1,500,1);
						ltx->DrawLatex(0.5,0.8,Form("Cent: %s, %2.2f < pT < %2.2f, %2.2f < y < %2.2f, file: %i", num_to_cent(30, i_cent_cuts).c_str(), jet_pT_lo, jet_pT_hi, eta_lo, eta_hi, i_files ));

						if (i_jetpt_cuts == 0 && i_eta_cuts == 0 && i_cent_cuts == 0 && i_files == 0) name = Form("eff.pdf(");
						else if (i_jetpt_cuts == n_jetpt_cuts -1 && i_eta_cuts == n_eta_cuts -1  && i_cent_cuts == n_cent_cuts - 1 && i_files == nFiles - 1) name = Form("eff.pdf)");
						else name = Form("eff.pdf");

						canvas1->Print(name.c_str(),Form("Title: eff_pt%ic%ie%if%i",i_jetpt_cuts, i_cent_cuts,i_eta_cuts,i_files));

						canvas1->Clear();
					}

					int pnt_number = 0;
					double x_tmp = 0;
					double y_tmp = 0;


					for (int i_bin=0; i_bin<h_efficiency_injet->GetNbinsX(); i_bin++)
					{
						double tmp_a = h1_eff_injet_matched_tmp->GetBinContent(i_bin+1);
						double tmp_b = h1_eff_injet_tmp->GetBinContent(i_bin+1);
						int entries = h1_eff_injet_entries_tmp->GetBinContent(i_bin+1);

						if (tmp_a > tmp_b) continue;
						if (entries < 50) continue;

						vec_entries.at(i_files).at(i_bin) = h1_eff_injet_entries_tmp->GetBinContent(i_bin+1);
						vec_e.at(i_files).at(i_bin) = h_efficiency_injet->GetBinContent(i_bin+1);
						vec_g.at(i_files).at(i_bin) = h1_eff_injet_tmp->GetBinContent(i_bin+1);
						vec_lo_e.at(i_files).at(i_bin) = h_efficiency_injet->GetBinError(i_bin+1);
						vec_hi_e.at(i_files).at(i_bin) = h_efficiency_injet->GetBinError(i_bin+1);
					}

					h_eff_injet_tmp->Reset();
					h_eff_injet_matched_tmp->Reset();;

					h1_eff_injet_tmp->Reset();;
					h1_eff_injet_matched_tmp->Reset();;

					h_efficiency_injet->Reset();
				}

				std::vector<double> px;
				std::vector<double> py;
				std::vector<double> d_l_py;
				std::vector<double> d_h_py;
				std::vector<double> d_l_px;
				std::vector<double> d_h_px;

				for (int i_bin=0; i_bin<h_final_efficiency_injet[0][0][0]->GetNbinsX(); i_bin++)
				{
					double num = 0;
					double den = 0;
					double rat = 0;

					double num_injet = 0;
					double den_injet = 0;
					double rat_injet = 0;

					double d_l_num_injet = 0;
					double d_l_rat_injet = 0;
					double d_h_num_injet = 0;
					double d_h_rat_injet = 0;

					for (int i_files=0; i_files<nFiles; i_files++)
					{
						if (vec_g.at(i_files).at(i_bin) < 0) continue;
						double num_tmp_injet = w[i_files]*vec_g.at(i_files).at(i_bin)*vec_e.at(i_files).at(i_bin);
						num_injet = num_injet + num_tmp_injet;

						double d_l_num_tmp_injet = pow(w[i_files]*vec_g.at(i_files).at(i_bin)*vec_lo_e.at(i_files).at(i_bin),2);
						d_l_num_injet = d_l_num_injet + d_l_num_tmp_injet;

						double d_h_num_tmp_injet = pow(w[i_files]*vec_g.at(i_files).at(i_bin)*vec_hi_e.at(i_files).at(i_bin),2);
						d_h_num_injet = d_h_num_injet + d_h_num_tmp_injet;

						double den_tmp_injet = w[i_files]*vec_g.at(i_files).at(i_bin);
						den_injet = den_injet + den_tmp_injet;
					}

					if (den_injet==0) continue;
					rat_injet = num_injet/den_injet;

					d_l_rat_injet = sqrt(d_l_num_injet/(den_injet*den_injet));
					d_h_rat_injet = sqrt(d_h_num_injet/(den_injet*den_injet));

					px.push_back(h_final_efficiency_injet[i_jetpt_cuts][i_cent_cuts][i_eta_cuts]->GetBinCenter(i_bin+1));
					py.push_back(rat_injet);
					d_l_py.push_back(d_l_rat_injet);
					d_h_py.push_back(d_h_rat_injet);

					h_final_efficiency_injet[i_jetpt_cuts][i_cent_cuts][i_eta_cuts]->SetBinContent(i_bin+1, rat_injet);
					h_final_efficiency_injet[i_jetpt_cuts][i_cent_cuts][i_eta_cuts]->SetBinError(i_bin+1, d_l_rat_injet);
				}


				f_eff->cd();
				h_final_efficiency_injet[i_jetpt_cuts][i_cent_cuts][i_eta_cuts]->SetTitle(Form("Efficiency: %s, %4.2f < p_{T} < %4.2f, %4.2f < #eta < %4.2f",num_to_cent(30,i_cent_cuts).c_str(), jet_pT_lo, jet_pT_hi, eta_lo, eta_hi));

				h_final_efficiency_injet[i_jetpt_cuts][i_cent_cuts][i_eta_cuts]->Write(Form("histo_eff_eta%i_cent%i_pt%i",i_eta_cuts, i_cent_cuts, i_jetpt_cuts+1));

				vec_g.clear();
				vec_e.clear();
				vec_hi_e.clear();
				vec_lo_e.clear();
				vec_entries.clear();
			}
		}
	}

	cout << Form("Done! Saving and closing") << endl;
	f_eff->Close();

	cout << "*********End EFFICIENCY*********" << endl;


}

