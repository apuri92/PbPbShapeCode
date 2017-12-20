#include "extras/global_variables.h"
static const int n_coarse_y = 4;
double coarse_y[n_coarse_y+1] = {0, 0.3, 0.8, 1.2, 2.1};


void combine_eff_jetpt_jety(vector<TFile*>& theFiles, double w[], string cut = "ppTight")
{
	gErrorIgnoreLevel = 3001;
	gStyle->SetOptTitle(0);
	SetAtlasStyle();

	string name;

	name = Form("mc_efficiency_jetpt_jety_%s.root", cut.c_str());
	TFile *output_file = new TFile(name.c_str(),"recreate");
	cout << Form("Creating output root file: %s", output_file->GetName()) << endl;

	TCanvas *canvas1 = new TCanvas("C1", "C1",0.,0.,800,600);
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

	TAxis* jet_pt_binning = (TAxis*)((TH3*)theFiles[0]->Get("h_eff_Injet_cent0"))->GetXaxis();
	int n_jetpt_cuts = jet_pt_binning->GetNbins();

	TAxis* trk_pt_binning = (TAxis*)((TH3*)theFiles[0]->Get("h_eff_Injet_cent0"))->GetYaxis();
	int n_trk_pt_bins = trk_pt_binning->GetNbins();

	TAxis* jet_y_binning_coarse = new TAxis(n_coarse_y, coarse_y);
	int n_jet_y_bins_coarse = jet_y_binning_coarse->GetNbins();

	output_file->cd();
	jet_pt_binning->Write("jet_pt_binning");
	trk_pt_binning->Write("trk_pt_binning");
	jet_y_binning_coarse->Write("jet_y_binning_coarse");

	vector<vector<vector<TH1*>>> h_final_efficiency_injet(n_jetpt_cuts, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (n_jet_y_bins_coarse)));

	cout << Form("Doing in jet y bins:") << endl;
	for (int i = 1; i <= n_jet_y_bins_coarse; i++)cout << Form("%2.2f - %2.2f", jet_y_binning_coarse->GetBinLowEdge(i), jet_y_binning_coarse->GetBinUpEdge(i)) << endl;
	cout << endl;

	cout << Form("Looping of jet pt, centrality, eta cuts, and files") << endl;
	for (int i_jetpt_cuts = 0; i_jetpt_cuts < n_jetpt_cuts; i_jetpt_cuts++)
	{
		double jet_pT_lo = jet_pt_binning->GetBinLowEdge(i_jetpt_cuts+1);
		double jet_pT_hi = jet_pt_binning->GetBinUpEdge(i_jetpt_cuts+1);

		for (int i_cent_cuts = 0; i_cent_cuts < n_cent_cuts; i_cent_cuts++)
		{
			for (int i_eta_cuts = 0; i_eta_cuts < n_jet_y_bins_coarse; i_eta_cuts++)
			{
				double eta_lo = jet_y_binning_coarse->GetBinLowEdge(i_eta_cuts+1);
				double eta_hi = jet_y_binning_coarse->GetBinUpEdge(i_eta_cuts+1);

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
						h_final_efficiency_injet.at(i_jetpt_cuts).at(i_cent_cuts).at(i_eta_cuts) = (TH1*)h1_eff_injet_tmp->Clone(Form("h_final_eff_injet_c%i_j%i_e%i_f%i",i_cent_cuts, i_jetpt_cuts, i_eta_cuts, i_files));
						h_final_efficiency_injet.at(i_jetpt_cuts).at(i_cent_cuts).at(i_eta_cuts)->Reset();
					}

					TH1* h_efficiency_injet = (TH1*)h1_eff_injet_matched_tmp->Clone(Form("h_efficiency_injet_c%i_j%i_e%i_f%i",i_cent_cuts, i_jetpt_cuts, i_eta_cuts, i_files));
					h_efficiency_injet->Reset();

					h_efficiency_injet->Divide(h1_eff_injet_matched_tmp,h1_eff_injet_tmp,1,1,"B");

					if (0)
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

						if (i_jetpt_cuts == 0 && i_eta_cuts == 0 && i_cent_cuts == 0 && i_files == 0) name = "(";
						else if (i_jetpt_cuts == n_jetpt_cuts -1 && i_eta_cuts == n_jet_y_bins_coarse -1  && i_cent_cuts == n_cent_cuts - 1 && i_files == nFiles - 1) name = ")";
						else name = "";
						canvas1->Print(Form("raw_eff_cent_jetpt_jety_%s.pdf%s",cut.c_str(), name.c_str()),Form("Title: eff_pt%i_c%i_e%i_f%i",i_jetpt_cuts, i_cent_cuts,i_eta_cuts,i_files));
						canvas1->Clear();
					}

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

				for (int i_bin=0; i_bin<h_final_efficiency_injet.at(i_jetpt_cuts).at(i_cent_cuts).at(i_eta_cuts)->GetNbinsX(); i_bin++)
				{
					double eff_numerator = 0;
					double eff_denominator = 0;

					double err_lo_numerator = 0;
					double err_hi_numerator = 0;

					double tmp = 0;

					bool no_files_used = true;
					for (int i_files=0; i_files<nFiles; i_files++)
					{
							double e = vec_e.at(i_files).at(i_bin);
							double g = vec_g.at(i_files).at(i_bin);
							double lo_e = vec_lo_e.at(i_files).at(i_bin);
							double hi_e = vec_hi_e.at(i_files).at(i_bin);
							double weight = w[i_files];

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
						cout << Form("e>1 e:%f - cent:%i	eta:%i	bin:%f - %f",efficiency, i_cent_cuts, i_eta_cuts, h_final_efficiency_injet.at(i_jetpt_cuts).at(i_cent_cuts).at(i_eta_cuts)->GetXaxis()->GetBinLowEdge(i_bin+1), h_final_efficiency_injet.at(i_jetpt_cuts).at(i_cent_cuts).at(i_eta_cuts)->GetXaxis()->GetBinUpEdge(i_bin+1)) << endl;
					}

					h_final_efficiency_injet.at(i_jetpt_cuts).at(i_cent_cuts).at(i_eta_cuts)->SetBinContent(i_bin+1, efficiency);
					h_final_efficiency_injet.at(i_jetpt_cuts).at(i_cent_cuts).at(i_eta_cuts)->SetBinError(i_bin+1, efficiency_lo);
				}
				vec_g.clear();
				vec_e.clear();
				vec_hi_e.clear();
				vec_lo_e.clear();
				vec_entries.clear();


				output_file->cd();
				h_final_efficiency_injet.at(i_jetpt_cuts).at(i_cent_cuts).at(i_eta_cuts)->SetTitle(Form("Efficiency: %s, %4.2f < p_{T} < %4.2f, %4.2f < #eta < %4.2f",num_to_cent(30,i_cent_cuts).c_str(), jet_pT_lo, jet_pT_hi, eta_lo, eta_hi));
				h_final_efficiency_injet.at(i_jetpt_cuts).at(i_cent_cuts).at(i_eta_cuts)->SetName(Form("histo_eff_eta%i_cent%i_pt%i",i_eta_cuts, i_cent_cuts, i_jetpt_cuts+1));
				h_final_efficiency_injet.at(i_jetpt_cuts).at(i_cent_cuts).at(i_eta_cuts)->Write(Form("histo_eff_eta%i_cent%i_pt%i",i_eta_cuts, i_cent_cuts, i_jetpt_cuts+1));


			}
		}
	}

	output_file->Close();
}
