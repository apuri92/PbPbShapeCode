#include "extras/global_variables.h"

void combine_eff_trketa(vector<TFile*>& input_files, double w[], string cut)
{
	cout << "*********Begin EFFICIENCY*********" << endl;

	string name;

	gStyle->SetOptTitle(0);

	TH3D* h_total_3d;
	TH3D* h_entries_3d;
	TH3D* h_matched_3d;

	TH1* h_total_1d;
	TH1* h_entries_1d;
	TH1* h_matched_1d;

	TH1* h_ratio;
	TGraphAsymmErrors* g_ratio;

	TH1* h_tmp;
	h_tmp = (TH1*)input_files[0]->Get("h_reco_trk_map");
	TAxis* axis_trk_pt = (TAxis*) h_tmp->GetXaxis();
	TAxis* axis_trk_eta = (TAxis*) h_tmp->GetYaxis();
	int n_trk_pt_bins = axis_trk_pt->GetNbins();
	int n_jet_pt_bins = axis_trk_pt->GetNbins();


	vector<vector<vector<<TH1*>>> h_efficiency(n_cent_cuts, vector<vector<TH1*>> (n_jet_pt_bins, vector<TH1*> (n_eta_sl));
	vector<vector<TGraphAsymmErrors*>> g_efficiency(n_cent_cuts, vector<TGraphAsymmErrors*> (eta_range.size() - 1));

	name = Form("mc_raw_eff_ptinclusive_%s.root", cut.c_str());
	TFile *output_file = new TFile(name.c_str(),"recreate");

	TCanvas *canvas1 = new TCanvas("C1", "C1",0.,0.,800,600);
	TCanvas *canvas2 = new TCanvas("C2", "C2",0.,0.,800,1400);
	TLine *line = new TLine();
	TLatex *ltx = new TLatex();
	ltx->SetNDC();
	ltx->SetTextFont(43);
	ltx->SetTextSize(18);
	TLegend *legend = new TLegend(0.15,0.15,0.45,0.35,NULL,"brNDC");
	legend->SetBorderSize(0);
	legend->SetNColumns(1);
	legend->SetTextFont(43);
	legend->SetTextSize(17);


	cout << Form("Looping of centrality, eta cuts, and files") << endl;
	bool first_pass = true;
	for (int i_cent_cuts = 0; i_cent_cuts < n_cent_cuts; i_cent_cuts++)
	{
		string centrality = num_to_cent(centrality_scheme,i_cent_cuts);

		for (int i_eta_cuts = 0; i_eta_cuts < axis_trk_eta_new->GetNbins(); i_eta_cuts++)
		{
//			double eta_lo = eta_range.at(i_eta_cuts);
//			double eta_hi = eta_range.at(i_eta_cuts+1);

			double eta_lo = axis_trk_eta_new->GetBinLowEdge(i_eta_cuts+1);
			double eta_hi = axis_trk_eta_new->GetBinLowEdge(i_eta_cuts+2);

			vector<vector<double> > vec_e(nFiles, vector<double> (n_trk_pt_bins,-1));
			vector<vector<double> > vec_g(nFiles, vector<double> (n_trk_pt_bins,-1));
			vector<vector<double> > vec_lo_e(nFiles, vector<double> (n_trk_pt_bins,-1));
			vector<vector<double> > vec_hi_e(nFiles, vector<double> (n_trk_pt_bins,-1));

			for (int i = 0; i < nFiles; i++)
			{
				for (int j = 0; j < n_trk_pt_bins; j++)
				{
					vec_e.at(i).at(j) = -1;
					vec_g.at(i).at(j) = -1;
					vec_lo_e.at(i).at(j) = -1;
					vec_hi_e.at(i).at(j) = -1;
				}

			}
			for (int i_files = 0; i_files < nFiles; i_files++)
			{
				name = Form("total_3d_c%i_e%i_f%i",i_cent_cuts, i_eta_cuts, i_files);
//				h_total_3d = (TH3D*)input_files.at(i_files)->Get(Form("h_eff_total_cent%i",i_cent_cuts))->Clone(name.c_str());
				h_total_3d = (TH3D*)input_files.at(i_files)->Get(Form("h_trk_foreff_full_cent%i",i_cent_cuts))->Clone(name.c_str());
				h_total_3d->SetName(name.c_str());

				name = Form("entries_3d_c%i_e%i_f%i",i_cent_cuts, i_eta_cuts, i_files);
				h_entries_3d = (TH3D*)input_files.at(i_files)->Get(Form("h_eff_entries_cent%i",i_cent_cuts))->Clone(name.c_str());
				h_entries_3d->SetName(name.c_str());

				name = Form("matched_3d_c%i_e%i_f%i",i_cent_cuts, i_eta_cuts, i_files);
//				h_matched_3d = (TH3D*)input_files.at(i_files)->Get(Form("h_eff_matched_cent%i",i_cent_cuts))->Clone(name.c_str());
				h_matched_3d = (TH3D*)input_files.at(i_files)->Get(Form("h_trk_foreff_matched_cent%i",i_cent_cuts))->Clone(name.c_str());
				h_matched_3d->SetName(name.c_str());

				int pt_lo_bin, pt_hi_bin, eta_lo_bin, eta_hi_bin;

				eta_lo_bin = h_total_3d->GetZaxis()->FindBin(eta_lo+0.001);
				eta_hi_bin = h_total_3d->GetZaxis()->FindBin(eta_hi-0.001); //to make sure I get the bin which has upper edge as eta hi

//				cout << Form("%f - %f : %i - %i",eta_lo, eta_hi, eta_lo_bin, eta_hi_bin) << endl;


				h_total_3d->GetZaxis()->SetRange(eta_lo_bin, eta_hi_bin);
				h_matched_3d->GetZaxis()->SetRange(eta_lo_bin, eta_hi_bin);
				h_entries_3d->GetZaxis()->SetRange(eta_lo_bin, eta_hi_bin);

				name = Form("total_1d_c%i_e%i_f%i",i_cent_cuts, i_eta_cuts, i_files);
				h_total_1d = (TH1*)h_total_3d->Project3D("y")->Clone(name.c_str());
				h_total_1d->SetName(name.c_str());

				name = Form("matched_1d_c%i_e%i_f%i",i_cent_cuts, i_eta_cuts, i_files);
				h_matched_1d = (TH1*)h_matched_3d->Project3D("y")->Clone(name.c_str());
				h_matched_1d->SetName(name.c_str());

				name = Form("entries_1d_c%i_e%i_f%i",i_cent_cuts, i_eta_cuts, i_files);
				h_entries_1d = (TH1*)h_matched_3d->Project3D("y")->Clone(name.c_str());
				h_entries_1d->SetName(name.c_str());

				if (i_files == 0)
				{
					name = Form("h_efficiency_c%i_e%i",i_cent_cuts, i_eta_cuts);
					h_efficiency.at(i_cent_cuts).at(i_eta_cuts) = (TH1*)h_matched_1d->Clone(name.c_str());
					h_efficiency.at(i_cent_cuts).at(i_eta_cuts)->Reset();
				}

				bool skip_histo = false;
				for (int i_bin=1; i_bin <= h_matched_1d->GetNbinsX(); i_bin++)
				{
					double matched = h_matched_1d->GetBinContent(i_bin);
					double total = h_total_1d->GetBinContent(i_bin);
					if (matched > total)
					{
						skip_histo = true;
						cout << Form("Warning - Matched > Truth Entries cent%i_eta%i_file%i_trkbin%i",i_cent_cuts, i_eta_cuts, i_files, i_bin) << endl;}
				}

				//if (skip_histo) continue;


				name = Form("ratio_1d_c%i_e%i_f%i",i_cent_cuts, i_eta_cuts, i_files);
				h_ratio = (TH1*)h_matched_1d->Clone(name.c_str());
				h_ratio->SetName(name.c_str());
				h_ratio->Reset();

				g_ratio = new TGraphAsymmErrors();
				h_ratio->Divide(h_matched_1d,h_total_1d,1,1,"B");
				g_ratio->Divide(h_matched_1d,h_total_1d,"pois");
//				g_ratio->Divide(h_matched_1d,h_total_1d,"cl=0.683 b(1,1) mode");


				if (1)
				{
					canvas1->Clear();

					canvas1->Divide(1,2);
					canvas1->cd(1);

					gPad->SetLogx();
					gPad->SetLogy();

					h_matched_1d->SetLineColor(kRed);
					h_total_1d->SetLineColor(kBlack);

					h_total_1d->Draw();
					h_matched_1d->Draw("same");

					name = Form("All vs Matched: %s, %4.2f < eta < %4.2f, JZ%i", num_to_cent(centrality_scheme, i_cent_cuts).c_str(), eta_lo, eta_hi ,i_files+1);
					ltx->SetTextAlign(22);
					ltx->DrawLatex(0.5,0.95,name.c_str());

					canvas1->cd(2);
					gPad->SetLogx();

					h_ratio->SetLineColor(kRed);
					h_ratio->GetYaxis()->SetRangeUser(0,1.2);

					h_ratio->Draw("hist");
					g_ratio->Draw("same");
					name = Form("Ratio Matched/All: cent%i_eta%i_file%i", i_cent_cuts,i_eta_cuts,i_files);
					ltx->DrawLatex(0.5,0.95,name.c_str());

					line->SetLineColor(kRed);
					line->SetLineStyle(3);
					line->DrawLine(1E-1,1,500,1);

					if (i_eta_cuts == 0 && i_cent_cuts == 0 && i_files == 0) name = Form("eff_ptinclusive_%s.pdf(",cut.c_str());
					else if (i_eta_cuts == eta_range.size() - 2  && i_cent_cuts == n_cent_cuts - 1 && i_files == nFiles - 1) name = Form("eff_ptinclusive_%s.pdf)",cut.c_str());
					else name = Form("eff_ptinclusive_%s.pdf",cut.c_str());

					canvas1->Print(name.c_str(),Form("Title: eff_c%ie%if%i", i_cent_cuts,i_eta_cuts,i_files));

					canvas1->Clear();
				}

				for (int i_bin = 0; i_bin < h_matched_1d->GetNbinsX(); i_bin++)
				{
					if (h_entries_1d->GetBinContent(i_bin+1) <= 0) continue;

					vec_e.at(i_files).at(i_bin) = h_ratio->GetBinContent(i_bin+1);
					vec_g.at(i_files).at(i_bin) = h_total_1d->GetBinContent(i_bin+1);
					vec_lo_e.at(i_files).at(i_bin) = h_ratio->GetBinError(i_bin+1);
					vec_hi_e.at(i_files).at(i_bin) = h_ratio->GetBinError(i_bin+1);

					//						vec_lo_e.at(i_files).at(i_bin) = g_ratio->GetErrorYlow(pnt_number);
					//						vec_hi_e.at(i_files).at(i_bin) = g_ratio->GetErrorYhigh(pnt_number);
					//						cout << Form("bin: %i | lo: %f | hi: %f", i_bin, g_ratio->GetErrorYlow(pnt_number),g_ratio->GetErrorYhigh(pnt_number) ) << endl;
				}

				delete h_ratio;
				delete g_ratio;
				delete h_matched_1d;
				delete h_total_1d;
				delete h_matched_3d;
				delete h_total_3d;
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

					if (g <= 0) continue;

					tmp = w[i_files] * g * e;
					eff_numerator = eff_numerator + tmp;

					tmp = w[i_files] * g;
					eff_denominator = eff_denominator + tmp;

					tmp = pow(w[i_files] * g * lo_e, 2);
					err_lo_numerator = err_lo_numerator + tmp;

					tmp = pow(w[i_files] * g * hi_e, 2);
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

			g_efficiency.at(i_cent_cuts).at(i_eta_cuts)=new TGraphAsymmErrors(h_efficiency.at(i_cent_cuts).at(i_eta_cuts));




			vec_e.clear();
			vec_g.clear();
			vec_lo_e.clear();
			vec_hi_e.clear();

			output_file->cd();
			g_efficiency.at(i_cent_cuts).at(i_eta_cuts)->SetTitle(Form("Efficiency: %s, %4.2f < #eta < %4.2f",centrality.c_str(), eta_lo, eta_hi));
			g_efficiency.at(i_cent_cuts).at(i_eta_cuts)->SetName(Form("g_efficiency_cent%i_eta%i", i_cent_cuts, i_eta_cuts));
			g_efficiency.at(i_cent_cuts).at(i_eta_cuts)->Write(Form("g_efficiency_cent%i_eta%i", i_cent_cuts, i_eta_cuts));


			h_efficiency.at(i_cent_cuts).at(i_eta_cuts)->SetTitle(Form("Efficiency: %s, %4.2f < #eta < %4.2f",centrality.c_str(), eta_lo, eta_hi));
			h_efficiency.at(i_cent_cuts).at(i_eta_cuts)->SetName(Form("h_efficiency_cent%i_eta%i", i_cent_cuts, i_eta_cuts));
			h_efficiency.at(i_cent_cuts).at(i_eta_cuts)->Write(Form("h_efficiency_cent%i_eta%i", i_cent_cuts, i_eta_cuts));
			
		}
	}

	
	output_file->cd();
	axis_trk_eta_new->Write("new_trk_eta_binning");

	output_file->Close();
	
	cout << "*********End EFFICIENCY*********" << endl;

	
}
