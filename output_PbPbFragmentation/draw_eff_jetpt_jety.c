//#include "combine_eff_dev.c"
#include "extras/global_variables.h"
#include "TVirtualFitter.h"
void draw_eff_jetpt_jety(string cut = "ppTight")
{
	cout << "Draw_eff_smooth.c" << endl;
	string name;

	gStyle->SetOptTitle(0);
	SetAtlasStyle();

	string pt_binning;
	pt_binning = "_pt_exclusive";

	//reading from input file
	name = Form("mc_efficiency_jetpt_jety_%s.root", cut.c_str());
	TFile *input_file = new TFile(name.c_str());

	cout << "********NOT FITTING IF PPTIGHT_TIGHT********" << endl;
	//creatign output file
	name = Form("mc_eff_fits_jetpt_jety_%s.root",cut.c_str());
	TFile *output_file = new TFile(name.c_str(),"recreate");
	TCanvas *canvas1 = new TCanvas("C1", "C1",0.,0.,900,600);
	TCanvas *canvas2 = new TCanvas("C3", "C2",0.,0.,900,600);

	TLine *line = new TLine();
	TLatex *ltx = new TLatex();
	ltx->SetNDC();
	ltx->SetTextFont(43);
	ltx->SetTextSize(15);
	ltx->SetTextAlign(11);
	TLegend *legend = new TLegend(0.10,0.3,0.9,0.70,NULL,"brNDC");
	legend->SetBorderSize(0);
	legend->SetNColumns(1);
	legend->SetTextFont(43);
	legend->SetTextSize(12);

	TLegend *legend1 = new TLegend(0.20,0.59,0.9,0.87,NULL,"brNDC");
	legend1->SetBorderSize(0);
	legend1->SetNColumns(1);
	legend1->SetTextFont(43);
	legend1->SetTextSize(12);

	TLegend *legend2 = new TLegend(0.20,0.59,0.9,0.87,NULL,"brNDC");
	legend2->SetBorderSize(0);
	legend2->SetNColumns(1);
	legend2->SetTextFont(43);
	legend2->SetTextSize(10);


	TAxis* jet_pt_binning = (TAxis*)input_file->Get("jet_pt_binning");
	TAxis* trk_pt_binning = (TAxis*)input_file->Get("trk_pt_binning");
	int n_jetpt_cuts = jet_pt_binning->GetNbins();
	int n_trk_pt_bins = trk_pt_binning->GetNbins();

	TGraphAsymmErrors* g_efficiency[n_cent_cuts][n_eta_cuts][n_jetpt_cuts];
	TH1* h_efficiency[n_cent_cuts][n_eta_cuts][n_jetpt_cuts];
	TH1* h_efficiency_fit[n_cent_cuts][n_eta_cuts][n_jetpt_cuts];
	TH1* h_fit_to_hist[n_cent_cuts][n_eta_cuts][n_jetpt_cuts];

	TF1 *fit_eff[n_cent_cuts][n_eta_cuts][n_jetpt_cuts];
	//	TH1* ConfidenceIntervals[n_cent_cuts][n_eta_cuts][n_jetpt_cuts];

	TH1* h_tmp;
	h_tmp = new TH1D("h_tmp","h_tmp",n_eta_cuts,eta_Slices);
	TAxis *jet_y_binning = (TAxis*)h_tmp->GetXaxis();

	output_file->cd();
	jet_pt_binning->Write("jet_pt_binning");
	jet_y_binning->Write("jet_y_binning");
	h_tmp->Write("h_tmp");

	canvas1->cd();
	legend->AddEntry(jet_pt_binning,"x","lp");
	legend->Draw();
	legend1->AddEntry(jet_pt_binning,"x","lp");
	legend1->Draw();
	legend2->AddEntry(jet_pt_binning,"x","lp");
	legend2->Draw();
	canvas1->Print("tmp.pdf");
	remove("tmp.pdf");

	legend->Clear();
	legend1->Clear();
	legend2->Clear();
	canvas1->Clear();
	bool smoothing = 0;

	for (int i_cent_cuts = 0; i_cent_cuts < n_cent_cuts; i_cent_cuts++)
	{
		string centrality = num_to_cent(centrality_scheme,i_cent_cuts);

		for (int i_eta_cuts = 0; i_eta_cuts < n_eta_cuts; i_eta_cuts++)
		{
			double eta_lo = eta_Slices[i_eta_cuts];
			double eta_hi = eta_Slices[i_eta_cuts+1];

			for (int i_jetpt_cuts = 0; i_jetpt_cuts < n_jetpt_cuts; i_jetpt_cuts++)
			{
				double jet_pT_lo = jet_pt_binning->GetBinLowEdge(i_jetpt_cuts+1);
				double jet_pT_hi = jet_pt_binning->GetBinLowEdge(i_jetpt_cuts+2);

				name = Form("histo_eff_eta%i_cent%i_pt%i",i_eta_cuts, i_cent_cuts, i_jetpt_cuts+1);
				h_efficiency[i_cent_cuts][i_eta_cuts][i_jetpt_cuts] = (TH1*)input_file->Get(name.c_str());

				//getting low edge of first bin with z = 1, this is the high fit range (bin with z = 1 is not included if low edge is used). all tracks up to the high fit range have z < 1
				double fit_range_lo = 9.98;
				double fit_range_hi =  h_efficiency[i_cent_cuts][i_eta_cuts][i_jetpt_cuts]->GetXaxis()->GetBinUpEdge(n_trk_pt_bins);
				int first_z_1_bin = n_trk_pt_bins;
				for (int i_bin=1; i_bin<=n_trk_pt_bins; i_bin++)
				{
					double trk_pT_hi = h_efficiency[i_cent_cuts][i_eta_cuts][i_jetpt_cuts]->GetXaxis()->GetBinUpEdge(i_bin);
					double trk_pT_lo = h_efficiency[i_cent_cuts][i_eta_cuts][i_jetpt_cuts]->GetXaxis()->GetBinLowEdge(i_bin);
					if (trk_pT_hi > jet_pT_lo && cut != "ppTight_tight")
					{
						//						cout << Form("%2.4f - %2.4f : %2.4f - %2.4f", jet_pT_lo, jet_pT_hi, trk_pT_lo, trk_pT_hi) << endl;

						fit_range_hi = h_efficiency[i_cent_cuts][i_eta_cuts][i_jetpt_cuts]->GetBinLowEdge(i_bin);
						first_z_1_bin = i_bin;
						break;
					}
					else if (trk_pT_lo > jet_pT_hi && cut == "ppTight_tight")
					{
						fit_range_hi = h_efficiency[i_cent_cuts][i_eta_cuts][i_jetpt_cuts]->GetBinLowEdge(i_bin);
						first_z_1_bin = i_bin;
						break;
					}

				}



				name = Form("fit_eff_eta%i_cent%i_pt%i",i_eta_cuts, i_cent_cuts, i_jetpt_cuts+1);
				fit_eff[i_cent_cuts][i_eta_cuts][i_jetpt_cuts] = new TF1(name.c_str(),"[0] + [1]*log(x) + [2]*pow(log(x),2) + [3]*pow(log(x),3)");
				fit_eff[i_cent_cuts][i_eta_cuts][i_jetpt_cuts]->SetRange(fit_range_lo,fit_range_hi);


				//fit within special range
				h_efficiency[i_cent_cuts][i_eta_cuts][i_jetpt_cuts]->Fit(fit_eff[i_cent_cuts][i_eta_cuts][i_jetpt_cuts],"RQ0","");


				//making histogram with fit. Using original points below low fit range (this is the start_change_bin), otherwise using fit->Eval() up to high fit range (last bin with z < 1). For the remaining bins use eff in previous bin. Relative errors are set in the new histogram. For drawing purposes (see loops below) the bins with trk_pt_hi > jet_pt_low are suppressed.

				name = Form("h_eff_fit_eff_eta%i_cent%i_pt%i",i_eta_cuts, i_cent_cuts, i_jetpt_cuts+1);
				h_efficiency_fit[i_cent_cuts][i_eta_cuts][i_jetpt_cuts] = (TH1*)h_efficiency[i_cent_cuts][i_eta_cuts][i_jetpt_cuts]->Clone(name.c_str());
				int start_change_bin = h_efficiency_fit[i_cent_cuts][i_eta_cuts][i_jetpt_cuts]->FindBin(fit_range_lo);
				for (int i_bin=start_change_bin; i_bin<=n_trk_pt_bins; i_bin++)
				{

					float e_original = h_efficiency[i_cent_cuts][i_eta_cuts][i_jetpt_cuts]->GetBinContent(i_bin);
					float e_err_original = h_efficiency[i_cent_cuts][i_eta_cuts][i_jetpt_cuts]->GetBinError(i_bin);

					if (e_original == 0) continue;

					float e_fit = fit_eff[i_cent_cuts][i_eta_cuts][i_jetpt_cuts]->Eval(h_efficiency_fit[i_cent_cuts][i_eta_cuts][i_jetpt_cuts]->GetBinCenter(i_bin));
					if (cut != "ppTight_tight")
					{
						h_efficiency_fit[i_cent_cuts][i_eta_cuts][i_jetpt_cuts]->SetBinContent(i_bin, e_fit);
						h_efficiency_fit[i_cent_cuts][i_eta_cuts][i_jetpt_cuts]->SetBinError(i_bin, e_err_original * e_fit/e_original);
					}
					if (i_bin >= first_z_1_bin)
					{
						double e_prev = h_efficiency_fit[i_cent_cuts][i_eta_cuts][i_jetpt_cuts]->GetBinContent(i_bin-1);
						double e_prev_error = h_efficiency_fit[i_cent_cuts][i_eta_cuts][i_jetpt_cuts]->GetBinError(i_bin-1);

						if (cut != "ppTight_tight")
						{
							h_efficiency_fit[i_cent_cuts][i_eta_cuts][i_jetpt_cuts]->SetBinContent(i_bin, e_prev);
							h_efficiency_fit[i_cent_cuts][i_eta_cuts][i_jetpt_cuts]->SetBinError(i_bin, e_err_original * e_prev/e_original);
						}

						if (cut == "ppTight_tight")
						{
							h_efficiency_fit[i_cent_cuts][i_eta_cuts][i_jetpt_cuts]->SetBinContent(i_bin, 0);
							h_efficiency_fit[i_cent_cuts][i_eta_cuts][i_jetpt_cuts]->SetBinError(i_bin, 0);
						}

					}



				}



				//comparing histogram from fit to original points
				name = Form("h_fit_to_hist_eta%i_cent%i_pt%i",i_eta_cuts, i_cent_cuts, i_jetpt_cuts+1);
				h_fit_to_hist[i_cent_cuts][i_eta_cuts][i_jetpt_cuts] = (TH1*)h_efficiency_fit[i_cent_cuts][i_eta_cuts][i_jetpt_cuts]->Clone(name.c_str());
				h_fit_to_hist[i_cent_cuts][i_eta_cuts][i_jetpt_cuts]->Divide(h_efficiency[i_cent_cuts][i_eta_cuts][i_jetpt_cuts]);

				output_file->cd();
				name = Form("histo_eff_eta%i_cent%i_pt%i",i_eta_cuts,i_cent_cuts,i_jetpt_cuts+1);
				h_efficiency_fit[i_cent_cuts][i_eta_cuts][i_jetpt_cuts]->SetName(name.c_str());
				h_efficiency_fit[i_cent_cuts][i_eta_cuts][i_jetpt_cuts]->Write(name.c_str());

				//Making final efficiency graph from new histogram from fit values
				g_efficiency[i_cent_cuts][i_eta_cuts][i_jetpt_cuts] = new TGraphAsymmErrors(h_efficiency_fit[i_cent_cuts][i_eta_cuts][i_jetpt_cuts]);
				name = Form("graph_eff_eta%i_cent%i_pt%i",i_eta_cuts,i_cent_cuts,i_jetpt_cuts+1);
				g_efficiency[i_cent_cuts][i_eta_cuts][i_jetpt_cuts]->Write(name.c_str());
			}
		}
	}


	for (int i_eta_cuts = 0; i_eta_cuts < n_eta_cuts; i_eta_cuts++)
	{
		double eta_lo = eta_Slices[i_eta_cuts];
		double eta_hi = eta_Slices[i_eta_cuts+1];

		canvas1->cd();
		canvas1->Clear();
		canvas1->Divide(3,2);

		for (int i_cent_cuts = 0; i_cent_cuts < n_cent_cuts ; i_cent_cuts++)
		{
			canvas1->cd(i_cent_cuts+1);
			gPad->SetLogx();
			gPad->SetLogy(0);

			canvas2->cd();
			canvas2->Clear();
			canvas2->Divide(3,2);


			string centrality = num_to_cent(centrality_scheme,i_cent_cuts);

			int start = 7;
			int end = 13;
			int canvas_number = 1;
			for (int i_jetpt_cuts = start; i_jetpt_cuts < end ; i_jetpt_cuts++)
			{
				double jet_pT_lo = jet_pt_binning->GetBinLowEdge(i_jetpt_cuts+1);
				double jet_pT_hi = jet_pt_binning->GetBinLowEdge(i_jetpt_cuts+2);

				for (int i_bin=1; i_bin<=n_trk_pt_bins; i_bin++)
				{
					double trk_pT_hi = h_efficiency[i_cent_cuts][i_eta_cuts][i_jetpt_cuts]->GetXaxis()->GetBinUpEdge(i_bin);
					double trk_pT_lo = h_efficiency[i_cent_cuts][i_eta_cuts][i_jetpt_cuts]->GetXaxis()->GetBinLowEdge(i_bin);

					if (trk_pT_lo >= jet_pT_hi) h_efficiency[i_cent_cuts][i_eta_cuts][i_jetpt_cuts]->SetBinContent(i_bin,-1);
				}


				SetHStyle(h_efficiency[i_cent_cuts][i_eta_cuts][i_jetpt_cuts], i_jetpt_cuts-start+1);
				smallify(h_efficiency[i_cent_cuts][i_eta_cuts][i_jetpt_cuts]);

				SetHStyle(h_efficiency_fit[i_cent_cuts][i_eta_cuts][i_jetpt_cuts], i_jetpt_cuts-start+1);
				smallify(h_efficiency_fit[i_cent_cuts][i_eta_cuts][i_jetpt_cuts]);

				h_efficiency_fit[i_cent_cuts][i_eta_cuts][i_jetpt_cuts]->GetYaxis()->SetTitle("Efficiency");
				h_efficiency_fit[i_cent_cuts][i_eta_cuts][i_jetpt_cuts]->GetXaxis()->SetTitle("p_{T}^{trk} [GeV]");

//				h_efficiency[i_cent_cuts][i_eta_cuts][i_jetpt_cuts]->GetXaxis()->SetLabelSize(12);
//				h_efficiency[i_cent_cuts][i_eta_cuts][i_jetpt_cuts]->GetXaxis()->SetTitleOffset(2.5);
//				h_efficiency[i_cent_cuts][i_eta_cuts][i_jetpt_cuts]->GetXaxis()->SetTitleSize(14);
//				h_efficiency[i_cent_cuts][i_eta_cuts][i_jetpt_cuts]->GetYaxis()->SetLabelSize(12);
//				h_efficiency[i_cent_cuts][i_eta_cuts][i_jetpt_cuts]->GetYaxis()->SetTitleOffset(3.);
//				h_efficiency[i_cent_cuts][i_eta_cuts][i_jetpt_cuts]->GetYaxis()->SetTitleSize(14);
//				h_efficiency[i_cent_cuts][i_eta_cuts][i_jetpt_cuts]->SetMarkerSize(0.9);

				h_efficiency_fit[i_cent_cuts][i_eta_cuts][i_jetpt_cuts]->GetYaxis()->SetRangeUser(0.2,1.);
				h_efficiency_fit[i_cent_cuts][i_eta_cuts][i_jetpt_cuts]->GetXaxis()->SetRangeUser(1.,1e3);
				h_efficiency[i_cent_cuts][i_eta_cuts][i_jetpt_cuts]->GetXaxis()->SetRangeUser(0.5,1e3);

				canvas1->cd(i_cent_cuts+1);
				if (i_jetpt_cuts == start) h_efficiency_fit[i_cent_cuts][i_eta_cuts][i_jetpt_cuts]->Draw("p");
				else h_efficiency_fit[i_cent_cuts][i_eta_cuts][i_jetpt_cuts]->Draw("same p");

				name = Form("%2.0f < p_{T}^{Jet} < %2.0f GeV",jet_pT_lo, jet_pT_hi);

				if (i_eta_cuts == 0 && i_cent_cuts == 0)
				{
					legend1->AddEntry(h_efficiency[i_cent_cuts][i_eta_cuts][i_jetpt_cuts],name.c_str(),"lp");
				}



				canvas2->cd(canvas_number);
				gPad->SetLogx();
				gPad->SetLogy(0);

				SetHStyle(h_fit_to_hist[i_cent_cuts][i_eta_cuts][i_jetpt_cuts], i_jetpt_cuts-start+1);
				smallify(h_fit_to_hist[i_cent_cuts][i_eta_cuts][i_jetpt_cuts]);
				h_fit_to_hist[i_cent_cuts][i_eta_cuts][i_jetpt_cuts]->GetYaxis()->SetRangeUser(0.8,1.2);
				h_fit_to_hist[i_cent_cuts][i_eta_cuts][i_jetpt_cuts]->GetXaxis()->SetRangeUser(1.,1e3);
				h_fit_to_hist[i_cent_cuts][i_eta_cuts][i_jetpt_cuts]->Draw("p");

				line->SetLineColor(kRed);
				line->SetLineStyle(3);
				line->DrawLine(1,1,500,1);

				name = Form("%s, %1.1f < y^{jet} < %1.1f",centrality.c_str(), eta_lo, eta_hi);
				ltx->SetTextAlign(11);
				ltx->DrawLatex(0.19,0.88,name.c_str());

				name = Form("#bf{%1.0f < p_{T}^{jet} < %1.0f}", jet_pT_lo, jet_pT_hi);
				ltx->SetTextAlign(11);
				ltx->DrawLatex(0.19,0.82,name.c_str());


				canvas_number++;

			}

			canvas2->cd(2);
			legend1->SetX1NDC(0.21);
			legend1->SetY1NDC(0.19);
			legend1->SetX2NDC(0.7);
			legend1->SetY2NDC(0.53);
			legend1->Draw();

			canvas1->cd(i_cent_cuts+1);
			name = Form("%s",centrality.c_str());
			ltx->SetTextAlign(11);
			ltx->DrawLatex(0.19,0.88,name.c_str());


			if (i_eta_cuts == 0 && i_cent_cuts == 0) name = "(";
			else if (i_eta_cuts == n_eta_cuts - 1 && i_cent_cuts == n_cent_cuts - 1) name = ")";
			else name = "";
			canvas2->Print(Form("eff_fitQual_cent_jetpt_jety_%s.pdf%s", cut.c_str(), name.c_str()), Form("Title: cent%i_jety%i",i_cent_cuts, i_eta_cuts));
//			canvas2->Print(Form("eff_fitQual_cent%i_jetpt_jety%i_%s.pdf", i_cent_cuts, i_eta_cuts, cut.c_str()));

		}

		canvas1->cd(1);
		PlotLabels_AtlasSim_q2_5(0.19, 0.455, 18, 11, true);

		ltx->SetTextAlign(11);
		name = Form("%4.2f < y^{jet} < %4.2f", eta_lo, eta_hi);
		ltx->DrawLatex(0.19,0.270,name.c_str());

		name = Form("%s", cut.c_str());
		if (cut == "ppTight_tight") name = "Tight";
		if (cut == "ppTight") name = "Default";
		if (cut == "manual") name = "Default (no N^{sh}_{Hits} cut)";
		ltx->DrawLatex(0.19,0.210,name.c_str());

		canvas1->cd(2);
		legend1->SetX1NDC(0.21);
		legend1->SetY1NDC(0.19);
		legend1->SetX2NDC(0.7);
		legend1->SetY2NDC(0.53);
		legend1->SetTextSize(12);
		legend1->Draw();


		if (i_eta_cuts == 0) name = Form("(");
		else if (i_eta_cuts == n_eta_cuts - 1) name = Form(")");
		else name = Form("");
		canvas1->Print(Form("eff_cent_jetpt_jety_%s.pdf%s",cut.c_str(), name.c_str()),Form("Title: jety%i",i_eta_cuts));
//		canvas1->Print(Form("eff_centrality_jetpt_jety%i_%s.pdf",i_eta_cuts, cut.c_str()));
		canvas1->Clear();
	}
	
	output_file->Close();
}


