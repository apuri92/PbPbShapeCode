#include "../functions/global_variables.h"
void draw_b2s()
{
	SetAtlasStyle();
	TFile *f_output = new TFile(Form("output_pdf_nominal/root/final_ChPS_MC_PbPb.root"));
	double r_max_range = 0.8;


	TAxis* dR_binning = (TAxis*)f_output->Get("dR_binning");
	TAxis* jetpT_binning = (TAxis*)f_output->Get("jetpT_binning");
	TAxis* trkpT_binning = (TAxis*)f_output->Get("trkpT_binning");
	int N_dR = dR_binning->GetNbins();
	int N_jetpt = jetpT_binning->GetNbins();
	int N_trkpt = trkpT_binning->GetNbins();



	string name;
	string pdf_label;

	TLine *line = new TLine();
	line->SetLineColor(kBlack);

	TLatex *ltx = new TLatex();
	ltx->SetTextFont(43);
	ltx->SetTextSize(12);

	double trk_pt_lo = 1.;
	double trk_pt_hi = 150.;

	double ratio_lo = 0;
	double ratio_hi = 2;

	int jet_pt_start = 7;
	int jet_pt_end = 11;

	int trk_pt_start = 2;
	int trk_pt_end = trkpT_binning->GetNbins()-2;
	vector<vector<vector<TH1*>>> h_ChPS_ratio_B2S_indR (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));


	{
		cout << "Doing B2S plots (as function of r, for trk pT)" << endl;

		gStyle->SetErrorX(0);
		TCanvas *c_B2S_dR = new TCanvas("c_B2S_dR","c_B2S_dR",900,600);
		TLegend *legend_B2S_dR = new TLegend(0.02, 0.62, 0.48, 0.87, "","brNDC");
		legend_B2S_dR->SetTextFont(43);
		legend_B2S_dR->SetBorderSize(0);
		legend_B2S_dR->SetTextSize(18);

		TLegend *legend_B2S_dR_2 = new TLegend(0.02, 0.62, 0.48, 0.87, "","brNDC");
		legend_B2S_dR_2->SetTextFont(43);
		legend_B2S_dR_2->SetBorderSize(0);
		legend_B2S_dR_2->SetTextSize(18);

		int jet_itr = 0;
		for (int i_jet = jet_pt_start; i_jet < jet_pt_start+1; i_jet++)
		{
			string jet_label = Form("%1.0f < #it{p}_{T}^{jet} < %1.0f GeV", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));
			c_B2S_dR->cd();
			c_B2S_dR->Clear();
			c_B2S_dR->Divide(3,2, 0.0, 0.0, 3);


			bool first_pass_cent = true;
			for (int i_cent = 0; i_cent < n_cent_cuts; i_cent++)
			{
				string centrality = num_to_cent(31,i_cent);
				if (i_cent == 6) continue;
				int trk_itr = 0;
				double max = 0., min = 999., tmp;

				c_B2S_dR->cd(i_cent+1);
				double x1, x2, y1, y2;
				gPad->GetPadPar(x1, y1, x2, y2);
				if (i_cent < 3) y1+=0.002;
				else y2-=0.002;
				if (i_cent == 0 || i_cent == 3) x2-=0.002;
				else {x1+=0.002; x2-=0.002;}
				gPad->SetPad(x1, y1, x2, y2);
//				cout << Form("%i: x1, y1, x2, y2: %f, %f, %f, %f",i_cent, x1, y1, x2, y2) << endl;
//				gPad->SetFillColor(i_cent+3);

				x1 = gPad->GetLeftMargin();
				x2 = gPad->GetRightMargin();
				y1 = gPad->GetBottomMargin();
				y2 = gPad->GetTopMargin();
//				if (x1 < 0.00001) x1=0.015;
				if (x2 < 0.00001) x2+=0.007;
//				if (y1 < 0.00001) y1=0.015;
				if (y2 < 0.00001) y2+=0.007;

				cout << Form("%i: %f, %f, %f, %f",i_cent, x1, x2, y1, y2) << endl;

//				if (i_cent < 3) y1+=0.005;
//				else y2-=0.005;
//				if (i_cent == 0 || i_cent == 3) x2-=0.005;
//				else {x1+=0.005; x2-=0.005;}

				gPad->SetMargin(x1, x2, y1, y2);


				for (int i_trk = trk_pt_start; i_trk < trk_pt_end; i_trk++)
				{
					if (i_trk < 2 || i_trk > 8) continue;
					name = Form("h_ChPS_ratio_B2S_indR_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);

					h_ChPS_ratio_B2S_indR.at(i_trk).at(i_cent).at(i_jet) = (TH1*)f_output->Get(name.c_str());

					string trk_label = Form("%1.1f < #it{p}_{T}^{ch} < %1.1f GeV", trkpT_binning->GetBinLowEdge(i_trk+1), trkpT_binning->GetBinUpEdge(i_trk+1));
//					cout << trk_label << " " << centrality << " " << jet_label << endl;

					SetHStyle_smallify(h_ChPS_ratio_B2S_indR.at(i_trk).at(i_cent).at(i_jet), trk_itr, 1);
					h_ChPS_ratio_B2S_indR.at(i_trk).at(i_cent).at(i_jet)->GetXaxis()->SetNdivisions(505);
					h_ChPS_ratio_B2S_indR.at(i_trk).at(i_cent).at(i_jet)->GetXaxis()->SetLabelFont(43);
					h_ChPS_ratio_B2S_indR.at(i_trk).at(i_cent).at(i_jet)->GetXaxis()->SetLabelSize(20);
					h_ChPS_ratio_B2S_indR.at(i_trk).at(i_cent).at(i_jet)->GetYaxis()->SetLabelFont(43);
					h_ChPS_ratio_B2S_indR.at(i_trk).at(i_cent).at(i_jet)->GetYaxis()->SetLabelSize(14);

					h_ChPS_ratio_B2S_indR.at(i_trk).at(i_cent).at(i_jet)->GetXaxis()->SetTitleFont(43);
					h_ChPS_ratio_B2S_indR.at(i_trk).at(i_cent).at(i_jet)->GetXaxis()->SetTitleSize(20);
					h_ChPS_ratio_B2S_indR.at(i_trk).at(i_cent).at(i_jet)->GetYaxis()->SetTitleFont(43);
					h_ChPS_ratio_B2S_indR.at(i_trk).at(i_cent).at(i_jet)->GetYaxis()->SetTitleSize(20);

					h_ChPS_ratio_B2S_indR.at(i_trk).at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(0, r_max_range);

					h_ChPS_ratio_B2S_indR.at(i_trk).at(i_cent).at(i_jet)->GetYaxis()->SetTickSize(0.02);
					if (i_cent < 3) h_ChPS_ratio_B2S_indR.at(i_trk).at(i_cent).at(i_jet)->GetYaxis()->SetRangeUser(8E-1,8E2);
					if (i_cent >=3 ) h_ChPS_ratio_B2S_indR.at(i_trk).at(i_cent).at(i_jet)->GetYaxis()->SetRangeUser(8E-1,2E2);

					TGraph *g_tmp = new TGraph(h_ChPS_ratio_B2S_indR.at(i_trk).at(i_cent).at(i_jet));
					g_tmp->GetYaxis()->SetTitle(h_ChPS_ratio_B2S_indR.at(i_trk).at(i_cent).at(i_jet)->GetYaxis()->GetTitle());
					g_tmp->GetXaxis()->SetTitle(h_ChPS_ratio_B2S_indR.at(i_trk).at(i_cent).at(i_jet)->GetXaxis()->GetTitle());


					//x axis
					g_tmp->GetXaxis()->SetTitleFont(43);
					g_tmp->GetXaxis()->SetTitleOffset(2);
					g_tmp->GetXaxis()->SetTitleSize(22);
					g_tmp->GetXaxis()->SetLabelFont(43);
					g_tmp->GetXaxis()->SetLabelSize(20);
					g_tmp->GetXaxis()->SetLimits(-0.02, 0.83);
					g_tmp->GetXaxis()->SetNdivisions(505);
					g_tmp->GetXaxis()->SetTickSize(0.02);

					//y axis
					g_tmp->GetYaxis()->SetTitleFont(43);
					g_tmp->GetYaxis()->SetTitleOffset(2.5);
					g_tmp->GetYaxis()->SetTitleSize(20);
					g_tmp->GetYaxis()->SetLabelFont(43);
					g_tmp->GetYaxis()->SetLabelSize(20);
					g_tmp->GetYaxis()->SetTickSize(0.02);


					if (i_cent < 3)  g_tmp->GetYaxis()->SetRangeUser(8E-1,2E3);
					if (i_cent >=3 ) g_tmp->GetYaxis()->SetRangeUser(8E-1,2E2);

					if (jet_itr == 0 && first_pass_cent && trk_itr < 3) legend_B2S_dR->AddEntry(h_ChPS_ratio_B2S_indR.at(i_trk).at(i_cent).at(i_jet),trk_label.c_str(),"p");
					if (jet_itr == 0 && first_pass_cent && trk_itr >= 3) legend_B2S_dR_2->AddEntry(h_ChPS_ratio_B2S_indR.at(i_trk).at(i_cent).at(i_jet),trk_label.c_str(),"p");

					c_B2S_dR->cd(i_cent+1);
					if (trk_itr == 0) g_tmp->Draw("AP");
					else g_tmp->Draw("same P");

					gPad->SetLogx(0);
					gPad->SetLogy();


					trk_itr++;

				} // end trk loop

				c_B2S_dR->cd(i_cent+1);
				ltx->SetTextAlign(32);
				ltx->SetTextSize(20);
				ltx->DrawLatexNDC(0.96, 0.92, Form("%s", centrality.c_str()));

				first_pass_cent = false;
			} //end cent loop

			c_B2S_dR->cd(1);
			double x_left = 0.20, x_right = 0.93, y = 0.88, y_diff = 0.09;
			ltx->SetTextAlign(11);
			ltx->SetTextSize(18);
			ltx->DrawLatexNDC(x_left, y, "#scale[1.3]{#font[72]{ATLAS} Internal}");
			ltx->DrawLatexNDC(x_left, y=y-y_diff, "Pb+Pb #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}");
			ltx->DrawLatexNDC(x_left, y=y-y_diff, Form("anti-#font[12]{k}_{#font[12]{t}} R=0.4"));
			ltx->DrawLatexNDC(x_left, y=y-y_diff, Form("%s", jet_label.c_str()));

			c_B2S_dR->cd(5);
			legend_B2S_dR->Draw();
			c_B2S_dR->cd(6);
			legend_B2S_dR_2->Draw();

			pdf_label = "";
			if (i_jet == jet_pt_start)
			{
				c_B2S_dR->Print(Form("output_pdf_nominal/PbPb/UE_B2S_single_0.pdf"));
			}

			jet_itr++;

		} //end jet loop

		delete c_B2S_dR;
		delete legend_B2S_dR;

	}
}
