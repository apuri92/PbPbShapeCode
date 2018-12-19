#include "output_dev/functions/global_variables.h"

void compare_final_UE()
{
	SetAtlasStyle();
	TFile *file_MC_method = new TFile("final_ChPS_data_PbPb_c21_UEMCmethod.root");
	TFile *file_cone_method = new TFile("final_ChPS_data_PbPb_c21_UEconemethod.root");

	TAxis* dR_binning = (TAxis*)file_MC_method->Get("dR_binning");
	TAxis* jetpT_binning = (TAxis*)file_MC_method->Get("jetpT_binning");
	TAxis* trkpT_binning = (TAxis*)file_MC_method->Get("trkpT_binning");

	double r_max_range = 0.8;

	{
		cout << "Comparing Final ChPS plots (as function of r, for trk pT)" << endl;

		TCanvas *c_ChPS_dR = new TCanvas("c_ChPS_dR","c_ChPS_dR",900,600);
		TLegend *legend_ChPS_dR = new TLegend(0.19, 0.700, 0.40, 0.93, "","brNDC");
		legend_ChPS_dR->SetTextFont(43);
		legend_ChPS_dR->SetBorderSize(0);
		legend_ChPS_dR->SetTextSize(13);
		int jet_pt_start = 7;
		int jet_pt_end = 11;

		int trk_pt_start = 1;
		int trk_pt_end = 9;

		int jet_itr = 0;

		TLine *line = new TLine();
		line->SetLineColor(kBlack);

		TLatex *ltx = new TLatex();
		ltx->SetTextFont(43);
		ltx->SetTextSize(12);
		string name;

		for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
		{
			string jet_label = Form("%1.0f < #it{p}_{T}^{jet} < %1.0f GeV", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));
			c_ChPS_dR->cd();
			c_ChPS_dR->Clear();
			c_ChPS_dR->Divide(3,2);

			bool first_pass_cent = true;
			for (int i_cent = 0; i_cent < 6; i_cent++)
			{
				string centrality = num_to_cent(31,i_cent);

				int trk_itr = 0;
				double max = 0., min = 999., tmp;
				for (int i_trk = 2; i_trk < 7; i_trk++)
				{
					string trk_label = Form("%1.1f < #it{p}_{T} < %1.1f GeV", trkpT_binning->GetBinLowEdge(i_trk+1), trkpT_binning->GetBinUpEdge(i_trk+1));

//					name = Form("h_ChPS_UE_indR_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
//					name = Form("h_ChPS_raw_indR_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
					name = Form("h_ChPS_raw_subtr_indR_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
//					name = Form("h_ChPS_raw_subtr_unf_indR_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
//					name = Form("h_ChPS_final_indR_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
					TH1* h_MC_meth = (TH1*)file_MC_method->Get(name.c_str());
					h_MC_meth->SetName(Form("%s_MCMeth", name.c_str()));

					TH1* h_cone_meth = (TH1*)file_cone_method->Get(name.c_str());
					h_cone_meth->SetName(Form("%s_ConeMeth", name.c_str()));

					TH1* h_ratio = (TH1*)h_cone_meth->Clone(Form("%s_ratio", name.c_str()));
					h_ratio->Divide(h_MC_meth);
					h_ratio->GetYaxis()->SetTitle("ChPS_{Cone Method} / ChPS_{MC Method}");
					SetHStyle_smallify(h_ratio, trk_itr, 1);

					for (int i = 0; i < h_ratio->GetXaxis()->GetNbins(); i++)
					{
						h_ratio->SetBinError(i+1,0.);
					}
					if (jet_itr == 0 && first_pass_cent) legend_ChPS_dR->AddEntry(h_ratio,trk_label.c_str(),"lp");

					c_ChPS_dR->cd(i_cent+1);

					h_ratio->GetYaxis()->SetRangeUser(0,2);
					h_ratio->GetXaxis()->SetRangeUser(0,r_max_range);

					if (trk_itr == 0) h_ratio->DrawCopy("p");
					else h_ratio->DrawCopy("p same");
					gPad->SetLogx(0);
					gPad->SetLogy(0);

					line->DrawLine(0, 1, r_max_range, 1);


					trk_itr++;

				} // end trk loop


				c_ChPS_dR->cd(i_cent+1);
				ltx->SetTextAlign(12);
				ltx->SetTextSize(14);
				ltx->DrawLatexNDC(0.19, 0.21, Form("%s", jet_label.c_str()));
//				ltx->SetTextAlign(32);
				ltx->DrawLatexNDC(0.19, 0.27, Form("%s", centrality.c_str()));

				first_pass_cent = false;
			} //end cent loop



			c_ChPS_dR->cd(1);
			ltx->SetTextAlign(32);
			ltx->SetTextSize(12);
			ATLASLabel(0.59, 0.83, "     Internal", "", kBlack);
			ltx->SetTextAlign(32);
			ltx->SetTextSize(15);
			ltx->DrawLatexNDC(0.93, 0.90, "Pb+Pb #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}");

			c_ChPS_dR->cd(2);
			legend_ChPS_dR->Draw();

			string pdf_label = "";
			if (i_jet == jet_pt_start) pdf_label = "(";
			if (i_jet == jet_pt_end-1) pdf_label = ")";
			c_ChPS_dR->Print(Form("ChPS_UE_Comparison.pdf%s", pdf_label.c_str()), Form("Title:jetpt%i", i_jet));
//			if (i_jet == jet_pt_start) c_ChPS_dR->Print("tmp.root");

			jet_itr++;

		} //end jet loop
	}








}
