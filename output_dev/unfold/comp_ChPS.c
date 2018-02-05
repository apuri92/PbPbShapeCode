#include "../functions/global_variables.h"
#include "TEnv.h"

void comp_ChPS()
{

	SetAtlasStyle();
	gErrorIgnoreLevel = 3001;

	cout << "Drawing ChPS..." << endl;

	TFile *f_data_PbPb = new TFile("final_ChPS_data_PbPb.root");
	TFile *f_data_pp = new TFile("final_ChPS_data_pp.root");

	TAxis* dR_binning = (TAxis*)f_data_PbPb->Get("dR_binning");
	TAxis* jetpT_binning = (TAxis*)f_data_PbPb->Get("jetpT_binning");
	TAxis* trkpT_binning = (TAxis*)f_data_PbPb->Get("trkpT_binning");

	int N_dR = dR_binning->GetNbins();
	int N_jetpt = jetpT_binning->GetNbins();
	int N_trkpt = trkpT_binning->GetNbins();

	double array_dr_bins[N_dR+1];
	for (int i_dR = 0; i_dR <= N_dR; i_dR++) array_dr_bins[i_dR] = dR_binning->GetBinLowEdge(i_dR+1);

	vector<vector<vector<TH1*>>> h_ChPS_final_data_PbPb (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_final_data_pp (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_final_data_PbPb_dR (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_final_data_pp_dR (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));

	vector<vector<vector<TH1*>>> h_ChPS_ratio_data_PbPb_data_pp (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_ratio_data_PbPb_data_pp_dR (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));

	vector<vector<vector<TH1*>>> h_ChPS_noUnf_data_PbPb (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_noUnf_data_pp (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_noUnf_data_PbPb_dR (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_noUnf_data_pp_dR (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));

	vector<vector<vector<TH1*>>> h_ChPS_noUnf_ratio_data_PbPb_data_pp (N_dR, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_noUnf_ratio_data_PbPb_data_pp_dR (N_trkpt, vector<vector<TH1*>> (n_cent_cuts, vector<TH1*> (N_jetpt)));

	string name;
	string pdf_label;

	TLine *line = new TLine();
	line->SetLineColor(kBlack);

	TLatex *ltx = new TLatex();
	ltx->SetTextFont(43);
	ltx->SetTextSize(12);

	double trk_pt_lo = 1.;
	double trk_pt_hi = 150.;
	double ratio_lo = 0.;
	double ratio_hi = 5.;

	int jet_pt_start = 7;
	int jet_pt_end = 11;

	for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
	{
		for (int i_cent = 0; i_cent < 6; i_cent++)
		{


			for (int i_trk = 0; i_trk < N_trkpt; i_trk++)
			{
				if (i_cent == 0) //get the pp plots from cent6
				{
					name = Form("h_ChPS_final_dR_trk%i_cent%i_jetpt%i", i_trk, 6, i_jet);
					h_ChPS_final_data_pp_dR.at(i_trk).at(6).at(i_jet) = (TH1*)f_data_pp->Get(name.c_str())->Clone(Form("pp_data_%s",name.c_str()));

					name = Form("h_ChPS_noUnf_dR_trk%i_cent%i_jetpt%i", i_trk, 6, i_jet);
					h_ChPS_noUnf_data_pp_dR.at(i_trk).at(6).at(i_jet) = (TH1*)f_data_pp->Get(name.c_str())->Clone(Form("pp_data_%s",name.c_str()));
				}

				//final
				name = Form("h_ChPS_final_dR_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_ChPS_final_data_PbPb_dR.at(i_trk).at(i_cent).at(i_jet) = (TH1*)f_data_PbPb->Get(name.c_str())->Clone(Form("PbPb_data_%s",name.c_str()));

				name = Form("h_ChPS_final_ratio_data_PbPb_data_pp_dR_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_ChPS_ratio_data_PbPb_data_pp_dR.at(i_trk).at(i_cent).at(i_jet) = (TH1*)h_ChPS_final_data_PbPb_dR.at(i_trk).at(i_cent).at(i_jet)->Clone(name.c_str());
				h_ChPS_ratio_data_PbPb_data_pp_dR.at(i_trk).at(i_cent).at(i_jet)->Divide(h_ChPS_final_data_pp_dR.at(i_trk).at(6).at(i_jet));


				//not unfolded
				name = Form("h_ChPS_noUnf_dR_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_ChPS_noUnf_data_PbPb_dR.at(i_trk).at(i_cent).at(i_jet) = (TH1*)f_data_PbPb->Get(name.c_str())->Clone(Form("PbPb_data_%s",name.c_str()));

				name = Form("h_ChPS_noUnf_ratio_data_PbPb_data_pp_dR_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_ChPS_noUnf_ratio_data_PbPb_data_pp_dR.at(i_trk).at(i_cent).at(i_jet) = (TH1*)h_ChPS_noUnf_data_PbPb_dR.at(i_trk).at(i_cent).at(i_jet)->Clone(name.c_str());
				h_ChPS_noUnf_ratio_data_PbPb_data_pp_dR.at(i_trk).at(i_cent).at(i_jet)->Divide(h_ChPS_noUnf_data_pp_dR.at(i_trk).at(6).at(i_jet));

			}


			for (int i_dR = 0; i_dR < N_dR; i_dR++)
			{

				if (i_cent == 0) //get the pp plots from cent6
				{
					name = Form("h_ChPS_final_dR%i_cent%i_jetpt%i", i_dR, 6, i_jet);
					h_ChPS_final_data_pp.at(i_dR).at(6).at(i_jet) = (TH1*)f_data_pp->Get(name.c_str())->Clone(Form("pp_data_%s",name.c_str()));

					name = Form("h_ChPS_noUnf_dR%i_cent%i_jetpt%i", i_dR, 6, i_jet);
					h_ChPS_noUnf_data_pp.at(i_dR).at(6).at(i_jet) = (TH1*)f_data_pp->Get(name.c_str())->Clone(Form("pp_data_%s",name.c_str()));
				}

				//final
				name = Form("h_ChPS_final_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet);
				h_ChPS_final_data_PbPb.at(i_dR).at(i_cent).at(i_jet) = (TH1*)f_data_PbPb->Get(name.c_str())->Clone(Form("PbPb_data_%s",name.c_str()));

				name = Form("h_ChPS_final_ratio_data_PbPb_data_pp_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet);
				h_ChPS_ratio_data_PbPb_data_pp.at(i_dR).at(i_cent).at(i_jet) = (TH1*)h_ChPS_final_data_PbPb.at(i_dR).at(i_cent).at(i_jet)->Clone(name.c_str());
				h_ChPS_ratio_data_PbPb_data_pp.at(i_dR).at(i_cent).at(i_jet)->Divide(h_ChPS_final_data_pp.at(i_dR).at(6).at(i_jet));


				//not unfolded
				name = Form("h_ChPS_noUnf_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet);
				h_ChPS_noUnf_data_PbPb.at(i_dR).at(i_cent).at(i_jet) = (TH1*)f_data_PbPb->Get(name.c_str())->Clone(Form("PbPb_data_%s",name.c_str()));

				name = Form("h_ChPS_noUnf_ratio_data_PbPb_data_pp_dR%i_cent%i_jetpt%i", i_dR, i_cent, i_jet);
				h_ChPS_noUnf_ratio_data_PbPb_data_pp.at(i_dR).at(i_cent).at(i_jet) = (TH1*)h_ChPS_noUnf_data_PbPb.at(i_dR).at(i_cent).at(i_jet)->Clone(name.c_str());
				h_ChPS_noUnf_ratio_data_PbPb_data_pp.at(i_dR).at(i_cent).at(i_jet)->Divide(h_ChPS_noUnf_data_pp.at(i_dR).at(6).at(i_jet));
			}
		}
	}


	//final ratio
	{
		cout << "Doing Final ChPS plots in jet pT" << endl;
		
		TCanvas *c_ChPS_final_ratio = new TCanvas("c_ChPS_final_ratio","c_ChPS_final_ratio",900,600);
		TLegend *legend_ChPS_final_ratio = new TLegend(0.20, 0.43, 0.40, 0.60, "","brNDC");
		legend_ChPS_final_ratio->SetTextFont(43);
		legend_ChPS_final_ratio->SetBorderSize(0);
		legend_ChPS_final_ratio->SetTextSize(8);
		
		for (int i_dR = 0; i_dR < N_dR; i_dR++)
		{
			string dr_label = Form("%1.2f < dR < %1.2f", dR_binning->GetBinLowEdge(i_dR+1), dR_binning->GetBinUpEdge(i_dR+1));

			c_ChPS_final_ratio->cd();
			c_ChPS_final_ratio->Clear();
			c_ChPS_final_ratio->Divide(3,2);
			
			bool first_pass_cent = true;
			for (int i_cent = 0; i_cent < 6; i_cent++)
			{
				string centrality = num_to_cent(31,i_cent);
				
				int jet_itr = 0;
				for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
				{
					string jet_label = Form("%1.0f < p_{T}^{Jet} < %1.0f", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));

					SetHStyle_smallify(h_ChPS_ratio_data_PbPb_data_pp.at(i_dR).at(i_cent).at(i_jet), jet_itr, 1);
					if (i_dR == 0 && first_pass_cent) legend_ChPS_final_ratio->AddEntry(h_ChPS_ratio_data_PbPb_data_pp.at(i_dR).at(i_cent).at(i_jet),jet_label.c_str(),"lp");

					h_ChPS_ratio_data_PbPb_data_pp.at(i_dR).at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(trk_pt_lo, trk_pt_hi);
					h_ChPS_ratio_data_PbPb_data_pp.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetRangeUser(ratio_lo, ratio_hi);
					h_ChPS_ratio_data_PbPb_data_pp.at(i_dR).at(i_cent).at(i_jet)->GetYaxis()->SetNdivisions(504);

					c_ChPS_final_ratio->cd(i_cent+1);
					if (i_jet == jet_pt_start) h_ChPS_ratio_data_PbPb_data_pp.at(i_dR).at(i_cent).at(i_jet)->Draw("");
					else h_ChPS_ratio_data_PbPb_data_pp.at(i_dR).at(i_cent).at(i_jet)->Draw("same");

					gPad->SetLogx();
					gPad->SetLogy(0);

					jet_itr++;
				} // end jet loop
				
				c_ChPS_final_ratio->cd(i_cent+1);
				ltx->SetTextAlign(32);
				ltx->SetTextSize(12);
				ltx->DrawLatexNDC(0.93, 0.90, Form("%s", dr_label.c_str()));
				ltx->DrawLatexNDC(0.93, 0.85, Form("%s", centrality.c_str()));
				line->DrawLine(trk_pt_lo, 1, trk_pt_hi, 1);

				legend_ChPS_final_ratio->Draw();
				first_pass_cent = false;
			} //end cent loop
			
			pdf_label = "";
			if (i_dR == 0) pdf_label = "(";
			if (i_dR == N_dR-1) pdf_label = ")";
			c_ChPS_final_ratio->Print(Form("ChPS_final_ratio.pdf%s", pdf_label.c_str()), Form("Title:dR%i", i_dR));
			
		} //end dR loop
	}


	{
		cout << "Doing Final ChPS ratio (PbPb/pp) plots in dR" << endl;

		TCanvas *c_ChPS_dR = new TCanvas("c_ChPS_dR","c_ChPS_dR",900,600);
		TLegend *legend_ChPS_dR = new TLegend(0.19, 0.80, 0.40, 0.92, "","brNDC");
		legend_ChPS_dR->SetTextFont(43);
		legend_ChPS_dR->SetBorderSize(0);
		legend_ChPS_dR->SetTextSize(10);

		int jet_itr = 0;
		for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
		{
			string jet_label = Form("%1.0f < p_{T}^{Jet} < %1.0f", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));
			c_ChPS_dR->cd();
			c_ChPS_dR->Clear();
			c_ChPS_dR->Divide(3,2);

			bool first_pass_cent = true;
			for (int i_cent = 0; i_cent < 6; i_cent++)
			{
				string centrality = num_to_cent(31,i_cent);

				c_ChPS_dR->cd(i_cent+1);

				int trk_itr = 0;
				for (int i_trk = 0; i_trk < N_trkpt; i_trk++)
				{
					if (i_trk != 2 &&
						i_trk != 5 &&
						i_trk != 7 &&
						i_trk != 8) continue;
					//1, 5, 10, 40

					string trk_label = Form("%1.2f < p_{T}^{Trk} < %1.2f", trkpT_binning->GetBinLowEdge(i_trk+1), trkpT_binning->GetBinUpEdge(i_trk+1));
//					cout << jet_itr << " " << trk_itr << " " << i_cent << endl;
					SetHStyle_smallify(h_ChPS_ratio_data_PbPb_data_pp_dR.at(i_trk).at(i_cent).at(i_jet), trk_itr, 1);
					if (jet_itr == 0 && first_pass_cent) legend_ChPS_dR->AddEntry(h_ChPS_ratio_data_PbPb_data_pp_dR.at(i_trk).at(i_cent).at(i_jet),trk_label.c_str(),"lp");

					h_ChPS_ratio_data_PbPb_data_pp_dR.at(i_trk).at(i_cent).at(i_jet)->Print("all");
//					h_ChPS_ratio_data_PbPb_data_pp_dR.at(i_trk).at(i_cent).at(i_jet)->GetXaxis()->SetRangeUser(0, 1.2);
					h_ChPS_ratio_data_PbPb_data_pp_dR.at(i_trk).at(i_cent).at(i_jet)->GetYaxis()->SetRangeUser(ratio_lo, ratio_hi);
					h_ChPS_ratio_data_PbPb_data_pp_dR.at(i_trk).at(i_cent).at(i_jet)->GetYaxis()->SetNdivisions(504);


					if (trk_itr == 0) h_ChPS_ratio_data_PbPb_data_pp_dR.at(i_trk).at(i_cent).at(i_jet)->Draw("");
					else h_ChPS_ratio_data_PbPb_data_pp_dR.at(i_trk).at(i_cent).at(i_jet)->Draw("same");

					gPad->SetLogx(0);
					gPad->SetLogy(0);

					trk_itr++;

				} // end trk loop

				c_ChPS_dR->cd(i_cent+1);
				ltx->SetTextAlign(32);
				ltx->SetTextSize(12);
				ltx->DrawLatexNDC(0.93, 0.90, Form("%s", jet_label.c_str()));
				ltx->DrawLatexNDC(0.93, 0.85, Form("%s", centrality.c_str()));
				line->DrawLine(0, 1, 1.2, 1);

				legend_ChPS_dR->Draw();
				first_pass_cent = false;
			} //end cent loop

			pdf_label = "";
			if (i_jet == jet_pt_start) pdf_label = "(";
			if (i_jet == jet_pt_end-1) pdf_label = ")";
			c_ChPS_dR->Print(Form("ChPS_final_ratio_dR.pdf%s", pdf_label.c_str()), Form("Title:jetpt%i", i_jet));
			//
			jet_itr++;

		} //end jet loop
	}


}
