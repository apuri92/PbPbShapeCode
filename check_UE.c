#include "output_dev/functions/global_variables.h"

void check_UE_updated()
{
	SetAtlasStyle();
	gErrorIgnoreLevel = 3001;

	TFile *input_file = new TFile("output_dev/raw_results/FF_MC_out_histo_PbPb_5p02_r001.root");
	//	TFile *input_file = new TFile("hist-local_mc.root");
	TAxis* dR_binning = (TAxis*)((TH3*)input_file->Get("h_dR_change_jetpt0_cent0"))->GetXaxis();
	TAxis* jetpT_binning = (TAxis*)((TH3*)input_file->Get("ChPS_raw_0_dR0_cent0"))->GetYaxis();
	TAxis* trkpT_binning = (TAxis*)((TH3*)input_file->Get("ChPS_raw_0_dR0_cent0"))->GetXaxis();


	string name;
	TCanvas *c_TM_MB_2D = new TCanvas("c_TM_MB_2D","c_TM_MB_2D",900,600);
	TCanvas *c_FS_MB_2D = new TCanvas("c_FS_MB_2D","c_FS_MB_2D",900,600);
	TCanvas *c_FNS_MB_2D = new TCanvas("c_FNS_MB_2D","c_FNS_MB_2D",900,600);

	TCanvas *c_TM_MB_1D = new TCanvas("c_TM_MB_1D","c_TM_MB_1D",900,600);
	TCanvas *c_FS_MB_1D = new TCanvas("c_FS_MB_1D","c_FS_MB_1D",900,600);
	TCanvas *c_FNS_MB_1D = new TCanvas("c_FNS_MB_1D","c_FNS_MB_1D",900,600);

	TLatex *ltx = new TLatex();
	ltx->SetTextFont(43);
	ltx->SetTextSize(12);
	ltx->SetTextAlign(12);

	TLegend *legend_TM_MB = new TLegend(0.19, 0.70, 0.40, 0.92, "","brNDC");
	legend_TM_MB->SetTextFont(43);
	legend_TM_MB->SetBorderSize(0);
	legend_TM_MB->SetTextSize(10);

	TLegend *legend_FS_MB = new TLegend(0.19, 0.70, 0.40, 0.92, "","brNDC");
	legend_FS_MB->SetTextFont(43);
	legend_FS_MB->SetBorderSize(0);
	legend_FS_MB->SetTextSize(10);

	TLegend *legend_FNS_MB = new TLegend(0.19, 0.70, 0.40, 0.92, "","brNDC");
	legend_FNS_MB->SetTextFont(43);
	legend_FNS_MB->SetBorderSize(0);
	legend_FNS_MB->SetTextSize(10);

	int jet_pt_start = 7;
	int jet_pt_end = 11;


	int N_dR = dR_binning->GetNbins();
	int N_jetpt = jetpT_binning->GetNbins();
	int N_trkpt = trkpT_binning->GetNbins();

	vector<vector<TH2*>> h_TM_MB_2D = vector<vector<TH2*>> (N_dR, vector<TH2*> (n_cent_cuts));
	vector<vector<TH2*>> h_FS_MB_2D = vector<vector<TH2*>> (N_dR, vector<TH2*> (n_cent_cuts));
	vector<vector<TH2*>> h_FNS_MB_2D = vector<vector<TH2*>> (N_dR, vector<TH2*> (n_cent_cuts));

	vector<vector<vector<TH1*>>> h_TM_MB_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_dR, vector<TH1*> (n_cent_cuts)));
	vector<vector<vector<TH1*>>> h_FS_MB_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_dR, vector<TH1*> (n_cent_cuts)));
	vector<vector<vector<TH1*>>> h_FNS_MB_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_dR, vector<TH1*> (n_cent_cuts)));

	vector<vector<vector<TH1*>>> h_TM_MB_r_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_trkpt, vector<TH1*> (n_cent_cuts)));
	vector<vector<vector<TH1*>>> h_FS_MB_r_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_trkpt, vector<TH1*> (n_cent_cuts)));
	vector<vector<vector<TH1*>>> h_FNS_MB_r_1D =  vector<vector<vector<TH1*>>> (N_jetpt, vector<vector<TH1*>> (N_trkpt, vector<TH1*> (n_cent_cuts)));


	double array_dr_bins[N_dR+1];
	for (int i_dR = 0; i_dR <= N_dR; i_dR++) array_dr_bins[i_dR] = dR_binning->GetBinLowEdge(i_dR+1);

	for (int i_cent = 0; i_cent < n_cent_cuts; i_cent++)
	{
		for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
		{
			for (int i_trk = 0; i_trk < N_trkpt; i_trk++)
			{
				name = Form("h_TM_MB_r_1D_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_TM_MB_r_1D.at(i_jet).at(i_trk).at(i_cent) = new TH1D(name.c_str(), name.c_str(), N_dR, array_dr_bins);

				name = Form("h_FS_MB_r_1D_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_FS_MB_r_1D.at(i_jet).at(i_trk).at(i_cent) = new TH1D(name.c_str(), name.c_str(), N_dR, array_dr_bins);

				name = Form("h_FNS_MB_r_1D_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				h_FNS_MB_r_1D.at(i_jet).at(i_trk).at(i_cent) = new TH1D(name.c_str(), name.c_str(), N_dR, array_dr_bins);
			}
		}
	}

	for (int i_cent = 0; i_cent < 6; i_cent++)
	{
		for (int i_dR = 0; i_dR < N_dR; i_dR++)
		{
			name = Form("h_reco_jet_spectrum_unW_y4_cent%i", i_cent);
			TH1* h_jet_spectra_unW = (TH1*)((TH1*)input_file->Get(name.c_str()))->Clone(Form("reco_jet_unW_y4_c%i", i_cent));
			h_jet_spectra_unW->Sumw2();

			name = Form("h_reco_jet_spectrum_y4_cent%i", i_cent);
			TH1* h_jet_spectra = (TH1*)((TH1*)input_file->Get(name.c_str()))->Clone(Form("reco_jet_y4_c%i", i_cent));
			h_jet_spectra->Sumw2();

			name = Form("h_true_jet_spectrum_y4_cent%i", i_cent);
			TH1* h_true_jet_spectra = (TH1*)((TH1*)input_file->Get(name.c_str()))->Clone(Form("true_jet_y4_c%i", i_cent));
			h_true_jet_spectra->Sumw2();


			name = Form("ChPS_raw_1_dR%i_cent%i", i_dR, i_cent);
			TH2* h_MB_method = (TH2*)input_file->Get(name.c_str());

			name = Form("ff_UE_pT_dR%i_cent%i", i_dR, i_cent);
			TH2* h_TM_method = (TH2*)input_file->Get(name.c_str());

			name = Form("ChPS_raw_2_dR%i_cent%i",i_dR,i_cent);
			TH2* h_FS = (TH2*)input_file->Get(name.c_str());

			name = Form("ChPS_raw_2_noSec_dR%i_cent%i",i_dR,i_cent);
			TH2* h_FNS = (TH2*)input_file->Get(name.c_str());


			for (int i_jet_bin = 1; i_jet_bin < N_jetpt; i_jet_bin++)
			{

				double n_jets_unw = h_jet_spectra_unW->GetBinContent(i_jet_bin);
				double n_jets = h_jet_spectra->GetBinContent(i_jet_bin);
				double n_jets_true = h_true_jet_spectra->GetBinContent(i_jet_bin);

				if (n_jets == 0) continue;

				for (int i_trk_bin = 1; i_trk_bin <= N_trkpt; i_trk_bin++)
				{
					double updated_UE_MB = h_MB_method->GetBinContent(i_trk_bin, i_jet_bin) / n_jets_unw;
					double updated_UE_TM = h_TM_method->GetBinContent(i_trk_bin, i_jet_bin) / n_jets;
					double updated_UE_FS = h_FS->GetBinContent(i_trk_bin, i_jet_bin) / n_jets_true;
					double updated_UE_FNS = h_FNS->GetBinContent(i_trk_bin, i_jet_bin) / n_jets_true;

					double updated_UE_MB_err = h_MB_method->GetBinError(i_trk_bin, i_jet_bin) / n_jets_unw;
					double updated_UE_TM_err = h_TM_method->GetBinError(i_trk_bin, i_jet_bin) / n_jets;
					double updated_UE_FS_err = h_FS->GetBinError(i_trk_bin, i_jet_bin) / n_jets_true;
					double updated_UE_FNS_err = h_FNS->GetBinError(i_trk_bin, i_jet_bin) / n_jets_true;

					h_MB_method->SetBinContent(i_trk_bin, i_jet_bin, updated_UE_MB);
					h_TM_method->SetBinContent(i_trk_bin, i_jet_bin, updated_UE_TM);
					h_FS->SetBinContent(i_trk_bin, i_jet_bin, updated_UE_FS);
					h_FNS->SetBinContent(i_trk_bin, i_jet_bin, updated_UE_FNS);

					h_MB_method->SetBinError(i_trk_bin, i_jet_bin, updated_UE_MB_err);
					h_TM_method->SetBinError(i_trk_bin, i_jet_bin, updated_UE_TM_err);
					h_FS->SetBinError(i_trk_bin, i_jet_bin, updated_UE_FS_err);
					h_FNS->SetBinError(i_trk_bin, i_jet_bin, updated_UE_FNS_err);
				}
			}


			name = Form("ratio_TM_MB_dR%i_cent%i", i_dR, i_cent);
			h_TM_MB_2D[i_dR][i_cent] = (TH2*)h_TM_method->Clone(name.c_str());
			h_TM_MB_2D[i_dR][i_cent]->Divide(h_MB_method);

			name = Form("ratio_FS_MB_dR%i_cent%i", i_dR, i_cent);
			h_FS_MB_2D[i_dR][i_cent] = (TH2*)h_FS->Clone(name.c_str());
			h_FS_MB_2D[i_dR][i_cent]->Divide(h_MB_method);

			name = Form("ratio_FNS_MB_dR%i_cent%i", i_dR, i_cent);
			h_FNS_MB_2D[i_dR][i_cent] = (TH2*)h_FNS->Clone(name.c_str());
			h_FNS_MB_2D[i_dR][i_cent]->Divide(h_MB_method);

			for (int i_jet_bin = jet_pt_start; i_jet_bin < jet_pt_end; i_jet_bin++)
			{
				h_TM_MB_1D[i_jet_bin][i_dR][i_cent] = (TH1*)h_TM_MB_2D[i_dR][i_cent]->ProjectionX(Form("TM_MB_%i_%i_%i", i_jet_bin, i_dR, i_cent), i_jet_bin+1, i_jet_bin+1);
				h_FS_MB_1D[i_jet_bin][i_dR][i_cent] = (TH1*)h_FS_MB_2D[i_dR][i_cent]->ProjectionX(Form("TFS_MB_%i_%i_%i", i_jet_bin, i_dR, i_cent), i_jet_bin+1, i_jet_bin+1);
				h_FNS_MB_1D[i_jet_bin][i_dR][i_cent] = (TH1*)h_FNS_MB_2D[i_dR][i_cent]->ProjectionX(Form("FNS_MB_%i_%i_%i", i_jet_bin, i_dR, i_cent), i_jet_bin+1, i_jet_bin+1);
			}

		}

//		for (int i_jet_bin = jet_pt_start; i_jet_bin < jet_pt_end; i_jet_bin++)
//		{
//			for (int i_trk_bin = 0; i_trk_bin < N_trkpt; i_trk_bin++)
//			{
//				for (int i_dR = 0; i_dR < 13; i_dR++)
//				{
//					h_TM_MB_r_1D[i_jet_bin][i_trk_bin][i_cent]->SetBinContent(i_dR+1, h_TM_MB_1D[i_jet_bin][i_dR][i_cent]->GetBinContent(i_trk_bin+1));
//					h_FS_MB_r_1D[i_jet_bin][i_trk_bin][i_cent]->SetBinContent(i_dR+1, h_FS_MB_1D[i_jet_bin][i_dR][i_cent]->GetBinContent(i_trk_bin+1));
//					h_FNS_MB_r_1D[i_jet_bin][i_trk_bin][i_cent]->SetBinContent(i_dR+1, h_FNS_MB_1D[i_jet_bin][i_dR][i_cent]->GetBinContent(i_trk_bin+1));
//				}
//			}
//		}
	}





	for (int i_dR = 0; i_dR < N_dR; i_dR++)
	{
		string dr_label = Form("%1.2f < r < %1.2f", dR_binning->GetBinLowEdge(i_dR+1), dR_binning->GetBinUpEdge(i_dR+1));

		c_TM_MB_2D->Clear();
		c_TM_MB_2D->Divide(3,2);

		c_FS_MB_2D->Clear();
		c_FS_MB_2D->Divide(3,2);

		c_FNS_MB_2D->Clear();
		c_FNS_MB_2D->Divide(3,2);

		for (int i_cent = 0; i_cent < 6; i_cent++)
		{

			h_TM_MB_2D[i_dR][i_cent]->GetYaxis()->SetTitle("p_{T}^{Jet}");
			h_FS_MB_2D[i_dR][i_cent]->GetYaxis()->SetTitle("p_{T}^{Jet}");
			h_FNS_MB_2D[i_dR][i_cent]->GetYaxis()->SetTitle("p_{T}^{Jet}");

			h_TM_MB_2D[i_dR][i_cent]->GetXaxis()->SetTitle("p_{T}^{Trk}");
			h_FS_MB_2D[i_dR][i_cent]->GetXaxis()->SetTitle("p_{T}^{Trk}");
			h_FNS_MB_2D[i_dR][i_cent]->GetXaxis()->SetTitle("p_{T}^{Trk}");

			h_TM_MB_2D[i_dR][i_cent]->GetXaxis()->SetRangeUser(1,10);
			h_FS_MB_2D[i_dR][i_cent]->GetXaxis()->SetRangeUser(1,10);
			h_FNS_MB_2D[i_dR][i_cent]->GetXaxis()->SetRangeUser(1,10);

			double low_range = 90;
			double hi_range = 400;

			h_TM_MB_2D[i_dR][i_cent]->GetYaxis()->SetRangeUser(low_range, hi_range);
			h_FS_MB_2D[i_dR][i_cent]->GetYaxis()->SetRangeUser(low_range, hi_range);
			h_FNS_MB_2D[i_dR][i_cent]->GetYaxis()->SetRangeUser(low_range, hi_range);

			h_TM_MB_2D[i_dR][i_cent]->GetZaxis()->SetRangeUser(0,2);
			h_FS_MB_2D[i_dR][i_cent]->GetZaxis()->SetRangeUser(0,2);
			h_FNS_MB_2D[i_dR][i_cent]->GetZaxis()->SetRangeUser(0,2);


			c_TM_MB_2D->cd(i_cent+1);
			h_TM_MB_2D[i_dR][i_cent]->Draw("colz text");

			ltx->SetTextAlign(12);
			ltx->DrawLatexNDC(0.19,0.98,dr_label.c_str());
			ltx->SetTextAlign(32);
			ltx->DrawLatexNDC(0.90,0.98,num_to_cent(31,i_cent).c_str());
			gPad->SetLogx();
			gPad->SetLogy();
//			gPad->SetLogz();

			c_FS_MB_2D->cd(i_cent+1);
			h_FS_MB_2D[i_dR][i_cent]->Draw("colz text");
			ltx->SetTextAlign(12);
			ltx->DrawLatexNDC(0.19,0.98,dr_label.c_str());
			ltx->SetTextAlign(32);
			ltx->DrawLatexNDC(0.90,0.98,num_to_cent(31,i_cent).c_str());
			gPad->SetLogx();
			gPad->SetLogy();
//			gPad->SetLogz();

			c_FNS_MB_2D->cd(i_cent+1);
			h_FNS_MB_2D[i_dR][i_cent]->Draw("colz text");
			ltx->SetTextAlign(12);
			ltx->DrawLatexNDC(0.19,0.98,dr_label.c_str());
			ltx->SetTextAlign(32);
			ltx->DrawLatexNDC(0.90,0.98,num_to_cent(31,i_cent).c_str());
			gPad->SetLogx();
			gPad->SetLogy();
//			gPad->SetLogz();

		}

		if (i_dR == 0) name = "(";
		else if (i_dR == 12) name = ")";
		else name = "";
		c_TM_MB_2D->Print(Form("UE_TM_MB_2D_up.pdf%s", name.c_str()), Form("Title: dR%i - %s", i_dR, dr_label.c_str()));
		c_FS_MB_2D->Print(Form("UE_FS_MB_2D_up.pdf%s", name.c_str()), Form("Title: dR%i - %s", i_dR, dr_label.c_str()));
		c_FNS_MB_2D->Print(Form("UE_FNS_MB_2D_up.pdf%s", name.c_str()), Form("Title: dR%i - %s", i_dR, dr_label.c_str()));
	}


	for (int i_dR = 0; i_dR < N_dR; i_dR++)
	{
		string dr_label = Form("%1.2f < r < %1.2f", dR_binning->GetBinLowEdge(i_dR+1), dR_binning->GetBinUpEdge(i_dR+1));

		c_TM_MB_1D->Clear();
		c_TM_MB_1D->Divide(3,2);
		
		c_FS_MB_1D->Clear();
		c_FS_MB_1D->Divide(3,2);
		
		c_FNS_MB_1D->Clear();
		c_FNS_MB_1D->Divide(3,2);

		for (int i_cent = 0; i_cent < 6; i_cent++)
		{

			int jet_itr = 0;
			for (int i_jet_bin = jet_pt_start; i_jet_bin < jet_pt_end; i_jet_bin++)
			{
				string jet_label = Form("%1.0f < p_{T}^{Jet} < %1.0f", jetpT_binning->GetBinLowEdge(i_jet_bin+1), jetpT_binning->GetBinUpEdge(i_jet_bin+1));


				SetHStyle_smallify(h_TM_MB_1D[i_jet_bin][i_dR][i_cent], jet_itr, 1);
				SetHStyle_smallify(h_FS_MB_1D[i_jet_bin][i_dR][i_cent], jet_itr, 1);
				SetHStyle_smallify(h_FNS_MB_1D[i_jet_bin][i_dR][i_cent], jet_itr, 1);


				if (i_cent == 0 && i_dR == 0) legend_TM_MB->AddEntry(h_TM_MB_1D[i_jet_bin][i_dR][i_cent],jet_label.c_str(),"lp");
				if (i_cent == 0 && i_dR == 0) legend_FS_MB->AddEntry(h_FS_MB_1D[i_jet_bin][i_dR][i_cent],jet_label.c_str(),"lp");
				if (i_cent == 0 && i_dR == 0) legend_FNS_MB->AddEntry(h_FNS_MB_1D[i_jet_bin][i_dR][i_cent],jet_label.c_str(),"lp");


				h_TM_MB_1D[i_jet_bin][i_dR][i_cent]->GetXaxis()->SetRangeUser(1,10);
				h_FS_MB_1D[i_jet_bin][i_dR][i_cent]->GetXaxis()->SetRangeUser(1,10);
				h_FNS_MB_1D[i_jet_bin][i_dR][i_cent]->GetXaxis()->SetRangeUser(1,10);

				double low_range = 0.5;
				double hi_range = 1.5;
				h_TM_MB_1D[i_jet_bin][i_dR][i_cent]->GetYaxis()->SetRangeUser(low_range, hi_range);
				h_FS_MB_1D[i_jet_bin][i_dR][i_cent]->GetYaxis()->SetRangeUser(low_range, hi_range);
				h_FNS_MB_1D[i_jet_bin][i_dR][i_cent]->GetYaxis()->SetRangeUser(low_range, hi_range);


				c_TM_MB_1D->cd(i_cent+1);
				if (i_jet_bin == 7) h_TM_MB_1D[i_jet_bin][i_dR][i_cent]->Draw("");
				else h_TM_MB_1D[i_jet_bin][i_dR][i_cent]->Draw("same");
				gPad->SetLogx();
//				gPad->SetLogy();

				c_FS_MB_1D->cd(i_cent+1);
				if (i_jet_bin == 7) h_FS_MB_1D[i_jet_bin][i_dR][i_cent]->Draw("");
				else h_FS_MB_1D[i_jet_bin][i_dR][i_cent]->Draw("same");
				gPad->SetLogx();
//				gPad->SetLogy();

				c_FNS_MB_1D->cd(i_cent+1);
				if (i_jet_bin == 7) h_FNS_MB_1D[i_jet_bin][i_dR][i_cent]->Draw("");
				else h_FNS_MB_1D[i_jet_bin][i_dR][i_cent]->Draw("same");
				gPad->SetLogx();
//				gPad->SetLogy();

				jet_itr++;
			}

			c_TM_MB_1D->cd(i_cent+1);
			ltx->SetTextAlign(12);
			ltx->DrawLatexNDC(0.19,0.98,dr_label.c_str());
			ltx->SetTextAlign(32);
			ltx->DrawLatexNDC(0.90,0.98,num_to_cent(31,i_cent).c_str());
			legend_TM_MB->Draw();

			c_FS_MB_1D->cd(i_cent+1);
			ltx->SetTextAlign(12);
			ltx->DrawLatexNDC(0.19,0.98,dr_label.c_str());
			ltx->SetTextAlign(32);
			ltx->DrawLatexNDC(0.90,0.98,num_to_cent(31,i_cent).c_str());
			legend_FS_MB->Draw();

			c_FNS_MB_1D->cd(i_cent+1);
			ltx->SetTextAlign(12);
			ltx->DrawLatexNDC(0.19,0.98,dr_label.c_str());
			ltx->SetTextAlign(32);
			ltx->DrawLatexNDC(0.90,0.98,num_to_cent(31,i_cent).c_str());
			legend_FNS_MB->Draw();
		}
		if (i_dR == 0) name = "(";
		else if (i_dR == 12) name = ")";
		else name = "";
		c_TM_MB_1D->Print(Form("c_TM_MB_1D.pdf%s", name.c_str()), Form("Title: dR%i - %s", i_dR, dr_label.c_str()));
		c_FS_MB_1D->Print(Form("c_FS_MB_1D.pdf%s", name.c_str()), Form("Title: dR%i - %s", i_dR, dr_label.c_str()));
		c_FNS_MB_1D->Print(Form("c_FNS_MB_1D.pdf%s", name.c_str()), Form("Title: dR%i - %s", i_dR, dr_label.c_str()));
	}

}
//

//			c_TM_MB_1D->Clear();
//			c_TM_MB_1D->Divide(3,2);
//
//			c_FS_MB_1D->Clear();
//			c_FS_MB_1D->Divide(3,2);
//
//			c_FNS_MB_1D->Clear();
//			c_FNS_MB_1D->Divide(3,2);

//	int jet_itr = 0;
//	for (int i_jet_bin = jet_pt_start; i_jet_bin < jet_pt_end; i_jet_bin++)
//	{
//
//		string jet_label = Form("%1.0f < p_{T}^{Jet} < %1.0f", jetpT_binning->GetBinLowEdge(i_jet_bin+1), jetpT_binning->GetBinUpEdge(i_jet_bin+1));
//
//
//		SetHStyle_smallify(h_TM_MB_1D[i_jet_bin][i_dR][i_cent], jet_itr, 1);
//		SetHStyle_smallify(h_FS_MB_1D[i_jet_bin][i_dR][i_cent], jet_itr, 1);
//		SetHStyle_smallify(h_FNS_MB_1D[i_jet_bin][i_dR][i_cent], jet_itr, 1);
//
//
//		if (i_cent == 0 && i_dR == 0) legend_TM_MB->AddEntry(h_TM_MB_1D[i_jet_bin][i_dR][i_cent],jet_label.c_str(),"lp");
//		if (i_cent == 0 && i_dR == 0) legend_FS_MB->AddEntry(h_FS_MB_1D[i_jet_bin][i_dR][i_cent],jet_label.c_str(),"lp");
//		if (i_cent == 0 && i_dR == 0) legend_FNS_MB->AddEntry(h_FNS_MB_1D[i_jet_bin][i_dR][i_cent],jet_label.c_str(),"lp");
//
//
//		h_TM_MB_1D[i_jet_bin][i_dR][i_cent]->GetXaxis()->SetRangeUser(1,10);
//		h_FS_MB_1D[i_jet_bin][i_dR][i_cent]->GetXaxis()->SetRangeUser(1,10);
//		h_FNS_MB_1D[i_jet_bin][i_dR][i_cent]->GetXaxis()->SetRangeUser(1,10);
//
//
//		h_TM_MB_1D[i_jet_bin][i_dR][i_cent]->GetYaxis()->SetRangeUser(0.95,1.05);
//		h_FS_MB_1D[i_jet_bin][i_dR][i_cent]->GetYaxis()->SetRangeUser(0.98,1.02);
//		h_FNS_MB_1D[i_jet_bin][i_dR][i_cent]->GetYaxis()->SetRangeUser(0.98,1.02);
//
//		c_TM_MB_1D->cd(i_cent+1);
//		if (i_jet_bin == 7) h_TM_MB_1D[i_jet_bin][i_dR][i_cent]->Draw("");
//		else h_TM_MB_1D[i_jet_bin][i_dR][i_cent]->Draw("same");
//		gPad->SetLogx();
//
//		c_FS_MB_1D->cd(i_cent+1);
//		if (i_jet_bin == 7) h_FS_MB_1D[i_jet_bin][i_dR][i_cent]->Draw("");
//		else h_FS_MB_1D[i_jet_bin][i_dR][i_cent]->Draw("same");
//		gPad->SetLogx();
//
//		c_FNS_MB_1D->cd(i_cent+1);
//		if (i_jet_bin == 7) h_FNS_MB_1D[i_jet_bin][i_dR][i_cent]->Draw("");
//		else h_FNS_MB_1D[i_jet_bin][i_dR][i_cent]->Draw("same");
//		gPad->SetLogx();
//
//		jet_itr++;
//	}
//
//	c_TM_MB_1D->cd(i_cent+1);
//	ltx->SetTextAlign(12);
//	ltx->DrawLatexNDC(0.19,0.98,dr_label.c_str());
//	ltx->SetTextAlign(32);
//	ltx->DrawLatexNDC(0.90,0.98,num_to_cent(31,i_cent).c_str());
//	legend_TM_MB->Draw();
//
//	c_FS_MB_1D->cd(i_cent+1);
//	ltx->SetTextAlign(12);
//	ltx->DrawLatexNDC(0.19,0.98,dr_label.c_str());
//	ltx->SetTextAlign(32);
//	ltx->DrawLatexNDC(0.90,0.98,num_to_cent(31,i_cent).c_str());
//	legend_FS_MB->Draw();
//
//	c_FNS_MB_1D->cd(i_cent+1);
//	ltx->SetTextAlign(12);
//	ltx->DrawLatexNDC(0.19,0.98,dr_label.c_str());
//	ltx->SetTextAlign(32);
//	ltx->DrawLatexNDC(0.90,0.98,num_to_cent(31,i_cent).c_str());
//	legend_FNS_MB->Draw();
//}
//	if (i_dR == 0) name = "(";
//	else if (i_dR == 12) name = ")";
//	else name = "";
//	c_TM_MB_2D->Print(Form("UE_TM_MB_2D.pdf%s", name.c_str()), Form("Title: dR%i - %s", i_dR, dr_label.c_str()));
//	c_FS_MB_2D->Print(Form("UE_FS_MB_2D.pdf%s", name.c_str()), Form("Title: dR%i - %s", i_dR, dr_label.c_str()));
//	c_FNS_MB_2D->Print(Form("UE_FNS_MB_2D.pdf%s", name.c_str()), Form("Title: dR%i - %s", i_dR, dr_label.c_str()));
//
//c_TM_MB_1D->Print(Form("c_TM_MB_1D.pdf%s", name.c_str()), Form("Title: dR%i - %s", i_dR, dr_label.c_str()));
//c_FS_MB_1D->Print(Form("c_FS_MB_1D.pdf%s", name.c_str()), Form("Title: dR%i - %s", i_dR, dr_label.c_str()));
//c_FNS_MB_1D->Print(Form("c_FNS_MB_1D.pdf%s", name.c_str()), Form("Title: dR%i - %s", i_dR, dr_label.c_str()));
//}
//}
//}

