#include "../functions/global_variables.h"

void UE_factors()
{
	SetAtlasStyle();
	gErrorIgnoreLevel = 3001;

	TFile *input_file = new TFile("../raw_results/FF_MC_out_histo_PbPb_5p02_r001.root");
	TFile *UE_factors = new TFile(Form("UE_factors.root"), "recreate");
	TAxis* dR_binning = (TAxis*)((TH3*)input_file->Get("h_dR_change_jetpt0_cent0"))->GetXaxis();
	TAxis* jetpT_binning = (TAxis*)((TH3*)input_file->Get("ChPS_raw_0_dR0_cent0"))->GetYaxis();
	TAxis* trkpT_binning = (TAxis*)((TH3*)input_file->Get("ChPS_raw_0_dR0_cent0"))->GetXaxis();

	vector<vector<TH2*>> h_ratio = vector<vector<TH2*>> (13, vector<TH2*> (6));

	string name;
	TCanvas *c1 = new TCanvas("c1","c1",800,400);
	TLatex *ltx = new TLatex();
	ltx->SetTextFont(43);
	ltx->SetTextSize(11);
	ltx->SetTextAlign(12);

	int N_dR = dR_binning->GetNbins();
	int N_jetpt = jetpT_binning->GetNbins();
	int N_trkpt = trkpT_binning->GetNbins();


	for (int i_dR = 0; i_dR < 13; i_dR++)
	{
		string dr_label = Form("%1.2f < r < %1.2f", dR_binning->GetBinLowEdge(i_dR+1), dR_binning->GetBinUpEdge(i_dR+1));

		c1->Clear();
		c1->Divide(3,2);
		for (int i_cent = 0; i_cent < 6; i_cent++)
		{
			name = Form("ChPS_raw_1_dR%i_cent%i", i_dR, i_cent);
			TH2* h_MB_method = (TH2*)input_file->Get(name.c_str());
			name = Form("h_reco_jet_spectrum_unW_y4_cent%i", i_cent);
			TH1* h_jet_spectra_unW = (TH1*)((TH1*)input_file->Get(name.c_str()))->Clone(Form("reco_uw_jet_y4_c%i", i_cent));
			h_jet_spectra_unW->Sumw2();


			name = Form("ff_UE_pT_dR%i_cent%i", i_dR, i_cent);
			TH2* h_TM_method = (TH2*)input_file->Get(name.c_str());
			name = Form("h_reco_jet_spectrum_y4_cent%i", i_cent);
			TH1* h_jet_spectra = (TH1*)((TH1*)input_file->Get(name.c_str()))->Clone(Form("reco_w_jet_y4_c%i", i_cent));
			h_jet_spectra->Sumw2();


			for (int i_jet_bin = 1; i_jet_bin <= N_jetpt; i_jet_bin++)
			{
				double n_jets_w = h_jet_spectra->GetBinContent(i_jet_bin);
				double n_jets_unw = h_jet_spectra_unW->GetBinContent(i_jet_bin);

				if (n_jets_w == 0) continue;

				for (int i_trk_bin = 1; i_trk_bin <= N_jetpt; i_trk_bin++)
				{
					double updated_UE_MB = h_MB_method->GetBinContent(i_trk_bin, i_jet_bin) / n_jets_unw;
					double updated_UE_TM = h_TM_method->GetBinContent(i_trk_bin, i_jet_bin) / n_jets_w;

					double updated_UE_MB_err = h_MB_method->GetBinError(i_trk_bin, i_jet_bin) / n_jets_unw;
					double updated_UE_TM_err = h_TM_method->GetBinError(i_trk_bin, i_jet_bin) / n_jets_w;

					h_MB_method->SetBinContent(i_trk_bin, i_jet_bin, updated_UE_MB);
					h_TM_method->SetBinContent(i_trk_bin, i_jet_bin, updated_UE_TM);

					h_MB_method->SetBinError(i_trk_bin, i_jet_bin, updated_UE_MB_err);
					h_TM_method->SetBinError(i_trk_bin, i_jet_bin, updated_UE_TM_err);
				}
			}


			name = Form("ratio_dR%i_cent%i", i_dR, i_cent);
			h_ratio[i_dR][i_cent] = (TH2*)h_TM_method->Clone(name.c_str());
			h_ratio[i_dR][i_cent]->Divide(h_MB_method);

			UE_factors->cd();
			name = Form("UE_ratio_dR%i_cent%i", i_dR, i_cent);
			h_ratio[i_dR][i_cent]->Write(name.c_str());

			c1->cd(i_cent+1);
			h_ratio[i_dR][i_cent]->Draw("colz text");
			h_ratio[i_dR][i_cent]->GetYaxis()->SetTitle("p_{T}^{Jet}");
			h_ratio[i_dR][i_cent]->GetXaxis()->SetTitle("p_{T}^{Trk}");

			h_ratio[i_dR][i_cent]->GetXaxis()->SetRangeUser(1,10);
			h_ratio[i_dR][i_cent]->GetYaxis()->SetRangeUser(90,500);
			h_ratio[i_dR][i_cent]->GetZaxis()->SetRangeUser(0,2);
			ltx->SetTextAlign(12);
			ltx->DrawLatexNDC(0.19,0.98,dr_label.c_str());
			ltx->SetTextAlign(32);
			ltx->DrawLatexNDC(0.90,0.98,num_to_cent(31,i_cent).c_str());
			gPad->SetLogx();
			gPad->SetLogy();
		}
		if (i_dR == 0) name = "(";
		else if (i_dR == 12) name = ")";
		else name = "";
		c1->Print(Form("UE_factors.pdf%s", name.c_str()), Form("Title: dR%i - %s", i_dR, dr_label.c_str()));
	}

}
