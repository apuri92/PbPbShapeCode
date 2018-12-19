#include "output_dev/functions/global_variables.h"

void grid_phi_dep()
{
	SetAtlasStyle();
	gErrorIgnoreLevel = 3001;

	TFile *f_fixed00_corr = new TFile("UE_c15_Corr.root");
	TFile *f_fixed01_corr = new TFile("UE_c14_Corr.root");
	TFile *f_float_corr = new TFile("UE_c16_Corr.root");
	TFile *f_fixed06_corr = new TFile("UE_c17_Corr.root");

	TFile *f_fixed00_nocorr = new TFile("UE_c15_noCorr.root");
	TFile *f_fixed01_nocorr = new TFile("UE_c14_noCorr.root");
	TFile *f_float_nocorr = new TFile("UE_c16_noCorr.root");
	TFile *f_fixed06_nocorr = new TFile("UE_c17_noCorr.root");


	TCanvas *c = new TCanvas("c","c",900,600);

	TAxis* dR_binning = (TAxis*)f_fixed00_corr->Get("dR_binning");
	TAxis* jetpT_binning = (TAxis*)f_fixed00_corr->Get("jetpT_binning");
	TAxis* trkpT_binning = (TAxis*)f_fixed00_corr->Get("trkpT_binning");


	string name;
	TLegend *legend = new TLegend(0.4, 0.20, 0.70, 0.30, "","brNDC");
	legend->SetTextFont(43);
	legend->SetBorderSize(0);
	legend->SetTextSize(10);

	TLine *line = new TLine();
	TLatex *ltx = new TLatex();

	bool first_pass = true;
	for (int i_trk = 2; i_trk < 7; i_trk++)
	{
		string trk_label = Form("%1.1f < p_{T}^{Trk} < %1.1f GeV", trkpT_binning->GetBinLowEdge(i_trk+1), trkpT_binning->GetBinUpEdge(i_trk+1));

		for (int i_jet = 7; i_jet < 11; i_jet++)
		{
			c->Clear();
			c->Divide(3,2);
			string jet_label = Form("%1.0f < p_{T}^{Jet} < %1.0f GeV", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));

			for (int i_cent = 0; i_cent < 6; i_cent++)
			{

				string name = Form("UE_cone_data_indR_jet%i_trk%i_cent%i", i_jet, i_trk, i_cent);


				TH1* h_fixed00_corr = (TH1*)f_fixed00_corr->Get(name.c_str());
				h_fixed00_corr->SetName("h_fixed00_corr");
				TH1* h_fixed00_nocorr = (TH1*)f_fixed00_nocorr->Get(name.c_str());
				h_fixed00_nocorr->SetName("h_fixed00_nocorr");

				TH1* h_fixed01_corr = (TH1*)f_fixed01_corr->Get(name.c_str());
				h_fixed01_corr->SetName("h_fixed01_corr");
				TH1* h_fixed01_nocorr = (TH1*)f_fixed01_nocorr->Get(name.c_str());
				h_fixed01_nocorr->SetName("h_fixed01_nocorr");

				TH1* h_fixed06_corr = (TH1*)f_fixed06_corr->Get(name.c_str());
				h_fixed06_corr->SetName("h_fixed06_corr");
				TH1* h_fixed06_nocorr = (TH1*)f_fixed06_nocorr->Get(name.c_str());
				h_fixed06_nocorr->SetName("h_fixed06_nocorr");

				TH1* h_float_corr = (TH1*)f_float_corr->Get(name.c_str());
				h_float_corr->SetName("h_float_corr");
				TH1* h_float_nocorr = (TH1*)f_float_nocorr->Get(name.c_str());
				h_float_nocorr->SetName("h_float_nocorr");

				TH1* h_ratio_00_corr = (TH1*)h_fixed00_corr->Clone(Form("h_00_corr_%s", name.c_str()));
				TH1* h_ratio_00_nocorr = (TH1*)h_fixed00_nocorr->Clone(Form("h_00_nocorr_%s", name.c_str()));

				TH1* h_ratio_01_corr = (TH1*)h_fixed01_corr->Clone(Form("h_01_corr_%s", name.c_str()));
				TH1* h_ratio_01_nocorr = (TH1*)h_fixed01_nocorr->Clone(Form("h_01_nocorr_%s", name.c_str()));

				TH1* h_ratio_06_corr = (TH1*)h_fixed06_corr->Clone(Form("h_06_corr_%s", name.c_str()));
				TH1* h_ratio_06_nocorr = (TH1*)h_fixed06_nocorr->Clone(Form("h_06_nocorr_%s", name.c_str()));

				TH1* h_ratio_00_01_corr = (TH1*)h_fixed01_corr->Clone(Form("h_00_01_corr_%s", name.c_str()));
				h_ratio_00_01_corr->Divide(h_fixed00_corr);

				TH1* h_ratio_00_01_nocorr = (TH1*)h_fixed01_nocorr->Clone(Form("h_00_01_nocorr_%s", name.c_str()));
				h_ratio_00_01_nocorr->Divide(h_fixed00_nocorr);

				h_ratio_00_corr->Divide(h_float_corr);
				h_ratio_00_nocorr->Divide(h_float_nocorr);
				h_ratio_01_corr->Divide(h_float_corr);
				h_ratio_01_nocorr->Divide(h_float_nocorr);
				h_ratio_06_corr->Divide(h_float_corr);
				h_ratio_06_nocorr->Divide(h_float_nocorr);

				SetHStyle_smallify(h_ratio_00_corr,0,1);
				SetHStyle_smallify(h_ratio_00_nocorr,1,1);

				SetHStyle_smallify(h_ratio_01_corr,2,1);
				SetHStyle_smallify(h_ratio_01_nocorr,3,1);

				SetHStyle_smallify(h_ratio_06_corr,4,1);
				SetHStyle_smallify(h_ratio_06_nocorr,5,1);

				SetHStyle_smallify(h_ratio_00_01_corr,6,1);
				SetHStyle_smallify(h_ratio_00_01_nocorr,7,1);

				if (first_pass)
				{
					legend->AddEntry(h_ratio_00_corr,"(#delta#phi = 0.0/#delta#phi = Rndm) (Corr.)");
					legend->AddEntry(h_ratio_01_corr,"(#delta#phi = 0.1/#delta#phi = Rndm) (Corr.)");
					legend->AddEntry(h_ratio_06_corr,"(#delta#phi = 0.6/#delta#phi = Rndm) (Corr.)");

//					legend->AddEntry(h_ratio_00_nocorr,"(#delta#phi = 0.0/#delta#phi = Rndm) (Not Corr.)");
//					legend->AddEntry(h_ratio_01_nocorr,"(#delta#phi = 0.1/#delta#phi = Rndm) (Not Corr.)");
//					legend->AddEntry(h_ratio_06_nocorr,"(#delta#phi = 0.6/#delta#phi = Rndm) (Not Corr.)");

//					legend->AddEntry(h_ratio_00_01_corr,"(#delta#phi = 0.1/#delta#phi = 0.0) (Corr.)");
//					legend->AddEntry(h_ratio_00_01_nocorr,"(#delta#phi = 0.1/#delta#phi = 0.0) (Not Corr.)");
					first_pass = false;
				}

				c->cd(i_cent+1);
				h_ratio_00_corr->GetYaxis()->SetRangeUser(0.98,1.02);
				h_ratio_00_corr->GetXaxis()->SetRangeUser(0.,0.6);

				h_ratio_00_corr->DrawCopy("");
//				h_ratio_00_nocorr->DrawCopy("same");
				h_ratio_01_corr->DrawCopy("same");
//				h_ratio_01_nocorr->DrawCopy("same");
				h_ratio_06_corr->DrawCopy("same");
//				h_ratio_06_nocorr->DrawCopy("same");
//				h_ratio_00_01_corr->DrawCopy("same");
//				h_ratio_00_01_nocorr->DrawCopy("same");


				ltx->SetTextAlign(32);
				ltx->DrawLatexNDC(0.94,0.92,num_to_cent(31,i_cent).c_str());
				ltx->DrawLatexNDC(0.94,0.86,jet_label.c_str());
				ltx->DrawLatexNDC(0.94,0.80,trk_label.c_str());

				line->DrawLine(0, 1, 0.6, 1);

			}
			c->cd(1);
			legend->Draw();
			if (i_trk == 2 && i_jet == 7) name = "(";
			else if (i_trk == 6 && i_jet == 10) name = ")";
			else name = "";
			c->Print(Form("phi_shift.pdf%s", name.c_str()), Form("Title: trk%i_jet%i", i_trk, i_jet));
		}
	}
}
