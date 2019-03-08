#include "output_dev/functions/global_variables.h"

void draw_run_dep()
{
	SetAtlasStyle();
	gErrorIgnoreLevel = 3001;

	TFile *f_nominal_full = new TFile("UE_c47.root");
	TFile *f_nominal = new TFile("UE_c67.root");
	TFile *f_UE_sys = new TFile("UE_RunDependentSys.root","recreate");

	vector<TFile*> f_run_UE;

	f_run_UE.push_back(new TFile("UE_c68.root"));
	f_run_UE.push_back(new TFile("UE_c69.root"));
	f_run_UE.push_back(new TFile("UE_c70.root"));

	double r_max_range = 0.8;

	vector<TH1*> h_run_UE = vector<TH1*> (f_run_UE.size());

	string name;
	TCanvas *c = new TCanvas("c","c",900,600);
	TCanvas *c1 = new TCanvas("c1","c1",900,600);
	TCanvas *c2 = new TCanvas("c2","c2",900,600);

	TLegend *legend_MB = new TLegend(0.19, 0.70, 0.40, 0.82, "","brNDC");
	legend_MB->SetTextFont(43);
	legend_MB->SetBorderSize(0);
	legend_MB->SetTextSize(10);
	TLine *line = new TLine();

	std::map<std::string, std::string> m;
	m["UE_c68.root"] = "Period I";
	m["UE_c69.root"] = "Period II";
	m["UE_c70.root"] = "Period III";

	TAxis* dR_binning = (TAxis*)f_nominal->Get("dR_binning");
	TAxis* jetpT_binning = (TAxis*)f_nominal->Get("jetpT_binning");
	TAxis* trkpT_binning = (TAxis*)f_nominal->Get("trkpT_binning");

	int N_dR = dR_binning->GetNbins();
	int N_jetpt = jetpT_binning->GetNbins();
	int N_trkpt = trkpT_binning->GetNbins();

	TLatex *ltx = new TLatex();
	ltx->SetTextFont(43);
	ltx->SetTextSize(12);
	ltx->SetTextAlign(12);

	//i_jet is the vector index, corresponds to bin = i_jet+1, i_jet = 6 -> Bin 7 for 100-126, i_jet = 11 -> bin 12 for 316 - 398
	for (int i_jet = 6; i_jet < 12; i_jet++)
	{

		//i_trk is the vector index, corresponds to bin = i_trk+1, i_trk = 1 -> Bin 2 for 0.9-1, i_trk = 6 -> bin 7 for 6.3 - 10
		for (int i_trk = 1; i_trk < 7; i_trk++)
		{

			string trk_label = Form("%1.1f < p_{T}^{Trk} < %1.1f GeV", trkpT_binning->GetBinLowEdge(i_trk+1), trkpT_binning->GetBinUpEdge(i_trk+1));

			string jet_label = Form("%1.0f < p_{T}^{Jet} < %1.0f GeV", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));


			c->cd();
			c->Clear();
			c->Divide(3,2);

			c1->cd();
			c1->Clear();
			c1->Divide(3,2);

			c2->cd();
			c2->Clear();
			c2->Divide(3,2);

			for (int i_cent = 0; i_cent < 6; i_cent++)
			{
				string which_UE = Form("UE_MB_data_indR_jet%i_trk%i_cent%i", i_jet, i_trk, i_cent);

				TH1* h_nominal = (TH1*)f_nominal->Get(which_UE.c_str());

				TH1* h_combined = (TH1*)h_nominal->Clone("h_combined");
				h_combined->Reset();

				TH1* h_systematic = (TH1*)h_nominal->Clone("h_systematic");
				h_systematic->Reset();

				double total = 0, weight = 0;
				for (int i_run = 0; i_run < f_run_UE.size(); i_run++)
				{
					name = Form("h_reco_jet_spectrum_y4_cent%i", i_cent);
					weight = ((TH1*)f_run_UE[i_run]->Get(name.c_str()))->GetBinContent(i_jet+1);
//					weight = ((TH1*)f_run_UE[i_run]->Get("Centrality"))->GetBinContent(i_cent+1);

					h_run_UE[i_run] = (TH1*)f_run_UE[i_run]->Get(which_UE.c_str());
					total = total+weight;
					h_combined->Add(h_run_UE[i_run], weight);

					f_UE_sys->cd();
					SetHStyle(h_run_UE[i_run],i_run+1);
					h_run_UE[i_run]->Write(Form("%s_run%i",which_UE.c_str(),i_run));
				}
				h_combined->Scale(1./total);

				TH1* h_nominal_full = (TH1*)f_nominal_full->Get(which_UE.c_str());
				h_nominal_full->SetTitle("UE_{All Runs (Full map)}");
				h_nominal->SetTitle("UE_{All Runs (Reduced map)}");
				TH1* h_up = (TH1*)h_nominal_full->Clone("h_up");
				h_up->SetTitle("UE_{Shifted Up}");
				TH1* h_down = (TH1*)h_nominal_full->Clone("h_down");
				h_down->SetTitle("UE_{Shifted Down}");

				for (int i_dR = 1; i_dR <= h_systematic->GetXaxis()->GetNbins() ; i_dR++)
				{
					double max = 0, delta = 0;
					for (int i_run = 0; i_run < f_run_UE.size(); i_run++)
					{
						//max delta % = (max deviation from nominal/ nominal, do 1+maxDelta% or 1-maxDelta%
						delta = abs(h_run_UE[i_run]->GetBinContent(i_dR) - h_nominal->GetBinContent(i_dR)) / h_nominal->GetBinContent(i_dR);
						if (delta > max) max = delta;
					}
					h_up->SetBinContent(i_dR,h_nominal_full->GetBinContent(i_dR)*(1+max));
					h_down->SetBinContent(i_dR,h_nominal_full->GetBinContent(i_dR)*(1-max));
					h_systematic->SetBinContent(i_dR, 1+max);
				}


				c1->cd(i_cent+1);
				SetHStyle_smallify(h_nominal_full,0,1);
				SetHStyle_smallify(h_up,1,1);
				SetHStyle_smallify(h_down,2,1);
				h_nominal_full->GetYaxis()->SetRangeUser(h_nominal_full->GetMaximum()*0.85, h_nominal_full->GetMaximum()*1.08);
				h_nominal_full->GetXaxis()->SetRangeUser(0,0.8);
				h_nominal_full->DrawCopy();
				h_up->DrawCopy("same");
				h_down->DrawCopy("same");
				gPad->BuildLegend();
//				gPad->BuildLegend(0.5,0.67,0.88,0.88,"","NDC");


				c2->cd(i_cent+1);
				TH1* h_ratio_full_reduced = (TH1*)h_nominal->Clone("h_ratio_full_reduced");
				h_ratio_full_reduced->Divide(h_nominal_full);
				h_ratio_full_reduced->SetTitle("Reduced Maps / Full Maps");
				SetHStyle_smallify(h_ratio_full_reduced,1,1);
				h_ratio_full_reduced->GetYaxis()->SetTitle("Ratio");
				h_ratio_full_reduced->GetXaxis()->SetRangeUser(0,0.8);
				h_ratio_full_reduced->GetYaxis()->SetRangeUser(0.94,1.06);
				h_ratio_full_reduced->DrawCopy("");
				gPad->BuildLegend();
				line->DrawLine(0,1,0.8,1);



				f_UE_sys->cd();
				SetHStyle(h_nominal,0);
				name = Form("%s_nom",which_UE.c_str());
				h_nominal->SetName(name.c_str());
				h_nominal->Write(name.c_str());

				SetHStyle(h_systematic,f_run_UE.size());
				name = Form("%s_SYS",which_UE.c_str());
				h_systematic->SetName(name.c_str());
				h_systematic->Write(name.c_str());



				//Drawing
				double min = 0.95;//h_nominal->GetBinContent(7) * 0.6;// 0.5;
				double max = 1.05;//h_nominal->GetBinContent(1) * 1.4;// 0.5;
//				double min = h_nominal->GetBinContent(7) * 0.6;// 0.5;
//				double max = h_nominal->GetBinContent(1) * 1.4;// 0.5;

				if (i_cent == 0 && i_trk == 2 && i_jet == 7)
				{
					for (int i_run = 0; i_run < f_run_UE.size(); i_run++) legend_MB->AddEntry(h_run_UE[i_run],m[f_run_UE[i_run]->GetName()].c_str(), "lp");
				}
				c->cd(i_cent+1);
				h_nominal->GetXaxis()->SetRangeUser(0, r_max_range);
				h_nominal->GetYaxis()->SetRangeUser(min, max);
				SetHStyle_smallify(h_nominal,1,1);
				SetHStyle_smallify(h_combined,0,1);

				h_nominal->GetYaxis()->SetTitle("Period_{i} / Nominal");
				h_nominal->DrawCopy("");
				for (int i_run = 0; i_run < f_run_UE.size(); i_run++)
				{
					SetHStyle_smallify(h_run_UE[i_run],i_run+2,1);
					h_run_UE[i_run]->Divide(h_nominal);
					h_run_UE[i_run]->DrawCopy("same");
				}
//				h_combined->Divide(h_nominal);
//				h_combined->Draw("same hist text");

				ltx->SetTextAlign(32);
				ltx->DrawLatexNDC(0.92,0.92,num_to_cent(31,i_cent).c_str());
				ltx->DrawLatexNDC(0.92,0.86,jet_label.c_str());
				ltx->DrawLatexNDC(0.92,0.80,trk_label.c_str());
				line->DrawLine(0,1,0.8,1);

			}

			c->cd(1);
			legend_MB->Draw();
			c->cd();
			string pdf_label = "";
			if (i_jet == 7 && i_trk == 2) pdf_label = "(";
			if (i_jet == 11 && i_trk == 6) pdf_label = ")";
			c->Print(Form("weightedRuns.pdf%s",pdf_label.c_str()),Form("Title: jet%i_trk%i", i_jet, i_trk));
			c1->Print(Form("weightedRuns_v_Nominal.pdf%s",pdf_label.c_str()),Form("Title: jet%i_trk%i", i_jet, i_trk));
			c2->Print(Form("Nominal_full_v_reducedMaps.pdf%s",pdf_label.c_str()),Form("Title: jet%i_trk%i", i_jet, i_trk));

		}
	}


}
