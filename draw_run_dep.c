#include "output_dev/functions/global_variables.h"

void draw_run_dep()
{
	gErrorIgnoreLevel = 3001;

	TFile *f_nominal = new TFile("UE_c67.root");
	TFile *f_UE_sys = new TFile("UE_sys.root","recreate");

	vector<TFile*> f_run_UE;

	f_run_UE.push_back(new TFile("UE_c68.root"));
	f_run_UE.push_back(new TFile("UE_c69.root"));
	f_run_UE.push_back(new TFile("UE_c70.root"));

	double r_max_range = 0.8;

	vector<TH1*> h_run_UE = vector<TH1*> (f_run_UE.size());

	string name;
	TCanvas *c = new TCanvas("c","c",800,600);

	TLegend *legend_MB = new TLegend(0.19, 0.70, 0.40, 0.82, "","brNDC");
	legend_MB->SetTextFont(43);
	legend_MB->SetBorderSize(0);
	legend_MB->SetTextSize(10);
	TLine *line = new TLine();

	for (int i_jet = 7; i_jet < 11; i_jet++)
	{
		for (int i_trk = 2; i_trk < 7; i_trk++)
		{
			c->cd();
			c->Clear();
			c->Divide(3,2);

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

				for (int i_dR = 1; i_dR <= h_systematic->GetXaxis()->GetNbins() ; i_dR++)
				{
					double max_delta = -999;
					for (int i_run = 0; i_run < f_run_UE.size(); i_run++)
					{
						double delta = abs(h_run_UE[i_run]->GetBinContent(i_dR)-h_nominal->GetBinContent(i_dR))  ;
						if (delta > max_delta)
						{
							max_delta = delta;

						}
					}
					h_systematic->SetBinContent(i_dR, max_delta);
				}


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
				double min = 0.5;
				double max = 1.5;

				if (i_cent == 0 && i_trk == 2 && i_jet == 7)
				{
					legend_MB->AddEntry(h_nominal,"Nominal", "lp");
					legend_MB->AddEntry(h_combined,"Combined", "lp");

					for (int i_run = 0; i_run < f_run_UE.size(); i_run++) legend_MB->AddEntry(h_run_UE[i_run],f_run_UE[i_run]->GetName(), "lp");
				}
				c->cd(i_cent+1);
				h_combined->GetXaxis()->SetRangeUser(0, r_max_range);
				h_combined->GetYaxis()->SetRangeUser(min, max);
				SetHStyle_smallify(h_nominal,1,1);
				SetHStyle_smallify(h_combined,0,1);

//				h_nominal->DrawCopy("");
				for (int i_run = 0; i_run < f_run_UE.size(); i_run++)
				{
					SetHStyle_smallify(h_run_UE[i_run],i_run+2,1);
					h_run_UE[i_run]->Divide(h_nominal);
//					h_run_UE[i_run]->DrawCopy("same");
				}
				h_combined->Divide(h_nominal);
				h_combined->Draw("same hist text");

				line->DrawLine(0,1,0.8,1);



//				for (int i_run = 0; i_run < f_run_UE.size(); i_run++)
//				{
//
//					h_run_UE[i_run] = (TH1*)f_run_UE[i_run]->Get(name.c_str());
//					SetHStyle_smallify(h_run_UE[i_run],i_run+1,1);
//					h_run_UE[i_run]->Draw("same hist");
//
//					if (i_cent == 0 && i_trk == 2 && i_jet == 7)
//					{
//						legend_MB->AddEntry(h_run_UE[i_run],Form("%s",f_run_UE[i_run]->GetName()) , "lp");
//					}
//
//				}
			}

			c->cd(1);
			legend_MB->Draw();
			c->cd();
			string pdf_label = "";
			if (i_jet == 7 && i_trk == 2) pdf_label = "(";
			if (i_jet == 10 && i_trk == 6) pdf_label = ")";
			c->Print(Form("weightedRuns.pdf%s",pdf_label.c_str()),Form("Title: jet%i_trk%i", i_jet, i_trk));

		}
	}


}
