#include "../functions/global_variables.h"


void iteration_stability()
{
	SetAtlasStyle();
	gErrorIgnoreLevel = 3001;


	vector<TFile*> files_data;
	vector<TFile*> files_mc;
	string name;
	for (int i = 1; i < 10; i++)
	{
		name = Form("iteration_test/final_ChPS_iter%i_MC_PbPb.root",i);
		files_mc.push_back(new TFile(name.c_str()));

		name = Form("iteration_test/final_ChPS_iter%i_data_PbPb.root",i);
		files_data.push_back(new TFile(name.c_str()));
	}


	TAxis* dR_binning = (TAxis*)files_mc[0]->Get("dR_binning");
	TAxis* jetpT_binning = (TAxis*)files_mc[0]->Get("jetpT_binning");
	TAxis* trkpT_binning = (TAxis*)files_mc[0]->Get("trkpT_binning");


	TCanvas *c = new TCanvas();
	TLegend *legend = new TLegend(0.48, 0.25, 0.90, 0.45, "","brNDC");
	legend->SetTextFont(43);
	legend->SetBorderSize(0);
	legend->SetNColumns(2);
	TLine *line = new TLine();




	c->Print("iteration_test/figures.pdf(","Title: Start");
	bool first_pass = true;

	for (int i_jet = 7; i_jet < 12; i_jet++)
	{
		string jet_label = Form("%1.0f < #it{p}_{T}^{jet} < %1.0f GeV", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));

		for (int i_trk = 2; i_trk < 9; i_trk++)
		{
			string trk_label = Form("%1.1f < #it{p}_{T}^{ch} < %1.1f GeV", trkpT_binning->GetBinLowEdge(i_trk+1), trkpT_binning->GetBinUpEdge(i_trk+1));

			c->Clear();
			c->Divide(3,2);
			for (int i_cent = 0; i_cent < 6; i_cent++)
			{
				string centrality = num_to_cent(31,i_cent);

//				name = Form("h_ChPS_truth_indR_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				name = Form("h_ChPS_raw_subtr_unf_indR_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				TH1* h_nominal = (TH1*)(files_data[3]->Get(name.c_str())->Clone(Form("%s_nom", name.c_str()) )) ;


				name = Form("h_ChPS_raw_subtr_unf_indR_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				c->cd(i_cent+1);

				for (int i_iter = 0; i_iter < files_mc.size(); i_iter++)
				{
					TH1* h_ChPS_data = (TH1*)files_data[i_iter]->Get(name.c_str());
					TH1* h_ChPS_mc = (TH1*)files_mc[i_iter]->Get(name.c_str());

					h_ChPS_data->Divide(h_nominal);
					h_ChPS_mc->Divide(h_nominal);

					h_ChPS_data->SetName(Form("%s_iter%i_data",name.c_str(), i_iter));
					h_ChPS_mc->SetName(Form("%s_iter%i_mc",name.c_str(), i_iter));

					SetHStyle(h_ChPS_data, i_iter, 1);
					SetHStyle(h_ChPS_mc, i_iter, 1);
					if (first_pass)
					{
						legend->AddEntry(h_ChPS_mc,Form("Iterations: %i", i_iter+1),"lp");
					}

					h_ChPS_mc->GetXaxis()->SetRangeUser(0., 0.8);
					h_ChPS_mc->GetYaxis()->SetRangeUser(0.6, 1.4);
					h_ChPS_data->GetXaxis()->SetRangeUser(0., 0.8);
					h_ChPS_data->GetYaxis()->SetRangeUser(0.6, 1.4);

					if (i_iter == 1) h_ChPS_data->DrawCopy("");
					else h_ChPS_data->DrawCopy("same");

//					if (i_iter == 1) h_ChPS_mc->DrawCopy("");
//					else h_ChPS_mc->DrawCopy("same");

				}
				first_pass = false;
				legend->Draw();
				line->DrawLine(0,1,0.8,1);
			}

			c->Print("iteration_test/figures.pdf",Form("Title: %s %s", trk_label.c_str(), jet_label.c_str()));
		}
	}
	c->Print("iteration_test/figures.pdf)","Title: End");

}
