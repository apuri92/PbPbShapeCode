//#include "combine_eff_dev.c"
#include "extras/global_variables.h"
#include "TVirtualFitter.h"
bool do_fine_eta = 1;
void draw_eff_trketa(bool isPbPb = 1, string cut = "ppTight")
{
	gErrorIgnoreLevel = 3001;
	gStyle->SetOptTitle(0);
	SetAtlasStyle();

	string name;
    string dataset;
    if (isPbPb) dataset = "PbPb";
    else dataset = "pp";
    
	//reading from input file
	name = Form("mc_efficiency_trketa_%s_%s.root",dataset.c_str(), cut.c_str());
	TFile *input_file = new TFile(name.c_str());

	name = Form("mc_eff_%s_trketa_jetptinc_%s.root",dataset.c_str(), cut.c_str());
	TFile *output_file = new TFile(name.c_str(),"recreate");
	cout << Form("%s -> %s", input_file->GetName(), output_file->GetName()) << endl;

	TAxis* trk_pt_binning = (TAxis*)input_file->Get("trk_pt_binning");
	TAxis* trk_eta_binning_coarse = (TAxis*)input_file->Get("trk_eta_binning_coarse");
    TAxis* trk_eta_binning_fine = (TAxis*)input_file->Get("trk_eta_binning_fine");

    int n_trk_pt_bins = trk_pt_binning->GetNbins();
	int n_trk_eta_bins_coarse = trk_eta_binning_coarse->GetNbins();
    int n_trk_eta_bins_fine = trk_eta_binning_fine->GetNbins();

    output_file->cd();
    if (do_fine_eta) trk_eta_binning_fine->Write("new_trk_eta_binning");
    
	TCanvas *canvas1 = new TCanvas("C1", "C1",900,600);

	TLine *line = new TLine();
	TLatex *ltx = new TLatex();
	ltx->SetNDC();
	ltx->SetTextFont(43);
	ltx->SetTextSize(15);
	ltx->SetTextAlign(11);

	TLegend *legend = new TLegend(0.21,0.20,0.70,0.40,NULL,"brNDC");
	legend->SetBorderSize(0);
	legend->SetNColumns(1);
	legend->SetTextFont(43);
	legend->SetTextSize(12);

    TAxis *eta_axis;
    if (do_fine_eta) eta_axis = (TAxis*)trk_eta_binning_fine->Clone("eta_axis");
    else eta_axis = (TAxis*)trk_eta_binning_coarse->Clone("eta_axis");
    int n_eta_bins = eta_axis->GetNbins();

    
	vector<vector<TH1*>> h_efficiency(n_cent_cuts+1, vector<TH1*> (n_eta_bins));
	vector<vector<TGraphAsymmErrors*>> g_efficiency(n_cent_cuts+1, vector<TGraphAsymmErrors*> (n_eta_bins));

    int start_eta = 0, end_eta = n_eta_bins-1;

    if (do_fine_eta)
    {
        start_eta = 7;
        end_eta = 14;
    }
	for (int i_cent_cuts = 0; i_cent_cuts < n_cent_cuts+1; i_cent_cuts++)
	{
		string centrality = num_to_cent(centrality_scheme,i_cent_cuts);
        
		for (int i_eta_cuts = 0; i_eta_cuts < n_eta_bins; i_eta_cuts++)
		{
			double eta_lo = eta_axis->GetBinLowEdge(i_eta_cuts+1);
			double eta_hi = eta_axis->GetBinUpEdge(i_eta_cuts+1);

			name = Form("histo_eff_eta%i_cent%i", i_eta_cuts, i_cent_cuts);
			h_efficiency.at(i_cent_cuts).at(i_eta_cuts) = (TH1*)input_file->Get(name.c_str());
            

			output_file->cd();

			g_efficiency.at(i_cent_cuts).at(i_eta_cuts) = new TGraphAsymmErrors(h_efficiency.at(i_cent_cuts).at(i_eta_cuts));
			name = Form("graph_eff_eta%i_cent%i",i_eta_cuts,i_cent_cuts);
			g_efficiency.at(i_cent_cuts).at(i_eta_cuts)->Write(name.c_str());

			name = Form("hist_eff_eta%i_cent%i",i_eta_cuts,i_cent_cuts);
			h_efficiency.at(i_cent_cuts).at(i_eta_cuts)->SetName(name.c_str());
			h_efficiency.at(i_cent_cuts).at(i_eta_cuts)->Write(name.c_str());


            if (i_eta_cuts < start_eta || i_eta_cuts > end_eta) continue;
            
            SetHStyle(h_efficiency.at(i_cent_cuts).at(i_eta_cuts),i_eta_cuts-start_eta);
//            smallify(h_efficiency.at(i_cent_cuts).at(i_eta_cuts));
            
            h_efficiency.at(i_cent_cuts).at(i_eta_cuts)->GetYaxis()->SetRangeUser(0.55,1.1);
            h_efficiency.at(i_cent_cuts).at(i_eta_cuts)->GetXaxis()->SetRangeUser(0.1,200);
            h_efficiency.at(i_cent_cuts).at(i_eta_cuts)->GetYaxis()->SetTitle("Efficiency");
            h_efficiency.at(i_cent_cuts).at(i_eta_cuts)->GetXaxis()->SetTitle("#it{p}_{T}^{truth} [GeV]");
            h_efficiency.at(i_cent_cuts).at(i_eta_cuts)->SetTitle(Form("Efficiency: %s, %4.2f < #eta < %4.2f",centrality.c_str(), eta_lo, eta_hi));
            
            if (i_cent_cuts == 0) legend->AddEntry(h_efficiency.at(i_cent_cuts).at(i_eta_cuts),Form("%4.2f < #eta < %4.2f",eta_lo, eta_hi),"lp");
            if (i_eta_cuts == start_eta) h_efficiency.at(i_cent_cuts).at(i_eta_cuts)->Draw("p");
            else h_efficiency.at(i_cent_cuts).at(i_eta_cuts)->Draw("same p");
            gPad->SetLogx();

            
		}
		name = Form("%s",centrality.c_str());
		ltx->SetTextAlign(11);
		ltx->DrawLatex(0.19,0.88,name.c_str());
        legend->Draw();

        
        if (i_cent_cuts == 0) name = "(";
        else if (i_cent_cuts == n_cent_cuts ) name = ")";
        else name = "";
        name = Form("eff_cent_trketa_%s_%s.pdf%s", dataset.c_str(), cut.c_str(), name.c_str());
        canvas1->Print(name.c_str());

	}


}

