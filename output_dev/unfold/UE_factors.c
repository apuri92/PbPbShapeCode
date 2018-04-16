#include "../functions/global_variables.h"

void UE_factors(string config_file = "ff_config.cfg")
{
	cout << "######### GETTING UE FACTORS #########" << endl;

	SetAtlasStyle();
	gErrorIgnoreLevel = 3001;

	//	##############	Reading config	##############"
	TEnv *m_config = new TEnv();
	m_config->ReadFile(config_file.c_str(), EEnvLevel(1));
	int sys_mode = -1; sys_mode = m_config->GetValue("sys_mode", sys_mode);
	int verbose = 0; verbose = m_config->GetValue("verbose", verbose);
	//	##############	Config done	##############"

	std::string sys_path = "";
	if (sys_mode > 0) sys_path = Form("systematics/%i", sys_mode);
	if (verbose) m_config->Print();


	TFile *input_file = new TFile(Form("../raw_results/%s/FF_MC_out_histo_PbPb_5p02_r001.root", sys_path.c_str()));
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
			name = Form("ChPS_MB_UE_dR%i_cent%i", i_dR, i_cent);
			TH2* h_MB_method = (TH2*)input_file->Get(name.c_str());
			name = Form("ChPS_MB_UE_err_dR%i_cent%i", i_dR, i_cent);
			TH2* h_MB_method_err = (TH2*)input_file->Get(name.c_str());


			name = Form("ChPS_TM_UE_dR%i_cent%i", i_dR, i_cent);
			TH2* h_TM_method = (TH2*)input_file->Get(name.c_str());


			//SET ERRORS FROM MB_ERR HISTOGRAM
			for (int i_jet_bin = 1; i_jet_bin <= N_jetpt; i_jet_bin++)
			{
				for (int i_trk_bin = 1; i_trk_bin <= N_jetpt; i_trk_bin++)
				{
					double updated_UE_MB_err = h_MB_method_err->GetBinContent(i_trk_bin, i_jet_bin);
					h_MB_method->SetBinError(i_trk_bin, i_jet_bin, updated_UE_MB_err);
				}
			}


			name = Form("ratio_dR%i_cent%i", i_dR, i_cent);
			h_ratio[i_dR][i_cent] = (TH2*)h_TM_method->Clone(name.c_str());
			h_ratio[i_dR][i_cent]->Divide(h_MB_method); //errors propagated assuming they are uncorrelated

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
		c1->Print(Form("output_pdf/PbPb/UE_factors.pdf%s", name.c_str()), Form("Title: dR%i - %s", i_dR, dr_label.c_str()));
	}

	cout << "######### DONE GETTING UE FACTORS #########" << endl << endl;;

}
