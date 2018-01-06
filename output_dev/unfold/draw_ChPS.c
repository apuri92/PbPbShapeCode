#include "../extras/global_variables.h"
#include "TEnv.h"
#include "TGaxis.h"

void draw_ChPS()
{
	SetAtlasStyle();
	gErrorIgnoreLevel = 3001;

	cout << "Drawing ChPS..." << endl;
	//	##############	Reading config	##############"
	TEnv *m_config = new TEnv();
	m_config->ReadFile("ff_config.cfg", EEnvLevel(1));
	m_config->Print();

	std::string dataset_type = "PbPb"; dataset_type = m_config->GetValue("dataset_type", dataset_type.c_str());
	std::string tracking_cut = "ppTight"; tracking_cut = m_config->GetValue("tracking_cut", tracking_cut.c_str());
	int centrality_scheme = 31; centrality_scheme = m_config->GetValue("centrality_scheme", centrality_scheme);
	int isMC = 1; isMC = m_config->GetValue("isMC", isMC);
	int n_unfold = 4; n_unfold = m_config->GetValue("n_unfold", n_unfold);

	std::string did = "data";
	if (isMC) did = "MC_JZ_comb";

	//	##############	Config done	##############"

	TFile *f_input = new TFile(Form("unfolded_%s_%s.root",did.c_str(), dataset_type.c_str()),"recreate");
    //Configuration done

	TAxis* dR_binning = (TAxis*)f_input->Get("dR_binning");
	TAxis* jetpT_binning = (TAxis*)f_input->Get("jetpT_binning");
	TAxis* trkpT_binning = (TAxis*)f_input->Get("trkpT_binning");

	int N_dR = dR_binning->GetNbins();
	int N_jetpt = jetpT_binning->GetNbins();
	int N_trkpt = trkpT_binning->GetNbins();

	//raw_0
	vector<vector<vector<TH1*>>> h_ChPS_raw (n_cent_cuts, vector<vector<TH1*>> (N_dR, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_raw_subtr (n_cent_cuts, vector<vector<TH1*>> (N_dR, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_raw_subtr_unf (n_cent_cuts, vector<vector<TH1*>> (N_dR, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_raw_subtr_unf_bbb (n_cent_cuts, vector<vector<TH1*>> (N_dR, vector<TH1*> (N_jetpt)));

	vector<vector<vector<TH1*>>> h_ChPS_ratio_subtr_raw (n_cent_cuts, vector<vector<TH1*>> (N_dR, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_ratio_unf_subtr (n_cent_cuts, vector<vector<TH1*>> (N_dR, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_ratio_closure (n_cent_cuts, vector<vector<TH1*>> (N_dR, vector<TH1*> (N_jetpt)));

	vector<vector<vector<TH1*>>> h_ChPS_diff_subtr_raw (n_cent_cuts, vector<vector<TH1*>> (N_dR, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_diff_unf_subtr (n_cent_cuts, vector<vector<TH1*>> (N_dR, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_diff_closure (n_cent_cuts, vector<vector<TH1*>> (N_dR, vector<TH1*> (N_jetpt)));

	//raw_rr
	vector<vector<vector<TH1*>>> h_ChPS_raw_rr (n_cent_cuts, vector<vector<TH1*>> (N_dR, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_raw_rr_unf (n_cent_cuts, vector<vector<TH1*>> (N_dR, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_raw_subtr_unf_bbb (n_cent_cuts, vector<vector<TH1*>> (N_dR, vector<TH1*> (N_jetpt)));

	vector<vector<vector<TH1*>>> h_ChPS_ratio_unf_raw_rr (n_cent_cuts, vector<vector<TH1*>> (N_dR, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_ratio_closure_rr (n_cent_cuts, vector<vector<TH1*>> (N_dR, vector<TH1*> (N_jetpt)));

	vector<vector<vector<TH1*>>> h_ChPS_diff_unf_raw_rr (n_cent_cuts, vector<vector<TH1*>> (N_dR, vector<TH1*> (N_jetpt)));
	vector<vector<vector<TH1*>>> h_ChPS_diff_closure_rr (n_cent_cuts, vector<vector<TH1*>> (N_dR, vector<TH1*> (N_jetpt)));

	//truth
	vector<vector<vector<TH1*>>> h_ChPS_truth (n_cent_cuts, vector<vector<TH1*>> (N_dR, vector<TH1*> (N_jetpt)));

	//UE
	vector<vector<vector<TH1*>>> h_ChPS_UE (n_cent_cuts, vector<vector<TH1*>> (N_dR, vector<TH1*> (N_jetpt)));


	/**/

//	vector<vector<TH1*>> h_ChPS_raw_injet (n_cent_cuts, vector<TH1*> (N_jetpt));
//	vector<vector<TH1*>> h_ChPS_UE_injet (n_cent_cuts, vector<TH1*> (N_jetpt));
//	vector<vector<TH1*>> h_ChPS_UE_SB_injet (n_cent_cuts, vector<TH1*> (N_jetpt));

	string name;

	TCanvas *c1 = new TCanvas("c1","c1",900,600);
//	TCanvas *c2 = new TCanvas("c2","c2",900,600);
//	TCanvas *c3 = new TCanvas("c3","c3",900,600);
//	TCanvas *c4 = new TCanvas("c4","c4",900,600);
//	TCanvas *c5 = new TCanvas("c5","c5",900,600);
//	TCanvas *c6 = new TCanvas("c6","c6",900,600);
//	TCanvas *c7 = new TCanvas("c7","c7",900,600);
//	TCanvas *c8 = new TCanvas("c8","c8",900,600);

	TLine *line = new TLine();
	line->SetLineColor(kBlack);

	TLegend *legend = new TLegend(0.20, 0.20, 0.48, 0.40, "","brNDC");
	legend->SetTextFont(43);
	legend->SetBorderSize(0);
	legend->SetTextSize(12);

	TLatex *ltx = new TLatex();
	ltx->SetTextFont(43);
	ltx->SetTextSize(16);

	double max = 0, min = 0;

	bool doing_injet = true;
	for (int i_dR = 0; i_dR < N_dR; i_dR++)
	{
		double dR_lo = dR_binning->GetBinLowEdge(i_dR+1);
		double dR_hi = dR_binning->GetBinUpEdge(i_dR+1);
		double area = TMath::Pi() * ((dR_hi*dR_hi) - (dR_lo*dR_lo));
		double area_injet = TMath::Pi() * (0.4*0.4);

		c1->Divide(3,2);
//		c8->Divide(3,2);
//
//		c2->Divide(3,2);
//		c7->Divide(3,2);
//		c3->Divide(3,2);
//		if (doing_injet) c4->Divide(3,2);
//		if (doing_injet) c6->Divide(3,2);
//		c5->Divide(3,2);

		for (int i_cent_cuts = 0; i_cent_cuts < n_cent_cuts; i_cent_cuts++)
		{
//            if (i_cent_cuts < 5) continue;

//            if (dataset_type == "pp" && i_cent_cuts !=5) continue;
			string centrality = num_to_cent(centrality_scheme,i_cent_cuts);
			name = Form("h_%s_cent%i_dR%i",ChPS_raw_type.c_str(), i_cent_cuts, i_dR);
            cout << name << endl;
            h_raw.at(i_cent_cuts).at(i_dR) = (TH2*)f_input->Get(name.c_str())->Clone(Form("c_%s",name.c_str()));
			h_raw.at(i_cent_cuts).at(i_dR)->SetName(name.c_str());
			h_raw.at(i_cent_cuts).at(i_dR)->Sumw2();

            
			name = Form("h_true_ChPS_cent%i_dR%i",i_cent_cuts, i_dR);
			h_true.at(i_cent_cuts).at(i_dR) = (TH2*)f_input->Get(name.c_str())->Clone(Form("c_%s",name.c_str()));
			h_true.at(i_cent_cuts).at(i_dR)->SetName(name.c_str());
			h_true.at(i_cent_cuts).at(i_dR)->Sumw2();

			name = Form("h_unfolded_%s_cent%i_dR%i", ChPS_raw_type.c_str(), i_cent_cuts, i_dR);
			h_unfolded.at(i_cent_cuts).at(i_dR) = (TH2*)f_input->Get(name.c_str())->Clone(Form("c_%s",name.c_str()));
			h_unfolded.at(i_cent_cuts).at(i_dR)->SetName(name.c_str());
			h_unfolded.at(i_cent_cuts).at(i_dR)->Sumw2();

			name = Form("h_UE_cent%i_dR%i", i_cent_cuts, i_dR);
			h_UE.at(i_cent_cuts).at(i_dR) = (TH2*)f_input->Get(name.c_str())->Clone(Form("c_%s",name.c_str()));
			h_UE.at(i_cent_cuts).at(i_dR)->SetName(name.c_str());
			h_UE.at(i_cent_cuts).at(i_dR)->Sumw2();

			if (doing_injet)
			{
				name = Form("h_%s_injet_cent%i",ChPS_raw_type.c_str(), i_cent_cuts);
				h_raw_injet.at(i_cent_cuts) = (TH2*)f_input->Get(name.c_str())->Clone(Form("c_%s",name.c_str()));
				h_raw_injet.at(i_cent_cuts)->SetName(name.c_str());
				h_raw_injet.at(i_cent_cuts)->Sumw2();

				name = Form("h_UE_injet_cent%i", i_cent_cuts);
				h_UE_injet.at(i_cent_cuts) = (TH2*)f_input->Get(name.c_str())->Clone(Form("c_%s",name.c_str()));
				h_UE_injet.at(i_cent_cuts)->SetName(name.c_str());
				h_UE_injet.at(i_cent_cuts)->Sumw2();
			}

			int jet_pt_start = 7;
			int jet_pt_end = 10;
			int style = 0;

			c1->cd(i_cent_cuts+1);
			gPad->Divide(1,2);

			c2->cd(i_cent_cuts+1);
			gPad->Divide(1,2);

			c8->cd(i_cent_cuts+1);

			for (int i_jet_bin = 0; i_jet_bin < N_jetpt; i_jet_bin++)
			{
				name = Form("h_ChPS_raw_cent%i_dR%i_jetpt%i",i_cent_cuts, i_dR, i_jet_bin);
				h_ChPS_raw.at(i_cent_cuts).at(i_dR).at(i_jet_bin) = (TH1*)h_raw.at(i_cent_cuts).at(i_dR)->ProjectionX(name.c_str(), i_jet_bin+1,i_jet_bin+1);
				h_ChPS_raw.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->Scale(1.,"width");
				h_ChPS_raw.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->Scale(1./area);
				h_ChPS_raw.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->Sumw2();

				name = Form("h_ChPS_unf_cent%i_dR%i_jetpt%i",i_cent_cuts, i_dR, i_jet_bin);
				h_ChPS_unf.at(i_cent_cuts).at(i_dR).at(i_jet_bin) = (TH1*)h_unfolded.at(i_cent_cuts).at(i_dR)->ProjectionX(name.c_str(), i_jet_bin+1,i_jet_bin+1);
				h_ChPS_unf.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->Scale(1.,"width");
				h_ChPS_unf.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->Scale(1./area);
				h_ChPS_unf.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->Sumw2();

				name = Form("h_ChPS_tru_cent%i_dR%i_jetpt%i",i_cent_cuts, i_dR, i_jet_bin);
				h_ChPS_tru.at(i_cent_cuts).at(i_dR).at(i_jet_bin) = (TH1*)h_true.at(i_cent_cuts).at(i_dR)->ProjectionX(name.c_str(), i_jet_bin+1, i_jet_bin+1);
				h_ChPS_tru.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->Scale(1.,"width");
				h_ChPS_tru.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->Scale(1./area);
				h_ChPS_tru.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->Sumw2();

				name = Form("h_ChPS_UE_cent%i_dR%i_jetpt%i",i_cent_cuts, i_dR, i_jet_bin);
				h_ChPS_UE.at(i_cent_cuts).at(i_dR).at(i_jet_bin) = (TH1*)h_UE.at(i_cent_cuts).at(i_dR)->ProjectionX(name.c_str(), i_jet_bin+1, i_jet_bin+1);
				h_ChPS_UE.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->Scale(1.,"width");
				h_ChPS_UE.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->Scale(1./area);
				h_ChPS_UE.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->Sumw2();

				name = Form("h_ChPS_closure_cent%i_dR%i_jetpt%i",i_cent_cuts, i_dR, i_jet_bin);
				h_ChPS_closure.at(i_cent_cuts).at(i_dR).at(i_jet_bin) = (TH1*)h_ChPS_unf.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->Clone(name.c_str());
				h_ChPS_closure.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->Divide(h_ChPS_tru.at(i_cent_cuts).at(i_dR).at(i_jet_bin));

				name = Form("h_ChPS_closure_diff_cent%i_dR%i_jetpt%i",i_cent_cuts, i_dR, i_jet_bin);
				h_ChPS_closure_diff.at(i_cent_cuts).at(i_dR).at(i_jet_bin) = (TH1*)h_ChPS_unf.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->Clone(name.c_str());
				h_ChPS_closure_diff.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->Add(h_ChPS_tru.at(i_cent_cuts).at(i_dR).at(i_jet_bin), -1);

				name = Form("h_ChPS_ratio_cent%i_dR%i_jetpt%i",i_cent_cuts, i_dR, i_jet_bin);
				h_ChPS_ratio.at(i_cent_cuts).at(i_dR).at(i_jet_bin) = (TH1*)h_ChPS_unf.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->Clone(name.c_str());
				h_ChPS_ratio.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->Divide(h_ChPS_raw.at(i_cent_cuts).at(i_dR).at(i_jet_bin));

				name = Form("h_ChPS_UE_SB_cent%i_dR%i_jetpt%i",i_cent_cuts, i_dR, i_jet_bin);
				h_ChPS_UE_SB.at(i_cent_cuts).at(i_dR).at(i_jet_bin) = (TH1*)h_ChPS_UE.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->Clone(name.c_str());
				h_ChPS_UE_SB.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->Divide(h_ChPS_raw.at(i_cent_cuts).at(i_dR).at(i_jet_bin));

				if (doing_injet)
				{
					name = Form("h_ChPS_raw_injet_cent%i_jetpt%i",i_cent_cuts, i_jet_bin);
					h_ChPS_raw_injet.at(i_cent_cuts).at(i_jet_bin) = (TH1*)h_raw_injet.at(i_cent_cuts)->ProjectionX(name.c_str(), i_jet_bin+1,i_jet_bin+1);
					h_ChPS_raw_injet.at(i_cent_cuts).at(i_jet_bin)->Scale(1.,"width");
					h_ChPS_raw_injet.at(i_cent_cuts).at(i_jet_bin)->Scale(1./area_injet);
					h_ChPS_raw_injet.at(i_cent_cuts).at(i_jet_bin)->Sumw2();

					name = Form("h_ChPS_UE_injet_cent%i_jetpt%i",i_cent_cuts, i_jet_bin);
					h_ChPS_UE_injet.at(i_cent_cuts).at(i_jet_bin) = (TH1*)h_UE_injet.at(i_cent_cuts)->ProjectionX(name.c_str(), i_jet_bin+1, i_jet_bin+1);
					h_ChPS_UE_injet.at(i_cent_cuts).at(i_jet_bin)->Scale(1.,"width");
					h_ChPS_UE_injet.at(i_cent_cuts).at(i_jet_bin)->Scale(1./area_injet);
					h_ChPS_UE_injet.at(i_cent_cuts).at(i_jet_bin)->Sumw2();


					name = Form("h_ChPS_UE_SB_injet_cent%i_jetpt%i",i_cent_cuts, i_jet_bin);
					h_ChPS_UE_SB_injet.at(i_cent_cuts).at(i_jet_bin) = (TH1*)h_ChPS_UE_injet.at(i_cent_cuts).at(i_jet_bin)->Clone(name.c_str());
					h_ChPS_UE_SB_injet.at(i_cent_cuts).at(i_jet_bin)->Divide(h_ChPS_raw_injet.at(i_cent_cuts).at(i_jet_bin));
				}

				f_output->cd();

				name = Form("h_ChPS_raw_cent%i_dR%i_jetpt%i",i_cent_cuts, i_dR, i_jet_bin);
				h_ChPS_raw.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->SetName(name.c_str());
				h_ChPS_raw.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->Write(name.c_str());

				name = Form("h_ChPS_tru_cent%i_dR%i_jetpt%i",i_cent_cuts, i_dR, i_jet_bin);
				h_ChPS_tru.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->SetName(name.c_str());
				h_ChPS_tru.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->Write(name.c_str());

				name = Form("h_ChPS_unf_cent%i_dR%i_jetpt%i",i_cent_cuts, i_dR, i_jet_bin);
				h_ChPS_unf.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->SetName(name.c_str());
				h_ChPS_unf.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->Write(name.c_str());

				name = Form("h_ChPS_closure_cent%i_dR%i_jetpt%i",i_cent_cuts, i_dR, i_jet_bin);
				h_ChPS_closure.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->SetName(name.c_str());
				h_ChPS_closure.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->Write(name.c_str());

				name = Form("h_ChPS_closure_diff_cent%i_dR%i_jetpt%i",i_cent_cuts, i_dR, i_jet_bin);
				h_ChPS_closure_diff.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->SetName(name.c_str());
				h_ChPS_closure_diff.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->Write(name.c_str());

				name = Form("h_ChPS_ratio_cent%i_dR%i_jetpt%i",i_cent_cuts, i_dR, i_jet_bin);
				h_ChPS_ratio.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->SetName(name.c_str());
				h_ChPS_ratio.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->Write(name.c_str());

				name = Form("h_ChPS_UE_cent%i_dR%i_jetpt%i",i_cent_cuts, i_dR, i_jet_bin);
				h_ChPS_UE.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->SetName(name.c_str());
				h_ChPS_UE.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->Write(name.c_str());

				name = Form("h_ChPS_UE_SB_cent%i_dR%i_jetpt%i",i_cent_cuts, i_dR, i_jet_bin);
				h_ChPS_UE_SB.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->SetName(name.c_str());
				h_ChPS_UE_SB.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->Write(name.c_str());

				if (doing_injet)
				{
					name = Form("h_ChPS_UE_SB_injet_cent%i_dR%i_jetpt%i",i_cent_cuts, i_dR, i_jet_bin);
					h_ChPS_UE_injet.at(i_cent_cuts).at(i_jet_bin)->SetName(name.c_str());
					h_ChPS_UE_injet.at(i_cent_cuts).at(i_jet_bin)->Write(name.c_str());

					name = Form("h_ChPS_UE_SB_injet_cent%i_dR%i_jetpt%i",i_cent_cuts, i_dR, i_jet_bin);
					h_ChPS_UE_SB_injet.at(i_cent_cuts).at(i_jet_bin)->SetName(name.c_str());
					h_ChPS_UE_SB_injet.at(i_cent_cuts).at(i_jet_bin)->Write(name.c_str());
				}

				if (i_jet_bin < jet_pt_start || i_jet_bin > jet_pt_end) continue;

				SetHStyle(h_ChPS_raw.at(i_cent_cuts).at(i_dR).at(i_jet_bin),style);
				SetHStyle(h_ChPS_unf.at(i_cent_cuts).at(i_dR).at(i_jet_bin),style);
				SetHStyle(h_ChPS_tru.at(i_cent_cuts).at(i_dR).at(i_jet_bin),style);
				SetHStyle(h_ChPS_closure.at(i_cent_cuts).at(i_dR).at(i_jet_bin),style);
				SetHStyle(h_ChPS_closure_diff.at(i_cent_cuts).at(i_dR).at(i_jet_bin),style);
				SetHStyle(h_ChPS_ratio.at(i_cent_cuts).at(i_dR).at(i_jet_bin),style);
				SetHStyle(h_ChPS_UE.at(i_cent_cuts).at(i_dR).at(i_jet_bin),style);
				SetHStyle(h_ChPS_UE_SB.at(i_cent_cuts).at(i_dR).at(i_jet_bin),style);
				if (doing_injet) SetHStyle(h_ChPS_UE_injet.at(i_cent_cuts).at(i_jet_bin),style);
				if (doing_injet) SetHStyle(h_ChPS_UE_SB_injet.at(i_cent_cuts).at(i_jet_bin),style);

				smallify(h_ChPS_raw.at(i_cent_cuts).at(i_dR).at(i_jet_bin));
				smallify(h_ChPS_unf.at(i_cent_cuts).at(i_dR).at(i_jet_bin));
				smallify(h_ChPS_tru.at(i_cent_cuts).at(i_dR).at(i_jet_bin));
				smallify(h_ChPS_closure.at(i_cent_cuts).at(i_dR).at(i_jet_bin));
				smallify(h_ChPS_closure_diff.at(i_cent_cuts).at(i_dR).at(i_jet_bin));
				smallify(h_ChPS_ratio.at(i_cent_cuts).at(i_dR).at(i_jet_bin));
				smallify(h_ChPS_UE.at(i_cent_cuts).at(i_dR).at(i_jet_bin));
				smallify(h_ChPS_UE_SB.at(i_cent_cuts).at(i_dR).at(i_jet_bin));
				if (doing_injet) smallify(h_ChPS_UE_injet.at(i_cent_cuts).at(i_jet_bin));
				if (doing_injet) smallify(h_ChPS_UE_SB_injet.at(i_cent_cuts).at(i_jet_bin));


				double range_pt_lo = 1;
				double range_pt_hi = 158.;
				if (dR_lo >= 0.5) range_pt_hi = 60.;
				if (dR_lo >= 0.7) range_pt_hi = 20.;

				h_ChPS_raw.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->GetXaxis()->SetRangeUser(range_pt_lo, range_pt_hi);
				h_ChPS_unf.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->GetXaxis()->SetRangeUser(range_pt_lo, range_pt_hi);
				h_ChPS_tru.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->GetXaxis()->SetRangeUser(range_pt_lo, range_pt_hi);
				h_ChPS_ratio.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->GetXaxis()->SetRangeUser(range_pt_lo, range_pt_hi);
				h_ChPS_closure.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->GetXaxis()->SetRangeUser(range_pt_lo, range_pt_hi);
				h_ChPS_closure_diff.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->GetXaxis()->SetRangeUser(range_pt_lo, range_pt_hi);
				h_ChPS_UE.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->GetXaxis()->SetRangeUser(range_pt_lo, range_pt_hi);
				h_ChPS_UE_injet.at(i_cent_cuts).at(i_jet_bin)->GetXaxis()->SetRangeUser(range_pt_lo, range_pt_hi);
				h_ChPS_UE_SB.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->GetXaxis()->SetRangeUser(range_pt_lo, range_pt_hi);
				h_ChPS_UE_SB_injet.at(i_cent_cuts).at(i_jet_bin)->GetXaxis()->SetRangeUser(range_pt_lo, range_pt_hi);

				// *** C1 *** //
				h_ChPS_raw.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]");
				h_ChPS_raw.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->GetYaxis()->SetTitle("#frac{1}{Area N_{Jets}} #frac{dN}{dp_{T}}");
				h_ChPS_tru.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]");
				h_ChPS_tru.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->GetYaxis()->SetTitle("#frac{1}{Area N_{Jets}} #frac{dN}{dp_{T}}");


				h_ChPS_raw.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->GetYaxis()->SetNdivisions(504);
				h_ChPS_unf.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->GetYaxis()->SetNdivisions(504);
				h_ChPS_tru.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->GetYaxis()->SetNdivisions(504);
				h_ChPS_closure.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->GetYaxis()->SetNdivisions(504);
				h_ChPS_closure_diff.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->GetYaxis()->SetNdivisions(504);

				c1->cd(i_cent_cuts+1)->cd(1);
				gPad->SetPad(0,0.51,1,.95);
				gPad->SetBottomMargin(0);
				gPad->SetRightMargin(0);
				if (i_jet_bin == jet_pt_start) h_ChPS_raw.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->Draw();
				else h_ChPS_raw.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->Draw("same");
				gPad->SetLogx();
				gPad->SetLogy();

				c1->cd(i_cent_cuts+1)->cd(2);
				gPad->SetPad(0, 0.0, 1, 0.51);
				gPad->SetTopMargin(0);
				gPad->SetRightMargin(0);
				if (i_jet_bin == jet_pt_start) h_ChPS_tru.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->Draw();
				else h_ChPS_tru.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->Draw("same");
				gPad->SetLogx();
				gPad->SetLogy();
				// ***

				// *** C8 *** //
				h_ChPS_unf.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]");
				h_ChPS_unf.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->GetYaxis()->SetTitle("#frac{1}{Area N_{Jets}} #frac{dN}{dp_{T}}");

				h_ChPS_unf.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->GetYaxis()->SetNdivisions(504);

				c8->cd(i_cent_cuts+1);
				if (i_jet_bin == jet_pt_start) h_ChPS_unf.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->Draw();
				else h_ChPS_unf.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->Draw("same");
				gPad->SetLogx();
				gPad->SetLogy();
				// ***


				// *** C2 *** //
				h_ChPS_closure.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]");
				h_ChPS_closure.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->GetYaxis()->SetRangeUser(0.45,1.55);

				h_ChPS_closure_diff.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]");
				h_ChPS_closure_diff.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->GetYaxis()->SetRangeUser(-0.15,0.15);

				name = "Closure"; //default title
				if (ChPS_raw_type == "ChPS_raw_tt") name = "Eff. Corr. Closure";
				if (ChPS_raw_type == "ChPS_raw_0") name = "Unf/Truth";
				h_ChPS_closure.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->GetYaxis()->SetTitle(name.c_str());
				h_ChPS_closure_diff.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->GetYaxis()->SetTitle("Unf - Truth");

				c2->cd(i_cent_cuts+1);
				c2->cd(i_cent_cuts+1)->cd(1);
				gPad->SetPad(0,0.51,1,.95);
				gPad->SetBottomMargin(0);
				gPad->SetRightMargin(0);
				if (i_jet_bin == jet_pt_start) h_ChPS_closure.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->Draw();
				else h_ChPS_closure.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->Draw("same");
				line->SetLineStyle(3);
				line->DrawLine(1,1.02, 150, 1.02);
				line->DrawLine(1,0.98, 150, 0.98);

				line->SetLineStyle(2);
				line->DrawLine(1,1, 150, 1);
				gPad->SetLogx();
				gPad->SetLogy(0);



				c2->cd(i_cent_cuts+1)->cd(2);
				gPad->SetPad(0, 0.0, 1, 0.51);
				gPad->SetTopMargin(0);
				gPad->SetRightMargin(0);
				if (i_jet_bin == jet_pt_start) h_ChPS_closure_diff.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->Draw();
				else h_ChPS_closure_diff.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->Draw("same");
				line->DrawLine(1,0,158,0);
				gPad->SetLogx();
				gPad->SetLogy(0);
				// ***


				// *** C7 *** //
				h_ChPS_ratio.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]");
				h_ChPS_ratio.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->GetYaxis()->SetRangeUser(0,2);
				name = "(Unf+BbB)/(Raw-UE)"; //default title
				if (ChPS_raw_type == "ChPS_raw_tt") name = "Eff. Corr. Closure";
				if (ChPS_raw_type == "ChPS_raw_0" && do_UE_subtr && do_unfolding && !do_BbB) name = "Unf/(Raw-UE)";
				h_ChPS_ratio.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->GetYaxis()->SetTitle(name.c_str());

				c7->cd(i_cent_cuts+1);
				if (i_jet_bin == jet_pt_start) h_ChPS_ratio.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->Draw();
				else h_ChPS_ratio.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->Draw("same");
				gPad->SetLogx();
				gPad->SetLogy(0);
				// ***

				// *** C3 *** //
				h_ChPS_UE.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]");
				h_ChPS_UE.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->GetYaxis()->SetTitle("#frac{1}{N_{jets}^{raw} Area} #frac{dN_{UE}}{dp_{T}}");
				h_ChPS_UE.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->GetYaxis()->SetTitleFont(43);
				h_ChPS_UE.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->GetYaxis()->SetTitleSize(12);
				h_ChPS_UE.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->GetYaxis()->SetTitleOffset(3);

				max = h_ChPS_UE.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->GetMaximum();
				min = h_ChPS_UE.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->GetMinimum();
				h_ChPS_UE.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->GetYaxis()->SetRangeUser(0.5*min,1.5*max);
				h_ChPS_UE.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->GetYaxis()->SetNdivisions(504);

				c3->cd(i_cent_cuts+1);
				if (i_jet_bin == jet_pt_start) h_ChPS_UE.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->Draw();
				else h_ChPS_UE.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->Draw("same");
				gPad->SetLogx();
				gPad->SetLogy(0);
				// ***

				// *** C5 *** //
				h_ChPS_UE_SB.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]");
				h_ChPS_UE_SB.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->GetYaxis()->SetTitle("#frac{dN_{UE}}{dN_{ChPS}^{Subtr}}");
				h_ChPS_UE_SB.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->GetYaxis()->SetTitleFont(43);
				h_ChPS_UE_SB.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->GetYaxis()->SetTitleSize(12);
				h_ChPS_UE_SB.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->GetYaxis()->SetTitleOffset(3);
				h_ChPS_UE_SB.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->GetYaxis()->SetNdivisions(504);

				c5->cd(i_cent_cuts+1);
				if (i_jet_bin == jet_pt_start) h_ChPS_UE_SB.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->Draw();
				else h_ChPS_UE_SB.at(i_cent_cuts).at(i_dR).at(i_jet_bin)->Draw("same");
				gPad->SetLogx();
				gPad->SetLogy(0);
				// ***

				// *** C4 *** //
				if (doing_injet)
				{
					h_ChPS_UE_injet.at(i_cent_cuts).at(i_jet_bin)->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]");
					h_ChPS_UE_injet.at(i_cent_cuts).at(i_jet_bin)->GetYaxis()->SetTitle("#frac{1}{N_{jets}^{raw} Area^{Jet}} #frac{dN_{UE}}{dp_{T}} (In jet)");
					h_ChPS_UE_injet.at(i_cent_cuts).at(i_jet_bin)->GetYaxis()->SetTitleFont(43);
					h_ChPS_UE_injet.at(i_cent_cuts).at(i_jet_bin)->GetYaxis()->SetTitleSize(12);
					h_ChPS_UE_injet.at(i_cent_cuts).at(i_jet_bin)->GetYaxis()->SetTitleOffset(3);
					max = h_ChPS_UE_injet.at(i_cent_cuts).at(i_jet_bin)->GetMaximum();
					min = h_ChPS_UE_injet.at(i_cent_cuts).at(i_jet_bin)->GetMinimum();
					h_ChPS_UE_injet.at(i_cent_cuts).at(i_jet_bin)->GetYaxis()->SetRangeUser(0,1.5*max);
					h_ChPS_UE_injet.at(i_cent_cuts).at(i_jet_bin)->GetYaxis()->SetNdivisions(504);

					h_ChPS_UE_SB_injet.at(i_cent_cuts).at(i_jet_bin)->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]");
					h_ChPS_UE_SB_injet.at(i_cent_cuts).at(i_jet_bin)->GetYaxis()->SetTitle("#frac{dN_{UE}}{dN_{ChPS}^{Subtr}} (In jet)");
					h_ChPS_UE_SB_injet.at(i_cent_cuts).at(i_jet_bin)->GetYaxis()->SetTitleFont(43);
					h_ChPS_UE_SB_injet.at(i_cent_cuts).at(i_jet_bin)->GetYaxis()->SetTitleSize(12);
					h_ChPS_UE_SB_injet.at(i_cent_cuts).at(i_jet_bin)->GetYaxis()->SetTitleOffset(3);
					max = h_ChPS_UE_SB_injet.at(i_cent_cuts).at(i_jet_bin)->GetMaximum();
					min = h_ChPS_UE_SB_injet.at(i_cent_cuts).at(i_jet_bin)->GetMinimum();
					h_ChPS_UE_SB_injet.at(i_cent_cuts).at(i_jet_bin)->GetYaxis()->SetRangeUser(0,1.5*max);
					h_ChPS_UE_SB_injet.at(i_cent_cuts).at(i_jet_bin)->GetYaxis()->SetNdivisions(504);


					c4->cd(i_cent_cuts+1);
					if (i_jet_bin == jet_pt_start) h_ChPS_UE_injet.at(i_cent_cuts).at(i_jet_bin)->Draw();
					else h_ChPS_UE_injet.at(i_cent_cuts).at(i_jet_bin)->Draw("same");
					gPad->SetLogx();
					gPad->SetLogy(0);

					c6->cd(i_cent_cuts+1);
					if (i_jet_bin == jet_pt_start) h_ChPS_UE_SB_injet.at(i_cent_cuts).at(i_jet_bin)->Draw();
					else h_ChPS_UE_SB_injet.at(i_cent_cuts).at(i_jet_bin)->Draw("same");
					gPad->SetLogx();
					gPad->SetLogy(0);

				}
				// ***



				if (i_cent_cuts == 0)
				{
					double pt_lo = jetpT_binning->GetBinLowEdge(i_jet_bin+1);
					double pt_hi = jetpT_binning->GetBinLowEdge(i_jet_bin+2);

					name = Form("%2.0f < p_{T}^{Jet} < %2.0f GeV",pt_lo, pt_hi);
					legend->AddEntry(h_ChPS_tru.at(i_cent_cuts).at(i_dR).at(i_jet_bin),name.c_str(),"lp");
				}

				style++;
			}

			// *****
			c1->cd(i_cent_cuts+1);
			ltx->SetTextAlign(32);
			ltx->DrawLatexNDC(1,0.98,num_to_cent(31, i_cent_cuts).c_str());
			ltx->SetTextAlign(32);
			name = "Raw (No Subtr)";
			if (do_UE_subtr) name = "Raw+Subtr";
			if (do_UE_subtr && do_unfolding) name = "Raw+Subtr+Unf";
			if (do_UE_subtr && do_unfolding && do_BbB) name = "Raw+Subtr+Unf+BbB";
			ltx->DrawLatexNDC(0.93,0.89,name.c_str());
			ltx->DrawLatexNDC(0.93,0.48,"Truth");
			ltx->SetTextAlign(12);
			ltx->DrawLatexNDC(0.19,0.98,Form("%4.2f < dR < %4.2f", dR_lo, dR_hi));
			// *****

			// *****
			c8->cd(i_cent_cuts+1);
			ltx->SetTextAlign(32);
			ltx->DrawLatexNDC(1,0.98,num_to_cent(31, i_cent_cuts).c_str());
			ltx->SetTextAlign(32);

			name = "Raw (No Subtr)";
			if (do_UE_subtr) name = "Raw+Subtr";
			if (do_UE_subtr && do_unfolding) name = "Raw+Subtr+Unf";
			if (do_UE_subtr && do_unfolding && do_BbB) name = "Raw+Subtr+Unf+BbB";
			ltx->DrawLatexNDC(0.93,0.89,name.c_str());
			ltx->SetTextAlign(12);
			ltx->DrawLatexNDC(0.19,0.98,Form("%4.2f < dR < %4.2f", dR_lo, dR_hi));
			// *****


			// *****
			c2->cd(i_cent_cuts+1);
			ltx->SetTextAlign(32);
			ltx->DrawLatexNDC(0.92,0.98,num_to_cent(31, i_cent_cuts).c_str());

			ltx->SetTextAlign(12);
			ltx->DrawLatexNDC(0.19,0.98,Form("%4.2f < dR < %4.2f", dR_lo, dR_hi));
			// *****

			// *****
			c7->cd(i_cent_cuts+1);
			ltx->SetTextAlign(32);
			ltx->DrawLatexNDC(0.92,0.98,num_to_cent(31, i_cent_cuts).c_str());

			line->SetLineStyle(3);
			line->DrawLine(1,1.02, 150, 1.02);
			line->DrawLine(1,0.98, 150, 0.98);

			line->SetLineStyle(2);
			line->DrawLine(1,1, 150, 1);

			ltx->SetTextAlign(12);
			ltx->DrawLatexNDC(0.19,0.98,Form("%4.2f < dR < %4.2f", dR_lo, dR_hi));
			// *****


			// *****
			c3->cd(i_cent_cuts+1);
			ltx->SetTextAlign(32);
			ltx->DrawLatexNDC(0.92,0.98,num_to_cent(31, i_cent_cuts).c_str());

			ltx->SetTextAlign(12);
			ltx->DrawLatexNDC(0.19,0.92,Form("%4.2f < dR < %4.2f", dR_lo, dR_hi));
			// *****


			// *****
			c5->cd(i_cent_cuts+1);
			ltx->SetTextAlign(32);
			ltx->DrawLatexNDC(0.92,0.98,num_to_cent(31, i_cent_cuts).c_str());

			ltx->SetTextAlign(12);
			ltx->DrawLatexNDC(0.19,0.92,Form("%4.2f < dR < %4.2f", dR_lo, dR_hi));
			// *****


			// *****
			if (doing_injet)
			{
				c4->cd(i_cent_cuts+1);
				ltx->SetTextAlign(32);
				ltx->DrawLatexNDC(0.92,0.98,num_to_cent(31, i_cent_cuts).c_str());
				ltx->SetTextAlign(12);
				ltx->DrawLatexNDC(0.19,0.98,Form("injet"));

				c6->cd(i_cent_cuts+1);
				ltx->SetTextAlign(32);
				ltx->DrawLatexNDC(0.92,0.98,num_to_cent(31, i_cent_cuts).c_str());
				ltx->SetTextAlign(12);
				ltx->DrawLatexNDC(0.19,0.98,Form("injet"));
			}
			// *****


		} //end cent loop


		c1->cd(1);
		legend->SetX1NDC(0.20);
		legend->SetY1NDC(0.20);
		legend->SetX2NDC(0.48);
		legend->SetY2NDC(0.40);
		legend->Draw();
		if (i_dR == 0) name = "(";
		else if (i_dR == N_dR - 1) name = ")";
		else name = "";
		c1->Print(Form("%s_spectra_%s_%s_%s.pdf%s",ChPS_raw_type.c_str(), dataset.c_str(), dataset_type.c_str(), evol.c_str(), name.c_str()),Form("Title: dR%i",i_dR));
		c1->Clear();


		c8->cd(1);
		legend->SetX1NDC(0.20);
		legend->SetY1NDC(0.20);
		legend->SetX2NDC(0.48);
		legend->SetY2NDC(0.40);
		legend->Draw();
		if (i_dR == 0) name = "(";
		else if (i_dR == N_dR - 1) name = ")";
		else name = "";
		c8->Print(Form("final_ChPS_%s_%s_%s.pdf%s",dataset.c_str(), dataset_type.c_str(), evol.c_str(), name.c_str()),Form("Title: dR%i",i_dR));
		c8->Clear();


		c2->cd(1)->cd(1);
		legend->SetX1NDC(0.50);
		legend->SetX2NDC(0.80);
		legend->SetY1NDC(0.50);
		legend->SetY2NDC(0.85);
		legend->Draw();
		c2->cd(1);
		ltx->SetTextAlign(12);
		if (i_dR == 0) name = "(";
		else if (i_dR == N_dR - 1) name = ")";
		else name = "";
		c2->Print(Form("%s_closure_%s_%s_%s.pdf%s",ChPS_raw_type.c_str(),dataset.c_str(), dataset_type.c_str(), evol.c_str(), name.c_str()),Form("Title: dR%i",i_dR));
		c2->Print("tmp.root");
		c2->Clear();

		c7->cd(1);
		legend->SetX1NDC(0.20);
		legend->SetY1NDC(0.20);
		legend->SetX2NDC(0.40);
		legend->SetY2NDC(0.40);
		legend->Draw();
		c7->cd(1);
		ltx->SetTextAlign(12);
		if (i_dR == 0) name = "(";
		else if (i_dR == N_dR - 1) name = ")";
		else name = "";
		c7->Print(Form("%s_ratio_%s_%s_%s.pdf%s",ChPS_raw_type.c_str(), dataset.c_str(), dataset_type.c_str(), evol.c_str(),name.c_str()),Form("Title: dR%i",i_dR));
		c7->Clear();


		c3->cd(1);
		legend->SetX1NDC(0.48);
		legend->SetX2NDC(0.80);
		legend->SetY1NDC(0.60);
		legend->SetY2NDC(0.80);
		legend->Draw();
		c3->cd(1);
		ltx->SetTextAlign(12);
		if (i_dR == 0) name = "(";
		else if (i_dR == N_dR - 1) name = ")";
		else name = "";
		c3->Print(Form("UE_%s_%s_%s.pdf%s", dataset.c_str(), dataset_type.c_str(), evol.c_str(),name.c_str()),Form("Title: dR%i",i_dR));
		c3->Clear();


		c5->cd(1);
		legend->SetX1NDC(0.48);
		legend->SetX2NDC(0.80);
		legend->SetY1NDC(0.60);
		legend->SetY2NDC(0.80);
		legend->Draw();
		c5->cd(1);
		ltx->SetTextAlign(12);
		if (i_dR == 0) name = "(";
		else if (i_dR == N_dR - 1) name = ")";
		else name = "";
		c5->Print(Form("UE_SB_%s_%s_%s.pdf%s", dataset.c_str(), dataset_type.c_str(), evol.c_str(),name.c_str()),Form("Title: dR%i",i_dR));
		c5->Clear();

		if (doing_injet)
		{
			c4->cd(1);
			legend->SetX1NDC(0.48);
			legend->SetX2NDC(0.80);
			legend->SetY1NDC(0.60);
			legend->SetY2NDC(0.80);
			legend->Draw();
			c4->cd(1);
			ltx->SetTextAlign(12);

			name = Form("UE_injet.pdf");
			c4->Print(name.c_str(),Form("Title: inJet"));
			c4->Clear();

			c6->cd(1);
			legend->SetX1NDC(0.48);
			legend->SetX2NDC(0.80);
			legend->SetY1NDC(0.60);
			legend->SetY2NDC(0.80);
			legend->Draw();
			c6->cd(1);
			ltx->SetTextAlign(12);

			name = Form("UE_SB_injet.pdf");
			c6->Print(name.c_str(),Form("Title: inJet"));
			c6->Clear();
		}



		doing_injet = false;
		legend->Clear();

	} //end DR loop



}
