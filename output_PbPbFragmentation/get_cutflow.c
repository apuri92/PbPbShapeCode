//#include "combine_eff_dev.c"
#include "extras/global_variables.h"

void get_cutflow(bool isPerf, bool isMC, string cut)
{

	cout << "*********** Begin Cutflow... ***********" << endl;
	string filename;
	string name;
	if (isMC && isPerf) filename = Form("mc_perf");
	else if (isMC && !isPerf) filename = Form("mc");
	else if (!isMC) filename = Form("data");
	else {cout << "ERROR IN COMBINING HISTOGRAMS" << endl; return;}

	TFile *input = new TFile(Form("input/Perf_MC_JZ_comb_out_histo_PbPb_5p02_r001.root",filename.c_str(),cut.c_str()));
//	TFile *input = new TFile(Form("jz3_perf_pptight.root"));

	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);

	TCanvas *canvas1 = new TCanvas("C1", "C1",0.,0.,800,600);
	TCanvas *canvas2 = new TCanvas("C2", "C2",0.,0.,800,600);

	TLine *line = new TLine();
	TLatex *ltx = new TLatex();
	ltx->SetNDC();
	ltx->SetTextFont(43);
	ltx->SetTextSize(15);
	ltx->SetTextAlign(11);
	TLegend *legend = new TLegend(0.12,0.11,0.60,0.40,NULL,"brNDC");
	legend->SetBorderSize(0);
	legend->SetNColumns(2);
	legend->SetTextFont(43);
	legend->SetTextSize(20);

	TH3 *h_tmp_3d;
	TH3 *h_tmp_1d;
	TH3 *h_truth_3d;
	TH1 *h_truth;

	TH1 *h_cutflow[14];
	TH1 *h_frac[14];
	TH1 *h_wrt_truth[14];


	vector<string> CutsOn;
	CutsOn.push_back("Reco"); //Must be first
	CutsOn.push_back("AllCuts"); //Must be last
	CutsOn.push_back("d0");


//	CutsOn.push_back("Z0");
//	CutsOn.push_back("z0sintheta");
//	CutsOn.push_back("TrtHits");
//	CutsOn.push_back("Truth");

	if (cut == "HITight")
	{
		CutsOn.push_back("InnermostLayersHits");
		CutsOn.push_back("PixelHits");
		CutsOn.push_back("SctHits");
		CutsOn.push_back("Z0SinTheta");
//		CutsOn.push_back("z0sintheta");
		CutsOn.push_back("FitQuality");
	}
	if (cut == "ppTight")
	{
		CutsOn.push_back("SiHits");
		CutsOn.push_back("InnermostLayersHits");
//		CutsOn.push_back("Z0SinTheta");
		CutsOn.push_back("z0sintheta");
	}
	if (cut == "ppTight_tight")
	{
		CutsOn.push_back("SiHits");
		CutsOn.push_back("InnermostLayersHits");
//		CutsOn.push_back("Z0SinTheta");
		CutsOn.push_back("z0sintheta");
		CutsOn.push_back("d0Sign");
		CutsOn.push_back("z0sinthetaSign");
	}
	
	h_truth_3d = (TH3*)input->Get("h_eff_Injet_cent6");
	h_truth = (TH1*)h_truth_3d->Project3D("y");

	TH1* h_cut_truth;
	h_cut_truth = (TH1*)h_truth->Clone("cut_truth");
	

	for (int i = 0; i < CutsOn.size(); i++)
	{
		h_tmp_3d = (TH3*)input->Get(Form("h_cut_flow_%s_cent6",CutsOn.at(i).c_str()));
		h_cutflow[i] = (TH1*)h_tmp_3d->Project3D("x");
		h_cutflow[i]->SetName(Form("h_cutflow_%s",CutsOn.at(i).c_str()));

		hstyle(h_cutflow[i],i);
		h_cutflow[i]->GetXaxis()->SetRangeUser(1,500);

		h_frac[i] = (TH1*)h_cutflow[i]->Clone(CutsOn.at(i).c_str());
		h_frac[i]->Divide(h_cutflow[0]);

		if (CutsOn.at(i) == "Reco") continue;
		legend->AddEntry(h_cutflow[i],CutsOn.at(i).c_str(),"lp");

		canvas1->cd();
		h_frac[i]->GetYaxis()->SetRangeUser(0.5,1.3);
		gPad->SetTicks(1,1);
		h_frac[i]->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]" );
		h_frac[i]->GetYaxis()->SetTitle("Passed/Reco" );
		gPad->SetLogx();
		if (i==0) h_frac[i]->Draw();
		else h_frac[i]->Draw("same");

	}

	canvas1->cd();
	PlotLabels_AtlasSim_q2_5(0.11, 0.84, 25, 11, false);
	legend->Draw();
	name = Form("%s", cut.c_str());
	ltx->SetTextSize(23);
	ltx->SetTextAlign(31);
	ltx->DrawLatex(0.88,0.84,cut.c_str());

//	h_cut_truth->Divide(h_cutflow[0]);
//	h_cut_truth->Draw("same hist");

	line->SetLineColor(kRed);
	line->SetLineStyle(3);
	line->DrawLine(4,1,500,1);

	legend->Draw();
	canvas1->Print(Form("cutflow_%s.pdf",cut.c_str()));


	canvas1->Clear();
	legend->Clear();

	
	for (int i = 0; i < CutsOn.size(); i++)
	{
		h_wrt_truth[i] = (TH1*)h_cutflow[i]->Clone(CutsOn.at(i).c_str());
		h_wrt_truth[i]->Divide(h_truth);

		legend->AddEntry(h_wrt_truth[i],CutsOn.at(i).c_str(),"lp");

		canvas1->cd();
		h_wrt_truth[i]->GetYaxis()->SetRangeUser(0,5);

		h_wrt_truth[i]->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]" );
		h_wrt_truth[i]->GetYaxis()->SetTitle("Passed/Truth" );
		gPad->SetLogx();
		if (i==0) h_wrt_truth[i]->Draw();
		else h_wrt_truth[i]->Draw("same");
		
	}

	line->SetLineColor(kRed);
	line->SetLineStyle(3);
	line->DrawLine(4,1,500,1);

	legend->Draw();
	canvas1->Print(Form("cutflow_truth_%s.pdf",cut.c_str()));

	canvas1->Clear();
	legend->Clear();



	TH1* h_d0sign_wrt_truth;
	TH1* h_z0sign_wrt_truth;


	h_tmp_3d = (TH3*)input->Get("h_d0sign_wrt_truth");
	h_d0sign_wrt_truth = (TH1*)h_tmp_3d->Project3D("x");
	h_tmp_3d = (TH3*)input->Get("h_z0sign_wrt_truth");
	h_z0sign_wrt_truth = (TH1*)h_tmp_3d->Project3D("x");

	h_d0sign_wrt_truth->Divide(h_truth);
	h_z0sign_wrt_truth->Divide(h_truth);

	hstyle(h_d0sign_wrt_truth,0);
	hstyle(h_z0sign_wrt_truth,2);

	legend->AddEntry(h_d0sign_wrt_truth,"d0 Sign. wrt truth","lp");
	legend->AddEntry(h_z0sign_wrt_truth,"d0 Sign. wrt truth","lp");

	canvas1->cd();
	h_d0sign_wrt_truth->GetYaxis()->SetRangeUser(0,1.2);

	h_d0sign_wrt_truth->GetXaxis()->SetTitle("p_{T}^{Trk} [GeV]" );
	h_d0sign_wrt_truth->GetYaxis()->SetTitle("Passed (matched)/Truth" );
	gPad->SetLogx();
	h_d0sign_wrt_truth->Draw();
	h_z0sign_wrt_truth->Draw("same");

	legend->Draw();
	canvas1->Print(Form("sign_truth_%s.pdf",cut.c_str()));


}

