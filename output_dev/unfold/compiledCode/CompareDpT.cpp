#include "CompareDpT.h"

using namespace std;

void CompareDpT::getDpT()
{
	cout << "Getting pp and PbPb DpT" << endl;

	vector<vector<vector<TH1*>>>	vh_PbPb (NTrk, vector<vector<TH1*>> (NJet, vector<TH1*> (NCent)));
	vector<vector<TH1*>>			vh_pp (NTrk, vector<TH1*> (NJet));

	for (int i_cent = 0; i_cent < NCent; i_cent++)
	{
		for (int i_trk = 0; i_trk < NTrk; i_trk++)
		{
			for (int i_jet = 0; i_jet < NJet; i_jet++)
			{
				if (i_cent == 0)
				{
					vh_pp[i_trk][i_jet] = (TH1*)inputFile_pp->Get(Form(histName.c_str(), i_trk, 6, i_jet));
					vh_pp[i_trk][i_jet]->SetName(Form("%s_pp",vh_pp[i_trk][i_jet]->GetName() ));
				}

				vh_PbPb[i_trk][i_jet][i_cent] = (TH1*)inputFile_PbPb->Get(Form(histName.c_str(), i_trk, i_cent, i_jet));
				vh_PbPb[i_trk][i_jet][i_cent]->SetName(Form("%s_PbPb",vh_PbPb[i_trk][i_jet][i_cent]->GetName() ));
			}
		}
	}

	spect_PbPb = vh_PbPb;
	spect_pp = vh_pp;

	CompareDpT::labelvHist(spect_pp, r_title, dptr_title);
	CompareDpT::labelvHist(spect_PbPb, r_title, dptr_title);


}


void CompareDpT::getRDpT()
{
	cout << "Getting RDpT distributions" << endl;

	vector<vector<vector<TH1*>>> vh_result (NTrk, vector<vector<TH1*>> (NJet, vector<TH1*> (NCent)));

	for (int i_cent = 0; i_cent < NCent; i_cent++)
	{
		for (int i_trk = 0; i_trk < NTrk; i_trk++)
		{
			for (int i_jet = 0; i_jet < NJet; i_jet++)
			{
				string name = Form(histName.c_str(), i_trk, i_cent, i_jet);
				vh_result[i_trk][i_jet][i_cent] = (TH1*)spect_PbPb[i_trk][i_jet][i_cent]->Clone(Form("%s_ratio",name.c_str()));
				vh_result[i_trk][i_jet][i_cent]->Divide(spect_pp[i_trk][i_jet]);
			}
		}
	}

	spect_rdpt = vh_result;
	CompareDpT::labelvHist(spect_rdpt, r_title, rdptr_title);

}


void CompareDpT::getDeltaDpT()
{
	cout << "Getting Delta DpT" << endl;

	vector<vector<vector<TH1*>>> vh_result (NTrk, vector<vector<TH1*>> (NJet, vector<TH1*> (NCent)));

	for (int i_cent = 0; i_cent < NCent; i_cent++)
	{
		for (int i_trk = 0; i_trk < NTrk; i_trk++)
		{
			for (int i_jet = 0; i_jet < NJet; i_jet++)
			{
				string name = Form(histName.c_str(), i_trk, i_cent, i_jet);
				vh_result[i_trk][i_jet][i_cent] = (TH1*)spect_PbPb[i_trk][i_jet][i_cent]->Clone(Form("%s_ratio",name.c_str()));
				vh_result[i_trk][i_jet][i_cent]->Add(spect_pp[i_trk][i_jet],-1);
			}
		}
	}

	spect_deltadpt = vh_result;
	CompareDpT::labelvHist(spect_deltadpt, r_title, deltadptr_title);
}

TH1* CompareDpT::cumulateR(TH1* hist)
{

	/*
	 scale back up by area of ring, add up cumulatively
	 and normalize by area of disk
	 */

	TH1* hist_integrated = (TH1*)hist->Clone(Form("%s_cumulateR",hist->GetName()));

	double bin_cont = 0, bin_err = 0;
	for (int i_dR = 1; i_dR <= NdR; i_dR++)
	{
		double r1 = dR_binning->GetBinLowEdge(i_dR);
		double r2 = dR_binning->GetBinUpEdge(i_dR);
		double area = getArea(r1, r2);

		bin_cont = hist_integrated->GetBinContent(i_dR) * area;
		bin_err = hist_integrated->GetBinError(i_dR) * area;
		hist_integrated->SetBinContent(i_dR, bin_cont);
		hist_integrated->SetBinError(i_dR, bin_err);
	}


	for (int i_dR = 1; i_dR <= NdR; i_dR++)
	{
		bin_cont = hist_integrated->GetBinContent(i_dR-1) + hist_integrated->GetBinContent(i_dR);
		hist_integrated->SetBinContent(i_dR, bin_cont);

		bin_err = pow(hist_integrated->GetBinError(i_dR-1),2) + pow(hist_integrated->GetBinError(i_dR),2);
		hist_integrated->SetBinError(i_dR, sqrt(bin_err));
	}


	for (int i_dR = 1; i_dR <= NdR; i_dR++)
	{
		double r = dR_binning->GetBinUpEdge(i_dR);
		double area = getArea(r);

		bin_cont = hist_integrated->GetBinContent(i_dR) / area;
		bin_err = hist_integrated->GetBinError(i_dR) / area;
		hist_integrated->SetBinContent(i_dR, bin_cont);
		hist_integrated->SetBinError(i_dR, bin_err);
	}

	return hist_integrated;
}


void CompareDpT::getJetShape()
{
	cout << "Getting Getting jet shape distributions" << endl;

	vector<vector<TH1*>>	vh_lowpt_integ_rdpt (NJet, vector<TH1*> (NCent));
	vector<vector<TH1*>>	vh_lowpt_integ_deltadpt (NJet, vector<TH1*> (NCent));
	vector<vector<TH1*>>	vh_jetshape_rdpt (NJet, vector<TH1*> (NCent));
	vector<vector<TH1*>>	vh_jetshape_deltadpt (NJet, vector<TH1*> (NCent));

	vector<vector<TH1*>>	vh_lowpt_integ_PbPb  (NJet, vector<TH1*> (NCent));
	vector<TH1*>			vh_lowpt_integ_pp (NJet);
	vector<vector<TH1*>>	vh_jetshape_PbPb  (NJet, vector<TH1*> (NCent));
	vector<TH1*>			vh_jetshape_pp (NJet);




	TH1* h_PbPb;
	TH1* h_pp;
	TH1* h_lowpT_int_rdpt;
	TH1* h_lowpT_int_deltadpt;
	TH1* h_lowpT_jetshape_rdpt;
	TH1* h_lowpT_jetshape_deltadpt;

	int pt_lo_index = trkpT_binning->FindBin(1.01) - 1;
	int pt_hi_index = trkpT_binning->FindBin(3.95) - 1;

	for (int i_cent = 0; i_cent < NCent; i_cent++)
	{
		for (int i_jet = 0; i_jet < NJet; i_jet++)
		{
			double bin_width = 0, sum_bin_width = 0;
			for (int i_trk = pt_lo_index; i_trk <= pt_hi_index; i_trk++)
			{
				bin_width = trkpT_binning->GetBinWidth(i_trk+1);
				sum_bin_width +=bin_width;

				if (i_trk == pt_lo_index)
				{
					h_pp = (TH1*)spect_pp[i_trk][i_jet]->Clone(Form("pp_integ_jet%i_cent%i",i_jet, i_cent));
					h_pp->Scale(bin_width);

					h_PbPb = (TH1*)spect_PbPb[i_trk][i_jet][i_cent]->Clone(Form("PbPb_integ_jet%i_cent%i",i_jet, i_cent));
					h_PbPb->Scale(bin_width);

				}
				else
				{
					h_pp->Add(spect_pp[i_trk][i_jet], bin_width);
					h_PbPb->Add(spect_PbPb[i_trk][i_jet][i_cent], bin_width);
				}
			}

			h_pp->Scale(1./sum_bin_width);
			h_PbPb->Scale(1./sum_bin_width);

			vh_lowpt_integ_PbPb[i_jet][i_cent] = (TH1*)h_PbPb->Clone(Form("PbPb_lowpT_integrated_ChPS_jet%i_cent%i",i_jet, i_cent));
			vh_lowpt_integ_pp[i_jet] = (TH1*)h_pp->Clone(Form("pp_lowpT_integrated_ChPS_jet%i_cent%i",i_jet, 6));

			h_lowpT_int_rdpt = (TH1*)h_PbPb->Clone(Form("lowpT_integrated_rdpt_jet%i_cent%i",i_jet, i_cent));
			h_lowpT_int_rdpt->Divide(h_pp);

			h_lowpT_int_deltadpt = (TH1*)h_PbPb->Clone(Form("lowpT_integrated_deltadpt_jet%i_cent%i",i_jet, i_cent));
			h_lowpT_int_deltadpt->Add(h_pp,-1);

			h_pp = cumulateR(h_pp);
			h_PbPb = cumulateR(h_PbPb);

			vh_jetshape_PbPb[i_jet][i_cent] = (TH1*)h_PbPb->Clone(Form("PbPb_jetshape_ChPS_jet%i_cent%i",i_jet, i_cent));
			vh_jetshape_pp[i_jet] = (TH1*)h_pp->Clone(Form("pp_jetshape_ChPS_jet%i_cent%i",i_jet, 6));

			h_lowpT_jetshape_rdpt = (TH1*)h_PbPb->Clone(Form("lowpT_jetshape_rdpt_jet%i_cent%i",i_jet, i_cent));
			h_lowpT_jetshape_rdpt->Divide(h_pp);

			h_lowpT_jetshape_deltadpt = (TH1*)h_PbPb->Clone(Form("lowpT_jetshape_deltadpt_jet%i_cent%i",i_jet, i_cent));
			h_lowpT_jetshape_deltadpt->Add(h_pp,-1);


			vh_lowpt_integ_rdpt[i_jet][i_cent] = h_lowpT_int_rdpt;
			vh_lowpt_integ_deltadpt[i_jet][i_cent] = h_lowpT_int_deltadpt;
			vh_jetshape_rdpt[i_jet][i_cent] = h_lowpT_jetshape_rdpt;
			vh_jetshape_deltadpt[i_jet][i_cent] = h_lowpT_jetshape_deltadpt;
		}
	}



	spect_lowpt_integ_PbPb = vh_lowpt_integ_PbPb;
	spect_lowpt_integ_pp = vh_lowpt_integ_pp;
	spect_jetshape_PbPb = vh_jetshape_PbPb;
	spect_jetshape_pp = vh_jetshape_pp;


	spect_lowpt_integ_rdpt = vh_lowpt_integ_rdpt;
	spect_lowpt_integ_deltadpt = vh_lowpt_integ_deltadpt;
	spect_jetshape_rdpt = vh_jetshape_rdpt;
	spect_jetshape_deltadpt = vh_jetshape_deltadpt;



	CompareDpT::labelvHist(spect_lowpt_integ_rdpt, r_title, lowpt_integ_rdpt_title);
	CompareDpT::labelvHist(spect_lowpt_integ_deltadpt, r_title, lowpt_integ_deltadpt_title);
	CompareDpT::labelvHist(spect_jetshape_rdpt, r_title, jetshape_rdpt_title);
	CompareDpT::labelvHist(spect_jetshape_deltadpt, r_title, jetshape_deltadpt_title);

	CompareDpT::labelvHist(spect_lowpt_integ_PbPb, r_title, lowpt_integ_dptr_title);
	CompareDpT::labelvHist(spect_lowpt_integ_pp, r_title, lowpt_integ_dptr_title);
	CompareDpT::labelvHist(spect_jetshape_PbPb, r_title, jetshape_dptr_title);
	CompareDpT::labelvHist(spect_jetshape_pp, r_title, jetshape_dptr_title);

}


void CompareDpT::writeToFile()
{
	cout << "Writing to file: " << output_file->GetName() << endl;
	output_file->cd();
	string name;

	dR_binning->Write("dR_binning");
	jetpT_binning->Write("jetpT_binning");
	trkpT_binning->Write("trkpT_binning");

	for (int i_jet = 0; i_jet < NJet; i_jet++)
	{
		name = Form("h_lowpt_integ_pp_ChPS_final_indR_cent%i_jetpt%i", 6, i_jet);
		spect_lowpt_integ_pp[i_jet]->Write(name.c_str());

		name = Form("h_jetshape_pp_ChPS_final_indR_cent%i_jetpt%i", 6, i_jet);
		spect_jetshape_pp[i_jet]->Write(name.c_str());

		for (int i_cent = 0; i_cent < NCent; i_cent++)
		{
			name = Form("h_lowpt_integ_PbPb_ChPS_final_indR_cent%i_jetpt%i", i_cent, i_jet);
			spect_lowpt_integ_PbPb[i_jet][i_cent]->Write(name.c_str());

			name = Form("h_jetshape_PbPb_ChPS_final_indR_cent%i_jetpt%i", i_cent, i_jet);
			spect_jetshape_PbPb[i_jet][i_cent]->Write(name.c_str());

		}
	}



	for (int i_cent = 0; i_cent < NCent; i_cent++)
	{
		for (int i_jet = 0; i_jet < NJet; i_jet++)
		{
			for (int i_trk = 0; i_trk < NTrk; i_trk++)
			{

				if (i_cent == 0)
				{
					name = Form("h_pp_ChPS_final_indR_trk%i_cent%i_jetpt%i", i_trk, 6, i_jet);
					spect_pp[i_trk][i_jet]->Write(Form("%s",name.c_str()));
				}

				name = Form("h_PbPb_ChPS_final_indR_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				spect_PbPb[i_trk][i_jet][i_cent]->Write(Form("%s",name.c_str()));

				name = Form("h_RDpT_final_indR_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				spect_rdpt[i_trk][i_jet][i_cent]->Write(Form("%s",name.c_str()));

				name = Form("h_DeltaDpT_final_indR_trk%i_cent%i_jetpt%i", i_trk, i_cent, i_jet);
				spect_deltadpt[i_trk][i_jet][i_cent]->Write(Form("%s",name.c_str()));
			}

			name = Form("h_lowpt_integ_RDpT_final_indR_cent%i_jetpt%i",i_cent, i_jet);
			spect_lowpt_integ_rdpt[i_jet][i_cent]->Write(Form("%s",name.c_str()));

			name = Form("h_lowpt_integ_DeltaDpT_final_indR_cent%i_jetpt%i",i_cent, i_jet);
			spect_lowpt_integ_deltadpt[i_jet][i_cent]->Write(Form("%s",name.c_str()));

			name = Form("h_jetshape_RDpT_final_indR_cent%i_jetpt%i",i_cent, i_jet);
			spect_jetshape_rdpt[i_jet][i_cent]->Write(Form("%s",name.c_str()));

			name = Form("h_jetshape_DeltaDpT_final_indR_cent%i_jetpt%i",i_cent, i_jet);
			spect_jetshape_deltadpt[i_jet][i_cent]->Write(Form("%s",name.c_str()));


		}
	}


}

void CompareDpT::loadAxes()
{
	TFile *file = new TFile(Form("../output_pdf_nominal/root/final_ChPS_MC_PbPb.root"));

	dR_binning = (TAxis*)file->Get("dR_binning");
	jetpT_binning = (TAxis*)file->Get("jetpT_binning");
	trkpT_binning = (TAxis*)file->Get("trkpT_binning");

	NTrk = trkpT_binning->GetNbins();
	NJet = jetpT_binning->GetNbins();
	NdR = dR_binning->GetNbins();

	file->Close();
	delete (file);
}

void CompareDpT::printInfo()
{
	cout << "Data Type: " << dataType << endl;
	cout << "Dimensions: " << NTrk << " " << NCent << " " << NJet << endl;
}


void CompareDpT::labelvHist(vector<vector<vector<TH1*>>> vvv_hist, string x_label, string y_label)
{
	for (int i = 0; i < vvv_hist.size(); i++) CompareDpT::labelvHist(vvv_hist[i], x_label, y_label);
}


void CompareDpT::labelvHist(vector<vector<TH1*>> vv_hist, string x_label, string y_label)
{
	for (int i = 0; i < vv_hist.size(); i++) {CompareDpT::labelvHist(vv_hist[i], x_label, y_label);}
}

void CompareDpT::labelvHist(vector<TH1*> v_hist, string x_label, string y_label)
{
	for (int i = 0; i < v_hist.size(); i++) {CompareDpT::labelvHist(v_hist[i], x_label, y_label);}
}

void CompareDpT::labelvHist(TH1* hist, string x_label, string y_label)
{
	hist->GetXaxis()->SetTitle(x_label.c_str());
	hist->GetYaxis()->SetTitle(y_label.c_str());
	hist->GetYaxis()->SetNdivisions(504);
}

void CompareDpT::loadFiles(int config)
{
	string file_path = "";
	if (config == 0) file_path = Form("nominal");
	else file_path = Form("sys%i", config);

	inputFile_PbPb = new TFile(Form("../output_pdf_%s/root/final_ChPS_%s_PbPb.root", file_path.c_str(), dataType.c_str()));
	inputFile_pp = new TFile(Form("../output_pdf_%s/root/final_ChPS_%s_pp.root", file_path.c_str(), dataType.c_str()));
	output_file = new TFile(Form("root/output_%s_%s.root",file_path.c_str(), dataType.c_str()), "recreate");
}
