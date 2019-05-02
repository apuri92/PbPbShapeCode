#include "drawingClass.h"
#include "/Users/Akshat/Dropbox/.RootUtils/Label.C"

using namespace std;
void drawingClass::loadObjects()
{
	canvas = new TCanvas("canvas","canvas",canvas_x,canvas_y);

	legend = new TLegend(legend_x1, legend_x2, legend_y1, legend_y2, "","brNDC");
	legend->SetTextFont(43);
	legend->SetBorderSize(0);
	legend->SetTextSize(legend_text_size);

	ltx = new TLatex();
	ltx->SetTextFont(43);
	ltx->SetTextSize(12);

	line = new TLine();
}

void drawingClass::drawCentPanels(vector<vector<vector<TH1*>>> vvvhist, string name)
{
	drawingClass::setSpecifics(name);

	legend->Clear();

	int jet_itr = 0;
	for (int i_jet = jetStart; i_jet < jetEnd; i_jet++)
	{
		canvas->cd();
		canvas->Clear();
		canvas->Divide(3,2);

		int cent_itr = 0;
		for (int i_cent = 0; i_cent < NCent; i_cent++)
		{
			int trk_itr = 0;
			for (int i_trk = trkStart; i_trk < trkEnd; i_trk++)
			{
				SetHStyle_smallify(vvvhist[i_trk][i_jet][i_cent], trk_itr, 1);
				if (jet_itr == 0 && cent_itr == 0) legend->AddEntry(vvvhist[i_trk][i_jet][i_cent],trk_label(i_trk+1).c_str(),"lp");

				vvvhist[i_trk][i_jet][i_cent]->GetXaxis()->SetRangeUser(0, 0.8);
				vvvhist[i_trk][i_jet][i_cent]->GetYaxis()->SetRangeUser(y_lo, y_hi);
				vvvhist[i_trk][i_jet][i_cent]->GetYaxis()->SetNdivisions(504);

				canvas->cd(i_cent+1);
				if (trk_itr == 0) vvvhist[i_trk][i_jet][i_cent]->Draw("");
				else vvvhist[i_trk][i_jet][i_cent]->Draw("same");
				gPad->SetLogx(logx);
				gPad->SetLogy(logy);

				trk_itr++;
			}

			canvas->cd(i_cent+1);
			ltx->SetTextAlign(32);
			ltx->SetTextSize(12);
			ltx->DrawLatexNDC(0.93, 0.90, Form("%s", jet_label(i_jet+1).c_str()));
			ltx->DrawLatexNDC(0.93, 0.85, Form("%s", cent_label(i_cent).c_str()));
			line->DrawLine(line_x1, line_y1, line_x2, line_y2);
			legend->Draw();

			cent_itr++;
		}

		string pdf_label = "";
		if (i_jet == jetStart) pdf_label = "(";
		if (i_jet == jetEnd-1) pdf_label = ")";
		canvas->Print(Form("figures/%s_config%i.pdf%s", name.c_str(), config, pdf_label.c_str()), Form("Title:jetpt%i", i_jet));


		jet_itr++;
	}

}



void drawingClass::drawCentPanels(vector<vector<TH1*>> vvhist, string name)
{
	drawingClass::setSpecifics(name);

	legend->Clear();
	canvas->cd();
	canvas->Clear();
	canvas->Divide(3,2);

	int cent_itr = 0;
	for (int i_cent = 0; i_cent < NCent; i_cent++)
	{
		int jet_itr = 0;
		for (int i_jet = jetStart; i_jet < jetEnd; i_jet++)
		{
			SetHStyle_smallify(vvhist[i_jet][i_cent], jet_itr, 1);
			if (cent_itr == 0) legend->AddEntry(vvhist[i_jet][i_cent],jet_label(i_jet+1).c_str(),"lp");

			vvhist[i_jet][i_cent]->GetXaxis()->SetRangeUser(0, 0.8);
			vvhist[i_jet][i_cent]->GetYaxis()->SetRangeUser(y_lo, y_hi);
			vvhist[i_jet][i_cent]->GetYaxis()->SetNdivisions(504);

			canvas->cd(i_cent+1);
			if (jet_itr == 0) vvhist[i_jet][i_cent]->Draw("");
			else vvhist[i_jet][i_cent]->Draw("same");
			gPad->SetLogx(logx);
			gPad->SetLogy(logy);

			jet_itr++;
		}

		canvas->cd(i_cent+1);
		ltx->SetTextAlign(32);
		ltx->SetTextSize(12);
		ltx->DrawLatexNDC(0.93, 0.85, Form("%s", cent_label(i_cent).c_str()));
		line->DrawLine(line_x1, line_y1, line_x2, line_y2);
		legend->Draw();

		cent_itr++;
	}
	canvas->Print(Form("figures/%s_config%i.pdf", name.c_str(), config));

}


void drawingClass::loadAxes()
{
	TFile *file = new TFile(Form("../output_pdf_nominal/root/final_ChPS_MC_PbPb.root"));
	dR_binning = (TAxis*)file->Get("dR_binning");
	jetpT_binning = (TAxis*)file->Get("jetpT_binning");
	trkpT_binning = (TAxis*)file->Get("trkpT_binning");
	file->Close();
	delete file;
}


std::string drawingClass::trk_label(int bin)
{
	return Form("%1.1f < #it{p}_{T}^{trk} < %1.1f",
				trkpT_binning->GetBinLowEdge(bin),
				trkpT_binning->GetBinUpEdge(bin));
}

std::string drawingClass::jet_label(int bin)
{
	return Form("%1.0f < #it{p}_{T}^{jet} < %1.0f",
				jetpT_binning->GetBinLowEdge(bin),
				jetpT_binning->GetBinUpEdge(bin));
}

std::string drawingClass::dR_label(int bin)
{
	return Form("%1.2f < #it{r} < %1.2f",
				dR_binning->GetBinLowEdge(bin),
				dR_binning->GetBinUpEdge(bin));
}

std::string drawingClass::cent_label(int bin)
{
	if (bin == 0)		return " 0 - 10%";
	else if (bin == 1)	return "10 - 20%";
	else if (bin == 2) 	return "20 - 30%";
	else if (bin == 3) 	return "30 - 40%";
	else if (bin == 4) 	return "40s - 60%";
	else if (bin == 5) 	return "60 - 80%";
	else return "Incorrect bin";
}

void drawingClass::setSpecifics(string hist_name)
{
	if (hist_name.compare("rdpt") == 0)
	{
		y_lo = 0, y_hi = 5;
		line_y1 = 1, line_y2 = 1;

	}
	if (hist_name.compare("deltadpt") == 0)
	{
		y_lo = -1, y_hi = 9;
		line_y1 = 0, line_y2 = 0;
	}

	if (hist_name.compare("lowpt_integ_rdpt") == 0)
	{
		y_lo = 0, y_hi = 3;
		line_y1 = 1, line_y2 = 1;
	}

	if (hist_name.compare("jetshape_rdpt") == 0)
	{
		y_lo = 0, y_hi = 3;
		line_y1 = 1, line_y2 = 1;
	}

	if (hist_name.compare("jetshape_deltadpt") == 0)
	{
		y_lo = -1, y_hi = 5;
		line_y1 = 0, line_y2 = 0;
	}

	if (hist_name.compare("lowpt_integ_deltadpt") == 0)
	{
		y_lo = -1, y_hi = 5;
		line_y1 = 0, line_y2 = 0;
	}



}
