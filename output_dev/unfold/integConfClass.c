#include "integConfClass.h"


void integConfClass::initHist()
{

	dR_binning = (TAxis*)f_nom->Get("dR_binning");
	jetpT_binning = (TAxis*)f_nom->Get("jetpT_binning");
	trkpT_binning = (TAxis*)f_nom->Get("trkpT_binning");

	N_dR = dR_binning->GetNbins();
	N_jetpt = jetpT_binning->GetNbins();
	N_trkpt = trkpT_binning->GetNbins();

	vector<vector<TH1*>> tmp_h_nom (6, vector<TH1*> (11));
	vector<vector<TH1*>> tmp_h_sys_p (6, vector<TH1*> (11));
	vector<vector<TH1*>> tmp_h_sys_n (6, vector<TH1*> (11));
	vector<vector<TGraphAsymmErrors*>> tmp_g_sys (6, vector<TGraphAsymmErrors*> (11));
	vector<vector<TGraphAsymmErrors*>> tmp_g_stat (6, vector<TGraphAsymmErrors*> (11));


	for (int i_jet = 7; i_jet < 11; i_jet++)
	{
		for (int i_cent = 0; i_cent < 6; i_cent++)
		{

			name = Form("h_%s_%s_final_indR_cent%i_jetpt%i", integType.c_str(), mode.c_str(), i_cent, i_jet);
			tmp_h_nom[i_cent][i_jet] = (TH1*)f_nom->Get(name.c_str());

			name = Form("h_%s_sys_cent%i_jetpt%i_total_p",mode.c_str(), i_cent, i_jet);
			tmp_h_sys_p[i_cent][i_jet] = (TH1*)f_sys->Get(name.c_str());

			name = Form("h_%s_sys_cent%i_jetpt%i_total_n",mode.c_str(), i_cent, i_jet);
			tmp_h_sys_n[i_cent][i_jet] = (TH1*)f_sys->Get(name.c_str());

			tmp_g_sys[i_cent][i_jet] = new TGraphAsymmErrors(tmp_h_nom[i_cent][i_jet]);
			tmp_g_sys[i_cent][i_jet]->GetYaxis()->SetTitle(axis_label_y.c_str());
			tmp_g_sys[i_cent][i_jet]->GetXaxis()->SetTitle(axis_label_x.c_str());

			tmp_g_stat[i_cent][i_jet] = new TGraphAsymmErrors(tmp_h_nom[i_cent][i_jet]);
			tmp_g_stat[i_cent][i_jet]->GetYaxis()->SetTitle(axis_label_y.c_str());
			tmp_g_stat[i_cent][i_jet]->GetXaxis()->SetTitle(axis_label_x.c_str());
		}
	}

	h_nom = tmp_h_nom;
	h_sys_p = tmp_h_sys_p;
	h_sys_n = tmp_h_sys_n;
	g_sys = tmp_g_sys;
	g_stat = tmp_g_stat;

}



void integConfClass::makeGraph()
{
	double nom, sys_hi, sys_lo, r_position, r_width, stat_hi, stat_lo;

	for (int i_jet = 7; i_jet < 11; i_jet++)
	{
		for (int i_cent = 0; i_cent < 6; i_cent++)
		{
			for (int i_dR = 0 ; i_dR < N_dR; i_dR++)
			{
				r_position = h_nom[i_cent][i_jet]->GetBinCenter(i_dR+1) + (i_cent*0.001 + i_jet*0.001);
				r_width = 0.020; //h_nom[i_cent][i_jet]->GetBinWidth(i_dR+1)/4;

				nom = h_nom[i_cent][i_jet]->GetBinContent(i_dR+1);
				sys_hi = nom * (h_sys_p[i_cent][i_jet]->GetBinContent(i_dR+1));
				sys_lo = fabs(nom * (h_sys_n[i_cent][i_jet]->GetBinContent(i_dR+1)));
				stat_hi = h_nom[i_cent][i_jet]->GetBinError(i_dR+1);
				stat_lo = h_nom[i_cent][i_jet]->GetBinError(i_dR+1);

				g_stat[i_cent][i_jet]->SetPoint(i_dR, r_position, nom );
				g_stat[i_cent][i_jet]->SetPointError(i_dR, 0, 0, stat_lo, stat_hi);

				g_sys[i_cent][i_jet]->SetPoint(i_dR, r_position, nom );
				g_sys[i_cent][i_jet]->SetPointError(i_dR, r_width/2, r_width/2, sys_lo, sys_hi);
			}
		}
	}
}

void integConfClass::drawAll()
{
	cout << "Drawing..." << endl;

	canvas = new TCanvas("canvas","canvas",700,800);
	legend = new TLegend(legend_x1, legend_y1, legend_x2, legend_y2, "","brNDC");
	legend->SetTextFont(43);
	legend->SetBorderSize(0);
	legend->SetTextSize(10);
	legend->SetNColumns(legend_cols);
	ltx = new TLatex();
	ltx->SetTextFont(43);
	ltx->SetTextSize(12);
	line = new TLine();
	line->SetLineStyle(3);

	string trk_label = Form("%1.1f < #it{p}_{T} < %1.1f GeV", 1.0, 4.0);

	bool first_pass_cent = true;

	canvas->Clear();
	canvas->Divide(2,3);

	for (int i_cent = 0; i_cent < 6; i_cent++)
	{

		string centrality = num_to_cent(31,i_cent);

		canvas->cd(i_cent+1);

		int jet_itr = 0;
		for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
		{
			string jet_label = Form("%1.0f < #it{p}_{T}^{jet} < %1.0f GeV", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));

			SetHStyle_smallify(h_nom[i_cent][i_jet], jet_itr, 1);
			SetHStyle_graph_smallify(g_sys[i_cent][i_jet], jet_itr, 1);
			SetHStyle_graph_smallify(g_stat[i_cent][i_jet], jet_itr, 1);

			if (first_pass_cent) legend->AddEntry(g_sys[i_cent][i_jet],jet_label.c_str(),"p");

			g_sys[i_cent][i_jet]->GetXaxis()->SetRangeUser(0, 0.8);
			g_sys[i_cent][i_jet]->GetYaxis()->SetRangeUser(y_range_lo, y_range_hi);
			g_sys[i_cent][i_jet]->GetYaxis()->SetNdivisions(504);
			g_sys[i_cent][i_jet]->GetYaxis()->SetTitleOffset(4);
			g_sys[i_cent][i_jet]->GetYaxis()->SetTitleFont(43);
			g_sys[i_cent][i_jet]->GetYaxis()->SetTitleSize(12);
			g_sys[i_cent][i_jet]->SetLineColor(h_nom[i_cent][i_jet]->GetMarkerColor());
			g_sys[i_cent][i_jet]->SetFillColorAlpha(g_sys[i_cent][i_jet]->GetFillColor(),opacity);
			g_sys[i_cent][i_jet]->SetLineWidth(1.);
			g_stat[i_cent][i_jet]->SetLineWidth(2.);

			h_nom[i_cent][i_jet]->GetXaxis()->SetRangeUser(0, 0.8);
			h_nom[i_cent][i_jet]->GetYaxis()->SetRangeUser(y_range_lo, y_range_hi);
			h_nom[i_cent][i_jet]->GetYaxis()->SetNdivisions(505);

			if (jet_itr == 0) g_sys[i_cent][i_jet]->Draw("a E2");
			else g_sys[i_cent][i_jet]->Draw("E2 same");
			g_stat[i_cent][i_jet]->Draw("P E1");

			jet_itr++;

		} // end trk loop

		canvas->cd(i_cent+1);
		line->DrawLine(line_x1, line_y1, line_x2, line_y2);
		legend->Draw();

		double x_left = 0.19, x_right = 0.93, y = 0.87, y_diff = 0.045;
		ltx->SetTextAlign(11);
		ltx->DrawLatexNDC(x_left, y, Form("%s", centrality.c_str()));
		ltx->SetTextAlign(31);
		ltx->DrawLatexNDC(x_right, y, "#scale[1.5]{#font[72]{ATLAS} Internal}");
		ltx->DrawLatexNDC(x_right, y=y-y_diff, "Pb+Pb #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}");
		ltx->DrawLatexNDC(x_right, y=y-y_diff, "#it{pp} #sqrt{#font[12]{s}} = 5.02 TeV, 25 pb^{-1}");
		ltx->DrawLatexNDC(x_right+0.3, y, Form("anti-#font[12]{k}_{#font[12]{t}} R=0.4"));
		first_pass_cent = false;


//		canvas->Print(Form("compiledCode/finalFig/RDpT_dR_jet%i_cent%i.pdf", i_jet, i_cent));

	} //end cent loop

	canvas->Print(Form("output_pdf_nominal/conf/%s_%s_dR.pdf", mode.c_str(), integType.c_str()));


	delete canvas;
	delete legend;
	delete line;
	delete ltx;
}


void integConfClass::cleanUp()
{
	for (int i_cent = 0; i_cent < 6; i_cent++)
	{
		for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
		{
			delete g_sys[i_cent][i_jet];
			delete g_stat[i_cent][i_jet];
		}
	}

}


void integConfClass::setSpecifics()
{
	if (mode.compare("DeltaDpT") == 0 && integType.compare("jetshape") == 0)
	{
		legend_x1 = 0.55, legend_y1 = 0.50, legend_x2 = legend_x1+0.25, legend_y2 = legend_y1+0.25;
		line_x1 = 0.0, line_x2 = 0.8, line_y1 = 0.0, line_y2 = 0.0;
		y_range_lo = -1, y_range_hi = 5;
		axis_label_y = jetshape_deltadpt_title;
		axis_label_x = r_title;
	}
	if (mode.compare("DeltaDpT") == 0 && integType.compare("lowpt_integ") == 0)
	{
		legend_x1 = 0.55, legend_y1 = 0.50, legend_x2 = legend_x1+0.25, legend_y2 = legend_y1+0.25;
		line_x1 = 0.0, line_x2 = 0.8, line_y1 = 0.0, line_y2 = 0.0;
		y_range_lo = -1, y_range_hi = 5;
		axis_label_y = lowpt_integ_deltadpt_title;
		axis_label_x = r_title;
	}


	if (mode.compare("RDpT") == 0 && integType.compare("jetshape") == 0)
	{
		legend_x1 = 0.18, legend_y1 = 0.55, legend_x2 = legend_x1+0.25, legend_y2 = legend_y1+0.25;
		line_x1 = 0.0, line_x2 = 0.8, line_y1 = 1.0, line_y2 = 1.0;
		y_range_lo = 0.5, y_range_hi = 3;
		axis_label_y = jetshape_rdpt_title;
		axis_label_x = r_title;
	}
	if (mode.compare("RDpT") == 0 && integType.compare("lowpt_integ") == 0)
	{
		legend_x1 = 0.18, legend_y1 = 0.55, legend_x2 = legend_x1+0.25, legend_y2 = legend_y1+0.25;
		line_x1 = 0.0, line_x2 = 0.8, line_y1 = 1.0, line_y2 = 1.0;
		y_range_lo = 0.0, y_range_hi = 5;
		axis_label_y = lowpt_integ_rdpt_title;
		axis_label_x = r_title;
	}



}
