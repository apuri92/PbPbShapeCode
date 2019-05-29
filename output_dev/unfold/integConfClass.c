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
				r_position = h_nom[i_cent][i_jet]->GetBinCenter(i_dR+1);
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

	canvas = new TCanvas("canvas","canvas",800,600);
	legend = new TLegend(legend_x1, legend_y1, legend_x2, legend_y2, "","brNDC");
	legend->SetTextFont(43);
	legend->SetBorderSize(0);
	legend->SetTextSize(20);
	legend->SetNColumns(legend_cols);
	legend->SetFillStyle(0);
	ltx = new TLatex();
	ltx->SetTextFont(43);
	ltx->SetTextSize(20);
	line = new TLine();
	line->SetLineStyle(3);

	string trk_label = Form("%1.1f < #it{p}_{T} < %1.1f GeV", 1.0, 4.0);

	bool first_pass_cent = true;

	for (int i_cent = 0; i_cent < 6; i_cent++)
	{
		if (i_cent != 0 && i_cent != 3 && i_cent != 5) continue;

		canvas->Clear();
		string centrality = num_to_cent(31,i_cent);
		canvas->cd(i_cent+1);

		int jet_itr = 0;
		for (int i_jet = jet_pt_start; i_jet < jet_pt_end; i_jet++)
		{
			string jet_label = Form("%1.0f < #it{p}_{T}^{jet} < %1.0f GeV", jetpT_binning->GetBinLowEdge(i_jet+1), jetpT_binning->GetBinUpEdge(i_jet+1));

			SetHStyle_smallify(h_nom[i_cent][i_jet], jet_itr, 0);
			SetHStyle_graph_smallify(g_sys[i_cent][i_jet], jet_itr, 0);
			SetHStyle_graph_smallify(g_stat[i_cent][i_jet], jet_itr, 0);

			if (first_pass_cent) legend->AddEntry(g_sys[i_cent][i_jet],jet_label.c_str(),"p");

			g_sys[i_cent][i_jet]->SetLineColor(h_nom[i_cent][i_jet]->GetMarkerColor());
			g_sys[i_cent][i_jet]->SetFillColorAlpha(g_sys[i_cent][i_jet]->GetFillColor(),opacity);
			g_sys[i_cent][i_jet]->SetLineWidth(1.);
			g_stat[i_cent][i_jet]->SetLineWidth(1.);

			TGraphAsymmErrors* tmp_sys = (TGraphAsymmErrors*)g_sys[i_cent][i_jet]->Clone("tmp_gsys");
			TGraphAsymmErrors* tmp_stat = (TGraphAsymmErrors*)g_stat[i_cent][i_jet]->Clone("tmp_gsys");
			double shift_size = 0.008;
			tmp_sys = shift(tmp_sys, jet_itr, shift_size);
			tmp_stat = shift(tmp_stat, jet_itr, shift_size);

			if (jet_itr == 0) tmp_sys->Draw("a E2");
			else tmp_sys->Draw("E2 same");
			tmp_stat->Draw("P E1");

			jet_itr++;

		} // end trk loop

		line->DrawLine(line_x1, line_y1, line_x2, line_y2);
		legend->Draw();

		double x_left = 0.19, x_right = 0.92, y = 0.88, y_diff = 0.045;


		if (mode.compare("DeltaDpT") == 0)
		{
			ltx->SetTextAlign(31);
			x_left = x_right;
		}
		else ltx->SetTextAlign(11);
		ltx->DrawLatexNDC(x_left, y, "#scale[1.5]{#font[72]{ATLAS} Internal}");
		ltx->DrawLatexNDC(x_left, y=y-y_diff, "Pb+Pb #sqrt{#font[12]{s_{NN}}} = 5.02 TeV, 0.49 nb^{-1}");
		ltx->DrawLatexNDC(x_left, y=y-y_diff, "#it{pp} #sqrt{#font[12]{s}} = 5.02 TeV, 25 pb^{-1}");
		ltx->DrawLatexNDC(x_left, y=y-y_diff, Form("anti-#font[12]{k}_{#font[12]{t}} R=0.4"));
		ltx->DrawLatexNDC(x_left, y=y-y_diff, Form("%s", centrality.c_str()));
		first_pass_cent = false;

		canvas->Print(Form("output_pdf_nominal/conf/%s_%s_cent%i.pdf", mode.c_str(), integType.c_str(), i_cent));
	} //end cent loop



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
		legend_x1 = 0.63, legend_y1 = 0.45, legend_x2 = legend_x1+0.25, legend_y2 = legend_y1+0.25;
		line_x1 = 0.0, line_x2 = 0.8, line_y1 = 0.0, line_y2 = 0.0;
		y_range_lo = -1, y_range_hi = 5;
		axis_label_y = jetshape_deltadpt_title;
		axis_label_x = r_title;
	}
	if (mode.compare("DeltaDpT") == 0 && integType.compare("lowpt_integ") == 0)
	{
		legend_x1 = 0.63, legend_y1 = 0.45, legend_x2 = legend_x1+0.25, legend_y2 = legend_y1+0.25;
		line_x1 = 0.0, line_x2 = 0.8, line_y1 = 0.0, line_y2 = 0.0;
		y_range_lo = -1, y_range_hi = 5;
		axis_label_y = lowpt_integ_deltadpt_title;
		axis_label_x = r_title;
	}


	if (mode.compare("RDpT") == 0 && integType.compare("jetshape") == 0)
	{
		legend_x1 = 0.24, legend_y1 = 0.17, legend_x2 = legend_x1+0.600, legend_y2 = legend_y1+0.17;
		line_x1 = 0.0, line_x2 = 0.8, line_y1 = 1.0, line_y2 = 1.0;
		y_range_lo = 0., y_range_hi = 3.1;
//		y_range_lo = 0.6, y_range_hi = 1.9;
		axis_label_y = jetshape_rdpt_title;
		axis_label_x = r_title;
		legend_cols = 2;

	}
	if (mode.compare("RDpT") == 0 && integType.compare("lowpt_integ") == 0)
	{
		legend_x1 = 0.24, legend_y1 = 0.17, legend_x2 = legend_x1+0.600, legend_y2 = legend_y1+0.17;
		line_x1 = 0.0, line_x2 = 0.8, line_y1 = 1.0, line_y2 = 1.0;
		y_range_lo = 0., y_range_hi = 3.1;
//		y_range_lo = 0., y_range_hi = 2.9;
		axis_label_y = lowpt_integ_rdpt_title;
		axis_label_x = r_title;
		legend_cols = 2;
	}

}

TGraphAsymmErrors* integConfClass::shift(TGraphAsymmErrors* g, int variable, double shift_size = 0.0025)
{
	for (int i = 0; i < g->GetN(); i++)
	{
		double x, y, err_x;
		g->GetPoint(i, x, y);
		err_x = g->GetErrorX(i);
		double orig_pos = x+variable*shift_size;// * err_x/2;
		g->SetPoint(i, orig_pos, y);
	}

	g->GetXaxis()->SetLimits(0, 0.8);
	g->GetYaxis()->SetRangeUser(y_range_lo, y_range_hi);
	g->GetYaxis()->SetNdivisions(504);


	return g;
}
