
void PlotLabels_AtlasSim_q2_5(double x = 0.2, double y = 0.8, int size = 24, int align = 11, bool split=true)
{
	TLatex *   tex = new TLatex();


	if (split) tex = new TLatex(x,y,"#splitline{#bf{#it{ATLAS Simulation}}}{Internal}");
	if (!split) tex = new TLatex(x,y,"#bf{#it{ATLAS Simulation}} Internal");


	tex->SetNDC();
	tex->SetTextFont(43);

	tex->SetTextSize(size);
	tex->SetTextAlign(align);

	tex->Draw();


	double diff = 0;
	if (size>24 && split) diff = 0.08;
	else if (size>24 && !split) diff = 0.05;
	else if (size<24 && split) diff = 0.103;
	else if (size<24 && !split) diff = 0.08;


	tex = new TLatex(x,y-diff,"Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV");
	tex->SetNDC();
	tex->SetTextFont(43);
	tex->SetTextSize(size - 2);
	tex->SetTextAlign(align);
	tex->Draw();


//	tex = new TLatex(x+0.12,y-diff,"#sqrt{#font[12]{s_{NN}}} = 5.02 TeV");
//	tex->SetNDC();
//	tex->SetTextFont(43);
//	tex->SetTextSize(size);
//	tex->SetTextAlign(align);
//	tex->Draw();

}
void PlotLabels_AtlasPrel_q2_5(double x = 0.2, double y = 0.8, int size = 24, int align = 11, bool split=true)
{
	TLatex *   tex = new TLatex();

	if (split) tex = new TLatex(x,y,"#splitline{#bf{#it{ATLAS Simulation}}}{Preliminary}");
	if (!split) tex = new TLatex(x,y,"#bf{#it{ATLAS Simulation}} Preliminary");


	tex->SetNDC();
	tex->SetTextFont(43);

	tex->SetTextSize(size);
	tex->SetTextAlign(align);

	tex->Draw();


	double diff = 0;
	if (size>24 && split) diff = 0.08;
	else if (size>24 && !split) diff = 0.05;
	else if (size<24 && split) diff = 0.103;
	else if (size<24 && !split) diff = 0.08;


	tex = new TLatex(x,y-diff,"Pb+Pb  #sqrt{#font[12]{s_{NN}}} = 5.02 TeV");
	tex->SetNDC();
	tex->SetTextFont(43);
	tex->SetTextSize(size - 2);
	tex->SetTextAlign(align);
	tex->Draw();


	//	tex = new TLatex(x+0.12,y-diff,"#sqrt{#font[12]{s_{NN}}} = 5.02 TeV");
	//	tex->SetNDC();
	//	tex->SetTextFont(43);
	//	tex->SetTextSize(size);
	//	tex->SetTextAlign(align);
	//	tex->Draw();
	
}

void PlotLabels_AtlasInt_q2_5(double x = 0.2, double y = 0.8, int size = 24, int align = 11)
{
	TLatex *   tex = new TLatex();

	tex = new TLatex(x,y,"#it{#bf{ATLAS} Internal}");
	tex->SetNDC();
	tex->SetTextFont(43);

	tex->SetTextSize(size);
	tex->SetTextAlign(align);

	tex->Draw();

	double diff = 0;
	if (size>24) diff = 0.06;
	else diff = 0.08;

	tex = new TLatex(x,y-diff,"#sqrt{#font[12]{s_{NN}}} = 5.02 TeV");
	tex->SetNDC();
	tex->SetTextFont(43);
	tex->SetTextSize(size - 4);
	tex->SetTextAlign(align);
	tex->Draw();

}

TGraphAsymmErrors* smallify(TGraphAsymmErrors* graph)
{
	graph->GetXaxis()->SetLabelSize(12);
	graph->GetXaxis()->SetTitleOffset(2.1);
	graph->GetXaxis()->SetTitleSize(14);
	graph->GetYaxis()->SetLabelSize(12);
	graph->GetYaxis()->SetTitleOffset(2.0);
	graph->GetYaxis()->SetTitleSize(14);
	graph->SetMarkerSize(0.9);

	return graph;
}





