
void set_plot_style()
{

	gStyle->SetOptStat(0);
	const Int_t NRGBs = 5;
	const Int_t NCont = 255;

	Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
	Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
	Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
	Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
	TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
	gStyle->SetNumberContours(NCont);
}

void pad_size(TVirtualPad * pad, double top = 0, double bot = 0.01, double lef = 0.05, double rig = -0.08)
{

	double xlo, ylo, xhi, yhi;
	gPad->GetPadPar(xlo, ylo, xhi, yhi);
	gPad->SetPad(xlo-0.005,ylo-0.006,xhi+0.01,yhi+0.006);

	double top_margin = pad->GetTopMargin() + top;
	double bot_margin = pad->GetBottomMargin() + bot;
	double lef_margin = pad->GetLeftMargin() + lef;
	double rig_margin = pad->GetRightMargin() + rig;
//	pad->SetMargin(lef_margin, rig_margin, bot_margin, top_margin);
	pad->SetTicks(1,1);

}




void hstyle(TH1* h1, int i)
{

	int color;

	h1->GetXaxis()->SetTitleOffset(1.3);

	if (i==0)
	{
		color = kBlack;
		h1->SetMarkerColor(color);
		h1->SetLineColor(color);
		h1->SetMarkerStyle(20);
	}
	else if(i==1)
	{
		color = kRed;
		h1->SetMarkerColor(color);
		h1->SetLineColor(color);
		h1->SetMarkerStyle(21);
	}
	else if(i==2)
	{
		color = kAzure-2;
		h1->SetMarkerColor(color);
		h1->SetLineColor(color);
		h1->SetMarkerStyle(22);
	}
	else if(i==3)
	{
		color = kGreen+2;
		h1->SetMarkerColor(color);
		h1->SetLineColor(color);
		h1->SetMarkerStyle(23);
	}
	else if(i==4)
	{
		color = kMagenta+3;
		h1->SetMarkerColor(color);
		h1->SetLineColor(color);
		h1->SetMarkerStyle(33);
	}
	else if(i==5)
	{
		color = kOrange-1;
		h1->SetMarkerColor(color);
		h1->SetLineColor(color);
		h1->SetMarkerStyle(34);
	}
	else if(i==6)
	{
		color = kBlack;
		h1->SetMarkerColor(color);
		h1->SetLineColor(color);
		h1->SetMarkerStyle(24);
	}

	else if(i==7)
	{
		color = kRed;
		h1->SetMarkerColor(color);
		h1->SetLineColor(color);
		h1->SetMarkerStyle(25);
	}

	else if(i==8)
	{
		color = kAzure-2;
		h1->SetMarkerColor(color);
		h1->SetLineColor(color);
		h1->SetMarkerStyle(26);
	}
	else if(i==9)
	{
		color = kGreen+2;
		h1->SetMarkerColor(color);
		h1->SetLineColor(color);
		h1->SetMarkerStyle(32);
	}
	else if(i==10)
	{
		color = kMagenta+3;
		h1->SetMarkerColor(color);
		h1->SetLineColor(color);
		h1->SetMarkerStyle(27);
	}

	else if(i==11)
	{
		color = kOrange-1;
		h1->SetMarkerColor(color);
		h1->SetLineColor(color);
		h1->SetMarkerStyle(28);
	}


}
//
void hstyle(TProfile* h1, int i)
{

	int color;

	h1->GetXaxis()->SetTitleOffset(1.3);

	if (i==0)
	{
		color = kBlack;
		h1->SetMarkerColor(color);
		h1->SetLineColor(color);
		h1->SetMarkerStyle(20);
	}
	else if(i==1)
	{
		color = kRed;
		h1->SetMarkerColor(color);
		h1->SetLineColor(color);
		h1->SetMarkerStyle(21);
	}
	else if(i==2)
	{
		color = kAzure-2;
		h1->SetMarkerColor(color);
		h1->SetLineColor(color);
		h1->SetMarkerStyle(22);
	}
	else if(i==3)
	{
		color = kGreen+2;
		h1->SetMarkerColor(color);
		h1->SetLineColor(color);
		h1->SetMarkerStyle(23);
	}
	else if(i==4)
	{
		color = kMagenta+3;
		h1->SetMarkerColor(color);
		h1->SetLineColor(color);
		h1->SetMarkerStyle(33);
	}
	else if(i==5)
	{
		color = kOrange-1;
		h1->SetMarkerColor(color);
		h1->SetLineColor(color);
		h1->SetMarkerStyle(34);
	}
	else if(i==6)
	{
		color = kBlack;
		h1->SetMarkerColor(color);
		h1->SetLineColor(color);
		h1->SetMarkerStyle(24);
	}

	else if(i==7)
	{
		color = kRed;
		h1->SetMarkerColor(color);
		h1->SetLineColor(color);
		h1->SetMarkerStyle(25);
	}

	else if(i==8)
	{
		color = kAzure-2;
		h1->SetMarkerColor(color);
		h1->SetLineColor(color);
		h1->SetMarkerStyle(26);
	}
	else if(i==9)
	{
		color = kGreen+2;
		h1->SetMarkerColor(color);
		h1->SetLineColor(color);
		h1->SetMarkerStyle(32);
	}
	else if(i==10)
	{
		color = kMagenta+3;
		h1->SetMarkerColor(color);
		h1->SetLineColor(color);
		h1->SetMarkerStyle(27);
	}

	else if(i==11)
	{
		color = kOrange-1;
		h1->SetMarkerColor(color);
		h1->SetLineColor(color);
		h1->SetMarkerStyle(28);
	}
	
	
}
void hstyle(TGraphAsymmErrors* g1, int i)
{

	int color;

//	g1->GetXaxis()->SetTitleOffset(1.3);

	if (i==0)
	{
		color = kBlack;
		g1->SetMarkerColor(color);
		g1->SetLineColor(color);
		g1->SetMarkerStyle(20);
	}
	else if(i==1)
	{
		color = kRed;
		g1->SetMarkerColor(color);
		g1->SetLineColor(color);
		g1->SetMarkerStyle(21);
	}
	else if(i==2)
	{
		color = kAzure-2;
		g1->SetMarkerColor(color);
		g1->SetLineColor(color);
		g1->SetMarkerStyle(22);
	}
	else if(i==3)
	{
		color = kGreen+2;
		g1->SetMarkerColor(color);
		g1->SetLineColor(color);
		g1->SetMarkerStyle(23);
	}
	else if(i==4)
	{
		color = kMagenta+3;
		g1->SetMarkerColor(color);
		g1->SetLineColor(color);
		g1->SetMarkerStyle(29);
	}
	else if(i==5)
	{
		color = kOrange-1;
		g1->SetMarkerColor(color);
		g1->SetLineColor(color);
		g1->SetMarkerStyle(33);
	}
	else if(i==6)
	{
		color = kCyan+3;
		g1->SetMarkerColor(color);
		g1->SetLineColor(color);
		g1->SetMarkerStyle(34);
	}

	else if(i==7)
	{
		color = kBlue+3;
		g1->SetMarkerColor(color);
		g1->SetLineColor(color);
		g1->SetMarkerStyle(3);
	}
}
//
//


void styleLabel(TH1* h)
{
	h->GetXaxis()->SetTitleFont(43);
	h->GetYaxis()->SetTitleFont(43);

	h->GetXaxis()->SetTitleSize(14);
	h->GetYaxis()->SetTitleSize(14);

	h->GetXaxis()->SetTitleOffset(1.5);
	h->GetYaxis()->SetTitleOffset(2.2);
}

void PlotLabels_AtlasSim_q2_5(double x = 0.2, double y = 0.8, int size = 24, int align = 11, bool split=true)
{
	TLatex *   tex;
	TPaveText *pt  = new TPaveText(.05,.1,.95,.8);


	if (split) tex = new TLatex(x,y,"#splitline{#bf{#it{ATLAS Simulation}}}{Internal}");
	if (!split) tex = new TLatex(x,y,"#bf{#it{ATLAS Simulation}} Internal");


	tex->SetNDC();
	tex->SetTextFont(43);

	tex->SetTextSize(size);
	tex->SetTextAlign(align);

	tex->Draw();


	double diff;
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
	TLatex *   tex;

	if (split) tex = new TLatex(x,y,"#splitline{#bf{#it{ATLAS Simulation}}}{Preliminary}");
	if (!split) tex = new TLatex(x,y,"#bf{#it{ATLAS Simulation}} Preliminary");


	tex->SetNDC();
	tex->SetTextFont(43);

	tex->SetTextSize(size);
	tex->SetTextAlign(align);

	tex->Draw();


	double diff;
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
	TLatex *   tex;

	tex = new TLatex(x,y,"#it{#bf{ATLAS} Internal}");
	tex->SetNDC();
	tex->SetTextFont(43);

	tex->SetTextSize(size);
	tex->SetTextAlign(align);

	tex->Draw();

	double diff;
	if (size>24) diff = 0.06;
	else diff = 0.08;

	tex = new TLatex(x,y-diff,"#sqrt{#font[12]{s_{NN}}} = 5.02 TeV");
	tex->SetNDC();
	tex->SetTextFont(43);
	tex->SetTextSize(size - 4);
	tex->SetTextAlign(align);
	tex->Draw();

}





void eight_pad_style(TH1* h1)			//for make_jer_plots_new.c
{

	h1->GetYaxis()->SetTitleSize(0.1);
	h1->GetYaxis()->SetTitleOffset(0.5);
	h1->GetYaxis()->SetLabelSize(0.08);
	h1->GetYaxis()->SetNdivisions(504);

	h1->GetXaxis()->SetTitleSize(0.09);
	h1->GetXaxis()->SetTitleOffset(1);
	h1->GetXaxis()->SetLabelSize(0.08);

	h1->SetTitle(0);
}

void eight_pad_style(TGraphAsymmErrors* h1)			//for make_jer_plots_new.c
{

	h1->GetYaxis()->SetTitleSize(0.1);
	h1->GetYaxis()->SetTitleOffset(0.5);
	h1->GetYaxis()->SetLabelSize(0.08);
	h1->GetYaxis()->SetNdivisions(504);

	h1->GetXaxis()->SetTitleSize(0.09);
	h1->GetXaxis()->SetTitleOffset(0.7);
	h1->GetXaxis()->SetLabelSize(0.08);

	h1->SetTitle(0);
}


void pad_eta_labels(double eta_lo, double eta_hi)
{

	string name;
	TLatex *label = new TLatex();
	label->SetNDC();
	label->SetTextFont(43);
	label->SetTextSize(22);
	label->SetTextAlign(22);
	label->SetTextSize(23);

	name = Form("%4.1f < |#eta| < %4.1f",eta_lo,eta_hi);
	label->DrawLatex(0.70,0.85,name.c_str());

}


void smallify(TH1* hist)
{
	hist->GetXaxis()->SetLabelFont(43);
	hist->GetYaxis()->SetLabelFont(43);
	hist->GetXaxis()->SetTitleFont(43);
	hist->GetYaxis()->SetTitleFont(43);

	hist->GetXaxis()->SetTitleOffset(2.5);
	hist->GetYaxis()->SetTitleOffset(2.5);

	hist->GetXaxis()->SetTitleSize(16);
	hist->GetYaxis()->SetTitleSize(16);

	hist->GetXaxis()->SetLabelSize(14);
	hist->GetYaxis()->SetLabelSize(14);

	hist->SetMarkerSize(0.9);

//	return hist;
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





