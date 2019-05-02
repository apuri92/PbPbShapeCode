#ifndef DRAWINGCLASS_H
#define DRAWINGCLASS_H

#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TH2.h>
#include <TH3.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLine.h>

#include <string>
#include <iostream>
#include <vector>
#include <map>


class drawingClass
{
private:
	TLegend* legend;
	TCanvas* canvas;
	TLatex *ltx;
	TLine *line;
	double line_x1 = 0, line_y1 = 1, line_x2 = 0.8, line_y2 = 1;
	double legend_x1 = 0.19, legend_y1 = 0.70, legend_x2 = 0.70, legend_y2 = 0.92;
	double legend_text_size = 10;
	double y_lo = 0, y_hi = 2;
	double canvas_x = 1200, canvas_y = 600;
	bool logx = false, logy = false;


public:
	int jetStart = 7;
	int jetEnd = 11;
	int trkStart = 2;
	int trkEnd = 9;
	int NCent = 6;
	int config;

	TAxis* dR_binning;
	TAxis* jetpT_binning;
	TAxis* trkpT_binning;

	drawingClass(int _config)
	{
		loadAxes();
		loadObjects();
		config = _config;
	}

	~drawingClass()
	{
		delete canvas;
		delete legend;
		delete ltx;
		delete line;
	}


	void loadAxes();
	void loadObjects();
	void drawCentPanels(std::vector<std::vector<std::vector<TH1*>>>, std::string);
	void drawCentPanels(std::vector<std::vector<TH1*>>, std::string);
	std::string trk_label(int bin);
	std::string jet_label(int bin);
	std::string dR_label(int bin);
	std::string cent_label(int bin);
	void setSpecifics(std::string);

};
#endif
