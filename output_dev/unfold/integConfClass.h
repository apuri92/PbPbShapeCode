#ifndef MANYHIST_H
#define MANYHIST_H
#include <string>
#include <iostream>
#include <vector>
#include <map>


class integConfClass
{

public:


	vector<vector<TH1*>> h_nom;
	vector<vector<TH1*>> h_sys_p;
	vector<vector<TH1*>> h_sys_n;
	vector<vector<TGraphAsymmErrors*>> g_sys;
	vector<vector<TGraphAsymmErrors*>> g_stat;
	string name;

	std::string rdptr_title = "#it{R}_{#it{D} (#it{p}_{T}, #it{r})}";
	std::string deltadptr_title = "#it{#Delta} #it{D} (#it{p}_{T}, #it{r}) [GeV]";
	std::string dptr_title = "#it{D} (#it{p}_{T}, #it{r}) [GeV^{-1}]";

	std::string r_title = "#it{r}";
	std::string trk_title = "#it{p}_{T}^{trk} [GeV]";
	std::string jet_title = "#it{p}_{T}^{jet} [GeV]";

	std::string jetshape_rdpt_title = "#it{R}_{#it{P(#it{r})}}";
	std::string jetshape_deltadpt_title = "#Delta_{#it{P(#it{r})}} [GeV^{-1}]";

	std::string lowpt_integ_rdpt_title = "#it{R}_{#it{#Theta(#it{r})}}";
	std::string lowpt_integ_deltadpt_title = "#Delta_{#it{#Theta(#it{r})}} [GeV^{-1}]";



	std::string lowpt_integ_dptr_title = Form("#int_{1}^{4} %s d#it{p}_{T}", dptr_title.c_str());
	std::string jetshape_dptr_title = Form("#int_{0}^{r'} #int_{1}^{4} %s d#it{p}_{T} dr", dptr_title.c_str());



	int jet_pt_start = 7;
	int jet_pt_end = 11;

	double opacity = 0.7;


	TAxis* dR_binning;
	TAxis* jetpT_binning;
	TAxis* trkpT_binning;
	int N_dR;
	int N_jetpt;
	int N_trkpt;

	TCanvas *canvas;
	TLegend *legend;
	TLatex *ltx;
	TLine* line;

	double legend_x1 = 0.0, legend_x2 = 0.1, legend_y1 = 0.0, legend_y2 = 0.1;
	int legend_cols = 1;
	double line_x1 = 0.0, line_x2 = 0.8, line_y1 = 0.0, line_y2 = 0.0;
	double y_range_lo = 0;
	double y_range_hi = 4;
	string axis_label_x = "", axis_label_y = "";

	TFile *f_nom;
	TFile *f_sys;

	string integType;
	string mode;
	string yaxis_label;
	integConfClass(string _mode, string _integType)
	{
		integType = _integType;
		mode = _mode;

		f_nom = new TFile(Form("compiledCode/root/output_nominal_data.root"));
		f_sys = new TFile(Form("compiledCode/root/final_%s_%s_sys_data.root",mode.c_str(), integType.c_str() ));


		setSpecifics();
		initHist();
		makeGraph();
		drawAll();
	}


	~integConfClass()
	{
		cleanUp();
		delete f_nom;
		delete f_sys;
	}


	void initHist();
	void makeGraph();
	void drawAll();
	void cleanUp();
	void setSpecifics();
	TGraphAsymmErrors* shift(TGraphAsymmErrors* g, int variable, double shift_size = 0.0025);


};


#endif
