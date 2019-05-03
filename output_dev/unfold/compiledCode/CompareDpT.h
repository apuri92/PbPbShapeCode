#ifndef CompareDpT_H
#define CompareDpT_H

#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TH2.h>
#include <TH3.h>
#include <TAxis.h>
#include <TMath.h>

#include <string>
#include <iostream>
#include <vector>
#include <map>


class CompareDpT
{
private:
	int NTrk, NCent, NJet, NdR;

public:
	std::string dataType;
	TFile *inputFile_pp = nullptr;
	TFile *inputFile_PbPb = nullptr;
	TFile *output_file = nullptr;
	std::string histName;

	TAxis* dR_binning;
	TAxis* jetpT_binning;
	TAxis* trkpT_binning;

	std::vector< std::vector<std::vector<TH1*>>>	spect_PbPb;
	std::vector< std::vector<TH1*>>					spect_pp;
	std::vector< std::vector<std::vector<TH1*>>>	spect_rdpt;
	std::vector< std::vector<std::vector<TH1*>>>	spect_deltadpt;


	std::vector<std::vector<TH1*>>	spect_lowpt_integ_PbPb;
	std::vector<TH1*>				spect_lowpt_integ_pp;
	std::vector<std::vector<TH1*>>	spect_jetshape_PbPb;
	std::vector<TH1*>				spect_jetshape_pp;


	std::vector<std::vector<TH1*>> spect_lowpt_integ_rdpt;
	std::vector<std::vector<TH1*>> spect_lowpt_integ_deltadpt;
	std::vector<std::vector<TH1*>> spect_jetshape_rdpt;
	std::vector<std::vector<TH1*>> spect_jetshape_deltadpt;

	std::string rdptr_title = "#it{R}_{#it{D} (#it{p}_{T}, #it{r})}";
	std::string deltadptr_title = "#it{#Delta} #it{D} (#it{p}_{T}, #it{r}) [GeV]";
	std::string dptr_title = "#it{D} (#it{p}_{T}, #it{r}) [GeV^{-1}]";
	std::string r_title = "#it{r}";
	std::string trk_title = "#it{p}_{T}^{trk} [GeV]";
	std::string jet_title = "#it{p}_{T}^{jet} [GeV]";

	std::string jetshape_rdpt_title = "#int_{0}^{r'}#int_{1}^{4} D_{PbPb} d#it{p}_{T}dr / #int_{0}^{r'} #int_{1}^{4}D_{pp} d#it{p}_{T}dr";
	std::string jetshape_deltadpt_title = "#int_{0}^{r'}#int_{1}^{4}#DeltaD d#it{p}_{T} dr";

	std::string lowpt_integ_rdpt_title = "#int_{1}^{4} D_{PbPb} d#it{p}_{T} / #int_{1}^{4}D_{pp} d#it{p}_{T}";
	std::string lowpt_integ_deltadpt_title = "#int_{1}^{4}#DeltaD d#it{p}_{T}";


	std::string lowpt_integ_dptr_title = Form("#int_{1}^{4} %s d#it{p}_{T}", dptr_title.c_str());
	std::string jetshape_dptr_title = Form("#int_{0}^{r'} #int_{1}^{4} %s d#it{p}_{T} dr", dptr_title.c_str());


	CompareDpT(std::string _dataType, std::string _histName, int config)
	{
		dataType = _dataType;
		histName = _histName;

		loadAxes();
		loadFiles(config);
		NCent = 6;

		getDpT();
	};

	~CompareDpT()
	{
		inputFile_pp->Close();
		inputFile_PbPb->Close();
		output_file->Close();

		delete (inputFile_pp);
		delete (inputFile_PbPb);
		delete (output_file);
	}

	void printInfo();
	void loadAxes();
	void writeToFile();
	void loadFiles(int);
	
	void getDpT();
	void getRDpT();
	void getDeltaDpT();
	void getJetShape();

	void labelvHist(TH1* hist, std::string x_label, std::string y_label);
	void labelvHist(std::vector<TH1*> hist, std::string x_label, std::string y_label);
	void labelvHist(std::vector<std::vector<TH1*>> hist, std::string x_label, std::string y_label);
	void labelvHist(std::vector<std::vector<std::vector<TH1*>>> hist, std::string x_label, std::string y_label);


	TH1* cumulateR(TH1*);
	double getArea(double r) {return TMath::Pi()*pow(r,2);}
	double getArea(double r1, double r2) {return TMath::Pi()*(pow(r2,2)-pow(r1,2));}


};

#endif
