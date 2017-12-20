//Takes in 5 jz samples and combines them to give weighted_output.

#include "extras/global_variables.h"
#include "combine_eff_jetpt_jety.c"
#include "combine_eff_trketa.c"

bool JX[5]={false,false,false,false,false};
double w[5]={0,0,0,0,0};

void combine_samples(bool isPerf, bool isMC, string cut)
{
	string filename[5];

	filename[0] = "";
	filename[1] = "";
	filename[2] = "";
	filename[3] = "";
	filename[4] = "";

//	filename[1] = "input/jz2_perf_pptight.root";
//	filename[2] = "input/jz3_perf_pptight.root";
//	filename[3] = "input/jz4_perf_pptight.root";

	filename[1] = "input/Perf_MC_JZ2_out_histo_PbPb_5p02_r001_comb.root";
	filename[2] = "input/Perf_MC_JZ3_out_histo_PbPb_5p02_r001_comb.root";
	filename[3] = "input/Perf_MC_JZ4_out_histo_PbPb_5p02_r001_comb.root";



	std::vector<TFile*> theFiles;

	for (int i = 0; i<5; i++)
	{
		if( filename[i].find("root") != std::string::npos )
		{
			theFiles.push_back( new TFile(filename[i].c_str()) );
			cout << filename[i] << endl;
			nFiles++;
			JX[i]=true;
		}
	}
	printf("Number of Files: %i \n",nFiles);

	//Number of evetns
	TH1* h1_tmp;
	int TotalNEvents[5];
	for (int i=0;i<nFiles;i++)
	{
		h1_tmp = (TH1*)theFiles[i]->Get("EventLoop_EventCount");
		TotalNEvents[i] = h1_tmp->GetBinContent(1);
		cout << "TotalNEvents "<< " " << i << " " << TotalNEvents[i] <<endl;
		delete h1_tmp;
	}

	//calculate weights
	for (int i=0;i<nFiles;i++)
	{
		int tmp_i = i+1;
		if (isMC) w[i]=(filterEff[tmp_i]*CS[tmp_i])/(TotalNEvents[i]*sum_jet_weights[tmp_i]);
		else w[i]=(1./TotalNEvents[i]);
		cout << Form("%i -> CS: %f, N: %i, fe: %f", i, CS[tmp_i], TotalNEvents[i], filterEff[tmp_i]) << endl;

		cout << Form("w[%i]: %1.2E",i, w[i]) << endl;
	}
	cout << endl;

	combine_eff_jetpt_jety(theFiles, w, cut);
	combine_eff_trketa(theFiles, w, cut);

}
