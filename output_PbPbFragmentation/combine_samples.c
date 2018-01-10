//Takes in 5 jz samples and combines them to give weighted_output.

#include "extras/global_variables.h"
#include "combine_eff_jetpt_jety.c"
#include "combine_eff_trketa.c"

bool JX[5]={false,false,false,false,false};
double w[5]={0,0,0,0,0};

void combine_samples(bool isPbPb = 1, string cut = "ppTight")
{
	string filename[5];

	filename[0] = "";
	filename[1] = "";
	filename[2] = "";
	filename[3] = "";
	filename[4] = "";

    if (isPbPb)
    {
        filename[1] = "input/Perf_MC_JZ2_out_histo_PbPb_5p02_r001_drmax_1p2.root";
        filename[2] = "input/Perf_MC_JZ3_out_histo_PbPb_5p02_r001_drmax_1p2.root";
        filename[3] = "input/Perf_MC_JZ4_out_histo_PbPb_5p02_r001_drmax_1p2.root";

//        filename[1] = "input/Perf_MC_JZ2_out_histo_PbPb_5p02_r001_drmax_0p4.root";
//        filename[2] = "input/Perf_MC_JZ3_out_histo_PbPb_5p02_r001_drmax_0p4.root";
//        filename[3] = "input/Perf_MC_JZ4_out_histo_PbPb_5p02_r001_drmax_0p4.root";

//        filename[1] = "input/Perf_MC_JZ2_out_histo_PbPb_5p02_r001_comb.root";
//        filename[2] = "input/Perf_MC_JZ3_out_histo_PbPb_5p02_r001_comb.root";
//        filename[3] = "input/Perf_MC_JZ4_out_histo_PbPb_5p02_r001_comb.root";
}
    if (!isPbPb)
    {
        filename[1] = "input/Perf_MC_JZ2_out_histo_pp_5p02_r001.root";
        filename[2] = "input/Perf_MC_JZ3_out_histo_pp_5p02_r001.root";
        filename[3] = "input/Perf_MC_JZ4_out_histo_pp_5p02_r001.root";

//        filename[1] = "input/Perf_MC_JZ2_out_histo_pp_5p02_r001_drmax_1p2.root";
//        filename[2] = "input/Perf_MC_JZ3_out_histo_pp_5p02_r001_drmax_1p2.root";
//        filename[3] = "input/Perf_MC_JZ4_out_histo_pp_5p02_r001_drmax_1p2.root";

//        filename[1] = "input/Perf_MC_JZ2_out_histo_pp_5p02_r001_drmax_0p4.root";
//        filename[2] = "input/Perf_MC_JZ3_out_histo_pp_5p02_r001_drmax_0p4.root";
//        filename[3] = "input/Perf_MC_JZ4_out_histo_pp_5p02_r001_drmax_0p4.root";
    }
    
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

    double FilterEff[5], CrossSec[5], SumJetW[5];
    get_weights(isPbPb, FilterEff, CrossSec, SumJetW);
    
	//calculate weights
	for (int i=0;i<nFiles;i++)
	{
		int tmp_i = i+1;
        w[i]=(FilterEff[tmp_i]*CrossSec[tmp_i])/(TotalNEvents[i]*SumJetW[tmp_i]);
		cout << Form("%i -> CS: %f, N: %i, fe: %f", i, CrossSec[tmp_i], TotalNEvents[i], FilterEff[tmp_i]) << endl;
		cout << Form("w[%i]: %1.2E",i, w[i]) << endl;
	}
	cout << endl;

	combine_eff_jetpt_jety(theFiles, w, cut);
//	combine_eff_trketa(isPbPb, theFiles, w, cut);

}
