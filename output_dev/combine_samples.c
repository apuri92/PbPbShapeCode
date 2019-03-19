//Takes in 5 jz samples and combines them to give weighted_output.

#include "functions/global_variables.h"
//#include "combine_eff_dr_jety_reduced.c"
#include "combine_eff_trketa.c"
//#include "combine_eff_dr.c"


void combine_samples()
{
	cout << endl << "	##############	Reading config	##############"<< endl;
	TEnv *m_config = new TEnv();
	m_config->ReadFile("perf_config.cfg", EEnvLevel(1));
	m_config->Print();

	int config = 0; config = m_config->GetValue("config", config);
	std::string dataset_type = "PbPb"; dataset_type = m_config->GetValue("dataset_type", dataset_type.c_str());
	std::string id = ""; id = m_config->GetValue("id", id.c_str());
	std::string tracking_cut = "ppTight"; tracking_cut = m_config->GetValue("tracking_cut", tracking_cut.c_str());
	bool eff_jetpt_jety = 0; eff_jetpt_jety = m_config->GetValue("module_eff_jetpt_jety", eff_jetpt_jety);
	bool eff_trketa = 0; eff_trketa = m_config->GetValue("module_eff_trketa", eff_trketa);

	cout << "	##############	Config done	##############" << endl << endl;

	string filename[5] = {"","","","",""};

	double FilterEff[5], CrossSec[5], SumJetW[5];
	get_weights(dataset_type, FilterEff, CrossSec, SumJetW);

	filename[1] = Form("raw_results/perf_c%i/Perf_MC_JZ2_out_histo_%s_5p02_r001_%s.root", config, dataset_type.c_str(), tracking_cut.c_str());
	filename[2] = Form("raw_results/perf_c%i/Perf_MC_JZ3_out_histo_%s_5p02_r001_%s.root", config, dataset_type.c_str(), tracking_cut.c_str());
	filename[3] = Form("raw_results/perf_c%i/Perf_MC_JZ4_out_histo_%s_5p02_r001_%s.root", config, dataset_type.c_str(), tracking_cut.c_str());

	filename[1] = Form("../hist-testdir.root", config, dataset_type.c_str(), tracking_cut.c_str());
	filename[2] = Form("../hist-testdir.root", config, dataset_type.c_str(), tracking_cut.c_str());
	filename[3] = Form("../hist-testdir.root", config, dataset_type.c_str(), tracking_cut.c_str());


	//setting up vectors (note that vectors correspond to number of files (not 5)
	std::vector<TFile*> theFiles;
	std::vector<double> nEvents;
	std::vector<double> weights;

	int nFile = 0;
	for (int i = 0; i<5; i++)
	{
		if( filename[i].find("root") != std::string::npos )
		{
			theFiles.push_back( new TFile(filename[i].c_str()) );
			nEvents.push_back(((TH1*)theFiles.at(nFile)->Get("EventLoop_EventCount"))->GetBinContent(1));
			weights.push_back((FilterEff[i]*CrossSec[i])/(nEvents.at(nFile)*SumJetW[i]));

			cout << Form("%s		%1.2e	%1.2e",theFiles.at(nFile)->GetName(), nEvents.at(nFile), weights.at(nFile)) << endl;

			nFile++;
		}
	}

//	if (1) combine_eff_dr_jety_reduced(theFiles, weights);
	if (1) combine_eff_trketa(theFiles, weights);
//	if (0) combine_eff_dr(theFiles, weights);



}
