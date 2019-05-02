#include "ManyHist.h"
#include "drawingClass.h"

#include <iostream>
#include "TCanvas.h"
#include <map>
#include "/Users/Akshat/Dropbox/.RootUtils/AtlasStyle.C"

using namespace std;


int main()
{
	SetAtlasStyle();
	gErrorIgnoreLevel = 3001;


	vector<int> configs = {0};//, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 42, 43, 44, 45};

//	vector<int> configs = {0, 1, 2};


	string histName = "h_ChPS_final_indR_trk%i_cent%i_jetpt%i";

	for (int i_cfg = 0; i_cfg < configs.size(); i_cfg++)
	{
		int config = configs[i_cfg];

		ManyHist *MH_data = new ManyHist("data",histName, config);
		MH_data->getRDpT();
		MH_data->getDeltaDpT();
		MH_data->getJetShape();
		MH_data->writeToFile();


//		drawingClass *DC = new drawingClass(config);
//		DC->drawCentPanels(MH_data->spect_rdpt,"rdpt");
//		DC->drawCentPanels(MH_data->spect_deltadpt,"deltadpt");
//
//		DC->drawCentPanels(MH_data->spect_lowpt_integ_rdpt,"lowpt_integ_rdpt");
//		DC->drawCentPanels(MH_data->spect_lowpt_integ_deltadpt,"lowpt_integ_deltadpt");
//		DC->drawCentPanels(MH_data->spect_jetshape_rdpt,"jetshape_rdpt");
//		DC->drawCentPanels(MH_data->spect_jetshape_deltadpt,"jetshape_deltadpt");
//		delete DC;

		delete MH_data;

		cout << "Done with cfg: " << config << endl << endl;
	}


	return 0;
}
