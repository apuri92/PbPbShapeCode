/*
 Use output of ./analysis (raw_unfolded_*.root) and make all comparison plots
 output of this will be used by systematics.c and integratedSystematics.c to
 get the systematics and summary plots for the systematics.

 The int/paper level plots will be drawn separately.
 */


#include "CompareDpT.h"
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


	vector<int> configs = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 42, 43, 44, 45};
//	vector<int> configs = {45};

	string histName = "h_ChPS_final_indR_trk%i_cent%i_jetpt%i";

	for (int i_cfg = 0; i_cfg < configs.size(); i_cfg++)
	{
		int config = configs[i_cfg];

		CompareDpT *compare_data = new CompareDpT("data",histName, config);
		compare_data->getRDpT();
		compare_data->getDeltaDpT();
		compare_data->getJetShape();
		compare_data->writeToFile();


		drawingClass *DC = new drawingClass(config);
		DC->drawCentPanels(compare_data->spect_rdpt,"rdpt");
		DC->drawCentPanels(compare_data->spect_deltadpt,"deltadpt");

		DC->drawCentPanels(compare_data->spect_lowpt_integ_rdpt,"lowpt_integ_rdpt");
		DC->drawCentPanels(compare_data->spect_lowpt_integ_deltadpt,"lowpt_integ_deltadpt");
		DC->drawCentPanels(compare_data->spect_jetshape_rdpt,"jetshape_rdpt");
		DC->drawCentPanels(compare_data->spect_jetshape_deltadpt,"jetshape_deltadpt");
		delete DC;

		delete compare_data;

		cout << "Done with cfg: " << config << endl << endl;
	}


	return 0;
}
