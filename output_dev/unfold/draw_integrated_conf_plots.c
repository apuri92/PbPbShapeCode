#include "../functions/global_variables.h"
#include "integConfClass.c"

void draw_integrated_conf_plots()
{
	cout << "######### DOING INTCONF_Plots #########" << endl;

	SetAtlasStyle();
	gErrorIgnoreLevel = 3001;

	//indR
	integConfClass *DeltaDpT_lowptInteg = new integConfClass("DeltaDpT", "lowpt_integ");
	integConfClass *RDpT_lowptInteg = new integConfClass("RDpT", "lowpt_integ");
	integConfClass *DeltaDpT_jetshape = new integConfClass("DeltaDpT", "jetshape");
	integConfClass *RDpT_jetshape = new integConfClass("RDpT", "jetshape");

	cout << "######### DONE INTCONF PLOTS #########" << endl << endl;;


}

