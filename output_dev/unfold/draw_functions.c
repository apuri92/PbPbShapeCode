#include "../functions/global_variables.h"
#include "TEnv.h"
#include "TGaxis.h"



void setup_canvas(TCanvas* &c, string dataset_type, int itr)
{
//	c->Clear();
//	if (dataset_type == "pp") c->Divide(1,2);
//	if (dataset_type == "PbPb") c->cd(itr+1)->Divide(1,2);

	if (dataset_type == "pp") c->cd()->cd(1);
	if (dataset_type == "PbPb") c->cd(itr+1)->cd(1);
	gPad->SetPad(0,0.40,0.95,0.95);
	gPad->SetTopMargin(0.05);
	gPad->SetBottomMargin(0);
	gPad->SetRightMargin(0);


	if (dataset_type == "pp") c->cd()->cd(2);
	if (dataset_type == "PbPb") c->cd(itr+1)->cd(2);
	gPad->SetPad(0,0.0,0.95,0.40);
	gPad->SetTopMargin(0);
	gPad->SetBottomMargin(0.30);
	gPad->SetRightMargin(0);
}
