#ifndef GLOBAL_VARIABLES_H
#define GLOBAL_VARIABLES_H

#include "TFile.h"
#include "TKey.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include "TPad.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include <iostream>
#include <TStyle.h>
#include "TMarker.h"
#include "TProfile.h"
#include "TLine.h"
#include "TText.h"
#include "TLatex.h"
#include "string.h"
#include <stdio.h>
#include <vector>
#include "TError.h"
#include "TAxis.h"
#include "TEnv.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "functions.c"

int n_cent_cuts = 7; //all bins, including inclusive

static const int N_JET_Y = 4;
double jet_y_binning[N_JET_Y+1] = {0, 0.3, 0.8, 1.2, 1.7};

#endif
