//data
//double CS[5]={1,1,1,1,1};
//double filterEff[5] = {1,1,1,1,1};

//Data Overlay (DF)
//double CS[6]={1.2794E+08, 1.9648E+07, 5.7613E+05, 4.1522E+04, 8.4338E+02};
//double filterEff[6]={1.5857E-03, 1.2948E-04, 4.2129E-05, 2.8563E-06, 5.9854E-07};

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
#include "functions.c"
#include "hstyle.c"
#include "TError.h"
#include "TAxis.h"

static const int n_cent_cuts = 7 - 1;

double CS[5]={1.2794E+08, 1.9648E+07, 5.7613E+05, 4.1522E+04, 8.4338E+02};
double filterEff[5]={1.5857E-03, 1.2948E-04, 4.2129E-05, 2.8563E-06, 5.9854E-07};

double sum_jet_weights[5]={3.428918e+08, 4.656080e+06, 5.805497e+04, 7.924295e+02, 1.764195e+01};

int nFiles=0;
int vec_ptbins;

static const int centrality_scheme = 30;

#endif
