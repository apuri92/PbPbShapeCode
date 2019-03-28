void get_weights(std::string dataset_type, double *FilterEff, double *CrossSec, double *SumJetW)
{

	double CS[5], FE[5], SUM_JET_WEIGHT[5];

//	if (dataset_type == "pp") //5p02 TeV pp MC Pythia8
//	{
//		CS[0]=6.7890E+07;
//		CS[1]=6.3997E+05;
//		CS[2]=4.7194E+03;
//		CS[3]=2.6602E+01;
//		CS[4]=2.2476E-01;
//
//		FE[0]=2.8320E-03;
//		FE[1]=4.2765E-03;
//		FE[2]=5.2883E-03;
//		FE[3]=4.5858E-03;
//		FE[4]=2.1831E-03;
//
//		SUM_JET_WEIGHT[0]=1;
//		SUM_JET_WEIGHT[1]=1;
//		SUM_JET_WEIGHT[2]=1;
//		SUM_JET_WEIGHT[3]=1;
//		SUM_JET_WEIGHT[4]=1;
//	}
//
//	if (dataset_type == "PbPb") //5p02 TeV PbPb MC Pythia8Powheg Overlay
	{
		CS[0]=1.2794E+08;
		CS[1]=1.9648E+07;
		CS[2]=5.7613E+05;
		CS[3]=4.1522E+04;
		CS[4]=8.4338E+02;

		FE[0]=1.5857E-03;
		FE[1]=1.2948E-04;
		FE[2]=4.2129E-05;
		FE[3]=2.8563E-06;
		FE[4]=5.9854E-07;

		SUM_JET_WEIGHT[0]=3.428918e+08;
		SUM_JET_WEIGHT[1]=4.656080e+06;
		SUM_JET_WEIGHT[2]=5.805497e+04;
		SUM_JET_WEIGHT[3]=7.924295e+02;
		SUM_JET_WEIGHT[4]=1.764195e+01;
	}


	for (int i = 0; i<5; i++)
	{
		FilterEff[i] = FE[i];
		CrossSec[i] = CS[i];
		SumJetW[i] = SUM_JET_WEIGHT[i];
	}

}

std::string num_to_cent(int scheme, int i)
{
	std::string cent_label;

	if (scheme == 31)
	{
		if (i==0) cent_label = "0 - 10%";
		if (i==1) cent_label = "10 - 20%";
		if (i==2) cent_label = "20 - 30%";
		if (i==3) cent_label = "30 - 40%";
		if (i==4) cent_label = "40 - 60%";
		if (i==5) cent_label = "60 - 80%";
		if (i==6) cent_label = "Inclusive";
	}

	return cent_label;
}

