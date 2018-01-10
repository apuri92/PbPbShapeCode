string num_to_cent(int scheme, int i)
{
	gErrorIgnoreLevel = 2001;

	string name;

	if (scheme == 30)
	{
		if (i==0) name = " 0 - 10%";
		if (i==1) name = "10 - 20%";
		if (i==2) name = "20 - 30%";
		if (i==3) name = "30 - 40%";
		if (i==4) name = "40 - 60%";
		if (i==5) name = "60 - 80%";
		if (i==6) name = "Inclusive";
	}

	return name;
}

void get_weights(bool isPbPb, double *FilterEff, double *CrossSec, double *SumJetW)
{

    double CS[5], FE[5], SUM_JET_WEIGHT[5];
    
    if (!isPbPb) //5p02 TeV pp MC Pythia8
    {
        CS[0]=6.7890E+07;
        CS[1]=6.3997E+05;
        CS[2]=4.7194E+03;
        CS[3]=2.6602E+01;
        CS[4]=2.2476E-01;

        FE[0]=2.8320E-03;
        FE[1]=4.2765E-03;
        FE[2]=5.2883E-03;
        FE[3]=4.5858E-03;
        FE[4]=2.1831E-03;
        
        SUM_JET_WEIGHT[0]=1;
        SUM_JET_WEIGHT[1]=1;
        SUM_JET_WEIGHT[2]=1;
        SUM_JET_WEIGHT[3]=1;
        SUM_JET_WEIGHT[4]=1;
    }

    if (isPbPb) //5p02 TeV PbPb MC Pythia8Powheg Overlay
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

//string num_to_dr(int i)
//{
//	gErrorIgnoreLevel = 2001;
//
//	string name;
//
//	if (scheme == 30)
//	{
//		if (i==0) name = " 0 - 10%";
//		if (i==1) name = "10 - 20%";
//		if (i==2) name = "20 - 30%";
//		if (i==3) name = "30 - 40%";
//		if (i==4) name = "40 - 50%";
//		if (i==5) name = "50 - 60%";
//		if (i==6) name = "60 - 80%";
//		if (i==7) name = "Inclusive";
//	}
//
//	if (scheme == 31)
//	{
//		if (i==0) name = " 0 - 10%";
//		if (i==1) name = "10 - 20%";
//		if (i==2) name = "20 - 30%";
//		if (i==3) name = "30 - 40%";
//		if (i==4) name = "40 - 60%";
//		if (i==5) name = "60 - 80%";
//		if (i==6) name = "Inclusive";
//	}
//
//	return name;
//}


TH2* axis_cuts(TH3* orig, string name, double eta1, double eta2,string axis)
{
	gErrorIgnoreLevel = 2001;

	TH3* histo3D_1 = new TH3D(); //original
	TH3* histo3D_2 = new TH3D(); //positive
	TH3* histo3D_3 = new TH3D(); //negative
	TH3* histo3D_4 = new TH3D(); //final

	TH2* histo2D_2 = new TH2D(); //intermediate
	TH2* histo2D_3 = new TH2D();	//intermediate
	TH2* h_resp_yx = new TH2D();	//final

	TH2* sliced;
	
	histo3D_1=(TH3*)orig->Clone("histo3D_1");
	histo3D_2=(TH3*)histo3D_1->Clone("histo3D_2"); //positive
	histo3D_3=(TH3*)histo3D_1->Clone("histo3D_3"); //negative

	if (axis == "z")
	{
		histo3D_2->GetZaxis()->SetRangeUser(eta1,eta2);
		histo3D_3->GetZaxis()->SetRangeUser(-eta2,-eta1);
		histo2D_2 = (TH2*)histo3D_2->Project3D("yx");
		histo2D_3 = (TH2*)histo3D_3->Project3D("yx");

	}

	if (axis == "y")
	{
		histo3D_2->GetYaxis()->SetRangeUser(eta1,eta2);
		histo3D_3->GetYaxis()->SetRangeUser(-eta2,-eta1);
		histo2D_2 = (TH2*)histo3D_2->Project3D("zx");
		histo2D_3 = (TH2*)histo3D_3->Project3D("zx");
	}

	if (axis == "x")
	{
		histo3D_2->GetYaxis()->SetRangeUser(eta1,eta2);
		histo3D_3->GetYaxis()->SetRangeUser(-eta2,-eta1);
		histo2D_2 = (TH2*)histo3D_2->Project3D("zy");
		histo2D_3 = (TH2*)histo3D_3->Project3D("zy");
	}

	histo2D_3->Add(histo2D_2,1);

	sliced=(TH2*)histo2D_3->Clone(name.c_str());

	delete histo3D_1;
	delete histo3D_2;
	delete histo3D_3;
	delete histo3D_4;

	delete histo2D_2;
	delete histo2D_3;
	delete h_resp_yx;

	return sliced;
}

void apply_pt_ranges(TH1 *histo, double pt_low, double pt_high){
	for(int i=1; i<histo->GetNbinsX()+1; i++){
		double center = histo->GetBinCenter(i);
		double width = histo->GetBinWidth(i);
		double low_edge = center - width/2.0;
		double high_edge = center + width/2.0;
		if(center < pt_low)histo->SetBinContent(i,-1);
		if(high_edge > pt_high){
			histo->SetBinContent(i,-1);
			cout << high_edge << " " << pt_high << endl;
		}
	}
}

void apply_z_ranges(TH1 *histo, double z_low, double z_high = 1){
	for(int i=1; i<histo->GetNbinsX()+1; i++)
	{
		double center = histo->GetBinCenter(i);
		double width = histo->GetBinWidth(i);
		double low_edge = center - width/2.0;
		double high_edge = center + width/2.0;
		if(low_edge < z_low)
		{
			histo->SetBinContent(i,20);
		}
		if(high_edge > z_high) histo->SetBinContent(i,20);
	}
}
bool replace(std::string& str, const std::string& from, const std::string& to) {
	size_t start_pos = str.find(from);
	if(start_pos == std::string::npos)
		return false;
	str.replace(start_pos, from.length(), to);
	return true;
}

TH1* make_histo(TGraphAsymmErrors* g_efficiency)
{

	double bins[1000];

	int Nbins = g_efficiency->GetN();
	double x, y;



	for (int i = 0; i < Nbins; i++)
	{
		g_efficiency->GetPoint(i,x,y);
		bins[i] = x-g_efficiency->GetErrorXlow(i);
	}

	g_efficiency->GetPoint(Nbins-1,x,y);
	bins[Nbins] = x + g_efficiency->GetErrorXhigh(Nbins-1);

	TH1* histo = new TH1D("histo","histo",Nbins,bins);

	for (int i = 1; i <= histo->GetNbinsX(); i++)
	{
		g_efficiency->GetPoint(i-1,x,y);
		histo->SetBinContent(i,y);
	}

	string name = Form("Conf_Int_%s",g_efficiency->GetName());
	histo->SetName(name.c_str());
	histo->SetTitle(g_efficiency->GetTitle());

	return histo;
}

void SetUncertainties(TH1 *h)
{
	for (int i=1;i<=h->GetNbinsX();i++){
		h->SetBinContent(i,h->GetBinError(i));
		h->SetBinError(i,0.);
	}

}

void fine_eta_binning(double array[1000], int &num)
{
	num = 50;
	array[0] = -2.5;
	for (int i=1; i<=num; i++)
	{
		array[i] = array[0]+i*0.1;
	}
}

