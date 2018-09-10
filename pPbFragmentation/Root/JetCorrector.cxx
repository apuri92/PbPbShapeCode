#include "pPbFragmentation/JetCorrector.h"


int JetCorrector::GetJetYBin(float y){
	int yBin=0;
	if (fabs(y) < 2.1) yBin =3;
	if (fabs(y) < 1.2) yBin =2;
	if (fabs(y) < 0.8) yBin =1;
	if (fabs(y) < 0.3) yBin =0;
	return yBin;
}

Int_t JetCorrector::GetJetEtaBin(float eta) {
	if (fabs(eta)<0.4) return 1;
	if (fabs(eta)<0.8) return 2;
	if (fabs(eta)<1.2) return 3;
	if (fabs(eta)<1.6) return 4;
	if (fabs(eta)<2.0) return 5;
	if (fabs(eta)<2.4) return 6;
	if (fabs(eta)<2.8) return 7;  
	else return 0;
}

Int_t JetCorrector::GetJetpTBin(float pt, TAxis* pt_bins) {

	int xBin = -1;
	for (int i = 0; i < pt_bins->GetNbins(); i++)
	{
		if (pt >= pt_bins->GetBinLowEdge(i+1) &&
			pt < pt_bins->GetBinUpEdge(i+1)) xBin = i;
        
        if (pt < pt_bins->GetBinLowEdge(1)) xBin = 0;
	}
	if (xBin == -1)
	{
		xBin = pt_bins->GetNbins() - 1;
		cout << Form("JetCorrector Warning: did not find bin for variable at %f on axis %s. Using %i bin.",pt, pt_bins->GetName(), xBin) << endl;
	}

	return xBin;

}

float JetCorrector::GetJetWeight(double pt, double eta, double phi)
{
	int xb=event_weight_histo->GetXaxis()->FindBin(pt);
	int yb=event_weight_histo->GetYaxis()->FindBin(eta);
	int zb=event_weight_histo->GetZaxis()->FindBin(phi);
	float jet_weight=event_weight_histo->GetBinContent(xb,yb,zb);

	return jet_weight;

}


float JetCorrector::GetFCalWeight(float FCalEt) {

    int centile = -1;
	float event_weight_fcal=1;
	for (int i=0; i<centiles.size(); i++)
	{
		if ( FCalEt < range_hi.at(i) && FCalEt >= range_lo.at(i) ) centile = centiles.at(i);
	}
	event_weight_fcal = cent_weight_histo->GetBinContent(cent_weight_histo->FindBin(centile));
	if (!m_isMB) event_weight_fcal= event_weight_fcal * FCal_HP_v_MB_weights_histo->GetBinContent(FCal_HP_v_MB_weights_histo->GetXaxis()->FindBin(FCalEt));
	return event_weight_fcal;			
}

float JetCorrector::GetFCalWeight(float FCalEt, int sample) {
	if (is_pp) return 1.0; //is using pp MC so no weighting required
	
	if (FCalEt>5.0) FCalEt=5.0;
	//MC+MB (overlay)->HP
	if (sample == 1) return FCal_HP_v_MCOV_weights_histo->GetBinContent(FCal_HP_v_MB_weights_histo->GetXaxis()->FindBin(FCalEt));
	//MB->HP
	else if (sample == 2) return FCal_HP_v_MB_weights_histo->GetBinContent(FCal_HP_v_MB_weights_histo->GetXaxis()->FindBin(FCalEt));
	//MB->MC + MB (Overlay)
	else if (sample == 3) return FCal_MCOV_v_MB_weights_histo->GetBinContent(FCal_HP_v_MB_weights_histo->GetXaxis()->FindBin(FCalEt));
	//MBOV->HP
	else if (sample == 4) return FCal_HP_v_MBOV_weights_histo->GetBinContent(FCal_HP_v_MB_weights_histo->GetXaxis()->FindBin(FCalEt));
	//MBOV->MC + MB (Overlay)
	else if (sample == 5) return FCal_MCOV_v_MBOV_weights_histo->GetBinContent(FCal_HP_v_MB_weights_histo->GetXaxis()->FindBin(FCalEt));
}

/*
float JetCorrector::GetFCalHPWeight(float FCalEt) {
	float event_weight_fcal= FCal_HP_v_MB_weights_histo->GetBinContent(FCal_HP_v_MB_weights_histo->GetXaxis()->FindBin(FCalEt));
	return event_weight_fcal;			
}
*/
bool JetCorrector::MCJetJERClean(float truth_jet_pt,float reco_jet_pt, float truth_jet_eta, int cent){
	bool pass = true;
	float JER =  GetJER(truth_jet_pt, truth_jet_eta, cent);
	if (fabs (( truth_jet_pt - reco_jet_pt) /  truth_jet_pt ) > JERcut * JER ) pass = false;
	return pass;
	
}

float JetCorrector::GetJER(float truth_jet_pt, float truth_jet_eta, int cent){
	float JER =  _h_JER[GetJetEtaBin(truth_jet_eta)][cent]->GetBinContent(_h_JER[GetJetEtaBin(truth_jet_eta)][cent]->FindBin(truth_jet_pt));
	return JER;	
}

float JetCorrector::GetJetReweightingFactor(double pt, double eta, int cent){
	//int jet_eta_bin = GetJetYBin(eta);
	int jet_eta_bin = 4;
	//cout << " jet w: " << jet_spectra_weight[7][cent]->Eval(pt) << " pt: " << pt << " cent: " << cent << endl;
	return jet_spectra_weight[jet_eta_bin][cent]->Eval(pt);
}

float JetCorrector::GetFFReweightingFactor(double z, double jet_pt, double jet_eta, int cent, bool isFine){
	//int jet_eta_bin = GetJetYBin(jet_eta);
	int jet_eta_bin = 4;
	double z_loc = z; 
	if (z>1.) z_loc = 0.95;
	if (z<0.00221) z_loc = 0.00221;
	if (jet_pt < min_jet_pt) jet_pt = min_jet_pt; //Minimum reco pT
	if (jet_pt > max_jet_pt) jet_pt = max_jet_pt; //Maximum reco pT  
	int jet_pt_bin = jet_pt_binning->FindBin(jet_pt);
	//Use rebinned version in the last two fine bins
	if (isFine && z_loc>=FF_weight[jet_eta_bin][cent][jet_pt_bin]->GetXaxis()->GetBinLowEdge(16)) isFine = false;
	int z_bin;
	if (isFine) z_bin = FF_weight_fine[jet_eta_bin][cent][jet_pt_bin]->FindBin(z_loc);
	else z_bin = FF_weight[jet_eta_bin][cent][jet_pt_bin]->FindBin(z_loc);
	//if ( cent ==5 ) cout << " FF w: " << FF_weight[7][cent][jet_pt_bin]->GetBinContent(z_bin) << " z " << z <<" jet pt: " << jet_pt << " cent: " << cent << endl;
	if (isFine) return FF_weight_fine[jet_eta_bin][cent][jet_pt_bin]->GetBinContent(z_bin);
	else return FF_weight[jet_eta_bin][cent][jet_pt_bin]->GetBinContent(z_bin);
}

float JetCorrector::GetCHPSReweightingFactor(double pt, double jet_pt, double jet_eta, int cent, bool isFine){
	//int jet_eta_bin = GetJetYBin(jet_eta);
	int jet_eta_bin = 4;
	if (jet_pt < min_jet_pt) jet_pt = min_jet_pt; //Minimum reco pT
	if (jet_pt > max_jet_pt) jet_pt = max_jet_pt; //Maximum reco pT
	int jet_pt_bin = jet_pt_binning->FindBin(jet_pt);
	int pt_bin;
	if (isFine) pt_bin = CHPS_weight_fine[jet_eta_bin][cent][jet_pt_bin]->FindBin(pt);
	else pt_bin = CHPS_weight[jet_eta_bin][cent][jet_pt_bin]->FindBin(pt);
	//cout << " CHPS w: " << CHPS_weight[7][cent][jet_pt_bin]->GetBinContent(pt_bin) << " jet pt: " << jet_pt << " cent: " << cent << endl;
	if (isFine) return CHPS_weight_fine[jet_eta_bin][cent][jet_pt_bin]->GetBinContent(pt_bin);
	else return CHPS_weight[jet_eta_bin][cent][jet_pt_bin]->GetBinContent(pt_bin);
}

float JetCorrector::GetEtaReweightingFactor(double jet_pt, double jet_eta, int cent)
{
	int pt_bin = (((TAxis*)eta_factors->Get("jetpT_binning"))->FindBin(jet_pt)) - 1;
	if (pt_bin < 0) return 1;
	int eta_bin = eta_weight[cent][pt_bin]->FindBin(jet_eta);
	double factor = eta_weight[cent][pt_bin]->GetBinContent(eta_bin);
	if (std::isnan(factor)) factor = 1.;
//	cout << Form("%i, %1.2f, %1.2f -> %i, %i -> %s, %f", cent, jet_pt, jet_eta, pt_bin, eta_bin, eta_weight[cent][pt_bin]->GetName(), factor) << endl;
	return factor;
}


