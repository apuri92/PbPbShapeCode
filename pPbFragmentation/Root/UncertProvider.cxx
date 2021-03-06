#include "pPbFragmentation/UncertProvider.h"

std::vector<double> UncertProvider::UEE_Uncert(double nominalUE, double UE_err, int size)
{
	std::vector<double> UE_sys_values;
	while (UE_sys_values.size() < size)
	{
		double ue_sys = r.Gaus(nominalUE, UE_err);
		UE_sys_values.push_back(ue_sys);
	}
	return UE_sys_values;
}

void UncertProvider::CorrectJet(const xAOD::EventInfo *eInfo, xAOD::Jet * reco, xAOD::Jet * truth = 0, int cent = 0, float FCalEt=0){
	switch (uncert_class){
		case 1:
		{	 
			 UncerJESIntrinsic(reco, eInfo);
			 break;
		}	 	 	 	  	 
		case 2:
		{
			UncerHIJESIntrinsic(reco);
			break;
		}
		case 3:
		{
		 	 UncerJERIntrinsic_Gaus(reco, truth, cent);
		 	 break;
		}
		case 6:
		{
		 	 UncerHIJESCentrality(reco,FCalEt);
		 	 break;
		} 	 	 	 	  	 
	}
	return;
}

float UncertProvider::CorrectTrackEff(float jetpt,float jety,float pt, float eta, float dR, int centrality){
	if (uncert_index == 10 || uncert_index == 11) return UncerEffMaterial(pt, eta);
	if (uncert_index == 12 || uncert_index == 13) return UncerEffFit(jetpt, jety, pt, centrality);
	if (uncert_index == 14) return UncerDense(dR);	
	else return 0;
}  

int UncertProvider::GetEtaUJERBin(float eta){
	int yBin=0;
	if (fabs(eta) < 0.3) yBin =0;
	if (fabs(eta) > 0.3 && fabs(eta) < 0.8) yBin =1;
	if (fabs(eta) > 0.8 && fabs(eta) < 1.2) yBin =2;
	if (fabs(eta) > 1.2 && fabs(eta) < 2.1) yBin =3;
	if (fabs(eta) > 2.1 && fabs(eta) < 2.8) yBin =4;
	if (fabs(eta) > 2.8 && fabs(eta) < 3.6) yBin =5;
	if (fabs(eta) > 3.6) yBin =6;

	return yBin;
}

void UncertProvider::UncerJERIntrinsic_Gaus(xAOD::Jet* recon, xAOD::Jet* truth, int cent)
{
       Float_t jetPtRecon = recon->pt();
  	   Float_t jetPtTruth = truth->pt();
   	   Float_t jetEtaRecon = recon->eta();
   	   Float_t jetEtaTruth = truth->eta();
   	   Float_t jetPhiRecon = recon->phi();
   	   Float_t jetMRecon = recon->m();
            
      Float_t uncertainty = 0;
      uncertainty = h_sJER_eta[GetEtaUJERBin(jetEtaRecon)]->Interpolate(jetPtRecon/1000.);

      //Float_t JER = JERhelper.GetJER(jetPtRecon/1000.,jetEtaRecon,cent);  centrality dependent JER
      Float_t JER =h_JER_eta[GetEtaUJERBin(jetEtaRecon)]->Interpolate(jetPtRecon/1000.);  
      Float_t smearingFactorSyst = sqrt(pow(JER+uncertainty,2)-pow(JER,2));
	  	  
	  Float_t correction = r.Gaus(0., smearingFactorSyst); 
	  //cout << "jet pt" << jetPtRecon << " truth " <<  jetPtTruth << " uncert: " << uncertainty << " JER " << JER << " correction " << correction << endl;
      Float_t jetPtRecoNew =jetPtRecon + jetPtTruth * correction;
      if (jetPtRecoNew<10) {jetPtRecoNew=1; jetMRecon=0;}
      recon->setJetP4( xAOD::JetFourMom_t(jetPtRecoNew,jetEtaRecon,jetPhiRecon,jetMRecon) );
}


//JES uncert
void UncertProvider::UncerJESIntrinsic(xAOD::Jet* recon, const xAOD::EventInfo* eInfo)
//@brief: applies JES shift based on the JES uncertainty provider tool
//@note: significance should be +/- 1
{
   
   Float_t significance=GetSysShift(uncert_index);
   Int_t component = GetJESSysComponent(uncert_index); 
   Float_t uncertainty=0;
   Float_t jetPt = recon->pt();
   Float_t jetEta = recon->eta();
   Float_t jetPhi = recon->phi();
   Float_t jetM = recon->m();

	//added eventinfo to prevent memory issues
	uncertainty = 1+ significance * (jesProv_UncerProv.getUncertainty(component, *recon, *eInfo));

   recon->setJetP4( xAOD::JetFourMom_t(jetPt*uncertainty,jetEta,jetPhi,jetM) ) ;
}


void UncertProvider::UncerHIJESIntrinsic(xAOD::Jet* recon)
//@brief: applies HI JES shift based on the JES uncertainty provider tool
//@note: significance should be +/- 1
{
   Float_t significance=GetSysShift(uncert_index);
   Int_t component = GetJESSysComponent(uncert_index);
   Float_t jetPt = recon->pt();
   Float_t jetEta = recon->eta();
   Float_t jetPhi = recon->phi();
   Float_t jetM = recon->m();
   
   Float_t uncertainty=0;
   Float_t HIJESuncertainty=0;
   if (component==1) HIJESuncertainty = sqrt(pow(HIJESProvider.GetUncertaintyComponent("flav_composition",jetPt, jetEta),2)+pow(HIJESProvider.GetUncertaintyComponent("flav_response",jetPt, jetEta),2));
   if (component==2) HIJESuncertainty = h_sJES_eta[GetEtaUJERBin(jetEta)]->Interpolate(jetPt/1000.)-1;
   //cout << "jet pt" << jetPt << " jet eta " << jetEta << " uncert " << HIJESuncertainty << endl;
   uncertainty = 1 + significance * HIJESuncertainty; 
   recon->setJetP4( xAOD::JetFourMom_t(jetPt*uncertainty,jetEta,jetPhi,jetM) );
}

void UncertProvider::UncerHIJESCentrality(xAOD::Jet* recon, float FCalEt)
//@brief: applies HI JES shift based on the JES uncertainty provider tool
//@note: significance should be +/- 1
{
   Float_t significance=GetSysShift(uncert_index);
   Float_t jetPt = recon->pt();
   Float_t jetEta = recon->eta();
   Float_t jetPhi = recon->phi();
   Float_t jetM = recon->m();
   
   Float_t uncertainty=0;
   Float_t HIJESuncertainty;
   int centrality = GetFinCentrality(FCalEt);
     
   if (centrality>60) HIJESuncertainty = 0;
   else HIJESuncertainty = (60.-(float)centrality)/60.*0.005;
   
   //cout << " FCal " << FCalEt << " cent " << centrality << " uncert " << HIJESuncertainty << endl;

   uncertainty = 1 + significance * HIJESuncertainty; 
   recon->setJetP4( xAOD::JetFourMom_t(jetPt*uncertainty,jetEta,jetPhi,jetM) );
}

int UncertProvider::GetFinCentrality(float FCalEt){
	int bin=nFineCentBins;
	
	for (int i=0;i<nFineCentBins;i++){ 
		if (fine_cent_bins[i]<FCalEt) bin = nFineCentBins - (i + 1);
		else break;
	}
	return bin;
}

void UncertProvider::UncerTrackMomentum(float& pt, float eta, float phi, int charge)
//Aply uncertainty to track momentum
{
	float charge_tmp=charge;
	if (charge_tmp>0.) charge_tmp=1;
	else if (charge_tmp<0.) charge_tmp =-1;
	
	float eta_tmp=eta;
	if(eta_tmp>2.499 )  eta_tmp= 2.499;
	if(eta_tmp<-2.499)  eta_tmp=-2.499;
	pt/=(1+charge*pt*(h_sagitta->GetBinContent(h_sagitta->FindBin(eta_tmp, phi)))*1e-3 );
}

float UncertProvider::UncerEffMaterial(float pt, float eta){
	Float_t significance=GetSysShift(uncert_index);
	if (pt>20.) pt=19.;//Maximum pt provided
	float uncertainty =0;
	for (int i=0;i<4;i++){
		int ptbin = h_eff_uncert_2015_material[i]->GetXaxis()->FindBin(pt);
		int etabin = h_eff_uncert_2015_material[i]->GetYaxis()->FindBin(eta);
		uncertainty = uncertainty + pow(h_eff_uncert_2015_material[i]->GetBinContent(ptbin,etabin),2);
	}
	//cout << "Ucnert: " << significance*sqrt(uncertainty) << endl;
	return significance*sqrt(uncertainty);			
}


float UncertProvider::UncerEffFit(float jet_pt, float jety, float pt, int centrality){
	Float_t significance=GetSysShift(uncert_index);
	if (pt < 1.29) pt = 1.29;
	if (jet_pt>500) jet_pt = 500.;
	if (pt>350.) pt=350.;//Maximum pt provided
	int ybin = jety_bins->FindBin(fabs(jety))-1;
	float uncertainty = _th1_eff_unc[ybin][centrality][jetpt_bins->FindBin(jet_pt)]->GetBinError(_th1_eff_unc[ybin][centrality][jetpt_bins->FindBin(jet_pt)]->FindBin(pt));
	//cout << "Ucnert: " << significance*uncertainty << endl;
	if (fabs(uncertainty) > 0.5) uncertainty = 0.5; // protection, the ucnertainty is <<5% 
	return significance*uncertainty;			
}


/*
float UncertProvider::UncerFit(float pt,TH1 * uncert){
	if (pt>400.) pt=400.;
	float uncertainty = 0;
	int ptbin = uncert->GetXaxis()->FindBin(pt);
	uncertainty = uncert->GetBinContent(ptbin);
	//cout << "uncert: " << uncertainty << endl;
	return uncertainty;	
}
*/
float UncertProvider::UncerDense(float dR){
	if (dR<0.1) return 0.004; //0.4% now
	else return 0.;	
}

float UncertProvider::GetMCProb(){
	if (uncert_index==4) mcprobcut = 0.4;
	else if (uncert_index==5) mcprobcut = 0.5;
	return mcprobcut;
}

//Mapping of systematic
void UncertProvider::GetTrackUncert(){         
   if (uncert_index==1) uncert_class=3; //JER uncert
   else if (uncert_index>5 && uncert_index<10) uncert_class=2; //HI JES  
   else if (uncert_index>17 && uncert_index < 42) uncert_class=1; //intrincis JES [SR is done from 18-23]
   else if (uncert_index>9 && uncert_index < 15) uncert_class=4; //tracking efficiency
   else if (uncert_index ==15) uncert_class=5; //Trk resolution
   else if (uncert_index ==16 || uncert_index ==17) uncert_class=6; //Trk resolution
   else if (uncert_index == 42) uncert_class=7; //UE map stats
   //else if (uncert_index ==18) uncert_class=7; //Trk charge scale for preliminary only
   else uncert_class = 0;
   cout << "Uncertainty class... " << uncert_class << endl; 
} 

string UncertProvider::GetSysName(int uncert){   
   string UncertLabel;
   UncertLabel="Non";
   int lastNonJES=19;
   switch(uncert){
		break;
		case 1: UncertLabel="JER_Intrinsic";
	    break;
	    case 2: UncertLabel="Tracking_Sign";
	    break;
	    case 3: UncertLabel="Unfolding";
	    break;
	    case 4: UncertLabel="Tracking_mcprob0p4";
	    break;
	    case 5: UncertLabel="Tracking_mcprob0p5";
	    break;
	    case 6: UncertLabel="JES_HI_1_P";
		break;
		case 7: UncertLabel="JES_HI_1_N";
		break;
		case 8: UncertLabel="JES_HI_2_P";
		break;
		case 9: UncertLabel="JES_HI_2_N";
		break;
		case 10: UncertLabel="Material_P";
		break;
		case 11: UncertLabel="Material_N";
		break;
		case 12: UncertLabel="eff_P";
		break;
		case 13: UncertLabel="eff_N";
		break;
		case 14: UncertLabel="effjet_N";
		break;
		case 15: UncertLabel="trk_resolution";
		break;
		case 16: UncertLabel="JES_HIC_P";
		break;
		case 17: UncertLabel="JES_HIC_N";
		break;
		case 18: UncertLabel="JES_Intrinsic_P_comp101";
		break;
		case 19: UncertLabel="JES_Intrinsic_N_comp101";
		break;
	   	case 20: UncertLabel="JES_Intrinsic_P_comp101";
		break;
	   	case 21: UncertLabel="JES_Intrinsic_N_comp101";
		break;
	   	case 22: UncertLabel="JES_Intrinsic_P_comp101";
		break;
	   	case 23: UncertLabel="JES_Intrinsic_N_comp101";
		break;
	   	case 42: UncertLabel="UE Map Statistics";
		break;


	}
	if (uncert>lastNonJES) {
		if (uncert%2==1) UncertLabel=Form("JES_Intrinsic_%i_N",(int)(uncert-lastNonJES)/2);
		else UncertLabel=Form("JES_Intrinsic_%i_P",(int)(uncert-lastNonJES+1)/2);
	}
	return UncertLabel;	
}

int UncertProvider::GetSysShift(int uncert){   
    int shift=1;
    if (uncert%2==1) shift=-1; //odds are negativni
	return shift;
}

//Return JES components
int UncertProvider::GetJESSysComponent(int uncert){   
   int Component=0;
   int lastNonJES=17;
   if (uncert<=lastNonJES){
	   switch(uncert){
			case 6: Component=1;
			break;
			case 7: Component=1;
			break;
			case 8: Component=2;
			break;
			case 9: Component=2;
			break;
		}
	}

	//TODO do it in better way 
	if (uncert>lastNonJES){
	   switch(uncert){
		   //using SR1 configuration for JES Uncertainties: all pileup and flavor is in group 100 (see config file) this corresponds to component 0, so we skip component 0 and do only components 1-3
		   case 18: Component=1;//2
			   break;
		   case 19: Component=1;//2
			   break;
		   case 20: Component=2;//3
			   break;
		   case 21: Component=2;//3
			   break;
		   case 22: Component=3;//4
			   break;
		   case 23: Component=3;//4
			   break;

			   
/* //insitu components needed for HI with JES2015_19NP.config (no pileup an b-jets)
			case 20: Component=0;//1
			break;
			case 21: Component=0;//1
			break;
			case 22: Component=1;//2
			break;
			case 23: Component=1;//2
			break;
			case 24: Component=8;//9
			break;
			case 25: Component=8;//9
			break;
			case 26: Component=9;//10
			break;
			case 27: Component=9;//10
			break;
			case 28: Component=10;//11
			break;
			case 29: Component=10;//11
			break;
			case 30: Component=11;//12
			break;
			case 31: Component=11;//12
			break;
			case 32: Component=12;//13
			break;
			case 33: Component=12;//13
			break;
			case 34: Component=13;//14
			break;
			case 35: Component=13;//14
			break;
			case 36: Component=14;//15
			break;
			case 37: Component=14;//15
			break;
			case 38: Component=15;//16
			break;
			case 39: Component=15;//16
			break;
			case 40: Component=16;//17
			break;
			case 41: Component=16;//17
			break;
 */
		}
    }
   return Component;
}
