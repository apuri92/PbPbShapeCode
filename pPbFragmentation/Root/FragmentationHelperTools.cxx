#include "pPbFragmentation/FragmentationHelperTools.h"

#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <TMath.h>
using namespace std;


double MTCorrector::GetZ(double track_pt, double track_eta, double track_phi, double jet_pt, double jet_eta, double jet_phi, bool UseAltzDef){
	
	double z;
	if (!UseAltzDef) {
		float R = JetHelperTools::DeltaR(track_phi,track_eta,jet_phi,jet_eta);
		z = cos(R)*track_pt / jet_pt;	
	}	
	else {
		TVector3 jet3vector; jet3vector.SetPtEtaPhi(jet_pt,jet_eta,jet_phi);
		TVector3 particle3vector; particle3vector.SetPtEtaPhi(track_pt,track_eta,track_phi);
		z =  (jet3vector * particle3vector)/ pow(jet_pt*cosh(jet_eta),2);	
	}
	return z;

}

std::vector<float> MTCorrector::GetIsolation(std::vector<float>& jetEt, std::vector<float>& jetEta, std::vector<float>& jetPhi, int jetRadius)
//@brief: remove jets that are not isolated according to a given scheme
{

  Float_t deltaRMin_recon = 0.;
  //Float_t ptMin_recon = 0.;
  Int_t  njets_recon = 0.;

  deltaRMin_recon = 1.0; njets_recon = jetEt.size();
 /* if (jetRadius==4) ptMin_recon = 100./2.;
  if (jetRadius==3) ptMin_recon = 92./2.;
  if (jetRadius==2) ptMin_recon = 85./2.; */
  std::vector<float> Isolated;
  Isolated.clear();
  for (int i=0; i<njets_recon; i++){
     //bool IsIsolated = true;
     float maxpT = -9999.;
     //if (jetEt.at(i) <= ptMin_recon) continue;
     for (int j=0; j<njets_recon; j++){
        if (i==j) continue;
        //if (jetEt.at(j) < ptMin_recon) continue;
        float R = JetHelperTools::DeltaR(jetPhi.at(i),jetEta.at(i),jetPhi.at(j),jetEta.at(j));
        if  (R > deltaRMin_recon) continue;  
        //IsIsolated = false;
        if (jetEt.at(j) < maxpT)  continue;
        maxpT = jetEt.at(j);
     }
     /*
     if (IsIsolated) {
     	Isolated.push_back(1);
     }
     else {
     	Isolated.push_back(0);
     }
     */
     Isolated.push_back(maxpT);
  }
  return Isolated;
}


std::vector<bool> MTCorrector::GetIsolation(std::vector<float>& jetpT, std::vector<float>& jetEta, std::vector<float>& jetPhi, int iso_R, double iso_pT)
{
	bool isIso = true;
	std::vector<bool> Isolated;
	double pT_cut = iso_pT;

	for (int i = 0; i<jetpT.size(); i++)
	{
		if (iso_pT < 1 ) pT_cut = jetpT.at(i)*iso_pT; // if in range 0-1: isolating pT requirement defined in % of jet pT
		if (iso_pT < 0 ) pT_cut = jetpT.at(i); // -1 is for isolating pT requirement is 100% of jet pT

		for (int j = 0; j < jetpT.size(); j++)
		{
			if (j == i) continue;

			double R = JetHelperTools::DeltaR(jetPhi.at(i),jetEta.at(i),jetPhi.at(j),jetEta.at(j));
			if (R < iso_R && jetpT.at(j) > pT_cut) isIso = false;
		}

		Isolated.push_back(isIso);
	}

	return Isolated;

}

bool MTCorrector::GetFJR(xAOD::Jet* jet, InDet::InDetTrackSelectionTool *m_track_selection_tool){
	
	const xAOD::JetFourMom_t jet_4mom = jet->jetP4();
	float jetEta  = jet_4mom.eta();
	
	const xAOD::Vertex* jetorigin=nullptr;
	bool hasVertex=jet->getAssociatedObject<xAOD::Vertex>("OriginVertex", jetorigin); 
	  
	std::vector<const xAOD::TrackParticle*> tracks;
	bool havetracks = jet->getAssociatedObjects("GhostTrack", tracks);
	if ( ! havetracks ) std::cerr << "no tracks found for given jet" << std::endl;
	float sumpt = 0.; //this is the sum pT track variable to cut on for the FJR
	   
	bool isFake = true;
	//Jet is outside the tracker -> no FJR
	if (fabs(jetEta) > 2.5  ) {
		return false;
	}
	for(unsigned int itrack = 0; itrack < tracks.size(); itrack++){
		if (!tracks[itrack]) continue;
		const xAOD::TrackParticle* t = tracks[itrack];
		bool removetrk = false;
		if(hasVertex) {
		   if(!m_track_selection_tool->accept(*t,jetorigin)) removetrk = true;
		 }
		else {
		   if(!m_track_selection_tool->accept(*t)) removetrk = true;
		 }
		TLorentzVector p_trk;
		if (t->pt()*1e-3 < 4.) continue; //require the tracks to be greater than 4 GeV
		if (removetrk == true) continue; //only use jets that pass the track selector
		p_trk.SetPtEtaPhiM(t->pt()*1e-3, t->eta(), t->phi(), 0.14); //mass pion in GeV 
		sumpt+=(t->pt()*1e-3);
	}
	if (sumpt > 8.) isFake = false;
	return isFake;
	
}

bool MTCorrector::GetFJR(float jetEta, float jetPhi, const xAOD::JetContainer* track_jets)
{

	double single_track_pT_cut = 8.;
	double di_track_pT_cut = 4.;
	int single_track_count_cut = 1;
	int di_track_count_cut = 2;
	double dRcut = 0.2;
	
	bool isFake = true;
	//Jet is outside the tracker -> no FJR
	if (fabs(jetEta) > 2.5  ) {
		return false;
	}	
	//Jet is inside the tracker: 
	int single_track_counts=0;
	int di_track_counts=0;
	xAOD::JetContainer::const_iterator track_jet_itr = track_jets->begin();
	xAOD::JetContainer::const_iterator track_jet_end = track_jets->end();
	for( ; track_jet_itr != track_jet_end; ++track_jet_itr )
	{
		xAOD::JetFourMom_t track_jet_4mom = (*track_jet_itr)->jetP4();
		double pt    = (track_jet_4mom.pt() * 0.001 );
		double eta    = (track_jet_4mom.eta());
		double phi    = (track_jet_4mom.phi());
		double R = JetHelperTools::DeltaR(jetPhi,jetEta,phi,eta);
		if (R > dRcut) continue;
		if (pt > single_track_pT_cut) single_track_counts++;
		if (pt > di_track_pT_cut) di_track_counts++;			
	}
	if (single_track_counts>=single_track_count_cut || di_track_counts>=di_track_count_cut) isFake = false;
	
	return isFake;

}

std::vector<int> MTCorrector::TruthMatching(std::vector<float>& reco_jetpT, std::vector<float>& reco_jetEta, std::vector<float>& reco_jetPhi,
											 std::vector<float>& truth_jetpT, std::vector<float>& truth_jetEta, std::vector<float>& truth_jetPhi,
											 float match_R)
{
	std::vector<int> TruthMatchedIndeces;

//	cout << "Start" << endl << endl;
	for(unsigned int i=0; i<reco_jetpT.size(); i++)
	{
		double jet_pt = reco_jetpT.at(i) / 1.0;
		double jet_eta = reco_jetEta.at(i);
		double jet_phi = reco_jetPhi.at(i);

		int TruthJetIndex=-1;
		float dR_min=999;

		for(unsigned int j=0; j<truth_jetpT.size(); j++)
		{

			double truth_pt = truth_jetpT.at(j);
			double truth_eta = truth_jetEta.at(j);
			double truth_phi = truth_jetPhi.at(j);

			float R = JetHelperTools::DeltaR(truth_phi,truth_eta,jet_phi,jet_eta);

			//Filter out outliers
			//if (fabs(jet_pt-truth_pt)/truth_pt > 4) continue;

			//dR cut
			if(R < match_R && R < dR_min)
			{
				dR_min = R;
				TruthJetIndex = j;
			}
		}

		TruthMatchedIndeces.push_back(TruthJetIndex);

	}

	return TruthMatchedIndeces;
	
}


float MTCorrector::R2Matching(float &reco_jetEta, float &reco_jetPhi,
							 std::vector<float>& a2_jetEta, std::vector<float>& a2_jetPhi,
							 int match_R)
{
	double jet_a2_eta = -999, jet_a2_phi = -999;
	float dR_min=999;
	bool matched_a2 = false;
	for(unsigned int j=0; j<a2_jetEta.size(); j++)
	{
		float R = JetHelperTools::DeltaR(a2_jetPhi.at(j),a2_jetEta.at(j),reco_jetPhi,reco_jetEta);

		if(R < 0.3 && R<dR_min)
		{
			matched_a2 = true;
			jet_a2_eta=a2_jetEta.at(j);
			jet_a2_phi=a2_jetPhi.at(j);
			dR_min = R;
		}
	}
	if (matched_a2)
	{
		reco_jetEta = jet_a2_eta;
		reco_jetPhi = jet_a2_phi;
	}
	return dR_min;
}


void MTCorrector::SetupBinning(Int_t scheme, string variable, Double_t array[1000], Int_t &num)
//@brief: to setup binning for histograms
{
   Float_t a, c;
   Int_t k;

   if ((scheme==0) && (variable=='z'))
     {num = 26;		
      array[0] = 0.001;// to overflow
      //c = 0.01;
      c = 0.002511886; 	
      a = pow( (1./c), 1./num );
      k = 1; //1 to overflow
      printf("... z-binning : ");
      for (int i=0; i<=num+1; i++)
        {Float_t value = c * pow( a, i );
         array[k] = value;
         printf("%.4f, ", array[k]);
         k++;
        }
      num = k-1; 
      cout << " ... " << num << endl;
     }
	/*
   if ((scheme==0) && (variable=="pt-trk"))
     {num = 35; 
      array[0] = 0;
      //c = 0.634242388;
      c = 0.1;
      a = 4./3.;
      k = 0;
      printf("\n... pt-trk-binning : ");
      for (int i=0; i<=num+1; i++)
        {Float_t value = c * pow( a, i );
         //if (value<2) continue;
         if (value>500) continue;
         array[k] = value;
         printf("%.4f, ", array[k]);
         k++;
        }
      num = k-1;
      cout << " ... \n" << num << endl;
     }
	*/
	if ((scheme==0) && (variable=="pt-trk"))
	{ num = 29;
		// old double bins[18] = {0.50, 0.63, 1.00, 1.58, 2.51, 3.98, 6.31, 10.00, 15.85, 25.12, 39.81, 63.09, 100.00, 125.89, 158.49, 199.53, 314.98, 419.98};
		double bins[30] = {0.1000, 0.1333, 0.1778, 0.2370, 0.3160, 0.4214, 0.5619, 0.9, 1.0, 1.3318, 1.7758, 2.3677, 3.1569, 4.2092, 5.6123, 7.4831, 9.9775, 13.3033, 17.7377, 23.6503, 31.5337, 42.0449, 56.0599, 74.7466, 99.6621, 132.8828, 177.1771, 236.2361, 314.9815, 419.9753};
		printf("\n... pt-trk-binning rebbin : ");
		for (int i=0; i<=num; i++)
		{
			array[i]=bins[i];
			printf("%.2f, ", array[i]);
		}
		cout << " ... " << num << endl;
	}
	
	if ((scheme==0) && (variable=="pt-trk-sum"))
	{ num = 40;
		// old double bins[18] = {0.50, 0.63, 1.00, 1.58, 2.51, 3.98, 6.31, 10.00, 15.85, 25.12, 39.81, 63.09, 100.00, 125.89, 158.49, 199.53, 314.98, 419.98};
		double bins[41] = {0.1000, 0.1333, 0.1778, 0.2370, 0.3160, 0.4214, 0.5619, 0.9, 1.0, 1.3318, 1.7758, 2.3677, 3.1569, 4.2092, 5. ,5.6123, 6.2 ,7.4831, 8.3 ,9.9775, 11. ,13.3033, 14.5 , 17.7377, 20. , 23.6503, 26. ,31.5337, 35. ,42.0449, 46.0, 56.0599, 64., 74.7466, 85. , 99.6621, 132.8828, 177.1771, 236.2361, 314.9815, 419.9753};
		printf("\n... pt-trk-binning rebbin : ");
		for (int i=0; i<=num; i++)
		{
			array[i]=bins[i];
			printf("%.2f, ", array[i]);
		}
		cout << " ... " << num << endl;
	}

	if ((scheme==0) && (variable=="pt-jet-rebin"))
	{num = 14;
		//def
		array[0] = 20; array[1] = 32;array[2] = 45;array[3] = 60;array[4] = 80;array[5] = 110;array[6] = 160;array[7] = 210;array[8] = 260;array[9] = 310;array[10] = 400;array[11] = 500; array[12] = 600; array[13] = 800; array[14] = 1000;
		//array[0] = 32; array[1] = 40;array[2] = 45;array[3] = 60;array[4] = 80;array[5] = 110;array[6] = 160;array[7] = 210;array[8] = 260;array[9] = 310;array[10] = 400;array[11] = 500; array[12] = 600; array[13] = 800; array[14] = 1000;
		printf("... pt-jet-binning : ");
		//num = 100;
		Float_t value=0;
		for (int i=0; i<=num; i++)
		{
			//array[i] = value;
			//value = value + 2;
			printf("%.0f, ", array[i]);
		}
		cout << " ... " << num << endl;
	}
	
	if ((scheme==0) && (variable=="pt-jet-PbPb"))
	{num = 16;
		//def
		//double bins[13]={63.096, 80.0, 100.000, 125.892,  158.488,  199.525,  251.186,  316.224,  398.101,  501.178,  630.944, 794.308, 999.970};
		//double bins[17]={25.119, 31.623, 40.0, 50.119, 63.096, 79.433, 100.000, 125.892,  158.488,  199.525,  251.186,  316.224,  398.101,  500.,  630.944, 794.308, 999.970};
		double bins[17]={25.119, 31.623, 40.0, 50.119, 63.096, 82., 100.000, 125.892,  158.488,  199.525,  251.186,  316.224,  398.101,  500.,  630.944, 794.308, 999.970};
		
		
		// Binning changed from 501.178 to 500 GeV
		//double bins[12]={80.0, 100.000, 125.892,  158.488,  199.525,  251.186,  316.224,  398.101,  501.178,  630.944, 794.308, 999.970};
		printf("... pt-jet-binning : ");
		//num = 100;
		Float_t value=0;
		for (int i=0; i<=num; i++)
		{
			array[i]=bins[i];
			//value = value + 2;
			printf("%.0f, ", array[i]);
		}
		cout << " ... " << num << endl;
	}
	
	if ((scheme==0) && (variable=="m-jet"))
	{num = 12;
		double bins[13] = { -50, -30,-10,0,10,20,30,40,50,75,100,200,300};
		printf("... m-jet-binning : ");
		//num = 100;
		Float_t value=0;
		for (int i=0; i<=num; i++)
		{
			array[i]=bins[i];
			//value = value + 2;
			printf("%.0f, ", array[i]);
		}
		cout << " ... " << num << endl;
	}

	if ((scheme==0) && (variable=="z-rebin"))
	{ num = 19;
		//double bins[17] = {0.0050, 0.0063, 0.0079, 0.0100, 0.0158, 0.0251, 0.0398, 0.0631, 0.1000, 0.1585, 0.2512, 0.3981, 0.6310, 1.0000, 1.2589, 1.5849, 1.9953};
		double bins[20] =      {0.0010, 0.0015, 0.0022, 0.0032, 0.0046, 0.0068, 0.0100, 0.0158, 0.0251, 0.0398, 0.0631, 0.1000, 0.1585, 0.2512, 0.3981, 0.6310, 1.0000, 1.2589, 1.5849, 1.9953};
		printf("\n... z-binning rebbin: ");
		for (int i=0; i<=num; i++)
		{
			array[i]=bins[i];
			printf("%.4f, ", array[i]);
		}
		cout << " ... " << num << endl;
	}

	if ((scheme==0) && (variable=="pt-trk-rebin"))
	{ num = 16;
		// old double bins[18] = {0.50, 0.63, 1.00, 1.58, 2.51, 3.98, 6.31, 10.00, 15.85, 25.12, 39.81, 63.09, 100.00, 125.89, 158.49, 199.53, 314.98, 419.98};
		double bins[17] = {0.50, 0.9, 1.00, 1.58, 2.51, 3.98, 6.31, 10.00, 15.85, 25.12, 39.81, 63.09, 100.00, 158.49, 251.18, 398.11, 630.96};
		printf("\n... pt-trk-binning rebbin : ");
		for (int i=0; i<=num; i++)
		{
			array[i]=bins[i];
			printf("%.2f, ", array[i]);
		}
		cout << " ... " << num << endl;
	}

	if ((scheme==0) && (variable=="pt-trk-shape"))
	{
//		num = 9;
//		double bins[10] = {0.50, 0.63, 1.00, 1.58, 2.51, 3.98, 6.31, 15.85, 39.81, 100.00};

		num = 10;
		double bins[11] = {0.50, 0.9, 1.00, 1.58, 2.51, 3.98, 6.31, 10.00, 25.12, 63.09, 158.49};

		printf("\n... pt-trk-binning for shape : ");
		for (int i=0; i<=num; i++)
		{
			array[i]=bins[i];
			printf("%.2f, ", array[i]);
		}
		cout << " ... " << num << endl;
	}

   if ((scheme==0) && (variable=="pt-jet"))
     {/*num = 100;
      array[0] = 0;
      k = 0;
      printf("\n... pt-jet-binning : ");
      for (int i=0; i<=num; i++)
        {Float_t value = i*5;
         array[k] = value;
         printf("%.4f, ", array[k]);
         k++;
        }
      num = k-1;
      cout << " ... \n" << num << endl;*/
      num = 14;
      array[0] = 10; array[1] = 32;array[2] = 45;array[3] = 60;array[4] = 80;array[5] = 110;array[6] = 160;array[7] = 210;array[8] = 260;array[9] = 310;array[10] = 400;array[11] = 500;array[12] = 600;array[13] = 800;array[14] = 1000;
      printf("\n... pt-jet-binning : ");
      for (int i=0; i<=num; i++)
        {
         printf("%.0f, ", array[i]);
        }
       cout << " ... \n" << num << endl;
     } 

   if ((scheme==0) && (variable=="eta-jet"))
     {
      printf("\n... eta-jet-binning : four simple bins ");
      array[0] = 0;
      array[1] = 0.3;
      array[2] = 0.8;
      array[3] = 1.2;
      array[4] = 2.1;
      num = 4;
     }
    if ((scheme==0) && (variable=="eta-trk"))
     {
      printf("\n... eta-trk-binning : four simple bins ");
      num = 50;
      double value=-2.500;
      for (int i=0; i<=num; i++)
        {
         array[i] = value;
         printf("%f, ", array[i]);
         value = value + 0.100;
        }
     }
     if ((scheme==0) && (variable=="phi-trk"))
     {
      printf("\n... phi-trk-binning : four simple bins ");
      num = 64;
      Float_t value=-TMath::Pi(); 
      for (int i=0; i<=num; i++)
        {
         array[i] = value;
         printf("%.4f, ", array[i]);
         value = value + 0.1;
        }
     }
     if ((scheme==0) && (variable=="eta-trk-coars"))
     {
      printf("\n... eta-trk-binning : four simple bins ");
      num = 21;
      double value=-2.1;
      for (int i=0; i<=num; i++)
        {
         array[i] = value;
         printf("%f, ", array[i]);
         value = value + 0.200;
        }
     }
     if ((scheme==0) && (variable=="phi-trk-coars"))
     {
      printf("\n... phi-trk-binning : four simple bins ");
      num = 32;
      Float_t value=-TMath::Pi(); 
      for (int i=0; i<=num; i++)
        {
         array[i] = value;
         printf("%.4f, ", array[i]);
         value = value + TMath::Pi()/16.;
        }
     }
     if ((scheme==0) && (variable=="eta-trk-fine"))
     {
      printf("\n... eta-trk-binning : four simple bins ");
      num = 100;
      double value=-2.500;
      for (int i=0; i<=num; i++)
        {
         array[i] = value;
         printf("%f, ", array[i]);
         value = value + 0.05;
        }
     }
     if ((scheme==0) && (variable=="phi-trk-fine"))
     {
      printf("\n... phi-trk-binning : four simple bins ");
      num = 128;
      Float_t value=-TMath::Pi(); 
      for (int i=0; i<=num; i++)
        {
         array[i] = value;
         printf("%.4f, ", array[i]);
         value = value + TMath::Pi()*2/num;
        }
     }
     if ((scheme==0) && (variable=="density"))
     {
      printf("\n... density-binning : four simple bins ");
      num = 50;
      Float_t value=0; 
      for (int i=0; i<=num; i++)
        {
         array[i] = value;
         printf("%.4f, ", array[i]);
         value = value + 1;
        }
     }
     if ((scheme==0) && (variable=="hits"))
     {
      printf("\n... hits-binning : four simple bins ");
      num = 100;
      Float_t value=0; 
      for (int i=0; i<=num; i++)
        {
         array[i] = value;
         printf("%.4f, ", array[i]);
         value = value + 5;
        }
     }
     if ((scheme==0) && (variable=="hits_fine"))
     {
      printf("\n... fine-hits-binning : four simple bins ");
      num = 20;
      Float_t value=0; 
      for (int i=0; i<=num; i++)
        {
         array[i] = value;
         printf("%.4f, ", array[i]);
         value = value + 1;
        }
     }
     if ((scheme==0) && (variable=="d0z0"))
     {
      printf("\n... d0zo0: four simple bins ");
      num = 300;
      Float_t value=-1.5; 
      for (int i=0; i<=num; i++)
        {
         array[i] = value;
         printf("%.4f, ", array[i]);
         value = value + 0.01;
        }
     }
     if ((scheme==0) && (variable=="dR"))
     {
      printf("\n... dR : four simple bins ");
      num = 15;
      //Float_t value=0; 
	  array[0]=0.;array[1]=0.001;array[2]=0.002;array[3]=0.004;array[4]=0.006;array[5]=0.008;array[6]=0.01;array[7]=0.015;array[8]=0.02;array[9]=0.03;array[10]=0.04;array[11]=0.05;array[12]=0.07;array[13]=0.1;array[14]=0.2;;array[15]=0.4;
     }
     
     if ((scheme==0) && (variable=="dR-c"))
     {
      printf("\n... density-binning : four simple bins ");
      num = 9;
      //Float_t value=0; 
	  array[0]=0.;array[1]=0.01;array[2]=0.02;array[3]=0.05;array[4]=0.1;array[5]=0.5;array[6]=1;array[7]=2;array[8]=5;array[9]=10;
     }
     
     if ((scheme==0) && (variable=="dR-RdR"))
     {
      printf("\n... density-binning : four simple bins ");
      num = 14;
      //Float_t value=0; 
	  array[0]=0.;array[1]=0.2;array[2]=0.4;array[3]=0.5;array[4]=0.6;array[5]=0.7;array[6]=0.8;array[7]=0.9;array[8]=1.0;array[9]=1.1;array[10]=1.2;array[11]=1.3;array[12]=1.4;array[13]=1.5;array[14]=1.6;
     }
     
     if ((scheme==0) && (variable=="MC-prob"))
     {
      printf("\n... MC probability binning ");
      num = 10;
      //Float_t value=0; 
	  array[0]=0.0;array[1]=0.1;array[2]=0.2;array[3]=0.3;array[4]=0.4;array[5]=0.5;array[6]=0.6;array[7]=0.7;array[8]=0.8;array[9]=0.9;array[10]=1.0;
      printf("\n... MC-prob-binning : ");
      for (int i=0; i<=num; i++)
        {
         printf("%.0f, ", array[i]);
        }
       cout << " ... \n" << num << endl;	
     }
     if ((scheme==0) && (variable=="dpT_fine"))
     {
      printf("\n... delta pT bins");
      num = 21;
      //Float_t value=0; 
	  array[0]=-1.;array[1]=-0.1;array[2]=-0.05;array[3]=-0.03;array[4]=-0.02;array[5]=-0.01;array[6]=-0.005;array[7]=-0.001;array[8]=-0.0005;array[9]=-0.0002;array[10]=-0.0001;array[11]=0.0;array[12]=0.0001;array[13]=0.0002;array[14]=0.0005;array[15]=0.001;array[16]=0.005;array[17]=0.01;array[18]=0.02;array[19]=0.03;array[20]=0.05;array[21]=0.1;array[21]=1;
      printf("\n... delta pt-binning : ");
      for (int i=0; i<=num; i++)
        {
         printf("%.0f, ", array[i]);
        }
       cout << " ... \n" << num << endl;	
     }
     if ((scheme==0) && (variable=="dpT"))
     {
      printf("\n... delta pT bins ");
      num = 14;
      //Float_t value=0; 
	  array[0]=-100;array[1]=-10.;array[2]=-1.;array[3]=-0.1;array[4]=-0.01;array[5]=-0.005;array[6]=-0.001;array[7]=0.0;array[8]=0.001;array[9]=0.005;array[10]=0.01;array[11]=0.1;array[12]=1;array[13]=10;array[14]=100;
      printf("\n... delta pt-binning : ");
      for (int i=0; i<=num; i++)
        {
         printf("%.0f, ", array[i]);
        }
       cout << " ... \n" << num << endl;	
     }
     if ((scheme==0) && (variable=="trk_res"))
     {
      printf("\n... trk-res-binning : four simple bins ");
      num = 400;
      Float_t value=-1; 
      for (int i=0; i<=num; i++)
        {
         array[i] = value;
         printf("%.4f, ", array[i]);
         value = value + 0.005;
        }
     }

	if ((scheme==0) && (variable=="resp"))
	{
		printf("\n... resp-binning : ");
		num = 80;
		array[0] = -0.4;
//		array[0] = 0.499375;
		for (int i=1; i<=num; i++)
		{
			array[i] = array[0]+0.01*i;
//			array[i] = array[0]+i/800.;
			printf("%4.4f, ", array[i]);
		}
		cout << " ... " << num << endl;


	}

	if ((scheme==0) && (variable=="resp-fine"))
	{
		printf("\n... fine-resp-binning : ");
		num = 800;
		array[0] = -0.020-(0.06/2400);
		//		array[0] = 0.499375;
		for (int i=1; i<=num; i++)
		{
			array[i] = array[0]+(0.06/1200)*i;
			printf("%4.4f, ", array[i]);
		}
		cout << " ... " << num << endl;

	}

	if ((scheme==0) && (variable=="dr_fine"))
	{
		printf("\n... R position res-binning : ");
		num = 240;
		array[0] = 0.;
		array[num] = 1.2;
		double w = (array[num] - array[0])/num;
		//		array[0] = 0.499375;
		for (int i=1; i<=num; i++)
		{
			array[i] = array[0]+(w*i);
			printf("%4.4f, ", array[i]);
		}
		cout << " ... " << num << endl;

	}


	if ((scheme==0) && (variable=="eta-fine"))
	{
		printf("\n... fine eta binning : ");
		num = 50;
		array[0] = -2.5;
		for (int i=1; i<=num; i++)
		{
			array[i] = array[0]+i*0.1;
			printf("%4.4f, ", array[i]);
		}
		cout << " ... " << num << endl;

	}
	
	if ((scheme==0) && (variable=="PbPb_runs"))
	{

		num = 32;
		double bins[33] = {286711, 286717, 286748, 286767, 286834, 286854, 286908, 286990, 287038, 287044, 287068, 287222, 287224, 287259, 287270, 287281, 287321, 287330, 287334, 287378, 287380, 287382, 287560, 287594, 287632, 287706, 287728, 287827, 287843, 287866, 287924, 287931, 287950};

		printf("\n... Runs  : ");
		for (int i=0; i<=num; i++)
		{
			array[i]=bins[i];
			printf("%.2f, ", array[i]);
		}
		cout << " ... " << num << endl;
	}


	if ((scheme==0) && (variable=="PsiBins"))
	{
		num = 16;
		int i = 0;
		printf("\n... Psi bins  : ");
		for (float dPsi=0.;dPsi<=TMath::Pi()+0.2;dPsi=dPsi+0.2){
			array[i] = dPsi;
			printf("%.2f, ", array[i]);
			i++;
		}

		cout << " ... " << num << endl;


	}

}

int MTCorrector::GetTrkpTBin(float pt){
	if (pt>100.) return 9;
	if (pt>70.) return 8;
	if (pt>40.) return 7;
	if (pt>20.) return 6;
	if (pt>10.) return 5;
	if (pt>4.) return 4;
	if (pt>2.) return 3;
	if (pt>1.) return 2;
	if (pt>0.5) return 1;
	if (pt>0.2) return 0;
	
	return 0;
}

int MTCorrector::GetRunNumberBin(int RunNumber){
	int num = 32;
	int bins[32] = {286711, 286717, 286748, 286767, 286834, 286854, 286908, 286990, 287038, 287044, 287068, 287222, 287224, 287259, 287270, 287281, 287321, 287330, 287334, 287378, 287380, 287382, 287560, 287594, 287632, 287706, 287728, 287827, 287843, 287866, 287924, 287931};
	for (int i=0;i<num;i++){
		if (bins[i]==RunNumber) return i;
	}
}
