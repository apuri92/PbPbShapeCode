#include "pPbFragmentation/GlobalHelper.h"
#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include <TMath.h>
//#include <xAODCaloEvent/​CaloClusterContainer.h>
using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// there are following options currently implemented :
//
//  0 = do not distinguish centrality => out->iNum = 1 
//  2 = Pb+Pb centrality              => out->iNum = 7, 0-10%=b, 10-20%=c, 20-30%=m, 30-40%=n, 40-50%=o, 50-60%=p, 60-80%=j
//  20 = p+Pb centrality              => out->iNum = 7, 0-10%=b, 10-20%=c, 20-30%=m, 30-40%=n, 40-50%=o, 50-60%=p, 60-100%=j
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


float GetEventPlane(const xAOD::CaloClusterContainer *hiclus)
// @brief: returns a event plane angle
{
	float totalEt_P=0;
	float totalEt_N=0;
	float m_Qnx_P;
	float m_Qny_P;
	float m_Qnx_N;
	float m_Qny_N;
	for(auto cl: *hiclus){
      double eta = cl->eta0();
      double phi = cl->phi0();
      if(fabs(eta) > 3.2 and fabs(eta) < 4.8){
        double et = cl->altE()/cosh(eta);
   		if(eta > 0.) totalEt_P +=et;
        if(eta < 0.) totalEt_N +=et;
        if(eta > 0.){
            m_Qnx_P +=(et*TMath::Cos(2*phi));
            m_Qny_P +=(et*TMath::Sin(2*phi));
        }
        if(eta < 0.){
            m_Qnx_N +=(et*TMath::Cos(2*phi));
            m_Qny_N +=(et*TMath::Sin(2*phi));
        }
     }
	}

	float psiEP_N = std::atan2(m_Qny_N/totalEt_N, m_Qnx_N/totalEt_N);
	float psiEP_P =  std::atan2(m_Qny_P/totalEt_P, m_Qnx_P/totalEt_P);
	return GetAveragePsi(psiEP_N, psiEP_P);
}

float GetEventPlane(const xAOD::HIEventShapeContainer* calos){
	float psi_2;
	//for(const xAOD::HIEventShape calo_itr : calos)
		//{
		    //std::string summary;
		    //if(calos->isAvailable<std::string>("Summary")) summary=calos->auxdata<std::string>("Summary");
		    //if(summary.compare("FCal")==0)
		    {
		        float FCal_Et=(calos->at(5)->et()*0.001*0.001);		        
		        double qx_2=(calos->at(5)->etCos().at(1));
		        double qy_2=(calos->at(5)->etSin().at(1));
		        float N_psi_2 = std::atan2(qy_2,qx_2);
		        psi_2 = N_psi_2/2.0;        
		        //break;
		    }
		//}
	return psi_2;	
}

Float_t GetAveragePsi(Float_t psi1, Float_t psi2)
//@brief: calculate the average event plane
{
   Float_t phase = (fabs(psi1-psi2)<(TMath::Pi()/2.))? 0:(TMath::Pi()/2.);
   return ( (psi1+psi2)/2. + phase );
}

Int_t GetPsiBin(Float_t psi)
{
   int bin = 0;
   for (float dPsi=0.1;dPsi<=TMath::Pi();dPsi=dPsi+0.1){
   		//cout << "psi " << psi << " dPsi " << dPsi << endl;
   		if (psi<=dPsi) return bin;
   		bin++;
   }
   return bin;
}

int GetCentralityBin(Int_t centralityScheme, float FCal_Et, bool isMC)
// @brief: returns a centrality bin [0-9] based on centralityScheme and MTGlobalEvent
{
	return GetGlobalBin(centralityScheme,FCal_Et, isMC);
}

void SetRejectionHistogram(TH1D* h)
{
	h->GetXaxis()->SetBinLabel(1,"All"); //0.5
	h->GetXaxis()->SetBinLabel(2,"centrality"); // 1.5
	h->GetXaxis()->SetBinLabel(3,"GRL"); // 2.5
	h->GetXaxis()->SetBinLabel(4,"Vertex"); //3.5
	h->GetXaxis()->SetBinLabel(5,"LAr Quality"); //4.5
	// h->GetXaxis()->SetBinLabel(4,"Timing");
	// 
	h->GetXaxis()->SetBinLabel(6,"vx_z"); //5.5
	h->GetXaxis()->SetBinLabel(7,"Pileup"); //6.5
	h->GetXaxis()->SetBinLabel(8,"Accepted");
	h->GetXaxis()->SetBinLabel(9,"Passed trigger");
}

int GetGlobalBin(Int_t centralityScheme, float FCal_Et, bool isMC)
// @brief: returns a global bin [0-9] based on centralityScheme, MTGlobalEvent, and MTEvent (e.g. position of leading jet wrt RP)
{
	Float_t centrality = FCal_Et;

	if (centralityScheme==1)
	{
		return 0;
	}
	else if (centralityScheme==30) // Pb+Pb 2015
	{
		// nominal 85%, full Fcal
		if ( 2.98931 	<=centrality && centrality< 6.00  ) return 0;		// 0-10%
		if ( 2.04651	<=centrality && centrality< 2.98931  ) return 1;	// 10-20%
		if ( 1.36875	<=centrality && centrality< 2.04651  ) return 2;	// 20-30%
		if ( 0.87541	<=centrality && centrality< 1.36875  ) return 3;	// 30-40%
		if ( 0.525092	<=centrality && centrality< 0.87541  ) return 4;	// 40-50%
		if ( 0.289595	<=centrality && centrality< 0.525092 ) return 5;	// 50-60%
		if ( 0.063719	<=centrality && centrality< 0.289595 ) return 6;	// 60-80%
		
		// Hijing doesn't have the same FCal distribution as data, let's keep everything
		if(isMC && centrality>= 6.0) return 0;
		if(isMC && 0.0<=centrality && centrality<0.063719) return 6;
		
		// you cannot do this -- if you want 60-70% or 70-80% you need to setup a different centrality scheme
		//if ( 0.144140 	<=centrality && centrality< 0.289595 ) return 7;		// 60-70%
		//if ( 0.063719 	<=centrality && centrality< 0.144140 ) return 8; 		// 70-80%

		return -1;
	}
	else if (centralityScheme==31) // Pb+Pb 2015, merged 40-60%
	{
		// nominal 85%, full Fcal
		if ( 2.98931 	<=centrality && centrality< 6.00  ) return 0;		// 0-10%
		if ( 2.04651	<=centrality && centrality< 2.98931  ) return 1;	// 10-20%
		if ( 1.36875	<=centrality && centrality< 2.04651  ) return 2;	// 20-30%
		if ( 0.87541	<=centrality && centrality< 1.36875  ) return 3;	// 30-40%
		if ( 0.525092	<=centrality && centrality< 0.87541  ) return 4;	// 40-50%
		if ( 0.289595	<=centrality && centrality< 0.525092 ) return 4;	// 50-60%
		if ( 0.063719	<=centrality && centrality< 0.289595 ) return 5;	// 60-80%
		
		// Hijing doesn't have the same FCal distribution as data, let's keep everything
		if(isMC && centrality>= 6.0) return 0;
		if(isMC && 0.0<=centrality && centrality<0.063719) return 5;
		
		return -1;
	}
	else if (centralityScheme==32) // Pb+Pb 2015, fine peripheral bins
	{
		// nominal 85%, full Fcal
		if ( 2.98931 	<=centrality && centrality< 6.00  ) return 0;		// 0-10%
		if ( 2.04651	<=centrality && centrality< 2.98931  ) return 1;	// 10-20%
		if ( 1.36875	<=centrality && centrality< 2.04651  ) return 2;	// 20-30%
		if ( 0.87541	<=centrality && centrality< 1.36875  ) return 3;	// 30-40%
		if ( 0.525092	<=centrality && centrality< 0.87541  ) return 4;	// 40-50%
		if ( 0.289595	<=centrality && centrality< 0.525092 ) return 4;	// 50-60%
		if ( 0.14414	<=centrality && centrality< 0.289595 ) return 5;	// 60-70%
		if ( 0.063719	<=centrality && centrality< 0.14414 ) return 6;	// 70-80%
		
		// Hijing doesn't have the same FCal distribution as data, let's keep everything
		if(isMC && centrality>= 6.0) return 0;
		if(isMC && 0.0<=centrality && centrality<0.063719) return 5;
		
		return -1;
	}
	else if (centralityScheme==33) // Pb+Pb 2015, merged 40-60%, fine central bins
	{
		// nominal 85%, full Fcal
		if ( 3.61844 	<=centrality && centrality< 6.00  ) return 0;		// 0-5%
		if ( 2.98931 	<=centrality && centrality< 3.61844  ) return 1;	// 5-10%
		if ( 2.04651	<=centrality && centrality< 2.98931  ) return 2;	// 10-20%
		if ( 1.36875	<=centrality && centrality< 2.04651  ) return 3;	// 20-30%
		if ( 0.87541	<=centrality && centrality< 1.36875  ) return 4;	// 30-40%
		if ( 0.525092	<=centrality && centrality< 0.87541  ) return 5;	// 40-50%
		if ( 0.289595	<=centrality && centrality< 0.525092 ) return 5;	// 50-60%
		if ( 0.063719	<=centrality && centrality< 0.289595 ) return 6;	// 60-80%
		
		// Hijing doesn't have the same FCal distribution as data, let's keep everything
		if(isMC && centrality>= 6.0) return 0;
		if(isMC && 0.0<=centrality && centrality<0.063719) return 5;
		
		return -1;
	}
	else if (centralityScheme==34) // Pb+Pb 2015, 10% bins, I think this is correct for JER (P.B.)
	{
		// nominal 85%, full Fcal
		if ( 3.61844 	<=centrality && centrality< 6.00     ) return 0;	// 0-5%
		if ( 2.98931 	<=centrality && centrality< 3.61844  ) return 0;	// 5-10%
		if ( 2.04651	<=centrality && centrality< 2.98931  ) return 1;	// 10-20%
		if ( 1.36875	<=centrality && centrality< 2.04651  ) return 2;	// 20-30%
		if ( 0.87541	<=centrality && centrality< 1.36875  ) return 3;	// 30-40%
		if ( 0.525092	<=centrality && centrality< 0.87541  ) return 4;	// 40-50%
		if ( 0.289595	<=centrality && centrality< 0.525092 ) return 5;	// 50-60%
		if ( 0.14414	<=centrality && centrality< 0.289595 ) return 6;	// 60-70%
		if ( 0.063719	<=centrality && centrality< 0.14414  ) return 7;	// 70-80%
		
		// Hijing doesn't have the same FCal distribution as data, let's keep everything
		if(isMC && centrality>= 6.0) return 0;
		if(isMC && 0.0<=centrality && centrality<0.063719) return 7;
		
		return -1;
	}
	else if (centralityScheme==35) // Pb+Pb 2015, coarse bins for test of y dependence
	{
		if ( 2.98931 	<=centrality && centrality< 6.00  ) return 0;		// 0-10%
		if ( 2.04651	<=centrality && centrality< 2.98931  ) return 0;	// 10-20%
		if ( 1.36875	<=centrality && centrality< 2.04651  ) return 0;	// 20-30%
		if ( 0.87541	<=centrality && centrality< 1.36875  ) return 1;	// 30-40%
		if ( 0.525092	<=centrality && centrality< 0.87541  ) return 1;	// 40-50%
		if ( 0.289595	<=centrality && centrality< 0.525092 ) return 1;	// 50-60%
		if ( 0.14414	<=centrality && centrality< 0.289595 ) return 1;	// 60-70%
		if ( 0.063719	<=centrality && centrality< 0.14414 ) return 1;	// 70-80%
		
		// Hijing doesn't have the same FCal distribution as data, let's keep everything
		if(isMC && centrality>= 6.0) return 0;
		if(isMC && 0.0<=centrality && centrality<0.063719) return 1;
		
		return -1;
	}
	else if (centralityScheme==20)	// p+Pb centrality
	{
		centrality = FCal_Et;

		if (53.74 <= centrality && centrality < 10e9 ) return 0; //0-10%
		if (40.04 <= centrality && centrality < 53.74) return 1; //10-20%
		if (31.07 <= centrality && centrality < 40.04) return 2; //20-30%
		if (24.10 <= centrality && centrality < 31.07) return 3; //30-40%
		if (13.41 <= centrality && centrality < 24.10) return 4; //40-60%
		if (5.585 <= centrality && centrality < 13.41) return 5; //60-90%

		return -1;
	}
	else if (centralityScheme==40)	// p+Pb, 2016, 8.16GeV centrality
	{
		centrality = FCal_Et;
		
		// this is stolen from p+Pb 2013 and basically crap, need to be changed for any real analysis
		if (53.74 <= centrality && centrality < 10e9 ) return 0; //0-10%
		if (40.04 <= centrality && centrality < 53.74) return 1; //10-20%
		if (31.07 <= centrality && centrality < 40.04) return 2; //20-30%
		if (24.10 <= centrality && centrality < 31.07) return 3; //30-40%
		if (13.41 <= centrality && centrality < 24.10) return 4; //40-60%
		if (5.585 <= centrality && centrality < 13.41) return 5; //60-90%
		if (0.000 <= centrality && centrality < 5.585) return 6; //90-100%
		
		return -1;
	}
	else if (centralityScheme==50)	// XeXe, 2017, 5.44TeV centrality
	{
		centrality = FCal_Et;
		
		// this is stolen from p+Pb 2013 and basically crap, need to be changed for any real analysis
		if (1.88743 <= centrality && centrality < 10e9 ) return 0; //0-10%
		if (1.30248 <= centrality && centrality < 1.88743) return 1; //10-20%
		if (0.88055 <= centrality && centrality < 1.30248) return 2; //20-30%
		if (0.57184 <= centrality && centrality < 0.88055) return 3; //30-40%
		if (0.35113 <= centrality && centrality < 0.57184) return 4; //40-50%
		if (0.20038 <= centrality && centrality < 0.35113) return 5; //50-60%
		if (0.10461 <= centrality && centrality < 0.20038) return 6; //60-70%
		if (0.04900 <= centrality && centrality < 0.10461) return 7; //70-80%
		if (0.000 <= centrality && centrality < 0.04900) return 8; //80-100%
		
		return -1;
	}
	else {
		centrality = 0;
		return -2;
	}
}

int GetCentralityNBins(Int_t centralityScheme)
{
	//Number + 1 to  include the inclusive bin
	if (centralityScheme==1) return 1;
	if (centralityScheme==2) return 7;
	if (centralityScheme==20) return 7;
	if (centralityScheme==30) return 8;
	if (centralityScheme==31) return 7;
	if (centralityScheme==32) return 8;
	if (centralityScheme==33) return 8;
	if (centralityScheme==34) return 9;
	if (centralityScheme==35) return 3;
	if (centralityScheme==40) return 8;
	if (centralityScheme==50) return 9;
	
	else return 1;
}

