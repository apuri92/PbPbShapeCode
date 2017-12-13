#define EventSelectors_cxx
#include "pPbFragmentation/pPbFragmentation.h"

using namespace std;

//Function to check pileup
bool BaseClass::HasPileUp() const
{
  if (vx_n < 3) return false;

  int vertexCount = 0;

  for (int ivx = 0; ivx < vx_n - 1; ivx++) {
    if (vx_sumPt.at(ivx)/1000 > 6 || vx_nTracks.at(ivx) > 10) {
      vertexCount++;
    }
  }
  return vertexCount > 1;
}

//Function to check primary vertex
bool BaseClass::HasVertex() const
{
  return vx_n >  1;
}

//Function to check timing
bool BaseClass::HasGoodTiming() const
{
  if ( (mbtime_timeA  != 0 ) && (mbtime_timeC  != 0 ) &&
	    (std::abs(mbtime_timeA) != 75) && (std::abs(mbtime_timeC) != 75) &&
	    (std::abs(mbtime_timeA - mbtime_timeC) < 10) ) return true;
  return false;
}


std::vector<float> BaseClass::FilterJetsByMuons(std::vector<float>& jetEta, std::vector<float>& jetPhi, std::pair< xAOD::MuonContainer*, xAOD::ShallowAuxContainer* >& muon_container, Float_t deltaRMin)
{
	std::vector<float> muon_pT;
	muon_pT.clear();
     
	xAOD::MuonContainer::const_iterator muon_itr = ( muon_container.first)->begin();
	xAOD::MuonContainer::const_iterator muon_end = ( muon_container.first)->end();
     
	for(unsigned int i=0; i<jetEta.size(); i++)
	{
		float max_mu_pT = 0.;
		for( ; muon_itr != muon_end; ++muon_itr )
		{
			if(!m_muonSelection->accept(**muon_itr)) continue;
			float R = DeltaR(jetPhi.at(i),jetEta.at(i),(*muon_itr)->phi(),(*muon_itr)->eta());
			float mu_pT = ( (*muon_itr)->pt()/1000.);
			if (R > deltaRMin || max_mu_pT >  mu_pT) continue;
			max_mu_pT = mu_pT;
		}
		muon_pT.push_back(max_mu_pT);
	}
	return muon_pT;
}

std::vector<float> BaseClass::FilterJetsByElectrons(std::vector<float>& jetEta, std::vector<float>& jetPhi, std::pair< xAOD::ElectronContainer*, xAOD::ShallowAuxContainer* >& electron_container, Float_t deltaRMin){
     std::vector<float> electron_pT;
     electron_pT.clear();
     
     for (unsigned int i=0; i<jetEta.size(); i++){
      	 float max_el_pT = 0.;
      	 for( xAOD::Electron_v1* shallowElectron : *(electron_container.first) ){
      	   if ( m_LHToolTight2015->accept( shallowElectron ) == false ) continue;   // check if passes electron id tool
      	   // correction tool
		   if( m_EgammaCalibrationAndSmearingTool->applyCorrection( *shallowElectron ) != CP::CorrectionCode::Ok ) {
			  Warning("execute()", "Problem in CP::EgammaCalibrationAndSmearingTool::applyCorrection()");
		   }
 
		   // central electron 4 vector
		   TLorentzVector eCent4v = shallowElectron->p4();
      	   
      	   // cuts 
		   // check if inside kinematic acceptance
		   if ( fabs( shallowElectron->caloCluster()->etaBE(2) ) > 2.47 ) continue;
		   // check if inside "crack" region  
		   if ( fabs( shallowElectron->caloCluster()->etaBE(2) ) < 1.52 &&
			    fabs( shallowElectron->caloCluster()->etaBE(2) ) > 1.37 ) continue;
      	   
      	   
      	             
           float R = DeltaR(jetPhi.at(i),jetEta.at(i),eCent4v.Phi(),eCent4v.Eta());
           float el_pT = ( eCent4v.Pt()/1000.);
           if (R > deltaRMin || max_el_pT >  el_pT) continue;
           max_el_pT = el_pT;
        }
		electron_pT.push_back(max_el_pT);
     }
     return electron_pT;
}

