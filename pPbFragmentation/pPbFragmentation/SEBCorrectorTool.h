#include <vector>

#include <iostream>

#include "TMath.h"

class SEBCorrectorTool {

 private:
  // this keeps track of whether each tower is claimed by a seed (1) or not (0)
  int _seed_map[100][64];
  // this keeps track of whether each tower is claimed by a jet (jet inded #) or not (-1)
  int _jet_map[100][64];
  // for a tower that is owned by a jet _only_, this is the ET in the
  // tower as reported by the jet constituents. for a tower that is
  // owned by a seed or unclaimed, this is zero (because the energy
  // there is not needed);
  // values in GeV, at the uncalibrated/EM scale
  float _energy_map[100][64];
  
  float _dR(float eta1, float eta2, float phi1, float phi2) {
    // assumed phi1, phi2 within -TMath::Pi() to TMath::Pi()
    float dphi = phi1 - phi2;
    if (dphi > TMath::Pi()) dphi -= 2*TMath::Pi();
    else if (dphi < -TMath::Pi()) dphi += 2*TMath::Pi();
    return sqrt((eta1 - eta2) * (eta1 - eta2) + dphi * dphi);
  }
  
  int _etaToBin(float eta) {
    eta += 4.95;
    eta = eta * 10;
    return (int) (eta + 0.1);
  }
  
  int _phiToBin(float phi) {
    phi += 31.5/32.0 * TMath::Pi();
    phi = phi / (TMath::Pi() / 32.0);
    return (int) (phi + 0.1);
  }
  
 public:
  
  // partial self and mutual bias
  std::vector<float> GetCorrections(
				    std::vector<std::vector<float> > *jet_const_eta, 
				    std::vector<std::vector<float> > *jet_const_phi, 
				    // this procedure is agnostic about the
				    // choice of units (MeV vs GeV) and will
				    // return corrections in the same units as
				    // what it is given
				    std::vector<std::vector<float> > *jet_const_ET, 
				    std::vector<float> *seed_eta,
				    std::vector<float> *seed_phi
				    ) 
    {

      // initialize the vector to be returned (there is a correction for
      // each jet)
      unsigned int Njets = jet_const_eta->size();

      std::vector<float> PSMB_corrections(Njets, 0);

      // in general, indexing is 100 eta bins from -5 to +5 and 64 phi bins from -TMath::Pi() to +TMath::Pi()
      
      // this procedure is O(N_seed) but with 6400 steps per seed
      
      // initialize the map to be unclaimed or covered by a seed
      for (int etabin = 0; etabin < 100; etabin++) {
	for (int phibin = 0; phibin < 64; phibin++) {
	  float this_eta = -5 + 0.05 + 0.1*etabin;
	  float this_phi = -TMath::Pi() + TMath::Pi()/64.0 + TMath::Pi()/32.0*phibin;
	  
	  // unclaimed by default
	  _seed_map[etabin][phibin] = 0;
	  _jet_map[etabin][phibin] = -1;
	  _energy_map[etabin][phibin] = 0.0;
	  
	  // iterate through seeds
	  for (unsigned int iseed = 0; iseed < seed_eta->size(); iseed++) {
	    if (_dR((*seed_eta)[iseed], this_eta, (*seed_phi)[iseed], this_phi) < 0.4) {
	      _seed_map[etabin][phibin] = 1;
	    }
	  }
	}
      }
      
      // this next procedure is O(N_jet) (really, O(N_jet_constituents) )
      
      
      // now, iterate through each jet
      for (unsigned int ijet = 0; ijet < Njets; ijet++) {
	// iterate through each jet constituent
	for (unsigned int ic = 0; ic < (*jet_const_eta)[ijet].size(); ic++) {
	  int etabin = _etaToBin((*jet_const_eta)[ijet][ic]);
	  int phibin = _phiToBin((*jet_const_phi)[ijet][ic]);
	  
	  if (_jet_map[etabin][phibin] != -1) {
	    // this tower is already owned by another jet!!  this should
	    // never happen, fail loudly
	    std::cout << " constituent #" << ic << " in jet #" << ijet << " is already owned by jet #" << _jet_map[etabin][phibin] << "!" << std::endl;
	    std::cout << " exiting procedure; check your jet collections and try again" << std::endl;
	    return PSMB_corrections;
	  }
	  
	  _jet_map[etabin][phibin] = ijet;
	  _energy_map[etabin][phibin] = (*jet_const_ET)[ijet][ic];
	  
	}
	
      }
      
      // now we are ready to compute the corrections
      
      // now, perform the corrections eta strip by strip
      for (int etabin = 0; etabin < 100; etabin++) {
	// in my slides, the following variables are needed:
	
	// n_{excl}, the towers this eta strip that have been excluded
	// from the background determination by the nearby presence of a
	// seed
	int total_seed_towers = 0;
	// n_{i}, the towers this eta strip that belong to a given jet
	std::vector<int> jet_num_towers_affected(Njets, 0);
	// \Sum E_{t}, the total energy in non-excluded jet towers
	float total_affected_energy = 0.0;
	
	// now we go through the towers, counting 
	for (int phibin = 0; phibin < 64; phibin++) {
	  int jet_owner = _jet_map[etabin][phibin];
	  int seed_owner = _seed_map[etabin][phibin];
	  if (jet_owner > -1)
	    jet_num_towers_affected[jet_owner]++;
	  if (seed_owner == 1)
	    total_seed_towers++;
	  if (jet_owner > -1 && seed_owner == 0)	
	    total_affected_energy += _energy_map[etabin][phibin];
	}
	
	// compute the correction per jet
	
	for (unsigned int ijet = 0; ijet < Njets; ijet++) {
	  if (jet_num_towers_affected[ijet] > 0) {
	    float correction = total_affected_energy * jet_num_towers_affected[ijet] / ( 64.0 - total_seed_towers);
	    // accumulate the corrections in this eta strip into the total
	    // correction
	    PSMB_corrections[ijet] += correction;
	  }
	}
	
      }
      
      // now we have integrated over all eta strips
      // return total correction
      return PSMB_corrections;

    }

};

