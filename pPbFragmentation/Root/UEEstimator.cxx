#include "pPbFragmentation/UEEstimator.h"

void UEEstimator::ExcludeConesByTrk(vector<float> &trk_pt,vector<float> &trk_eta,vector<float> &trk_phi)
//@brief: this is to indentify which cones have to be excluded due to possibly containing a jet 
{
   for (int i=0; i<_maxNCones; i++) _bkgrCones[i] = 1;	// any cone can be used unless we prove otherwise
   for (int i=0; i<_maxNCones; i++) _bkgrCones_hpT[i] = 0.;	// reseting maximal pT

   for (unsigned int j1=0; j1<trk_pt.size(); j1++)
     {
      //if (trk.ApplyAllCuts(j1)==0) continue; //?????
      //if (trk_pt.at(j1)/1000. < ptBkgrThreshold ) continue;

      for (int iEta=0; iEta<_nEta; iEta++)
        {for (int iPhi=0; iPhi<_nPhi; iPhi++)
		{	Float_t thePhi = cone_phi[iPhi];
			Float_t theEta = cone_eta[iEta];
            Float_t deltaR = DeltaR( trk_phi.at(j1), trk_eta.at(j1), thePhi, theEta );
            Int_t idx = _nEta * iPhi + iEta;
            if (deltaR<cone_radius)
              {
               if (trk_pt.at(j1) > ptBkgrThreshold ) _bkgrCones[idx] = 0;   // disable this cone since there is a jet
               if (trk_pt.at(j1) > _bkgrCones_hpT[idx] ) _bkgrCones_hpT[idx]= trk_pt.at(j1);
               }
           }
        }
     }

   _w_ncones = _maxNCones;
   for (int i=0; i<_maxNCones; i++)
     {if (_bkgrCones[i]==0) _w_ncones -= 1;
     }
   //cout << "_w_ncones (track rejection) : " << _w_ncones << endl;
}

void UEEstimator::ExcludeConesByJet(vector<float> &jet_pt,vector<float> &jet_eta,vector<float> &jet_phi)
//@brief: this is to indentify which cones have to be excluded due to possibly containing a jet 
{
   for (int i=0; i<_maxNCones; i++) _bkgrCones[i] = 1;	// any cone can be used unless we prove otherwise

   for (unsigned int j1=0; j1<jet_pt.size(); j1++)
     {
      if (jet_pt.at(j1) < jetptBkgrThreshold ) continue;

      for (int iEta=0; iEta<_nEta; iEta++)
        {for (int iPhi=0; iPhi<_nPhi; iPhi++)
		{	Float_t thePhi = cone_phi[iPhi];
			Float_t theEta = cone_eta[iEta];
            Float_t deltaR = DeltaR( jet_phi.at(j1), jet_eta.at(j1), thePhi, theEta );
            Int_t idx = _nEta * iPhi + iEta;
            if (deltaR<m_maxjetdeltaR)
              {
               _bkgrCones[idx] = 0;                           // disable this cone since there is a jet
              }
           }
        }
     }

   _w_ncones = _maxNCones;
   for (int i=0; i<_maxNCones; i++)
     {if (_bkgrCones[i]==0) _w_ncones -= 1;
     }
   //cout << "_w_ncones (jet rejection) : " << _w_ncones << endl;
}

void UEEstimator::ExcludeConesByJetandTrack(vector<float> &trk_pt,vector<float> &trk_eta,vector<float> &trk_phi, vector<float> &jet_pt,vector<float> &jet_eta,vector<float> &jet_phi)
//@brief: this is to indentify which cones have to be excluded due to possibly containing a jet 
{
   for (int i=0; i<_maxNCones; i++) _bkgrCones[i] = 1;	// any cone can be used unless we prove otherwise

   for (unsigned int j1=0; j1<jet_pt.size(); j1++)
     {
      if (jet_pt.at(j1) < jetptBkgrThreshold ) continue;

      for (int iEta=0; iEta<_nEta; iEta++)
        {for (int iPhi=0; iPhi<_nPhi; iPhi++)
		{	Float_t thePhi = cone_phi[iPhi];
			Float_t theEta = cone_eta[iEta];
            Float_t deltaR = DeltaR( jet_phi.at(j1), jet_eta.at(j1), thePhi, theEta );
            Int_t idx = _nEta * iPhi + iEta;
            if (deltaR<m_maxjetdeltaR)
              {
               _bkgrCones[idx] = 0;                           // disable this cone since there is a jet
//				  cout << Form("jet_%i: %f, %f, %f	cone: %f, %f	excluded:%i",j1, jet_eta.at(j1), jet_phi.at(j1), jet_pt.at(j1), theEta, thePhi, idx) << endl;
              }
           }
        }
     }
   for (unsigned int j1=0; j1<trk_pt.size(); j1++)
     {
      if (trk_pt.at(j1) < ptBkgrThreshold ) continue;

      for (int iEta=0; iEta<_nEta; iEta++)
        {for (int iPhi=0; iPhi<_nPhi; iPhi++)
		{	Float_t thePhi = cone_phi[iPhi];
			Float_t theEta = cone_eta[iEta];
            Float_t deltaR = DeltaR( trk_phi.at(j1), trk_eta.at(j1), thePhi, theEta );
            Int_t idx = _nEta * iPhi + iEta;
            if (deltaR<cone_radius)
              {
               _bkgrCones[idx] = 0;                           // disable this cone since there is a jet
//				  cout << Form("trk_%i: %f, %f, %f	cone: %f, %f	excluded:%i",j1, trk_eta.at(j1), trk_phi.at(j1), trk_pt.at(j1), theEta, thePhi, idx) << endl;
               if (trk_pt.at(j1) > _bkgrCones_hpT[idx] ) _bkgrCones_hpT[idx]= trk_pt.at(j1);
              }
           }
        }
     } 
   _w_ncones = _maxNCones;
   for (int i=0; i<_maxNCones; i++)
     {if (_bkgrCones[i]==0) _w_ncones -= 1;
     }
   //cout << "_w_ncones (jet rejection and track) : " << _w_ncones << endl;
}

void UEEstimator::FindCone(float trk_pt,float trk_eta,float trk_phi)
//@brief: find a minimum deltaR between a track and some non-disabled cone, that is find the cone
//@       that will be associated with this background particle, also store this deltaR in _deltaRToConeAxis 
{
   Float_t deltaROrthMin=999.;
   for (int iEta=0; iEta<_nEta; iEta++)
     {for (int iPhi=0; iPhi<_nPhi; iPhi++)
        {if (! _bkgrCones[iPhi*_nEta+iEta]) continue;
			Float_t thePhi = cone_phi[iPhi];
			Float_t theEta = cone_eta[iEta];
			Float_t deltaR = DeltaR( trk_phi, trk_eta, thePhi, theEta );
         if (deltaR < deltaROrthMin)
           {deltaROrthMin = deltaR;
            _etaOfCone = theEta;
            _phiOfCone = thePhi;
            _maxConePt = _bkgrCones_hpT[iPhi*_nEta+iEta];
            _maxConeIndex = iPhi*_nEta+iEta;
           }
        }
     }

   _deltaRToConeAxis = deltaROrthMin;
//	cout << Form("trk: %f, %f-%f, cone_%i, %f-%f, r: %f", trk_pt, trk_eta, trk_phi, _maxConeIndex, _etaOfCone, _phiOfCone, _deltaRToConeAxis) << endl;

}


Float_t UEEstimator::CalculateEtaWeight(float trk_pT, float trk_eta, float jet_eta, Int_t icent)
//@brief: calculate a weight that is due to a difference in the yield(eta) between a jet 
//@       position and a position of a given particle that is used to estimate the UE
//@note:  naming of some variable is not optimal as we don't want to deviate from previous codes
{
   
   float pT_temp = trk_pT;
   if (pT_temp < 1.) pT_temp = 1.; // weight in bin bellow 1 GeV is taken as weight at 1 GeV 
   int pT_bin=ptaxis->FindBin(pT_temp);
   //cout << "pt: " << trk_pT <<  " Pt bin" << pT_bin << " eta trk " <<  trk_eta << " jet eta" << jet_eta <<  " cent " << icent <<  " Period " << period << endl;
   Float_t trkEta = trk_eta;
   Float_t nearJetEta = jet_eta;
   Float_t deltaEtaOrthMin = trkEta - _etaOfCone;	// (old version: trkEta - theEta)
   deltaEtaOrthMin *= ((_etaOfCone!=0)&&(nearJetEta!=0))? ( ((_etaOfCone/fabs(_etaOfCone))==(nearJetEta/fabs(nearJetEta)))? 1:-1 ):1;
   //cout << "trk_pT " << trk_pT << " pT_bin " << pT_bin << " cent " << icent << endl;
    Float_t w_eta = _h_eta_w[pT_bin][icent]->Interpolate( nearJetEta + deltaEtaOrthMin ) / _h_eta_w[pT_bin][icent]->Interpolate( trkEta );
   //Float_t w_eta = _f1_trkEta[icent]->Eval( nearJetEta + deltaEtaOrthMin ) / _f1_trkEta[icent]->Eval( trkEta );

   return w_eta;
}

Float_t UEEstimator::CalculateEtaWeight_circle(float trk_pT, float trk_eta, float jet_eta, float circle_eta, Int_t icent)
//@brief: calculate a weight that is due to a difference in the yield(eta) between a jet
//@       position and a position of a given particle that is used to estimate the UE
//@note:  naming of some variable is not optimal as we don't want to deviate from previous codes
{

	float pT_temp = trk_pT;
	if (pT_temp < 1.) pT_temp = 1.; // weight in bin bellow 1 GeV is taken as weight at 1 GeV
	int pT_bin=ptaxis->FindBin(pT_temp);
	//cout << "pt: " << trk_pT <<  " Pt bin" << pT_bin << " eta trk " <<  trk_eta << " jet eta" << jet_eta <<  " cent " << icent <<  " Period " << period << endl;
	Float_t trkEta = trk_eta;
	Float_t nearJetEta = jet_eta;
	Float_t deltaEtaOrthMin = trkEta - circle_eta;	// (old version: trkEta - theEta)
	deltaEtaOrthMin *= ((circle_eta!=0)&&(nearJetEta!=0))? ( ((circle_eta/fabs(circle_eta))==(nearJetEta/fabs(nearJetEta)))? 1:-1 ):1;
	//cout << "trk_pT " << trk_pT << " pT_bin " << pT_bin << " cent " << icent << endl;
	Float_t w_eta = _h_eta_w[pT_bin][icent]->Interpolate( nearJetEta + deltaEtaOrthMin ) / _h_eta_w[pT_bin][icent]->Interpolate( trkEta );
	//Float_t w_eta = _f1_trkEta[icent]->Eval( nearJetEta + deltaEtaOrthMin ) / _f1_trkEta[icent]->Eval( trkEta );

	return w_eta;
}


Float_t UEEstimator::CalculateFlowWeight(float trk_pt,float trk_eta,float trk_phi, float nearJetPhi, float FCalEt)
//@brief: calculate a weight that is due to a difference in the flow between a jet 
//@       position and a position of a given particle that is used to estimate the UE
//@note:  naming of some variable is not optimal as we don't want to deviate from previous codes
{


   Float_t v2 = _h_v2_EP->GetBinContent(_h_v2_EP->GetXaxis()->FindBin(FCalEt),_h_v2_EP->GetYaxis()->FindBin(trk_pt),_h_v2_EP->GetZaxis()->FindBin(trk_eta));

			// calculate the event plane  
   Float_t w_flow  = (v2*cos(2*GetDeltaPsi( trk_phi, Psi))!=-0.5)?
                 ( (1 + 2*v2*cos(2*GetDeltaPsi( nearJetPhi, Psi)) )
                  /(1 + 2*v2*cos(2*GetDeltaPsi( trk_phi, Psi)) ) ):0; 	// modulate/demodulate

   if (w_flow > 1.4) w_flow = 1.4;	// if the correction is too large ...
   if (w_flow < 0.6) w_flow = 0.6;      // ... or too small, then take some maximal/minimal value

   return w_flow;
}

Float_t UEEstimator::CalculateV3Weight(float trk_pt,float trk_eta,float trk_phi, float nearJetPhi, float FCalEt)
//@brief: calculate a weight that is due to a difference in the flow between a jet 
//@       position and a position of a given particle that is used to estimate the UE
//@note:  naming of some variable is not optimal as we don't want to deviate from previous codes
{

   Float_t v3 = _h_v3_EP->GetBinContent(_h_v3_EP->GetXaxis()->FindBin(FCalEt),_h_v3_EP->GetYaxis()->FindBin(trk_pt),_h_v3_EP->GetZaxis()->FindBin(trk_eta));
			// calculate the event plane  
   Float_t w_flow  = (v3*cos(3*GetDeltaPsi3( trk_phi, Psi3))!=-0.5)?
                 ( (1 + 2*v3*cos(3*GetDeltaPsi3( nearJetPhi, Psi3)) )
                  /(1 + 2*v3*cos(3*GetDeltaPsi3( trk_phi, Psi3)) ) ):0; 	// modulate/demodulate
   if (w_flow > 1.4) w_flow = 1.4;	// if the correction is too large ...
   if (w_flow < 0.6) w_flow = 0.6;      // ... or too small, then take some maximal/minimal value

   return w_flow;
}

Float_t UEEstimator::GetDeltaPsi(double phi, double psi)
//@brief: distance from the event plane in convention of flow calculations
{
    Double_t diff;
    diff = fabs(phi - psi);
    while (diff > TMath::Pi()/2. ) diff = TMath::Pi() - diff;
    return fabs(diff);
}

Float_t UEEstimator::GetDeltaPsi3(double phi, double psi)
//@brief: distance from the event plane in convention of flow calculations
{
    Double_t diff;
    diff = fabs(phi - psi);
    //cout << "diff " << diff;
    //if (diff>1.047) diff=diff-0.002;
    //if (diff<-1.047) diff=diff+0.002;
    //cout << " diff " << diff <<endl;
    while (diff > TMath::Pi()/3. ) diff = 2*TMath::Pi()/3. - diff;
    return fabs(diff);
}

void UEEstimator::initShapeUE(bool isMC)
{

	for (int i_dR = 0; i_dR < 13; i_dR++)
	{
		for (int i_dPsi = 0; i_dPsi < 16; i_dPsi++)
		{
			for (int i_pt = 0; i_pt < 10; i_pt++)
			{
				for (int i_cent = 0; i_cent < 6; i_cent++)
				{
					std::string name;
					//if (isMC) name = Form("h_UE_MC_dR%i_dPsi%i_pt%i_cent%i", i_dR, i_dPsi, i_pt+1, i_cent);
					//else
					name = Form("h_UE_HP_dR%i_dPsi%i_pt%i_cent%i", i_dR, i_dPsi, i_pt+1, i_cent);
					h_UE[i_dR][i_dPsi][i_pt][i_cent] = (TH2*)_f_ShapeUE->Get(name.c_str());
				}
			}
		}
	}
}

void UEEstimator::initShapeUE(bool isMC, int uncert)
{
	_f_ShapeUE = new TFile("$ROOTCOREBIN/../pPbFragmentation/data/UE_MC_comb.root","READ");
	if ((uncert >= 10 && uncert <= 15) || uncert == 2 || uncert == 5)
	{
		_f_ShapeUE = new TFile(Form("$ROOTCOREBIN/../pPbFragmentation/data/UE_MC_comb_sys%i.root",uncert),"READ");
	}

	cout << "Using UE file: " << _f_ShapeUE->GetName() << endl;

	//UE with only v2 initialization
	for (int i_dR = 0; i_dR < 13; i_dR++)
	{
		for (int i_dPsi = 0; i_dPsi < 10; i_dPsi++)
		{
			for (int i_pt = 0; i_pt < 7; i_pt++)
			{
				for (int i_cent = 0; i_cent < 6; i_cent++)
				{
					for (int i_jet = 6; i_jet < 12; i_jet++)
					{

						std::string name;
						name = Form("h_UE_new_MC_dR%i_dPsi%i_pt%i_cent%i_jet%i",i_dR, i_dPsi, i_pt+1, i_cent, i_jet);
						h_UE_eta_phi_maps[i_jet][i_pt][i_dPsi][i_cent][i_dR] = (TH2*)_f_ShapeUE->Get(name.c_str());
					}
				}
			}
		}
	}
}

double UEEstimator::getShapeUE(int i_dR, int i_dPsi, int i_pt, int i_cent, double jet_eta, double jet_phi, double &error)
{
	int bin_eta = h_UE[i_dR][i_dPsi][i_pt][i_cent]->GetXaxis()->FindBin(jet_eta);
	int bin_phi = h_UE[i_dR][i_dPsi][i_pt][i_cent]->GetYaxis()->FindBin(jet_phi);
	double val =  h_UE[i_dR][i_dPsi][i_pt][i_cent]->GetBinContent(bin_eta, bin_phi);
	error =  h_UE[i_dR][i_dPsi][i_pt][i_cent]->GetBinError(bin_eta, bin_phi);
	return val;
}


double UEEstimator::getShapeUE(bool UE_MC, int i_dR, int i_dPsi, int i_pt, int i_cent, double jet_eta, double jet_phi, int i_jet, double &error)
{
	int bin_eta = h_UE_eta_phi_maps[i_jet][i_pt][i_dPsi][i_cent][i_dR]->GetXaxis()->FindBin(jet_eta);
	int bin_phi = h_UE_eta_phi_maps[i_jet][i_pt][i_dPsi][i_cent][i_dR]->GetYaxis()->FindBin(jet_phi);
	double val =  h_UE_eta_phi_maps[i_jet][i_pt][i_dPsi][i_cent][i_dR]->GetBinContent(bin_eta, bin_phi);
	error =  h_UE_eta_phi_maps[i_jet][i_pt][i_dPsi][i_cent][i_dR]->GetBinError(bin_eta, bin_phi);
	return val;
}


void UEEstimator::InitCones()
{
	for (int iEta=0; iEta<_nEta; iEta++){ cone_eta[iEta] = ((cone_eta_start+cone_radius) + (2*cone_radius*iEta));}
	for (int iPhi=0; iPhi<_nPhi; iPhi++){ cone_phi[iPhi] = ((cone_phi_start+cone_radius) + (2*cone_radius*iPhi));}

	for (int i=0; i<_maxNCones; i++) _bkgrCones[i] = 1;	// all cones can be used by default
	for (int i=0; i<_maxNCones; i++) _bkgrCones_hpT[i] = 0;	// maximal track pT
}



