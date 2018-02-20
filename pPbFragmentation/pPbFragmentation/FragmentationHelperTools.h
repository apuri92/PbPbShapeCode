#ifndef __FRAGMENTATIONHELPERTOOLS_H__
#define __FRAGMENTATIONHELPERTOOLS_H__


#include <string>
#include <iostream>
#include <vector>
#include <TTree.h>
#include "pPbFragmentation/JetHelperTools.h"
#include "xAODJet/JetContainer.h"
#include "xAODTracking/TrackParticleContainer.h"
#include "xAODTracking/VertexContainer.h"
#include "InDetTrackSelectionTool/InDetTrackSelectionTool.h"
using namespace std;

namespace MTCorrector
{
	std::vector<float> GetIsolation(std::vector<float>& jetEt, std::vector<float>& jetEta, std::vector<float>& jetPhi, int jetRadius);
	std::vector<bool> GetIsolation(std::vector<float>& jetpT, std::vector<float>& jetEta, std::vector<float>& jetPhi, int iso_R, double iso_pT);
	void SetupBinning(Int_t scheme, string variable, Double_t array[1000], Int_t &num);
	int GetTrkpTBin(float pt);
	std::vector<int> TruthMatching(std::vector<float>& reco_jetpT, std::vector<float>& reco_jetEta, std::vector<float>& reco_jetPhi,
								   std::vector<float>& truth_jetpT, std::vector<float>& truth_jetEta, std::vector<float>& truth_jetPhi,
								   float match_R);
	float R2Matching(float &reco_jetEta, float &reco_jetPhi,
					 std::vector<float>& a2_jetEta, std::vector<float>& a2_jetPhi,
					 int match_R);

	double GetZ(double track_pt, double track_eta, double track_phi, double jet_pt, double jet_eta, double jet_phi, bool UseAltzDef);
	bool GetFJR(float jetEta, float jetPhi, const xAOD::JetContainer* track_jets);
	bool GetFJR(xAOD::Jet* jet, InDet::InDetTrackSelectionTool *m_track_selection_tool);
	int GetRunNumberBin(int RunNumber);
};
#endif
