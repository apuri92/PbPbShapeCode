#ifndef GET_GLOBALHELPERS_H
#define GET_GLOBALHELPERS_H
#include <TROOT.h>
#include <TH1D.h>
#include "xAODEventInfo/EventInfo.h"
//CaloCluster include
#include  "xAODCaloEvent/CaloCluster.h" 
#include  "xAODCaloEvent/CaloClusterContainer.h"
#define private public
#include "xAODHIEvent/HIEventShapeAuxContainer.h"
#undef private
#include "xAODHIEvent/HIEventShapeContainer.h" 

//@CODE_begin
int GetGlobalBin(Int_t centralityScheme, float FCal_Et, bool isMC=false);
int GetCentralityBin(Int_t centralityScheme, float FCal_Et, bool isMC=false);
void SetRejectionHistogram(TH1D* h);
int GetCentralityNBins(Int_t centralityScheme);
float GetEventPlane(const xAOD::CaloClusterContainer *hiclus);
double GetEventPlane(const xAOD::HIEventShapeContainer* calos);
Float_t GetAveragePsi(Float_t psi1, Float_t psi2);
Int_t GetPsiBin(Float_t psi);
double GetEventPlane(const xAOD::HIEventShapeContainer* calos, int order);
Float_t GetDeltaPsi3(double phi, double psi);
//@CODE_end

#endif
