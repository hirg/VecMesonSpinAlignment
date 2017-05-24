#include "StRoot/StZdcSmdAna/StZdcSmdCut.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "StMessMgr.h"

ClassImp(StZdcSmdCut)

//---------------------------------------------------------------------------------

StZdcSmdCut::StZdcSmdCut()
{
}

//---------------------------------------------------------------------------------

StZdcSmdCut::~StZdcSmdCut()
{
}

//---------------------------------------------------------------------------------
bool StZdcSmdCut::passPhiCut(TLorentzVector lTrack) 
{
  float eta = lTrack.Eta();
  
  // eta cut
  if(fabs(eta) > vmsa::mEtaMax) return kFALSE;

  return kTRUE;
}

bool StZdcSmdCut::passDipAngleCut(TLorentzVector lKplus, TLorentzVector lKminus)
{
  double KplusPt = lKplus.Pt();
  double KplusPz = lKplus.Pz();
  double KplusP  = lKplus.P();

  double KminusPt = lKminus.Pt(); 
  double KminusPz = lKminus.Pz(); 
  double KminusP  = lKminus.P();

  double costheta = (KplusPt*KminusPt+KplusPz*KminusPz)/(KplusP*KminusP);
  double theta = acos(costheta);
  if(theta < 0.04) return kFALSE;

  return kTRUE;
}
//---------------------------------------------------------------------------------
