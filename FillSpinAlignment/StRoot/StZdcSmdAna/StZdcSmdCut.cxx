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
//---------------------------------------------------------------------------------
