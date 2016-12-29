#include "StVecMesonCut.h"
#include "../../../Utility/StSpinAlignmentCons.h"
#include "StRoot/StPicoDstMaker/StPicoTrack.h"
#include "StMessMgr.h"

ClassImp(StVecMesonCut)

//---------------------------------------------------------------------------------

StVecMesonCut::StVecMesonCut(Int_t energy)
{
  mEnergy = energy;
}

//---------------------------------------------------------------------------------

StVecMesonCut::~StVecMesonCut()
{
}

//---------------------------------------------------------------------------------
bool StVecMesonCut::passTrackEP(TLorentzVector lTrack, Float_t dca)
{
  // only used for primary track
  // dca cut for event plane reconstruction
  if(fabs(dca) > vmsa::mDcaEPMax[mEnergy])
  {
    return kFALSE;
  }

  // pt cut 0.2 - 2.0 GeV/c
  Float_t pt = lTrack.Perp();
  Float_t p = lTrack.P();
  if(!(pt > vmsa::mPrimPtMin[mEnergy] && pt < vmsa::mPrimPtMax && p < vmsa::mPrimMomMax))
  {
    return kFALSE;
  }

  // eta cut -> only important for particles reconstructed by global tracks
  Float_t eta = lTrack.Eta();
  if(fabs(eta) > vmsa::mEtaMax)
  {
    return kFALSE;
  }

  return kTRUE;
}
//---------------------------------------------------------------------------------
bool StVecMesonCut::passTrackEtaEast(TLorentzVector lTrack) // neg
{
  Float_t eta = lTrack.Eta();
  
  // eta cut
  // eta_gap between two sub event plane is 2*mEta_Gap
  if(!(eta > -1.0*vmsa::mEtaMax && eta < -1.0*vmsa::mEta_Gap))
  {
    return kFALSE;
  }

  return kTRUE;
}

bool StVecMesonCut::passTrackEtaWest(TLorentzVector lTrack) // pos
{
  Float_t eta = lTrack.Eta();

  // eta cut
  // eta_gap between two sub event plane is 2*mEta_Gap
  if(!(eta > vmsa::mEta_Gap && eta < vmsa::mEtaMax))
  {
    return kFALSE;
  }

  return kTRUE;
}

bool StVecMesonCut::passPhiEtaEast(TLorentzVector lTrack) // neg
{
  Float_t eta = lTrack.Eta();
  
  // eta cut
  // eta_gap between mother particle and sub event plane is mEta_Gap
  if(!(eta > -1.0*vmsa::mEtaMax && eta < 0.0))
  {
    return kFALSE;
  }

  return kTRUE;
}

bool StVecMesonCut::passPhiEtaWest(TLorentzVector lTrack) // pos
{
  Float_t eta = lTrack.Eta();

  // eta cut
  // eta_gap between mother particle and sub event plane is mEta_Gap
  if(!(eta > 0.0 && eta < vmsa::mEtaMax))
  {
    return kFALSE;
  }

  return kTRUE;
}
//---------------------------------------------------------------------------------
