#include "StEffCut.h"
#include <math.h>

ClassImp(StEffCut)

StEffCut::StEffCut()
{
}

StEffCut::~StEffCut()
{
  /* */
}

//-----------------------------------------------------------------------------------------------------
bool StEffCut::passTrackCutPhi(McVecMeson McPhi)
{
  if(fabs(McPhi.McEta) > vmsa::mEtaMax) return kFALSE; // eta cut

  return kTRUE;
}

bool StEffCut::passTrackCutPhi(RcVecMeson RcPhi)
{
  if(fabs(RcPhi.RcEta) > vmsa::mEtaMax) return kFALSE; // eta cut

  return kTRUE;
}
//-----------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------
bool StEffCut::passTrackCut(McDecayDau McKaon)
{
  if(fabs(McKaon.McEta) > vmsa::mEtaMax) return kFALSE; // eta cut for some reason

  float pt = McKaon.McPt;
  float momentum = McKaon.McPt*cosh(McKaon.McEta);
  // primary pt and momentum cut: PtMin > 0.1, PMax < 10.0
  if(!(pt > vmsa::mGlobPtMin && momentum < vmsa::mPrimMomMax)) return kFALSE;
  
  return kTRUE;
}

bool StEffCut::passTrackCut(RcDecayDau RcKaon)
{
  if(RcKaon.RcTpc < 0.5) return kFALSE; // no matched track
  if(fabs(RcKaon.RcEta) > vmsa::mEtaMax) return kFALSE; // eta cut

  float pt = RcKaon.RcPt;
  float momentum = RcKaon.RcPt*cosh(RcKaon.RcEta);
  // primary pt and momentum cut: PtMin = 0.1
  if(!(pt > vmsa::mGlobPtMin && momentum < vmsa::mPrimMomMax)) return kFALSE;

  return kTRUE;
}
//-----------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------
bool StEffCut::passDipAngleCut(McDecayDau McKplus, McDecayDau McKminus)
{
  TLorentzVector lMcKplus;
  lMcKplus.SetPtEtaPhiM(McKplus.McPt,McKplus.McEta,McKplus.McPhi,vmsa::mMassKaon);
  double KplusPt = lMcKplus.Pt();
  double KplusPz = lMcKplus.Pz();
  double KplusP  = lMcKplus.P();

  TLorentzVector lMcKminus;
  lMcKminus.SetPtEtaPhiM(McKminus.McPt,McKminus.McEta,McKminus.McPhi,vmsa::mMassKaon);
  double KminusPt = lMcKminus.Pt();
  double KminusPz = lMcKminus.Pz();
  double KminusP  = lMcKminus.P();

  double costheta = (KplusPt*KminusPt+KplusPz*KminusPz)/(KplusP*KminusP);
  double theta = acos(costheta);
  if(theta < 0.04) return kFALSE;

  return kTRUE;
}

bool StEffCut::passDipAngleCut(RcDecayDau RcKplus, RcDecayDau RcKminus)
{
  TLorentzVector lRcKplus;
  lRcKplus.SetPtEtaPhiM(RcKplus.RcPt,RcKplus.RcEta,RcKplus.RcPhi,vmsa::mMassKaon);
  double KplusPt = lRcKplus.Pt();
  double KplusPz = lRcKplus.Pz();
  double KplusP  = lRcKplus.P();

  TLorentzVector lRcKminus;
  lMcKminus.SetPtEtaPhiM(RcKminus.RcPt,RcKminus.RcEta,RcKminus.RcPhi,vmsa::mMassKaon);
  double KminusPt = lRcKminus.Pt();
  double KminusPz = lRcKminus.Pz();
  double KminusP  = lRcKminus.P();

  double costheta = (KplusPt*KminusPt+KplusPz*KminusPz)/(KplusP*KminusP);
  double theta = acos(costheta);
  if(theta < 0.04) return kFALSE;

  return kTRUE;
}
//-----------------------------------------------------------------------------------------------------
