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
  double KplusPt = McKplus.McPt;
  double KplusPz = KplusPt*sinh(McKplus.McEta);
  double KplusP  = sqrt(KplusPt*KplusPt+KplusPz*KplusPz);

  double KminusPt = McKminus.McPt;
  double KminusPz = KminusPz*sinh(McKminus.McEta);
  double KminusP  = sqrt(KminusPt*KminusPt+KminusPz*KminusPz);

  double costheta = (KplusPt*KminusPt+KplusPz*KminusPz)/(KplusP*KminusP);
  double theta = acos(costheta);
  if(theta < 0.04) return kFALSE;

  return kTRUE;
}

bool StEffCut::passDipAngleCut(RcDecayDau RcKplus, RcDecayDau RcKminus)
{
  double KplusPt = RcKplus.RcPt;
  double KplusPz = KplusPt*sinh(RcKplus.RcEta);
  double KplusP  = sqrt(KplusPt*KplusPt+KplusPz*KplusPz);

  double KminusPt = RcKminus.RcPt;
  double KminusPz = KminusPz*sinh(RcKminus.RcEta);
  double KminusP  = sqrt(KminusPt*KminusPt+KminusPz*KminusPz);

  double costheta = (KplusPt*KminusPt+KplusPz*KminusPz)/(KplusP*KminusP);
  double theta = acos(costheta);
  if(theta < 0.04) return kFALSE;

  return kTRUE;
}
//-----------------------------------------------------------------------------------------------------
