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
