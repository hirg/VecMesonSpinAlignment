#include "StEffCut.h"
#include <math.h>

ClassImp(StEffCut)

StEffCut::StEffCut(int energy)
{
  mEnergy = energy;
}

StEffCut::~StEffCut()
{
  /* */
}

//-----------------------------------------------------------------------------------------------------
bool StEffCut::passEventCut(McEvent EventHeader)
{
  // vertex cut
  float vx = EventHeader.McVx;
  float vy = EventHeader.McVy;
  float vz = EventHeader.McVz;
  if(fabs(vz) > vmsa::mVzMaxMap[mEnergy]) // Vz cut
    return false;
  if(sqrt(vx*vx+vy*vy) > vmsa::mVrMax) //Vr cut
    return false;
  if(mEnergy == 6)
  {
    float vzVpd = EventHeader.vzVpd;
    if(fabs(vz-vzVpd) > vmsa::mVzVpdDiffMax) return false;
  }

  // centrality cut
  // if(EventHeader.centrality < vmsa::cent_low[0] || EventHeader.centrality > vmsa::cent_high[0]) // 20%-60%
  if(EventHeader.Centrality < 0) return false;

  return true;
}
//-----------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------
bool StEffCut::passTrackCut(McDecayDau McKaon)
{
  // eta cut for some reason
  if(fabs(McKaon.McEta) > vmsa::mEtaMax) return false;
  // if(fabs(McKaon.McY) > vmsa::mEtaMax) return false; // rapidity cut

  float pt = McKaon.McPt;
  float momentum = McKaon.McPt*cosh(McKaon.McEta);
  // primary pt and momentum cut: PtMin > 0.1, PMax < 10.0
  if(!(pt > vmsa::mGlobPtMin && momentum < vmsa::mPrimMomMax)) return false;
  
  return true;
}

bool StEffCut::passTrackCut(RcDecayDau RcKaon)
{
  // cut for global tracks
  // if(RcKaon.gRcPt < -900.) return false; // no matched track
  // if(RcKaon.gRcNfit < vmsa::mHitsFitTPCMin) return false; // nHitsFit cut
  // // if(RcKaon.gRcNdedx < vmsa::mHitsDedxMin) return false;
  // if(RcKaon.gRcNmax <= vmsa::mHitsMaxTPCMin) return false; // nHitsRatio cut
  // if((float)RcKaon.gRcNfit/(float)RcKaon.gRcNmax < vmsa::mHitsRatioTPCMin) return false;
  // // if(fabs(RcKaon.gRcEta) > vmsa::mEtaMax) return false; // eta cut
  // if(fabs(RcKaon.McEta) > vmsa::mEtaMax) return false; // eta cut
  // // float pt = RcKaon.gRcPt;
  // // float momentum = RcKaon.gRcPt*cosh(RcKaon.gRcEta);
  // float pt = RcKaon.McPt;
  // float momentum = RcKaon.McPt*cosh(RcKaon.McEta);
  // // if(!(pt > vmsa::mGlobPtMin && momentum < vmsa::mPrimMomMax)) return false; // primary pt and momentum cut: PtMin = 0.1
  // if(RcKaon.gRcDca > vmsa::mDcaTrMax_phi) return false; // dca cut: 2.0 => only use dca from global RcTracks
  // if(RcKaon.gRcNcom <= 10) return false;

  // cut for primary tracks
  if(RcKaon.pRcPt < -900.) return false; // no matched track
  if(RcKaon.pRcNfit < vmsa::mHitsFitTPCMin) return false; // nHitsFit cut
  // if(RcKaon.pRcNdedx < vmsa::mHitsDedxMin) return false;
  if(RcKaon.pRcNmax <= vmsa::mHitsMaxTPCMin) return false; // nHitsRatio cut
  if((float)RcKaon.pRcNfit/(float)RcKaon.pRcNmax < vmsa::mHitsRatioTPCMin) return false;
  // if(fabs(RcKaon.pRcEta) > vmsa::mEtaMax) return false; // eta cut
  if(fabs(RcKaon.McEta) > vmsa::mEtaMax) return false; // eta cut
  // float pt = RcKaon.pRcPt;
  // float momentum = RcKaon.pRcPt*cosh(RcKaon.pRcEta);
  float pt = RcKaon.McPt;
  float momentum = RcKaon.McPt*cosh(RcKaon.McEta);
  if(!(pt > vmsa::mGlobPtMin && momentum < vmsa::mPrimMomMax)) return false; // primary pt and momentum cut: PtMin = 0.1
  if(RcKaon.gRcDca > vmsa::mDcaTrMax_phi) return false; // dca cut: 2.0 => only use dca from global RcTracks
  if(RcKaon.pRcNcom <= 10) return false;

  // Guannan's global cut
  // if(RcKaon.gRcPt < -900.) return false; // no matched track
  // if(RcKaon.gRcNfit < 20) return false; 
  // if(RcKaon.gRcNcom <= 10) return false;
  // if(fabs(RcKaon.gRcEta) > vmsa::mEtaMax) return false; 
  // if(RcKaon.gRcDca >= 1.5) return false;

  // Guannan's primary cut
  // if(RcKaon.pRcPt < -900.) return false; // no matched track
  // if(RcKaon.pRcNfit < 20) return false;
  // if(RcKaon.pRcNcom <= 10) return false;
  // if(fabs(RcKaon.pRcEta) > vmsa::mEtaMax) return false; 
  // if(RcKaon.gRcDca >= 1.5) return false;

  return true;
}
//-----------------------------------------------------------------------------------------------------
