#include "StEffCut.h"
#include <math.h>

ClassImp(StEffCut)

StEffCut::StEffCut(Int_t energy)
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
  float vx = EventHeader.mcVx;
  float vy = EventHeader.mcVy;
  float vz = EventHeader.mcVz;
  if(fabs(vz) > Efficiency::mVzMaxMap[mEnergy]) // Vz cut
    return kFALSE;
  if(sqrt(vx*vx+vy*vy) > Efficiency::mVrMax) //Vr cut
    return kFALSE;
  if(mEnergy == 6)
  {
    float vzVpd = EventHeader.vzVpd;
    if(fabs(vz-vzVpd) > Efficiency::mVzVpdDiffMax)
      return kFALSE;
  }

  // centrality cut
  // if(EventHeader.centrality < Efficiency::cent_low[0] || EventHeader.centrality > Efficiency::cent_high[0]) // 20%-60%
  if(EventHeader.centrality < 0) return kFALSE;

  return kTRUE;
}
//-----------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------
bool StEffCut::passTrackCutPhi(McVecMeson McPhi)
{
  // eta cut
  if(fabs(McPhi.eta) > Efficiency::mEtaMax) return kFALSE;

  return kTRUE;
}

bool StEffCut::passTrackCutPhi(RcVecMeson RcPhi)
{
  // eta cut
  if(fabs(RcPhi.McEta) > Efficiency::mEtaMax) return kFALSE;

  return kTRUE;
}
//-----------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------
bool StEffCut::passTrackCut(McDecayDau McKaon)
{
  // eta cut for some reason
  if(fabs(McKaon.eta) > Efficiency::mEtaMax) return kFALSE;

  float pt = McKaon.pt;
  float momentum = McKaon.pt*cosh(McKaon.eta);
  if(!(pt > Efficiency::mGlobPtMin && momentum < Efficiency::mPrimMomMax)) return kFALSE;
  
  return kTRUE;
}

bool StEffCut::passTrackCut(RcDecayDau RcKaon)
{
  // no matched track
  if( RcKaon.pt < -900. ) return kFALSE; 

  if(RcKaon.nCom <= 10) return kFALSE;
  // nHitsFit cut
  if(RcKaon.nFit < Efficiency::mHitsFitTPCMin) return kFALSE;
  // if(RcKaon.nDedx < Efficiency::mHitsDedxMin) return kFALSE;

  // nHitsRatio cut
  if(RcKaon.nMax <= Efficiency::mHitsMaxTPCMin) return kFALSE;
  if((Float_t)RcKaon.nFit/(Float_t)RcKaon.nMax < Efficiency::mHitsRatioTPCMin) return kFALSE;

  // eta cut
  if(fabs(RcKaon.McEta) > Efficiency::mEtaMax) return kFALSE; // use McKaon cut

  // nSigmaKaon cut 2.5
  // if(fabs(RcKaon.nSigKP) > Efficiency::mSigKaon) return kFALSE;
  // if(fabs(RcKaon.nSigKM) > Efficiency::mSigKaon) return kFALSE;

  // dca cut for flow analysis: 2.0
  if(RcKaon.dca > Efficiency::mDcaTrMax_phi) return kFALSE;

  float pt = RcKaon.McPt;
  float momentum = RcKaon.McPt*cosh(RcKaon.McEta);
  // primary pt and momentum cut: PtMin = 0.1
  if(!(pt > Efficiency::mGlobPtMin && momentum < Efficiency::mPrimMomMax)) return kFALSE;

  return kTRUE;
}
//-----------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------
bool StEffCut::passPhiEtaEast(McDecayDau McKaon) // neg
{
  if(!(McKaon.eta > -1.0*Efficiency::mEtaMax && McKaon.eta < 0.0)) // eta cut
    return kFALSE;

  return kTRUE;
}

bool StEffCut::passPhiEtaWest(McDecayDau McKaon) // pos
{
  if(!(McKaon.eta >= 0.0 && McKaon.eta < Efficiency::mEtaMax)) // eta cut
    return kFALSE;

  return kTRUE;
}

bool StEffCut::passPhiEtaEast(RcDecayDau RcKaon) // neg
{
  if(!(RcKaon.eta > -1.0*Efficiency::mEtaMax && RcKaon.eta < 0.0)) // eta cut
    return kFALSE;

  return kTRUE;
}

bool StEffCut::passPhiEtaWest(RcDecayDau RcKaon) // pos
{
  if(!(RcKaon.eta >= 0.0 && RcKaon.eta < Efficiency::mEtaMax)) // eta cut
    return kFALSE;

  return kTRUE;
}
//-----------------------------------------------------------------------------------------------------
