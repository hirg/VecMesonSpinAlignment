#include "StRoot/StZdcSmdMaker/StZdcSmdCut.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "StRoot/StPicoDstMaker/StPicoDst.h"
#include "StRoot/StPicoDstMaker/StPicoEvent.h"
#include "StRoot/StPicoDstMaker/StPicoTrack.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"
#include "StMessMgr.h"

ClassImp(StZdcSmdCut)

StRefMultCorr* StZdcSmdCut::mRefMultCorr = NULL;
//---------------------------------------------------------------------------------

StZdcSmdCut::StZdcSmdCut(Int_t energy)
{
  mEnergy = energy;
  if(!mRefMultCorr)
  {
    mRefMultCorr = CentralityMaker::instance()->getRefMultCorr();
  }
}

//---------------------------------------------------------------------------------

StZdcSmdCut::~StZdcSmdCut()
{
  /* */
}

//---------------------------------------------------------------------------------

bool StZdcSmdCut::passEventCut(StPicoDst *pico)
{
  // initialize mMatchedToF
  mMatchedToF = 0;
  mN_prim = 0;
  mN_non_prim = 0;

  StPicoEvent *event = pico->event();
  if(!event)
  {
    return kFALSE;
  }

  // initialize StRefMultCorr
  const Int_t runId = event->runId();
  const Int_t refMult = event->refMult();
  const Float_t vx = event->primaryVertex().x();
  const Float_t vy = event->primaryVertex().y();
  const Float_t vz = event->primaryVertex().z();
  const Float_t zdcX = event->ZDCx();
  const Float_t vzVpd = event->vzVpd();
  const Bool_t isBES = (event->energy() < 200.);
  mRefMultCorr->init(runId);

  // StRefMultCorr bad run cut
  if(mRefMultCorr->isBadRun(runId))
  {
    return kFALSE;
  }

  // minBias event cut
  if(!event->isMinBias())
  {
    return kFALSE;
  }

  // event vertex cut
  // vz cut
  if(fabs(vz) > vmsa::mVzMaxMap[mEnergy])
  {
    return kFALSE;
  }
  // vr cut
  if(sqrt(vx*vx+vy*vy) > vmsa::mVrMax)
  {
    return kFALSE;
  }
  // vz-vzVpd cut for 200 GeV
  if(!isBES)
  {
    if(fabs(vz-vzVpd) > vmsa::mVzVpdDiffMax)
    {
      return kFALSE;
    }
  }

  // refMult (0-80%) cut
  if(!isBES) mRefMultCorr->initEvent(refMult,vz,zdcX); // 200GeV
  if(isBES) mRefMultCorr->initEvent(refMult,vz); // BES
  if(!mRefMultCorr->isRefMultOk())
  {
    return kFALSE;
  }

  // ToF matched points cut
  Int_t nMatchedToF = 0;
  Int_t nN_prim = 0;
  Int_t nN_non_prim = 0;
  const Int_t nTracks = pico->numberOfTracks();
  for(Int_t i = 0; i < nTracks; i++)
  {
    StPicoTrack *track = (StPicoTrack*)pico->track(i);
    if(!track)
    {
      continue;
    }
    // stop loop if already have enough TOF matched points
//    if(nMatchedToF >= vmsa::mMatchedToFMin)
//    {
//      return kTRUE;
//    }
    if(track->dca() > 3) // global track
    {
      nN_non_prim++;
    }
    else
    {
      nN_prim++;
      if(track->btofMatchFlag() > 0 && track->btof() != 0 && track->btofBeta() != 0)
      {
	nMatchedToF++;
      }
    }
  }

  mMatchedToF = nMatchedToF;
  mN_prim = nN_prim;
  mN_non_prim = nN_non_prim;


  if(nMatchedToF < vmsa::mMatchedToFMin)
  {
    return kFALSE;
  }

  return kTRUE;
}

//---------------------------------------------------------------------------------

Int_t StZdcSmdCut::getMatchedToF()
{
  return mMatchedToF;
}

Int_t StZdcSmdCut::getNpirm()
{
  return mN_prim;
}

Int_t StZdcSmdCut::getNnonprim()
{
  return mN_non_prim;
}
//---------------------------------------------------------------------------------

Float_t StZdcSmdCut::getMass2(StPicoTrack *track)
{
  Float_t Mass2 = -100.0;
  Float_t Beta = track->btofBeta();
  Float_t Momentum = track->pMom().mag(); // primary momentum

  if(track->btofMatchFlag() > 0 && track->btof() != 0 && Beta != 0)
  {
    Mass2 = Momentum*Momentum*(1.0/(Beta*Beta) - 1.0);
  }

  return Mass2;
}

bool StZdcSmdCut::passSigPionCut(StPicoTrack* track, Float_t scale_nSigma_factor)
{
  Float_t nSigmaPion = track->nSigmaPion();
  if(fabs(nSigmaPion*scale_nSigma_factor) > vmsa::mNSigmaPionMax)
  {
    return kFALSE;
  }
  return kTRUE;
}

bool StZdcSmdCut::passSigKaonCut(StPicoTrack* track, Float_t scale_nSigma_factor)
{
  Float_t nSigmaKaon = track->nSigmaKaon();
  if(fabs(nSigmaKaon*scale_nSigma_factor) > vmsa::mNSigmaKaonMax)
  {
    return kFALSE;
  }
  return kTRUE;
}

bool StZdcSmdCut::passSigProntonCut(StPicoTrack* track, Float_t scale_nSigma_factor)
{
  Float_t nSigmaProton = track->nSigmaProton();
  if(fabs(nSigmaProton*scale_nSigma_factor) > vmsa::mNSigmaProtonMax)
  {
    return kFALSE;
  }
  return kTRUE;
}

//---------------------------------------------------------------------------------

bool StZdcSmdCut::passTrackBasic(StPicoTrack *track)
{
  // nHitsFit cut
  if(track->nHitsFit() < vmsa::mHitsFitTPCMin) return kFALSE;

  // nHitsRatio cut
  if(track->nHitsMax() <= vmsa::mHitsMaxTPCMin) return kFALSE;
  if((Float_t)track->nHitsFit()/(Float_t)track->nHitsMax() < vmsa::mHitsRatioTPCMin)
    return kFALSE;

  // eta cut
  Float_t eta = track->pMom().pseudoRapidity();
  if(fabs(eta) > vmsa::mEtaMax) return kFALSE;

  return kTRUE;
}

//---------------------------------------------------------------------------------

bool StZdcSmdCut::passTrackPhi(StPicoTrack *track)
{
  if(!track) return kFALSE;

  if(!passTrackBasic(track)) return kFALSE;

  // dca cut for flow analysis: 2.0
  if(track->dca() > vmsa::mDcaTrMax_phi) return kFALSE;

  // primary pt and momentum cut: PtMin = 0.1
  if(!(track->pMom().perp() > vmsa::mGlobPtMin && track->pMom().mag() < vmsa::mPrimMomMax))
    return kFALSE;

  return kTRUE;
}

//---------------------------------------------------------------------------------

bool StZdcSmdCut::passTrackV1(StPicoTrack *track)
{
  if(!track) return kFALSE;

  if(!passTrackBasic(track)) return kFALSE;

  // dca cut for flow analysis: 2.0
  if(track->dca() > vmsa::mDcaTrMax_phi) return kFALSE;

  // primary pt and momentum cut: PtMin = 0.1
  if(!(track->pMom().perp() > vmsa::mGlobPtMin && track->pMom().mag() < vmsa::mPrimMomMax))
    return kFALSE;

  return kTRUE;
}

//---------------------------------------------------------------------------------
