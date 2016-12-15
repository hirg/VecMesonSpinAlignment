#include "StTriFlowCut.h"
#include "StTriFlowConstants.h"
#include "StRoot/StPicoDstMaker/StPicoDst.h"
#include "StRoot/StPicoDstMaker/StPicoEvent.h"
#include "StRoot/StPicoDstMaker/StPicoTrack.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"
#include "StMessMgr.h"

ClassImp(StTriFlowCut)

StRefMultCorr* StTriFlowCut::mRefMultCorr = NULL;
//---------------------------------------------------------------------------------

StTriFlowCut::StTriFlowCut(Int_t energy)
{
  mEnergy = energy;
  if(!mRefMultCorr)
  {
    mRefMultCorr = CentralityMaker::instance()->getRefMultCorr();
  }
}

//---------------------------------------------------------------------------------

StTriFlowCut::~StTriFlowCut()
{
  /* */
}

//---------------------------------------------------------------------------------

bool StTriFlowCut::passEventCut(StPicoDst *pico)
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
//  mRefMultCorr->print();
//  cout << "StTriFLowV0:" << endl;
//  cout << mRefMultCorr->getWeight() << endl;

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
  if(fabs(vz) > TriFlow::mVzMaxMap[event->energy()])
  {
    return kFALSE;
  }
  // vr cut
  if(sqrt(vx*vx+vy*vy) > TriFlow::mVrMax)
  {
    return kFALSE;
  }
  // vz-vzVpd cut for 200 GeV
  if(!isBES)
  {
    if(fabs(vz-vzVpd) > TriFlow::mVzVpdDiffMax)
    {
      return kFALSE;
    }
  }

  // refMult (0-80%) cut
  if(!isBES) mRefMultCorr->initEvent(refMult,vz,zdcX); // 200GeV
  if(isBES) mRefMultCorr->initEvent(refMult,vz); // BES
//  cout << "eventId from StTriFlowCut: " << event->eventId()  << endl;
//  cout << "Centrality from StTriFlowCut: " << mRefMultCorr->getCentralityBin9() << endl;
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
//    if(nMatchedToF >= TriFlow::mMatchedToFMin)
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


  if(nMatchedToF < TriFlow::mMatchedToFMin)
  {
    return kFALSE;
  }

  return kTRUE;
}

//---------------------------------------------------------------------------------

Int_t StTriFlowCut::getMatchedToF()
{
  return mMatchedToF;
}

Int_t StTriFlowCut::getNpirm()
{
  return mN_prim;
}

Int_t StTriFlowCut::getNnonprim()
{
  return mN_non_prim;
}
//---------------------------------------------------------------------------------

Float_t StTriFlowCut::getMass2(StPicoTrack *track)
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

Float_t StTriFlowCut::getV0Mass2(StPicoTrack *track)
{
  Float_t Mass2 = -100.0;
  Float_t Beta = track->btofBeta();
  Float_t Momentum = track->gMom().mag(); // global momentum

  if(track->btofMatchFlag() > 0 && track->btof() != 0 && Beta != 0)
  {
    Mass2 = Momentum*Momentum*(1.0/(Beta*Beta) - 1.0);
  }

  return Mass2;
}

bool StTriFlowCut::passPIDCut(StPicoTrack *track)
{
  // mass2 cut
  Float_t Mass2 = getMass2(track);
  if(Mass2 < TriFlow::mMass2Min)
  {
    return kFALSE;
  }

  // nHitsDedx cut
  Int_t nHitsDedx = track->nHitsDedx();
  if(nHitsDedx < TriFlow::mHitsDedxMin)
  {
    return kFALSE;
  }

  // ToFYLocal cut
  Float_t ToFYLocal = track->btofYLocal();
  if(fabs(ToFYLocal) > TriFlow::mToFYLocalMax)
  {
    return kFALSE;
  }

  return kTRUE;
}

bool StTriFlowCut::passSigElectronCut(StPicoTrack* track, Float_t scale_nSigma_factor)
{
  Float_t nSigmaElectron = track->nSigmaElectron();
  if(fabs(nSigmaElectron*scale_nSigma_factor) > TriFlow::mNSigmaElectronMax)
  {
    return kFALSE;
  }
  return kTRUE;
}

bool StTriFlowCut::passSigPionCut(StPicoTrack* track, Float_t scale_nSigma_factor)
{
  Float_t nSigmaPion = track->nSigmaPion();
  if(fabs(nSigmaPion*scale_nSigma_factor) > TriFlow::mNSigmaPionMax)
  {
    return kFALSE;
  }
  return kTRUE;
}

bool StTriFlowCut::passSigKaonCut(StPicoTrack* track, Float_t scale_nSigma_factor)
{
  Float_t nSigmaKaon = track->nSigmaKaon();
  if(fabs(nSigmaKaon*scale_nSigma_factor) > TriFlow::mNSigmaKaonMax)
  {
    return kFALSE;
  }
  return kTRUE;
}

bool StTriFlowCut::passSigProntonCut(StPicoTrack* track, Float_t scale_nSigma_factor)
{
  Float_t nSigmaProton = track->nSigmaProton();
  if(fabs(nSigmaProton*scale_nSigma_factor) > TriFlow::mNSigmaProtonMax)
  {
    return kFALSE;
  }
  return kTRUE;
}

bool StTriFlowCut::passSigProntonCutSys(StPicoTrack* track, Float_t scale_nSigma_factor, Int_t i_proton)
{
  Float_t nSigmaProton = track->nSigmaProton();
  if(fabs(nSigmaProton*scale_nSigma_factor) > TriFlow::mNSigmaProtonMaxSys[i_proton])
  {
    return kFALSE;
  }
  return kTRUE;
}
//---------------------------------------------------------------------------------

bool StTriFlowCut::passTrackBasic(StPicoTrack *track)
{
  // nHitsFit cut
  if(track->nHitsFit() < TriFlow::mHitsFitTPCMin)
  {
    return kFALSE;
  }

  // nHitsRatio cut
  if(track->nHitsMax() <= TriFlow::mHitsMaxTPCMin)
  {
    return kFALSE;
  }
  if((Float_t)track->nHitsFit()/(Float_t)track->nHitsMax() < TriFlow::mHitsRatioTPCMin)
  {
    return kFALSE;
  }

  // eta cut
  Float_t eta = track->pMom().pseudoRapidity();
  if(fabs(eta) > TriFlow::mEtaMax)
  {
    return kFALSE;
  }

  return kTRUE;
}

bool StTriFlowCut::passTrackEP(StPicoTrack *track)
{
  if(!track) return kFALSE;

  if(!passTrackBasic(track)) return kFALSE;

  // dca cut for event plane reconstruction: 200GeV = 3.0, BES = 1.0
  if(track->dca() > TriFlow::mDcaEPMax[mEnergy])
  {
    return kFALSE;
  }

  // pt cut 0.2 - 2.0 GeV/c
  Float_t pt = track->pMom().perp();
  Float_t p  = track->pMom().mag();
  if(!(pt > TriFlow::mPrimPtMin[mEnergy] && pt < TriFlow::mPrimPtMax && p < TriFlow::mPrimMomMax))
  {
    return kFALSE;
  }

  return kTRUE;
}

bool StTriFlowCut::passTrackCut(StPicoTrack *track)
{
  if(!track) return kFALSE;

  if(!passTrackBasic(track)) return kFALSE;

  // dca cut for flow analysis: 1.0, 1.5 and 2.0
  if(track->dca() > TriFlow::mDcaTrMax)
  {
    return kFALSE;
  }

  // primary pt and momentum cut: PtMin = 0.15 for 200 GeV, PtMin = 0.2 for BES
  if(!(track->pMom().perp() > TriFlow::mPrimPtMin[mEnergy] && track->pMom().mag() < TriFlow::mPrimMomMax))
  {
    return kFALSE;
  }

  return kTRUE;
}

bool StTriFlowCut::passTrackCutSys(StPicoTrack *track, Int_t i_dca, Int_t i_nHitsFit)
{
  if(!track) return kFALSE;

  if(!passTrackBasic(track)) return kFALSE;

  // dca cut for flow analysis: 1.0, 1.5 and 2.0
  if(track->dca() > TriFlow::mDcaTrMaxSys[i_dca])
  {
    return kFALSE;
  }

  // nHitsFit cut: 15, 20
  if(track->nHitsFit() < TriFlow::mHitsFitTPCMinSys[i_nHitsFit])
  {
    return kFALSE;
  }

  // primary pt and momentum cut: PtMin = 0.15 for 200 GeV, PtMin = 0.2 for BES
  if(!(track->pMom().perp() > TriFlow::mPrimPtMin[mEnergy] && track->pMom().mag() < TriFlow::mPrimMomMax))
  {
    return kFALSE;
  }

  return kTRUE;
}
//---------------------------------------------------------------------------------
bool StTriFlowCut::passTrackPhi(StPicoTrack *track)
{
  if(!track) return kFALSE;

  if(!passTrackBasic(track)) return kFALSE;

  // dca cut for flow analysis: 2.0
  if(track->dca() > TriFlow::mDcaTrMax_phi)
  {
    return kFALSE;
  }

  // primary pt and momentum cut: PtMin = 0.1
  if(!(track->pMom().perp() > TriFlow::mGlobPtMin && track->pMom().mag() < TriFlow::mPrimMomMax))
  {
    return kFALSE;
  }

  return kTRUE;
}

bool StTriFlowCut::passTrackV0(StPicoTrack *track)
{
  if(!track) return kFALSE;

  // nHitsFit cut
  if(track->nHitsFit() < TriFlow::mHitsFitTPCMin)
  {
    return kFALSE;
  }

  // nHitsRatio cut
  if(track->nHitsMax() <= TriFlow::mHitsMaxTPCMin)
  {
    return kFALSE;
  }
  if((Float_t)track->nHitsFit()/(Float_t)track->nHitsMax() < TriFlow::mHitsRatioTPCMin)
  {
    return kFALSE;
  }

  // global pt and momentum cut: PtMin = 0.1
  if(!(track->gMom().perp() > TriFlow::mGlobPtMin && track->gMom().mag() < TriFlow::mPrimMomMax))
  {
    return kFALSE;
  }

  /*
  // eta cut
  Float_t eta = track->gMom().pseudoRapidity();
  if(fabs(eta) > TriFlow::mEtaMax)
  {
    return kFALSE;
  }
  */

  return kTRUE;
}
