#include "TFile.h"
#include "TH3F.h"
#include "TNtuple.h"
#include "TSystem.h"

#include "StarClassLibrary/StParticleDefinition.hh"
#include "StarClassLibrary/SystemOfUnits.h"
#include "StarClassLibrary/StParticleTypes.hh"
#include "StBFChain/StBFChain.h"

#include "StEvent/StEventTypes.h"
#include "StEvent/StTrack.h"
#include "StEvent/StTrackGeometry.h"
#include "StEvent/StTrackNode.h"
#include "StEvent/StGlobalTrack.h"
#include "StEvent/StTrackTopologyMap.h"
#include "StEvent/StEventSummary.h"
#include "StEvent/StBTofCollection.h"
#include "StEvent/StBTofHeader.h"
#include "StEvent/StEnumerations.h"
#include "StEvent/StEvent.h"
#include "StEvent/StTpcDedxPidAlgorithm.h"
#include "StEventUtilities/StuRefMult.hh"
#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StMuDSTMaker/COMMON/StMuPrimaryVertex.h"

#include "StMcEvent/StMcEventTypes.hh"
#include "StMcEvent/StMcContainers.hh"
#include "StMcEvent/StMcVertex.hh"
#include "StMcEvent/StMcEvent.hh"

#include "StAssociationMaker/StAssociationMaker.h"
#include "StAssociationMaker/StTrackPairInfo.hh"

#include "StRefMultCorr/StRefMultCorr.h"
#include "StMcPhiCuts.h"
#include "StMcPhiMaker.h"

ClassImp(StMcPhiMaker);

StMcPhiMaker::StMcPhiMaker(const char *name, const char *title): StMaker(name),
   mRefMultCorrUtil(NULL), mMuDst(NULL), mField(-999), mCentrality(-999), mFillTpcHitsNtuple(false), 
   mFile(NULL), mTracks(NULL), mEventCount(NULL), mMcEvent(NULL), mEvent(NULL), mAssoc(NULL)
{
   LOG_INFO << "StMcPhiMaker() - DONE" << endm;
}
//__________________________________
int StMcPhiMaker::Init()
{
   if (!mOutfileName.Length())
   {
      // StBFChain* bfChain = (StBFChain *) StMaker::GetChain();
      //
      // if (!bfChain) return kStFatal;
      //
      // mOutfileName = bfChain->GetFileIn();

      if (mOutfileName.Length())
      {
         LOG_INFO << mOutfileName << endm;
         mOutfileName = gSystem->BaseName(mOutfileName.Data());
         mOutfileName = mOutfileName.ReplaceAll(".event.root", "");
         mOutfileName = mOutfileName.ReplaceAll(".geant.root", "");
         mOutfileName = mOutfileName.ReplaceAll(".MuDst.root", "");
      }
      else
      {
         mOutfileName = "mcAnalysis";
      }
   }

   mOutfileName = mOutfileName.ReplaceAll(".root", "");
   mFile = new TFile(Form("%s.McAna.root", mOutfileName.Data()), "recreate");
   assert(mFile && !mFile->IsZombie());

   mAssoc = (StAssociationMaker*)GetMaker("StAssociationMaker");
   if (!mAssoc)
   {
      LOG_ERROR << "Could not get StAssociationMaker" << endm;
      return kStErr;
   }

   for (int ii = 0; ii < McPhiCuts::maxNumberOfTriggers; ++ii)
   {
      firedTriggersIndices.push_back(-1);
   };

   const char* evtVarList = "runId:eventId:mcVx:mcVy:mcVz:vx:vy:vz:vzVpd:"
     "centrality:gRefMult:RefMult:posRefMult:negRefMult:zdc:bbc:nMcTracks:nRTracks:magField:t0:t1:t2:t3:t4:t5";

   const char* varlist = "pt:p:eta:y:phi:label:geantId:" // MC phi 
     "kpPt:kpEta:kpPhi:kpGeantId:kpStartVtxX:kpStartVtxY:kpStartVtxZ:kpStopVtxX:kpStopVtxY:kpStopVtxZ:" // MC K+ 
     "kpRpt:kpReta:kpRphi:kpNfit:kpNmax:kpNcom:kpNdedx:kpDedx:kpNsigKP:kpNsigKM:kpDca:kpDcaXY:kpDcaZ:" // primary RC K+
     "kmPt:kmEta:kmPhi:kmGeantId:kmStartVtxX:kmStartVtxY:kmStartVtxZ:kmStopVtxX:kmStopVtxY:kmStopVtxZ:" // MC K-
     "kmRpt:kmReta:kmRphi:kmNfit:kmNmax:kmNcom:kmNdedx:kmDedx:kmNsigKP:kmNsigKM:kmDca:kmDcaXY:kmDcaZ:" // primary RC K+
     "invMass:rPt:rEta:rY:rphi"; // reconstructed phi

   mEventCount = new TNtuple("eventCount","eventCount",evtVarList);

   mTracks = new TNtuple("vecMeson","",varlist);

   LOG_INFO << "Init() - DONE" << endm;

   return kStOk;
}

//__________________________________
int StMcPhiMaker::Make()
{
   StMuDstMaker* muDstMaker = (StMuDstMaker*)GetMaker("MuDst");

   if (!muDstMaker)
   {
      LOG_WARN << " No MuDstMaker, will try to take all event information from StEvent" << endm;
      mMuDst = NULL;
   }
   else
   {
     mMuDst = (StMuDst*)muDstMaker->muDst();
   }

   if(!mMuDst || !mMuDst->event())
   {
     LOG_WARN << "MuDst or mMuDst->event() is missing, will try to take all event information from StEvent" << endm;
     mMuDst = NULL;
   }

   mMcEvent = (StMcEvent*)GetDataSet("StMcEvent");

   if (!mMcEvent)
   {
      LOG_WARN << "No StMcEvent" << endm;
      return kStWarn;
   }

   mEvent = (StEvent*)GetDataSet("StEvent");
   if (!mEvent)
   {
      LOG_WARN << "No StEvent" << endm;
      return kStWarn;
   }

   mField = (float)mEvent->summary()->magneticField();

   if(mRefMultCorrUtil && mMuDst)
   {
     mRefMultCorrUtil->init(mEvent->runId());

     mRefMultCorrUtil->initEvent(mMuDst->event()->refMult(),
                                  mEvent->primaryVertex()->position().z(), 
                                  mEvent->runInfo()->zdcCoincidenceRate());

     mCentrality  = mRefMultCorrUtil->getCentralityBin9();

     if (mRefMultCorrUtil->isBadRun(mEvent->runId()))
     {
       LOG_INFO << "This is a bad run from mRefMultCorrUtil! Skip! " << endm;

       return kStSkip;
     }
   }
   else
   {
     mCentrality = -999;
   }

   // Fill
   int nRTracks = -1;
   int nMcTracks = -1;

   int fillTracksStatus = kStOk;
   // fillTracksStatus = fillTracks(nRTracks, nMcTracks);
   if (passTrigger())
   {
     fillTracksStatus = fillTracks(nRTracks, nMcTracks);
   }
   else
   {
      LOG_INFO << "No interesting triggers. Counting event then skipping." << endm;
   }

   int fillEventCountStatus = fillEventCounts((float)nRTracks, (float)nMcTracks);

   return fillTracksStatus && fillEventCountStatus;
}

int StMcPhiMaker::fillTracks(int& nRTracks, int& nMcTracks)
{
   nRTracks = 0;
   nMcTracks = 0;

   LOG_INFO << "Filling " << mMcEvent->tracks().size() << " mcTracks..." << "\n";

   for (unsigned int iTrk = 0;  iTrk < mMcEvent->tracks().size(); ++iTrk)
   {
      StMcTrack* const mcTrack = mMcEvent->tracks()[iTrk];

      if (!mcTrack)
      {
         LOG_WARN << "Empty mcTrack container" << endm;
         continue;
      }

      if (!isGoodMcPhi(mcTrack)) continue; // MC phi meson

      float pTrkSvx, pTrkSvy;
      if(mcTrack->startVertex()) 
      {
	pTrkSvx = mcTrack->startVertex()->position().x();
	pTrkSvy = mcTrack->startVertex()->position().y();
      }
      if(sqrt(pTrkSvx*pTrkSvx + pTrkSvy*pTrkSvy)>3) continue; // dca < 3

      StMcTrack* phi_mcTrk = 0;
      StMcTrack* KP_mcTrk = 0; // D0 daughter
      StMcTrack* KM_mcTrk = 0;
      if(!getDaughters(mcTrack,&KP_mcTrk,&KM_mcTrk)) continue;
      phi_mcTrk = mcTrack;
      ++nMcTracks;

      // get reconstructed tracks
      int KPNcom = -9999;
      int KMNcom = -9999;
      const StTrack* KP_gTrk = findPartner(KP_mcTrk, KPNcom);
      const StTrack* KM_gTrk = findPartner(KM_mcTrk, KMNcom);
      StPrimaryTrack* KP_pTrk = 0;
      StPrimaryTrack* KM_pTrk = 0;
      if(KP_gTrk) KP_pTrk = (StPrimaryTrack*)KP_gTrk->node()->track(primary); // primary
      if(KM_gTrk) KM_pTrk = (StPrimaryTrack*)KM_gTrk->node()->track(primary);

      // if pi_pTrk and k_pTrk exists, reconstruct D0
      StLorentzVectorF* recoPhi = 0;
      if(KP_pTrk && KM_pTrk) 
      {
	recoPhi = getRecoPhi(KP_pTrk,KM_pTrk);
	++nRTracks;
      }

      float array[220];
      for (int ii = 0; ii < 220; ++ii) array[ii] = -999;

      int idx = 0;
      
      // MC phi
      array[idx++] = phi_mcTrk->pt();
      array[idx++] = phi_mcTrk->momentum().mag();
      array[idx++] = phi_mcTrk->pseudoRapidity();
      array[idx++] = phi_mcTrk->rapidity();
      array[idx++] = phi_mcTrk->momentum().phi();
      array[idx++] = phi_mcTrk->eventGenLabel();
      array[idx++] = phi_mcTrk->geantId();
      fillMcTrack(array,idx,KP_mcTrk); // K+ mcTrack
      fillRcTrack(array,idx,KP_gTrk,KP_pTrk,KPNcom); // K+ primary rcTrack
      fillMcTrack(array,idx,KM_mcTrk); // K- mcTrack
      fillRcTrack(array,idx,KM_gTrk,KM_pTrk,KMNcom); // K- primary rcTrack
      array[idx++] = recoPhi ? recoPhi->m() : -999.;
      array[idx++] = recoPhi ? recoPhi->perp() : -999.;
      array[idx++] = recoPhi ? recoPhi->pseudoRapidity() : -999.;
      array[idx++] = recoPhi ? recoPhi->rapidity() : -999.;
      array[idx++] = recoPhi ? recoPhi->phi() : -999.;

      mTracks->Fill(array);
   }

   LOG_INFO << endm;

   return kStOk;
}

void StMcPhiMaker::fillMcTrack(float* array, int& idx, StMcTrack const* const mcTrk)
{
  array[idx++] = mcTrk->pt();
  array[idx++] = mcTrk->pseudoRapidity();
  array[idx++] = mcTrk->momentum().phi();
  array[idx++] = mcTrk->geantId();
  array[idx++] = mcTrk->startVertex()->position().x();
  array[idx++] = mcTrk->startVertex()->position().y();
  array[idx++] = mcTrk->startVertex()->position().z();
  if (mcTrk->stopVertex())
  {
    array[idx++] = mcTrk->stopVertex()->position().x();
    array[idx++] = mcTrk->stopVertex()->position().y();
    array[idx++] = mcTrk->stopVertex()->position().z();
  }
  else
  {
    idx += 3;
  }
}

void StMcPhiMaker::fillRcTrack(float* array, int& idx, StTrack const* const gTrk, StPrimaryTrack const* const rcTrack, int const ncom)
{
  array[idx++] = rcTrack ? rcTrack->geometry()->momentum().perp() : -999.;
  array[idx++] = rcTrack ? rcTrack->geometry()->momentum().pseudoRapidity() : -999.;
  array[idx++] = rcTrack ? rcTrack->geometry()->momentum().phi() : -999.;
  array[idx++] = rcTrack ? rcTrack->fitTraits().numberOfFitPoints(kTpcId) : -999.;
  array[idx++] = rcTrack ? rcTrack->numberOfPossiblePoints(kTpcId) : -999.;
  array[idx++] = rcTrack ? ncom : -999.;

  // dedx info
  float nDedxPts = -999.;
  float dedx = -999.;
  float nSigKP = -999.;
  float nSigKM = -999.;
  static StTpcDedxPidAlgorithm aplus(McPhiCuts::dedxMethod);
  static StKaonPlus* KPlus = StKaonPlus::instance();
  static StKaonMinus* KMinus = StKaonMinus::instance();
  StParticleDefinition const* prtcl = 0;
  if(gTrk) prtcl = gTrk->pidTraits(aplus);
  if (prtcl)
  {
    nDedxPts = aplus.traits()->numberOfPoints();
    dedx = aplus.traits()->mean();
    nSigKP = aplus.numberOfSigma(KPlus);
    nSigKM = aplus.numberOfSigma(KMinus);
  }

  // array[idx++] = getNHitsDedx(gTrk);
  array[idx++] = nDedxPts;
  array[idx++] = dedx;
  array[idx++] = nSigKP;
  array[idx++] = nSigKM;

  float dca = -999.;
  float dcaXY = -999.;
  float dcaZ = -999.;

  if(gTrk) getDca(gTrk, dca, dcaXY, dcaZ);

  array[idx++] = dca;
  array[idx++] = dcaXY;
  array[idx++] = dcaZ;
}

bool StMcPhiMaker::isGoodMcPhi(StMcTrack const* const mcTrack) const
{
   // return mcTrack->geantId() == McPhiCuts::geantId && mcTrack->startVertex()->position().perp() < McPhiCuts::mcTrackStartVtxR;
   return mcTrack->geantId() == McPhiCuts::geantId_phi;
}

bool StMcPhiMaker::isGoodMcKPlus(StMcTrack const* const mcTrack) const
{
   return mcTrack->geantId() == McPhiCuts::geantId_KP;
}

bool StMcPhiMaker::isGoodMcKMinus(StMcTrack const* const mcTrack) const
{
   return mcTrack->geantId() == McPhiCuts::geantId_KM;
}

int StMcPhiMaker::fillEventCounts(float nRTracks, float nMcTracks)
{
   float vars[50];

   float vpdVz = -999;
   StBTofHeader* tofheader = 0;
   if (mEvent->btofCollection())  tofheader = mEvent->btofCollection()->tofHeader();
   if (tofheader) vpdVz = tofheader->vpdVz();

   int iVar = 0;
   vars[iVar++] = (float)mEvent->runId();
   vars[iVar++] = (float)mEvent->id();
   vars[iVar++] = (float)mMcEvent->primaryVertex()->position().x();
   vars[iVar++] = (float)mMcEvent->primaryVertex()->position().y();
   vars[iVar++] = (float)mMcEvent->primaryVertex()->position().z();
   vars[iVar++] = (float)mEvent->primaryVertex()->position().x();
   vars[iVar++] = (float)mEvent->primaryVertex()->position().y();
   vars[iVar++] = (float)mEvent->primaryVertex()->position().z();
   vars[iVar++] = vpdVz;
   vars[iVar++] = mCentrality;
   vars[iVar++] = mMuDst? mMuDst->event()->grefmult() : -999;
   vars[iVar++] = mMuDst? mMuDst->event()->refMult() : -999;
   vars[iVar++] = (float)uncorrectedNumberOfPositivePrimaries(*mEvent);
   vars[iVar++] = (float)uncorrectedNumberOfNegativePrimaries(*mEvent);
   vars[iVar++] = (float)mEvent->runInfo()->zdcCoincidenceRate();
   vars[iVar++] = (float)mEvent->runInfo()->bbcCoincidenceRate();
   vars[iVar++] = nMcTracks;
   vars[iVar++] = nRTracks;
   vars[iVar++] = (float)mEvent->summary()->magneticField() / 10;
   vars[iVar++] = firedTriggersIndices.at(0);
   vars[iVar++] = firedTriggersIndices.at(1);
   vars[iVar++] = firedTriggersIndices.at(2);
   vars[iVar++] = firedTriggersIndices.at(3);
   vars[iVar++] = firedTriggersIndices.at(4);
   vars[iVar++] = firedTriggersIndices.at(5);

   mEventCount->Fill(vars);

   return kStOk;
}

bool StMcPhiMaker::passTrigger()
{
   LOG_INFO << "Checking triggers..." << endm;
   bool interesting_event = false;

   if (!mEvent)
   {
      LOG_FATAL << "mEvent doesn't exist" << endm;
   }

   if (McPhiCuts::interesting_triggers.size() == 0)
   {
      LOG_WARN << "No triggers in McPhiCuts::interesting_triggers ... accepting event anyway" << endm;
      return true;
   }

   const StTriggerId* st_trgid = mEvent->triggerIdCollection()->nominal();

   for (unsigned int ii = 0; ii < firedTriggersIndices.size(); ++ii)
   {
      firedTriggersIndices[ii] = -1;
   }

   // Fill interesting triggers
   LOG_INFO << "Interesting fired triggers: " << "\n";

   for (unsigned int ii = 0; ii < st_trgid->maxTriggerIds(); ++ii)
   {
      unsigned int id = st_trgid->triggerId(ii);

      int trgIndex = -1;

      for (unsigned int jj = 0; jj < McPhiCuts::interesting_triggers.size(); ++jj)
      {
         if (McPhiCuts::interesting_triggers[jj] == id)
         {
            trgIndex = jj;
            interesting_event = true;
            LOG_INFO << id << " ";
            break;
         }
      }

      if (trgIndex != -1) firedTriggersIndices[trgIndex] = 1.0;
   }

   LOG_INFO << endm;

   return interesting_event;
}

StTrack const* StMcPhiMaker::findPartner(StMcTrack* mcTrack, int& maxCommonTpcHits) const
{
   //..StMcTrack find partner from the StTracks
   pair<mcTrackMapIter, mcTrackMapIter> p = mAssoc->mcTrackMap()->equal_range(mcTrack);

   const StTrack* maxTrack = 0;
   maxCommonTpcHits = 0;
   for (mcTrackMapIter k = p.first; k != p.second; ++k)
   {
      int commonTpcHits = k->second->commonTpcHits();
      const StTrack* track = k->second->partnerTrack()->node()->track(global);//should be global
      if (track && commonTpcHits > maxCommonTpcHits)
      {
         maxTrack = track;
         maxCommonTpcHits = commonTpcHits;
      }
   }
   return maxTrack;
}

void StMcPhiMaker::getDca(StTrack const* const gTrk, float& dca, float& dcaXY, float& dcaZ) const
{
   StPhysicalHelixD helix = gTrk->geometry()->helix();
   dca = helix.distance(mEvent->primaryVertex()->position());
   dcaXY = helix.geometricSignedDistance(mEvent->primaryVertex()->position().x(), mEvent->primaryVertex()->position().y());

   StThreeVectorF dcaPoint = helix.at(helix.pathLength(mEvent->primaryVertex()->position()));
   dcaZ = dcaPoint.z() - mEvent->primaryVertex()->position().z();
}

StMcTrack const* StMcPhiMaker::findPartner(StGlobalTrack* rcTrack, int& maxCommonTpcHits) const
{
   //.. StGlobalTracks find partner from StMcTracks.
   //.. See example from StRoot/StMcPhiMaker
   pair<rcTrackMapIter, rcTrackMapIter> p = mAssoc->rcTrackMap()->equal_range(rcTrack);

   const StMcTrack* maxTrack = 0;
   maxCommonTpcHits = 0;
   for (rcTrackMapIter k = p.first; k != p.second; ++k)
   {
      int commonTpcHits = k->second->commonTpcHits();

      const StMcTrack* track = k->second->partnerMcTrack();

      if (track && commonTpcHits > maxCommonTpcHits)
      {
         maxTrack = track;
         maxCommonTpcHits = commonTpcHits;
      }
   }
   return maxTrack;
}

int StMcPhiMaker::Finish()
{
   mFile->cd();
   mFile->Write();
   mFile->Close();
   return kStOk;
}

StDedxPidTraits const* StMcPhiMaker::findDedxPidTraits(StTrack const* const rcTrack) const
{
   StDedxPidTraits* pid = 0;
   StPtrVecTrackPidTraits traits = rcTrack->pidTraits(kTpcId);

   for (unsigned int ii = 0; ii < traits.size(); ++ii)
   {
      pid = dynamic_cast<StDedxPidTraits*>(traits[ii]);
      if (pid && pid->method() == McPhiCuts::dedxMethod) break;
   }

   return pid;
}

int StMcPhiMaker::getNHitsDedx(StTrack const* const t) const
{
   int ndedx = -1;
   StPtrVecTrackPidTraits pidTraits = t->pidTraits(kTpcId);

   if (pidTraits.size())
   {
      StDedxPidTraits* pid;
      for (unsigned int ii = 0; ii < pidTraits.size(); ++ii)
      {
         pid = dynamic_cast<StDedxPidTraits*>(pidTraits[ii]);

         if (pid && (pid->method() == McPhiCuts::dedxMethod))
         {
            ndedx = pid->numberOfPoints();            //number of dedx hits
            break;
         }
      }
   }

   return ndedx;
}

bool StMcPhiMaker::getDaughters(StMcTrack* pTrk,StMcTrack** KPTrk,StMcTrack** KMTrk)
{
  if(!pTrk->stopVertex()) return false;

  int ndaughter = pTrk->stopVertex()->numberOfDaughters();

  if(ndaughter!=2) return false;

  for(int i=0;i<ndaughter;i++)
  {
    StMcTrack* trk = pTrk->stopVertex()->daughter(i);

    // phi
    if( isGoodMcPhi(pTrk) && isGoodMcKPlus(trk) ) *KPTrk=trk;
    if( isGoodMcPhi(pTrk) && isGoodMcKMinus(trk) ) *KMTrk=trk;
  }

  if(!*KPTrk || !*KMTrk) return false;

  return true;
}

StLorentzVectorF* StMcPhiMaker::getRecoPhi(StPrimaryTrack* KPTrk,StPrimaryTrack* KMTrk)
{
  if (!KPTrk || !KMTrk) return 0;

  StLorentzVectorF KPFourMom(KPTrk->geometry()->momentum(),KPTrk->geometry()->momentum().massHypothesis(.493677));
  StLorentzVectorF KMFourMom(KMTrk->geometry()->momentum(),KMTrk->geometry()->momentum().massHypothesis(.493677));

  StLorentzVectorF* phiFourMom = new StLorentzVectorF();
  *phiFourMom = KPFourMom + KMFourMom;

  return phiFourMom;
}
