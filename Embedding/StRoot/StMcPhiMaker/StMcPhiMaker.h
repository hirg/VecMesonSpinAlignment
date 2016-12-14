#ifndef StMcPhiMaker_h
#define StMcPhiMaker_h

#include <vector>
#include "TString.h"

#include "StChain/StMaker.h"
#include "StThreeVectorD.hh"
#include "StLorentzVectorF.hh"

class TFile;
class TNtuple;
class TH3F;

class StMcTrack;
class StTrack;
class StGlobalTrack;
class StAssociationMaker;
class StMcEvent;
class StEvent;
class StMuDst;
class StDedxPidTraits;
class StRefMultCorr;
class StPrimaryTrack;

class StMcPhiMaker : public StMaker
{
private:
   TString mOutfileName;
   StRefMultCorr* mRefMultCorrUtil;
   StMuDst*       mMuDst;
   std::vector<float> firedTriggersIndices;
   double mField; //.. magnetic field
   int    mCentrality;
   bool mFillTpcHitsNtuple;

   TFile* mFile;
   TNtuple* mTracks;
   TNtuple* mEventCount; //.. For counting purposes

   StMcEvent* mMcEvent;
   StEvent* mEvent;
   StAssociationMaker* mAssoc;

   StTrack const* findPartner(StMcTrack*, int&) const;
   StMcTrack const* findPartner(StGlobalTrack*, int&) const;
   StDedxPidTraits const* findDedxPidTraits(StTrack const*) const;
   int getNHitsDedx(StTrack const*) const;

   bool passTrigger();
   int  fillEventCounts(float nRTracks = -1, float nMcTracks = -1);
   int  fillTracks(int& nRTracks, int& nMcTracks);
   void fillMcTrack(float* array,int& idx,StMcTrack const*);
   void fillRcTrack(float* array,int& idx,StTrack const* gTrk,StPrimaryTrack const* pTrk,int const ncom);
   void getDca(StTrack const*,float& dca, float& dcaXY, float& dcaZ) const;

   bool isGoodMcPhi(StMcTrack const*) const;
   bool isGoodMcKPlus(StMcTrack const*) const;
   bool isGoodMcKMinus(StMcTrack const*) const;

   bool getDaughters(StMcTrack* pTrk,StMcTrack** KPTrk,StMcTrack** KMTrk);

   StLorentzVectorF* getRecoPhi(StPrimaryTrack* KPTrk,StPrimaryTrack* KMTrk);

public:
   StMcPhiMaker (const char *name="StMcPhiMaker", const char *title="event/StMcPhiMaker");

   int Init();
   int Make();
   int Finish();

   void setOutFileName(std::string);
   void fillTpcHitsNtuple(bool t=false);
   void setRefMultCorr(StRefMultCorr*);

   ClassDef(StMcPhiMaker, 0)
};

inline void StMcPhiMaker::setOutFileName(std::string s){ mOutfileName = s.c_str();}
inline void StMcPhiMaker::setRefMultCorr(StRefMultCorr* refmultCorrUtil){ mRefMultCorrUtil = refmultCorrUtil;}
#endif
