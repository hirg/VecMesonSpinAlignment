#ifndef StEffTPC_h
#define StEffTPC_h
#include "StMessMgr.h"
#include "TString.h"
#include <map>

class TChain;
class TNtuple;
class TBranch;
class TFile;
class StRefMultCorr;
class TH1F;
class StEffCut;
class StEffHistManger;

typedef std::map<Int_t,Float_t> FloatMap;
typedef std::map<Int_t,Int_t> IntMap;

class StEffTPC
{
  public:
    StEffTPC(Int_t Energy, Long64_t StartEvent, Long64_t StopEvent, Int_t PID);
    ~StEffTPC();

    void SetInPutList(const TString inputlist);
    void SetOutPutFile(const TString outputfile);
    void SetStartEvent(Long64_t StartEvent);
    void SetStopEvent(Long64_t StopEvent);

    void Init();
    void InitMap();
    void Make();
    void Finish();
  private:
    TString mInPutList;
    TString mOutPutFile;
    Long64_t mStartEvent;
    Long64_t mStopEvent;
    TFile *mFile_OutPut;

    static Int_t mInput_flag;

    Int_t mEnergy;
    Int_t mPID;

    FloatMap mPsiEast;
    FloatMap mPsiWest;
    IntMap mEventID;

    StEffCut *mEffCut;
    StEffHistManger *mEffHistManger;

    StRefMultCorr *mRefMultCorr; 

    TChain *mChain_Event;   
    float mRunId;
    float mEventId;
    float mMcVx;
    float mMcVy;
    float mMcVz;
    float mVx;
    float mVy;
    float mVz;
    float mVzVpd;
    float mCentrality;
    float mGRefMult;
    float mRefMult;
    float mPosRefMult;
    float mNegRefMult;
    float mZdc;
    float mBbc;
    float mNMcTracks;
    float mNRTracks;
    float mMagField;
    float mT0;
    float mT1;
    float mT2;
    float mT3;
    float mT4;
    float mT5;

    TChain *mChain_Track;   
    float mPt;
    float mP;
    float mEta;
    float mY;
    float mPhi;
    float mLabel;
    float mGeantId;
    float mKpPt;
    float mKpEta;
    float mKpPhi;
    float mKpGeantId;
    float mKpStartVtxX;
    float mKpStartVtxY;
    float mKpStartVtxZ;
    float mKpStopVtxX;
    float mKpStopVtxY;
    float mKpStopVtxZ;
    float mKpRpt;
    float mKpReta;
    float mKpRphi;
    float mKpNfit;
    float mKpNmax;
    float mKpNcom;
    float mKpNdedx;
    float mKpDedx;
    float mKpNsigKP;
    float mKpNsigKM;
    float mKpDca;
    float mKpDcaXY;
    float mKpDcaZ;
    float mKmPt;
    float mKmEta;
    float mKmPhi;
    float mKmGeantId;
    float mKmStartVtxX;
    float mKmStartVtxY;
    float mKmStartVtxZ;
    float mKmStopVtxX;
    float mKmStopVtxY;
    float mKmStopVtxZ;
    float mKmRpt;
    float mKmReta;
    float mKmRphi;
    float mKmNfit;
    float mKmNmax;
    float mKmNcom;
    float mKmNdedx;
    float mKmDedx;
    float mKmNsigKP;
    float mKmNsigKM;
    float mKmDca;
    float mKmDcaXY;
    float mKmDcaZ;
    float mInvMass;
    float mRPt;
    float mREta;
    float mRY;
    float mRphi;

    ClassDef(StEffTPC,1)
};
#endif
