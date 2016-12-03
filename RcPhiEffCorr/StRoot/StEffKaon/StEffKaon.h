#ifndef StEffKaon_h
#define StEffKaon_h
#include "StMessMgr.h"
#include <string>
#include <map>

class TChain;
class TNtuple;
class TBranch;
class TFile;
class StRefMultCorr;
class TH1F;
class StEffCut;
class StEffHistManger;

class StEffKaon
{
  public:
    StEffKaon(int Energy, long StartEvent, long StopEvent, int PID);
    ~StEffKaon();

    void SetInPutList(const string inputlist);
    void SetOutPutFile(const string outputfile);
    void SetStartEvent(long StartEvent);
    void SetStopEvent(long StopEvent);

    void Init();
    void Make();
    void Finish();
  private:
    string mInPutList;
    string mOutPutFile;
    long mStartEvent;
    long mStopEvent;
    TFile *mFile_OutPut;

    static int mInput_flag;

    int mEnergy;
    int mPID;

    StEffCut *mEffCut;
    StEffHistManger *mEffHistManger;

    StRefMultCorr *mRefMultCorr; 

    TChain *mChain_Event; // event header
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
    float mNRcTracks;
    float mMagField;
    float mT0;
    float mT1;
    float mT2;
    float mT3;
    float mT4;
    float mT5;

    TChain *mChain_Track; // McTracks, global RcTracks and primary RcTracks
    float mMcPt;
    float mMcP;
    float mMcEta;
    float mMcY;
    float mMcPhi;
    float mMcGeantId;
    float mMcEventGenLabel;
    float mMcStartVtxX;
    float mMcStartVtxY;
    float mMcStartVtxZ;
    float mMcStopVtxX;
    float mMcStopVtxY;
    float mMcStopVtxZ;

    float mGRcPt;
    float mGRcEta;
    float mGRcPhi;
    float mGRcNfit;
    float mGRcNmax;
    float mGRcNcom;
    float mGRcNdedx;
    float mGRcDedx;
    float mGRcNsigKP;
    float mGRcNsigKM;
    float mGRcDca;
    float mGRcDcaXY;
    float mGRcDcaZ;

    float mPRcPt;
    float mPRcEta;
    float mPRcPhi;
    float mPRcNfit;
    float mPRcNmax;
    float mPRcNcom;
    float mPRcNdedx;
    float mPRcDedx;
    float mPRcNsigKP;
    float mPRcNsigKM;
    float mPRcDca;
    float mPRcDcaXY;
    float mPRcDcaZ;

    ClassDef(StEffKaon,1)
};
#endif
