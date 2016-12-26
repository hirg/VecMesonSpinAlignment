#ifndef StEffMcPhi_h
#define StEffMcPhi_h
#include "StMessMgr.h"
#include <string>

using namespace std;

class TNtuple;
class TFile;
class StEffHistManger;
class StEffCut;

class StEffMcPhi
{
  public:
    StEffMcPhi(int Energy, long StartEvent, long StopEvent, int PID, int Year, int Cut);
    ~StEffMcPhi();

    void SetInPutFile(const string inputfile);
    void SetOutPutFile(const string outputfile);
    void SetStartEvent(long StartEvent);
    void SetStopEvent(long StopEvent);

    void Init();
    void InitMap();
    void Make();
    void Finish();

  private:
    string mInPutFile;
    string mOutPutFile;
    long mStartEvent;
    long mStopEvent;
    TFile *mFile_InPut;
    TFile *mFile_OutPut;

    int energy;
    int pid;

    static int mInput_flag;

    StEffCut *mEffCut;
    StEffHistManger *mEffHistManger;

    TNtuple *mNtuple; // event header
    float mCentrality;
    float mMcPt;
    float mMcP;
    float mMcEta;
    float mMcY;
    float mMcPhi;
    float mMcInvMass;
    float mMcPid;

    float mKpMcPt;
    float mKpMcEta;
    float mKpMcY;
    float mKpMcPhi;
    float mKpMcM;
    float mKpMcPid;

    float mKpRcPt;
    float mKpRcEta;
    float mKpRcY;
    float mKpRcPhi;
    float mKpRcM;
    float mKpRcTpc;

    float mKmMcPt;
    float mKmMcEta;
    float mKmMcY;
    float mKmMcPhi;
    float mKmMcM;
    float mKmMcPid;

    float mKmRcPt;
    float mKmRcEta;
    float mKmRcY;
    float mKmRcPhi;
    float mKmRcM;
    float mKmRcTpc;

    float mRcPt;
    float mRcP;
    float mRcEta;
    float mRcY;
    float mRcPhi;
    float mRcInvMass;

    ClassDef(StEffMcPhi,1)
};

#endif
