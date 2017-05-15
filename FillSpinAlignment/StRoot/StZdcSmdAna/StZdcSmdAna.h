#ifndef StZdcSmdAna_h
#define StZdcSmdAna_h

#include "TObject.h"
#include "TString.h"

class StRefMultCorr;
class TFile;
class TChain;
class StAlexPhiMesonEvent;
class StAlexPhiMesonTrack;
class StZdcSmdCorr;
// class StZdcSmdCut;
class StZdcSmdHistoManger;
class StRunIdEventsDb;

class TStopwatch;

class StZdcSmdAna : public TObject
{
  public:
    StZdcSmdAna(int energy, int X_flag, int List, Long64_t start_event, Long64_t stop_event, int mode); // X_flag: 0 for Same Event, 1 for Mixed Event | List: number of list to use | mode: 0 for phi, 1 for Kstar, 2 for K0S
    ~StZdcSmdAna();

    void setInputDir(const TString inputdir);
    void setOutputfile(const TString outputfile);
    void setInPutList(const TString iInPutList);
    void setStopEvent(const Long64_t StopEvent);
    void setStartEvent(const Long64_t StartEvent);

    void Init();
    void InitSE();
    void InitME();
    void Make();
    void MakePhi();
    void Finish();

  private:
    TString mInputdir;
    TString mOutputfile;
    TString mInPutList;
    Long64_t mStopEvent;
    Long64_t mStartEvent;

    Long64_t mStart_Event;
    Long64_t mStop_Event;

    TFile *mFile_OutPut;
    TChain *mInPut;
    int mEnergy;
    int mX_flag; // 0 for Same Event, 1 for Mixed Event
    int mList;
    int mMode; // 0 for phi, 1 for Kstar
    StAlexPhiMesonEvent *mPhiMeson_event;
    StAlexPhiMesonTrack *mPhiMeson_track;
    StZdcSmdCorr *mZdcSmdCorr;
    // StZdcSmdCut *mZdcSmdCut;
    StZdcSmdHistoManger *mZdcSmdHistoManger;
    StRunIdEventsDb *mRunIdEventsDb;

    static StRefMultCorr *mRefMultCorr;
    static int mInPut_flag;
    static char* VM_EVENT_TREE;
    static char* VM_EVENT_BRANCH;

    TStopwatch *mStopWatch;

  ClassDef(StZdcSmdAna,1)
};
#endif
