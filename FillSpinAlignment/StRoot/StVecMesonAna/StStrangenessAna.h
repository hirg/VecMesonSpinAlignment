#ifndef StStrangenessAna_h
#define StStrangenessAna_h

#include "TObject.h"
#include "TString.h"

class StRefMultCorr;
class TFile;
class TChain;
class StAlexPhiMesonEvent;
class StAlexPhiMesonTrack;
class StStrangenessCorr;
class StStrangenessCut;
class StStrangenessHistoManger;
class StRunIdEventsDb;
class StV0Event;
class StV0Track;

class StStrangenessAna : public TObject
{
  public:
    StStrangenessAna(Int_t energy, Int_t X_flag, Int_t List, Long64_t start_event, Long64_t stop_event, Int_t mode, Int_t flag_Embedding); // energy: 0 for 200GeV, 1 for 39GeV | X_flag: 0 for Same Event, 1 for Mixed Event | List: number of list to use | mode: 0 for phi, 1 for Kstar | flag_Embedding: 0 for production, 1 for embedding
    ~StStrangenessAna();

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
    Int_t mEnergy;
    Int_t mX_flag; // 0 for Same Event, 1 for Mixed Event
    Int_t mList;
    Int_t mMode; // 0 for phi, 1 for Kstar
    Int_t mFlag_Embedding; // 0 for production, 1 for embedding
    StAlexPhiMesonEvent *mXuPhiMeson_event;
    StAlexPhiMesonTrack *mXuPhiMeson_track;
    StStrangenessCorr *mStrangenessCorr;
    StStrangenessCut *mStrangenessCut;
    StStrangenessHistoManger *mStrangenessHistoManger;
    StRunIdEventsDb *mRunIdEventsDb;

    static StRefMultCorr *mRefMultCorr;
    static Int_t mInPut_flag;
    static char* XUV0_EVENT_TREE;
    static char* XUV0_EVENT_BRANCH;

  ClassDef(StStrangenessAna,1)
};
#endif
