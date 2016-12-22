#ifndef StVecMesonMaker_h
#define StVecMesonMaker_h

#include "StMaker.h"
#include "TString.h"

class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StRefMultCorr;
class StRunIdEventsDb;
class StCombPID;
class StVecMesonCut;
class StVecMesonProManger;
class StVecMesonCorrection;
class StVecMesonHistoManger;
class StVecMesonTree;

class StVecMesonMaker : public StMaker {
  public:
    StVecMesonMaker(const char *name, StPicoDstMaker *picoMaker, const Int_t jobCounter, const Int_t Mode, const Int_t Energy, const Int_t Flag_ME);
    virtual ~StVecMesonMaker();
    
    virtual Int_t Init();
    virtual Int_t Make();
    virtual void  Clear(Option_t *opt="");
    virtual Int_t Finish();
    
  private:
    StPicoDstMaker *mPicoDstMaker;
    StPicoDst      *mPicoDst;
    StPicoEvent *mPicoEvent;
    static StRefMultCorr *mRefMultCorr;
    StRunIdEventsDb *mRunIdEventsDb;
    StVecMesonCut *mVecMesonCut;
    StVecMesonProManger *mVecMesonProManger;
    StVecMesonCorrection *mVecMesonCorrection;
    StVecMesonHistoManger *mVecMesonHistoManger;
    StVecMesonTree *mVecMesonTree;
    
    Int_t mMode;
    Int_t mEnergy;
    Int_t mFlag_ME;

    TString mInPut_Corr_ReCenter;

    TString mOutPut_ReCenterPar;
    TString mOutPut_Corr_ReCenter;
    TString mOutPut_Corr_Shift;
    TString mOutPut_Phi;

    TFile *mFile_ReCenterPar;
    TFile *mFile_Corr_ReCenter;
    TFile *mFile_Corr_Shift;
    TFile *mFile_Phi;

    Int_t mUsedTrackCounter;

    ClassDef(StVecMesonMaker, 1)
};

#endif
