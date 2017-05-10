#ifndef StZdcSmdMaker_h
#define StZdcSmdMaker_h

#include "StMaker.h"
#include "TString.h"

class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StRefMultCorr;
class StRunIdEventsDb;
class StCombPID;
class StZdcSmdCut;
class StZdcSmdProManger;
class StZdcSmdCorrection;
class StZdcSmdHistoManger;
class StZdcSmdTree;

class StZdcSmdMaker : public StMaker {
  public:
    StZdcSmdMaker(const char *name, StPicoDstMaker *picoMaker, const Int_t jobCounter, const Int_t Mode, const Int_t Energy, const Int_t Flag_ME);
    virtual ~StZdcSmdMaker();
    
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
    StZdcSmdCut *mZdcSmdCut;
    StZdcSmdProManger *mZdcSmdProManger;
    StZdcSmdCorrection *mZdcSmdCorrection;
    StZdcSmdHistoManger *mZdcSmdHistoManger;
    StZdcSmdTree *mZdcSmdTree;
    
    Int_t mMode;
    Int_t mEnergy;
    Int_t mFlag_ME;

    TString mInPut_Corr_ReCenter;

    TString mOutPut_GainCorrPar;
    TString mOutPut_ReCenterPar;
    TString mOutPut_ShiftPar;
    TString mOutPut_ShiftParFull;
    TString mOutPut_Resolution;
    TString mOutPut_DirectedFlow;
    TString mOutPut_Phi;

    TFile *mFile_GainCorrPar;
    TFile *mFile_ReCenterPar;
    TFile *mFile_ShiftPar;
    TFile *mFile_ShiftParFull;
    TFile *mFile_Resolution;
    TFile *mFile_DirectedFlow;
    TFile *mFile_Phi;

    Int_t mUsedTrackCounter;

    ClassDef(StZdcSmdMaker, 1)
};

#endif
