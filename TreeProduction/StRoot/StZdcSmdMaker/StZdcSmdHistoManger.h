#ifndef StZdcSmdHistoManger_h
#define StZdcSmdHistoManger_h

#include "StMessMgr.h"

class TH1F;
class TH2F;

class StZdcSmdHistoManger
{
  public:
    StZdcSmdHistoManger();
    virtual ~StZdcSmdHistoManger();

    void InitQA();
    void FillQA_Event(Float_t vz, Float_t refMult);
    void WriteQA();

    void InitGainCorr();
    void FillGainCorr(int i_eastwest, int i_verthori, int i_slat, int runIndex, float zdcsmd);
    void WriteGainCorr();
    
  private:

    // QA plots
    TH1F *h_mVz;
    TH1F *h_mRefMult;

    TH2F *h_mGainCorr[2][2][8]; // 0: east/west | 1: vertical(x)/horizontal(y) | 2: 7 slats(x)/8 slats(y); | x-axis: runId | y-axis: ADC

  ClassDef(StZdcSmdHistoManger,1)
};
#endif
