#ifndef StZdcSmdHistoManger_h
#define StZdcSmdHistoManger_h

#include "StMessMgr.h"
#include "TVector2.h"

class TH1F;
class TH2F;

class StZdcSmdHistoManger
{
  public:
    StZdcSmdHistoManger();
    virtual ~StZdcSmdHistoManger();

    void InitQA();
    void FillQA_Event(float vz, float refMult);
    void WriteQA();

    // fill adc distribution for gain correction
    void InitGainCorr();
    void FillGainCorr(int i_eastwest, int i_verthori, int i_slat, int runIndex, float zdcsmd);
    void WriteGainCorr();

    // fill QA for raw event plane
    void InitRawEP();
    void FillRawEP(TVector2 QEast, TVector2 QWest, TVector2 QFull, int Cent9, int runIndex);
    void WriteRawEP();

    // fill QA for re-center event plane
    void InitReCenterEP();
    void FillReCenterEP(TVector2 QEast, TVector2 QWest, TVector2 QFull, int Cent9, int runIndex);
    void WriteReCenterEP();

    // fill QA for shift event plane
    void InitShiftEP();
    void FillShiftEP(TVector2 QEast, TVector2 QWest, TVector2 QFull, int Cent9, int runIndex);
    void WriteShiftEP();
    
  private:

    // QA plots
    TH1F *h_mVz;
    TH1F *h_mRefMult;

    TH2F *h_mGainCorr[2][2][8]; // 0: east/west | 1: vertical(x)/horizontal(y) | 2: 7 slats(x)/8 slats(y); | x-axis: runId | y-axis: ADC

    TH2F *h_mRawEast[9];
    TH2F *h_mRawWest[9];
    TH2F *h_mRawFull[9]; // Qwest-QEast

    TH2F *h_mReCenterEast[9];
    TH2F *h_mReCenterWest[9];
    TH2F *h_mReCenterFull[9]; // Qwest-QEast

    TH2F *h_mShiftEast[9];
    TH2F *h_mShiftWest[9];
    TH2F *h_mShiftFull[9]; // Qwest-QEast

  ClassDef(StZdcSmdHistoManger,1)
};
#endif
