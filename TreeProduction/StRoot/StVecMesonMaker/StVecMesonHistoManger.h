#ifndef StVecMesonHistoManger_h
#define StVecMesonHistoManger_h

#include "StMessMgr.h"

class TH1F;
class TH2F;

class StVecMesonHistoManger
{
  public:
    StVecMesonHistoManger();
    virtual ~StVecMesonHistoManger();

    void InitQA();
    void FillQA_Detector(Float_t dEdx, Float_t Mass2, Float_t p);
    void FillQA_Event(Float_t vz, Float_t refMult);
    void FillEP_Eta(Float_t Psi2_East, Float_t Psi2_West);
    void FillEP_Full(Float_t Psi2_Full);
    void WriteQA();

    void InitEP();
    void FillEP_Sub(Float_t Psi2East_ReCenter, Float_t Psi2East_Shift, Float_t Psi2West_ReCenter, Float_t Psi2West_Shift);
    void FillEP_Ran(Float_t Psi2RanA_ReCenter, Float_t Psi2RanA_Shift, Float_t Psi2RanB_ReCenter, Float_t Psi2RanB_Shift, Float_t Psi2Full_ReCenter, Float_t Psi2Full_Shift);
    void WriteEP();
    
  private:

    // QA plots
    TH2F *h_mDEdx;
    TH2F *h_mMass2;

    TH1F *h_mEastRaw;
    TH1F *h_mWestRaw;
    TH1F *h_mFullRaw;
    TH1F *h_mVz;
    TH1F *h_mRefMult;

    // event plane distribution
    TH1F *h_mEastReCenter;
    TH1F *h_mWestReCenter;
    TH1F *h_mRanAReCenter;
    TH1F *h_mRanBReCenter;
    TH1F *h_mFullReCenter;

    TH1F *h_mEastShift;
    TH1F *h_mWestShift;
    TH1F *h_mRanAShift;
    TH1F *h_mRanBShift;
    TH1F *h_mFullShift;

  ClassDef(StVecMesonHistoManger,1)
};
#endif
