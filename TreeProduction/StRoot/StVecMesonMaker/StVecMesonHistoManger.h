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
    void FillEP_Eta(Float_t Psi2_East, Float_t Psi2_West);
    void FillEP_Full(Float_t Psi2_Full);
    void WriteQA();
    
  private:

    // QA plots
    TH2F *h_mDEdx;
    TH2F *h_mMass2;

    TH1F *h_East;
    TH1F *h_West;
    TH1F *h_Full;

  ClassDef(StVecMesonHistoManger,1)
};
#endif
