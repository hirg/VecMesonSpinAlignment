#ifndef StTriFlowHistoManger_h
#define StTriFlowHistoManger_h

#include "StMessMgr.h"

class TH2F;

class StTriFlowHistoManger
{
  public:
    StTriFlowHistoManger(Int_t energy);
    virtual ~StTriFlowHistoManger();

    void InitQA_Detector();
    void FillQA_Detector(Float_t dEdx, Float_t Mass2, Float_t p);
    void WriteQA_Detector();
    
  private:

    // QA plots
    TH2F *h_mDEdx;
    TH2F *h_mMass2;

    Int_t mEnergy;

  ClassDef(StTriFlowHistoManger,1)
};
#endif
