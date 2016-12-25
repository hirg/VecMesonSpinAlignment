#ifndef StVecMesonProManger_h
#define StVecMesonProManger_h

#include "TVector2.h"
#include "TString.h"
#include "StMessMgr.h"

class TProfile2D;
class TProfile;

class StVecMesonProManger
{
  public:
    StVecMesonProManger();
    virtual ~StVecMesonProManger();

    // ReCenter Correction
    void InitReCenter();
    void FillTrackEast(TVector2 q2Vector, Int_t Cent9, Int_t RunIndex, Int_t i, Float_t pt); // i = vertex pos/neg
    void FillTrackWest(TVector2 q2Vector, Int_t Cent9, Int_t RunIndex, Int_t i, Float_t pt);
    void FillTrackFull(TVector2 q2Vector, Int_t Cent9, Int_t RunIndex, Int_t i, Float_t pt);
    void WriteReCenter();

    // Shift Correction
    void InitShift();
    // Event Plane method
    void FillEventEast_EP(TVector2 Psi2Vector, Int_t Cent9, Int_t RunIndex, Int_t i, Int_t k); // i = vertex pos/neg, k = ShiftOrder
    void FillEventWest_EP(TVector2 Psi2Vector, Int_t Cent9, Int_t RunIndex, Int_t i, Int_t k);
    void FillEventFull_EP(TVector2 Psi2Vector, Int_t Cent9, Int_t RunIndex, Int_t i, Int_t k);
    void WriteShift();

    void InitResolution();
    void FillRes_Sub(Int_t Centrality, Float_t Psi2_East, Float_t Psi2_West);
    void FillRes_Ran(Int_t Centrality, Float_t Psi2_RanA, Float_t Psi2_RanB);
    void WriteResolution();

  private:
    // ReCenter Correction | x axis is RunIndex, y axis is Centrality
    TProfile2D *p_mq2x_East_EP[2]; // 0 = vertex pos/neg
    TProfile2D *p_mq2y_East_EP[2]; // 0 = vertex pos/neg
    TProfile2D *p_mq2x_West_EP[2]; // 0 = vertex pos/neg
    TProfile2D *p_mq2y_West_EP[2]; // 0 = vertex pos/neg
    TProfile2D *p_mq2x_Full_EP[2]; // 0 = vertex pos/neg
    TProfile2D *p_mq2y_Full_EP[2]; // 0 = vertex pos/neg

    // Shift Correction | x axis is RunIndex, y axis is Centrality
    TProfile2D *p_mcos2_East_EP[2][5]; // 0 = vertex pos/neg, 1 = ShiftOrder
    TProfile2D *p_msin2_East_EP[2][5]; // 0 = vertex pos/neg, 1 = ShiftOrder
    TProfile2D *p_mcos2_West_EP[2][5]; // 0 = vertex pos/neg, 1 = ShiftOrder
    TProfile2D *p_msin2_West_EP[2][5]; // 0 = vertex pos/neg, 1 = ShiftOrder
    TProfile2D *p_mcos2_Full_EP[2][5]; // 0 = vertex pos/neg, 1 = ShiftOrder
    TProfile2D *p_msin2_Full_EP[2][5]; // 0 = vertex pos/neg, 1 = ShiftOrder

    TProfile *p_mRes2_Sub; // eta_sub event plane resolution
    TProfile *p_mRes2_Ran; // full event plane resolution <= random sub

    static TString mVStr[2];

    ClassDef(StVecMesonProManger,1)
};

#endif
