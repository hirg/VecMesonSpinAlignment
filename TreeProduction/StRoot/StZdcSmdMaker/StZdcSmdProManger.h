#ifndef StZdcSmdProManger_h
#define StZdcSmdProManger_h

#include "TVector2.h"
#include "StMessMgr.h"
#include <string>

class TProfile2D;
class TProfile;

class StZdcSmdProManger
{
  public:
    StZdcSmdProManger();
    virtual ~StZdcSmdProManger();

    // ReCenter Correction
    void InitReCenter();
    void FillReCenterEast(TVector2 qVector, int Cent9, int RunIndex, int vz_sign); // vz_sign = vertex pos/neg
    void FillReCenterWest(TVector2 qVector, int Cent9, int RunIndex, int vz_sign);
    void WriteReCenter();

    // Shift Correction for East/West
    void InitShift();
    void FillShiftEast(TVector2 qVector, int Cent9, int RunIndex, int vz_sign);
    void FillShiftWest(TVector2 qVector, int Cent9, int RunIndex, int vz_sign);
    void WriteShift();

  private:
    // ReCenter Correction | x axis is RunIndex, y axis is Centrality
    TProfile2D *p_mQEastVertical[2]; // 0 = vertex pos/neg
    TProfile2D *p_mQEastHorizontal[2];
    TProfile2D *p_mQWestVertical[2];
    TProfile2D *p_mQWestHorizontal[2];

    // Shift Correction for East/West
    TProfile2D *p_mQEastCos[2][20]; // 0 = vertex pos/neg | 1 = shift correction harmonics
    TProfile2D *p_mQEastSin[2][20];
    TProfile2D *p_mQWestCos[2][20];
    TProfile2D *p_mQWestSin[2][20];


    static string mVStr[2];

    ClassDef(StZdcSmdProManger,1)
};

#endif
