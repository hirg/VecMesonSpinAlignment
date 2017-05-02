#ifndef StZdcSmdCorrection_h
#define StZdcSmdCorrection_h

#include "TObject.h"
#include "TVector2.h"
#include "TString.h"

class StPicoDst;
class StPicoTrack;
class TProfile2D;
class TFile;

class StZdcSmdCorrection : public TObject
{
  public:
    StZdcSmdCorrection(Int_t energy);
    virtual ~StZdcSmdCorrection();
    void clear();

    void SetZdcSmd(int eastwest,int verthori,int strip,const float zdcsmd);
    float GetZdcSmd(int eastwest,int verthori,int strip);

  private:

    int mEnergy;

    float mZdcSmd[2][2][8]; // 0: east/west | 1: vertical(x)/horizontal(y) | 2: 7 slats(x)/8 slats(y);

  ClassDef(StZdcSmdCorrection,1)
};

#endif
