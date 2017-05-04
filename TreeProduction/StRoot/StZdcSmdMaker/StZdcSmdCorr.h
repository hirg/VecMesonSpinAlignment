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
    StZdcSmdCorrection(int energy);
    virtual ~StZdcSmdCorrection();
    void clear();

    void SetZdcSmd(int eastwest,int verthori,int strip,const float zdcsmd);
    float GetZdcSmd(int eastwest,int verthori,int strip);

    void InitGainCorr();
    void SetZdcSmdGainCorr(int eastwest,int verthori,int strip,const float zdcsmd);
    void GetZdcSmdGainCorr(int eastwest,int verthori,int strip);
    float GetPosition(int eastwest,int verthori,int strip, int mode);
    TVector2 GetQEast(int mode);
    TVector2 GetQWest(int mode);

  private:

    int mEnergy;
    float mZdcSmd[2][2][8]; // 0: east/west | 1: vertical(x)/horizontal(y) | 2: 7 slats(x)/8 slats(y);
    float mGainCorrFactor[2][2][8];

    TFile *mFile_GainCorr;

  ClassDef(StZdcSmdCorrection,1)
};

#endif
