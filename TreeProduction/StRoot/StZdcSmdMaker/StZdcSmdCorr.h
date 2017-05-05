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

    void ReadGainCorr();
    void SetZdcSmdGainCorr(int eastwest,int verthori,int strip,const float zdcsmd);
    float GetZdcSmdGainCorr(int eastwest,int verthori,int strip);
    float GetPosition(int eastwest,int verthori,int strip, int mode);
    TVector2 GetQEast(int mode);
    TVector2 GetQWest(int mode);

    void ReadReCenterCorr();
    void SetZdcSmdCenter(int Cent9, int RunIndex, int vz_sign);

  private:

    int mEnergy;
    float mZdcSmd[2][2][8]; // 0: east/west | 1: vertical(x)/horizontal(y) | 2: 7 slats(x)/8 slats(y);
    float mGainCorrFactor[2][2][8];
    float mCenterEastVertical, mCenterEastHorizontal, mCenterWestVertical, mCenterWestHorizontal;

    TProfile2D *p_mQEastVertical[2]; // vz_sign
    TProfile2D *p_mQEastHorizontal[2];
    TProfile2D *p_mQWestVertical[2];
    TProfile2D *p_mQWestHorizontal[2];

    TFile *mFile_GainCorrPar;
    TFile *mFile_ReCenterPar;

  ClassDef(StZdcSmdCorrection,1)
};

#endif
