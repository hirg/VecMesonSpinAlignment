#ifndef StZdcSmdCorrection_h
#define StZdcSmdCorrection_h

#include "TObject.h"
#include "TVector2.h"
#include "TString.h"

class StPicoDst;
class StPicoTrack;
class TProfile2D;
class TProfile;
class TFile;

class StZdcSmdCorrection : public TObject
{
  public:
    StZdcSmdCorrection(int energy);
    virtual ~StZdcSmdCorrection();
    void clear();
    void InitEvent(int Cent9, int RunIndex, int vz_sign);

    void SetZdcSmd(int eastwest,int verthori,int strip,const float zdcsmd);
    float GetZdcSmd(int eastwest,int verthori,int strip);

    void ReadGainCorr();
    void SetZdcSmdGainCorr(int eastwest,int verthori,int strip,const float zdcsmd);
    float GetZdcSmdGainCorr(int eastwest,int verthori,int strip);
    float GetPosition(int eastwest,int verthori,int strip, int mode);
    TVector2 GetQEast(int mode);
    TVector2 GetQWest(int mode);

    void ReadReCenterCorr();
    void SetZdcSmdCenter();

    void ReadShiftCorr();
    TVector2 ApplyZdcSmdShiftCorrEast(TVector2 qVector);
    TVector2 ApplyZdcSmdShiftCorrWest(TVector2 qVector);
    float AngleShift(float Psi_shifted);

    void ReadShiftCorrFull();
    TVector2 ApplyZdcSmdShiftCorrFull(TVector2 qVector);
    TVector2 GetQFull(TVector2 QEast, TVector2 QWest);

    void ReadResolution();
    void CalResolution();
    float GetResolution(int Cent9);

  private:

    int mEnergy;
    int mCent9;
    int mRunIndex;
    int mVz_sign;
    float mZdcSmd[2][2][8]; // 0: east/west | 1: vertical(x)/horizontal(y) | 2: 7 slats(x)/8 slats(y);
    float mGainCorrFactor[2][2][8];

    // ReCenter Correction | x axis is RunIndex, y axis is Centrality
    TProfile2D *p_mQEastVertical[2]; // vz_sign
    TProfile2D *p_mQEastHorizontal[2];
    TProfile2D *p_mQWestVertical[2];
    TProfile2D *p_mQWestHorizontal[2];
    float mCenterEastVertical, mCenterEastHorizontal, mCenterWestVertical, mCenterWestHorizontal;

    // Shift Correction for East/West
    TProfile2D *p_mQEastCos[2][20]; // 0 = vertex pos/neg | 1 = shift correction harmonics
    TProfile2D *p_mQEastSin[2][20];
    TProfile2D *p_mQWestCos[2][20];
    TProfile2D *p_mQWestSin[2][20];

    // Shift Correction for East/West
    TProfile2D *p_mQFullCos[2][20]; // 0 = vertex pos/neg | 1 = shift correction harmonics
    TProfile2D *p_mQFullSin[2][20];

    // charged hadron v1 calculation
    TProfile *p_mResolution;
    float mResolution[9];

    TFile *mFile_GainCorrPar;
    TFile *mFile_ReCenterPar;
    TFile *mFile_ShiftPar;
    TFile *mFile_ShiftParFull;
    TFile *mFile_Resolution;

  ClassDef(StZdcSmdCorrection,1)
};

#endif
