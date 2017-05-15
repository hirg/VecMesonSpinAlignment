#ifndef StZdcSmdCorr_h
#define StZdcSmdCorr_h

#include "TObject.h"
#include <string>

class TFile;

class StZdcSmdCorr : public TObject
{
  public:
    StZdcSmdCorr(Int_t energy);
    virtual ~StZdcSmdCorr();

    // Resolution Correction
    void ReadResolution();
    float CalResolution();
    float GetResolution(int Cent9);

  private:
    TFile *mFile_Resolution; // input file for Resolution Correction
    TProfile *p_mResolution;
    float mResolution[9];
    int mEnergy;

  ClassDef(StZdcSmdCorr,1)
};
#endif
