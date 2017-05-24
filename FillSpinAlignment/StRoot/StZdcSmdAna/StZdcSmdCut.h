#ifndef StZdcSmdCut_h
#define StZdcSmdCut_h

#include "TObject.h"
#include "TString.h"
#include "TLorentzVector.h"

class StZdcSmdCut : public TObject
{
  public:
    StZdcSmdCut();
    ~StZdcSmdCut();

    bool passPhiCut(TLorentzVector); // eta cut for phi-meson candidate
    bool passDipAngleCut(TLorentzVector,TLorentzVector); // dip angle cut for phi-meson

  private:

    ClassDef(StZdcSmdCut,1)
};
#endif
