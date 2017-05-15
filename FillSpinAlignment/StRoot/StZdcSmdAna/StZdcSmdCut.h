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

    bool passPhiCut(TLorentzVector); // eta cut for Phi candidate

  private:

    ClassDef(StZdcSmdCut,1)
};
#endif
