#ifndef StVecMesonCut_h
#define StVecMesonCut_h

#include "TObject.h"
#include "TString.h"
#include "TLorentzVector.h"

class StVecMesonCut : public TObject
{
  public:
    StVecMesonCut(Int_t energy);
    ~StVecMesonCut();

    bool passTrackEP(TLorentzVector, Float_t);
    bool passTrackEtaEast(TLorentzVector, Int_t); // different eta_gap
    bool passTrackEtaWest(TLorentzVector, Int_t);
    bool passPhiEtaEast(TLorentzVector); // eta cut for Phi candidate
    bool passPhiEtaWest(TLorentzVector);

  private:
    Int_t mEnergy;

    ClassDef(StVecMesonCut,1)
};
#endif
