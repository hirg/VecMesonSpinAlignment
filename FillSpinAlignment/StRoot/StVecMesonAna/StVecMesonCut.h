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
    bool passTrackDCA(Float_t, Float_t); // apply dca cuts to Kaon for Systematic errors => default value is 2.0
    bool passTrackEtaEast(TLorentzVector); // different eta_gap
    bool passTrackEtaWest(TLorentzVector);
    bool passPhiEtaEast(TLorentzVector); // eta cut for Phi candidate
    bool passPhiEtaWest(TLorentzVector);

  private:
    Int_t mEnergy;

    ClassDef(StVecMesonCut,1)
};
#endif
