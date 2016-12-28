#ifndef StStrangenessCut_h
#define StStrangenessCut_h

#include "TObject.h"
#include "TString.h"
#include "TLorentzVector.h"

class StStrangenessCut : public TObject
{
  public:
    StStrangenessCut(Int_t energy);
    ~StStrangenessCut();

    bool passTrackEP(TLorentzVector, Float_t);
    bool passTrackEtaEast(TLorentzVector, Int_t); // different eta_gap
    bool passTrackEtaWest(TLorentzVector, Int_t);
    bool passPhiEtaEast(TLorentzVector); // eta cut for Phi candidate
    bool passPhiEtaWest(TLorentzVector);

  private:
    Int_t mEnergy;

    ClassDef(StStrangenessCut,1)
};
#endif
