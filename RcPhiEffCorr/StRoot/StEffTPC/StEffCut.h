#ifndef StEffCut_h
#define StEffCut_h

#include "TString.h"
#include "StEffCons.h"

class StEffCut
{
  public:
    StEffCut(Int_t energy);
    virtual ~StEffCut();

    bool passEventCut(McEvent);

    bool passTrackCutPhi(McVecMeson);
    bool passTrackCutPhi(RcVecMeson);

    bool passTrackCut(McDecayDau);
    bool passTrackCut(RcDecayDau);

    bool passPhiEtaEast(McDecayDau); // eta cut
    bool passPhiEtaWest(McDecayDau);

    bool passPhiEtaEast(RcDecayDau); // eta cut
    bool passPhiEtaWest(RcDecayDau);

  private:
    Int_t mEnergy;

    ClassDef(StEffCut,0)
};
#endif
