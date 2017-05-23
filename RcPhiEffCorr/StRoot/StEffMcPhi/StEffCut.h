#ifndef StEffCut_h
#define StEffCut_h

#include "StMessMgr.h"
#include "StEffStruct.h"
#include "../../../Utility/StSpinAlignmentCons.h"

class StEffCut
{
  public:
    StEffCut();
    virtual ~StEffCut();

    bool passTrackCutPhi(McVecMeson);
    bool passTrackCutPhi(RcVecMeson);

    bool passTrackCut(McDecayDau);
    bool passTrackCut(RcDecayDau);

    bool passDipAngleCut(McDecayDau,McDecayDau);
    bool passDipAngleCut(RcDecayDau,RcDecayDau);

  private:
    int mEnergy;

    ClassDef(StEffCut,0)
};
#endif
