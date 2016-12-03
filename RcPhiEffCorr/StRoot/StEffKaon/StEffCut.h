#ifndef StEffCut_h
#define StEffCut_h

#include "StMessMgr.h"
#include "StEffStruct.h"
#include "../../../Utility/StSpinAlignmentCons.h"

class StEffCut
{
  public:
    StEffCut(int energy);
    virtual ~StEffCut();

    bool passEventCut(McEvent);

    bool passTrackCut(McDecayDau);
    bool passTrackCut(RcDecayDau);

  private:
    int mEnergy;

    ClassDef(StEffCut,0)
};
#endif
