#ifndef StZdcSmdCut_h
#define StZdcSmdCut_h

#include "TObject.h"
#include "TString.h"

class StPicoDst;
class StPicoTrack;
class StRefMultCorr;

class StZdcSmdCut : public TObject
{
  public:
    StZdcSmdCut(Int_t energy);
    virtual ~StZdcSmdCut();

    bool passEventCut(StPicoDst*);
    bool passTrackBasic(StPicoTrack*);
    bool passSigPionCut(StPicoTrack*, Float_t);
    bool passSigKaonCut(StPicoTrack*, Float_t);
    bool passSigProntonCut(StPicoTrack*, Float_t);
    bool passTrackPhi(StPicoTrack*);
    bool passTrackV1(StPicoTrack*);
    Int_t getMatchedToF();
    Int_t getNpirm();
    Int_t getNnonprim();
    Float_t getMass2(StPicoTrack*);

  private:
    static StRefMultCorr *mRefMultCorr;
    Int_t mMatchedToF;
    Int_t mN_prim;
    Int_t mN_non_prim;
    Int_t mEnergy;

    ClassDef(StZdcSmdCut,1)
};
#endif
