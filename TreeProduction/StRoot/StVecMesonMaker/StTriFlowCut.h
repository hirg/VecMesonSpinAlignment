#ifndef StTriFlowCut_h
#define StTriFlowCut_h

#include "TObject.h"
#include "TString.h"

class StPicoDst;
class StPicoTrack;
class StRefMultCorr;

class StTriFlowCut : public TObject
{
  public:
    StTriFlowCut(Int_t energy);
    virtual ~StTriFlowCut();

    bool passEventCut(StPicoDst*);
    bool passTrackBasic(StPicoTrack*);
    bool passTrackEP(StPicoTrack*);
    bool passTrackCut(StPicoTrack*);
    bool passTrackCutSys(StPicoTrack*,Int_t,Int_t);
    bool passPIDCut(StPicoTrack*);
    bool passSigElectronCut(StPicoTrack*, Float_t);
    bool passSigPionCut(StPicoTrack*, Float_t);
    bool passSigKaonCut(StPicoTrack*, Float_t);
    bool passSigProntonCut(StPicoTrack*, Float_t);
    bool passSigProntonCutSys(StPicoTrack*, Float_t, Int_t);
    bool passTrackPhi(StPicoTrack*);
    bool passTrackV0(StPicoTrack*);
    Int_t getMatchedToF();
    Int_t getNpirm();
    Int_t getNnonprim();
    Float_t getMass2(StPicoTrack*);
    Float_t getV0Mass2(StPicoTrack*);

  private:
    static StRefMultCorr *mRefMultCorr;
    Int_t mMatchedToF;
    Int_t mN_prim;
    Int_t mN_non_prim;
    Int_t mEnergy;

    ClassDef(StTriFlowCut,1)
};
#endif
