#ifndef StVecMesonTree_h
#define StVecMesonTree_h

#include "StMessMgr.h"
#include "StPhysicalHelixD.hh"
#include "StThreeVectorF.hh"
#include "StVecMesonMEKey.h"
#include <vector>
#include "TVector2.h"

class StPicoDst;
class StAlexPhiMesonEvent;
class StAlexPhiMesonTrack;
class StV0TofCorrection;
class StV0Event;
class StV0Track;
class StVecMesonCut;
class TH1F;
class TH2F;
class TTree;
class TVector2;

class StVecMesonTree
{
  public:
    StVecMesonTree(Int_t energy);
    virtual ~StVecMesonTree();

    void InitPhi();

    void doPhi(Int_t,Int_t,Int_t,Int_t);
    void MixEvent_Phi(Int_t,StPicoDst*,Int_t,Float_t,Float_t);
    void clear_phi(Int_t,Int_t,Int_t);
    void size_phi(Int_t,Int_t,Int_t);

    void WritePhiMass2();

    void clearEvent();
    void passEvent(Int_t,Int_t,Int_t); // N_prim,N_non_prim,N_Tof_match
    void passEventPlane(TVector2,TVector2,TVector2); // qVector ater re-center: east, west, full
    void passNumTrack(Int_t,Int_t,Int_t,Int_t,Int_t); // Number of Tracks: east, west, full, full_east, full_west

  private:
    StVecMesonCut *mVecMesonCut;
    StV0TofCorrection *mTofCorr;
    TH2F *h_Mass2;
    Int_t mEventCounter2[9][10][5]; // 0 = centrality bin, 1 = vertexZ bin, 2 = EP bin

    // 0 = centrality bin, 1 = vertexZ bin, 2 = EP bin, 3 = mixed event bin, 4 = charge bin(0 for pos, 1 for neg) || push_back->track
    vectorHelixMap mHelix_Kaon;
    vectorFloatMap mMomentum;
    vectorFloatMap mMass2;
    vectorFloatMap mDca;
    vectorFloatMap mNHitsFit;
    vectorFloatMap mNSigmaKaon;

    TTree *mTree_Phi;
    StAlexPhiMesonEvent *mPhiMesonEvent;
    StAlexPhiMesonTrack *mPhiMesonTrack;

    // event information | 0 = centrality bin, 1 = vertexZ bin, 2 = EP bin || push_back->event
    std::vector<StThreeVectorF> mPrimaryvertex[9][10][5];
    std::vector<Int_t> mRefMult[9][10][5];
    std::vector<Int_t> mCentrality[9][10][5];
    std::vector<Int_t> mRunId[9][10][5];
    std::vector<Int_t> mEventId[9][10][5];
    std::vector<Int_t> mN_prim[9][10][5];
    std::vector<Int_t> mN_non_prim[9][10][5];
    std::vector<Int_t> mN_Tof_match[9][10][5];
    std::vector<Float_t> mZDCx[9][10][5];
    std::vector<Float_t> mBBCx[9][10][5];
    std::vector<Float_t> mVzVpd[9][10][5];
    std::vector<Float_t> mField[9][10][5];
    std::vector<UShort_t> mNumTracks[9][10][5];
    std::vector<TVector2> mQ2East[9][10][5];
    std::vector<TVector2> mQ2West[9][10][5];
    std::vector<TVector2> mQ2Full[9][10][5];
    std::vector<Int_t> mNumTrackEast[9][10][5];
    std::vector<Int_t> mNumTrackWest[9][10][5];
    std::vector<Int_t> mNumTrackFull[9][10][5];
    std::vector<Int_t> mNumTrackFullEast[9][10][5];
    std::vector<Int_t> mNumTrackFullWest[9][10][5];

    // passing variable
    Int_t mNumber_prim, mNumber_non_prim, mNumber_Tof_match;
    TVector2 mQVector2East, mQVector2West, mQVector2Full;
    Int_t mTrackEtaEast, mTrackEtaWest, mTrackFull, mTrackFullEast, mTrackFullWest;
    Int_t mEnergy;

  ClassDef(StVecMesonTree,1)
};
#endif
