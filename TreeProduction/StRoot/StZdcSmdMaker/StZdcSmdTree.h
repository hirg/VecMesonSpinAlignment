#ifndef StZdcSmdTree_h
#define StZdcSmdTree_h

#include "StMessMgr.h"
#include "StPhysicalHelixD.hh"
#include "StThreeVectorF.hh"
#include "StZdcSmdMEKey.h"
#include <vector>
#include "TVector2.h"

class StPicoDst;
class StAlexPhiMesonEvent;
class StAlexPhiMesonTrack;
class StZdcSmdCut;
class TH2F;
class TTree;

class StZdcSmdTree
{
  public:
    StZdcSmdTree(int energy);
    virtual ~StZdcSmdTree();

    void InitPhi();

    void doPhi(int,int,int,int);
    void MixEvent_Phi(int,StPicoDst*,int,float,float);
    void clear_phi(int,int,int);

    void WritePhiMass2();

    void clearEvent();
    void passEvent(int,int,int); // N_prim,N_non_prim,N_Tof_match
    void passEventPlane(TVector2,TVector2,TVector2); // qVector ater shift: east, west, full

  private:
    StZdcSmdCut *mZdcSmdCut;
    TH2F *h_Mass2;
    int mEventCounter[9][10][5]; // 0 = centrality bin, 1 = vertexZ bin, 2 = EP bin

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
    std::vector<int> mRefMult[9][10][5];
    std::vector<int> mCentrality[9][10][5];
    std::vector<int> mRunId[9][10][5];
    std::vector<int> mEventId[9][10][5];
    std::vector<int> mN_prim[9][10][5];
    std::vector<int> mN_non_prim[9][10][5];
    std::vector<int> mN_Tof_match[9][10][5];
    std::vector<float> mZDCx[9][10][5];
    std::vector<float> mBBCx[9][10][5];
    std::vector<float> mVzVpd[9][10][5];
    std::vector<float> mField[9][10][5];
    std::vector<unsigned short> mNumTracks[9][10][5];
    std::vector<TVector2> mQEast[9][10][5];
    std::vector<TVector2> mQWest[9][10][5];
    std::vector<TVector2> mQFull[9][10][5];

    std::vector<int> mNumTrackEast[9][10][5];
    std::vector<int> mNumTrackWest[9][10][5];
    std::vector<int> mNumTrackFull[9][10][5];
    std::vector<int> mNumTrackFullEast[9][10][5];
    std::vector<int> mNumTrackFullWest[9][10][5];

    // passing variable
    int mNumber_prim, mNumber_non_prim, mNumber_Tof_match;
    TVector2 mQVectorEast, mQVectorWest, mQVectorFull;
    int mEnergy;

  ClassDef(StZdcSmdTree,1)
};
#endif
