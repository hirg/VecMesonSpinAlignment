#ifndef StStrangenessCons_h
#define StStrangenessCons_h

#include "TString.h"

class Strangeness
{
  public:
    static TString v0_tree[2]; // XuPhiMesonEvent, KstarEvent
    static TString v0_branch[2]; // phi_SpinAlignment_branch, Kstar_SpinAlignment_branch
    static TString Energy[3]; // 0 for 200GeV, 1 for 39GeV, 2 for 27GeV
    static TString Partype[2]; // 0 for Phi, 1 for Kstar
    static TString XFlag[2]; // 0 for SE, 1 for ME 
    static Int_t   mList_Delta; // number of list in one job
    //------------------------------------------------------------------
    // Cuts
    static Float_t mDcaEPMax[5]; // 0: 200 GeV, 1: 39 GeV, 2: 27 GeV, 3: 19.6 GeV, 4: 62.4 GeV
    static Int_t   mTrackMin; // 2
    static Float_t mShiftOrder2[5];
    static Float_t mShiftOrder3[5];
    static Int_t   mHitsFitTPCMin; // 15
    static Int_t   mHitsMaxTPCMin; // 0
    static Float_t mHitsRatioTPCMin; // 0.51
    static Float_t mPrimPtMin[5]; // 0: 0.15(200 GeV), 1-4: 0.2 (BES)
    static Float_t mPrimPtMax; // 2.0
    static Float_t mPrimPtWeight; // 2.0
    static Float_t mPrimMomMax;  // 10.0
    static Float_t mEtaMax; // 1.0
    static Float_t mEta_Gap[4]; // 0 = 0.05, 1 = 0.1, 2 = 0.2, 3 = 0.5
    static Float_t mBeamEnergy[3]; // 200GeV, 39GeV, 27GeV
    static Int_t   mBeamYear[3]; 
    //--------------------------------------------------------------------
    // Histogram
    // pt bin
    static Float_t pt_low[25];
    static Float_t pt_up[25];

    // Centrality bin
    static Int_t   cent_low[4];
    static Int_t   cent_up[4];
    static TString Centrality_01[4]; // {"0080","0010","1040","4080"};
    static TString Centrality_23[4]; // {"0070","0010","1040","4070"};

    // phi-Psi bin
    static Double_t CosThetaStar_low[7];
    static Double_t CosThetaStar_up[7];

    static Int_t pt_total_phi;

    static Int_t Centrality_total;
    static Int_t Centrality_start;
    static Int_t Centrality_stop;

    static Int_t EtaGap_total;
    static Int_t EtaGap_start;
    static Int_t EtaGap_stop;

    static Int_t CosThetaStar_total;

    static Float_t InvMass_low[4];
    static Float_t InvMass_high[4];

  ClassDef(Strangeness,1)
};
#endif
