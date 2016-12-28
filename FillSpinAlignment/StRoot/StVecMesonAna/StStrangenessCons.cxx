#include "StStrangenessCons.h"
#include "TMath.h"

ClassImp(Strangeness)

//-----------------------------------------------------------------
TString Strangeness::v0_tree[2]  = {"XuPhiMesonEvent","KstarEvent"};
TString Strangeness::v0_branch[2] = {"phi_flow_branch","Kstar_SpinAlignment_branch"};
TString Strangeness::Energy[3] = {"200GeV","39GeV","27GeV"};
TString Strangeness::Partype[2] = {"Phi","Kstar"};
TString Strangeness::XFlag[2] = {"SE","ME"}; // 0 for SE, 1 for ME 
Int_t   Strangeness::mList_Delta = 50;

//--------------------------------------------------------------------
// Cuts 
Float_t Strangeness::mDcaEPMax[5] = {3.0,1.0,1.0,1.0,1.0}; // for event plane reconstruction: 3.0 for 200GeV, 1.0 for BES
Int_t   Strangeness::mTrackMin = 2;
Float_t Strangeness::mShiftOrder2[5] = {2.0, 4.0, 6.0, 8.0, 10.0};
Float_t Strangeness::mShiftOrder3[5] = {3.0, 6.0, 9.0, 12.0, 15.0};
Int_t   Strangeness::mHitsFitTPCMin = 15;
Int_t   Strangeness::mHitsMaxTPCMin = 0;
Float_t Strangeness::mHitsRatioTPCMin = 0.51;
Float_t Strangeness::mPrimPtMin[5] = {0.15,0.2,0.2,0.2,0.2}; // for event plane reconstruction and for pion, kaon, proton: 0.15 for 200 GeV, 0.2 for BES
Float_t Strangeness::mPrimPtMax = 2.0;
Float_t Strangeness::mPrimPtWeight = 2.0;
Float_t Strangeness::mPrimMomMax = 10.0;
Float_t Strangeness::mEtaMax = 1.0;
Float_t Strangeness::mEta_Gap[4] = {0.05,0.10,0.20,0.50};

Float_t Strangeness::mBeamEnergy[3] = {200.0,39.0,27.0};
Int_t   Strangeness::mBeamYear[3] = {2011,2010,2011};

//--------------------------------------------------------------------
// Histogram
// pt bin
//                                       0 , 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 ,10 ,11 ,12 ,13 ,14 ,15 ,16 ,17 ,18 ,19 ,10 ,21 ,22, 23, 24
Float_t Strangeness::pt_low[25] = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.4,5.8,6.2,6.6,7.2};
Float_t Strangeness::pt_up[25]  = {0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.4,5.8,6.2,6.6,7.2,8.0};

// Centrality bin
Int_t Strangeness::cent_low[4] = {2,7,4,0}; // 0 = 0-80%, 1 = 0-10%, 2 = 10-40%, 3 = 40-80%
Int_t Strangeness::cent_up[4]  = {5,8,6,3}; // 0 = 0-80%, 1 = 0-10%, 2 = 10-40%, 3 = 40-80%
TString Strangeness::Centrality_01[4] = {"0080","0010","1040","4080"};
TString Strangeness::Centrality_23[4] = {"0070","0010","1040","4070"};

// phi-Psi bin
Double_t Strangeness::CosThetaStar_low[7] = {0.0/7.0,1.0/7.0,2.0/7.0,3.0/7.0,4.0/7.0,5.0/7.0,6.0/7.0};
Double_t Strangeness::CosThetaStar_up[7]  = {1.0/7.0,2.0/7.0,3.0/7.0,4.0/7.0,5.0/7.0,6.0/7.0,7.0/7.0};

Int_t Strangeness::pt_total_phi = 25;

Int_t Strangeness::Centrality_total = 4;
Int_t Strangeness::Centrality_start = 0;
Int_t Strangeness::Centrality_stop  = 4;

Int_t Strangeness::EtaGap_total = 4;
Int_t Strangeness::EtaGap_start = 0;
Int_t Strangeness::EtaGap_stop  = 1;

Int_t Strangeness::CosThetaStar_total = 7;

Float_t Strangeness::InvMass_low[4] = {0.98,1.06,1.06,0.4};
Float_t Strangeness::InvMass_high[4] = {1.05,1.20,1.20,0.6};
