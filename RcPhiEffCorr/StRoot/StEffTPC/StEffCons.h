#ifndef StEffCons_h
#define StEffCons_h

#include <string>
#include "StEffStruct.h"
#include "TMath.h"

namespace Efficiency
{
  // event cut
  float const mVzMaxMap[7] = {70.0,50.0,70.0,70.0,40.0,40.0,30.0};
  float const mVrMax = 2.0;
  float const mVzVpdDiffMax = 3.0;

  // track cut
  float const mDcaTrMax = 1.0; // for pion, kaon, proton mDcaTrMax = 1.0 for flow
  float const mDcaTrMax_phi = 2.0; // for phi meson mDcaTrMax = 2.0 to fill a tree and apply an additional cut
  int const mHitsDedxMin = 5;
  int const mHitsFitTPCMin = 15;
  int const mHitsMaxTPCMin = 0;
  float const mHitsRatioTPCMin = 0.51;
  float const mEtaMax = 0.5;
  // float const mEtaMax = 0.4;
  float const mPrimPtMin[5] = {0.15,0.2,0.2,0.2,0.2}; // for event plane reconstruction and for pion, kaon, proton: 0.15 for 200 GeV, 0.2 for BES
  float const mGlobPtMin = 0.1; // for phi, Lambda, K0s
  float const mPrimPtMax = 2.0;
  float const mPrimPtWeight = 2.0;
  float const mPrimMomMax = 10.0; // also use for gMom
  float const mMassPion = 0.13957;
  float const mMassKaon = 0.49368;
  float const mMassProton = 0.93827;
  float const mMassPhi = 1.01946;
  float const mSigKaon = 2.5;

  // used constant
  float const mEta_Gap[4] = {0.05,0.10,0.20,0.50};

  int const pt_total = 25; // pt bin
  //                                    0 , 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 ,10 ,11 ,12 ,13 ,14 ,15 ,16 ,17 ,18 ,19 ,10 ,21 ,22, 23, 24
  float const pt_low[pt_total] = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.4,5.8,6.2,6.6,7.2};
  float const pt_high[pt_total]  = {0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.4,5.8,6.2,6.6,7.2,8.0};
  // double const pt_bin[pt_total+1] = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.4,5.8,6.2,6.6,7.2,8.0};

  /*
  int const pt_total = 8; // pt bin
  float const  pt_low[pt_total] = {0.4,0.8,1.2,1.8,2.4,3.0,4.2,5.4};
  float const  pt_up[pt_total]  = {0.8,1.2,1.8,2.4,3.0,4.2,5.4,7.2};
  double const pt_bin[pt_total+1] = {0.4,0.8,1.2,1.8,2.4,3.0,4.2,5.4,7.2};
  */

  // Centrality bin
  int const cent_low[4] = {2,7,4,0}; // 0 = 20-60%, 1 = 0-10%, 2 = 10-40%, 3 = 40-80%
  int const cent_high[4]  = {5,8,6,3}; // 0 = 20-60%, 1 = 0-10%, 2 = 10-40%, 3 = 40-80%

  int const Centrality_total = 4;
  int const Centrality_start = 0;
  int const Centrality_stop  = 1;

  int const EtaGap_total = 4;
  int const EtaGap_start = 0;
  int const EtaGap_stop  = 4;

  int const Charge_total = 2;
  int const Charge_start = 0;
  int const Charge_stop  = 2;

  float const ptMax = 8.0;
  int const BinPt  = 20; // DeltaPt = 0.1 GeV/c
  int const BinEta = 10;
  int const BinPhi = 36;
  int const BinCos = 10;

  double const Psi2_low[3] = {-3.0*TMath::Pi()/2.0,-1.0*TMath::Pi()/2.0,1.0*TMath::Pi()/2.0};
  double const Psi2_up[3]  = {-1.0*TMath::Pi()/2.0, 1.0*TMath::Pi()/2.0,3.0*TMath::Pi()/2.0};

  std::string const mBeamEnergy[7] = {"7GeV","11GeV","19GeV","27GeV","39GeV","62GeV","200GeV"};
  std::string const mParType[2] = {"Phi","KStar"};
}

#endif
