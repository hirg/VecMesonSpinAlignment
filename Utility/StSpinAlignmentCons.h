#ifndef StSpinAlignmentCons_h
#define StSpinAlignmentCons_h

#include <string>

namespace vmsa
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
  float const mEtaMax = 1.0;
  // float const mEtaMax = 0.5;
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
  int const pt_start = 0;
  int const pt_stop  = 25;
  float const ptRawStart[pt_total] = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.4,5.8,6.2,6.6,7.2};
  float const ptRawStop[pt_total]  = {0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.4,5.8,6.2,6.6,7.2,8.0};
  float const pt_bin[pt_total+1] = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.4,5.8,6.2,6.6,7.2,8.0};

  string const Centrality[9] = {"70%-80%","60%-70%","50%-60%","40%-50%","30%-40%","20%-30%","10%-20%","5%-10%","0%-5%"}; // Centrality bin
  int const Cent_start = 0;
  int const Cent_stop  = 1;
  int const cent_low[5] = {2,0,7,4,0}; // 0 = 20%-60%, 1 = 0-80%, 2 = 0-10%, 3 = 10-40%, 4 = 40-80%
  int const cent_up[5]  = {5,8,8,6,3}; // 0 = 20%-60%, 1 = 0-80%, 2 = 0-10%, 3 = 10-40%, 4 = 40-80%

  int const CTS_total = 7; // cos(theta*) bin
  int const CTS_start = 0;
  int const CTS_stop  = 7;

  int const Eta_start = 0; // EtaGap bin
  int const Eta_stop  = 1;
  int const Eta_QA    = 0;

  int const Sys_start = 0; // SysError bin
  int const Sys_stop  = 1;
  int const Sys_QA    = 0;

  // shared constant
  std::string const mBeamEnergy[7] = {"7GeV","11GeV","19GeV","27GeV","39GeV","62GeV","200GeV"};
  std::string const mPID[2] = {"Phi","KStar"};
  float const Norm_Start[2] = {1.04,0.41};
  float const Norm_Stop[2]  = {1.05,0.46};
  float const BW_Start[2] = {0.994,1.0};
  float const BW_Stop[2]  = {1.050,1.0};
  float const InvMass[2] = {1.019,0.892};
  float const Width[2]   = {0.00426,0.0487};
  float const nSigVec = 2.0;

  // float const ptMax = 12.0;
  float const ptMin = 0.2;
  float const ptMax = 5.0;
  int const BinPt  = 20;
  int const BinEta = 10;
  int const BinY = 20;
  int const BinPhi = 36;

  // used in McPhiResCorr
  double const acceptanceRapidity = 1.0;
  float const rhoDelta = 0.01;
  int const decayChannelsFirst[2] = {656,123}; // 0 phi, 1 K*
  int const decayChannelsSecond[2] = {666,133};
  int const decayChannels[2] = {656,123}; // 0: phi->K+K-

  std::string const mParType[2] = {"Kplus","Kminus"};
}

#endif
