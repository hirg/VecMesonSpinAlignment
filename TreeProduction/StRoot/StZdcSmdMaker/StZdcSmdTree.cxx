#include "StRoot/StZdcSmdMaker/StZdcSmdTree.h"
#include "StRoot/StZdcSmdMaker/StZdcSmdCut.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "StRoot/StPicoDstMaker/StPicoDst.h"
#include "StRoot/StPicoDstMaker/StPicoEvent.h"
#include "StRoot/StPicoDstMaker/StPicoTrack.h"
#include "StRoot/StAlexPhiMesonEvent/StAlexPhiMesonEvent.h"
#include "TLorentzVector.h"
#include "StThreeVectorF.hh"
#include "TH2F.h"
#include "TTree.h"
#include "StarClassLibrary/StThreeVector.hh"
#include "StarClassLibrary/SystemOfUnits.h"
#include "TMath.h"
#include "TObject.h"
#include "TVector3.h"
#include <string>

ClassImp(StZdcSmdTree)

//------------------------------------------------------------------------------------------------------------------
StZdcSmdTree::StZdcSmdTree(int energy)
{
  mEnergy = energy;
}

StZdcSmdTree::~StZdcSmdTree()
{
  /* */
}

//------------------------------------------------------------------------------------------------------------------

void StZdcSmdTree::InitPhi()
{
  mZdcSmdCut = new StZdcSmdCut(mEnergy);
  string HistName = "Mass2_pt";
  h_Mass2 = new TH2F(HistName.c_str(),HistName.c_str(),20,0.2,5.0,200,0.98,1.08);

  for(int i_cent = 0; i_cent < vmsa::Bin_Centrality; ++i_cent)
  {
    for(int i_vz = 0; i_vz < vmsa::Bin_VertexZ; ++i_vz)
    {
      for(int i_Psi = 0; i_Psi < vmsa::Bin_Phi_Psi; ++i_Psi)
      {
        mEventCounter[i_cent][i_vz][i_Psi] = 0;
	clear_phi(i_cent,i_vz,i_Psi);
      }
    }
  }

  mPhiMesonEvent = new StAlexPhiMesonEvent();
  mTree_Phi = new TTree("PhiMesonEvent","PhiMesonEvent");
  mTree_Phi->Branch("phi_SpinAlignment_branch","StAlexPhiMesonEvent",&mPhiMesonEvent);
  mTree_Phi->SetAutoSave(5000000);
}

//------------------------------------------------------------------------------------------------------------------

void StZdcSmdTree::WritePhiMass2()
{
  h_Mass2->Write();
  mTree_Phi->Write("",TObject::kOverwrite);
}

//------------------------------------------------------------------------------------------------------------------

void StZdcSmdTree::clear_phi(int cent9, int Bin_vz, int Bin_Psi)
{
  mPrimaryvertex[cent9][Bin_vz][Bin_Psi].clear();
  mRefMult[cent9][Bin_vz][Bin_Psi].clear();
  mCentrality[cent9][Bin_vz][Bin_Psi].clear();
  mRunId[cent9][Bin_vz][Bin_Psi].clear();
  mEventId[cent9][Bin_vz][Bin_Psi].clear();
  mN_prim[cent9][Bin_vz][Bin_Psi].clear();
  mN_non_prim[cent9][Bin_vz][Bin_Psi].clear();
  mN_Tof_match[cent9][Bin_vz][Bin_Psi].clear();
  mZDCx[cent9][Bin_vz][Bin_Psi].clear();
  mBBCx[cent9][Bin_vz][Bin_Psi].clear();
  mVzVpd[cent9][Bin_vz][Bin_Psi].clear();
  mQEast[cent9][Bin_vz][Bin_Psi].clear();
  mQWest[cent9][Bin_vz][Bin_Psi].clear();
  mQFull[cent9][Bin_vz][Bin_Psi].clear();

  mNumTrackEast[cent9][Bin_vz][Bin_Psi].clear();
  mNumTrackWest[cent9][Bin_vz][Bin_Psi].clear();
  mNumTrackFull[cent9][Bin_vz][Bin_Psi].clear();
  mNumTrackFullEast[cent9][Bin_vz][Bin_Psi].clear();
  mNumTrackFullWest[cent9][Bin_vz][Bin_Psi].clear();

  for(int Bin_Event = 0; Bin_Event < vmsa::Buffer_depth; Bin_Event++)
  {
    for(int charge = 0; charge < 2; charge++)
    {
      MEKey key = MEKey(cent9,Bin_vz,Bin_Psi,Bin_Event,charge);
      mHelix_Kaon[key].clear();
      mMomentum[key].clear();
      mMass2[key].clear();
      mDca[key].clear();
      mNHitsFit[key].clear();
      mNSigmaKaon[key].clear();
    }
  }
  mEventCounter[cent9][Bin_vz][Bin_Psi] = 0;
}

//------------------------------------------------------------------------------------------------------------------

void StZdcSmdTree::doPhi(int Flag_ME, int cent9, int Bin_vz, int Bin_Psi) // 0: Same Event, 1: Mix Event
{
  if(Flag_ME == 0) // same event
  {
    for(int i_event = 0; i_event < mEventCounter[cent9][Bin_vz][Bin_Psi]; i_event++)
    {
      // event header
      mPhiMesonEvent->clearTrackList();
      mPhiMesonEvent->setPrimaryVertex(mPrimaryvertex[cent9][Bin_vz][Bin_Psi][i_event]);
      mPhiMesonEvent->setRunId(mRunId[cent9][Bin_vz][Bin_Psi][i_event]);
      mPhiMesonEvent->setEventId(mEventId[cent9][Bin_vz][Bin_Psi][i_event]);
      mPhiMesonEvent->setRefMult(mRefMult[cent9][Bin_vz][Bin_Psi][i_event]);
      mPhiMesonEvent->setCentrality(mCentrality[cent9][Bin_vz][Bin_Psi][i_event]);
      mPhiMesonEvent->setN_prim(mN_prim[cent9][Bin_vz][Bin_Psi][i_event]);
      mPhiMesonEvent->setN_non_prim(mN_non_prim[cent9][Bin_vz][Bin_Psi][i_event]);
      mPhiMesonEvent->setN_Tof_match(mN_Tof_match[cent9][Bin_vz][Bin_Psi][i_event]);

      // QVector
      mPhiMesonEvent->setQ2East(mQEast[cent9][Bin_vz][Bin_Psi][i_event]);
      mPhiMesonEvent->setQ2West(mQWest[cent9][Bin_vz][Bin_Psi][i_event]);
      mPhiMesonEvent->setQ2Full(mQFull[cent9][Bin_vz][Bin_Psi][i_event]);
      // Number of Tracks
      mPhiMesonEvent->setNumTrackEast(mNumTrackEast[cent9][Bin_vz][Bin_Psi][i_event]); // 0 for 1st event plane
      mPhiMesonEvent->setNumTrackWest(mNumTrackWest[cent9][Bin_vz][Bin_Psi][i_event]);
      mPhiMesonEvent->setNumTrackFull(mNumTrackFull[cent9][Bin_vz][Bin_Psi][i_event]);
      mPhiMesonEvent->setNumTrackFullEast(mNumTrackFullEast[cent9][Bin_vz][Bin_Psi][i_event]);
      mPhiMesonEvent->setNumTrackFullWest(mNumTrackFullWest[cent9][Bin_vz][Bin_Psi][i_event]);

      mPhiMesonEvent->setZDCx(mZDCx[cent9][Bin_vz][Bin_Psi][i_event]);
      mPhiMesonEvent->setBBCx(mBBCx[cent9][Bin_vz][Bin_Psi][i_event]);
      mPhiMesonEvent->setVzVpd(mVzVpd[cent9][Bin_vz][Bin_Psi][i_event]);

      // start to select phi candidate in a event
      MEKey key_plus  = MEKey(cent9,Bin_vz,Bin_Psi,i_event,0);
      MEKey key_minus = MEKey(cent9,Bin_vz,Bin_Psi,i_event,1);

      TLorentzVector ltrackA, ltrackB;
      for(unsigned int i_kplus = 0; i_kplus < mHelix_Kaon[key_plus].size(); i_kplus++) // first track loop over K+ candidates
      {
	StThreeVectorF p_vecA = mHelix_Kaon[key_plus][i_kplus].cat(mHelix_Kaon[key_plus][i_kplus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi][i_event]));  // primary momentum
	p_vecA *= mMomentum[key_plus][i_kplus];
	ltrackA.SetXYZM(p_vecA.x(),p_vecA.y(),p_vecA.z(),vmsa::mMassKaon);

	for(unsigned int i_kminus = 0; i_kminus < mHelix_Kaon[key_minus].size(); i_kminus++) // second track loop over K- candidates
	{
	  StThreeVectorF p_vecB = mHelix_Kaon[key_minus][i_kminus].cat(mHelix_Kaon[key_minus][i_kminus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi][i_event]));  // primary momentum
	  p_vecB *= mMomentum[key_minus][i_kminus];
	  ltrackB.SetXYZM(p_vecB.x(),p_vecB.y(),p_vecB.z(),vmsa::mMassKaon);

	  TLorentzVector trackAB = ltrackA+ltrackB;
	  double InvMassAB = trackAB.M();
	  double pt = trackAB.Perp();

	  // fill phi candidate into mTree_Phi
	  if(InvMassAB > vmsa::mMassKaon*2 && InvMassAB < 1.05) 
	  {
	    mPhiMesonTrack = mPhiMesonEvent->createTrack();
	    mPhiMesonTrack->setMass2A(mMass2[key_plus][i_kplus]); // K+
	    mPhiMesonTrack->setMass2B(mMass2[key_minus][i_kminus]); // K-
	    mPhiMesonTrack->setNSigKaonA(mNSigmaKaon[key_plus][i_kplus]); // K+
	    mPhiMesonTrack->setNSigKaonB(mNSigmaKaon[key_minus][i_kminus]); // K-
	    mPhiMesonTrack->setDcaA(mDca[key_plus][i_kplus]); // K+
	    mPhiMesonTrack->setDcaB(mDca[key_minus][i_kminus]); // K-
	    mPhiMesonTrack->setTrackA(ltrackA); // K+
	    mPhiMesonTrack->setTrackB(ltrackB); // K-
	    mPhiMesonTrack->setFlagA(i_event); // K+
	    mPhiMesonTrack->setFlagB(i_event); // K-
	  }

	  // Fill histogram with InvMassAB information
	  h_Mass2->Fill(pt,InvMassAB);
	}
      }
    }
    mTree_Phi->Fill();
  }

  if(Flag_ME == 1) // mixed event
  {
    for(int i_event_A = 0; i_event_A < mEventCounter[cent9][Bin_vz][Bin_Psi]-1; i_event_A++)
    {
      MEKey key_A_plus  = MEKey(cent9,Bin_vz,Bin_Psi,i_event_A,0);
      MEKey key_A_minus = MEKey(cent9,Bin_vz,Bin_Psi,i_event_A,1);
      for(int i_event_B = i_event_A+1; i_event_B < mEventCounter[cent9][Bin_vz][Bin_Psi]; i_event_B++)
      {
	MEKey key_B_plus  = MEKey(cent9,Bin_vz,Bin_Psi,i_event_B,0);
	MEKey key_B_minus = MEKey(cent9,Bin_vz,Bin_Psi,i_event_B,1);

	if(i_event_A == 0 && i_event_B == 1)
	{
	  int i_event = i_event_A;
	  // event header
	  mPhiMesonEvent->clearTrackList();
	  mPhiMesonEvent->setPrimaryVertex(mPrimaryvertex[cent9][Bin_vz][Bin_Psi][i_event]);
	  mPhiMesonEvent->setRunId(mRunId[cent9][Bin_vz][Bin_Psi][i_event]);
	  mPhiMesonEvent->setEventId(mEventId[cent9][Bin_vz][Bin_Psi][i_event]);
	  mPhiMesonEvent->setRefMult(mRefMult[cent9][Bin_vz][Bin_Psi][i_event]);
	  mPhiMesonEvent->setCentrality(mCentrality[cent9][Bin_vz][Bin_Psi][i_event]);
	  mPhiMesonEvent->setN_prim(mN_prim[cent9][Bin_vz][Bin_Psi][i_event]);
	  mPhiMesonEvent->setN_non_prim(mN_non_prim[cent9][Bin_vz][Bin_Psi][i_event]);
	  mPhiMesonEvent->setN_Tof_match(mN_Tof_match[cent9][Bin_vz][Bin_Psi][i_event]);

	  // QVector
	  mPhiMesonEvent->setQ2East(mQEast[cent9][Bin_vz][Bin_Psi][i_event]);
	  mPhiMesonEvent->setQ2West(mQWest[cent9][Bin_vz][Bin_Psi][i_event]);
	  mPhiMesonEvent->setQ2Full(mQFull[cent9][Bin_vz][Bin_Psi][i_event]);

	  // Number of Tracks
	  mPhiMesonEvent->setNumTrackEast(mNumTrackEast[cent9][Bin_vz][Bin_Psi][i_event]);
	  mPhiMesonEvent->setNumTrackWest(mNumTrackWest[cent9][Bin_vz][Bin_Psi][i_event]);
	  mPhiMesonEvent->setNumTrackFull(mNumTrackFull[cent9][Bin_vz][Bin_Psi][i_event]);
	  mPhiMesonEvent->setNumTrackFullEast(mNumTrackFullEast[cent9][Bin_vz][Bin_Psi][i_event]);
	  mPhiMesonEvent->setNumTrackFullWest(mNumTrackFullWest[cent9][Bin_vz][Bin_Psi][i_event]);

	  mPhiMesonEvent->setZDCx(mZDCx[cent9][Bin_vz][Bin_Psi][i_event]);
	  mPhiMesonEvent->setBBCx(mBBCx[cent9][Bin_vz][Bin_Psi][i_event]);
	  mPhiMesonEvent->setVzVpd(mVzVpd[cent9][Bin_vz][Bin_Psi][i_event]);
	}

	TLorentzVector ltrackA, ltrackB;

	// start to mix events
	// mix K+ candidates from A event with K- candidates from B event
	for(unsigned int i_kplus = 0; i_kplus < mHelix_Kaon[key_A_plus].size(); i_kplus++) // first track loop over K+ candidates from event A
	{
	  StThreeVectorF p_vecA(mHelix_Kaon[key_A_plus][i_kplus].cat(mHelix_Kaon[key_A_plus][i_kplus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi][i_event_A])));  // primary momentum
	  p_vecA *= mMomentum[key_A_plus][i_kplus];
	  ltrackA.SetXYZM(p_vecA.x(),p_vecA.y(),p_vecA.z(),vmsa::mMassKaon); // K+

	  for(unsigned int i_kminus = 0; i_kminus < mHelix_Kaon[key_B_minus].size(); i_kminus++) // second track loop over K- candidates from event B
	  {
	    StThreeVectorF p_vecB(mHelix_Kaon[key_B_minus][i_kminus].cat(mHelix_Kaon[key_B_minus][i_kminus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi][i_event_B])));  // primary momentum
	    p_vecB *= mMomentum[key_B_minus][i_kminus];
	    ltrackB.SetXYZM(p_vecB.x(),p_vecB.y(),p_vecB.z(),vmsa::mMassKaon); // K-

	    TLorentzVector trackAB      = ltrackA+ltrackB;
	    double InvMassAB          = trackAB.M();
	    double pt = trackAB.Perp();

	    // fill phi candidate background into mTree_Phi
	    if(InvMassAB > vmsa::mMassKaon*2 && InvMassAB < 1.05) 
	    {
	      mPhiMesonTrack = mPhiMesonEvent->createTrack();
	      mPhiMesonTrack->setMass2A(mMass2[key_A_plus][i_kplus]); // K+
	      mPhiMesonTrack->setMass2B(mMass2[key_B_minus][i_kminus]); // K-
	      mPhiMesonTrack->setNSigKaonA(mNSigmaKaon[key_A_plus][i_kplus]); // K+
	      mPhiMesonTrack->setNSigKaonB(mNSigmaKaon[key_B_minus][i_kminus]); // K-
	      mPhiMesonTrack->setDcaA(mDca[key_A_plus][i_kplus]); // K+
	      mPhiMesonTrack->setDcaB(mDca[key_B_minus][i_kminus]); // K-
	      mPhiMesonTrack->setTrackA(ltrackA); // K+
	      mPhiMesonTrack->setTrackB(ltrackB); // K-
	      mPhiMesonTrack->setFlagA(i_event_A); // K+
	      mPhiMesonTrack->setFlagB(i_event_B); // K-
	    }

	    // Fill histogram with InvMassAB information
	    h_Mass2->Fill(pt,InvMassAB);
	  }
	}

	// mix K- candidates from A event with K+ candidates from B event
	for(unsigned int i_kminus = 0; i_kminus < mHelix_Kaon[key_A_minus].size(); i_kminus++) // first track loop over K- candidates from event A
	{
	  StThreeVectorF p_vecA(mHelix_Kaon[key_A_minus][i_kminus].cat(mHelix_Kaon[key_A_minus][i_kminus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi][i_event_A])));  // primary momentum
	  p_vecA *= mMomentum[key_A_minus][i_kminus];
	  ltrackA.SetXYZM(p_vecA.x(),p_vecA.y(),p_vecA.z(),vmsa::mMassKaon); // K-

	  for(unsigned int i_kplus = 0; i_kplus < mHelix_Kaon[key_B_plus].size(); i_kplus++) // second track loop over K+ candidates from event B
	  {
	    StThreeVectorF p_vecB(mHelix_Kaon[key_B_plus][i_kplus].cat(mHelix_Kaon[key_B_plus][i_kplus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi][i_event_B])));  // primary momentum
	    p_vecB *= mMomentum[key_B_plus][i_kplus];
	    ltrackB.SetXYZM(p_vecB.x(),p_vecB.y(),p_vecB.z(),vmsa::mMassKaon); // K+

	    TLorentzVector trackAB = ltrackA+ltrackB;
	    double InvMassAB = trackAB.M();
	    double pt = trackAB.Perp();

	    // fill phi candidate background into mTree_Phi
	    if(InvMassAB > vmsa::mMassKaon*2 && InvMassAB < 1.05) 
	    {
	      mPhiMesonTrack = mPhiMesonEvent->createTrack();
	      mPhiMesonTrack->setMass2A(mMass2[key_B_plus][i_kplus]); // K+
	      mPhiMesonTrack->setMass2B(mMass2[key_A_minus][i_kminus]); // K-
	      mPhiMesonTrack->setNSigKaonA(mNSigmaKaon[key_B_plus][i_kplus]); // K+
	      mPhiMesonTrack->setNSigKaonB(mNSigmaKaon[key_A_minus][i_kminus]); // K-
	      mPhiMesonTrack->setDcaA(mDca[key_B_plus][i_kplus]); // K+
	      mPhiMesonTrack->setDcaB(mDca[key_A_minus][i_kminus]); // K-
	      mPhiMesonTrack->setTrackA(ltrackB); // K+
	      mPhiMesonTrack->setTrackB(ltrackA); // K-
	      mPhiMesonTrack->setFlagA(i_event_B); // K+
	      mPhiMesonTrack->setFlagB(i_event_A); // K-
	    }

	    // Fill histogram with InvMassAB information
	    h_Mass2->Fill(pt,InvMassAB);
	  }
	}
      }
    }
    mTree_Phi->Fill();
  }
}

//------------------------------------------------------------------------------------------------------------------

void StZdcSmdTree::MixEvent_Phi(int Flag_ME, StPicoDst *pico, int cent9, float vz, float Psi)
{
  StPicoEvent *event = (StPicoEvent*)pico->event();

  int Bin_vz = -1;
  float vz_start = vmsa::mVzMaxMap[mEnergy];
  float vz_bin = 2*vz_start/vmsa::Bin_VertexZ;
  for(int i_vz = 0; i_vz < vmsa::Bin_VertexZ; ++i_vz)
  {
    if((vz > -1.0*vz_start+i_vz*vz_bin) && (vz <= -1.0*vz_start+(i_vz+1)*vz_bin))
    {
      Bin_vz = i_vz;
    }
  }

  int Bin_Psi = -1;
  float psi_start = TMath::Pi();
  float psi_bin = 2*psi_start/vmsa::Bin_Phi_Psi;
  for(int i_Psi = 0; i_Psi < vmsa::Bin_Phi_Psi; ++i_Psi)
  {
    if((Psi > -1.0*psi_start+i_Psi*psi_bin) && (Psi <= -1.0*psi_start+(i_Psi+1)*psi_bin))
    {
      Bin_Psi = i_Psi;
    }
  }

  int Bin_Event = mEventCounter[cent9][Bin_vz][Bin_Psi];

  const double MAGFIELDFACTOR = kilogauss;
  const unsigned int nTracks = pico->numberOfTracks();

  // store Enent Information
  mPrimaryvertex[cent9][Bin_vz][Bin_Psi].push_back(static_cast<StThreeVectorF>(event->primaryVertex()));
  mRefMult[cent9][Bin_vz][Bin_Psi].push_back(static_cast<int>(event->refMult()));
  mCentrality[cent9][Bin_vz][Bin_Psi].push_back(static_cast<int>(cent9));
  mRunId[cent9][Bin_vz][Bin_Psi].push_back(static_cast<int>(event->runId()));
  mEventId[cent9][Bin_vz][Bin_Psi].push_back(static_cast<int>(event->eventId()));
  mN_prim[cent9][Bin_vz][Bin_Psi].push_back(static_cast<int>(mNumber_prim));
  mN_non_prim[cent9][Bin_vz][Bin_Psi].push_back(static_cast<int>(mNumber_non_prim));
  mN_Tof_match[cent9][Bin_vz][Bin_Psi].push_back(static_cast<int>(mNumber_Tof_match));
  mZDCx[cent9][Bin_vz][Bin_Psi].push_back(static_cast<float>(event->ZDCx()));
  mBBCx[cent9][Bin_vz][Bin_Psi].push_back(static_cast<float>(event->BBCx()));
  mVzVpd[cent9][Bin_vz][Bin_Psi].push_back(static_cast<float>(event->vzVpd()));
  mNumTracks[cent9][Bin_vz][Bin_Psi].push_back(static_cast<unsigned int>(pico->numberOfTracks()));
  mQEast[cent9][Bin_vz][Bin_Psi].push_back(static_cast<TVector2>(mQVectorEast));
  mQWest[cent9][Bin_vz][Bin_Psi].push_back(static_cast<TVector2>(mQVectorWest));
  mQFull[cent9][Bin_vz][Bin_Psi].push_back(static_cast<TVector2>(mQVectorFull));

  mNumTrackEast[cent9][Bin_vz][Bin_Psi].push_back(static_cast<int>(0)); // not used in 1st EP, set to 0
  mNumTrackWest[cent9][Bin_vz][Bin_Psi].push_back(static_cast<int>(0));
  mNumTrackFull[cent9][Bin_vz][Bin_Psi].push_back(static_cast<int>(0));
  mNumTrackFullEast[cent9][Bin_vz][Bin_Psi].push_back(static_cast<int>(0));
  mNumTrackFullWest[cent9][Bin_vz][Bin_Psi].push_back(static_cast<int>(0));

  // store Track Information
  for(unsigned int i_track = 0; i_track < nTracks; ++i_track) // loop over all particles in event
  {
    StPicoTrack *track = pico->track(i_track);

    if(mZdcSmdCut->passTrackPhi(track))
    {
      float Mass2 = mZdcSmdCut->getMass2(track);
      float scale_nSigma_factor = vmsa::mSigScaleMap[mEnergy];
      float Polarity = static_cast<float>(track->charge());
      float momentum = track->pMom().mag();
      float Mass2_low;
      float Mass2_up;
      if(momentum < 0.5)
      {
        Mass2_low = 0.4*0.4;
	Mass2_up = 0.6*0.6;
      }
      if(momentum >= 0.5)
      {
	Mass2_low = 0.277205 - 0.0812931*momentum;
	Mass2_up = 0.215517 + 0.076801*momentum;
      }

      int charge = 0; // k+
      if(Polarity < 0) charge = 1; // k-


      if(mZdcSmdCut->passSigKaonCut(track,scale_nSigma_factor))
      {
	if(
	    (momentum < 0.65 && ((Mass2 > Mass2_low && Mass2 < Mass2_up) || Mass2 < -10.0)) // dE/dx + ToF
	    || (momentum >= 0.65 && (Mass2 > Mass2_low && Mass2 < Mass2_up)) // dE/dx + ToF(always)
	  )
	{
	  MEKey key = MEKey(cent9,Bin_vz,Bin_Psi,Bin_Event,charge);
	  mMass2[key].push_back(static_cast<float>(mZdcSmdCut->getMass2(track))); // mass2
	  mDca[key].push_back(static_cast<float>(track->dca()*track->charge())); // dca*charge 
	  mNHitsFit[key].push_back(static_cast<float>(track->nHitsFit())); // nHitsFit
	  mNSigmaKaon[key].push_back(static_cast<float>((track->nSigmaKaon())*scale_nSigma_factor)); // nSigmaKaon
	  mHelix_Kaon[key].push_back(static_cast<StPhysicalHelixD>(StPhysicalHelixD(track->pMom(),event->primaryVertex(),event->bField()*MAGFIELDFACTOR,track->charge())));// get helix from the pMom 
	  mMomentum[key].push_back(static_cast<float>(track->pMom().mag()));// get helix from the pMom 
	}
      }
    }
  }

  mEventCounter[cent9][Bin_vz][Bin_Psi]++;

  if(Flag_ME == 0) // same event
  {
    doPhi(Flag_ME,cent9,Bin_vz,Bin_Psi);
    clear_phi(cent9,Bin_vz,Bin_Psi);
  }

  if(Flag_ME == 1) // mix event
  {
    if(mEventCounter[cent9][Bin_vz][Bin_Psi] == vmsa::Buffer_depth)
    {
      doPhi(Flag_ME,cent9,Bin_vz,Bin_Psi);
      clear_phi(cent9,Bin_vz,Bin_Psi);
    }
  }
}

//------------------------------------------------------------------------------------------------------------------
// pass event information from Maker
void StZdcSmdTree::clearEvent()
{
  mNumber_prim = 0;
  mNumber_non_prim = 0;
  mNumber_Tof_match = 0;

  mQVectorEast.Set(0.0);
  mQVectorWest.Set(0.0);
  mQVectorFull.Set(0.0);
}

void StZdcSmdTree::passEvent(int N_prim, int N_non_prim, int N_Tof_match)
{
  mNumber_prim = N_prim;
  mNumber_non_prim = N_non_prim;
  mNumber_Tof_match = N_Tof_match;
}

void StZdcSmdTree::passEventPlane(TVector2 QEast, TVector2 QWest, TVector2 QFull)
{
  mQVectorEast = QEast;
  mQVectorWest = QWest;
  mQVectorFull = QFull;
}

//------------------------------------------------------------------------------------------------------------------
