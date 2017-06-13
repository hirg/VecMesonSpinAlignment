#include "StRoot/StVecMesonMaker/StVecMesonTree.h"
#include "StRoot/StVecMesonMaker/StVecMesonCut.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "StRoot/StPicoDstMaker/StPicoDst.h"
#include "StRoot/StPicoDstMaker/StPicoEvent.h"
#include "StRoot/StPicoDstMaker/StPicoTrack.h"
#include "StRoot/StAlexPhiMesonEvent/StAlexPhiMesonEvent.h"
#include <vector>
#include "TLorentzVector.h"
#include "StThreeVectorF.hh"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "StarClassLibrary/StThreeVector.hh"
#include "StarClassLibrary/SystemOfUnits.h"
#include "TMath.h"
#include "TObject.h"
#include "TVector3.h"

ClassImp(StVecMesonTree)

//------------------------------------------------------------------------------------------------------------------
StVecMesonTree::StVecMesonTree(Int_t energy)
{
  mEnergy = energy;
}

StVecMesonTree::~StVecMesonTree()
{
  /* */
}

//------------------------------------------------------------------------------------------------------------------

void StVecMesonTree::InitPhi()
{
  mVecMesonCut = new StVecMesonCut(mEnergy);
  TString HistName = "Mass2_pt";
  h_Mass2 = new TH2F(HistName.Data(),HistName.Data(),20,0.2,5.0,200,0.98,1.08);

  for(Int_t cent = 0; cent < vmsa::Bin_Centrality; cent++)
  {
    for(Int_t vz = 0; vz < vmsa::Bin_VertexZ; vz++)
    {
      for(Int_t phi_psi = 0; phi_psi < vmsa::Bin_Phi_Psi; phi_psi++)
      {
        mEventCounter2[cent][vz][phi_psi] = 0;
	clear_phi(cent,vz,phi_psi);
      }
    }
  }

  mPhiMesonEvent = new StAlexPhiMesonEvent();
  mTree_Phi = new TTree("PhiMesonEvent","PhiMesonEvent");
  mTree_Phi->Branch("phi_SpinAlignment_branch","StAlexPhiMesonEvent",&mPhiMesonEvent);
  mTree_Phi->SetAutoSave(5000000);
}

//------------------------------------------------------------------------------------------------------------------

void StVecMesonTree::WritePhiMass2()
{
  h_Mass2->Write();
  mTree_Phi->Write("",TObject::kOverwrite);
}

//------------------------------------------------------------------------------------------------------------------

void StVecMesonTree::clear_phi(Int_t cent9, Int_t Bin_vz, Int_t Bin_Psi2)
{
  mPrimaryvertex[cent9][Bin_vz][Bin_Psi2].clear();
  mRefMult[cent9][Bin_vz][Bin_Psi2].clear();
  mCentrality[cent9][Bin_vz][Bin_Psi2].clear();
  mRunId[cent9][Bin_vz][Bin_Psi2].clear();
  mEventId[cent9][Bin_vz][Bin_Psi2].clear();
  mN_prim[cent9][Bin_vz][Bin_Psi2].clear();
  mN_non_prim[cent9][Bin_vz][Bin_Psi2].clear();
  mN_Tof_match[cent9][Bin_vz][Bin_Psi2].clear();
  mZDCx[cent9][Bin_vz][Bin_Psi2].clear();
  mBBCx[cent9][Bin_vz][Bin_Psi2].clear();
  mVzVpd[cent9][Bin_vz][Bin_Psi2].clear();
  mQ2East[cent9][Bin_vz][Bin_Psi2].clear();
  mQ2West[cent9][Bin_vz][Bin_Psi2].clear();
  mQ2Full[cent9][Bin_vz][Bin_Psi2].clear();
  mNumTrackEast[cent9][Bin_vz][Bin_Psi2].clear();
  mNumTrackWest[cent9][Bin_vz][Bin_Psi2].clear();
  mNumTrackFull[cent9][Bin_vz][Bin_Psi2].clear();
  mNumTrackFullEast[cent9][Bin_vz][Bin_Psi2].clear();
  mNumTrackFullWest[cent9][Bin_vz][Bin_Psi2].clear();

  for(Int_t Bin_Event = 0; Bin_Event < vmsa::Buffer_depth; Bin_Event++)
  {
    for(Int_t charge = 0; charge < 2; charge++)
    {
      MEKey key = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,charge);
      mHelix_Kaon[key].clear();
      mMomentum[key].clear();
      mMass2[key].clear();
      mDca[key].clear();
      mNHitsFit[key].clear();
      mNSigmaKaon[key].clear();
    }
  }
  mEventCounter2[cent9][Bin_vz][Bin_Psi2] = 0;
}

void StVecMesonTree::size_phi(Int_t cent9, Int_t Bin_vz, Int_t Bin_Psi2)
{
  LOG_INFO << "Event Buffer: Centrality = " << cent9 << ", VertexZ = " << Bin_vz << ", Psi2 = " << Bin_Psi2 << endm;
  LOG_INFO << "Buffer_depth = " << mEventCounter2[cent9][Bin_vz][Bin_Psi2] << endm;

  LOG_INFO << "Size of primaryVertex = " << mPrimaryvertex[cent9][Bin_vz][Bin_Psi2].size() << endm;;
  LOG_INFO << "Size of refMult       = " << mRefMult[cent9][Bin_vz][Bin_Psi2].size() << endm;;
  LOG_INFO << "---------------------------------------------------------------------------" << endm;

  for(Int_t Bin_Event = 0; Bin_Event < mEventCounter2[cent9][Bin_vz][Bin_Psi2]; Bin_Event++)
  {
    MEKey key = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,0);
    LOG_INFO << "Event Number " << Bin_Event << ":" << endm; 
    LOG_INFO << "Positive Particle:" << endm;
    LOG_INFO << "  Size of Helix_Kplus  = " << mHelix_Kaon[key].size() << endm;;
    LOG_INFO << "  Size of Momentum     = " << mMomentum[key].size() << endm;
    LOG_INFO << "  Size of Mass2        = " << mMass2[key].size() << endm;
    LOG_INFO << "  Size of dca          = " << mDca[key].size() << endm;
    LOG_INFO << "  Size of nHitsFit     = " << mNHitsFit[key].size() << endm;
    LOG_INFO << "  Size of nSigmaKaon   = " << mNSigmaKaon[key].size() << endm;

    key = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,1);
    LOG_INFO << "Negative Particle:" << endm;
    LOG_INFO << "  Size of Helix_Kminus = " << mHelix_Kaon[key].size() << endm;
    LOG_INFO << "  Size of Momentum     = " << mMomentum[key].size() << endm;
    LOG_INFO << "  Size of Mass2        = " << mMass2[key].size() << endm;
    LOG_INFO << "  Size of dca          = " << mDca[key].size() << endm;
    LOG_INFO << "  Size of nHitsFit     = " << mNHitsFit[key].size() << endm;
    LOG_INFO << "  Size of nSigmaKaon   = " << mNSigmaKaon[key].size() << endm;
    LOG_INFO << "---------------------------------------------------------------------------" << endm;
  }
}

//------------------------------------------------------------------------------------------------------------------

void StVecMesonTree::doPhi(Int_t Flag_ME, Int_t cent9, Int_t Bin_vz, Int_t Bin_Psi2) // 0: Same Event, 1: Mix Event
{
  if(Flag_ME == 0) // same event
  {
    for(Int_t Bin_Event = 0; Bin_Event < mEventCounter2[cent9][Bin_vz][Bin_Psi2]; Bin_Event++)
    {
      // event header
      mPhiMesonEvent->clearTrackList();
      mPhiMesonEvent->setPrimaryVertex(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mPhiMesonEvent->setRunId(mRunId[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mPhiMesonEvent->setEventId(mEventId[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mPhiMesonEvent->setRefMult(mRefMult[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mPhiMesonEvent->setCentrality(mCentrality[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mPhiMesonEvent->setN_prim(mN_prim[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mPhiMesonEvent->setN_non_prim(mN_non_prim[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mPhiMesonEvent->setN_Tof_match(mN_Tof_match[cent9][Bin_vz][Bin_Psi2][Bin_Event]);

      // QVector
      mPhiMesonEvent->setQ2East(mQ2East[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mPhiMesonEvent->setQ2West(mQ2West[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mPhiMesonEvent->setQ2Full(mQ2Full[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      // Number of Tracks
      mPhiMesonEvent->setNumTrackEast(mNumTrackEast[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mPhiMesonEvent->setNumTrackWest(mNumTrackWest[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mPhiMesonEvent->setNumTrackFull(mNumTrackFull[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mPhiMesonEvent->setNumTrackFullEast(mNumTrackFullEast[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mPhiMesonEvent->setNumTrackFullWest(mNumTrackFullWest[cent9][Bin_vz][Bin_Psi2][Bin_Event]);

      mPhiMesonEvent->setZDCx(mZDCx[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mPhiMesonEvent->setBBCx(mBBCx[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mPhiMesonEvent->setVzVpd(mVzVpd[cent9][Bin_vz][Bin_Psi2][Bin_Event]);

      // start to select phi candidate in a event
      MEKey key_plus  = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,0);
      MEKey key_minus = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,1);

      TLorentzVector ltrackA, ltrackB;
      for(Int_t n_kplus = 0; n_kplus < mHelix_Kaon[key_plus].size(); n_kplus++) // first track loop over K+ candidates
      {
	StThreeVectorF p_vecA = mHelix_Kaon[key_plus][n_kplus].cat(mHelix_Kaon[key_plus][n_kplus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event]));  // primary momentum
	p_vecA *= mMomentum[key_plus][n_kplus];
	ltrackA.SetXYZM(p_vecA.x(),p_vecA.y(),p_vecA.z(),vmsa::mMassKaon);

	for(Int_t n_kminus = 0; n_kminus < mHelix_Kaon[key_minus].size(); n_kminus++) // second track loop over K- candidates
	{
	  StThreeVectorF p_vecB = mHelix_Kaon[key_minus][n_kminus].cat(mHelix_Kaon[key_minus][n_kminus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event]));  // primary momentum
	  p_vecB *= mMomentum[key_minus][n_kminus];
	  ltrackB.SetXYZM(p_vecB.x(),p_vecB.y(),p_vecB.z(),vmsa::mMassKaon);

	  TLorentzVector trackAB      = ltrackA+ltrackB;
	  Double_t InvMassAB          = trackAB.M();
	  Double_t pt = trackAB.Perp();

	  // fill phi candidate into mTree_Phi
	  if(InvMassAB > vmsa::mMassKaon*2 && InvMassAB < 1.05) 
	  {
	    mPhiMesonTrack = mPhiMesonEvent->createTrack();
	    mPhiMesonTrack->setMass2A(mMass2[key_plus][n_kplus]); // K+
	    mPhiMesonTrack->setMass2B(mMass2[key_minus][n_kminus]); // K-
	    mPhiMesonTrack->setNSigKaonA(mNSigmaKaon[key_plus][n_kplus]); // K+
	    mPhiMesonTrack->setNSigKaonB(mNSigmaKaon[key_minus][n_kminus]); // K-
	    mPhiMesonTrack->setDcaA(mDca[key_plus][n_kplus]); // K+
	    mPhiMesonTrack->setDcaB(mDca[key_minus][n_kminus]); // K-
	    mPhiMesonTrack->setTrackA(ltrackA); // K+
	    mPhiMesonTrack->setTrackB(ltrackB); // K-
	    mPhiMesonTrack->setFlagA(Bin_Event); // K+
	    mPhiMesonTrack->setFlagB(Bin_Event); // K-
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
    for(Int_t Bin_Event_A = 0; Bin_Event_A < mEventCounter2[cent9][Bin_vz][Bin_Psi2]-1; Bin_Event_A++)
    {
      MEKey key_A_plus  = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_A,0);
      MEKey key_A_minus = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_A,1);
      for(Int_t Bin_Event_B = Bin_Event_A+1; Bin_Event_B < mEventCounter2[cent9][Bin_vz][Bin_Psi2]; Bin_Event_B++)
      {
	MEKey key_B_plus  = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_B,0);
	MEKey key_B_minus = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_B,1);

	if(Bin_Event_A == 0 && Bin_Event_B == 1)
	{
	  Int_t Bin_Event = Bin_Event_A;
	  // event header
	  mPhiMesonEvent->clearTrackList();
	  mPhiMesonEvent->setPrimaryVertex(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mPhiMesonEvent->setRunId(mRunId[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mPhiMesonEvent->setEventId(mEventId[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mPhiMesonEvent->setRefMult(mRefMult[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mPhiMesonEvent->setCentrality(mCentrality[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mPhiMesonEvent->setN_prim(mN_prim[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mPhiMesonEvent->setN_non_prim(mN_non_prim[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mPhiMesonEvent->setN_Tof_match(mN_Tof_match[cent9][Bin_vz][Bin_Psi2][Bin_Event]);

	  // QVector
	  mPhiMesonEvent->setQ2East(mQ2East[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mPhiMesonEvent->setQ2West(mQ2West[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mPhiMesonEvent->setQ2Full(mQ2Full[cent9][Bin_vz][Bin_Psi2][Bin_Event]);

	  // Number of Tracks
	  mPhiMesonEvent->setNumTrackEast(mNumTrackEast[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mPhiMesonEvent->setNumTrackWest(mNumTrackWest[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mPhiMesonEvent->setNumTrackFull(mNumTrackFull[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mPhiMesonEvent->setNumTrackFullEast(mNumTrackFullEast[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mPhiMesonEvent->setNumTrackFullWest(mNumTrackFullWest[cent9][Bin_vz][Bin_Psi2][Bin_Event]);

	  mPhiMesonEvent->setZDCx(mZDCx[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mPhiMesonEvent->setBBCx(mBBCx[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mPhiMesonEvent->setVzVpd(mVzVpd[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	}

	TLorentzVector ltrackA, ltrackB;

	// start to mix events
	// mix K+ candidates from A event with K- candidates from B event
	for(Int_t n_kplus = 0; n_kplus < mHelix_Kaon[key_A_plus].size(); n_kplus++) // first track loop over K+ candidates from event A
	{
	  StThreeVectorF p_vecA(mHelix_Kaon[key_A_plus][n_kplus].cat(mHelix_Kaon[key_A_plus][n_kplus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_A])));  // primary momentum
	  p_vecA *= mMomentum[key_A_plus][n_kplus];
	  ltrackA.SetXYZM(p_vecA.x(),p_vecA.y(),p_vecA.z(),vmsa::mMassKaon); // K+

	  for(Int_t n_kminus = 0; n_kminus < mHelix_Kaon[key_B_minus].size(); n_kminus++) // second track loop over K- candidates from event B
	  {
	    StThreeVectorF p_vecB(mHelix_Kaon[key_B_minus][n_kminus].cat(mHelix_Kaon[key_B_minus][n_kminus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_B])));  // primary momentum
	    p_vecB *= mMomentum[key_B_minus][n_kminus];
	    ltrackB.SetXYZM(p_vecB.x(),p_vecB.y(),p_vecB.z(),vmsa::mMassKaon); // K-

	    TLorentzVector trackAB      = ltrackA+ltrackB;
	    Double_t InvMassAB          = trackAB.M();
	    Double_t pt = trackAB.Perp();

	    // fill phi candidate background into mTree_Phi
	    if(InvMassAB > vmsa::mMassKaon*2 && InvMassAB < 1.05) 
	    {
	      mPhiMesonTrack = mPhiMesonEvent->createTrack();
	      mPhiMesonTrack->setMass2A(mMass2[key_A_plus][n_kplus]); // K+
	      mPhiMesonTrack->setMass2B(mMass2[key_B_minus][n_kminus]); // K-
	      mPhiMesonTrack->setNSigKaonA(mNSigmaKaon[key_A_plus][n_kplus]); // K+
	      mPhiMesonTrack->setNSigKaonB(mNSigmaKaon[key_B_minus][n_kminus]); // K-
	      mPhiMesonTrack->setDcaA(mDca[key_A_plus][n_kplus]); // K+
	      mPhiMesonTrack->setDcaB(mDca[key_B_minus][n_kminus]); // K-
	      mPhiMesonTrack->setTrackA(ltrackA); // K+
	      mPhiMesonTrack->setTrackB(ltrackB); // K-
	      mPhiMesonTrack->setFlagA(Bin_Event_A); // K+
	      mPhiMesonTrack->setFlagB(Bin_Event_B); // K-
	    }

	    // Fill histogram with InvMassAB information
	    h_Mass2->Fill(pt,InvMassAB);
	  }
	}

	// mix K- candidates from A event with K+ candidates from B event
	for(Int_t n_kminus = 0; n_kminus < mHelix_Kaon[key_A_minus].size(); n_kminus++) // first track loop over K- candidates from event A
	{
	  StThreeVectorF p_vecA(mHelix_Kaon[key_A_minus][n_kminus].cat(mHelix_Kaon[key_A_minus][n_kminus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_A])));  // primary momentum
	  p_vecA *= mMomentum[key_A_minus][n_kminus];
	  ltrackA.SetXYZM(p_vecA.x(),p_vecA.y(),p_vecA.z(),vmsa::mMassKaon); // K-

	  for(Int_t n_kplus = 0; n_kplus < mHelix_Kaon[key_B_plus].size(); n_kplus++) // second track loop over K+ candidates from event B
	  {
	    StThreeVectorF p_vecB(mHelix_Kaon[key_B_plus][n_kplus].cat(mHelix_Kaon[key_B_plus][n_kplus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_B])));  // primary momentum
	    p_vecB *= mMomentum[key_B_plus][n_kplus];
	    ltrackB.SetXYZM(p_vecB.x(),p_vecB.y(),p_vecB.z(),vmsa::mMassKaon); // K+

	    TLorentzVector trackAB      = ltrackA+ltrackB;
	    Double_t InvMassAB          = trackAB.M();
	    Double_t pt = trackAB.Perp();

	    // fill phi candidate background into mTree_Phi
	    if(InvMassAB > vmsa::mMassKaon*2 && InvMassAB < 1.05) 
	    {
	      mPhiMesonTrack = mPhiMesonEvent->createTrack();
	      mPhiMesonTrack->setMass2A(mMass2[key_B_plus][n_kplus]); // K+
	      mPhiMesonTrack->setMass2B(mMass2[key_A_minus][n_kminus]); // K-
	      mPhiMesonTrack->setNSigKaonA(mNSigmaKaon[key_B_plus][n_kplus]); // K+
	      mPhiMesonTrack->setNSigKaonB(mNSigmaKaon[key_A_minus][n_kminus]); // K-
	      mPhiMesonTrack->setDcaA(mDca[key_B_plus][n_kplus]); // K+
	      mPhiMesonTrack->setDcaB(mDca[key_A_minus][n_kminus]); // K-
	      mPhiMesonTrack->setTrackA(ltrackB); // K+
	      mPhiMesonTrack->setTrackB(ltrackA); // K-
	      mPhiMesonTrack->setFlagA(Bin_Event_B); // K+
	      mPhiMesonTrack->setFlagB(Bin_Event_A); // K-
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

void StVecMesonTree::MixEvent_Phi(Int_t Flag_ME, StPicoDst *pico, Int_t cent9, Float_t vz, Float_t Psi2)
{
  StPicoEvent *event = (StPicoEvent*)pico->event();

  Int_t Bin_vz, Bin_Psi2;

  Float_t vz_start = vmsa::mVzMaxMap[mEnergy];
  Float_t vz_bin = 2*vz_start/vmsa::Bin_VertexZ;

  Float_t psi2_start = TMath::Pi()/2.0;
  Float_t psi2_bin = 2*psi2_start/vmsa::Bin_Phi_Psi;

  for(Int_t i = 0; i < vmsa::Bin_VertexZ; i++)
  {
    if((vz > -1.0*vz_start+i*vz_bin) && (vz <= -1.0*vz_start+(i+1)*vz_bin))
    {
      Bin_vz = i;
    }
  }
  for(Int_t i = 0; i < vmsa::Bin_Phi_Psi; i++)
  {
    if((Psi2 > -1.0*psi2_start+i*psi2_bin) && (Psi2 <= -1.0*psi2_start+(i+1)*psi2_bin))
    {
      Bin_Psi2 = i;
    }
  }

  Int_t Bin_Event = mEventCounter2[cent9][Bin_vz][Bin_Psi2];

  const Double_t MAGFIELDFACTOR = kilogauss;
  const Int_t nTracks = pico->numberOfTracks();

  // store Enent Information
  mPrimaryvertex[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<StThreeVectorF>(event->primaryVertex()));
  mRefMult[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(event->refMult()));
  mCentrality[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(cent9));
  mRunId[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(event->runId()));
  mEventId[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(event->eventId()));
  mN_prim[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mNumber_prim));
  mN_non_prim[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mNumber_non_prim));
  mN_Tof_match[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mNumber_Tof_match));
  mZDCx[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Float_t>(event->ZDCx()));
  mBBCx[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Float_t>(event->BBCx()));
  mVzVpd[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Float_t>(event->vzVpd()));
  mNumTracks[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<UShort_t>(pico->numberOfTracks()));
  mQ2East[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<TVector2>(mQVector2East));
  mQ2West[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<TVector2>(mQVector2West));
  mQ2Full[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<TVector2>(mQVector2Full));
  mNumTrackEast[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mTrackEtaEast));
  mNumTrackWest[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mTrackEtaWest));
  mNumTrackFull[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mTrackFull));
  mNumTrackFullEast[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mTrackFullEast));
  mNumTrackFullWest[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mTrackFullWest));

  // store Track Information
  for(Int_t i = 0; i < nTracks; i ++) // loop over all particles in event
  {
    StPicoTrack *track = pico->track(i);

    if(mVecMesonCut->passTrackPhi(track))
    {
      Float_t Mass2 = mVecMesonCut->getMass2(track);
      Float_t scale_nSigma_factor = vmsa::mSigScaleMap[mEnergy];
      Float_t Polarity = static_cast<Float_t>(track->charge());
      Float_t momentum = track->pMom().mag();
      Float_t Mass2_low;
      Float_t Mass2_up;
      if(momentum < 0.5)
      {
        Mass2_low = 0.4*0.4;
	Mass2_up = 0.6*0.6;
      }
      if(momentum >= 0.5)
      {
	// Mass2_low = 0.277205 - 0.0812931*momentum;
	// Mass2_up = 0.215517 + 0.076801*momentum;
	Mass2_low = 0.19 - 0.06*momentum;
	Mass2_up = 0.34 + 0.04*momentum;
      }

      Int_t charge = 0; // k+
      if(Polarity < 0) charge = 1; // k-


      if(mVecMesonCut->passSigKaonCut(track,scale_nSigma_factor))
      {
	if(
	    (momentum < 0.65 && ((Mass2 > Mass2_low && Mass2 < Mass2_up) || Mass2 < -10.0)) // dE/dx + ToF
	    || (momentum >= 0.65 && (Mass2 > Mass2_low && Mass2 < Mass2_up)) // dE/dx + ToF(always)
	  )
	{
	  MEKey key = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,charge);
	  mMass2[key].push_back(static_cast<Float_t>(mVecMesonCut->getMass2(track))); // mass2
	  mDca[key].push_back(static_cast<Float_t>(track->dca()*track->charge())); // dca*charge 
	  mNHitsFit[key].push_back(static_cast<Float_t>(track->nHitsFit())); // nHitsFit
	  mNSigmaKaon[key].push_back(static_cast<Float_t>((track->nSigmaKaon())*scale_nSigma_factor)); // nSigmaKaon
	  mHelix_Kaon[key].push_back(static_cast<StPhysicalHelixD>(StPhysicalHelixD(track->pMom(),event->primaryVertex(),event->bField()*MAGFIELDFACTOR,track->charge())));// get helix from the pMom 
	  mMomentum[key].push_back(static_cast<Float_t>(track->pMom().mag()));// get helix from the pMom 
	}
      }
    }
  }

  mEventCounter2[cent9][Bin_vz][Bin_Psi2]++;

  if(Flag_ME == 0) // same event
  {
    doPhi(Flag_ME,cent9,Bin_vz,Bin_Psi2);
    clear_phi(cent9,Bin_vz,Bin_Psi2);
  }

  if(Flag_ME == 1) // mix event
  {
    if(mEventCounter2[cent9][Bin_vz][Bin_Psi2] == vmsa::Buffer_depth)
    {
      doPhi(Flag_ME,cent9,Bin_vz,Bin_Psi2);
      clear_phi(cent9,Bin_vz,Bin_Psi2);
    }
  }
}

//------------------------------------------------------------------------------------------------------------------
// pass event information from Maker
void StVecMesonTree::clearEvent()
{
  mNumber_prim = 0;
  mNumber_non_prim = 0;
  mNumber_Tof_match = 0;

  mQVector2East.Set(-999.9,-999.9);
  mQVector2West.Set(-999.9,-999.9);
  mQVector2Full.Set(-999.9,-999.9);

  mTrackEtaEast  = 0;
  mTrackEtaWest  = 0;
  mTrackFull     = 0;
  mTrackFullEast = 0;
  mTrackFullWest = 0;
}

void StVecMesonTree::passEvent(Int_t N_prim, Int_t N_non_prim, Int_t N_Tof_match)
{
  mNumber_prim = N_prim;
  mNumber_non_prim = N_non_prim;
  mNumber_Tof_match = N_Tof_match;
}

void StVecMesonTree::passEventPlane(TVector2 Q2East, TVector2 Q2West, TVector2 Q2Full)
{
  mQVector2East = Q2East;
  mQVector2West = Q2West;
  mQVector2Full = Q2Full;
}

void StVecMesonTree::passNumTrack(Int_t NumTrackEast, Int_t NumTrackWest, Int_t NumTrackFull, Int_t NumTrackFullEast, Int_t NumTrackFullWest)
{
  mTrackEtaEast  = NumTrackEast;
  mTrackEtaWest  = NumTrackWest;
  mTrackFull     = NumTrackFull;
  mTrackFullEast = NumTrackFullEast;
  mTrackFullWest = NumTrackFullWest;
}
//------------------------------------------------------------------------------------------------------------------
