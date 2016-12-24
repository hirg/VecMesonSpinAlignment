#include "StVecMesonMaker.h"
#include "../../../Utility/StSpinAlignmentCons.h"
#include "StVecMesonCut.h"
#include "StVecMesonProManger.h"
#include "StVecMesonCorr.h"
#include "StVecMesonHistoManger.h"
#include "StVecMesonTree.h"
#include "StRoot/StPicoDstMaker/StPicoDst.h"
#include "StRoot/StPicoDstMaker/StPicoEvent.h"
#include "StRoot/StPicoDstMaker/StPicoTrack.h"
#include "StRoot/StPicoDstMaker/StPicoV0.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"
#include "StRoot/StRunIdEventsDb/StRunIdEventsDb.h"
#include "StThreeVectorF.hh"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TVector3.h"
#include "TMath.h"
#include "StMessMgr.h"
#include <algorithm>

ClassImp(StVecMesonMaker)

StRefMultCorr* StVecMesonMaker::mRefMultCorr = NULL;
//-----------------------------------------------------------------------------
StVecMesonMaker::StVecMesonMaker(const char* name, StPicoDstMaker *picoMaker, const Int_t jobCounter, const Int_t Mode, const Int_t energy, const Int_t flag_ME)
  : StMaker(name)
{
  mPicoDstMaker = picoMaker;
  mPicoDst = 0;
  mMode = Mode;
  mEnergy = energy;
  mFlag_ME = flag_ME;

  if(mMode == 0)
  {
    mOutPut_ReCenterPar = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/RecenterParameter/file_%s_ReCenterPar_%d.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str(),jobCounter);
  }
  if(mMode == 1)
  {
    mOutPut_Corr_Shift = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/Correction/Shift/file_%s_Corr_Shift_%d.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str(),jobCounter); 

    mOutPut_Corr_ReCenter = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/Correction/ReCenter/file_%s_Corr_ReCenter_%d.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str(),jobCounter); 
  }
  if(mMode == 3)
  {
    mOutPut_Phi = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/Phi/file_%s_Phi_%s_%d.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::MixEvent[mFlag_ME].Data(),jobCounter); 
  }
}

//----------------------------------------------------------------------------- 
StVecMesonMaker::~StVecMesonMaker()
{ /*  */ }

//----------------------------------------------------------------------------- 
Int_t StVecMesonMaker::Init() 
{
  if(!mRefMultCorr)
  {
    mRefMultCorr = CentralityMaker::instance()->getRefMultCorr();
  }
  mVecMesonCut = new StVecMesonCut(mEnergy);
  mVecMesonProManger = new StVecMesonProManger();
  mVecMesonCorrection = new StVecMesonCorrection(mEnergy);
  mVecMesonHistoManger = new StVecMesonHistoManger();

  if(mMode == 0)
  {
    mFile_ReCenterPar = new TFile(mOutPut_ReCenterPar.Data(),"RECREATE");
    mFile_ReCenterPar->cd();
    mVecMesonProManger->InitReCenter();
    mVecMesonHistoManger->InitQA();
  }

  if(mMode == 1)
  {
    mUsedTrackCounter = 0;
    mVecMesonCorrection->InitReCenterCorrection(mEnergy);
    mFile_Corr_Shift = new TFile(mOutPut_Corr_Shift.Data(),"RECREATE");
    mVecMesonProManger->InitShift();
  }
  if(mMode == 3)
  {
    mVecMesonTree = new StVecMesonTree(mEnergy);
    mFile_Phi = new TFile(mOutPut_Phi.Data(),"RECREATE");
    mFile_Phi->cd();
    mVecMesonTree->InitPhi();
    mVecMesonCorrection->InitReCenterCorrection(mEnergy);
    mVecMesonCorrection->InitShiftCorrection(mEnergy);
  }

  return kStOK;
}

//----------------------------------------------------------------------------- 
Int_t StVecMesonMaker::Finish() 
{
  if(mMode == 0)
  {
    if(mOutPut_ReCenterPar != "")
    {
      mFile_ReCenterPar->cd();
      mVecMesonProManger->WriteReCenter();
      mVecMesonHistoManger->WriteQA();
      mFile_ReCenterPar->Close();
    }
  }
  if(mMode == 1)
  {
    if(mOutPut_Corr_Shift != "")
    {
      mFile_Corr_Shift->cd();
      mVecMesonProManger->WriteShift();
      mFile_Corr_Shift->Close();
    }
  }
  if(mMode == 3)
  {
    if(mOutPut_Phi != "")
    {
      mFile_Phi->cd();
      mVecMesonTree->WritePhiMass2();
      mFile_Phi->Close();
    }
  }

  return kStOK;
}

//----------------------------------------------------------------------------- 
void StVecMesonMaker::Clear(Option_t *opt) {
}

//----------------------------------------------------------------------------- 
Int_t StVecMesonMaker::Make() 
{
  if(!mPicoDstMaker) 
  {
    LOG_WARN << " No PicoDstMaker! Skip! " << endm;
    return kStWarn;
  }

  mPicoDst = mPicoDstMaker->picoDst();
  if(!mPicoDst) 
  {
    LOG_WARN << " No PicoDst! Skip! " << endm;
    return kStWarn;
  }

  mPicoEvent = (StPicoEvent*)mPicoDst->event();
  if(!mPicoEvent)
  {
    LOG_WARN << " No PicoEvent! Skip! " << endm;
    return kStWarn;
  }

  // RefMult
  Int_t runId = mPicoEvent->runId();
  Int_t refMult = mPicoEvent->refMult();
  Float_t vz = mPicoEvent->primaryVertex().z();
  Float_t zdcX = mPicoEvent->ZDCx();
  mRefMultCorr->init(runId);
  if(mEnergy == 0) mRefMultCorr->initEvent(refMult,vz,zdcX); // for 200 GeV
  if(mEnergy != 0) mRefMultCorr->initEvent(refMult,vz); // for BES Energy

  // vz sign
  Int_t vz_sign;
  if(vz > 0.0)
  {
    vz_sign = 0;
  }
  else
  {
    vz_sign = 1;
  }

  // runIndex
  mRunIdEventsDb = StRunIdEventsDb::Instance((Float_t)mPicoEvent->energy(),(Float_t)mPicoEvent->year());
  const Int_t runIndex = mRunIdEventsDb->getRunIdIndex(runId); // expensive
//  cout << runIndex << endl;
//  cout << mRunIdEventsDb->getTotalNrRunIds() << endl;

  // Event Cut
  if(mVecMesonCut->passEventCut(mPicoDst)) // event cut
  {
    const Int_t nTracks = mPicoDst->numberOfTracks();
    const Int_t cent9 = mRefMultCorr->getCentralityBin9();
//    if(cent9 < 0) cout << cent9 << endl;
    const Double_t reweight = mRefMultCorr->getWeight();
    const Int_t nToFMatched = mVecMesonCut->getMatchedToF();

    for(Int_t i = 0; i < nTracks; i++) // track loop
    {
      StPicoTrack *track = (StPicoTrack*)mPicoDst->track(i);

      if(mMode == 0)
      {
	Float_t eta = track->pMom().pseudoRapidity();
	if(fabs(eta) < vmsa::mEtaMax && track->dca() < 3.0)
	{
	  Float_t Mass2 = mVecMesonCut->getMass2(track);
	  Float_t dEdx = track->dEdx();
	  Float_t p = track->pMom().mag();
	  mVecMesonHistoManger->FillQA_Detector(dEdx,Mass2,p);
	}
      }
      if(mVecMesonCut->passTrackEP(track)) // track cut
      {
	if(mMode == 0) // fill re-center parameter
	{
	  Float_t pt = track->pMom().perp();

	  if(mVecMesonCorrection->passTrackFull(track))
	  {
	    TVector2 q2Vector_Full = mVecMesonCorrection->calq2Vector(track);
	    mVecMesonProManger->FillTrackFull(q2Vector_Full,cent9,runIndex,vz_sign,pt);
	    mVecMesonCorrection->addTrack_FullRaw(track,cent9,runIndex);
	  }

	  if(mVecMesonCorrection->passTrackEtaEast(track)) // neg eta sub
	  {
	    TVector2 q2Vector_East = mVecMesonCorrection->calq2Vector(track);
	    mVecMesonProManger->FillTrackEast(q2Vector_East,cent9,runIndex,vz_sign,pt);
	    mVecMesonCorrection->addTrack_EastRaw(track,cent9,runIndex);
	  }
	  if(mVecMesonCorrection->passTrackEtaWest(track)) // pos eta sub
	  {
	    TVector2 q2Vector_West = mVecMesonCorrection->calq2Vector(track);
	    mVecMesonProManger->FillTrackWest(q2Vector_West,cent9,runIndex,vz_sign,pt);
	    mVecMesonCorrection->addTrack_WestRaw(track,cent9,runIndex);
	  }
	}
	else // calculate Q Vector after recentering for full event and eta sub event
	{
	  if(mVecMesonCorrection->passTrackFull(track))
	  {
	    mVecMesonCorrection->addTrack_Full(track,cent9,runIndex,vz_sign);
	    mUsedTrackCounter++;
	  }
	  if(mVecMesonCorrection->passTrackEtaEast(track)) // neg eta sub
	  {
	    mVecMesonCorrection->addTrack_East(track,cent9,runIndex,vz_sign);
	  }
	  if(mVecMesonCorrection->passTrackEtaWest(track)) // pos eta sub
	  {
	    mVecMesonCorrection->addTrack_West(track,cent9,runIndex,vz_sign);
	  }
	}
      }
    }

    if(mMode == 0) // fill raw EP
    {
      if(mVecMesonCorrection->passTrackEtaNumCut())
      {
	TVector2 Q2East = mVecMesonCorrection->getQVectorRaw(j,0); // 0 = eta_gap, 1 = east/west
	Float_t Psi2_East = TMath::ATan2(Q2East.Y(),Q2East.X())/2.0;
	TVector2 Q2West = mVecMesonCorrection->getQVectorRaw(j,1); // 0 = eta_gap, 1 = east/west
	Float_t Psi2_West = TMath::ATan2(Q2West.Y(),Q2West.X())/2.0;
	mVecMesonHistoManger->FillEP_Eta(Psi2_East,Psi2_West,j);
      }
      if(mVecMesonCorrection->passTrackFullNumCut())
      {
	TVector2 Q2Full = mVecMesonCorrection->getQVectorRaw(-1,2);
	Float_t Psi2_Full = TMath::ATan2(Q2Full.Y(),Q2Full.X())/2.0;
	mVecMesonHistoManger->FillEP_Full(Psi2_Full);
	mVecMesonCorrection->clear();
      }
    }
    if(mMode == 1) // calculate Q vector after recentering for Random Sub Event
    {
      Int_t iTrack[mUsedTrackCounter];
      Float_t ranCounter = (Float_t)mUsedTrackCounter/2.0 - 1;
      for(Int_t i = 0; i < mUsedTrackCounter; i++)
      {
        iTrack[i] = i;
      }
      std::srand(time(0));
      std::random_shuffle(iTrack,iTrack+mUsedTrackCounter);
      mUsedTrackCounter = 0;
      for(Int_t i = 0; i < nTracks; i++) // track loop
      {
	StPicoTrack *track = (StPicoTrack*)mPicoDst->track(i);
	if(mVecMesonCut->passTrackEP(track)) // track cut
	{
	  if(mVecMesonCorrection->passTrackFull(track))
	  {
	    if((Float_t)iTrack[mUsedTrackCounter] > ranCounter) // Sub Event A
	    {
	      mVecMesonCorrection->addTrack_A(track,cent9,runIndex,vz_sign);
	    }
	    else // Sub Event B
	    {
	      mVecMesonCorrection->addTrack_B(track,cent9,runIndex,vz_sign);
	    }
	    mUsedTrackCounter++;
	  }
	}
      }

      // full event shift parameter
      if(mVecMesonCorrection->passTrackFullNumCut())
      {
	for(Int_t k = 0; k < 5; k++) // ShiftOrder loop
	{
	  TVector2 Psi2Vector_Full_EP = mVecMesonCorrection->calPsi2_Full_EP(k);
	  mVecMesonProManger->FillEventFull_EP(Psi2Vector_Full_EP,cent9,runIndex,vz_sign,k);
	}
      }

      // eta sub event shift parameter
      if(mVecMesonCorrection->passTrackEtaNumCut(j))
      {
	for(Int_t k = 0; k < 5; k++)
	{
	  TVector2 Psi2Vector_East_EP = mVecMesonCorrection->calPsi2_East_EP(k);
	  mVecMesonProManger->FillEventEast_EP(Psi2Vector_East_EP,cent9,runIndex,vz_sign,k);

	  TVector2 Psi2Vector_West_EP = mVecMesonCorrection->calPsi2_West_EP(k);
	  mVecMesonProManger->FillEventWest_EP(Psi2Vector_West_EP,cent9,runIndex,vz_sign,k);
	}
      }
      mVecMesonCorrection->clear();
      mUsedTrackCounter = 0;
    }

    if(mMode == 3)
    { // phi meson
      Float_t Psi2_East;
      TVector2 Q2East, Q2West;
      Int_t NumTrackEast, NumTrackWest; 
      Q2East.Set(-999.9,-999.9); // initialize Q Vector to unreasonable value
      Q2West.Set(-999.9,-999.9);
      NumTrackEast = 0;
      NumTrackWest = 0;
      if(mVecMesonCorrection->passTrackEtaNumCut())
      {
	// get QVector of sub event
	Q2East = mVecMesonCorrection->getQVector(0); // east
	Q2West = mVecMesonCorrection->getQVector(1); // west
	NumTrackEast = mVecMesonCorrection->getNumTrack(0);
	NumTrackWest = mVecMesonCorrection->getNumTrack(1);
      }

      if(mVecMesonCorrection->passTrackEtaNumCut())
      {
	Psi2_East = mVecMesonCorrection->calShiftAngle2East_EP(runIndex,cent9,vz_sign);

	// get N_prim, N_non_prim, N_Tof_match
	Int_t N_prim = mVecMesonCut->getNpirm();
	Int_t N_non_prim = mVecMesonCut->getNnonprim();
	Int_t N_Tof_match = mVecMesonCut->getMatchedToF();

	// pass the event information to StVecMesonTree
	mVecMesonTree->clearEvent();
	mVecMesonTree->passEvent(N_prim, N_non_prim, N_Tof_match);

	// 2nd sub event plane
	mVecMesonTree->passEventPlane2East(Q2East);
	mVecMesonTree->passEventPlane2West(Q2West);

	// Number of Track in East and West part of TPC
	mVecMesonTree->passNumTrackEast(NumTrackEast);
	mVecMesonTree->passNumTrackWest(NumTrackWest);

	mVecMesonTree->MixEvent_Phi(mFlag_ME,mPicoDst,cent9,vz,Psi2_East);
      }
      mVecMesonCorrection->clear();
    }
  }

  return kStOK;
}

