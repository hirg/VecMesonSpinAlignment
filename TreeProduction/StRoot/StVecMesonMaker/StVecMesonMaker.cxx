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
    mOutPut_ShiftPar = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/ShiftParameter/file_%s_Corr_Shift_%d.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str(),jobCounter); 
  }
  if(mMode == 2)
  {
    mOutPut_Resolution = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/Resolution/file_%s_Resolution_%d.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str(),jobCounter); 
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
    mFile_ShiftPar = new TFile(mOutPut_ShiftPar.Data(),"RECREATE");
    mVecMesonProManger->InitShift();
  }

  if(mMode == 2)
  {
    mVecMesonCorrection->InitReCenterCorrection(mEnergy);
    mVecMesonCorrection->InitShiftCorrection(mEnergy);
    mVecMesonProManger->InitResolution();
    mVecMesonHistoManger->InitEP();
    mFile_Resolution = new TFile(mOutPut_Resolution.Data(),"RECREATE");
  }

  if(mMode == 3)
  {
    mVecMesonTree = new StVecMesonTree(mEnergy);
    mFile_Phi = new TFile(mOutPut_Phi.Data(),"RECREATE");
    mFile_Phi->cd();
    mVecMesonTree->InitPhi();
    mVecMesonCorrection->InitReCenterCorrection(mEnergy);
    mVecMesonCorrection->InitShiftCorrection(mEnergy);
    mVecMesonCorrection->InitResolutionCorr(mEnergy);
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
    if(mOutPut_ShiftPar != "")
    {
      mFile_ShiftPar->cd();
      mVecMesonProManger->WriteShift();
      mFile_ShiftPar->Close();
    }
  }
  if(mMode == 2)
  {
    if(mOutPut_Resolution != "")
    {
      mFile_Resolution->cd();
      mVecMesonHistoManger->WriteEP();
      mVecMesonProManger->WriteResolution();
      mFile_Resolution->Close();
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
      mVecMesonHistoManger->FillQA_Event(vz,refMult);
      if(mVecMesonCorrection->passTrackEtaNumRawCut())
      {
	TVector2 Q2East = mVecMesonCorrection->getQVectorRaw(0); // 0 = eta_gap, 1 = east/west
	Float_t Psi2_East = 0.5*TMath::ATan2(Q2East.Y(),Q2East.X());
	TVector2 Q2West = mVecMesonCorrection->getQVectorRaw(1); // 0 = eta_gap, 1 = east/west
	Float_t Psi2_West = 0.5*TMath::ATan2(Q2West.Y(),Q2West.X());
	mVecMesonHistoManger->FillEP_Eta(Psi2_East,Psi2_West);
      }
      if(mVecMesonCorrection->passTrackFullNumRawCut())
      {
	TVector2 Q2Full = mVecMesonCorrection->getQVectorRaw(2);
	Float_t Psi2_Full = 0.5*TMath::ATan2(Q2Full.Y(),Q2Full.X());
	mVecMesonHistoManger->FillEP_Full(Psi2_Full);
      }
    }

    if(mMode == 1)
    {
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
      if(mVecMesonCorrection->passTrackEtaNumCut())
      {
	for(Int_t k = 0; k < 5; k++)
	{
	  TVector2 Psi2Vector_East_EP = mVecMesonCorrection->calPsi2_East_EP(k);
	  mVecMesonProManger->FillEventEast_EP(Psi2Vector_East_EP,cent9,runIndex,vz_sign,k);

	  TVector2 Psi2Vector_West_EP = mVecMesonCorrection->calPsi2_West_EP(k);
	  mVecMesonProManger->FillEventWest_EP(Psi2Vector_West_EP,cent9,runIndex,vz_sign,k);
	}
      }
    }

    if(mMode == 2) // calculate resolution for eta_sub and random sub event plane
    {
      // calculate Q vector after recentering for Random Sub Event
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
      mUsedTrackCounter = 0;

      // calculate resolution
      TVector2 QVecEast = mVecMesonCorrection->getQVector(0);
      Float_t Psi2East_ReCenter = 0.5*TMath::ATan2(QVecEast.Y(),QVecEast.X());
      Float_t Psi2East_Shift = mVecMesonCorrection->calShiftAngle2East_EP(runIndex,cent9,vz_sign);

      TVector2 QVecWest = mVecMesonCorrection->getQVector(1);
      Float_t Psi2West_ReCenter = 0.5*TMath::ATan2(QVecWest.Y(),QVecWest.X());
      Float_t Psi2West_Shift = mVecMesonCorrection->calShiftAngle2West_EP(runIndex,cent9,vz_sign);

      TVector2 QVecFull = mVecMesonCorrection->getQVector(2);
      Float_t Psi2Full_ReCenter = 0.5*TMath::ATan2(QVecFull.Y(),QVecFull.X());
      Float_t Psi2Full_Shift = mVecMesonCorrection->calShiftAngle2Full_EP(runIndex,cent9,vz_sign);

      TVector2 QVecRanA = mVecMesonCorrection->getQVector(3);
      Float_t Psi2RanA_ReCenter = 0.5*TMath::ATan2(QVecRanA.Y(),QVecRanA.X());
      Float_t Psi2RanA_Shift = mVecMesonCorrection->calShiftAngle2A_EP(runIndex,cent9,vz_sign);

      TVector2 QVecRanB = mVecMesonCorrection->getQVector(4);
      Float_t Psi2RanB_ReCenter = 0.5*TMath::ATan2(QVecRanB.Y(),QVecRanB.X());
      Float_t Psi2RanB_Shift = mVecMesonCorrection->calShiftAngle2B_EP(runIndex,cent9,vz_sign);

      if(mVecMesonCorrection->passTrackEtaNumCut())
      {
	mVecMesonHistoManger->FillEP_Sub(Psi2East_ReCenter,Psi2West_ReCenter,Psi2RanA_ReCenter,Psi2RanB_ReCenter,Psi2Full_ReCenter);
	mVecMesonProManger->FillRes_Sub(cent9,Psi2East_Shift,Psi2West_Shift);
      }

      if(mVecMesonProManger->passTrackFullNumCut())
      {
	mVecMesonHistoManger->FillEP_Ran(Psi2East_Shift,Psi2West_Shift,Psi2RanA_Shift,Psi2RanB_Shift,Psi2Full_Shift);
	mVecMesonProManger->FillRes_Ran(cent9,Psi2RanA_Shift,Psi2RanB_Shift);
      }
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
    }
    mVecMesonCorrection->clear();
  }

  return kStOK;
}

