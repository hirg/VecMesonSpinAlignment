#include "StTriFlowMaker.h"
#include "StTriFlowConstants.h"
#include "StTriFlowCut.h"
#include "StTriFlowProManger.h"
#include "StTriFlowCorrection.h"
#include "StTriFlowHistoManger.h"
#include "StTriFlowV0.h"
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

ClassImp(StTriFlowMaker)

StRefMultCorr* StTriFlowMaker::mRefMultCorr = NULL;
//-----------------------------------------------------------------------------
StTriFlowMaker::StTriFlowMaker(const char* name, StPicoDstMaker *picoMaker, const Int_t jobCounter, const Int_t Mode, const Int_t energy, const Int_t flag_ME, const Int_t flag_Embedding)
  : StMaker(name)
{
  mPicoDstMaker = picoMaker;
  mPicoDst = 0;
  mMode = Mode;
  mEnergy = energy;
  mFlag_ME = flag_ME;
  mFlag_Embedding = flag_Embedding;


  if(mMode == 0)
  {
    mOutPut_ReCenterPar = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/RecenterParameter/file_%s_ReCenterPar_",TriFlow::Energy[energy].Data(),TriFlow::Energy[energy].Data());
    mOutPut_ReCenterPar += jobCounter;
    mOutPut_ReCenterPar += ".root";
  }
  if(mMode == 1)
  {
    mOutPut_Corr_Shift = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/Correction/Shift/file_%s_Corr_Shift_",TriFlow::Energy[energy].Data(),TriFlow::Energy[energy].Data()); 
    mOutPut_Corr_Shift += jobCounter;
    mOutPut_Corr_Shift += ".root";

    mOutPut_Corr_ReCenter = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/Correction/ReCenter/file_%s_Corr_ReCenter_",TriFlow::Energy[energy].Data(),TriFlow::Energy[energy].Data()); 
    mOutPut_Corr_ReCenter += jobCounter;
    mOutPut_Corr_ReCenter += ".root";
  }
  if(mMode == 2)
  {
    if(mFlag_Embedding == 0) mOutPut_Phi = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/Phi/file_%s_Phi_%s_",TriFlow::Energy[energy].Data(),TriFlow::Energy[energy].Data(),TriFlow::MixEvent[mFlag_ME].Data()); 
    if(mFlag_Embedding == 1) mOutPut_Phi = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/Embedding/Phi/file_%s_Phi_%s_",TriFlow::Energy[energy].Data(),TriFlow::Energy[energy].Data(),TriFlow::MixEvent[mFlag_ME].Data()); 
    mOutPut_Phi += jobCounter;
    mOutPut_Phi += ".root";
  }
}

//----------------------------------------------------------------------------- 
StTriFlowMaker::~StTriFlowMaker()
{ /*  */ }

//----------------------------------------------------------------------------- 
Int_t StTriFlowMaker::Init() 
{
  if(!mRefMultCorr)
  {
    mRefMultCorr = CentralityMaker::instance()->getRefMultCorr();
  }
  mTriFlowCut = new StTriFlowCut(mEnergy);
  mTriFlowProManger = new StTriFlowProManger();
  mTriFlowCorrection = new StTriFlowCorrection(mEnergy);
  mTriFlowHistoManger = new StTriFlowHistoManger(mEnergy);

  if(mMode == 0)
  {
    mFile_ReCenterPar = new TFile(mOutPut_ReCenterPar.Data(),"RECREATE");
    mFile_ReCenterPar->cd();
    mTriFlowProManger->InitReCenter();
    mTriFlowHistoManger->InitQA_Detector();
  }

  if(mMode == 1)
  {
    mUsedTrackCounter = 0;
    mTriFlowCorrection->InitReCenterCorrection(mEnergy);
    mFile_Corr_Shift = new TFile(mOutPut_Corr_Shift.Data(),"RECREATE");
    mTriFlowProManger->InitShift();
    mFile_Corr_ReCenter = new TFile(mOutPut_Corr_ReCenter.Data(),"RECREATE");
    mFile_Corr_ReCenter->cd();
    mTriFlowCorrection->InitNtuple();
  }
  if(mMode == 2)
  {
    mTriFlowV0 = new StTriFlowV0(mEnergy);
    mFile_Phi = new TFile(mOutPut_Phi.Data(),"RECREATE");
    mFile_Phi->cd();
    mTriFlowV0->InitPhi();
    mTriFlowCorrection->InitReCenterCorrection(mEnergy);
    mTriFlowCorrection->InitShiftCorrection(mEnergy);
  }

  return kStOK;
}

//----------------------------------------------------------------------------- 
Int_t StTriFlowMaker::Finish() 
{
  if(mMode == 0)
  {
    if(mOutPut_ReCenterPar != "")
    {
      mFile_ReCenterPar->cd();
      mTriFlowProManger->WriteReCenter();
      mTriFlowHistoManger->WriteQA_Detector();
      mFile_ReCenterPar->Close();
    }
  }
  if(mMode == 1)
  {
    if(mOutPut_Corr_ReCenter != "")
    {
      mFile_Corr_ReCenter->cd();
      mTriFlowCorrection->writeNtuple();
      mFile_Corr_ReCenter->Close();
    }
    if(mOutPut_Corr_Shift != "")
    {
      mFile_Corr_Shift->cd();
      mTriFlowProManger->WriteShift();
      mFile_Corr_Shift->Close();
    }
  }
  if(mMode == 2)
  {
    if(mOutPut_Phi != "")
    {
      mFile_Phi->cd();
      mTriFlowV0->WritePhiMass2();
      mFile_Phi->Close();
    }
  }

  return kStOK;
}

//----------------------------------------------------------------------------- 
void StTriFlowMaker::Clear(Option_t *opt) {
}

//----------------------------------------------------------------------------- 
Int_t StTriFlowMaker::Make() 
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
//  cout << "StTriFLowMaker:" << endl;
//  cout << mRefMultCorr->getWeight() << endl;

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
  if(mTriFlowCut->passEventCut(mPicoDst)) // event cut
  {
    const Int_t nTracks = mPicoDst->numberOfTracks();
    const Int_t cent9 = mRefMultCorr->getCentralityBin9();
//    if(cent9 < 0) cout << cent9 << endl;
    const Double_t reweight = mRefMultCorr->getWeight();
    const Int_t nToFMatched = mTriFlowCut->getMatchedToF();

//    cout << "nTracks = " << nTracks << endl;
    for(Int_t i = 0; i < nTracks; i++) // track loop
    {
      StPicoTrack *track = (StPicoTrack*)mPicoDst->track(i);

      if(mMode == 0)
      {
	Float_t eta = track->pMom().pseudoRapidity();
	if(fabs(eta) < TriFlow::mEtaMax && track->dca() < 3.0)
	{
	  Float_t Mass2 = mTriFlowCut->getMass2(track);
	  Float_t dEdx = track->dEdx();
	  Float_t p = track->pMom().mag();
	  mTriFlowHistoManger->FillQA_Detector(dEdx,Mass2,p);
	}
      }
      if(mTriFlowCut->passTrackEP(track)) // track cut
      {
	if(mMode == 0) // fill re-center parameter
	{
	  Float_t pt = track->pMom().perp();

	  if(mTriFlowCorrection->passTrackFull(track))
	  {
	    TVector2 q2Vector_Full = mTriFlowCorrection->calq2Vector(track);
	    TVector2 q3Vector_Full = mTriFlowCorrection->calq3Vector(track);
	    mTriFlowProManger->FillTrackFull(q2Vector_Full,q3Vector_Full,cent9,runIndex,vz_sign,pt);
	  }

	  for(Int_t j = 0; j < 4; j++) // eta_gap loop
	  {
	    if(mTriFlowCorrection->passTrackEtaEast(track,j,0)) // neg eta sub
	    {
	      TVector2 q2Vector_East = mTriFlowCorrection->calq2Vector(track);
	      TVector2 q3Vector_East = mTriFlowCorrection->calq3Vector(track);
	      mTriFlowProManger->FillTrackEast(q2Vector_East,q3Vector_East,cent9,runIndex,vz_sign,j,pt);
	    }
	    if(mTriFlowCorrection->passTrackEtaWest(track,j,0)) // pos eta sub
	    {
	      TVector2 q2Vector_West = mTriFlowCorrection->calq2Vector(track);
	      TVector2 q3Vector_West = mTriFlowCorrection->calq3Vector(track);
	      mTriFlowProManger->FillTrackWest(q2Vector_West,q3Vector_West,cent9,runIndex,vz_sign,j,pt);
	    }
	  }
	}
	else // calculate Q Vector after recentering for full event and eta sub event
//	if(mMode == 1 || mMode == 2)
	{
	  if(mTriFlowCorrection->passTrackFull(track))
	  {
	    mTriFlowCorrection->addTrack_Full(track,cent9,runIndex,vz_sign);
	    mUsedTrackCounter++;
	  }
	  for(Int_t j = 0; j < 4; j++) // eta_gap loop
	  {
	    if(mTriFlowCorrection->passTrackEtaEast(track,j,0)) // neg eta sub
	    {
	      mTriFlowCorrection->addTrack_East(track,cent9,runIndex,vz_sign,j);
	    }
	    if(mTriFlowCorrection->passTrackEtaWest(track,j,0)) // pos eta sub
	    {
	      mTriFlowCorrection->addTrack_West(track,cent9,runIndex,vz_sign,j);
	    }
	  }
	}
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
      random_shuffle(iTrack,iTrack+mUsedTrackCounter);
      mUsedTrackCounter = 0;
      for(Int_t i = 0; i < nTracks; i++) // track loop
      {
	StPicoTrack *track = (StPicoTrack*)mPicoDst->track(i);
	if(mTriFlowCut->passTrackEP(track)) // track cut
	{
	  if(mTriFlowCorrection->passTrackFull(track))
	  {
	    if((Float_t)iTrack[mUsedTrackCounter] > ranCounter) // Sub Event A
	    {
	      mTriFlowCorrection->addTrack_A(track,cent9,runIndex,vz_sign);
	    }
	    else // Sub Event B
	    {
	      mTriFlowCorrection->addTrack_B(track,cent9,runIndex,vz_sign);
	    }
	    mUsedTrackCounter++;
	  }
	}
      }
    }
    if(mMode == 1) // re-center and calculate shift parameter for EP and SP
    {
      mTriFlowCorrection->fillNtuple(mPicoDst,cent9,nToFMatched,runIndex);

      // full event shift parameter
      if(mTriFlowCorrection->passTrackFullNumCut())
      {
	for(Int_t k = 0; k < 5; k++) // ShiftOrder loop
	{
	  // Event Plane method
	  TVector2 Psi2Vector_Full_EP = mTriFlowCorrection->calPsi2_Full_EP(k);
	  TVector2 Psi3Vector_Full_EP = mTriFlowCorrection->calPsi3_Full_EP(k);
	  mTriFlowProManger->FillEventFull_EP(Psi2Vector_Full_EP,Psi3Vector_Full_EP,cent9,runIndex,vz_sign,k);

	  // Scalor Product method
	  TVector2 Psi2Vector_Full_SP = mTriFlowCorrection->calPsi2_Full_SP(k);
	  TVector2 Psi3Vector_Full_SP = mTriFlowCorrection->calPsi3_Full_SP(k);
	  mTriFlowProManger->FillEventFull_SP(Psi2Vector_Full_SP,Psi3Vector_Full_SP,cent9,runIndex,vz_sign,k);
	}
      }

      // eta sub event shift parameter
      for(Int_t j = 0; j < 4; j ++)
      {
        if(mTriFlowCorrection->passTrackEtaNumCut(j))
	{
	  for(Int_t k = 0; k < 5; k++)
	  {
	    // Event Plane method
	    TVector2 Psi2Vector_East_EP = mTriFlowCorrection->calPsi2_East_EP(j,k);
	    TVector2 Psi3Vector_East_EP = mTriFlowCorrection->calPsi3_East_EP(j,k);
	    mTriFlowProManger->FillEventEast_EP(Psi2Vector_East_EP,Psi3Vector_East_EP,cent9,runIndex,vz_sign,j,k);

	    TVector2 Psi2Vector_West_EP = mTriFlowCorrection->calPsi2_West_EP(j,k);
	    TVector2 Psi3Vector_West_EP = mTriFlowCorrection->calPsi3_West_EP(j,k);
	    mTriFlowProManger->FillEventWest_EP(Psi2Vector_West_EP,Psi3Vector_West_EP,cent9,runIndex,vz_sign,j,k);

	    // Scalor Product method
	    TVector2 Psi2Vector_East_SP = mTriFlowCorrection->calPsi2_East_SP(j,k);
	    TVector2 Psi3Vector_East_SP = mTriFlowCorrection->calPsi3_East_SP(j,k);
	    mTriFlowProManger->FillEventEast_SP(Psi2Vector_East_SP,Psi3Vector_East_SP,cent9,runIndex,vz_sign,j,k);

	    TVector2 Psi2Vector_West_SP = mTriFlowCorrection->calPsi2_West_SP(j,k);
	    TVector2 Psi3Vector_West_SP = mTriFlowCorrection->calPsi3_West_SP(j,k);
	    mTriFlowProManger->FillEventWest_SP(Psi2Vector_West_SP,Psi3Vector_West_SP,cent9,runIndex,vz_sign,j,k);
	  }
	}
      }
      mTriFlowCorrection->clear();
      mUsedTrackCounter = 0;
    }

    if(mMode == 2)
    { // phi meson
      Float_t Psi2_East;
      TVector2 Q2East[TriFlow::EtaGap_total], Q2West[TriFlow::EtaGap_total];
      TVector2 Q3East[TriFlow::EtaGap_total], Q3West[TriFlow::EtaGap_total]; 
      Int_t NumTrackEast[TriFlow::EtaGap_total], NumTrackWest[TriFlow::EtaGap_total]; 
      for(Int_t j = 0; j < TriFlow::EtaGap_total; j++)
      {
	Q2East[j].Set(-999.9,-999.9); // initialize Q Vector to unreasonable value
	Q2West[j].Set(-999.9,-999.9);
	Q3East[j].Set(-999.9,-999.9);
	Q3West[j].Set(-999.9,-999.9);
	NumTrackEast[j] = 0;
	NumTrackWest[j] = 0;
	if(mTriFlowCorrection->passTrackEtaNumCut(j))
	{
	  // get QVector of sub event
	  Q2East[j] = mTriFlowCorrection->getQVector(j,0,0); // 0 = eta_gap, 1 = flow type, 2 = east/west
	  Q2West[j] = mTriFlowCorrection->getQVector(j,0,1); // 0 = eta_gap, 1 = flow type, 2 = east/west
	  Q3East[j] = mTriFlowCorrection->getQVector(j,1,0); // 0 = eta_gap, 1 = flow type, 2 = east/west
	  Q3West[j] = mTriFlowCorrection->getQVector(j,1,1); // 0 = eta_gap, 1 = flow type, 2 = east/west
	  NumTrackEast[j] = mTriFlowCorrection->getNumTrack(j,0); // 0 = eta_gap, 1 = east/west
	  NumTrackWest[j] = mTriFlowCorrection->getNumTrack(j,1); // 0 = eta_gap, 1 = east/west
	}
      }

      if(mTriFlowCorrection->passTrackEtaNumCut(0))
      {
	// Event Plane method
	Psi2_East = mTriFlowCorrection->calShiftAngle2East_EP(runIndex,cent9,vz_sign,0);

	// get N_prim, N_non_prim, N_Tof_match
	Int_t N_prim = mTriFlowCut->getNpirm();
	Int_t N_non_prim = mTriFlowCut->getNnonprim();
	Int_t N_Tof_match = mTriFlowCut->getMatchedToF();

	// pass the event information to StTriFlowV0
	mTriFlowV0->clearEvent();
	mTriFlowV0->passEvent(N_prim, N_non_prim, N_Tof_match);

	// 2nd sub event plane
	mTriFlowV0->passEventPlane2East(Q2East[0],Q2East[1],Q2East[2],Q2East[3]);
	mTriFlowV0->passEventPlane2West(Q2West[0],Q2West[1],Q2West[2],Q2West[3]);

	// 3rd sub event plane
	mTriFlowV0->passEventPlane3East(Q3East[0],Q3East[1],Q3East[2],Q3East[3]);
	mTriFlowV0->passEventPlane3West(Q3West[0],Q3West[1],Q3West[2],Q3West[3]);

	// Number of Track in East and West part of TPC
	mTriFlowV0->passNumTrackEast(NumTrackEast[0],NumTrackEast[1],NumTrackEast[2],NumTrackEast[3]);
	mTriFlowV0->passNumTrackWest(NumTrackWest[0],NumTrackWest[1],NumTrackWest[2],NumTrackWest[3]);

	mTriFlowV0->MixEvent_Phi(mFlag_ME,mPicoDst,cent9,vz,Psi2_East);
      }
      mTriFlowCorrection->clear();
    }
  }

  return kStOK;
}

