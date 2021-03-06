#include "StRoot/StZdcSmdMaker/StZdcSmdMaker.h"
#include "StRoot/StZdcSmdMaker/StZdcSmdCut.h"
#include "StRoot/StZdcSmdMaker/StZdcSmdProManger.h"
#include "StRoot/StZdcSmdMaker/StZdcSmdCorr.h"
#include "StRoot/StZdcSmdMaker/StZdcSmdHistoManger.h"
#include "StRoot/StZdcSmdMaker/StZdcSmdTree.h"
#include "../Utility/StSpinAlignmentCons.h"
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
#include "TVector2.h"
#include "TMath.h"
#include "StMessMgr.h"
#include <algorithm>

ClassImp(StZdcSmdMaker)

StRefMultCorr* StZdcSmdMaker::mRefMultCorr = NULL;
//-----------------------------------------------------------------------------
StZdcSmdMaker::StZdcSmdMaker(const char* name, StPicoDstMaker *picoMaker, const Int_t jobCounter, const Int_t Mode, const Int_t energy, const Int_t flag_ME)
  : StMaker(name)
{
  mPicoDstMaker = picoMaker;
  mPicoDst = 0;
  mMode = Mode;
  mEnergy = energy;
  mFlag_ME = flag_ME;
  mStopWatch = new TStopwatch();

  if(mMode == 0)
  { // fill zdc-smd QA and gain correction factor
    mOutPut_GainCorrPar = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/ZDCSMD/GainCorrPar/file_%s_GainCorrPar_%d.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str(),jobCounter);
  }
  if(mMode == 1)
  { // apply gain correction and fill re-center parameter
    mOutPut_ReCenterPar = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/ZDCSMD/ReCenterPar/file_%s_ReCenterPar_%d.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str(),jobCounter);
  }
  if(mMode == 2)
  { // apply gian and re-center correction and fill shift parameter for East/West
    mOutPut_ShiftPar = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/ZDCSMD/ShiftPar/file_%s_ShiftPar_%d.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str(),jobCounter);
  }
  if(mMode == 3)
  { // apply gian, re-center and shift correction and fill shift parameter for Full 
    mOutPut_ShiftParFull = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/ZDCSMD/ShiftPar/file_%s_ShiftParFull_%d.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str(),jobCounter);
  }
  if(mMode == 4)
  { // apply gian, re-center and shift correction and fill shift parameter for Full 
    mOutPut_Resolution = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/ZDCSMD/Resolution/file_%s_Resolution_%d.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str(),jobCounter);
  }
  if(mMode == 5)
  { // apply gian, re-center and shift correction and fill shift parameter for Full 
    mOutPut_DirectedFlow = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/ZDCSMD/DirectedFlow/file_%s_DirectedFlow_%d.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str(),jobCounter);
  }
  if(mMode == 6)
  { // apply gian, re-center and shift correction and fill shift parameter for Full 
    mOutPut_Phi = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/ZDCSMD/Phi/Forest/file_%s_Phi_%s_%d.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::MixEvent[mFlag_ME].Data(),jobCounter);
  }
}

//----------------------------------------------------------------------------- 
StZdcSmdMaker::~StZdcSmdMaker()
{ /*  */ }

//----------------------------------------------------------------------------- 
Int_t StZdcSmdMaker::Init() 
{
  if(!mRefMultCorr)
  {
    mRefMultCorr = CentralityMaker::instance()->getRefMultCorr();
  }
  mZdcSmdCut = new StZdcSmdCut(mEnergy);
  mZdcSmdCorrection = new StZdcSmdCorrection(mEnergy);
  mZdcSmdHistoManger = new StZdcSmdHistoManger();
  mZdcSmdProManger = new StZdcSmdProManger();
  mStopWatch->Start();

  if(mMode == 0)
  {
    mFile_GainCorrPar = new TFile(mOutPut_GainCorrPar.Data(),"RECREATE");
    mFile_GainCorrPar->cd();
    mZdcSmdHistoManger->InitQA();
    mZdcSmdHistoManger->InitGainCorr();
  }
  if(mMode == 1)
  {
    mFile_ReCenterPar= new TFile(mOutPut_ReCenterPar.Data(),"RECREATE");
    mZdcSmdProManger->InitReCenter();
    mZdcSmdCorrection->ReadGainCorr();
    mZdcSmdHistoManger->InitRawEP();
    mFile_ReCenterPar->cd();
  }
  if(mMode == 2)
  {
    mFile_ShiftPar = new TFile(mOutPut_ShiftPar.Data(),"RECREATE");
    mZdcSmdProManger->InitShift();
    mZdcSmdCorrection->ReadGainCorr();
    mZdcSmdCorrection->ReadReCenterCorr();
    mZdcSmdHistoManger->InitReCenterEP();
    mFile_ShiftPar->cd();
  }
  if(mMode == 3)
  {
    mFile_ShiftParFull = new TFile(mOutPut_ShiftParFull.Data(),"RECREATE");
    mZdcSmdProManger->InitShiftFull();
    mZdcSmdCorrection->ReadGainCorr();
    mZdcSmdCorrection->ReadReCenterCorr();
    mZdcSmdCorrection->ReadShiftCorr();
    mZdcSmdHistoManger->InitShiftEP();
    mFile_ShiftParFull->cd();
  }
  if(mMode == 4)
  {
    mFile_Resolution = new TFile(mOutPut_Resolution.Data(),"RECREATE");
    mZdcSmdProManger->InitResolution();
    mZdcSmdCorrection->ReadGainCorr();
    mZdcSmdCorrection->ReadReCenterCorr();
    mZdcSmdCorrection->ReadShiftCorr();
    mZdcSmdCorrection->ReadShiftCorrFull();
    mZdcSmdHistoManger->InitShiftEPFull();
    mFile_Resolution->cd();
  }
  if(mMode == 5)
  {
    mFile_DirectedFlow = new TFile(mOutPut_DirectedFlow.Data(),"RECREATE");
    mZdcSmdCorrection->ReadGainCorr();
    mZdcSmdCorrection->ReadReCenterCorr();
    mZdcSmdCorrection->ReadShiftCorr();
    mZdcSmdCorrection->ReadShiftCorrFull();
    mZdcSmdCorrection->ReadResolution();
    mZdcSmdCorrection->CalResolution(); // get full event plane resolution
    mZdcSmdProManger->InitDirectedFlow();
    mFile_DirectedFlow->cd();
  }
  if(mMode == 6)
  {
    mZdcSmdTree = new StZdcSmdTree(mEnergy);
    mFile_Phi = new TFile(mOutPut_Phi.Data(),"RECREATE");
    mZdcSmdCorrection->ReadGainCorr();
    mZdcSmdCorrection->ReadReCenterCorr();
    mZdcSmdCorrection->ReadShiftCorr();
    mZdcSmdCorrection->ReadShiftCorrFull();
    mFile_Phi->cd();
    mZdcSmdTree->InitPhi();
  }

  return kStOK;
}

//----------------------------------------------------------------------------- 
Int_t StZdcSmdMaker::Finish() 
{
  if(mMode == 0)
  {
    if(mOutPut_GainCorrPar != "")
    {
      mFile_GainCorrPar->cd();
      mZdcSmdHistoManger->WriteQA();
      mZdcSmdHistoManger->WriteGainCorr();
      mFile_GainCorrPar->Close();
    }
  }
  if(mMode == 1)
  {
    if(mOutPut_ReCenterPar != "")
    {
      mFile_ReCenterPar->cd();
      mZdcSmdProManger->WriteReCenter();
      mZdcSmdHistoManger->WriteRawEP();
      mFile_ReCenterPar->Close();
    }
  }
  if(mMode == 2)
  {
    if(mOutPut_ShiftPar != "")
    {
      mFile_ShiftPar->cd();
      mZdcSmdProManger->WriteShift();
      mZdcSmdHistoManger->WriteReCenterEP();
      mFile_ShiftPar->Close();
    }
  }
  if(mMode == 3)
  {
    if(mOutPut_ShiftParFull != "")
    {
      mFile_ShiftParFull->cd();
      mZdcSmdProManger->WriteShiftFull();
      mZdcSmdHistoManger->WriteShiftEP();
      mFile_ShiftParFull->Close();
    }
  }
  if(mMode == 4)
  {
    if(mOutPut_Resolution != "")
    {
      mFile_Resolution->cd();
      mZdcSmdProManger->WriteResolution();
      mZdcSmdHistoManger->WriteShiftEPFull();
      mFile_Resolution->Close();
    }
  }
  if(mMode == 5)
  {
    if(mOutPut_DirectedFlow != "")
    {
      mFile_DirectedFlow->cd();
      mZdcSmdProManger->WriteDirectedFlow();
      mFile_DirectedFlow->Close();
    }
  }
  if(mMode == 6)
  {
    if(mOutPut_Phi != "")
    {
      mFile_Phi->cd();
      mZdcSmdTree->WritePhiMass2();
      mFile_Phi->Close();
    }
  }

  mStopWatch->Stop();
  mStopWatch->Print();

  return kStOK;
}

//----------------------------------------------------------------------------- 
void StZdcSmdMaker::Clear(Option_t *opt) {
}

//----------------------------------------------------------------------------- 
Int_t StZdcSmdMaker::Make() 
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
  if(mEnergy == 6) mRefMultCorr->initEvent(refMult,vz,zdcX); // for 200 GeV
  if(mEnergy != 6) mRefMultCorr->initEvent(refMult,vz); // for BES Energy

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

  if(mZdcSmdCut->passEventCut(mPicoDst)) // event cut
  {
    const unsigned int nTracks = mPicoDst->numberOfTracks();
    const int cent9 = mRefMultCorr->getCentralityBin9();
    // if(cent9 < 0) cout << cent9 << endl;
    const float reweight = mRefMultCorr->getWeight();
    const int nToFMatched = mZdcSmdCut->getMatchedToF();

    mZdcSmdCorrection->InitEvent(cent9,runIndex,vz_sign);

    // set ADC for each slats
    if(mMode == 0)
    {
      for(int i_slat = 0; i_slat < 8; ++i_slat) // read in raw ADC value from ZDC-SMD
      {
	mZdcSmdCorrection->SetZdcSmd(0,0,i_slat,mPicoEvent->ZdcSmdEastVertical(i_slat));
	mZdcSmdCorrection->SetZdcSmd(0,1,i_slat,mPicoEvent->ZdcSmdEastHorizontal(i_slat));
	mZdcSmdCorrection->SetZdcSmd(1,0,i_slat,mPicoEvent->ZdcSmdWestVertical(i_slat));
	mZdcSmdCorrection->SetZdcSmd(1,1,i_slat,mPicoEvent->ZdcSmdWestHorizontal(i_slat));
      }
    }
    if(mMode > 0)
    {
      for(int i_slat = 0; i_slat < 8; ++i_slat) // read in raw ADC value from ZDC-SMD
      {
	mZdcSmdCorrection->SetZdcSmdGainCorr(0,0,i_slat,mPicoEvent->ZdcSmdEastVertical(i_slat));
	mZdcSmdCorrection->SetZdcSmdGainCorr(0,1,i_slat,mPicoEvent->ZdcSmdEastHorizontal(i_slat));
	mZdcSmdCorrection->SetZdcSmdGainCorr(1,0,i_slat,mPicoEvent->ZdcSmdWestVertical(i_slat));
	mZdcSmdCorrection->SetZdcSmdGainCorr(1,1,i_slat,mPicoEvent->ZdcSmdWestHorizontal(i_slat));
      }
    }

    if(mMode == 0) // fill zdc-smd QA and gain correction fator
    {
      mZdcSmdHistoManger->FillQA_Event(vz,refMult);

      for(int i_eastwest = 0; i_eastwest < 2; ++i_eastwest)
      {
	for(int i_verthori = 0; i_verthori < 2; ++i_verthori)
	{
	  for(int i_slat = 0; i_slat < 8; ++i_slat)
	  {
	    mZdcSmdHistoManger->FillGainCorr(i_eastwest,i_verthori,i_slat,runIndex,mZdcSmdCorrection->GetZdcSmd(i_eastwest,i_verthori,i_slat));
	    // cout << "i_eastwest = " << i_eastwest << ", i_verthori = " << i_verthori << ", i_slat = " << i_slat << ", zdc = " << mZdcSmdCorrection->GetZdcSmd(i_eastwest,i_verthori,i_slat) << endl;
	  }
	}
      }
    }
    if(mMode == 1) // apply gain correction and fill recenter correction parameter
    {
      TVector2 QEast = mZdcSmdCorrection->GetQEast(mMode);
      TVector2 QWest = mZdcSmdCorrection->GetQWest(mMode);
      TVector2 QFull = QWest-QEast;
      if( !(QEast.Mod() < 1e-10 || QWest.Mod() < 1e-10 || QFull.Mod() < 1e-10) )
      {
	mZdcSmdProManger->FillReCenterEast(QEast,cent9,runIndex,vz_sign);
	mZdcSmdProManger->FillReCenterWest(QWest,cent9,runIndex,vz_sign);
	mZdcSmdHistoManger->FillRawEP(QEast,QWest,QFull,cent9,runIndex);
      }
    }
    if(mMode == 2) // apply gain and re-center correction and fill shift correction parameter
    {
      TVector2 QEast = mZdcSmdCorrection->GetQEast(mMode);
      TVector2 QWest = mZdcSmdCorrection->GetQWest(mMode);
      TVector2 QFull = QWest-QEast;
      if( !(QEast.Mod() < 1e-10 || QWest.Mod() < 1e-10 || QFull.Mod() < 1e-10) )
      {
	mZdcSmdProManger->FillShiftEast(QEast,cent9,runIndex,vz_sign);
	mZdcSmdProManger->FillShiftWest(QWest,cent9,runIndex,vz_sign);
	mZdcSmdHistoManger->FillReCenterEP(QEast,QWest,QFull,cent9,runIndex);
      }
    }
    if(mMode == 3) // apply gain, re-center and shift correction and fill shift correction parameter for full event plane
    {
      TVector2 QEast = mZdcSmdCorrection->GetQEast(mMode);
      TVector2 QWest = mZdcSmdCorrection->GetQWest(mMode);
      TVector2 QFull = QWest-QEast;
      if( !(QEast.Mod() < 1e-10 || QWest.Mod() < 1e-10 || QFull.Mod() < 1e-10) )
      {
	mZdcSmdProManger->FillShiftFull(QFull,cent9,runIndex,vz_sign);
	mZdcSmdHistoManger->FillShiftEP(QEast,QWest,QFull,cent9,runIndex);
      }
    }
    if(mMode == 4) // apply gain, re-center and shift correction and fill event plane resolution for East/West
    {
      TVector2 QEast = mZdcSmdCorrection->GetQEast(mMode);
      TVector2 QWest = mZdcSmdCorrection->GetQWest(mMode);
      TVector2 QFull = mZdcSmdCorrection->GetQFull(QEast,QWest);
      if( !(QEast.Mod() < 1e-10 || QWest.Mod() < 1e-10 || QFull.Mod() < 1e-10) )
      {
	mZdcSmdProManger->FillResolution(QEast,QWest,cent9);
	mZdcSmdHistoManger->FillShiftEPFull(QFull,cent9,runIndex);
      }
    }
    if(mMode == 5) // calculate v1 vs. eta for charged hadrons
    {
      TVector2 QEast = mZdcSmdCorrection->GetQEast(mMode);
      TVector2 QWest = mZdcSmdCorrection->GetQWest(mMode);
      TVector2 QFull = mZdcSmdCorrection->GetQFull(QEast,QWest);
      if(QEast.Mod() < 1e-10 || QWest.Mod() < 1e-10 || QFull.Mod() < 1e-10) return kStOK;
      float resolution = mZdcSmdCorrection->GetResolution(cent9);
      float Psi = TMath::ATan2(QFull.Y(),QFull.X());
      for(unsigned int i_track = 0; i_track < nTracks; ++i_track) // track loop
      {
	StPicoTrack *track = (StPicoTrack*)mPicoDst->track(i_track);
	if(mZdcSmdCut->passTrackV1(track)) // track cut
	{
	  float pt  = track->pMom().perp();
	  float phi = track->pMom().phi();
	  float eta = track->pMom().pseudoRapidity();
	  float v1 = TMath::Cos(phi-Psi);
	  mZdcSmdProManger->FillDirectedFlow(cent9,eta,pt,v1,resolution,reweight);
	}
      }
    }
    if(mMode == 6) // calculate v1 vs. eta for charged hadrons
    {
      TVector2 QEast = mZdcSmdCorrection->GetQEast(mMode);
      TVector2 QWest = mZdcSmdCorrection->GetQWest(mMode);
      TVector2 QFull = mZdcSmdCorrection->GetQFull(QEast,QWest);
      if(QEast.Mod() < 1e-10 || QWest.Mod() < 1e-10 || QFull.Mod() < 1e-10) return kStOK;

      // get N_prim, N_non_prim, N_Tof_match
      int N_prim = mZdcSmdCut->getNpirm();
      int N_non_prim = mZdcSmdCut->getNnonprim();
      int N_Tof_match = mZdcSmdCut->getMatchedToF();

      // pass the event information to StZdcSmdTree
      mZdcSmdTree->clearEvent();
      mZdcSmdTree->passEvent(N_prim,N_non_prim,N_Tof_match);

      // pass shifted event plane to StZdcSmdTree
      mZdcSmdTree->passEventPlane(QEast,QWest,QFull);

      float Psi = TMath::ATan2(QFull.Y(),QFull.X());
      mZdcSmdTree->MixEvent_Phi(mFlag_ME,mPicoDst,cent9,vz,Psi);
    }

    mZdcSmdCorrection->clear();
  }

  return kStOK;
}

