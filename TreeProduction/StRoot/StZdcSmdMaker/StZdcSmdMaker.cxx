#include "StRoot/StZdcSmdMaker/StZdcSmdMaker.h"
#include "StRoot/StZdcSmdMaker/StZdcSmdCut.h"
// #include "StRoot/StZdcSmdMaker/StZdcSmdProManger.h"
#include "StRoot/StZdcSmdMaker/StZdcSmdCorr.h"
#include "StRoot/StZdcSmdMaker/StZdcSmdHistoManger.h"
// #include "StRoot/StZdcSmdMaker/StZdcSmdTree.h"
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

  if(mMode == 0)
  { // fill zdc-smd QA and gain correction factor
    mOutPut_GainCorrPar = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/ZDCSMD/GainCorrPar/file_%s_GainCorrPar_%d.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str(),jobCounter);
  }
  if(mMode == 1)
  { // fill zdc-smd QA and gain correction factor
    mOutPut_ReCenterPar = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/ZDCSMD/ReCenterPar/file_%s_ReCenterPar_%d.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str(),jobCounter);
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
    mZdcSmdCorrection->InitGainCorr();
    mFile_ReCenterPar->cd();
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
      mFile_ReCenterPar->Close();
    }
  }

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
    const Int_t nTracks = mPicoDst->numberOfTracks();
    const Int_t cent9 = mRefMultCorr->getCentralityBin9();
    //    if(cent9 < 0) cout << cent9 << endl;
    const Double_t reweight = mRefMultCorr->getWeight();
    const Int_t nToFMatched = mZdcSmdCut->getMatchedToF();

    if(mMode == 0) // fill zdc-smd QA and gain correction fator
    {
      mZdcSmdHistoManger->FillQA_Event(vz,refMult);
      for(int i_slat = 0; i_slat < 8; ++i_slat) // read in raw ADC value from ZDC-SMD
      {
	mZdcSmdCorrection->SetZdcSmd(0,0,i_slat,mPicoEvent->ZdcSmdEastVertical(i_slat));
	mZdcSmdCorrection->SetZdcSmd(0,1,i_slat,mPicoEvent->ZdcSmdEastHorizontal(i_slat));
	mZdcSmdCorrection->SetZdcSmd(1,0,i_slat,mPicoEvent->ZdcSmdWestVertical(i_slat));
	mZdcSmdCorrection->SetZdcSmd(1,1,i_slat,mPicoEvent->ZdcSmdWestHorizontal(i_slat));
      }

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
      for(int i_slat = 0; i_slat < 8; ++i_slat) // read in raw ADC value from ZDC-SMD
      {
	mZdcSmdCorrection->SetZdcSmdGainCorr(0,0,i_slat,mPicoEvent->ZdcSmdEastVertical(i_slat));
	mZdcSmdCorrection->SetZdcSmdGainCorr(0,1,i_slat,mPicoEvent->ZdcSmdEastHorizontal(i_slat));
	mZdcSmdCorrection->SetZdcSmdGainCorr(1,0,i_slat,mPicoEvent->ZdcSmdWestVertical(i_slat));
	mZdcSmdCorrection->SetZdcSmdGainCorr(1,1,i_slat,mPicoEvent->ZdcSmdWestHorizontal(i_slat));
      }
      mZdcSmdProManger->FillReCenterEast(mZdcSmdCorrection->GetQEast(mMode),cent9,runIndex,vz_sign);
      mZdcSmdProManger->FillReCenterWest(mZdcSmdCorrection->GetQWest(mMode),cent9,runIndex,vz_sign);
    }

    mZdcSmdCorrection->clear();
  }

  return kStOK;
}

