#include "StVecMesonHistoManger.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TMath.h"
#include "TString.h"

ClassImp(StVecMesonHistoManger)

//-------------------------------------------------------------------------------------------

StVecMesonHistoManger::StVecMesonHistoManger()
{
}

//-------------------------------------------------------------------------------------------

StVecMesonHistoManger::~StVecMesonHistoManger()
{
  /* */
}

//-------------------------------------------------------------------------------------------
void StVecMesonHistoManger::InitQA()
{
  h_mDEdx = new TH2F("h_mDEdx","h_mDEdx",1000,0,4.0,1000,0,40);
  h_mMass2 = new TH2F("h_mMass2","h_mMass2",1000,0,4.0,1000,-0.3,1.7);
  h_mFullRaw = new TH1F("h_mFullRaw","h_mFullRaw",360,-TMath::Pi(),TMath::Pi());
  h_mEastRaw = new TH1F("h_mEastRaw","h_mEastRaw",360,-TMath::Pi(),TMath::Pi());
  h_mWestRaw = new TH1F("h_mWestRaw","h_mWestRaw",360,-TMath::Pi(),TMath::Pi());
  h_mVz   = new TH1F("h_mVz","h_mVz",201,-100.5,100.5);
  h_mRefMult = new TH1F("h_mRefMult","h_mRefMult",1000,-0.5,999.5);
}

void StVecMesonHistoManger::FillQA_Detector(Float_t dEdx, Float_t Mass2, Float_t p)
{
  h_mDEdx->Fill(p,dEdx);
  h_mMass2->Fill(p,Mass2);
}

void StVecMesonHistoManger::FillQA_Event(Float_t vz, Float_t refMult)
{
  h_mVz->Fill(vz);
  h_mRefMult->Fill(refMult);
}

void StVecMesonHistoManger::FillEP_Eta(Float_t Psi2_East, Float_t Psi2_West)
{
  h_mEastRaw->Fill(Psi2_East);
  h_mWestRaw->Fill(Psi2_West);
}

void StVecMesonHistoManger::FillEP_Full(Float_t Psi2_Full)
{
  h_mFullRaw->Fill(Psi2_Full);
}

void StVecMesonHistoManger::WriteQA()
{
  h_mDEdx->Write();
  h_mMass2->Write();
  h_mVz->Write();
  h_mRefMult->Write();
  h_mFullRaw->Write();
  h_mEastRaw->Write();
  h_mWestRaw->Write();
}
//-------------------------------------------------------------------------------------------
void StVecMesonHistoManger::InitEP()
{
  h_mEastReCenter = new TH1F("h_mEastReCenter","h_mEastReCenter",360,-TMath::Pi(),TMath::Pi());
  h_mWestReCenter = new TH1F("h_mWestReCenter","h_mWestReCenter",360,-TMath::Pi(),TMath::Pi());
  h_mRanAReCenter = new TH1F("h_mRanAReCenter","h_mRanAReCenter",360,-TMath::Pi(),TMath::Pi());
  h_mRanBReCenter = new TH1F("h_mRanBReCenter","h_mRanBReCenter",360,-TMath::Pi(),TMath::Pi());
  h_mFullReCenter = new TH1F("h_mFullReCenter","h_mFullReCenter",360,-TMath::Pi(),TMath::Pi());

  h_mEastShift = new TH1F("h_mEastShift","h_mEastShift",360,-TMath::Pi(),TMath::Pi());
  h_mWestShift = new TH1F("h_mWestShift","h_mWestShift",360,-TMath::Pi(),TMath::Pi());
  h_mRanAShift = new TH1F("h_mRanAShift","h_mRanAShift",360,-TMath::Pi(),TMath::Pi());
  h_mRanBShift = new TH1F("h_mRanBShift","h_mRanBShift",360,-TMath::Pi(),TMath::Pi());
  h_mFullShift = new TH1F("h_mFullShift","h_mFullShift",360,-TMath::Pi(),TMath::Pi());
}

void StVecMesonHistoManger::FillEP_Sub(Float_t Psi2East_ReCenter, Float_t Psi2East_Shift, Float_t Psi2West_ReCenter, Float_t Psi2West_Shift)
{
  h_mEastReCenter->Fill(Psi2East_ReCenter);
  h_mEastShift->Fill(Psi2East_Shift);

  h_mWestReCenter->Fill(Psi2West_ReCenter);
  h_mWestShift->Fill(Psi2West_Shift);
}

void StVecMesonHistoManger::FillEP_Ran(Float_t Psi2RanA_ReCenter, Float_t Psi2RanA_Shift, Float_t Psi2RanB_ReCenter, Float_t Psi2RanB_Shift, Float_t Psi2Full_ReCenter, Float_t Psi2Full_Shift)
{
  h_mRanAReCenter->Fill(Psi2RanA_ReCenter);
  h_mRanAShift->Fill(Psi2RanA_Shift);

  h_mRanBReCenter->Fill(Psi2RanB_ReCenter);
  h_mRanBShift->Fill(Psi2RanB_Shift);

  h_mFullReCenter->Fill(Psi2Full_ReCenter);
  h_mFullShift->Fill(Psi2Full_Shift);
}

void StVecMesonHistoManger::WriteEP()
{
  h_mEastReCenter->Write();
  h_mWestReCenter->Write();
  h_mRanAReCenter->Write();
  h_mRanBReCenter->Write();
  h_mFullReCenter->Write();

  h_mEastShift->Write();
  h_mWestShift->Write();
  h_mRanAShift->Write();
  h_mRanBShift->Write();
  h_mFullShift->Write();
}
