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
  h_mFull = new TH1F("h_mFull","h_mFull",360,-TMath::Pi(),TMath::Pi());
  h_mEast = new TH1F("h_mEast","h_mEast",360,-TMath::Pi(),TMath::Pi());
  h_mWest = new TH1F("h_mWest","h_mWest",360,-TMath::Pi(),TMath::Pi());
  h_mVz   = new TH1F("h_mVz","h_mVz",201,-100.5,100.5);
  h_mRefMult = new TH1F("h_mRefMult","h_mRefMult",1000,-0.5,999.5);
}
//-------------------------------------------------------------------------------------------
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
  h_East->Fill(Psi2_East);
  h_West->Fill(Psi2_West);
}

void StVecMesonHistoManger::FillEP_Full(Float_t Psi2_Full)
{
  h_Full->Fill(Psi2_Full);
}
//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
void StVecMesonHistoManger::WriteQA()
{
  h_mDEdx->Write();
  h_mMass2->Write();
  h_mVz->Write();
  h_mRefMult->Write();
  h_Full->Write();
  h_East->Write();
  h_West->Write();
}
//-------------------------------------------------------------------------------------------
