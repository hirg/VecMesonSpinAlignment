#include "StTriFlowHistoManger.h"
#include "TH2F.h"

ClassImp(StTriFlowHistoManger)

//-------------------------------------------------------------------------------------------

StTriFlowHistoManger::StTriFlowHistoManger(Int_t energy)
{
  mEnergy = energy;
}

//-------------------------------------------------------------------------------------------

StTriFlowHistoManger::~StTriFlowHistoManger()
{
  /* */
}

//-------------------------------------------------------------------------------------------
void StTriFlowHistoManger::InitQA_Detector()
{
  h_mDEdx = new TH2F("h_mDEdx","h_mDEdx",1000,0,4.0,1000,0,40);
  h_mMass2 = new TH2F("h_mMass2","h_mMass2",1000,0,4.0,1000,-0.3,1.7);
}
//-------------------------------------------------------------------------------------------
void StTriFlowHistoManger::FillQA_Detector(Float_t dEdx, Float_t Mass2, Float_t p)
{
  h_mDEdx->Fill(p,dEdx);
  h_mMass2->Fill(p,Mass2);
}
//-------------------------------------------------------------------------------------------
void StTriFlowHistoManger::WriteQA_Detector()
{
  h_mDEdx->Write();
  h_mMass2->Write();
}
//-------------------------------------------------------------------------------------------
