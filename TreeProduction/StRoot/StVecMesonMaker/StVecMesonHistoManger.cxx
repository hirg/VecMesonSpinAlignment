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
  h_Full = new TH1F("h_Full","h_Full",360,-TMath::Pi(),TMath::Pi());
  for(Int_t i = 0; i < 4; ++i)
  {
    TString HistName;
    HistName = Form("h_East_%d",i);
    h_East[i] = new TH1F(HistName.Data(),HistName.Data(),360,-TMath::Pi(),TMath::Pi());
    HistName = Form("h_West_%d",i);
    h_West[i] = new TH1F(HistName.Data(),HistName.Data(),360,-TMath::Pi(),TMath::Pi());
  }
}
//-------------------------------------------------------------------------------------------
void StVecMesonHistoManger::FillQA_Detector(Float_t dEdx, Float_t Mass2, Float_t p)
{
  h_mDEdx->Fill(p,dEdx);
  h_mMass2->Fill(p,Mass2);
}

void StVecMesonHistoManger::FillEP_Eta(Float_t Psi2_East, Float_t Psi2_West, Int_t j)
{
  h_East[j]->Fill(Psi2_East);
  h_West[j]->Fill(Psi2_West);
}

void StVecMesonHistoManger::FillEP_Full(Float_t Psi2_Full)
{
  h_Full->Fill(Psi2);
}
//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
void StVecMesonHistoManger::WriteQA()
{
  h_mDEdx->Write();
  h_mMass2->Write();
  h_Full->Write();
  for(Int_t i = 0; i < 4; ++i)
  {
    h_East[i]->Write();
    h_West[i]->Write();
  }
}
//-------------------------------------------------------------------------------------------
