#include "StRoot/StZdcSmdMaker/StZdcSmdHistoManger.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TMath.h"
#include "TString.h"
#include "TMath.h"

ClassImp(StZdcSmdHistoManger)

//-------------------------------------------------------------------------------------------

StZdcSmdHistoManger::StZdcSmdHistoManger()
{
}

//-------------------------------------------------------------------------------------------

StZdcSmdHistoManger::~StZdcSmdHistoManger()
{
  /* */
}

//-------------------------------------------------------------------------------------------
void StZdcSmdHistoManger::InitQA()
{
  h_mVz   = new TH1F("h_mVz","h_mVz",201,-100.5,100.5);
  h_mRefMult = new TH1F("h_mRefMult","h_mRefMult",1000,-0.5,999.5);
}

void StZdcSmdHistoManger::FillQA_Event(float vz, float refMult)
{
  h_mVz->Fill(vz);
  h_mRefMult->Fill(refMult);
}

void StZdcSmdHistoManger::WriteQA()
{
  h_mVz->Write();
  h_mRefMult->Write();
}

//-------------------------------------------------------------------------------------------

void StZdcSmdHistoManger::InitGainCorr()
{
  for(int i_eastwest = 0; i_eastwest < 2; ++i_eastwest)
  {
    for(int i_verthori = 0; i_verthori < 2; ++i_verthori)
    {
      for(int i_slat = 0; i_slat < 8; ++i_slat)
      {
	string HistName = Form("h_mGainCorr%s%s_%d",vmsa::mEastWest[i_eastwest].c_str(),vmsa::mVertHori[i_verthori].c_str(),i_slat);
	h_mGainCorr[i_eastwest][i_verthori][i_slat] = new TH2F(HistName.c_str(),HistName.c_str(),2000,-0.5,1999.5,1000,-5.0,995.0);
      }
    }
  }
}

void StZdcSmdHistoManger::FillGainCorr(int i_eastwest, int i_verthori, int i_slat, int runIndex, float zdcsmd)
{
  h_mGainCorr[i_eastwest][i_verthori][i_slat]->Fill((float)runIndex,zdcsmd);
}

void StZdcSmdHistoManger::WriteGainCorr()
{
  for(int i_eastwest = 0; i_eastwest < 2; ++i_eastwest)
  {
    for(int i_verthori = 0; i_verthori < 2; ++i_verthori)
    {
      for(int i_slat = 0; i_slat < 8; ++i_slat)
      {
	h_mGainCorr[i_eastwest][i_verthori][i_slat]->Write();
      }
    }
  }
}

//-------------------------------------------------------------------------------------------

void StZdcSmdHistoManger::InitRawEP()
{
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    string HistName = Form("h_mRawEast_%d",i_cent);
    h_mRawEast[i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),360,-1.0*TMath::Pi(),TMath::Pi(),1600,-0.5,1599.5);
    HistName = Form("h_mRawWest_%d",i_cent);
    h_mRawWest[i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),360,-1.0*TMath::Pi(),TMath::Pi(),1600,-0.5,1599.5);
    HistName = Form("h_mRawFull_%d",i_cent);
    h_mRawFull[i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),360,-1.0*TMath::Pi(),TMath::Pi(),1600,-0.5,1599.5);
  }
}

void StZdcSmdHistoManger::FillRawEP(TVector2 QEast, TVector2 QWest, TVector2 QFull, int Cent9, int runIndex)
{
  float PsiEast = TMath::ATan2(QEast.Y(),QEast.X()); h_mRawEast[Cent9]->Fill(PsiEast,runIndex);
  float PsiWest = TMath::ATan2(QWest.Y(),QWest.X()); h_mRawWest[Cent9]->Fill(PsiWest,runIndex);
  float PsiFull = TMath::ATan2(QFull.Y(),QFull.X()); h_mRawFull[Cent9]->Fill(PsiFull,runIndex);
}

void StZdcSmdHistoManger::WriteRawEP()
{
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    h_mRawEast[i_cent]->Write();
    h_mRawWest[i_cent]->Write();
    h_mRawFull[i_cent]->Write();
  }
}

//-------------------------------------------------------------------------------------------

void StZdcSmdHistoManger::InitReCenterEP()
{
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    string HistName = Form("h_mReCenterEast_%d",i_cent);
    h_mReCenterEast[i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),360,-1.0*TMath::Pi(),TMath::Pi(),1600,-0.5,1599.5);
    HistName = Form("h_mReCenterWest_%d",i_cent);
    h_mReCenterWest[i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),360,-1.0*TMath::Pi(),TMath::Pi(),1600,-0.5,1599.5);
    HistName = Form("h_mReCenterFull_%d",i_cent);
    h_mReCenterFull[i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),360,-1.0*TMath::Pi(),TMath::Pi(),1600,-0.5,1599.5);
  }
}

void StZdcSmdHistoManger::FillReCenterEP(TVector2 QEast, TVector2 QWest, TVector2 QFull, int Cent9, int runIndex)
{
  float PsiEast = TMath::ATan2(QEast.Y(),QEast.X()); h_mReCenterEast[Cent9]->Fill(PsiEast,runIndex);
  float PsiWest = TMath::ATan2(QWest.Y(),QWest.X()); h_mReCenterWest[Cent9]->Fill(PsiWest,runIndex);
  float PsiFull = TMath::ATan2(QFull.Y(),QFull.X()); h_mReCenterFull[Cent9]->Fill(PsiFull,runIndex);
}

void StZdcSmdHistoManger::WriteReCenterEP()
{
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    h_mReCenterEast[i_cent]->Write();
    h_mReCenterWest[i_cent]->Write();
    h_mReCenterFull[i_cent]->Write();
  }
}

//-------------------------------------------------------------------------------------------

void StZdcSmdHistoManger::InitShiftEP()
{
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    string HistName = Form("h_mShiftEast_%d",i_cent);
    h_mShiftEast[i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),360,-1.0*TMath::Pi(),TMath::Pi(),1600,-0.5,1599.5);
    HistName = Form("h_mShiftWest_%d",i_cent);
    h_mShiftWest[i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),360,-1.0*TMath::Pi(),TMath::Pi(),1600,-0.5,1599.5);
    HistName = Form("h_mShiftFull_%d",i_cent);
    h_mShiftFull[i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),360,-1.0*TMath::Pi(),TMath::Pi(),1600,-0.5,1599.5);
  }
}

void StZdcSmdHistoManger::FillShiftEP(TVector2 QEast, TVector2 QWest, TVector2 QFull, int Cent9, int runIndex)
{
  float PsiEast = TMath::ATan2(QEast.Y(),QEast.X()); h_mShiftEast[Cent9]->Fill(PsiEast,runIndex);
  float PsiWest = TMath::ATan2(QWest.Y(),QWest.X()); h_mShiftWest[Cent9]->Fill(PsiWest,runIndex);
  float PsiFull = TMath::ATan2(QFull.Y(),QFull.X()); h_mShiftFull[Cent9]->Fill(PsiFull,runIndex);
}

void StZdcSmdHistoManger::WriteShiftEP()
{
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    h_mShiftEast[i_cent]->Write();
    h_mShiftWest[i_cent]->Write();
    h_mShiftFull[i_cent]->Write();
  }
}

//-------------------------------------------------------------------------------------------
