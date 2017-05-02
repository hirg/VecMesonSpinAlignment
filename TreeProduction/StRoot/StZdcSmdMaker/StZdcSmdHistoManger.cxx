#include "StRoot/StZdcSmdMaker/StZdcSmdHistoManger.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TMath.h"
#include "TString.h"

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
	cout << HistName.c_str() << endl;
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
