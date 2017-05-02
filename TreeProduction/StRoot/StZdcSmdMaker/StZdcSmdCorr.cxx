#include "StRoot/StZdcSmdMaker/StZdcSmdCorr.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "StRoot/StPicoDstMaker/StPicoDst.h"
#include "StRoot/StPicoDstMaker/StPicoEvent.h"
#include "StRoot/StPicoDstMaker/StPicoTrack.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"
#include "StMessMgr.h"
#include "TMath.h"
#include "TString.h"
#include "TFile.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TRandom3.h"
#include "TF1.h"

Double_t Resolution_Full(Double_t *x_val, Double_t *par)
{
  Double_t y;
  Double_t chi = x_val[0];
  Double_t arg = chi*chi/4.0;
  Double_t norm = TMath::Sqrt(TMath::Pi()/2.0)/2.0;

  y = norm*chi*TMath::Exp(-1.0*arg)*(TMath::BesselI0(arg)+TMath::BesselI1(arg));

  return y;
}

ClassImp(StZdcSmdCorrection)

//---------------------------------------------------------------------------------

StZdcSmdCorrection::StZdcSmdCorrection(Int_t energy)
{
  mEnergy = energy;
}

StZdcSmdCorrection::~StZdcSmdCorrection()
{
  /* */
}

void StZdcSmdCorrection::clear()
{
  for(int i_eastwest = 0; i_eastwest < 2; ++i_eastwest)
  {
    for(int i_verthori = 0; i_verthori < 2; ++i_verthori)
    {
      for(int i_slat = 0; i_slat < 8; ++i_slat)
      {
	mZdcSmd[i_eastwest][i_verthori][i_slat] = 0.0;
      }
    }
  }
}

//---------------------------------------------------------------------------------

void StZdcSmdCorrection::SetZdcSmd(int i_eastwest,int i_verthori,int i_slat,const Float_t zdcsmd) 
{
  mZdcSmd[i_eastwest][i_verthori][i_slat] = (zdcsmd > 0.) ? zdcsmd : 0.;
}

float StZdcSmdCorrection::GetZdcSmd(int i_eastwest,int i_verthori,int i_slat)
{
  return mZdcSmd[i_eastwest][i_verthori][i_slat];
}

//---------------------------------------------------------------------------------

