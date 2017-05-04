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

StZdcSmdCorrection::StZdcSmdCorrection(int energy)
{
  mEnergy = energy;

  for(int i_eastwest = 0; i_eastwest < 2; ++i_eastwest)
  {
    for(int i_verthori = 0; i_verthori < 2; ++i_verthori)
    {
      for(int i_slat = 0; i_slat < 8; ++i_slat)
      {
	mZdcSmd[i_eastwest][i_verthori][i_slat] = 0.0;
	mGainCorrFactor[i_eastwest][i_verthori][i_slat] = 0.0;
      }
    }
  }
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
	mGainCorrFactor[i_eastwest][i_verthori][i_slat] = 0.0;
      }
    }
  }
}

//---------------------------------------------------------------------------------

void StZdcSmdCorrection::SetZdcSmd(int eastwest, int verthori, int slat, const float zdcsmd) 
{
  mZdcSmd[eastwest][verthori][slat] = (zdcsmd > 0.) ? zdcsmd : 0.;
}

float StZdcSmdCorrection::GetZdcSmd(int eastwest, int verthori, int slat)
{
  return mZdcSmd[eastwest][verthori][slat];
}

//---------------------------------------------------------------------------------

void StZdcSmdCorrection::InitGainCorr()
{
  string InPutFile = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/ZDCSMD/GainCorrPar/merged_file/file_%s_GainCorrFac.root",vmsa::mBeamEnergy[mEnergy].c_str(),vmsa::mBeamEnergy[mEnergy].c_str());
  mFile_GainCorrPar = TFile::Open(InPutFile.c_str());
  for(int i_eastwest = 0; i_eastwest < 2; ++i_eastwest)
  {
    for(int i_verthori = 0; i_verthori < 2; ++i_verthori)
    {
      for(int i_slat = 0; i_slat < 8; ++i_slat)
      {
	string HistName = Form("h_mGainCorrFactor%s%s_%d",vmsa::mEastWest[i_eastwest].c_str(),vmsa::mVertHori[i_verthori].c_str(),i_slat);
	TH1F *h_GainCorrFac = (TH1F*)mFile_GainCorrPar->Get(HistName.c_str());
	mGainCorrFactor[i_eastwest][i_verthori][i_slat] = h_GainCorrFac->GetBinContent(1);
	cout << "i_eastwest = " << i_eastwest << ", i_verthori = " << i_verthori << ", i_slat = " << i_slat << ", mGainCorrFactor = " << mGainCorrFactor[i_eastwest][i_verthori][i_slat] << endl;
      }
    }
  }
}

void StZdcSmdCorrection::SetZdcSmdGainCorr(int eastwest, int verthori, int slat, const float zdcsmd)
{
  mZdcSmd[eastwest][verthori][slat] = (zdcsmd > 0.) ? zdcsmd/mGainCorrFactor[eastwest][verthori][slat] : 0.;
}

float StZdcSmdCorrection::GetZdcSmdGainCorr(int eastwest, int verthori, int slat)
{
  return mZdcSmd[eastwest][verthori][slat];
}

float StZdcSmdCorrection::GetPosition(int eastwest, int verthori, int slat, int mode)
{
  //get position of each slat

  float zdcsmd_vert[7] = {0.5,2,3.5,5,6.5,8,9.5};
  float zdcsmd_hori[8] = {1.25,3.25,5.25,7.25,9.25,11.25,13.25,15.25};

  double mZDCSMDCenterEastX = 0.0;
  double mZDCSMDCenterEastY = 0.0;
  double mZDCSMDCenterWestX = 0.0;
  double mZDCSMDCenterWestY = 0.0;


  if(mode > 1) // with beam center corrected
  {
    if(eastwest == 0 && verthori == 0) return zdcsmd_vert[slat]-mZDCSMDCenterEastX;
    if(eastwest == 1 && verthori == 0) return mZDCSMDCenterWestX-zdcsmd_vert[slat];
    if(eastwest == 0 && verthori == 1) return zdcsmd_hori[slat]/sqrt(2.)-mZDCSMDCenterEastY;
    if(eastwest == 1 && verthori == 1) return zdcsmd_hori[slat]/sqrt(2.)-mZDCSMDCenterWestY;
  }
  else // raw beam center returned
  {
    if(eastwest == 0 && verthori == 0) return zdcsmd_vert[slat];
    if(eastwest == 1 && verthori == 0) return -zdcsmd_vert[slat];
    if(eastwest == 0 && verthori == 1) return zdcsmd_hori[slat]/sqrt(2.);
    if(eastwest == 1 && verthori == 1) return zdcsmd_hori[slat]/sqrt(2.);
  }

  return 0;
}

TVector2 StZdcSmdCorrection::GetQEast(int mode) 
{
  TVector2 qVector(0.0,0.0);
  float qXsum = 0.; float qYsum = 0.;
  float qXwgt = 0.; float qYwgt = 0.;

  for(int i_vert = 0; i_vert < 7; ++i_vert) // vertical
  {
    qXsum += GetPosition(0,0,i_vert,mode)*GetZdcSmdGainCorr(0,0,i_vert);
    qXwgt += GetZdcSmdGainCorr(0,0,i_vert);
  }
  for(int i_hori = 0; i_hori < 8; ++i_hori) // horizontal
  {
    qYsum += GetPosition(0,1,i_hori,mode)*GetZdcSmdGainCorr(0,1,i_hori);
    qYwgt += GetZdcSmdGainCorr(0,1,i_hori);
  }

  if(qXwgt > 0.0 && qYwgt > 0.0)
    qVector.Set(qXsum/qXwgt,qYsum/qYwgt);

  return qVector;
}

TVector2 StZdcSmdCorrection::GetQWest(int mode)
{
  TVector2 qVector(0.0,0.0);
  float qXsum = 0.; float qYsum = 0.;
  float qXwgt = 0.; float qYwgt = 0.;

  for(int i_vert = 0; i_vert < 7; ++i_vert) // vertical
  {
    qXsum += GetPosition(1,0,i_vert,mode)*GetZdcSmdGainCorr(1,0,i_vert);
    qXwgt += GetZdcSmdGainCorr(1,0,i_vert);
  }
  for(int i_hori = 0; i_hori < 8; ++i_hori) // horizontal
  {
    qYsum += GetPosition(1,1,i_hori,mode)*GetZdcSmdGainCorr(1,1,i_hori);
    qYwgt += GetZdcSmdGainCorr(1,1,i_hori);
  }

  if(qXwgt > 0.0 && qYwgt > 0.0)
    qVector.Set(qXsum/qXwgt,qYsum/qYwgt);

  return qVector;
}
//---------------------------------------------------------------------------------
