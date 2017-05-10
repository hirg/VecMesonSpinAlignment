#include "StRoot/StZdcSmdMaker/StZdcSmdProManger.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TMath.h"

ClassImp(StZdcSmdProManger)

string StZdcSmdProManger::mVStr[2] = {"pos","neg"};

//---------------------------------------------------------------------------------

StZdcSmdProManger::StZdcSmdProManger()
{
}

//---------------------------------------------------------------------------------

StZdcSmdProManger::~StZdcSmdProManger()
{
  /* */
}

//---------------------------------------------------------------------------------

void StZdcSmdProManger::InitReCenter()
{
  for(Int_t i_vz = 0; i_vz < 2; ++i_vz) // vertex pos/neg
  {
    string ProName;

    ProName = Form("p_mQEastVertical_%s",mVStr[i_vz].c_str());
    p_mQEastVertical[i_vz] = new TProfile2D(ProName.c_str(),ProName.c_str(),1600,-0.5,1599.5,9,-0.5,8.5);// neg eta || x axis is RunIndex, y axis is Centrality
    ProName = Form("p_mQEastHorizontal_%s",mVStr[i_vz].c_str());
    p_mQEastHorizontal[i_vz] = new TProfile2D(ProName.c_str(),ProName.c_str(),1600,-0.5,1599.5,9,-0.5,8.5);

    ProName = Form("p_mQWestVertical_%s",mVStr[i_vz].c_str());
    p_mQWestVertical[i_vz] = new TProfile2D(ProName.c_str(),ProName.c_str(),1600,-0.5,1599.5,9,-0.5,8.5);
    ProName = Form("p_mQWestHorizontal_%s",mVStr[i_vz].c_str());
    p_mQWestHorizontal[i_vz] = new TProfile2D(ProName.c_str(),ProName.c_str(),1600,-0.5,1599.5,9,-0.5,8.5);
  }
}

void StZdcSmdProManger::FillReCenterEast(TVector2 qVector, int Cent9, int RunIndex, int vz_sign) // vz_sign = vertex pos/neg 
{
  const double qx = qVector.X();
  const double qy = qVector.Y();

  // Event Plane method
  p_mQEastVertical[vz_sign]->Fill((double)RunIndex,(double)Cent9,(double)qx);
  p_mQEastHorizontal[vz_sign]->Fill((double)RunIndex,(double)Cent9,(double)qy);
}

void StZdcSmdProManger::FillReCenterWest(TVector2 qVector, int Cent9, int RunIndex, int vz_sign) // vz_sign = vertex pos/neg 
{
  const double qx = qVector.X();
  const double qy = qVector.Y();

  p_mQWestVertical[vz_sign]->Fill((double)RunIndex,(double)Cent9,(double)qx);
  p_mQWestHorizontal[vz_sign]->Fill((double)RunIndex,(double)Cent9,(double)qy);
}

void StZdcSmdProManger::WriteReCenter()
{
  for(Int_t i_vz = 0; i_vz < 2; ++i_vz) // vertex pos/neg
  {
    p_mQEastVertical[i_vz]->Write();
    p_mQEastHorizontal[i_vz]->Write();
    p_mQWestVertical[i_vz]->Write();
    p_mQWestHorizontal[i_vz]->Write();
  }
}

//---------------------------------------------------------------------------------

void StZdcSmdProManger::InitShift()
{
  for(int i_vz = 0; i_vz < 2; ++i_vz) // vertex pos/neg
  {
    for(int i_shift = 0; i_shift < 20; ++i_shift) // Shift Order
    {
      string ProName;

      ProName = Form("p_mQEastCos_%s_%d",mVStr[i_vz].c_str(),i_shift);
      p_mQEastCos[i_vz][i_shift] = new TProfile2D(ProName.c_str(),ProName.c_str(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality
      ProName = Form("p_mQEastSin_%s_%d",mVStr[i_vz].c_str(),i_shift);
      p_mQEastSin[i_vz][i_shift] = new TProfile2D(ProName.c_str(),ProName.c_str(),1600,-0.5,1599.5,9,-0.5,8.5);

      ProName = Form("p_mQWestCos_%s_%d",mVStr[i_vz].c_str(),i_shift);
      p_mQWestCos[i_vz][i_shift] = new TProfile2D(ProName.c_str(),ProName.c_str(),1600,-0.5,1599.5,9,-0.5,8.5);
      ProName = Form("p_mQWestSin_%s_%d",mVStr[i_vz].c_str(),i_shift);
      p_mQWestSin[i_vz][i_shift] = new TProfile2D(ProName.c_str(),ProName.c_str(),1600,-0.5,1599.5,9,-0.5,8.5);
    }
  }
}

void StZdcSmdProManger::FillShiftEast(TVector2 qVector, int Cent9, int RunIndex, int vz_sign)
{
  float Psi = TMath::ATan2(qVector.Y(),qVector.X());
  for(int i_shift = 0; i_shift < 20; ++i_shift)
  {
    p_mQEastCos[vz_sign][i_shift]->Fill((double)RunIndex,(double)Cent9,TMath::Cos((i_shift+1)*Psi));
    p_mQEastSin[vz_sign][i_shift]->Fill((double)RunIndex,(double)Cent9,TMath::Sin((i_shift+1)*Psi));
  }
}

void StZdcSmdProManger::FillShiftWest(TVector2 qVector, int Cent9, int RunIndex, int vz_sign)
{
  float Psi = TMath::ATan2(qVector.Y(),qVector.X());
  for(int i_shift = 0; i_shift < 20; ++i_shift)
  {
    p_mQWestCos[vz_sign][i_shift]->Fill((double)RunIndex,(double)Cent9,TMath::Cos((i_shift+1)*Psi));
    p_mQWestSin[vz_sign][i_shift]->Fill((double)RunIndex,(double)Cent9,TMath::Sin((i_shift+1)*Psi));
  }
}

void StZdcSmdProManger::WriteShift()
{
  for(int i_vz = 0; i_vz < 2; ++i_vz) // vertex pos/neg
  {
    for(int i_shift = 0; i_shift < 20; ++i_shift) // Shift Order
    {
      p_mQEastCos[i_vz][i_shift]->Write();
      p_mQEastSin[i_vz][i_shift]->Write();
      p_mQWestCos[i_vz][i_shift]->Write();
      p_mQWestSin[i_vz][i_shift]->Write();
    }
  }
}

//---------------------------------------------------------------------------------

void StZdcSmdProManger::InitShiftFull()
{
  for(int i_vz = 0; i_vz < 2; ++i_vz) // vertex pos/neg
  {
    for(int i_shift = 0; i_shift < 20; ++i_shift) // Shift Order
    {
      string ProName;

      ProName = Form("p_mQFullCos_%s_%d",mVStr[i_vz].c_str(),i_shift);
      p_mQFullCos[i_vz][i_shift] = new TProfile2D(ProName.c_str(),ProName.c_str(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality
      ProName = Form("p_mQFullSin_%s_%d",mVStr[i_vz].c_str(),i_shift);
      p_mQFullSin[i_vz][i_shift] = new TProfile2D(ProName.c_str(),ProName.c_str(),1600,-0.5,1599.5,9,-0.5,8.5);
    }
  }
}

void StZdcSmdProManger::FillShiftFull(TVector2 qVector, int Cent9, int RunIndex, int vz_sign)
{
  float Psi = TMath::ATan2(qVector.Y(),qVector.X());
  for(int i_shift = 0; i_shift < 20; ++i_shift)
  {
    p_mQFullCos[vz_sign][i_shift]->Fill((double)RunIndex,(double)Cent9,TMath::Cos((i_shift+1)*Psi));
    p_mQFullSin[vz_sign][i_shift]->Fill((double)RunIndex,(double)Cent9,TMath::Sin((i_shift+1)*Psi));
  }
}

void StZdcSmdProManger::WriteShiftFull()
{
  for(int i_vz = 0; i_vz < 2; ++i_vz) // vertex pos/neg
  {
    for(int i_shift = 0; i_shift < 20; ++i_shift) // Shift Order
    {
      p_mQFullCos[i_vz][i_shift]->Write();
      p_mQFullSin[i_vz][i_shift]->Write();
    }
  }
}
//---------------------------------------------------------------------------------

void StZdcSmdProManger::InitResolution()
{
  p_mResolution = new TProfile("p_mResolution","p_mResolution",9,-0.5,8.5);
}

void StZdcSmdProManger::FillResolution(TVector2 QEast, TVector2 QWest, int Cent9)
{
  float Psi_East = TMath::ATan2(QEast.Y(),QEast.X());
  float Psi_West = TMath::ATan2(QWest.Y(),QWest.X());
  float resolution = TMath::Cos(Psi_West-Psi_East+TMath::Pi());
  p_mResolution->Fill((double)Cent9,resolution);
}

void StZdcSmdProManger::WriteResolution()
{
  p_mResolution->Write();
}

//---------------------------------------------------------------------------------

void StZdcSmdProManger::InitDirectedFlow()
{
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    string ProName = Form("p_mDirectedFlow_%d",i_cent);
    p_mDirectedFlow[i_cent] = new TProfile(ProName.c_str(),ProName.c_str(),10,-1.0,1.0);
  }
  p_mDirectedFlowCom = new TProfile("p_mDirectedFlowCom","p_mDirectedFlowCom",10,-1.0,1.0);
}

void StZdcSmdProManger::FillDirectedFlow(int Cent9, float eta, float pt, float v1, float resolution, float reweight)
{
  if(pt > 0.15 && pt < 2.0)
  {
    p_mDirectedFlow[Cent9]->Fill(eta,v1/resolution);
    if(Cent9 >= 2 && Cent9 <= 4) p_mDirectedFlowCom->Fill(eta,v1/resolution,reweight); // 30-60%
  }
}

void StZdcSmdProManger::WriteDirectedFlow()
{
  for(int i_cent = 0; i_cent < 9; ++i_cent) p_mDirectedFlow[i_cent]->Write();
  p_mDirectedFlowCom->Write();
}
//---------------------------------------------------------------------------------
