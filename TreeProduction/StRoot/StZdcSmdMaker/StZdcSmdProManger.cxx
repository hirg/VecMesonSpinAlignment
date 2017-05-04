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

    ProName = Form("p_mQEastVertical_%s",mVStr[i].c_str());
    p_mQEastVertical[i] = new TProfile2D(ProName.c_str(),ProName.c_str(),1600,-0.5,1599.5,9,-0.5,8.5);// neg eta || x axis is RunIndex, y axis is Centrality
    ProName = Form("p_mQEastHorizontal_%s",mVStr[i].c_str());
    p_mQEastHorizontal[i] = new TProfile2D(ProName.c_str(),ProName.c_str(),1600,-0.5,1599.5,9,-0.5,8.5);

    ProName = Form("p_mQWestVertical_%s",mVStr[i].c_str());
    p_mQWestVertical[i] = new TProfile2D(ProName.c_str(),ProName.c_str(),1600,-0.5,1599.5,9,-0.5,8.5);
    ProName = Form("p_mQWestHorizontal_%s",mVStr[i].c_str());
    p_mQWestHorizontal[i] = new TProfile2D(ProName.c_str(),ProName.c_str(),1600,-0.5,1599.5,9,-0.5,8.5);
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
