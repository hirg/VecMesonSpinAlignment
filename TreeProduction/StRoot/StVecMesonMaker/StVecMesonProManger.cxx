#include "StVecMesonProManger.h"
#include "../../../Utility/StSpinAlignmentCons.h"
#include "TProfile2D.h"
#include "TProfile.h"

ClassImp(StVecMesonProManger)

TString StVecMesonProManger::mVStr[2] = {"pos","neg"};

//---------------------------------------------------------------------------------

StVecMesonProManger::StVecMesonProManger()
{
}

//---------------------------------------------------------------------------------

StVecMesonProManger::~StVecMesonProManger()
{
  /* */
}

//---------------------------------------------------------------------------------

void StVecMesonProManger::InitReCenter()
{
  for(Int_t i = 0; i < 2; i++) // vertex pos/neg
  {
    for(Int_t j = 0; j < 4; j++) // eta_gap
    {
      TString ProName;
      // Event Plane method
      ProName = Form("qx_2nd_Vertex_%s_EtaGap_%d_East_EP",mVStr[i].Data(),j);
      p_mq2x_East_EP[i][j] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5);// neg eta || x axis is RunIndex, y axis is Centrality
      ProName = Form("qy_2nd_Vertex_%s_EtaGap_%d_East_EP",mVStr[i].Data(),j);
      p_mq2y_East_EP[i][j] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5);// neg eta || x axis is RunIndex, y axis is Centrality
      ProName = Form("qx_2nd_Vertex_%s_EtaGap_%d_West_EP",mVStr[i].Data(),j);
      p_mq2x_West_EP[i][j] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5);// pos eta || x axis is RunIndex, y axis is Centrality
      ProName = Form("qy_2nd_Vertex_%s_EtaGap_%d_West_EP",mVStr[i].Data(),j);
      p_mq2y_West_EP[i][j] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5);// pos eta || x axis is RunIndex, y axis is Centrality
    }

    TString ProName_Full;
    // Event Plane method
    ProName_Full = Form("qx_2nd_Vertex_%s_Full_EP",mVStr[i].Data());
    p_mq2x_Full_EP[i] = new TProfile2D(ProName_Full.Data(),ProName_Full.Data(),1600,-0.5,1599.5,9,-0.5,8.5);// x axis is RunIndex, y axis is Centrality
    ProName_Full = Form("qy_2nd_Vertex_%s_Full_EP",mVStr[i].Data());
    p_mq2y_Full_EP[i] = new TProfile2D(ProName_Full.Data(),ProName_Full.Data(),1600,-0.5,1599.5,9,-0.5,8.5);// x axis is RunIndex, y axis is Centrality
    ProName_Full = Form("qx_3rd_Vertex_%s_Full_EP",mVStr[i].Data());
  }
}

//----------------------------------------------------------------------------

void StVecMesonProManger::InitShift()
{
  for(Int_t i = 0; i < 2; i++) // vertex pos/neg
  {
    for(Int_t j = 0; j < 4; j++) // eta_gap
    {
      for(Int_t k = 0; k < 5; k++) // Shift Order
      {
        TString ProName;
	// Event Plane method
	ProName = Form("CosPsi2_Vertex_%s_EtaGap_%d_Order_%d_East_EP",mVStr[i].Data(),j,k);
	p_mcos2_East_EP[i][j][k] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality
	ProName = Form("SinPsi2_Vertex_%s_EtaGap_%d_Order_%d_East_EP",mVStr[i].Data(),j,k);
	p_msin2_East_EP[i][j][k] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality
	ProName = Form("CosPsi2_Vertex_%s_EtaGap_%d_Order_%d_West_EP",mVStr[i].Data(),j,k);
	p_mcos2_West_EP[i][j][k] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality
	ProName = Form("SinPsi2_Vertex_%s_EtaGap_%d_Order_%d_West_EP",mVStr[i].Data(),j,k);
	p_msin2_West_EP[i][j][k] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality
      }
    }

    for(Int_t k = 0; k < 5; k++) // Shift Order
    {
      TString ProName_Full;
      // Event Plane method
      ProName_Full = Form("CosPsi2_Vertex_%s_Order_%d_Full_EP",mVStr[i].Data(),k);
      p_mcos2_Full_EP[i][k] = new TProfile2D(ProName_Full.Data(),ProName_Full.Data(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality
      ProName_Full = Form("SinPsi2_Vertex_%s_Order_%d_Full_EP",mVStr[i].Data(),k);
      p_msin2_Full_EP[i][k] = new TProfile2D(ProName_Full.Data(),ProName_Full.Data(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality
    }
  }
}

//----------------------------------------------------------------------------

void StVecMesonProManger::InitChargedFlow()
{
  // eta sub
  for(Int_t i = 0; i < 9; i++) // Centrality
  {
    for(Int_t j = 0; j < 4; j++) // eta_gap
    {
      TString ProName;
      ProName = Form("v2_Centrality_%d_EtaGap_%d_EP",i,j);
      p_mV2_EP[i][j] = new TProfile(ProName.Data(),ProName.Data(),30,0.15,5.0);
    }
  }
  for(Int_t j = 0; j < 4; j++) // eta_gap
  {
    TString ProName;
    ProName = Form("v2_minBias_EtaGap_%d_EP",j);
    p_mV2_EP[9][j] = new TProfile(ProName.Data(),ProName.Data(),30,0.15,5.0);
  }
  // random sub
  for(Int_t i = 0; i < 9; i++)
  {
    TString ProName_Ran;
    ProName_Ran = Form("v2_Centrality_%d_Ran_EP",i);
    p_mV2_Ran_EP[i] = new TProfile(ProName_Ran.Data(),ProName_Ran.Data(),30,0.15,5.0);
  }
  TString ProName_Ran_minBias;
  ProName_Ran_minBias = "v2_Centrality_minBias_Ran_EP";
  p_mV2_Ran_EP[9] = new TProfile(ProName_Ran_minBias.Data(),ProName_Ran_minBias.Data(),30,0.15,5.0);
}

//----------------------------------------------------------------------------

void StVecMesonProManger::FillTrackEast(TVector2 q2Vector, Int_t Cent9, Int_t RunIndex, Int_t i, Int_t j, Float_t pt) // i = vertex pos/neg, j = eta_gap
{
  const Float_t q2x = q2Vector.X();
  const Float_t q2y = q2Vector.Y();

  Float_t w;
  if(pt <= vmsa::mPrimPtWeight)
  {
    w = pt;
  }
  if(pt > vmsa::mPrimPtWeight)
  {
    w = vmsa::mPrimPtWeight;
  }

  // Event Plane method
  p_mq2x_East_EP[i][j]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)q2x,(Double_t)w);
  p_mq2y_East_EP[i][j]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)q2y,(Double_t)w);
}

void StVecMesonProManger::FillTrackWest(TVector2 q2Vector, Int_t Cent9, Int_t RunIndex, Int_t i, Int_t j, Float_t pt) // i = vertex pos/neg, j = eta_gap
{
  const Float_t q2x = q2Vector.X();
  const Float_t q2y = q2Vector.Y();

  Float_t w;
  if(pt <= vmsa::mPrimPtWeight)
  {
    w = pt;
  }
  if(pt > vmsa::mPrimPtWeight)
  {
    w = vmsa::mPrimPtWeight;
  }

  // Event Plane method
  p_mq2x_West_EP[i][j]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)q2x,(Double_t)w);
  p_mq2y_West_EP[i][j]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)q2y,(Double_t)w);
}

void StVecMesonProManger::FillTrackFull(TVector2 q2Vector, Int_t Cent9, Int_t RunIndex, Int_t i, Float_t pt) // i = vertex pos/neg
{
  const Float_t q2x = q2Vector.X();
  const Float_t q2y = q2Vector.Y();

  Float_t w;
  if(pt <= vmsa::mPrimPtWeight)
  {
    w = pt;
  }
  if(pt > vmsa::mPrimPtWeight)
  {
    w = vmsa::mPrimPtWeight;
  }

  // Event Plane method
  p_mq2x_Full_EP[i]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)q2x,(Double_t)w);
  p_mq2y_Full_EP[i]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)q2y,(Double_t)w);
}

//----------------------------------------------------------------------------
// Event Plane method
void StVecMesonProManger::FillEventEast_EP(TVector2 Psi2Vector, Int_t Cent9, Int_t RunIndex, Int_t i, Int_t j, Int_t k) // i = vertex pos/neg, j = eta_gap, k = ShiftOrder
{
  const Float_t cos2 = Psi2Vector.X();
  const Float_t sin2 = Psi2Vector.Y();

  p_mcos2_East_EP[i][j][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)cos2);
  p_msin2_East_EP[i][j][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)sin2);
}

void StVecMesonProManger::FillEventWest_EP(TVector2 Psi2Vector, Int_t Cent9, Int_t RunIndex, Int_t i, Int_t j, Int_t k) // i = vertex pos/neg, j = eta_gap, k = ShiftOrder
{
  const Float_t cos2 = Psi2Vector.X();
  const Float_t sin2 = Psi2Vector.Y();

  p_mcos2_West_EP[i][j][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)cos2);
  p_msin2_West_EP[i][j][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)sin2);
}

void StVecMesonProManger::FillEventFull_EP(TVector2 Psi2Vector, Int_t Cent9, Int_t RunIndex, Int_t i, Int_t k) // i = vertex pos/neg, k = ShiftOrder
{
  const Float_t cos2 = Psi2Vector.X();
  const Float_t sin2 = Psi2Vector.Y();

  p_mcos2_Full_EP[i][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)cos2);
  p_msin2_Full_EP[i][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)sin2);
}

//----------------------------------------------------------------------------
// eta sub
void StVecMesonProManger::FillEtaCharged2Flow(Float_t pt, Float_t v2_EP, Float_t Res2_EP, Float_t v2_SP, Float_t Res2_SP, Int_t Cent9, Int_t eta_gap, Double_t reweight)
{
  if(Res2_EP > 0.0)
  {
    p_mV2_EP[Cent9][eta_gap]->Fill(pt, v2_EP, reweight);
    p_mV2_EP[9][eta_gap]->Fill(pt, v2_EP, reweight);
  }
}

// random sub
void StVecMesonProManger::FillRanCharged2Flow(Float_t pt, Float_t v2_Ran_EP, Float_t Res2_Ran_EP, Float_t v2_Ran_SP, Float_t Res2_Ran_SP, Int_t Cent9, Double_t reweight)
{
  if(Res2_Ran_EP > 0.0)
  {
    p_mV2_Ran_EP[Cent9]->Fill(pt, v2_Ran_EP, reweight);
    p_mV2_Ran_EP[9]->Fill(pt, v2_Ran_EP, reweight);
  }
}

//----------------------------------------------------------------------------

void StVecMesonProManger::WriteReCenter()
{
  for(Int_t i = 0; i < 2; i++) // vertex pos/neg
  {
    for(Int_t j = 0; j < 4; j++) // eta_gap
    {
      // Event Plane method
      p_mq2x_East_EP[i][j]->Write();
      p_mq2y_East_EP[i][j]->Write();
      p_mq2x_West_EP[i][j]->Write();
      p_mq2y_West_EP[i][j]->Write();
    }
    // Event Plane method
    p_mq2x_Full_EP[i]->Write();
    p_mq2y_Full_EP[i]->Write();
  }
}

//----------------------------------------------------------------------------

void StVecMesonProManger::WriteShift()
{
  for(Int_t i = 0; i < 2; i++) // vertex pos/neg
  {
    for(Int_t j = 0; j < 4; j++) // eta_gap
    {
      for(Int_t k = 0; k < 5; k++) // Shift Order
      {
        // Event Plane method
	p_mcos2_East_EP[i][j][k]->Write();
	p_msin2_East_EP[i][j][k]->Write();
	p_mcos2_West_EP[i][j][k]->Write();
	p_msin2_West_EP[i][j][k]->Write();
      }
    }
    for(Int_t k = 0; k < 5; k++) // Shift Order
    {
      // Event Plane method
      p_mcos2_Full_EP[i][k]->Write();
      p_msin2_Full_EP[i][k]->Write();
    }
  }
}

//----------------------------------------------------------------------------

void StVecMesonProManger::WriteChargedFlow()
{
  for(Int_t i = 0; i < 10; i++)
  {
    for(Int_t j = 0; j < 4; j++)
    {
      p_mV2_EP[i][j]->Write();
    }
    p_mV2_Ran_EP[i]->Write();
  }
}

//----------------------------------------------------------------------------
