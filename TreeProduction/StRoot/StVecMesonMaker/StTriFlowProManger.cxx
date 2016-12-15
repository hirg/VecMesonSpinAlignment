#include "StTriFlowProManger.h"
#include "StTriFlowConstants.h"
#include "TProfile2D.h"
#include "TProfile.h"

ClassImp(StTriFlowProManger)

TString StTriFlowProManger::mVStr[2] = {"pos","neg"};

//---------------------------------------------------------------------------------

StTriFlowProManger::StTriFlowProManger()
{
}

//---------------------------------------------------------------------------------

StTriFlowProManger::~StTriFlowProManger()
{
  /* */
}

//---------------------------------------------------------------------------------

void StTriFlowProManger::InitReCenter()
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

      ProName = Form("qx_3rd_Vertex_%s_EtaGap_%d_East_EP",mVStr[i].Data(),j);
      p_mq3x_East_EP[i][j] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5);// neg eta || x axis is RunIndex, y axis is Centrality
      ProName = Form("qy_3rd_Vertex_%s_EtaGap_%d_East_EP",mVStr[i].Data(),j);
      p_mq3y_East_EP[i][j] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5);// neg eta || x axis is RunIndex, y axis is Centrality
      ProName = Form("qx_3rd_Vertex_%s_EtaGap_%d_West_EP",mVStr[i].Data(),j);
      p_mq3x_West_EP[i][j] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5);// pos eta || x axis is RunIndex, y axis is Centrality
      ProName = Form("qy_3rd_Vertex_%s_EtaGap_%d_West_EP",mVStr[i].Data(),j);
      p_mq3y_West_EP[i][j] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5);// pos eta || x axis is RunIndex, y axis is Centrality

      // Scalar Product method
      ProName = Form("qx_2nd_Vertex_%s_EtaGap_%d_East_SP",mVStr[i].Data(),j);
      p_mq2x_East_SP[i][j] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5);// neg eta || x axis is RunIndex, y axis is Centrality
      ProName = Form("qy_2nd_Vertex_%s_EtaGap_%d_East_SP",mVStr[i].Data(),j);
      p_mq2y_East_SP[i][j] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5);// neg eta || x axis is RunIndex, y axis is Centrality
      ProName = Form("qx_2nd_Vertex_%s_EtaGap_%d_West_SP",mVStr[i].Data(),j);
      p_mq2x_West_SP[i][j] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5);// pos eta || x axis is RunIndex, y axis is Centrality
      ProName = Form("qy_2nd_Vertex_%s_EtaGap_%d_West_SP",mVStr[i].Data(),j);
      p_mq2y_West_SP[i][j] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5);// pos eta || x axis is RunIndex, y axis is Centrality

      ProName = Form("qx_3rd_Vertex_%s_EtaGap_%d_East_SP",mVStr[i].Data(),j);
      p_mq3x_East_SP[i][j] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5);// neg eta || x axis is RunIndex, y axis is Centrality
      ProName = Form("qy_3rd_Vertex_%s_EtaGap_%d_East_SP",mVStr[i].Data(),j);
      p_mq3y_East_SP[i][j] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5);// neg eta || x axis is RunIndex, y axis is Centrality
      ProName = Form("qx_3rd_Vertex_%s_EtaGap_%d_West_SP",mVStr[i].Data(),j);
      p_mq3x_West_SP[i][j] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5);// pos eta || x axis is RunIndex, y axis is Centrality
      ProName = Form("qy_3rd_Vertex_%s_EtaGap_%d_West_SP",mVStr[i].Data(),j);
      p_mq3y_West_SP[i][j] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5);// pos eta || x axis is RunIndex, y axis is Centrality
    }

    TString ProName_Full;
    // Event Plane method
    ProName_Full = Form("qx_2nd_Vertex_%s_Full_EP",mVStr[i].Data());
    p_mq2x_Full_EP[i] = new TProfile2D(ProName_Full.Data(),ProName_Full.Data(),1600,-0.5,1599.5,9,-0.5,8.5);// x axis is RunIndex, y axis is Centrality
    ProName_Full = Form("qy_2nd_Vertex_%s_Full_EP",mVStr[i].Data());
    p_mq2y_Full_EP[i] = new TProfile2D(ProName_Full.Data(),ProName_Full.Data(),1600,-0.5,1599.5,9,-0.5,8.5);// x axis is RunIndex, y axis is Centrality
    ProName_Full = Form("qx_3rd_Vertex_%s_Full_EP",mVStr[i].Data());
    p_mq3x_Full_EP[i] = new TProfile2D(ProName_Full.Data(),ProName_Full.Data(),1600,-0.5,1599.5,9,-0.5,8.5);// x axis is RunIndex, y axis is Centrality
    ProName_Full = Form("qy_3rd_Vertex_%s_Full_EP",mVStr[i].Data());
    p_mq3y_Full_EP[i] = new TProfile2D(ProName_Full.Data(),ProName_Full.Data(),1600,-0.5,1599.5,9,-0.5,8.5);// x axis is RunIndex, y axis is Centrality

    // Scalar Product method
    ProName_Full = Form("qx_2nd_Vertex_%s_Full_SP",mVStr[i].Data());
    p_mq2x_Full_SP[i] = new TProfile2D(ProName_Full.Data(),ProName_Full.Data(),1600,-0.5,1599.5,9,-0.5,8.5);// x axis is RunIndex, y axis is Centrality
    ProName_Full = Form("qy_2nd_Vertex_%s_Full_SP",mVStr[i].Data());
    p_mq2y_Full_SP[i] = new TProfile2D(ProName_Full.Data(),ProName_Full.Data(),1600,-0.5,1599.5,9,-0.5,8.5);// x axis is RunIndex, y axis is Centrality
    ProName_Full = Form("qx_3rd_Vertex_%s_Full_SP",mVStr[i].Data());
    p_mq3x_Full_SP[i] = new TProfile2D(ProName_Full.Data(),ProName_Full.Data(),1600,-0.5,1599.5,9,-0.5,8.5);// x axis is RunIndex, y axis is Centrality
    ProName_Full = Form("qy_3rd_Vertex_%s_Full_SP",mVStr[i].Data());
    p_mq3y_Full_SP[i] = new TProfile2D(ProName_Full.Data(),ProName_Full.Data(),1600,-0.5,1599.5,9,-0.5,8.5);// x axis is RunIndex, y axis is Centrality
  }
}

//----------------------------------------------------------------------------

void StTriFlowProManger::InitShift()
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

	ProName = Form("CosPsi3_Vertex_%s_EtaGap_%d_Order_%d_East_EP",mVStr[i].Data(),j,k);
	p_mcos3_East_EP[i][j][k] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality
	ProName = Form("SinPsi3_Vertex_%s_EtaGap_%d_Order_%d_East_EP",mVStr[i].Data(),j,k);
	p_msin3_East_EP[i][j][k] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality
	ProName = Form("CosPsi3_Vertex_%s_EtaGap_%d_Order_%d_West_EP",mVStr[i].Data(),j,k);
	p_mcos3_West_EP[i][j][k] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality
	ProName = Form("SinPsi3_Vertex_%s_EtaGap_%d_Order_%d_West_EP",mVStr[i].Data(),j,k);
	p_msin3_West_EP[i][j][k] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality

	// Scalor Product method
	ProName = Form("CosPsi2_Vertex_%s_EtaGap_%d_Order_%d_East_SP",mVStr[i].Data(),j,k);
	p_mcos2_East_SP[i][j][k] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality
	ProName = Form("SinPsi2_Vertex_%s_EtaGap_%d_Order_%d_East_SP",mVStr[i].Data(),j,k);
	p_msin2_East_SP[i][j][k] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality
	ProName = Form("CosPsi2_Vertex_%s_EtaGap_%d_Order_%d_West_SP",mVStr[i].Data(),j,k);
	p_mcos2_West_SP[i][j][k] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality
	ProName = Form("SinPsi2_Vertex_%s_EtaGap_%d_Order_%d_West_SP",mVStr[i].Data(),j,k);
	p_msin2_West_SP[i][j][k] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality

	ProName = Form("CosPsi3_Vertex_%s_EtaGap_%d_Order_%d_East_SP",mVStr[i].Data(),j,k);
	p_mcos3_East_SP[i][j][k] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality
	ProName = Form("SinPsi3_Vertex_%s_EtaGap_%d_Order_%d_East_SP",mVStr[i].Data(),j,k);
	p_msin3_East_SP[i][j][k] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality
	ProName = Form("CosPsi3_Vertex_%s_EtaGap_%d_Order_%d_West_SP",mVStr[i].Data(),j,k);
	p_mcos3_West_SP[i][j][k] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality
	ProName = Form("SinPsi3_Vertex_%s_EtaGap_%d_Order_%d_West_SP",mVStr[i].Data(),j,k);
	p_msin3_West_SP[i][j][k] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality
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
      
      ProName_Full = Form("CosPsi3_Vertex_%s_Order_%d_Full_EP",mVStr[i].Data(),k);
      p_mcos3_Full_EP[i][k] = new TProfile2D(ProName_Full.Data(),ProName_Full.Data(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality
      ProName_Full = Form("SinPsi3_Vertex_%s_Order_%d_Full_EP",mVStr[i].Data(),k);
      p_msin3_Full_EP[i][k] = new TProfile2D(ProName_Full.Data(),ProName_Full.Data(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality

      // Scalor Product method
      ProName_Full = Form("CosPsi2_Vertex_%s_Order_%d_Full_SP",mVStr[i].Data(),k);
      p_mcos2_Full_SP[i][k] = new TProfile2D(ProName_Full.Data(),ProName_Full.Data(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality
      ProName_Full = Form("SinPsi2_Vertex_%s_Order_%d_Full_SP",mVStr[i].Data(),k);
      p_msin2_Full_SP[i][k] = new TProfile2D(ProName_Full.Data(),ProName_Full.Data(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality
      
      ProName_Full = Form("CosPsi3_Vertex_%s_Order_%d_Full_SP",mVStr[i].Data(),k);
      p_mcos3_Full_SP[i][k] = new TProfile2D(ProName_Full.Data(),ProName_Full.Data(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality
      ProName_Full = Form("SinPsi3_Vertex_%s_Order_%d_Full_SP",mVStr[i].Data(),k);
      p_msin3_Full_SP[i][k] = new TProfile2D(ProName_Full.Data(),ProName_Full.Data(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality
    }
  }
}

//----------------------------------------------------------------------------

void StTriFlowProManger::InitChargedFlow()
{
  // eta sub
  for(Int_t i = 0; i < 9; i++) // Centrality
  {
    for(Int_t j = 0; j < 4; j++) // eta_gap
    {
      TString ProName;
      ProName = Form("v2_Centrality_%d_EtaGap_%d_EP",i,j);
      p_mV2_EP[i][j] = new TProfile(ProName.Data(),ProName.Data(),30,0.15,5.0);
      ProName = Form("v3_Centrality_%d_EtaGap_%d_EP",i,j);
      p_mV3_EP[i][j] = new TProfile(ProName.Data(),ProName.Data(),30,0.15,5.0);
      ProName = Form("v2_Centrality_%d_EtaGap_%d_SP",i,j);
      p_mV2_SP[i][j] = new TProfile(ProName.Data(),ProName.Data(),30,0.15,5.0);
      ProName = Form("v3_Centrality_%d_EtaGap_%d_SP",i,j);
      p_mV3_SP[i][j] = new TProfile(ProName.Data(),ProName.Data(),30,0.15,5.0);
    }
  }
  for(Int_t j = 0; j < 4; j++) // eta_gap
  {
    TString ProName;
    ProName = Form("v2_minBias_EtaGap_%d_EP",j);
    p_mV2_EP[9][j] = new TProfile(ProName.Data(),ProName.Data(),30,0.15,5.0);
    ProName = Form("v3_minBias_EtaGap_%d_EP",j);
    p_mV3_EP[9][j] = new TProfile(ProName.Data(),ProName.Data(),30,0.15,5.0);
    ProName = Form("v2_minBias_EtaGap_%d_SP",j);
    p_mV2_SP[9][j] = new TProfile(ProName.Data(),ProName.Data(),30,0.15,5.0);
    ProName = Form("v3_minBias_EtaGap_%d_SP",j);
    p_mV3_SP[9][j] = new TProfile(ProName.Data(),ProName.Data(),30,0.15,5.0);
  }
  // random sub
  for(Int_t i = 0; i < 9; i++)
  {
    TString ProName_Ran;
    ProName_Ran = Form("v2_Centrality_%d_Ran_EP",i);
    p_mV2_Ran_EP[i] = new TProfile(ProName_Ran.Data(),ProName_Ran.Data(),30,0.15,5.0);
    ProName_Ran = Form("v3_Centrality_%d_Ran_EP",i);
    p_mV3_Ran_EP[i] = new TProfile(ProName_Ran.Data(),ProName_Ran.Data(),30,0.15,5.0);
    ProName_Ran = Form("v2_Centrality_%d_Ran_SP",i);
    p_mV2_Ran_SP[i] = new TProfile(ProName_Ran.Data(),ProName_Ran.Data(),30,0.15,5.0);
    ProName_Ran = Form("v3_Centrality_%d_Ran_SP",i);
    p_mV3_Ran_SP[i] = new TProfile(ProName_Ran.Data(),ProName_Ran.Data(),30,0.15,5.0);
  }
  TString ProName_Ran_minBias;
  ProName_Ran_minBias = "v2_Centrality_minBias_Ran_EP";
  p_mV2_Ran_EP[9] = new TProfile(ProName_Ran_minBias.Data(),ProName_Ran_minBias.Data(),30,0.15,5.0);
  ProName_Ran_minBias = "v3_Centrality_minBias_Ran_EP";
  p_mV3_Ran_EP[9] = new TProfile(ProName_Ran_minBias.Data(),ProName_Ran_minBias.Data(),30,0.15,5.0);
  ProName_Ran_minBias = "v2_Centrality_minBias_Ran_SP";
  p_mV2_Ran_SP[9] = new TProfile(ProName_Ran_minBias.Data(),ProName_Ran_minBias.Data(),30,0.15,5.0);
  ProName_Ran_minBias = "v3_Centrality_minBias_Ran_SP";
  p_mV3_Ran_SP[9] = new TProfile(ProName_Ran_minBias.Data(),ProName_Ran_minBias.Data(),30,0.15,5.0);
}

//----------------------------------------------------------------------------

void StTriFlowProManger::FillTrackEast(TVector2 q2Vector, TVector2 q3Vector, Int_t Cent9, Int_t RunIndex, Int_t i, Int_t j, Float_t pt) // i = vertex pos/neg, j = eta_gap
{
  const Float_t q2x = q2Vector.X();
  const Float_t q2y = q2Vector.Y();
  const Float_t q3x = q3Vector.X();
  const Float_t q3y = q3Vector.Y();

  Float_t w;
  if(pt <= TriFlow::mPrimPtWeight)
  {
    w = pt;
  }
  if(pt > TriFlow::mPrimPtWeight)
  {
    w = TriFlow::mPrimPtWeight;
  }

  // Event Plane method
  p_mq2x_East_EP[i][j]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)q2x,(Double_t)w);
  p_mq2y_East_EP[i][j]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)q2y,(Double_t)w);

  p_mq3x_East_EP[i][j]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)q3x,(Double_t)w);
  p_mq3y_East_EP[i][j]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)q3y,(Double_t)w);
  
  // Scalar Product method
  p_mq2x_East_SP[i][j]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)q2x);
  p_mq2y_East_SP[i][j]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)q2y);

  p_mq3x_East_SP[i][j]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)q3x);
  p_mq3y_East_SP[i][j]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)q3y);
}

void StTriFlowProManger::FillTrackWest(TVector2 q2Vector, TVector2 q3Vector, Int_t Cent9, Int_t RunIndex, Int_t i, Int_t j, Float_t pt) // i = vertex pos/neg, j = eta_gap
{
  const Float_t q2x = q2Vector.X();
  const Float_t q2y = q2Vector.Y();
  const Float_t q3x = q3Vector.X();
  const Float_t q3y = q3Vector.Y();

  Float_t w;
  if(pt <= TriFlow::mPrimPtWeight)
  {
    w = pt;
  }
  if(pt > TriFlow::mPrimPtWeight)
  {
    w = TriFlow::mPrimPtWeight;
  }

  // Event Plane method
  p_mq2x_West_EP[i][j]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)q2x,(Double_t)w);
  p_mq2y_West_EP[i][j]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)q2y,(Double_t)w);

  p_mq3x_West_EP[i][j]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)q3x,(Double_t)w);
  p_mq3y_West_EP[i][j]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)q3y,(Double_t)w);

  // Scalar Product method
  p_mq2x_West_SP[i][j]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)q2x);
  p_mq2y_West_SP[i][j]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)q2y);

  p_mq3x_West_SP[i][j]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)q3x);
  p_mq3y_West_SP[i][j]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)q3y);
}

void StTriFlowProManger::FillTrackFull(TVector2 q2Vector, TVector2 q3Vector, Int_t Cent9, Int_t RunIndex, Int_t i, Float_t pt) // i = vertex pos/neg
{
  const Float_t q2x = q2Vector.X();
  const Float_t q2y = q2Vector.Y();
  const Float_t q3x = q3Vector.X();
  const Float_t q3y = q3Vector.Y();

  Float_t w;
  if(pt <= TriFlow::mPrimPtWeight)
  {
    w = pt;
  }
  if(pt > TriFlow::mPrimPtWeight)
  {
    w = TriFlow::mPrimPtWeight;
  }

  // Event Plane method
  p_mq2x_Full_EP[i]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)q2x,(Double_t)w);
  p_mq2y_Full_EP[i]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)q2y,(Double_t)w);

  p_mq3x_Full_EP[i]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)q3x,(Double_t)w);
  p_mq3y_Full_EP[i]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)q3y,(Double_t)w);

  // Scalar Product method
  p_mq2x_Full_SP[i]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)q2x);
  p_mq2y_Full_SP[i]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)q2y);

  p_mq3x_Full_SP[i]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)q3x);
  p_mq3y_Full_SP[i]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)q3y);
}

//----------------------------------------------------------------------------
// Event Plane method
void StTriFlowProManger::FillEventEast_EP(TVector2 Psi2Vector, TVector2 Psi3Vector, Int_t Cent9, Int_t RunIndex, Int_t i, Int_t j, Int_t k) // i = vertex pos/neg, j = eta_gap, k = ShiftOrder
{
  const Float_t cos2 = Psi2Vector.X();
  const Float_t sin2 = Psi2Vector.Y();
  const Float_t cos3 = Psi3Vector.X();
  const Float_t sin3 = Psi3Vector.Y();

  p_mcos2_East_EP[i][j][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)cos2);
  p_msin2_East_EP[i][j][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)sin2);

  p_mcos3_East_EP[i][j][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)cos3);
  p_msin3_East_EP[i][j][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)sin3);
}

void StTriFlowProManger::FillEventWest_EP(TVector2 Psi2Vector, TVector2 Psi3Vector, Int_t Cent9, Int_t RunIndex, Int_t i, Int_t j, Int_t k) // i = vertex pos/neg, j = eta_gap, k = ShiftOrder
{
  const Float_t cos2 = Psi2Vector.X();
  const Float_t sin2 = Psi2Vector.Y();
  const Float_t cos3 = Psi3Vector.X();
  const Float_t sin3 = Psi3Vector.Y();

  p_mcos2_West_EP[i][j][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)cos2);
  p_msin2_West_EP[i][j][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)sin2);

  p_mcos3_West_EP[i][j][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)cos3);
  p_msin3_West_EP[i][j][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)sin3);
}

void StTriFlowProManger::FillEventFull_EP(TVector2 Psi2Vector, TVector2 Psi3Vector, Int_t Cent9, Int_t RunIndex, Int_t i, Int_t k) // i = vertex pos/neg, k = ShiftOrder
{
  const Float_t cos2 = Psi2Vector.X();
  const Float_t sin2 = Psi2Vector.Y();
  const Float_t cos3 = Psi3Vector.X();
  const Float_t sin3 = Psi3Vector.Y();

  p_mcos2_Full_EP[i][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)cos2);
  p_msin2_Full_EP[i][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)sin2);

  p_mcos3_Full_EP[i][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)cos3);
  p_msin3_Full_EP[i][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)sin3);
}

//----------------------------------------------------------------------------
// Scalor Product method
void StTriFlowProManger::FillEventEast_SP(TVector2 Psi2Vector, TVector2 Psi3Vector, Int_t Cent9, Int_t RunIndex, Int_t i, Int_t j, Int_t k) // i = vertex pos/neg, j = eta_gap, k = ShiftOrder
{
  const Float_t cos2 = Psi2Vector.X();
  const Float_t sin2 = Psi2Vector.Y();
  const Float_t cos3 = Psi3Vector.X();
  const Float_t sin3 = Psi3Vector.Y();

  p_mcos2_East_SP[i][j][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)cos2);
  p_msin2_East_SP[i][j][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)sin2);

  p_mcos3_East_SP[i][j][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)cos3);
  p_msin3_East_SP[i][j][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)sin3);
}

void StTriFlowProManger::FillEventWest_SP(TVector2 Psi2Vector, TVector2 Psi3Vector, Int_t Cent9, Int_t RunIndex, Int_t i, Int_t j, Int_t k) // i = vertex pos/neg, j = eta_gap, k = ShiftOrder
{
  const Float_t cos2 = Psi2Vector.X();
  const Float_t sin2 = Psi2Vector.Y();
  const Float_t cos3 = Psi3Vector.X();
  const Float_t sin3 = Psi3Vector.Y();

  p_mcos2_West_SP[i][j][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)cos2);
  p_msin2_West_SP[i][j][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)sin2);

  p_mcos3_West_SP[i][j][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)cos3);
  p_msin3_West_SP[i][j][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)sin3);
}

void StTriFlowProManger::FillEventFull_SP(TVector2 Psi2Vector, TVector2 Psi3Vector, Int_t Cent9, Int_t RunIndex, Int_t i, Int_t k) // i = vertex pos/neg, k = ShiftOrder
{
  const Float_t cos2 = Psi2Vector.X();
  const Float_t sin2 = Psi2Vector.Y();
  const Float_t cos3 = Psi3Vector.X();
  const Float_t sin3 = Psi3Vector.Y();

  p_mcos2_Full_SP[i][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)cos2);
  p_msin2_Full_SP[i][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)sin2);

  p_mcos3_Full_SP[i][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)cos3);
  p_msin3_Full_SP[i][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)sin3);
}

//----------------------------------------------------------------------------
// eta sub
void StTriFlowProManger::FillEtaCharged2Flow(Float_t pt, Float_t v2_EP, Float_t Res2_EP, Float_t v2_SP, Float_t Res2_SP, Int_t Cent9, Int_t eta_gap, Double_t reweight)
{
  if(Res2_EP > 0.0)
  {
    p_mV2_EP[Cent9][eta_gap]->Fill(pt, v2_EP, reweight);
    p_mV2_EP[9][eta_gap]->Fill(pt, v2_EP, reweight);
  }
  if(Res2_SP > 0.0)
  {
    p_mV2_SP[Cent9][eta_gap]->Fill(pt, v2_SP, reweight);
    p_mV2_SP[9][eta_gap]->Fill(pt, v2_SP, reweight);
  }
}

void StTriFlowProManger::FillEtaCharged3Flow(Float_t pt, Float_t v3_EP, Float_t Res3_EP, Float_t v3_SP, Float_t Res3_SP, Int_t Cent9, Int_t eta_gap, Double_t reweight)
{
  if(Res3_EP > 0.0)
  {
    p_mV3_EP[Cent9][eta_gap]->Fill(pt, v3_EP, reweight);
    p_mV3_EP[9][eta_gap]->Fill(pt, v3_EP, reweight);
  }
  if(Res3_SP > 0.0)
  {
    p_mV3_SP[Cent9][eta_gap]->Fill(pt, v3_SP, reweight);
    p_mV3_SP[9][eta_gap]->Fill(pt, v3_SP, reweight);
  }
}

// random sub
void StTriFlowProManger::FillRanCharged2Flow(Float_t pt, Float_t v2_Ran_EP, Float_t Res2_Ran_EP, Float_t v2_Ran_SP, Float_t Res2_Ran_SP, Int_t Cent9, Double_t reweight)
{
  if(Res2_Ran_EP > 0.0)
  {
    p_mV2_Ran_EP[Cent9]->Fill(pt, v2_Ran_EP, reweight);
    p_mV2_Ran_EP[9]->Fill(pt, v2_Ran_EP, reweight);
  }
  if(Res2_Ran_SP > 0.0)
  {
    p_mV2_Ran_SP[Cent9]->Fill(pt, v2_Ran_SP, reweight);
    p_mV2_Ran_SP[9]->Fill(pt, v2_Ran_SP, reweight);
  }
}

void StTriFlowProManger::FillRanCharged3Flow(Float_t pt, Float_t v3_Ran_EP, Float_t Res3_Ran_EP, Float_t v3_Ran_SP, Float_t Res3_Ran_SP, Int_t Cent9, Double_t reweight)
{
  if(Res3_Ran_EP > 0.0)
  {
    p_mV3_Ran_EP[Cent9]->Fill(pt, v3_Ran_EP, reweight);
    p_mV3_Ran_EP[9]->Fill(pt, v3_Ran_EP, reweight);
  }
  if(Res3_Ran_SP > 0.0)
  {
    p_mV3_Ran_SP[Cent9]->Fill(pt, v3_Ran_SP, reweight);
    p_mV3_Ran_SP[9]->Fill(pt, v3_Ran_SP, reweight);
  }
}

//----------------------------------------------------------------------------

void StTriFlowProManger::WriteReCenter()
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

      p_mq3x_East_EP[i][j]->Write();
      p_mq3y_East_EP[i][j]->Write();
      p_mq3x_West_EP[i][j]->Write();
      p_mq3y_West_EP[i][j]->Write();

      // Scalar Product method
      p_mq2x_East_SP[i][j]->Write();
      p_mq2y_East_SP[i][j]->Write();
      p_mq2x_West_SP[i][j]->Write();
      p_mq2y_West_SP[i][j]->Write();

      p_mq3x_East_SP[i][j]->Write();
      p_mq3y_East_SP[i][j]->Write();
      p_mq3x_West_SP[i][j]->Write();
      p_mq3y_West_SP[i][j]->Write();
    }
    // Event Plane method
    p_mq2x_Full_EP[i]->Write();
    p_mq2y_Full_EP[i]->Write();
    p_mq3x_Full_EP[i]->Write();
    p_mq3y_Full_EP[i]->Write();

    // Scalar Product method
    p_mq2x_Full_SP[i]->Write();
    p_mq2y_Full_SP[i]->Write();
    p_mq3x_Full_SP[i]->Write();
    p_mq3y_Full_SP[i]->Write();
  }
}

//----------------------------------------------------------------------------

void StTriFlowProManger::WriteShift()
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

	p_mcos3_East_EP[i][j][k]->Write();
	p_msin3_East_EP[i][j][k]->Write();
	p_mcos3_West_EP[i][j][k]->Write();
	p_msin3_West_EP[i][j][k]->Write();

        // Scalor Product method
	p_mcos2_East_SP[i][j][k]->Write();
	p_msin2_East_SP[i][j][k]->Write();
	p_mcos2_West_SP[i][j][k]->Write();
	p_msin2_West_SP[i][j][k]->Write();

	p_mcos3_East_SP[i][j][k]->Write();
	p_msin3_East_SP[i][j][k]->Write();
	p_mcos3_West_SP[i][j][k]->Write();
	p_msin3_West_SP[i][j][k]->Write();
      }
    }
    for(Int_t k = 0; k < 5; k++) // Shift Order
    {
      // Event Plane method
      p_mcos2_Full_EP[i][k]->Write();
      p_msin2_Full_EP[i][k]->Write();
      p_mcos3_Full_EP[i][k]->Write();
      p_msin3_Full_EP[i][k]->Write();

      // Scalor Product method
      p_mcos2_Full_SP[i][k]->Write();
      p_msin2_Full_SP[i][k]->Write();
      p_mcos3_Full_SP[i][k]->Write();
      p_msin3_Full_SP[i][k]->Write();
    }
  }
}

//----------------------------------------------------------------------------

void StTriFlowProManger::WriteChargedFlow()
{
  for(Int_t i = 0; i < 10; i++)
  {
    for(Int_t j = 0; j < 4; j++)
    {
      p_mV2_EP[i][j]->Write();
      p_mV3_EP[i][j]->Write();
      p_mV2_SP[i][j]->Write();
      p_mV3_SP[i][j]->Write();
    }
    p_mV2_Ran_EP[i]->Write();
    p_mV3_Ran_EP[i]->Write();
    p_mV2_Ran_SP[i]->Write();
    p_mV3_Ran_SP[i]->Write();
  }
}

//----------------------------------------------------------------------------
