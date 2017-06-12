#include "StRoot/StZdcSmdAna/StZdcSmdHistoManger.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "TString.h"
#include "TMath.h"
#include "TH1F.h"

ClassImp(StZdcSmdHistoManger)

//-------------------------------------------------------------
StZdcSmdHistoManger::StZdcSmdHistoManger()
{
}

StZdcSmdHistoManger::~StZdcSmdHistoManger()
{
}
//-------------------------------------------------------------

void StZdcSmdHistoManger::Init(int X_flag, int mode) // 0 for Same Event, 1 for Mixed Event
{
  TString Mode[2] = {"SE","ME"};
  // spin alignment analysis
  for(int i_pt = 0; i_pt < vmsa::pt_total; i_pt++) // pt bin 
  {
    for(int i_theta = 0; i_theta < vmsa::CTS_total; i_theta++) // phi-psi bin
    {
      for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // centrality bin
      {
	TString KEY_Mass2 = Form("pt_%d_Centrality_%d_CosThetaStar_%d_1st_%s_%s",i_pt,i_cent,i_theta,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
	h_mMass2[KEY_Mass2] = new TH1F(KEY_Mass2.Data(),KEY_Mass2.Data(),200,vmsa::InvMass_low[mode],vmsa::InvMass_high[mode]);
	h_mMass2[KEY_Mass2]->Sumw2();
      }
      for(int i_cent = 0; i_cent < 9; i_cent++) // centrality bin
      {
	TString KEY_Mass2_Cent = Form("pt_%d_Cent_%d_CosThetaStar_%d_1st_%s_%s",i_pt,i_cent,i_theta,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
	h_mMass_Cent[KEY_Mass2_Cent] = new TH1F(KEY_Mass2_Cent.Data(),KEY_Mass2_Cent.Data(),200,vmsa::InvMass_low[mode],vmsa::InvMass_high[mode]);
	h_mMass_Cent[KEY_Mass2_Cent]->Sumw2();
      }
    }
  }

  // raw pt spectra 
  for(int i_pt = 0; i_pt < vmsa::pt_total; i_pt++) // pt bin
  {
    for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // centrality bin
    {
      for(int i_pT = 0; i_pT < 2; i_pT++)
      {
	TString pT[2] = {"low","high"};
	TString KEY_Mass2_Spec = Form("Spec_pt_%d_%s_Centrality_%d_%s_%s",i_pt,pT[i_pT].Data(),i_cent,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
	h_mMass_Spec[KEY_Mass2_Spec] = new TH1F(KEY_Mass2_Spec.Data(),KEY_Mass2_Spec.Data(),200,vmsa::InvMass_low[mode],vmsa::InvMass_high[mode]);
	h_mMass_Spec[KEY_Mass2_Spec]->Sumw2();
      }
    }
  }
  
  // Yields
  for(int i_cent = 0; i_cent < 9; i_cent++) // centrality bin
  {
    TString KEY_Mass2_Yields = Form("Yields_Centrality_%d_%s_%s",i_cent,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
    h_mMass_Yields[KEY_Mass2_Yields] = new TH1F(KEY_Mass2_Yields.Data(),KEY_Mass2_Yields.Data(),200,vmsa::InvMass_low[mode],vmsa::InvMass_high[mode]);
    h_mMass_Yields[KEY_Mass2_Yields]->Sumw2();
  }
}
//-------------------------------------------------------------
void StZdcSmdHistoManger::Fill(float pt, int Cent9, float CosThetaStar, float Res, float InvMass, double reweight, int X_flag, int mode)
{
  TString Mode[2] = {"SE","ME"};
  if(Res > 0.0)
  {
    for(int i_pt = 0; i_pt < vmsa::pt_total; i_pt++) // pt_bin
    {
      if(pt >= vmsa::ptRawStart[i_pt] && pt < vmsa::ptRawStop[i_pt])
      {
	for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // centrality bin
	{
	  if(Cent9 >= vmsa::cent_low[i_cent] && Cent9 <= vmsa::cent_up[i_cent])
	  {
	    for(int i_theta = 0; i_theta < vmsa::CTS_total; i_theta++) // phi-psi2 bin
	    {
	      if(TMath::Abs(CosThetaStar) >= vmsa::CTS_low[i_theta] && TMath::Abs(CosThetaStar) < vmsa::CTS_up[i_theta])
	      {
		TString KEY_Mass2 = Form("pt_%d_Centrality_%d_CosThetaStar_%d_1st_%s_%s",i_pt,i_cent,i_theta,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
		h_mMass2[KEY_Mass2]->Fill(InvMass,reweight);
		// raw pt spectra
		if(pt < 0.5*(vmsa::ptRawStart[i_pt]+vmsa::ptRawStop[i_pt])) 
		{
		  TString KEY_Mass2_Spec = Form("Spec_pt_%d_low_Centrality_%d_%s_%s",i_pt,i_cent,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
		  h_mMass_Spec[KEY_Mass2_Spec]->Fill(InvMass,reweight);
		}
		else
		{
		  TString KEY_Mass2_Spec = Form("Spec_pt_%d_high_Centrality_%d_%s_%s",i_pt,i_cent,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
		  h_mMass_Spec[KEY_Mass2_Spec]->Fill(InvMass,reweight);
		}
	      }
	    }
	  }
	}
      }
    }
    for(int i_pt = 0; i_pt < vmsa::pt_total; i_pt++) // pt_bin
    {
      if(pt >= vmsa::ptRawStart[i_pt] && pt < vmsa::ptRawStop[i_pt])
      {
	for(int i_theta = 0; i_theta < vmsa::CTS_total; i_theta++) // phi-psi2 bin
	{
	  if(TMath::Abs(CosThetaStar) >= vmsa::CTS_low[i_theta] && TMath::Abs(CosThetaStar) < vmsa::CTS_up[i_theta])
	  {
	    TString KEY_Mass2_Cent = Form("pt_%d_Cent_%d_CosThetaStar_%d_1st_%s_%s",i_pt,Cent9,i_theta,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
	    h_mMass_Cent[KEY_Mass2_Cent]->Fill(InvMass,reweight);
	  }
	}
      }
    }
  }
  TString KEY_Mass2_Yields = Form("Yields_Centrality_%d_%s_%s",Cent9,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
  h_mMass_Yields[KEY_Mass2_Yields]->Fill(InvMass,reweight);
}
//-------------------------------------------------------------
void StZdcSmdHistoManger::Write(int X_flag, int mode)
{
  TString Mode[2] = {"SE","ME"};
  // flow
  for(int i_pt = 0; i_pt < vmsa::pt_total; i_pt++) // pt bin
  {
    for(int i_theta = 0; i_theta < vmsa::CTS_total; i_theta ++) // cos(theta*) bin
    {
      for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // centrality bin
      {
	TString KEY_Mass2 = Form("pt_%d_Centrality_%d_CosThetaStar_%d_1st_%s_%s",i_pt,i_cent,i_theta,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
	h_mMass2[KEY_Mass2]->Write();
      }
      for(int i_cent = 0; i_cent < 9; i_cent++) // centrality bin
      {
	TString KEY_Mass2_Cent = Form("pt_%d_Cent_%d_CosThetaStar_%d_1st_%s_%s",i_pt,i_cent,i_theta,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
	h_mMass_Cent[KEY_Mass2_Cent]->Write();
      }
    }
  }

  // Yields
  for(int i_cent = 0; i_cent < 9; i_cent++) // centrality bin
  {
    TString KEY_Mass2_Yields = Form("Yields_Centrality_%d_%s_%s",i_cent,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
    h_mMass_Yields[KEY_Mass2_Yields]->Write();
  }

  // raw pt spectra
  for(int i_pt = 0; i_pt < vmsa::pt_total; i_pt++) // pt bin
  {
    for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // centrality bin
    {
      for(int i_pT = 0; i_pT < 2; i_pT++)
      {
	TString pT[2] = {"low","high"};
	TString KEY_Mass2_Spec = Form("Spec_pt_%d_%s_Centrality_%d_%s_%s",i_pt,pT[i_pT].Data(),i_cent,vmsa::mPID[mode].c_str(),Mode[X_flag].Data());
	h_mMass_Spec[KEY_Mass2_Spec]->Write();
      }
    }
  }
}

//-------------------------------------------------------------
