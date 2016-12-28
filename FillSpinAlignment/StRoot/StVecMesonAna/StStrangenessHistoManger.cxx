#include "StStrangenessHistoManger.h"
#include "StStrangenessCons.h"
#include "TString.h"
#include "TMath.h"
#include "TH1F.h"

ClassImp(StStrangenessHistoManger)

//-------------------------------------------------------------
StStrangenessHistoManger::StStrangenessHistoManger()
{
}

StStrangenessHistoManger::~StStrangenessHistoManger()
{
}
//-------------------------------------------------------------

void StStrangenessHistoManger::Init(Int_t X_flag, Int_t mode) // 0 for Same Event, 1 for Mixed Event
{
  // flow analysis
  for(Int_t i_pt = 0; i_pt < Strangeness::pt_total_phi; i_pt++) // pt bin 
  {
    for(Int_t i_cent = Strangeness::Centrality_start; i_cent < Strangeness::Centrality_stop; i_cent++) // centrality bin
    {
      for(Int_t i_eta = Strangeness::EtaGap_start; i_eta < Strangeness::EtaGap_stop; i_eta++) // eta gap bin
      {
	for(Int_t i_theta = 0; i_theta < Strangeness::CosThetaStar_total; i_theta++) // phi-psi bin
	{
	  TString Mode[2] = {"SE","ME"};
	  TString KEY_Mass2 = Form("pt_%d_Centrality_%d_EtaGap_%d_CosThetaStar_%d_2nd_%s_%s",i_pt,i_cent,i_eta,i_theta,Strangeness::Partype[mode].Data(),Mode[X_flag].Data());
	  h_mMass2[KEY_Mass2] = new TH1F(KEY_Mass2.Data(),KEY_Mass2.Data(),200,Strangeness::InvMass_low[mode],Strangeness::InvMass_high[mode]);
	  h_mMass2[KEY_Mass2]->Sumw2();
	}
      }
    }
  }

  // raw pt spectra 
  for(Int_t i_pt = 0; i_pt < Strangeness::pt_total_phi; i_pt++) // pt bin
  {
    for(Int_t i_cent = Strangeness::Centrality_start; i_cent < Strangeness::Centrality_stop; i_cent++) // centrality bin
    {
      for(Int_t i_eta = Strangeness::EtaGap_start; i_eta < Strangeness::EtaGap_stop; i_eta++) // eta gap bin
      {
	for(Int_t i_pT = 0; i_pT < 2; i_pT++)
	{
	  TString Mode[2] = {"SE","ME"};
	  TString pT[2] = {"low","high"};
	  TString KEY_Mass2_Spec = Form("Spec_pt_%d_%s_Centrality_%d_EtaGap_%d_%s_%s",i_pt,pT[i_pT].Data(),i_cent,i_eta,Strangeness::Partype[mode].Data(),Mode[X_flag].Data());
	  h_mMass_Spec[KEY_Mass2_Spec] = new TH1F(KEY_Mass2_Spec.Data(),KEY_Mass2_Spec.Data(),200,Strangeness::InvMass_low[mode],Strangeness::InvMass_high[mode]);
	  h_mMass_Spec[KEY_Mass2_Spec]->Sumw2();
	}
      }
    }
  }
  
  // Yields
  for(Int_t i_cent = 0; i_cent < 9; i_cent++) // centrality bin
  {
    for(Int_t i_eta = Strangeness::EtaGap_start; i_eta < Strangeness::EtaGap_stop; i_eta++) // eta gap bin
    {
      TString Mode[2] = {"SE","ME"};
      TString KEY_Mass2_Yields = Form("Yields_Centrality_%d_EtaGap_%d_%s_%s",i_cent,i_eta,Strangeness::Partype[mode].Data(),Mode[X_flag].Data());
      h_mMass_Yields[KEY_Mass2_Yields] = new TH1F(KEY_Mass2_Yields.Data(),KEY_Mass2_Yields.Data(),200,Strangeness::InvMass_low[mode],Strangeness::InvMass_high[mode]);
      h_mMass_Yields[KEY_Mass2_Yields]->Sumw2();
    }
  }
}
//-------------------------------------------------------------
void StStrangenessHistoManger::Fill(Float_t pt, Int_t Cent9, Int_t eta_gap, Float_t CosThetaStar, Float_t Res2, Float_t InvMass, Double_t reweight, Int_t X_flag, Int_t mode)
{
  TString Mode[2] = {"SE","ME"};
  if(Res2 > 0.0)
  {
    for(Int_t i_pt = 0; i_pt < Strangeness::pt_total_phi; i_pt++) // pt_bin
    {
      if(pt >= Strangeness::pt_low[i_pt] && pt < Strangeness::pt_up[i_pt])
      {
	for(Int_t i_cent = Strangeness::Centrality_start; i_cent < Strangeness::Centrality_stop; i_cent++) // centrality bin
	{
	  if(Cent9 >= Strangeness::cent_low[i_cent] && Cent9 <= Strangeness::cent_up[i_cent])
	  {
	    for(Int_t i_theta = 0; i_theta < Strangeness::CosThetaStar_total; i_theta++) // phi-psi2 bin
	    {
	      if(TMath::Abs(CosThetaStar) >= Strangeness::CosThetaStar_low[i_theta] && TMath::Abs(CosThetaStar) < Strangeness::CosThetaStar_up[i_theta])
	      {
		// flow
		TString KEY_Mass2 = Form("pt_%d_Centrality_%d_EtaGap_%d_CosThetaStar_%d_2nd_%s_%s",i_pt,i_cent,eta_gap,i_theta,Strangeness::Partype[mode].Data(),Mode[X_flag].Data());
		h_mMass2[KEY_Mass2]->Fill(InvMass,reweight);
		// raw pt spectra
		if(pt < 0.5*(Strangeness::pt_low[i_pt]+Strangeness::pt_up[i_pt])) 
		{
		  TString KEY_Mass2_Spec = Form("Spec_pt_%d_low_Centrality_%d_EtaGap_%d_%s_%s",i_pt,i_cent,eta_gap,Strangeness::Partype[mode].Data(),Mode[X_flag].Data());
		  h_mMass_Spec[KEY_Mass2_Spec]->Fill(InvMass,reweight);
		}
		else
		{
		  TString KEY_Mass2_Spec = Form("Spec_pt_%d_high_Centrality_%d_EtaGap_%d_%s_%s",i_pt,i_cent,eta_gap,Strangeness::Partype[mode].Data(),Mode[X_flag].Data());
		  h_mMass_Spec[KEY_Mass2_Spec]->Fill(InvMass,reweight);
		}
	      }
	    }
	  }
	}
      }
    }
  }
  TString KEY_Mass2_Yields = Form("Yields_Centrality_%d_EtaGap_%d_%s_%s",Cent9,eta_gap,Strangeness::Partype[mode].Data(),Mode[X_flag].Data());
  h_mMass_Yields[KEY_Mass2_Yields]->Fill(InvMass,reweight);
}
//-------------------------------------------------------------
void StStrangenessHistoManger::Write(Int_t X_flag, Int_t mode)
{
  TString Mode[2] = {"SE","ME"};
  // flow
  for(Int_t i_pt = 0; i_pt < Strangeness::pt_total_phi; i_pt++) // pt bin
  {
    for(Int_t i_cent = Strangeness::Centrality_start; i_cent < Strangeness::Centrality_stop; i_cent++) // centrality bin
    {
      for(Int_t i_eta = Strangeness::EtaGap_start; i_eta < Strangeness::EtaGap_stop; i_eta++) // eta gap bin
      {
	for(Int_t i_theta = 0; i_theta < Strangeness::CosThetaStar_total; i_theta ++) // cos(theta*) bin
	{
	  TString KEY_Mass2 = Form("pt_%d_Centrality_%d_EtaGap_%d_CosThetaStar_%d_2nd_%s_%s",i_pt,i_cent,i_eta,i_theta,Strangeness::Partype[mode].Data(),Mode[X_flag].Data());
	  h_mMass2[KEY_Mass2]->Write();
	}
      }
    }
  }

  // Yields
  for(Int_t i_cent = 0; i_cent < 9; i_cent++) // centrality bin
  {
    for(Int_t i_eta = Strangeness::EtaGap_start; i_eta < Strangeness::EtaGap_stop; i_eta++) // eta gap bin
    {
      TString KEY_Mass2_Yields = Form("Yields_Centrality_%d_EtaGap_%d_%s_%s",i_cent,i_eta,Strangeness::Partype[mode].Data(),Mode[X_flag].Data());
      h_mMass_Yields[KEY_Mass2_Yields]->Write();
    }
  }

  // raw pt spectra
  for(Int_t i_pt = 0; i_pt < Strangeness::pt_total_phi; i_pt++) // pt bin
  {
    for(Int_t i_cent = Strangeness::Centrality_start; i_cent < Strangeness::Centrality_stop; i_cent++) // centrality bin
    {
      for(Int_t i_eta = Strangeness::EtaGap_start; i_eta < Strangeness::EtaGap_stop; i_eta++) // eta gap bin
      {
	for(Int_t i_pT = 0; i_pT < 2; i_pT++)
	{
	  TString pT[2] = {"low","high"};
	  TString KEY_Mass2_Spec = Form("Spec_pt_%d_%s_Centrality_%d_EtaGap_%d_%s_%s",i_pt,pT[i_pT].Data(),i_cent,i_eta,Strangeness::Partype[mode].Data(),Mode[X_flag].Data());
	  h_mMass_Spec[KEY_Mass2_Spec]->Write();
	}
      }
    }
  }
}
