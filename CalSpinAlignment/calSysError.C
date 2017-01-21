#include <string>
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TF1.h"
#include "TStyle.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "../Utility/type.h"
#include "../Utility/draw.h"
#include "../Utility/functions.h"

#ifndef _PlotQA_
#define _PlotQA_ 1
#endif

void calSysError(Int_t energy = 6, Int_t pid = 0)
{
  string inputfile = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/%s/rho00/ResRhoPtSys.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  TGraMap g_mRho;
  for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; ++i_cent) // Centrality loop
  {
    for(int i_eta = vmsa::Eta_start; i_eta < vmsa::Eta_stop; ++i_eta) // EtaGap loop
    {
      for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
      {
	for(int i_func = vmsa::Func_start; i_func < vmsa::Func_stop; ++i_func)
	{
	  for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
	  {
	    for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
	    {
	      string KEY_rho = Form("rhoRaw_Centrality_%d_EtaGap_%d_2nd_%s_Norm_%d_Func_%d_Sigma_%d_%s",i_cent,i_eta,vmsa::mPID[pid].c_str(),i_norm,i_func,i_sigma,vmsa::mInteMethod[i_method].c_str());
	      g_mRho[KEY_rho] = (TGraphAsymmErrors*)File_InPut->Get(KEY_rho.c_str());
	    }
	  }
	}
      }
    }
  }
  TH1F *h_frame = (TH1F*)File_InPut->Get("h_frame");

#if _PlotQA_
  TCanvas *c_rho = new TCanvas("c_rho","c_rho",10,10,800,800);
  c_rho->cd();
  c_rho->cd()->SetLeftMargin(0.15);
  c_rho->cd()->SetBottomMargin(0.15);
  c_rho->cd()->SetTicks(1,1);
  c_rho->cd()->SetGrid(0,0);
  h_frame->DrawCopy("pE");
  PlotLine(0.0,5.0,1.0/3.0,1.0/3.0,1,2,2);

  for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; ++i_cent) // Centrality loop
  {
    for(int i_eta = vmsa::Eta_start; i_eta < vmsa::Eta_stop; ++i_eta) // EtaGap loop
    {
      for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
      {
	for(int i_func = vmsa::Func_start; i_func < vmsa::Func_stop; ++i_func)
	{
	  for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
	  {
	    for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
	    {
	      string KEY_rho = Form("rhoRaw_Centrality_%d_EtaGap_%d_2nd_%s_Norm_%d_Func_%d_Sigma_%d_%s",i_cent,i_eta,vmsa::mPID[pid].c_str(),i_norm,i_func,i_sigma,vmsa::mInteMethod[i_method].c_str());
	      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mRho[KEY_rho],24,i_sigma+10*i_method+1,1.1);
	    }
	  }
	}
      }
    }
  }
#endif

  string KEY_Default = Form("rhoRaw_Centrality_0_EtaGap_0_2nd_%s_Norm_1_Func_0_Sigma_2_Count",vmsa::mPID[pid].c_str());

  TGraphAsymmErrors *g_SysErrors = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_StatErrors = new TGraphAsymmErrors();
  for(Int_t i_point = 0; i_point < g_mRho[KEY_Default]->GetN(); ++i_point)
  {
    int counter = 0;
    float total_rho = 0.0; 
    TGraph *g_diff = new TGraph();
    for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; ++i_cent) // Centrality loop
    {
      for(int i_eta = vmsa::Eta_start; i_eta < vmsa::Eta_stop; ++i_eta) // EtaGap loop
      {
	for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
	{
	  for(int i_func = vmsa::Func_start; i_func < vmsa::Func_stop; ++i_func)
	  {
	    for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
	    {
	      for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
	      {
		string KEY_rho = Form("rhoRaw_Centrality_%d_EtaGap_%d_2nd_%s_Norm_%d_Func_%d_Sigma_%d_%s",i_cent,i_eta,vmsa::mPID[pid].c_str(),i_norm,i_func,i_sigma,vmsa::mInteMethod[i_method].c_str());
		Double_t pt_val, rho_val;
		g_mRho[KEY_rho]->GetPoint(i_point,pt_val,rho_val);
		g_diff->SetPoint(counter,pt_val-0.1,rho_val);
		total_rho += rho_val;
		counter++;
	      }
	    }
	  }
	}
      }
    }
    Double_t rho_min = TMath::MinElement(g_diff->GetN(),g_diff->GetY());
    Double_t rho_max = TMath::MaxElement(g_diff->GetN(),g_diff->GetY());
    cout << "rho_min = " << rho_min << ", rho_max = " << rho_max << endl;
    Double_t pt, rho;
    g_mRho[KEY_Default]->GetPoint(i_point,pt,rho);

    float SysError_rho = (rho_max-rho_min)/TMath::Sqrt(12.0);
    float mean_rho = total_rho/(float)counter;
    g_SysErrors->SetPoint(i_point,pt,mean_rho);
    g_SysErrors->SetPointError(i_point,0.0,0.0,SysError_rho,SysError_rho);

    float StatError_rho = g_mRho[KEY_Default]->GetErrorYhigh(i_point);
    g_StatErrors->SetPoint(i_point,pt,mean_rho);
    g_StatErrors->SetPointError(i_point,0.0,0.0,StatError_rho,StatError_rho);
    cout << "number of combinations is: " << counter << endl;
    delete g_diff;
  }

  TCanvas *c_rho_SysError = new TCanvas("c_rho_SysError","c_rho_SysError",600,10,800,800);
  c_rho_SysError->cd();
  c_rho_SysError->cd()->SetLeftMargin(0.15);
  c_rho_SysError->cd()->SetBottomMargin(0.15);
  c_rho_SysError->cd()->SetTicks(1,1);
  c_rho_SysError->cd()->SetGrid(0,0);
  h_frame->Draw("pE");
  g_StatErrors->SetMarkerStyle(20);
  g_StatErrors->SetMarkerColor(kGray+2);
  g_StatErrors->SetLineColor(2);
  g_StatErrors->Draw("pE same");
  g_SysErrors->SetMarkerStyle(20);
  g_SysErrors->SetMarkerColor(kGray+2);
  g_SysErrors->SetLineColor(kGray+2);
  g_SysErrors->Draw("pE same");
  PlotLine(0.0,5.0,1.0/3.0,1.0/3.0,1,2,2);

  string OutPutFile = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/%s/rho00/Rho_SysErrors.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  cout << "OutPutFile set to: " << OutPutFile.c_str() << endl;
  TFile *File_OutPut = new TFile(OutPutFile.c_str(),"RECREATE");
  File_OutPut->cd();
  h_frame->Write();
  TString StatErrorRho = Form("g_rho00_%s_%s_StatError",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  g_StatErrors->SetName(StatErrorRho.Data());
  g_StatErrors->Write();
  TString SysErrorRho = Form("g_rho00_%s_%s_SysError",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  g_SysErrors->SetName(SysErrorRho.Data());
  g_SysErrors->Write();
  File_OutPut->Close();
}
