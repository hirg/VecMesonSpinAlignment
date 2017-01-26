#include <iostream>
#include <map>
#include <vector>
#include <string>
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TMath.h"
#include "TLegend.h"
#include "TGaxis.h"
#include "TProfile.h"
#include "TGraphAsymmErrors.h"
#include "TProfile2D.h"
#include "../Utility/functions.h"
#include "../Utility/draw.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "../Utility/type.h"

#ifndef _PlotQA_
#define _PlotQA_  1
#endif

TF1* getSgFunc(int nFunc, int pid);
TF1* getBgFunc(int nFunc, int pid);

static float const FuncPar[3][4][3] = { // 0 for norm | 1 for funcs | 2 for parameters
  {{-1,1,0.0},{0.1,-1.4,1.4},{0,-1,1},{0,1.5,0.2}}, // left
  {{-1,1,0.0},{0.0,-1.0,1.5},{0,-1,1},{0,1.5,0.2}}, // right
  {{-1,1,0.0},{0.0,-1.0,1.5},{0,-1,1},{0,1.5,0.2}}  // both
}; // 19GeV
// static float const FuncPar[3][4][3] = { // 0 for norm | 1 for funcs | 2 for parameters
//   {{-1,1,0.0},{0.0,-1.4,1.4},{0,-1,1},{0,1.5,0.2}}, // left
//   {{-1,1,0.0},{0.0,-1.5,1.5},{0,-1,1},{0,1.5,0.2}}, // right
//   {{-1,1,0.0},{1.0,-1.5,1.5},{0,-1,1},{0,1.5,0.2}}  // both
// }; // 27GeV
// static float const FuncPar[3][4][3] = { // 0 for norm | 1 for funcs | 2 for parameters
//   {{-1,1,0.0},{0.2,-1.45,1.35},{0,-1,1},{0,1.5,0.2}}, // left
//   {{-1,1,0.0},{0.0,-1.5,1.5},{0,-1,1},{0,1.5,0.2}}, // right
//   {{-1,1,0.0},{1.0,-1.4,1.4},{0,-1,1},{0,1.5,0.2}}  // both
// }; // 39GeV
// static float const FuncPar[3][4][3] = { // 0 for norm | 1 for funcs | 2 for parameters
//   {{-1,1,0.0},{0.0,-1.40,1.5},{0,-1,1},{0,1.5,0.2}}, // left
//   {{-1,1,0.0},{0.0,-1.50,1.5},{0,-1,1},{0,1.5,0.2}}, // right
//   {{-1,1,0.0},{0.0,-1.45,1.6},{0,-1,1},{0,1.5,0.2}}  // both
// }; // 62GeV
// static float const FuncPar[3][4][3] = { // 0 for norm | 1 for funcs | 2 for parameters
//   {{-1,1,0.0},{0.1,-1.6,1.3},{0,-1,1},{0,1.5,0.2}}, // left
//   {{-1,1,0.0},{0.2,-1.6,1.3},{0,-1,1},{0,1.5,0.2}}, // right
//   {{-1,1,0.0},{0.0,-1.6,1.9},{0,-1,1},{0,1.5,0.2}}  // both
// }; // 200GeV

static float const normL = 5.0;
static float const normR = 5.0;
void subBackGround(int energy = 3, int pid = 0, int year = 0)
{
  TGaxis::SetMaxDigits(4);
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(50000);

  string InPutFile_SE = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/%s/Yields/merged_file/Yields_SE_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str());
  
  string InPutFile_ME = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/%s/Yields/merged_file/Yields_ME_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str());
  TFile *File_SE = TFile::Open(InPutFile_SE.c_str());
  TFile *File_ME = TFile::Open(InPutFile_ME.c_str());

  // read in histogram for same event and mixed event
  // calculate SE - ME
  TH1FMap h_mInPut_SE, h_mInPut_ME;
  TH1FMap h_mMass_SE, h_mMass_ME, h_mMass_SM;
  for(int i_pt = vmsa::pt_start; i_pt < vmsa::pt_stop; i_pt++) // pt loop
  {
    for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // Centrality loop
    {
      for(int i_eta = vmsa::Eta_start; i_eta < vmsa::Eta_stop; i_eta++) // EtaGap loop
      {
	for(int i_theta = vmsa::CTS_start; i_theta < vmsa::CTS_stop; i_theta++) // cos(theta*) loop
	{
	  string KEY_InPutSE = Form("pt_%d_Centrality_%d_EtaGap_%d_CosThetaStar_%d_2nd_%s_SE",i_pt,i_cent,i_eta,i_theta,vmsa::mPID[pid].c_str());
	  h_mInPut_SE[KEY_InPutSE] = (TH1F*)File_SE->Get(KEY_InPutSE.c_str())->Clone(); 

	  string KEY_InPutME = Form("pt_%d_Centrality_%d_EtaGap_%d_CosThetaStar_%d_2nd_%s_ME",i_pt,i_cent,i_eta,i_theta,vmsa::mPID[pid].c_str());
	  h_mInPut_ME[KEY_InPutME] = (TH1F*)File_ME->Get(KEY_InPutME.c_str())->Clone(); 

	  for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
	  {
	    string KEY_SE = Form("pt_%d_Centrality_%d_EtaGap_%d_CosThetaStar_%d_2nd_%s_Norm_%d_SE",i_pt,i_cent,i_eta,i_theta,vmsa::mPID[pid].c_str(),i_norm);
	    h_mMass_SE[KEY_SE] = (TH1F*)h_mInPut_SE[KEY_InPutSE]->Clone(KEY_SE.c_str());
	    string KEY_ME = Form("pt_%d_Centrality_%d_EtaGap_%d_CosThetaStar_%d_2nd_%s_Norm_%d_ME",i_pt,i_cent,i_eta,i_theta,vmsa::mPID[pid].c_str(),i_norm);
	    h_mMass_ME[KEY_ME] = (TH1F*)h_mInPut_ME[KEY_InPutME]->Clone(KEY_ME.c_str());
	    string KEY_SM = Form("pt_%d_Centrality_%d_EtaGap_%d_CosThetaStar_%d_2nd_%s_Norm_%d_SM",i_pt,i_cent,i_eta,i_theta,vmsa::mPID[pid].c_str(),i_norm);
	    h_mMass_SM[KEY_SM] = (TH1F*)h_mMass_SE[KEY_SE]->Clone();
	    if(i_norm < 2)
	    {
	      int Norm_bin_start = h_mMass_SE[KEY_SE]->FindBin(vmsa::Norm_Start[pid][i_norm]);
	      int Norm_bin_stop  = h_mMass_SE[KEY_SE]->FindBin(vmsa::Norm_Stop[pid][i_norm]);

	      float Inte_SE = h_mMass_SE[KEY_SE]->Integral(Norm_bin_start,Norm_bin_stop);
	      float Inte_ME = h_mMass_ME[KEY_ME]->Integral(Norm_bin_start,Norm_bin_stop);

	      h_mMass_ME[KEY_ME]->Scale(Inte_SE/Inte_ME);
	      h_mMass_SM[KEY_SM]->Add(h_mMass_ME[KEY_ME],-1.0);
	    }
	    else
	    {
	      float Inte_SE = 0.0;
	      float Inte_ME = 0.0;

	      for(int i_inte = 0; i_inte < 2; ++i_inte)
	      {
		int Norm_bin_start = h_mMass_SE[KEY_SE]->FindBin(vmsa::Norm_Start[pid][i_inte]);
		int Norm_bin_stop  = h_mMass_SE[KEY_SE]->FindBin(vmsa::Norm_Stop[pid][i_inte]);
		Inte_SE += h_mMass_SE[KEY_SE]->Integral(Norm_bin_start,Norm_bin_stop);
		Inte_ME += h_mMass_ME[KEY_ME]->Integral(Norm_bin_start,Norm_bin_stop);
	      }

	      h_mMass_ME[KEY_ME]->Scale(Inte_SE/Inte_ME);
	      h_mMass_SM[KEY_SM]->Add(h_mMass_ME[KEY_ME],-1.0);
	    }
	  }
	}
      }
    }
  }

#if _PlotQA_
  // QA Plots for SE vs. ME
  TCanvas *c_peak = new TCanvas("c_peak","c_peak",10,10,800,800);
  c_peak->cd();
  c_peak->cd()->SetLeftMargin(0.15);
  c_peak->cd()->SetBottomMargin(0.15);
  c_peak->cd()->SetTicks(1,1);
  c_peak->cd()->SetGrid(0,0);
  string KEY_SE_QA = Form("pt_%d_Centrality_%d_EtaGap_%d_CosThetaStar_%d_2nd_%s_Norm_%d_SE",vmsa::pt_RawQA[energy],vmsa::Cent_start,vmsa::Eta_start,vmsa::CTS_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
  h_mMass_SE[KEY_SE_QA]->GetYaxis()->SetRangeUser(-0.1*h_mMass_SE[KEY_SE_QA]->GetMaximum(),1.1*h_mMass_SE[KEY_SE_QA]->GetMaximum());
  h_mMass_SE[KEY_SE_QA]->DrawCopy("PE");

  string KEY_ME_QA = Form("pt_%d_Centrality_%d_EtaGap_%d_CosThetaStar_%d_2nd_%s_Norm_%d_ME",vmsa::pt_RawQA[energy],vmsa::Cent_start,vmsa::Eta_start,vmsa::CTS_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
  h_mMass_ME[KEY_ME_QA]->SetLineColor(2);
  h_mMass_ME[KEY_ME_QA]->SetFillColor(2);
  h_mMass_ME[KEY_ME_QA]->SetFillStyle(3002);
  h_mMass_ME[KEY_ME_QA]->DrawCopy("h same");

  string KEY_SM_QA = Form("pt_%d_Centrality_%d_EtaGap_%d_CosThetaStar_%d_2nd_%s_Norm_%d_SM",vmsa::pt_RawQA[energy],vmsa::Cent_start,vmsa::Eta_start,vmsa::CTS_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
  h_mMass_SM[KEY_SM_QA]->SetLineColor(4);
  h_mMass_SM[KEY_SM_QA]->SetFillColor(4);
  h_mMass_SM[KEY_SM_QA]->SetFillStyle(3004);
  h_mMass_SM[KEY_SM_QA]->DrawCopy("h same");

  if(vmsa::Norm_QA == 0 || vmsa::Norm_QA == 2)
  {
    PlotLine(vmsa::Norm_Start[pid][0],vmsa::Norm_Start[pid][0],0,h_mMass_ME[KEY_ME_QA]->GetMaximum(),4,2,2);
    PlotLine(vmsa::Norm_Stop[pid][0],vmsa::Norm_Stop[pid][0],0,h_mMass_ME[KEY_ME_QA]->GetMaximum(),4,2,2);
  }
  if(vmsa::Norm_QA == 1 || vmsa::Norm_QA == 2)
  {
    PlotLine(vmsa::Norm_Start[pid][1],vmsa::Norm_Start[pid][1],0,h_mMass_ME[KEY_ME_QA]->GetMaximum(),4,2,2);
    PlotLine(vmsa::Norm_Stop[pid][1],vmsa::Norm_Stop[pid][1],0,h_mMass_ME[KEY_ME_QA]->GetMaximum(),4,2,2);
  }


  // QA Plots for pT bins
  TCanvas *c_pT = new TCanvas("c_pT","c_pT",10,10,1400,1400);
  c_pT->Divide(5,5);
  for(int i_pt = vmsa::pt_start; i_pt < vmsa::pt_stop; i_pt++) // pt loop
  {
    c_pT->cd(i_pt+1);
    c_pT->cd(i_pt+1)->SetLeftMargin(0.15);
    c_pT->cd(i_pt+1)->SetBottomMargin(0.15);
    c_pT->cd(i_pt+1)->SetTicks(1,1);
    c_pT->cd(i_pt+1)->SetGrid(0,0);
    string KEY_SE_QA = Form("pt_%d_Centrality_%d_EtaGap_%d_CosThetaStar_%d_2nd_%s_Norm_%d_SE",i_pt,vmsa::Cent_start,vmsa::Eta_start,vmsa::CTS_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
    h_mMass_SE[KEY_SE_QA]->GetYaxis()->SetRangeUser(-0.1*h_mMass_SE[KEY_SE_QA]->GetMaximum(),1.1*h_mMass_SE[KEY_SE_QA]->GetMaximum());
    h_mMass_SE[KEY_SE_QA]->DrawCopy();

    string KEY_ME_QA = Form("pt_%d_Centrality_%d_EtaGap_%d_CosThetaStar_%d_2nd_%s_Norm_%d_ME",i_pt,vmsa::Cent_start,vmsa::Eta_start,vmsa::CTS_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
    h_mMass_ME[KEY_ME_QA]->SetLineColor(2);
    h_mMass_ME[KEY_ME_QA]->SetFillColor(2);
    h_mMass_ME[KEY_ME_QA]->SetFillStyle(3002);
    h_mMass_ME[KEY_ME_QA]->DrawCopy("h same");

    string KEY_SM_QA = Form("pt_%d_Centrality_%d_EtaGap_%d_CosThetaStar_%d_2nd_%s_Norm_%d_SM",i_pt,vmsa::Cent_start,vmsa::Eta_start,vmsa::CTS_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
    h_mMass_SM[KEY_SM_QA]->SetLineColor(4);
    h_mMass_SM[KEY_SM_QA]->SetFillColor(4);
    h_mMass_SM[KEY_SM_QA]->SetFillStyle(3004);
    h_mMass_SM[KEY_SM_QA]->DrawCopy("h same");

    string pT_range = Form("[%.2f,%.2f]",vmsa::ptRawStart[i_pt],vmsa::ptRawStop[i_pt]);
    plotTopLegend((char*)pT_range.c_str(),0.2,0.7,0.08,1,0.0,42,1);

    if(vmsa::Norm_QA == 0 || vmsa::Norm_QA == 2)
    {
      PlotLine(vmsa::Norm_Start[pid][0],vmsa::Norm_Start[pid][0],0,h_mMass_ME[KEY_ME_QA]->GetMaximum(),4,2,2);
      PlotLine(vmsa::Norm_Stop[pid][0],vmsa::Norm_Stop[pid][0],0,h_mMass_ME[KEY_ME_QA]->GetMaximum(),4,2,2);
    }
    if(vmsa::Norm_QA == 1 || vmsa::Norm_QA == 2)
    {
      PlotLine(vmsa::Norm_Start[pid][1],vmsa::Norm_Start[pid][1],0,h_mMass_ME[KEY_ME_QA]->GetMaximum(),4,2,2);
      PlotLine(vmsa::Norm_Stop[pid][1],vmsa::Norm_Stop[pid][1],0,h_mMass_ME[KEY_ME_QA]->GetMaximum(),4,2,2);
    }
  }
#endif

  // pT rebin
  TH1FMap h_mMass; // rebinned InvMass distribution, SE-ME

  for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // Centrality loop
  {
    for(int i_eta = vmsa::Eta_start; i_eta < vmsa::Eta_stop; i_eta++) // EtaGap loop
    {
      for(int i_theta = vmsa::CTS_start; i_theta < vmsa::CTS_stop; i_theta++) // cos(theta*) loop
      {
	for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
	{
	  for(int pt_bin = vmsa::pt_rebin_first[energy]; pt_bin < vmsa::pt_rebin_last[energy]; pt_bin++) // pt loop
	  {
	    string KEY = Form("pt_%d_Centrality_%d_EtaGap_%d_CosThetaStar_%d_2nd_%s_Norm_%d",pt_bin,i_cent,i_eta,i_theta,vmsa::mPID[pid].c_str(),i_norm);
	    for(int i_pt = vmsa::pt_rebin_start[energy][pt_bin]; i_pt <= vmsa::pt_rebin_stop[energy][pt_bin]; i_pt++)
	    {
	      string KEY_SM = Form("pt_%d_Centrality_%d_EtaGap_%d_CosThetaStar_%d_2nd_%s_Norm_%d_SM",i_pt,i_cent,i_eta,i_theta,vmsa::mPID[pid].c_str(),i_norm);
	      //	      cout << "KEY= " << KEY.c_str() << ", KEY_SM = " << KEY_SM.c_str() << endl;
	      if(i_pt == vmsa::pt_rebin_start[energy][pt_bin])
	      {
		h_mMass[KEY] = (TH1F*)h_mMass_SM[KEY_SM]->Clone();
	      }
	      else
	      {
		h_mMass[KEY]->Add(h_mMass_SM[KEY_SM],1.0);
	      }
	    }
	  }
	}
      }
    }
  }

#if _PlotQA_
  // QA Plots for pT rebins
  TCanvas *c_pT_rebin = new TCanvas("c_pT_rebin","c_pT_rebin",10,10,1400,1400);
  c_pT_rebin->Divide(5,5);
  for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; i_pt++) // pt loop
  {
    c_pT_rebin->cd(vmsa::pt_rebin_start[energy][i_pt]+1);
    c_pT_rebin->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetLeftMargin(0.15);
    c_pT_rebin->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetBottomMargin(0.15);
    c_pT_rebin->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetTicks(1,1);
    c_pT_rebin->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetGrid(0,0);
    string KEY_QA = Form("pt_%d_Centrality_%d_EtaGap_%d_CosThetaStar_%d_2nd_%s_Norm_%d",i_pt,vmsa::Cent_start,vmsa::Eta_start,vmsa::CTS_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
    h_mMass[KEY_QA]->SetMarkerStyle(24);
    h_mMass[KEY_QA]->SetMarkerSize(1.2);
    h_mMass[KEY_QA]->SetLineColor(kGray+2);
    h_mMass[KEY_QA]->DrawCopy("pE");
    // cout << "KEY_QA = " << KEY_QA.c_str() << endl;

    string pT_range = Form("[%.2f,%.2f]",vmsa::pt_low[energy][i_pt],vmsa::pt_up[energy][i_pt]);
    plotTopLegend((char*)pT_range.c_str(),0.2,0.7,0.08,1,0.0,42,1);
    PlotLine(vmsa::InvMass_low[pid],vmsa::InvMass_high[pid],0.0,0.0,1,2,2);
  }
#endif

  TH1FMap h_mMass_func;
  if(pid == 0 || pid == 1) // Polynomial fit subtraction is only needed for phi meson
  {
    // Poly + Breit Wignar fit to phi integrated InvMass
    TH1FMap h_mMass_theta;
    vecFMap ParFit_theta;
    for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; i_pt++) // pt loop
    {
      for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // Centrality loop
      {
	for(int i_eta = vmsa::Eta_start; i_eta < vmsa::Eta_stop; i_eta++) // EtaGap loop
	{
	  for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
	  {
	    for(int i_func = vmsa::Func_start; i_func < vmsa::Func_stop; ++i_func)
	    {
	      string KEY_theta = Form("pt_%d_Centrality_%d_EtaGap_%d_2nd_%s_Norm_%d_Func_%d",i_pt,i_cent,i_eta,vmsa::mPID[pid].c_str(),i_norm,i_func);
	      for(int i_theta = vmsa::CTS_start; i_theta < vmsa::CTS_stop; i_theta++) // cos(theta*) loop
	      {
		string KEY = Form("pt_%d_Centrality_%d_EtaGap_%d_CosThetaStar_%d_2nd_%s_Norm_%d",i_pt,i_cent,i_eta,i_theta,vmsa::mPID[pid].c_str(),i_norm);
		if(i_theta == vmsa::CTS_start) h_mMass_theta[KEY_theta] = (TH1F*)h_mMass[KEY]->Clone(KEY_theta.c_str());
		else h_mMass_theta[KEY_theta]->Add(h_mMass[KEY],1.0);
	      }

	      TF1 *f_sg = getSgFunc(i_func,pid);
	      f_sg->SetParameter(0,vmsa::InvMass[pid]);
	      f_sg->SetParLimits(0,vmsa::InvMass[pid]-1.5*vmsa::Width[pid],vmsa::InvMass[pid]+1.5*vmsa::Width[pid]);
	      f_sg->SetParameter(1,vmsa::Width[pid]);
	      f_sg->SetParLimits(1,0.003,0.006);
	      f_sg->SetParameter(2,1);
	      for(int i_par = 3; i_par < vmsa::FuncParNum[i_func]; ++i_par)
	      {
		f_sg->SetParameter(i_par,0.0);
	      }
	      f_sg->SetParameter(2,h_mMass_theta[KEY_theta]->GetMaximum()/f_sg->GetMaximum());
	      for(int i_par = 3; i_par < vmsa::FuncParNum[i_func]; ++i_par)
	      {
		f_sg->SetParameter(i_par,FuncPar[i_norm][i_func][i_par-3]);
		// cout << "i_par = " << i_par << ", FuncPar = " << FuncPar[i_norm][i_func][i_par-3] << endl;
	      }
	      // if(energy == 4 && i_pt == 2 && i_norm == 0 && i_func == 1) f_sg->SetParameter(3,0.0);
	      // if(energy == 6 && i_pt > 0 && i_pt < 3 && i_norm == 0 && i_func == 1) f_sg->SetParameter(4,-1.5);
	      // if(energy == 6 && i_pt > 0 && i_pt < 3 && i_norm == 0 && i_func == 1) f_sg->SetParameter(5,1.5);
	      float Fit_start = vmsa::InvMass[pid]-normL*vmsa::Width[pid];
	      float Fit_stop  = vmsa::InvMass[pid]+normR*vmsa::Width[pid];
	      f_sg->SetRange(Fit_start,Fit_stop);
	      // f_sg->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
	      cout << "i_pt = " << i_pt << ", i_norm = " << i_norm << ", i_func = " << i_func << endl;
	      h_mMass_theta[KEY_theta]->Fit(f_sg,"LNRI");

	      ParFit_theta[KEY_theta].clear();
	      for(int i_par = 0; i_par < vmsa::FuncParNum[i_func]; ++i_par)
	      {
		ParFit_theta[KEY_theta].push_back(static_cast<float>(f_sg->GetParameter(i_par)));
	      }
	    }
	  }
	}
      }
    }

#if _PlotQA_
    // QA plots for Poly+Breit_Wignar fits for phi integrated InvMass
    TCanvas *c_mMass_theta = new TCanvas("c_mMass_theta","c_mMass_theta",10,10,1400,1400);
    c_mMass_theta->Divide(5,5);
    for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; i_pt++) // pt loop
    {
      c_mMass_theta->cd(vmsa::pt_rebin_start[energy][i_pt]+1);
      c_mMass_theta->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetLeftMargin(0.15);
      c_mMass_theta->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetBottomMargin(0.15);
      c_mMass_theta->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetTicks(1,1);
      c_mMass_theta->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetGrid(0,0);
      string KEY_theta_QA = Form("pt_%d_Centrality_%d_EtaGap_%d_2nd_%s_Norm_%d_Func_%d",i_pt,vmsa::Cent_start,vmsa::Eta_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA,vmsa::Func_QA);
      h_mMass_theta[KEY_theta_QA]->SetMarkerColor(1);
      h_mMass_theta[KEY_theta_QA]->SetMarkerStyle(24);
      h_mMass_theta[KEY_theta_QA]->SetMarkerSize(0.8);
      h_mMass_theta[KEY_theta_QA]->DrawCopy("PE");

      float Fit_start = vmsa::InvMass[pid]-normL*vmsa::Width[pid];
      float Fit_stop  = vmsa::InvMass[pid]+normR*vmsa::Width[pid];
      TF1 *f_sg = getSgFunc(vmsa::Func_QA,pid);
      for(int i_par = 0; i_par < vmsa::FuncParNum[vmsa::Func_QA]; ++i_par)
      {
	f_sg->SetParameter(i_par,ParFit_theta[KEY_theta_QA][i_par]);
      }
      f_sg->SetLineColor(2);
      f_sg->SetLineStyle(1);
      f_sg->SetLineWidth(2);
      f_sg->SetRange(Fit_start,Fit_stop);
      f_sg->DrawCopy("l same");

      TF1 *f_bg = getBgFunc(vmsa::Func_QA,pid);
      for(int i_par = 0; i_par < vmsa::FuncParNum[vmsa::Func_QA]-3;++i_par)
      {
	f_bg->SetParameter(i_par,ParFit_theta[KEY_theta_QA][i_par+3]);
      }
      f_bg->SetLineColor(4);
      f_bg->SetLineStyle(2);
      f_bg->SetLineWidth(4);
      f_bg->SetRange(Fit_start,Fit_stop);
      f_bg->DrawCopy("l same");

      string pT_range = Form("[%.2f,%.2f]",vmsa::pt_low[energy][i_pt],vmsa::pt_up[energy][i_pt]);
      plotTopLegend((char*)pT_range.c_str(),0.2,0.7,0.08,1,0.0,42,1);
      PlotLine(vmsa::InvMass_low[pid],vmsa::InvMass_high[pid],0.0,0.0,1,2,2);
    }
#endif

    // Poly+bw fits for phi differential InvMass
    vecFMap ParFit;
    for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; i_pt++) // pt loop
    {
      for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // Centrality loop
      {
	for(int i_eta = vmsa::Eta_start; i_eta < vmsa::Eta_stop; i_eta++) // EtaGap loop
	{
	  for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
	  {
	    for(int i_theta = vmsa::CTS_start; i_theta < vmsa::CTS_stop; i_theta++) // cos(theta*) loop
	    {
	      string KEY = Form("pt_%d_Centrality_%d_EtaGap_%d_CosThetaStar_%d_2nd_%s_Norm_%d",i_pt,i_cent,i_eta,i_theta,vmsa::mPID[pid].c_str(),i_norm);
	      for(int i_func = vmsa::Func_start; i_func < vmsa::Func_stop; ++i_func)
	      {
		string KEY_func = Form("pt_%d_Centrality_%d_EtaGap_%d_CosThetaStar_%d_2nd_%s_Norm_%d_Func_%d",i_pt,i_cent,i_eta,i_theta,vmsa::mPID[pid].c_str(),i_norm,i_func);
		h_mMass_func[KEY_func] = (TH1F*)h_mMass[KEY]->Clone(KEY_func.c_str());

		string KEY_theta = Form("pt_%d_Centrality_%d_EtaGap_%d_2nd_%s_Norm_%d_Func_%d",i_pt,i_cent,i_eta,vmsa::mPID[pid].c_str(),i_norm,i_func);
		TF1 *f_sg = getSgFunc(i_func,pid);
		f_sg->FixParameter(0,ParFit_theta[KEY_theta][0]);
		f_sg->FixParameter(1,ParFit_theta[KEY_theta][1]);
		f_sg->SetParameter(2,ParFit_theta[KEY_theta][2]/7.0);
		for(int i_par = 3; i_par < vmsa::FuncParNum[i_func]; ++i_par)
		{
		  f_sg->SetParameter(i_par,ParFit_theta[KEY_theta][i_par]/7.0);
		  // f_sg->FixParameter(i_par,ParFit_theta[KEY_theta][i_par]/7.0);
		}
		float Fit_start = vmsa::InvMass[pid]-normL*vmsa::Width[pid];
		float Fit_stop  = vmsa::InvMass[pid]+normR*vmsa::Width[pid];
		f_sg->SetRange(Fit_start,Fit_stop);
		// f_sg->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
		// cout << "i_pt = " << i_pt << ", i_theta = " << i_theta << ", i_norm = " << i_norm << ", i_func = " << i_func << endl;
		h_mMass_func[KEY_func]->Fit(f_sg,"LNQRI");
		ParFit[KEY_func].clear();
		for(int i_par = 0; i_par < vmsa::FuncParNum[i_func]; ++i_par)
		{
		  ParFit[KEY_func].push_back(static_cast<float>(f_sg->GetParameter(i_par)));
		}

		TF1 *f_bg = getBgFunc(i_func,pid);
		for(int i_par = 0; i_par < vmsa::FuncParNum[vmsa::Func_QA]-3;++i_par)
		{
		  f_bg->SetParameter(i_par,ParFit[KEY_func][i_par+3]);
		}

		h_mMass_func[KEY_func]->Add(f_bg,-1.0);
	      }
	    }
	  }
	}
      }
    }

#if _PlotQA_
    // QA plots for phi differential InvMass after linear background subtraction
    TCanvas *c_mMass_phi_diff = new TCanvas("c_mMass_phi_diff","c_mMass_phi_diff",10,10,1400,1400);
    c_mMass_phi_diff->Divide(5,5);
    for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; i_pt++) // pt loop
    {
      c_mMass_phi_diff->cd(vmsa::pt_rebin_start[energy][i_pt]+1);
      c_mMass_phi_diff->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetLeftMargin(0.15);
      c_mMass_phi_diff->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetBottomMargin(0.15);
      c_mMass_phi_diff->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetTicks(1,1);
      c_mMass_phi_diff->cd(vmsa::pt_rebin_start[energy][i_pt]+1)->SetGrid(0,0);
      string KEY_func_QA = Form("pt_%d_Centrality_%d_EtaGap_%d_CosThetaStar_%d_2nd_%s_Norm_%d_Func_%d",i_pt,vmsa::Cent_start,vmsa::Eta_start,vmsa::CTS_start,vmsa::mPID[pid].c_str(),vmsa::Norm_QA,vmsa::Func_QA);
      h_mMass_func[KEY_func_QA]->SetMarkerColor(1);
      h_mMass_func[KEY_func_QA]->SetMarkerStyle(24);
      h_mMass_func[KEY_func_QA]->SetMarkerSize(0.8);
      h_mMass_func[KEY_func_QA]->DrawCopy("PE");

      string pT_range = Form("[%.2f,%.2f]",vmsa::pt_low[energy][i_pt],vmsa::pt_up[energy][i_pt]);
      plotTopLegend((char*)pT_range.c_str(),0.2,0.7,0.08,1,0.0,42,1);
      PlotLine(vmsa::InvMass_low[pid],vmsa::InvMass_high[pid],0.0,0.0,1,2,2);
    }
#endif
  }

  // write background subtracted histograms to output file
  string outputfile = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/%s/rho00/InvMassSubBg.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  TFile *File_OutPut = new TFile(outputfile.c_str(),"RECREATE");
  File_OutPut->cd();
  for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; i_pt++) // pt loop
  {
    for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; i_cent++) // Centrality loop
    {
      for(int i_eta = vmsa::Eta_start; i_eta < vmsa::Eta_stop; i_eta++) // EtaGap loop
      {
	for(int i_theta = vmsa::CTS_start; i_theta < vmsa::CTS_stop; i_theta++) // cos(theta*) loop
	{
	  for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
	  {
	    for(int i_func = vmsa::Func_start; i_func < vmsa::Func_stop; ++i_func)
	    {
	      string KEY_func = Form("pt_%d_Centrality_%d_EtaGap_%d_CosThetaStar_%d_2nd_%s_Norm_%d_Func_%d",i_pt,i_cent,i_eta,i_theta,vmsa::mPID[pid].c_str(),i_norm,i_func);
	      h_mMass_func[KEY_func]->Write();
	    }
	  }
	}
      }
    }
  }
  File_OutPut->Close();
}

TF1* getSgFunc(int nFunc, int pid)
{
  TF1 *f_func;
  if(nFunc == 0) f_func = new TF1("f_func",PolyBreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],vmsa::FuncParNum[nFunc]);
  if(nFunc == 1) f_func = new TF1("f_func",ChebBreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],vmsa::FuncParNum[nFunc]);
  if(nFunc == 2) f_func = new TF1("f_func",LegeBreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],vmsa::FuncParNum[nFunc]);
  if(nFunc == 3) f_func = new TF1("f_func",SqRtBreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],vmsa::FuncParNum[nFunc]);
  for(int i_par = 0; i_par < vmsa::FuncParNum[nFunc]; ++i_par)
  {
    f_func->ReleaseParameter(i_par);
  }

  return f_func;
}

TF1* getBgFunc(int nFunc, int pid)
{
  TF1 *f_func;
  if(nFunc == 0) f_func = new TF1("f_func",Poly,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],vmsa::FuncParNum[nFunc]-3);
  if(nFunc == 1) f_func = new TF1("f_func",Cheb,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],vmsa::FuncParNum[nFunc]-3);
  if(nFunc == 2) f_func = new TF1("f_func",Lege,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],vmsa::FuncParNum[nFunc]-3);
  if(nFunc == 3) f_func = new TF1("f_func",SqRt,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],vmsa::FuncParNum[nFunc]-3);
  for(int i_par = 0; i_par < vmsa::FuncParNum[nFunc]-3; ++i_par)
  {
    f_func->ReleaseParameter(i_par);
  }

  return f_func;
}
